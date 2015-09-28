PROGRAM VARASSIM
! -----------------------------------------------------------------------------
!
! Land Data Assimilation System based on an Extended Kalman Filter
!
! Revised version : JFM (15 September 2008)
!
! Possibility to evolve the B matrix in the cycling - otherwise SEKF
!
! First version including patches (15 October 2008)
! Modifications to allow for LAI assimilation by C. Rüdiger (8 Jan 2009)
!               - added Biomass as control and LAI as observation variable

! Modifications to allow for SWI assimilation by A.Barbu (25 Sept 2009)
! The new state vector can be any element of (WG2, WG1, LAI, Biomass) - Choice in namelist
! The new observation vector can be any element of (WG1, LAI, SWI, FAPAR) - Choice in namelist
! Bug correction regarding the Kalman gain computation (otherwise underestimation of the gain)
! New formulation of the B matrix to ensure its positive definiteness
! -----------------------------------------------------------------------------
!
USE MODD_OFF_SURFEX_n
!
USE MODD_SURFEX_MPI, ONLY : NCOMM, NPROC, NRANK, NPIO, WLOG_MPI, PREP_LOG_MPI,   &
                            NINDEX, NSIZE_TASK, END_LOG_MPI
USE MODD_SURFEX_OMP, ONLY : NINDX2SFX, NWORK, XWORK, XWORK2, XWORK3, NBLOCK, NBLOCKTOT,&
                            NWORK_FULL, XWORK_FULL, XWORK2_FULL, INIT_DIM, RESET_DIM
USE MODD_MASK, ONLY: NMASK_FULL
                          
 USE MODD_TYPE_DATE_SURF
 USE MODD_SURF_PAR,ONLY : XUNDEF
 USE MODD_SURF_CONF, ONLY : CPROGNAME
 USE MODD_OL_FILEID
 USE MODD_IO_SURF_ASC,ONLY : CFILEIN, CFILEIN_SAVE, CFILEPGD, CFILEOUT 
!
 USE MODI_OPEN_NAMELIST
 USE MODI_CLOSE_NAMELIST
 USE MODE_POS_SURF 
 USE MODI_ABOR1_SFX
 USE MODI_INI_CSTS
 USE MODI_INI_DATA_COVER
 USE MODI_INIT_IO_SURF_n
 USE MODI_SET_SURFEX_FILEIN
 USE MODI_END_IO_SURF_n 
 USE MODI_READ_SURF
 USE MODI_GET_SIZE_FULL_n
 USE MODI_READ_COVER_n
 USE MODI_CONVERT_COVER_FRAC
 USE MODI_IO_BUFF_CLEAN
 USE MODI_WRITE_SURF
 USE MODI_ADD_FORECAST_TO_DATE_SURF
 USE MODI_CHOLDC
 USE MODI_CHOLSL
 USE MODI_INVERSE_MATRIX
 USE MODI_INIT_INDEX_MPI
!
 USE MODE_EKF, ONLY : GET_FILE_NAME
!
 USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
 USE PARKIND1  ,ONLY : JPRB
!
USE MODD_SURFEX_n
USE MODD_OFF_SURFEX_n
!
! -----------------------------------------------------------
!
 IMPLICIT NONE
!
#ifdef SFX_MPI
INCLUDE 'mpif.h'
#endif
!
 INTEGER                                  :: NOBSTYPE, NOBS
 CHARACTER(LEN=28)                        :: YNAMELIST = 'OPTIONS.nam '
 INTEGER, PARAMETER                       :: NVARMAX = 4      ! WG2, WG1, Biomass (order convention)
 INTEGER, PARAMETER                       :: NOBSMAX = 3      ! WG1, LAI (order convention)
 REAL,DIMENSION(NOBSMAX)                  :: XERROBS_M        ! Observational standard deviation (max dimension)
 REAL,DIMENSION(:),ALLOCATABLE            :: ERROBS          ! Observational standard deviation (observation dimension)
!
!    Declarations of local variables
!
 CHARACTER(LEN=3), PARAMETER                :: YINIT     = 'ALL'
 LOGICAL                                    :: LBFIXED 
 LOGICAL,DIMENSION(:,:),ALLOCATABLE         :: LPATCH                     ! logical vector for nb points, 12 patches and nb var control
 CHARACTER(LEN=6)                           :: YPROGRAM  = 'ASCII '       !!!!
 INTEGER                                    :: OBSCOUNT
 CHARACTER(LEN=1)                           :: LCHAR
 INTEGER                                    :: IYEAR                      ! current year (UTC)
 INTEGER                                    :: IMONTH                     ! current month (UTC)
 INTEGER                                    :: IDAY                       ! current day (UTC)
 INTEGER                                    :: IHOUR                      ! current hour (UTC)
 REAL                                       :: ZTIME                      ! current time since start of the run (s)
 REAL                                       :: PTSTEP_OUTPUT,PTSTEP_FORC  ! OUPUT step, Duration and FORCING Step
 INTEGER                                    :: IRESP,PATCH_NUMBER         ! return code
 INTEGER                                    :: ILOBS                      ! Namelist unit number
 TYPE (DATE_TIME)                           :: TTIME                      ! Current date and time
 INTEGER                                    :: INB_FORC                   ! NB Forcing step
 INTEGER                                    :: NBOUTPUT                   ! Number of time step
 INTEGER                                    :: ISTEP                      ! 
 REAL,DIMENSION(:,:,:),ALLOCATABLE          :: YF                         ! Vector of model observations (averaged)
 REAL,DIMENSION(:,:,:,:),ALLOCATABLE        :: YF_PATCH                   ! Vector of model observations (for each patch)
 REAL,DIMENSION(:,:,:,:),ALLOCATABLE        :: XF                         ! Vector of forecast control variables
 REAL,DIMENSION(:,:),ALLOCATABLE            :: LAI_PASS                   ! Storing the ini LAI value at the beginning of a simulation
 REAL,DIMENSION(:,:,:),ALLOCATABLE          :: XI                         ! Vector of control variables (at beginning of timestep)
 REAL,DIMENSION(:,:),ALLOCATABLE            :: ZBIO_PASS                  ! Vector of Biomass (not actual control var, though needs to be updated
 REAL,DIMENSION(:,:),ALLOCATABLE            :: ZBIO_OUT                   ! along with the assimilation of the LAI
 REAL,DIMENSION(12)                         :: XALPH
 REAL,DIMENSION(:,:,:),ALLOCATABLE          :: HO                         ! Jacobian of observation operator
 REAL,DIMENSION(:,:,:),ALLOCATABLE          :: HOWR                       ! copy of HO for writing out
 REAL,DIMENSION(:,:,:),ALLOCATABLE          :: HOT                        ! Transpose of HO
 REAL,DIMENSION(:,:,:),ALLOCATABLE          :: R                          ! covariance matrix of observation errors
 REAL,DIMENSION(:,:,:),ALLOCATABLE          :: B                          ! background error covariance matrix
 REAL,DIMENSION(:,:),ALLOCATABLE            :: IDENT                      ! identitiy matrix, used for Ba
 REAL, DIMENSION(:,:,:), ALLOCATABLE        :: GAIN                       ! Kalman gain (used explicitly for Ba)
 REAL, DIMENSION(:,:,:), ALLOCATABLE        :: LTM                        ! linear tangent matrix for the f'ward model
 REAL, DIMENSION(:,:,:), ALLOCATABLE        :: Q                          ! model error matrix
 REAL,DIMENSION(:,:),ALLOCATABLE            :: YO                         ! vector of observations
 REAL,DIMENSION(:),ALLOCATABLE              :: ZCLAY                      ! Percentage of clay (varies from 0 to 1)
 REAL,DIMENSION(:),ALLOCATABLE              :: ZSAND
 REAL,DIMENSION(:),ALLOCATABLE              :: COFSWI                     ! dynamic range (Wfc - Wwilt)
 REAL,DIMENSION(:),ALLOCATABLE              :: SMSAT                      ! saturation
 REAL,DIMENSION(:),ALLOCATABLE              :: WILT                       ! wilting point
 REAL,DIMENSION(:,:,:),ALLOCATABLE          :: ZEPS                       ! The perturbation amplitude
 REAL,DIMENSION(:,:),ALLOCATABLE            :: XINCR                      ! Analysis increment
 REAL,DIMENSION(:,:,:),ALLOCATABLE          :: SIMOBS                     ! Simulated temp, relative hum and LAI(available per patch?)
 REAL, DIMENSION (:,:),ALLOCATABLE          :: XPATCH                     ! Fraction covered by each patch
 REAL,DIMENSION(:,:,:),ALLOCATABLE          :: SIMMOD                     ! Control variables (for B propagation)
 REAL,DIMENSION(:,:),ALLOCATABLE            :: VECT                       ! The analysed variable
 
 CHARACTER(LEN=200)                         :: NMFILE_CANARI              ! Name of the observation, perturbed or reference file!
 CHARACTER(LEN=9)                           :: HFNAME
 CHARACTER(LEN=17)                          :: LFNAME
 CHARACTER(LEN=200)                         :: MINMAX
 REAL,DIMENSION(:,:),ALLOCATABLE            :: cpt                       ! index of abnormal LAI values
!
! Local Matrix for Analysis calculation
!
 REAL,DIMENSION(:,:,:),ALLOCATABLE          :: K1
 REAL,DIMENSION(:,:,:),ALLOCATABLE          :: KH,KRK
 REAL,DIMENSION(:,:,:),ALLOCATABLE          :: IdKH
 REAL,DIMENSION(:,:),ALLOCATABLE            :: ZX,ZB,ZP
 REAL,DIMENSION(:,:),ALLOCATABLE            :: YOWR
 INTEGER                                    :: JI,J,JJ,K,KK,L
 LOGICAL                                    :: LPRT                      ! Running VARASSIM in a perturbation mode
 LOGICAL                                    :: LSIM                      ! Running VARASSIM in a reading mode
 LOGICAL                                    :: LBEV                      ! Running VARASSIM to evolve B
 REAL,DIMENSION(NVARMAX)                    :: XTPRT_M                    ! The perturbation amplitude (max dim)
 REAL,DIMENSION(NVARMAX)                    :: XSIGMA_M                  ! covariance of background errors if B is fixed (max dim)
!                                                                        ! covariance of model errors if B evolving (max dim)
 CHARACTER(LEN=3),DIMENSION(NVARMAX)        :: CVAR_M !                  ! Name of contr var (syntax of surfex in PREP.txt file)(maxdim)
 CHARACTER(LEN=12)                          :: CBIO                      ! Name of Biomass variable
 CHARACTER(LEN=100),DIMENSION(NVARMAX)      :: CPREFIX_M                  ! The prefix of the control var. (in PREP.txt file) (max dim)
 CHARACTER(LEN=100)                         :: CPREFIX_BIO                ! The prefix of the Biomass variable
 INTEGER,DIMENSION(NVARMAX)                 :: NNCV                      ! Select the control variables to be used
 INTEGER,DIMENSION(NOBSMAX)                 :: NNCO                      ! Select the type of observations to be assimilated
! 
 REAL,DIMENSION(:),ALLOCATABLE              :: TPRT                      ! The perturbation amplitude
 REAL,DIMENSION(:),ALLOCATABLE              :: XSIGMA                    ! covariance of background errors if B is fixed
!                                                                        ! covariance of model errors if B evolving 
 CHARACTER(LEN=3),DIMENSION(:),ALLOCATABLE  :: XVAR !                    ! Name of control variables (syntax of surfex in PREP.txt file )
 CHARACTER(LEN=100),DIMENSION(:),ALLOCATABLE :: PREFIX                   ! The prefix of the control variables (in PREP.txt file) 
 CHARACTER(LEN=10),DIMENSION(NOBSMAX)        :: XOBS_M                   ! Identifier for simulated observations (max dimension)
 CHARACTER(LEN=10),DIMENSION(:),ALLOCATABLE :: XOBS                      ! Identifier for simulated observations
 REAL,DIMENSION(:,:),ALLOCATABLE            :: FROZ                      ! WGI1
 INTEGER                                    :: NIVAR                      ! counter for cntrl vars
 INTEGER                                    :: NVAR                      ! number of cntrl vars
 INTEGER                                    :: NDIM
 CHARACTER(LEN=100)                         :: HCOMMENT
 INTEGER                                    :: ILUOUT                    ! ascii output unit number
 INTEGER                                    :: ILUNAM                    ! namelist unit number
 INTEGER                                    :: ISTAT                    
 LOGICAL                                    :: GFOUND                    ! return logical when reading namelist
 REAL                                       :: XSCALE_Q,XSCALE_QLAI        ! scaling factor of Q matrix w.r.t. the initial B
 REAL,DIMENSION(:,:),ALLOCATABLE            :: WM,MWG2,XWG2              ! min and max value for wg/or swi2
 REAL                                       :: SCALE_SWI                 ! scaling factor of SWI depending on w2 
 REAL                                       :: WG2MIN                    ! minimal value of simulated wg2 in open-loop case
 REAL,DIMENSION(:),ALLOCATABLE              :: XVLAIMIN
 CHARACTER(LEN=50)                          :: OFORMAT
 REAL(KIND=JPRB)                            :: ZHOOK_HANDLE
 LOGICAL                                    :: LVARASSIM
!
!
! MPI variables
!
 CHARACTER(LEN=100) :: YNAME
INTEGER :: ILEVEL, INFOMPI, INKPROMA, JBLOCK
INTEGER, DIMENSION(:), ALLOCATABLE :: ISIZE_OMP

!
 NAMELIST/NAM_IO_VARASSIM/LPRT, LSIM, LBEV, LBFIXED
 NAMELIST/NAM_OBS/NOBSTYPE, XERROBS_M, NNCO
 NAMELIST/NAM_VAR/NIVAR, NVAR, CVAR_M, CPREFIX_M, XSIGMA_M, XTPRT_M, NNCV, XSCALE_Q, XSCALE_QLAI, CBIO, CPREFIX_BIO, XALPH
!
!     0.1.   MPI and dr_hook initializations
!
#ifdef SFX_MPI
 CALL MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,ILEVEL,INFOMPI)
#endif

 IF (LHOOK) CALL DR_HOOK('VARASSIM',0,ZHOOK_HANDLE)
!
#ifdef SFX_MPI
NCOMM = MPI_COMM_WORLD
 CALL MPI_COMM_SIZE(NCOMM,NPROC,INFOMPI)
 CALL MPI_COMM_RANK(NCOMM,NRANK,INFOMPI)
#endif
!
 CALL PREP_LOG_MPI

!--------------------------------------

!-------------------------------
!
!*     0.2.    Allocations of Surfex Types
!
! Allocation of Surfex Types
 CALL SURFEX_ALLOC_LIST(1)
!
 PRINT *
 PRINT *,'   --------------------------'
 PRINT *,'   |   ENTERING  VARASSIM   |'
 PRINT *,'   --------------------------'
 PRINT *
PRINT *,'V8 varassim '

!
!      0.3   Open namelist
!
 ILUOUT = 111
 CPROGNAME='ASCII '
!
NNCV(:) = 0
 CALL OPEN_NAMELIST(CPROGNAME,ILUNAM,YNAMELIST)
 CALL POSNAM(ILUNAM,'NAM_IO_VARASSIM',GFOUND,ILUOUT)
 IF (GFOUND) READ (UNIT=ILUNAM,NML=NAM_IO_VARASSIM)
 CALL CLOSE_NAMELIST(CPROGNAME,ILUNAM)
!
 CALL OPEN_NAMELIST(CPROGNAME,ILUNAM,YNAMELIST)
 CALL POSNAM(ILUNAM,'NAM_OBS',GFOUND,ILUOUT)
 IF (GFOUND) READ (UNIT=ILUNAM,NML=NAM_OBS)
 CALL CLOSE_NAMELIST(CPROGNAME,ILUNAM)
!
 CALL OPEN_NAMELIST(CPROGNAME,ILUNAM,YNAMELIST)
 CALL POSNAM(ILUNAM,'NAM_VAR',GFOUND,ILUOUT)
 IF (GFOUND) READ (UNIT=ILUNAM,NML=NAM_VAR)
 CALL CLOSE_NAMELIST(CPROGNAME,ILUNAM)
!
 CALL GOTO_MODEL(1)
!
 CALL INIT_INDEX_MPI(YSURF_CUR, &
                     CPROGNAME, 'LIN ', 0., .FALSE.)
!
!      0.4 Allocate arrays depending on the control vector dimension
!
 ALLOCATE (TPRT(NVAR))
 ALLOCATE (XSIGMA(NVAR))
 ALLOCATE (XVAR(NVAR))
 ALLOCATE (PREFIX(NVAR))
!
 IF (SUM(NNCV) /= NVAR) THEN
  PRINT *,' INCONSISTENCY in set-up of CONTROL VARIABLES',SUM(NNCV)
  PRINT *,'vars specified in namelist: ', sum(NNCV)
  PRINT *,'NVAR from runsh: ', NVAR
  CALL ABOR1_SFX("VARASSIM:INCONSISTENCY in set-up of CONTROL VARIABLES")
 ENDIF
!
! select relevant control vars
 J = 1
 DO JI = 1,NVARMAX
  IF (NNCV(JI) == 1 .AND. J <= NVAR ) THEN
   TPRT(J) = XTPRT_M(JI)
   XSIGMA(J) = XSIGMA_M(JI)
   XVAR(J) = CVAR_M(JI)
   PREFIX(J) = CPREFIX_M(JI)
   J = J + 1
  ENDIF
 ENDDO
!
!      0.5 Allocate arrays depending on the number of observation type
!
 ALLOCATE (XOBS(NOBSTYPE))
 ALLOCATE (ERROBS(NOBSTYPE))
!
!       assigning of observation variable names farther down
!
!      1.    Initializations
!
 ILOBS=55
!
 CALL INI_CSTS
!
 CALL INI_DATA_COVER(YSURF_CUR%DTCO, YSURF_CUR%U)

 CFILEIN = 'PREP.txt'
 CFILEPGD = 'PGD.txt'
 CFILEIN_SAVE = CFILEIN
!
!   Read grid dimension for allocation
!
CALL INIT_IO_SURF_n(YSURF_CUR%DTCO, YSURF_CUR%DGU, YSURF_CUR%U, &
                        YPROGRAM,'FULL  ','SURF  ','READ ')
!
!   Find current time
!
 CALL READ_SURF(&
                YPROGRAM,'DTCUR',TTIME,IRESP)
print *, 'Time from PREP file : ', TTIME 
!
 CALL END_IO_SURF_n(YPROGRAM)
 CALL SET_SURFEX_FILEIN(YPROGRAM,'PGD ') ! change input file name to pgd name
CALL INIT_IO_SURF_n(YSURF_CUR%DTCO, YSURF_CUR%DGU, YSURF_CUR%U, &
                        YPROGRAM,'FULL  ','SURF  ','READ ')
!
!   Reading grid characteristics to perform nature mask
!
 CALL READ_SURF(&
                YPROGRAM,'SEA   ',YSURF_CUR%U%CSEA   ,IRESP)
 CALL READ_SURF(&
                YPROGRAM,'WATER ',YSURF_CUR%U%CWATER ,IRESP)
 CALL READ_SURF(&
                YPROGRAM,'NATURE',YSURF_CUR%U%CNATURE,IRESP)
 CALL READ_SURF(&
                YPROGRAM,'TOWN  ',YSURF_CUR%U%CTOWN  ,IRESP)

!
 CALL READ_SURF(&
                YPROGRAM,'DIM_FULL  ',YSURF_CUR%U%NDIM_FULL,  IRESP)
 CALL READ_SURF(&
                YPROGRAM,'DIM_SEA   ',YSURF_CUR%U%NDIM_SEA,   IRESP)
 CALL READ_SURF(&
                YPROGRAM,'DIM_NATURE',YSURF_CUR%U%NDIM_NATURE,IRESP)
 CALL READ_SURF(&
                YPROGRAM,'DIM_WATER ',YSURF_CUR%U%NDIM_WATER, IRESP)
 CALL READ_SURF(&
                YPROGRAM,'DIM_TOWN  ',YSURF_CUR%U%NDIM_TOWN,  IRESP)

NINDX2SFX = YSURF_CUR%U%NDIM_FULL
ALLOCATE(NWORK(YSURF_CUR%U%NDIM_FULL))
ALLOCATE(XWORK(YSURF_CUR%U%NDIM_FULL))
ALLOCATE(XWORK2(YSURF_CUR%U%NDIM_FULL,10))
ALLOCATE(XWORK3(YSURF_CUR%U%NDIM_FULL,10,10))
IF (NRANK==NPIO) THEN
  ALLOCATE(NWORK_FULL(YSURF_CUR%U%NDIM_FULL))
  ALLOCATE(XWORK_FULL(YSURF_CUR%U%NDIM_FULL))
  ALLOCATE(XWORK2_FULL(YSURF_CUR%U%NDIM_FULL,10))
ELSE
  ALLOCATE(NWORK_FULL(0))
  ALLOCATE(XWORK_FULL(0))
  ALLOCATE(XWORK2_FULL(0,0))
ENDIF
!
!   Get number of points on this proc
!
 CALL GET_SIZE_FULL_n(YSURF_CUR%U, &
                      YPROGRAM,YSURF_CUR%U%NDIM_FULL,YSURF_CUR%U%NSIZE_FULL)
 CALL READ_COVER_n(YSURF_CUR%DTCO, YSURF_CUR%U, &
                   YPROGRAM)
 CALL END_IO_SURF_n(YPROGRAM)
!
!   Perform masks (only nature used)
!
 ALLOCATE(YSURF_CUR%U%XSEA(YSURF_CUR%U%NDIM_FULL))
 ALLOCATE(YSURF_CUR%U%XNATURE(YSURF_CUR%U%NDIM_FULL))
 ALLOCATE(YSURF_CUR%U%XWATER(YSURF_CUR%U%NDIM_FULL))
 ALLOCATE(YSURF_CUR%U%XTOWN(YSURF_CUR%U%NDIM_FULL))
!
 CALL CONVERT_COVER_FRAC(YSURF_CUR%DTCO, &
                         YSURF_CUR%U%XCOVER,YSURF_CUR%U%LCOVER,YSURF_CUR%U%XSEA,&
                         YSURF_CUR%U%XNATURE,YSURF_CUR%U%XTOWN,YSURF_CUR%U%XWATER)
!
 YSURF_CUR%U%NSIZE_NATURE = COUNT(YSURF_CUR%U%XNATURE(:) > 0.0)
 YSURF_CUR%U%NSIZE_TOWN   = COUNT(YSURF_CUR%U%XTOWN(:)   > 0.0)
 YSURF_CUR%U%NSIZE_WATER  = COUNT(YSURF_CUR%U%XWATER(:)  > 0.0)
 YSURF_CUR%U%NSIZE_SEA    = COUNT(YSURF_CUR%U%XSEA(:)    > 0.0)
!
!   Read number of patches
!
 CALL SET_SURFEX_FILEIN(YPROGRAM,'PGD ')
CALL INIT_IO_SURF_n(YSURF_CUR%DTCO, YSURF_CUR%DGU, YSURF_CUR%U, &
                        YPROGRAM,'NATURE','ISBA  ','READ ')
 CALL READ_SURF(&
                YPROGRAM,'PATCH_NUMBER',PATCH_NUMBER,IRESP)
!
 ALLOCATE(ZSAND(YSURF_CUR%U%NSIZE_NATURE))
 ALLOCATE(ZCLAY(YSURF_CUR%U%NSIZE_NATURE))
!
!   Read CLAY fraction to compute the SWI range (Wfc - Wwilt)
!   (XSIGMA is defined in terms of SWI), need to convert to equivalent v/v
!   using same clay fraction in both layers
!   Read SAND fraction to compute the saturation for conversion of ERS SWI
!
 CALL READ_SURF(&
                YPROGRAM,'CLAY',ZCLAY(:),IRESP)
 CALL READ_SURF(&
                YPROGRAM,'SAND',ZSAND(:),IRESP)
!
!   Define prefixes for simulated observations
! 
 CALL END_IO_SURF_n(YPROGRAM)
 CALL SET_SURFEX_FILEIN(YPROGRAM,'PREP') ! change input file name to pgd name
CALL INIT_IO_SURF_n(YSURF_CUR%DTCO, YSURF_CUR%DGU, YSURF_CUR%U, &
                        YPROGRAM,'NATURE','ISBA  ','READ ') 
!
 XOBS_M(1) = 'WG1'
 XOBS_M(2) = 'LAI'
!
 J = 1
 DO JI = 1,NOBSMAX
  IF (NNCO(JI) == 1 .AND. J <= NOBSTYPE ) THEN
   ERROBS(J) = XERROBS_M(JI)
   XOBS(J) = XOBS_M(JI)
   J = J + 1
  ENDIF
 ENDDO
!
 ALLOCATE(XI(YSURF_CUR%U%NSIZE_NATURE,PATCH_NUMBER,NVAR))
 ALLOCATE(VECT(YSURF_CUR%U%NSIZE_NATURE,PATCH_NUMBER))
 ALLOCATE(LAI_PASS(YSURF_CUR%U%NSIZE_NATURE,PATCH_NUMBER))
 ALLOCATE(ZBIO_PASS(YSURF_CUR%U%NSIZE_NATURE,PATCH_NUMBER))
 ALLOCATE(ZBIO_OUT(YSURF_CUR%U%NSIZE_NATURE,PATCH_NUMBER))
!
 ! Avoiding unsuitable LAI values
 ALLOCATE(XVLAIMIN(PATCH_NUMBER))
 XVLAIMIN(:) = 0.3
 XVLAIMIN(5:6) = 1.0
!
! Read in control variables
!
 DO L = 1,NVAR
   CALL READ_SURF(&
                YPROGRAM,XVAR(L),XI(:,:,L),IRESP)
   IF (XVAR(L) .EQ. 'LAI') THEN
      LAI_PASS=XI(:,:,L)
      IF (CBIO .NE. 'LAI') THEN
        CALL READ_SURF(&
                YPROGRAM,CBIO,ZBIO_PASS(:,:),IRESP)
      ELSE
        ZBIO_PASS(:,:) = XI(:,:,L)
      ENDIF
   ENDIF
 ENDDO
 WHERE (LAI_PASS(:,:)==XUNDEF) LAI_PASS(:,:)=0.
!
ALLOCATE(XPATCH(YSURF_CUR%U%NSIZE_NATURE,PATCH_NUMBER))
ALLOCATE(LPATCH(YSURF_CUR%U%NSIZE_NATURE,PATCH_NUMBER))
!
! Read fraction of each patch (=> LPGD=T in OPTIONS.nam)
!
 CALL END_IO_SURF_n(YPROGRAM)
 CALL IO_BUFF_CLEAN
 CALL SET_SURFEX_FILEIN(YPROGRAM,'PREP')
CALL INIT_IO_SURF_n(YSURF_CUR%DTCO, YSURF_CUR%DGU, YSURF_CUR%U, &
                        YPROGRAM,'NATURE','ISBA  ','READ ')  
 IF (PATCH_NUMBER > 1) THEN
   CALL READ_SURF(&
                YPROGRAM,'PATCH',XPATCH,IRESP)
 ELSE
   XPATCH(:,1) = 1.0 
 ENDIF
!
LPATCH(:,:) = .FALSE.
!
   DO JI=1,YSURF_CUR%U%NSIZE_NATURE
     DO J=1,PATCH_NUMBER
         IF (XPATCH(JI,J) > 0.0) THEN
           LPATCH(JI,J) = .TRUE.
         ENDIF
     ENDDO
   ENDDO
!
 ALLOCATE(COFSWI(YSURF_CUR%U%NSIZE_NATURE))
 ALLOCATE(WILT(YSURF_CUR%U%NSIZE_NATURE))
 ALLOCATE(SMSAT(YSURF_CUR%U%NSIZE_NATURE)) 
 ALLOCATE(B(YSURF_CUR%U%NSIZE_NATURE,PATCH_NUMBER*NVAR,PATCH_NUMBER*NVAR))
!
 DO JI=1,YSURF_CUR%U%NSIZE_NATURE
  COFSWI(JI)=0.001*(89.0467*((100.*ZCLAY(JI))**0.3496)-37.1342*((100.*ZCLAY(JI))**0.5))
  SMSAT(JI)=0.001*(-1.08*100*ZSAND(JI)+494.305) 
  WILT(JI)=0.001*37.1342*((100.*ZCLAY(JI))**0.5)
 ENDDO
!
!   Frequency of assimilation cycling and data availability
!
 IHOUR = 6
 NBOUTPUT = 1
 PTSTEP_OUTPUT = 24
!
!   Time
!
 IYEAR  = TTIME%TDATE%YEAR
 IMONTH = TTIME%TDATE%MONTH
 IDAY   = TTIME%TDATE%DAY
 ZTIME  = TTIME%TIME
!
! ====================================================================
!
! ----------------------
! VARASSIM OPTION : LPRT 
! ----------------------
!
! ====================================================================
!   Compute perturbations for control variables
!   (write to PREP.* then exit)
!   Element of the control variable chosen by NIVAR in namelist
!
 IF ( LPRT ) THEN

   PRINT *,'   ------------------------------------'
   PRINT *,'   |   ENTERING  VARASSIM WITH LPRT   |'
   PRINT *,'   ------------------------------------'
! read in control variable
   CFILEOUT='SURFOUT.txt'
   DO JI=1,YSURF_CUR%U%NSIZE_NATURE
     DO J=1,PATCH_NUMBER
       IF (XI(JI,J,NIVAR).NE.XUNDEF) THEN                                 ! check whether values are undefined
          VECT(JI,J) = XI(JI,J,NIVAR) + TPRT(NIVAR)*XI(JI,J,NIVAR) 
       ELSE
          VECT(JI,J) = XI(JI,J,NIVAR)
       ENDIF
       IF (XI(JI,J,NIVAR).NE.XUNDEF .AND. XVAR(NIVAR).EQ.'LAI') THEN
          ZBIO_OUT(JI,J) = ZBIO_PASS(JI,J) + TPRT(NIVAR)*ZBIO_PASS(JI,J)
       ELSE
          ZBIO_OUT(JI,J) = ZBIO_PASS(JI,J)
       ENDIF
     ENDDO
   ENDDO

   CALL END_IO_SURF_n(YPROGRAM)
   CALL IO_BUFF_CLEAN
CALL INIT_IO_SURF_n(YSURF_CUR%DTCO, YSURF_CUR%DGU, YSURF_CUR%U, &
                        YPROGRAM,'NATURE','ISBA  ','WRITE')     
   CALL WRITE_SURF(YSURF_CUR%DGU, YSURF_CUR%U, &
                   YPROGRAM,XVAR(NIVAR),VECT,IRESP,HCOMMENT=PREFIX(NIVAR))
   IF (XVAR(NIVAR).EQ.'LAI') THEN
     CALL WRITE_SURF(YSURF_CUR%DGU, YSURF_CUR%U, &
                   YPROGRAM,CBIO,ZBIO_OUT,IRESP,HCOMMENT=CPREFIX_BIO)
   ENDIF
   CALL END_IO_SURF_n(YPROGRAM)
   CALL IO_BUFF_CLEAN
!
! write out perturbation
   OPEN (unit=111,file='PERTURB',status='unknown',IOSTAT=istat)
   DO JI=1,YSURF_CUR%U%NSIZE_NATURE
     DO J=1,PATCH_NUMBER
       IF (XI(JI,J,NIVAR).NE.XUNDEF) THEN
         WRITE (111,*) TPRT(NIVAR)*XI(JI,J,NIVAR)
       ELSE
         WRITE (111,*) 1.0
       ENDIF
     ENDDO
   ENDDO
   CLOSE(111)
!
! Initialisation of B matrix
!
   B(:,:,:) = 0.0
   DO L=1,NVAR
     DO JI=1,YSURF_CUR%U%NSIZE_NATURE
       DO J=1,PATCH_NUMBER
         IF (XVAR(L) .EQ. 'WG1' .OR. XVAR(L) .EQ. 'WG2') THEN
          B(JI,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(L-1)) = XSIGMA(L)*XSIGMA(L)*COFSWI(JI)*COFSWI(JI)
         ELSEIF (XVAR(L) .EQ. 'LAI') THEN 
          B(JI,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(L-1)) = XSIGMA(L)*XSIGMA(L)!*LAI_PASS(I,J)*LAI_PASS(I,J)
         ELSE
           B(JI,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(L-1)) = XSIGMA(L)*XSIGMA(L)
         ENDIF
       ENDDO
     ENDDO
   ENDDO
!
   !write out initial B matrix
   OPEN (unit=111,file='BGROUNDin0',status='new',IOSTAT=istat)
   IF (istat .NE. 0) THEN
     STOP 'BGROUNDin0 already written'
   ELSE
     DO JI=1,YSURF_CUR%U%NSIZE_NATURE
       DO J=1,PATCH_NUMBER
         DO JJ=1,PATCH_NUMBER
           DO L=1,NVAR
             DO K=1,NVAR
               WRITE (111,*) B(JI,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(L-1))
             ENDDO
           ENDDO
         ENDDO
       ENDDO
     ENDDO
     CLOSE(111)
   ENDIF
!
   PRINT *
   PRINT *,'   -----------------------------------'
   PRINT *,'   |   EXITING VARASSIM AFTER LPRT   |'
   PRINT *,'   -----------------------------------'
   PRINT *
 ENDIF
! ====================================================================
!
! ----------------------
! VARASSIM OPTION : LSIM
! ----------------------
!
! ====================================================================
!   Store simulated observations and evolved prognostic variables. 
!   This option is called (NVAR+1) times:
!   1) with reference initial state  (once)
!   2) with perturbed initial states (NVAR times)
!
 
 IF ( LSIM ) THEN
!
   PRINT *,'   ----------------------------------'
   PRINT *,'   |   ENTERING VARASSIM WITH LSIM   |'
   PRINT *,'   ----------------------------------'
   CALL END_IO_SURF_n(YPROGRAM)
   CALL IO_BUFF_CLEAN
   CALL SET_SURFEX_FILEIN(YPROGRAM,'PREP')
CALL INIT_IO_SURF_n(YSURF_CUR%DTCO, YSURF_CUR%DGU, YSURF_CUR%U, &
                        YPROGRAM,'NATURE','ISBA  ','READ ')     
!
   ALLOCATE(SIMOBS(YSURF_CUR%U%NSIZE_NATURE,PATCH_NUMBER,NOBSTYPE))
   ALLOCATE(SIMMOD(YSURF_CUR%U%NSIZE_NATURE,PATCH_NUMBER,NVAR))
!
   !read observation from file
   DO K = 1,NOBSTYPE
     CALL READ_SURF(&
                YPROGRAM,XOBS(K),SIMOBS(:,:,K),IRESP) 
     IF( XOBS(K) == 'SWI2' ) THEN
       !read file containing wg2 minmax values -> use for normalization
       ALLOCATE(MWG2(YSURF_CUR%U%NSIZE_NATURE,PATCH_NUMBER))
       ALLOCATE(XWG2(YSURF_CUR%U%NSIZE_NATURE,PATCH_NUMBER))
       OPEN(UNIT=111,FILE='mxwg2.dat',FORM='FORMATTED',STATUS='OLD',IOSTAT=istat)
       IF (ISTAT.EQ.0) THEN
         DO JI=1,YSURF_CUR%U%NSIZE_NATURE
           DO J=1,PATCH_NUMBER
             READ(111,*) MWG2(JI,J),XWG2(JI,J)
             WG2MIN = MWG2(JI,J)
             SCALE_SWI = XWG2(JI,J) - MWG2(JI,J)
             IF( LPATCH(JI,J).AND. SCALE_SWI > 1.0E-6 ) THEN
               SIMOBS(JI,J,K)=(SIMOBS(JI,J,K)*COFSWI(JI)+WILT(JI)-WG2MIN)/SCALE_SWI
             ENDIF
           ENDDO
         ENDDO
         CLOSE(111)
       ELSE
         PRINT *, 'No mxwg2.dat file, use namelist given WG2min and scale values'
         DO JI=1,YSURF_CUR%U%NSIZE_NATURE
           DO J=1,PATCH_NUMBER
             SIMOBS(JI,J,K)=(SIMOBS(JI,J,K)*COFSWI(JI)+WILT(JI)-WG2MIN)/SCALE_SWI
           ENDDO
         ENDDO
       ENDIF
       DEALLOCATE(MWG2)
       DEALLOCATE(XWG2)
     ENDIF
   ENDDO

   CALL END_IO_SURF_n(YPROGRAM)
   CALL IO_BUFF_CLEAN
   CALL SET_SURFEX_FILEIN(YPROGRAM,'PREP')
CALL INIT_IO_SURF_n(YSURF_CUR%DTCO, YSURF_CUR%DGU, YSURF_CUR%U, &
                        YPROGRAM,'NATURE','ISBA  ','READ ')     
!
   OPEN (unit=111,file='OBSIMU',status='unknown')
   DO JI=1, YSURF_CUR%U%NSIZE_NATURE
     DO J=1,PATCH_NUMBER
       WRITE (111,*) (SIMOBS(JI,J,K),K = 1,NOBSTYPE)
     ENDDO
   ENDDO
   CLOSE(111)
!
   DO JI = 1, NVAR
     CALL READ_SURF( &
                YPROGRAM,XVAR(JI),SIMMOD(:,:,JI),IRESP)
   ENDDO
   OPEN (unit=111,file='MDSIMU',status='unknown',IOSTAT=istat)
   DO JI=1, YSURF_CUR%U%NSIZE_NATURE
     DO J=1,PATCH_NUMBER
       WRITE (111,*) (SIMMOD(JI,J,L),L=1,NVAR)
     ENDDO  
   ENDDO
   CLOSE(111)
! 
   PRINT *
   PRINT *,'   -----------------------------------'
   PRINT *,'   |   EXITING VARASSIM AFTER LSIM   |'
   PRINT *,'   -----------------------------------'
   PRINT *
 ENDIF
!
! ====================================================================
!
! ----------------------
! VARASSIM OPTION : LBEV
! ----------------------
!
! ====================================================================
!   Calculate the LTM, and evolve B. 
!
 IF (LBEV) THEN 
!
   PRINT *,'   -----------------------------------'
   PRINT *,'   |   ENTERING VARASSIM WITH LBEV   |'
   PRINT *,'   -----------------------------------'

   ALLOCATE(LTM( YSURF_CUR%U%NSIZE_NATURE,PATCH_NUMBER*NVAR,PATCH_NUMBER*NVAR))
   ALLOCATE(ZEPS( YSURF_CUR%U%NSIZE_NATURE,PATCH_NUMBER,NVAR))
   ALLOCATE(XF( YSURF_CUR%U%NSIZE_NATURE,PATCH_NUMBER,NVAR+1,NVAR))

!
   PRINT *, 'evolving B to time ', IYEAR,IMONTH,IDAY,IHOUR
!
 !read B matrix from previous analysis step: B(t-1)
  OPEN (unit=111,file='BGROUNDin',status='old',IOSTAT=istat)
  DO L=1,NVAR   ! control variable (x at previous time step)
     DO K=1,NVAR
       DO JI=1, YSURF_CUR%U%NSIZE_NATURE
         DO J=1,PATCH_NUMBER
           DO JJ=1,PATCH_NUMBER
             READ (111,*) B(JI,J+PATCH_NUMBER*(L-1),JJ+PATCH_NUMBER*(K-1))
           ENDDO
         ENDDO
       ENDDO
     ENDDO
   ENDDO
   CLOSE(111)
!
! Calculate the TLM of the forecast model
! 
! a) read in perturbed forecasts
!
   DO L=1,NVAR
     WRITE(LCHAR,'(I1) ') L
     NMFILE_CANARI='MDSIMU_PERT_'//LCHAR//'_'
     CALL GET_FILE_NAME(IYEAR,IMONTH,IDAY,IHOUR,NMFILE_CANARI)
     OPEN(UNIT=111,FILE=TRIM(NMFILE_CANARI)//".DAT",FORM='FORMATTED',STATUS='OLD')
     DO JI=1, YSURF_CUR%U%NSIZE_NATURE
       DO J=1,PATCH_NUMBER 
         READ(111,*) (XF(JI,J,L+1,K),K=1,NVAR)
       ENDDO
     ENDDO
     CLOSE(111)
     PRINT *, 'read in ptbd forecasts for L: ', NMFILE_CANARI(1:LEN_TRIM(NMFILE_CANARI))
   ENDDO
!
! b) read in reference forecasts
!
   NMFILE_CANARI='MDSIMU_REFR_'
   CALL GET_FILE_NAME(IYEAR,IMONTH,IDAY,IHOUR,NMFILE_CANARI)
   OPEN(UNIT=111,FILE=TRIM(NMFILE_CANARI)//".DAT",FORM='FORMATTED',STATUS='OLD')
   DO JI=1, YSURF_CUR%U%NSIZE_NATURE
     DO J=1,PATCH_NUMBER
       READ(111,*) (XF(JI,J,1,K),K=1,NVAR)
     ENDDO
   ENDDO
   CLOSE(111)
!
! c) read initial perturbation
!
   DO L=1,NVAR
     WRITE(LCHAR,'(I1) ') L
     NMFILE_CANARI='PERTURB_'//LCHAR//'_'
     CALL GET_FILE_NAME(IYEAR,IMONTH,IDAY,IHOUR,NMFILE_CANARI)
     PRINT *, 'read in PERTURBation: ', NMFILE_CANARI(1:LEN_TRIM(NMFILE_CANARI))
     OPEN(UNIT=111,FILE=TRIM(NMFILE_CANARI)//".DAT",FORM='FORMATTED',STATUS='OLD')
     DO JI=1, YSURF_CUR%U%NSIZE_NATURE
       DO J=1,PATCH_NUMBER
         READ(111,*) ZEPS(JI,J,L)
       ENDDO
     ENDDO
     CLOSE(111)
   ENDDO
!
! d) calculate LTM
!
LTM(:,:,:) = 0.0
DO L=1,NVAR    ! control variable (x at previous time step)
  DO K=1,NVAR
    DO JI=1, YSURF_CUR%U%NSIZE_NATURE 
      DO J=1,PATCH_NUMBER 
           IF (LPATCH(JI,J) .AND. XF(JI,J,L+1,K).NE.XUNDEF .AND. XF(JI,J,1,K).NE.XUNDEF ) THEN
             ! Jacobian of fwd model
             LTM(JI,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(K-1)) = (XF(JI,J,L+1,K) - XF(JI,J,1,K))/ZEPS(JI,J,L)
             LTM(JI,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(K-1)) = max(-0.1, LTM(JI,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(K-1)))! impose upper/lower limits 
             LTM(JI,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(K-1)) = min(1.0, LTM(JI,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(K-1)))
           ENDIF
         ENDDO
       ENDDO
     ENDDO
   ENDDO
!
! e) evolve B 
!
   DO JI=1, YSURF_CUR%U%NSIZE_NATURE
     B(JI,:,:)=MATMUL(LTM(JI,:,:),MATMUL(B(JI,:,:),TRANSPOSE(LTM(JI,:,:))))
   ENDDO
!
! write out the LTM for the forward model
!
   DO L=1,NVAR
     DO K=1,NVAR 
       LFNAME='LTM_del'//XVAR(K)//'_del'//XVAR(L)
       OPEN(UNIT=111,FILE=LFNAME,FORM='FORMATTED',STATUS='UNKNOWN',ACCESS='APPEND')
       DO JI=1, YSURF_CUR%U%NSIZE_NATURE
         DO J=1,PATCH_NUMBER
           WRITE (111,*) LTM(JI,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(K-1))
         ENDDO
       ENDDO
       CLOSE(111)
     ENDDO
   ENDDO
!
! write out current B (for use in next cycle)
!
  PRINT *,'writing out B in BGROUNDout'
  OPEN (unit=111,file='BGROUNDout',status='unknown')
  DO L=1,NVAR
    DO K=1,NVAR
      DO JI=1, YSURF_CUR%U%NSIZE_NATURE
        DO J=1,PATCH_NUMBER
           DO JJ=1,PATCH_NUMBER
             WRITE (111,*)  B(JI,J+PATCH_NUMBER*(L-1),JJ+PATCH_NUMBER*(K-1))
           ENDDO
         ENDDO
       ENDDO
     ENDDO
   ENDDO
   CLOSE(111)
   DEALLOCATE(B)
   DEALLOCATE(LTM)
   DEALLOCATE(ZEPS)
   DEALLOCATE(XF)
!
   PRINT *,'   -----------------------------------'
   PRINT *,'   |   EXITING VARASSIM AFTER LBEV   |'
   PRINT *,'   -----------------------------------'
 ENDIF
! ====================================================================
!
! if not LSIM, LPTR, or LBEV proceed with analysis
!
!   Number of available observation files
 PRINT *,'--> through to analysis section'
 PRINT *,'   ---------------------------------------'
 PRINT *,'   |   ENTERING VARASSIM WITH ANALYSIS   |'
 PRINT *,'   ---------------------------------------'
!
 NOBS=0
 ZTIME=PTSTEP_OUTPUT
!
 DO ISTEP=1,NBOUTPUT
   CALL ADD_FORECAST_TO_DATE_SURF(IYEAR,IMONTH,IDAY,ZTIME)
   ZTIME = ZTIME + PTSTEP_OUTPUT
   NMFILE_CANARI='CANARI_NATURE_'
   CALL GET_FILE_NAME(IYEAR,IMONTH,IDAY,IHOUR,NMFILE_CANARI)
   print *,'observation file name: ',NMFILE_CANARI(1:LEN_TRIM(NMFILE_CANARI))//".DAT"
   OPEN(UNIT=ILOBS,FILE=TRIM(NMFILE_CANARI)//".DAT",FORM='UNFORMATTED',STATUS='OLD',ERR=22)
   NOBS=NOBS+NOBSTYPE
   22 CONTINUE
   CLOSE(ILOBS)
 ENDDO
!
 IF (NOBS == 0) THEN
   PRINT *,'No observations read in OBS file - stop'
   STOP
 ELSE
   PRINT *,'Nobs=',NOBS
 ENDIF
!
!   Time reinitialization 
!
 IYEAR  = TTIME%TDATE%YEAR
 IMONTH = TTIME%TDATE%MONTH
 IDAY   = TTIME%TDATE%DAY
 ZTIME  = TTIME%TIME
!
!  Allocation
!
!  Perturbed simulations
!
 ALLOCATE(YF( YSURF_CUR%U%NSIZE_NATURE,NVAR+1,NOBS))
 ALLOCATE(YF_PATCH( YSURF_CUR%U%NSIZE_NATURE,PATCH_NUMBER,NVAR+1,NOBS))
 ALLOCATE(ZEPS( YSURF_CUR%U%NSIZE_NATURE,PATCH_NUMBER,NVAR))
 ALLOCATE(XF( YSURF_CUR%U%NSIZE_NATURE,PATCH_NUMBER,NVAR+1,NVAR))
!
! Read fraction of each patch (=> LPGD=T in OPTIONS.nam)
!
 CALL END_IO_SURF_n(YPROGRAM)
 CALL IO_BUFF_CLEAN
 CALL SET_SURFEX_FILEIN(YPROGRAM,'PREP')
CALL INIT_IO_SURF_n(YSURF_CUR%DTCO, YSURF_CUR%DGU, YSURF_CUR%U, &
                        YPROGRAM,'NATURE','ISBA  ','READ ')  
 IF (PATCH_NUMBER > 1) THEN
   CALL READ_SURF( &
                YPROGRAM,'PATCH',XPATCH,IRESP)
 ELSE
   XPATCH(:,1) = 1.0 
 ENDIF
!
!   Initial values (to be analysed)
!   Observations
!
 ALLOCATE(YO(YSURF_CUR%U%NSIZE_NATURE,NOBS))
 ALLOCATE(YOWR(YSURF_CUR%U%NDIM_FULL,NOBS))
!
!   Temporary vectors used by the EKF approach
!
 ALLOCATE(HO(YSURF_CUR%U%NSIZE_NATURE,NOBS,PATCH_NUMBER*NVAR))
 ALLOCATE(HOWR(YSURF_CUR%U%NSIZE_NATURE,NOBS,PATCH_NUMBER*NVAR))
 ALLOCATE(HOT(YSURF_CUR%U%NSIZE_NATURE,PATCH_NUMBER*NVAR,NOBS))
 ALLOCATE(R(YSURF_CUR%U%NSIZE_NATURE,NOBS,NOBS))
 ALLOCATE(IDENT(PATCH_NUMBER*NVAR,PATCH_NUMBER*NVAR))
 ALLOCATE(GAIN(YSURF_CUR%U%NSIZE_NATURE,PATCH_NUMBER*NVAR,NOBS))
 ALLOCATE(ZB(YSURF_CUR%U%NSIZE_NATURE,NOBS))
 ALLOCATE(ZX(YSURF_CUR%U%NSIZE_NATURE,NOBS))
 ALLOCATE(ZP(YSURF_CUR%U%NSIZE_NATURE,NOBS))
 ALLOCATE(XINCR(YSURF_CUR%U%NSIZE_NATURE,PATCH_NUMBER*NVAR))
 ALLOCATE(K1(YSURF_CUR%U%NSIZE_NATURE,NOBS,NOBS))
 ALLOCATE(KH(YSURF_CUR%U%NSIZE_NATURE,PATCH_NUMBER*NVAR,PATCH_NUMBER*NVAR))
 ALLOCATE(KRK(YSURF_CUR%U%NSIZE_NATURE,PATCH_NUMBER*NVAR,PATCH_NUMBER*NVAR))
 ALLOCATE(IdKH(YSURF_CUR%U%NSIZE_NATURE,PATCH_NUMBER*NVAR,PATCH_NUMBER*NVAR))
 ALLOCATE(FROZ(YSURF_CUR%U%NSIZE_NATURE,PATCH_NUMBER))
!
!   Initialisations
!
 HO(:,:,:)    = 0.0     ! Linearized observation matrix
 HOWR(:,:,:)  = 0.0   
 R(:,:,:)     = 0.0     ! Observation error matrix
 B(:,:,:)     = 0.0     ! Background error matrix
 YOWR(:,:)    = 0.0     ! Observation vector on the full grid to be written on file
 ZB(:,:)      = 0.0     ! Innovation vector
 YF(:,:,:)    = 0.0     ! Tile averaged simulated observation vector
 OBSCOUNT     = 0
!
 IDENT=0.0                  ! identity matrix
 DO L = 1,NVAR
   DO J=1,PATCH_NUMBER
     IDENT(J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(L-1)) = 1.0
   ENDDO
 ENDDO
!
!   Calculate delta for control variables 

   DO L=1,NVAR
     WRITE(LCHAR,'(I1) ') L
     NMFILE_CANARI='PERTURB_'//LCHAR//'_'
     CALL GET_FILE_NAME(IYEAR,IMONTH,IDAY,IHOUR,NMFILE_CANARI)
     PRINT *, 'read in PERTURBation for H: ', NMFILE_CANARI(1:LEN_TRIM(NMFILE_CANARI))//".DAT"
     OPEN(UNIT=111,FILE=TRIM(NMFILE_CANARI)//".DAT",FORM='FORMATTED',STATUS='OLD')
     DO JI=1,YSURF_CUR%U%NSIZE_NATURE
       DO J=1,PATCH_NUMBER
         READ(111,*) ZEPS(JI,J,L)
       ENDDO
     ENDDO
     CLOSE(111)
   ENDDO
!
 NOBS=0
 ZTIME=PTSTEP_OUTPUT
!
!   BEGINNING OF TIME LOOP
!
 TIMELOOP : DO ISTEP=1,NBOUTPUT
!
!   Update date
!
   CALL ADD_FORECAST_TO_DATE_SURF(IYEAR, IMONTH, IDAY, ZTIME)
   ZTIME = ZTIME + PTSTEP_OUTPUT
   NMFILE_CANARI='CANARI_NATURE_'
   CALL GET_FILE_NAME(IYEAR,IMONTH,IDAY,IHOUR,NMFILE_CANARI)
!
!   Open observation files
!
   OPEN(UNIT=ILOBS,FILE=TRIM(NMFILE_CANARI)//".DAT",FORM='FORMATTED',STATUS='OLD',ERR=10)
!
   NOBS=NOBS+NOBSTYPE 
!
!   If it exists, read observations
   DO JI=1,YSURF_CUR%U%NSIZE_NATURE
     READ (ILOBS,*)  (YO(JI,J),J=1,NOBSTYPE)
   ENDDO
   PRINT *,'read in obs: ', YO(1,:), NOBS
!
!   Calculate perturbations at the date of observations and innovation
!   (2nd dim of YF: index 1 is reference, 2,3... perturbed by cntrl var 1,2... 
!
   DO L=1,NVAR
     WRITE(LCHAR,'(I1) ') L
     NMFILE_CANARI='OBSIMU_PERT_'//LCHAR//'_'
     CALL GET_FILE_NAME(IYEAR,IMONTH,IDAY,IHOUR,NMFILE_CANARI)
     PRINT *, '--> read in ptbd forecasts for H: ', NMFILE_CANARI(1:LEN_TRIM(NMFILE_CANARI))//".DAT"
     OPEN(UNIT=111,FILE=TRIM(NMFILE_CANARI)//".DAT",FORM='FORMATTED',STATUS='OLD')
     DO JI=1,YSURF_CUR%U%NSIZE_NATURE
       DO J=1,PATCH_NUMBER
         READ(111,*) (YF_PATCH(JI,J,L+1,K+NOBS-NOBSTYPE),K=1,NOBSTYPE)
       ENDDO
     ENDDO
     CLOSE(111)
   ENDDO
!
   NMFILE_CANARI='OBSIMU_REFR_'
   CALL GET_FILE_NAME(IYEAR,IMONTH,IDAY,IHOUR,NMFILE_CANARI)
   PRINT *, '--> read in refr forecasts for H: ', NMFILE_CANARI(1:LEN_TRIM(NMFILE_CANARI))//".DAT"
   OPEN(UNIT=111,FILE=TRIM(NMFILE_CANARI)//".DAT",FORM='FORMATTED',STATUS='OLD')
!
   DO JI=1,YSURF_CUR%U%NSIZE_NATURE
     DO J=1,PATCH_NUMBER
       READ(111,*) (YF_PATCH(JI,J,1,K+NOBS-NOBSTYPE),K=1,NOBSTYPE)
     ENDDO
   ENDDO
   CLOSE(111)
   PRINT *, 'read in ref forecasts for H', YF_PATCH(1,10,1,:)
   10  CONTINUE ! if T2m HU2m observations does not exist
   CLOSE(ILOBS)
!
!  Mean simulated obs averaged over tiles
!
   DO JI=1,YSURF_CUR%U%NSIZE_NATURE
     DO J=1,PATCH_NUMBER
       DO K=1,NOBS
         DO L=1,NVAR+1
           IF (YF_PATCH(JI,J,L,K).NE.XUNDEF) THEN
             YF(JI,L,K) = YF(JI,L,K) + XPATCH(JI,J)*YF_PATCH(JI,J,L,K)
           ENDIF
         ENDDO
       ENDDO
     ENDDO
   ENDDO
 PRINT *, 'read in sim obs yf', YF(1,:,1)
 ENDDO TIMELOOP
!
! SET OBSERVATION ERROR 
!
!   Time reinitialization 
!
 IYEAR  = TTIME%TDATE%YEAR
 IMONTH = TTIME%TDATE%MONTH
 IDAY   = TTIME%TDATE%DAY
 ZTIME  = TTIME%TIME

 DO JI=1,YSURF_CUR%U%NSIZE_NATURE
   DO K=1,NOBS
     IF (XOBS(K) .EQ. 'LAI') THEN
            R(JI,K,K) = ERROBS(K)*ERROBS(K)*YO(JI,K)*YO(JI,K)
     ELSEIF (XOBS(K) .EQ. 'WG1') THEN
!  convert R for wg1 from SWI to abs value
       R(JI,K,K) = ERROBS(K)*ERROBS(K)*COFSWI(JI)*COFSWI(JI) 
     ELSE
       R(JI,K,K) = ERROBS(K)*ERROBS(K)
     ENDIF
   ENDDO
 ENDDO
!
! WRITE OUT OBS AND YERROR FOR DIAGNOSTIC PURPOSES
!
 OPEN (unit=111,file='OBSERRORout',status='unknown',IOSTAT=istat)
 WRITE (111,*) (R(JI,:,:), JI=1,YSURF_CUR%U%NSIZE_NATURE)
 CLOSE(111)
 OPEN (unit=111,file='OBSout',status='unknown',IOSTAT=istat)
 WRITE (111,*) (YO(JI,:), JI=1,YSURF_CUR%U%NSIZE_NATURE)
 CLOSE(111)
!--------------------------------------------------------------------
!
!   Background error matrix
!
!--------------------------------------------------------------------
 IF (LBFIXED) THEN 
   DO L=1,NVAR
     DO JI=1,YSURF_CUR%U%NSIZE_NATURE
       DO J=1,PATCH_NUMBER
         IF (XVAR(L) .EQ. 'WG1' .OR. XVAR(L) .EQ. 'WG2') THEN
           B(JI,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(L-1)) = XSIGMA(L)*XSIGMA(L)*COFSWI(JI)*COFSWI(JI)
         ELSEIF (XVAR(L) .EQ. 'LAI') THEN
           B(JI,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(L-1)) = XSIGMA(L)*XSIGMA(L)*LAI_PASS(JI,J)*LAI_PASS(JI,J)
         ELSE
           B(JI,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(L-1)) = XSIGMA(L)*XSIGMA(L)
         ENDIF
       ENDDO
     ENDDO
   ENDDO
 ELSE

ALLOCATE(Q(YSURF_CUR%U%NSIZE_NATURE,PATCH_NUMBER*NVAR,PATCH_NUMBER*NVAR))
   Q(:,:,:) = 0.0
   DO L=1,NVAR
     DO JI=1,YSURF_CUR%U%NSIZE_NATURE
       DO J=1,PATCH_NUMBER
           IF (XVAR(L) .EQ. 'WG2' .OR. XVAR(L) .EQ. 'WG1') THEN
             Q(JI,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(L-1)) = XSCALE_Q*XSCALE_Q*XSIGMA(L)*XSIGMA(L)*COFSWI(JI)*COFSWI(JI)
           ELSEIF (XVAR(L) .EQ. 'LAI') THEN
             Q(JI,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(L-1)) = XSCALE_QLAI*XSCALE_QLAI*XSIGMA(L)*XSIGMA(L)
          ELSE
             Q(JI,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(L-1)) = XSCALE_Q*XSCALE_Q*XSIGMA(L)*XSIGMA(L)
           ENDIF
       ENDDO
     ENDDO
   ENDDO
   print *
   PRINT *,'*** reading in B, then updating with Q ***'

 OPEN (unit=111,file='BGROUNDin',status='old',IOSTAT=istat)
 DO L=1,NVAR
   DO K=1,NVAR
     DO JI=1,YSURF_CUR%U%NSIZE_NATURE
       DO J=1,PATCH_NUMBER
         DO JJ=1,PATCH_NUMBER
            READ (111,*) B(JI,J+PATCH_NUMBER*(L-1),JJ+PATCH_NUMBER*(K-1))
         ENDDO
       ENDDO
     ENDDO
   ENDDO
 ENDDO
 CLOSE(111)
!
! B is the forecast matrix - need to add Q
!
   B = B + Q
!
   DEALLOCATE(Q)
 ENDIF! LBFIXED

!
!  check for frozen conditions
!
   CALL READ_SURF(&
                YPROGRAM,'WGI1', FROZ(:,:),IRESP)
   DO JI=1,YSURF_CUR%U%NSIZE_NATURE
     IF ( minval(FROZ(JI,:)) .GT. 0) THEN 
         YO(JI,:)=999.0
     ENDIF
   ENDDO
!
! Data type selection before assimilation (only if NOBS = NOBSTYPE)
!
 print *
 !PRINT *,'*** Calculating jacobians',NOBS,' with obs', YO(:,:), ' ***'
 !PRINT *,'*** Calculating Jacobians',NOBS, ' ***'
 print *
 DO K=1,NOBS
 DO L=1,NVAR
   DO JI=1,YSURF_CUR%U%NSIZE_NATURE
     DO J=1,PATCH_NUMBER  
           HOWR(JI,K,J+PATCH_NUMBER*(L-1)) = XPATCH(JI,J)*(YF_PATCH(JI,J,L+1,K) - YF_PATCH(JI,J,1,K))/ZEPS(JI,J,L) 
            IF(YO(JI,K) .NE. 999.0) THEN
              HO(JI,K,J+PATCH_NUMBER*(L-1)) = XPATCH(JI,J)*(YF_PATCH(JI,J,L+1,K) - YF_PATCH(JI,J,1,K))/ZEPS(JI,J,L)  ! Jacobian of obs operator
              HO(JI,K,J+PATCH_NUMBER*(L-1)) = max(-0.1,HO(JI,K,J+PATCH_NUMBER*(L-1)))
              HO(JI,K,J+PATCH_NUMBER*(L-1)) = min(1.0, HO(JI,K,J+PATCH_NUMBER*(L-1)))
              ZB(JI,K) = YO(JI,K) - YF(JI,1,K)                                                    ! innovation vector
              OBSCOUNT = OBSCOUNT + 1
            ELSE  
             HO(JI,K,J+PATCH_NUMBER*(L-1)) = 0.0
             ZB(JI,K) = 0.0 
            ENDIF
          ENDDO
     ENDDO
   ENDDO
 ENDDO
!
! *** Write innovations in ASCII file ***
!
 OPEN (unit=111,file='INNOV',status='unknown',IOSTAT=istat)
 DO JI=1,YSURF_CUR%U%NSIZE_NATURE
     WRITE(111,*) (ZB(JI,K),K=1,NOBS)
 ENDDO
 CLOSE(UNIT=111)

 OBSCOUNT = OBSCOUNT / PATCH_NUMBER / NVAR
!-----------------------------------------------------
!
!            ******  SOIL ANALYSIS *******
!
!-----------------------------------------------------
 PRINT *,'PERFORMING ANALYSIS'


 DO JI=1,YSURF_CUR%U%NSIZE_NATURE 
    HOT(JI,:,:) = TRANSPOSE(HO(JI,:,:)) 
    K1(JI,:,:) = MATMUL(HO(JI,:,:),MATMUL(B(JI,:,:),HOT(JI,:,:))) + R(JI,:,:)
    CALL CHOLDC(NOBS,K1(JI,:,:),ZP(JI,:))                                       ! Cholesky decomposition (1)
    CALL CHOLSL(NOBS,K1(JI,:,:),ZP(JI,:),ZB(JI,:),ZX(JI,:))                       ! Cholesky decomposition (2)
    XINCR(JI,:) = MATMUL(B(JI,:,:),MATMUL(HOT(JI,:,:),ZX(JI,:)))
    DO L=1,NVAR 
       IF (XVAR(L).EQ.'LAI') THEN
         DO J=1,PATCH_NUMBER
            XINCR(JI,J+PATCH_NUMBER*(L-1)) = MAX( XINCR(JI,J+PATCH_NUMBER*(L-1)), XVLAIMIN(J) - XI(JI,J,L) )
            XI(JI,J,L) = XI(JI,J,L) + XINCR(JI,J+PATCH_NUMBER*(L-1))
            ZBIO_PASS(JI,J) = ZBIO_PASS(JI,J) + XINCR(JI,J+PATCH_NUMBER*(L-1))*XALPH(J)
         ENDDO 
       ELSE
         DO J=1,PATCH_NUMBER           
            if (XI(JI,J,L)+ XINCR(JI,J+PATCH_NUMBER*(L-1)) .LT. 0) then
                 XINCR(JI,J+PATCH_NUMBER*(L-1))=0 
                !XINCR(I,J+PATCH_NUMBER*(L-1))= MAX( XINCR(I,J+PATCH_NUMBER*(L-1)), 0 - XI(I,J,L) )
            endif
                XI(JI,J,L) = XI(JI,J,L) + XINCR(JI,J+PATCH_NUMBER*(L-1))
          ENDDO
       ENDIF
     ENDDO
   ENDDO
!
! *** Write analysis results in ASCII file ***
!
 OPEN (unit=111,file='ANAL_INCR',status='unknown',IOSTAT=istat)
 DO JI=1,YSURF_CUR%U%NSIZE_NATURE
   DO J=1,PATCH_NUMBER
     WRITE(111,*) (XI(JI,J,L),L=1,NVAR), (XINCR(JI,J+PATCH_NUMBER*(L-1)),L=1,NVAR)
   ENDDO
 ENDDO
 CLOSE(UNIT=111)
!
!   Write analysis in PREP file
!
 CALL END_IO_SURF_n(YPROGRAM)
 CALL IO_BUFF_CLEAN
! CFILEOUT='PREP'
CALL INIT_IO_SURF_n(YSURF_CUR%DTCO, YSURF_CUR%DGU, YSURF_CUR%U, &
                        YPROGRAM,'NATURE','ISBA  ','WRITE')
 CALL IO_BUFF_CLEAN
 DO L = 1,NVAR
   CALL WRITE_SURF(YSURF_CUR%DGU, YSURF_CUR%U, &
                   YPROGRAM,XVAR(L),XI(:,:,L),IRESP,HCOMMENT=PREFIX(L))
   IF (XVAR(NVAR).EQ.'LAI') THEN
     CALL WRITE_SURF(YSURF_CUR%DGU, YSURF_CUR%U, &
                   YPROGRAM,CBIO,ZBIO_PASS(:,:),IRESP,HCOMMENT=CPREFIX_BIO)
   ENDIF
 ENDDO

 CALL END_IO_SURF_n(YPROGRAM)
!
! Analysis of B (for use in next cycle)
! Ba = (I-KH)Bf
! K = BfHT{K1}**-1
!
 DO JI=1,YSURF_CUR%U%NSIZE_NATURE
       !K1 = (R+H.B.HT) (calculate inverse -> output goes to K1)
       CALL INVERSE_MATRIX(NOBS,K1(JI,:,:),ZP(JI,:))
       GAIN(JI,:,:) = MATMUL(B(JI,:,:),MATMUL(HOT(JI,:,:),K1(JI,:,:)))
       KH(JI,:,:) = MATMUL(GAIN(JI,:,:),HO(JI,:,:))
       IdKH(JI,:,:) = IDENT - KH(JI,:,:) 
       KRK(JI,:,:) = MATMUL(GAIN(JI,:,:),MATMUL(R(JI,:,:),TRANSPOSE(GAIN(JI,:,:))))
       IF (.NOT.LBFIXED)  B(JI,:,:) = MATMUL(IdKH(JI,:,:),MATMUL(B(JI,:,:),TRANSPOSE(IdKH(JI,:,:))))+KRK(JI,:,:)
 ENDDO
!
! Write out analysed B (for use in next cycle)
!
 OPEN (unit=111,file='BGROUNDout',status='unknown',IOSTAT=istat)
 DO L=1,NVAR
   DO K=1,NVAR
     DO JI=1,YSURF_CUR%U%NSIZE_NATURE
       DO J=1,PATCH_NUMBER
         DO JJ=1,PATCH_NUMBER
           WRITE (111,*) B(JI,J+PATCH_NUMBER*(L-1),JJ+PATCH_NUMBER*(K-1))
         ENDDO
       ENDDO
     ENDDO
   ENDDO
 ENDDO
 CLOSE(111)
!
! **** Write out the observation operator + Gain matrix ****
!
 DO L=1,NVAR
   DO K=1, NOBSTYPE
     WRITE(LCHAR,'(I1)') K
     HFNAME='HO_'//XVAR(L)//'_v'//LCHAR
     OPEN(UNIT=111,FILE=HFNAME,FORM='FORMATTED',STATUS='NEW',IOSTAT=istat)
     DO JI=1,YSURF_CUR%U%NSIZE_NATURE
       DO J=1,PATCH_NUMBER
         WRITE(111,*) HOWR(JI,K,J+PATCH_NUMBER*(L-1)),GAIN(JI,J+PATCH_NUMBER*(L-1),K)
       ENDDO
     ENDDO
     CLOSE(111)
   ENDDO
 ENDDO

 DEALLOCATE(YOWR)
 DEALLOCATE(VECT)
 DEALLOCATE(YF)
 DEALLOCATE(YF_PATCH)
 DEALLOCATE(ZB)
 DEALLOCATE(HO)
 DEALLOCATE(R)
 DEALLOCATE(B)
 DEALLOCATE(GAIN)
 DEALLOCATE(ZP)
 DEALLOCATE(ZEPS)
 DEALLOCATE(YO)
 DEALLOCATE(HOT)
 DEALLOCATE(K1)
 DEALLOCATE(ZX)
 DEALLOCATE(XINCR)
 DEALLOCATE(HOWR)
 DEALLOCATE(LAI_PASS)
 DEALLOCATE(ZBIO_PASS)
 DEALLOCATE(ZBIO_OUT)
 DEALLOCATE(KH)
 DEALLOCATE(KRK)
 DEALLOCATE(IdKH)
 DEALLOCATE(XPATCH)
 DEALLOCATE(LPATCH)
 DEALLOCATE(COFSWI)
 DEALLOCATE(WILT)
 DEALLOCATE(FROZ)

DEALLOCATE(NWORK)
DEALLOCATE(NWORK_FULL)
DEALLOCATE(XWORK)
DEALLOCATE(XWORK_FULL)
DEALLOCATE(XWORK2)
DEALLOCATE(XWORK2_FULL)
DEALLOCATE(XWORK3)
!
 PRINT *
 PRINT *,'   ---------------------------------------'
 PRINT *,'   |   EXITING VARASSIM AFTER ANALYSIS   |'
 PRINT *,'   ---------------------------------------'
 PRINT *
 PRINT *,'Number of assimilated observations =',OBSCOUNT
 PRINT *
!
 CALL SURFEX_DEALLO_LIST
!
CALL END_LOG_MPI

IF (LHOOK) CALL DR_HOOK('VARASSIM',1,ZHOOK_HANDLE)
!
#ifdef SFX_MPI
 CALL MPI_FINALIZE(INFOMPI)
#endif
!-------------------------------------------------------------------------------
!
END PROGRAM VARASSIM
