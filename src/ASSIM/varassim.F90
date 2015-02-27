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
!		- added Biomass as control and LAI as observation variable

! Modifications to allow for SWI assimilation by A.Barbu (25 Sept 2009)
! The new state vector can be any element of (WG2, WG1, LAI, Biomass) - Choice in namelist
! The new observation vector can be any element of (WG1, LAI, SWI, FAPAR) - Choice in namelist
! Bug correction regarding the Kalman gain computation (otherwise underestimation of the gain)
! New formulation of the B matrix to ensure its positive definiteness
! -----------------------------------------------------------------------------
!
USE MODD_SURFEX_MPI, ONLY : NCOMM, NPROC, NRANK, NPIO, WLOG_MPI, PREP_LOG_MPI,   &
                            NINDEX, NSIZE_TASK, END_LOG_MPI
USE MODD_SURFEX_OMP, ONLY : NINDX2SFX, NWORK, XWORK, XWORK2, XWORK3, NBLOCK, NBLOCKTOT,&
                            NWORK_FULL, XWORK_FULL, XWORK2_FULL, INIT_DIM, RESET_DIM
USE MODD_MASK, ONLY: NMASK_FULL
                          
 USE MODD_TYPE_DATE_SURF
 USE MODD_SURF_PAR,ONLY : XUNDEF
 USE MODD_SURF_ATM_n,      ONLY :  CSEA,        CWATER,      CTOWN,      CNATURE,      &
                                   XSEA,        XWATER,      XTOWN,      XNATURE,      &
                                   NSIZE_SEA,   NSIZE_WATER, NSIZE_TOWN, NSIZE_NATURE, &
                                   XCOVER,      NDIM_FULL,   NSIZE_FULL, LCOVER      , &
                                   NDIM_NATURE, NDIM_SEA,    NDIM_WATER, NDIM_TOWN
 USE MODD_SURF_CONF, ONLY : CPROGNAME
 USE MODD_OL_FILEID
 USE MODD_IO_SURF_ASC,ONLY : CFILEIN, CFILEIN_SAVE, CFILEPGD, CFILEOUT 
!
 USE MODI_ALLOC_SURFEX
 USE MODI_OPEN_NAMELIST
 USE MODI_CLOSE_NAMELIST
 USE MODE_POS_SURF 
 USE MODI_GOTO_SURFEX
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
 USE MODI_IO_BUFF_CLEAN_n
 USE MODI_WRITE_SURF
 USE MODI_ADD_FORECAST_TO_DATE_SURF
 USE MODI_CHOLDC
 USE MODI_CHOLSL
 USE MODI_INVERSE_MATRIX
 USE MODI_DEALLOC_SURFEX
 USE MODI_INIT_INDEX_MPI
!
 USE MODE_EKF, ONLY : GET_FILE_NAME
!
 USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
 USE PARKIND1  ,ONLY : JPRB
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
 INTEGER                                    :: I,J,JJ,K,KK,L
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
 CALL ALLOC_SURFEX(1)
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
 CALL GOTO_SURFEX(1,.TRUE.)
 CALL INIT_INDEX_MPI(CPROGNAME, 'LIN ', 0., .FALSE.)
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
 DO I = 1,NVARMAX
  IF (NNCV(I) == 1 .AND. J <= NVAR ) THEN
   TPRT(J) = XTPRT_M(I)
   XSIGMA(J) = XSIGMA_M(I)
   XVAR(J) = CVAR_M(I)
   PREFIX(J) = CPREFIX_M(I)
   J = J + 1
  ENDIF
 ENDDO
!
!      0.5 Allocate arrays depending on the number of observation type
!
 ALLOCATE (XOBS(NOBSTYPE))
 ALLOCATE (ERROBS(NOBSTYPE))
!
!	assigning of observation variable names farther down
!
!      1.    Initializations
!
 ILOBS=55
!
 CALL INI_CSTS
!
 CALL INI_DATA_COVER

 CFILEIN = 'PREP.txt'          
 CFILEPGD = 'PGD.txt'
 CFILEIN_SAVE = CFILEIN
!
!   Read grid dimension for allocation
!
 CALL INIT_IO_SURF_n(YPROGRAM,'FULL  ','SURF  ','READ ')
!
!   Find current time
!
 CALL READ_SURF(YPROGRAM,'DTCUR',TTIME,IRESP)
print *, 'Time from PREP file : ', TTIME 
!
 CALL END_IO_SURF_n(YPROGRAM)
 CALL SET_SURFEX_FILEIN(YPROGRAM,'PGD ') ! change input file name to pgd name
 CALL INIT_IO_SURF_n(YPROGRAM,'FULL  ','SURF  ','READ ')
!
!   Reading grid characteristics to perform nature mask
!
 CALL READ_SURF(YPROGRAM,'SEA   ',CSEA   ,IRESP)
 CALL READ_SURF(YPROGRAM,'WATER ',CWATER ,IRESP)
 CALL READ_SURF(YPROGRAM,'NATURE',CNATURE,IRESP)
 CALL READ_SURF(YPROGRAM,'TOWN  ',CTOWN  ,IRESP)

!
 CALL READ_SURF(YPROGRAM,'DIM_FULL  ',NDIM_FULL,  IRESP)
 CALL READ_SURF(YPROGRAM,'DIM_SEA   ',NDIM_SEA,   IRESP)
 CALL READ_SURF(YPROGRAM,'DIM_NATURE',NDIM_NATURE,IRESP)
 CALL READ_SURF(YPROGRAM,'DIM_WATER ',NDIM_WATER, IRESP)
 CALL READ_SURF(YPROGRAM,'DIM_TOWN  ',NDIM_TOWN,  IRESP)

NINDX2SFX = NDIM_FULL
ALLOCATE(NWORK(NDIM_FULL))
ALLOCATE(XWORK(NDIM_FULL))
ALLOCATE(XWORK2(NDIM_FULL,10))
ALLOCATE(XWORK3(NDIM_FULL,10,10))
IF (NRANK==NPIO) THEN
  ALLOCATE(NWORK_FULL(NDIM_FULL))
  ALLOCATE(XWORK_FULL(NDIM_FULL))
  ALLOCATE(XWORK2_FULL(NDIM_FULL,10))
ELSE
  ALLOCATE(NWORK_FULL(0))
  ALLOCATE(XWORK_FULL(0))
  ALLOCATE(XWORK2_FULL(0,0))
ENDIF
!
!   Get number of points on this proc
!
 CALL GET_SIZE_FULL_n(YPROGRAM,NDIM_FULL,NSIZE_FULL)
 CALL READ_COVER_n(YPROGRAM)
 CALL END_IO_SURF_n(YPROGRAM)
!
!   Perform masks (only nature used)
!
 ALLOCATE(XSEA(NDIM_FULL))
 ALLOCATE(XNATURE(NDIM_FULL))
 ALLOCATE(XWATER(NDIM_FULL))
 ALLOCATE(XTOWN(NDIM_FULL))
!
 CALL CONVERT_COVER_FRAC(XCOVER,LCOVER,XSEA,XNATURE,XTOWN,XWATER)
!
 NSIZE_NATURE = COUNT(XNATURE(:) > 0.0)
 NSIZE_TOWN   = COUNT(XTOWN(:)   > 0.0)
 NSIZE_WATER  = COUNT(XWATER(:)  > 0.0)
 NSIZE_SEA    = COUNT(XSEA(:)    > 0.0)
!
!   Read number of patches
!
 CALL SET_SURFEX_FILEIN(YPROGRAM,'PGD ')
 CALL INIT_IO_SURF_n(YPROGRAM,'NATURE','ISBA  ','READ ')
 CALL READ_SURF(YPROGRAM,'PATCH_NUMBER',PATCH_NUMBER,IRESP)
!
 ALLOCATE(ZSAND(NSIZE_NATURE))
 ALLOCATE(ZCLAY(NSIZE_NATURE))
!
!   Read CLAY fraction to compute the SWI range (Wfc - Wwilt)
!   (XSIGMA is defined in terms of SWI), need to convert to equivalent v/v
!   using same clay fraction in both layers
!   Read SAND fraction to compute the saturation for conversion of ERS SWI
!
 CALL READ_SURF(YPROGRAM,'CLAY',ZCLAY(:),IRESP)
 CALL READ_SURF(YPROGRAM,'SAND',ZSAND(:),IRESP)
!
!   Define prefixes for simulated observations
! 
 CALL END_IO_SURF_n(YPROGRAM)
 CALL SET_SURFEX_FILEIN(YPROGRAM,'PREP') ! change input file name to pgd name
 CALL INIT_IO_SURF_n(YPROGRAM,'NATURE','ISBA  ','READ ')
!
 XOBS_M(1) = 'WG1'
 XOBS_M(2) = 'LAI'
!
 J = 1
 DO I = 1,NOBSMAX
  IF (NNCO(I) == 1 .AND. J <= NOBSTYPE ) THEN
   ERROBS(J) = XERROBS_M(I)
   XOBS(J) = XOBS_M(I)
   J = J + 1
  ENDIF
 ENDDO
!
 ALLOCATE(XI(NSIZE_NATURE,PATCH_NUMBER,NVAR))
 ALLOCATE(VECT(NSIZE_NATURE,PATCH_NUMBER))
 ALLOCATE(LAI_PASS(NSIZE_NATURE,PATCH_NUMBER))
 ALLOCATE(ZBIO_PASS(NSIZE_NATURE,PATCH_NUMBER))
 ALLOCATE(ZBIO_OUT(NSIZE_NATURE,PATCH_NUMBER))
!
 ! Avoiding unsuitable LAI values
 ALLOCATE(XVLAIMIN(PATCH_NUMBER))
 XVLAIMIN(:) = 0.3
 XVLAIMIN(5:6) = 1.0
!
! Read in control variables
!
 DO L = 1,NVAR
   CALL READ_SURF(YPROGRAM,XVAR(L),XI(:,:,L),IRESP)
   IF (XVAR(L) .EQ. 'LAI') THEN
      LAI_PASS=XI(:,:,L)
      IF (CBIO .NE. 'LAI') THEN
        CALL READ_SURF(YPROGRAM,CBIO,ZBIO_PASS(:,:),IRESP)        
      ELSE
        ZBIO_PASS(:,:) = XI(:,:,L)
      ENDIF
   ENDIF
 ENDDO
 WHERE (LAI_PASS(:,:)==XUNDEF) LAI_PASS(:,:)=0.
!
ALLOCATE(XPATCH(NSIZE_NATURE,PATCH_NUMBER))
ALLOCATE(LPATCH(NSIZE_NATURE,PATCH_NUMBER))
!
! Read fraction of each patch (=> LPGD=T in OPTIONS.nam)
!
 CALL END_IO_SURF_n(YPROGRAM)
 CALL IO_BUFF_CLEAN_n
 CALL SET_SURFEX_FILEIN(YPROGRAM,'PREP')
 CALL INIT_IO_SURF_n(YPROGRAM,'NATURE','ISBA  ','READ ')
 IF (PATCH_NUMBER > 1) THEN
   CALL READ_SURF(YPROGRAM,'PATCH',XPATCH,IRESP)
 ELSE
   XPATCH(:,1) = 1.0 
 ENDIF
!
LPATCH(:,:) = .FALSE.
!
   DO I=1,NSIZE_NATURE
     DO J=1,PATCH_NUMBER
         IF (XPATCH(I,J) > 0.0) THEN
           LPATCH(I,J) = .TRUE.
         ENDIF
     ENDDO
   ENDDO
!
 ALLOCATE(COFSWI(NSIZE_NATURE))
 ALLOCATE(WILT(NSIZE_NATURE))
 ALLOCATE(SMSAT(NSIZE_NATURE)) 
 ALLOCATE(B(NSIZE_NATURE,PATCH_NUMBER*NVAR,PATCH_NUMBER*NVAR))
!
 DO I=1,NSIZE_NATURE
  COFSWI(I)=0.001*(89.0467*((100.*ZCLAY(I))**0.3496)-37.1342*((100.*ZCLAY(I))**0.5))
  SMSAT(I)=0.001*(-1.08*100*ZSAND(I)+494.305) 
  WILT(I)=0.001*37.1342*((100.*ZCLAY(I))**0.5)   
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
   DO I=1,NSIZE_NATURE
     DO J=1,PATCH_NUMBER
       IF (XI(I,J,NIVAR).NE.XUNDEF) THEN                                 ! check whether values are undefined
          VECT(I,J) = XI(I,J,NIVAR) + TPRT(NIVAR)*XI(I,J,NIVAR) 
       ELSE                                                             
          VECT(I,J) = XI(I,J,NIVAR)
       ENDIF                 
       IF (XI(I,J,NIVAR).NE.XUNDEF .AND. XVAR(NIVAR).EQ.'LAI') THEN
          ZBIO_OUT(I,J) = ZBIO_PASS(I,J) + TPRT(NIVAR)*ZBIO_PASS(I,J)
       ELSE                                                     
          ZBIO_OUT(I,J) = ZBIO_PASS(I,J)
       ENDIF      
     ENDDO
   ENDDO

   CALL END_IO_SURF_n(YPROGRAM)
   CALL IO_BUFF_CLEAN_n  
   CALL INIT_IO_SURF_n(YPROGRAM,'NATURE','ISBA  ','WRITE')
   CALL WRITE_SURF(YPROGRAM,XVAR(NIVAR),VECT,IRESP,HCOMMENT=PREFIX(NIVAR))
   IF (XVAR(NIVAR).EQ.'LAI') THEN
     CALL WRITE_SURF(YPROGRAM,CBIO,ZBIO_OUT,IRESP,HCOMMENT=CPREFIX_BIO)
   ENDIF  
   CALL END_IO_SURF_n(YPROGRAM)
   CALL IO_BUFF_CLEAN_n
!
! write out perturbation
   OPEN (unit=111,file='PERTURB',status='unknown',IOSTAT=istat)
   DO I=1,NSIZE_NATURE
     DO J=1,PATCH_NUMBER
       IF (XI(I,J,NIVAR).NE.XUNDEF) THEN
         WRITE (111,*) TPRT(NIVAR)*XI(I,J,NIVAR)
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
     DO I=1,NSIZE_NATURE
       DO J=1,PATCH_NUMBER
         IF (XVAR(L) .EQ. 'WG1' .OR. XVAR(L) .EQ. 'WG2') THEN
          B(I,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(L-1)) = XSIGMA(L)*XSIGMA(L)*COFSWI(I)*COFSWI(I)
         ELSEIF (XVAR(L) .EQ. 'LAI') THEN 
          B(I,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(L-1)) = XSIGMA(L)*XSIGMA(L)!*LAI_PASS(I,J)*LAI_PASS(I,J)
         ELSE								
           B(I,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(L-1)) = XSIGMA(L)*XSIGMA(L)                    		     
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
     DO I=1,NSIZE_NATURE
       DO J=1,PATCH_NUMBER
         DO JJ=1,PATCH_NUMBER
           DO L=1,NVAR
             DO K=1,NVAR
               WRITE (111,*) B(I,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(L-1))
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
   CALL IO_BUFF_CLEAN_n
   CALL SET_SURFEX_FILEIN(YPROGRAM,'PREP')    
   CALL INIT_IO_SURF_n(YPROGRAM,'NATURE','ISBA  ','READ ')
!
   ALLOCATE(SIMOBS(NSIZE_NATURE,PATCH_NUMBER,NOBSTYPE))
   ALLOCATE(SIMMOD(NSIZE_NATURE,PATCH_NUMBER,NVAR))
!
   !read observation from file
   DO K = 1,NOBSTYPE
     CALL READ_SURF(YPROGRAM,XOBS(K),SIMOBS(:,:,K),IRESP) 
     IF( XOBS(K) == 'SWI2' ) THEN
       !read file containing wg2 minmax values -> use for normalization
       ALLOCATE(MWG2(NSIZE_NATURE,PATCH_NUMBER))
       ALLOCATE(XWG2(NSIZE_NATURE,PATCH_NUMBER))
       OPEN(UNIT=111,FILE='mxwg2.dat',FORM='FORMATTED',STATUS='OLD',IOSTAT=istat)
       IF (ISTAT.EQ.0) THEN
         DO I=1,NSIZE_NATURE
           DO J=1,PATCH_NUMBER
             READ(111,*) MWG2(I,J),XWG2(I,J)
             WG2MIN = MWG2(I,J)
             SCALE_SWI = XWG2(I,J) - MWG2(I,J)
             IF( LPATCH(I,J).AND. SCALE_SWI > 1.0E-6 ) THEN
               SIMOBS(I,J,K)=(SIMOBS(I,J,K)*COFSWI(I)+WILT(I)-WG2MIN)/SCALE_SWI
             ENDIF
           ENDDO
         ENDDO
         CLOSE(111)
       ELSE
         PRINT *, 'No mxwg2.dat file, use namelist given WG2min and scale values'
         DO I=1,NSIZE_NATURE
           DO J=1,PATCH_NUMBER
             SIMOBS(I,J,K)=(SIMOBS(I,J,K)*COFSWI(I)+WILT(I)-WG2MIN)/SCALE_SWI
           ENDDO
         ENDDO
       ENDIF
       DEALLOCATE(MWG2)
       DEALLOCATE(XWG2)
     ENDIF
  ENDDO
                                                                                                                                                                                                                                                                                       
   CALL END_IO_SURF_n(YPROGRAM)
   CALL IO_BUFF_CLEAN_n
   CALL SET_SURFEX_FILEIN(YPROGRAM,'PREP')
   CALL INIT_IO_SURF_n(YPROGRAM,'NATURE','ISBA  ','READ ')
!
   OPEN (unit=111,file='OBSIMU',status='unknown')
   DO I=1,NSIZE_NATURE
     DO J=1,PATCH_NUMBER
       WRITE (111,*) (SIMOBS(I,J,K),K = 1,NOBSTYPE)
     ENDDO
   ENDDO
   CLOSE(111)
!
   DO I = 1, NVAR
     CALL READ_SURF(YPROGRAM,XVAR(I),SIMMOD(:,:,I),IRESP)
   ENDDO
   OPEN (unit=111,file='MDSIMU',status='unknown',IOSTAT=istat)
   DO I=1,NSIZE_NATURE
     DO J=1,PATCH_NUMBER
       WRITE (111,*) (SIMMOD(I,J,L),L=1,NVAR)
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

   ALLOCATE(LTM(NSIZE_NATURE,PATCH_NUMBER*NVAR,PATCH_NUMBER*NVAR))
   ALLOCATE(ZEPS(NSIZE_NATURE,PATCH_NUMBER,NVAR))
   ALLOCATE(XF(NSIZE_NATURE,PATCH_NUMBER,NVAR+1,NVAR))

!
   PRINT *, 'evolving B to time ', IYEAR,IMONTH,IDAY,IHOUR  
!
 !read B matrix from previous analysis step: B(t-1)
  OPEN (unit=111,file='BGROUNDin',status='old',IOSTAT=istat)
  DO L=1,NVAR   ! control variable (x at previous time step)
     DO K=1,NVAR
       DO I=1,NSIZE_NATURE
         DO J=1,PATCH_NUMBER   
           DO JJ=1,PATCH_NUMBER   
             READ (111,*) B(I,J+PATCH_NUMBER*(L-1),JJ+PATCH_NUMBER*(K-1))
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
     DO I=1,NSIZE_NATURE
       DO J=1,PATCH_NUMBER 
         READ(111,*) (XF(I,J,L+1,K),K=1,NVAR)
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
   DO I=1,NSIZE_NATURE
     DO J=1,PATCH_NUMBER
       READ(111,*) (XF(I,J,1,K),K=1,NVAR)
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
     DO I=1,NSIZE_NATURE
       DO J=1,PATCH_NUMBER
         READ(111,*) ZEPS(I,J,L)
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
    DO I=1,NSIZE_NATURE 
      DO J=1,PATCH_NUMBER 
           IF (LPATCH(I,J) .AND. XF(I,J,L+1,K).NE.XUNDEF .AND. XF(I,J,1,K).NE.XUNDEF ) THEN
             ! Jacobian of fwd model
             LTM(I,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(K-1)) = (XF(I,J,L+1,K) - XF(I,J,1,K))/ZEPS(I,J,L)
             LTM(I,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(K-1)) = max(-0.1, LTM(I,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(K-1)))! impose upper/lower limits 
             LTM(I,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(K-1)) = min(1.0, LTM(I,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(K-1)))
           ENDIF
         ENDDO
       ENDDO
     ENDDO
   ENDDO
!
! e) evolve B 
!
   DO I=1,NSIZE_NATURE
     B(I,:,:)=MATMUL(LTM(I,:,:),MATMUL(B(I,:,:),TRANSPOSE(LTM(I,:,:))))     
   ENDDO
!
! write out the LTM for the forward model
!
   DO L=1,NVAR
     DO K=1,NVAR 
       LFNAME='LTM_del'//XVAR(K)//'_del'//XVAR(L)
       OPEN(UNIT=111,FILE=LFNAME,FORM='FORMATTED',STATUS='UNKNOWN',ACCESS='APPEND')
       DO I=1,NSIZE_NATURE
         DO J=1,PATCH_NUMBER
           WRITE (111,*) LTM(I,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(K-1))
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
      DO I=1,NSIZE_NATURE
        DO J=1,PATCH_NUMBER
           DO JJ=1,PATCH_NUMBER
             WRITE (111,*)  B(I,J+PATCH_NUMBER*(L-1),JJ+PATCH_NUMBER*(K-1))
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
 ALLOCATE(YF(NSIZE_NATURE,NVAR+1,NOBS))
 ALLOCATE(YF_PATCH(NSIZE_NATURE,PATCH_NUMBER,NVAR+1,NOBS))
 ALLOCATE(ZEPS(NSIZE_NATURE,PATCH_NUMBER,NVAR))
 ALLOCATE(XF(NSIZE_NATURE,PATCH_NUMBER,NVAR+1,NVAR))
!
! Read fraction of each patch (=> LPGD=T in OPTIONS.nam)
!
 CALL END_IO_SURF_n(YPROGRAM)
 CALL IO_BUFF_CLEAN_n
 CALL SET_SURFEX_FILEIN(YPROGRAM,'PREP')
 CALL INIT_IO_SURF_n(YPROGRAM,'NATURE','ISBA  ','READ ')
 IF (PATCH_NUMBER > 1) THEN
   CALL READ_SURF(YPROGRAM,'PATCH',XPATCH,IRESP)
 ELSE
   XPATCH(:,1) = 1.0 
 ENDIF
!
!   Initial values (to be analysed)
!   Observations
!
 ALLOCATE(YO(NSIZE_NATURE,NOBS))
 ALLOCATE(YOWR(NDIM_FULL,NOBS))
!
!   Temporary vectors used by the EKF approach
!
 ALLOCATE(HO(NSIZE_NATURE,NOBS,PATCH_NUMBER*NVAR))
 ALLOCATE(HOWR(NSIZE_NATURE,NOBS,PATCH_NUMBER*NVAR))
 ALLOCATE(HOT(NSIZE_NATURE,PATCH_NUMBER*NVAR,NOBS))
 ALLOCATE(R(NSIZE_NATURE,NOBS,NOBS))
 ALLOCATE(IDENT(PATCH_NUMBER*NVAR,PATCH_NUMBER*NVAR))
 ALLOCATE(GAIN(NSIZE_NATURE,PATCH_NUMBER*NVAR,NOBS))
 ALLOCATE(ZB(NSIZE_NATURE,NOBS))
 ALLOCATE(ZX(NSIZE_NATURE,NOBS))
 ALLOCATE(ZP(NSIZE_NATURE,NOBS))
 ALLOCATE(XINCR(NSIZE_NATURE,PATCH_NUMBER*NVAR))
 ALLOCATE(K1(NSIZE_NATURE,NOBS,NOBS))
 ALLOCATE(KH(NSIZE_NATURE,PATCH_NUMBER*NVAR,PATCH_NUMBER*NVAR))
 ALLOCATE(KRK(NSIZE_NATURE,PATCH_NUMBER*NVAR,PATCH_NUMBER*NVAR))
 ALLOCATE(IdKH(NSIZE_NATURE,PATCH_NUMBER*NVAR,PATCH_NUMBER*NVAR))
 ALLOCATE(FROZ(NSIZE_NATURE,PATCH_NUMBER))
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
     DO I=1,NSIZE_NATURE
       DO J=1,PATCH_NUMBER
         READ(111,*) ZEPS(I,J,L)
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
   DO I=1,NSIZE_NATURE
     READ (ILOBS,*)  (YO(I,J),J=1,NOBSTYPE)
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
     DO I=1,NSIZE_NATURE
       DO J=1,PATCH_NUMBER
         READ(111,*) (YF_PATCH(I,J,L+1,K+NOBS-NOBSTYPE),K=1,NOBSTYPE)
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
   DO I=1,NSIZE_NATURE
     DO J=1,PATCH_NUMBER
       READ(111,*) (YF_PATCH(I,J,1,K+NOBS-NOBSTYPE),K=1,NOBSTYPE)
     ENDDO
   ENDDO
   CLOSE(111)
   PRINT *, 'read in ref forecasts for H', YF_PATCH(1,10,1,:)				
   10  CONTINUE ! if T2m HU2m observations does not exist
   CLOSE(ILOBS)
!
!  Mean simulated obs averaged over tiles
!
   DO I=1,NSIZE_NATURE
     DO J=1,PATCH_NUMBER
       DO K=1,NOBS
         DO L=1,NVAR+1
           IF (YF_PATCH(I,J,L,K).NE.XUNDEF) THEN
             YF(I,L,K) = YF(I,L,K) + XPATCH(I,J)*YF_PATCH(I,J,L,K)
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

 DO I=1,NSIZE_NATURE
   DO K=1,NOBS
     IF (XOBS(K) .EQ. 'LAI') THEN
            R(I,K,K) = ERROBS(K)*ERROBS(K)*YO(I,K)*YO(I,K)                                                  	                                                        
     ELSEIF (XOBS(K) .EQ. 'WG1') THEN					        
!  convert R for wg1 from SWI to abs value
       R(I,K,K) = ERROBS(K)*ERROBS(K)*COFSWI(I)*COFSWI(I) 
     ELSE								         
       R(I,K,K) = ERROBS(K)*ERROBS(K)	                                
     ENDIF
   ENDDO
 ENDDO
!
! WRITE OUT OBS AND YERROR FOR DIAGNOSTIC PURPOSES
!
 OPEN (unit=111,file='OBSERRORout',status='unknown',IOSTAT=istat)
 WRITE (111,*) (R(I,:,:), I=1,NSIZE_NATURE)
 CLOSE(111)
 OPEN (unit=111,file='OBSout',status='unknown',IOSTAT=istat)
 WRITE (111,*) (YO(I,:), I=1,NSIZE_NATURE)
 CLOSE(111)
!--------------------------------------------------------------------
!
!   Background error matrix
!
!--------------------------------------------------------------------
 IF (LBFIXED) THEN 
   DO L=1,NVAR
     DO I=1,NSIZE_NATURE
       DO J=1,PATCH_NUMBER
         IF (XVAR(L) .EQ. 'WG1' .OR. XVAR(L) .EQ. 'WG2') THEN		
           B(I,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(L-1)) = XSIGMA(L)*XSIGMA(L)*COFSWI(I)*COFSWI(I)
         ELSEIF (XVAR(L) .EQ. 'LAI') THEN                                  	
           B(I,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(L-1)) = XSIGMA(L)*XSIGMA(L)*LAI_PASS(I,J)*LAI_PASS(I,J)
         ELSE								
           B(I,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(L-1)) = XSIGMA(L)*XSIGMA(L)                    		     
         ENDIF
       ENDDO
     ENDDO
   ENDDO
 ELSE

ALLOCATE(Q(NSIZE_NATURE,PATCH_NUMBER*NVAR,PATCH_NUMBER*NVAR))
   Q(:,:,:) = 0.0
   DO L=1,NVAR
     DO I=1,NSIZE_NATURE
       DO J=1,PATCH_NUMBER
           IF (XVAR(L) .EQ. 'WG2' .OR. XVAR(L) .EQ. 'WG1') THEN                     
             Q(I,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(L-1)) = XSCALE_Q*XSCALE_Q*XSIGMA(L)*XSIGMA(L)*COFSWI(I)*COFSWI(I)
           ELSEIF (XVAR(L) .EQ. 'LAI') THEN
             Q(I,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(L-1)) = XSCALE_QLAI*XSCALE_QLAI*XSIGMA(L)*XSIGMA(L)
          ELSE
             Q(I,J+PATCH_NUMBER*(L-1),J+PATCH_NUMBER*(L-1)) = XSCALE_Q*XSCALE_Q*XSIGMA(L)*XSIGMA(L)
           ENDIF
       ENDDO
     ENDDO
   ENDDO
   print *
   PRINT *,'*** reading in B, then updating with Q ***'

 OPEN (unit=111,file='BGROUNDin',status='old',IOSTAT=istat)
 DO L=1,NVAR
   DO K=1,NVAR
     DO I=1,NSIZE_NATURE
       DO J=1,PATCH_NUMBER
         DO JJ=1,PATCH_NUMBER
            READ (111,*) B(I,J+PATCH_NUMBER*(L-1),JJ+PATCH_NUMBER*(K-1))
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
   CALL READ_SURF(YPROGRAM,'WGI1', FROZ(:,:),IRESP)
   DO I=1,NSIZE_NATURE     
     IF ( minval(FROZ(I,:)) .GT. 0) THEN 
         YO(I,:)=999.0
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
   DO I=1,NSIZE_NATURE
     DO J=1,PATCH_NUMBER  
           HOWR(I,K,J+PATCH_NUMBER*(L-1)) = XPATCH(I,J)*(YF_PATCH(I,J,L+1,K) - YF_PATCH(I,J,1,K))/ZEPS(I,J,L) 
            IF(YO(I,K) .NE. 999.0) THEN         
              HO(I,K,J+PATCH_NUMBER*(L-1)) = XPATCH(I,J)*(YF_PATCH(I,J,L+1,K) - YF_PATCH(I,J,1,K))/ZEPS(I,J,L)  ! Jacobian of obs operator
              HO(I,K,J+PATCH_NUMBER*(L-1)) = max(-0.1,HO(I,K,J+PATCH_NUMBER*(L-1)))             
              HO(I,K,J+PATCH_NUMBER*(L-1)) = min(1.0, HO(I,K,J+PATCH_NUMBER*(L-1)))
              ZB(I,K) = YO(I,K) - YF(I,1,K)                                                    ! innovation vector
              OBSCOUNT = OBSCOUNT + 1
            ELSE  
             HO(I,K,J+PATCH_NUMBER*(L-1)) = 0.0
             ZB(I,K) = 0.0 
            ENDIF
          ENDDO
     ENDDO
   ENDDO
 ENDDO
!
! *** Write innovations in ASCII file ***
!
 OPEN (unit=111,file='INNOV',status='unknown',IOSTAT=istat)
 DO I=1,NSIZE_NATURE
     WRITE(111,*) (ZB(I,K),K=1,NOBS)
 ENDDO
 CLOSE(UNIT=111)

 OBSCOUNT = OBSCOUNT / PATCH_NUMBER / NVAR        
!-----------------------------------------------------
!
!            ******  SOIL ANALYSIS *******
!
!-----------------------------------------------------
 PRINT *,'PERFORMING ANALYSIS'


 DO I=1,NSIZE_NATURE 
    HOT(I,:,:) = TRANSPOSE(HO(I,:,:)) 
    K1(I,:,:) = MATMUL(HO(I,:,:),MATMUL(B(I,:,:),HOT(I,:,:))) + R(I,:,:)
    CALL CHOLDC(NOBS,K1(I,:,:),ZP(I,:))                                       ! Cholesky decomposition (1)
    CALL CHOLSL(NOBS,K1(I,:,:),ZP(I,:),ZB(I,:),ZX(I,:))                       ! Cholesky decomposition (2)     
    XINCR(I,:) = MATMUL(B(I,:,:),MATMUL(HOT(I,:,:),ZX(I,:)))
    DO L=1,NVAR 
       IF (XVAR(L).EQ.'LAI') THEN
         DO J=1,PATCH_NUMBER
            XINCR(I,J+PATCH_NUMBER*(L-1)) = MAX( XINCR(I,J+PATCH_NUMBER*(L-1)), XVLAIMIN(J) - XI(I,J,L) )
            XI(I,J,L) = XI(I,J,L) + XINCR(I,J+PATCH_NUMBER*(L-1))
            ZBIO_PASS(I,J) = ZBIO_PASS(I,J) + XINCR(I,J+PATCH_NUMBER*(L-1))*XALPH(J)
         ENDDO 
       ELSE
         DO J=1,PATCH_NUMBER           
            if (XI(I,J,L)+ XINCR(I,J+PATCH_NUMBER*(L-1)) .LT. 0) then
                 XINCR(I,J+PATCH_NUMBER*(L-1))=0 
                !XINCR(I,J+PATCH_NUMBER*(L-1))= MAX( XINCR(I,J+PATCH_NUMBER*(L-1)), 0 - XI(I,J,L) )         
            endif
                XI(I,J,L) = XI(I,J,L) + XINCR(I,J+PATCH_NUMBER*(L-1))
          ENDDO    
       ENDIF
     ENDDO
   ENDDO
!
! *** Write analysis results in ASCII file ***
!
 OPEN (unit=111,file='ANAL_INCR',status='unknown',IOSTAT=istat)
 DO I=1,NSIZE_NATURE
   DO J=1,PATCH_NUMBER
     WRITE(111,*) (XI(I,J,L),L=1,NVAR), (XINCR(I,J+PATCH_NUMBER*(L-1)),L=1,NVAR)
   ENDDO
 ENDDO
 CLOSE(UNIT=111)
!
!   Write analysis in PREP file
!
 CALL END_IO_SURF_n(YPROGRAM)
 CALL IO_BUFF_CLEAN_n
! CFILEOUT='PREP'
 CALL INIT_IO_SURF_n(YPROGRAM,'NATURE','ISBA  ','WRITE')
 CALL IO_BUFF_CLEAN_n
 DO L = 1,NVAR
   CALL WRITE_SURF(YPROGRAM,XVAR(L),XI(:,:,L),IRESP,HCOMMENT=PREFIX(L))
   IF (XVAR(NVAR).EQ.'LAI') THEN
     CALL WRITE_SURF(YPROGRAM,CBIO,ZBIO_PASS(:,:),IRESP,HCOMMENT=CPREFIX_BIO)         
   ENDIF   
 ENDDO

 CALL END_IO_SURF_n(YPROGRAM)
!
! Analysis of B (for use in next cycle)
! Ba = (I-KH)Bf
! K = BfHT{K1}**-1
!
 DO I=1,NSIZE_NATURE
       !K1 = (R+H.B.HT) (calculate inverse -> output goes to K1)
       CALL INVERSE_MATRIX(NOBS,K1(I,:,:),ZP(I,:))
       GAIN(I,:,:) = MATMUL(B(I,:,:),MATMUL(HOT(I,:,:),K1(I,:,:)))
       KH(I,:,:) = MATMUL(GAIN(I,:,:),HO(I,:,:))
       IdKH(I,:,:) = IDENT - KH(I,:,:) 
       KRK(I,:,:) = MATMUL(GAIN(I,:,:),MATMUL(R(I,:,:),TRANSPOSE(GAIN(I,:,:))))
       IF (.NOT.LBFIXED)  B(I,:,:) = MATMUL(IdKH(I,:,:),MATMUL(B(I,:,:),TRANSPOSE(IdKH(I,:,:))))+KRK(I,:,:)
 ENDDO
!
! Write out analysed B (for use in next cycle)
!
 OPEN (unit=111,file='BGROUNDout',status='unknown',IOSTAT=istat)
 DO L=1,NVAR
   DO K=1,NVAR
     DO I=1,NSIZE_NATURE
       DO J=1,PATCH_NUMBER
         DO JJ=1,PATCH_NUMBER
           WRITE (111,*) B(I,J+PATCH_NUMBER*(L-1),JJ+PATCH_NUMBER*(K-1))
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
     DO I=1,NSIZE_NATURE
       DO J=1,PATCH_NUMBER
         WRITE(111,*) HOWR(I,K,J+PATCH_NUMBER*(L-1)),GAIN(I,J+PATCH_NUMBER*(L-1),K)
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
CALL DEALLOC_SURFEX

CALL END_LOG_MPI

IF (LHOOK) CALL DR_HOOK('VARASSIM',1,ZHOOK_HANDLE)
!
#ifdef SFX_MPI
 CALL MPI_FINALIZE(INFOMPI)
#endif
!-------------------------------------------------------------------------------
!
END PROGRAM VARASSIM
