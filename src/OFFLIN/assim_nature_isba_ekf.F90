SUBROUTINE ASSIM_NATURE_ISBA_EKF(HPROGRAM, KI,   &
                                 PT2M,     PHU2M,&
                                 HTEST )

! -----------------------------------------------------------------------------
!
! Land Data Assimilation System based on an Extended Kalman Filter
!
! Revised version : JFM (15 September 2008)
!
! The control vector can be any element of (TG1,TG2,WG1,WG2) - Choice in namelist
!
! The observations can be any element of (T2M,HU2M,WG1) - Choice in namelist
!
! Possibility to evolve the B matrix in the cycling - otherwise SEKF
!
! First version including patches (15 October 2008)
!
! -----------------------------------------------------------------------------
!
 USE MODD_TYPE_DATE_SURF,ONLY : DATE_TIME
 USE MODD_ASSIM,         ONLY : LPRT,LSIM,LBEV,LBFIXED,NOBSTYPE,YERROBS,INCO,&
                                IVAR,NVAR,NVARMAX,TPRT_M,XSIGMA_M,XVAR_M,PREFIX_M,   &
                                INCV,SCALE_Q
 USE MODN_IO_OFFLINE,    ONLY : CPGDFILE,CPREPFILE
!
#ifdef LFI
 USE MODD_IO_SURF_LFI,   ONLY : CFILEIN_LFI, CFILEOUT_LFI
#endif
!
 USE MODD_SURF_ATM_n,    ONLY : NSIZE_NATURE,NDIM_FULL
!
 USE YOMHOOK,            ONLY : LHOOK,DR_HOOK
 USE PARKIND1,           ONLY : JPRB
!
 USE MODI_ABOR1_SFX
 USE MODI_INIT_IO_SURF_n 
 USE MODI_READ_SURF
 USE MODI_END_IO_SURF_n 
 USE MODI_IO_BUFF_CLEAN_n
 USE MODI_WRITE_SURF
 USE MODI_ADD_FORECAST_TO_DATE_SURF
!
! -----------------------------------------------------------
!
 IMPLICIT NONE
 CHARACTER(LEN=6),    INTENT(IN)            :: HPROGRAM     ! program calling surf. schemes
 INTEGER,             INTENT(IN)            :: KI
 REAL, DIMENSION(KI), INTENT(IN)            :: PT2M
 REAL, DIMENSION(KI), INTENT(IN)            :: PHU2M
 CHARACTER(LEN=2),    INTENT(IN)            :: HTEST        ! must be equal to 'OK'

 INTEGER                                    :: NOBS
 CHARACTER(LEN=28)                          :: YNAMELIST = 'OPTIONS.nam                 '
!
!    Declarations of local variables
!
 CHARACTER(LEN=3), PARAMETER                :: YINIT     = 'ALL'
 CHARACTER(LEN=6)                           :: YPROGRAM  = 'LFI   '
 CHARACTER(LEN=6)                           :: YPROGRAM2 = 'FA    '
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
 INTEGER                                    :: NBOUTPUT                   ! Number of time step
 INTEGER                                    :: ISTEP                      ! 
 REAL,DIMENSION(:,:,:),ALLOCATABLE          :: YF                         ! Vector of model observations (averaged)
 REAL,DIMENSION(:,:,:,:),ALLOCATABLE        :: YF_PATCH                   ! vector of model observations (for each pacth)
 REAL,DIMENSION(:,:,:,:),ALLOCATABLE        :: XF                         ! Vector of forecast control variables
 REAL,DIMENSION(:,:,:),ALLOCATABLE          :: XI                         ! Vector of control variables (at beginning of timestep)
 REAL,DIMENSION(:,:,:,:),ALLOCATABLE        :: HO                         ! Jacobian of observation operator
 REAL,DIMENSION(:,:,:,:),ALLOCATABLE        :: HOWR                       ! copy of HO for writing out
 REAL,DIMENSION(:,:,:,:),ALLOCATABLE        :: HOT                        ! Transpose of HO
 REAL,DIMENSION(:,:,:),ALLOCATABLE          :: R                          ! covariance matrix of observation errors
 REAL,DIMENSION(:,:,:,:),ALLOCATABLE        :: B                          ! background error covariance matrix
 REAL,DIMENSION(:,:),ALLOCATABLE            :: IDENT                      ! identitiy matrix, used for Ba
 REAL, DIMENSION(:,:,:,:), ALLOCATABLE      :: GAIN                       ! Kalman gain (used expicitly for Ba)               
 REAL, DIMENSION(:,:,:,:), ALLOCATABLE      :: LTM                        ! linear tangent matrix for teh f'ward model
 REAL, DIMENSION(:,:,:,:), ALLOCATABLE      :: Q                          ! model error matrix
 REAL,DIMENSION(:,:),ALLOCATABLE            :: YO                         ! vector of observations
 REAL,DIMENSION(:),ALLOCATABLE              :: ZCLAY                      ! Pourcentage of clay (varies from 0 to 1)
 REAL,DIMENSION(:),ALLOCATABLE              :: ZSAND
 REAL,DIMENSION(:),ALLOCATABLE              :: COFSWI                     ! The difference (Wfc - Wwilt)
 REAL,DIMENSION(:),ALLOCATABLE              :: SMSAT                      !  
 REAL,DIMENSION(:,:,:),ALLOCATABLE          :: ZEPS                       ! The perturbation amplitude
 REAL,DIMENSION(:,:,:),ALLOCATABLE          :: XINCR                      ! Analysis increment
 REAL,DIMENSION(:,:,:),ALLOCATABLE          :: SIMOBS                     ! Simulated temperature and relative humidity (available per patch ?)
 REAL, DIMENSION (:,:),ALLOCATABLE          :: XPATCH                     ! Fraction covered by each patch
 REAL,DIMENSION(:,:,:),ALLOCATABLE          :: SIMMOD                     ! Control variables (for B propagation)
 REAL,DIMENSION(:,:),ALLOCATABLE            :: VECT                       ! The analysed variable
 CHARACTER(LEN=200)                         :: NMFILE_CANARI              ! Name of the observation, perturbed or reference file
 CHARACTER(LEN=9)                           :: HFNAME
 CHARACTER(LEN=17)                          :: LFNAME
 INTEGER                                    :: IND                        
 INTEGER                                    :: LTEST                       
!
! Local Matrix for Analysis calculation
!
 REAL,DIMENSION(:,:,:),ALLOCATABLE          :: K1
 REAL,DIMENSION(:,:),ALLOCATABLE            :: ZX,ZB,ZP
 REAL,DIMENSION(:,:),ALLOCATABLE            :: YOWR
 INTEGER                                    :: I,J,K,KK,L
! 
 REAL,DIMENSION(:),ALLOCATABLE              :: TPRT                      ! The perturbation amplitude
 REAL,DIMENSION(:),ALLOCATABLE              :: XSIGMA                    ! covariance of background errors if B is fixed
!                                                                        ! covariance of model errors if B evolving 
 CHARACTER(LEN=3),DIMENSION(:),ALLOCATABLE  :: XVAR ! X is ctrl          ! Name of control variables (syntax of surfex in PREP.txt file )
 CHARACTER(LEN=100),DIMENSION(:),ALLOCATABLE:: PREFIX                    ! The prefix of the control variables (in PREP.txt file) 
 CHARACTER(LEN=10),DIMENSION(:),ALLOCATABLE :: XOBS                     ! Identifier for simulated observations
 CHARACTER(LEN=3)                           :: YREAD
!
 INTEGER                                    :: NDIM

 INTEGER                                    :: ILUOUT                    ! ascii output unit number
 INTEGER                                    :: ILUNAM                    ! namelist unit number
 INTEGER                                    :: ISTAT                    
 LOGICAL                                    :: GFOUND                    ! return logical when reading namelist
 
 REAL,DIMENSION(:),ALLOCATABLE              :: PT2M_O
 REAL,DIMENSION(:),ALLOCATABLE              :: PHU2M_O
 REAL,DIMENSION(:),ALLOCATABLE              :: PEPES

 INTEGER                                    :: COMPT 
 REAL(KIND=JPRB)                            :: ZHOOK_HANDLE

 IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_ISBA_EKF',0,ZHOOK_HANDLE)

 IF (HTEST/='OK') THEN
   CALL ABOR1_SFX('ASSIM_NATURE_ISBA_EKF: FATAL ERROR DURING ARGUMENT TRANSFER')
 END IF

 PRINT *
 PRINT *,'   --------------------------'
 PRINT *,'   |   ENTERING  VARASSIM   |'
 PRINT *,'   --------------------------'
 PRINT *
!
!      0.2 Allocate arrays depending on the control vector dimension
!
 ALLOCATE (TPRT(NVAR))
 ALLOCATE (XSIGMA(NVAR))
 ALLOCATE (XVAR(NVAR))
 ALLOCATE (PREFIX(NVAR))
!
 IF (SUM(INCV) /= NVAR) THEN
  PRINT *,' INCONSISTENCY in set-up of CONTROL VARIABLES',SUM(INCV),NVAR
  STOP
 ENDIF
!
 J = 1
 DO I = 1,NVARMAX
  IF (INCV(I) == 1 .AND. J <= NVAR ) THEN
   TPRT(J) = TPRT_M(I)
   XSIGMA(J) = XSIGMA_M(I)
   XVAR(J) = XVAR_M(I)
   PREFIX(J) = PREFIX_M(I)
   J = J + 1
  ENDIF
 ENDDO
!
!      0.3 Allocate arrays depending on the number of observation type
!
 ALLOCATE (XOBS(NOBSTYPE))
!
!      1.    Initializations
!
 ILOBS=55
!
!
!   LFI file handling
!
#ifdef LFI
 CFILEIN_LFI = CPREPFILE        ! output of PREP
#endif
!
!   Find current time
!
 CALL INIT_IO_SURF_n(YPROGRAM,'FULL  ','SURF  ','READ ')
 CALL READ_SURF(YPROGRAM,'DTCUR',TTIME,IRESP)
 CALL END_IO_SURF_n(YPROGRAM)
 CALL IO_BUFF_CLEAN_n
!
!   Read number of patches
!
#ifdef LFI
 CFILEIN_LFI = CPGDFILE        ! output of PGD
#endif
 CALL INIT_IO_SURF_n(YPROGRAM,'NATURE','ISBA  ','READ ')
 CALL READ_SURF(YPROGRAM,'PATCH_NUMBER',PATCH_NUMBER,IRESP)
 CALL END_IO_SURF_n(YPROGRAM)
 CALL IO_BUFF_CLEAN_n

!
!   Define prefixes for simulated observations
! 
 print *,'number of patches =',PATCH_NUMBER
 IF (PATCH_NUMBER > 1) THEN
   XOBS(1) = 'T2M_PATCH'
   XOBS(2) = 'HU2M_PATCH'
 ELSE
   XOBS(1) = 'T2M_ISBA'
   XOBS(2) = 'HU2M_ISBA'
 ENDIF
 IF (NOBSTYPE > 2) XOBS(3) = 'WG1'
!
 ALLOCATE(XI(NSIZE_NATURE,PATCH_NUMBER,NVAR))
 ALLOCATE(VECT(NSIZE_NATURE,PATCH_NUMBER))
 ALLOCATE(ZSAND(NSIZE_NATURE))
 ALLOCATE(SMSAT(NSIZE_NATURE))
!
! Read in control variables
!
#ifdef LFI
 CFILEIN_LFI = CPREPFILE        ! output of PREP
#endif
 CALL INIT_IO_SURF_n(YPROGRAM,'NATURE','ISBA  ','READ ')
 DO L = 1,NVAR
   CALL READ_SURF(YPROGRAM,XVAR(L),XI(:,:,L),IRESP)
   PRINT *,XVAR(L),' - initial ',XI(1,1,L)
 ENDDO
 CALL END_IO_SURF_n(YPROGRAM)
 CALL IO_BUFF_CLEAN_n

!
 ALLOCATE(COFSWI(NSIZE_NATURE))
 ALLOCATE(ZCLAY(NSIZE_NATURE))
 ALLOCATE(B(NSIZE_NATURE,PATCH_NUMBER,NVAR,NVAR))
!
!   Read CLAY fraction to  compute the SWI range (Wfc - Wwilt)
!   (XSIGMA is defined in terms of SWI), need to convert to equivalent v/v
!   using same clay fraction in both layers
!   Read SAND fraction to compute the saturation for conversion of ERS SWI
!
#ifdef LFI
 CFILEIN_LFI = CPGDFILE        ! output of PGD
#endif
 CALL INIT_IO_SURF_n(YPROGRAM,'NATURE','ISBA  ','READ ')
 CALL READ_SURF(YPROGRAM,'CLAY',ZCLAY(:),IRESP)
 CALL READ_SURF(YPROGRAM,'SAND',ZSAND(:),IRESP)
 CALL END_IO_SURF_n(YPROGRAM)
 CALL IO_BUFF_CLEAN_n

!
 DO I=1,NSIZE_NATURE
   COFSWI(I)=0.001*(89.0467*((100.*ZCLAY(I))**0.3496)-37.1342*((100.*ZCLAY(I))**0.5))
   SMSAT(I)=0.001*(-1.08*100*ZSAND(I)+494.305)
 ENDDO
!
!   Frequency of assimilation cycling and data availability
!
 IHOUR = 6
 NBOUTPUT = 1
 PTSTEP_OUTPUT = 6
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
!   (write to CPREPFILE then exit)
!   Element of the control variable chosen by IVAR in namelist
!

 IF ( LPRT ) THEN
   ! read in control variable
#ifdef LFI
   CFILEOUT_LFI=CPREPFILE
#endif
   DO I=1,NSIZE_NATURE
     DO J=1,PATCH_NUMBER
       VECT(I,J) = XI(I,J,IVAR) + TPRT(IVAR)*XI(I,J,IVAR) 
     ENDDO
   ENDDO
   PRINT *,XVAR(IVAR),' - ptrbd', vect(1,1) 
   CALL INIT_IO_SURF_n(YPROGRAM,'NATURE','ISBA  ','WRITE')
   CALL WRITE_SURF(YPROGRAM,XVAR(IVAR),VECT,IRESP,PREFIX(IVAR))
   CALL END_IO_SURF_n(YPROGRAM)
   CALL IO_BUFF_CLEAN_n
!
! Initialisation of B matrix
!
   B(:,:,:,:) = 0.0
   DO L=1,NVAR
     DO I=1,NSIZE_NATURE
       DO J=1,PATCH_NUMBER
         IF (XVAR(L) == 'WG2') B(I,J,L,L) = XSIGMA(L)*XSIGMA(L)*COFSWI(I)*COFSWI(I)
         IF (XVAR(L) == 'WG1') B(I,J,L,L) = XSIGMA(L)*XSIGMA(L)*COFSWI(I)*COFSWI(I)
         IF (XVAR(L) == 'TG2') B(I,J,L,L) = XSIGMA(L)*XSIGMA(L)
         IF (XVAR(L) == 'TG1') B(I,J,L,L) = XSIGMA(L)*XSIGMA(L)
       ENDDO
     ENDDO
   ENDDO
!
   PRINT *,'writing out initial B'
   OPEN (unit=111,file='BGROUNDin0',status='new',IOSTAT=istat)
   IF (istat .NE. 0) THEN
     STOP 'BGROUNDin0 already written'
   ELSE
     WRITE (111,*) ((B(I,J,:,:),J=1,PATCH_NUMBER),I=1,NSIZE_NATURE)
     CLOSE(111)
   ENDIF
!
   PRINT *
   PRINT *,'   -----------------------------------'
   PRINT *,'   |   EXITING VARASSIM AFTER LPRT   |'
   PRINT *,'   -----------------------------------'
   PRINT *
   STOP
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
#ifdef LFI
   CFILEIN_LFI = CPREPFILE        ! output of PREP
#endif
   CALL INIT_IO_SURF_n(YPROGRAM,'NATURE','ISBA  ','READ ')
!
   ALLOCATE(SIMOBS(NSIZE_NATURE,PATCH_NUMBER,NOBSTYPE))
   ALLOCATE(SIMMOD(NSIZE_NATURE,PATCH_NUMBER,NVAR))
!

   
   DO I = 1,NOBSTYPE
     CALL READ_SURF(YPROGRAM,XOBS(I),SIMOBS(:,:,I),IRESP) 
   ENDDO
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
   CALL END_IO_SURF_n(YPROGRAM)
   CALL IO_BUFF_CLEAN_n

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
   STOP
 ENDIF
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
   ALLOCATE(LTM(NSIZE_NATURE,PATCH_NUMBER,NVAR,NVAR))
   ALLOCATE(ZEPS(NSIZE_NATURE,PATCH_NUMBER,NVAR))
   ALLOCATE(XF(NSIZE_NATURE,PATCH_NUMBER,NVAR+1,NVAR))

!
   PRINT *, 'evolving B to time ', IYEAR,IMONTH,IDAY,IHOUR  
!
   OPEN (unit=111,file='BGROUNDin',status='unknown',IOSTAT=istat)
   READ (111,*) ((B(I,J,:,:),J=1,PATCH_NUMBER),I=1,NSIZE_NATURE)
   CLOSE(111)
   print *,'read previous B matrix  ==>',B(1,1,1,1),NVAR
!
! Calculate the TLM of the forecast model
! 
! a) read in perturbed forecasts
!
   DO L=1,NVAR
     WRITE(LCHAR,'(I1) ') L
     NMFILE_CANARI='MDSIMU_PERT_'//LCHAR//'_'
     CALL GET_FILE_NAME(IYEAR,IMONTH,IDAY,IHOUR,NMFILE_CANARI)
     PRINT *,'--> reading in ptb file: ',NMFILE_CANARI
     OPEN(UNIT=111,FILE=NMFILE_CANARI,FORM='FORMATTED',STATUS='OLD')
     DO I=1,NSIZE_NATURE
       DO J=1,PATCH_NUMBER 
         READ(111,*) (XF(I,J,L+1,K),K=1,NVAR)
       ENDDO
     ENDDO
     CLOSE(111)
     PRINT *, 'read in ptbd forecasts for L ',L, XF(1,1,L+1,1)
   ENDDO
!
! b) read in reference forecasts
!
   NMFILE_CANARI='MDSIMU_REFR_'
   PRINT *,'--> reading in reference file'
   CALL GET_FILE_NAME(IYEAR,IMONTH,IDAY,IHOUR,NMFILE_CANARI)
   OPEN(UNIT=111,FILE=NMFILE_CANARI,FORM='FORMATTED',STATUS='OLD')
!
   DO I=1,NSIZE_NATURE
     DO J=1,PATCH_NUMBER
       READ(111,*) (XF(I,J,1,K),K=1,NVAR)
     ENDDO
   ENDDO
   CLOSE(111)
   PRINT *, 'read in ref forecasts for L',L, XF(1,1,1,1)
!
! c) calculate initial perturbation
! 
   DO L=1,NVAR
     WHERE(XI(:,:,L).NE.999.0)
       ZEPS(:,:,L) = TPRT(L)*XI(:,:,L)
     ENDWHERE
   ENDDO
!
! d) calculate LTM - one patch only
!
   DO L=1,NVAR ! control variable (x at previous time step)
     DO I=1,NSIZE_NATURE 
       DO J=1,PATCH_NUMBER 
         DO K=1,NVAR     
           LTM(I,J,K,L) = (XF(I,J,L+1,K) - XF(I,J,1,K))/ZEPS(I,J,L) ! Jacobian of fwd model
         ENDDO
       ENDDO
     ENDDO
   ENDDO
!
! evolve B - one patch only
!
   DO I=1,NSIZE_NATURE
     DO J=1,PATCH_NUMBER
       B(I,J,:,:) = MATMUL(LTM(I,J,:,:),MATMUL(B(I,J,:,:),TRANSPOSE(LTM(I,J,:,:))))
     ENDDO
   ENDDO
!
   PRINT *,'LTM d(wg2)/d(wg2)', LTM(1,1,1,1)
!
! write out the LTM for the forward model
!
   DO L=1,NVAR
     DO K=1,NVAR 
       WRITE(LCHAR,'(I1)') K
       LFNAME='LTM_del'//XVAR(K)//'_del'//XVAR(L)
       OPEN(UNIT=111,FILE=LFNAME,FORM='FORMATTED',STATUS='UNKNOWN',POSITION='APPEND')
       DO I=1,NSIZE_NATURE
         DO J=1,PATCH_NUMBER
           WRITE (111,*) LTM(I,J,K,L)
         ENDDO
       ENDDO
       CLOSE(111)
     ENDDO
   ENDDO
!
! write out current B (for use in next cycle)
!
   print *,'store B matrix after TL evolution ==>',B(1,1,1,1)
   PRINT *,'writing out B'
   OPEN (unit=111,file='BGROUNDout',status='unknown')
   WRITE (111,*) ((B(I,J,:,:),J=1,PATCH_NUMBER),I=1,NSIZE_NATURE)
   CLOSE(111)
!
   DEALLOCATE(B)
   DEALLOCATE(LTM)
   DEALLOCATE(ZEPS)
   DEALLOCATE(XF)
!
   PRINT *
   PRINT *,'   -----------------------------------'
   PRINT *,'   |   EXITING VARASSIM AFTER LBEV   |'
   PRINT *,'   -----------------------------------'
   PRINT *
   STOP
 ENDIF
! ====================================================================
!
! if not LSIM,LPTR, or LBEV proceed with analysis
!
!   Number of available observation files
 PRINT *,'--> through to analysis section'



! NOBS=2


!
! NOBS=0
!R ZTIME=PTSTEP_OUTPUT
!
! DO ISTEP=1,NBOUTPUT
!R   CALL ADD_FORECAST_TO_DATE_SURF(IYEAR,IMONTH,IDAY,ZTIME)
!R   ZTIME = ZTIME + PTSTEP_OUTPUT
!R   NMFILE_CANARI='CANARI_NATURE_'
!R   CALL GET_FILE_NAME(IYEAR,IMONTH,IDAY,IHOUR,NMFILE_CANARI)
!R   print *,'observation file name =>',NMFILE_CANARI
!R   OPEN(UNIT=ILOBS,FILE=NMFILE_CANARI,FORM='UNFORMATTED',STATUS='OLD',ERR=22)
!    NOBS=NOBS+NOBSTYPE
!   22 CONTINUE
!R   CLOSE(ILOBS)
! ENDDO
!
!R IF (NOBS == 0) THEN
!R   PRINT *,'No observations read in OBS file - stop'
!R   STOP
!R ENDIF
!

 


 ALLOCATE (PT2M_O(NSIZE_NATURE))
 ALLOCATE (PHU2M_O(NSIZE_NATURE))

 COMPT=0
 DO I=1,KI
    COMPT=COMPT+1
    PT2M_O(COMPT)=PT2M(I) 
    PHU2M_O(COMPT)=PHU2M(I)
 ENDDO


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
 ALLOCATE(YF(NSIZE_NATURE,NVAR+1,NOBSTYPE))
 ALLOCATE(YF_PATCH(NSIZE_NATURE,PATCH_NUMBER,NVAR+1,NOBSTYPE))
 ALLOCATE(ZEPS(NSIZE_NATURE,PATCH_NUMBER,NVAR))
 ALLOCATE(XF(NSIZE_NATURE,PATCH_NUMBER,NVAR+1,NVAR))
 ALLOCATE(XPATCH(NSIZE_NATURE,PATCH_NUMBER))
!
! Read fraction of each patch (=> LPGD=T in OPTIONS.nam)
!
#ifdef LFI
 CFILEIN_LFI = CPGDFILE        ! output of PGD
#endif
 CALL INIT_IO_SURF_n(YPROGRAM,'NATURE','ISBA  ','READ ')
 IF (PATCH_NUMBER > 1) THEN
   CALL READ_SURF(YPROGRAM,'PATCH',XPATCH,IRESP)
 ELSE
   XPATCH(:,1) = 1.0 
 ENDIF
 CALL END_IO_SURF_n(YPROGRAM)
 CALL IO_BUFF_CLEAN_n

!


!   Initial values (to be analysed)
!   Observations
!
 ALLOCATE(YO(NSIZE_NATURE,NOBSTYPE))
 ALLOCATE(YOWR(NDIM_FULL,NOBSTYPE))
!
!   Temporary vectors used by the EKF approach
!
 ALLOCATE(ZB(NSIZE_NATURE,NOBSTYPE))
 ALLOCATE(HO(NSIZE_NATURE,PATCH_NUMBER,NOBSTYPE,NVAR))
 ALLOCATE(HOWR(NSIZE_NATURE,PATCH_NUMBER,NOBSTYPE,NVAR))
 ALLOCATE(HOT(NSIZE_NATURE,PATCH_NUMBER,NVAR,NOBSTYPE))
 ALLOCATE(R(NSIZE_NATURE,NOBSTYPE,NOBSTYPE))
 ALLOCATE(IDENT(NVAR,NVAR))
 ALLOCATE(GAIN(NSIZE_NATURE,PATCH_NUMBER,NVAR,NOBSTYPE))
 ALLOCATE(ZX(NSIZE_NATURE,NOBSTYPE))
 ALLOCATE(ZP(NSIZE_NATURE,NOBSTYPE))
 ALLOCATE(XINCR(NSIZE_NATURE,PATCH_NUMBER,NVAR))
 ALLOCATE(K1(NSIZE_NATURE,NOBSTYPE,NOBSTYPE))
!
!   Initialisations
!
 HO(:,:,:,:)    = 999.0     ! Linearized observation matrix
 HOWR(:,:,:,:)  = 999.0  
 R(:,:,:)       = 0.0       ! Observation error matrix
 B(:,:,:,:)     = 0.0       ! Background error matrix
 YOWR(:,:)      = 999.0     ! Observation vector on the full grid to be written on file
! YO(:,:)        = 999.0
 ZB(:,:)        = 999.0     ! Innovation vector
 YF(:,:,:)      = 0.0       ! Tile averaged simulated observation vector
 OBSCOUNT       = 0
!
 IDENT=0                    ! identity matrix
 DO I = 1,NVAR
   IDENT(I,I) = 1
 ENDDO
!
!   Calculate delta for control variables 
!
 DO L=1,NVAR
   ZEPS(:,:,L) = TPRT(L)*XI(:,:,L)
 ENDDO
!
! NOBS=2
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
!R   NMFILE_CANARI='CANARI_NATURE_'
!R   CALL GET_FILE_NAME(IYEAR,IMONTH,IDAY,IHOUR,NMFILE_CANARI)
!
!   Open observation files
!
!R   OPEN(UNIT=ILOBS,FILE=NMFILE_CANARI,FORM='FORMATTED',STATUS='OLD',ERR=10)
!
!   NOBS=NOBS+NOBSTYPE 
!
!   If it exists, read T2m and HU2m observations
!
!R     DO I=1,NSIZE_NATURE
!R     READ (ILOBS,*)  (YO(I,J),J=1,NOBSTYPE)        
!R     IF (YO(I,3) == 110.0) YO(I,3) = 999.0 ! frozen soils
!R     ENDDO



      YO(:,1)=PT2M_O(:)
      YO(:,2)=PHU2M_O(:)

   PRINT *,'read in obs: ', YO(1,1), YO(1,2), NOBS



!
!   Calculate perturbations at the date of observations and innovation
!   (2nd dim of YF: index 1 is reference, 2,3... perturbed by cntrl var 1,2... 
!
   DO L=1,NVAR
     WRITE(LCHAR,'(I1) ') L
     NMFILE_CANARI='OBSIMU_PERT_'//LCHAR//'_'
     CALL GET_FILE_NAME(IYEAR,IMONTH,IDAY,IHOUR,NMFILE_CANARI)
     PRINT *,'--> reading in ptb file: ',NMFILE_CANARI
     OPEN(UNIT=111,FILE=NMFILE_CANARI,FORM='FORMATTED',STATUS='OLD')
     DO I=1,NSIZE_NATURE
       DO J=1,PATCH_NUMBER
!         READ(111,*) (YF_PATCH(I,J,L+1,K+NOBS-NOBSTYPE),K=1,NOBSTYPE)
         READ(111,*) (YF_PATCH(I,J,L+1,K),K=1,NOBSTYPE)
       ENDDO
     ENDDO
     CLOSE(111)
     PRINT *, 'read in ptbd forecasts for H',L, YF_PATCH(1,1,L+1,1), YF_PATCH(1,1,L+1,2)
   ENDDO



!
   NMFILE_CANARI='OBSIMU_REFR_'
   PRINT *,'--> reading in reference file'
   CALL GET_FILE_NAME(IYEAR,IMONTH,IDAY,IHOUR,NMFILE_CANARI)
   OPEN(UNIT=111,FILE=NMFILE_CANARI,FORM='FORMATTED',STATUS='OLD')
!
   DO I=1,NSIZE_NATURE
     DO J=1,PATCH_NUMBER
!       READ(111,*) (YF_PATCH(I,J,1,K+NOBS-NOBSTYPE),K=1,NOBSTYPE)
       READ(111,*) (YF_PATCH(I,J,1,K),K=1,NOBSTYPE)
     ENDDO
   ENDDO
   CLOSE(111)
   PRINT *, 'read in ref forecasts for', YF_PATCH(1,1,1,1),YF_PATCH(1,1,1,2) 
   10  CONTINUE ! if T2m HU2m observations does not exist
   CLOSE(ILOBS)



!
!  Mean simulated obs averaged over tiles
!
   DO I=1,NSIZE_NATURE
     DO J=1,PATCH_NUMBER
       IF (XPATCH(I,J) > 0.0) THEN
         YF(I,:,:) = YF(I,:,:) + XPATCH(I,J)*YF_PATCH(I,J,:,:)
       ENDIF
     ENDDO
   ENDDO
!
 ENDDO TIMELOOP



!
! Rescale the satellite derived soil moisture observations
!
!R DO I=1,NSIZE_NATURE
!R   IF (YO(I,3) /= 999.0) THEN 
! CDF matching on ERS (revised formula)
!     YO(I,3) = (1.16117E-7*YO(I,3)**5 - 2.85771E-5*YO(I,3)**4 + 0.00235291*YO(I,3)**3 &
! &            - 0.070169*YO(I,3)**2   + 1.261382*YO(I,3) + 12.2607)*SMSAT(I)/100.
! CDF matching for ASCAT soil moisture
!R     YO(I,3) = (8.80461E-8*YO(I,3)**5 - 2.21598E-5*YO(I,3)**4 + 0.00188043*YO(I,3)**3 &
!R &            - 0.0575883*YO(I,3)**2   + 1.0249301*YO(I,3) + 15.7502)*SMSAT(I)/100.
!R   ENDIF
!R ENDDO
!




! SET OBSERVATION ERROR 
!
 DO I=1,NSIZE_NATURE
   DO K=1,NOBSTYPE
     IF (K .LT. 3) THEN
       R(I,K,K) = YERROBS(K)*YERROBS(K)
     ELSE
!  convert R for wg1 from SWI  to abs value
       R(I,K,K) = YERROBS(K)*YERROBS(K)*COFSWI(I)*COFSWI(I) 
     ENDIF
   ENDDO
 ENDDO
!
! WRITE OUT OBS AND YERROR FOR DIAGNOSTIC PURPOSES
!
 OPEN (unit=111,file='OBSERRORout',status='unknown',IOSTAT=istat)
 WRITE (111,*) (R(I,2,2), I=1,NSIZE_NATURE)
 CLOSE(111)
 OPEN (unit=111,file='OBSout',status='unknown',IOSTAT=istat)
 WRITE (111,*) (YO(I,2), I=1,NSIZE_NATURE)
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
         IF (XVAR(L) == 'WG2') B(I,J,L,L) = XSIGMA(L)*XSIGMA(L)*COFSWI(I)*COFSWI(I)
         IF (XVAR(L) == 'WG1') B(I,J,L,L) = XSIGMA(L)*XSIGMA(L)*COFSWI(I)*COFSWI(I)
         IF (XVAR(L) == 'TG2') B(I,J,L,L) = XSIGMA(L)*XSIGMA(L)
         IF (XVAR(L) == 'TG1') B(I,J,L,L) = XSIGMA(L)*XSIGMA(L)
       ENDDO
     ENDDO
   ENDDO
 ELSE
   ALLOCATE(Q(NSIZE_NATURE,PATCH_NUMBER,NVAR,NVAR))
   Q(:,:,:,:) = 0.0
   DO L=1,NVAR
     DO I=1,NSIZE_NATURE
       DO J=1,PATCH_NUMBER
         IF (XVAR(L) == 'WG2') Q(I,J,L,L) = SCALE_Q*SCALE_Q*XSIGMA(L)*XSIGMA(L)*COFSWI(I)*COFSWI(I)
         IF (XVAR(L) == 'WG1') Q(I,J,L,L) = SCALE_Q*SCALE_Q*XSIGMA(L)*XSIGMA(L)*COFSWI(I)*COFSWI(I)
         IF (XVAR(L) == 'TG2') Q(I,J,L,L) = SCALE_Q*SCALE_Q*XSIGMA(L)*XSIGMA(L)
         IF (XVAR(L) == 'TG1') Q(I,J,L,L) = SCALE_Q*SCALE_Q*XSIGMA(L)*XSIGMA(L)
       ENDDO
     ENDDO
   ENDDO
   PRINT *,'reading in B, then updating with Q'
   OPEN (unit=111,file='BGROUNDin',status='unknown',IOSTAT=istat)
   READ (111,*) ((B(I,J,:,:),J=1,PATCH_NUMBER),I=1,NSIZE_NATURE)
   CLOSE(111)
!
! B is the forecast matrix - need to add Q
!
   print *,'B before wg2 wg2 ==> ',sqrt(B(1,1,1,1))/COFSWI(1),B(1,1,1,1)
   print *,'Q value wg2 wg2 ==> ',sqrt(Q(1,1,1,1))/COFSWI(1),Q(1,1,1,1)
   B = B + Q
   print *,'B after wg2 wg2 ==>',sqrt(B(1,1,1,1))/COFSWI(1),B(1,1,1,1)
!
   DEALLOCATE(Q)
 ENDIF
!
 


! Data type selection before assimilation (only if NOBS = NOBSTYPE)
!
 IF (NOBSTYPE == NOBSTYPE) THEN
   DO I=1,NOBSTYPE
     IF (INCO(I) == 0) THEN 
       YO (:,I) = 999.0
       PRINT *,'OBSERVATION TYPE ',XOBS(I),' REMOVED'
     ENDIF
   ENDDO
 ENDIF
! YO(:,:) = 999.0 ! WARNING this instruction remove all observations
!
 PRINT *,'calculating jacobians',NOBS
 DO L=1,NVAR
   DO I=1,NSIZE_NATURE
     DO J=1,PATCH_NUMBER
       DO K=1,NOBSTYPE
         HOWR(I,J,K,L) = XPATCH(I,J)*(YF_PATCH(I,J,L+1,K) - YF_PATCH(I,J,1,K))/ZEPS(I,J,L) 
         IF(YO(I,K) .NE. 999.0) THEN
           HO(I,J,K,L) = XPATCH(I,J)*(YF_PATCH(I,J,L+1,K) - YF_PATCH(I,J,1,K))/ZEPS(I,J,L) ! Jacobian of obs operator
           ZB(I,K) = YO(I,K) - YF(I,1,K)                       ! innovation vector
           OBSCOUNT = OBSCOUNT + 1
         ELSE
! set obs op and obs innovation to zero if no obs are present
           HO(I,J,K,:) = 0.0
           ZB(I,K) = 0.0
         ENDIF
       ENDDO
     ENDDO
   ENDDO
 ENDDO




!-----------------------------------------------------
!
!            ******  SOIL ANALYSIS *******
!
!-----------------------------------------------------
 PRINT *,'PERFORMING ANALYSIS'
 DO I=1,NSIZE_NATURE
   DO J=1,PATCH_NUMBER
     HOT(I,J,:,:) = TRANSPOSE(HO(I,J,:,:))
     K1(I,:,:) = MATMUL(HO(I,J,:,:),MATMUL(B(I,J,:,:),HOT(I,J,:,:))) + R(I,:,:)
     CALL CHOLDC(NOBSTYPE,K1(I,:,:),ZP(I,:))                     ! Cholesky decomposition (1)
     CALL CHOLSL(NOBSTYPE,K1(I,:,:),ZP(I,:),ZB(I,:),ZX(I,:))     ! Cholesky decomposition (2)
     XINCR(I,J,:) = MATMUL(B(I,J,:,:),MATMUL(HOT(I,J,:,:),ZX(I,:))) 
     !TA: A possible safety check for WG2 and WG1. Is something needed? I have seen negative values for soil water
     DO L=1,NVAR
       !IF ( L <= 2 ) THEN
       !  XI(I,J,L) = MAX(0.000001,(XI(I,J,L) + XINCR(I,J,L)))
       !ELSE
         XI(I,J,L) = XI(I,J,L) + XINCR(I,J,L)
       !ENDIF
     ENDDO
   ENDDO
 ENDDO
!
! *** Write analysis results in ASCII file ***
!
 OPEN (unit=111,file='ANAL_INCR',status='unknown',IOSTAT=istat)
 DO I=1,NSIZE_NATURE
   WRITE(111,*) (XI(I,1,J), XINCR(I,1,J),J=1,NVAR),ZB(I,1),ZB(I,2)
 ENDDO
 CLOSE(UNIT=111)

! **********************************************
 DO L=1,NVAR
   PRINT *,'Sum and mean increments for ',TRIM(XVAR(l)),'=',SUM(XINCR(:,1,L)),SUM(XINCR(:,1,L))/REAL(NSIZE_NATURE)
 ENDDO
!
!   Write analysis in LFI file  (PREP.lfi) 
!
#ifdef LFI
 CFILEOUT_LFI=CPREPFILE
#endif
 CALL INIT_IO_SURF_n(YPROGRAM,'NATURE','ISBA  ','WRITE')
 DO L = 1,NVAR
   CALL WRITE_SURF(YPROGRAM,XVAR(L),XI(:,:,L),IRESP,PREFIX(L))
 ENDDO
 CALL END_IO_SURF_n(YPROGRAM)
 CALL IO_BUFF_CLEAN_n

!
! Analysis of B (for use in next cycle)
! Ba = (I-KH)Bf
! K = BfHT{K1}**-1
! K1 needs PATCH dim added
!
 DO I=1,NSIZE_NATURE
   DO J=1,PATCH_NUMBER
!  K1 = (R+H.B.HT) (calculate inverse -> output goes to K1)
     CALL INVERSE_MATRIX(NOBSTYPE,K1(I,:,:),ZP(I,:))
     GAIN(I,J,:,:) = MATMUL(B(I,J,:,:),MATMUL(HOT(I,J,:,:),K1(I,:,:))) 
     IF (.NOT.LBFIXED)  B(I,J,:,:) = MATMUL((IDENT - MATMUL(GAIN(I,J,:,:),HO(I,J,:,:))),B(I,J,:,:)) 
   ENDDO
 ENDDO
!
! Write out analysed B (for use in next cycle)
!
 OPEN (unit=111,file='BGROUNDout',status='unknown',IOSTAT=istat)
 WRITE (111,*) ((B(I,J,:,:),J=1,PATCH_NUMBER),I=1,NSIZE_NATURE)
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
         WRITE(111,*) HOWR(I,J,K,L),GAIN(I,J,L,K)
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
 DEALLOCATE(PT2M_O)
 DEALLOCATE(PHU2M_O)



!
 PRINT *
 PRINT *,'   ---------------------------------------'
 PRINT *,'   |   EXITING VARASSIM AFTER ANALYSIS   |'
 PRINT *,'   ---------------------------------------'
 PRINT *
 PRINT *,'Number of assimilated observations =',OBSCOUNT
 PRINT *
!

 IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_ISBA_EKF',1,ZHOOK_HANDLE)

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Subroutine  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GET_FILE_NAME(IYEAR,IMONTH,IDAY,IHOUR,NMFILE)
USE MODI_TRANS_CHAINE
 IMPLICIT NONE
 INTEGER, INTENT (IN) :: IYEAR, IMONTH, IDAY, IHOUR
 CHARACTER(LEN=50) :: CMONTH,CDAY,CYEAR,CHOUR
 CHARACTER(LEN=200), INTENT(INOUT) :: NMFILE
 REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
 IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_ISBA_EKF:GET_FILE_NAME',0,ZHOOK_HANDLE) 
!
 CALL TRANS_CHAINE(CMONTH,IMONTH,2)
 CALL TRANS_CHAINE(CDAY,IDAY,2)
 CALL TRANS_CHAINE(CYEAR,IYEAR,4)
 CALL TRANS_CHAINE(CHOUR,IHOUR,2)
 IF (IHOUR.eq.0.or.IHOUR.eq.24.or.IHOUR.eq.48) CHOUR='00'
!
 NMFILE = TRIM(NMFILE)//CYEAR(3:4)
 NMFILE = TRIM(NMFILE)//CMONTH
 NMFILE = TRIM(NMFILE)//CDAY
 NMFILE = TRIM(NMFILE)//'H'
 NMFILE = TRIM(NMFILE)//CHOUR
 NMFILE = TRIM(NMFILE)//'.DAT'
!
 IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_ISBA_EKF:GET_FILE_NAME',1,ZHOOK_HANDLE)
!
END SUBROUTINE GET_FILE_NAME

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE INVERSE_MATRIX(N,A,P)
!--------------------------------------------------------
!
! Explicit inversed matrix after Cholesky decomposition
!
!--------------------------------------------------------
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: N
 REAL, DIMENSION (N,N),  INTENT(INOUT) :: A
 REAL, DIMENSION (N),  INTENT(IN)      :: P
 REAL ZSUM
 INTEGER :: I, J, K
 REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
 IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_ISBA_EKF:INVERSE_MATRIX',0,ZHOOK_HANDLE)
!
 DO I=1,N
   A(I,I)=1./P(I)
   DO J=I+1,N
     ZSUM = 0.
     DO K=I,J-1
       ZSUM = ZSUM - A(J,K)*A(K,I)
     ENDDO
     A(J,I) = ZSUM/P(J)
   ENDDO
 ENDDO  
 DO I=1,N
   DO J=I+1,N
      A(I,J) =0.0
   ENDDO
 ENDDO
 A = MATMUL(TRANSPOSE(A),A)
!
 IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_ISBA_EKF:INVERSE_MATRIX',1,ZHOOK_HANDLE)
!
END SUBROUTINE INVERSE_MATRIX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CHOLDC(N,A,P)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: N
 REAL, DIMENSION(N,N), INTENT(INOUT) :: A
 REAL, DIMENSION(N), INTENT(OUT) :: P
 INTEGER :: I
 REAL :: SUMM
 REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
 IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_ISBA_EKF:CHOLDC',0,ZHOOK_HANDLE)
!
 DO I=1,N
   SUMM = A(I,I)- DOT_PRODUCT(A(I,1:I-1),A(I,1:I-1))
   IF (SUMM <= 0.0) THEN
     PRINT*,'ERROR IN CHOLDC'
     STOP
   ENDIF
   P(I) = SQRT(SUMM) 
   A(I+1:N,I) = (A(I,I+1:N) - MATMUL(A(I+1:N,1:I-1),A(I,1:I-1)))/P(I)
 ENDDO
!
 IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_ISBA_EKF:CHOLDC',1,ZHOOK_HANDLE)
!
END SUBROUTINE CHOLDC 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CHOLSL(N,A,P,B,X)
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: N
 REAL, DIMENSION (N,N),  INTENT(IN) :: A
 REAL, DIMENSION (N),  INTENT(IN)   :: P,B
 REAL, DIMENSION (N), INTENT(INOUT) :: X
 INTEGER :: I
 REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
 IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_ISBA_EKF:CHOLSL',0,ZHOOK_HANDLE)
!
 DO I=1,N
   X(I) = (B(I) - DOT_PRODUCT(A(I,1:I-1),X(1:I-1)))/P(I)
 ENDDO
 DO I=N,1,-1
   X(I) = (X(I) - DOT_PRODUCT(A(I+1:N,I),X(I+1:N)))/P(I)
 ENDDO
!
 IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_ISBA_EKF:CHOLSL',1,ZHOOK_HANDLE)
!
END SUBROUTINE CHOLSL

END SUBROUTINE ASSIM_NATURE_ISBA_EKF

