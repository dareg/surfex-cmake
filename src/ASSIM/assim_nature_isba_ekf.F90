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
! Trygve Aspelien, Separating IO  06/2013
!
! -----------------------------------------------------------------------------
!
 USE MODD_TYPE_DATE_SURF,ONLY : DATE_TIME
 USE MODD_ASSIM,         ONLY : LPRT,LSIM,LBEV,LBFIXED,NOBSTYPE,YERROBS,INCO,         &
                              & IVAR,NVAR,NVARMAX,INCV,SCALE_Q,NPRINTLEV,             &
                              & XAT2M_ISBA,XAHU2M_ISBA,                               &
                              & NBOUTPUT, PTSTEP_OUTPUT,YF_PATCH,XF,XOBS,XVAR,TPRT,   &
                              & XSIGMA
 USE MODD_ISBA_n,        ONLY : XTG,XWG,XSAND,XCLAY,TTIME,NPATCH,XPATCH
 USE MODN_IO_OFFLINE,    ONLY : CPGDFILE,CPREPFILE
 USE MODD_SURF_ATM_n,    ONLY : NSIZE_NATURE,NDIM_FULL
 USE YOMHOOK,            ONLY : LHOOK,DR_HOOK
 USE PARKIND1,           ONLY : JPRB
#ifdef ARO   
 USE YOMMP,              ONLY : MYPROC 
#endif
 USE MODI_ABOR1_SFX
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
!
!    Declarations of local variables
!
 INTEGER                                    :: OBSCOUNT
 CHARACTER(LEN=1)                           :: LCHAR
 INTEGER                                    :: IYEAR                      ! current year (UTC)
 INTEGER                                    :: IMONTH                     ! current month (UTC)
 INTEGER                                    :: IDAY                       ! current day (UTC)
 REAL                                       :: ZTIME                      ! current time since start of the run (s)
 INTEGER                                    :: IRESP,PATCH_NUMBER         ! return code
 INTEGER                                    :: ISTEP                      ! 
 INTEGER                                    :: ZMYPROC
 CHARACTER(LEN=7)                           :: CMYPROC
 REAL,DIMENSION(:,:,:),ALLOCATABLE          :: YF                         ! Vector of model observations (averaged) 
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
 REAL,DIMENSION(:,:),ALLOCATABLE            :: VECT                       ! The analysed variable
 CHARACTER(LEN=30)                          :: BGFILE
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
 LOGICAL                                    :: LBEXISTS

 IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_ISBA_EKF',0,ZHOOK_HANDLE)

 IF (HTEST/='OK') THEN
   CALL ABOR1_SFX('ASSIM_NATURE_ISBA_EKF: FATAL ERROR DURING ARGUMENT TRANSFER')
 END IF

 IF ( NPRINTLEV > 0 ) THEN
   WRITE(*,*)
   WRITE(*,*) '   --------------------------'
   WRITE(*,*) '   |   ENTERING  VARASSIM   |'
   WRITE(*,*) '   --------------------------'
   WRITE(*,*)
 ENDIF

#ifdef ARO
 IF ( MYPROC > 0 ) THEN 
   ZMYPROC=MYPROC
 ELSE
   ZMYPROC=1
 ENDIF
#else
 ZMYPROC=1
#endif
 WRITE(CMYPROC(1:7),'(I7.7)') ZMYPROC

 IF (SUM(INCV) /= NVAR) THEN
   WRITE(*,*) 'INCONSISTENCY in set-up of CONTROL VARIABLES',SUM(INCV),NVAR
   CALL ABOR1_SFX('INCONSISTENCY in set-up of CONTROL VARIABLES')
 ENDIF

 PATCH_NUMBER=NPATCH
 IF ( NPRINTLEV > 0 ) WRITE(*,*) 'number of patches =',PATCH_NUMBER
 NOBS=NOBSTYPE

 ALLOCATE(XI(NSIZE_NATURE,PATCH_NUMBER,NVAR))
 ALLOCATE(VECT(NSIZE_NATURE,PATCH_NUMBER))
 ALLOCATE(ZSAND(NSIZE_NATURE))
 ALLOCATE(SMSAT(NSIZE_NATURE))

 ! Set control variables
 DO L=1,NVAR
   SELECT CASE (TRIM(XVAR(L)))
     CASE("TG1")
       XI(:,:,L)=XTG(:,1,:)
     CASE("TG2")
       XI(:,:,L)=XTG(:,2,:)
     CASE("WG1")
       XI(:,:,L)=XWG(:,1,:)
     CASE("WG2")
       XI(:,:,L)=XWG(:,2,:)
     CASE DEFAULT
       CALL ABOR1_SFX("Mapping of "//XVAR(L)//" is not defined in EKF!")
   END SELECT
   IF ( NPRINTLEV > 0 ) WRITE(*,*) XVAR(L),' - initial ',XI(1,1,L)
 ENDDO

 ALLOCATE(COFSWI(NSIZE_NATURE))
 ALLOCATE(ZCLAY(NSIZE_NATURE))
 ALLOCATE(B(NSIZE_NATURE,PATCH_NUMBER,NVAR,NVAR))
 ALLOCATE(ZEPS(NSIZE_NATURE,PATCH_NUMBER,NVAR))

 ! Calculate delta for control variables 
 DO L=1,NVAR
   ZEPS(:,:,L) = TPRT(L)*XI(:,:,L)
 ENDDO

 !   Read CLAY fraction to  compute the SWI range (Wfc - Wwilt)
 !   (XSIGMA is defined in terms of SWI), need to convert to equivalent v/v
 !   using same clay fraction in both layers
 !   Read SAND fraction to compute the saturation for conversion of ERS SWI
 ZSAND=XSAND(:,1)
 ZCLAY=XCLAY(:,1)

 DO I=1,KI
   COFSWI(I)=0.001*(89.0467*((100.*ZCLAY(I))**0.3496)-37.1342*((100.*ZCLAY(I))**0.5))
   SMSAT(I)=0.001*(-1.08*100*ZSAND(I)+494.305)
 ENDDO

 ! Time
 IYEAR  = TTIME%TDATE%YEAR
 IMONTH = TTIME%TDATE%MONTH
 IDAY   = TTIME%TDATE%DAY
 ZTIME  = TTIME%TIME

 ! ----------------------
 ! VARASSIM OPTION : LBEV
 ! ----------------------
 !   Calculate the LTM, and evolve B. 
 !
 IF (LBEV) THEN
   ALLOCATE(LTM(NSIZE_NATURE,PATCH_NUMBER,NVAR,NVAR))

   ! Set the B input file depending of an existing B was found or not
   BGFILE="BGROUNDin."//CMYPROC
   INQUIRE (FILE=TRIM(BGFILE),EXIST=LBEXISTS)
   IF (LBEXISTS) THEN
     OPEN (unit=111,file=BGFILE,status='unknown',IOSTAT=istat)
     READ (111,*) ((B(I,J,:,:),J=1,PATCH_NUMBER),I=1,NSIZE_NATURE)
     CLOSE(111)
     IF ( NPRINTLEV > 0 ) WRITE(*,*) 'read previous B matrix  ==>',B(1,1,1,1),NVAR
   ELSE
     ! Initialization of B
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
     IF ( NPRINTLEV > 0 ) WRITE(*,*) 'Initialized B'
   ENDIF

   ! calculate LTM - one patch only
   DO L=1,NVAR ! control variable (x at previous time step)
     DO I=1,KI 
       DO J=1,PATCH_NUMBER 
         DO K=1,NVAR     
           LTM(I,J,K,L) = (XF(I,J,L+1,K) - XF(I,J,1,K))/ZEPS(I,J,L) ! Jacobian of fwd model
         ENDDO
       ENDDO
     ENDDO
   ENDDO

   ! evolve B - one patch only
   DO I=1,NSIZE_NATURE
     DO J=1,PATCH_NUMBER
       B(I,J,:,:) = MATMUL(LTM(I,J,:,:),MATMUL(B(I,J,:,:),TRANSPOSE(LTM(I,J,:,:))))
     ENDDO
   ENDDO

   IF ( NPRINTLEV > 0 ) WRITE(*,*) 'LTM d(wg2)/d(wg2)', LTM(1,1,1,1)

   ! write out the LTM for the forward model
   DO L=1,NVAR
     DO K=1,NVAR 
       WRITE(LCHAR,'(I1)') K
       LFNAME='LTM_del'//TRIM(XVAR(K))//'_del'//TRIM(XVAR(L))//"."//CMYPROC
       OPEN(UNIT=111,FILE=LFNAME,FORM='FORMATTED',STATUS='UNKNOWN',POSITION='APPEND')
       DO I=1,NSIZE_NATURE
         DO J=1,PATCH_NUMBER
           WRITE (111,*) LTM(I,J,K,L)
         ENDDO
       ENDDO
       CLOSE(111)
     ENDDO
   ENDDO

   ! Write out current B
   IF ( NPRINTLEV > 0 ) THEN
     WRITE(*,*) 'store B matrix after TL evolution ==>',B(1,1,1,1)
     WRITE(*,*) 'writing out B'
   ENDIF
   BGFILE="BGROUNDout_LBEV."//CMYPROC
   OPEN (unit=111,file=BGFILE,status='unknown')
   WRITE (111,*) ((B(I,J,:,:),J=1,PATCH_NUMBER),I=1,NSIZE_NATURE)
   CLOSE(111)

   DEALLOCATE(LTM)
   !--------------------------------------------------------------------
   !
   !   Adding model error to background error matrix 
   !
   !--------------------------------------------------------------------
   ALLOCATE(Q(NSIZE_NATURE,PATCH_NUMBER,NVAR,NVAR))
   Q(:,:,:,:) = 0.0
   DO L=1,NVAR
     DO I=1,NSIZE_NATURE
       DO J=1,PATCH_NUMBER
         IF (TRIM(XVAR(L)) == 'WG2') Q(I,J,L,L) = SCALE_Q*SCALE_Q*XSIGMA(L)*XSIGMA(L)*COFSWI(I)*COFSWI(I)
         IF (TRIM(XVAR(L)) == 'WG1') Q(I,J,L,L) = SCALE_Q*SCALE_Q*XSIGMA(L)*XSIGMA(L)*COFSWI(I)*COFSWI(I)
         IF (TRIM(XVAR(L)) == 'TG2') Q(I,J,L,L) = SCALE_Q*SCALE_Q*XSIGMA(L)*XSIGMA(L)
         IF (TRIM(XVAR(L)) == 'TG1') Q(I,J,L,L) = SCALE_Q*SCALE_Q*XSIGMA(L)*XSIGMA(L)
       ENDDO
     ENDDO
   ENDDO

   ! B is the forecast matrix - need to add Q
   IF ( NPRINTLEV > 0 ) THEN
     WRITE(*,*) 'B before wg2 wg2 ==> ',sqrt(B(1,1,1,1))/COFSWI(1),B(1,1,1,1)
     WRITE(*,*) 'Q value wg2 wg2 ==> ',sqrt(Q(1,1,1,1))/COFSWI(1),Q(1,1,1,1)
   ENDIF
   B = B + Q
   IF ( NPRINTLEV > 0 ) WRITE(*,*) 'B after wg2 wg2 ==>',sqrt(B(1,1,1,1))/COFSWI(1),B(1,1,1,1)
   DEALLOCATE(Q)

 ELSEIF (LBFIXED) THEN
   DO L=1,NVAR
     DO I=1,NSIZE_NATURE
       DO J=1,PATCH_NUMBER
         IF (TRIM(XVAR(L)) == 'WG2') B(I,J,L,L) = XSIGMA(L)*XSIGMA(L)*COFSWI(I)*COFSWI(I)
         IF (TRIM(XVAR(L)) == 'WG1') B(I,J,L,L) = XSIGMA(L)*XSIGMA(L)*COFSWI(I)*COFSWI(I)
         IF (TRIM(XVAR(L)) == 'TG2') B(I,J,L,L) = XSIGMA(L)*XSIGMA(L)
         IF (TRIM(XVAR(L)) == 'TG1') B(I,J,L,L) = XSIGMA(L)*XSIGMA(L)
       ENDDO
     ENDDO
   ENDDO
 ELSE
   CALL ABOR1_SFX("LBEV or LBFIXED should be .TRUE.!")
 ENDIF

 ! ====================================================================
 !
 ! Analysis
 !
 ! ====================================================================

 ALLOCATE (PT2M_O(NSIZE_NATURE))
 ALLOCATE (PHU2M_O(NSIZE_NATURE))

 COMPT=0
 DO I=1,KI
   COMPT=COMPT+1
   PT2M_O(COMPT)=PT2M(I) 
   PHU2M_O(COMPT)=PHU2M(I)
 ENDDO

 !   Time reinitialization 
 IYEAR  = TTIME%TDATE%YEAR
 IMONTH = TTIME%TDATE%MONTH
 IDAY   = TTIME%TDATE%DAY
 ZTIME  = TTIME%TIME

 !  Allocation
 !  Perturbed simulations
 ALLOCATE(YF(NSIZE_NATURE,NVAR+1,NOBSTYPE))

 ! Initial values (to be analysed)
 ! Observations
 ALLOCATE(YO(NSIZE_NATURE,NOBSTYPE))
 ALLOCATE(YOWR(NDIM_FULL,NOBSTYPE))

 ! Temporary vectors used by the EKF approach
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

 ! Initialisations
 HO(:,:,:,:)    = 999.0     ! Linearized observation matrix
 HOWR(:,:,:,:)  = 999.0  
 R(:,:,:)       = 0.0       ! Observation error matrix
 YOWR(:,:)      = 999.0     ! Observation vector on the full grid to be written on file
 ZB(:,:)        = 999.0     ! Innovation vector
 YF(:,:,:)      = 0.0       ! Tile averaged simulated observation vector
 OBSCOUNT       = 0

 IDENT=0                    ! identity matrix
 DO I = 1,NVAR
   IDENT(I,I) = 1
 ENDDO

 ZTIME=PTSTEP_OUTPUT
 ! BEGINNING OF TIME LOOP
 TIMELOOP : DO ISTEP=1,NBOUTPUT
   ! Update date
   CALL ADD_FORECAST_TO_DATE_SURF(IYEAR, IMONTH, IDAY, ZTIME)
   ZTIME = ZTIME + PTSTEP_OUTPUT

      YO(:,1)=PT2M_O(:)
      YO(:,2)=PHU2M_O(:)

   IF ( NPRINTLEV > 0 ) WRITE(*,*) 'read in obs: ', YO(1,1), YO(1,2), NOBS

   ! Mean simulated obs averaged over tiles
   DO I=1,NSIZE_NATURE
     DO J=1,PATCH_NUMBER
       IF (XPATCH(I,J) > 0.0) THEN
         YF(I,:,:) = YF(I,:,:) + XPATCH(I,J)*YF_PATCH(I,J,:,:)
       ENDIF
     ENDDO
   ENDDO
 ENDDO TIMELOOP

 ! SET OBSERVATION ERROR 
 DO I=1,NSIZE_NATURE
   DO K=1,NOBSTYPE
     IF (K .LT. 3) THEN
       R(I,K,K) = YERROBS(K)*YERROBS(K)
     ELSE
       ! convert R for wg1 from SWI  to abs value
       R(I,K,K) = YERROBS(K)*YERROBS(K)*COFSWI(I)*COFSWI(I) 
     ENDIF
   ENDDO
 ENDDO

 ! WRITE OUT OBS AND YERROR FOR DIAGNOSTIC PURPOSES
 OPEN (unit=111,file='OBSERRORout.'//CMYPROC,status='unknown',IOSTAT=istat)
 WRITE (111,*) (R(I,2,2), I=1,NSIZE_NATURE)
 CLOSE(111)
 OPEN (unit=111,file='OBSout.'//CMYPROC,status='unknown',IOSTAT=istat)
 WRITE (111,*) (YO(I,2), I=1,NSIZE_NATURE)
 CLOSE(111)

 ! Data type selection before assimilation (only if NOBS = NOBSTYPE)
 IF (NOBSTYPE == NOBSTYPE) THEN
   DO I=1,NOBSTYPE
     IF (INCO(I) == 0) THEN 
       YO (:,I) = 999.0
       IF ( NPRINTLEV > 0 ) WRITE(*,*) 'OBSERVATION TYPE ',XOBS(I),' REMOVED'
     ENDIF
   ENDDO
 ENDIF

 IF ( NPRINTLEV > 0 ) WRITE(*,*) 'calculating jacobians',NOBS
 DO L=1,NVAR
   DO I=1,NSIZE_NATURE
     DO J=1,PATCH_NUMBER
       DO K=1,NOBSTYPE
         HOWR(I,J,K,L) = XPATCH(I,J)*(YF_PATCH(I,J,L+1,K) - YF_PATCH(I,J,1,K))/ZEPS(I,J,L) 
         IF(YO(I,K) .NE. 999.0) THEN
           HO(I,J,K,L) = XPATCH(I,J)*(YF_PATCH(I,J,L+1,K) - YF_PATCH(I,J,1,K))/ZEPS(I,J,L) ! Jacobian of obs operator
           ZB(I,K) = YO(I,K) - YF(I,1,K)                                                   ! innovation vector
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
 IF ( NPRINTLEV > 0 ) WRITE(*,*) 'PERFORMING ANALYSIS'
 DO I=1,NSIZE_NATURE
   DO J=1,PATCH_NUMBER
     HOT(I,J,:,:) = TRANSPOSE(HO(I,J,:,:))
     K1(I,:,:) = MATMUL(HO(I,J,:,:),MATMUL(B(I,J,:,:),HOT(I,J,:,:))) + R(I,:,:)
     CALL CHOLDC(NOBSTYPE,K1(I,:,:),ZP(I,:))                     ! Cholesky decomposition (1)
     CALL CHOLSL(NOBSTYPE,K1(I,:,:),ZP(I,:),ZB(I,:),ZX(I,:))     ! Cholesky decomposition (2)
     XINCR(I,J,:) = MATMUL(B(I,J,:,:),MATMUL(HOT(I,J,:,:),ZX(I,:))) 
     DO L=1,NVAR

       ! Update the modified values
       XI(I,J,L) = XI(I,J,L) + XINCR(I,J,L)

       ! A possible sanity check
       SELECT CASE (TRIM(XVAR(L)))
         CASE("TG1")
         CASE("TG2")
         CASE("WG1")
           ! Preserve positive values
           !XI(I,J,L) = MAX(0.,(XI(I,J,L)))
         CASE("WG2")
           ! Preserve positive values
           !XI(I,J,L) = MAX(0.,(XI(I,J,L)))
         CASE DEFAULT
           CALL ABOR1_SFX("Mapping of "//TRIM(XVAR(L))//" is not defined in EKF!")
       END SELECT
       ! For no only warn if we have negative values.
       IF ( XI(I,J,L) < 0. ) WRITE(*,*) "WARNING X<0. for ",I,J," for variable ",TRIM(XVAR(L))
     ENDDO
   ENDDO
 ENDDO
 
 ! Write analysis results in ASCII file
 OPEN (unit=111,file='ANAL_INCR.'//CMYPROC,status='unknown',IOSTAT=istat)
 DO I=1,NSIZE_NATURE
   WRITE(111,*) (XI(I,1,J), XINCR(I,1,J),J=1,NVAR),ZB(I,1),ZB(I,2)
 ENDDO
 CLOSE(UNIT=111)

 ! **********************************************
 DO L=1,NVAR
   WRITE(*,*) 'Sum and mean increments for ',TRIM(XVAR(l)),'=',SUM(XINCR(:,1,L)),SUM(XINCR(:,1,L))/REAL(NSIZE_NATURE)

   ! Update the modified values
   SELECT CASE (TRIM(XVAR(L)))
     CASE("TG1")
       XTG(:,1,:)=XI(:,:,L)
     CASE("TG2")
       XTG(:,2,:)=XI(:,:,L)
     CASE("WG1")
       XWG(:,1,:)=XI(:,:,L)
     CASE("WG2")
       XWG(:,2,:)=XI(:,:,L)
     CASE DEFAULT
       CALL ABOR1_SFX("Mapping of "//TRIM(XVAR(L))//" is not defined in EKF!")
   END SELECT
 ENDDO

 ! Analysis of B (for use in next cycle)
 ! Ba = (I-KH)Bf
 ! K = BfHT{K1}**-1
 ! K1 needs PATCH dim added
 DO I=1,NSIZE_NATURE
   DO J=1,PATCH_NUMBER
     ! K1 = (R+H.B.HT) (calculate inverse -> output goes to K1)
     CALL INVERSE_MATRIX(NOBSTYPE,K1(I,:,:),ZP(I,:))
     GAIN(I,J,:,:) = MATMUL(B(I,J,:,:),MATMUL(HOT(I,J,:,:),K1(I,:,:))) 
     IF (.NOT.LBFIXED)  B(I,J,:,:) = MATMUL((IDENT - MATMUL(GAIN(I,J,:,:),HO(I,J,:,:))),B(I,J,:,:)) 
   ENDDO
 ENDDO

 ! Write out analysed B (for use in next cycle)
 BGFILE="BGROUNDout_ASSIM."//CMYPROC
 OPEN (unit=111,file=BGFILE,status='unknown',IOSTAT=istat)
 WRITE (111,*) ((B(I,J,:,:),J=1,PATCH_NUMBER),I=1,NSIZE_NATURE)
 CLOSE(111)

 ! **** Write out the observation operator + Gain matrix ****
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
 DEALLOCATE(XF)
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

 IF ( NPRINTLEV > 0 ) THEN
   WRITE(*,*)
   WRITE(*,*) '   ---------------------------------------'
   WRITE(*,*) '   |   EXITING VARASSIM AFTER ANALYSIS   |'
   WRITE(*,*) '   ---------------------------------------'
   WRITE(*,*)
   WRITE(*,*) 'Number of assimilated observations =',OBSCOUNT
   WRITE(*,*)
 ENDIF
 IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_ISBA_EKF',1,ZHOOK_HANDLE)

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Subroutine  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
   IF (SUMM <= 0.0) CALL ABOR1_SFX('ASSIM_NATURE_ISBA_EKF: ERROR IN CHOLDC')
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

