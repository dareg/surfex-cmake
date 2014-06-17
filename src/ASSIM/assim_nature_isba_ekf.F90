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
! Alina Barbu: bug correction of B matrix, otherwise understimation of the gain matrix (11/2013) 
! Alina Barbu: equivalent analysis of B matrix to ensure its symetric and positiv definiteness properties (11/2013) 
  
! -----------------------------------------------------------------------------
!
USE MODD_ASSIM,         ONLY : LBEV, LBFIXED, NOBSTYPE, XERROBS, NNCO, NVAR, NNCV, &
                               XSCALE_Q, NPRINTLEV, CVAR, XTPRT, XSIGMA, NBOUTPUT, &
                               XTSTEP_OUTPUT, XF_PATCH, XF, COBS
! 
USE MODD_SURF_PAR,      ONLY : XUNDEF
!
USE MODD_SURF_ATM_n,    ONLY : NDIM_FULL
USE MODD_ISBA_n,        ONLY : XTG, XWG, XSAND, XCLAY, TTIME, NPATCH, XPATCH
!
#ifdef ARO   
USE YOMMP,              ONLY : MYPROC 
#endif
!
USE YOMHOOK,            ONLY : LHOOK,DR_HOOK
USE PARKIND1,           ONLY : JPRB
!
USE MODI_ABOR1_SFX
USE MODI_ADD_FORECAST_TO_DATE_SURF
!
! -----------------------------------------------------------
!
IMPLICIT NONE
!
CHARACTER(LEN=6),    INTENT(IN)            :: HPROGRAM     ! program calling surf. schemes
INTEGER,             INTENT(IN)            :: KI
REAL, DIMENSION(KI), INTENT(IN)            :: PT2M
REAL, DIMENSION(KI), INTENT(IN)            :: PHU2M
CHARACTER(LEN=2),    INTENT(IN)            :: HTEST        ! must be equal to 'OK'
!
!    Declarations of local variables
!
CHARACTER(LEN=30)  :: YBGFILE
CHARACTER(LEN=17)  :: YLFNAME
CHARACTER(LEN=9)   :: YFNAME
CHARACTER(LEN=7)   :: YMYPROC
CHARACTER(LEN=1)   :: YCHAR
!
! Local Matrix for Analysis calculation
!
!  Allocation
!  Perturbed simulations
!
! Initial values (to be analysed)
! Observations
!
! Temporary vectors used by the EKF approach

REAL,DIMENSION(KI,NOBSTYPE,NOBSTYPE) :: ZK1
REAL,DIMENSION(KI,NOBSTYPE,NOBSTYPE) :: ZR                 ! covariance matrix of observation errors
REAL,DIMENSION(KI,NPATCH*NVAR,NPATCH*NVAR) :: ZKH, ZKRK
REAL,DIMENSION(KI,NPATCH*NVAR,NPATCH*NVAR) :: ZIDKH
REAL,DIMENSION(KI,NPATCH*NVAR,NPATCH*NVAR) :: ZB           ! background error covariance matrix
REAL,DIMENSION(KI,NPATCH*NVAR,NPATCH*NVAR) :: ZLTM         ! linear tangent matrix for the f'ward model
REAL,DIMENSION(KI,NPATCH*NVAR,NPATCH*NVAR) :: ZQ           ! model error matrix
REAL,DIMENSION(KI,NOBSTYPE,NPATCH*NVAR) :: ZHO             ! Jacobian of observation operator
REAL,DIMENSION(KI,NOBSTYPE,NPATCH*NVAR) :: ZHOWR           ! copy of HO for writing out
REAL,DIMENSION(KI,NPATCH*NVAR,NOBSTYPE) :: ZHOT            ! Transpose of HO
REAL,DIMENSION(KI,NPATCH*NVAR,NOBSTYPE) :: ZGAIN           ! Kalman gain (used explicitly for Ba) 
REAL,DIMENSION(KI,NPATCH,NVAR) :: ZXI                      ! Vector of control variables (at beginning of timestep)
REAL,DIMENSION(KI,NPATCH,NVAR) :: ZEPS                     ! The perturbation amplitude
REAL,DIMENSION(KI,NVAR+1,NOBSTYPE) :: ZYF                  ! Vector of model observations (averaged) 
REAL,DIMENSION(NPATCH*NVAR,NPATCH*NVAR) :: ZIDENT          ! identitiy matrix, used for Ba
REAL,DIMENSION(KI,NOBSTYPE) :: ZYO                         ! vector of observations
REAL,DIMENSION(KI,NOBSTYPE) :: ZX,ZB2,ZP
REAL,DIMENSION(KI,NPATCH*NVAR) :: ZXINCR                   ! Analysis increment
REAL,DIMENSION(KI,NPATCH) :: ZVECT                         ! The analysed variable
REAL,DIMENSION(NDIM_FULL,NOBSTYPE) :: ZYOWR
REAL,DIMENSION(KI) :: ZCOFSWI                              ! dynamic range (Wfc - Wwilt)
REAL,DIMENSION(KI) :: ZSMSAT                               ! saturation  
!
REAL :: ZTIME                      ! current time since start of the run (s)

INTEGER :: IOBSCOUNT
INTEGER :: IYEAR                      ! current year (UTC)
INTEGER :: IMONTH                     ! current month (UTC)
INTEGER :: IDAY                       ! current day (UTC)
INTEGER :: IRESP                      ! return code
INTEGER :: ISTEP                      ! 
INTEGER :: IMYPROC
INTEGER :: IOBS
INTEGER :: ISTAT
!
INTEGER :: I,J,K,JJ,L
!
LOGICAL :: GBEXISTS
!
REAL(KIND=JPRB)                            :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_ISBA_EKF',0,ZHOOK_HANDLE)
!
IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('ASSIM_NATURE_ISBA_EKF: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF
!
IF ( NPRINTLEV>0 ) THEN
  WRITE(*,*)
  WRITE(*,*) '   --------------------------'
  WRITE(*,*) '   |   ENTERING  VARASSIM   |'
  WRITE(*,*) '   --------------------------'
  WRITE(*,*)
ENDIF
!
#ifdef ARO
IF ( MYPROC > 0 ) THEN 
  IMYPROC = MYPROC
ELSE
  IMYPROC = 1
ENDIF
#else
IMYPROC = 1
#endif
!
WRITE(YMYPROC(1:7),'(I7.7)') IMYPROC
!
IF (SUM(NNCV) /= NVAR) THEN
  WRITE(*,*) 'INCONSISTENCY in set-up of CONTROL VARIABLES',SUM(NNCV),NVAR
  CALL ABOR1_SFX('INCONSISTENCY in set-up of CONTROL VARIABLES')
ENDIF
!
IF ( NPRINTLEV > 0 ) WRITE(*,*) 'number of patches =',NPATCH
!
IOBS = NOBSTYPE
!
! Set control variables
DO L = 1,NVAR
  !
  SELECT CASE (TRIM(CVAR(L)))
    CASE("TG1")
      ZXI(:,:,L) = XTG(:,1,:)
    CASE("TG2")
      ZXI(:,:,L) = XTG(:,2,:)
    CASE("WG1")
      ZXI(:,:,L) = XWG(:,1,:)
    CASE("WG2")
      ZXI(:,:,L) = XWG(:,2,:)
    CASE DEFAULT
      CALL ABOR1_SFX("Mapping of "//CVAR(L)//" is not defined in EKF!")
  END SELECT
  !
  IF ( NPRINTLEV > 0 ) WRITE(*,*) CVAR(L),' - initial ',ZXI(1,1,L)
  !
ENDDO
!
!
! Calculate delta for control variables 
DO L = 1,NVAR
  ZEPS(:,:,L) = XTPRT(L) * ZXI(:,:,L)
ENDDO
!
!   Read CLAY fraction to  compute the SWI range (Wfc - Wwilt)
!   (XSIGMA is defined in terms of SWI), need to convert to equivalent v/v
!   using same clay fraction in both layers
!   Read SAND fraction to compute the saturation for conversion of ERS SWI
!
DO I=1,KI
  ZCOFSWI(I) = 0.001 * (89.0467 * ((100.*XCLAY(I,1))**0.3496) - 37.1342*((100.*XCLAY(I,1))**0.5))
  ZSMSAT (I) = 0.001 * (-1.08*100.*XSAND(I,1) + 494.305) 
ENDDO
!
IF ( NPRINTLEV > 0 ) WRITE(*,*) KI,'ki ', KI,'nsize '
!
! ----------------------
! VARASSIM OPTION : LBEV
! ----------------------
!   Calculate the LTM, and evolve B. 
!
IF ( LBEV ) THEN
  !
  ! Set the B input file depending of an existing B was found or not
  YBGFILE = "BGROUNDin."//YMYPROC
  INQUIRE (FILE=TRIM(YBGFILE),EXIST=GBEXISTS)
  !
  IF (GBEXISTS) THEN
    !
    OPEN (unit=111,file=YBGFILE,status='unknown',IOSTAT=ISTAT)
    DO JJ = 1,NVAR   ! control variable (x at previous time step)
      DO L = 1,NVAR
        DO I = 1,KI
          DO J = 1,NPATCH   
            DO K = 1,NPATCH   
              READ (111,*) ZB(I, J+NPATCH*(JJ-1), K+NPATCH*(L-1))
            ENDDO
          ENDDO                 
        ENDDO
      ENDDO
    ENDDO
    CLOSE(111)
    IF ( NPRINTLEV > 0 ) WRITE(*,*) 'read previous B matrix  ==>',ZB(1,1,1),NVAR
    !
  ELSE
    ! Initialization of B 
    ZB(:,:,:) = 0.0
    DO L = 1,NVAR
      DO I = 1,KI
        DO J = 1,NPATCH
          IF ( CVAR(L)=='WG2') &
                  ZB(I,J+NPATCH*(L-1),J+NPATCH*(L-1)) = XSIGMA(L)*XSIGMA(L)*ZCOFSWI(I)*ZCOFSWI(I)
          IF ( CVAR(L)=='WG1') &
                  ZB(I,J+NPATCH*(L-1),J+NPATCH*(L-1)) = XSIGMA(L)*XSIGMA(L)*ZCOFSWI(I)*ZCOFSWI(I)
          IF ( CVAR(L)=='TG2') &
                  ZB(I,J+NPATCH*(L-1),J+NPATCH*(L-1)) = XSIGMA(L)*XSIGMA(L)
          IF ( CVAR(L)=='TG1') &
                  ZB(I,J+NPATCH*(L-1),J+NPATCH*(L-1)) = XSIGMA(L)*XSIGMA(L)       
        ENDDO
      ENDDO
    ENDDO
    IF ( NPRINTLEV > 0 ) WRITE(*,*) 'Initialized B'
    !
  ENDIF

  ! calculate LTM

  DO L = 1,NVAR    ! control variable (x at previous time step)
    DO K = 1,NVAR
      DO I = 1,KI 
        DO J = 1,NPATCH 
          
          IF ( XF(I,J,L+1,K).NE.XUNDEF .AND. XF(I,J,1,K).NE.XUNDEF ) THEN
            !
            ! Jacobian of fwd model
            ZLTM(I, J+NPATCH*(L-1), J+NPATCH*(K-1)) = (XF(I,J,L+1,K) - XF(I,J,1,K)) / ZEPS(I,J,L)
            ! impose upper/lower limits 
            !LTM(I,J+NPATCH*(L-1),J+NPATCH*(K-1)) = max(-0.1, LTM(I,J+NPATCH*(L-1),J+NPATCH*(K-1)))
            !LTM(I,J+NPATCH*(L-1),J+NPATCH*(K-1)) = min(1.0, LTM(I,J+NPATCH*(L-1),J+NPATCH*(K-1)))
            !
          ENDIF

        ENDDO
      ENDDO
    ENDDO
  ENDDO
  !
  IF ( NPRINTLEV > 0 ) WRITE(*,*) 'LTM d(wg2)/d(wg2)', ZLTM(1,1,1)
  !
  ! evolve B 
  !
  DO I=1,KI
    ZB(I,:,:) = MATMUL(ZLTM(I,:,:),MATMUL(ZB(I,:,:),TRANSPOSE(ZLTM(I,:,:))))     
  ENDDO
  !
  ! write out the LTM for the forward model
  DO L=1,NVAR
    DO K=1,NVAR 
      !
      WRITE(YCHAR,'(I1)') K
      YLFNAME='LTM_del'//TRIM(CVAR(K))//'_del'//TRIM(CVAR(L))//"."//YMYPROC
      !
      OPEN(UNIT=111,FILE=YLFNAME,FORM='FORMATTED',STATUS='UNKNOWN',POSITION='APPEND')
      DO I=1,KI
        DO J=1,NPATCH
          WRITE (111,*) ZLTM(I, J+NPATCH*(L-1), J+NPATCH*(K-1))
        ENDDO
      ENDDO
      !
      CLOSE(111)
      !
    ENDDO
  ENDDO
  !
  ! Write out current B
  IF ( NPRINTLEV > 0 ) THEN
    WRITE(*,*) 'store B matrix after TL evolution ==>',ZB(1,1,1)
    WRITE(*,*) 'writing out B'
  ENDIF
  !
  YBGFILE="BGROUNDout_LBEV."//YMYPROC
  !
  OPEN (unit=111,file=YBGFILE,status='unknown')
  DO L=1,NVAR
    DO K=1,NVAR
      DO I=1,KI
        DO J=1,NPATCH
          DO JJ=1,NPATCH
            WRITE (111,*)  ZB(I, J+NPATCH*(L-1), JJ+NPATCH*(K-1))
          ENDDO
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  CLOSE(111)
  !
  !--------------------------------------------------------------------
  !
  !   Adding model error to background error matrix 
  !
  !--------------------------------------------------------------------
  ! 
  ZQ(:,:,:) = 0.0
  DO L=1,NVAR
    DO I=1,KI
      DO J=1,NPATCH
        !
        IF (TRIM(CVAR(L)) == 'WG2') &
          ZQ(I, J+NPATCH*(L-1), J+NPATCH*(L-1)) = XSCALE_Q*XSCALE_Q*XSIGMA(L)*XSIGMA(L)*ZCOFSWI(I)*ZCOFSWI(I)
        IF (TRIM(CVAR(L)) == 'WG1') &
          ZQ(I, J+NPATCH*(L-1), J+NPATCH*(L-1)) = XSCALE_Q*XSCALE_Q*XSIGMA(L)*XSIGMA(L)*ZCOFSWI(I)*ZCOFSWI(I)
        IF (TRIM(CVAR(L)) == 'TG2') &
          ZQ(I, J+NPATCH*(L-1), J+NPATCH*(L-1)) = XSCALE_Q*XSCALE_Q*XSIGMA(L)*XSIGMA(L)
        IF (TRIM(CVAR(L)) == 'TG1') &
          ZQ(I, J+NPATCH*(L-1), J+NPATCH*(L-1)) = XSCALE_Q*XSCALE_Q*XSIGMA(L)*XSIGMA(L)
        !
      ENDDO
    ENDDO
  ENDDO
  !
  ! B is the forecast matrix - need to add Q
  IF ( NPRINTLEV > 0 ) THEN
    WRITE(*,*) 'B before wg2 wg2 ==> ',sqrt(ZB(1,1,1))/ZCOFSWI(1),ZB(1,1,1)
    WRITE(*,*) 'Q value wg2 wg2 ==> ',sqrt(ZQ(1,1,1))/ZCOFSWI(1),ZQ(1,1,1)
  ENDIF
  !
  ZB = ZB + ZQ
  !
  IF ( NPRINTLEV > 0 ) WRITE(*,*) 'B after wg2 wg2 ==>',sqrt(ZB(1,1,1))/ZCOFSWI(1),ZB(1,1,1)
  !
ELSEIF (LBFIXED) THEN
  !
  DO L=1,NVAR
    DO I=1,KI
      DO J=1,NPATCH
        IF (TRIM(CVAR(L)) == 'WG2') ZB(I,J+NPATCH*(L-1),J+NPATCH*(L-1)) = XSIGMA(L)*XSIGMA(L)*ZCOFSWI(I)*ZCOFSWI(I)
        IF (TRIM(CVAR(L)) == 'WG1') ZB(I,J+NPATCH*(L-1),J+NPATCH*(L-1)) = XSIGMA(L)*XSIGMA(L)*ZCOFSWI(I)*ZCOFSWI(I)
        IF (TRIM(CVAR(L)) == 'TG2') ZB(I,J+NPATCH*(L-1),J+NPATCH*(L-1)) = XSIGMA(L)*XSIGMA(L)
        IF (TRIM(CVAR(L)) == 'TG1') ZB(I,J+NPATCH*(L-1),J+NPATCH*(L-1)) = XSIGMA(L)*XSIGMA(L)       
      ENDDO
    ENDDO
  ENDDO
  !
ELSE
  CALL ABOR1_SFX("LBEV or LBFIXED should be .TRUE.!")
ENDIF
!
! ====================================================================
!
! Analysis
!
! ====================================================================
!
!   Time reinitialization 
IYEAR  = TTIME%TDATE%YEAR
IMONTH = TTIME%TDATE%MONTH
IDAY   = TTIME%TDATE%DAY
ZTIME  = TTIME%TIME
!
! Initialisations
ZB2(:,:)      = 999.0     ! Innovation vector
ZR(:,:,:)     = 0.0       ! Observation error matrix
ZHO(:,:,:)    = 999.0     ! Linearized observation matrix
ZHOWR(:,:,:)  = 999.0  
!
ZYOWR(:,:)    = 999.0     ! Observation vector on the full grid to be written on file
ZYF(:,:,:)    = 0.0       ! Tile averaged simulated observation vector
!
IOBSCOUNT     = 0
!
ZIDENT(:,:) = 0                    ! identity matrix
DO L = 1,NVAR
  DO J = 1,NPATCH
    ZIDENT(J+NPATCH*(L-1),J+NPATCH*(L-1)) = 1.0
  ENDDO
ENDDO
! 
ZTIME = XTSTEP_OUTPUT
! BEGINNING OF TIME LOOP
TIMELOOP : DO ISTEP=1,NBOUTPUT
  ! Update date
  CALL ADD_FORECAST_TO_DATE_SURF(IYEAR, IMONTH, IDAY, ZTIME)
  ZTIME = ZTIME + XTSTEP_OUTPUT
  !
  ZYO(:,1) = PT2M (:)
  ZYO(:,2) = PHU2M(:) 
  !
  IF ( NPRINTLEV > 0 ) WRITE(*,*) 'read in obs: ', ZYO(1,1), ZYO(1,2),IOBS
  !
  ! Mean simulated obs averaged over tiles
  DO I=1,KI
    DO J=1,NPATCH
      IF (XPATCH(I,J) > 0.0) THEN
        ZYF(I,:,:) = ZYF(I,:,:) + XPATCH(I,J)*XF_PATCH(I,J,:,:)
      ENDIF
    ENDDO
  ENDDO
  !
ENDDO TIMELOOP
!
! SET OBSERVATION ERROR 
DO I=1,KI
  DO K=1,NOBSTYPE
    IF (K .LT. 3) THEN
      ZR(I,K,K) = XERROBS(K)*XERROBS(K)
    ELSE
      ! convert R for wg1 from SWI  to abs value
      ZR(I,K,K) = XERROBS(K)*XERROBS(K)*ZCOFSWI(I)*ZCOFSWI(I) 
    ENDIF
  ENDDO
ENDDO
!
! WRITE OUT OBS AND YERROR FOR DIAGNOSTIC PURPOSES
OPEN (unit=111,file='OBSERRORout.'//YMYPROC,status='unknown',IOSTAT=istat)
WRITE (111,*) (ZR(I,:,:), I=1,KI)
CLOSE(111)
!
OPEN (unit=111,file='OBSout.'//YMYPROC,status='unknown',IOSTAT=istat)
WRITE (111,*) (ZYO(I,:), I=1,KI)
CLOSE(111)
!
! Data type selection before assimilation (only if NOBS = NOBSTYPE)
IF ( IOBS == NOBSTYPE ) THEN
  DO I = 1,NOBSTYPE
    IF (NNCO(I) == 0) THEN 
      ZYO (:,I) = 999.0
      IF ( NPRINTLEV > 0 ) WRITE(*,*) 'OBSERVATION TYPE ',COBS(I),' REMOVED'
    ENDIF
  ENDDO
ENDIF
!
IF ( NPRINTLEV > 0 ) WRITE(*,*) 'calculating jacobians',IOBS
!
DO L=1,NVAR
  DO I=1,KI
    DO J=1,NPATCH  
      DO K=1,IOBS
        !
        ZHOWR(I,K,J+NPATCH*(L-1)) = XPATCH(I,J)*(XF_PATCH(I,J,L+1,K) - XF_PATCH(I,J,1,K))/ZEPS(I,J,L) 
        IF(ZYO(I,K) .NE. 999.0) THEN         !if obs available
          ! Jacobian of obs operator
          ZHO(I,K,J+NPATCH*(L-1)) = XPATCH(I,J)*(XF_PATCH(I,J,L+1,K) - XF_PATCH(I,J,1,K))/ZEPS(I,J,L)
          ! impose limits  
          !HO(I,K,J+NPATCH*(L-1)) = max(-0.1, HO(I,K,J+NPATCH*(L-1)))             
          !HO(I,K,J+NPATCH*(L-1)) = min(1.0, HO(I,K,J+NPATCH*(L-1)))
          ! innovation vector
          ZB2(I,K) = ZYO(I,K) - ZYF(I,1,K)                 
          IOBSCOUNT = IOBSCOUNT + 1
        ELSE  !if no obs available
          ! set obs operator and innovation to zero if no obs available
          ZHO(I,K,J+NPATCH*(L-1)) = 0.0
          ZB2(I,K) = 0.0 
        ENDIF
      ENDDO
    ENDDO
  ENDDO
ENDDO
!
! *** Write innovations in ASCII file ***
!
OPEN (unit=111,file='INNOV.'//YMYPROC,status='unknown',IOSTAT=istat)
DO I=1,KI
  WRITE(111,*) (ZB2(I,K),K=1,IOBS)
ENDDO
CLOSE(UNIT=111)
!
!-----------------------------------------------------
!
!            ******  SOIL ANALYSIS *******
!
!-----------------------------------------------------
IF ( NPRINTLEV > 0 ) WRITE(*,*) 'PERFORMING ANALYSIS'
!
DO I=1,KI
  !
  ZHOT(I,:,:) = TRANSPOSE(ZHO(I,:,:))
  ZK1 (I,:,:) = MATMUL(ZHO(I,:,:),MATMUL(ZB(I,:,:),ZHOT(I,:,:))) + ZR(I,:,:)
  CALL CHOLDC(NOBSTYPE,ZK1(I,:,:),ZP(I,:))                         ! Cholesky decomposition (1)
  CALL CHOLSL(NOBSTYPE,ZK1(I,:,:),ZP(I,:),ZB2(I,:),ZX(I,:))         ! Cholesky decomposition (2)
  ZXINCR(I,:) = MATMUL(ZB(I,:,:),MATMUL(ZHOT(I,:,:),ZX(I,:))) 
  !
  DO L=1,NVAR
    DO J=1,NPATCH
      ! Update the modified values
      ZXI(I,J,L) = ZXI(I,J,L) + ZXINCR(I,J+NPATCH*(L-1))
      ! A possible sanity check
      SELECT CASE (TRIM(CVAR(L)))
        CASE("TG1")
        CASE("TG2")
        CASE("WG1")
          ! Preserve positive values
          !XI(I,J,L) = MAX(0.,(XI(I,J,L)))
        CASE("WG2")
          ! Preserve positive values
          !XI(I,J,L) = MAX(0.,(XI(I,J,L)))
        CASE DEFAULT
          CALL ABOR1_SFX("Mapping of "//TRIM(CVAR(L))//" is not defined in EKF!")
      END SELECT
      ! For no only warn if we have negative values.
      IF ( ZXI(I,J,L) < 0. ) WRITE(*,*) "WARNING X<0. for ",I,J," for variable ",TRIM(CVAR(L))
    ENDDO
  ENDDO
  !
ENDDO
!
! Write analysis results and increments in ASCII file
OPEN (unit=111,file='ANAL_INCR.'//YMYPROC,status='unknown',IOSTAT=istat)
DO I=1,KI
  DO J=1,NPATCH
    WRITE(111,*) (ZXI(I,J,L),L=1,NVAR), (ZXINCR(I,J+NPATCH*(L-1)),L=1,NVAR)
  ENDDO
ENDDO
CLOSE(UNIT=111)
!
! **********************************************
DO L=1,NVAR
  !
  WRITE(*,*) 'Sum and mean increments for ',TRIM(CVAR(l)),'=',SUM(ZXINCR(:,L)),SUM(ZXINCR(:,L))/REAL(KI)
  !
  ! Update the modified values
  SELECT CASE (TRIM(CVAR(L)))
    CASE("TG1")
      XTG(:,1,:) = ZXI(:,:,L)
    CASE("TG2")
      XTG(:,2,:) = ZXI(:,:,L)
    CASE("WG1")
      XWG(:,1,:) = ZXI(:,:,L)
    CASE("WG2")
      XWG(:,2,:) = ZXI(:,:,L)
    CASE DEFAULT
      CALL ABOR1_SFX("Mapping of "//TRIM(CVAR(L))//" is not defined in EKF!")
  END SELECT
ENDDO
!
! Analysis of B (for use in next cycle)
! Ba = (I-KH)Bf(I-KH)t+KRKt
! K = BfHt{K1}**-1
! K1 needs PATCH dim added
DO I=1,KI 
  !
  ! K1 = (R+H.B.HT) (calculate inverse -> output goes to K1)
  CALL INVERSE_MATRIX(NOBSTYPE,ZK1(I,:,:),ZP(I,:))
  ZGAIN(I,:,:) = MATMUL(ZB(I,:,:),MATMUL(ZHOT(I,:,:),ZK1(I,:,:)))
  ZKH(I,:,:) = MATMUL(ZGAIN(I,:,:),ZHO(I,:,:))
  ZIDKH(I,:,:) = ZIDENT(:,:) - ZKH(I,:,:) 
  ZKRK(I,:,:) = MATMUL(ZGAIN(I,:,:),MATMUL(ZR(I,:,:),TRANSPOSE(ZGAIN(I,:,:))))
  IF (.NOT.LBFIXED)  ZB(I,:,:) = MATMUL(ZIDKH(I,:,:),MATMUL(ZB(I,:,:),TRANSPOSE(ZIDKH(I,:,:))))+ZKRK(I,:,:)
  !
ENDDO
!
! Write out analysed B (for use in next cycle)
YBGFILE = "BGROUNDout_ASSIM."//YMYPROC
OPEN (unit=111,file=YBGFILE,status='unknown',IOSTAT=istat)
DO L=1,NVAR
  DO K=1,NVAR
    DO I=1,KI
      DO J=1,NPATCH
        DO JJ=1,NPATCH
          WRITE (111,*) ZB(I,J+NPATCH*(L-1),JJ+NPATCH*(K-1))
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDDO
CLOSE(111)
!
! **** Write out the observation operator + Gain matrix ****
DO L=1,NVAR
  DO K=1, NOBSTYPE
    WRITE(YCHAR,'(I1)') K
    YFNAME='HO_'//CVAR(L)//'_v'//YCHAR
    OPEN(UNIT=111,FILE=YFNAME,FORM='FORMATTED',STATUS='NEW',IOSTAT=ISTAT)
    DO I=1,KI
      DO J=1,NPATCH
        WRITE(111,*) ZHOWR(I,K,J+NPATCH*(L-1)),ZGAIN(I,J+NPATCH*(L-1),K)
      ENDDO
    ENDDO
    CLOSE(111)
  ENDDO
ENDDO
!
IF ( NPRINTLEV > 0 ) THEN
  WRITE(*,*)
  WRITE(*,*) '   ---------------------------------------'
  WRITE(*,*) '   |   EXITING VARASSIM AFTER ANALYSIS   |'
  WRITE(*,*) '   ---------------------------------------'
  WRITE(*,*)
  WRITE(*,*) 'Number of assimilated observations =',IOBSCOUNT
  WRITE(*,*)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_ISBA_EKF',1,ZHOOK_HANDLE)
!
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Subroutines of Linear Algebra  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
!
DO I=1,N
  DO J=I+1,N
    A(I,J) =0.0
  ENDDO
ENDDO
!
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

