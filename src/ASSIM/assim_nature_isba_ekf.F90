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
                               XSCALE_Q, NPRINTLEV, CVAR, XTPRT, XSIGMA, CBIO,     &
                               XF_PATCH, XF, COBS, XSCALE_QLAI,  LOBSFILE, XALPH,  &
                               NECHGU, NBOUTPUT, XEPS
! 
USE MODD_SURF_PAR,      ONLY : XUNDEF
!
USE MODD_SURF_ATM_n,    ONLY : NDIM_FULL
USE MODD_ISBA_n,        ONLY : XTG, XWG, XWGI, XLAI, XSAND, XCLAY, TTIME, &
                               NPATCH, XPATCH, XBIOMASS, XRESP_BIOMASS
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
USE MODI_GET_FILE_NAME
!
! -----------------------------------------------------------
!
IMPLICIT NONE
!
CHARACTER(LEN=6),    INTENT(IN)           :: HPROGRAM     ! program calling surf. schemes
INTEGER,             INTENT(IN)           :: KI
REAL, DIMENSION(:), INTENT(IN)            :: PT2M
REAL, DIMENSION(:), INTENT(IN)            :: PHU2M
CHARACTER(LEN=2),    INTENT(IN)           :: HTEST        ! must be equal to 'OK'
!
!    Declarations of local variables
!
 CHARACTER(LEN=200) :: YMFILE     ! Name of the observation, perturbed or reference file!  
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

REAL,DIMENSION(KI,NOBSTYPE*NBOUTPUT,NOBSTYPE*NBOUTPUT) :: ZK1
REAL,DIMENSION(KI,NOBSTYPE*NBOUTPUT,NOBSTYPE*NBOUTPUT) :: ZR                 ! covariance matrix of observation errors
REAL,DIMENSION(KI,NPATCH*NVAR,NPATCH*NVAR) :: ZKH, ZKRK
REAL,DIMENSION(KI,NPATCH*NVAR,NPATCH*NVAR) :: ZIDKH
REAL,DIMENSION(KI,NPATCH*NVAR,NPATCH*NVAR) :: ZB           ! background error covariance matrix
REAL,DIMENSION(KI,NPATCH*NVAR,NPATCH*NVAR) :: ZLTM         ! linear tangent matrix for the f'ward model
REAL,DIMENSION(KI,NPATCH*NVAR,NPATCH*NVAR) :: ZQ           ! model error matrix
REAL,DIMENSION(KI,NOBSTYPE*NBOUTPUT,NPATCH*NVAR) :: ZHO             ! Jacobian of observation operator
REAL,DIMENSION(KI,NOBSTYPE*NBOUTPUT,NPATCH*NVAR) :: ZHOWR           ! copy of HO for writing out
REAL,DIMENSION(KI,NPATCH*NVAR,NOBSTYPE*NBOUTPUT) :: ZHOT            ! Transpose of HO
REAL,DIMENSION(KI,NPATCH*NVAR,NOBSTYPE*NBOUTPUT) :: ZGAIN           ! Kalman gain (used explicitly for Ba) 
REAL,DIMENSION(KI,NPATCH,NVAR) :: ZEPS                     ! The perturbation amplitude
REAL,DIMENSION(KI,NVAR+1,NOBSTYPE*NBOUTPUT) :: ZYF                  ! Vector of model observations (averaged) 
REAL,DIMENSION(NPATCH*NVAR,NPATCH*NVAR) :: ZIDENT          ! identitiy matrix, used for Ba
REAL,DIMENSION(KI,NOBSTYPE*NBOUTPUT) :: ZYO                         ! vector of observations
REAL,DIMENSION(KI,NOBSTYPE*NBOUTPUT) :: ZX,ZB2,ZP
REAL,DIMENSION(KI,NPATCH*NVAR) :: ZXINCR                   ! Analysis increment
REAL,DIMENSION(KI,NPATCH) :: ZVECT                         ! The analysed variable
REAL,DIMENSION(KI,NPATCH) :: ZBIO_PASS
REAL,DIMENSION(KI,NPATCH) :: ZBIO_OUT
REAL,DIMENSION(KI) :: ZCOFSWI                              ! dynamic range (Wfc - Wwilt)
REAL,DIMENSION(KI) :: ZSMSAT                               ! saturation  
REAL,DIMENSION(KI) :: ZWILT
!
REAL,DIMENSION(NPATCH) :: ZVLAIMIN
!
REAL :: ZTIME                      ! current time since start of the run (s)

INTEGER :: IOBSCOUNT
INTEGER :: IYEAR                      ! current year (UTC)
INTEGER :: IMONTH                     ! current month (UTC)
INTEGER :: IDAY                       ! current day (UTC)
INTEGER :: IHOUR
INTEGER :: IRESP                      ! return code
INTEGER :: ISTEP                      ! 
INTEGER :: IMYPROC
INTEGER :: INOBS, IOBS
INTEGER :: ISTAT
!
INTEGER :: I,J,K,JJ,L,K1,L1
!
LOGICAL :: GBEXISTS
!
REAL(KIND=JPRB)                            :: ZHOOK_HANDLE
!
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
IF ( NPRINTLEV > 0 ) WRITE(*,*) 'number of patches =',NPATCH
!
! Set control variables
DO L = 1,NVAR
  IF ( NPRINTLEV > 0 ) WRITE(*,*) CVAR(L),' - initial ',XF(1,1,1,L)
ENDDO
!
SELECT CASE (TRIM(CBIO))
  CASE("BIOMA1","BIOMASS1")
    ZBIO_PASS(:,:) = XBIOMASS(:,1,:)
  CASE("BIOMA2","BIOMASS2")
    ZBIO_PASS(:,:) = XBIOMASS(:,2,:)
  CASE("RESPI1","RESP_BIOM1")
    ZBIO_PASS(:,:) = XRESP_BIOMASS(:,1,:)
  CASE("RESPI2","RESP_BIOM2")
    ZBIO_PASS(:,:) = XRESP_BIOMASS(:,2,:)
  CASE("LAI")
    ZBIO_PASS(:,:) = XLAI(:,:)
  CASE DEFAULT
    CALL ABOR1_SFX("Mapping of "//CBIO//" is not defined in EKF!")
END SELECT
!
!
!   Read CLAY fraction to  compute the SWI range (Wfc - Wwilt)
!   (XSIGMA is defined in terms of SWI), need to convert to equivalent v/v
!   using same clay fraction in both layers
!   Read SAND fraction to compute the saturation for conversion of ERS SWI
!
DO I=1,KI
  ZCOFSWI(I) = 0.001 * (89.0467 * ((100.*XCLAY(I,1))**0.3496) - 37.1342*((100.*XCLAY(I,1))**0.5))
  ZSMSAT (I) = 0.001 * (-1.08*100.*XSAND(I,1) + 494.305)
  ZWILT  (I) = 0.001 * 37.1342 * ((100.*XCLAY(I,1))**0.5) 
ENDDO
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
          L1 = J+NPATCH*(L-1)
          IF ( TRIM(CVAR(L))=='WG2' .OR. TRIM(CVAR(L))=='WG1' ) &
            ZB(I,L1,L1) = XSIGMA(L)*XSIGMA(L)*ZCOFSWI(I)*ZCOFSWI(I)
          IF ( TRIM(CVAR(L))=='TG2' .OR. TRIM(CVAR(L))=='TG1' ) &
            ZB(I,L1,L1) = XSIGMA(L)*XSIGMA(L)
          IF ( TRIM(CVAR(L))=='LAI' ) &
            ZB(I,L1,L1) = XSIGMA(L)*XSIGMA(L)*XLAI(I,J)*XLAI(I,J) 
        ENDDO
      ENDDO
    ENDDO
    IF ( NPRINTLEV > 0 ) WRITE(*,*) 'Initialized B'
    !
  ENDIF

  ! calculate LTM
  ZLTM(:,:,:) = 0.0
  DO L = 1,NVAR    ! control variable (x at previous time step)
    DO K = 1,NVAR
      DO I = 1,KI 
        DO J = 1,NPATCH 
          
          IF ( XPATCH(I,J)>0.0 .AND. XF(I,J,L+1,K).NE.XUNDEF .AND. XF(I,J,1,K).NE.XUNDEF ) THEN
            !
            L1 = J + NPATCH*(L-1)
            K1 = J + NPATCH*(K-1)
            ! Jacobian of fwd model
            ZLTM(I,L1,K1) = ( XF(I,J,L+1,K) - XF(I,J,1,K) ) / XEPS(I,J,L)
            ! impose upper/lower limits 
            ZLTM(I,L1,K1) = MAX(-0.1, ZLTM(I,L1,K1))
            ZLTM(I,L1,K1) = MIN( 1.0, ZLTM(I,L1,K1))
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
  ! Write out current B
  YBGFILE="BGROUNDout_LBEV."//YMYPROC
  OPEN (unit=112,file=YBGFILE,status='unknown')
  !
  DO L=1,NVAR
    DO K=1,NVAR
      !
      WRITE(YCHAR,'(I1)') K
      YLFNAME='LTM_del'//TRIM(CVAR(K))//'_del'//TRIM(CVAR(L))//"."//YMYPROC
      !
      OPEN(UNIT=111,FILE=YLFNAME,FORM='FORMATTED',STATUS='UNKNOWN',POSITION='APPEND')
      DO I=1,KI
        DO J=1,NPATCH
          !
          WRITE (111,*) ZLTM(I,J+NPATCH*(L-1),J+NPATCH*(K-1))
          !
          DO JJ=1,NPATCH
            WRITE (112,*)  ZB(I,J+NPATCH*(L-1),JJ+NPATCH*(K-1))
          ENDDO
          !     
        ENDDO
      ENDDO
      !
      CLOSE(111)
      !
    ENDDO
  ENDDO
  !
  CLOSE(112)
  !
  IF ( NPRINTLEV > 0 ) THEN
    WRITE(*,*) 'store B matrix after TL evolution ==>',ZB(1,1,1)
    WRITE(*,*) 'writing out B'
  ENDIF
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
        L1 = J+NPATCH*(L-1)
        !
        IF (TRIM(CVAR(L)) == 'WG2' .OR. TRIM(CVAR(L)) == 'WG1') &
          ZQ(I,L1,L1) = XSCALE_Q*XSCALE_Q*XSIGMA(L)*XSIGMA(L)*ZCOFSWI(I)*ZCOFSWI(I)
        IF (TRIM(CVAR(L)) == 'TG2' .OR. TRIM(CVAR(L)) == 'TG1') &
          ZQ(I,L1,L1) = XSCALE_Q*XSCALE_Q*XSIGMA(L)*XSIGMA(L)
        IF (TRIM(CVAR(L)) == 'LAI') &
          ZQ(I,L1,L1) = XSCALE_QLAI*XSCALE_QLAI*XSIGMA(L)*XSIGMA(L)
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
        L1 = J+NPATCH*(L-1)
        IF ( TRIM(CVAR(L)) == 'WG2' .OR. TRIM(CVAR(L)) == 'WG1' ) &
          ZB(I,L1,L1) = XSIGMA(L)*XSIGMA(L)*ZCOFSWI(I)*ZCOFSWI(I)
        IF ( TRIM(CVAR(L)) == 'TG2' .OR. TRIM(CVAR(L)) == 'TG1' ) &
          ZB(I,L1,L1) = XSIGMA(L)*XSIGMA(L)
        IF (TRIM(CVAR(L)) == 'LAI') &
          ZB(I,L1,L1) = XSIGMA(L)*XSIGMA(L)*XLAI(I,J)*XLAI(I,J)      
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
!
! Initialisations
ZB2(:,:)      = XUNDEF    ! Innovation vector
ZR(:,:,:)     = 0.0       ! Observation error matrix
ZHO(:,:,:)    = XUNDEF    ! Linearized observation matrix
ZHOWR(:,:,:)  = XUNDEF 
!
ZX(:,:) = 0.0
!
ZYF(:,:,:)    = 0.0       ! Tile averaged simulated observation vector
!
ZIDENT(:,:) = 0                    ! identity matrix
DO L = 1,NVAR
  DO J = 1,NPATCH
    ZIDENT(J+NPATCH*(L-1),J+NPATCH*(L-1)) = 1.0
  ENDDO
ENDDO
!
!
INOBS = 0
IHOUR = 0
ZTIME = FLOAT(NECHGU) * 3600.
! BEGINNING OF TIME LOOP
TIMELOOP : DO ISTEP=1,NBOUTPUT
  ! Update date
  CALL ADD_FORECAST_TO_DATE_SURF(IYEAR, IMONTH, IDAY, ZTIME)
  ZTIME = ZTIME + FLOAT(NECHGU) * 3600.
  IHOUR = IHOUR + NECHGU
  !
  IF ( LOBSFILE ) THEN
    !
    INOBS = 0
    YMFILE = 'CANARI_NATURE_'
    CALL GET_FILE_NAME(IYEAR,IMONTH,IDAY,IHOUR,YMFILE)
    OPEN(UNIT=55,FILE=TRIM(YMFILE)//".DAT",FORM='FORMATTED',STATUS='OLD',IOSTAT=ISTAT) 
    IF ( ISTAT==0 ) THEN
      !   If it exists, read observations
      DO I = 1,KI
        READ (55,*)  (ZYO(I,INOBS+J),J=1,NOBSTYPE)
      ENDDO
      INOBS = INOBS + NOBSTYPE      
      IF ( NPRINTLEV > 0 ) WRITE(*,*) 'read in obs: ', ZYO(1,:), INOBS
      CLOSE(55)
    ENDIF
    !
    IF ( INOBS == 0 ) THEN
      IF ( NPRINTLEV > 0 ) WRITE(*,*) 'No observations read for LAI in OBS file - stop'
      CALL ABOR1_SFX("ASSIM_NATURE_ISBA_EKF: No observations read for LAI in OBS file - stop")
    ENDIF
    !
  ELSE
    !
    INOBS = NOBSTYPE
    DO IOBS = 1,INOBS
      SELECT CASE (TRIM(COBS(IOBS)))   
        CASE("T2M") 
          ZYO(:,IOBS) = PT2M(:)
        CASE("HU2M")   
          ZYO(:,IOBS) = PHU2M(:) 
        CASE("WG1","LAI")  
          CALL ABOR1_SFX("Mapping of "//COBS(IOBS)//" is not defined in ASSIM_NATURE_ISBA_EKF!")
      END SELECT                 
    ENDDO
    !
  ENDIF
  !
ENDDO TIMELOOP
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
IF ( NPRINTLEV > 0 ) WRITE(*,*) 'read in sim obs yf', ZYF(1,:,1)
!
!
! SET OBSERVATION ERROR 
DO I=1,KI
  DO K=1,NOBSTYPE
    IF ( COBS(K) .EQ. "LAI" ) THEN
      ZR(I,K,K) = XERROBS(K)*XERROBS(K) * ZYO(I,K)*ZYO(I,K)
    ELSEIF (COBS(K) .EQ. "WG1") THEN
      ! convert R for wg1 from SWI  to abs value
      ZR(I,K,K) = XERROBS(K)*XERROBS(K) * ZCOFSWI(I)*ZCOFSWI(I)
    ELSE
      ZR(I,K,K) = XERROBS(K)*XERROBS(K)
    ENDIF
  ENDDO
ENDDO
!
! WRITE OUT OBS AND YERROR FOR DIAGNOSTIC PURPOSES
OPEN (UNIT=111,FILE='OBSERRORout.'//YMYPROC,STATUS='unknown',IOSTAT=ISTAT)
WRITE (111,*) (ZR(I,:,:), I=1,KI)
CLOSE(111)
!
OPEN (UNIT=111,FILE='OBSout.'//YMYPROC,STATUS='unknown',IOSTAT=ISTAT)
WRITE (111,*) (ZYO(I,:), I=1,KI)
CLOSE(111)
!
!
IOBS = NOBSTYPE
! Data type selection before assimilation 
DO I = 1,NOBSTYPE
  IF ( NNCO(I)==0 ) THEN 
    ZYO (:,I) = XUNDEF
    IF ( NPRINTLEV > 0 ) WRITE(*,*) 'OBSERVATION TYPE ',COBS(I),' REMOVED'
  ENDIF
ENDDO
!
DO I = 1,KI
  IF ( MINVAL(XWGI(I,1,:))>0. ) THEN
    ZYO (I,:) = XUNDEF
    IF ( NPRINTLEV > 0 ) WRITE(*,*) 'OBSERVATION FOR POINT ',I,' REMOVED'
  ENDIF
ENDDO
!
IF ( NPRINTLEV > 0 ) WRITE(*,*) 'calculating jacobians',IOBS
!
IOBSCOUNT = 0
DO L=1,NVAR
  DO I=1,KI
    DO J=1,NPATCH
      L1 = J + NPATCH*(L-1)
      DO K=1,IOBS
        !
        ZHOWR(I,K,L1) = XPATCH(I,J)*(XF_PATCH(I,J,L+1,K) - XF_PATCH(I,J,1,K))/XEPS(I,J,L) 
        IF( ZYO(I,K).NE.XUNDEF ) THEN         !if obs available
          ! Jacobian of obs operator
          ZHO(I,K,L1) = ZHOWR(I,K,J+NPATCH*(L-1))
          ! impose limits  
          ZHO(I,K,L1) = MAX(-0.1, ZHO(I,K,L1))             
          ZHO(I,K,L1) = MIN( 1.0, ZHO(I,K,L1))
          ! innovation vector
          ZB2(I,K) = ZYO(I,K) - ZYF(I,1,K)
          IOBSCOUNT = IOBSCOUNT + 1
        ELSE  !if no obs available
          ! set obs operator and innovation to zero if no obs available
          ZHO(I,K,L1) = 0.0
          ZB2(I,K) = 0.0 
        ENDIF
      ENDDO
    ENDDO
  ENDDO
ENDDO
!
!-----------------------------------------------------
!
!            ******  SOIL ANALYSIS *******
!
!-----------------------------------------------------
IF ( NPRINTLEV > 0 ) WRITE(*,*) 'PERFORMING ANALYSIS'
!
ZVLAIMIN=(/0.3,0.3,0.3,0.3,1.0,1.0,0.3,0.3,0.3,0.3,0.3,0.3/)
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
      L1 = J+NPATCH*(L-1)
      ! Update the modified values
      IF ( TRIM(CVAR(L))=="LAI" ) THEN
        ZXINCR(I,L1) = MAX( ZXINCR(I,L1), ZVLAIMIN(J)-XF(I,J,1,L) )
        ZBIO_PASS(I,J) = ZBIO_PASS(I,J) + ZXINCR(I,L1)*XALPH(J)
      ELSEIF ( XF(I,J,1,L)+ZXINCR(I,L1)<0. ) THEN
        ZXINCR(I,L1) = 0.         
      ENDIF
      XF(I,J,1,L) = XF(I,J,1,L) + ZXINCR(I,L1)
      ! For no only warn if we have negative values.
      IF ( XF(I,J,1,L) < 0. ) WRITE(*,*) "WARNING X<0. for ",I,J," for variable ",TRIM(CVAR(L))
    ENDDO
  ENDDO
  !
ENDDO
!
!
! **********************************************
DO L=1,NVAR
  !
  WRITE(*,*) 'Sum and mean increments for ',TRIM(CVAR(l)),'=',SUM(ZXINCR(:,L)),SUM(ZXINCR(:,L))/REAL(KI)
  !
  ! Update the modified values
  SELECT CASE (TRIM(CVAR(L)))
    CASE("TG1")
      XTG(:,1,:) = XF(:,:,1,L)
    CASE("TG2")
      XTG(:,2,:) = XF(:,:,1,L)
    CASE("WG1")
      XWG(:,1,:) = XF(:,:,1,L)
    CASE("WG2")
      XWG(:,2,:) = XF(:,:,1,L)
    CASE("LAI") 
      XLAI(:,:) = XF(:,:,1,L)
      SELECT CASE (TRIM(CBIO))
        CASE("BIOMA1","BIOMASS1")
          XBIOMASS(:,1,:) = ZBIO_PASS(:,:)
        CASE("BIOMA2","BIOMASS2")
          XBIOMASS(:,2,:) = ZBIO_PASS(:,:)
        CASE("RESPI1","RESP_BIOM1")
          XRESP_BIOMASS(:,1,:) = ZBIO_PASS(:,:)
        CASE("RESPI2","RESP_BIOM2")
          XRESP_BIOMASS(:,2,:) = ZBIO_PASS(:,:)
        CASE("LAI")
          XLAI(:,:) = ZBIO_PASS(:,:)
        CASE DEFAULT
          CALL ABOR1_SFX("Mapping of "//CBIO//" is not defined in EKF!")
      END SELECT
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
!
! *** Write innovations in ASCII file ***
!
OPEN (unit=111,file='INNOV.'//YMYPROC,status='unknown',IOSTAT=ISTAT)
DO I=1,KI
  WRITE(111,*) (ZB2(I,K),K=1,IOBS)
ENDDO
CLOSE(UNIT=111)
!
! Write analysis results and increments in ASCII file
OPEN (unit=111,file='ANAL_INCR.'//YMYPROC,status='unknown',IOSTAT=istat)
DO I=1,KI
  DO J=1,NPATCH
    WRITE(111,*) (XF(I,J,1,L),L=1,NVAR), (ZXINCR(I,J+NPATCH*(L-1)),L=1,NVAR)
  ENDDO
ENDDO
CLOSE(UNIT=111)
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
    OPEN(UNIT=111,FILE=YFNAME,FORM='FORMATTED',STATUS='UNKNOWN',IOSTAT=ISTAT)
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
  IOBSCOUNT = IOBSCOUNT / NPATCH / NVAR
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

SUBROUTINE INVERSE_MATRIX(KN,PA,PP)
!--------------------------------------------------------
!
! Explicit inversed matrix after Cholesky decomposition
!
!--------------------------------------------------------
IMPLICIT NONE
INTEGER, INTENT(IN) :: KN
REAL, DIMENSION (KN,KN), INTENT(INOUT) :: PA
REAL, DIMENSION (KN),    INTENT(IN)    :: PP
REAL :: ZSUM
INTEGER :: I, J, K
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_ISBA_EKF:INVERSE_MATRIX',0,ZHOOK_HANDLE)
!
DO I = 1,KN
  PA(I,I)=1./PP(I)
  DO J = I+1,KN
    ZSUM = 0.
    DO K = I,J-1
      ZSUM = ZSUM - PA(J,K)*PA(K,I)
    ENDDO
    PA(J,I) = ZSUM/PP(J)
  ENDDO
ENDDO  
!
DO I = 1,KN
  DO J = I+1,KN
    PA(I,J) = 0.0
  ENDDO
ENDDO
!
PA = MATMUL(TRANSPOSE(PA),PA)
!
IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_ISBA_EKF:INVERSE_MATRIX',1,ZHOOK_HANDLE)
!
END SUBROUTINE INVERSE_MATRIX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CHOLDC(KN,PA,PP)
IMPLICIT NONE
INTEGER, INTENT(IN) :: KN
REAL, DIMENSION(KN,KN), INTENT(INOUT) :: PA
REAL, DIMENSION(KN), INTENT(OUT) :: PP
INTEGER :: I
REAL :: ZSUM
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_ISBA_EKF:CHOLDC',0,ZHOOK_HANDLE)
!
DO I = 1,KN
  ZSUM = PA(I,I)- DOT_PRODUCT(PA(I,1:I-1),PA(I,1:I-1))
  IF ( ZSUM<=0.0 ) CALL ABOR1_SFX('ASSIM_NATURE_ISBA_EKF: ERROR IN CHOLDC')
  PP(I) = SQRT(ZSUM) 
  PA(I+1:KN,I) = ( PA(I,I+1:KN) - MATMUL(PA(I+1:KN,1:I-1),PA(I,1:I-1)) )/PP(I)
ENDDO
!
IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_ISBA_EKF:CHOLDC',1,ZHOOK_HANDLE)
!
END SUBROUTINE CHOLDC 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE CHOLSL(KN,PA,PP,PB,PX)
IMPLICIT NONE
INTEGER, INTENT(IN) :: KN
REAL, DIMENSION (KN,KN), INTENT(IN) :: PA
REAL, DIMENSION (KN), INTENT(IN)    :: PP,PB
REAL, DIMENSION (KN), INTENT(INOUT) :: PX
INTEGER :: I
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_ISBA_EKF:CHOLSL',0,ZHOOK_HANDLE)
!
DO I = 1,KN
  PX(I) = (PB(I) - DOT_PRODUCT(PA(I,1:I-1),PX(1:I-1)))/PP(I)
ENDDO
DO I = KN,1,-1
  PX(I) = (PX(I) - DOT_PRODUCT(PA(I+1:KN,I),PX(I+1:KN)))/PP(I)
ENDDO
!
IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_ISBA_EKF:CHOLSL',1,ZHOOK_HANDLE)
!
END SUBROUTINE CHOLSL
!
!
END SUBROUTINE ASSIM_NATURE_ISBA_EKF
!

