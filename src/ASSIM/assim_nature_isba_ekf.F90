SUBROUTINE ASSIM_NATURE_ISBA_EKF(HPROGRAM, KI, PT2M, PHU2M, HTEST)

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
                               XSCALE_Q, NPRINTLEV, CVAR, XSIGMA, CBIO, XI,        &
                               XF_PATCH, XF, COBS, XSCALE_QLAI,  LOBSFILE, XALPH,  &
                               NECHGU, NBOUTPUT, XTPRT, XLAI_PASS, XBIO_PASS
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
CHARACTER(LEN=6),   INTENT(IN) :: HPROGRAM     ! program calling surf. schemes
INTEGER,            INTENT(IN) :: KI
REAL, DIMENSION(:), INTENT(IN) :: PT2M
REAL, DIMENSION(:), INTENT(IN) :: PHU2M
CHARACTER(LEN=2),   INTENT(IN) :: HTEST        ! must be equal to 'OK'
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
REAL,DIMENSION(KI) :: ZCOFSWI                     ! dynamic range (Wfc - Wwilt)
!REAL,DIMENSION(KI) :: ZSMSAT                      ! saturation  
!REAL,DIMENSION(KI) :: ZWILT
!
REAL,DIMENSION(KI,NPATCH,NVAR) :: ZCOEF
REAL,DIMENSION(KI,NPATCH,NVAR) :: ZEPS            ! The perturbation amplitude
!
REAL,DIMENSION(KI,NOBSTYPE*NBOUTPUT) :: ZYO       ! vector of observations
REAL,DIMENSION(NVAR+1,NOBSTYPE) :: ZYF            ! Vector of model observations (averaged) 
!
REAL,DIMENSION(KI*NPATCH*NVAR*NPATCH*NVAR) :: ZBLONG
REAL,DIMENSION(KI,NPATCH*NVAR,NPATCH*NVAR) :: ZB           ! background error covariance matrix
!
REAL,DIMENSION(NPATCH*NVAR,NPATCH*NVAR) :: ZLTM         ! linear tangent matrix for the f'ward model
REAL,DIMENSION(NPATCH*NVAR,NPATCH*NVAR) :: ZQ           ! model error matrix
!
REAL,DIMENSION(NOBSTYPE*NBOUTPUT,NPATCH*NVAR) :: ZHO             ! Jacobian of observation operator
REAL,DIMENSION(NPATCH*NVAR,NOBSTYPE*NBOUTPUT) :: ZHOT            ! Transpose of HO
REAL,DIMENSION(NPATCH*NVAR,NOBSTYPE*NBOUTPUT) :: ZGAIN           ! Kalman gain (used explicitly for Ba) 
!
REAL,DIMENSION(NOBSTYPE*NBOUTPUT,NOBSTYPE*NBOUTPUT) :: ZR        ! covariance matrix of observation errors
REAL,DIMENSION(NOBSTYPE*NBOUTPUT,NOBSTYPE*NBOUTPUT) :: ZK1
REAL,DIMENSION(NOBSTYPE*NBOUTPUT) :: ZX,ZB2,ZP
!
REAL,DIMENSION(NPATCH*NVAR,NPATCH*NVAR) :: ZKRK
REAL,DIMENSION(NPATCH*NVAR,NPATCH*NVAR) :: ZIDKH
REAL,DIMENSION(NPATCH*NVAR,NPATCH*NVAR) :: ZIDENT          ! identitiy matrix, used for Ba
REAL,DIMENSION(NPATCH*NVAR) :: ZXINCR
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
INTEGER :: ISTAT, ICPT, IUNIT
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
!############################# BEGINNING ###############################
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
!############################# INITIALISATIONS ###############################
!
!   Read CLAY fraction to  compute the SWI range (Wfc - Wwilt)
!   (XSIGMA is defined in terms of SWI), need to convert to equivalent v/v
!   using same clay fraction in both layers
!   Read SAND fraction to compute the saturation for conversion of ERS SWI
!
DO I=1,KI
  ZCOFSWI(I) = 0.001 * (89.0467 * ((100.*XCLAY(I,1))**0.3496) - 37.1342*((100.*XCLAY(I,1))**0.5))
  !ZSMSAT (I) = 0.001 * (-1.08*100.*XSAND(I,1) + 494.305)
  !ZWILT  (I) = 0.001 * 37.1342 * ((100.*XCLAY(I,1))**0.5) 
ENDDO
!
! Set control variables
ZIDENT(:,:) = 0.                   ! identity matrix
DO L = 1,NVAR
  !
  DO J = 1,NPATCH
    ZIDENT(J+NPATCH*(L-1),J+NPATCH*(L-1)) = 1.0
  ENDDO
  !
  WHERE ( XI(:,:,L)/=XUNDEF )
    ZEPS(:,:,L) = XTPRT(L) * XI(:,:,L)
  ELSEWHERE
    ZEPS(:,:,L) = 1.
  ENDWHERE
  !
  IF ( TRIM(CVAR(L))=='WG2' .OR. TRIM(CVAR(L))=='WG1' ) THEN
    !
    DO I = 1,KI
      ZCOEF(I,:,L) = ZCOFSWI(I)*ZCOFSWI(I)
    ENDDO
    !
  ELSEIF ( TRIM(CVAR(L))=='LAI' .AND. LBFIXED ) THEN
    !
    DO I = 1,KI
      DO J = 1,NPATCH
        IF ( XLAI_PASS(I,J)/=XUNDEF ) THEN
          ZCOEF(I,J,L) = XLAI_PASS(I,J)*XLAI_PASS(I,J)
        ELSE 
          ZCOEF(I,J,L) = 0.
        ENDIF
      ENDDO
    ENDDO
    !
  ELSE
    !
    ZCOEF(:,:,L) = 1.
    !
  ENDIF
  !
ENDDO
!
!############################# B CALCULATION ###############################
!
! ----------------------
! VARASSIM OPTION : LBEV
! ----------------------
!   Calculate the LTM, and evolve B. 
!
! Set the B input file depending of an existing B was found or not
YBGFILE = "BGROUNDin."//YMYPROC
INQUIRE (FILE=TRIM(YBGFILE),EXIST=GBEXISTS)
!
IF ( LBEV .AND. GBEXISTS ) THEN
  !
  ZB(:,:,:) = 0.
  CALL B_BIG_LOOP("READ",YBGFILE,ZB)
  IF ( NPRINTLEV > 0 ) WRITE(*,*) 'read previous B matrix  ==>',ZB(1,1,1),NVAR
  !
ELSEIF ( LBEV .OR. LBFIXED ) THEN
  !
  ! Initialization of B 
  ZB(:,:,:) = 0.
  DO I = 1,KI
    DO L = 1,NVAR
      DO J = 1,NPATCH   
        !
        L1 = J + NPATCH *(L-1)
        ZB(I,L1,L1) = XSIGMA(L)*XSIGMA(L) * ZCOEF(I,J,L)
        !
      ENDDO
    ENDDO
    !
  ENDDO
  !
  IF ( LBEV ) THEN
    !
    ZBLONG(:) = 0.
    ICPT = 0
    !
    DO I = 1,KI
      DO J = 1,NPATCH
        DO JJ = 1,NPATCH
          DO L = 1,NVAR
            DO K = 1,NVAR
              !
              L1 = J + NPATCH * (L-1)
              !
              ICPT = ICPT + 1
              ZBLONG(ICPT) = ZB(I,L1,L1)
              !
            ENDDO
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !
    ZB(:,:,:) = 0.
    CALL B_BIG_LOOP("BUIL","",ZB,ZBLONG)
    IF ( NPRINTLEV > 0 ) WRITE(*,*) 'Initialized B'
    !
  ENDIF
  !
ELSE
  !
  CALL ABOR1_SFX("LBEV or LBFIXED should be .TRUE.!")
  !
ENDIF
!
IF ( LBEV ) THEN
  !
!//////////////////////TO WRITE LTM/////////////////////////////////////
  IUNIT = 110
  DO L=1,NVAR
    DO K=1,NVAR
      IUNIT = IUNIT + 1
      WRITE(YCHAR,'(I1)') K
      YLFNAME='LTM_del'//TRIM(CVAR(K))//'_del'//TRIM(CVAR(L))//"."//YMYPROC
      OPEN(UNIT=IUNIT,FILE=YLFNAME,FORM='FORMATTED',STATUS='UNKNOWN',POSITION='APPEND')
    ENDDO
  ENDDO
!/////////////////////TO WRITE LTM////////////////////////////////////////
  DO I = 1,KI
    !
    ! calculate LTM
    ZLTM(:,:) = 0.0
    IUNIT = 110
    DO L = 1,NVAR    ! control variable (x at previous time step)
      DO K = 1,NVAR 
        IUNIT = IUNIT + 1
        DO J = 1,NPATCH 
          !
          L1 = J + NPATCH*(L-1)
          K1 = J + NPATCH*(K-1)
          !
          IF ( XPATCH(I,J)>0.0 .AND. XF(I,J,L+1,K).NE.XUNDEF .AND. XF(I,J,1,K).NE.XUNDEF ) THEN
            !
            ! Jacobian of fwd model
            ZLTM(L1,K1) = ( XF(I,J,L+1,K) - XF(I,J,1,K) ) / ZEPS(I,J,L)
            ! impose upper/lower limits 
            ZLTM(L1,K1) = MAX(-0.1, ZLTM(L1,K1))
            ZLTM(L1,K1) = MIN( 1.0, ZLTM(L1,K1))
            !
          ENDIF
          !
          WRITE (IUNIT,*) ZLTM(L1,K1)
          !
        ENDDO
      ENDDO
    ENDDO
    !
!//////////////////////TO WRITE LTM/////////////////////////////////////      
    IUNIT = 110
    DO L=1,NVAR
      DO K=1,NVAR
        IUNIT = IUNIT + 1
        CLOSE(IUNIT)
      ENDDO
    ENDDO
!//////////////////////TO WRITE LTM/////////////////////////////////////      
    IF ( NPRINTLEV > 0 ) WRITE(*,*) 'LTM d(wg2)/d(wg2)', ZLTM(1,1)
    !
    ! evolve B 
    ZB(I,:,:) = MATMUL(ZLTM(:,:),MATMUL(ZB(I,:,:),TRANSPOSE(ZLTM(:,:))))
    !
    !
    !   Adding model error to background error matrix 
    ZQ(:,:) = 0.0
    DO L=1,NVAR
      DO J=1,NPATCH
        !
        L1 = J+NPATCH*(L-1)
        !
        ZQ(L1,L1) = XSIGMA(L)*XSIGMA(L)
        !
        IF (TRIM(CVAR(L)) == 'LAI') THEN
          ZQ(L1,L1) = XSCALE_QLAI*XSCALE_QLAI * ZQ(L1,L1)
        ELSE
          ZQ(L1,L1) = XSCALE_Q*XSCALE_Q * ZQ(L1,L1) * ZCOEF(I,J,L)
        ENDIF
        !
      ENDDO
    ENDDO
    !
    ! B is the forecast matrix - need to add Q
    IF ( NPRINTLEV > 0 ) THEN
      WRITE(*,*) 'B before wg2 wg2 ==> ',ZB(I,1,1)/ZCOFSWI(I),ZB(I,1,1)
      WRITE(*,*) 'Q value wg2 wg2 ==> ',ZQ(1,1)/ZCOFSWI(I),ZQ(1,1)
    ENDIF
    !
    ZB(I,:,:) = ZB(I,:,:) + ZQ(:,:)
    !
    IF ( NPRINTLEV > 0 ) WRITE(*,*) 'B after wg2 wg2 ==>',ZB(I,1,1)/ZCOFSWI(1),ZB(I,1,1)
    !
  ENDDO
  !
  ! write out the LTM for the forward model
  ! Write out current B
  YBGFILE="BGROUNDout_LBEV."//YMYPROC
  CALL B_BIG_LOOP("WRIT",YBGFILE,ZB)
  IF ( NPRINTLEV > 0 ) THEN
    WRITE(*,*) 'store B matrix after TL evolution ==>',ZB(1,1,1)
    WRITE(*,*) 'writing out B'
  ENDIF
  !
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
INOBS = 0
IHOUR = 0
ZTIME = FLOAT(NECHGU) * 3600.
!
!############################# READS OBSERVATIONS ###############################
!
! BEGINNING OF TIME LOOP
DO ISTEP=1,NBOUTPUT
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
ENDDO 
!
!//////////////////////TO WRITE OBS/////////////////////////////////////
OPEN (UNIT=111,FILE='OBSout.'//YMYPROC,STATUS='unknown',IOSTAT=ISTAT)
DO I = 1,KI
  IF ( MINVAL(XWGI(I,1,:))>0. ) THEN
    ZYO (I,:) = XUNDEF
    IF ( NPRINTLEV > 0 ) WRITE(*,*) 'OBSERVATION FOR POINT ',I,' REMOVED'
  ENDIF
  WRITE (111,*) ZYO(I,:)
ENDDO
CLOSE(111)
!//////////////////////TO WRITE OBS/////////////////////////////////////
!
!############################# ANALYSIS ###############################
!
IF ( NPRINTLEV > 0 ) THEN
  WRITE(*,*) 'calculating jacobians',INOBS
  WRITE(*,*) ' and then PERFORMING ANALYSIS'
ENDIF
!
!//////////////////////TO WRITE ANALYSIS ARRAYS/////////////////////////////////////
! WRITE OUT OBS AND YERROR FOR DIAGNOSTIC PURPOSES
OPEN (UNIT=111,FILE='OBSERRORout.'//YMYPROC,STATUS='unknown',IOSTAT=ISTAT)
! *** Write innovations in ASCII file ***
OPEN (unit=112,file='INNOV.'//YMYPROC,status='unknown',IOSTAT=ISTAT)
! Write analysis results and increments in ASCII file
OPEN (unit=113,file='ANAL_INCR.'//YMYPROC,status='unknown',IOSTAT=ISTAT)
! **** Write out the observation operator + Gain matrix ****
IUNIT = 113
DO L = 1,NVAR
  DO K=1,NOBSTYPE
    IUNIT = IUNIT + 1
    WRITE(YCHAR,'(I1)') K
    YFNAME='HO_'//CVAR(L)//'_v'//YCHAR
    OPEN(UNIT=IUNIT,FILE=YFNAME,FORM='FORMATTED',STATUS='UNKNOWN',IOSTAT=ISTAT)
  ENDDO
ENDDO
!//////////////////////TO WRITE ANALYSIS ARRAYS/////////////////////////////////////
!
IF (NPATCH==12) THEN
  ZVLAIMIN = (/0.3,0.3,0.3,0.3,1.0,1.0,0.3,0.3,0.3,0.3,0.3,0.3/)
ELSE
  ZVLAIMIN = (/0.3/)
ENDIF
!
IOBSCOUNT = 0
DO I=1,KI
  !
!---------------- MEAN SIMULATED OBS AVERAGED OVER TILES-----------------------
  ZYF(:,:) = 0. 
  DO J=1,NPATCH
    IF (XPATCH(I,J) > 0.0) THEN
      WHERE ( XF_PATCH(I,J,:,:)/=XUNDEF ) 
        ZYF(:,:) = ZYF(:,:) + XPATCH(I,J)*XF_PATCH(I,J,:,:)
      ENDWHERE
    ENDIF
  ENDDO
  IF ( NPRINTLEV > 0 ) WRITE(*,*) 'read in sim obs yf', ZYF(:,1)
  !
  !
  ZR   (:,:) = 0. ! Observation error matrix
  !
  ZHO  (:,:) = XUNDEF  ! Linearized observation matrix
  ZB2  (:)   = XUNDEF  ! Innovation vector
  
  DO ISTEP=1,NBOUTPUT
    !
    DO K = 1,NOBSTYPE
      !
      K1 = (ISTEP-1)*NOBSTYPE + K
      !
!--------------------- SET OBSERVATION ERROR ------------------      
      ZR(K1,K1) = XERROBS(K)*XERROBS(K)
      IF ( COBS(K) .EQ. "LAI" ) THEN
        ZR(K1,K1) = ZR(K1,K1) * ZYO(I,K1)*ZYO(I,K1)
      ELSEIF (COBS(K) .EQ. "WG1") THEN
        ! convert R for wg1 from SWI  to abs value
        ZR(K1,K1) = ZR(K1,K1) * ZCOFSWI(I)*ZCOFSWI(I)
      ENDIF
      !
!--------------------- CALCULATE JACOBIANS ------------------         
      DO L=1,NVAR
        DO J=1,NPATCH
          !
          L1 = J + NPATCH*(L-1)
          !
          IF( ZYO(I,K1).NE.XUNDEF .AND. ZYO(I,K1).NE.999.0 ) THEN         !if obs available
            ! Jacobian of obs operator
            ZHO(K1,L1) = XPATCH(I,J)*(XF_PATCH(I,J,L+1,K) - XF_PATCH(I,J,1,K))/ZEPS(I,J,L)
            ! impose limits  
            ZHO(K1,L1) = MAX(-0.1, ZHO(K1,L1))
            ZHO(K1,L1) = MIN( 1.0, ZHO(K1,L1))
            ! innovation vector
            ZB2(K1) = ZYO(I,K1) - ZYF(1,K)
            IOBSCOUNT = IOBSCOUNT + 1
          ELSE  !if no obs available
            ! set obs operator and innovation to zero if no obs available
            ZHO(K1,L1) = 0.0
            ZB2(K1) = 0.0 
          ENDIF
          !
        ENDDO
      ENDDO
      !
    ENDDO
    !
  ENDDO
  
  WRITE(111,*) ZR(:,:)
  WRITE(112,*) ZB2(:)
  
!---------------******  SOIL ANALYSIS *******--------------------------
  ZHOT(:,:) = 0.
  ZK1(:,:) = 0.
  ZP(:) = 0.
  ZX(:) = 0.
  ZXINCR(:) = 0.
  !
  ZHOT(:,:) = TRANSPOSE(ZHO(:,:))
  ZK1 (:,:) = MATMUL(ZHO(:,:),MATMUL(ZB(I,:,:),ZHOT(:,:))) + ZR(:,:)
  CALL CHOLDC(NOBSTYPE,ZK1(:,:),ZP(:))                         ! Cholesky decomposition (1)
  CALL CHOLSL(NOBSTYPE,ZK1(:,:),ZP(:),ZB2(:),ZX(:))            ! Cholesky decomposition (2)
  ZXINCR(:) = MATMUL(ZB(I,:,:),MATMUL(ZHOT(:,:),ZX(:)))
  
  DO L=1,NVAR
    DO J=1,NPATCH
      !
      L1 = J+NPATCH*(L-1)
      !
      ! Update the modified values
      IF ( TRIM(CVAR(L))=="LAI" ) THEN
        ZXINCR(L1) = MAX( ZXINCR(L1), ZVLAIMIN(J)-XF(I,J,1,L) )
        XBIO_PASS(I,J) = XBIO_PASS(I,J) + ZXINCR(L1)*XALPH(J)
      ELSEIF ( XF(I,J,1,L)+ZXINCR(L1)<0. ) THEN
        ZXINCR(L1) = 0.
      ENDIF
      !
      XF(I,J,1,L) = XF(I,J,1,L) + ZXINCR(L1)
      !
      ! For no only warn if we have negative values.
      IF ( NPRINTLEV > 0 ) THEN
        IF ( XF(I,J,1,L) < 0. ) WRITE(*,*) "WARNING X<0. for ",I,J," for variable ",TRIM(CVAR(L))
      ENDIF
      !
    ENDDO
  ENDDO
  
  DO J=1,NPATCH
    WRITE(113,*) (XF(I,J,1,L),L=1,NVAR), (ZXINCR(J+NPATCH*(L-1)),L=1,NVAR)
  ENDDO
  
  
!--------------------ANALYSIS OF B (FOR USE IN NEXT CYCLE)-------------------
  ! Ba = (I-KH)Bf(I-KH)t+KRKt
  ! K = BfHt{K1}**-1
  ! K1 needs PATCH dim added
  
  ZGAIN(:,:) = 0.
  ZIDKH(:,:) = 0.
  ZKRK(:,:) = 0.
  !
  ! K1 = (R+H.B.HT) (calculate inverse -> output goes to K1)
  CALL INVERSE_MATRIX(INOBS,ZK1(:,:),ZP(:))
  ZGAIN(:,:) = MATMUL(ZB(I,:,:),MATMUL(ZHOT(:,:),ZK1(:,:)))
  ZIDKH(:,:) = ZIDENT(:,:) - MATMUL(ZGAIN(:,:),ZHO(:,:))
  ZKRK (:,:) = MATMUL(ZGAIN(:,:),MATMUL(ZR(:,:),TRANSPOSE(ZGAIN(:,:))))
  IF (.NOT.LBFIXED)  ZB(I,:,:) = MATMUL(ZIDKH(:,:),MATMUL(ZB(I,:,:),TRANSPOSE(ZIDKH(:,:)))) + ZKRK(:,:)
  
  IUNIT = 113
  DO L = 1,NVAR
    DO K = 1,NOBSTYPE
      IUNIT = IUNIT + 1
      DO J=1,NPATCH
        WRITE(IUNIT,*) ZHO(K,J+NPATCH*(L-1)),ZGAIN(J+NPATCH*(L-1),K)
      ENDDO
    ENDDO
  ENDDO
  !  
ENDDO
!
!
!//////////////////////TO WRITE ANALYSIS ARRAYS/////////////////////////////////////
CLOSE(111)
CLOSE(112)
CLOSE(113)
IUNIT = 113
DO L = 1,NVAR
  DO K = 1,NOBSTYPE
    IUNIT = IUNIT + 1
    CLOSE(IUNIT)
  ENDDO
ENDDO
!//////////////////////TO WRITE ANALYSIS ARRAYS/////////////////////////////////////
!
! Write out analysed B (for use in next cycle)
YBGFILE = "BGROUNDout_ASSIM."//YMYPROC
CALL B_BIG_LOOP("WRIT",YBGFILE,ZB)
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
!############################# GET VARIABLES FOR OUTPUT WRITING ###############################
DO L=1,NVAR
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
          XBIOMASS(:,1,:) = XBIO_PASS(:,:)
        CASE("BIOMA2","BIOMASS2")
          XBIOMASS(:,2,:) = XBIO_PASS(:,:)
        CASE("RESPI1","RESP_BIOM1")
          XRESP_BIOMASS(:,1,:) = XBIO_PASS(:,:)
        CASE("RESPI2","RESP_BIOM2")
          XRESP_BIOMASS(:,2,:) = XBIO_PASS(:,:)
        CASE("LAI")
          XLAI(:,:) = XBIO_PASS(:,:)
        CASE DEFAULT
          CALL ABOR1_SFX("Mapping of "//CBIO//" is not defined in EKF!")
      END SELECT
    CASE DEFAULT
      CALL ABOR1_SFX("Mapping of "//TRIM(CVAR(L))//" is not defined in EKF!")
  END SELECT
ENDDO
!
!
IF (LHOOK) CALL DR_HOOK('ASSIM_NATURE_ISBA_EKF',1,ZHOOK_HANDLE)
!
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  Subroutines of Linear Algebra  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE B_BIG_LOOP(HACTION,HFILE,PTAB,PTAB_IN)
!
CHARACTER(LEN=4), INTENT(IN) :: HACTION
CHARACTER(LEN=*), INTENT(IN) :: HFILE
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PTAB
REAL, DIMENSION(:), INTENT(IN), OPTIONAL :: PTAB_IN
INTEGER :: I,J,K,L,JJ,L1,K1
!
IF (HACTION=="READ" .OR. HACTION=="WRIT") THEN
  OPEN (UNIT=111,FILE=HFILE,STATUS='unknown',IOSTAT=ISTAT)
ELSE
  ICPT = 0
ENDIF
!
DO L = 1,NVAR   ! control variable (x at previous time step)
  DO K = 1,NVAR
    DO I = 1,KI
      DO J = 1,NPATCH   
        DO JJ = 1,NPATCH
          !
          L1 = J+NPATCH*(L-1)
          K1 = JJ+NPATCH*(K-1)
          !
          IF ( HACTION=="READ" ) THEN
            READ (111,*) PTAB(I,L1,K1)
          ELSEIF ( HACTION=="WRIT" ) THEN
            WRITE(111,*) PTAB(I,L1,K1)
          ELSE
            ICPT = ICPT+1
            PTAB(I,L1,K1) = PTAB_IN(ICPT)
          ENDIF
          ! 
        ENDDO
      ENDDO
    ENDDO
  ENDDO
ENDDO
!
IF (HACTION=="READ" .OR. HACTION=="WRIT") CLOSE(111)
!
END SUBROUTINE B_BIG_LOOP

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

