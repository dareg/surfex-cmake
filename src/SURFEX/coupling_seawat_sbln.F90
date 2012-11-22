!     ###############################################################################
SUBROUTINE COUPLING_SEAWAT_SBL_n(HPROGRAM, HCOUPLING, HSURF,                               &
               PTSTEP, KYEAR, KMONTH, KDAY, PTIME, KI, KSV, KSW, PTSUN, PZENITH, PZENITH2, &
               PAZIM, PZREF, PUREF, PZS, PU, PV, PQA, PTA, PRHOA, PSV, PCO2, HSV,          &
               PRAIN, PSNOW, PLW, PDIR_SW, PSCA_SW, PSW_BANDS, PPS, PPA,                   &
               PSFTQ, PSFTH, PSFTS, PSFCO2, PSFU, PSFV, OSBL, PSST, PZ0,                   &
               PZ, PXU, KLVL, PTKE, PT, PQ, PLMO, PZF, PDZ, PDZF, PP,                      &
               K2M, PT2M, PQ2M, PHU2M, PZON10M, PMER10M, PWIND10M, PWIND10M_MAX,           &
               PT2M_MIN, PT2M_MAX, PHU2M_MIN, PHU2M_MAX,                                   &
               PTRAD, PDIR_ALB, PSCA_ALB, PEMIS,                                           &
               PPEW_A_COEF, PPEW_B_COEF,                                                   &
               PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,                         &
               HTEST                                                                       )
!     ###############################################################################
!
!!****  *COUPLING_SEAWAT_SBL_n * - Adds a SBL into SEAFLUX
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     V. Masson 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    09/2007
!!      V. Masson   05/2009 Implicitation of momentum fluxes
!!      S. Riette   06/2009 Initialisation of XT, PQ, XU and XTKE on canopy levels
!!      S. Riette   10/2009 Iterative computation of XZ0
!!      S. Riette   01/2010 Use of interpol_sbl to compute 10m wind diagnostic
!----------------------------------------------------------------
!
USE MODD_SURF_PAR,         ONLY : XUNDEF
USE MODD_CSTS,             ONLY : XCPD
! 
USE MODE_COUPLING_CANOPY
!
USE MODI_INIT_WATER_SBL
!
USE MODI_CANOPY_EVOL
USE MODI_CANOPY_GRID_UPDATE
!
USE MODI_COUPLING_SEAFLUX_n
USE MODI_COUPLING_WATFLUX_n
USE MODI_COUPLING_FLAKE_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=6),    INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
CHARACTER(LEN=1),    INTENT(IN)  :: HCOUPLING ! type of coupling
                                              ! 'E' : explicit
                                              ! 'I' : implicit
CHARACTER(LEN=1),    INTENT(IN)  :: HSURF     ! 
INTEGER,             INTENT(IN)  :: KYEAR     ! current year (UTC)
INTEGER,             INTENT(IN)  :: KMONTH    ! current month (UTC)
INTEGER,             INTENT(IN)  :: KDAY      ! current day (UTC)
REAL,                INTENT(IN)  :: PTIME     ! current time since midnight (UTC, s)
INTEGER,             INTENT(IN)  :: KI        ! number of points
INTEGER,             INTENT(IN)  :: KSV       ! number of scalars
INTEGER,             INTENT(IN)  :: KSW       ! number of short-wave spectral bands
REAL, DIMENSION(:), INTENT(IN)  :: PTSUN     ! solar time                    (s from midnight)
REAL,                INTENT(IN)  :: PTSTEP    ! atmospheric time-step                 (s)
REAL, DIMENSION(:), INTENT(IN)  :: PZREF     ! height of T,q forcing                 (m)
REAL, DIMENSION(:), INTENT(IN)  :: PUREF     ! height of wind forcing                (m)
!
REAL, DIMENSION(:), INTENT(IN)  :: PTA       ! air temperature forcing               (K)
REAL, DIMENSION(:), INTENT(IN)  :: PQA       ! air humidity forcing                  (kg/m3)
REAL, DIMENSION(:), INTENT(IN)  :: PRHOA     ! air density                           (kg/m3)
REAL, DIMENSION(:,:),INTENT(IN) :: PSV     ! scalar variables
!                                             ! chemistry:   first char. in HSV: '#'  (molecule/m3)
!                                             !
CHARACTER(LEN=6), DIMENSION(:),INTENT(IN):: HSV  ! name of all scalar variables
REAL, DIMENSION(:), INTENT(IN)  :: PU        ! zonal wind                            (m/s)
REAL, DIMENSION(:), INTENT(IN)  :: PV        ! meridian wind                         (m/s)
REAL, DIMENSION(:,:),INTENT(IN) :: PDIR_SW ! direct  solar radiation (on horizontal surf.)
!                                             !                                       (W/m2)
REAL, DIMENSION(:,:),INTENT(IN) :: PSCA_SW ! diffuse solar radiation (on horizontal surf.)
!                                             !                                       (W/m2)
REAL, DIMENSION(:),INTENT(IN)  :: PSW_BANDS ! mean wavelength of each shortwave band (m)
REAL, DIMENSION(:), INTENT(IN)  :: PZENITH   ! zenithal angle at t  (radian from the vertical)
REAL, DIMENSION(:), INTENT(IN)  :: PZENITH2  ! zenithal angle at t+1(radian from the vertical)
REAL, DIMENSION(:), INTENT(IN)  :: PAZIM     ! azimuthal angle      (radian from North, clockwise)
REAL, DIMENSION(:), INTENT(IN)  :: PLW       ! longwave radiation (on horizontal surf.)
!                                             !                                       (W/m2)
REAL, DIMENSION(:), INTENT(IN)  :: PPS       ! pressure at atmospheric model surface (Pa)
REAL, DIMENSION(:), INTENT(IN)  :: PPA       ! pressure at forcing level             (Pa)
REAL, DIMENSION(:), INTENT(IN)  :: PZS       ! atmospheric model CANOPY           (m)
REAL, DIMENSION(:), INTENT(IN)  :: PCO2      ! CO2 concentration in the air          (kg/m3)
REAL, DIMENSION(:), INTENT(IN)  :: PSNOW     ! snow precipitation                    (kg/m2/s)
REAL, DIMENSION(:), INTENT(IN)  :: PRAIN     ! liquid precipitation                  (kg/m2/s)
!
!
REAL, DIMENSION(:), INTENT(OUT) :: PSFTH     ! flux of heat                          (W/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PSFTQ     ! flux of water vapor                   (kg/m2/s)
REAL, DIMENSION(:), INTENT(OUT) :: PSFU      ! zonal momentum flux                   (Pa)
REAL, DIMENSION(:), INTENT(OUT) :: PSFV      ! meridian momentum flux                (Pa)
REAL, DIMENSION(:), INTENT(OUT) :: PSFCO2    ! flux of CO2                           (kg/m2/s)
REAL, DIMENSION(:,:),INTENT(OUT):: PSFTS   ! flux of scalar var.                   (kg/m2/s)
!
LOGICAL, INTENT(IN) :: OSBL
REAL, DIMENSION(:), INTENT(IN) :: PSST
REAL, DIMENSION(:), INTENT(INOUT) :: PZ0
REAL, DIMENSION(:,:), INTENT(INOUT) :: PZ
REAL, DIMENSION(:,:), INTENT(INOUT) :: PXU
INTEGER, INTENT(IN) :: KLVL
REAL, DIMENSION(:,:), INTENT(INOUT) :: PTKE
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT
REAL, DIMENSION(:,:), INTENT(INOUT) :: PQ
REAL, DIMENSION(:), INTENT(INOUT) :: PLMO
REAL, DIMENSION(:,:), INTENT(INOUT) :: PZF
REAL, DIMENSION(:,:), INTENT(INOUT) :: PDZ
REAL, DIMENSION(:,:), INTENT(INOUT) :: PDZF
REAL, DIMENSION(:,:), INTENT(INOUT) :: PP
INTEGER, INTENT(IN) :: K2M
REAL, DIMENSION(:), INTENT(INOUT) :: PT2M
REAL, DIMENSION(:), INTENT(INOUT) :: PQ2M
REAL, DIMENSION(:), INTENT(INOUT) :: PHU2M
REAL, DIMENSION(:), INTENT(INOUT) :: PZON10M
REAL, DIMENSION(:), INTENT(INOUT) :: PMER10M
REAL, DIMENSION(:), INTENT(INOUT) :: PWIND10M
REAL, DIMENSION(:), INTENT(INOUT) :: PWIND10M_MAX
REAL, DIMENSION(:), INTENT(INOUT) :: PT2M_MIN
REAL, DIMENSION(:), INTENT(INOUT) :: PT2M_MAX
REAL, DIMENSION(:), INTENT(INOUT) :: PHU2M_MIN
REAL, DIMENSION(:), INTENT(INOUT) :: PHU2M_MAX
!
REAL, DIMENSION(:), INTENT(OUT) :: PTRAD     ! radiative temperature                 (K)
REAL, DIMENSION(:,:),INTENT(OUT):: PDIR_ALB! direct albedo for each spectral band  (-)
REAL, DIMENSION(:,:),INTENT(OUT):: PSCA_ALB! diffuse albedo for each spectral band (-)
REAL, DIMENSION(:), INTENT(OUT) :: PEMIS     ! emissivity                            (-)
!
REAL, DIMENSION(:), INTENT(IN) :: PPEW_A_COEF! implicit coefficients
REAL, DIMENSION(:), INTENT(IN) :: PPEW_B_COEF! needed if HCOUPLING='I'
REAL, DIMENSION(:), INTENT(IN) :: PPET_A_COEF
REAL, DIMENSION(:), INTENT(IN) :: PPEQ_A_COEF
REAL, DIMENSION(:), INTENT(IN) :: PPET_B_COEF
REAL, DIMENSION(:), INTENT(IN) :: PPEQ_B_COEF
CHARACTER(LEN=2),    INTENT(IN) :: HTEST ! must be equal to 'OK'
!
!*      0.2    declarations of local variables
!
!* forcing variables
!
REAL, DIMENSION(SIZE(PTA))     :: ZWIND    ! lowest atmospheric level wind speed           (m/s)
REAL, DIMENSION(SIZE(PTA))     :: ZEXNA    ! Exner function at lowest SBL scheme level     (-)
REAL, DIMENSION(SIZE(PTA))     :: ZTA      ! temperature                                   (K)
REAL, DIMENSION(SIZE(PTA))     :: ZPA      ! pressure                                      (Pa)
REAL, DIMENSION(SIZE(PTA))     :: ZZREF    ! temperature forcing level                     (m)
REAL, DIMENSION(SIZE(PTA))     :: ZUREF    ! wind        forcing level                     (m)
REAL, DIMENSION(SIZE(PTA))     :: ZU       ! zonal wind                                    (m/s)
REAL, DIMENSION(SIZE(PTA))     :: ZV       ! meridian wind                                 (m/s)
REAL, DIMENSION(SIZE(PTA))     :: ZQA      ! specific humidity                             (kg/m3)
REAL, DIMENSION(SIZE(PTA))     :: ZPEQ_A_COEF ! specific humidity implicit
REAL, DIMENSION(SIZE(PTA))     :: ZPEQ_B_COEF ! coefficients (hum. in kg/kg)
!
!
! SBL turbulence scheme
!
REAL, DIMENSION(SIZE(PTA))        :: ZSFLUX_U  ! Surface flux u'w' (m2/s2)
REAL, DIMENSION(SIZE(PTA))        :: ZSFLUX_T  ! Surface flux w'T' (mK/s)
REAL, DIMENSION(SIZE(PTA))        :: ZSFLUX_Q  ! Surface flux w'q' (kgm2/s)
REAL, DIMENSION(SIZE(PTA),KLVL)   :: ZFORC_U   ! tendency due to drag force for wind
REAL, DIMENSION(SIZE(PTA),KLVL)   :: ZDFORC_UDU! formal derivative of
!                                              ! tendency due to drag force for wind
REAL, DIMENSION(SIZE(PTA),KLVL)   :: ZFORC_E   ! tendency due to drag force for TKE
REAL, DIMENSION(SIZE(PTA),KLVL)   :: ZDFORC_EDE! formal derivative of
!                                              ! tendency due to drag force for TKE
REAL, DIMENSION(SIZE(PTA),KLVL)   :: ZFORC_T   ! tendency due to drag force for Temp
REAL, DIMENSION(SIZE(PTA),KLVL)   :: ZDFORC_TDT! formal derivative of
!                                              ! tendency due to drag force for Temp
REAL, DIMENSION(SIZE(PTA),KLVL)   :: ZFORC_Q   ! tendency due to drag force for Temp
REAL, DIMENSION(SIZE(PTA),KLVL)   :: ZDFORC_QDQ! formal derivative of
!                                              ! tendency due to drag force for hum.
REAL, DIMENSION(SIZE(PTA),KLVL)   :: ZLMO      ! MO length
REAL, DIMENSION(SIZE(PTA),KLVL)   :: ZLM       ! mixing length
REAL, DIMENSION(SIZE(PTA),KLVL)   :: ZLEPS     ! dissipative length
REAL, DIMENSION(SIZE(PTA))     :: ZH           ! canopy height (m)
REAL, DIMENSION(SIZE(PTA))     :: ZUSTAR       ! friction velocity (m/s)
!
REAL, DIMENSION(SIZE(PTA))     :: ZPET_A_COEF ! temperature implicit
REAL, DIMENSION(SIZE(PTA))     :: ZPET_B_COEF ! coefficients (K)
REAL, DIMENSION(SIZE(PTA))     :: ZPEW_A_COEF ! wind implicit
REAL, DIMENSION(SIZE(PTA))     :: ZPEW_B_COEF ! coefficients (m/s)

REAL, DIMENSION(SIZE(PTA))   :: ZALFAU   ! V+(1) = - alfa rho u'w'(1) + beta
REAL, DIMENSION(SIZE(PTA))   :: ZBETAU   ! V+(1) = - alfa rho u'w'(1) + beta
REAL, DIMENSION(SIZE(PTA))   :: ZALFATH  ! Th+(1) = - alfa rho w'th'(1) + beta
REAL, DIMENSION(SIZE(PTA))   :: ZBETATH  ! Th+(1) = - alfa rho w'th'(1) + beta
REAL, DIMENSION(SIZE(PTA))   :: ZALFAQ   ! Q+(1) = - alfa rho w'q'(1) + beta
REAL, DIMENSION(SIZE(PTA))   :: ZBETAQ   ! Q+(1) = - alfa rho w'q'(1) + beta
!
CHARACTER(LEN=1) :: GCOUPLING
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
!
!*      1.     Preliminary computations of the SBL scheme
!              ------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('COUPLING_SEAWAT_SBL_N',0,ZHOOK_HANDLE)
IF (OSBL) THEN
!
!*      1.1    Updates SBL vertical grid as a function of forcing height
!              ---------------------------------------------------------
!
!* determines where is the forcing level and modifies the upper levels of the canopy grid
!
  ZH = 0.
  CALL CANOPY_GRID_UPDATE(KI,KLVL,ZH,PUREF,PZ,PZF,PDZ,PDZF)
!
!
!
!*     1.2     Initialisation at first time step
!              ---------------------------------
!
  IF(ANY(PT(:,:) == XUNDEF)) THEN
    CALL INIT_WATER_SBL(KLVL, PPA, PPS, PTA, PQA, PRHOA, PU, PV, PRAIN, PSNOW,     &
                        PSFTH, PSFTQ, PZREF, PUREF, PSST, PZ0, PZ,                 &
                        PT, PQ, PXU, PTKE, PP)
  ENDIF
!
!
!*      1.3    Allocations
!              -----------
!
  CALL INIT_FORC( ZFORC_U, ZDFORC_UDU, ZFORC_E, ZDFORC_EDE, &
                 ZFORC_T, ZDFORC_TDT, ZFORC_Q, ZDFORC_QDQ )
!
  ZSFLUX_U = 0.
  ZSFLUX_T = 0.
  ZSFLUX_Q = 0.
!
  ZLMO = SPREAD(PLMO,2,KLVL)
!
!*      1.3   Computes coefficients for implicitation
!             ---------------------------------------
!
  ZWIND = SQRT(PU**2+PV**2)
  CALL CANOPY_EVOL(KI,KLVL,PTSTEP,1,PZ,ZWIND,PTA,PQA,PPA,PRHOA,           &
                 ZSFLUX_U,ZSFLUX_T,ZSFLUX_Q,                              &
                 ZFORC_U,ZDFORC_UDU,ZFORC_E,ZDFORC_EDE,                   &
                 ZFORC_T,ZDFORC_TDT,ZFORC_Q,ZDFORC_QDQ,                   &
                 PZ,PZF,PDZ,PDZF,PXU,PTKE,PT,PQ,ZLMO,ZLM,ZLEPS,PP,ZUSTAR, &
                 ZALFAU,ZBETAU,ZALFATH,ZBETATH,ZALFAQ,ZBETAQ              )

!
!*     1.5     Goes from atmospheric forcing to canopy forcing height
!              ------------------------------------------------------
!
  IF (HSURF=='F') THEN
    GCOUPLING = 'E'
  ELSEIF (HSURF=='S' .OR. HSURF=='W') THEN
    GCOUPLING ='I'
  ENDIF
!
  CALL INIT_COUPLING_CANOPY( PP(:,1), PPA, PT(:,1), PQ(:,1), &
                           PU, PV, PZ(:,1), PXU(:,1),        &
                           PRHOA, ZALFAU, ZBETAU, ZALFATH,   &
                           ZBETATH, ZALFAQ, ZBETAQ,          &
                           ZPA, ZTA, ZQA, ZU, ZV,            &
                           ZUREF, ZZREF, ZEXNA,              &
                           ZPEW_A_COEF, ZPEW_B_COEF,         &
                           ZPET_A_COEF, ZPET_B_COEF,         &
                           ZPEQ_A_COEF, ZPEQ_B_COEF          )
!
!-------------------------------------------------------------------------------------
ELSE
!-------------------------------------------------------------------------------------
!
!*      2.     If no SBL scheme is used, forcing is not modified
!              -------------------------------------------------
!
  GCOUPLING = HCOUPLING
!
  CALL INIT_COUPLING( HCOUPLING,                  &
                      PPS, PPA, PTA, PQA, PU, PV, &
                      PUREF, PZREF,               &
                      PPEW_A_COEF, PPEW_B_COEF,   &
                      PPET_A_COEF, PPET_B_COEF,   &
                      PPEQ_A_COEF, PPEQ_B_COEF,   &
                      ZPA, ZTA, ZQA, ZU, ZV,      &
                      ZUREF, ZZREF,               &
                      ZPEW_A_COEF, ZPEW_B_COEF,   &
                      ZPET_A_COEF, ZPET_B_COEF,   &
                      ZPEQ_A_COEF, ZPEQ_B_COEF    ) 
!
END IF
!
!-------------------------------------------------------------------------------------
!
!*      2.     Call of SEAFLUX
!              ------------
!
IF (HSURF=='S') THEN
  CALL COUPLING_SEAFLUX_n(HPROGRAM, GCOUPLING,                                             &
             PTSTEP, KYEAR, KMONTH, KDAY, PTIME,                                           &
             KI, KSV, KSW,                                                                 &
             PTSUN, PZENITH, PZENITH2, PAZIM,                                              &
             ZZREF, ZUREF, PZS, ZU, ZV, ZQA, ZTA, PRHOA, PSV, PCO2, HSV,                   &
             PRAIN, PSNOW, PLW, PDIR_SW, PSCA_SW, PSW_BANDS, PPS, ZPA,                     &
             PSFTQ, PSFTH, PSFTS, PSFCO2, PSFU, PSFV,                                      &
             PTRAD, PDIR_ALB, PSCA_ALB, PEMIS,                                             &
             ZPEW_A_COEF, ZPEW_B_COEF,                                                     &
             ZPET_A_COEF, ZPEQ_A_COEF, ZPET_B_COEF, ZPEQ_B_COEF,                           &
             'OK'                                                                          )
ELSEIF (HSURF=='W') THEN
  CALL COUPLING_WATFLUX_n(HPROGRAM, GCOUPLING,                                             &
             PTSTEP, KYEAR, KMONTH, KDAY, PTIME,                                           &
             KI, KSV, KSW,                                                                 &
             PTSUN, PZENITH, PZENITH2, PAZIM,                                              &
             ZZREF, ZUREF, PZS, ZU, ZV, ZQA, ZTA, PRHOA, PSV, PCO2, HSV,                   &
             PRAIN, PSNOW, PLW, PDIR_SW, PSCA_SW, PSW_BANDS, PPS, ZPA,                     &
             PSFTQ, PSFTH, PSFTS, PSFCO2, PSFU, PSFV,                                      &
             PTRAD, PDIR_ALB, PSCA_ALB, PEMIS,                                             &
             ZPEW_A_COEF, ZPEW_B_COEF,                                                     &
             ZPET_A_COEF, ZPEQ_A_COEF, ZPET_B_COEF, ZPEQ_B_COEF,                           &
             'OK'                                                                          )
ELSEIF (HSURF=='F') THEN
  CALL COUPLING_FLAKE_n(HPROGRAM, GCOUPLING,                                                &
              PTSTEP, KYEAR, KMONTH, KDAY, PTIME,                                           &
              KI, KSV, KSW,                                                                 &
              PTSUN, PZENITH, PZENITH2, PAZIM,                                              &
              ZZREF, ZUREF, PZS, ZU, ZV, ZQA, ZTA, PRHOA, PSV, PCO2, HSV,                   &
              PRAIN, PSNOW, PLW, PDIR_SW, PSCA_SW, PSW_BANDS, PPS, ZPA,                     &
              PSFTQ, PSFTH, PSFTS, PSFCO2, PSFU, PSFV,                                      &
              PTRAD, PDIR_ALB, PSCA_ALB, PEMIS,                                             &
              ZPEW_A_COEF, ZPEW_B_COEF,                                                     &
              ZPET_A_COEF, ZPEQ_A_COEF, ZPET_B_COEF, ZPEQ_B_COEF,                           &
              'OK'                                                                          )
ENDIF
!
!-------------------------------------------------------------------------------------
!
!*      3.     End if no SBL is used
!              ---------------------
!
IF (.NOT. OSBL .AND. LHOOK) CALL DR_HOOK('COUPLING_SEAWAT_SBL_N',1,ZHOOK_HANDLE)
IF (.NOT. OSBL) RETURN
!
!-------------------------------------------------------------------------------------
!
!*      4.     Computes the impact of canopy and surfaces on air
!              -------------------------------------------------
!
CALL INIT_FORC( ZFORC_U, ZDFORC_UDU, ZFORC_E, ZDFORC_EDE, &
               ZFORC_T, ZDFORC_TDT, ZFORC_Q, ZDFORC_QDQ )
!
ZSFLUX_U = - SQRT(PSFU(:)**2+PSFV(:)**2) / PRHOA(:)
ZSFLUX_T(:) = PSFTH(:) / XCPD * ZEXNA(:) / PRHOA(:)
ZSFLUX_Q(:) = PSFTQ(:)
!
!-------------------------------------------------------------------------------------
!
!*      6.    Evolution of canopy air due to these impacts
!             --------------------------------------------
!
ZWIND = SQRT(PU**2+PV**2)
CALL CANOPY_EVOL(KI,KLVL,PTSTEP,2,PZ,ZWIND,PTA,PQA,PPA,PRHOA,                 &
                 ZSFLUX_U,ZSFLUX_T,ZSFLUX_Q,                                  &
                 ZFORC_U,ZDFORC_UDU,ZFORC_E,ZDFORC_EDE,                       &
                 ZFORC_T,ZDFORC_TDT,ZFORC_Q,ZDFORC_QDQ,                       &
                 PZ,PZF,PDZ,PDZF,PXU,PTKE,PT,PQ,ZLMO,ZLM,ZLEPS,PP,ZUSTAR,     &
                 ZALFAU,ZBETAU,ZALFATH,ZBETATH,ZALFAQ,ZBETAQ                  )
!
PLMO(:) = ZLMO(:,KLVL)
!
!-------------------------------------------------------------------------------------
!
!*      7.    2m and 10m diagnostics if canopy is used
!             ----------------------------------------
!
!
IF (OSBL .AND. K2M>=1) CALL INIT_2M_10M( PP(:,2), PT(:,2), PQ(:,2), PXU, PZ, &
                                         PU, PV, ZWIND, PRHOA,               &
                                         PT2M, PQ2M, PHU2M, PZON10M, PMER10M,&
                                         PWIND10M, PWIND10M_MAX, PT2M_MIN,   &
                                         PT2M_MAX, PHU2M_MIN, PHU2M_MAX      )
!
IF (LHOOK) CALL DR_HOOK('COUPLING_SEAWAT_SBL_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE COUPLING_SEAWAT_SBL_n
