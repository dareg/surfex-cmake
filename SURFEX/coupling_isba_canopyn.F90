!     ###############################################################################
SUBROUTINE COUPLING_ISBA_CANOPY_n(HPROGRAM, HCOUPLING,                                     &
               PTSTEP, KYEAR, KMONTH, KDAY, PTIME, KI, KSV, KSW, PTSUN, PZENITH, PZENITH2, &
               PAZIM, PZREF, PUREF, PZS, PU, PV, PQA, PTA, PRHOA, PSV, PCO2, HSV,          &
               PRAIN, PSNOW, PLW, PDIR_SW, PSCA_SW, PSW_BANDS, PPS, PPA,                   &
               PSFTQ, PSFTH, PSFTS, PSFCO2, PSFU, PSFV,                                    &
               PTRAD, PDIR_ALB, PSCA_ALB, PEMIS,                                           &
               PPEW_A_COEF, PPEW_B_COEF,                                                   &
               PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,                         &
               HTEST                                                                       )
!     ###############################################################################
!
!!****  *COUPLING_ISBA_CANOPY_n * - Adds a SBL into ISBA
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
!!      S. Riette   06/2009 Initialisation of XT, XQ, XU and XTKE on canopy levels
!!      S. Riette   01/2010 Use of interpol_sbl to compute 10m wind diagnostic
!----------------------------------------------------------------
!
USE MODE_THERMOS
USE MODD_CSTS,          ONLY : XCPD, XRD, XP00, XLVTT, XG
USE MODD_SURF_ATM,      ONLY : XWINDMIN
USE MODD_ISBA_n,        ONLY : LCANOPY, LCANOPY_DRAG, CROUGH, XZ0, XLAI, XPATCH, &
                               XSSO_STDEV, XSSO_SLOPE, XZ0_O_Z0H, XTG, CISBA,    &
                               TSNOW, CCPSURF, XWFC, XVEG, XLAI, XGAMMA, XRSMIN, &
                               XWR, XWRMAX_CF, XRESA, XRGL, XWSAT, XWG, XWGI
USE MODD_ISBA_CANOPY_n, ONLY : XZ, XU, NLVL, XTKE, XT, XQ, XLMO, XZF, XDZ, XDZF, XP
USE MODD_DIAG_ISBA_n,   ONLY : N2M, XAVG_T2M, XAVG_Q2M, XAVG_HU2M,                 &
                               XAVG_ZON10M, XAVG_MER10M, &
                               LSURF_BUDGET, XAVG_FMU, XAVG_FMV
USE MODD_SURF_PAR,      ONLY : XUNDEF
USE MODD_CTURBS,        ONLY : XALPSBL
! 
USE MODI_INIT_ISBA_SBL
USE MODI_ISBA_CANOPY
USE MODI_SSO_BELJAARS04
USE MODI_INTERPOL_SBL
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_CANOPY_EVOL
!
USE MODI_CANOPY_GRID_UPDATE
!
USE MODI_COUPLING_ISBA_n
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=6),    INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
CHARACTER(LEN=1),    INTENT(IN)  :: HCOUPLING ! type of coupling
                                              ! 'E' : explicit
                                              ! 'I' : implicit
INTEGER,             INTENT(IN)  :: KYEAR     ! current year (UTC)
INTEGER,             INTENT(IN)  :: KMONTH    ! current month (UTC)
INTEGER,             INTENT(IN)  :: KDAY      ! current day (UTC)
REAL,                INTENT(IN)  :: PTIME     ! current time since midnight (UTC, s)
INTEGER,             INTENT(IN)  :: KI        ! number of points
INTEGER,             INTENT(IN)  :: KSV       ! number of scalars
INTEGER,             INTENT(IN)  :: KSW       ! number of short-wave spectral bands
REAL, DIMENSION(KI), INTENT(IN)  :: PTSUN     ! solar time                    (s from midnight)
REAL,                INTENT(IN)  :: PTSTEP    ! atmospheric time-step                 (s)
REAL, DIMENSION(KI), INTENT(IN)  :: PZREF     ! height of T,q forcing                 (m)
REAL, DIMENSION(KI), INTENT(IN)  :: PUREF     ! height of wind forcing                (m)
!
REAL, DIMENSION(KI), INTENT(IN)  :: PTA       ! air temperature forcing               (K)
REAL, DIMENSION(KI), INTENT(IN)  :: PQA       ! air humidity forcing                  (kg/m3)
REAL, DIMENSION(KI), INTENT(IN)  :: PRHOA     ! air density                           (kg/m3)
REAL, DIMENSION(KI,KSV),INTENT(IN) :: PSV     ! scalar variables
!                                             ! chemistry:   first char. in HSV: '#'  (molecule/m3)
!                                             !
CHARACTER(LEN=6), DIMENSION(KSV),INTENT(IN):: HSV  ! name of all scalar variables
REAL, DIMENSION(KI), INTENT(IN)  :: PU        ! zonal wind                            (m/s)
REAL, DIMENSION(KI), INTENT(IN)  :: PV        ! meridian wind                         (m/s)
REAL, DIMENSION(KI,KSW),INTENT(IN) :: PDIR_SW ! direct  solar radiation (on horizontal surf.)
!                                             !                                       (W/m2)
REAL, DIMENSION(KI,KSW),INTENT(IN) :: PSCA_SW ! diffuse solar radiation (on horizontal surf.)
!                                             !                                       (W/m2)
REAL, DIMENSION(KSW),INTENT(IN)  :: PSW_BANDS ! mean wavelength of each shortwave band (m)
REAL, DIMENSION(KI), INTENT(IN)  :: PZENITH   ! zenithal angle at t       (radian from the vertical)
REAL, DIMENSION(KI), INTENT(IN)  :: PZENITH2  ! zenithal angle at t+1     (radian from the vertical)
REAL, DIMENSION(KI), INTENT(IN)  :: PAZIM     ! azimuthal angle      (radian from North, clockwise)
REAL, DIMENSION(KI), INTENT(IN)  :: PLW       ! longwave radiation (on horizontal surf.)
!                                             !                                       (W/m2)
REAL, DIMENSION(KI), INTENT(IN)  :: PPS       ! pressure at atmospheric model surface (Pa)
REAL, DIMENSION(KI), INTENT(IN)  :: PPA       ! pressure at forcing level             (Pa)
REAL, DIMENSION(KI), INTENT(IN)  :: PZS       ! atmospheric model CANOPY           (m)
REAL, DIMENSION(KI), INTENT(IN)  :: PCO2      ! CO2 concentration in the air          (kg/m3)
REAL, DIMENSION(KI), INTENT(IN)  :: PSNOW     ! snow precipitation                    (kg/m2/s)
REAL, DIMENSION(KI), INTENT(IN)  :: PRAIN     ! liquid precipitation                  (kg/m2/s)
!
!
REAL, DIMENSION(KI), INTENT(OUT) :: PSFTH     ! flux of heat                          (W/m2)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFTQ     ! flux of water vapor                   (kg/m2/s)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFU      ! zonal momentum flux                   (Pa)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFV      ! meridian momentum flux                (Pa)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFCO2    ! flux of CO2                           (kg/m2/s)
REAL, DIMENSION(KI,KSV),INTENT(OUT):: PSFTS   ! flux of scalar var.                   (kg/m2/s)
!
REAL, DIMENSION(KI), INTENT(OUT) :: PTRAD     ! radiative temperature                 (K)
REAL, DIMENSION(KI,KSW),INTENT(OUT):: PDIR_ALB! direct albedo for each spectral band  (-)
REAL, DIMENSION(KI,KSW),INTENT(OUT):: PSCA_ALB! diffuse albedo for each spectral band (-)
REAL, DIMENSION(KI), INTENT(OUT) :: PEMIS     ! emissivity                            (-)
!
REAL, DIMENSION(KI), INTENT(IN) :: PPEW_A_COEF! implicit coefficients
REAL, DIMENSION(KI), INTENT(IN) :: PPEW_B_COEF! needed if HCOUPLING='I'
REAL, DIMENSION(KI), INTENT(IN) :: PPET_A_COEF
REAL, DIMENSION(KI), INTENT(IN) :: PPEQ_A_COEF
REAL, DIMENSION(KI), INTENT(IN) :: PPET_B_COEF
REAL, DIMENSION(KI), INTENT(IN) :: PPEQ_B_COEF
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
! canopy turbulence scheme
!
REAL, DIMENSION(SIZE(PTA))        :: ZCANOPY   ! height of canopy   (m)
REAL, DIMENSION(SIZE(PTA))        :: ZSFLUX_U  ! Surface flux u'w' (m2/s2)
REAL, DIMENSION(SIZE(PTA))        :: ZSFLUX_T  ! Surface flux w'T' (mK/s)
REAL, DIMENSION(SIZE(PTA))        :: ZSFLUX_Q  ! Surface flux w'q' (kgm2/s)
REAL, DIMENSION(SIZE(PTA),NLVL)   :: ZFORC_U   ! tendency due to drag force for wind
REAL, DIMENSION(SIZE(PTA),NLVL)   :: ZDFORC_UDU! formal derivative of
!                                              ! tendency due to drag force for wind
REAL, DIMENSION(SIZE(PTA),NLVL)   :: ZFORC_E   ! tendency due to drag force for TKE
REAL, DIMENSION(SIZE(PTA),NLVL)   :: ZDFORC_EDE! formal derivative of
!                                              ! tendency due to drag force for TKE
REAL, DIMENSION(SIZE(PTA),NLVL)   :: ZFORC_T   ! tendency due to drag force for Temp
REAL, DIMENSION(SIZE(PTA),NLVL)   :: ZDFORC_TDT! formal derivative of
!                                              ! tendency due to drag force for Temp
REAL, DIMENSION(SIZE(PTA),NLVL)   :: ZFORC_Q   ! tendency due to drag force for Temp
REAL, DIMENSION(SIZE(PTA),NLVL)   :: ZDFORC_QDQ! formal derivative of
!                                              ! tendency due to drag force for hum.
REAL, DIMENSION(SIZE(PTA),NLVL)   :: ZLMO      ! MO length
REAL, DIMENSION(SIZE(PTA),NLVL)   :: ZLM       ! mixing length
REAL, DIMENSION(SIZE(PTA),NLVL)   :: ZLEPS     ! dissipative length
REAL, DIMENSION(SIZE(PTA))     :: ZH           ! canopy height (m)
REAL, DIMENSION(SIZE(PTA))     :: ZUSTAR       ! friction velocity including drag effect (m/s)
REAL, DIMENSION(SIZE(PTA))     :: ZUSTAR_GROUND! friction velocity at ground only (ISBA) (m/s)
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

CHARACTER(LEN=1) :: GCOUPLING
INTEGER          :: JLAYER    ! loop counter
INTEGER          :: JSWB      ! loop counter

REAL, DIMENSION(SIZE(PTA))   ::ZCANOPY_DENSITY
REAL, DIMENSION(SIZE(PTA))   ::ZUW_GROUND
REAL, DIMENSION(SIZE(PTA))   ::ZDUWDU_GROUND
REAL, DIMENSION(SIZE(PTA))   ::ZE_GROUND
REAL, DIMENSION(SIZE(PTA))   ::ZH_GROUND
!
REAL, DIMENSION(SIZE(PTA),NLVL)   :: ZZ        ! height above displacement height
!
INTEGER                      :: JJ
REAL, DIMENSION(SIZE(PTA))   :: ZU10
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------------------
!
!
!*      1.     Preliminary computations of the SBL scheme
!              ------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('COUPLING_ISBA_CANOPY_N',0,ZHOOK_HANDLE)
IF (LCANOPY) THEN
!
!*      1.1    Updates canopy vertical grid as a function of forcing height
!              ------------------------------------------------------------
!
!* determines where is the forcing level and modifies the upper levels of the canopy grid
!
  ZCANOPY = 0.
  CALL CANOPY_GRID_UPDATE(KI,NLVL,ZCANOPY,PUREF,XZ,XZF,XDZ,XDZF)
!
!
!
!*      1.2    Allocations and initialisations
!              -------------------------------
!
!
!       1.2.1  First time step canopy initialisation
!
  IF(ANY(XT(:,:) == XUNDEF)) THEN
    CALL INIT_ISBA_SBL(CISBA, CCPSURF, NLVL, PPA, PPS, PTA, PQA, PRHOA, PU, PV,           &
                       PDIR_SW, PSCA_SW, PSW_BANDS, PRAIN, PSNOW,                         &
                       PZREF, PUREF, XTG(:,1,:), XPATCH, XWG(:,1,:), XWGI(:,1,:),         &
                       XZ0, XSSO_SLOPE, XRESA, XVEG, XLAI, XWR, XRGL, XRSMIN,             &
                       XGAMMA, XWRMAX_CF, XZ0_O_Z0H, XWFC, XWSAT, TSNOW, XZ,              &
                       XT, XQ, XU, XTKE, XP)
  ENDIF
!
!*      1.3    Allocations
!              -----------
!
  ZFORC_U    = 0.
  ZDFORC_UDU = 0.
!
  ZFORC_E(:,:) = 0.
  ZDFORC_EDE(:,:) = 0.
!
  ZFORC_T(:,:) = 0.
  ZDFORC_TDT(:,:) = 0.
!
  ZFORC_Q(:,:) = 0.
  ZDFORC_QDQ(:,:) = 0.
!
  ZSFLUX_U = 0.
  ZSFLUX_T = 0.
  ZSFLUX_Q = 0.
!
  ZLMO = SPREAD(XLMO,2,NLVL)
!
!* default :
!* no canopy in ISBA scheme
!
  ZH = 0.
!
!
!* determine for each level the height above displacement height
!
  ZZ(:,:) = XZ(:,:)
!
!*      1.4   canopy for wind drag only
!             -------------------------
!
  IF (LCANOPY_DRAG) THEN
!* mean canopy height
!
!* in ecoclimap, height is set retrieved from roughness length (z0/0.13)
    ZH = SUM(XPATCH(:,:)*XZ0(:,:)/0.13,DIM=2)
    ZH = MIN(ZH, XZF(:,NLVL))
    WHERE (ZH<=XDZ(:,1)) ZH = 0.
!
!* canopy for wind drag only
    ZCANOPY_DENSITY=SUM(XPATCH(:,:)*XLAI(:,:),DIM=2)
    ZUW_GROUND=0.
    ZDUWDU_GROUND = 0.
!
!* computes tendencies on wind and Tke due to canopy
    CALL ISBA_CANOPY(KI,NLVL,XZ,XZF,XDZ,XDZF,ZH,ZCANOPY_DENSITY,XU,XTKE,    &
                    ZUW_GROUND, ZDUWDU_GROUND,                              &
                    ZFORC_U,ZDFORC_UDU,ZFORC_E,ZDFORC_EDE                   )
!
  ENDIF
!
!*      1.4   Subgrid-scale orographic drag (Beljaars et al 2004)
!             -----------------------------
!
  IF (CROUGH=='BE04') THEN
!
!* computes tendencies on wind and Tke due to subgridscale orography
    CALL SSO_BELJAARS04(KI,NLVL,XZ,XZF,XDZ,XDZF,XSSO_STDEV,XU,XTKE,    &
                    ZFORC_U,ZDFORC_UDU                                 )
!
  ENDIF
!
!
!*      1.5   Computes coefficients for implicitation
!             ---------------------------------------
!
  ZWIND = SQRT(PU**2+PV**2)
  CALL CANOPY_EVOL(KI,NLVL,PTSTEP,1,ZZ,ZWIND,PTA,PQA,PPA,PRHOA,             &
                   ZSFLUX_U,ZSFLUX_T,ZSFLUX_Q,                              &
                   ZFORC_U,ZDFORC_UDU,ZFORC_E,ZDFORC_EDE,                   &
                   ZFORC_T,ZDFORC_TDT,ZFORC_Q,ZDFORC_QDQ,                   &
                   XZ,XZF,XDZ,XDZF,XU,XTKE,XT,XQ,ZLMO,ZLM,ZLEPS,XP,ZUSTAR,  &
                   ZALFAU,ZBETAU,ZALFATH,ZBETATH,ZALFAQ,ZBETAQ              )
!
!*     1.6     Goes from atmospheric forcing to canopy forcing height
!              ------------------------------------------------------
!
  GCOUPLING ='I'
  ZTA = XT(:,1)
  ZQA = XQ(:,1)
  ZU  = PU / MAX(SQRT(PU**2+PV**2),XWINDMIN) * XU(:,1)
  ZV  = PV / MAX(SQRT(PU**2+PV**2),XWINDMIN) * XU(:,1)
  ZZREF = XZ(:,1)
  ZUREF = XZ(:,1)
  ZPA   = XP(:,1)
  ZEXNA(:)   = (XP(:,1)/XP00)**(XRD/XCPD)
!* ALMA conventions for implicit coefficients:
! U+  = - rho A u'w'  + B
! Th+ = - rho A w'th' + B
! q+  = - rho A w'q'  + B
!
  ZPEW_A_COEF = - ZALFAU / PRHOA
  ZPEW_B_COEF = ZBETAU
  ZPET_A_COEF = - ZALFATH / PRHOA
  ZPET_B_COEF = ZBETATH
  ZPEQ_A_COEF = - ZALFAQ / PRHOA
  ZPEQ_B_COEF = ZBETAQ
!
!-------------------------------------------------------------------------------------
ELSE
!-------------------------------------------------------------------------------------
!
!*      2.     If no canopy scheme is used, forcing is not modified
!              ----------------------------------------------------
!
  GCOUPLING = HCOUPLING
  ZTA = PTA
  ZQA = PQA
  ZU  = PU
  ZV  = PV
  ZUREF=PUREF
  ZZREF=PZREF
  ZPA = PPA
  ZPEQ_A_COEF(:) = PPEQ_A_COEF(:)
  ZPEQ_B_COEF(:) = PPEQ_B_COEF(:)
  ZPET_A_COEF(:) = PPET_A_COEF(:)
  ZPET_B_COEF(:) = PPET_B_COEF(:)
  ZPEW_A_COEF(:) = PPEW_A_COEF(:)
  ZPEW_B_COEF(:) = PPEW_B_COEF(:)
END IF
!
!-------------------------------------------------------------------------------------
!
!*      2.     Call of ISBA
!              ------------
!
CALL COUPLING_ISBA_n(HPROGRAM, GCOUPLING,                                                  &
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
!
!-------------------------------------------------------------------------------------
!
!*      3.     End if no canopy is used
!              ------------------------
!
IF (.NOT. LCANOPY .AND. LHOOK) CALL DR_HOOK('COUPLING_ISBA_CANOPY_N',1,ZHOOK_HANDLE)
IF (.NOT. LCANOPY) RETURN
!
!-------------------------------------------------------------------------------------
!
!*      4.     Computes the impact of surface on air
!              -------------------------------------
!
!
ZFORC_U    = 0.
ZDFORC_UDU = 0.
ZSFLUX_U = - SQRT(PSFU(:)**2+PSFV(:)**2) / PRHOA(:)
!
ZFORC_E(:,:) = 0.
ZDFORC_EDE(:,:) = 0.
!
ZFORC_T(:,:) = 0.
ZDFORC_TDT(:,:) = 0.
ZSFLUX_T(:) = PSFTH(:) / XCPD * ZEXNA(:) / PRHOA(:)
!
ZFORC_Q(:,:) = 0.
ZDFORC_QDQ(:,:) = 0.
ZSFLUX_Q(:) = PSFTQ(:)
!-------------------------------------------------------------------------------------
!
!*      5.     Computes the impact of canopy on air
!              ------------------------------------
!
IF (LCANOPY_DRAG) THEN
  ZUW_GROUND=-SQRT(PSFU**2+PSFV**2)/ PRHOA(:)
  ZDUWDU_GROUND = 0.
  WHERE (XU(:,1)   /=0.) ZDUWDU_GROUND = 2. * ZUW_GROUND / XU(:,1)

!* computes tendencies on wind and Tke due to canopy and surface
  CALL ISBA_CANOPY(KI,NLVL,XZ,XZF,XDZ,XDZF,ZH,ZCANOPY_DENSITY,XU,XTKE,  &
                  ZUW_GROUND, ZDUWDU_GROUND,                            &
                  ZFORC_U,ZDFORC_UDU,ZFORC_E,ZDFORC_EDE                 )

  ZSFLUX_U = 0.  ! surface friction is incorporated in ZFORC_U by ISBA_CANOPY routine
END IF
!
!
IF (CROUGH=='BE04') THEN
!
!* computes tendencies on wind and Tke due to subgridscale orography
  CALL SSO_BELJAARS04(KI,NLVL,XZ,XZF,XDZ,XDZF,XSSO_STDEV,XU,XTKE,    &
                      ZFORC_U,ZDFORC_UDU                             )
!
ENDIF
!
!-------------------------------------------------------------------------------------
!
!*      6.    Evolution of canopy air due to these impacts
!             --------------------------------------------
!
ZWIND = SQRT(PU**2+PV**2)
CALL CANOPY_EVOL(KI,NLVL,PTSTEP,2,ZZ,ZWIND,PTA,PQA,PPA,PRHOA,                 &
                 ZSFLUX_U,ZSFLUX_T,ZSFLUX_Q,                                  &
                 ZFORC_U,ZDFORC_UDU,ZFORC_E,ZDFORC_EDE,                       &
                 ZFORC_T,ZDFORC_TDT,ZFORC_Q,ZDFORC_QDQ,                       &
                 XZ,XZF,XDZ,XDZF,XU,XTKE,XT,XQ,ZLMO,ZLM,ZLEPS,XP,ZUSTAR,       &
                 ZALFAU,ZBETAU,ZALFATH,ZBETATH,ZALFAQ,ZBETAQ                  )
!
XLMO(:) = ZLMO(:,NLVL)
!
! Momentum fluxes if canopy is used
!
!* Total friction due to surface averaged friction and averaged canopy drag
IF (LCANOPY_DRAG .OR. CROUGH=='BE04') THEN
  ZUSTAR_GROUND=SQRT(SQRT(PSFU**2+PSFV**2)/PRHOA)
  WHERE (ZUSTAR_GROUND(:)>0.)
    PSFU(:) = PSFU(:) * ZUSTAR**2/ZUSTAR_GROUND**2
    PSFV(:) = PSFV(:) * ZUSTAR**2/ZUSTAR_GROUND**2
  END WHERE
END IF
!
!-------------------------------------------------------------------------------------
!
!*      7.    2m and 10m diagnostics if canopy is used
!             ----------------------------------------
!
!
IF (N2M>=1) THEN
  XAVG_T2M(:) = XT(:,2)
  XAVG_Q2M(:) = XQ(:,2) / PRHOA(:)
  XAVG_HU2M(:)= MIN( XAVG_Q2M(:) / QSAT(XAVG_T2M(:),XP(:,2)) , 1.)
  CALL INTERPOL_SBL(XZ(:,:),XU(:,:),10.,ZU10(:))
  DO JJ=1,SIZE(PTA)
    IF (ZU10(JJ)/=XUNDEF) THEN
      IF (ZWIND(JJ)>0.) THEN
        XAVG_ZON10M(JJ) = ZU10(JJ) * PU(JJ)/ZWIND(JJ)
        XAVG_MER10M(JJ) = ZU10(JJ) * PV(JJ)/ZWIND(JJ)
      ELSE
        XAVG_ZON10M(JJ) = 0.
        XAVG_MER10M(JJ) = 0.
      END IF
    ELSE
      XAVG_ZON10M(JJ) = XUNDEF
      XAVG_MER10M(JJ) = XUNDEF
    END IF
  END DO  
END IF
!
!* Total friction due to surface averaged friction and averaged canopy drag
IF ((LCANOPY_DRAG .OR. CROUGH=='BE04') .AND. LSURF_BUDGET) THEN
  XAVG_FMU = PSFU
  XAVG_FMV = PSFV
END IF
IF (LHOOK) CALL DR_HOOK('COUPLING_ISBA_CANOPY_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE COUPLING_ISBA_CANOPY_n
