!     #########
    SUBROUTINE GARDEN (DTCO, TG, T, TOP, DTGR, DTI, GB, VD, TV, TIR,            &
                       HIMPLICIT_WIND, TPTIME, PTSUN, PPEW_A_COEF, PPEW_B_COEF, &
                       PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,      &
                       PTSTEP, PZ_LOWCAN, PT_LOWCAN, PQ_LOWCAN, PEXNS, PRHOA,   &
                       PCO2, PPS, PRR, PSR, PZENITH, PSW, PLW, PU_LOWCAN,       &
                       PALBNIR_TVEG, PALBVIS_TVEG, PALBNIR_TSOIL, PALBVIS_TSOIL,&                
                       PRN_GARDEN, PH_GARDEN, PLE_GARDEN, PGFLUX_GARDEN, PSFCO2,&
                       PEVAP_GARDEN, PUW_GARDEN, PRUNOFF_GARDEN, PAC_GARDEN,    &
                       PQSAT_GARDEN, PAC_AGG_GARDEN, PHU_AGG_GARDEN,&
                       PDRAIN_GARDEN, PIRRIG_GARDEN         )  
!   ##########################################################################
!
!!****  *GARDEN*  
!!
!!    PURPOSE
!!    -------
!
!!call the vegetation scheme (ISBA) inside TEB
!     
!!**  METHOD
!     ------
!
!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!      A. Lemonsu          * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!    Original    05/2009
!     B. decharme 04/2013 : variables for surf/atm coupling
!                           dummy for water table / surface coupling
!!    P. Samuelsson  10/2014  Introduced dummy variables in call to ISBA for MEB
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_GRID_n, ONLY : GRID_t
USE MODD_SSO_n, ONLY : SSO_t, SSO_INIT
USE MODD_TEB_n, ONLY : TEB_1P_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
!
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_GR_BIOG_n, ONLY : GR_BIOG_t
USE MODD_SURFEX_n, ONLY : TEB_VEG_DIAG_t
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_t
USE MODD_TEB_IRRIG_n, ONLY : TEB_IRRIG_t
!
USE MODD_AGRI_n, ONLY : AGRI_t,AGRI_INIT
!
USE MODD_TYPE_DATE_SURF,    ONLY: DATE_TIME
USE MODD_SURF_PAR,          ONLY: XUNDEF
USE MODD_CSTS,              ONLY: XCPD
!
!
USE MODI_ISBA
USE MODI_VEGETATION_UPDATE
USE MODE_THERMOS
!
USE MODI_FLAG_TEB_GARDEN_n
USE MODI_CARBON_EVOL
USE MODI_VEGETATION_EVOL
USE MODI_TEB_IRRIG
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    Declarations of arguments
!
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(GRID_t), INTENT(INOUT) :: TG
TYPE(TEB_1P_t), INTENT(INOUT) :: T
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTGR
!
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTI
TYPE(GR_BIOG_t), INTENT(INOUT) :: GB
TYPE(TEB_VEG_DIAG_t), INTENT(INOUT) :: VD
TYPE(TEB_VEG_t), INTENT(INOUT) :: TV
TYPE(TEB_IRRIG_t), INTENT(INOUT) :: TIR
!
 CHARACTER(LEN=*),     INTENT(IN)  :: HIMPLICIT_WIND   ! wind implicitation option
!                                                     ! 'OLD' = direct
!                                                     ! 'NEW' = Taylor serie, order 1
TYPE(DATE_TIME)     , INTENT(IN)    :: TPTIME             ! current date and time from teb
REAL, DIMENSION(:)  , INTENT(IN)    :: PTSUN              ! solar time      (s from midnight)
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEW_A_COEF        ! implicit coefficients
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEW_B_COEF        ! for wind coupling
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEQ_A_COEF        ! implicit coefficients
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEQ_B_COEF        ! for humidity
REAL, DIMENSION(:)  , INTENT(IN)    :: PPET_A_COEF        ! implicit coefficients
REAL, DIMENSION(:)  , INTENT(IN)    :: PPET_B_COEF        ! for temperature
REAL                , INTENT(IN)    :: PTSTEP             ! time step
REAL, DIMENSION(:)  , INTENT(IN)    :: PZ_LOWCAN          ! height of atm. var. near the road
REAL, DIMENSION(:)  , INTENT(IN)    :: PT_LOWCAN          ! temp. near the road
REAL, DIMENSION(:)  , INTENT(IN)    :: PQ_LOWCAN          ! hum. near the road
REAL, DIMENSION(:)  , INTENT(IN)    :: PPS                ! pressure at the surface
REAL, DIMENSION(:)  , INTENT(IN)    :: PEXNS              ! surface exner function
REAL, DIMENSION(:)  , INTENT(IN)    :: PRHOA              ! air density at the lowest level
REAL, DIMENSION(:)  , INTENT(IN)    :: PCO2               ! CO2 concentration in the air    (kg/m3)
REAL, DIMENSION(:)  , INTENT(IN)    :: PRR                ! rain rate
REAL, DIMENSION(:)  , INTENT(IN)    :: PSR                ! snow rate
REAL, DIMENSION(:)  , INTENT(IN)    :: PZENITH            ! solar zenithal angle
REAL, DIMENSION(:)  , INTENT(IN)    :: PSW                ! incoming total solar rad on an horizontal surface
REAL, DIMENSION(:)  , INTENT(IN)    :: PLW                ! atmospheric infrared radiation
REAL, DIMENSION(:)  , INTENT(IN)    :: PU_LOWCAN          ! wind near the road
REAL, DIMENSION(:)  , INTENT(IN)    :: PALBNIR_TVEG       ! nearIR  veg tot albedo
REAL, DIMENSION(:)  , INTENT(IN)    :: PALBVIS_TVEG       ! visible veg tot albedo
REAL, DIMENSION(:)  , INTENT(IN)    :: PALBNIR_TSOIL      ! nearIR  soil tot albedo
REAL, DIMENSION(:)  , INTENT(IN)    :: PALBVIS_TSOIL      ! visible soil tot albedo
!
!
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN_GARDEN         ! net radiation over green areas
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_GARDEN          ! sensible heat flux over green areas
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_GARDEN         ! latent heat flux over green areas
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLUX_GARDEN      ! flux through the green areas
REAL, DIMENSION(:)  , INTENT(OUT)   :: PSFCO2             ! flux of CO2 positive toward the atmosphere (m/s*kg_CO2/kg_air)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PEVAP_GARDEN       ! total evaporation over gardens (kg/m2/s)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PUW_GARDEN         ! friction flux (m2/s2)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRUNOFF_GARDEN     ! runoff over garden (kg/m2/s)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_GARDEN         ! aerodynamical conductance
REAL, DIMENSION(:)  , INTENT(OUT)   :: PQSAT_GARDEN       ! saturation humidity
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_AGG_GARDEN     ! aggreg. aeodynamic resistance for green areas for latent heat flux
REAL, DIMENSION(:)  , INTENT(OUT)   :: PHU_AGG_GARDEN     ! aggreg. relative humidity for green areas for latent heat flux
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDRAIN_GARDEN      ! garden total (vertical) drainage
REAL, DIMENSION(:)  , INTENT(OUT)   :: PIRRIG_GARDEN      ! garden summer irrigation rate
!
!
!*      0.2    Declarations of local variables
!
TYPE(SSO_t) :: YGDSS
TYPE(AGRI_t) :: YAG
!
REAL, DIMENSION(SIZE(PPS)) :: ZDIRCOSZW           ! orography slope cosine (=1 in TEB)
REAL, DIMENSION(SIZE(PPS),TV%O%NNBIOMASS) :: ZRESP_BIOMASS_INST       ! instantaneous biomass respiration (kgCO2/kgair m/s)
REAL, DIMENSION(SIZE(PPS)) :: ZUSTAR
!
!  temperatures
!
REAL, DIMENSION(SIZE(PPS)) :: ZTA ! estimate of air temperature at future time
!                                 ! step as if modified by ISBA flux alone.
REAL, DIMENSION(SIZE(PPS)) :: ZDEEP_FLUX ! heat flux at base of the deep soil
!
!  surfaces relative fractions
!  for flood
REAL, DIMENSION(SIZE(PPS)) :: ZEMISF
!
!  variables for deep soil temperature
REAL, DIMENSION(SIZE(PPS)) :: ZTDEEP_A
!
! Dummy variables for MEB:
REAL, DIMENSION(SIZE(PPS)) :: ZP_MEB_SCA_SW, ZPALPHAN, ZZ0G_WITHOUT_SNOW, &
                              ZZ0_MEBV, ZZ0H_MEBV, ZZ0EFF_MEBV, ZZ0_MEBN, &
                              ZZ0H_MEBN, ZZ0EFF_MEBN
INTEGER                    :: ILU
LOGICAL :: GMASK, GALB
LOGICAL :: GUPDATED
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*      1.     various initialisations
!              -----------------------
!
IF (LHOOK) CALL DR_HOOK('GARDEN',0,ZHOOK_HANDLE)
ILU = SIZE(PPS)
!
ZDIRCOSZW = 1.
!
CALL SSO_INIT(YGDSS)
!
CALL AGRI_INIT(YAG)
!
!-------------------------------------------------------------------------------
!
!*      2.     Treatment of green areas
!              ------------------------
!*      2.1    Automatic irrigation
!              --------------------
!
CALL TEB_IRRIG(TIR%LPAR_GD_IRRIG, PTSTEP, TPTIME%TDATE%MONTH, PTSUN,        &
               TIR%XGD_START_MONTH, TIR%XGD_END_MONTH, TIR%XGD_START_HOUR,  &
               TIR%XGD_END_HOUR, TIR%XGD_24H_IRRIG, PIRRIG_GARDEN           ) 
!
! --------------------------------------------------------------------------------------
! Vegetation update (in case of non-interactive vegetation):
! --------------------------------------------------------------------------------------
!
TV%I%TTIME = TPTIME
TV%O%LECOCLIMAP = (.NOT. TV%O%LPAR)
!
GUPDATED=.FALSE.
GALB = .FALSE. 
IF (TV%O%CPHOTO=='LAI'.OR.TV%O%CPHOTO=='LST'.OR.TV%O%CPHOTO=='NIT'.OR.TV%O%CPHOTO=='NCB') GALB = .TRUE.
!
CALL VEGETATION_UPDATE(DTCO, DTI, TG%NDIM, TV%O, TV%M%T%CUR, TV%M%M, TV%M%I, TV%M%A, &
                       PTSTEP, TV%I%TTIME, TOP%XCOVER, TOP%LCOVER, .FALSE.,'GRD',  &
                       GALB, YGDSS, GUPDATED, OABSENT=(T%XGARDEN==0.)   )
!
!
VD%DP%XZ0(:) = TV%M%T%CUR%XZ0(:,1)
VD%DP%XZ0H(:) = TV%M%T%CUR%XZ0(:,1) / TV%M%X%XZ0_O_Z0H(:,1)
!
VD%DP%XZ0EFF(:) =  TV%M%T%CUR%XZ0(:,1)
!
!*      2.2    Call ISBA for green areas
!              -------------------------
!
ALLOCATE(GB%XIACAN(SIZE(PPS),SIZE(TV%IP%XABC),1))
!
 CALL ISBA(TV%O, TV%M%X, TV%M%T%CUR, TV%M%M, TV%M%I, TV%P, TV%IP, TV%I, TV%R%CUR, &
           TG, YAG, VD%D, VD%DP, VD%E, VD%EP, VD%M, TV%R%CUR%TSNOW%SCHEME, TPTIME,&
           TV%IP%XPOI, TV%IP%XABC, GB%XIACAN(:,:,1), .FALSE., PTSTEP,             &
           HIMPLICIT_WIND, PZ_LOWCAN, PZ_LOWCAN, ZDIRCOSZW, PT_LOWCAN, PQ_LOWCAN, &
           PEXNS, PRHOA, PPS, PEXNS,  PRR, PSR, PZENITH, ZP_MEB_SCA_SW, PSW, PLW, &
           PU_LOWCAN, PPEW_A_COEF, PPEW_B_COEF, PPET_A_COEF, PPEQ_A_COEF,         &
           PPET_B_COEF, PPEQ_B_COEF, PALBNIR_TVEG, PALBVIS_TVEG, PALBNIR_TSOIL,   &
           PALBVIS_TSOIL, ZPALPHAN, ZZ0G_WITHOUT_SNOW, ZZ0_MEBV, ZZ0H_MEBV,       &
           ZZ0EFF_MEBV, ZZ0_MEBN, ZZ0H_MEBN, ZZ0EFF_MEBN, ZTDEEP_A, PCO2,         &
           TV%I%XFFG(:,1), TV%I%XFFV(:,1), ZEMISF, ZUSTAR, PAC_AGG_GARDEN,        &
           PHU_AGG_GARDEN, ZRESP_BIOMASS_INST, ZDEEP_FLUX, PIRRIG_GARDEN )                                                  !
PRUNOFF_GARDEN(:) = VD%EP%XRUNOFF(:)
PDRAIN_GARDEN (:) = VD%EP%XDRAIN(:)
!
IF (TV%R%CUR%TSNOW%SCHEME=='3-L' .OR. TV%R%CUR%TSNOW%SCHEME=='CRO') &
    TV%R%CUR%TSNOW%TS(:,1)= VD%M%XSNOWTEMP(:,1)
!
IF (TV%O%LTR_ML) THEN
  GMASK = ( TPTIME%TIME - PTSTEP < 0. ) .AND. ( TPTIME%TIME >= 0. )
  IF (GMASK) THEN
    ALLOCATE(VD%M%XDFAPARC(ILU),VD%M%XDFAPIRC(ILU),VD%M%XDLAI_EFFC(ILU))
    VD%M%XDFAPARC  (:) = TV%R%CUR%XFAPARC   (:,1) / TV%R%CUR%XMUS (:,1)
    VD%M%XDFAPIRC  (:) = TV%R%CUR%XFAPIRC   (:,1) / TV%R%CUR%XMUS (:,1)
    VD%M%XDLAI_EFFC(:) = TV%R%CUR%XLAI_EFFC (:,1) / TV%R%CUR%XMUS (:,1)
  ENDIF
ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Vegetation evolution for interactive LAI
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
IF (TV%O%CPHOTO=='LAI' .OR. TV%O%CPHOTO=='LST' .OR. TV%O%CPHOTO=='NIT') THEN
  CALL VEGETATION_EVOL(TV%O, TV%IP, TV%M%X, TV%M%T%CUR, TV%M%A, TV%M%I, TV%R%CUR, &
                       .FALSE., PTSTEP, TPTIME%TDATE%MONTH, TPTIME%TDATE%DAY,     &
                       TPTIME%TIME, TG%XLAT, PRHOA, PCO2, YGDSS, ZRESP_BIOMASS_INST )         
END IF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Diagnostic of respiration carbon fluxes and soil carbon evolution
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
PSFCO2(:)=0.
VD%EP%XRESP_ECO (:)=0.
VD%EP%XRESP_AUTO(:)=0.
!
IF (TV%O%CPHOTO/='NON' .AND. TV%O%CRESPSL/='NON' .AND. ANY(TV%M%T%CUR%XLAI(:,1)/=XUNDEF)) THEN
  CALL CARBON_EVOL(TV%O, TV%P, TV%IP, TV%M%X, TV%M%T%CUR, TV%R%CUR, VD%EP, &
                   PTSTEP, PRHOA, ZRESP_BIOMASS_INST  )
  ! calculation of vegetation CO2 flux
  PSFCO2(:) = - VD%EP%XGPP(:) + VD%EP%XRESP_ECO(:)
END IF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!*      4.     Set undefined values for points where there is no garden
!              --------------------------------------------------------
!
! This way, these points are clearly flaged, and one will not try to interpret
! the values for those points
!
 CALL FLAG_TEB_GARDEN_n(TV%R%CUR, TV%O, TV%M%T%CUR%XLAI, T%XGARDEN, 2)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!*      9.     Fields required for TEB
!              -----------------------
!
WHERE (T%XGARDEN/=0.)
  !
  ! energy balance
  !
   PRN_GARDEN    (:) = VD%DP%XRN    (:)
   PH_GARDEN     (:) = VD%DP%XH     (:)
   PLE_GARDEN    (:) = TV%R%CUR%XLE(:,1)
   PGFLUX_GARDEN (:) = VD%DP%XGFLUX (:)
   PEVAP_GARDEN  (:) = VD%DP%XEVAP  (:)
  !
  !
  ! Estimate of green area aerodynamic conductance recomputed from heat flux,
  ! surface (radiative) temp. and forcing air temperature (estimated at future time step)
  ZTA = PPET_B_COEF + PPET_A_COEF * PH_GARDEN
  PAC_GARDEN = 0.
  WHERE (VD%DP%XTSRAD /= ZTA)
    PAC_GARDEN(:)   = MAX(PH_GARDEN(:) / XCPD / PRHOA(:) / (VD%DP%XTSRAD - ZTA) , 0.)
  ENDWHERE
  !
  ! Humidity of saturation for green areas
  PQSAT_GARDEN(:) = QSAT(TV%R%CUR%XTG(:,1,1),PPS(:))
  !
  !* friction flux
  PUW_GARDEN(:)    = -ZUSTAR(:)**2
  !
ELSEWHERE
  !
  PRN_GARDEN    (:) = XUNDEF
  PH_GARDEN     (:) = XUNDEF
  PLE_GARDEN    (:) = XUNDEF
  PGFLUX_GARDEN (:) = XUNDEF
  PEVAP_GARDEN  (:) = XUNDEF
  PAC_GARDEN    (:) = XUNDEF
  PQSAT_GARDEN  (:) = XUNDEF
  PUW_GARDEN    (:) = XUNDEF
  !
END WHERE
!
IF (LHOOK) CALL DR_HOOK('GARDEN',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE GARDEN
