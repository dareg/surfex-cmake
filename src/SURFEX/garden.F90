!     #########
    SUBROUTINE GARDEN (DTCO, TG, T, TOP, DTGR, GDM, &
                       HIMPLICIT_WIND, TPTIME, PTSUN, PPEW_A_COEF, PPEW_B_COEF,       &
                PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,                  &
                PTSTEP, PZ_LOWCAN,                                                   &
                PT_LOWCAN, PQ_LOWCAN, PEXNS, PRHOA, PCO2, PPS, PRR, PSR, PZENITH,    &
                PSW, PLW, PU_LOWCAN,                                                 &
                PALBNIR_TVEG, PALBVIS_TVEG, PALBNIR_TSOIL, PALBVIS_TSOIL,            &                
                PRN_GARDEN,PH_GARDEN,PLE_GARDEN,PGFLUX_GARDEN,PSFCO2,                &
                PEVAP_GARDEN, PUW_GARDEN,PRUNOFF_GARDEN,                             &
                PAC_GARDEN,PQSAT_GARDEN,PTS_GARDEN,                                  &
                PAC_AGG_GARDEN, PHU_AGG_GARDEN, PDRAIN_GARDEN, PIRRIG_GARDEN         )  
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
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_SURFEX_n, ONLY : TEB_GARDEN_MODEL_t
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
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTGR
TYPE(TEB_GARDEN_MODEL_t), INTENT(INOUT) :: GDM
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
REAL, DIMENSION(:)  , INTENT(INOUT) :: PTS_GARDEN         ! radiative surface temp. (snow free)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_AGG_GARDEN     ! aggreg. aeodynamic resistance for green areas for latent heat flux
REAL, DIMENSION(:)  , INTENT(OUT)   :: PHU_AGG_GARDEN     ! aggreg. relative humidity for green areas for latent heat flux
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDRAIN_GARDEN      ! garden total (vertical) drainage
REAL, DIMENSION(:)  , INTENT(OUT)   :: PIRRIG_GARDEN      ! garden summer irrigation rate
!
!
!*      0.2    Declarations of local variables
!
 CHARACTER(LEN=3)     :: HRAIN      ! Rainfall spatial distribution
                                   ! 'DEF' = No rainfall spatial distribution
                                   ! 'SGH' = Rainfall exponential spatial distribution
LOGICAL              :: OFLOOD     ! Activation of the flooding scheme
LOGICAL              :: OTEMP_ARP  ! True  = time-varying force-restore soil temperature (as in ARPEGE)
                                   ! False = No time-varying force-restore soil temperature (Default)
LOGICAL              :: OGLACIER   ! True = Over permanent snow and ice, 
!                                      initialise WGI=WSAT,
!                                      Hsnow>=10m and allow 0.8<SNOALB<0.85
                                   ! False = No specific treatment
REAL, DIMENSION(0)   ::  ZSODELX   ! Pulsation for each layer (Only used if LTEMP_ARP=True)     
!
REAL, DIMENSION(SIZE(PPS))               :: ZMUF    ! fraction of the grid cell reached by the rainfall
REAL, DIMENSION(SIZE(PPS))               :: ZFSAT   ! Topmodel (SGH not used in TEB) saturated fraction
REAL, DIMENSION(SIZE(PPS),GDM%TV%O%NGROUND_LAYER) :: ZTOPQS  ! Topmodel (SGH not used in TEB) lateral subsurface flow by layer
REAL, DIMENSION(SIZE(PPS))               :: ZQSB    ! Topmodel (SGH not used in TEB) output lateral subsurface
REAL, DIMENSION(SIZE(PPS))               :: ZFWTD   ! grid-cell fraction of water table to rise
REAL, DIMENSION(SIZE(PPS))               :: ZWTD    ! water table depth from Obs, TRIP or MODCOU
!
REAL, DIMENSION(SIZE(PPS)) :: ZDIRCOSZW           ! orography slope cosine (=1 in TEB)
REAL, DIMENSION(SIZE(PPS),GDM%TV%O%NNBIOMASS) :: ZRESP_BIOMASS_INST       ! instantaneous biomass respiration (kgCO2/kgair m/s)
!
!  temperatures
!
REAL, DIMENSION(SIZE(PPS)) :: ZTA ! estimate of air temperature at future time
!                                 ! step as if modified by ISBA flux alone.
REAL, DIMENSION(SIZE(PPS)) :: ZDEEP_FLUX ! heat flux at base of the deep soil
!
!   desactivated diag
!
REAL, DIMENSION(SIZE(PPS)) :: ZRN_ISBA      ! net radiative flux from snow-free surface 
REAL, DIMENSION(SIZE(PPS)) :: ZH_ISBA       ! sensible heat flux from snow-free surface 
REAL, DIMENSION(SIZE(PPS)) :: ZLEI_ISBA     ! baresoil evaporation from snow-free surface 
REAL, DIMENSION(SIZE(PPS)) :: ZLEG_ISBA     ! baresoil evaporation from snow-free surface 
REAL, DIMENSION(SIZE(PPS)) :: ZLEGI_ISBA    ! baresoil sublimation from snow-free surface 
REAL, DIMENSION(SIZE(PPS)) :: ZLEV_ISBA     ! total evapotranspiration from vegetation over 
REAL, DIMENSION(SIZE(PPS)) :: ZLETR_ISBA    ! transpiration from snow-free surface 
REAL, DIMENSION(SIZE(PPS)) :: ZUSTAR_ISBA   ! friction velocity from snow-free 
REAL, DIMENSION(SIZE(PPS)) :: ZLER_ISBA     ! evaporation from canopy water interception 
REAL, DIMENSION(SIZE(PPS)) :: ZLE_ISBA      ! latent heat flux from snow-free surface 
REAL, DIMENSION(SIZE(PPS)) :: ZGFLUX_ISBA   ! net energy flux into the snow-free surface 
REAL, DIMENSION(SIZE(PPS)) :: ZRNSNOW       ! net radiative flux from snow (ISBA-ES:3-L)    (W/m2)
REAL, DIMENSION(SIZE(PPS)) :: ZHSNOW        ! sensible heat flux from snow (ISBA-ES:3-L)    (W/m2)
REAL, DIMENSION(SIZE(PPS)) :: ZHPSNOW       ! heat release from rainfall (ISBA-ES:3-L)      (W/m2)
REAL, DIMENSION(SIZE(PPS)) :: ZGFLUXSNOW    ! net surface energy flux into snowpack (ISBA-ES:3-L)(W/m2)
REAL, DIMENSION(SIZE(PPS)) :: ZUSTARSNOW    ! friction velocity  over snow (ISBA-ES:3-L)    (m/s)
REAL, DIMENSION(SIZE(PPS)) :: ZGRNDFLUX     ! soil/snow interface heat flux (ISBA-ES:3-L)   (W/m2)
REAL, DIMENSION(SIZE(PPS)) :: ZSRSFC        ! snowfall over snowpack (ISBA-ES:3-L)          (kg/m2/s)
REAL, DIMENSION(SIZE(PPS)) :: ZRRSFC        ! rainfall over snowpack (ISBA-ES:3-L)          (kg/m2/s)
REAL, DIMENSION(SIZE(PPS)) :: ZLESL         ! snowpack evaporation (ISBA-ES:3-L)            (W/m2)
REAL, DIMENSION(SIZE(PPS)) :: ZCDSNOW       ! snow drag coefficient (ISBA-ES:3-L)           (-)
REAL, DIMENSION(SIZE(PPS)) :: ZCHSNOW       ! heat turbulent transfer coefficient           (-)
!
REAL, DIMENSION(SIZE(PPS)) :: ZCG
REAL, DIMENSION(SIZE(PPS)) :: ZC1
REAL, DIMENSION(SIZE(PPS)) :: ZC2
REAL, DIMENSION(SIZE(PPS)) :: ZWGEQ
REAL, DIMENSION(SIZE(PPS)) :: ZCT
REAL, DIMENSION(SIZE(PPS)) :: ZRS
REAL, DIMENSION(SIZE(PPS)) :: ZCH
REAL, DIMENSION(SIZE(PPS)) :: ZCD
REAL, DIMENSION(SIZE(PPS)) :: ZCDN
REAL, DIMENSION(SIZE(PPS)) :: ZRI
REAL, DIMENSION(SIZE(PPS)) :: ZHU
REAL, DIMENSION(SIZE(PPS)) :: ZHUG
REAL, DIMENSION(SIZE(PPS)) :: ZRN
REAL, DIMENSION(SIZE(PPS)) :: ZH
REAL, DIMENSION(SIZE(PPS)) :: ZLEI
REAL, DIMENSION(SIZE(PPS)) :: ZLEG
REAL, DIMENSION(SIZE(PPS)) :: ZLEGI
REAL, DIMENSION(SIZE(PPS)) :: ZLEV
REAL, DIMENSION(SIZE(PPS)) :: ZLES
REAL, DIMENSION(SIZE(PPS)) :: ZLER
REAL, DIMENSION(SIZE(PPS)) :: ZLETR
REAL, DIMENSION(SIZE(PPS)) :: ZEVAP
REAL, DIMENSION(SIZE(PPS)) :: ZSUBL !Sublimation
REAL, DIMENSION(SIZE(PPS)) :: ZGFLUX
REAL, DIMENSION(SIZE(PPS)) :: ZRESTORE
REAL, DIMENSION(SIZE(PPS)) :: ZUSTAR
REAL, DIMENSION(SIZE(PPS)) :: ZMELT
REAL, DIMENSION(SIZE(PPS),GDM%TV%R%CUR%TSNOW%NLAYER) :: ZSNOWTEMP
REAL, DIMENSION(SIZE(PPS),GDM%TV%R%CUR%TSNOW%NLAYER) :: ZSNOWLIQ
REAL, DIMENSION(SIZE(PPS),GDM%TV%R%CUR%TSNOW%NLAYER) :: ZSNOWDZ
REAL, DIMENSION(SIZE(PPS)) :: ZSNOWHMASS
REAL, DIMENSION(SIZE(PPS)) :: ZMELTADV
REAL, DIMENSION(SIZE(PPS),SIZE(GDM%TV%IP%XABC)) :: ZIACAN
REAL, DIMENSION(SIZE(PPS)) :: ZQS
REAL, DIMENSION(SIZE(PPS)) :: ZHV
REAL, DIMENSION(SIZE(PPS)) :: ZHORT
REAL, DIMENSION(SIZE(PPS)) :: ZDRIP
REAL, DIMENSION(SIZE(PPS)) :: ZTS
REAL, DIMENSION(SIZE(PPS)) :: ZRRVEG
REAL, DIMENSION(SIZE(PPS)) :: ZALBT
REAL, DIMENSION(SIZE(PPS)) :: ZEMIST
REAL, DIMENSION(SIZE(PPS)) :: ZGPP
REAL, DIMENSION(SIZE(PPS)) :: ZRESP_AUTO
REAL, DIMENSION(SIZE(PPS)) :: ZRESP_ECO
REAL, DIMENSION(SIZE(PPS)) :: ZFAPAR
REAL, DIMENSION(SIZE(PPS)) :: ZFAPIR
REAL, DIMENSION(SIZE(PPS)) :: ZFAPARC
REAL, DIMENSION(SIZE(PPS)) :: ZFAPIRC
REAL, DIMENSION(SIZE(PPS)) :: ZLAI_EFFC
REAL, DIMENSION(SIZE(PPS)) :: ZFAPAR_BS
REAL, DIMENSION(SIZE(PPS)) :: ZFAPIR_BS
REAL, DIMENSION(SIZE(PPS)) :: ZDFAPARC
REAL, DIMENSION(SIZE(PPS)) :: ZDFAPIRC
REAL, DIMENSION(SIZE(PPS)) :: ZDLAI_EFFC
REAL, DIMENSION(SIZE(PPS)) :: ZMUS
REAL, DIMENSION(SIZE(PPS)) :: ZIRRIG_FLUX
REAL, DIMENSION(0,0,0)     :: ZLITTER
REAL, DIMENSION(0,0)       :: ZLIGNIN_STRUC , ZSOILCARB, ZTURNOVER
!
REAL, DIMENSION(SIZE(PPS)) :: ZSNDRIFT
!
!  surfaces relative fractions
!  for flood
REAL, DIMENSION(SIZE(PPS)) :: ZFFG
REAL, DIMENSION(SIZE(PPS)) :: ZFFV
REAL, DIMENSION(SIZE(PPS)) :: ZFF
REAL, DIMENSION(SIZE(PPS)) :: ZALBF
REAL, DIMENSION(SIZE(PPS)) :: ZEMISF
REAL, DIMENSION(SIZE(PPS)) :: ZFFROZEN
REAL, DIMENSION(SIZE(PPS)) :: ZFFLOOD
REAL, DIMENSION(SIZE(PPS)) :: ZPIFLOOD
REAL, DIMENSION(SIZE(PPS)) :: ZIFLOOD
REAL, DIMENSION(SIZE(PPS)) :: ZPFLOOD
REAL, DIMENSION(SIZE(PPS)) :: ZLEFLOOD
REAL, DIMENSION(SIZE(PPS)) :: ZLEIFLOOD
REAL, DIMENSION(SIZE(PPS)) :: ZFFG_NOSNOW
REAL, DIMENSION(SIZE(PPS)) :: ZFFV_NOSNOW
!
!  variables for irrigation
REAL, DIMENSION(SIZE(PPS),1) :: ZIRRIG
REAL, DIMENSION(SIZE(PPS),1) :: ZWATSUP
REAL, DIMENSION(SIZE(PPS)) :: ZTHRESHOLDSPT
LOGICAL, DIMENSION(SIZE(PPS)) :: GIRRIGATE
LOGICAL, DIMENSION(SIZE(PPS)) :: GIRRIDAY
!
!
!  variables for deep soil temperature
REAL, DIMENSION(SIZE(PPS)) :: ZTDEEP_A
REAL, DIMENSION(SIZE(PPS)) :: ZTDEEP_B
REAL, DIMENSION(SIZE(PPS)) :: ZGAMMAT
!
REAL, DIMENSION(0) :: ZAOSIP  ! A/S for increasing x
REAL, DIMENSION(0) :: ZAOSIM  ! A/S for decreasing x
REAL, DIMENSION(0) :: ZAOSJP  ! A/S for increasing y
REAL, DIMENSION(0) :: ZAOSJM  ! A/S for decreasing y
REAL, DIMENSION(0) :: ZHO2IP  ! h/2 for increasing x
REAL, DIMENSION(0) :: ZHO2IM  ! h/2 for decreasing x
REAL, DIMENSION(0) :: ZHO2JP  ! h/2 for increasing y
REAL, DIMENSION(0) :: ZHO2JM  ! h/2 for decreasing y
REAL, DIMENSION(0,1) :: ZZ0EFFIP! roughness length for increasing x
REAL, DIMENSION(0,1) :: ZZ0EFFIM! roughness length for decreasing x
REAL, DIMENSION(0,1) :: ZZ0EFFJP! roughness length for increasing y
REAL, DIMENSION(0,1) :: ZZ0EFFJM! roughness length for decreasing y
REAL, DIMENSION(0) :: ZTAU_WOOD  ! residence time in wood (s)
REAL, DIMENSION(0,0) :: ZINCREASE
!
! Dummy variables for MEB:
LOGICAL,PARAMETER::OMEB=.FALSE.
LOGICAL,PARAMETER::OFORC_MEASURE=.FALSE.
LOGICAL,PARAMETER::OMEB_LITTER=.FALSE.
REAL, DIMENSION(SIZE(PPS)) :: ZP_MEB_SCA_SW,                     &
          ZP_RGLV, ZP_GAMMAV, ZP_RSMINV,                         &
          ZP_WRMAX_CFV, ZP_LAIV,                                 &
          ZP_BSLAI,ZP_LAIMIN,ZP_H_VEG,ZPALPHAN,                  &
          ZZ0G_WITHOUT_SNOW,                                     &
          ZZ0_MEBV,ZZ0H_MEBV,ZZ0EFF_MEBV,                        &
          ZZ0_MEBN,ZZ0H_MEBN,ZZ0EFF_MEBN,                        &
          ZP_ALBNIR_VEG, ZP_ALBVIS_VEG,                          &
          ZP_ALBNIR_SOIL, ZP_ALBVIS_SOIL, ZP_GNDLITTER
REAL, DIMENSION(SIZE(GDM%TV%M%X%XROOTFRAC,1),SIZE(GDM%TV%M%X%XROOTFRAC,2)) :: ZP_ROOTFRACV
REAL, DIMENSION(SIZE(PPS)) :: ZP_WRL,ZP_WRLI,ZP_WRVN,ZP_TV, ZP_TL
REAL, DIMENSION(SIZE(PPS)) :: ZP_TC,ZP_QC
REAL, DIMENSION(SIZE(PPS)) :: ZP_SWNET_V, ZP_SWNET_G, ZP_SWNET_N, ZP_SWNET_NS,    &
          ZP_LWNET_V, ZP_LWNET_G, ZP_LWNET_N,                                     &
          ZP_LEVCV, ZP_LESC, ZP_H_V_C, ZP_H_G_C,                                  &
          ZP_LETRGV, ZP_LETRCV, ZP_LERGV, ZP_LELITTER, ZP_LELITTERI,              &
          ZP_DRIPLIT,ZP_RRLIT, ZP_LERCV, ZP_H_C_A, ZP_H_N_C,                      &
          ZP_LE_C_A, ZP_LE_V_C, ZP_LE_G_C, ZP_LE_N_C,                             &
          ZP_EVAP_N_C, ZP_EVAP_G_C,                                               & 
          ZP_SR_GN, ZP_MELTCV, ZP_FRZCV,                                          &
          ZP_SWDOWN_GN, ZP_LWDOWN_GN
!
TYPE (DATE_TIME),   DIMENSION(0) :: TPSEED ! seeding date
TYPE (DATE_TIME),   DIMENSION(0) :: TPREAP ! reaping date
!
INTEGER                    :: ILU
!
LOGICAL :: GMASK, GAGRI_TO_GRASS
!
REAL, DIMENSION(0,0) :: ZGNDLITTER
REAL, DIMENSION(0,0) :: ZRGLGV
REAL, DIMENSION(0,0) :: ZGAMMAGV
REAL, DIMENSION(0,0) :: ZRSMINGV
REAL, DIMENSION(0,0) :: ZWRMAX_CFGV
REAL, DIMENSION(0,0) :: ZH_VEG
REAL, DIMENSION(0,0) :: ZLAIGV
REAL, DIMENSION(0,0) :: ZZ0LITTER
!
TYPE (DATE_TIME),  DIMENSION(0,0) :: TZSEED
TYPE (DATE_TIME), DIMENSION(0,0) :: TZREAP
!
!Snow options
LOGICAL :: GSNOWDRIFT,GSNOWDRIFT_SUBLIM,GSNOW_ABS_ZENITH
CHARACTER(3) :: YSNOWMETAMO,YSNOWRAD
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
HRAIN     = 'DEF'
OFLOOD    = .FALSE.
OTEMP_ARP = .FALSE.
OGLACIER  = .FALSE.
ZMUF      = 0.
ZFSAT     = 0.
ZTOPQS    = 0.
ZQSB      = 0.
ZFWTD     = 0.
ZWTD      = XUNDEF
ZSNDRIFT  = 0.
!
GAGRI_TO_GRASS = .FALSE.
!
!*      1.2    flood
!              -----
!
ZFFG          = 0.
ZFFV          = 0.
ZFF           = 0.
ZFFROZEN      = 0.
ZALBF         = 0.
ZEMISF        = 0.
ZIFLOOD       = 0.
ZPFLOOD       = 0.
ZFFLOOD       = 0.
ZPIFLOOD      = 0.
ZLEFLOOD      = 0.
ZLEIFLOOD     = 0.
ZFFG_NOSNOW   = 0.
ZFFV_NOSNOW   = 0.
!
!* ISBA like irrigation (not implemented)
!
ZIRRIG        = 0.
ZWATSUP       = 0.
ZTHRESHOLDSPT = 0.
GIRRIGATE     = .FALSE.
GIRRIDAY      = .FALSE.
!
ZTDEEP_A = XUNDEF
ZTDEEP_B = XUNDEF
ZGAMMAT  = XUNDEF
!
!-------------------------------------------------------------------------------
!
!* Variables required in TEB to allow coupling
!  with AROME/ALADIN/ARPEGE as LE or EVAP
!
ZLEI  = 0.0 ! sublimation heat flux (W/m2)
ZSUBL = 0.0 ! sublimation (kg/m2/s)
ZTS   = 0.0 ! surface temperature (K) (non-radiative)
!
!-------------------------------------------------------------------------------
! Snow options
GSNOWDRIFT=.TRUE.
GSNOWDRIFT_SUBLIM=.FALSE.
GSNOW_ABS_ZENITH=.FALSE.
YSNOWMETAMO="B92"
YSNOWRAD="B92"
!-------------------------------------------------------------------------------
!
!*      2.     Treatment of green areas
!              ------------------------
!
!
!
!*      2.1    Automatic irrigation
!              --------------------
!
CALL TEB_IRRIG(GDM%TIR%LPAR_GD_IRRIG, PTSTEP, TPTIME%TDATE%MONTH, PTSUN, &
               GDM%TIR%XGD_START_MONTH, GDM%TIR%XGD_END_MONTH, GDM%TIR%XGD_START_HOUR,   &
               GDM%TIR%XGD_END_HOUR, GDM%TIR%XGD_24H_IRRIG, PIRRIG_GARDEN        ) 
!
! --------------------------------------------------------------------------------------
! Vegetation update (in case of non-interactive vegetation):
! --------------------------------------------------------------------------------------
!
GUPDATED=.FALSE.
IF (GDM%TV%O%CPHOTO=='NON' .OR. GDM%TV%O%CPHOTO=='AGS' .OR. GDM%TV%O%CPHOTO=='AST') THEN
     CALL VEGETATION_UPDATE(DTCO, GDM%DTI, TG, GDM%TV%O,  &
                            PTSTEP,TPTIME,TOP%XCOVER,TOP%LCOVER,                 &
                         GDM%TV%O%CISBA,(.NOT. GDM%TV%O%LPAR), &
                         GDM%TV%O%CPHOTO, .FALSE.,     &
                         GDM%TV%O%LTR_ML, 'GRD',                                  &
                         GDM%TV%M%T%CUR%XLAI,GDM%TV%M%T%CUR%XVEG,GDM%TV%M%T%CUR%XZ0,         &
                         GDM%TV%M%T%CUR%XALBNIR,GDM%TV%M%T%CUR%XALBVIS,GDM%TV%M%T%CUR%XALBUV,&
                         GDM%TV%M%T%CUR%XEMIS,                   &
                         GDM%TV%M%T%CUR%XRSMIN,GDM%TV%M%T%CUR%XGAMMA,GDM%TV%M%T%CUR%XWRMAX_CF,   &
                         GDM%TV%M%T%CUR%XRGL,GDM%TV%M%T%CUR%XCV,                          &
                         GDM%TV%M%T%CUR%XGMES,GDM%TV%M%T%CUR%XBSLAI,GDM%TV%M%T%CUR%XLAIMIN,&
                         GDM%TV%M%T%CUR%XSEFOLD,GDM%TV%M%T%CUR%XGC,         &
                         GDM%TV%M%T%CUR%XF2I, GDM%TV%M%T%CUR%LSTRESS,                         &
                         ZAOSIP,ZAOSIM,ZAOSJP,ZAOSJM,                    &
                         ZHO2IP,ZHO2IM,ZHO2JP,ZHO2JM,                    &
                         ZZ0EFFIP,ZZ0EFFIM,ZZ0EFFJP,ZZ0EFFJM,            &
                         GDM%TV%O%CALBEDO, GDM%TV%M%T%CUR%XALBNIR_VEG, &
                         GDM%TV%M%T%CUR%XALBVIS_VEG, GDM%TV%M%T%CUR%XALBUV_VEG,  &
                         GDM%TV%M%A%XALBNIR_SOIL, GDM%TV%M%A%XALBVIS_SOIL, GDM%TV%M%A%XALBUV_SOIL,    &
                         GDM%TV%M%T%CUR%XCE_NITRO, GDM%TV%M%T%CUR%XCF_NITRO, GDM%TV%M%T%CUR%XCNA_NITRO,  &
                         TZSEED, TZREAP, ZWATSUP, ZIRRIG,                &
                         ZGNDLITTER,ZRGLGV,ZGAMMAGV,                     &
                         ZRSMINGV, ZWRMAX_CFGV,                          &
                         ZH_VEG, ZLAIGV, ZZ0LITTER,                      &
                         GUPDATED, OABSENT=(T%CUR%XGARDEN==0.)                 )
END IF
!
!*      2.2    Call ISBA for green areas
!              -------------------------
!
!
 CALL ISBA(GDM%TV%O%CISBA, GDM%TV%O%CPHOTO, GDM%TV%O%LTR_ML, GDM%TV%O%CRUNOFF, &
           GDM%TV%O%CKSAT, HRAIN, GDM%TV%O%CHORT,       &
          GDM%TV%O%CC1DRY, GDM%TV%O%CSCOND, GDM%TV%R%CUR%TSNOW%SCHEME, &
          GDM%TV%O%CSNOWRES, GDM%TV%O%CCPSURF, GDM%TV%O%CSOILFRZ,  &
          GDM%TV%O%CDIFSFCOND, TPTIME, OFLOOD, OTEMP_ARP, OGLACIER,        &
          OMEB, OFORC_MEASURE,OMEB_LITTER, PTSTEP, HIMPLICIT_WIND, GAGRI_TO_GRASS,&
          GSNOWDRIFT,GSNOWDRIFT_SUBLIM,GSNOW_ABS_ZENITH,              &
          YSNOWMETAMO,YSNOWRAD,                                       &
          GDM%TV%O%XCGMAX, PZ_LOWCAN, PZ_LOWCAN, ZDIRCOSZW, PT_LOWCAN,         &
          PQ_LOWCAN, PEXNS, PRHOA, PPS, PEXNS,  PRR, PSR, PZENITH,    &
          ZP_MEB_SCA_SW,                                              &
          PSW, PLW, PU_LOWCAN, PPEW_A_COEF, PPEW_B_COEF, PPET_A_COEF, &
          PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF, GDM%TV%M%T%CUR%XRSMIN(:,1), &
          GDM%TV%M%T%CUR%XRGL(:,1), GDM%TV%M%T%CUR%XGAMMA(:,1),&
          GDM%TV%M%T%CUR%XCV(:,1), GDM%TV%IP%XRUNOFFD(:,1), GDM%TV%IP%XSOILWGHT(:,:,1), &
          GDM%TV%O%NLAYER_HORT, GDM%TV%O%NLAYER_DUN,          &
          PALBNIR_TVEG, PALBVIS_TVEG, PALBNIR_TSOIL, PALBVIS_TSOIL,   &
          GDM%TV%R%CUR%XSNOWFREE_ALB(:,1), GDM%TV%M%T%CUR%XWRMAX_CF(:,1), GDM%TV%M%T%CUR%XVEG(:,1), &
          GDM%TV%M%T%CUR%XLAI(:,1), GDM%TV%M%T%CUR%XEMIS(:,1), GDM%TV%M%T%CUR%XZ0(:,1),       &
          GDM%TV%M%T%CUR%XZ0(:,1)/GDM%TV%M%X%XZ0_O_Z0H(:,1), GDM%TV%M%X%XVEGTYPE, GDM%TV%M%T%CUR%XZ0(:,1),  &
          ZP_RGLV, ZP_GAMMAV, ZP_RSMINV,                               &
          ZP_ROOTFRACV, ZP_WRMAX_CFV, ZP_LAIV,                        &
          ZP_BSLAI,ZP_LAIMIN,ZP_H_VEG,ZPALPHAN,                       &
          ZZ0G_WITHOUT_SNOW,                                          &
          ZZ0_MEBV,ZZ0H_MEBV,ZZ0EFF_MEBV,                             &
          ZZ0_MEBN,ZZ0H_MEBN,ZZ0EFF_MEBN, ZP_GNDLITTER,               &
          GDM%TV%P%XRUNOFFB, GDM%TV%IP%XCGSAT, GDM%TV%IP%XC1SAT(:,1),      &
          GDM%TV%IP%XC2REF(:,1), GDM%TV%IP%XC3(:,:,1), GDM%TV%IP%XC4B(:), GDM%TV%IP%XC4REF(:,1), &
          GDM%TV%IP%XACOEF(:), GDM%TV%IP%XPCOEF(:), GDM%TV%IP%XTAUICE, GDM%TV%P%XWDRAIN,&
          ZTDEEP_A, ZTDEEP_B, ZGAMMAT, GDM%TV%R%CUR%XPSN(:,1), GDM%TV%R%CUR%XPSNG(:,1), &
          GDM%TV%R%CUR%XPSNV(:,1), GDM%TV%R%CUR%XPSNV_A(:,1),   &
          GDM%TV%R%CUR%XSNOWFREE_ALB_VEG(:,1), GDM%TV%R%CUR%XSNOWFREE_ALB_SOIL(:,1), ZIRRIG(:,1), ZWATSUP(:,1),     &
          ZTHRESHOLDSPT, GIRRIGATE, GIRRIDAY, GDM%TV%M%T%CUR%LSTRESS(:,1), GDM%TV%M%T%CUR%XGC(:,1), GDM%TV%M%T%CUR%XF2I(:,1),     &
          GDM%TV%M%X%XDMAX(:,1), GDM%TV%IP%XAH(:,1), GDM%TV%IP%XBH(:,1), PCO2, GDM%TV%M%T%CUR%XGMES(:,1), GDM%TV%IP%XPOI, &
          GDM%TV%IP%XFZERO(:,1), GDM%TV%IP%XEPSO(:,1), GDM%TV%IP%XGAMM(:,1),   &
          GDM%TV%IP%XQDGAMM(:,1), GDM%TV%IP%XQDGMES(:,1), GDM%TV%IP%XT1GMES(:,1), GDM%TV%IP%XT2GMES(:,1), &
          GDM%TV%IP%XAMAX(:,1), GDM%TV%IP%XQDAMAX(:,1), GDM%TV%IP%XT1AMAX(:,1),&
          GDM%TV%IP%XT2AMAX(:,1), GDM%TV%IP%XABC, GDM%TV%M%X%XDG(:,:,1), GDM%TV%IP%XDZG(:,:,1), GDM%TV%IP%XDZDIF(:,:,1), &
          GDM%TV%M%X%NWG_LAYER(:,1), GDM%TV%M%X%XROOTFRAC(:,:,1),     &
          GDM%TV%IP%XWFC, GDM%TV%IP%XWWILT, GDM%TV%IP%XWSAT, GDM%TV%IP%XBCOEF,  &
          GDM%TV%IP%XCONDSAT(:,:,1), GDM%TV%IP%XMPOTSAT(:,:),           &
          GDM%TV%IP%XHCAPSOIL(:,:), GDM%TV%IP%XCONDDRY(:,:), GDM%TV%IP%XCONDSLD(:,:), GDM%TV%M%X%XD_ICE(:,1), &
          GDM%TV%IP%XKSAT_ICE(:,1), ZMUF, ZFF,&
          ZFFG, ZFFV, ZFFG_NOSNOW, ZFFV_NOSNOW, ZFFROZEN,  ZALBF,     &
          ZEMISF, ZFFLOOD, ZPIFLOOD, ZIFLOOD, ZPFLOOD, ZLEFLOOD,      &
          ZLEIFLOOD, ZSODELX, TG%XLAT, TG%XLON, GDM%TV%R%CUR%XTG(:,:,1), GDM%TV%R%CUR%XWG(:,:,1), &
          GDM%TV%R%CUR%XWGI(:,:,1), GDM%TV%IP%XPCPS(:,1),      &
          GDM%TV%IP%XPLVTT(:,1), GDM%TV%IP%XPLSTT(:,1), GDM%TV%R%CUR%XWR(:,1),                             &
          ZP_WRL,ZP_WRLI,ZP_WRVN,ZP_TV, ZP_TL,                                &
          GDM%TV%R%CUR%XRESA(:,1), GDM%TV%R%CUR%XANFM(:,1), ZFSAT, GDM%TV%R%CUR%TSNOW%ALB(:,1),          &
          GDM%TV%R%CUR%TSNOW%ALBVIS(:,1), GDM%TV%R%CUR%TSNOW%ALBNIR(:,1), GDM%TV%R%CUR%TSNOW%ALBFIR(:,1),    &
          GDM%TV%R%CUR%TSNOW%WSNOW(:,:,1), GDM%TV%R%CUR%TSNOW%HEAT(:,:,1), GDM%TV%R%CUR%TSNOW%RHO(:,:,1),    &
          GDM%TV%R%CUR%TSNOW%GRAN1(:,:,1), GDM%TV%R%CUR%TSNOW%GRAN2(:,:,1), GDM%TV%R%CUR%TSNOW%HIST(:,:,1),  &
          GDM%TV%R%CUR%TSNOW%AGE(:,:,1), ZGRNDFLUX, ZHPSNOW, ZSNOWHMASS,           &
          ZRNSNOW, ZHSNOW,  ZGFLUXSNOW, ZUSTARSNOW,                   &
          ZSRSFC, ZRRSFC, ZLESL, GDM%TV%R%CUR%TSNOW%EMIS(:,1), ZCDSNOW, ZCHSNOW,   &
          PTS_GARDEN, ZTS, ZHV, ZQS, ZSNOWTEMP, ZSNOWLIQ, ZSNOWDZ,    &
          ZCG, ZC1, ZC2, ZWGEQ, ZCT, ZCH, ZCD, ZCDN, ZRI, ZHU, ZHUG,  &
          ZEMIST, ZALBT, ZRS, GDM%TV%R%CUR%XLE(:,1),  ZRN, ZH, ZLEI, ZLEGI, ZLEG, ZLEV, &
          ZLES, ZLER, ZLETR, ZEVAP, ZGFLUX, ZRESTORE, ZUSTAR,         &
          PDRAIN_GARDEN, PRUNOFF_GARDEN,                              &
          ZMELT, ZMELTADV,                                            &
          ZP_TC,ZP_QC,                                                &
          ZRN_ISBA, ZH_ISBA, ZLEG_ISBA,                               &
          ZLEGI_ISBA, ZLEV_ISBA, ZLETR_ISBA, ZUSTAR_ISBA, ZLER_ISBA,  &
          ZLE_ISBA, ZLEI_ISBA, ZGFLUX_ISBA, ZHORT, ZDRIP, ZRRVEG,     &
          PAC_AGG_GARDEN, PHU_AGG_GARDEN, ZFAPARC, ZFAPIRC, ZMUS,     &
          ZLAI_EFFC, GDM%TV%R%CUR%XAN(:,1), GDM%TV%R%CUR%XANDAY(:,1), ZRESP_BIOMASS_INST, ZIACAN,  &
          ZGPP, ZFAPAR, ZFAPIR, ZFAPAR_BS, ZFAPIR_BS, ZIRRIG_FLUX,    &
          ZDEEP_FLUX,                                                 &
          ZP_SWNET_V, ZP_SWNET_G, ZP_SWNET_N, ZP_SWNET_NS,            &
          ZP_LWNET_V, ZP_LWNET_G, ZP_LWNET_N,                         &
          ZP_LEVCV, ZP_LESC, ZP_H_V_C, ZP_H_G_C,                      &
          ZP_LETRGV, ZP_LETRCV, ZP_LERGV, ZP_LELITTER,ZP_LELITTERI,   &
          ZP_DRIPLIT,ZP_RRLIT, ZP_LERCV, ZP_H_C_A, ZP_H_N_C,   &
          ZP_LE_C_A, ZP_LE_V_C, ZP_LE_G_C, ZP_LE_N_C,                 &
          ZP_EVAP_N_C, ZP_EVAP_G_C,                                   & 
          ZP_SR_GN, ZP_MELTCV, ZP_FRZCV,                              &
          ZP_SWDOWN_GN, ZP_LWDOWN_GN,                                 &
          PIRRIG_GARDEN, ZTOPQS, ZQSB, ZSUBL, ZFWTD, ZWTD, ZSNDRIFT   )                                                           
!
!
IF (GDM%TV%R%CUR%TSNOW%SCHEME=='3-L' .OR. GDM%TV%R%CUR%TSNOW%SCHEME=='CRO') &
        GDM%TV%R%CUR%TSNOW%TS(:,1)=ZSNOWTEMP(:,1)
!
IF (GDM%TV%O%LTR_ML) THEN
  GMASK = ( TPTIME%TIME - PTSTEP < 0. ) .AND. ( TPTIME%TIME >= 0. )
  IF (GMASK) THEN
    ZDFAPARC  (:) = ZFAPARC   (:) / ZMUS (:)
    ZDFAPIRC  (:) = ZFAPIRC   (:) / ZMUS (:)
    ZDLAI_EFFC(:) = ZLAI_EFFC (:) / ZMUS (:)
  ENDIF
ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Vegetation evolution for interactive LAI
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
IF (GDM%TV%O%CPHOTO=='LAI' .OR. GDM%TV%O%CPHOTO=='LST' .OR. GDM%TV%O%CPHOTO=='NIT') THEN
  CALL VEGETATION_EVOL(GDM%TV%O%CISBA, GDM%TV%O%CPHOTO, GDM%TV%O%CRESPSL, &
                       GDM%TV%O%CALBEDO, .FALSE., GDM%TV%O%LTR_ML,   &
                       GDM%TV%O%LNITRO_DILU, GAGRI_TO_GRASS,                        &
                       PTSTEP, TPTIME%TDATE%MONTH, TPTIME%TDATE%DAY, 1,    &
                       TPTIME%TIME, TG%XLAT, PRHOA, GDM%TV%M%X%XDG(:,:,1), GDM%TV%IP%XDZG(:,:,1), &
                       GDM%TV%M%X%NWG_LAYER(:,1),  & 
                       GDM%TV%R%CUR%XTG(:,:,1), GDM%TV%M%T%CUR%XALBNIR_VEG(:,1), GDM%TV%M%T%CUR%XALBVIS_VEG(:,1), &
                       GDM%TV%M%T%CUR%XALBUV_VEG(:,1), GDM%TV%M%A%XALBNIR_SOIL(:,1), GDM%TV%M%A%XALBVIS_SOIL(:,1), &
                       GDM%TV%M%A%XALBUV_SOIL(:,1), GDM%TV%M%X%XVEGTYPE, GDM%TV%M%T%CUR%XSEFOLD(:,1), GDM%TV%IP%XANMAX(:,1), &
                       GDM%TV%M%X%XH_TREE(:,1), GDM%TV%M%T%CUR%XBSLAI(:,1), GDM%TV%M%T%CUR%XLAIMIN(:,1), PCO2, &
                       GDM%TV%M%T%CUR%XCE_NITRO(:,1), GDM%TV%M%T%CUR%XCF_NITRO(:,1), GDM%TV%M%T%CUR%XCNA_NITRO(:,1),    &
                       GDM%TV%IP%XBSLAI_NITRO(:,1), GDM%TV%M%T%CUR%XGMES(:,1), ZTAU_WOOD, TPSEED,        &
                       TPREAP, ZAOSIP, ZAOSIM, ZAOSJP, ZAOSJM,             &
                       ZHO2IP, ZHO2IM, ZHO2JP, ZHO2JM, ZZ0EFFIP(:,1),           &
                       ZZ0EFFIM(:,1), ZZ0EFFJP(:,1), ZZ0EFFJM(:,1), GDM%TV%M%T%CUR%XLAI(:,1), GDM%TV%M%T%CUR%XVEG(:,1),   &
                       GDM%TV%M%T%CUR%XZ0(:,1), GDM%TV%M%T%CUR%XALBNIR(:,1), GDM%TV%M%T%CUR%XALBVIS(:,1), &
                       GDM%TV%M%T%CUR%XALBUV(:,1), GDM%TV%M%T%CUR%XEMIS(:,1), GDM%TV%R%CUR%XANFM(:,1), &
                       GDM%TV%R%CUR%XANDAY(:,1), GDM%TV%R%CUR%XBIOMASS(:,:,1), GDM%TV%R%CUR%XRESP_BIOMASS(:,:,1),        &
                       ZRESP_BIOMASS_INST, ZINCREASE, ZTURNOVER             )         
END IF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Diagnostic of respiration carbon fluxes and soil carbon evolution
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
PSFCO2(:)=0.
ZRESP_ECO (:)=0.
ZRESP_AUTO(:)=0.
!
IF (GDM%TV%O%CPHOTO/='NON' .AND. GDM%TV%O%CRESPSL/='NON' .AND. ANY(GDM%TV%M%T%CUR%XLAI(:,1)/=XUNDEF)) THEN
  CALL CARBON_EVOL(GDM%TV%O%CISBA, GDM%TV%O%CRESPSL, GDM%TV%O%CPHOTO, PTSTEP, 1,             &
                   PRHOA, GDM%TV%R%CUR%XTG(:,:,1), GDM%TV%R%CUR%XWG(:,:,1), GDM%TV%IP%XWFC(:,:), &
                   GDM%TV%IP%XWWILT(:,:), GDM%TV%IP%XWSAT(:,:), GDM%TV%P%XSAND,   &
                   GDM%TV%M%X%XDG(:,:,1), GDM%TV%IP%XDZG(:,:,1), GDM%TV%M%X%NWG_LAYER(:,1),           &                   
                   GDM%TV%M%X%XRE25(:,1), GDM%TV%M%T%CUR%XLAI(:,1), ZRESP_BIOMASS_INST, ZTURNOVER,    &
                   ZLITTER, ZLIGNIN_STRUC , ZSOILCARB,            &
                   ZRESP_AUTO, ZRESP_ECO                          )
  ! calculation of vegetation CO2 flux
  PSFCO2(:) = - ZGPP(:) + ZRESP_ECO(:)
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
 CALL FLAG_TEB_GARDEN_n(GDM%TV%R, GDM%TV%O, GDM%TV%M%T, T, 2)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!*      9.     Fields required for TEB
!              -----------------------
!
WHERE (T%CUR%XGARDEN/=0.)
  !
  ! energy balance
  !
   PRN_GARDEN    (:) = ZRN    (:)
   PH_GARDEN     (:) = ZH     (:)
   PLE_GARDEN    (:) = GDM%TV%R%CUR%XLE(:,1)
   PGFLUX_GARDEN (:) = ZGFLUX (:)
   PEVAP_GARDEN  (:) = ZEVAP  (:)
  !
  !
  ! Estimate of green area aerodynamic conductance recomputed from heat flux,
  ! surface (radiative) temp. and forcing air temperature (estimated at future time step)
  ZTA = PPET_B_COEF + PPET_A_COEF * PH_GARDEN
  PAC_GARDEN = 0.
  WHERE (PTS_GARDEN /= ZTA)
    PAC_GARDEN(:)   = MAX(PH_GARDEN(:) / XCPD / PRHOA(:) / (PTS_GARDEN - ZTA) , 0.)
  ENDWHERE
  !
  ! Humidity of saturation for green areas
  PQSAT_GARDEN(:) = QSAT(GDM%TV%R%CUR%XTG(:,1,1),PPS(:))
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
