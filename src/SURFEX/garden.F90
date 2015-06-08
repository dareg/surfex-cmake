!     #########
    SUBROUTINE GARDEN(HIMPLICIT_WIND, TPTIME, PTSUN, PPEW_A_COEF, PPEW_B_COEF,       &
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
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS
!
USE MODD_TYPE_DATE_SURF,    ONLY: DATE_TIME
USE MODD_SURF_PAR,          ONLY: XUNDEF
USE MODD_CSTS,              ONLY: XCPD
USE MODD_TEB_n, ONLY : T => TEB
USE MODD_TEB_GRID_n, ONLY : TG => TEB_GRID
USE MODD_TEB_VEG_n, ONLY : TVG => TEB_VEG_OPTIONS
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TGDO => TEB_GARDEN_OPTIONS
USE MODD_TEB_GARDEN_PGD_EVOL_n, ONLY : TGDPE => TEB_GARDEN_PGD_EVOL
USE MODD_TEB_GARDEN_PGD_n, ONLY : TGDP => TEB_GARDEN_PGD
USE MODD_TEB_GARDEN_n, ONLY : TGD => TEB_GARDEN
!
USE MODD_TEB_IRRIG_n, ONLY : TIR => TEB_IRRIG
!
USE MODI_ISBA
USE MODI_VEGETATION_UPDATE_GARDEN
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
REAL, DIMENSION(SIZE(PPS),TGDO%NGROUND_LAYER) :: ZTOPQS  ! Topmodel (SGH not used in TEB) lateral subsurface flow by layer
REAL, DIMENSION(SIZE(PPS))               :: ZQSB    ! Topmodel (SGH not used in TEB) output lateral subsurface
REAL, DIMENSION(SIZE(PPS))               :: ZFWTD   ! grid-cell fraction of water table to rise
REAL, DIMENSION(SIZE(PPS))               :: ZWTD    ! water table depth from Obs, TRIP or MODCOU
!
REAL, DIMENSION(SIZE(PPS)) :: ZDIRCOSZW           ! orography slope cosine (=1 in TEB)
REAL, DIMENSION(SIZE(PPS),TVG%NNBIOMASS) :: ZRESP_BIOMASS_INST       ! instantaneous biomass respiration (kgCO2/kgair m/s)
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
REAL, DIMENSION(SIZE(PPS),TGD%CUR%TSNOW%NLAYER) :: ZSNOWTEMP
REAL, DIMENSION(SIZE(PPS),TGD%CUR%TSNOW%NLAYER) :: ZSNOWLIQ
REAL, DIMENSION(SIZE(PPS),TGD%CUR%TSNOW%NLAYER) :: ZSNOWDZ
REAL, DIMENSION(SIZE(PPS)) :: ZSNOWHMASS
REAL, DIMENSION(SIZE(PPS)) :: ZMELTADV
REAL, DIMENSION(SIZE(PPS),SIZE(TGDP%XABC)) :: ZIACAN
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
REAL, DIMENSION(SIZE(PPS)) :: ZIRRIG
REAL, DIMENSION(SIZE(PPS)) :: ZWATSUP
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
REAL, DIMENSION(0) :: ZZ0EFFIP! roughness length for increasing x
REAL, DIMENSION(0) :: ZZ0EFFIM! roughness length for decreasing x
REAL, DIMENSION(0) :: ZZ0EFFJP! roughness length for increasing y
REAL, DIMENSION(0) :: ZZ0EFFJM! roughness length for decreasing y
REAL, DIMENSION(0) :: ZTAU_WOOD  ! residence time in wood (s)
REAL, DIMENSION(0,0) :: ZINCREASE
!
! Dummy variables for MEB:
LOGICAL,PARAMETER::OMEB=.FALSE.
LOGICAL,PARAMETER::OFORC_MEASURE=.FALSE.
REAL, DIMENSION(SIZE(PPS)) :: ZP_MEB_SCA_SW,                     &
          ZP_ZF_TALLVEG , ZP_RGLV, ZP_GAMMAV, ZP_RSMINV,         &
          ZP_WRMAX_CFV, ZP_LAIV,                                 &
          ZP_BSLAI,ZP_LAIMIN,ZP_H_VEG,ZPALPHAN,                  &
          ZZ0G_WITHOUT_SNOW,                                     &
          ZZ0_MEBV,ZZ0H_MEBV,ZZ0EFF_MEBV,                        &
          ZZ0_MEBN,ZZ0H_MEBN,ZZ0EFF_MEBN,                        &
          ZP_ALBNIR_VEG, ZP_ALBVIS_VEG,                          &
          ZP_ALBNIR_SOIL, ZP_ALBVIS_SOIL, ZP_GNDLITTER
REAL, DIMENSION(SIZE(TGDP%XROOTFRAC,1),SIZE(TGDP%XROOTFRAC,2)) :: ZP_ROOTFRACV
REAL, DIMENSION(SIZE(PPS)) :: ZP_WRV,ZP_WRVN,ZP_TV
REAL, DIMENSION(SIZE(PPS)) :: ZP_TC,ZP_QC
REAL, DIMENSION(SIZE(PPS)) :: ZP_SWNET_V, ZP_SWNET_G, ZP_SWNET_N, ZP_SWNET_NS,    &
          ZP_LWNET_V, ZP_LWNET_G, ZP_LWNET_N,                                     &
          ZP_LEVCV, ZP_LESC, ZP_H_V_C, ZP_H_G_C,                                  &
          ZP_LETRGV, ZP_LETRCV, ZP_LERGV, ZP_LERCV, ZP_H_C_A, ZP_H_N_C,           &
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
!Snow options
LOGICAL :: GSNOWDRIFT,GSNOWDRIFT_SUBLIM,GSNOW_ABS_ZENITH
CHARACTER(3) :: YSNOWMETAMO,YSNOWRAD
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
CALL TEB_IRRIG(TIR%LPAR_GD_IRRIG, PTSTEP, TPTIME%TDATE%MONTH, PTSUN, &
               TIR%XGD_START_MONTH, TIR%XGD_END_MONTH, TIR%XGD_START_HOUR,   &
               TIR%XGD_END_HOUR, TIR%XGD_24H_IRRIG, PIRRIG_GARDEN        ) 
!
!*      2.2    Call ISBA for green areas
!              -------------------------
!
!
 CALL ISBA(TVG%CISBA, TVG%CPHOTO, TVG%LTR_ML, TVG%CRUNOFF, TVG%CKSAT, HRAIN, TVG%CHORT,       &
          TVG%CC1DRY, TVG%CSCOND, TGD%CUR%TSNOW%SCHEME, TVG%CSNOWRES, TVG%CCPSURF, TVG%CSOILFRZ,  &
          TVG%CDIFSFCOND, TPTIME, OFLOOD, OTEMP_ARP, OGLACIER,            &
          OMEB, OFORC_MEASURE, PTSTEP, HIMPLICIT_WIND, GAGRI_TO_GRASS,&
          GSNOWDRIFT,GSNOWDRIFT_SUBLIM,GSNOW_ABS_ZENITH,              &
          YSNOWMETAMO,YSNOWRAD,                                       &
          TVG%XCGMAX, PZ_LOWCAN, PZ_LOWCAN, ZDIRCOSZW, PT_LOWCAN,         &
          PQ_LOWCAN, PEXNS, PRHOA, PPS, PEXNS,  PRR, PSR, PZENITH,    &
          ZP_MEB_SCA_SW,                                              &
          PSW, PLW, PU_LOWCAN, PPEW_A_COEF, PPEW_B_COEF, PPET_A_COEF, &
          PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF, TGDP%XRSMIN, TGDP%XRGL, TGDP%XGAMMA,&
          TGDP%XCV, TGDP%XRUNOFFD, TGDP%XSOILWGHT, TGDO%NLAYER_HORT, TGDO%NLAYER_DUN,          &
          PALBNIR_TVEG, PALBVIS_TVEG, PALBNIR_TSOIL, PALBVIS_TSOIL,   &
          TGD%CUR%XSNOWFREE_ALB, TGDP%XWRMAX_CF, TGDPE%CUR%XVEG, TGDPE%CUR%XLAI, TGDPE%CUR%XEMIS, TGDPE%CUR%XZ0,           &
          TGDPE%CUR%XZ0/TGDP%XZ0_O_Z0H, TGDP%XVEGTYPE, TGDPE%CUR%XZ0,                               &
          ZP_ZF_TALLVEG , ZP_RGLV, ZP_GAMMAV, ZP_RSMINV,              &
          ZP_ROOTFRACV, ZP_WRMAX_CFV, ZP_LAIV,                        &
          ZP_BSLAI,ZP_LAIMIN,ZP_H_VEG,ZPALPHAN,                       &
          ZZ0G_WITHOUT_SNOW,                                          &
          ZZ0_MEBV,ZZ0H_MEBV,ZZ0EFF_MEBV,                             &
          ZZ0_MEBN,ZZ0H_MEBN,ZZ0EFF_MEBN, ZP_GNDLITTER,               &
          TGDP%XRUNOFFB, TGDP%XCGSAT, TGDP%XC1SAT,                                   &
          TGDP%XC2REF, TGDP%XC3, TGDP%XC4B, TGDP%XC4REF, TGDP%XACOEF, TGDP%XPCOEF, TGDP%XTAUICE, TGDP%XWDRAIN,&
          ZTDEEP_A, ZTDEEP_B, ZGAMMAT, TGD%CUR%XPSN, TGD%CUR%XPSNG, TGD%CUR%XPSNV, TGD%CUR%XPSNV_A,   &
          TGD%CUR%XSNOWFREE_ALB_VEG, TGD%CUR%XSNOWFREE_ALB_SOIL, ZIRRIG, ZWATSUP,     &
          ZTHRESHOLDSPT, GIRRIGATE, GIRRIDAY, TGDP%LSTRESS, TGDP%XGC, TGDP%XF2I,     &
          TGDP%XDMAX, TGDP%XAH, TGDP%XBH, PCO2, TGDP%XGMES, TGDP%XPOI, TGDP%XFZERO, TGDP%XEPSO, TGDP%XGAMM,   &
          TGDP%XQDGAMM, TGDP%XQDGMES, TGDP%XT1GMES, TGDP%XT2GMES, TGDP%XAMAX, TGDP%XQDAMAX, TGDP%XT1AMAX,&
          TGDP%XT2AMAX, TGDP%XABC, TGDP%XDG, TGDP%XDZG, TGDP%XDZDIF, TGDP%NWG_LAYER, TGDP%XROOTFRAC,     &
          TGDP%XWFC, TGDP%XWWILT, TGDP%XWSAT, TGDP%XBCOEF,  TGDP%XCONDSAT, TGDP%XMPOTSAT,           &
          TGDP%XHCAPSOIL, TGDP%XCONDDRY, TGDP%XCONDSLD, TGDP%XD_ICE, TGDP%XKSAT_ICE, ZMUF, ZFF,&
          ZFFG, ZFFV, ZFFG_NOSNOW, ZFFV_NOSNOW, ZFFROZEN,  ZALBF,     &
          ZEMISF, ZFFLOOD, ZPIFLOOD, ZIFLOOD, ZPFLOOD, ZLEFLOOD,      &
          ZLEIFLOOD, ZSODELX, TG%XLAT, TG%XLON, TGD%CUR%XTG, TGD%CUR%XWG, TGD%CUR%XWGI, TGDP%XPCPS,      &
          TGDP%XPLVTT, TGDP%XPLSTT, TGD%CUR%XWR,                                        &
          ZP_WRV,ZP_WRVN,ZP_TV,                                       &
          TGD%CUR%XRESA, TGD%CUR%XANFM, ZFSAT, TGD%CUR%TSNOW%ALB(:,1),                        &
          TGD%CUR%TSNOW%ALBVIS(:,1), TGD%CUR%TSNOW%ALBNIR(:,1), TGD%CUR%TSNOW%ALBFIR(:,1),    &
          TGD%CUR%TSNOW%WSNOW(:,:,1), TGD%CUR%TSNOW%HEAT(:,:,1), TGD%CUR%TSNOW%RHO(:,:,1),    &
          TGD%CUR%TSNOW%GRAN1(:,:,1), TGD%CUR%TSNOW%GRAN2(:,:,1), TGD%CUR%TSNOW%HIST(:,:,1),  &
          TGD%CUR%TSNOW%AGE(:,:,1), ZGRNDFLUX, ZHPSNOW, ZSNOWHMASS,           &
          ZRNSNOW, ZHSNOW,  ZGFLUXSNOW, ZUSTARSNOW,                   &
          ZSRSFC, ZRRSFC, ZLESL, TGD%CUR%TSNOW%EMIS(:,1), ZCDSNOW, ZCHSNOW,   &
          PTS_GARDEN, ZTS, ZHV, ZQS, ZSNOWTEMP, ZSNOWLIQ, ZSNOWDZ,    &
          ZCG, ZC1, ZC2, ZWGEQ, ZCT, ZCH, ZCD, ZCDN, ZRI, ZHU, ZHUG,  &
          ZEMIST, ZALBT, ZRS, TGD%CUR%XLE,  ZRN, ZH, ZLEI, ZLEGI, ZLEG, ZLEV, &
          ZLES, ZLER, ZLETR, ZEVAP, ZGFLUX, ZRESTORE, ZUSTAR,         &
          PDRAIN_GARDEN, PRUNOFF_GARDEN,                              &
          ZMELT, ZMELTADV,                                            &
          ZP_TC,ZP_QC,                                                &
          ZRN_ISBA, ZH_ISBA, ZLEG_ISBA,                               &
          ZLEGI_ISBA, ZLEV_ISBA, ZLETR_ISBA, ZUSTAR_ISBA, ZLER_ISBA,  &
          ZLE_ISBA, ZLEI_ISBA, ZGFLUX_ISBA, ZHORT, ZDRIP, ZRRVEG,     &
          PAC_AGG_GARDEN, PHU_AGG_GARDEN, ZFAPARC, ZFAPIRC, ZMUS,     &
          ZLAI_EFFC, TGD%CUR%XAN, TGD%CUR%XANDAY, ZRESP_BIOMASS_INST, ZIACAN, TGDP%XANF,   &
          ZGPP, ZFAPAR, ZFAPIR, ZFAPAR_BS, ZFAPIR_BS, ZIRRIG_FLUX,    &
          ZDEEP_FLUX,                                                 &
          ZP_SWNET_V, ZP_SWNET_G, ZP_SWNET_N, ZP_SWNET_NS,            &
          ZP_LWNET_V, ZP_LWNET_G, ZP_LWNET_N,                         &
          ZP_LEVCV, ZP_LESC, ZP_H_V_C, ZP_H_G_C,                      &
          ZP_LETRGV, ZP_LETRCV, ZP_LERGV, ZP_LERCV, ZP_H_C_A, ZP_H_N_C,   &
          ZP_LE_C_A, ZP_LE_V_C, ZP_LE_G_C, ZP_LE_N_C,                 &
          ZP_EVAP_N_C, ZP_EVAP_G_C,                                   & 
          ZP_SR_GN, ZP_MELTCV, ZP_FRZCV,                              &
          ZP_SWDOWN_GN, ZP_LWDOWN_GN,                                 &
          PIRRIG_GARDEN, ZTOPQS, ZQSB, ZSUBL, ZFWTD, ZWTD, ZSNDRIFT   )                                                           
!
!
IF (TGD%CUR%TSNOW%SCHEME=='3-L' .OR. TGD%CUR%TSNOW%SCHEME=='CRO') TGD%CUR%TSNOW%TS(:,1)=ZSNOWTEMP(:,1)
!
IF (TVG%LTR_ML) THEN
  GMASK = ( TPTIME%TIME - PTSTEP < 0. ) .AND. ( TPTIME%TIME >= 0. )
  IF (GMASK) THEN
    ZDFAPARC  (:) = ZFAPARC   (:) / ZMUS (:)
    ZDFAPIRC  (:) = ZFAPIRC   (:) / ZMUS (:)
    ZDLAI_EFFC(:) = ZLAI_EFFC (:) / ZMUS (:)
  ENDIF
ENDIF
!
! --------------------------------------------------------------------------------------
! Vegetation update (in case of non-interactive vegetation):
! --------------------------------------------------------------------------------------
!
IF (TVG%CPHOTO=='NON' .OR. TVG%CPHOTO=='AGS' .OR. TVG%CPHOTO=='AST') THEN
     CALL VEGETATION_UPDATE_GARDEN(TGDO, TGDPE, TGDP, T, TOP, TVG, &
                                   TPTIME,PTSTEP,ILU)  
END IF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Vegetation evolution for interactive LAI
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
IF (TVG%CPHOTO=='LAI' .OR. TVG%CPHOTO=='LST' .OR. TVG%CPHOTO=='NIT') THEN
  CALL VEGETATION_EVOL(TVG%CISBA, TVG%CPHOTO, TVG%CRESPSL, TVG%CALBEDO, .FALSE., TVG%LTR_ML,   &
                       TVG%LNITRO_DILU, GAGRI_TO_GRASS,                        &
                       PTSTEP, TPTIME%TDATE%MONTH, TPTIME%TDATE%DAY, 1,    &
                       TPTIME%TIME, TG%XLAT, PRHOA, TGDP%XDG, TGDP%XDZG, TGDP%NWG_LAYER,     & 
                       TGD%CUR%XTG, TGDP%XALBNIR_VEG, TGDP%XALBVIS_VEG, TGDP%XALBUV_VEG,          &
                       TGDP%XALBNIR_SOIL, TGDP%XALBVIS_SOIL, TGDP%XALBUV_SOIL,            &
                       TGDP%XVEGTYPE, TGDP%XSEFOLD, TGDP%XANMAX, TGDP%XH_TREE, TGDP%XBSLAI,         &
                       TGDP%XLAIMIN, PCO2, TGDP%XCE_NITRO, TGDP%XCF_NITRO, TGDP%XCNA_NITRO,    &
                       TGDP%XBSLAI_NITRO, TGDP%XGMES, ZTAU_WOOD, TPSEED,             &
                       TPREAP, ZAOSIP, ZAOSIM, ZAOSJP, ZAOSJM,             &
                       ZHO2IP, ZHO2IM, ZHO2JP, ZHO2JM, ZZ0EFFIP,           &
                       ZZ0EFFIM, ZZ0EFFJP, ZZ0EFFJM, TGDPE%CUR%XLAI, TGDPE%CUR%XVEG,           &
                       TGDPE%CUR%XZ0, TGDPE%CUR%XALBNIR, TGDPE%CUR%XALBVIS, TGDPE%CUR%XALBUV, TGDPE%CUR%XEMIS,               &
                       TGD%CUR%XANFM, TGD%CUR%XANDAY, TGD%CUR%XBIOMASS, TGD%CUR%XRESP_BIOMASS,             &
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
IF (TVG%CPHOTO/='NON' .AND. TVG%CRESPSL/='NON' .AND. ANY(TGDPE%CUR%XLAI(:)/=XUNDEF)) THEN
  CALL CARBON_EVOL(TVG%CISBA, TVG%CRESPSL, TVG%CPHOTO, PTSTEP, 1,             &
                   PRHOA, TGD%CUR%XTG, TGD%CUR%XWG, TGDP%XWFC, TGDP%XWWILT, TGDP%XWSAT, TGDP%XSAND,   &
                   TGDP%XDG, TGDP%XDZG, TGDP%NWG_LAYER,                          &                   
                   TGDP%XRE25, TGDPE%CUR%XLAI, ZRESP_BIOMASS_INST, ZTURNOVER,    &
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
 CALL FLAG_TEB_GARDEN_n(TGD, TGDO, TGDPE, T, TVG, &
                        2)
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
   PLE_GARDEN    (:) = TGD%CUR%XLE    
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
  PQSAT_GARDEN(:) = QSAT(TGD%CUR%XTG(:,1),PPS(:))
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
