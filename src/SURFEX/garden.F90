!     #########
    SUBROUTINE GARDEN(HIMPLICIT_WIND, TPTIME, PPEW_A_COEF, PPEW_B_COEF,              &
                PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,                  &
                PTSTEP, PZ_LOWCAN,                                                   &
                PT_LOWCAN, PQ_LOWCAN, PEXNS, PRHOA, PCO2, PPS, PRR, PSR, PZENITH,    &
                PSW, PLW, PU_LOWCAN,                                                 &
                PRN_GARDEN,PH_GARDEN,PLE_GARDEN,PGFLUX_GARDEN,PSFCO2,                &
                PEVAP_GARDEN, PUW_GARDEN,                                            &
                PAC_GARDEN,PQSAT_GARDEN,PTS_GARDEN,                                  &
                PAC_AGG_GARDEN, PHU_AGG_GARDEN                                       )  
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
!!	A. Lemonsu          * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!    Original    05/2009
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_TYPE_DATE_SURF,    ONLY: DATE_TIME
USE MODD_SURF_PAR,          ONLY: XUNDEF
USE MODD_CSTS,              ONLY: XCPD
USE MODD_TEB_n,             ONLY: XCOVER, XGARDEN
USE MODD_TEB_GRID_n,        ONLY: XLAT, XLON
USE MODD_TEB_GARDEN_n,      ONLY: LPAR_GARDEN, &
                                  CISBA, CPHOTO, LTR_ML, CRUNOFF, CC1DRY,  &
                                  CSCOND, NNBIOMASS, CRESPSL, CALBEDO,     &
                                  CSOILFRZ, CDIFSFCOND, CCPSURF,           &
                                  CKSAT, CSOC, CHORT, CSNOWRES, TSNOW,     &
                                  XEMIS, XVEG, XLAI, XWRMAX_CF, XRSMIN,    &
                                  XGAMMA, XCV, XRGL, XRUNOFFD,             &
                                  XZ0, XZ0_O_Z0H, XRUNOFFB, XWDRAIN,       &
                                  XCGSAT, XC1SAT, XC2REF, XC3, XC4B,       &
                                  XC4REF, XACOEF, XPCOEF, XTAUICE,         &
                                  XTDEEP, XGAMMAT, XWR, XRESA, XAN,        &
                                  XANFM, XANDAY, XABC, XPOI,               &
                                  XFZERO, XEPSO, XGAMM, XQDGAMM,           &
                                  XGMES, XQDGMES, XT1GMES, XT2GMES,        &
                                  XRESP_BIOMASS, XBSLAI, XLAIMIN, XSEFOLD, &
                                  XLITTER, XSOILCARB, XLIGNIN_STRUC,       &
                                  XAMAX, XQDAMAX, XT1AMAX, XT2AMAX,        &
                                  LSTRESS, XF2I, XGC, XAH, XBH, XDMAX,     &
                                  XDG, XROOTFRAC, XDZG, XDZDIF, NWG_LAYER, &
                                  XTG, XWG, XWGI, XPCPS,                   &
                                  XPLVTT, XPLSTT, XWFC, XW33, XWWILT,XWSAT,&
                                  XBCOEF, XCONDSAT, XMPOTSAT, XHCAPSOIL,   &
                                  XZ0EFFIP, XZ0EFFIM, XZ0EFFJP, XZ0EFFJM,  &
                                  XAOSIP, XAOSIM, XAOSJP, XAOSJM,          &
                                  XHO2IP, XHO2IM, XHO2JP, XHO2JM,          &
                                  XCE_NITRO, XCF_NITRO, XCNA_NITRO,        &
                                  XCONDDRY, XCONDSLD, XRE25, TSEED, TREAP, &
                                  XIRRIG, XWATSUP, XCGMAX,                 &
                                  XKSAT_ICE, XD_ICE, XTURNOVER,            &
                                  XPATCH, XVEGTYPE_PATCH,                  &
                                  XALBNIR, XALBVIS, XALBUV,                &
                                  XALBNIR_VEG, XALBVIS_VEG, XALBUV_VEG,    &
                                  XALBNIR_SOIL, XALBVIS_SOIL, XALBUV_SOIL, &
                                  XLE, XANF, XSAND, XSOILWGHT,             &
                                  XPSN, XPSNV, XPSNG, XPSNV_A,             &
                                  XALBNIR_TVEG, XALBVIS_TVEG,              &
                                  XALBNIR_TSOIL, XLAI_EFFC,                &
                                  XALBVIS_TSOIL, XFAPARC, XFAPIRC, XMUS,   &
                                  NLAYER_HORT, NLAYER_DUN,                 &
                                  LSPINUPCARBS, LSPINUPCARBW, XSPINMAXS,   &
                                  XSPINMAXW, NNBYEARSPINS, NNBYEARSPINW,   &
                                  NNBYEARSOLD, NSPINS, NSPINW
!
USE MODD_AGRI_GARDEN, ONLY : LAGRIP
USE MODD_AGRI_GARDEN_n,     ONLY: LIRRIGATE, LIRRIDAY, XTHRESHOLDSPT
USE MODD_DIAG_TEB_GARDEN_n, ONLY: XCG, XC1, XC2, XWGEQ, XCT, XRS,           &
                                  XCH, XCD, XCDN, XRI, XHU, XHUG,           &
                                  XRN, XH, XLEI, XLEG, XLEGI, XLEV, XLES,   &
                                  XLER,  XLETR, XEVAP, XGFLUX,              &
                                  XRESTORE, XUSTAR, XDRAIN, XRUNOFF, XMELT, &
                                  XSNOWTEMP, XSNOWLIQ, XSNOWDZ, XSNOWHMASS, &
                                  XMELTADV, XIACAN, XQS, XHV, XHORT, XDRIP, &
                                  XSNOWFREE_ALB_VEG, XSNOWFREE_ALB_SOIL,    &
                                  XSNOWFREE_ALB, XTS, XTSRAD, XRRVEG, XALBT,&
                                  XEMIST, XGPP, XRESP_AUTO, XRESP_ECO,      &
                                  XFAPAR, XFAPIR, XFAPAR_BS, XFAPIR_BS,     &
                                  XDFAPARC, XDFAPIRC, XDLAI_EFFC,           &
                                  XIRRIG_FLUX  
!
USE MODI_ISBA
USE MODI_VEGETATION_UPDATE
USE MODE_THERMOS
!
USE MODI_FLAG_TEB_GARDEN_n
USE MODI_FLAG_DIAG_TEB_GARDEN
USE MODI_CARBON_EVOL
USE MODI_CARBON_SPINUP
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
!
TYPE(DATE_TIME)     , INTENT(IN)    :: TPTIME             ! current date and time from teb
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
REAL, DIMENSION(:),   INTENT(IN)    :: PSW                ! incoming total solar rad on an horizontal surface
REAL, DIMENSION(:)  , INTENT(IN)    :: PLW                ! atmospheric infrared radiation
REAL, DIMENSION(:)  , INTENT(IN)    :: PU_LOWCAN          ! wind near the road

REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN_GARDEN         ! net radiation over green areas
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_GARDEN          ! sensible heat flux over green areas
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_GARDEN         ! latent heat flux over green areas
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLUX_GARDEN      ! flux through the green areas
REAL, DIMENSION(:)  , INTENT(OUT)   :: PSFCO2             ! flux of CO2 positive toward the atmosphere (kg/m2/s)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PEVAP_GARDEN       ! total evaporation over gardens (kg/m2/s)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PUW_GARDEN         ! friction flux (m2/s2)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_GARDEN         ! aerodynamical conductance
REAL, DIMENSION(:)  , INTENT(OUT)   :: PQSAT_GARDEN       ! saturation humidity
REAL, DIMENSION(:)  , INTENT(INOUT) :: PTS_GARDEN         ! radiative surface temp. (snow free)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_AGG_GARDEN     ! aggreg. aeodynamic resistance for green areas for latent heat flux
REAL, DIMENSION(:)  , INTENT(OUT)   :: PHU_AGG_GARDEN     ! aggreg. relative humidity for green areas for latent heat flux
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
REAL, DIMENSION(0)   ::  PSODELX   ! Pulsation for each layer (Only used if LTEMP_ARP=True)                                    
REAL, DIMENSION(SIZE(PPS))  :: PMUF  ! fraction of the grid cell reached by the rainfall
REAL, DIMENSION(SIZE(PPS))  :: PFSAT ! Topmodel saturated fraction
REAL, DIMENSION(SIZE(PPS)) :: ZDIRCOSZW           ! orography slope cosine (=1 in TEB)
REAL, DIMENSION(SIZE(PPS),NNBIOMASS) :: ZRESP_BIOMASS_INST       ! instantaneous biomass respiration (kgCO2/kgair m/s)
!
!  temperatures
!
REAL, DIMENSION(SIZE(PPS)) :: ZTA ! estimate of air temperature at future time
!                                 ! step as if modified by ISBA flux alone.
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
REAL, DIMENSION(SIZE(PPS)) :: ZSMELTFLUX    ! energy removed from soil/vegetation surface
REAL, DIMENSION(SIZE(PPS)) :: ZGFLUXSNOW    ! net surface energy flux into snowpack (ISBA-ES:3-L)(W/m2)
REAL, DIMENSION(SIZE(PPS)) :: ZUSTARSNOW    ! friction velocity  over snow (ISBA-ES:3-L)    (m/s)
REAL, DIMENSION(SIZE(PPS)) :: ZGRNDFLUX     ! soil/snow interface heat flux (ISBA-ES:3-L)   (W/m2)
REAL, DIMENSION(SIZE(PPS)) :: ZSRSFC        ! snowfall over snowpack (ISBA-ES:3-L)          (kg/m2/s)
REAL, DIMENSION(SIZE(PPS)) :: ZRRSFC        ! rainfall over snowpack (ISBA-ES:3-L)          (kg/m2/s)
REAL, DIMENSION(SIZE(PPS)) :: ZLESL         ! snowpack evaporation (ISBA-ES:3-L)            (W/m2)
REAL, DIMENSION(SIZE(PPS)) :: ZCDSNOW       ! snow drag coefficient (ISBA-ES:3-L)           (-)
REAL, DIMENSION(SIZE(PPS)) :: ZCHSNOW       ! heat turbulent transfer coefficient           (-)
!
!  surfaces relative fractions
!
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
REAL :: ZTIME
!
LOGICAL :: GMASK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*      1.     various initialisations
!              -----------------------
!
IF (LHOOK) CALL DR_HOOK('GARDEN',0,ZHOOK_HANDLE)
ZDIRCOSZW = 1.
!
HRAIN     = 'DEF'
OFLOOD    = .FALSE.
OTEMP_ARP = .FALSE.
OGLACIER  = .FALSE.
PMUF      = 0.
PFSAT     = 0.
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
!-------------------------------------------------------------------------------
!
!*      2.     Treatment of green areas
!              ------------------------
!
!* Actualization of soil and wood carbon spinup
!
IF(LSPINUPCARBS.OR.LSPINUPCARBW)THEN
  ZTIME=TPTIME%TIME-PTSTEP
  CALL CARBON_SPINUP(TPTIME%TDATE%MONTH,TPTIME%TDATE%DAY,ZTIME,          &
                     LSPINUPCARBS, LSPINUPCARBW, XSPINMAXS, XSPINMAXW,   &
                     NNBYEARSPINS, NNBYEARSPINW, NNBYEARSOLD, CPHOTO,    &
                     CRESPSL, NSPINS, NSPINW                             )
ENDIF
!
!radiative temperature diagnostic
!-------------------------------
!
XTSRAD = PTS_GARDEN
!
!*      2.2    Call ISBA for green areas
!              -------------------------
!
  CALL ISBA(CISBA, CPHOTO, LTR_ML, CRUNOFF, CKSAT, CSOC, HRAIN, CHORT, &
            CC1DRY, CSCOND, TSNOW%SCHEME, CSNOWRES, CCPSURF, CSOILFRZ, &
            CDIFSFCOND, TPTIME, OFLOOD, OTEMP_ARP, OGLACIER, PTSTEP,   &
            HIMPLICIT_WIND,                                            &
            XCGMAX, PZ_LOWCAN, PZ_LOWCAN, ZDIRCOSZW, PT_LOWCAN,        &
            PQ_LOWCAN, PEXNS, PRHOA, PPS, PEXNS, PRR, PSR, PZENITH,    &
            PSW, PLW, PU_LOWCAN, PPEW_A_COEF, PPEW_B_COEF, PPET_A_COEF,&
            PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF, XRSMIN(:,1),        &
            XRGL(:,1), XGAMMA(:,1), XCV(:,1), XRUNOFFD(:,1),           &
            XSOILWGHT(:,:,1), NLAYER_HORT, NLAYER_DUN,                 &
            XALBNIR_TVEG, XALBVIS_TVEG, XALBNIR_TSOIL, XALBVIS_TSOIL,  &
            XSNOWFREE_ALB, XWRMAX_CF(:,1), XVEG(:,1), XLAI(:,1),       &
            XEMIS(:,1), XZ0(:,1), XZ0(:,1)/XZ0_O_Z0H(:,1),             &
            XVEGTYPE_PATCH(:,:,1), XZ0(:,1), XRUNOFFB, XCGSAT,         &
            XC1SAT(:,1), XC2REF(:,1), XC3(:,:,1), XC4B, XC4REF(:,1),   &
            XACOEF, XPCOEF, XTAUICE, XWDRAIN, XTDEEP, XGAMMAT,         &
            XPSN(:,1), XPSNG(:,1), XPSNV(:,1), XPSNV_A(:,1),           &
            XSNOWFREE_ALB_VEG, XSNOWFREE_ALB_SOIL, XIRRIG(:,1),        &
            XWATSUP(:,1), XTHRESHOLDSPT(:,1), LIRRIGATE(:,1),          &
            LIRRIDAY(:,1), LSTRESS(:,1), XGC(:,1), XF2I(:,1),          &
            XDMAX(:,1), XAH(:,1), XBH(:,1), PCO2, XGMES(:,1), XPOI,    &
            XFZERO(:,1), XEPSO(:,1), XGAMM(:,1), XQDGAMM(:,1),         &
            XQDGMES(:,1), XT1GMES(:,1), XT2GMES(:,1), XAMAX(:,1),      &
            XQDAMAX(:,1), XT1AMAX(:,1), XT2AMAX(:,1), XABC, XDG(:,:,1),&
            XDZG(:,:,1), XDZDIF(:,:,1), NWG_LAYER(:,1),                &
            XROOTFRAC(:,:,1), XWFC, XW33, XWWILT, XWSAT, XBCOEF,       &
            XCONDSAT(:,:,1), XMPOTSAT, XHCAPSOIL, XCONDDRY, XCONDSLD,  &
            XD_ICE(:,1), XKSAT_ICE(:,1), PMUF, ZFF, ZFFG, ZFFV,        &
            ZFFG_NOSNOW, ZFFV_NOSNOW, ZFFROZEN, ZALBF, ZEMISF, ZFFLOOD,&
            ZPIFLOOD, ZIFLOOD, ZPFLOOD, ZLEFLOOD, ZLEIFLOOD, PSODELX,  &
            XLAT, XLON, XTG(:,:,1), XWG(:,:,1), XWGI(:,:,1),           &
            XPCPS(:,1), XPLVTT(:,1), XPLSTT(:,1), XWR(:,1), XRESA(:,1),&
            XANFM(:,1), PFSAT, TSNOW%ALB(:,1), TSNOW%WSNOW(:,:,1),     &
            TSNOW%HEAT(:,:,1), TSNOW%RHO(:,:,1), TSNOW%GRAN1(:,:,1),   &
            TSNOW%GRAN2(:,:,1), TSNOW%HIST(:,:,1),TSNOW%AGE(:,:,1),    &
            ZGRNDFLUX, ZHPSNOW, XSNOWHMASS, ZSMELTFLUX, ZRNSNOW,       &
            ZHSNOW, ZGFLUXSNOW, ZUSTARSNOW, ZSRSFC, ZRRSFC, ZLESL,     &
            TSNOW%EMIS(:,1), ZCDSNOW, ZCHSNOW, PTS_GARDEN, XTS, XHV,   &
            XQS, XSNOWTEMP, XSNOWLIQ, XSNOWDZ, XCG, XC1, XC2, XWGEQ,   &
            XCT, XCH, XCD, XCDN, XRI, XHU, XHUG, XEMIST, XALBT, XRS,   &
            XLE(:,1), XRN, XH, XLEI, XLEGI, XLEG, XLEV, XLES, XLER,    &
            XLETR, XEVAP, XGFLUX, XRESTORE, XUSTAR, XDRAIN, XRUNOFF,   &
            XMELT, XMELTADV, ZRN_ISBA, ZH_ISBA, ZLEG_ISBA, ZLEGI_ISBA, &
            ZLEV_ISBA, ZLETR_ISBA, ZUSTAR_ISBA, ZLER_ISBA, ZLE_ISBA,   &
            ZLEI_ISBA, ZGFLUX_ISBA, XHORT, XDRIP, XRRVEG,              &
            PAC_AGG_GARDEN, PHU_AGG_GARDEN, XFAPARC(:,1), XFAPIRC(:,1),&
            XMUS(:,1), XLAI_EFFC(:,1), XAN(:,1), XANDAY(:,1),          &
            ZRESP_BIOMASS_INST, XIACAN, XANF(:,1), XGPP, XFAPAR,       &
            XFAPIR, XFAPAR_BS, XFAPIR_BS, XIRRIG_FLUX                  )
!
IF (TSNOW%SCHEME=='3-L' .OR. TSNOW%SCHEME=='CRO') TSNOW%TS(:,1)=XSNOWTEMP(:,1)
!
IF (LTR_ML) THEN
  GMASK = ( TPTIME%TIME - PTSTEP < 0. ) .AND. ( TPTIME%TIME >= 0. )
  IF (GMASK) THEN
    XDFAPARC  (:) = XFAPARC   (:,1) / XMUS (:,1)
    XDFAPIRC  (:) = XFAPIRC   (:,1) / XMUS (:,1)
    XDLAI_EFFC(:) = XLAI_EFFC (:,1) / XMUS (:,1)
  ENDIF
ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Diagnostic of respiration carbon fluxes and soil carbon evolution
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
PSFCO2    (:)=0.
XRESP_ECO (:)=0.
XRESP_AUTO(:)=0.
!
IF ( CPHOTO/='NON' .AND. CRESPSL/='NON' .AND. ANY(XLAI(:,1)/=XUNDEF) ) THEN
  CALL CARBON_EVOL(CISBA, CRESPSL, CPHOTO, PTSTEP, NSPINS,                      &
                   PRHOA, XTG(:,:,1), XWG(:,:,1), XWFC, XWWILT, XWSAT, XSAND,   &
                   XDG(:,:,1),XDZG(:,:,1), NWG_LAYER(:,1),                      &
                   XRE25(:,1), XLAI(:,1), ZRESP_BIOMASS_INST, XTURNOVER(:,:,1), &
                   XLITTER(:,:,:,1), XLIGNIN_STRUC(:,:,1) , XSOILCARB(:,:,1),   &
                   XRESP_AUTO, XRESP_ECO                       )  
  ! calculation of vegetation CO2 flux (Positive toward the atmosphere)
  PSFCO2(:) = XRESP_ECO(:) - XGPP(:)
END IF
!
! --------------------------------------------------------------------------------------
! Vegetation update (in case of non-interactive vegetation):
! --------------------------------------------------------------------------------------
!
IF (CPHOTO=='NON' .OR. CPHOTO=='AGS' .OR. CPHOTO=='AST') THEN
     CALL VEGETATION_UPDATE(PTSTEP,TPTIME,XCOVER,                        &
                         CISBA,(.NOT. LPAR_GARDEN), CPHOTO, LAGRIP, 'GRD',       &
                         XLAI,XVEG,XZ0,                                  &
                         XALBNIR,XALBVIS,XALBUV,XEMIS,                   &
                         XRSMIN,XGAMMA,XWRMAX_CF,                        &
                         XRGL,XCV,                                       &
                         XGMES,XBSLAI,XLAIMIN,XSEFOLD,XGC,XDMAX,         &
                         XF2I, LSTRESS,                                  &
                         XAOSIP,XAOSIM,XAOSJP,XAOSJM,                    &
                         XHO2IP,XHO2IM,XHO2JP,XHO2JM,                    &
                         XZ0EFFIP,XZ0EFFIM,XZ0EFFJP,XZ0EFFJM,            &
                         CALBEDO, XALBNIR_VEG, XALBVIS_VEG, XALBUV_VEG,  &
                         XALBNIR_SOIL, XALBVIS_SOIL, XALBUV_SOIL,        &
                         XCE_NITRO, XCF_NITRO, XCNA_NITRO,               &
                         TSEED, TREAP, XWATSUP, XIRRIG                   )  
END IF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!*      4.     Set undefined values for points where there is no garden
!              --------------------------------------------------------
!
! This way, these points are clearly flaged, and one will not try to interpret
! the values for those points
!
CALL FLAG_TEB_GARDEN_n(2)
CALL FLAG_DIAG_TEB_GARDEN
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!*      9.     Fields required for TEB
!              -----------------------
!
!
WHERE (XGARDEN/=0.)
  !
  ! energy balance
  !
   PRN_GARDEN    (:) = XRN    (:)
   PH_GARDEN     (:) = XH     (:)
   PLE_GARDEN    (:) = XLE    (:,1)
   PGFLUX_GARDEN (:) = XGFLUX (:)
   PEVAP_GARDEN  (:) = XEVAP  (:)
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
  PQSAT_GARDEN(:) = QSAT(XTG(:,1,1),PPS(:))
  !
  !* friction flux
  PUW_GARDEN(:)    = -XUSTAR(:)**2
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
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE GARDEN
