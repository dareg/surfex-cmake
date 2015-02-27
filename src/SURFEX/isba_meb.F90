!     #########
      SUBROUTINE ISBA_MEB(TPTIME, OMEB, OFORC_MEASURE, OGLACIER,               &
        OSNOWDRIFT, OSNOWDRIFT_SUBLIM, OSNOW_ABS_ZENITH, OIRRIGATE, OIRRIDAY,  &
        HSNOWMETAMO, HSNOWRAD,                                                 &           
        HISBA, HCPSURF, HRAIN, HSNOW_ISBA, HSNOWRES, HIMPLICIT_WIND,           &
        KWG_LAYER, PTSTEP, PVEGTYPE, PLAT, PLON,                               &
        PTHRESHOLD, PWATSUP, PIRRIG, PIRRIG_FLUX,                              &
        PSOILHCAPZ, PSOILCONDZ, PFROZEN1,                                      &
        PPS, PZENITH, PSCA_SW, PSW_RAD, PVMOD, PRR, PSR, PRHOA, PTA, PQA,      &
        PH_VEG, PDIRCOSZW,                                                     &
        PEXNS, PEXNA, PPET_A_COEF, PPET_B_COEF, PPEQ_A_COEF, PPEQ_B_COEF,      &
        PPEW_A_COEF, PPEW_B_COEF,                                              &
        PZREF, PUREF, PCH, PCD, PCDN, PRI, PRESA, PHUG, PHV, PHU, PQS,         &
        PZ0G_WITHOUT_SNOW,                                                     &
        PZ0_MEBV, PZ0H_MEBV, PZ0EFF_MEBV,                                      &
        PZ0_MEBN, PZ0H_MEBN, PZ0EFF_MEBN,                                      &
        PZ0_WITH_SNOW, PZ0H_WITH_SNOW, PZ0EFF,                                 &
        PTV, PTG, PTC, PQC, PWR, PWRVN, PWG, PWGI,                             &
        PWRMAX_CF, PRGL, PRSMIN, PGAMMA, PRS,                                  &
        PALBNIR_TVEG, PALBVIS_TVEG,PALBNIR_TSOIL, PALBVIS_TSOIL, PFALB,        &
        PSNOWALB, PSNOWALBVIS, PSNOWALBNIR, PSNOWALBFIR,                       &
        PGNDLITTER, PFF, PPSN, PPALPHAN, PZF_TALLVEG, PLAI, PROOTFRAC, PF2,    &
        PWSAT, PWFC, PWWILT,                                                   &
        PSNOWGRAN1, PSNOWGRAN2, PSNOWHIST,PSNOWAGE,                            &
        PSNOWRHO, PSNOWSWE, PSNOWHEAT, PSNOWTEMP, PSNOWDZ, PSNOWLIQ, PFEMIS,   &
        PSWNET_N, PSWNET_V, PSWNET_G, PSWNET_NS, PALBT, PSWDOWN_GN,            &
        PLW_RAD, PLWNET_N, PLWNET_V, PLWNET_G, PLWDOWN_GN,                     &
        PLEV_V_C, PLES_V_C, PH_V_C, PH_G_C, PLETR_V_C, PLER_V_C, PH_C_A,       &
        PH_N_C, PLE_V_C, PLE_G_C, PLE_C_A, PLE_N_C, PEVAP_N_C, PEVAP_G_C,      &
        PSR_GN, PMELTCV, PFRZCV, PMELT, PMELTADV,                              &
        PLE_FLOOD, PLEI_FLOOD,                                                 &
        PLE, PH, PRN, PLEI, PLEGI, PLEG, PLEV, PLES, PLER, PLETR, PEVAP,       &
        PSUBL, PGFLUX, PRESTORE, PGRNDFLUX, PFLSN_COR, PUSTAR,                 &
        PHPSNOW, PSNOWHMASS, PRNSNOW, PHSNOW, PGFLUXSNOW,                      &
        PUSTARSNOW, PSRSFC, PRRSFC, PLESL, PEMISNOW, PCDSNOW, PCHSNOW,         &
        PEMIST, PTS_RAD, PHU_AGG, PAC_AGG,                                     &
        PDELHEATV_SFC, PDELHEATG_SFC, PDELHEATG,                               &
        PDELHEATN, PDELHEATN_SFC, PRESTOREN,                                   &
        PD_G, PCPS, PLVTT, PLSTT, PCT, PCV, PCG, PFFROZEN,                     &
        PTDEEP_A, PTDEEP_B, PDEEP_FLUX, PMUF, PDRIP, PRRVEG,                   &
        PRISNOW, PSNOW_THRUFAL, PEVAPCOR, PSUBVCOR, PSNOWSFCH, PSNDRIFT, PQSNOW)
!     ##########################################################################
!
!                             
!!****  *isba_meb*  
!!
!!    PURPOSE
!!    -------
!       Monitor for the calculation of the surface fluxes and of the
!     prognostic variables of the surface over natural areas 
!     with an explicit vegetation layer
!
!     NOTE...currently MEB can be coupled with 
!     HISBA='DIF' or '3-L' soil options
!     HSNOW='3-L' snow scheme
!     Soon, HSNOW=CRO and HPHOTO/=NON (i.e. Ags will be added)
!     
!!**  METHOD
!!    ------
!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!    Noilhan and Planton (1989)
!!      
!!    AUTHOR
!!    ------
!!      A. Boone           * Meteo-France *
!!      P. Samuelsson      * SMHI *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2014
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_CSTS,           ONLY : XCPD, XDAY 
!
USE MODD_TYPE_DATE_SURF, ONLY : DATE_TIME
!
USE MODE_THERMOS
USE MODE_MEB,            ONLY : SNOW_INTERCEPT_EFF
!
USE MODI_WET_LEAVES_FRAC
USE MODI_VEG
USE MODI_SNOW_LEAVES_FRAC_MEB
USE MODI_PREPS_FOR_MEB_EBUD_RAD
USE MODI_ISBA_SWNET_MEB
USE MODI_ISBA_LWNET_MEB
USE MODI_DRAG_MEB
USE MODI_E_BUDGET_MEB
USE MODI_ISBA_FLUXES_MEB
USE MODI_SNOW_LOAD_MEB
USE MODI_HYDRO_VEG
USE MODI_SNOW3L_ISBA

!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!              -------------------------
!
!
!* general variables
!  -----------------
!
TYPE(DATE_TIME),      INTENT(IN)    :: TPTIME        ! current date and time
!
LOGICAL,              INTENT(IN)    :: OMEB          ! True = patch with multi-energy balance 
!                                                    ! False = patch with classical ISBA 
LOGICAL,              INTENT(IN)    :: OFORC_MEASURE ! switch for using measured data (drag scheme)
LOGICAL,              INTENT(IN)    :: OGLACIER      ! True = Over permanent snow and ice, 
!                                                    ! initialise WGI=WSAT,
!                                                    ! Hsnow>=10m and allow 0.8<SNOALB<0.85
!                                                    ! False = No specific treatment
LOGICAL,              INTENT(IN)    :: OSNOWDRIFT    ! if=T, activate snowdrift
LOGICAL,              INTENT(IN)    :: OSNOWDRIFT_SUBLIM ! if=T, activate snowdrift sublimation 
LOGICAL,              INTENT(IN)    :: OSNOW_ABS_ZENITH  ! if=T, activate parametrization of solar absorption 
!                                                        ! for polar regions
LOGICAL, DIMENSION(:),INTENT(IN)    :: OIRRIGATE     ! Irrigation FLAG
LOGICAL, DIMENSION(:),INTENT(INOUT) :: OIRRIDAY      ! Irrigation time 
!
CHARACTER(LEN=*),     INTENT(IN)    :: HISBA         ! type of ISBA version:
!                                                    ! '2-L' (default)
!                                                    ! '3-L'
!                                                    ! 'DIF'
CHARACTER(LEN=*),     INTENT(IN)    :: HCPSURF       ! Specific heat
!                                                    ! 'DRY' = dry Cp
!                                                    ! 'HUM' = humid Cp fct of qs
CHARACTER(LEN=*),     INTENT(IN)    :: HRAIN         ! Rainfall spatial distribution
                                                     ! 'DEF' = No rainfall spatial distribution
                                                     ! 'SGH' = Rainfall exponential spatial distribution
CHARACTER(LEN=*),     INTENT(IN)    :: HSNOW_ISBA    ! 'DEF' = Default F-R snow scheme
!                                                    !         (Douville et al. 1995)
!                                                    ! '3-L' = 3-L snow scheme (option)
!                                                    !         (Boone and Etchevers 2000)
CHARACTER(LEN=*),     INTENT(IN)    :: HSNOWRES      ! 'DEF' = Default: Louis (ISBA)
!                                                    ! 'RIL' = CROCUS (Martin) method
!                                                    !  ISBA-SNOW3L turbulant exchange option
CHARACTER(LEN=*),     INTENT(IN)    :: HIMPLICIT_WIND! wind implicitation option
!                                                    ! 'OLD' = direct
!                                                    ! 'NEW' = Taylor serie, order 1
CHARACTER(LEN=*),     INTENT(IN)    :: HSNOWMETAMO   ! Crocus metamorphism scheme:
!                                                    ! HSNOWMETAMO = B92 Brun et al 1992
!                                                    ! HSNOWMETAMO = C13 Carmagnola et al 2014
!                                                    ! HSNOWMETAMO = T07 Taillandier et al 2007
!                                                    ! HSNOWMETAMO = F06 Flanner et al 2006
CHARACTER(LEN=*),     INTENT(IN)    :: HSNOWRAD      ! Crocus radiative transfer scheme:
!                                                    ! HSNOWMETAMO = B92 Brun et al 1992
!                                                    ! HSNOWMETAMO = TAR TARTES (Libois et al 2013)
!                                                    ! HSNOWMETAMO = TA1 TARTES with constant impurities
!                                                    ! HSNOWMETAMO = TA2 TARTES with constant impurities as a 
!                                                    !                   function of ageing
!
INTEGER, DIMENSION(:),INTENT(IN)    :: KWG_LAYER     ! Number of soil moisture layers (DIF option)
!
REAL,                 INTENT(IN)    :: PTSTEP        ! Model time step (s)
REAL, DIMENSION(:,:), INTENT(IN)    :: PVEGTYPE      ! fraction of each vegetation (-)
REAL, DIMENSION(:),   INTENT(IN)    :: PLAT          ! Latitude (degrees North)
REAL, DIMENSION(:),   INTENT(IN)    :: PLON          ! Longitude (degrees East)
REAL, DIMENSION(:),   INTENT(IN)    :: PPS           ! Pressure [Pa]
REAL, DIMENSION(:),   INTENT(IN)    :: PZENITH       ! solar zenith angle
REAL, DIMENSION(:),   INTENT(IN)    :: PSW_RAD       ! solar (shortwave) incoming radiation [W/m2]
REAL, DIMENSION(:),   INTENT(IN)    :: PLW_RAD       ! thermal (longwave) incoming radiation [W/m2]
REAL, DIMENSION(:),   INTENT(IN)    :: PSCA_SW       ! solar diffuse incoming radiation [W/m2]
REAL, DIMENSION(:),   INTENT(IN)    :: PEXNA         ! Exner function: forcing level (-)
REAL, DIMENSION(:),   INTENT(IN)    :: PEXNS         ! Exner function: surface (-)
REAL, DIMENSION(:),   INTENT(IN)    :: PRR           ! Rain rate (kg/m2/s)
REAL, DIMENSION(:),   INTENT(IN)    :: PSR           ! Snow rate (kg/m2/s)
REAL, DIMENSION(:),   INTENT(IN)    :: PRHOA         ! air density (kg/m3)
REAL, DIMENSION(:),   INTENT(IN)    :: PVMOD         ! modulus of the wind
!                                                    ! parallel to the orography (m/s)
REAL, DIMENSION(:),   INTENT(IN)    :: PTA           ! Temperature of atmosphere (K)
REAL, DIMENSION(:),   INTENT(IN)    :: PQA           ! specific humidity of atmosphere (kg/kg)
REAL, DIMENSION(:),   INTENT(IN)    :: PH_VEG        ! height of vegetation
REAL, DIMENSION(:),   INTENT(IN)    :: PZREF         ! normal distance of the first
!                                                    ! atmospheric level to the
!                                                    ! orography (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PUREF         ! reference height of the wind (m)
!                                                    ! NOTE this is different from ZZREF
!                                                    ! ONLY in stand-alone/forced mode,
!                                                    ! NOT when coupled to a model (MesoNH)
REAL, DIMENSION(:),   INTENT(IN)    :: PDIRCOSZW     ! Director Cosinus along the z
!                                                    ! direction at the surface w-point
REAL, DIMENSION(:,:), INTENT(IN)    :: PSOILHCAPZ    ! ISBA-DF Soil heat capacity 
!                                                    ! profile [J/(m3 K)]
REAL, DIMENSION(:,:), INTENT(IN)    :: PSOILCONDZ    ! ISBA-DF Soil conductivity  
!                                                    ! profile  [W/(m K)]
REAL, DIMENSION(:),   INTENT(IN)    :: PFROZEN1      ! surface frozen fraction (-)
REAL, DIMENSION(:),   INTENT(IN)    :: PZF_TALLVEG   ! fraction of tall vegetation (compared to short veg) for
!                                                    !  shortwave radiative transfer (-)
REAL, DIMENSION(:),   INTENT(IN)    :: PLAI          ! vegetation Leaf Area Index (m2/m2)
REAL, DIMENSION(:),   INTENT(IN)    :: PGNDLITTER    ! fraction of litter on the surface
                                                     ! as==>0, baresoil below canopy,
                                                     ! as==>1, litter layer below canopy
REAL, DIMENSION(:),   INTENT(IN)    :: PRGL          ! maximum solar radiation
!                                                    ! usable in photosynthesis
REAL, DIMENSION(:),   INTENT(IN)    :: PRSMIN        ! minimum stomatal resistance (s/m)
REAL, DIMENSION(:),   INTENT(IN)    :: PGAMMA        ! coefficient for the calculation
!                                                    ! of the surface stomatal resistance
REAL, DIMENSION(:),   INTENT(IN)    :: PFF           ! Floodplain fraction at the surface
REAL, DIMENSION(:),   INTENT(IN)    :: PPSN          ! fraction of the grid covered
!                                                    ! by snow
REAL, DIMENSION(:),   INTENT(IN)    :: PPALPHAN      ! snow/canopy transition coefficient
REAL, DIMENSION(:),   INTENT(IN)    :: PFALB         ! Floodplain albedo
REAL, DIMENSION(:),   INTENT(IN)    :: PALBNIR_TVEG  ! albedo of vegetation in NIR 
!                                                    ! (needed for LM_TR or MEB)
REAL, DIMENSION(:),   INTENT(IN)    :: PALBVIS_TVEG  ! albedo of vegetation in VIS 
!                                                    ! (needed for LM_TR or MEB)
REAL, DIMENSION(:),   INTENT(IN)    :: PALBNIR_TSOIL ! albedo of bare soil in NIR 
!                                                    ! (needed for LM_TR or MEB)
REAL, DIMENSION(:),   INTENT(IN)    :: PALBVIS_TSOIL ! albedo of bare soil in VIS 
REAL, DIMENSION(:),   INTENT(IN)    :: PWRMAX_CF     ! maximum vegetation interception storage (kg/m2) 
REAL, DIMENSION(:),   INTENT(IN)    :: PFEMIS        ! Floodplain emissivity (-)
REAL, DIMENSION(:),   INTENT(IN)    :: PF2           ! Soil water stress factor for transpiration (-)
REAL, DIMENSION(:,:), INTENT(IN)    :: PROOTFRAC     ! vegetation root fraction (-)
REAL, DIMENSION(:,:), INTENT(IN)    :: PWFC          ! field capacity profile               (m3/m3)
REAL, DIMENSION(:,:), INTENT(IN)    :: PWWILT        ! wilting point profile                (m3/m3)
REAL, DIMENSION(:,:), INTENT(IN)    :: PWSAT         ! porosity profile                     (m3/m3)
REAL, DIMENSION(:,:), INTENT(IN)    :: PWG, PWGI     ! PWG  = soil liquid volumetric water content (m3/m3)
!                                                    ! PWGI = soil frozen volumetric water content (m3/m3)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0G_WITHOUT_SNOW ! roughness length for momentum at snow-free canopy floor (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0_MEBV      ! roughness length for momentum over MEB vegetation part of patch (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0H_MEBV     ! roughness length for heat over MEB vegetation part of path (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0EFF_MEBV   ! roughness length for momentum over MEB vegetation part of patch (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0_MEBN      ! roughness length for momentum over MEB snow part of patch (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0H_MEBN     ! roughness length for heat over MEB snow part of path (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0EFF_MEBN   ! roughness length for momentum over MEB snow part of patch (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0_WITH_SNOW ! roughness length for momentum
!                                                    ! (with snow taken into account) (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0H_WITH_SNOW ! roughness length for heat
!                                                    ! (with snow taken into account) (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0EFF        ! roughness length for momentum (m)
REAL, DIMENSION(:,:), INTENT(IN)    :: PD_G          ! Depth of Bottom of Soil layers       (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PCT           ! area-averaged effective inverse heat capacity [(K m2)/J]
REAL, DIMENSION(:),   INTENT(IN)    :: PCV           ! vegetation inverse heat capacity [(K m2)/J]
REAL, DIMENSION(:),   INTENT(IN)    :: PCG           ! soil inverse heat capacity [(K m2)/J]
REAL, DIMENSION(:),   INTENT(IN)    :: PFFROZEN      ! Fraction of frozen flood (-)
REAL, DIMENSION(:),   INTENT(IN)    :: PMUF          ! fraction of the grid cell reached by the rainfall (-)
!
! implicit atmospheric coupling coefficients:
!
REAL, DIMENSION(:),   INTENT(IN)    :: PPET_A_COEF, PPET_B_COEF, &
                                       PPEQ_A_COEF, PPEQ_B_COEF, &
                                       PPEW_A_COEF, PPEW_B_COEF  
!                                                    ! PPEW_A_COEF  A-wind coefficient
!                                                    ! PPEW_B_COEF  B-wind coefficient
!                                                    ! PPET_A_COEF  A-air temperature coefficient
!                                                    ! PPET_B_COEF  B-air temperature coefficient
!                                                    ! PPEQ_A_COEF  A-air specific humidity coefficient
!                                                    ! PPEQ_B_COEF  B-air specific humidity coefficient
REAL, DIMENSION(:),   INTENT(IN)    :: PTDEEP_A, PTDEEP_B ! Deep soil temperature boundary condition 
!                                                         ! (prescribed)     
!                                      PTDEEP_A = Deep soil temperature
!                                                 coefficient depending on flux
!                                      PTDEEP_B = Deep soil temperature (prescribed)
!                                                 which models heating/cooling from
!                                                 below the diurnal wave penetration
!                                                 (surface temperature) depth. If it
!                                                 is FLAGGED as undefined, then the zero
!                                                 flux lower BC is applied.
!                                                 Tdeep = PTDEEP_B + PTDEEP_A * PDEEP_FLUX
!                                                 (with PDEEP_FLUX in W/m2)
!
REAL, DIMENSION(:),   INTENT(IN)    :: PTHRESHOLD, PWATSUP, PIRRIG
!                                      PTHRESHOLD = threshold water level for irrigation (-)
!                                      PWATSUP    = irrigation water need to maintain a given moisture thresold (kg/m2)
!                                      PIRRIG     = irrigation mask (-)
!
! - - - - - - - - - - - - - - - - - - - - 
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PSNOWALB      ! Snow albedo
REAL, DIMENSION(:),   INTENT(INOUT) :: PSNOWALBVIS   ! Snow VIS albedo
REAL, DIMENSION(:),   INTENT(INOUT) :: PSNOWALBNIR   ! Snow NIR albedo
REAL, DIMENSION(:),   INTENT(INOUT) :: PSNOWALBFIR   ! Snow FIR albedo
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWSWE      ! Snow model layer liquid water equivalent or 
!                                                    ! SWE (kg m-2)  
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWHEAT     ! Snow layer heat content (J/m3) 
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWRHO      ! Snow layer average density (kg/m3)
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWGRAN1    ! Snow grain parameter 1 
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWGRAN2    ! Snow grain parameter 2 
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWHIST     ! Snow grain historical parameter
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWAGE      ! Snow grain age
!                                                    ! NOTE : methamorphism is only activated if the flag
!                                                    ! OSNOW_METAMO=TRUE
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PTG           ! Soil layer average temperature (K)
REAL, DIMENSION(:),   INTENT(INOUT) :: PTV           ! Canopy vegetation temperature (K)
REAL, DIMENSION(:),   INTENT(INOUT) :: PTC           ! Canopy air temperature [K]
REAL, DIMENSION(:),   INTENT(INOUT) :: PQC           ! Canopy air specific humidity [kg/kg]
REAL, DIMENSION(:),   INTENT(INOUT) :: PWR           ! liquid water retained on the foliage
!                                                    ! of the canopy vegetation [kg/m2]
REAL, DIMENSION(:),   INTENT(INOUT) :: PWRVN         ! liquid water equiv of snow retained on the foliage
!                                                    ! of the canopy vegetation [kg/m2]
REAL, DIMENSION(:),   INTENT(INOUT) :: PRESA         ! aerodynamic resistance (s/m)
REAL, DIMENSION(:),   INTENT(INOUT) :: PLE           ! total latent heat flux (W/m2)
REAL, DIMENSION(:),   INTENT(INOUT) :: PLE_FLOOD     ! Floodplains latent heat flux: liquid part [W/m2]
REAL, DIMENSION(:),   INTENT(INOUT) :: PLEI_FLOOD    ! Floodplains latent heat flux: frozen part [W/m2]
!
! - - - - - - - - - - - - - - - - - - - - 
!
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWTEMP     ! Snow layer average temperature (K)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWDZ       ! Snow layer thickness (m)
REAL, DIMENSION(:),   INTENT(OUT)   :: PEMISNOW      ! Snow surface emissivity (-)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSWNET_N      ! net snow shortwave radiation [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PSWNET_NS     ! net snow shortwave radiation for 
!                                                    ! the *surface* snow layer 
REAL, DIMENSION(:),   INTENT(OUT)   :: PSWNET_V      ! net vegetation canopy shortwave radiation 
!                                                    ! [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PSWNET_G      ! net surface (ground) shortwave radiation [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PALBT         ! total surface albedo
REAL, DIMENSION(:),   INTENT(OUT)   :: PSWDOWN_GN    ! total shortwave radiation transmitted through 
                                                     ! the vegetation canopy
REAL, DIMENSION(:),   INTENT(OUT)   :: PLWNET_V      ! net vegetation canopy longwave radiation [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PLWNET_G      ! net ground longwave radiation [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PLWNET_N      ! net snow longwave radiation [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PLWDOWN_GN    ! total shortwave radiation transmitted through and emitted by 
!                                                    !  the canopy reaching the snowpack/ground (explicit part) [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PRS           ! surface stomatal resistance (s/m)
REAL, DIMENSION(:),   INTENT(OUT)   :: PCH           ! drag coefficient for heat
REAL, DIMENSION(:),   INTENT(OUT)   :: PCD           ! drag coefficient for momentum
REAL, DIMENSION(:),   INTENT(OUT)   :: PCDN          ! neutral drag coefficient for momentum
REAL, DIMENSION(:),   INTENT(OUT)   :: PRI           ! Richardson number
REAL, DIMENSION(:),   INTENT(OUT)   :: PHV           ! Total effective Halstead coefficient
REAL, DIMENSION(:),   INTENT(OUT)   :: PHU           ! grid-area humidity of the soil
REAL, DIMENSION(:),   INTENT(OUT)   :: PHUG          ! ground relative humidity
REAL, DIMENSION(:),   INTENT(OUT)   :: PQS           ! surface humidity (kg/kg)
REAL, DIMENSION(:),   INTENT(OUT)   :: PRN           ! net radiation
REAL, DIMENSION(:),   INTENT(OUT)   :: PH            ! sensible heat flux
REAL, DIMENSION(:),   INTENT(OUT)   :: PLEI          ! sublimation latent heat flux
REAL, DIMENSION(:),   INTENT(OUT)   :: PLEGI         ! latent heat of sublimation over frozen soil
REAL, DIMENSION(:),   INTENT(OUT)   :: PLEG          ! latent heat of evaporation
!                                                    ! over the ground
REAL, DIMENSION(:),   INTENT(OUT)   :: PLEV          ! latent heat of evaporation
!                                                    ! over the vegetation
REAL, DIMENSION(:),   INTENT(OUT)   :: PLES          ! latent heat of sublimation
!                                                    ! over the snow
REAL, DIMENSION(:),   INTENT(OUT)   :: PLER          ! latent heat of the fraction
!                                                    ! delta of water retained on the
!                                                    ! foliage of the vegetation
REAL, DIMENSION(:),   INTENT(OUT)   :: PLETR         ! evapotranspiration of the rest
!                                                    ! of the vegetation
REAL, DIMENSION(:),   INTENT(OUT)   :: PEVAP         ! total evaporative flux (kg/m2/s)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSUBL         ! sublimation flux (kg/m2/s)
REAL, DIMENSION(:),   INTENT(OUT)   :: PGFLUX        ! flux through the ground
REAL, DIMENSION(:),   INTENT(OUT)   :: PRESTORE      ! surface restore flux for Force-Restore, diffusive flux between uppermost and second soil layers
!                                                    ! when using the DIF soil option (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PUSTAR        ! friction velocity
REAL, DIMENSION(:),   INTENT(OUT)   :: PMELT         ! melting rate of the snow (kg/m2/s)
REAL, DIMENSION(:),   INTENT(OUT)   :: PMELTADV      ! advection heat flux from snowmelt (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PCPS          ! heat capacity of air (J/kg/K)
REAL, DIMENSION(:),   INTENT(OUT)   :: PLVTT         ! latent heat of vaporization (J/kg)
REAL, DIMENSION(:),   INTENT(OUT)   :: PLSTT         ! latent heat of sublimation (J/kg)
REAL, DIMENSION(:),   INTENT(OUT)   :: PLEV_V_C      ! MEB: total evapotranspiration (no sublim) from vegetation canopy overstory [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PLES_V_C      ! MEB: total snow sublimation from vegetation canopy overstory [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PH_V_C        ! MEB: sensible heat flux from vegetation canopy overstory [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PH_G_C        ! MEB: sensible heat flux from ground [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PLETR_V_C     ! MEB: transpiration from overstory canopy vegetation [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PLER_V_C      ! MEB: interception evaporation from overstory canopy vegetation [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PH_C_A        ! MEB: sensible heat flux from canopy air space to the atmosphere [W/m2] 
                                                     !      NOTE total sensible heat flux to the atmosphere also possibly 
                                                     !      includes a contribution from snow covering the canopy
REAL, DIMENSION(:),   INTENT(OUT)   :: PH_N_C        ! MEB: sensible heat flux from the snow on the ground [W/m2]
                                                     !      NOTE total sensible heat flux from the snowpack
                                                     !      possibly includes a contribution from snow covering the canopy
REAL, DIMENSION(:),   INTENT(OUT)   :: PLE_V_C       ! MEB: latent heat flux from vegetation canopy overstory [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PLE_G_C       ! MEB: latent heat flux from ground [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PLE_C_A       ! MEB: latent heat flux from canopy air space to the atmosphere [W/m2] 
                                                     !      NOTE total latent heat flux to the atmosphere also possibly 
                                                     !      includes a contribution from snow covering the canopy
REAL, DIMENSION(:),   INTENT(OUT)   :: PLE_N_C       ! MEB: latent heat flux from the snow on the ground [W/m2]
                                                     !      NOTE total latent heat flux from the snowpack
                                                     !      possibly includes a contribution from snow covering the canopy
REAL, DIMENSION(:),   INTENT(OUT)   :: PEVAP_N_C     ! MEB: Total evap from snow on the ground to canopy air space  [kg/m2/s]
REAL, DIMENSION(:),   INTENT(OUT)   :: PEVAP_G_C     ! MEB: Total evap from ground to canopy air space [kg/m2/s]
REAL, DIMENSION(:),   INTENT(OUT)   :: PSR_GN        ! MEB: total snow reaching the ground snow [kg/m2/s]
REAL, DIMENSION(:),   INTENT(OUT)   :: PMELTCV       ! MEB: snow melt rate from the overstory snow reservoir [kg/m2/s]
REAL, DIMENSION(:),   INTENT(OUT)   :: PFRZCV        ! MEB: snow refreeze rate from the overstory snow reservoir [kg/m2/s]
REAL, DIMENSION(:),   INTENT(OUT)   :: PGRNDFLUX     ! snow/soil-biomass interface flux (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PFLSN_COR     ! soil/snow interface correction flux to conserve energy (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PHPSNOW       ! heat release from rainfall (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOWHMASS    ! snow heat content change from mass changes (J/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PRNSNOW       ! net radiative flux from snow (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PHSNOW        ! sensible heat flux from snow (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PGFLUXSNOW    ! net heat flux from snow (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PUSTARSNOW    ! friction velocity
REAL, DIMENSION(:),   INTENT(OUT)   :: PSRSFC        ! Snow rate falling outside of snow
!                                                    !  covered grid area [kg/(m2 s)]
REAL, DIMENSION(:),   INTENT(OUT)   :: PRRSFC        ! Rain rate falling outside of snow and flood
!                                                    !  covered grid area [kg/(m2 s)]
REAL, DIMENSION(:),   INTENT(OUT)   :: PLESL         ! Evaporation (liquid) from wet snow (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PCDSNOW       ! drag coefficient for momentum over snow
REAL, DIMENSION(:),   INTENT(OUT)   :: PCHSNOW       ! drag coefficient for heat over snow
REAL, DIMENSION(:),   INTENT(OUT)   :: PEMIST        ! total effective surface emissivity...LWUP = EMIST*TS_RAD**4 (-)
REAL, DIMENSION(:),   INTENT(OUT)   :: PTS_RAD       ! effective radiative temperature 
!                                                    !  of the natural surface (K)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWLIQ      ! snow layer liquid water content (m)
REAL, DIMENSION(:),   INTENT(OUT)   :: PAC_AGG       ! aggregated aerodynamic conductance
                                                     ! for evaporative flux calculations
REAL, DIMENSION(:),   INTENT(OUT)   :: PHU_AGG       ! aggregated relative humidity
                                                     ! for evaporative flux calculations
REAL, DIMENSION(:),   INTENT(OUT)   :: PDELHEATV_SFC ! change in heat storage of the vegetation canopy layer over the current time step (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PDELHEATG_SFC ! change in heat storage of the ground sfc layer over the current time step (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PDELHEATG     ! change in heat storage of the entire soil column over the current time step (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PRESTOREN     ! conductive heat flux between the surface and sub-surface soil layers 
!                                                    ! for the multi-layer snow schemes..for composite snow, it is 
!                                                    ! equal to PRESTORE (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PDELHEATN     ! change in heat storage of the entire snow column over the current time step (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PDELHEATN_SFC ! change in heat storage of the surface snow layer over the current time step (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PDEEP_FLUX    ! Heat flux at bottom of ISBA (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PDRIP         ! Water dripping from the vegetation canopy (kg/m2/s)
REAL, DIMENSION(:),   INTENT(OUT)   :: PRRVEG        ! Water intercepted by the vegetation canopy (kg/m2/s)
REAL, DIMENSION(:),   INTENT(OUT)   :: PRISNOW       ! Richarson number over ground-based snowpack (-)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOW_THRUFAL ! rate that liquid water leaves (explicit) snow pack: 
!                                                    ! ISBA-ES or CROCUS [kg/(m2 s)]
REAL, DIMENSION(:),   INTENT(OUT)   :: PEVAPCOR      !  evaporation correction as last traces of snow
!                                                    ! cover ablate..if sublimation exceeds trace amounts
                                                     ! of snow during time step, required residual mass taken 
                                                     ! from sfc soil layer [kg/(m2 s)]
REAL, DIMENSION(:),   INTENT(OUT)   :: PSUBVCOR      ! A possible snow mass correction (to be potentially    
!                                                    !  removed from soil)  (kg/m2/s)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOWSFCH     ! snow surface layer pseudo-heating term owing to
!                                                    !  changes in grid thickness            (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNDRIFT      ! blowing snow sublimation (kg/m2/s)
REAL, DIMENSION(:),   INTENT(OUT)   :: PQSNOW        ! snow surface specific humidity (kg/kg)
REAL, DIMENSION(:),   INTENT(OUT)   :: PIRRIG_FLUX   ! (kg/m2/s) irrigation flux (water need)
!
!
!*      0.2    declarations of local variables
!
!
REAL, PARAMETER                                    :: ZTSTEP_EB     = 300. ! s Minimum time tstep required 
!                                                                          !   to time-split MEB energy budget
INTEGER                                            :: JTSPLIT_EB           ! number of time splits
INTEGER                                            :: JDT                  ! time split loop index
!
REAL                                               :: ZTSTEP               ! Local time split timestep (s)
REAL, DIMENSION(SIZE(PPS))                         :: ZWORK                ! Working variable [*]
REAL, DIMENSION(SIZE(PSNOWSWE,1),SIZE(PSNOWSWE,2)) :: ZSNOWCOND            ! snow thermal conductivity  [W/(m K)] 
REAL, DIMENSION(SIZE(PSNOWSWE,1),SIZE(PSNOWSWE,2)) :: ZSNOWHCAP            ! snow heat capacity [J/(m3 K)]
REAL, DIMENSION(SIZE(PPS))                         :: ZCHIP                ! 
REAL, DIMENSION(SIZE(PPS))                         :: ZTAU_N               ! 
REAL, DIMENSION(SIZE(PPS))                         :: ZALBG                ! Effective ground albedo
REAL, DIMENSION(SIZE(PPS))                         :: ZSIGMA_F             ! LW transmission factor
REAL, DIMENSION(SIZE(PPS))                         :: ZSIGMA_FN            ! LW transmission factor - including buried (snow) 
!                                                                          ! vegetation effect
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_V_DTV        ! LW Jacobian: flux derrivative d LWnet_v/dTv [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_V_DTG        ! LW Jacobian: flux derrivative d LWnet_v/dTg [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_V_DTN        ! LW Jacobian: flux derrivative d LWnet_v/dTn [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_G_DTV        ! LW Jacobian: flux derrivative d LWnet_g/dTv [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_G_DTG        ! LW Jacobian: flux derrivative d LWnet_g/dTg [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_G_DTN        ! LW Jacobian: flux derrivative d LWnet_g/dTn [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_N_DTV        ! LW Jacobian: flux derrivative d LWnet_n/dTv [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_N_DTG        ! LW Jacobian: flux derrivative d LWnet_n/dTg [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_N_DTN        ! LW Jacobian: flux derrivative d LWnet_n/dTn [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZWRMAX               ! maximum canopy water equivalent interception capacity  [kg/m2]
REAL, DIMENSION(SIZE(PPS))                         :: ZLAIN                ! reduced (by burying by ground-based snow) LAI (m2/m2)
REAL, DIMENSION(SIZE(PPS))                         :: ZRS                  ! stomatal resistance (s/m)
REAL, DIMENSION(SIZE(PPS))                         :: ZRSN                 ! stomatal resistance of non-snow-buried canopy (s/m)
!                                                                          ! Etv=>0 as F2=>0 (-)  
REAL, DIMENSION(SIZE(PPS))                         :: ZWRVNMAX             ! maximum snow water equivalent interception capacity (kg/m2)
REAL, DIMENSION(SIZE(PPS))                         :: ZPSNCV               ! intercepted canopy snow fraction (-) NOTE! Not the same as the
!                                                                          ! ground-based snowpack
REAL, DIMENSION(SIZE(PPS))                         :: ZMELTVN              ! intercepted canopy snow net freeze/melt rate (kg/m2/s)
!                                                                          ! (if it is < 0, this signifies freezing)
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMA_TA            ! linear transform energy budget coefficient for Ta
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMB_TA            ! linear transform energy budget coefficient for Ta
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMA_TC            ! linear transform energy budget coefficient for Tc
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMB_TC            ! linear transform energy budget coefficient for Tc
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMA_TN            ! linear transform energy budget coefficient for Tn
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMB_TN            ! linear transform energy budget coefficient for Tn
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMA_TG            ! linear transform energy budget coefficient for Tg
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMB_TG            ! linear transform energy budget coefficient for Tg
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMA_TV            ! linear transform energy budget coefficient for Tv
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMB_TV            ! linear transform energy budget coefficient for Tv
REAL, DIMENSION(SIZE(PPS))                         :: ZPET_A_COEF          ! atmospheric coupling coefficient: Ta
REAL, DIMENSION(SIZE(PPS))                         :: ZPET_B_COEF          ! atmospheric coupling coefficient: Ta
REAL, DIMENSION(SIZE(PPS))                         :: ZKVN                 ! snow interception efficiency
REAL, DIMENSION(SIZE(PPS))                         :: ZVELC                ! wind speed at the top of the canopy (m/s)
REAL, DIMENSION(SIZE(PPS))                         :: ZDELTA               ! fraction of the foliage
!                                                                          ! covered with intercepted water (-)
REAL, DIMENSION(SIZE(PPS))                         :: ZHUGI                ! humidity over frozen bare ground (-)
REAL, DIMENSION(SIZE(PPS))                         :: ZHVN                 ! Halstead coefficient vegetation canopy above snow (-) 
REAL, DIMENSION(SIZE(PPS))                         :: ZHVG                 ! Halstead coefficient vegetation canopy above snow-free ground (-) 
REAL, DIMENSION(SIZE(PPS))                         :: ZLEG_DELTA           ! soil evaporation delta fn (-) 
REAL, DIMENSION(SIZE(PPS))                         :: ZLEGI_DELTA          ! soil sublimation delta fn (-) 
REAL, DIMENSION(SIZE(PPS))                         :: ZHSGL                ! surface halstead cofficient for bare soil (-)
REAL, DIMENSION(SIZE(PPS))                         :: ZHSGF                ! surface halstead cofficient for bare soil ice  (-)
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_C_A            ! turb transfer coef between vegetation canopy air and atmosphere (kg/m2/s) 
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_N_A            ! ...between the snow on the ground and atmosphere    (kg/m2/s) 
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_G_C            ! ...between snow-free ground and canopy air     (kg/m2/s)    
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_N_C            ! ...between snow on the ground and canopy air   (kg/m2/s)     
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_VG_C           ! ...between vegetation canopy over snow-free ground and canopy air   (kg/m2/s) 
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_VN_C           ! ...between vegetation canopy over the snow on the ground and canopy air  (kg/m2/s)  
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_V_C            ! ...between vegetation canopy and canopy air  (kg/m2/s)               
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_MOM            ! Effective drag coefficient for momentum [kg/(m2 s)]    
REAL, DIMENSION(SIZE(PPS))                         :: ZQSATG               ! saturation specific humidity for PTG (ground surface: kg kg-1)    
REAL, DIMENSION(SIZE(PPS))                         :: ZQSATV               ! saturation specific humidity for PTV (vegetation canopy: kg kg-1) 
REAL, DIMENSION(SIZE(PPS))                         :: ZQSATC               ! saturation specific humidity for PTC (canopy air: kg kg-1)      
REAL, DIMENSION(SIZE(PPS))                         :: ZQSATN               ! saturation specific humidity for PSNOWTEMP (snow surface: kg kg-1) 
REAL, DIMENSION(SIZE(PPS))                         :: ZDELTAVK             ! canopy interception capacity fraction including K-factor (-)  
REAL, DIMENSION(SIZE(PPS))                         :: ZCHEATV              ! Vegetation canopy *effective surface* heat capacity    (J m-2 K-1)
REAL, DIMENSION(SIZE(PPS))                         :: ZCHEATG              ! Understory-ground *effective surface* heat capacity    (J m-2 K-1)
REAL, DIMENSION(SIZE(PPS))                         :: ZCHEATN              ! Ground-based snow *effective surface* heat capacity    (J m-2 K-1)
REAL, DIMENSION(SIZE(PPS))                         :: ZHVGS                ! Dimensionless pseudo humidity factor for computing 
!                                                                          !  vapor fluxes from the non-buried part of the canopy 
!                                                                          !  to the canopy air                                     (-)
REAL, DIMENSION(SIZE(PPS))                         :: ZHVNS                ! Dimensionless pseudo humidity factor for computing 
!                                                                          !  vapor fluxes from the partly-buried part of the canopy 
!                                                                          !  to the canopy air                                     (-)
REAL, DIMENSION(SIZE(PPS))                         :: ZDQSAT_G             ! saturation specific humidity derivative for understory (kg kg-1 K-1)
REAL, DIMENSION(SIZE(PPS))                         :: ZDQSAT_V             ! saturation specific humidity derivative for the  
!                                                                          !  vegetation canopy                                     (kg kg-1 K-1)
REAL, DIMENSION(SIZE(PPS))                         :: ZDQSATI_N            ! saturation specific humidity derivative over ice for 
!                                                                          !  the ground-based snowpack                             (kg kg-1 K-1)
REAL, DIMENSION(SIZE(PPS))                         :: ZDELTAT_G            ! Time change in soil surface temperature                (K)
REAL, DIMENSION(SIZE(PPS))                         :: ZDELTAT_V            ! Time change in vegetation canopy temperature           (K)
REAL, DIMENSION(SIZE(PPS))                         :: ZDELTAT_N            ! Time change in snowpack surface temperature            (K)
REAL, DIMENSION(SIZE(PPS))                         :: ZRNET_V              ! Net vegetation canopy radiation                        (W/m2)
REAL, DIMENSION(SIZE(PPS))                         :: ZRNET_G              ! Net understory-ground radiation                        (W/m2)
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_C_A_F          ! Exchange coefficient between the snow on the ground and 
!                                                                          !  atmosphere modified by a partially to fully buried 
!                                                                          !  vegetation canopy                                     [kg/(m2 s)]
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_N_A_F          ! Exchange coefficient between vegetation canopy air and 
!                                                                          !  atmosphere modified by a partially to fully buried 
!                                                                          !  vegetation canopy                                     [kg/(m2 s)]
REAL, DIMENSION(SIZE(PPS))                         :: ZEVAP_C_A            ! Total canopy evapotranspiration and sublimation
!                                                                          !  of intercepted snow                                    (kg/m2/s)
REAL, DIMENSION(SIZE(PPS))                         :: ZEVAP_N_A            ! Vapor flux from the ground-based snowpack (part burying 
!                                                                          !  the canopy vegetation) to the atmosphere              (kg/m2/s)
REAL, DIMENSION(SIZE(PPS))                         :: ZH_N_A               ! Sensible heat flux from the ground-based snowpack (part 
!                                                                          !  burying the canopy vegetation) to the atmosphere      (W/m2)
REAL, DIMENSION(SIZE(PPS))                         :: ZVEGFACT             ! Fraction of canopy vegetation possibly receiving 
!                                                                          !  rainfall                                              (-)
REAL, DIMENSION(SIZE(PPS))                         :: ZRRSFC               ! The sum of all non-intercepted rain and canopy drip    (kg/m2/s)
REAL, DIMENSION(SIZE(PPS))                         :: ZLES3L               ! latent heat flux - sublimation of ice from the ground 
!                                                                          !  based snowpack (W/m2)
REAL, DIMENSION(SIZE(PPS))                         :: ZLEL3L               ! latent heat flux - evaporation of liquid water from the 
!                                                                          !  ground based snowpack (W/m2))
REAL, DIMENSION(SIZE(PPS))                         :: ZEVAP3L              ! total mass loss via evap & sublm from the ground based snowpack (kg/m2/s)
REAL, DIMENSION(SIZE(PPS))                         :: ZUSTAR2_IC           ! friction velocity (possibly implicitly coupled) (m/s)
REAL, DIMENSION(SIZE(PPS))                         :: ZTA_IC               ! atmospheric temperature (possibly implicitly coupled) (m/s)
REAL, DIMENSION(SIZE(PPS))                         :: ZQA_IC               ! atmospheric specific humidity (possibly implicitly coupled) (m/s)
REAL, DIMENSION(SIZE(PPS))                         :: ZSWUP                ! net upwelling shortwave radiation [W/m2]
REAL, DIMENSION(SIZE(PPS))                         :: ZLWUP                ! net upwelling longwave radiation [W/m2]
REAL, DIMENSION(SIZE(PPS))                         :: ZUSTAR2SNOW          ! snow fraciton velocity squared (m2/s2)
REAL, DIMENSION(SIZE(PPS))                         :: ZTA                  ! lowest level atmospheric temperature update estimate (K)
REAL, DIMENSION(SIZE(PPS))                         :: ZQA                  ! lowest level atmospheric spec. humidity update estimate (K)
REAL, DIMENSION(SIZE(PPS))                         :: ZVMOD                ! lowest level atmospheric wind speed update estimate (K)
REAL, DIMENSION(SIZE(PPS))                         :: ZRR                  ! combined rain rate (above canopy) and irrigation need (kg/m2/s)
REAL, DIMENSION(SIZE(PPS))                         :: ZFLSN_COR            ! snow/soil-biomass correction flux (W/m2) (not MEB)
!
!
! Working sums for flux averaging over MEB time split
!
REAL, DIMENSION(SIZE(PPS))   :: ZH_SUM, ZH_C_A_SUM, ZH_N_A_SUM, ZH_V_C_SUM, ZH_G_C_SUM, &
                                ZH_N_C_SUM, ZHSNOW_SUM, ZHPSNOW_SUM, ZGFLUXSNOW_SUM
REAL, DIMENSION(SIZE(PPS))   :: ZHU_AGG_SUM, ZAC_AGG_SUM

REAL, DIMENSION(SIZE(PPS))   :: ZLE_SUM, ZLE_C_A_SUM, ZLE_V_C_SUM, ZLE_G_C_SUM,           &
                                ZLE_N_C_SUM, ZLETR_V_C_SUM, ZLEG_SUM,                     &
                                ZLER_V_C_SUM, ZLE_FLOOD_SUM, ZLEI_FLOOD_SUM,              &
                                ZLES_V_C_SUM, ZLETR_SUM, ZLER_SUM, ZLEV_SUM,              &
                                ZLEI_SUM, ZLES3L_SUM, ZLEL3L_SUM, ZEVAP3L_SUM,            &
                                ZUSTAR2_SUM, ZUSTAR2SNOW_SUM, ZCDSNOW_SUM,                &
                                ZCHSNOW_SUM, ZRISNOW_SUM

REAL, DIMENSION(SIZE(PPS))   :: ZGRNDFLUX_SUM, ZRESTORE_SUM

REAL, DIMENSION(SIZE(PPS))   :: ZSWNET_V_SUM, ZSWNET_G_SUM, ZSWNET_N_SUM, ZLWNET_V_SUM, &
                                ZLWNET_G_SUM, ZLWNET_N_SUM, ZEMIST_SUM, ZSWUP_SUM,      &
                                ZLWUP_SUM
REAL, DIMENSION(SIZE(PPS))   :: ZDELHEATG_SFC_SUM, ZDELHEATV_SFC_SUM, ZDELHEATG_SUM
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*      1.0    Preliminaries
!              -------------
!
IF (LHOOK) CALL DR_HOOK('ISBA_MEB',0,ZHOOK_HANDLE)
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      2.0    Preliminaries for energy and radiation budget
!              ---------------------------------------------
!
CALL PREPS_FOR_MEB_EBUD_RAD(PPS,                                     &
        PLAI,PSNOWRHO,PSNOWSWE,PSNOWHEAT,                            &
        PSNOWTEMP,PSNOWDZ,ZSNOWCOND,ZSNOWHCAP,PEMISNOW,              &
        ZSIGMA_F,ZCHIP)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      3.0    Shortwave radiative transfer
!              ----------------------------  
!
! Calculate snow albedo: split into spectral bands:
!
CALL SNOWALB_SPECTRAL_BANDS_MEB(PVEGTYPE,PSNOWALB,PSNOWRHO(:,1),PSNOWAGE(:,1),PPS, &
                                PSNOWDZ(:,1), PZENITH,                             &
                                PSNOWALBVIS,PSNOWALBNIR,PSNOWALBFIR,ZTAU_N)
!
!
! NOTE, currently MEB only uses 2 of 3 potential snow albedo spectral bands
!
CALL ISBA_SWNET_MEB(PLAI,PZF_TALLVEG,                                        &
        PSNOWALBVIS,PSNOWALBNIR,                                             &
        PALBNIR_TVEG, PALBVIS_TVEG,PALBNIR_TSOIL, PALBVIS_TSOIL, PFALB,      &
        PFF,PPSN,PPALPHAN,ZTAU_N,PZENITH,PSCA_SW,PSW_RAD,                    &
        ZSWUP,PSWNET_N,PSWNET_V,PSWNET_G,PSWNET_NS,                          &
        PALBT,ZALBG,PSWDOWN_GN                                               )
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      4.0    Longwave radiative transfer
!              ---------------------------  
!
CALL ISBA_LWNET_MEB(PLAI,PPSN,PPALPHAN,                                 &
        PEMISNOW,PFEMIS,PFF,                                            &
        PTV,PTG(:,1),PSNOWTEMP(:,1),                                    &
        PLW_RAD,PLWNET_N,PLWNET_V,PLWNET_G,                             &
        ZDLWNET_V_DTV,ZDLWNET_V_DTG,ZDLWNET_V_DTN,                      &
        ZDLWNET_G_DTV,ZDLWNET_G_DTG,ZDLWNET_G_DTN,                      &
        ZDLWNET_N_DTV,ZDLWNET_N_DTG,ZDLWNET_N_DTN,                      &
        ZSIGMA_F,ZSIGMA_FN,PLWDOWN_GN                                   )
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      5.0    Fraction of leaves occupied by intercepted water
!              ------------------------------------------------
!
! Vegetation canopy:
!
! First, compute an effective veg fraction: it can only be < unity if vegetation is buried by snowpack...
!
ZWORK(:) = (1.0 - PPSN(:) + PPSN(:)*(1.0 - PPALPHAN(:))) 
! 
CALL WET_LEAVES_FRAC(PWR, ZWORK, PWRMAX_CF, PZ0_MEBV, PLAI, ZWRMAX, ZDELTA) 
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      6.0    Plant stress, stomatal resistance and, possibly, CO2 assimilation
!              --------------------------------------------------------------------
!
! MEB-NOTE here assumed HPHOTO=='DEF' for now, Ags options to be added later 
!
!
! Canopy vegetation (no snow, or snow below the main part of the canopy):
!
CALL VEG(PSW_RAD, PTC, PQC, PPS, PRGL, PLAI, PRSMIN, PGAMMA, PF2, ZRS)
!
! Possibly snow-buried canopy vegetation:
!
ZLAIN(:) = PLAI(:)*(1.0-PPALPHAN(:)) ! LAI of possibly buried vegetation canopy
!
CALL VEG(PSW_RAD, PTC, PQC, PPS, PRGL, ZLAIN, PRSMIN, PGAMMA, PF2, ZRSN)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      6.0    Canopy snow (intercepted) needed diagnostics:
!              ---------------------------------------------
!
CALL SNOW_LEAVES_FRAC_MEB(PPSN,PPALPHAN,PWRVN,PTV,ZCHIP,PLAI,        &
                             ZWRVNMAX,ZPSNCV,ZMELTVN)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      7.0    Aerodynamic drag and heat/mass transfer/fluxes 
!              and energy budget solution
!              ----------------------------------------------
!
! NOTE, this assumes thermodynamic variable herein is potential T

ZPET_A_COEF(:)  = -PPET_A_COEF(:)*XCPD 
ZPET_B_COEF(:)  =  PPET_B_COEF(:)*XCPD
ZTHRMA_TA(:)    =  XCPD/PEXNA(:)
ZTHRMB_TA(:)    =  0.0
ZWORK(:)        =  XCPD/PEXNS(:)
ZTHRMA_TC(:)    =  ZWORK(:)
ZTHRMB_TC(:)    =  0.0
ZTHRMA_TN(:)    =  ZWORK(:)
ZTHRMB_TN(:)    =  0.0
ZTHRMA_TG(:)    =  ZWORK(:)
ZTHRMB_TG(:)    =  0.0
ZTHRMA_TV(:)    =  ZWORK(:)
ZTHRMB_TV(:)    =  0.0
!
!
! Possibly split time step if large: 
! Although the energy budget is fully implicit, a very small canopy heat capacity 
! (and neglect of canopy air space heat capacity) can possibly lead to
! numerical shocks, especially during transition periods between stable and unstable 
! regimes. Thus, for relatively large steps, a simple time split scheme is activated.
! Note that soil moisture is held constant, while turbulent exchange coefficients are updated during the split.
! Also, experience shows that splitting at least once for moderately sized time steps 
! is quite effective in removing any lingering small but possible oscillations.
! Finally, for *very* small time steps (such as those for high res runs), no split is performed.
! Fluxes are averaged over the time split for conservation.
!
JTSPLIT_EB      = 1 + INT(PTSTEP/ZTSTEP_EB)  ! number of split-time steps
ZTSTEP          = PTSTEP/JTSPLIT_EB          ! split time step...for relatively small time steps, no split
!
! initialize time split sums for fluxes:
!
CALL INIT_SUM_FLUXES_MEB_TSPLIT 
!
!
! Note, when implicitly coupled to the atmosphere, these
! 3 variables will evolve during the split...we used updated
! values for turbulent exchange computations (drag_meb). 
! NOTE...when explicit coupling used, these 3 variables do NOT vary
! during the split.
!
ZVMOD(:) = PVMOD(:)
ZTA(:)   = PTA(:)
ZQA(:)   = PQA(:)
!
!
LOOP_TIME_SPLIT_EB: DO JDT=1,JTSPLIT_EB
!*      7.1    Aerodynamic drag and heat transfer coefficients
!              -----------------------------------------------
!
   CALL DRAG_MEB(OFORC_MEASURE,                                         &
              PTG(:,1), PTC, PTV, PSNOWTEMP(:,1), ZTA, PQC, ZQA, ZVMOD, &
              PWG(:,1), PWGI(:,1), PWSAT(:,1), PWFC(:,1),               &
              PEXNS, PEXNA, PPS,                                        &
              PRR, PSR, PRHOA, PZ0G_WITHOUT_SNOW,                       &
              PZ0_MEBV, PZ0H_MEBV, PZ0EFF_MEBV,                         &
              PZ0_MEBN, PZ0H_MEBN, PZ0EFF_MEBN,                         &
              PZ0_WITH_SNOW, PZ0H_WITH_SNOW, PZ0EFF,                    &
              PSNOWSWE(:,1),                                            &
              PWR, ZCHIP, ZTSTEP, ZRS, ZRSN,                            &
              PPSN, PPALPHAN, PZREF, PUREF, PH_VEG, PDIRCOSZW,          &
              ZPSNCV, ZDELTA, PLAI, PGNDLITTER,                         &
              PCH, PCD, PCDN, PRI, PRESA, ZVELC,                        &
              PCDSNOW, PCHSNOW, PRISNOW, ZUSTAR2SNOW,                   &
              PHUG, ZHUGI, PHV, ZHVG, ZHVN, PHU, PQS, PRS,              &
              ZLEG_DELTA, ZLEGI_DELTA, ZHSGL, ZHSGF,                    &
              ZFLXC_C_A, ZFLXC_N_A, ZFLXC_G_C, ZFLXC_N_C,               &    
              ZFLXC_VG_C, ZFLXC_VN_C, ZFLXC_MOM,                        &
              ZQSATG, ZQSATV, ZQSATC, ZQSATN, ZDELTAVK                  )
!
   ZKVN(:) = SNOW_INTERCEPT_EFF(ZCHIP,ZVELC,ZWRVNMAX)

!*      7.2    Resolution of the surface energy budgets
!              ----------------------------------------
!
   CALL E_BUDGET_MEB(HISBA,HCPSURF,ZTSTEP,                                             &
              PPS,PCG,PCT,PCV,PWRVN,PWR,                                               &
              PTDEEP_A,PTDEEP_B,PD_G,PSOILCONDZ,PSOILHCAPZ,                            &
              PSNOWDZ,ZSNOWCOND,ZSNOWHCAP,                                             &
              PSWNET_V,PSWNET_G,PSWNET_NS,                                             &
              PLWNET_V,PLWNET_G,PLWNET_N,                                              &
              ZDLWNET_V_DTV,ZDLWNET_V_DTG,ZDLWNET_V_DTN,                               &
              ZDLWNET_G_DTV,ZDLWNET_G_DTG,ZDLWNET_G_DTN,                               &
              ZDLWNET_N_DTV,ZDLWNET_N_DTG,ZDLWNET_N_DTN,                               &
              PPEW_A_COEF,PPEW_B_COEF,ZPET_A_COEF,PPEQ_A_COEF,ZPET_B_COEF,PPEQ_B_COEF, &
              ZTHRMA_TA,ZTHRMB_TA,ZTHRMA_TC,ZTHRMB_TC,                                 &
              ZTHRMA_TG,ZTHRMB_TG,ZTHRMA_TV,ZTHRMB_TV,ZTHRMA_TN,ZTHRMB_TN,             &
              ZQSATG,ZQSATV,ZQSATN,                                                    &
              PFF,PFFROZEN,PPSN,PPALPHAN,ZPSNCV,                                       &
              ZCHEATV,ZCHEATG,ZCHEATN,                                                 &
              ZLEG_DELTA,ZLEGI_DELTA,PHUG,ZHUGI,ZHVG,ZHVN,PFROZEN1,                    &
              ZFLXC_C_A,ZFLXC_G_C,ZFLXC_VG_C,ZFLXC_VN_C,ZFLXC_N_C,ZFLXC_N_A,           &
              ZFLXC_MOM,                                                               &
              PTG,PTV,PSNOWTEMP,                                                       &
              ZFLXC_V_C,ZHVGS,ZHVNS,                                                   &
              ZDQSAT_G,ZDQSAT_V,ZDQSATI_N,                                             &
              PTC,PQC,ZTA_IC,ZQA_IC,ZUSTAR2_IC,ZVMOD,                                  &
              ZDELTAT_G,ZDELTAT_V,ZDELTAT_N,PGRNDFLUX,PCPS,PLVTT,PLSTT,                &
              PHPSNOW,PMELTADV,PRESTORE,PDEEP_FLUX,                                    &
              PDELHEATV_SFC,PDELHEATG_SFC,PDELHEATG                                    )
!
!*      7.3    Energy and momentum fluxes and radiative temperature and emissivity
!              -------------------------------------------------------------------
!
   CALL ISBA_FLUXES_MEB(PRHOA,                                                         &
              ZSIGMA_F,ZSIGMA_FN,PEMISNOW,                                             &
              ZRNET_V,ZRNET_G,PRNSNOW,                                                 & 
              PSWNET_V,PSWNET_G,PSWNET_N,                                              &
              PLWNET_V,PLWNET_G,PLWNET_N,                                              &
              ZDLWNET_V_DTV,ZDLWNET_V_DTG,ZDLWNET_V_DTN,                               &
              ZDLWNET_G_DTV,ZDLWNET_G_DTG,ZDLWNET_G_DTN,                               &
              ZDLWNET_N_DTV,ZDLWNET_N_DTG,ZDLWNET_N_DTN,                               &
              ZTHRMA_TA,ZTHRMB_TA,ZTHRMA_TC,ZTHRMB_TC,                                 &
              ZTHRMA_TG,ZTHRMB_TG,ZTHRMA_TV,ZTHRMB_TV,ZTHRMA_TN,ZTHRMB_TN,             &
              ZQSATG,ZQSATV,ZQSATN,                                                    &
              PFF,PPSN,PPALPHAN,ZPSNCV,PFROZEN1,PFFROZEN,                              &
              ZLEG_DELTA,ZLEGI_DELTA,PHUG,ZHUGI,ZHVG,ZHVN,                             &
              ZFLXC_C_A,ZFLXC_G_C,ZFLXC_VG_C,ZFLXC_VN_C,ZFLXC_N_C,ZFLXC_N_A,           &
              ZFLXC_MOM,ZFLXC_V_C,ZHVGS,ZHVNS,                                         &
              PTG,PTV,PSNOWTEMP,                                                       &
              ZDQSAT_G,ZDQSAT_V,ZDQSATI_N,                                             &
              PTC,PQC,ZTA_IC,ZQA_IC,                                                   &
              ZDELTAVK,                                                                &
              ZDELTAT_G,ZDELTAT_V,ZDELTAT_N,                                           &
              ZSWUP,PSW_RAD,PLW_RAD,                                                   &
              PRN,ZLWUP,                                                               &
              PH_C_A,PH_V_C,PH_G_C,PH_N_C,ZH_N_A,PHSNOW,PH,                            &
              PLE_C_A,PLE_V_C,PLE_G_C,PLE_N_C,                                         &
              ZEVAP_C_A,PLEV_V_C,PEVAP_G_C,PEVAP_N_C,ZEVAP_N_A,                        &
              PEVAP,PSUBL,PLETR_V_C,PLER_V_C,PLEG,PLEGI,                               &
              PLE_FLOOD,PLEI_FLOOD,ZLES3L,ZLEL3L,                                      &
              ZEVAP3L,PLES_V_C,PLETR,PLER,PLEV,PLE,PLEI,                               &
              PTS_RAD,PEMIST                                                           )
!
! Compute aggregated coefficients for evaporation
! Sum(LEC+LES+LEL) = ACagg * Lv * RHOA * (HUagg.Qsat - Qa)
!
   ZFLXC_C_A_F(:) = ZFLXC_C_A(:)*(1.0-PPSN(:)*PPALPHAN(:))
   ZFLXC_N_A_F(:) = ZFLXC_N_A(:)*     PPSN(:)*PPALPHAN(:)

   PHU_AGG(:)     = (ZFLXC_C_A_F(:)*PQC(:)    + ZFLXC_N_A_F(:)*ZQSATN(:))/   &
                    (ZFLXC_C_A_F(:)*ZQSATC(:) + ZFLXC_N_A_F(:)*ZQSATN(:))

   PAC_AGG(:)     = ZFLXC_C_A_F(:) + ZFLXC_N_A_F(:) ! kg/m2/s
!
! Sum fluxes over time split:

   CALL SUM_FLUXES_MEB_TSPLIT  

ENDDO LOOP_TIME_SPLIT_EB
!
CALL AVG_FLUXES_MEB_TSPLIT     ! average fluxes over time split
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*     8.0    Snow explicit canopy loading/interception 
!             ------------------------------------------
!
CALL SNOW_LOAD_MEB(PTSTEP,PSR,PTV,ZWRVNMAX,ZKVN,ZCHEATV,PLER_V_C,PLES_V_C,ZMELTVN, &
     ZVELC,PMELTCV,PFRZCV,PSR_GN,PWR,PWRVN,PSUBVCOR)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*     9.0    Snow explicit canopy loading/interception 
!             ------------------------------------------
!
ZRR(:)         = PRR(:)
PIRRIG_FLUX(:) = 0.0
!
!* add irrigation over vegetation to liquid precipitation (rr)
!  Water "need" treated as if sprayed from above (over vegetation and soil):
!
IF (SIZE(OIRRIGATE)>0) THEN
   WHERE (OIRRIGATE(:) .AND. PIRRIG(:)>0. .AND. PIRRIG(:) /= XUNDEF .AND. (PF2(:)<PTHRESHOLD(:)) )
      PIRRIG_FLUX(:) = PWATSUP(:) / XDAY           
      ZRR        (:) = PRR(:) + PWATSUP(:)/XDAY
      OIRRIDAY   (:) = .TRUE.           
   END WHERE
ENDIF
!
! Call canopy interception...here because meltwater should be allowed to fall
! on understory snowpack
!
! Fraction of canopy vegetation possibly receiving rainfall/irrigation
!
ZVEGFACT(:) = ZSIGMA_F(:)*(1.0-PPALPHAN(:)*PPSN(:)) 
!
! The sum of all non-intercepted rain and drip is "ZRRSFC" (kg/m2/s):
! this is then partitioned by snow scheme into part falling on
! snowpack and part falling onto snow-free understory.
!
!
CALL HYDRO_VEG(HRAIN, PTSTEP, PMUF,                      &
        ZRR, PLEV_V_C, PLETR_V_C, ZVEGFACT, ZPSNCV,      &
        PWR, ZWRMAX, ZRRSFC, PDRIP, PRRVEG               )
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      10.0    Explicit snow scheme (MEB: impose surface fluxes as upper BC)
!              ----------------------------------------------------------------
!
CALL SNOW3L_ISBA(HISBA, HSNOW_ISBA, HSNOWRES, OMEB, OGLACIER, HIMPLICIT_WIND,          &
           TPTIME, PTSTEP, PVEGTYPE,                                                   &
           PSNOWSWE, PSNOWHEAT, PSNOWRHO, PSNOWALB,                                    &
           PSNOWGRAN1, PSNOWGRAN2, PSNOWHIST,PSNOWAGE,                                 &
           PTG(:,1), PCG, PCT, PSOILCONDZ(:,1),                                        &
           PPS, PTA, PSW_RAD, PQA, PVMOD, PLW_RAD, ZRRSFC, PSR_GN,                     &
           PRHOA, PUREF, PEXNS, PEXNA, PDIRCOSZW,                                      &
           PZREF, PZ0_WITH_SNOW, PZ0EFF, PZ0H_WITH_SNOW, ZALBG, PD_G(:,1),             &
           PPEW_A_COEF, PPEW_B_COEF,                                                   &
           PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,                         &
           PSNOW_THRUFAL, PGRNDFLUX, PFLSN_COR, PRESTOREN, PEVAPCOR,                   &
           PSWNET_N, PSWNET_NS, PLWNET_N,                                              &
           PRNSNOW, PHSNOW, PGFLUXSNOW, PHPSNOW, ZLES3L, ZLEL3L, ZEVAP3L,              &
           PSNDRIFT, PUSTARSNOW,                                                       & 
           PPSN, PSRSFC, PRRSFC, PSNOWSFCH, PDELHEATN, PDELHEATN_SFC,                  &
           PEMISNOW, PCDSNOW, PCHSNOW, PSNOWTEMP, PSNOWLIQ, PSNOWDZ,                   &
           PSNOWHMASS, PRISNOW, PZENITH, PDELHEATG, PDELHEATG_SFC, PLAT, PLON, PQSNOW, &
           OSNOWDRIFT, OSNOWDRIFT_SUBLIM, OSNOW_ABS_ZENITH,                            &
           HSNOWMETAMO, HSNOWRAD                                                       )  
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!
IF (LHOOK) CALL DR_HOOK('ISBA_MEB',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
CONTAINS
!
!===============================================================================
!
SUBROUTINE INIT_SUM_FLUXES_MEB_TSPLIT 
!
IMPLICIT NONE
!
!
!*      0.2    declarations of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ISBA_MEB:INIT_SUM_FLUXES_MEB_TSPLIT ',0,ZHOOK_HANDLE)
!
! sensible heat fluxes:
!
ZH_SUM(:)        = 0.0
ZH_C_A_SUM(:)    = 0.0
ZH_N_A_SUM(:)    = 0.0
ZH_V_C_SUM(:)    = 0.0
ZH_G_C_SUM(:)    = 0.0
ZH_N_C_SUM(:)    = 0.0
ZHSNOW_SUM(:)    = 0.0
!
! latent heat/water vapor fluxes:
!
ZLE_SUM(:)       = 0.0
ZLE_C_A_SUM(:)   = 0.0
ZLE_V_C_SUM(:)   = 0.0
ZLE_G_C_SUM(:)   = 0.0
ZLE_N_C_SUM(:)   = 0.0
ZLETR_V_C_SUM(:) = 0.0
ZLEG_SUM(:)      = 0.0
ZLER_V_C_SUM(:)  = 0.0
ZLE_FLOOD_SUM(:) = 0.0
ZLEI_FLOOD_SUM(:)= 0.0
ZLES_V_C_SUM(:)  = 0.0
ZLETR_SUM(:)     = 0.0
ZLER_SUM(:)      = 0.0
ZLEV_SUM(:)      = 0.0
ZLEI_SUM(:)      = 0.0
ZLES3L_SUM(:)    = 0.0
ZLEL3L_SUM(:)    = 0.0
ZEVAP3L_SUM(:)   = 0.0
!
ZHU_AGG_SUM(:)   = 0.0
ZAC_AGG_SUM(:)   = 0.0
!
! momentum/turb:
!
ZUSTAR2_SUM(:)     = 0.0
ZUSTAR2SNOW_SUM(:) = 0.
ZCDSNOW_SUM(:)     = 0.
ZCHSNOW_SUM(:)     = 0.
ZRISNOW_SUM(:)     = 0.
!
! surface interfacial/sub-surface heat fluxes:
!
ZGRNDFLUX_SUM(:) = 0.0
ZRESTORE_SUM(:)  = 0.0
ZGFLUXSNOW_SUM(:)= 0.0
ZHPSNOW_SUM(:)   = 0.0
!
! radiative fluxes:
!
ZSWNET_V_SUM(:)  = 0.0
ZSWNET_G_SUM(:)  = 0.0
ZSWNET_N_SUM(:)  = 0.0
ZLWNET_V_SUM(:)  = 0.0
ZLWNET_G_SUM(:)  = 0.0
ZLWNET_N_SUM(:)  = 0.0
ZEMIST_SUM(:)    = 0.0
ZSWUP_SUM(:)     = 0.0
ZLWUP_SUM(:)     = 0.0
!
ZDELHEATV_SFC_SUM(:) = 0.0 
ZDELHEATG_SFC_SUM(:) = 0.0 
ZDELHEATG_SUM(:)     = 0.0 
!
IF (LHOOK) CALL DR_HOOK('ISBA_MEB:INIT_SUM_FLUXES_MEB_TSPLIT ',1,ZHOOK_HANDLE)
!
END SUBROUTINE INIT_SUM_FLUXES_MEB_TSPLIT
!
!===============================================================================
!
SUBROUTINE SUM_FLUXES_MEB_TSPLIT 
!
IMPLICIT NONE
!
!*      0.2    declarations of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ISBA_MEB:SUM_FLUXES_MEB_TSPLIT ',0,ZHOOK_HANDLE)
!
! Sum fluxes over MEB TIME SPLIT:
!
! sensible heat fluxes:

ZH_SUM(:)        = ZH_SUM(:)        + PH(:)
ZH_C_A_SUM(:)    = ZH_C_A_SUM(:)    + PH_C_A(:)
ZH_N_A_SUM(:)    = ZH_N_A_SUM(:)    + ZH_N_A(:)
ZH_V_C_SUM(:)    = ZH_V_C_SUM(:)    + PH_V_C(:)
ZH_G_C_SUM(:)    = ZH_G_C_SUM(:)    + PH_G_C(:)
ZH_N_C_SUM(:)    = ZH_N_C_SUM(:)    + PH_N_C(:)
ZHSNOW_SUM(:)    = ZHSNOW_SUM(:)    + PHSNOW(:)
!
! latent heat/water vapor fluxes:
!
ZLE_SUM(:)       = ZLE_SUM(:)       + PLE(:)
ZLE_C_A_SUM(:)   = ZLE_C_A_SUM(:)   + PLE_C_A(:)
ZLE_V_C_SUM(:)   = ZLE_V_C_SUM(:)   + PLE_V_C(:) 
ZLE_G_C_SUM(:)   = ZLE_G_C_SUM(:)   + PLE_G_C(:) 
ZLE_N_C_SUM(:)   = ZLE_N_C_SUM(:)   + PLE_N_C(:) 
ZLETR_V_C_SUM(:) = ZLETR_V_C_SUM(:) + PLETR_V_C(:) 
ZLEG_SUM(:)      = ZLEG_SUM(:)      + PLEG(:) 
ZLER_V_C_SUM(:)  = ZLER_V_C_SUM(:)  + PLER_V_C(:) 
ZLE_FLOOD_SUM(:) = ZLE_FLOOD_SUM(:) + PLE_FLOOD(:)
ZLEI_FLOOD_SUM(:)= ZLEI_FLOOD_SUM(:)+ PLEI_FLOOD(:) 
ZLES_V_C_SUM(:)  = ZLES_V_C_SUM(:)  + PLES_V_C(:)
ZLETR_SUM(:)     = ZLETR_SUM(:)     + PLETR(:) 
ZLER_SUM(:)      = ZLER_SUM(:)      + PLER(:)
ZLEV_SUM(:)      = ZLEV_SUM(:)      + PLEV(:) 
ZLEI_SUM(:)      = ZLEI_SUM(:)      + PLEI(:)
ZLES3L_SUM(:)    = ZLES3L_SUM(:)    + ZLES3L(:) 
ZLEL3L_SUM(:)    = ZLEL3L_SUM(:)    + ZLEL3L(:) 
ZEVAP3L_SUM(:)   = ZEVAP3L_SUM(:)   + ZEVAP3L(:) 
!
ZHU_AGG_SUM(:)   = ZHU_AGG_SUM(:)   + PHU_AGG(:)
ZAC_AGG_SUM(:)   = ZAC_AGG_SUM(:)   + PAC_AGG(:)
!
! momentum/turb:
!
ZUSTAR2_SUM(:)     = ZUSTAR2_SUM(:)     + ZUSTAR2_IC(:)
ZUSTAR2SNOW_SUM(:) = ZUSTAR2SNOW_SUM(:) + ZUSTAR2SNOW(:)
ZCDSNOW_SUM(:)     = ZCDSNOW_SUM(:)     + PCDSNOW(:)
ZCHSNOW_SUM(:)     = ZCHSNOW_SUM(:)     + PCHSNOW(:)
ZRISNOW_SUM(:)     = ZRISNOW_SUM(:)     + PRISNOW(:)
!
! surface interfacial/sub-surface heat fluxes:
!
ZGRNDFLUX_SUM(:) = ZGRNDFLUX_SUM(:) + PGRNDFLUX(:) 
ZRESTORE_SUM(:)  = ZRESTORE_SUM(:)  + PRESTORE(:) 
ZGFLUXSNOW_SUM(:)= ZGFLUXSNOW_SUM(:)+ PGFLUXSNOW(:) 
ZHPSNOW_SUM(:)   = ZHPSNOW_SUM(:)   + PHPSNOW(:)
!
! radiative fluxes:
!
ZSWNET_V_SUM(:)  = ZSWNET_V_SUM(:)  +   PSWNET_V(:)
ZSWNET_G_SUM(:)  = ZSWNET_G_SUM(:)  +   PSWNET_G(:) 
ZSWNET_N_SUM(:)  = ZSWNET_N_SUM(:)  +   PSWNET_N(:) 
ZLWNET_V_SUM(:)  = ZLWNET_V_SUM(:)  +   PLWNET_V(:) 
ZLWNET_G_SUM(:)  = ZLWNET_G_SUM(:)  +   PLWNET_G(:)
ZLWNET_N_SUM(:)  = ZLWNET_N_SUM(:)  +   PLWNET_N(:) 
ZEMIST_SUM(:)    = ZEMIST_SUM(:)    +   PEMIST(:) 
ZSWUP_SUM(:)     = ZSWUP_SUM(:)     +   ZSWUP(:)
ZLWUP_SUM(:)     = ZLWUP_SUM(:)     +   ZLWUP(:)
!
ZDELHEATV_SFC_SUM(:) = ZDELHEATV_SFC_SUM(:) +   PDELHEATV_SFC(:) 
ZDELHEATG_SFC_SUM(:) = ZDELHEATG_SFC_SUM(:) +   PDELHEATG_SFC(:) 
ZDELHEATG_SUM(:)     = ZDELHEATG_SUM(:)     +   PDELHEATG(:) 
!
IF (LHOOK) CALL DR_HOOK('ISBA_MEB:SUM_FLUXES_MEB_TSPLIT ',1,ZHOOK_HANDLE)
!
END SUBROUTINE SUM_FLUXES_MEB_TSPLIT
!
!===============================================================================
!
SUBROUTINE AVG_FLUXES_MEB_TSPLIT 
!
USE MODD_CSTS, ONLY : XSTEFAN
!
IMPLICIT NONE
!
!*      0.2    declarations of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ISBA_MEB:AVG_FLUXES_MEB_TSPLIT ',0,ZHOOK_HANDLE)
!
! Average fluxes over MEB TIME SPLIT:
!
! sensible heat fluxes:
!
PH(:)        = ZH_SUM(:)        /JTSPLIT_EB
PH_C_A(:)    = ZH_C_A_SUM(:)    /JTSPLIT_EB
ZH_N_A(:)    = ZH_N_A_SUM(:)    /JTSPLIT_EB
PH_V_C(:)    = ZH_V_C_SUM(:)    /JTSPLIT_EB
PH_G_C(:)    = ZH_G_C_SUM(:)    /JTSPLIT_EB
PH_N_C(:)    = ZH_N_C_SUM(:)    /JTSPLIT_EB
PHSNOW(:)    = ZHSNOW_SUM(:)    /JTSPLIT_EB
!
! latent heat/water vapor fluxes:
!
PLE(:)       = ZLE_SUM(:)       /JTSPLIT_EB
PLE_C_A(:)   = ZLE_C_A_SUM(:)   /JTSPLIT_EB
PLE_V_C(:)   = ZLE_V_C_SUM(:)   /JTSPLIT_EB
PLE_G_C(:)   = ZLE_G_C_SUM(:)   /JTSPLIT_EB
PLE_N_C(:)   = ZLE_N_C_SUM(:)   /JTSPLIT_EB
PLETR_V_C(:) = ZLETR_V_C_SUM(:) /JTSPLIT_EB
PLEG(:)      = ZLEG_SUM(:)      /JTSPLIT_EB
PLER_V_C(:)  = ZLER_V_C_SUM(:)  /JTSPLIT_EB
PLE_FLOOD(:) = ZLE_FLOOD_SUM(:) /JTSPLIT_EB
PLEI_FLOOD(:)= ZLEI_FLOOD_SUM(:)/JTSPLIT_EB
PLES_V_C(:)  = ZLES_V_C_SUM(:)  /JTSPLIT_EB
PLETR(:)     = ZLETR_SUM(:)     /JTSPLIT_EB
PLER(:)      = ZLER_SUM(:)      /JTSPLIT_EB
PLEV(:)      = ZLEV_SUM(:)      /JTSPLIT_EB
PLEI(:)      = ZLEI_SUM(:)      /JTSPLIT_EB
ZLES3L(:)    = ZLES3L_SUM(:)    /JTSPLIT_EB
ZLEL3L(:)    = ZLEL3L_SUM(:)    /JTSPLIT_EB
ZEVAP3L(:)   = ZEVAP3L_SUM(:)   /JTSPLIT_EB
!
PHU_AGG(:)   = ZHU_AGG_SUM(:)   /JTSPLIT_EB
PAC_AGG(:)   = ZAC_AGG_SUM(:)   /JTSPLIT_EB
!
! momentum/turb:
!
PUSTAR(:)     = SQRT( ZUSTAR2_SUM(:)    /JTSPLIT_EB )
PUSTARSNOW(:) = SQRT( ZUSTAR2SNOW_SUM(:)/JTSPLIT_EB )
PCDSNOW(:)    = ZCDSNOW_SUM(:)          /JTSPLIT_EB 
PCHSNOW(:)    = ZCHSNOW_SUM(:)          /JTSPLIT_EB 
PRISNOW(:)    = ZRISNOW_SUM(:)          /JTSPLIT_EB 
!
! surface interfacial/sub-surface heat fluxes:
!
PGRNDFLUX(:) = ZGRNDFLUX_SUM(:) /JTSPLIT_EB
PRESTORE(:)  = ZRESTORE_SUM(:)  /JTSPLIT_EB
PGFLUXSNOW(:)= ZGFLUXSNOW_SUM(:)/JTSPLIT_EB
PHPSNOW(:)   = ZHPSNOW_SUM(:)   /JTSPLIT_EB
!
! radiative fluxes:
!
PSWNET_V(:)  = ZSWNET_V_SUM(:)  /JTSPLIT_EB
PSWNET_G(:)  = ZSWNET_G_SUM(:)  /JTSPLIT_EB
PSWNET_N(:)  = ZSWNET_N_SUM(:)  /JTSPLIT_EB
PLWNET_V(:)  = ZLWNET_V_SUM(:)  /JTSPLIT_EB
PLWNET_G(:)  = ZLWNET_G_SUM(:)  /JTSPLIT_EB
PLWNET_N(:)  = ZLWNET_N_SUM(:)  /JTSPLIT_EB
PEMIST(:)    = ZEMIST_SUM(:)    /JTSPLIT_EB
ZSWUP(:)     = ZSWUP_SUM(:)     /JTSPLIT_EB
ZLWUP(:)     = ZLWUP_SUM(:)     /JTSPLIT_EB
!
PDELHEATV_SFC(:) = ZDELHEATV_SFC_SUM(:) /JTSPLIT_EB 
PDELHEATG_SFC(:) = ZDELHEATG_SFC_SUM(:) /JTSPLIT_EB 
PDELHEATG(:)     = ZDELHEATG_SUM(:)     /JTSPLIT_EB 
!
! Additional diagnostics depending on AVG quantities:
!
PTS_RAD(:)   = ((ZLWUP(:) - PLW_RAD(:)*(1.0-PEMIST(:)))/(XSTEFAN*PEMIST(:)))**0.25
!
ZRNET_V(:)   = PSWNET_V(:) + PLWNET_V(:)
!
ZRNET_G(:)   = PSWNET_G(:) + PLWNET_G(:)
!
PRNSNOW(:)   = PSWNET_N(:) + PLWNET_N(:)
!
PRN(:)       = ZRNET_V(:) + ZRNET_G(:) + PRNSNOW(:) 
!
IF (LHOOK) CALL DR_HOOK('ISBA_MEB:AVG_FLUXES_MEB_TSPLIT ',1,ZHOOK_HANDLE)
!
END SUBROUTINE AVG_FLUXES_MEB_TSPLIT 
!
!===============================================================================
!
SUBROUTINE SNOWALB_SPECTRAL_BANDS_MEB(PVEGTYPE,PSNOWALB,PSNOWRHO,PSNOWAGE,PPS, &
                                      PSNOWDZ, PZENITH,                        &
                                      PSNOWALBVIS,PSNOWALBNIR,PSNOWALBFIR,     &
                                      PTAU_N)
!
! Split Total snow albedo into N-spectral bands
! NOTE currently MEB only uses 2 bands of the 3 possible.
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_DATA_COVER_PAR, ONLY : NVT_SNOW
USE MODD_MEB_PAR,        ONLY : XSW_WGHT_VIS, XSW_WGHT_NIR
USE MODD_SNOW_PAR,       ONLY : XDSGRAIN_MAX, XSNOW_AGRAIN, XSNOW_BGRAIN,  &
                                XVSPEC1, XVSPEC2, XVSPEC3,                 &
                                XES_CVEXT, XVBETA3, XVBETA5, XMINCOSZEN
!
USE MODE_SNOW3L,         ONLY : SNOW3LALB, SNOW3LRADABS
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PVEGTYPE      ! fraction of each vegetation (-)
REAL, DIMENSION(:),   INTENT(IN)    :: PSNOWALB      ! Snow albedo (total)
REAL, DIMENSION(:),   INTENT(IN)    :: PSNOWRHO      ! Snow layer average density (kg/m3)
REAL, DIMENSION(:),   INTENT(IN)    :: PSNOWDZ       ! Snow layer thickness (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PZENITH       ! Zenith angle (rad)
REAL, DIMENSION(:),   INTENT(IN)    :: PSNOWAGE      ! Snow grain age
REAL, DIMENSION(:),   INTENT(IN)    :: PPS           ! Pressure [Pa]
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOWALBVIS   ! Snow VIS albedo
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOWALBNIR   ! Snow NIR albedo
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOWALBFIR   ! Snow FIR (UV) albedo
REAL, DIMENSION(:),   INTENT(OUT)   :: PTAU_N        ! SW absorption (coef) in uppermost snow layer (-)
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PPS))          :: ZWORK, ZWORKA
REAL, DIMENSION(SIZE(PPS))          :: ZPROJLAT, ZDSGRAIN, ZBETA1, ZBETA2, ZBETA3, &
                                       ZOPTICALPATH1, ZOPTICALPATH2, ZOPTICALPATH3
REAL, DIMENSION(SIZE(PPS))          :: ZPERMSNOWFRAC
REAL, DIMENSION(SIZE(PPS),3)        :: ZSPECTRALALBEDO
!                                      ZSPECTRALALBEDO = spectral albedo (3 bands in algo: 
!                                                        MEB currently uses 2)
!                                                        1=VIS, 2=NIR, 3=UV
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ISBA_MEB:SNOWALB_SPECTRAL_BANDS_MEB',0,ZHOOK_HANDLE)
!
! 1) Spectral albedo
! ------------------
!
ZWORK(:)         = 0.0
ZWORKA(:)        = PSNOWALB(:)
ZPERMSNOWFRAC(:) = PVEGTYPE(:,NVT_SNOW)
!
CALL SNOW3LALB(ZWORKA,ZSPECTRALALBEDO,PSNOWRHO,PSNOWAGE,ZPERMSNOWFRAC,PPS)
!
! Since we only consider VIS and NIR bands for soil and veg in MEB currently:
! (also note, here PSNOWALB doesn't evolve...we just diagnose spectral components).
!
WHERE(PSNOWALB(:)/=XUNDEF)
!
   PSNOWALBVIS(:) = ZSPECTRALALBEDO(:,1)
!
! We diagnose NIR albedo such that total albedo is conserved
! (using just 2 spectral bands in MEB)
!
   PSNOWALBNIR(:) = (PSNOWALB(:) - XSW_WGHT_VIS*PSNOWALBVIS(:))/XSW_WGHT_NIR
!
! currently NOT used by MEB
!
   PSNOWALBFIR(:) = XUNDEF                                     
!
! For the surface layer absorbtion computation:
!
   ZSPECTRALALBEDO(:,1) = PSNOWALBVIS(:)
   ZSPECTRALALBEDO(:,2) = PSNOWALBNIR(:)
   ZSPECTRALALBEDO(:,3) = PSNOWALBFIR(:)
!
ELSEWHERE
!
   PSNOWALBVIS(:) = XUNDEF
   PSNOWALBNIR(:) = XUNDEF
   PSNOWALBFIR(:) = XUNDEF
!
END WHERE
!
!
! 2) SW absorption in uppermost snow layer 
! ----------------------------------------
!
PTAU_N(:) = SNOW3LRADABS(PSNOWRHO,PSNOWDZ,ZSPECTRALALBEDO,PZENITH,ZPERMSNOWFRAC)
!
IF (LHOOK) CALL DR_HOOK('ISBA_MEB:SNOWALB_SPECTRAL_BANDS_MEB',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWALB_SPECTRAL_BANDS_MEB
!===============================================================================

END SUBROUTINE ISBA_MEB
