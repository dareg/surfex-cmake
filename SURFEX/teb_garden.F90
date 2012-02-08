!     #########
    SUBROUTINE TEB_GARDEN (HZ0H, HCOUPLING, TPTIME,                           &
                     PT_CANYON, PQ_CANYON, PU_CANYON,                         &
                     PT_LOWCAN, PQ_LOWCAN, PU_LOWCAN, PZ_LOWCAN,              &
                     PTI_BLD,                                                 &
                     PT_ROOF, PT_ROAD, PT_WALL, PWS_ROOF,PWS_ROAD,            &
                     HSNOW_ROOF,                                              &
                     PWSNOW_ROOF, PTSNOW_ROOF, PRSNOW_ROOF, PASNOW_ROOF,      &
                     PTSSNOW_ROOF, PESNOW_ROOF,                               &
                     HSNOW_ROAD,                                              &
                     PWSNOW_ROAD, PTSNOW_ROAD, PRSNOW_ROAD, PASNOW_ROAD,      &
                     PTSSNOW_ROAD, PESNOW_ROAD,                               &
                     PPEW_A_COEF, PPEW_B_COEF,                                &
                     PPEW_A_COEF_LOWCAN, PPEW_B_COEF_LOWCAN,                  &
                     PPS, PPA, PEXNS, PEXNA,                                  &
                     PTA, PQA, PRHOA, PCO2,                                   &
                     PLW_RAD, PDIR_SW, PSCA_SW, PSW_BANDS, KSW,               &
                     PZENITH,                                                 &
                     PRR, PSR,                                                &
                     PZREF, PUREF, PVMOD,                                     &
                     PH_TRAFFIC, PLE_TRAFFIC, PH_INDUSTRY, PLE_INDUSTRY,      &
                     PTSTEP,                                                  &
                     PZ0_TOWN,                                                &
                     PBLD,PGARDEN,PROAD,                                      &
                     PBLD_HEIGHT,PWALL_O_HOR,PCAN_HW_RATIO,                   &
                     PALB_ROOF, PEMIS_ROOF,                                   &
                     PHC_ROOF,PTC_ROOF,PD_ROOF,                               &
                     PALB_ROAD, PEMIS_ROAD, PSVF_ROAD,                        &
                     PHC_ROAD,PTC_ROAD,PD_ROAD,                               &
                     PALB_WALL, PEMIS_WALL, PSVF_WALL,                        &
                     PSVF_GARDEN,                                             &
                     PHC_WALL,PTC_WALL,PD_WALL,                               &
                     PRN_ROOF, PH_ROOF, PLE_ROOF, PLEW_ROOF, PGFLUX_ROOF,     &
                     PRUNOFF_ROOF,                                            &
                     PRN_ROAD, PH_ROAD, PLE_ROAD, PLEW_ROAD, PGFLUX_ROAD,     &
                     PRUNOFF_ROAD,                                            &
                     PRN_WALL, PH_WALL, PLE_WALL, PGFLUX_WALL,                &
                     PRN_GARDEN,PH_GARDEN,PLE_GARDEN, PGFLUX_GARDEN,          &
                     PRN_BLT,PH_BLT,PLE_BLT, PGFLUX_BLT,                      &
                     PRNSNOW_ROOF, PHSNOW_ROOF, PLESNOW_ROOF, PGSNOW_ROOF,    &
                     PMELT_ROOF,                                              &
                     PRNSNOW_ROAD, PHSNOW_ROAD, PLESNOW_ROAD, PGSNOW_ROAD,    &
                     PMELT_ROAD,                                              &
                     PRN_GRND, PH_GRND, PLE_GRND, PGFLUX_GRND,                &
                     PRN_TOWN, PH_TOWN, PLE_TOWN, PGFLUX_TOWN, PEVAP_TOWN,    &
                     PRUNOFF_TOWN, PSFCO2,                                    &
                     PUW_GRND, PUW_ROOF, PDUWDU_GRND, PDUWDU_ROOF,            &
                     PUSTAR_TOWN, PCD, PCDN, PCH_TOWN, PRI_TOWN,              &
                     PTS_TOWN, PEMIS_TOWN, PDIR_ALB_TOWN, PSCA_ALB_TOWN,      &
                     PRESA_TOWN, PDQS_TOWN, PQF_TOWN, PQF_BLD, PQF_BLDWFR,    &
                     PTI_BLD_EQ, PTI_BLDWFR, PFLX_BLD, PAC_ROAD, PAC_GARDEN,  &
                     PAC_ROAD_WAT, PAC_GARDEN_WAT,                            &
                     PABS_SW_ROOF,PABS_LW_ROOF,                               &
                     PABS_SW_SNOW_ROOF,PABS_LW_SNOW_ROOF,                     &
                     PABS_SW_ROAD,PABS_LW_ROAD,                               &
                     PABS_SW_SNOW_ROAD,PABS_LW_SNOW_ROAD,                     &
                     PABS_SW_WALL,PABS_LW_WALL,                               &
                     PABS_SW_GARDEN,PABS_LW_GARDEN                            )
!   ##########################################################################
!
!!****  *TEB_GARDEN*  
!!
!!    PURPOSE
!!    -------
!
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
USE MODD_CSTS,              ONLY: XTT, XSTEFAN, XCPD, XLVTT
USE MODD_SURF_PAR,          ONLY: XUNDEF
USE MODD_SNOW_PAR,          ONLY: XEMISSN, XANSMAX
USE MODD_ISBA_PAR,          ONLY: XWGMIN
USE MODD_TEB_n,             ONLY: LGARDEN
!
USE MODE_THERMOS
USE MODE_SURF_SNOW_FRAC
!
USE MODI_URBAN_SOLAR_ABS
USE MODI_URBAN_LW_COEF
USE MODI_GARDEN_PROPERTIES
USE MODI_GARDEN
USE MODI_TEB
USE MODI_AVG_URBAN_FLUXES
USE MODI_FLAG_TEB_GARDEN_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    Declarations of arguments
!
CHARACTER(LEN=6)    , INTENT(IN)    :: HZ0H               ! TEB option for z0h roof & road
!                                                         ! 'MASC95' : Mascart et al 1995
!                                                         ! 'BRUT82' : Brustaert     1982
!                                                         ! 'KAND07' : Kanda         2007
CHARACTER(LEN=*)    , INTENT(IN)    :: HCOUPLING          ! type of coupling
                                                          ! 'E' : explicit
                                                          ! 'I' : implicit
TYPE(DATE_TIME)     , INTENT(IN)    :: TPTIME             ! current date and time from teb
!                                                         
REAL, DIMENSION(:)  , INTENT(INOUT) :: PT_CANYON          ! canyon air temperature
REAL, DIMENSION(:)  , INTENT(INOUT) :: PQ_CANYON          ! canyon air specific humidity
REAL, DIMENSION(:)  , INTENT(IN)    :: PU_CANYON          ! canyon hor. wind
REAL, DIMENSION(:)  , INTENT(IN)    :: PU_LOWCAN          ! wind near the road
REAL, DIMENSION(:)  , INTENT(IN)    :: PT_LOWCAN          ! temp. near the road
REAL, DIMENSION(:)  , INTENT(IN)    :: PQ_LOWCAN          ! hum. near the road
REAL, DIMENSION(:)  , INTENT(IN)    :: PZ_LOWCAN          ! height of atm. var. near the road
REAL, DIMENSION(:)  , INTENT(INOUT) :: PTI_BLD            ! inside building temperature
REAL, DIMENSION(:)  , INTENT(INOUT) :: PTI_BLD_EQ         ! inside building temperature
                                                          ! computed by its own evolution equation
REAL, DIMENSION(:)  , INTENT(INOUT) :: PTI_BLDWFR         ! inside building temperature
                                                          ! computed by its own evolution equation
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_ROOF            ! roof layers temperatures
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_ROAD            ! road layers temperatures
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_WALL            ! wall layers temperatures
REAL, DIMENSION(:)  , INTENT(INOUT) :: PWS_ROOF           ! roof water reservoir
REAL, DIMENSION(:)  , INTENT(INOUT) :: PWS_ROAD           ! road water reservoir
CHARACTER(LEN=*)    , INTENT(IN)    :: HSNOW_ROOF         ! snow roof scheme 'NONE', 'D95 ', '1-L '
CHARACTER(LEN=*)    , INTENT(IN)    :: HSNOW_ROAD         ! snow road scheme 'NONE', 'D95 ', '1-L '
REAL, DIMENSION(:,:), INTENT(INOUT) :: PWSNOW_ROOF        ! snow layers reservoir
REAL, DIMENSION(:,:), INTENT(INOUT) :: PTSNOW_ROOF        ! snow layers temperature
REAL, DIMENSION(:,:), INTENT(INOUT) :: PRSNOW_ROOF        ! snow layers density
REAL, DIMENSION(:)  , INTENT(INOUT) :: PASNOW_ROOF        ! snow albedo
REAL, DIMENSION(:)  , INTENT(INOUT) :: PESNOW_ROOF        ! snow emissivity
REAL, DIMENSION(:)  , INTENT(INOUT) :: PTSSNOW_ROOF       ! snow surface temperature
REAL, DIMENSION(:,:), INTENT(INOUT) :: PWSNOW_ROAD        ! snow layers reservoir
REAL, DIMENSION(:,:), INTENT(INOUT) :: PTSNOW_ROAD        ! snow layers temperature
REAL, DIMENSION(:,:), INTENT(INOUT) :: PRSNOW_ROAD        ! snow layers density
REAL, DIMENSION(:)  , INTENT(INOUT) :: PASNOW_ROAD        ! snow albedo
REAL, DIMENSION(:)  , INTENT(INOUT) :: PESNOW_ROAD        ! snow emissivity
REAL, DIMENSION(:)  , INTENT(INOUT) :: PTSSNOW_ROAD       ! snow surface temperature
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEW_A_COEF        ! implicit coefficients
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEW_B_COEF        ! for wind coupling
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEW_A_COEF_LOWCAN ! implicit coefficients for wind coupling
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEW_B_COEF_LOWCAN ! between low canyon wind and road
REAL, DIMENSION(:)  , INTENT(IN)    :: PPS                ! pressure at the surface
REAL, DIMENSION(:)  , INTENT(IN)    :: PPA                ! pressure at the first atmospheric level
REAL, DIMENSION(:)  , INTENT(IN)    :: PEXNS              ! surface exner function
REAL, DIMENSION(:)  , INTENT(IN)    :: PTA                ! temperature at the lowest level
REAL, DIMENSION(:)  , INTENT(IN)    :: PQA                ! specific humidity at the lowest level
REAL, DIMENSION(:)  , INTENT(IN)    :: PVMOD              ! module of the horizontal wind
REAL, DIMENSION(:)  , INTENT(IN)    :: PEXNA              ! exner function at the lowest level
REAL, DIMENSION(:)  , INTENT(IN)    :: PRHOA              ! air density at the lowest level
REAL, DIMENSION(:)  , INTENT(IN)    :: PCO2               ! CO2 concentration in the air    (kg/m3)
REAL, DIMENSION(:)  , INTENT(IN)    :: PLW_RAD            ! atmospheric infrared radiation
REAL, DIMENSION(:,:), INTENT(IN)    :: PDIR_SW            ! incoming direct solar rad on an horizontal surface
REAL, DIMENSION(:,:), INTENT(IN)    :: PSCA_SW            ! scattered incoming solar rad.
REAL, DIMENSION(:)  , INTENT(IN)    :: PSW_BANDS          ! mean wavelength of each shortwave band (m)
INTEGER,              INTENT(IN)    :: KSW                ! number of short-wave spectral bands
REAL, DIMENSION(:)  , INTENT(IN)    :: PZENITH            ! solar zenithal angle
REAL, DIMENSION(:)  , INTENT(IN)    :: PRR                ! rain rate
REAL, DIMENSION(:)  , INTENT(IN)    :: PSR                ! snow rate
REAL, DIMENSION(:)  , INTENT(IN)    :: PH_TRAFFIC         ! anthropogenic sensible heat fluxes due to traffic
REAL, DIMENSION(:)  , INTENT(IN)    :: PLE_TRAFFIC        ! anthropogenic latent heat fluxes due to traffic
REAL, DIMENSION(:)  , INTENT(IN)    :: PH_INDUSTRY        ! anthropogenic sensible heat fluxes due to factories
REAL, DIMENSION(:)  , INTENT(IN)    :: PLE_INDUSTRY       ! anthropogenic latent heat fluxes due to factories
REAL, DIMENSION(:)  , INTENT(IN)    :: PZREF              ! reference height of the first atm level (temperature)
REAL, DIMENSION(:)  , INTENT(IN)    :: PUREF              ! reference height of the first atm level (wind)
REAL                , INTENT(IN)    :: PTSTEP             ! time step
REAL, DIMENSION(:)  , INTENT(IN)    :: PZ0_TOWN           ! town roughness length for momentum
REAL, DIMENSION(:)  , INTENT(IN)    :: PBLD               ! fraction of buildings
REAL, DIMENSION(:)  , INTENT(IN)    :: PGARDEN            ! fraction of green areas
REAL, DIMENSION(:)  , INTENT(IN)    :: PROAD              ! fraction of roads
REAL, DIMENSION(:)  , INTENT(IN)    :: PBLD_HEIGHT        ! buildings h
REAL, DIMENSION(:)  , INTENT(IN)    :: PWALL_O_HOR        ! wall surf. / hor. surf.
REAL, DIMENSION(:)  , INTENT(IN)    :: PCAN_HW_RATIO      ! canyon    h/W
REAL, DIMENSION(:)  , INTENT(IN)    :: PALB_ROOF          ! roof albedo
REAL, DIMENSION(:)  , INTENT(IN)    :: PEMIS_ROOF         ! roof emissivity
REAL, DIMENSION(:,:), INTENT(IN)    :: PHC_ROOF           ! heat capacity for roof layers
REAL, DIMENSION(:,:), INTENT(IN)    :: PTC_ROOF           ! thermal conductivity for roof layers
REAL, DIMENSION(:,:), INTENT(IN)    :: PD_ROOF            ! depth of roof layers
REAL, DIMENSION(:)  , INTENT(IN)    :: PALB_ROAD          ! road albedo
REAL, DIMENSION(:)  , INTENT(IN)    :: PEMIS_ROAD         ! road emissivity
REAL, DIMENSION(:,:), INTENT(IN)    :: PHC_ROAD           ! heat capacity for road layers
REAL, DIMENSION(:,:), INTENT(IN)    :: PTC_ROAD           ! thermal conductivity for road layers
REAL, DIMENSION(:,:), INTENT(IN)    :: PD_ROAD            ! depth of road layers
REAL, DIMENSION(:)  , INTENT(IN)    :: PSVF_ROAD          ! road sky view factor
REAL, DIMENSION(:)  , INTENT(IN)    :: PALB_WALL          ! wall albedo
REAL, DIMENSION(:)  , INTENT(IN)    :: PEMIS_WALL         ! wall emissivity
REAL, DIMENSION(:,:), INTENT(IN)    :: PHC_WALL           ! heat capacity for wall layers
REAL, DIMENSION(:,:), INTENT(IN)    :: PTC_WALL           ! thermal conductivity for wall layers
REAL, DIMENSION(:,:), INTENT(IN)    :: PD_WALL            ! depth of wall layers
REAL, DIMENSION(:)  , INTENT(IN)    :: PSVF_WALL          ! wall sky view factor
REAL, DIMENSION(:)  , INTENT(IN)    :: PSVF_GARDEN        ! green area sky view factor
     !
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN_ROOF           ! net radiation over roof
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_ROOF            ! sensible heat flux over roof
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_ROOF           ! latent heat flux over roof
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLEW_ROOF          ! latent heat flux over roof (snow)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLUX_ROOF        ! flux through the roof
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRUNOFF_ROOF       ! runoff over the ground
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN_ROAD           ! net radiation over road
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_ROAD            ! sensible heat flux over road
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_ROAD           ! latent heat flux over road
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLEW_ROAD          ! latent heat flux over road (snow)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLUX_ROAD        ! flux through the road
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRUNOFF_ROAD       ! runoff over the ground
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN_WALL           ! net radiation over wall
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_WALL            ! sensible heat flux over wall
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_WALL           ! latent heat flux over wall
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLUX_WALL        ! flux through the wall
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN_GARDEN         ! net radiation over green areas
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_GARDEN          ! sensible heat flux over green areas
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_GARDEN         ! latent heat flux over green areas
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLUX_GARDEN      ! flux through the green areas
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN_BLT            ! net radiation over built surf 
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_BLT             ! sensible heat flux over built surf 
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_BLT            ! latent heat flux over built surf 
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLUX_BLT         ! flux through the built surf 
!
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRNSNOW_ROOF       ! net radiation over snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PHSNOW_ROOF        ! sensible heat flux over snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLESNOW_ROOF       ! latent heat flux over snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGSNOW_ROOF        ! flux under the snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PMELT_ROOF         ! snow melt
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRNSNOW_ROAD       ! net radiation over snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PHSNOW_ROAD        ! sensible heat flux over snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLESNOW_ROAD       ! latent heat flux over snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGSNOW_ROAD        ! flux under the snow
REAL, DIMENSION(:)  , INTENT(OUT)   :: PMELT_ROAD         ! snow melt
!
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN_GRND           ! net radiation over ground
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_GRND            ! sensible heat flux over ground
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_GRND           ! latent heat flux over ground
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLUX_GRND        ! flux through the ground
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN_TOWN           ! net radiation over town
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_TOWN            ! sensible heat flux over town
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_TOWN           ! latent heat flux over town
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLUX_TOWN        ! flux through the ground
REAL, DIMENSION(:)  , INTENT(OUT)   :: PEVAP_TOWN         ! evaporation flux (kg/m2/s)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRUNOFF_TOWN       ! runoff over the ground
REAL, DIMENSION(:)  , INTENT(OUT)   :: PSFCO2             ! flux of CO2       (kg/m2/s)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PUW_GRND           ! momentum flux for ground built surf
REAL, DIMENSION(:)  , INTENT(OUT)   :: PUW_ROOF           ! momentum flux for roofs
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDUWDU_GRND        !
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDUWDU_ROOF        !
REAL, DIMENSION(:)  , INTENT(OUT)   :: PUSTAR_TOWN        ! friciton velocity over town
REAL, DIMENSION(:)  , INTENT(OUT)   :: PCD                ! town averaged drag coefficient
REAL, DIMENSION(:)  , INTENT(OUT)   :: PCDN               ! town averaged neutral drag coefficient
REAL, DIMENSION(:)  , INTENT(OUT)   :: PCH_TOWN           ! town averaged heat transfer coefficient
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRI_TOWN           ! town averaged Richardson number
REAL, DIMENSION(:)  , INTENT(OUT)   :: PTS_TOWN           ! town surface temperature
REAL, DIMENSION(:)  , INTENT(OUT)   :: PEMIS_TOWN         ! town equivalent emissivity
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDIR_ALB_TOWN      ! town equivalent direct albedo
REAL, DIMENSION(:)  , INTENT(OUT)   :: PSCA_ALB_TOWN      ! town equivalent diffuse albedo
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRESA_TOWN         ! town aerodynamical resistance
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDQS_TOWN          ! heat storage inside town
REAL, DIMENSION(:)  , INTENT(OUT)   :: PQF_TOWN           ! total anthropogenic heat
REAL, DIMENSION(:)  , INTENT(OUT)   :: PQF_BLD            ! anthropogenic heat flux of domestic heating
REAL, DIMENSION(:)  , INTENT(OUT)   :: PQF_BLDWFR 
REAL, DIMENSION(:)  , INTENT(OUT)   :: PFLX_BLD           ! heat flx from inside bld through its structure
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_ROAD           ! road conductance
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_GARDEN         ! green area conductance
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_ROAD_WAT       ! road conductance for latent heat
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_GARDEN_WAT     ! green area conductance for latent heat
!
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_SW_ROOF       ! absorbed solar rad by roof
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_SW_SNOW_ROOF  ! absorbed solar rad by snow on roof
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_LW_ROOF       ! absorbed IR rad by roof
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_LW_SNOW_ROOF  ! absorbed IR rad by snow on roof
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_SW_ROAD       ! absorbed solar rad by road
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_SW_SNOW_ROAD  ! absorbed solar rad by snow on road
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_LW_ROAD       ! absorbed IR rad by road
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_LW_SNOW_ROAD  ! absorbed IR rad by snow on road
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_SW_WALL       ! absorbed solar rad by wall
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_LW_WALL       ! absorbed IR rad by wall
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_SW_GARDEN     ! absorbed solar rad by green areas
REAL, DIMENSION(:)  , INTENT(OUT)   :: PABS_LW_GARDEN     ! absorbed IR rad by green areas
!
!*      0.2    Declarations of local variables
!
REAL, DIMENSION(SIZE(PTA)) :: ZTA            ! air temperature extrapolated at roof level
REAL, DIMENSION(SIZE(PTA)) :: ZQA            ! air humidity extrapolated at roof level
!
REAL, DIMENSION(SIZE(PTA)) :: ZDN_ROOF       ! snow fraction on roofs
REAL, DIMENSION(SIZE(PTA)) :: ZDN_ROAD       ! snow fraction on roads
REAL, DIMENSION(SIZE(PTA)) :: ZDF_ROOF       ! free-snow fraction on roofs
REAL, DIMENSION(SIZE(PTA)) :: ZDF_ROAD       ! free-snow fraction on roads
REAL, DIMENSION(SIZE(PTA)) :: ZDELT_ROAD     ! fraction of water on roads
REAL, DIMENSION(SIZE(PTA)) :: ZDELT_ROOF     ! fraction of water on roofs
REAL, DIMENSION(SIZE(PTA)) :: ZAC_ROOF       ! roof conductance
REAL, DIMENSION(SIZE(PTA)) :: ZAC_ROOF_WAT   ! roof water conductance
REAL, DIMENSION(SIZE(PTA)) :: ZAC_WALL       ! wall conductance
REAL, DIMENSION(SIZE(PTA)) :: ZAC_TOP        ! top conductance
REAL, DIMENSION(SIZE(PTA)) :: ZQSAT_ROAD     ! hum of saturation for roads
REAL, DIMENSION(SIZE(PTA)) :: ZQSAT_GARDEN   ! hum of saturation for green areas
REAL, DIMENSION(SIZE(PTA)) :: ZQSAT_ROOF     ! hum of saturation for roofs
!
! coefficients for LW computations over snow (from previous time-step)
!
REAL, DIMENSION(SIZE(PTA)) :: ZTSSNOW_ROOF   ! roof snow temp at previous time-step
REAL, DIMENSION(SIZE(PTA)) :: ZTSSNOW_ROAD   ! road snow temp at previous time-step
REAL, DIMENSION(SIZE(PTA)) :: ZESNOW_ROOF    ! snow emissivity at previous time-step
REAL, DIMENSION(SIZE(PTA)) :: ZESNOW_ROAD    ! snow emissivity at previous time-step
!
! incoming shortwave radiation
!
REAL, DIMENSION(SIZE(PTA)) :: ZDIR_SW             ! direct  solar rad
REAL, DIMENSION(SIZE(PTA)) :: ZSCA_SW             ! diffuse solar rad
INTEGER                    :: JSWB
!
! albedo & emissivity
!
REAL, DIMENSION(SIZE(PTA)) :: ZALB_GARDEN    ! albedo     for green areas
REAL, DIMENSION(SIZE(PTA)) :: ZEMIS_GARDEN   ! emissivity for green areas
!
! radiation received by surfaces
!
REAL, DIMENSION(SIZE(PTA)) :: ZREC_SW_ROAD        ! solar rad received by roads
REAL, DIMENSION(SIZE(PTA)) :: ZREC_SW_WALL        ! solar rad received by walls
REAL, DIMENSION(SIZE(PTA)) :: ZREC_SW_GARDEN      ! solar rad received by gardens
REAL, DIMENSION(SIZE(PTA)) :: ZREC_SW_SNOW_ROAD   ! solar rad received by snow on roads
!
REAL, DIMENSION(SIZE(PTA)) :: ZREC_LW_GARDEN      ! IR rad received by gardens
!
REAL, DIMENSION(SIZE(PTA)) :: ZSW_RAD_GARDEN      ! solar radiation reaching urban green areas
!
! coefficients for LW contributions
!
REAL, DIMENSION(SIZE(PTA)) :: ZLW_W_TO_W          ! LW contrib. wall       -> wall
REAL, DIMENSION(SIZE(PTA)) :: ZLW_W_TO_R          ! LW contrib. wall       -> road
REAL, DIMENSION(SIZE(PTA)) :: ZLW_W_TO_G          ! LW contrib. wall       -> green
REAL, DIMENSION(SIZE(PTA)) :: ZLW_W_TO_NR         ! LW contrib. wall       -> road(snow)
REAL, DIMENSION(SIZE(PTA)) :: ZLW_R_TO_W          ! LW contrib. road       -> wall
REAL, DIMENSION(SIZE(PTA)) :: ZLW_R_TO_R          ! LW contrib. road       -> road
REAL, DIMENSION(SIZE(PTA)) :: ZLW_R_TO_G          ! LW contrib. road       -> green
REAL, DIMENSION(SIZE(PTA)) :: ZLW_R_TO_NR         ! LW contrib. road       -> road(snow)
REAL, DIMENSION(SIZE(PTA)) :: ZLW_G_TO_W          ! LW contrib. green      -> wall
REAL, DIMENSION(SIZE(PTA)) :: ZLW_G_TO_R          ! LW contrib. green      -> road
REAL, DIMENSION(SIZE(PTA)) :: ZLW_G_TO_G          ! LW contrib. green      -> green
REAL, DIMENSION(SIZE(PTA)) :: ZLW_G_TO_NR         ! LW contrib. green      -> road(snow)
REAL, DIMENSION(SIZE(PTA)) :: ZLW_NR_TO_W         ! LW contrib. road(snow) -> wall
REAL, DIMENSION(SIZE(PTA)) :: ZLW_NR_TO_R         ! LW contrib. road(snow) -> road
REAL, DIMENSION(SIZE(PTA)) :: ZLW_NR_TO_G         ! LW contrib. road(snow) -> green
REAL, DIMENSION(SIZE(PTA)) :: ZLW_NR_TO_NR        ! LW contrib. road(snow) -> road(snow)
REAL, DIMENSION(SIZE(PTA)) :: ZLW_S_TO_W          ! LW contrib. sky        -> wall
REAL, DIMENSION(SIZE(PTA)) :: ZLW_S_TO_R          ! LW contrib. sky        -> road
REAL, DIMENSION(SIZE(PTA)) :: ZLW_S_TO_G          ! LW contrib. sky        -> green
REAL, DIMENSION(SIZE(PTA)) :: ZLW_S_TO_NR         ! LW contrib. sky        -> road(snow)
!
! local variable at previous time-step
!
REAL, DIMENSION(SIZE(PTA))  :: ZPET_A_COEF          
REAL, DIMENSION(SIZE(PTA))  :: ZPET_B_COEF          
REAL, DIMENSION(SIZE(PTA))  :: ZPEQ_A_COEF          
REAL, DIMENSION(SIZE(PTA))  :: ZPEQ_B_COEF          
!
REAL, DIMENSION(SIZE(PTA))  :: ZUW_ROAD           ! momentum flux for roads
REAL, DIMENSION(SIZE(PTA))  :: ZUW_GARDEN         ! momentum flux for green areas
REAL, DIMENSION(SIZE(PTA))  :: ZDUWDU_ROAD        !
!
REAL, DIMENSION(SIZE(PTA))  :: ZAC_AGG_GARDEN     ! aggreg. aeodynamic resistance for green areas
REAL, DIMENSION(SIZE(PTA))  :: ZHU_AGG_GARDEN     ! aggreg. relative humidity for green areas
!
!  surfaces relative fractions
!
REAL, DIMENSION(SIZE(PTA)) :: ZROOF_FRAC          ! roof, wall and
REAL, DIMENSION(SIZE(PTA)) :: ZWALL_FRAC          ! road fractions
REAL, DIMENSION(SIZE(PTA)) :: ZROAD_FRAC          ! of exchange surf.
REAL, DIMENSION(SIZE(PTA)) :: ZGARDEN_FRAC        !                  
REAL, DIMENSION(SIZE(PTA)) :: ZTOTS_O_HORS        ! total canyon+roof surface
!                                                 ! over horizontal surface
REAL, DIMENSION(SIZE(PTA)) :: ZWALL_O_ROAD        ! wall surface over road surface
REAL, DIMENSION(SIZE(PTA)) :: ZWALL_O_GRND        ! wall surface over (road+green area) surface
!
! garden temperatures
!
REAL, DIMENSION(SIZE(PTA)) :: ZTS_GARDEN          ! surface temperature of urban green areas at t
!
! fluxes from green surfaces
!
REAL, DIMENSION(SIZE(PTA)) :: ZEVAP_GARDEN        ! evaporation (kg/m2/s)
REAL, DIMENSION(SIZE(PTA)) :: ZSFCO2_GARDEN       ! CO2 fluxes (kg/m2/s)
!
INTEGER :: JJ
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!*      1.     Initializations
!              ---------------
!
!*      1.0    broadband radiative fluxes
!              --------------------------
!
IF (LHOOK) CALL DR_HOOK('TEB_GARDEN',0,ZHOOK_HANDLE)
ZDIR_SW(:) = 0.
ZSCA_SW(:) = 0.
!
DO JSWB=1,KSW
  DO JJ=1,SIZE(PDIR_SW,1)
    ZDIR_SW(JJ) = ZDIR_SW(JJ) + PDIR_SW(JJ,JSWB)
    ZSCA_SW(JJ) = ZSCA_SW(JJ) + PSCA_SW(JJ,JSWB)
  ENDDO
END DO
!
!
!*      1.1    surfaces relative fractions
!              ---------------------------
!
DO JJ=1,SIZE(PROAD)
  IF (PROAD(JJ) .GT. 0.) THEN
    ZTOTS_O_HORS(JJ) = 1. + PWALL_O_HOR(JJ)
    ZROOF_FRAC(JJ)   = PBLD(JJ)        / ZTOTS_O_HORS(JJ)
    ZWALL_FRAC(JJ)   = PWALL_O_HOR(JJ) / ZTOTS_O_HORS(JJ)
    ZROAD_FRAC(JJ)   = PROAD(JJ)       / ZTOTS_O_HORS(JJ)
    ZGARDEN_FRAC(JJ) = PGARDEN(JJ)     / ZTOTS_O_HORS(JJ)
    ZWALL_O_ROAD(JJ) = ZWALL_FRAC(JJ) /  ZROAD_FRAC(JJ)
    ZWALL_O_GRND(JJ) = ZWALL_FRAC(JJ) / (ZROAD_FRAC(JJ)+ZGARDEN_FRAC(JJ))
  ELSEIF (PBLD(JJ) .EQ. 0. .AND. PROAD(JJ) .EQ. 0.) THEN
    ZTOTS_O_HORS(JJ) = 1. + PWALL_O_HOR(JJ)
    ZROOF_FRAC(JJ)   = 0.
    ZWALL_FRAC(JJ)   = PWALL_O_HOR(JJ) / ZTOTS_O_HORS(JJ)
    ZROAD_FRAC(JJ)   = 0.
    ZGARDEN_FRAC(JJ) = PGARDEN(JJ)     / ZTOTS_O_HORS(JJ)
    ZWALL_O_ROAD(JJ) = 0.
    ZWALL_O_GRND(JJ) = ZWALL_FRAC(JJ) / ZGARDEN_FRAC(JJ)
  ENDIF
ENDDO
!
!-------------------------------------------------------------------------------
!
!*      2.     Snow-covered surfaces relative effects
!              --------------------------------------
!
!*      2.1    Snow-covered surfaces relative fractions (at previous time-step)
!              ----------------------------------------
CALL SNOW_FRAC_ROAD(PWSNOW_ROAD(:,1),PSR(:)>0.,ZDN_ROAD,ZDF_ROAD)
CALL SNOW_FRAC_ROOF(PWSNOW_ROOF(:,1),PSR(:)>0.,ZDN_ROOF,ZDF_ROOF)
!
!* new snow albedo
!
WHERE (PWSNOW_ROAD(:,1)==0. .AND. PSR(:)>0.) PASNOW_ROAD(:) = XANSMAX
WHERE (PWSNOW_ROOF(:,1)==0. .AND. PSR(:)>0.) PASNOW_ROOF(:) = XANSMAX
!
!*      2.2    If snow was not present at previous time-step but is falling
!              ------------------------------------------------------------
!
WHERE (PWSNOW_ROAD(:,1)==0. .AND. PSR(:)>0.)
  PASNOW_ROAD(:) = XANSMAX
  PESNOW_ROAD(:) = XEMISSN
  PTSSNOW_ROAD(:)= MIN(PT_ROAD(:,1), XTT)
END WHERE
WHERE (PWSNOW_ROOF(:,1)==0. .AND. PSR(:)>0.)
  PASNOW_ROOF(:) = XANSMAX
  PESNOW_ROOF(:) = XEMISSN
  PTSSNOW_ROOF(:)= MIN(PT_ROOF(:,1), XTT)
END WHERE
!
!*      2.3    Radiative snow variables at previous time-step
!              ----------------------------------------------
!
ZESNOW_ROOF  (:) = PESNOW_ROOF  (:)
ZESNOW_ROAD  (:) = PESNOW_ROAD  (:)
ZTSSNOW_ROOF (:) = PTSSNOW_ROOF (:)
ZTSSNOW_ROAD (:) = PTSSNOW_ROAD (:)
!
!-------------------------------------------------------------------------------
!
!*      4.     Extrapolation of atmospheric T and q at roof level (for fluxes computation)
!              --------------------------------------------------
!
ZTA(:) = PTA(:) * PEXNS(:) / PEXNA(:)
ZQA(:) = PQA(:) * QSAT(PTA(:),PPS(:)) / QSAT(ZTA(:),PPA(:))
!
!-------------------------------------------------------------------------------
!
!*      1.3    Set physical values for points where there is no garden
!              -------------------------------------------------------
!
! This way, ISBA can run without problem for these points
!
CALL FLAG_TEB_GARDEN_n(1)

!-------------------------------------------------------------------------------
!
!*      5.     Grid-averaged albedo and emissivity of green areas
!              --------------------------------------------------
!
!
ZALB_GARDEN   = XUNDEF
ZEMIS_GARDEN  = XUNDEF
ZTS_GARDEN    = XUNDEF
!
IF (LGARDEN) THEN
 CALL GARDEN_PROPERTIES(PDIR_SW, PSCA_SW, PSW_BANDS, KSW,     &
                        ZTS_GARDEN, ZEMIS_GARDEN, ZALB_GARDEN )
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      6.     Solar radiation
!              ---------------
!
!* for sake of simplicity, the radiative properties of the vegetation (if any)
!  situated in the canyon (and then influencing h/w) is replaced by the road
!  properties. This influences only the radiative part due to reflections,
!  which are less important when h/w is low (i.e. large streets with garden),
!  than when h/w is large (from 0.5 and above approximately). This case
!  (h/w large) occurs in city downtowns, and then does not occur in presence
!  of significant vegetation area.
!
CALL URBAN_SOLAR_ABS(ZDIR_SW, ZSCA_SW, PZENITH,                    &
                     PBLD, PGARDEN, PROAD,                         &
                     PWALL_O_HOR, PCAN_HW_RATIO,                   &
                     PALB_ROOF,                                    &
                     PALB_ROAD, PSVF_ROAD, PALB_WALL, PSVF_WALL,   &
                     ZALB_GARDEN, PSVF_GARDEN,                     &
                     PASNOW_ROOF, PASNOW_ROAD,                     &
                     ZDN_ROOF, ZDF_ROOF, ZDN_ROAD, ZDF_ROAD,       &
                     PABS_SW_ROOF, PABS_SW_ROAD, PABS_SW_WALL,     &
                     PABS_SW_GARDEN,                               &
                     PABS_SW_SNOW_ROOF, PABS_SW_SNOW_ROAD,         &
                     ZREC_SW_ROAD,  ZREC_SW_SNOW_ROAD,             &
                     ZREC_SW_WALL,                                 &
                     ZREC_SW_GARDEN,                               &
                     PDIR_ALB_TOWN,PSCA_ALB_TOWN,                  &
                     ZSW_RAD_GARDEN                                )
!
!-------------------------------------------------------------------------------
!
!*      7.     LW properties
!              -------------
!
!*      7.1    Snow-free roads
!              ---------------
!
CALL URBAN_LW_COEF(PEMIS_ROAD,  PSVF_ROAD,                         &
                   PEMIS_WALL,  PSVF_WALL,                         &
                   ZEMIS_GARDEN,  PSVF_GARDEN,                     &
                   PROAD, PGARDEN,                                 &
                   ZDN_ROAD,   ZDF_ROAD,   PESNOW_ROAD,            &
                   ZLW_W_TO_W , ZLW_W_TO_R , ZLW_W_TO_G ,          &
                   ZLW_R_TO_W , ZLW_R_TO_R , ZLW_R_TO_G ,          &
                   ZLW_G_TO_W , ZLW_G_TO_R , ZLW_G_TO_G ,          &
                   ZLW_S_TO_W , ZLW_S_TO_R , ZLW_S_TO_G ,          &
                   ZLW_NR_TO_W, ZLW_NR_TO_R, ZLW_NR_TO_G           )
!
!*      7.2    Snow-covered roads
!              ------------------
!
  CALL URBAN_LW_COEF(PESNOW_ROAD,  PSVF_ROAD,                            &
                     PEMIS_WALL,   PSVF_WALL,                            &
                     ZEMIS_GARDEN, PSVF_GARDEN,                          &
                     PROAD, PGARDEN,                                     &
                     ZDF_ROAD,   ZDN_ROAD,   PEMIS_ROAD,                 &
                     ZLW_W_TO_W , ZLW_W_TO_NR , ZLW_W_TO_G ,             &
                     ZLW_NR_TO_W, ZLW_NR_TO_NR, ZLW_NR_TO_G,             &
                     ZLW_G_TO_W,  ZLW_G_TO_NR,  ZLW_G_TO_G,              &
                     ZLW_S_TO_W , ZLW_S_TO_NR , ZLW_S_TO_G ,             &
                     ZLW_R_TO_W , ZLW_R_TO_NR , ZLW_R_TO_G               )
!
!-------------------------------------------------------------------------------
!
!*      8.     Terms of radiation absorption
!              -----------------------------
!
!*      8.2    IR rad received by gardens (snow free and snow covered separately)
!              --------------------------
!
IF (LGARDEN) THEN
  ZREC_LW_GARDEN   (:) = ZLW_S_TO_G (:)*PLW_RAD       (:)            &
                       + ZLW_W_TO_G (:)*PT_WALL       (:,1)**4       &
                       + ZLW_R_TO_G (:)*PT_ROAD       (:,1)**4       &
                       + ZLW_G_TO_G (:)*ZTS_GARDEN    (:)**4         &
                       + ZLW_NR_TO_G(:)*ZTSSNOW_ROAD  (:)**4         &
                       + XSTEFAN*ZTS_GARDEN(:)**4
!
ELSE
  ZREC_LW_GARDEN      (:) = XUNDEF
END IF
!
!-------------------------------------------------------------------------------
!
!*      9.     Treatment of green areas
!              ------------------------
!
!*      9.1    Implicit coeefs for T and Q
!              ---------------------------
!
!* explicit coupling for the time being.
!  canopy may need implicitation if there is a lot a garden in the grid mesh
!
ZPET_A_COEF(:) = 0.
ZPET_B_COEF(:) = PT_LOWCAN(:) / PEXNS(:)
ZPEQ_A_COEF(:) = 0.
ZPEQ_B_COEF(:) = PQ_LOWCAN(:)
!
!
IF (LGARDEN) THEN
!
!*      9.2    Call ISBA for green areas
!              -------------------------
!
  CALL GARDEN(TPTIME, PPEW_A_COEF_LOWCAN, PPEW_B_COEF_LOWCAN,                      &
              ZPET_A_COEF, ZPEQ_A_COEF, ZPET_B_COEF, ZPEQ_B_COEF,                  &
              PTSTEP, PZ_LOWCAN,                                                   &
              PT_LOWCAN, PQ_LOWCAN, PEXNS, PRHOA, PCO2, PPS, PRR, PSR, PZENITH,    &
              ZREC_SW_GARDEN, ZREC_LW_GARDEN, PU_LOWCAN,                           &
              PRN_GARDEN,PH_GARDEN,PLE_GARDEN,PGFLUX_GARDEN, ZSFCO2_GARDEN,        &
              ZEVAP_GARDEN, ZUW_GARDEN,                                            &
              PAC_GARDEN,ZQSAT_GARDEN,ZTS_GARDEN,                                  &
              ZAC_AGG_GARDEN, ZHU_AGG_GARDEN                                       )

  PAC_GARDEN_WAT(:) = PAC_GARDEN(:)
  PABS_SW_GARDEN(:) = (1.-ZALB_GARDEN(:)) * ZREC_SW_GARDEN
  PABS_LW_GARDEN(:) = ZEMIS_GARDEN * ZREC_LW_GARDEN - XSTEFAN * ZEMIS_GARDEN * ZTS_GARDEN**4

ELSE
!
 PRN_GARDEN    (:) = 0.
 PH_GARDEN     (:) = 0.
 PLE_GARDEN    (:) = 0.
 PGFLUX_GARDEN (:) = 0.
 ZUW_GARDEN    (:) = 0.
 PAC_GARDEN    (:) = 0.
 PGFLUX_GARDEN (:) = 0.
 ZEVAP_GARDEN  (:) = 0.
 ZSFCO2_GARDEN (:) = 0.
 ZQSAT_GARDEN  (:) = XUNDEF
 ZTS_GARDEN    (:) = XUNDEF
 ZAC_AGG_GARDEN(:) = XUNDEF
 ZHU_AGG_GARDEN(:) = XUNDEF
 PABS_SW_GARDEN(:) = XUNDEF
 PABS_LW_GARDEN(:) = XUNDEF
!
ENDIF

!-------------------------------------------------------------------------------
!
!*     10.     Treatment of built covers
!              -------------------------
!
  CALL TEB  (HZ0H, HCOUPLING, PT_CANYON, PQ_CANYON, PU_CANYON,        &
             PT_LOWCAN, PQ_LOWCAN, PU_LOWCAN, PZ_LOWCAN,              &
             PTI_BLD,                                                 &
             PT_ROOF, PT_ROAD, PT_WALL, PWS_ROOF, PWS_ROAD,           &
             HSNOW_ROOF,                                              &
             PWSNOW_ROOF, PTSNOW_ROOF, PRSNOW_ROOF, PASNOW_ROOF,      &
             PTSSNOW_ROOF, PESNOW_ROOF,                               &
             HSNOW_ROAD,                                              &
             PWSNOW_ROAD, PTSNOW_ROAD, PRSNOW_ROAD, PASNOW_ROAD,      &
             PTSSNOW_ROAD, PESNOW_ROAD,                               &
             PPEW_A_COEF, PPEW_B_COEF,                                &
             PPEW_A_COEF_LOWCAN, PPEW_B_COEF_LOWCAN,                  &
             PPS, PPA, PEXNS, PEXNA,                                  &
             PTA, PQA, PRHOA,                                         &
             PLW_RAD, ZDIR_SW, ZSCA_SW, PZENITH,                      &
             PRR, PSR,                                                &
             PZREF, PUREF, PVMOD,                                     &
             PH_TRAFFIC, PLE_TRAFFIC, PH_INDUSTRY, PLE_INDUSTRY,      &
             PTSTEP,                                                  &
             PZ0_TOWN,                                                &
             PBLD, PGARDEN, PROAD,                                    &
             PBLD_HEIGHT, PWALL_O_HOR, PCAN_HW_RATIO,                 &
             ZDF_ROOF, ZDN_ROOF, ZDF_ROAD, ZDN_ROAD,                  &
             ZQSAT_ROOF, ZQSAT_ROAD, ZDELT_ROOF, ZDELT_ROAD,          &
             PALB_ROOF, PEMIS_ROOF,                                   &
             PHC_ROOF,PTC_ROOF,PD_ROOF,                               &
             PALB_ROAD, PEMIS_ROAD, PSVF_ROAD,                        &
             PHC_ROAD,PTC_ROAD,PD_ROAD,                               &
             PALB_WALL, PEMIS_WALL, PSVF_WALL,                        &
             ZTS_GARDEN,                                              &
             PHC_WALL,PTC_WALL,PD_WALL,                               &
             PRN_ROOF, PH_ROOF, PLE_ROOF, PLEW_ROOF, PGFLUX_ROOF,     &
             PRUNOFF_ROOF,                                            &
             PRN_ROAD, PH_ROAD, PLE_ROAD, PLEW_ROAD, PGFLUX_ROAD,     &
             PRUNOFF_ROAD,                                            &
             PRN_WALL, PH_WALL, PLE_WALL, PGFLUX_WALL,                &
             PRN_BLT,PH_BLT,PLE_BLT,PGFLUX_BLT,                       &
             PRNSNOW_ROOF, PHSNOW_ROOF, PLESNOW_ROOF, PGSNOW_ROOF,    &
             PMELT_ROOF,                                              &
             PRNSNOW_ROAD, PHSNOW_ROAD, PLESNOW_ROAD, PGSNOW_ROAD,    &
             PMELT_ROAD,                                              &
             PRUNOFF_TOWN,                                            &
             ZUW_ROAD, PUW_ROOF, ZDUWDU_ROAD, PDUWDU_ROOF,            &
             PUSTAR_TOWN, PCD, PCDN, PCH_TOWN, PRI_TOWN,              &
             PRESA_TOWN, PDQS_TOWN, PQF_TOWN, PQF_BLD, PQF_BLDWFR,    &
             PTI_BLD_EQ, PTI_BLDWFR, PFLX_BLD,                        &
             ZAC_ROOF, PAC_ROAD, ZAC_WALL, ZAC_TOP, PAC_GARDEN,       &
             ZAC_ROOF_WAT, PAC_ROAD_WAT,                              &
             PABS_SW_ROOF, PABS_LW_ROOF,                              &
             PABS_SW_SNOW_ROOF, PABS_LW_SNOW_ROOF,                    &
             PABS_SW_ROAD, PABS_LW_ROAD,                              &
             PABS_SW_SNOW_ROAD, PABS_LW_SNOW_ROAD,                    &
             PABS_SW_WALL, PABS_LW_WALL,                              &
             ZLW_W_TO_W ,  ZLW_W_TO_R  , ZLW_W_TO_G ,                 &
             ZLW_W_TO_NR ,                                            &
             ZLW_R_TO_W ,  ZLW_R_TO_R  , ZLW_R_TO_G ,                 &
             ZLW_R_TO_NR ,                                            &
             ZLW_G_TO_W ,  ZLW_G_TO_R  , ZLW_G_TO_G ,                 &
             ZLW_G_TO_NR ,                                            &
             ZLW_S_TO_W ,  ZLW_S_TO_R  , ZLW_S_TO_G ,                 &
             ZLW_S_TO_NR ,                                            &
             ZLW_NR_TO_W,  ZLW_NR_TO_R , ZLW_NR_TO_G,                 &
             ZLW_NR_TO_NR                                             )
!
!-------------------------------------------------------------------------------
!
!*     11.     Aggregation
!              -----------
!
CALL AVG_URBAN_FLUXES(PTS_TOWN, PEMIS_TOWN,                                    &
                     PT_CANYON, PQ_CANYON,                                     &
                     PT_LOWCAN, PQ_LOWCAN,                                     &
                     PT_ROOF(:,1),PT_ROAD(:,1),PT_WALL(:,1), ZTS_GARDEN,       &
                     ZTA, ZQA, PRHOA, PPS,                                     &
                     PH_TRAFFIC, PLE_TRAFFIC, PH_INDUSTRY, PLE_INDUSTRY,       &
                     PBLD, PROAD, PGARDEN, PWALL_O_HOR, ZWALL_O_GRND,          &
                     PEMIS_ROOF,                                               &
                     ZESNOW_ROOF,                                              &
                     PLW_RAD, ZLW_S_TO_W, ZLW_S_TO_R, ZLW_S_TO_G,              &
                     ZLW_S_TO_NR,                                              &
                     PABS_LW_ROOF, PABS_LW_WALL, PABS_LW_ROAD, PABS_LW_GARDEN, &
                     PABS_LW_SNOW_ROOF, PABS_LW_SNOW_ROAD,                     &
                     ZAC_ROOF, ZAC_ROOF_WAT,                                   &
                     ZAC_WALL, PAC_ROAD, PAC_ROAD_WAT, ZAC_TOP,                &
                     PAC_GARDEN,                                               &
                     ZQSAT_GARDEN, ZAC_AGG_GARDEN, ZHU_AGG_GARDEN,             &
                     ZQSAT_ROOF, ZQSAT_ROAD,                                   &
                     ZDELT_ROOF, ZDELT_ROAD,                                   &
                     ZROOF_FRAC, ZWALL_FRAC, ZROAD_FRAC, ZGARDEN_FRAC,         &
                     ZTOTS_O_HORS,                                             &
                     ZDF_ROOF, ZDN_ROOF, ZDF_ROAD, ZDN_ROAD,                   &
                     PRN_ROOF, PH_ROOF, PLE_ROOF, PGFLUX_ROOF,                 &
                     PRN_ROAD, PH_ROAD, PLE_ROAD, PGFLUX_ROAD,                 &
                     PRN_GARDEN, PH_GARDEN, PLE_GARDEN, PGFLUX_GARDEN,         &
                     PRN_WALL, PH_WALL, PLE_WALL, PGFLUX_WALL,                 &
                     PLEW_ROOF, PLESNOW_ROOF,                                  &
                     PLEW_ROAD, PLESNOW_ROAD, PHSNOW_ROAD,                     &
                     ZEVAP_GARDEN,                                             &
                     PRN_GRND, PH_GRND, PLE_GRND, PGFLUX_GRND,                 &
                     PRN_TOWN, PH_TOWN, PLE_TOWN, PGFLUX_TOWN, PEVAP_TOWN      )

PSFCO2(:) = PGARDEN(:) * ZSFCO2_GARDEN(:)  ! no CO2 flux from built and road yet.
!
!-------------------------------------------------------------------------------
!
!*     12.     Momentum flux for ground built surfaces
!              ---------------------------------------
!
WHERE (PROAD(:)+PGARDEN(:).NE.0.) 
        PUW_GRND (:)     = (PROAD(:)*ZUW_ROAD(:) + PGARDEN(:)*ZUW_GARDEN(:)) / (PROAD(:)+PGARDEN(:))
ELSEWHERE
        PUW_GRND (:)     = 0.
ENDWHERE
!
PDUWDU_GRND (:)  = 0.
IF (LHOOK) CALL DR_HOOK('TEB_GARDEN',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE TEB_GARDEN
