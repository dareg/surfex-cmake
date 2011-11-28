!     #########
    SUBROUTINE TEB  (HZ0H, HCOUPLING, PT_CANYON, PQ_CANYON, PU_CANYON,        &
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
                     PTA, PQA, PRHOA,                                         &
                     PLW_RAD, PDIR_SW_RAD, PSCA_SW_RAD, PZENITH,              &
                     PRR, PSR,                                                &
                     PZREF, PUREF, PVMOD,                                     &
                     PH_TRAFFIC, PLE_TRAFFIC, PH_INDUSTRY, PLE_INDUSTRY,      &
                     PTSTEP,                                                  &
                     PZ0_TOWN,                                                &
                     PBLD,PGARDEN,PROAD,                                      &
                     PBLD_HEIGHT,PWALL_O_HOR,PCAN_HW_RATIO,                   &
                     PDF_ROOF, PDN_ROOF, PDF_ROAD, PDN_ROAD,                  &
                     PQSAT_ROOF, PQSAT_ROAD, PDELT_ROOF, PDELT_ROAD,          &
                     PALB_ROOF, PEMIS_ROOF,                                   &
                     PHC_ROOF,PTC_ROOF,PD_ROOF,                               &
                     PALB_ROAD, PEMIS_ROAD, PSVF_ROAD,                        &
                     PHC_ROAD,PTC_ROAD,PD_ROAD,                               &
                     PALB_WALL, PEMIS_WALL, PSVF_WALL,                        &
                     PTS_GARDEN,                                              &
                     PHC_WALL,PTC_WALL,PD_WALL,                               &
                     PRN_ROOF, PH_ROOF, PLE_ROOF, PLEW_ROOF, PGFLUX_ROOF,     &
                     PRUNOFF_ROOF,                                            &
                     PRN_ROAD, PH_ROAD, PLE_ROAD, PLEW_ROAD, PGFLUX_ROAD,     &
                     PRUNOFF_ROAD,                                            &
                     PRN_WALL, PH_WALL, PLE_WALL, PGFLUX_WALL,                &
                     PRN_BLT,PH_BLT,PLE_BLT, PGFLUX_BLT,                      &
                     PRNSNOW_ROOF, PHSNOW_ROOF, PLESNOW_ROOF, PGSNOW_ROOF,    &
                     PMELT_ROOF,                                              &
                     PRNSNOW_ROAD, PHSNOW_ROAD, PLESNOW_ROAD, PGSNOW_ROAD,    &
                     PMELT_ROAD,                                              &
                     PRUNOFF_TOWN,                                            &
                     PUW_ROAD, PUW_ROOF, PDUWDU_ROAD, PDUWDU_ROOF,            &
                     PUSTAR_TOWN, PCD, PCDN, PCH_TOWN, PRI_TOWN,              &
                     PRESA_TOWN, PDQS_TOWN, PQF_TOWN, PQF_BLD, PQF_BLDWFR,    &
                     PTI_BLD_EQ, PTI_BLDWFR, PFLX_BLD,                        &
                     PAC_ROOF, PAC_ROAD, PAC_WALL, PAC_TOP, PAC_GARDEN,       &
                     PAC_ROOF_WAT, PAC_ROAD_WAT,                              &
                     PABS_SW_ROOF,PABS_LW_ROOF,                               &
                     PABS_SW_SNOW_ROOF, PABS_LW_SNOW_ROOF,                    &
                     PABS_SW_ROAD,PABS_LW_ROAD,                               &
                     PABS_SW_SNOW_ROAD, PABS_LW_SNOW_ROAD,                    &
                     PABS_SW_WALL,PABS_LW_WALL,                               &
                     PLW_W_TO_W  , PLW_W_TO_R  , PLW_W_TO_G ,                 &
                     PLW_W_TO_NR ,                                            &
                     PLW_R_TO_W  , PLW_R_TO_R  , PLW_R_TO_G ,                 &
                     PLW_R_TO_NR ,                                            &
                     PLW_G_TO_W  , PLW_G_TO_R  , PLW_G_TO_G ,                 &
                     PLW_G_TO_NR ,                                            &
                     PLW_S_TO_W  , PLW_S_TO_R  , PLW_S_TO_G ,                 &
                     PLW_S_TO_NR ,                                            &
                     PLW_NR_TO_W , PLW_NR_TO_R , PLW_NR_TO_G,                 &
                     PLW_NR_TO_NR                                             )
!   ##########################################################################
!
!!****  *TEB*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the evoultion of prognostic variables and the fluxes
!     over artificial surfaces as towns, taking into account the canyon like
!     geometry of urbanized areas.
!         
!     
!!**  METHOD
!     ------
!
!     The prognostic variables are:
!       - the surface temperature for roofs, roads, and walls
!       - the water reservoir, whose maximum value is 10mm
!
!
!    1 : Warning about snow
!        ******************
!
!     Except for snow mantel evolution, all other computation with snow
!   variables must be performed with these variables at previous time-step,
!   and NOT new time-step. This insure coherence between snow fractions
!   (computed at the begining) and other snow characteristics (albedo, Ts).
!
!
!    2 : computation of input solar radiation on each surface
!        ****************************************************
!
!    direct fluxes:
!    -------------
!
!    dir_Rg_road (Wm-2) =   S * 2*theta0/pi
!                         - S *2/tan(zen) * h/W /pi * (1-cos(theta0))
!
!    dir_Rg_wall (Wm-2) =   S / tan(zen) /pi * (1-cos(theta0))
!                         + S * W/h * (1/2 -theta0/pi)
!
!   where zen      is the zenithal angle, from horizon
!         h/W      is the aspect ratio of the canyon
!         S        is the direct solar radiation flux on a horizontal surface
!
!         theta0 = arcsin(min(W/h * tan(zen),1))
!
!   The surfaces will keep (1-a) times these fluxes, and reflect the
!   remaining
!
!    scattered fluxes:
!    ----------------
!
!   sca_Rg_road = sca_Rg * SVF_road
!
!   sca_Rg_wall = sca_Rg * SVF_wall
!
!
!    solar flux and isotropic reflections :
!    ------------------------------------
!
!  after 0 reflection, the absorbed part of the flux is:
!
!      ARg_r(0) = (1-a_r) (sca_Rg_road + dir_Rg_road)
!
!      ARg_w(0) = (1-a_w) (sca_Rg_wall + dir_Rg_wall)
!  
!    and the reflected parts are
!
!      RRg_r(0) = a_r (sca_Rg_road + dir_Rg_road)
!
!      RRg_w(0) = a_w (sca_Rg_wall + dir_Rg_wall)
!
!  after n reflection:
!
!      ARg_r(n) = ARg_r(n-1) + RRg_w(n-1) * (1-  SVF_r)(1-a_r)
!
!      ARg_w(n) = ARg_w(n-1) + RRg_r(n-1) *      SVF_w (1-a_w)
!                            + RRg_w(n-1) * (1-2*SVF_w)(1-a_w)
!
!      RRg_r(n) = (1- SVF_r) a_r RRg_w(n-1)
!
!      RRg_w(n) =     SVF_w  a_w RRg_r(n-1)
!                +(1-2SVF_w) a_w RRg_w(n-1)
!
!
!   i.e.
!                                               n-1
!      ARg_r(n) = ARg_r(0) + (1-  SVF_r)(1-a_r) SUM RRg_w(k)
!                                               k=0
!
!                                               n-1
!      ARg_w(n) = ARg_w(0) +      SVF_w (1-a_w) SUM RRg_r(k)
!                                               k=0
!                                               n-1
!                          + (1-2*SVF_w)(1-a_w) SUM RRg_w(k)
!                                               k=0
!
! with
!
!     n                             n-1
!    SUM RRg_r(k) = (1-  SVF_r) a_r SUM RRg_w(k)      +  RRg_r(0)
!    k=0                            k=0
!
!     n                             n-1
!    SUM RRg_w(k) =      SVF_w  a_w SUM RRg_r(k) 
!    k=0                            k=0
!                                   n-1
!                  +(1-2*SVF_w) a_w SUM RRg_w(k)      +  RRg_w(0)
!                                   k=0
!
!
!   Then
!
!     n                                        n-1
!    SUM RRg_w(k) =  (1-2*SVF_w)       a_w     SUM RRg_w(k)
!    k=0                                       k=0
!                                              n-2
!                  + (1-  SVF_r) SVF_w a_w a_r SUM RRg_w(k) 
!                                              k=0
!
!                  + RRg_w(0) + SVF_w a_w RRg_r(0)
!
!
!
!
!  solving this system, lead after an infinity of reflections/absorptions:
!
!    inf                      RRg_w(0) + SVF_w a_w RRg_r(0)
!    SUM RRg_w(k) = ----------------------------------------------------
!    k=0             1 - (1-2*SVF_w) a_w + (1-  SVF_r) SVF_w a_w a_r
!
!
!    inf            (1-  SVF_r) a_r ( a_w SVF_w RRg_r(0) + RRg_w(0) ) + RRg_r(0)
!    SUM RRg_r(k) = ------------------------------------------------------------
!    k=0             1 - (1-2*SVF_w) a_w + (1-  SVF_r) SVF_w a_w a_r
!
!
! ARg_r(n) and ARg_w(n) follow
!
!
!    3 : drag coefficient for momentum 
!        *****************************
!
!
!    4 : aerodynamical resistance for heat transfers
!        *******************************************
!
!
!    5 : equation for evolution of Ts_roof
!        *********************************
!
!
!       Rn = (dir_Rg + sca_Rg) (1-a) + emis * ( Rat - sigma Ts**4 (t+dt) )
!
!       H  = rho Cp CH V ( Ts (t+dt) - Tas )
!
!       LE = rho Lv CH V ( qs (t+dt) - qas )
!
!      where the as subscript denotes atmospheric values at ground level
!      (and not at first half level)
!
!
!    6 : equations for evolution of Ts_road and Ts_wall simultaneously
!        *************************************************************
!
!
!
!   Rn_w = abs_Rg_w 
!  - sigma * emis_w                                                   * Ts_w**4 (t+dt)
!  +         emis_w                       *      SVF_w                * Rat
!  + sigma * emis_w * emis_r              *      SVF_w                * Ts_r**4 (t+dt)
!  + sigma * emis_w * emis_w              * (1-2*SVF_w)               * Ts_w**4 (t+dt)
!  + sigma * emis_w * emis_w * (1-emis_r) *      SVF_w  * (1-  SVF_r) * Ts_w**4 (t+dt)
!  + sigma * emis_w * emis_w * (1-emis_w) * (1-2*SVF_w) * (1-2*SVF_w) * Ts_w**4 (t+dt)
!  + sigma * emis_w * emis_r * (1-emis_w) *      SVF_w  * (1-2*SVF_w) * Ts_r**4 (t+dt)
!
!   Rn_r = abs_Rg_r
!  - sigma * emis_r                                                   * Ts_r**4 (t+dt)
!  +         emis_r                       *    SVF_r                  * Rat
!  + sigma * emis_r * emis_w              * (1-SVF_r)                 * Ts_w**4 (t+dt)
!  + sigma * emis_r * emis_w * (1-emis_w) * (1-SVF_r)   * (1-2*SVF_w) * Ts_w**4 (t+dt)
!  + sigma * emis_r * emis_r * (1-emis_w) * (1-SVF_r)   *      SVF_w  * Ts_r**4 (t+dt)
!
!  H_w  = rho Cp CH V ( Ts_w (t+dt) - Ta_canyon )
!
!  LE_w = rho Lv CH V ( qs_w (t+dt) - qa_canyon )
!
!  H_r  = rho Cp CH V ( Ts_r (t+dt) - Ta_canyon )
!
!  LE_r = rho Lv CH V ( qs_r (t+dt) - qa_canyon )
!
! with again
!                AC_can * Swall/Sroad * Twall + AC_can * Troad + AC_top * Ta + H_traffic/Cp/rho/Sroad
!   Ta_canyon = -------------------------------------------------------------------------------------
!                AC_can * Swall/Sroad         + AC_can         + AC_top
!
!
!                 AC_can * delt_road * Hu_road * qsat(Troad) + AC_top * qa + LE_traffic/Lv/rho/Sroad
!   qa_canyon = ------------------------------------------------------------------------------------
!                 AC_can * delt_road                        + AC_top
!
!
!
!
!    7 : computation of fluxes for each surface type
!        *******************************************
!
!
!    8 : averaging of the fluxes
!        ***********************
!
!   This is done on the total exchange surface (roof + wall + road),
!  which is bigger than the horizontal surface (roof+road), leading
!  to bigger fluxes.
!
!   The fluxes due to industrial activity are directly added into the 
!  atmosphere
!
!
!    9 : road reservoir evolution
!        ************************
!
!   The roof reservoir runoff goes directly into the road reservoir.
!
!   Runoff occurs for road reservoir (too much water), as well as drainage
!   (evacuation system, typical time scale: 1 day)
!
!
!
!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    MODD_CST
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!	V. Masson           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    23/01/98 
!!     21 / 10 / 2003   P. Tulet    output aerodynamical resistance
!!     01 / 07 / 2005   P.Le Moigne Exner functions as arguments to urban_fluxes
!!     17 / 10 / 2005   (G. Pigeon) computation of anthropogenic heat from domestic heating
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,         ONLY : XTT, XSTEFAN, XCPD, XLVTT
USE MODD_SURF_PAR,     ONLY : XUNDEF
USE MODD_SNOW_PAR,     ONLY : XEMISSN, XANSMAX_ROOF, &
                          XANSMAX_ROAD,XWCRN_ROOF,XWCRN_ROAD
!
USE MODE_THERMOS
USE MODE_SURF_SNOW_FRAC
!
USE MODI_URBAN_DRAG
USE MODI_URBAN_SNOW_EVOL
USE MODI_ROOF_LAYER_E_BUDGET
USE MODI_ROAD_LAYER_E_BUDGET
USE MODI_WALL_LAYER_E_BUDGET
USE MODI_URBAN_FLUXES
USE MODI_URBAN_HYDRO
USE MODI_BLD_E_BUDGET
USE MODI_WIND_THRESHOLD
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    Declarations of arguments
!
!
CHARACTER(LEN=6), INTENT(IN)      :: HZ0H          ! TEB option for z0h roof & road
!                                                  ! 'MASC95' : Mascart et al 1995
!                                                  ! 'BRUT82' : Brustaert     1982
!                                                  ! 'KAND07' : Kanda         2007
CHARACTER(LEN=*)  , INTENT(IN)    :: HCOUPLING     ! type of coupling
                                                   ! 'E' : explicit
                                                   ! 'I' : implicit
REAL, DIMENSION(:), INTENT(INOUT) :: PT_CANYON     ! canyon air temperature
REAL, DIMENSION(:), INTENT(INOUT) :: PQ_CANYON     ! canyon air specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PU_CANYON     ! canyon hor. wind
REAL, DIMENSION(:), INTENT(IN)    :: PU_LOWCAN     ! wind near the road
REAL, DIMENSION(:), INTENT(IN)    :: PT_LOWCAN     ! temp. near the road
REAL, DIMENSION(:), INTENT(IN)    :: PQ_LOWCAN     ! hum. near the road
REAL, DIMENSION(:), INTENT(IN)    :: PZ_LOWCAN     ! height of atm. var. near the road
REAL, DIMENSION(:), INTENT(INOUT) :: PTI_BLD       ! inside building temperature
REAL, DIMENSION(:), INTENT(INOUT) :: PTI_BLD_EQ    ! inside building temperature
                                                   ! computed by its own evolution equation
REAL, DIMENSION(:), INTENT(INOUT) :: PTI_BLDWFR    ! inside building temperature
                                                   ! computed by its own evolution equation
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_ROOF     ! roof layers temperatures
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_ROAD     ! road layers temperatures
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_WALL     ! wall layers temperatures
REAL, DIMENSION(:), INTENT(INOUT) :: PWS_ROOF      ! roof water reservoir
REAL, DIMENSION(:), INTENT(INOUT) :: PWS_ROAD      ! road water reservoir
CHARACTER(LEN=*),   INTENT(IN)    :: HSNOW_ROOF    ! snow roof scheme
!                                                  ! 'NONE'
!                                                  ! 'D95 '
!                                                  ! '1-L '
CHARACTER(LEN=*),   INTENT(IN)    :: HSNOW_ROAD    ! snow road scheme
!                                                  ! 'NONE'
!                                                  ! 'D95 '
!                                                  ! '1-L '
REAL, DIMENSION(:,:), INTENT(INOUT) :: PWSNOW_ROOF ! snow layers reservoir
REAL, DIMENSION(:,:), INTENT(INOUT) :: PTSNOW_ROOF ! snow layers temperature
REAL, DIMENSION(:,:), INTENT(INOUT) :: PRSNOW_ROOF ! snow layers density
REAL, DIMENSION(:),   INTENT(INOUT) :: PASNOW_ROOF ! snow albedo
REAL, DIMENSION(:),   INTENT(INOUT) :: PESNOW_ROOF ! snow emissivity
REAL, DIMENSION(:),   INTENT(INOUT) :: PTSSNOW_ROOF! snow surface temperature
REAL, DIMENSION(:,:), INTENT(INOUT) :: PWSNOW_ROAD ! snow layers reservoir
REAL, DIMENSION(:,:), INTENT(INOUT) :: PTSNOW_ROAD ! snow layers temperature
REAL, DIMENSION(:,:), INTENT(INOUT) :: PRSNOW_ROAD ! snow layers density
REAL, DIMENSION(:),   INTENT(INOUT) :: PASNOW_ROAD ! snow albedo
REAL, DIMENSION(:),   INTENT(INOUT) :: PESNOW_ROAD ! snow emissivity
REAL, DIMENSION(:),   INTENT(INOUT) :: PTSSNOW_ROAD! snow surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_A_COEF   ! implicit coefficients
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_B_COEF   ! for wind coupling
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_A_COEF_LOWCAN ! implicit coefficients for wind coupling
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_B_COEF_LOWCAN ! between low canyon wind and road
REAL, DIMENSION(:), INTENT(IN)    :: PPS           ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PPA           ! pressure at the first atmospheric level
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS         ! surface exner function
REAL, DIMENSION(:), INTENT(IN)    :: PTA           ! temperature at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PQA           ! specific humidity
                                                   ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD         ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PEXNA         ! exner function
                                                   ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA         ! air density
                                                   ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD       ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PDIR_SW_RAD   ! incoming direct solar radiation
                                                   ! on an horizontal surface
REAL, DIMENSION(:), INTENT(IN)    :: PZENITH       ! solar zenithal angle
REAL, DIMENSION(:), INTENT(IN)    :: PSCA_SW_RAD   ! scattered incoming solar rad.
REAL, DIMENSION(:), INTENT(IN)    :: PRR           ! rain rate
REAL, DIMENSION(:), INTENT(IN)    :: PSR           ! snow rate
REAL, DIMENSION(:), INTENT(IN)    :: PH_TRAFFIC    ! anthropogenic sensible
!                                                  ! heat fluxes due to traffic
REAL, DIMENSION(:), INTENT(IN)    :: PLE_TRAFFIC   ! anthropogenic latent
!                                                  ! heat fluxes due to traffic
REAL, DIMENSION(:), INTENT(IN)    :: PH_INDUSTRY   ! anthropogenic sensible
!                                                  ! heat fluxes due to factories
REAL, DIMENSION(:), INTENT(IN)    :: PLE_INDUSTRY  ! anthropogenic latent
!                                                  ! heat fluxes due to factories
REAL, DIMENSION(:), INTENT(IN)    :: PZREF         ! reference height of the first
                                                   ! atmospheric level (temperature)
REAL, DIMENSION(:), INTENT(IN)    :: PUREF         ! reference height of the first
                                                   ! atmospheric level (wind)
REAL,               INTENT(IN)    :: PTSTEP        ! time step
REAL, DIMENSION(:), INTENT(IN)    :: PZ0_TOWN      ! town roughness length
                                                   ! for momentum
REAL, DIMENSION(:), INTENT(IN)    :: PBLD          ! fraction of buildings
REAL, DIMENSION(:), INTENT(IN)    :: PGARDEN       ! fraction of GARDEN areas
REAL, DIMENSION(:), INTENT(IN)    :: PROAD         ! fraction of roads
REAL, DIMENSION(:), INTENT(IN)    :: PBLD_HEIGHT   ! buildings h
REAL, DIMENSION(:), INTENT(IN)    :: PWALL_O_HOR   ! wall surf. / hor. surf.
REAL, DIMENSION(:), INTENT(IN)    :: PCAN_HW_RATIO ! canyon    h/W
!
REAL, DIMENSION(:), INTENT(INOUT) :: PDF_ROOF      ! snow-free    fraction on roofs
REAL, DIMENSION(:), INTENT(INOUT) :: PDN_ROOF      ! snow-covered fraction on roofs
REAL, DIMENSION(:), INTENT(INOUT) :: PDF_ROAD      ! snow-free    fraction on roads
REAL, DIMENSION(:), INTENT(INOUT) :: PDN_ROAD      ! snow-covered fraction on roads
REAL, DIMENSION(:), INTENT(OUT)   :: PQSAT_ROOF    ! hum at saturation over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PQSAT_ROAD    ! hum at saturation over road
REAL, DIMENSION(:), INTENT(OUT)   :: PDELT_ROOF    ! water fraction on roof
REAL, DIMENSION(:), INTENT(OUT)   :: PDELT_ROAD    ! water fraction on road
!
REAL, DIMENSION(:), INTENT(IN)    :: PALB_ROOF     ! roof albedo
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_ROOF    ! roof emissivity
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_ROOF      ! heat capacity for roof layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_ROOF      ! thermal conductivity for roof layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_ROOF       ! depth of roof layers
REAL, DIMENSION(:), INTENT(IN)    :: PALB_ROAD     ! road albedo
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_ROAD    ! road emissivity
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_ROAD      ! heat capacity for road layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_ROAD      ! thermal conductivity for road layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_ROAD       ! depth of road layers
REAL, DIMENSION(:), INTENT(IN)    :: PSVF_ROAD     ! road sky view factor
REAL, DIMENSION(:), INTENT(IN)    :: PALB_WALL     ! wall albedo
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_WALL    ! wall emissivity
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_WALL      ! heat capacity for wall layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_WALL      ! thermal conductivity for wall layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_WALL       ! depth of wall layers
REAL, DIMENSION(:), INTENT(IN)    :: PSVF_WALL     ! wall sky view factor
REAL, DIMENSION(:), INTENT(IN)    :: PTS_GARDEN     ! GARDEN area surf temp.

REAL, DIMENSION(:), INTENT(OUT)   :: PRN_ROOF     ! net radiation over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PH_ROOF      ! sensible heat flux over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_ROOF     ! latent heat flux over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PLEW_ROOF    ! latent heat flux over roof (snow)
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_ROOF  ! flux through the roof
REAL, DIMENSION(:), INTENT(OUT)   :: PRUNOFF_ROOF ! runoff over the ground
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_ROAD     ! net radiation over road
REAL, DIMENSION(:), INTENT(OUT)   :: PH_ROAD      ! sensible heat flux over road
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_ROAD     ! latent heat flux over road
REAL, DIMENSION(:), INTENT(OUT)   :: PLEW_ROAD    ! latent heat flux over road (snow)
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_ROAD  ! flux through the road
REAL, DIMENSION(:), INTENT(OUT)   :: PRUNOFF_ROAD ! runoff over the ground
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_WALL     ! net radiation over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PH_WALL      ! sensible heat flux over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_WALL     ! latent heat flux over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_WALL  ! flux through the wall
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_BLT       ! net radiation over built surf 
REAL, DIMENSION(:), INTENT(OUT)   :: PH_BLT        ! sensible heat flux over built surf 
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_BLT       ! latent heat flux over built surf 
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_BLT    ! flux through the built surf 
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRNSNOW_ROOF ! net radiation over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PHSNOW_ROOF  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PLESNOW_ROOF ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PGSNOW_ROOF  ! flux under the snow
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT_ROOF   ! snow melt
REAL, DIMENSION(:), INTENT(OUT)   :: PRNSNOW_ROAD ! net radiation over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PHSNOW_ROAD  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PLESNOW_ROAD ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(OUT)   :: PGSNOW_ROAD  ! flux under the snow
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT_ROAD   ! snow melt
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRUNOFF_TOWN ! runoff over the ground
REAL, DIMENSION(:), INTENT(OUT)   :: PUW_ROAD     ! Momentum flux for roads
REAL, DIMENSION(:), INTENT(OUT)   :: PUW_ROOF     ! Momentum flux for roofs
REAL, DIMENSION(:), INTENT(OUT)   :: PDUWDU_ROAD  !
REAL, DIMENSION(:), INTENT(OUT)   :: PDUWDU_ROOF  !
REAL, DIMENSION(:), INTENT(OUT)   :: PUSTAR_TOWN  ! friciton velocity over town
REAL, DIMENSION(:), INTENT(OUT)   :: PCD          ! town averaged drag coefficient
REAL, DIMENSION(:), INTENT(OUT)   :: PCDN         ! town averaged neutral drag coefficient
REAL, DIMENSION(:), INTENT(OUT)   :: PCH_TOWN     ! town averaged heat transfer
!                                                 ! coefficient
REAL, DIMENSION(:), INTENT(OUT)   :: PRI_TOWN      ! town averaged Richardson number
REAL, DIMENSION(:), INTENT(OUT)   :: PRESA_TOWN    ! town aerodynamical resistance
REAL, DIMENSION(:), INTENT(OUT)   :: PDQS_TOWN     ! heat storage inside town
REAL, DIMENSION(:), INTENT(OUT)   :: PQF_TOWN      ! total anthropogenic heat
REAL, DIMENSION(:), INTENT(OUT)   :: PQF_BLD       ! anthropogenic heat flux of
                                                   ! domestic heating  
REAL, DIMENSION(:), INTENT(OUT)   :: PQF_BLDWFR 
REAL, DIMENSION(:), INTENT(OUT)   :: PFLX_BLD      ! heat flx from inside bld through its structure
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_ROOF      ! roof conductance
REAL, DIMENSION(:), INTENT(INOUT) :: PAC_ROAD      ! road conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_WALL      ! wall conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_TOP       ! top conductance
REAL, DIMENSION(:), INTENT(IN)    :: PAC_GARDEN    ! garden conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_ROOF_WAT  ! roof water conductance
REAL, DIMENSION(:), INTENT(OUT)   :: PAC_ROAD_WAT  ! roof water conductance
!
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_ROOF      ! absorbed solar rad by roof
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_SNOW_ROOF ! absorbed solar rad by snow on roof
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_ROAD      ! absorbed solar rad by road
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_SNOW_ROAD ! absorbed solar rad by snow on road
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_WALL      ! absorbed solar rad by wall
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_ROOF      ! absorbed IR rad by roof
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_SNOW_ROOF ! absorbed IR rad by snow on roof
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_ROAD      ! absorbed IR rad by road
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_SNOW_ROAD ! absorbed IR rad by snow on road
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_WALL      ! absorbed IR rad by wall
!
REAL, DIMENSION(:), INTENT(IN)    :: PLW_W_TO_W          ! LW contrib. wall       -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_W_TO_R          ! LW contrib. wall       -> road
REAL, DIMENSION(:), INTENT(IN)    :: PLW_W_TO_G          ! LW contrib. wall       -> GARDEN
REAL, DIMENSION(:), INTENT(IN)    :: PLW_W_TO_NR         ! LW contrib. wall       -> road(snow)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_R_TO_W          ! LW contrib. road       -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_R_TO_R          ! LW contrib. road       -> road
REAL, DIMENSION(:), INTENT(IN)    :: PLW_R_TO_G          ! LW contrib. road       -> GARDEN
REAL, DIMENSION(:), INTENT(IN)    :: PLW_R_TO_NR         ! LW contrib. road       -> road(snow)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_G_TO_W          ! LW contrib. GARDEN     -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_G_TO_R          ! LW contrib. GARDEN     -> road
REAL, DIMENSION(:), INTENT(IN)    :: PLW_G_TO_G          ! LW contrib. GARDEN     -> GARDEN
REAL, DIMENSION(:), INTENT(IN)    :: PLW_G_TO_NR         ! LW contrib. GARDEN     -> road(snow)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_NR_TO_W         ! LW contrib. road(snow) -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_NR_TO_R         ! LW contrib. road(snow) -> road
REAL, DIMENSION(:), INTENT(IN)    :: PLW_NR_TO_G         ! LW contrib. road(snow) -> GARDEN
REAL, DIMENSION(:), INTENT(IN)    :: PLW_NR_TO_NR        ! LW contrib. road(snow) -> road(snow)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_W          ! LW contrib. sky        -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_R          ! LW contrib. sky        -> road
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_G          ! LW contrib. sky        -> GARDEN
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_NR         ! LW contrib. sky        -> road(snow)
!
!*      0.2    Declarations of local variables
!
LOGICAL                    :: GTI_EVOL       ! true --> internal temperature
!                                            !          of buildings evolves
!                                            ! false--> it is fixed
!
REAL, DIMENSION(SIZE(PTA)) :: ZVMOD          ! wind
REAL, DIMENSION(SIZE(PTA)) :: ZWS_ROOF_MAX   ! maximum deepness of roof
REAL, DIMENSION(SIZE(PTA)) :: ZWS_ROAD_MAX   ! and road water reservoirs
!
REAL, DIMENSION(SIZE(PTA)) :: ZAC_BLD        ! surface conductance inside
                                             ! the building itself
REAL, DIMENSION(SIZE(PTA)) :: ZTA            ! air temperature extrapolated at roof level
REAL, DIMENSION(SIZE(PTA)) :: ZQA            ! air humidity extrapolated at roof level
!
REAL, DIMENSION(SIZE(PTA)) :: ZROOF_FRAC     ! roof, wall and
REAL, DIMENSION(SIZE(PTA)) :: ZWALL_FRAC     ! road fractions
REAL, DIMENSION(SIZE(PTA)) :: ZROAD_FRAC     ! of exchange surf.
REAL, DIMENSION(SIZE(PTA)) :: ZGARDEN_FRAC   !                  
REAL, DIMENSION(SIZE(PTA)) :: ZTOTS_O_HORS   ! total canyon+roof surface
!                                            ! over horizontal surface
REAL, DIMENSION(SIZE(PTA)) :: ZWALL_O_ROAD   ! wall surface over road surface
REAL, DIMENSION(SIZE(PTA)) :: ZWALL_O_GRND   ! wall surface over (road+GARDEN area) surface
REAL, DIMENSION(SIZE(PTA)) :: ZDQS_ROAD      ! heat storage inside road
REAL, DIMENSION(SIZE(PTA)) :: ZDQS_ROOF      ! heat storage inside roof
REAL, DIMENSION(SIZE(PTA)) :: ZDQS_WALL      ! heat storage inside wall
REAL, DIMENSION(SIZE(PTA)) :: ZQF_WALL       ! 
REAL, DIMENSION(SIZE(PTA)) :: ZQF_ROOF       !heating through the roof 
REAL, DIMENSION(SIZE(PTA)) :: ZFLX_BLD_ROOF  !heat flux from inside through roof
REAL, DIMENSION(SIZE(PTA)) :: ZFLX_BLD_WALL  !heat flux from inside through wall
!
REAL, DIMENSION(SIZE(PTA)) :: ZDQS_SNOW_ROOF ! heat storage inside roof snowpack
REAL, DIMENSION(SIZE(PTA)) :: ZDQS_SNOW_ROAD ! heat storage inside road snowpack
REAL, DIMENSION(SIZE(PTA)) :: ZMELT_BLT      ! Snow melt for built & impervious part
!
! coefficients for LW computations over snow (from previous time-step)
!
REAL, DIMENSION(SIZE(PTA)) :: ZTSSNOW_ROAD ! road snow temperature
!                                          ! at previous time-step
REAL, DIMENSION(SIZE(PTA)) :: ZTS_ROAD     ! road surface temperature 
!                                          ! at previous time-step
REAL, DIMENSION(SIZE(PTA)) :: ZTS_WALL     ! wall surface temperature 
!                                          ! at previous time-step
REAL, DIMENSION(SIZE(PTA)) :: ZTS_ROOF     ! roof surface temperature 
!                                          ! at previous time-step
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!*      1.     Initializations
!              ---------------
!
!*      1.1    Water reservoirs
!              ----------------
!
IF (LHOOK) CALL DR_HOOK('TEB',0,ZHOOK_HANDLE)
ZWS_ROOF_MAX =  1. ! (1mm) maximum deepness of roof water reservoir
ZWS_ROAD_MAX =  1. ! (1mm) maximum deepness of road water reservoir
!
!*      1.2    surfaces relative fractions
!              ---------------------------
!
ZTOTS_O_HORS = 1. + PWALL_O_HOR(:)
ZROOF_FRAC   = PBLD        / ZTOTS_O_HORS
ZWALL_FRAC   = PWALL_O_HOR / ZTOTS_O_HORS
ZROAD_FRAC   = PROAD       / ZTOTS_O_HORS
ZGARDEN_FRAC = PGARDEN     / ZTOTS_O_HORS
WHERE (ZROAD_FRAC(:)+ZGARDEN_FRAC(:).GT. 0.) 
!
  ZWALL_O_GRND = ZWALL_FRAC  / (ZROAD_FRAC+ZGARDEN_FRAC)
!
ELSEWHERE
!
  ZWALL_O_GRND = 0.
!
END WHERE
!
WHERE (PROAD(:) .GT. 0.)
!
 ZWALL_O_ROAD = ZWALL_FRAC /  ZROAD_FRAC
!
ELSEWHERE (PBLD(:) .EQ. 0. .AND. PROAD(:) .EQ. 0.)
!
 ZWALL_O_ROAD = 0.
!
ENDWHERE
!
!*      1.3    radiative temperatures variables at previous time-step
!              ------------------------------------------------------
!
ZTSSNOW_ROAD(:)=PTSSNOW_ROAD(:)
ZTS_ROAD    (:)=PT_ROAD     (:,1)
ZTS_WALL    (:)=PT_WALL     (:,1)
ZTS_ROOF    (:)=PT_ROOF     (:,1)
!
!*       1.4   Option for internal temperature of buildings
!              --------------------------------------------
!
GTI_EVOL = .TRUE.
!
!-------------------------------------------------------------------------------
!
!*      2.     Snow-covered surfaces relative effects
!              --------------------------------------
!
!*      2.1    Effects on water reservoirs
!              ---------------------------
!
ZWS_ROOF_MAX(:) = ZWS_ROOF_MAX(:) * PDF_ROOF(:)
ZWS_ROAD_MAX(:) = ZWS_ROAD_MAX(:) * PDF_ROAD(:)
!
!-------------------------------------------------------------------------------
!
!*      3.     Surface drag
!              ------------
!
CALL URBAN_DRAG(HZ0H, PTSTEP, PT_CANYON, PQ_CANYON, PU_CANYON,       &
                PT_LOWCAN, PQ_LOWCAN, PU_LOWCAN, PZ_LOWCAN,          &
                ZTS_ROOF, ZTS_ROAD, ZTS_WALL,                        &
                PTS_GARDEN, PDN_ROOF, PDN_ROAD,                      &
                PEXNS, PEXNA, PTA, PQA, PPS, PRHOA,                  &
                PZREF, PUREF, PVMOD,                                 &
                PZ0_TOWN,                                            &
                PBLD, PGARDEN, PROAD,                                &
                PBLD_HEIGHT, PCAN_HW_RATIO,                          &
                ZWALL_O_GRND,                                        &
                PWS_ROOF, PWS_ROAD,                                  &
                ZWS_ROOF_MAX, ZWS_ROAD_MAX,                          &
                PPEW_A_COEF, PPEW_B_COEF,                            &
                PPEW_A_COEF_LOWCAN, PPEW_B_COEF_LOWCAN,              &
                PQSAT_ROOF, PQSAT_ROAD, PDELT_ROOF, PDELT_ROAD,      &
                PCD, PCDN, PAC_ROOF, PAC_ROOF_WAT,                   &
                PAC_WALL, PAC_ROAD, PAC_ROAD_WAT, PAC_TOP,           &
                PAC_GARDEN, PRI_TOWN,                                &
                PUW_ROAD, PUW_ROOF, PDUWDU_ROAD, PDUWDU_ROOF,        &
                PUSTAR_TOWN                                          )
!
!* area-averaged heat transfer coefficient
!
ZVMOD(:) = WIND_THRESHOLD(PVMOD(:),PUREF(:))
!
PCH_TOWN(:) = (PBLD(:) * PAC_ROOF(:) + (1.-PBLD(:)) * PAC_TOP (:)) / ZVMOD(:)
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
!*      5.     Evolution of interior building air temperature
!              ----------------------------------------------
!
CALL BLD_E_BUDGET(GTI_EVOL, PTSTEP, PBLD, PWALL_O_HOR,             &
                  PRHOA, PT_ROOF, PT_WALL, PTI_BLD, ZAC_BLD,       &
                  PTI_BLD_EQ, PTI_BLDWFR                           )
!
!-------------------------------------------------------------------------------
!
!*      6.     Snow mantel model
!              -----------------
!
CALL URBAN_SNOW_EVOL(PT_LOWCAN, PQ_LOWCAN, PU_LOWCAN,                         &
                     ZTS_ROOF, ZTS_ROAD, ZTS_WALL, PTS_GARDEN,                &
                     HSNOW_ROOF,                                              &
                     PWSNOW_ROOF, PTSNOW_ROOF, PRSNOW_ROOF, PASNOW_ROOF,      &
                     PTSSNOW_ROOF, PESNOW_ROOF,                               &
                     HSNOW_ROAD,                                              &
                     PWSNOW_ROAD, PTSNOW_ROAD, PRSNOW_ROAD, PASNOW_ROAD,      &
                     PTSSNOW_ROAD, PESNOW_ROAD,                               &
                     PPS, ZTA, ZQA, PRHOA,                                    &
                     PLW_RAD,                                                 &
                     PSR, PZREF, PUREF, PVMOD,                                &
                     PTSTEP,                                                  &
                     PZ_LOWCAN,                                               &
                     PDN_ROOF, PABS_SW_SNOW_ROOF, PABS_LW_SNOW_ROOF,          &
                     PDN_ROAD, PABS_SW_SNOW_ROAD, PABS_LW_SNOW_ROAD,          &
                     PRNSNOW_ROOF, PHSNOW_ROOF, PLESNOW_ROOF, PGSNOW_ROOF,    &
                     PMELT_ROOF,                                              &
                     PRNSNOW_ROAD, PHSNOW_ROAD, PLESNOW_ROAD, PGSNOW_ROAD,    &
                     PMELT_ROAD,                                              &
                     PLW_W_TO_NR , PLW_R_TO_NR , PLW_G_TO_NR ,                &
                     PLW_NR_TO_NR, PLW_S_TO_NR,                               &
                     ZDQS_SNOW_ROOF, ZDQS_SNOW_ROAD                           )
!
!-------------------------------------------------------------------------------
!
!*      10.    LW properties
!              -------------
!
PDF_ROAD (:) = 1. - PDN_ROAD (:)
!
!-------------------------------------------------------------------------------
!
!*      12.    Roof Ts computation
!              -------------------
!
!* ts_roof and qsat_roof are updated
!
CALL ROOF_LAYER_E_BUDGET(PT_ROOF, PQSAT_ROOF,                              &
                         ZTA, ZQA, PPS,                                    &
                         PLW_RAD, PTSTEP,                                  &
                         PEMIS_ROOF, PHC_ROOF, PTC_ROOF, PD_ROOF,          &
                         PTI_BLD, PTI_BLDWFR, ZAC_BLD, PDELT_ROOF,         &
                         PDN_ROOF, PGSNOW_ROOF,                            &
                         PRHOA, PAC_ROOF, PAC_ROOF_WAT,                    &
                         PABS_SW_ROOF, PABS_LW_ROOF,                       &
                         ZDQS_ROOF,ZQF_ROOF, ZFLX_BLD_ROOF                 )
!
!-------------------------------------------------------------------------------
!
!*      13.    Wall temperatures computations
!              ------------------------------
!
CALL WALL_LAYER_E_BUDGET(ZTS_ROAD, PT_WALL, PT_CANYON, PQ_CANYON,      &
                              PTS_GARDEN,                              &
                              ZTA, ZQA, PPS,                           &
                              PLW_RAD,  PTSTEP,                        &
                              PHC_WALL, PTC_WALL, PD_WALL,             &
                              PTI_BLD, PTI_BLDWFR, ZAC_BLD,            &
                              ZTSSNOW_ROAD, PRHOA, PAC_WALL,           &
                              PABS_SW_WALL, PABS_LW_WALL,              &
                              PLW_W_TO_W, PLW_R_TO_W, PLW_G_TO_W,      &
                              PLW_S_TO_W, PLW_NR_TO_W,                 &
                              ZDQS_WALL, ZQF_WALL, ZFLX_BLD_WALL       )
!
!-------------------------------------------------------------------------------
!
!*      14.    Road temperatures computations
!              ------------------------------
!
CALL ROAD_LAYER_E_BUDGET(PT_ROAD, ZTS_WALL, PQSAT_ROAD,                &
                              PT_LOWCAN, PQ_LOWCAN,                    &
                              PTS_GARDEN,                              &
                              ZTA, ZQA, PPS,                           &
                              PLW_RAD,  PTSTEP,                        &
                              PHC_ROAD, PTC_ROAD, PD_ROAD,             &
                              PDELT_ROAD, PDN_ROAD,                    &
                              ZTSSNOW_ROAD, PGSNOW_ROAD,               &
                              PRHOA,                                   &
                              PAC_ROAD, PAC_ROAD_WAT,                  &
                              PABS_SW_ROAD, PABS_LW_ROAD,              &
                              PLW_W_TO_R, PLW_R_TO_R, PLW_G_TO_R,      &
                              PLW_S_TO_R, PLW_NR_TO_R, ZDQS_ROAD       )
!
!-------------------------------------------------------------------------------
!
!*      15.    Fluxes over built surfaces
!              --------------------------
!
CALL URBAN_FLUXES   (PT_CANYON, PQ_CANYON,                                    &
                     PT_LOWCAN, PQ_LOWCAN,                                    &
                     PT_ROOF(:,1),PT_ROAD(:,1),PT_WALL(:,1),                  &
                     PPEW_A_COEF, PPEW_B_COEF,                                &
                     PEXNS, PEXNA, ZTA, ZQA, PRHOA,                           &
                     PVMOD,                                                   &
                     PH_TRAFFIC, PLE_TRAFFIC, PH_INDUSTRY, PLE_INDUSTRY,      &
                     PBLD, PROAD, PWALL_O_HOR,                                &
                     PEMIS_ROOF,                                              &
                     PABS_SW_ROOF, PABS_LW_ROOF, PABS_SW_WALL, PABS_LW_WALL,  &
                     PABS_SW_ROAD, PABS_LW_ROAD,                              &
                     PAC_ROOF, PAC_ROOF_WAT,                                  &
                     PAC_WALL, PAC_ROAD, PAC_ROAD_WAT,                        &
                     PCD,                                                     &
                     PQSAT_ROOF, PQSAT_ROAD,                                  &
                     PDELT_ROOF, PDELT_ROAD,                                  &
                     ZROOF_FRAC, ZWALL_FRAC, ZROAD_FRAC,                      &
                     ZTOTS_O_HORS,                                            &
                     PDF_ROOF, PDN_ROOF, PDF_ROAD, PDN_ROAD,                  &
                     PRN_ROOF, PH_ROOF, PLE_ROOF, PLEW_ROOF, PGFLUX_ROOF,     &
                     PRN_ROAD, PH_ROAD, PLE_ROAD, PLEW_ROAD, PGFLUX_ROAD,     &
                     PRN_WALL, PH_WALL, PLE_WALL, PGFLUX_WALL,                &
                     PRNSNOW_ROOF, PHSNOW_ROOF, PLESNOW_ROOF, PGSNOW_ROOF,    &
                     PRNSNOW_ROAD, PHSNOW_ROAD, PLESNOW_ROAD, PGSNOW_ROAD,    &
                     PRN_BLT,  PH_BLT,  PLE_BLT,  PGFLUX_BLT,                 &
                     PUSTAR_TOWN, ZDQS_WALL, ZDQS_ROAD, ZDQS_ROOF, PQF_BLD,   &
                     PQF_BLDWFR, PDQS_TOWN, PQF_TOWN, ZQF_WALL, ZQF_ROOF,     &
                     ZFLX_BLD_WALL, ZFLX_BLD_ROOF, PFLX_BLD, PMELT_ROOF,      &
                     ZDQS_SNOW_ROOF, PMELT_ROAD, ZMELT_BLT                    )
!
! Water transfer from snow reservoir to water reservoir in case of snow melt
!
WHERE (PMELT_ROOF(:) .GT. 0.)
 PWS_ROOF(:) = MIN(ZWS_ROOF_MAX,PWS_ROOF(:) + PMELT_ROOF(:)*PTSTEP)
ENDWHERE
!
WHERE (PMELT_ROAD(:) .GT. 0.)
 PWS_ROAD(:) = MIN(ZWS_ROAD_MAX,PWS_ROAD(:) + PMELT_ROAD(:)*PTSTEP)
ENDWHERE
!
!-------------------------------------------------------------------------------
!
!*      16.    Roof ans road reservoirs evolution
!              ----------------------------------
!
CALL URBAN_HYDRO(ZWS_ROOF_MAX,ZWS_ROAD_MAX, PWS_ROOF, PWS_ROAD,        &
                 PRR, PTSTEP, PBLD, PLE_ROOF, PLE_ROAD,                &
                 PRUNOFF_ROOF,                                         &
                 PRUNOFF_ROAD,                                         &
                 PRUNOFF_TOWN                                          )
!
!-------------------------------------------------------------------------------
!
!*      17.    Compute aerodynamical resistance 
!              --------------------------------
!
PRESA_TOWN(:) = 1. / ( PBLD(:) * PAC_ROOF(:)  + ( 1. - PBLD(:)) * PAC_TOP (:))
IF (LHOOK) CALL DR_HOOK('TEB',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE TEB
