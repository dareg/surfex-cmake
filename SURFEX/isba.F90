!     #########
      SUBROUTINE ISBA( HISBA, HPHOTO, HRUNOFF, HC1DRY, HSCOND, HSOILFRZ,         &
                         HDIFSFCOND, HSNOW_ISBA, HSNOWRES, HCPSURF,              &
                         HDIF, TPTIME, OFLOOD, OTEMP_ARP, OGLACIER,              &
                         PPEW_A_COEF, PPEW_B_COEF,                               &
                         PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,     &
                         PTSTEP, PZREF, PUREF, PDIRCOSZW,                        &
                         PTA, PQA, PEXNA, PRHOA, PPS, PEXNS,                     &
                         PRR, PSR, PZENITH, PSW_RAD, PLW_RAD,                    &
                         PVMOD,                                                  &
                         PVEG,  PLAI, PWRMAX_CF, PRSMIN, PGAMMA, PCV, PRGL,      &
                         PRUNOFFD, PEMIS, PALB, PZ0_WITH_SNOW, PZ0H_WITH_SNOW,   &
                         PZ0EFF,                                                 &
                         PRUNOFFB, PWDRAIN,                                      &
                         PCGSAT, PC1SAT, PC2REF, PC3,                            &
                         PC4B, PC4REF,                                           &
                         PACOEF, PPCOEF,                                         &
                         PTAUICE,                                                &
                         PCG, PC1, PC2, PWGEQ, PCT, PRS,                         &
                         PTDEEP, PGAMMAT,                                        &
                         PWR,                                                    &
                         PRESA, PCH, PCD, PCDN, PRI, PHU, PHUG,                  &
                         PPSN, PPSNG, PPSNV,                                     &
                         PRN, PH, PLE, PLEI, PLEG, PLEGI, PLEV, PLES, PLER,      &
                         PLETR, PEVAP, PGFLUX, PRESTORE,                         &
                         PUSTAR, PDRAIN, PRUNOFF, PMELT,                         &
                         PVEGTYPE, PAN, PANF,                                    &
                         PANFM, PANDAY, PRESP_BIOMASS_INST, PRESP_BIOMASS,       &
                         PABC, PPOI, PCSP,                                       &
                         PFZERO, PEPSO, PGAMM, PQDGAMM, PGMES, PQDGMES,          &
                         PT1GMES, PT2GMES, PAMAX, PQDAMAX, PT1AMAX, PT2AMAX,     &
                         OSTRESSDEF, PF2I, PGC, PAH, PBH, PDMAX,                 &
                         PSNOWSWE, PSNOWHEAT, PSNOWRHO, PSNOWALB,                &
                         PSNOWGRAN1, PSNOWGRAN2, PSNOWHIST,PSNOWAGE,             &
                         PRNSNOW, PHSNOW,  PHPSNOW, PSMELTFLUX,                  &
                         PGFLUXSNOW, PUSTARSNOW,                                 &
                         PGRNDFLUX, PSRSFC, PRRSFC, PLESL,                       &
                         PEMISNOW, PCDSNOW, PCHSNOW,                             &
                         PSNOWTEMP, PSNOWLIQ, PSNOWDZ, PTS_RAD, PTS, PSNOWHMASS, &
                         PRN_ISBA, PH_ISBA, PLEG_ISBA, PLEGI_ISBA, PLEV_ISBA,    &
                         PLETR_ISBA, PUSTAR_ISBA, PLER_ISBA, PLE_ISBA,           &
                         PLEI_ISBA, PGFLUX_ISBA, PMELTADV,                       &
                         PD_G, PROOTFRAC,                                        &
                         PTG, PWG, PWGI,                                         &
                         PCPS, PLVTT, PLSTT,                                     &
                         PWFC, PWWILT, PWSAT, PBCOEF, PCONDSAT, PMPOTSAT,        &
                         PHCAPSOIL, PCONDDRY, PCONDSLD, PEMIST, PALBT, PIACAN,   &
                         PHV, PSNOWFREE_ALB_VEG, PSNOWFREE_ALB_SOIL,             &
                         PPSNV_A, PQS,                                           &
                         PIRRIG, PWATSUP, PTHRESHOLD, LIRRIGATE, LIRRIDAY,       &
                         PCGMAX, HKSAT, HRAIN, HHORT, PMUF, PFSAT, PKSAT_ICE,    &
                         PD_ICE, PHORT, PDRIP, PGPP, PFFLOOD, PPIFLOOD,          &
                         PIFLOOD, PPFLOOD, PLE_FLOOD, PLEI_FLOOD, PFFG, PFFV,    &
                         PFF, PFALB, PFEMIS, PFFG_NOSNOW, PFFV_NOSNOW, PFFROZEN, &
                         PRRVEG, PCONDSAT_EXP, PEXP_DIF, PALPHA, PN, PM, PL,     &
                         PWRES, PAC_AGG, PHU_AGG, PSODELX, PLAT, PLON            )                         
!     ##########################################################################
!
!
!!****  *ISBA*  
!!
!!    PURPOSE
!!    -------
!       Monitor for the calculation of the surface fluxes and of the
!     prognostic variables of the surface over natural areas
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
!!	S. Belair           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/03/95 
!!      (J.Stein)   25/10/95  add the rain flux computation at the ground
!!                            and the lbc
!!      (J.Stein)   15/11/95  include the strong slopes cases
!!      (J.Stein)   06/02/96  bug correction for the precipitation flux writing 
!!      (J.Stein)   20/05/96  set the right IGRID value for the rain rate
!!      (J.Viviand) 04/02/97  add cold and convective precipitation rate
!!      (J.Stein)   22/06/97  use the absolute pressure    
!!      (V.Masson)  09/07/97  add directional z0 computations and RESA correction     
!!      (V.Masson)  13/02/98  simplify the routine: only vegetation computation
!!                            are now made here.
!!      (A.Boone)   05/10/98  add: Boone et al. (1999) 3 soil-water Layers version
!!      (V.Masson)                 Dumenil and Todini (1992) runoff
!!                                 Calvet (1998) biomass and CO2 assimilation
!!                                 Calvet (1998) LAI evolution
!!      (A.Boone)  03/15/99   Soil ice scheme: modify CG, C1, C2, WSAT, WFC, WILT,
!!                            LEG (add soil ice sublimation); Can modify TS and T2.
!!                            New variables WGI1, WGI2
!!      (A.Boone)  18/01/00   ISBA-ES (3-layer explicit snow scheme option)
!!                            (Boone and Etchevers 2000)
!!                            New variable PSNOWHEAT
!!      (V. Masson) 01/2004   wet leaves fraction computed in separate routine
!!                            all vegetation stress (ISBA, AGS, AST) routines
!!                            called at the same point
!!      (P. LeMoigne) 03/2004 computation of QSAT 
!!      (P. LeMoigne) 10/2004 halstead coefficient as diagnostic for isba
!!      (A. Bogatchev)09/2005 EBA snow option
!!      (P. LeMoigne) 02/2006 z0h and snow
!!      (B. Decharme) 05/2008 Add floodplains scheme
!!      (R. Hamdi)    01/09   Cp and L are not constants (As in ALADIN)
!!      (A.L. Gibelin)  03/2009 : Add respiration diagnostics
!!      A.L. Gibelin   06/09 : move calculations of CO2 fluxes
!!      A.L. Gibelin 07/2009 : Suppress PPST and PPSTF as outputs
!!      (A. Boone)    11/2009 Add local variable: total soil temperature change (before
!!                            phase change) for use by LWT scheme in ISBA-DIF. 
!!      (A. Boone)    03/2010 Add local variable: delta functions for LEG and LEGI
!!                            to numerically correct for when they should be
!!                            zero when hug(i) Qsat < Qa and Qsat > Qa
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
USE MODD_CO2V_PAR,   ONLY : XMC, XMCO2, XPCCO2
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_CSTS,       ONLY : XLVTT, XLSTT
!
USE MODD_TYPE_DATE_SURF, ONLY: DATE_TIME
!
USE MODI_SOIL
USE MODI_SOILDIF
USE MODI_SOILSTRESS
USE MODI_WET_LEAVES_FRAC
USE MODI_VEG
USE MODI_DRAG
USE MODI_SNOW3L_ISBA
USE MODI_E_BUDGET
USE MODI_HYDRO
USE MODI_ISBA_SNOW_AGR
!
USE MODI_COTWORES
USE MODI_COTWORESTRESS
USE MODI_ISBA_FLUXES
!
USE MODE_THERMOS
!
!
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
CHARACTER(LEN=*),     INTENT(IN)  :: HISBA      ! type of ISBA version:
!                                               ! '2-L' (default)
!                                               ! '3-L'
!                                               ! 'DIF'
CHARACTER(LEN=*),     INTENT(IN)  :: HPHOTO     ! Kind of photosynthesis
!                                               ! 'NON'
!                                               ! 'AGS'
!                                               ! 'LAI'
!                                               ! 'AST'
!                                               ! 'LST'
CHARACTER(LEN=*),     INTENT(IN)  :: HRUNOFF    ! surface runoff formulation
!                                               ! 'WSAT'
!                                               ! 'DT92'
!                                               ! 'SGH ' Topmodel
CHARACTER(LEN=*),     INTENT(IN)  :: HC1DRY     ! C1 for dry soil formulation
!                                               ! 'DEF' Default: Giard and Bazile
!                                               ! 'GB93' Giordani 1993, Braud 1993
!                                               ! (discontinuous at WILT)
CHARACTER(LEN=*),     INTENT(IN)  :: HSCOND     ! Thermal conductivity
!                                               ! 'NP89' = NP89 implicit method
!                                               ! 'PL98' = Peters-Lidard et al. 1998 used
!                                               ! for explicit computation of CG
!
CHARACTER(LEN=*),     INTENT(IN)  :: HSNOW_ISBA ! 'DEF' = Default F-R snow scheme
!                                               !         (Douville et al. 1995)
!                                               ! '3-L' = 3-L snow scheme (option)
!                                               !         (Boone and Etchevers 2000)
CHARACTER(LEN=*),     INTENT(IN)  :: HSNOWRES   ! 'DEF' = Default: Louis (ISBA)
!                                               ! 'RIL' = CROCUS (Martin) method
!                                               !  ISBA-SNOW3L turbulant exchange option
CHARACTER(LEN=*),      INTENT(IN) :: HCPSURF    ! Specific heat
!                                               ! 'DRY' = dry Cp
!                                               ! 'HUM' = humid Cp fct of qs
!
CHARACTER(LEN=*),   INTENT(IN) :: HSOILFRZ      ! soil freezing-physics option
!                                               ! 'DEF'   Default (Boone et al. 2000; Giard and Bazile 2000)
!                                               ! 'LWT'   phase changes as above, but relation between unfrozen 
!                                                     water and temperature considered
!
CHARACTER(LEN=*),     INTENT(IN)  :: HDIFSFCOND ! NOTE: Only used when HISBA = DIF
!                                               ! MLCH' = include the insulating effect of leaf
!                                               !         litter/mulch on the surface thermal cond.
!                                               ! 'DEF' = no mulch effect
CHARACTER(LEN=*),     INTENT(IN)  :: HDIF       ! NOTE: Only used when HISBA = DIF
!                                               ! 'BC' = Brook and Corey
!                                               ! 'VG' = Van Genuchten
!
TYPE(DATE_TIME), INTENT(IN)       :: TPTIME     ! current date and time
LOGICAL, INTENT(IN)               :: OFLOOD     ! Activation of the flooding scheme
LOGICAL, INTENT(IN)               :: OTEMP_ARP  ! True  = time-varying force-restore soil temperature (as in ARPEGE)
                                                ! False = No time-varying force-restore soil temperature (Default)
LOGICAL, INTENT(IN)               :: OGLACIER   ! True = Over permanent snow and ice, 
!                                                     initialise WGI=WSAT,
!                                                     Hsnow>=10m and allow 0.8<SNOALB<0.85
                                                ! False = No specific treatment
REAL,                 INTENT(IN) :: PTSTEP      ! timestep of the integration
REAL, DIMENSION(:), INTENT(IN)   :: PZREF       ! normal distance of the first
!                                               ! atmospheric level to the
!                                               ! orography
REAL, DIMENSION(:), INTENT(IN)   :: PUREF       ! reference height of the wind
!                                               ! NOTE this is different from ZZREF
!                                               ! ONLY in stand-alone/forced mode,
!                                               ! NOT when coupled to a model (MesoNH)
REAL, DIMENSION(:), INTENT(IN)   ::  PDIRCOSZW  ! Director Cosinus along z
!                                               ! directions at surface w-point
REAL,                 INTENT(IN) :: PCGMAX      ! maximum soil heat capacity
!
!* atmospheric variables
!  ---------------------
!
!            suffix 'A' stands for atmospheric variable at first model level
!            suffix 'S' stands for atmospheric variable at ground level
!
REAL, DIMENSION(:), INTENT(IN)  :: PTA        ! Temperature
REAL, DIMENSION(:), INTENT(IN)  :: PQA        ! specific humidity
REAL, DIMENSION(:), INTENT(IN)  :: PEXNA      ! Exner function
REAL, DIMENSION(:), INTENT(IN)  :: PRHOA      ! air density
!
REAL, DIMENSION(:), INTENT(IN)  :: PPS        ! Pressure
REAL, DIMENSION(:), INTENT(IN)  :: PEXNS      ! Exner function
!
REAL, DIMENSION(:), INTENT(IN)  :: PRR        ! Rain rate (in kg/m2/s)
REAL, DIMENSION(:), INTENT(IN)  :: PSR        ! Snow rate (in kg/m2/s)
!
REAL, DIMENSION(:), INTENT(IN)  :: PZENITH    ! solar zenith angle
REAL, DIMENSION(:), INTENT(IN)  :: PSW_RAD    ! solar   incoming radiation
REAL, DIMENSION(:), INTENT(IN)  :: PLW_RAD    ! thermal incoming radiation
!
REAL, DIMENSION(:), INTENT(IN)  :: PVMOD      ! modulus of the wind
!                                             ! parallel to the orography
!
! implicit coupling coefficients:
!
REAL, DIMENSION(:), INTENT(IN)  :: PPEW_A_COEF, PPEW_B_COEF,           &
                             PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF,      &
                             PPEQ_B_COEF  
!                                  PPEW_A_COEF ! A-wind coefficient
!                                  PPEW_B_COEF ! B-wind coefficient
!                                  PPET_A_COEF ! A-air temperature coefficient
!                                  PPET_B_COEF ! B-air temperature coefficient
!                                  PPEQ_A_COEF ! A-air specific humidity coefficient
!                                  PPEQ_B_COEF ! B-air specific humidity coefficient
!
!* vegetation parameters
!  ---------------------
!
REAL, DIMENSION(:), INTENT(IN)  :: PRSMIN     ! minimum stomatal resistance
REAL, DIMENSION(:), INTENT(IN)  :: PRGL       ! maximum solar radiation
!                                             ! usable in photosynthesis
REAL, DIMENSION(:), INTENT(IN)  :: PGAMMA     ! coefficient for the calculation
!                                             ! of the surface stomatal
!                                             ! resistance
REAL, DIMENSION(:), INTENT(IN)  :: PCV        ! 2*sqrt(pi/day)/sqrt(Cveg*hveg)
!                                             ! where Cveg and hveg are the
!                                             ! heat capacity and conductivity
!                                             ! of the vegetation
REAL, DIMENSION(:), INTENT(IN)  :: PRUNOFFD   ! depth over which sub-grid runoff computed (m)
REAL, DIMENSION(:), INTENT(IN)  :: PALB       ! albedo of vegetation
REAL, DIMENSION(:), INTENT(IN)  :: PWRMAX_CF  ! coefficient for maximum water interception 
!                                             ! storage capacity on the vegetation (-)
REAL, DIMENSION(:), INTENT(IN)    :: PVEG     ! fraction of vegetation of the
!                                             ! mesh covered by natural or
!                                             ! agricultural areas
!                                             ! 1-PVEG --> bare soil
REAL, DIMENSION(:), INTENT(IN)    :: PLAI     ! LAI as a function of time:
!                                             ! as a function of growth,
!                                             ! decay, assimilation.
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS    ! emissivity of natural surfaces
!                                             ! (without prognostic snow)
REAL, DIMENSION(:), INTENT(IN)    :: PZ0_WITH_SNOW ! roughness length for momentum
!                                             ! (with snow taken into account)
REAL, DIMENSION(:), INTENT(IN)    :: PZ0H_WITH_SNOW ! roughness length for heat
!                                             ! (with snow taken into account)
!
!* ISBA-Ags (with LAI evolution) parameters
!  ----------------------------------------
!
REAL, DIMENSION(:,:), INTENT(IN) :: PVEGTYPE ! fraction of each vegetation
!
!* ISBA-Ags parameters
!  -------------------
!
LOGICAL, DIMENSION(:), INTENT(IN) :: OSTRESSDEF ! vegetation response type to water
!                                         ! stress (true:defensive false:offensive)
REAL, DIMENSION(:), INTENT(IN) :: PGC     ! cuticular conductance (m s-1)
REAL, DIMENSION(:), INTENT(IN) :: PF2I    ! critical normilized soil water 
!                                         ! content for stress parameterisation
REAL, DIMENSION(:), INTENT(IN) :: PDMAX   ! maximum air saturation deficit
!                                         ! tolerate by vegetation
REAL, DIMENSION(:), INTENT(IN) :: PAH,PBH ! coefficients for herbaceous water stress 
!                                         ! response (offensive or defensive) 
!
REAL, DIMENSION(:), INTENT(IN) :: PCSP    ! atmospheric CO2 concentration
REAL, DIMENSION(:), INTENT(IN) :: PGMES   ! mesophyll conductance (m s-1)
!
REAL, DIMENSION(:), INTENT(IN) :: PABC    ! abscissa needed for integration
!                                         ! of net assimilation and stomatal
!                                         ! conductance over canopy depth
REAL, DIMENSION(:), INTENT(IN) :: PPOI    ! Gaussian weights (as above)
!
REAL, DIMENSION(:), INTENT(IN) :: PFZERO  ! ideal value of F, no photo- 
!                                         ! respiration or saturation deficit
REAL, DIMENSION(:), INTENT(IN) :: PEPSO   ! maximum initial quantum use
!                                         ! efficiency (mg J-1 PAR)
REAL, DIMENSION(:), INTENT(IN) :: PGAMM   ! CO2 conpensation concentration (ppmv)
REAL, DIMENSION(:), INTENT(IN) :: PQDGAMM ! Q10 function for CO2 conpensation 
!                                         ! concentration
REAL, DIMENSION(:), INTENT(IN) :: PQDGMES ! Q10 function for mesophyll conductance 
REAL, DIMENSION(:), INTENT(IN) :: PT1GMES ! reference temperature for computing 
!                                         ! compensation concentration function for 
!                                         ! mesophyll conductance: minimum
!                                         ! temperature 
REAL, DIMENSION(:), INTENT(IN) :: PT2GMES ! reference temperature for computing 
!                                         ! compensation concentration function for 
!                                         ! mesophyll conductance: maximum
!                                         ! temperature
REAL, DIMENSION(:), INTENT(IN) :: PAMAX   ! leaf photosynthetic capacity (kgCO2 m-2 s-1)
REAL, DIMENSION(:), INTENT(IN) :: PQDAMAX ! Q10 function for leaf photosynthetic capacity
REAL, DIMENSION(:), INTENT(IN) :: PT1AMAX ! reference temperature for computing 
!                                         ! compensation concentration function for leaf 
!                                         ! photosynthetic capacity: minimum
!                                         ! temperature
REAL, DIMENSION(:), INTENT(IN) :: PT2AMAX ! reference temperature for computing 
!                                         ! compensation concentration function for leaf 
!                                         ! photosynthetic capacity: maximum
!                                         ! temperature
!
!
!* subgrid-scale orography parameters
!  ----------------------------------
!
REAL, DIMENSION(:), INTENT(IN)  :: PZ0EFF     ! roughness length for momentum
!
REAL, DIMENSION(:), INTENT(IN)  :: PRUNOFFB   ! slope of the runoff curve
!
!* soil parameters
!  ---------------
!
REAL, DIMENSION(:), INTENT(IN)  :: PCGSAT     ! thermal coefficient at
!                                             ! saturation
REAL, DIMENSION(:), INTENT(IN)  :: PC1SAT     ! C1 coefficient at saturation
REAL, DIMENSION(:), INTENT(IN)  :: PC2REF     ! reference value of C2
REAL, DIMENSION(:,:), INTENT(IN):: PC3        ! C3 coefficient
REAL, DIMENSION(:), INTENT(IN)  :: PC4B       ! fiiting soil paramater for vertical diffusion (C4)
REAL, DIMENSION(:), INTENT(IN)  :: PC4REF     !         "
REAL, DIMENSION(:), INTENT(IN)  :: PACOEF     ! a and p coefficients for
REAL, DIMENSION(:), INTENT(IN)  :: PPCOEF     ! the wgeq calculations.
!
REAL, DIMENSION(:), INTENT(IN)  :: PTAUICE    ! characteristic time scale for phase change
!                                             ! within the soil
!
REAL, DIMENSION(:), INTENT(IN)  :: PWDRAIN    ! minimum Wg for drainage (m3/m3)
!
!
REAL, DIMENSION(:), INTENT(IN)  :: PTDEEP     ! Deep soil temperature (prescribed)
!                                             ! which models heating/cooling from
!                                             ! below the diurnal wave penetration
!                                             ! (surface temperature) depth.
REAL, DIMENSION(:), INTENT(IN)  :: PGAMMAT    ! Deep soil heat transfer coefficient:
!                                             ! assuming homogeneous soil so that
!                                             ! this can be prescribed in units of 
!                                             ! (1/days): associated time scale with
!                                             ! PTDEEP.
!
REAL, DIMENSION(:), INTENT(IN)  :: PPSN       ! fraction of the grid covered
!                                             ! by snow
REAL, DIMENSION(:), INTENT(IN)  :: PPSNG      ! fraction of the the bare
!                                             ! ground covered by snow
REAL, DIMENSION(:), INTENT(IN)  :: PPSNV      ! fraction of the the veg.
!                                             ! covered by snow
REAL, DIMENSION(:), INTENT(IN)  :: PPSNV_A    ! snow free albedo of vegetation 
                                      ! for EBA
REAL, DIMENSION(:), INTENT(IN)  :: PSNOWFREE_ALB_VEG  ! snow free albedo of vegetation 
REAL, DIMENSION(:), INTENT(IN)  :: PSNOWFREE_ALB_SOIL ! snow free albedo of soil 
!
REAL   ,DIMENSION(:),INTENT(IN)    :: PIRRIG
REAL   ,DIMENSION(:),INTENT(IN)    :: PWATSUP
REAL   ,DIMENSION(:),INTENT(IN)    :: PTHRESHOLD
LOGICAL,DIMENSION(:),INTENT(INOUT) :: LIRRIDAY
LOGICAL,DIMENSION(:),INTENT(IN)    :: LIRRIGATE
!
!
!* ISBA-SNOW3L variables/parameters:
!  ---------------------------------
!
! Prognostic variables:
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PSNOWALB   ! Snow albedo
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWSWE   ! Snow model layer liquid water equivalent or SWE (kg m-2)  
!                                                 ! NOTE for 'DEF' snow option, only uppermost element
!                                                 ! of this array is non-zero (as it's a one layer scheme)
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWHEAT  ! Snow layer heat content (J/m3) 
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWRHO   ! Snow layer average density (kg/m3)
!                                                 ! NOTE for 'DEF' snow option, only uppermost element
!                                                 ! of this array is used (as it's a one layer scheme)
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWGRAN1 ! Snow grain parameter 1 
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWGRAN2 ! Snow grain parameter 2 
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWHIST  ! Snow grain historical parameter
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWAGE  ! Snow grain age
!                                                 ! NOTE : methamorphism is only activated if the flag
!                                                 OSNOW_METAMO=TRUE
! 
! Diagnostics:
!
!                                                  (ISBA) surface
REAL, DIMENSION(:), INTENT(OUT)     :: PGRNDFLUX  ! snow/soil-biomass interface flux (W/m2)
!
REAL, DIMENSION(:), INTENT(OUT)     :: PHPSNOW    ! heat release from rainfall (W/m2)
REAL, DIMENSION(:), INTENT(OUT)     :: PSNOWHMASS ! snow heat content change from mass changes (J/m2)
REAL, DIMENSION(:), INTENT(OUT)     :: PSMELTFLUX ! energy removed from soil/vegetation surface
!                                                 ! when last traces of snow melted (W/m2)
REAL, DIMENSION(:), INTENT(OUT)     :: PRNSNOW    ! net radiative flux from snow (W/m2)
REAL, DIMENSION(:), INTENT(OUT)     :: PHSNOW     ! sensible heat flux from snow (W/m2)
REAL, DIMENSION(:), INTENT(OUT)     :: PGFLUXSNOW ! net heat flux from snow (W/m2)
REAL, DIMENSION(:), INTENT(OUT)     :: PUSTARSNOW ! friction velocity
REAL, DIMENSION(:), INTENT(OUT)     :: PSRSFC     ! Snow rate falling outside of snow
!                                                   covered grid area [kg/(m2 s)]
REAL, DIMENSION(:), INTENT(OUT)     :: PRRSFC     ! Rain rate falling outside of snow and flood
!                                                   covered grid area [kg/(m2 s)]
REAL, DIMENSION(:), INTENT(OUT)     :: PLESL      ! Evaporation (liquid) from wet snow (W/m2)
REAL, DIMENSION(:), INTENT(OUT)     :: PEMISNOW   ! snow surface emissivity
REAL, DIMENSION(:), INTENT(OUT)     :: PCDSNOW    ! drag coefficient for momentum over snow
REAL, DIMENSION(:), INTENT(OUT)     :: PCHSNOW    ! drag coefficient for heat over snow
REAL, DIMENSION(:), INTENT(OUT)     :: PTS_RAD    ! effective radiative temperature 
!                                                   of the natural surface (K)
REAL, DIMENSION(:), INTENT(OUT)     :: PTS        ! effective surface temperature (K)
REAL, DIMENSION(:), INTENT(OUT)     :: PHV        ! Halstead coefficient
REAL, DIMENSION(:), INTENT(OUT)     :: PQS        ! surface humidity (kg/kg)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWTEMP  ! snow layer temperatures (K)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWLIQ   ! snow layer liquid water content (m)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWDZ    ! snow layer thickness (m)
!
!
!* ISBA-DF variables/parameters:                  
!  ---------------------------------
! Parameters:
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PD_G       ! Depth of Bottom of Soil layers       (m)
REAL, DIMENSION(:,:), INTENT(IN)    :: PROOTFRAC  ! root fraction                        (-)
REAL, DIMENSION(:,:), INTENT(IN)    :: PWFC       ! field capacity profile               (m3/m3)
REAL, DIMENSION(:,:), INTENT(IN)    :: PWWILT     ! wilting point profile                (m3/m3)
REAL, DIMENSION(:,:), INTENT(IN)    :: PWSAT      ! porosity profile                     (m3/m3)
REAL, DIMENSION(:,:), INTENT(IN)    :: PBCOEF     ! soil water CH78 b-parameter          (-)
REAL, DIMENSION(:,:), INTENT(IN)    :: PCONDSAT   ! hydraulic conductivity at saturation (m/s)
REAL, DIMENSION(:,:), INTENT(IN)    :: PMPOTSAT   ! matric potential at saturation       (m)
REAL, DIMENSION(:,:), INTENT(IN)    :: PHCAPSOIL  ! soil heat capacity                   [J/(K m3)]
REAL, DIMENSION(:,:), INTENT(IN)    :: PCONDDRY   ! soil dry thermal conductivity        [W/(m K)]
REAL, DIMENSION(:,:), INTENT(IN)    :: PCONDSLD   ! soil solids thermal conductivity     [W/(m K)]
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PCONDSAT_EXP! exponential hydraulic conductivity
REAL, DIMENSION(:,:), INTENT(IN)    :: PEXP_DIF    ! exponential coef for hydraulic conductivity
REAL, DIMENSION(:,:), INTENT(IN)    :: PALPHA      ! Van Genuchten parameter
REAL, DIMENSION(:,:), INTENT(IN)    :: PN          ! Van Genuchten parameter
REAL, DIMENSION(:,:), INTENT(IN)    :: PM          ! Van Genuchten parameter
REAL, DIMENSION(:,:), INTENT(IN)    :: PL          ! Van Genuchten parameter
REAL, DIMENSION(:,:), INTENT(IN)    :: PWRES       ! Van Genuchten parameter
!
! Prognostic variables:
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PTG, PWG, PWGI
!                                      PTG   ! soil layer average temperatures        (K)
!                                      PWG   ! soil liquid volumetric water content   (m3/m3)
!                                      PWGI  ! soil frozen volumetric water content   (m3/m3)
!
REAL, DIMENSION(:), INTENT(INOUT)  :: PCPS, PLVTT, PLSTT
!
!* output soil parameters
!  ----------------------
!
REAL, DIMENSION(:), INTENT(OUT) :: PCG        ! heat capacity of the ground
REAL, DIMENSION(:), INTENT(OUT) :: PC1        ! coefficients for the moisure
REAL, DIMENSION(:), INTENT(OUT) :: PC2        ! equation.
REAL, DIMENSION(:), INTENT(OUT) :: PWGEQ      ! equilibrium volumetric water
!                                             ! content
REAL, DIMENSION(:), INTENT(OUT) :: PCT        ! area-averaged heat capacity
!
!* prognostic variables
!  --------------------
!
REAL, DIMENSION(:), INTENT(INOUT) :: PWR      ! liquid water retained on the
!                                             ! foliage of the vegetation
!                                             ! canopy
REAL, DIMENSION(:), INTENT(INOUT) :: PRESA    ! aerodynamic resistance
!
REAL, DIMENSION(:), INTENT(INOUT) :: PANFM    ! maximum leaf assimilation
!
!* diagnostic variables for Carbon assimilation
!  --------------------------------------------
!
REAL, DIMENSION(:), INTENT(INOUT) :: PAN      ! net CO2 assimilation
REAL, DIMENSION(:), INTENT(INOUT) :: PANDAY   ! daily net CO2 assimilation
REAL, DIMENSION(:), INTENT(OUT)   :: PANF     ! total assimilation over canopy
REAL, DIMENSION(:,:), INTENT(OUT) :: PIACAN   ! PAR in the canopy at different gauss level
REAL, DIMENSION(:,:),INTENT(OUT)  :: PRESP_BIOMASS_INST  ! instantaneous biomass respiration (kgCO2/kgair m/s)
REAL, DIMENSION(:,:), INTENT(INOUT) :: PRESP_BIOMASS     ! daily cumulated biomass respiration (kg/m2/day)
REAL, DIMENSION(:),INTENT(OUT)    :: PGPP     ! Gross Primary Production
!
!* diagnostic variables
!  --------------------
!
REAL, DIMENSION(:), INTENT(OUT) :: PCH        ! drag coefficient for heat
REAL, DIMENSION(:), INTENT(OUT) :: PCD        ! drag coefficient for momentum
REAL, DIMENSION(:), INTENT(OUT) :: PCDN       ! neutral drag coefficient for momentum
REAL, DIMENSION(:), INTENT(OUT) :: PRI        ! Richardson number
REAL, DIMENSION(:), INTENT(OUT) :: PHU        ! grid-area humidity of the soil
REAL, DIMENSION(:), INTENT(OUT) :: PHUG       ! ground relative humidity
REAL, DIMENSION(:), INTENT(OUT) :: PEMIST     ! total surface emissivity
REAL, DIMENSION(:), INTENT(OUT) :: PALBT      ! total surface albedo
REAL, DIMENSION(:), INTENT(OUT) :: PRS        ! surface stomatal resistance
!
!* surface fluxes
!  --------------
!
REAL, DIMENSION(:), INTENT(OUT) :: PRN        ! net radiation
REAL, DIMENSION(:), INTENT(OUT) :: PH         ! sensible heat flux
REAL, DIMENSION(:), INTENT(INOUT) :: PLE      ! total latent heat flux
REAL, DIMENSION(:), INTENT(OUT) :: PLEI       ! sublimation latent heat flux
REAL, DIMENSION(:), INTENT(OUT) :: PLEGI      ! latent heat of sublimation over frozen soil
REAL, DIMENSION(:), INTENT(OUT) :: PLEG       ! latent heat of evaporation
!                                             ! over the ground
REAL, DIMENSION(:), INTENT(OUT) :: PLEV       ! latent heat of evaporation
!                                             ! over the vegetation
REAL, DIMENSION(:), INTENT(OUT) :: PLES       ! latent heat of sublimation
!                                             ! over the snow
REAL, DIMENSION(:), INTENT(OUT) :: PLER       ! latent heat of the fraction
!                                             ! delta of water retained on the
!                                             ! foliage of the vegetation
REAL, DIMENSION(:), INTENT(OUT) :: PLETR      ! evapotranspiration of the rest
!                                             ! of the vegetation
REAL, DIMENSION(:), INTENT(OUT) :: PEVAP      ! total evaporative flux (kg/m2/s)
REAL, DIMENSION(:), INTENT(OUT) :: PGFLUX     ! flux through the ground
REAL, DIMENSION(:), INTENT(OUT) :: PRESTORE   ! surface restore flux (W m-2)
REAL, DIMENSION(:), INTENT(OUT) :: PUSTAR     ! friction velocity
REAL, DIMENSION(:), INTENT(OUT) :: PDRAIN     ! drainage
REAL, DIMENSION(:), INTENT(OUT) :: PRUNOFF    ! runoff
REAL, DIMENSION(:), INTENT(OUT) :: PMELT      ! melting rate of the snow (kg/m2/s)
REAL, DIMENSION(:), INTENT(OUT) :: PMELTADV   ! advection heat flux from snowmelt (W/m2)
!
! The following surface fluxes are from snow-free portion of grid
! box when the ISBA-ES option is ON. Otherwise, they are equal
! to the same variables without the _ISBA extension.
!
REAL, DIMENSION(:), INTENT(OUT) :: PRN_ISBA   ! net radiation
REAL, DIMENSION(:), INTENT(OUT) :: PH_ISBA    ! sensible heat flux
REAL, DIMENSION(:), INTENT(OUT) :: PLEG_ISBA  ! latent heat of evaporation (ground)
REAL, DIMENSION(:), INTENT(OUT) :: PLEGI_ISBA ! latent heat of sublimation (ground)
REAL, DIMENSION(:), INTENT(OUT) :: PLEV_ISBA  ! latent heat of evaporation (vegetation)
REAL, DIMENSION(:), INTENT(OUT) :: PLETR_ISBA ! latent heat of evaporation (transpiration)
REAL, DIMENSION(:), INTENT(OUT) :: PUSTAR_ISBA! friction velocity
REAL, DIMENSION(:), INTENT(OUT) :: PLER_ISBA  ! latent heat of evaporation (plant interception)
REAL, DIMENSION(:), INTENT(OUT) :: PLE_ISBA   ! total latent heat flux 
REAL, DIMENSION(:), INTENT(OUT) :: PLEI_ISBA  ! sublimation latent heat flux 
REAL, DIMENSION(:), INTENT(OUT) :: PGFLUX_ISBA! flux through the ground
!
CHARACTER(LEN=*),   INTENT(IN)  :: HKSAT      ! soil hydraulic profil option
!                                             ! 'DEF'  = ISBA homogenous soil
!                                             ! 'SGH'  = ksat exponential decay
!
CHARACTER(LEN=*),   INTENT(IN)  :: HRAIN      ! Rainfall spatial distribution
                                      ! 'DEF' = No rainfall spatial distribution
                                      ! 'SGH' = Rainfall exponential spatial distribution
                                      ! 
!
CHARACTER(LEN=*),   INTENT(IN)  :: HHORT      ! Horton runoff
                                      ! 'DEF' = no Horton runoff
                                      ! 'SGH' = Horton runoff
!                                        
REAL, DIMENSION(:), INTENT(IN)  :: PD_ICE     !depth of the soil column for the calculation
!                                              of the frozen soil fraction (m)
REAL, DIMENSION(:), INTENT(IN)  :: PKSAT_ICE  !hydraulic conductivity at saturation (m/s)
!                                            
REAL, DIMENSION(:), INTENT(IN)  :: PMUF       !fraction of the grid cell reached by the rainfall
REAL, DIMENSION(:), INTENT(IN)  :: PFSAT      !Topmodel saturated fraction
!
REAL, DIMENSION(:), INTENT(OUT) :: PHORT      !Horton runoff (kg/m2/s)
!
REAL, DIMENSION(:), INTENT(OUT) :: PDRIP      !Dripping from the vegetation (kg/m2/s)
REAL, DIMENSION(:), INTENT(OUT) :: PRRVEG     !Precip. intercepted by vegetation (kg/m2/s)
!
REAL, DIMENSION(:), INTENT(IN)   :: PFFG      !Floodplain fraction over the ground
REAL, DIMENSION(:), INTENT(IN)   :: PFFG_NOSNOW ! Without snow (ES)
REAL, DIMENSION(:), INTENT(IN)   :: PFFV      !Floodplain fraction over vegetation
REAL, DIMENSION(:), INTENT(IN)   :: PFFV_NOSNOW ! Without snow (ES)
REAL, DIMENSION(:), INTENT(IN)   :: PFFROZEN  !Fraction of frozen flood
REAL, DIMENSION(:), INTENT(IN)   :: PFF       !Floodplain fraction at the surface
REAL, DIMENSION(:), INTENT(IN)   :: PFALB     !Floodplain albedo
REAL, DIMENSION(:), INTENT(IN)   :: PFEMIS    !Floodplain emis
REAL, DIMENSION(:), INTENT(IN)   :: PFFLOOD   !Efective floodplain fraction
REAL, DIMENSION(:), INTENT(IN)   :: PPIFLOOD  !Floodplains potential infiltration           [kg/m²/s]
REAL, DIMENSION(:), INTENT(INOUT):: PIFLOOD   !Floodplains infiltration                     [kg/m²/s]
REAL, DIMENSION(:), INTENT(INOUT):: PPFLOOD   !Floodplains direct precipitation             [kg/m²/s]
REAL, DIMENSION(:), INTENT(INOUT):: PLE_FLOOD, PLEI_FLOOD !Floodplains latent heat flux     [W/m²]
!
REAL, DIMENSION(:),  INTENT(OUT) :: PAC_AGG  ! aggregated aerodynamic conductance
                                     ! for evaporative flux calculations
REAL, DIMENSION(:),  INTENT(OUT) :: PHU_AGG  ! aggregated relative humidity
                                     ! for evaporative flux calculations
!
REAL, DIMENSION(:), INTENT(IN)   ::  PSODELX  ! Pulsation for each layer (Only used if LTEMP_ARP=True)
REAL, DIMENSION(:), INTENT(IN)   ::  PLAT
REAL, DIMENSION(:), INTENT(IN)   ::  PLON
!
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PWR)) :: ZVMOD     ! modulus of the wind
REAL, DIMENSION(SIZE(PWR)) :: ZCS       ! heat capacity of the snow
REAL, DIMENSION(SIZE(PWR)) :: ZFROZEN1  ! ice fraction in superficial soil
REAL, DIMENSION(SIZE(PWR)) :: ZDELTA    ! fraction of the foliage
!                                       ! covered with intercepted
!                                       ! water
REAL, DIMENSION(SIZE(PWR)) :: ZQSAT     ! expression for the saturation 
!                                       ! specific humidity 
REAL, DIMENSION(SIZE(PWR)) :: ZDQSAT    ! expression for the saturation
!                                       ! specific humidity derivative
REAL, DIMENSION(SIZE(PWR)) :: ZTS_RAD   ! effective radiative temperature 
!                                         of the natural surface (K)
!                                         (snow free part in case of 3-L snow scheme)
!
REAL, DIMENSION(SIZE(PWR)) :: ZWRMAX    ! maximum canopy water interception
!
REAL, DIMENSION(SIZE(PWR)) :: ZF2       ! water stress coefficient
!
REAL, DIMENSION(SIZE(PWR)) :: ZF5       ! water stress coefficient (based on F2)
!                                       ! to enforce Etv=>0 as F2=>0
!
REAL, DIMENSION(SIZE(PWR)) :: ZDWGI1, ZDWGI2 ! Liquid equivalent volumetric soil
!                                              ice content time tendencies (m3/m3)
!
REAL, DIMENSION(SIZE(PWR)) :: ZHUGI    ! humidity over frozen bare ground
!
REAL, DIMENSION(SIZE(PWR)) :: ZEVAPCOR ! evaporation correction as last traces of snow
!                                      ! cover ablate
REAL, DIMENSION(SIZE(PWR)) :: ZLES3L   ! sublimation from ISBA-ES(3L)
REAL, DIMENSION(SIZE(PWR)) :: ZLEL3L   ! evaporation heat flux of water in the snow (W/m2)
REAL, DIMENSION(SIZE(PWR)) :: ZEVAP3L  ! evaporation flux over snow from ISBA-ES (kg/m2/s)
REAL, DIMENSION(SIZE(PWR)) :: ZSNOW_THRUFAL ! rate that liquid water leaves snow pack: 
!                                      ! ISBA-ES [kg/(m2 s)]
REAL, DIMENSION(SIZE(PWR)) :: ZT2M     ! restore temperature before time integration (K)
REAL, DIMENSION(SIZE(PWR)) :: ZTSM     ! surface temperature before time integration (K)
!
REAL, DIMENSION(SIZE(PWR)) :: ZLEG_DELTA  ! soil evaporation delta fn
REAL, DIMENSION(SIZE(PWR)) :: ZLEGI_DELTA ! soil sublimation delta fn
!
REAL, DIMENSION(SIZE(PTG,1),SIZE(PTG,2)) :: ZDELTAT
!                                      ! change in temperature over the time
!                                      ! step before adjustment owing to phase 
!                                      ! changes (K)
! ISBA-DF:
!                                                              
REAL, DIMENSION(SIZE(PWG,1),SIZE(PWG,2)) :: ZSOILHCAPZ ! ISBA-DF Soil heat capacity 
!                                                      ! profile [J/(m3 K)]
REAL, DIMENSION(SIZE(PWG,1),SIZE(PWG,2)) :: ZSOILCONDZ ! ISBA-DF Soil conductivity  
!                                                      ! profile  [W/(m K)]
!
REAL, DIMENSION(SIZE(PWG,1),SIZE(PWG,2)) :: ZF2WGHT    ! water stress factor
!
REAL, DIMENSION(SIZE(PWR)) :: ZTA_IC, ZQA_IC, ZVMOD_IC ! TA, QA and Wind speed updated values
!                                                      ! if implicit coupling with atmosphere used.
REAL, DIMENSION(SIZE(PWR)) :: ZTDIURN ! Ice maximum penetration depth for restore (m)
!
! Necessary to close the energy budget between surfex and the atmosphere:
!
REAL, DIMENSION(SIZE(PWR))   :: ZEMIST
REAL, DIMENSION(SIZE(PWR))   :: ZALBT
!
LOGICAL, DIMENSION(SIZE(PTG,1))  :: GMASK         ! mask where evolution occurs
LOGICAL                          :: GTMPMASK
INTEGER                          :: I
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*      1.0    Preliminaries
!              -------------
!
IF (LHOOK) CALL DR_HOOK('ISBA',0,ZHOOK_HANDLE)
!

PC1(:)          = XUNDEF
PC2(:)          = XUNDEF
PWGEQ(:)        = XUNDEF
ZCS(:)          = XUNDEF           
!
ZTA_IC(:)       = XUNDEF
ZQA_IC(:)       = XUNDEF
ZVMOD_IC(:)     = XUNDEF
!
ZTDIURN     (:) = 0.0
!
ZEMIST      (:) = XUNDEF
ZALBT       (:) = XUNDEF
!
ZSOILHCAPZ(:,:) = XUNDEF
ZSOILCONDZ(:,:) = XUNDEF
ZF2WGHT   (:,:) = XUNDEF
ZF5         (:) = XUNDEF
!
! Save surface and sub-surface temperature values at beginning of time step for 
! budget and flux calculations:
!
ZTSM(:)         = PTG(:,1)
ZT2M(:)         = PTG(:,2)
!
PRS         (:)   = 0.0
!
!-------------------------------------------------------------------------------
!
!*      2.0    Soil parameters
!              ---------------
!
IF(HISBA =='2-L' .OR. HISBA == '3-L')THEN
!
   CALL SOIL (HC1DRY, HSCOND, HSNOW_ISBA,                                       &
     PSNOWRHO(:,1), PVEG,                                                       &
     PCGSAT, PCGMAX,                                                            &
     PC1SAT, PC2REF, PACOEF, PPCOEF, PCV,                                       &
     PPSN, PPSNG, PPSNV, PFFG, PFFV, PFF,                                       &
     PCG, PC1, PC2, PWGEQ, PCT, ZCS, ZFROZEN1,                                  &
     PTG(:,1), PWG, PWGI,                                                       &
     PHCAPSOIL(:,1), PCONDDRY(:,1), PCONDSLD(:,1),                              &
     PBCOEF(:,1), PWSAT(:,1), PWWILT(:,1),                                      &
     HKSAT,PCONDSAT,PFFG_NOSNOW,PFFV_NOSNOW                                     )  
!
ELSE
!
   CALL SOILDIF (HSCOND, HDIFSFCOND, HDIF,                                      &
     PVEG, PCV, PFFG_NOSNOW, PFFV_NOSNOW,                                       &
     PCG, PCGMAX, PCT, ZFROZEN1,                                                &
     PWG, PWGI,                                                                 &
     PHCAPSOIL, PCONDDRY, PCONDSLD,                                             &
     PBCOEF, PWSAT, PMPOTSAT, PALPHA, PN, PM,                                   &
     PWRES, ZSOILCONDZ, ZSOILHCAPZ                                              )  
ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      3.0    Explicit snow scheme
!              --------------------
!
CALL SNOW3L_ISBA(HISBA, HSNOW_ISBA, HSNOWRES, OGLACIER, TPTIME, PTSTEP,         &
           PVEGTYPE, PSNOWSWE, PSNOWHEAT, PSNOWRHO, PSNOWALB,                   &
           PSNOWGRAN1, PSNOWGRAN2, PSNOWHIST,PSNOWAGE,                          &
           PTG(:,1), PCG, PCT, ZSOILCONDZ(:,1),                                 &
           PPS, PTA, PSW_RAD, PQA, PVMOD, PLW_RAD, PRR, PSR,                    &
           PRHOA, PUREF, PEXNS, PEXNA, PDIRCOSZW,                               &
           PZREF, PZ0_WITH_SNOW, PZ0EFF, PZ0H_WITH_SNOW, PALB, PD_G(:,1),       &
           PPEW_A_COEF, PPEW_B_COEF,                                            &
           PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,                  &
           ZSNOW_THRUFAL, PGRNDFLUX, ZEVAPCOR,                                  &
           PRNSNOW, PHSNOW, PGFLUXSNOW, PHPSNOW, ZLES3L, ZLEL3L, ZEVAP3L,       &
           PUSTARSNOW, PPSN, PSRSFC, PRRSFC, PSMELTFLUX,                        &
           PEMISNOW, PCDSNOW, PCHSNOW, PSNOWTEMP, PSNOWLIQ, PSNOWDZ,            &
           PSNOWHMASS, PZENITH, PLAT, PLON                                      )  
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      4.0    Fraction of leaves occupied by intercepted water
!              ------------------------------------------------
!
CALL WET_LEAVES_FRAC(PWR, PVEG, PWRMAX_CF, PZ0_WITH_SNOW, PLAI, ZWRMAX, ZDELTA)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      5.0    Plant stress due to soil water deficit
!              --------------------------------------
!
CALL SOILSTRESS(HISBA, ZF2,                &
         PROOTFRAC, PWSAT, PWFC, PWWILT,   &
         PWG, PWGI, ZF2WGHT, ZF5           )  
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      6.0    Plant stress, stomatal resistance and, possibly, CO2 assimilation
!              --------------------------------------------------------------------
!
ZQSAT=QSAT(PTG(:,1),PPS(:))
!
IF (HPHOTO=='NON') THEN
   CALL VEG(PSW_RAD, PTA, PQA, PPS, PRGL, PLAI, PRSMIN, PGAMMA, ZF2, PRS)
!
ELSE IF (MAXVAL(PGMES).NE.XUNDEF .OR. MINVAL(PGMES).NE.XUNDEF) THEN
  IF (HPHOTO=='AGS' .OR. HPHOTO=='LAI') THEN
   CALL COTWORES(PTSTEP, PANF, PABC, PPOI, PAN, PANDAY,                          &
          PANFM, PCSP, PTG(:,1),                                                   &
         PSW_RAD, PRESA, PQA, ZQSAT, PLE, PPSNV, ZDELTA, ZF2,                      &
         PLAI, PRS, PZENITH, PGC, PRHOA,                                           &
         PFZERO, PEPSO, PGAMM, PQDGAMM, PGMES, PQDGMES,                            &
         PT1GMES, PT2GMES, PAMAX, PQDAMAX, PT1AMAX, PT2AMAX, PIACAN,               &
         PGPP, PRESP_BIOMASS_INST(:,1)                                             )  
!
  ELSE IF (HPHOTO=='AST' .OR. HPHOTO=='LST' .OR. HPHOTO=='NIT' .OR. HPHOTO=='NCB') THEN

   CALL COTWORESTRESS(PTSTEP, PVEGTYPE, OSTRESSDEF, PAH, PBH, PF2I,           &
         PANF, PABC, PPOI, PAN, PANDAY, PANFM, PCSP, PTG(:,1),                  &
         PSW_RAD, PRESA, PQA, ZQSAT, PLE, PPSNV, ZDELTA, ZF2,                   &
         PLAI, PRS, PZENITH, PGC, PRHOA, PDMAX,                                 &
         PFZERO, PEPSO, PGAMM, PQDGAMM, PGMES, PQDGMES,                         &
         PT1GMES, PT2GMES, PAMAX, PQDAMAX, PT1AMAX, PT2AMAX, PIACAN,            &
         PGPP, PRESP_BIOMASS_INST(:,1)                                          )  
  ENDIF
ELSE
  PRESP_BIOMASS_INST(:,1) = 0.0
  PGPP(:) = 0.0   
ENDIF
!
! Mask where vegetation evolution is performed (just before solar midnight)
GTMPMASK = ( TPTIME%TIME - PTSTEP < 0. ) .AND. ( TPTIME%TIME >= 0. )
GMASK(:) = GTMPMASK
       
IF (HPHOTO/='NON') THEN
  !        
  ! Cumulate daily leaf respiration in RESP_BIOMASS(:,1)
  ! compiler doesnt' like mask in array of different size.
  DO I=1,SIZE(GMASK,1)
    IF   (GMASK(I)) PRESP_BIOMASS(I,1) = 0.0
  END DO
  !
  PRESP_BIOMASS(:,1) = PRESP_BIOMASS(:,1) + PRESP_BIOMASS_INST(:,1) * (PTSTEP*PRHOA(:)*XMC)/(XPCCO2*XMCO2)
!
ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      7.0    Aerodynamic drag and heat transfer coefficients
!              -----------------------------------------------
!
CALL DRAG(HISBA, HSNOW_ISBA, HCPSURF,                                         &
    PTG(:,1), PWG(:,1), PWGI(:,1), PEXNS, PEXNA, PTA, PVMOD, PQA, PRR, PSR,     &
    PPS, PRS, PVEG, PZ0_WITH_SNOW, PZ0EFF, PZ0H_WITH_SNOW,                      &
    PWFC(:,1), PWSAT(:,1), PPSNG, PPSNV, PZREF, PUREF,                          &
    PDIRCOSZW, ZDELTA, ZF5, PRESA,                                              &
    PCH, PCD, PCDN, PRI, PHUG, ZHUGI, PHV, PHU, PCPS, PQS,                      &
    PFFG, PFFV, PFF, PFFG_NOSNOW, PFFV_NOSNOW,                                  &
    ZLEG_DELTA, ZLEGI_DELTA                                                     )  
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      8.0    Resolution of the surface and soil energy budget
!              ------------------------------------------------
!
CALL E_BUDGET(HISBA, HSNOW_ISBA, OFLOOD, OTEMP_ARP, PSODELX,                  &
        PUREF, PPEW_A_COEF, PPEW_B_COEF,                                        &
        PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,                     &
        PVMOD, PCD,                                                             &
        PTG, PTSTEP, PSNOWALB,                                                  &
        PSW_RAD, PLW_RAD, PTA, PQA, PPS, PRHOA,                                 &
        PEXNS,PEXNA, PCPS, PLVTT, PLSTT,                                        &
        PVEG, PHUG, ZHUGI, PHV,                                                 &
        ZLEG_DELTA, ZLEGI_DELTA,                                                &
        PEMIS, PALB, PRESA,                                                     &
        PCT, PPSN, PPSNV, PPSNG,                                                &
        PGRNDFLUX, PSMELTFLUX, ZSNOW_THRUFAL,                                   &
        PD_G, ZSOILCONDZ, ZSOILHCAPZ,                                           &
        ZALBT, ZEMIST, ZQSAT, ZDQSAT,                                           &
        ZFROZEN1, PTDEEP, PGAMMAT, ZTA_IC, ZQA_IC, ZVMOD_IC,                    &
        PSNOWFREE_ALB_VEG, PPSNV_A, PSNOWFREE_ALB_SOIL,                         &
        PFFG, PFFV, PFF, PFFROZEN, PFALB, PFEMIS, ZDELTAT)  
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      9.0    Energy and momentum fluxes
!              --------------------------
!
!*******************************************************************************
! WARNING: at this stage, ZALBT and ZEMIST have two different meanings according
!          to the ISBA snow-scheme option:
!  'D95' : they represent aggregated (snow + flood + snow-flood-free) albedo and emissivity
!  '3-L' : they represent                    flood + snow-flood-free  albedo and emissivity
!*******************************************************************************
!
CALL ISBA_FLUXES(HISBA, HSNOW_ISBA, HSOILFRZ, HDIF, OTEMP_ARP, PTSTEP, PSODELX, &
           PCD, PVMOD, PSW_RAD, PLW_RAD, ZTA_IC, ZQA_IC,                        &
           ZVMOD_IC, PRHOA, PEXNS, PEXNA, PCPS, PLVTT, PLSTT,                   &
           PLAI, PVEG, PHUG, ZHUGI, PHV, ZLEG_DELTA, ZLEGI_DELTA, ZDELTA, PRESA,&
           ZF5, PRS, ZCS, PCG, PCT, PSNOWSWE(:,1), ZT2M, ZTSM,                  &
           PPSN, PPSNV, PPSNG, ZFROZEN1, PTAUICE,                               &
           ZALBT, ZEMIST, ZQSAT, ZDQSAT, ZSNOW_THRUFAL,                         &
           PRN, PH, PLE, PLEG, PLEGI, PLEV,                                     &
           PLES, PLER, PLETR, PEVAP, PGFLUX, PMELTADV, PMELT, PRESTORE,         &
           PUSTAR, ZTS_RAD,                                                     &
           ZDWGI1, ZDWGI2, ZDELTAT,                                             &
           ZSOILCONDZ, ZSOILHCAPZ, PWSAT, PMPOTSAT,                             &
           PBCOEF, PD_G, PTG, PWGI, PWG, PSRSFC, PPSNV_A,                       &
           PFFG, PFFV, PFF, PFFROZEN, PLE_FLOOD, PLEI_FLOOD, PSNOWTEMP(:,1),    &
           PALPHA, PN, PM, PWRES, ZTDIURN                                       )  
!
! Compute aggregated coefficients for evaporation
! Sum(LEV+LEG+LEGI+LES) = ACagg * Lv * RHOA * (HUagg.Qsat - Qa)
!
PAC_AGG(:) =   1. / PRESA(:) / XLVTT                              &
     * ( XLVTT*    PVEG(:) *(1.-PPSNV(:))                 *PHV(:)   &
       + XLVTT*(1.-PVEG(:))*(1.-PPSNG(:))*(1.-ZFROZEN1(:))          &
       + XLSTT*(1.-PVEG(:))*(1.-PPSNG(:))*    ZFROZEN1(:)           &
       + XLSTT*                 PPSN (:)                            )  
!
PHU_AGG(:) =   1. / (PRESA(:) * PAC_AGG(:)) / XLVTT               &
     * ( XLVTT*    PVEG(:) *(1.-PPSNV(:))                 *PHV(:)   &
       + XLVTT*(1.-PVEG(:))*(1.-PPSNG(:))*(1.-ZFROZEN1(:))*PHUG(:)  &
       + XLSTT*(1.-PVEG(:))*(1.-PPSNG(:))*    ZFROZEN1(:) *ZHUGI(:) &
       + XLSTT*                 PPSN (:)                            )  
!
!*******************************************************************************
! WARNING: at this stage, all fluxes have two different meanings according
!          to the ISBA snow-scheme option:
!  'D95' : they represent aggregated (snow + flood + snow-flood-free) fluxes
!  '3-L' : they represent                    flood + snow-flood-free  fluxes
!
! The variables concerned by this are: PRN, PH, PLE, PLEI, PLEG, PLEGI, PLEV, PLES, 
!                                      PLER, PLETR, PEVAP, PUSTAR, PGFLUX
!*******************************************************************************
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*     10.0    Water transfers and phase change in the soil
!              --------------------------------------------
!
CALL HYDRO(HISBA, HDIF, HSNOW_ISBA, HRUNOFF, OGLACIER,                          &
     PTSTEP, PVEGTYPE,                                                          &
     PRRSFC, PSRSFC, PLEV, PLETR, PLEG, PLES,                                   &
     PRUNOFFB, PWDRAIN,                                                         &
     PC1, PC2, PC3, PC4B, PC4REF, PWGEQ, PCT,                                   &
     PVEG, ZWRMAX, PRUNOFFD, PMELT,                                             &
     ZDWGI1, ZDWGI2, PLEGI,                                                     &
     PPSNV, PPSNG, ZSNOW_THRUFAL, ZEVAPCOR,                                     &
     PWR,                                                                       &
     PSNOWSWE(:,1), PSNOWALB, PSNOWRHO(:,1),                                    &
     PBCOEF, PWSAT, PCONDSAT, PMPOTSAT, PWFC,                                   &
     PWWILT, ZF2WGHT, ZF2, PD_G, PPS,                                           &
     PWG, PWGI, PTG,                                                            &
     PDRAIN, PRUNOFF,                                                           &
     PIRRIG, PWATSUP, PTHRESHOLD, LIRRIDAY, LIRRIGATE,                          &
     HKSAT, HRAIN, HHORT, PMUF, PFSAT, PKSAT_ICE, PD_ICE, PHORT, PDRIP,         &
     PFFG, PFFV, PFFLOOD, PPIFLOOD, PIFLOOD, PPFLOOD, PRRVEG,                   &
   PCONDSAT_EXP, PEXP_DIF, PALPHA, PN, PM, PL, PWRES, ZTDIURN                   )  
!
!-------------------------------------------------------------------------------
!
!*     11.0    Aggregated output fluxes and diagnostics
!              -----------------------------------------
!
!* add snow component to output radiative parameters and fluxes in case 
!  of 3-L snow scheme
!

CALL ISBA_SNOW_AGR( HSNOW_ISBA,                                 &
          ZEMIST, ZALBT,                                          &
          PPSN, PPSNG, PPSNV,                                     &
          PRN, PH, PLE, PLEI, PLEG, PLEGI, PLEV, PLES, PLER,      &
          PLETR, PEVAP, PGFLUX, PLVTT, PLSTT,                     &
          PUSTAR,                                                 &
          ZLES3L, ZLEL3L, ZEVAP3L,                                &
          PSNOWALB,                                               &
          PRNSNOW, PHSNOW,  PHPSNOW,                              &
          PGFLUXSNOW, PUSTARSNOW,                                 &
          PGRNDFLUX, PLESL,                                       &
          PEMISNOW,                                               &
          PSNOWTEMP, PTS_RAD, PTS, PSNOWHMASS,                    &
          PRN_ISBA, PH_ISBA, PLEG_ISBA, PLEGI_ISBA, PLEV_ISBA,    &
          PLETR_ISBA, PUSTAR_ISBA, PLER_ISBA, PLE_ISBA,           &
          PLEI_ISBA, PGFLUX_ISBA, PMELTADV,                       &
          ZTS_RAD, PTG,                                           &
          PEMIST, PALBT, PLE_FLOOD, PLEI_FLOOD, PFFG, PFFV, PFF   )  
!
!***************************************************************************
! All output fluxes and radiative variables have recovered the same physical
! meaning, that is they are aggregated quantities (snow + snow-free)
!***************************************************************************
!
IF (LHOOK) CALL DR_HOOK('ISBA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ISBA
