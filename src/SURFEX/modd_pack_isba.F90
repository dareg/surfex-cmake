!######################
MODULE MODD_PACK_ISBA
!######################
!
!!****  *MODD_PACK_ISBA - declaration of packed surface parameters for ISBA scheme
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      A. Boone   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       20/09/02
!!      A.L. Gibelin    04/2009 : BIOMASS and RESP_BIOMASS arrays 
!!      A.L. Gibelin    04/2009 : TAU_WOOD for NCB option 
!!      A.L. Gibelin    05/2009 : Add carbon spinup
!!      A.L. Gibelin    06/2009 : Soil carbon variables for CNT option
!!      A.L. Gibelin    07/2009 : Suppress RDK and transform GPP as a diagnostic
!!      A.L. Gibelin    07/2009 : Suppress PPST and PPSTF as outputs
!!      P. Samuelsson   10/2014 : MEB and additional snow albedos
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TYPE_DATE_SURF      
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!-------------------------------------------------------------------------------
!
TYPE PACK_ISBA_t

  INTEGER :: NSIZE_LSIMPLE
  INTEGER :: NSIZE_L0
  INTEGER :: NSIZE_NSIMPLE
  INTEGER :: NSIZE_N0
  INTEGER :: NSIZE_TSIMPLE
  INTEGER :: NSIZE_T0
  INTEGER :: NSIZE_SIMPLE
  INTEGER :: NSIZE_GROUND
  INTEGER :: NSIZE_VEGTYPE
  INTEGER :: NSIZE_TG
  INTEGER :: NSIZE_SNOW
  INTEGER :: NSIZE_ALB
  INTEGER :: NSIZE_2
  INTEGER :: NSIZE_BIOMASS
  INTEGER :: NSIZE_SOILCARB
  INTEGER :: NSIZE_LITTLEVS
  INTEGER :: NSIZE_LITTER
  INTEGER :: NSIZE_0
  INTEGER :: NSIZE_00
  INTEGER :: NSIZE_000
  INTEGER :: NSIZE_01
  LOGICAL, POINTER, DIMENSION(:,:) :: LBLOCK_SIMPLE 
  LOGICAL, POINTER, DIMENSION(:,:) :: LBLOCK_0
  INTEGER, POINTER, DIMENSION(:,:) :: NBLOCK_SIMPLE 
  INTEGER, POINTER, DIMENSION(:,:) :: NBLOCK_0
  TYPE(DATE_TIME), POINTER, DIMENSION(:,:) :: TBLOCK_SIMPLE
  TYPE(DATE_TIME), POINTER, DIMENSION(:,:) :: TBLOCK_0
  REAL, POINTER, DIMENSION(:,:) :: XBLOCK_SIMPLE
  REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_GROUND
  REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_VEGTYPE
  REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_TG
  REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_SNOW
  REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_ALB
  REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_2
  REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_BIOMASS
  REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_SOILCARB
  REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_LITTLEVS
  REAL, POINTER, DIMENSION(:,:,:,:) :: XBLOCK_LITTER
  REAL, POINTER, DIMENSION(:,:) :: XBLOCK_0
  REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_00
  REAL, POINTER, DIMENSION(:,:,:,:) :: XBLOCK_000
  REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_01
!
! Mask and number of grid elements containing patches/tiles:
!
  REAL, POINTER, DIMENSION(:,:)  :: XP_VEGTYPE_PATCH ! fraction of each vegetation type for
!                                                      ! each vegetation unit/patch              (-)
!
! General surface parameters:
!
  REAL, POINTER, DIMENSION(:)    :: XP_SSO_SLOPE     ! orography slope within the grid mesh    (-)
  REAL, POINTER, DIMENSION(:)    :: XP_LAT           ! latitude    (-)
  REAL, POINTER, DIMENSION(:)    :: XP_LON           ! latitude    (-)
!
! Subgrid orography parameters
!
  REAL, DIMENSION(:), POINTER :: XP_AOSIP,XP_AOSIM,XP_AOSJP,XP_AOSJM
! directional A/S quantities in 4 coordinate directions
! (IP: i index up;  IM: i index down;  JP: j index up;  JM: j index down)
! They are used in soil routines to compute effective roughness length
!
  REAL, DIMENSION(:), POINTER :: XP_HO2IP,XP_HO2IM,XP_HO2JP,XP_HO2JM
! directional h/2 quantities in 4 coordinate directions
! (IP: i index up;  IM: i index down;  JP: j index up;  JM: j index down)
! They are used in soil routines to compute effective roughness length
!
  REAL, DIMENSION(:), POINTER :: XP_Z0EFFIP,XP_Z0EFFIM,XP_Z0EFFJP,XP_Z0EFFJM
! directional total roughness lenghts in 4 coordinate directions
! (IP: i index up;  IM: i index down;  JP: j index up;  JM: j index down)
!
  REAL, POINTER, DIMENSION(:) :: XP_Z0REL         ! relief roughness length                 (m)

!
! Input Parameters:
!
! - bare soil:
!
  REAL, POINTER, DIMENSION(:,:) :: XP_CLAY         ! clay fraction profile                   (-)
  REAL, POINTER, DIMENSION(:,:) :: XP_SAND         ! sand fraction profile                   (-)

  REAL, POINTER, DIMENSION(:) :: XP_ALBNIR_DRY     ! near-infra-red albedo of dry soil       (-)
  REAL, POINTER, DIMENSION(:) :: XP_ALBVIS_DRY     ! visible albedo of dry soil              (-)
  REAL, POINTER, DIMENSION(:) :: XP_ALBUV_DRY      ! UV albedo of dry soil                   (-)
  REAL, POINTER, DIMENSION(:) :: XP_ALBNIR_WET     ! near-infra-red albedo of wet soil       (-)
  REAL, POINTER, DIMENSION(:) :: XP_ALBVIS_WET     ! visible albedo of wet soil              (-)
  REAL, POINTER, DIMENSION(:) :: XP_ALBUV_WET      ! UV albedo of wet soil                   (-)
  REAL, POINTER, DIMENSION(:) :: XP_ALBNIR_SOIL    ! near-infra-red albedo of wet soil       (-)
  REAL, POINTER, DIMENSION(:) :: XP_ALBVIS_SOIL    ! visible albedo of soil                  (-)
  REAL, POINTER, DIMENSION(:) :: XP_ALBUV_SOIL     ! UV albedo of soil                       (-)
!
! - vegetation + bare soil:
!
  REAL, POINTER, DIMENSION(:) :: XP_Z0_O_Z0H       ! ratio of surface roughness lengths
!                                                    ! (momentum to heat)                      (-)
  REAL, POINTER, DIMENSION(:) :: XP_ALBNIR         ! near-infra-red albedo                   (-)
  REAL, POINTER, DIMENSION(:) :: XP_ALBVIS         ! visible albedo                          (-)
  REAL, POINTER, DIMENSION(:) :: XP_ALBUV          ! UV albedo                               (-)
  REAL, POINTER, DIMENSION(:) :: XP_EMIS           ! snow-free surface emissivity                      (-)
  REAL, POINTER, DIMENSION(:) :: XP_Z0             ! snow-free surface roughness length                (m)
!
! - vegetation :
!
  REAL, POINTER, DIMENSION(:) :: XP_ALBNIR_VEG     ! near-infra-red albedo of vegetation     (-)
  REAL, POINTER, DIMENSION(:) :: XP_ALBVIS_VEG     ! visible albedo of vegetation            (-)
  REAL, POINTER, DIMENSION(:) :: XP_ALBUV_VEG      ! UV albedo of vegetation                 (-)
!
!
! - vegetation: default option (Jarvis) and general parameters:
!
  REAL, POINTER, DIMENSION(:) :: XP_VEG            ! vegetation cover fraction               (-)
  REAL, POINTER, DIMENSION(:) :: XP_WRMAX_CF       ! coefficient for maximum water 
!                                                    ! interception 
!                                                    ! storage capacity on the vegetation      (-)
  REAL, POINTER, DIMENSION(:) :: XP_RSMIN          ! minimum stomatal resistance             (s/m)
  REAL, POINTER, DIMENSION(:) :: XP_GAMMA          ! coefficient for the calculation
!                                                    ! of the surface stomatal
!                                                    ! resistance
  REAL, POINTER, DIMENSION(:) :: XP_CV             ! vegetation thermal inertia coefficient  (K m2/J)
  REAL, POINTER, DIMENSION(:) :: XP_RGL            ! maximum solar radiation
!                                                    ! usable in photosynthesis                (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XP_ROOTFRAC     ! root fraction profile ('DIF' option)
!
! - vegetation: Ags parameters ('AGS', 'LAI', 'AST', 'LST', 'NIT', 'NCB' options)
!
  REAL, POINTER, DIMENSION(:)    :: XP_BSLAI      ! ratio d(biomass)/d(lai)                 (kg/m2)
  REAL, POINTER, DIMENSION(:)    :: XP_LAIMIN     ! minimum LAI (Leaf Area Index)           (m2/m2)
  REAL, POINTER, DIMENSION(:)    :: XP_SEFOLD     ! e-folding time for senescence           (s)
  REAL, POINTER, DIMENSION(:)    :: XP_H_TREE     ! height of trees                         (m)
  REAL, POINTER, DIMENSION(:)    :: XP_ANF        ! total assimilation over canopy          (
  REAL, POINTER, DIMENSION(:)    :: XP_ANMAX      ! maximum photosynthesis rate             (
  REAL, POINTER, DIMENSION(:)    :: XP_FZERO      ! ideal value of F, no photo- 
!                                                   ! respiration or saturation deficit       (
  REAL, POINTER, DIMENSION(:)    :: XP_EPSO       ! maximum initial quantum use             
!                                                   ! efficiency                              (mg J-1 PAR)
  REAL, POINTER, DIMENSION(:)    :: XP_GAMM       ! CO2 conpensation concentration          (ppmv)
  REAL, POINTER, DIMENSION(:)    :: XP_QDGAMM     ! Log of Q10 function for CO2 conpensation 
!                                                   ! concentration                           (-)
  REAL, POINTER, DIMENSION(:)    :: XP_GMES       ! mesophyll conductance                   (m s-1)
  REAL, POINTER, DIMENSION(:)    :: XP_RE25       ! Ecosystem respiration parameter         (kg m-2 s-1)
  REAL, POINTER, DIMENSION(:)    :: XP_QDGMES     ! Log of Q10 function for mesophyll conductance  (-)
  REAL, POINTER, DIMENSION(:)    :: XP_T1GMES     ! reference temperature for computing 
!                                                   ! compensation concentration function for 
!                                                   ! mesophyll conductance: minimum
!                                                   ! temperature                             (K)
  REAL, POINTER, DIMENSION(:)    :: XP_T2GMES     ! reference temperature for computing 
!                                                   ! compensation concentration function for 
!                                                   ! mesophyll conductance: maximum
!                                                   ! temperature                             (K)
  REAL, POINTER, DIMENSION(:)    :: XP_AMAX       ! leaf photosynthetic capacity            (kg m-2 s-1)
  REAL, POINTER, DIMENSION(:)    :: XP_QDAMAX     ! Log of Q10 function for leaf photosynthetic 
!                                                   ! capacity                                (-)
  REAL, POINTER, DIMENSION(:)    :: XP_T1AMAX     ! reference temperature for computing 
!                                                   ! compensation concentration function for 
!                                                   ! leaf photosynthetic capacity: minimum
!                                                   ! temperature                             (K)
  REAL, POINTER, DIMENSION(:)    :: XP_T2AMAX     ! reference temperature for computing 
!                                                   ! compensation concentration function for 
!                                                   ! leaf photosynthetic capacity: maximum
!                                                   ! temperature                             (K)
!
! - vegetation: Ags Stress parameters ('AST', 'LST', 'NIT', 'NCB' options)
!
  LOGICAL, POINTER, DIMENSION(:) :: LP_STRESS      ! vegetation response type to water
!                                                    ! stress (true:defensive false:offensive) (-)
  REAL, POINTER, DIMENSION(:)    :: XP_F2I         ! critical normilized soil water 
!                                                    ! content for stress parameterisation
  REAL, POINTER, DIMENSION(:)    :: XP_GC          ! cuticular conductance                   (m s-1)
  REAL, POINTER, DIMENSION(:)    :: XP_AH          ! coefficients for herbaceous water stress 
!                                                    ! response (offensive or defensive)       (log(mm/s))
  REAL, POINTER, DIMENSION(:)    :: XP_BH          ! coefficients for herbaceous water stress 
!                                                    ! response (offensive or defensive)       (-)
  REAL, POINTER, DIMENSION(:)    :: XP_TAU_WOOD    ! residence time in woody biomass         (s)
  REAL, POINTER, DIMENSION(:)    :: XP_DMAX        ! maximum air saturation deficit
!                                                    ! tolerate by vegetation                  (kg/kg)
!
! - vegetation: Ags Nitrogen-model parameters ('NIT', 'NCB' option)
!
  REAL, POINTER, DIMENSION(:)    :: XP_CE_NITRO      ! leaf aera ratio sensitivity to 
!                                                      ! nitrogen concentration                (m2/kg)
  REAL, POINTER, DIMENSION(:)    :: XP_CF_NITRO      ! lethal minimum value of leaf area
!                                                      ! ratio                                 (m2/kg)
  REAL, POINTER, DIMENSION(:)    :: XP_CNA_NITRO     ! nitrogen concentration of active 
!                                                      ! biomass                               (kg/kg)
  REAL, POINTER, DIMENSION(:)    :: XP_BSLAI_NITRO! biomass/LAI ratio from nitrogen 
!                                                      ! decline theory                        (kg/m2)
!
! - soil: primary parameters
!
  REAL, POINTER, DIMENSION(:)      :: XP_RUNOFFB     ! sub-grid surface runoff slope parameter (-)
  REAL, POINTER, DIMENSION(:)      :: XP_WDRAIN      ! continuous drainage parameter           (-)
  REAL, POINTER, DIMENSION(:)      :: XP_TAUICE      ! soil freezing characteristic timescale  (s)
  REAL, POINTER, DIMENSION(:)      :: XP_GAMMAT      ! 'Force-Restore' timescale when using a
!                                                      ! prescribed lower boundary temperature   (1/days)
  REAL, POINTER, DIMENSION(:,:)    :: XP_DG              ! soil layer depth                  (m)
!                                                      ! NOTE: in Force-Restore mode, the 
!                                                      ! uppermost layer thickness is superficial
!                                                      ! and is only explicitly used for soil 
!                                                      ! water phase changes                     (m)
  REAL, POINTER, DIMENSION(:,:)    :: XP_DZG             ! soil layer thicknesses                  (m)
  REAL, POINTER, DIMENSION(:,:)    :: XP_DZDIF           ! distance between consecuative layer mid-points(m)
  INTEGER, POINTER, DIMENSION(:)   :: NK_WG_LAYER        ! Number of soil moisture layers for DIF

  REAL, POINTER, DIMENSION(:)      :: XP_RUNOFFD       ! depth over which sub-grid runoff is
!                                                      ! computed: in Force-Restore this is the
!                                                      ! total soil column ('2-L'), or root zone
!                                                      ! ('3-L'). For the 'DIF' option, it can
!                                                      ! be any depth within soil column         (m)
!
  REAL, POINTER, DIMENSION(:,:)  :: XP_SOILWGHT      ! ISBA-DIF: weights for vertical
!                                                  ! integration of soil water and properties
!
! - soil: Secondary parameters: hydrology
!
  REAL, POINTER, DIMENSION(:)    :: XP_C1SAT       ! 'Force-Restore' C1 coefficient at 
!                                                    ! saturation                              (-)
  REAL, POINTER, DIMENSION(:)    :: XP_C2REF       ! 'Force-Restore' reference value of C2   (-)
  REAL, POINTER, DIMENSION(:,:)  :: XP_C3          ! 'Force-Restore' C3 drainage coefficient (m)
  REAL, POINTER, DIMENSION(:)    :: XP_C4B         ! 'Force-Restore' sub-surface vertical 
!                                                    ! diffusion coefficient (slope parameter) (-)
  REAL, POINTER, DIMENSION(:)    :: XP_C4REF       ! 'Force-Restore' sub-surface vertical 
!                                                    ! diffusion coefficient                   (-)
  REAL, POINTER, DIMENSION(:)    :: XP_ACOEF       ! 'Force-Restore' surface vertical 
!                                                    ! diffusion coefficient                   (-)
  REAL, POINTER, DIMENSION(:)    :: XP_PCOEF       ! 'Force-Restore' surface vertical 
!                                                    ! diffusion coefficient                   (-)
  REAL, POINTER, DIMENSION(:,:)  :: XP_WFC         ! field capacity volumetric water content
!                                                    ! profile                                 (m3/m3)
  REAL, POINTER, DIMENSION(:,:)  :: XP_WWILT       ! wilting point volumetric water content 
!                                                    ! profile                                 (m3/m3)
  REAL, POINTER, DIMENSION(:,:)  :: XP_WSAT        ! porosity profile                        (m3/m3) 
  REAL, POINTER, DIMENSION(:,:)  :: XP_BCOEF       ! soil water CH78 b-parameter             (-)
  REAL, POINTER, DIMENSION(:,:)  :: XP_CONDSAT     ! hydraulic conductivity at saturation    (m/s)
  REAL, POINTER, DIMENSION(:,:)  :: XP_MPOTSAT     ! matric potential at saturation          (m)
!
! - soil: Secondary parameters: thermal 
!
  REAL, POINTER, DIMENSION(:)    :: XP_CGSAT       ! soil thermal inertia coefficient at 
!                                                    ! saturation                              (K m2/J)
  REAL, POINTER, DIMENSION(:,:)  :: XP_HCAPSOIL    ! soil heat capacity                      (J/K/m3)
  REAL, POINTER, DIMENSION(:,:)  :: XP_CONDDRY     ! soil dry thermal conductivity           (W/m/K)
  REAL, POINTER, DIMENSION(:,:)  :: XP_CONDSLD     ! soil solids thermal conductivity        (W/m/K)
  REAL, POINTER, DIMENSION(:)    :: XP_TDEEP       ! prescribed deep soil temperature 
!                                                    ! (optional)                              (K)
! Prognostic variables:
!
! - Snow Cover:
!
  REAL,  POINTER, DIMENSION(:,:) :: XP_SNOWSWE     ! snow (& liq. water) content             (kg/m2)
  REAL,  POINTER, DIMENSION(:,:) :: XP_SNOWHEAT    ! heat content                            (J/m2)
  REAL,  POINTER, DIMENSION(:,:) :: XP_SNOWRHO     ! density                                 (kg m-3)
  REAL,  POINTER, DIMENSION(:,:) :: XP_SNOWGRAN1   ! grain parameter 1                       (-)
  REAL,  POINTER, DIMENSION(:,:) :: XP_SNOWGRAN2   ! grain parameter 2                       (-)
  REAL,  POINTER, DIMENSION(:,:) :: XP_SNOWHIST    ! historical parameter                    (-)
  REAL, POINTER, DIMENSION(:,:)  :: XP_SNOWAGE     ! Snow grain age                          (days)
  REAL,  POINTER, DIMENSION(:)   :: XP_SNOWALB     ! snow tot albedo                         (-)
  REAL,  POINTER, DIMENSION(:)   :: XP_SNOWALBVIS     ! snow VIS albedo                                  (-)
  REAL,  POINTER, DIMENSION(:)   :: XP_SNOWALBNIR     ! snow NIR albedo                                  (-)
  REAL,  POINTER, DIMENSION(:)   :: XP_SNOWALBFIR     ! snow FIR albedo                                  (-)
  REAL,  POINTER, DIMENSION(:)   :: XP_SNOWEMIS    ! snow emissivity (ISBA-ES:3-L)           (-)
!
  REAL,  POINTER, DIMENSION(:)   :: XP_ICE_STO
!
! - Soil and vegetation heat and water:
!
  REAL, POINTER, DIMENSION(:)     :: XP_WR         ! liquid water retained on the
!                                                    ! foliage of the vegetation
!                                                    ! canopy                                  (kg/m2)
  REAL, POINTER, DIMENSION(:,:)   :: XP_TG         ! surface and sub-surface soil 
!                                                    ! temperature profile                     (K)
  REAL, POINTER, DIMENSION(:,:)   :: XP_WG         ! soil volumetric water content profile   (m3/m3)
  REAL, POINTER, DIMENSION(:,:)   :: XP_WGI        ! soil liquid water equivalent volumetric 
!                                                    ! ice content profile                     (m3/m3)
  REAL, POINTER, DIMENSION(:)     :: XP_RESA       ! aerodynamic resistance                  (s/m)
! - For multi-energy balance:
!
  REAL, POINTER, DIMENSION(:)     :: XP_WRV        ! liquid water retained on the foliage
!                                                    ! of the canopy vegetation                (kg/m2)
  REAL, POINTER, DIMENSION(:)     :: XP_WRVN       ! snow retained on the foliage
!                                                    ! of the canopy vegetation                (kg/m2)
  REAL, POINTER, DIMENSION(:)     :: XP_TV         ! canopy vegetation temperature           (K)
  REAL, POINTER, DIMENSION(:)     :: XP_TC         ! canopy air temperature                  (K)
  REAL, POINTER, DIMENSION(:)     :: XP_QC         ! canopy air specific humidity            (kg/kg)
!
  REAL, POINTER, DIMENSION(:)     :: XP_ZF_TALLVEG   ! MEB tall vegetation binary              (-)
  REAL, POINTER, DIMENSION(:)     :: XP_RGLV         ! canopy veg maximum solar radiation
  REAL, POINTER, DIMENSION(:)     :: XP_GAMMAV       ! coefficient for the calculation
!                                                      ! of the canopy veg surface stomatal
  REAL, POINTER, DIMENSION(:)     :: XP_RSMINV       ! canopy veg minimum stomatal resistance  (s/m)
  REAL, POINTER, DIMENSION(:,:)   :: XP_ROOTFRACV    ! canopy veg root fraction profile ('DIF' option)
  REAL, POINTER, DIMENSION(:)     :: XP_WRMAX_CFV    ! canopy veg coefficient for maximum water 
!                                                      ! interception
  REAL, POINTER, DIMENSION(:)     :: XP_LAIV         ! canopy veg Leaf Area Index              (m2/m2)
  REAL, POINTER, DIMENSION(:)     :: XP_Z0V          ! canopy veg roughness length             (m)
  REAL, POINTER, DIMENSION(:)     :: XP_H_VEG        ! vegetation height                       (m)
  REAL, POINTER, DIMENSION(:)     :: XP_GNDLITTER    ! ground litter cover                     (-)
  REAL, POINTER, DIMENSION(:)     :: XP_Z0LITTER     ! ground litter roughness length          (m)
!
! - Vegetation: Ags Prognos
!
  REAL, POINTER, DIMENSION(:)     :: XP_FWTD       ! grid-cell fraction of water table to rise
  REAL, POINTER, DIMENSION(:)     :: XP_WTD        ! water table depth                  (m)
!
! - Vegetation: Ags Prognostic (YPHOTO = 'LAI', 'LST', 'NIT', 'NCB') or prescribed (YPHOTO = 'NON', 'AGS', 'AST')
!
  REAL, POINTER, DIMENSION(:)     :: XP_LAI        ! Leaf Area Index                         (m2/m2)
!
! - Vegetation: Ags Prognostic (YPHOTO = 'AGS', 'LAI', 'AST', 'LST', 'NIT', 'NCB')
!
  REAL, POINTER, DIMENSION(:)     :: XP_AN         ! net CO2 assimilation                    (mg/m2/s)
  REAL, POINTER, DIMENSION(:)     :: XP_ANDAY      ! daily net CO2 assimilation              (mg/m2)
  REAL, POINTER, DIMENSION(:)     :: XP_ANFM       ! maximum leaf assimilation               (mg/m2/s)
  REAL, POINTER, DIMENSION(:)     :: XP_LE         ! evapotranspiration                      (W/m2)
  REAL, POINTER, DIMENSION(:)     :: XP_LEI        ! sublimation                             (W/m2)
  REAL, POINTER, DIMENSION(:)     :: XP_FAPARC     ! FAPAR of vegetation (cumul)
  REAL, POINTER, DIMENSION(:)     :: XP_FAPIRC     ! FAPIR of vegetation (cumul)
  REAL, POINTER, DIMENSION(:)     :: XP_LAI_EFFC   ! effective LAI (cumul)
  REAL, POINTER, DIMENSION(:)     :: XP_MUS        ! 
!
! - Vegetation: Ags Prognostic (YPHOTO = 'NIT', 'NCB')
!
  REAL, POINTER, DIMENSION(:,:)   :: XP_RESP_BIOMASS    ! daily cumulated respiration of 
!                                                         ! biomass                            (kg/m2/s)
  REAL, POINTER, DIMENSION(:,:)   :: XP_BIOMASS         ! biomass of previous day            (kg/m2) 
  REAL, POINTER, DIMENSION(:,:)   :: XP_INCREASE        ! biomass increase                   (kg/m2/day)
!
! - Soil carbon (ISBA-CC, YRESPSL = 'CNT')
!
  REAL, POINTER, DIMENSION(:,:,:)   :: XP_LITTER            ! litter pools                       (gC/m2)
  REAL, POINTER, DIMENSION(:,:)     :: XP_SOILCARB          ! soil carbon pools                  (gC/m2) 
  REAL, POINTER, DIMENSION(:,:)     :: XP_LIGNIN_STRUC      ! ratio Lignin/Carbon in structural
!                                                         litter                               (gC/m2)
!
  REAL, POINTER, DIMENSION(:,:)   :: XP_TURNOVER        ! turnover rates from biomass to litter (gC/m2/s)
!
! - Irrigation
!
  LOGICAL, POINTER, DIMENSION(:)     :: XP_LIRRIGATE    ! high level switch for irrigation
  LOGICAL, POINTER, DIMENSION(:)     :: XP_LIRRIDAY     ! flag used for daily irrigation stage
  REAL, POINTER, DIMENSION(:)        :: XP_THRESHOLD    ! threshold for stages
  REAL, POINTER, DIMENSION(:)        :: XP_WATSUP       ! water supply
  REAL, POINTER, DIMENSION(:)        :: XP_IRRIG        ! fraction of irrigated vegetation
  TYPE(DATE_TIME), DIMENSION(:), POINTER :: TP_SEED         ! seeding date
  TYPE(DATE_TIME), DIMENSION(:), POINTER :: TP_REAP         ! reaping date
!                                                         ! previous day                         (kg/m2)
! - SGH scheme
!
  REAL, POINTER, DIMENSION(:)        :: XP_D_ICE     !depth of the soil column for the calculation
!                                                       of the frozen soil fraction (m)
  REAL, POINTER, DIMENSION(:)        :: XP_KSAT_ICE  !hydraulic conductivity at saturation
!                                                       over frozen area (m s-1)
  REAL, POINTER, DIMENSION(:)        :: XP_FSAT      !Topmodel saturated fraction
  REAL, POINTER, DIMENSION(:)        :: XP_MUF       !Rainfall surface fraction 
  REAL, POINTER, DIMENSION(:,:)      :: XP_TOPQS     !Topmodel baseflow by layer (m s-1)
!
! - Courant time step properties
!
  REAL, POINTER, DIMENSION(:)        :: XP_PSN       ! fraction of the grid covered by snow          (-)
  REAL, POINTER, DIMENSION(:)        :: XP_PSNG      ! fraction of the the bare ground covered by snow (-)
  REAL, POINTER, DIMENSION(:)        :: XP_PSNV      ! fraction of the the vegetation covered by snow(-)
  REAL, POINTER, DIMENSION(:)        :: XP_PSNV_A    ! fraction of the the vegetation covered by snow for EBA scheme(-)
  REAL, POINTER, DIMENSION(:,:)      :: XP_DIR_ALB_WITH_SNOW ! Total direct albedo
  REAL, POINTER, DIMENSION(:,:)      :: XP_SCA_ALB_WITH_SNOW ! Total diffuse albedo
!
! - Flood scheme
!
  REAL, POINTER, DIMENSION(:)        :: XP_ALBF
  REAL, POINTER, DIMENSION(:)        :: XP_EMISF
!
  REAL, POINTER, DIMENSION(:)        :: XP_FF        ! flood fraction over the surface
  REAL, POINTER, DIMENSION(:)        :: XP_FFG       ! flood fraction over the ground
  REAL, POINTER, DIMENSION(:)        :: XP_FFV       ! flood fraction over the vegetation
  REAL, POINTER, DIMENSION(:)        :: XP_FFROZEN   ! fraction of frozen flood
  REAL, POINTER, DIMENSION(:)        :: XP_FFLOOD  ! Grdi-cell flood fraction           (-)
  REAL, POINTER, DIMENSION(:)        :: XP_PIFLOOD ! Floodplains potential infiltration (kg/m2/s)
!
  REAL, POINTER, DIMENSION(:)        :: XP_CPS, XP_LVTT, XP_LSTT

END TYPE PACK_ISBA_t
!
!-------------------------------------------------------------------------------
!
TYPE(PACK_ISBA_t), ALLOCATABLE, TARGET, SAVE :: PACK_ISBA_MODEL(:)
!
INTEGER, POINTER :: NSIZE_LSIMPLE
!$OMP THREADPRIVATE(NSIZE_LSIMPLE)
INTEGER, POINTER :: NSIZE_L0
!$OMP THREADPRIVATE(NSIZE_L0)
INTEGER, POINTER :: NSIZE_NSIMPLE
!$OMP THREADPRIVATE(NSIZE_NSIMPLE)
INTEGER, POINTER :: NSIZE_N0
!$OMP THREADPRIVATE(NSIZE_N0)
INTEGER, POINTER :: NSIZE_TSIMPLE
!$OMP THREADPRIVATE(NSIZE_TSIMPLE)
INTEGER, POINTER :: NSIZE_T0
!$OMP THREADPRIVATE(NSIZE_T0)
INTEGER, POINTER :: NSIZE_SIMPLE
!$OMP THREADPRIVATE(NSIZE_SIMPLE)
INTEGER, POINTER :: NSIZE_GROUND
!$OMP THREADPRIVATE(NSIZE_GROUND)
INTEGER, POINTER :: NSIZE_VEGTYPE
!$OMP THREADPRIVATE(NSIZE_VEGTYPE)
INTEGER, POINTER :: NSIZE_TG
!$OMP THREADPRIVATE(NSIZE_TG)
INTEGER, POINTER :: NSIZE_SNOW
!$OMP THREADPRIVATE(NSIZE_SNOW)
INTEGER, POINTER :: NSIZE_ALB
!$OMP THREADPRIVATE(NSIZE_ALB)
INTEGER, POINTER :: NSIZE_2
!$OMP THREADPRIVATE(NSIZE_2)
INTEGER, POINTER :: NSIZE_BIOMASS
!$OMP THREADPRIVATE(NSIZE_BIOMASS)
INTEGER, POINTER :: NSIZE_SOILCARB
!$OMP THREADPRIVATE(NSIZE_SOILCARB)
INTEGER, POINTER :: NSIZE_LITTLEVS
!$OMP THREADPRIVATE(NSIZE_LITTLEVS)
INTEGER, POINTER :: NSIZE_LITTER
!$OMP THREADPRIVATE(NSIZE_LITTER)
INTEGER, POINTER :: NSIZE_0
!$OMP THREADPRIVATE(NSIZE_0)
INTEGER, POINTER :: NSIZE_00
!$OMP THREADPRIVATE(NSIZE_00)
INTEGER, POINTER :: NSIZE_000
!$OMP THREADPRIVATE(NSIZE_000)
INTEGER, POINTER :: NSIZE_01
!$OMP THREADPRIVATE(NSIZE_01)
LOGICAL, POINTER, DIMENSION(:,:) :: LBLOCK_SIMPLE 
!$OMP THREADPRIVATE(LBLOCK_SIMPLE)
LOGICAL, POINTER, DIMENSION(:,:) :: LBLOCK_0
!$OMP THREADPRIVATE(LBLOCK_0)
INTEGER, POINTER, DIMENSION(:,:) :: NBLOCK_SIMPLE 
!$OMP THREADPRIVATE(NBLOCK_SIMPLE)
INTEGER, POINTER, DIMENSION(:,:) :: NBLOCK_0
!$OMP THREADPRIVATE(NBLOCK_0)
TYPE(DATE_TIME), POINTER, DIMENSION(:,:) :: TBLOCK_SIMPLE
!$OMP THREADPRIVATE(TBLOCK_SIMPLE)
TYPE(DATE_TIME), POINTER, DIMENSION(:,:) :: TBLOCK_0
!$OMP THREADPRIVATE(TBLOCK_0)
REAL, POINTER, DIMENSION(:,:) :: XBLOCK_SIMPLE
!$OMP THREADPRIVATE(XBLOCK_SIMPLE)
REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_GROUND
!$OMP THREADPRIVATE(XBLOCK_GROUND)
REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_VEGTYPE
!$OMP THREADPRIVATE(XBLOCK_VEGTYPE)
REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_TG
!$OMP THREADPRIVATE(XBLOCK_TG)
REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_SNOW
!$OMP THREADPRIVATE(XBLOCK_SNOW)
REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_ALB
!$OMP THREADPRIVATE(XBLOCK_ALB)
REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_2
!$OMP THREADPRIVATE(XBLOCK_2)
REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_BIOMASS
!$OMP THREADPRIVATE(XBLOCK_BIOMASS)
REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_SOILCARB
!$OMP THREADPRIVATE(XBLOCK_SOILCARB)
REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_LITTLEVS
!$OMP THREADPRIVATE(XBLOCK_LITTLEVS)
REAL, POINTER, DIMENSION(:,:,:,:) :: XBLOCK_LITTER
!$OMP THREADPRIVATE(XBLOCK_LITTER)
REAL, POINTER, DIMENSION(:,:) :: XBLOCK_0
!$OMP THREADPRIVATE(XBLOCK_0)
REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_00
!$OMP THREADPRIVATE(XBLOCK_00)
REAL, POINTER, DIMENSION(:,:,:,:) :: XBLOCK_000
!$OMP THREADPRIVATE(XBLOCK_000)
REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_01
!$OMP THREADPRIVATE(XBLOCK_01)
!
REAL, POINTER, DIMENSION(:,:)  :: XP_VEGTYPE_PATCH 
!$OMP THREADPRIVATE(XP_VEGTYPE_PATCH)
!
REAL, POINTER, DIMENSION(:)    :: XP_SSO_SLOPE     
!$OMP THREADPRIVATE(XP_SSO_SLOPE)
REAL, POINTER, DIMENSION(:)    :: XP_LAT           
!$OMP THREADPRIVATE(XP_LAT)
REAL, POINTER, DIMENSION(:)    :: XP_LON           
!$OMP THREADPRIVATE(XP_LON)
!
REAL, DIMENSION(:), POINTER :: XP_AOSIP,XP_AOSIM,XP_AOSJP,XP_AOSJM
!$OMP THREADPRIVATE(XP_AOSIP, XP_AOSIM, XP_AOSJP, XP_AOSJM)
!
REAL, DIMENSION(:), POINTER :: XP_HO2IP,XP_HO2IM,XP_HO2JP,XP_HO2JM
!$OMP THREADPRIVATE(XP_HO2IP, XP_HO2IM, XP_HO2JP, XP_HO2JM)
!
REAL, DIMENSION(:), POINTER :: XP_Z0EFFIP,XP_Z0EFFIM,XP_Z0EFFJP,XP_Z0EFFJM
!$OMP THREADPRIVATE(XP_Z0EFFIP, XP_Z0EFFIM, XP_Z0EFFJP, XP_Z0EFFJM)
!
REAL, POINTER, DIMENSION(:) :: XP_Z0REL         
!$OMP THREADPRIVATE(XP_Z0REL)
!
REAL, POINTER, DIMENSION(:,:) :: XP_CLAY         
!$OMP THREADPRIVATE(XP_CLAY)
REAL, POINTER, DIMENSION(:,:) :: XP_SAND         
!$OMP THREADPRIVATE(XP_SAND)
!
REAL, POINTER, DIMENSION(:) :: XP_ALBNIR_DRY     
!$OMP THREADPRIVATE(XP_ALBNIR_DRY)
REAL, POINTER, DIMENSION(:) :: XP_ALBVIS_DRY     
!$OMP THREADPRIVATE(XP_ALBVIS_DRY)
REAL, POINTER, DIMENSION(:) :: XP_ALBUV_DRY      
!$OMP THREADPRIVATE(XP_ALBUV_DRY)
REAL, POINTER, DIMENSION(:) :: XP_ALBNIR_WET     
!$OMP THREADPRIVATE(XP_ALBNIR_WET)
REAL, POINTER, DIMENSION(:) :: XP_ALBVIS_WET     
!$OMP THREADPRIVATE(XP_ALBVIS_WET)
REAL, POINTER, DIMENSION(:) :: XP_ALBUV_WET      
!$OMP THREADPRIVATE(XP_ALBUV_WET)
REAL, POINTER, DIMENSION(:) :: XP_ALBNIR_SOIL    
!$OMP THREADPRIVATE(XP_ALBNIR_SOIL)
REAL, POINTER, DIMENSION(:) :: XP_ALBVIS_SOIL    
!$OMP THREADPRIVATE(XP_ALBVIS_SOIL)
REAL, POINTER, DIMENSION(:) :: XP_ALBUV_SOIL     
!$OMP THREADPRIVATE(XP_ALBUV_SOIL)
!
REAL, POINTER, DIMENSION(:) :: XP_Z0_O_Z0H       
!$OMP THREADPRIVATE(XP_Z0_O_Z0H)
REAL, POINTER, DIMENSION(:) :: XP_ALBNIR         
!$OMP THREADPRIVATE(XP_ALBNIR)
REAL, POINTER, DIMENSION(:) :: XP_ALBVIS         
!$OMP THREADPRIVATE(XP_ALBVIS)
REAL, POINTER, DIMENSION(:) :: XP_ALBUV          
!$OMP THREADPRIVATE(XP_ALBUV)
REAL, POINTER, DIMENSION(:) :: XP_EMIS           
!$OMP THREADPRIVATE(XP_EMIS)
REAL, POINTER, DIMENSION(:) :: XP_Z0             
!$OMP THREADPRIVATE(XP_Z0)
!
REAL, POINTER, DIMENSION(:) :: XP_ALBNIR_VEG     
!$OMP THREADPRIVATE(XP_ALBNIR_VEG)
REAL, POINTER, DIMENSION(:) :: XP_ALBVIS_VEG     
!$OMP THREADPRIVATE(XP_ALBVIS_VEG)
REAL, POINTER, DIMENSION(:) :: XP_ALBUV_VEG      
!$OMP THREADPRIVATE(XP_ALBUV_VEG)
!
REAL, POINTER, DIMENSION(:) :: XP_VEG            
!$OMP THREADPRIVATE(XP_VEG)
REAL, POINTER, DIMENSION(:) :: XP_WRMAX_CF       
!$OMP THREADPRIVATE(XP_WRMAX_CF)
REAL, POINTER, DIMENSION(:) :: XP_RSMIN          
!$OMP THREADPRIVATE(XP_RSMIN)
REAL, POINTER, DIMENSION(:) :: XP_GAMMA          
!$OMP THREADPRIVATE(XP_GAMMA)
REAL, POINTER, DIMENSION(:) :: XP_CV             
!$OMP THREADPRIVATE(XP_CV)
REAL, POINTER, DIMENSION(:) :: XP_RGL            
!$OMP THREADPRIVATE(XP_RGL)
REAL, POINTER, DIMENSION(:,:) :: XP_ROOTFRAC     
!$OMP THREADPRIVATE(XP_ROOTFRAC)
!
REAL, POINTER, DIMENSION(:)    :: XP_BSLAI      
!$OMP THREADPRIVATE(XP_BSLAI)
REAL, POINTER, DIMENSION(:)    :: XP_LAIMIN     
!$OMP THREADPRIVATE(XP_LAIMIN)
REAL, POINTER, DIMENSION(:)    :: XP_SEFOLD     
!$OMP THREADPRIVATE(XP_SEFOLD)
REAL, POINTER, DIMENSION(:)    :: XP_H_TREE     
!$OMP THREADPRIVATE(XP_H_TREE)
REAL, POINTER, DIMENSION(:)    :: XP_ANF        
!$OMP THREADPRIVATE(XP_ANF)
REAL, POINTER, DIMENSION(:)    :: XP_ANMAX      
!$OMP THREADPRIVATE(XP_ANMAX)
REAL, POINTER, DIMENSION(:)    :: XP_FZERO      
!$OMP THREADPRIVATE(XP_FZERO)
REAL, POINTER, DIMENSION(:)    :: XP_EPSO       
!$OMP THREADPRIVATE(XP_EPSO)
REAL, POINTER, DIMENSION(:)    :: XP_GAMM       
!$OMP THREADPRIVATE(XP_GAMM)
REAL, POINTER, DIMENSION(:)    :: XP_QDGAMM     
!$OMP THREADPRIVATE(XP_QDGAMM)
REAL, POINTER, DIMENSION(:)    :: XP_GMES       
!$OMP THREADPRIVATE(XP_GMES)
REAL, POINTER, DIMENSION(:)    :: XP_RE25       
!$OMP THREADPRIVATE(XP_RE25)
REAL, POINTER, DIMENSION(:)    :: XP_QDGMES     
!$OMP THREADPRIVATE(XP_QDGMES)
REAL, POINTER, DIMENSION(:)    :: XP_T1GMES     
!$OMP THREADPRIVATE(XP_T1GMES)
REAL, POINTER, DIMENSION(:)    :: XP_T2GMES     
!$OMP THREADPRIVATE(XP_T2GMES)
REAL, POINTER, DIMENSION(:)    :: XP_AMAX       
!$OMP THREADPRIVATE(XP_AMAX)
REAL, POINTER, DIMENSION(:)    :: XP_QDAMAX     
!$OMP THREADPRIVATE(XP_QDAMAX)
REAL, POINTER, DIMENSION(:)    :: XP_T1AMAX     
!$OMP THREADPRIVATE(XP_T1AMAX)
REAL, POINTER, DIMENSION(:)    :: XP_T2AMAX     
!$OMP THREADPRIVATE(XP_T2AMAX)
!
LOGICAL, POINTER, DIMENSION(:) :: LP_STRESS      
!$OMP THREADPRIVATE(LP_STRESS)
REAL, POINTER, DIMENSION(:)    :: XP_F2I         
!$OMP THREADPRIVATE(XP_F2I)
REAL, POINTER, DIMENSION(:)    :: XP_GC          
!$OMP THREADPRIVATE(XP_GC)
REAL, POINTER, DIMENSION(:)    :: XP_AH          
!$OMP THREADPRIVATE(XP_AH)
REAL, POINTER, DIMENSION(:)    :: XP_BH          
!$OMP THREADPRIVATE(XP_BH)
REAL, POINTER, DIMENSION(:)    :: XP_TAU_WOOD    
!$OMP THREADPRIVATE(XP_TAU_WOOD)
REAL, POINTER, DIMENSION(:)    :: XP_DMAX        
!$OMP THREADPRIVATE(XP_DMAX)
!
REAL, POINTER, DIMENSION(:)    :: XP_CE_NITRO      
!$OMP THREADPRIVATE(XP_CE_NITRO)
REAL, POINTER, DIMENSION(:)    :: XP_CF_NITRO      
!$OMP THREADPRIVATE(XP_CF_NITRO)
REAL, POINTER, DIMENSION(:)    :: XP_CNA_NITRO     
!$OMP THREADPRIVATE(XP_CNA_NITRO)
REAL, POINTER, DIMENSION(:)    :: XP_BSLAI_NITRO
!$OMP THREADPRIVATE(XP_BSLAI_NITRO)
!
REAL, POINTER, DIMENSION(:)      :: XP_RUNOFFB     
!$OMP THREADPRIVATE(XP_RUNOFFB)
REAL, POINTER, DIMENSION(:)      :: XP_WDRAIN      
!$OMP THREADPRIVATE(XP_WDRAIN)
REAL, POINTER, DIMENSION(:)      :: XP_TAUICE      
!$OMP THREADPRIVATE(XP_TAUICE)
REAL, POINTER, DIMENSION(:)      :: XP_GAMMAT      
!$OMP THREADPRIVATE(XP_GAMMAT)
REAL, POINTER, DIMENSION(:,:)    :: XP_DG              
!$OMP THREADPRIVATE(XP_DG)
REAL, POINTER, DIMENSION(:,:)    :: XP_DZG             
!$OMP THREADPRIVATE(XP_DZG)
REAL, POINTER, DIMENSION(:,:)    :: XP_DZDIF           
!$OMP THREADPRIVATE(XP_DZDIF)
INTEGER, POINTER, DIMENSION(:)   :: NK_WG_LAYER        
!$OMP THREADPRIVATE(NK_WG_LAYER)
!
REAL, POINTER, DIMENSION(:)      :: XP_RUNOFFD       
!$OMP THREADPRIVATE(XP_RUNOFFD)
!
REAL, POINTER, DIMENSION(:,:)  :: XP_SOILWGHT      
!$OMP THREADPRIVATE(XP_SOILWGHT)
!
REAL, POINTER, DIMENSION(:)    :: XP_C1SAT       
!$OMP THREADPRIVATE(XP_C1SAT)
REAL, POINTER, DIMENSION(:)    :: XP_C2REF       
!$OMP THREADPRIVATE(XP_C2REF)
REAL, POINTER, DIMENSION(:,:)  :: XP_C3          
!$OMP THREADPRIVATE(XP_C3)
REAL, POINTER, DIMENSION(:)    :: XP_C4B         
!$OMP THREADPRIVATE(XP_C4B)
REAL, POINTER, DIMENSION(:)    :: XP_C4REF       
!$OMP THREADPRIVATE(XP_C4REF)
REAL, POINTER, DIMENSION(:)    :: XP_ACOEF       
!$OMP THREADPRIVATE(XP_ACOEF)
REAL, POINTER, DIMENSION(:)    :: XP_PCOEF       
!$OMP THREADPRIVATE(XP_PCOEF)
REAL, POINTER, DIMENSION(:,:)  :: XP_WFC         
!$OMP THREADPRIVATE(XP_WFC)
!
REAL, POINTER, DIMENSION(:,:)  :: XP_WWILT       
!$OMP THREADPRIVATE(XP_WWILT)
REAL, POINTER, DIMENSION(:,:)  :: XP_WSAT        
!$OMP THREADPRIVATE(XP_WSAT)
REAL, POINTER, DIMENSION(:,:)  :: XP_BCOEF       
!$OMP THREADPRIVATE(XP_BCOEF)
REAL, POINTER, DIMENSION(:,:)  :: XP_CONDSAT     
!$OMP THREADPRIVATE(XP_CONDSAT)
REAL, POINTER, DIMENSION(:,:)  :: XP_MPOTSAT     
!$OMP THREADPRIVATE(XP_MPOTSAT)
!
REAL, POINTER, DIMENSION(:)    :: XP_CGSAT       
!$OMP THREADPRIVATE(XP_CGSAT)
REAL, POINTER, DIMENSION(:,:)  :: XP_HCAPSOIL    
!$OMP THREADPRIVATE(XP_HCAPSOIL)
REAL, POINTER, DIMENSION(:,:)  :: XP_CONDDRY     
!$OMP THREADPRIVATE(XP_CONDDRY)
REAL, POINTER, DIMENSION(:,:)  :: XP_CONDSLD     
!$OMP THREADPRIVATE(XP_CONDSLD)
REAL, POINTER, DIMENSION(:)    :: XP_TDEEP       
!$OMP THREADPRIVATE(XP_TDEEP)
!
REAL,  POINTER, DIMENSION(:,:) :: XP_SNOWSWE     
!$OMP THREADPRIVATE(XP_SNOWSWE)
REAL,  POINTER, DIMENSION(:,:) :: XP_SNOWHEAT    
!$OMP THREADPRIVATE(XP_SNOWHEAT)
REAL,  POINTER, DIMENSION(:,:) :: XP_SNOWRHO     
!$OMP THREADPRIVATE(XP_SNOWRHO)
REAL,  POINTER, DIMENSION(:,:) :: XP_SNOWGRAN1   
!$OMP THREADPRIVATE(XP_SNOWGRAN1)
REAL,  POINTER, DIMENSION(:,:) :: XP_SNOWGRAN2   
!$OMP THREADPRIVATE(XP_SNOWGRAN2)
REAL,  POINTER, DIMENSION(:,:) :: XP_SNOWHIST    
!$OMP THREADPRIVATE(XP_SNOWHIST)
REAL, POINTER, DIMENSION(:,:)  :: XP_SNOWAGE     
!$OMP THREADPRIVATE(XP_SNOWAGE)
REAL,  POINTER, DIMENSION(:)   :: XP_SNOWALB     
!$OMP THREADPRIVATE(XP_SNOWALB)
REAL,  POINTER, DIMENSION(:)   :: XP_SNOWALBVIS     
!$OMP THREADPRIVATE(XP_SNOWALBVIS)
REAL,  POINTER, DIMENSION(:)   :: XP_SNOWALBNIR     
!$OMP THREADPRIVATE(XP_SNOWALBNIR)
REAL,  POINTER, DIMENSION(:)   :: XP_SNOWALBFIR     
!$OMP THREADPRIVATE(XP_SNOWALBFIR)
REAL,  POINTER, DIMENSION(:)   :: XP_SNOWEMIS    
!$OMP THREADPRIVATE(XP_SNOWEMIS)
!
REAL,  POINTER, DIMENSION(:)   :: XP_ICE_STO
!$OMP THREADPRIVATE(XP_ICE_STO)
!
REAL, POINTER, DIMENSION(:)     :: XP_WR         
!$OMP THREADPRIVATE(XP_WR)
REAL, POINTER, DIMENSION(:,:)   :: XP_TG         
!$OMP THREADPRIVATE(XP_TG)
REAL, POINTER, DIMENSION(:,:)   :: XP_WG         
!$OMP THREADPRIVATE(XP_WG)
REAL, POINTER, DIMENSION(:,:)   :: XP_WGI        
!$OMP THREADPRIVATE(XP_WGI)
REAL, POINTER, DIMENSION(:)     :: XP_RESA       
!$OMP THREADPRIVATE(XP_RESA)
!
REAL, POINTER, DIMENSION(:)     :: XP_WRV        
!$OMP THREADPRIVATE(XP_WRV)
REAL, POINTER, DIMENSION(:)     :: XP_WRVN       
!$OMP THREADPRIVATE(XP_WRVN)
REAL, POINTER, DIMENSION(:)     :: XP_TV         
!$OMP THREADPRIVATE(XP_TV)
REAL, POINTER, DIMENSION(:)     :: XP_TC         
!$OMP THREADPRIVATE(XP_TC)
REAL, POINTER, DIMENSION(:)     :: XP_QC         
!$OMP THREADPRIVATE(XP_QC)
!
REAL, POINTER, DIMENSION(:)     :: XP_ZF_TALLVEG   
!$OMP THREADPRIVATE(XP_ZF_TALLVEG)
REAL, POINTER, DIMENSION(:)     :: XP_RGLV         
!$OMP THREADPRIVATE(XP_RGLV)
REAL, POINTER, DIMENSION(:)     :: XP_GAMMAV       
!$OMP THREADPRIVATE(XP_GAMMAV)
REAL, POINTER, DIMENSION(:)     :: XP_RSMINV       
!$OMP THREADPRIVATE(XP_RSMINV)
REAL, POINTER, DIMENSION(:,:)   :: XP_ROOTFRACV    
!$OMP THREADPRIVATE(XP_ROOTFRACV)
REAL, POINTER, DIMENSION(:)     :: XP_WRMAX_CFV    
!$OMP THREADPRIVATE(XP_WRMAX_CFV)
REAL, POINTER, DIMENSION(:)     :: XP_LAIV         
!$OMP THREADPRIVATE(XP_LAIV)
REAL, POINTER, DIMENSION(:)     :: XP_Z0V          
!$OMP THREADPRIVATE(XP_Z0V)
REAL, POINTER, DIMENSION(:)     :: XP_H_VEG        
!$OMP THREADPRIVATE(XP_H_VEG)
REAL, POINTER, DIMENSION(:)     :: XP_GNDLITTER    
!$OMP THREADPRIVATE(XP_GNDLITTER)
REAL, POINTER, DIMENSION(:)     :: XP_Z0LITTER     
!$OMP THREADPRIVATE(XP_Z0LITTER)
!
REAL, POINTER, DIMENSION(:)     :: XP_FWTD       
!$OMP THREADPRIVATE(XP_FWTD)
REAL, POINTER, DIMENSION(:)     :: XP_WTD        
!$OMP THREADPRIVATE(XP_WTD)
!
REAL, POINTER, DIMENSION(:)     :: XP_LAI        
!$OMP THREADPRIVATE(XP_LAI)
!
REAL, POINTER, DIMENSION(:)     :: XP_AN         
!$OMP THREADPRIVATE(XP_AN)
REAL, POINTER, DIMENSION(:)     :: XP_ANDAY      
!$OMP THREADPRIVATE(XP_ANDAY)
REAL, POINTER, DIMENSION(:)     :: XP_ANFM       
!$OMP THREADPRIVATE(XP_ANFM)
REAL, POINTER, DIMENSION(:)     :: XP_LE         
!$OMP THREADPRIVATE(XP_LE)
REAL, POINTER, DIMENSION(:)     :: XP_LEI        
!$OMP THREADPRIVATE(XP_LEI)
REAL, POINTER, DIMENSION(:)     :: XP_FAPARC     
!$OMP THREADPRIVATE(XP_FAPARC)
REAL, POINTER, DIMENSION(:)     :: XP_FAPIRC     
!$OMP THREADPRIVATE(XP_FAPIRC)
REAL, POINTER, DIMENSION(:)     :: XP_LAI_EFFC   
!$OMP THREADPRIVATE(XP_LAI_EFFC)
REAL, POINTER, DIMENSION(:)     :: XP_MUS        
!$OMP THREADPRIVATE(XP_MUS)
!
REAL, POINTER, DIMENSION(:,:)   :: XP_RESP_BIOMASS    
!$OMP THREADPRIVATE(XP_RESP_BIOMASS)
REAL, POINTER, DIMENSION(:,:)   :: XP_BIOMASS         
!$OMP THREADPRIVATE(XP_BIOMASS)
REAL, POINTER, DIMENSION(:,:)   :: XP_INCREASE        
!$OMP THREADPRIVATE(XP_INCREASE)
!
REAL, POINTER, DIMENSION(:,:,:)   :: XP_LITTER            
!$OMP THREADPRIVATE(XP_LITTER)
REAL, POINTER, DIMENSION(:,:)     :: XP_SOILCARB          
!$OMP THREADPRIVATE(XP_SOILCARB)
REAL, POINTER, DIMENSION(:,:)     :: XP_LIGNIN_STRUC      
!$OMP THREADPRIVATE(XP_LIGNIN_STRUC)
!
REAL, POINTER, DIMENSION(:,:)   :: XP_TURNOVER        
!$OMP THREADPRIVATE(XP_TURNOVER)
!
LOGICAL, POINTER, DIMENSION(:)     :: XP_LIRRIGATE    
!$OMP THREADPRIVATE(XP_LIRRIGATE)
LOGICAL, POINTER, DIMENSION(:)     :: XP_LIRRIDAY     
!$OMP THREADPRIVATE(XP_LIRRIDAY)
REAL, POINTER, DIMENSION(:)        :: XP_THRESHOLD    
!$OMP THREADPRIVATE(XP_THRESHOLD)
REAL, POINTER, DIMENSION(:)        :: XP_WATSUP       
!$OMP THREADPRIVATE(XP_WATSUP)
REAL, POINTER, DIMENSION(:)        :: XP_IRRIG        
!$OMP THREADPRIVATE(XP_IRRIG)
TYPE(DATE_TIME), DIMENSION(:), POINTER :: TP_SEED         
!$OMP THREADPRIVATE(TP_SEED)
TYPE(DATE_TIME), DIMENSION(:), POINTER :: TP_REAP         
!$OMP THREADPRIVATE(TP_REAP)
!
REAL, POINTER, DIMENSION(:)        :: XP_D_ICE     !depth of the soil column for the calculation
!$OMP THREADPRIVATE(XP_D_ICE)
REAL, POINTER, DIMENSION(:)        :: XP_KSAT_ICE  !hydraulic conductivity at saturation
!$OMP THREADPRIVATE(XP_KSAT_ICE)
REAL, POINTER, DIMENSION(:)        :: XP_FSAT      !Topmodel saturated fraction
!$OMP THREADPRIVATE(XP_FSAT)
REAL, POINTER, DIMENSION(:)        :: XP_MUF       !Rainfall surface fraction 
!$OMP THREADPRIVATE(XP_MUF)
REAL, POINTER, DIMENSION(:,:)      :: XP_TOPQS     !Topmodel baseflow by layer (m s-1)
!$OMP THREADPRIVATE(XP_TOPQS)
!
REAL, POINTER, DIMENSION(:)        :: XP_PSN       
!$OMP THREADPRIVATE(XP_PSN)
REAL, POINTER, DIMENSION(:)        :: XP_PSNG      
!$OMP THREADPRIVATE(XP_PSNG)
REAL, POINTER, DIMENSION(:)        :: XP_PSNV      
!$OMP THREADPRIVATE(XP_PSNV)
REAL, POINTER, DIMENSION(:)        :: XP_PSNV_A    
!$OMP THREADPRIVATE(XP_PSNV_A)
REAL, POINTER, DIMENSION(:,:)      :: XP_DIR_ALB_WITH_SNOW 
!$OMP THREADPRIVATE(XP_DIR_ALB_WITH_SNOW)
REAL, POINTER, DIMENSION(:,:)      :: XP_SCA_ALB_WITH_SNOW 
!$OMP THREADPRIVATE(XP_SCA_ALB_WITH_SNOW)
!
REAL, POINTER, DIMENSION(:)        :: XP_ALBF
!$OMP THREADPRIVATE(XP_ALBF)
REAL, POINTER, DIMENSION(:)        :: XP_EMISF
!$OMP THREADPRIVATE(XP_EMISF)
!
REAL, POINTER, DIMENSION(:)        :: XP_FF        
!$OMP THREADPRIVATE(XP_FF)
REAL, POINTER, DIMENSION(:)        :: XP_FFG       
!$OMP THREADPRIVATE(XP_FFG)
REAL, POINTER, DIMENSION(:)        :: XP_FFV       
!$OMP THREADPRIVATE(XP_FFV)
REAL, POINTER, DIMENSION(:)        :: XP_FFROZEN   
!$OMP THREADPRIVATE(XP_FFROZEN)
REAL, POINTER, DIMENSION(:)        :: XP_FFLOOD  
!$OMP THREADPRIVATE(XP_FFLOOD)
REAL, POINTER, DIMENSION(:)        :: XP_PIFLOOD 
!$OMP THREADPRIVATE(XP_PIFLOOD)
!
REAL, POINTER, DIMENSION(:)        :: XP_CPS, XP_LVTT, XP_LSTT
!$OMP THREADPRIVATE(XP_CPS, XP_LVTT, XP_LSTT)
!
!------------------------------------------------------------------------------
!
CONTAINS

SUBROUTINE PACK_ISBA_GOTO_MODEL(KFROM,KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
!
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Save current state for allocated arrays
IF (LKFROM) THEN
  PACK_ISBA_MODEL(KFROM)%LBLOCK_SIMPLE=>LBLOCK_SIMPLE
  PACK_ISBA_MODEL(KFROM)%LBLOCK_0=>LBLOCK_0
  PACK_ISBA_MODEL(KFROM)%NBLOCK_SIMPLE=>NBLOCK_SIMPLE
  PACK_ISBA_MODEL(KFROM)%NBLOCK_0=>NBLOCK_0
  PACK_ISBA_MODEL(KFROM)%TBLOCK_SIMPLE=>TBLOCK_SIMPLE
  PACK_ISBA_MODEL(KFROM)%TBLOCK_0=>TBLOCK_0
  PACK_ISBA_MODEL(KFROM)%XBLOCK_SIMPLE=>XBLOCK_SIMPLE
  PACK_ISBA_MODEL(KFROM)%XBLOCK_GROUND=>XBLOCK_GROUND
  PACK_ISBA_MODEL(KFROM)%XBLOCK_VEGTYPE=>XBLOCK_VEGTYPE
  PACK_ISBA_MODEL(KFROM)%XBLOCK_TG=>XBLOCK_TG
  PACK_ISBA_MODEL(KFROM)%XBLOCK_SNOW=>XBLOCK_SNOW
  PACK_ISBA_MODEL(KFROM)%XBLOCK_ALB=>XBLOCK_ALB
  PACK_ISBA_MODEL(KFROM)%XBLOCK_2=>XBLOCK_2
  PACK_ISBA_MODEL(KFROM)%XBLOCK_BIOMASS=>XBLOCK_BIOMASS
  PACK_ISBA_MODEL(KFROM)%XBLOCK_SOILCARB=>XBLOCK_SOILCARB
  PACK_ISBA_MODEL(KFROM)%XBLOCK_LITTLEVS=>XBLOCK_LITTLEVS
  PACK_ISBA_MODEL(KFROM)%XBLOCK_LITTER=>XBLOCK_LITTER
  PACK_ISBA_MODEL(KFROM)%XBLOCK_0=>XBLOCK_0
  PACK_ISBA_MODEL(KFROM)%XBLOCK_00=>XBLOCK_00
  PACK_ISBA_MODEL(KFROM)%XBLOCK_000=>XBLOCK_000
  PACK_ISBA_MODEL(KFROM)%XBLOCK_01=>XBLOCK_01
  PACK_ISBA_MODEL(KFROM)%XP_VEGTYPE_PATCH=>XP_VEGTYPE_PATCH
  PACK_ISBA_MODEL(KFROM)%XP_SSO_SLOPE=>XP_SSO_SLOPE
  PACK_ISBA_MODEL(KFROM)%XP_LAT=>XP_LAT
  PACK_ISBA_MODEL(KFROM)%XP_LON=>XP_LON
  PACK_ISBA_MODEL(KFROM)%XP_AOSIP=>XP_AOSIP
  PACK_ISBA_MODEL(KFROM)%XP_AOSIM=>XP_AOSIM
  PACK_ISBA_MODEL(KFROM)%XP_AOSJP=>XP_AOSJP
  PACK_ISBA_MODEL(KFROM)%XP_AOSJM=>XP_AOSJM
  PACK_ISBA_MODEL(KFROM)%XP_HO2IP=>XP_HO2IP
  PACK_ISBA_MODEL(KFROM)%XP_HO2IM=>XP_HO2IM
  PACK_ISBA_MODEL(KFROM)%XP_HO2JP=>XP_HO2JP
  PACK_ISBA_MODEL(KFROM)%XP_HO2JM=>XP_HO2JM
  PACK_ISBA_MODEL(KFROM)%XP_Z0EFFIP=>XP_Z0EFFIP
  PACK_ISBA_MODEL(KFROM)%XP_Z0EFFIM=>XP_Z0EFFIM
  PACK_ISBA_MODEL(KFROM)%XP_Z0EFFJP=>XP_Z0EFFJP
  PACK_ISBA_MODEL(KFROM)%XP_Z0EFFJM=>XP_Z0EFFJM
  PACK_ISBA_MODEL(KFROM)%XP_Z0REL=>XP_Z0REL
  PACK_ISBA_MODEL(KFROM)%XP_CLAY=>XP_CLAY
  PACK_ISBA_MODEL(KFROM)%XP_SAND=>XP_SAND
  PACK_ISBA_MODEL(KFROM)%XP_ALBNIR_DRY=>XP_ALBNIR_DRY
  PACK_ISBA_MODEL(KFROM)%XP_ALBVIS_DRY=>XP_ALBVIS_DRY
  PACK_ISBA_MODEL(KFROM)%XP_ALBUV_DRY=>XP_ALBUV_DRY
  PACK_ISBA_MODEL(KFROM)%XP_ALBNIR_WET=>XP_ALBNIR_WET
  PACK_ISBA_MODEL(KFROM)%XP_ALBVIS_WET=>XP_ALBVIS_WET
  PACK_ISBA_MODEL(KFROM)%XP_ALBUV_WET=>XP_ALBUV_WET
  PACK_ISBA_MODEL(KFROM)%XP_ALBNIR_SOIL=>XP_ALBNIR_SOIL
  PACK_ISBA_MODEL(KFROM)%XP_ALBVIS_SOIL=>XP_ALBVIS_SOIL
  PACK_ISBA_MODEL(KFROM)%XP_ALBUV_SOIL=>XP_ALBUV_SOIL
  PACK_ISBA_MODEL(KFROM)%XP_Z0_O_Z0H=>XP_Z0_O_Z0H
  PACK_ISBA_MODEL(KFROM)%XP_ALBNIR=>XP_ALBNIR
  PACK_ISBA_MODEL(KFROM)%XP_ALBVIS=>XP_ALBVIS
  PACK_ISBA_MODEL(KFROM)%XP_ALBUV=>XP_ALBUV
  PACK_ISBA_MODEL(KFROM)%XP_EMIS=>XP_EMIS
  PACK_ISBA_MODEL(KFROM)%XP_Z0=>XP_Z0
  PACK_ISBA_MODEL(KFROM)%XP_ALBNIR_VEG=>XP_ALBNIR_VEG
  PACK_ISBA_MODEL(KFROM)%XP_ALBVIS_VEG=>XP_ALBVIS_VEG
  PACK_ISBA_MODEL(KFROM)%XP_ALBUV_VEG=>XP_ALBUV_VEG
  PACK_ISBA_MODEL(KFROM)%XP_VEG=>XP_VEG
  PACK_ISBA_MODEL(KFROM)%XP_WRMAX_CF=>XP_WRMAX_CF
  PACK_ISBA_MODEL(KFROM)%XP_RSMIN=>XP_RSMIN
  PACK_ISBA_MODEL(KFROM)%XP_GAMMA=>XP_GAMMA
  PACK_ISBA_MODEL(KFROM)%XP_CV=>XP_CV
  PACK_ISBA_MODEL(KFROM)%XP_RGL=>XP_RGL
  PACK_ISBA_MODEL(KFROM)%XP_ROOTFRAC=>XP_ROOTFRAC
  PACK_ISBA_MODEL(KFROM)%XP_BSLAI=>XP_BSLAI
  PACK_ISBA_MODEL(KFROM)%XP_LAIMIN=>XP_LAIMIN
  PACK_ISBA_MODEL(KFROM)%XP_SEFOLD=>XP_SEFOLD
  PACK_ISBA_MODEL(KFROM)%XP_H_TREE=>XP_H_TREE
  PACK_ISBA_MODEL(KFROM)%XP_ANF=>XP_ANF
  PACK_ISBA_MODEL(KFROM)%XP_ANMAX=>XP_ANMAX
  PACK_ISBA_MODEL(KFROM)%XP_FZERO=>XP_FZERO
  PACK_ISBA_MODEL(KFROM)%XP_EPSO=>XP_EPSO
  PACK_ISBA_MODEL(KFROM)%XP_GAMM=>XP_GAMM
  PACK_ISBA_MODEL(KFROM)%XP_QDGAMM=>XP_QDGAMM
  PACK_ISBA_MODEL(KFROM)%XP_GMES=>XP_GMES
  PACK_ISBA_MODEL(KFROM)%XP_RE25=>XP_RE25
  PACK_ISBA_MODEL(KFROM)%XP_QDGMES=>XP_QDGMES
  PACK_ISBA_MODEL(KFROM)%XP_T1GMES=>XP_T1GMES
  PACK_ISBA_MODEL(KFROM)%XP_T2GMES=>XP_T2GMES
  PACK_ISBA_MODEL(KFROM)%XP_AMAX=>XP_AMAX
  PACK_ISBA_MODEL(KFROM)%XP_QDAMAX=>XP_QDAMAX
  PACK_ISBA_MODEL(KFROM)%XP_T1AMAX=>XP_T1AMAX
  PACK_ISBA_MODEL(KFROM)%XP_T2AMAX=>XP_T2AMAX
  PACK_ISBA_MODEL(KFROM)%LP_STRESS=>LP_STRESS
  PACK_ISBA_MODEL(KFROM)%XP_F2I=>XP_F2I
  PACK_ISBA_MODEL(KFROM)%XP_GC=>XP_GC
  PACK_ISBA_MODEL(KFROM)%XP_AH=>XP_AH
  PACK_ISBA_MODEL(KFROM)%XP_BH=>XP_BH
  PACK_ISBA_MODEL(KFROM)%XP_TAU_WOOD=>XP_TAU_WOOD
  PACK_ISBA_MODEL(KFROM)%XP_DMAX=>XP_DMAX
  PACK_ISBA_MODEL(KFROM)%XP_CE_NITRO=>XP_CE_NITRO
  PACK_ISBA_MODEL(KFROM)%XP_CF_NITRO=>XP_CF_NITRO
  PACK_ISBA_MODEL(KFROM)%XP_CNA_NITRO=>XP_CNA_NITRO
  PACK_ISBA_MODEL(KFROM)%XP_BSLAI_NITRO=>XP_BSLAI_NITRO
  PACK_ISBA_MODEL(KFROM)%XP_RUNOFFB=>XP_RUNOFFB
  PACK_ISBA_MODEL(KFROM)%XP_WDRAIN=>XP_WDRAIN
  PACK_ISBA_MODEL(KFROM)%XP_TAUICE=>XP_TAUICE
  PACK_ISBA_MODEL(KFROM)%XP_GAMMAT=>XP_GAMMAT
  PACK_ISBA_MODEL(KFROM)%XP_DG=>XP_DG
  PACK_ISBA_MODEL(KFROM)%XP_DZG=>XP_DZG
  PACK_ISBA_MODEL(KFROM)%XP_DZDIF=>XP_DZDIF
  PACK_ISBA_MODEL(KFROM)%NK_WG_LAYER=>NK_WG_LAYER
  PACK_ISBA_MODEL(KFROM)%XP_RUNOFFD=>XP_RUNOFFD
  PACK_ISBA_MODEL(KFROM)%XP_SOILWGHT=>XP_SOILWGHT
  PACK_ISBA_MODEL(KFROM)%XP_C1SAT=>XP_C1SAT
  PACK_ISBA_MODEL(KFROM)%XP_C2REF=>XP_C2REF
  PACK_ISBA_MODEL(KFROM)%XP_C3=>XP_C3
  PACK_ISBA_MODEL(KFROM)%XP_C4B=>XP_C4B
  PACK_ISBA_MODEL(KFROM)%XP_C4REF=>XP_C4REF
  PACK_ISBA_MODEL(KFROM)%XP_ACOEF=>XP_ACOEF
  PACK_ISBA_MODEL(KFROM)%XP_PCOEF=>XP_PCOEF
  PACK_ISBA_MODEL(KFROM)%XP_WFC=>XP_WFC
  PACK_ISBA_MODEL(KFROM)%XP_WWILT=>XP_WWILT
  PACK_ISBA_MODEL(KFROM)%XP_WSAT=>XP_WSAT
  PACK_ISBA_MODEL(KFROM)%XP_BCOEF=>XP_BCOEF
  PACK_ISBA_MODEL(KFROM)%XP_CONDSAT=>XP_CONDSAT
  PACK_ISBA_MODEL(KFROM)%XP_MPOTSAT=>XP_MPOTSAT
  PACK_ISBA_MODEL(KFROM)%XP_CGSAT=>XP_CGSAT
  PACK_ISBA_MODEL(KFROM)%XP_HCAPSOIL=>XP_HCAPSOIL
  PACK_ISBA_MODEL(KFROM)%XP_CONDDRY=>XP_CONDDRY
  PACK_ISBA_MODEL(KFROM)%XP_CONDSLD=>XP_CONDSLD
  PACK_ISBA_MODEL(KFROM)%XP_TDEEP=>XP_TDEEP
  PACK_ISBA_MODEL(KFROM)%XP_SNOWSWE=>XP_SNOWSWE
  PACK_ISBA_MODEL(KFROM)%XP_SNOWHEAT=>XP_SNOWHEAT
  PACK_ISBA_MODEL(KFROM)%XP_SNOWRHO=>XP_SNOWRHO
  PACK_ISBA_MODEL(KFROM)%XP_SNOWGRAN1=>XP_SNOWGRAN1
  PACK_ISBA_MODEL(KFROM)%XP_SNOWGRAN2=>XP_SNOWGRAN2
  PACK_ISBA_MODEL(KFROM)%XP_SNOWHIST=>XP_SNOWHIST
  PACK_ISBA_MODEL(KFROM)%XP_SNOWAGE=>XP_SNOWAGE
  PACK_ISBA_MODEL(KFROM)%XP_SNOWALB=>XP_SNOWALB
  PACK_ISBA_MODEL(KFROM)%XP_SNOWALBVIS=>XP_SNOWALBVIS
  PACK_ISBA_MODEL(KFROM)%XP_SNOWALBNIR=>XP_SNOWALBNIR
  PACK_ISBA_MODEL(KFROM)%XP_SNOWALBFIR=>XP_SNOWALBFIR
  PACK_ISBA_MODEL(KFROM)%XP_SNOWEMIS=>XP_SNOWEMIS
  PACK_ISBA_MODEL(KFROM)%XP_ICE_STO=>XP_ICE_STO
  PACK_ISBA_MODEL(KFROM)%XP_WR=>XP_WR
  PACK_ISBA_MODEL(KFROM)%XP_TG=>XP_TG
  PACK_ISBA_MODEL(KFROM)%XP_WG=>XP_WG
  PACK_ISBA_MODEL(KFROM)%XP_WGI=>XP_WGI
  PACK_ISBA_MODEL(KFROM)%XP_RESA=>XP_RESA
  PACK_ISBA_MODEL(KFROM)%XP_WRV=>XP_WRV
  PACK_ISBA_MODEL(KFROM)%XP_WRVN=>XP_WRVN
  PACK_ISBA_MODEL(KFROM)%XP_TV=>XP_TV
  PACK_ISBA_MODEL(KFROM)%XP_TC=>XP_TC
  PACK_ISBA_MODEL(KFROM)%XP_QC=>XP_QC
  PACK_ISBA_MODEL(KFROM)%XP_ZF_TALLVEG=>XP_ZF_TALLVEG
  PACK_ISBA_MODEL(KFROM)%XP_RGLV=>XP_RGLV
  PACK_ISBA_MODEL(KFROM)%XP_GAMMAV=>XP_GAMMAV
  PACK_ISBA_MODEL(KFROM)%XP_RSMINV=>XP_RSMINV
  PACK_ISBA_MODEL(KFROM)%XP_ROOTFRACV=>XP_ROOTFRACV
  PACK_ISBA_MODEL(KFROM)%XP_WRMAX_CFV=>XP_WRMAX_CFV
  PACK_ISBA_MODEL(KFROM)%XP_LAIV=>XP_LAIV
  PACK_ISBA_MODEL(KFROM)%XP_Z0V=>XP_Z0V
  PACK_ISBA_MODEL(KFROM)%XP_H_VEG=>XP_H_VEG
  PACK_ISBA_MODEL(KFROM)%XP_GNDLITTER=>XP_GNDLITTER
  PACK_ISBA_MODEL(KFROM)%XP_Z0LITTER=>XP_Z0LITTER
  PACK_ISBA_MODEL(KFROM)%XP_FWTD=>XP_FWTD
  PACK_ISBA_MODEL(KFROM)%XP_WTD=>XP_WTD
  PACK_ISBA_MODEL(KFROM)%XP_LAI=>XP_LAI
  PACK_ISBA_MODEL(KFROM)%XP_AN=>XP_AN
  PACK_ISBA_MODEL(KFROM)%XP_ANDAY=>XP_ANDAY
  PACK_ISBA_MODEL(KFROM)%XP_ANFM=>XP_ANFM
  PACK_ISBA_MODEL(KFROM)%XP_LE=>XP_LE
  PACK_ISBA_MODEL(KFROM)%XP_LEI=>XP_LEI
  PACK_ISBA_MODEL(KFROM)%XP_FAPARC=>XP_FAPARC
  PACK_ISBA_MODEL(KFROM)%XP_FAPIRC=>XP_FAPIRC
  PACK_ISBA_MODEL(KFROM)%XP_LAI_EFFC=>XP_LAI_EFFC
  PACK_ISBA_MODEL(KFROM)%XP_MUS=>XP_MUS
  PACK_ISBA_MODEL(KFROM)%XP_RESP_BIOMASS=>XP_RESP_BIOMASS
  PACK_ISBA_MODEL(KFROM)%XP_BIOMASS=>XP_BIOMASS
  PACK_ISBA_MODEL(KFROM)%XP_INCREASE=>XP_INCREASE
  PACK_ISBA_MODEL(KFROM)%XP_LITTER=>XP_LITTER
  PACK_ISBA_MODEL(KFROM)%XP_SOILCARB=>XP_SOILCARB
  PACK_ISBA_MODEL(KFROM)%XP_LIGNIN_STRUC=>XP_LIGNIN_STRUC
  PACK_ISBA_MODEL(KFROM)%XP_TURNOVER=>XP_TURNOVER
  PACK_ISBA_MODEL(KFROM)%XP_LIRRIGATE=>XP_LIRRIGATE
  PACK_ISBA_MODEL(KFROM)%XP_LIRRIDAY=>XP_LIRRIDAY
  PACK_ISBA_MODEL(KFROM)%XP_THRESHOLD=>XP_THRESHOLD
  PACK_ISBA_MODEL(KFROM)%XP_WATSUP=>XP_WATSUP
  PACK_ISBA_MODEL(KFROM)%XP_IRRIG=>XP_IRRIG
  PACK_ISBA_MODEL(KFROM)%TP_SEED=>TP_SEED
  PACK_ISBA_MODEL(KFROM)%TP_REAP=>TP_REAP
  PACK_ISBA_MODEL(KFROM)%XP_D_ICE=>XP_D_ICE
  PACK_ISBA_MODEL(KFROM)%XP_KSAT_ICE=>XP_KSAT_ICE
  PACK_ISBA_MODEL(KFROM)%XP_FSAT=>XP_FSAT
  PACK_ISBA_MODEL(KFROM)%XP_MUF=>XP_MUF
  PACK_ISBA_MODEL(KFROM)%XP_TOPQS=>XP_TOPQS
  PACK_ISBA_MODEL(KFROM)%XP_PSN=>XP_PSN
  PACK_ISBA_MODEL(KFROM)%XP_PSNG=>XP_PSNG
  PACK_ISBA_MODEL(KFROM)%XP_PSNV=>XP_PSNV
  PACK_ISBA_MODEL(KFROM)%XP_PSNV_A=>XP_PSNV_A
  PACK_ISBA_MODEL(KFROM)%XP_DIR_ALB_WITH_SNOW=>XP_DIR_ALB_WITH_SNOW
  PACK_ISBA_MODEL(KFROM)%XP_SCA_ALB_WITH_SNOW=>XP_SCA_ALB_WITH_SNOW
  PACK_ISBA_MODEL(KFROM)%XP_ALBF=>XP_ALBF
  PACK_ISBA_MODEL(KFROM)%XP_EMISF=>XP_EMISF
  PACK_ISBA_MODEL(KFROM)%XP_FF=>XP_FF
  PACK_ISBA_MODEL(KFROM)%XP_FFG=>XP_FFG
  PACK_ISBA_MODEL(KFROM)%XP_FFV=>XP_FFV
  PACK_ISBA_MODEL(KFROM)%XP_FFROZEN=>XP_FFROZEN
  PACK_ISBA_MODEL(KFROM)%XP_FFLOOD=>XP_FFLOOD
  PACK_ISBA_MODEL(KFROM)%XP_PIFLOOD=>XP_PIFLOOD
  PACK_ISBA_MODEL(KFROM)%XP_CPS=>XP_CPS
  PACK_ISBA_MODEL(KFROM)%XP_LVTT=>XP_LVTT
  PACK_ISBA_MODEL(KFROM)%XP_LSTT=>XP_LSTT
ENDIF
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_PACK_ISBA_N:PACK_ISBA_GOTO_MODEL',0,ZHOOK_HANDLE)
NSIZE_LSIMPLE=>PACK_ISBA_MODEL(KTO)%NSIZE_LSIMPLE
NSIZE_L0=>PACK_ISBA_MODEL(KTO)%NSIZE_L0
NSIZE_NSIMPLE=>PACK_ISBA_MODEL(KTO)%NSIZE_NSIMPLE
NSIZE_N0=>PACK_ISBA_MODEL(KTO)%NSIZE_N0
NSIZE_TSIMPLE=>PACK_ISBA_MODEL(KTO)%NSIZE_TSIMPLE
NSIZE_T0=>PACK_ISBA_MODEL(KTO)%NSIZE_T0
NSIZE_SIMPLE=>PACK_ISBA_MODEL(KTO)%NSIZE_SIMPLE
NSIZE_GROUND=>PACK_ISBA_MODEL(KTO)%NSIZE_GROUND
NSIZE_VEGTYPE=>PACK_ISBA_MODEL(KTO)%NSIZE_VEGTYPE
NSIZE_TG=>PACK_ISBA_MODEL(KTO)%NSIZE_TG
NSIZE_SNOW=>PACK_ISBA_MODEL(KTO)%NSIZE_SNOW
NSIZE_ALB=>PACK_ISBA_MODEL(KTO)%NSIZE_ALB
NSIZE_2=>PACK_ISBA_MODEL(KTO)%NSIZE_2
NSIZE_BIOMASS=>PACK_ISBA_MODEL(KTO)%NSIZE_BIOMASS
NSIZE_SOILCARB=>PACK_ISBA_MODEL(KTO)%NSIZE_SOILCARB
NSIZE_LITTLEVS=>PACK_ISBA_MODEL(KTO)%NSIZE_LITTLEVS
NSIZE_LITTER=>PACK_ISBA_MODEL(KTO)%NSIZE_LITTER
NSIZE_0=>PACK_ISBA_MODEL(KTO)%NSIZE_0
NSIZE_00=>PACK_ISBA_MODEL(KTO)%NSIZE_00
NSIZE_000=>PACK_ISBA_MODEL(KTO)%NSIZE_000
NSIZE_01=>PACK_ISBA_MODEL(KTO)%NSIZE_01
LBLOCK_SIMPLE=>PACK_ISBA_MODEL(KTO)%LBLOCK_SIMPLE
LBLOCK_0=>PACK_ISBA_MODEL(KTO)%LBLOCK_0
NBLOCK_SIMPLE=>PACK_ISBA_MODEL(KTO)%NBLOCK_SIMPLE
NBLOCK_0=>PACK_ISBA_MODEL(KTO)%NBLOCK_0
TBLOCK_SIMPLE=>PACK_ISBA_MODEL(KTO)%TBLOCK_SIMPLE
TBLOCK_0=>PACK_ISBA_MODEL(KTO)%TBLOCK_0
XBLOCK_SIMPLE=>PACK_ISBA_MODEL(KTO)%XBLOCK_SIMPLE
XBLOCK_GROUND=>PACK_ISBA_MODEL(KTO)%XBLOCK_GROUND
XBLOCK_VEGTYPE=>PACK_ISBA_MODEL(KTO)%XBLOCK_VEGTYPE
XBLOCK_TG=>PACK_ISBA_MODEL(KTO)%XBLOCK_TG
XBLOCK_SNOW=>PACK_ISBA_MODEL(KTO)%XBLOCK_SNOW
XBLOCK_ALB=>PACK_ISBA_MODEL(KTO)%XBLOCK_ALB
XBLOCK_2=>PACK_ISBA_MODEL(KTO)%XBLOCK_2
XBLOCK_BIOMASS=>PACK_ISBA_MODEL(KTO)%XBLOCK_BIOMASS
XBLOCK_SOILCARB=>PACK_ISBA_MODEL(KTO)%XBLOCK_SOILCARB
XBLOCK_LITTLEVS=>PACK_ISBA_MODEL(KTO)%XBLOCK_LITTLEVS
XBLOCK_LITTER=>PACK_ISBA_MODEL(KTO)%XBLOCK_LITTER
XBLOCK_0=>PACK_ISBA_MODEL(KTO)%XBLOCK_0
XBLOCK_00=>PACK_ISBA_MODEL(KTO)%XBLOCK_00
XBLOCK_000=>PACK_ISBA_MODEL(KTO)%XBLOCK_000
XBLOCK_01=>PACK_ISBA_MODEL(KTO)%XBLOCK_01
XP_VEGTYPE_PATCH=>PACK_ISBA_MODEL(KTO)%XP_VEGTYPE_PATCH
XP_SSO_SLOPE=>PACK_ISBA_MODEL(KTO)%XP_SSO_SLOPE
XP_LAT=>PACK_ISBA_MODEL(KTO)%XP_LAT
XP_LON=>PACK_ISBA_MODEL(KTO)%XP_LON
XP_AOSIP=>PACK_ISBA_MODEL(KTO)%XP_AOSIP
XP_AOSIM=>PACK_ISBA_MODEL(KTO)%XP_AOSIM
XP_AOSJP=>PACK_ISBA_MODEL(KTO)%XP_AOSJP
XP_AOSJM=>PACK_ISBA_MODEL(KTO)%XP_AOSJM
XP_HO2IP=>PACK_ISBA_MODEL(KTO)%XP_HO2IP
XP_HO2IM=>PACK_ISBA_MODEL(KTO)%XP_HO2IM
XP_HO2JP=>PACK_ISBA_MODEL(KTO)%XP_HO2JP
XP_HO2JM=>PACK_ISBA_MODEL(KTO)%XP_HO2JM
XP_Z0EFFIP=>PACK_ISBA_MODEL(KTO)%XP_Z0EFFIP
XP_Z0EFFIM=>PACK_ISBA_MODEL(KTO)%XP_Z0EFFIM
XP_Z0EFFJP=>PACK_ISBA_MODEL(KTO)%XP_Z0EFFJP
XP_Z0EFFJM=>PACK_ISBA_MODEL(KTO)%XP_Z0EFFJM
XP_Z0REL=>PACK_ISBA_MODEL(KTO)%XP_Z0REL
XP_CLAY=>PACK_ISBA_MODEL(KTO)%XP_CLAY
XP_SAND=>PACK_ISBA_MODEL(KTO)%XP_SAND
XP_ALBNIR_DRY=>PACK_ISBA_MODEL(KTO)%XP_ALBNIR_DRY
XP_ALBVIS_DRY=>PACK_ISBA_MODEL(KTO)%XP_ALBVIS_DRY
XP_ALBUV_DRY=>PACK_ISBA_MODEL(KTO)%XP_ALBUV_DRY
XP_ALBNIR_WET=>PACK_ISBA_MODEL(KTO)%XP_ALBNIR_WET
XP_ALBVIS_WET=>PACK_ISBA_MODEL(KTO)%XP_ALBVIS_WET
XP_ALBUV_WET=>PACK_ISBA_MODEL(KTO)%XP_ALBUV_WET
XP_ALBNIR_SOIL=>PACK_ISBA_MODEL(KTO)%XP_ALBNIR_SOIL
XP_ALBVIS_SOIL=>PACK_ISBA_MODEL(KTO)%XP_ALBVIS_SOIL
XP_ALBUV_SOIL=>PACK_ISBA_MODEL(KTO)%XP_ALBUV_SOIL
XP_Z0_O_Z0H=>PACK_ISBA_MODEL(KTO)%XP_Z0_O_Z0H
XP_ALBNIR=>PACK_ISBA_MODEL(KTO)%XP_ALBNIR
XP_ALBVIS=>PACK_ISBA_MODEL(KTO)%XP_ALBVIS
XP_ALBUV=>PACK_ISBA_MODEL(KTO)%XP_ALBUV
XP_EMIS=>PACK_ISBA_MODEL(KTO)%XP_EMIS
XP_Z0=>PACK_ISBA_MODEL(KTO)%XP_Z0
XP_ALBNIR_VEG=>PACK_ISBA_MODEL(KTO)%XP_ALBNIR_VEG
XP_ALBVIS_VEG=>PACK_ISBA_MODEL(KTO)%XP_ALBVIS_VEG
XP_ALBUV_VEG=>PACK_ISBA_MODEL(KTO)%XP_ALBUV_VEG
XP_VEG=>PACK_ISBA_MODEL(KTO)%XP_VEG
XP_WRMAX_CF=>PACK_ISBA_MODEL(KTO)%XP_WRMAX_CF
XP_RSMIN=>PACK_ISBA_MODEL(KTO)%XP_RSMIN
XP_GAMMA=>PACK_ISBA_MODEL(KTO)%XP_GAMMA
XP_CV=>PACK_ISBA_MODEL(KTO)%XP_CV
XP_RGL=>PACK_ISBA_MODEL(KTO)%XP_RGL
XP_ROOTFRAC=>PACK_ISBA_MODEL(KTO)%XP_ROOTFRAC
XP_BSLAI=>PACK_ISBA_MODEL(KTO)%XP_BSLAI
XP_LAIMIN=>PACK_ISBA_MODEL(KTO)%XP_LAIMIN
XP_SEFOLD=>PACK_ISBA_MODEL(KTO)%XP_SEFOLD
XP_H_TREE=>PACK_ISBA_MODEL(KTO)%XP_H_TREE
XP_ANF=>PACK_ISBA_MODEL(KTO)%XP_ANF
XP_ANMAX=>PACK_ISBA_MODEL(KTO)%XP_ANMAX
XP_FZERO=>PACK_ISBA_MODEL(KTO)%XP_FZERO
XP_EPSO=>PACK_ISBA_MODEL(KTO)%XP_EPSO
XP_GAMM=>PACK_ISBA_MODEL(KTO)%XP_GAMM
XP_QDGAMM=>PACK_ISBA_MODEL(KTO)%XP_QDGAMM
XP_GMES=>PACK_ISBA_MODEL(KTO)%XP_GMES
XP_RE25=>PACK_ISBA_MODEL(KTO)%XP_RE25
XP_QDGMES=>PACK_ISBA_MODEL(KTO)%XP_QDGMES
XP_T1GMES=>PACK_ISBA_MODEL(KTO)%XP_T1GMES
XP_T2GMES=>PACK_ISBA_MODEL(KTO)%XP_T2GMES
XP_AMAX=>PACK_ISBA_MODEL(KTO)%XP_AMAX
XP_QDAMAX=>PACK_ISBA_MODEL(KTO)%XP_QDAMAX
XP_T1AMAX=>PACK_ISBA_MODEL(KTO)%XP_T1AMAX
XP_T2AMAX=>PACK_ISBA_MODEL(KTO)%XP_T2AMAX
LP_STRESS=>PACK_ISBA_MODEL(KTO)%LP_STRESS
XP_F2I=>PACK_ISBA_MODEL(KTO)%XP_F2I
XP_GC=>PACK_ISBA_MODEL(KTO)%XP_GC
XP_AH=>PACK_ISBA_MODEL(KTO)%XP_AH
XP_BH=>PACK_ISBA_MODEL(KTO)%XP_BH
XP_TAU_WOOD=>PACK_ISBA_MODEL(KTO)%XP_TAU_WOOD
XP_DMAX=>PACK_ISBA_MODEL(KTO)%XP_DMAX
XP_CE_NITRO=>PACK_ISBA_MODEL(KTO)%XP_CE_NITRO
XP_CF_NITRO=>PACK_ISBA_MODEL(KTO)%XP_CF_NITRO
XP_CNA_NITRO=>PACK_ISBA_MODEL(KTO)%XP_CNA_NITRO
XP_BSLAI_NITRO=>PACK_ISBA_MODEL(KTO)%XP_BSLAI_NITRO
XP_RUNOFFB=>PACK_ISBA_MODEL(KTO)%XP_RUNOFFB
XP_WDRAIN=>PACK_ISBA_MODEL(KTO)%XP_WDRAIN
XP_TAUICE=>PACK_ISBA_MODEL(KTO)%XP_TAUICE
XP_GAMMAT=>PACK_ISBA_MODEL(KTO)%XP_GAMMAT
XP_DG=>PACK_ISBA_MODEL(KTO)%XP_DG
XP_DZG=>PACK_ISBA_MODEL(KTO)%XP_DZG
XP_DZDIF=>PACK_ISBA_MODEL(KTO)%XP_DZDIF
NK_WG_LAYER=>PACK_ISBA_MODEL(KTO)%NK_WG_LAYER
XP_RUNOFFD=>PACK_ISBA_MODEL(KTO)%XP_RUNOFFD
XP_SOILWGHT=>PACK_ISBA_MODEL(KTO)%XP_SOILWGHT
XP_C1SAT=>PACK_ISBA_MODEL(KTO)%XP_C1SAT
XP_C2REF=>PACK_ISBA_MODEL(KTO)%XP_C2REF
XP_C3=>PACK_ISBA_MODEL(KTO)%XP_C3
XP_C4B=>PACK_ISBA_MODEL(KTO)%XP_C4B
XP_C4REF=>PACK_ISBA_MODEL(KTO)%XP_C4REF
XP_ACOEF=>PACK_ISBA_MODEL(KTO)%XP_ACOEF
XP_PCOEF=>PACK_ISBA_MODEL(KTO)%XP_PCOEF
XP_WFC=>PACK_ISBA_MODEL(KTO)%XP_WFC
XP_WWILT=>PACK_ISBA_MODEL(KTO)%XP_WWILT
XP_WSAT=>PACK_ISBA_MODEL(KTO)%XP_WSAT
XP_BCOEF=>PACK_ISBA_MODEL(KTO)%XP_BCOEF
XP_CONDSAT=>PACK_ISBA_MODEL(KTO)%XP_CONDSAT
XP_MPOTSAT=>PACK_ISBA_MODEL(KTO)%XP_MPOTSAT
XP_CGSAT=>PACK_ISBA_MODEL(KTO)%XP_CGSAT
XP_HCAPSOIL=>PACK_ISBA_MODEL(KTO)%XP_HCAPSOIL
XP_CONDDRY=>PACK_ISBA_MODEL(KTO)%XP_CONDDRY
XP_CONDSLD=>PACK_ISBA_MODEL(KTO)%XP_CONDSLD
XP_TDEEP=>PACK_ISBA_MODEL(KTO)%XP_TDEEP
XP_SNOWSWE=>PACK_ISBA_MODEL(KTO)%XP_SNOWSWE
XP_SNOWHEAT=>PACK_ISBA_MODEL(KTO)%XP_SNOWHEAT
XP_SNOWRHO=>PACK_ISBA_MODEL(KTO)%XP_SNOWRHO
XP_SNOWGRAN1=>PACK_ISBA_MODEL(KTO)%XP_SNOWGRAN1
XP_SNOWGRAN2=>PACK_ISBA_MODEL(KTO)%XP_SNOWGRAN2
XP_SNOWHIST=>PACK_ISBA_MODEL(KTO)%XP_SNOWHIST
XP_SNOWAGE=>PACK_ISBA_MODEL(KTO)%XP_SNOWAGE
XP_SNOWALB=>PACK_ISBA_MODEL(KTO)%XP_SNOWALB
XP_SNOWALBVIS=>PACK_ISBA_MODEL(KTO)%XP_SNOWALBVIS
XP_SNOWALBNIR=>PACK_ISBA_MODEL(KTO)%XP_SNOWALBNIR
XP_SNOWALBFIR=>PACK_ISBA_MODEL(KTO)%XP_SNOWALBFIR
XP_SNOWEMIS=>PACK_ISBA_MODEL(KTO)%XP_SNOWEMIS
XP_ICE_STO=>PACK_ISBA_MODEL(KTO)%XP_ICE_STO
XP_WR=>PACK_ISBA_MODEL(KTO)%XP_WR
XP_TG=>PACK_ISBA_MODEL(KTO)%XP_TG
XP_WG=>PACK_ISBA_MODEL(KTO)%XP_WG
XP_WGI=>PACK_ISBA_MODEL(KTO)%XP_WGI
XP_RESA=>PACK_ISBA_MODEL(KTO)%XP_RESA
XP_WRV=>PACK_ISBA_MODEL(KTO)%XP_WRV
XP_WRVN=>PACK_ISBA_MODEL(KTO)%XP_WRVN
XP_TV=>PACK_ISBA_MODEL(KTO)%XP_TV
XP_TC=>PACK_ISBA_MODEL(KTO)%XP_TC
XP_QC=>PACK_ISBA_MODEL(KTO)%XP_QC
XP_ZF_TALLVEG=>PACK_ISBA_MODEL(KTO)%XP_ZF_TALLVEG
XP_RGLV=>PACK_ISBA_MODEL(KTO)%XP_RGLV
XP_GAMMAV=>PACK_ISBA_MODEL(KTO)%XP_GAMMAV
XP_RSMINV=>PACK_ISBA_MODEL(KTO)%XP_RSMINV
XP_ROOTFRACV=>PACK_ISBA_MODEL(KTO)%XP_ROOTFRACV
XP_WRMAX_CFV=>PACK_ISBA_MODEL(KTO)%XP_WRMAX_CFV
XP_LAIV=>PACK_ISBA_MODEL(KTO)%XP_LAIV
XP_Z0V=>PACK_ISBA_MODEL(KTO)%XP_Z0V
XP_H_VEG=>PACK_ISBA_MODEL(KTO)%XP_H_VEG
XP_GNDLITTER=>PACK_ISBA_MODEL(KTO)%XP_GNDLITTER
XP_Z0LITTER=>PACK_ISBA_MODEL(KTO)%XP_Z0LITTER
XP_FWTD=>PACK_ISBA_MODEL(KTO)%XP_FWTD
XP_WTD=>PACK_ISBA_MODEL(KTO)%XP_WTD
XP_LAI=>PACK_ISBA_MODEL(KTO)%XP_LAI
XP_AN=>PACK_ISBA_MODEL(KTO)%XP_AN
XP_ANDAY=>PACK_ISBA_MODEL(KTO)%XP_ANDAY
XP_ANFM=>PACK_ISBA_MODEL(KTO)%XP_ANFM
XP_LE=>PACK_ISBA_MODEL(KTO)%XP_LE
XP_LEI=>PACK_ISBA_MODEL(KTO)%XP_LEI
XP_FAPARC=>PACK_ISBA_MODEL(KTO)%XP_FAPARC
XP_FAPIRC=>PACK_ISBA_MODEL(KTO)%XP_FAPIRC
XP_LAI_EFFC=>PACK_ISBA_MODEL(KTO)%XP_LAI_EFFC
XP_MUS=>PACK_ISBA_MODEL(KTO)%XP_MUS
XP_RESP_BIOMASS=>PACK_ISBA_MODEL(KTO)%XP_RESP_BIOMASS
XP_BIOMASS=>PACK_ISBA_MODEL(KTO)%XP_BIOMASS
XP_INCREASE=>PACK_ISBA_MODEL(KTO)%XP_INCREASE
XP_LITTER=>PACK_ISBA_MODEL(KTO)%XP_LITTER
XP_SOILCARB=>PACK_ISBA_MODEL(KTO)%XP_SOILCARB
XP_LIGNIN_STRUC=>PACK_ISBA_MODEL(KTO)%XP_LIGNIN_STRUC
XP_TURNOVER=>PACK_ISBA_MODEL(KTO)%XP_TURNOVER
XP_LIRRIGATE=>PACK_ISBA_MODEL(KTO)%XP_LIRRIGATE
XP_LIRRIDAY=>PACK_ISBA_MODEL(KTO)%XP_LIRRIDAY
XP_THRESHOLD=>PACK_ISBA_MODEL(KTO)%XP_THRESHOLD
XP_WATSUP=>PACK_ISBA_MODEL(KTO)%XP_WATSUP
XP_IRRIG=>PACK_ISBA_MODEL(KTO)%XP_IRRIG
TP_SEED=>PACK_ISBA_MODEL(KTO)%TP_SEED
TP_REAP=>PACK_ISBA_MODEL(KTO)%TP_REAP
XP_D_ICE=>PACK_ISBA_MODEL(KTO)%XP_D_ICE
XP_KSAT_ICE=>PACK_ISBA_MODEL(KTO)%XP_KSAT_ICE
XP_FSAT=>PACK_ISBA_MODEL(KTO)%XP_FSAT
XP_MUF=>PACK_ISBA_MODEL(KTO)%XP_MUF
XP_TOPQS=>PACK_ISBA_MODEL(KTO)%XP_TOPQS
XP_PSN=>PACK_ISBA_MODEL(KTO)%XP_PSN
XP_PSNG=>PACK_ISBA_MODEL(KTO)%XP_PSNG
XP_PSNV=>PACK_ISBA_MODEL(KTO)%XP_PSNV
XP_PSNV_A=>PACK_ISBA_MODEL(KTO)%XP_PSNV_A
XP_DIR_ALB_WITH_SNOW=>PACK_ISBA_MODEL(KTO)%XP_DIR_ALB_WITH_SNOW
XP_SCA_ALB_WITH_SNOW=>PACK_ISBA_MODEL(KTO)%XP_SCA_ALB_WITH_SNOW
XP_ALBF=>PACK_ISBA_MODEL(KTO)%XP_ALBF
XP_EMISF=>PACK_ISBA_MODEL(KTO)%XP_EMISF
XP_FF=>PACK_ISBA_MODEL(KTO)%XP_FF
XP_FFG=>PACK_ISBA_MODEL(KTO)%XP_FFG
XP_FFV=>PACK_ISBA_MODEL(KTO)%XP_FFV
XP_FFROZEN=>PACK_ISBA_MODEL(KTO)%XP_FFROZEN
XP_FFLOOD=>PACK_ISBA_MODEL(KTO)%XP_FFLOOD
XP_PIFLOOD=>PACK_ISBA_MODEL(KTO)%XP_PIFLOOD
XP_CPS=>PACK_ISBA_MODEL(KTO)%XP_CPS
XP_LVTT=>PACK_ISBA_MODEL(KTO)%XP_LVTT
XP_LSTT=>PACK_ISBA_MODEL(KTO)%XP_LSTT
IF (LHOOK) CALL DR_HOOK('MODD_PACK_ISBA_N:PACK_ISBA_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE PACK_ISBA_GOTO_MODEL

SUBROUTINE PACK_ISBA_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_PACK_ISBA_N:PACK_ISBA_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(PACK_ISBA_MODEL(KMODEL))
DO J=1,KMODEL
  NULLIFY(PACK_ISBA_MODEL(J)%LBLOCK_SIMPLE)
  NULLIFY(PACK_ISBA_MODEL(J)%LBLOCK_0)
  NULLIFY(PACK_ISBA_MODEL(J)%NBLOCK_SIMPLE)
  NULLIFY(PACK_ISBA_MODEL(J)%NBLOCK_0)
  NULLIFY(PACK_ISBA_MODEL(J)%TBLOCK_SIMPLE)
  NULLIFY(PACK_ISBA_MODEL(J)%TBLOCK_0)
  NULLIFY(PACK_ISBA_MODEL(J)%XBLOCK_SIMPLE)
  NULLIFY(PACK_ISBA_MODEL(J)%XBLOCK_GROUND)
  NULLIFY(PACK_ISBA_MODEL(J)%XBLOCK_VEGTYPE)
  NULLIFY(PACK_ISBA_MODEL(J)%XBLOCK_TG)
  NULLIFY(PACK_ISBA_MODEL(J)%XBLOCK_SNOW)
  NULLIFY(PACK_ISBA_MODEL(J)%XBLOCK_ALB)
  NULLIFY(PACK_ISBA_MODEL(J)%XBLOCK_2)
  NULLIFY(PACK_ISBA_MODEL(J)%XBLOCK_BIOMASS)
  NULLIFY(PACK_ISBA_MODEL(J)%XBLOCK_SOILCARB)
  NULLIFY(PACK_ISBA_MODEL(J)%XBLOCK_LITTLEVS)
  NULLIFY(PACK_ISBA_MODEL(J)%XBLOCK_LITTER)
  NULLIFY(PACK_ISBA_MODEL(J)%XBLOCK_0)
  NULLIFY(PACK_ISBA_MODEL(J)%XBLOCK_00)
  NULLIFY(PACK_ISBA_MODEL(J)%XBLOCK_000)
  NULLIFY(PACK_ISBA_MODEL(J)%XBLOCK_01)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_VEGTYPE_PATCH)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_SSO_SLOPE)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_LAT)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_LON)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_AOSIP)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_AOSIM)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_AOSJP)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_AOSJM)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_HO2IP)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_HO2IM)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_HO2JP)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_HO2JM)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_Z0EFFIP)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_Z0EFFIM)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_Z0EFFJP)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_Z0EFFJM)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_Z0REL)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_CLAY)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_SAND)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_ALBNIR_DRY)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_ALBVIS_DRY)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_ALBUV_DRY)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_ALBNIR_WET)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_ALBVIS_WET)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_ALBUV_WET)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_ALBNIR_SOIL)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_ALBVIS_SOIL)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_ALBUV_SOIL)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_Z0_O_Z0H)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_ALBNIR)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_ALBVIS)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_ALBUV)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_EMIS)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_Z0)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_ALBNIR_VEG)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_ALBVIS_VEG)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_ALBUV_VEG)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_VEG)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_WRMAX_CF)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_RSMIN)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_GAMMA)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_CV)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_RGL)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_ROOTFRAC)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_BSLAI)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_LAIMIN)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_SEFOLD)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_H_TREE)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_ANF)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_ANMAX)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_FZERO)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_EPSO)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_GAMM)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_QDGAMM)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_GMES)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_RE25)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_QDGMES)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_T1GMES)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_T2GMES)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_AMAX)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_QDAMAX)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_T1AMAX)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_T2AMAX)
  NULLIFY(PACK_ISBA_MODEL(J)%LP_STRESS)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_F2I)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_GC)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_AH)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_BH)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_TAU_WOOD)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_DMAX)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_CE_NITRO)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_CF_NITRO)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_CNA_NITRO)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_BSLAI_NITRO)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_RUNOFFB)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_WDRAIN)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_TAUICE)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_GAMMAT)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_DG)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_DZG)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_DZDIF)
  NULLIFY(PACK_ISBA_MODEL(J)%NK_WG_LAYER)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_RUNOFFD)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_SOILWGHT)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_C1SAT)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_C2REF)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_C3)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_C4B)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_C4REF)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_ACOEF)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_PCOEF)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_WFC)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_WWILT)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_WSAT)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_BCOEF)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_CONDSAT)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_MPOTSAT)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_CGSAT)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_HCAPSOIL)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_CONDDRY)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_CONDSLD)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_TDEEP)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_SNOWSWE)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_SNOWHEAT)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_SNOWRHO)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_SNOWGRAN1)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_SNOWGRAN2)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_SNOWHIST)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_SNOWAGE)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_SNOWALB)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_SNOWALBVIS)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_SNOWALBNIR)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_SNOWALBFIR)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_SNOWEMIS)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_ICE_STO)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_WR)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_TG)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_WG)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_WGI)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_RESA)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_WRV)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_WRVN)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_TV)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_TC)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_QC)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_ZF_TALLVEG)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_RGLV)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_GAMMAV)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_RSMINV)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_ROOTFRACV)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_WRMAX_CFV)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_LAIV)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_Z0V)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_H_VEG)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_GNDLITTER)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_Z0LITTER)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_FWTD)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_WTD)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_LAI)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_AN)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_ANDAY)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_ANFM)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_LE)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_LEI)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_FAPARC)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_FAPIRC)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_LAI_EFFC)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_MUS)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_RESP_BIOMASS)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_BIOMASS)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_INCREASE)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_LITTER)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_SOILCARB)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_LIGNIN_STRUC)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_TURNOVER)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_LIRRIGATE)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_LIRRIDAY)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_THRESHOLD)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_WATSUP)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_IRRIG)
  NULLIFY(PACK_ISBA_MODEL(J)%TP_SEED)
  NULLIFY(PACK_ISBA_MODEL(J)%TP_REAP)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_D_ICE)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_KSAT_ICE)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_FSAT)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_MUF)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_TOPQS)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_PSN)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_PSNG)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_PSNV)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_PSNV_A)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_DIR_ALB_WITH_SNOW)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_SCA_ALB_WITH_SNOW)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_ALBF)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_EMISF)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_FF)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_FFG)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_FFV)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_FFROZEN)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_FFLOOD)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_PIFLOOD)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_CPS)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_LVTT)
  NULLIFY(PACK_ISBA_MODEL(J)%XP_LSTT)
ENDDO
IF (LHOOK) CALL DR_HOOK("MODD_PACK_ISBA_N:PACK_ISBA_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE PACK_ISBA_ALLOC

SUBROUTINE PACK_ISBA_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_PACK_ISBA_N:PACK_ISBA_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(PACK_ISBA_MODEL)) DEALLOCATE(PACK_ISBA_MODEL)
IF (LHOOK) CALL DR_HOOK("MODD_PACK_ISBA_N:PACK_ISBA_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE PACK_ISBA_DEALLO

END MODULE MODD_PACK_ISBA
