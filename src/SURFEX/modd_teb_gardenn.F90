!##################
MODULE MODD_TEB_GARDEN_n
!##################
!
!!****  *MODD_TEB_GARDEN - declaration of packed surface parameters for ISBA scheme
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
!!	A. Boone   *Meteo France*
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
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TYPE_DATE_SURF
USE MODD_TYPE_SNOW
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE TEB_GARDEN_t
!-------------------------------------------------------------------------------
!
! ISBA Scheme Options:
!
  CHARACTER(LEN=4)               :: CROUGH   ! type of roughness length
                                           ! 'Z01D'
                                           ! 'Z04D'
  CHARACTER(LEN=3)               :: CISBA    ! type of ISBA version:
!                                          ! '2-L' (default)
!                                          ! '3-L'
!                                          ! 'DIF'
!
  CHARACTER(LEN=4)               :: CPEDOTF! NOTE: Only used when HISBA = DIF
!                                          ! 'CH78' = Clapp and Hornberger 1978 for BC (Default)
!                                          ! 'CO84' = Cosby et al. 1988 for BC
!
  CHARACTER(LEN=3)               :: CPHOTO   ! type of photosynthesis
!                                          ! 'NON'
!                                          ! 'AGS'
!                                          ! 'LAI'
!                                          ! 'LST'
!                                          ! 'AST'
!                                          ! 'NIT'
!                                          ! 'NCB'
  LOGICAL                        :: LTR_ML
  CHARACTER(LEN=4)               :: CALBEDO  ! albedo type
!                                          ! 'DRY ' 
!                                          ! 'EVOL' 
!                                          ! 'WET ' 
!                                          ! 'USER' 
  CHARACTER(LEN=4)               :: CSCOND   ! Thermal conductivity
!                                          ! 'DEF ' = DEFault: NP89 implicit method
!                                          ! 'PL98' = Peters-Lidard et al. 1998 used
!                                          ! for explicit computation of CG
  CHARACTER(LEN=4)               :: CC1DRY   ! C1 formulation for dry soils
!                                          ! 'DEF ' = DEFault: Giard-Bazile formulation
!                                          ! 'GB93' = Giordani 1993, Braud 1993 
!                                          !discontinuous at WILT
  CHARACTER(LEN=3)               :: CSOILFRZ ! soil freezing-physics option
!                                          ! 'DEF' = Default (Boone et al. 2000; 
!                                          !        Giard and Bazile 2000)
!                                          ! 'LWT' = Phase changes as above,
!                                          !         but relation between unfrozen 
!                                          !         water and temperature considered
!                            NOTE that when using the YISBA='DIF' multi-layer soil option,
!                            the 'LWT' method is used. It is only an option
!                            when using the force-restore soil method ('2-L' or '3-L')
!
  CHARACTER(LEN=4)               :: CDIFSFCOND ! Mulch effects
!                                          ! 'MLCH' = include the insulating effect of
!                                          ! leaf litter/mulch on the surf. thermal cond.
!                                          ! 'DEF ' = no mulch effect
!                           NOTE: Only used when YISBA = DIF
!
  CHARACTER(LEN=3)               :: CSNOWRES ! Turbulent exchanges over snow
!	                                   ! 'DEF' = Default: Louis (ISBA)
!       	                           ! 'RIL' = Maximum Richardson number limit
!                                          !         for stable conditions ISBA-SNOW3L
!                                          !         turbulent exchange option
!                                           
  CHARACTER(LEN=3)               :: CRESPSL  ! Soil respiration
!                                          ! 'DEF' = Default: Norman (1992)
!                                          ! 'PRM' = New Parameterization
!                                          ! 'CNT' = CENTURY model (Gibelin 2007)
!                                           
  CHARACTER(LEN=3)               :: CCPSURF! specific heat at surface
!                                          ! 'DRY' = default value (dry Cp)
!                                          ! 'HUM' = Cp as a fct of specific humidity
!
!-------------------------------------------------------------------------------
!
  LOGICAL                        :: LCANOPY ! T: SBL scheme within the canopy
!                                           ! F: no atmospheric layers below forcing level
  LOGICAL                        :: LCANOPY_DRAG ! T: drag activated in SBL scheme within the canopy
!                                                ! F: no drag activated in SBL atmospheric layers
!-------------------------------------------------------------------------------
!
! type of initialization of vegetation: from cover types (ecoclimap) or parameters prescribed
!
  LOGICAL                        :: LPAR_GARDEN ! T: parameters computed from ecoclimap
!                                               ! F: they are read in the file
!
!-------------------------------------------------------------------------------
!
! Mask and number of grid elements containing patches/tiles:
!
  INTEGER, POINTER, DIMENSION(:)   :: NSIZE_NATURE_P ! number of sub-patchs/tiles              (-)
  INTEGER, POINTER, DIMENSION(:,:) :: NR_NATURE_P    ! patch/tile mask                         (-)
  REAL, POINTER, DIMENSION(:,:)    :: XPATCH         ! fraction of each tile/patch             (-)
  REAL, POINTER, DIMENSION(:,:)    :: XVEGTYPE       ! fraction of each vegetation type for
!                                                      ! each grid mesh                          (-)
  REAL, POINTER, DIMENSION(:,:,:)  :: XVEGTYPE_PATCH ! fraction of each vegetation type for
!                                                      ! each vegetation unit/patch              (-)
  INTEGER                          :: NPATCH           ! maximum number of sub-tiles (patches)
!                                                      ! used at any grid point within a 
!                                                      ! natural surface fraction
  INTEGER                          :: NGROUND_LAYER    ! number of ground layers
!
  REAL, POINTER, DIMENSION(:)      :: XSOILGRID      ! Soil layer grid as reference for DIF
  INTEGER                              :: NNBIOMASS    ! number of biomass pools
  INTEGER                              :: NNLITTER     ! number of litter pools
  INTEGER                              :: NNLITTLEVS   ! number of litter levels
  INTEGER                              :: NNSOILCARB   ! number of soil carbon pools  
!
!-------------------------------------------------------------------------------
!
! General surface parameters:
!
  REAL, POINTER, DIMENSION(:)   :: XZS               ! relief                                  (m)
  REAL, POINTER, DIMENSION(:,:) :: XCOVER            ! fraction of each ecosystem              (-)
  LOGICAL, POINTER, DIMENSION(:):: LCOVER            ! GCOVER(i)=T --> ith cover field is not 0.
!
! Averaged Surface radiative parameters:
!
  REAL, POINTER, DIMENSION(:)   :: XALBNIR_DRY       ! dry soil near-infra-red albedo          (-)
  REAL, POINTER, DIMENSION(:)   :: XALBVIS_DRY       ! dry soil visible albedo                 (-)
  REAL, POINTER, DIMENSION(:)   :: XALBUV_DRY        ! dry soil UV albedo                      (-)
  REAL, POINTER, DIMENSION(:)   :: XALBNIR_WET       ! wet soil near-infra-red albedo          (-)
  REAL, POINTER, DIMENSION(:)   :: XALBVIS_WET       ! wet soil visible albedo                 (-)
  REAL, POINTER, DIMENSION(:)   :: XALBUV_WET        ! wet soil UV albedo                      (-)
  REAL, POINTER, DIMENSION(:,:) :: XALBNIR_SOIL      ! soil near-infra-red albedo              (-)
  REAL, POINTER, DIMENSION(:,:) :: XALBVIS_SOIL      ! soil visible albedo                     (-)
  REAL, POINTER, DIMENSION(:,:) :: XALBUV_SOIL       ! soil UV albedo                          (-)
  REAL, POINTER, DIMENSION(:)   :: XALBNIR_TSOIL     ! total near-infra-red albedo of wet soil (-)
  REAL, POINTER, DIMENSION(:)   :: XALBVIS_TSOIL     ! total visible albedo of soil            (-)
!
!
! Subgrid orography parameters
!
  REAL, DIMENSION(:), POINTER :: XAOSIP,XAOSIM,XAOSJP,XAOSJM
! directional A/S quantities in 4 coordinate directions
! (IP: i index up;  IM: i index down;  JP: j index up;  JM: j index down)
! They are used in soil routines to compute effective roughness length
!
  REAL, DIMENSION(:), POINTER :: XHO2IP,XHO2IM,XHO2JP,XHO2JM
! directional h/2 quantities in 4 coordinate directions
! (IP: i index up;  IM: i index down;  JP: j index up;  JM: j index down)
! They are used in soil routines to compute effective roughness length
!
  REAL, DIMENSION(:,:), POINTER :: XZ0EFFIP,XZ0EFFIM,XZ0EFFJP,XZ0EFFJM
! directional total roughness lenghts in 4 coordinate directions
! (IP: i index up;  IM: i index down;  JP: j index up;  JM: j index down)
!
  REAL, DIMENSION(:), POINTER   :: XZ0EFFJPDIR    ! heading of J direction (deg from N clockwise)
!
  REAL, DIMENSION(:), POINTER   :: XZ0REL         ! relief roughness length                 (m)
!
  REAL, DIMENSION(:), POINTER   :: XSSO_SLOPE     ! slope of S.S.O.                         (-)
  REAL, DIMENSION(:), POINTER   :: XSSO_STDEV     ! relief  standard deviation              (m)
!-------------------------------------------------------------------------------
!
! Input Parameters, per patch:
!
! - vegetation + bare soil:
!
  REAL, POINTER, DIMENSION(:,:) :: XZ0_O_Z0H         ! ratio of surface roughness lengths
!                                                    ! (momentum to heat)                      (-)
  REAL, POINTER, DIMENSION(:,:) :: XALBNIR           ! near-infra-red albedo                   (-)
  REAL, POINTER, DIMENSION(:,:) :: XALBVIS           ! visible albedo                          (-)
  REAL, POINTER, DIMENSION(:,:) :: XALBUV            ! UV albedo                               (-)
  REAL, POINTER, DIMENSION(:,:) :: XEMIS             ! surface emissivity                      (-)
  REAL, POINTER, DIMENSION(:,:) :: XZ0               ! surface roughness length                (m)
!
! - vegetation:
!
  REAL, POINTER, DIMENSION(:,:) :: XALBNIR_VEG       ! vegetation near-infra-red albedo        (-)
  REAL, POINTER, DIMENSION(:,:) :: XALBVIS_VEG       ! vegetation visible albedo               (-)
  REAL, POINTER, DIMENSION(:,:) :: XALBUV_VEG        ! vegetation UV albedo                    (-)
  REAL, POINTER, DIMENSION(:)   :: XALBNIR_TVEG      ! total near-infra-red albedo of vegetation (-)
  REAL, POINTER, DIMENSION(:)   :: XALBVIS_TVEG      ! total visible albedo of vegetation        (-)  
!
! - vegetation: default option (Jarvis) and general parameters:
!
  REAL, POINTER, DIMENSION(:,:) :: XVEG              ! vegetation cover fraction               (-)
  REAL, POINTER, DIMENSION(:,:) :: XWRMAX_CF         ! coefficient for maximum water 
!                                                      ! interception 
!                                                      ! storage capacity on the vegetation      (-)
  REAL, POINTER, DIMENSION(:,:) :: XRSMIN            ! minimum stomatal resistance             (s/m)
  REAL, POINTER, DIMENSION(:,:) :: XGAMMA            ! coefficient for the calculation
!                                                      ! of the surface stomatal
!                                                      ! resistance
  REAL, POINTER, DIMENSION(:,:) :: XCV               ! vegetation thermal inertia coefficient  (K m2/J)
  REAL, POINTER, DIMENSION(:,:) :: XRGL              ! maximum solar radiation
!                                                      ! usable in photosynthesis                (W/m2)
  REAL, POINTER, DIMENSION(:,:,:) :: XROOTFRAC       ! root fraction profile ('DIF' option)
!
!-------------------------------------------------------------------------------
!
! - vegetation: Ags parameters ('AGS', 'LAI', 'AST', 'LST', 'NIT', 'NCB' options)
!
  REAL, POINTER, DIMENSION(:)      :: XABC          ! abscissa needed for integration
!                                                   ! of net assimilation and stomatal
!                                                   ! conductance over canopy depth           (-)
  REAL, POINTER, DIMENSION(:)      :: XPOI          ! Gaussian weights for integration
!                                                   ! of net assimilation and stomatal
!                                                   ! conductance over canopy depth           (-)
  REAL, POINTER, DIMENSION(:,:)    :: XBSLAI        ! ratio d(biomass)/d(lai)                 (kg/m2)
  REAL, POINTER, DIMENSION(:,:)    :: XLAIMIN       ! minimum LAI (Leaf Area Index)           (m2/m2)
  REAL, POINTER, DIMENSION(:,:)    :: XSEFOLD       ! e-folding time for senescence           (s)
  REAL, POINTER, DIMENSION(:,:)    :: XTAU_WOOD     ! residence time in woody biomass         (s)  
  REAL, POINTER, DIMENSION(:,:)    :: XH_TREE       ! height of trees                         (m)
  REAL, POINTER, DIMENSION(:,:)    :: XANF          ! total assimilation over canopy          (
  REAL, POINTER, DIMENSION(:,:)    :: XANMAX        ! maximum photosynthesis rate             (
  REAL, POINTER, DIMENSION(:,:)    :: XFZERO        ! ideal value of F, no photo- 
!                                                     ! respiration or saturation deficit       (
  REAL, POINTER, DIMENSION(:,:)    :: XEPSO         ! maximum initial quantum use             
!                                                     ! efficiency                              (mg J-1 PAR)
  REAL, POINTER, DIMENSION(:,:)    :: XGAMM         ! CO2 conpensation concentration          (ppm)
  REAL, POINTER, DIMENSION(:,:)    :: XQDGAMM       ! Q10 function for CO2 conpensation 
!                                                     ! concentration                           (-)
  REAL, POINTER, DIMENSION(:,:)    :: XGMES         ! mesophyll conductance                   (m s-1)
  REAL, POINTER, DIMENSION(:,:)    :: XRE25         ! Ecosystem respiration parameter         (kg/kg.m.s-1)
  REAL, POINTER, DIMENSION(:,:)    :: XQDGMES       ! Q10 function for mesophyll conductance  (-)
  REAL, POINTER, DIMENSION(:,:)    :: XT1GMES       ! reference temperature for computing 
!                                                     ! compensation concentration function for 
!                                                     ! mesophyll conductance: minimum
!                                                     ! temperature                             (K)
  REAL, POINTER, DIMENSION(:,:)    :: XT2GMES       ! reference temperature for computing 
!                                                     ! compensation concentration function for 
!                                                     ! mesophyll conductance: maximum
!                                                     ! temperature                             (K)
  REAL, POINTER, DIMENSION(:,:)    :: XAMAX         ! leaf photosynthetic capacity            (mg m-2 s-1)
  REAL, POINTER, DIMENSION(:,:)    :: XQDAMAX       ! Q10 function for leaf photosynthetic 
!                                                     ! capacity                                (-)
  REAL, POINTER, DIMENSION(:,:)    :: XT1AMAX       ! reference temperature for computing 
!                                                     ! compensation concentration function for 
!                                                     ! leaf photosynthetic capacity: minimum
!                                                     ! temperature                             (K)
  REAL, POINTER, DIMENSION(:,:)    :: XT2AMAX       ! reference temperature for computing 
!                                                     ! compensation concentration function for 
!                                                     ! leaf photosynthetic capacity: maximum
!                                                     ! temperature                             (K)
!                                      
!-------------------------------------------------------------------------------
!
! - vegetation: Ags Stress parameters ('AST', 'LST', 'NIT', 'NCB' options)
!
  LOGICAL, POINTER, DIMENSION(:,:) :: LSTRESS       ! vegetation response type to water
!                                                     ! stress (true:defensive false:offensive) (-)
  REAL, POINTER, DIMENSION(:,:)    :: XF2I          ! critical normilized soil water 
!                                                     ! content for stress parameterisation
  REAL, POINTER, DIMENSION(:,:)    :: XGC           ! cuticular conductance                   (m s-1)
  REAL, POINTER, DIMENSION(:,:)    :: XAH           ! coefficients for herbaceous water stress 
!                                                     ! response (offensive or defensive)       (log(mm/s))
  REAL, POINTER, DIMENSION(:,:)    :: XBH           ! coefficients for herbaceous water stress 
!                                                     ! response (offensive or defensive)       (-)
  REAL, POINTER, DIMENSION(:,:)    :: XDMAX         ! maximum air saturation deficit
!                                                     ! tolerate by vegetation                  (kg/kg)
!
!-------------------------------------------------------------------------------
!
! - vegetation: Ags Nitrogen-model parameters ('NIT', 'NCB' option)
!
  REAL, POINTER, DIMENSION(:,:)    :: XCE_NITRO       ! leaf aera ratio sensitivity to 
!                                                       ! nitrogen concentration                (m2/kg)
  REAL, POINTER, DIMENSION(:,:)    :: XCF_NITRO       ! lethal minimum value of leaf area
!                                                       ! ratio                                 (m2/kg)
  REAL, POINTER, DIMENSION(:,:)    :: XCNA_NITRO      ! nitrogen concentration of active 
!                                                       ! biomass                               (kg/kg)
  REAL, POINTER, DIMENSION(:,:)    :: XBSLAI_NITRO   ! biomass/LAI ratio from nitrogen 
!                                                       ! decline theory                        (kg/m2)
!
!-------------------------------------------------------------------------------
!
! - soil: primary parameters
!
  REAL, POINTER, DIMENSION(:,:)    :: XSAND          ! sand fraction                           (-)
  REAL, POINTER, DIMENSION(:,:)    :: XCLAY          ! clay fraction                           (-)
  REAL, POINTER, DIMENSION(:)      :: XWDRAIN        ! continuous drainage parameter           (-)
  REAL, POINTER, DIMENSION(:)      :: XTAUICE        ! soil freezing characteristic timescale  (s)
  REAL, POINTER, DIMENSION(:)      :: XGAMMAT        ! 'Force-Restore' timescale when using a
!                                                      ! prescribed lower boundary temperature   (1/days)
  REAL, POINTER, DIMENSION(:,:,:)  :: XDG            ! soil layer thicknesses                  (m)
!                                                      ! NOTE: in Force-Restore mode, the 
!                                                      ! uppermost layer thickness is superficial
!                                                      ! and is only explicitly used for soil 
!                                                      ! water phase changes                     (m)
!
  REAL, POINTER, DIMENSION(:,:,:)  :: XDZG           ! soil layers thicknesses (DIF option)
  REAL, POINTER, DIMENSION(:,:,:)  :: XDZDIF         ! distance between consecuative layer mid-points (DIF option)
!
  INTEGER, POINTER, DIMENSION(:,:) :: NWG_LAYER      ! Number of soil moisture layers for DIF
  REAL, POINTER, DIMENSION(:,:)    :: XDROOT         ! effective root depth for DIF (m)
  REAL, POINTER, DIMENSION(:,:)    :: XDG2           ! root depth for DIF as 3-L (m)
!
!-------------------------------------------------------------------------------
!
! - soil: Secondary parameters: hydrology
!
  REAL, POINTER, DIMENSION(:,:)    :: XC1SAT         ! 'Force-Restore' C1 coefficient at 
!                                                    ! saturation                              (-)
  REAL, POINTER, DIMENSION(:,:)    :: XC2REF         ! 'Force-Restore' reference value of C2   (-)
  REAL, POINTER, DIMENSION(:,:,:)  :: XC3            ! 'Force-Restore' C3 drainage coefficient (m)
  REAL, POINTER, DIMENSION(:)      :: XC4B           ! 'Force-Restore' sub-surface vertical 
!                                                    ! diffusion coefficient (slope parameter) (-)
  REAL, POINTER, DIMENSION(:,:)    :: XC4REF         ! 'Force-Restore' sub-surface vertical 
!                                                    ! diffusion coefficient                   (-)
  REAL, POINTER, DIMENSION(:)      :: XACOEF         ! 'Force-Restore' surface vertical 
!                                                    ! diffusion coefficient                   (-)
  REAL, POINTER, DIMENSION(:)      :: XPCOEF         ! 'Force-Restore' surface vertical 
!                                                    ! diffusion coefficient                   (-)
  REAL, POINTER, DIMENSION(:,:)    :: XWFC           ! field capacity volumetric water content
!                                                      ! profile                                 (m3/m3)
  REAL, POINTER, DIMENSION(:,:)    :: XWWILT         ! wilting point volumetric water content 
!                                                      ! profile                                 (m3/m3)
  REAL, POINTER, DIMENSION(:,:)    :: XWSAT          ! porosity profile                        (m3/m3) 
  REAL, POINTER, DIMENSION(:,:)    :: XBCOEF         ! soil water CH78 b-parameter             (-)
  REAL, POINTER, DIMENSION(:,:,:)  :: XCONDSAT       ! hydraulic conductivity at saturation    (m/s)
  REAL, POINTER, DIMENSION(:,:)    :: XMPOTSAT       ! matric potential at saturation          (m)
!
!-------------------------------------------------------------------------------
!
! - soil: Secondary parameters: thermal 
!
  REAL, POINTER, DIMENSION(:)      :: XCGSAT         ! soil thermal inertia coefficient at 
!                                                      ! saturation                              (K m2/J)
  REAL, POINTER, DIMENSION(:,:)    :: XHCAPSOIL      ! soil heat capacity                      (J/K/m3)
  REAL, POINTER, DIMENSION(:,:)    :: XCONDDRY       ! soil dry thermal conductivity           (W/m/K)
  REAL, POINTER, DIMENSION(:,:)    :: XCONDSLD       ! soil solids thermal conductivity        (W/m/K)
  REAL, POINTER, DIMENSION(:)      :: XTDEEP         ! prescribed deep soil temperature 
!                                                      ! (optional)                              (K)
!-------------------------------------------------------------------------------
!
! Prognostic variables:
!
! - Snow Cover:
!
  TYPE(SURF_SNOW)                       :: TSNOW         ! snow state: 
!                                                      ! scheme type/option                      (-)
!                                                      ! number of layers                        (-)
!                                                      ! snow (& liq. water) content             (kg/m2)
!                                                      ! heat content                            (J/m2)
!                                                      ! temperature                             (K)
!                                                      ! density                                 (kg m-3)
!
!-------------------------------------------------------------------------------
!
! - Soil and vegetation heat and water:
!
  REAL, POINTER, DIMENSION(:,:)     :: XWR           ! liquid water retained on the
!                                                      ! foliage of the vegetation
!                                                      ! canopy                                  (kg/m2)
  REAL, POINTER, DIMENSION(:,:,:)   :: XTG           ! surface and sub-surface soil 
!                                                      ! temperature profile                     (K)
  REAL, POINTER, DIMENSION(:,:,:)   :: XWG           ! soil volumetric water content profile   (m3/m3)
  REAL, POINTER, DIMENSION(:,:,:)   :: XWGI          ! soil liquid water equivalent volumetric 
!                                                      ! ice content profile                     (m3/m3)
  REAL, POINTER, DIMENSION(:,:)     :: XRESA         ! aerodynamic resistance                  (s/m)

  REAL, POINTER, DIMENSION(:,:)     :: XPCPS
  REAL, POINTER, DIMENSION(:,:)     :: XPLVTT
  REAL, POINTER, DIMENSION(:,:)     :: XPLSTT 
!
!-------------------------------------------------------------------------------
!
! - Vegetation: Ags Prognostic (YPHOTO = ('LAI', 'LST', or 'NIT') or prescribed (YPHOTO='NON', 'AGS' or 'LST')
!
  REAL, POINTER, DIMENSION(:,:)     :: XLAI          ! Leaf Area Index                         (m2/m2)
!
!-------------------------------------------------------------------------------
!
! - Vegetation: Ags Prognostic (YPHOTO = 'AGS', 'LAI', 'AST', 'LST', 'NIT', 'NCB')
!
  REAL, POINTER, DIMENSION(:,:)     :: XAN           ! net CO2 assimilation                    (mg/m2/s)
  REAL, POINTER, DIMENSION(:,:)     :: XANDAY        ! daily net CO2 assimilation              (mg/m2)
  REAL, POINTER, DIMENSION(:,:)     :: XANFM         ! maximum leaf assimilation               (mg/m2/s)
  REAL, POINTER, DIMENSION(:,:)     :: XLE           ! evapotranspiration                      (W/m2)
  REAL, POINTER, DIMENSION(:,:)     :: XFAPARC       ! Fapar of vegetation (cumul)
  REAL, POINTER, DIMENSION(:,:)     :: XFAPIRC       ! Fapir of vegetation (cumul)
  REAL, POINTER, DIMENSION(:,:)     :: XLAI_EFFC     ! Effective LAI (cumul)
  REAL, POINTER, DIMENSION(:,:)     :: XMUS          ! cos zenithal angle (cumul)     
!
!-------------------------------------------------------------------------------
!
! - Vegetation: Ags Prognostic (YPHOTO = 'NIT', 'NCB')
!
  REAL, POINTER, DIMENSION(:,:,:)   :: XRESP_BIOMASS    ! daily cumulated respiration of 
!                                                       ! biomass                              (kg/m2/s)
  REAL, POINTER, DIMENSION(:,:,:)   :: XBIOMASS         ! biomass of previous day              (kg/m2) 
  REAL, POINTER, DIMENSION(:,:,:)   :: XINCREASE        ! biomass increase                     (kg/m2/day)
!
!
!-------------------------------------------------------------------------------
!
! - Soil carbon (ISBA-CC, YRESPSL = 'CNT')
!
  REAL, POINTER, DIMENSION(:,:,:,:) :: XLITTER          ! litter pools                         (gC/m2)
  REAL, POINTER, DIMENSION(:,:,:)   :: XSOILCARB        ! soil carbon pools                    (gC/m2) 
  REAL, POINTER, DIMENSION(:,:,:)   :: XLIGNIN_STRUC    ! ratio Lignin/Carbon in structural
!                                                         litter                               (gC/m2)
!
  REAL, POINTER, DIMENSION(:,:,:)   :: XTURNOVER        ! turnover rates from biomass to litter (gC/m2/s)
!
!-------------------------------------------------------------------------------
!
! - Irrigation, seeding and reaping
!
  TYPE (DATE_TIME), POINTER, DIMENSION(:,:)  :: TSEED          ! date of seeding
  TYPE (DATE_TIME), POINTER, DIMENSION(:,:)  :: TREAP          ! date of reaping
  REAL, POINTER, DIMENSION(:,:)         :: XWATSUP        ! water supply during irrigation process (mm)
  REAL, POINTER, DIMENSION(:,:)         :: XIRRIG         ! flag for irrigation (irrigation if >0.)
!-------------------------------------------------------------------------------
!
! - Adjustable physical parameters
!
  REAL                                  :: XCGMAX           ! maximum soil heat capacity
!
  REAL                                  :: XCDRAG           ! drag coefficient in canopy
!-------------------------------------------------------------------------------
!
! - SGH scheme
!                                                     
  CHARACTER(LEN=4)               :: CRUNOFF! surface runoff formulation
!                                          ! 'WSAT'
!                                          ! 'DT92'
!                                          ! 'SGH ' Topmodel
!                                                     
  CHARACTER(LEN=3)               :: CTOPREG! Wolock and McCabe (2000) linear regression for Topmodel
                                           ! 'DEF' = Reg
                                           ! 'NON' = no Reg  
!                                           
  CHARACTER(LEN=3)               :: CKSAT  ! ksat
!                                          ! 'DEF' = default value 
!                                          ! 'SGH' = profil exponentiel
!                                           
  CHARACTER(LEN=3)               :: CSOM   ! soil organic matter effect
!                                          ! 'DEF' = default value 
!                                          ! 'SGH' = soil organic matter profil
!
  CHARACTER(LEN=3)               :: CHORT  ! Horton runoff
                                           ! 'DEF' = no Horton runoff
                                           ! 'SGH' = Horton runoff
!
!
  INTEGER                          :: NLAYER_HORT
  INTEGER                          :: NLAYER_DUN
!
  REAL, POINTER, DIMENSION(:)      :: XRUNOFFB       ! sub-grid dt92 surface runoff slope parameter (-)
  REAL, POINTER, DIMENSION(:,:)    :: XRUNOFFD       ! depth over which sub-grid runoff is
!                                                      ! computed: in Force-Restore this is the
!                                                      ! total soil column ('2-L'), or root zone
!                                                      ! ('3-L'). For the 'DIF' option, it can
!                                                      ! be any depth within soil column         (m)

!
  REAL, POINTER, DIMENSION(:,:,:)  :: XSOILWGHT      ! ISBA-DIF: weights for vertical
!                                        
  REAL, POINTER, DIMENSION(:,:)  :: XD_ICE    !depth of the soil column for the calculation
!                                              of the frozen soil fraction (m)
  REAL, POINTER, DIMENSION(:,:)  :: XKSAT_ICE !hydraulic conductivity at saturation
!                                              over frozen area (m s-1)                                     
!-------------------------------------------------------------------------------
!
! - Snow and flood fractions and total albedo at time t:
!
  REAL, POINTER, DIMENSION(:,:) :: XPSNG         ! Snow fraction over ground
  REAL, POINTER, DIMENSION(:,:) :: XPSNV         ! Snow fraction over vegetation
  REAL, POINTER, DIMENSION(:,:) :: XPSNV_A       ! Snow fraction over vegetation
  REAL, POINTER, DIMENSION(:,:) :: XPSN          ! Total Snow fraction
! 
!-------------------------------------------------------------------------------
!                                 
! Type of vegetation (simplification of vegetation charaterization)
CHARACTER(LEN=4)             :: CTYPE_HVEG   ! type of high vegetation
CHARACTER(LEN=4)             :: CTYPE_LVEG   ! type of low vegetation
CHARACTER(LEN=4)             :: CTYPE_NVEG   ! type of bare soil (no vegetation)
!-------------------------------------------------------------------------------
!
END TYPE TEB_GARDEN_t

TYPE(TEB_GARDEN_t), ALLOCATABLE, TARGET, SAVE :: TEB_GARDEN_MODEL(:)

CHARACTER(LEN=4), POINTER :: CROUGH=>NULL()
!$OMP THREADPRIVATE(CROUGH)
CHARACTER(LEN=3), POINTER :: CISBA=>NULL()
!$OMP THREADPRIVATE(CISBA)
CHARACTER(LEN=4), POINTER :: CPEDOTF=>NULL()
!$OMP THREADPRIVATE(CPEDOTF)
CHARACTER(LEN=3), POINTER :: CPHOTO=>NULL()
!$OMP THREADPRIVATE(CPHOTO)
LOGICAL,          POINTER :: LTR_ML=>NULL()
!$OMP THREADPRIVATE(LTR_ML)
CHARACTER(LEN=4), POINTER :: CALBEDO=>NULL()
!$OMP THREADPRIVATE(CALBEDO)
CHARACTER(LEN=4), POINTER :: CRUNOFF=>NULL()
!$OMP THREADPRIVATE(CRUNOFF)
CHARACTER(LEN=4), POINTER :: CSCOND=>NULL()
!$OMP THREADPRIVATE(CSCOND)
CHARACTER(LEN=4), POINTER :: CC1DRY=>NULL()
!$OMP THREADPRIVATE(CC1DRY)
CHARACTER(LEN=3), POINTER :: CSOILFRZ=>NULL()
!$OMP THREADPRIVATE(CSOILFRZ)
CHARACTER(LEN=4), POINTER :: CDIFSFCOND=>NULL()
!$OMP THREADPRIVATE(CDIFSFCOND)
CHARACTER(LEN=3), POINTER :: CSNOWRES=>NULL()
!$OMP THREADPRIVATE(CSNOWRES)
CHARACTER(LEN=3), POINTER :: CRESPSL=>NULL()
!$OMP THREADPRIVATE(CRESPSL)
CHARACTER(LEN=3), POINTER :: CCPSURF=>NULL()
!$OMP THREADPRIVATE(CCPSURF)
LOGICAL, POINTER :: LCANOPY=>NULL()
!$OMP THREADPRIVATE(LCANOPY)
LOGICAL, POINTER :: LCANOPY_DRAG=>NULL()
!$OMP THREADPRIVATE(LCANOPY_DRAG)
LOGICAL, POINTER :: LPAR_GARDEN=>NULL()
!$OMP THREADPRIVATE(LPAR_GARDEN)
INTEGER, POINTER, DIMENSION(:)   :: NSIZE_NATURE_P=>NULL()
!$OMP THREADPRIVATE(NSIZE_NATURE_P)
INTEGER, POINTER, DIMENSION(:,:) :: NR_NATURE_P=>NULL()
!$OMP THREADPRIVATE(NR_NATURE_P)
REAL, POINTER, DIMENSION(:,:)    :: XPATCH=>NULL()
!$OMP THREADPRIVATE(XPATCH)
REAL, POINTER, DIMENSION(:,:)    :: XVEGTYPE=>NULL()
!$OMP THREADPRIVATE(XVEGTYPE)
REAL, POINTER, DIMENSION(:,:,:)  :: XVEGTYPE_PATCH=>NULL()
!$OMP THREADPRIVATE(XVEGTYPE_PATCH)
INTEGER, POINTER :: NPATCH=>NULL()
!$OMP THREADPRIVATE(NPATCH)
INTEGER, POINTER :: NGROUND_LAYER=>NULL()
!$OMP THREADPRIVATE(NGROUND_LAYER)
REAL, POINTER, DIMENSION(:)   :: XSOILGRID=>NULL()
!$OMP THREADPRIVATE(XSOILGRID)
INTEGER, POINTER :: NNBIOMASS=>NULL()
!$OMP THREADPRIVATE(NNBIOMASS)
INTEGER, POINTER :: NNLITTER=>NULL()
!$OMP THREADPRIVATE(NNLITTER)
INTEGER, POINTER :: NNLITTLEVS=>NULL()
!$OMP THREADPRIVATE(NNLITTLEVS)
INTEGER, POINTER :: NNSOILCARB=>NULL()
!$OMP THREADPRIVATE(NNSOILCARB)
REAL, POINTER, DIMENSION(:)   :: XZS=>NULL()
!$OMP THREADPRIVATE(XZS)
REAL, POINTER, DIMENSION(:,:) :: XCOVER=>NULL()
!$OMP THREADPRIVATE(XCOVER)
LOGICAL, POINTER, DIMENSION(:):: LCOVER=>NULL()
!$OMP THREADPRIVATE(LCOVER)
REAL, POINTER, DIMENSION(:)   :: XALBNIR_DRY=>NULL()
!$OMP THREADPRIVATE(XALBNIR_DRY)
REAL, POINTER, DIMENSION(:)   :: XALBVIS_DRY=>NULL()
!$OMP THREADPRIVATE(XALBVIS_DRY)
REAL, POINTER, DIMENSION(:)   :: XALBUV_DRY=>NULL()
!$OMP THREADPRIVATE(XALBUV_DRY)
REAL, POINTER, DIMENSION(:)   :: XALBNIR_WET=>NULL()
!$OMP THREADPRIVATE(XALBNIR_WET)
REAL, POINTER, DIMENSION(:)   :: XALBVIS_WET=>NULL()
!$OMP THREADPRIVATE(XALBVIS_WET)
REAL, POINTER, DIMENSION(:)   :: XALBUV_WET=>NULL()
!$OMP THREADPRIVATE(XALBUV_WET)
REAL, POINTER, DIMENSION(:,:) :: XALBNIR_SOIL=>NULL()
!$OMP THREADPRIVATE(XALBNIR_SOIL)
REAL, POINTER, DIMENSION(:,:) :: XALBVIS_SOIL=>NULL()
!$OMP THREADPRIVATE(XALBVIS_SOIL)
REAL, POINTER, DIMENSION(:,:) :: XALBUV_SOIL=>NULL()
!$OMP THREADPRIVATE(XALBUV_SOIL)
REAL, POINTER, DIMENSION(:)   :: XALBNIR_TSOIL=>NULL()
!$OMP THREADPRIVATE(XALBNIR_TSOIL)
REAL, POINTER, DIMENSION(:)   :: XALBVIS_TSOIL=>NULL()
!$OMP THREADPRIVATE(XALBVIS_TSOIL)
REAL, DIMENSION(:), POINTER :: XAOSIP=>NULL(),XAOSIM=>NULL(),XAOSJP=>NULL(),XAOSJM=>NULL()
!$OMP THREADPRIVATE(XAOSIP,XAOSIM,XAOSJP,XAOSJM)
REAL, DIMENSION(:), POINTER :: XHO2IP=>NULL(),XHO2IM=>NULL(),XHO2JP=>NULL(),XHO2JM=>NULL()
!$OMP THREADPRIVATE(XHO2IP,XHO2IM,XHO2JP,XHO2JM)
REAL, DIMENSION(:,:), POINTER :: XZ0EFFIP=>NULL(),XZ0EFFIM=>NULL(),XZ0EFFJP=>NULL(),XZ0EFFJM=>NULL()
!$OMP THREADPRIVATE(XZ0EFFIP,XZ0EFFIM,XZ0EFFJP,XZ0EFFJM)
REAL, DIMENSION(:), POINTER   :: XZ0EFFJPDIR=>NULL()
!$OMP THREADPRIVATE(XZ0EFFJPDIR)
REAL, DIMENSION(:), POINTER   :: XZ0REL=>NULL()
!$OMP THREADPRIVATE(XZ0REL)
REAL, DIMENSION(:), POINTER   :: XSSO_SLOPE=>NULL()
!$OMP THREADPRIVATE(XSSO_SLOPE)
REAL, DIMENSION(:), POINTER   :: XSSO_STDEV=>NULL()
!$OMP THREADPRIVATE(XSSO_STDEV)
REAL, POINTER, DIMENSION(:,:) :: XZ0_O_Z0H=>NULL()
!$OMP THREADPRIVATE(XZ0_O_Z0H)
REAL, POINTER, DIMENSION(:,:) :: XALBNIR=>NULL()
!$OMP THREADPRIVATE(XALBNIR)
REAL, POINTER, DIMENSION(:,:) :: XALBVIS=>NULL()
!$OMP THREADPRIVATE(XALBVIS)
REAL, POINTER, DIMENSION(:,:) :: XALBUV=>NULL()
!$OMP THREADPRIVATE(XALBUV)
REAL, POINTER, DIMENSION(:,:) :: XEMIS=>NULL()
!$OMP THREADPRIVATE(XEMIS)
REAL, POINTER, DIMENSION(:,:) :: XZ0=>NULL()
!$OMP THREADPRIVATE(XZ0)
REAL, POINTER, DIMENSION(:,:) :: XALBNIR_VEG=>NULL()
!$OMP THREADPRIVATE(XALBNIR_VEG)
REAL, POINTER, DIMENSION(:,:) :: XALBVIS_VEG=>NULL()
!$OMP THREADPRIVATE(XALBVIS_VEG)
REAL, POINTER, DIMENSION(:,:) :: XALBUV_VEG=>NULL()
!$OMP THREADPRIVATE(XALBUV_VEG)
REAL, POINTER, DIMENSION(:)   :: XALBNIR_TVEG=>NULL()
!$OMP THREADPRIVATE(XALBNIR_TVEG)
REAL, POINTER, DIMENSION(:)   :: XALBVIS_TVEG=>NULL()
!$OMP THREADPRIVATE(XALBVIS_TVEG)
REAL, POINTER, DIMENSION(:,:) :: XVEG=>NULL()
!$OMP THREADPRIVATE(XVEG)
REAL, POINTER, DIMENSION(:,:) :: XWRMAX_CF=>NULL()
!$OMP THREADPRIVATE(XWRMAX_CF)
REAL, POINTER, DIMENSION(:,:) :: XRSMIN=>NULL()
!$OMP THREADPRIVATE(XRSMIN)
REAL, POINTER, DIMENSION(:,:) :: XGAMMA=>NULL()
!$OMP THREADPRIVATE(XGAMMA)
REAL, POINTER, DIMENSION(:,:) :: XCV=>NULL()
!$OMP THREADPRIVATE(XCV)
REAL, POINTER, DIMENSION(:,:) :: XRGL=>NULL()
!$OMP THREADPRIVATE(XRGL)
REAL, POINTER, DIMENSION(:,:,:) :: XROOTFRAC=>NULL()
!$OMP THREADPRIVATE(XROOTFRAC)
REAL, POINTER, DIMENSION(:) :: XABC=>NULL()
!$OMP THREADPRIVATE(XABC)
REAL, POINTER, DIMENSION(:) :: XPOI=>NULL()
!$OMP THREADPRIVATE(XPOI)
REAL, POINTER, DIMENSION(:,:)    :: XBSLAI=>NULL()
!$OMP THREADPRIVATE(XBSLAI)
REAL, POINTER, DIMENSION(:,:)    :: XLAIMIN=>NULL()
!$OMP THREADPRIVATE(XLAIMIN)
REAL, POINTER, DIMENSION(:,:)    :: XSEFOLD=>NULL()
!$OMP THREADPRIVATE(XSEFOLD)
REAL, POINTER, DIMENSION(:,:)    :: XTAU_WOOD=>NULL()
!$OMP THREADPRIVATE(XTAU_WOOD)
REAL, POINTER, DIMENSION(:,:)    :: XH_TREE=>NULL()
!$OMP THREADPRIVATE(XH_TREE)
REAL, POINTER, DIMENSION(:,:)    :: XANF=>NULL()
!$OMP THREADPRIVATE(XANF)
REAL, POINTER, DIMENSION(:,:)    :: XANMAX=>NULL()
!$OMP THREADPRIVATE(XANMAX)
REAL, POINTER, DIMENSION(:,:)    :: XFZERO=>NULL()
!$OMP THREADPRIVATE(XFZERO)
REAL, POINTER, DIMENSION(:,:)    :: XEPSO=>NULL()
!$OMP THREADPRIVATE(XEPSO)
REAL, POINTER, DIMENSION(:,:)    :: XGAMM=>NULL()
!$OMP THREADPRIVATE(XGAMM)
REAL, POINTER, DIMENSION(:,:)    :: XQDGAMM=>NULL()
!$OMP THREADPRIVATE(XQDGAMM)
REAL, POINTER, DIMENSION(:,:)    :: XGMES=>NULL()
!$OMP THREADPRIVATE(XGMES)
REAL, POINTER, DIMENSION(:,:)    :: XRE25=>NULL()
!$OMP THREADPRIVATE(XRE25)
REAL, POINTER, DIMENSION(:,:)    :: XQDGMES=>NULL()
!$OMP THREADPRIVATE(XQDGMES)
REAL, POINTER, DIMENSION(:,:)    :: XT1GMES=>NULL()
!$OMP THREADPRIVATE(XT1GMES)
REAL, POINTER, DIMENSION(:,:)    :: XT2GMES=>NULL()
!$OMP THREADPRIVATE(XT2GMES)
REAL, POINTER, DIMENSION(:,:)    :: XAMAX=>NULL()
!$OMP THREADPRIVATE(XAMAX)
REAL, POINTER, DIMENSION(:,:)    :: XQDAMAX=>NULL()
!$OMP THREADPRIVATE(XQDAMAX)
REAL, POINTER, DIMENSION(:,:)    :: XT1AMAX=>NULL()
!$OMP THREADPRIVATE(XT1AMAX)
REAL, POINTER, DIMENSION(:,:)    :: XT2AMAX=>NULL()
!$OMP THREADPRIVATE(XT2AMAX)
LOGICAL, POINTER, DIMENSION(:,:) :: LSTRESS=>NULL()
!$OMP THREADPRIVATE(LSTRESS)
REAL, POINTER, DIMENSION(:,:)    :: XF2I=>NULL()
!$OMP THREADPRIVATE(XF2I)
REAL, POINTER, DIMENSION(:,:)    :: XGC=>NULL()
!$OMP THREADPRIVATE(XGC)
REAL, POINTER, DIMENSION(:,:)    :: XAH=>NULL()
!$OMP THREADPRIVATE(XAH)
REAL, POINTER, DIMENSION(:,:)    :: XBH=>NULL()
!$OMP THREADPRIVATE(XBH)
REAL, POINTER, DIMENSION(:,:)    :: XDMAX=>NULL()
!$OMP THREADPRIVATE(XDMAX)
REAL, POINTER, DIMENSION(:,:)    :: XCE_NITRO=>NULL()
!$OMP THREADPRIVATE(XCE_NITRO)
REAL, POINTER, DIMENSION(:,:)    :: XCF_NITRO=>NULL()
!$OMP THREADPRIVATE(XCF_NITRO)
REAL, POINTER, DIMENSION(:,:)    :: XCNA_NITRO=>NULL()
!$OMP THREADPRIVATE(XCNA_NITRO)
REAL, POINTER, DIMENSION(:,:)    :: XBSLAI_NITRO=>NULL()
!$OMP THREADPRIVATE(XBSLAI_NITRO)
REAL, POINTER, DIMENSION(:,:)    :: XSAND=>NULL()
!$OMP THREADPRIVATE(XSAND)
REAL, POINTER, DIMENSION(:,:)    :: XCLAY=>NULL()
!$OMP THREADPRIVATE(XCLAY)
REAL, POINTER, DIMENSION(:)      :: XRUNOFFB=>NULL()
!$OMP THREADPRIVATE(XRUNOFFB)
REAL, POINTER, DIMENSION(:)      :: XWDRAIN=>NULL()
!$OMP THREADPRIVATE(XWDRAIN)
REAL, POINTER, DIMENSION(:)      :: XTAUICE=>NULL()
!$OMP THREADPRIVATE(XTAUICE)
REAL, POINTER, DIMENSION(:)      :: XGAMMAT=>NULL()
!$OMP THREADPRIVATE(XGAMMAT)
REAL, POINTER, DIMENSION(:,:,:)  :: XDG=>NULL()
!$OMP THREADPRIVATE(XDG)
REAL, POINTER, DIMENSION(:,:,:)  :: XDZG=>NULL()
!$OMP THREADPRIVATE(XDZG)
REAL, POINTER, DIMENSION(:,:,:)  :: XDZDIF=>NULL()
!$OMP THREADPRIVATE(XDZDIF)
INTEGER, POINTER, DIMENSION(:,:) :: NWG_LAYER=>NULL()
!$OMP THREADPRIVATE(NWG_LAYER)
REAL, POINTER, DIMENSION(:,:)    :: XDROOT=>NULL()
!$OMP THREADPRIVATE(XDROOT)
REAL, POINTER, DIMENSION(:,:)    :: XDG2=>NULL()
!$OMP THREADPRIVATE(XDG2)
REAL, POINTER, DIMENSION(:,:)    :: XRUNOFFD=>NULL()
!$OMP THREADPRIVATE(XRUNOFFD)
REAL, POINTER, DIMENSION(:,:,:)  :: XSOILWGHT=>NULL()
!$OMP THREADPRIVATE(XSOILWGHT)
REAL, POINTER, DIMENSION(:,:)    :: XC1SAT=>NULL()
!$OMP THREADPRIVATE(XC1SAT)
REAL, POINTER, DIMENSION(:,:)    :: XC2REF=>NULL()
!$OMP THREADPRIVATE(XC2REF)
REAL, POINTER, DIMENSION(:,:,:)  :: XC3=>NULL()
!$OMP THREADPRIVATE(XC3)
REAL, POINTER, DIMENSION(:)      :: XC4B=>NULL()
!$OMP THREADPRIVATE(XC4B)
REAL, POINTER, DIMENSION(:,:)    :: XC4REF=>NULL()
!$OMP THREADPRIVATE(XC4REF)
REAL, POINTER, DIMENSION(:)      :: XACOEF=>NULL()
!$OMP THREADPRIVATE(XACOEF)
REAL, POINTER, DIMENSION(:)      :: XPCOEF=>NULL()
!$OMP THREADPRIVATE(XPCOEF)
REAL, POINTER, DIMENSION(:,:)    :: XWFC=>NULL()
!$OMP THREADPRIVATE(XWFC)
REAL, POINTER, DIMENSION(:,:)    :: XWWILT=>NULL()
!$OMP THREADPRIVATE(XWWILT)
REAL, POINTER, DIMENSION(:,:)    :: XWSAT=>NULL()
!$OMP THREADPRIVATE(XWSAT)
REAL, POINTER, DIMENSION(:,:)    :: XBCOEF=>NULL()
!$OMP THREADPRIVATE(XBCOEF)
REAL, POINTER, DIMENSION(:,:,:)  :: XCONDSAT=>NULL()
!$OMP THREADPRIVATE(XCONDSAT)
REAL, POINTER, DIMENSION(:,:)    :: XMPOTSAT=>NULL()
!$OMP THREADPRIVATE(XMPOTSAT)
REAL, POINTER, DIMENSION(:)      :: XCGSAT=>NULL()
!$OMP THREADPRIVATE(XCGSAT)
REAL, POINTER, DIMENSION(:,:)    :: XHCAPSOIL=>NULL()
!$OMP THREADPRIVATE(XHCAPSOIL)
REAL, POINTER, DIMENSION(:,:)    :: XCONDDRY=>NULL()
!$OMP THREADPRIVATE(XCONDDRY)
REAL, POINTER, DIMENSION(:,:)    :: XCONDSLD=>NULL()
!$OMP THREADPRIVATE(XCONDSLD)
REAL, POINTER, DIMENSION(:)      :: XTDEEP=>NULL()
!$OMP THREADPRIVATE(XTDEEP)
TYPE(SURF_SNOW), POINTER :: TSNOW=>NULL()
!$OMP THREADPRIVATE(TSNOW)
REAL, POINTER, DIMENSION(:,:)     :: XWR=>NULL()
!$OMP THREADPRIVATE(XWR)
REAL, POINTER, DIMENSION(:,:,:)   :: XTG=>NULL()
!$OMP THREADPRIVATE(XTG)
REAL, POINTER, DIMENSION(:,:,:)   :: XWG=>NULL()
!$OMP THREADPRIVATE(XWG)
REAL, POINTER, DIMENSION(:,:,:)   :: XWGI=>NULL()
!$OMP THREADPRIVATE(XWGI)
REAL, POINTER, DIMENSION(:,:)     :: XRESA=>NULL()
!$OMP THREADPRIVATE(XRESA)
REAL, POINTER, DIMENSION(:,:)     :: XPCPS=>NULL()
!$OMP THREADPRIVATE(XPCPS)
REAL, POINTER, DIMENSION(:,:)     :: XPLVTT=>NULL()
!$OMP THREADPRIVATE(XPLVTT)
REAL, POINTER, DIMENSION(:,:)     :: XPLSTT=>NULL()
!$OMP THREADPRIVATE(XPLSTT)
REAL, POINTER, DIMENSION(:,:)     :: XLAI=>NULL()
!$OMP THREADPRIVATE(XLAI)
REAL, POINTER, DIMENSION(:,:)     :: XAN=>NULL()
!$OMP THREADPRIVATE(XAN)
REAL, POINTER, DIMENSION(:,:)     :: XANDAY=>NULL()
!$OMP THREADPRIVATE(XANDAY)
REAL, POINTER, DIMENSION(:,:)     :: XANFM=>NULL()
!$OMP THREADPRIVATE(XANFM)
REAL, POINTER, DIMENSION(:,:)     :: XLE=>NULL()
!$OMP THREADPRIVATE(XLE)
REAL, POINTER, DIMENSION(:,:)     :: XFAPARC=>NULL()
!$OMP THREADPRIVATE(XFAPARC)
REAL, POINTER, DIMENSION(:,:)     :: XFAPIRC=>NULL()
!$OMP THREADPRIVATE(XFAPIRC)
REAL, POINTER, DIMENSION(:,:)     :: XLAI_EFFC=>NULL()
!$OMP THREADPRIVATE(XLAI_EFFC)
REAL, POINTER, DIMENSION(:,:)     :: XMUS=>NULL()
!$OMP THREADPRIVATE(XMUS)
REAL, POINTER, DIMENSION(:,:,:)   :: XRESP_BIOMASS=>NULL()
!$OMP THREADPRIVATE(XRESP_BIOMASS)
REAL, POINTER, DIMENSION(:,:,:)   :: XBIOMASS=>NULL()
!$OMP THREADPRIVATE(XBIOMASS)
REAL, POINTER, DIMENSION(:,:,:)   :: XINCREASE=>NULL()
!$OMP THREADPRIVATE(XINCREASE)
REAL, POINTER, DIMENSION(:,:,:,:) :: XLITTER=>NULL()
!$OMP THREADPRIVATE(XLITTER)
REAL, POINTER, DIMENSION(:,:,:)   :: XSOILCARB=>NULL()
!$OMP THREADPRIVATE(XSOILCARB)
REAL, POINTER, DIMENSION(:,:,:)   :: XLIGNIN_STRUC=>NULL()
!$OMP THREADPRIVATE(XLIGNIN_STRUC)
REAL, POINTER, DIMENSION(:,:,:)   :: XTURNOVER=>NULL()
!$OMP THREADPRIVATE(XTURNOVER)
TYPE (DATE_TIME), POINTER, DIMENSION(:,:) :: TSEED=>NULL()
!$OMP THREADPRIVATE(TSEED)
TYPE (DATE_TIME), POINTER, DIMENSION(:,:) :: TREAP=>NULL()
!$OMP THREADPRIVATE(TREAP)
REAL, POINTER, DIMENSION(:,:)        :: XWATSUP=>NULL()
!$OMP THREADPRIVATE(XWATSUP)
REAL, POINTER, DIMENSION(:,:)        :: XIRRIG=>NULL()
!$OMP THREADPRIVATE(XIRRIG)
REAL, POINTER :: XCGMAX=>NULL()
!$OMP THREADPRIVATE(XCGMAX)
REAL, POINTER :: XCDRAG=>NULL()
!$OMP THREADPRIVATE(XCDRAG)
!
REAL, POINTER, DIMENSION(:,:) :: XPSNG=>NULL()
!$OMP THREADPRIVATE(XPSNG)
REAL, POINTER, DIMENSION(:,:) :: XPSNV=>NULL()
!$OMP THREADPRIVATE(XPSNV)
REAL, POINTER, DIMENSION(:,:) :: XPSNV_A=>NULL()
!$OMP THREADPRIVATE(XPSNV_A)
REAL, POINTER, DIMENSION(:,:) :: XPSN=>NULL()
!$OMP THREADPRIVATE(XPSN)
!
!SGH scheme
!
CHARACTER(LEN=3), POINTER      :: CTOPREG=>NULL()
!$OMP THREADPRIVATE(CTOPREG)
CHARACTER(LEN=3), POINTER      :: CKSAT=>NULL()
!$OMP THREADPRIVATE(CKSAT)
CHARACTER(LEN=3), POINTER      :: CSOM=>NULL()
!$OMP THREADPRIVATE(CSOM)
CHARACTER(LEN=3), POINTER      :: CHORT=>NULL()
!$OMP THREADPRIVATE(CHORT)
!
REAL, POINTER, DIMENSION(:,:)  :: XD_ICE=>NULL()
!$OMP THREADPRIVATE(XD_ICE)
REAL, POINTER, DIMENSION(:,:)  :: XKSAT_ICE=>NULL()
!$OMP THREADPRIVATE(XKSAT_ICE)
!
INTEGER, POINTER :: NLAYER_HORT=>NULL()
!$OMP THREADPRIVATE(NLAYER_HORT)
INTEGER, POINTER :: NLAYER_DUN=>NULL()
!$OMP THREADPRIVATE(NLAYER_DUN)
!
! Type of vegetation (simplification of veg characterization)
CHARACTER(LEN=4), POINTER :: CTYPE_HVEG=>NULL()
!$OMP THREADPRIVATE(CTYPE_HVEG)
CHARACTER(LEN=4), POINTER :: CTYPE_LVEG=>NULL()
!$OMP THREADPRIVATE(CTYPE_LVEG)
CHARACTER(LEN=4), POINTER :: CTYPE_NVEG=>NULL()
!$OMP THREADPRIVATE(CTYPE_NVEG)
!
CONTAINS

SUBROUTINE TEB_GARDEN_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Save current state for allocated arrays
IF (LKFROM) THEN
TEB_GARDEN_MODEL(KFROM)%NSIZE_NATURE_P=>NSIZE_NATURE_P
TEB_GARDEN_MODEL(KFROM)%NR_NATURE_P=>NR_NATURE_P
TEB_GARDEN_MODEL(KFROM)%XPATCH=>XPATCH
TEB_GARDEN_MODEL(KFROM)%XVEGTYPE=>XVEGTYPE
TEB_GARDEN_MODEL(KFROM)%XVEGTYPE_PATCH=>XVEGTYPE_PATCH
TEB_GARDEN_MODEL(KFROM)%XSOILGRID=>XSOILGRID
TEB_GARDEN_MODEL(KFROM)%XZS=>XZS
TEB_GARDEN_MODEL(KFROM)%XCOVER=>XCOVER
TEB_GARDEN_MODEL(KFROM)%LCOVER=>LCOVER
TEB_GARDEN_MODEL(KFROM)%XALBNIR_DRY=>XALBNIR_DRY
TEB_GARDEN_MODEL(KFROM)%XALBVIS_DRY=>XALBVIS_DRY
TEB_GARDEN_MODEL(KFROM)%XALBUV_DRY=>XALBUV_DRY
TEB_GARDEN_MODEL(KFROM)%XALBNIR_WET=>XALBNIR_WET
TEB_GARDEN_MODEL(KFROM)%XALBVIS_WET=>XALBVIS_WET
TEB_GARDEN_MODEL(KFROM)%XALBUV_WET=>XALBUV_WET
TEB_GARDEN_MODEL(KFROM)%XALBNIR_SOIL=>XALBNIR_SOIL
TEB_GARDEN_MODEL(KFROM)%XALBVIS_SOIL=>XALBVIS_SOIL
TEB_GARDEN_MODEL(KFROM)%XALBUV_SOIL=>XALBUV_SOIL
TEB_GARDEN_MODEL(KFROM)%XALBNIR_TSOIL=>XALBNIR_TSOIL
TEB_GARDEN_MODEL(KFROM)%XALBVIS_TSOIL=>XALBVIS_TSOIL
TEB_GARDEN_MODEL(KFROM)%XAOSIP=>XAOSIP
TEB_GARDEN_MODEL(KFROM)%XAOSIM=>XAOSIM
TEB_GARDEN_MODEL(KFROM)%XAOSJP=>XAOSJP
TEB_GARDEN_MODEL(KFROM)%XAOSJM=>XAOSJM
TEB_GARDEN_MODEL(KFROM)%XHO2IP=>XHO2IP
TEB_GARDEN_MODEL(KFROM)%XHO2IM=>XHO2IM
TEB_GARDEN_MODEL(KFROM)%XHO2JP=>XHO2JP
TEB_GARDEN_MODEL(KFROM)%XHO2JM=>XHO2JM
TEB_GARDEN_MODEL(KFROM)%XZ0EFFIP=>XZ0EFFIP
TEB_GARDEN_MODEL(KFROM)%XZ0EFFIM=>XZ0EFFIM
TEB_GARDEN_MODEL(KFROM)%XZ0EFFJP=>XZ0EFFJP
TEB_GARDEN_MODEL(KFROM)%XZ0EFFJM=>XZ0EFFJM
TEB_GARDEN_MODEL(KFROM)%XZ0EFFJPDIR=>XZ0EFFJPDIR
TEB_GARDEN_MODEL(KFROM)%XZ0REL=>XZ0REL
TEB_GARDEN_MODEL(KFROM)%XSSO_SLOPE=>XSSO_SLOPE
TEB_GARDEN_MODEL(KFROM)%XSSO_STDEV=>XSSO_STDEV
TEB_GARDEN_MODEL(KFROM)%XZ0_O_Z0H=>XZ0_O_Z0H
TEB_GARDEN_MODEL(KFROM)%XALBNIR=>XALBNIR
TEB_GARDEN_MODEL(KFROM)%XALBVIS=>XALBVIS
TEB_GARDEN_MODEL(KFROM)%XALBUV=>XALBUV
TEB_GARDEN_MODEL(KFROM)%XEMIS=>XEMIS
TEB_GARDEN_MODEL(KFROM)%XZ0=>XZ0
TEB_GARDEN_MODEL(KFROM)%XALBNIR_VEG=>XALBNIR_VEG
TEB_GARDEN_MODEL(KFROM)%XALBVIS_VEG=>XALBVIS_VEG
TEB_GARDEN_MODEL(KFROM)%XALBUV_VEG=>XALBUV_VEG
TEB_GARDEN_MODEL(KFROM)%XALBNIR_TVEG=>XALBNIR_TVEG
TEB_GARDEN_MODEL(KFROM)%XALBVIS_TVEG=>XALBVIS_TVEG
TEB_GARDEN_MODEL(KFROM)%XVEG=>XVEG
TEB_GARDEN_MODEL(KFROM)%XWRMAX_CF=>XWRMAX_CF
TEB_GARDEN_MODEL(KFROM)%XRSMIN=>XRSMIN
TEB_GARDEN_MODEL(KFROM)%XGAMMA=>XGAMMA
TEB_GARDEN_MODEL(KFROM)%XCV=>XCV
TEB_GARDEN_MODEL(KFROM)%XRGL=>XRGL
TEB_GARDEN_MODEL(KFROM)%XROOTFRAC=>XROOTFRAC
TEB_GARDEN_MODEL(KFROM)%XABC=>XABC
TEB_GARDEN_MODEL(KFROM)%XPOI=>XPOI
TEB_GARDEN_MODEL(KFROM)%XBSLAI=>XBSLAI
TEB_GARDEN_MODEL(KFROM)%XLAIMIN=>XLAIMIN
TEB_GARDEN_MODEL(KFROM)%XSEFOLD=>XSEFOLD
TEB_GARDEN_MODEL(KFROM)%XTAU_WOOD=>XTAU_WOOD
TEB_GARDEN_MODEL(KFROM)%XH_TREE=>XH_TREE
TEB_GARDEN_MODEL(KFROM)%XANF=>XANF
TEB_GARDEN_MODEL(KFROM)%XANMAX=>XANMAX
TEB_GARDEN_MODEL(KFROM)%XFZERO=>XFZERO
TEB_GARDEN_MODEL(KFROM)%XEPSO=>XEPSO
TEB_GARDEN_MODEL(KFROM)%XGAMM=>XGAMM
TEB_GARDEN_MODEL(KFROM)%XQDGAMM=>XQDGAMM
TEB_GARDEN_MODEL(KFROM)%XGMES=>XGMES
TEB_GARDEN_MODEL(KFROM)%XRE25=>XRE25
TEB_GARDEN_MODEL(KFROM)%XQDGMES=>XQDGMES
TEB_GARDEN_MODEL(KFROM)%XT1GMES=>XT1GMES
TEB_GARDEN_MODEL(KFROM)%XT2GMES=>XT2GMES
TEB_GARDEN_MODEL(KFROM)%XAMAX=>XAMAX
TEB_GARDEN_MODEL(KFROM)%XQDAMAX=>XQDAMAX
TEB_GARDEN_MODEL(KFROM)%XT1AMAX=>XT1AMAX
TEB_GARDEN_MODEL(KFROM)%XT2AMAX=>XT2AMAX
TEB_GARDEN_MODEL(KFROM)%LSTRESS=>LSTRESS
TEB_GARDEN_MODEL(KFROM)%XF2I=>XF2I
TEB_GARDEN_MODEL(KFROM)%XGC=>XGC
TEB_GARDEN_MODEL(KFROM)%XAH=>XAH
TEB_GARDEN_MODEL(KFROM)%XBH=>XBH
TEB_GARDEN_MODEL(KFROM)%XDMAX=>XDMAX
TEB_GARDEN_MODEL(KFROM)%XCE_NITRO=>XCE_NITRO
TEB_GARDEN_MODEL(KFROM)%XCF_NITRO=>XCF_NITRO
TEB_GARDEN_MODEL(KFROM)%XCNA_NITRO=>XCNA_NITRO
TEB_GARDEN_MODEL(KFROM)%XBSLAI_NITRO=>XBSLAI_NITRO
TEB_GARDEN_MODEL(KFROM)%XSAND=>XSAND
TEB_GARDEN_MODEL(KFROM)%XCLAY=>XCLAY
TEB_GARDEN_MODEL(KFROM)%XRUNOFFB=>XRUNOFFB
TEB_GARDEN_MODEL(KFROM)%XWDRAIN=>XWDRAIN
TEB_GARDEN_MODEL(KFROM)%XTAUICE=>XTAUICE
TEB_GARDEN_MODEL(KFROM)%XGAMMAT=>XGAMMAT
TEB_GARDEN_MODEL(KFROM)%XDG=>XDG
TEB_GARDEN_MODEL(KFROM)%XDZG=>XDZG
TEB_GARDEN_MODEL(KFROM)%XDZDIF=>XDZDIF
TEB_GARDEN_MODEL(KFROM)%NWG_LAYER=>NWG_LAYER
TEB_GARDEN_MODEL(KFROM)%XDROOT=>XDROOT
TEB_GARDEN_MODEL(KFROM)%XDG2=>XDG2
TEB_GARDEN_MODEL(KFROM)%XRUNOFFD=>XRUNOFFD
TEB_GARDEN_MODEL(KFROM)%XSOILWGHT=>XSOILWGHT
TEB_GARDEN_MODEL(KFROM)%XC1SAT=>XC1SAT
TEB_GARDEN_MODEL(KFROM)%XC2REF=>XC2REF
TEB_GARDEN_MODEL(KFROM)%XC3=>XC3
TEB_GARDEN_MODEL(KFROM)%XC4B=>XC4B
TEB_GARDEN_MODEL(KFROM)%XC4REF=>XC4REF
TEB_GARDEN_MODEL(KFROM)%XACOEF=>XACOEF
TEB_GARDEN_MODEL(KFROM)%XPCOEF=>XPCOEF
TEB_GARDEN_MODEL(KFROM)%XWFC=>XWFC
TEB_GARDEN_MODEL(KFROM)%XWWILT=>XWWILT
TEB_GARDEN_MODEL(KFROM)%XWSAT=>XWSAT
TEB_GARDEN_MODEL(KFROM)%XBCOEF=>XBCOEF
TEB_GARDEN_MODEL(KFROM)%XCONDSAT=>XCONDSAT
TEB_GARDEN_MODEL(KFROM)%XMPOTSAT=>XMPOTSAT
TEB_GARDEN_MODEL(KFROM)%XCGSAT=>XCGSAT
TEB_GARDEN_MODEL(KFROM)%XHCAPSOIL=>XHCAPSOIL
TEB_GARDEN_MODEL(KFROM)%XCONDDRY=>XCONDDRY
TEB_GARDEN_MODEL(KFROM)%XCONDSLD=>XCONDSLD
TEB_GARDEN_MODEL(KFROM)%XTDEEP=>XTDEEP
TEB_GARDEN_MODEL(KFROM)%XWR=>XWR
TEB_GARDEN_MODEL(KFROM)%XTG=>XTG
TEB_GARDEN_MODEL(KFROM)%XWG=>XWG
TEB_GARDEN_MODEL(KFROM)%XWGI=>XWGI
TEB_GARDEN_MODEL(KFROM)%XRESA=>XRESA
TEB_GARDEN_MODEL(KFROM)%XPCPS=>XPCPS
TEB_GARDEN_MODEL(KFROM)%XPLVTT=>XPLVTT
TEB_GARDEN_MODEL(KFROM)%XPLSTT=>XPLSTT
TEB_GARDEN_MODEL(KFROM)%XLAI=>XLAI
TEB_GARDEN_MODEL(KFROM)%XAN=>XAN
TEB_GARDEN_MODEL(KFROM)%XANDAY=>XANDAY
TEB_GARDEN_MODEL(KFROM)%XANFM=>XANFM
TEB_GARDEN_MODEL(KFROM)%XLE=>XLE
TEB_GARDEN_MODEL(KFROM)%XFAPARC=>XFAPARC
TEB_GARDEN_MODEL(KFROM)%XFAPIRC=>XFAPIRC
TEB_GARDEN_MODEL(KFROM)%XLAI_EFFC=>XLAI_EFFC
TEB_GARDEN_MODEL(KFROM)%XMUS=>XMUS
TEB_GARDEN_MODEL(KFROM)%XRESP_BIOMASS=>XRESP_BIOMASS
TEB_GARDEN_MODEL(KFROM)%XBIOMASS=>XBIOMASS
TEB_GARDEN_MODEL(KFROM)%XINCREASE=>XINCREASE
TEB_GARDEN_MODEL(KFROM)%XLITTER=>XLITTER
TEB_GARDEN_MODEL(KFROM)%XSOILCARB=>XSOILCARB
TEB_GARDEN_MODEL(KFROM)%XLIGNIN_STRUC=>XLIGNIN_STRUC
TEB_GARDEN_MODEL(KFROM)%XTURNOVER=>XTURNOVER
TEB_GARDEN_MODEL(KFROM)%TSEED=>TSEED
TEB_GARDEN_MODEL(KFROM)%TREAP=>TREAP
TEB_GARDEN_MODEL(KFROM)%XWATSUP=>XWATSUP
TEB_GARDEN_MODEL(KFROM)%XIRRIG=>XIRRIG
!
TEB_GARDEN_MODEL(KFROM)%XPSNG=>XPSNG
TEB_GARDEN_MODEL(KFROM)%XPSNV_A=>XPSNV_A
TEB_GARDEN_MODEL(KFROM)%XPSNV=>XPSNV
TEB_GARDEN_MODEL(KFROM)%XPSN=>XPSN
!
!SGH scheme
!                    
TEB_GARDEN_MODEL(KFROM)%XD_ICE=>XD_ICE
TEB_GARDEN_MODEL(KFROM)%XKSAT_ICE=>XKSAT_ICE
ENDIF
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_TEB_GARDEN_N:TEB_GARDEN_GOTO_MODEL',0,ZHOOK_HANDLE)
CROUGH=>TEB_GARDEN_MODEL(KTO)%CROUGH
CISBA=>TEB_GARDEN_MODEL(KTO)%CISBA
CPEDOTF=>TEB_GARDEN_MODEL(KTO)%CPEDOTF
CPHOTO=>TEB_GARDEN_MODEL(KTO)%CPHOTO
LTR_ML=>TEB_GARDEN_MODEL(KTO)%LTR_ML
CALBEDO=>TEB_GARDEN_MODEL(KTO)%CALBEDO
CRUNOFF=>TEB_GARDEN_MODEL(KTO)%CRUNOFF
CSCOND=>TEB_GARDEN_MODEL(KTO)%CSCOND
CC1DRY=>TEB_GARDEN_MODEL(KTO)%CC1DRY
CSOILFRZ=>TEB_GARDEN_MODEL(KTO)%CSOILFRZ
CDIFSFCOND=>TEB_GARDEN_MODEL(KTO)%CDIFSFCOND
CSNOWRES=>TEB_GARDEN_MODEL(KTO)%CSNOWRES
CRESPSL=>TEB_GARDEN_MODEL(KTO)%CRESPSL
CCPSURF=>TEB_GARDEN_MODEL(KTO)%CCPSURF
LCANOPY=>TEB_GARDEN_MODEL(KTO)%LCANOPY
LCANOPY_DRAG=>TEB_GARDEN_MODEL(KTO)%LCANOPY_DRAG
LPAR_GARDEN=>TEB_GARDEN_MODEL(KTO)%LPAR_GARDEN
NSIZE_NATURE_P=>TEB_GARDEN_MODEL(KTO)%NSIZE_NATURE_P
NR_NATURE_P=>TEB_GARDEN_MODEL(KTO)%NR_NATURE_P
XPATCH=>TEB_GARDEN_MODEL(KTO)%XPATCH
XVEGTYPE=>TEB_GARDEN_MODEL(KTO)%XVEGTYPE
XVEGTYPE_PATCH=>TEB_GARDEN_MODEL(KTO)%XVEGTYPE_PATCH
NPATCH=>TEB_GARDEN_MODEL(KTO)%NPATCH
NGROUND_LAYER=>TEB_GARDEN_MODEL(KTO)%NGROUND_LAYER
XSOILGRID=>TEB_GARDEN_MODEL(KTO)%XSOILGRID
NNBIOMASS=>TEB_GARDEN_MODEL(KTO)%NNBIOMASS
NNLITTER=>TEB_GARDEN_MODEL(KTO)%NNLITTER
NNLITTLEVS=>TEB_GARDEN_MODEL(KTO)%NNLITTLEVS
NNSOILCARB=>TEB_GARDEN_MODEL(KTO)%NNSOILCARB
XZS=>TEB_GARDEN_MODEL(KTO)%XZS
XCOVER=>TEB_GARDEN_MODEL(KTO)%XCOVER
LCOVER=>TEB_GARDEN_MODEL(KTO)%LCOVER
XALBNIR_DRY=>TEB_GARDEN_MODEL(KTO)%XALBNIR_DRY
XALBVIS_DRY=>TEB_GARDEN_MODEL(KTO)%XALBVIS_DRY
XALBUV_DRY=>TEB_GARDEN_MODEL(KTO)%XALBUV_DRY
XALBNIR_WET=>TEB_GARDEN_MODEL(KTO)%XALBNIR_WET
XALBVIS_WET=>TEB_GARDEN_MODEL(KTO)%XALBVIS_WET
XALBUV_WET=>TEB_GARDEN_MODEL(KTO)%XALBUV_WET
XALBNIR_SOIL=>TEB_GARDEN_MODEL(KTO)%XALBNIR_SOIL
XALBVIS_SOIL=>TEB_GARDEN_MODEL(KTO)%XALBVIS_SOIL
XALBUV_SOIL=>TEB_GARDEN_MODEL(KTO)%XALBUV_SOIL
XALBNIR_TSOIL=>TEB_GARDEN_MODEL(KTO)%XALBNIR_TSOIL
XALBVIS_TSOIL=>TEB_GARDEN_MODEL(KTO)%XALBVIS_TSOIL
XAOSIP=>TEB_GARDEN_MODEL(KTO)%XAOSIP
XAOSIM=>TEB_GARDEN_MODEL(KTO)%XAOSIM
XAOSJP=>TEB_GARDEN_MODEL(KTO)%XAOSJP
XAOSJM=>TEB_GARDEN_MODEL(KTO)%XAOSJM
XHO2IP=>TEB_GARDEN_MODEL(KTO)%XHO2IP
XHO2IM=>TEB_GARDEN_MODEL(KTO)%XHO2IM
XHO2JP=>TEB_GARDEN_MODEL(KTO)%XHO2JP
XHO2JM=>TEB_GARDEN_MODEL(KTO)%XHO2JM
XZ0EFFIP=>TEB_GARDEN_MODEL(KTO)%XZ0EFFIP
XZ0EFFIM=>TEB_GARDEN_MODEL(KTO)%XZ0EFFIM
XZ0EFFJP=>TEB_GARDEN_MODEL(KTO)%XZ0EFFJP
XZ0EFFJM=>TEB_GARDEN_MODEL(KTO)%XZ0EFFJM
XZ0EFFJPDIR=>TEB_GARDEN_MODEL(KTO)%XZ0EFFJPDIR
XZ0REL=>TEB_GARDEN_MODEL(KTO)%XZ0REL
XSSO_SLOPE=>TEB_GARDEN_MODEL(KTO)%XSSO_SLOPE
XSSO_STDEV=>TEB_GARDEN_MODEL(KTO)%XSSO_STDEV
XZ0_O_Z0H=>TEB_GARDEN_MODEL(KTO)%XZ0_O_Z0H
XALBNIR=>TEB_GARDEN_MODEL(KTO)%XALBNIR
XALBVIS=>TEB_GARDEN_MODEL(KTO)%XALBVIS
XALBUV=>TEB_GARDEN_MODEL(KTO)%XALBUV
XEMIS=>TEB_GARDEN_MODEL(KTO)%XEMIS
XZ0=>TEB_GARDEN_MODEL(KTO)%XZ0
XALBNIR_VEG=>TEB_GARDEN_MODEL(KTO)%XALBNIR_VEG
XALBVIS_VEG=>TEB_GARDEN_MODEL(KTO)%XALBVIS_VEG
XALBUV_VEG=>TEB_GARDEN_MODEL(KTO)%XALBUV_VEG
XALBNIR_TVEG=>TEB_GARDEN_MODEL(KTO)%XALBNIR_TVEG
XALBVIS_TVEG=>TEB_GARDEN_MODEL(KTO)%XALBVIS_TVEG
XVEG=>TEB_GARDEN_MODEL(KTO)%XVEG
XWRMAX_CF=>TEB_GARDEN_MODEL(KTO)%XWRMAX_CF
XRSMIN=>TEB_GARDEN_MODEL(KTO)%XRSMIN
XGAMMA=>TEB_GARDEN_MODEL(KTO)%XGAMMA
XCV=>TEB_GARDEN_MODEL(KTO)%XCV
XRGL=>TEB_GARDEN_MODEL(KTO)%XRGL
XROOTFRAC=>TEB_GARDEN_MODEL(KTO)%XROOTFRAC
XABC=>TEB_GARDEN_MODEL(KTO)%XABC
XPOI=>TEB_GARDEN_MODEL(KTO)%XPOI
XBSLAI=>TEB_GARDEN_MODEL(KTO)%XBSLAI
XLAIMIN=>TEB_GARDEN_MODEL(KTO)%XLAIMIN
XSEFOLD=>TEB_GARDEN_MODEL(KTO)%XSEFOLD
XTAU_WOOD=>TEB_GARDEN_MODEL(KTO)%XTAU_WOOD
XH_TREE=>TEB_GARDEN_MODEL(KTO)%XH_TREE
XANF=>TEB_GARDEN_MODEL(KTO)%XANF
XANMAX=>TEB_GARDEN_MODEL(KTO)%XANMAX
XFZERO=>TEB_GARDEN_MODEL(KTO)%XFZERO
XEPSO=>TEB_GARDEN_MODEL(KTO)%XEPSO
XGAMM=>TEB_GARDEN_MODEL(KTO)%XGAMM
XQDGAMM=>TEB_GARDEN_MODEL(KTO)%XQDGAMM
XGMES=>TEB_GARDEN_MODEL(KTO)%XGMES
XRE25=>TEB_GARDEN_MODEL(KTO)%XRE25
XQDGMES=>TEB_GARDEN_MODEL(KTO)%XQDGMES
XT1GMES=>TEB_GARDEN_MODEL(KTO)%XT1GMES
XT2GMES=>TEB_GARDEN_MODEL(KTO)%XT2GMES
XAMAX=>TEB_GARDEN_MODEL(KTO)%XAMAX
XQDAMAX=>TEB_GARDEN_MODEL(KTO)%XQDAMAX
XT1AMAX=>TEB_GARDEN_MODEL(KTO)%XT1AMAX
XT2AMAX=>TEB_GARDEN_MODEL(KTO)%XT2AMAX
LSTRESS=>TEB_GARDEN_MODEL(KTO)%LSTRESS
XF2I=>TEB_GARDEN_MODEL(KTO)%XF2I
XGC=>TEB_GARDEN_MODEL(KTO)%XGC
XAH=>TEB_GARDEN_MODEL(KTO)%XAH
XBH=>TEB_GARDEN_MODEL(KTO)%XBH
XDMAX=>TEB_GARDEN_MODEL(KTO)%XDMAX
XCE_NITRO=>TEB_GARDEN_MODEL(KTO)%XCE_NITRO
XCF_NITRO=>TEB_GARDEN_MODEL(KTO)%XCF_NITRO
XCNA_NITRO=>TEB_GARDEN_MODEL(KTO)%XCNA_NITRO
XBSLAI_NITRO=>TEB_GARDEN_MODEL(KTO)%XBSLAI_NITRO
XSAND=>TEB_GARDEN_MODEL(KTO)%XSAND
XCLAY=>TEB_GARDEN_MODEL(KTO)%XCLAY
XRUNOFFB=>TEB_GARDEN_MODEL(KTO)%XRUNOFFB
XWDRAIN=>TEB_GARDEN_MODEL(KTO)%XWDRAIN
XTAUICE=>TEB_GARDEN_MODEL(KTO)%XTAUICE
XGAMMAT=>TEB_GARDEN_MODEL(KTO)%XGAMMAT
XDG=>TEB_GARDEN_MODEL(KTO)%XDG
XDZG=>TEB_GARDEN_MODEL(KTO)%XDZG
XDZDIF=>TEB_GARDEN_MODEL(KTO)%XDZDIF
NWG_LAYER=>TEB_GARDEN_MODEL(KTO)%NWG_LAYER
XDROOT=>TEB_GARDEN_MODEL(KTO)%XDROOT
XDG2=>TEB_GARDEN_MODEL(KTO)%XDG2
XRUNOFFD=>TEB_GARDEN_MODEL(KTO)%XRUNOFFD
XSOILWGHT=>TEB_GARDEN_MODEL(KTO)%XSOILWGHT
XC1SAT=>TEB_GARDEN_MODEL(KTO)%XC1SAT
XC2REF=>TEB_GARDEN_MODEL(KTO)%XC2REF
XC3=>TEB_GARDEN_MODEL(KTO)%XC3
XC4B=>TEB_GARDEN_MODEL(KTO)%XC4B
XC4REF=>TEB_GARDEN_MODEL(KTO)%XC4REF
XACOEF=>TEB_GARDEN_MODEL(KTO)%XACOEF
XPCOEF=>TEB_GARDEN_MODEL(KTO)%XPCOEF
XWFC=>TEB_GARDEN_MODEL(KTO)%XWFC
XWWILT=>TEB_GARDEN_MODEL(KTO)%XWWILT
XWSAT=>TEB_GARDEN_MODEL(KTO)%XWSAT
XBCOEF=>TEB_GARDEN_MODEL(KTO)%XBCOEF
XCONDSAT=>TEB_GARDEN_MODEL(KTO)%XCONDSAT
XMPOTSAT=>TEB_GARDEN_MODEL(KTO)%XMPOTSAT
XCGSAT=>TEB_GARDEN_MODEL(KTO)%XCGSAT
XHCAPSOIL=>TEB_GARDEN_MODEL(KTO)%XHCAPSOIL
XCONDDRY=>TEB_GARDEN_MODEL(KTO)%XCONDDRY
XCONDSLD=>TEB_GARDEN_MODEL(KTO)%XCONDSLD
XTDEEP=>TEB_GARDEN_MODEL(KTO)%XTDEEP
TSNOW=>TEB_GARDEN_MODEL(KTO)%TSNOW
XWR=>TEB_GARDEN_MODEL(KTO)%XWR
XTG=>TEB_GARDEN_MODEL(KTO)%XTG
XWG=>TEB_GARDEN_MODEL(KTO)%XWG
XWGI=>TEB_GARDEN_MODEL(KTO)%XWGI
XRESA=>TEB_GARDEN_MODEL(KTO)%XRESA
XPCPS=>TEB_GARDEN_MODEL(KTO)%XPCPS
XPLVTT=>TEB_GARDEN_MODEL(KTO)%XPLVTT
XPLSTT=>TEB_GARDEN_MODEL(KTO)%XPLSTT
XLAI=>TEB_GARDEN_MODEL(KTO)%XLAI
XAN=>TEB_GARDEN_MODEL(KTO)%XAN
XANDAY=>TEB_GARDEN_MODEL(KTO)%XANDAY
XANFM=>TEB_GARDEN_MODEL(KTO)%XANFM
XLE=>TEB_GARDEN_MODEL(KTO)%XLE
XFAPARC=>TEB_GARDEN_MODEL(KTO)%XFAPARC
XFAPIRC=>TEB_GARDEN_MODEL(KTO)%XFAPIRC
XLAI_EFFC=>TEB_GARDEN_MODEL(KTO)%XLAI_EFFC
XMUS=>TEB_GARDEN_MODEL(KTO)%XMUS
XRESP_BIOMASS=>TEB_GARDEN_MODEL(KTO)%XRESP_BIOMASS
XBIOMASS=>TEB_GARDEN_MODEL(KTO)%XBIOMASS
XINCREASE=>TEB_GARDEN_MODEL(KTO)%XINCREASE
XLITTER=>TEB_GARDEN_MODEL(KTO)%XLITTER
XSOILCARB=>TEB_GARDEN_MODEL(KTO)%XSOILCARB
XLIGNIN_STRUC=>TEB_GARDEN_MODEL(KTO)%XLIGNIN_STRUC
XTURNOVER=>TEB_GARDEN_MODEL(KTO)%XTURNOVER
TSEED=>TEB_GARDEN_MODEL(KTO)%TSEED
TREAP=>TEB_GARDEN_MODEL(KTO)%TREAP
XWATSUP=>TEB_GARDEN_MODEL(KTO)%XWATSUP
XIRRIG=>TEB_GARDEN_MODEL(KTO)%XIRRIG
XCGMAX=>TEB_GARDEN_MODEL(KTO)%XCGMAX
XCDRAG=>TEB_GARDEN_MODEL(KTO)%XCDRAG
!
XPSNG=>TEB_GARDEN_MODEL(KTO)%XPSNG
XPSNV_A=>TEB_GARDEN_MODEL(KTO)%XPSNV_A
XPSNV=>TEB_GARDEN_MODEL(KTO)%XPSNV
XPSN=>TEB_GARDEN_MODEL(KTO)%XPSN
!
!SGH scheme
!
CTOPREG=>TEB_GARDEN_MODEL(KTO)%CTOPREG
CKSAT=>TEB_GARDEN_MODEL(KTO)%CKSAT
CSOM=>TEB_GARDEN_MODEL(KTO)%CSOM
CHORT=>TEB_GARDEN_MODEL(KTO)%CHORT
!
XD_ICE=>TEB_GARDEN_MODEL(KTO)%XD_ICE
XKSAT_ICE=>TEB_GARDEN_MODEL(KTO)%XKSAT_ICE
!
NLAYER_HORT=>TEB_GARDEN_MODEL(KTO)%NLAYER_HORT
NLAYER_DUN=>TEB_GARDEN_MODEL(KTO)%NLAYER_DUN
!
! Types of vegetation (simplification of veg characterization)
!
CTYPE_HVEG=>TEB_GARDEN_MODEL(KTO)%CTYPE_HVEG
CTYPE_LVEG=>TEB_GARDEN_MODEL(KTO)%CTYPE_LVEG
CTYPE_NVEG=>TEB_GARDEN_MODEL(KTO)%CTYPE_NVEG
IF (LHOOK) CALL DR_HOOK('MODD_TEB_GARDEN_N:TEB_GARDEN_GOTO_MODEL',1,ZHOOK_HANDLE)
!
END SUBROUTINE TEB_GARDEN_GOTO_MODEL

SUBROUTINE TEB_GARDEN_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GARDEN_N:TEB_GARDEN_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(TEB_GARDEN_MODEL(KMODEL))
DO J=1,KMODEL
  NULLIFY(TEB_GARDEN_MODEL(J)%NSIZE_NATURE_P)
  NULLIFY(TEB_GARDEN_MODEL(J)%NR_NATURE_P)
  NULLIFY(TEB_GARDEN_MODEL(J)%XPATCH)
  NULLIFY(TEB_GARDEN_MODEL(J)%XVEGTYPE)
  NULLIFY(TEB_GARDEN_MODEL(J)%XVEGTYPE_PATCH)
  NULLIFY(TEB_GARDEN_MODEL(J)%XSOILGRID)
  NULLIFY(TEB_GARDEN_MODEL(J)%XZS)
  NULLIFY(TEB_GARDEN_MODEL(J)%XCOVER)
  NULLIFY(TEB_GARDEN_MODEL(J)%LCOVER)
  NULLIFY(TEB_GARDEN_MODEL(J)%XALBNIR_DRY)
  NULLIFY(TEB_GARDEN_MODEL(J)%XALBVIS_DRY)
  NULLIFY(TEB_GARDEN_MODEL(J)%XALBUV_DRY)
  NULLIFY(TEB_GARDEN_MODEL(J)%XALBNIR_WET)
  NULLIFY(TEB_GARDEN_MODEL(J)%XALBVIS_WET)
  NULLIFY(TEB_GARDEN_MODEL(J)%XALBUV_WET)
  NULLIFY(TEB_GARDEN_MODEL(J)%XALBNIR_SOIL)
  NULLIFY(TEB_GARDEN_MODEL(J)%XALBVIS_SOIL)
  NULLIFY(TEB_GARDEN_MODEL(J)%XALBUV_SOIL)
  NULLIFY(TEB_GARDEN_MODEL(J)%XALBNIR_TSOIL)
  NULLIFY(TEB_GARDEN_MODEL(J)%XALBVIS_TSOIL)    
  NULLIFY(TEB_GARDEN_MODEL(J)%XAOSIP)
  NULLIFY(TEB_GARDEN_MODEL(J)%XAOSIM)
  NULLIFY(TEB_GARDEN_MODEL(J)%XAOSJP)
  NULLIFY(TEB_GARDEN_MODEL(J)%XAOSJM)
  NULLIFY(TEB_GARDEN_MODEL(J)%XHO2IP)
  NULLIFY(TEB_GARDEN_MODEL(J)%XHO2IM)
  NULLIFY(TEB_GARDEN_MODEL(J)%XHO2JP)
  NULLIFY(TEB_GARDEN_MODEL(J)%XHO2JM)
  NULLIFY(TEB_GARDEN_MODEL(J)%XZ0EFFIP)
  NULLIFY(TEB_GARDEN_MODEL(J)%XZ0EFFIM)
  NULLIFY(TEB_GARDEN_MODEL(J)%XZ0EFFJP)
  NULLIFY(TEB_GARDEN_MODEL(J)%XZ0EFFJM)
  NULLIFY(TEB_GARDEN_MODEL(J)%XZ0EFFJPDIR)
  NULLIFY(TEB_GARDEN_MODEL(J)%XZ0REL)
  NULLIFY(TEB_GARDEN_MODEL(J)%XSSO_SLOPE)
  NULLIFY(TEB_GARDEN_MODEL(J)%XSSO_STDEV)
  NULLIFY(TEB_GARDEN_MODEL(J)%XZ0_O_Z0H)
  NULLIFY(TEB_GARDEN_MODEL(J)%XALBNIR)
  NULLIFY(TEB_GARDEN_MODEL(J)%XALBVIS)
  NULLIFY(TEB_GARDEN_MODEL(J)%XALBUV)
  NULLIFY(TEB_GARDEN_MODEL(J)%XEMIS)
  NULLIFY(TEB_GARDEN_MODEL(J)%XZ0)
  NULLIFY(TEB_GARDEN_MODEL(J)%XALBNIR_VEG)
  NULLIFY(TEB_GARDEN_MODEL(J)%XALBVIS_VEG)
  NULLIFY(TEB_GARDEN_MODEL(J)%XALBUV_VEG)
  NULLIFY(TEB_GARDEN_MODEL(J)%XALBNIR_TVEG)
  NULLIFY(TEB_GARDEN_MODEL(J)%XALBVIS_TVEG)    
  NULLIFY(TEB_GARDEN_MODEL(J)%XVEG)
  NULLIFY(TEB_GARDEN_MODEL(J)%XWRMAX_CF)
  NULLIFY(TEB_GARDEN_MODEL(J)%XRSMIN)
  NULLIFY(TEB_GARDEN_MODEL(J)%XGAMMA)
  NULLIFY(TEB_GARDEN_MODEL(J)%XCV)
  NULLIFY(TEB_GARDEN_MODEL(J)%XRGL)
  NULLIFY(TEB_GARDEN_MODEL(J)%XROOTFRAC)
  NULLIFY(TEB_GARDEN_MODEL(J)%XABC)
  NULLIFY(TEB_GARDEN_MODEL(J)%XPOI)  
  NULLIFY(TEB_GARDEN_MODEL(J)%XBSLAI)
  NULLIFY(TEB_GARDEN_MODEL(J)%XLAIMIN)
  NULLIFY(TEB_GARDEN_MODEL(J)%XSEFOLD)
  NULLIFY(TEB_GARDEN_MODEL(J)%XTAU_WOOD)
  NULLIFY(TEB_GARDEN_MODEL(J)%XH_TREE)
  NULLIFY(TEB_GARDEN_MODEL(J)%XANF)
  NULLIFY(TEB_GARDEN_MODEL(J)%XANMAX)
  NULLIFY(TEB_GARDEN_MODEL(J)%XFZERO)
  NULLIFY(TEB_GARDEN_MODEL(J)%XEPSO)
  NULLIFY(TEB_GARDEN_MODEL(J)%XGAMM)
  NULLIFY(TEB_GARDEN_MODEL(J)%XQDGAMM)
  NULLIFY(TEB_GARDEN_MODEL(J)%XGMES)
  NULLIFY(TEB_GARDEN_MODEL(J)%XRE25)
  NULLIFY(TEB_GARDEN_MODEL(J)%XQDGMES)
  NULLIFY(TEB_GARDEN_MODEL(J)%XT1GMES)
  NULLIFY(TEB_GARDEN_MODEL(J)%XT2GMES)
  NULLIFY(TEB_GARDEN_MODEL(J)%XAMAX)
  NULLIFY(TEB_GARDEN_MODEL(J)%XQDAMAX)
  NULLIFY(TEB_GARDEN_MODEL(J)%XT1AMAX)
  NULLIFY(TEB_GARDEN_MODEL(J)%XT2AMAX)
  NULLIFY(TEB_GARDEN_MODEL(J)%LSTRESS)
  NULLIFY(TEB_GARDEN_MODEL(J)%XF2I)
  NULLIFY(TEB_GARDEN_MODEL(J)%XGC)
  NULLIFY(TEB_GARDEN_MODEL(J)%XAH)
  NULLIFY(TEB_GARDEN_MODEL(J)%XBH)
  NULLIFY(TEB_GARDEN_MODEL(J)%XDMAX)
  NULLIFY(TEB_GARDEN_MODEL(J)%XCE_NITRO)
  NULLIFY(TEB_GARDEN_MODEL(J)%XCF_NITRO)
  NULLIFY(TEB_GARDEN_MODEL(J)%XCNA_NITRO)
  NULLIFY(TEB_GARDEN_MODEL(J)%XBSLAI_NITRO)
  NULLIFY(TEB_GARDEN_MODEL(J)%XSAND)
  NULLIFY(TEB_GARDEN_MODEL(J)%XCLAY)
  NULLIFY(TEB_GARDEN_MODEL(J)%XRUNOFFB)
  NULLIFY(TEB_GARDEN_MODEL(J)%XWDRAIN)
  NULLIFY(TEB_GARDEN_MODEL(J)%XTAUICE)
  NULLIFY(TEB_GARDEN_MODEL(J)%XGAMMAT)
  NULLIFY(TEB_GARDEN_MODEL(J)%XDG)
  NULLIFY(TEB_GARDEN_MODEL(J)%XDZG)
  NULLIFY(TEB_GARDEN_MODEL(J)%XDZDIF)
  NULLIFY(TEB_GARDEN_MODEL(J)%NWG_LAYER)
  NULLIFY(TEB_GARDEN_MODEL(J)%XDROOT)
  NULLIFY(TEB_GARDEN_MODEL(J)%XDG2)
  NULLIFY(TEB_GARDEN_MODEL(J)%XRUNOFFD)
  NULLIFY(TEB_GARDEN_MODEL(J)%XSOILWGHT)
  NULLIFY(TEB_GARDEN_MODEL(J)%XC1SAT)
  NULLIFY(TEB_GARDEN_MODEL(J)%XC2REF)
  NULLIFY(TEB_GARDEN_MODEL(J)%XC3)
  NULLIFY(TEB_GARDEN_MODEL(J)%XC4B)
  NULLIFY(TEB_GARDEN_MODEL(J)%XC4REF)
  NULLIFY(TEB_GARDEN_MODEL(J)%XACOEF)
  NULLIFY(TEB_GARDEN_MODEL(J)%XPCOEF)
  NULLIFY(TEB_GARDEN_MODEL(J)%XWFC)
  NULLIFY(TEB_GARDEN_MODEL(J)%XWWILT)
  NULLIFY(TEB_GARDEN_MODEL(J)%XWSAT)
  NULLIFY(TEB_GARDEN_MODEL(J)%XBCOEF)
  NULLIFY(TEB_GARDEN_MODEL(J)%XCONDSAT)
  NULLIFY(TEB_GARDEN_MODEL(J)%XMPOTSAT)
  NULLIFY(TEB_GARDEN_MODEL(J)%XCGSAT)
  NULLIFY(TEB_GARDEN_MODEL(J)%XHCAPSOIL)
  NULLIFY(TEB_GARDEN_MODEL(J)%XCONDDRY)
  NULLIFY(TEB_GARDEN_MODEL(J)%XCONDSLD)
  NULLIFY(TEB_GARDEN_MODEL(J)%XTDEEP)
  NULLIFY(TEB_GARDEN_MODEL(J)%XWR)
  NULLIFY(TEB_GARDEN_MODEL(J)%XTG)
  NULLIFY(TEB_GARDEN_MODEL(J)%XWG)
  NULLIFY(TEB_GARDEN_MODEL(J)%XWGI)
  NULLIFY(TEB_GARDEN_MODEL(J)%XRESA)
  NULLIFY(TEB_GARDEN_MODEL(J)%XPCPS)
  NULLIFY(TEB_GARDEN_MODEL(J)%XPLVTT)
  NULLIFY(TEB_GARDEN_MODEL(J)%XPLSTT)
  NULLIFY(TEB_GARDEN_MODEL(J)%XLAI)
  NULLIFY(TEB_GARDEN_MODEL(J)%XAN)
  NULLIFY(TEB_GARDEN_MODEL(J)%XANDAY)
  NULLIFY(TEB_GARDEN_MODEL(J)%XANFM)
  NULLIFY(TEB_GARDEN_MODEL(J)%XLE)
  NULLIFY(TEB_GARDEN_MODEL(J)%XFAPARC)
  NULLIFY(TEB_GARDEN_MODEL(J)%XFAPIRC)
  NULLIFY(TEB_GARDEN_MODEL(J)%XLAI_EFFC)  
  NULLIFY(TEB_GARDEN_MODEL(J)%XMUS)    
  NULLIFY(TEB_GARDEN_MODEL(J)%XRESP_BIOMASS)
  NULLIFY(TEB_GARDEN_MODEL(J)%XBIOMASS)
  NULLIFY(TEB_GARDEN_MODEL(J)%XINCREASE)
  NULLIFY(TEB_GARDEN_MODEL(J)%XLITTER)
  NULLIFY(TEB_GARDEN_MODEL(J)%XSOILCARB)
  NULLIFY(TEB_GARDEN_MODEL(J)%XLIGNIN_STRUC)
  NULLIFY(TEB_GARDEN_MODEL(J)%XTURNOVER)
  NULLIFY(TEB_GARDEN_MODEL(J)%XWATSUP)
  NULLIFY(TEB_GARDEN_MODEL(J)%XIRRIG)
  NULLIFY(TEB_GARDEN_MODEL(J)%XD_ICE)
  NULLIFY(TEB_GARDEN_MODEL(J)%XKSAT_ICE)
  NULLIFY(TEB_GARDEN_MODEL(J)%XPSNG)
  NULLIFY(TEB_GARDEN_MODEL(J)%XPSNV)
  NULLIFY(TEB_GARDEN_MODEL(J)%XPSNV_A)
  NULLIFY(TEB_GARDEN_MODEL(J)%XPSN)
ENDDO
TEB_GARDEN_MODEL(:)%CROUGH=' '
TEB_GARDEN_MODEL(:)%CISBA=' '
TEB_GARDEN_MODEL(:)%CPEDOTF=' '
TEB_GARDEN_MODEL(:)%CPHOTO=' '
TEB_GARDEN_MODEL(:)%LTR_ML=.FALSE.
TEB_GARDEN_MODEL(:)%CALBEDO=' '
TEB_GARDEN_MODEL(:)%CSCOND=' '
TEB_GARDEN_MODEL(:)%CC1DRY=' '
TEB_GARDEN_MODEL(:)%CSOILFRZ=' '
TEB_GARDEN_MODEL(:)%CDIFSFCOND=' '
TEB_GARDEN_MODEL(:)%CSNOWRES=' '
TEB_GARDEN_MODEL(:)%CRESPSL=' '
TEB_GARDEN_MODEL(:)%CCPSURF=' '
TEB_GARDEN_MODEL(:)%LCANOPY=.FALSE.
TEB_GARDEN_MODEL(:)%LCANOPY_DRAG=.FALSE.
TEB_GARDEN_MODEL(:)%LPAR_GARDEN=.TRUE.
TEB_GARDEN_MODEL(:)%NPATCH=0
TEB_GARDEN_MODEL(:)%NGROUND_LAYER=0
TEB_GARDEN_MODEL(:)%NNBIOMASS=0
TEB_GARDEN_MODEL(:)%NNLITTER=0
TEB_GARDEN_MODEL(:)%NNLITTLEVS=0
TEB_GARDEN_MODEL(:)%NNSOILCARB=0
TEB_GARDEN_MODEL(:)%XCGMAX=0.
TEB_GARDEN_MODEL(:)%XCDRAG=0.
TEB_GARDEN_MODEL(:)%CRUNOFF=' '
TEB_GARDEN_MODEL(:)%CTOPREG=' '
TEB_GARDEN_MODEL(:)%CKSAT=' '
TEB_GARDEN_MODEL(:)%CSOM=' '
TEB_GARDEN_MODEL(:)%CHORT=' '
TEB_GARDEN_MODEL(:)%CTYPE_HVEG=' '
TEB_GARDEN_MODEL(:)%CTYPE_LVEG=' '
TEB_GARDEN_MODEL(:)%CTYPE_NVEG=' '
TEB_GARDEN_MODEL(:)%NLAYER_HORT=0
TEB_GARDEN_MODEL(:)%NLAYER_DUN=0
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GARDEN_N:TEB_GARDEN_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_GARDEN_ALLOC

SUBROUTINE TEB_GARDEN_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GARDEN_N:TEB_GARDEN_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(TEB_GARDEN_MODEL)) DEALLOCATE(TEB_GARDEN_MODEL)
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GARDEN_N:TEB_GARDEN_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_GARDEN_DEALLO

END MODULE MODD_TEB_GARDEN_n
