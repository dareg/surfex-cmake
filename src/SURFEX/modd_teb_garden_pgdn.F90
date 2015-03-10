!##################
MODULE MODD_TEB_GARDEN_PGD_n
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
!!      A. Lemonsu   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2011
!!      V. Masson      06/2013 splits module in 4
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

!-------------------------------------------------------------------------------
TYPE TEB_GARDEN_PGD_t
!-------------------------------------------------------------------------------
!
! Mask and number of grid elements containing patches/tiles:
!
  REAL, POINTER, DIMENSION(:,:) :: XVEGTYPE          ! fraction of each vegetation type for
!                                                    ! each grid mesh                          (-)
!-------------------------------------------------------------------------------
!
! Averaged Surface radiative parameters:
!
  REAL, POINTER, DIMENSION(:)   :: XALBNIR_DRY       ! dry soil near-infra-red albedo          (-)
  REAL, POINTER, DIMENSION(:)   :: XALBVIS_DRY       ! dry soil visible albedo                 (-)
  REAL, POINTER, DIMENSION(:)   :: XALBUV_DRY        ! dry soil UV albedo                      (-)
  REAL, POINTER, DIMENSION(:)   :: XALBNIR_WET       ! wet soil near-infra-red albedo          (-)
  REAL, POINTER, DIMENSION(:)   :: XALBVIS_WET       ! wet soil visible albedo                 (-)
  REAL, POINTER, DIMENSION(:)   :: XALBUV_WET        ! wet soil UV albedo                      (-)
  REAL, POINTER, DIMENSION(:)   :: XALBNIR_SOIL      ! soil near-infra-red albedo              (-)
  REAL, POINTER, DIMENSION(:)   :: XALBVIS_SOIL      ! soil visible albedo                     (-)
  REAL, POINTER, DIMENSION(:)   :: XALBUV_SOIL       ! soil UV albedo                          (-)
!
!-------------------------------------------------------------------------------
!
! Input Parameters, per patch:
!
! - vegetation + bare soil:
!
  REAL, POINTER, DIMENSION(:) :: XZ0_O_Z0H         ! ratio of surface roughness lengths
!                                                    ! (momentum to heat)                      (-)
!
! - vegetation:
!
  REAL, POINTER, DIMENSION(:) :: XALBNIR_VEG       ! vegetation near-infra-red albedo        (-)
  REAL, POINTER, DIMENSION(:) :: XALBVIS_VEG       ! vegetation visible albedo               (-)
  REAL, POINTER, DIMENSION(:) :: XALBUV_VEG        ! vegetation UV albedo                    (-)
!
! - vegetation: default option (Jarvis) and general parameters:
!
  REAL, POINTER, DIMENSION(:) :: XWRMAX_CF         ! coefficient for maximum water 
!                                                      ! interception 
!                                                      ! storage capacity on the vegetation      (-)
  REAL, POINTER, DIMENSION(:) :: XRSMIN            ! minimum stomatal resistance             (s/m)
  REAL, POINTER, DIMENSION(:) :: XGAMMA            ! coefficient for the calculation
!                                                      ! of the surface stomatal
!                                                      ! resistance
  REAL, POINTER, DIMENSION(:) :: XCV               ! vegetation thermal inertia coefficient  (K m2/J)
  REAL, POINTER, DIMENSION(:) :: XRGL              ! maximum solar radiation
!                                                      ! usable in photosynthesis                (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XROOTFRAC       ! root fraction profile ('DIF' option)
!
!-------------------------------------------------------------------------------
!
! - vegetation: Ags parameters ('AGS', 'LAI', 'AST', 'LST', 'NIT', 'NCB' options)
!
!  REAL,              DIMENSION(3)   :: XABC             ! abscissa needed for integration
  REAL, POINTER,     DIMENSION(:)   :: XABC       ! abscissa needed for integration
!                                                 ! of net assimilation and stomatal
!                                                 ! conductance over canopy depth           (-)
!  REAL,              DIMENSION(3)   :: XPOI             ! Gaussian weights for integration
  REAL, POINTER,     DIMENSION(:)   :: XPOI       ! Gaussian weights for integration
!                                                 ! of net assimilation and stomatal
!                                                 ! conductance over canopy depth           (-)
  REAL, POINTER, DIMENSION(:)    :: XBSLAI        ! ratio d(biomass)/d(lai)                 (kg/m2)
  REAL, POINTER, DIMENSION(:)    :: XLAIMIN       ! minimum LAI (Leaf Area Index)           (m2/m2)
  REAL, POINTER, DIMENSION(:)    :: XSEFOLD       ! e-folding time for senescence           (s)
  REAL, POINTER, DIMENSION(:)    :: XH_TREE       ! height of trees                         (m)
  REAL, POINTER, DIMENSION(:)    :: XANF          ! total assimilation over canopy          (
  REAL, POINTER, DIMENSION(:)    :: XANMAX        ! maximum photosynthesis rate             (
  REAL, POINTER, DIMENSION(:)    :: XFZERO        ! ideal value of F, no photo- 
!                                                     ! respiration or saturation deficit       (
  REAL, POINTER, DIMENSION(:)    :: XEPSO         ! maximum initial quantum use             
!                                                     ! efficiency                              (mg J-1 PAR)
  REAL, POINTER, DIMENSION(:)    :: XGAMM         ! CO2 conpensation concentration          (ppm)
  REAL, POINTER, DIMENSION(:)    :: XQDGAMM       ! Log of Q10 function for CO2 conpensation 
!                                                 ! concentration                           (-)
  REAL, POINTER, DIMENSION(:)    :: XGMES         ! mesophyll conductance                   (m s-1)
  REAL, POINTER, DIMENSION(:)    :: XRE25         ! Ecosystem respiration parameter         (kg/kg.m.s-1)
  REAL, POINTER, DIMENSION(:)    :: XQDGMES       ! Log of Q10 function for mesophyll conductance  (-)
  REAL, POINTER, DIMENSION(:)    :: XT1GMES       ! reference temperature for computing 
!                                                     ! compensation concentration function for 
!                                                     ! mesophyll conductance: minimum
!                                                     ! temperature                             (K)
  REAL, POINTER, DIMENSION(:)    :: XT2GMES       ! reference temperature for computing 
!                                                     ! compensation concentration function for 
!                                                     ! mesophyll conductance: maximum
!                                                     ! temperature                             (K)
  REAL, POINTER, DIMENSION(:)    :: XAMAX         ! leaf photosynthetic capacity            (mg m-2 s-1)
  REAL, POINTER, DIMENSION(:)    :: XQDAMAX       ! Log of Q10 function for leaf photosynthetic 
!                                                     ! capacity                                (-)
  REAL, POINTER, DIMENSION(:)    :: XT1AMAX       ! reference temperature for computing 
!                                                     ! compensation concentration function for 
!                                                     ! leaf photosynthetic capacity: minimum
!                                                     ! temperature                             (K)
  REAL, POINTER, DIMENSION(:)    :: XT2AMAX       ! reference temperature for computing 
!                                                     ! compensation concentration function for 
!                                                     ! leaf photosynthetic capacity: maximum
!                                                     ! temperature                             (K)
!                                      

!-------------------------------------------------------------------------------
!
! - vegetation: Ags Stress parameters ('AST', 'LST', 'NIT', 'NCB' options)
!
  LOGICAL, POINTER, DIMENSION(:) :: LSTRESS       ! vegetation response type to water
!                                                     ! stress (true:defensive false:offensive) (-)
  REAL, POINTER, DIMENSION(:)    :: XF2I          ! critical normilized soil water 
!                                                     ! content for stress parameterisation
  REAL, POINTER, DIMENSION(:)    :: XGC           ! cuticular conductance                   (m s-1)
  REAL, POINTER, DIMENSION(:)    :: XAH           ! coefficients for herbaceous water stress 
!                                                     ! response (offensive or defensive)       (log(mm/s))
  REAL, POINTER, DIMENSION(:)    :: XBH           ! coefficients for herbaceous water stress 
!                                                     ! response (offensive or defensive)       (-)
  REAL, POINTER, DIMENSION(:)    :: XDMAX         ! maximum air saturation deficit
!                                                     ! tolerate by vegetation                  (kg/kg)
!
!-------------------------------------------------------------------------------
!
! - vegetation: Ags Nitrogen-model parameters ('NIT', 'NCB' option)
!
  REAL, POINTER, DIMENSION(:)    :: XCE_NITRO       ! leaf aera ratio sensitivity to 
!                                                       ! nitrogen concentration                (m2/kg)
  REAL, POINTER, DIMENSION(:)    :: XCF_NITRO       ! lethal minimum value of leaf area
!                                                       ! ratio                                 (m2/kg)
  REAL, POINTER, DIMENSION(:)    :: XCNA_NITRO      ! nitrogen concentration of active 
!                                                       ! biomass                               (kg/kg)
  REAL, POINTER, DIMENSION(:)    :: XBSLAI_NITRO   ! biomass/LAI ratio from nitrogen 
!                                                       ! decline theory                        (kg/m2)
!
!-------------------------------------------------------------------------------
!
! - soil: primary parameters
!
  REAL, POINTER, DIMENSION(:,:)    :: XSAND          ! sand fraction                           (-)
  REAL, POINTER, DIMENSION(:,:)    :: XCLAY          ! clay fraction                           (-)
  REAL, POINTER, DIMENSION(:)      :: XRUNOFFB       ! sub-grid surface runoff slope parameter (-)
  REAL, POINTER, DIMENSION(:)      :: XWDRAIN        ! continuous drainage parameter           (-)
  REAL, POINTER, DIMENSION(:)      :: XTAUICE        ! soil freezing characteristic timescale  (s)
  REAL, POINTER, DIMENSION(:)      :: XGAMMAT        ! 'Force-Restore' timescale when using a
!                                                    ! prescribed lower boundary temperature   (1/days)
  REAL, POINTER, DIMENSION(:,:)    :: XDG            ! soil layer thicknesses                  (m)
!                                                    ! NOTE: in Force-Restore mode, the 
!                                                    ! uppermost layer thickness is superficial
!                                                    ! and is only explicitly used for soil 
!                                                    ! water phase changes                     (m)
  REAL, POINTER, DIMENSION(:)      :: XRUNOFFD       ! depth over which sub-grid runoff is
!                                                    ! computed: in Force-Restore this is the
!                                                    ! total soil column ('2-L'), or root zone
!                                                    ! ('3-L'). For the 'DIF' option, it can
!                                                    ! be any depth within soil column         (m)
!
  REAL, POINTER, DIMENSION(:,:)  :: XSOILWGHT      ! ISBA-DIF: weights for vertical
  REAL, POINTER, DIMENSION(:,:)  :: XDZG           ! soil layers thicknesses (DIF option)
  REAL, POINTER, DIMENSION(:,:)  :: XDZDIF         ! distance between consecuative layer mid-points (DIF option)
!
  INTEGER, POINTER, DIMENSION(:) :: NWG_LAYER      ! Number of soil moisture layers for DIF
  REAL, POINTER, DIMENSION(:)    :: XDROOT         ! effective root depth for DIF (m)
  REAL, POINTER, DIMENSION(:)    :: XDG2           ! root depth for DIF as 3-L (m)
!-------------------------------------------------------------------------------
!
! - soil: Secondary parameters: hydrology
!
  REAL, POINTER, DIMENSION(:)    :: XC1SAT         ! 'Force-Restore' C1 coefficient at 
!                                                    ! saturation                              (-)
  REAL, POINTER, DIMENSION(:)    :: XC2REF         ! 'Force-Restore' reference value of C2   (-)
  REAL, POINTER, DIMENSION(:,:)  :: XC3            ! 'Force-Restore' C3 drainage coefficient (m)
  REAL, POINTER, DIMENSION(:)      :: XC4B           ! 'Force-Restore' sub-surface vertical 
!                                                    ! diffusion coefficient (slope parameter) (-)
  REAL, POINTER, DIMENSION(:)    :: XC4REF         ! 'Force-Restore' sub-surface vertical 
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
  REAL, POINTER, DIMENSION(:,:)  :: XCONDSAT       ! hydraulic conductivity at saturation    (m/s)
  REAL, POINTER, DIMENSION(:,:)  :: XMPOTSAT       ! matric potential at saturation          (m)
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
!                                                    ! (optional)                              (K) 
  REAL, POINTER, DIMENSION(:)     :: XPCPS
  REAL, POINTER, DIMENSION(:)     :: XPLVTT
  REAL, POINTER, DIMENSION(:)     :: XPLSTT 
!-------------------------------------------------------------------------------
!
! - SGH scheme
!                                                     
  REAL, POINTER, DIMENSION(:)  :: XD_ICE    !depth of the soil column for the calculation
!                                              of the frozen soil fraction (m)
  REAL, POINTER, DIMENSION(:)  :: XKSAT_ICE !hydraulic conductivity at saturation
!                                              over frozen area (m s-1)                                     
!-------------------------------------------------------------------------------
!
! Type of vegetation (simplification of vegetation charaterization)
 CHARACTER(LEN=4)             :: CTYPE_HVEG   ! type of high vegetation
 CHARACTER(LEN=4)             :: CTYPE_LVEG   ! type of low vegetation
 CHARACTER(LEN=4)             :: CTYPE_NVEG   ! type of bare soil (no vegetation)
!-------------------------------------------------------------------------------
!
END TYPE TEB_GARDEN_PGD_t
!

TYPE(TEB_GARDEN_PGD_t),     ALLOCATABLE, TARGET, SAVE :: TEB_GARDEN_PGD_MODEL(:)

REAL, POINTER, DIMENSION(:,:) :: XVEGTYPE=>NULL()
!$OMP THREADPRIVATE(XVEGTYPE)
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
REAL, POINTER, DIMENSION(:) :: XALBNIR_SOIL=>NULL()
!$OMP THREADPRIVATE(XALBNIR_SOIL)
REAL, POINTER, DIMENSION(:) :: XALBVIS_SOIL=>NULL()
!$OMP THREADPRIVATE(XALBVIS_SOIL)
REAL, POINTER, DIMENSION(:) :: XALBUV_SOIL=>NULL()
!$OMP THREADPRIVATE(XALBUV_SOIL)
REAL, POINTER, DIMENSION(:) :: XZ0_O_Z0H=>NULL()
!$OMP THREADPRIVATE(XZ0_O_Z0H)
REAL, POINTER, DIMENSION(:) :: XALBNIR_VEG=>NULL()
!$OMP THREADPRIVATE(XALBNIR_VEG)
REAL, POINTER, DIMENSION(:) :: XALBVIS_VEG=>NULL()
!$OMP THREADPRIVATE(XALBVIS_VEG)
REAL, POINTER, DIMENSION(:) :: XALBUV_VEG=>NULL()
!$OMP THREADPRIVATE(XALBUV_VEG)
REAL, POINTER, DIMENSION(:) :: XWRMAX_CF=>NULL()
!$OMP THREADPRIVATE(XWRMAX_CF)
REAL, POINTER, DIMENSION(:) :: XRSMIN=>NULL()
!$OMP THREADPRIVATE(XRSMIN)
REAL, POINTER, DIMENSION(:) :: XGAMMA=>NULL()
!$OMP THREADPRIVATE(XGAMMA)
REAL, POINTER, DIMENSION(:) :: XCV=>NULL()
!$OMP THREADPRIVATE(XCV)
REAL, POINTER, DIMENSION(:) :: XRGL=>NULL()
!$OMP THREADPRIVATE(XRGL)
REAL, POINTER, DIMENSION(:,:) :: XROOTFRAC=>NULL()
!$OMP THREADPRIVATE(XROOTFRAC)
REAL, DIMENSION(:), POINTER :: XABC=>NULL()
!$OMP THREADPRIVATE(XABC)
REAL, DIMENSION(:), POINTER :: XPOI=>NULL()
!$OMP THREADPRIVATE(XPOI)
REAL, POINTER, DIMENSION(:)    :: XBSLAI=>NULL()
!$OMP THREADPRIVATE(XBSLAI)
REAL, POINTER, DIMENSION(:)    :: XLAIMIN=>NULL()
!$OMP THREADPRIVATE(XLAIMIN)
REAL, POINTER, DIMENSION(:)    :: XSEFOLD=>NULL()
!$OMP THREADPRIVATE(XSEFOLD)
REAL, POINTER, DIMENSION(:)    :: XH_TREE=>NULL()
!$OMP THREADPRIVATE(XH_TREE)
REAL, POINTER, DIMENSION(:)    :: XANF=>NULL()
!$OMP THREADPRIVATE(XANF)
REAL, POINTER, DIMENSION(:)    :: XANMAX=>NULL()
!$OMP THREADPRIVATE(XANMAX)
REAL, POINTER, DIMENSION(:)    :: XFZERO=>NULL()
!$OMP THREADPRIVATE(XFZERO)
REAL, POINTER, DIMENSION(:)    :: XEPSO=>NULL()
!$OMP THREADPRIVATE(XEPSO)
REAL, POINTER, DIMENSION(:)    :: XGAMM=>NULL()
!$OMP THREADPRIVATE(XGAMM)
REAL, POINTER, DIMENSION(:)    :: XQDGAMM=>NULL()
!$OMP THREADPRIVATE(XQDGAMM)
REAL, POINTER, DIMENSION(:)    :: XGMES=>NULL()
!$OMP THREADPRIVATE(XGMES)
REAL, POINTER, DIMENSION(:)    :: XRE25=>NULL()
!$OMP THREADPRIVATE(XRE25)
REAL, POINTER, DIMENSION(:)    :: XQDGMES=>NULL()
!$OMP THREADPRIVATE(XQDGMES)
REAL, POINTER, DIMENSION(:)    :: XT1GMES=>NULL()
!$OMP THREADPRIVATE(XT1GMES)
REAL, POINTER, DIMENSION(:)    :: XT2GMES=>NULL()
!$OMP THREADPRIVATE(XT2GMES)
REAL, POINTER, DIMENSION(:)    :: XAMAX=>NULL()
!$OMP THREADPRIVATE(XAMAX)
REAL, POINTER, DIMENSION(:)    :: XQDAMAX=>NULL()
!$OMP THREADPRIVATE(XQDAMAX)
REAL, POINTER, DIMENSION(:)    :: XT1AMAX=>NULL()
!$OMP THREADPRIVATE(XT1AMAX)
REAL, POINTER, DIMENSION(:)    :: XT2AMAX=>NULL()
!$OMP THREADPRIVATE(XT2AMAX)
LOGICAL, POINTER, DIMENSION(:) :: LSTRESS=>NULL()
!$OMP THREADPRIVATE(LSTRESS)
REAL, POINTER, DIMENSION(:)    :: XF2I=>NULL()
!$OMP THREADPRIVATE(XF2I)
REAL, POINTER, DIMENSION(:)    :: XGC=>NULL()
!$OMP THREADPRIVATE(XGC)
REAL, POINTER, DIMENSION(:)    :: XAH=>NULL()
!$OMP THREADPRIVATE(XAH)
REAL, POINTER, DIMENSION(:)    :: XBH=>NULL()
!$OMP THREADPRIVATE(XBH)
REAL, POINTER, DIMENSION(:)    :: XDMAX=>NULL()
!$OMP THREADPRIVATE(XDMAX)
REAL, POINTER, DIMENSION(:)    :: XCE_NITRO=>NULL()
!$OMP THREADPRIVATE(XCE_NITRO)
REAL, POINTER, DIMENSION(:)    :: XCF_NITRO=>NULL()
!$OMP THREADPRIVATE(XCF_NITRO)
REAL, POINTER, DIMENSION(:)    :: XCNA_NITRO=>NULL()
!$OMP THREADPRIVATE(XCNA_NITRO)
REAL, POINTER, DIMENSION(:)    :: XBSLAI_NITRO=>NULL()
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
REAL, POINTER, DIMENSION(:,:)  :: XDG=>NULL()
!$OMP THREADPRIVATE(XDG)
REAL, POINTER, DIMENSION(:)    :: XRUNOFFD=>NULL()
!$OMP THREADPRIVATE(XRUNOFFD)
REAL, POINTER, DIMENSION(:,:)  :: XSOILWGHT=>NULL()
!$OMP THREADPRIVATE(XSOILWGHT)
REAL, POINTER, DIMENSION(:,:)  :: XDZG=>NULL()
!$OMP THREADPRIVATE(XDZG)
REAL, POINTER, DIMENSION(:,:)  :: XDZDIF=>NULL()
!$OMP THREADPRIVATE(XDZDIF)
INTEGER, POINTER, DIMENSION(:) :: NWG_LAYER=>NULL()
!$OMP THREADPRIVATE(NWG_LAYER)
REAL, POINTER, DIMENSION(:)    :: XDROOT=>NULL()
!$OMP THREADPRIVATE(XDROOT)
REAL, POINTER, DIMENSION(:)    :: XDG2=>NULL()
!$OMP THREADPRIVATE(XDG2)
REAL, POINTER, DIMENSION(:)    :: XC1SAT=>NULL()
!$OMP THREADPRIVATE(XC1SAT)
REAL, POINTER, DIMENSION(:)    :: XC2REF=>NULL()
!$OMP THREADPRIVATE(XC2REF)
REAL, POINTER, DIMENSION(:,:)  :: XC3=>NULL()
!$OMP THREADPRIVATE(XC3)
REAL, POINTER, DIMENSION(:)      :: XC4B=>NULL()
!$OMP THREADPRIVATE(XC4B)
REAL, POINTER, DIMENSION(:)    :: XC4REF=>NULL()
!$OMP THREADPRIVATE(XC4REF)
REAL, POINTER, DIMENSION(:)      :: XACOEF=>NULL()
!$OMP THREADPRIVATE(XACOEF)
REAL, POINTER, DIMENSION(:)      :: XPCOEF=>NULL()
!$OMP THREADPRIVATE(XPCOEF)
REAL, POINTER, DIMENSION(:,:)  :: XWFC=>NULL()
!$OMP THREADPRIVATE(XWFC)
REAL, POINTER, DIMENSION(:,:)  :: XWWILT=>NULL()
!$OMP THREADPRIVATE(XWWILT)
REAL, POINTER, DIMENSION(:,:)  :: XWSAT=>NULL()
!$OMP THREADPRIVATE(XWSAT)
REAL, POINTER, DIMENSION(:,:)  :: XBCOEF=>NULL()
!$OMP THREADPRIVATE(XBCOEF)
REAL, POINTER, DIMENSION(:,:)  :: XCONDSAT=>NULL()
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
REAL, POINTER, DIMENSION(:)     :: XPCPS=>NULL()
!$OMP THREADPRIVATE(XPCPS)
REAL, POINTER, DIMENSION(:)     :: XPLVTT=>NULL()
!$OMP THREADPRIVATE(XPLVTT)
REAL, POINTER, DIMENSION(:)     :: XPLSTT=>NULL()
!$OMP THREADPRIVATE(XPLSTT)
!
!SGH scheme
!
REAL, POINTER, DIMENSION(:)  :: XD_ICE=>NULL()
!$OMP THREADPRIVATE(XD_ICE)
REAL, POINTER, DIMENSION(:)  :: XKSAT_ICE=>NULL()
!$OMP THREADPRIVATE(XKSAT_ICE)
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

SUBROUTINE TEB_GARDEN_PGD_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Save current state for allocated arrays
IF (LKFROM) THEN
TEB_GARDEN_PGD_MODEL(KFROM)%XVEGTYPE=>XVEGTYPE
TEB_GARDEN_PGD_MODEL(KFROM)%XALBNIR_DRY=>XALBNIR_DRY
TEB_GARDEN_PGD_MODEL(KFROM)%XALBVIS_DRY=>XALBVIS_DRY
TEB_GARDEN_PGD_MODEL(KFROM)%XALBUV_DRY=>XALBUV_DRY
TEB_GARDEN_PGD_MODEL(KFROM)%XALBNIR_WET=>XALBNIR_WET
TEB_GARDEN_PGD_MODEL(KFROM)%XALBVIS_WET=>XALBVIS_WET
TEB_GARDEN_PGD_MODEL(KFROM)%XALBUV_WET=>XALBUV_WET
TEB_GARDEN_PGD_MODEL(KFROM)%XALBNIR_SOIL=>XALBNIR_SOIL
TEB_GARDEN_PGD_MODEL(KFROM)%XALBVIS_SOIL=>XALBVIS_SOIL
TEB_GARDEN_PGD_MODEL(KFROM)%XALBUV_SOIL=>XALBUV_SOIL
TEB_GARDEN_PGD_MODEL(KFROM)%XZ0_O_Z0H=>XZ0_O_Z0H
TEB_GARDEN_PGD_MODEL(KFROM)%XALBNIR_VEG=>XALBNIR_VEG
TEB_GARDEN_PGD_MODEL(KFROM)%XALBVIS_VEG=>XALBVIS_VEG
TEB_GARDEN_PGD_MODEL(KFROM)%XALBUV_VEG=>XALBUV_VEG
TEB_GARDEN_PGD_MODEL(KFROM)%XWRMAX_CF=>XWRMAX_CF
TEB_GARDEN_PGD_MODEL(KFROM)%XRSMIN=>XRSMIN
TEB_GARDEN_PGD_MODEL(KFROM)%XGAMMA=>XGAMMA
TEB_GARDEN_PGD_MODEL(KFROM)%XCV=>XCV
TEB_GARDEN_PGD_MODEL(KFROM)%XRGL=>XRGL
TEB_GARDEN_PGD_MODEL(KFROM)%XROOTFRAC=>XROOTFRAC
TEB_GARDEN_PGD_MODEL(KFROM)%XABC=>XABC
TEB_GARDEN_PGD_MODEL(KFROM)%XPOI=>XPOI
TEB_GARDEN_PGD_MODEL(KFROM)%XBSLAI=>XBSLAI
TEB_GARDEN_PGD_MODEL(KFROM)%XLAIMIN=>XLAIMIN
TEB_GARDEN_PGD_MODEL(KFROM)%XSEFOLD=>XSEFOLD
TEB_GARDEN_PGD_MODEL(KFROM)%XH_TREE=>XH_TREE
TEB_GARDEN_PGD_MODEL(KFROM)%XANF=>XANF
TEB_GARDEN_PGD_MODEL(KFROM)%XANMAX=>XANMAX
TEB_GARDEN_PGD_MODEL(KFROM)%XFZERO=>XFZERO
TEB_GARDEN_PGD_MODEL(KFROM)%XEPSO=>XEPSO
TEB_GARDEN_PGD_MODEL(KFROM)%XGAMM=>XGAMM
TEB_GARDEN_PGD_MODEL(KFROM)%XQDGAMM=>XQDGAMM
TEB_GARDEN_PGD_MODEL(KFROM)%XGMES=>XGMES
TEB_GARDEN_PGD_MODEL(KFROM)%XRE25=>XRE25
TEB_GARDEN_PGD_MODEL(KFROM)%XQDGMES=>XQDGMES
TEB_GARDEN_PGD_MODEL(KFROM)%XT1GMES=>XT1GMES
TEB_GARDEN_PGD_MODEL(KFROM)%XT2GMES=>XT2GMES
TEB_GARDEN_PGD_MODEL(KFROM)%XAMAX=>XAMAX
TEB_GARDEN_PGD_MODEL(KFROM)%XQDAMAX=>XQDAMAX
TEB_GARDEN_PGD_MODEL(KFROM)%XT1AMAX=>XT1AMAX
TEB_GARDEN_PGD_MODEL(KFROM)%XT2AMAX=>XT2AMAX
TEB_GARDEN_PGD_MODEL(KFROM)%LSTRESS=>LSTRESS
TEB_GARDEN_PGD_MODEL(KFROM)%XF2I=>XF2I
TEB_GARDEN_PGD_MODEL(KFROM)%XGC=>XGC
TEB_GARDEN_PGD_MODEL(KFROM)%XAH=>XAH
TEB_GARDEN_PGD_MODEL(KFROM)%XBH=>XBH
TEB_GARDEN_PGD_MODEL(KFROM)%XDMAX=>XDMAX
TEB_GARDEN_PGD_MODEL(KFROM)%XCE_NITRO=>XCE_NITRO
TEB_GARDEN_PGD_MODEL(KFROM)%XCF_NITRO=>XCF_NITRO
TEB_GARDEN_PGD_MODEL(KFROM)%XCNA_NITRO=>XCNA_NITRO
TEB_GARDEN_PGD_MODEL(KFROM)%XBSLAI_NITRO=>XBSLAI_NITRO
TEB_GARDEN_PGD_MODEL(KFROM)%XSAND=>XSAND
TEB_GARDEN_PGD_MODEL(KFROM)%XCLAY=>XCLAY
TEB_GARDEN_PGD_MODEL(KFROM)%XRUNOFFB=>XRUNOFFB
TEB_GARDEN_PGD_MODEL(KFROM)%XWDRAIN=>XWDRAIN
TEB_GARDEN_PGD_MODEL(KFROM)%XTAUICE=>XTAUICE
TEB_GARDEN_PGD_MODEL(KFROM)%XGAMMAT=>XGAMMAT
TEB_GARDEN_PGD_MODEL(KFROM)%XDG=>XDG
TEB_GARDEN_PGD_MODEL(KFROM)%XRUNOFFD=>XRUNOFFD
TEB_GARDEN_PGD_MODEL(KFROM)%XSOILWGHT=>XSOILWGHT
TEB_GARDEN_PGD_MODEL(KFROM)%XDZG=>XDZG
TEB_GARDEN_PGD_MODEL(KFROM)%XDZDIF=>XDZDIF
TEB_GARDEN_PGD_MODEL(KFROM)%NWG_LAYER=>NWG_LAYER
TEB_GARDEN_PGD_MODEL(KFROM)%XDROOT=>XDROOT
TEB_GARDEN_PGD_MODEL(KFROM)%XDG2=>XDG2
TEB_GARDEN_PGD_MODEL(KFROM)%XPCPS=>XPCPS
TEB_GARDEN_PGD_MODEL(KFROM)%XPLVTT=>XPLVTT
TEB_GARDEN_PGD_MODEL(KFROM)%XPLSTT=>XPLSTT
TEB_GARDEN_PGD_MODEL(KFROM)%XC1SAT=>XC1SAT
TEB_GARDEN_PGD_MODEL(KFROM)%XC2REF=>XC2REF
TEB_GARDEN_PGD_MODEL(KFROM)%XC3=>XC3
TEB_GARDEN_PGD_MODEL(KFROM)%XC4B=>XC4B
TEB_GARDEN_PGD_MODEL(KFROM)%XC4REF=>XC4REF
TEB_GARDEN_PGD_MODEL(KFROM)%XACOEF=>XACOEF
TEB_GARDEN_PGD_MODEL(KFROM)%XPCOEF=>XPCOEF
TEB_GARDEN_PGD_MODEL(KFROM)%XWFC=>XWFC
TEB_GARDEN_PGD_MODEL(KFROM)%XWWILT=>XWWILT
TEB_GARDEN_PGD_MODEL(KFROM)%XWSAT=>XWSAT
TEB_GARDEN_PGD_MODEL(KFROM)%XBCOEF=>XBCOEF
TEB_GARDEN_PGD_MODEL(KFROM)%XCONDSAT=>XCONDSAT
TEB_GARDEN_PGD_MODEL(KFROM)%XMPOTSAT=>XMPOTSAT
TEB_GARDEN_PGD_MODEL(KFROM)%XCGSAT=>XCGSAT
TEB_GARDEN_PGD_MODEL(KFROM)%XHCAPSOIL=>XHCAPSOIL
TEB_GARDEN_PGD_MODEL(KFROM)%XCONDDRY=>XCONDDRY
TEB_GARDEN_PGD_MODEL(KFROM)%XCONDSLD=>XCONDSLD
TEB_GARDEN_PGD_MODEL(KFROM)%XTDEEP=>XTDEEP
!
!
!SGH scheme
!                    
TEB_GARDEN_PGD_MODEL(KFROM)%XD_ICE=>XD_ICE
TEB_GARDEN_PGD_MODEL(KFROM)%XKSAT_ICE=>XKSAT_ICE
ENDIF
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_TEB_GARDEN_PGD_N:TEB_GARDEN_PGD_GOTO_MODEL',0,ZHOOK_HANDLE)
XVEGTYPE=>TEB_GARDEN_PGD_MODEL(KTO)%XVEGTYPE
XALBNIR_DRY=>TEB_GARDEN_PGD_MODEL(KTO)%XALBNIR_DRY
XALBVIS_DRY=>TEB_GARDEN_PGD_MODEL(KTO)%XALBVIS_DRY
XALBUV_DRY=>TEB_GARDEN_PGD_MODEL(KTO)%XALBUV_DRY
XALBNIR_WET=>TEB_GARDEN_PGD_MODEL(KTO)%XALBNIR_WET
XALBVIS_WET=>TEB_GARDEN_PGD_MODEL(KTO)%XALBVIS_WET
XALBUV_WET=>TEB_GARDEN_PGD_MODEL(KTO)%XALBUV_WET
XALBNIR_SOIL=>TEB_GARDEN_PGD_MODEL(KTO)%XALBNIR_SOIL
XALBVIS_SOIL=>TEB_GARDEN_PGD_MODEL(KTO)%XALBVIS_SOIL
XALBUV_SOIL=>TEB_GARDEN_PGD_MODEL(KTO)%XALBUV_SOIL
XZ0_O_Z0H=>TEB_GARDEN_PGD_MODEL(KTO)%XZ0_O_Z0H
XALBNIR_VEG=>TEB_GARDEN_PGD_MODEL(KTO)%XALBNIR_VEG
XALBVIS_VEG=>TEB_GARDEN_PGD_MODEL(KTO)%XALBVIS_VEG
XALBUV_VEG=>TEB_GARDEN_PGD_MODEL(KTO)%XALBUV_VEG
XWRMAX_CF=>TEB_GARDEN_PGD_MODEL(KTO)%XWRMAX_CF
XRSMIN=>TEB_GARDEN_PGD_MODEL(KTO)%XRSMIN
XGAMMA=>TEB_GARDEN_PGD_MODEL(KTO)%XGAMMA
XCV=>TEB_GARDEN_PGD_MODEL(KTO)%XCV
XRGL=>TEB_GARDEN_PGD_MODEL(KTO)%XRGL
XROOTFRAC=>TEB_GARDEN_PGD_MODEL(KTO)%XROOTFRAC
XABC=>TEB_GARDEN_PGD_MODEL(KTO)%XABC
XPOI=>TEB_GARDEN_PGD_MODEL(KTO)%XPOI
XBSLAI=>TEB_GARDEN_PGD_MODEL(KTO)%XBSLAI
XLAIMIN=>TEB_GARDEN_PGD_MODEL(KTO)%XLAIMIN
XSEFOLD=>TEB_GARDEN_PGD_MODEL(KTO)%XSEFOLD
XH_TREE=>TEB_GARDEN_PGD_MODEL(KTO)%XH_TREE
XANF=>TEB_GARDEN_PGD_MODEL(KTO)%XANF
XANMAX=>TEB_GARDEN_PGD_MODEL(KTO)%XANMAX
XFZERO=>TEB_GARDEN_PGD_MODEL(KTO)%XFZERO
XEPSO=>TEB_GARDEN_PGD_MODEL(KTO)%XEPSO
XGAMM=>TEB_GARDEN_PGD_MODEL(KTO)%XGAMM
XQDGAMM=>TEB_GARDEN_PGD_MODEL(KTO)%XQDGAMM
XGMES=>TEB_GARDEN_PGD_MODEL(KTO)%XGMES
XRE25=>TEB_GARDEN_PGD_MODEL(KTO)%XRE25
XQDGMES=>TEB_GARDEN_PGD_MODEL(KTO)%XQDGMES
XT1GMES=>TEB_GARDEN_PGD_MODEL(KTO)%XT1GMES
XT2GMES=>TEB_GARDEN_PGD_MODEL(KTO)%XT2GMES
XAMAX=>TEB_GARDEN_PGD_MODEL(KTO)%XAMAX
XQDAMAX=>TEB_GARDEN_PGD_MODEL(KTO)%XQDAMAX
XT1AMAX=>TEB_GARDEN_PGD_MODEL(KTO)%XT1AMAX
XT2AMAX=>TEB_GARDEN_PGD_MODEL(KTO)%XT2AMAX
LSTRESS=>TEB_GARDEN_PGD_MODEL(KTO)%LSTRESS
XF2I=>TEB_GARDEN_PGD_MODEL(KTO)%XF2I
XGC=>TEB_GARDEN_PGD_MODEL(KTO)%XGC
XAH=>TEB_GARDEN_PGD_MODEL(KTO)%XAH
XBH=>TEB_GARDEN_PGD_MODEL(KTO)%XBH
XDMAX=>TEB_GARDEN_PGD_MODEL(KTO)%XDMAX
XCE_NITRO=>TEB_GARDEN_PGD_MODEL(KTO)%XCE_NITRO
XCF_NITRO=>TEB_GARDEN_PGD_MODEL(KTO)%XCF_NITRO
XCNA_NITRO=>TEB_GARDEN_PGD_MODEL(KTO)%XCNA_NITRO
XBSLAI_NITRO=>TEB_GARDEN_PGD_MODEL(KTO)%XBSLAI_NITRO
XSAND=>TEB_GARDEN_PGD_MODEL(KTO)%XSAND
XCLAY=>TEB_GARDEN_PGD_MODEL(KTO)%XCLAY
XRUNOFFB=>TEB_GARDEN_PGD_MODEL(KTO)%XRUNOFFB
XWDRAIN=>TEB_GARDEN_PGD_MODEL(KTO)%XWDRAIN
XTAUICE=>TEB_GARDEN_PGD_MODEL(KTO)%XTAUICE
XGAMMAT=>TEB_GARDEN_PGD_MODEL(KTO)%XGAMMAT
XDG=>TEB_GARDEN_PGD_MODEL(KTO)%XDG
XRUNOFFD=>TEB_GARDEN_PGD_MODEL(KTO)%XRUNOFFD
XSOILWGHT=>TEB_GARDEN_PGD_MODEL(KTO)%XSOILWGHT
XDZG=>TEB_GARDEN_PGD_MODEL(KTO)%XDZG
XDZDIF=>TEB_GARDEN_PGD_MODEL(KTO)%XDZDIF
NWG_LAYER=>TEB_GARDEN_PGD_MODEL(KTO)%NWG_LAYER
XDROOT=>TEB_GARDEN_PGD_MODEL(KTO)%XDROOT
XDG2=>TEB_GARDEN_PGD_MODEL(KTO)%XDG2
XPCPS=>TEB_GARDEN_PGD_MODEL(KTO)%XPCPS
XPLVTT=>TEB_GARDEN_PGD_MODEL(KTO)%XPLVTT
XPLSTT=>TEB_GARDEN_PGD_MODEL(KTO)%XPLSTT
XC1SAT=>TEB_GARDEN_PGD_MODEL(KTO)%XC1SAT
XC2REF=>TEB_GARDEN_PGD_MODEL(KTO)%XC2REF
XC3=>TEB_GARDEN_PGD_MODEL(KTO)%XC3
XC4B=>TEB_GARDEN_PGD_MODEL(KTO)%XC4B
XC4REF=>TEB_GARDEN_PGD_MODEL(KTO)%XC4REF
XACOEF=>TEB_GARDEN_PGD_MODEL(KTO)%XACOEF
XPCOEF=>TEB_GARDEN_PGD_MODEL(KTO)%XPCOEF
XWFC=>TEB_GARDEN_PGD_MODEL(KTO)%XWFC
XWWILT=>TEB_GARDEN_PGD_MODEL(KTO)%XWWILT
XWSAT=>TEB_GARDEN_PGD_MODEL(KTO)%XWSAT
XBCOEF=>TEB_GARDEN_PGD_MODEL(KTO)%XBCOEF
XCONDSAT=>TEB_GARDEN_PGD_MODEL(KTO)%XCONDSAT
XMPOTSAT=>TEB_GARDEN_PGD_MODEL(KTO)%XMPOTSAT
XCGSAT=>TEB_GARDEN_PGD_MODEL(KTO)%XCGSAT
XHCAPSOIL=>TEB_GARDEN_PGD_MODEL(KTO)%XHCAPSOIL
XCONDDRY=>TEB_GARDEN_PGD_MODEL(KTO)%XCONDDRY
XCONDSLD=>TEB_GARDEN_PGD_MODEL(KTO)%XCONDSLD
XTDEEP=>TEB_GARDEN_PGD_MODEL(KTO)%XTDEEP
!
!SGH scheme
!
XD_ICE=>TEB_GARDEN_PGD_MODEL(KTO)%XD_ICE
XKSAT_ICE=>TEB_GARDEN_PGD_MODEL(KTO)%XKSAT_ICE
!
! Types of vegetation (simplification of veg characterization)
!
CTYPE_HVEG=>TEB_GARDEN_PGD_MODEL(KTO)%CTYPE_HVEG
CTYPE_LVEG=>TEB_GARDEN_PGD_MODEL(KTO)%CTYPE_LVEG
CTYPE_NVEG=>TEB_GARDEN_PGD_MODEL(KTO)%CTYPE_NVEG
IF (LHOOK) CALL DR_HOOK('MODD_TEB_GARDEN_PGD_N:TEB_GARDEN_PGD_GOTO_MODEL',1,ZHOOK_HANDLE)
!
END SUBROUTINE TEB_GARDEN_PGD_GOTO_MODEL

SUBROUTINE TEB_GARDEN_PGD_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GARDEN_PGD_N:TEB_GARDEN_PGD_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(TEB_GARDEN_PGD_MODEL(KMODEL))
DO J=1,KMODEL
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XVEGTYPE)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XALBNIR_DRY)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XALBVIS_DRY)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XALBUV_DRY)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XALBNIR_WET)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XALBVIS_WET)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XALBUV_WET)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XALBNIR_SOIL)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XALBVIS_SOIL)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XALBUV_SOIL)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XZ0_O_Z0H)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XALBNIR_VEG)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XALBVIS_VEG)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XALBUV_VEG)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XWRMAX_CF)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XRSMIN)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XGAMMA)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XCV)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XRGL)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XROOTFRAC)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XBSLAI)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XLAIMIN)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XSEFOLD)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XH_TREE)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XANF)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XANMAX)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XFZERO)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XEPSO)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XGAMM)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XQDGAMM)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XGMES)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XRE25)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XQDGMES)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XT1GMES)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XT2GMES)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XAMAX)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XQDAMAX)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XT1AMAX)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XT2AMAX)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%LSTRESS)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XF2I)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XGC)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XAH)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XBH)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XDMAX)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XCE_NITRO)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XCF_NITRO)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XCNA_NITRO)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XBSLAI_NITRO)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XSAND)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XCLAY)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XRUNOFFB)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XWDRAIN)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XTAUICE)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XGAMMAT)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XDG)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XRUNOFFD)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XSOILWGHT)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XDZG)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XDZDIF)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%NWG_LAYER)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XDROOT)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XDG2)  
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XPCPS)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XPLVTT)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XPLSTT)  
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XC1SAT)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XC2REF)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XC3)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XC4B)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XC4REF)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XACOEF)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XPCOEF)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XWFC)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XWWILT)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XWSAT)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XBCOEF)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XCONDSAT)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XMPOTSAT)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XCGSAT)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XHCAPSOIL)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XCONDDRY)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XCONDSLD)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XTDEEP)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XD_ICE)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XKSAT_ICE)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XABC)
  NULLIFY(TEB_GARDEN_PGD_MODEL(J)%XPOI)
ENDDO
TEB_GARDEN_PGD_MODEL(:)%CTYPE_HVEG=' '
TEB_GARDEN_PGD_MODEL(:)%CTYPE_LVEG=' '
TEB_GARDEN_PGD_MODEL(:)%CTYPE_NVEG=' '
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GARDEN_PGD_N:TEB_GARDEN_PGD_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_GARDEN_PGD_ALLOC

SUBROUTINE TEB_GARDEN_PGD_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GARDEN_PGD_N:TEB_GARDEN_PGD_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(TEB_GARDEN_PGD_MODEL)) DEALLOCATE(TEB_GARDEN_PGD_MODEL)
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GARDEN_PGD_N:TEB_GARDEN_PGD_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_GARDEN_PGD_DEALLO

END MODULE MODD_TEB_GARDEN_PGD_n
