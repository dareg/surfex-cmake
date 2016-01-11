!##################
MODULE MODD_ISBA_INIT_n
!##################
!
!!****  *MODD_ISBA - declaration of packed surface parameters for ISBA scheme
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
!!      P. Samuelsson   02/2012 : MEB
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
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
TYPE ISBA_INIT_PGD_t
!
! Mask and number of grid elements containing patches/tiles:
!
REAL, POINTER, DIMENSION(:,:)    :: XPATCH         ! fraction of each tile/patch             (-)
REAL, POINTER, DIMENSION(:,:,:)  :: XVEGTYPE_PATCH ! fraction of each vegetation type for
INTEGER, POINTER, DIMENSION(:)   :: NSIZE_NATURE_P ! number of sub-patchs/tiles              (-)
INTEGER, POINTER, DIMENSION(:,:) :: NR_NATURE_P    ! patch/tile mask                         (-)
!
REAL, POINTER, DIMENSION(:,:)    :: XPATCH_OLD     ! fraction of each tile/patchfor land use (-)
!
! * DEEPSOIL = T
!
REAL, POINTER, DIMENSION(:)      :: XTDEEP         ! prescribed deep soil temperature 
!                                                  ! (optional)
REAL, POINTER, DIMENSION(:)      :: XGAMMAT        ! 'Force-Restore' timescale when using a
!                                                  ! prescribed lower boundary temperature   (1/days)
! * HPHOTO /= NON , co2_init_n
!
REAL, POINTER, DIMENSION(:)      :: XABC           ! abscissa needed for integration
!                                                  ! of net assimilation and stomatal
!                                                  ! conductance over canopy depth           (-)
REAL, POINTER, DIMENSION(:)      :: XPOI           ! Gaussian weights for integration
!                                                  ! of net assimilation and stomatal
!                                                  ! conductance over canopy depth           (-)
REAL, POINTER, DIMENSION(:,:)    :: XANMAX         ! maximum photosynthesis rate             (
REAL, POINTER, DIMENSION(:,:)    :: XFZERO         ! ideal value of F, no photo- 
!                                                  ! respiration or saturation deficit  
REAL, POINTER, DIMENSION(:,:)    :: XEPSO          ! maximum initial quantum use             
!                                                  ! efficiency                              (mg J-1 PAR)
REAL, POINTER, DIMENSION(:,:)    :: XGAMM          ! CO2 conpensation concentration          (ppm)
REAL, POINTER, DIMENSION(:,:)    :: XQDGAMM        ! Log of Q10 function for CO2 conpensation 
!                                                 ! concentration                           (-)
REAL, POINTER, DIMENSION(:,:)    :: XQDGMES        ! Log of Q10 function for mesophyll conductance  (-)
REAL, POINTER, DIMENSION(:,:)    :: XT1GMES        ! reference temperature for computing 
!                                                  ! compensation concentration function for 
!                                                  ! mesophyll conductance: minimum
!                                                  ! temperature                             (K)
REAL, POINTER, DIMENSION(:,:)    :: XT2GMES        ! reference temperature for computing 
!                                                  ! compensation concentration function for 
!                                                  ! mesophyll conductance: maximum
!                                                  ! temperature                             (K)
REAL, POINTER, DIMENSION(:,:)    :: XAMAX          ! leaf photosynthetic capacity            (mg m-2 s-1)
REAL, POINTER, DIMENSION(:,:)    :: XQDAMAX        ! Log of Q10 function for leaf photosynthetic 
!                                                  ! capacity                                (-)
REAL, POINTER, DIMENSION(:,:)    :: XT1AMAX        ! reference temperature for computing 
!                                                  ! compensation concentration function for 
!                                                  ! leaf photosynthetic capacity: minimum
!                                                  ! temperature                             (K)
REAL, POINTER, DIMENSION(:,:)    :: XT2AMAX        ! reference temperature for computing 
!                                                  ! compensation concentration function for 
!                                                  ! leaf photosynthetic capacity: maximum
!                                                  ! temperature                             (K)
REAL, POINTER, DIMENSION(:,:)    :: XAH            ! coefficients for herbaceous water stress 
!                                                  ! response (offensive or defensive)       (log(mm/s))
REAL, POINTER, DIMENSION(:,:)    :: XBH            ! coefficients for herbaceous water stress 
!                                                  ! response (offensive or defensive)
REAL, POINTER, DIMENSION(:,:)    :: XTAU_WOOD      ! residence time in woody biomass         (s)
REAL, POINTER, DIMENSION(:,:,:)   :: XINCREASE     ! biomass increase                     (kg/m2/day)
REAL, POINTER, DIMENSION(:,:,:)   :: XTURNOVER     ! turnover rates from biomass to litter (gC/m2/s)
!
! directional total roughness lenghts in 4 coordinate directions
! (IP: i index up;  IM: i index down;  JP: j index up;  JM: j index down)
!
REAL, DIMENSION(:,:), POINTER :: XZ0EFFIP,XZ0EFFIM,XZ0EFFJP,XZ0EFFJM
REAL, DIMENSION(:), POINTER   :: XZ0REL         ! relief roughness length                 (m)
!
! *Soil hydraulic characteristics
!
REAL, POINTER, DIMENSION(:,:,:)  :: XCONDSAT       ! hydraulic conductivity at saturation    (m/s)
REAL, POINTER, DIMENSION(:,:)    :: XMPOTSAT       ! matric potential at saturation          (m)
REAL, POINTER, DIMENSION(:,:)    :: XBCOEF         ! soil water CH78 b-parameter             (-)
REAL, POINTER, DIMENSION(:,:)    :: XWWILT         ! wilting point volumetric water content 
!                                                  ! profile                                 (m3/m3)
REAL, POINTER, DIMENSION(:,:)    :: XWFC           ! field capacity volumetric water content
!                                                  ! profile                                 (m3/m3)
REAL, POINTER, DIMENSION(:,:)    :: XWSAT          ! porosity profile                        (m3/m3) 
REAL, POINTER, DIMENSION(:)      :: XTAUICE        ! soil freezing characteristic timescale  (s)
!
REAL, POINTER, DIMENSION(:)      :: XCGSAT         ! soil thermal inertia coefficient at 
!                                                  ! saturation   
REAL, POINTER, DIMENSION(:,:)    :: XC1SAT         ! 'Force-Restore' C1 coefficient at 
!                                                  ! saturation                              (-)
REAL, POINTER, DIMENSION(:,:)    :: XC2REF         ! 'Force-Restore' reference value of C2   (-)
REAL, POINTER, DIMENSION(:,:,:)  :: XC3            ! 'Force-Restore' C3 drainage coefficient (m)
REAL, POINTER, DIMENSION(:)      :: XC4B           ! 'Force-Restore' sub-surface vertical 
!                                                  ! diffusion coefficient (slope parameter) (-)
REAL, POINTER, DIMENSION(:)      :: XACOEF         ! 'Force-Restore' surface vertical 
!                                                  ! diffusion coefficient                   (-)
REAL, POINTER, DIMENSION(:)      :: XPCOEF         ! 'Force-Restore' surface vertical 
!                                                  ! diffusion coefficient                   (-)
REAL, POINTER, DIMENSION(:,:)    :: XC4REF         ! 'Force-Restore' sub-surface vertical 
!                                                  ! for lateral drainage ('DIF' option)
!
! * SGH initializations
!
REAL, POINTER, DIMENSION(:,:)  :: XWD0     ! water content equivalent to TOPMODEL maximum deficit
REAL, POINTER, DIMENSION(:,:)  :: XKANISO  ! Anisotropy coeficient for hydraulic conductivity 
!
! * Soil thermal characteristics
!
REAL, POINTER, DIMENSION(:,:)    :: XPCPS
REAL, POINTER, DIMENSION(:,:)    :: XPLVTT
REAL, POINTER, DIMENSION(:,:)    :: XPLSTT 
!
REAL, POINTER, DIMENSION(:,:)    :: XHCAPSOIL      ! soil heat capacity                      (J/K/m3)
REAL, POINTER, DIMENSION(:,:)    :: XCONDDRY       ! soil dry thermal conductivity           (W/m/K)
REAL, POINTER, DIMENSION(:,:)    :: XCONDSLD       ! soil solids thermal conductivity        (W/m/K)
!
! * Initialize hydrology
!
  REAL, POINTER, DIMENSION(:,:)    :: XRUNOFFD     ! depth over which sub-grid runoff is
!                                                  ! computed: in Force-Restore this is the
!                                                  ! total soil column ('2-L'), or root zone
!                                                  ! ('3-L'). For the 'DIF' option, it can
!                                                  ! be any depth within soil column         (m)
!
REAL, POINTER, DIMENSION(:,:,:)  :: XDZG           ! soil layers thicknesses (DIF option)
REAL, POINTER, DIMENSION(:,:,:)  :: XDZDIF         ! distance between consecuative layer mid-points (DIF option)
REAL, POINTER, DIMENSION(:,:,:)  :: XSOILWGHT      ! ISBA-DIF: weights for vertical
!                                                  ! integration of soil water and properties
! - Water table depth coupling
!  
REAL, POINTER, DIMENSION(:)  :: XFWTD         ! grid-cell fraction of water table rise
REAL, POINTER, DIMENSION(:)  :: XWTD          ! water table depth (m)
!
REAL, POINTER, DIMENSION(:,:)  :: XKSAT_ICE   ! hydraulic conductivity at saturation
!                                               over frozen area (m s-1)
!
! * Physiographic radiative fields
!
REAL, POINTER, DIMENSION(:)   :: XALBNIR_DRY       ! dry soil near-infra-red albedo          (-)
REAL, POINTER, DIMENSION(:)   :: XALBVIS_DRY       ! dry soil visible albedo                 (-)
REAL, POINTER, DIMENSION(:)   :: XALBUV_DRY        ! dry soil UV albedo                      (-)
REAL, POINTER, DIMENSION(:)   :: XALBNIR_WET       ! wet soil near-infra-red albedo          (-)
REAL, POINTER, DIMENSION(:)   :: XALBVIS_WET       ! wet soil visible albedo                 (-)
REAL, POINTER, DIMENSION(:)   :: XALBUV_WET        ! wet soil UV albedo                      (-)
!
REAL, POINTER, DIMENSION(:,:)    :: XBSLAI_NITRO   ! biomass/LAI ratio from nitrogen 
!                                                  ! decline theory                        (kg/m2)
!
END TYPE ISBA_INIT_PGD_t
!
TYPE ISBA_INIT_PGD_PATCH_t
!
TYPE(ISBA_INIT_PGD_t), ALLOCATABLE :: AL(:) 
!
END TYPE ISBA_INIT_PGD_PATCH_t
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
TYPE ISBA_INIT_t
!
TYPE (DATE_TIME) :: TTIME            ! current date and time
!
REAL, POINTER, DIMENSION(:,:) :: XFRACSOC ! Fraction of organic carbon in each soil layer
!
REAL, POINTER, DIMENSION(:,:) :: XTAB_FSAT !Satured fraction array
REAL, POINTER, DIMENSION(:,:) :: XTAB_WTOP !Active TOPMODEL-layer array
REAL, POINTER, DIMENSION(:,:) :: XTAB_QTOP !Subsurface flow TOPMODEL array
!
REAL, POINTER, DIMENSION(:,:,:) :: XTOPQS  ! Topmodel subsurface flow by layer (m/s)
!
REAL, POINTER, DIMENSION(:) :: XF_PARAM
REAL, POINTER, DIMENSION(:) :: XC_DEPTH_RATIO
!
REAL, POINTER, DIMENSION(:)    :: XMUF     ! fraction of the grid cell reached by the rainfall
REAL, POINTER, DIMENSION(:)    :: XFSAT    ! Topmodel or dt92 saturated fracti
!
! - Coupling with river routing model
!  
REAL, POINTER, DIMENSION(:)  :: XCPL_DRAIN   ! Surface runoff
REAL, POINTER, DIMENSION(:)  :: XCPL_RUNOFF  ! Deep drainage or gourdwater recharge
REAL, POINTER, DIMENSION(:)  :: XCPL_ICEFLUX ! Calving flux
REAL, POINTER, DIMENSION(:)  :: XCPL_RECHARGE! Groundwater recharge
REAL, POINTER, DIMENSION(:)  :: XCPL_EFLOOD  ! floodplains evaporation
REAL, POINTER, DIMENSION(:)  :: XCPL_PFLOOD  ! floodplains precipitation interception
REAL, POINTER, DIMENSION(:)  :: XCPL_IFLOOD  ! floodplains infiltration
!
! - Flood scheme
!
REAL, POINTER, DIMENSION(:)  :: XFFLOOD      ! Grid-cell flood fraction
REAL, POINTER, DIMENSION(:)  :: XPIFLOOD     ! flood potential infiltration (kg/m2)
REAL, POINTER, DIMENSION(:,:) :: XFF         ! Total Flood fraction  
REAL, POINTER, DIMENSION(:,:) :: XFFG        ! Flood fraction over ground
REAL, POINTER, DIMENSION(:,:) :: XFFV        ! Flood fraction over vegetation
REAL, POINTER, DIMENSION(:,:) :: XFFROZEN    ! Fraction of frozen floodplains
REAL, POINTER, DIMENSION(:,:) :: XALBF       ! Flood albedo
REAL, POINTER, DIMENSION(:,:) :: XEMISF      ! Flood emissivity
!
! - Snow and flood fractions and total albedo at time t:                             (-)
!
REAL, POINTER, DIMENSION(:,:,:) :: XDIR_ALB_WITH_SNOW ! total direct albedo by bands
REAL, POINTER, DIMENSION(:,:,:) :: XSCA_ALB_WITH_SNOW ! total diffuse albedo by bands
!
REAL, POINTER, DIMENSION(:)   :: XEMIS_NAT         ! patch averaged emissivity               (-)
!
!  - Random perturbations
!
REAL, POINTER, DIMENSION(:)     :: XPERTVEG
REAL, POINTER, DIMENSION(:)     :: XPERTLAI
REAL, POINTER, DIMENSION(:)     :: XPERTCV
REAL, POINTER, DIMENSION(:)     :: XPERTALB
REAL, POINTER, DIMENSION(:)     :: XPERTZ0
!
END TYPE ISBA_INIT_t
!
TYPE ISBA_INIT_PATCH_t
!
TYPE(ISBA_INIT_t), ALLOCATABLE :: AL(:) 
!
END TYPE ISBA_INIT_PATCH_t
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
CONTAINS
!
SUBROUTINE ISBA_INIT_PGD_PATCH_INIT(YISBA_INIT_PGD_PATCH,KPATCH)
TYPE(ISBA_INIT_PGD_PATCH_t), INTENT(INOUT) :: YISBA_INIT_PGD_PATCH 
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_INIT_N:ISBA_INIT_PGD_PATCH_INIT",0,ZHOOK_HANDLE)
 ALLOCATE(YISBA_INIT_PGD_PATCH%AL(KPATCH))
DO JP=1,KPATCH
  CALL ISBA_INIT_PGD_INIT(YISBA_INIT_PGD_PATCH%AL(JP))
ENDDO         
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_INIT_N:ISBA_INIT_PGD_PATCH_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_INIT_PGD_PATCH_INIT
!
SUBROUTINE ISBA_INIT_PATCH_INIT(YISBA_INIT_PATCH,KPATCH)
TYPE(ISBA_INIT_PATCH_t), INTENT(INOUT) :: YISBA_INIT_PATCH 
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_INIT_N:ISBA_INIT_PATCH_INIT",0,ZHOOK_HANDLE)
 ALLOCATE(YISBA_INIT_PATCH%AL(KPATCH))
DO JP=1,KPATCH
  CALL ISBA_INIT_INIT(YISBA_INIT_PATCH%AL(JP))
ENDDO         
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_INIT_N:ISBA_INIT_PATCH_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_INIT_PATCH_INIT
!
SUBROUTINE ISBA_INIT_PGD_INIT(YISBA_INIT_PGD)
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: YISBA_INIT_PGD
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_INIT_N:ISBA_INIT_PGD_INIT",0,ZHOOK_HANDLE)
!
NULLIFY(YISBA_INIT_PGD%XPATCH)
NULLIFY(YISBA_INIT_PGD%XVEGTYPE_PATCH)
NULLIFY(YISBA_INIT_PGD%NSIZE_NATURE_P)
NULLIFY(YISBA_INIT_PGD%NR_NATURE_P)
NULLIFY(YISBA_INIT_PGD%XPATCH_OLD)  
NULLIFY(YISBA_INIT_PGD%XTDEEP)  
NULLIFY(YISBA_INIT_PGD%XGAMMAT)
NULLIFY(YISBA_INIT_PGD%XABC)
NULLIFY(YISBA_INIT_PGD%XPOI) 
NULLIFY(YISBA_INIT_PGD%XANMAX)
NULLIFY(YISBA_INIT_PGD%XFZERO)
NULLIFY(YISBA_INIT_PGD%XEPSO)
NULLIFY(YISBA_INIT_PGD%XGAMM)
NULLIFY(YISBA_INIT_PGD%XQDGAMM)
NULLIFY(YISBA_INIT_PGD%XQDGMES)
NULLIFY(YISBA_INIT_PGD%XT1GMES)
NULLIFY(YISBA_INIT_PGD%XT2GMES)
NULLIFY(YISBA_INIT_PGD%XAMAX)
NULLIFY(YISBA_INIT_PGD%XQDAMAX)
NULLIFY(YISBA_INIT_PGD%XT1AMAX)
NULLIFY(YISBA_INIT_PGD%XT2AMAX)
NULLIFY(YISBA_INIT_PGD%XAH)
NULLIFY(YISBA_INIT_PGD%XBH)
NULLIFY(YISBA_INIT_PGD%XTAU_WOOD)
NULLIFY(YISBA_INIT_PGD%XINCREASE)
NULLIFY(YISBA_INIT_PGD%XTURNOVER)
NULLIFY(YISBA_INIT_PGD%XZ0EFFIP)
NULLIFY(YISBA_INIT_PGD%XZ0EFFIM)
NULLIFY(YISBA_INIT_PGD%XZ0EFFJP)
NULLIFY(YISBA_INIT_PGD%XZ0EFFJM)  
NULLIFY(YISBA_INIT_PGD%XZ0REL)  
NULLIFY(YISBA_INIT_PGD%XCONDSAT)
NULLIFY(YISBA_INIT_PGD%XMPOTSAT)
NULLIFY(YISBA_INIT_PGD%XBCOEF)
NULLIFY(YISBA_INIT_PGD%XWWILT)
NULLIFY(YISBA_INIT_PGD%XWFC)  
NULLIFY(YISBA_INIT_PGD%XWSAT)
NULLIFY(YISBA_INIT_PGD%XTAUICE)
NULLIFY(YISBA_INIT_PGD%XCGSAT)
NULLIFY(YISBA_INIT_PGD%XC1SAT)
NULLIFY(YISBA_INIT_PGD%XC2REF)
NULLIFY(YISBA_INIT_PGD%XC3)
NULLIFY(YISBA_INIT_PGD%XC4B)
NULLIFY(YISBA_INIT_PGD%XACOEF)
NULLIFY(YISBA_INIT_PGD%XPCOEF)  
NULLIFY(YISBA_INIT_PGD%XC4REF)
NULLIFY(YISBA_INIT_PGD%XWD0)
NULLIFY(YISBA_INIT_PGD%XKANISO)
NULLIFY(YISBA_INIT_PGD%XPCPS)
NULLIFY(YISBA_INIT_PGD%XPLVTT)
NULLIFY(YISBA_INIT_PGD%XPLSTT)
NULLIFY(YISBA_INIT_PGD%XHCAPSOIL)
NULLIFY(YISBA_INIT_PGD%XCONDDRY)
NULLIFY(YISBA_INIT_PGD%XCONDSLD)
NULLIFY(YISBA_INIT_PGD%XRUNOFFD)
NULLIFY(YISBA_INIT_PGD%XDZG)
NULLIFY(YISBA_INIT_PGD%XDZDIF)
NULLIFY(YISBA_INIT_PGD%XSOILWGHT)
NULLIFY(YISBA_INIT_PGD%XFWTD)
NULLIFY(YISBA_INIT_PGD%XWTD)    
NULLIFY(YISBA_INIT_PGD%XKSAT_ICE)
NULLIFY(YISBA_INIT_PGD%XALBNIR_DRY)
NULLIFY(YISBA_INIT_PGD%XALBVIS_DRY)
NULLIFY(YISBA_INIT_PGD%XALBUV_DRY)
NULLIFY(YISBA_INIT_PGD%XALBNIR_WET)
NULLIFY(YISBA_INIT_PGD%XALBVIS_WET)
NULLIFY(YISBA_INIT_PGD%XALBUV_WET)
NULLIFY(YISBA_INIT_PGD%XBSLAI_NITRO)
!
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_INIT_N:ISBA_INIT_PGD_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_INIT_PGD_INIT


SUBROUTINE ISBA_INIT_INIT(YISBA_INIT)
TYPE(ISBA_INIT_t), INTENT(INOUT) :: YISBA_INIT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_INIT_N:ISBA_INIT_INIT",0,ZHOOK_HANDLE)
!
NULLIFY(YISBA_INIT%XFRACSOC)
NULLIFY(YISBA_INIT%XTAB_FSAT)
NULLIFY(YISBA_INIT%XTAB_WTOP)
NULLIFY(YISBA_INIT%XTAB_QTOP)
NULLIFY(YISBA_INIT%XTOPQS)
NULLIFY(YISBA_INIT%XF_PARAM)
NULLIFY(YISBA_INIT%XC_DEPTH_RATIO)
NULLIFY(YISBA_INIT%XMUF)
NULLIFY(YISBA_INIT%XFSAT)
NULLIFY(YISBA_INIT%XCPL_DRAIN)
NULLIFY(YISBA_INIT%XCPL_RUNOFF)
NULLIFY(YISBA_INIT%XCPL_ICEFLUX)
NULLIFY(YISBA_INIT%XCPL_RECHARGE)
NULLIFY(YISBA_INIT%XCPL_EFLOOD)
NULLIFY(YISBA_INIT%XCPL_PFLOOD)
NULLIFY(YISBA_INIT%XCPL_IFLOOD)
NULLIFY(YISBA_INIT%XFFLOOD)
NULLIFY(YISBA_INIT%XPIFLOOD)  
NULLIFY(YISBA_INIT%XFF)
NULLIFY(YISBA_INIT%XFFG)
NULLIFY(YISBA_INIT%XFFV)
NULLIFY(YISBA_INIT%XFFROZEN)
NULLIFY(YISBA_INIT%XALBF)
NULLIFY(YISBA_INIT%XEMISF)
NULLIFY(YISBA_INIT%XDIR_ALB_WITH_SNOW)
NULLIFY(YISBA_INIT%XSCA_ALB_WITH_SNOW)
NULLIFY(YISBA_INIT%XEMIS_NAT)
NULLIFY(YISBA_INIT%XPERTVEG)
NULLIFY(YISBA_INIT%XPERTLAI)
NULLIFY(YISBA_INIT%XPERTCV)
NULLIFY(YISBA_INIT%XPERTALB)
NULLIFY(YISBA_INIT%XPERTZ0)
!
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_INIT_N:ISBA_INIT_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_INIT_INIT

END MODULE MODD_ISBA_INIT_n
