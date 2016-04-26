!##################
MODULE MODD_ISBA_INIT1_n
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
TYPE ISBA_INIT_PGD_1P_t
!
INTEGER, POINTER   :: NSIZE_NATURE_P ! number of sub-patchs/tiles              (-)
!
! Mask and number of grid elements containing patches/tiles:
!
REAL, POINTER, DIMENSION(:,:)  :: XVEGTYPE_PATCH ! fraction of each vegetation type for
!
! * DEEPSOIL = T
!
REAL, POINTER, DIMENSION(:)      :: XTDEEP         ! prescribed deep soil temperature 
!                                                  ! (optional)
REAL, POINTER, DIMENSION(:)      :: XGAMMAT        ! 'Force-Restore' timescale when using a
!                                                  ! prescribed lower boundary temperature   (1/days)
REAL, POINTER, DIMENSION(:,:)    :: XMPOTSAT       ! matric potential at saturation          (m)
REAL, POINTER, DIMENSION(:,:)    :: XBCOEF         ! soil water CH78 b-parameter             (-)
REAL, POINTER, DIMENSION(:,:)    :: XWWILT         ! wilting point volumetric water content 
!                                                  ! profile                                 (m3/m3)
REAL, POINTER, DIMENSION(:,:)    :: XWFC           ! field capacity volumetric water content
!                                                  ! profile                                 (m3/m3)
REAL, POINTER, DIMENSION(:,:)    :: XWSAT          ! porosity profile                        (m3/m3) 
REAL, POINTER, DIMENSION(:)      :: XTAUICE        ! soil freezing characteristic timescale  (s)
REAL, POINTER, DIMENSION(:)      :: XC4B           ! 'Force-Restore' sub-surface vertical 
!                                                  ! diffusion coefficient (slope parameter) (-)
REAL, POINTER, DIMENSION(:)      :: XACOEF         ! 'Force-Restore' surface vertical 
!                                                  ! diffusion coefficient                   (-)
REAL, POINTER, DIMENSION(:)      :: XPCOEF         ! 'Force-Restore' surface vertical 
!                                                  ! diffusion coefficient                   (-)
REAL, POINTER, DIMENSION(:)      :: XCGSAT         ! soil thermal inertia coefficient at 
!                                                  ! saturation   
REAL, POINTER, DIMENSION(:,:)    :: XHCAPSOIL      ! soil heat capacity                      (J/K/m3)
REAL, POINTER, DIMENSION(:,:)    :: XCONDDRY       ! soil dry thermal conductivity           (W/m/K)
REAL, POINTER, DIMENSION(:,:)    :: XCONDSLD       ! soil solids thermal conductivity        (W/m/K)
!
! - Water table depth coupling
!  
REAL, POINTER, DIMENSION(:)  :: XFWTD         ! grid-cell fraction of water table rise
REAL, POINTER, DIMENSION(:)  :: XWTD          ! water table depth (m)
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
!
! * HPHOTO /= NON , co2_init_n
!
REAL, POINTER, DIMENSION(:)      :: XABC           ! abscissa needed for integration
!                                                  ! of net assimilation and stomatal
!                                                  ! conductance over canopy depth           (-)
REAL, POINTER, DIMENSION(:)      :: XPOI           ! Gaussian weights for integration
!                                                  ! of net assimilation and stomatal
!                                                  ! conductance over canopy depth           (-)
!
REAL, POINTER, DIMENSION(:)    :: XANMAX         ! maximum photosynthesis rate             (
REAL, POINTER, DIMENSION(:)    :: XFZERO         ! ideal value of F, no photo- 
!                                                ! respiration or saturation deficit  
REAL, POINTER, DIMENSION(:)    :: XEPSO          ! maximum initial quantum use             
!                                                ! efficiency                              (mg J-1 PAR)
REAL, POINTER, DIMENSION(:)    :: XGAMM          ! CO2 conpensation concentration          (ppm)
REAL, POINTER, DIMENSION(:)    :: XQDGAMM        ! Log of Q10 function for CO2 conpensation 
!                                               ! concentration                           (-)
REAL, POINTER, DIMENSION(:)    :: XQDGMES        ! Log of Q10 function for mesophyll conductance  (-)
REAL, POINTER, DIMENSION(:)    :: XT1GMES        ! reference temperature for computing 
!                                                ! compensation concentration function for 
!                                                ! mesophyll conductance: minimum
!                                                ! temperature                             (K)
REAL, POINTER, DIMENSION(:)    :: XT2GMES        ! reference temperature for computing 
!                                                ! compensation concentration function for 
!                                                ! mesophyll conductance: maximum
!                                                ! temperature                             (K)
REAL, POINTER, DIMENSION(:)    :: XAMAX          ! leaf photosynthetic capacity            (mg m-2 s-1)
REAL, POINTER, DIMENSION(:)    :: XQDAMAX        ! Log of Q10 function for leaf photosynthetic 
!                                                ! capacity                                (-)
REAL, POINTER, DIMENSION(:)    :: XT1AMAX        ! reference temperature for computing 
!                                                ! compensation concentration function for 
!                                                ! leaf photosynthetic capacity: minimum
!                                                ! temperature                             (K)
REAL, POINTER, DIMENSION(:)    :: XT2AMAX        ! reference temperature for computing 
!                                                ! compensation concentration function for 
!                                                ! leaf photosynthetic capacity: maximum
!                                                ! temperature                             (K)
REAL, POINTER, DIMENSION(:)    :: XAH            ! coefficients for herbaceous water stress 
!                                                ! response (offensive or defensive)       (log(mm/s))
REAL, POINTER, DIMENSION(:)    :: XBH            ! coefficients for herbaceous water stress 
!                                                ! response (offensive or defensive)
REAL, POINTER, DIMENSION(:)    :: XTAU_WOOD      ! residence time in woody biomass         (s)
REAL, POINTER, DIMENSION(:,:)   :: XINCREASE     ! biomass increase                     (kg/m2/day)
REAL, POINTER, DIMENSION(:,:)   :: XTURNOVER     ! turnover rates from biomass to litter (gC/m2/s)
!
! *Soil hydraulic characteristics
!
REAL, POINTER, DIMENSION(:,:)  :: XCONDSAT       ! hydraulic conductivity at saturation    (m/s)
!
REAL, POINTER, DIMENSION(:)    :: XC1SAT         ! 'Force-Restore' C1 coefficient at 
!                                                ! saturation                              (-)
REAL, POINTER, DIMENSION(:)    :: XC2REF         ! 'Force-Restore' reference value of C2   (-)
REAL, POINTER, DIMENSION(:,:)  :: XC3            ! 'Force-Restore' C3 drainage coefficient (m)
REAL, POINTER, DIMENSION(:)    :: XC4REF         ! 'Force-Restore' sub-surface vertical 
!                                                ! for lateral drainage ('DIF' option)
!
! * SGH initializations
!
REAL, POINTER, DIMENSION(:,:)  :: XWD0     ! water content equivalent to TOPMODEL maximum deficit
REAL, POINTER, DIMENSION(:,:)  :: XKANISO  ! Anisotropy coeficient for hydraulic conductivity
!
! * Soil thermal characteristics
!
REAL, POINTER, DIMENSION(:)    :: XPCPS
REAL, POINTER, DIMENSION(:)    :: XPLVTT
REAL, POINTER, DIMENSION(:)    :: XPLSTT 
!
! * Initialize hydrology
!
  REAL, POINTER, DIMENSION(:)    :: XRUNOFFD     ! depth over which sub-grid runoff is
!                                                  ! computed: in Force-Restore this is the
!                                                  ! total soil column ('2-L'), or root zone
!                                                  ! ('3-L'). For the 'DIF' option, it can
!                                                  ! be any depth within soil column         (m)
!
REAL, POINTER, DIMENSION(:,:)  :: XDZG           ! soil layers thicknesses (DIF option)
REAL, POINTER, DIMENSION(:,:)  :: XDZDIF         ! distance between consecuative layer mid-points (DIF option)
REAL, POINTER, DIMENSION(:,:)  :: XSOILWGHT      ! ISBA-DIF: weights for vertical
!                                                  ! integration of soil water and properties
!
REAL, POINTER, DIMENSION(:)  :: XKSAT_ICE        ! hydraulic conductivity at saturation
!                                                    over frozen area (m s-1)
!
REAL, POINTER, DIMENSION(:)    :: XBSLAI_NITRO   ! biomass/LAI ratio from nitrogen 
!                                                  ! decline theory                        (kg/m2)
!
END TYPE ISBA_INIT_PGD_1P_t
!
TYPE ISBA_INIT_PGD_NP_t
!
TYPE(ISBA_INIT_PGD_1P_t), ALLOCATABLE :: AL(:) 
!
END TYPE ISBA_INIT_PGD_NP_t
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
TYPE ISBA_INIT_1P_t
!
TYPE (DATE_TIME) :: TTIME
!
REAL, POINTER, DIMENSION(:)    :: XMUF     ! fraction of the grid cell reached by the rainfall
REAL, POINTER, DIMENSION(:)    :: XFSAT    ! Topmodel or dt92 saturated fracti
!
REAL, POINTER, DIMENSION(:)  :: XFFLOOD      ! Grid-cell flood fraction
REAL, POINTER, DIMENSION(:)  :: XPIFLOOD     ! flood potential infiltration (kg/m2)
!
!
! - Snow and flood fractions and total albedo at time t:                             (-)
!
REAL, POINTER, DIMENSION(:,:) :: XDIR_ALB_WITH_SNOW ! total direct albedo by bands
REAL, POINTER, DIMENSION(:,:) :: XSCA_ALB_WITH_SNOW ! total diffuse albedo by bands
!
! - Flood scheme
!
REAL, POINTER, DIMENSION(:) :: XFF         ! Total Flood fraction  
REAL, POINTER, DIMENSION(:) :: XFFG        ! Flood fraction over ground
REAL, POINTER, DIMENSION(:) :: XFFV        ! Flood fraction over vegetation
REAL, POINTER, DIMENSION(:) :: XFFROZEN    ! Fraction of frozen floodplains
REAL, POINTER, DIMENSION(:) :: XALBF       ! Flood albedo
REAL, POINTER, DIMENSION(:) :: XEMISF      ! Flood emissivity
!
REAL, POINTER, DIMENSION(:,:) :: XTOPQS  ! Topmodel subsurface flow by layer (m/s)
!
END TYPE ISBA_INIT_1P_t
!
TYPE ISBA_INIT_NP_t
!
TYPE(ISBA_INIT_1P_t), ALLOCATABLE :: AL(:) 
!
END TYPE ISBA_INIT_NP_t
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
CONTAINS
!
SUBROUTINE ISBA_INIT_PGD_1P_INIT(IP)
TYPE(ISBA_INIT_PGD_1P_t), INTENT(INOUT) :: IP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_INIT1_N:ISBA_INIT_PGD_1P_INIT",0,ZHOOK_HANDLE)
  NULLIFY(IP%XVEGTYPE_PATCH)
  NULLIFY(IP%XTDEEP)  
  NULLIFY(IP%XGAMMAT)
  NULLIFY(IP%XABC)
  NULLIFY(IP%XPOI)   
  NULLIFY(IP%XANMAX)
  NULLIFY(IP%XFZERO)
  NULLIFY(IP%XEPSO)
  NULLIFY(IP%XGAMM)
  NULLIFY(IP%XQDGAMM)
  NULLIFY(IP%XQDGMES)
  NULLIFY(IP%XT1GMES)
  NULLIFY(IP%XT2GMES)
  NULLIFY(IP%XAMAX)
  NULLIFY(IP%XQDAMAX)
  NULLIFY(IP%XT1AMAX)
  NULLIFY(IP%XT2AMAX)
  NULLIFY(IP%XAH)
  NULLIFY(IP%XBH)
  NULLIFY(IP%XTAU_WOOD)
  NULLIFY(IP%XINCREASE)
  NULLIFY(IP%XTURNOVER) 
  NULLIFY(IP%XCONDSAT)
  NULLIFY(IP%XMPOTSAT)
  NULLIFY(IP%XBCOEF)
  NULLIFY(IP%XWWILT)
  NULLIFY(IP%XWFC)  
  NULLIFY(IP%XWSAT)
  NULLIFY(IP%XTAUICE)
  NULLIFY(IP%XCGSAT)
  NULLIFY(IP%XC1SAT)
  NULLIFY(IP%XC2REF)
  NULLIFY(IP%XC3)
  NULLIFY(IP%XC4B)
  NULLIFY(IP%XACOEF)
  NULLIFY(IP%XPCOEF)  
  NULLIFY(IP%XC4REF)
  NULLIFY(IP%XWD0)
  NULLIFY(IP%XKANISO)  
  NULLIFY(IP%XPCPS)
  NULLIFY(IP%XPLVTT)
  NULLIFY(IP%XPLSTT)
  NULLIFY(IP%XHCAPSOIL)
  NULLIFY(IP%XCONDDRY)
  NULLIFY(IP%XCONDSLD)
  NULLIFY(IP%XRUNOFFD)
  NULLIFY(IP%XDZG)
  NULLIFY(IP%XDZDIF)
  NULLIFY(IP%XSOILWGHT)
  NULLIFY(IP%XFWTD)
  NULLIFY(IP%XWTD)    
  NULLIFY(IP%XKSAT_ICE)
  NULLIFY(IP%XALBNIR_DRY)
  NULLIFY(IP%XALBVIS_DRY)
  NULLIFY(IP%XALBUV_DRY)
  NULLIFY(IP%XALBNIR_WET)
  NULLIFY(IP%XALBVIS_WET)
  NULLIFY(IP%XALBUV_WET)
  NULLIFY(IP%XBSLAI_NITRO)
  IP%NSIZE_NATURE_P = 0
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_INIT1_N:ISBA_INIT_PGD_1P_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_INIT_PGD_1P_INIT
!
SUBROUTINE ISBA_INIT_PGD_NP_INIT(IP,KPATCH)
TYPE(ISBA_INIT_PGD_NP_t), INTENT(INOUT) :: IP
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_INIT1_N:ISBA_INIT_PGD_NP_INIT",0,ZHOOK_HANDLE)
 ALLOCATE(IP%AL(KPATCH))
DO JP=1,KPATCH
  CALL ISBA_INIT_PGD_1P_INIT(IP%AL(JP))
ENDDO    
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_INIT1_N:ISBA_INIT_PGD_NP_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_INIT_PGD_NP_INIT
!
SUBROUTINE ISBA_INIT_1P_INIT(I)
TYPE(ISBA_INIT_1P_t), INTENT(INOUT) :: I
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_INIT1_N:ISBA_INIT_1P_INIT",0,ZHOOK_HANDLE)
  NULLIFY(I%XMUF)
  NULLIFY(I%XFSAT)
  NULLIFY(I%XFFLOOD)
  NULLIFY(I%XPIFLOOD)  
  NULLIFY(I%XDIR_ALB_WITH_SNOW)
  NULLIFY(I%XSCA_ALB_WITH_SNOW)
  NULLIFY(I%XFF)
  NULLIFY(I%XFFG)
  NULLIFY(I%XFFV)
  NULLIFY(I%XFFROZEN)
  NULLIFY(I%XALBF)
  NULLIFY(I%XEMISF)
  NULLIFY(I%XTOPQS)
        
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_INIT1_N:ISBA_INIT_1P_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_INIT_1P_INIT
!
SUBROUTINE ISBA_INIT_NP_INIT(I,KPATCH)
TYPE(ISBA_INIT_NP_t), INTENT(INOUT) :: I
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_INIT1_N:ISBA_INIT_NP_INIT",0,ZHOOK_HANDLE)
 ALLOCATE(I%AL(KPATCH))
DO JP=1,KPATCH
  CALL ISBA_INIT_1P_INIT(I%AL(JP))
ENDDO         
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_INIT1_N:ISBA_INIT_NP_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_INIT_NP_INIT
!
END MODULE MODD_ISBA_INIT1_n
