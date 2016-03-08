!##################
MODULE MODD_ISBA_PARAM_n
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
TYPE ISBA_PARAM_FIX_t
!
REAL, POINTER, DIMENSION(:,:)  :: XVEGTYPE       ! fraction of each vegetation type for
!                                                ! each grid mesh
!
REAL, POINTER, DIMENSION(:,:,:) :: XDG           ! soil layer depth                  (m)
!                                                ! NOTE: in Force-Restore mode, the 
!                                                ! uppermost layer depth is superficial
!                                                ! and is only explicitly used for soil 
!                                                ! water phase changes                     (m)
REAL, POINTER, DIMENSION(:,:,:)  :: XDG_OLD      ! For land use
!
INTEGER, POINTER, DIMENSION(:,:) :: NWG_LAYER    ! Number of soil moisture layers for DIF
REAL, POINTER, DIMENSION(:,:)    :: XDROOT       ! effective root depth for DIF (m)
REAL, POINTER, DIMENSION(:,:)    :: XDG2         ! root depth for DIF as 3-L (m)
REAL, POINTER, DIMENSION(:,:,:)  :: XROOTFRAC    ! root fraction profile ('DIF' option)
!
REAL, POINTER, DIMENSION(:,:) :: XD_ICE          ! depth of the soil column for the calculation
!                                                 of the frozen soil fraction (m)
!
REAL, POINTER, DIMENSION(:,:) :: XH_TREE         ! height of trees                         (m)
!
REAL, POINTER, DIMENSION(:,:) :: XZ0_O_Z0H       ! ratio of surface roughness lengths
!                                                ! (momentum to heat)                      (-)
!
REAL, POINTER, DIMENSION(:,:) :: XRE25           ! Ecosystem respiration parameter         (kg/kg.m.s-1)
!
REAL, POINTER, DIMENSION(:,:) :: XDMAX           ! maximum air saturation deficit
!                                                ! tolerate by vegetation
!                                                (kg/kg)
!
END TYPE ISBA_PARAM_FIX_t
!
!
!TYPE ISBA_PARAM_FIX_t
!!
!REAL, POINTER, DIMENSION(:,:)  :: XVEGTYPE       ! fraction of each vegetation type for
!!                                                ! each grid mesh
!!
!TYPE(ISBA_PARAM_FIX_1P_t), ALLOCATABLE :: AL(:) ! => NULL()
!!TYPE(ISBA_PARAM_FIX_1P_t), POINTER :: CUR => NULL()
!!
!END TYPE ISBA_PARAM_FIX_t
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
TYPE ISBA_PARAM_TIME_t
!
REAL, POINTER, DIMENSION(:,:) :: XVEG            ! vegetation cover fraction               (-)
!
REAL, POINTER, DIMENSION(:,:) :: XLAI          ! Leaf Area Index                         (m2/m2)
!
REAL, POINTER, DIMENSION(:,:) :: XEMIS         ! surface emissivity                      (-)
REAL, POINTER, DIMENSION(:,:) :: XZ0           ! surface roughness length                (m)
!
REAL, POINTER, DIMENSION(:,:) :: XRSMIN        ! minimum stomatal resistance             (s/m)
REAL, POINTER, DIMENSION(:,:) :: XGAMMA        ! coefficient for the calculation
!                                              ! of the surface stomatal
!                                              ! resistance
REAL, POINTER, DIMENSION(:,:) :: XWRMAX_CF     ! coefficient for maximum water 
!                                              ! interception 
!                                              ! storage capacity on the vegetation      (-)
REAL, POINTER, DIMENSION(:,:) :: XRGL          ! maximum solar radiation
!                                              ! usable in photosynthesis      
REAL, POINTER, DIMENSION(:,:) :: XCV           ! vegetation thermal inertia coefficient  (K m2/J)
REAL, POINTER, DIMENSION(:,:)    :: XLAIMIN    ! minimum LAI (Leaf Area Index)           (m2/m2)
REAL, POINTER, DIMENSION(:,:)    :: XSEFOLD    ! e-folding time for senescence           (s)
REAL, POINTER, DIMENSION(:,:)    :: XGMES      ! mesophyll conductance                   (m s-1)
REAL, POINTER, DIMENSION(:,:)    :: XGC        ! cuticular conductance                   (m s-1)
REAL, POINTER, DIMENSION(:,:)    :: XF2I       ! critical normilized soil water 
!                                              ! content for stress parameterisation
REAL, POINTER, DIMENSION(:,:)    :: XBSLAI     ! ratio d(biomass)/d(lai)                 (kg/m2)
!
REAL, POINTER, DIMENSION(:,:)    :: XCE_NITRO  ! leaf aera ratio sensitivity to 
!                                              ! nitrogen concentration                (m2/kg)
REAL, POINTER, DIMENSION(:,:)    :: XCF_NITRO  ! lethal minimum value of leaf area
!                                              ! ratio                                 (m2/kg)
REAL, POINTER, DIMENSION(:,:)    :: XCNA_NITRO ! nitrogen concentration of active 
!                                              ! biomass       
LOGICAL, POINTER, DIMENSION(:,:) :: LSTRESS    ! vegetation response type to water
!                                              ! stress (true:defensive false:offensive) (-)
!
REAL, POINTER, DIMENSION(:,:) :: XALBNIR_VEG       ! vegetation near-infra-red albedo        (-)
REAL, POINTER, DIMENSION(:,:) :: XALBVIS_VEG       ! vegetation visible albedo               (-)
REAL, POINTER, DIMENSION(:,:) :: XALBUV_VEG        ! vegetation UV albedo                    (-)
!
REAL, POINTER, DIMENSION(:,:) :: XALBNIR       ! near-infra-red albedo                   (-)
REAL, POINTER, DIMENSION(:,:) :: XALBVIS       ! visible albedo                          (-)
REAL, POINTER, DIMENSION(:,:) :: XALBUV        ! UV albedo
!
END TYPE ISBA_PARAM_TIME_t
!
!TYPE ISBA_PARAM_TIME_t
!!
!TYPE(ISBA_PARAM_TIME_1P_t), ALLOCATABLE :: AL(:) !=> NULL()
!!TYPE(ISBA_PARAM_TIME_1P_t), POINTER :: CUR => NULL()
!!
!END TYPE ISBA_PARAM_TIME_t
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
TYPE ISBA_PARAM_MEB_t
!
! - Multi-energy balance (MEB) parameters.
! - Postfix GV denotes understory ground vegetation
!
REAL, POINTER, DIMENSION(:,:) :: XGNDLITTER        ! ground litter fraction                  (-)
REAL, POINTER, DIMENSION(:,:) :: XLAIGV            ! understory veg Leaf Area Index                 (m2/m2)
REAL, POINTER, DIMENSION(:,:) :: XH_VEG            ! height of vegetation                           (m)
REAL, POINTER, DIMENSION(:,:) :: XZ0LITTER         ! ground litter roughness length                 (m)
REAL, POINTER, DIMENSION(:,:) :: XRSMINGV          ! understory veg minimum stomatal resistance     (s/m)
REAL, POINTER, DIMENSION(:,:) :: XGAMMAGV          ! understory veg coefficient for the calculation
!                                                    ! of the surface stomatal resistance
REAL, POINTER, DIMENSION(:,:) :: XWRMAX_CFGV       ! understory veg coefficient for maximum water 
!                                                    ! interception
REAL, POINTER, DIMENSION(:,:) :: XRGLGV            ! understory veg maximum solar radiation
!                                                    ! usable in photosynthesis                       (W/m2)
REAL, POINTER, DIMENSION(:,:) :: XZ0GV
REAL, POINTER, DIMENSION(:,:,:) :: XROOTFRACGV     ! understory veg root fraction profile
!
END TYPE ISBA_PARAM_MEB_t
!
!TYPE ISBA_PARAM_MEB_t
!!
!TYPE(ISBA_PARAM_MEB_1P_t), ALLOCATABLE :: AL(:) !=> NULL()
!!TYPE(ISBA_PARAM_MEB_1P_t), POINTER :: CUR => NULL()
!!
!END TYPE ISBA_PARAM_MEB_t
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
TYPE ISBA_PARAM_ALB_t
!
REAL, POINTER, DIMENSION(:,:) :: XALBNIR_SOIL      ! soil near-infra-red albedo              (-)
REAL, POINTER, DIMENSION(:,:) :: XALBVIS_SOIL      ! soil visible albedo                     (-)
REAL, POINTER, DIMENSION(:,:) :: XALBUV_SOIL       ! soil UV albedo
!
END TYPE ISBA_PARAM_ALB_t
!
!TYPE ISBA_PARAM_ALB_t
!!
!TYPE(ISBA_PARAM_ALB_1P_t), ALLOCATABLE :: AL(:) !=> NULL()
!!TYPE(ISBA_PARAM_ALB_1P_t), POINTER :: CUR => NULL()
!!
!END TYPE ISBA_PARAM_ALB_t
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
TYPE ISBA_PARAM_IRRIG_t
!
! - Irrigation, seeding and reaping
!
TYPE (DATE_TIME), POINTER, DIMENSION(:,:)  :: TSEED          ! date of seeding
TYPE (DATE_TIME), POINTER, DIMENSION(:,:)  :: TREAP          ! date of reaping
REAL, POINTER, DIMENSION(:,:)         :: XWATSUP        ! water supply during irrigation process (mm)
REAL, POINTER, DIMENSION(:,:)         :: XIRRIG         ! flag for irrigation (irrigation if >0.)
!
END TYPE ISBA_PARAM_IRRIG_t
!
!TYPE ISBA_PARAM_IRRIG_t
!!
!TYPE(ISBA_PARAM_IRRIG_1P_t), ALLOCATABLE :: AL(:) !=> NULL()
!!TYPE(ISBA_PARAM_IRRIG_1P_t), POINTER :: CUR => NULL()
!!
!END TYPE ISBA_PARAM_IRRIG_t
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
TYPE ISBA_PARAM_t
!
TYPE(ISBA_PARAM_FIX_t) :: X
TYPE(ISBA_PARAM_TIME_t) :: T
TYPE(ISBA_PARAM_MEB_t) :: M
TYPE(ISBA_PARAM_ALB_t) :: A
TYPE(ISBA_PARAM_IRRIG_t) :: I
!
END TYPE ISBA_PARAM_t
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
CONTAINS
!
SUBROUTINE ISBA_PARAM_TIME_INIT(YISBA_PARAM_TIME)
TYPE(ISBA_PARAM_TIME_t), INTENT(INOUT) :: YISBA_PARAM_TIME
!INTEGER, INTENT(IN) :: KPATCH
!INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_PARAM_N:ISBA_PARAM_TIME_INIT",0,ZHOOK_HANDLE)
!
!ALLOCATE(YISBA_PARAM_TIME%AL(KPATCH))
!YISBA_PARAM_TIME%CUR => YISBA_PARAM_TIME%AL(1)
!DO JP=1,KPATCH
  NULLIFY(YISBA_PARAM_TIME%XLAI)  
  NULLIFY(YISBA_PARAM_TIME%XVEG)
  NULLIFY(YISBA_PARAM_TIME%XEMIS)
  NULLIFY(YISBA_PARAM_TIME%XZ0)
  NULLIFY(YISBA_PARAM_TIME%XRSMIN)
  NULLIFY(YISBA_PARAM_TIME%XGAMMA)
  NULLIFY(YISBA_PARAM_TIME%XWRMAX_CF)
  NULLIFY(YISBA_PARAM_TIME%XRGL)  
  NULLIFY(YISBA_PARAM_TIME%XCV)
  NULLIFY(YISBA_PARAM_TIME%XLAIMIN)
  NULLIFY(YISBA_PARAM_TIME%XSEFOLD)
  NULLIFY(YISBA_PARAM_TIME%XGMES)
  NULLIFY(YISBA_PARAM_TIME%XGC)
  NULLIFY(YISBA_PARAM_TIME%XF2I)
  NULLIFY(YISBA_PARAM_TIME%XBSLAI)
  NULLIFY(YISBA_PARAM_TIME%XCE_NITRO)
  NULLIFY(YISBA_PARAM_TIME%XCF_NITRO)
  NULLIFY(YISBA_PARAM_TIME%XCNA_NITRO)
  NULLIFY(YISBA_PARAM_TIME%LSTRESS)
  NULLIFY(YISBA_PARAM_TIME%XALBNIR_VEG)
  NULLIFY(YISBA_PARAM_TIME%XALBVIS_VEG)
  NULLIFY(YISBA_PARAM_TIME%XALBUV_VEG)
  NULLIFY(YISBA_PARAM_TIME%XALBNIR)
  NULLIFY(YISBA_PARAM_TIME%XALBVIS)
  NULLIFY(YISBA_PARAM_TIME%XALBUV)
!ENDDO
!
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_PARAM_N:ISBA_PARAM_TIME_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_PARAM_TIME_INIT
!
SUBROUTINE ISBA_PARAM_INIT(YISBA_PARAM_FIX, YISBA_PARAM_MEB, YISBA_PARAM_ALB,  YISBA_PARAM_IRR)
TYPE(ISBA_PARAM_FIX_t), INTENT(INOUT) :: YISBA_PARAM_FIX
TYPE(ISBA_PARAM_MEB_t), INTENT(INOUT) :: YISBA_PARAM_MEB
TYPE(ISBA_PARAM_ALB_t), INTENT(INOUT) :: YISBA_PARAM_ALB
TYPE(ISBA_PARAM_IRRIG_t), INTENT(INOUT) :: YISBA_PARAM_IRR
!INTEGER, INTENT(IN) :: KPATCH
!INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_PARAM_N:ISBA_PARAM_INIT",0,ZHOOK_HANDLE)
!
NULLIFY(YISBA_PARAM_FIX%XVEGTYPE)
!
!ALLOCATE(YISBA_PARAM_FIX%AL(KPATCH))
!YISBA_PARAM_FIX%CUR => YISBA_PARAM_FIX%AL(1)
!DO JP=1,KPATCH
  NULLIFY(YISBA_PARAM_FIX%XDG)
  NULLIFY(YISBA_PARAM_FIX%XDG_OLD)
  NULLIFY(YISBA_PARAM_FIX%NWG_LAYER)
  NULLIFY(YISBA_PARAM_FIX%XDROOT)
  NULLIFY(YISBA_PARAM_FIX%XDG2)
  NULLIFY(YISBA_PARAM_FIX%XROOTFRAC)
  NULLIFY(YISBA_PARAM_FIX%XD_ICE)
  NULLIFY(YISBA_PARAM_FIX%XH_TREE)
  NULLIFY(YISBA_PARAM_FIX%XZ0_O_Z0H)
  NULLIFY(YISBA_PARAM_FIX%XRE25)  
  NULLIFY(YISBA_PARAM_FIX%XDMAX)
!ENDDO
!
!ALLOCATE(YISBA_PARAM_MEB%AL(KPATCH))
!YISBA_PARAM_MEB%CUR => YISBA_PARAM_MEB%AL(1)
!DO JP=1,KPATCH
  NULLIFY(YISBA_PARAM_MEB%XGNDLITTER)
  NULLIFY(YISBA_PARAM_MEB%XLAIGV)
  NULLIFY(YISBA_PARAM_MEB%XH_VEG)
  NULLIFY(YISBA_PARAM_MEB%XZ0LITTER)
  NULLIFY(YISBA_PARAM_MEB%XRSMINGV)
  NULLIFY(YISBA_PARAM_MEB%XGAMMAGV)
  NULLIFY(YISBA_PARAM_MEB%XWRMAX_CFGV)  
  NULLIFY(YISBA_PARAM_MEB%XRGLGV)
  NULLIFY(YISBA_PARAM_MEB%XZ0GV)
  NULLIFY(YISBA_PARAM_MEB%XROOTFRACGV)
!ENDDO
!
!ALLOCATE(YISBA_PARAM_ALB%AL(KPATCH))
!YISBA_PARAM_ALB%CUR => YISBA_PARAM_ALB_%AL(1)
!DO JP=1,KPATCH
  NULLIFY(YISBA_PARAM_ALB%XALBNIR_SOIL)
  NULLIFY(YISBA_PARAM_ALB%XALBVIS_SOIL)
  NULLIFY(YISBA_PARAM_ALB%XALBUV_SOIL)
!ENDDO
!
!ALLOCATE(YISBA_PARAM_IRR%AL(KPATCH))
!YISBA_PARAM_IRR%CUR => YISBA_PARAM_ALB%AL(1)
!DO JP=1,KPATCH
  NULLIFY(YISBA_PARAM_IRR%XWATSUP)
  NULLIFY(YISBA_PARAM_IRR%XIRRIG)
!ENDDO
!
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_PARAM_N:ISBA_PARAM_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_PARAM_INIT

END MODULE MODD_ISBA_PARAM_n
