!##################
MODULE MODD_ISBA_PARAM1_n
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
TYPE ISBA_PARAM_FIX_1P_t
!
REAL, POINTER, DIMENSION(:,:) :: XVEGTYPE
!
REAL, POINTER, DIMENSION(:,:) :: XDG           ! soil layer depth                  (m)
!                                                ! NOTE: in Force-Restore mode, the 
!                                                ! uppermost layer depth is superficial
!                                                ! and is only explicitly used for soil 
!                                                ! water phase changes                     (m)
!
REAL, POINTER, DIMENSION(:) :: XDG2
INTEGER, POINTER, DIMENSION(:) :: NWG_LAYER    ! Number of soil moisture layers for DIF
REAL, POINTER, DIMENSION(:)    :: XDROOT       ! effective root depth for DIF (m)
REAL, POINTER, DIMENSION(:,:)  :: XROOTFRAC    ! root fraction profile ('DIF' option)
!
REAL, POINTER, DIMENSION(:) :: XD_ICE          ! depth of the soil column for the calculation
!                                                 of the frozen soil fraction (m)
!
REAL, POINTER, DIMENSION(:) :: XH_TREE         ! height of trees                         (m)
!
REAL, POINTER, DIMENSION(:) :: XZ0_O_Z0H       ! ratio of surface roughness lengths
!                                                ! (momentum to heat)                      (-)
!
REAL, POINTER, DIMENSION(:) :: XRE25           ! Ecosystem respiration parameter         (kg/kg.m.s-1)
!
REAL, POINTER, DIMENSION(:) :: XDMAX           ! maximum air saturation deficit
!                                                ! tolerate by vegetation
!                                                (kg/kg)
!
END TYPE ISBA_PARAM_FIX_1P_t
!
!
TYPE ISBA_PARAM_FIX_NP_t
!
TYPE(ISBA_PARAM_FIX_1P_t), ALLOCATABLE :: AL(:)
!
END TYPE ISBA_PARAM_FIX_NP_t
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
TYPE ISBA_PARAM_TIME_1P_t
!
REAL, POINTER, DIMENSION(:) :: XVEG            ! vegetation cover fraction               (-)
!
REAL, POINTER, DIMENSION(:) :: XLAI          ! Leaf Area Index                         (m2/m2)
!
REAL, POINTER, DIMENSION(:) :: XEMIS         ! surface emissivity                      (-)
REAL, POINTER, DIMENSION(:) :: XZ0           ! surface roughness length                (m)
!
REAL, POINTER, DIMENSION(:) :: XRSMIN        ! minimum stomatal resistance             (s/m)
REAL, POINTER, DIMENSION(:) :: XGAMMA        ! coefficient for the calculation
!                                              ! of the surface stomatal
!                                              ! resistance
REAL, POINTER, DIMENSION(:) :: XWRMAX_CF     ! coefficient for maximum water 
!                                              ! interception 
!                                              ! storage capacity on the vegetation      (-)
REAL, POINTER, DIMENSION(:) :: XRGL          ! maximum solar radiation
!                                              ! usable in photosynthesis      
REAL, POINTER, DIMENSION(:) :: XCV           ! vegetation thermal inertia coefficient  (K m2/J)
REAL, POINTER, DIMENSION(:)    :: XLAIMIN    ! minimum LAI (Leaf Area Index)           (m2/m2)
REAL, POINTER, DIMENSION(:)    :: XSEFOLD    ! e-folding time for senescence           (s)
REAL, POINTER, DIMENSION(:)    :: XGMES      ! mesophyll conductance                   (m s-1)
REAL, POINTER, DIMENSION(:)    :: XGC        ! cuticular conductance                   (m s-1)
REAL, POINTER, DIMENSION(:)    :: XF2I       ! critical normilized soil water 
!                                              ! content for stress parameterisation
REAL, POINTER, DIMENSION(:)    :: XBSLAI     ! ratio d(biomass)/d(lai)                 (kg/m2)
!
REAL, POINTER, DIMENSION(:)    :: XCE_NITRO  ! leaf aera ratio sensitivity to 
!                                            ! nitrogen concentration                (m2/kg)
REAL, POINTER, DIMENSION(:)    :: XCF_NITRO  ! lethal minimum value of leaf area
!                                              ! ratio                                 (m2/kg)
REAL, POINTER, DIMENSION(:)    :: XCNA_NITRO ! nitrogen concentration of active 
!                                              ! biomass       
LOGICAL, POINTER, DIMENSION(:) :: LSTRESS    ! vegetation response type to water
!                                              ! stress (true:defensive false:offensive) (-)
!
REAL, POINTER, DIMENSION(:) :: XALBNIR_VEG       ! vegetation near-infra-red albedo        (-)
REAL, POINTER, DIMENSION(:) :: XALBVIS_VEG       ! vegetation visible albedo               (-)
REAL, POINTER, DIMENSION(:) :: XALBUV_VEG        ! vegetation UV albedo                    (-)
!
REAL, POINTER, DIMENSION(:) :: XALBNIR       ! near-infra-red albedo                   (-)
REAL, POINTER, DIMENSION(:) :: XALBVIS       ! visible albedo                          (-)
REAL, POINTER, DIMENSION(:) :: XALBUV        ! UV albedo
!
END TYPE ISBA_PARAM_TIME_1P_t
!
!
TYPE ISBA_PARAM_TIME_NP_t
!
TYPE(ISBA_PARAM_TIME_1P_t), ALLOCATABLE :: AL(:)
!
END TYPE ISBA_PARAM_TIME_NP_t
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
TYPE ISBA_PARAM_MEB_1P_t
!
! - Multi-energy balance (MEB) parameters.
! - Postfix GV denotes understory ground vegetation
!
REAL, POINTER, DIMENSION(:) :: XGNDLITTER        ! ground litter fraction                  (-)
REAL, POINTER, DIMENSION(:) :: XH_VEG            ! height of vegetation                           (m)
REAL, POINTER, DIMENSION(:) :: XZ0LITTER         ! ground litter roughness length                 (m)
!
END TYPE ISBA_PARAM_MEB_1P_t
!
TYPE ISBA_PARAM_MEB_NP_t
!
TYPE(ISBA_PARAM_MEB_1P_t), ALLOCATABLE :: AL(:)
!
END TYPE ISBA_PARAM_MEB_NP_t
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
TYPE ISBA_PARAM_ALB_1P_t
!
REAL, POINTER, DIMENSION(:) :: XALBNIR_SOIL      ! soil near-infra-red albedo              (-)
REAL, POINTER, DIMENSION(:) :: XALBVIS_SOIL      ! soil visible albedo                     (-)
REAL, POINTER, DIMENSION(:) :: XALBUV_SOIL       ! soil UV albedo
!
END TYPE ISBA_PARAM_ALB_1P_t
!
TYPE ISBA_PARAM_ALB_NP_t
!
TYPE(ISBA_PARAM_ALB_1P_t), ALLOCATABLE :: AL(:)
!
END TYPE ISBA_PARAM_ALB_NP_t
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
TYPE ISBA_PARAM_IRRIG_1P_t
!
! - Irrigation, seeding and reaping
!
TYPE (DATE_TIME), POINTER, DIMENSION(:)  :: TSEED          ! date of seeding
TYPE (DATE_TIME), POINTER, DIMENSION(:)  :: TREAP          ! date of reaping
REAL, POINTER, DIMENSION(:)         :: XWATSUP        ! water supply during irrigation process (mm)
REAL, POINTER, DIMENSION(:)         :: XIRRIG         ! flag for irrigation (irrigation if >0.)
!
END TYPE ISBA_PARAM_IRRIG_1P_t
!
TYPE ISBA_PARAM_IRRIG_NP_t
!
TYPE(ISBA_PARAM_IRRIG_1P_t), ALLOCATABLE :: AL(:)
!
END TYPE ISBA_PARAM_IRRIG_NP_t
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
TYPE ISBA_PARAM_NP_t
!
TYPE(ISBA_PARAM_FIX_NP_t) :: X
TYPE(ISBA_PARAM_TIME_NP_t) :: T
TYPE(ISBA_PARAM_MEB_NP_t) :: M
TYPE(ISBA_PARAM_ALB_NP_t) :: A
TYPE(ISBA_PARAM_IRRIG_NP_t) :: I
!
END TYPE ISBA_PARAM_NP_t
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
CONTAINS
!
SUBROUTINE ISBA_PARAM_TIME_1P_INIT(MT)
TYPE(ISBA_PARAM_TIME_1P_t), INTENT(INOUT) :: MT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_PARAM1_N:ISBA_PARAM_TIME_1P_INIT",0,ZHOOK_HANDLE)
!
  NULLIFY(MT%XLAI)  
  NULLIFY(MT%XVEG)
  NULLIFY(MT%XEMIS)
  NULLIFY(MT%XZ0)
  NULLIFY(MT%XRSMIN)
  NULLIFY(MT%XGAMMA)
  NULLIFY(MT%XWRMAX_CF)
  NULLIFY(MT%XRGL)  
  NULLIFY(MT%XCV)
  NULLIFY(MT%XLAIMIN)
  NULLIFY(MT%XSEFOLD)
  NULLIFY(MT%XGMES)
  NULLIFY(MT%XGC)
  NULLIFY(MT%XF2I)
  NULLIFY(MT%XBSLAI)
  NULLIFY(MT%XCE_NITRO)
  NULLIFY(MT%XCF_NITRO)
  NULLIFY(MT%XCNA_NITRO)
  NULLIFY(MT%LSTRESS)
  NULLIFY(MT%XALBNIR_VEG)
  NULLIFY(MT%XALBVIS_VEG)
  NULLIFY(MT%XALBUV_VEG)
  NULLIFY(MT%XALBNIR)
  NULLIFY(MT%XALBVIS)
  NULLIFY(MT%XALBUV)
!
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_PARAM1_N:ISBA_PARAM_TIME_1P_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_PARAM_TIME_1P_INIT
!
SUBROUTINE ISBA_PARAM_TIME_NP_INIT(MT,KPATCH)
TYPE(ISBA_PARAM_TIME_NP_t), INTENT(INOUT) :: MT
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_PARAM1_N:ISBA_PARAM_TIME_NP_INIT",0,ZHOOK_HANDLE)
!
ALLOCATE(MT%AL(KPATCH))
!
DO JP=1,KPATCH
  CALL ISBA_PARAM_TIME_1P_INIT(MT%AL(JP))
ENDDO
!
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_PARAM1_N:ISBA_PARAM_TIME_NP_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_PARAM_TIME_NP_INIT
!
SUBROUTINE ISBA_PARAM_FIX_1P_INIT(MX)
TYPE(ISBA_PARAM_FIX_1P_t), INTENT(INOUT) :: MX
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_PARAM1_N:ISBA_PARAM_FIX_1P_INIT",0,ZHOOK_HANDLE)
!
  NULLIFY(MX%XVEGTYPE)
  NULLIFY(MX%XDG)
  NULLIFY(MX%XDG2)
  NULLIFY(MX%NWG_LAYER)
  NULLIFY(MX%XDROOT)  
  NULLIFY(MX%XROOTFRAC)
  NULLIFY(MX%XD_ICE)
  NULLIFY(MX%XH_TREE)
  NULLIFY(MX%XZ0_O_Z0H)
  NULLIFY(MX%XRE25)  
  NULLIFY(MX%XDMAX)
!
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_PARAM1_N:ISBA_PARAM_FIX_1P_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_PARAM_FIX_1P_INIT
!
SUBROUTINE ISBA_PARAM_MEB_1P_INIT(MM)
TYPE(ISBA_PARAM_MEB_1P_t), INTENT(INOUT) :: MM
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_PARAM1_N:ISBA_PARAM_MEB_1P_INIT",0,ZHOOK_HANDLE)
!
  NULLIFY(MM%XGNDLITTER)
  NULLIFY(MM%XH_VEG)
  NULLIFY(MM%XZ0LITTER)
!
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_PARAM1_N:ISBA_PARAM_MEB_1P_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_PARAM_MEB_1P_INIT
!
SUBROUTINE ISBA_PARAM_ALB_1P_INIT(MA)
TYPE(ISBA_PARAM_ALB_1P_t), INTENT(INOUT) :: MA
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_PARAM1_N:ISBA_PARAM_ALB_1P_INIT",0,ZHOOK_HANDLE)
!
  NULLIFY(MA%XALBNIR_SOIL)
  NULLIFY(MA%XALBVIS_SOIL)
  NULLIFY(MA%XALBUV_SOIL)
!
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_PARAM1_N:ISBA_PARAM_ALB_1P_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_PARAM_ALB_1P_INIT
!
SUBROUTINE ISBA_PARAM_IRRIG_1P_INIT(MI)
TYPE(ISBA_PARAM_IRRIG_1P_t), INTENT(INOUT) :: MI
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_PARAM1_N:ISBA_PARAM_IRRIG_1P_INIT",0,ZHOOK_HANDLE)
!
  NULLIFY(MI%XWATSUP)
  NULLIFY(MI%XIRRIG)
!
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_PARAM1_N:ISBA_PARAM_IRRIG_1P_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_PARAM_IRRIG_1P_INIT
!
SUBROUTINE ISBA_PARAM_1P_INIT(MX, MM, MA, MI)
TYPE(ISBA_PARAM_FIX_1P_t), INTENT(INOUT) :: MX
TYPE(ISBA_PARAM_MEB_1P_t), INTENT(INOUT) :: MM
TYPE(ISBA_PARAM_ALB_1P_t), INTENT(INOUT) :: MA
TYPE(ISBA_PARAM_IRRIG_1P_t), INTENT(INOUT) :: MI
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_PARAM1_N:ISBA_PARAM_1P_INIT",0,ZHOOK_HANDLE)
!
 CALL ISBA_PARAM_FIX_1P_INIT(MX)
 CALL ISBA_PARAM_MEB_1P_INIT(MM)
 CALL ISBA_PARAM_ALB_1P_INIT(MA)
 CALL ISBA_PARAM_IRRIG_1P_INIT(MI)
!
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_PARAM1_N:ISBA_PARAM_1P_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_PARAM_1P_INIT
!
SUBROUTINE ISBA_PARAM_NP_INIT(MX, MM, MA, MI, KPATCH)
TYPE(ISBA_PARAM_FIX_NP_t), INTENT(INOUT) :: MX
TYPE(ISBA_PARAM_MEB_NP_t), INTENT(INOUT) :: MM
TYPE(ISBA_PARAM_ALB_NP_t), INTENT(INOUT) :: MA
TYPE(ISBA_PARAM_IRRIG_NP_t), INTENT(INOUT) :: MI
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_PARAM1_N:ISBA_PARAM_NP_INIT",0,ZHOOK_HANDLE)
!
ALLOCATE(MX%AL(KPATCH))
!
DO JP=1,KPATCH
  CALL ISBA_PARAM_FIX_1P_INIT(MX%AL(JP))
ENDDO
!
ALLOCATE(MM%AL(KPATCH))
!
DO JP=1,KPATCH
  CALL ISBA_PARAM_MEB_1P_INIT(MM%AL(JP))
ENDDO
!
ALLOCATE(MA%AL(KPATCH))
!
DO JP=1,KPATCH
  CALL ISBA_PARAM_ALB_1P_INIT(MA%AL(JP))
ENDDO
!
ALLOCATE(MI%AL(KPATCH))
!
DO JP=1,KPATCH
  CALL ISBA_PARAM_IRRIG_1P_INIT(MI%AL(JP))
ENDDO
!
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_PARAM1_N:ISBA_PARAM_NP_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_PARAM_NP_INIT

END MODULE MODD_ISBA_PARAM1_n
