!     ################
      MODULE MODD_TEB_n
!     ################
!
!!****  *MODD_TEB_n - declaration of surface parameters for urban surface
!!
!!    PURPOSE
!!    -------
!     Declaration of surface parameters
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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2004
!!      A. Lemonsu      07/2012         Key for urban hydrology
!!      V. Masson       06/2013         splits module in two
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TYPE_SNOW
USE MODD_TYPE_DATE_SURF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!--------------------------------------------------------------------------

TYPE TEB_t
! TEB scheme option
!
! Geometric Parameters:
!
  REAL, POINTER, DIMENSION(:)   :: XROAD_DIR     ! Road direction (deg from North, clockwise)
  REAL, POINTER, DIMENSION(:)   :: XGARDEN       ! fraction of veg in the streets   (-)
  REAL, POINTER, DIMENSION(:)   :: XGREENROOF    ! fraction of greenroofs on roofs  (-)
  REAL, POINTER, DIMENSION(:)   :: XBLD          ! fraction of buildings            (-)
  REAL, POINTER, DIMENSION(:)   :: XROAD         ! fraction of roads                (-)
  REAL, POINTER, DIMENSION(:)   :: XCAN_HW_RATIO ! canyon    h/W                    (-)
  REAL, POINTER, DIMENSION(:)   :: XBLD_HEIGHT   ! buildings height 'h'             (m)
  REAL, POINTER, DIMENSION(:)   :: XWALL_O_HOR   ! wall surf. / hor. surf.          (-)
  REAL, POINTER, DIMENSION(:)   :: XROAD_O_GRND  ! road surf. / (road + garden surf.) (-)
  REAL, POINTER, DIMENSION(:)   :: XGARDEN_O_GRND! gard. surf. / (road + garden surf.)(-)
  REAL, POINTER, DIMENSION(:)   :: XWALL_O_GRND  ! wall surf. / (road + garden surf.) (-)
  REAL, POINTER, DIMENSION(:)   :: XWALL_O_BLD   ! wall surf. / bld surf. (-)
  REAL, POINTER, DIMENSION(:)   :: XZ0_TOWN      ! roughness length for momentum    (m)
  REAL, POINTER, DIMENSION(:)   :: XSVF_ROAD     ! road sky view factor             (-)
  REAL, POINTER, DIMENSION(:)   :: XSVF_GARDEN   ! green area sky view factor       (-)
  REAL, POINTER, DIMENSION(:)   :: XSVF_WALL     ! wall sky view factor             (-)
!
! Roof parameters
!
  REAL, POINTER, DIMENSION(:)   :: XALB_ROOF     ! roof albedo                      (-)
  REAL, POINTER, DIMENSION(:)   :: XEMIS_ROOF    ! roof emissivity                  (-)
  REAL, POINTER, DIMENSION(:,:) :: XHC_ROOF      ! roof layers heat capacity        (J/K/m3)
  REAL, POINTER, DIMENSION(:,:) :: XTC_ROOF      ! roof layers thermal conductivity (W/K/m)
  REAL, POINTER, DIMENSION(:,:) :: XD_ROOF       ! depth of roof layers             (m)
  REAL, POINTER, DIMENSION(:)   :: XROUGH_ROOF   ! roof roughness coef
!
!
! Road parameters
!
  REAL, POINTER, DIMENSION(:)   :: XALB_ROAD     ! road albedo                      (-)
  REAL, POINTER, DIMENSION(:)   :: XEMIS_ROAD    ! road emissivity                  (-)
  REAL, POINTER, DIMENSION(:,:) :: XHC_ROAD      ! road layers heat capacity        (J/K/m3)
  REAL, POINTER, DIMENSION(:,:) :: XTC_ROAD      ! road layers thermal conductivity (W/K/m)
  REAL, POINTER, DIMENSION(:,:) :: XD_ROAD       ! depth of road layers             (m)
!
! Wall parameters
!
  REAL, POINTER, DIMENSION(:)   :: XALB_WALL     ! wall albedo                      (-)
  REAL, POINTER, DIMENSION(:)   :: XEMIS_WALL    ! wall emissivity                  (-)
  REAL, POINTER, DIMENSION(:,:) :: XHC_WALL      ! wall layers heat capacity        (J/K/m3)
  REAL, POINTER, DIMENSION(:,:) :: XTC_WALL      ! wall layers thermal conductivity (W/K/m)
  REAL, POINTER, DIMENSION(:,:) :: XD_WALL       ! depth of wall layers             (m)
  REAL, POINTER, DIMENSION(:)   :: XROUGH_WALL   ! wall roughness coef
!
! Building's use type
!
  REAL, POINTER, DIMENSION(:)   :: XRESIDENTIAL  ! fraction of Residential use      (-)
  REAL                          :: XDT_RES       ! target temperature change when unoccupied (K) (residential buildings)
  REAL                          :: XDT_OFF       ! target temperature change when unoccupied (K) (offices buildings)
  
!
! anthropogenic fluxes
!
  REAL, POINTER, DIMENSION(:)   :: XH_TRAFFIC    ! anthropogenic sensible
!                                                  ! heat fluxes due to traffic       (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XLE_TRAFFIC   ! anthropogenic latent
!                                                  ! heat fluxes due to traffic       (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XH_INDUSTRY   ! anthropogenic sensible                   
!                                                  ! heat fluxes due to factories     (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XLE_INDUSTRY  ! anthropogenic latent
!                                                  ! heat fluxes due to factories     (W/m2)
!
! temperatures for boundary conditions
!
  REAL, POINTER, DIMENSION(:)   :: XTI_ROAD      ! road interior temperature        (K)
!
! Prognostic variables:
!
  REAL, POINTER, DIMENSION(:)   :: XWS_ROOF      ! roof water reservoir             (kg/m2)
  REAL, POINTER, DIMENSION(:)   :: XWS_ROAD      ! road water reservoir             (kg/m2)
  REAL, POINTER, DIMENSION(:,:) :: XT_ROOF       ! roof layer temperatures          (K)
  REAL, POINTER, DIMENSION(:,:) :: XT_ROAD       ! road layer temperatures          (K)
  REAL, POINTER, DIMENSION(:,:) :: XT_WALL_A     ! wall layer temperatures          (K)
  REAL, POINTER, DIMENSION(:,:) :: XT_WALL_B     ! wall layer temperatures          (K)
!
  REAL, POINTER, DIMENSION(:)   :: XAC_ROOF      ! roof aerodynamic conductance     ()
  REAL, POINTER, DIMENSION(:)   :: XAC_ROAD      ! road aerodynamic conductance     ()
  REAL, POINTER, DIMENSION(:)   :: XAC_WALL      ! wall aerodynamic conductance     ()
  REAL, POINTER, DIMENSION(:)   :: XAC_TOP       ! top  aerodynamic conductance     ()
  REAL, POINTER, DIMENSION(:)   :: XAC_ROOF_WAT  ! water aerodynamic conductance    ()
  REAL, POINTER, DIMENSION(:)   :: XAC_ROAD_WAT  ! water aerodynamic conductance    ()
!
  REAL, POINTER, DIMENSION(:)   :: XQSAT_ROOF    ! humidity of saturation for roofs (kg/kg)
  REAL, POINTER, DIMENSION(:)   :: XQSAT_ROAD    ! humidity of saturation for roads (kg/kg)
!
  REAL, POINTER, DIMENSION(:)   :: XDELT_ROOF    ! humidity of saturation for roofs (-)
  REAL, POINTER, DIMENSION(:)   :: XDELT_ROAD    ! humidity of saturation for roads (-)
!
! Semi-prognostic variables:
!
  REAL, POINTER, DIMENSION(:)   :: XT_CANYON     ! canyon air temperature           (K)
  REAL, POINTER, DIMENSION(:)   :: XQ_CANYON     ! canyon air specific humidity     (kg/kg)
!
!
! Prognostic snow:
!
  TYPE(SURF_SNOW)                 :: TSNOW_ROOF      ! snow state on roofs: 
!                                                  ! scheme type/option               (-)
!                                                  ! number of layers                 (-)
!                                                  ! snow (& liq. water) content      (kg/m2)
!                                                  ! heat content                     (J/m2)
!                                                  ! temperature                      (K)
!                                                  ! density                          (kg m-3)
  TYPE(SURF_SNOW)                 :: TSNOW_ROAD      ! snow state on roads: 
!                                                  ! scheme type/option               (-)
!                                                  ! number of layers                 (-)
!                                                  ! snow (& liq. water) content      (kg/m2)
!                                                  ! heat content                     (J/m2)
!                                                  ! temperature                      (K)
!                                                  ! density                          (kg m-3)
!                                                  ! density                          (kg m-3)
  TYPE(SURF_SNOW)                 :: TSNOW_GARDEN    ! snow state on green areas:
!                                                  ! scheme type/option               (-)
!                                                  ! number of layers                 (-)
!                                                  ! snow (& liq. water) content      (kg/m2)
!                                                  ! heat content                     (J/m2)
!                                                  ! temperature                      (K)
!                                                  ! density                          (kg m-3)
!
END TYPE TEB_t

TYPE(TEB_t), ALLOCATABLE, TARGET, SAVE :: TEB_MODEL(:,:)

TYPE(TEB_t), POINTER :: TEB => NULL()
!$OMP THREADPRIVATE(TEB)

CONTAINS
!----------------------------------------------------------------------------

SUBROUTINE TEB_GOTO_MODEL(KFROM, KTO, LKFROM, KFROM_PATCH, KTO_PATCH)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
INTEGER, INTENT(IN) :: KFROM_PATCH, KTO_PATCH
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!


! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_TEB_N:TEB_GOTO_MODEL',0,ZHOOK_HANDLE)

TEB => TEB_MODEL(KTO,KTO_PATCH)

IF (LHOOK) CALL DR_HOOK('MODD_TEB_N:TEB_GOTO_MODEL',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE TEB_GOTO_MODEL

SUBROUTINE TEB_ALLOC(KMODEL,KPATCH)
INTEGER, INTENT(IN) :: KMODEL
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: J,JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_N:TEB_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(TEB_MODEL(KMODEL,KPATCH))
DO J=1,KMODEL
 DO JP=1,KPATCH
  NULLIFY(TEB_MODEL(J,JP)%XROAD_DIR)
  NULLIFY(TEB_MODEL(J,JP)%XGARDEN)
  NULLIFY(TEB_MODEL(J,JP)%XGREENROOF)
  NULLIFY(TEB_MODEL(J,JP)%XBLD)
  NULLIFY(TEB_MODEL(J,JP)%XROAD)
  NULLIFY(TEB_MODEL(J,JP)%XCAN_HW_RATIO)
  NULLIFY(TEB_MODEL(J,JP)%XBLD_HEIGHT)
  NULLIFY(TEB_MODEL(J,JP)%XWALL_O_HOR)
  NULLIFY(TEB_MODEL(J,JP)%XROAD_O_GRND)
  NULLIFY(TEB_MODEL(J,JP)%XGARDEN_O_GRND)
  NULLIFY(TEB_MODEL(J,JP)%XWALL_O_GRND)
  NULLIFY(TEB_MODEL(J,JP)%XWALL_O_BLD)
  NULLIFY(TEB_MODEL(J,JP)%XZ0_TOWN)
  NULLIFY(TEB_MODEL(J,JP)%XSVF_ROAD)
  NULLIFY(TEB_MODEL(J,JP)%XSVF_GARDEN)
  NULLIFY(TEB_MODEL(J,JP)%XSVF_WALL)
  NULLIFY(TEB_MODEL(J,JP)%XALB_ROOF)
  NULLIFY(TEB_MODEL(J,JP)%XEMIS_ROOF)
  NULLIFY(TEB_MODEL(J,JP)%XHC_ROOF)
  NULLIFY(TEB_MODEL(J,JP)%XTC_ROOF)
  NULLIFY(TEB_MODEL(J,JP)%XD_ROOF)
  NULLIFY(TEB_MODEL(J,JP)%XALB_ROAD)
  NULLIFY(TEB_MODEL(J,JP)%XEMIS_ROAD)
  NULLIFY(TEB_MODEL(J,JP)%XHC_ROAD)
  NULLIFY(TEB_MODEL(J,JP)%XTC_ROAD)
  NULLIFY(TEB_MODEL(J,JP)%XD_ROAD)
  NULLIFY(TEB_MODEL(J,JP)%XALB_WALL)
  NULLIFY(TEB_MODEL(J,JP)%XEMIS_WALL)
  NULLIFY(TEB_MODEL(J,JP)%XHC_WALL)
  NULLIFY(TEB_MODEL(J,JP)%XTC_WALL)
  NULLIFY(TEB_MODEL(J,JP)%XD_WALL)
  NULLIFY(TEB_MODEL(J,JP)%XH_TRAFFIC)
  NULLIFY(TEB_MODEL(J,JP)%XLE_TRAFFIC)
  NULLIFY(TEB_MODEL(J,JP)%XH_INDUSTRY)
  NULLIFY(TEB_MODEL(J,JP)%XLE_INDUSTRY)
  NULLIFY(TEB_MODEL(J,JP)%XTI_ROAD)
  NULLIFY(TEB_MODEL(J,JP)%XWS_ROOF)
  NULLIFY(TEB_MODEL(J,JP)%XWS_ROAD)
  NULLIFY(TEB_MODEL(J,JP)%XT_ROOF)
  NULLIFY(TEB_MODEL(J,JP)%XT_ROAD)
  NULLIFY(TEB_MODEL(J,JP)%XT_WALL_A)
  NULLIFY(TEB_MODEL(J,JP)%XT_WALL_B)
  NULLIFY(TEB_MODEL(J,JP)%XAC_ROOF)
  NULLIFY(TEB_MODEL(J,JP)%XAC_ROAD)
  NULLIFY(TEB_MODEL(J,JP)%XAC_WALL)
  NULLIFY(TEB_MODEL(J,JP)%XAC_TOP)
  NULLIFY(TEB_MODEL(J,JP)%XAC_ROOF_WAT)
  NULLIFY(TEB_MODEL(J,JP)%XAC_ROAD_WAT)
  NULLIFY(TEB_MODEL(J,JP)%XQSAT_ROOF)
  NULLIFY(TEB_MODEL(J,JP)%XQSAT_ROAD)
  NULLIFY(TEB_MODEL(J,JP)%XDELT_ROOF)
  NULLIFY(TEB_MODEL(J,JP)%XDELT_ROAD)
  NULLIFY(TEB_MODEL(J,JP)%XT_CANYON)
  NULLIFY(TEB_MODEL(J,JP)%XQ_CANYON)
  NULLIFY(TEB_MODEL(J,JP)%XROUGH_ROOF)
  NULLIFY(TEB_MODEL(J,JP)%XROUGH_WALL)
  NULLIFY(TEB_MODEL(J,JP)%XRESIDENTIAL)
 ENDDO
ENDDO
TEB_MODEL(:,:)%XDT_RES=0.
TEB_MODEL(:,:)%XDT_OFF=0.

IF (LHOOK) CALL DR_HOOK("MODD_TEB_N:TEB_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_ALLOC

SUBROUTINE TEB_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_N:TEB_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(TEB_MODEL)) DEALLOCATE(TEB_MODEL)
IF (ASSOCIATED(TEB)) NULLIFY(TEB)
IF (LHOOK) CALL DR_HOOK("MODD_TEB_N:TEB_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_DEALLO

END MODULE MODD_TEB_n
