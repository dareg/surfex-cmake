!     ################
      MODULE MODD_BLD_DESCRIPTION_n
!     ################
!
!!****  *MODD_BLD_DESCRIPTION_n - declaration of surface parameters for typical
!                               buildings
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
!!      G. Pigeon   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       08/2011
!!       V. Masson     08/2013 adds solar panels
!!       V. Masson     10/2013 adds residential fraction
!!----------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE BLD_DESC_t
!
! Number of layers
!
  INTEGER                       :: NDESC_BLD          ! number of types of buildings
  INTEGER                       :: NDESC_AGE          ! number of building's construction dates ranges
  INTEGER                       :: NDESC_CODE         ! number of codes for buildings (merges type & age)
  INTEGER                       :: NDESC_USE          ! number of types of building's uses
  INTEGER                       :: NDESC_ROOF_LAYER   ! number of layers in roofs
  INTEGER                       :: NDESC_ROAD_LAYER   ! number of layers in roads
  INTEGER                       :: NDESC_WALL_LAYER   ! number of layers in walls
  INTEGER                       :: NDESC_FLOOR_LAYER  ! number of layers in floor
  INTEGER, POINTER, DIMENSION(:):: NDESC_BLD_LIST     ! list of the types for buildings
  INTEGER, POINTER, DIMENSION(:):: NDESC_AGE_DATE     ! list of the contruction dates for buildings
  INTEGER, POINTER, DIMENSION(:):: NDESC_CODE_LIST    ! list of the codes for buildings
  INTEGER, POINTER, DIMENSION(:):: NDESC_AGE_LIST     ! list of the contruction dates' codes
  INTEGER, POINTER, DIMENSION(:):: NDESC_USE_LIST     ! list of the codes for building's uses
  !
  REAL, POINTER, DIMENSION(:)   :: XDESC_ALB_ROOF     ! Roof albedo
  REAL, POINTER, DIMENSION(:)   :: XDESC_ALB_ROAD     ! Road albedo
  REAL, POINTER, DIMENSION(:)   :: XDESC_ALB_WALL     ! Wall albedo
  REAL, POINTER, DIMENSION(:)   :: XDESC_EMIS_ROOF    ! Roof emissivity
  REAL, POINTER, DIMENSION(:)   :: XDESC_EMIS_ROAD    ! Road emissivity
  REAL, POINTER, DIMENSION(:)   :: XDESC_EMIS_WALL    ! Wall emissivity
  REAL, POINTER, DIMENSION(:)   :: XDESC_TCOOL_TARGET ! cooling setpoint of indoor air
  REAL, POINTER, DIMENSION(:)   :: XDESC_THEAT_TARGET ! heating setpoint of indoor air
  REAL, POINTER, DIMENSION(:)   :: XDESC_F_WASTE_CAN  ! fraction of waste heat into the canyon
  REAL, POINTER, DIMENSION(:)   :: XDESC_EFF_HEAT     ! efficiency of the heating system
  REAL, POINTER, DIMENSION(:,:) :: XDESC_HC_FLOOR     ! heat capacity of floor layers [J m-3 K-1]
  REAL, POINTER, DIMENSION(:,:) :: XDESC_TC_FLOOR     ! thermal conductivity of floor layers [W m-1 K-1]
  REAL, POINTER, DIMENSION(:,:) :: XDESC_D_FLOOR      ! thickness of floor layers [m]
  REAL, POINTER, DIMENSION(:,:) :: XDESC_HC_ROOF      ! heat capacity of roof layers [J m-3 K-1]
  REAL, POINTER, DIMENSION(:,:) :: XDESC_TC_ROOF      ! thermal conductivity of roof layers [W m-1 K-1]
  REAL, POINTER, DIMENSION(:,:) :: XDESC_D_ROOF       ! thickness of roof layers [m]
  REAL, POINTER, DIMENSION(:,:) :: XDESC_HC_ROAD      ! heat capacity of road layers [J m-3 K-1]
  REAL, POINTER, DIMENSION(:,:) :: XDESC_TC_ROAD      ! thermal conductivity of road layers [W m-1 K-1]
  REAL, POINTER, DIMENSION(:,:) :: XDESC_D_ROAD       ! thickness of road layers [m]
  REAL, POINTER, DIMENSION(:,:) :: XDESC_HC_WALL      ! heat capacity of wall layers [J m-3 K-1]
  REAL, POINTER, DIMENSION(:,:) :: XDESC_TC_WALL      ! thermal conductivity of wall layers [W m-1 K-1]
  REAL, POINTER, DIMENSION(:,:) :: XDESC_D_WALL       ! thickness of wall layers [m]
  REAL, POINTER, DIMENSION(:)   :: XDESC_QIN          ! internal heat gains [W m-2(floor)]
  REAL, POINTER, DIMENSION(:)   :: XDESC_QIN_FRAD     ! radiant fraction of internal heat gains
  REAL, POINTER, DIMENSION(:)   :: XDESC_SHGC         ! solar transmitance of windows
  REAL, POINTER, DIMENSION(:)   :: XDESC_U_WIN        ! glazing thermal resistance [K m W-2]
  REAL, POINTER, DIMENSION(:)   :: XDESC_GR           ! glazing ratio
  REAL, POINTER, DIMENSION(:)   :: XDESC_SHGC_SH      ! solar transmitance of windows + shading
  REAL, POINTER, DIMENSION(:)   :: XDESC_FLOOR_HEIGHT ! building floor height [m]
  REAL, POINTER, DIMENSION(:)   :: XDESC_INF          ! infiltration/ventilation flow rate [AC/H]
!
  REAL, POINTER, DIMENSION(:)   :: XDESC_F_WATER_COND ! fraction of evaporation for condensers
  REAL, POINTER, DIMENSION(:)   :: XDESC_SHADE        ! Flag to activate shading devices 0->inactivated , 1->activated
  REAL, POINTER, DIMENSION(:)   :: XDESC_NATVENT      ! Flag to describe bld surventilation solution 0-> NONE ; 1 -> MANU ; 2-> AUTO
  REAL, POINTER, DIMENSION(:)   :: XDESC_QIN_FLAT     ! Latent franction of internal heat gains
  REAL, POINTER, DIMENSION(:)   :: XDESC_HR_TARGET    ! Relative humidity setpoint
  REAL, POINTER, DIMENSION(:)   :: XDESC_V_VENT       ! Ventilation flow rate [AC/H]
  REAL, POINTER, DIMENSION(:)   :: XDESC_COP_RAT      ! Rated COP of the cooling system
  REAL, POINTER, DIMENSION(:)   :: XDESC_GREENROOF    ! Greenroof fraction
  REAL, POINTER, DIMENSION(:)   :: XDESC_EMIS_PANEL   ! Emissivity of Solar panels
  REAL, POINTER, DIMENSION(:)   :: XDESC_ALB_PANEL    ! Albedo     of Solar panels
  REAL, POINTER, DIMENSION(:)   :: XDESC_EFF_PANEL    ! Efficiency of Solar panels
  REAL, POINTER, DIMENSION(:)   :: XDESC_FRAC_PANEL   ! Fraction   of Solar panels on roofs
  REAL, POINTER, DIMENSION(:)   :: XDESC_RESIDENTIAL  ! Fraction of residential use
!
END TYPE BLD_DESC_t

TYPE(BLD_DESC_t), ALLOCATABLE, TARGET, SAVE :: BLD_DESC_MODEL(:)

TYPE(BLD_DESC_t), POINTER :: BLD_DESC => NULL()
!$OMP THREADPRIVATE(BLD_DESC)

CONTAINS

SUBROUTINE BLD_DESC_GOTO_MODEL(KFROM, KTO, LKFROM)
INTEGER, INTENT(IN) :: KFROM, KTO
LOGICAL, INTENT(IN) :: LKFROM
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
!
IF (LHOOK) CALL DR_HOOK('MODD_BLD_DESCRIPTION_n:BLD_DESC_GOTO_MODEL',0,ZHOOK_HANDLE)

BLD_DESC => BLD_DESC_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_BLD_DESCRIPTION_n:BLD_DESC_GOTO_MODEL',1,ZHOOK_HANDLE)
!
END SUBROUTINE BLD_DESC_GOTO_MODEL
!
SUBROUTINE BLD_DESC_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_BLD_DESCRIPTION_n:BLD_DESC_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(BLD_DESC_MODEL(KMODEL))
BLD_DESC => BLD_DESC_MODEL(KMODEL)
DO J=1,KMODEL
  NULLIFY(BLD_DESC_MODEL(J)%NDESC_BLD_LIST)
  NULLIFY(BLD_DESC_MODEL(J)%NDESC_CODE_LIST)
  NULLIFY(BLD_DESC_MODEL(J)%NDESC_AGE_LIST)
  NULLIFY(BLD_DESC_MODEL(J)%NDESC_AGE_DATE)
  NULLIFY(BLD_DESC_MODEL(J)%NDESC_USE_LIST)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_ALB_ROOF)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_ALB_ROAD)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_ALB_WALL)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_EMIS_ROOF)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_EMIS_ROAD)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_EMIS_WALL)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_TCOOL_TARGET)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_THEAT_TARGET)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_F_WASTE_CAN)  
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_EFF_HEAT)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_HC_FLOOR)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_TC_FLOOR)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_D_FLOOR)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_HC_ROOF)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_TC_ROOF)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_D_ROOF)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_HC_ROAD)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_TC_ROAD)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_D_ROAD)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_HC_WALL)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_TC_WALL)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_D_WALL)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_QIN)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_QIN_FRAD)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_SHGC) 
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_U_WIN)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_GR)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_FLOOR_HEIGHT)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_INF)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_F_WATER_COND)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_QIN_FLAT)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_HR_TARGET)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_V_VENT)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_COP_RAT)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_GREENROOF)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_SHADE)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_NATVENT)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_EMIS_PANEL)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_ALB_PANEL)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_EFF_PANEL)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_FRAC_PANEL)
  NULLIFY(BLD_DESC_MODEL(J)%XDESC_RESIDENTIAL)
ENDDO
BLD_DESC_MODEL(:)%NDESC_BLD=0
BLD_DESC_MODEL(:)%NDESC_AGE=0
BLD_DESC_MODEL(:)%NDESC_CODE=0
BLD_DESC_MODEL(:)%NDESC_USE=0
BLD_DESC_MODEL(:)%NDESC_ROOF_LAYER=0
BLD_DESC_MODEL(:)%NDESC_ROAD_LAYER=0
BLD_DESC_MODEL(:)%NDESC_WALL_LAYER=0
BLD_DESC_MODEL(:)%NDESC_FLOOR_LAYER=0
IF (LHOOK) CALL DR_HOOK("MODD_BLD_DESCRIPTION_n:BLD_DESC_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE BLD_DESC_ALLOC
!
SUBROUTINE BLD_DESC_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_BLD_DESCRIPTION_n:BLD_DESC_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(BLD_DESC_MODEL)) DEALLOCATE(BLD_DESC_MODEL)
IF (ASSOCIATED(BLD_DESC)) NULLIFY(BLD_DESC)
IF (LHOOK) CALL DR_HOOK("MODD_BLD_DESCRIPTION_n:BLD_DESC_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE BLD_DESC_DEALLO
!
END MODULE MODD_BLD_DESCRIPTION_n
