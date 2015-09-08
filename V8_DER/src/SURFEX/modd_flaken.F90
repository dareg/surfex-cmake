!     ####################
MODULE MODD_FLAKE_n
!     ####################
!
!!****  *MODD_FLAKE_n - declaration of surface parameters for the FLake model 
!!                      for inland water surfaces
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
!!      V. Masson  *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2004
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
TYPE FLAKE_t 
!
!-------------------------------------------------------------------------------------
! General surface: 
!-------------------------------------------------------------------------------------
!
  REAL, POINTER, DIMENSION(:) :: XZS       ! orography                     (m)
  REAL, POINTER, DIMENSION(:) :: XZ0       ! roughness length              (m)
  REAL, POINTER, DIMENSION(:) :: XUSTAR    ! air friction velocity         (m/s)
  REAL, POINTER, DIMENSION(:) :: XEMIS     ! water surface emissivity (NOT USED BY FLAKE)
!
  REAL, POINTER, DIMENSION(:,:) :: XCOVER  ! fraction of each ecosystem    (-)
!                                          ! F: no atmospheric layers below forcing level  
!
  LOGICAL, POINTER, DIMENSION(:) :: LCOVER ! GCOVER(i)=T --> ith cover field is not 0.
  LOGICAL                        :: LSBL   ! T: SBL scheme within the Surface Boundary Layer
!
!-------------------------------------------------------------------------------------
! Date and time:
!-------------------------------------------------------------------------------------
!
  TYPE (DATE_TIME)                  :: TTIME         ! current date and time
!
  REAL                              :: XTSTEP        ! time step
!
  REAL                              :: XOUT_TSTEP    ! output writing time step
!
!-------------------------------------------------------------------------------------
! FLake switches
!-------------------------------------------------------------------------------------
!
  LOGICAL            :: LSEDIMENTS  ! flag to use or not the bottom sediments
  LOGICAL            :: LSKINTEMP   ! flag to use or not the skin temperature computation
  CHARACTER(LEN=3)   :: CSNOW_FLK   ! FLake snow scheme
  CHARACTER(LEN=5)   :: CFLK_FLUX   ! Type of flux computation
  CHARACTER(LEN=4)   :: CFLK_ALB    ! Type of albedo
!
!-------------------------------------------------------------------------------------
! FLake parameters and variables
!-------------------------------------------------------------------------------------
!
  REAL, POINTER, DIMENSION(:) :: XWATER_DEPTH  ! Lake depth (m)
  REAL, POINTER, DIMENSION(:) :: XWATER_FETCH  ! Lake fetch (m)
  REAL, POINTER, DIMENSION(:) :: XT_BS         ! Temperature at the outer edge of the thermally 
                                               !       active layer of the bottom sediments [K]
  REAL, POINTER, DIMENSION(:) :: XDEPTH_BS     ! Depth of the thermally active layer of the
                                               !       bottom sediments [m]
  REAL, POINTER, DIMENSION(:) :: XCORIO        ! The Coriolis parameter [s^{-1}]
  REAL, POINTER, DIMENSION(:) :: XDIR_ALB      ! Water surface direct albedo
  REAL, POINTER, DIMENSION(:) :: XSCA_ALB      ! Water surface diffuse albedo
  REAL, POINTER, DIMENSION(:) :: XICE_ALB      ! Ice surface albedo (for ESM coupling)
  REAL, POINTER, DIMENSION(:) :: XSNOW_ALB     ! Snow surface albedo
  REAL, POINTER, DIMENSION(:) :: XEXTCOEF_WATER ! Extinction coefficient for the water [m^{-1}]
  REAL, POINTER, DIMENSION(:) :: XEXTCOEF_ICE   ! Extinction coefficient for the ice [m^{-1}]
  REAL, POINTER, DIMENSION(:) :: XEXTCOEF_SNOW  ! Extinction coefficient for the snow [m^{-1}] 
  REAL, POINTER, DIMENSION(:) :: XT_SNOW       ! Temperature at the air-snow interface [K]    
  REAL, POINTER, DIMENSION(:) :: XT_ICE        ! Temperature at the snow-ice or air-ice 
                                               !        interface [K]
  REAL, POINTER, DIMENSION(:) :: XT_MNW        ! Mean temperature of the water column [K]
  REAL, POINTER, DIMENSION(:) :: XT_WML        ! Mixed-layer temperature [K]
  REAL, POINTER, DIMENSION(:) :: XT_BOT        ! Temperature at the water-bottom sediment 
                                               !        interface [K]
  REAL, POINTER, DIMENSION(:) :: XT_B1         ! Temperature at the bottom of the upper 
                                               !        layer of the sediments [K]
  REAL, POINTER, DIMENSION(:) :: XCT           ! Shape factor (thermocline)
  REAL, POINTER, DIMENSION(:) :: XH_SNOW       ! Snow thickness [m]
  REAL, POINTER, DIMENSION(:) :: XH_ICE        ! Ice thickness [m]
  REAL, POINTER, DIMENSION(:) :: XH_ML         ! Thickness of the mixed-layer [m]
  REAL, POINTER, DIMENSION(:) :: XH_B1         ! Thickness of the upper layer of bottom sediments [m]                                    
!
  REAL, POINTER, DIMENSION(:) :: XTS  ! surface temperature  (K)
                                      ! (water or ice or snow)
!
!-------------------------------------------------------------------------------------
! Coupling field for Earth system model
!-------------------------------------------------------------------------------------
!
  REAL, POINTER, DIMENSION(:) :: XCPL_FLAKE_EVAP ! Evaporation for ESM coupling
  REAL, POINTER, DIMENSION(:) :: XCPL_FLAKE_RAIN ! Rainfall for ESM coupling
  REAL, POINTER, DIMENSION(:) :: XCPL_FLAKE_SNOW ! Snowfall for ESM coupling
!
END TYPE FLAKE_t
!
!-------------------------------------------------------------------------------------
!
TYPE(FLAKE_t), ALLOCATABLE, TARGET, SAVE :: FLAKE_MODEL(:)

TYPE(FLAKE_t), POINTER :: FLAKE => NULL()
!$OMP THREADPRIVATE(FLAKE)

CONTAINS
!
!-------------------------------------------------------------------------------------
!
SUBROUTINE FLAKE_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
!
IF (LHOOK) CALL DR_HOOK('MODD_FLAKE_N:FLAKE_GOTO_MODEL',0,ZHOOK_HANDLE)

FLAKE => FLAKE_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_FLAKE_N:FLAKE_GOTO_MODEL',1,ZHOOK_HANDLE)
!
END SUBROUTINE FLAKE_GOTO_MODEL
!
!-------------------------------------------------------------------------------------
!
SUBROUTINE FLAKE_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_FLAKE_N:FLAKE_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(FLAKE_MODEL(KMODEL))
FLAKE => FLAKE_MODEL(KMODEL)
DO J=1,KMODEL
  NULLIFY(FLAKE_MODEL(J)%XZS)
  NULLIFY(FLAKE_MODEL(J)%XZ0)
  NULLIFY(FLAKE_MODEL(J)%XUSTAR)
  NULLIFY(FLAKE_MODEL(J)%XCOVER)
  NULLIFY(FLAKE_MODEL(J)%LCOVER)
  NULLIFY(FLAKE_MODEL(J)%XEMIS)
  NULLIFY(FLAKE_MODEL(J)%XWATER_DEPTH)
  NULLIFY(FLAKE_MODEL(J)%XWATER_FETCH)
  NULLIFY(FLAKE_MODEL(J)%XT_BS)
  NULLIFY(FLAKE_MODEL(J)%XDEPTH_BS)
  NULLIFY(FLAKE_MODEL(J)%XCORIO)
  NULLIFY(FLAKE_MODEL(J)%XDIR_ALB)
  NULLIFY(FLAKE_MODEL(J)%XSCA_ALB)
  NULLIFY(FLAKE_MODEL(J)%XICE_ALB)
  NULLIFY(FLAKE_MODEL(J)%XSNOW_ALB)
  NULLIFY(FLAKE_MODEL(J)%XEXTCOEF_WATER)
  NULLIFY(FLAKE_MODEL(J)%XEXTCOEF_ICE)
  NULLIFY(FLAKE_MODEL(J)%XEXTCOEF_SNOW)
  NULLIFY(FLAKE_MODEL(J)%XT_SNOW)
  NULLIFY(FLAKE_MODEL(J)%XT_ICE)
  NULLIFY(FLAKE_MODEL(J)%XT_MNW)
  NULLIFY(FLAKE_MODEL(J)%XT_WML)
  NULLIFY(FLAKE_MODEL(J)%XT_BOT)
  NULLIFY(FLAKE_MODEL(J)%XT_B1)
  NULLIFY(FLAKE_MODEL(J)%XCT)
  NULLIFY(FLAKE_MODEL(J)%XH_SNOW)
  NULLIFY(FLAKE_MODEL(J)%XH_ICE)
  NULLIFY(FLAKE_MODEL(J)%XH_ML)
  NULLIFY(FLAKE_MODEL(J)%XH_B1)
  NULLIFY(FLAKE_MODEL(J)%XTS)
  NULLIFY(FLAKE_MODEL(J)%XCPL_FLAKE_EVAP)
  NULLIFY(FLAKE_MODEL(J)%XCPL_FLAKE_RAIN)
  NULLIFY(FLAKE_MODEL(J)%XCPL_FLAKE_SNOW)
ENDDO
FLAKE_MODEL(:)%LSBL=.FALSE.
FLAKE_MODEL(:)%XTSTEP=0.
FLAKE_MODEL(:)%XOUT_TSTEP=0.
FLAKE_MODEL(:)%LSEDIMENTS=.FALSE.
FLAKE_MODEL(:)%LSKINTEMP=.FALSE.
FLAKE_MODEL(:)%CSNOW_FLK='   '
FLAKE_MODEL(:)%CFLK_ALB='    '
FLAKE_MODEL(:)%CFLK_FLUX='     '
IF (LHOOK) CALL DR_HOOK("MODD_FLAKE_N:FLAKE_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE FLAKE_ALLOC
!
!-------------------------------------------------------------------------------------
!
SUBROUTINE FLAKE_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_FLAKE_N:FLAKE_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(FLAKE_MODEL)) DEALLOCATE(FLAKE_MODEL)
IF (ASSOCIATED(FLAKE)) NULLIFY(FLAKE)
IF (LHOOK) CALL DR_HOOK("MODD_FLAKE_N:FLAKE_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE FLAKE_DEALLO
!
!-------------------------------------------------------------------------------------
!
END MODULE MODD_FLAKE_n
