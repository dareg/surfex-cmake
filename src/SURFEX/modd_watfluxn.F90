!     ####################
      MODULE MODD_WATFLUX_n
!     ####################
!
!!****  *MODD_WATFLUX_n - declaration of surface parameters for an inland water surface
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
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE WATFLUX_t
!
! General surface: 
!
  REAL, POINTER, DIMENSION(:)   :: XZS     ! orography                     (m)
  REAL, POINTER, DIMENSION(:,:) :: XCOVER  ! fraction of each ecosystem    (-)
  LOGICAL, POINTER, DIMENSION(:):: LCOVER  ! GCOVER(i)=T --> ith cover field is not 0.
  LOGICAL                       :: LSBL    ! T: SBL scheme within the Surface Boundary Layer
!                                          ! F: no atmospheric layers below forcing level
  CHARACTER(LEN=4)              :: CWAT_ALB ! type of albedo
!
  LOGICAL                       :: LINTERPOL_TS ! Interpotalation of monthly TS
  CHARACTER(LEN=6)              :: CINTERPOL_TS ! Interpotalation of monthly TS
!
! Inland water:
!
  REAL, POINTER, DIMENSION(:) :: XTS       ! water surface temperature               (K)
  REAL, POINTER, DIMENSION(:) :: XTICE     ! water ice temperature
  REAL, POINTER, DIMENSION(:) :: XZ0       ! water surface roughness length          (-)
  REAL, POINTER, DIMENSION(:) :: XEMIS     ! water surface emissivity                (-)
  REAL, POINTER, DIMENSION(:) :: XDIR_ALB  ! water surface direct albedo             (-)
  REAL, POINTER, DIMENSION(:) :: XSCA_ALB  ! water surface diffuse albedo            (-)
  REAL, POINTER, DIMENSION(:) :: XICE_ALB  ! water ice albedo (for ESM coupling)     (-)
!
  REAL, POINTER, DIMENSION(:,:) :: XTS_MTH   ! Monthly water surface temperature               (K)
!
  REAL, POINTER, DIMENSION(:) :: XCPL_WATER_WIND ! 10m wind speed for ESM coupling
  REAL, POINTER, DIMENSION(:) :: XCPL_WATER_FWSU ! zonal wind stress for ESM coupling
  REAL, POINTER, DIMENSION(:) :: XCPL_WATER_FWSV ! meridian wind stress for ESM coupling
  REAL, POINTER, DIMENSION(:) :: XCPL_WATER_SNET ! Solar net heat flux
  REAL, POINTER, DIMENSION(:) :: XCPL_WATER_HEAT ! Non solar net heat flux
  REAL, POINTER, DIMENSION(:) :: XCPL_WATER_EVAP ! Evaporation for ESM coupling
  REAL, POINTER, DIMENSION(:) :: XCPL_WATER_RAIN ! Rainfall for ESM coupling
  REAL, POINTER, DIMENSION(:) :: XCPL_WATER_SNOW ! Snowfall for ESM coupling
  REAL, POINTER, DIMENSION(:) :: XCPL_WATER_FWSM ! wind stress module for ESM coupling
!
  REAL, POINTER, DIMENSION(:) :: XCPL_WATERICE_SNET ! solar net heat flux
  REAL, POINTER, DIMENSION(:) :: XCPL_WATERICE_HEAT ! Non solar net heat flux
  REAL, POINTER, DIMENSION(:) :: XCPL_WATERICE_EVAP ! Sublimation for ESM coupling
!
! Date:
!
  TYPE (DATE_TIME)                  :: TTIME         ! current date and time
  TYPE (DATE_TIME)                  :: TZTIME  
!
! Time-step:
!
  REAL                              :: XTSTEP        ! time step
!
  REAL                              :: XOUT_TSTEP    ! output writing time step
!
!
END TYPE WATFLUX_t

TYPE(WATFLUX_t), ALLOCATABLE, TARGET, SAVE :: WATFLUX_MODEL(:)

TYPE(WATFLUX_t), POINTER :: WATFLUX => NULL()
!$OMP THREADPRIVATE(WATFLUX)

CONTAINS

SUBROUTINE WATFLUX_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_WATFLUX_N:WATFLUX_GOTO_MODEL',0,ZHOOK_HANDLE)

WATFLUX => WATFLUX_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_WATFLUX_N:WATFLUX_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE WATFLUX_GOTO_MODEL

SUBROUTINE WATFLUX_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_WATFLUX_N:WATFLUX_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(WATFLUX_MODEL(KMODEL))
WATFLUX => WATFLUX_MODEL(KMODEL)
DO J=1,KMODEL
  NULLIFY(WATFLUX_MODEL(J)%XZS)
  NULLIFY(WATFLUX_MODEL(J)%XCOVER)
  NULLIFY(WATFLUX_MODEL(J)%LCOVER)
  NULLIFY(WATFLUX_MODEL(J)%XTS)
  NULLIFY(WATFLUX_MODEL(J)%XTICE)
  NULLIFY(WATFLUX_MODEL(J)%XZ0)
  NULLIFY(WATFLUX_MODEL(J)%XEMIS)
  NULLIFY(WATFLUX_MODEL(J)%XDIR_ALB)
  NULLIFY(WATFLUX_MODEL(J)%XSCA_ALB)
  NULLIFY(WATFLUX_MODEL(J)%XICE_ALB)
  NULLIFY(WATFLUX_MODEL(J)%XTS_MTH)
  NULLIFY(WATFLUX_MODEL(J)%XCPL_WATER_WIND)
  NULLIFY(WATFLUX_MODEL(J)%XCPL_WATER_FWSU)
  NULLIFY(WATFLUX_MODEL(J)%XCPL_WATER_FWSV)
  NULLIFY(WATFLUX_MODEL(J)%XCPL_WATER_SNET)
  NULLIFY(WATFLUX_MODEL(J)%XCPL_WATER_HEAT)
  NULLIFY(WATFLUX_MODEL(J)%XCPL_WATER_EVAP)
  NULLIFY(WATFLUX_MODEL(J)%XCPL_WATER_RAIN)
  NULLIFY(WATFLUX_MODEL(J)%XCPL_WATER_SNOW)
  NULLIFY(WATFLUX_MODEL(J)%XCPL_WATER_FWSM)
  NULLIFY(WATFLUX_MODEL(J)%XCPL_WATERICE_SNET)
  NULLIFY(WATFLUX_MODEL(J)%XCPL_WATERICE_HEAT)
  NULLIFY(WATFLUX_MODEL(J)%XCPL_WATERICE_EVAP)
ENDDO
WATFLUX_MODEL(:)%LSBL=.FALSE.
WATFLUX_MODEL(:)%CWAT_ALB=' '
WATFLUX_MODEL(:)%LINTERPOL_TS=.FALSE.
WATFLUX_MODEL(:)%CINTERPOL_TS=' '
WATFLUX_MODEL(:)%XTSTEP=0.
WATFLUX_MODEL(:)%XOUT_TSTEP=0.
IF (LHOOK) CALL DR_HOOK("MODD_WATFLUX_N:WATFLUX_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE WATFLUX_ALLOC

SUBROUTINE WATFLUX_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_WATFLUX_N:WATFLUX_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(WATFLUX_MODEL)) DEALLOCATE(WATFLUX_MODEL)
IF (ASSOCIATED(WATFLUX)) NULLIFY(WATFLUX)
IF (LHOOK) CALL DR_HOOK("MODD_WATFLUX_N:WATFLUX_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE WATFLUX_DEALLO

END MODULE MODD_WATFLUX_n
