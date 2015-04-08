!     ################
      MODULE MODD_TEB_IRRIG_n
!     ################
!
!!****  *MODD_TEB_IRRIG_n - declaration of surface parameters for urban canopy
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
!!      Original       07/2006
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE TEB_IRRIG_t
!
  LOGICAL                        :: LTEB_IRRIG        ! flag to use irrigation for gardens or greenroofs
  LOGICAL                        :: LPAR_GD_IRRIG     ! flag to use prescribed irrigation for gardens
  LOGICAL                        :: LPAR_GR_IRRIG     ! flag to use prescribed irrigation for greenroofs
  LOGICAL                        :: LPAR_RD_IRRIG     ! flag to use prescribed irrigation for roads
  REAL,    POINTER, DIMENSION(:) :: XGD_START_MONTH   ! gardens : start month for irrigation (included)
  REAL,    POINTER, DIMENSION(:) :: XGD_END_MONTH     ! gardens : end   month for irrigation (included)
  REAL,    POINTER, DIMENSION(:) :: XGD_START_HOUR    ! gardens : start solar hour for irrigation (included, hour)
  REAL,    POINTER, DIMENSION(:) :: XGD_END_HOUR      ! gardens : end   solar hour for irrigation (excluded, hour)
  REAL,    POINTER, DIMENSION(:) :: XGD_24H_IRRIG     ! gardens : total irrigation over 24 hours (kg/m2)
  REAL,    POINTER, DIMENSION(:) :: XGR_START_MONTH   ! greenroofs : start month for irrigation (included)
  REAL,    POINTER, DIMENSION(:) :: XGR_END_MONTH     ! greenroofs : end   month for irrigation (included)
  REAL,    POINTER, DIMENSION(:) :: XGR_START_HOUR    ! greenroofs : start solar hour for irrigation (included, hour)
  REAL,    POINTER, DIMENSION(:) :: XGR_END_HOUR      ! greenroofs : end   solar hour for irrigation (excluded, hour)
  REAL,    POINTER, DIMENSION(:) :: XGR_24H_IRRIG     ! greenroofs : total irrigation over 24 hours (kg/m2)
  REAL,    POINTER, DIMENSION(:) :: XRD_START_MONTH   ! roads : start month for irrigation (included)
  REAL,    POINTER, DIMENSION(:) :: XRD_END_MONTH     ! roads : end   month for irrigation (included)
  REAL,    POINTER, DIMENSION(:) :: XRD_START_HOUR    ! roads : start solar hour for irrigation (included, hour)
  REAL,    POINTER, DIMENSION(:) :: XRD_END_HOUR      ! roads : end   solar hour for irrigation (excluded, hour)
  REAL,    POINTER, DIMENSION(:) :: XRD_24H_IRRIG     ! roads : total irrigation over 24 hours (kg/m2)
!
END TYPE TEB_IRRIG_t

TYPE(TEB_IRRIG_t), ALLOCATABLE, TARGET, SAVE :: TEB_IRRIG_MODEL(:)

TYPE(TEB_IRRIG_t), POINTER :: TEB_IRRIG => NULL()
!$OMP THREADPRIVATE(TEB_IRRIG)

CONTAINS

SUBROUTINE TEB_IRRIG_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_TEB_IRRIG_N:TEB_IRRIG_GOTO_MODEL',0,ZHOOK_HANDLE)

TEB_IRRIG => TEB_IRRIG_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_TEB_IRRIG_N:TEB_IRRIG_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE TEB_IRRIG_GOTO_MODEL

SUBROUTINE TEB_IRRIG_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_IRRIG_N:TEB_IRRIG_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(TEB_IRRIG_MODEL(KMODEL))
DO J=1,KMODEL
  NULLIFY(TEB_IRRIG_MODEL(J)%XGD_START_MONTH)
  NULLIFY(TEB_IRRIG_MODEL(J)%XGD_END_MONTH)
  NULLIFY(TEB_IRRIG_MODEL(J)%XGD_START_HOUR)
  NULLIFY(TEB_IRRIG_MODEL(J)%XGD_END_HOUR)
  NULLIFY(TEB_IRRIG_MODEL(J)%XGD_24H_IRRIG)
  NULLIFY(TEB_IRRIG_MODEL(J)%XGR_START_MONTH)
  NULLIFY(TEB_IRRIG_MODEL(J)%XGR_END_MONTH)
  NULLIFY(TEB_IRRIG_MODEL(J)%XGR_START_HOUR)
  NULLIFY(TEB_IRRIG_MODEL(J)%XGR_END_HOUR)
  NULLIFY(TEB_IRRIG_MODEL(J)%XGR_24H_IRRIG)
  NULLIFY(TEB_IRRIG_MODEL(J)%XRD_START_MONTH)
  NULLIFY(TEB_IRRIG_MODEL(J)%XRD_END_MONTH)
  NULLIFY(TEB_IRRIG_MODEL(J)%XRD_START_HOUR)
  NULLIFY(TEB_IRRIG_MODEL(J)%XRD_END_HOUR)
  NULLIFY(TEB_IRRIG_MODEL(J)%XRD_24H_IRRIG)
ENDDO
TEB_IRRIG_MODEL(:)%LTEB_IRRIG = .FALSE.
TEB_IRRIG_MODEL(:)%LPAR_GD_IRRIG = .FALSE.
TEB_IRRIG_MODEL(:)%LPAR_GR_IRRIG = .FALSE.
TEB_IRRIG_MODEL(:)%LPAR_RD_IRRIG = .FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_TEB_IRRIG_N:TEB_IRRIG_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_IRRIG_ALLOC

SUBROUTINE TEB_IRRIG_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_IRRIG_N:TEB_IRRIG_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(TEB_IRRIG_MODEL)) DEALLOCATE(TEB_IRRIG_MODEL)
IF (ASSOCIATED(TEB_IRRIG)) NULLIFY(TEB_IRRIG)
IF (LHOOK) CALL DR_HOOK("MODD_TEB_IRRIG_N:TEB_IRRIG_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_IRRIG_DEALLO

END MODULE MODD_TEB_IRRIG_n
