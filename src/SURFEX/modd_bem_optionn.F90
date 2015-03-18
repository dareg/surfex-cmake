!     ################
      MODULE MODD_BEM_OPTION_n
!     ################
!
!!****  *MODD_BEM_n - declaration of parameters and option for BEM
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
!!      B. Bueno   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       10/2010
!!      G. Pigeon      06/2011 add LSHAD_DAY
!!      G. Pigeon      07/2011 add LNATVENT_NIGHT
!!      G. Pigeon      08/2011 change from MODD_BLD -> MODD_BEM
!!      G. Pigeon      10/2011 add indoor relative surf. and view factors
!!      G. Pigeon      09/2012 add TRAN_WIN
!!      G. Pigeon      10/2012 add XF_WIN_WIN
!!      V. Masson      06/2013 splits module in two
!
!*       0.   DECLARATIONS
!             ------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE BEM_OPTIONS_t
! BLD scheme option
!
! Number of layers
!
  INTEGER                       :: NFLOOR_LAYER   ! number of layers in walls
  CHARACTER(LEN=6)              :: CCOOL_COIL    ! type of cooling coil
  CHARACTER(LEN=6)              :: CHEAT_COIL    ! type of heating coil
  LOGICAL                       :: LAUTOSIZE     ! Flag to activate autosize calculations
!
END TYPE BEM_OPTIONS_t
!
TYPE(BEM_OPTIONS_t), ALLOCATABLE, TARGET, SAVE :: BEM_OPTIONS_MODEL(:)
!
INTEGER, POINTER :: NFLOOR_LAYER=>NULL()
!$OMP THREADPRIVATE(NFLOOR_LAYER)
 CHARACTER(LEN=6), POINTER     :: CCOOL_COIL=>NULL()
!$OMP THREADPRIVATE(CCOOL_COIL)
 CHARACTER(LEN=6), POINTER     :: CHEAT_COIL=>NULL()
!$OMP THREADPRIVATE(CHEAT_COIL)
LOGICAL, POINTER :: LAUTOSIZE=>NULL()
!$OMP THREADPRIVATE(LAUTOSIZE)
!
!
!--------------------------------------------------------------------------
!
CONTAINS

SUBROUTINE BEM_OPTIONS_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Save current state for allocated arrays
IF (LKFROM) THEN
ENDIF
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_BEM_N:BEM_OPTIONS_GOTO_MODEL',0,ZHOOK_HANDLE)
NFLOOR_LAYER=>BEM_OPTIONS_MODEL(KTO)%NFLOOR_LAYER
CCOOL_COIL=>BEM_OPTIONS_MODEL(KTO)%CCOOL_COIL
CHEAT_COIL=>BEM_OPTIONS_MODEL(KTO)%CHEAT_COIL
LAUTOSIZE=>BEM_OPTIONS_MODEL(KTO)%LAUTOSIZE
IF (LHOOK) CALL DR_HOOK('MODD_BEM_N:BEM_OPTIONS_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE BEM_OPTIONS_GOTO_MODEL

SUBROUTINE BEM_OPTIONS_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_BEM_N:BEM_OPTIONS_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(BEM_OPTIONS_MODEL(KMODEL))
BEM_OPTIONS_MODEL(:)%NFLOOR_LAYER = 0
BEM_OPTIONS_MODEL(:)%CCOOL_COIL   = '      '
BEM_OPTIONS_MODEL(:)%CHEAT_COIL   = '      '
BEM_OPTIONS_MODEL(:)%LAUTOSIZE    = .FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_BEM_N:BEM_OPTIONS_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE BEM_OPTIONS_ALLOC

SUBROUTINE BEM_OPTIONS_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_BEM_N:BEM_OPTIONS_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(BEM_OPTIONS_MODEL)) DEALLOCATE(BEM_OPTIONS_MODEL)
IF (LHOOK) CALL DR_HOOK("MODD_BEM_N:BEM_OPTIONS_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE BEM_OPTIONS_DEALLO

!----------------------------------------------------------------------------
!
END MODULE MODD_BEM_OPTION_n
