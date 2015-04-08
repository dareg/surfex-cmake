!##################
MODULE MODD_TEB_GARDEN_OPTION_n
!##################
!
!!****  *MODD_TEB_GARDEN - declaration of packed surface parameters for ISBA scheme
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
!!      A. Lemonsu   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2011
!!      V. Masson      06/2013 splits module in two
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TYPE_SNOW
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE TEB_GARDEN_OPTIONS_t
!-------------------------------------------------------------------------------
!
! type of initialization of vegetation: from cover types (ecoclimap) or parameters prescribed
!
  LOGICAL                              :: LPAR_GARDEN      ! T: parameters computed from ecoclimap
!                                                          ! F: they are read in the file
!
! Number of inside garden vegetation (not TEB) patches and of layers
! 
!
  INTEGER                              :: NGROUND_LAYER    ! number of ground layers
!
  INTEGER                              :: NLAYER_HORT
  INTEGER                              :: NLAYER_DUN
!
  REAL, POINTER, DIMENSION(:)          :: XSOILGRID        ! Soil layer grid as reference for DIF
!
END TYPE TEB_GARDEN_OPTIONS_t
!-------------------------------------------------------------------------------

TYPE(TEB_GARDEN_OPTIONS_t), ALLOCATABLE, TARGET, SAVE :: TEB_GARDEN_OPTIONS_MODEL(:)

TYPE(TEB_GARDEN_OPTIONS_t), POINTER :: TEB_GARDEN_OPTIONS => NULL()
!$OMP THREADPRIVATE(TEB_GARDEN_OPTIONS)

CONTAINS

SUBROUTINE TEB_GARDEN_OPTIONS_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_TEB_GARDEN_N:TEB_GARDEN_OPTIONS_GOTO_MODEL',0,ZHOOK_HANDLE)

TEB_GARDEN_OPTIONS => TEB_GARDEN_OPTIONS_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_TEB_GARDEN_N:TEB_GARDEN_OPTIONS_GOTO_MODEL',1,ZHOOK_HANDLE)
!
END SUBROUTINE TEB_GARDEN_OPTIONS_GOTO_MODEL

SUBROUTINE TEB_GARDEN_OPTIONS_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GARDEN_N:TEB_GARDEN_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(TEB_GARDEN_OPTIONS_MODEL(KMODEL))
DO J=1,KMODEL
  NULLIFY(TEB_GARDEN_OPTIONS_MODEL(J)%XSOILGRID)
ENDDO
TEB_GARDEN_OPTIONS_MODEL(:)%LPAR_GARDEN=.TRUE.
TEB_GARDEN_OPTIONS_MODEL(:)%NGROUND_LAYER=0
TEB_GARDEN_OPTIONS_MODEL(:)%NLAYER_HORT=0
TEB_GARDEN_OPTIONS_MODEL(:)%NLAYER_DUN=0
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GARDEN_N:TEB_GARDEN_OPTIONS_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_GARDEN_OPTIONS_ALLOC

SUBROUTINE TEB_GARDEN_OPTIONS_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GARDEN_N:TEB_GARDEN_OPTIONS_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(TEB_GARDEN_OPTIONS_MODEL)) DEALLOCATE(TEB_GARDEN_OPTIONS_MODEL)
IF (ASSOCIATED(TEB_GARDEN_OPTIONS)) NULLIFY(TEB_GARDEN_OPTIONS)
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GARDEN_N:TEB_GARDEN_OPTIONS_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_GARDEN_OPTIONS_DEALLO

END MODULE MODD_TEB_GARDEN_OPTION_n
