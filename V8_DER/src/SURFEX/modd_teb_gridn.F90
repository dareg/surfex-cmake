!     ##################
      MODULE MODD_TEB_GRID_n
!     ##################
!
!!****  *MODD_TEB_GRID - declaration of TEB grid
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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2004
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE TEB_GRID_t
!-------------------------------------------------------------------------------
!
! Grid definition
!
  INTEGER                         :: NDIM        ! number of points
  CHARACTER(LEN=10)               :: CGRID       ! grid type
!                                              ! "NONE        " : no grid computations
!                                              ! "CONF PROJ   " : conformal projection
!                                              ! "SURF ATM    " : town points of surf. atm. grid
!
  REAL, POINTER,     DIMENSION(:) :: XGRID_PAR   ! lits of parameters used to define the grid
!                                              ! (depends on value of CGRID)
!
!-------------------------------------------------------------------------------
!
! General surface parameters:
!
  REAL, POINTER, DIMENSION(:) :: XLAT        ! latitude (degrees +North)               (-)
  REAL, POINTER, DIMENSION(:) :: XLON        ! longitude (degrees +East)               (-)
  REAL, POINTER, DIMENSION(:) :: XMESH_SIZE  ! mesh size                               (m2)
!-------------------------------------------------------------------------------
!

END TYPE TEB_GRID_t

TYPE(TEB_GRID_t), ALLOCATABLE, TARGET, SAVE :: TEB_GRID_MODEL(:)

TYPE(TEB_GRID_t), POINTER :: TEB_GRID => NULL()
!$OMP THREADPRIVATE(TEB_GRID)

CONTAINS

SUBROUTINE TEB_GRID_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_TEB_GRID_N:TEB_GRID_GOTO_MODEL',0,ZHOOK_HANDLE)

TEB_GRID => TEB_GRID_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_TEB_GRID_N:TEB_GRID_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE TEB_GRID_GOTO_MODEL

SUBROUTINE TEB_GRID_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GRID_N:TEB_GRID_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(TEB_GRID_MODEL(KMODEL))
TEB_GRID => TEB_GRID_MODEL(KMODEL)
DO J=1,KMODEL
  NULLIFY(TEB_GRID_MODEL(J)%XGRID_PAR)
  NULLIFY(TEB_GRID_MODEL(J)%XLAT)
  NULLIFY(TEB_GRID_MODEL(J)%XLON)
  NULLIFY(TEB_GRID_MODEL(J)%XMESH_SIZE)
ENDDO
TEB_GRID_MODEL(:)%NDIM=0
TEB_GRID_MODEL(:)%CGRID=' '
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GRID_N:TEB_GRID_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_GRID_ALLOC

SUBROUTINE TEB_GRID_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GRID_N:TEB_GRID_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(TEB_GRID_MODEL)) DEALLOCATE(TEB_GRID_MODEL)
IF (ASSOCIATED(TEB_GRID)) NULLIFY(TEB_GRID)
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GRID_N:TEB_GRID_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_GRID_DEALLO

END MODULE MODD_TEB_GRID_n
