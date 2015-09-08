!     ##################
MODULE MODD_FLAKE_GRID_n
!     ##################
!
!!****  *MODD_FLAKE_GRID - declaration of FLAKE grid
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

TYPE FLAKE_GRID_t
!-------------------------------------------------------------------------------
!
! Grid definition
!
  INTEGER                         :: NDIM        ! number of points
  CHARACTER(LEN=10)               :: CGRID       ! grid type
!                                              ! "NONE        " : no grid computations
!                                              ! "CONF PROJ   " : conformal projection
!                                              ! "SURF ATM    " : water points of surf. atm. grid
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

END TYPE FLAKE_GRID_t

TYPE(FLAKE_GRID_t), ALLOCATABLE, TARGET, SAVE :: FLAKE_GRID_MODEL(:)

TYPE(FLAKE_GRID_t), POINTER :: FLAKE_GRID => NULL()
!$OMP THREADPRIVATE(FLAKE_GRID)

CONTAINS

SUBROUTINE FLAKE_GRID_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_FLAKE_GRID_N:FLAKE_GRID_GOTO_MODEL',0,ZHOOK_HANDLE)

FLAKE_GRID => FLAKE_GRID_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_FLAKE_GRID_N:FLAKE_GRID_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE FLAKE_GRID_GOTO_MODEL

SUBROUTINE FLAKE_GRID_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_FLAKE_GRID_N:FLAKE_GRID_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(FLAKE_GRID_MODEL(KMODEL))
FLAKE_GRID => FLAKE_GRID_MODEL(KMODEL)
DO J=1,KMODEL
  NULLIFY(FLAKE_GRID_MODEL(J)%XGRID_PAR)
  NULLIFY(FLAKE_GRID_MODEL(J)%XLAT)
  NULLIFY(FLAKE_GRID_MODEL(J)%XLON)
  NULLIFY(FLAKE_GRID_MODEL(J)%XMESH_SIZE)
ENDDO
FLAKE_GRID_MODEL(:)%NDIM=0
FLAKE_GRID_MODEL(:)%CGRID=' '
IF (LHOOK) CALL DR_HOOK("MODD_FLAKE_GRID_N:FLAKE_GRID_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE FLAKE_GRID_ALLOC

SUBROUTINE FLAKE_GRID_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_FLAKE_GRID_N:FLAKE_GRID_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(FLAKE_GRID_MODEL)) DEALLOCATE(FLAKE_GRID_MODEL)
IF (ASSOCIATED(FLAKE_GRID)) NULLIFY(FLAKE_GRID)
IF (LHOOK) CALL DR_HOOK("MODD_FLAKE_GRID_N:FLAKE_GRID_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE FLAKE_GRID_DEALLO

END MODULE MODD_FLAKE_GRID_n
