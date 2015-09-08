!     ##################
      MODULE MODD_ISBA_GRID_n
!     ##################
!
!!****  *MODD_ISBA - declaration of grid for ISBA scheme
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

TYPE ISBA_GRID_t
!-------------------------------------------------------------------------------
!
! Grid definition
!
  INTEGER                         :: NDIM        ! number of points
  CHARACTER(LEN=10)               :: CGRID       ! grid type
!                                              ! "NONE        " : no grid computations
!                                              ! "CONF PROJ   " : conformal projection
!                                              ! "SURF ATM    " : nature points of surf. atm. grid
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

END TYPE ISBA_GRID_t

TYPE(ISBA_GRID_t), ALLOCATABLE, TARGET, SAVE :: ISBA_GRID_MODEL(:)

TYPE(ISBA_GRID_t), POINTER :: ISBA_GRID => NULL()
!$OMP THREADPRIVATE(ISBA_GRID)

CONTAINS

SUBROUTINE ISBA_GRID_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_ISBA_GRID_N:ISBA_GRID_GOTO_MODEL',0,ZHOOK_HANDLE)

ISBA_GRID => ISBA_GRID_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_ISBA_GRID_N:ISBA_GRID_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE ISBA_GRID_GOTO_MODEL

SUBROUTINE ISBA_GRID_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_GRID_N:ISBA_GRID_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(ISBA_GRID_MODEL(KMODEL))
ISBA_GRID => ISBA_GRID_MODEL(KMODEL)
DO J=1,KMODEL
  NULLIFY(ISBA_GRID_MODEL(J)%XGRID_PAR)
  NULLIFY(ISBA_GRID_MODEL(J)%XLAT)
  NULLIFY(ISBA_GRID_MODEL(J)%XLON)
  NULLIFY(ISBA_GRID_MODEL(J)%XMESH_SIZE)
ENDDO
ISBA_GRID_MODEL(:)%NDIM=0
ISBA_GRID_MODEL(:)%CGRID=' '
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_GRID_N:ISBA_GRID_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_GRID_ALLOC

SUBROUTINE ISBA_GRID_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_GRID_N:ISBA_GRID_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(ISBA_GRID_MODEL)) DEALLOCATE(ISBA_GRID_MODEL)
IF (ASSOCIATED(ISBA_GRID)) NULLIFY(ISBA_GRID)
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_GRID_N:ISBA_GRID_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_GRID_DEALLO

END MODULE MODD_ISBA_GRID_n
