!     ##################
      MODULE MODD_WATFLUX_GRID_n
!     ##################
!
!!****  *MODD_WATFLUX_GRID - declaration of WATFLUX grid
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

TYPE WATFLUX_GRID_t
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

END TYPE WATFLUX_GRID_t

TYPE(WATFLUX_GRID_t), ALLOCATABLE, TARGET, SAVE :: WATFLUX_GRID_MODEL(:)

TYPE(WATFLUX_GRID_t), POINTER :: WATFLUX_GRID => NULL()
!$OMP THREADPRIVATE(WATFLUX_GRID)

CONTAINS

SUBROUTINE WATFLUX_GRID_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_WATFLUX_GRID_N:WATFLUX_GRID_GOTO_MODEL',0,ZHOOK_HANDLE)

WATFLUX_GRID => WATFLUX_GRID_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_WATFLUX_GRID_N:WATFLUX_GRID_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE WATFLUX_GRID_GOTO_MODEL

SUBROUTINE WATFLUX_GRID_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_WATFLUX_GRID_N:WATFLUX_GRID_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(WATFLUX_GRID_MODEL(KMODEL))
DO J=1,KMODEL
  NULLIFY(WATFLUX_GRID_MODEL(J)%XGRID_PAR)
  NULLIFY(WATFLUX_GRID_MODEL(J)%XLAT)
  NULLIFY(WATFLUX_GRID_MODEL(J)%XLON)
  NULLIFY(WATFLUX_GRID_MODEL(J)%XMESH_SIZE)
ENDDO
WATFLUX_GRID_MODEL(:)%NDIM=0
WATFLUX_GRID_MODEL(:)%CGRID=' '
IF (LHOOK) CALL DR_HOOK("MODD_WATFLUX_GRID_N:WATFLUX_GRID_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE WATFLUX_GRID_ALLOC

SUBROUTINE WATFLUX_GRID_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_WATFLUX_GRID_N:WATFLUX_GRID_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(WATFLUX_GRID_MODEL)) DEALLOCATE(WATFLUX_GRID_MODEL)
IF (ASSOCIATED(WATFLUX_GRID)) NULLIFY(WATFLUX_GRID)
IF (LHOOK) CALL DR_HOOK("MODD_WATFLUX_GRID_N:WATFLUX_GRID_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE WATFLUX_GRID_DEALLO

END MODULE MODD_WATFLUX_GRID_n
