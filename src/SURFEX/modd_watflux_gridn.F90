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



CONTAINS

!




SUBROUTINE WATFLUX_GRID_INIT(YWATFLUX_GRID)
TYPE(WATFLUX_GRID_t), INTENT(INOUT) :: YWATFLUX_GRID
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_WATFLUX_GRID_N:WATFLUX_GRID_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YWATFLUX_GRID%XGRID_PAR)
  NULLIFY(YWATFLUX_GRID%XLAT)
  NULLIFY(YWATFLUX_GRID%XLON)
  NULLIFY(YWATFLUX_GRID%XMESH_SIZE)
YWATFLUX_GRID%NDIM=0
YWATFLUX_GRID%CGRID=' '
IF (LHOOK) CALL DR_HOOK("MODD_WATFLUX_GRID_N:WATFLUX_GRID_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE WATFLUX_GRID_INIT


END MODULE MODD_WATFLUX_GRID_n
