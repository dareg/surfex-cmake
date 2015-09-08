!     ##################
      MODULE MODD_SEAFLUX_GRID_n
!     ##################
!
!!****  *MODD_SEAFLUX_GRID - declaration of SEAFLUX grid
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
!!      V. Masson  *Meteo France*
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

TYPE SEAFLUX_GRID_t
!-------------------------------------------------------------------------------
!
! Grid definition
!
  INTEGER                         :: NDIM        ! number of points
  CHARACTER(LEN=10)               :: CGRID       ! grid type
!                                              ! "NONE        " : no grid computations
!                                              ! "CONF PROJ   " : conformal projection
!                                              ! "SURF ATM    " : SEA points of surf. atm. grid
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

END TYPE SEAFLUX_GRID_t



CONTAINS

!




SUBROUTINE SEAFLUX_GRID_INIT(YSEAFLUX_GRID)
TYPE(SEAFLUX_GRID_t), INTENT(INOUT) :: YSEAFLUX_GRID
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_SEAFLUX_GRID_N:SEAFLUX_GRID_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YSEAFLUX_GRID%XGRID_PAR)
  NULLIFY(YSEAFLUX_GRID%XLAT)
  NULLIFY(YSEAFLUX_GRID%XLON)
  NULLIFY(YSEAFLUX_GRID%XMESH_SIZE)
YSEAFLUX_GRID%NDIM=0
YSEAFLUX_GRID%CGRID=' '
IF (LHOOK) CALL DR_HOOK("MODD_SEAFLUX_GRID_N:SEAFLUX_GRID_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE SEAFLUX_GRID_INIT


END MODULE MODD_SEAFLUX_GRID_n
