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



CONTAINS

!




SUBROUTINE TEB_GRID_INIT(YTEB_GRID)
TYPE(TEB_GRID_t), INTENT(INOUT) :: YTEB_GRID
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GRID_N:TEB_GRID_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YTEB_GRID%XGRID_PAR)
  NULLIFY(YTEB_GRID%XLAT)
  NULLIFY(YTEB_GRID%XLON)
  NULLIFY(YTEB_GRID%XMESH_SIZE)
YTEB_GRID%NDIM=0
YTEB_GRID%CGRID=' '
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GRID_N:TEB_GRID_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_GRID_INIT


END MODULE MODD_TEB_GRID_n
