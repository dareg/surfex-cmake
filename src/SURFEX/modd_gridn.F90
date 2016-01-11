!     ##################
      MODULE MODD_GRID_n
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

TYPE GRID_t
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
  INTEGER                         :: NGRID_PAR   ! size of XGRID_PAR
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

END TYPE GRID_t



CONTAINS

!




SUBROUTINE GRID_INIT(YGRID)
TYPE(GRID_t), INTENT(INOUT) :: YGRID
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_GRID_N:IGRID_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YGRID%XGRID_PAR)
  NULLIFY(YGRID%XLAT)
  NULLIFY(YGRID%XLON)
  NULLIFY(YGRID%XMESH_SIZE)
YGRID%NDIM=0
YGRID%CGRID=' '
IF (LHOOK) CALL DR_HOOK("MODD_GRID_N:GRID_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE GRID_INIT


END MODULE MODD_GRID_n
