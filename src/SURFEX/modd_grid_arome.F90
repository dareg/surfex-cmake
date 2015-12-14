!     ################
      MODULE MODD_GRID_AROME
!     ################
!
!!****  *MODD_GRID_AROME - declaration of Arome grid characteristics
!!
!!    PURPOSE
!!    -------
!     Used if CINGRID_TYPE = 'AROME '
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
!!      V. Masson    *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2004
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
REAL    :: XX  ! total physical size in X direction (meters)
REAL    :: XY  ! total physical size in Y direction (meters)
INTEGER :: NX  ! dimension in X direction           (number of points)
INTEGER :: NY  ! dimension in Y direction           (number of points)
!
REAL    :: XLAT0  ! reference latitude
REAL    :: XLON0  ! reference longitude
REAL    :: XLATOR ! origin latitude
REAL    :: XLONOR ! origin longitude
REAL    :: XRPK   ! projection parameter for the conformal projection
REAL    :: XBETA  ! rotation   parameter for the conformal projection
!
REAL, DIMENSION(:), ALLOCATABLE :: XZX       ! X coordinate
REAL, DIMENSION(:), ALLOCATABLE :: XZY       ! Y coordinate
INTEGER, DIMENSION(:), ALLOCATABLE  :: NIX  ! number of points on each line
!
END MODULE MODD_GRID_AROME
