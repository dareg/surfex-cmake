!     ##################
      MODULE MODD_OCEAN_GRID
!     ##################
!
!!****  *MODD_OCEAN_GRID - declaration of grid for oceanic model
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
!!      C. Lebeaupin Brossier   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2008
!       11/2014 : NOCKMAX not parameter
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
INTEGER, SAVE :: NOCKMIN  !first ocean level indice
INTEGER, SAVE :: NOCKMAX  ! last ocean level indice
!
REAL, POINTER, DIMENSION(:) :: XK1
REAL, POINTER, DIMENSION(:) :: XK2
REAL, POINTER, DIMENSION(:) :: XK3
REAL, POINTER, DIMENSION(:) :: XK4
REAL, POINTER, DIMENSION(:) :: XZHOC
REAL, POINTER, DIMENSION(:) :: XZ2
REAL, POINTER, DIMENSION(:) :: XDZ1
REAL, POINTER, DIMENSION(:) :: XDZ2
REAL, POINTER, DIMENSION(:) :: XRAY
!
END MODULE MODD_OCEAN_GRID
