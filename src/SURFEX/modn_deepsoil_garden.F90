!     ##################
      MODULE MODN_DEEPSOIL_GARDEN
!     ##################
!
!!****  *MODN_DEEPSOIL_GARDEN - deep soil characteristics
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
!!	P. Le Moigne   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original   05/2008
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_DEEPSOIL_GARDEN, ONLY:  LDEEPSOIL, LPHYSDOMC
!
IMPLICIT NONE
!
NAMELIST/NAM_DEEPSOIL/LDEEPSOIL, LPHYSDOMC
!
END MODULE MODN_DEEPSOIL_GARDEN
