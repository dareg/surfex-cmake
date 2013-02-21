!     ####################
      MODULE MODD_SELECT
!     ####################
!
!!****  *MODD_SELECT - declaration of surface ATM
!!
!!    PURPOSE
!!    -------
!     Declaration of surface parameters
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
!!	P. Le Moigne *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       3/2009
!
!*       0.   DECLARATIONS
!             ------------
!
!
IMPLICIT NONE
!
!-----------------------------------------------------------------------------------------------------
LOGICAL    :: LSELECT = .FALSE.
              ! activates output selection from namelist
!
 CHARACTER(LEN=12), DIMENSION(200)    :: CNAME_SELECT
              ! name of output fields in namelist
!
LOGICAL    :: LSELECT_USER
 CHARACTER(LEN=12), DIMENSION(:), POINTER    :: CNAME_USER
!
!-----------------------------------------------------------------------------------------------------
!
END MODULE MODD_SELECT
