!####################
MODULE MODD_SURF_PAR
!####################
!
!!****  *MODD_SURF_PAR - declaration of surface parameters
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
!!	V. Masson *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       02/2004
!
!*       0.   DECLARATIONS
!             ------------
!
!
IMPLICIT NONE
!
!-----------------------------------------------------------------------------------------------------
INTEGER :: NVERSION  ! surface version
INTEGER :: NBUGFIX   ! bugfix number of this version
!
REAL,    PARAMETER :: XUNDEF = 1.E+20! undefined value
INTEGER, PARAMETER :: NUNDEF = 1E+9  ! undefined value
!-----------------------------------------------------------------------------------------------------
!
END MODULE MODD_SURF_PAR
