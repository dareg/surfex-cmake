!     ################
      MODULE MODD_HORIBL
!     ################
!
!!****  *MODD_HORIBL - declaration for field interpolations
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
LOGICAL :: LGLOBLON=.FALSE., LGLOBS=.FALSE., LGLOBN=.FALSE.
REAL :: XILO1H=0., XILO2H=0.

INTEGER, DIMENSION(:,:), ALLOCATABLE :: NO
INTEGER, DIMENSION(:), ALLOCATABLE :: NINLOH
REAL, DIMENSION(:,:), ALLOCATABLE :: XLA
REAL, DIMENSION(:), ALLOCATABLE :: XOLA, XOLO
INTEGER, DIMENSION(:,:), ALLOCATABLE :: NP
REAL, DIMENSION(:,:), ALLOCATABLE :: XLOPH
!
END MODULE MODD_HORIBL
