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
!$OMP THREADPRIVATE(NO)
INTEGER, DIMENSION(:), ALLOCATABLE :: NINLOH
!$OMP THREADPRIVATE(NINLOH)
REAL, DIMENSION(:,:), ALLOCATABLE :: XLA
!$OMP THREADPRIVATE(XLA)
REAL, DIMENSION(:), ALLOCATABLE :: XOLA, XOLO
!$OMP THREADPRIVATE(XOLA, XOLO)
INTEGER, DIMENSION(:,:), ALLOCATABLE :: NP
!$OMP THREADPRIVATE(NP)
REAL, DIMENSION(:,:), ALLOCATABLE :: XLOPH
!$OMP THREADPRIVATE(XLOPH)
!
END MODULE MODD_HORIBL
