!     ##################
      MODULE MODD_IO_BUFF
!     ##################
!
!!****  *MODD_IO_IO_BUFF - 
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
!!      V. Masson    *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!    
!
!*       0.   DECLARATIONS
!
!
IMPLICIT NONE

 CHARACTER(LEN=12), DIMENSION(3000) :: CREC   ! list of records already read/written
INTEGER                            :: NREC   ! number of records read/written

!
!
END MODULE MODD_IO_BUFF
