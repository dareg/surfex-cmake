!     ##################
      MODULE MODD_WRITE_TXT
!     ##################
!
!!****  *MODD_WRITE_TXT - 
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
!!	P. LE MOIGNE    *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!    
!
!*       0.   DECLARATIONS
!
IMPLICIT NONE
!      
INTEGER, PARAMETER                  :: JPVAR = 600             ! maximum number of fields to write      
 CHARACTER(LEN=12), DIMENSION(JPVAR) :: CVAR='                ' ! names of fields to write
 CHARACTER(LEN=12), DIMENSION(JPVAR) :: CVARN='                ' ! names of fields to write
INTEGER, DIMENSION(JPVAR)           :: NVAR                    ! unit number associated to CVAR elements
INTEGER                             :: NIND                    ! current unit number
INTEGER                             :: NUNIT0=33
!
END MODULE MODD_WRITE_TXT
