!     ##################
      MODULE MODD_WRITE_BIN
!     ##################
!
!!****  *MODD_WRITE_BIN - 
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
!!      A. LEMONSU      *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!    
!
!*       0.   DECLARATIONS
!
IMPLICIT NONE
!      
INTEGER, PARAMETER                  :: JPVAR = 700             ! maximum number of fields to write      
 CHARACTER(LEN=12), DIMENSION(JPVAR) :: CVAR='                ' ! names of fields to write
INTEGER, DIMENSION(JPVAR)           :: NVAR                    ! unit number associated to CVAR elements
INTEGER                             :: NIND                    ! current unit number
INTEGER                             :: NUNIT0=33
INTEGER                             :: NWRITE                  ! counter for writing
!
!
END MODULE MODD_WRITE_BIN
