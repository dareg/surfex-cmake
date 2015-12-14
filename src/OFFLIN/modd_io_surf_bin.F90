!     ##################
      MODULE MODD_IO_SURF_BIN
!     ##################
!
!!****  *MODD_IO_SURF_BIN - 
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
!!      A. LEMONSU   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!    
!
!*       0.   DECLARATIONS
!
IMPLICIT NONE
INTEGER, DIMENSION(:), POINTER :: NMASK                   ! 1D mask to read only interesting
 CHARACTER(LEN=6)               :: CMASK                   ! surface mask type
 CHARACTER(LEN=28),SAVE         :: CFILEIN  ='SURFIN.txt'  ! Name of the input
 CHARACTER(LEN=28),SAVE         :: CFILEOUT ='SURFOUT.txt' ! Name of the input
INTEGER                        :: NFULL                   ! total number for points of surface
INTEGER                        :: NUNIT                   ! logical unit of surface file
INTEGER                        :: NLUOUT                  ! logical unit of output file
!
END MODULE MODD_IO_SURF_BIN
