!     ##################
      MODULE MODD_IO_SURF_TXT
!     ##################
!
!!****  *MODD_IO_SURF_TXT - 
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
!!	V. Masson    *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!    
!!     P. LeMoigne 04/2004 : distinguish in and out file name
!
!*       0.   DECLARATIONS
!
IMPLICIT NONE
INTEGER, DIMENSION(:), POINTER :: NMASK ! 1D mask to read only interesting
!$OMP THREADPRIVATE(NMASK)
CHARACTER(LEN=6)               :: CMASK ! surface mask type
!$OMP THREADPRIVATE(CMASK)
INTEGER                        :: NFULL ! total number for points of surface
!$OMP THREADPRIVATE(NFULL)
!
END MODULE MODD_IO_SURF_TXT
