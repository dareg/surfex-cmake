!     ##################
      MODULE MODN_PGD_GRID
!     ##################
!
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!       
!!    AUTHOR
!!    ------
!!	V. Masson    *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2003                    
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PGD_GRID, ONLY : CGRID
!
USE MODD_POINT_OVERLAY, ONLY : NOVMX
!
IMPLICIT NONE
!
 CHARACTER(LEN=28):: YINIFILE ! name of input file
 CHARACTER(LEN=6) :: YINIFILETYPE! type of input file
!
!
NAMELIST/NAM_PGD_GRID/CGRID,NOVMX,YINIFILE,YINIFILETYPE
!
END MODULE MODN_PGD_GRID
