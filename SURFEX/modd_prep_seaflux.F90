!     ################
      MODULE MODD_PREP_SEAFLUX
!     ################
!
!!****  *MODD_PREP_SEAFLUX - declaration for field interpolations
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
!!	S.Malardel    *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       10/03
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
!--------------------------------------------------------------------------
!
CHARACTER(LEN=28) :: CFILE_SEAFLX   ! input file name
CHARACTER(LEN=6)  :: CTYPE          ! input file type
!
REAL              :: XSST_UNIF   !  uniform prescribed SST
!
!--------------------------------------------------------------------------
!
END MODULE MODD_PREP_SEAFLUX


