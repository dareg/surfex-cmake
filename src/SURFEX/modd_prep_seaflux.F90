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
!!      Modified     09/2013 : S. Senesi : introduce variables for sea-ice model
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
SAVE
!--------------------------------------------------------------------------
!
 CHARACTER(LEN=28) :: CFILE_SEAFLX   ! input file name
 CHARACTER(LEN=6)  :: CTYPE_SEAFLX   ! input file type
 CHARACTER(LEN=28) :: CFILEPGD_SEAFLX   ! input file name
 CHARACTER(LEN=6)  :: CTYPEPGD          ! input file type
!
REAL              :: XSST_UNIF   !  uniform prescribed SST
REAL              :: XSSS_UNIF   !  uniform prescribed SSS
REAL              :: XSIC_UNIF   !  uniform prescribed Seaice cover
!
!--------------------------------------------------------------------------
!
END MODULE MODD_PREP_SEAFLUX


