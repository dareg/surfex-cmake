!     ######################
      MODULE MODD_IO_SURF_OL
!     ######################
!
!!****  *MODD_IO_SURF_OL* Keep in memory the netcdf ID of the output files
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
!!	F. Habets   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      modified 04/04 by P. LeMoigne: add logical for town, sea and water
!
!*       0.   DECLARATIONS
!             ------------
!
!
IMPLICIT NONE
!------------------------------------------------------------------------------
!
!* variables for each patch
!
INTEGER, DIMENSION(:),POINTER :: NMASK
INTEGER                       :: NSTEP_OUTPUT
LOGICAL                       :: LMASK = .FALSE.
LOGICAL                       :: LPARTR,LPARTW
LOGICAL                       :: LDEFINED_NATURE
LOGICAL                       :: LDEFINED_SEA
LOGICAL                       :: LDEFINED_WATER
LOGICAL                       :: LDEFINED_TOWN
LOGICAL                       :: LDEFINED_SURF_ATM
LOGICAL,DIMENSION(5)          :: LTIME_WRITTEN
INTEGER                       :: XTYPE
INTEGER                       :: XSTART,XCOUNT,XSTRIDE
INTEGER                       :: XSTARTW,XCOUNTW



!------------------------------------------------------------------------------
!
END MODULE MODD_IO_SURF_OL

