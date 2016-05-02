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
!!      F. Habets   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      modified 04/04 by P. LeMoigne: add logical for town, sea and water
!!      modified 07/11 by B. Decharme: add mask for IGN grid
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
!$OMP THREADPRIVATE(NMASK)
INTEGER                       :: NSTEP_OUTPUT
LOGICAL                       :: LMASK = .FALSE.
LOGICAL                       :: LPARTR,LPARTW
LOGICAL,DIMENSION(5)          :: LTIME_WRITTEN
INTEGER                       :: XTYPE
INTEGER                       :: XSTART,XCOUNT,XSTRIDE
INTEGER                       :: XSTARTW,XCOUNTW

INTEGER, DIMENSION(:),ALLOCATABLE :: NMASK_IGN

!------------------------------------------------------------------------------
!
END MODULE MODD_IO_SURF_OL

