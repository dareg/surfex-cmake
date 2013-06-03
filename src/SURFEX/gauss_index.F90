!     #########
      SUBROUTINE GAUSS_INDEX(HPROGRAM,HFILE,HFLAG,OINDEX_STORE,HINDEX_1KM,&
                               HINDEX_10KM,HINDEX_100KM,HRES_COMP)  
!     #####################################################################
!!
!!    PURPOSE
!!    -------
!!     
!!    Read or calculate and store gauss indexes
!!
!!    METHOD
!!    ------
!!   
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    B. Decharme                   Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original     02/2010
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_PGD_GRID,   ONLY : XMESHLENGTH
USE MODD_PGDWORK,    ONLY : NSSO
!
USE MODI_GET_LUOUT
USE MODI_INI_SSOWORK
!
USE MODE_GAUSS_INDEX
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of dummy arguments
!            ------------------------------
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM
 CHARACTER(LEN=28), INTENT(IN)  :: HFILE
 CHARACTER(LEN=6),  INTENT(IN)  :: HFLAG
LOGICAL,           INTENT(IN)  :: OINDEX_STORE      ! Store index in a binary file
 CHARACTER(LEN=28), INTENT(IN)  :: HINDEX_1KM
 CHARACTER(LEN=28), INTENT(IN)  :: HINDEX_10KM
 CHARACTER(LEN=28), INTENT(IN)  :: HINDEX_100KM
 CHARACTER(LEN=5), INTENT(INOUT):: HRES_COMP
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
 CHARACTER(LEN=5)            :: YRES         ! Resolution in file 
REAL, DIMENSION(:), POINTER :: ZLAT         ! latitude of data points
REAL, DIMENSION(:), POINTER :: ZLON         ! longitude of data points
REAL                        :: ZDLAT        ! latitude mesh in the data file
REAL                        :: ZDLON        ! longitude mesh in the data file
INTEGER                     :: INDIM        ! number of grid point in file
INTEGER                     :: INLON        ! number of longitude rows in file
INTEGER                     :: INLAT        ! number of latitude  rows in file
INTEGER                     :: ILUOUT       ! output listing logical unit
 CHARACTER(LEN=28)           :: YINDEX       ! file name for gauss index           
!
INTEGER, DIMENSION(:), ALLOCATABLE :: IINDEX ! mesh index of all input points
                                             ! 0 indicates the point is out of the domain
INTEGER, DIMENSION(:), ALLOCATABLE :: ISSOX  ! X submesh index in their mesh of all input points
INTEGER, DIMENSION(:), ALLOCATABLE :: ISSOY  ! Y submesh index in their mesh of all input points
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!  
!-------------------------------------------------------------------------------
!
END SUBROUTINE GAUSS_INDEX
