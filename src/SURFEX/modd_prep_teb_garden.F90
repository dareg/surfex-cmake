!     ################
      MODULE MODD_PREP_TEB_GARDEN
!     ################
!
!!****  *MODD_PREP - declaration for field interpolations
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
!!	V. Masson    *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2004
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
SAVE
!--------------------------------------------------------------------------
!
 CHARACTER(LEN=28) :: CFILE_GD     ! input file name
 CHARACTER(LEN=6)  :: CTYPE          ! input file type
 CHARACTER(LEN=28) :: CFILEPGD_GD  ! input file name
 CHARACTER(LEN=6)  :: CTYPEPGD       ! input file type
 CHARACTER(LEN=28) :: CFILE_SNOW_GD     ! input file name for Snow
 CHARACTER(LEN=6)  :: CTYPE_SNOW     ! input file type for Snow
 CHARACTER(LEN=28) :: CFILEPGD_SNOW_GD     ! input file name for Snow
 CHARACTER(LEN=6)  :: CTYPEPGD_SNOW     ! input file type for Snow 
 CHARACTER(LEN=28) :: CFILE_HUG_GD      ! input file name for Wg, Wgi
 CHARACTER(LEN=6)  :: CTYPE_HUG      ! input file type for Wg, Wgi
 CHARACTER(LEN=28) :: CFILE_TG_GD       ! input file name for Tg
 CHARACTER(LEN=6)  :: CTYPE_TG       ! input file type for Tg
 CHARACTER(LEN=28) :: CFILE_HUG_SURF_GD ! input file name for HUG_SURF
 CHARACTER(LEN=28) :: CFILE_HUG_ROOT_GD ! input file name for HUG_ROOT
 CHARACTER(LEN=28) :: CFILE_HUG_DEEP_GD ! input file name for HUG_DEEP
 CHARACTER(LEN=28) :: CFILE_TG_SURF_GD  ! input file name for TG_SURF
 CHARACTER(LEN=28) :: CFILE_TG_ROOT_GD  ! input file name for TG_ROOT
 CHARACTER(LEN=28) :: CFILE_TG_DEEP_GD  ! input file name for TG_DEEP
!
REAL              :: XHUG_SURF_GD      ! surface relative soil humidity
REAL              :: XHUG_ROOT_GD      ! root layer relative soil humidity
REAL              :: XHUG_DEEP_GD      ! deep layer relative soil humidity
REAL              :: XHUGI_SURF_GD     ! surf layer relative ice content
REAL              :: XHUGI_ROOT_GD     ! root layer relative ice content
REAL              :: XHUGI_DEEP_GD     ! deep layer relative ice content
REAL              :: XTG_SURF_GD       ! surface temperature
REAL              :: XTG_ROOT_GD       ! root layer temperature
REAL              :: XTG_DEEP_GD       ! deep layer temperature
!
LOGICAL :: LSNOW_IDEAL_GD 
!
REAL, DIMENSION(:), POINTER :: XWSNOW_GD         ! Snow reservoir
REAL, DIMENSION(:), POINTER :: XRSNOW_GD         ! snow density
REAL, DIMENSION(:), POINTER :: XTSNOW_GD         ! snow temperature
REAL              :: XASNOW_GD         ! snow albedo
!
REAL              :: XWR_DEF        ! default for leaves interception reservoir
!--------------------------------------------------------------------------
!
!* normalized dimensions for interpolation grids for soil
INTEGER, PARAMETER           :: NGRID_LEVEL = 20
REAL, DIMENSION(NGRID_LEVEL) :: XGRID_SOIL = (/ 0., 0.01, 0.02, 0.05, 0.1, 0.2,&
    0.4, 0.7, 1., 1.3, 1.6, 2., 2.5, 3., 3.5, 4., 5., 6., 8., 10.  /)  
!
!--------------------------------------------------------------------------
!
END MODULE MODD_PREP_TEB_GARDEN


