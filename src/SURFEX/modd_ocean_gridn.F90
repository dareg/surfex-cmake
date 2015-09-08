!     ##################
      MODULE MODD_OCEAN_GRID_n
!     ##################
!
!!****  *MODD_OCEAN_GRID_n - declaration of grid for oceanic model
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
!!      C. Lebeaupin Brossier   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2008
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!------------------------------------------------------------------------------
!
!TYPE OCEAN_GRID_t
!!-------------------------------------------------------------------------------
!!
!! Grid definition
!!
!  REAL,POINTER,     DIMENSION(:)  :: XK1,XK2,XK3,XK4 !oceanic lengths between levels
!  REAL,POINTER,     DIMENSION(:)  :: XZHOC           !oceanic levels depth
!  REAL,POINTER,     DIMENSION(:)  :: XZ2,XDZ1,XDZ2   !oceanic levels and inter-levels
!  REAL,POINTER,     DIMENSION(:)  :: XRAY            !solar penetration coef. 
!!
!!-------------------------------------------------------------------------------
!!
!
!END TYPE OCEAN_GRID_t



REAL, POINTER, DIMENSION(:) :: XK1
REAL, POINTER, DIMENSION(:) :: XK2
REAL, POINTER, DIMENSION(:) :: XK3
REAL, POINTER, DIMENSION(:) :: XK4
REAL, POINTER, DIMENSION(:) :: XZHOC
REAL, POINTER, DIMENSION(:) :: XZ2
REAL, POINTER, DIMENSION(:) :: XDZ1
REAL, POINTER, DIMENSION(:) :: XDZ2
REAL, POINTER, DIMENSION(:) :: XRAY

CONTAINS

!!
!! Save current state for allocated arrays
!!
!
!SUBROUTINE OCEAN_GRID_INIT(Y!SUBROUTINE OCEAN_GRID)
!TYPE(SUBROUTINE OCEAN_GRID_t), INTENT(INOUT) :: YSUBROUTINE OCEAN_GRID
!REAL(KIND=JPRB) :: ZHOOK_HANDLE
!IF (LHOOK) CALL DR_HOOK("MODD_OCEAN_GRID_N:OCEAN_GRID_INIT",0,ZHOOK_HANDLE)
!  NULLIFY(YOCEAN_GRID%XK1)
!  NULLIFY(YOCEAN_GRID%XK2)
!  NULLIFY(YOCEAN_GRID%XK3)
!  NULLIFY(YOCEAN_GRID%XK4)
!  NULLIFY(YOCEAN_GRID%XZHOC)
!  NULLIFY(YOCEAN_GRID%XZ2)
!  NULLIFY(YOCEAN_GRID%XDZ1)
!  NULLIFY(YOCEAN_GRID%XDZ2)
!  NULLIFY(YOCEAN_GRID%XRAY)
!IF (LHOOK) CALL DR_HOOK("MODD_OCEAN_GRID_N:OCEAN_GRID_INIT",1,ZHOOK_HANDLE)
!END SUBROUTINE OCEAN_GRID_INIT
!

END MODULE MODD_OCEAN_GRID_n
