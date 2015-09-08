!     #############################################################
      SUBROUTINE ALLOC_TRIP(KMODEL)
!     #############################################################
!
!!    AUTHOR
!!    ------
!!      R. El Khatib   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2010
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_TRIP_DIAG, ONLY : TRIP_DIAG_ALLOC
USE MODD_TRIP_GRID, ONLY : TRIP_GRID_ALLOC
USE MODD_TRIP,      ONLY : TRIP_ALLOC
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
INTEGER, INTENT(IN) :: KMODEL    ! actual number of surfaces
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*       0.2   Declarations of local variables
!              -------------------------------

!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ALLOC_TRIP',0,ZHOOK_HANDLE)
CALL TRIP_DIAG_ALLOC(KMODEL)
CALL TRIP_GRID_ALLOC(KMODEL)
CALL TRIP_ALLOC(KMODEL)
IF (LHOOK) CALL DR_HOOK('ALLOC_TRIP',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ALLOC_TRIP
