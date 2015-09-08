!     #################################################################################
SUBROUTINE DEALLOC_TRIP
!     #################################################################################
!
!!****  *DEALLOC_TRIP * - Deallocate  the surfaces
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     R. El Khatib
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2010
!!------------------------------------------------------------------
!
USE MODD_TRIP_DIAG,   ONLY : TRIP_DIAG_DEALLO
USE MODD_TRIP_GRID,   ONLY : TRIP_GRID_DEALLO
USE MODD_TRIP,        ONLY : TRIP_DEALLO

USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('DEALLOC_TRIP',0,ZHOOK_HANDLE)
CALL TRIP_DIAG_DEALLO
CALL TRIP_GRID_DEALLO
CALL TRIP_DEALLO
IF (LHOOK) CALL DR_HOOK('DEALLOC_TRIP',1,ZHOOK_HANDLE)
!
!--------------------------------------------------------------------------------------
!
END SUBROUTINE DEALLOC_TRIP
