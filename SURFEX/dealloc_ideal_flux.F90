!     #################################################################################
SUBROUTINE DEALLOC_IDEAL_FLUX
!     #################################################################################
!
!!****  *DEALLOC_IDEAL_FLUX * - Deallocate all arrays
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
!!     V. Masson 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!------------------------------------------------------------------
!
USE MODD_IDEAL_FLUX, ONLY : XSFTS
!
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

IF (LHOOK) CALL DR_HOOK('DEALLOC_IDEAL_FLUX',0,ZHOOK_HANDLE)
IF (ALLOCATED(XSFTS))  DEALLOCATE(XSFTS)
IF (LHOOK) CALL DR_HOOK('DEALLOC_IDEAL_FLUX',1,ZHOOK_HANDLE)
!
!--------------------------------------------------------------------------------------
!
END SUBROUTINE DEALLOC_IDEAL_FLUX


