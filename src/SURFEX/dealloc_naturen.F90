!     #########
SUBROUTINE DEALLOC_NATURE_n
!     ###############################################################################
!
!!****  *DEALLOC_NATURE_n * - Deallocate all arrays
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

!
USE MODD_CH_ISBA_n, ONLY : CHI => CH_ISBA
USE MODD_DATA_ISBA_n, ONLY : DTI => DATA_ISBA
USE MODD_GR_BIOG_n, ONLY : GB => GR_BIOG
USE MODD_ISBA_GRID_n, ONLY : IG => ISBA_GRID
USE MODD_ISBA_n, ONLY : I => ISBA
!
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_DEALLOC_IDEAL_FLUX
!
USE MODI_DEALLOC_ISBA_n
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

IF (LHOOK) CALL DR_HOOK('DEALLOC_NATURE_N',0,ZHOOK_HANDLE)
IF (U%CNATURE=='ISBA  ' .OR. U%CNATURE=='TSZ0  ') THEN
  CALL DEALLOC_ISBA_n(CHI, DTI, GB, IG, I)
ELSE IF (U%CNATURE=='FLUX  ') THEN
  CALL DEALLOC_IDEAL_FLUX
END IF
IF (LHOOK) CALL DR_HOOK('DEALLOC_NATURE_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DEALLOC_NATURE_n
