!     #########
SUBROUTINE DEALLOC_TOWN_n
!     ###############################################################################
!
!!****  *DEALLOC_TOWN_n * - Deallocate all arrays
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
USE MODD_BEM_n, ONLY : B => BEM
USE MODD_BEM_OPTION_n, ONLY : BOP => BEM_OPTIONS
USE MODD_CH_TEB_n, ONLY : CHT => CH_TEB
USE MODD_DATA_TEB_n, ONLY : DTT => DATA_TEB
USE MODD_TEB_GRID_n, ONLY : TG => TEB_GRID
USE MODD_TEB_n, ONLY : T => TEB
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS
USE MODD_TEB_PANEL_n, ONLY : TPN => TEB_PANEL
!
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_DEALLOC_IDEAL_FLUX
!
USE MODI_DEALLOC_TEB_n
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

IF (LHOOK) CALL DR_HOOK('DEALLOC_TOWN_N',0,ZHOOK_HANDLE)
IF (U%CTOWN=='TEB   ') THEN
  CALL DEALLOC_TEB_n(B, BOP, CHT, DTT, TG, T, TOP, TPN)
ELSE IF (U%CTOWN=='FLUX  ') THEN
  CALL DEALLOC_IDEAL_FLUX
END IF
IF (LHOOK) CALL DR_HOOK('DEALLOC_TOWN_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DEALLOC_TOWN_n
