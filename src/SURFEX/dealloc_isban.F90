!     #################################################################################
SUBROUTINE DEALLOC_ISBA_n (IM)
!     #################################################################################
!
!!****  *DEALLOC_ISBA_n * - Deallocate all arrays
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
!!      P Samuelsson 10/2014  MEB
!!------------------------------------------------------------------
!
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_INIT
USE MODD_ISBA_n, ONLY : ISBA_S_INIT, ISBA_K_INIT, ISBA_NK_INIT, ISBA_NP_INIT, ISBA_NPE_INIT
USE MODD_SSO_n, ONLY : SSO_INIT, SSO_NP_INIT
USE MODD_GRID_n, ONLY : GRID_INIT, GRID_NP_INIT
USE MODD_CH_ISBA_n, ONLY : CH_ISBA_INIT, CH_ISBA_NP_INIT
USE MODD_AGRI_n, ONLY : AGRI_NP_INIT
USE MODD_GR_BIOG_n, ONLY : GR_BIOG_NP_INIT
!
USE MODD_SURFEX_n, ONLY : ISBA_MODEL_t
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
TYPE(ISBA_MODEL_t), INTENT(INOUT) :: IM
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('DEALLOC_ISBA_N',0,ZHOOK_HANDLE)
!
CALL ISBA_S_INIT(IM%S)
CALL ISBA_K_INIT(IM%K)
CALL SSO_INIT(IM%ISS)
!
CALL ISBA_NK_INIT(IM%NK,IM%O%NPATCH)
CALL SSO_NP_INIT(IM%NISS,IM%O%NPATCH)
CALL ISBA_NPE_INIT(IM%NPE,IM%O%NPATCH)
!
CALL ISBA_NP_INIT(IM%NP,IM%O%NPATCH)
!
CALL GRID_INIT(IM%G)
CALL GRID_NP_INIT(IM%NG,IM%O%NPATCH)
!
CALL CH_ISBA_INIT(IM%CHI)
CALL CH_ISBA_NP_INIT(IM%NCHI,IM%O%NPATCH)
!
CALL DATA_ISBA_INIT(IM%DTI)
!
CALL AGRI_NP_INIT(IM%NAG, IM%O%NPATCH)
!
CALL GR_BIOG_NP_INIT(IM%NGB, IM%O%NPATCH)
!
!-------------------------------------------------------------------------------------
!
IF (ASSOCIATED(IM%GB%XMONOPOT)) DEALLOCATE(IM%GB%XMONOPOT)
IF (ASSOCIATED(IM%GB%XISOPOT))  DEALLOCATE(IM%GB%XISOPOT)
!
!-------------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DEALLOC_ISBA_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DEALLOC_ISBA_n
