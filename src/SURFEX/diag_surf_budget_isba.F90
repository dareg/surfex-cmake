!     #########
       SUBROUTINE DIAG_SURF_BUDGET_ISBA (PDIR_SW, PSCA_SW, PLW, INI, DGIP )  
!     ###############################################################################
!
!!****  *DIAG_SURF_BUDGET_ISBA * - Computes diagnostics over ISBA
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
!!     P. Le Moigne 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2006
!!      Modified    08/2008 (B. Decharme) LWU diag
!!------------------------------------------------------------------
!
USE MODD_DIAG_n, ONLY : DIAG_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_t
!
USE MODD_CSTS,           ONLY : XSTEFAN
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:),INTENT(IN)  :: PDIR_SW   ! direct  solar radiation (on horizontal surf.)
!                                             !                                       (W/m2)
REAL, DIMENSION(:,:),INTENT(IN)  :: PSCA_SW   ! diffuse solar radiation (on horizontal surf.)
!                                             !                                       (W/m2)
REAL, DIMENSION(:), INTENT(IN)   :: PLW       ! longwave radiation (on horizontal surf.)
!
TYPE(ISBA_INIT_t), INTENT(INOUT) :: INI
TYPE(DIAG_t), INTENT(INOUT) :: DGIP
!
!*      0.2    declarations of local variables
!
INTEGER                          :: ISWB      ! number of SW bands
INTEGER                          :: JSWB      ! loop counter on number of SW bands
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_BUDGET_ISBA',0,ZHOOK_HANDLE)
ISWB = SIZE(PDIR_SW,2)
! 
!* total incoming and outgoing SW
!
DO JSWB=1,ISWB
  DGIP%XSWBD(:,JSWB) = PDIR_SW(:,JSWB) + PSCA_SW(:,JSWB)
  DGIP%XSWBU(:,JSWB) = PDIR_SW(:,JSWB) * INI%XDIR_ALB_WITH_SNOW(:,JSWB,1) + &
                       PSCA_SW(:,JSWB) * INI%XSCA_ALB_WITH_SNOW(:,JSWB,1) 
ENDDO
!
DGIP%XSWD(:) = 0.
DGIP%XSWU(:) = 0.
DO JSWB=1,ISWB
   DGIP%XSWD(:) = DGIP%XSWD(:) + DGIP%XSWBD(:,JSWB)
   DGIP%XSWU(:) = DGIP%XSWU(:) + DGIP%XSWBU(:,JSWB)
ENDDO
!
!*incoming outgoing LW
!
!Wrong old diag : LWU=EMIS*STEFAN*Ts**4 + (1.-EMIS)*LW
!Due to e_budget.f90 linearization, LWU can not be calculated using actual Ts
!
DGIP%XLWD(:) = PLW(:)
DGIP%XLWU(:) = DGIP%XSWD(:) - DGIP%XSWU(:) + DGIP%XLWD(:) - DGIP%XRN(:)
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_BUDGET_ISBA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_SURF_BUDGET_ISBA
