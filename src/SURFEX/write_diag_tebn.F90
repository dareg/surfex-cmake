!     #########
SUBROUTINE WRITE_DIAG_TEB_n(HPROGRAM,HWRITE)
!     ###############################################################################
!
!!****  *WRITE_DIAG_TEB_n * - diagnostics for TEB
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
USE MODD_CH_TEB_n, ONLY : CHT => CH_TEB
USE MODD_DIAG_UTCI_TEB_n, ONLY : DGUT => DIAG_UTCI_TEB
!
USE MODD_BEM_n, ONLY : B => BEM
USE MODD_BEM_OPTION_n, ONLY : BOP => BEM_OPTIONS
USE MODD_TEB_PANEL_n, ONLY : TPN => TEB_PANEL
!
USE MODD_TEB_GARDEN_PGD_EVOL_n, ONLY : TGDPE => TEB_GARDEN_PGD_EVOL
USE MODD_TEB_GARDEN_PGD_n, ONLY : TGDP => TEB_GARDEN_PGD
USE MODD_TEB_VEG_n, ONLY : TVG => TEB_VEG_OPTIONS
!
USE MODD_DIAG_CUMUL_TEB_n, ONLY : DGCT => DIAG_CUMUL_TEB
USE MODD_DIAG_MISC_TEB_n, ONLY : DGMT => DIAG_MISC_TEB
USE MODD_DIAG_MISC_TEB_OPTION_n, ONLY : DGMTO => DIAG_MISC_TEB_OPTIONS
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS
USE MODD_TEB_n, ONLY : T => TEB
USE MODD_DIAG_TEB_n, ONLY : DGT => DIAG_TEB
!
USE MODI_GOTO_TEB
USE MODI_WRITE_DIAG_SEB_TEB_n
USE MODI_WRITE_DIAG_MISC_TEB_n
USE MODI_WRITE_DIAG_PGD_TEB_n
USE MODI_WRITE_DIAG_PGD_GRDN_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM ! program calling surf. schemes
 CHARACTER(LEN=3),   INTENT(IN)  :: HWRITE   ! 'PGD' : only physiographic fields are written
!                                           ! 'ALL' : all fields are written
!
!*      0.2    declarations of local variables
!
INTEGER         :: JTEB_PATCH
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_TEB_N',0,ZHOOK_HANDLE)
IF (HWRITE/='PGD') THEN
!        
   IF (DGT%XDIAG_TSTEP==XUNDEF .OR. ABS(NINT(TOP%TTIME%TIME/DGT%XDIAG_TSTEP)*DGT%XDIAG_TSTEP-TOP%TTIME%TIME)<1.E-3 ) THEN
      CALL WRITE_DIAG_SEB_TEB_n(CHT, DGT, DGUT, &
                                HPROGRAM)
      DO JTEB_PATCH=1,TOP%NTEB_PATCH
        CALL GOTO_TEB(JTEB_PATCH)
        CALL WRITE_DIAG_MISC_TEB_n(DGCT, DGMT, DGMTO, T, TOP, &
                                   HPROGRAM,JTEB_PATCH)
      END DO      
   END IF
!
ENDIF
!
IF (DGT%LPGD) THEN
  IF (DGT%XDIAG_TSTEP==XUNDEF .OR. ABS(NINT(TOP%TTIME%TIME/DGT%XDIAG_TSTEP)*DGT%XDIAG_TSTEP-TOP%TTIME%TIME)<1.E-3 ) THEN
    IF (ASSOCIATED(T%XBLD)) THEN
      CALL WRITE_DIAG_PGD_TEB_n(B, BOP, T, TOP, TPN, &
                                HPROGRAM)
      IF (TOP%LGARDEN) CALL WRITE_DIAG_PGD_GRDN_n(DGMTO, TGDPE, TGDP, TVG, &
                                                  HPROGRAM)
    ENDIF
  END IF
END IF
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_TEB_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_DIAG_TEB_n
