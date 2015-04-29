!     #########
SUBROUTINE WRITE_DIAG_FLAKE_n(HPROGRAM,HWRITE)
!     ###############################################################################
!
!!****  *WRITE_DIAG_FLAKE_n * - diagnostics for lakes
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
USE MODD_CH_FLAKE_n, ONLY : CHF => CH_FLAKE
USE MODD_DIAG_SURF_ATM_n, ONLY : DGU => DIAG_SURF_ATM
!
USE MODD_DIAG_MISC_FLAKE_n, ONLY : DGMF => DIAG_MISC_FLAKE
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_FLAKE_n, ONLY : F => FLAKE
USE MODD_DIAG_FLAKE_n, ONLY : DGF => DIAG_FLAKE
!
USE MODI_WRITE_DIAG_SEB_FLAKE_n
! 
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_WRITE_DIAG_MISC_FLAKE_n
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM ! program calling surf. schemes
 CHARACTER(LEN=3),   INTENT(IN)  :: HWRITE   ! 'PGD' : only physiographic fields are written
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!                                           ! 'ALL' : all fields are written
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_FLAKE_N',0,ZHOOK_HANDLE)
IF (HWRITE/='PGD') THEN
!        
   IF (DGF%XDIAG_TSTEP==XUNDEF .OR. ABS(NINT(F%TTIME%TIME/DGF%XDIAG_TSTEP)*DGF%XDIAG_TSTEP-F%TTIME%TIME)<1.E-3 ) THEN
      CALL WRITE_DIAG_SEB_FLAKE_n(CHF, DGF, DGU, &
                                  HPROGRAM)
      CALL WRITE_DIAG_MISC_FLAKE_n(DGMF, &
                                   HPROGRAM)
   END IF
!        
ENDIF
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_FLAKE_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_DIAG_FLAKE_n
