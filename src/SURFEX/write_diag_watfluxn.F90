!     #########
SUBROUTINE WRITE_DIAG_WATFLUX_n(HPROGRAM,HWRITE)
!     ###############################################################################
!
!!****  *WRITE_DIAG_WATFLUX_n * - diagnostics for lakes
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
USE MODD_CH_WATFLUX_n, ONLY : CHW => CH_WATFLUX
USE MODD_DIAG_SURF_ATM_n, ONLY : DGU => DIAG_SURF_ATM
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_WATFLUX_n, ONLY : W => WATFLUX
USE MODD_DIAG_WATFLUX_n, ONLY : DGW => DIAG_WATFLUX
!
USE MODI_WRITE_DIAG_SEB_WATFLUX_n
! 
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
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!                                           ! 'ALL' : all fields are written
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_WATFLUX_N',0,ZHOOK_HANDLE)
IF (HWRITE/='PGD') THEN
!        
   IF (DGW%XDIAG_TSTEP==XUNDEF .OR. ABS(NINT(W%TTIME%TIME/DGW%XDIAG_TSTEP)*DGW%XDIAG_TSTEP-W%TTIME%TIME)<1.E-3 ) THEN
      CALL WRITE_DIAG_SEB_WATFLUX_n(CHW, DGU, DGW, &
                                    HPROGRAM)
   END IF
!        
ENDIF
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_WATFLUX_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_DIAG_WATFLUX_n
