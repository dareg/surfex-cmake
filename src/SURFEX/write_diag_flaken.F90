!     #########
SUBROUTINE WRITE_DIAG_FLAKE_n (DTCO, DGU, IOB, U, FM, &
                               HPROGRAM,HWRITE)
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
USE MODD_SURFEX_n, ONLY : FLAKE_MODEL_t
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DIAG_SURF_ATM_n, ONLY : DIAG_SURF_ATM_t
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
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
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DIAG_SURF_ATM_t), INTENT(INOUT) :: DGU
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(FLAKE_MODEL_t), INTENT(INOUT) :: FM
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
   IF (FM%DGF%XDIAG_TSTEP==XUNDEF .OR. &
           ABS(NINT(FM%F%TTIME%TIME/FM%DGF%XDIAG_TSTEP)*FM%DGF%XDIAG_TSTEP-FM%F%TTIME%TIME)<1.E-3 ) THEN
      CALL WRITE_DIAG_SEB_FLAKE_n(DTCO, DGU, IOB, U, FM%CHF, FM%DGF, &
                                  HPROGRAM)
      CALL WRITE_DIAG_MISC_FLAKE_n(DTCO, DGU, IOB, U, FM%DGMF, &
                                   HPROGRAM)
   END IF
!        
ENDIF
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_FLAKE_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_DIAG_FLAKE_n
