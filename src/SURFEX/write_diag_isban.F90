!     #########
SUBROUTINE WRITE_DIAG_ISBA_n (DTCO, DGU, IOB, U, CHI, DGEI, DGI, DGMI, DST, GB, I, &
                              HPROGRAM,HWRITE)
!     ###############################################################################
!
!!****  *WRITE_DIAG_ISBA_n * - Stores ISBA diagnostics
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
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DIAG_SURF_ATM_n, ONLY : DIAG_SURF_ATM_t
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_CH_ISBA_n, ONLY : CH_ISBA_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_ISBA_n, ONLY : DIAG_ISBA_t
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t
USE MODD_DST_n, ONLY : DST_t
USE MODD_GR_BIOG_n, ONLY : GR_BIOG_t
USE MODD_ISBA_n, ONLY : ISBA_t
!
USE MODD_SURF_PAR,    ONLY : XUNDEF
! 
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_WRITE_DIAG_MISC_ISBA_n
USE MODI_WRITE_DIAG_PGD_ISBA_n
USE MODI_WRITE_DIAG_SEB_ISBA_n
!
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
TYPE(CH_ISBA_t), INTENT(INOUT) :: CHI
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEI
TYPE(DIAG_ISBA_t), INTENT(INOUT) :: DGI
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DGMI
TYPE(DST_t), INTENT(INOUT) :: DST
TYPE(GR_BIOG_t), INTENT(INOUT) :: GB
TYPE(ISBA_t), INTENT(INOUT) :: I
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM ! program calling surf. schemes
 CHARACTER(LEN=3),   INTENT(IN)  :: HWRITE    ! 'PGD' : only physiographic fields are written
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!                                            ! 'ALL' : all fields are written
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_ISBA_N',0,ZHOOK_HANDLE)
IF (HWRITE/='PGD') THEN
  IF (DGI%XDIAG_TSTEP==XUNDEF .OR. ABS(NINT(I%TTIME%TIME/DGI%XDIAG_TSTEP)*DGI%XDIAG_TSTEP-I%TTIME%TIME)<1.E-3 ) THEN
    CALL WRITE_DIAG_SEB_ISBA_n(DTCO, DGU, IOB, U, CHI, DGEI, DGI, DST, GB, I, &
                               HPROGRAM)
    CALL WRITE_DIAG_MISC_ISBA_n(DTCO, DGU, IOB, U, DGI, DGMI, I, &
                                HPROGRAM)
  END IF
END IF
!
IF (DGI%LPGD) THEN
  IF (DGI%XDIAG_TSTEP==XUNDEF .OR. ABS(NINT(I%TTIME%TIME/DGI%XDIAG_TSTEP)*DGI%XDIAG_TSTEP-I%TTIME%TIME)<1.E-3 ) THEN
    CALL WRITE_DIAG_PGD_ISBA_n(DTCO, DGU, IOB, U, CHI, DGMI, I, &
                               HPROGRAM)
  END IF
END IF
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_ISBA_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_DIAG_ISBA_n
