!     #########
SUBROUTINE WRITE_DIAG_SEAFLUX_n (DTCO, DGU, IOB, U, CHS, DGS, DGSI, S, DGO, &   
                                 HPROGRAM,HWRITE)
!     ###############################################################################
!
!!****  *WRITE_DIAG_SEAFLUX_n * - diagnostics for SEAFLUX
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
!!      Modified    09/2013 : S. Senesi : call WRITE_DIAG_SEB_SEAICE_n
!!------------------------------------------------------------------
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DIAG_SURF_ATM_n, ONLY : DIAG_SURF_ATM_t
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_CH_SEAFLUX_n, ONLY : CH_SEAFLUX_t
USE MODD_DIAG_SEAFLUX_n, ONLY : DIAG_SEAFLUX_t
USE MODD_DIAG_SEAICE_n, ONLY : DIAG_SEAICE_t
USE MODD_SEAFLUX_n, ONLY : SEAFLUX_t
USE MODD_DIAG_OCEAN_n, ONLY : DIAG_OCEAN_t
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
!
USE MODI_WRITE_DIAG_SEB_SEAFLUX_n
USE MODI_WRITE_DIAG_SEB_OCEAN_n
USE MODI_WRITE_DIAG_SEB_SEAICE_n
! 
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DIAG_SURF_ATM_t), INTENT(INOUT) :: DGU
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(CH_SEAFLUX_t), INTENT(INOUT) :: CHS
TYPE(DIAG_SEAFLUX_t), INTENT(INOUT) :: DGS
TYPE(DIAG_SEAICE_t), INTENT(INOUT) :: DGSI
TYPE(SEAFLUX_t), INTENT(INOUT) :: S
TYPE(DIAG_OCEAN_t), INTENT(INOUT) :: DGO
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
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_SEAFLUX_N',0,ZHOOK_HANDLE)
IF (HWRITE/='PGD') THEN
!        
   IF (DGS%XDIAG_TSTEP==XUNDEF .OR. ABS(NINT(S%TTIME%TIME/DGS%XDIAG_TSTEP)*DGS%XDIAG_TSTEP-S%TTIME%TIME)<1.E-3 ) THEN
      CALL WRITE_DIAG_SEB_SEAFLUX_n(DTCO, DGU, IOB, U, CHS, DGS, DGSI, S, &
                                    HPROGRAM)
      IF (DGO%LDIAG_OCEAN)  CALL WRITE_DIAG_SEB_OCEAN_n(DTCO, DGU, IOB, U, DGO, &
                                                        HPROGRAM)
      IF (DGSI%LDIAG_SEAICE) CALL WRITE_DIAG_SEB_SEAICE_n(DTCO, DGU, IOB, U, DGS, DGSI, S, &
                                                          HPROGRAM)
   END IF
!        
ENDIF
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_SEAFLUX_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_DIAG_SEAFLUX_n
