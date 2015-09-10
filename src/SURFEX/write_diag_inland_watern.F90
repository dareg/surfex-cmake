!     #########
SUBROUTINE WRITE_DIAG_INLAND_WATER_n (DTCO, DGU, IOB, U, CHW, DGW, W, CHF, DGF, DGMF, F, &
                                      HPROGRAM,HWRITE)
!     ###############################################################################
!
!!****  *WRITE_DIAG_INLAND_WATER_n * - Chooses the surface schemes for lakes diagnostics
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
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DIAG_SURF_ATM_n, ONLY : DIAG_SURF_ATM_t
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_CH_WATFLUX_n, ONLY : CH_WATFLUX_t
USE MODD_DIAG_WATFLUX_n, ONLY : DIAG_WATFLUX_t
USE MODD_WATFLUX_n, ONLY : WATFLUX_t
USE MODD_CH_FLAKE_n, ONLY : CH_FLAKE_t
USE MODD_DIAG_FLAKE_n, ONLY : DIAG_FLAKE_t
USE MODD_DIAG_MISC_FLAKE_n, ONLY : DIAG_MISC_FLAKE_t
USE MODD_FLAKE_n, ONLY : FLAKE_t
!
USE MODD_SURF_PAR,   ONLY : XUNDEF

USE MODI_WRITE_DIAG_WATFLUX_n
USE MODI_WRITE_DIAG_FLAKE_n
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
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DIAG_SURF_ATM_t), INTENT(INOUT) :: DGU
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(CH_WATFLUX_t), INTENT(INOUT) :: CHW
TYPE(DIAG_WATFLUX_t), INTENT(INOUT) :: DGW
TYPE(WATFLUX_t), INTENT(INOUT) :: W
TYPE(CH_FLAKE_t), INTENT(INOUT) :: CHF
TYPE(DIAG_FLAKE_t), INTENT(INOUT) :: DGF
TYPE(DIAG_MISC_FLAKE_t), INTENT(INOUT) :: DGMF
TYPE(FLAKE_t), INTENT(INOUT) :: F
!
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM ! program calling surf. schemes
 CHARACTER(LEN=3),   INTENT(IN)  :: HWRITE   ! 'PGD' : only physiographic fields are written
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!                                           ! 'ALL' : all fields are written
!
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_INLAND_WATER_N',0,ZHOOK_HANDLE)
IF (U%CWATER=='WATFLX') THEN
  CALL WRITE_DIAG_WATFLUX_n(DTCO, DGU, IOB, U, CHW, DGW, W, &
                            HPROGRAM,HWRITE)
END IF
IF (U%CWATER=='FLAKE ') THEN
  CALL WRITE_DIAG_FLAKE_n(DTCO, DGU, IOB, U, CHF, DGF, DGMF, F, &
                          HPROGRAM,HWRITE)
END IF
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_INLAND_WATER_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_DIAG_INLAND_WATER_n
