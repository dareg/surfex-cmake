!     #########
SUBROUTINE WRITE_DIAG_SEA_n (CHS, DGO, DGS, DGSI, DGU, S, U, &
                             HPROGRAM,HWRITE)
!     ###############################################################################
!
!!****  *WRITE_DIAG_SEA_n * - Chooses the surface schemes for sea diagnostics
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
!
!
USE MODD_CH_SEAFLUX_n, ONLY : CH_SEAFLUX_t
USE MODD_DIAG_OCEAN_n, ONLY : DIAG_OCEAN_t
USE MODD_DIAG_SEAFLUX_n, ONLY : DIAG_SEAFLUX_t
USE MODD_DIAG_SEAICE_n, ONLY : DIAG_SEAICE_t
USE MODD_DIAG_SURF_ATM_n, ONLY : DIAG_SURF_ATM_t
USE MODD_SEAFLUX_n, ONLY : SEAFLUX_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE MODI_WRITE_DIAG_SEAFLUX_n
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
TYPE(CH_SEAFLUX_t), INTENT(INOUT) :: CHS
TYPE(DIAG_OCEAN_t), INTENT(INOUT) :: DGO
TYPE(DIAG_SEAFLUX_t), INTENT(INOUT) :: DGS
TYPE(DIAG_SEAICE_t), INTENT(INOUT) :: DGSI
TYPE(DIAG_SURF_ATM_t), INTENT(INOUT) :: DGU
TYPE(SEAFLUX_t), INTENT(INOUT) :: S
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
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
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_SEA_N',0,ZHOOK_HANDLE)
IF (U%CSEA=='SEAFLX') THEN
  CALL WRITE_DIAG_SEAFLUX_n(CHS, DGO, DGS, DGSI, DGU, S, &
                            HPROGRAM,HWRITE)
END IF
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_SEA_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_DIAG_SEA_n
