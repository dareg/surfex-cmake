!     #########
SUBROUTINE WRITE_DIAG_NATURE_n (CHI, DGEI, DGI, DGMI, DGU, DST, GB, I, U, &
                                HPROGRAM,HWRITE)
!     ###############################################################################
!
!!****  *WRITE_DIAG_NATURE_n * - Chooses the surface schemes for diagnostics over
!!    natural continental parts
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
USE MODD_CH_ISBA_n, ONLY : CH_ISBA_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_ISBA_n, ONLY : DIAG_ISBA_t
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t
USE MODD_DIAG_SURF_ATM_n, ONLY : DIAG_SURF_ATM_t
USE MODD_DST_n, ONLY : DST_t
USE MODD_GR_BIOG_n, ONLY : GR_BIOG_t
USE MODD_ISBA_n, ONLY : ISBA_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE MODI_WRITE_DIAG_ISBA_n
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
TYPE(CH_ISBA_t), INTENT(INOUT) :: CHI
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEI
TYPE(DIAG_ISBA_t), INTENT(INOUT) :: DGI
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DGMI
TYPE(DIAG_SURF_ATM_t), INTENT(INOUT) :: DGU
TYPE(DST_t), INTENT(INOUT) :: DST
TYPE(GR_BIOG_t), INTENT(INOUT) :: GB
TYPE(ISBA_t), INTENT(INOUT) :: I
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
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_NATURE_N',0,ZHOOK_HANDLE)
IF (U%CNATURE=='ISBA  ') THEN
  CALL WRITE_DIAG_ISBA_n(CHI, DGEI, DGI, DGMI, DGU, DST, GB, I, &
                         HPROGRAM,HWRITE)
END IF
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_NATURE_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_DIAG_NATURE_n
