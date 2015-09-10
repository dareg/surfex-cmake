!     #########
      SUBROUTINE WRITE_PGD_TEB_n (DTCO, DGU, IOB, U, &
                                  B, BOP, BDD, DTB, DTGD, DTGR, DTT, DGCT, DGMT, &
                                  T, TGD, TGDO, TGDP, TGDPE, TGR, TGRPE, TGRO, TGRP, &
                                  TG, TIR, TOP, TVG, &
                                  HPROGRAM)
!     ####################################
!
!!****  *WRITE_PGD_TEB_n* - routine to write pgd surface variables in their respective files
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      B. Decharme   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/2011
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DIAG_SURF_ATM_n, ONLY : DIAG_SURF_ATM_t
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_BEM_n, ONLY : BEM_t
USE MODD_BEM_OPTION_n, ONLY : BEM_OPTIONS_t
USE MODD_BLD_DESCRIPTION_n, ONLY : BLD_DESC_t
USE MODD_DATA_BEM_n, ONLY : DATA_BEM_t
USE MODD_DATA_TEB_GARDEN_n, ONLY : DATA_TEB_GARDEN_t
USE MODD_DATA_TEB_GREENROOF_n, ONLY : DATA_TEB_GREENROOF_t
USE MODD_DATA_TEB_n, ONLY : DATA_TEB_t
USE MODD_DIAG_CUMUL_TEB_n, ONLY : DIAG_CUMUL_TEB_t
USE MODD_DIAG_MISC_TEB_n, ONLY : DIAG_MISC_TEB_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_TEB_GARDEN_n, ONLY : TEB_GARDEN_t
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TEB_GARDEN_OPTIONS_t
USE MODD_TEB_GARDEN_PGD_n, ONLY : TEB_GARDEN_PGD_t
USE MODD_TEB_GARDEN_PGD_EVOL_n, ONLY : TEB_GARDEN_PGD_EVOL_t
USE MODD_TEB_GREENROOF_n, ONLY : TEB_GREENROOF_t
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TEB_GREENROOF_OPTIONS_t
USE MODD_TEB_GREENROOF_PGD_n, ONLY : TEB_GREENROOF_PGD_t
USE MODD_TEB_GREENROOF_PGD_EVOL_n, ONLY : TEB_GREENROOF_PGD_EVOL_t
USE MODD_TEB_GRID_n, ONLY : TEB_GRID_t
USE MODD_TEB_IRRIG_n, ONLY : TEB_IRRIG_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_OPTIONS_t
!
USE MODI_INIT_IO_SURF_n
USE MODI_WRITESURF_PGD_TEB_n
USE MODI_END_IO_SURF_n
USE MODI_GOTO_WRAPPER_TEB_PATCH
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DIAG_SURF_ATM_t), INTENT(INOUT) :: DGU
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(BEM_t), INTENT(INOUT) :: B
TYPE(BEM_OPTIONS_t), INTENT(INOUT) :: BOP
TYPE(BLD_DESC_t), INTENT(INOUT) :: BDD
TYPE(DATA_BEM_t), INTENT(INOUT) :: DTB
TYPE(DATA_TEB_GARDEN_t), INTENT(INOUT) :: DTGD
TYPE(DATA_TEB_GREENROOF_t), INTENT(INOUT) :: DTGR
TYPE(DATA_TEB_t), INTENT(INOUT) :: DTT
TYPE(DIAG_CUMUL_TEB_t), INTENT(INOUT) :: DGCT
TYPE(DIAG_MISC_TEB_t), INTENT(INOUT) :: DGMT
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(TEB_GARDEN_t), INTENT(INOUT) :: TGD
TYPE(TEB_GARDEN_OPTIONS_t), INTENT(INOUT) :: TGDO
TYPE(TEB_GARDEN_PGD_t), INTENT(INOUT) :: TGDP
TYPE(TEB_GARDEN_PGD_EVOL_t), INTENT(INOUT) :: TGDPE
TYPE(TEB_GREENROOF_t), INTENT(INOUT) :: TGR
TYPE(TEB_GREENROOF_PGD_EVOL_t), INTENT(INOUT) :: TGRPE
TYPE(TEB_GREENROOF_OPTIONS_t), INTENT(INOUT) :: TGRO
TYPE(TEB_GREENROOF_PGD_t), INTENT(INOUT) :: TGRP
TYPE(TEB_GRID_t), INTENT(INOUT) :: TG
TYPE(TEB_IRRIG_t), INTENT(INOUT) :: TIR
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_VEG_OPTIONS_t), INTENT(INOUT) :: TVG
!
 CHARACTER(LEN=6),    INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!                                             ! 'ALL' : all fields are written
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
!
!         Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('WRITE_PGD_TEB_N',0,ZHOOK_HANDLE)
CALL INIT_IO_SURF_n(DTCO, DGU, IOB, U, &
                     HPROGRAM,'TOWN  ','TEB   ','WRITE')
!
!*       1.     Selection of surface scheme
!               ---------------------------
!
 CALL GOTO_WRAPPER_TEB_PATCH(B, DGCT, DGMT, T, TGD, TGDPE, TGR, TGRPE, 1)
 CALL WRITESURF_PGD_TEB_n(BOP, BDD, DTB, DTGD, DTGR, DTT, DGU, IOB, U, TGDO, TGDP, &
                                      TGRO, TGRP, TG, TIR, TOP, TVG, &
                          HPROGRAM)
!
!-------------------------------------------------------------------------------
!
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('WRITE_PGD_TEB_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_PGD_TEB_n
