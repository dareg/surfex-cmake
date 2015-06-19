!     #########
      SUBROUTINE WRITE_TOWN_n(HPROGRAM,HWRITE)
!     ####################################
!
!!****  *WRITE_TOWN_n* - routine to write surface variables in their respective files
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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2003
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_BEM_n, ONLY : B => BEM
USE MODD_BEM_OPTION_n, ONLY : BOP => BEM_OPTIONS
USE MODD_CH_TEB_n, ONLY : CHT => CH_TEB
USE MODD_DATA_TEB_n, ONLY : DTT => DATA_TEB
USE MODD_DIAG_MISC_TEB_OPTION_n, ONLY : DGMTO => DIAG_MISC_TEB_OPTIONS
USE MODD_DIAG_SURF_ATM_n, ONLY : DGU => DIAG_SURF_ATM
USE MODD_DIAG_TEB_n, ONLY : DGT => DIAG_TEB
USE MODD_DIAG_UTCI_TEB_n, ONLY : DGUT => DIAG_UTCI_TEB
USE MODD_TEB_CANOPY_n, ONLY : TCP => TEB_CANOPY
USE MODD_TEB_GARDEN_n, ONLY : TGD => TEB_GARDEN
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TGDO => TEB_GARDEN_OPTIONS
USE MODD_TEB_GARDEN_PGD_EVOL_n, ONLY : TGDPE => TEB_GARDEN_PGD_EVOL
USE MODD_TEB_GREENROOF_n, ONLY : TGR => TEB_GREENROOF
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TGRO => TEB_GREENROOF_OPTIONS
USE MODD_TEB_GREENROOF_PGD_EVOL_n, ONLY : TGRPE => TEB_GREENROOF_PGD_EVOL
USE MODD_TEB_n, ONLY : T => TEB
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS
USE MODD_TEB_PANEL_n, ONLY : TPN => TEB_PANEL
USE MODD_TEB_VEG_n, ONLY : TVG => TEB_VEG_OPTIONS
!
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
USE MODI_WRITE_TEB_n
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=6),    INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=3),    INTENT(IN)  :: HWRITE    ! 'PREP' : does not write SBL XUNDEF fields
!                                             ! 'ALL' : all fields are written
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!*       1.     Selection of surface scheme
!               ---------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITE_TOWN_N',0,ZHOOK_HANDLE)
IF (U%CTOWN=='TEB   ') THEN
  CALL WRITE_TEB_n(B, BOP, CHT, DTT, DGMTO, DGU, DGT, DGUT, TCP, TGD, TGDO, TGDPE, TGR, TGRO, TGRPE, T, TOP, TPN, TVG, &
                   HPROGRAM,HWRITE)
END IF
IF (LHOOK) CALL DR_HOOK('WRITE_TOWN_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_TOWN_n
