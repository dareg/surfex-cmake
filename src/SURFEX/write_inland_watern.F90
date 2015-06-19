!     #########
      SUBROUTINE WRITE_INLAND_WATER_n(HPROGRAM,HWRITE)
!     ####################################
!
!!****  *WRITE_INLAND_WATER_n* - routine to write surface variables in their respective files
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
USE MODD_CH_WATFLUX_n, ONLY : CHW => CH_WATFLUX
USE MODD_WATFLUX_n, ONLY : W => WATFLUX
USE MODD_WATFLUX_SBL_n, ONLY : WSB => WATFLUX_SBL
!
USE MODD_CH_FLAKE_n, ONLY : CHF => CH_FLAKE
USE MODD_DIAG_MISC_FLAKE_n, ONLY : DGMF => DIAG_MISC_FLAKE
USE MODD_DIAG_SURF_ATM_n, ONLY : DGU => DIAG_SURF_ATM
USE MODD_FLAKE_n, ONLY : F => FLAKE
USE MODD_FLAKE_SBL_n, ONLY : FSB => FLAKE_SBL
!
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
USE MODI_WRITE_WATFLUX_n
USE MODI_WRITE_FLAKE_n
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
IF (LHOOK) CALL DR_HOOK('WRITE_INLAND_WATER_N',0,ZHOOK_HANDLE)
IF (U%CWATER=='WATFLX') THEN
  CALL WRITE_WATFLUX_n(CHW, DGU, W, WSB, &
                       HPROGRAM,HWRITE)
ELSE IF (U%CWATER=='FLAKE ') THEN
  CALL WRITE_FLAKE_n(CHF, DGMF, DGU, F, FSB, &
                     HPROGRAM,HWRITE)
END IF
IF (LHOOK) CALL DR_HOOK('WRITE_INLAND_WATER_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_INLAND_WATER_n
