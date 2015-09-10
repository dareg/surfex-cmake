!     #########
      SUBROUTINE WRITE_NATURE_n (DTCO, DGU, IOB, U, CHI, DGEI, DGI, DGMI, I, DST, ICP, &
                                 HPROGRAM,HWRITE,OLAND_USE)
!     ####################################
!
!!****  *WRITE_NATURE_n* - routine to write surface variables in their respective files
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
!       Modified    07/2011, B.Decharme: land_use key for writing semi-prognostic variables
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DIAG_SURF_ATM_n, ONLY : DIAG_SURF_ATM_t
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_CH_ISBA_n, ONLY : CH_ISBA_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_ISBA_n, ONLY : DIAG_ISBA_t
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t
USE MODD_ISBA_n, ONLY : ISBA_t
USE MODD_DST_n, ONLY : DST_t
USE MODD_ISBA_CANOPY_n, ONLY : ISBA_CANOPY_t
!
USE MODI_WRITE_ISBA_n
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
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DIAG_SURF_ATM_t), INTENT(INOUT) :: DGU
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(CH_ISBA_t), INTENT(INOUT) :: CHI
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEI
TYPE(DIAG_ISBA_t), INTENT(INOUT) :: DGI
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DGMI
TYPE(ISBA_t), INTENT(INOUT) :: I
TYPE(DST_t), INTENT(INOUT) :: DST
TYPE(ISBA_CANOPY_t), INTENT(INOUT) :: ICP
!
 CHARACTER(LEN=6),    INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=3),    INTENT(IN)  :: HWRITE    ! 'PREP' : does not write SBL XUNDEF fields
!                                             ! 'ALL' : all fields are written
LOGICAL,             INTENT(IN)  :: OLAND_USE !
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!*       1.     Selection of surface scheme
!               ---------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITE_NATURE_N',0,ZHOOK_HANDLE)
IF (U%CNATURE=='ISBA' .OR. U%CNATURE=='TSZ0') THEN
  CALL WRITE_ISBA_n(DTCO, DGU, IOB, U, CHI, DGEI, DGI, DGMI, I, DST, ICP, &
                    HPROGRAM,HWRITE,OLAND_USE)
END IF
IF (LHOOK) CALL DR_HOOK('WRITE_NATURE_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_NATURE_n
