!     #########
      SUBROUTINE WRITE_PGD_NATURE_n (U, &
                                     HPROGRAM)
!     ####################################
!
!!****  *WRITE_PGD_NATURE_n* - routine to write pgd surface variables in their respective files
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
USE MODD_DATA_ISBA_n, ONLY : DTI => DATA_ISBA
USE MODD_DATA_TSZ0_n, ONLY : DTZ => DATA_TSZ0
USE MODD_ISBA_GRID_n, ONLY : IG => ISBA_GRID
USE MODD_ISBA_n, ONLY : I => ISBA
!
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODI_WRITE_PGD_ISBA_n
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
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
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
!*       1.     Selection of surface scheme
!               ---------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITE_PGD_NATURE_N',0,ZHOOK_HANDLE)
IF (U%CNATURE=='ISBA' .OR. U%CNATURE=='TSZ0') THEN
  CALL WRITE_PGD_ISBA_n(DTI, DTZ, IG, I, U, &
                        HPROGRAM)
END IF
IF (LHOOK) CALL DR_HOOK('WRITE_PGD_NATURE_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_PGD_NATURE_n
