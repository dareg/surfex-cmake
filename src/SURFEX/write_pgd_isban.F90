!     #########
      SUBROUTINE WRITE_PGD_ISBA_n (DTCO, HSELECT, U, DTI, DTZ, IG, ISS, I, HPROGRAM)
!     ####################################
!
!!****  *WRITE_PGD_ISBA_n* - routine to write pgd surface variables in their respective files
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
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_DATA_TSZ0_n, ONLY : DATA_TSZ0_t
USE MODD_GRID_n, ONLY : GRID_t
USE MODD_SSO_n, ONLY : SSO_t
USE MODD_ISBA_n, ONLY : ISBA_t
!
USE MODI_INIT_IO_SURF_n
USE MODI_WRITESURF_PGD_ISBA_n
USE MODI_END_IO_SURF_n
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
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
 CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: HSELECT
!
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTI
TYPE(DATA_TSZ0_t), INTENT(INOUT) :: DTZ
TYPE(GRID_t), INTENT(INOUT) :: IG
TYPE(SSO_t), INTENT(INOUT) :: ISS
TYPE(ISBA_t), INTENT(INOUT) :: I
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
!         Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('WRITE_PGD_ISBA_N',0,ZHOOK_HANDLE)
CALL INIT_IO_SURF_n(DTCO, U, HPROGRAM,'NATURE','ISBA  ','WRITE')
!
!*       1.     Selection of surface scheme
!               ---------------------------
!
 CALL WRITESURF_PGD_ISBA_n(HSELECT, U%CNATURE, DTI, DTZ, IG, ISS, I%O, I%P, &
                           I%M%X%XVEGTYPE, HPROGRAM)
!
!-------------------------------------------------------------------------------
!
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('WRITE_PGD_ISBA_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_PGD_ISBA_n
