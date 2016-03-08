!     #########
      SUBROUTINE PGD_NATURE (DTCO, DTI, DTZ, IG, ISS, I, UG, U, USS, HPROGRAM)
!     #############################################################
!
!!****  *PGD_NATURE* - routine to choose initialization of vegetation scheme
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
!!      Original    03/2004
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_DATA_TSZ0_n, ONLY : DATA_TSZ0_t
USE MODD_GRID_n, ONLY : GRID_t
USE MODD_SSO_n, ONLY : SSO_t
USE MODD_ISBA_n, ONLY : ISBA_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SSO_n, ONLY : SSO_t

!
USE MODI_PGD_ISBA
USE MODI_PGD_TSZ0_PAR
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTI
TYPE(DATA_TSZ0_t), INTENT(INOUT) :: DTZ
TYPE(GRID_t), INTENT(INOUT) :: IG
TYPE(SSO_t), INTENT(INOUT) :: ISS
TYPE(ISBA_t), INTENT(INOUT) :: I
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SSO_t), INTENT(INOUT) :: USS
!
 REAL, DIMENSION(:,:), POINTER :: PVEGTYPE
!
 CHARACTER(LEN=6), INTENT(IN)  :: HPROGRAM   ! program calling surf. schemes
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!                                           ! F if all parameters must be specified
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
!*       2.     Selection of surface scheme
!               ---------------------------
!
IF (LHOOK) CALL DR_HOOK('PGD_NATURE',0,ZHOOK_HANDLE)
IF (U%CNATURE=='NONE  ') THEN
  IF (LHOOK) CALL DR_HOOK('PGD_NATURE',1,ZHOOK_HANDLE)
  RETURN
ELSE IF (U%CNATURE=='FLUX  ') THEN
  IF (LHOOK) CALL DR_HOOK('PGD_NATURE',1,ZHOOK_HANDLE)
  RETURN
ELSE IF (U%CNATURE=='ISBA  ' .OR. U%CNATURE=='TSZ0') THEN
  CALL PGD_ISBA(DTCO, DTI, IG, I%O, I%P, ISS, UG, U, USS, I%M%X%XVEGTYPE, HPROGRAM)
  IF (U%CNATURE=='TSZ0') CALL PGD_TSZ0_PAR(DTZ, &
                                           HPROGRAM)
END IF
IF (LHOOK) CALL DR_HOOK('PGD_NATURE',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE PGD_NATURE
