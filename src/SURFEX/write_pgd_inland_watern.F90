!     #########
      SUBROUTINE WRITE_PGD_INLAND_WATER_n (U, &
                                           HPROGRAM)
!     ####################################
!
!!****  *WRITE_PGD_INLAND_WATER_n* - routine to write pgd surface variables in their respective files
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
!
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODI_WRITE_PGD_WATFLUX_n
USE MODI_WRITE_PGD_FLAKE_n
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
IF (LHOOK) CALL DR_HOOK('WRITE_PGD_INLAND_WATER_N',0,ZHOOK_HANDLE)
IF (U%CWATER=='WATFLX') THEN
  CALL WRITE_PGD_WATFLUX_n(HPROGRAM)
ELSE IF (U%CWATER=='FLAKE ') THEN
  CALL WRITE_PGD_FLAKE_n(HPROGRAM)
END IF
IF (LHOOK) CALL DR_HOOK('WRITE_PGD_INLAND_WATER_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_PGD_INLAND_WATER_n
