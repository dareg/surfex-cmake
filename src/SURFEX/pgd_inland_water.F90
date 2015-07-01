!     #########
      SUBROUTINE PGD_INLAND_WATER(HPROGRAM,OECOCLIMAP,ORM_RIVER)
!     #############################################################
!
!!****  *PGD_INLAND_WATER* - routine to choose initialization of lake scheme
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
!!     B. Decharme  02/2014  Add LRM_RIVER
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
USE MODD_WATFLUX_GRID_n, ONLY : WG => WATFLUX_GRID
USE MODD_WATFLUX_n, ONLY : W => WATFLUX
!
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
USE MODI_PGD_WATFLUX
USE MODI_PGD_FLAKE
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
CHARACTER(LEN=6),    INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
LOGICAL,             INTENT(IN)  :: OECOCLIMAP
LOGICAL,             INTENT(IN)  :: ORM_RIVER ! delete river coverage (default = false)
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*       2.     Selection of surface scheme
!               ---------------------------
!
IF (LHOOK) CALL DR_HOOK('PGD_INLAND_WATER',0,ZHOOK_HANDLE)
!
IF (U%CWATER=='NONE  ') THEN
  IF (LHOOK) CALL DR_HOOK('PGD_INLAND_WATER',1,ZHOOK_HANDLE)
  RETURN
ELSE IF (U%CWATER=='FLUX  ') THEN
  IF (LHOOK) CALL DR_HOOK('PGD_INLAND_WATER',1,ZHOOK_HANDLE)
  RETURN
ELSE IF (U%CWATER=='WATFLX') THEN
  CALL PGD_WATFLUX(DTCO, U, WG, W, &
                   HPROGRAM)
ELSE IF (U%CWATER=='FLAKE ') THEN
  CALL PGD_FLAKE(HPROGRAM,OECOCLIMAP,ORM_RIVER)
END IF
!
IF (LHOOK) CALL DR_HOOK('PGD_INLAND_WATER',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE PGD_INLAND_WATER
