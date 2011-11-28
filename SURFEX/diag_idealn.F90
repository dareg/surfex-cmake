!     #########
SUBROUTINE DIAG_IDEAL_n(HPROGRAM, PH, PLE, PRN, PGFLUX)
!     ###############################################################################
!
!!****  *DIAG_IDEAL_n * - Stores IDEAL_n diagnostics
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
!!     P. Le Moigne 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2009
!!------------------------------------------------------------------
!

!
USE MODD_SURF_PAR,    ONLY : XUNDEF
USE MODD_DIAG_IDEAL_n, ONLY : LSURF_BUDGET, XH, XLE, XRN, XGFLUX
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM ! program calling surf. schemes
!
REAL, DIMENSION(:), INTENT(OUT) :: PH       ! Sensible heat flux  (W/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PLE      ! Latent heat flux    (W/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PRN      ! net flux    (W/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PGFLUX   ! net flux    (W/m2)
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DIAG_IDEAL_N',0,ZHOOK_HANDLE)
IF (LSURF_BUDGET) THEN
  PH       = XH
  PLE      = XLE
  PRN      = XRN
  PGFLUX   = XGFLUX
END IF
IF (LHOOK) CALL DR_HOOK('DIAG_IDEAL_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_IDEAL_n
