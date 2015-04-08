!     #########
      SUBROUTINE DIAG_IDEAL_INIT_n(KLU,KSW)
!     #####################
!
!!****  *DIAG_IDEAL_INIT_n* - routine to initialize IDEAL diagnostic variables
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
!!      P. Le Moigne   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2009 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_DIAG_IDEAL_n, ONLY : DGL => DIAG_IDEAL
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
INTEGER, INTENT(IN) :: KLU   ! size of arrays
INTEGER, INTENT(IN) :: KSW   ! spectral bands
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
!* surface energy budget
!
IF (LHOOK) CALL DR_HOOK('DIAG_IDEAL_INIT_N',0,ZHOOK_HANDLE)
IF (DGL%LSURF_BUDGET) THEN
  ALLOCATE(DGL%XRN     (KLU))
  ALLOCATE(DGL%XH      (KLU))
  ALLOCATE(DGL%XLE     (KLU))
  ALLOCATE(DGL%XGFLUX  (KLU))
  ALLOCATE(DGL%XSWD    (KLU))
  ALLOCATE(DGL%XSWU    (KLU))
  ALLOCATE(DGL%XSWBD   (KLU,KSW))
  ALLOCATE(DGL%XSWBU   (KLU,KSW))
  ALLOCATE(DGL%XLWD    (KLU))
  ALLOCATE(DGL%XLWU    (KLU))
  ALLOCATE(DGL%XFMU    (KLU))
  ALLOCATE(DGL%XFMV    (KLU))
  !
  DGL%XRN      = XUNDEF
  DGL%XH       = XUNDEF
  DGL%XLE      = XUNDEF
  DGL%XGFLUX   = XUNDEF
  DGL%XSWD     = XUNDEF
  DGL%XSWU     = XUNDEF
  DGL%XSWBD    = XUNDEF
  DGL%XSWBU    = XUNDEF
  DGL%XLWD     = XUNDEF
  DGL%XLWU     = XUNDEF
  DGL%XFMU     = XUNDEF
  DGL%XFMV     = XUNDEF
ELSE
  ALLOCATE(DGL%XRN     (0))
  ALLOCATE(DGL%XH      (0))
  ALLOCATE(DGL%XLE     (0))
  ALLOCATE(DGL%XGFLUX  (0))
  ALLOCATE(DGL%XSWD    (0))
  ALLOCATE(DGL%XSWU    (0))
  ALLOCATE(DGL%XLWD    (0))
  ALLOCATE(DGL%XLWU    (0))
  ALLOCATE(DGL%XSWBD   (0,0))
  ALLOCATE(DGL%XSWBU   (0,0))
  ALLOCATE(DGL%XFMU    (0))
  ALLOCATE(DGL%XFMV    (0))
END IF
!
!* parameters at 2m
!
IF (DGL%N2M>=1) THEN
  ALLOCATE(DGL%XRI     (KLU))
  ALLOCATE(DGL%XT2M    (KLU))
  ALLOCATE(DGL%XQ2M    (KLU))
  ALLOCATE(DGL%XHU2M   (KLU))
  ALLOCATE(DGL%XZON10M (KLU))
  ALLOCATE(DGL%XMER10M (KLU))
  !
  DGL%XRI      = XUNDEF
  DGL%XT2M     = XUNDEF
  DGL%XQ2M     = XUNDEF
  DGL%XHU2M    = XUNDEF
  DGL%XZON10M  = XUNDEF
  DGL%XMER10M  = XUNDEF
ELSE
  ALLOCATE(DGL%XRI     (0))
  ALLOCATE(DGL%XT2M    (0))
  ALLOCATE(DGL%XQ2M    (0))
  ALLOCATE(DGL%XHU2M   (0))
  ALLOCATE(DGL%XZON10M (0))
  ALLOCATE(DGL%XMER10M (0))
END IF
!
!* transfer coefficients
!
IF (DGL%LCOEF) THEN
  ALLOCATE(DGL%XCD     (KLU))
  ALLOCATE(DGL%XCH     (KLU))
  ALLOCATE(DGL%XCE     (KLU))
  ALLOCATE(DGL%XZ0     (KLU))
  ALLOCATE(DGL%XZ0H    (KLU))
  !
  DGL%XCD      = XUNDEF
  DGL%XCH      = XUNDEF
  DGL%XCE      = XUNDEF
  DGL%XZ0      = XUNDEF
  DGL%XZ0H     = XUNDEF
ELSE
  ALLOCATE(DGL%XCD     (0))
  ALLOCATE(DGL%XCH     (0))
  ALLOCATE(DGL%XCE     (0))
  ALLOCATE(DGL%XZ0     (0))
  ALLOCATE(DGL%XZ0H    (0))
END IF
!
!
!* surface humidity
!
IF (DGL%LSURF_VARS) THEN
  ALLOCATE(DGL%XQS     (KLU))
  !
  DGL%XQS      = XUNDEF
ELSE
  ALLOCATE(DGL%XQS     (0))
END IF
IF (LHOOK) CALL DR_HOOK('DIAG_IDEAL_INIT_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_IDEAL_INIT_n
