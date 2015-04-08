!     #########
      SUBROUTINE CPL_GCM_n(KI,PRAIN,PSNOW,PZ0,PZ0H,PQSURF)
!     ####################################################
!
!!****  *CPL_GCM_n* - Save some phisical fields for ARPEGE/ALADIN run
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!
!!    EXTERNAL
!!    --------
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
!!      Original    04/2013
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER,             INTENT(IN)           :: KI       ! number of points
!
REAL, DIMENSION(KI), INTENT(IN), OPTIONAL :: PRAIN    ! total rainfall rate (kg/m2/s)
REAL, DIMENSION(KI), INTENT(IN), OPTIONAL :: PSNOW    ! total snowfall rate (kg/m2/s)
REAL, DIMENSION(KI), INTENT(IN), OPTIONAL :: PZ0      ! roughness length for momentum         (m)
REAL, DIMENSION(KI), INTENT(IN), OPTIONAL :: PZ0H     ! roughness length for heat             (m)
REAL, DIMENSION(KI), INTENT(IN), OPTIONAL :: PQSURF   ! specific humidity at surface          (kg/kg)
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CPL_GCM_N',0,ZHOOK_HANDLE)
!
IF(PRESENT(PRAIN )) U%XRAIN (:) = PRAIN (:)
IF(PRESENT(PSNOW )) U%XSNOW (:) = PSNOW (:)
IF(PRESENT(PZ0   )) U%XZ0   (:) = PZ0   (:)
IF(PRESENT(PZ0H  )) U%XZ0H  (:) = PZ0H  (:)
IF(PRESENT(PQSURF)) U%XQSURF(:) = PQSURF(:)
!
IF (LHOOK) CALL DR_HOOK('CPL_GCM_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE CPL_GCM_n
