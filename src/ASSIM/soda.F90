! *****************************************************************************************
PROGRAM SODA
! ******************************************************************************************
USE MODI_SODA_CONTROL
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
REAL (KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SODA',0,ZHOOK_HANDLE)
!
 CALL SODA_CONTROL (ODINLINE = .FALSE.)
!
IF (LHOOK) CALL DR_HOOK('SODA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END PROGRAM SODA
