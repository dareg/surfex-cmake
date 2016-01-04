!     #########
    SUBROUTINE ALLOC_DIAG_TEB_GREENROOF (DGGR, DGMGR, DGEGR, TGR, &
                                         KLU,KLAYER_GR,KSW)
!   ##########################################################################
!
!
!
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_ISBA_n, ONLY : DIAG_ISBA_t
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_PROG_t
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!
TYPE(DIAG_ISBA_t), INTENT(INOUT) :: DGGR
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEGR
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DGMGR
TYPE(TEB_VEG_PROG_t), INTENT(INOUT) :: TGR
!
INTEGER, INTENT(IN) :: KLU
INTEGER, INTENT(IN) :: KLAYER_GR
INTEGER, INTENT(IN) :: KSW
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! Diagnostic variables:
!
IF (LHOOK) CALL DR_HOOK('ALLOC_DIAG_TEB_GREENROOF',0,ZHOOK_HANDLE)
ALLOCATE(DGGR%XRI                     (KLU,  1)) 
ALLOCATE(DGGR%XCD                     (KLU,  1)) 
ALLOCATE(DGGR%XCH                     (KLU,  1)) 
ALLOCATE(DGGR%XRN                     (KLU,  1)) 
ALLOCATE(DGGR%XH                      (KLU,  1)) 
ALLOCATE(DGGR%XGFLUX                  (KLU,  1)) 
ALLOCATE(DGGR%XQS                     (KLU,  1)) 
!
ALLOCATE(DGGR%XLEI                    (KLU,  1))
!
ALLOCATE(DGEGR%XLEG                    (KLU,  1)) 
ALLOCATE(DGEGR%XLEGI                   (KLU,  1)) 
ALLOCATE(DGEGR%XLEV                    (KLU,  1)) 
ALLOCATE(DGEGR%XLES                    (KLU,  1)) 
ALLOCATE(DGEGR%XLER                    (KLU,  1)) 
ALLOCATE(DGEGR%XLETR                   (KLU,  1)) 
ALLOCATE(DGEGR%XEVAP                   (KLU,  1)) 
ALLOCATE(DGEGR%XDRAIN                  (KLU,  1)) 
ALLOCATE(DGEGR%XRUNOFF                 (KLU,  1)) 
ALLOCATE(DGEGR%XHORT                   (KLU,  1)) 
ALLOCATE(DGEGR%XDRIP                   (KLU,  1)) 
ALLOCATE(DGEGR%XRRVEG                  (KLU,  1)) 
ALLOCATE(DGEGR%XMELT                   (KLU,  1)) 
!
ALLOCATE(DGMGR%XSNOWTEMP               (KLU,TGR%CUR%TSNOW%NLAYER, 1)) 
ALLOCATE(DGMGR%XSNOWLIQ                (KLU,TGR%CUR%TSNOW%NLAYER, 1)) 
!
ALLOCATE(DGMGR%XHV                     (KLU,  1))
ALLOCATE(DGMGR%XALBT                   (KLU,  1)) 
IF (LHOOK) CALL DR_HOOK('ALLOC_DIAG_TEB_GREENROOF',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE ALLOC_DIAG_TEB_GREENROOF
