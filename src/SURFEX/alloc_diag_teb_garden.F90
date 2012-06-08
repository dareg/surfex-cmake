!     #########
    SUBROUTINE ALLOC_DIAG_TEB_GARDEN(KLU,KGROUND_LAYER,KSW)
!   ##########################################################################
!
USE MODD_TEB_GARDEN_n
USE MODD_DIAG_TEB_GARDEN_n
USE MODD_AGRI_GARDEN_n
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: KLU
INTEGER, INTENT(IN) :: KGROUND_LAYER
INTEGER, INTENT(IN) :: KSW
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! Diagnostic variables:
!
IF (LHOOK) CALL DR_HOOK('ALLOC_DIAG_TEB_GARDEN',0,ZHOOK_HANDLE)
ALLOCATE(XRI                     (KLU                     )) 
ALLOCATE(XCD                     (KLU                     )) 
ALLOCATE(XCH                     (KLU                     )) 
ALLOCATE(XRN                     (KLU                     )) 
ALLOCATE(XH                      (KLU                     )) 
ALLOCATE(XGFLUX                  (KLU                     )) 
ALLOCATE(XQS                     (KLU                     )) 
!
ALLOCATE(XLEI                    (KLU                     )) 
ALLOCATE(XLEG                    (KLU                     )) 
ALLOCATE(XLEGI                   (KLU                     )) 
ALLOCATE(XLEV                    (KLU                     )) 
ALLOCATE(XLES                    (KLU                     )) 
ALLOCATE(XLER                    (KLU                     )) 
ALLOCATE(XLETR                   (KLU                     )) 
ALLOCATE(XEVAP                   (KLU                     )) 
ALLOCATE(XDRAIN                  (KLU                     )) 
ALLOCATE(XRUNOFF                 (KLU                     )) 
ALLOCATE(XHORT                   (KLU                     )) 
ALLOCATE(XDRIP                   (KLU                     )) 
ALLOCATE(XRRVEG                  (KLU                     )) 
ALLOCATE(XMELT                   (KLU                     )) 
!
ALLOCATE(XCG                     (KLU                     )) 
ALLOCATE(XC1                     (KLU                     )) 
ALLOCATE(XC2                     (KLU                     )) 
ALLOCATE(XWGEQ                   (KLU                     )) 
ALLOCATE(XCT                     (KLU                     )) 
ALLOCATE(XRS                     (KLU                     )) 
ALLOCATE(XCDN                    (KLU                     )) 
ALLOCATE(XHU                     (KLU                     )) 
ALLOCATE(XHUG                    (KLU                     )) 
ALLOCATE(XRESTORE                (KLU                     )) 
ALLOCATE(XUSTAR                  (KLU                     )) 
ALLOCATE(XIACAN                  (KLU,SIZE(XABC)          )) 
!
ALLOCATE(XFAPAR                  (KLU                     ))
ALLOCATE(XFAPIR                  (KLU                     ))
ALLOCATE(XFAPAR_BS               (KLU                     ))
ALLOCATE(XFAPIR_BS               (KLU                     ))
ALLOCATE(XDFAPARC                (KLU                     ))
ALLOCATE(XDFAPIRC                (KLU                     ))
ALLOCATE(XDLAI_EFFC              (KLU                     ))
!
ALLOCATE(XSNOWTEMP               (KLU,TSNOW%NLAYER        )) 
ALLOCATE(XSNOWLIQ                (KLU,TSNOW%NLAYER        )) 
ALLOCATE(XSNOWDZ                 (KLU,TSNOW%NLAYER        )) 
ALLOCATE(XSNOWHMASS              (KLU                     )) 
ALLOCATE(XMELTADV                (KLU                     )) 
!
ALLOCATE(XHV                     (KLU                     ))
ALLOCATE(XALBT                   (KLU                     )) 
ALLOCATE(XEMIST                  (KLU                     )) 
ALLOCATE(XSNOWFREE_ALB           (KLU                     )) 
ALLOCATE(XSNOWFREE_ALB_VEG       (KLU                     )) 
ALLOCATE(XSNOWFREE_ALB_SOIL      (KLU                     )) 
IF (LHOOK) CALL DR_HOOK('ALLOC_DIAG_TEB_GARDEN',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE ALLOC_DIAG_TEB_GARDEN
