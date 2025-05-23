!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE DIAG_TEB_VEG_INIT_n(DK, DEK, DECK, DMK, KLU, KSNOW_LAYER)
              !     #####################
!
USE MODD_SURF_PAR, ONLY : XUNDEF
!
USE MODE_DIAG
!
USE MODD_DIAG_n, ONLY : DIAG_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
TYPE(DIAG_t), INTENT(INOUT) :: DK
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DEK
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DECK
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DMK
INTEGER, INTENT(IN) :: KLU
INTEGER, INTENT(IN) :: KSNOW_LAYER
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('DIAG_TEB_VEG_INIT',0,ZHOOK_HANDLE)
!
! Main diag
!
!general
ALLOCATE(DK%XRN   (KLU))
ALLOCATE(DK%XH    (KLU))
ALLOCATE(DK%XLE   (KLU))
ALLOCATE(DK%XLEI  (KLU))
ALLOCATE(DK%XGFLUX(KLU))
ALLOCATE(DK%XEVAP (KLU))
ALLOCATE(DK%XSUBL (KLU))
!
DK%XRN   (:) = XUNDEF
DK%XH    (:) = XUNDEF
DK%XLE   (:) = XUNDEF
DK%XLEI  (:) = XUNDEF
DK%XGFLUX(:) = XUNDEF
DK%XEVAP (:) = XUNDEF
DK%XSUBL (:) = XUNDEF
!
!
ALLOCATE(DK%XRI(KLU))
!
DK%XRI(:) = XUNDEF
!
!
ALLOCATE(DK%XCD (KLU))
ALLOCATE(DK%XCH (KLU))
ALLOCATE(DK%XZ0 (KLU))
ALLOCATE(DK%XZ0H(KLU))
!
DK%XCD (:) = XUNDEF
DK%XCH (:) = XUNDEF
DK%XZ0 (:) = XUNDEF
DK%XZ0H(:) = XUNDEF
!
!
ALLOCATE(DK%XZ0EFF(KLU))
ALLOCATE(DK%XCDN  (KLU))
ALLOCATE(DK%XHUG  (KLU))
ALLOCATE(DK%XHU   (KLU))
ALLOCATE(DK%XQS   (KLU))
!
DK%XZ0EFF(:) = XUNDEF
DK%XCDN  (:) = XUNDEF
DK%XHUG  (:) = XUNDEF
DK%XHU   (:) = XUNDEF
DK%XQS   (:) = XUNDEF
!
!
ALLOCATE(DK%XTS   (KLU))
ALLOCATE(DK%XTSRAD(KLU))
ALLOCATE(DK%XALBT (KLU))
!
DK%XTS   (:) = XUNDEF
DK%XTSRAD(:) = XUNDEF
DK%XALBT (:) = XUNDEF
!
!---------------------------
!
! Evaporation diag
!
! general
ALLOCATE(DEK%XLEG (KLU))
ALLOCATE(DEK%XLEGI(KLU))
ALLOCATE(DEK%XLEV (KLU))
ALLOCATE(DEK%XLES (KLU))
ALLOCATE(DEK%XLESL(KLU))
ALLOCATE(DEK%XLER (KLU))
ALLOCATE(DEK%XLETR(KLU))
!
DEK%XLEG (:) = XUNDEF
DEK%XLEGI(:) = XUNDEF
DEK%XLEV (:) = XUNDEF
DEK%XLES (:) = XUNDEF
DEK%XLESL(:) = XUNDEF
DEK%XLER (:) = XUNDEF
DEK%XLETR(:) = XUNDEF
!
!
ALLOCATE(DEK%XDRAIN (KLU))
ALLOCATE(DEK%XRUNOFF(KLU))
ALLOCATE(DEK%XHORT  (KLU))
ALLOCATE(DEK%XQSB   (KLU))
ALLOCATE(DEK%XIRRIG_FLUX(KLU))
!
DEK%XDRAIN (:) = XUNDEF
DEK%XRUNOFF(:) = XUNDEF
DEK%XHORT  (:) = XUNDEF
DEK%XQSB   (:) = XUNDEF
DEK%XIRRIG_FLUX(:) = XUNDEF
!
!
ALLOCATE(DEK%XDRIP (KLU))
ALLOCATE(DEK%XRRVEG(KLU))
!
DEK%XDRIP (:) = XUNDEF
DEK%XRRVEG(:) = XUNDEF
!
!
ALLOCATE(DEK%XMELT    (KLU))
ALLOCATE(DEK%XMELTADV (KLU))
ALLOCATE(DEK%XRESTORE (KLU))
ALLOCATE(DEK%XSNDRIFT (KLU))
ALLOCATE(DEK%XSWNET_N (KLU))
ALLOCATE(DEK%XSWNET_NS(KLU))
ALLOCATE(DEK%XLWNET_N (KLU))
!
DEK%XMELT    (:) = XUNDEF
DEK%XMELTADV (:) = XUNDEF
DEK%XRESTORE (:) = XUNDEF
DEK%XSNDRIFT (:) = XUNDEF
DEK%XSWNET_N (:) = XUNDEF
DEK%XSWNET_NS(:) = XUNDEF
DEK%XLWNET_N (:) = XUNDEF
!
!
ALLOCATE(DEK%XLE_FLOOD (KLU))
ALLOCATE(DEK%XLEI_FLOOD(KLU))
ALLOCATE(DEK%XPFLOOD   (KLU))
ALLOCATE(DEK%XIFLOOD   (KLU))
!
DEK%XLE_FLOOD (:) = XUNDEF
DEK%XLEI_FLOOD(:) = XUNDEF
DEK%XPFLOOD   (:) = XUNDEF
DEK%XIFLOOD   (:) = XUNDEF
!
!
ALLOCATE(DEK%XGPP      (KLU))
ALLOCATE(DEK%XRESP_ECO (KLU))
ALLOCATE(DEK%XRESP_AUTO(KLU))
!
DEK%XGPP      (:) = XUNDEF
DEK%XRESP_ECO (:) = XUNDEF
DEK%XRESP_AUTO(:) = XUNDEF
!
!
ALLOCATE(DEK%XRN_SN_FR   (KLU))
ALLOCATE(DEK%XH_SN_FR    (KLU))
ALLOCATE(DEK%XLEI_SN_FR  (KLU))
ALLOCATE(DEK%XLE_SN_FR   (KLU))
ALLOCATE(DEK%XGFLUX_SN_FR(KLU))
ALLOCATE(DEK%XLEG_SN_FR  (KLU))
ALLOCATE(DEK%XLEGI_SN_FR (KLU))
ALLOCATE(DEK%XLEV_SN_FR  (KLU))
ALLOCATE(DEK%XLETR_SN_FR (KLU))
ALLOCATE(DEK%XUSTAR_SN_FR(KLU))
ALLOCATE(DEK%XLER_SN_FR  (KLU))
!
DEK%XRN_SN_FR   (:) = XUNDEF
DEK%XH_SN_FR    (:) = XUNDEF
DEK%XLEI_SN_FR  (:) = XUNDEF
DEK%XLE_SN_FR   (:) = XUNDEF
DEK%XGFLUX_SN_FR(:) = XUNDEF
DEK%XLEG_SN_FR  (:) = XUNDEF
DEK%XLEGI_SN_FR (:) = XUNDEF
DEK%XLEV_SN_FR  (:) = XUNDEF
DEK%XLETR_SN_FR (:) = XUNDEF
DEK%XUSTAR_SN_FR(:) = XUNDEF
DEK%XLER_SN_FR  (:) = XUNDEF
!
!
ALLOCATE(DECK%XDRAIN (KLU))
ALLOCATE(DECK%XRUNOFF(KLU))
!
DECK%XDRAIN (:) = 0.
DECK%XRUNOFF(:) = 0.
!
!--------------------------------------
!
! Misc diag
!
ALLOCATE(DMK%XC1  (KLU))
ALLOCATE(DMK%XC2  (KLU))
ALLOCATE(DMK%XWGEQ(KLU))
ALLOCATE(DMK%XCG  (KLU))
ALLOCATE(DMK%XCT  (KLU))
ALLOCATE(DMK%XRS  (KLU))
ALLOCATE(DMK%XHV  (KLU))
ALLOCATE(DMK%XGRNDFLUX (KLU))
ALLOCATE(DMK%XSNOWTEMP (KLU,KSNOW_LAYER))
ALLOCATE(DMK%XSNOWHMASS(KLU))
ALLOCATE(DMK%XSNOWLIQ  (KLU,KSNOW_LAYER))
ALLOCATE(DMK%XSNOWDZ   (KLU,KSNOW_LAYER))
ALLOCATE(DMK%XSRSFC    (KLU))
ALLOCATE(DMK%XRRSFC    (KLU))
ALLOCATE(DMK%XRNSNOW   (KLU))
ALLOCATE(DMK%XHSNOW    (KLU))
ALLOCATE(DMK%XGFLUXSNOW(KLU))
ALLOCATE(DMK%XHPSNOW   (KLU))
ALLOCATE(DMK%XUSTARSNOW(KLU))
ALLOCATE(DMK%XCDSNOW   (KLU))
ALLOCATE(DMK%XCHSNOW   (KLU))
!
DMK%XC1        = XUNDEF
DMK%XC2        = XUNDEF
DMK%XWGEQ      = XUNDEF
DMK%XCG        = XUNDEF
DMK%XCT        = XUNDEF
DMK%XRS        = XUNDEF
DMK%XHV        = XUNDEF
DMK%XGRNDFLUX  = XUNDEF
DMK%XSNOWTEMP  = XUNDEF
DMK%XSNOWHMASS = XUNDEF
DMK%XSNOWLIQ   = XUNDEF
DMK%XSNOWDZ    = XUNDEF
DMK%XSRSFC     = XUNDEF
DMK%XRRSFC     = XUNDEF
DMK%XRNSNOW    = XUNDEF
DMK%XHSNOW     = XUNDEF
DMK%XGFLUXSNOW = XUNDEF
DMK%XHPSNOW    = XUNDEF
DMK%XUSTARSNOW = XUNDEF
DMK%XCDSNOW    = XUNDEF
DMK%XCHSNOW    = XUNDEF
!
! LPROSNOW = .FALSE. for GARDEN GREENROOF
ALLOCATE(DMK%XSNOWDEND(0,0))
!
IF (LHOOK) CALL DR_HOOK('DIAG_TEB_VEG_INIT',1,ZHOOK_HANDLE)
!
END SUBROUTINE DIAG_TEB_VEG_INIT_n
