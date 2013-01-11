!     #########
    SUBROUTINE ALLOCATE_TEB_GREENROOF_PGD(KLU,KVEGTYPE,KLAYER_GR, KDIMTAB)  
!   ##########################################################################
!
USE MODD_TEB_GREENROOF_n
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: KLU
INTEGER, INTENT(IN) :: KVEGTYPE
INTEGER, INTENT(IN) :: KLAYER_GR
INTEGER, INTENT(IN) :: KDIMTAB
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! Mask and number of grid elements containing tiles:
!
IF (LHOOK) CALL DR_HOOK('ALLOCATE_TEB_GREENROOF_PGD',0,ZHOOK_HANDLE)
ALLOCATE(XVEGTYPE                (KLU,KVEGTYPE            ))
!
!-------------------------------------------------------------------------------
!
! Input Parameters:
!
! - vegetation + bare soil:
!
ALLOCATE(XZ0_O_Z0H               (KLU                     )) 
ALLOCATE(XEMIS                   (KLU                     )) 
ALLOCATE(XZ0                     (KLU                     )) 
!
! - vegetation:
!
ALLOCATE(XALBNIR_VEG             (KLU                     )) 
ALLOCATE(XALBVIS_VEG             (KLU                     )) 
ALLOCATE(XALBUV_VEG              (KLU                     )) 
!
! - vegetation: default option (Jarvis) and general parameters:
!
ALLOCATE(XVEG                    (KLU                     )) 
ALLOCATE(XWRMAX_CF               (KLU                     )) 
ALLOCATE(XGAMMA                  (KLU                     )) 
ALLOCATE(XCV                     (KLU                     )) 
ALLOCATE(XRGL                    (KLU                     )) 
ALLOCATE(XRSMIN                  (KLU                     )) 
ALLOCATE(XROOTFRAC               (KLU,KLAYER_GR           ))
ALLOCATE(NWG_LAYER               (KLU                     ))
ALLOCATE(XDROOT                  (KLU                     ))
ALLOCATE(XDG2                    (KLU                     ))
!
!-------------------------------------------------------------------------------
!
! - LAI: Physiographic prescribed field (YPHOTO='NON', 'AGS' or 'LST')
!
ALLOCATE(XLAI                    (KLU                    ))
!
!-------------------------------------------------------------------------------
!
! - vegetation: Ags parameters ('AGS', 'LAI', 'AST', 'LST', 'NIT' options)
!
ALLOCATE(XBSLAI                  (KLU                     )) 
ALLOCATE(XLAIMIN                 (KLU                     )) 
ALLOCATE(XSEFOLD                 (KLU                     )) 
ALLOCATE(XH_TREE                 (KLU                     )) 
ALLOCATE(XANF                    (KLU                     ))
ALLOCATE(XGMES                   (KLU                     ))
ALLOCATE(XRE25                   (KLU                     ))
!
!-------------------------------------------------------------------------------
!
! - vegetation: Ags Stress parameters ('AST', 'LST', 'NIT' options)
!
ALLOCATE(LSTRESS                 (KLU                     )) 
ALLOCATE(XF2I                    (KLU                     )) 
ALLOCATE(XGC                     (KLU                     )) 
ALLOCATE(XAH                     (KLU                     )) 
ALLOCATE(XBH                     (KLU                     )) 
ALLOCATE(XDMAX                   (KLU                     )) 
!
!-------------------------------------------------------------------------------
!
! - vegetation: Ags Nitrogen-model parameters ('NIT' option)
!
ALLOCATE(XCE_NITRO               (KLU                     )) 
ALLOCATE(XCF_NITRO               (KLU                     )) 
ALLOCATE(XCNA_NITRO              (KLU                     )) 
!
!-------------------------------------------------------------------------------
!
! - soil: primary parameters
!
ALLOCATE(XOM_GR                  (KLU,KLAYER_GR           ))  
ALLOCATE(XSAND_GR                (KLU,KLAYER_GR           ))  
ALLOCATE(XCLAY_GR                (KLU,KLAYER_GR           ))  
ALLOCATE(XRUNOFFB_GR             (KLU                     ))  
ALLOCATE(XWDRAIN_GR              (KLU                     ))  
ALLOCATE(XTAUICE                 (KLU                     )) 
ALLOCATE(XGAMMAT                 (KLU                     )) 
ALLOCATE(XDG                     (KLU,KLAYER_GR           )) 
ALLOCATE(XRUNOFFD                (KLU                     )) 
!
!-------------------------------------------------------------------------------
!
! - SGH scheme
!                                   
ALLOCATE(XD_ICE                  (KLU                     )) 
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ALLOCATE_TEB_GREENROOF_PGD',1,ZHOOK_HANDLE)
!
END SUBROUTINE ALLOCATE_TEB_GREENROOF_PGD
