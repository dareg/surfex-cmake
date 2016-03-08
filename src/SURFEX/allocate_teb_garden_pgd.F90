!     #########
    SUBROUTINE ALLOCATE_TEB_GARDEN_PGD (VMX, VMT, TVP, TVIP, OALLOC,KLU,KVEGTYPE,KGROUND_LAYER)  
!   ##########################################################################
!
!
!
!
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_FIX_t, ISBA_PARAM_TIME_t
USE MODD_ISBA_PGD_n, ONLY : ISBA_PGD_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!
TYPE(ISBA_PARAM_FIX_t), INTENT(INOUT) :: VMX
TYPE(ISBA_PARAM_TIME_t), INTENT(INOUT) :: VMT
TYPE(ISBA_PGD_t), INTENT(INOUT) :: TVP
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: TVIP
!
LOGICAL, INTENT(IN) :: OALLOC ! True if constant PGD fields must be allocated
INTEGER, INTENT(IN) :: KLU
INTEGER, INTENT(IN) :: KVEGTYPE
INTEGER, INTENT(IN) :: KGROUND_LAYER
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('ALLOCATE_TEB_GARDEN_PGD',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
! - Physiographic field that can evolve prognostically
!
ALLOCATE(VMT%XLAI                    (KLU                     ,1))
ALLOCATE(VMT%XVEG                    (KLU                     ,1)) 
ALLOCATE(VMT%XEMIS                   (KLU                     ,1)) 
ALLOCATE(VMT%XZ0                     (KLU                     ,1)) 
!
! - vegetation: default option (Jarvis) and general parameters:
!
ALLOCATE(VMT%XRSMIN                  (KLU                     ,1))
ALLOCATE(VMT%XGAMMA                  (KLU                     ,1)) 
ALLOCATE(VMT%XWRMAX_CF               (KLU                     ,1))
ALLOCATE(VMT%XRGL                    (KLU                     ,1))
ALLOCATE(VMT%XCV                     (KLU                     ,1)) 
!
ALLOCATE(VMT%XLAIMIN                 (KLU                     ,1)) 
ALLOCATE(VMT%XSEFOLD                 (KLU                     ,1)) 
ALLOCATE(VMT%XGMES                   (KLU                     ,1))
ALLOCATE(VMT%XGC                     (KLU                     ,1)) 
ALLOCATE(VMT%XF2I                    (KLU                     ,1)) 
ALLOCATE(VMT%XBSLAI                  (KLU                     ,1)) 
!
! - vegetation:
!
ALLOCATE(VMT%XALBNIR_VEG             (KLU                     ,1)) 
ALLOCATE(VMT%XALBVIS_VEG             (KLU                     ,1)) 
ALLOCATE(VMT%XALBUV_VEG              (KLU                     ,1)) 
!
ALLOCATE(VMT%LSTRESS                 (KLU                     ,1)) 
!
!-------------------------------------------------------------------------------
!
! - vegetation: Ags Nitrogen-model parameters ('NIT' option)
!
ALLOCATE(VMT%XCE_NITRO               (KLU                     ,1)) 
ALLOCATE(VMT%XCF_NITRO               (KLU                     ,1)) 
ALLOCATE(VMT%XCNA_NITRO              (KLU                     ,1)) 
!
IF (.NOT. OALLOC) THEN
  IF (LHOOK) CALL DR_HOOK('ALLOCATE_TEB_GARDEN_PGD',1,ZHOOK_HANDLE)
  RETURN
END IF
!-------------------------------------------------------------------------------
!
! Mask and number of grid elements containing tiles:
ALLOCATE(VMX%XVEGTYPE                (KLU,KVEGTYPE            ))
!
!-------------------------------------------------------------------------------
!
! Input Parameters:
!
! - vegetation + bare soil:
!
ALLOCATE(VMX%XZ0_O_Z0H               (KLU                     ,1)) 
!
ALLOCATE(VMX%XROOTFRAC               (KLU,KGROUND_LAYER       ,1))
ALLOCATE(VMX%NWG_LAYER               (KLU                     ,1))
ALLOCATE(VMX%XDROOT                  (KLU                     ,1))
ALLOCATE(VMX%XDG2                    (KLU                     ,1))
!
!-------------------------------------------------------------------------------
!
! - vegetation: Ags parameters ('AGS', 'LAI', 'AST', 'LST', 'NIT' options)
!
!
ALLOCATE(VMX%XH_TREE                 (KLU                     ,1)) 
!
!
ALLOCATE(VMX%XRE25                   (KLU                     ,1))
!
!-------------------------------------------------------------------------------
!
! - vegetation: Ags Stress parameters ('AST', 'LST', 'NIT' options)
!
!
ALLOCATE(TVIP%XAH                     (KLU                     ,1)) 
ALLOCATE(TVIP%XBH                     (KLU                     ,1)) 
!
ALLOCATE(VMX%XDMAX                   (KLU                     ,1)) 
!
!-------------------------------------------------------------------------------
!
! - soil: primary parameters
!
ALLOCATE(TVP%XSOC                 (KLU,KGROUND_LAYER          ))  
ALLOCATE(TVP%XSAND                   (KLU,KGROUND_LAYER       )) 
ALLOCATE(TVP%XCLAY                   (KLU,KGROUND_LAYER       )) 
ALLOCATE(TVP%XRUNOFFB                (KLU                     )) 
ALLOCATE(TVP%XWDRAIN                 (KLU                     )) 
ALLOCATE(TVIP%XTAUICE                 (KLU                    ))
!
ALLOCATE(TVIP%XGAMMAT                 (KLU                     )) 
!
ALLOCATE(VMX%XDG                     (KLU,KGROUND_LAYER       ,1)) 
!
ALLOCATE(TVIP%XRUNOFFD                (KLU                     ,1)) 
!
!-------------------------------------------------------------------------------
!
! - SGH scheme
!                                   
ALLOCATE(VMX%XD_ICE                  (KLU                     ,1)) 
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ALLOCATE_TEB_GARDEN_PGD',1,ZHOOK_HANDLE)
!
END SUBROUTINE ALLOCATE_TEB_GARDEN_PGD
