!     #########
    SUBROUTINE ALLOCATE_TEB_GREENROOF_PGD (TVM, TVP, TVIP, &
                                           OALLOC,KLU,KVEGTYPE,KLAYER_GR, KDIMTAB)  
!   ##########################################################################
!
!
!
!
USE MODD_TEB_VEG_PARAM_n, ONLY : TEB_VEG_PARAM_t
USE MODD_ISBA_PGD_n, ONLY : ISBA_PGD_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!
TYPE(TEB_VEG_PARAM_t), INTENT(INOUT) :: TVM
TYPE(ISBA_PGD_t), INTENT(INOUT) :: TVP
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: TVIP
!
LOGICAL, INTENT(IN) :: OALLOC ! True if constant PGD fields must be allocated
INTEGER, INTENT(IN) :: KLU
INTEGER, INTENT(IN) :: KVEGTYPE
INTEGER, INTENT(IN) :: KLAYER_GR
INTEGER, INTENT(IN) :: KDIMTAB
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ALLOCATE_TEB_GREENROOF_PGD',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
! - Physiographic field that can evolve prognostically
!
ALLOCATE(TVM%T%CUR%XLAI                    (KLU                     ,1))
ALLOCATE(TVM%T%CUR%XVEG                    (KLU                     ,1)) 
ALLOCATE(TVM%T%CUR%XEMIS                   (KLU                     ,1)) 
ALLOCATE(TVM%T%CUR%XZ0                     (KLU                     ,1)) 
!
ALLOCATE(TVM%T%CUR%XRSMIN                  (KLU                     ,1)) 
ALLOCATE(TVM%T%CUR%XGAMMA                  (KLU                     ,1)) 
ALLOCATE(TVM%T%CUR%XWRMAX_CF               (KLU                     ,1)) 
ALLOCATE(TVM%T%CUR%XRGL                    (KLU                     ,1)) 
ALLOCATE(TVM%T%CUR%XCV                     (KLU                     ,1)) 
!
ALLOCATE(TVM%T%CUR%XLAIMIN                 (KLU                     ,1)) 
ALLOCATE(TVM%T%CUR%XSEFOLD                 (KLU                     ,1)) 
ALLOCATE(TVM%T%CUR%XGMES                   (KLU                     ,1))
ALLOCATE(TVM%T%CUR%XGC                     (KLU                     ,1))
ALLOCATE(TVM%T%CUR%XF2I                    (KLU                     ,1)) 
ALLOCATE(TVM%T%CUR%XBSLAI                  (KLU                     ,1))
!
! - vegetation:
!
ALLOCATE(TVM%T%CUR%XALBNIR_VEG             (KLU                     ,1)) 
ALLOCATE(TVM%T%CUR%XALBVIS_VEG             (KLU                     ,1)) 
ALLOCATE(TVM%T%CUR%XALBUV_VEG              (KLU                     ,1)) 
!
!
ALLOCATE(TVM%T%CUR%LSTRESS                 (KLU                     ,1)) 
!
!-------------------------------------------------------------------------------
!
! - vegetation: Ags Nitrogen-model parameters ('NIT' option)
!
ALLOCATE(TVM%T%CUR%XCE_NITRO               (KLU                     ,1)) 
ALLOCATE(TVM%T%CUR%XCF_NITRO               (KLU                     ,1)) 
ALLOCATE(TVM%T%CUR%XCNA_NITRO              (KLU                     ,1)) 
!
IF (.NOT. OALLOC) THEN
  IF (LHOOK) CALL DR_HOOK('ALLOCATE_TEB_GARDEN_PGD',1,ZHOOK_HANDLE)
  RETURN
END IF
!-------------------------------------------------------------------------------
! Mask and number of grid elements containing tiles:
!
ALLOCATE(TVM%X%XVEGTYPE                (KLU,KVEGTYPE            ))
!
!-------------------------------------------------------------------------------
!
! Input Parameters:
!
! - vegetation + bare soil:
!
ALLOCATE(TVM%X%XZ0_O_Z0H               (KLU                     ,1)) 
!
! - vegetation: default option (Jarvis) and general parameters:
!
ALLOCATE(TVM%X%XROOTFRAC               (KLU,KLAYER_GR           ,1))
ALLOCATE(TVM%X%NWG_LAYER               (KLU                     ,1))
ALLOCATE(TVM%X%XDROOT                  (KLU                     ,1))
ALLOCATE(TVM%X%XDG2                    (KLU                     ,1))
!
!-------------------------------------------------------------------------------
!
! - vegetation: Ags parameters ('AGS', 'LAI', 'AST', 'LST', 'NIT' options)
!
ALLOCATE(TVM%X%XH_TREE                 (KLU                     ,1)) 
!
ALLOCATE(TVM%X%XRE25                   (KLU                     ,1))
!
!-------------------------------------------------------------------------------
!
! - vegetation: Ags Stress parameters ('AST', 'LST', 'NIT' options)
!
ALLOCATE(TVIP%XAH                     (KLU                     ,1)) 
ALLOCATE(TVIP%XBH                     (KLU                     ,1)) 

ALLOCATE(TVM%X%XDMAX                   (KLU                     ,1)) 
!
!-------------------------------------------------------------------------------
!
! - soil: primary parameters
!
ALLOCATE(TVIP%XTAUICE                 (KLU                     )) 

ALLOCATE(TVIP%XGAMMAT                 (KLU                     )) 

ALLOCATE(TVM%X%XDG                     (KLU,KLAYER_GR           ,1)) 

ALLOCATE(TVIP%XRUNOFFD                (KLU                     ,1)) 
!
!-------------------------------------------------------------------------------
!
! - SGH scheme
!                                   
ALLOCATE(TVM%X%XD_ICE                  (KLU                     ,1)) 
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ALLOCATE_TEB_GREENROOF_PGD',1,ZHOOK_HANDLE)
!
END SUBROUTINE ALLOCATE_TEB_GREENROOF_PGD
