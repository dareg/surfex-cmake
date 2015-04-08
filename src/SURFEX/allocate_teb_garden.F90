!     #########
    SUBROUTINE ALLOCATE_TEB_GARDEN(KLU,KGROUND_LAYER)  
!   ##########################################################################
!
USE MODD_TEB_VEG_n, ONLY : TVG => TEB_VEG_OPTIONS
USE MODD_TEB_GARDEN_n, ONLY : TGD => TEB_GARDEN
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: KLU
INTEGER, INTENT(IN) :: KGROUND_LAYER
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! Mask and number of grid elements containing patches/tiles:
!
IF (LHOOK) CALL DR_HOOK('ALLOCATE_TEB_GARDEN',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
! Averaged Surface radiative parameters:
!
ALLOCATE(TGD%XSNOWFREE_ALB           (KLU))
ALLOCATE(TGD%XSNOWFREE_ALB_VEG       (KLU))
ALLOCATE(TGD%XSNOWFREE_ALB_SOIL      (KLU))
!
!-------------------------------------------------------------------------------
!
! Prognostic variables:
!
!
! - Soil and vegetation heat and water:
!
ALLOCATE(TGD%XWR                     (KLU                    )) 
ALLOCATE(TGD%XTG                     (KLU,KGROUND_LAYER      )) 
ALLOCATE(TGD%XWG                     (KLU,KGROUND_LAYER      )) 
ALLOCATE(TGD%XWGI                    (KLU,KGROUND_LAYER      )) 
ALLOCATE(TGD%XRESA                   (KLU                    )) 
!
! - Vegetation: Ags Prognostic (YPHOTO = 'LAI', 'LST', 'AGS' or 'LST')
!
ALLOCATE(TGD%XAN                     (KLU                    )) 
ALLOCATE(TGD%XANDAY                  (KLU                    )) 
ALLOCATE(TGD%XANFM                   (KLU                    )) 
ALLOCATE(TGD%XLE                     (KLU                    ))
!
! - Vegetation (Ags 'NIT' 'NCB' option):
!
ALLOCATE(TGD%XBIOMASS                (KLU,TVG%NNBIOMASS          ))
ALLOCATE(TGD%XRESP_BIOMASS           (KLU,TVG%NNBIOMASS          ))
!
IF (LHOOK) CALL DR_HOOK('ALLOCATE_TEB_GARDEN',1,ZHOOK_HANDLE)
!
END SUBROUTINE ALLOCATE_TEB_GARDEN
