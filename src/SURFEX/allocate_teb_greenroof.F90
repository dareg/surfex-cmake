!     #########
    SUBROUTINE ALLOCATE_TEB_GREENROOF (TGR, TVG, &
                                       KLU,KLAYER_GR)
!   ##########################################################################
!
!
!
!
USE MODD_TEB_GREENROOF_n, ONLY : TEB_GREENROOF_t
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_OPTIONS_t
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!
TYPE(TEB_GREENROOF_t), INTENT(INOUT) :: TGR
TYPE(TEB_VEG_OPTIONS_t), INTENT(INOUT) :: TVG
!
INTEGER, INTENT(IN) :: KLU
INTEGER, INTENT(IN) :: KLAYER_GR
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! Mask and number of grid elements containing patches/tiles:
!
IF (LHOOK) CALL DR_HOOK('ALLOCATE_TEB_GREENROOF',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
! Averaged Surface radiative parameters:
!
ALLOCATE(TGR%XSNOWFREE_ALB           (KLU))
ALLOCATE(TGR%XSNOWFREE_ALB_VEG       (KLU))
ALLOCATE(TGR%XSNOWFREE_ALB_SOIL      (KLU))
!
!-------------------------------------------------------------------------------
!
! Prognostic variables:
!
! - Soil and vegetation heat and water:
!
ALLOCATE(TGR%XWR                     (KLU                     )) 
ALLOCATE(TGR%XTG                     (KLU,KLAYER_GR       )) 
ALLOCATE(TGR%XWG                     (KLU,KLAYER_GR       )) 
ALLOCATE(TGR%XWGI                    (KLU,KLAYER_GR       )) 
ALLOCATE(TGR%XRESA                   (KLU                     )) 
!
! - Vegetation: Ags Prognostic (YPHOTO = 'LAI', 'LST', 'AGS' or 'LST')
!
ALLOCATE(TGR%XAN                     (KLU                     )) 
ALLOCATE(TGR%XANDAY                  (KLU                     )) 
ALLOCATE(TGR%XANFM                   (KLU                     )) 
ALLOCATE(TGR%XLE                     (KLU                     ))
!
! - Vegetation (Ags 'NIT' 'NCB' option):
!
ALLOCATE(TGR%XBIOMASS                (KLU,TVG%NNBIOMASS           ))
ALLOCATE(TGR%XRESP_BIOMASS           (KLU,TVG%NNBIOMASS           ))
!
IF (LHOOK) CALL DR_HOOK('ALLOCATE_TEB_GREENROOF',1,ZHOOK_HANDLE)
!
END SUBROUTINE ALLOCATE_TEB_GREENROOF
