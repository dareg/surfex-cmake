!     #########
    SUBROUTINE ALLOCATE_TEB_GREENROOF (TVR, KLU,KLAYER_GR,KNBIOMASS)
!   ##########################################################################
!
!
!
!
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!
TYPE(ISBA_PROG_t), INTENT(INOUT) :: TVR
!
INTEGER, INTENT(IN) :: KLU
INTEGER, INTENT(IN) :: KLAYER_GR
INTEGER, INTENT(IN) :: KNBIOMASS
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
ALLOCATE(TVR%XSNOWFREE_ALB           (KLU,1))
ALLOCATE(TVR%XSNOWFREE_ALB_VEG       (KLU,1))
ALLOCATE(TVR%XSNOWFREE_ALB_SOIL      (KLU,1))
!
!-------------------------------------------------------------------------------
!
! Prognostic variables:
!
! - Soil and vegetation heat and water:
!
ALLOCATE(TVR%XWR                     (KLU                     ,1)) 
ALLOCATE(TVR%XTG                     (KLU,KLAYER_GR       ,1)) 
ALLOCATE(TVR%XWG                     (KLU,KLAYER_GR       ,1)) 
ALLOCATE(TVR%XWGI                    (KLU,KLAYER_GR       ,1)) 
ALLOCATE(TVR%XRESA                   (KLU                     ,1)) 
!
! - Vegetation: Ags Prognostic (YPHOTO = 'LAI', 'LST', 'AGS' or 'LST')
!
ALLOCATE(TVR%XAN                     (KLU                     ,1)) 
ALLOCATE(TVR%XANDAY                  (KLU                     ,1)) 
ALLOCATE(TVR%XANFM                   (KLU                     ,1)) 
ALLOCATE(TVR%XLE                     (KLU                     ,1))
!
! - Vegetation (Ags 'NIT' 'NCB' option):
!
ALLOCATE(TVR%XBIOMASS                (KLU,KNBIOMASS           ,1))
ALLOCATE(TVR%XRESP_BIOMASS           (KLU,KNBIOMASS           ,1))
!
IF (LHOOK) CALL DR_HOOK('ALLOCATE_TEB_GREENROOF',1,ZHOOK_HANDLE)
!
END SUBROUTINE ALLOCATE_TEB_GREENROOF
