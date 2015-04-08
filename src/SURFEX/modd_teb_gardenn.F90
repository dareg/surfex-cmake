!##################
MODULE MODD_TEB_GARDEN_n
!##################
!
!!****  *MODD_TEB_GARDEN - declaration of packed surface parameters for ISBA scheme
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      A. Lemonsu   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2011
!!      V. Masson      06/2013 splits module in 4
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TYPE_SNOW
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

!-------------------------------------------------------------------------------
TYPE TEB_GARDEN_t
!-------------------------------------------------------------------------------
!
! Prognostic variables:
!
! - Snow Cover:
!
  TYPE(SURF_SNOW)                       :: TSNOW         ! snow state: 
!                                                      ! scheme type/option                      (-)
!                                                      ! number of layers                        (-)
!                                                      ! snow (& liq. water) content             (kg/m2)
!                                                      ! heat content                            (J/m2)
!                                                      ! temperature                             (K)
!                                                      ! density                                 (kg m-3)
!
!-------------------------------------------------------------------------------
!
! - Soil and vegetation heat and water:
!
  REAL, POINTER, DIMENSION(:)     :: XWR           ! liquid water retained on the
!                                                      ! foliage of the vegetation
!                                                      ! canopy                                  (kg/m2)
  REAL, POINTER, DIMENSION(:,:)   :: XTG           ! surface and sub-surface soil 
!                                                      ! temperature profile                     (K)
  REAL, POINTER, DIMENSION(:,:)   :: XWG           ! soil volumetric water content profile   (m3/m3)
  REAL, POINTER, DIMENSION(:,:)   :: XWGI          ! soil liquid water equivalent volumetric 
!                                                      ! ice content profile                     (m3/m3)
  REAL, POINTER, DIMENSION(:)     :: XRESA         ! aerodynamic resistance                  (s/m)
!
!-------------------------------------------------------------------------------
!
! - Vegetation: Ags Prognostic (YPHOTO = 'AGS', 'LAI', 'AST', 'LST', 'NIT', 'NCB')
!
  REAL, POINTER, DIMENSION(:)     :: XAN           ! net CO2 assimilation                    (mg/m2/s)
  REAL, POINTER, DIMENSION(:)     :: XANDAY        ! daily net CO2 assimilation              (mg/m2)
  REAL, POINTER, DIMENSION(:)     :: XANFM         ! maximum leaf assimilation               (mg/m2/s)
  REAL, POINTER, DIMENSION(:)     :: XLE           ! evapotranspiration                      (W/m2)
  REAL, POINTER, DIMENSION(:)     :: XFAPARC       ! Fapar of vegetation (cumul)
  REAL, POINTER, DIMENSION(:)     :: XFAPIRC       ! Fapir of vegetation (cumul)
  REAL, POINTER, DIMENSION(:)     :: XLAI_EFFC     ! Effective LAI (cumul)
  REAL, POINTER, DIMENSION(:)     :: XMUS          ! cos zenithal angle (cumul)     
!
!-------------------------------------------------------------------------------
!
! - Vegetation: Ags Prognostic (YPHOTO = 'NIT', 'NCB')
!
  REAL, POINTER, DIMENSION(:,:)   :: XRESP_BIOMASS    ! daily cumulated respiration of 
!                                                       ! biomass                              (kg/m2/s)
  REAL, POINTER, DIMENSION(:,:)   :: XBIOMASS         ! biomass of previous day              (kg/m2) 
!
!
!-------------------------------------------------------------------------------
!
! - Snow and flood fractions and total albedo at time t:
!
  REAL, POINTER, DIMENSION(:) :: XPSNG         ! Snow fraction over ground
  REAL, POINTER, DIMENSION(:) :: XPSNV         ! Snow fraction over vegetation
  REAL, POINTER, DIMENSION(:) :: XPSNV_A       ! Snow fraction over vegetation
  REAL, POINTER, DIMENSION(:) :: XPSN          ! Total Snow fraction
! 
  REAL, POINTER, DIMENSION(:)   :: XSNOWFREE_ALB     ! snow free albedo                        (-)
  REAL, POINTER, DIMENSION(:)   :: XSNOWFREE_ALB_VEG ! snow free albedo for vegetation         (-)
  REAL, POINTER, DIMENSION(:)   :: XSNOWFREE_ALB_SOIL! snow free albedo for soil
!-------------------------------------------------------------------------------
!                                 
END TYPE TEB_GARDEN_t


TYPE(TEB_GARDEN_t),         ALLOCATABLE, TARGET, SAVE :: TEB_GARDEN_MODEL(:,:)

TYPE(TEB_GARDEN_t), POINTER :: TEB_GARDEN => NULL()
!$OMP THREADPRIVATE(TEB_GARDEN)

CONTAINS


SUBROUTINE TEB_GARDEN_GOTO_MODEL(KFROM, KTO, LKFROM, KFROM_PATCH, KTO_PATCH)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
INTEGER, INTENT(IN) :: KFROM_PATCH, KTO_PATCH
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_TEB_GARDEN_N:TEB_GARDEN_GOTO_MODEL',0,ZHOOK_HANDLE)

TEB_GARDEN => TEB_GARDEN_MODEL(KTO,KTO_PATCH)

IF (LHOOK) CALL DR_HOOK('MODD_TEB_GARDEN_N:TEB_GARDEN_GOTO_MODEL',1,ZHOOK_HANDLE)
!
END SUBROUTINE TEB_GARDEN_GOTO_MODEL

SUBROUTINE TEB_GARDEN_ALLOC(KMODEL,KPATCH)
INTEGER, INTENT(IN) :: KMODEL, KPATCH
INTEGER :: J, JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GARDEN_N:TEB_GARDEN_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(TEB_GARDEN_MODEL(KMODEL,KPATCH))
DO J=1,KMODEL
 DO JP=1,KPATCH
  NULLIFY(TEB_GARDEN_MODEL(J,JP)%XWR)
  NULLIFY(TEB_GARDEN_MODEL(J,JP)%XTG)
  NULLIFY(TEB_GARDEN_MODEL(J,JP)%XWG)
  NULLIFY(TEB_GARDEN_MODEL(J,JP)%XWGI)
  NULLIFY(TEB_GARDEN_MODEL(J,JP)%XRESA)
  NULLIFY(TEB_GARDEN_MODEL(J,JP)%XAN)
  NULLIFY(TEB_GARDEN_MODEL(J,JP)%XANDAY)
  NULLIFY(TEB_GARDEN_MODEL(J,JP)%XANFM)
  NULLIFY(TEB_GARDEN_MODEL(J,JP)%XLE)
  NULLIFY(TEB_GARDEN_MODEL(J,JP)%XFAPARC)
  NULLIFY(TEB_GARDEN_MODEL(J,JP)%XFAPIRC)
  NULLIFY(TEB_GARDEN_MODEL(J,JP)%XLAI_EFFC)  
  NULLIFY(TEB_GARDEN_MODEL(J,JP)%XMUS)    
  NULLIFY(TEB_GARDEN_MODEL(J,JP)%XRESP_BIOMASS)
  NULLIFY(TEB_GARDEN_MODEL(J,JP)%XBIOMASS)
  NULLIFY(TEB_GARDEN_MODEL(J,JP)%XPSNG)
  NULLIFY(TEB_GARDEN_MODEL(J,JP)%XPSNV)
  NULLIFY(TEB_GARDEN_MODEL(J,JP)%XPSNV_A)
  NULLIFY(TEB_GARDEN_MODEL(J,JP)%XPSN)
  NULLIFY(TEB_GARDEN_MODEL(J,JP)%XSNOWFREE_ALB)
  NULLIFY(TEB_GARDEN_MODEL(J,JP)%XSNOWFREE_ALB_VEG)
  NULLIFY(TEB_GARDEN_MODEL(J,JP)%XSNOWFREE_ALB_SOIL)
 ENDDO
ENDDO
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GARDEN_N:TEB_GARDEN_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_GARDEN_ALLOC

SUBROUTINE TEB_GARDEN_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GARDEN_N:TEB_GARDEN_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(TEB_GARDEN_MODEL)) DEALLOCATE(TEB_GARDEN_MODEL)
IF (ASSOCIATED(TEB_GARDEN)) NULLIFY(TEB_GARDEN)
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GARDEN_N:TEB_GARDEN_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_GARDEN_DEALLO


END MODULE MODD_TEB_GARDEN_n
