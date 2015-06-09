!##################
MODULE MODD_TEB_GREENROOF_n
!##################
!
!!****  *MODD_TEB_GREENROOF - declaration of ISBA scheme packed surface parameters for urban green roofs
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
!!      A. Lemonsu *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       09/2009
!!      C. de Munck     06/2011 
!!      V. Masson       06/2013  splits module in 4
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

TYPE TEB_GREENROOF_1P_t
!-------------------------------------------------------------------------------
!
! Prognostic variables:
!
! - Snow Cover:
!
  TYPE(SURF_SNOW)                :: TSNOW            ! snow state: 
                                                     ! scheme type/option                      (-)
                                                     ! number of layers                        (-)
                                                     ! snow (& liq. water) content             (kg/m2)
                                                     ! heat content                            (J/m2)
                                                     ! temperature                             (K)
                                                     ! density                                 (kg m-3)
!
! - Soil and vegetation heat and water:
!
  REAL, POINTER, DIMENSION(:)    :: XWR              ! liquid water retained on the
                                                     ! foliage of the vegetation
                                                     ! canopy                                  (kg/m2)
  REAL, POINTER, DIMENSION(:,:)  :: XTG              ! surface and sub-surface soil 
                                                     ! temperature profile                     (K)
  REAL, POINTER, DIMENSION(:,:)  :: XWG              ! soil volumetric water content profile   (m3/m3)
  REAL, POINTER, DIMENSION(:,:)  :: XWGI             ! soil liquid water equivalent volumetric 
                                                     ! ice content profile                     (m3/m3)
  REAL, POINTER, DIMENSION(:)    :: XRESA            ! aerodynamic resistance                  (s/m)
!

! - Vegetation: Ags Prognostic (YPHOTO = 'AGS', 'LAI', 'AST', 'LST', 'NIT', 'NCB')
!
  REAL, POINTER, DIMENSION(:)    :: XAN              ! net CO2 assimilation                    (mg/m2/s)
  REAL, POINTER, DIMENSION(:)    :: XANDAY           ! daily net CO2 assimilation              (mg/m2)
  REAL, POINTER, DIMENSION(:)    :: XANFM            ! maximum leaf assimilation               (mg/m2/s)
  REAL, POINTER, DIMENSION(:)    :: XLE              ! evapotranspiration                      (W/m2)
  REAL, POINTER, DIMENSION(:)     :: XFAPARC       ! Fapar of vegetation (cumul)
  REAL, POINTER, DIMENSION(:)     :: XFAPIRC       ! Fapir of vegetation (cumul)
  REAL, POINTER, DIMENSION(:)     :: XLAI_EFFC     ! Effective LAI (cumul)
  REAL, POINTER, DIMENSION(:)     :: XMUS          ! cos zenithal angle (cumul)    
!
! - Vegetation: Ags Prognostic (YPHOTO = 'NIT', 'NCB')
!
  REAL, POINTER, DIMENSION(:,:)  :: XRESP_BIOMASS    ! daily cumulated respiration of 
                                                     ! biomass                              (kg/m2/s)
  REAL, POINTER, DIMENSION(:,:)  :: XBIOMASS         ! biomass of previous day              (kg/m2) 
!
! - SGH scheme
!                                                     
  REAL, POINTER, DIMENSION(:)    :: XKSAT_ICE       ! hydraulic conductivity at saturation
                                                    ! over frozen area (m s-1)                                     
!-------------------------------------------------------------------------------
!
! - Snow and flood fractions and total albedo at time t:
!
  REAL, POINTER, DIMENSION(:)    :: XPSNG              ! Snow fraction over ground
  REAL, POINTER, DIMENSION(:)    :: XPSNV              ! Snow fraction over vegetation
  REAL, POINTER, DIMENSION(:)    :: XPSNV_A            ! Snow fraction over vegetation
  REAL, POINTER, DIMENSION(:)    :: XPSN               ! Total Snow fraction
! 
  REAL, POINTER, DIMENSION(:)    :: XSNOWFREE_ALB      ! snow free albedo                        (-)
  REAL, POINTER, DIMENSION(:)    :: XSNOWFREE_ALB_VEG  ! snow free albedo for vegetation         (-)
  REAL, POINTER, DIMENSION(:)    :: XSNOWFREE_ALB_SOIL ! snow free albedo for soil               (-)
!
!-------------------------------------------------------------------------------
!
END TYPE TEB_GREENROOF_1P_t
!
TYPE TEB_GREENROOF_t
  !
  TYPE(TEB_GREENROOF_1P_t), POINTER :: ALP(:) => NULL()
  TYPE(TEB_GREENROOF_1P_t), POINTER :: CUR => NULL()
  !
END TYPE TEB_GREENROOF_t
!
TYPE(TEB_GREENROOF_t),         ALLOCATABLE, TARGET, SAVE :: TEB_GREENROOF_MODEL(:)

TYPE(TEB_GREENROOF_t), POINTER :: TEB_GREENROOF => NULL()
!$OMP THREADPRIVATE(TEB_GREENROOF)

CONTAINS

SUBROUTINE TEB_GREENROOF_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_TEB_GREENROOF_N:TEB_GREENROOF_GOTO_MODEL',0,ZHOOK_HANDLE)

TEB_GREENROOF => TEB_GREENROOF_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_TEB_GREENROOF_N:TEB_GREENROOF_GOTO_MODEL',1,ZHOOK_HANDLE)
!
END SUBROUTINE TEB_GREENROOF_GOTO_MODEL
!
SUBROUTINE TEB_GREENROOF_GOTO_PATCH(KTO_PATCH)
INTEGER, INTENT(IN) :: KTO_PATCH
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current patch is set to patch KTO_PATCH
IF (LHOOK) CALL DR_HOOK('MODD_TEB_GREENROOF_N:TEB_GREENROOF_GOTO_PATCH',0,ZHOOK_HANDLE)

TEB_GREENROOF%CUR => TEB_GREENROOF%ALP(KTO_PATCH)

IF (LHOOK) CALL DR_HOOK('MODD_TEB_GREENROOF_N:TEB_GREENROOF_GOTO_PATCH',1,ZHOOK_HANDLE)
!
END SUBROUTINE TEB_GREENROOF_GOTO_PATCH
!
SUBROUTINE TEB_GREENROOF_ALLOC(KMODEL,KPATCH)
INTEGER, INTENT(IN) :: KMODEL, KPATCH
INTEGER :: J, JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GREENROOF_N:TEB_GREENROOF_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(TEB_GREENROOF_MODEL(KMODEL))
TEB_GREENROOF => TEB_GREENROOF_MODEL(KMODEL)
DO J=1,KMODEL
 ALLOCATE(TEB_GREENROOF_MODEL(J)%ALP(KPATCH))
 DO JP=1,KPATCH
  NULLIFY(TEB_GREENROOF_MODEL(J)%ALP(JP)%XWR)
  NULLIFY(TEB_GREENROOF_MODEL(J)%ALP(JP)%XTG)
  NULLIFY(TEB_GREENROOF_MODEL(J)%ALP(JP)%XWG)
  NULLIFY(TEB_GREENROOF_MODEL(J)%ALP(JP)%XWGI)
  NULLIFY(TEB_GREENROOF_MODEL(J)%ALP(JP)%XRESA)
  NULLIFY(TEB_GREENROOF_MODEL(J)%ALP(JP)%XAN)
  NULLIFY(TEB_GREENROOF_MODEL(J)%ALP(JP)%XANDAY)
  NULLIFY(TEB_GREENROOF_MODEL(J)%ALP(JP)%XANFM)
  NULLIFY(TEB_GREENROOF_MODEL(J)%ALP(JP)%XLE)
  NULLIFY(TEB_GREENROOF_MODEL(J)%ALP(JP)%XFAPARC)
  NULLIFY(TEB_GREENROOF_MODEL(J)%ALP(JP)%XFAPIRC)
  NULLIFY(TEB_GREENROOF_MODEL(J)%ALP(JP)%XLAI_EFFC)  
  NULLIFY(TEB_GREENROOF_MODEL(J)%ALP(JP)%XMUS)     
  NULLIFY(TEB_GREENROOF_MODEL(J)%ALP(JP)%XRESP_BIOMASS)
  NULLIFY(TEB_GREENROOF_MODEL(J)%ALP(JP)%XBIOMASS)
  NULLIFY(TEB_GREENROOF_MODEL(J)%ALP(JP)%XKSAT_ICE)
  NULLIFY(TEB_GREENROOF_MODEL(J)%ALP(JP)%XPSNG)
  NULLIFY(TEB_GREENROOF_MODEL(J)%ALP(JP)%XPSNV)
  NULLIFY(TEB_GREENROOF_MODEL(J)%ALP(JP)%XPSNV_A)
  NULLIFY(TEB_GREENROOF_MODEL(J)%ALP(JP)%XPSN)
  NULLIFY(TEB_GREENROOF_MODEL(J)%ALP(JP)%XSNOWFREE_ALB)
  NULLIFY(TEB_GREENROOF_MODEL(J)%ALP(JP)%XSNOWFREE_ALB_VEG)
  NULLIFY(TEB_GREENROOF_MODEL(J)%ALP(JP)%XSNOWFREE_ALB_SOIL)
 ENDDO
ENDDO
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GREENROOF_N:TEB_GREENROOF_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_GREENROOF_ALLOC

SUBROUTINE TEB_GREENROOF_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GREENROOF_N:TEB_GREENROOF_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(TEB_GREENROOF_MODEL)) DEALLOCATE(TEB_GREENROOF_MODEL)
IF (ASSOCIATED(TEB_GREENROOF)) NULLIFY(TEB_GREENROOF)
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GREENROOF_N:TEB_GREENROOF_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_GREENROOF_DEALLO

END MODULE MODD_TEB_GREENROOF_n
