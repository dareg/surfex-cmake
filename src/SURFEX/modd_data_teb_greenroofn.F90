!     ##################
      MODULE MODD_DATA_TEB_GREENROOF_n
!     ##################
!
!!****  *MODD_DATA_ISBA - declaration of DATA surface parameters for ISBA scheme
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
!!      V. Masson  *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original                    05/2005
!!      A. Lemonsu / C. de Munck     04/2011  : TEB GreenRoof
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE DATA_TEB_GREENROOF_t
!-------------------------------------------------------------------------------
!

  REAL, POINTER, DIMENSION(:,:) :: XPAR_OM_GR           ! fraction of organic matter (OM) in green roof layer
  REAL, POINTER, DIMENSION(:,:) :: XPAR_CLAY_GR         ! fraction of clay for the non-OM part of the green roof layer
  REAL, POINTER, DIMENSION(:,:) :: XPAR_SAND_GR         ! fraction of sand for the non-OM part of the green roof layer
  REAL, POINTER, DIMENSION(:,:) :: XPAR_LAI_GR          ! LAI of green roof vegetation
!
!
! Mask and number of grid elements containing patches/tiles:
!
  REAL, POINTER, DIMENSION(:,:)  :: XPAR_VEGTYPE        ! fraction of each vegetation type for
                                                        ! each grid mesh                          (-)
!
!-------------------------------------------------------------------------------
!
! Input Parameters, per patch:
!
! - vegetation + bare soil:
!
  REAL, POINTER, DIMENSION(:)   :: XPAR_Z0_O_Z0H        ! ratio of surface roughness lengths
                                                        ! (momentum to heat)                      (-)
  REAL, POINTER, DIMENSION(:,:)   :: XPAR_EMIS          ! surface emissivity                      (-)
  REAL, POINTER, DIMENSION(:,:)   :: XPAR_Z0            ! surface roughness length                (m)
!
! - vegetation:
!
  REAL, POINTER, DIMENSION(:)   :: XPAR_ALBNIR_VEG       ! vegetation near-infra-red albedo        (-)
  REAL, POINTER, DIMENSION(:)   :: XPAR_ALBVIS_VEG       ! vegetation visible albedo               (-)
  REAL, POINTER, DIMENSION(:)   :: XPAR_ALBUV_VEG        ! vegetation UV albedo                    (-)
!
! - vegetation: default option (Jarvis) and general parameters:
!
  REAL, POINTER, DIMENSION(:,:) :: XPAR_VEG              ! vegetation cover fraction               (-)
  REAL, POINTER, DIMENSION(:)   :: XPAR_WRMAX_CF         ! coefficient for maximum water 
                                                         ! interception 
                                                         ! storage capacity on the vegetation      (-)
  REAL, POINTER, DIMENSION(:)   :: XPAR_RSMIN            ! minimum stomatal resistance             (s/m)
  REAL, POINTER, DIMENSION(:)   :: XPAR_GAMMA            ! coefficient for the calculation
                                                         ! of the surface stomatal
                                                         ! resistance
  REAL, POINTER, DIMENSION(:)   :: XPAR_CV               ! vegetation thermal inertia coefficient  (K m2/J)
  REAL, POINTER, DIMENSION(:)   :: XPAR_RGL              ! maximum solar radiation
                                                         ! usable in photosynthesis                (W/m2)
  REAL, POINTER, DIMENSION(:,:)   :: XPAR_ROOTFRAC       ! root fraction profile ('DIF' option)
!
!-------------------------------------------------------------------------------
!
! - vegetation: Ags parameters ('AGS', 'LAI', 'AST', 'LST', 'NIT', 'NCB' options)
!
  REAL, POINTER, DIMENSION(:)      :: XPAR_BSLAI         ! ratio d(biomass)/d(lai)                 (kg/m2)
  REAL, POINTER, DIMENSION(:)      :: XPAR_LAIMIN        ! minimum LAI (Leaf Area Index)           (m2/m2)
  REAL, POINTER, DIMENSION(:)      :: XPAR_SEFOLD        ! e-folding time for senescence           (s)
  REAL, POINTER, DIMENSION(:)      :: XPAR_H_TREE        ! height of trees                         (m)
  REAL, POINTER, DIMENSION(:)      :: XPAR_GMES          ! mesophyll conductance                   (m s-1)
  REAL, POINTER, DIMENSION(:)      :: XPAR_RE25          ! Ecosystem respiration parameter         (kg m2 s-1)
!
!-------------------------------------------------------------------------------
!
! - vegetation: Ags Stress parameters ('AST', 'LST', 'NIT', 'NCB' options)
!
  LOGICAL, POINTER, DIMENSION(:)   :: LDATA_STRESS       ! vegetation response type to water
                                                         ! stress (true:defensive false:offensive) (-)
  REAL, POINTER, DIMENSION(:)      :: XPAR_F2I           ! critical normilized soil water 
                                                         ! content for stress parameterisation
  REAL, POINTER, DIMENSION(:)      :: XPAR_GC            ! cuticular conductance                   (m s-1)
  REAL, POINTER, DIMENSION(:)      :: XPAR_DMAX          ! maximum air saturation deficit
                                                         ! tolerate by vegetation                  (kg/kg)
!
  REAL, POINTER, DIMENSION(:)      :: XPAR_BSLAI_ST      ! ratio d(biomass)/d(lai)                 (kg/m2)
  REAL, POINTER, DIMENSION(:)      :: XPAR_SEFOLD_ST     ! e-folding time for senescence           (s)
  REAL, POINTER, DIMENSION(:)      :: XPAR_GMES_ST       ! mesophyll conductance                   (m s-1)
  REAL, POINTER, DIMENSION(:)      :: XPAR_GC_ST         ! cuticular conductance                   (m s-1)
  REAL, POINTER, DIMENSION(:)      :: XPAR_DMAX_ST       ! maximum air saturation deficit
!-------------------------------------------------------------------------------
!
! - vegetation: Ags Nitrogen-model parameters ('NIT', 'NCB' option)
!
  REAL, POINTER, DIMENSION(:)      :: XPAR_CE_NITRO      ! leaf aera ratio sensitivity to 
                                                         ! nitrogen concentration                (m2/kg)
  REAL, POINTER, DIMENSION(:)      :: XPAR_CF_NITRO      ! lethal minimum value of leaf area
                                                         ! ratio                                 (m2/kg)
  REAL, POINTER, DIMENSION(:)      :: XPAR_CNA_NITRO     ! nitrogen concentration of active 
                                                         ! biomass                               (kg/kg)
!
!-------------------------------------------------------------------------------
!
! - soil: primary parameters
!
  REAL, POINTER, DIMENSION(:,:)    :: XPAR_DG            ! soil layer thicknesses                  (m)
                                                         ! NOTE: in Force-Restore mode, the 
                                                         ! uppermost layer thickness is superficial
                                                         ! and is only explicitly used for soil 
                                                         ! water phase changes                     (m)
!
  REAL, POINTER,DIMENSION(:)       :: XPAR_DICE          ! depth of the soil column for the calculation
                                                         ! of the frozen soil fraction (m)
!
! - bare soil albedo
!
  REAL, POINTER, DIMENSION(:)   :: XPAR_ALBNIR_SOIL      ! soil near-infra-red albedo        (-)
  REAL, POINTER, DIMENSION(:)   :: XPAR_ALBVIS_SOIL      ! soil visible albedo               (-)
  REAL, POINTER, DIMENSION(:)   :: XPAR_ALBUV_SOIL       ! soil UV albedo                    (-)
  REAL, POINTER, DIMENSION(:)   :: XPAR_ALBNIR_DRY       ! dry soil near-infra-red albedo    (-)
  REAL, POINTER, DIMENSION(:)   :: XPAR_ALBVIS_DRY       ! dry soil visible albedo           (-)
  REAL, POINTER, DIMENSION(:)   :: XPAR_ALBUV_DRY        ! dry soil UV albedo                (-)
  REAL, POINTER, DIMENSION(:)   :: XPAR_ALBNIR_WET       ! wet soil near-infra-red albedo    (-)
  REAL, POINTER, DIMENSION(:)   :: XPAR_ALBVIS_WET       ! wet soil visible albedo           (-)
  REAL, POINTER, DIMENSION(:)   :: XPAR_ALBUV_WET        ! wet soil UV albedo                (-)
!
!-------------------------------------------------------------------------------
!
!- Vegetation: Ags Prognostic (YPHOTO = ('LAI', 'LST', 'NIT', or 'NCB') or prescribed (YPHOTO='NON', 'AGS' or 'AST')
!
  REAL, POINTER, DIMENSION(:,:)     :: XPAR_LAI          ! Leaf Area Index                         (m2/m2)
!
!-------------------------------------------------------------------------------
!

END TYPE DATA_TEB_GREENROOF_t

TYPE(DATA_TEB_GREENROOF_t), ALLOCATABLE, TARGET, SAVE :: DATA_TEB_GREENROOF_MODEL(:)

TYPE(DATA_TEB_GREENROOF_t), POINTER :: DATA_TEB_GREENROOF => NULL()
!$OMP THREADPRIVATE(DATA_TEB_GREENROOF)

CONTAINS

SUBROUTINE DATA_TEB_GREENROOF_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_DATA_TEB_GREENROOF_N:DATA_TEB_GREENROOF_GOTO_MODEL',0,ZHOOK_HANDLE)

DATA_TEB_GREENROOF => DATA_TEB_GREENROOF_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_DATA_TEB_GREENROOF_N:DATA_TEB_GREENROOF_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE DATA_TEB_GREENROOF_GOTO_MODEL

SUBROUTINE DATA_TEB_GREENROOF_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DATA_TEB_GREENROOF_N:DATA_TEB_GREENROOF_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(DATA_TEB_GREENROOF_MODEL(KMODEL))
DO J=1,KMODEL
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_OM_GR)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_CLAY_GR)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_SAND_GR)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_LAI_GR)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_VEGTYPE)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_Z0_O_Z0H)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_EMIS)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_Z0)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_ALBNIR_VEG)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_ALBVIS_VEG)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_ALBUV_VEG)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_VEG)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_WRMAX_CF)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_RSMIN)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_GAMMA)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_CV)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_RGL)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_ROOTFRAC)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_BSLAI)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_LAIMIN)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_SEFOLD)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_H_TREE)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_GMES)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_RE25)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%LDATA_STRESS)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_F2I)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_GC)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_DMAX)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_BSLAI_ST)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_SEFOLD_ST)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_GMES_ST)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_GC_ST)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_DMAX_ST)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_CE_NITRO)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_CF_NITRO)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_CNA_NITRO)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_DG)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_DICE)  
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_ALBNIR_SOIL)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_ALBVIS_SOIL)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_ALBUV_SOIL)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_ALBNIR_DRY)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_ALBVIS_DRY)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_ALBUV_DRY)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_ALBNIR_WET)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_ALBVIS_WET)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_ALBUV_WET)
  NULLIFY(DATA_TEB_GREENROOF_MODEL(J)%XPAR_LAI)
ENDDO
IF (LHOOK) CALL DR_HOOK("MODD_DATA_TEB_GREENROOF_N:DATA_TEB_GREENROOF_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE DATA_TEB_GREENROOF_ALLOC

SUBROUTINE DATA_TEB_GREENROOF_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DATA_TEB_GREENROOF_N:DATA_TEB_GREENROOF_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(DATA_TEB_GREENROOF_MODEL)) DEALLOCATE(DATA_TEB_GREENROOF_MODEL)
IF (ASSOCIATED(DATA_TEB_GREENROOF)) NULLIFY(DATA_TEB_GREENROOF)
IF (LHOOK) CALL DR_HOOK("MODD_DATA_TEB_GREENROOF_N:DATA_TEB_GREENROOF_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE DATA_TEB_GREENROOF_DEALLO

END MODULE MODD_DATA_TEB_GREENROOF_n
