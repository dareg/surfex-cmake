!     ##################
      MODULE MODD_DATA_ISBA_n
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
!!      Original       05/2005
!!      P Samuelsson   02/2012  MEB
!
!*       0.   DECLARATIONS
!             ------------
!
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE DATA_ISBA_t
!-------------------------------------------------------------------------------
!
! Mask and number of grid elements containing patches/tiles:
!
  REAL, POINTER, DIMENSION(:,:)  :: XPAR_VEGTYPE       ! fraction of each vegetation type for
!                                                  ! each grid mesh                          (-)
!
!-------------------------------------------------------------------------------
!
  INTEGER                       :: NTIME               ! number of time data
!                                                      ! for VEG, LAI, EMIS, Z0
  LOGICAL                        :: LDATA_MIXPAR
!
  LOGICAL                        :: LDATA_VEGTYPE
  LOGICAL                        :: LDATA_LAI
  LOGICAL                        :: LDATA_H_TREE
  LOGICAL                        :: LDATA_DG
  LOGICAL                        :: LDATA_DICE  
  LOGICAL                        :: LDATA_ROOTFRAC
  LOGICAL                        :: LDATA_GROUND_DEPTH
  LOGICAL                        :: LDATA_ROOT_DEPTH
  LOGICAL                        :: LDATA_ROOT_EXTINCTION
  LOGICAL                        :: LDATA_ROOT_LIN
  LOGICAL                        :: LDATA_VEG
  LOGICAL                        :: LDATA_Z0
  LOGICAL                        :: LDATA_EMIS
  LOGICAL                        :: LDATA_ALBNIR_VEG
  LOGICAL                        :: LDATA_ALBVIS_VEG
  LOGICAL                        :: LDATA_ALBUV_VEG
  LOGICAL                        :: LDATA_RSMIN
  LOGICAL                        :: LDATA_GAMMA
  LOGICAL                        :: LDATA_WRMAX_CF
  LOGICAL                        :: LDATA_CV
  LOGICAL                        :: LDATA_Z0_O_Z0H
  LOGICAL                        :: LDATA_RGL
  LOGICAL                        :: LDATA_BSLAI
  LOGICAL                        :: LDATA_LAIMIN
  LOGICAL                        :: LDATA_SEFOLD
  LOGICAL                        :: LDATA_GMES
  LOGICAL                        :: LDATA_RE25
  LOGICAL                        :: LDATA_STRESS
  LOGICAL                        :: LDATA_F2I
  LOGICAL                        :: LDATA_GC
  LOGICAL                        :: LDATA_DMAX
  LOGICAL                        :: LDATA_CE_NITRO
  LOGICAL                        :: LDATA_CF_NITRO
  LOGICAL                        :: LDATA_CNA_NITRO
  LOGICAL                        :: LDATA_ALBNIR_SOIL
  LOGICAL                        :: LDATA_ALBVIS_SOIL
  LOGICAL                        :: LDATA_ALBUV_SOIL
  LOGICAL                        :: LDATA_IRRIG
  LOGICAL                        :: LDATA_WATSUP
! - For multi-energy balance (MEB)
!
  LOGICAL                        :: LDATA_GNDLITTER
  LOGICAL                        :: LDATA_LAIGV
  LOGICAL                        :: LDATA_Z0LITTER
  LOGICAL                        :: LDATA_RSMINGV
  LOGICAL                        :: LDATA_GAMMAGV
  LOGICAL                        :: LDATA_WRMAX_CFGV
  LOGICAL                        :: LDATA_RGLGV
  LOGICAL                        :: LDATA_ZF_TALLVEG
  LOGICAL                        :: LDATA_ROOTFRACGV
  LOGICAL                        :: LDATA_ROOT_DEPTHGV
  LOGICAL                        :: LDATA_ROOT_EXTINCTIONGV
  LOGICAL                        :: LDATA_H_VEG
!
!  
! Input Parameters, per patch:
!
! - vegetation + bare soil:
!
  REAL, POINTER, DIMENSION(:,:) :: XPAR_Z0_O_Z0H         ! ratio of surface roughness lengths
!                                                      ! (momentum to heat)                      (-)
  REAL, POINTER, DIMENSION(:,:,:) :: XPAR_EMIS             ! surface emissivity                      (-)
  REAL, POINTER, DIMENSION(:,:,:) :: XPAR_Z0               ! surface roughness length                (m)
!
! - vegetation:
!
  REAL, POINTER, DIMENSION(:,:) :: XPAR_ALBNIR_VEG       ! vegetation near-infra-red albedo        (-)
  REAL, POINTER, DIMENSION(:,:) :: XPAR_ALBVIS_VEG       ! vegetation visible albedo               (-)
  REAL, POINTER, DIMENSION(:,:) :: XPAR_ALBUV_VEG        ! vegetation UV albedo                    (-)
!
! - vegetation: default option (Jarvis) and general parameters:
!
  REAL, POINTER, DIMENSION(:,:,:) :: XPAR_VEG            ! vegetation cover fraction               (-)
  REAL, POINTER, DIMENSION(:,:) :: XPAR_WRMAX_CF         ! coefficient for maximum water 
!                                                      ! interception 
!                                                      ! storage capacity on the vegetation      (-)
  REAL, POINTER, DIMENSION(:,:) :: XPAR_RSMIN            ! minimum stomatal resistance             (s/m)
  REAL, POINTER, DIMENSION(:,:) :: XPAR_GAMMA            ! coefficient for the calculation
!                                                      ! of the surface stomatal
!                                                      ! resistance
  REAL, POINTER, DIMENSION(:,:) :: XPAR_CV               ! vegetation thermal inertia coefficient  (K m2/J)
  REAL, POINTER, DIMENSION(:,:) :: XPAR_RGL              ! maximum solar radiation
!                                                      ! usable in photosynthesis                (W/m2)
  REAL, POINTER, DIMENSION(:,:,:) :: XPAR_ROOTFRAC       ! root fraction profile ('DIF' option)
!
  REAL, POINTER, DIMENSION(:,:) :: XPAR_ROOT_DEPTH       ! root depth ('DIF' option)
  REAL, POINTER, DIMENSION(:,:) :: XPAR_ROOT_EXTINCTION  ! root extinction parameter ('DIF' option)
  REAL, POINTER, DIMENSION(:,:) :: XPAR_ROOT_LIN         ! root linear parameter ('DIF' option)
!
! - For multi-energy balance (MEB)
!
  REAL, POINTER, DIMENSION(:,:,:) :: XPAR_GNDLITTER      ! ground litter fraction
  REAL, POINTER, DIMENSION(:,:,:) :: XPAR_LAIGV          ! understory LAI
  REAL, POINTER, DIMENSION(:,:,:) :: XPAR_Z0LITTER       ! ground litter roughness length
  REAL, POINTER, DIMENSION(:,:)   :: XPAR_RSMINGV        ! understory minimum surface resistance
  REAL, POINTER, DIMENSION(:,:)   :: XPAR_GAMMAGV        !
  REAL, POINTER, DIMENSION(:,:)   :: XPAR_WRMAX_CFGV     !
  REAL, POINTER, DIMENSION(:,:)   :: XPAR_RGLGV          !
  REAL, POINTER, DIMENSION(:,:)   :: XPAR_ZF_TALLVEG     !
  REAL, POINTER, DIMENSION(:,:,:) :: XPAR_ROOTFRACGV     ! understory root fraction profile
  REAL, POINTER, DIMENSION(:,:,:) :: XPAR_H_VEG          ! height of canopy vegetation
  REAL, POINTER, DIMENSION(:,:) :: XPAR_ROOT_DEPTHGV       ! root depth ('DIF' option)
  REAL, POINTER, DIMENSION(:,:) :: XPAR_ROOT_EXTINCTIONGV  ! root extinction parameter ('DIF' option)
!
!  
!-------------------------------------------------------------------------------
!
! - vegetation: Ags parameters ('AGS', 'LAI', 'AST', 'LST', 'NIT', 'NCB' options)
!
  REAL, POINTER, DIMENSION(:,:)    :: XPAR_BSLAI        ! ratio d(biomass)/d(lai)                 (kg/m2)
  REAL, POINTER, DIMENSION(:,:)    :: XPAR_LAIMIN       ! minimum LAI (Leaf Area Index)           (m2/m2)
  REAL, POINTER, DIMENSION(:,:)    :: XPAR_SEFOLD       ! e-folding time for senescence           (s)
  REAL, POINTER, DIMENSION(:,:)    :: XPAR_H_TREE       ! height of trees                         (m)
  REAL, POINTER, DIMENSION(:,:)    :: XPAR_GMES         ! mesophyll conductance                   (m s-1)
  REAL, POINTER, DIMENSION(:,:)    :: XPAR_RE25         ! Ecosystem respiration parameter         (kg m2 s-1)
!
!-------------------------------------------------------------------------------
!
! - vegetation: Ags Stress parameters ('AST', 'LST', 'NIT', 'NCB' options)
!
  LOGICAL, POINTER, DIMENSION(:,:) :: LPAR_STRESS       ! vegetation response type to water
!                                                     ! stress (true:defensive false:offensive) (-)
  REAL, POINTER, DIMENSION(:,:)    :: XPAR_F2I          ! critical normilized soil water 
!                                                     ! content for stress parameterisation
  REAL, POINTER, DIMENSION(:,:)    :: XPAR_GC           ! cuticular conductance                   (m s-1)
  REAL, POINTER, DIMENSION(:,:)    :: XPAR_DMAX         ! maximum air saturation deficit
!                                                     ! tolerate by vegetation                  (kg/kg)
!
!-------------------------------------------------------------------------------
!
! - vegetation: Ags Nitrogen-model parameters ('NIT', 'NCB' option)
!
  REAL, POINTER, DIMENSION(:,:)    :: XPAR_CE_NITRO       ! leaf aera ratio sensitivity to 
!                                                       ! nitrogen concentration                (m2/kg)
  REAL, POINTER, DIMENSION(:,:)    :: XPAR_CF_NITRO       ! lethal minimum value of leaf area
!                                                       ! ratio                                 (m2/kg)
  REAL, POINTER, DIMENSION(:,:)    :: XPAR_CNA_NITRO      ! nitrogen concentration of active 
!                                                       ! biomass                               (kg/kg)
!
!-------------------------------------------------------------------------------
!
! - soil: primary parameters
!
  REAL, POINTER, DIMENSION(:,:,:)  :: XPAR_DG          ! soil layer depth                        (m)
!                                                      ! NOTE: in Force-Restore mode, the 
!                                                      ! uppermost layer thickness is superficial
!                                                      ! and is only explicitly used for soil 
!                                                      ! water phase changes                     (m)
!
  REAL, POINTER,DIMENSION(:,:)     :: XPAR_GROUND_DEPTH ! ground depth (DIF option)
!
  REAL, POINTER,DIMENSION(:,:)     :: XPAR_DICE        ! depth of the soil column for the calculation
!                                                        of the frozen soil fraction (m) (Force restore)
!
! - bare soil albedo
!
  REAL, POINTER, DIMENSION(:,:) :: XPAR_ALBNIR_SOIL      ! soil near-infra-red albedo        (-)
  REAL, POINTER, DIMENSION(:,:) :: XPAR_ALBVIS_SOIL      ! soil visible albedo               (-)
  REAL, POINTER, DIMENSION(:,:) :: XPAR_ALBUV_SOIL       ! soil UV albedo                    (-)
!
!-------------------------------------------------------------------------------
!
! - Vegetation: Ags Prognostic (YPHOTO = ('LAI', 'LST', 'NIT', or 'NCB') or prescribed (YPHOTO='NON', 'AGS' or 'AST')
!
  REAL, POINTER, DIMENSION(:,:,:)     :: XPAR_LAI          ! Leaf Area Index                         (m2/m2)
!
  REAL, POINTER, DIMENSION(:,:,:)     :: XPAR_IRRIG
  REAL, POINTER, DIMENSION(:,:,:)     :: XPAR_WATSUP
!
!-------------------------------------------------------------------------------
!

END TYPE DATA_ISBA_t

TYPE(DATA_ISBA_t), ALLOCATABLE, TARGET, SAVE :: DATA_ISBA_MODEL(:)

TYPE(DATA_ISBA_t), POINTER :: DATA_ISBA => NULL()
!$OMP THREADPRIVATE(DATA_ISBA)

CONTAINS

SUBROUTINE DATA_ISBA_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_DATA_ISBA_N:DATA_ISBA_GOTO_MODEL',0,ZHOOK_HANDLE)

DATA_ISBA => DATA_ISBA_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_DATA_ISBA_N:DATA_ISBA_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE DATA_ISBA_GOTO_MODEL

SUBROUTINE DATA_ISBA_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DATA_ISBA_N:DATA_ISBA_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(DATA_ISBA_MODEL(KMODEL))
DATA_ISBA => DATA_ISBA_MODEL(KMODEL)
DO J=1,KMODEL
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_VEGTYPE)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_Z0_O_Z0H)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_EMIS)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_Z0)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_ALBNIR_VEG)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_ALBVIS_VEG)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_ALBUV_VEG)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_VEG)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_WRMAX_CF)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_RSMIN)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_GAMMA)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_CV)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_RGL)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_ROOTFRAC)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_BSLAI)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_LAIMIN)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_SEFOLD)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_H_TREE)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_GMES)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_RE25)
  NULLIFY(DATA_ISBA_MODEL(J)%LPAR_STRESS)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_F2I)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_GC)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_DMAX)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_CE_NITRO)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_CF_NITRO)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_CNA_NITRO)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_DG)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_DICE)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_GROUND_DEPTH)  
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_ROOT_DEPTH)  
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_ROOT_EXTINCTION)  
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_ROOT_LIN)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_ALBNIR_SOIL)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_ALBVIS_SOIL)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_ALBUV_SOIL)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_LAI)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_IRRIG)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_WATSUP)  
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_GNDLITTER)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_LAIGV)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_Z0LITTER)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_RSMINGV)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_GAMMAGV)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_WRMAX_CFGV)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_RGLGV)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_ZF_TALLVEG)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_ROOTFRACGV)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_ROOT_DEPTHGV)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_ROOT_EXTINCTIONGV)
  NULLIFY(DATA_ISBA_MODEL(J)%XPAR_H_VEG)
ENDDO
DATA_ISBA_MODEL(:)%NTIME=0
DATA_ISBA_MODEL(:)%LDATA_MIXPAR=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_VEGTYPE=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_LAI=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_H_TREE=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_DG=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_DICE=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_GROUND_DEPTH=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_ROOT_DEPTH=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_ROOT_EXTINCTION=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_ROOT_LIN=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_ROOTFRAC=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_VEG=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_Z0=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_EMIS=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_ALBNIR_VEG=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_ALBVIS_VEG=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_ALBUV_VEG=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_RSMIN=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_GAMMA=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_WRMAX_CF=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_CV=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_RGL=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_Z0_O_Z0H=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_BSLAI=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_LAIMIN=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_SEFOLD=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_GMES=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_RE25=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_STRESS=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_F2I=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_GC=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_DMAX=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_CE_NITRO=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_CF_NITRO=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_CNA_NITRO=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_ALBNIR_SOIL=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_ALBVIS_SOIL=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_ALBUV_SOIL=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_IRRIG=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_WATSUP=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_GNDLITTER=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_LAIGV=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_Z0LITTER=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_RSMINGV=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_GAMMAGV=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_WRMAX_CFGV=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_RGLGV=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_ZF_TALLVEG=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_ROOTFRACGV=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_ROOT_DEPTHGV=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_ROOT_EXTINCTIONGV=.FALSE.
DATA_ISBA_MODEL(:)%LDATA_H_VEG=.FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_DATA_ISBA_N:DATA_ISBA_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE DATA_ISBA_ALLOC

SUBROUTINE DATA_ISBA_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DATA_ISBA_N:DATA_ISBA_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(DATA_ISBA_MODEL)) DEALLOCATE(DATA_ISBA_MODEL)
IF (ASSOCIATED(DATA_ISBA)) NULLIFY(DATA_ISBA)
IF (LHOOK) CALL DR_HOOK("MODD_DATA_ISBA_N:DATA_ISBA_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE DATA_ISBA_DEALLO

END MODULE MODD_DATA_ISBA_n
