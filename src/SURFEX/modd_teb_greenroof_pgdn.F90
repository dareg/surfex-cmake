!##################
MODULE MODD_TEB_GREENROOF_PGD_n
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
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE


TYPE TEB_GREENROOF_PGD_t
!-------------------------------------------------------------------------------
!
! Mask and number of grid elements containing patches/tiles:
!
  REAL, POINTER, DIMENSION(:,:) :: XVEGTYPE          ! fraction of each vegetation type for
!                                                    ! each grid mesh                          (-)
!
!-------------------------------------------------------------------------------
!
! Averaged Surface radiative parameters:
!
  REAL, POINTER, DIMENSION(:)   :: XALBNIR_DRY       ! dry soil near-infra-red albedo          (-)
  REAL, POINTER, DIMENSION(:)   :: XALBVIS_DRY       ! dry soil visible albedo                 (-)
  REAL, POINTER, DIMENSION(:)   :: XALBUV_DRY        ! dry soil UV albedo                      (-)
  REAL, POINTER, DIMENSION(:)   :: XALBNIR_WET       ! wet soil near-infra-red albedo          (-)
  REAL, POINTER, DIMENSION(:)   :: XALBVIS_WET       ! wet soil visible albedo                 (-)
  REAL, POINTER, DIMENSION(:)   :: XALBUV_WET        ! wet soil UV albedo                      (-)
  REAL, POINTER, DIMENSION(:)   :: XALBNIR_SOIL      ! soil near-infra-red albedo              (-)
  REAL, POINTER, DIMENSION(:)   :: XALBVIS_SOIL      ! soil visible albedo                     (-)
  REAL, POINTER, DIMENSION(:)   :: XALBUV_SOIL       ! soil UV albedo                          (-)
!
!-------------------------------------------------------------------------------
!
! Input Parameters, per patch:
!
! - vegetation + bare soil:
!
  REAL, POINTER, DIMENSION(:)   :: XZ0_O_Z0H         ! ratio of surface roughness lengths
!                                                    ! (momentum to heat)                      (-)

!
! - vegetation:
!
  REAL, POINTER, DIMENSION(:)   :: XALBNIR_VEG       ! vegetation near-infra-red albedo        (-)
  REAL, POINTER, DIMENSION(:)   :: XALBVIS_VEG       ! vegetation visible albedo               (-)
  REAL, POINTER, DIMENSION(:)   :: XALBUV_VEG        ! vegetation UV albedo                    (-)
!
! - vegetation: default option (Jarvis) and general parameters:
!
  REAL, POINTER, DIMENSION(:)   :: XWRMAX_CF         ! coefficient for maximum water 
                                                     ! interception 
                                                     ! storage capacity on the vegetation      (-)
  REAL, POINTER, DIMENSION(:)   :: XRSMIN            ! minimum stomatal resistance             (s/m)
  REAL, POINTER, DIMENSION(:)   :: XGAMMA            ! coefficient for the calculation
                                                     ! of the surface stomatal
                                                     ! resistance
  REAL, POINTER, DIMENSION(:)   :: XCV               ! vegetation thermal inertia coefficient  (K m2/J)
  REAL, POINTER, DIMENSION(:)   :: XRGL              ! maximum solar radiation
                                                     ! usable in photosynthesis                (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XROOTFRAC         ! root fraction profile ('DIF' option)
!
!-------------------------------------------------------------------------------
!
! - vegetation: Ags parameters ('AGS', 'LAI', 'AST', 'LST', 'NIT', 'NCB' options)
!
  REAL, POINTER, DIMENSION(:)   :: XABC              ! abscissa needed for integration
                                                     ! of net assimilation and stomatal
                                                     ! conductance over canopy depth           (-)
  REAL, POINTER, DIMENSION(:)   :: XPOI              ! Gaussian weights for integration
                                                     ! of net assimilation and stomatal
                                                     ! conductance over canopy depth           (-)
  REAL, POINTER, DIMENSION(:)   :: XBSLAI            ! ratio d(biomass)/d(lai)                 (kg/m2)
  REAL, POINTER, DIMENSION(:)   :: XLAIMIN           ! minimum LAI (Leaf Area Index)           (m2/m2)
  REAL, POINTER, DIMENSION(:)   :: XSEFOLD           ! e-folding time for senescence           (s)
  REAL, POINTER, DIMENSION(:)   :: XH_TREE           ! height of trees                         (m)
  REAL, POINTER, DIMENSION(:)   :: XANF              ! total assimilation over canopy          (
  REAL, POINTER, DIMENSION(:)   :: XANMAX            ! maximum photosynthesis rate             (
  REAL, POINTER, DIMENSION(:)   :: XFZERO            ! ideal value of F, no photo- 
                                                     ! respiration or saturation deficit       (
  REAL, POINTER, DIMENSION(:)   :: XEPSO             ! maximum initial quantum use             
                                                     ! efficiency                              (mg J-1 PAR)
  REAL, POINTER, DIMENSION(:)   :: XGAMM             ! CO2 conpensation concentration          (ppm)
  REAL, POINTER, DIMENSION(:)   :: XQDGAMM           ! Q10 function for CO2 conpensation 
                                                     ! concentration                           (-)
  REAL, POINTER, DIMENSION(:)   :: XGMES             ! mesophyll conductance                   (m s-1)
  REAL, POINTER, DIMENSION(:)   :: XRE25             ! Ecosystem respiration parameter         (kg/kg.m.s-1)
  REAL, POINTER, DIMENSION(:)   :: XQDGMES           ! Q10 function for mesophyll conductance  (-)
  REAL, POINTER, DIMENSION(:)   :: XT1GMES           ! reference temperature for computing 
                                                     ! compensation concentration function for 
                                                     ! mesophyll conductance: minimum
                                                     ! temperature                             (K)
  REAL, POINTER, DIMENSION(:)   :: XT2GMES           ! reference temperature for computing 
                                                     ! compensation concentration function for 
                                                     ! mesophyll conductance: maximum
                                                     ! temperature                             (K)
  REAL, POINTER, DIMENSION(:)   :: XAMAX             ! leaf photosynthetic capacity            (mg m-2 s-1)
  REAL, POINTER, DIMENSION(:)   :: XQDAMAX           ! Q10 function for leaf photosynthetic 
                                                     ! capacity                                (-)
  REAL, POINTER, DIMENSION(:)   :: XT1AMAX           ! reference temperature for computing 
                                                     ! compensation concentration function for 
                                                     ! leaf photosynthetic capacity: minimum
                                                     ! temperature                             (K)
  REAL, POINTER, DIMENSION(:)   :: XT2AMAX           ! reference temperature for computing 
                                                     ! compensation concentration function for 
                                                     ! leaf photosynthetic capacity: maximum
                                                     ! temperature                             (K)
!                                      
!-------------------------------------------------------------------------------
!
! - vegetation: Ags Stress parameters ('AST', 'LST', 'NIT', 'NCB' options)
!
  LOGICAL, POINTER, DIMENSION(:):: LSTRESS           ! vegetation response (offensive/defensive)
  REAL, POINTER, DIMENSION(:)   :: XF2I              ! critical normilized soil water 
                                                     ! content for stress parameterisation
  REAL, POINTER, DIMENSION(:)   :: XGC               ! cuticular conductance                   (m s-1)
  REAL, POINTER, DIMENSION(:)   :: XAH               ! coefficients for herbaceous water stress 
                                                     ! response (offensive or defensive)       (log(mm/s))
  REAL, POINTER, DIMENSION(:)   :: XBH               ! coefficients for herbaceous water stress 
                                                     ! response (offensive or defensive)       (-)
  REAL, POINTER, DIMENSION(:)   :: XDMAX             ! maximum air saturation deficit
                                                     ! tolerate by vegetation                  (kg/kg)
!
!-------------------------------------------------------------------------------
!
! - vegetation: Ags Nitrogen-model parameters ('NIT', 'NCB' option)
!
  REAL, POINTER, DIMENSION(:)    :: XCE_NITRO        ! leaf aera ratio sensitivity to 
                                                     ! nitrogen concentration                (m2/kg)
  REAL, POINTER, DIMENSION(:)    :: XCF_NITRO        ! lethal minimum value of leaf area
                                                     ! ratio                                 (m2/kg)
  REAL, POINTER, DIMENSION(:)    :: XCNA_NITRO       ! nitrogen concentration of active 
                                                     ! biomass                               (kg/kg)
  REAL, POINTER, DIMENSION(:)    :: XBSLAI_NITRO     ! biomass/LAI ratio from nitrogen 
                                                     ! decline theory                        (kg/m2)
!
!-------------------------------------------------------------------------------
!
! - soil: primary parameters
!
  REAL, POINTER, DIMENSION(:,:)  :: XOM_GR           ! green roof OM fraction (-)
  REAL, POINTER, DIMENSION(:,:)  :: XSAND_GR         ! green roof sand fraction of the non-OM part (-)
  REAL, POINTER, DIMENSION(:,:)  :: XCLAY_GR         ! green roof clay fraction of the non-OM part (-)
  REAL, POINTER, DIMENSION(:)    :: XRUNOFFB_GR      ! green roof sub-grid surface runoff slope parameter (-)
  REAL, POINTER, DIMENSION(:)    :: XWDRAIN_GR       ! green roof continuous drainage parameter           (-)
  REAL, POINTER, DIMENSION(:)    :: XTAUICE          ! soil freezing characteristic timescale  (s)
  REAL, POINTER, DIMENSION(:)    :: XGAMMAT          ! 'Force-Restore' timescale when using a
                                                     ! prescribed lower boundary temperature   (1/days)
  REAL, POINTER, DIMENSION(:,:)  :: XDG              ! soil layer thicknesses                  (m)
                                                     ! NOTE: in Force-Restore mode, the 
                                                     ! uppermost layer thickness is superficial
                                                     ! and is only explicitly used for soil 
                                                     ! water phase changes                     (m)
  REAL, POINTER, DIMENSION(:)    :: XRUNOFFD         ! depth over which sub-grid runoff is
                                                     ! computed: in Force-Restore this is the
                                                     ! total soil column ('2-L'), or root zone
                                                     ! ('3-L'). For the 'DIF' option, it can
                                                     ! be any depth within soil column         (m)
!
  REAL, POINTER, DIMENSION(:,:)  :: XSOILWGHT      ! ISBA-DIF: weights for vertical
  REAL, POINTER, DIMENSION(:,:)  :: XDZG           ! soil layers thicknesses (DIF option)
  REAL, POINTER, DIMENSION(:,:)  :: XDZDIF         ! distance between consecuative layer mid-points (DIF option)
!
  INTEGER, POINTER, DIMENSION(:) :: NWG_LAYER      ! Number of soil moisture layers for DIF
  REAL, POINTER, DIMENSION(:)    :: XDROOT         ! effective root depth for DIF (m)
  REAL, POINTER, DIMENSION(:)    :: XDG2           ! root depth for DIF as 3-L (m)
!-------------------------------------------------------------------------------
!
! - soil: Secondary parameters: hydrology
!
  REAL, POINTER, DIMENSION(:)    :: XC1SAT           ! 'Force-Restore' C1 coefficient at 
                                                     ! saturation                              (-)
  REAL, POINTER, DIMENSION(:)    :: XC2REF           ! 'Force-Restore' reference value of C2   (-)
  REAL, POINTER, DIMENSION(:,:)  :: XC3              ! 'Force-Restore' C3 drainage coefficient (m)
  REAL, POINTER, DIMENSION(:)    :: XC4B             ! 'Force-Restore' sub-surface vertical 
                                                     ! diffusion coefficient (slope parameter) (-)
  REAL, POINTER, DIMENSION(:)    :: XC4REF           ! 'Force-Restore' sub-surface vertical 
                                                     ! diffusion coefficient                   (-)
  REAL, POINTER, DIMENSION(:)    :: XACOEF           ! 'Force-Restore' surface vertical 
                                                     ! diffusion coefficient                   (-)
  REAL, POINTER, DIMENSION(:)    :: XPCOEF           ! 'Force-Restore' surface vertical 
                                                     ! diffusion coefficient                   (-)
  REAL, POINTER, DIMENSION(:,:)  :: XWFC             ! field capacity volumetric water content
                                                     ! profile                             (m3/m3)
  REAL, POINTER, DIMENSION(:,:)  :: XWWILT           ! wilting point volumetric water content 
                                                     ! profile         
  REAL, POINTER, DIMENSION(:,:)  :: XWSAT            ! porosity profile                      (m3/m3) 
  REAL, POINTER, DIMENSION(:,:)  :: XBCOEF           ! soil water CH78 b-parameter             (-)
  REAL, POINTER, DIMENSION(:,:)  :: XCONDSAT         ! hydraulic conductivity at saturation    (m/s)
  REAL, POINTER, DIMENSION(:,:)  :: XMPOTSAT         ! matric potential at saturation          (m)
!
  REAL, POINTER, DIMENSION(:)    :: XPCPS
  REAL, POINTER, DIMENSION(:)    :: XPLVTT
  REAL, POINTER, DIMENSION(:)    :: XPLSTT 
!
!-------------------------------------------------------------------------------
!
! - soil: Secondary parameters: thermal 
!
  REAL, POINTER, DIMENSION(:)    :: XCGSAT           ! soil thermal inertia coefficient at 
                                                     ! saturation                              (K m2/J)
  REAL, POINTER, DIMENSION(:,:)  :: XHCAPSOIL        ! soil heat capacity                      (J/K/m3)
  REAL, POINTER, DIMENSION(:,:)  :: XCONDDRY         ! soil dry thermal conductivity           (W/m/K)
  REAL, POINTER, DIMENSION(:,:)  :: XCONDSLD         ! soil solids thermal conductivity        (W/m/K)
  REAL, POINTER, DIMENSION(:)    :: XTDEEP           ! prescribed deep soil temperature 
                                                     ! (optional)                              (K)
!
! - SGH scheme
!
  REAL, POINTER, DIMENSION(:)    :: XD_ICE          ! depth of the soil column for the calculation
                                                    ! of the frozen soil fraction (m)
  REAL, POINTER, DIMENSION(:)    :: XKSAT_ICE       ! hydraulic conductivity at saturation
                                                    ! over frozen area (m s-1)
!-------------------------------------------------------------------------------
!
END TYPE TEB_GREENROOF_PGD_t


TYPE(TEB_GREENROOF_PGD_t),     ALLOCATABLE, TARGET, SAVE :: TEB_GREENROOF_PGD_MODEL(:)

TYPE(TEB_GREENROOF_PGD_t), POINTER :: TEB_GREENROOF_PGD => NULL()
!$OMP THREADPRIVATE(TEB_GREENROOF_PGD)

CONTAINS



SUBROUTINE TEB_GREENROOF_PGD_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_TEB_GREENROOF_N:TEB_GREENROOF_PGD_GOTO_MODEL',0,ZHOOK_HANDLE)

TEB_GREENROOF_PGD => TEB_GREENROOF_PGD_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_TEB_GREENROOF_N:TEB_GREENROOF_PGD_GOTO_MODEL',1,ZHOOK_HANDLE)
!
END SUBROUTINE TEB_GREENROOF_PGD_GOTO_MODEL

SUBROUTINE TEB_GREENROOF_PGD_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GREENROOF_N:TEB_GREENROOF_PGD_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(TEB_GREENROOF_PGD_MODEL(KMODEL))
DO J=1,KMODEL
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XVEGTYPE)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XALBNIR_DRY)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XALBVIS_DRY)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XALBUV_DRY)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XALBNIR_WET)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XALBVIS_WET)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XALBUV_WET)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XALBNIR_SOIL)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XALBVIS_SOIL)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XALBUV_SOIL)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XZ0_O_Z0H)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XALBNIR_VEG)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XALBVIS_VEG)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XALBUV_VEG)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XWRMAX_CF)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XRSMIN)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XGAMMA)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XCV)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XRGL)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XROOTFRAC)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XABC)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XPOI)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XBSLAI)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XLAIMIN)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XSEFOLD)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XH_TREE)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XANF)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XANMAX)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XFZERO)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XEPSO)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XGAMM)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XQDGAMM)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XGMES)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XRE25)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XQDGMES)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XT1GMES)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XT2GMES)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XAMAX)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XQDAMAX)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XT1AMAX)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XT2AMAX)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%LSTRESS)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XF2I)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XGC)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XAH)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XBH)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XDMAX)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XCE_NITRO)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XCF_NITRO)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XCNA_NITRO)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XBSLAI_NITRO)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XOM_GR)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XSAND_GR)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XCLAY_GR)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XRUNOFFB_GR)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XWDRAIN_GR)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XTAUICE)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XGAMMAT)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XDG)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XRUNOFFD)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XSOILWGHT)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XDZG)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XDZDIF)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%NWG_LAYER)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XDROOT)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XDG2)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XPCPS)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XPLVTT)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XPLSTT)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XC1SAT)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XC2REF)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XC3)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XC4B)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XC4REF)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XACOEF)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XPCOEF)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XWFC)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XWWILT)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XWSAT)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XBCOEF)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XCONDSAT)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XMPOTSAT)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XCGSAT)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XHCAPSOIL)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XCONDDRY)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XCONDSLD)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XTDEEP)
NULLIFY(TEB_GREENROOF_PGD_MODEL(J)%XD_ICE)
ENDDO
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GREENROOF_N:TEB_GREENROOF_PGD_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_GREENROOF_PGD_ALLOC

SUBROUTINE TEB_GREENROOF_PGD_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GREENROOF_N:TEB_GREENROOF_PGD_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(TEB_GREENROOF_PGD_MODEL)) DEALLOCATE(TEB_GREENROOF_PGD_MODEL)
IF (ASSOCIATED(TEB_GREENROOF_PGD)) NULLIFY(TEB_GREENROOF_PGD)
IF (LHOOK) CALL DR_HOOK("MODD_TEB_GREENROOF_N:TEB_GREENROOF_PGD_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_GREENROOF_PGD_DEALLO

END MODULE MODD_TEB_GREENROOF_PGD_n
