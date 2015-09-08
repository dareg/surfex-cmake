!     ################
      MODULE MODD_TEB_VEG_n
!     ################
!
!!****  *MODD_TEB_VEG_n - declaration of options and parameters for urban vegetation
!!                        (for parameters common to all types of urban vegetation)
!!
!!    PURPOSE
!!    -------
!     Declaration of options and surface parameters
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
!!      C. de Munck & A. Lemonsu   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       07/2012
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE


TYPE TEB_VEG_OPTIONS_t
! ISBA options common of all types of urban vegetation
!
!
  LOGICAL                          :: LCANOPY_DRAG ! T: drag activated in SBL scheme within the canopy
                                                   ! F: no drag activated in SBL atmospheric layers
!
  LOGICAL                          :: LVEGUPD      ! T = update vegetation parameters every decade
                                                   ! F = keep vegetation parameters constant in time
!
  LOGICAL                          :: LTR_ML
!  
  LOGICAL                          :: LNITRO_DILU  ! nitrogen dilution fct of CO2 (Calvet et al. 2008)
!-------------------------------------------------------------------------------
!
  CHARACTER(LEN=3)                 :: CISBA       ! type of ISBA version:
                                                  ! '2-L' (default)
                                                  ! '3-L'
                                                  ! 'DIF'
!
  CHARACTER(LEN=4)                 :: CROUGH      ! type of roughness length
                                                  ! 'Z01D'
                                                  ! 'Z04D'
!
  CHARACTER(LEN=4)                 :: CPEDOTF     ! NOTE: Only used when HISBA = DIF
                                                  ! 'CH78' = Clapp and Hornberger 1978 for BC (Default)
                                                  ! 'CO84' = Cosby et al. 1988 for BC
                                                  ! 'CP88' = Carsel and Parrish 1988 for VG
                                                  ! 'WO99' = Wosten et al. 1999 for VG
!
  CHARACTER(LEN=3)                 :: CPHOTO      ! type of photosynthesis
                                                  ! 'NON'
                                                  ! 'AGS'
                                                  ! 'LAI'
                                                  ! 'LST'
                                                  ! 'AST'
                                                  ! 'NIT'
                                                  ! 'NCB'
!
  CHARACTER(LEN=4)                 :: CALBEDO     ! albedo type
                                                  ! 'DRY ' 
                                                  ! 'EVOL' 
                                                  ! 'WET ' 
                                                  ! 'USER' 
!
  CHARACTER(LEN=4)                 :: CSCOND      ! Thermal conductivity
                                                  ! 'DEF ' = DEFault: NP89 implicit method
                                                  ! 'PL98' = Peters-Lidard et al. 1998 used
                                                  ! for explicit computation of CG
!
  CHARACTER(LEN=4)                 :: CC1DRY      ! C1 formulation for dry soils
                                                  ! 'DEF ' = DEFault: Giard-Bazile formulation
                                                  ! 'GB93' = Giordani 1993, Braud 1993 
                                                  !discontinuous at WILT
!
  CHARACTER(LEN=3)                 :: CSOILFRZ    ! soil freezing-physics option
                                                  ! 'DEF' = Default (Boone et al. 2000; 
                                                  !        Giard and Bazile 2000)
                                                  ! 'LWT' = Phase changes as above,
                                                  !         but relation between unfrozen 
                                                  !         water and temperature considered
!                            NOTE that when using the YISBA='DIF' multi-layer soil option,
!                            the 'LWT' method is used. It is only an option
!                            when using the force-restore soil method ('2-L' or '3-L')
!
  CHARACTER(LEN=4)                 :: CDIFSFCOND  ! Mulch effects
                                                  ! 'MLCH' = include the insulating effect of
                                                  ! leaf litter/mulch on the surf. thermal cond.
                                                  ! 'DEF ' = no mulch effect
!                           NOTE: Only used when YISBA = DIF
!
  CHARACTER(LEN=3)                 :: CSNOWRES    ! Turbulent exchanges over snow
                                                  ! 'DEF' = Default: Louis (ISBA)
                                                  ! 'RIL' = Maximum Richardson number limit
                                                  !         for stable conditions ISBA-SNOW3L
                                                  !         turbulent exchange option
!                                           
  CHARACTER(LEN=3)                 :: CRESPSL     ! Soil respiration
                                                  ! 'DEF' = Default: Norman (1992)
                                                  ! 'PRM' = New Parameterization
                                                  ! 'CNT' = CENTURY model (Gibelin 2007)
!                                           
  CHARACTER(LEN=3)                 :: CCPSURF     ! specific heat at surface
                                                  ! 'DRY' = default value (dry Cp)
                                                  ! 'HUM' = Cp as a fct of specific humidity
! - SGH scheme and vertical hydrology
!                                                     
  CHARACTER(LEN=4)                 :: CRUNOFF     ! surface runoff formulation
                                                  ! 'WSAT'
                                                  ! 'DT92'
                                                  ! 'SGH ' Topmodel
!                                                     
  CHARACTER(LEN=3)                 :: CKSAT       ! ksat
                                                  ! 'DEF' = default value 
                                                  ! 'SGH' = profil exponentiel
!
  LOGICAL                          :: LSOC        ! soil organic carbon effect
!                                                 ! FALSE = default value 
!                                                 ! TRUE  = SOC profil
!
  CHARACTER(LEN=3)                 :: CRAIN       ! Rainfall spatial distribution
                                                  ! 'DEF' = No rainfall spatial distribution
                                                  ! 'SGH' = Rainfall exponential spatial distribution
!
  CHARACTER(LEN=3)                 :: CHORT       ! Horton runoff
                                                  ! 'DEF' = no Horton runoff
                                                  ! 'SGH' = Horton runoff
!
! -----------------------------------------------------------------------------------------------------------
!
  INTEGER                          :: NNBIOMASS   ! number of biomass pools
  REAL                             :: XCGMAX      ! maximum soil heat capacity (=2.E-5)
  REAL                             :: XCDRAG      ! drag coefficient in canopy
  REAL                             :: XTSTEP      !  time step  
!
END TYPE TEB_VEG_OPTIONS_t

TYPE(TEB_VEG_OPTIONS_t), ALLOCATABLE, TARGET, SAVE :: TEB_VEG_OPTIONS_MODEL(:)

TYPE(TEB_VEG_OPTIONS_t), POINTER :: TEB_VEG_OPTIONS => NULL()
!$OMP THREADPRIVATE(TEB_VEG_OPTIONS)

CONTAINS
!----------------------------------------------------------------------------

SUBROUTINE TEB_VEG_OPTIONS_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_TEB_VEG_N:TEB_VEG_OPTIONS_GOTO_MODEL',0,ZHOOK_HANDLE)

TEB_VEG_OPTIONS => TEB_VEG_OPTIONS_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_TEB_VEG_N:TEB_VEG_OPTIONS_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE TEB_VEG_OPTIONS_GOTO_MODEL

SUBROUTINE TEB_VEG_OPTIONS_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J,JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_VEG_N:TEB_VEG_OPTIONS_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(TEB_VEG_OPTIONS_MODEL(KMODEL))
TEB_VEG_OPTIONS => TEB_VEG_OPTIONS_MODEL(KMODEL)
TEB_VEG_OPTIONS_MODEL(:)%LCANOPY_DRAG=.FALSE.
TEB_VEG_OPTIONS_MODEL(:)%LVEGUPD=.FALSE. 
TEB_VEG_OPTIONS_MODEL(:)%LNITRO_DILU=.FALSE. 
TEB_VEG_OPTIONS_MODEL(:)%LTR_ML=.FALSE.
TEB_VEG_OPTIONS_MODEL(:)%CISBA=' '
TEB_VEG_OPTIONS_MODEL(:)%CROUGH=' '
TEB_VEG_OPTIONS_MODEL(:)%CSCOND=' '
TEB_VEG_OPTIONS_MODEL(:)%CPEDOTF=' '
TEB_VEG_OPTIONS_MODEL(:)%CPHOTO=' '
TEB_VEG_OPTIONS_MODEL(:)%CALBEDO=' '
TEB_VEG_OPTIONS_MODEL(:)%CC1DRY=' '
TEB_VEG_OPTIONS_MODEL(:)%CSOILFRZ=' '
TEB_VEG_OPTIONS_MODEL(:)%CDIFSFCOND=' '
TEB_VEG_OPTIONS_MODEL(:)%CSNOWRES=' '
TEB_VEG_OPTIONS_MODEL(:)%CRESPSL=' '
TEB_VEG_OPTIONS_MODEL(:)%CCPSURF=' '
TEB_VEG_OPTIONS_MODEL(:)%CRUNOFF=' '
TEB_VEG_OPTIONS_MODEL(:)%CKSAT=' '
TEB_VEG_OPTIONS_MODEL(:)%LSOC=.FALSE.
TEB_VEG_OPTIONS_MODEL(:)%CRAIN=' '
TEB_VEG_OPTIONS_MODEL(:)%CHORT=' '
TEB_VEG_OPTIONS_MODEL(:)%NNBIOMASS=0
TEB_VEG_OPTIONS_MODEL(:)%XCGMAX=0.
TEB_VEG_OPTIONS_MODEL(:)%XCDRAG=0.
TEB_VEG_OPTIONS_MODEL(:)%XTSTEP=0.
IF (LHOOK) CALL DR_HOOK("MODD_TEB_VEG_N:TEB_VEG_OPTIONS_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_VEG_OPTIONS_ALLOC

SUBROUTINE TEB_VEG_OPTIONS_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_VEG_N:TEB_VEG_OPTIONS_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(TEB_VEG_OPTIONS_MODEL)) DEALLOCATE(TEB_VEG_OPTIONS_MODEL)
IF (ASSOCIATED(TEB_VEG_OPTIONS)) NULLIFY(TEB_VEG_OPTIONS)
IF (LHOOK) CALL DR_HOOK("MODD_TEB_VEG_N:TEB_VEG_OPTIONS_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_VEG_OPTIONS_DEALLO

!----------------------------------------------------------------------------

END MODULE MODD_TEB_VEG_n
