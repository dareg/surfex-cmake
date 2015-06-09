!######################
MODULE MODD_PACK_DIAG_ISBA
!######################
!
!!****  *MODD_PACK_DIAG_ISBA - declaration of packed diagnostics for ISBA scheme
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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2004
!!      Modified       10/2004 by P. Le Moigne: add Halstead coefficient
!!      Modified       11/2009 by S. Senesi: add precipitation intercepted by the vegetation (XP_RRVEG)
!!      Modified       04-09 by A.L. Gibelin  : Add carbon diagnostics
!!      Modified       10-14 by P. Samuelsson: MEB
!
!*       0.   DECLARATIONS
!             ------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!------------------------------------------------------------------------------
!
TYPE PACK_DIAG_ISBA_t

  INTEGER :: NSIZE_SIMPLE
  INTEGER :: NSIZE_GROUND
  INTEGER :: NSIZE_SNOW
  INTEGER :: NSIZE_KSW
  INTEGER :: NSIZE_ABC
  INTEGER :: NSIZE_0
  INTEGER :: NSIZE_00
  REAL, POINTER, DIMENSION(:,:) :: XBLOCK_SIMPLE
  REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_GROUND
  REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_SNOW
  REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_KSW
  REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_ABC
  REAL, POINTER, DIMENSION(:,:) :: XBLOCK_0
  REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_00
!
  REAL, POINTER, DIMENSION(:) :: XP_RNSNOW    ! net radiative flux from snow (ISBA-ES:3-L)    (W/m2)
  REAL, POINTER, DIMENSION(:) :: XP_HSNOW     ! sensible heat flux from snow (ISBA-ES:3-L)    (W/m2)
  REAL, POINTER, DIMENSION(:) :: XP_HPSNOW    ! heat release from rainfall (ISBA-ES:3-L)      (W/m2)
  REAL, POINTER, DIMENSION(:) :: XP_GFLUXSNOW ! net surface energy flux into snowpack      
!                                               ! (ISBA-ES:3-L)                                 (W/m2)
  REAL, POINTER, DIMENSION(:) :: XP_USTARSNOW ! friction velocity  over snow 
!                                               ! (ISBA-ES:3-L)                                 (m/s)
  REAL, POINTER, DIMENSION(:) :: XP_GRNDFLUX  ! soil/snow interface heat flux (ISBA-ES:3-L)   (W/m2)
  REAL, POINTER, DIMENSION(:) :: XP_SRSFC     ! snowfall over snowpack (ISBA-ES:3-L)          (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XP_RRSFC     ! rainfall over snowpack (ISBA-ES:3-L)          (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XP_LESL      ! snowpack evaporation (ISBA-ES:3-L)            (W/m2)
  REAL, POINTER, DIMENSION(:) :: XP_CDSNOW    ! snow drag coefficient (ISBA-ES:3-L)           (-)
  REAL, POINTER, DIMENSION(:) :: XP_CHSNOW    ! heat turbulent transfer coefficient 
!                                               ! (ISBA-ES:3-L)                                 (-)
  REAL, POINTER, DIMENSION(:,:)::XP_SNOWTEMP  ! snow temperature profile (ISBA-ES:3-L)        (K)
  REAL, POINTER, DIMENSION(:,:)::XP_SNOWLIQ   ! snow liquid water profile (ISBA-ES:3-L)       (m)
  REAL, POINTER, DIMENSION(:,:)::XP_SNOWDZ    ! snow layer thicknesses                        (m)
  REAL, POINTER, DIMENSION(:) :: XP_SNOWHMASS ! heat content change due to mass
!                                               ! changes in snowpack: for budget
!                                               ! calculations only. (ISBA-ES:3-L)              (J/m2)
  REAL, POINTER, DIMENSION(:) :: XP_RN_ISBA   ! net radiative flux from snow-free surface 
!                                               ! (ISBA-ES:3-L)                                 (W/m2) 
  REAL, POINTER, DIMENSION(:) :: XP_H_ISBA    ! sensible heat flux from snow-free surface 
!                                               ! (ISBA-ES:3-L)                                 (W/m2) 
  REAL, POINTER, DIMENSION(:) :: XP_LEG_ISBA  ! baresoil evaporation from snow-free surface 
!                                               ! (ISBA-ES:3-L)                                 (W/m2) 
  REAL, POINTER, DIMENSION(:) :: XP_LEGI_ISBA ! baresoil sublimation from snow-free surface 
!                                               ! (ISBA-ES:3-L)                                 (W/m2) 
  REAL, POINTER, DIMENSION(:) :: XP_LEV_ISBA  ! total evapotranspiration from vegetation over 
!                                               ! snow-free surface (ISBA-ES:3-L)               (W/m2) 
  REAL, POINTER, DIMENSION(:) :: XP_LETR_ISBA ! transpiration from snow-free surface 
!                                               ! (ISBA-ES:3-L)                                 (W/m2)
  REAL, POINTER, DIMENSION(:) :: XP_USTAR_ISBA! friction velocity from snow-free 
!                                               ! surface (ISBA-ES:3-L)                         (m/s)
  REAL, POINTER, DIMENSION(:) :: XP_LER_ISBA  ! evaporation from canopy water interception 
!                                               ! store over snow-free surface (ISBA-ES:3-L)    (W/m2)
  REAL, POINTER, DIMENSION(:) :: XP_LE_ISBA   ! total latent heat flux from snow-free surface 
  REAL, POINTER, DIMENSION(:) :: XP_LEI_ISBA  ! sublimation latent heat flux from snow-free surface 
!                                               ! (ISBA-ES:3-L)                                 (W/m2) 
  REAL, POINTER, DIMENSION(:) :: XP_GFLUX_ISBA! net energy flux into the snow-free surface 
!                                               ! (ISBA-ES:3-L)                                 (W/m2) 
  REAL, POINTER, DIMENSION(:) :: XP_MELTADV   ! advective energy from snow melt water 
!                                               ! (ISBA-ES:3-L)                                 (W/m2)
  REAL, POINTER, DIMENSION(:) :: XP_CH        ! thermal diffusion coefficient                 (W/s)
  REAL, POINTER, DIMENSION(:) :: XP_CE        ! transfer coefficient for vapor                (W/s/K)
  REAL, POINTER, DIMENSION(:) :: XP_CD        ! drag coefficient                              (-)
  REAL, POINTER, DIMENSION(:) :: XP_CDN       ! neutral drag coefficient                      (-)
  REAL, POINTER, DIMENSION(:) :: XP_RI        ! Bulk-Richardson number                        (-)
  REAL, POINTER, DIMENSION(:) :: XP_HU        ! area averaged surface humidity coefficient    (-)
  REAL, POINTER, DIMENSION(:) :: XP_HUG       ! baresoil surface humidity coefficient         (-)
  REAL, POINTER, DIMENSION(:) :: XP_HV        ! Halstead coefficient                          (-)
!
  REAL, POINTER, DIMENSION(:) :: XP_ALBT      ! Total Albedo                                  (-)
!
  REAL, POINTER, DIMENSION(:) :: XP_RN        ! net radiation at surface                      (W/m2)
  REAL, POINTER, DIMENSION(:) :: XP_H         ! sensible heat flux                            (W/m2)
  REAL, POINTER, DIMENSION(:) :: XP_LEG       ! baresoil evaporation                          (W/m2)
  REAL, POINTER, DIMENSION(:) :: XP_LEGI      ! baresoil sublimation                          (W/m2)
  REAL, POINTER, DIMENSION(:) :: XP_LEV       ! total evapotranspiration from vegetation      (W/m2)
  REAL, POINTER, DIMENSION(:) :: XP_LES       ! snow sublimation                              (W/m2)
  REAL, POINTER, DIMENSION(:) :: XP_LER       ! evaporation from canopy water interception 
!                                               ! store                                         (W/m2)
  REAL, POINTER, DIMENSION(:) :: XP_LETR      ! transpiration                                 (W/m2)
  REAL, POINTER, DIMENSION(:) :: XP_EVAP      ! evapotranspiration                            (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XP_SUBL      ! sublimation                                   (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XP_SNDRIFT   ! blowing snow sublimation (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XP_LEI       ! sublimation latent heat flux                  (W/m2)
  REAL, POINTER, DIMENSION(:) :: XP_GFLUX     ! net soil-vegetation flux                      (W/m2)
  REAL, POINTER, DIMENSION(:) :: XP_RESTORE   ! surface energy budget restore term            (W/m2)
  REAL, POINTER, DIMENSION(:) :: XP_DRAIN     ! soil drainage flux                            (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XP_QSB       ! lateral subsurface flux                       (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XP_RUNOFF    ! sub-grid and supersaturation runoff           (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XP_MELT      ! melting rate of the snow                      (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XP_SNOWFREE_ALB ! snow-free global albedo                    (-)
  REAL, POINTER, DIMENSION(:) :: XP_SNOWFREE_ALB_VEG ! snow-free global  albedo of vegetation
  REAL, POINTER, DIMENSION(:) :: XP_SNOWFREE_ALB_SOIL! snow-free soil albedo
  REAL, POINTER, DIMENSION(:) :: XP_Z0_WITH_SNOW ! total roughness length (including snow)    (m)
  REAL, POINTER, DIMENSION(:) :: XP_Z0H_WITH_SNOW! roughness length for heat (including snow) (m)
  REAL, POINTER, DIMENSION(:) :: XP_Z0EFF     ! effective roughness length (with relief added)(m)
!
  REAL, POINTER, DIMENSION(:,:)::XP_IACAN     ! PAR in the canopy at different gauss level    (micmolphot/m2/s)
!
  REAL, POINTER, DIMENSION(:) :: XP_CG        ! heat capacity of the ground
  REAL, POINTER, DIMENSION(:) :: XP_C1        ! coefficients for the moisure
  REAL, POINTER, DIMENSION(:) :: XP_C2        ! equation.
  REAL, POINTER, DIMENSION(:) :: XP_WGEQ      ! equilibrium volumetric water
!                                               ! content
  REAL, POINTER, DIMENSION(:) :: XP_CT        ! area-averaged heat capacity
  REAL, POINTER, DIMENSION(:) :: XP_RS        ! stomatal resistance                            (s/m)
!
!------------------------------------------------------------------------------
!
  REAL, POINTER, DIMENSION(:) :: XP_TS        ! Surface temperature                            (K)
  REAL, POINTER, DIMENSION(:) :: XP_TSRAD     ! Radiative surface temperature                  (K)
  REAL, POINTER, DIMENSION(:) :: XP_T2M       ! Air temperature       at 2 meters              (K)
  REAL, POINTER, DIMENSION(:) :: XP_Q2M       ! Air spec. humidity    at 2 meters              (kg/kg)
  REAL, POINTER, DIMENSION(:) :: XP_HU2M      ! Air rela. humidity    at 2 meters              (-)
  REAL, POINTER, DIMENSION(:) :: XP_ZON10M    ! zonal Wind at 10 meters                        (m/s)
  REAL, POINTER, DIMENSION(:) :: XP_MER10M    ! meridian Wind at 10 meters                     (m/s)
!
!------------------------------------------------------------------------------
!
  REAL, POINTER, DIMENSION(:)   :: XP_QS      ! humidity at surface                            (Kg/kg)
  REAL, POINTER, DIMENSION(:,:) :: XP_SWI     ! soil wetness index profile                     (-)
  REAL, POINTER, DIMENSION(:,:) :: XP_TSWI    ! total soil wetness index profile               (-)
!
!------------------------------------------------------------------------------
!
  REAL, POINTER, DIMENSION(:) :: XP_TWSNOW     ! total snow reservoir (kg/m2)
  REAL, POINTER, DIMENSION(:) :: XP_TDSNOW     ! total snow height (m)
!
!------------------------------------------------------------------------------
!
  REAL, POINTER, DIMENSION(:) :: XP_SWD       ! downward short wave radiation    (W/m2)
  REAL, POINTER, DIMENSION(:) :: XP_SWU       ! upward short wave radiation      (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XP_SWBD    ! downward short wave radiation by spectral band   (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XP_SWBU    ! upward short wave radiation by spectral band (W/m2)
  REAL, POINTER, DIMENSION(:) :: XP_LWD       ! downward long wave radiation     (W/m2)
  REAL, POINTER, DIMENSION(:) :: XP_LWU       ! upward long wave radiation       (W/m2)
  REAL, POINTER, DIMENSION(:) :: XP_FMU       ! horizontal momentum flux zonal   (m2/s2)
  REAL, POINTER, DIMENSION(:) :: XP_FMV       ! horizontal momentum flux meridian (m2/s2)
!
!------------------------------------------------------------------------------
!
  REAL, POINTER, DIMENSION(:) :: XP_HORT      ! sub-grid Horton runoff from the SGH scheme   (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XP_DRIP      ! dripping from the vegetation reservoir       (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XP_IFLOOD    ! flood infiltration                           (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XP_PFLOOD    ! precipitation intercepted by the floodplains (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XP_LE_FLOOD  ! flood evaporation                            (W/m2)
  REAL, POINTER, DIMENSION(:) :: XP_LEI_FLOOD ! frozen flood evaporation                     (W/m2)
  REAL, POINTER, DIMENSION(:) :: XP_ICEFLUX
  REAL, POINTER, DIMENSION(:) :: XP_RRVEG     ! precipitation intercepted by the vegetation   (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XP_IRRIG_FLUX! irrigation rate                               (kg/m2/s)
!
!------------------------------------------------------------------------------
!
  REAL, POINTER, DIMENSION(:) :: XP_GPP         ! Gross primary production (kgCO2/m2/s)
  REAL, POINTER, DIMENSION(:) :: XP_RESP_AUTO   ! Autotrophic respiration  (kgCO2/m2/s)
  REAL, POINTER, DIMENSION(:) :: XP_RESP_ECO    ! Ecosystem respiration    (kgCO2/m2/s)
  REAL, POINTER, DIMENSION(:) :: XP_FAPAR       ! Fapar of vegetation
  REAL, POINTER, DIMENSION(:) :: XP_FAPIR       ! Fapir of vegetation 
  REAL, POINTER, DIMENSION(:) :: XP_FAPAR_BS    ! Fapar of bare soil
  REAL, POINTER, DIMENSION(:) :: XP_FAPIR_BS    ! Fapir of bare soil
!
!------------------------------------------------------------------------------
!
!* diagnostic variables for multi-energy balance (MEB)
!  ---------------------------------------------------
!
  REAL, POINTER, DIMENSION(:) :: XP_SWUP          ! MEB: net *total* (surface) upwelling shortwave radiation to atmosphere [W/m2]
  REAL, POINTER, DIMENSION(:) :: XP_SWNET_V       ! MEB: net vegetation canopy shortwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XP_SWNET_G       ! MEB: net ground shortwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XP_SWNET_N       ! MEB: net snow shortwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XP_SWNET_NS      ! MEB: net snow shortwave radiation for *surface* layer 
                                                !     (i.e. net snow shortwave radiation less absorbed radiation) [W/m2]
  REAL, POINTER, DIMENSION(:) :: XP_LWUP          ! MEB: net *total* (surface) upwelling longwave radiation to atmosphere [W/m2]
  REAL, POINTER, DIMENSION(:) :: XP_LWNET_V       ! MEB: net vegetation canopy longwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XP_LWNET_G       ! MEB: net ground longwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XP_LWNET_N       ! MEB: net snow longwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XP_LEVCV         ! MEB: total evapotranspiration from vegetation canopy overstory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XP_LESC          ! MEB: total snow sublimation from vegetation canopy overstory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XP_H_V_C         ! MEB: sensible heat flux from vegetation canopy overstory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XP_H_G_C         ! MEB: sensible heat flux from understory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XP_LETRGV        ! MEB: transpiration from understory vegetation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XP_LETRCV        ! MEB: transpiration from overstory canopy vegetation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XP_LERGV         ! MEB: interception evaporation from understory vegetation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XP_LERCV         ! MEB: interception evaporation from overstory canopy vegetation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XP_H_C_A         ! MEB: sensible heat flux from canopy air space to the atmosphere [W/m2] 
  REAL, POINTER, DIMENSION(:) :: XP_H_N_C         ! MEB: sensible heat flux from the snow on the ground [W/m2]
                                                !      NOTE total sensible heat flux to the atmosphere also possibly 
                                                !      includes a contribution from snow covering the canopy
  REAL, POINTER, DIMENSION(:) :: XP_LE_V_C        ! MEB: latent heat flux from vegetation canopy overstory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XP_LE_G_C        ! MEB: latent heat flux from understory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XP_LE_C_A        ! MEB: latent heat flux from canopy air space to the atmosphere [W/m2] 
                                                !      NOTE total latent heat flux to the atmosphere also possibly 
                                                !      includes a contribution from snow covering the canopy
  REAL, POINTER, DIMENSION(:) :: XP_LE_N_C        ! MEB: latent heat flux from the snow on the ground [W/m2]
                                                !      NOTE total latent heat flux from the snowpack
                                                !      possibly includes a contribution from snow covering the canopy
  REAL, POINTER, DIMENSION(:) :: XP_EVAP_N_C      ! MEB: Total evap from snow on the ground to canopy air space  [kg/m2/s]
  REAL, POINTER, DIMENSION(:) :: XP_EVAP_G_C      ! MEB: Total evap from ground to canopy air space [kg/m2/s]
  REAL, POINTER, DIMENSION(:) :: XP_SR_GN         ! MEB: total snow reacing the ground snow [kg/m2/s]
  REAL, POINTER, DIMENSION(:) :: XP_MELTCV        ! MEB: snow melt rate from the overstory snow reservoir [kg/m2/s]
  REAL, POINTER, DIMENSION(:) :: XP_FRZCV         ! MEB: snow refreeze rate from the overstory snow reservoir [kg/m2/s]
  REAL, POINTER, DIMENSION(:) :: XP_SWDOWN_GN     ! MEB: total shortwave radiation transmitted through the canopy
                                                !      reaching the snowpack/ground understory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XP_LWDOWN_GN     ! MEB: total shortwave radiation transmitted through and emitted by the canopy
!                                               !      reaching the snowpack/ground understory (explicit part) [W/m2]
!
!------------------------------------------------------------------------------
!
  REAL, POINTER, DIMENSION(:) :: XP_DWG         ! liquid soil moisture time tendencies  (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XP_DWGI        ! solid soil moisture time tendencies   (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XP_DWR         ! canopy water time tendencies          (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XP_DSWE        ! snow water equivalent time tendencies (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XP_WATBUD      ! ISBA water budget                     (kg/m2/s)

END TYPE PACK_DIAG_ISBA_t
!
!-------------------------------------------------------------------------------
!
TYPE(PACK_DIAG_ISBA_t), ALLOCATABLE, TARGET, SAVE :: PACK_DIAG_ISBA_MODEL(:)

TYPE(PACK_DIAG_ISBA_t), POINTER :: PACK_DIAG_ISBA => NULL()
!$OMP THREADPRIVATE(PACK_DIAG_ISBA)

CONTAINS

SUBROUTINE PACK_DIAG_ISBA_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
!
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_PACK_DIAG_ISBA_N:PACK_DIAG_ISBA_GOTO_MODEL',0,ZHOOK_HANDLE)

PACK_DIAG_ISBA => PACK_DIAG_ISBA_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_PACK_DIAG_ISBA_N:PACK_DIAG_ISBA_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE PACK_DIAG_ISBA_GOTO_MODEL

SUBROUTINE PACK_DIAG_ISBA_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_PACK_DIAG_ISBA_N:PACK_DIAG_ISBA_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(PACK_DIAG_ISBA_MODEL(KMODEL))
PACK_DIAG_ISBA => PACK_DIAG_ISBA_MODEL(KMODEL)
DO J=1,KMODEL
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XBLOCK_SIMPLE)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XBLOCK_GROUND)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XBLOCK_SNOW)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XBLOCK_KSW)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XBLOCK_ABC)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XBLOCK_0)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XBLOCK_00)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_RNSNOW)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_HSNOW)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_HPSNOW)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_GFLUXSNOW)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_USTARSNOW)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_GRNDFLUX)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_SRSFC)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_RRSFC)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LESL)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_CDSNOW)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_CHSNOW)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_SNOWTEMP)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_SNOWLIQ)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_SNOWDZ)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_SNOWHMASS)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_RN_ISBA)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_H_ISBA)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LEG_ISBA)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LEGI_ISBA)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LEV_ISBA)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LETR_ISBA)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_USTAR_ISBA)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LER_ISBA)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LE_ISBA)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LEI_ISBA)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_GFLUX_ISBA)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_MELTADV)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_CH)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_CE)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_CD)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_CDN)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_RI)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_HU)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_HUG)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_HV)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_ALBT)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_RN)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_H)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LEG)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LEGI)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LEV)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LES)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LER)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LETR)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_EVAP)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_SUBL)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_SNDRIFT)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LEI)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_GFLUX)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_RESTORE)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_DRAIN)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_QSB)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_RUNOFF)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_MELT)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_SNOWFREE_ALB)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_SNOWFREE_ALB_VEG)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_SNOWFREE_ALB_SOIL)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_Z0_WITH_SNOW)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_Z0H_WITH_SNOW)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_Z0EFF)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_IACAN)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_CG)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_C1)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_C2)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_WGEQ)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_CT)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_RS)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_TS)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_TSRAD)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_T2M)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_Q2M)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_HU2M)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_ZON10M)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_MER10M)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_QS)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_SWI)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_TSWI)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_TWSNOW)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_TDSNOW)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_SWD)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_SWU)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_SWBD)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_SWBU)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LWD)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LWU)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_FMU)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_FMV)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_HORT)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_DRIP)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_IFLOOD)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_PFLOOD)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LE_FLOOD)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LEI_FLOOD)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_ICEFLUX)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_RRVEG)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_IRRIG_FLUX)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_GPP)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_RESP_AUTO)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_RESP_ECO)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_FAPAR)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_FAPIR)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_FAPAR_BS)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_FAPIR_BS)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_SWUP)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_SWNET_V)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_SWNET_G)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_SWNET_N)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_SWNET_NS)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LWUP)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LWNET_V)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LWNET_G)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LWNET_N)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LEVCV)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LESC)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_H_V_C)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_H_G_C)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LETRGV)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LETRCV)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LERGV)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LERCV)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_H_C_A)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_H_N_C)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LE_V_C)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LE_G_C)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LE_C_A)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LE_N_C)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_EVAP_N_C)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_EVAP_G_C)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_SR_GN)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_MELTCV)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_FRZCV)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_SWDOWN_GN)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_LWDOWN_GN)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_DWG)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_DWGI)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_DWR)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_DSWE)
  NULLIFY(PACK_DIAG_ISBA_MODEL(J)%XP_WATBUD)
ENDDO
IF (LHOOK) CALL DR_HOOK("MODD_PACK_DIAG_ISBA_N:PACK_DIAG_ISBA_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE PACK_DIAG_ISBA_ALLOC

SUBROUTINE PACK_DIAG_ISBA_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_PACK_DIAG_ISBA_N:PACK_DIAG_ISBA_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(PACK_DIAG_ISBA_MODEL)) DEALLOCATE(PACK_DIAG_ISBA_MODEL)
IF (ASSOCIATED(PACK_DIAG_ISBA)) NULLIFY(PACK_DIAG_ISBA)
IF (LHOOK) CALL DR_HOOK("MODD_PACK_DIAG_ISBA_N:PACK_DIAG_ISBA_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE PACK_DIAG_ISBA_DEALLO

END MODULE MODD_PACK_DIAG_ISBA
