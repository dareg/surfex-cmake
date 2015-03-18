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
!
INTEGER, POINTER :: NSIZE_SIMPLE
!$OMP THREADPRIVATE(NSIZE_SIMPLE)
INTEGER, POINTER :: NSIZE_GROUND
!$OMP THREADPRIVATE(NSIZE_GROUND)
INTEGER, POINTER :: NSIZE_SNOW
!$OMP THREADPRIVATE(NSIZE_SNOW)
INTEGER, POINTER :: NSIZE_KSW
!$OMP THREADPRIVATE(NSIZE_KSW)
INTEGER, POINTER :: NSIZE_ABC
!$OMP THREADPRIVATE(NSIZE_ABC)
INTEGER, POINTER :: NSIZE_0
!$OMP THREADPRIVATE(NSIZE_0)
INTEGER, POINTER :: NSIZE_00
!$OMP THREADPRIVATE(NSIZE_00)
REAL, POINTER, DIMENSION(:,:) :: XBLOCK_SIMPLE
!$OMP THREADPRIVATE(XBLOCK_SIMPLE)
REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_GROUND
!$OMP THREADPRIVATE(XBLOCK_GROUND)
REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_SNOW
!$OMP THREADPRIVATE(XBLOCK_SNOW)
REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_KSW
!$OMP THREADPRIVATE(XBLOCK_KSW)
REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_ABC
!$OMP THREADPRIVATE(XBLOCK_ABC)
REAL, POINTER, DIMENSION(:,:) :: XBLOCK_0
!$OMP THREADPRIVATE(XBLOCK_0)
REAL, POINTER, DIMENSION(:,:,:) :: XBLOCK_00
!$OMP THREADPRIVATE(XBLOCK_00)
!
REAL, POINTER, DIMENSION(:) :: XP_RNSNOW    
!$OMP THREADPRIVATE(XP_RNSNOW)
REAL, POINTER, DIMENSION(:) :: XP_HSNOW     
!$OMP THREADPRIVATE(XP_HSNOW)
REAL, POINTER, DIMENSION(:) :: XP_HPSNOW    
!$OMP THREADPRIVATE(XP_HPSNOW)
REAL, POINTER, DIMENSION(:) :: XP_GFLUXSNOW
!$OMP THREADPRIVATE(XP_GFLUXSNOW)
REAL, POINTER, DIMENSION(:) :: XP_USTARSNOW
!$OMP THREADPRIVATE(XP_USTARSNOW)
REAL, POINTER, DIMENSION(:) :: XP_GRNDFLUX
!$OMP THREADPRIVATE(XP_GRNDFLUX)
REAL, POINTER, DIMENSION(:) :: XP_SRSFC
!$OMP THREADPRIVATE(XP_SRSFC)
REAL, POINTER, DIMENSION(:) :: XP_RRSFC
!$OMP THREADPRIVATE(XP_RRSFC)
REAL, POINTER, DIMENSION(:) :: XP_LESL
!$OMP THREADPRIVATE(XP_LESL)
REAL, POINTER, DIMENSION(:) :: XP_CDSNOW
!$OMP THREADPRIVATE(XP_CDSNOW)
REAL, POINTER, DIMENSION(:) :: XP_CHSNOW
!$OMP THREADPRIVATE(XP_CHSNOW)
REAL, POINTER, DIMENSION(:,:)::XP_SNOWTEMP
!$OMP THREADPRIVATE(XP_SNOWTEMP)
REAL, POINTER, DIMENSION(:,:)::XP_SNOWLIQ
!$OMP THREADPRIVATE(XP_SNOWLIQ)
REAL, POINTER, DIMENSION(:,:)::XP_SNOWDZ
!$OMP THREADPRIVATE(XP_SNOWDZ)
REAL, POINTER, DIMENSION(:) :: XP_SNOWHMASS
!$OMP THREADPRIVATE(XP_SNOWHMASS)
!
REAL, POINTER, DIMENSION(:) :: XP_RN_ISBA
!$OMP THREADPRIVATE(XP_RN_ISBA)
REAL, POINTER, DIMENSION(:) :: XP_H_ISBA
!$OMP THREADPRIVATE(XP_H_ISBA)
REAL, POINTER, DIMENSION(:) :: XP_LEG_ISBA
!$OMP THREADPRIVATE(XP_LEG_ISBA)
REAL, POINTER, DIMENSION(:) :: XP_LEGI_ISBA 
!$OMP THREADPRIVATE(XP_LEGI_ISBA)
REAL, POINTER, DIMENSION(:) :: XP_LEV_ISBA  
!$OMP THREADPRIVATE(XP_LEV_ISBA)
REAL, POINTER, DIMENSION(:) :: XP_LETR_ISBA 
!$OMP THREADPRIVATE(XP_LETR_ISBA)
REAL, POINTER, DIMENSION(:) :: XP_USTAR_ISBA
!$OMP THREADPRIVATE(XP_USTAR_ISBA)
REAL, POINTER, DIMENSION(:) :: XP_LER_ISBA  
!$OMP THREADPRIVATE(XP_LER_ISBA)
REAL, POINTER, DIMENSION(:) :: XP_LE_ISBA   
!$OMP THREADPRIVATE(XP_LE_ISBA)
REAL, POINTER, DIMENSION(:) :: XP_LEI_ISBA  
!$OMP THREADPRIVATE(XP_LEI_ISBA)
REAL, POINTER, DIMENSION(:) :: XP_GFLUX_ISBA
!$OMP THREADPRIVATE(XP_GFLUX_ISBA)
REAL, POINTER, DIMENSION(:) :: XP_MELTADV   
!$OMP THREADPRIVATE(XP_MELTADV)
REAL, POINTER, DIMENSION(:) :: XP_CH        
!$OMP THREADPRIVATE(XP_CH)
REAL, POINTER, DIMENSION(:) :: XP_CE        
!$OMP THREADPRIVATE(XP_CE)
REAL, POINTER, DIMENSION(:) :: XP_CD        
!$OMP THREADPRIVATE(XP_CD)
REAL, POINTER, DIMENSION(:) :: XP_CDN       
!$OMP THREADPRIVATE(XP_CDN)
REAL, POINTER, DIMENSION(:) :: XP_RI        
!$OMP THREADPRIVATE(XP_RI)
REAL, POINTER, DIMENSION(:) :: XP_HU        
!$OMP THREADPRIVATE(XP_HU)
REAL, POINTER, DIMENSION(:) :: XP_HUG       
!$OMP THREADPRIVATE(XP_HUG)
REAL, POINTER, DIMENSION(:) :: XP_HV        
!$OMP THREADPRIVATE(XP_HV)
!
REAL, POINTER, DIMENSION(:) :: XP_ALBT      
!$OMP THREADPRIVATE(XP_ALBT)
!
REAL, POINTER, DIMENSION(:) :: XP_RN        
!$OMP THREADPRIVATE(XP_RN)
REAL, POINTER, DIMENSION(:) :: XP_H         
!$OMP THREADPRIVATE(XP_H)
REAL, POINTER, DIMENSION(:) :: XP_LEG       
!$OMP THREADPRIVATE(XP_LEG)
REAL, POINTER, DIMENSION(:) :: XP_LEGI      
!$OMP THREADPRIVATE(XP_LEGI)
REAL, POINTER, DIMENSION(:) :: XP_LEV       
!$OMP THREADPRIVATE(XP_LEV)
REAL, POINTER, DIMENSION(:) :: XP_LES       
!$OMP THREADPRIVATE(XP_LES)
REAL, POINTER, DIMENSION(:) :: XP_LER       
!$OMP THREADPRIVATE(XP_LER)

REAL, POINTER, DIMENSION(:) :: XP_LETR      
!$OMP THREADPRIVATE(XP_LETR)
REAL, POINTER, DIMENSION(:) :: XP_EVAP      
!$OMP THREADPRIVATE(XP_EVAP)
REAL, POINTER, DIMENSION(:) :: XP_SUBL      
!$OMP THREADPRIVATE(XP_SUBL)
REAL, POINTER, DIMENSION(:) :: XP_SNDRIFT   
!$OMP THREADPRIVATE(XP_SNDRIFT)
REAL, POINTER, DIMENSION(:) :: XP_LEI       
!$OMP THREADPRIVATE(XP_LEI)
REAL, POINTER, DIMENSION(:) :: XP_GFLUX     
!$OMP THREADPRIVATE(XP_GFLUX)
REAL, POINTER, DIMENSION(:) :: XP_RESTORE   
!$OMP THREADPRIVATE(XP_RESTORE)
REAL, POINTER, DIMENSION(:) :: XP_DRAIN     
!$OMP THREADPRIVATE(XP_DRAIN)
REAL, POINTER, DIMENSION(:) :: XP_QSB       
!$OMP THREADPRIVATE(XP_QSB)
REAL, POINTER, DIMENSION(:) :: XP_RUNOFF    
!$OMP THREADPRIVATE(XP_RUNOFF)
REAL, POINTER, DIMENSION(:) :: XP_MELT      
!$OMP THREADPRIVATE(XP_MELT)
REAL, POINTER, DIMENSION(:) :: XP_SNOWFREE_ALB 
!$OMP THREADPRIVATE(XP_SNOWFREE_ALB)
REAL, POINTER, DIMENSION(:) :: XP_SNOWFREE_ALB_VEG 
!$OMP THREADPRIVATE(XP_SNOWFREE_ALB_VEG)
REAL, POINTER, DIMENSION(:) :: XP_SNOWFREE_ALB_SOIL
!$OMP THREADPRIVATE(XP_SNOWFREE_ALB_SOIL)
REAL, POINTER, DIMENSION(:) :: XP_Z0_WITH_SNOW 
!$OMP THREADPRIVATE(XP_Z0_WITH_SNOW)
REAL, POINTER, DIMENSION(:) :: XP_Z0H_WITH_SNOW
!$OMP THREADPRIVATE(XP_Z0H_WITH_SNOW)
REAL, POINTER, DIMENSION(:) :: XP_Z0EFF     
!$OMP THREADPRIVATE(XP_Z0EFF)
!
REAL, POINTER, DIMENSION(:,:)::XP_IACAN     
!$OMP THREADPRIVATE(XP_IACAN)
!
REAL, POINTER, DIMENSION(:) :: XP_CG        
!$OMP THREADPRIVATE(XP_CG)
REAL, POINTER, DIMENSION(:) :: XP_C1        
!$OMP THREADPRIVATE(XP_C1)
REAL, POINTER, DIMENSION(:) :: XP_C2        
!$OMP THREADPRIVATE(XP_C2)
REAL, POINTER, DIMENSION(:) :: XP_WGEQ      
!$OMP THREADPRIVATE(XP_WGEQ)

REAL, POINTER, DIMENSION(:) :: XP_CT        
!$OMP THREADPRIVATE(XP_CT)
REAL, POINTER, DIMENSION(:) :: XP_RS        
!$OMP THREADPRIVATE(XP_RS)
!
!------------------------------------------------------------------------------
!
REAL, POINTER, DIMENSION(:) :: XP_TS        
!$OMP THREADPRIVATE(XP_TS)
REAL, POINTER, DIMENSION(:) :: XP_TSRAD     
!$OMP THREADPRIVATE(XP_TSRAD)
REAL, POINTER, DIMENSION(:) :: XP_T2M       
!$OMP THREADPRIVATE(XP_T2M)
REAL, POINTER, DIMENSION(:) :: XP_Q2M       
!$OMP THREADPRIVATE(XP_Q2M)
REAL, POINTER, DIMENSION(:) :: XP_HU2M      
!$OMP THREADPRIVATE(XP_HU2M)
REAL, POINTER, DIMENSION(:) :: XP_ZON10M    
!$OMP THREADPRIVATE(XP_ZON10M)
REAL, POINTER, DIMENSION(:) :: XP_MER10M    
!$OMP THREADPRIVATE(XP_MER10M)
!
!------------------------------------------------------------------------------
!
REAL, POINTER, DIMENSION(:)   :: XP_QS      
!$OMP THREADPRIVATE(XP_QS)
REAL, POINTER, DIMENSION(:,:) :: XP_SWI     
!$OMP THREADPRIVATE(XP_SWI)
REAL, POINTER, DIMENSION(:,:) :: XP_TSWI    
!$OMP THREADPRIVATE(XP_TSWI)
!
!------------------------------------------------------------------------------
!
REAL, POINTER, DIMENSION(:) :: XP_TWSNOW     
!$OMP THREADPRIVATE(XP_TWSNOW)
REAL, POINTER, DIMENSION(:) :: XP_TDSNOW     
!$OMP THREADPRIVATE(XP_TDSNOW)
!
!------------------------------------------------------------------------------
!
REAL, POINTER, DIMENSION(:) :: XP_SWD       
!$OMP THREADPRIVATE(XP_SWD)
REAL, POINTER, DIMENSION(:) :: XP_SWU       
!$OMP THREADPRIVATE(XP_SWU)
REAL, POINTER, DIMENSION(:,:) :: XP_SWBD    
!$OMP THREADPRIVATE(XP_SWBD)
REAL, POINTER, DIMENSION(:,:) :: XP_SWBU    
!$OMP THREADPRIVATE(XP_SWBU)
REAL, POINTER, DIMENSION(:) :: XP_LWD       
!$OMP THREADPRIVATE(XP_LWD)
REAL, POINTER, DIMENSION(:) :: XP_LWU       
!$OMP THREADPRIVATE(XP_LWU)
REAL, POINTER, DIMENSION(:) :: XP_FMU       
!$OMP THREADPRIVATE(XP_FMU)
REAL, POINTER, DIMENSION(:) :: XP_FMV       
!$OMP THREADPRIVATE(XP_FMV)
!
!------------------------------------------------------------------------------
!
REAL, POINTER, DIMENSION(:) :: XP_HORT      
!$OMP THREADPRIVATE(XP_HORT)
REAL, POINTER, DIMENSION(:) :: XP_DRIP      
!$OMP THREADPRIVATE(XP_DRIP)
REAL, POINTER, DIMENSION(:) :: XP_IFLOOD    
!$OMP THREADPRIVATE(XP_IFLOOD)
REAL, POINTER, DIMENSION(:) :: XP_PFLOOD    
!$OMP THREADPRIVATE(XP_PFLOOD)
REAL, POINTER, DIMENSION(:) :: XP_LE_FLOOD  
!$OMP THREADPRIVATE(XP_LE_FLOOD)
REAL, POINTER, DIMENSION(:) :: XP_LEI_FLOOD 
!$OMP THREADPRIVATE(XP_LEI_FLOOD)
REAL, POINTER, DIMENSION(:) :: XP_ICEFLUX
!$OMP THREADPRIVATE(XP_ICEFLUX)
REAL, POINTER, DIMENSION(:) :: XP_RRVEG     
!$OMP THREADPRIVATE(XP_RRVEG)
REAL, POINTER, DIMENSION(:) :: XP_IRRIG_FLUX
!$OMP THREADPRIVATE(XP_IRRIG_FLUX)
!
!------------------------------------------------------------------------------
!
REAL, POINTER, DIMENSION(:) :: XP_GPP         
!$OMP THREADPRIVATE(XP_GPP)
REAL, POINTER, DIMENSION(:) :: XP_RESP_AUTO   
!$OMP THREADPRIVATE(XP_RESP_AUTO)
REAL, POINTER, DIMENSION(:) :: XP_RESP_ECO    
!$OMP THREADPRIVATE(XP_RESP_ECO)
REAL, POINTER, DIMENSION(:) :: XP_FAPAR       
!$OMP THREADPRIVATE(XP_FAPAR)
REAL, POINTER, DIMENSION(:) :: XP_FAPIR       
!$OMP THREADPRIVATE(XP_FAPIR)
REAL, POINTER, DIMENSION(:) :: XP_FAPAR_BS    
!$OMP THREADPRIVATE(XP_FAPAR_BS)
REAL, POINTER, DIMENSION(:) :: XP_FAPIR_BS    
!$OMP THREADPRIVATE(XP_FAPIR_BS)
!
!------------------------------------------------------------------------------
!
!* diagnostic variables for multi-energy balance (MEB)
!  ---------------------------------------------------
!
REAL, POINTER, DIMENSION(:) :: XP_SWUP          
!$OMP THREADPRIVATE(XP_SWUP)
REAL, POINTER, DIMENSION(:) :: XP_SWNET_V       
!$OMP THREADPRIVATE(XP_SWNET_V)
REAL, POINTER, DIMENSION(:) :: XP_SWNET_G       
!$OMP THREADPRIVATE(XP_SWNET_G)
REAL, POINTER, DIMENSION(:) :: XP_SWNET_N       
!$OMP THREADPRIVATE(XP_SWNET_N)
REAL, POINTER, DIMENSION(:) :: XP_SWNET_NS      
                                                
!$OMP THREADPRIVATE(XP_SWNET_NS)
REAL, POINTER, DIMENSION(:) :: XP_LWUP          
!$OMP THREADPRIVATE(XP_LWUP)
REAL, POINTER, DIMENSION(:) :: XP_LWNET_V       
!$OMP THREADPRIVATE(XP_LWNET_V)
REAL, POINTER, DIMENSION(:) :: XP_LWNET_G       
!$OMP THREADPRIVATE(XP_LWNET_G)
REAL, POINTER, DIMENSION(:) :: XP_LWNET_N       
!$OMP THREADPRIVATE(XP_LWNET_N)
REAL, POINTER, DIMENSION(:) :: XP_LEVCV         
!$OMP THREADPRIVATE(XP_LEVCV)
REAL, POINTER, DIMENSION(:) :: XP_LESC          
!$OMP THREADPRIVATE(XP_LESC)
REAL, POINTER, DIMENSION(:) :: XP_H_V_C         
!$OMP THREADPRIVATE(XP_H_V_C)
REAL, POINTER, DIMENSION(:) :: XP_H_G_C         
!$OMP THREADPRIVATE(XP_H_G_C)
REAL, POINTER, DIMENSION(:) :: XP_LETRGV        
!$OMP THREADPRIVATE(XP_LETRGV)
REAL, POINTER, DIMENSION(:) :: XP_LETRCV        
!$OMP THREADPRIVATE(XP_LETRCV)
REAL, POINTER, DIMENSION(:) :: XP_LERGV         
!$OMP THREADPRIVATE(XP_LERGV)
REAL, POINTER, DIMENSION(:) :: XP_LERCV         
!$OMP THREADPRIVATE(XP_LERCV)
REAL, POINTER, DIMENSION(:) :: XP_H_C_A         
!$OMP THREADPRIVATE(XP_H_C_A)
REAL, POINTER, DIMENSION(:) :: XP_H_N_C         
                                                
                                                
!$OMP THREADPRIVATE(XP_H_N_C)
REAL, POINTER, DIMENSION(:) :: XP_LE_V_C        
!$OMP THREADPRIVATE(XP_LE_V_C)
REAL, POINTER, DIMENSION(:) :: XP_LE_G_C        
!$OMP THREADPRIVATE(XP_LE_G_C)
REAL, POINTER, DIMENSION(:) :: XP_LE_C_A        
                                                
                                                
!$OMP THREADPRIVATE(XP_LE_C_A)
REAL, POINTER, DIMENSION(:) :: XP_LE_N_C        
                                                
                                                
!$OMP THREADPRIVATE(XP_LE_N_C)
REAL, POINTER, DIMENSION(:) :: XP_EVAP_N_C      
!$OMP THREADPRIVATE(XP_EVAP_N_C)
REAL, POINTER, DIMENSION(:) :: XP_EVAP_G_C      
!$OMP THREADPRIVATE(XP_EVAP_G_C)
REAL, POINTER, DIMENSION(:) :: XP_SR_GN         
!$OMP THREADPRIVATE(XP_SR_GN)
REAL, POINTER, DIMENSION(:) :: XP_MELTCV        
!$OMP THREADPRIVATE(XP_MELTCV)
REAL, POINTER, DIMENSION(:) :: XP_FRZCV         
!$OMP THREADPRIVATE(XP_FRZCV)
REAL, POINTER, DIMENSION(:) :: XP_SWDOWN_GN     
                                                
!$OMP THREADPRIVATE(XP_SWDOWN_GN)
REAL, POINTER, DIMENSION(:) :: XP_LWDOWN_GN     

!$OMP THREADPRIVATE(XP_LWDOWN_GN)
!
!------------------------------------------------------------------------------
!
REAL, POINTER, DIMENSION(:) :: XP_DWG         
!$OMP THREADPRIVATE(XP_DWG)
REAL, POINTER, DIMENSION(:) :: XP_DWGI        
!$OMP THREADPRIVATE(XP_DWGI)
REAL, POINTER, DIMENSION(:) :: XP_DWR         
!$OMP THREADPRIVATE(XP_DWR)
REAL, POINTER, DIMENSION(:) :: XP_DSWE        
!$OMP THREADPRIVATE(XP_DSWE)
REAL, POINTER, DIMENSION(:) :: XP_WATBUD      
!$OMP THREADPRIVATE(XP_WATBUD)
!
!------------------------------------------------------------------------------
!
CONTAINS

SUBROUTINE PACK_DIAG_ISBA_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
!
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Save current state for allocated arrays
IF (LKFROM) THEN
  PACK_DIAG_ISBA_MODEL(KFROM)%XBLOCK_SIMPLE=>XBLOCK_SIMPLE
  PACK_DIAG_ISBA_MODEL(KFROM)%XBLOCK_GROUND=>XBLOCK_GROUND
  PACK_DIAG_ISBA_MODEL(KFROM)%XBLOCK_SNOW=>XBLOCK_SNOW
  PACK_DIAG_ISBA_MODEL(KFROM)%XBLOCK_KSW=>XBLOCK_KSW
  PACK_DIAG_ISBA_MODEL(KFROM)%XBLOCK_ABC=>XBLOCK_ABC
  PACK_DIAG_ISBA_MODEL(KFROM)%XBLOCK_0=>XBLOCK_0
  PACK_DIAG_ISBA_MODEL(KFROM)%XBLOCK_00=>XBLOCK_00
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_RNSNOW=>XP_RNSNOW
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_HSNOW=>XP_HSNOW
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_HPSNOW=>XP_HPSNOW
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_GFLUXSNOW=>XP_GFLUXSNOW
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_USTARSNOW=>XP_USTARSNOW
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_GRNDFLUX=>XP_GRNDFLUX
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_SRSFC=>XP_SRSFC
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_RRSFC=>XP_RRSFC
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LESL=>XP_LESL
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_CDSNOW=>XP_CDSNOW
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_CHSNOW=>XP_CHSNOW
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_SNOWTEMP=>XP_SNOWTEMP
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_SNOWLIQ=>XP_SNOWLIQ
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_SNOWDZ=>XP_SNOWDZ
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_SNOWHMASS=>XP_SNOWHMASS
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_RN_ISBA=>XP_RN_ISBA
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_H_ISBA=>XP_H_ISBA
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LEG_ISBA=>XP_LEG_ISBA
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LEGI_ISBA=>XP_LEGI_ISBA
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LEV_ISBA=>XP_LEV_ISBA
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LETR_ISBA=>XP_LETR_ISBA
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_USTAR_ISBA=>XP_USTAR_ISBA
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LER_ISBA=>XP_LER_ISBA
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LE_ISBA=>XP_LE_ISBA
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LEI_ISBA=>XP_LEI_ISBA
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_GFLUX_ISBA=>XP_GFLUX_ISBA
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_MELTADV=>XP_MELTADV
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_CH=>XP_CH
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_CE=>XP_CE
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_CD=>XP_CD
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_CDN=>XP_CDN
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_RI=>XP_RI
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_HU=>XP_HU
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_HUG=>XP_HUG
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_HV=>XP_HV
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_ALBT=>XP_ALBT
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_RN=>XP_RN
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_H=>XP_H
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LEG=>XP_LEG
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LEGI=>XP_LEGI
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LEV=>XP_LEV
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LES=>XP_LES
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LER=>XP_LER
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LETR=>XP_LETR
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_EVAP=>XP_EVAP
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_SUBL=>XP_SUBL
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_SNDRIFT=>XP_SNDRIFT
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LEI=>XP_LEI
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_GFLUX=>XP_GFLUX
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_RESTORE=>XP_RESTORE
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_DRAIN=>XP_DRAIN
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_QSB=>XP_QSB
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_RUNOFF=>XP_RUNOFF
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_MELT=>XP_MELT
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_SNOWFREE_ALB=>XP_SNOWFREE_ALB
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_SNOWFREE_ALB_VEG=>XP_SNOWFREE_ALB_VEG
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_SNOWFREE_ALB_SOIL=>XP_SNOWFREE_ALB_SOIL
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_Z0_WITH_SNOW=>XP_Z0_WITH_SNOW
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_Z0H_WITH_SNOW=>XP_Z0H_WITH_SNOW
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_Z0EFF=>XP_Z0EFF
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_IACAN=>XP_IACAN
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_CG=>XP_CG
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_C1=>XP_C1
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_C2=>XP_C2
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_WGEQ=>XP_WGEQ
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_CT=>XP_CT
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_RS=>XP_RS
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_TS=>XP_TS
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_TSRAD=>XP_TSRAD
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_T2M=>XP_T2M
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_Q2M=>XP_Q2M
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_HU2M=>XP_HU2M
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_ZON10M=>XP_ZON10M
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_MER10M=>XP_MER10M
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_QS=>XP_QS
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_SWI=>XP_SWI
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_TSWI=>XP_TSWI
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_TWSNOW=>XP_TWSNOW
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_TDSNOW=>XP_TDSNOW
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_SWD=>XP_SWD
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_SWU=>XP_SWU
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_SWBD=>XP_SWBD
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_SWBU=>XP_SWBU
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LWD=>XP_LWD
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LWU=>XP_LWU
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_FMU=>XP_FMU
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_FMV=>XP_FMV
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_HORT=>XP_HORT
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_DRIP=>XP_DRIP
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_IFLOOD=>XP_IFLOOD
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_PFLOOD=>XP_PFLOOD
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LE_FLOOD=>XP_LE_FLOOD
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LEI_FLOOD=>XP_LEI_FLOOD
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_ICEFLUX=>XP_ICEFLUX
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_RRVEG=>XP_RRVEG
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_IRRIG_FLUX=>XP_IRRIG_FLUX
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_GPP=>XP_GPP
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_RESP_AUTO=>XP_RESP_AUTO
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_RESP_ECO=>XP_RESP_ECO
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_FAPAR=>XP_FAPAR
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_FAPIR=>XP_FAPIR
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_FAPAR_BS=>XP_FAPAR_BS
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_FAPIR_BS=>XP_FAPIR_BS
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_SWUP=>XP_SWUP
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_SWNET_V=>XP_SWNET_V
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_SWNET_G=>XP_SWNET_G
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_SWNET_N=>XP_SWNET_N
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_SWNET_NS=>XP_SWNET_NS
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LWUP=>XP_LWUP
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LWNET_V=>XP_LWNET_V
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LWNET_G=>XP_LWNET_G
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LWNET_N=>XP_LWNET_N
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LEVCV=>XP_LEVCV
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LESC=>XP_LESC
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_H_V_C=>XP_H_V_C
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_H_G_C=>XP_H_G_C
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LETRGV=>XP_LETRGV
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LETRCV=>XP_LETRCV
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LERGV=>XP_LERGV
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LERCV=>XP_LERCV
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_H_C_A=>XP_H_C_A
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_H_N_C=>XP_H_N_C
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LE_V_C=>XP_LE_V_C
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LE_G_C=>XP_LE_G_C
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LE_C_A=>XP_LE_C_A
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LE_N_C=>XP_LE_N_C
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_EVAP_N_C=>XP_EVAP_N_C
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_EVAP_G_C=>XP_EVAP_G_C
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_SR_GN=>XP_SR_GN
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_MELTCV=>XP_MELTCV
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_FRZCV=>XP_FRZCV
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_SWDOWN_GN=>XP_SWDOWN_GN
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_LWDOWN_GN=>XP_LWDOWN_GN
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_DWG=>XP_DWG
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_DWGI=>XP_DWGI
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_DWR=>XP_DWR
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_DSWE=>XP_DSWE
  PACK_DIAG_ISBA_MODEL(KFROM)%XP_WATBUD=>XP_WATBUD
ENDIF
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_PACK_DIAG_ISBA_N:PACK_DIAG_ISBA_GOTO_MODEL',0,ZHOOK_HANDLE)
NSIZE_SIMPLE=>PACK_DIAG_ISBA_MODEL(KTO)%NSIZE_SIMPLE
NSIZE_GROUND=>PACK_DIAG_ISBA_MODEL(KTO)%NSIZE_GROUND
NSIZE_SNOW=>PACK_DIAG_ISBA_MODEL(KTO)%NSIZE_SNOW
NSIZE_KSW=>PACK_DIAG_ISBA_MODEL(KTO)%NSIZE_KSW
NSIZE_ABC=>PACK_DIAG_ISBA_MODEL(KTO)%NSIZE_ABC
NSIZE_0=>PACK_DIAG_ISBA_MODEL(KTO)%NSIZE_0
NSIZE_00=>PACK_DIAG_ISBA_MODEL(KTO)%NSIZE_00
XBLOCK_SIMPLE=>PACK_DIAG_ISBA_MODEL(KTO)%XBLOCK_SIMPLE
XBLOCK_GROUND=>PACK_DIAG_ISBA_MODEL(KTO)%XBLOCK_GROUND
XBLOCK_SNOW=>PACK_DIAG_ISBA_MODEL(KTO)%XBLOCK_SNOW
XBLOCK_KSW=>PACK_DIAG_ISBA_MODEL(KTO)%XBLOCK_KSW
XBLOCK_ABC=>PACK_DIAG_ISBA_MODEL(KTO)%XBLOCK_ABC
XBLOCK_0=>PACK_DIAG_ISBA_MODEL(KTO)%XBLOCK_0
XBLOCK_00=>PACK_DIAG_ISBA_MODEL(KTO)%XBLOCK_00
XP_RNSNOW=>PACK_DIAG_ISBA_MODEL(KTO)%XP_RNSNOW
XP_HSNOW=>PACK_DIAG_ISBA_MODEL(KTO)%XP_HSNOW
XP_HPSNOW=>PACK_DIAG_ISBA_MODEL(KTO)%XP_HPSNOW
XP_GFLUXSNOW=>PACK_DIAG_ISBA_MODEL(KTO)%XP_GFLUXSNOW
XP_USTARSNOW=>PACK_DIAG_ISBA_MODEL(KTO)%XP_USTARSNOW
XP_GRNDFLUX=>PACK_DIAG_ISBA_MODEL(KTO)%XP_GRNDFLUX
XP_SRSFC=>PACK_DIAG_ISBA_MODEL(KTO)%XP_SRSFC
XP_RRSFC=>PACK_DIAG_ISBA_MODEL(KTO)%XP_RRSFC
XP_LESL=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LESL
XP_CDSNOW=>PACK_DIAG_ISBA_MODEL(KTO)%XP_CDSNOW
XP_CHSNOW=>PACK_DIAG_ISBA_MODEL(KTO)%XP_CHSNOW
XP_SNOWTEMP=>PACK_DIAG_ISBA_MODEL(KTO)%XP_SNOWTEMP
XP_SNOWLIQ=>PACK_DIAG_ISBA_MODEL(KTO)%XP_SNOWLIQ
XP_SNOWDZ=>PACK_DIAG_ISBA_MODEL(KTO)%XP_SNOWDZ
XP_SNOWHMASS=>PACK_DIAG_ISBA_MODEL(KTO)%XP_SNOWHMASS
XP_RN_ISBA=>PACK_DIAG_ISBA_MODEL(KTO)%XP_RN_ISBA
XP_H_ISBA=>PACK_DIAG_ISBA_MODEL(KTO)%XP_H_ISBA
XP_LEG_ISBA=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LEG_ISBA
XP_LEGI_ISBA=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LEGI_ISBA
XP_LEV_ISBA=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LEV_ISBA
XP_LETR_ISBA=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LETR_ISBA
XP_USTAR_ISBA=>PACK_DIAG_ISBA_MODEL(KTO)%XP_USTAR_ISBA
XP_LER_ISBA=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LER_ISBA
XP_LE_ISBA=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LE_ISBA
XP_LEI_ISBA=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LEI_ISBA
XP_GFLUX_ISBA=>PACK_DIAG_ISBA_MODEL(KTO)%XP_GFLUX_ISBA
XP_MELTADV=>PACK_DIAG_ISBA_MODEL(KTO)%XP_MELTADV
XP_CH=>PACK_DIAG_ISBA_MODEL(KTO)%XP_CH
XP_CE=>PACK_DIAG_ISBA_MODEL(KTO)%XP_CE
XP_CD=>PACK_DIAG_ISBA_MODEL(KTO)%XP_CD
XP_CDN=>PACK_DIAG_ISBA_MODEL(KTO)%XP_CDN
XP_RI=>PACK_DIAG_ISBA_MODEL(KTO)%XP_RI
XP_HU=>PACK_DIAG_ISBA_MODEL(KTO)%XP_HU
XP_HUG=>PACK_DIAG_ISBA_MODEL(KTO)%XP_HUG
XP_HV=>PACK_DIAG_ISBA_MODEL(KTO)%XP_HV
XP_ALBT=>PACK_DIAG_ISBA_MODEL(KTO)%XP_ALBT
XP_RN=>PACK_DIAG_ISBA_MODEL(KTO)%XP_RN
XP_H=>PACK_DIAG_ISBA_MODEL(KTO)%XP_H
XP_LEG=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LEG
XP_LEGI=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LEGI
XP_LEV=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LEV
XP_LES=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LES
XP_LER=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LER
XP_LETR=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LETR
XP_EVAP=>PACK_DIAG_ISBA_MODEL(KTO)%XP_EVAP
XP_SUBL=>PACK_DIAG_ISBA_MODEL(KTO)%XP_SUBL
XP_SNDRIFT=>PACK_DIAG_ISBA_MODEL(KTO)%XP_SNDRIFT
XP_LEI=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LEI
XP_GFLUX=>PACK_DIAG_ISBA_MODEL(KTO)%XP_GFLUX
XP_RESTORE=>PACK_DIAG_ISBA_MODEL(KTO)%XP_RESTORE
XP_DRAIN=>PACK_DIAG_ISBA_MODEL(KTO)%XP_DRAIN
XP_QSB=>PACK_DIAG_ISBA_MODEL(KTO)%XP_QSB
XP_RUNOFF=>PACK_DIAG_ISBA_MODEL(KTO)%XP_RUNOFF
XP_MELT=>PACK_DIAG_ISBA_MODEL(KTO)%XP_MELT
XP_SNOWFREE_ALB=>PACK_DIAG_ISBA_MODEL(KTO)%XP_SNOWFREE_ALB
XP_SNOWFREE_ALB_VEG=>PACK_DIAG_ISBA_MODEL(KTO)%XP_SNOWFREE_ALB_VEG
XP_SNOWFREE_ALB_SOIL=>PACK_DIAG_ISBA_MODEL(KTO)%XP_SNOWFREE_ALB_SOIL
XP_Z0_WITH_SNOW=>PACK_DIAG_ISBA_MODEL(KTO)%XP_Z0_WITH_SNOW
XP_Z0H_WITH_SNOW=>PACK_DIAG_ISBA_MODEL(KTO)%XP_Z0H_WITH_SNOW
XP_Z0EFF=>PACK_DIAG_ISBA_MODEL(KTO)%XP_Z0EFF
XP_IACAN=>PACK_DIAG_ISBA_MODEL(KTO)%XP_IACAN
XP_CG=>PACK_DIAG_ISBA_MODEL(KTO)%XP_CG
XP_C1=>PACK_DIAG_ISBA_MODEL(KTO)%XP_C1
XP_C2=>PACK_DIAG_ISBA_MODEL(KTO)%XP_C2
XP_WGEQ=>PACK_DIAG_ISBA_MODEL(KTO)%XP_WGEQ
XP_CT=>PACK_DIAG_ISBA_MODEL(KTO)%XP_CT
XP_RS=>PACK_DIAG_ISBA_MODEL(KTO)%XP_RS
XP_TS=>PACK_DIAG_ISBA_MODEL(KTO)%XP_TS
XP_TSRAD=>PACK_DIAG_ISBA_MODEL(KTO)%XP_TSRAD
XP_T2M=>PACK_DIAG_ISBA_MODEL(KTO)%XP_T2M
XP_Q2M=>PACK_DIAG_ISBA_MODEL(KTO)%XP_Q2M
XP_HU2M=>PACK_DIAG_ISBA_MODEL(KTO)%XP_HU2M
XP_ZON10M=>PACK_DIAG_ISBA_MODEL(KTO)%XP_ZON10M
XP_MER10M=>PACK_DIAG_ISBA_MODEL(KTO)%XP_MER10M
XP_QS=>PACK_DIAG_ISBA_MODEL(KTO)%XP_QS
XP_SWI=>PACK_DIAG_ISBA_MODEL(KTO)%XP_SWI
XP_TSWI=>PACK_DIAG_ISBA_MODEL(KTO)%XP_TSWI
XP_TWSNOW=>PACK_DIAG_ISBA_MODEL(KTO)%XP_TWSNOW
XP_TDSNOW=>PACK_DIAG_ISBA_MODEL(KTO)%XP_TDSNOW
XP_SWD=>PACK_DIAG_ISBA_MODEL(KTO)%XP_SWD
XP_SWU=>PACK_DIAG_ISBA_MODEL(KTO)%XP_SWU
XP_SWBD=>PACK_DIAG_ISBA_MODEL(KTO)%XP_SWBD
XP_SWBU=>PACK_DIAG_ISBA_MODEL(KTO)%XP_SWBU
XP_LWD=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LWD
XP_LWU=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LWU
XP_FMU=>PACK_DIAG_ISBA_MODEL(KTO)%XP_FMU
XP_FMV=>PACK_DIAG_ISBA_MODEL(KTO)%XP_FMV
XP_HORT=>PACK_DIAG_ISBA_MODEL(KTO)%XP_HORT
XP_DRIP=>PACK_DIAG_ISBA_MODEL(KTO)%XP_DRIP
XP_IFLOOD=>PACK_DIAG_ISBA_MODEL(KTO)%XP_IFLOOD
XP_PFLOOD=>PACK_DIAG_ISBA_MODEL(KTO)%XP_PFLOOD
XP_LE_FLOOD=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LE_FLOOD
XP_LEI_FLOOD=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LEI_FLOOD
XP_ICEFLUX=>PACK_DIAG_ISBA_MODEL(KTO)%XP_ICEFLUX
XP_RRVEG=>PACK_DIAG_ISBA_MODEL(KTO)%XP_RRVEG
XP_IRRIG_FLUX=>PACK_DIAG_ISBA_MODEL(KTO)%XP_IRRIG_FLUX
XP_GPP=>PACK_DIAG_ISBA_MODEL(KTO)%XP_GPP
XP_RESP_AUTO=>PACK_DIAG_ISBA_MODEL(KTO)%XP_RESP_AUTO
XP_RESP_ECO=>PACK_DIAG_ISBA_MODEL(KTO)%XP_RESP_ECO
XP_FAPAR=>PACK_DIAG_ISBA_MODEL(KTO)%XP_FAPAR
XP_FAPIR=>PACK_DIAG_ISBA_MODEL(KTO)%XP_FAPIR
XP_FAPAR_BS=>PACK_DIAG_ISBA_MODEL(KTO)%XP_FAPAR_BS
XP_FAPIR_BS=>PACK_DIAG_ISBA_MODEL(KTO)%XP_FAPIR_BS
XP_SWUP=>PACK_DIAG_ISBA_MODEL(KTO)%XP_SWUP
XP_SWNET_V=>PACK_DIAG_ISBA_MODEL(KTO)%XP_SWNET_V
XP_SWNET_G=>PACK_DIAG_ISBA_MODEL(KTO)%XP_SWNET_G
XP_SWNET_N=>PACK_DIAG_ISBA_MODEL(KTO)%XP_SWNET_N
XP_SWNET_NS=>PACK_DIAG_ISBA_MODEL(KTO)%XP_SWNET_NS
XP_LWUP=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LWUP
XP_LWNET_V=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LWNET_V
XP_LWNET_G=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LWNET_G
XP_LWNET_N=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LWNET_N
XP_LEVCV=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LEVCV
XP_LESC=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LESC
XP_H_V_C=>PACK_DIAG_ISBA_MODEL(KTO)%XP_H_V_C
XP_H_G_C=>PACK_DIAG_ISBA_MODEL(KTO)%XP_H_G_C
XP_LETRGV=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LETRGV
XP_LETRCV=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LETRCV
XP_LERGV=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LERGV
XP_LERCV=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LERCV
XP_H_C_A=>PACK_DIAG_ISBA_MODEL(KTO)%XP_H_C_A
XP_H_N_C=>PACK_DIAG_ISBA_MODEL(KTO)%XP_H_N_C
XP_LE_V_C=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LE_V_C
XP_LE_G_C=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LE_G_C
XP_LE_C_A=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LE_C_A
XP_LE_N_C=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LE_N_C
XP_EVAP_N_C=>PACK_DIAG_ISBA_MODEL(KTO)%XP_EVAP_N_C
XP_EVAP_G_C=>PACK_DIAG_ISBA_MODEL(KTO)%XP_EVAP_G_C
XP_SR_GN=>PACK_DIAG_ISBA_MODEL(KTO)%XP_SR_GN
XP_MELTCV=>PACK_DIAG_ISBA_MODEL(KTO)%XP_MELTCV
XP_FRZCV=>PACK_DIAG_ISBA_MODEL(KTO)%XP_FRZCV
XP_SWDOWN_GN=>PACK_DIAG_ISBA_MODEL(KTO)%XP_SWDOWN_GN
XP_LWDOWN_GN=>PACK_DIAG_ISBA_MODEL(KTO)%XP_LWDOWN_GN
XP_DWG=>PACK_DIAG_ISBA_MODEL(KTO)%XP_DWG
XP_DWGI=>PACK_DIAG_ISBA_MODEL(KTO)%XP_DWGI
XP_DWR=>PACK_DIAG_ISBA_MODEL(KTO)%XP_DWR
XP_DSWE=>PACK_DIAG_ISBA_MODEL(KTO)%XP_DSWE
XP_WATBUD=>PACK_DIAG_ISBA_MODEL(KTO)%XP_WATBUD
IF (LHOOK) CALL DR_HOOK('MODD_PACK_DIAG_ISBA_N:PACK_DIAG_ISBA_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE PACK_DIAG_ISBA_GOTO_MODEL

SUBROUTINE PACK_DIAG_ISBA_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_PACK_DIAG_ISBA_N:PACK_DIAG_ISBA_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(PACK_DIAG_ISBA_MODEL(KMODEL))
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
IF (LHOOK) CALL DR_HOOK("MODD_PACK_DIAG_ISBA_N:PACK_DIAG_ISBA_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE PACK_DIAG_ISBA_DEALLO

END MODULE MODD_PACK_DIAG_ISBA
