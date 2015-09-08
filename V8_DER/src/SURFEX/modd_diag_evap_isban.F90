!########################
MODULE MODD_DIAG_EVAP_ISBA_n
!########################
!
!!****  *MODD_DIAG_ISBA - declaration of packed surface parameters for ISBA scheme
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
!!      P. Le Moigne   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       07/11/03
!!      P. Samuelsson  04/2012   MEB
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE DIAG_EVAP_ISBA_t
!------------------------------------------------------------------------------
!
  LOGICAL :: LSURF_EVAP_BUDGET   ! flag for all terms of evaporation
  LOGICAL :: LSURF_BUDGETC       ! flag for surface cumulated energy budget
  LOGICAL :: LRESET_BUDGETC      ! flag for surface cumulated energy budget
  LOGICAL :: LWATER_BUDGET       ! flag for isba water budget including input  
                                 ! fluxes (rain and snow) and reservoir tendencies
!
!* variables for each patch
!
  REAL, POINTER, DIMENSION(:,:) :: XLEG          ! latent heat of evaporation over the ground   (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XLEGI         ! surface soil ice sublimation                 (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XLEV          ! latent heat of evaporation over vegetation   (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XLES          ! latent heat of sublimation over the snow     (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XLESL         ! latent heat of evaporation over the snow     (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XLER          ! evaporation from canopy water interception   (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XLETR         ! evapotranspiration of the vegetation         (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XEVAP         ! evapotranspiration                           (kg/m2/s)
  REAL, POINTER, DIMENSION(:,:) :: XSUBL         ! sublimation                                  (kg/m2/s)
  REAL, POINTER, DIMENSION(:,:) :: XSNDRIFT      ! blowing snow sublimation (ES or Crocus)      (kg/m2/s)
  REAL, POINTER, DIMENSION(:,:) :: XDRAIN        ! soil drainage flux                           (kg/m2/s)
  REAL, POINTER, DIMENSION(:,:) :: XQSB          ! lateral subsurface flux (dif option)         (kg/m2/s)
  REAL, POINTER, DIMENSION(:,:) :: XRUNOFF       ! sub-grid and supersaturation runoff          (kg/m2/s)
  REAL, POINTER, DIMENSION(:,:) :: XHORT         ! sub-grid Horton runoff from the SGH scheme   (kg/m2/s)
  REAL, POINTER, DIMENSION(:,:) :: XRRVEG        !  precipitation intercepted by the vegetation (kg/m2/s)
  REAL, POINTER, DIMENSION(:,:) :: XMELT         ! snow melt                                    (kg/m2/s)      
  REAL, POINTER, DIMENSION(:,:) :: XIFLOOD       ! Floodplains infiltration                     (kg/m2/s)      
  REAL, POINTER, DIMENSION(:,:) :: XPFLOOD       ! Precipitation intercepted by the floodplains (kg/m2/s)      
  REAL, POINTER, DIMENSION(:,:) :: XLE_FLOOD     ! Floodplains evapotration                     (W/m2)      
  REAL, POINTER, DIMENSION(:,:) :: XLEI_FLOOD    ! Floodplains evapotration                     (W/m2)      
  REAL, POINTER, DIMENSION(:,:) :: XDRIP         ! dripping from the vegetation reservoir       (kg/m2/s)
  REAL, POINTER, DIMENSION(:,:) :: XIRRIG_FLUX   ! irrigation rate (as soil input)              (kg/m2/s)
!  
  REAL, POINTER, DIMENSION(:,:) :: XGPP          ! Gross Primary Production                     (kgCO2/m2/s)
  REAL, POINTER, DIMENSION(:,:) :: XRESP_AUTO    ! Autotrophic respiration                      (kgCO2/m2/s)
  REAL, POINTER, DIMENSION(:,:) :: XRESP_ECO     ! Ecosystem respiration                        (kgCO2/m2/s)
!
  REAL, POINTER, DIMENSION(:,:) :: XLEVCV        ! MEB: total evapotranspiration from vegetation canopy overstory [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XLESC         ! MEB: total snow sublimation from vegetation canopy overstory [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XLETRGV       ! MEB: transpiration from understory vegetation [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XLETRCV       ! MEB: transpiration from overstory canopy vegetation [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XLERGV        ! MEB: interception evaporation from understory vegetation [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XLERCV        ! MEB: interception evaporation from overstory canopy vegetation [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XLE_V_C       ! MEB: latent heat flux from vegetation canopy overstory [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XLE_G_C       ! MEB: latent heat flux from understory [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XLE_C_A       ! MEB: latent heat flux from canopy air space to the atmosphere [W/m2] 
                                                 !      NOTE total latent heat flux to the atmosphere also possibly 
                                                 !      includes a contribution from snow covering the canopy
  REAL, POINTER, DIMENSION(:,:) :: XLE_N_C       ! MEB: latent heat flux from the snow on the ground [W/m2]
                                                 !      NOTE total latent heat flux from the snowpack
                                                 !      possibly includes a contribution from snow covering the canopy
  REAL, POINTER, DIMENSION(:,:) :: XSWNET_V      ! MEB: net vegetation canopy shortwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XSWNET_G      ! MEB: net ground shortwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XSWNET_N      ! MEB: net snow shortwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XSWNET_NS     ! MEB: net snow shortwave radiation for *surface* layer 
                                                 !     (i.e. net snow shortwave radiation less absorbed radiation) [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XLWNET_V      ! MEB: net vegetation canopy longwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XLWNET_G      ! MEB: net ground longwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XLWNET_N      ! MEB: net snow longwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XH_V_C        ! MEB: sensible heat flux from vegetation canopy overstory [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XH_G_C        ! MEB: sensible heat flux from understory [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XH_C_A        ! MEB: sensible heat flux from canopy air space to the atmosphere [W/m2] 
                                                 !      NOTE total sensible heat flux to the atmosphere also possibly 
                                                 !      includes a contribution from snow covering the canopy
  REAL, POINTER, DIMENSION(:,:) :: XH_N_C        ! MEB: sensible heat flux from the snow on the ground [W/m2]
                                                 !      NOTE total sensible heat flux from the snowpack
                                                 !      possibly includes a contribution from snow covering the canopy
  REAL, POINTER, DIMENSION(:,:) :: XSR_GN        ! MEB: snow unloading rate from the overstory reservoir [kg/m2/s]
  REAL, POINTER, DIMENSION(:,:) :: XMELTCV       ! MEB: snow melt rate from the overstory snow reservoir [kg/m2/s]
  REAL, POINTER, DIMENSION(:,:) :: XFRZCV        ! MEB: snow refreeze rate from the overstory snow reservoir [kg/m2/s]
  REAL, POINTER, DIMENSION(:,:) :: XSWDOWN_GN    ! MEB: total shortwave radiation transmitted through the canopy
                                                 !      reaching the snowpack/ground understory [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XLWDOWN_GN    ! MEB: total shortwave radiation transmitted through and emitted by the canopy
                                                 !      reaching the snowpack/ground understory (explicit part) [W/m2]
!
  REAL, POINTER, DIMENSION(:,:) :: XDWG          ! liquid soil moisture time tendencies         (kg/m2/s)
  REAL, POINTER, DIMENSION(:,:) :: XDWGI         ! solid soil moisture time tendencies          (kg/m2/s)
  REAL, POINTER, DIMENSION(:,:) :: XDWR          ! canopy water time tendencies                 (kg/m2/s)
  REAL, POINTER, DIMENSION(:,:) :: XDSWE         ! snow water equivalent time tendencies        (kg/m2/s)
  REAL, POINTER, DIMENSION(:,:) :: XWATBUD       ! ISBA water budget                            (kg/m2/s)
!
!* average variables
!
  REAL, POINTER, DIMENSION(:)   :: XAVG_LEG      ! latent heat of evaporation over the ground   (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LEGI     ! surface soil ice sublimation                 (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LEV      ! latent heat of evaporation over vegetation   (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LES      ! latent heat of sublimation over the snow     (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LESL     ! latent heat of evaporation over the snow     (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LER      ! evaporation from canopy water interception   (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LETR     ! evapotranspiration of the vegetation         (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_EVAP     ! evapotranspiration                           (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_SUBL     ! sublimation                                  (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_SNDRIFT  ! blowing snow sublimation (ES or Crocus)      (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_DRAIN    ! soil drainage flux                           (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_QSB      ! lateral subsurface flux (dif option)         (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_RUNOFF   ! sub-grid and supersaturation runoff          (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_HORT     ! sub-grid Horton runoff from the SGH scheme   (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_DRIP     ! dripping from the vegetation reservoir       (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_MELT     ! snow melt                                    (kg/m2/s)      
  REAL, POINTER, DIMENSION(:)   :: XAVG_IFLOOD   ! Floodplains infiltration                     (kg/m2/s)      
  REAL, POINTER, DIMENSION(:)   :: XAVG_PFLOOD   ! Precipitation intercepted by the floodplains (kg/m2/s)      
  REAL, POINTER, DIMENSION(:)   :: XAVG_LE_FLOOD ! Floodplains evapotration                     (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LEI_FLOOD! Floodplains evapotration                     (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_RRVEG    ! precipitation intercepted by the vegetation  (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_IRRIG_FLUX! irrigation rate (as soil input)             (kg/m2/s)
!  
  REAL, POINTER, DIMENSION(:)   :: XAVG_GPP      ! Gross Primary Production                     (kgCO2/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_RESP_AUTO! Autotrophic respiration                      (kgCO2/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_RESP_ECO ! Ecosystem respiration                        (kgCO2/m2/s)
!
  REAL, POINTER, DIMENSION(:) :: XAVG_LEVCV        ! MEB: total evapotranspiration from vegetation canopy overstory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_LESC         ! MEB: total snow sublimation from vegetation canopy overstory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_LETRGV       ! MEB: transpiration from understory vegetation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_LETRCV       ! MEB: transpiration from overstory canopy vegetation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_LERGV        ! MEB: interception evaporation from understory vegetation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_LERCV        ! MEB: interception evaporation from overstory canopy vegetation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_LE_V_C       ! MEB: latent heat flux from vegetation canopy overstory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_LE_G_C       ! MEB: latent heat flux from understory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_LE_C_A       ! MEB: latent heat flux from canopy air space to the atmosphere [W/m2] 
                                                 !      NOTE total latent heat flux to the atmosphere also possibly 
                                                 !      includes a contribution from snow covering the canopy
  REAL, POINTER, DIMENSION(:) :: XAVG_LE_N_C       ! MEB: latent heat flux from the snow on the ground [W/m2]
                                                 !      NOTE total latent heat flux from the snowpack
                                                 !      possibly includes a contribution from snow covering the canopy
  REAL, POINTER, DIMENSION(:) :: XAVG_SWNET_V      ! MEB: net vegetation canopy shortwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_SWNET_G      ! MEB: net ground shortwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_SWNET_N      ! MEB: net snow shortwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_SWNET_NS     ! MEB: net snow shortwave radiation for *surface* layer 
                                                 !     (i.e. net snow shortwave radiation less absorbed radiation) [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_LWNET_V      ! MEB: net vegetation canopy longwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_LWNET_G      ! MEB: net ground longwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_LWNET_N      ! MEB: net snow longwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_H_V_C        ! MEB: sensible heat flux from vegetation canopy overstory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_H_G_C        ! MEB: sensible heat flux from understory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_H_C_A        ! MEB: sensible heat flux from canopy air space to the atmosphere [W/m2] 
                                                 !      NOTE total sensible heat flux to the atmosphere also possibly 
                                                 !      includes a contribution from snow covering the canopy
  REAL, POINTER, DIMENSION(:) :: XAVG_H_N_C        ! MEB: sensible heat flux from the snow on the ground [W/m2]
                                                 !      NOTE total sensible heat flux from the snowpack
                                                 !      possibly includes a contribution from snow covering the canopy
  REAL, POINTER, DIMENSION(:) :: XAVG_SR_GN        ! MEB: snow unloading rate from the overstory reservoir [kg/m2/s]
  REAL, POINTER, DIMENSION(:) :: XAVG_MELTCV       ! MEB: snow melt rate from the overstory snow reservoir [kg/m2/s]
  REAL, POINTER, DIMENSION(:) :: XAVG_FRZCV        ! MEB: snow refreeze rate from the overstory snow reservoir [kg/m2/s]
  REAL, POINTER, DIMENSION(:) :: XAVG_SWDOWN_GN    ! MEB: total shortwave radiation transmitted through the canopy
                                                 !      reaching the snowpack/ground understory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_LWDOWN_GN    ! MEB: total shortwave radiation transmitted through and emitted by the canopy
                                                 !      reaching the snowpack/ground understory (explicit part) [W/m2]
!
  REAL, POINTER, DIMENSION(:)   :: XRAINFALL     ! input rainfall rate for LWATER_BUDGET        (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XSNOWFALL     ! input snowfall rate for LWATER_BUDGET        (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_DWG      ! liquid soil moisture time tendencies         (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_DWGI     ! solid soil moisture time tendencies          (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_DWR      ! canopy water time tendencies                 (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_DSWE     ! snow water equivalent time tendencies        (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_WATBUD   ! ISBA water budget                            (kg/m2/s)
!  
!* budget summation variables for each patch
!
  REAL, POINTER, DIMENSION(:,:) :: XRNC          ! net radiation at surface                     (J/m2)
  REAL, POINTER, DIMENSION(:,:) :: XHC           ! sensible heat flux                           (J/m2)
  REAL, POINTER, DIMENSION(:,:) :: XLEC          ! total latent heat flux                       (J/m2)
  REAL, POINTER, DIMENSION(:,:) :: XLEIC         ! sublimation latent heat flux                 (J/m2)
  REAL, POINTER, DIMENSION(:,:) :: XGFLUXC       ! net soil-vegetation flux                     (J/m2)
  REAL, POINTER, DIMENSION(:,:) :: XLEGC         ! latent heat of evaporation over the ground   (J/m2)
  REAL, POINTER, DIMENSION(:,:) :: XLEGIC        ! surface soil ice sublimation                 (J/m2)
  REAL, POINTER, DIMENSION(:,:) :: XLEVC         ! latent heat of evaporation over vegetation   (J/m2)
  REAL, POINTER, DIMENSION(:,:) :: XLESAC         ! latent heat of sublimation over the snow     (J/m2)
  REAL, POINTER, DIMENSION(:,:) :: XLESLC        ! latent heat of evaporation over the snow     (J/m2)
  REAL, POINTER, DIMENSION(:,:) :: XLERC         ! evaporation from canopy water interception   (J/m2)
  REAL, POINTER, DIMENSION(:,:) :: XLETRC        ! evapotranspiration of the vegetation         (J/m2)
  REAL, POINTER, DIMENSION(:,:) :: XEVAPC        ! evapotranspiration                           (kg/m2)
  REAL, POINTER, DIMENSION(:,:) :: XSUBLC        ! sublimation                                  (kg/m2)
  REAL, POINTER, DIMENSION(:,:) :: XSNDRIFTC     ! blowing snow sublimation (ES or Crocus)      (kg/m2)
  REAL, POINTER, DIMENSION(:,:) :: XDRAINC       ! soil drainage flux                           (kg/m2) 
  REAL, POINTER, DIMENSION(:,:) :: XQSBC         ! lateral subsurface flux (dif option)         (kg/m2) 
  REAL, POINTER, DIMENSION(:,:) :: XRUNOFFC      ! sub-grid and supersaturation runoff          (kg/m2)
  REAL, POINTER, DIMENSION(:,:) :: XHORTC        ! sub-grid Horton runoff from the SGH scheme   (kg/m2)
  REAL, POINTER, DIMENSION(:,:) :: XDRIPC        ! dripping from the vegetation reservoir       (kg/m2/s)
  REAL, POINTER, DIMENSION(:,:) :: XMELTC        ! snow melt                                    (kg/m2)      
  REAL, POINTER, DIMENSION(:,:) :: XIFLOODC      ! Floodplains infiltration                     (kg/m2)      
  REAL, POINTER, DIMENSION(:,:) :: XPFLOODC      ! Precipitation intercepted by the floodplains (kg/m2)      
  REAL, POINTER, DIMENSION(:,:) :: XLE_FLOODC    ! Floodplains evapotration                     (J/m2)  
  REAL, POINTER, DIMENSION(:,:) :: XLEI_FLOODC   ! Floodplains evapotration                     (J/m2)  
  REAL, POINTER, DIMENSION(:,:) :: XICEFLUXC     ! ice calving flux                             (kg/m2)
  REAL, POINTER, DIMENSION(:,:) :: XRRVEGC       ! precipitation intercepted by the vegetation  (kg/m2)
  REAL, POINTER, DIMENSION(:,:) :: XIRRIG_FLUXC  ! irrigation rate (as soil input)              (kg/m2)
!  
  REAL, POINTER, DIMENSION(:,:) :: XGPPC         ! Gross Primary Production                     (kgCO2/m2)
  REAL, POINTER, DIMENSION(:,:) :: XRESPC_AUTO   ! Autotrophic respiration                      (kgCO2/m2)
  REAL, POINTER, DIMENSION(:,:) :: XRESPC_ECO    ! Ecosystem respiration                        (kgCO2/m2)
!
  REAL, POINTER, DIMENSION(:,:) :: XLEVCVC        ! MEB: total evapotranspiration from vegetation canopy overstory [J/m2]
  REAL, POINTER, DIMENSION(:,:) :: XLESCC         ! MEB: total snow sublimation from vegetation canopy overstory [J/m2]
  REAL, POINTER, DIMENSION(:,:) :: XLETRGVC       ! MEB: transpiration from understory vegetation [J/m2]
  REAL, POINTER, DIMENSION(:,:) :: XLETRCVC       ! MEB: transpiration from overstory canopy vegetation [J/m2]
  REAL, POINTER, DIMENSION(:,:) :: XLERGVC        ! MEB: interception evaporation from understory vegetation [J/m2]
  REAL, POINTER, DIMENSION(:,:) :: XLERCVC        ! MEB: interception evaporation from overstory canopy vegetation [J/m2]
  REAL, POINTER, DIMENSION(:,:) :: XLE_V_CC       ! MEB: latent heat flux from vegetation canopy overstory [J/m2]
  REAL, POINTER, DIMENSION(:,:) :: XLE_G_CC       ! MEB: latent heat flux from understory [J/m2]
  REAL, POINTER, DIMENSION(:,:) :: XLE_C_AC       ! MEB: latent heat flux from canopy air space to the atmosphere [J/m2] 
                                                  !      NOTE total latent heat flux to the atmosphere also possibly 
                                                  !      includes a contribution from snow covering the canopy
  REAL, POINTER, DIMENSION(:,:) :: XLE_N_CC       ! MEB: latent heat flux from the snow on the ground [J/m2]
                                                  !      NOTE total latent heat flux from the snowpack
                                                  !      possibly includes a contribution from snow covering the canopy
  REAL, POINTER, DIMENSION(:,:) :: XSWNET_VC      ! MEB: net vegetation canopy shortwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XSWNET_GC      ! MEB: net ground shortwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XSWNET_NC      ! MEB: net snow shortwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XSWNET_NSC     ! MEB: net snow shortwave radiation for *surface* layer 
                                                  !     (i.e. net snow shortwave radiation less absorbed radiation) [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XLWNET_VC      ! MEB: net vegetation canopy longwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XLWNET_GC      ! MEB: net ground longwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XLWNET_NC      ! MEB: net snow longwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XH_V_CC        ! MEB: sensible heat flux from vegetation canopy overstory [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XH_G_CC        ! MEB: sensible heat flux from understory [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XH_C_AC        ! MEB: sensible heat flux from canopy air space to the atmosphere [W/m2] 
                                                  !      NOTE total sensible heat flux to the atmosphere also possibly 
                                                  !      includes a contribution from snow covering the canopy
  REAL, POINTER, DIMENSION(:,:) :: XH_N_CC        ! MEB: sensible heat flux from the snow on the ground [W/m2]
                                                  !      NOTE total sensible heat flux from the snowpack
                                                  !      possibly includes a contribution from snow covering the canopy
  REAL, POINTER, DIMENSION(:,:) :: XSR_GNC        ! MEB: snow unloading rate from the overstory reservoir [kg/m2/s]
  REAL, POINTER, DIMENSION(:,:) :: XMELTCVC       ! MEB: snow melt rate from the overstory snow reservoir [kg/m2/s]
  REAL, POINTER, DIMENSION(:,:) :: XFRZCVC        ! MEB: snow refreeze rate from the overstory snow reservoir [kg/m2/s]
  REAL, POINTER, DIMENSION(:,:) :: XSWDOWN_GNC    ! MEB: total shortwave radiation transmitted through the canopy
                                                  !      reaching the snowpack/ground understory [W/m2]
  REAL, POINTER, DIMENSION(:,:) :: XLWDOWN_GNC    ! MEB: total shortwave radiation transmitted through and emitted by the canopy
                                                  !      reaching the snowpack/ground understory (explicit part) [W/m2]
!
  REAL, POINTER, DIMENSION(:,:) :: XDWGC         ! liquid soil moisture time tendencies         (kg/m2)
  REAL, POINTER, DIMENSION(:,:) :: XDWGIC        ! solid soil moisture time tendencies          (kg/m2)
  REAL, POINTER, DIMENSION(:,:) :: XDWRC         ! canopy water time tendencies                 (kg/m2)
  REAL, POINTER, DIMENSION(:,:) :: XDSWEC        ! snow water equivalent time tendencies        (kg/m2)
  REAL, POINTER, DIMENSION(:,:) :: XWATBUDC      ! ISBA water budget                            (kg/m2)
!
!* average budget summation variables
!
  REAL, POINTER, DIMENSION(:)   :: XAVG_RNC       ! net radiation at surface                     (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_HC        ! sensible heat flux                           (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LEC       ! total latent heat flux                       (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LEIC      ! sublimation latent heat flux                 (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_GFLUXC    ! net soil-vegetation flux                     (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LEGC      ! latent heat of evaporation over the ground   (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LEGIC     ! surface soil ice sublimation                 (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LEVC      ! latent heat of evaporation over vegetation   (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LESAC      ! latent heat of sublimation over the snow     (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LESLC     ! latent heat of evaporation over the snow     (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LERC      ! evaporation from canopy water interception   (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LETRC     ! evapotranspiration of the vegetation         (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_EVAPC     ! evapotranspiration                           (kg/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_SUBLC     ! sublimation                                  (kg/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_SNDRIFTC  ! blowing snow sublimation (ES or Crocus)      (kg/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_DRAINC    ! soil drainage flux                           (kg/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_QSBC      ! lateral subsurface flux (dif option)         (kg/m2) 
  REAL, POINTER, DIMENSION(:)   :: XAVG_RUNOFFC   ! sub-grid and supersaturation runoff          (kg/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_HORTC     ! sub-grid Horton runoff from the SGH scheme   (kg/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_DRIPC     ! dripping from the vegetation reservoir       (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_MELTC     ! snow melt                                    (kg/m2)   
  REAL, POINTER, DIMENSION(:)   :: XAVG_IFLOODC   ! Floodplains infiltration                     (kg/m2)      
  REAL, POINTER, DIMENSION(:)   :: XAVG_PFLOODC   ! Precipitation intercepted by the floodplains (kg/m2)      
  REAL, POINTER, DIMENSION(:)   :: XAVG_LE_FLOODC ! Floodplains evapotration                     (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LEI_FLOODC! Floodplains evapotration                     (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_ICEFLUXC  ! ice calving flux                             (kg/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_RRVEGC    ! precipitation intercepted by the vegetation  (kg/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_IRRIG_FLUXC! irrigation rate (as soil input)             (kg/m2)
!  
  REAL, POINTER, DIMENSION(:)   :: XAVG_GPPC      ! Gross Primary Production                     (kgCO2/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_RESPC_AUTO! Autotrophic respiration                      (kgCO2/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_RESPC_ECO ! Ecosystem respiration                        (kgCO2/m2)  
!
  REAL, POINTER, DIMENSION(:) :: XAVG_LEVCVC        ! MEB: total evapotranspiration from vegetation canopy overstory [J/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_LESCC         ! MEB: total snow sublimation from vegetation canopy overstory [J/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_LETRGVC       ! MEB: transpiration from understory vegetation [J/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_LETRCVC       ! MEB: transpiration from overstory canopy vegetation [J/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_LERGVC        ! MEB: interception evaporation from understory vegetation [J/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_LERCVC        ! MEB: interception evaporation from overstory canopy vegetation [J/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_LE_V_CC       ! MEB: latent heat flux from vegetation canopy overstory [J/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_LE_G_CC       ! MEB: latent heat flux from understory [J/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_LE_C_AC       ! MEB: latent heat flux from canopy air space to the atmosphere [J/m2] 
                                                  !      NOTE total latent heat flux to the atmosphere also possibly 
                                                  !      includes a contribution from snow covering the canopy
  REAL, POINTER, DIMENSION(:) :: XAVG_LE_N_CC       ! MEB: latent heat flux from the snow on the ground [J/m2]
                                                  !      NOTE total latent heat flux from the snowpack
                                                  !      possibly includes a contribution from snow covering the canopy
  REAL, POINTER, DIMENSION(:) :: XAVG_SWNET_VC      ! MEB: net vegetation canopy shortwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_SWNET_GC      ! MEB: net ground shortwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_SWNET_NC      ! MEB: net snow shortwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_SWNET_NSC     ! MEB: net snow shortwave radiation for *surface* layer 
                                                  !     (i.e. net snow shortwave radiation less absorbed radiation) [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_LWNET_VC      ! MEB: net vegetation canopy longwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_LWNET_GC      ! MEB: net ground longwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_LWNET_NC      ! MEB: net snow longwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_H_V_CC        ! MEB: sensible heat flux from vegetation canopy overstory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_H_G_CC        ! MEB: sensible heat flux from understory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_H_C_AC        ! MEB: sensible heat flux from canopy air space to the atmosphere [W/m2] 
                                                  !      NOTE total sensible heat flux to the atmosphere also possibly 
                                                  !      includes a contribution from snow covering the canopy
  REAL, POINTER, DIMENSION(:) :: XAVG_H_N_CC        ! MEB: sensible heat flux from the snow on the ground [W/m2]
                                                  !      NOTE total sensible heat flux from the snowpack
                                                  !      possibly includes a contribution from snow covering the canopy
  REAL, POINTER, DIMENSION(:) :: XAVG_SR_GNC        ! MEB: snow unloading rate from the overstory reservoir [kg/m2/s]
  REAL, POINTER, DIMENSION(:) :: XAVG_MELTCVC       ! MEB: snow melt rate from the overstory snow reservoir [kg/m2/s]
  REAL, POINTER, DIMENSION(:) :: XAVG_FRZCVC        ! MEB: snow refreeze rate from the overstory snow reservoir [kg/m2/s]
  REAL, POINTER, DIMENSION(:) :: XAVG_SWDOWN_GNC    ! MEB: total shortwave radiation transmitted through the canopy
                                                  !      reaching the snowpack/ground understory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XAVG_LWDOWN_GNC    ! MEB: total shortwave radiation transmitted through and emitted by the canopy
                                                  !      reaching the snowpack/ground understory (explicit part) [W/m2]
!
  REAL, POINTER, DIMENSION(:)   :: XRAINFALLC     ! input rainfall rate for LWATER_BUDGET        (kg/m2)
  REAL, POINTER, DIMENSION(:)   :: XSNOWFALLC     ! input snowfall rate for LWATER_BUDGET        (kg/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_DWGC      ! liquid soil moisture time tendencies         (kg/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_DWGIC     ! solid soil moisture time tendencies          (kg/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_DWRC      ! canopy water time tendencies                 (kg/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_DSWEC     ! snow water equivalent time tendencies        (kg/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_WATBUDC   ! ISBA water budget                            (kg/m2)
! 
!------------------------------------------------------------------------------
!

END TYPE DIAG_EVAP_ISBA_t

TYPE(DIAG_EVAP_ISBA_t), ALLOCATABLE, TARGET, SAVE :: DIAG_EVAP_ISBA_MODEL(:)

TYPE(DIAG_EVAP_ISBA_t), POINTER :: DIAG_EVAP_ISBA => NULL()
!$OMP THREADPRIVATE(DIAG_EVAP_ISBA)

CONTAINS
!
SUBROUTINE DIAG_EVAP_ISBA_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODD_DIAG_EVAP_ISBA_N:DIAG_EVAP_ISBA_GOTO_MODEL',0,ZHOOK_HANDLE)

DIAG_EVAP_ISBA => DIAG_EVAP_ISBA_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_DIAG_EVAP_ISBA_N:DIAG_EVAP_ISBA_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE DIAG_EVAP_ISBA_GOTO_MODEL

SUBROUTINE DIAG_EVAP_ISBA_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_EVAP_ISBA_N:DIAG_EVAP_ISBA_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(DIAG_EVAP_ISBA_MODEL(KMODEL))
DIAG_EVAP_ISBA => DIAG_EVAP_ISBA_MODEL(KMODEL)
DO J=1,KMODEL
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLEG)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLEGI)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLEV)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLES)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLESL)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLER)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLETR)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XEVAP)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XSUBL)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XSNDRIFT)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XDRAIN)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XQSB)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XRUNOFF)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XHORT)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XRRVEG)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XMELT)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XIFLOOD)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XPFLOOD)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLE_FLOOD)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLEI_FLOOD)
!
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLEVCV)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLESC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLETRGV)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLETRCV)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLERGV)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLERCV)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLE_V_C)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLE_G_C)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLE_C_A)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLE_N_C)
!
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XSWNET_V)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XSWNET_G)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XSWNET_N)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XSWNET_NS)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLWNET_V)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLWNET_G)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLWNET_N)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XSWDOWN_GN)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLWDOWN_GN)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XH_V_C)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XH_G_C)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XH_C_A)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XH_N_C)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XSR_GN)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XMELTCV)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XFRZCV)
!
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XDRIP)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XIRRIG_FLUX)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XGPP)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XRESP_AUTO)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XRESP_ECO)  
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XDWG)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XDWGI)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XDWR)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XDSWE)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XWATBUD)  
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LEG)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LEGI)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LEV)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LES)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LESL)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LER)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LETR)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_EVAP)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_SUBL)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_SNDRIFT)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_DRAIN)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_QSB)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_RUNOFF)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_HORT)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_DRIP)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_MELT)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_IFLOOD)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_PFLOOD)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LE_FLOOD)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LEI_FLOOD)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_RRVEG)
!
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LEVCV)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LESC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LETRGV)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LETRCV)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LERGV)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LERCV)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LE_V_C)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LE_G_C)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LE_C_A)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LE_N_C)
!
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_SWNET_V)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_SWNET_G)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_SWNET_N)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_SWNET_NS)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LWNET_V)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LWNET_G)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LWNET_N)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_SWDOWN_GN)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LWDOWN_GN)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_H_V_C)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_H_G_C)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_H_C_A)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_H_N_C)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_SR_GN)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_MELTCV)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_FRZCV)
!
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_IRRIG_FLUX)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_GPP)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_RESP_AUTO)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_RESP_ECO)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XRAINFALL)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XSNOWFALL)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_DWG)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_DWGI)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_DWR)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_DSWE)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_WATBUD)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XRNC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XHC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLEC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLEIC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XGFLUXC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLEGC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLEGIC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLEVC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLESAC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLESLC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLERC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLETRC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XEVAPC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XSUBLC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XSNDRIFTC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XDRAINC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XQSBC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XRUNOFFC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XHORTC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XDRIPC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XMELTC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XIFLOODC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XPFLOODC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLE_FLOODC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLEI_FLOODC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XICEFLUXC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XRRVEGC)
!
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLEVCVC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLESCC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLETRGVC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLETRCVC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLERGVC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLERCVC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLE_V_CC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLE_G_CC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLE_C_AC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLE_N_CC)
!
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XSWNET_VC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XSWNET_GC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XSWNET_NC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XSWNET_NSC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLWNET_VC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLWNET_GC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLWNET_NC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XSWDOWN_GNC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XLWDOWN_GNC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XH_V_CC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XH_G_CC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XH_C_AC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XH_N_CC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XSR_GNC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XMELTCVC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XFRZCVC)
!
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XIRRIG_FLUXC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XGPPC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XRESPC_AUTO)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XRESPC_ECO) 
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XDWGC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XDWGIC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XDWRC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XDSWEC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XWATBUDC) 
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_RNC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_HC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LEC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LEIC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_GFLUXC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LEGC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LEGIC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LEVC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LESAC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LESLC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LERC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LETRC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_EVAPC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_SUBLC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_SNDRIFTC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_DRAINC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_QSBC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_RUNOFFC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_HORTC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_DRIPC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_MELTC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_IFLOODC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_PFLOODC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LE_FLOODC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LEI_FLOODC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_ICEFLUXC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_RRVEGC)
!
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LEVCVC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LESCC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LETRGVC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LETRCVC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LERGVC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LERCVC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LE_V_CC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LE_G_CC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LE_C_AC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LE_N_CC)
!
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_SWNET_VC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_SWNET_GC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_SWNET_NC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_SWNET_NSC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LWNET_VC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LWNET_GC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LWNET_NC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_SWDOWN_GNC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_LWDOWN_GNC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_H_V_CC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_H_G_CC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_H_C_AC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_H_N_CC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_SR_GNC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_MELTCVC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_FRZCVC)
!
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_IRRIG_FLUXC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_GPPC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_RESPC_AUTO)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_RESPC_ECO)  
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XRAINFALLC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XSNOWFALLC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_DWGC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_DWGIC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_DWRC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_DSWEC)
  NULLIFY(DIAG_EVAP_ISBA_MODEL(J)%XAVG_WATBUDC)  
ENDDO
DIAG_EVAP_ISBA_MODEL(:)%LSURF_EVAP_BUDGET=.FALSE.
DIAG_EVAP_ISBA_MODEL(:)%LSURF_BUDGETC=.FALSE.
DIAG_EVAP_ISBA_MODEL(:)%LRESET_BUDGETC=.FALSE.
DIAG_EVAP_ISBA_MODEL(:)%LWATER_BUDGET=.FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_EVAP_ISBA_N:DIAG_EVAP_ISBA_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_EVAP_ISBA_ALLOC

SUBROUTINE DIAG_EVAP_ISBA_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_EVAP_ISBA_N:DIAG_EVAP_ISBA_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(DIAG_EVAP_ISBA_MODEL)) DEALLOCATE(DIAG_EVAP_ISBA_MODEL)
IF (ASSOCIATED(DIAG_EVAP_ISBA)) NULLIFY(DIAG_EVAP_ISBA)
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_EVAP_ISBA_N:DIAG_EVAP_ISBA_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_EVAP_ISBA_DEALLO

END MODULE MODD_DIAG_EVAP_ISBA_n
