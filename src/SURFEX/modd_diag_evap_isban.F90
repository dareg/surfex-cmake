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
  REAL, POINTER, DIMENSION(:) :: XLEG          ! latent heat of evaporation over the ground   (W/m2)
  REAL, POINTER, DIMENSION(:) :: XLEGI         ! surface soil ice sublimation                 (W/m2)
  REAL, POINTER, DIMENSION(:) :: XLEV          ! latent heat of evaporation over vegetation   (W/m2)
  REAL, POINTER, DIMENSION(:) :: XLES          ! latent heat of sublimation over the snow     (W/m2)
  REAL, POINTER, DIMENSION(:) :: XLESL         ! latent heat of evaporation over the snow     (W/m2)
  REAL, POINTER, DIMENSION(:) :: XLER          ! evaporation from canopy water interception   (W/m2)
  REAL, POINTER, DIMENSION(:) :: XLETR         ! evapotranspiration of the vegetation         (W/m2)
  REAL, POINTER, DIMENSION(:) :: XUSTAR       ! friction velocity from snow-free 
!                                               ! surface (ISBA-ES:3-L)                         (m/s  
  REAL, POINTER, DIMENSION(:) :: XSNDRIFT      ! blowing snow sublimation (ES or Crocus)      (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XRESTORE   ! surface energy budget restore term            (W/m2)  
  REAL, POINTER, DIMENSION(:) :: XDRAIN        ! soil drainage flux                           (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XQSB          ! lateral subsurface flux (dif option)         (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XRUNOFF       ! sub-grid and supersaturation runoff          (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XHORT         ! sub-grid Horton runoff from the SGH scheme   (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XRRVEG        !  precipitation intercepted by the vegetation (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XMELT         ! snow melt                                    (kg/m2/s)  
  REAL, POINTER, DIMENSION(:) :: XMELTADV   ! advective energy from snow melt water 
!                                               ! (ISBA-ES:3-L)                                 (W/m2)
  REAL, POINTER, DIMENSION(:) :: XIFLOOD       ! Floodplains infiltration                     (kg/m2/s)      
  REAL, POINTER, DIMENSION(:) :: XPFLOOD       ! Precipitation intercepted by the floodplains (kg/m2/s)      
  REAL, POINTER, DIMENSION(:) :: XLE_FLOOD     ! Floodplains evapotration                     (W/m2)      
  REAL, POINTER, DIMENSION(:) :: XLEI_FLOOD    ! Floodplains evapotration                     (W/m2) 
  REAL, POINTER, DIMENSION(:) :: XICEFLUX    ! ice calving flux                             (kg/m2)  
  REAL, POINTER, DIMENSION(:) :: XDRIP         ! dripping from the vegetation reservoir       (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XIRRIG_FLUX   ! irrigation rate (as soil input)              (kg/m2/s)
!  
  REAL, POINTER, DIMENSION(:) :: XGPP          ! Gross Primary Production                     (kgCO2/m2/s)
  REAL, POINTER, DIMENSION(:) :: XRESP_AUTO    ! Autotrophic respiration                      (kgCO2/m2/s)
  REAL, POINTER, DIMENSION(:) :: XRESP_ECO     ! Ecosystem respiration                        (kgCO2/m2/s)

  REAL, POINTER, DIMENSION(:) :: XLEVCV        ! MEB: total evapotranspiration from vegetation canopy overstory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XLESC         ! MEB: total snow sublimation from vegetation canopy overstory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XLETRGV       ! MEB: transpiration from understory vegetation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XLETRCV       ! MEB: transpiration from overstory canopy vegetation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XLERGV        ! MEB: interception evaporation from understory vegetation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XLELITTER        ! MEB: interception evaporation from understory vegetation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XLELITTERI        ! MEB: interception evaporation from understory vegetation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XDRIPLIT       ! 
  REAL, POINTER, DIMENSION(:) :: XRRLIT         ! 
  REAL, POINTER, DIMENSION(:) :: XLERCV        ! MEB: interception evaporation from overstory canopy vegetation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XLE_V_C       ! MEB: latent heat flux from vegetation canopy overstory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XLE_G_C       ! MEB: latent heat flux from understory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XLE_C_A       ! MEB: latent heat flux from canopy air space to the atmosphere [W/m2] 
                                                 !      NOTE total latent heat flux to the atmosphere also possibly 
                                                 !      includes a contribution from snow covering the canopy
  REAL, POINTER, DIMENSION(:) :: XLE_N_C       ! MEB: latent heat flux from the snow on the ground [W/m2]
                                                 !      NOTE total latent heat flux from the snowpack
                                                 !      possibly includes a contribution from snow covering the canopy
  REAL, POINTER, DIMENSION(:) :: XEVAP_N_C      ! MEB: Total evap from snow on the ground to canopy air space  [kg/m2/s]
  REAL, POINTER, DIMENSION(:) :: XEVAP_G_C      ! MEB: Total evap from ground to canopy air space [kg/m2/s]                                                 
  REAL, POINTER, DIMENSION(:) :: XSWUP       ! MEB: net *total* (surface) upwelling shortwave radiation to atmosphere [W/m2]

  REAL, POINTER, DIMENSION(:) :: XSWNET_V      ! MEB: net vegetation canopy shortwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XSWNET_G      ! MEB: net ground shortwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XSWNET_N      ! MEB: net snow shortwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XSWNET_NS     ! MEB: net snow shortwave radiation for *surface* layer 
                                                 !     (i.e. net snow shortwave radiation less absorbed radiation) [W/m2]
  REAL, POINTER, DIMENSION(:) :: XLWUP          ! MEB: net *total* (surface) upwelling longwave radiation to atmosphere [W/m2]    
  REAL, POINTER, DIMENSION(:) :: XLWNET_V      ! MEB: net vegetation canopy longwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XLWNET_G      ! MEB: net ground longwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XLWNET_N      ! MEB: net snow longwave radiation [W/m2]
  REAL, POINTER, DIMENSION(:) :: XH_V_C        ! MEB: sensible heat flux from vegetation canopy overstory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XH_G_C        ! MEB: sensible heat flux from understory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XH_C_A        ! MEB: sensible heat flux from canopy air space to the atmosphere [W/m2] 
                                                 !      NOTE total sensible heat flux to the atmosphere also possibly 
                                                 !      includes a contribution from snow covering the canopy
  REAL, POINTER, DIMENSION(:) :: XH_N_C        ! MEB: sensible heat flux from the snow on the ground [W/m2]
                                                 !      NOTE total sensible heat flux from the snowpack
                                                 !      possibly includes a contribution from snow covering the canopy
  REAL, POINTER, DIMENSION(:) :: XSR_GN        ! MEB: snow unloading rate from the overstory reservoir [kg/m2/s]
  REAL, POINTER, DIMENSION(:) :: XMELTCV       ! MEB: snow melt rate from the overstory snow reservoir [kg/m2/s]
  REAL, POINTER, DIMENSION(:) :: XFRZCV        ! MEB: snow refreeze rate from the overstory snow reservoir [kg/m2/s]
  REAL, POINTER, DIMENSION(:) :: XSWDOWN_GN    ! MEB: total shortwave radiation transmitted through the canopy
                                                 !      reaching the snowpack/ground understory [W/m2]
  REAL, POINTER, DIMENSION(:) :: XLWDOWN_GN    ! MEB: total shortwave radiation transmitted through and emitted by the canopy
                                                 !      reaching the snowpack/ground understory (explicit part) [W/m2]
!
  REAL, POINTER, DIMENSION(:) :: XDWG          ! liquid soil moisture time tendencies         (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XDWGI         ! solid soil moisture time tendencies          (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XDWR          ! canopy water time tendencies                 (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XDSWE         ! snow water equivalent time tendencies        (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XWATBUD       ! ISBA water budget                            (kg/m2/s)
!
  REAL, POINTER, DIMENSION(:)   :: XRAINFALL     ! input rainfall rate for LWATER_BUDGET        (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XSNOWFALL     ! input snowfall rate for LWATER_BUDGET        (kg/m2/s)
!  
!------------------------------------------------------------------------------
!

END TYPE DIAG_EVAP_ISBA_t
!
TYPE DIAG_EVAP_ISBA_PATCH_t
!
TYPE(DIAG_EVAP_ISBA_t), ALLOCATABLE :: AL(:) 
!
END TYPE DIAG_EVAP_ISBA_PATCH_t
!
CONTAINS
!
SUBROUTINE DIAG_EVAP_ISBA_PATCH_INIT(YDIAG_EVAP_ISBA_PATCH,KPATCH)
TYPE(DIAG_EVAP_ISBA_PATCH_t), INTENT(INOUT) :: YDIAG_EVAP_ISBA_PATCH 
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_EVAP_ISBA_N:DIAG_EVAP_ISBA_PATCH_INIT",0,ZHOOK_HANDLE)
 ALLOCATE(YDIAG_EVAP_ISBA_PATCH%AL(KPATCH))
DO JP=1,KPATCH
  CALL DIAG_EVAP_ISBA_INIT(YDIAG_EVAP_ISBA_PATCH%AL(JP))
ENDDO         
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_EVAP_ISBA_N:DIAG_EVAP_ISBA_PATCH_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_EVAP_ISBA_PATCH_INIT
!
SUBROUTINE DIAG_EVAP_ISBA_INIT(YDIAG_EVAP_ISBA)
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: YDIAG_EVAP_ISBA
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_EVAP_ISBA_N:DIAG_EVAP_ISBA_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YDIAG_EVAP_ISBA%XLEG)
  NULLIFY(YDIAG_EVAP_ISBA%XLEGI)
  NULLIFY(YDIAG_EVAP_ISBA%XLEV)
  NULLIFY(YDIAG_EVAP_ISBA%XLES)
  NULLIFY(YDIAG_EVAP_ISBA%XLESL)
  NULLIFY(YDIAG_EVAP_ISBA%XLER)
  NULLIFY(YDIAG_EVAP_ISBA%XLETR)
  NULLIFY(YDIAG_EVAP_ISBA%XUSTAR)
  NULLIFY(YDIAG_EVAP_ISBA%XSNDRIFT)
  NULLIFY(YDIAG_EVAP_ISBA%XRESTORE)
  NULLIFY(YDIAG_EVAP_ISBA%XDRAIN)
  NULLIFY(YDIAG_EVAP_ISBA%XQSB)
  NULLIFY(YDIAG_EVAP_ISBA%XRUNOFF)
  NULLIFY(YDIAG_EVAP_ISBA%XHORT)
  NULLIFY(YDIAG_EVAP_ISBA%XRRVEG)
  NULLIFY(YDIAG_EVAP_ISBA%XMELT)
  NULLIFY(YDIAG_EVAP_ISBA%XMELTADV)
  NULLIFY(YDIAG_EVAP_ISBA%XIFLOOD)
  NULLIFY(YDIAG_EVAP_ISBA%XPFLOOD)
  NULLIFY(YDIAG_EVAP_ISBA%XLE_FLOOD)
  NULLIFY(YDIAG_EVAP_ISBA%XLEI_FLOOD)
!
  NULLIFY(YDIAG_EVAP_ISBA%XICEFLUX)  
!
  NULLIFY(YDIAG_EVAP_ISBA%XLEVCV)
  NULLIFY(YDIAG_EVAP_ISBA%XLESC)
  NULLIFY(YDIAG_EVAP_ISBA%XLETRGV)
  NULLIFY(YDIAG_EVAP_ISBA%XLETRCV)
  NULLIFY(YDIAG_EVAP_ISBA%XLERGV)

  NULLIFY(YDIAG_EVAP_ISBA%XLELITTER)
  NULLIFY(YDIAG_EVAP_ISBA%XLELITTERI)
  NULLIFY(YDIAG_EVAP_ISBA%XDRIPLIT)
  NULLIFY(YDIAG_EVAP_ISBA%XRRLIT)

  NULLIFY(YDIAG_EVAP_ISBA%XLERCV)
  NULLIFY(YDIAG_EVAP_ISBA%XLE_V_C)
  NULLIFY(YDIAG_EVAP_ISBA%XLE_G_C)
  NULLIFY(YDIAG_EVAP_ISBA%XLE_C_A)
  NULLIFY(YDIAG_EVAP_ISBA%XLE_N_C)
!
  NULLIFY(YDIAG_EVAP_ISBA%XEVAP_N_C)
  NULLIFY(YDIAG_EVAP_ISBA%XEVAP_G_C)

  NULLIFY(YDIAG_EVAP_ISBA%XSWUP)
  NULLIFY(YDIAG_EVAP_ISBA%XSWNET_V)
  NULLIFY(YDIAG_EVAP_ISBA%XSWNET_G)
  NULLIFY(YDIAG_EVAP_ISBA%XSWNET_N)
  NULLIFY(YDIAG_EVAP_ISBA%XSWNET_NS)
  NULLIFY(YDIAG_EVAP_ISBA%XLWUP)
  NULLIFY(YDIAG_EVAP_ISBA%XLWNET_V)
  NULLIFY(YDIAG_EVAP_ISBA%XLWNET_G)
  NULLIFY(YDIAG_EVAP_ISBA%XLWNET_N)
  NULLIFY(YDIAG_EVAP_ISBA%XSWDOWN_GN)
  NULLIFY(YDIAG_EVAP_ISBA%XLWDOWN_GN)
  NULLIFY(YDIAG_EVAP_ISBA%XH_V_C)
  NULLIFY(YDIAG_EVAP_ISBA%XH_G_C)
  NULLIFY(YDIAG_EVAP_ISBA%XH_C_A)
  NULLIFY(YDIAG_EVAP_ISBA%XH_N_C)
  NULLIFY(YDIAG_EVAP_ISBA%XSR_GN)
  NULLIFY(YDIAG_EVAP_ISBA%XMELTCV)
  NULLIFY(YDIAG_EVAP_ISBA%XFRZCV)
!
  NULLIFY(YDIAG_EVAP_ISBA%XDRIP)
  NULLIFY(YDIAG_EVAP_ISBA%XIRRIG_FLUX)
  NULLIFY(YDIAG_EVAP_ISBA%XGPP)
  NULLIFY(YDIAG_EVAP_ISBA%XRESP_AUTO)
  NULLIFY(YDIAG_EVAP_ISBA%XRESP_ECO)  
  NULLIFY(YDIAG_EVAP_ISBA%XDWG)
  NULLIFY(YDIAG_EVAP_ISBA%XDWGI)
  NULLIFY(YDIAG_EVAP_ISBA%XDWR)
  NULLIFY(YDIAG_EVAP_ISBA%XDSWE)
  NULLIFY(YDIAG_EVAP_ISBA%XWATBUD)  
!
  NULLIFY(YDIAG_EVAP_ISBA%XRAINFALL)
  NULLIFY(YDIAG_EVAP_ISBA%XSNOWFALL)
!
YDIAG_EVAP_ISBA%LSURF_EVAP_BUDGET=.FALSE.
YDIAG_EVAP_ISBA%LSURF_BUDGETC=.FALSE.
YDIAG_EVAP_ISBA%LRESET_BUDGETC=.FALSE.
YDIAG_EVAP_ISBA%LWATER_BUDGET=.FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_EVAP_ISBA_N:DIAG_EVAP_ISBA_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_EVAP_ISBA_INIT
!
END MODULE MODD_DIAG_EVAP_ISBA_n
