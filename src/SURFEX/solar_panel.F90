!     #########
    SUBROUTINE SOLAR_PANEL(TPN, DGMT, PTSTEP, PTSUN, PRESIDENTIAL, PEMIT_LW_ROOF, &
                           PEMIT_LWDN_PANEL, PLW_RAD, PTA, PN_FLOOR, PPROD_BLD )
!   ##########################################################################
!
!!****  *SOLAR_PANEL*  
!!
!!    PURPOSE
!!    -------
!
! Computes the energy budget of the solar panels.
!         
!     
!!**  METHOD
!     ------
!
!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!      V. Masson           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    08/2013
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_TEB_PANEL_n, ONLY : TEB_PANEL_t
USE MODD_DIAG_MISC_TEB_n, ONLY : DIAG_MISC_TEB_1P_t
!
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_CSTS,     ONLY : XSTEFAN
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(TEB_PANEL_t), INTENT(INOUT) :: TPN
TYPE(DIAG_MISC_TEB_1P_t), INTENT(INOUT) :: DGMT
!
REAL,               INTENT(IN)  :: PTSTEP          ! time step  (s)
REAL, DIMENSION(:), INTENT(IN)  :: PTSUN           ! solar time (s since solar midnight)
REAL, DIMENSION(:), INTENT(IN)  :: PRESIDENTIAL    ! Buildings Residential use fraction        (-)
REAL, DIMENSION(:), INTENT(IN)  :: PEMIT_LW_ROOF   ! Upwards   LW flux from roof               (W/m2)
REAL, DIMENSION(:), INTENT(IN)  :: PEMIT_LWDN_PANEL! Downwards LW flux from panel              (W/m2)
REAL, DIMENSION(:), INTENT(IN)  :: PLW_RAD         ! Incoming Longwave radiation               (W/m2)
REAL, DIMENSION(:), INTENT(IN)  :: PTA             ! Air temperature                           (K)
REAL, DIMENSION(:), INTENT(IN)  :: PN_FLOOR        ! number of floors                          (-)
!
REAL, DIMENSION(:), INTENT(OUT)  :: PPROD_BLD
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*      0.2    declarations of parameters
!
! coefficient to take into account irradiance to estimate surface temperature of panel
REAL, PARAMETER :: XKT    = 0.05   ! (Km2/W)
!
! coefficient to take into account optimum orientation of the solar panel
REAL, PARAMETER :: XFT    = 1.10   ! (-)
!
! optimum panel temperature for photovoltaic production
REAL, PARAMETER :: XT_OPT = 298.15 ! (K)
!
! coefficient to take into account decrease of production for warm temperatures
REAL, PARAMETER :: XT_LOSS= 0.005  ! (K-1)
!
! difference of temperature between cold and hot water
REAL, PARAMETER :: XWATER_DT = 45. ! (K)
!
! density of m2 of panel per m2 of floor necessary for hot water production
REAL, PARAMETER :: XTHER_FLOOR= 1./30. ! (-)
!
! Annual thermal production for an increase of 1K of the water temperature
REAL, PARAMETER :: XTHER_RATE= 0.5 * 1.16*32. ! (kWh/year/m2 panel)
!
!! Maximum irradiance above which there is no more heating necessary (target temperature reached quickly)
!REAL, PARAMETER :: XTHER_IRR = 500.           ! (W/m2)
! Efficiency for thermal panels
REAL, PARAMETER :: XTHER_EFF = 0.60           ! (-)
!
!
!
!
!*      0.3    declarations of local variables
!
REAL                       :: ZTHER_DAILY_TARGET ! daily target for thermal production         (J/m2)
REAL, DIMENSION(SIZE(PTA)) :: ZIRRADIANCE ! incoming solar radiation normal to the solar panel (W/m2)
REAL, DIMENSION(SIZE(PTA)) :: ZTS_PANEL   ! Surface temperature of the upwards face of the solar panel (K)
REAL, DIMENSION(SIZE(PTA)) :: ZLWU_PANEL  ! Upwards longwave radiation from the solar panel (W/m2)
REAL, DIMENSION(SIZE(PTA)) :: ZTHER_FRAC  ! Fraction of thermal      panels per m2 of roof  (m2 panel/m2 bld)
REAL, DIMENSION(SIZE(PTA)) :: ZPHOT_FRAC  ! Fraction of photovoltaic panels per m2 of roof  (m2 panel/m2 bld)
REAL, DIMENSION(SIZE(PTA)) :: ZTHER_PRODC_DAY ! guess of daily production by thermal panels (J/m2)
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SOLAR_PANEL',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!*      0.4    Default values for output variables
!              -----------------------------------
!
DGMT%XABS_LW_PANEL   = XUNDEF ! Absorbed LW enerby by solar panel         (W/m2)
DGMT%XH_PANEL        = XUNDEF ! Sensible heat released by the solar panel (W/m2)
DGMT%XRN_PANEL       = XUNDEF ! Net radiation     of the solar panel      (W/m2)
DGMT%XTHER_PROD_PANEL= XUNDEF ! Thermal      Energy production of the solar panel      (W/m2)
DGMT%XPHOT_PROD_PANEL= XUNDEF ! Photovoltaic Energy production of the solar panel      (W/m2)
DGMT%XPROD_PANEL     = XUNDEF ! Averaged     Energy production of the solar panel      (W/m2)
DGMT%XTHER_PROD_BLD  = 0.     ! Thermal      Energy production of the solar panel      (W/m2)
DGMT%XPHOT_PROD_BLD  = 0.     ! Photovoltaic Energy production of the solar panel      (W/m2)
PPROD_BLD       = 0.     ! Averaged     Energy production of the solar panel      (W/m2)
!
!-------------------------------------------------------------------------------
!
!*      1.     Initializations
!              ---------------
!
!*      1.1    Checks hour for daily thermal production reset
!              ----------------------------------------------
!
!* energy reset between midnight and 1AM.
WHERE (PTSUN(:)>=0. .AND. PTSUN(:)<=3600.)  TPN%XTHER_PRODC_DAY(:) = 0.
!
!
!*      1.2    Daily target production for thermal panels
!              ------------------------------------------
!
ZTHER_DAILY_TARGET = (2.*XTHER_RATE) * XWATER_DT * (1000. / 365. * 3600. ) ! (J/m2)
! the factor 2 is to remove the assumption of cloudy days in the annual mean production
!
!-------------------------------------------------------------------------------
!* Note that computations are done only where solar panels are present
WHERE (TPN%XFRAC_PANEL(:)>0.)
!-------------------------------------------------------------------------------
!
!*      2.     Irradiance on panel
!              -------------------
!
  ZIRRADIANCE(:) = XFT * DGMT%XABS_SW_PANEL(:) / (1.-TPN%XALB_PANEL(:) ) 
!
!-------------------------------------------------------------------------------
!
!*      3.     Solar panel temperature
!              -----------------------
!
  ZTS_PANEL  (:) = PTA(:) + XKT * ZIRRADIANCE(:)
!
!-------------------------------------------------------------------------------
!
!*      4.     Upwards solar panel LW radiation
!              --------------------------------
!
  ZLWU_PANEL  (:) =       TPN%XEMIS_PANEL(:)  * XSTEFAN * ZTS_PANEL(:)**4 &
                    + (1.-TPN%XEMIS_PANEL(:)) * PLW_RAD(:)
!
!-------------------------------------------------------------------------------
!
!*      5.     Solar panel LW budget
!              ---------------------
!
  DGMT%XABS_LW_PANEL(:)= PLW_RAD(:) + PEMIT_LW_ROOF(:) - PEMIT_LWDN_PANEL(:) - ZLWU_PANEL(:)
!
!-------------------------------------------------------------------------------
!
!*      6.     Solar panel Net radiation
!              -------------------------
!
  DGMT%XRN_PANEL(:)    = DGMT%XABS_SW_PANEL(:) + DGMT%XABS_LW_PANEL(:)
!
!-------------------------------------------------------------------------------
!
!*      7.     Fraction of panel surface types
!              -------------------------------
!
!*      7.1    Panel dedicated to thermal production of hot water 
!              --------------------------------------------------
!
  ZTHER_FRAC(:) = MIN( XTHER_FLOOR * PN_FLOOR(:) * PRESIDENTIAL(:), TPN%XFRAC_PANEL(:) )   ! (m2 thermal       panel / m2 roof)
!
!*      7.2    Photovoltaic panel
!              ------------------
!
  ZPHOT_FRAC(:) =  TPN%XFRAC_PANEL(:) - ZTHER_FRAC(:)                      ! (m2 photovoltaic panel / m2 roof)
!
!-------------------------------------------------------------------------------
!
!*      8.     Thermal Production of hot water (W/m2 thermal panel)
!              -------------------------------
!
!*      8.1    Instantaneous production
!              ------------------------
!
  DGMT%XTHER_PROD_PANEL(:)= XTHER_EFF * ZIRRADIANCE(:)                       ! (W/m2)
!  DGMT%XTHER_PROD_PANEL(:) = XTHER_RATE * XWATER_DT * (1000. / 24. / 365.)
!
!*      8.2    Integrated daily production
!              ---------------------------
!
  ZTHER_PRODC_DAY(:) = TPN%XTHER_PRODC_DAY(:) + PTSTEP * DGMT%XTHER_PROD_PANEL(:) ! (J/m2)
!
!*      8.3    Daily production limited by daily target
!              ----------------------------------------
!
  ZTHER_PRODC_DAY(:) = MIN( ZTHER_PRODC_DAY(:) , ZTHER_DAILY_TARGET )
!
!*      8.4     Instantaneous production taking into account target limit if reached
!               --------------------------------------------------------------------
!
  DGMT%XTHER_PROD_PANEL(:)= ( ZTHER_PRODC_DAY(:) - TPN%XTHER_PRODC_DAY(:) ) / PTSTEP
!
!*      8.5    Updates daily production
!              ------------------------
!
  TPN%XTHER_PRODC_DAY(:) = ZTHER_PRODC_DAY(:)

!-------------------------------------------------------------------------------
!
!*      9.     Photovoltaic Production  (W/m2 photovoltaic panel)
!              -----------------------
!
  DGMT%XPHOT_PROD_PANEL(:) = TPN%XEFF_PANEL(:) * ZIRRADIANCE(:) * MIN(1.,1.-XT_LOSS*(ZTS_PANEL(:)-XT_OPT))
!
!-------------------------------------------------------------------------------
!
!*     10.     Averaged Production  (W/m2 panel)
!              -------------------
!
  DGMT%XPROD_PANEL(:) =  DGMT%XTHER_PROD_PANEL(:) * (ZTHER_FRAC(:) / TPN%XFRAC_PANEL(:)) &
                  + DGMT%XPHOT_PROD_PANEL(:) * (ZPHOT_FRAC(:) / TPN%XFRAC_PANEL(:))
!
!-------------------------------------------------------------------------------
!
!*     11.     Sensible heat flux  (W/m2 panel)
!              ------------------
!
  DGMT%XH_PANEL(:) = DGMT%XRN_PANEL(:) - DGMT%XPROD_PANEL(:)
!
!-------------------------------------------------------------------------------
!
!*     12.     Productions per building (W/m2 bld)
!              ------------------------
!
  DGMT%XTHER_PROD_BLD(:) = DGMT%XTHER_PROD_PANEL(:) * ZTHER_FRAC(:)
  DGMT%XPHOT_PROD_BLD(:) = DGMT%XPHOT_PROD_PANEL(:) * ZPHOT_FRAC(:)
  PPROD_BLD     (:) = DGMT%XTHER_PROD_BLD  (:) + DGMT%XPHOT_PROD_PANEL(:)
!
!-------------------------------------------------------------------------------
END WHERE
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SOLAR_PANEL',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE SOLAR_PANEL
