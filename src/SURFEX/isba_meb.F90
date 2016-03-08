!     #########
      SUBROUTINE ISBA_MEB(IO, IMX, IMT, IMM, IMI, IP, INI, IR, DGIP, DGEIP, DGMI, G, AG, &
                          TPTIME, OMEB, OSHADE, HSNOW_ISBA, HIMPLICIT_WIND, PTSTEP, &
                          PSOILHCAPZ, PSOILCONDZ, PFROZEN1, PPS, PZENITH,     &
                          PSCA_SW, PSW_RAD, PVMOD, PRR, PSR, PRHOA, PTA, PQA, &
                          PDIRCOSZW, PEXNS, PEXNA, PPET_A_COEF, PPET_B_COEF,  &
                          PPEQ_A_COEF, PPEQ_B_COEF, PPEW_A_COEF, PPEW_B_COEF, &
                          PZREF, PUREF, PZ0G_WITHOUT_SNOW, PZ0_MEBV, PZ0H_MEBV,&
                          PZ0EFF_MEBV, PZ0_MEBN, PZ0H_MEBN, PZ0EFF_MEBN,      &
                          PALBNIR_TVEG, PALBVIS_TVEG, PALBNIR_TSOIL, PALBVIS_TSOIL, &
                          PABC, PIACAN, PPOI, PCSP, PRESP_BIOMASS_INST, PPALPHAN,   &
                          PF2, PLW_RAD, PGRNDFLUX, PFLSN_COR, PUSTAR, PEMIST,       &
                          PHU_AGG, PAC_AGG, PDELHEATV_SFC, PDELHEATG_SFC, PDELHEATG,&
                          PDELHEATN, PDELHEATN_SFC, PRESTOREN, PTDEEP_A, PDEEP_FLUX, &
                          PRISNOW, PSNOW_THRUFAL, PEVAPCOR, PSUBVCOR, PSNOWSFCH, PQSNOW)
!     ##########################################################################
!
!                             
!!****  *isba_meb*  
!!
!!    PURPOSE
!!    -------
!       Monitor for the calculation of the surface fluxes and of the
!     prognostic variables of the surface over natural areas 
!     with an explicit vegetation layer
!
!     NOTE...currently MEB can be coupled with 
!     IO%CISBA='DIF' or '3-L' soil options
!     HSNOW='3-L' snow scheme
!     Soon, HSNOW=CRO and IO%CPHOTO/=NON (i.e. Ags will be added)
!     
!!**  METHOD
!!    ------
!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!    Noilhan and Planton (1989)
!!      
!!    AUTHOR
!!    ------
!!	A. Boone           * Meteo-France *
!!      P. Samuelsson      * SMHI *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       10/2014
!!      (A. Napoly)    09/2015  Add Litter layer option code
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_FIX_t, ISBA_PARAM_TIME_t, ISBA_PARAM_MEB_t, &
                              ISBA_PARAM_IRRIG_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t, ISBA_INIT_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
USE MODD_GRID_n, ONLY : GRID_t
USE MODD_AGRI_n, ONLY : AGRI_t
USE MODD_DIAG_n, ONLY : DIAG_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_CSTS,           ONLY : XCPD, XDAY, XRHOLW 
USE MODD_MEB_PAR,        ONLY : XSW_WGHT_VIS, XSW_WGHT_NIR
USE MODD_ISBA_PAR,       ONLY : XRS_MAX 
USE MODD_DATA_COVER_PAR, ONLY : NVT_SNOW
!
USE MODD_TYPE_DATE_SURF, ONLY : DATE_TIME
!
USE MODE_THERMOS
USE MODE_MEB,            ONLY : SNOW_INTERCEPT_EFF
!
USE MODI_WET_LEAVES_FRAC
USE MODI_VEG
USE MODI_SNOW_LEAVES_FRAC_MEB
USE MODI_PREPS_FOR_MEB_EBUD_RAD
USE MODI_ISBA_LWNET_MEB
USE MODI_DRAG_MEB
USE MODI_E_BUDGET_MEB
USE MODI_ISBA_FLUXES_MEB
USE MODI_SNOW_LOAD_MEB
USE MODI_HYDRO_VEG
USE MODI_SNOW3L_ISBA
USE MODI_RADIATIVE_TRANSFERT
USE MODI_COTWORES
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!              -------------------------
!
!
!* general variables
!  -----------------
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_PARAM_FIX_t), INTENT(INOUT) :: IMX
TYPE(ISBA_PARAM_TIME_t), INTENT(INOUT) :: IMT
TYPE(ISBA_PARAM_MEB_t), INTENT(INOUT) :: IMM
TYPE(ISBA_PARAM_IRRIG_t), INTENT(INOUT) :: IMI
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
TYPE(ISBA_INIT_t), INTENT(INOUT) :: INI
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
TYPE(GRID_t), INTENT(INOUT) :: G
TYPE(AGRI_t), INTENT(INOUT) :: AG
TYPE(DIAG_t), INTENT(INOUT) :: DGIP
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEIP
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DGMI
!
TYPE(DATE_TIME),      INTENT(IN)    :: TPTIME        ! current date and time
!
LOGICAL,              INTENT(IN)    :: OMEB          ! True = patch with multi-energy balance 
!                                                    ! False = patch with classical ISBA 
LOGICAL, DIMENSION(:),INTENT(INOUT) :: OSHADE        ! where vegetation evolution occurs
CHARACTER(LEN=*),     INTENT(IN)    :: HSNOW_ISBA    ! 'DEF' = Default F-R snow scheme
!                                                    !         (Douville et al. 1995)
!                                                    ! '3-L' = 3-L snow scheme (option)
!                                                    !         (Boone and Etchevers 2000)
CHARACTER(LEN=*),     INTENT(IN)    :: HIMPLICIT_WIND! wind implicitation option
!                                                    ! 'OLD' = direct
!                                                    ! 'NEW' = Taylor serie, order 1
REAL,                 INTENT(IN)    :: PTSTEP        ! Model time step (s)
!
REAL, DIMENSION(:),   INTENT(IN)    :: PPS           ! Pressure [Pa]
REAL, DIMENSION(:),   INTENT(IN)    :: PZENITH       ! solar zenith angle
REAL, DIMENSION(:),   INTENT(IN)    :: PSW_RAD       ! solar (shortwave) incoming radiation [W/m2]
REAL, DIMENSION(:),   INTENT(IN)    :: PLW_RAD       ! thermal (longwave) incoming radiation [W/m2]
REAL, DIMENSION(:),   INTENT(IN)    :: PSCA_SW       ! solar diffuse incoming radiation [W/m2]
REAL, DIMENSION(:),   INTENT(IN)    :: PEXNA         ! Exner function: forcing level (-)
REAL, DIMENSION(:),   INTENT(IN)    :: PEXNS         ! Exner function: surface (-)
REAL, DIMENSION(:),   INTENT(IN)    :: PRR           ! Rain rate (kg/m2/s)
REAL, DIMENSION(:),   INTENT(IN)    :: PSR           ! Snow rate (kg/m2/s)
REAL, DIMENSION(:),   INTENT(IN)    :: PRHOA         ! air density (kg/m3)
REAL, DIMENSION(:),   INTENT(IN)    :: PVMOD         ! modulus of the wind
!                                                    ! parallel to the orography (m/s)
REAL, DIMENSION(:),   INTENT(IN)    :: PTA           ! Temperature of atmosphere (K)
REAL, DIMENSION(:),   INTENT(IN)    :: PQA           ! specific humidity of atmosphere (kg/kg)
REAL, DIMENSION(:),   INTENT(IN)    :: PZREF         ! normal distance of the first
!                                                    ! atmospheric level to the
!                                                    ! orography (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PUREF         ! reference height of the wind (m)
!                                                    ! NOTE this is different from ZZREF
!                                                    ! ONLY in stand-alone/forced mode,
!                                                    ! NOT when coupled to a model (MesoNH)
REAL, DIMENSION(:),   INTENT(IN)    :: PDIRCOSZW     ! Director Cosinus along the z
!                                                    ! direction at the surface w-point
REAL, DIMENSION(:,:), INTENT(IN)    :: PSOILHCAPZ    ! ISBA-DF Soil heat capacity 
!                                                    ! profile [J/(m3 K)]
REAL, DIMENSION(:,:), INTENT(IN)    :: PSOILCONDZ    ! ISBA-DF Soil conductivity  
!                                                    ! profile  [W/(m K)]
REAL, DIMENSION(:),   INTENT(IN)    :: PFROZEN1      ! surface frozen fraction (-)
!
REAL, DIMENSION(:),   INTENT(IN)    :: PPALPHAN      ! snow/canopy transition coefficient
REAL, DIMENSION(:),   INTENT(IN)    :: PALBNIR_TVEG  ! albedo of vegetation in NIR 
!                                                    ! (needed for LM_TR or MEB)
REAL, DIMENSION(:),   INTENT(IN)    :: PALBVIS_TVEG  ! albedo of vegetation in VIS 
!                                                    ! (needed for LM_TR or MEB)
REAL, DIMENSION(:),   INTENT(IN)    :: PALBNIR_TSOIL ! albedo of bare soil in NIR 
!                                                    ! (needed for LM_TR or MEB)
REAL, DIMENSION(:),   INTENT(IN)    :: PALBVIS_TSOIL ! albedo of bare soil in VIS 
REAL, DIMENSION(:),   INTENT(IN)    :: PF2           ! Soil water stress factor for transpiration (-)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0G_WITHOUT_SNOW ! roughness length for momentum at snow-free canopy floor (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0_MEBV      ! roughness length for momentum over MEB vegetation part of patch (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0H_MEBV     ! roughness length for heat over MEB vegetation part of path (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0EFF_MEBV   ! roughness length for momentum over MEB vegetation part of patch (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0_MEBN      ! roughness length for momentum over MEB snow part of patch (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0H_MEBN     ! roughness length for heat over MEB snow part of path (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0EFF_MEBN   ! roughness length for momentum over MEB snow part of patch (m)
!
! implicit atmospheric coupling coefficients:
!
REAL, DIMENSION(:),   INTENT(IN)    :: PPET_A_COEF, PPET_B_COEF, &
                                       PPEQ_A_COEF, PPEQ_B_COEF, &
                                       PPEW_A_COEF, PPEW_B_COEF  
!                                                    ! PPEW_A_COEF  A-wind coefficient
!                                                    ! PPEW_B_COEF  B-wind coefficient
!                                                    ! PPET_A_COEF  A-air temperature coefficient
!                                                    ! PPET_B_COEF  B-air temperature coefficient
!                                                    ! PPEQ_A_COEF  A-air specific humidity coefficient
!                                                    ! PPEQ_B_COEF  B-air specific humidity coefficient
REAL, DIMENSION(:),   INTENT(IN)    :: PTDEEP_A          ! Deep soil temperature boundary condition 
!                                                         ! (prescribed)     
!                                      PTDEEP_A = Deep soil temperature
!                                                 coefficient depending on flux
!
! ISBA-Ags parameters
! (see also parameters with 'Ags:' in comments)
!
REAL, DIMENSION(:),   INTENT(IN) :: PCSP       ! atmospheric CO2 concentration
!                                                 [ppmm]=[kg CO2 / kg air]
REAL, DIMENSION(:),   INTENT(IN) :: PPOI       ! Gaussian weights (as above)
!                                              ! temperature
! - - - - - - - - - - - - - - - - - - - - 
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PABC          ! Ags: abscissa needed for integration
!                                                    ! of net assimilation and stomatal
!                                                    ! conductance over canopy depth
REAL, DIMENSION(:,:), INTENT(OUT)   :: PIACAN        ! PAR in the canopy at different gauss levels
!                                                    ! when using the DIF soil option (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PUSTAR        ! friction velocity
!
REAL, DIMENSION(:),   INTENT(OUT)   :: PGRNDFLUX     ! snow/soil-biomass interface flux (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PFLSN_COR     ! soil/snow interface correction flux to conserve energy (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PEMIST        ! total effective surface emissivity...LWUP = EMIST*TS_RAD**4 (-)
REAL, DIMENSION(:),   INTENT(OUT)   :: PAC_AGG       ! aggregated aerodynamic conductance
                                                     ! for evaporative flux calculations
REAL, DIMENSION(:),   INTENT(OUT)   :: PHU_AGG       ! aggregated relative humidity
                                                     ! for evaporative flux calculations
REAL, DIMENSION(:),   INTENT(OUT)   :: PDELHEATV_SFC ! change in heat storage of the vegetation canopy layer over the current time step (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PDELHEATG_SFC ! change in heat storage of the ground sfc layer over the current time step (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PDELHEATG     ! change in heat storage of the entire soil column over the current time step (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PRESTOREN     ! conductive heat flux between the surface and sub-surface soil layers 
!                                                    ! for the multi-layer snow schemes..for composite snow, it is 
!                                                    ! equal to DGEIP%XRESTORE (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PDELHEATN     ! change in heat storage of the entire snow column over the current time step (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PDELHEATN_SFC ! change in heat storage of the surface snow layer over the current time step (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PDEEP_FLUX    ! Heat flux at bottom of ISBA (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PRISNOW       ! Richarson number over ground-based snowpack (-)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOW_THRUFAL ! rate that liquid water leaves (explicit) snow pack: 
!                                                    ! ISBA-ES or CROCUS [kg/(m2 s)]
REAL, DIMENSION(:),   INTENT(OUT)   :: PEVAPCOR      !  evaporation correction as last traces of snow
!                                                    ! cover ablate..if sublimation exceeds trace amounts
                                                     ! of snow during time step, required residual mass taken 
                                                     ! from sfc soil layer [kg/(m2 s)]
REAL, DIMENSION(:),   INTENT(OUT)   :: PSUBVCOR      ! A possible snow mass correction (to be potentially    
!                                                    !  removed from soil)  (kg/m2/s)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSNOWSFCH     ! snow surface layer pseudo-heating term owing to
!                                                    !  changes in grid thickness            (W/m2)
REAL, DIMENSION(:),   INTENT(OUT)   :: PQSNOW        ! snow surface specific humidity (kg/kg)
!
! diagnostic variables for Carbon assimilation:
!
REAL, DIMENSION(:,:), INTENT(OUT)   :: PRESP_BIOMASS_INST ! instantaneous biomass respiration (kgCO2/kgair m/s)
!
!*      0.2    declarations of local variables
!
!
REAL, PARAMETER                                    :: ZTSTEP_EB     = 300. ! s Minimum time tstep required 
!                                                                          !   to time-split MEB energy budget
INTEGER                                            :: JTSPLIT_EB           ! number of time splits
INTEGER                                            :: JDT                  ! time split loop index
!
REAL                                               :: ZTSTEP               ! Local time split timestep (s)
REAL, DIMENSION(SIZE(PPS))                         :: ZWORK,ZWORK2,ZWORK3,ZWORK4  ! Working variables [*]
REAL, DIMENSION(SIZE(IR%TSNOW%WSNOW,1),SIZE(IR%TSNOW%WSNOW,2)) :: ZSNOWCOND            ! snow thermal conductivity  [W/(m K)] 
REAL, DIMENSION(SIZE(IR%TSNOW%WSNOW,1),SIZE(IR%TSNOW%WSNOW,2)) :: ZSNOWHCAP            ! snow heat capacity [J/(m3 K)]
REAL, DIMENSION(SIZE(IR%TSNOW%WSNOW,1),SIZE(IR%TSNOW%WSNOW,2)) :: ZSNOWRHO             ! snow layer density (kg/m3)
REAL, DIMENSION(SIZE(IR%TSNOW%WSNOW,1),SIZE(IR%TSNOW%WSNOW,2)) :: ZSNOWAGE             ! snow layer grain age
REAL, DIMENSION(SIZE(IR%TSNOW%WSNOW,1),SIZE(IR%TSNOW%WSNOW,2)) :: ZSNOWSWE             ! snow layer liquid water equivalent (kg/m2)
REAL, DIMENSION(SIZE(IR%TSNOW%WSNOW,1),SIZE(IR%TSNOW%WSNOW,2)) :: ZTAU_N               ! snow rad transmission coef at layer base (-)
REAL, DIMENSION(SIZE(PPS))                         :: ZCHIP                ! 
REAL, DIMENSION(SIZE(PPS))                         :: ZALBG                ! Effective ground albedo
REAL, DIMENSION(SIZE(PPS))                         :: ZSIGMA_F             ! LW transmission factor
REAL, DIMENSION(SIZE(PPS))                         :: ZSIGMA_FN            ! LW transmission factor - including buried (snow) 
!                                                                          ! vegetation effect
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_V_DTV        ! LW Jacobian: flux derrivative d LWnet_v/dTv [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_V_DTG        ! LW Jacobian: flux derrivative d LWnet_v/dTg [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_V_DTN        ! LW Jacobian: flux derrivative d LWnet_v/dTn [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_G_DTV        ! LW Jacobian: flux derrivative d LWnet_g/dTv [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_G_DTG        ! LW Jacobian: flux derrivative d LWnet_g/dTg [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_G_DTN        ! LW Jacobian: flux derrivative d LWnet_g/dTn [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_N_DTV        ! LW Jacobian: flux derrivative d LWnet_n/dTv [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_N_DTG        ! LW Jacobian: flux derrivative d LWnet_n/dTg [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_N_DTN        ! LW Jacobian: flux derrivative d LWnet_n/dTn [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZWRMAX               ! maximum canopy water equivalent interception capacity  [kg/m2]
REAL, DIMENSION(SIZE(PPS))                         :: ZWRLMAX              ! maximum litter water equivalent interception capacity  [kg/m2]
REAL, DIMENSION(SIZE(PPS))                         :: ZRS                  ! stomatal resistance (s/m)
REAL, DIMENSION(SIZE(PPS))                         :: ZRSN                 ! stomatal resistance of non-snow-buried canopy (s/m)
!                                                                          ! Etv=>0 as F2=>0 (-)  
REAL, DIMENSION(SIZE(PPS))                         :: ZWRVNMAX             ! maximum snow water equivalent interception capacity (kg/m2)
REAL, DIMENSION(SIZE(PPS))                         :: ZPSNCV               ! intercepted canopy snow fraction (-) NOTE! Not the same as the
!                                                                          ! ground-based snowpack
REAL, DIMENSION(SIZE(PPS))                         :: ZMELTVN              ! intercepted canopy snow net freeze/melt rate (kg/m2/s)
!                                                                          ! (if it is < 0, this signifies freezing)
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMA_TA            ! linear transform energy budget coefficient for Ta
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMB_TA            ! linear transform energy budget coefficient for Ta
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMA_TC            ! linear transform energy budget coefficient for Tc
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMB_TC            ! linear transform energy budget coefficient for Tc
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMA_TN            ! linear transform energy budget coefficient for Tn
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMB_TN            ! linear transform energy budget coefficient for Tn
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMA_TG            ! linear transform energy budget coefficient for Tg
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMB_TG            ! linear transform energy budget coefficient for Tg
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMA_TV            ! linear transform energy budget coefficient for Tv
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMB_TV            ! linear transform energy budget coefficient for Tv
REAL, DIMENSION(SIZE(PPS))                         :: ZPET_A_COEF          ! atmospheric coupling coefficient: Ta
REAL, DIMENSION(SIZE(PPS))                         :: ZPET_B_COEF          ! atmospheric coupling coefficient: Ta
REAL, DIMENSION(SIZE(PPS))                         :: ZKVN                 ! snow interception efficiency
REAL, DIMENSION(SIZE(PPS))                         :: ZVELC                ! wind speed at the top of the canopy (m/s)
REAL, DIMENSION(SIZE(PPS))                         :: ZDELTA               ! fraction of the foliage
!                                                                          ! covered with intercepted water (-)
REAL, DIMENSION(SIZE(PPS))                         :: ZHUGI                ! humidity over frozen bare ground (-)
REAL, DIMENSION(SIZE(PPS))                         :: ZHVN                 ! Halstead coefficient vegetation canopy above snow (-) 
REAL, DIMENSION(SIZE(PPS))                         :: ZHVG                 ! Halstead coefficient vegetation canopy above snow-free ground (-) 
REAL, DIMENSION(SIZE(PPS))                         :: ZLEG_DELTA           ! soil evaporation delta fn (-) 
REAL, DIMENSION(SIZE(PPS))                         :: ZLEGI_DELTA          ! soil sublimation delta fn (-) 
REAL, DIMENSION(SIZE(PPS))                         :: ZHSGL                ! surface halstead cofficient for bare soil (-)
REAL, DIMENSION(SIZE(PPS))                         :: ZHSGF                ! surface halstead cofficient for bare soil ice  (-)
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_C_A            ! turb transfer coef between vegetation canopy air and atmosphere (kg/m2/s) 
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_N_A            ! ...between the snow on the ground and atmosphere    (kg/m2/s) 
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_G_C            ! ...between snow-free ground and canopy air     (kg/m2/s)    
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_N_C            ! ...between snow on the ground and canopy air   (kg/m2/s)     
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_VG_C           ! ...between vegetation canopy over snow-free ground and canopy air   (kg/m2/s) 
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_VN_C           ! ...between vegetation canopy over the snow on the ground and canopy air  (kg/m2/s)  
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_V_C            ! ...between vegetation canopy and canopy air  (kg/m2/s)               
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_MOM            ! Effective drag coefficient for momentum [kg/(m2 s)]    
REAL, DIMENSION(SIZE(PPS))                         :: ZQSATG               ! saturation specific humidity for IR%XTG (ground surface: kg kg-1)    
REAL, DIMENSION(SIZE(PPS))                         :: ZQSATV               ! saturation specific humidity for IR%XTV (vegetation canopy: kg kg-1) 
REAL, DIMENSION(SIZE(PPS))                         :: ZQSATC               ! saturation specific humidity for IR%XTC (canopy air: kg kg-1)      
REAL, DIMENSION(SIZE(PPS))                         :: ZQSATN               ! saturation specific humidity for DGMI%XSNOWTEMP (snow surface: kg kg-1) 
REAL, DIMENSION(SIZE(PPS))                         :: ZDELTAVK             ! canopy interception capacity fraction including K-factor (-)  
REAL, DIMENSION(SIZE(PPS))                         :: ZCHEATV              ! Vegetation canopy *effective surface* heat capacity    (J m-2 K-1)
REAL, DIMENSION(SIZE(PPS))                         :: ZCHEATG              ! Understory-ground *effective surface* heat capacity    (J m-2 K-1)
REAL, DIMENSION(SIZE(PPS))                         :: ZCHEATN              ! Ground-based snow *effective surface* heat capacity    (J m-2 K-1)
REAL, DIMENSION(SIZE(PPS))                         :: ZHVGS                ! Dimensionless pseudo humidity factor for computing 
!                                                                          !  vapor fluxes from the non-buried part of the canopy 
!                                                                          !  to the canopy air                                     (-)
REAL, DIMENSION(SIZE(PPS))                         :: ZHVNS                ! Dimensionless pseudo humidity factor for computing 
!                                                                          !  vapor fluxes from the partly-buried part of the canopy 
!                                                                          !  to the canopy air                                     (-)
REAL, DIMENSION(SIZE(PPS))                         :: ZDQSAT_G             ! saturation specific humidity derivative for understory (kg kg-1 K-1)
REAL, DIMENSION(SIZE(PPS))                         :: ZDQSAT_V             ! saturation specific humidity derivative for the  
!                                                                          !  vegetation canopy                                     (kg kg-1 K-1)
REAL, DIMENSION(SIZE(PPS))                         :: ZDQSATI_N            ! saturation specific humidity derivative over ice for 
!                                                                          !  the ground-based snowpack                             (kg kg-1 K-1)
REAL, DIMENSION(SIZE(PPS))                         :: ZDELTAT_G            ! Time change in soil surface temperature                (K)
REAL, DIMENSION(SIZE(PPS))                         :: ZDELTAT_V            ! Time change in vegetation canopy temperature           (K)
REAL, DIMENSION(SIZE(PPS))                         :: ZDELTAT_N            ! Time change in snowpack surface temperature            (K)
REAL, DIMENSION(SIZE(PPS))                         :: ZRNET_V              ! Net vegetation canopy radiation                        (W/m2)
REAL, DIMENSION(SIZE(PPS))                         :: ZRNET_G              ! Net understory-ground radiation                        (W/m2)
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_C_A_F          ! Exchange coefficient between the snow on the ground and 
!                                                                          !  atmosphere modified by a partially to fully buried 
!                                                                          !  vegetation canopy                                     [kg/(m2 s)]
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_N_A_F          ! Exchange coefficient between vegetation canopy air and 
!                                                                          !  atmosphere modified by a partially to fully buried 
!                                                                          !  vegetation canopy                                     [kg/(m2 s)]
REAL, DIMENSION(SIZE(PPS))                         :: ZEVAP_C_A            ! Total canopy evapotranspiration and sublimation
!                                                                          !  of intercepted snow                                    (kg/m2/s)
REAL, DIMENSION(SIZE(PPS))                         :: ZEVAP_N_A            ! Vapor flux from the ground-based snowpack (part burying 
!                                                                          !  the canopy vegetation) to the atmosphere              (kg/m2/s)
REAL, DIMENSION(SIZE(PPS))                         :: ZH_N_A               ! Sensible heat flux from the ground-based snowpack (part 
!                                                                          !  burying the canopy vegetation) to the atmosphere      (W/m2)
REAL, DIMENSION(SIZE(PPS))                         :: ZVEGFACT             ! Fraction of canopy vegetation possibly receiving 
!                                                                          !  rainfall                                              (-)
REAL, DIMENSION(SIZE(PPS))                         :: ZRRSFC               ! The sum of all non-intercepted rain and canopy drip    (kg/m2/s)
REAL, DIMENSION(SIZE(PPS))                         :: ZRRSFCL              ! The sum of all non-intercepted rain and litter drip    (kg/m2/s)
REAL, DIMENSION(SIZE(PPS))                         :: ZLES3L               ! latent heat flux - sublimation of ice from the ground 
!                                                                          !  based snowpack (W/m2)
REAL, DIMENSION(SIZE(PPS))                         :: ZLEL3L               ! latent heat flux - evaporation of liquid water from the 
!                                                                          !  ground based snowpack (W/m2))
REAL, DIMENSION(SIZE(PPS))                         :: ZEVAP3L              ! total mass loss via evap & sublm from the ground based snowpack (kg/m2/s)
REAL, DIMENSION(SIZE(PPS))                         :: ZUSTAR2_IC           ! friction velocity (possibly implicitly coupled) (m/s)
REAL, DIMENSION(SIZE(PPS))                         :: ZTA_IC               ! atmospheric temperature (possibly implicitly coupled) (m/s)
REAL, DIMENSION(SIZE(PPS))                         :: ZQA_IC               ! atmospheric specific humidity (possibly implicitly coupled) (m/s)
REAL, DIMENSION(SIZE(PPS))                         :: ZSWUP                ! net upwelling shortwave radiation [W/m2]
REAL, DIMENSION(SIZE(PPS))                         :: ZLWUP                ! net upwelling longwave radiation [W/m2]
REAL, DIMENSION(SIZE(PPS))                         :: ZUSTAR2SNOW          ! snow fraciton velocity squared (m2/s2)
REAL, DIMENSION(SIZE(PPS))                         :: ZVMOD                ! lowest level atmospheric wind speed update estimate (K)
REAL, DIMENSION(SIZE(PPS))                         :: ZRR                  ! combined rain rate (above canopy) and irrigation need (kg/m2/s)
REAL, DIMENSION(SIZE(PPS))                         :: ZFLSN_COR            ! snow/soil-biomass correction flux (W/m2) (not MEB)
REAL, DIMENSION(SIZE(PPS))                         :: ZWSFC                ! surface liquid water content for resistances  (m3/m3)
REAL, DIMENSION(SIZE(PPS))                         :: ZWISFC               ! surface frozen water content for resistances  (m3/m3)
REAL, DIMENSION(SIZE(PPS))                         :: ZLESFC               ! evaporation from the surface (soil or litter) (W/m2)
REAL, DIMENSION(SIZE(PPS))                         :: ZLESFCI              ! sublimation from the surface (soil or litter) (W/m2)
REAL, DIMENSION(SIZE(PPS))                         :: ZPERMSNOWFRAC        ! fraction of permanent snow/ice
!
! - TR_ML radiation option: NOTE...always used by MEB
!
REAL, DIMENSION(SIZE(PPS),SIZE(PABC))              :: ZIACAN_SUNLIT        ! Absorbed PAR of each level within the
REAL, DIMENSION(SIZE(PPS),SIZE(PABC))              :: ZIACAN_SHADE         !    canopy - Split into SHADEd and SUNLIT
REAL, DIMENSION(SIZE(PPS),SIZE(PABC))              :: ZFRAC_SUN            !    fraction of sunlit leaves
!
REAL, DIMENSION(SIZE(PPS))                         :: ZLAI                 ! Potentially covered/buried canopy LAI (m2/m2)
REAL, DIMENSION(SIZE(PPS))                         :: ZALBVIS_TSOIL        ! average snow-free ground VIS albedo (soil plus flooded fraction) 
REAL, DIMENSION(SIZE(PPS))                         :: ZALBNIR_TSOIL        ! average snow-free ground NIR albedo (soil plus flooded fraction)
REAL, DIMENSION(SIZE(PPS))                         :: ZSWNET_S             ! Net SW radiation at the surface (below canopy snow/ground/flooded zone)
!
!
! - CPHOTO/=NON (Ags Option(s)):
!
REAL, DIMENSION(SIZE(PPS))                         :: ZQSAT                ! CPHOTO/=NON (Ags Option(s))diagnosed (past time step) Qsat relative to canopy (for Ags)
REAL, DIMENSION(SIZE(PPS))                         :: ZFFV                 ! submerged vegetation (by flooding) fraction (-)
REAL, DIMENSION(SIZE(PPS),SIZE(PABC))              :: ZIACAN               ! PAR in the canopy at different gauss levels: local working needed if
!                                                                          ! Ags if off (i.e. CPHOTO==NON)
!
REAL, DIMENSION(:,:), ALLOCATABLE                  :: ZTGL                 ! Temporary temperature of litter + soil
REAL, DIMENSION(:,:), ALLOCATABLE                  :: ZSOILHCAPZ           ! Temporary heat capacity of litter + soil
REAL, DIMENSION(:,:), ALLOCATABLE                  :: ZSOILCONDZ           ! Temporary heat conductivity of litter + soil
REAL, DIMENSION(:,:), ALLOCATABLE                  :: ZD_G                 ! Temporary depth of bottom litter + soil layers
REAL, DIMENSION(:,:), ALLOCATABLE                  :: ZDZG                 ! Temporary thickness of litter + soil layers
REAL, DIMENSION(:,:), ALLOCATABLE                  :: ZWFC                 ! Temporary Wfc of bottom litter + soil layers
REAL, DIMENSION(:,:), ALLOCATABLE                  :: ZWSAT                ! Temporary Wsat of bottom litter + soil layers
!
! Working sums for flux averaging over MEB time split
!
REAL, DIMENSION(SIZE(PPS))   :: ZH_SUM, ZH_C_A_SUM, ZH_N_A_SUM, ZH_V_C_SUM, ZH_G_C_SUM, &
                                ZH_N_C_SUM, ZHSNOW_SUM, ZHPSNOW_SUM
REAL, DIMENSION(SIZE(PPS))   :: ZHU_AGG_SUM, ZAC_AGG_SUM

REAL, DIMENSION(SIZE(PPS))   :: ZLE_SUM, ZLE_C_A_SUM, ZLE_V_C_SUM, ZLE_G_C_SUM,           &
                                ZLE_N_C_SUM, ZLETR_V_C_SUM, ZLEG_SUM,ZLEGI_SUM,ZLESFC_SUM,&
                                ZLESFCI_SUM,                                              &
                                ZLER_V_C_SUM, ZLE_FLOOD_SUM, ZLEI_FLOOD_SUM,              &
                                ZLES_V_C_SUM, ZLETR_SUM, ZLER_SUM, ZLEV_SUM,              &
                                ZLEI_SUM, ZLES3L_SUM, ZLEL3L_SUM, ZEVAP3L_SUM,            &
                                ZUSTAR2_SUM, ZUSTAR2SNOW_SUM, ZCDSNOW_SUM,                &
                                ZCHSNOW_SUM, ZRISNOW_SUM

REAL, DIMENSION(SIZE(PPS))   :: ZGRNDFLUX_SUM, ZRESTORE_SUM

REAL, DIMENSION(SIZE(PPS))   :: ZSWNET_V_SUM, ZSWNET_G_SUM, ZSWNET_N_SUM, ZLWNET_V_SUM, &
                                ZLWNET_G_SUM, ZLWNET_N_SUM, ZEMIST_SUM, ZSWUP_SUM,      &
                                ZLWUP_SUM
REAL, DIMENSION(SIZE(PPS))   :: ZDELHEATG_SFC_SUM, ZDELHEATV_SFC_SUM, ZDELHEATG_SUM
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
INTEGER :: INJ, INL, JJ, JL
REAL, DIMENSION(SIZE(IR%XWR,1))         :: ZPHASEL  ! Phase changement in litter (W/m2)
REAL, DIMENSION(SIZE(IR%XWR,1))         :: ZCTSFC 
!-------------------------------------------------------------------------------
!
!*      1.0    Preliminaries
!              -------------
!
IF (LHOOK) CALL DR_HOOK('ISBA_MEB',0,ZHOOK_HANDLE)
!
!
PIACAN(:,:)        = 0.
DGMI%XFAPAR(:)          = 0.
DGMI%XFAPIR(:)          = 0.
DGMI%XFAPAR_BS(:)       = 0.
DGMI%XFAPIR_BS(:)       = 0.
DGEIP%XRRLIT(:)          =0.0
DGEIP%XDRIPLIT(:)        =0.0
!
DGEIP%XLEGI(:)  = 0.
DGEIP%XLEG(:)   = 0.
ZLESFCI(:)= 0.
ZLESFC(:) = 0.
!
ZIACAN_SUNLIT(:,:) = XUNDEF
ZIACAN_SHADE(:,:)  = XUNDEF
ZFRAC_SUN (:,:)    = XUNDEF
ZLAI (:)           = XUNDEF
ZALBVIS_TSOIL(:)   = XUNDEF
ZALBNIR_TSOIL(:)   = XUNDEF
ZSWNET_S(:)        = XUNDEF
ZQSAT(:)           = XUNDEF
!
!*      1.1    Preliminaries for litter parameters
!              -----------------------------------
!
INJ=SIZE(IR%XWG,1)
INL=SIZE(IR%XWG,2)
!
CALL ALLOCATE_LOCAL_VARS_PREP_GRID_SOIL
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      1.2    Preliminaries for litter temperature
!              ------------------------------------
!
! Concatenate IR%XTL and IR%XTG and the parameters linked to heat transfer into the soil
!

CALL PREP_MEB_SOIL(IO%LMEB_LITTER, PSOILHCAPZ, PSOILCONDZ, IMX%XDG(:,:,1), IP, IR, &
                   IMM%XGNDLITTER(:,1) ,ZD_G, ZDZG,ZTGL, ZSOILHCAPZ, ZSOILCONDZ,&
                   ZWSAT, ZWFC, ZWSFC, ZWISFC, ZCTSFC, DGMI%XCT      )
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      2.0    Preliminaries for energy and radiation budget
!              ---------------------------------------------
!
ZPERMSNOWFRAC(:) = IP%XVEGTYPE_PATCH(:,NVT_SNOW,1)
!
! Local working:
! - possibly adjust these prognostic variables locally, but do not save
!
ZSNOWRHO(:,:)    = IR%TSNOW%RHO(:,:,1)
ZSNOWAGE(:,:)    = IR%TSNOW%AGE(:,:,1)
ZSNOWSWE(:,:)    = IR%TSNOW%WSNOW(:,:,1)
!
CALL PREPS_FOR_MEB_EBUD_RAD(PPS, IMT%XLAI(:,1), ZSNOWRHO, ZSNOWSWE, IR%TSNOW%HEAT(:,:,1), &
                            DGMI%XSNOWTEMP, DGMI%XSNOWDZ, ZSNOWCOND, ZSNOWHCAP, IR%TSNOW%EMIS(:,1), &
                            ZSIGMA_F, ZCHIP, PTSTEP, PSR, PTA, PVMOD, ZSNOWAGE, ZPERMSNOWFRAC  )
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      3.0    Shortwave radiative transfer
!              ----------------------------  
!
! Calculate snow albedo: split into spectral bands:
!
CALL SNOWALB_SPECTRAL_BANDS_MEB(IP%XVEGTYPE_PATCH(:,:,1), IR, ZSNOWRHO, ZSNOWAGE, PPS,&
                                DGMI%XSNOWDZ,PZENITH, ZTAU_N)
!
!
! NOTE, currently MEB only uses 2 of 3 potential snow albedo spectral bands
!
!
WHERE(IR%TSNOW%ALB(:,1) /= XUNDEF)
   ZLAI(:)          = IMT%XLAI(:,1)*(1.0-PPALPHAN(:))
   ZALBVIS_TSOIL(:) = PALBVIS_TSOIL(:)*(1.-IR%XPSN(:,1)) + IR%XPSN(:,1)*IR%TSNOW%ALBVIS(:,1)
   ZALBNIR_TSOIL(:) = PALBNIR_TSOIL(:)*(1.-IR%XPSN(:,1)) + IR%XPSN(:,1)*IR%TSNOW%ALBNIR(:,1)
ELSEWHERE
   ZLAI(:)          = IMT%XLAI(:,1)
   ZALBVIS_TSOIL(:) = PALBVIS_TSOIL(:)
   ZALBNIR_TSOIL(:) = PALBNIR_TSOIL(:)
END WHERE
!
  CALL RADIATIVE_TRANSFERT(IO%LAGRI_TO_GRASS, IP%XVEGTYPE_PATCH(:,:,1), &
                           PALBVIS_TVEG, ZALBVIS_TSOIL, PALBNIR_TVEG, ZALBNIR_TSOIL,  &
                           PSW_RAD, ZLAI, PZENITH, PABC, IR%XFAPARC(:,1), IR%XFAPIRC(:,1), &
                           IR%XMUS(:,1), IR%XLAI_EFFC(:,1), OSHADE, ZIACAN,  ZIACAN_SUNLIT, &
                           ZIACAN_SHADE, ZFRAC_SUN, DGMI%XFAPAR, DGMI%XFAPIR, DGMI%XFAPAR_BS, &
                           DGMI%XFAPIR_BS )    

! Total effective surface (canopy, ground/flooded zone, snow) all-wavelength
! albedo: diagnosed from shortwave energy budget closure

DGIP%XALBT(:)      = 1. - (XSW_WGHT_VIS*(DGMI%XFAPAR(:)+DGMI%XFAPAR_BS(:)) +     &
                           XSW_WGHT_NIR*(DGMI%XFAPIR(:)+DGMI%XFAPIR_BS(:)))
ZSWUP(:)      = PSW_RAD(:)*DGIP%XALBT(:)
DGIP%XALBT(:)      = ZSWUP(:)/MAX(1.E-5, PSW_RAD(:))

! Diagnose all-wavelength SW radiative budget components:

DGEIP%XSWNET_V(:)   = PSW_RAD(:)*(XSW_WGHT_VIS*DGMI%XFAPAR(:)    +         &
                                 XSW_WGHT_NIR*DGMI%XFAPIR(:)   )
ZSWNET_S(:)   = PSW_RAD(:)*(XSW_WGHT_VIS*DGMI%XFAPAR_BS(:) +                   &
                            XSW_WGHT_NIR*DGMI%XFAPIR_BS(:))
DGEIP%XSWNET_N(:)   = ZSWNET_S(:)*    IR%XPSN(:,1)
DGEIP%XSWNET_G(:)   = ZSWNET_S(:)*(1.-IR%XPSN(:,1))

! Quantity of net shortwave radiation absorbed in surface snow layer 

DGEIP%XSWNET_NS(:)  = DGEIP%XSWNET_N(:)*(1.0 - ZTAU_N(:,1))

! Compute all-wavelength effective ground albedo

ZALBG(:)      = XSW_WGHT_NIR*ZALBNIR_TSOIL(:) +                           & 
                XSW_WGHT_VIS*ZALBVIS_TSOIL(:)

! Any SW radiation reaching the base of the lowest snow layer can pass
! into the soil:

ZTAU_N(:,SIZE(IR%TSNOW%WSNOW,2)) = ZTAU_N(:,SIZE(IR%TSNOW%WSNOW,2))*(1.-ZALBG(:))

! Downwelling SW radiation arriving at ground/snow surface

DGEIP%XSWDOWN_GN(:) = ZSWNET_S(:)/(1.-ZALBG(:))
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      4.0    Longwave radiative transfer
!              ---------------------------  
!
CALL ISBA_LWNET_MEB(IMT%XLAI(:,1), IR%XPSN(:,1), PPALPHAN,IR%TSNOW%EMIS(:,1),&
                    INI%XEMISF(:,1), INI%XFF(:,1), IR%XTV(:,1), ZTGL(:,1), &
                    DGMI%XSNOWTEMP(:,1),  PLW_RAD, DGEIP%XLWNET_N, DGEIP%XLWNET_V, &
                    DGEIP%XLWNET_G, ZDLWNET_V_DTV, ZDLWNET_V_DTG, ZDLWNET_V_DTN,   &
                    ZDLWNET_G_DTV, ZDLWNET_G_DTG, ZDLWNET_G_DTN, ZDLWNET_N_DTV, &
                    ZDLWNET_N_DTG, ZDLWNET_N_DTN, ZSIGMA_F, ZSIGMA_FN, DGEIP%XLWDOWN_GN )
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      5.0    Fraction of leaves occupied by intercepted water
!              ------------------------------------------------
!
! Vegetation canopy:
!
! First, compute an effective veg fraction: it can only be < unity if vegetation is buried by snowpack...
!
ZWORK(:) = (1.0 - IR%XPSN(:,1) + IR%XPSN(:,1)*(1.0 - PPALPHAN(:))) 
! 
CALL WET_LEAVES_FRAC(IR%XWR(:,1), ZWORK, IMT%XWRMAX_CF(:,1), PZ0_MEBV, IMT%XLAI(:,1), ZWRMAX, ZDELTA) 
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      6.0    Plant stress, stomatal resistance and, possibly, CO2 assimilatio
!              --------------------------------------------------------------------
!
!              MEB-NOTE here assumed IO%CPHOTO=='DEF' or 'AST' for now
!              More Ags options to be added later 
!
IF (IO%CPHOTO=='NON') THEN
!
! Canopy vegetation (no snow, or snow below the main part of the canopy):
!
   CALL VEG(PSW_RAD, IR%XTC(:,1), IR%XQC(:,1), PPS, IMT%XRGL(:,1), &
            IMT%XLAI(:,1), IMT%XRSMIN(:,1), IMT%XGAMMA(:,1), PF2, ZRS)
!
!
ELSE IF (MAXVAL(IMT%XGMES) /= XUNDEF .OR. MINVAL(IMT%XGMES) /= XUNDEF) THEN
!
! NOTE: For now we assume that forest canopy can be flooded.
! However, we need to likely compute a fraction like PALPHAN (for snow vertical extent)
! for floods for grasses/crops/shrubs...i.e. low vegetation

   ZFFV(:)  = 0.0

   ZQSAT(:) = QSAT(IR%XTV(:,1),PPS)  
   CALL COTWORES(PTSTEP, IO, OSHADE,  IP, IMT, IR,  IMX%XDMAX(:,1),  &
                 PPOI, PCSP, IR%XTV(:,1), PF2, PSW_RAD, IR%XQC(:,1), &
                 ZQSAT, PPALPHAN, ZDELTA, PRHOA, PZENITH, ZFFV,      &
                 ZIACAN_SUNLIT, ZIACAN_SHADE, ZFRAC_SUN, ZIACAN,     &
                 PABC, ZRS, DGEIP%XGPP, PRESP_BIOMASS_INST(:,1))
!
   PIACAN(:,:)             = ZIACAN(:,:)
!
ELSE
   PRESP_BIOMASS_INST(:,1) = 0.0
   DGEIP%XGPP(:)                 = 0.0
ENDIF
!
! Additional resistance for possibly snow-buried canopy vegetation:
!
ZRSN(:) = ZRS(:)/( 1.0 - MIN(PPALPHAN(:), 1.0 - (ZRS(:)/XRS_MAX)) ) 
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      6.0    Canopy snow (intercepted) needed diagnostics:
!              ---------------------------------------------
!
CALL SNOW_LEAVES_FRAC_MEB(IR%XPSN(:,1), PPALPHAN, IR%XWRVN(:,1), IR%XTV(:,1), ZCHIP, &
                          IMT%XLAI(:,1), ZWRVNMAX, ZPSNCV, ZMELTVN)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      7.0    Aerodynamic drag and heat/mass transfer/fluxes 
!              and energy budget solution
!              ----------------------------------------------
!
! NOTE, this assumes thermodynamic variable herein is potential T

ZPET_A_COEF(:)  = -PPET_A_COEF(:)*XCPD 
ZPET_B_COEF(:)  =  PPET_B_COEF(:)*XCPD
ZTHRMA_TA(:)    =  XCPD/PEXNA(:)
ZTHRMB_TA(:)    =  0.0
ZWORK(:)        =  XCPD/PEXNS(:)
ZTHRMA_TC(:)    =  ZWORK(:)
ZTHRMB_TC(:)    =  0.0
ZTHRMA_TN(:)    =  ZWORK(:)
ZTHRMB_TN(:)    =  0.0
ZTHRMA_TG(:)    =  ZWORK(:)
ZTHRMB_TG(:)    =  0.0
ZTHRMA_TV(:)    =  ZWORK(:)
ZTHRMB_TV(:)    =  0.0
!
!
! Possibly split time step if large: 
! Although the energy budget is fully implicit, a very small canopy heat capacity 
! (and neglect of canopy air space heat capacity) can possibly lead to
! numerical shocks, especially during transition periods between stable and unstable 
! regimes. Thus, for relatively large steps, a simple time split scheme is activated.
! Note that soil moisture is held constant, while turbulent exchange coefficients are updated during the split.
! Also, experience shows that splitting at least once for moderately sized time steps 
! is quite effective in removing any lingering small but possible oscillations.
! Finally, for *very* small time steps (such as those for high res runs), no split is performed.
! Fluxes are averaged over the time split for conservation.
!
JTSPLIT_EB      = 1 + INT(PTSTEP/ZTSTEP_EB)  ! number of split-time steps
ZTSTEP          = PTSTEP/JTSPLIT_EB          ! split time step...for relatively small time steps, no split
!
! initialize time split sums for fluxes:
!
CALL INIT_SUM_FLUXES_MEB_TSPLIT 
!
!
! Note, when implicitly coupled to the atmosphere, these
! 3 variables will evolve during the split...we used updated
! values for turbulent exchange computations (drag_meb). 
! NOTE...when explicit coupling used, these 3 variables do NOT vary
! during the split.
!
ZVMOD(:)  = PVMOD(:)
ZTA_IC(:) = PTA(:)
ZQA_IC(:) = PQA(:)
!
!
LOOP_TIME_SPLIT_EB: DO JDT=1,JTSPLIT_EB
!*      7.1    Aerodynamic drag and heat transfer coefficients
!              -----------------------------------------------
!
   CALL DRAG_MEB(IO, IR, DGMI, DGIP, IMT%XLAI(:,1), IMM%XH_VEG(:,1), &
                 ZTGL(:,1), ZTA_IC,  ZQA_IC, ZVMOD, ZWSFC, ZWISFC,   &
                 ZWSAT(:,1), ZWFC(:,1), PEXNS, PEXNA, PPS, PRR, PSR, &
                 PRHOA, PZ0G_WITHOUT_SNOW, PZ0_MEBV, PZ0H_MEBV,      &
                 PZ0EFF_MEBV, PZ0_MEBN, PZ0H_MEBN, PZ0EFF_MEBN,      &
                 ZSNOWSWE(:,1), ZCHIP, ZTSTEP, ZRS, ZRSN, PPALPHAN,  &
                 PZREF, PUREF, PDIRCOSZW, ZPSNCV, ZDELTA, ZVELC,     &
                 PRISNOW, ZUSTAR2SNOW, ZHUGI, ZHVG,                  &
                 ZHVN, ZLEG_DELTA, ZLEGI_DELTA, ZHSGL, ZHSGF,        &
                 ZFLXC_C_A, ZFLXC_N_A, ZFLXC_G_C, ZFLXC_N_C,         &
                 ZFLXC_VG_C, ZFLXC_VN_C, ZFLXC_MOM, ZQSATG, ZQSATV,  &
                 ZQSATC, ZQSATN, ZDELTAVK )
!
   ZKVN(:) = SNOW_INTERCEPT_EFF(ZCHIP,ZVELC,ZWRVNMAX)

!*      7.2    Resolution of the surface energy budgets
!              ----------------------------------------
!
   CALL E_BUDGET_MEB(IO, IP, INI, IR, DGIP, DGEIP, DGMI, IMT%XCV(:,1), &
                     ZTSTEP, PPS, ZCTSFC, PTDEEP_A, ZD_G, ZSOILCONDZ, ZSOILHCAPZ,&
                     ZSNOWCOND, ZSNOWHCAP, ZTAU_N, ZDLWNET_V_DTV, ZDLWNET_V_DTG, &
                     ZDLWNET_V_DTN, ZDLWNET_G_DTV, ZDLWNET_G_DTG, ZDLWNET_G_DTN, &
                     ZDLWNET_N_DTV, ZDLWNET_N_DTG, ZDLWNET_N_DTN, PPEW_A_COEF,   &
                     PPEW_B_COEF, ZPET_A_COEF, PPEQ_A_COEF, ZPET_B_COEF,         &
                     PPEQ_B_COEF, ZTHRMA_TA, ZTHRMB_TA, ZTHRMA_TC, ZTHRMB_TC,    &
                     ZTHRMA_TG, ZTHRMB_TG, ZTHRMA_TV, ZTHRMB_TV, ZTHRMA_TN,      &
                     ZTHRMB_TN, ZQSATG, ZQSATV, ZQSATN, PPALPHAN, ZPSNCV,        &
                     ZCHEATV, ZCHEATG, ZCHEATN, ZLEG_DELTA, ZLEGI_DELTA, ZHUGI,  &
                     ZHVG, ZHVN, PFROZEN1, ZFLXC_C_A, ZFLXC_G_C, ZFLXC_VG_C,     &
                     ZFLXC_VN_C, ZFLXC_N_C, ZFLXC_N_A, ZFLXC_MOM, ZTGL,          &
                     ZFLXC_V_C, ZHVGS, ZHVNS, ZDQSAT_G,ZDQSAT_V,ZDQSATI_N,       &
                     ZTA_IC, ZQA_IC, ZUSTAR2_IC, ZVMOD, ZDELTAT_G, ZDELTAT_V,    &
                     ZDELTAT_N, PGRNDFLUX, PDEEP_FLUX, PDELHEATV_SFC,            &
                     PDELHEATG_SFC, PDELHEATG                              )
!
!*      7.3    Energy and momentum fluxes and radiative temperature and emissivity
!              -------------------------------------------------------------------
!
   CALL ISBA_FLUXES_MEB(INI, IR, DGIP, DGEIP, DGMI, PRHOA, ZSIGMA_F, ZSIGMA_FN, &
                        ZRNET_V, ZRNET_G, ZDLWNET_V_DTV, ZDLWNET_V_DTG,         &
                        ZDLWNET_V_DTN, ZDLWNET_G_DTV, ZDLWNET_G_DTG,            &
                        ZDLWNET_G_DTN, ZDLWNET_N_DTV, ZDLWNET_N_DTG,            &
                        ZDLWNET_N_DTN, ZTHRMA_TA, ZTHRMB_TA, ZTHRMA_TC,         &
                        ZTHRMB_TC, ZTHRMA_TG, ZTHRMB_TG, ZTHRMA_TV, ZTHRMB_TV,  &
                        ZTHRMA_TN, ZTHRMB_TN,  ZQSATG, ZQSATV, ZQSATN, PPALPHAN,&
                        ZPSNCV, PFROZEN1, ZLEG_DELTA, ZLEGI_DELTA, ZHUGI, ZHVG, &
                        ZHVN, ZFLXC_C_A, ZFLXC_G_C, ZFLXC_VG_C, ZFLXC_VN_C,     &
                        ZFLXC_N_C, ZFLXC_N_A, ZFLXC_MOM, ZFLXC_V_C, ZHVGS,      &
                        ZHVNS, ZTGL, ZDQSAT_G, ZDQSAT_V, ZDQSATI_N, ZTA_IC,     &
                        ZQA_IC, ZDELTAVK, ZDELTAT_G, ZDELTAT_V, ZDELTAT_N,      &
                        ZSWUP, PSW_RAD, PLW_RAD, ZLWUP, ZH_N_A, ZEVAP_C_A,      &
                        ZEVAP_N_A, ZLESFC, ZLESFCI, ZLES3L, ZLEL3L, ZEVAP3L,    &
                        PEMIST                                   )
!
! Compute aggregated coefficients for evaporation
! Sum(LEC+LES+LEL) = ACagg * Lv * RHOA * (HUagg.Qsat - Qa)
!
   ZFLXC_C_A_F(:) = ZFLXC_C_A(:)*(1.0-IR%XPSN(:,1)*PPALPHAN(:))
   ZFLXC_N_A_F(:) = ZFLXC_N_A(:)*     IR%XPSN(:,1)*PPALPHAN(:)

   PHU_AGG(:)     = (ZFLXC_C_A_F(:)*IR%XQC(:,1)    + ZFLXC_N_A_F(:)*ZQSATN(:))/   &
                    (ZFLXC_C_A_F(:)*ZQSATC(:) + ZFLXC_N_A_F(:)*ZQSATN(:))

   PAC_AGG(:)     = ZFLXC_C_A_F(:) + ZFLXC_N_A_F(:) ! kg/m2/s
!
! Sum fluxes over time split:

   CALL SUM_FLUXES_MEB_TSPLIT  

ENDDO LOOP_TIME_SPLIT_EB
!
CALL AVG_FLUXES_MEB_TSPLIT     ! average fluxes over time split
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*     8.0    Snow explicit canopy loading/interception 
!             ------------------------------------------
!
CALL SNOW_LOAD_MEB(IR, DGEIP, PTSTEP, PSR, ZWRVNMAX, ZKVN, ZCHEATV, ZMELTVN, &
                   ZVELC, PSUBVCOR)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*     9.0    Snow explicit canopy loading/interception 
!             ------------------------------------------
!
ZRR(:)         = PRR(:)
DGEIP%XIRRIG_FLUX(:) = 0.0
!
!* add irrigation over vegetation to liquid precipitation (rr)
!  Water "need" treated as if sprayed from above (over vegetation and soil):
!
IF (SIZE(AG%LIRRIGATE,1)>0) THEN
   WHERE (AG%LIRRIGATE(:,1) .AND. IMI%XIRRIG(:,1)>0. .AND. IMI%XIRRIG(:,1) /= XUNDEF &
                            .AND. (PF2(:)<AG%XTHRESHOLDSPT(:,1)) )
      DGEIP%XIRRIG_FLUX(:) = IMI%XWATSUP(:,1) / XDAY           
      ZRR              (:) = PRR(:) + IMI%XWATSUP(:,1)/XDAY
      AG%LIRRIDAY   (:,1) = .TRUE.           
   END WHERE
ENDIF
!
! Call canopy interception...here because meltwater should be allowed to fall
! on understory snowpack
!
! Fraction of canopy vegetation possibly receiving rainfall/irrigation
!
ZVEGFACT(:) = ZSIGMA_F(:)*(1.0-PPALPHAN(:)*IR%XPSN(:,1)) 
!
! The sum of all non-intercepted rain and drip is "ZRRSFC" (kg/m2/s):
! this is then partitioned by snow scheme into part falling on
! snowpack and part falling onto snow-free understory.
!
!
CALL HYDRO_VEG(IO%CRAIN, PTSTEP, INI%XMUF, ZRR, DGEIP%XLEVCV, DGEIP%XLETRCV, &
               ZVEGFACT, ZPSNCV, IR%XWR(:,1), ZWRMAX, ZRRSFC, DGEIP%XDRIP, DGEIP%XRRVEG  )
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      10.0    Explicit snow scheme (MEB: impose surface fluxes as upper BC)
!              ----------------------------------------------------------------
!
 CALL SNOW3L_ISBA(IO, G, IR, DGIP, DGEIP, DGMI, HSNOW_ISBA, OMEB, HIMPLICIT_WIND,  &
                  TPTIME, PTSTEP, IP%XVEGTYPE_PATCH(:,:,1),  ZTGL, ZCTSFC,  &
                  ZSOILHCAPZ, ZSOILCONDZ(:,1), PPS, PTA, PSW_RAD, PQA,      &
                  PVMOD, PLW_RAD, ZRRSFC, DGEIP%XSR_GN, PRHOA, PUREF, PEXNS,&
                  PEXNA, PDIRCOSZW, PZREF, ZALBG, ZD_G, ZDZG, PPEW_A_COEF,  &
                  PPEW_B_COEF, PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF,       &
                  PPEQ_B_COEF, PSNOW_THRUFAL, PGRNDFLUX, PFLSN_COR,         &
                  PRESTOREN, PEVAPCOR, ZLES3L, ZLEL3L, ZEVAP3L, PSNOWSFCH,  &
                  PDELHEATN, PDELHEATN_SFC, PRISNOW, PZENITH, PDELHEATG,    &
                  PDELHEATG_SFC, PQSNOW     )  
!
! If a litter layer exists, compute hydrology:
!
IF(IO%LMEB_LITTER)THEN
!
   ZWORK(:)   = 0.
   ZWORK2(:)  = IR%XWRL(:,1)
   ZWORK3(:)  = 1.
   ZWORK4(:)  = PSNOW_THRUFAL(:) + DGMI%XRRSFC(:)
   ZWRLMAX(:) = IMM%XGNDLITTER(:,1)*ZWFC(:,1)*XRHOLW

   CALL HYDRO_VEG(IO%CRAIN, PTSTEP, INI%XMUF, ZWORK4(:), ZLESFC, ZWORK, &
                  ZWORK3, ZWORK, IR%XWRL(:,1) , ZWRLMAX, ZRRSFCL, DGEIP%XDRIPLIT, DGEIP%XRRLIT  )

ELSE

   ZRRSFCL(:) = ZRRSFC(:)

ENDIF
!
!*      11.0    Separate litter and soil temperature
!              ------------------------------------
!
CALL RESHIFT_MEB_SOIL(IO%LMEB_LITTER, ZTGL, ZLESFC, ZLESFCI, IR, DGEIP)              
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
CALL DEALLOCATE_LOCAL_VARS_PREP_GRID_SOIL
!
IF(IO%LMEB_LITTER)THEN
!
CALL ICE_LITTER(PTSTEP, DGEIP%XLELITTERI, PSOILHCAPZ, IR, IMX%NWG_LAYER(:,1), &
                IP%XDZG(:,:,1), IMM%XGNDLITTER(:,1), ZPHASEL,ZCTSFC   )
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('ISBA_MEB',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
CONTAINS
!
!===============================================================================
!
SUBROUTINE INIT_SUM_FLUXES_MEB_TSPLIT 
!
IMPLICIT NONE
!
!
!*      0.2    declarations of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ISBA_MEB:INIT_SUM_FLUXES_MEB_TSPLIT ',0,ZHOOK_HANDLE)
!
! sensible heat fluxes:
!
ZH_SUM(:)        = 0.0
ZH_C_A_SUM(:)    = 0.0
ZH_N_A_SUM(:)    = 0.0
ZH_V_C_SUM(:)    = 0.0
ZH_G_C_SUM(:)    = 0.0
ZH_N_C_SUM(:)    = 0.0
ZHSNOW_SUM(:)    = 0.0
!
! latent heat/water vapor fluxes:
!
ZLE_SUM(:)       = 0.0
ZLE_C_A_SUM(:)   = 0.0
ZLE_V_C_SUM(:)   = 0.0
ZLE_G_C_SUM(:)   = 0.0
ZLE_N_C_SUM(:)   = 0.0
ZLETR_V_C_SUM(:) = 0.0
ZLEG_SUM(:)      = 0.0
ZLEGI_SUM(:)     = 0.0
ZLESFC_SUM(:)    = 0.0
ZLESFCI_SUM(:)   = 0.0
ZLER_V_C_SUM(:)  = 0.0
ZLE_FLOOD_SUM(:) = 0.0
ZLEI_FLOOD_SUM(:)= 0.0
ZLES_V_C_SUM(:)  = 0.0
ZLETR_SUM(:)     = 0.0
ZLER_SUM(:)      = 0.0
ZLEV_SUM(:)      = 0.0
ZLEI_SUM(:)      = 0.0
ZLES3L_SUM(:)    = 0.0
ZLEL3L_SUM(:)    = 0.0
ZEVAP3L_SUM(:)   = 0.0
!
ZHU_AGG_SUM(:)   = 0.0
ZAC_AGG_SUM(:)   = 0.0
!
! momentum/turb:
!
ZUSTAR2_SUM(:)     = 0.0
ZUSTAR2SNOW_SUM(:) = 0.
ZCDSNOW_SUM(:)     = 0.
ZCHSNOW_SUM(:)     = 0.
ZRISNOW_SUM(:)     = 0.
!
! surface interfacial/sub-surface heat fluxes:
!
ZGRNDFLUX_SUM(:) = 0.0
ZRESTORE_SUM(:)  = 0.0
ZHPSNOW_SUM(:)   = 0.0
!
! radiative fluxes:
!
ZSWNET_V_SUM(:)  = 0.0
ZSWNET_G_SUM(:)  = 0.0
ZSWNET_N_SUM(:)  = 0.0
ZLWNET_V_SUM(:)  = 0.0
ZLWNET_G_SUM(:)  = 0.0
ZLWNET_N_SUM(:)  = 0.0
ZEMIST_SUM(:)    = 0.0
ZSWUP_SUM(:)     = 0.0
ZLWUP_SUM(:)     = 0.0
!
ZDELHEATV_SFC_SUM(:) = 0.0 
ZDELHEATG_SFC_SUM(:) = 0.0 
ZDELHEATG_SUM(:)     = 0.0 
!
IF (LHOOK) CALL DR_HOOK('ISBA_MEB:INIT_SUM_FLUXES_MEB_TSPLIT ',1,ZHOOK_HANDLE)
!
END SUBROUTINE INIT_SUM_FLUXES_MEB_TSPLIT
!
!===============================================================================
!
SUBROUTINE SUM_FLUXES_MEB_TSPLIT 
!
IMPLICIT NONE
!
!*      0.2    declarations of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ISBA_MEB:SUM_FLUXES_MEB_TSPLIT ',0,ZHOOK_HANDLE)
!
! Sum fluxes over MEB TIME SPLIT:
!
! sensible heat fluxes:

ZH_SUM(:)        = ZH_SUM(:)        + DGIP%XH(:)
ZH_C_A_SUM(:)    = ZH_C_A_SUM(:)    + DGEIP%XH_C_A(:)
ZH_N_A_SUM(:)    = ZH_N_A_SUM(:)    + ZH_N_A(:)
ZH_V_C_SUM(:)    = ZH_V_C_SUM(:)    + DGEIP%XH_V_C(:)
ZH_G_C_SUM(:)    = ZH_G_C_SUM(:)    + DGEIP%XH_G_C(:)
ZH_N_C_SUM(:)    = ZH_N_C_SUM(:)    + DGEIP%XH_N_C(:)
ZHSNOW_SUM(:)    = ZHSNOW_SUM(:)    + DGMI%XHSNOW(:)
!
! latent heat/water vapor fluxes:
!
ZLE_SUM(:)       = ZLE_SUM(:)       + IR%XLE(:,1)
ZLE_C_A_SUM(:)   = ZLE_C_A_SUM(:)   + DGEIP%XLE_C_A(:)
ZLE_V_C_SUM(:)   = ZLE_V_C_SUM(:)   + DGEIP%XLE_V_C(:) 
ZLE_G_C_SUM(:)   = ZLE_G_C_SUM(:)   + DGEIP%XLE_G_C(:) 
ZLE_N_C_SUM(:)   = ZLE_N_C_SUM(:)   + DGEIP%XLE_N_C(:) 
ZLETR_V_C_SUM(:) = ZLETR_V_C_SUM(:) + DGEIP%XLETRCV(:) 
ZLEG_SUM(:)      = ZLEG_SUM(:)      + DGEIP%XLEG(:) 
ZLEGI_SUM(:)     = ZLEGI_SUM(:)     + DGEIP%XLEGI(:) 
ZLESFC_SUM(:)    = ZLESFC_SUM(:)    + ZLESFC(:) 
ZLESFCI_SUM(:)   = ZLESFCI_SUM(:)   + ZLESFCI(:) 
ZLER_V_C_SUM(:)  = ZLER_V_C_SUM(:)  + DGEIP%XLERCV(:) 
ZLE_FLOOD_SUM(:) = ZLE_FLOOD_SUM(:) + DGEIP%XLE_FLOOD(:)
ZLEI_FLOOD_SUM(:)= ZLEI_FLOOD_SUM(:)+ DGEIP%XLEI_FLOOD(:) 
ZLES_V_C_SUM(:)  = ZLES_V_C_SUM(:)  + DGEIP%XLESC(:)
ZLETR_SUM(:)     = ZLETR_SUM(:)     + DGEIP%XLETR(:) 
ZLER_SUM(:)      = ZLER_SUM(:)      + DGEIP%XLER(:)
ZLEV_SUM(:)      = ZLEV_SUM(:)      + DGEIP%XLEV(:) 
ZLEI_SUM(:)      = ZLEI_SUM(:)      + DGIP%XLEI(:)
ZLES3L_SUM(:)    = ZLES3L_SUM(:)    + ZLES3L(:) 
ZLEL3L_SUM(:)    = ZLEL3L_SUM(:)    + ZLEL3L(:) 
ZEVAP3L_SUM(:)   = ZEVAP3L_SUM(:)   + ZEVAP3L(:) 
!
ZHU_AGG_SUM(:)   = ZHU_AGG_SUM(:)   + PHU_AGG(:)
ZAC_AGG_SUM(:)   = ZAC_AGG_SUM(:)   + PAC_AGG(:)
!
! momentum/turb:
!
ZUSTAR2_SUM(:)     = ZUSTAR2_SUM(:)     + ZUSTAR2_IC(:)
ZUSTAR2SNOW_SUM(:) = ZUSTAR2SNOW_SUM(:) + ZUSTAR2SNOW(:)
ZCDSNOW_SUM(:)     = ZCDSNOW_SUM(:)     + DGMI%XCDSNOW(:)
ZCHSNOW_SUM(:)     = ZCHSNOW_SUM(:)     + DGMI%XCHSNOW(:)
ZRISNOW_SUM(:)     = ZRISNOW_SUM(:)     + PRISNOW(:)
!
! surface interfacial/sub-surface heat fluxes:
!
ZGRNDFLUX_SUM(:) = ZGRNDFLUX_SUM(:) + PGRNDFLUX(:) 
ZRESTORE_SUM(:)  = ZRESTORE_SUM(:)  + DGEIP%XRESTORE(:) 
ZHPSNOW_SUM(:)   = ZHPSNOW_SUM(:)   + DGMI%XHPSNOW(:)
!
! radiative fluxes:
!
ZSWNET_V_SUM(:)  = ZSWNET_V_SUM(:)  +   DGEIP%XSWNET_V(:)
ZSWNET_G_SUM(:)  = ZSWNET_G_SUM(:)  +   DGEIP%XSWNET_G(:) 
ZSWNET_N_SUM(:)  = ZSWNET_N_SUM(:)  +   DGEIP%XSWNET_N(:) 
ZLWNET_V_SUM(:)  = ZLWNET_V_SUM(:)  +   DGEIP%XLWNET_V(:) 
ZLWNET_G_SUM(:)  = ZLWNET_G_SUM(:)  +   DGEIP%XLWNET_G(:)
ZLWNET_N_SUM(:)  = ZLWNET_N_SUM(:)  +   DGEIP%XLWNET_N(:) 
ZEMIST_SUM(:)    = ZEMIST_SUM(:)    +   PEMIST(:) 
ZSWUP_SUM(:)     = ZSWUP_SUM(:)     +   ZSWUP(:)
ZLWUP_SUM(:)     = ZLWUP_SUM(:)     +   ZLWUP(:)
!
ZDELHEATV_SFC_SUM(:) = ZDELHEATV_SFC_SUM(:) +   PDELHEATV_SFC(:) 
ZDELHEATG_SFC_SUM(:) = ZDELHEATG_SFC_SUM(:) +   PDELHEATG_SFC(:) 
ZDELHEATG_SUM(:)     = ZDELHEATG_SUM(:)     +   PDELHEATG(:) 
!
IF (LHOOK) CALL DR_HOOK('ISBA_MEB:SUM_FLUXES_MEB_TSPLIT ',1,ZHOOK_HANDLE)
!
END SUBROUTINE SUM_FLUXES_MEB_TSPLIT
!
!===============================================================================
!
SUBROUTINE AVG_FLUXES_MEB_TSPLIT 
!
USE MODD_CSTS, ONLY : XSTEFAN
!
IMPLICIT NONE
!
!*      0.2    declarations of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ISBA_MEB:AVG_FLUXES_MEB_TSPLIT ',0,ZHOOK_HANDLE)
!
! Average fluxes over MEB TIME SPLIT:
!
! sensible heat fluxes:
!
DGIP%XH(:)      = ZH_SUM(:)     /JTSPLIT_EB
DGEIP%XH_C_A(:) = ZH_C_A_SUM(:) /JTSPLIT_EB
ZH_N_A(:)       = ZH_N_A_SUM(:) /JTSPLIT_EB
DGEIP%XH_V_C(:) = ZH_V_C_SUM(:) /JTSPLIT_EB
DGEIP%XH_G_C(:) = ZH_G_C_SUM(:) /JTSPLIT_EB
DGEIP%XH_N_C(:) = ZH_N_C_SUM(:) /JTSPLIT_EB
DGMI%XHSNOW(:)  = ZHSNOW_SUM(:) /JTSPLIT_EB
!
! latent heat/water vapor fluxes:
!
IR%XLE(:,1)         = ZLE_SUM(:)       /JTSPLIT_EB
DGEIP%XLE_C_A(:)    = ZLE_C_A_SUM(:)   /JTSPLIT_EB
DGEIP%XLE_V_C(:)    = ZLE_V_C_SUM(:)   /JTSPLIT_EB
DGEIP%XLE_G_C(:)    = ZLE_G_C_SUM(:)   /JTSPLIT_EB
DGEIP%XLE_N_C(:)    = ZLE_N_C_SUM(:)   /JTSPLIT_EB
DGEIP%XLETRCV(:)    = ZLETR_V_C_SUM(:) /JTSPLIT_EB
DGEIP%XLEG(:)       = ZLEG_SUM(:)      /JTSPLIT_EB
DGEIP%XLEGI(:)      = ZLEGI_SUM(:)     /JTSPLIT_EB
ZLESFC(:)           = ZLESFC_SUM(:)    /JTSPLIT_EB
ZLESFCI(:)          = ZLESFCI_SUM(:)   /JTSPLIT_EB
DGEIP%XLERCV(:)     = ZLER_V_C_SUM(:)  /JTSPLIT_EB
DGEIP%XLE_FLOOD(:)  = ZLE_FLOOD_SUM(:) /JTSPLIT_EB
DGEIP%XLEI_FLOOD(:) = ZLEI_FLOOD_SUM(:)/JTSPLIT_EB
DGEIP%XLESC(:)      = ZLES_V_C_SUM(:)  /JTSPLIT_EB
DGEIP%XLETR(:)      = ZLETR_SUM(:)     /JTSPLIT_EB
DGEIP%XLER(:)       = ZLER_SUM(:)      /JTSPLIT_EB
DGEIP%XLEV(:)       = ZLEV_SUM(:)      /JTSPLIT_EB
DGIP%XLEI(:)        = ZLEI_SUM(:)      /JTSPLIT_EB
ZLES3L(:)           = ZLES3L_SUM(:)    /JTSPLIT_EB
ZLEL3L(:)           = ZLEL3L_SUM(:)    /JTSPLIT_EB
ZEVAP3L(:)          = ZEVAP3L_SUM(:)   /JTSPLIT_EB
!
PHU_AGG(:)   = ZHU_AGG_SUM(:)   /JTSPLIT_EB
PAC_AGG(:)   = ZAC_AGG_SUM(:)   /JTSPLIT_EB
!
! momentum/turb:
!
PUSTAR(:)          = SQRT( ZUSTAR2_SUM(:)    /JTSPLIT_EB )
DGMI%XUSTARSNOW(:) = SQRT( ZUSTAR2SNOW_SUM(:)/JTSPLIT_EB )
DGMI%XCDSNOW(:)    = ZCDSNOW_SUM(:)          /JTSPLIT_EB 
DGMI%XCHSNOW(:)    = ZCHSNOW_SUM(:)          /JTSPLIT_EB 
PRISNOW(:)         = ZRISNOW_SUM(:)          /JTSPLIT_EB 
!
! surface interfacial/sub-surface heat fluxes:
!
PGRNDFLUX(:)      = ZGRNDFLUX_SUM(:) /JTSPLIT_EB
DGEIP%XRESTORE(:) = ZRESTORE_SUM(:)  /JTSPLIT_EB
DGMI%XHPSNOW(:)   = ZHPSNOW_SUM(:)   /JTSPLIT_EB
!
! radiative fluxes:
!
DGEIP%XSWNET_V(:)  = ZSWNET_V_SUM(:)  /JTSPLIT_EB
DGEIP%XSWNET_G(:)  = ZSWNET_G_SUM(:)  /JTSPLIT_EB
DGEIP%XSWNET_N(:)  = ZSWNET_N_SUM(:)  /JTSPLIT_EB
DGEIP%XLWNET_V(:)  = ZLWNET_V_SUM(:)  /JTSPLIT_EB
DGEIP%XLWNET_G(:)  = ZLWNET_G_SUM(:)  /JTSPLIT_EB
DGEIP%XLWNET_N(:)  = ZLWNET_N_SUM(:)  /JTSPLIT_EB
PEMIST(:)    = ZEMIST_SUM(:)    /JTSPLIT_EB
ZSWUP(:)     = ZSWUP_SUM(:)     /JTSPLIT_EB
ZLWUP(:)     = ZLWUP_SUM(:)     /JTSPLIT_EB
!
PDELHEATV_SFC(:) = ZDELHEATV_SFC_SUM(:) /JTSPLIT_EB 
PDELHEATG_SFC(:) = ZDELHEATG_SFC_SUM(:) /JTSPLIT_EB 
PDELHEATG(:)     = ZDELHEATG_SUM(:)     /JTSPLIT_EB 
!
! Additional diagnostics depending on AVG quantities:
!
DGIP%XTSRAD(:)  = ((ZLWUP(:) - PLW_RAD(:)*(1.0-PEMIST(:)))/(XSTEFAN*PEMIST(:)))**0.25
!
ZRNET_V(:)      = DGEIP%XSWNET_V(:) + DGEIP%XLWNET_V(:)
!
ZRNET_G(:)      = DGEIP%XSWNET_G(:) + DGEIP%XLWNET_G(:)
!
DGMI%XRNSNOW(:) = DGEIP%XSWNET_N(:) + DGEIP%XLWNET_N(:)
!
DGIP%XRN(:)     = ZRNET_V(:) + ZRNET_G(:) + DGMI%XRNSNOW(:) 
!
IF (LHOOK) CALL DR_HOOK('ISBA_MEB:AVG_FLUXES_MEB_TSPLIT ',1,ZHOOK_HANDLE)
!
END SUBROUTINE AVG_FLUXES_MEB_TSPLIT 
!
!===============================================================================
SUBROUTINE SNOWALB_SPECTRAL_BANDS_MEB(PVEGTYPE,IR,PSNOWRHO,PSNOWAGE,PPS, &
                                      PSNOWDZ,PZENITH,PTAU_N)
!
! Split Total snow albedo into N-spectral bands
! NOTE currently MEB only uses 2 bands of the 3 possible.
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_DATA_COVER_PAR, ONLY : NVT_SNOW
USE MODD_MEB_PAR,        ONLY : XSW_WGHT_VIS, XSW_WGHT_NIR
USE MODD_SNOW_PAR,       ONLY : NSPEC_BAND_SNOW
USE MODD_SNOW_METAMO,    ONLY : XSNOWDZMIN
!
USE MODE_SNOW3L,         ONLY : SNOW3LALB, SNOW3LDOPT 
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PVEGTYPE      ! fraction of each vegetation (-)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWRHO      ! Snow layer average density (kg/m3)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWDZ       ! Snow layer thickness (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PZENITH       ! Zenith angle (rad)
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWAGE      ! Snow grain age
REAL, DIMENSION(:),   INTENT(IN)    :: PPS           ! Pressure [Pa]
REAL, DIMENSION(:,:), INTENT(OUT)   :: PTAU_N        ! SW absorption (coef) in uppermost snow layer (-)
!
!*      0.2    declarations of local variables
!
INTEGER                             :: JJ, JI, INJ, INLVLS
REAL, DIMENSION(SIZE(PPS))          :: ZWORK, ZWORKA, ZAGE
REAL, DIMENSION(SIZE(PPS))          :: ZPROJLAT, ZDSGRAIN, ZBETA1, ZBETA2, ZBETA3, &
                                       ZOPTICALPATH1, ZOPTICALPATH2, ZOPTICALPATH3
REAL, DIMENSION(SIZE(PPS))          :: ZPERMSNOWFRAC
REAL, DIMENSION(SIZE(PSNOWDZ,1),SIZE(PSNOWDZ,2)) :: ZSNOWDZ
REAL, DIMENSION(SIZE(PPS),NSPEC_BAND_SNOW)       :: ZSPECTRALALBEDO
!                                      ZSPECTRALALBEDO = spectral albedo (3 bands in algo: 
!                                                        MEB currently uses 2)
!                                                        1=VIS, 2=NIR, 3=UV
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ISBA_MEB:SNOWALB_SPECTRAL_BANDS_MEB',0,ZHOOK_HANDLE)
!
INJ    = SIZE(PSNOWDZ,1)
INLVLS = SIZE(PSNOWDZ,2)
!
! 1) Spectral albedo
! ------------------
!
ZWORK(:)         = 0.0
ZWORKA(:)        = IR%TSNOW%ALB(:,1)
ZPERMSNOWFRAC(:) = PVEGTYPE(:,NVT_SNOW)
!
CALL SNOW3LALB(ZWORKA,ZSPECTRALALBEDO,PSNOWRHO(:,1),PSNOWAGE(:,1),ZPERMSNOWFRAC,PPS)
!
! Since we only consider VIS and NIR bands for soil and veg in MEB currently:
! (also note, here PSNOWALB doesn't evolve...we just diagnose spectral components).
!
WHERE(IR%TSNOW%ALB(:,1)/=XUNDEF)
!
   IR%TSNOW%ALBVIS(:,1) = ZSPECTRALALBEDO(:,1)
!
! We diagnose NIR albedo such that total albedo is conserved
! (using just 2 spectral bands in MEB)
!
   IR%TSNOW%ALBNIR(:,1) = (IR%TSNOW%ALB(:,1) - XSW_WGHT_VIS*IR%TSNOW%ALBVIS(:,1))/XSW_WGHT_NIR
!
! currently NOT used by MEB
!
   IR%TSNOW%ALBFIR(:,1) = XUNDEF                                     
!
! For the surface layer absorbtion computation:
!
   ZSPECTRALALBEDO(:,1) = IR%TSNOW%ALBVIS(:,1)
   ZSPECTRALALBEDO(:,2) = IR%TSNOW%ALBNIR(:,1)
   ZSPECTRALALBEDO(:,3) = IR%TSNOW%ALBFIR(:,1)
!
ELSEWHERE
!
   IR%TSNOW%ALBVIS(:,1) = XUNDEF
   IR%TSNOW%ALBNIR(:,1) = XUNDEF
   IR%TSNOW%ALBFIR(:,1) = XUNDEF
!
END WHERE
!
! Snow optical grain diameter (no age dependency over polar regions):
!
ZAGE(:) = (1.0-ZPERMSNOWFRAC(:))*PSNOWAGE(:,1)
!
ZDSGRAIN(:) = SNOW3LDOPT(PSNOWRHO(:,1),ZAGE)
!
! 2) SW absorption in uppermost snow layer 
! ----------------------------------------
! For now, consider just 2 bands with MEB, so renormalize:

ZSPECTRALALBEDO(:,1) = ZSPECTRALALBEDO(:,1)
ZSPECTRALALBEDO(:,2) = (IR%TSNOW%ALB(:,1) - XSW_WGHT_VIS*ZSPECTRALALBEDO(:,1))/XSW_WGHT_NIR
!
! Adjust thickness to be as in snow computations:
!
DO JJ=1,INLVLS
   DO JI=1,INJ
      ZSNOWDZ(JI,JJ) = PSNOWDZ(JI,JJ)/MAX(1.E-4,IR%XPSN(JI,1))
   ENDDO
ENDDO
!
CALL SNOW3LRADTRANS(XSNOWDZMIN, ZSPECTRALALBEDO, ZSNOWDZ, PSNOWRHO, &
                           ZPERMSNOWFRAC, PZENITH,  PSNOWAGE, PTAU_N)
!
IF (LHOOK) CALL DR_HOOK('ISBA_MEB:SNOWALB_SPECTRAL_BANDS_MEB',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWALB_SPECTRAL_BANDS_MEB
!===============================================================================
      SUBROUTINE SNOW3LRADTRANS(PSNOWDZMIN, PSPECTRALALBEDO, PSNOWDZ, PSNOWRHO, &
                                PPERMSNOWFRAC, PZENITH,  PSNOWAGE, PRADTRANS)
!
!!    PURPOSE
!!    -------
!     Calculate the transmission of shortwave (solar) radiation
!     through the snowpack (using a form of Beer's Law: exponential
!     decay of radiation with increasing snow depth).
!
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_SNOW_PAR, ONLY : XVSPEC1,XVSPEC2,XVSPEC3,XVBETA1,XVBETA2, &
                          XVBETA4,XVBETA3,XVBETA5, XMINCOSZEN
USE MODD_MEB_PAR,  ONLY : XSW_WGHT_VIS, XSW_WGHT_NIR
!
USE MODE_SNOW3L,   ONLY : SNOW3LDOPT
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,                 INTENT(IN)    :: PSNOWDZMIN
!
REAL, DIMENSION(:),   INTENT(IN)    :: PPERMSNOWFRAC
REAL, DIMENSION(:),   INTENT(IN)    :: PZENITH
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWRHO, PSNOWDZ, PSNOWAGE
REAL, DIMENSION(:,:), INTENT(IN)    :: PSPECTRALALBEDO
!
REAL, DIMENSION(:,:), INTENT(OUT)   :: PRADTRANS
!
!
!*      0.2    declarations of local variables
!
INTEGER                              :: JJ, JI
!
INTEGER                              :: INJ
INTEGER                              :: INLVLS
!
REAL, DIMENSION(SIZE(PSNOWRHO,1))    :: ZRADTOT, ZPROJLAT, ZCOSZEN
REAL, DIMENSION(SIZE(PSNOWRHO,1))    :: ZOPTICALPATH1, ZOPTICALPATH2, ZOPTICALPATH3
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZDSGRAIN, ZCOEF, ZSNOWDZ, ZAGE
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZBETA1, ZBETA2, ZBETA3, ZWORK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
! 0. Initialization:
! ------------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3LRADTRANS',0,ZHOOK_HANDLE)
!
INJ    = SIZE(PSNOWDZ(:,:),1)
INLVLS = SIZE(PSNOWDZ(:,:),2)
!
!
! 1. Vanishingly thin snowpack check:
! -----------------------------------
!    For vanishingly thin snowpacks, much of the radiation
!    can pass through snowpack into underlying soil, making
!    a large (albeit temporary) thermal gradient: by imposing
!    a minimum thickness, this increases the radiation absorbtion
!    for vanishingly thin snowpacks.
!
ZSNOWDZ(:,:) = MAX(PSNOWDZMIN, PSNOWDZ(:,:))
!
!
! 2. Extinction of net shortwave radiation
! ----------------------------------------
! Fn of snow depth and density (Loth and Graf 1993:
! SNOWCVEXT => from Bohren and Barkstrom 1974
! SNOWAGRAIN and SNOWBGRAIN=> from Jordan 1976)
!
! Coefficient for taking into account the increase of path length of rays
! in snow due to zenithal angle
!
ZCOSZEN(:)=MAX(XMINCOSZEN,COS(PZENITH(:)))
!
! This formulation is incorrect but it compensate partly the fact that 
! the albedo formulation does not account for zenithal angle.
! Only for polar or glacier regions
!
ZPROJLAT(:)=(1.0-PPERMSNOWFRAC(:))+PPERMSNOWFRAC(:)/ZCOSZEN(:)
!
! Snow optical grain diameter (no age dependency over polar regions):
!
ZAGE(:,:) = 0.
DO JJ=1,INLVLS
   DO JI=1,INJ
      IF(PSNOWAGE(JI,JJ)/=XUNDEF)THEN
         ZAGE(JI,JJ) = (1.0-PPERMSNOWFRAC(JI))*PSNOWAGE(JI,JJ)
      ENDIF
   ENDDO
ENDDO
!
ZDSGRAIN(:,:) = SNOW3LDOPT(PSNOWRHO,ZAGE)
!
! Extinction coefficient from Brun et al. (1989):
!
ZWORK(:,:)=SQRT(ZDSGRAIN(:,:))
!
ZBETA1(:,:)=MAX(XVBETA1*PSNOWRHO(:,:)/ZWORK(:,:),XVBETA2)
ZBETA2(:,:)=MAX(XVBETA3*PSNOWRHO(:,:)/ZWORK(:,:),XVBETA4)
ZBETA3(:,:)=XVBETA5
!
ZOPTICALPATH1(:) = 0.0
ZOPTICALPATH2(:) = 0.0
ZOPTICALPATH3(:) = 0.0
!
DO JJ=1,INLVLS
   DO JI=1,INJ
      !
         ZOPTICALPATH1(JI) = ZOPTICALPATH1(JI) + ZBETA1(JI,JJ)*ZSNOWDZ(JI,JJ)
         ZOPTICALPATH2(JI) = ZOPTICALPATH2(JI) + ZBETA2(JI,JJ)*ZSNOWDZ(JI,JJ)

         ZCOEF (JI,JJ) = XSW_WGHT_VIS*(1.0-PSPECTRALALBEDO(JI,1))*EXP(-ZOPTICALPATH1(JI)*ZPROJLAT(JI)) &
                       + XSW_WGHT_NIR*(1.0-PSPECTRALALBEDO(JI,2))*EXP(-ZOPTICALPATH2(JI)*ZPROJLAT(JI)) 

   ENDDO
ENDDO
!
! 3. Radiation trans at base of each layer
! ----------------------------------
! NOTE, at level=0, rad = Swnet*(1-alb)  so radcoef(0)=1 implicitly
!
PRADTRANS(:,:)  = ZCOEF(:,:)
!
IF (LHOOK) CALL DR_HOOK('SNOW3LRADTRANS',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOW3LRADTRANS
!===============================================================================

SUBROUTINE ALLOCATE_LOCAL_VARS_PREP_GRID_SOIL
!
IMPLICIT NONE
!
!*      0.2    declarations of local variables
!
INTEGER         :: INLL
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ISBA_MEB:ALLOCATE_LOCAL_VARS_PREP_GRID_SOIL ',0,ZHOOK_HANDLE)

INLL = INL
IF(IO%LMEB_LITTER)INLL = INL + 1

ALLOCATE ( ZTGL        (INJ, INLL ))
ALLOCATE ( ZSOILHCAPZ (INJ, INLL ))
ALLOCATE ( ZSOILCONDZ (INJ, INLL ))
ALLOCATE ( ZD_G       (INJ, INLL ))
ALLOCATE ( ZDZG       (INJ, INLL ))
ALLOCATE ( ZWFC       (INJ, INLL ))
ALLOCATE ( ZWSAT       (INJ, INLL ))

IF (LHOOK) CALL DR_HOOK('ISBA_MEB:ALLOCATE_LOCAL_VARS_PREP_GRID_SOIL ',1,ZHOOK_HANDLE)

END SUBROUTINE ALLOCATE_LOCAL_VARS_PREP_GRID_SOIL
!===============================================================================
SUBROUTINE DEALLOCATE_LOCAL_VARS_PREP_GRID_SOIL
!
IMPLICIT NONE
!
!*      0.2    declarations of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ISBA_MEB:DEALLOCATE_LOCAL_VARS_PREP_GRID_SOIL ',0,ZHOOK_HANDLE)

DEALLOCATE ( ZTGL        )
DEALLOCATE ( ZSOILHCAPZ )
DEALLOCATE ( ZSOILCONDZ )
DEALLOCATE ( ZD_G       )
DEALLOCATE ( ZDZG       )
DEALLOCATE ( ZWSAT       )
DEALLOCATE ( ZWFC       )

IF (LHOOK) CALL DR_HOOK('ISBA_MEB:DEALLOCATE_LOCAL_VARS_PREP_GRID_SOIL ',1,ZHOOK_HANDLE)

END SUBROUTINE DEALLOCATE_LOCAL_VARS_PREP_GRID_SOIL
!===============================================================================
SUBROUTINE RESHIFT_MEB_SOIL(OMEB_LITTER,PTGL,PLESFC,PLESFCI,IR,DGEIP)
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
LOGICAL,                INTENT(IN)    :: OMEB_LITTER
REAL,   DIMENSION(:,:), INTENT(IN)    :: PTGL
REAL,   DIMENSION(:),   INTENT(IN)    :: PLESFC
REAL,   DIMENSION(:),   INTENT(IN)    :: PLESFCI
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEIP
!
!*      0.2    declarations of local variables
!
INTEGER                               :: JJ, JL
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------

INJ  = SIZE(IR%XTG(:,:,1),1)
INL  = SIZE(IR%XTG(:,:,1),2)

IF (LHOOK) CALL DR_HOOK('ISBA_MEB:FINISH_MEB_SOIL ',0,ZHOOK_HANDLE)

IF (OMEB_LITTER)THEN

   IR%XTL(:,1)           = PTGL(:,1)

   DO JL=1,INL
      DO JJ=1,INJ
         IR%XTG(JJ,JL,1) = PTGL(JJ,JL+1) 
      ENDDO
   ENDDO

   DGEIP%XLEG(:)          = 0.0
   DGEIP%XLEGI(:)         = 0.0
   DGEIP%XLELITTER(:)     = PLESFC(:)
   DGEIP%XLELITTERI(:)    = PLESFCI(:)
ELSE

   IR%XTG(:,:,1)          = PTGL(:,:) 

   DGEIP%XLEG(:)          = PLESFC(:)
   DGEIP%XLEGI(:)         = PLESFCI(:)
   DGEIP%XLELITTER(:)     = 0.
   DGEIP%XLELITTERI(:)    = 0.

ENDIF


IF (LHOOK) CALL DR_HOOK('ISBA_MEB:FINISH_MEB_SOIL ',1,ZHOOK_HANDLE)

END SUBROUTINE RESHIFT_MEB_SOIL
!===============================================================================
SUBROUTINE PREP_MEB_SOIL(OMEB_LITTER,PSOILHCAPZ,PSOILCONDZ,PD_G,IP,IR,PGNDLITTER,PD_GL,&
                         PDZGL,PTGL,PSOILHCAPL,PSOILCONDL,PWSATL,PWFCL,PWSFC,PWISFC,PCTSFC,PCT)
!
USE MODD_CSTS,       ONLY : XRHOLW
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
LOGICAL,                INTENT(IN)    :: OMEB_LITTER
REAL,   DIMENSION(:,:), INTENT(IN)    :: PSOILHCAPZ
REAL,   DIMENSION(:,:), INTENT(IN)    :: PSOILCONDZ
REAL,   DIMENSION(:,:), INTENT(IN)    :: PD_G
!
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
!
REAL,   DIMENSION(:),   INTENT(IN)    :: PCT
REAL,   DIMENSION(:),   INTENT(IN)    :: PGNDLITTER
REAL,   DIMENSION(:,:), INTENT(OUT)   :: PD_GL
REAL,   DIMENSION(:,:), INTENT(OUT)   :: PDZGL
REAL,   DIMENSION(:,:), INTENT(OUT)   :: PTGL
REAL,   DIMENSION(:,:), INTENT(OUT)   :: PSOILHCAPL
REAL,   DIMENSION(:,:), INTENT(OUT)   :: PSOILCONDL
REAL,   DIMENSION(:,:), INTENT(OUT)   :: PWSATL
REAL,   DIMENSION(:,:), INTENT(OUT)   :: PWFCL
REAL,   DIMENSION(:),   INTENT(OUT)   :: PWSFC
REAL,   DIMENSION(:),   INTENT(OUT)   :: PWISFC
REAL,   DIMENSION(:),   INTENT(OUT)   :: PCTSFC
!
!*      0.2    declarations of local variables
!
INTEGER                               :: INJ, INL, JJ, JL
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*      0.3    declarations of local parameters
!
REAL, PARAMETER                       :: Z1 = 1900.0     !massic organic matter heat capacity (J/kg/K)
REAL, PARAMETER                       :: Z2 = 45.0       !litter bulk density (kg/m3)
REAL, PARAMETER                       :: Z3 = 4180.0     !massic water heat capacity (J/kg/K)
REAL, PARAMETER                       :: Z4 = 0.1        !coeff for litter conductivity (K/m)
REAL, PARAMETER                       :: Z5 = 0.03        !coeff for litter conductivity
REAL, PARAMETER                       :: Z6 = 0.95       !litter porosity       (m3/m3)
REAL, PARAMETER                       :: Z7 = 0.12       !litter field capacity (m3/m3)

!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ISBA_MEB:PREP_MEB_SOIL',0,ZHOOK_HANDLE)
!
INJ  = SIZE(IR%XTG,1)
INL  = SIZE(IR%XTG,2)
!
ZWORK(:) = 0.0
IF(OMEB_LITTER)THEN
   PTGL(:,1)                  = IR%XTL(:,1)
   ZWORK(:)                   = IR%XWRL(:,1)/PGNDLITTER(:)
   PSOILHCAPL(:,1)            = Z1*Z2 + Z3*1000/XRHOLW*ZWORK(:)
   PSOILCONDL(:,1)            = Z4 +  (Z5/XRHOLW)     *ZWORK(:)
   PWSATL(:,1)                = Z6
   PWFCL(:,1)                 = Z7
   PD_GL(:,1)                 = PGNDLITTER(:)
   PDZGL(:,1)                 = PGNDLITTER(:)
   PCTSFC(:)                  = 1. / (PSOILHCAPL(:,1) * PGNDLITTER(:))

   DO JL=1,INL
      DO JJ=1,INJ
         PTGL(JJ,JL+1)        = IR%XTG(JJ,JL,1) 
         PSOILHCAPL(JJ,JL+1)  = PSOILHCAPZ(JJ,JL)
         PSOILCONDL(JJ,JL+1)  = PSOILCONDZ(JJ,JL)
         PWSATL(JJ,JL+1)      = IP%XWSAT(JJ,JL)
         PWFCL(JJ,JL+1)       = IP%XWFC(JJ,JL)
         PD_GL(JJ,JL+1)       = PGNDLITTER(JJ) + PD_G(JJ,JL)
         PDZGL(JJ,JL+1)       = IP%XDZG(JJ,JL,1)
      ENDDO
   ENDDO
   PWSFC(:)                   = IR%XWRL(:,1) /(XRHOLW*PGNDLITTER(:)) ! (m3/m3)
   PWISFC(:)                  = IR%XWRLI(:,1)/(XRHOLW*PGNDLITTER(:)) ! (m3/m3)

ELSE
   PTGL(:,:)                  = IR%XTG(:,:,1)
   PSOILHCAPL(:,:)            = PSOILHCAPZ(:,:)
   PSOILCONDL(:,:)            = PSOILCONDZ(:,:)
   PWSATL(:,:)                = IP%XWSAT(:,:)
   PWFCL(:,:)                 = IP%XWFC(:,:)
   PD_GL(:,:)                 = PD_G(:,:)
   PDZGL(:,:)                 = IP%XDZG(:,:,1)
   PCTSFC(:)                  = PCT(:)
   PWSFC(:)                   = IR%XWG(:,1,1)
   PWISFC(:)                  = IR%XWGI(:,1,1)
ENDIF
IF (LHOOK) CALL DR_HOOK('ISBA_MEB:PREP_MEB_SOIL',1,ZHOOK_HANDLE)

END SUBROUTINE PREP_MEB_SOIL
!===============================================================================
SUBROUTINE ICE_LITTER(PTSTEP, PLELITTERI, PSOILHCAPZ, IR, &
                      KWG_LAYER, PDZG, PGNDLITTER, PPHASEL, PCTSFC   )
!
USE MODD_CSTS,     ONLY : XLMTT, XTT, XCI, XRHOLI, XRHOLW, XLSTT
USE MODD_ISBA_PAR, ONLY : XWGMIN
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
!
REAL, INTENT(IN)                    :: PTSTEP  
!                                      PTSTEP     = Model time step (s)
!
REAL, DIMENSION(:), INTENT(IN)      :: PLELITTERI
!                                      PLELITTERI = ice sublimation (m s-1)
REAL, DIMENSION(:), INTENT(IN)      :: PCTSFC
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSOILHCAPZ
!                                      PSOILHCAPZ = soil heat capacity [J/(m3 K)]
REAL, DIMENSION(:,:), INTENT(IN)    :: PDZG
!                                      PDZG       = Layer thickness (DIF option)
REAL, DIMENSION(:), INTENT(IN)      :: PGNDLITTER
!                                      PGNDLITTER        = Litter thickness (m)
!
INTEGER, DIMENSION(:), INTENT(IN)   :: KWG_LAYER  
!                                      KWG_LAYER = Number of soil moisture layers (DIF option)
!
REAL, DIMENSION(:), INTENT(OUT)     :: PPHASEL
!                                      PPHASEL = Phase changement in litter (W/m2)
!
!*      0.2    declarations of local variables
!
INTEGER                             :: JL     ! loop control
!
INTEGER                             :: INL    ! Number of explicit soil layers
!
REAL, DIMENSION(SIZE(IR%XTG(:,:,1),1))             :: ZEXCESS, ZK,ZTAUICE,ZHCAPL,ZELITTERI, &
                                            ZDELTAT,ZPHASE,ZPHASEM,ZPHASEF,ZPHASEX,&
                                            ZWRL,ZWRLI,Z0,ZWRLSAT,ZPHASEC
!
REAL                                     :: ZPSIMAX, ZPSI,  &
                                            ZTGM,ZWGIM,   &
                                            ZEFFIC, ZWORK, &
                                            ZAPPHEATCAP
!                                            !
!
!*      0.3    declaration of local parameters
!
REAL, PARAMETER                     :: ZERTOL = 1.E-6 ! error tolerance
!
!-------------------------------------------------------------------------------
!
! Initialization:
! ---------------
!
!
INL = MAXVAL(KWG_LAYER(:))
!
ZEXCESS(:)  =0.0
ZPHASEC(:)  =0.0

!
ZTAUICE = 3300.
ZWRLSAT(:) = 0.95
ZHCAPL(:)  = 1/(PCTSFC*PGNDLITTER)

!-------------------------------------------------------------------------------
!
! 1. Surface layer vegetation insulation coefficient (-)
!    ---------------------------------------------------
!
      ZK(:) = 1.0 
!
! 1.1 Convert to m3/m3
!    -----------------
!
ZWRL(:)= IR%XWRL(:,1) /(XRHOLW*PGNDLITTER(:))    
ZWRLI(:)=IR%XWRLI(:,1)/(XRHOLW*PGNDLITTER(:)) 
!
! 2. Litter ice evolution computation:
!    --------------------------------
!  
ZDELTAT(:) = IR%XTL(:,1) - XTT
!
!     
!     *Melt* ice if energy and ice available:
ZPHASEM(:)  = (PTSTEP/ZTAUICE(:))*MIN((XCI*XRHOLI)*MAX(0.0,ZDELTAT),ZWRLI(:)*(XLMTT*XRHOLW))
!
!     *Freeze* liquid water if energy and water available and do not exceed porosity:
ZPHASEF(:)  = (PTSTEP/ZTAUICE(:))*MIN(ZK(:)*(XCI*XRHOLI)*MAX(0.0,-ZDELTAT),ZWRL(:)*(XLMTT*XRHOLW))
ZPHASEF(:)  = min(ZPHASEF(:) , (ZWRLSAT(:) - 0.1 - ZWRLI(:)) * (XLMTT*XRHOLW) ) !!!!! LOOK!!!!!!!!
!
ZPHASE(:) = ZPHASEF(:) - ZPHASEM(:)

!     Update heat content if melting or freezing
IR%XTL(:,1) = IR%XTL(:,1) + ZPHASE(:)/ZHCAPL(:)                                    
!
!     Get estimate of actual total phase change (J/m3) for equivalent litter water changes:

ZPHASEX(:) = ZPHASE(:)
!
!     Adjust ice and liquid water conents (m3/m3) accordingly :
!   
ZWRL (:) = ZWRL (:) - ZPHASEX/(XLMTT*XRHOLW)
ZWRLI(:) = ZWRLI(:) + ZPHASEX/(XLMTT*XRHOLW)
!
! 2.1 Convert to Kg/m2
!    -----------------
!
IR%XWRL(:,1) = ZWRL(:)  * PGNDLITTER(:) * XRHOLW
IR%XWRLI(:,1)= ZWRLI(:) * PGNDLITTER(:) * XRHOLW
!
! 3. Adjust litter ice content for sublimation
!    -----------------------------------------
!
!ZELITTERI  = PLELITTERI /XLSTT * PTSTEP                    
ZELITTERI  = PLELITTERI * (PTSTEP/XLSTT)                    
ZEXCESS(:) = MAX( 0.0 , ZELITTERI - IR%XWRLI(:,1) )       
IR%XWRLI  (:,1) = IR%XWRLI(:,1) - ( ZELITTERI - ZEXCESS )                    
!
! 4. Prevent some possible problems
!    ------------------------------
!
IR%XWGI (:,1,1) = IR%XWGI(:,1,1)- ZEXCESS / (XRHOLW * PDZG(:,1))             
!
ZEXCESS(:) = max( 0.0, - IR%XWGI(:,1,1) )
IR%XWGI(:,1,1)  = IR%XWGI(:,1,1) + ZEXCESS(:)                                
IR%XWG (:,1,1)  = IR%XWG (:,1,1) - ZEXCESS(:)                                
IR%XTG (:,1,1)  = IR%XTG (:,1,1) + ZEXCESS(:) * (XLMTT*XRHOLW)/PSOILHCAPZ(:,1) 
!
DO JL=1,INL-1                 
   ZEXCESS = max(0.0,-IR%XWG(:,JL,1))
   IR%XWG(:,JL+1,1) = IR%XWG(:,JL+1,1) - ZEXCESS*PDZG(:,JL)/PDZG(:,JL+1)
   IR%XWG(:,JL,1)   = IR%XWG(:,JL,1)   + ZEXCESS
ENDDO
!
! 5. Prevent from keeping track of ice in litter
!    -------------------------------------------
!
DO JJ=1,INJ
   IF (IR%XWRLI(JJ,1) < ZERTOL ) THEN 
      IR%XWRL(JJ,1)= IR%XWRL(JJ,1) + IR%XWRLI(JJ,1) 
      IR%XTL(JJ,1) = IR%XTL(JJ,1) + IR%XWRLI(JJ,1) * XLMTT / PGNDLITTER(JJ) / ZHCAPL(JJ)
      ZPHASEC(:)   = IR%XWRLI(JJ,1) * XLMTT / PGNDLITTER(JJ)
      IR%XWRLI(JJ,1) = 0.0
   ELSE
      ZPHASEC(:) = 0.0
   ENDIF
ENDDO
!
PPHASEL(:)=(ZPHASE(:) + ZPHASEC(:))/PTSTEP*PGNDLITTER
!
END SUBROUTINE ICE_LITTER

!===============================================================================

END SUBROUTINE ISBA_MEB
