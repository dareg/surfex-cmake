!     #########
      SUBROUTINE ISBA(IO, IMX, IMT, IMM, IMI, P, IP, INI, IR, G, AG, DGI, DGIP,  &
                      DGEI, DGEIP, DGMI, HSNOW_ISBA, TPTIME, PPOI, PABC, PIACAN, &
                      OMEB, PTSTEP, HIMPLICIT_WIND, PZREF, PUREF, PDIRCOSZW,     &
                      PTA, PQA, PEXNA, PRHOA, PPS, PEXNS, PRR, PSR, PZENITH,     &
                      PSCA_SW, PSW_RAD, PLW_RAD, PVMOD, PPEW_A_COEF, PPEW_B_COEF,&
                      PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,        &
                      PALBNIR_TVEG, PALBVIS_TVEG, PALBNIR_TSOIL, PALBVIS_TSOIL,  &
                      PPALPHAN, PZ0G_WITHOUT_SNOW, PZ0_MEBV, PZ0H_MEBV,          &
                      PZ0EFF_MEBV, PZ0_MEBN, PZ0H_MEBN, PZ0EFF_MEBN, PTDEEP_A,   &
                      PCSP, PFFG_NOSNOW, PFFV_NOSNOW, PEMIST, PUSTAR, PAC_AGG,   &
                      PHU_AGG, PRESP_BIOMASS_INST, PDEEP_FLUX, PIRRIG_GR     )
!     ##########################################################################
!
!
!!****  *ISBA*  
!!
!!    PURPOSE
!!    -------
!       Monitor for the calculation of the surface fluxes and of the
!     prognostic variables of the surface over natural areas
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
!!      S. Belair           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/03/95 
!!      (J.Stein)   25/10/95  add the rain flux computation at the ground
!!                            and the lbc
!!      (J.Stein)   15/11/95  include the strong slopes cases
!!      (J.Stein)   06/02/96  bug correction for the precipitation flux writing 
!!      (J.Stein)   20/05/96  set the right IGRID value for the rain rate
!!      (J.Viviand) 04/02/97  add cold and convective precipitation rate
!!      (J.Stein)   22/06/97  use the absolute pressure    
!!      (V.Masson)  09/07/97  add directional z0 computations and RESA correction     
!!      (V.Masson)  13/02/98  simplify the routine: only vegetation computation
!!                            are now made here.
!!      (A.Boone)   05/10/98  add: Boone et al. (1999) 3 soil-water Layers version
!!      (V.Masson)                 Dumenil and Todini (1992) runoff
!!                                 Calvet (1998) biomass and CO2 assimilation
!!                                 Calvet (1998) LAI evolution
!!      (A.Boone)  03/15/99   Soil ice scheme: modify CG, C1, C2, WSAT, WFC, WILT,
!!                            LEG (add soil ice sublimation); Can modify TS and T2.
!!                            New variables WGI1, WGI2
!!      (A.Boone)  18/01/00   ISBA-ES (3-layer explicit snow scheme option)
!!                            (Boone and Etchevers 2000)
!!                            New variable IR%TSNOW%HEAT(:,:,1)
!!      (V. Masson) 01/2004   wet leaves fraction computed in separate routine
!!                            all vegetation stress (ISBA, AGS, AST) routines
!!                            called at the same point
!!      (P. LeMoigne) 03/2004 computation of QSAT 
!!      (P. LeMoigne) 10/2004 halstead coefficient as diagnostic for isba
!!      (A. Bogatchev)09/2005 EBA snow option
!!      (P. LeMoigne) 02/2006 z0h and snow
!!      (B. Decharme) 05/2008 Add floodplains scheme
!!      (R. Hamdi)    01/09   Cp and L are not constants (As in ALADIN)
!!      (A.L. Gibelin)  03/2009 : Add respiration diagnostics
!!      A.L. Gibelin   06/09 : move calculations of CO2 fluxes
!!      A.L. Gibelin 07/2009 : Suppress PPST and PPSTF as outputs
!!      (A. Boone)    11/2009 Add local variable: total soil temperature change (before
!!                            phase change) for use by LWT scheme in ISBA-DIF. 
!!      (A. Boone)    03/2010 Add local variable: delta functions for LEG and LEGI
!!                            to numerically correct for when they should be
!!                            zero when hug(i) Qsat < Qa and Qsat > Qa
!!     (A. Carrer)    04/2011 : new radiative transfert (AGS)
!!      (B. Decharme) 09/2012 Bug : Save snow albedo values at beginning
!!                                  of time step for total albedo calculation
!!                            Bug : flood fraction in COTWORES
!!                            new wind implicitation
!!                            Irrigation rate diag
!!     (C. de Munck) 03/2013  Specified irrigation for ground
!!      (B. Decharme) 04/2013 Bug : Wrong radiative temperature
!!                            DIF lateral subsurface drainage
!!                            Sublimation diag flux
!!                            Qs for 3l or crocus (needed for coupling with atm)
!!                            water table / surface coupling
!!                            Routines drag, e_budget and isba_fluxes now in isba_ceb
!!      (A. Boone & P. Samuelsson) (10/2014) Added MEB v1
!!      (P. LeMoigne) 12/2014 EBA scheme update
!!      (A. Boone)    02/2015 Consider spectral band dependence of snow for IO%LTR_ML radiation option 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_FIX_t, ISBA_PARAM_TIME_t, ISBA_PARAM_MEB_t, &
                              ISBA_PARAM_IRRIG_t
USE MODD_ISBA_PGD_n, ONLY : ISBA_PGD_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t, ISBA_INIT_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
USE MODD_GRID_n, ONLY : GRID_t
USE MODD_AGRI_n, ONLY : AGRI_t
USE MODD_DIAG_n, ONLY : DIAG_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t
!
USE MODD_CO2V_PAR,   ONLY : XMC, XMCO2, XPCCO2
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE MODD_CSTS,           ONLY : XTT
USE MODD_CO2V_PAR,       ONLY : XMC, XMCO2, XPCCO2
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_DATA_COVER_PAR, ONLY : NVT_SNOW
USE MODD_MEB_PAR,        ONLY : XSW_WGHT_VIS, XSW_WGHT_NIR
!
USE MODD_TYPE_DATE_SURF, ONLY : DATE_TIME
!
USE MODI_SOIL
USE MODI_SOILDIF
USE MODI_SOILSTRESS
USE MODI_WET_LEAVES_FRAC
USE MODI_VEG
USE MODI_SNOW3L_ISBA
USE MODI_HYDRO
USE MODI_ISBA_SNOW_AGR
!
USE MODI_RADIATIVE_TRANSFERT
USE MODI_COTWORES
!
!
USE MODI_ISBA_CEB
USE MODI_ISBA_MEB
!
USE MODE_THERMOS
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
TYPE(ISBA_PGD_t), INTENT(INOUT) :: P
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
TYPE(ISBA_INIT_t), INTENT(INOUT) :: INI
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
TYPE(GRID_t), INTENT(INOUT) :: G
TYPE(AGRI_t), INTENT(INOUT) :: AG
TYPE(DIAG_t), INTENT(INOUT) :: DGI
TYPE(DIAG_t), INTENT(INOUT) :: DGIP
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEI
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEIP
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DGMI
!
 CHARACTER(LEN=*),     INTENT(IN)  :: HSNOW_ISBA ! 'DEF' = Default F-R snow scheme
!                                               !         (Douville et al. 1995)
!                                               ! '3-L' = 3-L snow scheme (option)
!                                               !         (Boone and Etchevers 2000)
!                                               ! 'DEF' = no mulch effect
TYPE(DATE_TIME), INTENT(IN)       :: TPTIME     ! current date and time
!
REAL, DIMENSION(:),    INTENT(IN) :: PPOI       ! Gaussian weights (as above)
REAL, DIMENSION(:), INTENT(INOUT) :: PABC       ! abscissa needed for integration
!                                               ! of net assimilation and stomatal
!                                               ! conductance over canopy depth
REAL, DIMENSION(:,:),   INTENT(OUT) :: PIACAN   ! PAR in the canopy at different gauss level
LOGICAL, INTENT(IN)                 :: OMEB     ! True = patch with multi-energy balance 
!                                               ! False = patch with classical ISBA 
REAL,                 INTENT(IN)  :: PTSTEP     ! timestep of the integration
CHARACTER(LEN=*),     INTENT(IN)  :: HIMPLICIT_WIND   ! wind implicitation option
!                                                     ! 'OLD' = direct
!                                                     ! 'NEW' = Taylor serie, order 1
REAL, DIMENSION(:),   INTENT(IN) :: PZREF       ! normal distance of the first
!                                               ! atmospheric level to the
!                                               ! orography
REAL, DIMENSION(:),   INTENT(IN) :: PUREF       ! reference height of the wind
!                                               ! NOTE this is different from ZZREF
!                                               ! ONLY in stand-alone/forced mode,
!                                               ! NOT when coupled to a model (MesoNH)
REAL, DIMENSION(:),   INTENT(IN) ::  PDIRCOSZW  ! Director Cosinus along z
!                                               ! directions at surface w-point
!
!* atmospheric variables
!  ---------------------
!
!            suffix 'A' stands for atmospheric variable at first model level
!            suffix 'S' stands for atmospheric variable at ground level
!
REAL, DIMENSION(:), INTENT(IN)  :: PTA        ! Temperature
REAL, DIMENSION(:), INTENT(IN)  :: PQA        ! specific humidity
REAL, DIMENSION(:), INTENT(IN)  :: PEXNA      ! Exner function
REAL, DIMENSION(:), INTENT(IN)  :: PRHOA      ! air density
!
REAL, DIMENSION(:), INTENT(IN)  :: PPS        ! Pressure
REAL, DIMENSION(:), INTENT(IN)  :: PEXNS      ! Exner function
!
REAL, DIMENSION(:), INTENT(IN)  :: PRR        ! Rain rate (in kg/m2/s)
REAL, DIMENSION(:), INTENT(IN)  :: PSR        ! Snow rate (in kg/m2/s)
!
REAL, DIMENSION(:), INTENT(IN)  :: PZENITH    ! solar zenith angle
REAL, DIMENSION(:), INTENT(IN)  :: PSW_RAD    ! solar   incoming radiation
REAL, DIMENSION(:), INTENT(IN)  :: PSCA_SW    ! solar diffuse incoming radiation
REAL, DIMENSION(:), INTENT(IN)  :: PLW_RAD    ! thermal incoming radiation
!
REAL, DIMENSION(:), INTENT(IN)  :: PVMOD      ! modulus of the wind
!                                             ! parallel to the orography
!
! implicit coupling coefficients:
!
REAL, DIMENSION(:), INTENT(IN)  :: PPEW_A_COEF, PPEW_B_COEF, &
                                   PPET_A_COEF, PPEQ_A_COEF, &
                                   PPET_B_COEF, PPEQ_B_COEF  
!                                  PPEW_A_COEF ! A-wind coefficient
!                                  PPEW_B_COEF ! B-wind coefficient
!                                  PPET_A_COEF ! A-air temperature coefficient
!                                  PPET_B_COEF ! B-air temperature coefficient
!                                  PPEQ_A_COEF ! A-air specific humidity coefficient
!                                  PPEQ_B_COEF ! B-air specific humidity coefficient
!
!* vegetation parameters
!  ---------------------
!
REAL, DIMENSION(:), INTENT(IN)  :: PALBNIR_TVEG  ! tot albedo of vegetation in NIR (needed for LM_TR)
REAL, DIMENSION(:), INTENT(IN)  :: PALBVIS_TVEG  ! tot albedo of vegetation in VIS
REAL, DIMENSION(:), INTENT(IN)  :: PALBNIR_TSOIL ! tot albedo of bare soil in NIR
REAL, DIMENSION(:), INTENT(IN)  :: PALBVIS_TSOIL ! tot albedo of bare soil in VIS
!
! For multi-energy balance
!
REAL, DIMENSION(:), INTENT(IN)    :: PPALPHAN           ! snow/canopy transition coefficient
REAL, DIMENSION(:), INTENT(IN)    :: PZ0G_WITHOUT_SNOW  ! roughness length for momentum at snow-free canopy floor
REAL, DIMENSION(:), INTENT(IN)    :: PZ0_MEBV           ! roughness length for momentum over MEB vegetation part of patch
REAL, DIMENSION(:), INTENT(IN)    :: PZ0H_MEBV          ! roughness length for heat over MEB vegetation part of path
REAL, DIMENSION(:), INTENT(IN)    :: PZ0EFF_MEBV        ! roughness length for momentum over MEB vegetation part of patch
REAL, DIMENSION(:), INTENT(IN)    :: PZ0_MEBN           ! roughness length for momentum over MEB snow part of patch
REAL, DIMENSION(:), INTENT(IN)    :: PZ0H_MEBN          ! roughness length for heat over MEB snow part of path
REAL, DIMENSION(:), INTENT(IN)    :: PZ0EFF_MEBN        ! roughness length for momentum over MEB snow part of patch
!
!* soil parameters
!  ---------------
!
REAL, DIMENSION(:), INTENT(IN)  :: PTDEEP_A       ! Deep soil temperature (prescribed)
!                                                 PTDEEP_A = Deep soil temperature
!                                                 coefficient depending on flux
!* ISBA-Ags parameters
!  -------------------
!
REAL, DIMENSION(:),    INTENT(IN) :: PCSP       ! atmospheric CO2 concentration
!                                                 [ppmm]=[kg CO2 / kg air]
!
!* ISBA-DF variables/parameters:                  
!  ------------------------------
!
REAL, DIMENSION(:), INTENT(IN)   :: PFFG_NOSNOW ! Without snow (ES)
REAL, DIMENSION(:), INTENT(IN)   :: PFFV_NOSNOW ! Without snow (ES)
!
!* diagnostic variables
!  --------------------
!
REAL, DIMENSION(:), INTENT(OUT) :: PEMIST     ! grid-area surface emissivity
!
!* surface fluxes
!  --------------
!
REAL, DIMENSION(:), INTENT(OUT) :: PUSTAR     ! friction velocity
!
! The following surface fluxes are from snow-free portion of grid
! box when the ISBA-ES option is ON. Otherwise, they are equal
! to the same variables without the _ISBA extension.
!
REAL, DIMENSION(:),  INTENT(OUT) :: PAC_AGG  ! aggregated aerodynamic conductance
                                     ! for evaporative flux calculations
REAL, DIMENSION(:),  INTENT(OUT) :: PHU_AGG  ! aggregated relative humidity
                                     ! for evaporative flux calculations
!
!* diagnostic variables for Carbon assimilation
!  --------------------------------------------
!
REAL, DIMENSION(:,:),   INTENT(OUT) :: PRESP_BIOMASS_INST  ! instantaneous biomass respiration (kgCO2/kgair m/s)
!
!* diagnostic variables for multi-energy balance (MEB)
!  ---------------------------------------------------
!
REAL, DIMENSION(:),     INTENT(OUT) :: PDEEP_FLUX ! Heat flux at bottom of ISBA (W/m2)
!
REAL   ,DIMENSION(:),INTENT(IN)    :: PIRRIG_GR ! ground irrigation rate (kg/m2/s)
!
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(IR%XWR)) :: ZCS       ! heat capacity of the snow
REAL, DIMENSION(SIZE(IR%XWR)) :: ZFROZEN1  ! ice fraction in superficial soil
REAL, DIMENSION(SIZE(IR%XWR)) :: ZDELTA    ! fraction of the foliage
!                                       ! covered with intercepted
!                                       ! water
REAL, DIMENSION(SIZE(IR%XWR)) :: ZQSAT     ! expression for the saturation 
!                                       ! specific humidity 
!
REAL, DIMENSION(SIZE(IR%XWR)) :: ZWRMAX    ! maximum canopy water interception
!
REAL, DIMENSION(SIZE(IR%XWR)) :: ZF2       ! water stress coefficient
!
REAL, DIMENSION(SIZE(IR%XWR)) :: ZF5       ! water stress coefficient (based on F2)
!                                       ! to enforce Etv=>0 as F2=>0
!
REAL, DIMENSION(SIZE(IR%XWR)) :: ZHUGI    ! humidity over frozen bare ground
!
REAL, DIMENSION(SIZE(IR%XWR)) :: ZEVAPCOR ! evaporation correction as last traces of snow
!                                      ! cover ablate
REAL, DIMENSION(SIZE(IR%XWR)) :: ZLES3L   ! sublimation from ISBA-ES(3L)
REAL, DIMENSION(SIZE(IR%XWR)) :: ZLEL3L   ! evaporation heat flux of water in the snow (W/m2)
REAL, DIMENSION(SIZE(IR%XWR)) :: ZEVAP3L  ! evaporation flux over snow from ISBA-ES (kg/m2/s)
REAL, DIMENSION(SIZE(IR%XWR)) :: ZSNOW_THRUFAL ! rate that liquid water leaves snow pack: 
!                                           ! ISBA-ES [kg/(m2 s)]
REAL, DIMENSION(SIZE(IR%XWR)) :: ZALB3L   !Snow albedo at t-dt for total albedo calculation (ES/CROCUS)
REAL, DIMENSION(SIZE(IR%XWR)) :: ZRI3L    !Snow Ridcharson number (ES/CROCUS)
REAL, DIMENSION(SIZE(IR%XWR)) :: ZQS3L    ! surface humidity (kg/kg) (ES/CROCUS)
!
REAL, DIMENSION(SIZE(IR%XWR)) :: ZVEG
!
REAL, DIMENSION(SIZE(IR%XWR),SIZE(PABC)) :: ZIACAN_SHADE, ZIACAN_SUNLIT
!                                      ! absorbed PAR of each level within the
!                                      ! canopy - Split into shaded and SUNLIT
REAL, DIMENSION(SIZE(IR%XWR),SIZE(PABC)) :: ZFRAC_SUN  ! fraction of sunlit leaves
!
! ISBA-DF:
!                                                              
REAL, DIMENSION(SIZE(IR%XWG,1),SIZE(IR%XWG,2)) :: ZSOILHCAPZ ! ISBA-DF Soil heat capacity 
!                                                      ! profile [J/(m3 K)]
REAL, DIMENSION(SIZE(IR%XWG,1),SIZE(IR%XWG,2)) :: ZSOILCONDZ ! ISBA-DF Soil conductivity  
!                                                      ! profile  [W/(m K)]
!
REAL, DIMENSION(SIZE(IR%XWG,1),SIZE(IR%XWG,2)) :: ZF2WGHT    ! water stress factor
!
REAL, DIMENSION(SIZE(IR%XWR))               :: ZGRNDFLUX  ! snow/soil-biomass interface flux (W/m2)
REAL, DIMENSION(SIZE(IR%XWR))               :: ZFLSN_COR  ! snow/soil-biomass correction flux (W/m2)
!
! MEB:
!
REAL, DIMENSION(SIZE(IR%XWR))               :: ZSUBVCOR   ! A possible snow (intercepted by the canopy) mass correction 
!                                                       (to be potentially removed from soil) when MEB activated (kg/m2/s)
!
! Misc :
!
! -----------------------------------------------------------------------------------------------------------------------------------------------------
! Budget: Add to arguments, diags

REAL, DIMENSION(SIZE(IR%XWR))                   :: ZDELHEATV_SFC  ! Change in heat storage of the explicit vegetation (MEB) layer over the current time step (W m-2)
REAL, DIMENSION(SIZE(IR%XWR))                   :: ZDELHEATG      ! change in heat storage of the entire soil column over the current time step (W m-2) 
REAL, DIMENSION(SIZE(IR%XWR))                   :: ZDELHEATG_SFC  ! change in heat storage of the surface soil layer over the current time step (W m-2)
REAL, DIMENSION(SIZE(IR%XWR))                   :: ZDELPHASEG     ! latent heating due to soil freeze-thaw in the entire soil column            (W m-2) 
REAL, DIMENSION(SIZE(IR%XWR))                   :: ZDELPHASEG_SFC ! latent heating due to soil freeze-thaw in the surface soil layer            (W m-2) 
REAL, DIMENSION(SIZE(IR%XWR))                   :: ZDELHEATN      ! change in heat storage of the entire snow column over the current time step (W m-2)
REAL, DIMENSION(SIZE(IR%XWR))                   :: ZDELHEATN_SFC  ! change in heat storage of the surface snow layer over the current time step (W m-2)
REAL, DIMENSION(SIZE(IR%XWR))                   :: ZSNOWSFCH      ! snow surface layer pseudo-heating term owing to
!                                                              !  changes in grid thickness            (W m-2)
REAL, DIMENSION(SIZE(IR%XWR))                   :: ZGSFCSNOW      ! conductive heat flux between the surface and sub-surface soil layers 
!                                                              ! for the multi-layer snow schemes..for composite snow, it is 
!                                                              ! equal to DGEIP%XRESTORE (W m-2)
!
!
! Necessary to close the energy budget between surfex and the atmosphere:
!
REAL, DIMENSION(SIZE(IR%XWR))   :: ZEMIST, ZZHV
REAL, DIMENSION(SIZE(IR%XWR))   :: ZALBT, ZEV, ZETR, ZER
!
LOGICAL, DIMENSION(SIZE(IR%XTG,1))  :: GSHADE         ! mask where evolution occurs
!
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*      1.0    Preliminaries
!              -------------
!
IF (LHOOK) CALL DR_HOOK('ISBA',0,ZHOOK_HANDLE)
!
!print*,'tg1 ',IR%XTG(1,1,1)
DGMI%XC1(:)          = XUNDEF
DGMI%XC2(:)          = XUNDEF
DGMI%XWGEQ(:)        = XUNDEF
ZCS(:)          = XUNDEF
!
ZEMIST      (:) = XUNDEF
ZALBT       (:) = XUNDEF
ZRI3L       (:) = XUNDEF
!
ZSOILHCAPZ(:,:) = XUNDEF
ZSOILCONDZ(:,:) = XUNDEF
ZF2WGHT   (:,:) = XUNDEF
!
DGMI%XRS    (:)   = 0.0
PAC_AGG     (:)   = 0.0
PHU_AGG     (:)   = 0.0
DGMI%XSNOWTEMP   (:,:) = XTT
DGEIP%XMELT       (:)   = 0.0
!
!
!
! MEB:
!
ZDELHEATV_SFC (:) = 0.0
ZDELHEATG     (:) = 0.0 
ZDELHEATG_SFC (:) = 0.0
ZDELPHASEG    (:) = 0.0 
ZDELPHASEG_SFC(:) = 0.0 
ZDELHEATN     (:) = 0.0
ZDELHEATN_SFC (:) = 0.0
ZSNOWSFCH     (:) = 0.0
ZGSFCSNOW     (:) = 0.0
!
ZSUBVCOR(:)     = 0.0
!
IF(OMEB)THEN
   ZVEG(:) = 0.0
ELSE
   ZVEG(:) = IMT%XVEG(:,1)
ENDIF
!
! Save snow albedo values at beginning of time step for total albedo calculation
!
ZALB3L(:)=IR%TSNOW%ALB(:,1)
!
!-------------------------------------------------------------------------------
!
!*      2.0    Soil parameters
!              ---------------
!
IF(IO%CISBA =='2-L' .OR. IO%CISBA == '3-L')THEN
!
   CALL SOIL (IO, HSNOW_ISBA, IP, INI, IR, DGMI, ZVEG, IMT%XCV(:,1), &
              ZCS, ZFROZEN1, PFFG_NOSNOW, PFFV_NOSNOW  )  
!
ELSE
!
   CALL SOILDIF (IO, IP, IMX, IR, DGMI, ZVEG, IMT%XCV(:,1), INI%XPIFLOOD, &
                 ZFROZEN1, PFFG_NOSNOW, PFFV_NOSNOW, ZSOILCONDZ, ZSOILHCAPZ  )
!
ENDIF
!
!print*,'tg2 ',IR%XTG(1,1,1)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      3.0    Plant stress due to soil water deficit
!              --------------------------------------
!
CALL SOILSTRESS(IO%CISBA, ZF2, IP, IMX, IR, ZF2WGHT, ZF5 )  
!
!print*,'tg3 ',IR%XTG(1,1,1)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      4.0    Explicit Canopy Vegetation Option
!              ---------------------------------
!
IF(OMEB)THEN
   CALL ISBA_MEB(IO, IMX, IMT, IMM, IMI, IP, INI, IR, DGIP, DGEIP, DGMI, G, AG, &
                 TPTIME, OMEB, GSHADE, HSNOW_ISBA, HIMPLICIT_WIND, PTSTEP, &
                 ZSOILHCAPZ, ZSOILCONDZ, ZFROZEN1, PPS, PZENITH,      &
                 PSCA_SW, PSW_RAD, PVMOD, PRR, PSR, PRHOA, PTA, PQA,  &
                 PDIRCOSZW, PEXNS, PEXNA, PPET_A_COEF, PPET_B_COEF,   &
                 PPEQ_A_COEF, PPEQ_B_COEF, PPEW_A_COEF, PPEW_B_COEF,  &
                 PZREF, PUREF, PZ0G_WITHOUT_SNOW, PZ0_MEBV, PZ0H_MEBV,&
                 PZ0EFF_MEBV, PZ0_MEBN, PZ0H_MEBN, PZ0EFF_MEBN,       & 
                 PALBNIR_TVEG, PALBVIS_TVEG,PALBNIR_TSOIL, PALBVIS_TSOIL, &
                 PABC, PIACAN, PPOI, PCSP, PRESP_BIOMASS_INST,  PPALPHAN, &
                 ZF2, PLW_RAD, ZGRNDFLUX, ZFLSN_COR, PUSTAR, ZEMIST,  &
                 PHU_AGG, PAC_AGG, ZDELHEATV_SFC, ZDELHEATG_SFC, ZDELHEATG, &
                 ZDELHEATN, ZDELHEATN_SFC, ZGSFCSNOW, PTDEEP_A, PDEEP_FLUX, &
                 ZRI3L, ZSNOW_THRUFAL, ZEVAPCOR, ZSUBVCOR, ZSNOWSFCH, ZQS3L   )

ELSE
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      5.0    Radiative transfert
!              -------------------
!
  IF (IO%LTR_ML) THEN
    CALL RADIATIVE_TRANSFERT(IO%LAGRI_TO_GRASS, IP%XVEGTYPE_PATCH(:,:,1), PALBVIS_TVEG, &
                             PALBVIS_TSOIL, PALBNIR_TVEG, PALBNIR_TSOIL, PSW_RAD,       &
                             IMT%XLAI(:,1), PZENITH, PABC, IR%XFAPARC(:,1),             &
                             IR%XFAPIRC(:,1), IR%XMUS(:,1), IR%XLAI_EFFC(:,1), GSHADE,  &
                             PIACAN, ZIACAN_SUNLIT, ZIACAN_SHADE, ZFRAC_SUN,            &
                             DGMI%XFAPAR, DGMI%XFAPIR, DGMI%XFAPAR_BS, DGMI%XFAPIR_BS  )
   ENDIF
!
!print*,'tg4 ',IR%XTG(1,1,1)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      6.0    Fraction of leaves occupied by intercepted water
!              ------------------------------------------------
!
   CALL WET_LEAVES_FRAC(IR%XWR(:,1), IMT%XVEG(:,1), IMT%XWRMAX_CF(:,1), DGIP%XZ0, &
                        IMT%XLAI(:,1), ZWRMAX, ZDELTA)
!
!print*,'tg5 ',IR%XTG(1,1,1)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      7.0    Explicit snow scheme
!              --------------------
!
   CALL SNOW3L_ISBA(IO, G, IR, DGIP, DGEIP, DGMI, HSNOW_ISBA, OMEB, HIMPLICIT_WIND,    &
                    TPTIME, PTSTEP, IP%XVEGTYPE_PATCH(:,:,1), IR%XTG(:,:,1), DGMI%XCT, &
                    ZSOILHCAPZ, ZSOILCONDZ(:,1), PPS, PTA, PSW_RAD, PQA, PVMOD,        &
                    PLW_RAD, PRR, PSR, PRHOA, PUREF, PEXNS, PEXNA, PDIRCOSZW, PZREF,   &
                    IR%XSNOWFREE_ALB(:,1), IMX%XDG(:,:,1), IP%XDZG(:,:,1), PPEW_A_COEF,&
                    PPEW_B_COEF, PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,   &
                    ZSNOW_THRUFAL, ZGRNDFLUX, ZFLSN_COR, ZGSFCSNOW, ZEVAPCOR, ZLES3L,  &
                    ZLEL3L, ZEVAP3L, ZSNOWSFCH, ZDELHEATN, ZDELHEATN_SFC, ZRI3L,       &
                    PZENITH, ZDELHEATG, ZDELHEATG_SFC, ZQS3L                   ) 
!print*,'tg6 ',IR%XTG(1,1,1)            
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      8.0    Plant stress, stomatal resistance and, possibly, CO2 assimilation
!              --------------------------------------------------------------------
!
!print*,'HPHOTO ',IO%CPHOTO
   IF (IO%CPHOTO=='NON') THEN
      CALL VEG(PSW_RAD, PTA, PQA, PPS, IMT%XRGL(:,1), IMT%XLAI(:,1), &
                IMT%XRSMIN(:,1), IMT%XGAMMA(:,1), ZF2, DGMI%XRS)
   ELSE IF (MAXVAL(IMT%XGMES(:,1)).NE.XUNDEF .OR. MINVAL(IMT%XGMES(:,1)).NE.XUNDEF) THEN
      ZQSAT(:)=QSAT(IR%XTG(:,1,1),PPS(:))  
      CALL COTWORES(PTSTEP, IO, GSHADE, IP, IMT, IR, IMX%XDMAX(:,1), PPOI, PCSP,    &
                    IR%XTG(:,1,1), ZF2, PSW_RAD, PQA, ZQSAT, IR%XPSNV(:,1), ZDELTA, &
                    PRHOA, PZENITH, INI%XFFV(:,1), ZIACAN_SUNLIT, ZIACAN_SHADE,     &
                    ZFRAC_SUN, PIACAN, PABC, DGMI%XRS, DGEIP%XGPP, PRESP_BIOMASS_INST(:,1))
   ELSE
      PRESP_BIOMASS_INST(:,1) = 0.0
      DGEIP%XGPP(:) = 0.0
   ENDIF
!
!print*,'tg7 ',IR%XTG(1,1,1)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      9.0    ISBA Composit Energy Budget
!              -----------------------------------------------
!
  CALL ISBA_CEB(IO, IP, INI, IMX, IMT, IR, DGIP, DGEIP, DGMI,       &
                HSNOW_ISBA, HIMPLICIT_WIND, PTSTEP, PPEW_A_COEF,    &
                PPEW_B_COEF, PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, &
                PPEQ_B_COEF, PSW_RAD, PLW_RAD, PEXNS, PEXNA, PTA,   &
                PVMOD, PQA, PRR, PSR, PPS, PZREF, PUREF, PDIRCOSZW, &
                ZF5, PFFG_NOSNOW, PFFV_NOSNOW,  PRHOA, ZCS,         &
                ZSOILCONDZ, ZSOILHCAPZ, ZFROZEN1, PTDEEP_A,         &
                ZGRNDFLUX, ZFLSN_COR, ZSNOW_THRUFAL, ZDELTA, ZHUGI, &
                ZALBT, ZEMIST, PDEEP_FLUX, PUSTAR, PAC_AGG, PHU_AGG )
!
!print*,'tg8 ',IR%XTG(1,1,1)
ENDIF
!
!*******************************************************************************
! WARNING: at this stage, all fluxes have two different meanings according
!          to the ISBA snow-scheme option:
!  'D95' : they represent aggregated (snow + flood + snow-flood-free) fluxes
!  '3-L' : they represent                    flood + snow-flood-free  fluxes
!
! The variables concerned by this are: PRN, PH, PLE, PLEI, DGEIP%XLEG, DGEIP%XLEGI, DGEIP%XLEV, DGEIP%XLES, 
!                                      DGEIP%XLER, DGEIP%XLETR, PEVAP, PUSTAR, PGFLUX
!*******************************************************************************
!
!*     12.0    Water transfers and phase change in the soil
!              --------------------------------------------
!
CALL HYDRO(IO, P, IP, INI, IMX, IMT, IMI, AG, IR, DGEIP, DGMI, &
           HSNOW_ISBA, OMEB, PTSTEP, ZVEG, ZWRMAX, ZSNOW_THRUFAL, &
           ZEVAPCOR, ZSUBVCOR, ZSOILHCAPZ, ZF2WGHT, ZF2, PPS,  &
           PIRRIG_GR, ZDELHEATG, ZDELHEATG_SFC,  ZDELPHASEG,   &
           ZDELPHASEG_SFC                                )
!print*,'tg9 ',IR%XTG(1,1,1)   
!-------------------------------------------------------------------------------
!
!*     13.0    Aggregated output fluxes and diagnostics
!              -----------------------------------------
!
!* add snow component to output radiative parameters and fluxes in case 
!  of ES or CROCUS snow schemes
!
CALL ISBA_SNOW_AGR(IP, INI, IR, DGI, DGEI, DGMI, DGIP, DGEIP, &
                   HSNOW_ISBA, OMEB, PEXNS, PEXNA, PTA, PQA,  &
                   PZREF, PUREF, PDIRCOSZW, PVMOD, PRR, PSR,  &
                   ZEMIST, ZALBT, PUSTAR, ZLES3L, ZLEL3L,     &
                   ZEVAP3L, ZQS3L, ZALB3L, ZGSFCSNOW,         &
                   ZGRNDFLUX, ZFLSN_COR, PEMIST, PPALPHAN    )  
!
!print*,'tg10 ',IR%XTG(1,1,1)
!***************************************************************************
! All output fluxes and radiative variables have recovered the same physical
! meaning, that is they are aggregated quantities (snow + snow-free)
!***************************************************************************
!
IF (LHOOK) CALL DR_HOOK('ISBA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ISBA
