!     #########
SUBROUTINE SNOW3L_ISBA(HISBA, HSNOW_ISBA, HSNOWRES, OMEB, OGLACIER, HIMPLICIT_WIND,          &
                         TPTIME, PTSTEP, PVEGTYPE,                                           &
                         PSNOWSWE, PSNOWHEAT, PSNOWRHO, PSNOWALB,                            &
                         PSNOWGRAN1, PSNOWGRAN2, PSNOWHIST,PSNOWAGE,PSNOWIMPUR,PSNOWIMPURV2, &
                         PTG, PCG, PCT, PSOILHCAPZ, PSOILCONDZ,                              &
                         PPS, PTA, PSW_RAD, PQA, PVMOD, PVDIR, PLW_RAD, PRR, PSR,            &
                         PRHOA, PUREF, PEXNS, PEXNA, PDIRCOSZW,PSLOPEDIR,                    &
                         PZREF, PZ0NAT, PZ0EFF, PZ0HNAT, PALB, PD_G, PDZG,                   &
                         PPEW_A_COEF, PPEW_B_COEF,                                           &
                         PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,                 &
                         PTHRUFAL, PGRNDFLUX, PFLSN_COR, PGSFCSNOW, PEVAPCOR,                &
                         PSWNETSNOW, PSWNETSNOWS, PLWNETSNOW,                                &
                         PRNSNOW, PHSNOW, PGFLUXSNOW, PHPSNOW, PLES3L, PLEL3L, PEVAP,        &
                         PSNDRIFT, PUSTARSNOW, PPSN, PSRSFC, PRRSFC, PSNOWSFCH,              &
                         PDELHEATN, PDELHEATN_SFC,                                           &
                         PEMISNOW, PCDSNOW, PCHSNOW, PSNOWTEMP, PSNOWLIQ, PSNOWDZ,           &
                         PSNOWHMASS, PRI, PZENITH, PDELHEATG, PDELHEATG_SFC, PLAT, PLON, PQS,&
                         OSNOWDRIFT,OSNOWDRIFT_SUBLIM,OSNOW_ABS_ZENITH,                      &
                         HSNOWMETAMO, HSNOWRAD,OATMORAD,OSNOWSYTRON,KTAB_SYT,PSYTMASS,                &
                         PSNOWDEND,PSNOWSPHER,PSNOWSIZE,PSNOWSSA,PSNOWTYPEMEPRA,PSNOWRAM,    &
                         PSNOWSHEAR,PSNOWDEPTH_1DAYS,PSNOWDEPTH_3DAYS,PSNOWDEPTH_5DAYS,      &
                         PSNOWDEPTH_7DAYS,PSNOWSWE_1DAYS,PSNOWSWE_3DAYS,PSNOWSWE_5DAYS,      &
                         PSNOWSWE_7DAYS,PSNOWRAM_SONDE,PSNOW_WETTHICKNESS,PSNOW_REFROZENTHICKNESS,&
                         P_DIR_SW, P_SCA_SW, PSPEC_ALB, PDIFF_RATIO)                               
!     ######################################################################################
!
!!****  *SNOW3L_ISBA*  
!!
!!    PURPOSE
!!    -------
!
!     3-Layer snow scheme option (Boone and Etchevers 1999)
!     This routine is NOT called as snow depth goes below
!     a critical threshold which is vanishingly small.
!     This routine acts as an interface between SNOW3L and ISBA.
!     
!!**  METHOD
!!    ------
!
!     Direct calculation
!
!!    EXTERNAL
!!    --------
!
!     None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!    Boone and Etchevers (1999)
!!    Belair (1995)
!!    Noilhan and Planton (1989)
!!    Noilhan and Mahfouf (1996)
!!      
!!    AUTHOR
!!    ------
!!	A. Boone           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        7/99  Boone
!!      Packing added   4/00  Masson & Boone
!!      z0h and snow    2/06  LeMoigne
!!
!!      Modified by B. Decharme  (03/2009): Consistency with Arpege permanent
!!                                          snow/ice treatment
!!      Modified by A. Boone     (04/2010): Implicit coupling with atmosphere permitted.
!!
!!      Modified by B. Decharme  (04/2010): check suspicious low temperature for ES and CROCUS
!!      Modified by B. Decharme  (08/2013): Qsat as argument (needed for coupling with atm)
!!      Modified by A. Boone     (10/2014): MEB: pass in fluxes when using MEB
!!      Modified by M. Lafaysse (08/2015): MEB-Crocus coupling
!!
!-------------------------------------------------------------------------------
!
USE MODD_CSTS,       ONLY : XTT, XPI, XDAY, XLMTT, XLSTT
USE MODD_SNOW_PAR,   ONLY : XRHOSMAX_ES, XSNOWDMIN, XRHOSMIN_ES, XEMISSN
USE MODD_PREP_SNOW,   ONLY : NIMPUR
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_TYPE_DATE_SURF, ONLY: DATE_TIME
!
USE MODD_DATA_COVER_PAR, ONLY : NVT_SNOW
!
USE MODI_SNOW3L
USE MODI_SNOWCRO
USE MODI_SNOWCRO_DIAG
USE MODI_SNOW_SYTRON
!
#ifdef SFX_OL
USE MODN_IO_OFFLINE, ONLY : XTSTEP_OUTPUT
#endif
!
USE MODI_ABOR1_SFX
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP
!                                      PTSTEP    = time step of the integration
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PVEGTYPE ! fraction of each vegetation
!
 CHARACTER(LEN=*),     INTENT(IN)    :: HISBA
!                                      HISBA     = FLAG to use Force-Restore or DIFfusion
!                                      soil heat and mass transfer method
!
 CHARACTER(LEN=*),     INTENT(IN)    :: HSNOW_ISBA
!                                      HSNOW_ISBA = FLAG to use SNOW3L or not 
!                                      (or default FR method)
!
 CHARACTER(LEN=*),     INTENT(IN)    :: HSNOWRES
!                                      HSNOWRES  = ISBA-SNOW3L turbulant exchange option
!                                      'DEF' = Default: Louis (ISBA: Noilhan and Mahfouf 1996)
!                                      'RIL' = Limit Richarson number under very stable 
!                                              conditions (currently testing)
!
LOGICAL, INTENT(IN)                 :: OMEB       ! True = coupled to MEB. This means surface fluxes ae IMPOSED
!                                                 ! as an upper boundary condition to the explicit snow schemes. 
!                                                 ! If = False, then energy
!                                                 ! budget and fluxes are computed herein.
!
LOGICAL, INTENT(IN)                 :: OGLACIER   ! True = Over permanent snow and ice, 
!                                                     initialise WGI=WSAT,
!                                                     Hsnow>=10m and allow 0.8<SNOALB<0.85
!
CHARACTER(LEN=*),     INTENT(IN)    :: HIMPLICIT_WIND   ! wind implicitation option
!                                                       ! 'OLD' = direct
!                                                       ! 'NEW' = Taylor serie, order 1
!
TYPE(DATE_TIME), INTENT(IN)         :: TPTIME     ! current date and time
!
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSOILHCAPZ, PD_G, PDZG
REAL, DIMENSION(:),   INTENT(IN)    :: PCG, PCT, PSOILCONDZ  
!                                      PD_G      = Depth to bottom of each soil layer (m)
!                                      PDZG      = Soil layer thicknesses (m)
!                                      PCG       = area-averaged soil heat capacity [(K m2)/J]
!                                      PCT       = area-averaged surface heat capacity [(K m2)/J]
!                                      PSOILCONDZ= soil thermal conductivity (W m-1 K-1)
!                                      PSOILHCAPZ= soil heat capacity (J m-3 K-1)
!
REAL, DIMENSION(:), INTENT(IN)      :: PPS, PTA, PSW_RAD, PQA,                       &
                                       PVMOD, PVDIR, PLW_RAD, PSR, PRR  
!                                      PSW_RAD = incoming solar radiation (W/m2)
!                                      PLW_RAD = atmospheric infrared radiation (W/m2)
!                                      PRR     = rain rate [kg/(m2 s)]
!                                      PSR     = snow rate (SWE) [kg/(m2 s)]
!                                      PTA     = atmospheric temperature at level za (K)
!                                      PVMOD   = modulus of the wind parallel to the orography (m/s)
!                                      PVDIR   = wind direction (rad)
!                                      PPS     = surface pressure
!                                      PQA     = atmospheric specific humidity
!                                                at level za
REAL, DIMENSION(:,:), INTENT(IN)    :: P_DIR_SW, P_SCA_SW
!					P_DIR_SW, P_SCA_SW = direct and diffuse spectral solar irradiance (W/m2/um)
!
REAL, DIMENSION(:), INTENT(IN)      :: PZREF, PUREF, PEXNS, PEXNA, PDIRCOSZW, PRHOA, PZ0NAT, &
                                       PZ0EFF, PZ0HNAT, PALB,PSLOPEDIR
!                                      PZ0EFF    = roughness length for momentum 
!                                      PZ0NAT    = grid box average roughness length
!                                      PZ0HNAT   = grid box average roughness length
!                                      PZREF     = reference height of the first
!                                                  atmospheric level
!                                      PUREF     = reference height of the wind
!                                      PRHOA     = air density
!                                      PEXNS     = Exner function at surface
!                                      PEXNA     = Exner function at lowest atmos level
!                                      PDIRCOSZW = Cosinus of the angle between the 
!                                                  normal to the surface and the vertical
!                                      PALB      = soil/vegetation albedo
!                                      PSLOPE_DIR = direction of the slope
!
REAL, DIMENSION(:), INTENT(IN)      :: PPSN
!                                      PPSN     = Snow cover fraction (total) 
!
REAL, DIMENSION(:), INTENT(IN)      :: PPEW_A_COEF, PPEW_B_COEF,                   &
                                       PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF,      &
                                       PPEQ_B_COEF  
!                                      PPEW_A_COEF = wind coefficient
!                                      PPEW_B_COEF = wind coefficient
!                                      PPET_A_COEF = A-air temperature coefficient
!                                      PPET_B_COEF = B-air temperature coefficient
!                                      PPEQ_A_COEF = A-air specific humidity coefficient
!                                      PPEQ_B_COEF = B-air specific humidity coefficient                         !
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PTG
!                                      PTG       = Soil temperature profile (K)
!
INTEGER , DIMENSION(:), INTENT(IN)   ::  KTAB_SYT    ! Array of index defining
!                                             opposite points for Sytron
REAL, DIMENSION(:), INTENT(INOUT)   :: PSNOWALB
!                                      PSNOWALB = Prognostic surface snow albedo
!                                                 (does not include anything but
!                                                 the actual snow cover)
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWHEAT, PSNOWRHO, PSNOWSWE
!                                      PSNOWHEAT = Snow layer(s) heat content (J/m3)
!                                      PSNOWRHO  = Snow layer(s) averaged density (kg/m3)
!                                      PSNOWSWE  = Snow layer(s) liquid Water Equivalent (SWE:kg m-2)
!
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWGRAN1, PSNOWGRAN2, PSNOWHIST
!                                      PSNOWGRAN1 = Snow layer(s) grain parameter 1
!                                      PSNOWGRAN2 = Snow layer(s) grain parameter 2
!                                      PSNOWHIST  = Snow layer(s) grain historical parameter
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWAGE  ! Snow grain age
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWIMPUR  ! Snow impurity content N6K
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PSNOWIMPURV2  ! Snow impurity content N6K
!
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PRNSNOW, PHSNOW, PLES3L, PLEL3L,     &
                                       PHPSNOW, PEMISNOW, PEVAP, PGRNDFLUX, PSWNETSNOW, &
                                       PLWNETSNOW, PSWNETSNOWS, PDELHEATG, PDELHEATG_SFC
!                                      PLEL3L        = evaporation heat flux from snow (W/m2)
!                                      PLES3L        = sublimation (W/m2)
!                                      PHPSNOW       = heat release from rainfall (W/m2)
!                                      PRNSNOW       = net radiative flux from snow (W/m2)
!                                      PHSNOW        = sensible heat flux from snow (W/m2)
!                                      PEMISNOW      = snow surface emissivity
!                                      PEVAP         = total evaporative flux from snow (kg/m2/s)
!                                      PGRNDFLUX     = soil/snow interface heat flux (W/m2)
!                                      PSWNETSNOW    = net shortwave radiation entering top of snowpack (W/m2)
!                                      PSWNETSNOWS   = net shortwave radiation in uppermost layer of snowpack (W/m2)
!                                                      (for surface energy budget closure diagnostics)
!                                      PLWNETSNOW    = net longwave radiation entering top of snowpack (W/m2)
!                                      PDELHEATG     = ground heat content change (diagnostic) (W/m2)
!                                                      note, modified if ground-snow flux adjusted
!                                      PDELHEATG_SFC = ground heat content change in sfc only (diagnostic) (W/m2)
!                                                      note, modified if ground-snow flux adjusted
!
!
REAL, DIMENSION(:), INTENT(OUT)     :: PGFLUXSNOW
!                                      PGFLUXSNOW    = net heat flux from snow (W/m2)
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PUSTARSNOW, PCDSNOW, PCHSNOW, PRI
!                                      PCDSNOW    = drag coefficient for momentum over snow (-)
!                                      PUSTARSNOW = friction velocity over snow (m/s)
!                                      PCHSNOW    = drag coefficient for heat over snow (-)
!                                      PRI        = Richardson number (-)
!
REAL, DIMENSION(:), INTENT(OUT)     :: PTHRUFAL, PFLSN_COR, PEVAPCOR, PSNOWHMASS, PGSFCSNOW
!                                      PTHRUFAL  = rate that liquid water leaves snow pack: 
!                                                  paritioned into soil infiltration/runoff 
!                                                  by ISBA [kg/(m2 s)]
!                                      PFLSN_COR = soil/snow correction heat flux (W/m2) (not MEB)
!                                      PEVAPCOR  = evaporation/sublimation correction term:
!                                                  extract any evaporation exceeding the
!                                                  actual snow cover (as snow vanishes)
!                                                  and apply it as a surface soil water
!                                                  sink. [kg/(m2 s)]
!                                      PSNOWHMASS = heat content change due to mass
!                                                   changes in snowpack (J/m2): for budget
!                                                   calculations only.
!                                      PGSFCSNOW  = heat flux between the surface and sub-surface 
!                                                   snow layers (for energy budget diagnostics) (W/m2)
!
REAL, DIMENSION(:), INTENT(OUT)     :: PSNDRIFT
!                                      PSNDRIFT    = blowing snow sublimation (kg/m2/s)
REAL, DIMENSION(:), INTENT(OUT)     :: PSYTMASS 
!                                      PSYTMASS    = eroded/cumulated snow mass in SYTRON (kg/m2/s)
!
REAL, DIMENSION(:,:), INTENT(OUT) :: PSNOWDEND
REAL, DIMENSION(:,:), INTENT(OUT) :: PSNOWSPHER
REAL, DIMENSION(:,:), INTENT(OUT) :: PSNOWSIZE
REAL, DIMENSION(:,:), INTENT(OUT) :: PSNOWSSA
REAL, DIMENSION(:,:), INTENT(OUT) :: PSNOWTYPEMEPRA
REAL, DIMENSION(:,:), INTENT(OUT) :: PSNOWRAM
REAL, DIMENSION(:,:), INTENT(OUT) :: PSNOWSHEAR
REAL, DIMENSION(:), INTENT(OUT) :: PSNOWDEPTH_1DAYS
REAL, DIMENSION(:), INTENT(OUT) :: PSNOWDEPTH_3DAYS
REAL, DIMENSION(:), INTENT(OUT) :: PSNOWDEPTH_5DAYS
REAL, DIMENSION(:), INTENT(OUT) :: PSNOWDEPTH_7DAYS
REAL, DIMENSION(:), INTENT(OUT) :: PSNOWSWE_1DAYS
REAL, DIMENSION(:), INTENT(OUT) :: PSNOWSWE_3DAYS
REAL, DIMENSION(:), INTENT(OUT) :: PSNOWSWE_5DAYS
REAL, DIMENSION(:), INTENT(OUT) :: PSNOWSWE_7DAYS
REAL, DIMENSION(:), INTENT(OUT) :: PSNOWRAM_SONDE
REAL, DIMENSION(:), INTENT(OUT) :: PSNOW_WETTHICKNESS
REAL, DIMENSION(:), INTENT(OUT) :: PSNOW_REFROZENTHICKNESS
REAL, DIMENSION(:,:), INTENT(OUT) :: PSPEC_ALB, PDIFF_RATIO ! spectral albedo and diffuse to total irradiance ratio



REAL, DIMENSION(:), INTENT(OUT)     :: PSRSFC, PRRSFC, PSNOWSFCH, PDELHEATN, PDELHEATN_SFC
!                                      PSRSFC = snow rate on soil/veg surface when SNOW3L in use
!                                      PRRSFC = rain rate on soil/veg surface when SNOW3L in use
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWTEMP
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWLIQ, PSNOWDZ
!                                      PSNOWLIQ  = Snow layer(s) liquid water content (m)
!                                      PSNOWTEMP = Snow layer(s) temperature (m)
!                                      PSNOWDZ   = Snow layer(s) thickness (m)
!
REAL, DIMENSION(:), INTENT(OUT)     :: PQS
!                                      PQS = surface humidity (kg/kg)
!
! ajout_EB pour prendre en compte angle zenithal du soleil dans LRAD
! puis plus tard dans LALB
REAL, DIMENSION(:), INTENT(IN)      :: PZENITH    ! solar zenith angle
REAL, DIMENSION(:), INTENT(IN)      :: PLAT
REAL, DIMENSION(:), INTENT(IN)      :: PLON
!
LOGICAL, INTENT(IN)                 :: OSNOWDRIFT, OSNOWDRIFT_SUBLIM ! activate snowdrift, sublimation during drift
LOGICAL, INTENT(IN)                 :: OSNOW_ABS_ZENITH ! activate parametrization of solar absorption for polar regions
CHARACTER(3), INTENT(IN)            :: HSNOWMETAMO, HSNOWRAD
                                         !-----------------------
                                         ! Crocus metamorphism scheme
                                         ! HSNOWMETAMO=B92 Brun et al 1992
                                         ! HSNOWMETAMO=C13 Carmagnola et al 2014
                                         ! HSNOWMETAMO=T07 Taillandier et al 2007
                                         ! HSNOWMETAMO=F06 Flanner et al 2006
                                         !-----------------------
                                         ! Crocus radiative transfer scheme
                                         ! HSNOWMETAMO=B92 Brun et al 1992
                                         ! HSNOWMETAMO=TAR TARTES (Libois et al 2013)
                                         ! HSNOWMETAMO=TA1 TARTES with constant impurities
                                         ! HSNOWMETAMO=TA2 TARTES with constant impurities as function of ageing
LOGICAL, INTENT(IN)                 :: OSNOWSYTRON ! activate SYTROn snow redistribution scheme
LOGICAL, INTENT(IN)                 :: OATMORAD ! activate atmotartes scheme

!*      0.2    declarations of local variables
!
REAL, PARAMETER                     :: ZCHECK_TEMP = 100.0 
!                                      Limit to check suspicious low temperature (K)
!
INTEGER                             :: JWRK, JJ,JIMP ! Loop control
!
INTEGER                             :: INLVLS   ! maximum number of snow layers
INTEGER                             :: INLVLG   ! number of ground layers
INTEGER                             :: IBLOWSNW     ! number of blowing snow variables
!
REAL, DIMENSION(SIZE(PTA))          :: ZRRSNOW, ZSOILCOND, ZSNOW, ZSNOWFALL,  &
                                       ZSNOWABLAT_DELTA, ZSNOWSWE_1D, ZSNOWD, & 
                                       ZSNOWH, ZSNOWH1, ZGRNDFLUXN, ZPSN,     &
                                       ZSOILCOR, ZSNOWSWE_OUT, ZTHRUFAL,      &
                                       ZSNOW_MASS_BUDGET
!                                      ZSOILCOND    = soil thermal conductivity [W/(m K)]
!                                      ZRRSNOW      = rain rate over snow [kg/(m2 s)]
!                                      ZSNOW        = snow depth (m) 
!                                      ZSNOWFALL    = minimum equivalent snow depth
!                                                     for snow falling during the
!                                                     current time step (m)
!                                      ZSNOWABLAT_DELTA = FLAG =1 if snow ablates completely
!                                                     during current time step, else=0
!                                      ZSNOWSWE_1D  = TOTAL snowpack SWE (kg m-2)
!                                      ZSNOWD       = snow depth
!                                      ZSNOWH       = snow total heat content (J m-2)
!                                      ZSNOWH1      = snow surface layer heat content (J m-2)
!                                      ZGRNDFLUXN   = corrected snow-ground flux (if snow fully ablated during timestep)
!                                      ZPSN         = snow fraction working array
!                                      ZSOILCOR = for vanishingy thin snow cover,
!                                                 allow any excess evaporation
!                                                 to be extracted from the soil
!                                                 to maintain an accurate water
!                                                 balance [kg/(m2 s)]
!                                      ZSNOW_MASS_BUDGET = snow water equivalent budget (kg/m2/s)

REAL, DIMENSION(SIZE(PTA),4)        :: ZBLOWSNW   ! Properties of deposited blowing snow
                                      !    1 : Deposition flux (kg/m2/s)
                                      !    2 : Density of deposited snow (kg/m3)
                                      !    3 : SGRA1 of deposited snow
                                      !    4 : SGRA2 of deposited snow
REAL, DIMENSION(SIZE(PTA))          :: ZBLOWSNW_ACC
!                                      ZBLOWSNW_ACC  = minimum equivalent snow depth
!                                                      for deposition of blown snow particles
!                                                      during the current time step (m)
!
!*      0.3    declarations of packed  variables
!
INTEGER                            :: ISIZE_SNOW ! number of points where computations are done
INTEGER, DIMENSION(SIZE(PTA))      :: NMASK      ! indices correspondance between arrays
!
LOGICAL, DIMENSION(SIZE(PTA))      :: LREMOVE_SNOW

LOGICAL :: GCOMPUTECRODIAG ! flag to compute Crocus-MEPRA diagnostics

!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! - - ---------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3L_ISBA',0,ZHOOK_HANDLE)
!
!*       0.     Initialize variables:
!               ---------------------
!
IF (SIZE(PSNOWDEND)>1) THEN
  PSNOWDEND(:,:)=XUNDEF
  PSNOWSPHER(:,:)=XUNDEF
  PSNOWSIZE(:,:)=XUNDEF
  PSNOWSSA(:,:)=XUNDEF
  PSNOWTYPEMEPRA(:,:)=XUNDEF
  PSNOWRAM(:,:)=XUNDEF
  PSNOWSHEAR(:,:)=XUNDEF
ENDIF

PFLSN_COR(:)   = 0.0
PTHRUFAL(:)    = 0.0
PEVAPCOR(:)    = 0.0
PSNDRIFT(:)    = 0.0
PSNOWHMASS(:)  = 0.0
PSRSFC(:)      = PSR(:)         ! these are snow and rain rates passed to ISBA,
PRRSFC(:)      = PRR(:)         ! so initialize here if SNOW3L not used:
PQS(:)         = XUNDEF
PSYTMASS(:)    = 0.0
!
ZSNOW(:)       = 0.0
ZSNOWD(:)      = 0.0
ZGRNDFLUXN(:)  = 0.0
ZSNOWH(:)      = 0.0
ZSNOWH1(:)     = 0.0
ZSNOWSWE_1D(:) = 0.0
ZSNOWSWE_OUT(:)= 0.0
ZSOILCOND(:)   = 0.0
ZRRSNOW(:)     = 0.0
ZSNOWFALL(:)   = 0.0
ZSNOWABLAT_DELTA(:) = 0.0
PSNOWLIQ(:,:)  = 0.0
PSNOWDZ(:,:)   = 0.0
ZBLOWSNW(:,:)  = 0.0
ZBLOWSNW_ACC(:)  = 0.0
!
INLVLS          = SIZE(PSNOWSWE(:,:),2)                         
INLVLG          = MIN(SIZE(PD_G(:,:),2),SIZE(PTG(:,:),2))                         
IBLOWSNW       = SIZE(ZBLOWSNW(:,:),2)
! PSPEC_ALB(:,:)=0.
! PDIFF_RATIO(:,:)=0.
!
IF(.NOT.OMEB)THEN 
!
! If MEB activated, these values are input, else initialize here:
!
   PGRNDFLUX(:)   = 0.0
   PLES3L(:)      = 0.0
   PLEL3L(:)      = 0.0
   PEVAP(:)       = 0.0
   PRNSNOW(:)     = 0.0
   PHSNOW(:)      = 0.0
   PGFLUXSNOW(:)  = 0.0
   PHPSNOW(:)     = 0.0
   PSWNETSNOW(:)  = 0.0
   PSWNETSNOWS(:) = 0.0
   PLWNETSNOW(:)  = 0.0
   PEMISNOW(:)    = XEMISSN
   PUSTARSNOW(:)  = 0.0
   PCDSNOW(:)     = 0.0
   PCHSNOW(:)     = 0.0
   PRI(:)         = XUNDEF
END IF
!
!
!GSYTRON = .NOT.( ALL (KTAB_SYT(:)==-999))
!
! Use ISBA-SNOW3L or NOT: NOTE that if explicit soil diffusion method in use,
! then *must* use explicit snow model:
!
IF (HSNOW_ISBA=='3-L' .OR. HISBA == 'DIF' .OR. HSNOW_ISBA == 'CRO') THEN
!
! - Snow and rain falling onto the 3-L grid space:
!
   PSRSFC(:)=0.0
!
   DO JJ=1,SIZE(PSR)
      ZRRSNOW(JJ)        = PPSN(JJ)*PRR(JJ)
      PRRSFC(JJ)         = PRR(JJ) - ZRRSNOW(JJ)
      ZSNOWFALL(JJ)      = PSR(JJ)*PTSTEP/XRHOSMAX_ES    ! maximum possible snowfall depth (m)
   ENDDO
!
! Calculate preliminary snow depth (m)
! 
!    ZSNOW(:)=0.
!    ZSNOWH(:)=0.
!    ZSNOWSWE_1D(:)=0.
!    ZSNOWH1(:)              = PSNOWHEAT(:,1)*PSNOWSWE(:,1)/PSNOWRHO(:,1) ! sfc layer only
! !
!    DO JWRK=1,SIZE(PSNOWSWE,2)
!       DO JJ=1,SIZE(PSNOWSWE,1)
!          ZSNOWSWE_1D(JJ)     = ZSNOWSWE_1D(JJ) + PSNOWSWE(JJ,JWRK)
!          ZSNOW(JJ)           = ZSNOW(JJ)       + PSNOWSWE(JJ,JWRK)/PSNOWRHO(JJ,JWRK)
!          ZSNOWH(JJ)          = ZSNOWH(JJ)      + PSNOWHEAT(JJ,JWRK)*PSNOWSWE(JJ,JWRK)/PSNOWRHO(JJ,JWRK)
!       END DO
!    ENDDO
! !
   IF(HISBA == 'DIF')THEN
      ZSOILCOND(:)   = PSOILCONDZ(:)
   ELSE
!
! - Soil thermal conductivity
!   is implicit in Force-Restore soil method, so it
!   must be backed-out of surface thermal coefficients
!   (Etchevers and Martin 1997):
!
      ZSOILCOND(:)    = 4.*XPI/( PCG(:)*PCG(:)*XDAY/(PD_G(:,1)*PCT(:)) )
!
   ENDIF
!
! ===============================================================
!        Snow redistribution scheme Sytron
! 
! 
IF (HSNOW_ISBA=='CRO' .AND. OSNOWSYTRON) THEN

  CALL SNOW_SYTRON(PTSTEP,PPS,PTA,PQA,PVMOD,PVDIR,PSLOPEDIR,PDIRCOSZW,     &
                        PSNOWHEAT,PSNOWSWE,PSNOWRHO,                       &
                        PSNOWGRAN1,PSNOWGRAN2,PSNOWHIST,PSNOWAGE,KTAB_SYT, &
                        ZBLOWSNW,PSYTMASS)
!
! Calculate maximum snow depth (m) of deposited blown snow particles
!
  WHERE(ZBLOWSNW(:,1)> 0.)
     ZBLOWSNW_ACC(:)=ZBLOWSNW(:,1)*PTSTEP/ZBLOWSNW(:,2)
  END WHERE
ENDIF
!
! Calculate preliminary snow depth (m)
! Now after SNOW_SYTRON to account for cases of total snowpack erosion 
 ZSNOW(:)=0.
 ZSNOWH(:)=0.
 ZSNOWSWE_1D(:)=0.
 ZSNOWH1(:)              = PSNOWHEAT(:,1)*PSNOWSWE(:,1)/PSNOWRHO(:,1) ! sfc layer only




  DO JWRK=1,SIZE(PSNOWSWE,2)
    DO JJ=1,SIZE(PSNOWSWE,1)
     ZSNOWSWE_1D(JJ)     = ZSNOWSWE_1D(JJ) + PSNOWSWE(JJ,JWRK)
     ZSNOW(JJ)           = ZSNOW(JJ)       + PSNOWSWE(JJ,JWRK)/PSNOWRHO(JJ,JWRK)
     ZSNOWH(JJ)          = ZSNOWH(JJ)      + PSNOWHEAT(JJ,JWRK)*PSNOWSWE(JJ,JWRK)/PSNOWRHO(JJ,JWRK)
    END DO
  ENDDO

! ===============================================================
! === Packing: Only call snow model when there is snow on the surface
!              exceeding a minimum threshold OR if the equivalent
!              snow depth falling during the current time step exceeds 
!              this limit.
!
! counts the number of points where the computations will be made
!
!
   ISIZE_SNOW = 0
   NMASK(:) = 0
!
   DO JJ=1,SIZE(ZSNOW)
      IF (ZSNOW(JJ) >= XSNOWDMIN .OR. ZSNOWFALL(JJ) >= XSNOWDMIN .OR.ZBLOWSNW_ACC(JJ) >= XSNOWDMIN) THEN
         ISIZE_SNOW = ISIZE_SNOW + 1
         NMASK(ISIZE_SNOW) = JJ
      ENDIF
   ENDDO
!  
   IF (ISIZE_SNOW>0) CALL CALL_MODEL(ISIZE_SNOW,INLVLS,INLVLG,IBLOWSNW,NMASK)
!
! ===============================================================
!
! Remove trace amounts of snow and reinitialize snow prognostic variables
! if snow cover is ablated.
! If MEB used, soil T already computed, therefore correct heating/cooling
! effect of updated snow-soil flux
!
   ZSNOWD(:) = 0.
   ZSNOWSWE_OUT(:) = 0.
   DO JWRK=1,SIZE(PSNOWSWE,2)
      DO JJ=1,SIZE(PSNOWSWE,1)
         ZSNOWD      (JJ) = ZSNOWD      (JJ) + PSNOWSWE(JJ,JWRK)/PSNOWRHO(JJ,JWRK)
         ZSNOWSWE_OUT(JJ) = ZSNOWSWE_OUT(JJ) + PSNOWSWE(JJ,JWRK)
      ENDDO
   END DO
!
   LREMOVE_SNOW(:)=(ZSNOWD(:)<XSNOWDMIN*1.1)
!
!
   IF(OMEB)THEN
     ZPSN(:)=1.0
   ELSE
!    To Conserve mass in ISBA without MEB, 
!    EVAP must be weignted by the snow fraction
!    in the calulation of THRUFAL
     ZPSN(:)=PPSN(:)
   ENDIF
!
   ZSNOWABLAT_DELTA(:) = 0.0
   ZTHRUFAL        (:) = PTHRUFAL(:)
!
   WHERE(LREMOVE_SNOW(:))
      ZSNOWSWE_OUT(:)     = 0.0
      PLES3L(:)           = MIN(PLES3L(:), XLSTT*(ZSNOWSWE_1D(:)/PTSTEP + PSR(:)))
      PLEL3L(:)           = 0.0
      PEVAP(:)            = PLES3L(:)/XLSTT
      PTHRUFAL(:)         = MAX(0.0, ZSNOWSWE_1D(:)/PTSTEP + PSR(:) - PEVAP(:)*ZPSN(:) + ZRRSNOW(:)) ! kg m-2 s-1
      ZTHRUFAL(:)         = MAX(0.0, ZSNOWSWE_1D(:)/PTSTEP + PSR(:) - PEVAP(:)         + ZRRSNOW(:)) ! kg m-2 s-1
      PSRSFC(:)           = 0.0
      PRRSFC(:)           = PRRSFC(:)
      ZSNOWABLAT_DELTA(:) = 1.0
      PSNOWALB(:)         = XUNDEF
      PEVAPCOR(:)         = 0.0
      ZSOILCOR(:)         = 0.0
      PGFLUXSNOW(:)       = PRNSNOW(:) - PHSNOW(:) - PLES3L(:) - PLEL3L(:)
      PSNOWHMASS(:)       = -PSR(:)*(XLMTT*PTSTEP)
      PGSFCSNOW(:)        = 0.0
      PDELHEATN(:)        = -ZSNOWH(:) /PTSTEP
      PDELHEATN_SFC(:)    = -ZSNOWH1(:)/PTSTEP
      PSNOWSFCH(:)        = PDELHEATN_SFC(:) - (PSWNETSNOWS(:) + PLWNETSNOW(:)    &
                          - PHSNOW(:) - PLES3L(:) - PLEL3L(:)) + PGSFCSNOW(:)     &
                          - PSNOWHMASS(:)/PTSTEP 
      ZGRNDFLUXN(:)       = (ZSNOWH(:)+PSNOWHMASS(:))/PTSTEP + PGFLUXSNOW(:)
      PTG(:,1)            = PTG(:,1) + PTSTEP*PCT(:)*ZPSN(:)*(ZGRNDFLUXN(:) - PGRNDFLUX(:) - PFLSN_COR(:))
      PDELHEATG(:)        = PDELHEATG(:)     + ZPSN(:)*(ZGRNDFLUXN(:) - PGRNDFLUX(:) - PFLSN_COR(:))
      PDELHEATG_SFC(:)    = PDELHEATG_SFC(:) + ZPSN(:)*(ZGRNDFLUXN(:) - PGRNDFLUX(:) - PFLSN_COR(:))
      PGRNDFLUX(:)        = ZGRNDFLUXN(:)
      PFLSN_COR(:)        = 0.0
     END WHERE
!
!
   DO JWRK=1,INLVLS
      DO JJ=1,SIZE(PSNOWSWE,1)
         PSNOWSWE (JJ,JWRK)  = (1.0-ZSNOWABLAT_DELTA(JJ))*PSNOWSWE(JJ,JWRK)
         PSNOWHEAT(JJ,JWRK)  = (1.0-ZSNOWABLAT_DELTA(JJ))*PSNOWHEAT(JJ,JWRK)
         PSNOWRHO (JJ,JWRK)  = (1.0-ZSNOWABLAT_DELTA(JJ))*PSNOWRHO(JJ,JWRK)  + &
                                    ZSNOWABLAT_DELTA(JJ)*XRHOSMIN_ES  
         PSNOWTEMP(JJ,JWRK)  = (1.0-ZSNOWABLAT_DELTA(JJ))*PSNOWTEMP(JJ,JWRK) + &
                                    ZSNOWABLAT_DELTA(JJ)*XTT  
         PSNOWLIQ (JJ,JWRK)  = (1.0-ZSNOWABLAT_DELTA(JJ))*PSNOWLIQ(JJ,JWRK)        
         PSNOWDZ  (JJ,JWRK)  = (1.0-ZSNOWABLAT_DELTA(JJ))*PSNOWDZ(JJ,JWRK)
         PSNOWAGE (JJ,JWRK)  = (1.0-ZSNOWABLAT_DELTA(JJ))*PSNOWAGE (JJ,JWRK)
         PSNOWIMPUR (JJ,JWRK)=(1.0-ZSNOWABLAT_DELTA(JJ))*PSNOWIMPUR(JJ,JWRK) !N6K
      ENDDO
   ENDDO
   DO JIMP=1,NIMPUR
     DO JWRK=1,INLVLS
        DO JJ=1,SIZE(PSNOWSWE,1)
           PSNOWIMPURV2 (JJ,JWRK,NIMPUR)=(1.0-ZSNOWABLAT_DELTA(JJ))*PSNOWIMPURV2(JJ,JWRK,NIMPUR) !N6K
        ENDDO
     ENDDO
   ENDDO
!  
   IF (HSNOW_ISBA=='CRO') THEN
      DO JWRK=1,INLVLS
         DO JJ=1,SIZE(PSNOWGRAN1,1)
            PSNOWGRAN1(JJ,JWRK)  = (1.0-ZSNOWABLAT_DELTA(JJ))*PSNOWGRAN1(JJ,JWRK) 
            PSNOWGRAN2(JJ,JWRK)  = (1.0-ZSNOWABLAT_DELTA(JJ))*PSNOWGRAN2(JJ,JWRK)
            PSNOWHIST(JJ,JWRK)   = (1.0-ZSNOWABLAT_DELTA(JJ))*PSNOWHIST(JJ,JWRK)
         ENDDO
      ENDDO
   ENDIF
!
!  ===============================================================
!
!  Compute snow mass budget 
!
   ZSNOW_MASS_BUDGET(:) = (ZSNOWSWE_1D(:)-ZSNOWSWE_OUT(:))/PTSTEP + PSR     (:)+ZRRSNOW (:) &
                                                                  - PEVAP   (:)-ZTHRUFAL(:) &
                                                                  + PEVAPCOR(:)+ZSOILCOR(:)
!
!
!  ===============================================================
!
!  To Conserve mass in ISBA, the latent heat flux part of 
!  the EVAPCOR term must be weignted by the snow fraction 
!
   PEVAPCOR (:) = PEVAPCOR(:)*ZPSN(:) + ZSOILCOR(:)
!
! ===============================================================
!
! check suspicious low temperature
!
   DO JWRK=1,INLVLS
      DO JJ=1,SIZE(PSNOWSWE,1)

         IF(PSNOWSWE(JJ,JWRK)>0.0) THEN
           IF (PSNOWTEMP(JJ,JWRK)<ZCHECK_TEMP)THEN
              WRITE(*,*) 'Suspicious low temperature :',PSNOWTEMP(JJ,JWRK)
              WRITE(*,*) 'At point and location      :',JJ,'LAT=',PLAT(JJ),'LON=',PLON(JJ)
              WRITE(*,*) 'At snow level / total layer:',JWRK,'/',INLVLS
              WRITE(*,*) 'SNOW MASS BUDGET (kg/m2/s) :',ZSNOW_MASS_BUDGET(JJ)
              WRITE(*,*) 'SWE BY LAYER      (kg/m2)  :',PSNOWSWE (JJ,1:INLVLS)
              WRITE(*,*) 'DEPTH BY LAYER      (m)    :',PSNOWDZ  (JJ,1:INLVLS)
              WRITE(*,*) 'DENSITY BY LAYER   (kg/m3) :',PSNOWRHO (JJ,1:INLVLS)
              WRITE(*,*) 'TEMPERATURE BY LAYER (K)   :',PSNOWTEMP(JJ,1:INLVLS)
              CALL ABOR1_SFX('SNOW3L_ISBA: Suspicious low temperature')                
            ENDIF
            
         ELSE
           !Prognostic variables forced to XUNDEF for correct outputs
           PSNOWDZ(JJ,JWRK)=XUNDEF
           ! Careful : to compute average surface temperature in ISBA_SNOW_AGR
           ! PSNOWTEMP(JJ,1) is required when PPSN(JJ)>0 even if PSNOWSWE(JJ,1)==0
           ! (vanishing snowpack)
           IF (.NOT.((PPSN(JJ)>0.0).AND.(JWRK==1))) THEN
               PSNOWTEMP(JJ,JWRK)=XUNDEF
           ENDIF
           PSNOWLIQ(JJ,JWRK)=XUNDEF
           PSNOWHEAT(JJ,JWRK)=XUNDEF
           PSNOWRHO(JJ,JWRK)=XUNDEF
           PSNOWAGE(JJ,JWRK)=XUNDEF

           IF (HSNOW_ISBA=='CRO') THEN
               PSNOWGRAN1(JJ,JWRK)=XUNDEF
               PSNOWGRAN2(JJ,JWRK)=XUNDEF
               PSNOWHIST(JJ,JWRK)=XUNDEF
           END IF
         ENDIF
         
      ENDDO
   ENDDO
!
! ===============================================================
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('SNOW3L_ISBA',1,ZHOOK_HANDLE)
!
CONTAINS
!
!================================================================
SUBROUTINE CALL_MODEL(KSIZE1,KSIZE2,KSIZE3,KSIZE4,KMASK)
!
IMPLICIT NONE
!
INTEGER, INTENT(IN) :: KSIZE1
INTEGER, INTENT(IN) :: KSIZE2
INTEGER, INTENT(IN) :: KSIZE3
INTEGER, INTENT(IN) :: KSIZE4
INTEGER, DIMENSION(:), INTENT(IN) :: KMASK
!
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWSWE
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWDZ
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWRHO
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWHEAT
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWTEMP
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWLIQ
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWGRAN1
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWGRAN2
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWHIST
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWAGE
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWIMPUR !N6K
REAL, DIMENSION(KSIZE1,KSIZE2,NIMPUR) :: ZP_SNOWIMPURV2 !N6K
REAL, DIMENSION(KSIZE1,KSIZE4) :: ZP_BLOWSNW
REAL, DIMENSION(KSIZE1)        :: ZP_SNOWALB
REAL, DIMENSION(KSIZE1)        :: ZP_SWNETSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_SWNETSNOWS
REAL, DIMENSION(KSIZE1)        :: ZP_LWNETSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_PS
REAL, DIMENSION(KSIZE1)        :: ZP_SRSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_RRSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_PSN3L
REAL, DIMENSION(KSIZE1)        :: ZP_TA
REAL, DIMENSION(KSIZE1)        :: ZP_CT
REAL, DIMENSION(KSIZE1,KSIZE3) :: ZP_TG
REAL, DIMENSION(KSIZE1,KSIZE3) :: ZP_D_G
REAL, DIMENSION(KSIZE1,KSIZE3) :: ZP_DZG
REAL, DIMENSION(KSIZE1,KSIZE3) :: ZP_SOILHCAPZ
REAL, DIMENSION(KSIZE1)        :: ZP_SOILD
REAL, DIMENSION(KSIZE1)        :: ZP_DELHEATG
REAL, DIMENSION(KSIZE1)        :: ZP_DELHEATG_SFC
REAL, DIMENSION(KSIZE1)        :: ZP_SW_RAD
REAL, DIMENSION(KSIZE1)        :: ZP_QA
REAL, DIMENSION(KSIZE1)        :: ZP_VMOD
REAL, DIMENSION(KSIZE1)        :: ZP_LW_RAD
REAL, DIMENSION(KSIZE1)        :: ZP_RHOA
REAL, DIMENSION(KSIZE1)        :: ZP_UREF
REAL, DIMENSION(KSIZE1)        :: ZP_EXNS
REAL, DIMENSION(KSIZE1)        :: ZP_EXNA
REAL, DIMENSION(KSIZE1)        :: ZP_DIRCOSZW
REAL, DIMENSION(KSIZE1)        :: ZP_ZREF
REAL, DIMENSION(KSIZE1)        :: ZP_Z0NAT
REAL, DIMENSION(KSIZE1)        :: ZP_Z0HNAT
REAL, DIMENSION(KSIZE1)        :: ZP_Z0EFF
REAL, DIMENSION(KSIZE1)        :: ZP_ALB
REAL, DIMENSION(KSIZE1)        :: ZP_SOILCOND
REAL, DIMENSION(KSIZE1)        :: ZP_THRUFAL
REAL, DIMENSION(KSIZE1)        :: ZP_GRNDFLUX
REAL, DIMENSION(KSIZE1)        :: ZP_FLSN_COR
REAL, DIMENSION(KSIZE1)        :: ZP_GSFCSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_EVAPCOR
REAL, DIMENSION(KSIZE1)        :: ZP_SOILCOR
REAL, DIMENSION(KSIZE1)        :: ZP_GFLXCOR
REAL, DIMENSION(KSIZE1)        :: ZP_RNSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_HSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_GFLUXSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_DELHEATN
REAL, DIMENSION(KSIZE1)        :: ZP_DELHEATN_SFC
REAL, DIMENSION(KSIZE1)        :: ZP_SNOWSFCH
REAL, DIMENSION(KSIZE1)        :: ZP_HPSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_LES3L
REAL, DIMENSION(KSIZE1)        :: ZP_LEL3L
REAL, DIMENSION(KSIZE1)        :: ZP_EVAP
REAL, DIMENSION(KSIZE1)        :: ZP_SNDRIFT
REAL, DIMENSION(KSIZE1)        :: ZP_RI
REAL, DIMENSION(KSIZE1)        :: ZP_QS
REAL, DIMENSION(KSIZE1)        :: ZP_EMISNOW
REAL, DIMENSION(KSIZE1)        :: ZP_CDSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_USTARSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_CHSNOW
REAL, DIMENSION(KSIZE1)        :: ZP_SNOWHMASS
REAL, DIMENSION(KSIZE1)        :: ZP_VEGTYPE
REAL, DIMENSION(KSIZE1)        :: ZP_PEW_A_COEF
REAL, DIMENSION(KSIZE1)        :: ZP_PEW_B_COEF
REAL, DIMENSION(KSIZE1)        :: ZP_PET_A_COEF
REAL, DIMENSION(KSIZE1)        :: ZP_PET_B_COEF
REAL, DIMENSION(KSIZE1)        :: ZP_PEQ_A_COEF
REAL, DIMENSION(KSIZE1)        :: ZP_PEQ_B_COEF
REAL, DIMENSION(KSIZE1)        :: ZP_ZENITH
REAL, DIMENSION(KSIZE1)        :: ZP_LAT,ZP_LON
REAL, DIMENSION(KSIZE1)        :: ZP_PSN_INV
REAL, DIMENSION(KSIZE1)        :: ZP_PSN
REAL, DIMENSION(KSIZE1)        :: ZP_PSN_GFLXCOR
REAL, DIMENSION(KSIZE1)        :: ZP_WORK

REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWDEND
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWSPHER
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWSIZE
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWSSA
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWTYPEMEPRA
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWRAM
REAL, DIMENSION(KSIZE1,KSIZE2) :: ZP_SNOWSHEAR
REAL, DIMENSION(KSIZE1) :: ZP_SNOWDEPTH_1DAYS
REAL, DIMENSION(KSIZE1) :: ZP_SNOWDEPTH_3DAYS
REAL, DIMENSION(KSIZE1) :: ZP_SNOWDEPTH_5DAYS
REAL, DIMENSION(KSIZE1) :: ZP_SNOWDEPTH_7DAYS
REAL, DIMENSION(KSIZE1) :: ZP_SNOWSWE_1DAYS
REAL, DIMENSION(KSIZE1) :: ZP_SNOWSWE_3DAYS
REAL, DIMENSION(KSIZE1) :: ZP_SNOWSWE_5DAYS
REAL, DIMENSION(KSIZE1) :: ZP_SNOWSWE_7DAYS
REAL, DIMENSION(KSIZE1) :: ZP_SNOWRAM_SONDE
REAL, DIMENSION(KSIZE1) :: ZP_SNOW_WETTHICKNESS
REAL, DIMENSION(KSIZE1) :: ZP_SNOW_REFROZENTHICKNESS

!
REAL, PARAMETER :: ZDEPTHABS = 0.60 ! m
!
INTEGER :: JWRK, JJ, JI
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SNOW3L_ISBA:CALL_MODEL',0,ZHOOK_HANDLE)

!
! Initialize:
!
ZP_PSN_GFLXCOR(:)  = 0.
ZP_WORK(:)         = 0.
ZP_SOILD(:)        = 0.
!
! pack the variables
!
DO JWRK=1,KSIZE2
   DO JJ=1,KSIZE1
      JI = KMASK(JJ)
      ZP_SNOWSWE (JJ,JWRK) = PSNOWSWE (JI,JWRK)
      ZP_SNOWRHO (JJ,JWRK) = PSNOWRHO (JI,JWRK)
      ZP_SNOWHEAT(JJ,JWRK) = PSNOWHEAT(JI,JWRK)
      ZP_SNOWTEMP(JJ,JWRK) = PSNOWTEMP(JI,JWRK)
      ZP_SNOWLIQ (JJ,JWRK) = PSNOWLIQ (JI,JWRK)
      ZP_SNOWDZ  (JJ,JWRK) = PSNOWDZ  (JI,JWRK)
      ZP_SNOWAGE (JJ,JWRK) = PSNOWAGE (JI,JWRK)

   ENDDO
ENDDO
!

IF (HSNOW_ISBA=='CRO') THEN
   DO JWRK=1,KSIZE2
      DO JJ=1,KSIZE1
         JI = KMASK(JJ)
         ZP_SNOWGRAN1(JJ,JWRK) = PSNOWGRAN1 (JI,JWRK)
         ZP_SNOWGRAN2(JJ,JWRK) = PSNOWGRAN2 (JI,JWRK)
         ZP_SNOWHIST (JJ,JWRK) = PSNOWHIST  (JI,JWRK)
         ZP_SNOWIMPUR(JJ,JWRK) =PSNOWIMPUR(JI,JWRK)  !N6K
      ENDDO
   ENDDO
   
   DO JIMP=1,NIMPUR
     DO JWRK=1,KSIZE2
        DO JJ=1,KSIZE1
           JI = KMASK(JJ)
           ZP_SNOWIMPURV2(JJ,JWRK,JIMP) =PSNOWIMPURV2(JI,JWRK,JIMP)  !N6K
        ENDDO
     ENDDO
   ENDDO
   
   DO JWRK=1,KSIZE4
      DO JJ=1,KSIZE1
        JI = KMASK(JJ)
        ZP_BLOWSNW(JJ,JWRK) = ZBLOWSNW(JI,JWRK)
      ENDDO
   ENDDO

ELSE
   DO JWRK=1,KSIZE2
      DO JJ=1,KSIZE1
         ZP_SNOWGRAN1(JJ,JWRK) = XUNDEF
         ZP_SNOWGRAN2(JJ,JWRK) = XUNDEF
         ZP_SNOWHIST (JJ,JWRK) = XUNDEF
         ZP_SNOWIMPUR(JJ,JWRK) = XUNDEF         !N6K
      ENDDO
   ENDDO
   
   
   DO JIMP=1,NIMPUR
     DO JWRK=1,KSIZE2
        DO JJ=1,KSIZE1
           ZP_SNOWIMPURV2(JJ,JWRK,JIMP) = XUNDEF         !N6K
        ENDDO
     ENDDO
   ENDDO

  DO JWRK=1,KSIZE4
     DO JJ=1,KSIZE1
        ZP_BLOWSNW(JJ,JWRK) = XUNDEF
     ENDDO
  ENDDO

ENDIF
!  
DO JWRK=1,KSIZE3
   DO JJ=1,KSIZE1
      JI                    = KMASK           (JJ)
      ZP_TG       (JJ,JWRK) = PTG        (JI,JWRK)
      ZP_D_G      (JJ,JWRK) = PD_G       (JI,JWRK)
      ZP_SOILHCAPZ(JJ,JWRK) = PSOILHCAPZ (JI,JWRK)
   ENDDO
ENDDO
!
IF (OMEB) THEN
  DO JWRK=1,KSIZE3
    DO JJ=1,KSIZE1
      JI                    = KMASK           (JJ)
      ZP_DZG      (JJ,JWRK) = PDZG       (JI,JWRK)
    ENDDO
  ENDDO
ENDIF
!
DO JJ=1,KSIZE1
   JI = KMASK(JJ)
   ZP_SNOWALB (JJ) = PSNOWALB (JI)    
   ZP_PS      (JJ) = PPS      (JI)
   ZP_SRSNOW  (JJ) = PSR      (JI)
   ZP_RRSNOW  (JJ) = ZRRSNOW  (JI)
   ZP_PSN3L   (JJ) = PPSN     (JI)
   ZP_CT      (JJ) = PCT      (JI)
   ZP_TA      (JJ) = PTA      (JI)
   ZP_DELHEATG(JJ) = PDELHEATG(JI)
   ZP_DELHEATG_SFC(JJ) = PDELHEATG_SFC(JI)
   ZP_SW_RAD  (JJ) = PSW_RAD  (JI)
   ZP_QA      (JJ) = PQA      (JI)
   ZP_VMOD    (JJ) = PVMOD    (JI)
   ZP_LW_RAD  (JJ) = PLW_RAD  (JI)
   ZP_RHOA    (JJ) = PRHOA    (JI)
   ZP_UREF    (JJ) = PUREF    (JI)
   ZP_EXNS    (JJ) = PEXNS    (JI)
   ZP_EXNA    (JJ) = PEXNA    (JI)
   ZP_DIRCOSZW(JJ) = PDIRCOSZW(JI)
   ZP_ZREF    (JJ) = PZREF    (JI)
   ZP_Z0NAT   (JJ) = PZ0NAT   (JI)
   ZP_Z0HNAT  (JJ) = PZ0HNAT  (JI)
   ZP_Z0EFF   (JJ) = PZ0EFF   (JI)
   ZP_ALB     (JJ) = PALB     (JI)
   ZP_SOILCOND(JJ) = ZSOILCOND(JI)
   !  
   ZP_PEW_A_COEF(JJ) = PPEW_A_COEF(JI)
   ZP_PEW_B_COEF(JJ) = PPEW_B_COEF(JI)
   ZP_PET_A_COEF(JJ) = PPET_A_COEF(JI)
   ZP_PEQ_A_COEF(JJ) = PPEQ_A_COEF(JI)      
   ZP_PET_B_COEF(JJ) = PPET_B_COEF(JI)
   ZP_PEQ_B_COEF(JJ) = PPEQ_B_COEF(JI)
   !
   ZP_LAT  (JJ)      = PLAT(JI)
   ZP_LON  (JJ)      = PLON(JI)
   ZP_ZENITH(JJ)     = PZENITH  (JI)
!
   ZP_GRNDFLUX    (JJ) = PGRNDFLUX    (JI)
   ZP_RNSNOW      (JJ) = PRNSNOW      (JI)
   ZP_HSNOW       (JJ) = PHSNOW       (JI)
   ZP_DELHEATN    (JJ) = PDELHEATN    (JI)
   ZP_DELHEATN_SFC(JJ) = PDELHEATN_SFC(JI)
   ZP_SNOWSFCH    (JJ) = PSNOWSFCH    (JI)
   ZP_HPSNOW      (JJ) = PHPSNOW      (JI)
   ZP_LES3L       (JJ) = PLES3L       (JI) 
   ZP_LEL3L       (JJ) = PLEL3L       (JI)  
   ZP_EVAP        (JJ) = PEVAP        (JI)
   ZP_EMISNOW     (JJ) = PEMISNOW     (JI) 
   ZP_SWNETSNOW   (JJ) = PSWNETSNOW   (JI) 
   ZP_SWNETSNOWS  (JJ) = PSWNETSNOWS  (JI) 
   ZP_LWNETSNOW   (JJ) = PLWNETSNOW   (JI) 
ENDDO
!
DO JJ=1,KSIZE1
   JI = KMASK(JJ)
   ZP_VEGTYPE (JJ) = PVEGTYPE (JI,NVT_SNOW)
ENDDO
!
!
! ===============================================================
! conversion of snow heat from J/m3 into J/m2
WHERE(ZP_SNOWSWE(:,:)>0.) &
  ZP_SNOWHEAT(:,:) = ZP_SNOWHEAT(:,:) / ZP_SNOWRHO (:,:) * ZP_SNOWSWE (:,:)  
! ===============================================================
!
ZP_PSN_INV(:)       = 0.0
ZP_PSN(:)           = ZP_PSN3L(:)
!
IF(OMEB)THEN
!
!   MEB (case of imposed surface fluxes)
!   - Prepare inputs for explicit snow scheme(s):
!     If using MEB, these are INPUTs ONLY:
!     divide fluxes by snow fraction to make "snow-relative"
!
   ZP_PSN(:)         = MAX(1.E-4, ZP_PSN3L(:))
   ZP_PSN_INV(:)     = 1.0/ZP_PSN(:)
!
   ZP_RNSNOW(:)      = ZP_RNSNOW(:)      *ZP_PSN_INV(:)
   ZP_SWNETSNOW(:)   = ZP_SWNETSNOW(:)   *ZP_PSN_INV(:)
   ZP_SWNETSNOWS(:)  = ZP_SWNETSNOWS(:)  *ZP_PSN_INV(:)
   ZP_LWNETSNOW(:)   = ZP_LWNETSNOW(:)   *ZP_PSN_INV(:)
   ZP_HSNOW(:)       = ZP_HSNOW(:)       *ZP_PSN_INV(:)
   ZP_LES3L(:)       = ZP_LES3L(:)       *ZP_PSN_INV(:)
   ZP_LEL3L(:)       = ZP_LEL3L(:)       *ZP_PSN_INV(:)
   ZP_GRNDFLUX(:)    = ZP_GRNDFLUX(:)    *ZP_PSN_INV(:)
   ZP_EVAP(:)        = ZP_EVAP(:)        *ZP_PSN_INV(:)
   ZP_HPSNOW(:)      = ZP_HPSNOW(:)      *ZP_PSN_INV(:)
   ZP_DELHEATN(:)    = ZP_DELHEATN(:)    *ZP_PSN_INV(:)
   ZP_DELHEATN_SFC(:)= ZP_DELHEATN_SFC(:)*ZP_PSN_INV(:)
   ZP_SNOWSFCH(:)    = ZP_SNOWSFCH(:)    *ZP_PSN_INV(:)

   ZP_SRSNOW(:)      = ZP_SRSNOW(:)      *ZP_PSN_INV(:)
   ZP_RRSNOW(:)      = ZP_RRSNOW(:)      *ZP_PSN_INV(:)

   DO JJ=1,KSIZE2
      DO JI=1,KSIZE1
         ZP_SNOWSWE(JI,JJ)  = ZP_SNOWSWE(JI,JJ) *ZP_PSN_INV(JI)
         ZP_SNOWHEAT(JI,JJ) = ZP_SNOWHEAT(JI,JJ)*ZP_PSN_INV(JI)
         ZP_SNOWDZ(JI,JJ)   = ZP_SNOWDZ(JI,JJ)  *ZP_PSN_INV(JI)
      ENDDO
   ENDDO
   !
ENDIF
!
! Call ISBA-SNOW3L model:  
!  
IF (HSNOW_ISBA=='CRO') THEN 

   CALL SNOWCRO(HSNOWRES, TPTIME, OMEB, OGLACIER, HIMPLICIT_WIND,          &
             ZP_PEW_A_COEF, ZP_PEW_B_COEF,                                 &
             ZP_PET_A_COEF, ZP_PEQ_A_COEF, ZP_PET_B_COEF, ZP_PEQ_B_COEF,   &
             ZP_SNOWSWE,ZP_SNOWRHO, ZP_SNOWHEAT, ZP_SNOWALB,               &
             ZP_SNOWGRAN1, ZP_SNOWGRAN2, ZP_SNOWHIST, ZP_SNOWAGE,          &
             ZP_SNOWIMPUR,PTSTEP,  &
             ZP_PS, ZP_SRSNOW, ZP_RRSNOW ,ZP_PSN3L, ZP_TA, ZP_TG(:,1),     &
             ZP_SW_RAD, ZP_QA, ZP_VMOD, ZP_LW_RAD, ZP_RHOA, ZP_UREF,       &
             ZP_EXNS, ZP_EXNA, ZP_DIRCOSZW, ZP_ZREF, ZP_Z0NAT, ZP_Z0EFF,   &
             ZP_Z0HNAT, ZP_ALB, ZP_SOILCOND, ZP_D_G(:,1), ZP_SNOWLIQ,      &
             ZP_SNOWTEMP, ZP_SNOWDZ, ZP_THRUFAL, ZP_GRNDFLUX, ZP_EVAPCOR,  &
             ZP_GFLXCOR, ZP_SWNETSNOW, ZP_SWNETSNOWS, ZP_LWNETSNOW,        &             
             ZP_RNSNOW, ZP_HSNOW, ZP_GFLUXSNOW, ZP_HPSNOW, ZP_LES3L,       &
             ZP_LEL3L, ZP_EVAP, ZP_SNDRIFT, ZP_RI,                         &
             ZP_EMISNOW, ZP_CDSNOW, ZP_USTARSNOW,                          &
             ZP_CHSNOW, ZP_SNOWHMASS, ZP_QS, ZP_VEGTYPE, ZP_ZENITH,        &
             ZP_LAT, ZP_LON, ZP_BLOWSNW, OSNOWDRIFT,OSNOWDRIFT_SUBLIM,     &
             OSNOW_ABS_ZENITH, HSNOWMETAMO,HSNOWRAD,OATMORAD,P_DIR_SW, P_SCA_SW,&
             PSPEC_ALB, PDIFF_RATIO)
!
  ZP_GFLXCOR (:) = 0.0
  ZP_FLSN_COR(:) = 0.0
  ZP_SOILCOR (:) = 0.0
!
   !SIZE(PSNOWDEND)>1 is equivalent to test the value of DGMI%LPROSNOW which does not enter in ISBA
#ifdef SFX_OL
  ! En couplé il faudra voir si on veut virer les diagnostics, les calculer tout le temps, ou trouver une autre solution
  GCOMPUTECRODIAG=(SIZE(PSNOWDEND)>1)
#else
  GCOMPUTECRODIAG=(SIZE(PSNOWDEND)>1).AND.(MOD(TPTIME%TIME,XTSTEP_OUTPUT) == 0.)
#endif

   !Ajout test sur pas de temps de sortie
   IF (GCOMPUTECRODIAG) THEN
     CALL SNOWCRO_DIAG(HSNOWMETAMO,&
                     ZP_SNOWDZ,ZP_SNOWSWE,ZP_SNOWRHO,ZP_SNOWGRAN1,ZP_SNOWGRAN2,ZP_SNOWAGE,ZP_SNOWHIST,ZP_SNOWTEMP,&
                     ZP_SNOWLIQ,ZP_DIRCOSZW,ZP_SNOWDEND,ZP_SNOWSPHER,ZP_SNOWSIZE,ZP_SNOWSSA,ZP_SNOWTYPEMEPRA,     &
                     ZP_SNOWRAM, ZP_SNOWSHEAR,ZP_SNOWDEPTH_1DAYS,ZP_SNOWDEPTH_3DAYS,ZP_SNOWDEPTH_5DAYS,           &
                     ZP_SNOWDEPTH_7DAYS,ZP_SNOWSWE_1DAYS,ZP_SNOWSWE_3DAYS,ZP_SNOWSWE_5DAYS,                       &
                     ZP_SNOWSWE_7DAYS,ZP_SNOWRAM_SONDE,ZP_SNOW_WETTHICKNESS,ZP_SNOW_REFROZENTHICKNESS)
   ENDIF


ELSE 
!
  CALL SNOW3L(HSNOWRES, TPTIME, OMEB, HIMPLICIT_WIND,                      &
             ZP_PEW_A_COEF, ZP_PEW_B_COEF,                                 &
             ZP_PET_A_COEF, ZP_PEQ_A_COEF,ZP_PET_B_COEF, ZP_PEQ_B_COEF,    &
             ZP_SNOWSWE, ZP_SNOWRHO, ZP_SNOWHEAT, ZP_SNOWALB,              &
             ZP_SNOWGRAN1, ZP_SNOWGRAN2, ZP_SNOWHIST, ZP_SNOWAGE, PTSTEP,  &
             ZP_PS, ZP_SRSNOW, ZP_RRSNOW, ZP_PSN3L, ZP_TA, ZP_TG(:,1),     &
             ZP_SW_RAD, ZP_QA, ZP_VMOD, ZP_LW_RAD, ZP_RHOA, ZP_UREF,       &
             ZP_EXNS, ZP_EXNA, ZP_DIRCOSZW, ZP_ZREF, ZP_Z0NAT, ZP_Z0EFF,   &
             ZP_Z0HNAT, ZP_ALB, ZP_SOILCOND, ZP_D_G(:,1), ZP_SNOWLIQ,      &
             ZP_SNOWTEMP, ZP_SNOWDZ, ZP_THRUFAL, ZP_GRNDFLUX ,             &
             ZP_EVAPCOR, ZP_SOILCOR, ZP_GFLXCOR, ZP_SNOWSFCH,              &
             ZP_DELHEATN, ZP_DELHEATN_SFC,                                 &
             ZP_SWNETSNOW, ZP_SWNETSNOWS, ZP_LWNETSNOW, ZP_GSFCSNOW,       &
             ZP_RNSNOW, ZP_HSNOW, ZP_GFLUXSNOW, ZP_HPSNOW, ZP_LES3L,       &
             ZP_LEL3L, ZP_EVAP, ZP_SNDRIFT, ZP_RI,                         &
             ZP_EMISNOW, ZP_CDSNOW, ZP_USTARSNOW,                          &
             ZP_CHSNOW, ZP_SNOWHMASS, ZP_QS, ZP_VEGTYPE, ZP_ZENITH,        &
             ZP_LAT, ZP_LON, OSNOWDRIFT, OSNOWDRIFT_SUBLIM                 )
!
  IF(OMEB)THEN
!
! - reverse transform: back to surface-relative
!
     ZP_RNSNOW(:)      = ZP_RNSNOW(:)      /ZP_PSN_INV(:)
     ZP_SWNETSNOW(:)   = ZP_SWNETSNOW(:)   /ZP_PSN_INV(:)
     ZP_SWNETSNOWS(:)  = ZP_SWNETSNOWS(:)  /ZP_PSN_INV(:)
     ZP_LWNETSNOW(:)   = ZP_LWNETSNOW(:)   /ZP_PSN_INV(:)
     ZP_HSNOW(:)       = ZP_HSNOW(:)       /ZP_PSN_INV(:)
     ZP_LES3L(:)       = ZP_LES3L(:)       /ZP_PSN_INV(:)
     ZP_LEL3L(:)       = ZP_LEL3L(:)       /ZP_PSN_INV(:)
     ZP_GRNDFLUX(:)    = ZP_GRNDFLUX(:)    /ZP_PSN_INV(:)
     ZP_EVAP(:)        = ZP_EVAP(:)        /ZP_PSN_INV(:)
     ZP_HPSNOW(:)      = ZP_HPSNOW(:)      /ZP_PSN_INV(:)
     ZP_GFLUXSNOW(:)   = ZP_GFLUXSNOW(:)   /ZP_PSN_INV(:) 
     ZP_DELHEATN(:)    = ZP_DELHEATN(:)    /ZP_PSN_INV(:) 
     ZP_DELHEATN_SFC(:)= ZP_DELHEATN_SFC(:)/ZP_PSN_INV(:) 
     ZP_SNOWSFCH(:)    = ZP_SNOWSFCH(:)    /ZP_PSN_INV(:) 
     ZP_GSFCSNOW(:)    = ZP_GSFCSNOW(:)    /ZP_PSN_INV(:) 

     ZP_SRSNOW(:)      = ZP_SRSNOW(:)      /ZP_PSN_INV(:)
     ZP_RRSNOW(:)      = ZP_RRSNOW(:)      /ZP_PSN_INV(:)
     DO JJ=1,KSIZE2
        DO JI=1,KSIZE1
           ZP_SNOWSWE(JI,JJ)  = ZP_SNOWSWE(JI,JJ) /ZP_PSN_INV(JI)
           ZP_SNOWHEAT(JI,JJ) = ZP_SNOWHEAT(JI,JJ)/ZP_PSN_INV(JI)
           ZP_SNOWDZ(JI,JJ)   = ZP_SNOWDZ(JI,JJ)  /ZP_PSN_INV(JI)
        ENDDO
     ENDDO
     
     ZP_SNOWHMASS(:)  = ZP_SNOWHMASS(:)/ZP_PSN_INV(:)
     ZP_THRUFAL(:)    = ZP_THRUFAL(:)  /ZP_PSN_INV(:)
!
!    Final Adjustments:
!    ------------------
!    Add cooling/heating flux correction to underlying soil.
!    This term is usually active for vanishingly thin snowpacks..
!    it is put outside of the snow scheme owing to it's dependence on
!    snow fraction. It is related to a possible correction to the ground-snow
!    heat flux when it is imposed (using MEB).
!    Also, it is added as a heat sink/source here since
!    fluxes have already be computed and should not be adjusted at this point:
!    applying it to the soil has the same impact as soil freeze-thaw, in the
!    sense it is computed after the fluxes have been updated.
!    (and update heat storage diagnostic in a consistent manner)
!
!    Energy is thickness weighted, thus thicker layers receive more energy and energy
!    is evenly distributed to depth ZDEPTHABS. An
!    alternate method is to weight near surface layers more and diminish weights
!    (thus eenrgy received by each layer) with depth. Both methods conserve energy as 
!    long as vertical weights are normalized.

!    i) Determine soil depth for energy absorption:

     ZP_SOILD(:) = ZP_DZG(:,1)
     DO JJ=2,KSIZE3
        DO JI=1,KSIZE1
           IF(ZP_DZG(JI,JJ) <= ZDEPTHABS)THEN
              ZP_SOILD(JI) = ZP_DZG(JI,JJ)
           ENDIF
        ENDDO
     ENDDO

!    ii) Distribute (possible) energy to absorb vertically over some layer (defined above):

     ZP_PSN_GFLXCOR(:)  = ZP_PSN(:)*ZP_GFLXCOR(:)                                ! (W/m2)
     ZP_WORK(:)         = ZP_PSN_GFLXCOR(:)*PTSTEP/ZP_SOILD(:)

     ZP_TG(:,1)         = ZP_TG(:,1)         + ZP_WORK(:)*ZP_CT(:)*ZP_D_G(:,1)   ! (K)
     DO JJ=2,KSIZE3
        DO JI=1,KSIZE1
           IF(ZP_SOILD(JI) <= ZDEPTHABS)THEN
              ZP_TG(JI,JJ) = ZP_TG(JI,JJ)    + ZP_WORK(JI)/ZP_SOILHCAPZ(JI,JJ)   ! K
           ENDIF
        ENDDO
     ENDDO

     ZP_DELHEATG(:)     = ZP_DELHEATG(:)     + ZP_PSN_GFLXCOR(:)                 ! (W/m2)
     ZP_DELHEATG_SFC(:) = ZP_DELHEATG_SFC(:) + ZP_PSN_GFLXCOR(:)                 ! (W/m2)
!
     ZP_FLSN_COR(:)     = 0.0
!
  ELSE
!
!    To conserve energy in ISBA, the correction flux must be distributed at least
!    over the first 60cm depth. This method prevent numerical oscillations
!    especially when explicit snow vanishes. Final Adjustments are done in ISBA_CEB
!
     ZP_FLSN_COR(:) = ZP_GFLXCOR(:) ! (W/m2)
!
  ENDIF
!
ENDIF
!
!===============================================================
!conversion of snow heat from J/m2 into J/m3
WHERE(ZP_SNOWSWE (:,:)>0.)
      ZP_SNOWHEAT(:,:)=ZP_SNOWHEAT(:,:)*ZP_SNOWRHO(:,:)/ZP_SNOWSWE(:,:)  
ENDWHERE
!===============================================================
!
! === Packing:
!
! unpack variables
!
DO JWRK=1,KSIZE2
  DO JJ=1,KSIZE1
    JI = KMASK(JJ)
    PSNOWSWE  (JI,JWRK) = ZP_SNOWSWE  (JJ,JWRK)
    PSNOWRHO  (JI,JWRK) = ZP_SNOWRHO  (JJ,JWRK)
    PSNOWHEAT (JI,JWRK) = ZP_SNOWHEAT (JJ,JWRK)
    PSNOWTEMP (JI,JWRK) = ZP_SNOWTEMP (JJ,JWRK)
    PSNOWLIQ  (JI,JWRK) = ZP_SNOWLIQ  (JJ,JWRK)
    PSNOWDZ   (JI,JWRK) = ZP_SNOWDZ   (JJ,JWRK)
    PSNOWAGE  (JI,JWRK) = ZP_SNOWAGE  (JJ,JWRK)
  ENDDO
ENDDO
!
IF (HSNOW_ISBA=='CRO') THEN
  DO JWRK=1,KSIZE2
    DO JJ=1,KSIZE1
      JI = KMASK(JJ)
      PSNOWGRAN1(JI,JWRK) = ZP_SNOWGRAN1(JJ,JWRK)
      PSNOWGRAN2(JI,JWRK) = ZP_SNOWGRAN2(JJ,JWRK)
      PSNOWHIST (JI,JWRK) = ZP_SNOWHIST (JJ,JWRK)
      PSNOWIMPUR(JI,JWRK) = ZP_SNOWIMPUR(JJ,JWRK)
    ENDDO
  ENDDO
  
  DO JIMP=1,NIMPUR
    DO JWRK=1,INLVLS
      DO JJ=1,SIZE(PSNOWSWE,1)
        JI = KMASK(JJ)
        PSNOWIMPURV2(JI,JWRK,JIMP) = ZP_SNOWIMPURV2(JJ,JWRK,JIMP) !N6K
      ENDDO
    ENDDO
  ENDDO      
  
  IF (SIZE(PSNOWDEND)>1) THEN
  ! This is equivalent to test the value of DGMI%LPROSNOW which does not enter in ISBA
    DO JWRK=1,KSIZE2
      DO JJ=1,KSIZE1
        JI = KMASK(JJ)
        PSNOWDEND(JI,JWRK) = ZP_SNOWDEND(JJ,JWRK)
        PSNOWSPHER(JI,JWRK) = ZP_SNOWSPHER(JJ,JWRK) 
        PSNOWSIZE(JI,JWRK) = ZP_SNOWSIZE(JJ,JWRK)
        PSNOWSSA(JI,JWRK) = ZP_SNOWSSA(JJ,JWRK)
        PSNOWTYPEMEPRA(JI,JWRK) = ZP_SNOWTYPEMEPRA(JJ,JWRK)
        PSNOWRAM(JI,JWRK) =  ZP_SNOWRAM(JJ,JWRK)
        PSNOWSHEAR(JI,JWRK) =  ZP_SNOWSHEAR(JJ,JWRK)
      ENDDO
    ENDDO
  ENDIF
ENDIF
!
DO JWRK=1,KSIZE3
   DO JJ=1,KSIZE1
      JI              = KMASK          (JJ)
      PTG    (JI,JWRK)= ZP_TG        (JJ,JWRK)
   ENDDO
ENDDO
!
DO JJ=1,KSIZE1
  JI                  = KMASK          (JJ)
  PDELHEATG    (JI)   = ZP_DELHEATG    (JJ)
  PDELHEATG_SFC(JI)   = ZP_DELHEATG_SFC(JJ)
  PSNOWALB     (JI)   = ZP_SNOWALB     (JJ)
  PTHRUFAL     (JI)   = ZP_THRUFAL     (JJ)
  PEVAPCOR     (JI)   = ZP_EVAPCOR     (JJ)
  ZSOILCOR     (JI)   = ZP_SOILCOR     (JJ)
  PRI          (JI)   = ZP_RI          (JJ)
  PQS          (JI)   = ZP_QS          (JJ)
  PCDSNOW      (JI)   = ZP_CDSNOW      (JJ)
  PUSTARSNOW   (JI)   = ZP_USTARSNOW   (JJ)
  PCHSNOW      (JI)   = ZP_CHSNOW      (JJ)
  PSNOWHMASS   (JI)   = ZP_SNOWHMASS   (JJ)
  PGRNDFLUX    (JI)   = ZP_GRNDFLUX    (JJ)
  PFLSN_COR    (JI)   = ZP_FLSN_COR    (JJ)
  PRNSNOW      (JI)   = ZP_RNSNOW      (JJ)
  PHSNOW       (JI)   = ZP_HSNOW       (JJ)
  PGFLUXSNOW   (JI)   = ZP_GFLUXSNOW   (JJ)
  PDELHEATN    (JI)   = ZP_DELHEATN    (JJ)
  PDELHEATN_SFC(JI)   = ZP_DELHEATN_SFC(JJ)
  PSNOWSFCH    (JI)   = ZP_SNOWSFCH    (JJ)
  PGSFCSNOW    (JI)   = ZP_GSFCSNOW    (JJ)
  PHPSNOW      (JI)   = ZP_HPSNOW      (JJ)
  PLES3L       (JI)   = ZP_LES3L       (JJ)
  PLEL3L       (JI)   = ZP_LEL3L       (JJ)
  PEVAP        (JI)   = ZP_EVAP        (JJ)
  PEMISNOW     (JI)   = ZP_EMISNOW     (JJ)
  PSWNETSNOW   (JI)   = ZP_SWNETSNOW   (JJ)
  PSWNETSNOWS  (JI)   = ZP_SWNETSNOWS  (JJ)
  PLWNETSNOW   (JI)   = ZP_LWNETSNOW   (JJ)
ENDDO

IF ((HSNOW_ISBA=='CRO')  .AND. (SIZE(PSNOWDEND)>1)) THEN
      ! This is equivalent to test the value of DGMI%LPROSNOW which does not enter in ISBATHEN
  PSNOWDEPTH_1DAYS(:)=XUNDEF
  PSNOWDEPTH_3DAYS(:)=XUNDEF
  PSNOWDEPTH_5DAYS(:)=XUNDEF
  PSNOWDEPTH_7DAYS(:)=XUNDEF  
  PSNOWSWE_1DAYS(:)=XUNDEF 
  PSNOWSWE_3DAYS(:)=XUNDEF 
  PSNOWSWE_5DAYS(:)=XUNDEF 
  PSNOWSWE_7DAYS(:)=XUNDEF 
  PSNOWRAM_SONDE(:)=XUNDEF
  PSNOW_WETTHICKNESS(:)=XUNDEF
  PSNOW_REFROZENTHICKNESS(:)=XUNDEF   
  DO JJ=1,KSIZE1
    JI = KMASK(JJ)
    PSNOWDEPTH_1DAYS(JI)=ZP_SNOWDEPTH_1DAYS(JJ)
    PSNOWDEPTH_3DAYS(JI)=ZP_SNOWDEPTH_3DAYS(JJ)
    PSNOWDEPTH_5DAYS(JI)=ZP_SNOWDEPTH_5DAYS(JJ)
    PSNOWDEPTH_7DAYS(JI)=ZP_SNOWDEPTH_7DAYS(JJ)
    PSNOWSWE_1DAYS(JI)=ZP_SNOWSWE_1DAYS(JJ)
    PSNOWSWE_3DAYS(JI)=ZP_SNOWSWE_3DAYS(JJ)
    PSNOWSWE_5DAYS(JI)=ZP_SNOWSWE_5DAYS(JJ)
    PSNOWSWE_7DAYS(JI)=ZP_SNOWSWE_7DAYS(JJ)     
    PSNOWRAM_SONDE(JI)=ZP_SNOWRAM_SONDE(JJ)
    PSNOW_WETTHICKNESS(JI)=ZP_SNOW_WETTHICKNESS(JJ)
    PSNOW_REFROZENTHICKNESS(JI)=ZP_SNOW_REFROZENTHICKNESS(JJ)    
  ENDDO
ENDIF
!
IF (LHOOK) CALL DR_HOOK('SNOW3L_ISBA:CALL_MODEL',1,ZHOOK_HANDLE)
!
END SUBROUTINE CALL_MODEL
!
END SUBROUTINE SNOW3L_ISBA
