!     ######spl
      SUBROUTINE ISBA_FLUXES(HISBA, HSNOW_ISBA, HSOILFRZ, HDIF, OTEMP_ARP,     &
                          PTSTEP, PSODELX,                                     &
                          PCD, PVMOD, PSW_RAD, PLW_RAD, PTA, PQA,              &
                          PVMODP, PRHOA, PEXNS, PEXNA, PCPS, PLVTT, PLSTT,     &
                          PLAI, PVEG, PHUG, PHUI, PHV,                         &
                          PLEG_DELTA, PLEGI_DELTA, PDELTA, PRA,                &
                          PF5, PRS, PCS, PCG, PCT, PSNOWSWE, PT2M, PTSM,       &
                          PPSN, PPSNV, PPSNG, PFROZEN1,                        &
                          PTAUICE,                                             &
                          PALBT, PEMIST, PQSAT, PDQSAT, PSNOW_THRUFAL,         &
                          PRN, PH, PLE, PLEG, PLEGI, PLEV,                     &
                          PLES, PLER, PLETR, PEVAP,                            &
                          PGFLUX, PMELTADV, PMELT, PRESTORE,                   &
                          PUSTAR, PTS_RAD,                                     &
                          PDWGI1, PDWGI2, PDELTAT,                             &
                          PSOILCONDZ, PSOILHCAPZ, PWSATZ, PMPOTSATZ,           &
                          PBCOEFZ, PD_G, PTG, PWGI, PWG, PSR, PPSNV_A,         &
                          PFFG, PFFV, PFF, PFFROZEN, PLE_FLOOD, PLEI_FLOOD,    &
                          PSNOWTEMP, PALPHA, PN, PM, PWRES, PTDIURN            )
!     ##########################################################################
!
!!****  *ISBA_FLUXES*  
!!
!!    PURPOSE
!!    -------
!
!     Calculates the 
!     1.) evolution of the surface and deep-soil temperature
!     (i.e., Ts and T2) due to snow melt (DEF option) and soil ice phase changes, 
!     2.) the surface fluxes.
!         
!     
!!**  METHOD
!!    ------
!
!     1- snow melt latent heat, liquid rate (DEF option)
!     2- latent heating from soil ice phase changes
!     3- derive the surface fluxes.
!
!!    EXTERNAL
!!    --------
!!
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!    Noilhan and Planton (1989)
!!    Belair (1995)
!!    Douville et al. (1995)
!!    Boone et al. (2000)
!!      
!!    AUTHOR
!!    ------
!!
!!	S. Belair           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    14/03/95 
!!      (J.Stein)   15/11/95 use the wind components in the flux computation
!!      (J.Noilhan) 15/03/96 use the potential temperature instead of the
!!                           temperature for the heat flux computation 
!!      (J.Stein)   27/03/96 use only H and LE in the soil scheme
!!      (A.Boone V.Masson) 05/10/98 splits e_budget in two for CO2
!!      (A.Boone)   03/10/99 explicit latent heat of sublimation calculated 
!!      (A.Boone)   08/11/00 snow melt and soil ice phase changes herein
!!      (A.Boone)   06/05/02 Updates, ordering. Addition of 'HSOILFRZ' option
!!                           and introduction of snow melt timescale to 'DEF' snow option
!!      (P.LeMoigne) 01/07/05 expression of latent heat flux as a function of
!!                            w'theta' instead of w'T' (divison by surface exner)
!!      (P.LeMoigne) 28/07/05 dependence on qs for cp
!!      (A. Dziedzic and PLM) 10/2006 EBA snow option
!!      (B. Decharme)01/2009  Floodplains
!!      (B. Decharme)03/2009  BUG : effect of insolation due to vegetation cover
!!                                  at 1 for bare soil
!!      (R. Hamdi)   01/09    Cp and L are not constants (As in ALADIN)
!!      (B. Decharme)09/2009  Close the energy budget with the D95 snow scheme
!!      (A.Boone)    03/2010  Add delta fnctions to force LEG ans LEGI=0
!!                            when hug(i)Qsat < Qa and Qsat > Qa
!!      (A.Boone)    11/2011  Add RS_max limit to Etv
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,       ONLY : XSTEFAN, XCPD, XLSTT, XLVTT, XCL, XTT, XPI, XDAY, &
                            XCI, XRHOLI, XLMTT, XRHOLW, XG, XCL, XCONDI
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_ISBA_PAR,   ONLY : XWGMIN, XSPHSOIL, XDRYWGHT, XRS_MAX
!
USE MODE_THERMOS
!
USE MODI_ICE_SOILDIF              
!
USE MODE_SURF_SNOW_FRAC
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=*),   INTENT(IN)   :: HISBA   ! type of soil (Force-Restore OR Diffusion)
!                                           ! '2-L'
!                                           ! '3-L'
!                                           ! 'DIF'   ISBA-DF
!
CHARACTER(LEN=*), INTENT(IN)  :: HSNOW_ISBA ! 'DEF' = Default F-R snow scheme
!                                           !         (Douville et al. 1995)
!                                           ! '3-L' = 3-L snow scheme (option)
!                                           !         (Boone and Etchevers 2001)
!
CHARACTER(LEN=*),   INTENT(IN) :: HSOILFRZ  ! soil freezing-physics option
!                                           ! 'DEF'   Default (Boone et al. 2000; Giard and Bazile 2000)
!                                           ! 'LWT'   phase changes as above, but relation between unfrozen 
!                                                     water and temperature considered
!
CHARACTER(LEN=*),  INTENT(IN)  :: HDIF      ! NOTE: Only used when HISBA = DIF
!                                           ! 'BC' = Brook and Corey
!                                           ! 'VG' = Van Genuchten
!
LOGICAL, INTENT(IN)              :: OTEMP_ARP ! True  = time-varying force-restore soil temperature (as in ARPEGE)
                                               ! False = No time-varying force-restore soil temperature (Default)
!
!
REAL, INTENT (IN)                :: PTSTEP ! model time step (s)
!
!
REAL, DIMENSION(:), INTENT(IN)   :: PSODELX  ! Pulsation for each layer (Only used if LTEMP_ARP=True)
!
                                       
REAL, DIMENSION(:), INTENT (IN)  :: PSW_RAD, PLW_RAD, PTA, PQA, PVMODP, PRHOA
!                                     PSW_RAD = incoming solar radiation
!                                     PLW_RAD = atmospheric infrared radiation
!                                     PTA = near-ground air temperature
!                                     PQA = near-ground air specific humidity
!                                     PVMODP = near-ground wind speed
!                                     PRHOA = near-ground air density
!
REAL, DIMENSION(:), INTENT(IN)   :: PEXNS, PEXNA
REAL, DIMENSION(:), INTENT(IN)   :: PVEG, PLAI
REAL, DIMENSION(:), INTENT(IN)   :: PHUG, PHUI, PHV, PDELTA, PRA, PRS, PF5
REAL, DIMENSION(:), INTENT(IN)   :: PPSN, PPSNV, PPSNG, PFROZEN1, PTAUICE
REAL, DIMENSION(:), INTENT(IN)   :: PALBT, PEMIST
REAL, DIMENSION(:), INTENT(IN)   :: PQSAT, PDQSAT
REAL, DIMENSION(:), INTENT(IN)   :: PLEG_DELTA, PLEGI_DELTA
!                                     PVEG = fraction of vegetation
!                                     PHUG = relative humidity of the soil
!                                     PHV = Halstead coefficient
!                                     PF5 = water stress numerical correction factor (based on F2)
!                                     PDELTA = fraction of the foliage covered
!                                              by intercepted water
!                                     PRA = aerodynamic surface resistance for
!                                           heat transfers
!                                     PRS = surface stomatal resistance
!                                     PPSN = grid fraction covered by snow
!                                     PPSNV = fraction of the vegetation covered
!                                             by snow
!                                     PPSNG = fraction of the ground covered by
!                                             snow 
!                                     PFROZEN1 = fraction of ice in near-surface
!                                                ground
!                                     PALBT = area averaged albedo
!                                     PEMIST = area averaged emissivity
!                                     PQSAT = stauration vapor humidity at 't'
!                                     PDQSAT= stauration vapor humidity derivative at 't'
!                                     PLAI  = Leaf Area Index (m2 m-2)
!                                     PTAUICE = characteristic time scale for soil water
!                                               phase changes (s)
!                                     PLEG_DELTA = soil evaporation delta fn
!                                     PLEGI_DELTA = soil evaporation delta fn
!
REAL, DIMENSION(:), INTENT (IN)  :: PCS, PCG, PCT, PT2M, PTSM, PSNOWSWE
!                                     PCT    = area-averaged heat capacity (K m2 J-1)
!                                     PCS    = heat capacity of the snow (K m2 J-1)
!                                     PCG    = heat capacity of the soil (K m2 J-1)
!                                     PT2M   = mean surface (or restore) temperature at start 
!                                              of time step (K)
!                                     PTSM   = surface temperature at start 
!                                              of time step (K)
!                                     PSNOWSWE = equivalent water content of
!                                              the snow reservoir (kg m-2)
!
REAL, DIMENSION(:), INTENT(IN)   :: PCD, PVMOD
!                                     PCD = drag coefficient for momentum
!                                     PVMOD = module of the surface tangential
!                                             wind
!
REAL, DIMENSION(:), INTENT(IN)   :: PSNOW_THRUFAL
!                                     PSNOW_THRUFAL = rate that liquid water leaves snow pack: 
!                                                ISBA-ES [kg/(m2 s)]
REAL, DIMENSION(:), INTENT(IN)   :: PSR 
!                                     PSR = snow precipitation rate [kg/(m2 s)]
REAL, DIMENSION(:), INTENT(IN)   :: PPSNV_A
!                                     PPSNV_A = vegetation covered by snow EBA scheme
!
REAL, DIMENSION(:), INTENT(OUT)  :: PRN, PH, PLE, PLEG, PLEV, PLES
REAL, DIMENSION(:), INTENT(OUT)  :: PLER, PLETR, PEVAP, PGFLUX, PMELTADV, PMELT, PRESTORE
!                                     PRN = net radiation at the surface
!                                     PH = sensible heat flux
!                                     PLE = latent heat flux
!                                     PLEG = latent heat flux from the soil surface
!                                     PLEV = latent heat flux from the vegetation
!                                     PLES = latent heat flux from the snow
!                                     PLER = direct evaporation from the fraction
!                                            delta of the foliage
!                                     PLETR = transpiration of the remaining
!                                             part of the leaves
!                                     PEVAP = total evaporative flux (kg/m2/s)
!                                     PGFLUX = ground flux
!                                     PMELTADV = heat advection by melting snow
!                                                (acts to restore temperature to
!                                                 melting point) (W/m2)
!                                     PMELT = melting rate of snow (kg m-2 s-1)
!                                     PRESTORE = surface restore flux (W m-2)
!
REAL, DIMENSION(:), INTENT(OUT)  :: PLEGI
!                                     PLEGI = sublimation component of the 
!                                             latent heat flux from the soil surface
!
REAL, DIMENSION(:), INTENT(OUT)  :: PUSTAR
!                                     friction velocity
REAL, DIMENSION(:), INTENT(OUT)  :: PTS_RAD    ! effective radiative temperature 
!                                                of the natural surface (K)
!
REAL, DIMENSION(:), INTENT(OUT) ::    PDWGI1, PDWGI2
!                                     PDWGI1   = near-surface liquid water equivalent
!                                                volumetric ice content tendency
!                                     PDWGI2   = deep ground liquid water equivalent
!                                                volumetric ice content tendency
!
REAL, DIMENSION(:,:), INTENT(IN) ::   PDELTAT
!                                      PDELTAT = change in temperature over the time
!                                                step before adjustment owing to phase 
!                                                changes (K)
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PD_G, PWSATZ, PSOILHCAPZ, PSOILCONDZ
!                                      PD_G   = Depth of bottom of Soil layers (m)
!                                      PWSATZ    = porosity (m3/m3)
!                                      PSOILHCAPZ=ISBA-DF Soil heat capacity profile [J/(m3 K)]
!                                      PSOILCONDZ= ISBA-DF Soil conductivity profile  [W/(m K)]
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PMPOTSATZ, PBCOEFZ
!                                      PMPOTSATZ = matric potential at saturation (m)
!                                      PBCOEFZ   = slope of the water retention curve (-)
!
REAL, DIMENSION(:), INTENT(IN)      :: PCPS, PLVTT, PLSTT
!                                      PCPS  = heat capacity at surface
!
REAL, DIMENSION(:,:), INTENT(INOUT):: PWG, PWGI, PTG 
!                                     PWGI   = soil frozen volumetric water content (m3/m3)
!                                     PWG    = soil liquid volumetric water content (m3/m3)
!                                     PTG    = soil temperature profile (K)
!
REAL, DIMENSION(:), INTENT(IN)   :: PFFV      !Floodplain fraction over vegetation
REAL, DIMENSION(:), INTENT(IN)   :: PFF       !Floodplain fraction at the surface
REAL, DIMENSION(:), INTENT(IN)   :: PFFG   !Efective floodplain fraction
REAL, DIMENSION(:), INTENT(IN)   :: PFFROZEN  !fraction of frozen flood
REAL, DIMENSION(:), INTENT(OUT)  :: PLE_FLOOD, PLEI_FLOOD !Floodplains latent heat flux           [W/m²]
REAL, DIMENSION(:), INTENT(OUT)  :: PSNOWTEMP  ! snow layer temperatures (K)
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PALPHA      ! Van Genuchten parameter
REAL, DIMENSION(:,:), INTENT(IN)    :: PN          ! Van Genuchten parameter
REAL, DIMENSION(:,:), INTENT(IN)    :: PM          ! Van Genuchten parameter
REAL, DIMENSION(:,:), INTENT(IN)    :: PWRES       ! Van Genuchten parameter
!
REAL, DIMENSION(:), INTENT(OUT)     :: PTDIURN
!                                      PTDIURN      = penetration depth for restore (m)
!
!*      0.2    declarations of local variables
!
REAL                        ::   ZKSOIL     ! coefficient for soil freeze/thaw
!
REAL, DIMENSION(SIZE(PLAI)) :: ZZHV,                  &
!                                             for the calculation of the latent
!                                             heat of evapotranspiration
!
                               ZTN
!                                             average temperature used in 
!                                             the calculation of the 
!                                             melting effect
!
REAL, DIMENSION(SIZE(PLAI)) ::   ZKSFC_FRZ, ZFREEZING, ZICE_MELT, ZWIM,       &
                                 ZWIT, ZWGI1, ZWGI2, ZWM, ZSOILHEATCAP,       &
                                 ZICEEFF, ZEFFIC, ZDT, ZKSFC_IVEG,            &
                                 ZWGMIN, ZTGMAX, ZMATPOT, ZDELTAT, ZDELTATI
!                                ZKSFC_FRZ    = surface insolation coefficient (kg m-2 K-1)
!                                ZKSFC_IVEG   = non-dimensional vegetation insolation coefficient
!                                ZFREEZING    = rate for freezing soil water (kg m-2)
!                                ZICE_MELT    = rate for melting soil ice (kg m-2)
!                                ZWIM,ZWIT    = available ice content (m3 m-3)
!                                ZWGI1,ZWGI2  = volumetric ice contents (m3 m-3)
!                                ZWM          = available liquid water content (m3 m-3)
!                                ZSOILHEATCAP = soil heat capacity (J  m-3 K-1)
!                                ZICEEFF      = effective soil ice penetration depth (m)
!                                ZEFFIC       = phase change efficiency
!                                ZDT          = initial temperature change over time step (K)
!                                ZMATPOT      = soil matric potential (m)
!                                ZWGMIN       = volumetric water content above which soil water can
!                                               be unfrozen (if energy and mass available)(m3 m-3)
!                                ZTGMAX       = temperature below which liquid water 
!                                               can be frozen (if energy and mass available)(K)
!                                ZDELTA       = Freezing or melting temperature depression (K) after 
!                                               possible flux correction
!                                ZDELTAI      = Initial Freezing or melting temperature depression (K)
!
REAL, DIMENSION(SIZE(PLAI)) ::  ZPSN, ZPSNV, ZPSNG, ZFRAC
!                               ZPSN, ZPSNV, ZPSNG = snow fractions corresponding to
!                                                    dummy arguments PPSN, PPSNG, PPSNV
!                                                    if HSNOW_ISBA = 'DEF' (composite
!                                                    or Force-Restore snow scheme), else
!                                                    they are zero for explicit snow case
!                                                    as snow fluxes calculated outside of
!                                                    this routine using the 
!                                                    HSNOW_ISBA = '3-L' option.
!
REAL, DIMENSION(SIZE(PLAI)) ::  ZNEXTSNOW
!                               ZNEXTSNOW = Futur snow reservoir to close the
!                                           energy budget (see hydro_snow.f90)
!
REAL, DIMENSION(SIZE(PLAI)) ::  ZCONDAVG, ZWSAT_AVGZ
!                               ZCONDAVG   = average thermal conductivity of surface
!                                            and sub-surface layers (W m-1 K-1)
!                               ZWSAT_AVGZ = soil column average porosity (m3 m-3)
!
!
!*      0.2    local arrays for EBA scheme
!
REAL                            :: ZEPS1
!
!*      0.3    declarations of local parameters
!
REAL, PARAMETER             :: ZINSOLFRZ_VEG = 0.20  ! (-)       Vegetation insolation coefficient
!
REAL, PARAMETER             :: ZINSOLFRZ_LAI = 30.0  ! (m2 m-2)  Vegetation insolation coefficient
!
REAL, PARAMETER             :: ZTAU_SNOWMELT = 300.  ! (s)       Snow Melt timescale: needed to
!                                                                prevent time step dependence of melt
!                                                                when snow fraction < unity. 
!
REAL, PARAMETER             :: ZTWGHT     = 0.50 ! (-)   (0 < ZTWGHT <= 1/2)
!                                                                Weight for averaging the actual and flux corrected
!                                                                temperature depressions. Default is 1/2
!                                                                If ZTWGHT=0, then flux correction is OFF.
REAL, PARAMETER             :: ZEFFIC_MIN = 0.01 ! (-)   (0 <= ZEFFIC_MIN << 1)
!                                                                This parameter ensures
!                                                                a small minimum melt or freeze efficiency...
!                                                                It is numerical. When it is small, it has
!                                                                a only small impact on results, except
!                                                                that it keeps very small values of ice from persisting
!                                                                over long periods of time as they approach zero.
!                                                                If it is zero, then this effect off.
INTEGER         :: JJ
REAL, DIMENSION(SIZE(PLAI))          :: ZWORK1, ZWORK2, ZWORK3
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!
!*       0.     Initialization
!               --------------
IF (LHOOK) CALL DR_HOOK('ISBA_FLUXES',0,ZHOOK_HANDLE)
IF (HSNOW_ISBA == 'EBA') ZEPS1=1.0E-8
!
PMELT(:)        = 0.0
PLER(:)         = 0.0 
!
ZTN(:)          = 0.0
!
ZFREEZING(:)    = 0.0
ZKSFC_FRZ(:)    = 0.0
ZKSFC_IVEG(:)   = 0.0
ZEFFIC(:)       = 0.0
ZICE_MELT(:)    = 0.0
ZWGI1(:)        = 0.0
ZWIM(:)         = 0.0
ZSOILHEATCAP(:) = 0.0
PTDIURN(:)      = 0.0
ZWIT(:)         = 0.0
ZWGI2(:)        = 0.0
ZCONDAVG(:)     = 0.0
ZTGMAX(:)       = 0.0
ZWGMIN(:)       = 0.0
ZMATPOT(:)      = 0.0
ZWSAT_AVGZ(:)   = XUNDEF
!
! If ISBA-ES option in use, then snow covered surface
! fluxes calculated outside of this routine, so set
! the local snow fractions here to zero:
! 
IF(HSNOW_ISBA == '3-L' .OR. HSNOW_ISBA == 'CRO' .OR. HISBA == 'DIF')THEN
   ZPSN(:)      = 0.0
   ZPSNG(:)     = 0.0
   ZPSNV(:)     = 0.0
ELSE
   ZPSN(:)      = PPSN(:)
   ZPSNG(:)     = PPSNG(:)+PFFG(:)
   ZPSNV(:)     = PPSNV(:)+PFFV(:)
   ZFRAC(:)     = PPSNG(:)
ENDIF
!
!
DO JJ=1,SIZE(PTG,1)
!-------------------------------------------------------------------------------
!*       1.     FLUX CALCULATIONS
!               -----------------
!                                            temperature change
  ZDT(JJ)      = PTG(JJ,1) - PTSM(JJ)
!
!                                            net radiation
!
  PRN(JJ)      = (1. - PALBT(JJ)) * PSW_RAD(JJ) + PEMIST(JJ) *      &
              (PLW_RAD(JJ) - XSTEFAN * (PTSM(JJ)** 3)*(4.*PTG(JJ,1) - 3.*PTSM(JJ)))
!
  PTS_RAD(JJ)=((PTSM(JJ)** 3)*(4.*PTG(JJ,1) - 3.*PTSM(JJ)))**0.25
!                                            sensible heat flux
!
  PH(JJ)       = PRHOA(JJ) * PCPS(JJ) * (PTG(JJ,1) - PTA(JJ)*PEXNS(JJ)/PEXNA(JJ)) &
              / PRA(JJ) / PEXNS(JJ)
!
  ZWORK1(JJ) = PRHOA(JJ) * (1.-PVEG(JJ))*(1.-ZPSNG(JJ)) / PRA(JJ)
  ZWORK2(JJ) = PQSAT(JJ)+PDQSAT(JJ)*ZDT(JJ) 
!                                            latent heat of sublimation from
!                                            the ground
!
  PLEGI(JJ)    = ZWORK1(JJ) * PLSTT(JJ) * ( PHUI(JJ) * ZWORK2(JJ) - PQA(JJ)) * PFROZEN1(JJ) * PLEGI_DELTA(JJ)
!
!                                            total latent heat of evaporation from
!                                            the ground
!
  PLEG(JJ)     = ZWORK1(JJ) * PLVTT(JJ) * ( PHUG(JJ) * ZWORK2(JJ) - PQA(JJ)) * (1.-PFROZEN1(JJ)) * PLEG_DELTA(JJ)
!
  ZWORK2(JJ) = PRHOA(JJ) * (ZWORK2(JJ) - PQA(JJ))
  ZWORK3(JJ) = ZWORK2(JJ) / PRA(JJ)
!                                            latent heat of evaporation from 
!                                            the snow canopy
!
  PLES(JJ)     = PLSTT(JJ) * ZPSN(JJ) * ZWORK3(JJ)
!
!                                            latent heat of evaporation from
!                                            evaporation
!
  PLEV(JJ)     = PLVTT(JJ) * PVEG(JJ)*(1.-ZPSNV(JJ)) * PHV(JJ) * ZWORK3(JJ)
!
!                                            latent heat of evapotranspiration
!                                            
  ZZHV(JJ)     = MAX(0., SIGN(1.,PQSAT(JJ) - PQA(JJ)))
  PLETR(JJ)    = ZZHV(JJ) * (1. - PDELTA(JJ)) * PLVTT(JJ) * PVEG(JJ)*(1-ZPSNV(JJ))          &
               * ZWORK2(JJ) *( (1/(PRA(JJ) + PRS(JJ))) - ((1.-PF5(JJ))/(PRA(JJ) + XRS_MAX)) )

!
  PLER(JJ)     = PLEV(JJ) - PLETR(JJ)
!
!                                            latent heat of free water (floodplains)
!
  PLE_FLOOD(JJ)  = PLVTT(JJ) * (1.-PFFROZEN(JJ)) * PFF(JJ) * ZWORK3(JJ) 
!
  PLEI_FLOOD(JJ) = PLSTT(JJ) * PFFROZEN(JJ) * PFF(JJ) * ZWORK3(JJ) 
!
!                                            total latent heat of evaporation
!                                            without flood
!
  PLE(JJ)      = PLEG(JJ) + PLEV(JJ) + PLES(JJ) + PLEGI(JJ)
!
!                                            heat flux into the ground
!                                            without flood
!
  PGFLUX(JJ)   = PRN(JJ) - PH(JJ) - PLE(JJ)
!
!                                            heat flux due to snow melt
!                                            (ISBA-ES/SNOW3L)
!
  PMELTADV(JJ) = PSNOW_THRUFAL(JJ)*XCL*(XTT - PTG(JJ,1))
!
!                                            restore heat flux in FR mode,
!                                            or surface to sub-surface heat
!                                            flux using the DIF mode.
!
!
  PEVAP(JJ)    = ((PLEV(JJ) + PLEG(JJ))/PLVTT(JJ)) + ((PLEGI(JJ) + PLES(JJ))/PLSTT(JJ))
!                                            total evaporative flux (kg/m2/s)
!                                            without flood
!-------------------------------------------------------------------------------
!
!*       2.     FRICTION VELOCITY
!               -----------------
!
!
  PUSTAR(JJ) = SQRT( PCD(JJ) * PVMOD(JJ) * PVMODP(JJ) )
!
ENDDO

IF(HSNOW_ISBA == 'D95')THEN
!cdir nodep
  DO JJ=1,SIZE(PTG,1)
    PLE    (JJ)  = PLE    (JJ) + PLE_FLOOD(JJ) + PLEI_FLOOD(JJ)
    PGFLUX (JJ)  = PGFLUX (JJ) - PLE_FLOOD(JJ) - PLEI_FLOOD(JJ)
    PEVAP  (JJ)  = PEVAP  (JJ) + PLE_FLOOD(JJ)/PLVTT(JJ) + PLEI_FLOOD(JJ)/PLSTT(JJ)
  ENDDO
ENDIF
!
IF (HISBA /= 'DIF') THEN
!
! "Restore" flux between surface and deep layer:
!
   PRESTORE(:) = 2.0*XPI*(PTG(:,1)-PT2M(:))/(PCT(:)*XDAY)  
!
   IF(OTEMP_ARP)PRESTORE(:)=PRESTORE(:)/(PSODELX(1)*(PSODELX(1)+PSODELX(2)))
!
ELSE
!
! "Restore" flux here is actually the heat flux between the surface
! and sub-surface layers:
!
  DO JJ=1,SIZE(PTG,1)
    ZCONDAVG(JJ) = (PD_G(JJ,1)*PSOILCONDZ(JJ,1) + (PD_G(JJ,2)-PD_G(JJ,1))*PSOILCONDZ(JJ,2))/     &
               PD_G(JJ,2)
!
    PRESTORE(JJ) = 2.*ZCONDAVG(JJ)*(PTG(JJ,1)-PT2M(JJ))/PD_G(JJ,2)
  ENDDO
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       3.     SNOWMELT LATENT HEATING EFFECTS ('DEF' option)
!               ----------------------------------------------
!
IF( (HSNOW_ISBA == 'D95' .OR. HSNOW_ISBA == 'EBA') .AND. HISBA /= 'DIF' )THEN
!                                            temperature tn
!
    IF (HSNOW_ISBA == 'D95') THEN
!           
      ZTN       (:) = (1.-PVEG(:))*PTG(:,1) + PVEG(:)*PT2M(:)
!
!     Only diag
      PSNOWTEMP (:) = ZTN (:)
!
!
!                                            melting rate
!                                            there is melting only if T > T0 and
!                                            of course when SNOWSWE > 0.
!
      WHERE ( ZTN(:) > XTT .AND. PSNOWSWE(:) > 0.0 )
        PMELT(:) = ZPSN(:)*(ZTN(:)-XTT) / (PCS(:)*XLMTT*ZTAU_SNOWMELT)
      END WHERE
!
!                                            close the energy budget: cannot melt 
!                                            more than the futur available snow
!      
      ZNEXTSNOW(:) = PSNOWSWE(:) + PTSTEP * (PSR(:) - PLES(:) / PLSTT(:))
!
      WHERE ( PMELT(:) > 0.0 )
!              
              PMELT(:)=MIN(PMELT(:),ZNEXTSNOW(:)/PTSTEP)      
              ZNEXTSNOW(:) = ZNEXTSNOW(:) - PTSTEP * PMELT
!              
!             removes very small fraction
              ZFRAC(:) = SNOW_FRAC_GROUND(ZNEXTSNOW(:))
              WHERE(ZFRAC(:)<1.0E-4)
                    PMELT(:)     = PMELT(:) + ZNEXTSNOW(:) / PTSTEP       
              ENDWHERE   
!       
      ENDWHERE   
!    
    ELSEIF (HSNOW_ISBA == 'EBA') THEN
!    
      PMELT(:)=MIN( PSNOWSWE(:)/PTSTEP + PSR(:) - PLES(:)/PLSTT(:) , &
                    MAX(0.0,(PTG(:,1)-XTT))  / MAX(ZEPS1,PCT*PTSTEP) / XLMTT )
!
    ENDIF
!
!                                            new temperature Ts(t) after melting
!                                            (cooling due to the melting of the
!                                            snow)
!
  PTG(:,1) = PTG(:,1) - PCT(:)*XLMTT*PMELT(:)*PTSTEP
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       4.     EFFECT OF THE MELTING/FREEZING 
!               ON THE SOIL HEAT AND ICE CONTENTS
!               ('DIF' soil option)
!               ---------------------------------
!
! The values for the two coefficients (multiplied by VEG and LAI) 
! in the expression below are from 
! Giard and Bazile (2000), Mon. Wea. Rev.: they model the effect of insolation due to
! vegetation cover. This used by both 'DEF' (code blocks 3.-4.) and 'DIF' options.
!
WHERE(PLAI(:)/=XUNDEF .AND. PVEG(:)/=0.)
    ZKSFC_IVEG(:) = (1.0-ZINSOLFRZ_VEG*PVEG(:)) * (1.0-(PLAI(:)/ZINSOLFRZ_LAI))
ELSEWHERE
    ZKSFC_IVEG(:) = 1.0 ! No vegetation
ENDWHERE
!
!
IF (HISBA == 'DIF') THEN
!
! Determine soil phase changes and ice evolution:
!
!
   CALL ICE_SOILDIF(HDIF, PTSTEP, PTAUICE, ZKSFC_IVEG,          &
                 PSOILHCAPZ, PWSATZ, PMPOTSATZ, PBCOEFZ,        &
                 PTG, PWGI, PWG, PDELTAT, PALPHA, PN, PM,       &
                 PWRES                                          )
!
!
ELSE
!
!*       5.1    Melting/freezing normalized coefficient
!               ---------------------------------------
!
  ZKSOIL       = (0.5 * SQRT(XCONDI*XCI*XRHOLI*XDAY/XPI))/XLMTT

  DO JJ=1,SIZE(PTG,1)
!-------------------------------------------------------------------------------
!*       5.     EFFECT OF THE MELTING/FREEZING 
!               ON THE SURFACE-SOIL HEAT AND ICE CONTENTS
!               ('DEF' or Force-Restore soil option)
!               -----------------------------------------
!
!        5.0    Average soil-column porosity
!               ----------------------------
!               if Force-Restore option in use, then vertical
!               profiles of soil hydrological parameters are constant,
!               so use the values in uppermost element (arbitrary)
!
    ZWSAT_AVGZ(JJ) = PWSATZ(JJ,1)
!
!
!               Influence of vegetation insolation on surface:
!
    ZKSFC_FRZ(JJ) = ZKSOIL * ZKSFC_IVEG(JJ)
!
  ENDDO
!*       5.2    Water freezing
!               --------------
!
  IF(HSOILFRZ == 'LWT')THEN
!
! use option to control phase changes based on a relationship
! between the unfrozen liquid water content and temperature.
! Uses the Clapp and Hornberger model for water potential.
! The energy-limit method used by Boone et al. 2000 and
! Giard and Bazile (2000) is the default. 
!
    DO JJ=1,SIZE(PTG,1)
      ZMATPOT(JJ)   = MIN(PMPOTSATZ(JJ,1), XLMTT*(PTG(JJ,1)-XTT)/(XG*PTG(JJ,1)) )
      ZWGMIN(JJ)    = ZWSAT_AVGZ(JJ)*( (ZMATPOT(JJ)/PMPOTSATZ(JJ,1))**(-1./PBCOEFZ(JJ,1)) )

      ZMATPOT(JJ)   = PMPOTSATZ(JJ,1)*( (PWG(JJ,1)/ZWSAT_AVGZ(JJ))**(-PBCOEFZ(JJ,1)) )
      ZTGMAX(JJ)    = XLMTT*XTT/(XLMTT - XG* ZMATPOT(JJ))
    ENDDO
  ELSE
    ZWGMIN(:)    = XWGMIN
    ZTGMAX(:)    = XTT
  ENDIF
!
!
   ZDELTATI(:)  = PTG(:,1) - ZTGMAX(:) ! initial temperature depression 
!
! Limit phase change energy by energy actually available for phase changes:
! (this has little effect generally but was found to prevent oscillations
! between soil T and water/ice in rare circumstances, notably for large time steps).
! We take an average of this flux-limited value and actual
! (one could eventually put in a time dependence on the weights used in averaging,
! perhaps giving larger weight to limited flux for larger time steps...for now
! just a constant factor assumed)
!
   ZDELTAT(:)  = ZDELTATI(:)
   ZDELTAT(:)  = MIN(ZDELTAT(:), MAX(0.0,ZDT(:)))
   ZDELTAT(:)  = MAX(ZDELTAT(:), MIN(0.0,ZDT(:)))
   ZDELTAT(:)  = ZTWGHT*ZDELTAT(:) + (1.0-ZTWGHT)*ZDELTATI(:)   
   WHERE(ZDT(:)*ZDELTAT(:) < 0.0)ZDELTAT(:) = 0.0
!

  DO JJ=1,SIZE(PTG,1)
!    
    ZWORK2(JJ) = XRHOLW*PD_G(JJ,1)
    ZEFFIC(JJ)    = MAX(ZEFFIC_MIN,(PWG(JJ,1)-XWGMIN)/ZWSAT_AVGZ(JJ))
    ZFREEZING(JJ) = MIN( MAX(0.0,PWG(JJ,1)-ZWGMIN(JJ))*ZWORK2(JJ),    &  
                  ZKSFC_FRZ(JJ)*ZEFFIC(JJ)*MAX( -ZDELTAT(JJ), 0.) )
!
!*       5.3    Ground Ice melt
!               ---------------
!
    ZEFFIC(JJ)    =  MAX(ZEFFIC_MIN,PWGI(JJ,1)/(ZWSAT_AVGZ(JJ)-XWGMIN))
    ZICE_MELT(JJ) = MIN( PWGI(JJ,1)*ZWORK2(JJ),                      &
                  ZKSFC_FRZ(JJ)*ZEFFIC(JJ)*MAX( ZDELTAT(JJ), 0. ) )
!
!*       5.4    Ice reservoir evolution
!               -----------------------
!
! Melting of ice/freezing of water:
!
    ZWGI1(JJ) = PWGI(JJ,1) + (PTSTEP/PTAUICE(JJ))*(1.0-ZPSNG(JJ))*        &
              (ZFREEZING(JJ) - ZICE_MELT(JJ))/ZWORK2(JJ) 
!
!
    ZWGI1(JJ)  = MAX( ZWGI1(JJ) , 0.             )
    ZWGI1(JJ)  = MIN( ZWGI1(JJ) , ZWSAT_AVGZ(JJ)-XWGMIN)
!
! Time tendency:
!
    PDWGI1(JJ) = ZWGI1(JJ) - PWGI(JJ,1)
!
!
!*       5.5    Effect on temperature
!               ---------------------
!
    PTG(JJ,1)   = PTG(JJ,1) + PDWGI1(JJ)*XLMTT*PCT(JJ)*ZWORK2(JJ)
!
!-------------------------------------------------------------------------------
!
!*       6.     EFFECT OF THE MELTING/FREEZING 
!               ON THE DEEP-SOIL HEAT AND ICE CONTENTS
!               ('DEF' or Force-Restore soil option)
!               --------------------------------------
!
    ZWORK1(JJ) = PD_G(JJ,1)/PD_G(JJ,2)
!*       6.1  Available Deep ice content
!             --------------------------
!
    ZWIM(JJ) = ( PWGI(JJ,2) - ZWORK1(JJ) * PWGI(JJ,1) )  / ( 1. - ZWORK1(JJ) )
!
    ZWIM(JJ) = MAX(0.,ZWIM(JJ))  ! Just in case of round-off errors
!
!*       6.2  Deep liquid water content
!             -------------------------
!
    ZWM(JJ)  = ( PWG(JJ,2) - ZWORK1(JJ) * PWG(JJ,1) )  / ( 1. - ZWORK1(JJ) )
!
!*       6.3    Water freezing
!               --------------
!
! Total soil volumetric heat capacity [J/(m3 K)]:
!
    ZSOILHEATCAP(JJ) = XCL*XRHOLW*PWG(JJ,2)  +                           &
                     XCI*XRHOLI*PWGI(JJ,2) +                           &
                     XSPHSOIL*XDRYWGHT*(1.0-ZWSAT_AVGZ(JJ))*(1.0-ZWSAT_AVGZ(JJ))
!
! Soil thickness which corresponds to T2 (m): 2 times the diurnal
! surface temperature wave penetration depth as T2 is the average
! temperature for this layer:
!
    PTDIURN(JJ)   = MIN(PD_G(JJ,2), 4./(ZSOILHEATCAP(JJ)*PCG(JJ)))
!
! Effective soil ice penetration depth (m):
!
    ZICEEFF(JJ)   = (PWGI(JJ,2)/(PWGI(JJ,2)+PWG(JJ,2)))*PD_G(JJ,2)
!
    ZDT(JJ)       = PTG(JJ,2) - PT2M(JJ) ! temperature depression (K)
!
  ENDDO

  IF(HSOILFRZ == 'LWT')THEN
!
! as for the surface layer (above)JJ 
! Note also that if the 'DIF'
! soil option is not in force, then the soil parameters are assumed
! to be homogeneous (in the verticalJJ thus we use 1st element of 2nd dimension
! of the 2D-soil parameter arrays).
!
    DO JJ=1,SIZE(PTG,1)

       ZMATPOT(JJ)   = MIN(PMPOTSATZ(JJ,1), XLMTT*(PTG(JJ,2)-XTT)/(XG*PTG(JJ,2)) )
       ZWGMIN(JJ)    = ZWSAT_AVGZ(JJ)*( (ZMATPOT(JJ)/PMPOTSATZ(JJ,1))**(-1./PBCOEFZ(JJ,1)) )

       ZMATPOT(JJ)   = PMPOTSATZ(JJ,1)*( (PWG(JJ,2)/ZWSAT_AVGZ(JJ))**(-PBCOEFZ(JJ,1)) )
       ZTGMAX(JJ)    = XLMTT*XTT/(XLMTT - XG* ZMATPOT(JJ))
    ENDDO
  ELSE
    ZWGMIN(:)    = XWGMIN
    ZTGMAX(:)    = XTT
  ENDIF
!
   ZDELTATI(:)  = PTG(:,2) - ZTGMAX(:) ! initial temperature depression 
!
! Limit phase change energy by energy actually available for phase changes:
! (this has little effect generally but was found to prevent oscillations
! between soil T and water/ice in rare circumstances, notably for large time steps).
! We take an average of this flux-limited value and actual
! (one could eventually put in a time dependence on the weights used in averaging,
! perhaps giving larger weight to limited flux for larger time steps...for now
! just a constant factor assumed)
!
   ZDELTAT(:)  = ZDELTATI(:)
   ZDELTAT(:)  = MIN(ZDELTAT(:), MAX(0.0,ZDT(:)))
   ZDELTAT(:)  = MAX(ZDELTAT(:), MIN(0.0,ZDT(:)))
   ZDELTAT(:)  = ZTWGHT*ZDELTAT(:) + (1.0-ZTWGHT)*ZDELTATI(:)   
   WHERE(ZDT(:)*ZDELTAT(:) < 0.0)ZDELTAT(:) = 0.0
!
! Allow freezing by T2 up to a certain depth so that
! T2 energy can not be used to freeze soil water
! at levels sufficiently deep in the soil.
!
  DO JJ=1,SIZE(PTG,1)
  
    ZWORK1(JJ) = PD_G(JJ,1)/PD_G(JJ,2)
    ZWORK2(JJ) = XRHOLW*(PD_G(JJ,2)-PD_G(JJ,1))

    ZFREEZING(JJ) = 0.0
    IF (ZICEEFF(JJ) <= PTDIURN(JJ)) THEN
!
      ZEFFIC(JJ)    = MAX(ZEFFIC_MIN, MAX(0.0,ZWM(JJ) - XWGMIN)/ZWSAT_AVGZ(JJ))
      ZFREEZING(JJ) = MIN( MAX(0.0, ZWM(JJ) - ZWGMIN(JJ))* ZWORK2(JJ),            &
                     ZKSOIL*ZEFFIC(JJ)*MAX( -ZDELTAT(JJ) , 0. ) )
    ENDIF
!
!
!*       6.4    Ground Ice melt
!               ---------------
!
    ZEFFIC(JJ)    = MAX(ZEFFIC_MIN, ZWIM(JJ)/(ZWSAT_AVGZ(JJ)-XWGMIN))
    ZICE_MELT(JJ) = MIN( ZWIM(JJ)*ZWORK2(JJ),             &
                  ZKSOIL*ZEFFIC(JJ)*MAX( ZDELTAT(JJ) , 0. ) )
!
!
!*       6.5    Deep-part of deep-soil Ice reservoir evolution
!               ----------------------------------------------
!
    ZWIT(JJ)   = ZWIM(JJ) + (PTSTEP/PTAUICE(JJ))*(1.0-ZPSNG(JJ))*       &
               ((ZFREEZING(JJ) - ZICE_MELT(JJ))/ ZWORK2(JJ))
!
    ZWIT(JJ)   = MAX( ZWIT(JJ) , 0.             )
    ZWIT(JJ)   = MIN( ZWIT(JJ) , ZWSAT_AVGZ(JJ)-XWGMIN)
!
!
!*       6.6    Effect on temperature
!               ---------------------
!
    PTG(JJ,2) = PTG(JJ,2) + (ZWIT(JJ)-ZWIM(JJ))*XLMTT*PCG(JJ)*ZWORK2(JJ)
!
!*       6.7    Add reservoir evolution from surface freezing (WI2 contains WI1)
!               ----------------------------------------------------------------
!
    ZWGI2(JJ)  = (1.-ZWORK1(JJ))*ZWIT(JJ) +  ZWORK1(JJ)*ZWGI1(JJ)
!
    PDWGI2(JJ) = ZWGI2(JJ) - PWGI(JJ,2)
  ENDDO
ENDIF
!
IF (LHOOK) CALL DR_HOOK('ISBA_FLUXES',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE ISBA_FLUXES











