!     #########
      SUBROUTINE SNOW3L(HSNOWRES, TPTIME, OMEB, HIMPLICIT_WIND,           &
                PPALPHAN, PLAIV,                                          &
                PPEW_A_COEF, PPEW_B_COEF,                                 &
                PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,       &
                PSNOWSWE,PSNOWRHO,PSNOWHEAT,PSNOWALB,                     &
                PSNOWGRAN1,PSNOWGRAN2,PSNOWHIST,PSNOWAGE,                 &                
                PTSTEP,PPS,PSR,PRR,PPSN3L,                                &
                PTA,PTG,PSW_RAD,PQA,PVMOD,PLW_RAD, PRHOA,                 &
                PUREF,PEXNS,PEXNA,PDIRCOSZW,                              &
                PZREF,PZ0,PZ0EFF,PZ0H,PALB,                               &
                PCT,PSOILCOND,PD_G,                                       &
                PSNOWLIQ,PSNOWTEMP,PSNOWDZ,                               &
                PTHRUFAL,PGRNDFLUX,PEVAPCOR,PGFLXCOR,PSNOWSFCH,           &
                PDELHEATN, PDELHEATN_SFC,                                 &
                PSWNETSNOW,PSWNETSNOWS,PLWNETSNOW,PSNOWFLUX,              &
                PRNSNOW,PHSNOW,PGFLUXSNOW,                                &
                PHPSNOW,PLES3L,PLEL3L,PEVAP,PSNDRIFT,PRI,                 &
                PEMISNOW,PCDSNOW,PUSTAR,PCHSNOW,PSNOWHMASS,PQS,           &
                PPERMSNOWFRAC,PZENITH,PXLAT,PXLON,                        &
                OSNOWDRIFT,OSNOWDRIFT_SUBLIM                              )  
!     ##########################################################################
!
!!****  *SNOW3L*
!!
!!    PURPOSE
!!    -------
!
!     3-Layer snow scheme option (Boone and Etchevers 1999)
!     For shallow snow cover, Default method of Douville et al. (1995)
!     used with this option: Model "turns on" when snow sufficiently
!     deep/above some preset critical snow depth.
!
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
!!
!!    REFERENCE
!!    ---------
!!
!!    ISBA-ES: Boone and Etchevers (2001)
!!    ISBA: Belair (1995)
!!    ISBA: Noilhan and Planton (1989)
!!    ISBA: Noilhan and Mahfouf (1996)
!!
!!    AUTHOR
!!    ------
!!      A. Boone           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    7/99
!!      Modified by A.Boone 05/02 (code, not physics)
!!      Modified by A.Boone 11/04 i) maximum density limit imposed (although
!!                                rarely if ever reached), ii) check to
!!                                see if upermost layer completely sublimates
!!                                during a timestep (as snowpack becomes vanishly
!!                                thin), iii) impose maximum grain size limit
!!                                in radiation transmission computation.
!!
!!      Modified by B. Decharme  (03/2009): Consistency with Arpege permanent
!!                                          snow/ice treatment (LGLACIER for alb)
!!      Modified by A. Boone     (04/2010): Implicit coupling and replace Qsat and DQsat
!!                                          by Qsati and DQsati, respectively.
!!      Modified by E. Brun      (08/2012): Mass conservation in SNOW3LEVAPGONE
!!      Modified by B. Decharme  (08/2012): Loop optimization
!!      Modified by B. Decharme  (09/2012): New wind implicitation
!!      Modified by E. Brun      (10/2012): Bug in vertical snow redistribution
!!      Modified by B. Decharme  (08/2013): Qsat as argument (needed for coupling with atm)
!!                                          Add snow drift effect like in CROCUS
!!                                          Improve snow albedo
!!      Modified by A. Boone     (10/2014): MEB modifs: permit option to impose fluxes at sfc  
!!      Modified by P. Samuelsson(10/2014): MEB added an optional snow albedo calculation
!!      Modified by A. Boone     (10/2014): SNOW3LREFRZ and SNOW3LEVAPN edited to give
!!                                          better enthalpy conservation.
!!
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,     ONLY : XTT, XRHOLW, XLMTT, XCL, XDAY, XLVTT
USE MODD_SURF_PAR, ONLY : XUNDEF
!
USE MODE_SNOW3L
!
USE MODD_TYPE_DATE_SURF, ONLY: DATE_TIME
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP
!                                      PTSTEP    = time step of the integration
TYPE(DATE_TIME), INTENT(IN)         :: TPTIME      ! current date and time
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
 CHARACTER(LEN=*),     INTENT(IN)  :: HIMPLICIT_WIND   ! wind implicitation option
!                                                     ! 'OLD' = direct
!                                                     ! 'NEW' = Taylor serie, order 1
!
REAL, DIMENSION(:), INTENT(IN)    :: PPS, PTA, PSW_RAD, PQA,                       &
                                         PVMOD, PLW_RAD, PSR, PRR  
!                                      PSW_RAD = incoming solar radiation (W/m2)
!                                      PLW_RAD = atmospheric infrared radiation (W/m2)
!                                      PRR     = rain rate [kg/(m2 s)]
!                                      PSR     = snow rate (SWE) [kg/(m2 s)]
!                                      PTA     = atmospheric temperature at level za (K)
!                                      PVMOD   = modulus of the wind parallel to the orography (m/s)
!                                      PPS     = surface pressure
!                                      PQA     = atmospheric specific humidity
!                                                at level za
!
REAL, DIMENSION(:), INTENT(IN)    :: PSOILCOND, PD_G, PPSN3L, PCT, PPALPHAN, PLAIV
!                                      PCT       = thermal inertia of underlying soil (m2 K J-1)
!                                      PSOILCOND = soil thermal conductivity [W/(m K)]
!                                      PD_G      = Assumed first soil layer thickness (m)
!                                                  Used to calculate ground/snow heat flux
!                                      PPSN3L    = snow fraction
!                                      PPALPHAN  = Snow/canopy transition coefficient
!                                      PLAIV     = LAI of canopy
!
REAL, DIMENSION(:), INTENT(IN)    :: PZREF, PUREF, PEXNS, PEXNA, PDIRCOSZW, PRHOA, PZ0, PZ0EFF, &
                                       PALB, PZ0H, PPERMSNOWFRAC  
!                                      PZ0EFF    = roughness length for momentum
!                                      PZ0       = grid box average roughness length
!                                      PZ0H      = grid box average roughness length for heat
!                                      PZREF     = reference height of the first
!                                                  atmospheric level
!                                      PUREF     = reference height of the wind
!                                      PRHOA     = air density
!                                      PEXNS     = Exner function at surface
!                                      PEXNA     = Exner function at lowest atmos level
!                                      PDIRCOSZW = Cosinus of the angle between the
!                                                  normal to the surface and the vertical
!                                      PALB      = soil/vegetation albedo
!                                      PPERMSNOWFRAC  = fraction of permanet snow/ice
!
REAL, DIMENSION(:), INTENT(IN)      :: PPEW_A_COEF, PPEW_B_COEF,                   &
                                         PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF,      &
                                         PPEQ_B_COEF  
!                                      PPEW_A_COEF = wind coefficient (m2s/kg)
!                                      PPEW_B_COEF = wind coefficient (m/s)
!                                      PPET_A_COEF = A-air temperature coefficient
!                                      PPET_B_COEF = B-air temperature coefficient
!                                      PPEQ_A_COEF = A-air specific humidity coefficient
!                                      PPEQ_B_COEF = B-air specific humidity coefficient
!
REAL, DIMENSION(:), INTENT(IN)    :: PTG
!                                      PTG       = Surface soil temperature (effective
!                                                  temperature the of layer lying below snow)
REAL, DIMENSION(:), INTENT(INOUT) :: PSNOWALB
!                                      PSNOWALB = Prognostic surface snow albedo
!                                                 (does not include anything but
!                                                 the actual snow cover)
!
REAL, DIMENSION(:,:), INTENT(INOUT):: PSNOWHEAT, PSNOWRHO, PSNOWSWE
!                                      PSNOWHEAT = Snow layer(s) heat content (J/m2)
!                                      PSNOWRHO  = Snow layer(s) averaged density (kg/m3)
!                                      PSNOWSWE  = Snow layer(s) liquid Water Equivalent (SWE:kg m-2)
!
REAL, DIMENSION(:,:), INTENT(INOUT):: PSNOWGRAN1, PSNOWGRAN2, PSNOWHIST
!                                      PSNOWGRAN1 = Snow layers grain feature 1
!                                      PSNOWGRAN2 = Snow layer grain feature 2
!                                      PSNOWHIST  = Snow layer grain historical
!                                                   parameter (only for non
!                                                   dendritic snow)
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWAGE  ! Snow grain age
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PRNSNOW, PHSNOW, PGFLUXSNOW, PLES3L, PLEL3L, &
                                       PHPSNOW, PEVAP, PGRNDFLUX, PEMISNOW
!                                      PGRNDFLUX = soil/snow interface heat flux (W/m2)
!                                      PLES3L      = evaporation heat flux from snow (W/m2)
!                                      PLEL3L      = sublimation (W/m2)
!                                      PHPSNOW     = heat release from rainfall (W/m2)
!                                      PRNSNOW     = net radiative flux from snow (W/m2)
!                                      PHSNOW      = sensible heat flux from snow (W/m2)
!                                      PGFLUXSNOW  = net heat flux from snow (W/m2)
!                                      PEVAP       = total evaporative flux (kg/m2/s)
!                                      PEMISNOW    = snow surface emissivity
!
REAL, DIMENSION(:), INTENT(INOUT) :: PSWNETSNOW, PLWNETSNOW, PSWNETSNOWS
!                                      PSWNETSNOW = net shortwave radiation entering top of snowpack 
!                                                  (W m-2) Imposed if MEB=T, diagnosed herein if MEB=F
!                                      PSWNETSNOWS= net shortwave radiation in uppermost layer of snowpack 
!                                                  (W m-2) Imposed if MEB=T, diagnosed herein if MEB=F
!                                                   Used for surface energy budget diagnostics
!                                      PLWNETSNOW = net longwave radiation entering top of snowpack 
!                                                  (W m-2) Imposed if MEB=T, diagnosed herein if MEB=F
!
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWLIQ, PSNOWTEMP, PSNOWDZ
!                                      PSNOWLIQ  = Snow layer(s) liquid water content (m)
!                                      PSNOWTEMP = Snow layer(s) temperature (m)
!                                      PSNOWDZ   = Snow layer(s) thickness (m)
!
REAL, DIMENSION(:), INTENT(OUT)     :: PTHRUFAL, PEVAPCOR, PGFLXCOR, PSNOWFLUX, PSNOWSFCH, &
                                       PDELHEATN, PDELHEATN_SFC
!                                      PTHRUFAL  = rate that liquid water leaves snow pack:
!                                                  paritioned into soil infiltration/runoff
!                                                  by ISBA [kg/(m2 s)]
!                                      PGRNDFLUX = soil/snow interface heat flux (W/m2)
!                                      PEVAPCOR  = evaporation/sublimation correction term:
!                                                  extract any evaporation exceeding the
!                                                  actual snow cover (as snow vanishes)
!                                                  and apply it as a surface soil water
!                                                  sink. [kg/(m2 s)]
!                                      PGFLXCOR  = flux correction to underlying soil for vanishing snowpack
!                                                  (to put any energy excess from snow to soil) (W/m2)
!                                      PSNOWFLUX = heat flux between the surface and sub-surface 
!                                                  snow layers (W/m2)
!                                      PSNOWSFCH = snow surface layer pseudo-heating term owing to
!                                                  changes in grid thickness            (W m-2)
!                                      PDELHEATN = total snow heat content change in the surface layer (W m-2)
!                                      PDELHEATN_SFC = total snow heat content change during the timestep (W m-2)
!
REAL, DIMENSION(:), INTENT(OUT)     :: PCDSNOW, PUSTAR, PCDSNOW, PSNDRIFT
!                                      PCDSNOW     = drag coefficient for momentum over snow
!                                      PUSTAR      = friction velocity over snow (m/s)
!                                      PEVAP       = total evaporative flux (kg/m2/s)
!                                      PSNDRIFT    = blowing snow sublimation (kg/m2/s)
!
REAL, DIMENSION(:), INTENT(OUT)   :: PCHSNOW, PSNOWHMASS
!                                      PCHSNOW     = drag coefficient for heat over snow
!                                      PSNOWHMASS  = heat content change due to mass
!                                                    changes in snowpack (J/m2): for budget
!                                                    calculations only.
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRI, PQS
!                                    PRI = Ridcharson number
!                                    PQS = surface humidity
!
REAL, DIMENSION(:), INTENT(IN)    :: PZENITH ! solar zenith angle
REAL, DIMENSION(:), INTENT(IN)    :: PXLAT,PXLON ! LAT/LON after packing
!
LOGICAL, INTENT(IN)               :: OSNOWDRIFT, OSNOWDRIFT_SUBLIM ! activate snowdrift, sublimation during drift
!
!*      0.2    declarations of local variables
!
INTEGER                            :: JJ, JI     ! Loop control
!
INTEGER                            :: INI        ! number of point
INTEGER                            :: INLVLS     ! number of snow layers
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWTEMP, ZSCAP, ZSNOWDZN, ZSCOND,    &
                                                        ZRADSINK  
!                                      ZSNOWTEMP  = Snow layer(s) averaged temperature (K)
!                                      ZSCAP      = Snow layer(s) heat capacity [J/(K m3)]
!                                      ZSNOWDZN   = Updated snow layer thicknesses (m)
!                                      ZSCOND     = Snow layer(s) thermal conducivity [W/(m K)]
!                                      ZRADSINK   = Snow solar Radiation source terms (W/m2)
!
REAL, DIMENSION(SIZE(PTA))          :: ZSNOW, ZSFCFRZ, ZTSTERM1, ZTSTERM2,                   &
                                       ZCT, ZRA, ZSNOWTEMPO1  
!                                      ZSNOW      = Total snow depth (m)
!                                      ZCT        = inverse of the product of snow heat capacity
!                                                   and layer thickness [(m2 K)/J]
!                                      ZRA        = Surface aerodynamic resistance
!                                      ZTSTERM1,ZTSTERM2 = Surface energy budget coefficients
!                                      ZSNOWTEMPO1= value of uppermost snow temperature
!                                                   before time integration (K)
!
LOGICAL, DIMENSION(SIZE(PTA))       :: GSFCMELT
!                                      GSFCMELT   = FLAG if surface melt is occurring, used
!                                                   for surface albedo calculation.
!
REAL, DIMENSION(SIZE(PTA))          :: ZRSRA, ZDQSAT, ZQSAT, ZRADXS, ZMELTXS, ZLIQHEATXS, &
                                         ZLWUPSNOW, ZGRNDFLUXO, ZGRNDFLUXI, ZPSN3L
!                                      ZRSRA    = air density over aerodynamic resistance
!                                      ZDQSAT   = derrivative of saturation specific humidity
!                                      ZQSAT    = saturation specific humidity
!                                      ZRADXS   = shortwave radiation absorbed by soil surface
!                                                 (for thin snow sover) (W m-2)
!                                      ZMELTXS  = excess energy for snowmelt for vanishingly
!                                                 thin snowpacks: used to heat underlying surface
!                                                 in order to maintain energy conservation (W m-2)
!                                      ZLIQHEATXS = excess snowpack heating for vanishingly thin
!                                                 snow cover: add energy to snow/ground heat
!                                                 flux (W m-2)
!                                      ZLWUPSNOW = upwelling longwave raaditive flux (W m-2)
!                                      ZGRNDFLUXO= snow-ground flux before correction (W m-2)
!                                                  This is used simply to test if snow
!                                                  completely melts during a timestep.
!                                      ZGRNDFLUXI= for the case where the ground flux is imposed,
!                                                  this is the actual imposed value.
!                                      ZPSN3L    = snow fraction: different use if MEB "on".
!                                                  In this case, it is only used for Tg update
!                                                  since only this variable has a sub-grid relevance.
!
REAL, DIMENSION(SIZE(PTA))          :: ZUSTAR2_IC, ZTA_IC, ZQA_IC, ZWORK, ZWORK2,                  &
                                       ZPET_A_COEF_T, ZPEQ_A_COEF_T, ZPET_B_COEF_T, ZPEQ_B_COEF_T  
!                                      ZUSTAR2_IC    = implicit lowest atmospheric level friction (m2/s2)
!                                      ZTA_IC        = implicit lowest atmospheric level air temperature
!                                      ZQA_IC        = implicit lowest atmospheric level specific humidity
!                                      ZPET_A_COEF_T = transformed A-air temperature coefficient
!                                      ZPET_B_COEF_T = transformed B-air temperature coefficient
!                                      ZPEQ_A_COEF_T = transformed A-air specific humidity coefficient
!                                      ZPEQ_B_COEF_T = transformed B-air specific humidity coefficient
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),3)  :: ZSPECTRALALBEDO, ZSPECTRALWORK
!                                       ZSPECTRALALBEDO=spectral albedo (3 bands for now)
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWHEAT0
!                                      ZSNOWHEAT0    = initial snow heat content (saved for diagnostics) (J m-2)
!
! *      0.3    declarations of local parameters
!
REAL, PARAMETER                      :: ZSNOWDZMIN   = 0.0001   ! (m)
!                                       ZSNOWDZMIN   = minimum snow layer thickness for
!                                                      thermal calculations. Used to prevent
!                                                      numerical problems as snow becomes
!                                                      vanishingly thin.
!
REAL, PARAMETER                      :: ZSNOW_AGRAIN = 1.6e-4   ! (m)
REAL, PARAMETER                      :: ZSNOW_BGRAIN = 1.1e-13  ! (m13/kg4)
REAL, PARAMETER                      :: ZDSGRAIN_MAX = 2.796e-3 ! m
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! - - ---------------------------------------------------
!
!       0.     Initialization
!               --------------
! NOTE that snow layer thickness is used throughout this code: SWE
! is only used to diagnose the thickness at the beginning of this routine
! and it is updated at the end of this routine.
!
IF (LHOOK) CALL DR_HOOK('SNOW3L',0,ZHOOK_HANDLE)
PSNOWDZ(:,:)  = PSNOWSWE(:,:)/PSNOWRHO(:,:)
!
INI           = SIZE(PSNOWSWE(:,:),1)
INLVLS        = SIZE(PSNOWSWE(:,:),2)    ! total snow layers
!
ZUSTAR2_IC    = 0.0
ZTA_IC        = 0.0
ZQA_IC        = 0.0
!
ZSNOWHEAT0(:,:) = PSNOWHEAT(:,:) ! save initial heat content
!
!*       1.     Snow total depth
!               ----------------
!
ZSNOW(:) = 0.
DO JJ=1,INLVLS
   DO JI=1,INI
      ZSNOW(JI) = ZSNOW(JI) + PSNOWDZ(JI,JJ)  ! m
   ENDDO
END DO
!
ZWORK(:)=ZSNOW(:)
!
! Caluclate new snow albedo at time t
!
ZWORK2(:)=PSNOWALB(:)
!
CALL SNOW3LALB(OMEB,ZWORK2,ZSPECTRALALBEDO,PSNOWRHO(:,1),PSNOWAGE(:,1),PPERMSNOWFRAC,PPS)
ZSPECTRALALBEDO(:,1)=ZSPECTRALALBEDO(:,1)*PSNOWALB(:)/ZWORK2(:)
ZSPECTRALALBEDO(:,2)=ZSPECTRALALBEDO(:,2)*PSNOWALB(:)/ZWORK2(:)
ZSPECTRALALBEDO(:,3)=ZSPECTRALALBEDO(:,3)*PSNOWALB(:)/ZWORK2(:)
!
!*       2.     Snowfall
!               --------
!
! Caluclate uppermost density and thickness changes due to snowfall,
! and add heat content of falling snow
!
 CALL SNOW3LFALL(PTSTEP,PSR,PTA,PVMOD,ZSNOW,PSNOWRHO,PSNOWDZ,          &
                 PSNOWHEAT,PSNOWHMASS,PSNOWALB,PSNOWAGE,PPERMSNOWFRAC  )
!
! Caluclate new snow albedo at time t if snowfall
!
CALL SNOW3LALB(OMEB,ZWORK2,ZSPECTRALWORK,PSNOWRHO(:,1),PSNOWAGE(:,1),PPERMSNOWFRAC,PPS)
!
WHERE(ZWORK(:)==0.0.AND.PSR(:)>0.0)
     ZSPECTRALALBEDO(:,1)=ZSPECTRALWORK(:,1)
     ZSPECTRALALBEDO(:,2)=ZSPECTRALWORK(:,2)
     ZSPECTRALALBEDO(:,3)=ZSPECTRALWORK(:,3)
     PSNOWALB(:)=ZWORK2(:)
ENDWHERE
!
!*       3.     Update grid/discretization
!               --------------------------
! Reset grid to conform to model specifications:
!
 CALL SNOW3LGRID(ZSNOWDZN,ZSNOW,PSNOWDZ_OLD=PSNOWDZ)
!
! Mass/Heat redistribution:
!
 CALL SNOW3LTRANSF(ZSNOW,PSNOWDZ,ZSNOWDZN,PSNOWRHO,PSNOWHEAT,PSNOWAGE)
!
!
!*       4.     Liquid water content and snow temperature
!               -----------------------------------------
!
! First diagnose snow temperatures and liquid
! water portion of the snow from snow heat content:
!
ZSCAP(:,:)     = SNOW3LSCAP(PSNOWRHO)
!
ZSNOWTEMP(:,:) = XTT + ( ((PSNOWHEAT(:,:)/PSNOWDZ(:,:))                   &
                   + XLMTT*PSNOWRHO(:,:))/ZSCAP(:,:) )  
!
PSNOWLIQ(:,:)  = MAX(0.0,ZSNOWTEMP(:,:)-XTT)*ZSCAP(:,:)*                  &
                   PSNOWDZ(:,:)/(XLMTT*XRHOLW)  
!
ZSNOWTEMP(:,:) = MIN(XTT,ZSNOWTEMP(:,:))
!
!
!*       5.     Snow Compaction
!               ---------------
!
! Calculate snow density: compaction/aging: density increases
!
 CALL SNOW3LCOMPACTN(PTSTEP,PSNOWRHO,PSNOWDZ,ZSNOWTEMP,ZSNOW)
!
! Snow compaction and metamorphism due to drift
!
PSNDRIFT(:) = 0.0
IF (OSNOWDRIFT) THEN
   CALL SNOW3LDRIFT(PTSTEP,PVMOD,PTA,PQA,PPS,PRHOA,PSNOWRHO,&
                    PSNOWDZ,ZSNOW,OSNOWDRIFT_SUBLIM,PSNDRIFT)
ENDIF
!
! Update snow heat content (J/m2):
!
ZSCAP(:,:)     = SNOW3LSCAP(PSNOWRHO)
PSNOWHEAT(:,:) = PSNOWDZ(:,:)*( ZSCAP(:,:)*(ZSNOWTEMP(:,:)-XTT)        &
                   - XLMTT*PSNOWRHO(:,:) ) + XLMTT*XRHOLW*PSNOWLIQ(:,:)  
!
!
!*       6.     Solar radiation transmission
!               -----------------------------
!
! Heat source (-sink) term due to shortwave
! radiation transmission within the snowpack:
!
 CALL SNOW3LRAD(OMEB,ZSNOWDZMIN,PSW_RAD,PSNOWALB,      &
                ZSPECTRALALBEDO,PSNOWDZ,PSNOWRHO,PALB, &
                PPERMSNOWFRAC,PZENITH,PSWNETSNOW,      &
                PSWNETSNOWS,ZRADSINK,ZRADXS)  
!
!
!*       7.     Heat transfer and surface energy budget
!               ---------------------------------------
! Snow thermal conductivity:
!
 CALL SNOW3LTHRM(PSNOWRHO,ZSCOND,ZSNOWTEMP,PPS,PSNOWLIQ)
!
! Precipitation heating term:
! Rainfall renders it's heat to the snow when it enters
! the snowpack:
! NOTE: for now set to zero because we'd need to remove this heat from
!       the atmosphere to conserve energy. This can be done, but should
!       be tested within an atmos model first probably...NOTE, this
!       could be actiavted for use in OFFLINE mode however.
!
!PHPSNOW(:)     = PRR(:)*XCL*(MAX(XTT,PTA(:))-XTT)    ! (W/m2)
PHPSNOW(:) = 0.0
!
! Surface Energy Budget calculations using ISBA linearized form
! and standard ISBA turbulent transfer formulation
!
IF(OMEB)THEN 

! - surface fluxes prescribed:

   CALL SNOW3LEBUDMEB(PTSTEP,ZSNOWDZMIN,                                     &                
                ZSNOWTEMP(:,1),PSNOWDZ(:,1),PSNOWDZ(:,2),                    & 
                ZSCOND(:,1),ZSCOND(:,2),ZSCAP(:,1),                          &
                PRNSNOW,PHSNOW,PLES3L,PLEL3L,ZRADSINK(:,1),PHPSNOW,          & 
                ZCT,ZTSTERM1,ZTSTERM2,PGFLUXSNOW)

! - for mass fluxes, snow fraction notion removed

   ZPSN3L(:) = 1.0

! - as prescribed, we don't need output from here, but must initialize

   PCDSNOW(:) = 0.
   PCHSNOW(:) = 0.
   PUSTAR(:)  = 0.

ELSE

   ZPSN3L(:) = PPSN3L(:)

! - compute surface energy budget and fluxes:

   CALL SNOW3LEBUD(HSNOWRES, HIMPLICIT_WIND,                                   &
                  PPEW_A_COEF, PPEW_B_COEF,                                    &
                  PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,          &
                  ZSNOWDZMIN,                                                  &
                  PZREF,ZSNOWTEMP(:,1),PSNOWRHO(:,1),PSNOWLIQ(:,1),ZSCAP(:,1), &
                  ZSCOND(:,1),ZSCOND(:,2),                                     &
                  PUREF,PEXNS,PEXNA,PDIRCOSZW,PVMOD,                           &
                  PLW_RAD,PSW_RAD,PTA,PQA,PPS,PTSTEP,                          &
                  PSNOWDZ(:,1),PSNOWDZ(:,2),PSNOWALB,PZ0,PZ0EFF,PZ0H,          &
                  ZSFCFRZ,ZRADSINK(:,1),PHPSNOW,                               &
                  ZCT,PEMISNOW,PRHOA,ZTSTERM1,ZTSTERM2,ZRA,PCDSNOW,PCHSNOW,    &
                  ZQSAT, ZDQSAT, ZRSRA, ZUSTAR2_IC, PRI,                       &
                  ZPET_A_COEF_T,ZPEQ_A_COEF_T,ZPET_B_COEF_T,ZPEQ_B_COEF_T      )  
!
ENDIF
!
! Heat transfer: simple diffusion along the thermal gradient
!
ZSNOWTEMPO1(:) = ZSNOWTEMP(:,1) ! save surface snow temperature before update
!
ZGRNDFLUXI(:)  = PGRNDFLUX(:)
!
CALL SNOW3LSOLVT(OMEB,PTSTEP,ZSNOWDZMIN,PSNOWDZ,ZSCOND,ZSCAP,PTG,              &
                   PSOILCOND,PD_G,ZRADSINK,ZCT,ZTSTERM1,ZTSTERM2,              &
                   ZPET_A_COEF_T,ZPEQ_A_COEF_T,ZPET_B_COEF_T,ZPEQ_B_COEF_T,    &
                   ZTA_IC,ZQA_IC,PGRNDFLUX,ZGRNDFLUXO,ZSNOWTEMP,PSNOWFLUX      )  
!
!
!*       8.     Surface fluxes
!               --------------
! 
! Since surface fluxes already computed under MEB option (and already
! recomputed for the case when T>Tf), Only call if MEB not in use:
!
IF(.NOT.OMEB)THEN 

   CALL SNOW3LFLUX(ZSNOWTEMP(:,1),PSNOWDZ(:,1),PEXNS,PEXNA,           &
                  ZUSTAR2_IC,                                         &
                  PTSTEP,PSNOWALB,PSW_RAD,                            &
                  PEMISNOW,ZLWUPSNOW,PLW_RAD,PLWNETSNOW,              &
                  ZTA_IC,ZSFCFRZ,ZQA_IC,PHPSNOW,                      &
                  ZSNOWTEMPO1,PSNOWFLUX,ZCT,ZRADSINK(:,1),            &
                  ZQSAT,ZDQSAT,ZRSRA,                                 &
                  PRNSNOW,PHSNOW,PGFLUXSNOW,PLES3L,PLEL3L,PEVAP,      &
                  PUSTAR,GSFCMELT                                     )
!
ENDIF                  
!
!
!*       9.     Snow melt
!               ---------
!
! First Test to see if snow pack vanishes during this time step:
!
 CALL SNOW3LGONE(PTSTEP,PPSN3L,PLEL3L,PLES3L,PSNOWRHO,                     &
                  PSNOWHEAT,ZRADSINK(:,INLVLS),PEVAPCOR,PTHRUFAL,PGRNDFLUX, &
                  PGFLUXSNOW,ZGRNDFLUXO,PSNOWDZ,PSNOWLIQ,ZSNOWTEMP,ZRADXS)  
!
! Add radiation not absorbed by snow to soil/vegetation interface flux
! (for thin snowpacks):
!
PGRNDFLUX(:) = PGRNDFLUX(:) + ZRADXS(:)
!
! For "normal" melt: transform excess heat content into snow liquid:
!
 CALL SNOW3LMELT(PTSTEP,ZSCAP,ZSNOWTEMP,PSNOWDZ,PSNOWRHO,  &
                  PSNOWLIQ,ZMELTXS)  
!
! Add any excess heat for melting to underlying surface
! (for vanishingly thin snowpacks):
!
PGRNDFLUX(:) = PGRNDFLUX(:) + ZMELTXS(:)
!
!
!*      10.     Snow water flow and refreezing
!               ------------------------------
! Liquid water vertical transfer and possible snowpack runoff
! And refreezing/freezing of meltwater/rainfall (ripening of the snow)
!
 CALL SNOW3LREFRZ(PTSTEP,PRR,PSNOWRHO,ZSNOWTEMP,PSNOWDZ,PSNOWLIQ,PTHRUFAL)
!
ZSCAP(:,:)        = SNOW3LSCAP(PSNOWRHO)
PSNOWHEAT(:,:)    = PSNOWDZ(:,:)*( ZSCAP(:,:)*(ZSNOWTEMP(:,:)-XTT)        &
                      - XLMTT*PSNOWRHO(:,:) ) + XLMTT*XRHOLW*PSNOWLIQ(:,:) 
!
!
!*      11.     Snow Evaporation/Sublimation mass updates:
!               ------------------------------------------
!
CALL SNOW3LEVAPN(PPSN3L,PLES3L,PLEL3L,PTSTEP,ZSNOWTEMP(:,1),PSNOWRHO(:,1),  &
                   PSNOWDZ,PSNOWLIQ(:,1),PTA,PSNOWHEAT,                     &
                   PEVAPCOR,PGFLXCOR)
!
! If all snow in uppermost layer evaporates/sublimates, re-distribute
! grid (below could be evoked for vanishingly thin snowpacks):
!
 CALL SNOW3LEVAPGONE(PSNOWHEAT,PSNOWDZ,PSNOWRHO,ZSNOWTEMP,PSNOWLIQ)
!
!
!*      12.     Update surface albedo:
!               ----------------------
! Snow clear sky albedo:
!
 CALL SNOW3LALB(OMEB,PSNOWALB,ZSPECTRALALBEDO,PSNOWRHO(:,1),PSNOWAGE(:,1),PPERMSNOWFRAC,PPS)
!
!
!*      13.     Update snow heat content:
!               -------------------------
! Update the heat content (variable stored each time step)
! using current snow temperature and liquid water content:
!
! First, make check to make sure heat content not too large
! (this can result due to signifigant heating of thin snowpacks):
! add any excess heat to ground flux:
!
DO JJ=1,INLVLS
   DO JI=1,INI
      ZLIQHEATXS(JI) = MAX(0.0,PSNOWLIQ(JI,JJ)*XRHOLW-PSNOWDZ(JI,JJ)*PSNOWRHO(JI,JJ))*XLMTT/PTSTEP  
      PSNOWLIQ(JI,JJ)= PSNOWLIQ(JI,JJ) - ZLIQHEATXS(JI)*PTSTEP/(XRHOLW*XLMTT)
      PSNOWLIQ(JI,JJ)= MAX(0.0, PSNOWLIQ(JI,JJ))
      PGRNDFLUX(JI)  = PGRNDFLUX(JI)   + ZLIQHEATXS(JI)
   ENDDO
ENDDO
!
PSNOWTEMP(:,:)    = ZSNOWTEMP(:,:)
!
ZSCAP(:,:)        = SNOW3LSCAP(PSNOWRHO)
!
PSNOWHEAT(:,:)    = PSNOWDZ(:,:)*( ZSCAP(:,:)*(PSNOWTEMP(:,:)-XTT)        &
                      - XLMTT*PSNOWRHO(:,:) ) + XLMTT*XRHOLW*PSNOWLIQ(:,:)  
!
PSNOWSWE(:,:)     = PSNOWDZ(:,:)*PSNOWRHO(:,:)
!
!*      14.     Update snow mass and age:
!               -------------------------
!
PSNOWSWE(:,:) = PSNOWDZ(:,:)*PSNOWRHO(:,:)
!
WHERE (PSNOWSWE(:,:)>0)
  PSNOWAGE(:,:)=PSNOWAGE(:,:)+(PTSTEP/XDAY)
ELSEWHERE
  PSNOWAGE(:,:)= XUNDEF
ENDWHERE
!
!
!*      15.     Update surface specific humidity:
!               ---------------------------------
!
PQS(:) = ZQSAT(:)
!
!
!*      16.     Energy Budget Diagnostics:
!               --------------------------
! Final adjustment: If using MEB, then all excess heat (above initial guess ground flux)
! is used herein to adjust the ground temperature since for the MEB case, boundary
! fluxes are imposed (the energy budget of the ground has already been computed).
! This correction is needed since snow-soil flux here is implicit, while guess
! used in surface energy soil budget was semi-implicit. This forces the flux
! seen by the ground to be *the same* as that leaving the snow (energy conservation).
! Also, add any excessing cooling from sublimation as snowpack becomes vanishingly thin.
! This is added to underlying soil as a heating/cooling term (W/m2)
!
PGFLXCOR(:)     = (PGRNDFLUX(:) - ZGRNDFLUXI(:)) + PGFLXCOR(:)
!
!
! NOTE: To check enthalpy conservation, the error in W m-2 is defined HERE as
!
! error = (SUM(PSNOWHEAT(:,:),2)-SUM(ZSNOWHEAT0(:,:),2))/PTSTEP         &
!          - PSNOWHMASS(:)/PTSTEP - (PGFLUXSNOW(:) - PGRNDFLUX(:))
!
! where ZSNOWHEAT0 is the value of PSNOWHEAT that enters this routine before any adjustments/processes.
! Errors should be VERY small ( ~0) at each time step.
! The above is "snow relative"...the actual tile or grid box error is multiplied by PPSN3L.
!
! Snow heat storage change over the time step (W m-2):

PDELHEATN(:) = 0.0
DO JJ=1,INLVLS
   DO JI=1,INI
      PDELHEATN(JI) = PDELHEATN(JI) + (PSNOWHEAT(JI,JJ)-ZSNOWHEAT0(JI,JJ))
   ENDDO
END DO
PDELHEATN(:)     =  PDELHEATN(:)                   /PTSTEP
PDELHEATN_SFC(:) = (PSNOWHEAT(:,1)-ZSNOWHEAT0(:,1))/PTSTEP

!
! The pseudo-heating (W m-2) is derrived as a diagnostic from the surface energy budget
! here...NOTE that the total snow energy budget is explicitly computed, and since
! it balances, then simply compute this term as a residue (which ensures sfc energy balance):
!
PSNOWSFCH(:)     = PDELHEATN_SFC(:) - (PSWNETSNOWS(:) +PLWNETSNOW(:) - PHSNOW(:) -PLES3L(:)-PLEL3L(:))    &
                      + PSNOWFLUX(:) - PSNOWHMASS(:)/PTSTEP
!

IF (LHOOK) CALL DR_HOOK('SNOW3L',1,ZHOOK_HANDLE)
CONTAINS
!
!
!
!####################################################################
!####################################################################
!####################################################################
SUBROUTINE SNOW3LFALL(PTSTEP,PSR,PTA,PVMOD,PSNOW,PSNOWRHO,PSNOWDZ,        &
                      PSNOWHEAT,PSNOWHMASS,PSNOWALB,PSNOWAGE,PPERMSNOWFRAC)  
!
!!    PURPOSE
!!    -------
!     Calculate changes to snowpack resulting from snowfall.
!     Update mass and heat content of uppermost layer.
!
!
USE MODD_CSTS,     ONLY : XLMTT, XTT, XCI
USE MODD_SNOW_PAR, ONLY : XRHOSMIN_ES, XSNOWDMIN
!                   
USE MODE_SNOW3L
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP
!
REAL, DIMENSION(:), INTENT(IN)      :: PSR, PTA, PVMOD, PPERMSNOWFRAC
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PSNOW, PSNOWALB
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWRHO, PSNOWDZ, PSNOWHEAT, PSNOWAGE
!
REAL, DIMENSION(:), INTENT(OUT)     :: PSNOWHMASS
!
!
!*      0.2    declarations of local variables
!
INTEGER                             :: JJ, JI
!
INTEGER                             :: INI
INTEGER                             :: INLVLS
!
REAL, DIMENSION(SIZE(PTA))          :: ZSNOWFALL, ZRHOSNEW,        &
                                       ZSNOW, ZSNOWTEMP,           &
                                       ZSNOWFALL_DELTA, ZSCAP,     &
                                       ZAGENEW
!
! ISBA-ES CROCUS (Pahaut 1976): snowfall density coefficients:
!
REAL, PARAMETER                      :: ZSNOWFALL_A_SN = 109.0  ! kg/m3
REAL, PARAMETER                      :: ZSNOWFALL_B_SN =   6.0  ! kg/(m3 K)
REAL, PARAMETER                      :: ZSNOWFALL_C_SN =  26.0  ! kg/(m7/2 s1/2)
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! 0. Initialize:
! ------------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3LFALL',0,ZHOOK_HANDLE)
!
INI             = SIZE(PSNOWDZ(:,:),1)
INLVLS          = SIZE(PSNOWDZ(:,:),2)
!
ZRHOSNEW(:)     = XRHOSMIN_ES
ZAGENEW (:)     = 0.0
ZSNOWFALL(:)    = 0.0
ZSCAP(:)        = 0.0
ZSNOW(:)        = PSNOW(:)
!
PSNOWHMASS(:)   = 0.0
!
! 1. Incorporate snowfall into snowpack:
! --------------------------------------
!
!
! Heat content of newly fallen snow (J/m2):
! NOTE for now we assume the snowfall has
! the temperature of the snow surface upon reaching the snow.
! This is done as opposed to using the air temperature since
! this flux is quite small and has little to no impact
! on the time scales of interest. If we use the above assumption
! then, then the snowfall advective heat flux is zero.
!
ZSNOWTEMP(:)  = XTT
!
WHERE (PSR(:) > 0.0 .AND. PSNOWDZ(:,1)>0.)
  ZSCAP(:)      = SNOW3LSCAP(PSNOWRHO(:,1))
  ZSNOWTEMP(:)  = XTT + (PSNOWHEAT(:,1) +                              &
                    XLMTT*PSNOWRHO(:,1)*PSNOWDZ(:,1))/                   &
                    (ZSCAP(:)*MAX(XSNOWDMIN/INLVLS,PSNOWDZ(:,1)))  
  ZSNOWTEMP(:)  = MIN(XTT, ZSNOWTEMP(:))
END WHERE
!
WHERE (PSR(:) > 0.0)
!
  PSNOWHMASS(:) = PSR(:)*(XCI*(ZSNOWTEMP(:)-XTT)-XLMTT)*PTSTEP
!
! Snowfall density: Following CROCUS (Pahaut 1976)
!
   ZRHOSNEW(:)   = MAX(XRHOSMIN_ES, ZSNOWFALL_A_SN + ZSNOWFALL_B_SN*(PTA(:)-XTT)+         &
                     ZSNOWFALL_C_SN*SQRT(PVMOD(:)))  
!
!
! Fresh snowfall changes the snowpack age,
! decreasing in uppermost snow layer (mass weighted average):
!
   PSNOWAGE(:,1) = (PSNOWAGE(:,1)*PSNOWDZ(:,1)*PSNOWRHO(:,1)+ZAGENEW(:)*PSR(:)*PTSTEP) / &
                   (PSNOWDZ(:,1)*PSNOWRHO(:,1)+PSR(:)*PTSTEP)
!
! Augment total pack depth:
!
   ZSNOWFALL(:)  = PSR(:)*PTSTEP/ZRHOSNEW(:)    ! snowfall thickness (m)
!
   PSNOW(:)      = PSNOW(:) + ZSNOWFALL(:)
!
! Fresh snowfall changes the snowpack
! density, increases the total liquid water
! equivalent: in uppermost snow layer:
!
   PSNOWRHO(:,1) = (PSNOWDZ(:,1)*PSNOWRHO(:,1) + ZSNOWFALL(:)*ZRHOSNEW(:))/     &
                   (PSNOWDZ(:,1)+ZSNOWFALL(:))  
!
   PSNOWDZ(:,1)  = PSNOWDZ(:,1) + ZSNOWFALL(:)
!
! Add energy of snowfall to snowpack:
! Update heat content (J/m2) (therefore the snow temperature
! and liquid content):
!
   PSNOWHEAT(:,1)  = PSNOWHEAT(:,1) + PSNOWHMASS(:)
!
END WHERE
!
!
! 2. Case of new snowfall on a previously snow-free surface:
! ----------------------------------------------------------
!
! When snow first falls on a surface devoid of snow,
! redistribute the snow mass throughout the 3 layers:
! (temperature already set in the calling routine
! for this case)
!
ZSNOWFALL_DELTA(:)    = 0.0
WHERE(ZSNOW(:) == 0.0 .AND. PSR(:) > 0.0)
   ZSNOWFALL_DELTA(:) = 1.0
END WHERE
!
DO JJ=1,INLVLS
   DO JI=1,INI
!
      PSNOWDZ(JI,JJ)   = ZSNOWFALL_DELTA(JI)*(ZSNOWFALL(JI) /INLVLS) + &
                        (1.0-ZSNOWFALL_DELTA(JI))*PSNOWDZ(JI,JJ)  
!
      PSNOWHEAT(JI,JJ) = ZSNOWFALL_DELTA(JI)*(PSNOWHMASS(JI)/INLVLS) + &
                       (1.0-ZSNOWFALL_DELTA(JI))*PSNOWHEAT(JI,JJ)  
!
      PSNOWRHO(JI,JJ)  = ZSNOWFALL_DELTA(JI)*ZRHOSNEW(JI)            + &
                       (1.0-ZSNOWFALL_DELTA(JI))*PSNOWRHO(JI,JJ)
!
      PSNOWAGE(JI,JJ)  = ZSNOWFALL_DELTA(JI)*(ZAGENEW(JI)/INLVLS)    + &
                       (1.0-ZSNOWFALL_DELTA(JI))*PSNOWAGE(JI,JJ)  
!
   ENDDO
ENDDO
!
IF (LHOOK) CALL DR_HOOK('SNOW3LFALL',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE SNOW3LFALL
!####################################################################
!####################################################################
!####################################################################
SUBROUTINE SNOW3LCOMPACTN(PTSTEP,PSNOWRHO,PSNOWDZ,PSNOWTEMP,PSNOW)  
!
!!    PURPOSE
!!    -------
!     Snow compaction due to overburden and settling.
!     Mass is unchanged: layer thickness is reduced
!     in proportion to density increases. Method
!     of Anderson (1976): see Loth and Graf, 1993,
!     J. of Geophys. Res., 98, 10,451-10,464.
!
!
USE MODD_CSTS,     ONLY : XTT, XG
USE MODD_SNOW_PAR, ONLY : XRHOSMAX_ES
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWTEMP
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWRHO, PSNOWDZ
!
REAL, DIMENSION(:), INTENT(OUT)     :: PSNOW
!
!
!*      0.2    declarations of local variables
!
INTEGER                             :: JJ, JI
!
INTEGER                             :: INI
INTEGER                             :: INLVLS
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWRHO2, ZSETTLE, &
                                                        ZVISCOCITY, ZSMASS, &
                                                        ZWSNOWDZ  
!
REAL, DIMENSION(SIZE(PSNOW))         :: ZSMASSC
!
! ISBA-ES Compaction/Settling Coefficients from Loth and Graf (1993)
! following Mellor (1964) and Anderson (1976):
! see Boone, Meteo-France/CNRM Note de Centre No. 70 (2002)
!
REAL, PARAMETER                      :: ZSNOWCMPCT_RHOD  = 150.0  ! (kg/m3)
REAL, PARAMETER                      :: ZSNOWCMPCT_ACM   = 2.8e-6 ! (1/s)
REAL, PARAMETER                      :: ZSNOWCMPCT_BCM   = 0.04   ! (1/K)
REAL, PARAMETER                      :: ZSNOWCMPCT_CCM   = 0.460  ! (m3/kg)
REAL, PARAMETER                      :: ZSNOWCMPCT_V0    = 3.7e7  ! (Pa/s)
REAL, PARAMETER                      :: ZSNOWCMPCT_VT    = 0.081  ! (1/K)
REAL, PARAMETER                      :: ZSNOWCMPCT_VR    = 0.018  ! (m3/kg)
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! 0. Initialization:
! ------------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3LCOMPACTN',0,ZHOOK_HANDLE)
!
ZSNOWRHO2(:,:)  = PSNOWRHO(:,:)
ZSETTLE(:,:)    = ZSNOWCMPCT_ACM
ZVISCOCITY(:,:) = ZSNOWCMPCT_V0
!
INI             = SIZE(PSNOWDZ(:,:),1)
INLVLS          = SIZE(PSNOWDZ(:,:),2)
!
!
! 1. Cumulative snow mass (kg/m2):
! --------------------------------
!
ZWSNOWDZ(:,:) = PSNOWDZ(:,:)*PSNOWRHO(:,:)
ZSMASSC(:)    = 0.0
DO JJ=1,INLVLS
   DO JI=1,INI
      ZSMASS(JI,JJ)    = ZSMASSC(JI) + ZWSNOWDZ(JI,JJ)
      ZSMASSC(JI)      = ZSMASSC(JI) + ZWSNOWDZ(JI,JJ)
   ENDDO
ENDDO
!
! 2. Compaction/Settling
! ----------------------
! Compaction/settling if density below upper limit
! (compaction is generally quite small above ~ 500 kg m-3):
!
DO JJ=1,INLVLS
   DO JI=1,INI
!   
      IF(PSNOWRHO(JI,JJ) < XRHOSMAX_ES)THEN
!              
!       First calculate settling due to freshly fallen snow:
        ZSETTLE(JI,JJ) = ZSNOWCMPCT_ACM*EXP( -ZSNOWCMPCT_BCM*(XTT-PSNOWTEMP(JI,JJ))  &
                       - ZSNOWCMPCT_CCM*MAX(0.0,PSNOWRHO(JI,JJ)-ZSNOWCMPCT_RHOD)  )  
!                       
!       Snow viscocity:
        ZVISCOCITY(JI,JJ) = ZSNOWCMPCT_V0*EXP( ZSNOWCMPCT_VT*(XTT-PSNOWTEMP(JI,JJ))  &
                          + ZSNOWCMPCT_VR*PSNOWRHO(JI,JJ) )                          
!                          
!       Calculate snow density: compaction from weight/over-burden Anderson 1976 method:
        ZSNOWRHO2(JI,JJ) = PSNOWRHO(JI,JJ) + PSNOWRHO(JI,JJ)*PTSTEP &
                         * ( (XG*ZSMASS(JI,JJ)/ZVISCOCITY(JI,JJ))+ ZSETTLE(JI,JJ) )
!                         
!       Conserve mass by decreasing grid thicknesses in response to density increases
        PSNOWDZ(JI,JJ) = PSNOWDZ(JI,JJ)*(PSNOWRHO(JI,JJ)/ZSNOWRHO2(JI,JJ))  
!        
      ENDIF
   ENDDO
ENDDO
!
! 3. Update total snow depth and density profile:
! -----------------------------------------------
!
! Compaction/augmentation of total snowpack depth
!
PSNOW(:) = 0.
DO JJ=1,INLVLS
   DO JI=1,INI
      PSNOW(JI) = PSNOW(JI) + PSNOWDZ(JI,JJ)
   ENDDO
ENDDO
!
! Update density (kg m-3):
!
PSNOWRHO(:,:)  = ZSNOWRHO2(:,:)
!
IF (LHOOK) CALL DR_HOOK('SNOW3LCOMPACTN',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOW3LCOMPACTN
!####################################################################
!####################################################################
!####################################################################
SUBROUTINE SNOW3LDRIFT(PTSTEP,PVMOD,PTA,PQA,PPS,PRHOA,PSNOWRHO,&
                       PSNOWDZ,PSNOW,OSNOWDRIFT_SUBLIM,PSNDRIFT)
!
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODE_THERMOS
!
!
!!    PURPOSE
!!    -------
!     Snow compaction  and metamorphism due to drift
!     Mass is unchanged: layer thickness is reduced
!     in proportion to density increases. Method inspired from
!     Brun et al. (1997) and Guyomarch
!     
!     - computes a mobility index of each snow layer from its  density
!     - computes a drift index of each layer from its mobility and wind speed
!     - computes a transport index 
!     - increases density in case of transport
!
!     HISTORY:
!     Basic parameterization from Crocus/ARPEGE Coupling (1997) 
!     Adaptation from Crocus SNOWDRIFT subroutine
!
!!    PURPOSE
!!    -------
!     Snow compaction due to overburden and settling.
!     Mass is unchanged: layer thickness is reduced
!     in proportion to density increases. Method
!     of Anderson (1976): see Loth and Graf, 1993,
!     J. of Geophys. Res., 98, 10,451-10,464.
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP
!
REAL, DIMENSION(:), INTENT(IN)      :: PVMOD
REAL, DIMENSION(:), INTENT(IN)      :: PTA
REAL, DIMENSION(:), INTENT(IN)      :: PQA
REAL, DIMENSION(:), INTENT(IN)      :: PPS
REAL, DIMENSION(:), INTENT(IN)      :: PRHOA
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWRHO, PSNOWDZ
!
REAL, DIMENSION(:), INTENT(OUT)     :: PSNOW
!
LOGICAL,            INTENT(IN)      :: OSNOWDRIFT_SUBLIM
REAL, DIMENSION(:), INTENT(OUT)     :: PSNDRIFT !blowing snow sublimation (kg/m2/s)
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWRHO2
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZRMOB
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZRDRIFT
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZRT
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZDRO
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZQS_EFFECT
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZDRIFT_EFFECT
REAL, DIMENSION(SIZE(PSNOWRHO,1)                 ) :: ZPROFEQU
REAL, DIMENSION(SIZE(PSNOWRHO,1)                 ) :: ZWIND
REAL, DIMENSION(SIZE(PSNOWRHO,1)                 ) :: ZQSATI
REAL, DIMENSION(SIZE(PSNOWRHO,1)                 ) :: ZVT
REAL, DIMENSION(SIZE(PSNOWRHO,1)                 ) :: ZQS
REAL, DIMENSION(SIZE(PSNOWRHO,1)                 ) :: ZW
REAL, DIMENSION(SIZE(PSNOWRHO,1)                 ) :: ZT
REAL, DIMENSION(SIZE(PSNOWRHO,1)                 ) :: ZSNOWDZ1
!
INTEGER                             :: JJ, JI
!
INTEGER                             :: INI
INTEGER                             :: INLVLS
!
!*      0.3    declarations of Calibration coefficients
!
REAL, PARAMETER                     ::ZVTIME    = 48.*3600. ! characteristic time for compaction and metamorphism by wind drift
REAL, PARAMETER                     ::ZVROMAX   = 350.      ! maximum density for drift compaction (kg m-3)
REAL, PARAMETER                     ::ZVROMIN   = 50.       ! minimum density for mobility computation (kg m-3)
REAL, PARAMETER                     ::ZVMOB1    = 0.295     ! coef for computing the mobility index
REAL, PARAMETER                     ::ZVDRIFT1  = 2.868     ! coef for computing the drift index
REAL, PARAMETER                     ::ZVDRIFT2  = 0.085     ! coef for computing the drift index
REAL, PARAMETER                     ::ZVDRIFT3  = 3.25      ! coef for computing the drift index
!
! modif_EB pour sublim 
! a pour but de tenir compte du fait que le vent moyen est > rafales
! on en tient compte egalement pour diminuer la duree de l'effet
REAL, PARAMETER :: ZCOEF_FF    =1.25 ! coefficient for gust diagnosis from average wind 
REAL, PARAMETER :: ZCOEF_EFFECT=1.0  ! coefficient for impact on density du drift
REAL, PARAMETER :: ZQS_REF     =2.E-5! valeur de reference de ZQS pour effet neige
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3LDRIFT',0,ZHOOK_HANDLE)
!
! 0. Initialization:
! ------------------
!
INI    = SIZE(PSNOWDZ(:,:),1)
INLVLS = SIZE(PSNOWDZ(:,:),2)
!
ZSNOWRHO2    (:,:) = PSNOWRHO(:,:)
ZRDRIFT      (:,:) = XUNDEF
ZRT          (:,:) = XUNDEF
ZDRO         (:,:) = XUNDEF
ZQS_EFFECT   (:,:) = 0.0
ZDRIFT_EFFECT(:,:) = 0.0
!
ZSNOWDZ1(:) = PSNOWDZ(:,1)
!
! 1. Initialazation of drift and induced settling
! --------------------------------
!
ZWIND    (:) = ZCOEF_FF * PVMOD(:)
ZPROFEQU (:) = 0.0
!
! mobility index computation of a layer as a function of its density 
!
ZRMOB(:,:)= 1.25-1.25E-3*(MAX(PSNOWRHO(:,:),ZVROMIN)-ZVROMIN) / ZVMOB1
!
! 2. Specific case for blowing snow sublimation
! ---------------------------------------------
!
IF(OSNOWDRIFT_SUBLIM)THEN
!
  ZQSATI(:)= QSATI(PTA(:),PPS(:))
!
! computation of wind speed threshold QSATI and RH withe respect to ice
!
  ZW (:)=MAX(-0.99,ZRMOB(:,1))
  ZVT(:)=LOG(ZVDRIFT1/(1.0+ZW(:)))/ZVDRIFT2
!
  ZW(:)=LOG(ZWIND(:)/ZVT(:))
  ZW(:)=EXP(3.6*ZW(:))
!
! computation of sublimation rate according to Gordon's PhD
!
  ZT(:)=XTT/PTA(:)
  ZT(:)=0.0018*(ZT(:)**4)
!
  ZQS(:)=ZT(:)*ZVT(:)*PRHOA(:)*ZQSATI(:)*(1.-PQA(:)/ZQSATI(:))*ZW(:)
!
! surface depth decrease in case of blowing snow sublimation
!
  PSNOWDZ(:,1)=MAX(0.5*PSNOWDZ(:,1),PSNOWDZ(:,1)-MAX(0.,ZQS(:))*PTSTEP/(ZCOEF_FF*PSNOWRHO(:,1)))
!
  PSNDRIFT(:) = (ZSNOWDZ1(:)-PSNOWDZ(:,1))*PSNOWRHO(:,1)/PTSTEP
!
  ZQS_EFFECT(:,1)=MIN(3.,MAX(0.,ZQS(:))/ZQS_REF)
!
ENDIF
!
! 3. Computation of drift and induced settling
! --------------------------------
!
DO JJ=1,INLVLS
   DO JI=1,INI
!      
!     computation of the drift index inclunding the decay by overburden snow 
      ZRDRIFT(JI,JJ) =  ZRMOB(JI,JJ)-(ZVDRIFT1*EXP(-ZVDRIFT2*ZWIND(JI))-1.)

!     update the decay coeff by half the current layer
      ZPROFEQU(JI) = ZPROFEQU(JI)  + 0.5 * PSNOWDZ(JI,JJ) * 0.1 * (ZVDRIFT3-ZRDRIFT(JI,JJ))
!      
!     computation of the drift index inclunding the decay by overburden snow 
      ZRT(JI,JJ) = MAX(0.0,ZRDRIFT(JI,JJ)) * EXP(-ZPROFEQU(JI)*100.0)         
!     
      ZDRIFT_EFFECT(JI,JJ) = (ZQS_EFFECT(JI,JJ)+ZCOEF_EFFECT)*ZRT(JI,JJ)/(ZVTIME*ZCOEF_FF)
!
      IF(PSNOWRHO(JI,JJ) < ZVROMAX )THEN
!
!       settling by wind transport only in case of not too dense snow
        ZDRO(JI,JJ) = (ZVROMAX - PSNOWRHO(JI,JJ)) * ZDRIFT_EFFECT(JI,JJ) * PTSTEP
!          
!       Calculate new snow density:
        ZSNOWRHO2(JI,JJ) = MIN(ZVROMAX,PSNOWRHO(JI,JJ)+ZDRO(JI,JJ))

!       Conserve mass by decreasing grid thicknesses in response to density increases
        PSNOWDZ(JI,JJ) = PSNOWDZ(JI,JJ)*(PSNOWRHO(JI,JJ)/ZSNOWRHO2(JI,JJ))    
!
      ENDIF
!       
!     update the decay coeff by half the current layer
      ZPROFEQU(JI) = ZPROFEQU(JI)  + 0.5 * PSNOWDZ(JI,JJ) * 0.1 * (ZVDRIFT3-ZRDRIFT(JI,JJ))
!      
   ENDDO
ENDDO
!
! 4. Update total snow depth and density profile:
! -----------------------------------------------
!
! Compaction/augmentation of total snowpack depth
!
PSNOW(:) = 0.
DO JJ=1,INLVLS
   DO JI=1,INI
      PSNOW(JI) = PSNOW(JI) + PSNOWDZ(JI,JJ)
   ENDDO
ENDDO
!
! Update density (kg m-3):
!
PSNOWRHO(:,:)  = ZSNOWRHO2(:,:)
!
IF (LHOOK) CALL DR_HOOK('SNOW3LDRIFT',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOW3LDRIFT
!####################################################################
!####################################################################
!####################################################################
        SUBROUTINE SNOW3LTRANSF(PSNOW,PSNOWDZ,PSNOWDZN,    &
                                PSNOWRHO,PSNOWHEAT,PSNOWAGE)  
!
!!    PURPOSE
!!    -------
!     Snow mass,heat and characteristics redistibution in case of
!     grid resizing. Total mass and heat content of the overall snowpack
!     unchanged/conserved within this routine.
!     Same method as in Crocus
!
USE MODD_SNOW_PAR, ONLY : XSNOWCRITD
!
IMPLICIT NONE
!
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:  ), INTENT(IN)    :: PSNOW
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWDZN  
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWHEAT
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWRHO
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWDZ
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWAGE
!
!*      0.2    declarations of local variables
!
INTEGER                             :: JI, JL, JLO
!
INTEGER                             :: INI
INTEGER                             :: INLVLS
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWRHON
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWHEATN
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWAGEN
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWZTOP_NEW
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWZBOT_NEW                                       
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWRHOO
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWHEATO
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWAGEO
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWDZO
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWZTOP_OLD
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWZBOT_OLD  
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWHEAN
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWAGN
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZMASTOTN
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZMASSDZO
!
REAL, DIMENSION(SIZE(PSNOW)) :: ZPSNOW_OLD, ZPSNOW_NEW
REAL, DIMENSION(SIZE(PSNOW)) :: ZSUMHEAT, ZSUMSWE, ZSUMAGE, ZSNOWMIX_DELTA
!
REAL :: ZPROPOR
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! 0. Initialization:
! ------------------
!
!
IF (LHOOK) CALL DR_HOOK('SNOW3LTRANSF',0,ZHOOK_HANDLE)
!
INI        = SIZE(PSNOWRHO,1)
INLVLS     = SIZE(PSNOWRHO,2)
!
ZPSNOW_NEW(:) = 0.0
ZPSNOW_OLD(:) = PSNOW(:)
!
DO JL=1,INLVLS
   DO JI=1,INI
      ZPSNOW_NEW(JI)=ZPSNOW_NEW(JI)+PSNOWDZN(JI,JL)
   ENDDO
ENDDO
!
! initialization of variables describing the initial snowpack 
!
ZSNOWDZO  (:,:) = PSNOWDZ  (:,:)
ZSNOWRHOO (:,:) = PSNOWRHO (:,:)
ZSNOWHEATO(:,:) = PSNOWHEAT(:,:)
ZSNOWAGEO (:,:) = PSNOWAGE (:,:)
ZMASSDZO  (:,:) = XUNDEF
!
! 1. Calculate vertical grid limits (m):
! --------------------------------------
!
ZSNOWZTOP_OLD(:,1) = ZPSNOW_OLD(:)
ZSNOWZTOP_NEW(:,1) = ZPSNOW_NEW(:)
ZSNOWZBOT_OLD(:,1) = ZSNOWZTOP_OLD(:,1)-ZSNOWDZO(:,1)
ZSNOWZBOT_NEW(:,1) = ZSNOWZTOP_NEW(:,1)-PSNOWDZN(:,1)
!
DO JL=2,INLVLS
   DO JI=1,INI
      ZSNOWZTOP_OLD(JI,JL) = ZSNOWZBOT_OLD(JI,JL-1)
      ZSNOWZTOP_NEW(JI,JL) = ZSNOWZBOT_NEW(JI,JL-1)
      ZSNOWZBOT_OLD(JI,JL) = ZSNOWZTOP_OLD(JI,JL  )-ZSNOWDZO(JI,JL)
      ZSNOWZBOT_NEW(JI,JL) = ZSNOWZTOP_NEW(JI,JL  )-PSNOWDZN(JI,JL)
   ENDDO
ENDDO
ZSNOWZBOT_OLD(:,INLVLS)=0.0
ZSNOWZBOT_NEW(:,INLVLS)=0.0
!
! 3. Calculate mass, heat, charcateristics mixing due to vertical grid resizing:
! --------------------------------------------------------------------
!
! loop over the new snow layers
! Summ or avergage of the constituting quantities of the old snow layers
! which are totally or partially inserted in the new snow layer
! For snow age, mass weighted average is used.
!
ZSNOWHEAN(:,:)=0.0
ZMASTOTN (:,:)=0.0
ZSNOWAGN (:,:)=0.0
!
DO JL=1,INLVLS
   DO JLO=1, INLVLS   
      DO JI=1,INI
        IF((ZSNOWZTOP_OLD(JI,JLO)>ZSNOWZBOT_NEW(JI,JL)).AND.(ZSNOWZBOT_OLD(JI,JLO)<ZSNOWZTOP_NEW(JI,JL)))THEN
!                
          ZPROPOR = (MIN(ZSNOWZTOP_OLD(JI,JLO), ZSNOWZTOP_NEW(JI,JL)) &
                  -  MAX(ZSNOWZBOT_OLD(JI,JLO), ZSNOWZBOT_NEW(JI,JL)))&
                  / ZSNOWDZO(JI,JLO) 
!
          ZMASSDZO (JI,JLO)=ZSNOWRHOO(JI,JLO)*ZSNOWDZO(JI,JLO)*ZPROPOR
!
          ZMASTOTN (JI,JL)=ZMASTOTN (JI,JL)+ZMASSDZO  (JI,JLO)
          ZSNOWAGN (JI,JL)=ZSNOWAGN (JI,JL)+ZSNOWAGEO (JI,JLO)*ZMASSDZO(JI,JLO)
!
          ZSNOWHEAN(JI,JL)=ZSNOWHEAN(JI,JL)+ZSNOWHEATO(JI,JLO)*ZPROPOR
!          
        ENDIF
      ENDDO 
    ENDDO 
ENDDO  
!
! the new layer inherits from the weighted average properties of the old ones
! heat and mass
!
ZSNOWHEATN(:,:)= ZSNOWHEAN(:,:)
ZSNOWAGEN (:,:)= ZSNOWAGN (:,:)/ZMASTOTN(:,:)
ZSNOWRHON (:,:)= ZMASTOTN (:,:)/PSNOWDZN(:,:)
!
! 4. Vanishing or very thin snowpack check:
! -----------------------------------------
!
! NOTE: ONLY for very shallow snowpacks, mix properties (homogeneous):
! this avoids problems related to heat and mass exchange for
! thin layers during heavy snowfall or signifigant melt: one
! new/old layer can exceed the thickness of several old/new layers.
! Therefore, mix (conservative):
!
ZSUMHEAT(:)       = 0.0
ZSUMSWE(:)        = 0.0
ZSUMAGE(:)        = 0.0
ZSNOWMIX_DELTA(:) = 0.0
!
DO JL=1,INLVLS
   DO JI=1,INI
      IF(PSNOW(JI) < XSNOWCRITD)THEN
         ZSUMHEAT      (JI) = ZSUMHEAT(JI) + PSNOWHEAT(JI,JL)
         ZSUMSWE       (JI) = ZSUMSWE (JI) + PSNOWRHO (JI,JL)*PSNOWDZ(JI,JL)
         ZSUMAGE       (JI) = ZSUMAGE (JI) + PSNOWAGE (JI,JL)
         ZSNOWMIX_DELTA(JI) = 1.0
      ENDIF
   ENDDO
ENDDO
!
! Heat and mass are evenly distributed vertically:
! heat and mass (density and thickness) are constant
! in profile:
!
DO JL=1,INLVLS
   DO JI=1,INI
!
      ZSNOWHEATN(JI,JL) = ZSNOWMIX_DELTA(JI)*(ZSUMHEAT(JI)/INLVLS)  + &
                         (1.0-ZSNOWMIX_DELTA(JI))*ZSNOWHEATN(JI,JL)  
!
      PSNOWDZN(JI,JL)   = ZSNOWMIX_DELTA(JI)*(PSNOW(JI)/INLVLS)     + &
                        (1.0-ZSNOWMIX_DELTA(JI))*PSNOWDZN(JI,JL)  
!
      ZSNOWRHON(JI,JL)  = ZSNOWMIX_DELTA(JI)*(ZSUMSWE(JI)/PSNOW(JI)) + &
                        (1.0-ZSNOWMIX_DELTA(JI))*ZSNOWRHON(JI,JL)
!
      ZSNOWAGEN(JI,JL)  = ZSNOWMIX_DELTA(JI)*(ZSUMAGE(JI)/INLVLS)  + &
                         (1.0-ZSNOWMIX_DELTA(JI))*ZSNOWAGEN(JI,JL)  
!
   ENDDO
ENDDO
!
! 5. Update mass (density and thickness) and heat:
! ------------------------------------------------
!
PSNOWDZ  (:,:) = PSNOWDZN  (:,:)
PSNOWRHO (:,:) = ZSNOWRHON (:,:)
PSNOWHEAT(:,:) = ZSNOWHEATN(:,:)
PSNOWAGE (:,:) = ZSNOWAGEN (:,:)
!
IF (LHOOK) CALL DR_HOOK('SNOW3LTRANSF',1,ZHOOK_HANDLE)
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOW3LTRANSF
!####################################################################
!####################################################################
!####################################################################
SUBROUTINE SNOW3LALB(OMEB,PALBEDOSC,PSPECTRALALBEDO,PSNOWRHO,PSNOWAGE,PPERMSNOWFRAC,PPS )  
!
!!    PURPOSE
!!    -------
!     Calculate the snow surface albedo. Use the method of
!     CROCUS with 3 spectral albedo depending on snow density 
!     and age
!
!
USE MODD_SNOW_PAR, ONLY : XVAGING_GLACIER, XVAGING_NOGLACIER,  &
                          XVALB2,XVALB3,XVALB4,XVALB5,XVALB6,  &
                          XVALB7,XVALB8,XVALB9,XVALB10,XVALB11,&
                          XVDIOP1,XVRPRE1,XVRPRE2,XVPRES1,     &
                          XVW1,XVW2,XVSPEC1,XVSPEC2,XVSPEC3
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
LOGICAL,            INTENT(IN)      :: OMEB
!
REAL, DIMENSION(:), INTENT(IN)      :: PSNOWRHO
REAL, DIMENSION(:), INTENT(IN)      :: PSNOWAGE
REAL, DIMENSION(:), INTENT(IN)      :: PPERMSNOWFRAC
REAL, DIMENSION(:), INTENT(IN)      :: PPS
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PALBEDOSC
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSPECTRALALBEDO
!
!
!*      0.2    declarations of local variables
!
REAL, PARAMETER                 :: ZALBIR = 0.3
REAL, PARAMETER                 :: ZALBUV = 0.0
!
REAL, DIMENSION(SIZE(PSNOWRHO)) :: ZNVAGE, ZDIAM
!
REAL, DIMENSION(SIZE(PSNOWRHO)) :: ZALB1, ZALB2, ZALB3
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! 0. Initialize:
! ------------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3LALB',0,ZHOOK_HANDLE)
!
ZNVAGE(:)=XVAGING_GLACIER * PPERMSNOWFRAC(:) + XVAGING_NOGLACIER * (1.0-PPERMSNOWFRAC(:))
!
! 1. Snow grain size :
! --------------------
!
ZDIAM(:) = MIN(ZDSGRAIN_MAX, ZSNOW_AGRAIN + ZSNOW_BGRAIN*(PSNOWRHO(:)**4))
!
! 3. spectral albedo over 3 bands :
! ---------------------------------
!
ZALB1(:)=MIN(XVALB4,XVALB2-XVALB3*SQRT(ZDIAM(:)))
ZALB1(:)=MAX(XVALB11,ZALB1(:)-MIN(MAX(PPS(:)/XVPRES1,XVRPRE1),XVRPRE2)*XVALB10*PSNOWAGE(:)/ZNVAGE(:))
!
ZALB2(:)=MAX(ZALBIR,XVALB5-XVALB6*SQRT(ZDIAM(:)))
!
ZDIAM(:)=MIN(XVDIOP1,ZDIAM(:))
ZALB3(:)=MAX(ZALBUV,XVALB7*ZDIAM(:)-XVALB8*SQRT(ZDIAM(:))+XVALB9)
!
PSPECTRALALBEDO(:,1)=ZALB1(:)
PSPECTRALALBEDO(:,2)=ZALB2(:)
PSPECTRALALBEDO(:,3)=ZALB3(:)
!
! 4. total albedo :
! -----------------
!
PALBEDOSC(:)=XVSPEC1*ZALB1(:)+XVSPEC2*ZALB2(:)+XVSPEC3*ZALB3(:) 
!
!
IF (LHOOK) CALL DR_HOOK('SNOW3LALB',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOW3LALB
!####################################################################
!####################################################################
!####################################################################
      SUBROUTINE SNOW3LRAD(OMEB, PSNOWDZMIN, PSW_RAD, PSNOWALB,      &
                           PSPECTRALALBEDO, PSNOWDZ, PSNOWRHO, PALB, &
                           PPERMSNOWFRAC, PZENITH,  PSWNETSNOW,      &
                           PSWNETSNOWS, PRADSINK, PRADXS )  
!
!!    PURPOSE
!!    -------
!     Calculate the transmission of shortwave (solar) radiation
!     through the snowpack (using a form of Beer's Law: exponential
!     decay of radiation with increasing snow depth).
!
USE MODD_SNOW_PAR, ONLY : XVSPEC1,XVSPEC2,XVSPEC3
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
LOGICAL,            INTENT(IN)      :: OMEB ! if=T, then uppermost abs is diagnosed
!                                           !       since fluxes known
!
REAL,               INTENT(IN)      :: PSNOWDZMIN
!
REAL, DIMENSION(:), INTENT(IN)      :: PSW_RAD
REAL, DIMENSION(:), INTENT(IN)      :: PSNOWALB
REAL, DIMENSION(:), INTENT(IN)      :: PALB
REAL, DIMENSION(:), INTENT(IN)      :: PPERMSNOWFRAC
REAL, DIMENSION(:), INTENT(IN)      :: PZENITH
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWRHO, PSNOWDZ
REAL, DIMENSION(:,:), INTENT(IN)    :: PSPECTRALALBEDO
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PSWNETSNOW, PSWNETSNOWS
!
REAL, DIMENSION(:), INTENT(OUT)     :: PRADXS
!
REAL, DIMENSION(:,:), INTENT(OUT)   :: PRADSINK
!
!
!*      0.2    declarations of local variables
!
INTEGER                              :: JJ, JI
!
INTEGER                              :: INI
INTEGER                              :: INLVLS
!
REAL, DIMENSION(SIZE(PSNOWRHO,1))    :: ZRADTOT, ZPROJLAT
REAL, DIMENSION(SIZE(PSNOWRHO,1))    :: ZOPTICALPATH1, ZOPTICALPATH2, ZOPTICALPATH3
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZDSGRAIN, ZCOEF, ZSNOWDZ
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZBETA1, ZBETA2, ZBETA3
!
!
! ISBA-ES Radiation extinction coefficients: (see Loth and Graf 1993):
! see Boone, Meteo-France/CNRM Note de Centre No. 70 (2002)
!
REAL, PARAMETER                      :: ZSNOWRAD_CVEXT  = 3.8e-3   ! [(m5/2)/kg]
!
! calibration coefficients for exctinction computation from CROCUS
REAL, PARAMETER :: ZPVBETA3=1.098E-2
REAL, PARAMETER :: ZPVBETA5=2000.
!
REAL, PARAMETER :: ZMINCOSZEN      = 0.01
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
! 0. Initialization:
! ------------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3LRAD',0,ZHOOK_HANDLE)
!
INI    = SIZE(PSNOWDZ(:,:),1)
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
ZPROJLAT(:)=(1.0-PPERMSNOWFRAC(:))+PPERMSNOWFRAC(:)/MAX(ZMINCOSZEN,COS(PZENITH(:)))
!
! Snow grain size:
!
ZDSGRAIN(:,:) = MIN(ZDSGRAIN_MAX, ZSNOW_AGRAIN + ZSNOW_BGRAIN*(PSNOWRHO(:,:)**4))
!
! extinction coefficient:
!
ZBETA1(:,:)=ZSNOWRAD_CVEXT*PSNOWRHO(:,:)/SQRT(ZDSGRAIN(:,:))
ZBETA2(:,:)=ZPVBETA3*PSNOWRHO(:,:)/SQRT(ZDSGRAIN(:,:))
ZBETA3(:,:)=ZPVBETA5
!
ZOPTICALPATH1(:) = 0.0
ZOPTICALPATH2(:) = 0.0
ZOPTICALPATH3(:) = 0.0
!
DO JJ=1,INLVLS
   DO JI=1,INI
      !
      ZOPTICALPATH1(JI   ) = ZOPTICALPATH1(JI) + ZBETA1(JI,JJ)*ZSNOWDZ(JI,JJ)
      ZOPTICALPATH2(JI   ) = ZOPTICALPATH2(JI) + ZBETA2(JI,JJ)*ZSNOWDZ(JI,JJ)
      ZOPTICALPATH3(JI   ) = ZOPTICALPATH3(JI) + ZBETA3(JI,JJ)*ZSNOWDZ(JI,JJ)
      !
      ZCOEF (JI,JJ) = XVSPEC1*(1.0-PSPECTRALALBEDO(JI,1))*EXP(-ZOPTICALPATH1(JI)*ZPROJLAT(JI)) &
                    + XVSPEC2*(1.0-PSPECTRALALBEDO(JI,2))*EXP(-ZOPTICALPATH2(JI)*ZPROJLAT(JI)) &
                    + XVSPEC3*(1.0-PSPECTRALALBEDO(JI,3))*EXP(-ZOPTICALPATH3(JI)*ZPROJLAT(JI))
      !
   ENDDO
ENDDO
!
! 3. Radiation at each level: (W/m2)
! ----------------------------------
!
IF(OMEB)THEN

   ZCOEF(:,1)          = 1.0 - PSWNETSNOWS(:)/MAX(1.E-4,PSWNETSNOW(:))

ELSE

   PSWNETSNOW(:)       = PSW_RAD(:)*(1.-PSNOWALB(:))
   PSWNETSNOWS(:)      = PSWNETSNOW(:)*(1.0-ZCOEF(:,1)) ! SW rad absorbed in uppermost layer

ENDIF
!
! NOTE: as this is a "sink" term for heat, it
! as assigned a negative value.
DO JJ=1,INLVLS
   DO JI=1,INI
      PRADSINK(JI,JJ)  = -PSW_RAD(JI)*ZCOEF(JI,JJ)
   ENDDO
ENDDO
!
! For thin snow packs, radiation might reach base of
! snowpack...so we influence this amount with sfc albedo
! and (outside of this routine) add any excess heat
! to underlying soil:
!
PRADSINK(:,INLVLS) = PRADSINK(:,INLVLS)*(1.0-PALB(:))
!
! 4. Excess radiation not absorbed by snowpack (W/m2):
! ----------------------------------------------------
!
ZRADTOT(:)    = PRADSINK(:,1)+(1.-PSNOWALB(:))*PSW_RAD(:)
DO JJ=2,INLVLS
   DO JI=1,INI
      ZRADTOT(JI) = ZRADTOT(JI) + PRADSINK(JI,JJ)-PRADSINK(JI,JJ-1)
   ENDDO
ENDDO
!
PRADXS(:)     = (1.-PSNOWALB(:))*PSW_RAD(:) - ZRADTOT(:)
!
IF (LHOOK) CALL DR_HOOK('SNOW3LRAD',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOW3LRAD
!####################################################################
!####################################################################
!####################################################################
      SUBROUTINE SNOW3LTHRM(PSNOWRHO,PSCOND,PSNOWTEMP,PPS,PSNOWLIQ)
!
!!    PURPOSE
!!    -------
!     Calculate snow thermal conductivity from
!     Sun et al. 1999, J. of Geophys. Res., 104, 19587-19579
!     (vapor) and Anderson, 1976, NOAA Tech. Rep. NWS 19 (snow).
!
!
USE MODD_CSTS,ONLY : XP00, XCONDI, XRHOLW
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:), INTENT(IN)      :: PPS
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWTEMP, PSNOWRHO, PSNOWLIQ
!
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSCOND
!
!
!*      0.2    declarations of local variables
!
INTEGER                              :: JJ, JI
!
INTEGER                              :: INI
INTEGER                              :: INLVLS
!
! ISBA-ES Thermal conductivity coefficients from Anderson (1976):
! see Boone, Meteo-France/CNRM Note de Centre No. 70 (2002)
!
REAL, PARAMETER                      :: ZSNOWTHRMCOND1 = 0.02    ! [W/(m K)]
REAL, PARAMETER                      :: ZSNOWTHRMCOND2 = 2.5E-6  ! [W m5/(kg2 K)]
!
! ISBA-ES Thermal conductivity coefficients from Yen (1981):
!
REAL, PARAMETER                      :: ZSNOWTHRMEXP = 1.88
!
! ISBA-ES Thermal conductivity: Implicit vapor diffn effects
! (sig only for new snow OR high altitudes)
! from Sun et al. (1999): based on data from Jordan (1991)
! see Boone, Meteo-France/CNRM Note de Centre No. 70 (2002)
!
REAL, PARAMETER                      :: ZSNOWTHRMCOND_AVAP  = -0.06023 ! [W/(m K)]
REAL, PARAMETER                      :: ZSNOWTHRMCOND_BVAP  = -2.5425  ! (W/m)
REAL, PARAMETER                      :: ZSNOWTHRMCOND_CVAP  = -289.99  ! (K)
!
CHARACTER(LEN=5)                     :: YSNOWCOND !should be in namelist
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3LTHRM',0,ZHOOK_HANDLE)
!
INI    = SIZE(PSNOWRHO(:,:),1)
INLVLS = SIZE(PSNOWRHO(:,:),2)
!
! 1. Snow thermal conductivity
! ----------------------------
!
YSNOWCOND='YEN81' !should be in namelist
!
IF(YSNOWCOND=='AND76')THEN
!  Thermal conductivity coefficients from Anderson (1976)
  PSCOND(:,:) = (ZSNOWTHRMCOND1 + ZSNOWTHRMCOND2*PSNOWRHO(:,:)*PSNOWRHO(:,:))
ELSE
! Thermal conductivity coefficients from Yen (1981)
  PSCOND(:,:) = XCONDI * EXP(ZSNOWTHRMEXP * LOG(PSNOWRHO(:,:)/XRHOLW))
ENDIF
!
! 2. Implicit vapor diffn effects
! -------------------------------
!
DO JJ=1,INLVLS
   DO JI=1,INI
    PSCOND(JI,JJ) = PSCOND(JI,JJ) + MAX(0.0,(ZSNOWTHRMCOND_AVAP+(ZSNOWTHRMCOND_BVAP/(PSNOWTEMP(JI,JJ) &
                                  + ZSNOWTHRMCOND_CVAP)))*(XP00/PPS(JI)))
   ENDDO
ENDDO
!
IF (LHOOK) CALL DR_HOOK('SNOW3LTHRM',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOW3LTHRM
!####################################################################
!####################################################################
!####################################################################
      SUBROUTINE SNOW3LEBUD(HSNOWRES, HIMPLICIT_WIND,                                     &
                              PPEW_A_COEF, PPEW_B_COEF,                                   &
                              PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,         &
                              PSNOWDZMIN,                                                 &
                              PZREF,PTS,PSNOWRHO,PSNOWLIQ,PSCAP,PSCOND1,PSCOND2,          &
                              PUREF,PEXNS,PEXNA,PDIRCOSZW,PVMOD,                          &
                              PLW_RAD,PSW_RAD,PTA,PQA,PPS,PTSTEP,                         &
                              PSNOWDZ1,PSNOWDZ2,PALBT,PZ0,PZ0EFF,PZ0H,                    &
                              PSFCFRZ,PRADSINK,PHPSNOW,                                   &
                              PCT,PEMIST,PRHOA,PTSTERM1,PTSTERM2,PRA,PCDSNOW,PCHSNOW,     &
                              PQSAT,PDQSAT,PRSRA,PUSTAR2_IC,PRI,                            &
                              PPET_A_COEF_T,PPEQ_A_COEF_T,PPET_B_COEF_T,PPEQ_B_COEF_T     )  
!
!!    PURPOSE
!!    -------
!     Calculate surface energy budget linearization (terms) and turbulent
!     exchange coefficients/resistance between surface and atmosphere.
!     (Noilhan and Planton 1989; Giordani 1993; Noilhan and Mahfouf 1996)
!
USE MODD_CSTS,     ONLY : XCPD, XRHOLW, XSTEFAN, XLVTT, XLSTT
USE MODD_SNOW_PAR, ONLY : X_RI_MAX, XEMISSN
!
USE MODE_THERMOS
!
USE MODI_SURFACE_RI
USE MODI_SURFACE_AERO_COND
USE MODI_SURFACE_CD
!
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP, PSNOWDZMIN
!
 CHARACTER(LEN=*),     INTENT(IN)    :: HSNOWRES ! type of sfc resistance
!                                      DEFAULT=Louis (1979), standard ISBA
!                                      method. Option to limit Ri number
!                                      for very srtable conditions
!
 CHARACTER(LEN=*),     INTENT(IN)  :: HIMPLICIT_WIND   ! wind implicitation option
!                                                     ! 'OLD' = direct
!                                                     ! 'NEW' = Taylor serie, order 1
!
REAL, DIMENSION(:), INTENT(IN)      :: PPEW_A_COEF, PPEW_B_COEF,                   &
                                         PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF,      &
                                         PPEQ_B_COEF  
!                                      PPEW_A_COEF = wind coefficient (m2s/kg)
!                                      PPEW_B_COEF = wind coefficient (m/s)
!                                      PPET_A_COEF = A-air temperature coefficient
!                                      PPET_B_COEF = B-air temperature coefficient
!                                      PPEQ_A_COEF = A-air specific humidity coefficient
!                                      PPEQ_B_COEF = B-air specific humidity coefficient
!
REAL, DIMENSION(:), INTENT(IN)    :: PZREF, PTS, PSNOWDZ1, PSNOWDZ2,          &
                                         PRADSINK, PSNOWRHO, PSNOWLIQ, PSCAP,   &
                                         PSCOND1, PSCOND2,                      &
                                         PZ0, PHPSNOW,                          &
                                         PALBT, PZ0EFF, PZ0H  
!
REAL, DIMENSION(:), INTENT(IN)    :: PSW_RAD, PLW_RAD, PTA, PQA, PPS, PRHOA
!
REAL, DIMENSION(:), INTENT(IN)    :: PUREF, PEXNS, PEXNA, PDIRCOSZW, PVMOD
!
REAL, DIMENSION(:), INTENT(OUT)   :: PTSTERM1, PTSTERM2, PEMIST, PRA,         &
                                       PCT, PSFCFRZ, PCDSNOW, PCHSNOW,          &
                                       PQSAT, PDQSAT, PRSRA  
!
REAL, DIMENSION(:), INTENT(OUT)   :: PUSTAR2_IC,                        &
                                       PPET_A_COEF_T, PPEQ_A_COEF_T,    &
                                       PPET_B_COEF_T, PPEQ_B_COEF_T 
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRI
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PTS))        :: ZAC, ZRI,                          &
                                       ZSCONDA, ZA, ZB, ZC,             &
                                       ZCDN, ZSNOWDZM1, ZSNOWDZM2,      &
                                       ZVMOD, ZUSTAR2, ZTS3, ZLVT,      &
                                       Z_CCOEF  
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! 1. New saturated specific humidity and derrivative:
! ---------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3LEBUD',0,ZHOOK_HANDLE)
!
ZRI   (:) = XUNDEF
!
PQSAT (:) =  QSATI(PTS(:),PPS(:))
PDQSAT(:) = DQSATI(PTS(:),PPS(:),PQSAT(:))
!
!
! 2. Surface properties:
! ----------------------
! It might be of interest to use snow-specific roughness
! or a temperature dependence on emissivity:
! but for now, use ISBA defaults.
!
PEMIST(:) = XEMISSN
!
! 2. Computation of resistance and drag coefficient
! -------------------------------------------------
!
 CALL SURFACE_RI(PTS, PQSAT, PEXNS, PEXNA, PTA, PQA,                  &
                PZREF, PUREF, PDIRCOSZW, PVMOD, ZRI                  )  
!
! Simple adaptation of method by Martin and Lejeune (1998)
! to apply a lower limit to turbulent transfer coef
! by defining a maximum Richarson number for stable
! conditions:
!
IF(HSNOWRES=='RIL') THEN
  DO JJ = 1,SIZE(ZRI)
    ZRI(JJ) = MIN(X_RI_MAX,ZRI(JJ))
  ENDDO
ENDIF
!
PRI(:)=ZRI(:)
!
! Surface aerodynamic resistance for heat transfers
!
 CALL SURFACE_AERO_COND(ZRI, PZREF, PUREF, PVMOD, PZ0, PZ0H, ZAC, PRA, PCHSNOW)
!
! For atmospheric model coupling:
!
 CALL SURFACE_CD(ZRI, PZREF, PUREF, PZ0EFF, PZ0H, PCDSNOW, ZCDN)
!
PRSRA(:) = PRHOA(:) / PRA(:)
!
!
! Modify flux-form implicit coupling coefficients:
! - wind components:
!
IF(HIMPLICIT_WIND=='OLD')THEN
! old implicitation
  ZUSTAR2(:) = ( PCDSNOW(:)*PVMOD(:)*PPEW_B_COEF(:)) /        &
               (1.0-PRHOA(:)*PCDSNOW(:)*PVMOD(:)*PPEW_A_COEF(:))  
ELSE
! new implicitation
  ZUSTAR2(:) = (PCDSNOW(:)*PVMOD(:)*(2.*PPEW_B_COEF(:)-PVMOD(:)))  &
               / (1.0-2.0*PRHOA(:)*PCDSNOW(:)*PVMOD(:)*PPEW_A_COEF(:))
ENDIF               
!
ZVMOD(:)       = PRHOA(:)*PPEW_A_COEF(:)*ZUSTAR2(:) + PPEW_B_COEF(:)
ZVMOD(:)       = MAX(ZVMOD(:),0.)
!
WHERE(PPEW_A_COEF(:)/= 0.)
      ZUSTAR2(:) = MAX( ( ZVMOD(:) - PPEW_B_COEF(:) ) / (PRHOA(:)*PPEW_A_COEF(:)), 0.)
ENDWHERE
!
! implicit wind friction
ZUSTAR2(:) = MAX(ZUSTAR2(:),0.)
!
PUSTAR2_IC(:) =  ZUSTAR2(:)
!
! 3. Calculate linearized surface energy budget components:
! ---------------------------------------------------------
! To prevent numerical difficulties for very thin snow
! layers, limit the grid "thinness": this is important as
! layers become vanishing thin:
!
ZSNOWDZM1(:) = MAX(PSNOWDZ1(:), PSNOWDZMIN)
ZSNOWDZM2(:) = MAX(PSNOWDZ2(:), PSNOWDZMIN)
!
! Surface thermal inertia:
!
PCT(:)      = 1.0/(PSCAP(:)*ZSNOWDZM1(:))
!
! Fraction of surface frozen (sublimation) with the remaining
! fraction being liquid (evaporation):
!
PSFCFRZ(:)  = 1.0 - PSNOWLIQ(:)*XRHOLW/(ZSNOWDZM1(:)*PSNOWRHO(:))
!
! Thermal conductivity between uppermost and lower snow layers:
!
ZSCONDA(:)  = (ZSNOWDZM1(:)*PSCOND1(:) + ZSNOWDZM2(:)*PSCOND2(:))/           &
                (ZSNOWDZM1(:)+ZSNOWDZM2(:))  
!
! Transform implicit coupling coefficients: 
! Note, surface humidity is 100% over snow.
!
! - specific humidity:
!
Z_CCOEF(:)       = 1.0 - PPEQ_A_COEF(:)*PRSRA(:)
!
PPEQ_A_COEF_T(:) = - PPEQ_A_COEF(:)*PRSRA(:)*PDQSAT(:)/Z_CCOEF(:)
!
PPEQ_B_COEF_T(:) = ( PPEQ_B_COEF(:) - PPEQ_A_COEF(:)*PRSRA(:)*(PQSAT(:) - &
                       PDQSAT(:)*PTS(:)) )/Z_CCOEF(:)  
!
! - air temperature:
!   (assumes A and B correspond to potential T):
!
Z_CCOEF(:)       = (1.0 - PPET_A_COEF(:)*PRSRA(:))/PEXNA(:)
!
PPET_A_COEF_T(:) = - PPET_A_COEF(:)*PRSRA(:)/(PEXNS(:)*Z_CCOEF(:))
!
PPET_B_COEF_T(:) = PPET_B_COEF(:)/Z_CCOEF(:)
!
!
! Energy budget solution terms:
!
ZTS3(:) = PEMIST(:) * XSTEFAN * PTS(:)**3
ZLVT(:) = (1.-PSFCFRZ(:))*XLVTT + PSFCFRZ(:)*XLSTT
!
ZA(:)   = 1. / PTSTEP + PCT(:) * (4. * ZTS3(:) +                               &
            PRSRA(:) *  ZLVT(:) * (PDQSAT(:) - PPEQ_A_COEF_T(:))                 &
            + PRSRA(:) * XCPD * ( (1./PEXNS(:))-(PPET_A_COEF_T(:)/PEXNA(:)) )    &
            + (2*ZSCONDA(:)/(ZSNOWDZM2(:)+ZSNOWDZM1(:))) )  
!
ZB(:)   = 1. / PTSTEP + PCT(:) * (3. * ZTS3(:) +                               &
            PRSRA(:) * PDQSAT(:) *  ZLVT(:) )  
!
ZC(:)   = PCT(:) * (PRSRA(:) * XCPD * PPET_B_COEF_T(:)/PEXNA(:) + PSW_RAD(:) * &
            (1. - PALBT(:)) + PEMIST(:)*PLW_RAD(:) - PRSRA(:) *                  &
            ZLVT(:) *  (PQSAT(:)-PPEQ_B_COEF_T(:))                               &
            + PHPSNOW(:) + PRADSINK(:) )  
!
!
! Coefficients needed for implicit solution
! of linearized surface energy budget:
!
PTSTERM2(:) = 2*ZSCONDA(:)*PCT(:)/(ZA(:)*(ZSNOWDZM2(:)+ZSNOWDZM1(:)))
!
PTSTERM1(:) = (PTS(:)*ZB(:) + ZC(:))/ZA(:)
IF (LHOOK) CALL DR_HOOK('SNOW3LEBUD',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOW3LEBUD
!####################################################################
!####################################################################
!####################################################################
      SUBROUTINE SNOW3LSOLVT(OMEB,PTSTEP,PSNOWDZMIN,                  &
                               PSNOWDZ,PSCOND,PSCAP,PTG,              &
                               PSOILCOND,PD_G,                        &
                               PRADSINK,PCT,PTERM1,PTERM2,            &
                               PPET_A_COEF_T,PPEQ_A_COEF_T,           &
                               PPET_B_COEF_T,PPEQ_B_COEF_T,           &
                               PTA_IC, PQA_IC,                        &
                               PGRNDFLUX,PGRNDFLUXO,PSNOWTEMP,        &
                               PSNOWFLUX                              )  
!
!!    PURPOSE
!!    -------
!     This subroutine solves the 1-d diffusion of 'ZSNOWTEMP' using a
!     layer averaged set of equations which are time differenced
!     using the backward-difference scheme (implicit).
!     The eqs are solved rapidly by taking advantage of the
!     fact that the matrix is tridiagonal. This is a very
!     general routine and can be used for the 1-d diffusion of any
!     quantity as long as the diffusity is not a function of the
!     quantity being diffused. Aaron Boone 8-98. Soln to the eq:
!
!                 c  dQ    d  k dQ    dS
!                    -- =  --   -- -  --
!                    dt    dx   dx    dx
!
!     where k = k(x) (thermal conductivity), c = c(x) (heat capacity)
!     as an eg. for temperature/heat eq. S is a sink (-source) term.
!     Diffusivity is k/c
!
!
USE MODD_CSTS,ONLY : XTT
!
USE MODI_TRIDIAG_GROUND
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
LOGICAL,            INTENT(IN)      :: OMEB
!
REAL,               INTENT(IN)      :: PTSTEP, PSNOWDZMIN
!
REAL, DIMENSION(:), INTENT(IN)      :: PTG, PSOILCOND, PD_G,         &
                                         PCT, PTERM1, PTERM2  

!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWDZ, PSCOND, PSCAP,       &
                                         PRADSINK  
!
REAL, DIMENSION(:), INTENT(IN)      :: PPET_A_COEF_T, PPEQ_A_COEF_T, &
                                         PPET_B_COEF_T, PPEQ_B_COEF_T  
!                                       
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWTEMP
!
REAL, DIMENSION(:), INTENT(OUT)     :: PGRNDFLUX, PGRNDFLUXO, PSNOWFLUX,     &
                                         PTA_IC, PQA_IC   
!
!
!*      0.2    declarations of local variables
!
!
INTEGER                             :: JJ, JI
!
INTEGER                             :: INI
INTEGER                             :: INLVLS
!
REAL, DIMENSION(SIZE(PTG))                     :: ZSNOWTEMP_DELTA
!
REAL, DIMENSION(SIZE(PSNOWDZ,1),SIZE(PSNOWDZ,2)) :: ZSNOWTEMP, ZDTERM, ZCTERM, &
                                                      ZFRCV, ZAMTRX, ZBMTRX,     &
                                                      ZCMTRX  
!
REAL, DIMENSION(SIZE(PSNOWDZ,1),SIZE(PSNOWDZ,2)) :: ZWORK1, ZWORK2, ZDZDIF,    &
                                                      ZSNOWDZM  
!
REAL, DIMENSION(SIZE(PSNOWDZ,1),SIZE(PSNOWDZ,2)-1) :: ZSNOWTEMP_M,             &
                                                      ZFRCV_M, ZAMTRX_M,         &
                                                      ZBMTRX_M, ZCMTRX_M  
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
! 0. Initialize:
! ------------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3LSOLVT',0,ZHOOK_HANDLE)
ZSNOWTEMP(:,:)  = PSNOWTEMP(:,:)
INI             = SIZE(PSNOWDZ(:,:),1)
INLVLS          = SIZE(PSNOWDZ(:,:),2)
!
!
! 1. Calculate tri-diagnoal matrix coefficients:
! ----------------------------------------------
! For heat transfer, assume a minimum grid
! thickness (to prevent numerical
! problems for very thin snow cover):
!
ZSNOWDZM(:,:)  = MAX(PSNOWDZ(:,:), PSNOWDZMIN)
!
ZWORK1(:,:)      = ZSNOWDZM(:,:)*PSCOND(:,:)
DO JJ=1,INLVLS-1
   DO JI=1,INI
      ZDZDIF(JI,JJ)  = ZSNOWDZM(JI,JJ) + ZSNOWDZM(JI,JJ+1)
      ZWORK2(JI,JJ)  = ZSNOWDZM(JI,JJ+1)*PSCOND(JI,JJ+1)
   ENDDO
ENDDO
ZDZDIF(:,INLVLS) = ZSNOWDZM(:,INLVLS) + PD_G(:)
ZWORK2(:,INLVLS) = PD_G(:)*PSOILCOND(:)
!
ZDTERM(:,:)      = 2.0*(ZWORK1(:,:)+ZWORK2(:,:))/(ZDZDIF(:,:)*ZDZDIF(:,:))
!
ZCTERM(:,:)      = PSCAP(:,:)*ZSNOWDZM(:,:)/PTSTEP
!
! 2. Set up tri-diagonal matrix
! -----------------------------
!
! Upper BC
!
ZAMTRX(:,1) =  0.0
ZBMTRX(:,1) =  1./(PCT(:)*PTSTEP)
ZCMTRX(:,1) = -PTERM2(:)*ZBMTRX(:,1)
ZFRCV(:,1)  =  PTERM1(:)*ZBMTRX(:,1)
!
!
! Interior Grid
!
DO JJ=2,INLVLS-1
   DO JI=1,INI
      ZAMTRX(JI,JJ) = -ZDTERM(JI,JJ-1)
      ZBMTRX(JI,JJ) =  ZCTERM(JI,JJ) + ZDTERM(JI,JJ-1) + ZDTERM(JI,JJ)
      ZCMTRX(JI,JJ) = -ZDTERM(JI,JJ)
      ZFRCV (JI,JJ) =  ZCTERM(JI,JJ)*PSNOWTEMP(JI,JJ) - (PRADSINK(JI,JJ-1)-PRADSINK(JI,JJ))  
   ENDDO
ENDDO
!
!Lower BC
!
ZAMTRX(:,INLVLS) = -ZDTERM(:,INLVLS-1)
ZBMTRX(:,INLVLS) =  ZCTERM(:,INLVLS) + ZDTERM(:,INLVLS-1) +                   &
                      ZDTERM(:,INLVLS)  
ZCMTRX(:,INLVLS) =  0.0
ZFRCV(:,INLVLS)  =  ZCTERM(:,INLVLS)*PSNOWTEMP(:,INLVLS) +                    &
                      ZDTERM(:,INLVLS)*PTG(:)                                   &
                      - (PRADSINK(:,INLVLS-1)-PRADSINK(:,INLVLS))  
!
! - - -------------------------------------------------
!
! 4. Compute solution vector
! --------------------------
!
 CALL TRIDIAG_GROUND(ZAMTRX,ZBMTRX,ZCMTRX,ZFRCV,ZSNOWTEMP)
!
! Heat flux between surface and 2nd snow layers: (W/m2)
!
PSNOWFLUX(:)      = ZDTERM(:,1)*(ZSNOWTEMP(:,1) - ZSNOWTEMP(:,2))
!
!
! 5. Snow melt case
! -----------------
! If melting in uppermost layer, assume surface layer
! temperature at freezing point and re-evaluate lower
! snowpack temperatures. This is done as it is most likely
! most signigant melting will occur within a time step in surface layer.
! Surface energy budget (and fluxes) will
! be re-calculated (outside of this routine):
!
IF(.NOT.OMEB)THEN
!
   ZAMTRX_M(:,1)  =  0.0
   ZBMTRX_M(:,1)  =  ZCTERM(:,2) + ZDTERM(:,1) + ZDTERM(:,2)
   ZCMTRX_M(:,1)  = -ZDTERM(:,2)
   ZFRCV_M(:,1)   =  ZCTERM(:,2)*PSNOWTEMP(:,2) + XTT*ZDTERM(:,1)  -  &
                    (PRADSINK(:,1)-PRADSINK(:,2))  
!
   DO JJ=2,INLVLS-1
      DO JI=1,INI
         ZAMTRX_M   (JI,JJ) = ZAMTRX   (JI,JJ+1)
         ZBMTRX_M   (JI,JJ) = ZBMTRX   (JI,JJ+1)
         ZCMTRX_M   (JI,JJ) = ZCMTRX   (JI,JJ+1)
         ZFRCV_M    (JI,JJ) = ZFRCV    (JI,JJ+1)
         ZSNOWTEMP_M(JI,JJ) = PSNOWTEMP(JI,JJ+1)
      ENDDO
   ENDDO
!
   CALL TRIDIAG_GROUND(ZAMTRX_M,ZBMTRX_M,ZCMTRX_M,ZFRCV_M,ZSNOWTEMP_M)
!
! If melting for 2 consecuative time steps, then replace current T-profile
! with one assuming T=Tf in surface layer:
!
   ZSNOWTEMP_DELTA(:)    = 0.0
   !
   WHERE(ZSNOWTEMP(:,1) > XTT .AND. PSNOWTEMP(:,1) == XTT)
      PSNOWFLUX(:)       = ZDTERM(:,1)*(XTT - ZSNOWTEMP_M(:,1))
      ZSNOWTEMP_DELTA(:) = 1.0
   END WHERE
!
   DO JJ=2,INLVLS
      DO JI=1,INI
         ZSNOWTEMP(JI,JJ) = ZSNOWTEMP_DELTA(JI)*ZSNOWTEMP_M(JI,JJ-1)   &
              + (1.0-ZSNOWTEMP_DELTA(JI))*ZSNOWTEMP(JI,JJ)  
      ENDDO
   ENDDO
!
ENDIF
!
!
! 6. Lower boundary flux:
! -----------------------
! NOTE: evaluate this term assuming the snow layer
! can't exceed the freezing point as this adjustment
! is made in melting routine. Then must adjust temperature
! to conserve energy:
!
PGRNDFLUXO(:)       = ZDTERM(:,INLVLS)*(ZSNOWTEMP(:,INLVLS)         -PTG(:))
PGRNDFLUX(:)        = ZDTERM(:,INLVLS)*(MIN(XTT,ZSNOWTEMP(:,INLVLS))-PTG(:))
!
ZSNOWTEMP(:,INLVLS) = ZSNOWTEMP(:,INLVLS) + (PGRNDFLUXO(:)-PGRNDFLUX(:))/ZCTERM(:,INLVLS)
!
! 7. Update temperatute profile in time:
! --------------------------------------
!
PSNOWTEMP(:,:)      = ZSNOWTEMP(:,:)
!
!
! 8. Compute new (implicit) air T and specific humidity
! -----------------------------------------------------
!
IF(.NOT.OMEB)THEN
!
   PTA_IC(:) = PPET_B_COEF_T(:) + PPET_A_COEF_T(:)* PSNOWTEMP(:,1)
!
   PQA_IC(:) = PPEQ_B_COEF_T(:) + PPEQ_A_COEF_T(:)* PSNOWTEMP(:,1)
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('SNOW3LSOLVT',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE SNOW3LSOLVT
!####################################################################
!####################################################################
!####################################################################
      SUBROUTINE SNOW3LMELT(PTSTEP,PSCAP,PSNOWTEMP,PSNOWDZ,         &
                              PSNOWRHO,PSNOWLIQ,PMELTXS               )  
!
!
!!    PURPOSE
!!    -------
!     Calculate snow melt (resulting from surface fluxes, ground fluxes,
!     or internal shortwave radiation absorbtion). It is used to
!     augment liquid water content, maintain temperatures
!     at or below freezing, and possibly reduce the mass
!     or compact the layer(s).
!
!
USE MODD_CSTS,ONLY : XTT, XLMTT, XRHOLW, XRHOLI
!
USE MODE_SNOW3L
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSCAP
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWDZ, PSNOWTEMP, PSNOWRHO,   &
                                           PSNOWLIQ  
!
REAL, DIMENSION(:), INTENT(OUT)     :: PMELTXS
!
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZPHASE, ZCMPRSFACT,   &
                                                        ZSNOWLWE, ZWHOLDMAX,  &
                                                        ZSNOWMELT, ZSNOWTEMP, &
                                                        ZMELTXS  
!
INTEGER :: JWRK, JI ! loop counter
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
! 0. Initialize:
! ---------------------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3LMELT',0,ZHOOK_HANDLE)
ZPHASE(:,:)     = 0.0
ZCMPRSFACT(:,:) = 0.0
ZSNOWLWE(:,:)   = 0.0
ZWHOLDMAX(:,:)  = 0.0
ZSNOWMELT(:,:)  = 0.0
ZSNOWTEMP(:,:)  = 0.0
ZMELTXS(:,:)    = 0.0
!
! 1. Determine amount of melt in each layer:
! ------------------------------------------
!
WHERE(PSNOWDZ > 0.0)
!
! Total Liquid equivalent water content of snow (m):
!
   ZSNOWLWE(:,:) = PSNOWRHO(:,:)*PSNOWDZ(:,:)/XRHOLW
!
! Melt snow if excess energy and snow available:
! Phase change (J/m2)
!
   ZPHASE(:,:)  = MIN(PSCAP(:,:)*MAX(0.0, PSNOWTEMP(:,:) - XTT)*      &
                    PSNOWDZ(:,:),                                       &
                    MAX(0.0,ZSNOWLWE(:,:)-PSNOWLIQ(:,:))*XLMTT*XRHOLW)  
!
!
! Update snow liq water content and temperature if melting:
! liquid water available for next layer from melting of snow
! which is assumed to be leaving the current layer (m):
!
   ZSNOWMELT(:,:) = ZPHASE(:,:)/(XLMTT*XRHOLW)
!
! Cool off snow layer temperature due to melt:
!
   ZSNOWTEMP(:,:) = PSNOWTEMP(:,:) - ZPHASE(:,:)/(PSCAP(:,:)*PSNOWDZ(:,:))
!
   PSNOWTEMP(:,:) = MIN(XTT, ZSNOWTEMP(:,:))
!
   ZMELTXS(:,:)   = (ZSNOWTEMP(:,:)-PSNOWTEMP(:,:))*PSCAP(:,:)*PSNOWDZ(:,:)
!
! Loss of snowpack depth: (m) and liquid equiv (m):
! Compression factor for melt loss: this decreases
! layer thickness and increases density thereby leaving
! total SWE constant. NOTE: All melt water
! in excess of the holding capacity is NOT used
! for compression, rather it decreases the layer
! thickness ONLY, causing a loss of SWE (outside
! of this routine).
!
   ZWHOLDMAX(:,:)  = SNOW3LHOLD(PSNOWRHO,PSNOWDZ)
!
   ZCMPRSFACT(:,:) = (ZSNOWLWE(:,:)-MIN(PSNOWLIQ(:,:)+ZSNOWMELT(:,:),     &
                        ZWHOLDMAX(:,:)))/                                   &
                       (ZSNOWLWE(:,:)-MIN(PSNOWLIQ(:,:),ZWHOLDMAX(:,:)))  
!
   PSNOWDZ(:,:)    = PSNOWDZ(:,:)*ZCMPRSFACT(:,:)
   PSNOWRHO(:,:)   = ZSNOWLWE(:,:)*XRHOLW/PSNOWDZ(:,:)
!
! Make sure maximum density is not surpassed! If it is, lower the density
! and increase the snow thickness accordingly:

   ZCMPRSFACT(:,:) = MAX(XRHOLI, PSNOWRHO(:,:))/XRHOLI
   PSNOWDZ(:,:)    = PSNOWDZ(:,:)*ZCMPRSFACT(:,:)
   PSNOWRHO(:,:)   = ZSNOWLWE(:,:)*XRHOLW/PSNOWDZ(:,:)
!
!
! 2. Add snow melt to current snow liquid water content:
! ------------------------------------------------------
!
   PSNOWLIQ(:,:)   = PSNOWLIQ(:,:) + ZSNOWMELT(:,:)
!
END WHERE
!
! 3. Excess heat from melting
! ---------------------------
! use it to warm underlying ground/vegetation layer to conserve energy
!
PMELTXS(:) = 0.
DO JWRK = 1, SIZE(ZMELTXS,2)
   DO JI = 1, SIZE(ZMELTXS,1)
      PMELTXS(JI) = PMELTXS(JI) + ZMELTXS(JI,JWRK)
   ENDDO
ENDDO
PMELTXS(:) = PMELTXS(:) / PTSTEP   ! (W/m2)
IF (LHOOK) CALL DR_HOOK('SNOW3LMELT',1,ZHOOK_HANDLE)
!
!
!
END SUBROUTINE SNOW3LMELT
!####################################################################
!####################################################################
!####################################################################
      SUBROUTINE SNOW3LREFRZ(PTSTEP,PRR,                            &
                               PSNOWRHO,PSNOWTEMP,PSNOWDZ,PSNOWLIQ,   &
                               PTHRUFAL                               )  
!
!
!!    PURPOSE
!!    -------
!     Calculate any freezing/refreezing of liquid water in the snowpack.
!     Also, calculate liquid water transmission and snow runoff.
!     Water flow causes layer thickness reduction and can cause
!     rapid densification of a layer.
!
!
USE MODD_CSTS,     ONLY : XTT, XLMTT, XRHOLW
USE MODD_SNOW_PAR, ONLY : XSNOWDMIN
!
USE MODE_SNOW3L
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                      :: PTSTEP
!
REAL, DIMENSION(:), INTENT(IN)        :: PRR
!
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWDZ, PSNOWTEMP, PSNOWLIQ, PSNOWRHO
!
REAL, DIMENSION(:), INTENT(INOUT)     :: PTHRUFAL
!
!
!
!*      0.2    declarations of local variables
!
INTEGER                               :: JJ, JI
!
INTEGER                               :: INI
INTEGER                               :: INLVLS
!
REAL, DIMENSION(SIZE(PRR))            :: ZPCPXS, ZTOTWCAP, ZRAINFALL
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZPHASE, ZFLOWLIQ,    &
                                                        ZSNOWLIQ, ZSNOWRHO,  &
                                                        ZWHOLDMAX, ZSNOWDZ,  &
                                                        ZSNOWTEMP, ZSCAP  
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),0:SIZE(PSNOWRHO,2)) :: ZFLOWLIQT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! 0. Initialize:
! --------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3LREFRZ',0,ZHOOK_HANDLE)
ZSNOWRHO(:,:)  = PSNOWRHO(:,:)
ZSNOWLIQ(:,:)  = PSNOWLIQ(:,:)
ZSNOWTEMP(:,:) = PSNOWTEMP(:,:)
INI            = SIZE(PSNOWDZ(:,:),1)
INLVLS         = SIZE(PSNOWDZ(:,:),2)
!
!
! 1. Refreeze due to heat transfer
!    -----------------------------
! Freeze liquid water in any layers which cooled due
! to heat transfer.
!
ZSCAP(:,:)     = SNOW3LSCAP(ZSNOWRHO)
ZPHASE(:,:)    = MIN(ZSCAP(:,:)*                                          &
                   MAX(0.0, XTT - ZSNOWTEMP(:,:))*PSNOWDZ(:,:),             &
                   ZSNOWLIQ(:,:)*XLMTT*XRHOLW)  
!
!
! Warm layer and reduce liquid if freezing occurs:
!
ZSNOWDZ(:,:)   = MAX(XSNOWDMIN/INLVLS, PSNOWDZ(:,:))
!
PSNOWTEMP(:,:) = ZSNOWTEMP(:,:) + ZPHASE(:,:)/(ZSCAP(:,:)*ZSNOWDZ(:,:))
!
! Reduce liquid portion if freezing occurs:
!
ZSNOWLIQ(:,:)  = ZSNOWLIQ(:,:) - ( (PSNOWTEMP(:,:)-ZSNOWTEMP(:,:))*       &
                   ZSCAP(:,:)*ZSNOWDZ(:,:)/(XLMTT*XRHOLW) )  
!
!
! 2. Reduce thickness due to snowmelt in excess of holding capacity
!    --------------------------------------------------------------
! Any water in excess of the
! Maximum holding space for liquid water
! amount is drained into next layer down.
! Loss of water due to snowmelt causes a reduction
! in snow layer mass by a reduction in thickness.
!
ZWHOLDMAX(:,:) = SNOW3LHOLD(PSNOWRHO,PSNOWDZ)
!
ZFLOWLIQ(:,:)  = MAX(0.,ZSNOWLIQ(:,:)-ZWHOLDMAX(:,:))
!
ZSNOWLIQ(:,:)  = ZSNOWLIQ(:,:) - ZFLOWLIQ(:,:)
!
ZSNOWDZ(:,:)   = PSNOWDZ(:,:) - ZFLOWLIQ(:,:)*XRHOLW/ZSNOWRHO(:,:)
!
ZSNOWDZ(:,:)   = MAX(0.0, ZSNOWDZ(:,:))  ! to prevent possible very small
!                                          negative values (machine prescision
!                                          as snow vanishes
!
!
! 3. Liquid water flow: liquid precipitation and meltwater
!    -----------------------------------------------------
!
! Rainfall flowing into uppermost snow layer:
! If rainfall is excessive enough (or layers thin enough)
! it is simply routed directly to runoff: First calculate
! the total snow pack available liquid water holding capacity:
!
ZTOTWCAP(:)   = 0.
DO JJ=1,INLVLS
   DO JI=1,INI
      ZTOTWCAP(JI) = ZTOTWCAP(JI) + ZWHOLDMAX(JI,JJ)
   ENDDO
ENDDO
!
! Rain entering snow (m):
!
ZRAINFALL(:)  = PRR(:)*PTSTEP/XRHOLW                ! rainfall (m)
!
ZFLOWLIQT(:,0)= MIN(ZRAINFALL(:),ZTOTWCAP(:))
!
! Rain assumed to directly pass through the pack to runoff (m):
!
ZPCPXS(:)     = ZRAINFALL(:) - ZFLOWLIQT(:,0)
!
DO JJ=1,INLVLS
   DO JI=1,INI
      ZFLOWLIQT(JI,JJ) = ZFLOWLIQ(JI,JJ)
   ENDDO
ENDDO
!
!
! Thickness is maintained during water through-flow,
! so that mass transfer is represented by
! density changes: NOTE a maximum density
! is assumed (XRHOSMAX_ES) so that all flow
! which would result in densities exceeding
! this limit are routed to next layer down.
! First test for saturation, then
! rout excess water down to next layer down
! and repeat calculation. Net gain in liquid (mass) is
! translated into a density increase:
!
ZFLOWLIQ(:,:)  = 0.0                ! clear this array for work
PSNOWLIQ(:,:)  = ZSNOWLIQ(:,:)      ! reset liquid water content
!
DO JJ=1,INLVLS
   DO JI=1,INI
      ZSNOWLIQ(JI,JJ)  = ZSNOWLIQ(JI,JJ) + ZFLOWLIQT(JI,JJ-1)
      ZFLOWLIQ(JI,JJ)  = MAX(0.0, ZSNOWLIQ(JI,JJ)-ZWHOLDMAX(JI,JJ))
      ZSNOWLIQ(JI,JJ)  = ZSNOWLIQ(JI,JJ) - ZFLOWLIQ(JI,JJ)
      ZSNOWRHO(JI,JJ)  = ZSNOWRHO(JI,JJ)  + (ZSNOWLIQ(JI,JJ) - PSNOWLIQ(JI,JJ))*       &
                          XRHOLW/MAX(XSNOWDMIN/INLVLS,ZSNOWDZ(JI,JJ))  
      ZFLOWLIQT(JI,JJ) = ZFLOWLIQT(JI,JJ) + ZFLOWLIQ(JI,JJ)
   ENDDO
ENDDO
!
!
! Any remaining throughflow after freezing is available to
! the soil for infiltration or surface runoff (m).
! I.E. This is the amount of water leaving the snowpack:
! Rate water leaves the snowpack [kg/(m2 s)]:
!
PTHRUFAL(:)  = PTHRUFAL(:) + ZFLOWLIQT(:,INLVLS)
!
! Add excess rain (rain which flowed directly through the snow
! due to saturation):
!
PTHRUFAL(:)  = (PTHRUFAL(:) + ZPCPXS(:))*XRHOLW/PTSTEP
!
!
! 4. Update thickness and density and any freezing:
!    ----------------------------------------------
!
PSNOWDZ(:,:)  = ZSNOWDZ(:,:)
PSNOWRHO(:,:) = ZSNOWRHO(:,:)
PSNOWLIQ(:,:) = ZSNOWLIQ(:,:)
IF (LHOOK) CALL DR_HOOK('SNOW3LREFRZ',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOW3LREFRZ
!####################################################################
!####################################################################
!####################################################################
      SUBROUTINE SNOW3LFLUX(PSNOWTEMP,PSNOWDZ,PEXNS,PEXNA,            &
                              PUSTAR2_IC,                               &
                              PTSTEP,PALBT,PSW_RAD,PEMIST,PLWUPSNOW,  &
                              PLW_RAD,PTA,PSFCFRZ,PQA,PHPSNOW,        &
                              PSNOWTEMPO1,PSNOWFLUX,PCT,PRADSINK,     &
                              PQSAT,PDQSAT,PRSRA,                     &
                              PRN,PH,PGFLUX,PLES3L,PLEL3L,PEVAP,      &
                              PUSTAR,OSFCMELT                         )  
!
!
!!    PURPOSE
!!    -------
!     Calculate the surface fluxes (atmospheric/surface).
!     (Noilhan and Planton 1989; Noilhan and Mahfouf 1996)
!
!
USE MODD_CSTS,ONLY : XSTEFAN, XCPD, XLSTT, XLVTT, XTT
!
USE MODE_THERMOS
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP
!
REAL, DIMENSION(:), INTENT(IN)      :: PSNOWDZ, PSNOWTEMPO1, PSNOWFLUX, PCT, &
                                         PRADSINK, PEXNS, PEXNA  
!
REAL, DIMENSION(:), INTENT(IN)      :: PALBT, PSW_RAD, PEMIST, PLW_RAD,        &
                                         PTA, PSFCFRZ, PQA,                    &
                                         PHPSNOW, PQSAT, PDQSAT, PRSRA,        &
                                         PUSTAR2_IC  
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PSNOWTEMP
!
REAL, DIMENSION(:), INTENT(OUT)     :: PRN, PH, PGFLUX, PLES3L, PLEL3L,      &
                                         PEVAP, PLWUPSNOW, PUSTAR  
!
LOGICAL, DIMENSION(:), INTENT(OUT)  :: OSFCMELT
!
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSNOWDZ))      :: ZEVAPC, ZLE, ZSNOWTEMP, ZSMSNOW, ZGFLUX,  &
                                         ZDELTAT, ZSNOWTO3  
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! 0. Initialize:
! --------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3LFLUX',0,ZHOOK_HANDLE)
ZSNOWTEMP(:)  = PSNOWTEMP(:)
ZLE(:)        = 0.0
ZSMSNOW(:)    = 0.0
ZGFLUX(:)     = 0.0
!
OSFCMELT(:)   = .FALSE.
!
ZSNOWTO3(:)   = PSNOWTEMPO1(:) ** 3  ! to save some CPU time, store this
!
! 1. Flux calculations when melt not occuring at surface (W/m2):
! --------------------------------------------------------------
!
!
ZDELTAT(:)   = PSNOWTEMP(:) - PSNOWTEMPO1(:)   ! surface T time change:
!
PLWUPSNOW(:) = PEMIST(:) * XSTEFAN * ZSNOWTO3(:)*( PSNOWTEMPO1(:) + 4.* ZDELTAT(:) )
!
PRN(:)       = (1. - PALBT(:)) * PSW_RAD(:) + PEMIST(:) * PLW_RAD(:) -  PLWUPSNOW(:)
!
PH(:)        = PRSRA(:) * XCPD * (PSNOWTEMP(:)/PEXNS(:) - PTA(:)/PEXNA(:))
!
ZEVAPC(:)    = PRSRA(:) * ( (PQSAT(:) - PQA(:)) + PDQSAT(:)*ZDELTAT(:) )
!
PLES3L(:)    = PSFCFRZ(:)     * XLSTT * ZEVAPC(:)
!
PLEL3L(:)    = (1.-PSFCFRZ(:))* XLVTT * ZEVAPC(:)
!
ZLE(:)       = PLES3L(:) + PLEL3L(:)
!
PGFLUX(:)    = PRN(:) - PH(:) - ZLE(:) + PHPSNOW(:)
!
!
! 2. Initial melt adjustment
! --------------------------
! If energy avalabile to melt snow, then recalculate fluxes
! at the freezing point and add residual heat to layer
! average heat.
!
! A) If temperature change is > 0 and passes freezing point this timestep,
!    then recalculate fluxes at freezing point and excess energy
!    will be used outside of this routine to change snow heat content:
!
WHERE (PSNOWTEMP > XTT .AND. PSNOWTEMPO1 < XTT)
!
   OSFCMELT(:)= .TRUE.
!
   ZDELTAT(:)    = XTT - PSNOWTEMPO1(:)
!
   PLWUPSNOW(:) = PEMIST(:) * XSTEFAN * ZSNOWTO3(:)*( PSNOWTEMPO1(:) + 4.* ZDELTAT(:) ) 
!
   PRN(:)       = (1. - PALBT(:)) * PSW_RAD(:) + PEMIST(:) * PLW_RAD(:) - PLWUPSNOW(:)
!
   PH(:)        = PRSRA(:) * XCPD * (XTT/PEXNS(:) - PTA(:)/PEXNA(:))   
!
   ZEVAPC(:)    = PRSRA(:) * ( (PQSAT(:) - PQA(:)) + PDQSAT(:)*ZDELTAT(:) )
!
   PLES3L(:)    = PSFCFRZ(:)     * XLSTT * ZEVAPC(:)
!
   PLEL3L(:)    = (1.-PSFCFRZ(:))* XLVTT * ZEVAPC(:)
!
   ZLE(:)       = PLES3L(:) + PLEL3L(:)
!
   ZGFLUX(:)    = PRN(:) - PH(:) - ZLE(:) + PHPSNOW(:)
!
   ZSMSNOW(:)   = PGFLUX(:) - ZGFLUX(:)
!
   PGFLUX(:)  = ZGFLUX(:)
!
! This will be used to change heat content of snow:
!
   ZSNOWTEMP(:) = PSNOWTEMP(:) - ZSMSNOW(:)*PTSTEP*PCT(:)
!
END WHERE
!
! 3. Ongoing melt adjustment: explicit solution
! ---------------------------------------------
!    If temperature change is 0 and at freezing point this timestep,
!    then recalculate fluxes and surface temperature *explicitly*
!    as this is *exact* for snow at freezing point (Brun, Martin)
!
WHERE(PSNOWTEMP(:) > XTT .AND. PSNOWTEMPO1(:) >= XTT)
!
   OSFCMELT(:)  = .TRUE.   
!
   PLWUPSNOW(:) = PEMIST(:) * XSTEFAN * (XTT ** 4) 
!
   PRN(:)       = (1. - PALBT(:)) * PSW_RAD(:) + PEMIST(:) * PLW_RAD(:) - PLWUPSNOW(:)
!
   PH(:)        = PRSRA(:) * XCPD * (XTT/PEXNS(:) - PTA(:)/PEXNA(:))
!
   ZEVAPC(:)    = PRSRA(:) * (PQSAT(:) - PQA(:))
!
   PLES3L(:)    = PSFCFRZ(:)     * XLSTT * ZEVAPC(:)
!
   PLEL3L(:)    = (1.-PSFCFRZ(:))* XLVTT * ZEVAPC(:)
!
   ZLE(:)       = PLES3L(:) + PLEL3L(:)
!
   PGFLUX(:)    = PRN(:) - PH(:) - ZLE(:) + PHPSNOW(:)
!
   ZSNOWTEMP(:) = XTT + PTSTEP*PCT(:)*(PGFLUX(:) + PRADSINK(:) - PSNOWFLUX(:))
!
END WHERE
!
! 4. Update surface temperature:
! ------------------------------
!
PSNOWTEMP(:) = ZSNOWTEMP(:)
!
! 5. Final evaporative flux (kg/m2/s)
!
PEVAP(:)     = ZEVAPC(:)
!
! 6. Friction velocity
! --------------------
!
PUSTAR(:) = SQRT(PUSTAR2_IC(:)) 
!
IF (LHOOK) CALL DR_HOOK('SNOW3LFLUX',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOW3LFLUX
!####################################################################
!####################################################################
!####################################################################
      SUBROUTINE SNOW3LEVAPN(PPSN3L,PLES3L,PLEL3L,PTSTEP,PSNOWTEMP, &
                               PSNOWRHO,PSNOWDZ,PSNOWLIQ,PTA,       &
                               PSNOWHEAT,PEVAPCOR,PGFLXCOR)
!
!
!!    PURPOSE
!!    -------
!     Remove mass from uppermost snow layer in response to
!     evaporation (liquid) and sublimation.
!
!
USE MODD_CSTS,     ONLY : XLVTT, XRHOLW, XLSTT, XLMTT, XCI, XTT
USE MODD_SNOW_PAR, ONLY : XRHOSMIN_ES, XSNOWDMIN
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP
!
REAL, DIMENSION(:), INTENT(IN)      :: PPSN3L 
!
REAL, DIMENSION(:), INTENT(IN)      :: PLES3L, PLEL3L   ! (W/m2)
!
REAL, DIMENSION(:), INTENT(IN)      :: PTA
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWHEAT, PSNOWDZ
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PSNOWRHO, PSNOWLIQ, &
                                       PEVAPCOR, PSNOWTEMP
REAL, DIMENSION(:), INTENT(OUT)     :: PGFLXCOR 
!
!*      0.2    declarations of local variables
!
INTEGER                             :: INI, INLVLS, JJ, JI
!
REAL, DIMENSION(SIZE(PLES3L))       :: ZSNOWEVAPS, ZSNOWEVAP, ZSNOWEVAPX,          &
                                         ZSNOWDZ, ZEVAPCOR, ZSNOWHEAT, ZSCAP ,     &
                                         ZSNOWTEMP
!                                      ZEVAPCOR = for vanishingy thin snow cover,
!                                                 allow any excess evaporation
!                                                 to be extracted from the soil
!                                                 to maintain an accurate water
!                                                 balance [kg/(m2 s)]


REAL, DIMENSION(SIZE(PLES3L))       :: ZXSE, ZISNOWD

!*      0.3    declarations of local parameters

REAL, PARAMETER                     :: ZSNOWDEMIN = 1.E-4 ! m
REAL, PARAMETER                     :: ZTDIF      = 15.   ! K To prevent a possible
                                                          !   decoupling of sfc snow T
                                                          !   when vanishingly thin, impose
                                                          !   this max T-diff based on obs...
                                                          !   between Ta-Ts


REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!-------------------------------------------------------------------------------
!
! 0. Initialize:
! --------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3LEVAPN',0,ZHOOK_HANDLE)
ZEVAPCOR(:)    = 0.0
ZSNOWEVAPS(:)  = 0.0
ZSNOWEVAP(:)   = 0.0
ZSNOWEVAPX(:)  = 0.0
ZSNOWDZ(:)     = 0.0
ZSCAP(:)       = 0.0
ZSNOWHEAT(:)   = 0.0
ZSNOWTEMP(:)   = 0.0
ZXSE(:)        = 0.0
PGFLXCOR(:)    = 0.0
!
INI            = SIZE(PSNOWDZ(:,:),1)
INLVLS         = SIZE(PSNOWDZ(:,:),2)
!
WHERE(PSNOWDZ(:,1) > 0.0)
!
! 1. Evaporation of snow liquid water
! -----------------------------------
!     Evaporation reduces liq water equivalent content
!     of snow pack either by reducing the snow density
!     (evaporation) or the layer thickness (sublimation).
!     Condensation does the reverse.
!
!     Multiply evaporation components by snow fraction
!     to be consistent with fluxes from the snow covered
!     fraction of grid box
!
! Evaporation:
! Reduce density and liquid water content:
!
   ZSNOWEVAP(:)   = PPSN3L(:)*PLEL3L(:)*PTSTEP/(XLVTT*XRHOLW)
   ZSNOWEVAPX(:)  = MIN(ZSNOWEVAP(:),PSNOWLIQ(:))
!
!  This should not change enthalpy (since already accounted for 
!  in sfc e budget), so make density change insuring constant enthalpy:
!
   PSNOWLIQ(:)    = PSNOWLIQ(:) - ZSNOWEVAPX(:)
   PSNOWRHO(:) = (PSNOWHEAT(:,1)-XLMTT*XRHOLW*PSNOWLIQ(:))/ &
                 (PSNOWDZ(:,1)*(XCI*(PSNOWTEMP(:)-XTT)-XLMTT))
!
! Budget check: If density drops below minimum, then extract the
! corresponding water mass from soil (for vanishingly thin snow covers):
!
   ZEVAPCOR(:)    = MAX(0.0,XRHOSMIN_ES-PSNOWRHO(:))*PSNOWDZ(:,1)/PTSTEP
   PSNOWRHO(:)    = MAX(XRHOSMIN_ES,PSNOWRHO(:))
!
!
!! 2. Sublimation of snow ice
! ---------------------------

! Update heat capacity:

   ZSCAP(:)       = SNOW3LSCAP(PSNOWRHO)
!
! Budget check: as last traces of liquid in snow evaporates, it is possible
! evaporation could exceed liquid present in snow. If this is the case,
! remove residual mass from solid snow in order to maintain high
! accuracy water balance:

   ZSNOWEVAPX(:)  = MAX(0.0, ZSNOWEVAP(:) - ZSNOWEVAPX(:))
   ZSNOWDZ(:)     = PSNOWDZ(:,1) - ZSNOWEVAPX(:)*XRHOLW/PSNOWRHO(:)
   PSNOWDZ(:,1)   = MAX(0.0, ZSNOWDZ(:))
   ZEVAPCOR(:)    = ZEVAPCOR(:) + MAX(0.0,-ZSNOWDZ(:))*PSNOWRHO(:)/PTSTEP
   
! Sublimation: Reduce layer thickness and total snow depth
! if sublimation: add to correction term if potential
! sublimation exceeds available snow cover.
!
   ZSNOWEVAPS(:)  = PPSN3L(:)*PLES3L(:)*PTSTEP/(XLSTT*PSNOWRHO(:))
   ZSNOWDZ(:)     = PSNOWDZ(:,1) - ZSNOWEVAPS(:)
   PSNOWDZ(:,1)   = MAX(0.0, ZSNOWDZ(:))
   ZEVAPCOR(:)    = ZEVAPCOR(:) + MAX(0.0,-ZSNOWDZ(:))*PSNOWRHO(:)/PTSTEP

! Note that the effect on enthalpy is already accounted for via the surface energy
! budget, so the change in enthalpy
! here owing to a mass loss must be offset by possible adjustments to T and Wl
! to conserve Enthalpy.

! For vanishingly thin snow, a decoupling can occur between forcing level T
! and snow sfc T (as snowpack vanishes owing to sublimation), so a simple limit
! imposed on this gradient:

   PSNOWTEMP(:)   = XTT + ( ((PSNOWHEAT(:,1)/MAX(ZSNOWDEMIN,PSNOWDZ(:,1)))   &
                    + XLMTT*PSNOWRHO(:))/ZSCAP(:) )

   PSNOWTEMP(:)   = MAX(MIN(XTT,PTA(:)-ZTDIF), PSNOWTEMP(:))

! Update surface enthalpy and re-diagnose T and liq:

   PSNOWLIQ(:)    = MAX(0.0,PSNOWTEMP(:)-XTT)*ZSCAP(:)*                  &
                    PSNOWDZ(:,1)/(XLMTT*XRHOLW)  

   PSNOWTEMP(:)   = MIN(XTT,PSNOWTEMP(:)) 

   ZSNOWHEAT(:)   = PSNOWHEAT(:,1)
   PSNOWHEAT(:,1) = PSNOWDZ(:,1)*( ZSCAP(:)*(PSNOWTEMP(:)-XTT)             &
                    - XLMTT*PSNOWRHO(:) ) + XLMTT*XRHOLW*PSNOWLIQ(:) 

   ZXSE(:)        = PSNOWHEAT(:,1) - ZSNOWHEAT(:) ! excess cooling (removed from sfc)
 
END WHERE
!
! If the sfc T-Gradient was limited (and thus the enthalpy in the uppermost layer changed),
! distribute this energy correction (cooling) in the layers below (thickness-based weighting)
! to conserve total snowpack enthalpy. Normally this only occurs for this snowpacks,
! but it can rarely occur otherwise also.
!
ZISNOWD(:) = 0.
DO JJ=2,INLVLS
   DO JI=1,INI
      ZISNOWD(JI) = ZISNOWD(JI) + PSNOWDZ(JI,JJ)        ! m
   ENDDO
END DO
ZISNOWD(:)        = ZXSE(:)/MAX(ZISNOWD(:),ZSNOWDEMIN)  ! J kg-1 m-1 
!
DO JJ=2,INLVLS
   DO JI=1,INI
      PSNOWHEAT(JI,JJ) = PSNOWHEAT(JI,JJ) - PSNOWDZ(JI,JJ)*ZISNOWD(JI)
   ENDDO
ENDDO
!
!
! 3. Update evaporation correction term:
! --------------------------------------
!
PEVAPCOR(:) = PEVAPCOR(:) + ZEVAPCOR(:)
IF (LHOOK) CALL DR_HOOK('SNOW3LEVAPN',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOW3LEVAPN
!####################################################################
!####################################################################
!####################################################################
SUBROUTINE SNOW3LGONE(PTSTEP,PPSN3L,PLEL3L,PLES3L,PSNOWRHO,         &
                   PSNOWHEAT,PRADSINK,PEVAPCOR,PTHRUFAL,PGRNDFLUX,    &
                   PGFLUXSNOW,PGRNDFLUXO,PSNOWDZ,PSNOWLIQ,PSNOWTEMP,  &
                   PRADXS)  
!
!
!
!!    PURPOSE
!!    -------
!     Account for the case when the last trace of snow melts
!     during a time step: ensure mass and heat balance of
!     snow AND underlying surface.
!
!
USE MODD_CSTS,ONLY : XTT, XLSTT, XLVTT
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP
!
REAL, DIMENSION(:), INTENT(IN)      :: PPSN3L, PLEL3L, PLES3L, PGFLUXSNOW, &
                                         PRADSINK, PGRNDFLUXO  
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWRHO, PSNOWHEAT
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PGRNDFLUX, PRADXS
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWDZ, PSNOWLIQ, PSNOWTEMP
!
REAL, DIMENSION(:), INTENT(OUT)     :: PTHRUFAL   ! melt water [kg/(m2 s)]
!
REAL, DIMENSION(:), INTENT(OUT)     :: PEVAPCOR   ! [kg/(m2 s)]
!                                      PEVAPCOR = for vanishingy thin snow cover,
!                                                 allow any excess evaporation
!                                                 to be extracted from the soil
!                                                 to maintain an accurate water
!                                                 balance.
!
!*      0.2    declarations of local variables
!
INTEGER                             :: JJ, JI
!
INTEGER                             :: INI
INTEGER                             :: INLVLS
!
REAL, DIMENSION(SIZE(PLES3L))       :: ZSNOWHEATC, ZSNOWGONE_DELTA
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! 0. Initialize:
! --------------
!
IF (LHOOK) CALL DR_HOOK('SNOW3LGONE',0,ZHOOK_HANDLE)
PEVAPCOR(:)           = 0.0
PTHRUFAL(:)           = 0.0
INI                   = SIZE(PSNOWDZ(:,:),1)
INLVLS                = SIZE(PSNOWDZ(:,:),2)
ZSNOWHEATC(:)         = 0.
DO JJ=1,INLVLS
   DO JI=1,INI
      ZSNOWHEATC(JI) = ZSNOWHEATC(JI) + PSNOWHEAT(JI,JJ) ! total heat content (J m-2)
   ENDDO
ENDDO
ZSNOWGONE_DELTA(:)    = 1.0
!
! 1. Simple test to see if snow vanishes:
! ---------------------------------------
! If so, set thicknesses (and therefore mass and heat) and liquid content
! to zero, and adjust fluxes of water, evaporation and heat into underlying
! surface. Note, test with flux computed *before* correction since this represents
! actual inflow of heat from below (as heat content correction owing to a corrected
! flux has not yet been done: here we compare to pre-corrected heat content).
!
WHERE(PGFLUXSNOW(:) - PGRNDFLUXO(:) + PRADSINK(:) >= (-ZSNOWHEATC(:)/PTSTEP) )
   PGRNDFLUX(:)       = PGFLUXSNOW(:) + (ZSNOWHEATC(:)/PTSTEP)
   PEVAPCOR(:)        = PPSN3L(:)*((PLEL3L(:)/XLVTT) + (PLES3L(:)/(XLSTT*PSNOWRHO(:,1))))
   PRADXS(:)          = 0.0
   ZSNOWGONE_DELTA(:) = 0.0          ! FLAG...if=0 then snow vanishes, else=1
END WHERE
!
PTHRUFAL(:) = 0.
DO JJ=1,INLVLS
   DO JI=1,INI
      PTHRUFAL(JI) = PTHRUFAL(JI) + (1.0-ZSNOWGONE_DELTA(JI))*PSNOWRHO(JI,JJ)*PSNOWDZ(JI,JJ)/PTSTEP
   ENDDO
END DO
!
! 2. Final update of snow state
! -----------------------------
! (either still present or not):
!
DO JJ=1,INLVLS
   DO JI=1,INI
      PSNOWDZ  (JI,JJ) =                                 PSNOWDZ  (JI,JJ)*ZSNOWGONE_DELTA(JI)
      PSNOWLIQ (JI,JJ) =                                 PSNOWLIQ (JI,JJ)*ZSNOWGONE_DELTA(JI)
      PSNOWTEMP(JI,JJ) = (1.0-ZSNOWGONE_DELTA(JI))*XTT + PSNOWTEMP(JI,JJ)*ZSNOWGONE_DELTA(JI)
   ENDDO
ENDDO
IF (LHOOK) CALL DR_HOOK('SNOW3LGONE',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE SNOW3LGONE
!####################################################################
!####################################################################
!####################################################################
SUBROUTINE SNOW3LEVAPGONE(PSNOWHEAT,PSNOWDZ,PSNOWRHO,PSNOWTEMP,PSNOWLIQ)
!
!!    PURPOSE
!!    -------
!
!     If all snow in uppermost layer evaporates/sublimates, re-distribute
!     grid (below assumes very thin snowpacks so layer-thicknesses are
!     constant).
!
!
USE MODD_CSTS,     ONLY : XTT, XRHOLW, XLMTT
USE MODD_SNOW_PAR, ONLY : XRHOSMIN_ES, XSNOWDMIN, XRHOSMAX_ES
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWRHO   ! snow density profile                (kg/m3)
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWDZ    ! snow layer thickness profile        (m)
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWHEAT  ! snow heat content/enthalpy          (J/m2)
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWTEMP  ! snow temperature profile            (K)
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWLIQ   ! snow liquid water profile           (m)
!
!*      0.2    declarations of local variables
!
INTEGER                               :: JJ, JI
!
INTEGER                               :: INI
INTEGER                               :: INLVLS
!
REAL, DIMENSION(SIZE(PSNOWDZ,1))      :: ZSNOWHEAT_1D ! total heat content                (J/m2)
REAL, DIMENSION(SIZE(PSNOWDZ,1))      :: ZSNOW        ! total snow depth                  (m)
REAL, DIMENSION(SIZE(PSNOWDZ,1))      :: ZMASS        ! total mass
!
REAL, DIMENSION(SIZE(PSNOWDZ,1),SIZE(PSNOWDZ,2)) :: ZSCAP  ! Snow layer heat capacity          (J/K/m3)
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! Initialize:
!
IF (LHOOK) CALL DR_HOOK('SNOW3LEVAPGONE',0,ZHOOK_HANDLE)
INI             = SIZE(PSNOWDZ,1)
INLVLS          = SIZE(PSNOWDZ,2)
!
! First, determine where uppermost snow layer has completely
! evaporated/sublimated (as it becomes thin):
!
ZSNOWHEAT_1D(:) = 0.
ZSNOW(:)        = 0.
ZMASS(:)        = 0.
!
ZSCAP(:,:) = SNOW3LSCAP(PSNOWRHO(:,:))
!
DO JJ=2,INLVLS
   DO JI=1,INI
      IF(PSNOWDZ(JI,1) == 0.0)THEN
         ZSNOWHEAT_1D(JI) = ZSNOWHEAT_1D(JI) + XLMTT*XRHOLW*PSNOWLIQ(JI,JJ)  &
                          + PSNOWDZ(JI,JJ)*(ZSCAP(JI,JJ)*(PSNOWTEMP(JI,JJ)-XTT) &
                          - XLMTT*PSNOWRHO(JI,JJ) )           
         ZSNOW       (JI) = ZSNOW(JI) + PSNOWDZ(JI,JJ)
         ZMASS       (JI) = ZMASS(JI) + PSNOWDZ(JI,JJ)*PSNOWRHO(JI,JJ)    
       ENDIF
    ENDDO
ENDDO
!
! Where uppermost snow layer has vanished, redistribute vertical
! snow mass and heat profiles (and associated quantities):
!
DO JJ=1,INLVLS
   DO JI=1,INI
      IF(ZSNOW(JI)/= 0.0)THEN
        ZSNOW    (JI)    = MAX(0.5*XSNOWDMIN,ZSNOW(JI))
        PSNOWDZ  (JI,JJ) = ZSNOW(JI)/REAL(INLVLS)
        PSNOWHEAT(JI,JJ) = ZSNOWHEAT_1D(JI)/REAL(INLVLS)
        PSNOWRHO (JI,JJ) = ZMASS (JI)/ZSNOW(JI)
      ENDIF
    ENDDO
ENDDO        
!
ZSCAP(:,:) = SNOW3LSCAP(PSNOWRHO(:,:))
!
DO JJ=1,INLVLS
   DO JI=1,INI
      IF(ZSNOW(JI)/= 0.0)THEN
        PSNOWTEMP(JI,JJ) = XTT + ( ((PSNOWHEAT(JI,JJ)/PSNOWDZ(JI,JJ))   &
                               + XLMTT*PSNOWRHO(JI,JJ))/ZSCAP(JI,JJ) )  
        PSNOWLIQ (JI,JJ) = MAX(0.0,PSNOWTEMP(JI,JJ)-XTT)*ZSCAP(JI,JJ)*  &
                                   PSNOWDZ(JI,JJ)/(XLMTT*XRHOLW)  
        PSNOWTEMP(JI,JJ) = MIN(XTT,PSNOWTEMP(JI,JJ))
      ENDIF
    ENDDO
ENDDO
IF (LHOOK) CALL DR_HOOK('SNOW3LEVAPGONE',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOW3LEVAPGONE
!
!
!
END SUBROUTINE SNOW3L
