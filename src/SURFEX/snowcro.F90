!     ##########################################################################
      SUBROUTINE SNOWCRO(HSNOWRES, TPTIME, OGLACIER,                    &
               PPEW_A_COEF, PPEW_B_COEF,                                 &
               PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,       &
               PSNOWSWE,PSNOWRHO,PSNOWHEAT,PSNOWALB,                     &
               PSNOWGRAN1,PSNOWGRAN2,PSNOWHIST,PSNOWAGE,                 &
               PTSTEP,PPS,PSR,PRR,PPSN3L,                                &
               PTA,PTG,PSW_RAD,PQA,PVMOD,PLW_RAD, PRHOA,                 &
               PUREF,PEXNS,PEXNA,PDIRCOSZW,                              &
               PZREF,PZ0,PZ0EFF,PZ0H,PALB,                               &
               PSOILCOND,PD_G,                                           &
               PSNOWLIQ,PSNOWTEMP,PSNOWDZ,                               &
               PTHRUFAL,PGRNDFLUX,PEVAPCOR,PRNSNOW,PHSNOW,PGFLUXSNOW,    &
               PHPSNOW,PLES3L,PLEL3L,PEVAP,                              &
               PEMISNOW,PCDSNOW,PUSTAR,PCHSNOW,PSNOWHMASS,               &
               PPERMSNOWFRAC,PZENITH,  PXLAT, PXLON                      ) 
!     ##########################################################################
!
!!****  *SNOWCRO*
!!
!!    PURPOSE
!!    -------
!
!     Detailed snowpack scheme Crocus, computationnally based on the
!     3-Layer snow scheme option (Boone and Etchevers 1999)
!     For shallow snow cover, Default method of Douville et al. (1995)
!     used with this option: Model "turns on" when snow sufficiently
!     deep/above some preset critical snow depth.
!
!
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
!!    ISBA: Belair (1995)
!!    ISBA: Noilhan and Planton (1989)
!!    ISBA: Noilhan and Mahfouf (1996)
!!    ISBA-ES: Boone and Etchevers (2001)
!!    Crocus : Brun et al., 1989 (J. Glaciol.)
!!    Crocus : Brun et al., 1992 (J. Glaciol.)
!!    Crocus : Vionnet et al., in prep (Geosci. Mod. Devel. Discuss.)
!!
!!
!!    AUTHOR
!!    ------
!!      A. Boone           * Meteo-France *
!!      V. Vionnet         * Meteo-France *
!!      E. Brun            * Meteo-France *
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
!!      Modified by E. Brun, V. Vionnet, S. Morin (05/2011): 
!!                                          Addition of Crocus processes and 
!!                                          parametrizations to 
!!                                          the SNOW-3L code. This includes the dynamic handling
!!                                          of snow layers and the inclusion of snow metamorphism
!!                                          rules similar to the original Crocus implementation.
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
USE MODD_CSTS,     ONLY : XTT, XRHOLW, XLMTT,XLSTT,XLVTT, XCL, XCI
!
USE MODE_SNOW3L
!
USE MODD_SNOW_METAMO
!
USE MODD_TYPE_DATE_SURF, ONLY: DATE_TIME
!
USE MODI_ABOR1_SFX
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_GREGODSTRATI
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
LOGICAL, INTENT(IN)                 :: OGLACIER   ! True = Over permanent snow and ice, 
!                                                     initialise WGI=WSAT,
!                                                     Hsnow>=10m and allow 0.8<SNOALB<0.85
                                                  ! False = No specific treatment
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
REAL, DIMENSION(:), INTENT(IN)    :: PTG, PSOILCOND, PD_G, PPSN3L
!                                      PTG       = Surface soil temperature (effective
!                                                  temperature the of layer lying below snow)
!                                      PSOILCOND = soil thermal conductivity [W/(m K)]
!                                      PD_G      = Assumed first soil layer thickness (m)
!                                                  Used to calculate ground/snow heat flux
!                                      PPSN3L    = snow fraction
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
!                                      PPEW_A_COEF = wind coefficient
!                                      PPEW_B_COEF = wind coefficient
!                                      PPET_A_COEF = A-air temperature coefficient
!                                      PPET_B_COEF = B-air temperature coefficient
!                                      PPEQ_A_COEF = A-air specific humidity coefficient
!                                      PPEQ_B_COEF = B-air specific humidity coefficient
!
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
REAL, DIMENSION(:,:), INTENT(OUT) :: PSNOWLIQ, PSNOWTEMP, PSNOWDZ
!                                      PSNOWLIQ  = Snow layer(s) liquid water content (m)
!                                      PSNOWTEMP = Snow layer(s) temperature (m)
!                                      PSNOWDZ   = Snow layer(s) thickness (m)
!
REAL, DIMENSION(:), INTENT(OUT)   :: PTHRUFAL, PGRNDFLUX, PEVAPCOR
!                                      PTHRUFAL  = rate that liquid water leaves snow pack:
!                                                  paritioned into soil infiltration/runoff
!                                                  by ISBA [kg/(m2 s)]
!                                      PGRNDFLUX = soil/snow interface heat flux (W/m2)
!                                      PEVAPCOR  = evaporation/sublimation correction term:
!                                                  extract any evaporation exceeding the
!                                                  actual snow cover (as snow vanishes)
!                                                  and apply it as a surface soil water
!                                                  sink. [kg/(m2 s)]
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRNSNOW, PHSNOW, PGFLUXSNOW, PLES3L, PLEL3L, &
                                        PHPSNOW, PCDSNOW, PUSTAR, PEVAP 
!                                      PLES3L      = evaporation heat flux from snow (W/m2)
!                                      PLEL3L      = sublimation (W/m2)
!                                      PHPSNOW     = heat release from rainfall (W/m2)
!                                      PRNSNOW     = net radiative flux from snow (W/m2)
!                                      PHSNOW      = sensible heat flux from snow (W/m2)
!                                      PGFLUXSNOW  = net heat flux from snow (W/m2)
!                                      PCDSNOW     = drag coefficient for momentum over snow
!                                      PUSTAR      = friction velocity over snow (m/s)
!                                      PEVAP       = total evaporative flux (kg/m2/s)
!
REAL, DIMENSION(:), INTENT(OUT)   :: PCHSNOW, PEMISNOW, PSNOWHMASS
!                                      PEMISNOW    = snow surface emissivity
!                                      PCHSNOW     = drag coefficient for heat over snow
!                                      PSNOWHMASS  = heat content change due to mass
!                                                    changes in snowpack (J/m2): for budget
!                                                    calculations only.
REAL, DIMENSION(:), INTENT(IN)    :: PZENITH ! solar zenith angle
REAL, DIMENSION(:), INTENT(IN)    :: PXLAT,PXLON ! LAT/LON after packing
!
!
!*      0.2    declarations of local variables
!
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWTEMP, ZSCAP, ZSNOWDZN, ZSCOND,    &
                                                       ZRADSINK 
!                                      ZSNOWTEMP  = Snow layer(s) averaged temperature (K)
!                                      ZSCAP      = Snow layer(s) heat capacity [J/(K m3)]
!                                      ZSNOWDZN   = Updated snow layer thicknesses (m)
!                                      ZSCOND     = Snow layer(s) thermal conducivity [W/(m K)]
!                                      ZRADSINK   = Snow solar Radiation source terms (W/m2)
!
!  
REAL, DIMENSION(SIZE(PTA))          :: ZSNOWBIS
!                                     ZSNOWBIS      = Total snow depth after snowfall  
!
REAL, DIMENSION(SIZE(PTA))          :: ZSNOW, ZSFCFRZ, ZTSTERM1, ZTSTERM2,                   &
                                        ZCT, ZRA, ZSNOWFLUX, ZSNOWTEMPO1 
!                                      ZSNOW      = Total snow depth (m)
!                                      ZCT        = inverse of the product of snow heat capacity
!                                                   and layer thickness [(m2 K)/J]
!                                      ZRA        = Surface aerodynamic resistance
!                                      ZTSTERM1,ZTSTERM2 = Surface energy budget coefficients
!                                      ZSNOWFLUX  = heat flux between 1st and 2nd snow layers:
!                                                   used during surface melting (W/m2)
!                                      ZSNOWTEMPO1= value of uppermost snow temperature
!                                                   before time integration (K)
!
LOGICAL, DIMENSION(SIZE(PTA))       :: GSNOWFALL,GMODIF_MAILLAGE
!                                      GSNOWFALL  = FLAG if snowfall exceed PSNOW/10, used for 
!                                                   grid updating. 
!
REAL, DIMENSION(SIZE(PTA))          :: ZRSRA, ZDQSAT, ZQSAT, ZRADXS, ZLIQHEATXS, &
                                        ZLWUPSNOW 
!                                      ZRSRA    = air density over aerodynamic resistance
!                                      ZDQSAT   = derrivative of saturation specific humidity
!                                      ZQSAT    = saturation specific humidity
!                                      ZRADXS   = shortwave radiation absorbed by soil surface
!                                                 (for thin snow sover) (W m-2)
!                                      ZLIQHEATXS = excess snowpack heating for vanishingly thin
!                                                 snow cover: add energy to snow/ground heat
!                                                 flux (W m-2)
!                                      ZLWUPSNOW = upwelling longwave raaditive flux (W m-2)
!
REAL, DIMENSION(SIZE(PTA))          :: ZVMOD_IC, ZTA_IC, ZQA_IC,                                   &
                                        ZPET_A_COEF_T, ZPEQ_A_COEF_T, ZPET_B_COEF_T, ZPEQ_B_COEF_T  
!                                      ZVMOD_IC      = implicit lowest atmospheric level wind speed
!                                      ZTA_IC        = implicit lowest atmospheric level air temperature
!                                      ZQA_IC        = implicit lowest atmospheric level specific humidity
!                                      ZPET_A_COEF_T = transformed A-air temperature coefficient
!                                      ZPET_B_COEF_T = transformed B-air temperature coefficient
!                                      ZPEQ_A_COEF_T = transformed A-air specific humidity coefficient
!                                      ZPEQ_B_COEF_T = transformed B-air specific humidity coefficient
!
REAL, PARAMETER                     :: ZSNOWDZMIN = 0.0001 ! (m)
!                                      ZSNOWDZMIN = minimum snow layer thickness for
!                                                   thermal calculations. Used to prevent
!                                                   numerical problems as snow becomes
!                                                   vanishingly thin.
!
INTEGER                             :: JJ,JST   ! looping indexes
!
LOGICAL                             :: OCOND_GRAIN,ODRIFT_GRAIN, OCOND_YEN
REAL, DIMENSION(SIZE(PTA))          :: ZSNOWRHOF, ZSNOWDZF,ZSNOWGRAN1F,&
                                        ZSNOWGRAN2F, ZSNOWHISTF 
REAL, DIMENSION(SIZE(PTA))          :: ZSNOWAGEF
INTEGER, DIMENSION(SIZE(PTA)) :: INLVLS_USE ! varying number of effective layers 
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2))  :: ZWHOLDMAX 
INTEGER                         :: I  ! gridpoint number to be printed
LOGICAL                         ::LP ! flag  to print gridpoints diagnostics
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SNOWCRO',0,ZHOOK_HANDLE)
OCOND_GRAIN=.TRUE.
ODRIFT_GRAIN=.TRUE.
OCOND_YEN=.TRUE.!FALSE. !(if TRUE : use of the Yen (1981) thermal conductivity paramztrization ; 
!                   otherwise, use the default ISBA-ES thermal conductivity parametrization)
!
PGRNDFLUX=0.
PSNOWHMASS=0.
PHSNOW=0.
PRNSNOW=0.
PLES3L=0.
PLEL3L=0.
PHPSNOW=0.
PEVAPCOR=0.
PTHRUFAL=0.
! pour imprimer des diagnostics sur un des points
IF(SIZE(PTA)>= 1) THEN 
        LP=.FALSE.
        I=1
ELSE
        LP=.FALSE.
        I=1
ENDIF        
!
! - - ---------------------------------------------------
!
!       0.     Initialization
!               --------------
! NOTE that snow layer thickness is used throughout this code: SWE
! is only used to diagnose the thickness at the beginning of this routine
! and it is updated at the end of this routine.
!
! Initialization of the actual number of snow layers, total snow depth 
!  and layer thicknesses
!
GSNOWFALL(:)    = .FALSE.
INLVLS_USE(:)=0
DO JJ=1, SIZE(ZSNOW)
  DO JST=1,SIZE(PSNOWSWE(:,:),2) 
     IF (PSNOWSWE(JJ,JST)>0.) THEN 
        PSNOWDZ(JJ,JST) = PSNOWSWE(JJ,JST)/PSNOWRHO(JJ,JST)
        INLVLS_USE(JJ)= JST
     ELSE
        PSNOWDZ(JJ,JST) = 0.
     ENDIF
  ENDDO  !  end loop snow layers
ENDDO    ! end loop grid points
!
!  Print of some simulation characteristics
!
IF(TPTIME%TIME ==0.0) THEN
  WRITE(*,FMT="(I4,2I4,F6.2,A12,I3,A12,I4)")          &
     TPTIME%TDATE%YEAR,TPTIME%TDATE%MONTH,TPTIME%TDATE%DAY,TPTIME%TIME/3600.,&
     'nlayer(I):',INLVLS_USE(I), ' nbpoints:', SIZE(ZSNOW) 
  WRITE(*,*) 'PZ0H: ', PZ0H(I), 'drift ROf 2 x vent '
  write(*,*) 'Fraction neige =',PPSN3L(1)                   
ENDIF
!
!*       1.     Snow total depth
!               ----------------
!
ZSNOW(:) = 0.
DO JJ=1, SIZE(ZSNOW)
   ZSNOW(JJ)=SUM(PSNOWDZ(JJ,1:INLVLS_USE(JJ)))
ENDDO
!
ZSNOWBIS(:)        = ZSNOW(:)
!
!
!*       2.     Snowfall
!               --------
! Calculate uppermost density and thickness changes due to snowfall,
! and add heat content of falling snow
!
!
  CALL SNOWNLFALL_UPGRID(TPTIME, OGLACIER,  &
                  PTSTEP,PSR,PTA,PVMOD,ZSNOWBIS,PSNOWRHO,PSNOWDZ, &
                  PSNOWHEAT,PSNOWHMASS,PSNOWALB,PPERMSNOWFRAC, &
                  PSNOWGRAN1,PSNOWGRAN2,GSNOWFALL, ZSNOWDZN,&
                  ZSNOWRHOF, ZSNOWDZF,ZSNOWGRAN1F, ZSNOWGRAN2F, ZSNOWHISTF,    &
                  ZSNOWAGEF, GMODIF_MAILLAGE,INLVLS_USE) 
!
  ZSNOW(:)      = ZSNOWBIS(:)
!
!*       3.     Update grid/discretization
!               --------------------------
! Reset grid to conform to model specifications:
!
  DO JJ=1, SIZE(ZSNOW)
    IF(GMODIF_MAILLAGE(JJ)) THEN
!        
    CALL  SNOWNLGRIDFRESH_1D(ZSNOW(JJ),PSNOWDZ(JJ,:),ZSNOWDZN(JJ,:),   &
           PSNOWRHO(JJ,:),PSNOWHEAT(JJ,:),PSNOWGRAN1(JJ,:),PSNOWGRAN2(JJ,:), &
           PSNOWHIST(JJ,:),PSNOWAGE(JJ,:),GSNOWFALL(JJ), &
           ZSNOWRHOF(JJ),ZSNOWDZF(JJ),PSNOWHMASS(JJ),ZSNOWGRAN1F(JJ),&
           ZSNOWGRAN2F(JJ), ZSNOWHISTF(JJ),ZSNOWAGEF(JJ),   &
           INLVLS_USE(JJ))                                
    ENDIF
  ENDDO                                       
!  
!
!*       4.     Liquid water content and snow temperature
!               -----------------------------------------
!
! First diagnose snow temperatures and liquid
! water portion of the snow from snow heat content:
! update some snow layers parameters after new discretization 
!
DO JJ=1, SIZE(ZSNOW)
! active layers
  DO JST=1,INLVLS_USE(JJ)
      PSNOWSWE(JJ,JST)     = PSNOWDZ(JJ,JST)*PSNOWRHO(JJ,JST)
!
      ZSCAP(JJ,JST)     = SNOW3LSCAP(PSNOWRHO(JJ,JST))
!
      ZSNOWTEMP(JJ,JST) = XTT + ( ((PSNOWHEAT(JJ,JST)/PSNOWDZ(JJ,JST))   &
                + XLMTT*PSNOWRHO(JJ,JST))/ZSCAP(JJ,JST) ) 
!
      PSNOWLIQ(JJ,JST)  = MAX(0.0,ZSNOWTEMP(JJ,JST)-XTT)*ZSCAP(JJ,JST)*  &
                  PSNOWDZ(JJ,JST)/(XLMTT*XRHOLW) 
!
      ZSNOWTEMP(JJ,JST) = MIN(XTT,ZSNOWTEMP(JJ,JST))
  ENDDO  !  end loop active snow layers
! unactive layers  
  DO JST=INLVLS_USE(JJ)+1,SIZE(PSNOWSWE,2)
      PSNOWSWE(JJ,JST)=0.0
      PSNOWRHO(JJ,JST)=999.
      PSNOWDZ(JJ,JST)=0.
      PSNOWGRAN1(JJ,JST)=0.
      PSNOWGRAN2(JJ,JST)=0.
      PSNOWHIST(JJ,JST)=0.
      PSNOWAGE(JJ,JST)=0.
      PSNOWHEAT(JJ,JST)=0.
      ZSNOWTEMP(JJ,JST)=XTT
      PSNOWLIQ(JJ,JST)=0.
  ENDDO  !  end loop unactive snow layers
ENDDO    ! end loop grid points
!
!
!        4.BIS   Snow metamorphism
!                ----------------- 
!        
CALL SNOWCROMETAMO(PSNOWDZ,PSNOWGRAN1,PSNOWGRAN2,PSNOWHIST,ZSNOWTEMP,         &
                       PSNOWLIQ,PTSTEP,PSNOWSWE,INLVLS_USE) 
!

!*       5.     Snow Compaction
!               ---------------
! Calculate snow density: compaction/aging: density increases
!
CALL SNOWCROCOMPACTN(PTSTEP,PSNOWRHO,PSNOWDZ,ZSNOWTEMP,ZSNOW,     &
            PSNOWGRAN1, PSNOWGRAN2, PSNOWHIST, PSNOWLIQ,INLVLS_USE,&
             PDIRCOSZW) 
!       
!*       5.1    Snow Compaction and Metamorphism due to snow drift
!               ---------------
IF (ODRIFT_GRAIN) CALL SNOWDRIFT(PTSTEP, PVMOD, PSNOWRHO,PSNOWDZ, ZSNOW, &
                   PSNOWGRAN1,PSNOWGRAN2,PSNOWHIST,INLVLS_USE) 
!
! Update snow heat content (J/m2) using dry density instead of total density:
!
DO JJ=1, SIZE(ZSNOW)
  DO JST=1,INLVLS_USE(JJ)
  ZSCAP(JJ,JST)     = SNOW3LSCAP(PSNOWRHO(JJ,JST)-PSNOWLIQ(JJ,JST)* &
                         XRHOLW/MAX(PSNOWDZ(JJ,JST),ZSNOWDZMIN)) 
  PSNOWHEAT(JJ,JST) = PSNOWDZ(JJ,JST)*( ZSCAP(JJ,JST)*(ZSNOWTEMP(JJ,JST)-XTT)- &
                         XLMTT*PSNOWRHO(JJ,JST) ) + XLMTT*XRHOLW*PSNOWLIQ(JJ,JST) 
  ENDDO  !  end loop snow layers
ENDDO    ! end loop grid points
!
!*       6.     Solar radiation transmission
!               -----------------------------
!
! Heat source (-sink) term due to shortwave
! radiation transmission within the snowpack:
!
CALL SNOWCRORAD(TPTIME,OGLACIER,&
                PSW_RAD,PSNOWALB,PSNOWDZ,PSNOWRHO,  &
                PALB,ZRADSINK,ZRADXS,                          &
                PSNOWGRAN1, PSNOWGRAN2, PSNOWAGE,PPS,          &
                PZENITH, PPERMSNOWFRAC,INLVLS_USE              ) 
!
!*       7.     Heat transfer and surface energy budget
!               ---------------------------------------
! Snow thermal conductivity:
!
CALL SNOWCROTHRM(PSNOWRHO,ZSCOND,ZSNOWTEMP,PPS,PSNOWLIQ,OCOND_GRAIN,OCOND_YEN)
!
! Precipitation heating term:
! Rainfall renders it's heat to the snow when it enters
! the snowpack:
!
PHPSNOW(:)     = PRR(:)*XCL*(MAX(XTT,PTA(:))-XTT)    ! (W/m2)
!
! Surface Energy Budget calculations using ISBA linearized form
! and standard ISBA turbulent transfer formulation
!
CALL SNOWCROEBUD(HSNOWRES,                                                   &
                 PPEW_A_COEF, PPEW_B_COEF,                                    &
                 PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,          &
                 ZSNOWDZMIN,                                                  &
                 PZREF,ZSNOWTEMP(:,1),PSNOWRHO(:,1),PSNOWLIQ(:,1),ZSCAP(:,1), &
                 ZSCOND(:,1),ZSCOND(:,2),                                     &
                 PUREF,PEXNS,PEXNA,PDIRCOSZW,PVMOD,                           &
                 PLW_RAD,PSW_RAD,PTA,PQA,PPS,PTSTEP,                          &
                 PSNOWDZ(:,1),PSNOWDZ(:,2),PSNOWALB,PZ0,         PZ0EFF,PZ0H, &
                 ZSFCFRZ,ZRADSINK(:,1),PHPSNOW,                               &
                 ZCT,PEMISNOW,PRHOA,ZTSTERM1,ZTSTERM2,ZRA,PCDSNOW,PCHSNOW,    &
                 ZQSAT, ZDQSAT, ZRSRA, ZVMOD_IC,                              &
                 ZPET_A_COEF_T,ZPEQ_A_COEF_T,ZPET_B_COEF_T,ZPEQ_B_COEF_T      ) 
!
! Heat transfer: simple diffusion along the thermal gradient
!
ZSNOWTEMPO1(:) = ZSNOWTEMP(:,1) ! save surface snow temperature before update
!
CALL SNOWCROSOLVT(PTSTEP,ZSNOWDZMIN,PSNOWDZ,ZSCOND,ZSCAP,PTG,                &
                  PSOILCOND,PD_G,ZRADSINK,ZCT,ZTSTERM1,ZTSTERM2,              &
                  ZPET_A_COEF_T,ZPEQ_A_COEF_T,ZPET_B_COEF_T,ZPEQ_B_COEF_T,    &
                  ZTA_IC,ZQA_IC,PGRNDFLUX, ZSNOWTEMP ,ZSNOWFLUX,              &
                  INLVLS_USE                                                  ) 
!
!
!*       8.     Surface fluxes
!               --------------
!
CALL SNOWCROFLUX(ZSNOWTEMP(:,1),PSNOWDZ(:,1),PEXNS,PEXNA,           &
                 PCDSNOW,PVMOD,ZVMOD_IC,                             &
                 PTSTEP,PSNOWALB,PSW_RAD,PEMISNOW,ZLWUPSNOW,PLW_RAD, &
                 ZTA_IC,ZSFCFRZ,ZQA_IC,PHPSNOW,                      &
                 ZSNOWTEMPO1,ZSNOWFLUX,ZCT,ZRADSINK(:,1),            &
                 ZQSAT,ZDQSAT,ZRSRA,                                 &
                 PRNSNOW,PHSNOW,PGFLUXSNOW,PLES3L,PLEL3L,PEVAP,      &
                 PUSTAR                                     ) 
!
!*       9.     Snow melt
!               ---------
!
! First Test to see if snow pack vanishes during this time step:
!
CALL SNOWCROGONE(PTSTEP,PLEL3L,PLES3L,PSNOWRHO,                     &
                 PSNOWHEAT,ZRADSINK,PEVAPCOR,PTHRUFAL,PGRNDFLUX, &
                 PGFLUXSNOW,PSNOWDZ,PSNOWLIQ,ZSNOWTEMP,ZRADXS,    &
                 PRR,INLVLS_USE) 
!
! Add radiation not absorbed by snow to soil/vegetation interface flux
! (for thin snowpacks):
!
PGRNDFLUX(:) = PGRNDFLUX(:) + ZRADXS(:)
!
! Second Test to see if one or several snow layers vanishe during this time
! step. In such a case, the concerned snow layers are agregated to neighbours

 CALL SNOWCROLAYER_GONE(PTSTEP,ZSCAP,ZSNOWTEMP,PSNOWDZ,         &
                         PSNOWRHO,PSNOWLIQ,PSNOWGRAN1,PSNOWGRAN2,   &
                         PSNOWHIST,PSNOWAGE,PLES3L, INLVLS_USE) 
!
!
! For partial melt: transform excess heat content into snow liquid:
!
CALL SNOWCROMELT(ZSCAP,ZSNOWTEMP,PSNOWDZ,PSNOWRHO,PSNOWLIQ,INLVLS_USE)
!                
!
!*      10.     Snow water flow and refreezing
!               ------------------------------
! Liquid water vertical transfer and possible snowpack runoff
! And refreezing/freezing of meltwater/rainfall (ripening of the snow)
!
CALL SNOWCROREFRZ(PTSTEP,PRR,PSNOWRHO,ZSNOWTEMP,PSNOWDZ,PSNOWLIQ,PTHRUFAL, &
                  ZSCAP,PLEL3L,INLVLS_USE) 
!
!
!*      11.     Snow Evaporation/Sublimation mass updates:
!               ------------------------------------------
!
CALL SNOWCROEVAPN(PLES3L,PTSTEP,ZSNOWTEMP(:,1),PSNOWRHO(:,1),  &
                  PSNOWDZ(:,1),PEVAPCOR,PSNOWHMASS) 
!
! If all snow in uppermost layer evaporates/sublimates, re-distribute
! grid (below could be evoked for vanishingly thin snowpacks):
!
CALL SNOWCROEVAPGONE(PSNOWHEAT,PSNOWDZ,PSNOWRHO,ZSNOWTEMP,PSNOWLIQ,PSNOWGRAN1,&
                     PSNOWGRAN2,PSNOWHIST,PSNOWAGE,INLVLS_USE) 
!
!*      12.     Update surface albedo:
!               ----------------------
! Snow clear sky albedo:
!
CALL SNOWCROALB(TPTIME,OGLACIER,&
                PSNOWALB,PSNOWDZ(:,1),PSNOWRHO(:,1),   &
                PPERMSNOWFRAC,PSNOWGRAN1(:,1),PSNOWGRAN2(:,1), &
                PSNOWAGE(:,1),PSNOWGRAN1(:,2),PSNOWGRAN2(:,2),PSNOWAGE(:,2),&
                PPS, PZENITH, INLVLS_USE) 
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
DO JJ=1, SIZE(ZSNOW)
! active layers
  DO JST=1,INLVLS_USE(JJ)
    ZWHOLDMAX(JJ,JST) = SNOWCROHOLD(PSNOWRHO(JJ,JST),PSNOWLIQ(JJ,JST),PSNOWDZ(JJ,JST))
    ZLIQHEATXS(JJ)    = MAX(0.0, PSNOWLIQ(JJ,JST)*XRHOLW - &
                            ZWHOLDMAX(JJ,JST)*XRHOLW) *XLMTT/PTSTEP 
    PSNOWLIQ(JJ,JST)  = PSNOWLIQ(JJ,JST) -ZLIQHEATXS(JJ)*PTSTEP/(XRHOLW*XLMTT)
    PSNOWLIQ(JJ,JST)  = MAX(0.0, PSNOWLIQ(JJ,JST))
    PGRNDFLUX(JJ)     = PGRNDFLUX(JJ)   + ZLIQHEATXS(JJ)
    PSNOWTEMP(JJ,JST) = ZSNOWTEMP(JJ,JST)
!   Heat content using total density
    ZSCAP(JJ,JST)     = SNOW3LSCAP(PSNOWRHO(JJ,JST))
    PSNOWHEAT(JJ,JST) = PSNOWDZ(JJ,JST)*( ZSCAP(JJ,JST)*(PSNOWTEMP(JJ,JST)-XTT)&
                         - XLMTT*PSNOWRHO(JJ,JST)) +XLMTT*XRHOLW*PSNOWLIQ(JJ,JST) 
!
    PSNOWSWE(JJ,JST)  = PSNOWDZ(JJ,JST)*PSNOWRHO(JJ,JST)
  ENDDO  !  end loop active snow layers
!  
! unactive layers
  DO JST=INLVLS_USE(JJ)+1,SIZE(PSNOWSWE,2)
    PSNOWSWE(JJ,JST)=0.
    PSNOWHEAT(JJ,JST) =0.
    PSNOWRHO(JJ,JST)=999.
    PSNOWTEMP(JJ,JST)=0.
    PSNOWDZ(JJ,JST)=0.
  ENDDO  !  end loop unactive snow layers
!  
ENDDO    ! end loop grid points
!
! print some final diagnostics
IF (INLVLS_USE(I)>0) THEN
 write(*,FMT="(I4,2I4,F4.0,A7,F5.2,A10,F7.1,A11,F6.2,A13,F6.2)")  &
     TPTIME%TDATE%YEAR,TPTIME%TDATE%MONTH,TPTIME%TDATE%DAY,TPTIME%TIME/3600.,&
          'HTN= ',SUM(PSNOWDZ(I,1:INLVLS_USE(I))), 'FLUX Sol=', PGRNDFLUX(I),&
  'Tsurf_sol=',PTG(I)-273.16, 'Tbase_neige=',PSNOWTEMP(I,INLVLS_USE(I))-273.16 
ENDIF


! check suspect low temperature
DO JJ=1,SIZE(ZSNOW)
!IF(INLVLS_USE(JJ)>0) write(*,*) 'SOL:',JJ,INLVLS_USE(JJ),PGRNDFLUX(JJ),PTG(JJ),&
! PSNOWTEMP(jj,INLVLS_USE(JJ)),PSNOWTEMP(jj,1),PZENITH(JJ)           
       do JST=1,INLVLS_USE(JJ)        
       if (PSNOWTEMP(JJ,JST) < 100.) then
           write(*,*) 'pb tempe snow :',JJ,JST,PSNOWTEMP(JJ,JST)
           write(*,*) 'XLAT=',PXLAT(JJ),'XLON=',PXLON(JJ)
           write(*,*) 'INLVLS_USE(JJ):',INLVLS_USE(JJ)
           write(*,*) PSNOWDZ(JJ,1:INLVLS_USE(JJ))
           write(*,*) PSNOWRHO(JJ,1:INLVLS_USE(JJ))
           write(*,*) PSNOWTEMP(JJ,1:INLVLS_USE(JJ))
           CALL ABOR1_SFX('SNOWCRO: erreur tempe snow')
        endif
        enddo
ENDDO!
!
IF (LHOOK) CALL DR_HOOK('SNOWCRO',1,ZHOOK_HANDLE)
CONTAINS
!
!
!
!####################################################################
!####################################################################
!####################################################################
        SUBROUTINE SNOWCROCOMPACTN(PTSTEP,PSNOWRHO,PSNOWDZ,          &
                    PSNOWTEMP,PSNOW,PSNOWGRAN1,PSNOWGRAN2,PSNOWHIST, &
                    PSNOWLIQ,INLVLS_USE,PDIRCOSZW) 
!
!!    PURPOSE
!!    -------
!     Snow compaction due to overburden and settling.
!     Mass is unchanged: layer thickness is reduced
!     in proportion to density increases. Method
!     directly inherited from Crocus v2.4 and
!     coarsely described in Brun et al., J. Glac 1989 and 1992
!
!     de/e = -sigma/eta * dt
!
!     where e is layer thickness, sigma is the vertical stress, dt is the
!     time step and eta is the snow viscosity
!     * sigma is currently calculated taking into account only the overburden
!     (a term representing "metamorphism stress" in fresh snow may be added
!      in the future)
!     * eta is computed as a function of snowtype, density and temperature
!
!     The local slope is taken into account, through the variable PDIRCOSZW
!     which is directly the cosine of the local slope
!
!
!     HISTORY:
!     Basic structure from ISBA-ES model (Boone and Etchevers, 2001)
!     Implementation of Crocus laws : E. Brun, S. Morin, J.-M. Willemet July 2010.
!     Implementation of slope effect on settling : V. Vionnet, S. Morin May 2011
!
!
USE MODD_CSTS,     ONLY : XTT, XG
USE MODD_SNOW_PAR, ONLY : XRHOSMAX_ES
USE MODD_SNOW_METAMO
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP       ! Time step UNIT : s
REAL, DIMENSION(:), INTENT(IN)      :: PDIRCOSZW    ! cosine of local slope
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWTEMP    ! Snow temperature UNIT : K
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWRHO, PSNOWDZ   ! Density UNIT : kg m-3, Layer thickness UNIT : m
!
REAL, DIMENSION(:), INTENT(OUT)     :: PSNOW        ! Snowheight UNIT : m
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWGRAN1, PSNOWGRAN2, PSNOWHIST,     &!Snowtype variables
                                        PSNOWLIQ     ! Snow liquid water content UNIT ??? 
INTEGER, DIMENSION(:), INTENT(IN)   :: INLVLS_USE   ! Number of snow layers used

!
!*      0.2    declarations of local variables
!
INTEGER                             :: JJ,JST   ! looping indexes
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWRHO2,    &! work snow density UNIT : kg m-3
                                                       ZVISCOSITY,   &! Snow viscosity UNIT : N m-2 s (= Pa s)
                                                       ZSMASS !, &	     ! overburden mass for a given layer UNIT : kg m-2 
!                                                      ZWSNOWDZ       ! mass of each snow layer UNIT : kg m-2
!
! Compaction/Settling Coefficients from Crocus v2.4
!
REAL, PARAMETER                      ::  VVISC1=7.62237E6   ! pre-exponential viscosity factor (UNIT : N m-2 s)
REAL, PARAMETER                      ::  VVISC3=0.023       ! density adjustement in the exponential correction for viscosity (UNIT : m3 kg-1)
REAL, PARAMETER                      ::  VVISC4=.1          ! temperature adjustement in the exponential correction for viscosity (UNIT : K-1)
REAL, PARAMETER                      ::  VVISC5=1.          ! factor for viscosity adjustement to grain type - to be checked
REAL, PARAMETER                      ::  VVISC6=60.         ! factor for viscosity adjustement to grain type - to be checked 
!							      (especially this one ; inconsistency with Crocus v2.4)
REAL, PARAMETER                      ::  VVISC7=10.         ! factor for viscosity adjustement to grain type - to be checked
REAL, PARAMETER                      ::  VRO11=250.  ! normalization term for density dependence of the viscosity calculation (UNIT : kg m-3)
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! 1. Cumulative snow mass (kg/m2):
! --------------------------------
!

IF (LHOOK) CALL DR_HOOK('SNOWCROCOMPACTN',0,ZHOOK_HANDLE)
DO JJ=1,SIZE(PSNOW) 
     ZSMASS(JJ,1) = 0.0
     DO JST=2,INLVLS_USE(JJ)
        ZSMASS(JJ,JST) = ZSMASS(JJ,JST-1) + PSNOWDZ(JJ,JST-1) * PSNOWRHO(JJ,JST-1)
     ENDDO
     ZSMASS(JJ,1) = 0.5* PSNOWDZ(JJ,1) * PSNOWRHO(JJ,1)  ! overburden of half the mass of the uppermost layer applied to itself 
ENDDO

!
! 2. Compaction/Settling
! ----------------------

DO JJ=1, SIZE(PSNOW)
  DO JST=1, INLVLS_USE(JJ)
!
! Snow viscosity basic equation (depend on temperature and density only):
!
     ZVISCOSITY(JJ,JST)=VVISC1 *EXP(VVISC3*PSNOWRHO(JJ,JST)+       &
            VVISC4*ABS(XTT-PSNOWTEMP(JJ,JST)))*PSNOWRHO(JJ,JST)/VRO11 

! Equations below apply changes to the basic viscosity value, based on snow microstructure properties

     IF (PSNOWLIQ(JJ,JST)>0.0) THEN
        ZVISCOSITY(JJ,JST)=ZVISCOSITY(JJ,JST)/        &
             (VVISC5+VVISC6*PSNOWLIQ(JJ,JST)/PSNOWDZ(JJ,JST)) 
     ENDIF

     IF(PSNOWLIQ(JJ,JST)/PSNOWDZ(JJ,JST) <= 0.01       &
            .AND.PSNOWHIST(JJ,JST)  >= NVHIS2) THEN   
         ZVISCOSITY(JJ,JST)=ZVISCOSITY(JJ,JST)*VVISC7
     ENDIF

     IF ((PSNOWGRAN1(JJ,JST)>=0..AND.PSNOWGRAN1(JJ,JST) < VGRAN6)  &
         .AND. PSNOWHIST(JJ,JST) == NVHIS1) THEN    
         ZVISCOSITY(JJ,JST)=ZVISCOSITY(JJ,JST)*   &
               MIN(4.,EXP(MIN(VDIAM1,PSNOWGRAN2(JJ,JST)-VDIAM4)/VDIAM6)) 
     ENDIF

! 
! Calculate new snow snow density: compaction from weight/over-burden
!
     ZSNOWRHO2(JJ,JST) = PSNOWRHO(JJ,JST) + PSNOWRHO(JJ,JST)*PTSTEP* &
                         ((XG*PDIRCOSZW(JJ)*ZSMASS(JJ,JST)/ZVISCOSITY(JJ,JST)))                  
!     
! Calculate new grid thickness in response to density change
!
     PSNOWDZ(JJ,JST)   = PSNOWDZ(JJ,JST)*(PSNOWRHO(JJ,JST)/ZSNOWRHO2(JJ,JST))
!
!  Update density (kg m-3):
!
     PSNOWRHO(JJ,JST)  = ZSNOWRHO2(JJ,JST)
!
!
  ENDDO    ! end loop snow layers
ENDDO      ! end loop grid points
!
!
! 3. Update total snow depth:
! -----------------------------------------------
!
DO JJ=1, SIZE(PSNOWDZ,1)
   PSNOW(JJ)=SUM(PSNOWDZ(JJ,1:INLVLS_USE(JJ)))
ENDDO
IF (LHOOK) CALL DR_HOOK('SNOWCROCOMPACTN',1,ZHOOK_HANDLE)

!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOWCROCOMPACTN


!####################################################################
!####################################################################
!####################################################################
SUBROUTINE SNOWCROMETAMO(PSNOWDZ,PSNOWGRAN1, PSNOWGRAN2, &
                               PSNOWHIST, PSNOWTEMP, PSNOWLIQ, PTSTEP, &
                                PSNOWSWE,INLVLS_USE)  
!                               

!**** *METAMO* - METAMORPHOSE DES GRAINS
!              - SNOW METAMORPHISM
!     OBJET.
!     ------
!     METAMORPHOSE DU MANTEAU NEIGEUX.
!     EVOLUTION DU TYPE DE GRAINS 
!     MISE A JOUR DES VARIABLES HISTORIQUES.
!     METAMORPHISM OF THE SNOW GRAINS,
!     HISTORICAL VARIABLES

!**   INTERFACE.
!     ----------
!     FORMALISME ADOPTE POUR LA REPRESENTATION DES GRAINS :
!     FORMALISM FOR THE REPRESENTATION OF GRAINS
!     -----------------------------------------------------


!                    1       - -1                 NEIGE FRAICHE
!                   / \      |                    -------------
!                  /   \     |  DENDRICITE        DECRITE EN TERME
!                 /     \    |  DENDRICITY        DE DENDRICITE ET
!                /       \   |                    SPHERICITE
!               2---------3  -  0                 DESCRIBED WITH
!                                                 SPHERICITY AND
!               |---------|                       DENDRICITY
!               0         1
!               SPHERICITE
!               SPHERICITY

!               4---------5  -
!               |         |  |
!               |         |  | DIAMETRE  (OU TAILLE)
!               |         |  | DIAMETER  (OR SIZE  )
!               |         |  |                 
!               |         |  |                   NEIGE NON DENDRITIQUE
!               6---------7  -                   ---------------------

!                                                SPHERICITE ET TAILLE
!                                                SPHERICITY AND SIZE 

!              LES VARIABLES DU MODELE : 
!              -------------------------
!              CAS DENDRITIQUE             CAS NON DENDRITIQUE
!  
!            SGRAN1(JST) : DENDRICITE      SGRAN1(JST) : SPHERICITE
!            SGRAN2(JST) : SPHERICITE      SGRAN2(JST) : TAILLE (EN METRE)
!                                                        SIZE

! 
!    CAS DENDRITIQUE/ DENDRITIC CASE
!    -------------------------------
!    SGRAN1(JST) VARIE DE -VGRAN1 (-99 PAR DEFAUT) (ETOILE) A 0
!  (DENDRICITY)  >D OU LA DIVISION PAR -VGRAN1 POUR OBTENIR DES VALEURS
!                 ENTRE 1 ET 0
!                 VARIES FROM -VGRAN1 (DEFAULT -99) (FRESH SNOW) TO 0
!                 DIVISION BY -VGRAN1 TO OBTAIN VALUES BETWEEN 0 AND 1

!    SGRAN2(JST) VARIE DE 0 (CAS COMPLETEMENT ANGULEUX) A VGRAN1
!  (SPHERICITY)  (99 PAR DEFAUT)
!                >D OU LA DIVISION PAR VGRAN1 POUR OBTENIR DES VALEURS
!                 ENTRE 0 ET 1
!                 VARIES FROM 0 (SPHERICITY=0) TO VGRAN1


!    CAS NON DENDRITIQUE / NON DENDRITIC CASE
!    ---------------------------------------

!    SGRAN1(JST) VARIE DE 0 (CAS COMPLETEMENT ANGULEUX) A VGRAN1
!  (SPHERICITY)  (99 PAR DEFAUT) (CAS SPHERIQUE)
!                >D OU LA DIVISION PAR VGRAN1 POUR OBTENIR DES VALEURS
!                 ENTRE 0 ET 1
!                 VARIES FROM 0 TO 99

!    SGRAN2(JST) EST SUPERIEUR A VDIAM1-SPHERICITE (3.E-4 M) ET NE FAIT QUE CROITRE
!     (SIZE)     IS GREATER THAN VDIAM1-SPHERICITE (3.E-4 M) ALWAYS INCREASE


!    EXEMPLES : POINTS CARACTERISTIQUES DE LA FIGURE
!    --------

!                 SGRAN1     SGRAN2    DENDRICITE  SPHERICITE  TAILLE
!                                      DENDRICITY  SPHERICITY  SIZE
!      --------------------------------------------------------------
!                                                               (M)
!        1        -VGRAN1    VNSPH3        1           0.5
!        2           0         0           0            0
!        3           0       VGRAN1        0            1
!        4           0       VDIAM1                     0        4.E-4 
!        5         VGRAN1    VDIAM1-VSPHE1              1        3.E-4
!        6           0         --                       0        --
!        7         VGRAN1      --                       1        --

!     PAR DEFAUT : VGRAN1 =99   VNSPH3=50 VSPHE1=1. VDIAM1=4.E-4


!     METHODE.
!     --------
!     EVOLUTION DES TYPES DE GRAINS : SELON LES LOIS DECRITES 
!     DANS BRUN ET AL (1992)
!     PLUSIEURS CAS SONT A DISTINGUER
!      1.2 NEIGE HUMIDE
!      1.3 METAMORPHOSE NEIGE SECHE
!        1.3.1 FAIBLE GRADIENT
!        1.3.2 GRADIENT MOYEN
!        1.3.3 FORT GRADIENT
!     DANS CHAQUE CAS ON SEPARE NEIGE DENDRITIQUE ET NON DENDRITIQUE
!     LE PASSAGE DENDRITIQUE => NON DENDRITIQUE SE FAIT LORSQUE 
!     SGRAN1 DEVIENT > 0
         
!     TASSEMENT : LOIS DE VISCOSITE ADAPTEE SELON LE TYPE DE GRAINS

!     VARIABLES HISTORIQUES (CAS NON DENDRITIQUE SEULEMENT)

!     MSHIST DEFAUT
!        0           CAS NORMAL
!     NVHIS1   1     GRAINS ANGULEUX
!     NVHIS2   2     GRAINS AYANT ETE EN PRESENCE D EAU LIQUIDE 
!                    MAIS N'AYANT PAS EU DE CARATERE ANGULEUX 
!     NVHIS3   3     GRAINS AYANT ETE EN PRESENCE D EAU LIQUIDE
!                    AYANT EU AUPARAVANT UN CARACTERE ANGULEUX 

!     GRAIN METAMORPHISM ACCORDING TO BRUN ET AL (1992)
!     THE DIFFERENT CASES ARE :
!     1.2 WET SNOW
!     1.3 DRY SNOW
!       1.3.1. LOW      TEMPERATURE GRADIENT
!       1.3.2. MODERATE TEMPERATURE GRADIENT
!       1.3.3. HIGH     TEMPERATURE GRADIENT
!     THE CASE OF DENTRITIC OR NON DENDRITIC SNOW IS TREATED SEPARATELY
!     THE LIMIT DENTRITIC ==> NON DENDRITIC IS REACHED WHEN SGRAN1>0

!     SNOW SETTLING : VISCOSITY DEPENDS ON THE GRAIN TYPES

!     HISTORICAL VARIABLES (NON DENDRITIC CASE)
!     MSHIST DEFAUT
!        0           CAS NORMAL
!     NVHIS1   1     FACETED CRISTAL
!     NVHIS2   2     LIQUID WATER AND NO FACETED CRISTALS BEFORE
!     NVHIS3   3     LIQUID WATER AND FACETED CRISTALS BEFORE

!     EXTERNES.
!     ---------

!     REFERENCES.
!     -----------

!     AUTEURS.
!     --------
!        ERIC BRUN ET AL. - JOURNAL OF GLACIOLOGY 1989/1992.

!     MODIFICATIONS.
!     --------------
!        08/95: YANNICK DANIELOU - CODAGE A LA NORME DOCTOR.
!        09/96: ERIC MARTIN      - CORRECTION COMMENTAIRES
!        03/06: JM Willemet      - F90 and SI units
!        08/06: JM Willemet      - new formulation for TEL (Mwat/(Mice+Mwat) instead of Mwat/Mice.
!                                  Threshold on the diameter increasing of the wet grains.
!        01/07 : JM Willemet     - CORRECTION DES COUCHES SATUREES SUBISSANT DU TASSEMENT
!                                  CORRECTION ON THE SATURATED LAYERS WHICH ARE SETTLED
   
USE MODD_SNOW_METAMO
USE MODD_CSTS, ONLY : XTT, XPI,XRHOLW
USE MODE_SNOW3L
!
      IMPLICIT NONE 
!
!     0.1 declarations of arguments  
!      
      REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWDZ, PSNOWTEMP, &
                                              PSNOWLIQ, PSNOWSWE 
!
      REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWGRAN1, PSNOWGRAN2, &
                                              PSNOWHIST              
!   
      REAL, INTENT(IN)                    :: PTSTEP    
! 
INTEGER, DIMENSION(:), INTENT(IN)      :: INLVLS_USE      
!      
!     0.2 declaration of local variables      
!     
       REAL                              ::  ZGRADT,ZTELM, ZVDENT,      &
                                               ZDENT, ZSPHE,ZVAP, &
                                               ZDANGL 
       REAL :: diam                                       
      INTEGER                           :: INLVLS
      INTEGER JST,JJ                                !Loop controls 
      REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     INITIALISATION
!     --------------

IF (LHOOK) CALL DR_HOOK('SNOWCROMETAMO',0,ZHOOK_HANDLE)
INLVLS=SIZE(PSNOWGRAN1(:,:),2)  ! total snow layers

!*    1. METAMORPHOSES DANS LES STRATES.
!        METAMORPHISM
!        -------------------------------
DO JJ=1,SIZE(PSNOWRHO,1)
DO JST=1,INLVLS_USE(JJ)

!     1.1 INITIALISATION: GRADIENT DE TEMPERATURE / TEMPERATURE GRADIENT

  IF(JST == INLVLS_USE(JJ))THEN
    ZGRADT=ABS(PSNOWTEMP(JJ,INLVLS_USE(JJ))-PSNOWTEMP(JJ,INLVLS_USE(JJ)-1))*2. &
              /(PSNOWDZ(JJ,INLVLS_USE(JJ))+PSNOWDZ(JJ,INLVLS_USE(JJ)-1)) 
  ELSEIF( JST == 1)THEN
    ZGRADT=ABS(PSNOWTEMP(JJ,2)-PSNOWTEMP(JJ,1))*2./(PSNOWDZ(JJ,1)+PSNOWDZ(JJ,2))
  ELSE
    ZGRADT=ABS(PSNOWTEMP(JJ,JST+1)-PSNOWTEMP(JJ,JST-1))   &
           /(PSNOWDZ(JJ,JST-1)*.5+PSNOWDZ(JJ,JST)+PSNOWDZ(JJ,JST+1)*.5) 
  ENDIF

  IF(PSNOWLIQ(JJ,JST) > UEPSI)THEN


!     1.2 METAMORPHOSE HUMIDE.
!         WET SNOW METAMORPHISM

    ! TENEUR EN EAU LIQUIDE / LIQUID WATER CONTENT

    ZTELM=UPOURC*PSNOWLIQ(JJ,JST)*XRHOLW/PSNOWSWE(JJ,JST)

!   VITESSE DE DIMINUTION DE LA DENDRICITE.
!   RATE OF THE DENDRICITY DECREASE
    ZVDENT=MAX(VDENT2*ZTELM**NVDENT1,VDENT1*EXP(VVAP1/XTT))  

    IF(PSNOWGRAN1(JJ,JST) < -UEPSI)THEN

!     1.2.1 CAS DENDRITIQUE/DENDRITIC CASE.


!     VARIABLES DESCRIPTIVES DE LA DENDRICITE ET LA SPHERICITE.
      ZDENT=-PSNOWGRAN1(JJ,JST)/VGRAN1
      ZSPHE=PSNOWGRAN2(JJ,JST)/VGRAN1
!     CALCUL NOUVELLE DENDRICITE ET SPHERICITE.
      ZDENT=ZDENT-ZVDENT*PTSTEP
      ZSPHE=ZSPHE+ZVDENT*PTSTEP
      IF(ZDENT <= UEPSI)THEN
!       EVOLUTION DE SGRAN1 ET SGRAN2 ET TEST PASSAGE
!       DENDRITIQUE > NON DENDRITIQUE.
        PSNOWGRAN1(JJ,JST)=MIN(VGRAN1,ZSPHE*VGRAN1)
        PSNOWGRAN2(JJ,JST)=VDIAM1-VDIAM5*MIN(ZSPHE,VSPHE1)
      ELSE
        PSNOWGRAN1(JJ,JST)=-ZDENT*VGRAN1
        PSNOWGRAN2(JJ,JST)=MIN(VGRAN1,ZSPHE*VGRAN1)
      ENDIF
 
    ELSEIF(PSNOWGRAN1(JJ,JST) < VGRAN1-UEPSI)THEN

!     1.2.2 CAS NON DENDRITIQUE NON COMPLETEMENT SPHERIQUE.
!           EVOLUTION DE LA SPHERICITE SEULEMENT.
!           NON DENDRITIC AND NOT COMPLETELY SPHERIC CASE
!           EVOLUTION OF SPHERICITY ONLY (NOT SIZE)

            ZSPHE=PSNOWGRAN1(JJ,JST)/VGRAN1
            ZSPHE=ZSPHE+ZVDENT*PTSTEP
            PSNOWGRAN1(JJ,JST)=MIN(VGRAN1,ZSPHE*VGRAN1)
    ELSE

!     1.2.3 CAS NON DENDRITIQUE ET SPHERIQUE/NON DENDRITIC AND SPHERIC.
!           EVOLUTION DE LA TAILLE SEULEMENT/EVOLUTION OF SIZE ONLY

      diam=PSNOWGRAN2(JJ,JST)
      PSNOWGRAN2(JJ,JST)=2.                                               &
                    *(3./(4.*XPI)*(4.*XPI/3.*(PSNOWGRAN2(JJ,JST)/2.)**3   &
                    +(VTAIL1+VTAIL2*ZTELM**NVDENT1)*PTSTEP))**(1./3.) 

    ENDIF

  ELSE

!  1.3 METAMORPHOSES SECHES/DRY METAMORPHISM.

    IF(ZGRADT < VGRAT1)THEN

!     1.3.1 CALCUL METAMORPHOSE  FAIBLE/LOW GRADIENT (0-5 DEG/M).


      ZVAP=EXP(VVAP1/PSNOWTEMP(JJ,JST))
      IF(PSNOWGRAN1(JJ,JST) < -UEPSI)THEN

!       1.3.1.1 CAS DENDRITIQUE / DENDRITIC CASE.

        ZDENT=-PSNOWGRAN1(JJ,JST)/VGRAN1
        ZSPHE=PSNOWGRAN2(JJ,JST)/VGRAN1
        ZDENT=ZDENT-VDENT1*ZVAP*PTSTEP
        ZSPHE=ZSPHE+VSPHE2*ZVAP*PTSTEP
        IF(ZDENT < UEPSI)THEN
!         EVOLUTION DE SGRAN1 ET SGRAN2 ET TEST PASSAGE
!         DENDRITIQUE > NON DENDRITIQUE.
          PSNOWGRAN1(JJ,JST)=MIN(VGRAN1,ZSPHE*VGRAN1)
          PSNOWGRAN2(JJ,JST)=VDIAM1-VDIAM5*MIN(VSPHE1,ZSPHE)
        ELSE
          PSNOWGRAN1(JJ,JST)=-ZDENT*VGRAN1
          PSNOWGRAN2(JJ,JST)=MIN(VGRAN1,ZSPHE*VGRAN1)
        ENDIF
      ELSE

!       1.3.1.2 CAS NON DENDRITIQUE / NON DENDRITIC CASE.

        ZSPHE=PSNOWGRAN1(JJ,JST)/VGRAN1

        IF(PSNOWHIST(JJ,JST) /= NVHIS1 .OR. PSNOWGRAN2(JJ,JST) < VDIAM2)THEN
          ZSPHE=ZSPHE+VSPHE2*ZVAP*PTSTEP
        ELSE

!         CAS HISTORIQUE=2 OU 3 ET GROS GRAINS
!         SPHERICITE LIMITEE
!         CASE HISTORY=2 OR 3 AND BIG GRAINS
!         LIMITED SPHERICITY 

          ZSPHE=ZSPHE+     &
         VSPHE2*ZVAP*PTSTEP*EXP(MIN(0.,VDIAM3-PSNOWGRAN2(JJ,JST))/VDIAM6) 
          ZSPHE=MIN(VSPHE3,ZSPHE)
        ENDIF
        PSNOWGRAN1(JJ,JST)=MIN(VGRAN1,ZSPHE*VGRAN1)
      ENDIF
    ELSEIF (ZGRADT < VGRAT2)THEN

!     1.3.2 CALCUL METAMORPHOSE GRADIENT MOYEN/MODERATE (5-15).


      ZVAP=VDENT1*EXP(VVAP1/PSNOWTEMP(JJ,JST))*(ZGRADT)**VVAP2
      IF(PSNOWGRAN1(JJ,JST) < -UEPSI)THEN

!       1.3.2.1 CAS DENDRITIQUE / DENDRITIC CASE.

        ZDENT=-PSNOWGRAN1(JJ,JST)/VGRAN1
        ZSPHE=PSNOWGRAN2(JJ,JST)/VGRAN1
        ZDENT=ZDENT-ZVAP*PTSTEP
        ZSPHE=ZSPHE-ZVAP*PTSTEP
        IF(ZDENT < UEPSI)THEN
!         EVOLUTION DE ZSGRAN1 ET ZSGRAN2 ET TEST PASSAGE
!         DENDRITIQUE > NON DENDRITIQUE.
          PSNOWGRAN1(JJ,JST)=MAX(0.,ZSPHE*VGRAN1)
          PSNOWGRAN2(JJ,JST)=VDIAM1-VDIAM5*MAX(ZSPHE,0.)
        ELSE
          PSNOWGRAN1(JJ,JST)=-ZDENT*VGRAN1
          PSNOWGRAN2(JJ,JST)=MAX(0.,ZSPHE*VGRAN1)
        ENDIF
      ELSE

!       1.3.2.2 CAS NON DENDRITIQUE / NON DENDRITIC.

        ZSPHE=PSNOWGRAN1(JJ,JST)/VGRAN1
        ZSPHE=ZSPHE-ZVAP*PTSTEP
        PSNOWGRAN1(JJ,JST)=MAX(0.,ZSPHE*VGRAN1)
      ENDIF
    ELSE

!     1.3.3 CALCUL METAMORPHOSE FORT /HIGH GRADIENT


      ZVAP=VDENT1*EXP(VVAP1/PSNOWTEMP(JJ,JST))*(ZGRADT)**VVAP2
      IF(PSNOWGRAN1(JJ,JST) < -UEPSI)THEN

!       1.3.3.1 CAS DENDRITIQUE / DENDRITIC CASE.

        ZDENT=-PSNOWGRAN1(JJ,JST)/VGRAN1
        ZSPHE=PSNOWGRAN2(JJ,JST)/VGRAN1
        ZDENT=ZDENT-ZVAP*PTSTEP
!       CAS NON DENDRITIQUE ET ANGULEUX.
        ZSPHE=ZSPHE-ZVAP*PTSTEP
        IF(ZDENT < UEPSI)THEN
!         EVOLUTION DE SGRAN1 ET SGRAN2 ET TEST PASSAGE
!         DENDRITIQUE > NON DENDRITIQUE.
          PSNOWGRAN1(JJ,JST)=MAX(0.,ZSPHE*VGRAN1)
          PSNOWGRAN2(JJ,JST)=VDIAM1-VDIAM5*MAX(ZSPHE,0.)
        ELSE
          PSNOWGRAN1(JJ,JST)=-ZDENT*VGRAN1
          PSNOWGRAN2(JJ,JST)=MAX(0.,ZSPHE*VGRAN1)
        ENDIF
      ELSEIF(PSNOWGRAN1(JJ,JST) > UEPSI)THEN

!       1.3.3.2 CAS NON DENDRITIQUE NON COMPLETEMENT ANGULEUX.
!           NON DENDRITIC AND SPERICITY GT. 0

        ZSPHE=PSNOWGRAN1(JJ,JST)/VGRAN1
        ZSPHE=ZSPHE-ZVAP*PTSTEP
        diam=PSNOWGRAN1(JJ,JST)
        PSNOWGRAN1(JJ,JST)=MAX(0.,ZSPHE*VGRAN1)
      ELSE

!       1.3.3.3 CAS NON DENDRITIQUE ET ANGULEUX/DENDRITIC AND SPERICITY=0.
        ZDANGL = SNOW3L_MARBOUTY(PSNOWRHO(JJ,JST),PSNOWTEMP(JJ,JST),ZGRADT) 
        diam=PSNOWGRAN2(JJ,JST)
        PSNOWGRAN2(JJ,JST) = PSNOWGRAN2(JJ,JST) + ZDANGL*VFI*PTSTEP
      ENDIF
    ENDIF
  ENDIF
ENDDO
ENDDO      

!*    2. MISE A JOUR VARIABLES HISTORIQUES (CAS NON DENDRITIQUE).
!        UPDATE OF THE HISTORICAL VARIABLES
!        --------------------------------------------------------
DO JJ=1,SIZE(PSNOWRHO,1)
DO JST=1,INLVLS_USE(JJ)
  IF(PSNOWGRAN1(JJ,JST) >= 0.)THEN
    IF(PSNOWGRAN1(JJ,JST) < VSPHE4 .AND. PSNOWHIST(JJ,JST) == 0)THEN
      PSNOWHIST(JJ,JST)=NVHIS1
    ELSEIF(VGRAN1-PSNOWGRAN1(JJ,JST) < VSPHE4 .AND.        &
           PSNOWLIQ(JJ,JST)/(PSNOWDZ(JJ,JST)) > VTELV1)THEN 
      IF(PSNOWHIST(JJ,JST) == 0) PSNOWHIST(JJ,JST)=NVHIS2
      IF(PSNOWHIST(JJ,JST) == NVHIS1) PSNOWHIST(JJ,JST)=NVHIS3
    ELSEIF(PSNOWTEMP(JJ,JST) < XTT)THEN
      IF(PSNOWHIST(JJ,JST) == NVHIS2) PSNOWHIST(JJ,JST)=NVHIS4
      IF(PSNOWHIST(JJ,JST) == NVHIS3) PSNOWHIST(JJ,JST)=NVHIS5
    ENDIF
  ENDIF
ENDDO
ENDDO
IF (LHOOK) CALL DR_HOOK('SNOWCROMETAMO',1,ZHOOK_HANDLE)

END SUBROUTINE SNOWCROMETAMO  

!####################################################################
!####################################################################
!####################################################################
!  
!####################################################################
!####################################################################
!####################################################################
!
        SUBROUTINE SNOWCROALB(TPTIME,OGLACIER,&
                              PALBEDOSC,PSNOWDZ,&
                              PSNOWRHO,PPERMSNOWFRAC,   &
                              PSNOWGRAN1_TOP,PSNOWGRAN2_TOP,PSNOWAGE_TOP, &
                              PSNOWGRAN1_BOT,PSNOWGRAN2_BOT,PSNOWAGE_BOT, &
                              PPS, PZENITH, INLVLS_USE) 
!
!!    PURPOSE
!!    -------
!     Calculate the snow surface albedo. Use the method of original
!     Crocus which considers a specified spectral distribution of solar 
!     solar radiation (to be replaced by an input forcing when available)
!     In addition to original crocus, the top 2 surface snow layers are
!     considered in the calculation, using an arbitrary weighting, in order
!     to avoid time discontinuities due to layers agregation
!     Ageing depends on the presence of permanent snow cover 
!
USE MODD_CSTS,     ONLY : XDAY
USE MODD_SNOW_PAR, ONLY : XWCRN, XANSMAX, XANSMIN, XANS_TODRY, &
                           XSNOWDMIN, XANS_T, XAGLAMIN, XAGLAMAX 
USE MODD_TYPE_DATE_SURF, ONLY : DATE_TIME                                                    
!
!
USE MODE_SNOW3L
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(DATE_TIME), INTENT(IN)       :: TPTIME      ! current date and time
LOGICAL, INTENT(IN)               :: OGLACIER   ! True = Over permanent snow and ice, 
!                                                   initialise WGI=WSAT,
!                                                   Hsnow>=10m and allow 0.8<SNOALB<0.85
                                                ! False = No specific treatment
REAL, DIMENSION(:), INTENT(IN)    :: PSNOWDZ,PSNOWRHO,PPERMSNOWFRAC
!
REAL, DIMENSION(:), INTENT(INOUT) :: PALBEDOSC
!  
REAL, DIMENSION(:), INTENT(IN)  :: PSNOWGRAN1_TOP,PSNOWGRAN2_TOP,PSNOWAGE_TOP, &
                                PSNOWGRAN1_BOT,PSNOWGRAN2_BOT,PSNOWAGE_BOT, PPS 
INTEGER, DIMENSION(:), INTENT(IN)    :: INLVLS_USE                
! 
REAL, DIMENSION(:), INTENT(IN)    :: PZENITH ! solar zenith angle for future use
!
!*      0.2    declarations of local variables
!
INTEGER                             :: JJ   ! looping indexes
!
REAL, DIMENSION(SIZE(PSNOWRHO))   ::  ZANSMIN,ZANSMAX  
!
REAL, DIMENSION(SIZE(PSNOWRHO))   :: ZDIAM_TOP,ZDIAM_BOT
REAL, DIMENSION(SIZE(PSNOWRHO),3)   :: ZALB,ZALB_TOP,ZALB_BOT
!
! calibration coefficients for albedo computation
! for grains effects:
! REAL, PARAMETER :: D1=1., D2=3., D3=4., X=99.,VALB2=.96, VALB3=1.58,&
!                     VALB4=.94,VALB5=.95,VALB6=15.4,VALB7=346.3, VALB8=32.31,  &
!                     VALB9=.88, VALB10=.175,VALB11=.7,VDIOP1=2.3E-3, &
!                     VRPRE1=.5,VRPRE2=1. 
! ! for ageing effects:
! REAL, PARAMETER :: VAGING_SOIL=90. , VAGING_GLACIER=900. , VPRES1=87000.

! modifs SM 20110805 tests SIberie - albedo
REAL, PARAMETER :: D1=1., D2=3., D3=4., X=99.,VALB2=.96, VALB3=1.58,&
                   VALB4=.92,VALB5=.90,VALB6=15.4,VALB7=346.3, VALB8=32.31,&  
                   VALB9=.88, VALB10=.200,VALB11=.6,VDIOP1=2.3E-3, &
                   VRPRE1=.5,VRPRE2=1.5
! for ageing effects:
REAL, PARAMETER :: VAGING_SOIL=60., VAGING_GLACIER=900. , VPRES1=87000.
REAL, PARAMETER :: Z_SWE_ALB=5.
! end modifs SM 20110805 tests SIberie - albedo


! for spectral distribution and thickness effects
REAL, PARAMETER :: VSPEC1=.71, VSPEC2=.21 , VSPEC3=.08
! modif SM 20110519
!REAL, PARAMETER :: VSPEC1=.68, VSPEC2=.25 , VSPEC3=.07
! end modif SM 20110519
! for thickness effects
REAL, PARAMETER :: VW1=.80, VW2=.20 , VD1=.02, VD2=.01
!
REAL, DIMENSION(SIZE(PALBEDOSC))         :: NVAGE1                   
REAL            :: ZAGE_NOW
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------------
!
! 0. Initialize:
! ------------------
!
IF (LHOOK) CALL DR_HOOK('SNOWCROALB',0,ZHOOK_HANDLE)
IF(OGLACIER)THEN
   ZANSMIN(:) = XAGLAMIN * PPERMSNOWFRAC(:) + XANSMIN * (1.0-PPERMSNOWFRAC(:))
   ZANSMAX(:) = XAGLAMAX * PPERMSNOWFRAC(:) + XANSMAX * (1.0-PPERMSNOWFRAC(:))  
   NVAGE1(:)  = VAGING_GLACIER * PPERMSNOWFRAC(:) &
                  + VAGING_SOIL * (1.0-PPERMSNOWFRAC(:)) 
ELSE
   ZANSMIN(:) = XANSMIN
   ZANSMAX(:) = XANSMAX   
   NVAGE1(:)  = VAGING_SOIL
ENDIF
!  date computation for ageing effects
  CALL GREGODSTRATI(TPTIME%TDATE%YEAR,TPTIME%TDATE%MONTH,TPTIME%TDATE%DAY,   &
                      TPTIME%TIME,ZAGE_NOW) 
DO JJ=1, SIZE(PALBEDOSC)
  IF (INLVLS_USE(JJ)==0) THEN
! case with no snow on the ground 
      PALBEDOSC(JJ) = ZANSMIN(JJ)
 ELSE
! case with snow on the ground 
! top layer
    IF (PSNOWGRAN1_TOP(JJ)<0.) THEN
         ZDIAM_TOP(JJ)=-PSNOWGRAN1_TOP(JJ)*D1/X+(1.+PSNOWGRAN1_TOP(JJ)/X)* &
                   (PSNOWGRAN2_TOP(JJ)*D2/X+(1.-PSNOWGRAN2_TOP(JJ)/X)*D3) 
         ZDIAM_TOP(JJ)=ZDIAM_TOP(JJ)/10000.      
    ELSE 
         ZDIAM_TOP(JJ)=MAX(0.0003,0.5*(1.+PSNOWGRAN1_TOP(JJ)/X)* &
                          PSNOWGRAN2_TOP(JJ))       
    ENDIF
    ZALB_TOP(JJ,1)=MIN(VALB2-VALB3*SQRT(ZDIAM_TOP(JJ)),VALB4)
    ZALB_TOP(JJ,2)=MAX(0.3,VALB5-VALB6*SQRT(ZDIAM_TOP(JJ)))
    ZDIAM_TOP(JJ)=MIN(ZDIAM_TOP(JJ),VDIOP1)
    ZALB_TOP(JJ,3)=MAX(0.,VALB7*ZDIAM_TOP(JJ)-VALB8*SQRT(ZDIAM_TOP(JJ))+VALB9)
!
! AGE CORRECTION ONLY FOR VISIBLE BAND
!
    ZALB_TOP(JJ,1)=MAX(VALB11,ZALB_TOP(JJ,1)-MIN(MAX(PPS(JJ)/VPRES1,VRPRE1), &
                       VRPRE2)*VALB10*(ZAGE_NOW-PSNOWAGE_TOP(JJ))/NVAGE1(JJ))                      
     IF (INLVLS_USE(JJ)>=1) THEN
! second surface layer when it exists 
     
       IF(PSNOWGRAN1_BOT(JJ)<0.) THEN
         ZDIAM_BOT(JJ)=-PSNOWGRAN1_BOT(JJ)*D1/X+(1.+PSNOWGRAN1_BOT(JJ)/X)* &
                   (PSNOWGRAN2_BOT(JJ)*D2/X+(1.-PSNOWGRAN2_BOT(JJ)/X)*D3) 
         ZDIAM_BOT(JJ)=ZDIAM_BOT(JJ)/10000.      
       ELSE 
         ZDIAM_BOT(JJ)=MAX(0.0003,0.5*(1.+PSNOWGRAN1_BOT(JJ)/X)* &
                      PSNOWGRAN2_BOT(JJ))       
       ENDIF
              
       ZALB_BOT(JJ,1)=MIN(VALB2-VALB3*SQRT(ZDIAM_BOT(JJ)),VALB4)
       ZALB_BOT(JJ,2)=MAX(0.3,VALB5-VALB6*SQRT(ZDIAM_BOT(JJ)))
       ZDIAM_BOT(JJ)=MIN(ZDIAM_BOT(JJ),VDIOP1)
       ZALB_BOT(JJ,3)=MAX(0.,VALB7*ZDIAM_BOT(JJ)-VALB8*SQRT(ZDIAM_BOT(JJ))+ &
                        VALB9) 
!
! AGE CORRECTION ONLY FOR VISIBLE BAND

   ZALB_BOT(JJ,1)=MAX(VALB11,ZALB_BOT(JJ,1)-MIN(MAX(PPS(JJ)/VPRES1,VRPRE1), &
         VRPRE2)*VALB10*MIN(365.,ZAGE_NOW-PSNOWAGE_BOT(JJ))/NVAGE1(JJ))                      
      ELSE
! when it does not exist, the second surface layer gets top layer albedo   
        ZALB_BOT(JJ,1)= ZALB_TOP(JJ,1)
        ZALB_BOT(JJ,2)= ZALB_TOP(JJ,2)
        ZALB_BOT(JJ,3)= ZALB_TOP(JJ,3)
      ENDIF
! 
! computation of spectral albedo over 3 bands taking into account the respective
! depths of top layers 
!
     ZALB(JJ,1)=(VW1*min(1.,PSNOWDZ(JJ)/VD1)+ &
                 VW2*min(1.,max(0.,(PSNOWDZ(JJ)-VD1)/VD2))) *ZALB_TOP(JJ,1) + &
                 (VW1*(1.-min(1.,PSNOWDZ(JJ)/VD1))+ &
                 VW2*(1.- min(1.,max(0.,(PSNOWDZ(JJ)-VD1)/VD2))))*ZALB_BOT(JJ,1) 
     ZALB(JJ,2)=(VW1*min(1.,PSNOWDZ(JJ)/VD1)+ &
                 VW2*min(1.,max(0.,(PSNOWDZ(JJ)-VD1)/VD2))) *ZALB_TOP(JJ,2) + &
                 (VW1*(1.-min(1.,PSNOWDZ(JJ)/VD1))+ &
                 VW2*(1.- min(1.,max(0.,(PSNOWDZ(JJ)-VD1)/VD2))))*ZALB_BOT(JJ,2) 
     ZALB(JJ,3)=(VW1*min(1.,PSNOWDZ(JJ)/VD1)+ &
                 VW2*min(1.,max(0.,(PSNOWDZ(JJ)-VD1)/VD2))) *ZALB_TOP(JJ,3) + &
                 (VW1*(1.-min(1.,PSNOWDZ(JJ)/VD1))+ &
                 VW2*(1.- min(1.,max(0.,(PSNOWDZ(JJ)-VD1)/VD2))))*ZALB_BOT(JJ,3) 
!                
! arbitrarily specified spectral distribution  
! to be changed when solar radiation distribution is an input variable 
    PALBEDOSC(JJ)=VSPEC1*ZALB(JJ,1)+VSPEC2*ZALB(JJ,2)+VSPEC3*ZALB(JJ,3) 
!    
  ENDIF ! end case with snow on the ground
ENDDO ! end loop grid points
IF (LHOOK) CALL DR_HOOK('SNOWCROALB',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOWCROALB
!####################################################################
!####################################################################
!####################################################################
      SUBROUTINE SNOWCRORAD(TPTIME, OGLACIER,  &
                            PSW_RAD, PSNOWALB, PSNOWDZ,  &
                            PSNOWRHO, PALB, PRADSINK, PRADXS,        &
                            PSNOWGRAN1, PSNOWGRAN2, PSNOWAGE,PPS,    &
                            PZENITH, PPERMSNOWFRAC,INLVLS_USE        ) 
!
!!    PURPOSE
!!    -------
!     Calculate the transmission of shortwave (solar) radiation
!     through the snowpack (using a form of Beer's Law: exponential
!     decay of radiation with increasing snow depth).
!     Needs a first calculation of the albedo to stay coherent with
!     ISBA-ES ==> make sure to keep SNOWCRORAD coherent with SNOWCROALB
!
USE MODD_CSTS,     ONLY : XDAY
USE MODD_SNOW_PAR, ONLY : XWCRN, XANSMAX, XANSMIN, XANS_TODRY, &
                           XSNOWDMIN, XANS_T, XAGLAMIN, XAGLAMAX 
USE MODD_TYPE_DATE_SURF, ONLY: DATE_TIME                      
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(DATE_TIME), INTENT(IN)       :: TPTIME      ! current date and time
LOGICAL, INTENT(IN)               :: OGLACIER   ! True = Over permanent snow and ice, 
!                                                   initialise WGI=WSAT,
!                                                   Hsnow>=10m and allow 0.8<SNOALB<0.85
                                                ! False = No specific treatment
!
REAL, DIMENSION(:), INTENT(IN)      :: PSW_RAD, PSNOWALB, PALB,PPERMSNOWFRAC
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWRHO, PSNOWDZ
!
REAL, DIMENSION(:), INTENT(OUT)     :: PRADXS
!
REAL, DIMENSION(:,:), INTENT(OUT)   :: PRADSINK
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWGRAN1, PSNOWGRAN2, PSNOWAGE 
REAL, DIMENSION(:), INTENT(IN)      :: PPS 
INTEGER, DIMENSION(:), INTENT(IN)   :: INLVLS_USE
REAL, DIMENSION(:), INTENT(IN)      :: PZENITH 
!*      0.2    declarations of local variables
!
INTEGER                             :: JJ,JST   ! looping indexes
!
REAL, DIMENSION(SIZE(PSNOWRHO,1))    :: ZRADTOT, ZANSMIN, ZANSMAX
!
REAL     :: ZALB_NEW, ZPROJLAT
REAL, DIMENSION(SIZE(PSNOWRHO,2))   :: ZDIAM
REAL, DIMENSION(SIZE(PSNOWRHO,2),3)   :: ZBETA
REAL, DIMENSION(3)   :: ZALB,ZALB_TOP,ZALB_BOT,ZOPTICALPATH
! calibration coefficients for albedo computation
! for grains effects:
! REAL, PARAMETER :: D1=1., D2=3., D3=4., X=99.,VALB2=.96, VALB3=1.58,&
!                     VALB4=.94,VALB5=.95,VALB6=15.4,VALB7=346.3, VALB8=32.31,  &
!                     VALB9=.88, VALB10=.175,VALB11=.7,VDIOP1=2.3E-3, &
!                     VRPRE1=.5,VRPRE2=1. 
! ! for ageing effects:
! REAL, PARAMETER :: VAGING_SOIL=90. , VAGING_GLACIER=900. , VPRES1=87000.
! for spectral distribution and thickness effects
REAL, PARAMETER :: VSPEC1=.71, VSPEC2=.21 , VSPEC3=.08


! modifs SM 20110805 tests SIberie - albedo
REAL, PARAMETER :: D1=1., D2=3., D3=4., X=99.,VALB2=.96, VALB3=1.58,&
                   VALB4=.92,VALB5=.90,VALB6=15.4,VALB7=346.3, VALB8=32.31,&  
                   VALB9=.88, VALB10=.200,VALB11=.6,VDIOP1=2.3E-3, &
                   VRPRE1=.5,VRPRE2=1.5
! for ageing effects:
REAL, PARAMETER :: VAGING_SOIL=60., VAGING_GLACIER=900. , VPRES1=87000.
REAL, PARAMETER :: Z_SWE_ALB=5.
REAL  :: Z_RATE_SURF 
! end modifs SM 20110805 tests SIberie - albedo

! for thickness effects
REAL, PARAMETER :: VW1=.80, VW2=.20 , VD1=.02, VD2=.01
!
! calibration coefficients for exctinction computation
REAL, PARAMETER ::  VBETA1=1.92E-3,VBETA2=40.,VBETA3=1.098E-2,            &
                     VBETA4=100., VBETA5=2000. 
REAL, DIMENSION(SIZE(PALB))         :: NVAGE1                  
REAL                              :: ZAGE_NOW
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
! 0. Initialization:
! ------------------
!
  IF (LHOOK) CALL DR_HOOK('SNOWCRORAD',0,ZHOOK_HANDLE)
  PRADSINK(:,:)=0.
!  
IF(OGLACIER)THEN
   ZANSMIN(:) = XAGLAMIN * PPERMSNOWFRAC(:) + XANSMIN * (1.0-PPERMSNOWFRAC(:))
   ZANSMAX(:) = XAGLAMAX * PPERMSNOWFRAC(:) + XANSMAX * (1.0-PPERMSNOWFRAC(:))  
   NVAGE1(:)  = VAGING_GLACIER * PPERMSNOWFRAC(:) &
                  + VAGING_SOIL * (1.0-PPERMSNOWFRAC(:)) 
ELSE
   ZANSMIN(:) = XANSMIN
   ZANSMAX(:) = XANSMAX   
   NVAGE1(:)  = VAGING_SOIL
ENDIF
!
!
! 1. Computation of the new albedo (see SNOWCROALB):
! -----------------------------------
!
!  
  CALL GREGODSTRATI(TPTIME%TDATE%YEAR,TPTIME%TDATE%MONTH,TPTIME%TDATE%DAY,   &
                      TPTIME%TIME,ZAGE_NOW) 
DO JJ=1, SIZE(PSW_RAD)
! Coefficient for taking into account the increase of path length of rays
! in snow due to zenithal angle
ZPROJLAT=1./MAX(UEPSI,COS(PZENITH(JJ)))
!
  DO JST=1, INLVLS_USE(JJ)                     
  IF (PSNOWGRAN1(JJ,JST)<0.) THEN
     ZDIAM(JST)=-PSNOWGRAN1(JJ,JST)*D1/X+(1.+PSNOWGRAN1(JJ,JST)/X)* &
                     (PSNOWGRAN2(JJ,JST)*D2/X+(1.-PSNOWGRAN2(JJ,JST)/X)*D3) 
     ZDIAM(JST)=ZDIAM(JST)/10000.      
  ELSE 
     ZDIAM(JST)=MAX(0.0003,0.5*(1.+PSNOWGRAN1(JJ,JST)/X)*PSNOWGRAN2(JJ,JST))      
  ENDIF
  ENDDO    ! end loop snow layers
! surface layer
  ZALB_TOP(1)=MIN(VALB2-VALB3*SQRT(ZDIAM(1)),VALB4)
  ZALB_TOP(2)=MAX(0.3,VALB5-VALB6*SQRT(ZDIAM(1)))
  ZALB_TOP(3)=MAX(0.,VALB7*(MIN(ZDIAM(1),VDIOP1))-    &
             VALB8*SQRT(MIN(ZDIAM(1),VDIOP1))+VALB9) 
  ZALB_TOP(1)=MAX(VALB11,ZALB_TOP(1)-MIN(MAX(PPS(JJ)/VPRES1,VRPRE1), &
         VRPRE2)*VALB10*(ZAGE_NOW-PSNOWAGE(JJ,1))/NVAGE1(JJ)) 
! second surface layer (always exists)
  ZALB_BOT(1)=MIN(VALB2-VALB3*SQRT(ZDIAM(2)),VALB4)
  ZALB_BOT(2)=MAX(0.3,VALB5-VALB6*SQRT(ZDIAM(2)))
  ZALB_BOT(3)=MAX(0.,VALB7*(MIN(ZDIAM(2),VDIOP1))-    &
             VALB8*SQRT(MIN(ZDIAM(2),VDIOP1))+VALB9) 
  ZALB_BOT(1)=MAX(VALB11,ZALB_BOT(1)-MIN(MAX(PPS(JJ)/VPRES1,VRPRE1), &
         VRPRE2)*VALB10*MIN(365.,ZAGE_NOW-PSNOWAGE(JJ,2))/NVAGE1(JJ)) 
!
  ZALB(1)=(VW1*min(1.,PSNOWDZ(JJ,1)/VD1)+ &
               VW2*min(1.,max(0.,(PSNOWDZ(JJ,1)-VD1)/VD2))) *ZALB_TOP(1) + &
              (VW1*(1.-min(1.,PSNOWDZ(JJ,1)/VD1))+ &
               VW2*(1.- min(1.,max(0.,(PSNOWDZ(JJ,1)-VD1)/VD2))))*ZALB_BOT(1) 
  ZALB(2)=(VW1*min(1.,PSNOWDZ(JJ,1)/VD1)+ &
               VW2*min(1.,max(0.,(PSNOWDZ(JJ,1)-VD1)/VD2))) *ZALB_TOP(2) + &
              (VW1*(1.-min(1.,PSNOWDZ(JJ,1)/VD1))+ &
               VW2*(1.- min(1.,max(0.,(PSNOWDZ(JJ,1)-VD1)/VD2))))*ZALB_BOT(2) 
  ZALB(3)=(VW1*min(1.,PSNOWDZ(JJ,1)/VD1)+ &
               VW2*min(1.,max(0.,(PSNOWDZ(JJ,1)-VD1)/VD2))) *ZALB_TOP(3) + &
              (VW1*(1.-min(1.,PSNOWDZ(JJ,1)/VD1))+ &
               VW2*(1.- min(1.,max(0.,(PSNOWDZ(JJ,1)-VD1)/VD2))))*ZALB_BOT(3) 
!                
  ZALB_NEW=VSPEC1*ZALB(1)+VSPEC2*ZALB(2)+VSPEC3*ZALB(3)
  
! modif_EB test albedo ponderation basee sur SWE 5mm redondant avec ci-dessus
Z_RATE_SURF=MIN(1.,PSNOWDZ(JJ,1)*PSNOWRHO(JJ,1)/Z_SWE_ALB)
ZALB(1)=Z_RATE_SURF*ZALB_TOP(1) +(1.-Z_RATE_SURF)*ZALB_BOT(1)
ZALB(2)=Z_RATE_SURF*ZALB_TOP(2) +(1.-Z_RATE_SURF)*ZALB_BOT(2)
ZALB(3)=Z_RATE_SURF*ZALB_TOP(3) +(1.-Z_RATE_SURF)*ZALB_BOT(3)
ZALB_NEW=VSPEC1*ZALB(1)+VSPEC2*ZALB(2)+VSPEC3*ZALB(3)  
  
!
! 2. Extinction of net shortwave radiation
! ----------------------------------------
! First calculates extinction coefficients fn of grain size and density
! then calculates exctinction in the layer and increases optical path length
!
!   Initialize optical depth  
    ZOPTICALPATH(1)=0.
    ZOPTICALPATH(2)=0.
    ZOPTICALPATH(3)=0.
  DO JST=1, INLVLS_USE(JJ)                     
    ZBETA(JST,1)=MAX(VBETA1*PSNOWRHO(JJ,JST)/SQRT(ZDIAM(JST)),VBETA2)
    ZBETA(JST,2)=MAX(VBETA3*PSNOWRHO(JJ,JST)/SQRT(ZDIAM(JST)),VBETA4)
    ZBETA(JST,3)=VBETA5
    ZOPTICALPATH(1)=ZOPTICALPATH(1)+ ZBETA(JST,1)*PSNOWDZ(JJ,JST)
    ZOPTICALPATH(2)=ZOPTICALPATH(2)+ ZBETA(JST,2)*PSNOWDZ(JJ,JST)
    ZOPTICALPATH(3)=ZOPTICALPATH(3)+ ZBETA(JST,3)*PSNOWDZ(JJ,JST)
    PRADSINK(JJ,JST)  = -PSW_RAD(JJ)*(1.-PSNOWALB(JJ))/(1.-ZALB_NEW)*    &
         ( VSPEC1*(1.-ZALB(1))*EXP(-ZOPTICALPATH(1)*ZPROJLAT)  &
         + VSPEC2*(1.-ZALB(2))*EXP(-ZOPTICALPATH(2)*ZPROJLAT)  &
         + VSPEC3*(1.-ZALB(3))*EXP(-ZOPTICALPATH(3)*ZPROJLAT))       
  ENDDO    ! end loop snow layers
!
! For thin snow packs, radiation might reach base of
! snowpack...so we influence this amount with sfc albedo
! and (outside of this routine) add any excess heat
! to underlying soil:
!
  PRADSINK(JJ,INLVLS_USE(JJ)) = PRADSINK(JJ,INLVLS_USE(JJ))*(1.0-PALB(JJ))
!
! 4. Excess radiation not absorbed by snowpack (W/m2)JJ
! ----------------------------------------------------
!
  ZRADTOT(JJ)    = PRADSINK(JJ,1)+(1.-PSNOWALB(JJ))*PSW_RAD(JJ)
  DO JST=2, INLVLS_USE(JJ)                     
       ZRADTOT(JJ) = ZRADTOT(JJ) + PRADSINK(JJ,JST)-PRADSINK(JJ,JST-1)
  ENDDO    ! end loop snow layers
!
  PRADXS(JJ)     = (1.-PSNOWALB(JJ))*PSW_RAD(JJ) - ZRADTOT(JJ)
ENDDO    !end loop grid points
IF (LHOOK) CALL DR_HOOK('SNOWCRORAD',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOWCRORAD
!####################################################################
!####################################################################
!####################################################################
      SUBROUTINE SNOWCROTHRM(PSNOWRHO,PSCOND,PSNOWTEMP,PPS,PSNOWLIQ,  &
                             OCOND_GRAIN,OCOND_YEN) 
!
!!    PURPOSE
!!    -------
!     Calculate snow thermal conductivity from
!     Sun et al. 1999, J. of Geophys. Res., 104, 19587-19579
!     (vapor) and Anderson, 1976, NOAA Tech. Rep. NWS 19 (snow).
!
!     Upon activation of flag OCOND_YEN, use the Yen (1981) formula for thermal conductivity
!     This formula was originally used in Crocus.
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
LOGICAL, INTENT(IN)                 :: OCOND_GRAIN, OCOND_YEN
!
!*      0.2    declarations of local variables
!
INTEGER JJ, JST ! looping indexes
!
! ISBA-ES Thermal conductivity coefficients from Anderson (1976):
! see Boone, Meteo-France/CNRM Note de Centre No. 70 (2002)
!
REAL, PARAMETER                      :: ZSNOWTHRMCOND1 = 0.02    ! [W/(m K)]
REAL, PARAMETER                      :: ZSNOWTHRMCOND2 = 2.5E-6  ! [W m5/(kg2 K)]
!
! ISBA-ES Thermal conductivity: Implicit vapor diffn effects
! (sig only for new snow OR high altitudes)
! from Sun et al. (1999): based on data from Jordan (1991)
! see Boone, Meteo-France/CNRM Note de Centre No. 70 (2002)
!
REAL, PARAMETER                      :: ZSNOWTHRMCOND_AVAP  = -0.06023 ! [W/(m K)]
REAL, PARAMETER                      :: ZSNOWTHRMCOND_BVAP  = -2.5425  ! (W/m)
REAL, PARAMETER                      :: ZSNOWTHRMCOND_CVAP  = -289.99  ! (K)
!-------------------------------------------------------------------------------
! Crocus thermal conducitivity coefficient from Yen (1981)
REAL, PARAMETER                      :: VRKZ6 = 1.88
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! 1. Snow thermal conductivity
! ----------------------------
!

IF (LHOOK) CALL DR_HOOK('SNOWCROTHRM',0,ZHOOK_HANDLE)
DO JJ=1,SIZE(PSNOWRHO(:,:),1)
    DO JST=1,SIZE(PSNOWRHO(:,:),2)
!
           IF (OCOND_YEN) THEN
             PSCOND(JJ,JST) = XCONDI * EXP(VRKZ6 * LOG(PSNOWRHO(JJ,JST)/XRHOLW))
           ELSE
         PSCOND(JJ,JST) = (ZSNOWTHRMCOND1 + ZSNOWTHRMCOND2*PSNOWRHO(JJ,JST)* &
                     PSNOWRHO(JJ,JST))  +      &
                    MAX(0.0,(ZSNOWTHRMCOND_AVAP+  &
                   (ZSNOWTHRMCOND_BVAP/(PSNOWTEMP(JJ,JST)+  &
                   ZSNOWTHRMCOND_CVAP)))*(XP00/PPS(JJ))) 
           ENDIF
!
! Snow thermal conductivity is set to be above 0.04 W m-1 K-1
           IF(OCOND_GRAIN) THEN                  
               PSCOND(JJ,JST) = MAX(0.04,PSCOND(JJ,JST))
! Snow thermal conductivity is annihilated in presence of liquid water            
               IF(PSNOWLIQ(JJ,JST)>UEPSI) PSCOND(JJ,JST) =0.01*PSCOND(JJ,JST)  
           ENDIF
!
   ENDDO   !  end loop JST
ENDDO      ! end loop JST
IF (LHOOK) CALL DR_HOOK('SNOWCROTHRM',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOWCROTHRM
!####################################################################
!####################################################################
!####################################################################
      SUBROUTINE SNOWCROEBUD(HSNOWRES,                                                  &
                             PPEW_A_COEF, PPEW_B_COEF,                                   &
                             PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,         &
                             PSNOWDZMIN,                                                 &
                             PZREF,PTS,PSNOWRHO,PSNOWLIQ,PSCAP,PSCOND1,PSCOND2,          &
                             PUREF,PEXNS,PEXNA,PDIRCOSZW,PVMOD,                          &
                             PLW_RAD,PSW_RAD,PTA,PQA,PPS,PTSTEP,                         &
                             PSNOWDZ1,PSNOWDZ2,PALBT,PZ0,PZ0EFF,PZ0H,                    &
                             PSFCFRZ,PRADSINK,PHPSNOW,                                   &
                             PCT,PEMIST,PRHOA,PTSTERM1,PTSTERM2,PRA,PCDSNOW,PCHSNOW,     &
                             PQSAT,PDQSAT,PRSRA,PVMOD_IC,                                &
                             PPET_A_COEF_T,PPEQ_A_COEF_T,PPET_B_COEF_T,PPEQ_B_COEF_T     ) 
!
!!    PURPOSE
!!    -------
!     Calculate surface energy budget linearization (terms) and turbulent
!     exchange coefficients/resistance between surface and atmosphere.
!     (Noilhan and Planton 1989; Giordani 1993; Noilhan and Mahfouf 1996)
!
!!    MODIFICATIONS
!!    -------------
!!      Original A. Boone 
!!      Modified by E.Brun 05/11 The ratio PSFCFRZ describing liquid/solid
!                                vapor exchanges is equal to 1. as soon as
!                                there is condensation
!
USE MODD_CSTS,     ONLY : XCPD, XRHOLW, XSTEFAN, XLVTT, XLSTT, XRHOLW
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
!                                      for very stable conditions
!
REAL, DIMENSION(:), INTENT(IN)      :: PPEW_A_COEF, PPEW_B_COEF,                   &
                                        PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF,      &
                                        PPEQ_B_COEF  
!                                      PPEW_A_COEF = wind coefficient
!                                      PPEW_B_COEF = wind coefficient
!                                      PPET_A_COEF = A-air temperature coefficient
!                                      PPET_B_COEF = B-air temperature coefficient
!                                      PPEQ_A_COEF = A-air specific humidity coefficient
!                                      PPEQ_B_COEF = B-air specific humidity coefficient                                    
!
REAL, DIMENSION(:), INTENT(IN)      :: PZREF, PTS, PSNOWDZ1, PSNOWDZ2,        &
                                        PRADSINK, PSNOWRHO, PSNOWLIQ, PSCAP,   &
                                        PSCOND1, PSCOND2,                      &
                                        PZ0, PHPSNOW,                          &
                                        PALBT, PZ0EFF, PZ0H 
!
REAL, DIMENSION(:), INTENT(IN)      :: PSW_RAD, PLW_RAD, PTA, PQA, PPS, PRHOA
!
REAL, DIMENSION(:), INTENT(IN)      :: PUREF, PEXNS, PEXNA, PDIRCOSZW, PVMOD
!
REAL, DIMENSION(:), INTENT(OUT)   :: PTSTERM1, PTSTERM2, PEMIST, PRA,         &
                                      PCT, PSFCFRZ, PCDSNOW, PCHSNOW,          &
                                      PQSAT, PDQSAT, PRSRA 
!
REAL, DIMENSION(:), INTENT(OUT)   :: PVMOD_IC,                        &
                                      PPET_A_COEF_T, PPEQ_A_COEF_T,    &
                                      PPET_B_COEF_T, PPEQ_B_COEF_T 
!
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PTS))        :: ZAC, ZRI,                        &
                                      ZSCONDA, ZA, ZB, ZC,             &
                                      ZCDN, ZSNOWDZM1, ZSNOWDZM2,      &
                                      ZVMOD, ZUSTAR2, ZTS3, ZLVT,      &
                                      Z_CCOEF 
REAL, DIMENSION(SIZE(PTS))       ::  ZSNOWEVAPX                                 !   
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! 1. New saturated specific humidity and derrivative:
! ---------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SNOWCROEBUD',0,ZHOOK_HANDLE)
PQSAT(:)  = QSATI(PTS(:),PPS(:))
!
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
IF(HSNOWRES=='RIL') ZRI(:) = MIN(X_RI_MAX, ZRI(:))
!
! Surface aerodynamic resistance for heat transfers
!
CALL SURFACE_AERO_COND(ZRI, PZREF, PUREF, PVMOD, PZ0, PZ0H, ZAC, PRA, PCHSNOW)
!
PRSRA(:) = PRHOA(:) / PRA(:)
!
! For atmospheric model coupling:
!
CALL SURFACE_CD(ZRI, PZREF, PUREF, PZ0EFF, PZ0H, PCDSNOW, ZCDN)
!
!
! Modify flux-form implicit coupling coefficients:
! - wind components:
!
ZUSTAR2(:)     = (PRHOA(:)*PCDSNOW(:)*PVMOD(:)*PPEW_B_COEF(:))/        &
                  (1.0-PRHOA(:)*PCDSNOW(:)*PVMOD(:)*PPEW_A_COEF(:))  
!
ZVMOD(:)       = PPEW_A_COEF(:)*ZUSTAR2(:) + PPEW_B_COEF(:)
ZVMOD(:)       = MAX(ZVMOD(:),0.)
!
PVMOD_IC(:) =  ZVMOD(:) ! implicit wind speed: used to compute implicit wind stress
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
PSFCFRZ(:)  = 1.0 
WHERE (PSNOWLIQ(:)>UEPSI)
! computation of maximal evaporation in case  T = XTT
        ZSNOWEVAPX(:)  =  PRSRA(:) * (PQSAT(:) - PQA(:))*PTSTEP !(kg/m2)
        WHERE (ZSNOWEVAPX(:)>UEPSI)
! case when liquid water is initially present and evaporates
             PSFCFRZ(:)  = MAX(0.,1.- PSNOWLIQ(:)*XRHOLW/ZSNOWEVAPX(:))
        ELSEWHERE
! case when liquid water is initially present and liquid condensation occures 
             PSFCFRZ(:)  = 0.
        ENDWHERE        
ENDWHERE        
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
Z_CCOEF(:)       = 1.0 - PPET_A_COEF(:)*PRSRA(:)
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
IF (LHOOK) CALL DR_HOOK('SNOWCROEBUD',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOWCROEBUD
!####################################################################
!####################################################################
!####################################################################
      SUBROUTINE SNOWCROSOLVT(PTSTEP,PSNOWDZMIN,                    &
                              PSNOWDZ,PSCOND,PSCAP,PTG,              &
                              PSOILCOND,PD_G,                        &
                              PRADSINK,PCT,PTERM1,PTERM2,            &
                              PPET_A_COEF_T,PPEQ_A_COEF_T,           &
                              PPET_B_COEF_T,PPEQ_B_COEF_T,           &
                              PTA_IC, PQA_IC,                        &
                              PGBAS,PSNOWTEMP,PSNOWFLUX,             &
                              INLVLS_USE                             ) 
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
!!     MODIFICATIONS
!!     -------------
!!      Original A. Boone 
!!       05/2011: Brun  Special treatment to tackle the variable number
!!                      of snow layers
!
USE MODD_CSTS,ONLY : XTT
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP, PSNOWDZMIN
!
REAL, DIMENSION(:), INTENT(IN)      :: PTG, PSOILCOND, PD_G,        &
                                        PCT, PTERM1, PTERM2 

!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWDZ, PSCOND, PSCAP,      &
                                        PRADSINK 
!
REAL, DIMENSION(:), INTENT(IN)      :: PPET_A_COEF_T, PPEQ_A_COEF_T, &
                                        PPET_B_COEF_T, PPEQ_B_COEF_T 
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWTEMP
!
REAL, DIMENSION(:), INTENT(OUT)     :: PGBAS, PSNOWFLUX, PTA_IC, PQA_IC
!
INTEGER, DIMENSION(:), INTENT(IN)   :: INLVLS_USE
!
!*      0.2    declarations of local variables
!
!
INTEGER                               :: JJ   ! looping indexes
!
INTEGER                             :: INLVLS
!
REAL, DIMENSION(SIZE(PTG))                     :: ZGBAS, ZSNOWTEMP_DELTA
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
IF (LHOOK) CALL DR_HOOK('SNOWCROSOLVT',0,ZHOOK_HANDLE)
ZSNOWTEMP(:,:)  = PSNOWTEMP(:,:)
INLVLS          = SIZE(PSNOWDZ(:,:),2)
!
!
! 1. Calculate tri-diagnoal matrix coefficients:
! ----------------------------------------------
! For heat transfer, assume a minimum grid
! thickness (to prevent numerical
! problems for very thin snow cover):
!
DO JJ=1,SIZE(PTG)        
ZSNOWDZM(JJ,1:INLVLS_USE(JJ)) = MAX(PSNOWDZ(JJ,1:INLVLS_USE(JJ)), PSNOWDZMIN)
ZDZDIF(JJ,1:INLVLS_USE(JJ)-1) = ZSNOWDZM(JJ,1:INLVLS_USE(JJ)-1) + &
                                 ZSNOWDZM(JJ,2:INLVLS_USE(JJ)) 
ZDZDIF(JJ,INLVLS_USE(JJ))     = ZSNOWDZM(JJ,INLVLS_USE(JJ)) + PD_G(JJ)
ZWORK1(JJ,1:INLVLS_USE(JJ))   = ZSNOWDZM(JJ,1:INLVLS_USE(JJ))* &
                                 PSCOND(JJ,1:INLVLS_USE(JJ)) 
ZWORK2(JJ,1:INLVLS_USE(JJ)-1) = ZSNOWDZM(JJ,2:INLVLS_USE(JJ))*PSCOND(JJ,2:INLVLS_USE(JJ))
ZWORK2(JJ,INLVLS_USE(JJ))     = PD_G(JJ)*PSOILCOND(JJ)
!
ZDTERM(JJ,1:INLVLS_USE(JJ))   = 2.0*(ZWORK1(JJ,1:INLVLS_USE(JJ))+ &
                                ZWORK2(JJ,1:INLVLS_USE(JJ)))/ &
                                (ZDZDIF(JJ,1:INLVLS_USE(JJ))* &
                                ZDZDIF(JJ,1:INLVLS_USE(JJ))) 
!
ZCTERM(JJ,1:INLVLS_USE(JJ))   = PSCAP(JJ,1:INLVLS_USE(JJ))* &
                                 ZSNOWDZM(JJ,1:INLVLS_USE(JJ))/PTSTEP 
ENDDO
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
!
DO JJ=1,SIZE(PTG)        
! Interior Grid
ZAMTRX(JJ,2:INLVLS_USE(JJ)-1) = -ZDTERM(JJ,1:INLVLS_USE(JJ)-2)
ZBMTRX(JJ,2:INLVLS_USE(JJ)-1) =  ZCTERM(JJ,2:INLVLS_USE(JJ)-1) + &
                                      ZDTERM(JJ,1:INLVLS_USE(JJ)-2) + &
                                      ZDTERM(JJ,2:INLVLS_USE(JJ)-1) 
ZCMTRX(JJ,2:INLVLS_USE(JJ)-1) = -ZDTERM(JJ,2:INLVLS_USE(JJ)-1)
ZFRCV(JJ,2:INLVLS_USE(JJ)-1)  =  ZCTERM(JJ,2:INLVLS_USE(JJ)-1)* &
                                      PSNOWTEMP(JJ,2:INLVLS_USE(JJ)-1) - &
                                      (PRADSINK(JJ,1:INLVLS_USE(JJ)-2)- &
                                      PRADSINK(JJ,2:INLVLS_USE(JJ)-1)) 
!Lower BC
ZAMTRX(JJ,INLVLS_USE(JJ)) = -ZDTERM(JJ,INLVLS_USE(JJ)-1)
ZBMTRX(JJ,INLVLS_USE(JJ)) =  ZCTERM(JJ,INLVLS_USE(JJ)) +    &
                         ZDTERM(JJ,INLVLS_USE(JJ)-1) +ZDTERM(JJ,INLVLS_USE(JJ)) 
ZCMTRX(JJ,INLVLS_USE(JJ)) =  0.0
ZFRCV(JJ,INLVLS_USE(JJ))  =  ZCTERM(JJ,INLVLS_USE(JJ))*     &
          PSNOWTEMP(JJ,INLVLS_USE(JJ)) +  ZDTERM(JJ,INLVLS_USE(JJ))*PTG(JJ)   &
           - (PRADSINK(JJ,INLVLS_USE(JJ)-1)-PRADSINK(JJ,INLVLS_USE(JJ))) 
ENDDO
!
! - - -------------------------------------------------
!
! 4. Compute solution vector
! --------------------------
!
CALL TRIDIAG_GROUND_SNOWCRO(ZAMTRX,ZBMTRX,ZCMTRX,ZFRCV,ZSNOWTEMP,   &
            INLVLS_USE,0) 
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
ZAMTRX_M(:,1)  =  0.0
ZBMTRX_M(:,1)  =  ZCTERM(:,2) + ZDTERM(:,1) + ZDTERM(:,2)
ZCMTRX_M(:,1)  = -ZDTERM(:,2)
ZFRCV_M(:,1)   =  ZCTERM(:,2)*PSNOWTEMP(:,2) + XTT*ZDTERM(:,1)  -  &
                   (PRADSINK(:,1)-PRADSINK(:,2)) 
!
DO JJ=1,SIZE(PTG)        
   ZAMTRX_M(JJ,2:INLVLS_USE(JJ)-1)    = ZAMTRX(JJ,3:INLVLS_USE(JJ))
   ZBMTRX_M(JJ,2:INLVLS_USE(JJ)-1)    = ZBMTRX(JJ,3:INLVLS_USE(JJ))
   ZCMTRX_M(JJ,2:INLVLS_USE(JJ)-1)    = ZCMTRX(JJ,3:INLVLS_USE(JJ))
   ZFRCV_M(JJ,2:INLVLS_USE(JJ)-1)     = ZFRCV(JJ,3:INLVLS_USE(JJ))
   ZSNOWTEMP_M(JJ,2:INLVLS_USE(JJ)-1) = PSNOWTEMP(JJ,3:INLVLS_USE(JJ))
ENDDO
!
CALL TRIDIAG_GROUND_SNOWCRO(ZAMTRX_M,ZBMTRX_M,ZCMTRX_M,ZFRCV_M,ZSNOWTEMP_M,   &
                    INLVLS_USE,1) 
!
! If melting for 2 consecuative time steps, then replace current T-profile
! with one assuming T=Tf in surface layer:
!
ZSNOWTEMP_DELTA(:)    = 0.0
!
WHERE(ZSNOWTEMP(:,1) > XTT .AND. PSNOWTEMP(:,1) >= XTT)
   PSNOWFLUX(:)       = ZDTERM(:,1)*(XTT - ZSNOWTEMP_M(:,1))
   ZSNOWTEMP_DELTA(:) = 1.0
END WHERE
!
DO JJ=1,SIZE(PTG)        
   ZSNOWTEMP(JJ,2:INLVLS_USE(JJ)) = ZSNOWTEMP_DELTA(JJ)* &
                      ZSNOWTEMP_M(JJ,1:INLVLS_USE(JJ)-1)   + &
                      (1.0-ZSNOWTEMP_DELTA(JJ))*ZSNOWTEMP(JJ,2:INLVLS_USE(JJ)) 
ENDDO
!
!
! 6. Lower boundary flux:
! -----------------------
! NOTE: evaluate this term assuming the snow layer
! can't exceed the freezing point as this adjustment
! is made in melting routine. Then must adjust temperature
! to conserve energy:
!
DO JJ=1, SIZE(PTG)
 ZGBAS(JJ)            = ZDTERM(JJ,INLVLS_USE(JJ))*    &
                          (ZSNOWTEMP(JJ,INLVLS_USE(JJ))         -PTG(JJ)) 
 PGBAS(JJ)            = ZDTERM(JJ,INLVLS_USE(JJ))*    &
                          (MIN(XTT,ZSNOWTEMP(JJ,INLVLS_USE(JJ)))-PTG(JJ)) 
 ZSNOWTEMP(JJ,INLVLS_USE(JJ)) = ZSNOWTEMP(JJ,INLVLS_USE(JJ)) + &
                       (ZGBAS(JJ)-PGBAS(JJ))/ZCTERM(JJ,INLVLS_USE(JJ))      
ENDDO
!

!
! 7. Update temperatute profile in time:
! --------------------------------------
!
DO JJ=1, SIZE(PTG)
PSNOWTEMP(JJ,1:INLVLS_USE(JJ))      = ZSNOWTEMP(JJ,1:INLVLS_USE(JJ))
ENDDO
!
!
! 8. Compute new (implicit) air T and specific humidity
! -----------------------------------------------------
!
PTA_IC(:) = PPET_B_COEF_T(:) + PPET_A_COEF_T(:)* PSNOWTEMP(:,1)
!
PQA_IC(:) = PPEQ_B_COEF_T(:) + PPEQ_A_COEF_T(:)* PSNOWTEMP(:,1)
IF (LHOOK) CALL DR_HOOK('SNOWCROSOLVT',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE SNOWCROSOLVT
!####################################################################
!####################################################################
!####################################################################
      SUBROUTINE SNOWCROMELT(PSCAP,PSNOWTEMP,PSNOWDZ,         &
                             PSNOWRHO,PSNOWLIQ,INLVLS_USE     ) 
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
REAL, DIMENSION(:,:), INTENT(IN)    :: PSCAP
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWDZ, PSNOWTEMP, PSNOWRHO,   &
                                          PSNOWLIQ 
!
INTEGER, DIMENSION(:), INTENT(IN)   :: INLVLS_USE 
!
!*      0.2    declarations of local variables
!
INTEGER                             :: JJ, JST ! looping indexes
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZPHASE, ZCMPRSFACT,   &
                                                       ZSNOWLWE,   &
                                                       ZSNOWMELT, ZSNOWTEMP 
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
! 0. Initialize:
! ---------------------------
!
IF (LHOOK) CALL DR_HOOK('SNOWCROMELT',0,ZHOOK_HANDLE)
DO JJ=1,SIZE(PSNOWDZ,1)
   DO JST=1,INLVLS_USE(JJ)
ZPHASE(JJ,JST)     = 0.0
ZCMPRSFACT(JJ,JST) = 0.0
ZSNOWLWE(JJ,JST)   = 0.0
ZSNOWMELT(JJ,JST)  = 0.0
ZSNOWTEMP(JJ,JST)  = 0.0
   ENDDO
ENDDO
!
! 1. Determine amount of melt in each layer:
! ------------------------------------------
!
DO JJ=1,SIZE(PSNOWDZ,1)
   DO JST=1, INLVLS_USE(JJ)
!
! Total Liquid equivalent water content of snow (m):
!
   ZSNOWLWE(JJ,JST) = PSNOWRHO(JJ,JST)*PSNOWDZ(JJ,JST)/XRHOLW
!
! Melt snow if excess energy and snow available:
! Phase change (J/m2)
!
   ZPHASE(JJ,JST)  = MIN(PSCAP(JJ,JST)*MAX(0.0, PSNOWTEMP(JJ,JST) - XTT)*      &
                   PSNOWDZ(JJ,JST),                                       &
                   MAX(0.0,ZSNOWLWE(JJ,JST)-PSNOWLIQ(JJ,JST))*XLMTT*XRHOLW) 
!
! Update snow liq water content and temperature if melting:
! liquid water available for next layer from melting of snow
! which is assumed to be leaving the current layer (m):
!
   ZSNOWMELT(JJ,JST) = ZPHASE(JJ,JST)/(XLMTT*XRHOLW)
!
! Cool off snow layer temperature due to melt:
!
   ZSNOWTEMP(JJ,JST) = PSNOWTEMP(JJ,JST) - ZPHASE(JJ,JST)/(PSCAP(JJ,JST)*PSNOWDZ(JJ,JST))
!
! Difference with ISBA_ES: ZMELTXS should never be different of 0.
! because of the introduction of the tests in LLAYERGONE
!
   PSNOWTEMP(JJ,JST) =  ZSNOWTEMP(JJ,JST)
! The control below should be suppressed after further tests
IF (PSNOWTEMP(JJ,JST)-XTT > UEPSI) THEN
  WRITE(*,*) 'pb dans MELT PSNOWTEMP(JJ,JST) >XTT:', JJ,JST, PSNOWTEMP(JJ,JST)
  CALL ABOR1_SFX('SNOWCRO: pb dans MELT')
ENDIF
!
! Loss of snowpack depth: (m) and liquid equiv (m):
! Compression factor for melt loss: this decreases
! layer thickness and increases density thereby leaving
! total SWE constant. 
!
! Difference with ISBA_ES: All melt is considered to decrease the depth
! without consideration to the irreducible content
!
   ZCMPRSFACT(JJ,JST) = (ZSNOWLWE(JJ,JST)-(PSNOWLIQ(JJ,JST)+ZSNOWMELT(JJ,JST)))&
                       / (ZSNOWLWE(JJ,JST)-PSNOWLIQ(JJ,JST)) 
   PSNOWDZ(JJ,JST)    = PSNOWDZ(JJ,JST)*ZCMPRSFACT(JJ,JST)
   PSNOWRHO(JJ,JST)   = ZSNOWLWE(JJ,JST)*XRHOLW/PSNOWDZ(JJ,JST)
!
! 2. Add snow melt to current snow liquid water content:
! ------------------------------------------------------
!
   PSNOWLIQ(JJ,JST)   = PSNOWLIQ(JJ,JST) + ZSNOWMELT(JJ,JST)
!
   ENDDO   ! loop JST active snow layers
ENDDO   ! loop JJ grid points 
IF (LHOOK) CALL DR_HOOK('SNOWCROMELT',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE SNOWCROMELT
!####################################################################
!####################################################################
!####################################################################
      SUBROUTINE SNOWCROREFRZ(PTSTEP,PRR,                            &
                              PSNOWRHO,PSNOWTEMP,PSNOWDZ,PSNOWLIQ,   &
                              PTHRUFAL, PSCAP, PLEL3L,INLVLS_USE) 
!
!
!!    PURPOSE
!!    -------
!     Calculate any freezing/refreezing of liquid water in the snowpack.
!     Also, calculate liquid water transmission and snow runoff.
!     Refreezing causes densification of a layer.
!
!
USE MODD_CSTS,     ONLY : XTT, XLMTT, XRHOLW, XCI,XRHOLI
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
! modifs_EB layers
INTEGER, DIMENSION(:), INTENT(IN)      :: INLVLS_USE
REAL, DIMENSION(:,:), INTENT(IN)       :: PSCAP
REAL, DIMENSION(:), INTENT(IN)         :: PLEL3L

!
!*      0.2    declarations of local variables
!
INTEGER                               :: JJ, JST   ! looping indexes
!
INTEGER                               :: INLVLS     ! maximum snow layers number
!
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZPHASE,     &
                                                       ZSNOWLIQ, ZSNOWRHO,  &
                                                       ZWHOLDMAX, ZSNOWDZ,  &
                                                       ZSNOWTEMP 
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),0:SIZE(PSNOWRHO,2)) :: ZFLOWLIQ
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! 0. Initialize:
! --------------
!
IF (LHOOK) CALL DR_HOOK('SNOWCROREFRZ',0,ZHOOK_HANDLE)
INLVLS = SIZE(PSNOWDZ,2)
DO JJ=1,SIZE(PSNOWDZ,1)        
   DO JST=1,INLVLS_USE(JJ)
      ZSNOWRHO(JJ,JST)  = PSNOWRHO(JJ,JST)
      ZSNOWTEMP(JJ,JST) = PSNOWTEMP(JJ,JST)
      ZWHOLDMAX(JJ,JST) = SNOWCROHOLD(PSNOWRHO(JJ,JST),PSNOWLIQ(JJ,JST),PSNOWDZ(JJ,JST))
   ENDDO
ENDDO
!
DO JJ=1,SIZE(PSNOWDZ,1)  ! loop JJ grid points        
!
! 1. Increases Liquid Water of top layer from rain
!    ---------------------------------------------
!
!  Rainfall (m) initialises the liquid flow whih feeds the top layer 
!  and evaporation/condensation are taken into account
! 
IF(INLVLS_USE(JJ) > 0.) THEN
          ZFLOWLIQ(JJ,0)= PRR(JJ)*PTSTEP/XRHOLW 
          ZFLOWLIQ(JJ,0)= MAX(0.,ZFLOWLIQ(JJ,0)- PLEL3L(JJ)*PTSTEP/(XLVTT*XRHOLW)) 
ELSE
          ZFLOWLIQ(JJ,0)=0
ENDIF
!
!
DO JST=1,INLVLS_USE(JJ) ! loop JST active snow layers


! 2. Increases Liquid Water from the upper layers flow (or rain for top layer) 
!    -----------------------------
!
      PSNOWLIQ(JJ,JST)  = PSNOWLIQ(JJ,JST) + ZFLOWLIQ(JJ,JST-1)
!                        
! 3. Freezes liquid water in any cold layers
!    ---------------------------------------
!                        
! Calculate the maximum possible refreezing 
!
ZPHASE(JJ,JST)    = MIN(PSCAP(JJ,JST)*                                          &
                   MAX(0.0, XTT - ZSNOWTEMP(JJ,JST))*PSNOWDZ(JJ,JST),             &
                   PSNOWLIQ(JJ,JST)*XLMTT*XRHOLW) 
!
! Reduce liquid content if freezing occurs:
!
ZSNOWLIQ(JJ,JST)  = PSNOWLIQ(JJ,JST) - ZPHASE(JJ,JST)/(XLMTT*XRHOLW)       
!
! Warm layer and reduce liquid if freezing occurs:
!
ZSNOWDZ(JJ,JST)   = MAX(XSNOWDMIN/INLVLS, PSNOWDZ(JJ,JST))
!
!
! Difference with ISBA-ES: a possible cooling of current refreezing water
!                          is taken into account to calculate temperature change
!
PSNOWTEMP(JJ,JST)=XTT+(ZSNOWTEMP(JJ,JST)-XTT)* &
     (ZSNOWRHO(JJ,JST)*ZSNOWDZ(JJ,JST)-(PSNOWLIQ(JJ,JST)-ZFLOWLIQ(JJ,JST-1))*XRHOLW)/ &
     (ZSNOWRHO(JJ,JST)*ZSNOWDZ(JJ,JST)-(ZSNOWLIQ(JJ,JST)-ZFLOWLIQ(JJ,JST-1))*XRHOLW) &
     + ZPHASE(JJ,JST)/(XCI*(ZSNOWRHO(JJ,JST)*ZSNOWDZ(JJ,JST)-  &
       (ZSNOWLIQ(JJ,JST)-ZFLOWLIQ(JJ,JST-1))*XRHOLW)) 
     
!
! 4. Calculate flow from the excess of holding capacity
!    --------------------------------------------------------------
! Any water in excess of the maximum holding space for liquid water
! amount is drained into next layer down.
!
ZFLOWLIQ(JJ,JST)  = MAX(0.,ZSNOWLIQ(JJ,JST)-ZWHOLDMAX(JJ,JST))
!
ZSNOWLIQ(JJ,JST)  = ZSNOWLIQ(JJ,JST) - ZFLOWLIQ(JJ,JST)
!
! 5. Density is adjusted to conserve the mass
!    --------------------------------------------------------------
ZSNOWRHO(JJ,JST)  = (ZSNOWRHO(JJ,JST)*PSNOWDZ(JJ,JST)  &
                  -(ZFLOWLIQ(JJ,JST)-ZFLOWLIQ(JJ,JST-1))*XRHOLW)/ZSNOWDZ(JJ,JST) 
! keeps snow density below ice density
IF(ZSNOWRHO(JJ,JST)> XRHOLI) THEN
 PSNOWDZ(JJ,JST)=PSNOWDZ(JJ,JST)*ZSNOWRHO(JJ,JST)/XRHOLI
 ZSNOWRHO(JJ,JST)=XRHOLI
ENDIF

! 6. Update thickness and density and any freezing:
!    ----------------------------------------------
!
PSNOWRHO(JJ,JST) = ZSNOWRHO(JJ,JST)
PSNOWLIQ(JJ,JST) = ZSNOWLIQ(JJ,JST)
!
ENDDO ! loop JST active snow layers
!
! Any remaining throughflow after freezing is available to
! the soil for infiltration or surface runoff (m).
! I.E. This is the amount of water leaving the snowpack:
! Rate water leaves the snowpack [kg/(m2 s)]:
!
PTHRUFAL(JJ)  = PTHRUFAL(JJ) + ZFLOWLIQ(JJ,INLVLS_USE(JJ))*XRHOLW/PTSTEP
ENDDO ! loop JJ grid points
IF (LHOOK) CALL DR_HOOK('SNOWCROREFRZ',1,ZHOOK_HANDLE)
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOWCROREFRZ
!####################################################################
!####################################################################
!####################################################################
      SUBROUTINE SNOWCROFLUX(PSNOWTEMP,PSNOWDZ,PEXNS,PEXNA,         &
                             PCD,PVMOD,PVMOD_IC,                     &
                             PTSTEP,PALBT,PSW_RAD,PEMIST,PLWUPSNOW,  &
                             PLW_RAD,PTA,PSFCFRZ,PQA,PHPSNOW,        &
                             PSNOWTEMPO1,PSNOWFLUX,PCT,PRADSINK,     &
                             PQSAT,PDQSAT,PRSRA,                     &
                             PRN,PH,PGFLUX,PLES3L,PLEL3L,PEVAP,      &
                             PUSTAR                                  ) 
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
REAL, DIMENSION(:), INTENT(IN)      :: PALBT, PSW_RAD, PEMIST, PLW_RAD,      &
                                        PTA, PSFCFRZ, PQA,                    &
                                        PHPSNOW, PQSAT, PDQSAT, PRSRA,        &
                                        PCD, PVMOD, PVMOD_IC 
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PSNOWTEMP
!
REAL, DIMENSION(:), INTENT(OUT)     :: PRN, PH, PGFLUX, PLES3L, PLEL3L,      &
                                        PEVAP, PLWUPSNOW, PUSTAR 
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
IF (LHOOK) CALL DR_HOOK('SNOWCROFLUX',0,ZHOOK_HANDLE)
ZSNOWTEMP(:)  = PSNOWTEMP(:)
ZLE(:)        = 0.0
ZSMSNOW(:)    = 0.0
ZGFLUX(:)     = 0.0
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
PRN(:)       = (1. - PALBT(:)) * PSW_RAD(:) + PEMIST(:) * PLW_RAD(:) - PLWUPSNOW(:)
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
!write (*,*) 'attention test LFLUX traitement XTT supprime!'
WHERE (PSNOWTEMP > XTT .AND. PSNOWTEMPO1 < XTT)
!
   ZDELTAT(:) = XTT - PSNOWTEMPO1(:)
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
   PGFLUX(:)    = ZGFLUX(:)
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
WHERE(PSNOWTEMP(:) > XTT .AND. PSNOWTEMPO1(:) >= XTT )
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
 !write(*,*) 'Flux surface:',PGFLUX(1),PRN(1),PH(1), ZLE(1), PHPSNOW(1)
!
! 5. Friction velocity
! --------------------
!
PUSTAR(:) = SQRT( PCD(:) * PVMOD(:) * PVMOD_IC(:)) 
IF (LHOOK) CALL DR_HOOK('SNOWCROFLUX',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOWCROFLUX
!####################################################################
!####################################################################
!####################################################################
      SUBROUTINE SNOWCROEVAPN(PLES3L,PTSTEP,PSNOWTEMP, &
                              PSNOWRHO,PSNOWDZ,PEVAPCOR,    &
                              PSNOWHMASS ) 
!
!
!!    PURPOSE
!!    -------
!     Remove mass from uppermost snow layer in response to
!     evaporation (liquid) and sublimation.
!
!!     MODIFICATIONS
!!     -------------
!!      Original A. Boone 
!!      05/2011: E. Brun  Takes only into account sublimation and solid
!!                         condensation. Evaporation and liquid condensation
!!                         are taken into account in SNOWCROREFRZ
!
USE MODD_CSTS,     ONLY : XLVTT, XRHOLW, XLSTT, XLMTT, XCI, XTT
USE MODD_SNOW_PAR, ONLY : XRHOSMIN_ES
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP
!
REAL, DIMENSION(:), INTENT(IN)      :: PSNOWTEMP
!
REAL, DIMENSION(:), INTENT(IN)      :: PLES3L   ! (W/m2)
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PSNOWRHO, PSNOWDZ, PSNOWHMASS, &
                                        PEVAPCOR 
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PLES3L))       :: ZSNOWEVAPS, ZSNOWEVAP, ZSNOWEVAPX,          &
                                        ZSNOWDZ, ZEVAPCOR 
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!                                      ZEVAPCOR = for vanishingy thin snow cover,
!                                                 allow any excess evaporation
!                                                 to be extracted from the soil
!                                                 to maintain an accurate water
!                                                 balance [kg/(m2 s)]
!
!
!-------------------------------------------------------------------------------
!
! 0. Initialize:
! --------------
!
IF (LHOOK) CALL DR_HOOK('SNOWCROEVAPN',0,ZHOOK_HANDLE)
ZEVAPCOR(:)    = 0.0
ZSNOWEVAPS(:)  = 0.0
ZSNOWEVAP(:)   = 0.0
ZSNOWEVAPX(:)  = 0.0
ZSNOWDZ(:)     = 0.0
!
!
!
WHERE(PSNOWDZ > 0.0)
!
!! 1. Sublimation/condensation of snow ice
! ----------------------------------------
! Reduce layer thickness and total snow depth
! if sublimation: add to correction term if potential
! sublimation exceeds available snow cover.
!
   ZSNOWEVAPS(:)  = PLES3L(:)*PTSTEP/(XLSTT*PSNOWRHO(:))
   ZSNOWDZ(:)     = PSNOWDZ(:) - ZSNOWEVAPS(:)
   PSNOWDZ(:)     = MAX(0.0, ZSNOWDZ(:))
   ZEVAPCOR(:)    = ZEVAPCOR(:) + MAX(0.0,-ZSNOWDZ(:))*PSNOWRHO(:)/PTSTEP
!
! Total heat content change due to snowfall and sublimation (added here):
! (for budget calculations):
!
   PSNOWHMASS(:)  = PSNOWHMASS(:) - PLES3L(:)*(PTSTEP/XLSTT)*     &
                     ( XCI*(PSNOWTEMP(:)-XTT)- XLMTT) 
END WHERE
!
! 3. Update evaporation correction term:
! --------------------------------------
!
PEVAPCOR(:) = PEVAPCOR(:) + ZEVAPCOR(:)
IF (LHOOK) CALL DR_HOOK('SNOWCROEVAPN',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SNOWCROEVAPN
!####################################################################
!####################################################################
!####################################################################
SUBROUTINE SNOWCROGONE(PTSTEP,PLEL3L,PLES3L,PSNOWRHO,         &
                  PSNOWHEAT,PRADSINK_2D,PEVAPCOR,PTHRUFAL,PGRNDFLUX,    &
                  PGFLUXSNOW,PSNOWDZ,PSNOWLIQ,PSNOWTEMP,PRADXS,   &
                  PRR,INLVLS_USE) 
!
!
!!    PURPOSE
!!    -------
!     Account for the case when the last trace of snow melts
!     during a time step: ensure mass and heat balance of
!     snow AND underlying surface.
!     Original A. Boone 
!     05/2011: E. Brun  Takes into account sublimation and PGRNDFLUX 
!                       Adds rain and evaporation/liquid condensation
!                       in PTHRUFAL
!
USE MODD_CSTS,ONLY : XTT, XLSTT, XLVTT
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP
!
REAL, DIMENSION(:), INTENT(IN)      :: PLEL3L, PLES3L, PGFLUXSNOW,PRR
!  
REAL, DIMENSION(:,:), INTENT(IN)    :: PRADSINK_2D
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
INTEGER, DIMENSION(:), INTENT(INOUT)      :: INLVLS_USE
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PLES3L))       :: PRADSINK
!
INTEGER                             :: JJ
!
REAL, DIMENSION(SIZE(PLES3L))       :: ZSNOWHEATC
INTEGER, DIMENSION(SIZE(PLES3L))    :: ISNOWGONE_DELTA
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! 0. Initialize:
! --------------
!
IF (LHOOK) CALL DR_HOOK('SNOWCROGONE',0,ZHOOK_HANDLE)
PEVAPCOR(:)           = 0.0
PTHRUFAL(:)           = 0.0
!
DO JJ=1, SIZE(PRADSINK)
PRADSINK(JJ) = PRADSINK_2D(JJ,INLVLS_USE(JJ)) 
ZSNOWHEATC(JJ) = SUM(PSNOWHEAT(JJ,1:INLVLS_USE(JJ))) !total heat content (J m-2)
END DO
!
ISNOWGONE_DELTA(:)    = 1
!
! 1. Simple test to see if snow vanishes:
! ---------------------------------------
! If so, set thicknesses (and therefore mass and heat) and liquid content
! to zero, and adjust fluxes of water, evaporation and heat into underlying
! surface.
!
! takes into account the heat content corresponding to the occasional
! sublimation  and then PGRNDFLUX
!
 ZSNOWHEATC(:) = ZSNOWHEATC(:)+MAX(0.,PLES3L(:)*PTSTEP/XLSTT)* &
                  XLMTT 
!
WHERE(PGFLUXSNOW(:) + PRADSINK(:) - PGRNDFLUX(:)  >= (-ZSNOWHEATC(:)/PTSTEP) )
   PGRNDFLUX(:)       = PGFLUXSNOW(:) + (ZSNOWHEATC(:)/PTSTEP)
   PEVAPCOR(:)        = PLES3L(:)/XLSTT
   PRADXS(:)          = 0.0
   ISNOWGONE_DELTA(:) = 0          ! FLAG...if=0 then snow vanishes, else=1
END WHERE
!
! 2. Final update of snow state and computation of corresponding flow
!    Only if snow vanishes
! -----------------------------
PTHRUFAL(:) = 0.
DO JJ=1, SIZE(PRADSINK)
   IF(ISNOWGONE_DELTA(JJ) == 0 ) THEN
    PTHRUFAL(JJ) = PTHRUFAL(JJ) + SUM(PSNOWRHO(JJ,1:INLVLS_USE(JJ))* &
                                       PSNOWDZ(JJ,1:INLVLS_USE(JJ)))/PTSTEP 
! takes into account rain and condensation/evaporation
    PTHRUFAL(JJ)= PTHRUFAL(JJ) + PRR(JJ)-PLEL3L(JJ)/XLVTT
    PSNOWTEMP(JJ,:) = XTT
    PSNOWDZ(JJ,:)   = 0.
    PSNOWLIQ(JJ,:)  = 0.
    INLVLS_USE(JJ)=0
   ENDIF
ENDDO
IF (LHOOK) CALL DR_HOOK('SNOWCROGONE',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE SNOWCROGONE
!####################################################################
!####################################################################
!####################################################################
        SUBROUTINE SNOWCROEVAPGONE(PSNOWHEAT,PSNOWDZ,PSNOWRHO,PSNOWTEMP,PSNOWLIQ,&
                   PSNOWGRAN1,PSNOWGRAN2,PSNOWHIST,PSNOWAGE, INLVLS_USE) 
!
!!    PURPOSE
!!    -------
!
!     If all snow in uppermost layer evaporates/sublimates, re-distribute
!     grid (below assumes very thin snowpacks so layer-thicknesses are
!     constant).
!     Original A. Boone 
!     05/2011: E. Brun  Takes into account previous changes in the energy
!                       content
!
!
USE MODD_CSTS,     ONLY : XTT, XRHOLW, XLMTT
USE MODD_SNOW_PAR, ONLY : XRHOSMIN_ES, XSNOWDMIN, XRHOSMAX_ES
USE MODE_SNOW3L
USE MODD_SNOW_METAMO
USE MODD_TYPE_DATE_SURF
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWRHO   ! snow density profile                (kg/m3)
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWDZ    ! snow layer thickness profile        (m)
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWHEAT  ! snow heat content/enthalpy          (J/m2)
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWGRAN1 ! snow grain parameter 1              (-)
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWGRAN2 ! snow grain parameter 2              (-)
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWHIST  ! snow grain historical variable      (-)
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWAGE  ! Snow grain age
! 
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWTEMP  ! snow temperature profile            (K)
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSNOWLIQ   ! snow liquid water profile           (m)
! 
INTEGER, DIMENSION(:), INTENT(IN)      :: INLVLS_USE
!
!*      0.2    declarations of local variables
!
INTEGER                               :: JJ,JST          ! loop control
INTEGER                               :: JNLVLS
!
REAL, DIMENSION(SIZE(PSNOWDZ,1))      :: ZSNOWHEAT_1D ! total heat content                (J/m2)
REAL, DIMENSION(SIZE(PSNOWDZ,1))      :: ZSNOWRHO_1D  ! total snowpack average density    (kg/m3)
REAL, DIMENSION(SIZE(PSNOWDZ,1))      :: ZSNOW        ! total snow depth                  (m)
REAL, DIMENSION(SIZE(PSNOWDZ,1))      :: ZSCAP        ! Snow layer heat capacity          (J/K/m3)
REAL, DIMENSION(SIZE(PSNOWDZ,1))      :: ZNDENT       ! Number of dendritic layers        (-)
REAL, DIMENSION(SIZE(PSNOWDZ,1))      :: ZNVIEU       ! Number of non dendritic layers    (-)
REAL, DIMENSION(SIZE(PSNOWDZ,1))      :: ZSNOWAGE_1D  ! total snowpack average 
!age (days)
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) ::ZSNOWGRAN1N,   &
                                 ZSNOWGRAN2N,ZSNOWHISTN 
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!
!-------------------------------------------------------------------------------
!
! Initialize:
!
IF (LHOOK) CALL DR_HOOK('SNOWCROEVAPGONE',0,ZHOOK_HANDLE)
JNLVLS          = SIZE(PSNOWDZ,2)
ZSNOWHEAT_1D(:) = 0.0; ZSNOWRHO_1D(:) = 0.0; ZSNOW(:) = 0.0; ZSCAP(:) = 0.0
!
! First, determine where uppermost snow layer has completely
! evaporated/sublimated (as it becomes thin):
!
ZSNOWHEAT_1D(:) = 0.
ZSNOW(:)        = 0.
ZSNOWRHO_1D(:)  = 0.
ZNDENT(:)       = 0.
ZNVIEU(:)       = 0.  
ZSNOWAGE_1D(:)  = 0.
!           
DO JJ=1,size(PSNOWRHO,1)
   IF(PSNOWDZ(JJ,1) == 0.0)THEN   
      DO JST=2,INLVLS_USE(JJ)
         ZSNOWHEAT_1D(JJ) = ZSNOWHEAT_1D(JJ) +  PSNOWDZ(JJ,JST)* &
            (SNOW3LSCAP(PSNOWRHO(JJ,JST))*(ZSNOWTEMP(JJ,JST)-XTT) &
             - XLMTT*PSNOWRHO(JJ,JST) ) + XLMTT*XRHOLW*PSNOWLIQ(JJ,JST) 
         ZSNOW       (JJ) = ZSNOW       (JJ) + PSNOWDZ  (JJ,JST)
         ZSNOWRHO_1D (JJ) = ZSNOWRHO_1D (JJ) + &
                           PSNOWDZ  (JJ,JST)*PSNOWRHO(JJ,JST)     
         ZSNOWAGE_1D (JJ) = ZSNOWAGE_1D (JJ) + &
                           PSNOWDZ  (JJ,JST)*PSNOWRHO(JJ,JST)*PSNOWAGE(JJ,JST)   
!        snow grains
         IF(PSNOWGRAN1(JJ,JST)<-XEPSI)THEN   ! Dendritic snow
            ZNDENT(JJ) = ZNDENT(JJ)+ 1.0
         ELSE                                ! Non dendritic snow
            ZNVIEU(JJ) = ZNVIEU(JJ)+1.0
         ENDIF
      ENDDO
   ENDIF
END DO
ZSNOWRHO_1D (:) = ZSNOWRHO_1D (:) / MAX(XSNOWDMIN,ZSNOW(:))
ZSNOWAGE_1D (:)= ZSNOWAGE_1D (:)/ MAX(XSNOWDMIN,ZSNOW(:)*ZSNOWRHO_1D (:))
ZSNOWRHO_1D (:) = MAX(XRHOSMIN_ES,MIN(XRHOSMAX_ES,ZSNOWRHO_1D(:)))
!
! Where uppermost snow layer has vanished, redistribute vertical
! snow mass and heat profiles (and associated quantities):
!
CALL SNOW3LAVGRAIN(PSNOWGRAN1, PSNOWGRAN2, PSNOWHIST, &
                 ZSNOWGRAN1N, ZSNOWGRAN2N, ZSNOWHISTN,ZNDENT, ZNVIEU)         
!
DO JJ=1,SIZE(PSNOWRHO,1)
  IF(ZSNOW(JJ) /= 0.0) THEN
      PSNOWDZ(JJ,1:INLVLS_USE(JJ))   = ZSNOW(JJ)/INLVLS_USE(JJ)
      PSNOWHEAT(JJ,1:INLVLS_USE(JJ)) = ZSNOWHEAT_1D(JJ)/INLVLS_USE(JJ)
      PSNOWRHO(JJ,1:INLVLS_USE(JJ))  = ZSNOWRHO_1D(JJ)
!
      ZSCAP(JJ)        = SNOW3LSCAP(ZSNOWRHO_1D(JJ))
      PSNOWTEMP(JJ,1:INLVLS_USE(JJ)) = XTT + ( ((PSNOWHEAT(JJ,1:INLVLS_USE(JJ))/PSNOWDZ(JJ,1:INLVLS_USE(JJ)))        &
                         + XLMTT*PSNOWRHO(JJ,1:INLVLS_USE(JJ)))/ZSCAP(JJ) ) 
      PSNOWLIQ(JJ,1:INLVLS_USE(JJ))  = MAX(0.0,PSNOWTEMP(JJ,1:INLVLS_USE(JJ))-XTT)*ZSCAP(JJ)*       &
                         PSNOWDZ(JJ,1:INLVLS_USE(JJ))/(XLMTT*XRHOLW) 
      PSNOWTEMP(JJ,1:INLVLS_USE(JJ)) = MIN(XTT,PSNOWTEMP(JJ,1:INLVLS_USE(JJ)))
   ENDIF
ENDDO
IF (LHOOK) CALL DR_HOOK('SNOWCROEVAPGONE',1,ZHOOK_HANDLE)

!
END SUBROUTINE SNOWCROEVAPGONE
!
!
!

!####################################################################
!####################################################################
!####################################################################
SUBROUTINE SNOWNLFALL_UPGRID(TPTIME, OGLACIER,                 &
            PTSTEP,PSR,PTA,PVMOD,                               &
            PSNOW,PSNOWRHO,PSNOWDZ,PSNOWHEAT,PSNOWHMASS,PSNOWALB,PPERMSNOWFRAC, &
            PSNOWGRAN1,PSNOWGRAN2, GSNOWFALL,PSNOWDZN, &
            PSNOWRHOF,PSNOWDZF,PSNOWGRAN1F,PSNOWGRAN2F,PSNOWHISTF,PSNOWAGEF, &
            OMODIF_GRID,INLVLS_USE) 
!
!!    PURPOSE
!!    -------
! Adds new snowfall and updates the vertical grid in order to keep an
! optimal discertisation
!
!!    AUTHOR
!!    ------
!!      E. Brun           * Meteo-France *
!!
!
USE MODD_CSTS,     ONLY : XLMTT, XTT, XCI
USE MODD_SNOW_PAR, ONLY : XRHOSMIN_ES, XSNOWDMIN, XANSMAX,XAGLAMAX,XSNOWCRITD
USE MODD_SNOW_METAMO
!
USE MODE_SNOW3L
USE MODD_TYPE_DATE_SURF,  ONLY: DATE_TIME
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(DATE_TIME), INTENT(IN)       :: TPTIME      ! current date and time
LOGICAL, INTENT(IN)               :: OGLACIER   ! True = Over permanent snow and ice, 
!                                                   initialise WGI=WSAT,
!                                                   Hsnow>=10m and allow 0.8<SNOALB<0.85
                                                ! False = No specific treatment
!
REAL, INTENT(IN)                    :: PTSTEP
!
REAL, DIMENSION(:), INTENT(IN)      :: PSR, PTA, PVMOD, PPERMSNOWFRAC
!
REAL, DIMENSION(:), INTENT(INOUT)   :: PSNOW, PSNOWALB
! 
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWRHO, PSNOWDZ, PSNOWHEAT
!
REAL, DIMENSION(:), INTENT(OUT)     :: PSNOWHMASS
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PSNOWGRAN1, PSNOWGRAN2
!
LOGICAL, DIMENSION(:), INTENT(INOUT):: GSNOWFALL   
!
! Fresh snow characteristics
REAL, DIMENSION(:), INTENT(OUT)     :: PSNOWRHOF, PSNOWDZF 
REAL, DIMENSION(:), INTENT(OUT)     :: PSNOWGRAN1F, PSNOWGRAN2F, PSNOWHISTF 
REAL, DIMENSION(:), INTENT(OUT)     :: PSNOWAGEF  
! New vertical grid
REAL, DIMENSION(:,:), INTENT(OUT)   :: PSNOWDZN
! 
INTEGER, DIMENSION(:), INTENT(INOUT)       :: INLVLS_USE
!*      0.2    declarations of local variables
!
INTEGER                             :: INLVLS, INLVLSMIN, INLVLSMAX, JJ, JST
!
LOGICAL, DIMENSION(SIZE(PTA))       :: OAGREG_SURF,OMODIF_GRID 
!
REAL, DIMENSION(SIZE(PTA))          :: ZSNOWFALL,ZSNOWTEMP,ZSCAP,ZANSMAX
REAL                                :: ZAGE_NOW
!
! ISBA-ES CROCUS (Pahaut 1976): snowfall density coefficients:
!
REAL, PARAMETER                      :: ZSNOWFALL_A_SN = 109.0  ! kg/m3
REAL, PARAMETER                      :: ZSNOWFALL_B_SN =   6.0  ! kg/(m3 K)
REAL, PARAMETER                      :: ZSNOWFALL_C_SN =  26.0  ! kg/(m7/2 s1/2)
!
! Coefficients for the optimal vertical grid calculation
REAL, PARAMETER                      ::  DZ1=0.01
REAL, PARAMETER                      ::  DZ2=0.0125
REAL, PARAMETER                      ::  DZ3=0.015
REAL, PARAMETER                      ::  DZ3_bis=0.03
REAL, PARAMETER                      ::  DZ4=0.04
REAL, PARAMETER                      ::  DZ5=0.05
REAL, PARAMETER                      ::  DZ_base=0.02
REAL, PARAMETER                      ::  DZ_internal=0.07
REAL, PARAMETER                      ::  scale_cm=100.
!
! Coefficients for computing the difference in 2 snow layer characteristics
REAL, PARAMETER                      ::  DIFF_1=20.
REAL, PARAMETER                      ::  DIFF_MAX=200.
REAL, PARAMETER                      ::  scale_DIFF=25.
!
! Coeefficients for snow layer splitting
REAL, PARAMETER                      ::  DZMIN_TOP=0.01
REAL, PARAMETER                      ::  DZMIN_TOP_BIS=0.005
REAL, PARAMETER                      ::  DZMIN_BOT=0.02
REAL, PARAMETER                      ::  SPLIT_COEF=8.
!
! Coeefficients for snow layer agregation 
REAL, PARAMETER                      ::  AGREG_COEF_1=5.
REAL, PARAMETER                      ::  AGREG_COEF_2=4.5
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) ::ZDZOPT 

REAL        :: ZPENALTY, ZDIFTYPE_INF,ZDIFTYPE_SUP,ZCRITSIZE, &
                     ZCRITSIZE_INF, ZCRITSIZE_SUP 
!        
INTEGER JST_1, JJ_A_AGREG_SUP, JJ_A_AGREG_INF,JJ_A_DEDOUB
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*      1.0   Initialization and snowage calculation for the present date 
!
IF (LHOOK) CALL DR_HOOK('SNOWNLFALL_UPGRID',0,ZHOOK_HANDLE)
INLVLS         = SIZE (PSNOWRHO(:,:),2)
INLVLSMAX      = SIZE (PSNOWRHO(:,:),2)
INLVLSMIN      = 3 
ZSNOWTEMP(:)   = XTT
GSNOWFALL(:)   =.FALSE. 
PSNOWHMASS(:)  = 0.0 
PSNOWRHOF(:)   = 0.0 
PSNOWDZF(:)    = 0.0 
PSNOWGRAN1F(:) = 0.0 
PSNOWGRAN2F(:) = 0.0      
PSNOWHISTF(:)  = 0.0
PSNOWDZN(:,:)  = PSNOWDZ(:,:)
OAGREG_SURF(:)=.FALSE.
OMODIF_GRID(:)=.FALSE.
!
! 
  CALL GREGODSTRATI(TPTIME%TDATE%YEAR,TPTIME%TDATE%MONTH,TPTIME%TDATE%DAY,   &
                      TPTIME%TIME,ZAGE_NOW) 
!                     
!*      1.1   Calculation of the optimal vertical grid size 
!             as a function of maximum number of layers and of current 
!             snow depth
!
DO JJ=1, SIZE(PSNOW(:))
   ZDZOPT(JJ,1) = MIN(DZ1,PSNOW(JJ)/MAX(INLVLSMIN,INLVLS_USE(JJ))) 
   ZDZOPT(JJ,2) = MIN(DZ2,PSNOW(JJ)/MAX(INLVLSMIN,INLVLS_USE(JJ))) 
   ZDZOPT(JJ,3) = MIN(DZ3,PSNOW(JJ)/MAX(INLVLSMIN,INLVLS_USE(JJ))) 
   IF (INLVLS_USE(JJ)>3)  ZDZOPT(JJ,3) = MIN(DZ3_bis,PSNOW(JJ)/INLVLS_USE(JJ))
   IF (INLVLS_USE(JJ)>4)  ZDZOPT(JJ,4) = MIN(DZ4,PSNOW(JJ)/INLVLS_USE(JJ))
   IF (INLVLS_USE(JJ)>5)  ZDZOPT(JJ,5) = MIN(DZ5,PSNOW(JJ)/INLVLS_USE(JJ))
   IF (INLVLS_USE(JJ)>0) ZDZOPT(JJ,INLVLS_USE(JJ))=   &
              MIN(DZ_base,PSNOW(JJ)/MAX(INLVLSMIN,INLVLS_USE(JJ)))  
   DO JST=6,INLVLS_USE(JJ)-1,1
       ZDZOPT(JJ,JST) = MAX(DZ_internal,(PSNOW(JJ) - SUM(ZDZOPT(JJ,1:5))-    &
                         ZDZOPT(JJ,INLVLS_USE(JJ))) /(INLVLS_USE(JJ)-6)) 
   END DO
END DO
!
!*      2.0   Fresh snow characteristics
!
!
!
! Heat content of newly fallen snow (J/m2):
! NOTE for now we assume the snowfall has
! the temperature of the snow surface upon reaching the snow.
! This is done as opposed to using the air temperature since
! this flux is quite small and has little to no impact
! on the time scales of interest. If we use the above assumption
! then, then the snowfall advective heat flux is zero.
!!
DO JJ=1, SIZE(PSNOW(:)) 
IF (PSR(JJ) > UEPSI ) THEN  
! newly fallen snow characteristics:
  IF (INLVLS_USE(JJ)>0) THEN !Case of new snowfall on a previously snow-free surface 
    ZSCAP(JJ)      = XCI*PSNOWRHO(JJ,1)
    ZSNOWTEMP(JJ)  = XTT + (PSNOWHEAT(JJ,1) +                              &
                       XLMTT*PSNOWRHO(JJ,1)*PSNOWDZ(JJ,1))/                 &
                       (ZSCAP(JJ)*MAX(XSNOWDMIN/INLVLS,PSNOWDZ(JJ,1))) 
  ELSE  ! case with bare ground
    ZSNOWTEMP(JJ)=PTA(JJ)
  ENDIF
  ZSNOWTEMP(JJ)  = MIN(XTT, ZSNOWTEMP(JJ))
!
  PSNOWHMASS(JJ) = PSR(JJ)*(XCI*(ZSNOWTEMP(JJ)-XTT)-XLMTT)*PTSTEP
  PSNOWRHOF(JJ)  = MAX(XRHOSMIN_ES, ZSNOWFALL_A_SN + &
                    ZSNOWFALL_B_SN*(PTA(JJ)-XTT)+ &
                    ZSNOWFALL_C_SN*MIN(PVMOD(JJ),SQRT(PVMOD(JJ)))) 
  ZSNOWFALL(JJ)  = PSR(JJ)*PTSTEP/PSNOWRHOF(JJ)    ! snowfall thickness (m)
  PSNOW(JJ)      = PSNOW(JJ) + ZSNOWFALL(JJ)
  PSNOWDZF(JJ)   = ZSNOWFALL(JJ)
  PSNOWGRAN1F(JJ)= MAX(MIN(XNDEN1*PVMOD(JJ)-XNDEN2,XNDEN3),-XGRAN)
  PSNOWGRAN2F(JJ)= MIN(MAX(XNSPH1*PVMOD(JJ)+XNSPH2,XNSPH3),XNSPH4)     
  PSNOWHISTF(JJ) = 0.0
  PSNOWAGEF(JJ)  = ZAGE_NOW
  GSNOWFALL(JJ)  =.TRUE. 
  OMODIF_GRID(JJ)=.TRUE.
!      
ENDIF
ENDDO
!                   
!  intialize the albedo:
!  penser a changer 0.000001 par UEPSI
IF(OGLACIER)THEN
   ZANSMAX(:) = XAGLAMAX * PPERMSNOWFRAC(:) + XANSMAX * (1.0-PPERMSNOWFRAC(:))
ELSE
   ZANSMAX(:) = XANSMAX
ENDIF
WHERE(GSNOWFALL(:) .AND. ABS(PSNOW(:)-ZSNOWFALL(:))< 0.000001)
   PSNOWALB(:) = ZANSMAX(:)
END WHERE
!
! Computation of the new grid size
! It starts with successive exclusive cases
! Each case is described inside the corresponding condition
! 
! cases with fresh snow
!
DO JJ=1,SIZE(PSNOW(:)) ! grid point loop
IF((.NOT.GSNOWFALL(JJ)).AND.PSNOW(JJ)>=XSNOWCRITD.AND. &
                       INLVLS_USE(JJ)>=INLVLSMIN) THEN 
! no fresh snow + deep enough snowpack + enough snow layers ==> no change
!
ELSEIF(PSNOW(JJ)<XSNOWCRITD.OR.PSNOW(JJ)==ZSNOWFALL(JJ).OR. &
                              INLVLS_USE(JJ)<INLVLSMIN  ) THEN 
! too shallow snowpack or only fresh snow or too few layers
! ==> uniform grid and identical snow layers / number depends on snow depth 
  OMODIF_GRID(JJ)=.TRUE.
  INLVLS_USE(JJ)= MAX(INLVLSMIN,MIN(INLVLSMAX,INT(PSNOW(JJ)*scale_cm)))
  PSNOWDZN(JJ,1:INLVLS_USE(JJ))=PSNOW(JJ)/INLVLS_USE(JJ)
ELSE
! fresh snow over snow covered ground + enough snow layers 
  OMODIF_GRID(JJ)=.TRUE.
  ZDIFTYPE_SUP= SNOW3LDIFTYP(PSNOWGRAN1(JJ,1),PSNOWGRAN1F(JJ),&
                        PSNOWGRAN2(JJ,1),PSNOWGRAN2F(JJ)) 
  IF (ZDIFTYPE_SUP<DIFF_1.AND.PSNOWDZ(JJ,1)< ZDZOPT(JJ,1)) THEN
! fresh snow is agregated to the similar and shallow surface layer
       PSNOWDZN(JJ,1)=PSNOWDZ(JJ,1)+PSNOWDZF(JJ)
       DO JST=  INLVLS_USE(JJ), 2,-1
          PSNOWDZN(JJ, JST)=PSNOWDZ(JJ,JST)
       ENDDO
   ELSEIF(INLVLS_USE(JJ)<INLVLSMAX) THEN
! fresh snow is too different from the surface or the surface is too deep
! and there is room for extra layers ==> we create a new layer
       INLVLS_USE(JJ)=INLVLS_USE(JJ)+1
       PSNOWDZN(JJ,1)=PSNOWDZF(JJ)
       DO JST=  INLVLS_USE(JJ), 2,-1
          PSNOWDZN(JJ, JST)=PSNOWDZ(JJ,JST-1)
       ENDDO
   ELSE
! fresh snow is too different from the surface or the surface is too deep
! and there is no room for extra layers
! ==> we agregate internal most similar snowlayers and create a new surface layer
   JJ_A_AGREG_SUP=1
   JJ_A_AGREG_INF=2
   ZDIFTYPE_INF= SNOW3LDIFTYP(PSNOWGRAN1(JJ,1),PSNOWGRAN1F(JJ),&
                     PSNOWGRAN2(JJ,1),PSNOWGRAN2F(JJ)) 
   ZCRITSIZE_INF= scale_DIFF*(PSNOWDZ(JJ,1)/ZDZOPT(JJ,1)+ &
                     PSNOWDZ(JJ,2)/ZDZOPT(JJ,2)) 
   ZPENALTY=ZDIFTYPE_INF+ZCRITSIZE_INF
   DO JST=2, INLVLS_USE(JJ)-1
   ZDIFTYPE_SUP= SNOW3LDIFTYP(PSNOWGRAN1(JJ,JST-1),PSNOWGRAN1(JJ, JST),&
                        PSNOWGRAN2(JJ,JST-1),PSNOWGRAN2(JJ, JST)) 
   ZDIFTYPE_INF= SNOW3LDIFTYP(PSNOWGRAN1(JJ,JST+1),PSNOWGRAN1(JJ, JST),&
                        PSNOWGRAN2(JJ,JST+1),PSNOWGRAN2(JJ, JST)) 
   ZCRITSIZE_SUP= scale_DIFF*(PSNOWDZ(JJ,JST)/ZDZOPT(JJ,JST)+ &
                        PSNOWDZ(JJ,JST-1)/ZDZOPT(JJ,JST-1))  
   ZCRITSIZE_INF= scale_DIFF*(PSNOWDZ(JJ,JST)/ZDZOPT(JJ,JST)+ &
                        PSNOWDZ(JJ,JST+1)/ZDZOPT(JJ,JST+1))  
   IF(ZDIFTYPE_SUP+ZCRITSIZE_SUP < ZPENALTY) THEN
     ZPENALTY=ZDIFTYPE_SUP+ZCRITSIZE_SUP
     JJ_A_AGREG_SUP=JST-1
     JJ_A_AGREG_INF=JST
   ENDIF
   IF(ZDIFTYPE_INF+ZCRITSIZE_INF < ZPENALTY) THEN
     ZPENALTY=ZDIFTYPE_INF+ZCRITSIZE_INF
     JJ_A_AGREG_SUP=JST
     JJ_A_AGREG_INF=JST+1
   ENDIF
   ENDDO
   ZDIFTYPE_SUP= SNOW3LDIFTYP(PSNOWGRAN1(JJ,INLVLS_USE(JJ)-1),   &
        PSNOWGRAN1(JJ,INLVLS_USE(JJ)), PSNOWGRAN2(JJ,INLVLS_USE(JJ)-1),  &
        PSNOWGRAN2(JJ, INLVLS_USE(JJ))) 
   ZCRITSIZE_SUP=scale_DIFF*(PSNOWDZ(JJ,INLVLS_USE(JJ))/ZDZOPT(JJ,INLVLS_USE(JJ))+ &
                  PSNOWDZ(JJ,INLVLS_USE(JJ)-1)/ZDZOPT(JJ,INLVLS_USE(JJ)-1))  
   IF(ZDIFTYPE_SUP+ZCRITSIZE_SUP < ZPENALTY) THEN
      JJ_A_AGREG_SUP=INLVLS_USE(JJ)-1
      JJ_A_AGREG_INF=INLVLS_USE(JJ)
   ENDIF
! agregation of the similar layers and shift of upper layers
   PSNOWDZN(JJ,JJ_A_AGREG_INF)=PSNOWDZ(JJ,JJ_A_AGREG_INF) + &
                  PSNOWDZ(JJ,JJ_A_AGREG_SUP) 
   DO JST=  JJ_A_AGREG_SUP, 2,-1
     PSNOWDZN(JJ, JST)=PSNOWDZ(JJ,JST-1)
   ENDDO
     PSNOWDZN(JJ,1)=PSNOWDZF(JJ)
   ENDIF
ENDIF ! end of the case with fresh snow
ENDDO ! end loop grid points
!
! cases with no fresh snow and no previous grid resize
!
IF(INLVLSMIN==INLVLSMAX) THEN ! specific case with INLSVSMIN = INLVLSMAX  (INLVLS)
! check if surface layer depth is too small
! in such a case looks for an other layer to be split
 DO JJ=1,SIZE(PSNOW(:)) ! loop grid points
  IF(.NOT.OMODIF_GRID(JJ)) THEN 
   IF(PSNOWDZ(JJ,1)< DZMIN_TOP) THEN
     OMODIF_GRID(JJ)=.TRUE.
     OAGREG_SURF(JJ)=.TRUE.
     ZPENALTY=PSNOWDZ(JJ,2)/ZDZOPT(JJ,2)
     IF(PSNOWDZ(JJ,2)<DZMIN_TOP) ZPENALTY=0.
       JJ_A_DEDOUB=2
       DO JST=3, INLVLS
          ZCRITSIZE=PSNOWDZ(JJ,JST)/ZDZOPT(JJ,JST)
          IF(JST==INLVLS.AND.PSNOWDZ(JJ,JST)<DZMIN_BOT) ZCRITSIZE=0.
          IF (ZCRITSIZE > ZPENALTY) THEN
             ZPENALTY=ZCRITSIZE
             JJ_A_DEDOUB=JST
          ENDIF
       ENDDO
       IF(JJ_A_DEDOUB == 2) THEN
            PSNOWDZN(JJ,1)=0.5*(PSNOWDZ(JJ,1)+PSNOWDZ(JJ,2))
            PSNOWDZN(JJ,2)=0.5*(PSNOWDZ(JJ,1)+PSNOWDZ(JJ,2))
       ELSE
            PSNOWDZN(JJ,1)=(PSNOWDZ(JJ,1)+PSNOWDZ(JJ,2))
            DO JST=2, JJ_A_DEDOUB-2
               PSNOWDZN(JJ,JST)=PSNOWDZ(JJ,JST+1)
            ENDDO
            PSNOWDZN(JJ,JJ_A_DEDOUB-1)=0.5*PSNOWDZ(JJ,JJ_A_DEDOUB)
            PSNOWDZN(JJ,JJ_A_DEDOUB)=0.5*PSNOWDZ(JJ,JJ_A_DEDOUB)
       ENDIF
     ENDIF
  ENDIF   
 ENDDO ! end grid points loop
!
! check if bottom layer depth is too small
! in such a case agregation with upper layer and
! looks for an other layer to be splitted
 DO JJ=1,SIZE(PSNOW(:)) ! loop grid points
  IF(.NOT.OMODIF_GRID(JJ)) THEN 
   IF(PSNOWDZ(JJ,INLVLS)< DZMIN_TOP) THEN
      OMODIF_GRID(JJ)=.TRUE.
      ZPENALTY=PSNOWDZ(JJ,INLVLS-2)/ZDZOPT(JJ,INLVLS-2)
      JJ_A_DEDOUB=INLVLS-2
      DO JST= MAX(1,INLVLS-3),1,-1
       ZCRITSIZE=PSNOWDZ(JJ,JST)/ZDZOPT(JJ,JST)
       IF(JST==1.AND.PSNOWDZ(JJ,JST)<DZMIN_BOT) ZCRITSIZE=0.
       IF (ZCRITSIZE > ZPENALTY) THEN
         ZPENALTY=ZCRITSIZE
         JJ_A_DEDOUB=JST
       ENDIF
      ENDDO
      IF(JJ_A_DEDOUB == INLVLS-1) THEN
       PSNOWDZN(JJ,INLVLS)=0.5*(PSNOWDZ(JJ,INLVLS)+PSNOWDZ(JJ,INLVLS-1))
       PSNOWDZN(JJ,INLVLS-1)=0.5*(PSNOWDZ(JJ,INLVLS)+PSNOWDZ(JJ,INLVLS-1))
! write( *,*) '2 couches basales agregees pour creer 2 couches identiques'       
      ELSE
! write( *,*) 'couches basales agregees et on dedouble ',JJ_A_DEDOUB       
       PSNOWDZN(JJ,INLVLS)=(PSNOWDZ(JJ,INLVLS)+PSNOWDZ(JJ,INLVLS-1))
       DO JST=INLVLS-1, JJ_A_DEDOUB+1,-1
          PSNOWDZN(JJ,JST)=PSNOWDZ(JJ,JST-1)
       ENDDO
       PSNOWDZN(JJ,JJ_A_DEDOUB)=0.5*PSNOWDZ(JJ,JJ_A_DEDOUB)
       PSNOWDZN(JJ,JJ_A_DEDOUB+1)=0.5*PSNOWDZ(JJ,JJ_A_DEDOUB)
     ENDIF
   ENDIF
  ENDIF   
 ENDDO ! end grid points loop
ENDIF  ! end specific case INLSVSMIN = INLVLSMAX
!
! case without new snowfall and INVLS > INLVLSMIN
!
! check if surface layer depth is too small
! in such a case agregation with layer beneath
! in case of reaching INLVLSMIN, looks for an other layer to be splitted
DO JJ=1,SIZE(PSNOW(:))
IF((.NOT.GSNOWFALL(JJ)).AND.PSNOW(JJ)> XSNOWCRITD.AND. &
  (.NOT.OMODIF_GRID(JJ))) THEN  
  IF(PSNOWDZ(JJ,1)< DZMIN_TOP_BIS) THEN ! case shallow surface layer 
    OMODIF_GRID(JJ)=.TRUE.
    OAGREG_SURF(JJ)=.TRUE.
    IF(INLVLS_USE(JJ)>INLVLSMIN) THEN ! case minimum not reached 
     INLVLS_USE(JJ)=INLVLS_USE(JJ)-1
     PSNOWDZN(JJ,1)=PSNOWDZ(JJ,1)+PSNOWDZ(JJ,2)
     DO JST=2,INLVLS_USE(JJ) 
          PSNOWDZN(JJ,JST)=PSNOWDZ(JJ,JST+1)
     ENDDO
    ELSE ! case minimum reached         
     ZPENALTY=PSNOWDZ(JJ,2)/ZDZOPT(JJ,2)
     IF(PSNOWDZ(JJ,2)<DZMIN_TOP) ZPENALTY=0.
     JJ_A_DEDOUB=2
     DO JST=3, INLVLS_USE(JJ)
       ZCRITSIZE=PSNOWDZ(JJ,JST)/ZDZOPT(JJ,JST)
       IF(JST==INLVLS_USE(JJ).AND.PSNOWDZ(JJ,JST)<DZMIN_BOT) ZCRITSIZE=0.
       IF (ZCRITSIZE > ZPENALTY) THEN
         ZPENALTY=ZCRITSIZE
         JJ_A_DEDOUB=JST
       ENDIF
     ENDDO
     IF(JJ_A_DEDOUB == 2) THEN ! case splitted layer == 2
       PSNOWDZN(JJ,1)=0.5*(PSNOWDZ(JJ,1)+PSNOWDZ(JJ,2))
       PSNOWDZN(JJ,2)=0.5*(PSNOWDZ(JJ,1)+PSNOWDZ(JJ,2))
     ELSE ! case splitted layer =/ 2
       PSNOWDZN(JJ,1)=PSNOWDZ(JJ,1)+PSNOWDZ(JJ,2)
       DO JST=2, JJ_A_DEDOUB-2
          PSNOWDZN(JJ,JST)=PSNOWDZ(JJ,JST+1)
       ENDDO
       PSNOWDZN(JJ,JJ_A_DEDOUB-1)=0.5*PSNOWDZ(JJ,JJ_A_DEDOUB)
       PSNOWDZN(JJ,JJ_A_DEDOUB)=0.5*PSNOWDZ(JJ,JJ_A_DEDOUB)
     ENDIF ! end case splitted layer =/ 2
    ENDIF ! end case minimum reached
   ENDIF ! end case shallow surface layer
ENDIF   
ENDDO ! end grid points loop
!
! check if bottom layer depth is too small
! in such a case agregation with above layer 
! in case of reaching INLVLSMIN, looks for an other layer to be splitted
DO JJ=1,SIZE(PSNOW(:))
IF((.NOT.GSNOWFALL(JJ)).AND.(.NOT.OAGREG_SURF(JJ)).AND.  &
        (.NOT.OMODIF_GRID(JJ)).AND. PSNOW(JJ)> XSNOWCRITD) THEN  
   IF(PSNOWDZ(JJ,INLVLS_USE(JJ))< DZMIN_TOP) THEN ! case shallow bottom layer
      OMODIF_GRID(JJ)=.TRUE.
     IF(INLVLS_USE(JJ)>INLVLSMIN) THEN ! case minimum not reached
        OMODIF_GRID(JJ)=.TRUE.
        INLVLS_USE(JJ)=INLVLS_USE(JJ)-1
        PSNOWDZN(JJ,INLVLS_USE(JJ))=PSNOWDZ(JJ,INLVLS_USE(JJ))+  &
                  PSNOWDZ(JJ,INLVLS_USE(JJ)+1) 
     ELSE  ! case minimum reached
     ZPENALTY=PSNOWDZ(JJ,INLVLS_USE(JJ)-2)/ZDZOPT(JJ,INLVLS_USE(JJ)-2)
     JJ_A_DEDOUB=INLVLS_USE(JJ)-2
     DO JST= MAX(1,INLVLS_USE(JJ)-3),1,-1
       ZCRITSIZE=PSNOWDZ(JJ,JST)/ZDZOPT(JJ,JST)
       IF(JST==1.AND.PSNOWDZ(JJ,JST)<DZMIN_BOT) ZCRITSIZE=0.
       IF (ZCRITSIZE > ZPENALTY) THEN
         ZPENALTY=ZCRITSIZE
         JJ_A_DEDOUB=JST
       ENDIF
     ENDDO
     IF(JJ_A_DEDOUB == INLVLS_USE(JJ)-1) THEN
       PSNOWDZN(JJ,INLVLS_USE(JJ))=0.5*(PSNOWDZ(JJ,INLVLS_USE(JJ))+   &
                 PSNOWDZ(JJ,INLVLS_USE(JJ)-1)) 
       PSNOWDZN(JJ,INLVLS_USE(JJ)-1)=0.5*(PSNOWDZ(JJ,INLVLS_USE(JJ))+ &
                 PSNOWDZ(JJ,INLVLS_USE(JJ)-1)) 
     ELSE
       PSNOWDZN(JJ,INLVLS_USE(JJ))=(PSNOWDZ(JJ,INLVLS_USE(JJ))+    &
                 PSNOWDZ(JJ,INLVLS_USE(JJ)-1)) 
       DO JST=INLVLS_USE(JJ)-1, JJ_A_DEDOUB+1,-1
          PSNOWDZN(JJ,JST)=PSNOWDZ(JJ,JST-1)
       ENDDO
       PSNOWDZN(JJ,JJ_A_DEDOUB)=0.5*PSNOWDZ(JJ,JJ_A_DEDOUB)
       PSNOWDZN(JJ,JJ_A_DEDOUB+1)=0.5*PSNOWDZ(JJ,JJ_A_DEDOUB)
     ENDIF ! end case splitted layer =/ 2
    ENDIF ! end case minimum reached
   ENDIF ! end case shallow surface layer
ENDIF   
ENDDO ! end grid points loop       
!
! case whithout new snow fall and without a previous grid resize 
! looks for a shallow layer to be splitted according to its depth and to
! the optimal grid size
DO JJ=1,SIZE(PSNOW(:))
IF (.NOT.OMODIF_GRID(JJ))THEN
        DO JST=1,INLVLS-3
          IF (JST<=INLVLS_USE(JJ).AND.INLVLS_USE(JJ)<INLVLS-3.AND. &
            (.NOT.OMODIF_GRID(JJ))) THEN 
                IF(PSNOWDZ(JJ,JST)>(SPLIT_COEF-FLOAT(INLVLS-INLVLS_USE(JJ))/&
                  MAX(1.,FLOAT(INLVLS-INLVLSMIN)))*ZDZOPT(JJ,JST)) THEN 
                    DO JST_1=INLVLS_USE(JJ)+1, JST+2,-1
                        PSNOWDZN(JJ,JST_1)=PSNOWDZ(JJ,JST_1-1)
                        ZDZOPT(JJ,JST_1)=ZDZOPT(JJ,JST_1-1)
                    ENDDO
                    PSNOWDZN(JJ,JST+1)=0.5*PSNOWDZ(JJ,JST)
                    PSNOWDZN(JJ,JST)=0.5*PSNOWDZ(JJ,JST)
                    INLVLS_USE(JJ)=INLVLS_USE(JJ)+1
                    OMODIF_GRID(JJ)=.TRUE.
                ENDIF

          ENDIF
        ENDDO
ENDIF   
ENDDO
!
! case whithout new snow fall and without a previous grid resize 
! looks for a deep layer to be agregated to the layer beneath if similar
! according to its depth and to the optimal grid size
DO JJ=1,SIZE(PSNOW(:))
IF (.NOT.OMODIF_GRID(JJ))THEN
  DO JST=2,INLVLS
    IF (JST<=INLVLS_USE(JJ)-1.AND.INLVLS_USE(JJ)>INLVLSMIN+2.AND. &
            (.NOT.OMODIF_GRID(JJ))) THEN 
       ZDIFTYPE_INF= SNOW3LDIFTYP(PSNOWGRAN1(JJ,JST+1),PSNOWGRAN1(JJ, JST),&
                        PSNOWGRAN2(JJ,JST+1),PSNOWGRAN2(JJ, JST)) 
       ZDIFTYPE_INF=MAX(DIFF_1,MIN(DIFF_MAX,ZDIFTYPE_INF))
       IF(PSNOWDZ(JJ,JST)<ZDZOPT(JJ,JST)*AGREG_COEF_1/ZDIFTYPE_INF &
             .AND.PSNOWDZ(JJ,JST)+PSNOWDZ(JJ,JST+1)<AGREG_COEF_2*&
              MAX(ZDZOPT(JJ,JST),ZDZOPT(JJ,JST+1))) THEN 
           PSNOWDZN(JJ,JST)= PSNOWDZ(JJ,JST)+PSNOWDZ(JJ,JST+1)
           ZDZOPT(JJ,JST)=ZDZOPT(JJ,JST+1)
           DO JST_1=JST+1,INLVLS_USE(JJ)-1
              PSNOWDZN(JJ,JST_1)=PSNOWDZ(JJ,JST_1+1)
              ZDZOPT(JJ,JST_1)=ZDZOPT(JJ,JST_1+1)
           ENDDO
           INLVLS_USE(JJ)=INLVLS_USE(JJ)-1
           OMODIF_GRID(JJ)=.TRUE.
       ENDIF
    ENDIF
  ENDDO
ENDIF   
ENDDO
!
!final check of the consistensy of the new grid size
!
DO JJ=1,SIZE(PSNOW(:))
  IF(ABS(SUM(PSNOWDZN(JJ,1:INLVLS_USE(JJ)))-PSNOW(JJ)) > UEPSI) THEN
   write(*,*) 'error in grid resizing',JJ, SUM(PSNOWDZN(JJ,1:INLVLS_USE(JJ))), &
               PSNOW(JJ),INLVLS_USE(JJ) 
write( *,*) 'PSNOWDZ:',     PSNOWDZ
write( *,*) 'PSNOWDZN:',     PSNOWDZN
CALL ABOR1_SFX("SNOWCRO: error in grid resizing")
   ENDIF       
ENDDO
IF (LHOOK) CALL DR_HOOK('SNOWNLFALL_UPGRID',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWNLFALL_UPGRID

!###############################################################################
!###############################################################################
!################################################################################
!################################################################################
!
SUBROUTINE SNOWNLGRIDFRESH_1D (PSNOW,PSNOWDZ,PSNOWDZN, &
                                PSNOWRHO,PSNOWHEAT,PSNOWGRAN1,PSNOWGRAN2,   &
                                  PSNOWHIST,PSNOWAGE,GSNOWFALL, &
                                  PSNOWRHOF, PSNOWDZF,PSNOWHEATF,PSNOWGRAN1F,&
                                      PSNOWGRAN2F, PSNOWHISTF,PSNOWAGEF,     &
                                      INLVLS_USE) 
!
!!    PURPOSE
!!    -------
!     Snow mass,heat and characteristics redistibution in case of
!     grid resizing. Total mass and heat content of the overall snowpack
!     unchanged/conserved within this routine.
!     Grain size and type of mixed layers is deduced from the conservation
!     of the average optical size
!
!!    AUTHOR
!!    ------
!!      E. Brun           * Meteo-France *
!!
!
USE MODD_SNOW_PAR, ONLY : XSNOWCRITD
USE MODD_SNOW_METAMO
USE MODE_SNOW3L
!
!
IMPLICIT NONE
!
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                  :: PSNOW
!
REAL, DIMENSION(:), INTENT(INOUT) :: PSNOWHEAT, PSNOWRHO, PSNOWDZ, &
                                        PSNOWDZN, PSNOWGRAN1, PSNOWGRAN2, &
                                        PSNOWHIST 
REAL, DIMENSION(:), INTENT(INOUT) :: PSNOWAGE  
REAL,  INTENT(IN)                 :: PSNOWRHOF, PSNOWDZF,PSNOWHEATF,&
                                        PSNOWGRAN1F,PSNOWGRAN2F, PSNOWHISTF 
REAL, INTENT(IN)                  :: PSNOWAGEF                                 
!
INTEGER, INTENT(IN)               :: INLVLS_USE
!
LOGICAL, INTENT(IN)               :: GSNOWFALL
!
!*      0.2    declarations of local variables
!
INTEGER JST 
!
INTEGER                             :: INLVLS_OLD, INLVLS_NEW
!
REAL, DIMENSION(SIZE(PSNOWRHO,1)) :: ZSNOWRHON,ZSNOWGRAN1N,ZSNOWGRAN2N,     &
                                        ZSNOWHEATN,ZSNOWHISTN &
           , ZSNOWZTOP_NEW, ZSNOWZBOT_NEW                                       
REAL,DIMENSION(SIZE(PSNOWRHO,1)) ::ZSNOWAGEN          
REAL, DIMENSION(SIZE(PSNOWRHO,1)+1) :: ZSNOWRHOO,ZSNOWGRAN1O,ZSNOWGRAN2O,     &
                                        ZSNOWHEATO,ZSNOWHISTO,ZSNOWDZO &
           , ZSNOWZTOP_OLD,  ZSNOWZBOT_OLD                                        
REAL,DIMENSION(SIZE(PSNOWRHO,1)+1) ::ZSNOWAGEO          
!
                                       
REAL, PARAMETER :: D1=1., D2=3., D3=4., X=99., VALB5= .95, VALB6=15.4                                      
INTEGER :: JST_OLD,JST_NEW
REAL :: ZDENTMOYN ,ZSPHERMOYN, ZALBMOYN, ZMASTOTN,ZMASTOTO, ZSNOWHEAN, &
          ZSNOWHEAO,ZPROPOR,ZMASDZ_OLD, ZDIAM,ZHISTMOYN 
REAL :: ZAGEMOYN
REAL :: ZPSNOW_OLD, ZPSNOW_NEW
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------------
!
! 0. Initialization:
! ------------------
!
! starts by checking the consistency between both vertical grid sizes
!
IF (LHOOK) CALL DR_HOOK('SNOWNLGRIDFRESH_1D',0,ZHOOK_HANDLE)
INLVLS_NEW     = INLVLS_USE 
INLVLS_OLD = -1
ZPSNOW_NEW=0.
ZPSNOW_OLD=0.
!
DO JST_NEW=1,INLVLS_NEW
ZPSNOW_NEW=ZPSNOW_NEW+PSNOWDZN(JST_NEW)
ENDDO
IF (ABS(ZPSNOW_NEW-PSNOWDZF)<UEPSI) THEN
        INLVLS_OLD=0
ELSE
   DO JST_OLD=1,SIZE(PSNOWRHO)
     IF(PSNOWDZ(JST_OLD)>=UEPSI) THEN
        ZPSNOW_OLD=ZPSNOW_OLD+PSNOWDZ(JST_OLD)
        IF (ABS(ZPSNOW_NEW-PSNOWDZF-ZPSNOW_OLD)<UEPSI) THEN
           INLVLS_OLD=JST_OLD
     ENDIF
     ENDIF
   ENDDO
   IF (INLVLS_OLD==-1) THEN
       write(*,*)'pb INLVLS_OLD INLVLS_NEW=',INLVLS_NEW
       write(*,*)'pb INLVLS_OLD',PSNOWDZF
       write(*,*)'pb INLVLS_OLD',PSNOWDZN
       write(*,*)'pb INLVLS_OLD',PSNOWDZ
       CALL ABOR1_SFX('SNOWCRO: Error INLVLS_OLD')
   ENDIF
ENDIF
!
IF (GSNOWFALL) INLVLS_OLD=INLVLS_OLD+1
ZPSNOW_OLD=PSNOW
ZPSNOW_NEW=ZPSNOW_OLD

! initialization of variables describing the initial snowpack + new snowfall

IF(GSNOWFALL) THEN
        DO JST_OLD=2, INLVLS_OLD
            ZSNOWDZO(JST_OLD) =PSNOWDZ(JST_OLD-1)
            ZSNOWRHOO(JST_OLD) =PSNOWRHO(JST_OLD-1)
            ZSNOWHEATO(JST_OLD) = PSNOWHEAT(JST_OLD-1)
            ZSNOWGRAN1O(JST_OLD) = PSNOWGRAN1(JST_OLD-1)
            ZSNOWGRAN2O(JST_OLD) = PSNOWGRAN2(JST_OLD-1)
            ZSNOWHISTO(JST_OLD) = PSNOWHIST(JST_OLD-1)
            ZSNOWAGEO(JST_OLD) = PSNOWAGE(JST_OLD-1)
        ENDDO
            ZSNOWDZO(1) =PSNOWDZF
            ZSNOWRHOO(1) =PSNOWRHOF
            ZSNOWHEATO(1) = PSNOWHEATF
            ZSNOWGRAN1O(1) = PSNOWGRAN1F
            ZSNOWGRAN2O(1) = PSNOWGRAN2F
            ZSNOWHISTO(1) = PSNOWHISTF
            ZSNOWAGEO(1) = PSNOWAGEF
ELSE
        DO JST_OLD=1, INLVLS_OLD
            ZSNOWDZO(JST_OLD) =PSNOWDZ(JST_OLD)
            ZSNOWRHOO(JST_OLD) =PSNOWRHO(JST_OLD)
            ZSNOWHEATO(JST_OLD) = PSNOWHEAT(JST_OLD)
            ZSNOWGRAN1O(JST_OLD) = PSNOWGRAN1(JST_OLD)
            ZSNOWGRAN2O(JST_OLD) = PSNOWGRAN2(JST_OLD)
            ZSNOWHISTO(JST_OLD) = PSNOWHIST(JST_OLD)
            ZSNOWAGEO(JST_OLD) = PSNOWAGE(JST_OLD)
        ENDDO
ENDIF        

            

!
! 1. Calculate vertical grid limits (m):
! --------------------------------------
!
ZSNOWZTOP_OLD(1) = ZPSNOW_OLD
ZSNOWZTOP_NEW(1) = ZPSNOW_NEW
ZSNOWZBOT_OLD(1) = ZSNOWZTOP_OLD(1)-ZSNOWDZO(1)
ZSNOWZBOT_NEW(1) = ZSNOWZTOP_NEW(1)-PSNOWDZN(1)
!
DO JST_OLD=2,INLVLS_OLD
   ZSNOWZTOP_OLD(JST_OLD) = ZSNOWZBOT_OLD(JST_OLD-1)
   ZSNOWZBOT_OLD(JST_OLD) = ZSNOWZTOP_OLD(JST_OLD)-ZSNOWDZO(JST_OLD)
ENDDO
! Check consistency
IF (ABS(ZSNOWZBOT_OLD(INLVLS_OLD)) > 0.00001) write (*,*) 'Error bottom OLD'
!
ZSNOWZBOT_OLD(INLVLS_OLD)=0.
DO JST_NEW=2,INLVLS_NEW
   ZSNOWZTOP_NEW(JST_NEW) = ZSNOWZBOT_NEW(JST_NEW-1)
   ZSNOWZBOT_NEW(JST_NEW) = ZSNOWZTOP_NEW(JST_NEW)-PSNOWDZN(JST_NEW)
ENDDO
! Check consistency
if (ABS(ZSNOWZBOT_NEW(INLVLS_NEW)) > 0.00001) write (*,*) 'Error bottom NEW'
ZSNOWZBOT_NEW(INLVLS_NEW)=0.
!
! 3. Calculate mass, heat, charcateristics mixing due to vertical grid resizing:
! --------------------------------------------------------------------
!

! loop over the new snow layers
! Summ or avergage of the constituting quantities of the old snow layers
! which are totally or partially inserted in the new snow layer

DO JST_NEW=1,INLVLS_NEW
    ZSNOWHEAN=0.
    ZMASTOTN=0.
    ZDENTMOYN=0.
    ZSPHERMOYN=0.
    ZALBMOYN=0.
    ZDIAM=0.
    ZHISTMOYN=0.
    ZAGEMOYN=0.

    ! lopp over the ols snow layers 
    DO JST_OLD=1, INLVLS_OLD
        IF( ZSNOWZTOP_OLD(JST_OLD) <= ZSNOWZBOT_NEW(JST_NEW)) THEN
                ! JST_OLD lower than JJ_NEW ==> no contribution
        ELSEIF ( ZSNOWZBOT_OLD(JST_OLD) >= ZSNOWZTOP_NEW(JST_NEW)) THEN
                ! JST_OLD higher than JJ_NEW ==> no contribution
        ELSE
                ! old layer contributes to the new one
          ! computation of its contributing ratio and mass/heat 
          ZPROPOR= (MIN(ZSNOWZTOP_OLD(JST_OLD), ZSNOWZTOP_NEW(JST_NEW))&
                 - MAX(ZSNOWZBOT_OLD(JST_OLD), ZSNOWZBOT_NEW(JST_NEW))) &
                 / ZSNOWDZO(JST_OLD) 
          ZMASDZ_OLD= ZPROPOR*ZSNOWRHOO(JST_OLD)*ZSNOWDZO(JST_OLD)   
          ZMASTOTN=ZMASTOTN + ZMASDZ_OLD
          ZSNOWHEAN=ZSNOWHEAN+ZPROPOR*ZSNOWHEATO(JST_OLD)
          ! contribution to the grain optical size and then to the albedo
          IF(ZSNOWGRAN1O(JST_OLD)<0.) THEN
            ZDIAM=-ZSNOWGRAN1O(JST_OLD)*D1/X+(1.+ZSNOWGRAN1O(JST_OLD)/X)* &
                  (ZSNOWGRAN2O(JST_OLD)*D2/X+(1.-ZSNOWGRAN2O(JST_OLD)/X)*D3) 
            ZDIAM=ZDIAM/10000.      
            ZDENTMOYN= ZDENTMOYN-ZMASDZ_OLD*ZSNOWGRAN1O(JST_OLD)/X
            ZSPHERMOYN=ZSPHERMOYN+ZMASDZ_OLD*ZSNOWGRAN2O(JST_OLD)/X
          ELSE
            ZDIAM=ZSNOWGRAN2O(JST_OLD)
            ZDENTMOYN= ZDENTMOYN+ZMASDZ_OLD*0.
            ZSPHERMOYN=ZSPHERMOYN+ZMASDZ_OLD*ZSNOWGRAN1O(JST_OLD)/X
          ENDIF
          ZALBMOYN=ZALBMOYN+MAX(0.,ZMASDZ_OLD*(VALB5-VALB6*SQRT(ZDIAM)))
          ZHISTMOYN=ZHISTMOYN+ZMASDZ_OLD*ZSNOWHISTO(JST_OLD)
          ZAGEMOYN=ZAGEMOYN+ZMASDZ_OLD*ZSNOWAGEO(JST_OLD)
        ENDIF
    ENDDO 
! the new layer inherits from the weihted average properties of the old ones
     ! heat and mass
     ZSNOWHEATN(JST_NEW)= ZSNOWHEAN
     ZSNOWRHON(JST_NEW)= ZMASTOTN/PSNOWDZN(JST_NEW)
     ! grain type and size decuced from the average albedo
     ZALBMOYN=ZALBMOYN/ZMASTOTN
     ZSPHERMOYN=MAX(0.,ZSPHERMOYN/ZMASTOTN)
     ZDENTMOYN=MAX(0.,ZDENTMOYN/ZMASTOTN)
     ZDIAM=((VALB5-ZALBMOYN)/VALB6)**2
     IF (ZDIAM <D2/10000. - 0.0000001) THEN
         ! dendricity is preserved if possible and sphericity is adjusted
         ZSNOWGRAN1N(JST_NEW)=-X*ZDENTMOYN
               IF(ABS(ZSNOWGRAN1N(JST_NEW)+X)< 0.01) THEN
                 ZSNOWGRAN2N(JST_NEW)=X*ZSPHERMOYN                  
               ELSEIF  (ABS(ZSNOWGRAN1N(JST_NEW))< 0.0001) THEN
               ! dendritic snow               
                ZSNOWGRAN1N(JST_NEW)=X*ZSPHERMOYN
                ZSNOWGRAN2N(JST_NEW)=ZDIAM 
               ELSE ! non dendritic
                ZSNOWGRAN2N(JST_NEW)=X*((ZDIAM*10000.+ZSNOWGRAN1N(JST_NEW)*D1/X) &
                  / (1.+ZSNOWGRAN1N(JST_NEW)/X)-D3)/(D2-D3) 
                 IF (ZSNOWGRAN2N(JST_NEW)<0.) THEN
                   ZSNOWGRAN2N(JST_NEW)=0.
                 ENDIF  
                 IF (ZSNOWGRAN2N(JST_NEW)>  X + 0.0000001) THEN
                   ZSNOWGRAN2N(JST_NEW)=X
                 ENDIF  
               ENDIF
     ELSEIF  (ZDIAM >D3/10000.) THEN
                ZSNOWGRAN1N(JST_NEW)=X*ZSPHERMOYN
                ZSNOWGRAN2N(JST_NEW)=ZDIAM 
     ELSEIF(ZDENTMOYN<= 0.+0.0000001) THEN
               ! size between D2 and D3 and dendricity == 0         
               ZSNOWGRAN1N(JST_NEW)=X*ZSPHERMOYN
               ZSNOWGRAN2N(JST_NEW)=ZDIAM
     ELSE
               ! size between D2 and D3 and dendricity < 0         
               ! sphericity is firts preserved, if possible. If not,
               ! denditricity =0
               ZSNOWGRAN1N(JST_NEW)=-X*ZDENTMOYN
               ZSNOWGRAN2N(JST_NEW)=X*((ZDIAM*10000.+ZSNOWGRAN1N(JST_NEW)*D1/X) &
                / (1.+ZSNOWGRAN1N(JST_NEW)/X)-D3)/(D2-D3) 
               IF ( ZSNOWGRAN2N(JST_NEW)<0..OR. ZSNOWGRAN2N(JST_NEW)> X) THEN
               ! inconsistency with ZDIAM ==>  dendricity = 0  
                  ZSNOWGRAN1N(JST_NEW)=X*ZSPHERMOYN
                  ZSNOWGRAN2N(JST_NEW)=ZDIAM
               ENDIF
     ENDIF
        ZSNOWHISTN(JST_NEW)=NINT(ZHISTMOYN/ZMASTOTN)
        ZSNOWAGEN(JST_NEW)=ZAGEMOYN/ZMASTOTN
ENDDO        
!
! check of consistency between new and old snowpacks
    ZSNOWHEAN=0.
    ZMASTOTN=0.
    ZSNOWHEAO=0.
    ZMASTOTO=0.    
    ZPSNOW_NEW=0.
    ZPSNOW_OLD=0.

DO JST=1,INLVLS_NEW
ZSNOWHEAN=ZSNOWHEAN+ZSNOWHEATN(JST)
ZMASTOTN=ZMASTOTN+ZSNOWRHON(JST)*PSNOWDZN(JST)
ZPSNOW_NEW=ZPSNOW_NEW +PSNOWDZN(JST)
ENDDO
DO JST=1,INLVLS_OLD
ZSNOWHEAO=ZSNOWHEAO+ZSNOWHEATO(JST)
ZMASTOTO=ZMASTOTO+ZSNOWRHOO(JST)*ZSNOWDZO(JST)
ZPSNOW_OLD=ZPSNOW_OLD +ZSNOWDZO(JST)
ENDDO
    if (abs(ZSNOWHEAN-ZSNOWHEAO)>0.0001.OR.ABS(ZMASTOTN-ZMASTOTO)>0.0001 &
     .OR. ABS(ZPSNOW_NEW-ZPSNOW_OLD)> 0.0001) THEN 
write(*,*) 'Warning diff', ZSNOWHEAN-ZSNOWHEAO,ZMASTOTN-ZMASTOTO,ZPSNOW_NEW-ZPSNOW_OLD
ENDIF


!
! 5. Update mass (density and thickness) and heat:
! ------------------------------------------------
!
PSNOWDZ(:)    = PSNOWDZN(:)
PSNOWRHO(:)   = ZSNOWRHON(:)
PSNOWHEAT(:)  = ZSNOWHEATN(:)
PSNOWGRAN1(:) = ZSNOWGRAN1N(:)
PSNOWGRAN2(:) = ZSNOWGRAN2N(:)
PSNOWHIST(:)  =  ZSNOWHISTN(:)
PSNOWAGE(:)  =  ZSNOWAGEN(:)
        IF (LHOOK) CALL DR_HOOK('SNOWNLGRIDFRESH_1D',1,ZHOOK_HANDLE)
!
!
!--------------------------------------------------------------------------------
!
        END SUBROUTINE SNOWNLGRIDFRESH_1D
!####################################################################
!###################################################################
        SUBROUTINE SNOWDRIFT(PTSTEP, PVMOD, PSNOWRHO,PSNOWDZ,PSNOW,    &
                    PSNOWGRAN1, PSNOWGRAN2,PSNOWHIST,INLVLS_USE) 
!
!!    PURPOSE
!!    -------
!     Snow compaction  and metamorphism due to drift
!     Mass is unchanged: layer thickness is reduced
!     in proportion to density increases. Method inspired from
!     Brun et al. (1997) and Guyomarch
!     
!     - computes a mobility index of each snow layer from its grains, density
!                 and history
!     - computes a drift index of each layer from its mobility and wind speed
!     - computes a transport index with an exponential decay taking into
!                 account its depth and the mobility of upper layers
!     - increases density and changes grains in case of transport
!
!     HISTORY:
!     Basic parameterization from Crocus/ARPEGE Coupling (1997) 
!     Implementation in V5 
!     Insertion in V6 of grains type evolution in case of dendritic snow (V.
!     Vionnet) 
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP
!
REAL, DIMENSION(:), INTENT(IN)      :: PVMOD
!

REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWRHO, PSNOWDZ,PSNOWGRAN1,    &
                                        PSNOWGRAN2,PSNOWHIST 
REAL, DIMENSION(:), INTENT(OUT)     :: PSNOW
!
INTEGER, DIMENSION(:), INTENT(IN)   :: INLVLS_USE
!
!*      0.2    declarations of local variables
!
INTEGER                             :: JJ,JST   ! looping indexes
!
!
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWRHO2
!
REAL                                :: ZPROFEQU, ZRMOB, ZRDRIFT, ZRT, ZDRO, &
                                         ZDGR1, ZDGR2 
!
! Calibration coefficients
REAL, PARAMETER                      ::VTIME= 48*3600. ! characteristic time for
!compaction and metamorphism by wind drift
!
REAL, PARAMETER                      ::VROMAX= 350. !  maximum density for
! drift compaction     UNIT : kg m-3
REAL, PARAMETER                      ::VROMIN= 50.  !  minimum density for
! mobility computation UNIT : kg m-3
REAL, PARAMETER                      ::VMOB1= 0.295  !  coefficient for computing
! the mobility index
REAL, PARAMETER                      ::VMOB2= 0.833  !  coefficient for computing
! the mobility index
REAL, PARAMETER                      ::VMOB3= 0.583  !  coefficient for computing
! the mobility index
REAL, PARAMETER                      ::VMOB4= -0.0583 !  coefficient for computing
! the mobility index
REAL, PARAMETER                      ::VDRIFT1= 2.868 !  coefficient for computing
! the drift index
REAL, PARAMETER                      ::VDRIFT2= 0.085 !  coefficient for computing
! the drift index
REAL, PARAMETER                      ::VDRIFT3= 3.25  !  coefficient for computing
! the drift index
REAL, PARAMETER                      ::VSIZEMIN=3.E-4 !  minimum size decrease 
REAL(KIND=JPRB) :: ZHOOK_HANDLE
! by drift  UNIT = m
!
!-------------------------------------------------------------------------------
!
! 0. Initialization:
! ------------------
!
IF (LHOOK) CALL DR_HOOK('SNOWDRIFT',0,ZHOOK_HANDLE)
DO JJ=1, SIZE(PSNOW)
  DO JST=1, INLVLS_USE(JJ)
    ZSNOWRHO2(JJ,JST)  = PSNOWRHO(JJ,JST)
  ENDDO
ENDDO
!
! 1. Computation of drift and induced settling and metamorphism
! ------------------
!
DO JJ=1, SIZE(PSNOW)
!
!   initialization decay coeff
      ZPROFEQU=0.
  DO JST=1, INLVLS_USE(JJ)
!  mobility index computation of a layer as a function of its properties
      IF(PSNOWGRAN1(JJ,JST) < 0.) THEN        
!     dendritic case              
          ZRMOB=0.34*(0.5+0.75*(-PSNOWGRAN1(JJ,JST)/99.)-0.5*  &
           PSNOWGRAN2(JJ,JST)/99.)+   &
           0.66*(1.25-1.25*(MAX(PSNOWRHO(JJ,JST),VROMIN)/1000.-VROMIN/1000.) / VMOB1)  
      ELSE
!     non dendritic case         
          ZRMOB=0.34*(VMOB2-VMOB2*PSNOWGRAN1(JJ,JST)/99.-VMOB3*  &
           PSNOWGRAN2(JJ,JST)*10000./100.)+    &
           0.66*(1.25-1.25*(MAX(PSNOWRHO(JJ,JST),VROMIN)/1000.-VROMIN/1000.) / VMOB1) 
      ENDIF
!     correction in case of former wet snow
      IF (PSNOWHIST(JJ,JST) >= 2) ZRMOB = MIN(ZRMOB, VMOB4)
!      
!     computation of drift index supposing no overburden snow
      ZRDRIFT = ZRMOB - (VDRIFT1 *exp(-VDRIFT2*PVMOD(JJ))-1.)
!     update the decay coeff by half the current layer
      ZPROFEQU = ZPROFEQU  +0.5 * PSNOWDZ(JJ,JST) *0.1 * (VDRIFT3-ZRDRIFT)
!     computation of the drift index inclunding the decay by overburden snow 
      ZRT = max(0.,ZRDRIFT * exp(-ZPROFEQU*100))
!      
!     settling by wind transport only in case of not too dense snow
      IF(PSNOWRHO(JJ,JST) < VROMAX) THEN
          ZDRO = (VROMAX - PSNOWRHO(JJ,JST))*PTSTEP/ VTIME *ZRT
          PSNOWRHO(JJ,JST) = MIN(VROMAX , PSNOWRHO(JJ,JST)+ZDRO)
          PSNOWDZ(JJ,JST) = PSNOWDZ(JJ,JST)*ZSNOWRHO2(JJ,JST)/   &
                             PSNOWRHO(JJ,JST) 
      ENDIF
!      
!     metamorphism induced by snow drift
      IF(PSNOWGRAN1(JJ,JST) < 0.) THEN        
!     dendritic case              
          ZDGR1 = -PSNOWGRAN1(JJ,JST)*PTSTEP/ VTIME *ZRT*0.5
          ZDGR1 = MIN(ZDGR1,-0.99*PSNOWGRAN1(JJ,JST))
          PSNOWGRAN1(JJ,JST) = PSNOWGRAN1(JJ,JST)+ZDGR1
! modif_VV_140910          
          ZDGR2 = PTSTEP/ VTIME *ZRT*(99.-PSNOWGRAN2(JJ,JST))
          PSNOWGRAN2(JJ,JST) = MIN(99.,PSNOWGRAN2(JJ,JST)+ZDGR2)
! fin modif_VV_140910          
      ELSE
!     non dendritic case         
          ZDGR1 = PTSTEP/ VTIME *ZRT*(99.-PSNOWGRAN1(JJ,JST))
          ZDGR2 = PTSTEP/ VTIME *ZRT*5./10000.
          PSNOWGRAN1(JJ,JST) = MIN(99.,PSNOWGRAN1(JJ,JST)+ZDGR1)
          PSNOWGRAN2(JJ,JST) = MAX(VSIZEMIN,PSNOWGRAN2(JJ,JST)-ZDGR2)
      ENDIF

!     update the decay coeff by half the current layer
      ZPROFEQU = ZPROFEQU  +0.5 * PSNOWDZ(JJ,JST) *0.1 * (VDRIFT3-ZRDRIFT)
!      
  ENDDO  ! snow layers loop
ENDDO    ! grid points loop
!
! 2. Update total snow depth:
! -----------------------------------------------
!
! Compaction of total snowpack depth
!
DO JJ=1, SIZE(PSNOWDZ,1)
   PSNOW(JJ)=SUM(PSNOWDZ(JJ,1:INLVLS_USE(JJ)))
ENDDO
        IF (LHOOK) CALL DR_HOOK('SNOWDRIFT',1,ZHOOK_HANDLE)

!
        END SUBROUTINE SNOWDRIFT
!####################################################################
!###################################################################
       SUBROUTINE TRIDIAG_GROUND_SNOWCRO(PA,PB,PC,PY,PX,INLVLS_USE,IDIFLOOP)
!      #########################################
!
!
!!****   *TRIDIAG_GROUND* - routine to solve a time implicit scheme
!!
!!
!!     PURPOSE
!!     -------
!        The purpose of this routine is to resolve the linear system:
!
!       A.X = Y
!
!      where A is a tridiagonal matrix, and X and Y two vertical vectors.
!     However, the computations are performed at the same time for all
!     the verticals where an inversion of the system is necessary.
!     This explain the dimansion of the input variables.
!
!!**   METHOD
!!     ------
!!                      
!!        Then, the classical tridiagonal algorithm is used to invert the 
!!     implicit operator. Its matrix is given by:
!!
!!     (  b(1)      c(1)      0        0        0         0        0        0  )
!!     (  a(2)      b(2)     c(2)      0  ...    0        0        0        0  ) 
!!     (   0        a(3)     b(3)     c(3)       0        0        0        0  ) 
!!      .......................................................................
!!     (   0   ...   0      a(k)      b(k)     c(k)       0   ...  0        0  ) 
!!      .......................................................................
!!     (   0         0        0        0        0 ...  a(n-1)   b(n-1)   c(n-1))
!!     (   0         0        0        0        0 ...     0      a(n)     b(n) )
!!
!!
!!       All these computations are purely vertical and vectorizations are 
!!     easely achieved by processing all the verticals in parallel.
!!
!!     EXTERNAL
!!     --------
!!
!!       NONE
!!
!!     IMPLICIT ARGUMENTS
!!     ------------------
!!
!!     REFERENCE
!!     ---------
!!
!!     AUTHOR
!!     ------
!!       V. Masson
!! 
!!     MODIFICATIONS
!!     -------------
!!       Original        May 13, 1998
!!       05/2011: Brun  Special treatment to tackle the variable number
!!                      of snow layers 
!!                      In case of second call, a shift of 1 snow layer
!!                      is applied in the control loops.
!! ---------------------------------------------------------------------
!
!*       0. DECLARATIONS
!
IMPLICIT NONE
!
!
!*       0.1 declarations of arguments
!
REAL,    DIMENSION(:,:), INTENT(IN)  :: PA  ! lower diag. elements of A matrix
REAL,    DIMENSION(:,:), INTENT(IN)  :: PB  ! main  diag. elements of A matrix
REAL,    DIMENSION(:,:), INTENT(IN)  :: PC  ! upper diag. elements of A matrix
REAL,    DIMENSION(:,:), INTENT(IN)  :: PY  ! r.h.s. term   
!
REAL,    DIMENSION(:,:), INTENT(OUT) :: PX  ! solution of A.X = Y 
!
INTEGER,    DIMENSION(:), INTENT(IN) :: INLVLS_USE ! number of effective layers 
!
INTEGER           :: IDIFLOOP       ! shift in control loops: 0 or 1
!*       0.2 declarations of local variables
!
INTEGER           :: JK             ! vertical loop control
INTEGER           :: IN             ! number of vertical levels
!
REAL, DIMENSION(SIZE(PA,1)           ) :: ZDET ! work array
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2)) :: ZW   ! work array
REAL(KIND=JPRB) :: ZHOOK_HANDLE
! ---------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TRIDIAG_GROUND_SNOWCRO',0,ZHOOK_HANDLE)
IN=SIZE(PX,2)
!
!
!*       1.  levels going up
!            ---------------
!
!*       1.1 first level
!            -----------
!
ZDET(:)   = PB(:,1)
!
PX  (:,1) = PY(:,1) / ZDET(:)
!
!*       1.2 other levels
!            ------------
!
DO JK=2,IN
WHERE (JK<=INLVLS_USE(:)-IDIFLOOP)
  ZW(:,JK)    = PC(:,JK-1)/ZDET(:)
  ZDET(:)       = PB(:,JK  ) - PA(:,JK)*ZW(:,JK)

  PX  (:,JK)    = ( PY(:,JK) - PA(:,JK)*PX(:,JK-1) ) / ZDET(:)
ENDWHERE
END DO
!
!-------------------------------------------------------------------------------
!
!*       2.  levels going down
!            -----------------
!
DO JK=IN-1,1,-1
WHERE (JK<=INLVLS_USE(:)-1-IDIFLOOP)
  PX  (:,JK) = PX(:,JK) - ZW(:,JK+1)*PX(:,JK+1)
ENDWHERE
END DO
IF (LHOOK) CALL DR_HOOK('TRIDIAG_GROUND_SNOWCRO',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE TRIDIAG_GROUND_SNOWCRO


!####################################################################
!####################################################################
!####################################################################
      SUBROUTINE SNOWCROLAYER_GONE(PTSTEP,PSCAP,PSNOWTEMP,PSNOWDZ,      &
                             PSNOWRHO,PSNOWLIQ,PSNOWGRAN1,PSNOWGRAN2,   &
                             PSNOWHIST,PSNOWAGE,PLES3L, INLVLS_USE) 
!
!
!!    PURPOSE
!     Account for the case when one or several snow layers melt
!     during a time step:
!     in that case, merge these layers with the underlying layer
!     except for the bottom layer which is merged to the abovelying layer
!     energy and mass are conserved
!     a new merged layer keeps the grain, histo and age properties of the 
!     non-melted layer
!
USE MODD_CSTS,ONLY : XTT, XLMTT, XRHOLW, XRHOLI, XLVTT
!
USE MODE_SNOW3L
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, INTENT(IN)                    :: PTSTEP
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSCAP
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWDZ, PSNOWTEMP, PSNOWRHO,   &
                                          PSNOWLIQ 
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWGRAN1,PSNOWGRAN2,PSNOWHIST,PSNOWAGE
!
INTEGER, DIMENSION(:), INTENT(INOUT)  :: INLVLS_USE ! 
!
REAL, DIMENSION(:), INTENT(IN) :: PLES3L
!
!*      0.2    declarations of local variables
!
REAL                                :: ZHEAT, ZMASS, ZDZ, ZLIQ, ZSNOWLWE
!
INTEGER :: JJ,JST,JST_1, JST_2,JST_MAX,IDIFF_LAYER ! loop counter
!
REAL, PARAMETER                     :: ZSNOWDZMIN = 0.0001 ! (m)
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!                                      ZSNOWDZMIN = minimum snow layer thickness
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SNOWCROLAYER_GONE',0,ZHOOK_HANDLE)
DO JJ=1,SIZE(PSNOWRHO,1)  ! loop on gridpoints

     JST_MAX= INLVLS_USE(JJ)
     IDIFF_LAYER=0  ! used as shift counter of previously melted layers
!     
     DO JST_1=JST_MAX,1-JST_MAX,-1 ! loop on 2 x layers in case of multi melt
!   
       JST= JST_1+IDIFF_LAYER     
       ! Merge is possible only in case of 2 active layers or more
       IF(JST>=1.AND.INLVLS_USE(JJ)>1) THEN 
! Total Liquid equivalent water content of snow (m):
       ZSNOWLWE = PSNOWRHO(JJ,JST)*PSNOWDZ(JJ,JST)/XRHOLW
! Consideration of sublimation if any       
       IF(JST==1) ZSNOWLWE = ZSNOWLWE- MAX(0.,PLES3L(JJ)*PTSTEP/(XLSTT*XRHOLW))
!
! Test if avalaible energy exceeds total latent heat
       IF (PSCAP(JJ,JST)*MAX(0.0, PSNOWTEMP(JJ,JST) - XTT)* PSNOWDZ(JJ,JST)>=  &
                   (ZSNOWLWE-PSNOWLIQ(JJ,JST))*XLMTT*XRHOLW-UEPSI) THEN 
        IF(JST==INLVLS_USE(JJ)) THEN
! Case of a total melt of the bottom layer: merge with above layer 
!        which keeps its grain, histo and age properties
           ZHEAT= PSNOWDZ(JJ,JST)*(PSCAP(JJ,JST)*(PSNOWTEMP(JJ,JST)-XTT)  &
                 - XLMTT*PSNOWRHO(JJ,JST) ) + XLMTT*XRHOLW*PSNOWLIQ(JJ,JST) &
                 + PSNOWDZ(JJ,JST-1)*(PSCAP(JJ,JST-1)*(PSNOWTEMP(JJ,JST-1)-XTT)&
                 - XLMTT*PSNOWRHO(JJ,JST-1) ) + XLMTT*XRHOLW*PSNOWLIQ(JJ,JST-1) 
           ZMASS= PSNOWDZ(JJ,JST)*PSNOWRHO(JJ,JST)   &
                 + PSNOWDZ(JJ,JST-1)*PSNOWRHO(JJ,JST-1)    
           ZDZ  = PSNOWDZ(JJ,JST) + PSNOWDZ(JJ,JST-1)   
           ZLIQ = PSNOWLIQ(JJ,JST)+PSNOWLIQ(JJ,JST-1)
           PSNOWDZ(JJ,JST-1)=ZDZ
           PSNOWRHO(JJ,JST-1)=ZMASS/ZDZ
           PSNOWLIQ(JJ,JST-1)=ZLIQ
! Temperature of the merged layer is deduced from the heat content 
           PSCAP(JJ,JST-1)= SNOW3LSCAP(PSNOWRHO(JJ,JST-1)-PSNOWLIQ(JJ,JST-1) &
                           * XRHOLW / MAX(PSNOWDZ(JJ,JST-1),ZSNOWDZMIN)) 
           PSNOWTEMP(JJ,JST-1)=XTT +((((ZHEAT-XLMTT*XRHOLW*PSNOWLIQ(JJ,JST-1)) &
               /PSNOWDZ(JJ,JST-1)) + XLMTT*PSNOWRHO(JJ,JST-1))/PSCAP(JJ,JST-1) ) 
! Decrease the number of active snow layers              
           INLVLS_USE(JJ)=INLVLS_USE(JJ)-1           
!
        ELSE                
! Case of a total melt of the bottom layer: merge with beneath layer 
!        which keeps its grain, histo and age properties
! 
           ZHEAT= PSNOWDZ(JJ,JST)*(PSCAP(JJ,JST)*(PSNOWTEMP(JJ,JST)-XTT)  &
                 - XLMTT*PSNOWRHO(JJ,JST) ) + XLMTT*XRHOLW*PSNOWLIQ(JJ,JST) &
                 + PSNOWDZ(JJ,JST+1)*(PSCAP(JJ,JST+1)*(PSNOWTEMP(JJ,JST+1)-XTT)&
                 - XLMTT*PSNOWRHO(JJ,JST+1) ) + XLMTT*XRHOLW*PSNOWLIQ(JJ,JST+1) 
           ZMASS= PSNOWDZ(JJ,JST)*PSNOWRHO(JJ,JST)   &
                 + PSNOWDZ(JJ,JST+1)*PSNOWRHO(JJ,JST+1)    
           ZDZ  = PSNOWDZ(JJ,JST) + PSNOWDZ(JJ,JST+1)   
           ZLIQ = PSNOWLIQ(JJ,JST)+PSNOWLIQ(JJ,JST+1)
           PSNOWDZ(JJ,JST)=ZDZ
           PSNOWRHO(JJ,JST)=ZMASS/ZDZ
           PSNOWLIQ(JJ,JST)=ZLIQ
           PSNOWGRAN1(JJ,JST)=PSNOWGRAN1(JJ,JST+1)
           PSNOWGRAN2(JJ,JST)=PSNOWGRAN2(JJ,JST+1)
           PSNOWHIST(JJ,JST)=PSNOWHIST(JJ,JST+1)
           PSNOWAGE(JJ,JST)=PSNOWAGE(JJ,JST+1)
! Temperature of the merged layer is deduced from the heat content 
           PSCAP(JJ,JST)= SNOW3LSCAP(PSNOWRHO(JJ,JST)-PSNOWLIQ(JJ,JST) &
                           * XRHOLW / MAX(PSNOWDZ(JJ,JST),ZSNOWDZMIN)) 
           PSNOWTEMP(JJ,JST) = XTT + ( (((ZHEAT-XLMTT*XRHOLW*PSNOWLIQ(JJ,JST)) &
                   /PSNOWDZ(JJ,JST)) + XLMTT*PSNOWRHO(JJ,JST))/PSCAP(JJ,JST) ) 
! Shift the above layers
           DO JST_2=JST+1,INLVLS_USE(JJ)-1
              PSNOWTEMP(JJ,JST_2)=PSNOWTEMP(JJ,JST_2+1)
              PSCAP(JJ,JST_2)=PSCAP(JJ,JST_2+1)
              PSNOWDZ(JJ,JST_2)=PSNOWDZ(JJ,JST_2+1)
              PSNOWRHO(JJ,JST_2)=PSNOWRHO(JJ,JST_2+1)
              PSNOWLIQ(JJ,JST_2)=PSNOWLIQ(JJ,JST_2+1)
              PSNOWGRAN1(JJ,JST_2)=PSNOWGRAN1(JJ,JST_2+1)
              PSNOWGRAN2(JJ,JST_2)=PSNOWGRAN2(JJ,JST_2+1)
              PSNOWHIST(JJ,JST_2)=PSNOWHIST(JJ,JST_2+1)
              PSNOWAGE(JJ,JST_2)=PSNOWAGE(JJ,JST_2+1)
           ENDDO !  loop JST_2
! Decrease the number of active snow layers              
           INLVLS_USE(JJ)=INLVLS_USE(JJ)-1           
! Update the shift counter IDIFF_LAYER
           IDIFF_LAYER=IDIFF_LAYER+1
        ENDIF ! end test of bottom layer       
       ENDIF ! end test on availibility of energy      
      ENDIF ! end test on the number of remaining active layers
     ENDDO ! end loop on the snow layers
ENDDO ! end loop gridpoints
IF (LHOOK) CALL DR_HOOK('SNOWCROLAYER_GONE',1,ZHOOK_HANDLE)
!
END SUBROUTINE SNOWCROLAYER_GONE
!####################################################################
!###################################################################


!
END SUBROUTINE SNOWCRO
