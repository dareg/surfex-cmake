!     ###############################################################################
SUBROUTINE COUPLING_SEAFLUX_n(HPROGRAM, HCOUPLING, PTIMEC,                                   &
                 PTSTEP, KYEAR, KMONTH, KDAY, PTIME, KI, KSV, KSW, PTSUN, PZENITH, PZENITH2, &
                 PAZIM, PZREF, PUREF, PZS, PU, PV, PQA, PTA, PRHOA, PSV, PCO2, HSV,          &
                 PRAIN, PSNOW, PLW, PDIR_SW, PSCA_SW, PSW_BANDS, PPS, PPA,                   &
                 PSFTQ, PSFTH, PSFTS, PSFCO2, PSFU, PSFV,                                    &
                 PTRAD, PDIR_ALB, PSCA_ALB, PEMIS, PTSURF, PZ0, PZ0H, PQSURF,                &
                 PPEW_A_COEF, PPEW_B_COEF,                                                   &
                 PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,                         &
                 HTEST                                                                       )  
!     ###############################################################################
!
!!****  *COUPLING_SEAFLUX_n * - Driver of the WATER_FLUX scheme for sea   
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     V. Masson 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      Modified    01/2006 : sea flux parameterization.
!!      Modified    09/2006 : P. Tulet Introduce Sea salt aerosol Emission/Deposition
!!      Modified    03/2009 : B. Decharme SST could change during a run => ALB and EMIS 
!!      Modified    05/2009 : V. Masson : implicitation of momentum fluxes
!!      Modified    09/2009 : B. Decharme Radiative properties at time t+1
!!      Modified    01/2010 : B. Decharme Add XTTS
!!      Modified    09/2012 : B. Decharme New wind implicitation
!!      Modified    10/2012 : P. Le Moigne CMO1D update
!!      Modified    04/2013 : P. Le Moigne Wind implicitation and SST update displaced
!!      Modified    04/2013 : B. Decharme new coupling variables
!!      Modified    01/2014 : S. Senesi : handle sea-ice cover, sea-ice model interface, 
!!                               and apply to Gelato
!!      Modified    01/2014 : S. Belamari Remove MODE_THERMOS and XLVTT
!!      Modified    05/2014 : S. Belamari New ECUME : Include salinity & atm. pressure impact 
!!                                       
!!---------------------------------------------------------------------
!
USE MODD_REPROD_OPER, ONLY : CIMPLICIT_WIND
!
USE MODD_CSTS,       ONLY : XRD, XCPD, XP00, XTT, XTTS, XTTSI, XDAY
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_SFX_OASIS,  ONLY : LCPL_SEA, LCPL_SEAICE
USE MODD_WATER_PAR,  ONLY : XEMISWAT, XEMISWATICE
!
USE MODD_DATA_SEAFLUX_n,  ONLY : LSST_DATA
USE MODD_SEAFLUX_n,  ONLY : XSST, XSSS, XSIC, XTICE, XZ0, XDIR_ALB, XSCA_ALB, &
                              XEMIS, TTIME, CSEA_ALB, CSEA_FLUX, XUMER, XVMER,&
                              LINTERPOL_SST, LINTERPOL_SSS , LINTERPOL_SIC,   &
                              LINTERPOL_SIT, XFSIC, XFSIT, LPERTFLUX,         &
                              XPERTFLUX,                                      &
                              XICHCE, LPRECIP, LPWEBB , LPWG, XICE_ALB,       &
                              LHANDLE_SIC, CSEAICE_SCHEME, TGLT, NZ0
USE MODD_DIAG_SEAFLUX_n, ONLY : XEVAP
USE MODD_WATER_PAR, ONLY : XALBSEAICE
!
USE MODD_OCEAN_n, ONLY : LMERCATOR                            
USE MODD_CH_SEAFLUX_n, ONLY : CSV, CCH_DRY_DEP, XDEP, NBEQ, NSV_CHSBEG, NSV_CHSEND,&
                                NSV_DSTBEG, NSV_DSTEND, NAEREQ, NDSTEQ, NSLTEQ, &
                                NSV_AERBEG, NSV_AEREND, NSV_SLTBEG, NSV_SLTEND  
!
USE MODI_WATER_FLUX
USE MODI_MR98
USE MODI_ECUME_SEAFLUX
USE MODI_COARE30_SEAFLUX
USE MODI_ADD_FORECAST_TO_DATE_SURF
USE MODI_MOD1D_n
USE MODI_DIAG_INLINE_SEAFLUX_n
USE MODI_CH_AER_DEP
USE MODI_CH_DEP_WATER
USE MODI_DSLT_DEP
USE MODI_SST_UPDATE
USE MODI_INTERPOL_SST_MTH
USE MODI_UPDATE_RAD_SEA
!
USE MODE_DSLT_SURF
USE MODD_DST_SURF
USE MODD_SLT_SURF
USE MODD_DST_n,    ONLY: XEMISRADIUS_DST, XEMISSIG_DST
USE MODD_SLT_n,    ONLY: XEMISRADIUS_SLT, XEMISSIG_SLT
! 
USE MODD_SEAFLUX_GRID_n, ONLY : XLAT
USE MODD_OCEAN_CSTS,   ONLY : NOCKMIN
USE MODD_OCEAN_REL_n,      ONLY : XSEAT_REL
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
USE MODI_COUPLING_ICEFLUX_n
USE MODI_SEAICE_GELATO1D_n
!
USE MODI_COUPLING_SLT_n
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
 CHARACTER(LEN=6),    INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=1),    INTENT(IN)  :: HCOUPLING ! type of coupling
                                              ! 'E' : explicit
                                              ! 'I' : implicit
REAL,                INTENT(IN)  :: PTIMEC    ! current duration since start of the run (s)
INTEGER,             INTENT(IN)  :: KYEAR     ! current year (UTC)
INTEGER,             INTENT(IN)  :: KMONTH    ! current month (UTC)
INTEGER,             INTENT(IN)  :: KDAY      ! current day (UTC)
REAL,                INTENT(IN)  :: PTIME     ! current time since midnight (UTC, s)
INTEGER,             INTENT(IN)  :: KI        ! number of points
INTEGER,             INTENT(IN)  :: KSV       ! number of scalars
INTEGER,             INTENT(IN)  :: KSW       ! number of short-wave spectral bands
REAL, DIMENSION(KI), INTENT(IN)  :: PTSUN     ! solar time                    (s from midnight)
REAL,                INTENT(IN)  :: PTSTEP    ! atmospheric time-step                 (s)
REAL, DIMENSION(KI), INTENT(IN)  :: PZREF     ! height of T,q forcing                 (m)
REAL, DIMENSION(KI), INTENT(IN)  :: PUREF     ! height of wind forcing                (m)
!
REAL, DIMENSION(KI), INTENT(IN)  :: PTA       ! air temperature forcing               (K)
REAL, DIMENSION(KI), INTENT(IN)  :: PQA       ! air humidity forcing                  (kg/m3)
REAL, DIMENSION(KI), INTENT(IN)  :: PRHOA     ! air density                           (kg/m3)
REAL, DIMENSION(KI,KSV),INTENT(IN) :: PSV     ! scalar variables
!                                             ! chemistry:   first char. in HSV: '#'  (molecule/m3)
!                                             !
 CHARACTER(LEN=6), DIMENSION(KSV),INTENT(IN):: HSV  ! name of all scalar variables
REAL, DIMENSION(KI), INTENT(IN)  :: PU        ! zonal wind                            (m/s)
REAL, DIMENSION(KI), INTENT(IN)  :: PV        ! meridian wind                         (m/s)
REAL, DIMENSION(KI,KSW),INTENT(IN) :: PDIR_SW ! direct  solar radiation (on horizontal surf.)
!                                             !                                       (W/m2)
REAL, DIMENSION(KI,KSW),INTENT(IN) :: PSCA_SW ! diffuse solar radiation (on horizontal surf.)
!                                             !                                       (W/m2)
REAL, DIMENSION(KSW),INTENT(IN)  :: PSW_BANDS ! mean wavelength of each shortwave band (m)
REAL, DIMENSION(KI), INTENT(IN)  :: PZENITH   ! zenithal angle at t  (radian from the vertical)
REAL, DIMENSION(KI), INTENT(IN)  :: PZENITH2  ! zenithal angle at t+1(radian from the vertical)
REAL, DIMENSION(KI), INTENT(IN)  :: PAZIM     ! azimuthal angle      (radian from North, clockwise)
REAL, DIMENSION(KI), INTENT(IN)  :: PLW       ! longwave radiation (on horizontal surf.)
!                                             !                                       (W/m2)
REAL, DIMENSION(KI), INTENT(IN)  :: PPS       ! pressure at atmospheric model surface (Pa)
REAL, DIMENSION(KI), INTENT(IN)  :: PPA       ! pressure at forcing level             (Pa)
REAL, DIMENSION(KI), INTENT(IN)  :: PZS       ! atmospheric model orography           (m)
REAL, DIMENSION(KI), INTENT(IN)  :: PCO2      ! CO2 concentration in the air          (kg/m3)
REAL, DIMENSION(KI), INTENT(IN)  :: PSNOW     ! snow precipitation                    (kg/m2/s)
REAL, DIMENSION(KI), INTENT(IN)  :: PRAIN     ! liquid precipitation                  (kg/m2/s)
!
REAL, DIMENSION(KI), INTENT(OUT) :: PSFTH     ! flux of heat                          (W/m2)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFTQ     ! flux of water vapor                   (kg/m2/s)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFU      ! zonal momentum flux                   (Pa)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFV      ! meridian momentum flux                (Pa)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFCO2    ! flux of CO2                           (m/s*kg_CO2/kg_air)
REAL, DIMENSION(KI,KSV),INTENT(OUT):: PSFTS   ! flux of scalar var.                   (kg/m2/s)
!
REAL, DIMENSION(KI), INTENT(OUT) :: PTRAD     ! radiative temperature                 (K)
REAL, DIMENSION(KI,KSW),INTENT(OUT):: PDIR_ALB! direct albedo for each spectral band  (-)
REAL, DIMENSION(KI,KSW),INTENT(OUT):: PSCA_ALB! diffuse albedo for each spectral band (-)
REAL, DIMENSION(KI), INTENT(OUT) :: PEMIS     ! emissivity                            (-)
!
REAL, DIMENSION(KI), INTENT(OUT) :: PTSURF    ! surface effective temperature         (K)
REAL, DIMENSION(KI), INTENT(OUT) :: PZ0       ! roughness length for momentum         (m)
REAL, DIMENSION(KI), INTENT(OUT) :: PZ0H      ! roughness length for heat             (m)
REAL, DIMENSION(KI), INTENT(OUT) :: PQSURF    ! specific humidity at surface          (kg/kg)
!
REAL, DIMENSION(KI), INTENT(IN) :: PPEW_A_COEF! implicit coefficients   (m2s/kg)
REAL, DIMENSION(KI), INTENT(IN) :: PPEW_B_COEF! needed if HCOUPLING='I' (m/s)
REAL, DIMENSION(KI), INTENT(IN) :: PPET_A_COEF
REAL, DIMENSION(KI), INTENT(IN) :: PPEQ_A_COEF
REAL, DIMENSION(KI), INTENT(IN) :: PPET_B_COEF
REAL, DIMENSION(KI), INTENT(IN) :: PPEQ_B_COEF
 CHARACTER(LEN=2),    INTENT(IN) :: HTEST ! must be equal to 'OK'
!
!*      0.2    declarations of local variables
!     
REAL, DIMENSION(KI,KSW) :: ZDIR_ALB   ! Direct albedo at time t
REAL, DIMENSION(KI,KSW) :: ZSCA_ALB   ! Diffuse albedo at time t
!
REAL, DIMENSION(KI) :: ZEXNA      ! Exner function at forcing level
REAL, DIMENSION(KI) :: ZEXNS      ! Exner function at surface level
REAL, DIMENSION(KI) :: ZU         ! zonal wind
REAL, DIMENSION(KI) :: ZV         ! meridian wind
REAL, DIMENSION(KI) :: ZWIND      ! Wind
REAL, DIMENSION(KI) :: ZWORK      ! Work array
REAL, DIMENSION(KI) :: ZCD        ! Drag coefficient on open sea
REAL, DIMENSION(KI) :: ZCD_ICE    ! "     "          on seaice
REAL, DIMENSION(KI) :: ZCDN       ! Neutral Drag coefficient on open sea
REAL, DIMENSION(KI) :: ZCDN_ICE   ! "     "          on seaice
REAL, DIMENSION(KI) :: ZCH        ! Heat transfer coefficient on open sea
REAL, DIMENSION(KI) :: ZCH_ICE    ! "     "          on seaice
REAL, DIMENSION(KI) :: ZCE        ! Vaporization heat transfer coefficient on open sea
REAL, DIMENSION(KI) :: ZCE_ICE    ! "     "          on seaice
REAL, DIMENSION(KI) :: ZRI        ! Richardson number on open sea
REAL, DIMENSION(KI) :: ZRI_ICE    ! "     "          on seaice
REAL, DIMENSION(KI) :: ZRESA_SEA  ! aerodynamical resistance on open sea
REAL, DIMENSION(KI) :: ZRESA_SEA_ICE  ! "     "          on seaice
REAL, DIMENSION(KI) :: ZUSTAR     ! friction velocity (m/s) on open sea
REAL, DIMENSION(KI) :: ZUSTAR_ICE ! "     "          on seaice
REAL, DIMENSION(KI) :: ZZ0        ! roughness length over open sea
REAL, DIMENSION(KI) :: ZZ0_ICE    ! roughness length over seaice
REAL, DIMENSION(KI) :: ZZ0H       ! heat roughness length over open sea
REAL, DIMENSION(KI) :: ZZ0H_ICE   ! heat roughness length over seaice
REAL, DIMENSION(KI) :: ZZ0W       ! Work array for Z0 and Z0H computation
REAL, DIMENSION(KI) :: ZQSAT      ! humidity at saturation on open sea
REAL, DIMENSION(KI) :: ZQSAT_ICE  ! "     "          on seaice
!
REAL, DIMENSION(KI) :: ZSFTH      ! Heat flux for open sea (and for sea-ice points if merged)
REAL, DIMENSION(KI) :: ZSFTQ      ! Water vapor flux on open sea (and for sea-ice points if merged)
REAL, DIMENSION(KI) :: ZSFU       ! zonal      momentum flux on open sea (and for sea-ice points if merged)(Pa)
REAL, DIMENSION(KI) :: ZSFV       ! meridional momentum flux on open sea (and for sea-ice points if merged)(Pa)
!
REAL, DIMENSION(KI) :: ZSFTH_ICE  ! Heat flux on sea ice 
REAL, DIMENSION(KI) :: ZSFTQ_ICE  ! Sea-ice sublimation flux
REAL, DIMENSION(KI) :: ZSFU_ICE   ! zonal      momentum flux on seaice (Pa)
REAL, DIMENSION(KI) :: ZSFV_ICE   ! meridional momentum flux on seaice (Pa)

REAL, DIMENSION(KI) :: ZHU        ! Near surface relative humidity
REAL, DIMENSION(KI) :: ZQA        ! specific humidity (kg/kg)
REAL, DIMENSION(KI) :: ZEMIS      ! Emissivity at time t
REAL, DIMENSION(KI) :: ZTRAD      ! Radiative temperature at time t
!
REAL, DIMENSION(KI) :: ZSST       ! XSST corrected for anomalously low values (which actually are sea-ice temp)
REAL, DIMENSION(KI) :: ZMASK      ! A mask for diagnosing where seaice exists (or, for coupling_iceflux, may appear)
!
REAL                             :: ZCONVERTFACM0_SLT, ZCONVERTFACM0_DST
REAL                             :: ZCONVERTFACM3_SLT, ZCONVERTFACM3_DST
REAL                             :: ZCONVERTFACM6_SLT, ZCONVERTFACM6_DST
!
INTEGER                          :: ISIZE_WATER  ! number of points with some sea water 
INTEGER                          :: ISIZE_ICE    ! number of points with some sea ice
!
INTEGER                          :: ISWB       ! number of shortwave spectral bands
INTEGER                          :: JSWB       ! loop counter on shortwave spectral bands
INTEGER                          :: ISLT       ! number of sea salt variable
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
! Preliminaries:
!-------------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('COUPLING_SEAFLUX_N',0,ZHOOK_HANDLE)
IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('COUPLING_SEAFLUXN: FATAL ERROR DURING ARGUMENT TRANSFER')        
END IF
!-------------------------------------------------------------------------------------
!
ZEXNA    (:) = XUNDEF
ZEXNS    (:) = XUNDEF
ZU       (:) = XUNDEF
ZV       (:) = XUNDEF
ZWIND    (:) = XUNDEF
ZWORK    (:) = XUNDEF
ZSFTQ    (:) = XUNDEF
ZSFTH    (:) = XUNDEF
ZCD      (:) = XUNDEF    
ZCDN     (:) = XUNDEF
ZCH      (:) = XUNDEF
ZCE      (:) = XUNDEF
ZRI      (:) = XUNDEF
ZHU      (:) = XUNDEF
ZRESA_SEA(:) = XUNDEF
ZUSTAR   (:) = XUNDEF
ZZ0      (:) = XUNDEF
ZZ0H     (:) = XUNDEF
ZQSAT    (:) = XUNDEF
!
ZSFTQ_ICE(:) = XUNDEF
ZSFTH_ICE(:) = XUNDEF
ZCD_ICE  (:) = XUNDEF    
ZCDN_ICE (:) = XUNDEF
ZCH_ICE  (:) = XUNDEF
ZCE_ICE  (:) = XUNDEF
ZRI_ICE  (:) = XUNDEF
ZRESA_SEA_ICE= XUNDEF
ZUSTAR_ICE(:) = XUNDEF
ZZ0_ICE  (:) = XUNDEF
ZZ0H_ICE (:) = XUNDEF
ZQSAT_ICE(:) = XUNDEF
!
ZEMIS    (:) = XUNDEF
ZTRAD    (:) = XUNDEF
ZDIR_ALB (:,:) = XUNDEF
ZSCA_ALB (:,:) = XUNDEF
!
!-------------------------------------------------------------------------------------
!
ZEXNS(:)     = (PPS(:)/XP00)**(XRD/XCPD)
ZEXNA(:)     = (PPA(:)/XP00)**(XRD/XCPD)
!
IF(LCPL_SEA)THEN 
  !Sea currents are taken into account
  ZU(:)=PU(:)-XUMER(:)
  ZV(:)=PV(:)-XVMER(:)
ELSE
  ZU(:)=PU(:)
  ZV(:)=PV(:)        
ENDIF
!
ZWIND(:) = SQRT(ZU(:)**2+ZV(:)**2)
!
PSFTS(:,:) = 0.
!
ZHU = 1.
!
ZQA(:) = PQA(:) / PRHOA(:)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Time evolution
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
TTIME%TIME = TTIME%TIME + PTSTEP
 CALL ADD_FORECAST_TO_DATE_SURF(TTIME%TDATE%YEAR,TTIME%TDATE%MONTH,TTIME%TDATE%DAY,TTIME%TIME)
!
!--------------------------------------------------------------------------------------
! Fluxes over water according to Charnock formulae
!--------------------------------------------------------------------------------------
!
IF (LHANDLE_SIC) THEN 
   ! Flux for sea are computed everywhere
   ISIZE_WATER = SIZE(ZMASK)
   ! Ensure freezing SST values where XSST actually has very low (sea-ice) values (old habits)
   ZSST(:)=MAX(XSST(:), XTTSI) 
   ! Flux over sea-ice will not be computed by next calls, but by coupling_iceflux. Hence :
   ISIZE_ICE   = 0
   ! Flux over sea-ice will be computed by coupling_iceflux anywhere sea-ice could form in one 
   ! time-step (incl. under forcing). ZMASK value is set to 1. on these points
   ZMASK(:)=0.
   WHERE ( XSIC(:) > 0. ) ZMASK(:)=1. 
   ! To be large, assume that seaice may form where SST is < 10°C
   WHERE ( XSST(:) - XTTS <= 10. ) ZMASK(:)=1.
   IF (LINTERPOL_SIC) WHERE (XFSIC(:) > 0. ) ZMASK(:)=1. 
   IF (LINTERPOL_SIT) WHERE (XFSIT(:) > 0. ) ZMASK(:)=1.
ELSE
   ZSST (:) = XSST(:)
   ZMASK(:) = XSST(:) - XTTS
   ISIZE_WATER = COUNT(ZMASK(:)>=0.)
   ISIZE_ICE   = SIZE(XSST) - ISIZE_WATER
ENDIF
!
SELECT CASE (CSEA_FLUX)
  CASE ('DIRECT')
    CALL WATER_FLUX(XZ0,                                              &
                      PTA, ZEXNA, PRHOA, ZSST, ZEXNS, ZQA, PRAIN,     &
                      PSNOW, XTTS,                                    &
                      ZWIND, PZREF, PUREF,                            &
                      PPS, LHANDLE_SIC, ZQSAT,                        &
                      ZSFTH, ZSFTQ, ZUSTAR,                           &
                      ZCD, ZCDN, ZCH, ZRI, ZRESA_SEA, ZZ0H            )  
  CASE ('ITERAT')
    CALL MR98      (XZ0,                                              &
                      PTA, ZEXNA, PRHOA, XSST, ZEXNS, ZQA,            &
                      XTTS,                                           &
                      ZWIND, PZREF, PUREF,                            &
                      PPS, ZQSAT,                                     &
                      ZSFTH, ZSFTQ, ZUSTAR,                           &
                      ZCD, ZCDN, ZCH, ZRI, ZRESA_SEA, ZZ0H            )  

  CASE ('ECUME ','ECUME6')
    CALL ECUME_SEAFLUX(XZ0, ZMASK, ISIZE_WATER, ISIZE_ICE,            &
                      PTA, ZEXNA ,PRHOA, ZSST, XSSS, ZEXNS, ZQA,      &
                      PRAIN, PSNOW,                                   &
                      ZWIND, PZREF, PUREF, PPS, PPA,                  &
                      XICHCE, LPRECIP, LPWEBB, LPWG, NZ0,             &
                      LHANDLE_SIC, ZQSAT, ZSFTH, ZSFTQ, ZUSTAR,       &
                      ZCD, ZCDN, ZCH, ZCE, ZRI, ZRESA_SEA, ZZ0H,      &
                      LPERTFLUX, XPERTFLUX, CSEA_FLUX                 )
  CASE ('COARE3')
    CALL COARE30_SEAFLUX(XZ0, ZMASK, ISIZE_WATER, ISIZE_ICE,          &
                      PTA, ZEXNA ,PRHOA, ZSST, ZEXNS, ZQA, PRAIN,     &
                      PSNOW,                                          &
                      ZWIND, PZREF, PUREF,                            &
                      PPS, LHANDLE_SIC, ZQSAT,                        &
                      ZSFTH, ZSFTQ, ZUSTAR,                           &
                      ZCD, ZCDN, ZCH, ZCE, ZRI, ZRESA_SEA, ZZ0H       )  
END SELECT
!
!-------------------------------------------------------------------------------------
!radiative properties at time t
!-------------------------------------------------------------------------------------
!
ISWB = SIZE(PSW_BANDS)
!
DO JSWB=1,ISWB
  ZDIR_ALB(:,JSWB) = XDIR_ALB(:)
  ZSCA_ALB(:,JSWB) = XSCA_ALB(:)
END DO
!
IF (LHANDLE_SIC) THEN 
   ZEMIS(:) =   (1 - XSIC(:)) * XEMISWAT    + XSIC(:) * XEMISWATICE
   ZTRAD(:) = (((1 - XSIC(:)) * XEMISWAT    * XSST (:)**4 + &
                     XSIC(:)  * XEMISWATICE * XTICE(:)**4)/ &
                     ZEMIS(:)) ** 0.25
ELSE
   ZTRAD(:) = XSST (:)
   ZEMIS(:) = XEMIS(:)
END IF
!
!-------------------------------------------------------------------------------------
!Specific fields for seaice model (when using earth system model or embedded  
!seaice scheme)
!-------------------------------------------------------------------------------------
!
IF(LCPL_SEAICE.OR.LHANDLE_SIC)THEN
  CALL COUPLING_ICEFLUX_n(KI, PTA, ZEXNA, PRHOA, XTICE, ZEXNS, &
                     ZQA, PRAIN, PSNOW, ZWIND, PZREF, PUREF,   &
                     PPS, XSST, XTTS, ZSFTH_ICE, ZSFTQ_ICE,    &  
                     LHANDLE_SIC, ZMASK, ZQSAT_ICE, ZZ0_ICE,   &
                     ZUSTAR_ICE, ZCD_ICE, ZCDN_ICE, ZCH_ICE,   &
                     ZRI_ICE, ZRESA_SEA_ICE, ZZ0H_ICE          )
ENDIF
!
IF (LHANDLE_SIC) CALL COMPLEMENT_EACH_OTHER_FLUX
!
!-------------------------------------------------------------------------------------
! Momentum fluxes over sea or se-ice
!-------------------------------------------------------------------------------------
!
CALL SEA_MOMENTUM_FLUXES(ZCD, ZSFU, ZSFV)
!
! Momentum fluxes over sea-ice if embedded seaice scheme is used
!
IF (LHANDLE_SIC) CALL SEA_MOMENTUM_FLUXES(ZCD_ICE, ZSFU_ICE, ZSFV_ICE)
!
! CO2 flux
!
PSFCO2(:) = 0.0
!
!IF(LCPL_SEA.AND.CSEACO2=='NONE')THEN
!  PSFCO2(:) = XSEACO2(:)
!ELSEIF(CSEACO2=='CST ')THEN
!  PSFCO2 = E * deltapCO2 
!  According to Wanninkhof (medium hypothesis) : 
!  E = 1.13.10^-3 * WIND^2 CO2mol.m-2.yr-1.µatm-1 
!    = 1.13.10^-3 * WIND^2 * Mco2.10^-3 * (1/365*24*3600)
!  deltapCO2 = -8.7 µatm (Table 1 half hypothesis)
PSFCO2(:) = - ZWIND(:)**2 * 1.13E-3 * 8.7 * 44.E-3 / ( 365*24*3600 )
!ENDIF
!
!-------------------------------------------------------------------------------------
! Scalar fluxes:
!-------------------------------------------------------------------------------------
!
IF (NBEQ>0) THEN
  IF (CCH_DRY_DEP == "WES89") THEN

    CALL CH_DEP_WATER  (ZRESA_SEA, ZUSTAR, PTA, ZTRAD,      &
                          PSV(:,NSV_CHSBEG:NSV_CHSEND),       &
                          CSV(NSV_CHSBEG:NSV_CHSEND),         &
                          XDEP(:,1:NBEQ) )  

   PSFTS(:,NSV_CHSBEG:NSV_CHSEND) = - PSV(:,NSV_CHSBEG:NSV_CHSEND)  &
                                               * XDEP(:,1:NBEQ)  
     IF (NAEREQ > 0 ) THEN
        CALL CH_AER_DEP(PSV(:,NSV_AERBEG:NSV_AEREND),&
                          PSFTS(:,NSV_AERBEG:NSV_AEREND),&
                          ZUSTAR,ZRESA_SEA,PTA,PRHOA)     
      END IF

  ELSE
    PSFTS(:,NSV_CHSBEG:NSV_CHSEND) =0.
    IF (NSV_AEREND.GT.NSV_AERBEG)     PSFTS(:,NSV_AERBEG:NSV_AEREND) =0.
  ENDIF
ENDIF
!
IF (NSLTEQ>0) THEN
  ISLT = NSV_SLTEND - NSV_SLTBEG + 1

  CALL COUPLING_SLT_n(           &
       SIZE(ZUSTAR,1),           & !I [nbr] number of sea point
       ISLT,                     & !I [nbr] number of sea salt variables
       ZWIND,                    & !I [m/s] wind velocity
       PSFTS(:,NSV_SLTBEG:NSV_SLTEND) )   
ENDIF
!
IF (NDSTEQ>0) THEN
  CALL DSLT_DEP(PSV(:,NSV_DSTBEG:NSV_DSTEND), PSFTS(:,NSV_DSTBEG:NSV_DSTEND),   &
                ZUSTAR, ZRESA_SEA, PTA, PRHOA, XEMISSIG_DST, XEMISRADIUS_DST,   &
                JPMODE_DST, XDENSITY_DST, XMOLARWEIGHT_DST, ZCONVERTFACM0_DST,  &
                ZCONVERTFACM6_DST, ZCONVERTFACM3_DST, LVARSIG_DST, LRGFIX_DST,  &
                CVERMOD  )  

  CALL MASSFLUX2MOMENTFLUX(         &
    PSFTS(:,NSV_DSTBEG:NSV_DSTEND), & !I/O ![kg/m2/sec] In: flux of only mass, out: flux of moments
    PRHOA,                          & !I [kg/m3] air density
    XEMISRADIUS_DST,                &!I [um] emitted radius for the modes (max 3)
    XEMISSIG_DST,                   &!I [-] emitted sigma for the different modes (max 3)
    NDSTMDE,                        &
    ZCONVERTFACM0_DST,              &
    ZCONVERTFACM6_DST,              &
    ZCONVERTFACM3_DST,              &
    LVARSIG_DST, LRGFIX_DST         )  
ENDIF


IF (NSLTEQ>0) THEN
  CALL DSLT_DEP(PSV(:,NSV_SLTBEG:NSV_SLTEND), PSFTS(:,NSV_SLTBEG:NSV_SLTEND),   &
                ZUSTAR, ZRESA_SEA, PTA, PRHOA, XEMISSIG_SLT, XEMISRADIUS_SLT,   &
                JPMODE_SLT, XDENSITY_SLT, XMOLARWEIGHT_SLT, ZCONVERTFACM0_SLT,  &
                ZCONVERTFACM6_SLT, ZCONVERTFACM3_SLT, LVARSIG_SLT, LRGFIX_SLT,  &
                CVERMOD  )  

  CALL MASSFLUX2MOMENTFLUX(         &
    PSFTS(:,NSV_SLTBEG:NSV_SLTEND), & !I/O ![kg/m2/sec] In: flux of only mass, out: flux of moments
    PRHOA,                          & !I [kg/m3] air density
    XEMISRADIUS_SLT,                &!I [um] emitted radius for the modes (max 3)
    XEMISSIG_SLT,                   &!I [-] emitted sigma for the different modes (max 3)
    NSLTMDE,                        &
    ZCONVERTFACM0_SLT,              &
    ZCONVERTFACM6_SLT,              &
    ZCONVERTFACM3_SLT,              &
    LVARSIG_SLT, LRGFIX_SLT         ) 
ENDIF
!
!-------------------------------------------------------------------------------
! Inline diagnostics at time t for SST and TRAD
!-------------------------------------------------------------------------------
!
CALL DIAG_INLINE_SEAFLUX_n(PTSTEP, PTA, XSST, ZQA, PPA, PPS, PRHOA, PU, &
     PV, PZREF, PUREF, ZCD, ZCDN, ZCH, ZCE, ZRI, ZHU,       &
     XZ0, ZZ0H, ZQSAT, ZSFTH, ZSFTQ, ZSFU, ZSFV,            &
     PDIR_SW, PSCA_SW, PLW, ZDIR_ALB, ZSCA_ALB, XICE_ALB,   &
     ZEMIS, ZTRAD, PRAIN, PSNOW,                            &
     TGLT, XSIC, LHANDLE_SIC, XTICE,                        &
     ZCD_ICE, ZCDN_ICE, ZCH_ICE, ZCE_ICE, ZRI_ICE,          &
     ZZ0_ICE, ZZ0H_ICE, ZQSAT_ICE, ZSFTH_ICE, ZSFTQ_ICE,    &
     ZSFU_ICE, ZSFV_ICE)
!
!-------------------------------------------------------------------------------
! A kind of "average_flux"
!-------------------------------------------------------------------------------
!
IF (LHANDLE_SIC) THEN
   PSFTH  (:) = ZSFTH (:) * ( 1 - XSIC (:)) + ZSFTH_ICE(:) * XSIC(:)
   PSFTQ  (:) = ZSFTQ (:) * ( 1 - XSIC (:)) + ZSFTQ_ICE(:) * XSIC(:)
   PSFU   (:) = ZSFU  (:) * ( 1 - XSIC (:)) +  ZSFU_ICE(:) * XSIC(:)
   PSFV   (:) = ZSFV  (:) * ( 1 - XSIC (:)) +  ZSFV_ICE(:) * XSIC(:)
ELSE
   PSFTH  (:) = ZSFTH (:) 
   PSFTQ  (:) = ZSFTQ (:) 
   PSFU   (:) = ZSFU  (:) 
   PSFV   (:) = ZSFV  (:) 
ENDIF
!
!-------------------------------------------------------------------------------
! IMPOSED SSS OR INTERPOLATED SSS AT TIME t+1
!-------------------------------------------------------------------------------
!
! Daily update Sea surface salinity from monthly data
!
IF (LINTERPOL_SSS .AND. MOD(TTIME%TIME,XDAY) == 0.) THEN
      CALL INTERPOL_SST_MTH(TTIME%TDATE%YEAR,TTIME%TDATE%MONTH,TTIME%TDATE%DAY,'S',XSSS)
ENDIF
!
!-------------------------------------------------------------------------------
! SEA-ICE coupling at time t+1
!-------------------------------------------------------------------------------
!
IF (LHANDLE_SIC) THEN
   IF (LINTERPOL_SIC) THEN
      IF ((MOD(TTIME%TIME,XDAY) == 0.) .OR. (PTIMEC <= PTSTEP )) THEN
      ! Daily update Sea Ice Cover constraint from monthly data
         CALL INTERPOL_SST_MTH(TTIME%TDATE%YEAR,TTIME%TDATE%MONTH,TTIME%TDATE%DAY,'C',XFSIC)
      ENDIF
   ENDIF
   IF (LINTERPOL_SIT) THEN
      IF ((MOD(TTIME%TIME,XDAY) == 0.) .OR. (PTIMEC <= PTSTEP )) THEN
      ! Daily update Sea Ice Thickness constraint from monthly data
         CALL INTERPOL_SST_MTH(TTIME%TDATE%YEAR,TTIME%TDATE%MONTH,TTIME%TDATE%DAY,'H',XFSIT)
      ENDIF
   ENDIF
   IF (CSEAICE_SCHEME=='GELATO') THEN
      CALL SEAICE_GELATO1D_n(HPROGRAM,PTIMEC, PTSTEP, TGLT, XSST, XSSS, XFSIC, XFSIT, XSIC, XTICE, XICE_ALB)
   ELSE
      XTICE=XSST
      XSIC=XFSIC
      XICE_ALB=XALBSEAICE
   ENDIF
   ! Update of cell-averaged albedo, emissivity and radiative 
   ! temperature is done later
ENDIF
!
!-------------------------------------------------------------------------------
! OCEANIC COUPLING, IMPOSED SST OR INTERPOLATED SST AT TIME t+1
!-------------------------------------------------------------------------------
!
IF (LMERCATOR) THEN
   !
   ! Update SST reference profile for relaxation purpose
   IF (LSST_DATA) THEN
      CALL SST_UPDATE(XSEAT_REL(:,NOCKMIN+1), TTIME)
      !
      ! Convert to degree C for ocean model
      XSEAT_REL(:,NOCKMIN+1) = XSEAT_REL(:,NOCKMIN+1) - XTT
   ENDIF
   !
   CALL MOD1D_n(HPROGRAM,PTIME,ZEMIS(:),ZDIR_ALB(:,1:KSW),ZSCA_ALB(:,1:KSW),&
                PLW(:),PSCA_SW(:,1:KSW),PDIR_SW(:,1:KSW),PSFTH(:),          &
                PSFTQ(:),PSFU(:),PSFV(:),PRAIN(:),XSST(:))
   !
ELSEIF(LSST_DATA)THEN 
   !
   ! Imposed SST 
   !
   CALL SST_UPDATE(XSST, TTIME)
   !
ELSEIF (LINTERPOL_SST.AND.MOD(TTIME%TIME,XDAY) == 0.) THEN
   !
   ! Imposed monthly SST 
   !
   CALL INTERPOL_SST_MTH(TTIME%TDATE%YEAR,TTIME%TDATE%MONTH,TTIME%TDATE%DAY,'T',XSST)
   !
ENDIF
!
!-------------------------------------------------------------------------------
!Physical properties see by the atmosphere in order to close the energy budget 
!between surfex and the atmosphere. All variables should be at t+1 but very 
!difficult to do. Maybe it will be done later. However, Ts is at time t+1
!-------------------------------------------------------------------------------
!
IF (LHANDLE_SIC) THEN
   PTSURF (:) = XSST  (:) * ( 1 - XSIC (:)) +     XTICE(:) * XSIC(:)
   PQSURF (:) = ZQSAT (:) * ( 1 - XSIC (:)) + ZQSAT_ICE(:) * XSIC(:)
   ZZ0W   (:) = ( 1 - XSIC(:) ) * 1.0/(LOG(PUREF(:)/ZZ0(:))    **2)  +  &
                      XSIC(:)   * 1.0/(LOG(PUREF(:)/ZZ0_ICE(:))**2)  
   PZ0    (:) = PUREF (:) * EXP ( - SQRT ( 1./  ZZ0W(:) ))
   ZZ0W   (:) = ( 1 - XSIC(:) ) * 1.0/(LOG(PZREF(:)/ZZ0H(:))    **2)  +  &
                      XSIC(:)   * 1.0/(LOG(PZREF(:)/ZZ0H_ICE(:))**2)  
   PZ0H   (:) = PZREF (:) * EXP ( - SQRT ( 1./  ZZ0W(:) ))
ELSE
   PTSURF (:) = XSST  (:) 
   PQSURF (:) = ZQSAT (:) 
   PZ0    (:) = XZ0   (:) 
   PZ0H   (:) = ZZ0H  (:) 
ENDIF
!
!-------------------------------------------------------------------------------
!Radiative properties at time t+1 (see by the atmosphere) in order to close
!the energy budget between surfex and the atmosphere
!-------------------------------------------------------------------------------
!
 CALL UPDATE_RAD_SEA(CSEA_ALB,XSST,PZENITH2,XTTS,XEMIS,  &
                     XDIR_ALB,XSCA_ALB,PDIR_ALB,PSCA_ALB,&
                     PEMIS,PTRAD,LHANDLE_SIC,XTICE,XSIC, &
                     XICE_ALB)
!
!=======================================================================================
!
IF (LHOOK) CALL DR_HOOK('COUPLING_SEAFLUX_N',1,ZHOOK_HANDLE)
!
!=======================================================================================
!
CONTAINS
!
SUBROUTINE SEA_MOMENTUM_FLUXES(PCD, PSFU, PSFV)
!
IMPLICIT NONE
!
REAL, DIMENSION(KI), INTENT(IN)  :: PCD       ! Drag coefficient (on open sea or seaice)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFU      ! zonal momentum flux                   (Pa)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFV      ! meridian momentum flux                (Pa)
!
REAL, DIMENSION(KI) :: ZUSTAR2    ! square of friction velocity (m2/s2)
REAL, DIMENSION(KI) :: ZWORK      ! Work array
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('COUPLING_SEAFLUX_N: SEA_MOMENTUM_FLUXES',0,ZHOOK_HANDLE)
!
ZWORK  (:) = XUNDEF
ZUSTAR2(:) = XUNDEF
!
IF(CIMPLICIT_WIND=='OLD')THEN
! old implicitation (m2/s2)
  ZUSTAR2(:) = (PCD(:)*ZWIND(:)*PPEW_B_COEF(:)) /            &
              (1.0-PRHOA(:)*PCD(:)*ZWIND(:)*PPEW_A_COEF(:))
ELSE
! new implicitation (m2/s2)
  ZUSTAR2(:) = (PCD(:)*ZWIND(:)*(2.*PPEW_B_COEF(:)-ZWIND(:))) /&
              (1.0-2.0*PRHOA(:)*PCD(:)*ZWIND(:)*PPEW_A_COEF(:))
!                   
  ZWORK(:)  = PRHOA(:)*PPEW_A_COEF(:)*ZUSTAR2(:) + PPEW_B_COEF(:)
  ZWORK(:) = MAX(ZWORK(:),0.)
!
  WHERE(PPEW_A_COEF(:)/= 0.)
        ZUSTAR2(:) = MAX( ( ZWORK(:) - PPEW_B_COEF(:) ) / (PRHOA(:)*PPEW_A_COEF(:)), 0.)
  ENDWHERE
!              
ENDIF
!
PSFU = 0.
PSFV = 0.
WHERE (ZWIND(:)>0.)
  PSFU(:) = - PRHOA(:) * ZUSTAR2(:) * ZU(:) / ZWIND(:)
  PSFV(:) = - PRHOA(:) * ZUSTAR2(:) * ZV(:) / ZWIND(:)
END WHERE
!
IF (LHOOK) CALL DR_HOOK('COUPLING_SEAFLUX_N: SEA_MOMENTUM_FLUXES',1,ZHOOK_HANDLE)
!
END SUBROUTINE SEA_MOMENTUM_FLUXES
!
!=======================================================================================
!
SUBROUTINE COMPLEMENT_EACH_OTHER_FLUX
!
! Provide dummy fluxes on places with no open-sea or no sea-ice 
! Allows a smooth computing of CLS parameters in all cases while avoiding 
! having to pack arrays (in routines PARAM_CLS and CLS_TQ)
!
IMPLICIT NONE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('COUPLING_SEAFLUX_N: COMPLEMENT_EACH_OTHER_FLUX',0,ZHOOK_HANDLE)
!
  WHERE (XSIC(:) == 1.)
     ZSFTH=ZSFTH_ICE 
     ZSFTQ=ZSFTQ_ICE 
     ZSFU=ZSFU_ICE
     ZSFV=ZSFV_ICE
     ZQSAT=ZQSAT_ICE
     ZCD=ZCD_ICE
     ZCDN=ZCDN_ICE
     ZCH=ZCH_ICE
     ZCE=ZCE_ICE
     ZRI=ZRI_ICE 
     ZZ0H=ZZ0H_ICE
  END WHERE
  WHERE (XSIC(:) == 0.)
     ZSFTH_ICE=ZSFTH 
     ZSFTQ_ICE=ZSFTQ 
     ZSFU_ICE=ZSFU
     ZSFV_ICE=ZSFV
     ZQSAT_ICE=ZQSAT
     ZCD_ICE=ZCD
     ZCDN_ICE=ZCDN
     ZCH_ICE=ZCH
     ZCE_ICE=ZCE
     ZRI_ICE=ZRI 
     ZZ0H_ICE=ZZ0H
  END WHERE
!
IF (LHOOK) CALL DR_HOOK('COUPLING_SEAFLUX_N: COMPLEMENT_EACH_OTHER_FLUX',1,ZHOOK_HANDLE)
!
END SUBROUTINE COMPLEMENT_EACH_OTHER_FLUX
!
!=======================================================================================
!
END SUBROUTINE COUPLING_SEAFLUX_n
