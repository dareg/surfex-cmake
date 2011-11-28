!     ###############################################################################
SUBROUTINE COUPLING_FLAKE_n(HPROGRAM, HCOUPLING,                                         &
                 PTSTEP, KYEAR, KMONTH, KDAY, PTIME, KI, KSV, KSW, PTSUN, PZENITH, PAZIM,    &
                 PZREF, PUREF, PZS, PU, PV, PQA, PTA, PRHOA, PSV, PCO2, HSV,                 &
                 PRAIN, PSNOW, PLW, PDIR_SW, PSCA_SW, PSW_BANDS, PPS, PPA,                   &
                 PSFTQ, PSFTH, PSFTS, PSFCO2, PSFU, PSFV,                                    &
                 PTRAD, PDIR_ALB, PSCA_ALB, PEMIS,                                           &
                 PPEW_A_COEF, PPEW_B_COEF,                                                   &
                 PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,                         &
                 HTEST                                                                       )  
!     ###############################################################################

!
!!****  *COUPLING_FLAKE_n * - Driver for FLAKE scheme for lakes
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
!!      V. Masson   05/2009 Implicitation of momentum fluxes
!!      B. Decharme 01/2010 Add XTT in water_flux
!!      V. Masson   11/2011 Ch limited to 1.E-7 in all cases and Cd coming from
!!                          Flake_interface routine if computed by flake
!!------------------------------------------------------------------
!
!
USE MODD_CSTS,       ONLY : XRD, XCPD, XP00, XLVTT, XKARMAN, XTT
! USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE MODD_FLAKE_n,  ONLY : TTIME         ,  &
                          ! Flake parameters
                            XEMIS         ,  &! 
                            XWATER_DEPTH  ,  &! Lake depth (m)
                            XWATER_FETCH  ,  &! Lake fetch (m)
                            XT_BS         ,  &! Temperature at the outer edge of the thermally
                                            !       active layer of the bottom sediments [K]
                            XDEPTH_BS     ,  &! Depth of the thermally active layer of the 
                                            !       bottom sediments [m]
                            XCOR          ,  &! The Coriolis parameter [s^{-1}]
                            XALB_WATER    ,  &! Water surface albedo
                            XALB_ICE      ,  &! Ice surface albedo
                            XALB_SNOW     ,  &! Snow surface albedo
                            XEXTCOEF_WATER,  &! Extinction coefficient for the water [m^{-1}]
                            XEXTCOEF_ICE  ,  &! Extinction coefficient for the ice [m^{-1}]
                            XEXTCOEF_SNOW ,  &! Extinction coefficient for the snow [m^{-1}] 
                          ! FLake variables
                            XT_SNOW       ,  &! Temperature at the air-snow interface [K]    
                            XT_ICE        ,  &! Temperature at the snow-ice or air-ice 
                                            !         interface [K]
                            XT_MNW        ,  &! Mean temperature of the water column [K]
                            XT_WML        ,  &! Mixed-layer temperature [K]
                            XT_BOT        ,  &! Temperature at the water-bottom sediment 
                                            !         interface [K]
                            XT_B1         ,  &! Temperature at the bottom of the upper layer 
                                            !         of the sediments [K]
                            XCT           ,  &! Shape factor (thermocline)
                            XH_SNOW       ,  &! Snow thickness [m]
                            XH_ICE        ,  &! Ice thickness [m]
                            XH_ML         ,  &! Thickness of the mixed-layer [m]
                            XH_B1         ,  &! Thickness of the upper layer of bottom sediments [m]                                   
                            XTS           ,  &! surface temperature  (K)                            
                            XZ0           , &
                            XUSTAR        ,  &! air friction velocity (m/s)
                            LSEDIMENTS    ,  &! Switch, .TRUE. -> use the bottom-sediment scheme 
                            LWATFLX           ! Switch, .TRUE. -> compute the surface fluxes with water_flux  
!                          
!salgado - keep the same ch_ routines and modules used in watflux_n
USE MODD_CH_WATFLUX_n, ONLY : CSV, CCH_DRY_DEP, XDEP, NBEQ, NSV_CHSBEG, NSV_CHSEND,&
                                NSV_DSTBEG, NSV_DSTEND, NAEREQ, NDSTEQ, NSLTEQ, &
                                NSV_AERBEG, NSV_AEREND, NSV_SLTBEG, NSV_SLTEND  
!
! 
USE MODI_WATER_FLUX
USE MODI_ADD_FORECAST_TO_DATE_SURF
USE MODI_DIAG_INLINE_FLAKE_n 
USE MODI_DIAG_MISC_FLAKE_n
USE MODI_CH_AER_DEP
USE MODI_CH_DEP_WATER
USE MODI_DST_DEP
USE MODI_SLT_DEP
USE MODE_DST_SURF
USE MODE_SLT_SURF
USE MODD_SLT_n,       ONLY: XEMISRADIUS_SLT,XEMISSIG_SLT
USE MODD_DST_n,       ONLY: XEMISRADIUS,XEMISSIG

USE MODI_WIND_THRESHOLD
!
USE MODE_THERMOS
!


! 

!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
USE MODI_FLAKE_INTERFACE
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=6),    INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
CHARACTER(LEN=1),    INTENT(IN)  :: HCOUPLING ! type of coupling
                                              ! 'E' : explicit
                                              ! 'I' : implicit
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
REAL, DIMENSION(KI), INTENT(IN)  :: PZENITH   ! zenithal angle       (radian from the vertical)
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
!
REAL, DIMENSION(KI), INTENT(OUT) :: PSFTH     ! flux of heat                          (W/m2)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFTQ     ! flux of water vapor                   (kg/m2/s)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFU      ! zonal momentum flux                   (Pa)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFV      ! meridian momentum flux                (Pa)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFCO2    ! flux of CO2                           (kg/m2/s)
REAL, DIMENSION(KI,KSV),INTENT(OUT):: PSFTS   ! flux of scalar var.                   (kg/m2/s)
!
REAL, DIMENSION(KI), INTENT(OUT) :: PTRAD     ! radiative temperature                 (K)
REAL, DIMENSION(KI,KSW),INTENT(OUT):: PDIR_ALB! direct albedo for each spectral band  (-)
REAL, DIMENSION(KI,KSW),INTENT(OUT):: PSCA_ALB! diffuse albedo for each spectral band (-)
REAL, DIMENSION(KI), INTENT(OUT) :: PEMIS     ! emissivity                            (-)
!
REAL, DIMENSION(KI), INTENT(IN) :: PPEW_A_COEF! implicit coefficients
REAL, DIMENSION(KI), INTENT(IN) :: PPEW_B_COEF! needed if HCOUPLING='I'
REAL, DIMENSION(KI), INTENT(IN) :: PPET_A_COEF
REAL, DIMENSION(KI), INTENT(IN) :: PPEQ_A_COEF
REAL, DIMENSION(KI), INTENT(IN) :: PPET_B_COEF
REAL, DIMENSION(KI), INTENT(IN) :: PPEQ_B_COEF
CHARACTER(LEN=2),    INTENT(IN) :: HTEST ! must be equal to 'OK'
!
!*      0.2    declarations of local variables
!
INTEGER                     :: ISWB   ! number of shortwave spectral bands
INTEGER                     :: JSWB   ! loop counter on shortwave spectral bands
!         
REAL, DIMENSION(SIZE(PTA))  :: ZEXNA  ! Exner function at forcing level
REAL, DIMENSION(SIZE(PTA))  :: ZEXNS  ! Exner function at surface level
!
REAL, DIMENSION(SIZE(PTA))  :: ZWIND  ! Wind
REAL, DIMENSION(SIZE(PTA))  :: ZSW    ! Solar radiation flux at the surface (W/m2) 
REAL, DIMENSION(SIZE(PTA))  :: ZQA    ! Air specific humidity (kg/kg)
!
REAL, DIMENSION(SIZE(PTA))  :: ZUSTAR ! friction velocity (m/s)
REAL, DIMENSION(SIZE(PTA))  :: ZSFM   ! flux of momentum (Pa)
!?? Not used, for the moment
REAL, DIMENSION(SIZE(PTA))  :: ZRESA_WATER ! aerodynamical resistance
!
!salgado only for inline diagnostics - not used for the moment
!                                      flake don't have it
REAL, DIMENSION(SIZE(PTA))  :: ZCD    ! Drag coefficient
REAL, DIMENSION(SIZE(PTA))  :: ZCDN   ! Neutral Drag coefficient
REAL, DIMENSION(SIZE(PTA))  :: ZCH    ! Heat transfer coefficient
REAL, DIMENSION(SIZE(PTA))  :: ZRI    ! Richardson number
REAL, DIMENSION(SIZE(PTA))  :: ZHU    ! Near surface relative humidity
REAL, DIMENSION(SIZE(PTA))  :: ZZ0    ! roughness length
REAL, DIMENSION(SIZE(PTA))  :: ZZ0H   ! heat roughness length
REAL, DIMENSION(SIZE(PTA))  :: ZQSAT  ! humidity at saturation
REAL, DIMENSION(SIZE(PTA))  :: ZTSTEP ! time-step

INTEGER                     :: ILUOUT ! output logical unit
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
! Preliminaries:
!-------------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('COUPLING_FLAKE_N',0,ZHOOK_HANDLE)
IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('COUPLING_FLAKEN: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF
!-------------------------------------------------------------------------------------
! Variables needed by flake:
!-------------------------------------------------------------------------------------
!
ZTSTEP(:) = PTSTEP
!
ZEXNS(:)     = (PPS(:)/XP00)**(XRD/XCPD)
ZEXNA(:)     = (PPA(:)/XP00)**(XRD/XCPD)
!
!
ZWIND(:) = SQRT(PU(:)**2+PV(:)**2)
!
ZWIND(:) = WIND_THRESHOLD(ZWIND(:),PUREF)
!
ZSW(:) = 0.
DO JSWB=1,KSW
  ZSW(:) = ZSW(:) + PDIR_SW(:,JSWB)+PSCA_SW(:,JSWB)
END DO
ZQA(:) = PQA/PRHOA
!
PSFTS(:,:) = 0.
!
ZHU = 1.
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Time evolution
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
TTIME%TIME = TTIME%TIME + PTSTEP
CALL ADD_FORECAST_TO_DATE_SURF(TTIME%TDATE%YEAR,TTIME%TDATE%MONTH,TTIME%TDATE%DAY,TTIME%TIME)
!
!--------------------------------------------------------------------------------------
PSFU = 0.
PSFV = 0.
ZSFM = 0.
IF (LWATFLX) then ! Call water_flux to compute fluxes over water
    CALL WATER_FLUX(XZ0,                                        &
                  PTA, ZEXNA, PRHOA, XTS, ZEXNS, PQA,PRAIN, PSNOW,&
                  XTT, ZWIND, PZREF, PUREF,                       &
                  PPS, ZQSAT,                                     &
                  PSFTH, PSFTQ, ZUSTAR,                           &
                  ZCD, ZCDN, ZCH, ZRI, ZRESA_WATER, ZZ0H          )  
    WHERE (ZWIND(:)>0.)
!      ZSFM(:) = - PRHOA(:) * ZUSTAR(:)**2 
      ZUSTAR(:) = SQRT(  (ZCD(:)*ZWIND(:)*PPEW_B_COEF(:))/                &
                          (1.0-PRHOA(:)*ZCD(:)*ZWIND(:)*PPEW_A_COEF(:))     )  
      ZSFM(:) = - PRHOA(:) * ZUSTAR(:)**2 
      PSFU(:) = ZSFM(:) * PU(:) / ZWIND(:)
      PSFV(:) = ZSFM(:) * PV(:) / ZWIND(:)
    END WHERE
  ! PSFTQ become temporarly the flux of heat flux (W/m2)
    PSFTQ = PSFTQ * XLVTT
!
ENDIF
!
!--------------------------------------------------------------------------------------
! Call FLake 
! to compute Fluxes over water if LWATFLX = .FALSE.
! to actualize FLake variables, namely water surface temperature
!--------------------------------------------------------------------------------------
!
ZZ0 = XZ0
!
CALL FLAKE_INTERFACE( KI, &
! Atmospheric forcing:
                       PSNOW, ZSW, PLW, PUREF, PZREF, ZWIND, PTA, ZQA, PPS, &
! Constant parameters
                       XWATER_DEPTH, XWATER_FETCH, XDEPTH_BS, XT_BS, XCOR,  &
                       ZTSTEP,                                              &
! Parameters that may change (constants for the moment)
                       XALB_WATER, XALB_ICE, XALB_SNOW, XEXTCOEF_WATER,     &
                       XEXTCOEF_ICE, XEXTCOEF_SNOW,                         &
! Flake variables
                       XT_SNOW, XT_ICE, XT_MNW, XT_WML, XT_BOT, XT_B1, XCT, &
                       XH_SNOW, XH_ICE, XH_ML, XH_B1, XTS,                  &
! Surface heat and momentum fluxes
                       PSFTH, PSFTQ, ZSFM, ZZ0, ZZ0H, ZRI, XUSTAR, ZCD,     &
! Flags              
                       LSEDIMENTS, LWATFLX, PPEW_A_COEF, PPEW_B_COEF, PRHOA )   
!
IF (.NOT. LWATFLX) then 
    XZ0 = ZZ0
ENDIF
!
!-------------------------------------------------------------------------------------
! Outputs:
!-------------------------------------------------------------------------------------
!
! Momentum fluxes
!
IF (.NOT. LWATFLX) THEN
  WHERE (ZWIND(:)>0.)
    PSFU(:) = ZSFM(:) * PU(:) / ZWIND(:)
    PSFV(:) = ZSFM(:) * PV(:) / ZWIND(:)
  END WHERE
  !
  ! 
  ! ZUSTAR and ZRESA_WATER are not in Flake but are needed to the ch_* routines
  !
  ZUSTAR(:)       = SQRT (ABS(ZSFM(:))/PRHOA(:))
  ZEXNS (:)       = (PPS(:)/XP00)**(XRD/XCPD)
  ZEXNA (:)       = (PPA(:)/XP00)**(XRD/XCPD)
  ZRESA_WATER=2.E4
  WHERE (PSFTH/=0.)
  ZRESA_WATER (:) = XCPD * PRHOA(:) * (XTS(:) - PTA(:) * ZEXNS(:)/ZEXNA(:)) &
                     / (PSFTH(:) * ZEXNS(:))  
  END WHERE
  !                               
ENDIF
! flux of water vapor (kg/m2/s)
PSFTQ = PSFTQ / XLVTT
!
! CO2 flux
!
PSFCO2(:)       =  0.0    ! Assumes no CO2 emission over water bodies
!
!radiative properties
!
ISWB = SIZE(PSW_BANDS)
!
DO JSWB=1,ISWB 
  PDIR_ALB(:,JSWB) = XALB_WATER(:) !salgado - to be changed (to account for snow or ice)
  PSCA_ALB(:,JSWB) = XALB_WATER(:)
END DO
PEMIS  = XEMIS
PTRAD  = XTS
!
!-------------------------------------------------------------------------------------
! Scalar fluxes:
!-------------------------------------------------------------------------------------
!
!
!salgado The scalar fluxes are computed as in watflux
IF (NBEQ>0) THEN
  IF (CCH_DRY_DEP == "WES89") THEN
    CALL CH_DEP_WATER  (ZRESA_WATER, ZUSTAR, PTA, PTRAD,      &
                          PSV(:,NSV_CHSBEG:NSV_CHSEND),       &
                          CSV(NSV_CHSBEG:NSV_CHSEND),         &
                          XDEP(:,1:NBEQ) )  

   PSFTS(:,NSV_CHSBEG:NSV_CHSEND) = - PSV(:,NSV_CHSBEG:NSV_CHSEND)  &
                                               * XDEP(:,1:NBEQ)  
     IF (NAEREQ > 0 ) THEN
        CALL CH_AER_DEP(PSV(:,NSV_AERBEG:NSV_AEREND),&
                          PSFTS(:,NSV_AERBEG:NSV_AEREND),&
                          ZUSTAR,ZRESA_WATER,PTA,PRHOA)     
      END IF

  ELSE
    PSFTS(:,NSV_CHSBEG:NSV_CHSEND) =0.
    IF(NSV_AERBEG.LT.NSV_AEREND) PSFTS(:,NSV_AERBEG:NSV_AEREND) =0.
  ENDIF
ENDIF

IF (NDSTEQ>0) THEN
     CALL DST_DEP(PSV(:,NSV_DSTBEG:NSV_DSTEND),&
                    PSFTS(:,NSV_DSTBEG:NSV_DSTEND),&
                    ZUSTAR,ZRESA_WATER,PTA,PRHOA)     
     CALL MASSFLUX2MOMENTFLUX(          &
            PSFTS(:,NSV_DSTBEG:NSV_DSTEND), &!I/O ![kg/m2/sec] In: flux of only mass, out: flux of moments
            PRHOA,                         &!I [kg/m3] air density
            XEMISRADIUS,                   &!I [um] emitted radius for the modes (max 3)
            XEMISSIG                      &!I [-] emitted sigma for the different modes (max 3)
            )  
ENDIF


IF (NSLTEQ>0) THEN
     CALL SLT_DEP(PSV(:,NSV_SLTBEG:NSV_SLTEND),&
                    PSFTS(:,NSV_SLTBEG:NSV_SLTEND),&
                    ZUSTAR,ZRESA_WATER,PTA,PRHOA)     
     CALL MASSFLUX2MOMENTFLUX_SLT(      &
            PSFTS(:,NSV_SLTBEG:NSV_SLTEND), &!I/O ![kg/m2/sec] In: flux of only mass, out: flux of moments
            PRHOA,                         &!I [kg/m3] air density
            XEMISRADIUS_SLT,                &!I [um] emitted radius for the modes (max 3)
            XEMISSIG_SLT                   &!I [-] emitted sigma for the different modes (max 3)
            )  

ENDIF

!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Inline diagnostics
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
IF (.NOT.LWATFLX) THEN  !compute some variables not present in FLake code
  ZCH = 1.E-7
!
  WHERE (ABS((XTS(:) - PTA(:) * ZEXNS(:)/ZEXNA(:))) > 1.E-2)
     ZCH = MAX(1.E-7,PSFTH / (XCPD * PRHOA(:) * ZWIND(:) * (XTS(:) - PTA(:) * ZEXNS(:)/ZEXNA(:))) * ZEXNS(:))
  END WHERE
!
!
  ZCDN = (XKARMAN/LOG(PUREF(:)/XZ0(:)))**2
!
  ZQSAT(:) = QSAT(XTS(:),PPS(:))
ENDIF
!
CALL DIAG_INLINE_FLAKE_n(PTA, XTS, ZQA, PPA, PPS, PRHOA, PU, PV, PZREF, PUREF,  &
                             ZCD, ZCDN, ZCH, ZRI, ZHU,                           &
                             XZ0, ZZ0H, &
                             ZQSAT,            &
                             PSFTH, PSFTQ, PSFU, PSFV,                              &
                             PDIR_SW, PSCA_SW, PLW,                                 &
                             PDIR_ALB, PSCA_ALB, PEMIS, PTRAD                       )  
!
!-------------------------------------------------------------------------------------
!
CALL DIAG_MISC_FLAKE_n(XT_WML,XT_BOT,XH_ML,XCT,XWATER_DEPTH)
IF (LHOOK) CALL DR_HOOK('COUPLING_FLAKE_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE COUPLING_FLAKE_n
