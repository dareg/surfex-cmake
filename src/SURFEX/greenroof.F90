!     #########
    SUBROUTINE GREENROOF (DTCO, TG, T, TOP, DTGD, TIR, DTI, GB, VD, TV,            &
                          HIMPLICIT_WIND, TPTIME, PTSUN, PPEW_A_COEF, PPEW_B_COEF, &
                          PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,      &
                          PTSTEP, PZREF, PUREF, PTA, PQA, PEXNS, PEXNA, PRHOA,     &
                          PCO2, PPS, PRR, PSR, PZENITH, PSW, PLW, PVMOD,           &
                          PALBNIR_TVEG, PALBVIS_TVEG, PALBNIR_TSOIL, PALBVIS_TSOIL,&                
                          PRN_GREENROOF, PH_GREENROOF, PLE_GREENROOF,              &
                          PGFLUX_GREENROOF, PSFCO2, PEVAP_GREENROOF, PUW_GREENROOF,&
                          PAC_GREENROOF, PQSAT_GREENROOF, PTS_GREENROOF,           &
                          PAC_AGG_GREENROOF, PHU_AGG_GREENROOF, PDEEP_FLUX,        &
                          PRUNOFF_GREENROOF, PDRAIN_GREENROOF, PIRRIG_GREENROOF )  
!   ##################################################################################
!
!!****  *GREENROOF*  
!!
!!    PURPOSE
!!    -------
!!
!!    call the vegetation scheme (ISBA) inside TEB for greenroofs
!!     
!!**  METHOD
!!     ------
!!    based on subroutine "garden" 
!!
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
!!    Based on subroutine "garden"
!!      
!!    AUTHOR
!!    ------
!!
!!      C. de Munck & A. Lemonsu          * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!     Original    09/2011 
!     C. de Munck   02/2013  irrigation (drip irrigation)
!     B. decharme 04/2013 : Variables required in TEB to allow coupling with AROME/ALADIN/ARPEGE
!                           phasing call isba
!                           calculation of vegetation CO2 flux
!                           dummy for water table / surface coupling
!!    P. Samuelsson  10/2014  Introduced dummy variables in call to ISBA for MEB
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_GRID_n, ONLY : GRID_t
USE MODD_SSO_n, ONLY : SSO_t, SSO_INIT
USE MODD_TEB_n, ONLY : TEB_1P_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_IRRIG_n, ONLY : TEB_IRRIG_t
!
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_GR_BIOG_n, ONLY : GR_BIOG_t
USE MODD_SURFEX_n, ONLY : TEB_VEG_DIAG_t
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_t
!
USE MODD_AGRI_n, ONLY : AGRI_t, AGRI_INIT
!
USE MODD_SURF_PAR,             ONLY: XUNDEF
USE MODD_TYPE_DATE_SURF,       ONLY: DATE_TIME
USE MODD_CSTS,                 ONLY: XCPD
!
USE MODI_ISBA
USE MODI_VEGETATION_UPDATE
USE MODI_VEGETATION_EVOL
USE MODI_CARBON_EVOL
USE MODE_THERMOS
USE MODI_ROOF_IMPL_COEF
USE MODI_TEB_IRRIG
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    Declarations of arguments
!
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(GRID_t), INTENT(INOUT) :: TG
TYPE(TEB_1P_t), INTENT(INOUT) :: T
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTGD
TYPE(TEB_IRRIG_t), INTENT(INOUT) :: TIR
!
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTI
TYPE(GR_BIOG_t), INTENT(INOUT) :: GB
TYPE(TEB_VEG_DIAG_t), INTENT(INOUT) :: VD
TYPE(TEB_VEG_t), INTENT(INOUT) :: TV
!
 CHARACTER(LEN=*),     INTENT(IN)  :: HIMPLICIT_WIND   ! wind implicitation option
!                                                     ! 'OLD' = direct
!                                                     ! 'NEW' = Taylor serie, order 1
TYPE(DATE_TIME)     , INTENT(IN)    :: TPTIME             ! current date and time from teb
REAL, DIMENSION(:)  , INTENT(IN)    :: PTSUN              ! solar time      (s from midnight)
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEW_A_COEF        ! implicit coefficients
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEW_B_COEF        ! for wind coupling
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEQ_A_COEF        ! implicit coefficients
REAL, DIMENSION(:)  , INTENT(IN)    :: PPEQ_B_COEF        ! for humidity
REAL, DIMENSION(:)  , INTENT(IN)    :: PPET_A_COEF        ! implicit coefficients
REAL, DIMENSION(:)  , INTENT(IN)    :: PPET_B_COEF        ! for temperature
REAL                , INTENT(IN)    :: PTSTEP             ! time step
REAL, DIMENSION(:)  , INTENT(IN)    :: PZREF              ! height of the first atmospheric level
REAL, DIMENSION(:)  , INTENT(IN)    :: PUREF              ! reference height for the wind
REAL, DIMENSION(:)  , INTENT(IN)    :: PTA                ! temperature at first atm. level 
REAL, DIMENSION(:)  , INTENT(IN)    :: PQA                ! specific humidity at first atm. level
REAL, DIMENSION(:)  , INTENT(IN)    :: PPS                ! pressure at the surface
REAL, DIMENSION(:)  , INTENT(IN)    :: PEXNA              ! Exner function at first atm. level
REAL, DIMENSION(:)  , INTENT(IN)    :: PEXNS              ! surface Exner function
REAL, DIMENSION(:)  , INTENT(IN)    :: PRHOA              ! air density at the lowest level
REAL, DIMENSION(:)  , INTENT(IN)    :: PCO2               ! CO2 concentration in the air    (kg/m3)
REAL, DIMENSION(:)  , INTENT(IN)    :: PRR                ! rain rate
REAL, DIMENSION(:)  , INTENT(IN)    :: PSR                ! snow rate
REAL, DIMENSION(:)  , INTENT(IN)    :: PZENITH            ! solar zenithal angle
REAL, DIMENSION(:)  , INTENT(IN)    :: PSW                ! incoming total solar rad on an horizontal surface
REAL, DIMENSION(:)  , INTENT(IN)    :: PLW                ! atmospheric infrared radiation
REAL, DIMENSION(:)  , INTENT(IN)    :: PVMOD              ! module of horizontal wind near first atm. level
REAL, DIMENSION(:)  , INTENT(IN)    :: PALBNIR_TVEG       ! nearIR  veg tot albedo
REAL, DIMENSION(:)  , INTENT(IN)    :: PALBVIS_TVEG       ! visible veg tot albedo
REAL, DIMENSION(:)  , INTENT(IN)    :: PALBNIR_TSOIL      ! nearIR  soil tot albedo
REAL, DIMENSION(:)  , INTENT(IN)    :: PALBVIS_TSOIL      ! visible soil tot albedo
!
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN_GREENROOF         ! net radiation over greenroofs
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH_GREENROOF          ! sensible heat flux over greenroofs
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE_GREENROOF         ! latent heat flux over greenroofs
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLUX_GREENROOF      ! flux through the greenroofs
REAL, DIMENSION(:)  , INTENT(OUT)   :: PSFCO2                ! flux of greenroof CO2       (m/s*kg_CO2/kg_air)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PEVAP_GREENROOF       ! total evaporation over greenroofs (kg/m2/s)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PUW_GREENROOF         ! friction flux (m2/s2)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_GREENROOF         ! greenroof aerodynamical conductance
REAL, DIMENSION(:)  , INTENT(OUT)   :: PQSAT_GREENROOF       ! saturation humidity
REAL, DIMENSION(:)  , INTENT(INOUT) :: PTS_GREENROOF         ! greenroof radiative surface temp. (snow free)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_AGG_GREENROOF     ! aggreg. aeodynamic resistance for greenroofs for latent heat flux
REAL, DIMENSION(:)  , INTENT(OUT)   :: PHU_AGG_GREENROOF     ! aggreg. relative humidity for greenroofs for latent heat flux
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDEEP_FLUX            ! Heat Flux at the bottom layer of the greenroof
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRUNOFF_GREENROOF     ! greenroof surface runoff
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDRAIN_GREENROOF      ! greenroof total (vertical) drainage
REAL, DIMENSION(:)  , INTENT(OUT)   :: PIRRIG_GREENROOF      ! greenroof summer irrigation rate
!
!
!*      0.2    Declarations of local variables
!
TYPE(SSO_t) :: YGRSS
TYPE(AGRI_t) :: YAG
!
REAL, DIMENSION(SIZE(PPS)) :: ZDIRCOSZW           ! orography slope cosine (=1 in TEB)
REAL, DIMENSION(SIZE(PPS),TV%O%NNBIOMASS) :: ZRESP_BIOMASS_INST       ! instantaneous biomass respiration (kgCO2/kgair m/s)
REAL, DIMENSION(SIZE(PPS)) :: ZUSTAR
!
!  temperatures
!
REAL, DIMENSION(SIZE(PPS)) :: ZTA ! estimate of air temperature at future time
!                                 ! step as if modified by ISBA flux alone.
!
!  surfaces relative fractions
!  for flood
REAL, DIMENSION(SIZE(PPS)) :: ZEMISF
!
!  variables for deep soil temperature
REAL, DIMENSION(SIZE(PPS)) :: ZTDEEP_A
!
! Dummy variables for MEB:
REAL, DIMENSION(SIZE(PPS)) :: ZP_MEB_SCA_SW, ZPALPHAN, ZZ0G_WITHOUT_SNOW, &
                              ZZ0_MEBV, ZZ0H_MEBV, ZZ0EFF_MEBV, ZZ0_MEBN, &
                              ZZ0H_MEBN, ZZ0EFF_MEBN
!
!
INTEGER                    :: ILU
LOGICAL :: GUPDATED, GALB
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*      1.     various initialisations
!              -----------------------
!
IF (LHOOK) CALL DR_HOOK('GREENROOF',0,ZHOOK_HANDLE)
ILU = SIZE(PPS)
!
ZDIRCOSZW = 1.
!
 CALL SSO_INIT(YGRSS)
!
 CALL AGRI_INIT(YAG)
!
!* automatic summer irrigation 
!
PIRRIG_GREENROOF(:) = 0.
!
!* deep soil implicitation with roof
!
 CALL ROOF_IMPL_COEF(T, PTSTEP, ZTDEEP_A, TV%IP%XTDEEP)
!
!-------------------------------------------------------------------------------
!
!*      9.     Treatment of green areas
!              ------------------------
!
!radiative temperature diagnostic
!-------------------------------
!
!*      9.1    Summer irrigation 
!              ------------------
!
!* irrigation automatique de type goutte à goutte (arrosage du sol seulement)
!
CALL TEB_IRRIG(TIR%LPAR_GR_IRRIG, PTSTEP, TPTIME%TDATE%MONTH, PTSUN, &
               TIR%XGR_START_MONTH, TIR%XGR_END_MONTH, TIR%XGR_START_HOUR,   &
               TIR%XGR_END_HOUR, TIR%XGR_24H_IRRIG, PIRRIG_GREENROOF     )
!
! --------------------------------------------------------------------------------------
! Vegetation update (in case of non-interactive vegetation):
! --------------------------------------------------------------------------------------
!
TV%I%TTIME = TPTIME
TV%O%LECOCLIMAP = (.NOT. TV%O%LPAR)
!
GUPDATED=.FALSE.
GALB = .FALSE. 
!print*,TV%O%CPHOTO,TV%M%T%CUR%XLAI(1,1)
IF (TV%O%CPHOTO=='LAI'.OR.TV%O%CPHOTO=='LST'.OR.TV%O%CPHOTO=='NIT'.OR.TV%O%CPHOTO=='NCB') GALB = .TRUE.
!
CALL VEGETATION_UPDATE(DTCO, DTI, TG%NDIM, TV%O, TV%M%T%CUR, TV%M%M, TV%M%I, TV%M%A, &
                       PTSTEP, TV%I%TTIME, TOP%XCOVER, TOP%LCOVER, .FALSE.,'GNR',  &
                       GALB, YGRSS, GUPDATED, OABSENT=(T%XGREENROOF==0.)     )
!
!print*,'update',TV%M%T%CUR%XLAI(1,1)
!*      9.2    Call ISBA for greenroofs
!              ------------------------
!
VD%DP%XZ0(:) = TV%M%T%CUR%XZ0(:,1)
VD%DP%XZ0H(:) = TV%M%T%CUR%XZ0(:,1) / TV%M%X%XZ0_O_Z0H(:,1)
!
VD%DP%XZ0EFF(:) =  TV%M%T%CUR%XZ0(:,1)
!
ALLOCATE(GB%XIACAN(SIZE(PPS),SIZE(TV%IP%XABC),1))
!
 CALL ISBA(TV%O, TV%M%X, TV%M%T%CUR, TV%M%M, TV%M%I, TV%P, TV%IP, TV%I, TV%R%CUR, &
           TG, YAG, VD%D, VD%DP, VD%E, VD%EP, VD%M, TV%R%CUR%TSNOW%SCHEME, TPTIME,&
           TV%IP%XPOI, TV%IP%XABC, GB%XIACAN(:,:,1), .FALSE., PTSTEP,             &
           HIMPLICIT_WIND, PZREF, PUREF, ZDIRCOSZW, PTA, PQA, PEXNA, PRHOA, PPS,  &
           PEXNS,  PRR, PSR, PZENITH, ZP_MEB_SCA_SW, PSW, PLW, PVMOD, PPEW_A_COEF,&
           PPEW_B_COEF, PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,       &
           PALBNIR_TVEG, PALBVIS_TVEG, PALBNIR_TSOIL, PALBVIS_TSOIL, ZPALPHAN,    &
           ZZ0G_WITHOUT_SNOW, ZZ0_MEBV, ZZ0H_MEBV, ZZ0EFF_MEBV, ZZ0_MEBN,         &
           ZZ0H_MEBN, ZZ0EFF_MEBN, ZTDEEP_A, PCO2, TV%I%XFFG(:,1), TV%I%XFFV(:,1),&
           ZEMISF, ZUSTAR, PAC_AGG_GREENROOF, PHU_AGG_GREENROOF,                  &
           ZRESP_BIOMASS_INST, PDEEP_FLUX, PIRRIG_GREENROOF )
!
PTS_GREENROOF(:) = VD%DP%XTSRAD(:)
PRUNOFF_GREENROOF(:) = VD%EP%XRUNOFF(:)
PDRAIN_GREENROOF (:) = VD%EP%XDRAIN(:)
!
IF (TV%R%CUR%TSNOW%SCHEME=='3-L' .OR. TV%R%CUR%TSNOW%SCHEME=='CRO') &
        TV%R%CUR%TSNOW%TS(:,1) = VD%M%XSNOWTEMP(:,1)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Diagnostic of respiration carbon fluxes and soil carbon evolution
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Vegetation evolution for interactive LAI
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
IF (TV%O%CPHOTO=='LAI' .OR. TV%O%CPHOTO=='LST' .OR. TV%O%CPHOTO=='NIT') THEN
        !print*,'LAI evol ',TV%M%T%CUR%XLAI(1,1)
  CALL VEGETATION_EVOL(TV%O, TV%IP, TV%M%X, TV%M%T%CUR, TV%M%A, TV%M%I, TV%R%CUR, &
                       .FALSE., PTSTEP, TPTIME%TDATE%MONTH, TPTIME%TDATE%DAY,     &
                       TPTIME%TIME, TG%XLAT, PRHOA, PCO2, YGRSS, ZRESP_BIOMASS_INST )
        !print*,'LAI evol2 ',TV%M%T%CUR%XLAI(1,1)        
END IF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
PSFCO2    (:)=0.
VD%EP%XRESP_ECO (:)=0.
VD%EP%XRESP_AUTO(:)=0.
!
!
IF (TV%O%CPHOTO/='NON' .AND. TV%O%CRESPSL/='NON' .AND. ANY(TV%M%T%CUR%XLAI(:,1)/=XUNDEF)) THEN
  ! faire intervenir le type de vegetation du greenroof ? (CTYP_GR)
  CALL CARBON_EVOL(TV%O, TV%P, TV%IP, TV%M%X, TV%M%T%CUR, TV%R%CUR, VD%EP, &
                   PTSTEP, PRHOA, ZRESP_BIOMASS_INST  )
  ! calculation of vegetation CO2 flux
  ! Positive toward the atmosphere
  PSFCO2(:) = VD%EP%XRESP_ECO(:) - VD%EP%XGPP(:)
END IF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!*      9.     Fields required for TEB
!              -----------------------
!
! energy balance
!
 PRN_GREENROOF    (:) = VD%DP%XRN    (:)
 PH_GREENROOF     (:) = VD%DP%XH     (:)
 PLE_GREENROOF    (:) = TV%R%CUR%XLE (:,1)
 PGFLUX_GREENROOF (:) = VD%DP%XGFLUX (:)
 PEVAP_GREENROOF  (:) = VD%DP%XEVAP  (:)
!
!
! Estimate of green area aerodynamic conductance recomputed from heat flux,
! surface (radiative) temp. and forcing air temperature (estimated at future time step)
 ZTA = PPET_B_COEF + PPET_A_COEF * PH_GREENROOF
 PAC_GREENROOF = 0.
 WHERE (PTS_GREENROOF /= ZTA)
   PAC_GREENROOF(:)   = MAX(PH_GREENROOF(:) / XCPD / PRHOA(:) / (PTS_GREENROOF - ZTA) , 0.)
 ENDWHERE
!
! Humidity of saturation for green areas
 PQSAT_GREENROOF(:) = QSAT(TV%R%CUR%XTG(:,1,1),PPS(:))
!
!* friction flux
  PUW_GREENROOF(:)    = -ZUSTAR(:)**2
IF (LHOOK) CALL DR_HOOK('GREENROOF',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE GREENROOF
