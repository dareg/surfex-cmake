!     #########
    SUBROUTINE GARDEN (DTCO, G, T, TOP, TIR, DTI, GB, VD, GDO, S, K, P, PEK,  &
                       HIMPLICIT_WIND, TPTIME, PTSUN, PPEW_A_COEF, PPEW_B_COEF, &
                       PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,      &
                       PTSTEP, PZREF, PTA, PQA, PEXNS, PRHOA,                   &
                       PCO2, PPS, PRR, PSR, PZENITH, PSW, PLW, PVMOD,           &
                       PALBNIR_TVEG, PALBVIS_TVEG, PALBNIR_TSOIL, PALBVIS_TSOIL,&                
                       PRN, PH, PLE, PGFLUX, PSFCO2, PEVAP, PUW, PAC,PQSAT,     &
                       PAC_AGG, PHU_AGG, PRUNOFF, PDRAIN, PIRRIG         )  
!   ##########################################################################
!
!!****  *GARDEN*  
!!
!!    PURPOSE
!!    -------
!
!!call the vegetation scheme (ISBA) inside TEB
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
!!      A. Lemonsu          * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!    Original    05/2009
!     B. decharme 04/2013 : variables for surf/atm coupling
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
!
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_GR_BIOG_n, ONLY : GR_BIOG_t
USE MODD_SURFEX_n, ONLY : TEB_VEG_DIAG_t
USE MODD_TEB_IRRIG_n, ONLY : TEB_IRRIG_t
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n, ONLY : ISBA_S_t, ISBA_K_t, ISBA_P_t, ISBA_PE_t
!
USE MODD_AGRI_n, ONLY : AGRI_t,AGRI_INIT
!
USE MODD_TYPE_DATE_SURF,    ONLY: DATE_TIME
USE MODD_SURF_PAR,          ONLY: XUNDEF
USE MODD_CSTS,              ONLY: XCPD
!
!
USE MODI_ISBA
USE MODI_VEGETATION_UPDATE
USE MODE_THERMOS
!
USE MODI_FLAG_TEB_VEG_n
USE MODI_CARBON_EVOL
USE MODI_VEGETATION_EVOL
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
TYPE(GRID_t), INTENT(INOUT) :: G
TYPE(TEB_1P_t), INTENT(INOUT) :: T
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
!
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTI
TYPE(GR_BIOG_t), INTENT(INOUT) :: GB
TYPE(TEB_VEG_DIAG_t), INTENT(INOUT) :: VD
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: GDO
TYPE(ISBA_S_t), INTENT(INOUT) :: S
TYPE(ISBA_K_t), INTENT(INOUT) :: K
TYPE(ISBA_P_t), INTENT(INOUT) :: P
TYPE(ISBA_PE_t), INTENT(INOUT) :: PEK
!
TYPE(TEB_IRRIG_t), INTENT(INOUT) :: TIR
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
REAL, DIMENSION(:)  , INTENT(IN)    :: PZREF              ! height of atm. var. near the road
REAL, DIMENSION(:)  , INTENT(IN)    :: PTA                ! temp. near the road
REAL, DIMENSION(:)  , INTENT(IN)    :: PQA                ! hum. near the road
REAL, DIMENSION(:)  , INTENT(IN)    :: PPS                ! pressure at the surface
REAL, DIMENSION(:)  , INTENT(IN)    :: PEXNS              ! surface exner function
REAL, DIMENSION(:)  , INTENT(IN)    :: PRHOA              ! air density at the lowest level
REAL, DIMENSION(:)  , INTENT(IN)    :: PCO2               ! CO2 concentration in the air    (kg/m3)
REAL, DIMENSION(:)  , INTENT(IN)    :: PRR                ! rain rate
REAL, DIMENSION(:)  , INTENT(IN)    :: PSR                ! snow rate
REAL, DIMENSION(:)  , INTENT(IN)    :: PZENITH            ! solar zenithal angle
REAL, DIMENSION(:)  , INTENT(IN)    :: PSW                ! incoming total solar rad on an horizontal surface
REAL, DIMENSION(:)  , INTENT(IN)    :: PLW                ! atmospheric infrared radiation
REAL, DIMENSION(:)  , INTENT(IN)    :: PVMOD              ! wind near the road
REAL, DIMENSION(:)  , INTENT(IN)    :: PALBNIR_TVEG       ! nearIR  veg tot albedo
REAL, DIMENSION(:)  , INTENT(IN)    :: PALBVIS_TVEG       ! visible veg tot albedo
REAL, DIMENSION(:)  , INTENT(IN)    :: PALBNIR_TSOIL      ! nearIR  soil tot albedo
REAL, DIMENSION(:)  , INTENT(IN)    :: PALBVIS_TSOIL      ! visible soil tot albedo
!
!
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRN         ! net radiation over green areas
REAL, DIMENSION(:)  , INTENT(OUT)   :: PH          ! sensible heat flux over green areas
REAL, DIMENSION(:)  , INTENT(OUT)   :: PLE         ! latent heat flux over green areas
REAL, DIMENSION(:)  , INTENT(OUT)   :: PGFLUX      ! flux through the green areas
REAL, DIMENSION(:)  , INTENT(OUT)   :: PSFCO2             ! flux of CO2 positive toward the atmosphere (m/s*kg_CO2/kg_air)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PEVAP       ! total evaporation over gardens (kg/m2/s)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PUW         ! friction flux (m2/s2)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PRUNOFF     ! runoff over garden (kg/m2/s)
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC         ! aerodynamical conductance
REAL, DIMENSION(:)  , INTENT(OUT)   :: PQSAT       ! saturation humidity
REAL, DIMENSION(:)  , INTENT(OUT)   :: PAC_AGG     ! aggreg. aeodynamic resistance for green areas for latent heat flux
REAL, DIMENSION(:)  , INTENT(OUT)   :: PHU_AGG     ! aggreg. relative humidity for green areas for latent heat flux
REAL, DIMENSION(:)  , INTENT(OUT)   :: PDRAIN      ! garden total (vertical) drainage
REAL, DIMENSION(:)  , INTENT(OUT)   :: PIRRIG      ! garden summer irrigation rate
!
!
!*      0.2    Declarations of local variables
!
TYPE(SSO_t) :: YSS
TYPE(AGRI_t) :: YAG
!
REAL, DIMENSION(SIZE(PPS)) :: ZDIRCOSZW           ! orography slope cosine (=1 in TEB)
REAL, DIMENSION(SIZE(PPS),GDO%NNBIOMASS) :: ZRESP_BIOMASS_INST       ! instantaneous biomass respiration (kgCO2/kgair m/s)
REAL, DIMENSION(SIZE(PPS)) :: ZUSTAR
!
!  temperatures
!
REAL, DIMENSION(SIZE(PPS)) :: ZTA ! estimate of air temperature at future time
!                                 ! step as if modified by ISBA flux alone.
REAL, DIMENSION(SIZE(PPS)) :: ZDEEP_FLUX ! heat flux at base of the deep soil
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
INTEGER                    :: ILU
LOGICAL :: GMASK, GALB
LOGICAL :: GUPDATED
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*      1.     various initialisations
!              -----------------------
!
IF (LHOOK) CALL DR_HOOK('GARDEN',0,ZHOOK_HANDLE)
ILU = SIZE(PPS)
!
ZDIRCOSZW = 1.
!
CALL SSO_INIT(YSS)
!
CALL AGRI_INIT(YAG)
!
!-------------------------------------------------------------------------------
!
!*      2.     Treatment of green areas
!              ------------------------
!*      2.1    Automatic irrigation
!              --------------------
!
CALL TEB_IRRIG(TIR%LPAR_GD_IRRIG, PTSTEP, TPTIME%TDATE%MONTH, PTSUN,        &
               TIR%XGD_START_MONTH, TIR%XGD_END_MONTH, TIR%XGD_START_HOUR,  &
               TIR%XGD_END_HOUR, TIR%XGD_24H_IRRIG, PIRRIG           ) 
!
! --------------------------------------------------------------------------------------
! Vegetation update (in case of non-interactive vegetation):
! --------------------------------------------------------------------------------------
!
S%TTIME = TPTIME
!
GUPDATED=.FALSE.
GALB = .FALSE. 
IF (GDO%CPHOTO=='LAI'.OR.GDO%CPHOTO=='LST'.OR.GDO%CPHOTO=='NIT'.OR.GDO%CPHOTO=='NCB') GALB = .TRUE.
!
  CALL VEGETATION_UPDATE(DTCO, DTI, G%NDIM, GDO, K, P, PEK, 1, &
                         PTSTEP, S%TTIME, TOP%XCOVER, TOP%LCOVER,  .FALSE., &
                         'GRD', GALB, YSS, GUPDATED, OABSENT=(T%XGARDEN==0.)  )
!
!
VD%DP%XZ0 (:) = PEK%XZ0(:)
VD%DP%XZ0H(:) = PEK%XZ0(:) / P%XZ0_O_Z0H(:)
!
VD%DP%XZ0EFF(:) =  PEK%XZ0(:)
!
!*      2.2    Call ISBA for green areas
!              -------------------------
!
ALLOCATE(GB%XIACAN(SIZE(PPS),SIZE(S%XABC)))
!
 CALL ISBA(GDO, K, P, PEK, G, YAG, VD%DP, VD%DEP, VD%DM, &
           TPTIME, S%XPOI, S%XABC, GB%XIACAN, .FALSE., PTSTEP,             &
           HIMPLICIT_WIND, PZREF, PZREF, ZDIRCOSZW, PTA, PQA, PEXNS, PRHOA, PPS,  &
           PEXNS, PRR, PSR, PZENITH, ZP_MEB_SCA_SW, PSW, PLW, PVMOD, PPEW_A_COEF, &
           PPEW_B_COEF, PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,       &
           PALBNIR_TVEG, PALBVIS_TVEG, PALBNIR_TSOIL, PALBVIS_TSOIL, ZPALPHAN,    &
           ZZ0G_WITHOUT_SNOW, ZZ0_MEBV, ZZ0H_MEBV, ZZ0EFF_MEBV, ZZ0_MEBN,         &
           ZZ0H_MEBN, ZZ0EFF_MEBN, ZTDEEP_A, PCO2, K%XFFG(:), K%XFFV(:),&
           ZEMISF, ZUSTAR, PAC_AGG, PHU_AGG, ZRESP_BIOMASS_INST, ZDEEP_FLUX, PIRRIG )     
!
PRUNOFF(:) = VD%DEP%XRUNOFF(:)
PDRAIN (:) = VD%DEP%XDRAIN(:)
!
IF (PEK%TSNOW%SCHEME=='3-L' .OR. PEK%TSNOW%SCHEME=='CRO') &
    PEK%TSNOW%TS(:)= VD%DM%XSNOWTEMP(:,1)
!
IF (GDO%LTR_ML) THEN
  GMASK = ( TPTIME%TIME - PTSTEP < 0. ) .AND. ( TPTIME%TIME >= 0. )
  IF (GMASK) THEN
    ALLOCATE(VD%DM%XDFAPARC(ILU),VD%DM%XDFAPIRC(ILU),VD%DM%XDLAI_EFFC(ILU))
    VD%DM%XDFAPARC  (:) = PEK%XFAPARC   (:) / PEK%XMUS (:)
    VD%DM%XDFAPIRC  (:) = PEK%XFAPIRC   (:) / PEK%XMUS (:)
    VD%DM%XDLAI_EFFC(:) = PEK%XLAI_EFFC (:) / PEK%XMUS (:)
  ENDIF
ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Vegetation evolution for interactive LAI
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
IF (GDO%CPHOTO=='LAI' .OR. GDO%CPHOTO=='LST' .OR. GDO%CPHOTO=='NIT') THEN
  CALL VEGETATION_EVOL(GDO, P, PEK, .FALSE., PTSTEP, &
                       TPTIME%TDATE%MONTH, TPTIME%TDATE%DAY,     &
                       TPTIME%TIME, G%XLAT, PRHOA, PCO2, YSS, ZRESP_BIOMASS_INST )         
END IF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Diagnostic of respiration carbon fluxes and soil carbon evolution
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
PSFCO2(:)=0.
VD%DEP%XRESP_ECO (:)=0.
VD%DEP%XRESP_AUTO(:)=0.
!
IF (GDO%CPHOTO/='NON' .AND. GDO%CRESPSL/='NON' .AND. ANY(PEK%XLAI(:)/=XUNDEF)) THEN
  CALL CARBON_EVOL(GDO, K, P, PEK, VD%DEP, &
                   PTSTEP, PRHOA, ZRESP_BIOMASS_INST  )
  ! calculation of vegetation CO2 flux
  PSFCO2(:) = - VD%DEP%XGPP(:) + VD%DEP%XRESP_ECO(:)
END IF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!*      4.     Set undefined values for points where there is no garden
!              --------------------------------------------------------
!
! This way, these points are clearly flaged, and one will not try to interpret
! the values for those points
!
 CALL FLAG_TEB_VEG_n(PEK, GDO, T%XGARDEN, 2)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!*      9.     Fields required for TEB
!              -----------------------
!
WHERE (T%XGARDEN/=0.)
  !
  ! energy balance
  !
  PRN    (:) = VD%DP%XRN    (:)
  PH     (:) = VD%DP%XH     (:)
  PLE    (:) = PEK%XLE(:)
  PGFLUX (:) = VD%DP%XGFLUX (:)
  PEVAP  (:) = VD%DP%XEVAP  (:)
  !
  !
  ! Estimate of green area aerodynamic conductance recomputed from heat flux,
  ! surface (radiative) temp. and forcing air temperature (estimated at future time step)
  ZTA = PPET_B_COEF + PPET_A_COEF * PH
  PAC = 0.
  WHERE (VD%DP%XTSRAD /= ZTA)
    PAC(:)   = MAX(PH(:) / XCPD / PRHOA(:) / (VD%DP%XTSRAD - ZTA) , 0.)
  ENDWHERE
  !
  ! Humidity of saturation for green areas
  PQSAT(:) = QSAT(PEK%XTG(:,1),PPS(:))
  !
  !* friction flux
  PUW(:)    = -ZUSTAR(:)**2
  !
ELSEWHERE
  !
  PRN    (:) = XUNDEF
  PH     (:) = XUNDEF
  PLE    (:) = XUNDEF
  PGFLUX (:) = XUNDEF
  PEVAP  (:) = XUNDEF
  PAC    (:) = XUNDEF
  PQSAT  (:) = XUNDEF
  PUW    (:) = XUNDEF
  !
END WHERE
!
IF (LHOOK) CALL DR_HOOK('GARDEN',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE GARDEN
