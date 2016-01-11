!     ###############################################################################
SUBROUTINE COUPLING_ISBA_n (DTCO, UG, U, USS, IM, DTGD, DTGR, DST, SLT,   &
                             HPROGRAM, HCOUPLING,                                              &
                 PTSTEP, KYEAR, KMONTH, KDAY, PTIME, KI, KSV, KSW, PTSUN, PZENITH, PZENITH2, &
                 PZREF, PUREF, PZS, PU, PV, PQA, PTA, PRHOA, PSV, PCO2, HSV,                 &
                 PRAIN, PSNOW, PLW, PDIR_SW, PSCA_SW, PSW_BANDS, PPS, PPA,                   &
                 PSFTQ, PSFTH, PSFTS, PSFCO2, PSFU, PSFV,                                    &
                 PTRAD, PDIR_ALB, PSCA_ALB, PEMIS, PTSURF, PZ0, PZ0H, PQSURF,                &
                 PPEW_A_COEF, PPEW_B_COEF,                                                   &
                 PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,                         &
                 HTEST                                                                       )  
!     ###############################################################################
!
!!****  *COUPLING_ISBA_n * - Driver for ISBA time step   
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!! First, all actions dependant on each patch is donbe independantly
!!     (loop on patches)
!! Second, actions common to all patches (e.g. prescription of new vegetation)
!! Third, energy fluxes are averaged
!!
!! Nota that chemical fluxes are also treated.
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
!!      P Le Moigne 11/2004 add new diagnostics for isba
!!      A.Bogatchev 09/2005 EBA snow option
!!      P Le Moigne 09/2005 AGS modifs of L. Jarlan
!!      P Le Moigne 02/2006 z0h with snow
!!      P.Le Moigne 06/2006 seeding and irrigation
!!      B. Decharme   2008  reset the subgrid topographic effect on the forcing
!!                          PSNV allways <= PSNG
!!                          News diag
!!                          Flooding scheme and allows TRIP variables coupling
!!      A.L. Gibelin 04/2009 : Add respiration diagnostics
!!      A.L. Gibelin 04/2009 : BIOMASS and RESP_BIOMASS arrays 
!!      A.L. Gibelin 04/2009 : TAU_WOOD for NCB option 
!!      A.L. Gibelin 05/2009 : Add carbon spinup
!!      A.L. Gibelin 06/2009 : Soil carbon variables for CNT option
!!      A.L. Gibelin 07/2009 : Suppress RDK and transform GPP as a diagnostic
!!      A.L. Gibelin 07/2009 : Suppress PPST and PPSTF as outputs
!!        S.Lafont   01/2011 : add PTSTEP as arg of diag_misc
!!       B.Decharme  09/2012 : Bug in hydro_glacier calculation with ES or Crocus
!!                             New wind implicitation
!!                             New soil carbon spinup and diag
!!                             Isba budget
!!      F. Bouttier  01/2013 : Apply random perturbations for ensembles
!!      B. Decharme  04/2013 new coupling variables
!!                           Subsurface runoff if SGH (DIF option only)
!!                   07/2013 Surface / Water table depth coupling
!!      P Samuelsson 10/2014 : MEB
!!      P. LeMoigne  12/2014 EBA scheme update
!!      R. Seferian  05/2015 : Add coupling fiels to vegetation_evol call
!!-------------------------------------------------------------------
!
USE MODD_SURFEX_n, ONLY : ISBA_MODEL_t
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SURF_ATM_SSO_n, ONLY : SURF_ATM_SSO_t
USE MODD_DATA_TEB_GARDEN_n, ONLY : DATA_TEB_GARDEN_t
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_DST_n, ONLY : DST_t
USE MODD_SLT_n, ONLY : SLT_t
!
USE MODD_REPROD_OPER, ONLY : CIMPLICIT_WIND
!
USE MODD_CSTS,         ONLY : XRD, XRV, XP00, XCPD, XPI, XAVOGADRO, XMD
USE MODD_CO2V_PAR,     ONLY : XMCO2, XSPIN_CO2
!
USE MODD_SURF_PAR,     ONLY : XUNDEF
USE MODD_SNOW_PAR,     ONLY : XZ0SN
USE MODD_TYPE_DATE_SURF
!
USE MODD_SURF_ATM,    ONLY : LNOSOF
!
USE MODD_DST_SURF
USE MODD_SLT_SURF
USE MODE_DSLT_SURF
USE MODE_MEB
!
!
!                         

USE MODD_DATA_COVER_PAR, ONLY : NVT_NO, NVT_ROCK
!
USE MODD_AGRI,           ONLY : LAGRIP
USE MODD_DEEPSOIL,       ONLY : LDEEPSOIL
!
#ifdef TOPD
USE MODD_COUPLING_TOPD,  ONLY : LCOUPL_TOPD, NMASKT_PATCH
#endif
!
USE MODI_IRRIGATION_UPDATE
USE MODI_ADD_FORECAST_TO_DATE_SURF
USE MODI_Z0EFF
USE MODI_ISBA
USE MODI_AVERAGE_FLUX
USE MODI_AVERAGE_PHY
USE MODI_AVERAGE_RAD
USE MODI_AVERAGE_DIAG_ISBA_n
USE MODI_VEGETATION_EVOL
USE MODI_VEGETATION_UPDATE
USE MODI_ALBEDO_VEG_UPDATE
USE MODI_CARBON_EVOL
USE MODI_SUBSCALE_Z0EFF
USE MODI_SOIL_ALBEDO
USE MODI_ALBEDO
USE MODI_DIAG_INLINE_ISBA_n
USE MODI_DIAG_EVAP_ISBA_n
USE MODI_DIAG_MISC_ISBA_n
!
USE MODI_UPDATE_RAD_ISBA_n
USE MODI_DEEPSOIL_UPDATE
USE MODI_ISBA_SGH_UPDATE
USE MODI_ISBA_FLOOD_PROPERTIES
USE MODI_DIAG_CPL_ESM_ISBA
USE MODI_HYDRO_GLACIER
USE MODI_ISBA_ALBEDO
USE MODI_CARBON_SPINUP
USE MODI_PACK_ISBA_PATCH_n    
USE MODI_PACK_ISBA_PATCH_GET_SIZE_n
USE MODI_PACK_CH_ISBA_PATCH_n     
USE MODI_PACK_DIAG_PATCH_n
USE MODI_PACK_DIAG_PATCH_GET_SIZE_n
USE MODI_UNPACK_ISBA_PATCH_n     
USE MODI_UNPACK_CH_ISBA_PATCH_n     
USE MODI_UNPACK_DIAG_PATCH_n     
USE MODI_CH_AER_DEP
USE MODI_ABOR1_SFX
USE MODI_AVERAGE_DIAG_EVAP_ISBA_n
USE MODI_AVERAGE_DIAG_MISC_ISBA_n
USE MODI_CH_BVOCEM_n
USE MODI_SOILEMISNO_n
USE MODI_CH_DEP_ISBA
USE MODI_DSLT_DEP
USE MODI_COUPLING_DST_n
USE MODI_COUPLING_SURF_TOPD
USE MODI_ISBA_BUDGET_INIT
USE MODI_ISBA_BUDGET
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(ISBA_MODEL_t), INTENT(INOUT) :: IM
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SURF_ATM_SSO_t), INTENT(INOUT) :: USS
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTGD
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTGR
TYPE(DST_t), INTENT(INOUT) :: DST
TYPE(SLT_t), INTENT(INOUT) :: SLT
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
!   
 CHARACTER(LEN=6), DIMENSION(KSV),INTENT(IN):: HSV  ! name of all scalar variables!
REAL, DIMENSION(KI), INTENT(IN)  :: PU        ! zonal wind                            (m/s)
REAL, DIMENSION(KI), INTENT(IN)  :: PV        ! meridian wind                         (m/s)
REAL, DIMENSION(KI,KSW),INTENT(IN) :: PDIR_SW ! direct  solar radiation (on horizontal surf.)
!                                             !                                       (W/m2)
REAL, DIMENSION(KI,KSW),INTENT(IN) :: PSCA_SW ! diffuse solar radiation (on horizontal surf.)
!                                             !                                       (W/m2)
REAL, DIMENSION(KSW),INTENT(IN)  :: PSW_BANDS ! mean wavelength of each shortwave band (m)
REAL, DIMENSION(KI), INTENT(IN)  :: PZENITH   ! zenithal angle at t  (radian from the vertical)
REAL, DIMENSION(KI), INTENT(IN)  :: PZENITH2  ! zenithal angle at t+1(radian from the vertical)
REAL, DIMENSION(KI), INTENT(IN)  :: PLW       ! longwave radiation (on horizontal surf.)
!                                             !                                       (W/m2)
REAL, DIMENSION(KI), INTENT(IN)  :: PPS       ! pressure at atmospheric model surface (Pa)
REAL, DIMENSION(KI), INTENT(IN)  :: PPA       ! pressure at forcing level             (Pa)
REAL, DIMENSION(KI), INTENT(IN)  :: PZS       ! atmospheric model orography           (m)
REAL, DIMENSION(KI), INTENT(IN)  :: PCO2      ! CO2 concentration in the air          (kg_CO2/m3)
REAL, DIMENSION(KI), INTENT(IN)  :: PSNOW     ! snow precipitation                    (kg/m2/s)
REAL, DIMENSION(KI), INTENT(IN)  :: PRAIN     ! liquid precipitation                  (kg/m2/s)
!
!
REAL, DIMENSION(KI), INTENT(OUT) :: PSFTH     ! flux of heat                          (W/m2)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFTQ     ! flux of water vapor                   (kg/m2/s)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFU      ! zonal momentum flux                   (Pa)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFV      ! meridian momentum flux                (Pa)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFCO2    ! flux of CO2 positive toward the atmosphere (m/s*kg_CO2/kg_air)
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
REAL, DIMENSION(KI), INTENT(IN) :: PPEW_A_COEF! implicit coefficients
REAL, DIMENSION(KI), INTENT(IN) :: PPEW_B_COEF! needed if HCOUPLING='I'
REAL, DIMENSION(KI), INTENT(IN) :: PPET_A_COEF
REAL, DIMENSION(KI), INTENT(IN) :: PPEQ_A_COEF
REAL, DIMENSION(KI), INTENT(IN) :: PPET_B_COEF
REAL, DIMENSION(KI), INTENT(IN) :: PPEQ_B_COEF
 CHARACTER(LEN=2),    INTENT(IN) :: HTEST ! must be equal to 'OK'
!
!
!*      0.2    declarations of local variables
!
!* forcing variables
!
REAL, DIMENSION(KI)     :: ZWIND    ! lowest atmospheric level wind speed           (m/s)
REAL, DIMENSION(KI)     :: ZDIR     ! wind direction                        (rad from N clockwise)
REAL, DIMENSION(KI)     :: ZEXNA    ! Exner function at lowest atmospheric level    (-)
REAL, DIMENSION(KI)     :: ZEXNS    ! Exner function at surface                     (-)
REAL, DIMENSION(KI)     :: ZALFA    ! Wind direction                                (-)
REAL, DIMENSION(KI)     :: ZQA      ! specific humidity                             (kg/kg)
REAL, DIMENSION(KI)     :: ZCO2     ! CO2 concentration                             (kg/kg)
REAL                    :: ZSPINCO2 ! CO2 concentration                             (ppmv)
REAL, DIMENSION(KI)     :: ZPEQ_A_COEF ! specific humidity implicit
REAL, DIMENSION(KI)     :: ZPEQ_B_COEF ! coefficients (hum. in kg/kg)
!
INTEGER                 ::ISPINEND
!
! Patch outputs:
!
REAL, DIMENSION(KI,IM%I%O%NPATCH) :: ZSFTH_TILE     ! surface heat flux (W/m2)
REAL, DIMENSION(KI,IM%I%O%NPATCH) :: ZSFTQ_TILE     ! surface vapor flux (kg/m2/s)
REAL, DIMENSION(KI,IM%I%O%NPATCH) :: ZSFCO2_TILE    ! surface CO2 flux positive toward the atmosphere (m/s*kg_CO2/kg_air)
REAL, DIMENSION(KI,IM%I%O%NPATCH) :: ZSFU_TILE      ! zonal momentum flux
REAL, DIMENSION(KI,IM%I%O%NPATCH) :: ZSFV_TILE      ! meridian momentum flux
REAL, DIMENSION(KI,IM%I%O%NPATCH) :: ZTRAD_TILE     ! radiative surface temperature
REAL, DIMENSION(KI,IM%I%O%NPATCH) :: ZEMIS_TILE     ! emissivity
REAL, DIMENSION(KI,IM%I%O%NPATCH) :: ZTSURF_TILE    ! surface effective temperature
REAL, DIMENSION(KI,IM%I%O%NPATCH) :: ZZ0_TILE       ! roughness length for momentum
REAL, DIMENSION(KI,IM%I%O%NPATCH) :: ZZ0H_TILE      ! roughness length for heat
REAL, DIMENSION(KI,IM%I%O%NPATCH) :: ZQSURF_TILE    ! specific humidity at surface
REAL, DIMENSION(KI,KSW,IM%I%O%NPATCH) :: ZDIR_ALB_TILE  ! direct albedo
REAL, DIMENSION(KI,KSW,IM%I%O%NPATCH) :: ZSCA_ALB_TILE  ! diffuse albedo
REAL, DIMENSION(KI,KSV,IM%I%O%NPATCH) :: ZSFTS_TILE     ! scalar surface flux
!
REAL, DIMENSION(KI, IM%I%O%NPATCH) :: ZCPL_DRAIN     ! For the coupling with TRIP
REAL, DIMENSION(KI, IM%I%O%NPATCH) :: ZCPL_RUNOFF    ! For the coupling with TRIP
REAL, DIMENSION(KI, IM%I%O%NPATCH) :: ZCPL_EFLOOD    ! For the coupling with TRIP
REAL, DIMENSION(KI, IM%I%O%NPATCH) :: ZCPL_PFLOOD    ! For the coupling with TRIP
REAL, DIMENSION(KI, IM%I%O%NPATCH) :: ZCPL_IFLOOD    ! For the coupling with TRIP
REAL, DIMENSION(KI, IM%I%O%NPATCH) :: ZCPL_ICEFLUX
!
! for chemical computations
!
REAL, DIMENSION(KI, IM%I%O%NPATCH) :: ZSW_FORBIO
!
REAL                       :: ZCONVERTFACM0_SLT, ZCONVERTFACM0_DST
REAL                       :: ZCONVERTFACM3_SLT, ZCONVERTFACM3_DST
REAL                       :: ZCONVERTFACM6_SLT, ZCONVERTFACM6_DST
!
! dimensions and loop counters
!
INTEGER :: ISWB   ! number of spectral shortwave bands
INTEGER :: JSWB   ! loop on number of spectral shortwave bands
INTEGER :: JPATCH ! loop on patches
INTEGER :: JSV, IDST, IMOMENT, II
INTEGER :: JLAYER, JMODE, JSV_IDX
!
! logical units
!
INTEGER :: JJ
LOGICAL :: LUPDATED              ! T if VEGETATION_UPDATE has reset fields
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! --------------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('COUPLING_ISBA_N',0,ZHOOK_HANDLE)
IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('COUPLING_ISBAN: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF
! --------------------------------------------------------------------------------------
!
!*      1.     Initializations
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Allocations:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
ZSFTH_TILE   (:,:)   = XUNDEF
ZSFTQ_TILE   (:,:)   = XUNDEF
ZSFCO2_TILE  (:,:)   = XUNDEF
ZSFU_TILE    (:,:)   = XUNDEF
ZSFV_TILE    (:,:)   = XUNDEF
ZTRAD_TILE   (:,:)   = XUNDEF
ZEMIS_TILE   (:,:)   = XUNDEF
ZDIR_ALB_TILE(:,:,:) = XUNDEF
ZSCA_ALB_TILE(:,:,:) = XUNDEF
ZTSURF_TILE  (:,:)   = XUNDEF
ZZ0_TILE     (:,:)   = XUNDEF
ZZ0H_TILE    (:,:)   = XUNDEF
ZQSURF_TILE  (:,:)   = XUNDEF
!
ZSFTS_TILE(:,:,:) = 0.
!
ZCPL_DRAIN(:,:)   = 0.0
ZCPL_RUNOFF(:,:)  = 0.0
ZCPL_EFLOOD(:,:)  = 0.0
ZCPL_PFLOOD(:,:)  = 0.0
ZCPL_IFLOOD(:,:)  = 0.0
ZCPL_ICEFLUX(:,:) = 0.0
!
ZSW_FORBIO(:,:)   =  XUNDEF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Forcing Modifications:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
ZDIR=0.
!
DO JJ=1,SIZE(PQA) 
! specific humidity (conversion from kg/m3 to kg/kg)
!
  ZQA(JJ) = PQA(JJ) / PRHOA(JJ)
  ZPEQ_A_COEF(JJ) = PPEQ_A_COEF(JJ) / PRHOA(JJ)
  ZPEQ_B_COEF(JJ) = PPEQ_B_COEF(JJ) / PRHOA(JJ)
!
  ZCO2(JJ) = PCO2(JJ) / PRHOA(JJ)
!
! Other forcing variables depending on incoming forcing (argument list)JJ
!
  ZEXNS(JJ)   = (PPS(JJ)/XP00)**(XRD/XCPD)
  ZEXNA(JJ)   = (PPA(JJ)/XP00)**(XRD/XCPD)
!
!* wind strength
!
  ZWIND(JJ) = SQRT(PU(JJ)**2+PV(JJ)**2)
!
!* wind direction
!
  IF (ZWIND(JJ)>0.)  ZDIR(JJ)=ATAN2(PU(JJ),PV(JJ))
!
!* angle between z0eff J axis and wind direction (rad., clockwise)
!
  ZALFA(JJ) = ZDIR(JJ) - IM%I%P%XZ0EFFJPDIR(JJ) * XPI/180.

  IF (ZALFA(JJ)<-XPI) ZALFA(JJ) = ZALFA(JJ) + 2.*XPI
  IF (ZALFA(JJ)>=XPI) ZALFA(JJ) = ZALFA(JJ) - 2.*XPI
!
ENDDO
!
!* number of shortwave spectral bands
!
ISWB = KSW
!
!* irrigation
!
IF (LAGRIP .AND. (IM%I%O%CPHOTO=='LAI' .OR. IM%I%O%CPHOTO=='LST' .OR. &
                IM%I%O%CPHOTO=='NIT'.OR. IM%I%O%CPHOTO=='NCB') ) THEN
   CALL IRRIGATION_UPDATE(IM%AG, &
                          IM%I%M%I%XIRRIG,PTSTEP,KMONTH,KDAY,PTIME,               &
                            IM%I%M%I%TSEED(:,:)%TDATE%MONTH,IM%I%M%I%TSEED(:,:)%TDATE%DAY,   &
                            IM%I%M%I%TREAP(:,:)%TDATE%MONTH,IM%I%M%I%TREAP(:,:)%TDATE%DAY    )  
ENDIF
!
!* Actualization of the SGH variable (Fmu, Fsat)
!
 CALL ISBA_SGH_UPDATE(IM%IG, IM%I, &
                      IM%I%O%CISBA,IM%I%O%CRUNOFF,IM%I%O%CRAIN,PRAIN,IM%I%I%XMUF,IM%I%I%XFSAT,IM%I%I%XTOPQS)
!
!
!* Actualization of deep soil characteristics
!
IF (LDEEPSOIL) THEN
   CALL DEEPSOIL_UPDATE(IM%I, &
                        IM%I%I%TTIME%TDATE%MONTH)
ENDIF
!
!* Actualization of soil and wood carbon spinup
!
! During soil carbon spinup with ISBA-CC: 
!        (1) Atmospheric CO2 concentration fixed to Pre-industrial CO2 consentration XCO2_START
!        (2) Atmospheric CO2 concentration rampin up from XCO2_START to XCO2_END
!
IF(IM%I%O%LSPINUPCARBS.OR.IM%I%O%LSPINUPCARBW)THEN
!
  ISPINEND=IM%I%O%NNBYEARSPINS-NINT(IM%I%O%NNBYEARSPINS*XSPIN_CO2)
!  
  IM%I%O%LAGRI_TO_GRASS = .FALSE.
!
  IF ( IM%I%O%LSPINUPCARBS .AND. (IM%I%O%NNBYEARSOLD <= ISPINEND) ) THEN
!
   IM%I%O%LAGRI_TO_GRASS = .TRUE.
!
   ZCO2(:) = IM%I%O%XCO2_START * 1.E-6 * XMCO2 / XMD
!
  ELSEIF(IM%I%O%LSPINUPCARBS .AND. (IM%I%O%NNBYEARSOLD > ISPINEND) .AND. &
                  (IM%I%O%NNBYEARSOLD <= IM%I%O%NNBYEARSPINS) )THEN
!
   ZSPINCO2 = IM%I%O%XCO2_START + (IM%I%O%XCO2_END-IM%I%O%XCO2_START) * &
                REAL(IM%I%O%NNBYEARSOLD - ISPINEND) / REAL(IM%I%O%NNBYEARSPINS - ISPINEND)
!
   ZCO2 (:) = ZSPINCO2 * 1.E-6 * XMCO2 / XMD
!
  ENDIF
!
  CALL CARBON_SPINUP(IM%I%I%TTIME%TDATE%MONTH,IM%I%I%TTIME%TDATE%DAY,IM%I%I%TTIME%TIME,       &
                     IM%I%O%LSPINUPCARBS, IM%I%O%LSPINUPCARBW, IM%I%O%XSPINMAXS, IM%I%O%XSPINMAXW,   &
                     IM%I%O%NNBYEARSPINS, IM%I%O%NNBYEARSPINW, IM%I%O%NNBYEARSOLD, IM%I%O%CPHOTO,    &
                     IM%I%O%CRESPSL, IM%I%O%NSPINS, IM%I%O%NSPINW                             )
!
ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Time evolution
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
IM%I%I%TTIME%TIME = IM%I%I%TTIME%TIME + PTSTEP
 CALL ADD_FORECAST_TO_DATE_SURF(IM%I%I%TTIME%TDATE%YEAR,IM%I%I%TTIME%TDATE%MONTH,&
                IM%I%I%TTIME%TDATE%DAY,IM%I%I%TTIME%TIME)
!
! --------------------------------------------------------------------------------------
!
!*      2.     Physical evolution
!
! --------------------------------------------------------------------------------------
! Patch Dependent Calculations
! --------------------------------------------------------------------------------------
!
PATCH_LOOP: DO JPATCH=1,IM%I%O%NPATCH
!
  IF (IM%I%IP%NSIZE_NATURE_P(JPATCH) == 0 ) CYCLE
!
! Pack dummy arguments for each patch:
!
#ifdef TOPD
  IF (LCOUPL_TOPD)NMASKT_PATCH(:)=IM%I%IP%NR_NATURE_P(:,JPATCH)
#endif
  CALL TREAT_PATCH(IM%I%IP%NSIZE_NATURE_P(JPATCH),IM%I%IP%NR_NATURE_P(:,JPATCH))
!
ENDDO PATCH_LOOP
!
! --------------------------------------------------------------------------------------
! SFX - RRM coupling update if used :
! --------------------------------------------------------------------------------------
!
IF(IM%I%O%LCPL_RRM)THEN
  CALL DIAG_CPL_ESM_ISBA(IM%I, &
                         PTSTEP,ZCPL_DRAIN,ZCPL_RUNOFF,ZCPL_EFLOOD, &
                           ZCPL_PFLOOD,ZCPL_IFLOOD,ZCPL_ICEFLUX     )  
ENDIF
!
! --------------------------------------------------------------------------------------
! Vegetation update (in case of non-interactive vegetation):
! Or
! Vegetation albedo only update (in case of interactive vegetation):
! --------------------------------------------------------------------------------------
!
LUPDATED=.FALSE.
IF ((IM%I%O%CPHOTO=='NON' .OR. IM%I%O%CPHOTO=='AGS' .OR. IM%I%O%CPHOTO=='AST') .AND. IM%I%O%LVEGUPD) THEN
     CALL VEGETATION_UPDATE(DTCO, IM%DTI, IM%IG, IM%I%O, &
                            PTSTEP,IM%I%I%TTIME,IM%I%P%XCOVER, IM%I%P%LCOVER,                 &
                         IM%I%O%CISBA,IM%I%O%LECOCLIMAP,IM%I%O%CPHOTO,LAGRIP,IM%I%O%LTR_ML,'NAT',    &
                         IM%I%M%T%XLAI,IM%I%M%T%XVEG,IM%I%M%T%XZ0,                                  &
                         IM%I%M%T%XALBNIR,IM%I%M%T%XALBVIS,IM%I%M%T%XALBUV,IM%I%M%T%XEMIS,                   &
                         IM%I%M%T%XRSMIN,IM%I%M%T%XGAMMA,IM%I%M%T%XWRMAX_CF,                        &
                         IM%I%M%T%XRGL,IM%I%M%T%XCV,                                       &
                         IM%I%M%T%XGMES,IM%I%M%T%XBSLAI,IM%I%M%T%XLAIMIN,IM%I%M%T%XSEFOLD,IM%I%M%T%XGC,  &
                         IM%I%M%T%XF2I, IM%I%M%T%LSTRESS,                                  &
                         IM%I%P%XAOSIP,IM%I%P%XAOSIM,IM%I%P%XAOSJP,IM%I%P%XAOSJM,                    &
                         IM%I%P%XHO2IP,IM%I%P%XHO2IM,IM%I%P%XHO2JP,IM%I%P%XHO2JM,                    &
                         IM%I%IP%XZ0EFFIP,IM%I%IP%XZ0EFFIM,IM%I%IP%XZ0EFFJP,IM%I%IP%XZ0EFFJM,            &
                         IM%I%O%CALBEDO, IM%I%M%T%XALBNIR_VEG, IM%I%M%T%XALBVIS_VEG, IM%I%M%T%XALBUV_VEG,  &
                         IM%I%M%A%XALBNIR_SOIL, IM%I%M%A%XALBVIS_SOIL, IM%I%M%A%XALBUV_SOIL,        &
                         IM%I%M%T%XCE_NITRO, IM%I%M%T%XCF_NITRO, IM%I%M%T%XCNA_NITRO,               &
                         IM%I%M%I%TSEED, IM%I%M%I%TREAP, IM%I%M%I%XWATSUP, IM%I%M%I%XIRRIG,                  &
                         IM%I%M%M%XGNDLITTER, IM%I%M%M%XRGLGV,IM%I%M%M%XGAMMAGV,                     &
                         IM%I%M%M%XRSMINGV, IM%I%M%M%XWRMAX_CFGV,                          &
                         IM%I%M%M%XH_VEG, IM%I%M%M%XLAIGV, IM%I%M%M%XZ0LITTER, LUPDATED             )  
!
ELSEIF ((IM%I%O%CPHOTO=='LAI'.OR.IM%I%O%CPHOTO=='LST'.OR.IM%I%O%CPHOTO=='NIT'.OR.IM%I%O%CPHOTO=='NCB').AND.IM%I%O%LVEGUPD) THEN
!
  CALL ALBEDO_VEG_UPDATE(DTCO, IM%DTI, IM%IG, IM%I%O, &
                         PTSTEP,IM%I%I%TTIME,IM%I%P%XCOVER, IM%I%P%LCOVER,                    &
                         IM%I%O%CISBA,IM%I%O%LECOCLIMAP,IM%I%O%CPHOTO,LAGRIP,IM%I%O%LTR_ML,'NAT',    &
                         IM%I%M%T%XVEG,IM%I%M%T%XALBNIR,IM%I%M%T%XALBVIS,IM%I%M%T%XALBUV,                    &
                         IM%I%O%CALBEDO, IM%I%M%T%XALBNIR_VEG, IM%I%M%T%XALBVIS_VEG, IM%I%M%T%XALBUV_VEG,  &
                         IM%I%M%A%XALBNIR_SOIL, IM%I%M%A%XALBVIS_SOIL, IM%I%M%A%XALBUV_SOIL         )
END IF
!
IF(IM%I%O%LPERTSURF.AND.LUPDATED) THEN
  ! random perturbation for ensembles:
  ! reset these fields to their original values, as in compute_isba_parameters
  IM%I%M%T%XVEG(:,1) = IM%I%I%XPERTVEG(:)
  IM%I%M%T%XLAI(:,1) = IM%I%I%XPERTLAI(:)
  IM%I%M%T%XCV(:,1)  = IM%I%I%XPERTCV(:)
  ! reapply original perturbation patterns
  WHERE(IM%I%M%T%XALBNIR(:,1)/=XUNDEF)  IM%I%M%T%XALBNIR(:,1) =IM%I%M%T%XALBNIR(:,1) *( 1.+ IM%I%I%XPERTALB(:) )
  WHERE(IM%I%M%T%XALBVIS(:,1)/=XUNDEF)  IM%I%M%T%XALBVIS(:,1) =IM%I%M%T%XALBVIS(:,1) *( 1.+ IM%I%I%XPERTALB(:) )
  WHERE(IM%I%M%T%XALBUV(:,1)/=XUNDEF)   IM%I%M%T%XALBUV(:,1)  =IM%I%M%T%XALBUV(:,1)  *( 1.+ IM%I%I%XPERTALB(:) )
  WHERE(IM%I%M%T%XZ0(:,1)/=XUNDEF)      IM%I%M%T%XZ0(:,1)     =IM%I%M%T%XZ0(:,1)     *( 1.+ IM%I%I%XPERTZ0(:) )
  WHERE(IM%I%IP%XZ0EFFIP(:,1)/=XUNDEF) IM%I%IP%XZ0EFFIP(:,1)=IM%I%IP%XZ0EFFIP(:,1)*( 1.+ IM%I%I%XPERTZ0(:) )
  WHERE(IM%I%IP%XZ0EFFIM(:,1)/=XUNDEF) IM%I%IP%XZ0EFFIM(:,1)=IM%I%IP%XZ0EFFIM(:,1)*( 1.+ IM%I%I%XPERTZ0(:) )
  WHERE(IM%I%IP%XZ0EFFJP(:,1)/=XUNDEF) IM%I%IP%XZ0EFFJP(:,1)=IM%I%IP%XZ0EFFJP(:,1)*( 1.+ IM%I%I%XPERTZ0(:) )
  WHERE(IM%I%IP%XZ0EFFJM(:,1)/=XUNDEF) IM%I%IP%XZ0EFFJM(:,1)=IM%I%IP%XZ0EFFJM(:,1)*( 1.+ IM%I%I%XPERTZ0(:) )

ENDIF
!
! --------------------------------------------------------------------------------------
! Outputs for the atmospheric model or update the snow/flood fraction at time t+1
! --------------------------------------------------------------------------------------
! Grid box average fluxes/properties: Arguments and standard diagnostics at time t+1
!
 CALL AVERAGE_FLUX(IM%I%IP%XPATCH,                                             &
                  ZSFTH_TILE, ZSFTQ_TILE, ZSFTS_TILE, ZSFCO2_TILE,    &
                  ZSFU_TILE, ZSFV_TILE,                               &                   
                  PSFTH, PSFTQ, PSFTS, PSFCO2,                        &
                  PSFU, PSFV                                          )  
!
!
!-------------------------------------------------------------------------------
!Physical properties see by the atmosphere in order to close the energy budget 
!between surfex and the atmosphere. All variables should be at t+1 but very 
!difficult to do. Maybe it will be done later. However, Ts is at time t+1
!-------------------------------------------------------------------------------
!   
 CALL AVERAGE_PHY(IM%I%IP%XPATCH,                                         &
                  ZTSURF_TILE, ZZ0_TILE, ZZ0H_TILE, ZQSURF_TILE,  &    
                  PUREF, PZREF, PTSURF, PZ0, PZ0H, PQSURF         )
!
!-------------------------------------------------------------------------------------
!Radiative properties at time t+1 (see by the atmosphere) in order to close
!the energy budget between surfex and the atmosphere
!-------------------------------------------------------------------------------------
!
 CALL UPDATE_RAD_ISBA_n(IM%I, &
                        IM%I%O%LFLOOD, IM%I%R%TSNOW%SCHEME, PZENITH2, PSW_BANDS,      &
                       IM%I%M%T%XVEG, IM%I%M%T%XLAI, IM%I%M%T%XZ0,                                 &
                       IM%I%O%LMEB_PATCH,IM%I%M%M%XLAIGV,IM%I%M%M%XGNDLITTER,IM%I%M%M%XZ0LITTER,IM%I%M%M%XH_VEG,   &
                       IM%I%M%T%XALBNIR, IM%I%M%T%XALBVIS, IM%I%M%T%XALBUV, IM%I%M%T%XEMIS,                 &
                       ZDIR_ALB_TILE,ZSCA_ALB_TILE,ZEMIS_TILE,          &
                       PDIR_SW, PSCA_SW,                                &
                       IM%I%M%T%XALBNIR_VEG, IM%I%M%A%XALBNIR_SOIL,                       &
                       IM%I%M%T%XALBVIS_VEG, IM%I%M%A%XALBVIS_SOIL                        )
!
 CALL AVERAGE_RAD(IM%I%IP%XPATCH,                                              &
                 ZDIR_ALB_TILE, ZSCA_ALB_TILE, ZEMIS_TILE, ZTRAD_TILE, &
                 PDIR_ALB,      PSCA_ALB,      IM%I%I%XEMIS_NAT,  IM%I%R%XTSRAD_NAT  )  
!
PEMIS = IM%I%I%XEMIS_NAT
PTRAD = IM%I%R%XTSRAD_NAT
!
!-------------------------------------------------------------------------------------
!
! Any additional diagnostics (stored in MODD_DIAG_ISBA_n)
!
 CALL AVERAGE_DIAG_ISBA_n(IM%DGI, IM%DGIC, IM%DGIP, IM%DGIPC, &
              IM%DGEI%LSURF_BUDGETC, IM%I%O%LCANOPY, IM%I%IP%XPATCH, IM%I%R%XLE, &
                          PUREF,PZREF,PSFCO2,PTRAD)
!
! Cumulated diagnostics (stored in MODD_DIAG_EVAP_ISBA_n)
!
 CALL AVERAGE_DIAG_EVAP_ISBA_n(IM%DGEI, IM%DGEIC, IM%DGEIP, IM%DGEIPC, &
                               IM%I%O%LGLACIER, IM%I%O%LMEB_PATCH, IM%I%IP%XPATCH, &
                               PTSTEP,PRAIN,PSNOW)
!
! Miscellaneous diagnostics (stored in MODD_DIAG_MISC_ISBA_n)
!
 CALL AVERAGE_DIAG_MISC_ISBA_n(IM%DGMI, IM%DGMIP, IM%I)
!
!--------------------------------------------------------------------------------------
!
 CALL COUPLING_SURF_TOPD(IM%DGEI, IM%DGEIC, IM%DGIC, IM%DGMI, IM%IG, IM%I, UG, U, &
                         HPROGRAM,U%NDIM_FULL)
!
! --------------------------------------------------------------------------------------
! Snow/Flood fractions, albedo and emissivity update :
! --------------------------------------------------------------------------------------
!
! --------------------------------------------------------------------------------------
! Chemical fluxes :
! --------------------------------------------------------------------------------------
!
IF (IM%CHI%SVI%NBEQ>0 .AND. IM%CHI%LCH_BIO_FLUX) THEN
 CALL CH_BVOCEM_n(IM%CHI, IM%GB, IM%I, &
                  ZSW_FORBIO,PRHOA,PSFTS)
ENDIF
!
!SOILNOX
IF (IM%CHI%LCH_NO_FLUX) THEN
  CALL SOILEMISNO_n(IM%GB, IM%I, &
                    PU,PV)
ENDIF
!
!==========================================================================================
!
IF (LHOOK) CALL DR_HOOK('COUPLING_ISBA_N',1,ZHOOK_HANDLE)
CONTAINS
!
!=======================================================================================
SUBROUTINE TREAT_PATCH(KSIZE,KMASK)
!
IMPLICIT NONE
!
INTEGER, INTENT(IN)               :: KSIZE
INTEGER, INTENT(IN), DIMENSION(KI) :: KMASK
!
REAL, DIMENSION(KSIZE) :: ZP_ZREF    ! height of T,q forcing                 (m)
REAL, DIMENSION(KSIZE) :: ZP_UREF    ! height of wind forcing                (m)
REAL, DIMENSION(KSIZE) :: ZP_U       ! zonal wind                            (m/s)
REAL, DIMENSION(KSIZE) :: ZP_V       ! meridian wind                         (m/s)
REAL, DIMENSION(KSIZE) :: ZP_WIND    ! wind                                  (m/s)
REAL, DIMENSION(KSIZE) :: ZP_DIR     ! wind direction                        (rad from N clockwise)
REAL, DIMENSION(KSIZE) :: ZP_QA      ! air specific humidity forcing         (kg/kg)
REAL, DIMENSION(KSIZE) :: ZP_TA      ! air temperature forcing               (K)
REAL, DIMENSION(KSIZE) :: ZP_CO2     ! CO2 concentration in the air          (kg/kg)
REAL, DIMENSION(KSIZE,KSV) :: ZP_SV      ! scalar concentration in the air       (kg/kg)
REAL, DIMENSION(KSIZE) :: ZP_ZENITH  ! zenithal angle        radian from the vertical)
REAL, DIMENSION(KSIZE) :: ZP_PEW_A_COEF ! implicit coefficients
REAL, DIMENSION(KSIZE) :: ZP_PEW_B_COEF ! needed if HCOUPLING='I'
REAL, DIMENSION(KSIZE) :: ZP_PET_A_COEF
REAL, DIMENSION(KSIZE) :: ZP_PET_B_COEF
REAL, DIMENSION(KSIZE) :: ZP_PEQ_A_COEF
REAL, DIMENSION(KSIZE) :: ZP_PEQ_B_COEF
REAL, DIMENSION(KSIZE) :: ZP_RAIN    ! liquid precipitation                  (kg/m2/s)
REAL, DIMENSION(KSIZE) :: ZP_SNOW    ! snow precipitation                    (kg/m2/s)
REAL, DIMENSION(KSIZE) :: ZP_LW      ! longwave radiation (W/m2)
REAL, DIMENSION(KSIZE,ISWB) :: ZP_DIR_SW  ! direct  solar radiation (W/m2)
REAL, DIMENSION(KSIZE,ISWB) :: ZP_SCA_SW  ! diffuse solar radiation (W/m2)
REAL, DIMENSION(KSIZE) :: ZP_PS      ! pressure at atmospheric model surface (Pa)
REAL, DIMENSION(KSIZE) :: ZP_PA      ! pressure at forcing level             (Pa)
REAL, DIMENSION(KSIZE) :: ZP_ZS      ! atmospheric model orography           (m)
REAL, DIMENSION(KSIZE) :: ZP_SFTQ    ! flux of water vapor <w'q'>            (kg.m-2.s-1)
REAL, DIMENSION(KSIZE) :: ZP_SFTH    ! flux of temperature <w'T'>            (W/m2)
REAL, DIMENSION(KSIZE,KSV) :: ZP_SFTS    ! flux of scalar      <w'sv'>           (mkg/kg/s)
REAL, DIMENSION(KSIZE) :: ZP_SFCO2   ! flux of CO2 positive toward the atmosphere (m/s*kg_CO2/kg_air)
REAL, DIMENSION(KSIZE) :: ZP_USTAR   ! friction velocity                     (m/s)
REAL, DIMENSION(KSIZE) :: ZP_SFU     ! zonal momentum flux                   (pa)
REAL, DIMENSION(KSIZE) :: ZP_SFV     ! meridian momentum flux                (pa)
REAL, DIMENSION(KSIZE) :: ZP_TRAD    ! radiative temperature                 (K)
REAL, DIMENSION(KSIZE) :: ZP_TSURF   ! surface effective temperature (K)
REAL, DIMENSION(KSIZE) :: ZP_Z0      ! roughness length for momentum (m)
REAL, DIMENSION(KSIZE) :: ZP_Z0H     ! roughness length for heat     (m)
REAL, DIMENSION(KSIZE) :: ZP_QSURF   ! specific humidity at surface  (kg/kg)
!
!*  other forcing variables (packed for each patch)
!
REAL, DIMENSION(KSIZE) :: ZP_RHOA    ! lowest atmospheric level air density          (kg/m3)
REAL, DIMENSION(KSIZE) :: ZP_EXNA    ! Exner function at lowest atmospheric level    (-)
REAL, DIMENSION(KSIZE) :: ZP_EXNS    ! Exner function at surface                     (-)
REAL, DIMENSION(KSIZE) :: ZP_ALFA    ! Wind direction   (-)
!
!*  working variables (packed for each patch)
!
REAL, DIMENSION(KSIZE)      :: ZP_ALBNIR_TVEG         ! total vegetation albedo in ir
REAL, DIMENSION(KSIZE)      :: ZP_ALBNIR_TSOIL        ! total soil albedo in ir
REAL, DIMENSION(KSIZE)      :: ZP_ALBVIS_TVEG         ! total vegetation albedo in vis
REAL, DIMENSION(KSIZE)      :: ZP_ALBVIS_TSOIL        ! total soil albedo in vis
REAL, DIMENSION(KSIZE) :: ZP_EMIS                      ! emissivity
REAL, DIMENSION(KSIZE) :: ZP_GLOBAL_SW                 ! global incoming SW rad.
REAL, DIMENSION(KSIZE) :: ZP_SLOPE_COS                 ! typical slope in the grid cosine
!
REAL, DIMENSION(KSIZE) :: ZP_Z0FLOOD  !Floodplain 
REAL, DIMENSION(KSIZE) :: ZP_FFGNOS   !Floodplain fraction over the ground without snow
REAL, DIMENSION(KSIZE) :: ZP_FFVNOS   !Floodplain fraction over vegetation without snow
!
REAL, DIMENSION(KSIZE,IM%I%O%NNBIOMASS) :: ZP_RESP_BIOMASS_INST         ! instantaneous biomass respiration (kgCO2/kgair m/s)
!
!*  Aggregated coeffs for evaporative flux calculations
!
REAL, DIMENSION(KSIZE) :: ZP_AC_AGG      ! aggregated aerodynamic resistance
REAL, DIMENSION(KSIZE) :: ZP_HU_AGG      ! aggregated relative humidity
!
!*  For multi-energy balance
!
REAL, DIMENSION(KSIZE) :: ZPALPHAN                     ! snow/canopy transition coefficient
REAL, DIMENSION(KSIZE) :: ZSNOWDEPTH                   ! total snow depth
REAL, DIMENSION(KSIZE) :: ZZ0G_WITHOUT_SNOW            ! roughness length for momentum at snow-free canopy floor
REAL, DIMENSION(KSIZE) :: ZZ0_MEBV                     ! roughness length for momentum over MEB vegetation part of patch
REAL, DIMENSION(KSIZE) :: ZZ0H_MEBV                    ! roughness length for heat over MEB vegetation part of path
REAL, DIMENSION(KSIZE) :: ZZ0EFF_MEBV                  ! roughness length for momentum over MEB vegetation part of patch
REAL, DIMENSION(KSIZE) :: ZZ0_MEBN                     ! roughness length for momentum over MEB snow part of patch
REAL, DIMENSION(KSIZE) :: ZZ0H_MEBN                    ! roughness length for heat over MEB snow part of path
REAL, DIMENSION(KSIZE) :: ZZ0EFF_MEBN                  ! roughness length for momentum over MEB snow part of patch
! Temporary
REAL, DIMENSION(KSIZE) :: ZP_MEB_SCA_SW                ! diffuse incoming SW rad.
!
!*  ISBA water and energy budget
!
REAL, DIMENSION(KSIZE) :: ZP_WG_INI
REAL, DIMENSION(KSIZE) :: ZP_WGI_INI
REAL, DIMENSION(KSIZE) :: ZP_WR_INI
REAL, DIMENSION(KSIZE) :: ZP_SWE_INI
!
! miscellaneous
!
REAL, DIMENSION(KSIZE)               :: ZP_DEEP_FLUX ! Flux at the bottom of the soil
REAL, DIMENSION(KSIZE)               :: ZP_TDEEP_A   ! coefficient for implicitation of Tdeep
REAL, DIMENSION(KSIZE)               :: ZIRRIG_GR    ! green roof ground irrigation rate 
!
! For multi-energy balance
LOGICAL :: GMEB  ! True if multi-energy balance should be used for the specific patch
!
INTEGER :: JJ, JI, JK
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('COUPLING_ISBA_n:TREAT_PATCH',0,ZHOOK_HANDLE)
!
!--------------------------------------------------------------------------------------
!
! Pack isba forcing outputs
!
IF (IM%I%O%NPATCH==1) THEN
   ZP_ZENITH(:)     = PZENITH     (:)
   ZP_ZREF(:)       = PZREF       (:)
   ZP_UREF(:)       = PUREF       (:)
   ZP_WIND(:)       = ZWIND       (:)
   ZP_U(:)          = PU          (:)
   ZP_V(:)          = PV          (:)
   ZP_DIR(:)        = ZDIR        (:)
   ZP_QA(:)         = ZQA         (:)
   ZP_TA(:)         = PTA         (:)
   ZP_CO2(:)        = ZCO2        (:)
   ZP_SV(:,:)       = PSV         (:,:)
   ZP_PEW_A_COEF(:) = PPEW_A_COEF (:)
   ZP_PEW_B_COEF(:) = PPEW_B_COEF (:)
   ZP_PET_A_COEF(:) = PPET_A_COEF (:)
   ZP_PET_B_COEF(:) = PPET_B_COEF (:)
   ZP_PEQ_A_COEF(:) = ZPEQ_A_COEF (:)
   ZP_PEQ_B_COEF(:) = ZPEQ_B_COEF (:)
   ZP_RAIN(:)       = PRAIN       (:)
   ZP_SNOW(:)       = PSNOW       (:)
   ZP_LW(:)         = PLW         (:)
   ZP_DIR_SW(:,:)   = PDIR_SW     (:,:)
   ZP_SCA_SW(:,:)   = PSCA_SW     (:,:)
   ZP_PS(:)         = PPS         (:)
   ZP_PA(:)         = PPA         (:)
   ZP_ZS(:)         = PZS         (:)
!
   ZP_RHOA(:)       = PRHOA       (:)
   ZP_EXNA(:)       = ZEXNA       (:)
   ZP_EXNS(:)       = ZEXNS       (:)
   ZP_ALFA(:)       = ZALFA       (:)
ELSE
!cdir nodep
!cdir unroll=8
  DO JJ=1,KSIZE
   JI = KMASK(JJ)
   ZP_ZENITH(JJ)     = PZENITH     (JI)
   ZP_ZREF(JJ)       = PZREF       (JI)
   ZP_UREF(JJ)       = PUREF       (JI)
   ZP_WIND(JJ)       = ZWIND       (JI)
   ZP_U(JJ)          = PU          (JI)
   ZP_V(JJ)          = PV          (JI)
   ZP_DIR(JJ)        = ZDIR        (JI)
   ZP_QA(JJ)         = ZQA         (JI)
   ZP_TA(JJ)         = PTA         (JI)
   ZP_CO2(JJ)        = ZCO2        (JI)
   ZP_PEW_A_COEF(JJ) = PPEW_A_COEF (JI)
   ZP_PEW_B_COEF(JJ) = PPEW_B_COEF (JI)
   ZP_PET_A_COEF(JJ) = PPET_A_COEF (JI)
   ZP_PET_B_COEF(JJ) = PPET_B_COEF (JI)
   ZP_PEQ_A_COEF(JJ) = ZPEQ_A_COEF (JI)
   ZP_PEQ_B_COEF(JJ) = ZPEQ_B_COEF (JI)
   ZP_RAIN(JJ)       = PRAIN       (JI)
   ZP_SNOW(JJ)       = PSNOW       (JI)
   ZP_LW(JJ)         = PLW         (JI)
   ZP_PS(JJ)         = PPS         (JI)
   ZP_PA(JJ)         = PPA         (JI)
   ZP_ZS(JJ)         = PZS         (JI)
!
   ZP_RHOA(JJ)       = PRHOA       (JI)
   ZP_EXNA(JJ)       = ZEXNA       (JI)
   ZP_EXNS(JJ)       = ZEXNS       (JI)
   ZP_ALFA(JJ)       = ZALFA       (JI)
  ENDDO
!
  DO JK=1,KSV
!cdir nodep
!cdir unroll=8
    DO JJ=1,KSIZE
      JI=KMASK(JJ)
      ZP_SV(JJ,JK) = PSV(JI,JK)
    ENDDO
  ENDDO
!
  DO JK=1,SIZE(PDIR_SW,2)
!cdir nodep
!cdir unroll=8
    DO JJ=1,KSIZE
      JI=KMASK(JJ)
      ZP_DIR_SW(JJ,JK) = PDIR_SW (JI,JK)
      ZP_SCA_SW(JJ,JK) = PSCA_SW (JI,JK)
    ENDDO
  ENDDO
!
ENDIF
!
!--------------------------------------------------------------------------------------
!
! For multi-energy balance
   GMEB=IM%I%O%LMEB_PATCH(JPATCH)
!
! Pack ISBA input and prognostic variables (modd_isban) for each patch:
!
 CALL PACK_ISBA_PATCH_GET_SIZE_n(IM%I, IM%PK, &
                                 JPATCH)
!
 CALL PACK_DIAG_PATCH_GET_SIZE_n(IM%DGEI, IM%DGI, IM%DGMI, IM%I, IM%PKD, &
                                 JPATCH)
!
 CALL PACK_ISBA_PATCH_n(IM%AG, IM%IG, IM%I, IM%PK, &
                        KMASK,KSIZE,JPATCH)     
!
! Pack chemistry input and prognostic variables (modd_ch_isban) for each patch:
!
IF (IM%CHI%SVI%NBEQ>0) THEN
  IF( IM%CHI%CCH_DRY_DEP == "WES89") THEN
    CALL PACK_CH_ISBA_PATCH_n(IM%CHI, IM%PKCI, &
                              KMASK,KSIZE,IM%I%O%NPATCH,JPATCH)     
  END IF
END IF
!
! Allocate ISBA diagnostics for each patch:
!
 CALL PACK_DIAG_PATCH_n(IM%DGEI, IM%DGI, IM%DGMI, IM%I, IM%PK, IM%PKD, &
                        KSIZE,ISWB,JPATCH)     
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Cosine of the slope typically encoutered in the grid mesh (including subgrid orography)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
ZP_SLOPE_COS(:) = 1./SQRT(1.+IM%PK%I%P%XSSO_SLOPE(:)**2)
IF(LNOSOF)ZP_SLOPE_COS(:) = 1.0
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Snow fractions
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! now caculated at the initialization and at the end of the time step 
! (see update_frac_alb_emis_isban.f90) in order to close the energy budget
! between surfex and the atmosphere. This fact do not change the offline runs.
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! No implicitation of Tdeep
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ZP_TDEEP_A = 0.
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Flood properties 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
IF(IM%I%O%LFLOOD)THEN
  CALL ISBA_FLOOD_PROPERTIES(IM%PK%I%M%T%XLAI(:,1),IM%PK%I%I%XFFLOOD,IM%PK%I%I%XFFROZEN(:,1),  &
                             ZP_Z0FLOOD,ZP_FFGNOS,ZP_FFVNOS)  
ELSE
  ZP_Z0FLOOD = XUNDEF
  ZP_FFGNOS  = 0.0
  ZP_FFVNOS  = 0.0
ENDIF
!
! For multi-energy balance
   IF(GMEB)THEN
     ZSNOWDEPTH(:) = SUM(IM%PK%I%R%TSNOW%WSNOW(:,:,1)/IM%PK%I%R%TSNOW%RHO(:,:,1),2)
     ZPALPHAN(:)=MEBPALPHAN(ZSNOWDEPTH,IM%PK%I%M%M%XH_VEG(:,1))
   ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Surface Roughness lengths (m):
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!* effective roughness
!
 CALL Z0EFF(IM%I%R%TSNOW%SCHEME, &
            IM%I%O%CROUGH, GMEB, ZP_ALFA, ZP_ZREF, ZP_UREF, &
            IM%PK%I%M%T%XZ0(:,1), IM%PK%I%IP%XZ0REL, IM%PK%I%R%XPSN(:,1),   &
            ZPALPHAN,IM%PK%I%M%M%XZ0LITTER(:,1), IM%PK%I%R%TSNOW%WSNOW(:,1,1),           &
     IM%PK%I%IP%XZ0EFFIP(:,1),IM%PK%I%IP%XZ0EFFIM(:,1),IM%PK%I%IP%XZ0EFFJP(:,1), &
     IM%PK%I%IP%XZ0EFFJM(:,1), IM%PK%I%I%XFF(:,1), ZP_Z0FLOOD,     &
     IM%PK%I%P%XAOSIP,IM%PK%I%P%XAOSIM,IM%PK%I%P%XAOSJP,IM%PK%I%P%XAOSJM,         &
     IM%PK%I%P%XHO2IP,IM%PK%I%P%XHO2IM,IM%PK%I%P%XHO2JP,IM%PK%I%P%XHO2JM,        &
     IM%PK%I%M%X%XZ0_O_Z0H(:,1), IM%PKD%DGIP%XZ0, IM%PKD%DGIP%XZ0H, &
     IM%PKD%DGIP%XZ0EFF, ZZ0G_WITHOUT_SNOW,                                 &
     ZZ0_MEBV,ZZ0H_MEBV,ZZ0EFF_MEBV,                                     &
     ZZ0_MEBN,ZZ0H_MEBN,ZZ0EFF_MEBN                                      )
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Shortwave computations for outputs (albedo for radiative scheme)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! now caculated at the initialization and at the end of the time step 
! (see update_frac_alb_emis_isban.f90) in order to close the energy budget
! between surfex and the atmosphere. This fact do not change the offline runs.
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Shortwave computations for ISBA inputs (global snow-free albedo)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! ISBA needs global incoming solar radiation: it currently does
! not distinguish between the scattered and direct components,
! or between different wavelengths.
!
!
!* Snow-free surface albedo for each wavelength
!
 CALL ISBA_ALBEDO(IM%I%R%TSNOW%SCHEME, IM%I%O%LTR_ML, GMEB,                &
                   ZP_DIR_SW, ZP_SCA_SW, PSW_BANDS,ISWB,                 &
                   IM%PK%I%M%T%XALBNIR(:,1), IM%PK%I%M%T%XALBVIS(:,1), IM%PK%I%M%T%XALBUV(:,1),           &
                   IM%PK%I%M%T%XALBNIR_VEG(:,1), IM%PK%I%M%T%XALBVIS_VEG(:,1), IM%PK%I%M%T%XALBUV_VEG(:,1),           &
                   IM%PK%I%M%A%XALBNIR_SOIL(:,1), IM%PK%I%M%A%XALBVIS_SOIL(:,1), IM%PK%I%M%A%XALBUV_SOIL(:,1),        &
                   IM%PK%I%I%XALBF(:,1), IM%PK%I%I%XFFV(:,1), IM%PK%I%I%XFFG(:,1),    & 
                   ZP_GLOBAL_SW, IM%PK%I%R%XSNOWFREE_ALB(:,1), IM%PK%I%R%XSNOWFREE_ALB_VEG(:,1),   &
                   IM%PK%I%R%XSNOWFREE_ALB_SOIL(:,1), ZP_MEB_SCA_SW,                  &
                   ZP_ALBNIR_TVEG, ZP_ALBVIS_TVEG,                       &
                   ZP_ALBNIR_TSOIL, ZP_ALBVIS_TSOIL                      )  
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Intialize computation of ISBA water and energy budget
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
 CALL ISBA_BUDGET_INIT(IM%DGEI%LWATER_BUDGET, &
                       IM%I%O%CISBA,IM%I%R%TSNOW%SCHEME,            &
                      IM%PK%I%R%XWG(:,:,1),IM%PK%I%R%XWGI(:,:,1),IM%PK%I%R%XWR(:,1),IM%PK%I%R%TSNOW%WSNOW(:,:,1), &
                      IM%PK%I%M%X%XDG(:,:,1), IM%PK%I%IP%XDZG(:,:,1), ZP_WG_INI,      &
                      ZP_WGI_INI, ZP_WR_INI,         &
                      ZP_SWE_INI                     )
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Over Natural Land Surfaces:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ZIRRIG_GR(:)= 0.
!
 CALL ISBA(IM%I%O%CISBA, IM%I%O%CPHOTO, IM%I%O%LTR_ML, IM%I%O%CRUNOFF, IM%I%O%CKSAT, IM%I%O%CRAIN, &
           IM%I%O%CHORT, IM%I%O%CC1DRY, IM%I%O%CSCOND, IM%I%R%TSNOW%SCHEME, IM%I%O%CSNOWRES, &
           IM%I%O%CCPSURF, IM%I%O%CSOILFRZ, IM%I%O%CDIFSFCOND, IM%I%I%TTIME, IM%I%O%LFLOOD, &
           IM%I%O%LTEMP_ARP, IM%I%O%LGLACIER, GMEB, IM%I%O%LFORC_MEASURE, IM%I%O%LMEB_LITTER, PTSTEP, &
           CIMPLICIT_WIND, IM%I%O%LAGRI_TO_GRASS, IM%I%O%LSNOWDRIFT, IM%I%O%LSNOWDRIFT_SUBLIM, &
           IM%I%O%LSNOW_ABS_ZENITH, IM%I%O%CSNOWMETAMO, IM%I%O%CSNOWRAD, IM%I%O%XCGMAX, ZP_ZREF, &
           ZP_UREF, ZP_SLOPE_COS, ZP_TA, ZP_QA, ZP_EXNA, ZP_RHOA, ZP_PS, ZP_EXNS, ZP_RAIN, &
           ZP_SNOW, ZP_ZENITH, ZP_MEB_SCA_SW, ZP_GLOBAL_SW, ZP_LW, ZP_WIND, ZP_PEW_A_COEF, &
           ZP_PEW_B_COEF, ZP_PET_A_COEF, ZP_PEQ_A_COEF,  ZP_PET_B_COEF, ZP_PEQ_B_COEF, &
           IM%PK%I%M%T%XRSMIN(:,1), IM%PK%I%M%T%XRGL(:,1), IM%PK%I%M%T%XGAMMA(:,1), IM%PK%I%M%T%XCV(:,1), &
           IM%PK%I%IP%XRUNOFFD(:,1), &
           IM%PK%I%IP%XSOILWGHT(:,:,1), IM%I%O%NLAYER_HORT, IM%I%O%NLAYER_DUN, ZP_ALBNIR_TVEG, ZP_ALBVIS_TVEG,  &
           ZP_ALBNIR_TSOIL, ZP_ALBVIS_TSOIL, IM%PK%I%R%XSNOWFREE_ALB(:,1), IM%PK%I%M%T%XWRMAX_CF(:,1), &
           IM%PK%I%M%T%XVEG(:,1), IM%PK%I%M%T%XLAI(:,1), IM%PK%I%M%T%XEMIS(:,1), IM%PKD%DGIP%XZ0, &
           IM%PKD%DGIP%XZ0H, IM%PK%I%IP%XVEGTYPE_PATCH(:,:,1), IM%PKD%DGIP%XZ0EFF,   &
           IM%PK%I%M%M%XRGLGV(:,1), IM%PK%I%M%M%XGAMMAGV(:,1), IM%PK%I%M%M%XRSMINGV(:,1), &
           IM%PK%I%M%M%XROOTFRACGV(:,:,1), IM%PK%I%M%M%XWRMAX_CFGV(:,1), IM%PK%I%M%M%XLAIGV(:,1), &
           IM%PK%I%M%T%XBSLAI(:,1), &
           IM%PK%I%M%T%XLAIMIN(:,1),IM%PK%I%M%M%XH_VEG(:,1),ZPALPHAN, ZZ0G_WITHOUT_SNOW, ZZ0_MEBV,     &
           ZZ0H_MEBV,ZZ0EFF_MEBV, ZZ0_MEBN,ZZ0H_MEBN,ZZ0EFF_MEBN, IM%PK%I%M%M%XGNDLITTER(:,1),  &
           IM%PK%I%P%XRUNOFFB, IM%PK%I%IP%XCGSAT, IM%PK%I%IP%XC1SAT(:,1), IM%PK%I%IP%XC2REF(:,1), &
           IM%PK%I%IP%XC3(:,:,1), IM%PK%I%IP%XC4B, IM%PK%I%IP%XC4REF(:,1), IM%PK%I%IP%XACOEF, IM%PK%I%IP%XPCOEF, &
           IM%PK%I%IP%XTAUICE, IM%PK%I%P%XWDRAIN, ZP_TDEEP_A, IM%PK%I%IP%XTDEEP, IM%PK%I%IP%XGAMMAT,  &
           IM%PK%I%R%XPSN(:,1), IM%PK%I%R%XPSNG(:,1), IM%PK%I%R%XPSNV(:,1), IM%PK%I%R%XPSNV_A(:,1), &
           IM%PK%I%R%XSNOWFREE_ALB_VEG(:,1), IM%PK%I%R%XSNOWFREE_ALB_SOIL(:,1), IM%PK%I%M%I%XIRRIG(:,1), &
           IM%PK%I%M%I%XWATSUP(:,1), IM%PK%AG%XTHRESHOLDSPT(:,1), IM%PK%AG%LIRRIGATE(:,1), IM%PK%AG%LIRRIDAY(:,1), &
           IM%PK%I%M%T%LSTRESS(:,1), IM%PK%I%M%T%XGC(:,1), IM%PK%I%M%T%XF2I(:,1), IM%PK%I%M%X%XDMAX(:,1), &
           IM%PK%I%IP%XAH(:,1), &
           IM%PK%I%IP%XBH(:,1), ZP_CO2, IM%PK%I%M%T%XGMES(:,1), IM%I%IP%XPOI, IM%PK%I%IP%XFZERO(:,1), IM%PK%I%IP%XEPSO(:,1), &
           IM%PK%I%IP%XGAMM(:,1), IM%PK%I%IP%XQDGAMM(:,1), IM%PK%I%IP%XQDGMES(:,1), IM%PK%I%IP%XT1GMES(:,1), &
           IM%PK%I%IP%XT2GMES(:,1), &
           IM%PK%I%IP%XAMAX(:,1), IM%PK%I%IP%XQDAMAX(:,1), IM%PK%I%IP%XT1AMAX(:,1), IM%PK%I%IP%XT2AMAX(:,1), IM%I%IP%XABC, &
           IM%PK%I%M%X%XDG(:,:,1), IM%PK%I%IP%XDZG(:,:,1), IM%PK%I%IP%XDZDIF(:,:,1), &
           IM%PK%I%M%X%NWG_LAYER(:,1), IM%PK%I%M%X%XROOTFRAC(:,:,1), &
           IM%PK%I%IP%XWFC, IM%PK%I%IP%XWWILT, IM%PK%I%IP%XWSAT, IM%PK%I%IP%XBCOEF, IM%PK%I%IP%XCONDSAT(:,:,1), &
           IM%PK%I%IP%XMPOTSAT, IM%PK%I%IP%XHCAPSOIL, IM%PK%I%IP%XCONDDRY, IM%PK%I%IP%XCONDSLD, IM%PK%I%M%X%XD_ICE(:,1), &
           IM%PK%I%IP%XKSAT_ICE(:,1), IM%PK%I%I%XMUF, IM%PK%I%I%XFF(:,1), IM%PK%I%I%XFFG(:,1), IM%PK%I%I%XFFV(:,1), ZP_FFGNOS,  &
           ZP_FFVNOS, IM%PK%I%I%XFFROZEN(:,1), IM%PK%I%I%XALBF(:,1), IM%PK%I%I%XEMISF(:,1), IM%PK%I%I%XFFLOOD, IM%PK%I%I%XPIFLOOD, &
           IM%PKD%DGEIP%XIFLOOD, IM%PKD%DGEIP%XPFLOOD, IM%PKD%DGEIP%XLE_FLOOD, IM%PKD%DGEIP%XLEI_FLOOD, IM%I%O%XSODELX, &
           IM%PK%G%XLAT, IM%PK%G%XLON, IM%PK%I%R%XTG(:,:,1), IM%PK%I%R%XWG(:,:,1), IM%PK%I%R%XWGI(:,:,1), IM%PK%I%IP%XPCPS(:,1), &
           IM%PK%I%IP%XPLVTT(:,1), IM%PK%I%IP%XPLSTT(:,1), IM%PK%I%R%XWR(:,1), IM%PK%I%R%XWRL(:,1), &
           IM%PK%I%R%XWRLI(:,1), IM%PK%I%R%XWRVN(:,1), &
           IM%PK%I%R%XTV(:,1), IM%PK%I%R%XTL(:,1), &
           IM%PK%I%R%XRESA(:,1), IM%PK%I%R%XANFM(:,1), IM%PK%I%I%XFSAT, IM%PK%I%R%TSNOW%ALB(:,1), IM%PK%I%R%TSNOW%ALBVIS(:,1), &
           IM%PK%I%R%TSNOW%ALBNIR(:,1), IM%PK%I%R%TSNOW%ALBFIR(:,1), IM%PK%I%R%TSNOW%WSNOW(:,:,1), IM%PK%I%R%TSNOW%HEAT(:,:,1), &
           IM%PK%I%R%TSNOW%RHO(:,:,1), IM%PK%I%R%TSNOW%GRAN1(:,:,1), IM%PK%I%R%TSNOW%GRAN2(:,:,1), &
           IM%PK%I%R%TSNOW%HIST(:,:,1), IM%PK%I%R%TSNOW%AGE(:,:,1), &
           IM%PKD%DGMI%XGRNDFLUX, IM%PKD%DGMI%XHPSNOW, IM%PKD%DGMI%XSNOWHMASS, IM%PKD%DGMI%XRNSNOW, IM%PKD%DGMI%XHSNOW, &
           IM%PKD%DGMI%XGFLUXSNOW, IM%PKD%DGMI%XUSTARSNOW, IM%PKD%DGMI%XSRSFC, IM%PKD%DGMI%XRRSFC, IM%PKD%DGEIP%XLESL,   &
           IM%PK%I%R%TSNOW%EMIS(:,1), IM%PKD%DGMI%XCDSNOW, IM%PKD%DGMI%XCHSNOW, IM%PKD%DGIP%XTSRAD, IM%PKD%DGIP%XTS, &
           IM%PKD%DGMI%XHV, IM%PKD%DGIP%XQS, IM%PKD%DGMI%XSNOWTEMP, IM%PKD%DGMI%XSNOWLIQ, IM%PKD%DGMI%XSNOWDZ, &
           IM%PKD%DGMI%XCG, IM%PKD%DGMI%XC1, IM%PKD%DGMI%XC2, IM%PKD%DGMI%XWGEQ, IM%PKD%DGMI%XCT, IM%PKD%DGIP%XCH, &
           IM%PKD%DGIP%XCD, &
           IM%PKD%DGIP%XCDN, IM%PKD%DGIP%XRI, IM%PKD%DGIP%XHU, IM%PKD%DGIP%XHUG, ZP_EMIS, IM%PKD%DGIP%XALBT, IM%PKD%DGMI%XRS, &
           IM%PK%I%R%XLE(:,1), IM%PKD%DGIP%XRN, IM%PKD%DGIP%XH, IM%PKD%DGIP%XLEI, IM%PKD%DGEIP%XLEGI, IM%PKD%DGEIP%XLEG, &
           IM%PKD%DGEIP%XLEV, &
           IM%PKD%DGEIP%XLES, IM%PKD%DGEIP%XLER, IM%PKD%DGEIP%XLETR, IM%PKD%DGIP%XEVAP, &
           IM%PKD%DGIP%XGFLUX, IM%PKD%DGEIP%XRESTORE, &
           ZP_USTAR, IM%PKD%DGEIP%XDRAIN, IM%PKD%DGEIP%XRUNOFF, IM%PKD%DGEIP%XMELT, IM%PKD%DGEIP%XMELTADV, IM%PK%I%R%XTC(:,1), &
           IM%PK%I%R%XQC(:,1), IM%PKD%DGI%XRN, IM%PKD%DGI%XH, IM%PKD%DGEI%XLEG, IM%PKD%DGEI%XLEGI, &
           IM%PKD%DGEI%XLEV, IM%PKD%DGEI%XLETR, IM%PKD%DGEI%XUSTAR, IM%PKD%DGEI%XLER, IM%PKD%DGI%XLE, &
           IM%PKD%DGI%XLEI, IM%PKD%DGI%XGFLUX, IM%PKD%DGEIP%XHORT, IM%PKD%DGEIP%XDRIP, IM%PKD%DGEIP%XRRVEG, &
           ZP_AC_AGG, ZP_HU_AGG, IM%PK%I%R%XFAPARC(:,1), IM%PK%I%R%XFAPIRC(:,1), IM%PK%I%R%XMUS(:,1), &
           IM%PK%I%R%XLAI_EFFC(:,1), IM%PK%I%R%XAN(:,1),   &
           IM%PK%I%R%XANDAY(:,1), ZP_RESP_BIOMASS_INST, IM%PKD%GB%XIACAN(:,:,1), IM%PKD%DGEIP%XGPP, IM%PKD%DGMI%XFAPAR, &
           IM%PKD%DGMI%XFAPIR, IM%PKD%DGMI%XFAPAR_BS, IM%PKD%DGMI%XFAPIR_BS, IM%PKD%DGEIP%XIRRIG_FLUX, ZP_DEEP_FLUX,  &
           IM%PKD%DGEIP%XSWNET_V, IM%PKD%DGEIP%XSWNET_G, IM%PKD%DGEIP%XSWNET_N, IM%PKD%DGEIP%XSWNET_NS, IM%PKD%DGEIP%XLWNET_V, &
           IM%PKD%DGEIP%XLWNET_G, IM%PKD%DGEIP%XLWNET_N, IM%PKD%DGEIP%XLEVCV, IM%PKD%DGEIP%XLESC, IM%PKD%DGEIP%XH_V_C, &
           IM%PKD%DGEIP%XH_G_C, IM%PKD%DGEIP%XLETRGV, IM%PKD%DGEIP%XLETRCV, IM%PKD%DGEIP%XLERGV, IM%PKD%DGEIP%XLELITTER, &
           IM%PKD%DGEIP%XLELITTERI,IM%PKD%DGEIP%XDRIPLIT,IM%PKD%DGEIP%XRRLIT, IM%PKD%DGEIP%XLERCV, IM%PKD%DGEIP%XH_C_A, &
           IM%PKD%DGEIP%XH_N_C, IM%PKD%DGEIP%XLE_C_A, IM%PKD%DGEIP%XLE_V_C, IM%PKD%DGEIP%XLE_G_C,IM%PKD%DGEIP%XLE_N_C, &
           IM%PKD%DGEIP%XEVAP_N_C, IM%PKD%DGEIP%XEVAP_G_C, IM%PKD%DGEIP%XSR_GN, IM%PKD%DGEIP%XMELTCV, IM%PKD%DGEIP%XFRZCV,   &
           IM%PKD%DGEIP%XSWDOWN_GN, IM%PKD%DGEIP%XLWDOWN_GN, ZIRRIG_GR, IM%PK%I%I%XTOPQS(:,:,1), IM%PKD%DGEIP%XQSB, &
           IM%PKD%DGIP%XSUBL, &
           IM%PK%I%IP%XFWTD, IM%PK%I%IP%XWTD, IM%PKD%DGEIP%XSNDRIFT               )
!
ZP_TRAD=IM%PKD%DGIP%XTSRAD
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Glacier : ice runoff flux (especally for Earth System Model)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
IF(IM%I%O%LGLACIER)THEN
!           
  CALL HYDRO_GLACIER(IM%I%R%TSNOW%SCHEME, &
                     PTSTEP,ZP_SNOW,IM%PK%I%R%TSNOW%RHO(:,:,1),&
                     IM%PK%I%R%TSNOW%WSNOW(:,:,1),IM%PK%I%R%XICE_STO(:,1),IM%PKD%DGEIP%XICEFLUX)
!     
ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Calculation of ISBA water and energy budget (and time tendencies of each reservoir)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
CALL ISBA_BUDGET(IM%DGEI%LWATER_BUDGET, &
                 IM%I%O%CISBA,IM%I%R%TSNOW%SCHEME,IM%I%O%LGLACIER,PTSTEP,          &
                 IM%PK%I%R%XWG(:,:,1),IM%PK%I%R%XWGI(:,:,1),IM%PK%I%R%XWR(:,1),&
                 IM%PK%I%R%TSNOW%WSNOW(:,:,1),IM%PK%I%M%X%XDG(:,:,1),IM%PK%I%IP%XDZG(:,:,1),  & 
                 ZP_WG_INI,ZP_WGI_INI,ZP_WR_INI,ZP_SWE_INI,   &
                 ZP_RAIN,ZP_SNOW,IM%PKD%DGIP%XEVAP,IM%PKD%DGEIP%XDRAIN,IM%PKD%DGEIP%XRUNOFF,  &
                 IM%PKD%DGEIP%XIFLOOD,IM%PKD%DGEIP%XPFLOOD,IM%PKD%DGEIP%XICEFLUX,IM%PKD%DGEIP%XIRRIG_FLUX,&
                 IM%PKD%DGEIP%XSNDRIFT,                                  &
                 IM%PKD%DGEIP%XDWG,IM%PKD%DGEIP%XDWGI,IM%PKD%DGEIP%XDWR,IM%PKD%DGEIP%XDSWE,IM%PKD%DGEIP%XWATBUD      )
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Evolution of soil albedo, when depending on surface soil wetness:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
IF (IM%I%O%CALBEDO=='EVOL' .AND. IM%I%O%LECOCLIMAP) THEN
  CALL SOIL_ALBEDO(IM%I%O%CALBEDO,                                    &
                   IM%PK%I%IP%XWSAT(:,1),IM%PK%I%R%XWG(:,1,1),                    &
                   IM%PK%I%IP%XALBVIS_DRY,IM%PK%I%IP%XALBNIR_DRY,IM%PK%I%IP%XALBUV_DRY,   &
                   IM%PK%I%IP%XALBVIS_WET,IM%PK%I%IP%XALBNIR_WET,IM%PK%I%IP%XALBUV_WET,   &
                   IM%PK%I%M%A%XALBVIS_SOIL(:,1),IM%PK%I%M%A%XALBNIR_SOIL(:,1),IM%PK%I%M%A%XALBUV_SOIL(:,1) )  
  !
  CALL ALBEDO(IM%I%O%CALBEDO,                                          &
              IM%PK%I%M%T%XALBVIS_VEG,IM%PK%I%M%T%XALBNIR_VEG,IM%PK%I%M%T%XALBUV_VEG,IM%PK%I%M%T%XVEG,  &
              IM%PK%I%M%A%XALBVIS_SOIL,IM%PK%I%M%A%XALBNIR_SOIL,IM%PK%I%M%A%XALBUV_SOIL,      &
              IM%PK%I%M%T%XALBVIS,IM%PK%I%M%T%XALBNIR,IM%PK%I%M%T%XALBUV                      )  
END IF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Vegetation evolution for interactive LAI
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
IF (IM%I%O%CPHOTO=='LAI' .OR. IM%I%O%CPHOTO=='LST' .OR. IM%I%O%CPHOTO=='NIT' .OR. IM%I%O%CPHOTO=='NCB') THEN
  CALL VEGETATION_EVOL(IM%I%O%CISBA, IM%I%O%CPHOTO, IM%I%O%CRESPSL, IM%I%O%CALBEDO, LAGRIP, IM%I%O%LTR_ML,           &
                       IM%I%O%LNITRO_DILU, IM%I%O%LAGRI_TO_GRASS,                               &
                       PTSTEP, KMONTH, KDAY, IM%I%O%NSPINW, PTIME, IM%PK%G%XLAT, ZP_RHOA,      &
                       IM%PK%I%M%X%XDG(:,:,1), IM%PK%I%IP%XDZG(:,:,1), IM%PK%I%M%X%NWG_LAYER(:,1),                        &  
                       IM%PK%I%R%XTG(:,:,1), IM%PK%I%M%T%XALBNIR_VEG(:,1), &
                       IM%PK%I%M%T%XALBVIS_VEG(:,1), IM%PK%I%M%T%XALBUV_VEG(:,1),         &
                       IM%PK%I%M%A%XALBNIR_SOIL(:,1), IM%PK%I%M%A%XALBVIS_SOIL(:,1), IM%PK%I%M%A%XALBUV_SOIL(:,1),             &
                       IM%PK%I%IP%XVEGTYPE_PATCH(:,:,1), IM%PK%I%M%T%XSEFOLD(:,1), IM%PK%I%IP%XANMAX(:,1), &
                       IM%PK%I%M%X%XH_TREE(:,1), IM%PK%I%M%T%XBSLAI(:,1),&
                       IM%PK%I%M%T%XLAIMIN(:,1), ZP_CO2, IM%PK%I%M%T%XCE_NITRO(:,1), &
                       IM%PK%I%M%T%XCF_NITRO(:,1), IM%PK%I%M%T%XCNA_NITRO(:,1), &
                       IM%PK%I%IP%XBSLAI_NITRO(:,1), IM%PK%I%M%T%XGMES(:,1), IM%PK%I%IP%XTAU_WOOD(:,1), &
                       IM%PK%I%M%I%TSEED(:,1),             &
                       IM%PK%I%M%I%TREAP(:,1), IM%PK%I%P%XAOSIP, IM%PK%I%P%XAOSIM, IM%PK%I%P%XAOSJP, IM%PK%I%P%XAOSJM,           &
                       IM%PK%I%P%XHO2IP, IM%PK%I%P%XHO2IM, IM%PK%I%P%XHO2JP, IM%PK%I%P%XHO2JM, IM%PK%I%IP%XZ0EFFIP(:,1),        &
                       IM%PK%I%IP%XZ0EFFIM(:,1), IM%PK%I%IP%XZ0EFFJP(:,1), IM%PK%I%IP%XZ0EFFJM(:,1), IM%PK%I%M%T%XLAI(:,1), &
                       IM%PK%I%M%T%XVEG(:,1),        &
                       IM%PK%I%M%T%XZ0(:,1), IM%PK%I%M%T%XALBNIR(:,1), IM%PK%I%M%T%XALBVIS(:,1), IM%PK%I%M%T%XALBUV(:,1), &
                       IM%PK%I%M%T%XEMIS(:,1),            &
                       IM%PK%I%R%XANFM(:,1), IM%PK%I%R%XANDAY(:,1), IM%PK%I%R%XBIOMASS(:,:,1), IM%PK%I%R%XRESP_BIOMASS(:,:,1),  &
                       ZP_RESP_BIOMASS_INST, IM%PK%I%IP%XINCREASE(:,:,1), IM%PK%I%IP%XTURNOVER(:,:,1), &
                       ! add optional for accurate dependency to nitrogen
                       ! limitation
                        PSWDIR=ZP_GLOBAL_SW ) 
END IF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Diagnostic of respiration carbon fluxes and soil carbon evolution
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
ZP_SFCO2    (:)=0.
IM%PKD%DGEIP%XRESP_ECO (:)=0.
IM%PKD%DGEIP%XRESP_AUTO(:)=0.
!
IF ( IM%I%O%CPHOTO/='NON' .AND. IM%I%O%CRESPSL/='NON' .AND. ANY(IM%PK%I%M%T%XLAI(:,1)/=XUNDEF) ) THEN
  CALL CARBON_EVOL(IM%I%O%CISBA, IM%I%O%CRESPSL, IM%I%O%CPHOTO, PTSTEP, IM%I%O%NSPINS,                   &
                   ZP_RHOA, IM%PK%I%R%XTG(:,:,1), IM%PK%I%R%XWG(:,:,1), &
                   IM%PK%I%IP%XWFC, IM%PK%I%IP%XWWILT, IM%PK%I%IP%XWSAT, IM%PK%I%P%XSAND,&
                   IM%PK%I%M%X%XDG(:,:,1), IM%PK%I%IP%XDZG(:,:,1), IM%PK%I%M%X%NWG_LAYER(:,1),             &                   
                   IM%PK%I%M%X%XRE25(:,1), IM%PK%I%M%T%XLAI(:,1), ZP_RESP_BIOMASS_INST, IM%PK%I%IP%XTURNOVER(:,:,1),       &
                   IM%PK%I%R%XLITTER(:,:,:,1), IM%PK%I%R%XLIGNIN_STRUC(:,:,1) , IM%PK%I%R%XSOILCARB(:,:,1),                 &
                   IM%PKD%DGEIP%XRESP_AUTO, IM%PKD%DGEIP%XRESP_ECO                                 )  
  ! calculation of vegetation CO2 flux
  ! Positive toward the atmosphere
  ZP_SFCO2(:) = IM%PKD%DGEIP%XRESP_ECO(:) - IM%PKD%DGEIP%XGPP(:)  
END IF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Reset effecitve roughness lentgh to its nominal value when snow has just disappeared
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
 CALL SUBSCALE_Z0EFF(IM%PK%I%P%XAOSIP,IM%PK%I%P%XAOSIM,IM%PK%I%P%XAOSJP,IM%PK%I%P%XAOSJM,            &
                    IM%PK%I%P%XHO2IP,IM%PK%I%P%XHO2IM,IM%PK%I%P%XHO2JP,IM%PK%I%P%XHO2JM,IM%PK%I%M%T%XZ0(:,1),      &
                    IM%PK%I%IP%XZ0EFFIP(:,1),IM%PK%I%IP%XZ0EFFIM(:,1),IM%PK%I%IP%XZ0EFFJP(:,1),IM%PK%I%IP%XZ0EFFJM(:,1),    &
                    OMASK=(IM%PK%I%R%TSNOW%WSNOW(:,1,1)==0. .AND. IM%PK%I%R%XPSN(:,1)>0.)  )   
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Turbulent fluxes
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
ZP_SFTH(:) = IM%PKD%DGIP%XH(:)
ZP_SFTQ(:) = IM%PKD%DGIP%XEVAP(:)

ZP_SFU (:) = 0.
ZP_SFV (:) = 0.
WHERE (ZP_WIND>0.)
  ZP_SFU (:) = - ZP_U(:)/ZP_WIND(:) * ZP_USTAR(:)**2 * ZP_RHOA(:)
  ZP_SFV (:) = - ZP_V(:)/ZP_WIND(:) * ZP_USTAR(:)**2 * ZP_RHOA(:)
END WHERE
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Scalar fluxes
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
ZP_SFTS(:,:) = 0.
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
! --------------------------------------------------------------------------------------
! Chemical dry deposition :
! --------------------------------------------------------------------------------------
IF (IM%CHI%SVI%NBEQ>0) THEN
  IF( IM%CHI%CCH_DRY_DEP == "WES89") THEN

    CALL CH_DEP_ISBA    (ZP_USTAR, IM%PKD%DGIP%XHU, IM%PK%I%R%XPSN(:,1),             &
                        IM%PK%I%M%T%XVEG(:,1), IM%PK%I%M%T%XLAI(:,1), IM%PK%I%P%XSAND, IM%PK%I%P%XCLAY, &
                        IM%PK%I%R%XRESA(:,1), &
                        IM%PKD%DGMI%XRS(:),  IM%PK%I%M%T%XZ0(:,1),                       &
                        ZP_TA, ZP_PA, ZP_TRAD(:),                  &
                        IM%PK%I%IP%XVEGTYPE_PATCH(:,NVT_NO,1),                &
                        IM%PK%I%IP%XVEGTYPE_PATCH(:,NVT_ROCK,1),              &
                        IM%CHI%SVI%CSV(IM%CHI%SVI%NSV_CHSBEG:IM%CHI%SVI%NSV_CHSEND),                &
                        IM%PKCI%XP_SOILRC_SO2,  IM%PKCI%XP_SOILRC_O3 ,             &
                        IM%PKCI%XP_DEP(:,1:IM%CHI%SVI%NBEQ)                           )  
 
    ZP_SFTS(:,IM%CHI%SVI%NSV_CHSBEG:IM%CHI%SVI%NSV_CHSEND) = - ZP_SV(:,IM%CHI%SVI%NSV_CHSBEG:IM%CHI%SVI%NSV_CHSEND)  &
                                                    * IM%PKCI%XP_DEP(:,1:IM%CHI%SVI%NBEQ)  
    IF (IM%CHI%SVI%NAEREQ > 0 ) THEN
      CALL CH_AER_DEP(ZP_SV(:,IM%CHI%SVI%NSV_AERBEG:IM%CHI%SVI%NSV_AEREND),&
                           ZP_SFTS(:,IM%CHI%SVI%NSV_AERBEG:IM%CHI%SVI%NSV_AEREND),&
                           ZP_USTAR, IM%PK%I%R%XRESA(:,1),ZP_TA,ZP_RHOA)     
    END IF
  ELSE
    ZP_SFTS(:,IM%CHI%SVI%NSV_CHSBEG:IM%CHI%SVI%NSV_CHSEND) = 0.
    ZP_SFTS(:,IM%CHI%SVI%NSV_AERBEG:IM%CHI%SVI%NSV_AEREND) = 0.
  ENDIF
ENDIF
!
! --------------------------------------------------------------------------------------
! Dust deposition and emission:
! --------------------------------------------------------------------------------------
!
IF(IM%CHI%SVI%NDSTEQ>0)THEN
  IDST = IM%CHI%SVI%NSV_DSTEND - IM%CHI%SVI%NSV_DSTBEG + 1

  CALL COUPLING_DST_n(DST, IM%PK%I%IP%XVEGTYPE_PATCH(:,:,1), &
            HPROGRAM,                    &!I [char] Name of program
            KSIZE,      &!I [nbr] number of points in patch
            IDST,                        &!I [nbr] number of dust emissions variables
            JPATCH,                      &!I [idx] patch in question
            IM%PK%I%P%XCLAY(:,1),                &!I [frc] mass fraction clay in first soil layer
            ZP_PS,                       &!I [Pa] surface pressure
            IM%PKD%DGIP%XTS,                       &!I [K] surface temperature
            ZP_QA,                       &!I [kg/kg] specific humidity
            IM%PK%I%R%XRESA(:,1),                     &!I [s/m] atmospheric resistance
            ZP_RHOA,                     &!I [kg/m3] atmospheric density
            IM%PK%I%P%XSAND(:,1),                &!I [frc] mass fraction of sand in first soil layer
            ZP_PA,                       &!I [K] Atmospheric pressure
            ZP_TA,                       &!I [K] Atmospheric temperature
            IM%PK%I%R%XTG(:,1,1),                  &!I [K] Ground temperature
            ZP_U,                        &!I [m/s] zonal wind at atmospheric height 
            ZP_UREF,                     &!I [m] reference height of wind
            ZP_V,                        &!I [m/s] meridional wind at atmospheric height
            IM%PK%I%R%XWG(:,1,1),                  &!I [m3/m3] ground volumetric water content
            IM%PK%I%IP%XWSAT(:,1),                &!I [m3/m3] saturation volumetric water content
            ZP_ZREF,                     &!I [m] reference height of wind
            IM%PKD%DGIP%XCD,                       &!I [] Drag Coefficient for momentum
            IM%PKD%DGIP%XCDN,                      &!I [] Drag neutral Coefficient for momentum
            IM%PKD%DGIP%XCH,                       &!I [] drag coefficient for heat
            IM%PKD%DGIP%XRI,                       &!I [] Richardson number
            IM%PKD%DGIP%XZ0H,            &!I [frc] Z0 (heat) with snow
            ZP_SFTS(:,IM%CHI%SVI%NSV_DSTBEG:IM%CHI%SVI%NSV_DSTEND)  &!O [kg/m2/sec] flux of dust            
            )  
!
   IF (IM%CHI%SVI%NSV_AEREND > 0)  THEN ! case of dust/ anthropogenic aerosols coupling
     DO JMODE=1,NDSTMDE
       !
       !Make index which is 0 for first mode, 3 for second, 6 for third etc
       IF (LVARSIG_DST) THEN
         JSV_IDX = (JMODE-1)*3
       ELSE IF (LRGFIX_DST) THEN
         JSV_IDX = JMODE-2
       ELSE
         JSV_IDX = (JMODE-1)*2
       END IF
       !
       DO JSV=1, size(HSV)
         IF ((TRIM(HSV(JSV)) == "@DSTI").AND.(JMODE==3)) THEN 
           ! add dust flux and conversion kg/m2/s into molec.m2/s
           ZP_SFTS(:,JSV) = ZP_SFTS(:,JSV) + ZP_SFTS(:,IM%CHI%SVI%NSV_DSTBEG-1+JSV_IDX+2)*XAVOGADRO/XMOLARWEIGHT_DST
         END IF
         IF ( (TRIM(HSV(JSV)) == "@DSTJ").AND.(JMODE==2)) THEN 
           ! add dust flux and conversion kg/m2/sec into molec.m2/s
           ZP_SFTS(:,JSV) = ZP_SFTS(:,JSV) + ZP_SFTS(:,IM%CHI%SVI%NSV_DSTBEG-1+JSV_IDX+2)*XAVOGADRO/XMOLARWEIGHT_DST
         END IF
       END DO
       !
     END DO
    END IF
!    
!Modify fluxes due to dry deposition, we introduce a negative flux where dust is lost
  CALL DSLT_DEP(ZP_SV(:,IM%CHI%SVI%NSV_DSTBEG:IM%CHI%SVI%NSV_DSTEND), &
                ZP_SFTS(:,IM%CHI%SVI%NSV_DSTBEG:IM%CHI%SVI%NSV_DSTEND), &
                ZP_USTAR, IM%PK%I%R%XRESA(:,1), ZP_TA, ZP_RHOA, DST%XEMISSIG_DST, DST%XEMISRADIUS_DST, &
                JPMODE_DST, XDENSITY_DST, XMOLARWEIGHT_DST, ZCONVERTFACM0_DST,    &
                ZCONVERTFACM6_DST, ZCONVERTFACM3_DST, LVARSIG_DST, LRGFIX_DST,    &
                CVERMOD  )
!
!Transfer these fluxes to fluxes understandable by all moments
  CALL MASSFLUX2MOMENTFLUX(           &
    ZP_SFTS(:,IM%CHI%SVI%NSV_DSTBEG:IM%CHI%SVI%NSV_DSTEND), & !I/O ![kg/m2/sec] In: flux of only mass, out: flux of moments
    ZP_RHOA,                          & !I [kg/m3] air density
    DST%XEMISRADIUS_DST,                  &!I [um] emitted radius for the modes (max 3)
    DST%XEMISSIG_DST,                     &!I [-] emitted sigma for the different modes (max 3)
    NDSTMDE,                          &
    ZCONVERTFACM0_DST,                &
    ZCONVERTFACM6_DST,                &
    ZCONVERTFACM3_DST,                &
    LVARSIG_DST, LRGFIX_DST           )   
!
ENDIF !Check on CDSTYN
!
! --------------------------------------------------------------------------------------
! Sea Salt deposition
! --------------------------------------------------------------------------------------
!
IF (IM%CHI%SVI%NSLTEQ>0) THEN
  CALL DSLT_DEP(ZP_SV(:,IM%CHI%SVI%NSV_SLTBEG:IM%CHI%SVI%NSV_SLTEND), &
                ZP_SFTS(:,IM%CHI%SVI%NSV_SLTBEG:IM%CHI%SVI%NSV_SLTEND), &
                ZP_USTAR, IM%PK%I%R%XRESA(:,1), ZP_TA, ZP_RHOA, SLT%XEMISSIG_SLT, SLT%XEMISRADIUS_SLT, &
                JPMODE_SLT, XDENSITY_SLT, XMOLARWEIGHT_SLT, ZCONVERTFACM0_SLT,    &
                ZCONVERTFACM6_SLT, ZCONVERTFACM3_SLT, LVARSIG_SLT, LRGFIX_SLT,    &
                CVERMOD  )  

  CALL MASSFLUX2MOMENTFLUX(           &
    ZP_SFTS(:,IM%CHI%SVI%NSV_SLTBEG:IM%CHI%SVI%NSV_SLTEND), & !I/O ![kg/m2/sec] In: flux of only mass, out: flux of moments
    ZP_RHOA,                          & !I [kg/m3] air density
    SLT%XEMISRADIUS_SLT,                  &!I [um] emitted radius for the modes (max 3)
    SLT%XEMISSIG_SLT,                     &!I [-] emitted sigma for the different modes (max 3)
    NSLTMDE,                          &
    ZCONVERTFACM0_SLT,                &
    ZCONVERTFACM6_SLT,                &
    ZCONVERTFACM3_SLT,                &
    LVARSIG_SLT, LRGFIX_SLT         ) 
ENDIF !Check on CSLTYN
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Inline diagnostics
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
 CALL DIAG_INLINE_ISBA_n(IM%DGI, IM%PKD, IM%DGEI%LSURF_BUDGETC, IM%I%O%LCANOPY, &
                         ZP_TA, IM%PKD%DGIP%XTS, ZP_QA, ZP_PA, ZP_PS, ZP_RHOA, ZP_U, ZP_V,       &
                          ZP_ZREF, ZP_UREF,                                            &
                          IM%PKD%DGIP%XCD, IM%PKD%DGIP%XCDN, IM%PKD%DGIP%XCH, IM%PKD%DGIP%XRI, &
                          IM%PKD%DGIP%XHU, IM%PKD%DGIP%XZ0,         &
                          IM%PKD%DGIP%XZ0H, IM%PKD%DGIP%XZ0EFF,                                  &
                          ZP_SFTH, ZP_SFTQ, ZP_SFU, ZP_SFV, IM%PKD%DGIP%XQS,                     &
                          IM%PK%I%I%XDIR_ALB_WITH_SNOW(:,:,1), IM%PK%I%I%XSCA_ALB_WITH_SNOW(:,:,1),                  &
                          ZP_DIR_SW, ZP_SCA_SW, ZP_LW, IM%PKD%DGIP%XRN                           )  
!
!
!-------------------------------------------------------------------------------
!Physical properties see by the atmosphere in order to close the energy budget 
!between surfex and the atmosphere. All variables should be at t+1 but very 
!difficult to do. Maybe it will be done later. However, Ts can be at time t+1
!-------------------------------------------------------------------------------
!
ZP_TSURF (:) = IM%PKD%DGIP%XTS (:)
ZP_Z0    (:) = IM%PKD%DGIP%XZ0 (:)
ZP_Z0H   (:) = IM%PKD%DGIP%XZ0H(:)
ZP_QSURF (:) = IM%PKD%DGIP%XQS (:)
!
!-------------------------------------------------------------------------------
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Isba offline diagnostics for each patch
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
 CALL DIAG_EVAP_ISBA_n(IM%DGEI, IM%DGEIP, IM%DGEIPC, &
                       IM%DGIPC, IM%PKD, IM%PK, IM%I%O%LGLACIER, &
                       IM%I%O%LMEB_PATCH, IM%I%R%TSNOW%SCHEME, &
                       IM%I%O%CPHOTO,PTSTEP,KMASK,KSIZE,JPATCH,ZP_RHOA)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Isba offline diagnostics for miscellaneous terms over each patch
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
 CALL DIAG_MISC_ISBA_n(IM%DGMIP, IM%PKD, IM%DGMI%LSURF_MISC_BUDGET, &
                       PTSTEP, IM%I%O%CISBA, IM%I%O%CPHOTO, IM%I%R%TSNOW%SCHEME, LAGRIP, IM%I%O%LTR_ML,    &
                      PTIME, KSIZE, JPATCH, KMASK, IM%PK%AG%XTHRESHOLDSPT(:,1),              &
                      IM%PK%I%R%XPSN(:,1), IM%PK%I%R%XPSNG(:,1), IM%PK%I%R%XPSNV(:,1), &
                      IM%PK%I%I%XFF(:,1), IM%PK%I%I%XFFG(:,1), &
                      IM%PK%I%I%XFFV(:,1),        &
                      IM%PK%I%R%XWG(:,:,1), IM%PK%I%R%XWGI(:,:,1), IM%PK%I%IP%XWFC, IM%PK%I%IP%XWWILT, &
                      IM%PK%I%R%TSNOW%WSNOW(:,:,1), IM%PK%I%R%TSNOW%RHO(:,:,1),&
                      IM%PK%I%R%XFAPARC(:,1), IM%PK%I%R%XFAPIRC(:,1), &
                      IM%PK%I%R%XLAI_EFFC(:,1), IM%PK%I%R%XMUS(:,1), IM%PK%I%I%XFSAT,     &
                      IM%PK%I%M%X%XDG(:,:,1), IM%PK%I%R%XTG(:,:,1)       )                  
!
! Unpack ISBA diagnostics (modd_diag_isban) for each patch:ISIZE_MAX = MAXVAL(NSIZE_NATURE_P)

!  (MUST be done BEFORE UNPACK_ISBA_PATCH, because of XP_LE)
!
 CALL UNPACK_DIAG_PATCH_n(IM%DGI, IM%DGIP, IM%GB, IM%I, IM%PKD, IM%PK, &
                          KMASK,KSIZE,IM%I%O%NPATCH,JPATCH, &
                           ZCPL_DRAIN,ZCPL_RUNOFF,ZCPL_EFLOOD,ZCPL_PFLOOD,           &
                           ZCPL_IFLOOD, ZCPL_ICEFLUX)  
!
! for chemical deposition
!
IF (IM%CHI%SVI%NBEQ>0) THEN
  IF( IM%CHI%CCH_DRY_DEP == "WES89") THEN
    CALL UNPACK_CH_ISBA_PATCH_n(IM%CHI, IM%PKCI, &
                                KMASK,KSIZE,IM%I%O%NPATCH,JPATCH)     
  END IF
END IF
!
! Unpack ISBA variables (modd_isban) for each patch:
!
 CALL UNPACK_ISBA_PATCH_n(IM%AG, IM%I, IM%PK, &
                          KMASK,KSIZE,JPATCH)
!
!----------------------------------------------------------------------
!
! for further chemical biogenic emissions
!
IF (IM%CHI%SVI%NBEQ>0 .AND. IM%CHI%LCH_BIO_FLUX) THEN
  !
  DO JJ=1,KSIZE
    ZSW_FORBIO(KMASK(JJ),JPATCH) = 0.
  ENDDO
  !
  DO JSWB=1,ISWB
!cdir nodep
!cdir unroll=8
    DO JJ=1,KSIZE
      ZSW_FORBIO(KMASK(JJ),JPATCH) = ZSW_FORBIO(KMASK(JJ),JPATCH)              &
                                     + ZP_DIR_SW(JJ,JSWB) + ZP_SCA_SW(JJ,JSWB)  
    ENDDO
  ENDDO
  !
ENDIF
!----------------------------------------------------------------------
!
! Unpack output dummy arguments for each patch:
!
IF (IM%I%O%NPATCH==1) THEN
   ZSFTQ_TILE      (:,JPATCH)  = ZP_SFTQ      (:)
   ZSFTH_TILE      (:,JPATCH)  = ZP_SFTH      (:)
   ZSFTS_TILE      (:,:,JPATCH)= ZP_SFTS      (:,:)
   ZSFCO2_TILE     (:,JPATCH)  = ZP_SFCO2     (:)
   ZSFU_TILE       (:,JPATCH)  = ZP_SFU       (:)
   ZSFV_TILE       (:,JPATCH)  = ZP_SFV       (:)
   ZTRAD_TILE      (:,JPATCH)  = ZP_TRAD      (:)
   ZTSURF_TILE     (:,JPATCH)  = ZP_TSURF     (:)
   ZZ0_TILE        (:,JPATCH)  = ZP_Z0        (:)
   ZZ0H_TILE       (:,JPATCH)  = ZP_Z0H       (:)
   ZQSURF_TILE     (:,JPATCH)  = ZP_QSURF     (:)   
ELSE
!cdir nodep
!cdir unroll=8
 DO JJ=1,KSIZE
   JI = KMASK(JJ)
   ZSFTQ_TILE      (JI,JPATCH)  = ZP_SFTQ      (JJ)
   ZSFTH_TILE      (JI,JPATCH)  = ZP_SFTH      (JJ)
   ZSFCO2_TILE     (JI,JPATCH)  = ZP_SFCO2     (JJ)
   ZSFU_TILE       (JI,JPATCH)  = ZP_SFU       (JJ)
   ZSFV_TILE       (JI,JPATCH)  = ZP_SFV       (JJ)
   ZTRAD_TILE      (JI,JPATCH)  = ZP_TRAD      (JJ)
   ZTSURF_TILE     (JI,JPATCH)  = ZP_TSURF     (JJ)
   ZZ0_TILE        (JI,JPATCH)  = ZP_Z0        (JJ)
   ZZ0H_TILE       (JI,JPATCH)  = ZP_Z0H       (JJ)
   ZQSURF_TILE     (JI,JPATCH)  = ZP_QSURF     (JJ)   
 ENDDO
!
!cdir nodep
!cdir unroll=8
  DO JK=1,SIZE(ZP_SFTS,2)
    DO JJ=1,KSIZE
      JI=KMASK(JJ)    
      ZSFTS_TILE      (JI,JK,JPATCH)= ZP_SFTS      (JJ,JK)
    ENDDO
  ENDDO
ENDIF
!
!----------------------------------------------------------------------
!
! Get output dust flux if we are calculating dust
IF (NDSTMDE .GE. 1) IMOMENT = INT(IDST / NDSTMDE)
IF (IM%CHI%SVI%NDSTEQ>0) THEN
  DO JSV = 1,NDSTMDE
    IF (IMOMENT == 1) THEN
      DST%XSFDST(:,JSV,JPATCH)=ZSFTS_TILE(:,NDST_MDEBEG+JSV-1,JPATCH)
    ELSE
      DST%XSFDST(:,JSV,JPATCH)=ZSFTS_TILE(:,NDST_MDEBEG+(JSV-1)*IMOMENT+1,JPATCH)
    END IF

    DST%XSFDSTM(:,JSV,JPATCH)=DST%XSFDSTM(:,JSV,JPATCH) + DST%XSFDST(:,JSV,JPATCH) * PTSTEP
  ENDDO
ENDIF
!
IF (LHOOK) CALL DR_HOOK('COUPLING_ISBA_n:TREAT_PATCH',1,ZHOOK_HANDLE)
!
END SUBROUTINE TREAT_PATCH
!==========================================================================================
END SUBROUTINE COUPLING_ISBA_n
