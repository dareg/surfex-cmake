!     ###############################################################################
SUBROUTINE COUPLING_ISBA_n(HPROGRAM, HCOUPLING,                                              &
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
!!-------------------------------------------------------------------
!
USE MODD_REPROD_OPER, ONLY : CIMPLICIT_WIND
!
USE MODD_CSTS,         ONLY : XRD, XRV, XP00, XCPD, XPI, XAVOGADRO, XMD
USE MODD_CO2V_PAR,     ONLY : XMCO2, XSPIN_CO2
!
USE MODD_SURF_PAR,     ONLY : XUNDEF
USE MODD_SNOW_PAR,     ONLY : XZ0SN
USE MODD_TYPE_DATE_SURF
USE MODD_ISBA_n, ONLY : I => ISBA
!
USE MODD_SURF_ATM,    ONLY : LNOSOF
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
USE MODD_DST_n, ONLY : DST => DST
USE MODD_SLT_n, ONLY : SLT => SLT
USE MODD_DST_SURF
USE MODD_SLT_SURF
USE MODE_DSLT_SURF
USE MODE_MEB
!
USE MODD_PACK_ISBA, ONLY : PKI => PACK_ISBA
!
USE MODD_PACK_DIAG_ISBA, ONLY : PKDI => PACK_DIAG_ISBA
!                         
USE MODD_PACK_CH_ISBA, ONLY : PKCI => PACK_CH_ISBA

USE MODD_CH_ISBA_n, ONLY : CHI => CH_ISBA
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
REAL, DIMENSION(KI), INTENT(IN)  :: PCO2      ! CO2 concentration in the air          (kg/kg)
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
REAL, DIMENSION(KI,I%NPATCH) :: ZSFTH_TILE     ! surface heat flux (W/m2)
REAL, DIMENSION(KI,I%NPATCH) :: ZSFTQ_TILE     ! surface vapor flux (kg/m2/s)
REAL, DIMENSION(KI,I%NPATCH) :: ZSFCO2_TILE    ! surface CO2 flux positive toward the atmosphere (m/s*kg_CO2/kg_air)
REAL, DIMENSION(KI,I%NPATCH) :: ZSFU_TILE      ! zonal momentum flux
REAL, DIMENSION(KI,I%NPATCH) :: ZSFV_TILE      ! meridian momentum flux
REAL, DIMENSION(KI,I%NPATCH) :: ZTRAD_TILE     ! radiative surface temperature
REAL, DIMENSION(KI,I%NPATCH) :: ZEMIS_TILE     ! emissivity
REAL, DIMENSION(KI,I%NPATCH) :: ZTSURF_TILE    ! surface effective temperature
REAL, DIMENSION(KI,I%NPATCH) :: ZZ0_TILE       ! roughness length for momentum
REAL, DIMENSION(KI,I%NPATCH) :: ZZ0H_TILE      ! roughness length for heat
REAL, DIMENSION(KI,I%NPATCH) :: ZQSURF_TILE    ! specific humidity at surface
REAL, DIMENSION(KI,KSW,I%NPATCH) :: ZDIR_ALB_TILE  ! direct albedo
REAL, DIMENSION(KI,KSW,I%NPATCH) :: ZSCA_ALB_TILE  ! diffuse albedo
REAL, DIMENSION(KI,KSV,I%NPATCH) :: ZSFTS_TILE     ! scalar surface flux
!
REAL, DIMENSION(KI, I%NPATCH) :: ZCPL_DRAIN     ! For the coupling with TRIP
REAL, DIMENSION(KI, I%NPATCH) :: ZCPL_RUNOFF    ! For the coupling with TRIP
REAL, DIMENSION(KI, I%NPATCH) :: ZCPL_EFLOOD    ! For the coupling with TRIP
REAL, DIMENSION(KI, I%NPATCH) :: ZCPL_PFLOOD    ! For the coupling with TRIP
REAL, DIMENSION(KI, I%NPATCH) :: ZCPL_IFLOOD    ! For the coupling with TRIP
REAL, DIMENSION(KI, I%NPATCH) :: ZCPL_ICEFLUX
!
! for chemical computations
!
REAL, DIMENSION(KI, I%NPATCH) :: ZSW_FORBIO
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
  ZALFA(JJ) = ZDIR(JJ) - I%XZ0EFFJPDIR(JJ) * XPI/180.

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
IF (LAGRIP .AND. (I%CPHOTO=='LAI' .OR. I%CPHOTO=='LST' .OR. I%CPHOTO=='NIT'.OR. I%CPHOTO=='NCB') ) THEN
   CALL IRRIGATION_UPDATE(I%XIRRIG,PTSTEP,KMONTH,KDAY,PTIME,               &
                            I%TSEED(:,:)%TDATE%MONTH,I%TSEED(:,:)%TDATE%DAY,   &
                            I%TREAP(:,:)%TDATE%MONTH,I%TREAP(:,:)%TDATE%DAY    )  
ENDIF
!
!* Actualization of the SGH variable (Fmu, Fsat)
!
 CALL ISBA_SGH_UPDATE(I%CISBA,I%CRUNOFF,I%CRAIN,PRAIN,I%XMUF,I%XFSAT,I%XTOPQS)
!
!
!* Actualization of deep soil characteristics
!
IF (LDEEPSOIL) THEN
   CALL DEEPSOIL_UPDATE(I%TTIME%TDATE%MONTH)
ENDIF
!
!* Actualization of soil and wood carbon spinup
!
! During soil carbon spinup with ISBA-CC: 
!        (1) Atmospheric CO2 concentration fixed to Pre-industrial CO2 consentration XCO2_START
!        (2) Atmospheric CO2 concentration rampin up from XCO2_START to XCO2_END
!
IF(I%LSPINUPCARBS.OR.I%LSPINUPCARBW)THEN
!
  ISPINEND=I%NNBYEARSPINS-NINT(I%NNBYEARSPINS*XSPIN_CO2)
!  
  I%LAGRI_TO_GRASS = .FALSE.
!
  IF ( I%LSPINUPCARBS .AND. (I%NNBYEARSOLD <= ISPINEND) ) THEN
!
   I%LAGRI_TO_GRASS = .TRUE.
!
   ZCO2(:) = I%XCO2_START * 1.E-6 * XMCO2 / XMD
!
  ELSEIF(I%LSPINUPCARBS .AND. (I%NNBYEARSOLD > ISPINEND) .AND. (I%NNBYEARSOLD <= I%NNBYEARSPINS) )THEN
!
   ZSPINCO2 = I%XCO2_START + (I%XCO2_END-I%XCO2_START) * REAL(I%NNBYEARSOLD - ISPINEND) / REAL(I%NNBYEARSPINS - ISPINEND)
!
   ZCO2 (:) = ZSPINCO2 * 1.E-6 * XMCO2 / XMD
!
  ENDIF
!
  CALL CARBON_SPINUP(I%TTIME%TDATE%MONTH,I%TTIME%TDATE%DAY,I%TTIME%TIME,       &
                     I%LSPINUPCARBS, I%LSPINUPCARBW, I%XSPINMAXS, I%XSPINMAXW,   &
                     I%NNBYEARSPINS, I%NNBYEARSPINW, I%NNBYEARSOLD, I%CPHOTO,    &
                     I%CRESPSL, I%NSPINS, I%NSPINW                             )
!
ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Time evolution
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
I%TTIME%TIME = I%TTIME%TIME + PTSTEP
 CALL ADD_FORECAST_TO_DATE_SURF(I%TTIME%TDATE%YEAR,I%TTIME%TDATE%MONTH,I%TTIME%TDATE%DAY,I%TTIME%TIME)
!
! --------------------------------------------------------------------------------------
!
!*      2.     Physical evolution
!
! --------------------------------------------------------------------------------------
! Patch Dependent Calculations
! --------------------------------------------------------------------------------------
!
PATCH_LOOP: DO JPATCH=1,I%NPATCH
!
  IF (I%NSIZE_NATURE_P(JPATCH) == 0 ) CYCLE
!
! Pack dummy arguments for each patch:
!
#ifdef TOPD
  IF (LCOUPL_TOPD)NMASKT_PATCH(:)=I%NR_NATURE_P(:,JPATCH)
#endif
  CALL TREAT_PATCH(I%NSIZE_NATURE_P(JPATCH),I%NR_NATURE_P(:,JPATCH))
!
ENDDO PATCH_LOOP
!
! --------------------------------------------------------------------------------------
! SFX - RRM coupling update if used :
! --------------------------------------------------------------------------------------
!
IF(I%LCPL_RRM)THEN
  CALL DIAG_CPL_ESM_ISBA(PTSTEP,ZCPL_DRAIN,ZCPL_RUNOFF,ZCPL_EFLOOD, &
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
IF ((I%CPHOTO=='NON' .OR. I%CPHOTO=='AGS' .OR. I%CPHOTO=='AST') .AND. I%LVEGUPD) THEN
     CALL VEGETATION_UPDATE(PTSTEP,I%TTIME,I%XCOVER, I%LCOVER,                 &
                         I%CISBA,I%LECOCLIMAP,I%CPHOTO,LAGRIP,I%LTR_ML,'NAT',    &
                         I%XLAI,I%XVEG,I%XZ0,                                  &
                         I%XALBNIR,I%XALBVIS,I%XALBUV,I%XEMIS,                   &
                         I%XRSMIN,I%XGAMMA,I%XWRMAX_CF,                        &
                         I%XRGL,I%XCV,                                       &
                         I%XGMES,I%XBSLAI,I%XLAIMIN,I%XSEFOLD,I%XGC,I%XDMAX,         &
                         I%XF2I, I%LSTRESS,                                  &
                         I%XAOSIP,I%XAOSIM,I%XAOSJP,I%XAOSJM,                    &
                         I%XHO2IP,I%XHO2IM,I%XHO2JP,I%XHO2JM,                    &
                         I%XZ0EFFIP,I%XZ0EFFIM,I%XZ0EFFJP,I%XZ0EFFJM,            &
                         I%CALBEDO, I%XALBNIR_VEG, I%XALBVIS_VEG, I%XALBUV_VEG,  &
                         I%XALBNIR_SOIL, I%XALBVIS_SOIL, I%XALBUV_SOIL,        &
                         I%XCE_NITRO, I%XCF_NITRO, I%XCNA_NITRO,               &
                         I%TSEED, I%TREAP, I%XWATSUP, I%XIRRIG,                  &
                         I%XGNDLITTER,I%XZF_TALLVEG, I%XRGLGV,I%XGAMMAGV,        &
                         I%XRSMINGV, I%XWRMAX_CFGV,                          &
                         I%XH_VEG, I%XLAIGV, I%XZ0LITTER, LUPDATED             )  
!
ELSEIF ((I%CPHOTO=='LAI'.OR.I%CPHOTO=='LST'.OR.I%CPHOTO=='NIT'.OR.I%CPHOTO=='NCB').AND.I%LVEGUPD) THEN
!
  CALL ALBEDO_VEG_UPDATE(PTSTEP,I%TTIME,I%XCOVER, I%LCOVER,                    &
                         I%CISBA,I%LECOCLIMAP,I%CPHOTO,LAGRIP,I%LTR_ML,'NAT',    &
                         I%XVEG,I%XALBNIR,I%XALBVIS,I%XALBUV,                    &
                         I%CALBEDO, I%XALBNIR_VEG, I%XALBVIS_VEG, I%XALBUV_VEG,  &
                         I%XALBNIR_SOIL, I%XALBVIS_SOIL, I%XALBUV_SOIL         )
END IF
!
IF(I%LPERTSURF.AND.LUPDATED) THEN
  ! random perturbation for ensembles:
  ! reset these fields to their original values, as in compute_isba_parameters
  I%XVEG(:,1) = I%XPERTVEG(:)
  I%XLAI(:,1) = I%XPERTLAI(:)
  I%XCV(:,1)  = I%XPERTCV(:)
  ! reapply original perturbation patterns
  WHERE(I%XALBNIR(:,1)/=XUNDEF)  I%XALBNIR(:,1) =I%XALBNIR(:,1) *( 1.+ I%XPERTALB(:) )
  WHERE(I%XALBVIS(:,1)/=XUNDEF)  I%XALBVIS(:,1) =I%XALBVIS(:,1) *( 1.+ I%XPERTALB(:) )
  WHERE(I%XALBUV(:,1)/=XUNDEF)   I%XALBUV(:,1)  =I%XALBUV(:,1)  *( 1.+ I%XPERTALB(:) )
  WHERE(I%XZ0(:,1)/=XUNDEF)      I%XZ0(:,1)     =I%XZ0(:,1)     *( 1.+ I%XPERTZ0(:) )
  WHERE(I%XZ0EFFIP(:,1)/=XUNDEF) I%XZ0EFFIP(:,1)=I%XZ0EFFIP(:,1)*( 1.+ I%XPERTZ0(:) )
  WHERE(I%XZ0EFFIM(:,1)/=XUNDEF) I%XZ0EFFIM(:,1)=I%XZ0EFFIM(:,1)*( 1.+ I%XPERTZ0(:) )
  WHERE(I%XZ0EFFJP(:,1)/=XUNDEF) I%XZ0EFFJP(:,1)=I%XZ0EFFJP(:,1)*( 1.+ I%XPERTZ0(:) )
  WHERE(I%XZ0EFFJM(:,1)/=XUNDEF) I%XZ0EFFJM(:,1)=I%XZ0EFFJM(:,1)*( 1.+ I%XPERTZ0(:) )

ENDIF
!
! --------------------------------------------------------------------------------------
! Outputs for the atmospheric model or update the snow/flood fraction at time t+1
! --------------------------------------------------------------------------------------
! Grid box average fluxes/properties: Arguments and standard diagnostics at time t+1
!
 CALL AVERAGE_FLUX(I%XPATCH,                                             &
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
 CALL AVERAGE_PHY(I%XPATCH,                                         &
                  ZTSURF_TILE, ZZ0_TILE, ZZ0H_TILE, ZQSURF_TILE,  &    
                  PUREF, PZREF, PTSURF, PZ0, PZ0H, PQSURF         )
!
!-------------------------------------------------------------------------------------
!Radiative properties at time t+1 (see by the atmosphere) in order to close
!the energy budget between surfex and the atmosphere
!-------------------------------------------------------------------------------------
!
 CALL UPDATE_RAD_ISBA_n(I%LFLOOD, I%TSNOW%SCHEME, PZENITH2, PSW_BANDS,      &
                       I%XVEG, I%XLAI, I%XZ0,                                 &
                       I%LMEB_PATCH,I%XLAIGV,I%XGNDLITTER,I%XZ0LITTER,I%XH_VEG,   &
                       I%XALBNIR, I%XALBVIS, I%XALBUV, I%XEMIS,                 &
                       ZDIR_ALB_TILE,ZSCA_ALB_TILE,ZEMIS_TILE,          &
                       PDIR_SW, PSCA_SW,                                &
                       I%XZF_TALLVEG,                                     &
                       I%XALBNIR_VEG, I%XALBNIR_SOIL,                       &
                       I%XALBVIS_VEG, I%XALBVIS_SOIL                        )
!
 CALL AVERAGE_RAD(I%XPATCH,                                              &
                 ZDIR_ALB_TILE, ZSCA_ALB_TILE, ZEMIS_TILE, ZTRAD_TILE, &
                 PDIR_ALB,      PSCA_ALB,      I%XEMIS_NAT,  I%XTSRAD_NAT  )  
!
PEMIS = I%XEMIS_NAT
PTRAD = I%XTSRAD_NAT
!
!-------------------------------------------------------------------------------------
!
! Any additional diagnostics (stored in MODD_DIAG_ISBA_n)
!
 CALL AVERAGE_DIAG_ISBA_n(PUREF,PZREF,PSFCO2,PTRAD)
!
! Cumulated diagnostics (stored in MODD_DIAG_EVAP_ISBA_n)
!
 CALL AVERAGE_DIAG_EVAP_ISBA_n(PTSTEP,PRAIN,PSNOW)
!
! Miscellaneous diagnostics (stored in MODD_DIAG_MISC_ISBA_n)
!
 CALL AVERAGE_DIAG_MISC_ISBA_n
!
!--------------------------------------------------------------------------------------
!
 CALL COUPLING_SURF_TOPD(HPROGRAM,U%NDIM_FULL)
!
! --------------------------------------------------------------------------------------
! Snow/Flood fractions, albedo and emissivity update :
! --------------------------------------------------------------------------------------
!
! --------------------------------------------------------------------------------------
! Chemical fluxes :
! --------------------------------------------------------------------------------------
!
IF (CHI%NBEQ>0 .AND. CHI%LCH_BIO_FLUX) THEN
 CALL CH_BVOCEM_n(ZSW_FORBIO,PRHOA,PSFTS)
ENDIF
!
!SOILNOX
IF (CHI%LCH_NO_FLUX) THEN
  CALL SOILEMISNO_n(PU,PV)
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
REAL, DIMENSION(KSIZE,I%NNBIOMASS) :: ZP_RESP_BIOMASS_INST         ! instantaneous biomass respiration (kgCO2/kgair m/s)
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
IF (I%NPATCH==1) THEN
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
   GMEB=I%LMEB_PATCH(JPATCH)
!
! Pack ISBA input and prognostic variables (modd_isban) for each patch:
!
 CALL PACK_ISBA_PATCH_GET_SIZE_n(JPATCH)
!
 CALL PACK_DIAG_PATCH_GET_SIZE_n(JPATCH)
!
 CALL PACK_ISBA_PATCH_n(KMASK,KSIZE,JPATCH)     
!
! Pack chemistry input and prognostic variables (modd_ch_isban) for each patch:
!
IF (CHI%NBEQ>0) THEN
  IF( CHI%CCH_DRY_DEP == "WES89") THEN
    CALL PACK_CH_ISBA_PATCH_n(KMASK,KSIZE,I%NPATCH,JPATCH)     
  END IF
END IF
!
! Allocate ISBA diagnostics for each patch:
!
 CALL PACK_DIAG_PATCH_n(KSIZE,ISWB,JPATCH)     
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Cosine of the slope typically encoutered in the grid mesh (including subgrid orography)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
ZP_SLOPE_COS(:) = 1./SQRT(1.+PKI%XP_SSO_SLOPE(:)**2)
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
IF(I%LFLOOD)THEN
  CALL ISBA_FLOOD_PROPERTIES(PKI%XP_LAI,PKI%XP_FFLOOD,PKI%XP_FFROZEN,  &
                             ZP_Z0FLOOD,ZP_FFGNOS,ZP_FFVNOS)  
ELSE
  ZP_Z0FLOOD = XUNDEF
  ZP_FFGNOS  = 0.0
  ZP_FFVNOS  = 0.0
ENDIF
!
! For multi-energy balance
   IF(GMEB)THEN
     ZSNOWDEPTH(:) = SUM(PKI%XP_SNOWSWE(:,:)/PKI%XP_SNOWRHO(:,:),2)
     ZPALPHAN(:)=MEBPALPHAN(ZSNOWDEPTH,PKI%XP_H_VEG)
   ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Surface Roughness lengths (m):
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!* effective roughness
!
 CALL Z0EFF(I%CROUGH, GMEB, ZP_ALFA, ZP_ZREF, ZP_UREF, PKI%XP_Z0, PKI%XP_Z0REL, PKI%XP_PSN,   &
     ZPALPHAN,PKI%XP_Z0LITTER, PKI%XP_SNOWSWE(:,1),                              &
     PKI%XP_Z0EFFIP,PKI%XP_Z0EFFIM,PKI%XP_Z0EFFJP,PKI%XP_Z0EFFJM, PKI%XP_FF, ZP_Z0FLOOD,     &
     PKI%XP_AOSIP,PKI%XP_AOSIM,PKI%XP_AOSJP,PKI%XP_AOSJM,                                &
     PKI%XP_HO2IP,PKI%XP_HO2IM,PKI%XP_HO2JP,PKI%XP_HO2JM,                                &
     PKI%XP_Z0_O_Z0H, PKDI%XP_Z0_WITH_SNOW, PKDI%XP_Z0H_WITH_SNOW, PKDI%XP_Z0EFF,           &
     ZZ0G_WITHOUT_SNOW,                                                  &
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
 CALL ISBA_ALBEDO(I%TSNOW%SCHEME, I%LTR_ML, GMEB,                            &
                   ZP_DIR_SW, ZP_SCA_SW, PSW_BANDS,ISWB,                 &
                   PKI%XP_ALBNIR, PKI%XP_ALBVIS, PKI%XP_ALBUV,                       &
                   PKI%XP_ALBNIR_VEG, PKI%XP_ALBVIS_VEG, PKI%XP_ALBUV_VEG,           &
                   PKI%XP_ALBNIR_SOIL, PKI%XP_ALBVIS_SOIL, PKI%XP_ALBUV_SOIL,        &
                   PKI%XP_SNOWALB, PKI%XP_PSNV, PKI%XP_PSNG, PKI%XP_ALBF, PKI%XP_FFV, PKI%XP_FFG,& 
                   ZP_GLOBAL_SW, PKDI%XP_SNOWFREE_ALB, PKDI%XP_SNOWFREE_ALB_VEG,   &
                   PKDI%XP_SNOWFREE_ALB_SOIL, ZP_MEB_SCA_SW,                  &
                   ZP_ALBNIR_TVEG, ZP_ALBVIS_TVEG,                       &
                   ZP_ALBNIR_TSOIL, ZP_ALBVIS_TSOIL                      )  
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Intialize computation of ISBA water and energy budget
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
 CALL ISBA_BUDGET_INIT(I%CISBA,I%TSNOW%SCHEME,            &
                      PKI%XP_WG,PKI%XP_WGI,PKI%XP_WR,PKI%XP_SNOWSWE, &
                      PKI%XP_DG, PKI%XP_DZG, ZP_WG_INI,      &
                      ZP_WGI_INI, ZP_WR_INI,         &
                      ZP_SWE_INI                     )
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Over Natural Land Surfaces:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ZIRRIG_GR(:)= 0.
!
 CALL ISBA(I%CISBA, I%CPHOTO, I%LTR_ML, I%CRUNOFF, I%CKSAT, I%CRAIN, I%CHORT, I%CC1DRY, I%CSCOND,  &
 I%TSNOW%SCHEME, I%CSNOWRES, I%CCPSURF, I%CSOILFRZ, I%CDIFSFCOND, I%TTIME, I%LFLOOD, I%LTEMP_ARP,&
 I%LGLACIER, GMEB, I%LFORC_MEASURE, PTSTEP, CIMPLICIT_WIND, I%LAGRI_TO_GRASS,          &
 I%LSNOWDRIFT, I%LSNOWDRIFT_SUBLIM, I%LSNOW_ABS_ZENITH, I%CSNOWMETAMO, I%CSNOWRAD,         &
 I%XCGMAX, ZP_ZREF, ZP_UREF, ZP_SLOPE_COS, ZP_TA, ZP_QA, ZP_EXNA,                  &
 ZP_RHOA, ZP_PS, ZP_EXNS, ZP_RAIN, ZP_SNOW, ZP_ZENITH, ZP_MEB_SCA_SW,            &
 ZP_GLOBAL_SW, ZP_LW,                                                            &
 ZP_WIND, ZP_PEW_A_COEF, ZP_PEW_B_COEF, ZP_PET_A_COEF, ZP_PEQ_A_COEF,            &
 ZP_PET_B_COEF, ZP_PEQ_B_COEF,  PKI%XP_RSMIN, PKI%XP_RGL, PKI%XP_GAMMA, PKI%XP_CV, PKI%XP_RUNOFFD,   &
 PKI%XP_SOILWGHT, I%NLAYER_HORT, I%NLAYER_DUN, ZP_ALBNIR_TVEG, ZP_ALBVIS_TVEG,           &
 ZP_ALBNIR_TSOIL, ZP_ALBVIS_TSOIL, PKDI%XP_SNOWFREE_ALB, PKI%XP_WRMAX_CF, PKI%XP_VEG, PKI%XP_LAI, &
 PKI%XP_EMIS, PKDI%XP_Z0_WITH_SNOW, PKDI%XP_Z0H_WITH_SNOW, PKI%XP_VEGTYPE_PATCH, PKDI%XP_Z0EFF,         &
 PKI%XP_ZF_TALLVEG , PKI%XP_RGLV, PKI%XP_GAMMAV, PKI%XP_RSMINV, PKI%XP_ROOTFRACV, PKI%XP_WRMAX_CFV,      &
 PKI%XP_LAIV, PKI%XP_BSLAI,PKI%XP_LAIMIN,PKI%XP_H_VEG,ZPALPHAN, ZZ0G_WITHOUT_SNOW, ZZ0_MEBV,     &
 ZZ0H_MEBV,ZZ0EFF_MEBV, ZZ0_MEBN,ZZ0H_MEBN,ZZ0EFF_MEBN, PKI%XP_GNDLITTER,            &
 PKI%XP_RUNOFFB, PKI%XP_CGSAT, PKI%XP_C1SAT, PKI%XP_C2REF, PKI%XP_C3, PKI%XP_C4B, PKI%XP_C4REF, PKI%XP_ACOEF,    &
 PKI%XP_PCOEF, PKI%XP_TAUICE, PKI%XP_WDRAIN, ZP_TDEEP_A, PKI%XP_TDEEP, PKI%XP_GAMMAT,                &
 PKI%XP_PSN, PKI%XP_PSNG, PKI%XP_PSNV,                                                       &
 PKI%XP_PSNV_A, PKDI%XP_SNOWFREE_ALB_VEG, PKDI%XP_SNOWFREE_ALB_SOIL, PKI%XP_IRRIG, PKI%XP_WATSUP,      &
 PKI%XP_THRESHOLD, PKI%XP_LIRRIGATE, PKI%XP_LIRRIDAY, PKI%LP_STRESS, PKI%XP_GC, PKI%XP_F2I, PKI%XP_DMAX,     &
 PKI%XP_AH, PKI%XP_BH, ZP_CO2, PKI%XP_GMES, I%XPOI, PKI%XP_FZERO, PKI%XP_EPSO, PKI%XP_GAMM, PKI%XP_QDGAMM,     &
 PKI%XP_QDGMES, PKI%XP_T1GMES, PKI%XP_T2GMES, PKI%XP_AMAX, PKI%XP_QDAMAX,  PKI%XP_T1AMAX, PKI%XP_T2AMAX,     &
 I%XABC, PKI%XP_DG, PKI%XP_DZG, PKI%XP_DZDIF, PKI%NK_WG_LAYER, PKI%XP_ROOTFRAC, PKI%XP_WFC,                &
 PKI%XP_WWILT, PKI%XP_WSAT, PKI%XP_BCOEF, PKI%XP_CONDSAT, PKI%XP_MPOTSAT, PKI%XP_HCAPSOIL, PKI%XP_CONDDRY,   &
 PKI%XP_CONDSLD, PKI%XP_D_ICE, PKI%XP_KSAT_ICE, PKI%XP_MUF, PKI%XP_FF, PKI%XP_FFG, PKI%XP_FFV, ZP_FFGNOS,    &
 ZP_FFVNOS, PKI%XP_FFROZEN, PKI%XP_ALBF, PKI%XP_EMISF, PKI%XP_FFLOOD, PKI%XP_PIFLOOD, PKDI%XP_IFLOOD,     &
 PKDI%XP_PFLOOD, PKDI%XP_LE_FLOOD, PKDI%XP_LEI_FLOOD, I%XSODELX, PKI%XP_LAT, PKI%XP_LON, PKI%XP_TG, PKI%XP_WG,    &
 PKI%XP_WGI, PKI%XP_CPS, PKI%XP_LVTT, PKI%XP_LSTT, PKI%XP_WR, PKI%XP_WRV,PKI%XP_WRVN,PKI%XP_TV,                  &
 PKI%XP_RESA, PKI%XP_ANFM, PKI%XP_FSAT, PKI%XP_SNOWALB, PKI%XP_SNOWALBVIS, PKI%XP_SNOWALBNIR,            &
 PKI%XP_SNOWALBFIR, PKI%XP_SNOWSWE, PKI%XP_SNOWHEAT, PKI%XP_SNOWRHO, PKI%XP_SNOWGRAN1, PKI%XP_SNOWGRAN2, &
 PKI%XP_SNOWHIST, PKI%XP_SNOWAGE, PKDI%XP_GRNDFLUX, PKDI%XP_HPSNOW, PKDI%XP_SNOWHMASS,                  &
 PKDI%XP_RNSNOW, PKDI%XP_HSNOW, PKDI%XP_GFLUXSNOW, PKDI%XP_USTARSNOW, PKDI%XP_SRSFC, PKDI%XP_RRSFC, PKDI%XP_LESL,   &
 PKI%XP_SNOWEMIS, PKDI%XP_CDSNOW, PKDI%XP_CHSNOW, PKDI%XP_TSRAD, PKDI%XP_TS, PKDI%XP_HV, PKDI%XP_QS, PKDI%XP_SNOWTEMP,  &
 PKDI%XP_SNOWLIQ, PKDI%XP_SNOWDZ, PKDI%XP_CG, PKDI%XP_C1, PKDI%XP_C2, PKDI%XP_WGEQ, PKDI%XP_CT, PKDI%XP_CH, PKDI%XP_CD,       &
 PKDI%XP_CDN, PKDI%XP_RI, PKDI%XP_HU, PKDI%XP_HUG, ZP_EMIS, PKDI%XP_ALBT, PKDI%XP_RS, PKI%XP_LE, PKDI%XP_RN, PKDI%XP_H,      &
 PKDI%XP_LEI, PKDI%XP_LEGI, PKDI%XP_LEG, PKDI%XP_LEV, PKDI%XP_LES, PKDI%XP_LER, PKDI%XP_LETR, PKDI%XP_EVAP, PKDI%XP_GFLUX,    &
 PKDI%XP_RESTORE, ZP_USTAR, PKDI%XP_DRAIN, PKDI%XP_RUNOFF, PKDI%XP_MELT, PKDI%XP_MELTADV,                 &
 PKI%XP_TC,PKI%XP_QC, PKDI%XP_RN_ISBA,                                                        &
 PKDI%XP_H_ISBA, PKDI%XP_LEG_ISBA, PKDI%XP_LEGI_ISBA, PKDI%XP_LEV_ISBA, PKDI%XP_LETR_ISBA, PKDI%XP_USTAR_ISBA, &
 PKDI%XP_LER_ISBA, PKDI%XP_LE_ISBA, PKDI%XP_LEI_ISBA, PKDI%XP_GFLUX_ISBA, PKDI%XP_HORT, PKDI%XP_DRIP, PKDI%XP_RRVEG,&
 ZP_AC_AGG, ZP_HU_AGG, PKI%XP_FAPARC, PKI%XP_FAPIRC, PKI%XP_MUS, PKI%XP_LAI_EFFC, PKI%XP_AN,         &
 PKI%XP_ANDAY, ZP_RESP_BIOMASS_INST, PKDI%XP_IACAN, PKI%XP_ANF, PKDI%XP_GPP, PKDI%XP_FAPAR, PKDI%XP_FAPIR,   &
 PKDI%XP_FAPAR_BS, PKDI%XP_FAPIR_BS, PKDI%XP_IRRIG_FLUX, ZP_DEEP_FLUX,                          &
 PKDI%XP_SWNET_V, PKDI%XP_SWNET_G, PKDI%XP_SWNET_N, PKDI%XP_SWNET_NS, PKDI%XP_LWNET_V, PKDI%XP_LWNET_G,        &
 PKDI%XP_LWNET_N, PKDI%XP_LEVCV, PKDI%XP_LESC, PKDI%XP_H_V_C, PKDI%XP_H_G_C, PKDI%XP_LETRGV, PKDI%XP_LETRCV,        &
 PKDI%XP_LERGV, PKDI%XP_LERCV, PKDI%XP_H_C_A, PKDI%XP_H_N_C, PKDI%XP_LE_C_A, PKDI%XP_LE_V_C, PKDI%XP_LE_G_C,        &
 PKDI%XP_LE_N_C, PKDI%XP_EVAP_N_C, PKDI%XP_EVAP_G_C, PKDI%XP_SR_GN, PKDI%XP_MELTCV, PKDI%XP_FRZCV,             &
 PKDI%XP_SWDOWN_GN, PKDI%XP_LWDOWN_GN,                                                     &
 ZIRRIG_GR, PKI%XP_TOPQS, PKDI%XP_QSB, PKDI%XP_SUBL, PKI%XP_FWTD, PKI%XP_WTD, PKDI%XP_SNDRIFT               )
!
ZP_TRAD=PKDI%XP_TSRAD
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Glacier : ice runoff flux (especally for Earth System Model)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
IF(I%LGLACIER)THEN
!           
  CALL HYDRO_GLACIER(PTSTEP,ZP_SNOW,PKI%XP_SNOWRHO,PKI%XP_SNOWSWE,PKI%XP_ICE_STO,PKDI%XP_ICEFLUX)
!     
ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Calculation of ISBA water and energy budget (and time tendencies of each reservoir)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
CALL ISBA_BUDGET(I%CISBA,I%TSNOW%SCHEME,I%LGLACIER,PTSTEP,          &
                 PKI%XP_WG,PKI%XP_WGI,PKI%XP_WR,PKI%XP_SNOWSWE,PKI%XP_DG,PKI%XP_DZG,  & 
                 ZP_WG_INI,ZP_WGI_INI,ZP_WR_INI,ZP_SWE_INI,   &
                 ZP_RAIN,ZP_SNOW,PKDI%XP_EVAP,PKDI%XP_DRAIN,PKDI%XP_RUNOFF,  &
                 PKDI%XP_IFLOOD,PKDI%XP_PFLOOD,PKDI%XP_ICEFLUX,PKDI%XP_IRRIG_FLUX,&
                 PKDI%XP_SNDRIFT,                                  &
                 PKDI%XP_DWG,PKDI%XP_DWGI,PKDI%XP_DWR,PKDI%XP_DSWE,PKDI%XP_WATBUD      )
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Evolution of soil albedo, when depending on surface soil wetness:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
IF (I%CALBEDO=='EVOL' .AND. I%LECOCLIMAP) THEN
  CALL SOIL_ALBEDO(I%CALBEDO,                                    &
                   PKI%XP_WSAT(:,1),PKI%XP_WG(:,1),                    &
                   PKI%XP_ALBVIS_DRY,PKI%XP_ALBNIR_DRY,PKI%XP_ALBUV_DRY,   &
                   PKI%XP_ALBVIS_WET,PKI%XP_ALBNIR_WET,PKI%XP_ALBUV_WET,   &
                   PKI%XP_ALBVIS_SOIL,PKI%XP_ALBNIR_SOIL,PKI%XP_ALBUV_SOIL )  
  !
  CALL ALBEDO(I%CALBEDO,                                          &
              PKI%XP_ALBVIS_VEG,PKI%XP_ALBNIR_VEG,PKI%XP_ALBUV_VEG,PKI%XP_VEG,  &
              PKI%XP_ALBVIS_SOIL,PKI%XP_ALBNIR_SOIL,PKI%XP_ALBUV_SOIL,      &
              PKI%XP_ALBVIS,PKI%XP_ALBNIR,PKI%XP_ALBUV                      )  
END IF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Vegetation evolution for interactive LAI
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
IF (I%CPHOTO=='LAI' .OR. I%CPHOTO=='LST' .OR. I%CPHOTO=='NIT' .OR. I%CPHOTO=='NCB') THEN
  CALL VEGETATION_EVOL(I%CISBA, I%CPHOTO, I%CRESPSL, I%CALBEDO, LAGRIP, I%LTR_ML,           &
                       I%LNITRO_DILU, I%LAGRI_TO_GRASS,                               &
                       PTSTEP, KMONTH, KDAY, I%NSPINW, PTIME, PKI%XP_LAT, ZP_RHOA,      &
                       PKI%XP_DG, PKI%XP_DZG, PKI%NK_WG_LAYER,                                &                       
                       PKI%XP_TG, PKI%XP_ALBNIR_VEG, PKI%XP_ALBVIS_VEG, PKI%XP_ALBUV_VEG,         &
                       PKI%XP_ALBNIR_SOIL, PKI%XP_ALBVIS_SOIL, PKI%XP_ALBUV_SOIL,             &
                       PKI%XP_VEGTYPE_PATCH, PKI%XP_SEFOLD, PKI%XP_ANMAX, PKI%XP_H_TREE, PKI%XP_BSLAI,&
                       PKI%XP_LAIMIN, ZP_CO2, PKI%XP_CE_NITRO, PKI%XP_CF_NITRO, PKI%XP_CNA_NITRO, &
                       PKI%XP_BSLAI_NITRO, PKI%XP_GMES, PKI%XP_TAU_WOOD, PKI%TP_SEED,             &
                       PKI%TP_REAP, PKI%XP_AOSIP, PKI%XP_AOSIM, PKI%XP_AOSJP, PKI%XP_AOSJM,           &
                       PKI%XP_HO2IP, PKI%XP_HO2IM, PKI%XP_HO2JP, PKI%XP_HO2JM, PKI%XP_Z0EFFIP,        &
                       PKI%XP_Z0EFFIM, PKI%XP_Z0EFFJP, PKI%XP_Z0EFFJM, PKI%XP_LAI, PKI%XP_VEG,        &
                       PKI%XP_Z0, PKI%XP_ALBNIR, PKI%XP_ALBVIS, PKI%XP_ALBUV, PKI%XP_EMIS,            &
                       PKI%XP_ANFM, PKI%XP_ANDAY, PKI%XP_BIOMASS, PKI%XP_RESP_BIOMASS,            &
                       ZP_RESP_BIOMASS_INST, PKI%XP_INCREASE, PKI%XP_TURNOVER             )  
END IF
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Diagnostic of respiration carbon fluxes and soil carbon evolution
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
ZP_SFCO2    (:)=0.
PKDI%XP_RESP_ECO (:)=0.
PKDI%XP_RESP_AUTO(:)=0.
!
IF ( I%CPHOTO/='NON' .AND. I%CRESPSL/='NON' .AND. ANY(PKI%XP_LAI(:)/=XUNDEF) ) THEN
  CALL CARBON_EVOL(I%CISBA, I%CRESPSL, I%CPHOTO, PTSTEP, I%NSPINS,                   &
                   ZP_RHOA, PKI%XP_TG, PKI%XP_WG, PKI%XP_WFC, PKI%XP_WWILT, PKI%XP_WSAT, PKI%XP_SAND,&
                   PKI%XP_DG, PKI%XP_DZG, PKI%NK_WG_LAYER,                               &                   
                   PKI%XP_RE25, PKI%XP_LAI, ZP_RESP_BIOMASS_INST, PKI%XP_TURNOVER,       &
                   PKI%XP_LITTER, PKI%XP_LIGNIN_STRUC , PKI%XP_SOILCARB,                 &
                   PKDI%XP_RESP_AUTO, PKDI%XP_RESP_ECO                                 )  
  ! calculation of vegetation CO2 flux
  ! Positive toward the atmosphere
  ZP_SFCO2(:) = PKDI%XP_RESP_ECO(:) - PKDI%XP_GPP(:)  
END IF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Reset effecitve roughness lentgh to its nominal value when snow has just disappeared
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
 CALL SUBSCALE_Z0EFF(PKI%XP_AOSIP,PKI%XP_AOSIM,PKI%XP_AOSJP,PKI%XP_AOSJM,            &
                    PKI%XP_HO2IP,PKI%XP_HO2IM,PKI%XP_HO2JP,PKI%XP_HO2JM,PKI%XP_Z0,      &
                    PKI%XP_Z0EFFIP,PKI%XP_Z0EFFIM,PKI%XP_Z0EFFJP,PKI%XP_Z0EFFJM,    &
                    OMASK=(PKI%XP_SNOWSWE(:,1)==0. .AND. PKI%XP_PSN(:)>0.)  )   
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Turbulent fluxes
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
ZP_SFTH(:) = PKDI%XP_H(:)
ZP_SFTQ(:) = PKDI%XP_EVAP(:)

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
IF (CHI%NBEQ>0) THEN
  IF( CHI%CCH_DRY_DEP == "WES89") THEN

    CALL CH_DEP_ISBA         (ZP_USTAR, PKDI%XP_HU, PKI%XP_PSN,             &
                        PKI%XP_VEG, PKI%XP_LAI, PKI%XP_SAND, PKI%XP_CLAY, PKI%XP_RESA, &
                        PKDI%XP_RS(:),  PKI%XP_Z0(:),                       &
                        ZP_TA, ZP_PA, ZP_TRAD(:),                  &
                        PKI%XP_VEGTYPE_PATCH(:,NVT_NO),                &
                        PKI%XP_VEGTYPE_PATCH(:,NVT_ROCK),              &
                        CHI%CSV(CHI%NSV_CHSBEG:CHI%NSV_CHSEND),                &
                        PKCI%XP_SOILRC_SO2,  PKCI%XP_SOILRC_O3 ,             &
                        PKCI%XP_DEP(:,1:CHI%NBEQ)                           )  
 
    ZP_SFTS(:,CHI%NSV_CHSBEG:CHI%NSV_CHSEND) = - ZP_SV(:,CHI%NSV_CHSBEG:CHI%NSV_CHSEND)  &
                                                    * PKCI%XP_DEP(:,1:CHI%NBEQ)  
    IF (CHI%NAEREQ > 0 ) THEN
      CALL CH_AER_DEP(ZP_SV(:,CHI%NSV_AERBEG:CHI%NSV_AEREND),&
                           ZP_SFTS(:,CHI%NSV_AERBEG:CHI%NSV_AEREND),&
                           ZP_USTAR, PKI%XP_RESA,ZP_TA,ZP_RHOA)     
    END IF
  ELSE
    ZP_SFTS(:,CHI%NSV_CHSBEG:CHI%NSV_CHSEND) = 0.
    ZP_SFTS(:,CHI%NSV_AERBEG:CHI%NSV_AEREND) = 0.
  ENDIF
ENDIF
!
! --------------------------------------------------------------------------------------
! Dust deposition and emission:
! --------------------------------------------------------------------------------------
!
IF(CHI%NDSTEQ>0)THEN
  IDST = CHI%NSV_DSTEND - CHI%NSV_DSTBEG + 1

  CALL COUPLING_DST_n(                   &
            HPROGRAM,                    &!I [char] Name of program
            KSIZE,      &!I [nbr] number of points in patch
            IDST,                        &!I [nbr] number of dust emissions variables
            JPATCH,                      &!I [idx] patch in question
            PKI%XP_CLAY(:,1),                &!I [frc] mass fraction clay in first soil layer
            ZP_PS,                       &!I [Pa] surface pressure
            PKDI%XP_TS,                       &!I [K] surface temperature
            ZP_QA,                       &!I [kg/kg] specific humidity
            PKI%XP_RESA,                     &!I [s/m] atmospheric resistance
            ZP_RHOA,                     &!I [kg/m3] atmospheric density
            PKI%XP_SAND(:,1),                &!I [frc] mass fraction of sand in first soil layer
            ZP_PA,                       &!I [K] Atmospheric pressure
            ZP_TA,                       &!I [K] Atmospheric temperature
            PKI%XP_TG(:,1),                  &!I [K] Ground temperature
            ZP_U,                        &!I [m/s] zonal wind at atmospheric height 
            ZP_UREF,                     &!I [m] reference height of wind
            ZP_V,                        &!I [m/s] meridional wind at atmospheric height
            PKI%XP_WG(:,1),                  &!I [m3/m3] ground volumetric water content
            PKI%XP_WSAT(:,1),                &!I [m3/m3] saturation volumetric water content
            ZP_ZREF,                     &!I [m] reference height of wind
            PKDI%XP_CD,                       &!I [] Drag Coefficient for momentum
            PKDI%XP_CDN,                      &!I [] Drag neutral Coefficient for momentum
            PKDI%XP_CH,                       &!I [] drag coefficient for heat
            PKDI%XP_RI,                       &!I [] Richardson number
            PKDI%XP_Z0H_WITH_SNOW,            &!I [frc] Z0 (heat) with snow
            ZP_SFTS(:,CHI%NSV_DSTBEG:CHI%NSV_DSTEND)  &!O [kg/m2/sec] flux of dust            
            )  
!
   IF (CHI%NSV_AEREND > 0)  THEN ! case of dust/ anthropogenic aerosols coupling
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
           ZP_SFTS(:,JSV) = ZP_SFTS(:,JSV) + ZP_SFTS(:,CHI%NSV_DSTBEG-1+JSV_IDX+2)*XAVOGADRO/XMOLARWEIGHT_DST
         END IF
         IF ( (TRIM(HSV(JSV)) == "@DSTJ").AND.(JMODE==2)) THEN 
           ! add dust flux and conversion kg/m2/sec into molec.m2/s
           ZP_SFTS(:,JSV) = ZP_SFTS(:,JSV) + ZP_SFTS(:,CHI%NSV_DSTBEG-1+JSV_IDX+2)*XAVOGADRO/XMOLARWEIGHT_DST
         END IF
       END DO
       !
     END DO
    END IF
!    
!Modify fluxes due to dry deposition, we introduce a negative flux where dust is lost
  CALL DSLT_DEP(ZP_SV(:,CHI%NSV_DSTBEG:CHI%NSV_DSTEND), ZP_SFTS(:,CHI%NSV_DSTBEG:CHI%NSV_DSTEND), &
                ZP_USTAR, PKI%XP_RESA, ZP_TA, ZP_RHOA, DST%XEMISSIG_DST, DST%XEMISRADIUS_DST, &
                JPMODE_DST, XDENSITY_DST, XMOLARWEIGHT_DST, ZCONVERTFACM0_DST,    &
                ZCONVERTFACM6_DST, ZCONVERTFACM3_DST, LVARSIG_DST, LRGFIX_DST,    &
                CVERMOD  )
!
!Transfer these fluxes to fluxes understandable by all moments
  CALL MASSFLUX2MOMENTFLUX(           &
    ZP_SFTS(:,CHI%NSV_DSTBEG:CHI%NSV_DSTEND), & !I/O ![kg/m2/sec] In: flux of only mass, out: flux of moments
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
IF (CHI%NSLTEQ>0) THEN
  CALL DSLT_DEP(ZP_SV(:,CHI%NSV_SLTBEG:CHI%NSV_SLTEND), ZP_SFTS(:,CHI%NSV_SLTBEG:CHI%NSV_SLTEND), &
                ZP_USTAR, PKI%XP_RESA, ZP_TA, ZP_RHOA, SLT%XEMISSIG_SLT, SLT%XEMISRADIUS_SLT, &
                JPMODE_SLT, XDENSITY_SLT, XMOLARWEIGHT_SLT, ZCONVERTFACM0_SLT,    &
                ZCONVERTFACM6_SLT, ZCONVERTFACM3_SLT, LVARSIG_SLT, LRGFIX_SLT,    &
                CVERMOD  )  

  CALL MASSFLUX2MOMENTFLUX(           &
    ZP_SFTS(:,CHI%NSV_SLTBEG:CHI%NSV_SLTEND), & !I/O ![kg/m2/sec] In: flux of only mass, out: flux of moments
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
 CALL DIAG_INLINE_ISBA_n(ZP_TA, PKDI%XP_TS, ZP_QA, ZP_PA, ZP_PS, ZP_RHOA, ZP_U, ZP_V,       &
                          ZP_ZREF, ZP_UREF,                                            &
                          PKDI%XP_CD, PKDI%XP_CDN, PKDI%XP_CH, PKDI%XP_RI, PKDI%XP_HU, PKDI%XP_Z0_WITH_SNOW,         &
                          PKDI%XP_Z0H_WITH_SNOW, PKDI%XP_Z0EFF,                                  &
                          ZP_SFTH, ZP_SFTQ, ZP_SFU, ZP_SFV, PKDI%XP_QS,                     &
                          PKI%XP_DIR_ALB_WITH_SNOW, PKI%XP_SCA_ALB_WITH_SNOW,                  &
                          ZP_DIR_SW, ZP_SCA_SW, ZP_LW, PKDI%XP_RN                           )  
!
!
!-------------------------------------------------------------------------------
!Physical properties see by the atmosphere in order to close the energy budget 
!between surfex and the atmosphere. All variables should be at t+1 but very 
!difficult to do. Maybe it will be done later. However, Ts can be at time t+1
!-------------------------------------------------------------------------------
!
ZP_TSURF (:) = PKDI%XP_TS (:)
ZP_Z0    (:) = PKDI%XP_Z0_WITH_SNOW (:)
ZP_Z0H   (:) = PKDI%XP_Z0H_WITH_SNOW(:)
ZP_QSURF (:) = PKDI%XP_QS (:)
!
!-------------------------------------------------------------------------------
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Isba offline diagnostics for each patch
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
 CALL DIAG_EVAP_ISBA_n(I%CPHOTO,PTSTEP,KMASK,KSIZE,JPATCH,ZP_RHOA)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Isba offline diagnostics for miscellaneous terms over each patch
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
 CALL DIAG_MISC_ISBA_n(PTSTEP, I%CISBA, I%CPHOTO, I%TSNOW%SCHEME, LAGRIP, I%LTR_ML,    &
                      PTIME, KSIZE, JPATCH, KMASK, PKI%XP_THRESHOLD,              &
                      PKI%XP_PSN, PKI%XP_PSNG, PKI%XP_PSNV, PKI%XP_FF, PKI%XP_FFG, PKI%XP_FFV,        &
                      PKI%XP_WG, PKI%XP_WGI, PKI%XP_WFC, PKI%XP_WWILT, PKI%XP_SNOWSWE, PKI%XP_SNOWRHO,&
                      PKI%XP_FAPARC, PKI%XP_FAPIRC, PKI%XP_LAI_EFFC, PKI%XP_MUS, PKI%XP_FSAT,     &
                      PKI%XP_DG, PKI%XP_TG       )                  
!
! Unpack ISBA diagnostics (modd_diag_isban) for each patch:ISIZE_MAX = MAXVAL(NSIZE_NATURE_P)

!  (MUST be done BEFORE UNPACK_ISBA_PATCH, because of XP_LE)
!
 CALL UNPACK_DIAG_PATCH_n(KMASK,KSIZE,I%NPATCH,JPATCH, &
                           ZCPL_DRAIN,ZCPL_RUNOFF,ZCPL_EFLOOD,ZCPL_PFLOOD,           &
                           ZCPL_IFLOOD, ZCPL_ICEFLUX)  
!
! for chemical deposition
!
IF (CHI%NBEQ>0) THEN
  IF( CHI%CCH_DRY_DEP == "WES89") THEN
    CALL UNPACK_CH_ISBA_PATCH_n(KMASK,KSIZE,I%NPATCH,JPATCH)     
  END IF
END IF
!
! Unpack ISBA variables (modd_isban) for each patch:
!
 CALL UNPACK_ISBA_PATCH_n(KMASK,KSIZE,JPATCH)
!
!----------------------------------------------------------------------
!
! for further chemical biogenic emissions
!
IF (CHI%NBEQ>0 .AND. CHI%LCH_BIO_FLUX) THEN
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
IF (I%NPATCH==1) THEN
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
IF (CHI%NDSTEQ>0) THEN
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
