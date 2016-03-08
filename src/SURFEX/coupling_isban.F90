!     ###############################################################################
SUBROUTINE COUPLING_ISBA_n (DTCO, UG, U, USS, AG, CHI, DTI, DGI, GB, ISS, IG,           &
                            IO, IP, INP, INI, IPM, IR, PKCI, PKD, PK, DST, SLT,         &
                            HPROGRAM, HCOUPLING, PTSTEP,  KYEAR, KMONTH, KDAY, PTIME,   &
                            KI, KSV, KSW, PTSUN, PZENITH, PZENITH2, PZREF, PUREF, PZS,  &
                            PU, PV, PQA, PTA, PRHOA, PSV, PCO2, HSV, PRAIN, PSNOW, PLW, &
                            PDIR_SW, PSCA_SW, PSW_BANDS, PPS, PPA, PSFTQ, PSFTH, PSFTS, &
                            PSFCO2, PSFU, PSFV, PTRAD, PDIR_ALB, PSCA_ALB, PEMIS,       &
                            PTSURF, PZ0, PZ0H, PQSURF, PPEW_A_COEF, PPEW_B_COEF,        &
                            PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF, HTEST    )  
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
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_AGRI_n, ONLY : AGRI_t
USE MODD_CH_ISBA_n, ONLY : CH_ISBA_t
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_SURFEX_n, ONLY : ISBA_DIAG_t
USE MODD_GR_BIOG_n, ONLY : GR_BIOG_t
USE MODD_SSO_n, ONLY : SSO_t
USE MODD_GRID_n, ONLY : GRID_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_PGD_n, ONLY : ISBA_PGD_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t, ISBA_INIT_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
USE MODD_PACK_CH_ISBA, ONLY : PACK_CH_ISBA_t
USE MODD_PACK_DIAG_ISBA, ONLY : PACK_DIAG_ISBA_t
USE MODD_PACK_ISBA, ONLY : PACK_ISBA_t
!
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
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SSO_t), INTENT(INOUT) :: USS
!
TYPE(AGRI_t), INTENT(INOUT) :: AG
TYPE(CH_ISBA_t), INTENT(INOUT) :: CHI
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTI
TYPE(ISBA_DIAG_t), INTENT(INOUT) :: DGI
TYPE(GR_BIOG_t), INTENT(INOUT) :: GB
TYPE(SSO_t), INTENT(INOUT) :: ISS 
TYPE(GRID_t), INTENT(INOUT) :: IG
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_PGD_t), INTENT(INOUT) :: IP
TYPE(ISBA_PARAM_t), INTENT(INOUT) :: IPM
TYPE(ISBA_INIT_t), INTENT(INOUT) :: INI
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: INP
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
TYPE(PACK_CH_ISBA_t), INTENT(INOUT) :: PKCI
TYPE(PACK_DIAG_ISBA_t), INTENT(INOUT) :: PKD
TYPE(PACK_ISBA_t), INTENT(INOUT) :: PK
!
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
REAL, DIMENSION(KI,IO%NPATCH) :: ZSFTH_TILE     ! surface heat flux (W/m2)
REAL, DIMENSION(KI,IO%NPATCH) :: ZSFTQ_TILE     ! surface vapor flux (kg/m2/s)
REAL, DIMENSION(KI,IO%NPATCH) :: ZSFCO2_TILE    ! surface CO2 flux positive toward the atmosphere (m/s*kg_CO2/kg_air)
REAL, DIMENSION(KI,IO%NPATCH) :: ZSFU_TILE      ! zonal momentum flux
REAL, DIMENSION(KI,IO%NPATCH) :: ZSFV_TILE      ! meridian momentum flux
REAL, DIMENSION(KI,IO%NPATCH) :: ZTRAD_TILE     ! radiative surface temperature
REAL, DIMENSION(KI,IO%NPATCH) :: ZEMIS_TILE     ! emissivity
REAL, DIMENSION(KI,IO%NPATCH) :: ZTSURF_TILE    ! surface effective temperature
REAL, DIMENSION(KI,IO%NPATCH) :: ZZ0_TILE       ! roughness length for momentum
REAL, DIMENSION(KI,IO%NPATCH) :: ZZ0H_TILE      ! roughness length for heat
REAL, DIMENSION(KI,IO%NPATCH) :: ZQSURF_TILE    ! specific humidity at surface
REAL, DIMENSION(KI,KSW,IO%NPATCH) :: ZDIR_ALB_TILE  ! direct albedo
REAL, DIMENSION(KI,KSW,IO%NPATCH) :: ZSCA_ALB_TILE  ! diffuse albedo
REAL, DIMENSION(KI,KSV,IO%NPATCH) :: ZSFTS_TILE     ! scalar surface flux
!
REAL, DIMENSION(KI, IO%NPATCH) :: ZCPL_DRAIN     ! For the coupling with TRIP
REAL, DIMENSION(KI, IO%NPATCH) :: ZCPL_RUNOFF    ! For the coupling with TRIP
REAL, DIMENSION(KI, IO%NPATCH) :: ZCPL_EFLOOD    ! For the coupling with TRIP
REAL, DIMENSION(KI, IO%NPATCH) :: ZCPL_PFLOOD    ! For the coupling with TRIP
REAL, DIMENSION(KI, IO%NPATCH) :: ZCPL_IFLOOD    ! For the coupling with TRIP
REAL, DIMENSION(KI, IO%NPATCH) :: ZCPL_ICEFLUX
!
! for chemical computations
!
REAL, DIMENSION(KI, IO%NPATCH) :: ZSW_FORBIO
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
INTEGER :: JJ, IBEG, IEND, ISIZE
LOGICAL :: GUPDATED, GALB             ! T if VEGETATION_UPDATE has reset fields
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
  ZALFA(JJ) = ZDIR(JJ) - ISS%XZ0EFFJPDIR(JJ) * XPI/180.

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
IF (LAGRIP .AND. (IO%CPHOTO=='LAI' .OR. IO%CPHOTO=='LST' .OR. IO%CPHOTO=='NIT'.OR. IO%CPHOTO=='NCB') ) THEN
   CALL IRRIGATION_UPDATE(AG, IPM%I, PTSTEP, KMONTH, KDAY, PTIME  )  
ENDIF
!
!* Actualization of the SGH variable (Fmu, Fsat)
!
 CALL ISBA_SGH_UPDATE(IG%XMESH_SIZE, IO, IP, INP, INI, IR, PRAIN )
!
!
!* Actualization of deep soil characteristics
!
IF (LDEEPSOIL) THEN
   CALL DEEPSOIL_UPDATE(INP%XTDEEP, INP%XGAMMAT, INI%TTIME%TDATE%MONTH)
ENDIF
!
!* Actualization of soil and wood carbon spinup
!
! During soil carbon spinup with ISBA-CC: 
!        (1) Atmospheric CO2 concentration fixed to Pre-industrial CO2 consentration XCO2_START
!        (2) Atmospheric CO2 concentration rampin up from XCO2_START to XCO2_END
!
IF(IO%LSPINUPCARBS.OR.IO%LSPINUPCARBW)THEN
!
  ISPINEND = IO%NNBYEARSPINS-NINT(IO%NNBYEARSPINS*XSPIN_CO2)
!  
  IO%LAGRI_TO_GRASS = .FALSE.
!
  IF ( IO%LSPINUPCARBS .AND. (IO%NNBYEARSOLD <= ISPINEND) ) THEN
!
   IO%LAGRI_TO_GRASS = .TRUE.
!
   ZCO2(:) = IO%XCO2_START * 1.E-6 * XMCO2 / XMD
!
  ELSEIF(IO%LSPINUPCARBS .AND. (IO%NNBYEARSOLD > ISPINEND) .AND. (IO%NNBYEARSOLD <= IO%NNBYEARSPINS) )THEN
!
   ZSPINCO2 = IO%XCO2_START + (IO%XCO2_END-IO%XCO2_START) * &
                REAL(IO%NNBYEARSOLD - ISPINEND) / REAL(IO%NNBYEARSPINS - ISPINEND)
!
   ZCO2 (:) = ZSPINCO2 * 1.E-6 * XMCO2 / XMD
!
  ENDIF
!
  CALL CARBON_SPINUP( INI%TTIME, IO )
!
ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Time evolution
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
INI%TTIME%TIME = INI%TTIME%TIME + PTSTEP
 CALL ADD_FORECAST_TO_DATE_SURF(INI%TTIME%TDATE%YEAR,INI%TTIME%TDATE%MONTH,INI%TTIME%TDATE%DAY,INI%TTIME%TIME)
!
! --------------------------------------------------------------------------------------
!
!*      2.     Physical evolution
!
! --------------------------------------------------------------------------------------
! Patch Dependent Calculations
! --------------------------------------------------------------------------------------
!
PATCH_LOOP: DO JPATCH=1,IO%NPATCH
!
  IF (INP%NSIZE_NATURE_P(JPATCH) == 0 ) CYCLE
!
! Pack dummy arguments for each patch:
!
#ifdef TOPD
  IF (LCOUPL_TOPD) NMASKT_PATCH(:) = INP%NR_NATURE_P(:,JPATCH)
#endif
  CALL TREAT_PATCH(INP%NSIZE_NATURE_P(JPATCH),INP%NR_NATURE_P(1:INP%NSIZE_NATURE_P(JPATCH),JPATCH))
!
ENDDO PATCH_LOOP
!
! --------------------------------------------------------------------------------------
! SFX - RRM coupling update if used :
! --------------------------------------------------------------------------------------
!
IF(IO%LCPL_RRM)THEN
  CALL DIAG_CPL_ESM_ISBA(IO, IP, INP, INI, PTSTEP, ZCPL_DRAIN, ZCPL_RUNOFF, &
                         ZCPL_EFLOOD, ZCPL_PFLOOD, ZCPL_IFLOOD, ZCPL_ICEFLUX  )  
ENDIF
!
! --------------------------------------------------------------------------------------
! Vegetation update (in case of non-interactive vegetation):
! Or
! Vegetation albedo only update (in case of interactive vegetation):
! --------------------------------------------------------------------------------------
!
GUPDATED=.FALSE.
!
IF (IO%LVEGUPD) THEN
  GALB = .FALSE. 
  IF (IO%CPHOTO=='LAI'.OR.IO%CPHOTO=='LST'.OR.IO%CPHOTO=='NIT'.OR.IO%CPHOTO=='NCB') GALB = .TRUE.
  CALL VEGETATION_UPDATE(DTCO, DTI, IG%NDIM, IO, IPM%T, IPM%M, IPM%I, IPM%A, &
                         PTSTEP, INI%TTIME, IP%XCOVER, IP%LCOVER,  LAGRIP, &
                         'NAT', GALB, ISS, GUPDATED             )  
!
ENDIF
!
IF(IO%LPERTSURF.AND.GUPDATED) THEN
  ! random perturbation for ensembles:
  ! reset these fields to their original values, as in compute_isba_parameters
  IPM%T%XVEG(:,1) = INI%XPERTVEG(:)
  IPM%T%XLAI(:,1) = INI%XPERTLAI(:)
  IPM%T%XCV (:,1)  = INI%XPERTCV(:)
  ! reapply original perturbation patterns
  WHERE(IPM%T%XALBNIR(:,1)/=XUNDEF)  IPM%T%XALBNIR(:,1) = IPM%T%XALBNIR(:,1) *( 1.+ INI%XPERTALB(:) )
  WHERE(IPM%T%XALBVIS(:,1)/=XUNDEF)  IPM%T%XALBVIS(:,1) = IPM%T%XALBVIS(:,1) *( 1.+ INI%XPERTALB(:) )
  WHERE(IPM%T%XALBUV(:,1)/=XUNDEF)   IPM%T%XALBUV (:,1) = IPM%T%XALBUV (:,1) *( 1.+ INI%XPERTALB(:) )
  WHERE(IPM%T%XZ0(:,1)/=XUNDEF)      IPM%T%XZ0(:,1)     = IPM%T%XZ0(:,1)     *( 1.+ INI%XPERTZ0(:) )
  WHERE(ISS%XZ0EFFIP(:,1)/=XUNDEF) ISS%XZ0EFFIP(:,1) = ISS%XZ0EFFIP(:,1)*( 1.+ INI%XPERTZ0(:) )
  WHERE(ISS%XZ0EFFIM(:,1)/=XUNDEF) ISS%XZ0EFFIM(:,1) = ISS%XZ0EFFIM(:,1)*( 1.+ INI%XPERTZ0(:) )
  WHERE(ISS%XZ0EFFJP(:,1)/=XUNDEF) ISS%XZ0EFFJP(:,1) = ISS%XZ0EFFJP(:,1)*( 1.+ INI%XPERTZ0(:) )
  WHERE(ISS%XZ0EFFJM(:,1)/=XUNDEF) ISS%XZ0EFFJM(:,1) = ISS%XZ0EFFJM(:,1)*( 1.+ INI%XPERTZ0(:) )

ENDIF
!
! --------------------------------------------------------------------------------------
! Outputs for the atmospheric model or update the snow/flood fraction at time t+1
! --------------------------------------------------------------------------------------
! Grid box average fluxes/properties: Arguments and standard diagnostics at time t+1
!
 CALL AVERAGE_FLUX(INP%XPATCH, ZSFTH_TILE, ZSFTQ_TILE, ZSFTS_TILE, &
                   ZSFCO2_TILE, ZSFU_TILE, ZSFV_TILE, PSFTH, PSFTQ,&
                   PSFTS, PSFCO2, PSFU, PSFV    )  
!
!
!-------------------------------------------------------------------------------
!Physical properties see by the atmosphere in order to close the energy budget 
!between surfex and the atmosphere. All variables should be at t+1 but very 
!difficult to do. Maybe it will be done later. However, Ts is at time t+1
!-------------------------------------------------------------------------------
!   
 CALL AVERAGE_PHY(INP%XPATCH, ZTSURF_TILE, ZZ0_TILE, ZZ0H_TILE, &
                  ZQSURF_TILE, PUREF, PZREF, PTSURF, PZ0, PZ0H, PQSURF )
!
!-------------------------------------------------------------------------------------
!Radiative properties at time t+1 (see by the atmosphere) in order to close
!the energy budget between surfex and the atmosphere
!-------------------------------------------------------------------------------------
!
 CALL UPDATE_RAD_ISBA_n(IO, IPM%T, IPM%M, IPM%A, INP, INI, IR, IPM%X%XVEGTYPE, &
                        PZENITH2, PSW_BANDS, ZDIR_ALB_TILE, ZSCA_ALB_TILE, &
                        ZEMIS_TILE, PDIR_SW, PSCA_SW  )
!
 CALL AVERAGE_RAD(INP%XPATCH, ZDIR_ALB_TILE, ZSCA_ALB_TILE, ZEMIS_TILE, &
                  ZTRAD_TILE, PDIR_ALB, PSCA_ALB, INI%XEMIS_NAT,IR%XTSRAD_NAT  )  
!
PEMIS = INI%XEMIS_NAT
PTRAD = IR%XTSRAD_NAT
!
!-------------------------------------------------------------------------------------
!
! Any additional diagnostics (stored in MODD_DIAG_ISBA_n)
!
 CALL AVERAGE_DIAG_ISBA_n(DGI%O, DGI%D, DGI%DC, DGI%DP, DGI%DPC, DGI%O%LSURF_BUDGETC, &
                          IO%LCANOPY, INP%XPATCH, IR%XLE, PUREF, PZREF, PSFCO2, PTRAD)
!
! Cumulated diagnostics (stored in MODD_DIAG_EVAP_ISBA_n)
!
 CALL AVERAGE_DIAG_EVAP_ISBA_n(DGI%O%LSURF_BUDGETC, DGI%DE, DGI%DEC, DGI%DEP, DGI%DEPC, &
                               IO%LGLACIER, IO%LMEB_PATCH, INP%XPATCH, PTSTEP, PRAIN, PSNOW)
!
! Miscellaneous diagnostics (stored in MODD_DIAG_MISC_ISBA_n)
!
 CALL AVERAGE_DIAG_MISC_ISBA_n(DGI%DM, DGI%DMP, IO, INP, IPM%X, IPM%T%XLAI, IR)
!
!--------------------------------------------------------------------------------------
!
 CALL COUPLING_SURF_TOPD(DGI%DE, DGI%DEC, DGI%DC, DGI%DM, IG, &
                         IO, INP, INI, IPM%X, IR, UG, U, HPROGRAM, U%NDIM_FULL)
!
! --------------------------------------------------------------------------------------
! Snow/Flood fractions, albedo and emissivity update :
! --------------------------------------------------------------------------------------
!
! --------------------------------------------------------------------------------------
! Chemical fluxes :
! --------------------------------------------------------------------------------------
!
IF (CHI%SVI%NBEQ>0 .AND. CHI%LCH_BIO_FLUX) THEN
 CALL CH_BVOCEM_n(CHI%SVI, GB, IO, INP, IPM%X%XVEGTYPE, IR%XTG, ZSW_FORBIO, PRHOA, PSFTS)
ENDIF
!
!SOILNOX
IF (CHI%LCH_NO_FLUX) THEN
  CALL SOILEMISNO_n(GB, IP, IR, IPM%T%XLAI, PU, PV)
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
INTEGER, INTENT(IN), DIMENSION(:) :: KMASK
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
REAL, DIMENSION(KSIZE,IO%NNBIOMASS) :: ZP_RESP_BIOMASS_INST         ! instantaneous biomass respiration (kgCO2/kgair m/s)
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
IF (IO%NPATCH==1) THEN
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
   GMEB=IO%LMEB_PATCH(JPATCH)
!
! Pack ISBA input and prognostic variables (modd_isban) for each patch:
!
 CALL PACK_ISBA_PATCH_GET_SIZE_n(IO, IR%TSNOW%SCHEME, PK, JPATCH)
!
 CALL PACK_DIAG_PATCH_GET_SIZE_n(DGI%DE, DGI%O, DGI%DM, IO, IR%TSNOW%SCHEME, PKD, JPATCH)
!
 CALL PACK_ISBA_PATCH_n(AG, IG, ISS, IO, IP, INP, INI, IPM, IR, PK, &
                        KMASK,KSIZE,JPATCH)     
!
! Pack chemistry input and prognostic variables (modd_ch_isban) for each patch:
!
IF (CHI%SVI%NBEQ>0) THEN
  IF( CHI%CCH_DRY_DEP == "WES89") THEN
    CALL PACK_CH_ISBA_PATCH_n(CHI, PKCI, KMASK, KSIZE, JPATCH)     
  END IF
END IF
!
! Allocate ISBA diagnostics for each patch:
!
 CALL PACK_DIAG_PATCH_n(DGI%DE, DGI%O, DGI%D, DGI%DM, IO, IR%TSNOW, INP, PK, PKD, &
                        KSIZE,ISWB,JPATCH)     
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Cosine of the slope typically encoutered in the grid mesh (including subgrid orography)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
ZP_SLOPE_COS(:) = 1./SQRT(1.+PK%ISS%XSSO_SLOPE(:)**2)
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
IF(IO%LFLOOD)THEN
  CALL ISBA_FLOOD_PROPERTIES(PK%I%M%T%XLAI(:,1),PK%I%I%XFFLOOD,PK%I%I%XFFROZEN(:,1),  &
                             ZP_Z0FLOOD, ZP_FFGNOS, ZP_FFVNOS)  
ELSE
  ZP_Z0FLOOD = XUNDEF
  ZP_FFGNOS  = 0.0
  ZP_FFVNOS  = 0.0
ENDIF
!
! For multi-energy balance
   IF(GMEB)THEN
     ZSNOWDEPTH(:) = SUM(PK%I%R%TSNOW%WSNOW(:,:,1)/PK%I%R%TSNOW%RHO(:,:,1),2)
     ZPALPHAN(:)=MEBPALPHAN(ZSNOWDEPTH,PK%I%M%M%XH_VEG(:,1))
   ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Surface Roughness lengths (m):
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!* effective roughness
!
 CALL Z0EFF(IR%TSNOW%SCHEME, IO%CROUGH, GMEB, ZP_ALFA, ZP_ZREF, ZP_UREF,  &
            PK%I%M%T%XZ0(:,1), PK%ISS%XZ0REL, PK%I%R%XPSN(:,1), ZPALPHAN, &
            PK%I%M%M%XZ0LITTER(:,1), PK%I%R%TSNOW%WSNOW(:,1,1),  PK%ISS,  &
            PK%I%I%XFF(:,1), ZP_Z0FLOOD, PK%I%M%X%XZ0_O_Z0H(:,1),         &
            PKD%DGIP%XZ0, PKD%DGIP%XZ0H, PKD%DGIP%XZ0EFF, ZZ0G_WITHOUT_SNOW, &
            ZZ0_MEBV, ZZ0H_MEBV, ZZ0EFF_MEBV, ZZ0_MEBN, ZZ0H_MEBN, ZZ0EFF_MEBN )
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
 CALL ISBA_ALBEDO(IR%TSNOW%SCHEME, IO%LTR_ML, GMEB, ZP_DIR_SW, ZP_SCA_SW, PSW_BANDS,ISWB, &
                  PK%I%M%T, PK%I%M%A, PK%I%R, PK%I%I%XALBF(:,1), PK%I%I%XFFV(:,1), &
                  PK%I%I%XFFG(:,1), ZP_GLOBAL_SW,  ZP_MEB_SCA_SW, ZP_ALBNIR_TVEG, &
                  ZP_ALBVIS_TVEG, ZP_ALBNIR_TSOIL, ZP_ALBVIS_TSOIL                   )  
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Intialize computation of ISBA water and energy budget
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
 CALL ISBA_BUDGET_INIT(DGI%DE%LWATER_BUDGET, IO%CISBA, PK%I%R, PK%I%M%X%XDG(:,:,1), PK%I%IP%XDZG(:,:,1), &
                      ZP_WG_INI, ZP_WGI_INI, ZP_WR_INI, ZP_SWE_INI  )
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Over Natural Land Surfaces:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ZIRRIG_GR(:)= 0.
!
 CALL ISBA(IO, PK%I%M%X, PK%I%M%T, PK%I%M%M, PK%I%M%I, PK%I%P, PK%I%IP, PK%I%I, PK%I%R,       &
           PK%G, PK%AG, PKD%DGI, PKD%DGIP, PKD%DGEI, PKD%DGEIP, PKD%DGMI, IR%TSNOW%SCHEME,    &
           INI%TTIME, INP%XPOI, INP%XABC, PKD%GB%XIACAN(:,:,1), GMEB, PTSTEP, CIMPLICIT_WIND, &
           ZP_ZREF, ZP_UREF, ZP_SLOPE_COS, ZP_TA, ZP_QA, ZP_EXNA, ZP_RHOA, ZP_PS, ZP_EXNS,    &
           ZP_RAIN, ZP_SNOW, ZP_ZENITH, ZP_MEB_SCA_SW, ZP_GLOBAL_SW, ZP_LW, ZP_WIND,          &
           ZP_PEW_A_COEF, ZP_PEW_B_COEF, ZP_PET_A_COEF, ZP_PEQ_A_COEF, ZP_PET_B_COEF,         &
           ZP_PEQ_B_COEF, ZP_ALBNIR_TVEG, ZP_ALBVIS_TVEG, ZP_ALBNIR_TSOIL, ZP_ALBVIS_TSOIL,   &
           ZPALPHAN, ZZ0G_WITHOUT_SNOW, ZZ0_MEBV, ZZ0H_MEBV, ZZ0EFF_MEBV, ZZ0_MEBN, ZZ0H_MEBN,&
           ZZ0EFF_MEBN, ZP_TDEEP_A, ZP_CO2, ZP_FFGNOS, ZP_FFVNOS, ZP_EMIS, ZP_USTAR,          &
           ZP_AC_AGG, ZP_HU_AGG, ZP_RESP_BIOMASS_INST, ZP_DEEP_FLUX, ZIRRIG_GR             )
!
ZP_TRAD=PKD%DGIP%XTSRAD
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Glacier : ice runoff flux (especally for Earth System Model)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
IF(IO%LGLACIER)THEN
!           
  CALL HYDRO_GLACIER(IR%TSNOW%SCHEME, PTSTEP, ZP_SNOW, PK%I%R, PKD%DGEIP%XICEFLUX)
!     
ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Calculation of ISBA water and energy budget (and time tendencies of each reservoir)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
CALL ISBA_BUDGET(IO, PK%I%R, PKD%DGEIP, DGI%DE%LWATER_BUDGET, &
                 IR%TSNOW%SCHEME, PTSTEP, PK%I%M%X%XDG(:,:,1),&
                 PK%I%IP%XDZG(:,:,1), ZP_WG_INI, ZP_WGI_INI,  &
                 ZP_WR_INI, ZP_SWE_INI, ZP_RAIN, ZP_SNOW, PKD%DGIP%XEVAP   )
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Evolution of soil albedo, when depending on surface soil wetness:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
IF (IO%CALBEDO=='EVOL' .AND. IO%LECOCLIMAP) THEN
  CALL SOIL_ALBEDO(IO%CALBEDO, PK%I%IP%XWSAT(:,1),PK%I%R%XWG(:,1,:), PK%I%IP, PK%I%M%A, "ALL")  
  !
  CALL ALBEDO(IO%CALBEDO, PK%I%M%T, PK%I%M%A )  
END IF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Vegetation evolution for interactive LAI
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
IF (IO%CPHOTO=='LAI' .OR. IO%CPHOTO=='LST' .OR. IO%CPHOTO=='NIT' .OR. IO%CPHOTO=='NCB') THEN
  CALL VEGETATION_EVOL(IO, PK%I%IP, PK%I%M%X, PK%I%M%T, PK%I%M%A, PK%I%M%I, PK%I%R, &
                       LAGRIP, PTSTEP, KMONTH, KDAY, PTIME, PK%G%XLAT, ZP_RHOA, ZP_CO2, &
                       PK%ISS, ZP_RESP_BIOMASS_INST,  &
                       ! add optional for accurate dependency to nitrogen
                       ! limitation
                        PSWDIR=ZP_GLOBAL_SW ) 
END IF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Diagnostic of respiration carbon fluxes and soil carbon evolution
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!ii
ZP_SFCO2    (:)=0.
PKD%DGEIP%XRESP_ECO (:)=0.
PKD%DGEIP%XRESP_AUTO(:)=0.
!
IF ( IO%CPHOTO/='NON' .AND. IO%CRESPSL/='NON' .AND. ANY(PK%I%M%T%XLAI(:,1)/=XUNDEF) ) THEN
  CALL CARBON_EVOL(IO, PK%I%P, PK%I%IP, PK%I%M%X, PK%I%M%T, PK%I%R,  PKD%DGEIP, &
                   PTSTEP, ZP_RHOA, ZP_RESP_BIOMASS_INST )  
  ! calculation of vegetation CO2 flux
  ! Positive toward the atmosphere
  ZP_SFCO2(:) = PKD%DGEIP%XRESP_ECO(:) - PKD%DGEIP%XGPP(:)  
END IF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Reset effecitve roughness lentgh to its nominal value when snow has just disappeared
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
 CALL SUBSCALE_Z0EFF(PK%ISS,PK%I%M%T%XZ0(:,:),.FALSE.,&
                     OMASK=(PK%I%R%TSNOW%WSNOW(:,1,1)==0. .AND. PK%I%R%XPSN(:,1)>0.)  )   
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Turbulent fluxes
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
ZP_SFTH(:) = PKD%DGIP%XH(:)
ZP_SFTQ(:) = PKD%DGIP%XEVAP(:)

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
IF (CHI%SVI%NBEQ>0) THEN
  IF( CHI%CCH_DRY_DEP == "WES89") THEN

    IBEG = CHI%SVI%NSV_CHSBEG
    IEND = CHI%SVI%NSV_CHSEND 
    ISIZE = IEND - IBEG + 1 

    CALL CH_DEP_ISBA(PK%I%P, PK%I%IP, PK%I%M%T, PK%I%R, PKD%DGIP, PKD%DGMI, PKCI, &
                     ZP_USTAR, ZP_TA, ZP_PA, ZP_TRAD(:), ISIZE )  
 
    ZP_SFTS(:,IBEG:IEND) = - ZP_SV(:,IBEG:IEND) * PKCI%XDEP(:,1:CHI%SVI%NBEQ)  

    IF (CHI%SVI%NAEREQ > 0 ) THEN
       
      IBEG = CHI%SVI%NSV_AERBEG
      IEND = CHI%SVI%NSV_AEREND
      CALL CH_AER_DEP(ZP_SV(:,IBEG:IEND), ZP_SFTS(:,IBEG:IEND), &
                      ZP_USTAR, PK%I%R%XRESA(:,1), ZP_TA, ZP_RHOA)     
    END IF
  ELSE

    IBEG = CHI%SVI%NSV_AERBEG
    IEND = CHI%SVI%NSV_AEREND
    ZP_SFTS(:,IBEG:IEND) = 0.
    ZP_SFTS(:,IBEG:IEND) = 0.

  ENDIF
ENDIF
!
! --------------------------------------------------------------------------------------
! Dust deposition and emission:
! --------------------------------------------------------------------------------------
!
IF(CHI%SVI%NDSTEQ>0)THEN

  IBEG = CHI%SVI%NSV_DSTBEG
  IEND = CHI%SVI%NSV_DSTEND
  IDST = IEND - IBEG + 1

  CALL COUPLING_DST_n(DST, PK%I%P, PK%I%IP, PK%I%R, PKD%DGIP, &
            HPROGRAM,                    &!I [char] Name of program
            KSIZE,      &!I [nbr] number of points in patch
            IDST,                        &!I [nbr] number of dust emissions variables
            JPATCH,                      &!I [idx] patch in question
            ZP_PS,                       &!I [Pa] surface pressure
            ZP_QA,                       &!I [kg/kg] specific humidity
            ZP_RHOA,                     &!I [kg/m3] atmospheric density
            ZP_PA,                       &!I [K] Atmospheric pressure
            ZP_TA,                       &!I [K] Atmospheric temperature
            ZP_U,                        &!I [m/s] zonal wind at atmospheric height 
            ZP_UREF,                     &!I [m] reference height of wind
            ZP_V,                        &!I [m/s] meridional wind at atmospheric height
            ZP_ZREF,                     &!I [m] reference height of wind
            ZP_SFTS(:,IBEG:IEND)  &!O [kg/m2/sec] flux of dust            
            )  
!
   IF (CHI%SVI%NSV_AEREND > 0)  THEN ! case of dust/ anthropogenic aerosols coupling

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
           ZP_SFTS(:,JSV) = ZP_SFTS(:,JSV) + ZP_SFTS(:,IBEG-1+JSV_IDX+2)*XAVOGADRO/XMOLARWEIGHT_DST
         END IF
         IF ( (TRIM(HSV(JSV)) == "@DSTJ").AND.(JMODE==2)) THEN 
           ! add dust flux and conversion kg/m2/sec into molec.m2/s
           ZP_SFTS(:,JSV) = ZP_SFTS(:,JSV) + ZP_SFTS(:,IBEG-1+JSV_IDX+2)*XAVOGADRO/XMOLARWEIGHT_DST
         END IF
       END DO
       !
     END DO
    END IF
!    
!Modify fluxes due to dry deposition, we introduce a negative flux where dust is lost
  CALL DSLT_DEP(ZP_SV(:,IBEG:IEND), ZP_SFTS(:,IBEG:IEND), ZP_USTAR, PK%I%R%XRESA(:,1), &
                ZP_TA, ZP_RHOA, DST%XEMISSIG_DST, DST%XEMISRADIUS_DST, JPMODE_DST,     &
                XDENSITY_DST, XMOLARWEIGHT_DST, ZCONVERTFACM0_DST, ZCONVERTFACM6_DST,  &
                ZCONVERTFACM3_DST, LVARSIG_DST, LRGFIX_DST, CVERMOD  )
!
!Transfer these fluxes to fluxes understandable by all moments
  CALL MASSFLUX2MOMENTFLUX(           &
    ZP_SFTS(:,IBEG:IEND),             & !I/O ![kg/m2/sec] In: flux of only mass, out: flux of moments
    ZP_RHOA,                          & !I [kg/m3] air density
    DST%XEMISRADIUS_DST,              & !I [um] emitted radius for the modes (max 3)
    DST%XEMISSIG_DST,                 & !I [-] emitted sigma for the different modes (max 3)
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
IF (CHI%SVI%NSLTEQ>0) THEN
  !
  IBEG = CHI%SVI%NSV_SLTBEG
  IEND = CHI%SVI%NSV_SLTEND
  !
  CALL DSLT_DEP(ZP_SV(:,IBEG:IEND), ZP_SFTS(:,IBEG:IEND), ZP_USTAR, PK%I%R%XRESA(:,1), &
                ZP_TA, ZP_RHOA, SLT%XEMISSIG_SLT, SLT%XEMISRADIUS_SLT, JPMODE_SLT,     &
                XDENSITY_SLT, XMOLARWEIGHT_SLT, ZCONVERTFACM0_SLT, ZCONVERTFACM6_SLT,  &
                ZCONVERTFACM3_SLT, LVARSIG_SLT, LRGFIX_SLT, CVERMOD  )  

  CALL MASSFLUX2MOMENTFLUX(           &
    ZP_SFTS(:,IBEG:IEND),             & !I/O ![kg/m2/sec] In: flux of only mass, out: flux of moments
    ZP_RHOA,                          & !I [kg/m3] air density
    SLT%XEMISRADIUS_SLT,              & !I [um] emitted radius for the modes (max 3)
    SLT%XEMISSIG_SLT,                 & !I [-] emitted sigma for the different modes (max 3)
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
 CALL DIAG_INLINE_ISBA_n(DGI%O, PK%I%I, PKD%DGIP, IO%LCANOPY, ZP_TA, ZP_QA, ZP_PA, &
                         ZP_PS, ZP_RHOA, ZP_U, ZP_V, ZP_ZREF, ZP_UREF, ZP_SFTH, &
                         ZP_SFTQ, ZP_SFU, ZP_SFV, ZP_DIR_SW, ZP_SCA_SW, ZP_LW )  
!
!
!-------------------------------------------------------------------------------
!Physical properties see by the atmosphere in order to close the energy budget 
!between surfex and the atmosphere. All variables should be at t+1 but very 
!difficult to do. Maybe it will be done later. However, Ts can be at time t+1
!-------------------------------------------------------------------------------
!
ZP_TSURF (:) = PKD%DGIP%XTS (:)
ZP_Z0    (:) = PKD%DGIP%XZ0 (:)
ZP_Z0H   (:) = PKD%DGIP%XZ0H(:)
ZP_QSURF (:) = PKD%DGIP%XQS (:)
!
!-------------------------------------------------------------------------------
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Isba offline diagnostics for each patch
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
 CALL DIAG_EVAP_ISBA_n(DGI%O%LSURF_BUDGETC, DGI%DE, DGI%DEP, DGI%DEPC, DGI%DPC, PKD%DGEIP, &
                       PKD%DGIP, PK%I%R, IO, IR%TSNOW%SCHEME, PTSTEP, KMASK, KSIZE,JPATCH, ZP_RHOA)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Isba offline diagnostics for miscellaneous terms over each patch
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
 CALL DIAG_MISC_ISBA_n(DGI%DMP, PKD%DGMI, PK%I%IP, PK%I%I, PK%I%M%X, PK%I%R, PK%AG, IO, &
                       DGI%DM%LSURF_MISC_BUDGET, PTSTEP, IR%TSNOW%SCHEME, LAGRIP, PTIME, &
                       KSIZE, JPATCH, KMASK     )                  
!
! Unpack ISBA diagnostics (modd_diag_isban) for each patch:ISIZE_MAX = MAXVAL(NSIZE_NATURE_P)

!  (MUST be done BEFORE UNPACK_ISBA_PATCH, because of XP_LE)
!
 CALL UNPACK_DIAG_PATCH_n(DGI%O, DGI%DP, GB, IO, IR%TSNOW, PKD, PK,    &
                          KMASK, KSIZE, IO%NPATCH, JPATCH, ZCPL_DRAIN, &
                          ZCPL_RUNOFF, ZCPL_EFLOOD, ZCPL_PFLOOD,       &
                          ZCPL_IFLOOD, ZCPL_ICEFLUX)  
!
! for chemical deposition
!
IF (CHI%SVI%NBEQ>0) THEN
  IF( CHI%CCH_DRY_DEP == "WES89") THEN
    CALL UNPACK_CH_ISBA_PATCH_n(CHI, PKCI, KMASK, KSIZE, JPATCH)     
  END IF
END IF
!
! Unpack ISBA variables (modd_isban) for each patch:
!
 CALL UNPACK_ISBA_PATCH_n(AG, IO, INP, IPM%T, IPM%A, IR, ISS, PK, KMASK, KSIZE, JPATCH)
!
!----------------------------------------------------------------------
!
! for further chemical biogenic emissions
!
IF (CHI%SVI%NBEQ>0 .AND. CHI%LCH_BIO_FLUX) THEN
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
IF (IO%NPATCH==1) THEN
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
IF (CHI%SVI%NDSTEQ>0) THEN
  DO JSV = 1,NDSTMDE
    IF (IMOMENT == 1) THEN
      DST%XSFDST(:,JSV,JPATCH) = ZSFTS_TILE(:,NDST_MDEBEG+JSV-1,JPATCH)
    ELSE
      DST%XSFDST(:,JSV,JPATCH) = ZSFTS_TILE(:,NDST_MDEBEG+(JSV-1)*IMOMENT+1,JPATCH)
    END IF

    DST%XSFDSTM(:,JSV,JPATCH) = DST%XSFDSTM(:,JSV,JPATCH) + DST%XSFDST(:,JSV,JPATCH) * PTSTEP
  ENDDO
ENDIF
!
IF (LHOOK) CALL DR_HOOK('COUPLING_ISBA_n:TREAT_PATCH',1,ZHOOK_HANDLE)
!
END SUBROUTINE TREAT_PATCH
!==========================================================================================
END SUBROUTINE COUPLING_ISBA_n
