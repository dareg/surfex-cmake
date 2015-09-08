!#############################################################
SUBROUTINE INIT_ISBA_n (AG, BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTCO, &
                            DTI, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, &
                            DGSI, DGU, DGT, DGUT, DGW, DST, F, FSB, GB, IOB, ICP, &
                            IG, I, O, S, SSB, SLT, UG, U, SV, TCP, TGD, &
                            TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                            HPROGRAM,HINIT,OLAND_USE,                    &
                             KI,KSV,KSW,                                &
                             HSV,PCO2,PRHOA,                            &
                             PZENITH,PAZIM,PSW_BANDS,PDIR_ALB,PSCA_ALB, &
                             PEMIS,PTSRAD, PTSURF,                      &
                             KYEAR, KMONTH,KDAY, PTIME,                 &
                             HATMFILE,HATMFILETYPE,                     &
                             HTEST                                      )  
!#############################################################
!
!!****  *INIT_ISBA_n* - routine to initialize ISBA
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      Modified by P. Le Moigne (11/2004): miscellaneous diagnostics
!!      Modified by P. Le Moigne (06/2006): seeding and irrigation    
!!      Modified by B. Decharme    (2008) : SGH and Flooding scheme
!!      Modified by B. Decharme  (01/2009): optional deep soil temperature as in Arpege
!!      Modified by R. Hamdi     (01/2009): Cp and L
!!      Modified by B. Decharme  (06/2009): read topographic index statistics
!!      Modified by P. Le Moigne (01/2009): Beljaars sso
!!      Modified by B. Decharme  (08/2009): Active Trip coupling variable if Earth System Model
!!      A.L. Gibelin   04/09 : change BSLAI_NITRO initialisation
!!      A.L. Gibelin   04/09 : modifications for CENTURY model 
!!      A.L. Gibelin   06/09 : soil carbon initialisation
!!      B. Decharme    07/11 : read pgd+prep
!!      R. Alkama      05/12 : new carbon spinup
!!      J.Escobar      11/13 : add USE MODI_DEFAULT_CROCUS
!!      B. Decharme  04/2013 new coupling variables
!!      P. Samuelsson  10/14 : MEB
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
!
!
!
!
USE MODD_AGRI_n, ONLY : AGRI_t
USE MODD_BEM_OPTION_n, ONLY : BEM_OPTIONS_t
USE MODD_BLD_DESCRIPTION_n, ONLY : BLD_DESC_t
USE MODD_CH_EMIS_FIELD_n, ONLY : CH_EMIS_FIELD_t
USE MODD_CH_ISBA_n, ONLY : CH_ISBA_t
USE MODD_CH_SEAFLUX_n, ONLY : CH_SEAFLUX_t
USE MODD_CH_SNAP_n, ONLY : CH_EMIS_SNAP_t
USE MODD_CH_SURF_n, ONLY : CH_SURF_t
USE MODD_CH_TEB_n, ONLY : CH_TEB_t
USE MODD_CH_WATFLUX_n, ONLY : CH_WATFLUX_t
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_DATA_SEAFLUX_n, ONLY : DATA_SEAFLUX_t
USE MODD_DATA_TEB_n, ONLY : DATA_TEB_t
USE MODD_DATA_TSZ0_n, ONLY : DATA_TSZ0_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_FLAKE_n, ONLY : DIAG_FLAKE_t
USE MODD_DIAG_ISBA_n, ONLY : DIAG_ISBA_t
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t
USE MODD_DIAG_MISC_TEB_OPTION_n, ONLY : DIAG_MISC_TEB_OPTIONS_t
USE MODD_DIAG_OCEAN_n, ONLY : DIAG_OCEAN_t
USE MODD_DIAG_SEAFLUX_n, ONLY : DIAG_SEAFLUX_t
USE MODD_DIAG_SEAICE_n, ONLY : DIAG_SEAICE_t
USE MODD_DIAG_SURF_ATM_n, ONLY : DIAG_SURF_ATM_t
USE MODD_DIAG_TEB_n, ONLY : DIAG_TEB_t
USE MODD_DIAG_UTCI_TEB_n, ONLY : DIAG_UTCI_TEB_t
USE MODD_DIAG_WATFLUX_n, ONLY : DIAG_WATFLUX_t
USE MODD_DST_n, ONLY : DST_t
USE MODD_FLAKE_n, ONLY : FLAKE_t
USE MODD_FLAKE_SBL_n, ONLY : FLAKE_SBL_t
USE MODD_GR_BIOG_n, ONLY : GR_BIOG_t
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
USE MODD_ISBA_CANOPY_n, ONLY : ISBA_CANOPY_t
USE MODD_ISBA_GRID_n, ONLY : ISBA_GRID_t
USE MODD_ISBA_n, ONLY : ISBA_t
USE MODD_OCEAN_n, ONLY : OCEAN_t
USE MODD_SEAFLUX_n, ONLY : SEAFLUX_t
USE MODD_SEAFLUX_SBL_n, ONLY : SEAFLUX_SBL_t
USE MODD_SLT_n, ONLY : SLT_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SV_n, ONLY : SV_t
USE MODD_TEB_CANOPY_n, ONLY : TEB_CANOPY_t
USE MODD_TEB_GARDEN_n, ONLY : TEB_GARDEN_t
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TEB_GARDEN_OPTIONS_t
USE MODD_TEB_GREENROOF_n, ONLY : TEB_GREENROOF_t
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TEB_GREENROOF_OPTIONS_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_OPTIONS_t
USE MODD_WATFLUX_n, ONLY : WATFLUX_t
USE MODD_WATFLUX_SBL_n, ONLY : WATFLUX_SBL_t
!
USE MODD_DATA_COVER,     ONLY : XDATA_LAI, XDATA_H_TREE,                          &
                                XDATA_ALBNIR_VEG, XDATA_ALBVIS_VEG,               &
                                XDATA_ALBUV_VEG, XDATA_RSMIN,                     &
                                XDATA_ROOT_EXTINCTION,XDATA_ROOT_LIN,             &
                                XDATA_RGL, XDATA_CV, XDATA_GAMMA, XDATA_GMES,     &
                                XDATA_GC, XDATA_BSLAI, XDATA_SEFOLD, XDATA_LAIMIN,&
                                XDATA_DMAX, XDATA_STRESS, XDATA_F2I,              &
                                XDATA_VEG, XDATA_GREEN, XDATA_Z0, XDATA_Z0_O_Z0H, &
                                XDATA_EMIS_ECO, XDATA_WRMAX_CF,                   &
                                XDATA_CE_NITRO,XDATA_CF_NITRO,XDATA_CNA_NITRO,    &
                                XDATA_SOILRC_SO2, XDATA_SOILRC_O3, XDATA_RE25,    &
                                XDATA_GMES_ST, XDATA_BSLAI_ST, XDATA_SEFOLD_ST,   &
                                XDATA_GC_ST, XDATA_DMAX_ST
!
!

USE MODD_SURF_PAR,       ONLY : XUNDEF, NUNDEF
USE MODD_AGRI,           ONLY : LAGRIP
!
USE MODE_TARTES, ONLY : INIT_TARTES
USE MODE_SNOWCRO_FLANNER, ONLY : READ_FZ06
!
USE MODD_READ_NAMELIST,  ONLY : LNAM_READ
!
USE MODD_CO2V_PAR,  ONLY : XMCO2, XSPIN_CO2
USE MODD_CSTS,      ONLY : XMD
!
USE MODI_INIT_IO_SURF_n
!
USE MODI_GET_LUOUT
USE MODI_ABOR1_SFX
USE MODI_DEFAULT_ISBA
USE MODI_DEFAULT_CH_DEP
USE MODI_DEFAULT_CH_BIO_FLUX
USE MODI_DEFAULT_DIAG_ISBA
USE MODI_DEFAULT_CROCUS
USE MODI_READ_DEFAULT_ISBA_n
USE MODI_READ_ISBA_CONF_n
USE MODI_READ_PREP_ISBA_SNOW
USE MODI_READ_PREP_ISBA_CARBON
USE MODI_READ_SURF
USE MODI_PREP_CTRL_ISBA
USE MODI_READ_ISBA_DATE
USE MODI_READ_PGD_ISBA_n
USE MODI_COMPUTE_ISBA_PARAMETERS
USE MODI_READ_NAM_PREP_ISBA_n
USE MODI_INI_DATA_PARAM
!
USE MODI_SET_SURFEX_FILEIN
!
USE MODI_END_IO_SURF_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
TYPE(AGRI_t), INTENT(INOUT) :: AG
TYPE(BEM_OPTIONS_t), INTENT(INOUT) :: BOP
TYPE(BLD_DESC_t), INTENT(INOUT) :: BDD
TYPE(CH_EMIS_FIELD_t), INTENT(INOUT) :: CHE
TYPE(CH_ISBA_t), INTENT(INOUT) :: CHI
TYPE(CH_SEAFLUX_t), INTENT(INOUT) :: CHS
TYPE(CH_EMIS_SNAP_t), INTENT(INOUT) :: CHN
TYPE(CH_SURF_t), INTENT(INOUT) :: CHU
TYPE(CH_TEB_t), INTENT(INOUT) :: CHT
TYPE(CH_WATFLUX_t), INTENT(INOUT) :: CHW
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTI
TYPE(DATA_SEAFLUX_t), INTENT(INOUT) :: DTS
TYPE(DATA_TEB_t), INTENT(INOUT) :: DTT
TYPE(DATA_TSZ0_t), INTENT(INOUT) :: DTZ
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEI
TYPE(DIAG_FLAKE_t), INTENT(INOUT) :: DGF
TYPE(DIAG_ISBA_t), INTENT(INOUT) :: DGI
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DGMI
TYPE(DIAG_MISC_TEB_OPTIONS_t), INTENT(INOUT) :: DGMTO
TYPE(DIAG_OCEAN_t), INTENT(INOUT) :: DGO
TYPE(DIAG_SEAFLUX_t), INTENT(INOUT) :: DGS
TYPE(DIAG_SEAICE_t), INTENT(INOUT) :: DGSI
TYPE(DIAG_SURF_ATM_t), INTENT(INOUT) :: DGU
TYPE(DIAG_TEB_t), INTENT(INOUT) :: DGT
TYPE(DIAG_UTCI_TEB_t), INTENT(INOUT) :: DGUT
TYPE(DIAG_WATFLUX_t), INTENT(INOUT) :: DGW
TYPE(DST_t), INTENT(INOUT) :: DST
TYPE(FLAKE_t), INTENT(INOUT) :: F
TYPE(FLAKE_SBL_t), INTENT(INOUT) :: FSB
TYPE(GR_BIOG_t), INTENT(INOUT) :: GB
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
TYPE(ISBA_CANOPY_t), INTENT(INOUT) :: ICP
TYPE(ISBA_GRID_t), INTENT(INOUT) :: IG
TYPE(ISBA_t), INTENT(INOUT) :: I
TYPE(OCEAN_t), INTENT(INOUT) :: O
TYPE(SEAFLUX_t), INTENT(INOUT) :: S
TYPE(SEAFLUX_SBL_t), INTENT(INOUT) :: SSB
TYPE(SLT_t), INTENT(INOUT) :: SLT
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SV_t), INTENT(INOUT) :: SV
TYPE(TEB_CANOPY_t), INTENT(INOUT) :: TCP
TYPE(TEB_GARDEN_t), INTENT(INOUT) :: TGD
TYPE(TEB_GARDEN_OPTIONS_t), INTENT(INOUT) :: TGDO
TYPE(TEB_GREENROOF_t), INTENT(INOUT) :: TGR
TYPE(TEB_GREENROOF_OPTIONS_t), INTENT(INOUT) :: TGRO
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_VEG_OPTIONS_t), INTENT(INOUT) :: TVG
TYPE(WATFLUX_t), INTENT(INOUT) :: W
TYPE(WATFLUX_SBL_t), INTENT(INOUT) :: WSB
!
 CHARACTER(LEN=6),                 INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=3),                 INTENT(IN)  :: HINIT     ! choice of fields to initialize
LOGICAL,                          INTENT(IN)  :: OLAND_USE !
INTEGER,                          INTENT(IN)  :: KI        ! number of points
INTEGER,                          INTENT(IN)  :: KSV       ! number of scalars
INTEGER,                          INTENT(IN)  :: KSW       ! number of short-wave spectral bands
 CHARACTER(LEN=6), DIMENSION(KSV), INTENT(IN)  :: HSV       ! name of all scalar variables
REAL,             DIMENSION(KI),  INTENT(IN)  :: PCO2      ! CO2 concentration (kg/m3)
REAL,             DIMENSION(KI),  INTENT(IN)  :: PRHOA     ! air density
REAL,             DIMENSION(KI),  INTENT(IN)  :: PZENITH   ! solar zenithal angle
REAL,             DIMENSION(KI),  INTENT(IN)  :: PAZIM     ! solar azimuthal angle (rad from N, clock)
REAL,             DIMENSION(KSW), INTENT(IN)  :: PSW_BANDS ! middle wavelength of each band
REAL,             DIMENSION(KI,KSW),INTENT(OUT) :: PDIR_ALB  ! direct albedo for each band
REAL,             DIMENSION(KI,KSW),INTENT(OUT) :: PSCA_ALB  ! diffuse albedo for each band
REAL,             DIMENSION(KI),  INTENT(OUT) :: PEMIS     ! emissivity
REAL,             DIMENSION(KI),  INTENT(OUT) :: PTSRAD    ! radiative temperature
REAL,             DIMENSION(KI),  INTENT(OUT) :: PTSURF    ! surface effective temperature         (K)
!
INTEGER,                          INTENT(IN)  :: KYEAR     ! current year (UTC)
INTEGER,                          INTENT(IN)  :: KMONTH    ! current month (UTC)
INTEGER,                          INTENT(IN)  :: KDAY      ! current day (UTC)
REAL,                             INTENT(IN)  :: PTIME     ! current time since
                                                          !  midnight (UTC, s)
!
 CHARACTER(LEN=28),                INTENT(IN)  :: HATMFILE    ! atmospheric file name
 CHARACTER(LEN=6),                 INTENT(IN)  :: HATMFILETYPE! atmospheric file type
 CHARACTER(LEN=2),                 INTENT(IN)  :: HTEST       ! must be equal to 'OK'
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL, DIMENSION(KI) :: ZCO2     ! CO2 concentration  (kg/m3)
REAL                :: ZSPINCO2
INTEGER             :: ISPINEND
!
INTEGER             :: ILUOUT   ! unit of output listing file
INTEGER             :: IVERSION       ! surface version
INTEGER             :: IRESP   ! return code
INTEGER             :: ISIZE_LMEB_PATCH   ! Number of patches where multi-energy balance should be applied
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!               Initialisation for IO
!
!
IF (LHOOK) CALL DR_HOOK('INIT_ISBA_N',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('INIT_ISBAN: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF
!
!               Other little things
!
IF (LNAM_READ) THEN
 !
 !*       0.     Defaults
!               --------

 !        0.1. Hard defaults
 !      
 CALL DEFAULT_ISBA(I%XTSTEP, I%XOUT_TSTEP,                           &
                     I%CROUGH,I%CRUNOFF,I%CALBEDO,I%CSCOND,              &
                     I%CC1DRY, I%CSOILFRZ, I%CDIFSFCOND, I%CSNOWRES,     &
                     I%CCPSURF, I%XCGMAX, I%XCDRAG, I%CKSAT, I%LSOC,       &
                     I%CRAIN, I%CHORT, I%LGLACIER, I%LCANOPY_DRAG,       &
                     I%LVEGUPD, I%LSPINUPCARBS, I%LSPINUPCARBW,        &
                     I%XSPINMAXS, I%XSPINMAXW, I%XCO2_START, I%XCO2_END, &
                     I%NNBYEARSPINS, I%NNBYEARSPINW, I%LNITRO_DILU     )
 !                  
 CALL DEFAULT_CH_DEP(CHI%CCH_DRY_DEP)
 CALL DEFAULT_CH_BIO_FLUX(CHI%LCH_BIO_FLUX)                  
 CALL DEFAULT_DIAG_ISBA(DGI%N2M,DGI%LSURF_BUDGET,DGI%L2M_MIN_ZS,DGI%LRAD_BUDGET,   &
                        DGI%LCOEF,DGI%LSURF_VARS,DGEI%LSURF_EVAP_BUDGET,        &
                        DGMI%LSURF_MISC_BUDGET,DGMI%LSURF_DIAG_ALBEDO,       &
                        DGEI%LSURF_BUDGETC,DGMI%LSURF_MISC_DIF,DGI%LPATCH_BUDGET,&
                        DGI%LPGD,DGEI%LRESET_BUDGETC,DGEI%LWATER_BUDGET,         &
                        DGI%XDIAG_TSTEP                                )  
 !
 CALL DEFAULT_CROCUS(I%LSNOWDRIFT,I%LSNOWDRIFT_SUBLIM,I%LSNOW_ABS_ZENITH,&
                     I%CSNOWMETAMO,I%CSNOWRAD)
 ! 
ENDIF
!
!        0.2. Defaults from file header
!    
 CALL READ_DEFAULT_ISBA_n(CHI, DGEI, DGI, DGMI, I, &
                          HPROGRAM)
!
 CALL READ_ISBA_CONF_n(CHI, DGEI, DGI, DGMI, I, &
                       HPROGRAM)
!
 CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                     HPROGRAM,'FULL  ','ISBA  ','READ ')
 CALL READ_SURF(IOB, &
                HPROGRAM,'VERSION',IVERSION,IRESP)
 CALL END_IO_SURF_n(HPROGRAM)
!
!*       1.     Reading of configuration:
!               -------------------------
!
!* initialization of snow and carbon schemes
!
I%NNBYEARSOLD = 1
I%NSPINS      = 1
I%NSPINW      = 1
!
IF (HINIT=='PRE') THEN 
  CALL READ_PREP_ISBA_SNOW(HPROGRAM,I%TSNOW%SCHEME,I%TSNOW%NLAYER) 
!
!* initialization of soil carbon scheme
!
  CALL READ_PREP_ISBA_CARBON(HPROGRAM,I%CRESPSL)
!
  IF (I%CRESPSL=='CNT') THEN
    I%NNLITTER = 2
    I%NNLITTLEVS = 2
    I%NNSOILCARB = 3
  ELSE
    I%NNLITTER = 0
    I%NNLITTLEVS = 0
    I%NNSOILCARB = 0
  ENDIF

ELSEIF (HINIT=='ALL') THEN
!
  CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                     HPROGRAM,'NATURE','ISBA  ','READ ')
!
  IF (IVERSION<6) THEN
    I%CRESPSL='DEF'
  ELSE  
    CALL READ_SURF(IOB, &
                HPROGRAM,'RESPSL',I%CRESPSL,IRESP)
    CALL READ_SURF(IOB, &
                HPROGRAM,'NLITTER',I%NNLITTER,IRESP)
    CALL READ_SURF(IOB, &
                HPROGRAM,'NLITTLEVS',I%NNLITTLEVS,IRESP)
    CALL READ_SURF(IOB, &
                HPROGRAM,'NSOILCARB',I%NNSOILCARB,IRESP)
    IF(IVERSION>=7.AND.(I%LSPINUPCARBS.OR.I%LSPINUPCARBW))THEN
      CALL READ_SURF(IOB, &
                HPROGRAM,'NBYEARSOLD',I%NNBYEARSOLD,IRESP)
    ELSE
      I%NNBYEARSOLD=NUNDEF
    ENDIF
  ENDIF
!
  CALL END_IO_SURF_n(HPROGRAM)
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       2.     Physiographic fields
!               --------------------
!
!
!* date
!
SELECT CASE (HINIT)
  CASE ('PGD')
    I%TTIME%TDATE%YEAR = NUNDEF
    I%TTIME%TDATE%MONTH= NUNDEF
    I%TTIME%TDATE%DAY  = NUNDEF
    I%TTIME%TIME       = XUNDEF

  CASE ('PRE')
    CALL PREP_CTRL_ISBA(DGI%N2M,DGI%LSURF_BUDGET,DGI%L2M_MIN_ZS,DGI%LRAD_BUDGET,DGI%LCOEF,DGI%LSURF_VARS,&
                          DGEI%LSURF_EVAP_BUDGET,DGMI%LSURF_MISC_BUDGET,DGEI%LSURF_BUDGETC,     &
                          DGI%LPATCH_BUDGET,DGMI%LSURF_MISC_DIF,ILUOUT                    )    
    IF (LNAM_READ) CALL READ_NAM_PREP_ISBA_n(HPROGRAM)                        
    CALL READ_ISBA_DATE(IOB, &
                        HPROGRAM,HINIT,ILUOUT,HATMFILE,HATMFILETYPE,KYEAR,KMONTH,KDAY,PTIME,I%TTIME)

  CASE DEFAULT
    CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                     HPROGRAM,'NATURE','ISBA  ','READ ')
    CALL READ_SURF(IOB, &
                HPROGRAM,'DTCUR',I%TTIME,IRESP)
    CALL END_IO_SURF_n(HPROGRAM)
END SELECT
!
!-----------------------------------------------------------------------------------------------------
! READ PGD FILE
!-----------------------------------------------------------------------------------------------------
!
! initialization for I/O
!
 CALL SET_SURFEX_FILEIN(HPROGRAM,'PGD ') ! change input file name to pgd name
 CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                     HPROGRAM,'NATURE','ISBA  ','READ ')
!
!
!*       2.1    Cover, soil and orographic fields:
!               ---------------------------------
!
 CALL READ_PGD_ISBA_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTCO, DTI, &
                                  DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, &
                                  DGU, DGT, DGUT, DGW, F, FSB, GB, IOB, ICP, IG, I, &
                                  O, S, SSB, UG, U, SV, TCP, TGD, TGDO, TGR, TGRO, &
                                  T, TOP, TVG, W, WSB, &
                      HPROGRAM,OLAND_USE)
!
ISIZE_LMEB_PATCH=COUNT(I%LMEB_PATCH(:))
!
!
!*       2.2    Check:
!               ------
!
IF ( I%CPHOTO/='NON' .AND. I%NPATCH/=12 .AND. I%NPATCH/=19 )THEN
  CALL ABOR1_SFX('INIT_ISBAN: INCONSISTENCY BETWEEN CPHOTO AND NPATCH')
ENDIF
!
IF (HINIT=='PRE' .AND. I%TSNOW%SCHEME.NE.'3-L' .AND. I%TSNOW%SCHEME.NE.'CRO' .AND. I%CISBA=='DIF') THEN
    CALL ABOR1_SFX("INIT_ISBAN: WITH CISBA = DIF, CSNOW MUST BE 3-L OR CRO")
ENDIF
IF ( I%CPHOTO/='LAI' .AND. I%CPHOTO/='LST' .AND. I%CPHOTO/='NIT' .AND. I%CPHOTO/='NCB' .AND. LAGRIP) THEN
  CALL ABOR1_SFX('INIT_ISBAN: INCONSISTENCY BETWEEN CPHOTO AND LAGRIP')
ENDIF
IF ( I%CPHOTO/='NCB' .AND. I%CRESPSL=='CNT') THEN
  CALL ABOR1_SFX('INIT_ISBAN: INCONSISTENCY BETWEEN CPHOTO AND CRESPSL')
ENDIF
IF (HINIT=='PRE' .AND. ISIZE_LMEB_PATCH>0 .AND. I%TSNOW%SCHEME.NE.'3-L' .AND. I%TSNOW%SCHEME.NE.'CRO') THEN
    CALL ABOR1_SFX("INIT_ISBAN: WITH LMEB_PATCH = TRUE, CSNOW MUST BE 3-L OR CRO")
ENDIF
IF(I%CPHOTO/='NCB'.AND.I%LSPINUPCARBW)THEN
  CALL ABOR1_SFX('INIT_ISBAN: INCONSISTENCY BETWEEN CPHOTO AND LSPINUPCARBW (if not NCB must be false)')
ENDIF
IF(I%CRESPSL/='CNT'.AND.I%LSPINUPCARBS)THEN
  CALL ABOR1_SFX('INIT_ISBAN: INCONSISTENCY BETWEEN CRESPSL AND LSPINUPCARBS (if not CNT must be false)')
ENDIF
IF(I%LSPINUPCARBW.AND.REAL(I%NNBYEARSPINW)>REAL(I%NNBYEARSPINS)*0.5)THEN
  WRITE(ILUOUT,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  WRITE(ILUOUT,*)'INIT_ISBAN: INCONSISTENCY BETWEEN NNBYEARSPINW AND NNBYEARSPINS'
  WRITE(ILUOUT,*)'NNBYEARSPINW MUST BE < TO 0.5 * NNBYEARSPINS'
  WRITE(ILUOUT,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  CALL ABOR1_SFX('INIT_ISBAN: INCONSISTENCY BETWEEN NNBYEARSPINW AND NNBYEARSPINS')
ENDIF
IF(I%LSPINUPCARBS.AND.(I%XCO2_START==XUNDEF.OR.I%XCO2_END==XUNDEF))THEN
  WRITE(ILUOUT,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  WRITE(ILUOUT,*)'INIT_ISBAN: INCONSISTENCY BETWEEN LSPINUPCARBS AND XCO2_START OR XCO2_END'
  WRITE(ILUOUT,*)'FOR ISBA-CC SPINUP XCO2_START AND XCO2_END MUST BE DEFINED'
  WRITE(ILUOUT,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  CALL ABOR1_SFX('INIT_ISBAN: INCONSISTENCY BETWEEN LSPINUPCARBS AND XCO2_START OR XCO2_END')
ENDIF
!
 CALL END_IO_SURF_n(HPROGRAM)
 CALL SET_SURFEX_FILEIN(HPROGRAM,'PREP') ! restore input file name
!
!-------------------------------------------------------------------------------
!
! During soil carbon spinup with ISBA-CC: 
!        (1) grass parameters are attributed to all agricultural PFT with atmospheric CO2 concentration 
!            fixed to Pre-industrial CO2 consentration XCO2_START
!        (2) Atmospheric CO2 concentration rampin up from XCO2_START to XCO2_END
!
ISPINEND=I%NNBYEARSPINS-NINT(I%NNBYEARSPINS*XSPIN_CO2)
!
I%LAGRI_TO_GRASS = .FALSE.
!
IF ( I%LSPINUPCARBS .AND. (I%NNBYEARSOLD <= ISPINEND) ) THEN
!
   I%LAGRI_TO_GRASS = .TRUE.
!
   CALL INI_DATA_PARAM(DTCO%XDATA_VEGTYPE, PSURF=DTCO%XDATA_NATURE, PH_TREE=XDATA_H_TREE,PLAI=XDATA_LAI,                   &
                                  PALBNIR_VEG=XDATA_ALBNIR_VEG, PALBVIS_VEG=XDATA_ALBVIS_VEG,                    &
                                  PALBUV_VEG=XDATA_ALBUV_VEG, PRSMIN=XDATA_RSMIN,                                &
                                  PRGL=XDATA_RGL, PCV=XDATA_CV, PGAMMA=XDATA_GAMMA,                              &
                                  PGMES=XDATA_GMES, PGC=XDATA_GC, PBSLAI=XDATA_BSLAI,                            &
                                  PSEFOLD=XDATA_SEFOLD, PLAIMIN_OUT=XDATA_LAIMIN, PDMAX=XDATA_DMAX,              &
                                  PSTRESS=XDATA_STRESS, PF2I=XDATA_F2I, PVEG_OUT=XDATA_VEG,                      &
                                  PGREEN=XDATA_GREEN, PZ0=XDATA_Z0, PZ0_O_Z0H=XDATA_Z0_O_Z0H,                    &
                                  PEMIS_ECO=XDATA_EMIS_ECO, PWRMAX_CF=XDATA_WRMAX_CF,                            &
                                  PROOT_LIN=XDATA_ROOT_LIN, PROOT_EXTINCTION=XDATA_ROOT_EXTINCTION,              &
                                  PSOILRC_SO2=XDATA_SOILRC_SO2, PSOILRC_O3=XDATA_SOILRC_O3, PRE25=XDATA_RE25,    &
                                  PCE_NITRO=XDATA_CE_NITRO,PCF_NITRO=XDATA_CF_NITRO,PCNA_NITRO=XDATA_CNA_NITRO,  &
                                  PGMES_ST=XDATA_GMES_ST, PGC_ST=XDATA_GC_ST, PBSLAI_ST=XDATA_BSLAI_ST,          &
                                  PSEFOLD_ST=XDATA_SEFOLD_ST, PDMAX_ST=XDATA_DMAX_ST, OAGRI_TO_GRASS=I%LAGRI_TO_GRASS)
!
   ZCO2(:) = PRHOA(:) * I%XCO2_START * 1.E-6 * XMCO2 / XMD
!
ELSEIF(I%LSPINUPCARBS .AND. (I%NNBYEARSOLD > ISPINEND) .AND. (I%NNBYEARSOLD <= I%NNBYEARSPINS) )THEN
!
   ZSPINCO2 = I%XCO2_START + (I%XCO2_END-I%XCO2_START) * REAL(I%NNBYEARSOLD - ISPINEND) / REAL(I%NNBYEARSPINS - ISPINEND)
!
   ZCO2(:) = PRHOA(:) * ZSPINCO2 * 1.E-6 * XMCO2 / XMD
!
ELSE
!
   ZCO2(:) = PCO2(:)
!
ENDIF
!
!-----------------------------------------------------------------------------------------------------
! END READ PGD FILE
!-----------------------------------------------------------------------------------------------------
!
IF (OLAND_USE .OR. HINIT=='PGD') THEN
  IF (LHOOK) CALL DR_HOOK('INIT_ISBA_N',1,ZHOOK_HANDLE)
  RETURN
END IF
!
CALL COMPUTE_ISBA_PARAMETERS(AG, BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTCO, &
                                    DTI, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, &
                                    DGSI, DGU, DGT, DGUT, DGW, DST, F, FSB, GB, IOB, ICP, &
                                    IG, I, O, S, SSB, SLT, UG, U, SV, TCP, TGD, &
                                    TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                             HPROGRAM,HINIT,OLAND_USE,                  &
                             KI,KSV,KSW,                                &
                             HSV,ZCO2,PRHOA,                            &
                             PZENITH,PSW_BANDS,PDIR_ALB,PSCA_ALB,       &
                             PEMIS,PTSRAD,PTSURF,                       &
                             HTEST                                )
!
IF ( I%CSNOWMETAMO/="B92" ) THEN
  CALL READ_FZ06('drdt_bst_fit_60.nc')
ENDIF
!
IF ( I%CSNOWRAD=="TAR" .OR. I%CSNOWRAD=="TA1" .OR.  I%CSNOWRAD=="TA2" ) THEN
  CALL INIT_TARTES()
END IF
!
IF (LHOOK) CALL DR_HOOK('INIT_ISBA_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE INIT_ISBA_n
