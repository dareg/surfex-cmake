!     #############################################################
      SUBROUTINE INIT_SEAFLUX_n (BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTCO, DTS, &
                                 DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, &
                                 DGT, DGUT, DGW, F, FSB, GB, IOB, ICP, I, O, OR, &
                                 SG, S, SSB, UG, U, SV, TCP, TGD, TGDO, TGR, TGRO, &
                                 T, TOP, TVG, W, WSB, &
                                 HPROGRAM,HINIT,                            &
                                  KI,KSV,KSW,                                &
                                  HSV,PCO2,PRHOA,                            &
                                  PZENITH,PAZIM,PSW_BANDS,PDIR_ALB,PSCA_ALB, &
                                  PEMIS,PTSRAD,PTSURF,                       &
                                  KYEAR, KMONTH,KDAY, PTIME,                 &
                                  HATMFILE,HATMFILETYPE,                     &
                                  HTEST                                      )  
!     #############################################################
!
!!****  *INIT_SEAFLUX_n* - routine to initialize SEAFLUX
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
!!      Original    01/2003
!!      Modified    01/2006 : sea flux parameterization.
!!                  01/2008 : coupling with 1D ocean
!!      B. Decharme 08/2009 : specific treatment for sea/ice in the Earth System Model 
!!      B. Decharme 07/2011 : read pgd+prep 
!!      B. Decharme 04/2013 : new coupling variables
!!      S. Senesi   01/2014 : introduce sea-ice model 
!!      S. Belamari 03/2014 : add NZ0 (to choose PZ0SEA formulation)
!!      R. Séférian 01/2015 : introduce interactive ocean surface albedo
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
!
!
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
USE MODD_FLAKE_n, ONLY : FLAKE_t
USE MODD_FLAKE_SBL_n, ONLY : FLAKE_SBL_t
USE MODD_GR_BIOG_n, ONLY : GR_BIOG_t
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
USE MODD_ISBA_CANOPY_n, ONLY : ISBA_CANOPY_t
USE MODD_ISBA_n, ONLY : ISBA_t
USE MODD_OCEAN_n, ONLY : OCEAN_t
USE MODD_OCEAN_REL_n, ONLY : OCEAN_REL_t
USE MODD_SEAFLUX_GRID_n, ONLY : SEAFLUX_GRID_t
USE MODD_SEAFLUX_n, ONLY : SEAFLUX_t
USE MODD_SEAFLUX_SBL_n, ONLY : SEAFLUX_SBL_t
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
USE MODD_SFX_OASIS,      ONLY : LCPL_SEA, LCPL_SEAICE
!
USE MODD_READ_NAMELIST,  ONLY : LNAM_READ
USE MODD_CSTS,           ONLY : XTTS
USE MODD_SNOW_PAR,       ONLY : XZ0HSN
USE MODD_SURF_PAR,       ONLY : XUNDEF, NUNDEF
USE MODD_CHS_AEROSOL,    ONLY: LVARSIGI, LVARSIGJ
USE MODD_DST_SURF,       ONLY: LVARSIG_DST, NDSTMDE, NDST_MDEBEG, LRGFIX_DST
USE MODD_SLT_SURF,       ONLY: LVARSIG_SLT, NSLTMDE, NSLT_MDEBEG, LRGFIX_SLT
!
USE MODI_INIT_IO_SURF_n
USE MODI_DEFAULT_CH_DEP
!
USE MODI_DEFAULT_SEAFLUX
USE MODI_DEFAULT_DIAG_SEAFLUX
USE MODI_READ_DEFAULT_SEAFLUX_n
USE MODI_READ_SEAFLUX_CONF_n
USE MODI_READ_SEAFLUX_n
!
USE MODI_READ_OCEAN_n
!
USE MODI_DEFAULT_SEAICE
USE MODI_READ_SEAICE_n
USE MODI_GLTOOLS_READNAM
!
USE MODI_READ_PGD_SEAFLUX_n
USE MODI_DIAG_SEAFLUX_INIT_n
USE MODI_END_IO_SURF_n
USE MODI_GET_LUOUT
USE MODI_READ_SURF
USE MODI_READ_SEAFLUX_DATE
USE MODI_READ_NAM_PREP_SEAFLUX_n
USE MODI_INIT_CHEMICAL_n
USE MODI_PREP_CTRL_SEAFLUX
USE MODI_UPDATE_RAD_SEA
USE MODI_READ_SEAFLUX_SBL_n
USE MODI_ABOR1_SFX
!
USE MODI_SET_SURFEX_FILEIN
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
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
TYPE(FLAKE_t), INTENT(INOUT) :: F
TYPE(FLAKE_SBL_t), INTENT(INOUT) :: FSB
TYPE(GR_BIOG_t), INTENT(INOUT) :: GB
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
TYPE(ISBA_CANOPY_t), INTENT(INOUT) :: ICP
TYPE(ISBA_t), INTENT(INOUT) :: I
TYPE(OCEAN_t), INTENT(INOUT) :: O
TYPE(OCEAN_REL_t), INTENT(INOUT) :: OR
TYPE(SEAFLUX_GRID_t), INTENT(INOUT) :: SG
TYPE(SEAFLUX_t), INTENT(INOUT) :: S
TYPE(SEAFLUX_SBL_t), INTENT(INOUT) :: SSB
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
INTEGER           :: ILU    ! sizes of SEAFLUX arrays
INTEGER           :: ILUOUT ! unit of output listing file
INTEGER           :: IRESP  ! return code
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!         Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('INIT_SEAFLUX_N',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('INIT_SEAFLUXN: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF
!
!
!         Others litlle things
!
PDIR_ALB = XUNDEF
PSCA_ALB = XUNDEF
PEMIS    = XUNDEF
PTSRAD   = XUNDEF
PTSURF   = XUNDEF
!
O%LMERCATOR = .FALSE.
O%LCURRENT  = .FALSE.
!
IF (LNAM_READ) THEN
 !
 !*       0.     Defaults
 !               --------
 !
 !        0.1. Hard defaults
 !      
 
 CALL DEFAULT_SEAFLUX(S%XTSTEP,S%XOUT_TSTEP,S%CSEA_ALB,S%CSEA_FLUX,S%LPWG, &
                        S%LPRECIP,S%LPWEBB,S%NZ0,S%NGRVWAVES,O%LPROGSST,   &
                        O%NTIME_COUPLING,S%XICHCE,S%CINTERPOL_SST,     &
                        S%CINTERPOL_SSS                            )
 CALL DEFAULT_SEAICE(HPROGRAM,                                   &
                     S%CINTERPOL_SIC,S%CINTERPOL_SIT, S%XFREEZING_SST, &
                     S%XSEAICE_TSTEP, S%XSIC_EFOLDING_TIME,          &
                     S%XSIT_EFOLDING_TIME, S%XCD_ICE_CST, S%XSI_FLX_DRV)     
 !                     
 CALL DEFAULT_CH_DEP(CHS%CCH_DRY_DEP) 
 !            
 CALL DEFAULT_DIAG_SEAFLUX(DGS%N2M,DGS%LSURF_BUDGET,DGS%L2M_MIN_ZS,DGS%LRAD_BUDGET,DGS%LCOEF,DGS%LSURF_VARS,&
                           DGO%LDIAG_OCEAN,DGSI%LDIAG_SEAICE,DGS%LSURF_BUDGETC,DGS%LRESET_BUDGETC,DGS%XDIAG_TSTEP )  

ENDIF
!
!
!        0.2. Defaults from file header
!    
 CALL READ_DEFAULT_SEAFLUX_n(CHS, DGO, DGS, DGSI, O, S, &
                             HPROGRAM)
!
!*       1.1    Reading of configuration:
!               -------------------------
!
 CALL READ_SEAFLUX_CONF_n(CHS, DGO, DGS, DGSI, O, S, &
                          HPROGRAM)
!
S%LINTERPOL_SST=.FALSE.
S%LINTERPOL_SSS=.FALSE.
S%LINTERPOL_SIC=.FALSE.
S%LINTERPOL_SIT=.FALSE.
IF(LCPL_SEA)THEN 
  IF(DGS%N2M<1)THEN
     CALL ABOR1_SFX('INIT_SEAFLUX_n: N2M must be set >0 in case of LCPL_SEA')
  ENDIF
! No STT / SSS interpolation in Earth System Model
  S%CINTERPOL_SST='NONE  '
  S%CINTERPOL_SSS='NONE  '
  S%CINTERPOL_SIC='NONE  '
  S%CINTERPOL_SIT='NONE  '
ELSE
   IF(TRIM(S%CINTERPOL_SST)/='NONE')THEN
      S%LINTERPOL_SST=.TRUE.
   ENDIF
   IF(TRIM(S%CINTERPOL_SSS)/='NONE')THEN
      S%LINTERPOL_SSS=.TRUE.
   ENDIF
   IF(TRIM(S%CINTERPOL_SIC)/='NONE')THEN
      S%LINTERPOL_SIC=.TRUE.
   ENDIF
   IF(TRIM(S%CINTERPOL_SIT)/='NONE')THEN
      S%LINTERPOL_SIT=.TRUE.
   ENDIF
ENDIF
!
!*       1.     Cover fields and grid:
!               ---------------------
!* date
!
SELECT CASE (HINIT)
!
  CASE ('PGD')
!
    S%TTIME%TDATE%YEAR = NUNDEF
    S%TTIME%TDATE%MONTH= NUNDEF
    S%TTIME%TDATE%DAY  = NUNDEF
    S%TTIME%TIME       = XUNDEF
!
  CASE ('PRE')
!
    CALL PREP_CTRL_SEAFLUX(DGS%N2M,DGS%LSURF_BUDGET,DGS%L2M_MIN_ZS,DGS%LRAD_BUDGET,DGS%LCOEF,DGS%LSURF_VARS,&
                             DGO%LDIAG_OCEAN,DGSI%LDIAG_SEAICE,ILUOUT,DGS%LSURF_BUDGETC ) 
    IF (LNAM_READ) CALL READ_NAM_PREP_SEAFLUX_n(HPROGRAM)      
    CALL READ_SEAFLUX_DATE(IOB, O, &
                           HPROGRAM,HINIT,ILUOUT,HATMFILE,HATMFILETYPE,KYEAR,KMONTH,KDAY,PTIME,S%TTIME)
!
  CASE DEFAULT
!
    CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                        HPROGRAM,'SEA   ','SEAFLX','READ ')
    CALL READ_SURF(IOB, &
                   HPROGRAM,'DTCUR',S%TTIME,IRESP)
    CALL END_IO_SURF_n(HPROGRAM)
!
END SELECT
!
!-----------------------------------------------------------------------------------------------------
! READ PGD FILE
!-----------------------------------------------------------------------------------------------------
!
!         Initialisation for IO
!
 CALL SET_SURFEX_FILEIN(HPROGRAM,'PGD ') ! change input file name to pgd name
 CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                        HPROGRAM,'SEA   ','SEAFLX','READ ')
!
!         Reading of the fields
!
 CALL READ_PGD_SEAFLUX_n(DTCO, DTS, IOB, SG, S, U, &
                         HPROGRAM)
!
 CALL END_IO_SURF_n(HPROGRAM)
 CALL SET_SURFEX_FILEIN(HPROGRAM,'PREP') ! restore input file name
!-------------------------------------------------------------------------------
!
!* if only physiographic fields are to be initialized, stop here.
!
IF (HINIT/='ALL' .AND. HINIT/='SOD') THEN
  IF (LHOOK) CALL DR_HOOK('INIT_SEAFLUX_N',1,ZHOOK_HANDLE)
  RETURN
END IF
!
!-------------------------------------------------------------------------------
!
!         Initialisation for IO
!
 CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                        HPROGRAM,'SEA   ','SEAFLX','READ ')
!
!*       2.     Prognostic fields:
!               ----------------
!
 CALL READ_SEAFLUX_n(DTCO, IOB, SG, S, U, &
                     HPROGRAM,ILUOUT)
!
IF (HINIT/='ALL') THEN
  CALL END_IO_SURF_n(HPROGRAM)
  IF (LHOOK) CALL DR_HOOK('INIT_SEAFLUX_N',1,ZHOOK_HANDLE)
  RETURN
END IF
!-------------------------------------------------------------------------------
!
!*       2.1    Ocean fields:
!               -------------
!
 CALL READ_OCEAN_n(DTCO, IOB, O, OR, U, &
                   HPROGRAM)
!
!-------------------------------------------------------------------------------
!
ILU = SIZE(S%XCOVER,1)
!
ALLOCATE(S%XSST_INI    (ILU))
S%XSST_INI(:) = S%XSST(:)
!
ALLOCATE(S%XZ0H(ILU))
WHERE (S%XSST(:)>=XTTS)
  S%XZ0H(:) = S%XZ0(:)
ELSEWHERE
  S%XZ0H(:) = XZ0HSN
ENDWHERE
!
!-------------------------------------------------------------------------------
!
!*       3.     Specific fields when using earth system model or sea-ice scheme
!               (Sea current and Sea-ice temperature)
!               -----------------------------------------------------------------
!
IF(LCPL_SEA.OR.S%LHANDLE_SIC)THEN       
! 
  ALLOCATE(S%XUMER   (ILU))
  ALLOCATE(S%XVMER   (ILU))
!
  S%XUMER   (:)=XUNDEF
  S%XVMER   (:)=XUNDEF
!
  IF(LCPL_SEAICE.OR.S%LHANDLE_SIC)THEN       
    ALLOCATE(S%XTICE   (ILU))
    ALLOCATE(S%XICE_ALB(ILU))
    S%XTICE   (:)=XUNDEF
    S%XICE_ALB(:)=XUNDEF
  ENDIF
!
ELSE
! 
  ALLOCATE(S%XTICE   (0))
  ALLOCATE(S%XICE_ALB(0))
  ALLOCATE(S%XUMER   (0))
  ALLOCATE(S%XVMER   (0))
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       4.     Seaice prognostic variables and forcings :
!
CALL READ_SEAICE_n(IOB, &
                   SG, S, &
                   HPROGRAM,ILU,ILUOUT)
!
!-------------------------------------------------------------------------------
!
!*       5.     Albedo, emissivity and temperature fields on the mix (open sea + sea ice)
!               -----------------------------------------------------------------
!
ALLOCATE(S%XEMIS    (ILU))
S%XEMIS    = 0.0
!
CALL UPDATE_RAD_SEA(S%CSEA_ALB,S%XSST,PZENITH,XTTS,S%XEMIS,S%XDIR_ALB,&
                    S%XSCA_ALB,PDIR_ALB,PSCA_ALB,PEMIS,PTSRAD,  &
                    S%LHANDLE_SIC,S%XTICE,S%XSIC,S%XICE_ALB           )  
!
IF (S%LHANDLE_SIC) THEN
   PTSURF(:) = S%XSST(:) * ( 1 - S%XSIC(:)) + S%XTICE(:) * S%XSIC(:)
ELSE
   PTSURF(:) = S%XSST(:)
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       6.     SBL air fields:
!               --------------
!
 CALL READ_SEAFLUX_SBL_n(DTCO, IOB, S, SSB, U, &
                         HPROGRAM)
!
!-------------------------------------------------------------------------------
!
!*       7.     Chemistry /dust
!               ---------
!
 CALL INIT_CHEMICAL_n(ILUOUT, KSV, HSV, CHS%NBEQ, CHS%CSV, CHS%NAEREQ,            &
                     CHS%NSV_CHSBEG, CHS%NSV_CHSEND, CHS%NSV_AERBEG, CHS%NSV_AEREND, &
                     CHS%CCH_NAMES, CHS%CAER_NAMES, CHS%NDSTEQ, CHS%NSV_DSTBEG,      &
                     CHS%NSV_DSTEND, CHS%NSLTEQ, CHS%NSV_SLTBEG, CHS%NSV_SLTEND,     &
                     HDSTNAMES=CHS%CDSTNAMES, HSLTNAMES=CHS%CSLTNAMES        )
!
!* deposition scheme
!
IF (CHS%NBEQ>0 .AND. CHS%CCH_DRY_DEP=='WES89') THEN
  ALLOCATE(CHS%XDEP(ILU,CHS%NBEQ))
ELSE
  ALLOCATE(CHS%XDEP(0,0))
END IF
!
!-------------------------------------------------------------------------------
!
!*       8.     diagnostics initialization
!               --------------------------
!
IF(.NOT.(S%LHANDLE_SIC.OR.LCPL_SEAICE))THEN
  DGSI%LDIAG_SEAICE=.FALSE.
ENDIF
!
CALL DIAG_SEAFLUX_INIT_n(IOB, &
                         DGO, DGS, DGSI, DGU, S, &
                         HPROGRAM,ILU,KSW)
!
!-------------------------------------------------------------------------------
!
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('INIT_SEAFLUX_N',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE INIT_SEAFLUX_n
