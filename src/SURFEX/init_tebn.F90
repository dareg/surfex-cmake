!     #############################################################
      SUBROUTINE INIT_TEB_n (B, BOP, BDD, CHI, CHT, DTB, DTCO, DTI, DTGD, DTGR, &
                             DTT, DGCT, DGMT, DGMTO, DGU, DGT, DGUT, DST, GBGD, &
                             GBGR, I, IOB, SLT, UG, U, TCP, TGD, TGDO, TGDPE, TGDP, &
                             TGR, TGRO, TGRPE, TGRP, TG, TIR, T, TOP, TPN, TVG, &
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
!!****  *INIT_TEB_n* - routine to initialize TEB
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
!!      G. Pigeon   09/2012: add ROUGH_WALL/ROUGH_ROOF/CH_BEM for conv. coef.
!!      B. Decharme  04/2013 new coupling variables
!!                           delete CTOPREG option (never used)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_BEM_n, ONLY : BEM_t
USE MODD_BEM_OPTION_n, ONLY : BEM_OPTIONS_t
USE MODD_BLD_DESCRIPTION_n, ONLY : BLD_DESC_t
USE MODD_CH_ISBA_n, ONLY : CH_ISBA_t
USE MODD_CH_TEB_n, ONLY : CH_TEB_t
USE MODD_DATA_BEM_n, ONLY : DATA_BEM_t
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_DATA_TEB_GARDEN_n, ONLY : DATA_TEB_GARDEN_t
USE MODD_DATA_TEB_GREENROOF_n, ONLY : DATA_TEB_GREENROOF_t
USE MODD_DATA_TEB_n, ONLY : DATA_TEB_t
USE MODD_DIAG_CUMUL_TEB_n, ONLY : DIAG_CUMUL_TEB_t
USE MODD_DIAG_MISC_TEB_n, ONLY : DIAG_MISC_TEB_t
USE MODD_DIAG_MISC_TEB_OPTION_n, ONLY : DIAG_MISC_TEB_OPTIONS_t
USE MODD_DIAG_SURF_ATM_n, ONLY : DIAG_SURF_ATM_t
USE MODD_DIAG_TEB_n, ONLY : DIAG_TEB_t
USE MODD_DIAG_UTCI_TEB_n, ONLY : DIAG_UTCI_TEB_t
USE MODD_DST_n, ONLY : DST_t
USE MODD_GR_BIOG_GARDEN_n, ONLY : GR_BIOG_GARDEN_t
USE MODD_GR_BIOG_GREENROOF_n, ONLY : GR_BIOG_GREENROOF_t
USE MODD_ISBA_n, ONLY : ISBA_t
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
USE MODD_SLT_n, ONLY : SLT_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_TEB_CANOPY_n, ONLY : TEB_CANOPY_t
USE MODD_TEB_GARDEN_n, ONLY : TEB_GARDEN_t
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TEB_GARDEN_OPTIONS_t
USE MODD_TEB_GARDEN_PGD_EVOL_n, ONLY : TEB_GARDEN_PGD_EVOL_t
USE MODD_TEB_GARDEN_PGD_n, ONLY : TEB_GARDEN_PGD_t
USE MODD_TEB_GREENROOF_n, ONLY : TEB_GREENROOF_t
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TEB_GREENROOF_OPTIONS_t
USE MODD_TEB_GREENROOF_PGD_EVOL_n, ONLY : TEB_GREENROOF_PGD_EVOL_t
USE MODD_TEB_GREENROOF_PGD_n, ONLY : TEB_GREENROOF_PGD_t
USE MODD_TEB_GRID_n, ONLY : TEB_GRID_t
USE MODD_TEB_IRRIG_n, ONLY : TEB_IRRIG_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_PANEL_n, ONLY : TEB_PANEL_t
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_OPTIONS_t
!
USE MODD_IO_SURF_ASC,ONLY: CMASK
USE MODD_SNOW_PAR, ONLY : XEMISSN
!
USE MODD_READ_NAMELIST, ONLY : LNAM_READ


USE MODD_CHS_AEROSOL,     ONLY: LVARSIGI, LVARSIGJ
USE MODD_DST_SURF,        ONLY: LVARSIG_DST, NDSTMDE, NDST_MDEBEG, LRGFIX_DST 
USE MODD_SLT_SURF,        ONLY: LVARSIG_SLT, NSLTMDE, NSLT_MDEBEG, LRGFIX_SLT
USE MODD_SURF_PAR,        ONLY: XUNDEF, NUNDEF
!
USE MODI_INIT_IO_SURF_n
USE MODI_DEFAULT_CH_DEP
USE MODI_DEFAULT_TEB
USE MODI_DEFAULT_DIAG_TEB
USE MODI_READ_DEFAULT_TEB_n
USE MODI_READ_TEB_CONF_n
USE MODI_PREP_CTRL_TEB
USE MODI_READ_TEB_n
USE MODI_READ_PGD_TEB_n
USE MODI_CONVERT_TEB
USE MODI_CONVERT_PATCH_TEB
USE MODI_INIT_SNOW_LW
USE MODI_AVERAGED_TSRAD_TEB
USE MODI_AVERAGED_ALBEDO_TEB
USE MODI_DIAG_TEB_INIT_n
USE MODI_DIAG_MISC_TEB_INIT_n
USE MODI_END_IO_SURF_n
USE MODI_GET_LUOUT
USE MODI_READ_SURF
USE MODI_READ_PREP_TEB_SNOW
USE MODI_READ_TEB_DATE
USE MODI_READ_NAM_PREP_TEB_n
USE MODI_INIT_CHEMICAL_n
USE MODI_GARDEN_PROPERTIES
USE MODI_HVAC_AUTOSIZE
USE MODI_GOTO_WRAPPER_TEB_PATCH
!
USE MODI_INIT_TEB_GARDEN_n
USE MODI_INIT_TEB_GARDEN_PGD_n
USE MODI_INIT_TEB_VEG_OPTIONS_n
USE MODI_TEB_MORPHO
USE MODI_INIT_BEM_n
USE MODI_INIT_TEB_GREENROOF_n
USE MODI_INIT_TEB_GREENROOF_PGD_n
USE MODI_GREENROOF_PROPERTIES
USE MODI_READ_PGD_TEB_IRRIG_n
!
USE MODI_READ_COVER_GARDEN
USE MODI_ABOR1_SFX
USE MODI_READ_TEB_CANOPY_n
USE MODI_SET_SURFEX_FILEIN
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
TYPE(BEM_t), INTENT(INOUT) :: B
TYPE(BEM_OPTIONS_t), INTENT(INOUT) :: BOP
TYPE(BLD_DESC_t), INTENT(INOUT) :: BDD
TYPE(CH_ISBA_t), INTENT(INOUT) :: CHI
TYPE(CH_TEB_t), INTENT(INOUT) :: CHT
TYPE(DATA_BEM_t), INTENT(INOUT) :: DTB
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTI
TYPE(DATA_TEB_GARDEN_t), INTENT(INOUT) :: DTGD
TYPE(DATA_TEB_GREENROOF_t), INTENT(INOUT) :: DTGR
TYPE(DATA_TEB_t), INTENT(INOUT) :: DTT
TYPE(DIAG_CUMUL_TEB_t), INTENT(INOUT) :: DGCT
TYPE(DIAG_MISC_TEB_t), INTENT(INOUT) :: DGMT
TYPE(DIAG_MISC_TEB_OPTIONS_t), INTENT(INOUT) :: DGMTO
TYPE(DIAG_SURF_ATM_t), INTENT(INOUT) :: DGU
TYPE(DIAG_TEB_t), INTENT(INOUT) :: DGT
TYPE(DIAG_UTCI_TEB_t), INTENT(INOUT) :: DGUT
TYPE(DST_t), INTENT(INOUT) :: DST
TYPE(GR_BIOG_GARDEN_t), INTENT(INOUT) :: GBGD
TYPE(GR_BIOG_GREENROOF_t), INTENT(INOUT) :: GBGR
TYPE(ISBA_t), INTENT(INOUT) :: I
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
TYPE(SLT_t), INTENT(INOUT) :: SLT
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(TEB_CANOPY_t), INTENT(INOUT) :: TCP
TYPE(TEB_GARDEN_t), INTENT(INOUT) :: TGD
TYPE(TEB_GARDEN_OPTIONS_t), INTENT(INOUT) :: TGDO
TYPE(TEB_GARDEN_PGD_EVOL_t), INTENT(INOUT) :: TGDPE
TYPE(TEB_GARDEN_PGD_t), INTENT(INOUT) :: TGDP
TYPE(TEB_GREENROOF_t), INTENT(INOUT) :: TGR
TYPE(TEB_GREENROOF_OPTIONS_t), INTENT(INOUT) :: TGRO
TYPE(TEB_GREENROOF_PGD_EVOL_t), INTENT(INOUT) :: TGRPE
TYPE(TEB_GREENROOF_PGD_t), INTENT(INOUT) :: TGRP
TYPE(TEB_GRID_t), INTENT(INOUT) :: TG
TYPE(TEB_IRRIG_t), INTENT(INOUT) :: TIR
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_PANEL_t), INTENT(INOUT) :: TPN
TYPE(TEB_VEG_OPTIONS_t), INTENT(INOUT) :: TVG
!
 CHARACTER(LEN=6),                   INTENT(IN)  :: HPROGRAM    ! program calling surf. schemes
 CHARACTER(LEN=3),                   INTENT(IN)  :: HINIT       ! choice of fields to initialize
INTEGER,                            INTENT(IN)  :: KI          ! number of points
INTEGER,                            INTENT(IN)  :: KSV         ! number of scalars
INTEGER,                            INTENT(IN)  :: KSW         ! number of short-wave spectral bands
 CHARACTER(LEN=6), DIMENSION(KSV),   INTENT(IN)  :: HSV         ! name of all scalar variables
REAL,             DIMENSION(KI),    INTENT(IN)  :: PCO2        ! CO2 concentration (kg/m3)
REAL,             DIMENSION(KI),    INTENT(IN)  :: PRHOA       ! air density
REAL,             DIMENSION(KI),    INTENT(IN)  :: PZENITH     ! solar zenithal angle
REAL,             DIMENSION(KI),    INTENT(IN)  :: PAZIM       ! solar azimuthal angle (rad from N, clock)
REAL,             DIMENSION(KSW),   INTENT(IN)  :: PSW_BANDS   ! middle wavelength of each band
REAL,             DIMENSION(KI,KSW),INTENT(OUT) :: PDIR_ALB    ! direct albedo for each band
REAL,             DIMENSION(KI,KSW),INTENT(OUT) :: PSCA_ALB    ! diffuse albedo for each band
REAL,             DIMENSION(KI),    INTENT(OUT) :: PEMIS       ! emissivity
REAL,             DIMENSION(KI),    INTENT(OUT) :: PTSRAD      ! radiative temperature
REAL,             DIMENSION(KI),    INTENT(OUT) :: PTSURF      ! surface effective temperature         (K)
INTEGER,                            INTENT(IN)  :: KYEAR       ! current year (UTC)
INTEGER,                            INTENT(IN)  :: KMONTH      ! current month (UTC)
INTEGER,                            INTENT(IN)  :: KDAY        ! current day (UTC)
REAL,                               INTENT(IN)  :: PTIME       ! current time since
                                                               !  midnight (UTC, s)
!
 CHARACTER(LEN=28),                  INTENT(IN)  :: HATMFILE    ! atmospheric file name
 CHARACTER(LEN=6),                   INTENT(IN)  :: HATMFILETYPE! atmospheric file type
 CHARACTER(LEN=2),                   INTENT(IN)  :: HTEST       ! must be equal to 'OK'
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER                         :: ILU              ! sizes of TEB arrays
INTEGER                         :: ILUOUT           ! unit of output listing file
INTEGER                         :: IRESP            ! return code
!
INTEGER                         :: ISWB             ! number of shortwave spectral bands
INTEGER                         :: JSWB             ! loop on shortwave spectral bands
!
REAL                            :: ZDEF_ROAD_DIR    ! default raod direction
REAL, DIMENSION(:), ALLOCATABLE :: ZDIR_ALB         ! direct town albedo
REAL, DIMENSION(:), ALLOCATABLE :: ZSCA_ALB         ! diffuse town albedo
!
!              local variables for urban green areas
REAL, DIMENSION(KI,KSW)         :: ZDIR_SW          ! direct  SW for each band
REAL, DIMENSION(KI,KSW)         :: ZSCA_SW          ! diffuse SW for each band
REAL, DIMENSION(KI)             :: ZEMIS_GARDEN     ! emissivity
REAL, DIMENSION(KI)             :: ZALB_GARDEN      ! albedo
REAL, DIMENSION(KI)             :: ZTS_GARDEN       ! radiative temperature
!
!              local variables for urban greenroofs
REAL, DIMENSION(KI)             :: ZEMIS_GREENROOF     ! emissivity
REAL, DIMENSION(KI)             :: ZALB_GREENROOF      ! albedo
REAL, DIMENSION(KI)             :: ZTS_GREENROOF       ! radiative temperature
!
INTEGER                         :: JPATCH
INTEGER                         :: IVERSION, IBUGFIX

REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!         Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('INIT_TEB_N',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('INIT_TEBN: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF
!
!         Other little things
!
PDIR_ALB = XUNDEF
PSCA_ALB = XUNDEF
PEMIS    = XUNDEF
PTSRAD   = XUNDEF
PTSURF   = XUNDEF
!
DGMTO%LSURF_EVAP_BUDGET = .FALSE.
!
IF (LNAM_READ) THEN
 !
 !*       0.     Defaults
 !               --------
 !
 !        0.1. Hard defaults
 !      
 CALL DEFAULT_TEB(TOP%CZ0H,TOP%XTSTEP,TOP%XOUT_TSTEP, TOP%CCH_BEM, T%CUR%XDT_RES, T%CUR%XDT_OFF)
 CALL DEFAULT_CH_DEP(CHT%CCH_DRY_DEP)
 CALL DEFAULT_DIAG_TEB(DGT%N2M,DGT%LSURF_BUDGET,DGT%L2M_MIN_ZS,DGT%LRAD_BUDGET,DGT%LCOEF,DGT%LSURF_VARS, &
                       DGMTO%LSURF_MISC_BUDGET,DGMTO%LSURF_DIAG_ALBEDO,DGUT%LUTCI,DGT%LPGD,DGT%LPGD_FIX,DGT%XDIAG_TSTEP)  
!
ENDIF
!
!        0.2. Defaults from file header
!    
 CALL READ_DEFAULT_TEB_n(CHT, DGMTO, DGT, DGUT, TGRO, T, TOP, &
                         HPROGRAM)
!
!*       1.     Reading of configuration:
!               -------------------------
!
 CALL READ_TEB_CONF_n(CHT, DGMTO, DGT, DGUT, T, TOP, &
                      HPROGRAM)
!
!* initialization of snow scheme
!
IF (HINIT=='PRE') THEN
  DO JPATCH=1,TOP%NTEB_PATCH
    CALL GOTO_WRAPPER_TEB_PATCH(B, DGCT, DGMT, T, TGD, TGDPE, TGR, TGRPE, JPATCH)
    CALL READ_PREP_TEB_SNOW(HPROGRAM,T%CUR%TSNOW_ROOF%SCHEME,T%CUR%TSNOW_ROOF%NLAYER, &
                                     T%CUR%TSNOW_ROAD%SCHEME,T%CUR%TSNOW_ROAD%NLAYER)
  END DO
ENDIF
!
!*       2.     Cover fields and grid:
!               ---------------------
!* date
!
SELECT CASE (HINIT)
  CASE ('PGD')
    TOP%TTIME%TDATE%YEAR = NUNDEF
    TOP%TTIME%TDATE%MONTH= NUNDEF
    TOP%TTIME%TDATE%DAY  = NUNDEF
    TOP%TTIME%TIME       = XUNDEF

  CASE ('PRE')
    CALL PREP_CTRL_TEB(DGT%N2M,DGT%LSURF_BUDGET,DGT%L2M_MIN_ZS,DGT%LRAD_BUDGET,DGT%LCOEF,DGT%LSURF_VARS,&
                         DGMTO%LSURF_EVAP_BUDGET,DGMTO%LSURF_MISC_BUDGET,DGUT%LUTCI,ILUOUT )           
    IF (LNAM_READ) CALL READ_NAM_PREP_TEB_n(HPROGRAM)   
    CALL READ_TEB_DATE(IOB, &
                       HPROGRAM,HINIT,ILUOUT,HATMFILE,HATMFILETYPE,KYEAR,KMONTH,KDAY,PTIME,TOP%TTIME)

  CASE DEFAULT
CALL INIT_IO_SURF_n(DTCO, DGU, IOB, U, &
                        HPROGRAM,'TOWN  ','TEB   ','READ ')
    CALL READ_SURF(IOB, &
                   HPROGRAM,'DTCUR',TOP%TTIME,IRESP)
    CALL END_IO_SURF_n(HPROGRAM)
END SELECT
!
!-----------------------------------------------------------------------------------------------------
! READ PGD FILE
!-----------------------------------------------------------------------------------------------------
!
!         Initialisation for IO
!
 CALL SET_SURFEX_FILEIN(HPROGRAM,'PGD ') ! change input file name to pgd name
CALL INIT_IO_SURF_n(DTCO, DGU, IOB, U, &
                        HPROGRAM,'TOWN  ','TEB   ','READ ')
!
 CALL READ_SURF(IOB, &
                   HPROGRAM,'VERSION',IVERSION,IRESP)
 CALL READ_SURF(IOB, &
                   HPROGRAM,'BUG',IBUGFIX,IRESP)
!
!         Reading of the fields
!
 CALL READ_COVER_GARDEN(IOB, &
                        HPROGRAM,TOP%LGARDEN)
!
 CALL READ_PGD_TEB_n(BOP, BDD, DTB, DTCO, DTT, IOB, U, TG, TOP, &
                     HPROGRAM)
!
 CALL END_IO_SURF_n(HPROGRAM)
! 
!*        Fraction of each patch in the grid mesh
!
ILU = SIZE(TOP%XCOVER,1)
!
ALLOCATE(TOP%XTEB_PATCH(ILU,TOP%NTEB_PATCH))
 CALL CONVERT_TEB(TOP, &
                  TOP%XCOVER,TOP%XTEB_PATCH)
!
 CALL SET_SURFEX_FILEIN(HPROGRAM,'PREP') ! restore input file name
CALL INIT_IO_SURF_n(DTCO, DGU, IOB, U, &
                        HPROGRAM,'TOWN  ','TEB   ','READ ')
!
 CALL READ_SURF(IOB, &
                   HPROGRAM,'VERSION',IVERSION,IRESP)
 CALL READ_SURF(IOB, &
                   HPROGRAM,'BUG',IBUGFIX,IRESP)
!
!* reads what is the option defined for road orientations & walls
!
IF (HINIT=='ALL') THEN
  TOP%CROAD_DIR='UNIF'
  TOP%CWALL_OPT='UNIF'
  IF (IVERSION>7 .OR. (IVERSION==7 .AND. IBUGFIX>=3)) THEN
    CALL READ_SURF(IOB, &
                   HPROGRAM,'ROAD_DIR',TOP%CROAD_DIR,IRESP)
    CALL READ_SURF(IOB, &
                   HPROGRAM,'WALL_OPT',TOP%CWALL_OPT,IRESP)
  END IF
END IF
 CALL END_IO_SURF_n(HPROGRAM)
!-----------------------------------------------------------------------------------
!
!*              LOOP ON TEB PATCHES
!               -------------------
!
DO JPATCH=1,TOP%NTEB_PATCH

  CALL GOTO_WRAPPER_TEB_PATCH(B, DGCT, DGMT, T, TGD, TGDPE, TGR, TGRPE, JPATCH)
  !-----------------------------------------------------------------------------------
  !
  !*       3.     Physiographic data fields from land cover:
  !               -----------------------------------------
  !
  ALLOCATE(T%CUR%XZ0_TOWN     (ILU))
  ALLOCATE(T%CUR%XALB_ROOF    (ILU))
  ALLOCATE(T%CUR%XEMIS_ROOF   (ILU))
  ALLOCATE(T%CUR%XALB_ROAD    (ILU))
  ALLOCATE(T%CUR%XEMIS_ROAD   (ILU))
  ALLOCATE(T%CUR%XALB_WALL    (ILU))
  ALLOCATE(T%CUR%XEMIS_WALL   (ILU))
  ALLOCATE(T%CUR%XBLD         (ILU))
  ALLOCATE(T%CUR%XROAD_DIR    (ILU))
  ALLOCATE(T%CUR%XROAD        (ILU))
  ALLOCATE(T%CUR%XBLD_HEIGHT  (ILU))
  ALLOCATE(T%CUR%XWALL_O_HOR  (ILU))
  ALLOCATE(T%CUR%XCAN_HW_RATIO(ILU))
  ALLOCATE(T%CUR%XROAD_O_GRND(ILU))
  ALLOCATE(T%CUR%XGARDEN_O_GRND(ILU))
  ALLOCATE(T%CUR%XWALL_O_GRND(ILU))
  ALLOCATE(T%CUR%XWALL_O_BLD(ILU))
  ALLOCATE(T%CUR%XH_TRAFFIC   (ILU))
  ALLOCATE(T%CUR%XLE_TRAFFIC  (ILU))
  ALLOCATE(T%CUR%XH_INDUSTRY  (ILU))
  ALLOCATE(T%CUR%XLE_INDUSTRY (ILU))
  ALLOCATE(T%CUR%XHC_ROOF     (ILU,TOP%NROOF_LAYER))
  ALLOCATE(T%CUR%XTC_ROOF     (ILU,TOP%NROOF_LAYER))
  ALLOCATE(T%CUR%XD_ROOF      (ILU,TOP%NROOF_LAYER))
  ALLOCATE(T%CUR%XHC_ROAD     (ILU,TOP%NROAD_LAYER))
  ALLOCATE(T%CUR%XTC_ROAD     (ILU,TOP%NROAD_LAYER))
  ALLOCATE(T%CUR%XD_ROAD      (ILU,TOP%NROAD_LAYER))
  ALLOCATE(T%CUR%XHC_WALL     (ILU,TOP%NWALL_LAYER))
  ALLOCATE(T%CUR%XTC_WALL     (ILU,TOP%NWALL_LAYER))
  ALLOCATE(T%CUR%XD_WALL      (ILU,TOP%NWALL_LAYER))
  ALLOCATE(T%CUR%XROUGH_ROOF      (ILU))
  ALLOCATE(T%CUR%XROUGH_WALL      (ILU))
  ALLOCATE(T%CUR%XRESIDENTIAL     (ILU))
  ALLOCATE(T%CUR%XGREENROOF       (ILU))
  ALLOCATE(T%CUR%XGARDEN          (ILU))
  ALLOCATE(TPN%XEMIS_PANEL      (ILU))
  ALLOCATE(TPN%XALB_PANEL       (ILU))
  ALLOCATE(TPN%XEFF_PANEL       (ILU))
  ALLOCATE(TPN%XFRAC_PANEL      (ILU))
  !
  T%CUR%XROAD_DIR(:) = 0.
  T%CUR%XROAD    (:) = 0.
  !
  ZDEF_ROAD_DIR = 0.
  IF (TOP%CROAD_DIR/='UNIF') THEN
    !* road direction if not specified by the user depends on patch number
    !  First patch has a Notrh-South road. Other patches have roads spaced by
    !  regular angles
    ZDEF_ROAD_DIR = 180. * FLOAT(JPATCH-1) / FLOAT(TOP%NTEB_PATCH)
  END IF
  !
  CALL CONVERT_PATCH_TEB(BDD, DTB, DTCO, DTT, TOP, &
                         TOP%XCOVER, TOP%LCOVER, ZDEF_ROAD_DIR,                          &
                      PZ0_TOWN=T%CUR%XZ0_TOWN,                                         &
                      PALB_ROOF=T%CUR%XALB_ROOF,                                       &
                      PEMIS_ROOF=T%CUR%XEMIS_ROOF,PHC_ROOF=T%CUR%XHC_ROOF,PTC_ROOF=T%CUR%XTC_ROOF, &
                      PD_ROOF=T%CUR%XD_ROOF,                                           &
                      PALB_ROAD=T%CUR%XALB_ROAD,                                       &
                      PEMIS_ROAD=T%CUR%XEMIS_ROAD,PHC_ROAD=T%CUR%XHC_ROAD,PTC_ROAD=T%CUR%XTC_ROAD, &
                      PD_ROAD=T%CUR%XD_ROAD,                                           &
                      PALB_WALL=T%CUR%XALB_WALL,                                       &
                      PEMIS_WALL=T%CUR%XEMIS_WALL,PHC_WALL=T%CUR%XHC_WALL,PTC_WALL=T%CUR%XTC_WALL, &
                      PD_WALL=T%CUR%XD_WALL,                                           &
                      PBLD_HEIGHT=T%CUR%XBLD_HEIGHT,                                   &
                      PWALL_O_HOR=T%CUR%XWALL_O_HOR,PBLD=T%CUR%XBLD, PROAD_DIR=T%CUR%XROAD_DIR,    &
                      PGARDEN=T%CUR%XGARDEN,                                           &
                      PH_TRAFFIC=T%CUR%XH_TRAFFIC, PLE_TRAFFIC=T%CUR%XLE_TRAFFIC,            &
                      PH_INDUSTRY=T%CUR%XH_INDUSTRY, PLE_INDUSTRY=T%CUR%XLE_INDUSTRY,        &
                      PROUGH_ROOF = T%CUR%XROUGH_ROOF, PROUGH_WALL = T%CUR%XROUGH_WALL,      &
                      PRESIDENTIAL = T%CUR%XRESIDENTIAL,                               &
                      PGREENROOF = T%CUR%XGREENROOF,                                   &
                      PEMIS_PANEL=TPN%XEMIS_PANEL, PALB_PANEL=TPN%XALB_PANEL,            &
                      PEFF_PANEL=TPN%XEFF_PANEL, PFRAC_PANEL=TPN%XFRAC_PANEL             )
  !
  IF (.NOT. TOP%LGREENROOF .AND. MAXVAL(T%CUR%XGREENROOF)>0. ) THEN !<== A paralleliser pour un stop propre
    WRITE(ILUOUT,*) 'You choose NOT to have greenroofs, BUT your greenroof fraction is not zero'
    WRITE(ILUOUT,*) 'Please activate the greenroof option (and rerun the SURFEX suite from the PGD step)'
    WRITE(ILUOUT,*) 'Or be sure NOT to have any greenroofs in your area'
    CALL ABOR1_SFX('INIT_TEBN: GREENROOF OPTION NOT ACTIVATED WHILE GREENROOFS ARE PRESENT')
  ENDIF
  !
  IF (.NOT. TOP%LSOLAR_PANEL .AND. MAXVAL(TPN%XFRAC_PANEL)>0. ) THEN !<== A paralleliser pour un stop propre
    WRITE(ILUOUT,*) 'You choose NOT to have solar panels, BUT your solar panel fraction is not zero'
    WRITE(ILUOUT,*) 'Please activate the solar panel option (and rerun the SURFEX suite from the PGD step)'
    WRITE(ILUOUT,*) 'Or be sure NOT to have any solar panel in your area'
    CALL ABOR1_SFX('INIT_TEBN: SOLAR_PANEL OPTION NOT ACTIVATED WHILE SOLAR PANELS ARE PRESENT')
  ENDIF
  !
  !-------------------------------------------------------------------------------
  !
  !*       5.     Sky-view-factors:
  !               ----------------
  !
  ALLOCATE(T%CUR%XSVF_ROAD  (ILU))
  ALLOCATE(T%CUR%XSVF_GARDEN(ILU))
  ALLOCATE(T%CUR%XSVF_WALL  (ILU))
  !
  ALLOCATE(B%CUR%XGR          (ILU))
  ALLOCATE(B%CUR%XALB_WIN     (ILU))
  ALLOCATE(B%CUR%XF_WASTE_CAN (ILU))
  !
  !
  CALL TEB_MORPHO(HPROGRAM, T%CUR%XBLD, T%CUR%XWALL_O_HOR, T%CUR%XGARDEN, T%CUR%XBLD_HEIGHT, T%CUR%XROAD, T%CUR%XROAD_O_GRND, &
                T%CUR%XGARDEN_O_GRND, T%CUR%XWALL_O_GRND, T%CUR%XCAN_HW_RATIO, T%CUR%XSVF_ROAD, T%CUR%XSVF_GARDEN,    &
                T%CUR%XSVF_WALL, T%CUR%XZ0_TOWN, T%CUR%XWALL_O_BLD, T%CUR%XH_TRAFFIC, T%CUR%XLE_TRAFFIC               )
                !
  !-------------------------------------------------------------------------------
  !
  !*       6.     Building Energy Model
  !               ---------------------
  !
  CALL INIT_BEM_n(DGU, IOB, B, BOP, BDD, DTB, DTCO, DTT, UG, U, TG, T, TOP, &
                  ILUOUT)
  !
  !-------------------------------------------------------------------------------
  !
  !*      7.      Case of urban green areas
  !               -------------------------
  !
  IF (TOP%LGARDEN) THEN
  !
    CALL SET_SURFEX_FILEIN(HPROGRAM,'PGD ') ! change input file name to pgd name
CALL INIT_IO_SURF_n(DTCO, DGU, IOB, U, &
                        HPROGRAM,'TOWN  ','TEB   ','READ ')     
    IF (JPATCH==1) CALL INIT_TEB_VEG_OPTIONS_n(IOB, &
                                               CHT, DGMTO, TGDO, TVG, &
                                               HPROGRAM)
    CALL INIT_TEB_GARDEN_PGD_n(CHI, CHT, DTCO, DTI, DTGD, DST, GBGD, IOB, I, SLT, U, TGDO, TGDPE, TGDP, TG, T, TOP, TVG, &
                               HPROGRAM,HINIT,(JPATCH==1),KI,KSV,HSV,IVERSION,IBUGFIX,PCO2,PRHOA)
    ! Case of urban green roofs
    IF (TOP%LGREENROOF) CALL INIT_TEB_GREENROOF_PGD_n(CHI, CHT, DTCO, DTI, DTGR, DST, GBGR, IOB, I, SLT, U, &
                                                      TGRO, TGRPE, TGRP, TG, T, TOP, TVG, &
                                                      HPROGRAM,HINIT,(JPATCH==1),KI,KSV,HSV,IVERSION,PCO2,PRHOA)
    CALL END_IO_SURF_n(HPROGRAM)
    !
  ENDIF
!-------------------------------------------------------------------------------
END DO ! end of loop on TEB patches
!-------------------------------------------------------------------------------
!
!* Read irrigation parameters for TEB
!
 CALL SET_SURFEX_FILEIN(HPROGRAM,'PGD ') ! change input file name to pgd name
CALL INIT_IO_SURF_n(DTCO, DGU, IOB, U, &
                        HPROGRAM,'TOWN  ','TEB   ','READ ')     
 CALL READ_PGD_TEB_IRRIG_n(IOB, &
                           TG, TIR, &
                           HPROGRAM)
 CALL END_IO_SURF_n(HPROGRAM)
!
!-------------------------------------------------------------------------------
!
!* if only physiographic fields are to be initialized, stop here.
!
IF (HINIT/='ALL' .AND. HINIT/='SOD') THEN
  IF (LHOOK) CALL DR_HOOK('INIT_TEB_N',1,ZHOOK_HANDLE)
  RETURN
END IF
!
!-------------------------------------------------------------------------------
!
!         Initialisation for IO
!
 CALL SET_SURFEX_FILEIN(HPROGRAM,'PREP') ! restore input file name
CALL INIT_IO_SURF_n(DTCO, DGU, IOB, U, &
                        HPROGRAM,'TOWN  ','TEB   ','READ ')
!
!*       9.     Prognostic fields:
!               -----------------
!
!               -------------------------
!

!
!*              LOOP ON TEB PATCHES
!               -------------------
!
DO JPATCH=1,TOP%NTEB_PATCH
  CALL GOTO_WRAPPER_TEB_PATCH(B, DGCT, DGMT, T, TGD, TGDPE, TGR, TGRPE, JPATCH)
!
!* TEB fields
  CALL READ_TEB_n(B, BOP, DTCO, DGU, IOB, U, T, TOP, TPN, &
                  HPROGRAM,JPATCH)
!
  ALLOCATE(T%CUR%XAC_ROOF    (ILU))
  ALLOCATE(T%CUR%XAC_ROAD    (ILU))
  ALLOCATE(T%CUR%XAC_WALL    (ILU))
  ALLOCATE(T%CUR%XAC_TOP     (ILU))
  ALLOCATE(T%CUR%XAC_ROOF_WAT(ILU))
  ALLOCATE(T%CUR%XAC_ROAD_WAT(ILU))
  ALLOCATE(T%CUR%XQSAT_ROOF  (ILU))
  ALLOCATE(T%CUR%XQSAT_ROAD  (ILU))
  ALLOCATE(T%CUR%XDELT_ROOF  (ILU))
  ALLOCATE(T%CUR%XDELT_ROAD  (ILU))
!
!* Case of urban green areas
  IF (TOP%LGARDEN) THEN
!    CALL SET_SURFEX_FILEIN(HPROGRAM,'PREP') ! change input file name to pgd name
!    CALL INIT_IO_SURF_n(HPROGRAM,'TOWN  ','TEB   ','READ ')       
    CALL INIT_TEB_GARDEN_n(DTCO, DTGD, DGMTO, DGU, IOB, UG, &
                           U, TGD, TGDO, TGDPE, TGDP, TOP, TVG, &
                           HPROGRAM,HINIT,KI,KSW,PSW_BANDS)
  ! Case of urban green roofs
    IF (TOP%LGREENROOF) CALL INIT_TEB_GREENROOF_n(DTCO, DTGR, DGMTO, IOB, U, TGR, TGRO, TGRPE, TGRP, TOP, TVG, &
                                                  HPROGRAM,HINIT,KI,KSV,PSW_BANDS)
!    CALL END_IO_SURF_n(HPROGRAM)
  ENDIF
!-------------------------------------------------------------------------------
!
!*      10.     Infra-red Radiative fields:
!               --------------------------
!
!* snow long-wave properties (not initialized in read_gr_snow)
!
  CALL INIT_SNOW_LW(XEMISSN,T%CUR%TSNOW_ROOF)
  CALL INIT_SNOW_LW(XEMISSN,T%CUR%TSNOW_ROAD)
!
  IF (TOP%LGARDEN) THEN
    ZDIR_SW=0. ! night as first guess for albedo computation
    ZSCA_SW=0. !
    CALL GARDEN_PROPERTIES(TGD, TGDO, TGDPE, TGDP, T, TVG, &
                           ZDIR_SW, ZSCA_SW, PSW_BANDS, KSW,     &
                           ZTS_GARDEN, ZEMIS_GARDEN, ZALB_GARDEN )      
  ELSE
    ZALB_GARDEN = XUNDEF
    ZEMIS_GARDEN= XUNDEF
    ZTS_GARDEN  = XUNDEF
  END IF
  !
  IF (TOP%LGREENROOF) THEN
    ZDIR_SW=0. ! night as first guess for albedo computation
    ZSCA_SW=0. !
    CALL GREENROOF_PROPERTIES(TGR, TGRO, TGRPE, TGRP, T, TVG, &
                              ZDIR_SW, ZSCA_SW, PSW_BANDS, KSW,              &
                              ZTS_GREENROOF, ZEMIS_GREENROOF, ZALB_GREENROOF )  
  ELSE
    ZALB_GREENROOF  = XUNDEF
    ZEMIS_GREENROOF = XUNDEF
    ZTS_GREENROOF   = XUNDEF
  END IF
!
!* averaged albedo, emissivity and radiative temperature
!
  CALL AVERAGED_TSRAD_TEB(T%CUR%XEMIS_ROOF,T%CUR%XT_ROOF(:,1),       &
                        T%CUR%XEMIS_ROAD,T%CUR%XT_ROAD(:,1),       &
                        T%CUR%XEMIS_WALL,                    &
                        T%CUR%XT_WALL_A(:,1),                &
                        T%CUR%XT_WALL_B(:,1),                &
                        ZEMIS_GARDEN, ZTS_GARDEN,      &
                        ZEMIS_GREENROOF, ZTS_GREENROOF,&
                        T%CUR%TSNOW_ROOF,T%CUR%TSNOW_ROAD,         &
                        T%CUR%XROAD, T%CUR%XGREENROOF, T%CUR%XGARDEN,    &
                        T%CUR%XBLD,T%CUR%XWALL_O_HOR,              &
                        T%CUR%XSVF_ROAD,T%CUR%XSVF_WALL,           &
                        T%CUR%XSVF_GARDEN,                   &
                        PEMIS,PTSRAD, B%CUR%XT_WIN1,         &
                        B%CUR%XGR                            )
!
!
!*       9.     Visible and near-infra-red Radiative fields:
!               -------------------------------------------
!
  ALLOCATE(ZDIR_ALB(ILU))
  ALLOCATE(ZSCA_ALB(ILU))
!
  CALL AVERAGED_ALBEDO_TEB(TOP%CBEM,TOP%CROAD_DIR,TOP%CWALL_OPT,PZENITH,PAZIM, &
                       T%CUR%XBLD, T%CUR%XGARDEN, T%CUR%XROAD_DIR, T%CUR%XROAD, T%CUR%XGREENROOF,&
                       TPN%XFRAC_PANEL, TPN%XALB_PANEL,                    &
                       T%CUR%XWALL_O_HOR, T%CUR%XCAN_HW_RATIO,                 &
                       T%CUR%XALB_ROOF,                                  &
                       T%CUR%XALB_ROAD, T%CUR%XSVF_ROAD,                       &
                       T%CUR%XALB_WALL, T%CUR%XSVF_WALL,                       &
                       ZALB_GARDEN, T%CUR%XSVF_GARDEN,                   &
                       ZALB_GREENROOF,                             &
                       T%CUR%TSNOW_ROOF, T%CUR%TSNOW_ROAD,                     &
                       B%CUR%XGR, B%CUR%XSHGC, B%CUR%XSHGC_SH, B%CUR%XABS_WIN, B%CUR%XALB_WIN,   &
                       B%CUR%LSHAD_DAY,                                  &
                       ZDIR_ALB, ZSCA_ALB, B%CUR%XTRAN_WIN               )  

  ISWB=SIZE(PSW_BANDS)
  DO JSWB=1,ISWB
    PDIR_ALB(:,JSWB) = ZDIR_ALB(:)
    PSCA_ALB(:,JSWB) = ZSCA_ALB(:)
  END DO
  !
  DEALLOCATE(ZDIR_ALB)
  DEALLOCATE(ZSCA_ALB)
!-------------------------------------------------------------------------------
!
!*      10.     Chemistry /dust
!               ---------------
!
  CALL INIT_CHEMICAL_n(ILUOUT, KSV, HSV, CHT%NBEQ, CHT%CSV, CHT%NAEREQ,            &
                     CHT%NSV_CHSBEG, CHT%NSV_CHSEND, CHT%NSV_AERBEG, CHT%NSV_AEREND, &
                     CHT%CCH_NAMES, CHT%CAER_NAMES, CHT%NDSTEQ, CHT%NSV_DSTBEG,      &
                     CHT%NSV_DSTEND, CHT%NSLTEQ, CHT%NSV_SLTBEG, CHT%NSV_SLTEND,     &
                     HDSTNAMES=CHT%CDSTNAMES, HSLTNAMES=CHT%CSLTNAMES        )
!
!* Initialization of dry deposition scheme (chemistry)
!
  IF (CHT%NBEQ>0 .AND. CHT%CCH_DRY_DEP=='WES89') THEN
    ALLOCATE(CHT%XDEP(ILU,CHT%NBEQ))
  ELSE
    ALLOCATE(CHT%XDEP(0,0))
  END IF
!
!-------------------------------------------------------------------------------
END DO ! end of loop on patches
!
IF (HINIT/='ALL') THEN
  CALL END_IO_SURF_n(HPROGRAM)
  IF (LHOOK) CALL DR_HOOK('INIT_TEB_N',1,ZHOOK_HANDLE)
  RETURN
END IF
!-------------------------------------------------------------------------------
!
!*       7.     Canopy air fields:
!               ------------------
!
 CALL READ_TEB_CANOPY_n(DTCO, IOB, U, TCP, TOP, &
                        HPROGRAM)
!
!-------------------------------------------------------------------------------
!
!*      11.     Diagnostics:
!               -----------
!
 CALL DIAG_TEB_INIT_n(DGT, DGUT, &
                      HPROGRAM,ILU,ISWB)
DO JPATCH=1,TOP%NTEB_PATCH
  CALL GOTO_WRAPPER_TEB_PATCH(B, DGCT, DGMT, T, TGD, TGDPE, TGR, TGRPE, JPATCH)
  CALL DIAG_MISC_TEB_INIT_n(DGCT, DGMT, DGMTO, TOP, &
                            HPROGRAM,ILU,ISWB)
END DO ! end of loop on patches
!
!-------------------------------------------------------------------------------
!
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('INIT_TEB_N',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE INIT_TEB_n
