!     #############################################################
      SUBROUTINE INIT_TEB_n     (HPROGRAM,HINIT,                            &
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
USE MODD_IO_SURF_ASC,ONLY: CMASK
USE MODD_SNOW_PAR, ONLY : XEMISSN
!
USE MODD_READ_NAMELIST, ONLY : LNAM_READ

USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS
USE MODD_TEB_n, ONLY : T => TEB

USE MODD_TEB_PANEL_n, ONLY : TPN => TEB_PANEL

USE MODD_BEM_OPTION_n, ONLY : BOP => BEM_OPTIONS
USE MODD_BEM_n, ONLY : B => BEM

USE MODD_TEB_VEG_n, ONLY : TVG => TEB_VEG_OPTIONS

USE MODD_CH_TEB_n, ONLY : CHT => CH_TEB


USE MODD_CHS_AEROSOL,     ONLY: LVARSIGI, LVARSIGJ
USE MODD_DST_SURF,        ONLY: LVARSIG_DST, NDSTMDE, NDST_MDEBEG, LRGFIX_DST 
USE MODD_SLT_SURF,        ONLY: LVARSIG_SLT, NSLTMDE, NSLT_MDEBEG, LRGFIX_SLT
USE MODD_DIAG_TEB_n, ONLY : DGT => DIAG_TEB
USE MODD_DIAG_MISC_TEB_OPTION_n, ONLY : DGMTO => DIAG_MISC_TEB_OPTIONS
USE MODD_DIAG_UTCI_TEB_n, ONLY : DGUT => DIAG_UTCI_TEB
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
USE MODI_GOTO_TEB
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
 CALL DEFAULT_TEB(TOP%CZ0H,TOP%XTSTEP,TOP%XOUT_TSTEP, TOP%CCH_BEM, T%XDT_RES, T%XDT_OFF)
 CALL DEFAULT_CH_DEP(CHT%CCH_DRY_DEP)
 CALL DEFAULT_DIAG_TEB(DGT%N2M,DGT%LSURF_BUDGET,DGT%L2M_MIN_ZS,DGT%LRAD_BUDGET,DGT%LCOEF,DGT%LSURF_VARS, &
                       DGMTO%LSURF_MISC_BUDGET,DGMTO%LSURF_DIAG_ALBEDO,DGUT%LUTCI,DGT%LPGD,DGT%LPGD_FIX,DGT%XDIAG_TSTEP)  
!
ENDIF
!
!        0.2. Defaults from file header
!    
 CALL READ_DEFAULT_TEB_n(HPROGRAM)
!
!*       1.     Reading of configuration:
!               -------------------------
!
 CALL READ_TEB_CONF_n(HPROGRAM)
!
!* initialization of snow scheme
!
IF (HINIT=='PRE') THEN
  DO JPATCH=1,TOP%NTEB_PATCH
    CALL GOTO_TEB(JPATCH)
    CALL READ_PREP_TEB_SNOW(HPROGRAM,T%TSNOW_ROOF%SCHEME,T%TSNOW_ROOF%NLAYER, &
                                     T%TSNOW_ROAD%SCHEME,T%TSNOW_ROAD%NLAYER)
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
    CALL READ_TEB_DATE(HPROGRAM,HINIT,ILUOUT,HATMFILE,HATMFILETYPE,KYEAR,KMONTH,KDAY,PTIME,TOP%TTIME)

  CASE DEFAULT
    CALL INIT_IO_SURF_n(HPROGRAM,'TOWN  ','TEB   ','READ ')
    CALL READ_SURF(HPROGRAM,'DTCUR',TOP%TTIME,IRESP)
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
 CALL INIT_IO_SURF_n(HPROGRAM,'TOWN  ','TEB   ','READ ')
!
 CALL READ_SURF(HPROGRAM,'VERSION',IVERSION,IRESP)
 CALL READ_SURF(HPROGRAM,'BUG',IBUGFIX,IRESP)
!
!         Reading of the fields
!
 CALL READ_COVER_GARDEN(HPROGRAM,TOP%LGARDEN)
!
 CALL READ_PGD_TEB_n(HPROGRAM)
!
 CALL END_IO_SURF_n(HPROGRAM)
! 
!*        Fraction of each patch in the grid mesh
!
ILU = SIZE(TOP%XCOVER,1)
!
ALLOCATE(TOP%XTEB_PATCH(ILU,TOP%NTEB_PATCH))
 CALL CONVERT_TEB(TOP%XCOVER,TOP%XTEB_PATCH)
!
 CALL SET_SURFEX_FILEIN(HPROGRAM,'PREP') ! restore input file name
 CALL INIT_IO_SURF_n(HPROGRAM,'TOWN  ','TEB   ','READ ')
!
 CALL READ_SURF(HPROGRAM,'VERSION',IVERSION,IRESP)
 CALL READ_SURF(HPROGRAM,'BUG',IBUGFIX,IRESP)
!
!* reads what is the option defined for road orientations & walls
!
IF (HINIT=='ALL') THEN
  TOP%CROAD_DIR='UNIF'
  TOP%CWALL_OPT='UNIF'
  IF (IVERSION>7 .OR. (IVERSION==7 .AND. IBUGFIX>=3)) THEN
    CALL READ_SURF(HPROGRAM,'ROAD_DIR',TOP%CROAD_DIR,IRESP)
    CALL READ_SURF(HPROGRAM,'WALL_OPT',TOP%CWALL_OPT,IRESP)
  END IF
END IF
 CALL END_IO_SURF_n(HPROGRAM)
!-----------------------------------------------------------------------------------
!
!*              LOOP ON TEB PATCHES
!               -------------------
!
DO JPATCH=1,TOP%NTEB_PATCH

  CALL GOTO_TEB(JPATCH)
  !-----------------------------------------------------------------------------------
  !
  !*       3.     Physiographic data fields from land cover:
  !               -----------------------------------------
  !
  ALLOCATE(T%XZ0_TOWN     (ILU))
  ALLOCATE(T%XALB_ROOF    (ILU))
  ALLOCATE(T%XEMIS_ROOF   (ILU))
  ALLOCATE(T%XALB_ROAD    (ILU))
  ALLOCATE(T%XEMIS_ROAD   (ILU))
  ALLOCATE(T%XALB_WALL    (ILU))
  ALLOCATE(T%XEMIS_WALL   (ILU))
  ALLOCATE(T%XBLD         (ILU))
  ALLOCATE(T%XROAD_DIR    (ILU))
  ALLOCATE(T%XROAD        (ILU))
  ALLOCATE(T%XBLD_HEIGHT  (ILU))
  ALLOCATE(T%XWALL_O_HOR  (ILU))
  ALLOCATE(T%XCAN_HW_RATIO(ILU))
  ALLOCATE(T%XROAD_O_GRND(ILU))
  ALLOCATE(T%XGARDEN_O_GRND(ILU))
  ALLOCATE(T%XWALL_O_GRND(ILU))
  ALLOCATE(T%XWALL_O_BLD(ILU))
  ALLOCATE(T%XH_TRAFFIC   (ILU))
  ALLOCATE(T%XLE_TRAFFIC  (ILU))
  ALLOCATE(T%XH_INDUSTRY  (ILU))
  ALLOCATE(T%XLE_INDUSTRY (ILU))
  ALLOCATE(T%XHC_ROOF     (ILU,TOP%NROOF_LAYER))
  ALLOCATE(T%XTC_ROOF     (ILU,TOP%NROOF_LAYER))
  ALLOCATE(T%XD_ROOF      (ILU,TOP%NROOF_LAYER))
  ALLOCATE(T%XHC_ROAD     (ILU,TOP%NROAD_LAYER))
  ALLOCATE(T%XTC_ROAD     (ILU,TOP%NROAD_LAYER))
  ALLOCATE(T%XD_ROAD      (ILU,TOP%NROAD_LAYER))
  ALLOCATE(T%XHC_WALL     (ILU,TOP%NWALL_LAYER))
  ALLOCATE(T%XTC_WALL     (ILU,TOP%NWALL_LAYER))
  ALLOCATE(T%XD_WALL      (ILU,TOP%NWALL_LAYER))
  ALLOCATE(T%XROUGH_ROOF      (ILU))
  ALLOCATE(T%XROUGH_WALL      (ILU))
  ALLOCATE(T%XRESIDENTIAL     (ILU))
  ALLOCATE(T%XGREENROOF       (ILU))
  ALLOCATE(T%XGARDEN          (ILU))
  ALLOCATE(TPN%XEMIS_PANEL      (ILU))
  ALLOCATE(TPN%XALB_PANEL       (ILU))
  ALLOCATE(TPN%XEFF_PANEL       (ILU))
  ALLOCATE(TPN%XFRAC_PANEL      (ILU))
  !
  T%XROAD_DIR(:) = 0.
  T%XROAD    (:) = 0.
  !
  ZDEF_ROAD_DIR = 0.
  IF (TOP%CROAD_DIR/='UNIF') THEN
    !* road direction if not specified by the user depends on patch number
    !  First patch has a Notrh-South road. Other patches have roads spaced by
    !  regular angles
    ZDEF_ROAD_DIR = 180. * FLOAT(JPATCH-1) / FLOAT(TOP%NTEB_PATCH)
  END IF
  !
  CALL CONVERT_PATCH_TEB(TOP%XCOVER, TOP%LCOVER, ZDEF_ROAD_DIR,                          &
                      PZ0_TOWN=T%XZ0_TOWN,                                         &
                      PALB_ROOF=T%XALB_ROOF,                                       &
                      PEMIS_ROOF=T%XEMIS_ROOF,PHC_ROOF=T%XHC_ROOF,PTC_ROOF=T%XTC_ROOF, &
                      PD_ROOF=T%XD_ROOF,                                           &
                      PALB_ROAD=T%XALB_ROAD,                                       &
                      PEMIS_ROAD=T%XEMIS_ROAD,PHC_ROAD=T%XHC_ROAD,PTC_ROAD=T%XTC_ROAD, &
                      PD_ROAD=T%XD_ROAD,                                           &
                      PALB_WALL=T%XALB_WALL,                                       &
                      PEMIS_WALL=T%XEMIS_WALL,PHC_WALL=T%XHC_WALL,PTC_WALL=T%XTC_WALL, &
                      PD_WALL=T%XD_WALL,                                           &
                      PBLD_HEIGHT=T%XBLD_HEIGHT,                                   &
                      PWALL_O_HOR=T%XWALL_O_HOR,PBLD=T%XBLD, PROAD_DIR=T%XROAD_DIR,    &
                      PGARDEN=T%XGARDEN,                                           &
                      PH_TRAFFIC=T%XH_TRAFFIC, PLE_TRAFFIC=T%XLE_TRAFFIC,            &
                      PH_INDUSTRY=T%XH_INDUSTRY, PLE_INDUSTRY=T%XLE_INDUSTRY,        &
                      PROUGH_ROOF = T%XROUGH_ROOF, PROUGH_WALL = T%XROUGH_WALL,      &
                      PRESIDENTIAL = T%XRESIDENTIAL,                               &
                      PGREENROOF = T%XGREENROOF,                                   &
                      PEMIS_PANEL=TPN%XEMIS_PANEL, PALB_PANEL=TPN%XALB_PANEL,            &
                      PEFF_PANEL=TPN%XEFF_PANEL, PFRAC_PANEL=TPN%XFRAC_PANEL             )
  !
  IF (.NOT. TOP%LGREENROOF .AND. MAXVAL(T%XGREENROOF)>0. ) THEN !<== A paralleliser pour un stop propre
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
  ALLOCATE(T%XSVF_ROAD  (ILU))
  ALLOCATE(T%XSVF_GARDEN(ILU))
  ALLOCATE(T%XSVF_WALL  (ILU))
  !
  ALLOCATE(B%XGR          (ILU))
  ALLOCATE(B%XALB_WIN     (ILU))
  ALLOCATE(B%XF_WASTE_CAN (ILU))
  !
  !
  CALL TEB_MORPHO(HPROGRAM, T%XBLD, T%XWALL_O_HOR, T%XGARDEN, T%XBLD_HEIGHT, T%XROAD, T%XROAD_O_GRND, &
                T%XGARDEN_O_GRND, T%XWALL_O_GRND, T%XCAN_HW_RATIO, T%XSVF_ROAD, T%XSVF_GARDEN,    &
                T%XSVF_WALL, T%XZ0_TOWN, T%XWALL_O_BLD, T%XH_TRAFFIC, T%XLE_TRAFFIC               )
                !
  !-------------------------------------------------------------------------------
  !
  !*       6.     Building Energy Model
  !               ---------------------
  !
  CALL INIT_BEM_n(ILUOUT)
  !
  !-------------------------------------------------------------------------------
  !
  !*      7.      Case of urban green areas
  !               -------------------------
  !
  IF (TOP%LGARDEN) THEN
  !
    CALL SET_SURFEX_FILEIN(HPROGRAM,'PGD ') ! change input file name to pgd name
    CALL INIT_IO_SURF_n(HPROGRAM,'TOWN  ','TEB   ','READ ')     
    IF (JPATCH==1) CALL INIT_TEB_VEG_OPTIONS_n(HPROGRAM)
    CALL INIT_TEB_GARDEN_PGD_n(HPROGRAM,HINIT,(JPATCH==1),KI,KSV,HSV,IVERSION,IBUGFIX,PCO2,PRHOA)
    ! Case of urban green roofs
    IF (TOP%LGREENROOF) CALL INIT_TEB_GREENROOF_PGD_n(HPROGRAM,HINIT,(JPATCH==1), &
                                                  KI,KSV,HSV,IVERSION,PCO2,PRHOA)
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
 CALL INIT_IO_SURF_n(HPROGRAM,'TOWN  ','TEB   ','READ ')     
 CALL READ_PGD_TEB_IRRIG_n(HPROGRAM)
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
 CALL INIT_IO_SURF_n(HPROGRAM,'TOWN  ','TEB   ','READ ')
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
  CALL GOTO_TEB(JPATCH)
!
!* TEB fields
  CALL READ_TEB_n(HPROGRAM,JPATCH)
!
  ALLOCATE(T%XAC_ROOF    (ILU))
  ALLOCATE(T%XAC_ROAD    (ILU))
  ALLOCATE(T%XAC_WALL    (ILU))
  ALLOCATE(T%XAC_TOP     (ILU))
  ALLOCATE(T%XAC_ROOF_WAT(ILU))
  ALLOCATE(T%XAC_ROAD_WAT(ILU))
  ALLOCATE(T%XQSAT_ROOF  (ILU))
  ALLOCATE(T%XQSAT_ROAD  (ILU))
  ALLOCATE(T%XDELT_ROOF  (ILU))
  ALLOCATE(T%XDELT_ROAD  (ILU))
!
!* Case of urban green areas
  IF (TOP%LGARDEN) THEN
!    CALL SET_SURFEX_FILEIN(HPROGRAM,'PREP') ! change input file name to pgd name
!    CALL INIT_IO_SURF_n(HPROGRAM,'TOWN  ','TEB   ','READ ')       
    CALL INIT_TEB_GARDEN_n(HPROGRAM,HINIT,KI,KSW,PSW_BANDS)
  ! Case of urban green roofs
    IF (TOP%LGREENROOF) CALL INIT_TEB_GREENROOF_n(HPROGRAM,HINIT,KI,KSV,PSW_BANDS)
!    CALL END_IO_SURF_n(HPROGRAM)
  ENDIF
!-------------------------------------------------------------------------------
!
!*      10.     Infra-red Radiative fields:
!               --------------------------
!
!* snow long-wave properties (not initialized in read_gr_snow)
!
  CALL INIT_SNOW_LW(XEMISSN,T%TSNOW_ROOF)
  CALL INIT_SNOW_LW(XEMISSN,T%TSNOW_ROAD)
!
  IF (TOP%LGARDEN) THEN
    ZDIR_SW=0. ! night as first guess for albedo computation
    ZSCA_SW=0. !
    CALL GARDEN_PROPERTIES(ZDIR_SW, ZSCA_SW, PSW_BANDS, KSW,     &
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
    CALL GREENROOF_PROPERTIES(ZDIR_SW, ZSCA_SW, PSW_BANDS, KSW,              &
                              ZTS_GREENROOF, ZEMIS_GREENROOF, ZALB_GREENROOF )  
  ELSE
    ZALB_GREENROOF  = XUNDEF
    ZEMIS_GREENROOF = XUNDEF
    ZTS_GREENROOF   = XUNDEF
  END IF
!
!* averaged albedo, emissivity and radiative temperature
!
  CALL AVERAGED_TSRAD_TEB(T%XEMIS_ROOF,T%XT_ROOF(:,1),       &
                        T%XEMIS_ROAD,T%XT_ROAD(:,1),       &
                        T%XEMIS_WALL,                    &
                        T%XT_WALL_A(:,1),                &
                        T%XT_WALL_B(:,1),                &
                        ZEMIS_GARDEN, ZTS_GARDEN,      &
                        ZEMIS_GREENROOF, ZTS_GREENROOF,&
                        T%TSNOW_ROOF,T%TSNOW_ROAD,         &
                        T%XROAD, T%XGREENROOF, T%XGARDEN,    &
                        T%XBLD,T%XWALL_O_HOR,              &
                        T%XSVF_ROAD,T%XSVF_WALL,           &
                        T%XSVF_GARDEN,                   &
                        PEMIS,PTSRAD, B%XT_WIN1,         &
                        B%XGR                            )
!
!
!*       9.     Visible and near-infra-red Radiative fields:
!               -------------------------------------------
!
  ALLOCATE(ZDIR_ALB(ILU))
  ALLOCATE(ZSCA_ALB(ILU))
!
  CALL AVERAGED_ALBEDO_TEB(TOP%CBEM,TOP%CROAD_DIR,TOP%CWALL_OPT,PZENITH,PAZIM, &
                       T%XBLD, T%XGARDEN, T%XROAD_DIR, T%XROAD, T%XGREENROOF,&
                       TPN%XFRAC_PANEL, TPN%XALB_PANEL,                    &
                       T%XWALL_O_HOR, T%XCAN_HW_RATIO,                 &
                       T%XALB_ROOF,                                  &
                       T%XALB_ROAD, T%XSVF_ROAD,                       &
                       T%XALB_WALL, T%XSVF_WALL,                       &
                       ZALB_GARDEN, T%XSVF_GARDEN,                   &
                       ZALB_GREENROOF,                             &
                       T%TSNOW_ROOF, T%TSNOW_ROAD,                     &
                       B%XGR, B%XSHGC, B%XSHGC_SH, B%XABS_WIN, B%XALB_WIN,   &
                       B%LSHAD_DAY,                                  &
                       ZDIR_ALB, ZSCA_ALB, B%XTRAN_WIN               )  

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
 CALL READ_TEB_CANOPY_n(HPROGRAM)
!
!-------------------------------------------------------------------------------
!
!*      11.     Diagnostics:
!               -----------
!
 CALL DIAG_TEB_INIT_n(HPROGRAM,ILU,ISWB)
DO JPATCH=1,TOP%NTEB_PATCH
  CALL GOTO_TEB(JPATCH)
  CALL DIAG_MISC_TEB_INIT_n(HPROGRAM,ILU,ISWB)
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
