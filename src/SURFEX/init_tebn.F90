!     #############################################################
      SUBROUTINE INIT_TEB_n     (HPROGRAM,HINIT,                            &
                                   KI,KSV,KSW,                                &
                                   HSV,PCO2,PRHOA,                            &
                                   PZENITH,PAZIM,PSW_BANDS,PDIR_ALB,PSCA_ALB, &
                                   PEMIS,PTSRAD,                              &
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
!!	V. Masson   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2003
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SNOW_PAR, ONLY : XEMISSN
!
USE MODD_READ_NAMELIST, ONLY : LNAM_READ
!
USE MODD_TEB_n,           ONLY: XTSTEP, XOUT_TSTEP, TTIME, XCOVER,                       &
                                  XH_TRAFFIC, XLE_TRAFFIC, XH_INDUSTRY, XLE_INDUSTRY,      &
                                  XZ0_TOWN, XBLD, XGARDEN, XVEG_ROOF,                      &
                                  XROAD, XBLD_HEIGHT, XWALL_O_HOR, XCAN_HW_RATIO,          &
                                  XALB_ROOF, XEMIS_ROOF, XHC_ROOF,XTC_ROOF, XD_ROOF,       &
                                  XALB_ROAD, XEMIS_ROAD, XHC_ROAD,XTC_ROAD, XD_ROAD,       &
                                  XALB_WALL, XEMIS_WALL, XHC_WALL,XTC_WALL, XD_WALL,       &
                                  XSVF_ROAD, XSVF_GARDEN, XSVF_WALL,                       &
                                  TSNOW_ROOF, TSNOW_ROAD,                                  &
                                  NROOF_LAYER, NROAD_LAYER, NWALL_LAYER,                   &
                                  XT_ROOF, XT_ROAD, XT_WALL, CZ0H,                         &
                                  XT_CANYON, XQ_CANYON,                                    &
                                  XAC_ROOF, XAC_ROAD, XAC_WALL, XAC_TOP,                   &
                                  XAC_ROOF_WAT, XAC_ROAD_WAT,                              &
                                  XQSAT_ROOF, XQSAT_ROAD, XDELT_ROOF, XDELT_ROAD,          &
                                  LGARDEN 
USE MODD_CH_TEB_n,        ONLY: XDEP, CCH_DRY_DEP, CSV, CCH_NAMES,                       &
                                  NBEQ, NSV_CHSBEG, NSV_CHSEND,                            &
                                  NAEREQ, NSV_AERBEG, NSV_AEREND, CAER_NAMES,              &
                                  NSV_DSTBEG, NSV_DSTEND, NDSTEQ, CDSTNAMES,               &
                                  NSV_SLTBEG, NSV_SLTEND, NSLTEQ, CSLTNAMES  
USE MODD_CHS_AEROSOL,     ONLY: LVARSIGI, LVARSIGJ
USE MODD_DST_SURF,        ONLY: LVARSIG_DST, NDSTMDE, NDST_MDEBEG, LRGFIX_DST 
USE MODD_SLT_SURF,        ONLY: LVARSIG_SLT, NSLTMDE, NSLT_MDEBEG, LRGFIX_SLT
USE MODD_DIAG_TEB_n,      ONLY: N2M, LSURF_BUDGET, LRAD_BUDGET, XDIAG_TSTEP,             &
                                  LPGD, LPGD_FIX, L2M_MIN_ZS, LCOEF, LSURF_VARS  
USE MODD_DIAG_MISC_TEB_n, ONLY: LSURF_MISC_BUDGET, LSURF_BUDGETC, LRESET_BUDGETC,        &
                                  LSURF_DIAG_ALBEDO, LSURF_EVAP_BUDGET  
USE MODD_SURF_PAR,        ONLY: XUNDEF, NUNDEF
!
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
USE MODI_CONVERT_PATCH_TEB
USE MODI_INIT_SNOW_LW
USE MODI_AVERAGED_TSRAD_TEB
USE MODI_AVERAGED_ALBEDO_TEB
USE MODI_DIAG_TEB_INIT_n
USE MODI_END_IO_SURF_n
USE MODI_GET_LUOUT
USE MODI_READ_SURF
USE MODI_READ_PREP_TEB_SNOW
USE MODI_READ_TEB_DATE
USE MODI_READ_NAM_PREP_TEB_n
USE MODI_INIT_CHEMICAL_n
USE MODI_GARDEN_PROPERTIES
USE MODI_READ_TEB_CANOPY_n
USE MODI_INIT_TEB_GARDEN_n
USE MODI_ABOR1_SFX
USE MODI_READ_COVER_GARDEN
USE MODI_WRITE_COVER_TEX_TEB
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
CHARACTER(LEN=6),                   INTENT(IN)  :: HPROGRAM    ! program calling surf. schemes
CHARACTER(LEN=3),                   INTENT(IN)  :: HINIT       ! choice of fields to initialize
INTEGER,                            INTENT(IN)  :: KI          ! number of points
INTEGER,                            INTENT(IN)  :: KSV         ! number of scalars
INTEGER,                            INTENT(IN)  :: KSW         ! number of short-wave spectral bands
CHARACTER(LEN=6), DIMENSION(:),   INTENT(IN)  :: HSV         ! name of all scalar variables
REAL,             DIMENSION(:),    INTENT(IN)  :: PCO2        ! CO2 concentration (kg/m3)
REAL,             DIMENSION(:),    INTENT(IN)  :: PRHOA       ! air density
REAL,             DIMENSION(:),    INTENT(IN)  :: PZENITH     ! solar zenithal angle
REAL,             DIMENSION(:),    INTENT(IN)  :: PAZIM       ! solar azimuthal angle (rad from N, clock)
REAL,             DIMENSION(:),   INTENT(IN)  :: PSW_BANDS   ! middle wavelength of each band
REAL,             DIMENSION(:,:),INTENT(OUT) :: PDIR_ALB    ! direct albedo for each band
REAL,             DIMENSION(:,:),INTENT(OUT) :: PSCA_ALB    ! diffuse albedo for each band
REAL,             DIMENSION(:),    INTENT(OUT) :: PEMIS       ! emissivity
REAL,             DIMENSION(:),    INTENT(OUT) :: PTSRAD      ! radiative temperature
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
!              local variables for urban green areas
REAL, DIMENSION(SIZE(PZENITH),SIZE(PSW_BANDS))         :: ZDIR_ALB_GARDEN  ! direct  albedo for each band
REAL, DIMENSION(SIZE(PZENITH),SIZE(PSW_BANDS))         :: ZSCA_ALB_GARDEN  ! diffuse albedo for each band
REAL, DIMENSION(SIZE(PZENITH),SIZE(PSW_BANDS))         :: ZDIR_SW          ! direct  SW for each band
REAL, DIMENSION(SIZE(PZENITH),SIZE(PSW_BANDS))         :: ZSCA_SW          ! diffuse SW for each band
REAL, DIMENSION(SIZE(PZENITH))             :: ZEMIS_GARDEN     ! emissivity
REAL, DIMENSION(SIZE(PZENITH))             :: ZALB_GARDEN      ! albedo
REAL, DIMENSION(SIZE(PZENITH))             :: ZTS_GARDEN       ! radiative temperature
!
INTEGER                         :: ILU              ! sizes of TEB arrays
INTEGER                         :: ILUOUT           ! unit of output listing file
INTEGER                         :: IRESP            ! return code
!
INTEGER                         :: ISWB             ! number of shortwave spectral bands
INTEGER                         :: JSWB             ! loop on shortwave spectral bands
!
REAL, DIMENSION(:), ALLOCATABLE :: ZDIR_ALB         ! direct town albedo
REAL, DIMENSION(:), ALLOCATABLE :: ZSCA_ALB         ! diffuse town albedo
!
!
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
!         Other littel things
!
!
PDIR_ALB = XUNDEF
PSCA_ALB = XUNDEF
PEMIS    = XUNDEF
PTSRAD   = XUNDEF
!
LSURF_DIAG_ALBEDO = .FALSE.
LSURF_EVAP_BUDGET = .FALSE.
IF (LNAM_READ) THEN
 !
 !*       0.     Defaults
 !               --------
 !
 !        0.1. Hard defaults
 !      
 CALL DEFAULT_TEB(CZ0H,XTSTEP,XOUT_TSTEP)
 CALL DEFAULT_CH_DEP(CCH_DRY_DEP)
 CALL DEFAULT_DIAG_TEB(N2M,LSURF_BUDGET,L2M_MIN_ZS,LRAD_BUDGET,LCOEF,LSURF_VARS, &
                       LSURF_MISC_BUDGET,LSURF_BUDGETC,LRESET_BUDGETC,           &
                       LPGD,LPGD_FIX,XDIAG_TSTEP)  
!
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
!
CALL READ_TEB_CONF_n(HPROGRAM)
!
!* initialization of snow scheme
!
IF (HINIT=='PRE') THEN
  CALL READ_PREP_TEB_SNOW(HPROGRAM,TSNOW_ROOF%SCHEME,TSNOW_ROOF%NLAYER, &
                                   TSNOW_ROAD%SCHEME,TSNOW_ROAD%NLAYER)
ENDIF
! 
!
!
!*       2.     Cover fields and grid:
!               ---------------------
!* date
!
SELECT CASE (HINIT)
  CASE ('PGD')
    TTIME%TDATE%YEAR = NUNDEF
    TTIME%TDATE%MONTH= NUNDEF
    TTIME%TDATE%DAY  = NUNDEF
    TTIME%TIME       = XUNDEF

  CASE ('PRE')
    CALL PREP_CTRL_TEB(N2M,LSURF_BUDGET,L2M_MIN_ZS,LRAD_BUDGET,LCOEF,LSURF_VARS,&
                         LSURF_EVAP_BUDGET,LSURF_MISC_BUDGET,LSURF_BUDGETC,ILUOUT )           
    IF (LNAM_READ) CALL READ_NAM_PREP_TEB_n(HPROGRAM)   
    CALL READ_TEB_DATE(HPROGRAM,HINIT,ILUOUT,HATMFILE,HATMFILETYPE,KYEAR,KMONTH,KDAY,PTIME,TTIME)

  CASE DEFAULT
    CALL INIT_IO_SURF_n(HPROGRAM,'TOWN  ','TEB   ','READ ')
    CALL READ_SURF(HPROGRAM,'DTCUR',TTIME,IRESP)
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
!         Reading of the fields
!
CALL READ_COVER_GARDEN(HPROGRAM,LGARDEN)
CALL READ_PGD_TEB_n(HPROGRAM)
!
CALL END_IO_SURF_n(HPROGRAM)
!
!*      4.      Case of Urban green areas
!               -------------------------
!
ILU = SIZE(XCOVER,1)
!
IF (HINIT/='PGD') THEN
  ALLOCATE(XGARDEN      (ILU))
  XGARDEN(:) = 0.
  CALL CONVERT_PATCH_TEB(XCOVER, PGARDEN=XGARDEN  )
ENDIF
!
IF (LGARDEN) THEN     
 CALL INIT_TEB_GARDEN_n(HPROGRAM,HINIT,KI,KSV,KSW,HSV,PCO2,PRHOA,    &
                          PSW_BANDS,ZDIR_ALB_GARDEN,ZSCA_ALB_GARDEN, &
                          ZEMIS_GARDEN,ZTS_GARDEN,                   &
                          HATMFILE,HATMFILETYPE,HTEST                )  
ELSE
  ZEMIS_GARDEN  = XUNDEF
  ZTS_GARDEN    = XUNDEF
ENDIF
!
IF (HINIT=='PGD') THEN
  IF (LHOOK) CALL DR_HOOK('INIT_TEB_N',1,ZHOOK_HANDLE)
  RETURN
END IF
CALL SET_SURFEX_FILEIN(HPROGRAM,'PREP') ! restore input file name
!
!-----------------------------------------------------------------------------------------------------
! END READ PGD FILE
!-----------------------------------------------------------------------------------------------------
!
!*       3.     Physiographic data fields from land cover:
!               -----------------------------------------
!
ALLOCATE(XZ0_TOWN     (ILU))
ALLOCATE(XALB_ROOF    (ILU))
ALLOCATE(XEMIS_ROOF   (ILU))
ALLOCATE(XALB_ROAD    (ILU))
ALLOCATE(XEMIS_ROAD   (ILU))
ALLOCATE(XALB_WALL    (ILU))
ALLOCATE(XEMIS_WALL   (ILU))
ALLOCATE(XBLD         (ILU))
ALLOCATE(XVEG_ROOF    (ILU))
ALLOCATE(XROAD        (ILU))
ALLOCATE(XBLD_HEIGHT  (ILU))
ALLOCATE(XWALL_O_HOR  (ILU))
ALLOCATE(XCAN_HW_RATIO(ILU))
ALLOCATE(XH_TRAFFIC   (ILU))
ALLOCATE(XLE_TRAFFIC  (ILU))
ALLOCATE(XH_INDUSTRY  (ILU))
ALLOCATE(XLE_INDUSTRY (ILU))
ALLOCATE(XHC_ROOF     (ILU,NROOF_LAYER))
ALLOCATE(XTC_ROOF     (ILU,NROOF_LAYER))
ALLOCATE(XD_ROOF      (ILU,NROOF_LAYER))
ALLOCATE(XHC_ROAD     (ILU,NROAD_LAYER))
ALLOCATE(XTC_ROAD     (ILU,NROAD_LAYER))
ALLOCATE(XD_ROAD      (ILU,NROAD_LAYER))
ALLOCATE(XHC_WALL     (ILU,NWALL_LAYER))
ALLOCATE(XTC_WALL     (ILU,NWALL_LAYER))
ALLOCATE(XD_WALL      (ILU,NWALL_LAYER))
!
XROAD  (:) = 0.
!
CALL CONVERT_PATCH_TEB(XCOVER,       &
                      PZ0_TOWN=XZ0_TOWN,                   &
                      PALB_ROOF=XALB_ROOF,                 &
                      PEMIS_ROOF=XEMIS_ROOF,PHC_ROOF=XHC_ROOF,PTC_ROOF=XTC_ROOF, &
                      PD_ROOF=XD_ROOF,                     &
                      PALB_ROAD=XALB_ROAD,                 &
                      PEMIS_ROAD=XEMIS_ROAD,PHC_ROAD=XHC_ROAD,PTC_ROAD=XTC_ROAD, &
                      PD_ROAD=XD_ROAD,                     &
                      PALB_WALL=XALB_WALL,                 &
                      PEMIS_WALL=XEMIS_WALL,PHC_WALL=XHC_WALL,PTC_WALL=XTC_WALL, &
                      PD_WALL=XD_WALL,                    &
                      PBLD_HEIGHT=XBLD_HEIGHT,            &
                      PWALL_O_HOR=XWALL_O_HOR,PBLD=XBLD, PVEG_ROOF=XVEG_ROOF, &
                      PH_TRAFFIC=XH_TRAFFIC, PLE_TRAFFIC=XLE_TRAFFIC,          &
                      PH_INDUSTRY=XH_INDUSTRY, PLE_INDUSTRY=XLE_INDUSTRY  )
!
!-------------------------------------------------------------------------------
!
XROAD      (:) = 1.-(XGARDEN(:)+XBLD(:))
!
XCAN_HW_RATIO(:)    = 0.5 * XWALL_O_HOR(:) / (1.-XBLD(:))
!
!*       4.     User physiographic fields:
!               -------------------------
!
!
!
!*       5.     Sky-view-factors:
!               ----------------
!
ALLOCATE(XSVF_ROAD  (ILU))
ALLOCATE(XSVF_GARDEN(ILU))
ALLOCATE(XSVF_WALL  (ILU))
!
!
XSVF_ROAD  (:) = (SQRT(XCAN_HW_RATIO(:)**2+1.) - XCAN_HW_RATIO(:))
XSVF_GARDEN(:) = (SQRT(XCAN_HW_RATIO(:)**2+1.) - XCAN_HW_RATIO(:))
XSVF_WALL  (:) =  0.5*(XCAN_HW_RATIO(:)+1.-SQRT(XCAN_HW_RATIO(:)**2+1.))/XCAN_HW_RATIO(:)
!
WHERE (XBLD(:) .EQ. 0.)
 XBLD_HEIGHT(:)   = 0.
 XZ0_TOWN(:)      = 1.
 XWALL_O_HOR(:)   = 0.
 XCAN_HW_RATIO(:) = 0.
 XSVF_ROAD(:)     = 1.
 XSVF_WALL(:)     = 0.
 XSVF_GARDEN(:)   = 1.
ENDWHERE
!
!-------------------------------------------------------------------------------
!
!* if only physiographic fields are to be initialized, stop here.
!
CALL WRITE_COVER_TEX_TEB
!
IF (HINIT/='ALL') THEN
  IF (LHOOK) CALL DR_HOOK('INIT_TEB_N',1,ZHOOK_HANDLE)
  RETURN
END IF
!
!-------------------------------------------------------------------------------
!
!
!         Initialisation for IO
!
CALL INIT_IO_SURF_n(HPROGRAM,'TOWN  ','TEB   ','READ ')
!
!*       6.     Prognostic fields:
!               -----------------
!
CALL READ_TEB_n(HPROGRAM)
!
ALLOCATE(XAC_ROOF    (ILU))
ALLOCATE(XAC_ROAD    (ILU))
ALLOCATE(XAC_WALL    (ILU))
ALLOCATE(XAC_TOP     (ILU))
ALLOCATE(XAC_ROOF_WAT(ILU))
ALLOCATE(XAC_ROAD_WAT(ILU))
ALLOCATE(XQSAT_ROOF  (ILU))
ALLOCATE(XQSAT_ROAD  (ILU))
ALLOCATE(XDELT_ROOF  (ILU))
ALLOCATE(XDELT_ROAD  (ILU))
!
!-------------------------------------------------------------------------------
!
!*       7.     Canopy air fields:
!               ------------------
!
CALL READ_TEB_CANOPY_n(HPROGRAM)
!
!-------------------------------------------------------------------------------
!
!*       8.     Infra-red Radiative fields:
!               --------------------------
!
!* snow long-wave properties (not initialized in read_gr_snow)
!
CALL INIT_SNOW_LW(XEMISSN,TSNOW_ROOF)
CALL INIT_SNOW_LW(XEMISSN,TSNOW_ROAD)
!
IF (LGARDEN) THEN
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
!
!* averaged albedo, emissivity and radiative temperature
!
CALL AVERAGED_TSRAD_TEB(XEMIS_ROOF,XT_ROOF(:,1), &
                          XEMIS_ROAD,XT_ROAD(:,1), &
                          XEMIS_WALL,XT_WALL(:,1), &
                          ZEMIS_GARDEN, ZTS_GARDEN,&
                          TSNOW_ROOF,TSNOW_ROAD,   &
                          XROAD, XGARDEN,          &
                          XBLD,XWALL_O_HOR,        &
                          XSVF_ROAD,XSVF_WALL,     &
                          XSVF_GARDEN,             &
                          PEMIS,PTSRAD             )  
!
!
!*       9.     Visible and near-infra-red Radiative fields:
!               -------------------------------------------
!
ALLOCATE(ZDIR_ALB(ILU))
ALLOCATE(ZSCA_ALB(ILU))
!
CALL AVERAGED_ALBEDO_TEB(PZENITH,                                  &
                       XBLD, XGARDEN, XROAD,                         &
                       XWALL_O_HOR, XCAN_HW_RATIO,                   &
                       XALB_ROOF,                                    &
                       XALB_ROAD, XSVF_ROAD,                         &
                       XALB_WALL, XSVF_WALL,                         &
                       ZALB_GARDEN, XSVF_GARDEN,                     &
                       TSNOW_ROOF, TSNOW_ROAD,                       &
                       ZDIR_ALB, ZSCA_ALB                            )  
!
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
CALL INIT_CHEMICAL_n(ILUOUT, KSV, HSV, NBEQ, CSV, NAEREQ,            &
                     NSV_CHSBEG, NSV_CHSEND, NSV_AERBEG, NSV_AEREND, &
                     CCH_NAMES, CAER_NAMES, NDSTEQ, NSV_DSTBEG,      &
                     NSV_DSTEND, NSLTEQ, NSV_SLTBEG, NSV_SLTEND,     &
                     HDSTNAMES=CDSTNAMES, HSLTNAMES=CSLTNAMES        )
!
!* Initialization of dry deposition scheme (chemistry)
!
IF (NBEQ>0 .AND. CCH_DRY_DEP=='WES89') THEN
  ALLOCATE(XDEP(ILU,NBEQ))
ELSE
  ALLOCATE(XDEP(0,0))
END IF
!
!-------------------------------------------------------------------------------
!
!*      11.     Diagnostics:
!               -----------
!
CALL DIAG_TEB_INIT_n(HPROGRAM,ILU,ISWB)
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
