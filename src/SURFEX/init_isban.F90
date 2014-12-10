!#############################################################
SUBROUTINE INIT_ISBA_n    (HPROGRAM,HINIT,OLAND_USE,                    &
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
!!	V. Masson   *Meteo France*	
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
USE MODD_DATA_COVER_n, ONLY : XDATA_VEGTYPE, XDATA_NATURE
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
USE MODD_ISBA_n,   ONLY : CROUGH ,CISBA, CPHOTO, CRUNOFF, CALBEDO, CSCOND,    &
                          CC1DRY, CSOILFRZ, CDIFSFCOND, CSNOWRES, CRESPSL,    &
                          LMEB_PATCH,                                         &
                          NNLITTER, NNLITTLEVS, NNSOILCARB, NPATCH,           &
                          TSNOW, TTIME, XTSTEP, XOUT_TSTEP,                   &
                          LGLACIER, LVEGUPD, LCANOPY_DRAG,                    &
                          CCPSURF, CHORT, XCGMAX, XCDRAG, CKSAT,              &
                          LSOC, CRAIN, LSPINUPCARBS, LSPINUPCARBW,            &
                          NNBYEARSOLD, NSPINS, NSPINW, LAGRI_TO_GRASS,        &
                          XSPINMAXS, XSPINMAXW, NNBYEARSPINS, NNBYEARSPINW,   &
                          LSNOWDRIFT,LSNOWDRIFT_SUBLIM,LSNOW_ABS_ZENITH,      &
                          CSNOWMETAMO,CSNOWRAD,XCO2_START,XCO2_END,LNITRO_DILU
!
USE MODD_CH_ISBA_n,      ONLY : LCH_BIO_FLUX, CCH_DRY_DEP  

USE MODD_DIAG_ISBA_n,    ONLY : N2M, LSURF_BUDGET, LRAD_BUDGET,          &
                                  XDIAG_TSTEP, LPGD,  L2M_MIN_ZS, LCOEF, &
                                  LSURF_VARS, LPATCH_BUDGET  
USE MODD_DIAG_EVAP_ISBA_n, ONLY : LSURF_EVAP_BUDGET, LSURF_BUDGETC, LRESET_BUDGETC, &
                                  LWATER_BUDGET
USE MODD_DIAG_MISC_ISBA_n, ONLY : LSURF_MISC_BUDGET, LSURF_DIAG_ALBEDO, LSURF_MISC_DIF  
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
LSURF_DIAG_ALBEDO = .FALSE.
!
IF (LNAM_READ) THEN
 !
 !*       0.     Defaults
!               --------

 !        0.1. Hard defaults
 !      
 CALL DEFAULT_ISBA(XTSTEP, XOUT_TSTEP,                           &
                     CROUGH,CRUNOFF,CALBEDO,CSCOND,              &
                     CC1DRY, CSOILFRZ, CDIFSFCOND, CSNOWRES,     &
                     CCPSURF, XCGMAX, XCDRAG, CKSAT, LSOC,       &
                     CRAIN, CHORT, LGLACIER, LCANOPY_DRAG,       &
                     LVEGUPD, LSPINUPCARBS, LSPINUPCARBW,        &
                     XSPINMAXS, XSPINMAXW, XCO2_START, XCO2_END, &
                     NNBYEARSPINS, NNBYEARSPINW, LNITRO_DILU     )
 !                  
 CALL DEFAULT_CH_DEP(CCH_DRY_DEP)
 CALL DEFAULT_CH_BIO_FLUX(LCH_BIO_FLUX)                  
 CALL DEFAULT_DIAG_ISBA(N2M,LSURF_BUDGET,L2M_MIN_ZS,LRAD_BUDGET,   &
                        LCOEF,LSURF_VARS,LSURF_EVAP_BUDGET,        &
                        LSURF_MISC_BUDGET,LSURF_BUDGETC,           &
                        LSURF_MISC_DIF,LPATCH_BUDGET,              &
                        LPGD,LRESET_BUDGETC,LWATER_BUDGET,         &
                        XDIAG_TSTEP                                )  
 !
 CALL DEFAULT_CROCUS(LSNOWDRIFT,LSNOWDRIFT_SUBLIM,LSNOW_ABS_ZENITH,&
                     CSNOWMETAMO,CSNOWRAD)
 ! 
ENDIF
!
!        0.2. Defaults from file header
!    
 CALL READ_DEFAULT_ISBA_n(HPROGRAM)
!
 CALL READ_ISBA_CONF_n(HPROGRAM)
!
 CALL INIT_IO_SURF_n(HPROGRAM,'FULL  ','ISBA  ','READ ')
 CALL READ_SURF(HPROGRAM,'VERSION',IVERSION,IRESP)
 CALL END_IO_SURF_n(HPROGRAM)
!
!*       1.     Reading of configuration:
!               -------------------------
!
!* initialization of snow and carbon schemes
!
NNBYEARSOLD = 1
NSPINS      = 1
NSPINW      = 1
!
IF (HINIT=='PRE') THEN 
  CALL READ_PREP_ISBA_SNOW(HPROGRAM,TSNOW%SCHEME,TSNOW%NLAYER) 
!
!* initialization of soil carbon scheme
!
  CALL READ_PREP_ISBA_CARBON(HPROGRAM,CRESPSL)
!
  IF (CRESPSL=='CNT') THEN
    NNLITTER = 2
    NNLITTLEVS = 2
    NNSOILCARB = 3
  ELSE
    NNLITTER = 0
    NNLITTLEVS = 0
    NNSOILCARB = 0
  ENDIF

ELSEIF (HINIT=='ALL') THEN
!
  CALL INIT_IO_SURF_n(HPROGRAM,'NATURE','ISBA  ','READ ')
!
  IF (IVERSION<6) THEN
    CRESPSL='DEF'
  ELSE  
    CALL READ_SURF(HPROGRAM,'RESPSL',CRESPSL,IRESP)
    CALL READ_SURF(HPROGRAM,'NLITTER',NNLITTER,IRESP)
    CALL READ_SURF(HPROGRAM,'NLITTLEVS',NNLITTLEVS,IRESP)
    CALL READ_SURF(HPROGRAM,'NSOILCARB',NNSOILCARB,IRESP)
    IF(IVERSION>=7.AND.(LSPINUPCARBS.OR.LSPINUPCARBW))THEN
      CALL READ_SURF(HPROGRAM,'NBYEARSOLD',NNBYEARSOLD,IRESP)
    ELSE
      NNBYEARSOLD=NUNDEF
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
    TTIME%TDATE%YEAR = NUNDEF
    TTIME%TDATE%MONTH= NUNDEF
    TTIME%TDATE%DAY  = NUNDEF
    TTIME%TIME       = XUNDEF

  CASE ('PRE')
    CALL PREP_CTRL_ISBA(N2M,LSURF_BUDGET,L2M_MIN_ZS,LRAD_BUDGET,LCOEF,LSURF_VARS,&
                          LSURF_EVAP_BUDGET,LSURF_MISC_BUDGET,LSURF_BUDGETC,     &
                          LPATCH_BUDGET,LSURF_MISC_DIF,ILUOUT                    )    
    IF (LNAM_READ) CALL READ_NAM_PREP_ISBA_n(HPROGRAM)                        
    CALL READ_ISBA_DATE(HPROGRAM,HINIT,ILUOUT,HATMFILE,HATMFILETYPE,KYEAR,KMONTH,KDAY,PTIME,TTIME)

  CASE DEFAULT
    CALL INIT_IO_SURF_n(HPROGRAM,'NATURE','ISBA  ','READ ')
    CALL READ_SURF(HPROGRAM,'DTCUR',TTIME,IRESP)
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
 CALL INIT_IO_SURF_n(HPROGRAM,'NATURE','ISBA  ','READ ')
!
!
!*       2.1    Cover, soil and orographic fields:
!               ---------------------------------
!
 CALL READ_PGD_ISBA_n(HPROGRAM,OLAND_USE)
!
ISIZE_LMEB_PATCH=COUNT(LMEB_PATCH(:))
!
!
!*       2.2    Check:
!               ------
!
IF ( CPHOTO/='NON' .AND. NPATCH/=12) THEN
  CALL ABOR1_SFX('INIT_ISBAN: INCONSISTENCY BETWEEN CPHOTO AND NPATCH')
ENDIF
!
IF (HINIT=='PRE' .AND. TSNOW%SCHEME.NE.'3-L' .AND. TSNOW%SCHEME.NE.'CRO' .AND. CISBA=='DIF') THEN
    CALL ABOR1_SFX("INIT_ISBAN: WITH CISBA = DIF, CSNOW MUST BE 3-L OR CRO")
ENDIF
IF ( CPHOTO/='LAI' .AND. CPHOTO/='LST' .AND. CPHOTO/='NIT' .AND. CPHOTO/='NCB' .AND. LAGRIP) THEN
  CALL ABOR1_SFX('INIT_ISBAN: INCONSISTENCY BETWEEN CPHOTO AND LAGRIP')
ENDIF
IF ( CPHOTO/='NCB' .AND. CRESPSL=='CNT') THEN
  CALL ABOR1_SFX('INIT_ISBAN: INCONSISTENCY BETWEEN CPHOTO AND CRESPSL')
ENDIF
IF (HINIT=='PRE' .AND. ISIZE_LMEB_PATCH>0 .AND. TSNOW%SCHEME.NE.'3-L' .AND. TSNOW%SCHEME.NE.'CRO') THEN
    CALL ABOR1_SFX("INIT_ISBAN: WITH LMEB_PATCH = TRUE, CSNOW MUST BE 3-L OR CRO")
ENDIF
IF(CPHOTO/='NCB'.AND.LSPINUPCARBW)THEN
  CALL ABOR1_SFX('INIT_ISBAN: INCONSISTENCY BETWEEN CPHOTO AND LSPINUPCARBW (if not NCB must be false)')
ENDIF
IF(CRESPSL/='CNT'.AND.LSPINUPCARBS)THEN
  CALL ABOR1_SFX('INIT_ISBAN: INCONSISTENCY BETWEEN CRESPSL AND LSPINUPCARBS (if not CNT must be false)')
ENDIF
IF(LSPINUPCARBW.AND.REAL(NNBYEARSPINW)>REAL(NNBYEARSPINS)*0.5)THEN
  WRITE(ILUOUT,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  WRITE(ILUOUT,*)'INIT_ISBAN: INCONSISTENCY BETWEEN NNBYEARSPINW AND NNBYEARSPINS'
  WRITE(ILUOUT,*)'NNBYEARSPINW MUST BE < TO 0.5 * NNBYEARSPINS'
  WRITE(ILUOUT,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  CALL ABOR1_SFX('INIT_ISBAN: INCONSISTENCY BETWEEN NNBYEARSPINW AND NNBYEARSPINS')
ENDIF
IF(LSPINUPCARBS.AND.(XCO2_START==XUNDEF.OR.XCO2_END==XUNDEF))THEN
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
ISPINEND=NNBYEARSPINS-NINT(NNBYEARSPINS*XSPIN_CO2)
!
LAGRI_TO_GRASS = .FALSE.
!
IF ( LSPINUPCARBS .AND. (NNBYEARSOLD <= ISPINEND) ) THEN
!
   LAGRI_TO_GRASS = .TRUE.
!
   CALL INI_DATA_PARAM(XDATA_VEGTYPE, PSURF=XDATA_NATURE, PH_TREE=XDATA_H_TREE,PLAI=XDATA_LAI,                   &
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
                                  PSEFOLD_ST=XDATA_SEFOLD_ST, PDMAX_ST=XDATA_DMAX_ST, OAGRI_TO_GRASS=LAGRI_TO_GRASS)
!
   ZCO2(:) = PRHOA(:) * XCO2_START * 1.E-6 * XMCO2 / XMD
!
ELSEIF(LSPINUPCARBS .AND. (NNBYEARSOLD > ISPINEND) .AND. (NNBYEARSOLD <= NNBYEARSPINS) )THEN
!
   ZSPINCO2 = XCO2_START + (XCO2_END-XCO2_START) * REAL(NNBYEARSOLD - ISPINEND) / REAL(NNBYEARSPINS - ISPINEND)
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
CALL COMPUTE_ISBA_PARAMETERS(HPROGRAM,HINIT,OLAND_USE,                  &
                             KI,KSV,KSW,                                &
                             HSV,ZCO2,PRHOA,                            &
                             PZENITH,PSW_BANDS,PDIR_ALB,PSCA_ALB,       &
                             PEMIS,PTSRAD,PTSURF,                       &
                             HTEST                                )
!
IF ( CSNOWMETAMO/="B92" ) THEN
  CALL READ_FZ06('drdt_bst_fit_60.nc')
ENDIF
!
IF ( CSNOWRAD=="TAR" .OR. CSNOWRAD=="TA1" .OR.  CSNOWRAD=="TA2" ) THEN
  CALL INIT_TARTES()
END IF
!
IF (LHOOK) CALL DR_HOOK('INIT_ISBA_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE INIT_ISBA_n
