!#############################################################
SUBROUTINE COMPUTE_ISBA_PARAMETERS(HPROGRAM,HINIT,OLAND_USE,            &
                             KI,KSV,KSW,                                &
                             HSV,PCO2,PRHOA,                            &
                             PZENITH,PAZIM,PSW_BANDS,PDIR_ALB,PSCA_ALB, &
                             PEMIS,PTSRAD,                              &
                             KYEAR, KMONTH,KDAY, PTIME,                 &
                             HATMFILE,HATMFILETYPE,                     &
                             HTEST                                      )  
!#############################################################
!
!!****  *COMPUTE_ISBA_PARAMETERS_n* - routine to initialize ISBA
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
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_TYPE_DATE_SURF
USE MODD_TYPE_SNOW
!
USE MODD_SNOW_PAR, ONLY : XEMISSN
USE MODD_ISBA_PAR, ONLY : XTAU_ICE
USE MODD_ISBA_n,   ONLY : CROUGH,CISBA,CDIF,CPEDOTF,CPHOTO,CRUNOFF,CALBEDO,     &
                            CSCOND, CC1DRY, CSOILFRZ, CDIFSFCOND, CSNOWRES,     &
                            CRESPSL, NNBIOMASS, NNLITTER, NNLITTLEVS,           &
                            NNSOILCARB, XCLAY, XSAND, XORGMAT, XDENSITY,        &
                            XWWILT, XWFC, XWSAT,                                &
                            XCOVER, XVEG, XLAI, XRSMIN, XGAMMA, XRGL, XCV,      &
                            XDG, XZ0, XZ0_O_Z0H,                                &
                            XALBNIR_VEG, XALBVIS_VEG, XALBUV_VEG,               &
                            XEMIS, XVEGTYPE, XGMES, XRE25, XBSLAI, XLAIMIN, XGC,&
                            XDMAX, LSTRESS, XF2I,                               &
                            XSEFOLD, XH_TREE, XPATCH, NPATCH, XWRMAX_CF,        &
                            NR_NATURE_P, NSIZE_NATURE_P,                        &
                            XALBNIR_DRY, XALBVIS_DRY, XALBUV_DRY,               &
                            XALBNIR_WET, XALBVIS_WET, XALBUV_WET,               &
                            XALBNIR_SOIL, XALBVIS_SOIL, XALBUV_SOIL,            &
                            XWG, XTG, TSNOW, XALBNIR, XALBVIS, XALBUV,          &
                            XEMIS_NAT,                                          &
                            XAOSIP,XAOSIM,XAOSJP,XAOSJM,                        &
                            XHO2IP,XHO2IM,XHO2JP,XHO2JM,                        &
                            XZ0EFFIP,XZ0EFFIM,XZ0EFFJP,XZ0EFFJM, XZ0REL,        &
                            XVEGTYPE_PATCH,XROOTFRAC,XRUNOFFD,                  &
                            XCGSAT, XC1SAT, XC2REF, XC3, XC4B, XACOEF, XPCOEF,  &
                            XTAUICE, XACOEF, XPCOEF, XTAUICE, XBCOEF, XCONDSAT, &
                            XHCAPSOIL, XCONDDRY, XCONDSLD, XC4REF, XMPOTSAT,    &
                            XTDEEP, XGAMMAT, NGROUND_LAYER, TTIME,              &
                            XCE_NITRO, XCF_NITRO,                               &
                            XCNA_NITRO, XBSLAI_NITRO, CCPSURF, TSEED,           &
                            TREAP, XWATSUP, XIRRIG, XCGMAX,                     &
                            CKSAT, CTOPREG, CRAIN, CHORT,                       &
                            XTI_MIN, XTI_MAX, XTI_MEAN, XTI_STD, XTI_SKEW,      &
                            XTAB_FSAT, XTAB_WTOP, XD_ICE, XKSAT_ICE,            &
                            XFSAT, XMUF, LTRIP, LFLOOD, XFFLOOD, XFFROZEN,      &
                            XPIFLOOD, XCPL_EFLOOD, XCPL_PFLOOD, XCPL_IFLOOD,    &
                            XCPL_DRAIN, XCPL_RUNOFF, LGLACIER,                  &
                            LTEMP_ARP, NTEMPLAYER_ARP, XPSN, XPSNG, XPSNV,      &
                            XPSNV_A, XFF, XFFG, XFFV, XPCPS, XPLVTT, XPLSTT,    &
                            LCANOPY, LCANOPY_DRAG, XDIR_ALB_WITH_SNOW,          &
                            XSCA_ALB_WITH_SNOW, XALBF, XEMISF, XCPL_ICEFLUX,    &
                            XWR, XWGI, XAN, XANDAY, XANFM, XRESA, XICE_STO,     &
                            XCONDSAT_EXP, XEXP_DIF, XALPHA, XNVG, XMVG, XLVG,   &
                            XWRES
!
USE MODD_CH_ISBA_n,      ONLY : CSV, CCH_NAMES, NBEQ, NSV_CHSBEG, NSV_CHSEND, &
                                  CCHEM_SURF_FILE, NDSTEQ, NSV_DSTBEG, NSV_DSTEND,&
                                  NSV_AERBEG, NSV_AEREND, NAEREQ, CDSTNAMES, CAER_NAMES,&
                                  NSLTEQ, NSV_SLTBEG,  NSV_SLTEND, CSLTNAMES,   &
                                  LCH_BIO_FLUX, CCH_DRY_DEP  
USE MODD_CHS_AEROSOL,    ONLY: LVARSIGI, LVARSIGJ
USE MODD_DST_SURF,       ONLY: LVARSIG, CDSTYN, NDSTMDE,  NDST_MDEBEG, LRGFIX_DST
USE MODD_SLT_SURF,       ONLY: LVARSIG_SLT, CSLTYN, NSLTMDE,  NSLT_MDEBEG, LRGFIX_SLT

USE MODD_DST_n
USE MODD_SLT_n
USE MODD_DIAG_ISBA_n,    ONLY : LPATCH_BUDGET
USE MODD_DIAG_MISC_ISBA_n, ONLY : LSURF_DIAG_ALBEDO
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE
USE MODD_SURF_PAR,       ONLY : XUNDEF, NUNDEF
USE MODD_ASSIM
USE MODD_AGRI,           ONLY : LAGRIP, XTHRESHOLD
USE MODD_AGRI_n,         ONLY : NIRRINUM, XTHRESHOLDSPT, LIRRIDAY, LIRRIGATE
USE MODD_SGH_PAR,        ONLY : NDIMTAB, XICE_DEPH_MAX, XF_DECAY
USE MODD_DEEPSOIL,       ONLY : LPHYSDOMC, LDEEPSOIL, XTDEEP_CLI, XGAMMAT_CLI
USE MODD_CSTS,           ONLY : XCPD, XLVTT, XLSTT
USE MODD_SURF_ATM,       ONLY : LCPL_ARP, LCPL_ESM
!
!
USE MODI_INIT_IO_SURF_n
USE MODI_END_IO_SURF_n
!
USE MODI_READ_ISBA_n
USE MODI_INI_VAR_FROM_VEGTYPE_DATA
USE MODI_CONVERT_PATCH_ISBA
USE MODI_INIT_ISBA_MIXPAR
USE MODI_CH_INIT_NAMES
USE MODI_SUBSCALE_Z0EFF
USE MODI_DRY_WET_SOIL_ALBEDOS
USE MODI_INIT_SNOW_LW
USE MODI_AVERAGED_ALBEDO_EMIS_ISBA
USE MODI_CARBON_INIT
USE MODI_CO2_INIT_n
USE MODI_THRMCONDZ
USE MODI_HEATCAPZ
USE MODI_EMIS_FROM_VEG
USE MODI_VEG_FROM_LAI
USE MODI_Z0V_FROM_LAI
USE MODI_SURF_PATCH
USE MODI_WRITE_COVER_TEX_ISBA
USE MODI_WRITE_COVER_TEX_ISBA_PAR
USE MODI_DIAG_ISBA_INIT_n
USE MODI_CH_INIT_DEP_ISBA_n
USE MODI_OPEN_NAMELIST
USE MODI_CLOSE_NAMELIST
USE MODI_DST_INIT_NAMES
USE MODI_DST_INIT_MODES
USE MODI_SLT_INIT_NAMES
USE MODI_SLT_INIT_MODES
USE MODI_GET_LUOUT
USE MODI_INIT_TOP
USE MODI_EXP_DECAY_SOIL_FR
USE MODI_EXP_DECAY_SOIL_DIF
USE MODI_SET_ROUGH
USE MODI_READ_ISBA_CANOPY_n
!
USE MODI_SOILTEMP_ARP_PAR
!
USE MODE_SOIL
USE MODE_SOIL_VG
USE MODE_HYDRO_DIF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
USE MODI_GET_1D_MASK
USE MODI_INIT_DST_n
USE MODI_INIT_SLT_n
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
INTEGER           :: ILU      ! sizes of ISBA arrays
INTEGER           :: JILU     ! loop increment
INTEGER           :: ILUOUT   ! unit of output listing file
INTEGER           :: ICH      ! unit of input chemistry file
!
INTEGER           :: ISWB     ! number of shortwave spectral bands
!
INTEGER           :: IDECADE  ! decade of simulation
!
INTEGER :: JLAYER  ! loop counter on layers
INTEGER :: JPATCH  ! loop counter on tiles
INTEGER :: IRESP   ! return code
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZWG1 ! work array for surface water content
REAL, DIMENSION(:,:), ALLOCATABLE :: ZTG1 ! work array for surface temperature
REAL, DIMENSION(SIZE(PCO2))       :: ZCO2  ! CO2 concentration  (kg/kg)
REAL, DIMENSION(SIZE(PTSRAD))     :: ZTSRAD_NAT !radiative temperature
!
REAL, DIMENSION(:), ALLOCATABLE   :: ZM
REAL, DIMENSION(:,:), ALLOCATABLE :: ZF
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZTOPSOIL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!               Initialisation for IO
!
!
IF (LHOOK) CALL DR_HOOK('COMPUTE_ISBA_PARAMETERS',0,ZHOOK_HANDLE)
CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('COMPUTE_ISBA_PARAMETERS: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF
!
ISWB=SIZE(PSW_BANDS)
!
! initialization for I/O
!
CALL INIT_IO_SURF_n(HPROGRAM,'NATURE','ISBA  ','READ ')
!
!
!*       2.3    Physiographic data fields from land cover:
!               -----------------------------------------
!
!
ILU = SIZE(XCOVER,1)
ALLOCATE(XLAI       (ILU,NPATCH))
ALLOCATE(XVEG       (ILU,NPATCH))
ALLOCATE(XRSMIN     (ILU,NPATCH))
ALLOCATE(XGAMMA     (ILU,NPATCH))
ALLOCATE(XWRMAX_CF  (ILU,NPATCH))
ALLOCATE(XRGL       (ILU,NPATCH))
ALLOCATE(XCV        (ILU,NPATCH))
ALLOCATE(XDG        (ILU,NGROUND_LAYER,NPATCH))
ALLOCATE(XROOTFRAC  (ILU,NGROUND_LAYER,NPATCH))
ALLOCATE(XD_ICE     (ILU,NPATCH))
ALLOCATE(XZ0        (ILU,NPATCH))
ALLOCATE(XZ0_O_Z0H  (ILU,NPATCH))
ALLOCATE(XALBNIR_VEG(ILU,NPATCH))
ALLOCATE(XALBVIS_VEG(ILU,NPATCH))
ALLOCATE(XALBUV_VEG (ILU,NPATCH))
ALLOCATE(XEMIS      (ILU,NPATCH))
ALLOCATE(XVEGTYPE   (ILU,NVEGTYPE))
ALLOCATE(XGMES      (ILU,NPATCH))
ALLOCATE(XRE25      (ILU,NPATCH))
ALLOCATE(XBSLAI     (ILU,NPATCH))
ALLOCATE(XLAIMIN    (ILU,NPATCH))
ALLOCATE(XSEFOLD    (ILU,NPATCH))
ALLOCATE(XGC        (ILU,NPATCH))
ALLOCATE(XDMAX      (ILU,NPATCH))
ALLOCATE(XF2I       (ILU,NPATCH))
ALLOCATE(LSTRESS    (ILU,NPATCH))
ALLOCATE(XH_TREE    (ILU,NPATCH))
ALLOCATE(XCE_NITRO  (ILU,NPATCH))
ALLOCATE(XCF_NITRO  (ILU,NPATCH))
ALLOCATE(XCNA_NITRO (ILU,NPATCH))
ALLOCATE(XBSLAI_NITRO (ILU,NPATCH))
ALLOCATE(TSEED      (ILU,NPATCH))
ALLOCATE(TREAP      (ILU,NPATCH))
ALLOCATE(XWATSUP    (ILU,NPATCH))
ALLOCATE(XIRRIG     (ILU,NPATCH))
!
!
IF (TTIME%TDATE%MONTH /= NUNDEF) THEN
  IDECADE = 3 * ( TTIME%TDATE%MONTH - 1 ) + MIN(TTIME%TDATE%DAY-1,29) / 10 + 1
ELSE
  IDECADE = 1
END IF
!
CALL INIT_ISBA_MIXPAR(CISBA,XCOVER,CPHOTO,'NAT')
!
CALL CONVERT_PATCH_ISBA(CISBA,IDECADE,XCOVER,CPHOTO,                           &
                        'NAT',PVEG=XVEG,PLAI=XLAI,                             &
                        PRSMIN=XRSMIN,PGAMMA=XGAMMA,PWRMAX_CF=XWRMAX_CF,       &
                        PRGL=XRGL,PCV=XCV,PDG=XDG,PZ0=XZ0,PZ0_O_Z0H=XZ0_O_Z0H, &
                        PALBNIR_VEG=XALBNIR_VEG,PALBVIS_VEG=XALBVIS_VEG,       &
                        PALBUV_VEG=XALBUV_VEG,PEMIS_ECO=XEMIS,                 &
                        PVEGTYPE=XVEGTYPE,PROOTFRAC=XROOTFRAC,                 &
                        PGMES=XGMES,PBSLAI=XBSLAI,PLAIMIN=XLAIMIN,             &
                        PSEFOLD=XSEFOLD,PGC=XGC,                               &
                        PDMAX=XDMAX,PF2I=XF2I,OSTRESS=LSTRESS,PH_TREE=XH_TREE, &
                        PRE25=XRE25,PCE_NITRO=XCE_NITRO,PCF_NITRO=XCF_NITRO,   &
                        PCNA_NITRO=XCNA_NITRO,PD_ICE=XD_ICE,TPSEED=TSEED,      &
                        TPREAP=TREAP,PWATSUP=XWATSUP,PIRRIG=XIRRIG             ) 
!
IF (CISBA=='DIF') THEN
  DO JLAYER = 1, NGROUND_LAYER
    IF (ANY((XROOTFRAC(:,JLAYER,:)<0. .OR. XROOTFRAC(:,JLAYER,:)>1.) .AND. XPATCH(:,:).NE.0.)) &
      CALL ABOR1_SFX('COMPUTE_ISBA_PARAMETERS: WITH CISBA=DIF ROOTFRAC MUST BE DEFINED')
  ENDDO
ENDIF
!
!
!*       2.4    Fraction of each tile
!               ---------------------
!
ALLOCATE(XPATCH        (ILU,NPATCH))
ALLOCATE(XVEGTYPE_PATCH(ILU,NVEGTYPE,NPATCH))
CALL SURF_PATCH(NPATCH,XVEGTYPE,XPATCH,XVEGTYPE_PATCH)
!
!*       2.5    Masks for tiles
!               ---------------
!
ALLOCATE(NSIZE_NATURE_P (NPATCH))
DO JPATCH=1,NPATCH
  NSIZE_NATURE_P(JPATCH) = COUNT(XPATCH(:,JPATCH) > 0.0)
ENDDO
ALLOCATE(NR_NATURE_P(ILU,NPATCH))
NR_NATURE_P(:,:) = 0
DO JPATCH=1,NPATCH
  CALL GET_1D_MASK( NSIZE_NATURE_P(JPATCH), ILU, XPATCH(:,JPATCH), NR_NATURE_P(:NSIZE_NATURE_P(JPATCH),JPATCH))
ENDDO
!
!
!*       2.6    Miscellaneous fields for ISBA:
!               -----------------------------
!
!* default value for:
! lateral water flux, deep soil temperature climatology and its relaxation time-scale
!
ALLOCATE(XTDEEP (ILU))
ALLOCATE(XGAMMAT(ILU))
XTDEEP (:) = XUNDEF
XGAMMAT(:) = XUNDEF
!
IF (LDEEPSOIL) THEN
   DO JILU = 1, ILU
      !
      XTDEEP (JILU) = XTDEEP_CLI (TTIME%TDATE%MONTH)
      !
      XGAMMAT(JILU) = 1. / XGAMMAT_CLI(TTIME%TDATE%MONTH)
      !
   END DO
   !
   WRITE(ILUOUT,*)' LDEEPSOIL = ',LDEEPSOIL,' LPHYSDOMC = ',LPHYSDOMC
   WRITE(ILUOUT,*)' XTDEEP    = ',MINVAL(XTDEEP(:)),MAXVAL(XTDEEP(:))
   WRITE(ILUOUT,*)' XGAMMAT   = ',MINVAL(XGAMMAT(:)),MAXVAL(XGAMMAT(:))
   !
ENDIF
!
!
!*       2.7    Irrigation
!               ----------
!
IF (LAGRIP) THEN
   ALLOCATE(NIRRINUM(ILU,NPATCH))
   ALLOCATE(LIRRIDAY(ILU,NPATCH))
   ALLOCATE(LIRRIGATE(ILU,NPATCH))
   ALLOCATE(XTHRESHOLDSPT(ILU,NPATCH))
   !
   NIRRINUM (:,:) = 1
   LIRRIDAY (:,:) = .FALSE.                          
   LIRRIGATE(:,:) = .FALSE.                          
   !
   DO JILU = 1, ILU
      DO JPATCH = 1, NPATCH
         XTHRESHOLDSPT(JILU,JPATCH) = XTHRESHOLD(NIRRINUM(JILU,JPATCH))
      END DO
   END DO
END IF
!
!
!
!*       2.8    Additional fields for ISBA-AGS:
!               ------------------------------                        
!
IF(CPHOTO /= 'NON' .AND. HINIT == 'ALL') THEN
  ZCO2(:) = PCO2(:) / PRHOA(:)
  CALL CO2_INIT_n(ZCO2)
END IF
!
!
!*       2.9    Nitrogen version for isbaAgs
!               ------------------------------                        
!
IF (CPHOTO=='NIT' .OR. CPHOTO=='NCB') THEN
  WHERE ((XCE_NITRO (:,:)*XCNA_NITRO(:,:)+XCF_NITRO (:,:)) /= 0. )
      XBSLAI_NITRO(:,:) = 1. / (XCE_NITRO (:,:)*XCNA_NITRO(:,:)+XCF_NITRO (:,:))
  ELSEWHERE
      XBSLAI_NITRO(:,:) = XUNDEF
  ENDWHERE
ENDIF
!
!
!*       2.10   Soil carbon
!               -----------                        
!
IF (CRESPSL=='CNT' .AND. HINIT == 'ALL') THEN
  CALL CARBON_INIT(CPHOTO, CRESPSL, NNBIOMASS, NNLITTER, NNLITTLEVS, NNSOILCARB)
ENDIF
!
!-------------------------------------------------------------------------------
!
!        3.  Initialize Chemical Deposition
!            ------------------------------
!
!        3.1 Chemical gazes
!            --------------
!
IF (KSV /= 0) THEN
  ALLOCATE(CSV(KSV))
  CALL CH_INIT_NAMES(ILUOUT,HSV,NBEQ,NAEREQ,CSV,&
                       NSV_CHSBEG,NSV_CHSEND, &
                       NSV_AERBEG, NSV_AEREND,&
                       LVARSIGI, LVARSIGJ)  

  IF (NBEQ > 0 ) THEN
    ALLOCATE(CCH_NAMES(NBEQ))
    CCH_NAMES(:) = CSV(NSV_CHSBEG:NSV_CHSEND)
    CALL OPEN_NAMELIST(HPROGRAM,ICH,HFILE=CCHEM_SURF_FILE)
    CALL CH_INIT_DEP_ISBA_n(ICH,ILUOUT,CSV,ILU)
    CALL CLOSE_NAMELIST(HPROGRAM,ICH)

  ELSE
    ALLOCATE(CCH_NAMES(0))
  END IF

  IF (NAEREQ > 0 ) THEN
    ALLOCATE(CAER_NAMES(NAEREQ))
    CAER_NAMES(:) = CSV(NSV_AERBEG:NSV_AEREND)
  ELSE
    ALLOCATE(CAER_NAMES(0))
  END IF

   CALL DST_INIT_NAMES(         &
          ILUOUT,                   &!I [idx] index of writing unit
          HSV,                     &!I [char] list of scalar variables
          NDSTEQ,                  &!O [nbr] number of dust related tracers
          NSV_DSTBEG,              &!O [idx] first dust related scalar variable
          NSV_DSTEND,              &!O [idx] last dust related scalar variable
          LVARSIG,                 &!O type of standard deviation (fixed or variable)
          LRGFIX_DST,              &!O type of mean radius (fixed or variable)        
          CDSTYN                  &!O [char] Y/N, are dust related scalars present?
          )  

    IF (NDSTEQ >=1) THEN
   CALL DST_INIT_MODES(         &
          NDSTEQ,                   &!I [nbr] number of dust related variables in scalar list
          NSV_DSTBEG,              &!I [idx] index of first dust related variable in scalar list
          NSV_DSTEND,              &!I [idx] index of last dust related variable in scalar list
          LVARSIG,                 &!I type of standard deviation (fixed or variable)
          LRGFIX_DST,              &!O type of mean radius (fixed or variable)        
          NDST_MDEBEG,             &!O [idx] index of mass for first mode in scalar list
          NDSTMDE                 &!O [nbr] number of modes to be transported
          )  
     IF(.NOT. ASSOCIATED(CDSTNAMES)) ALLOCATE (CDSTNAMES(NDSTEQ))
     CDSTNAMES(:) = HSV(NSV_DSTBEG:NSV_DSTEND)
     ALLOCATE (XSFDST(ILU,NDSTEQ,NPATCH))  !Output array
     ALLOCATE (XSFDSTM(ILU,NDSTEQ,NPATCH))  !Output array
     XSFDST(:,:,:)  = 0.
     XSFDSTM(:,:,:) = 0.     
     CALL INIT_DST_n(HPROGRAM)    
    ELSE
     ALLOCATE(XSFDST(0,0,0))
     ALLOCATE(XSFDSTM(0,0,0))
    END IF

   CALL SLT_INIT_NAMES(         &
          ILUOUT,                   &!I [idx] index of writing unit
          HSV,                     &!I [char] list of scalar variables
          NSLTEQ,                  &!O [nbr] number of sea salt related tracers
          NSV_SLTBEG,              &!O [idx] first sea salt related scalar variable
          NSV_SLTEND,              &!O [idx] last sea salt related scalar variable
          LVARSIG_SLT,             &!O type of standard deviation (fixed or variable)
          LRGFIX_SLT,              &!O type of mean radius (fixed or variable)        
          CSLTYN                  &!O [char] Y/N, are sea salt related scalars present?
          )  

    IF (NSLTEQ >=1) THEN
   CALL SLT_INIT_MODES(         &
          NSLTEQ,                   &!I [nbr] number of sea salt related variables in scalar list
          NSV_SLTBEG,              &!I [idx] index of first sea salt related variable in scalar list
          NSV_SLTEND,              &!I [idx] index of last sea salt related variable in scalar list
          LVARSIG_SLT,             &!I type of standard deviation (fixed or variable)
          LRGFIX_SLT,              &!O type of mean radius (fixed or variable)
          NSLT_MDEBEG,             &!O [idx] index of mass for first mode in scalar list
          NSLTMDE                 &!O [nbr] number of modes to be transported
          )  
     IF(.NOT. ASSOCIATED(CSLTNAMES)) ALLOCATE (CSLTNAMES(NSLTEQ))
     CSLTNAMES(:) = HSV(NSV_SLTBEG:NSV_SLTEND)
     ALLOCATE (XSFSLT(ILU,NSLTEQ,NPATCH))  !Output array
     CALL INIT_SLT_n(HPROGRAM)    
    ELSE
     ALLOCATE(XSFSLT(0,0,0))
    END IF

ELSE
  ALLOCATE(CSV      (0))
  ALLOCATE(CDSTNAMES(0))
  ALLOCATE(CSLTNAMES(0))
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       4.     Orographic roughness length
!               ---------------------------
!
ALLOCATE(XZ0EFFIP(ILU,NPATCH))
ALLOCATE(XZ0EFFIM(ILU,NPATCH))
ALLOCATE(XZ0EFFJP(ILU,NPATCH))
ALLOCATE(XZ0EFFJM(ILU,NPATCH))
ALLOCATE(XZ0REL  (ILU))
!
CALL SUBSCALE_Z0EFF(XAOSIP,XAOSIM,XAOSJP,XAOSJM,         &
                      XHO2IP,XHO2IM,XHO2JP,XHO2JM,XZ0,     &
                      XZ0EFFIP,XZ0EFFIM,XZ0EFFJP,XZ0EFFJM, &
                      XZ0REL                               )  
!
!-------------------------------------------------------------------------------
!
!*       5.1     Soil hydraulic characteristics:
!                -------------------------------
!
ALLOCATE(XCONDSAT (ILU,NGROUND_LAYER,NPATCH))
!
IF(CDIF=='BC')THEN
!
! Note that if ISBA/=DIF, alway CDIF = 'BC' and CPEDOTF = 'CH78'
!        
  DO JLAYER=1,NGROUND_LAYER
     DO JPATCH=1,NPATCH
        XCONDSAT(:,JLAYER,JPATCH) = HYDCONDSAT_FUNC(XCLAY(:,JLAYER),XSAND(:,JLAYER),CPEDOTF)
     ENDDO
  END DO
!
ELSE  
!    
  ALLOCATE(ZTOPSOIL(ILU,NGROUND_LAYER))
  ZTOPSOIL(:,1)=1.0
  ZTOPSOIL(:,2)=1.0
  IF(NGROUND_LAYER>2)THEN
    DO JLAYER=3,NGROUND_LAYER
       ZTOPSOIL(:,JLAYER)=0.0
    ENDDO
  ENDIF
!
  DO JLAYER=1,NGROUND_LAYER
     DO JPATCH=1,NPATCH
        XCONDSAT(:,JLAYER,JPATCH) = HYDCONDSAT_VG(XCLAY(:,JLAYER),XSAND(:,JLAYER),XORGMAT(:,JLAYER), &
                                                  XDENSITY(:,JLAYER),ZTOPSOIL(:,JLAYER),CPEDOTF      )
     ENDDO
  END DO
!  
ENDIF
!
IF(CDIF=='BC')THEN
!        
! Note that if ISBA/=DIF, alway CDIF = 'BC' and CPEDOTF = 'CH78'
  ALLOCATE(XMPOTSAT(ILU,NGROUND_LAYER))
  ALLOCATE(XBCOEF  (ILU,NGROUND_LAYER))
  DO JLAYER=1,NGROUND_LAYER
     XBCOEF  (:,JLAYER) = BCOEF_FUNC     (XCLAY(:,JLAYER),XSAND(:,JLAYER),CPEDOTF)
     XMPOTSAT(:,JLAYER) = MATPOTSAT_FUNC (XCLAY(:,JLAYER),XSAND(:,JLAYER),CPEDOTF)
  END DO
  ALLOCATE(XALPHA(0,0))
  ALLOCATE(XNVG  (0,0))
  ALLOCATE(XMVG  (0,0))
  ALLOCATE(XLVG  (0,0)) 
  ALLOCATE(XWRES (0,0)) 
!  
ELSE 
!        
  ALLOCATE(XALPHA(ILU,NGROUND_LAYER))
  ALLOCATE(XNVG  (ILU,NGROUND_LAYER))
  ALLOCATE(XMVG  (ILU,NGROUND_LAYER))
  ALLOCATE(XLVG  (ILU,NGROUND_LAYER))
  ALLOCATE(XWRES (ILU,NGROUND_LAYER)) 
  DO JLAYER=1,NGROUND_LAYER
     XWRES (:,JLAYER) = WRES_VG (XCLAY(:,JLAYER),XSAND(:,JLAYER),CPEDOTF)
     XALPHA(:,JLAYER) = ALPHA_VG(XCLAY(:,JLAYER),XSAND(:,JLAYER),XORGMAT(:,JLAYER), &
                                 XDENSITY(:,JLAYER),ZTOPSOIL(:,JLAYER),CPEDOTF      )
     XNVG  (:,JLAYER) = NCOEF_VG(XCLAY(:,JLAYER),XSAND(:,JLAYER),XORGMAT(:,JLAYER), &
                                 XDENSITY(:,JLAYER),ZTOPSOIL(:,JLAYER),CPEDOTF      )      
     XLVG  (:,JLAYER) = LCOEF_VG(XCLAY(:,JLAYER),XSAND(:,JLAYER),XORGMAT(:,JLAYER), &
                                 XDENSITY(:,JLAYER),CPEDOTF                         )
     WHERE(XNVG(:,JLAYER)/=XUNDEF)
           XMVG  (:,JLAYER) = 1.0-1.0/XNVG(:,JLAYER)
     ELSEWHERE
           XMVG  (:,JLAYER) = XUNDEF
     ENDWHERE
     
  END DO
  ALLOCATE(XMPOTSAT(0,0))
  ALLOCATE(XBCOEF  (0,0))
!  
ENDIF
!
ALLOCATE(XWWILT(ILU,NGROUND_LAYER)) ! wilting point
ALLOCATE(XWFC  (ILU,NGROUND_LAYER)) ! field capacity
ALLOCATE(XWSAT (ILU,NGROUND_LAYER)) ! saturation
!
XCLAY=MAX(0.01,XCLAY)
!
IF(CDIF=='BC')THEN
!        
! Note that if ISBA/=DIF, alway CDIF = 'BC' and CPEDOTF = 'CH78'
  DO JLAYER=1,NGROUND_LAYER
     XWSAT (:,JLAYER) = WSAT_FUNC (XCLAY(:,JLAYER),XSAND(:,JLAYER),CPEDOTF)
     XWWILT(:,JLAYER) = WWILT_FUNC(XCLAY(:,JLAYER),XSAND(:,JLAYER),CPEDOTF)
     XWFC  (:,JLAYER) = WFC_FUNC  (XCLAY(:,JLAYER),XSAND(:,JLAYER),CPEDOTF)
  END DO   
!  
ELSE
!        
  DO JLAYER=1,NGROUND_LAYER
     XWSAT (:,JLAYER) = WSAT_VG (XCLAY(:,JLAYER),XSAND(:,JLAYER),XORGMAT(:,JLAYER), &
                                 XDENSITY(:,JLAYER),ZTOPSOIL(:,JLAYER),CPEDOTF      )
     XWWILT(:,JLAYER) = WWILT_VG(XWSAT(:,JLAYER),XALPHA(:,JLAYER),XNVG(:,JLAYER),   &
                                 XMVG(:,JLAYER),XWRES(:,JLAYER),CPEDOTF)
     XWFC  (:,JLAYER) = WFC_VG  (XWSAT(:,JLAYER),XCONDSAT(:,JLAYER,1),XNVG(:,JLAYER), &
                                 XWRES(:,JLAYER),CPEDOTF)
  END DO
!  
  DEALLOCATE(ZTOPSOIL)
!
ENDIF
!
ALLOCATE(XTAUICE  (ILU))
ALLOCATE(XRUNOFFD (ILU,NPATCH))
XTAUICE(:) = XTAU_ICE
!
IF (CISBA == 'DIF') THEN
  DO JPATCH=1,NPATCH
    XRUNOFFD(:,JPATCH) = XDG(:,NGROUND_LAYER-1,JPATCH)
  END DO
ELSE
  DO JPATCH=1,NPATCH
    XRUNOFFD(:,JPATCH) = XDG(:,2,JPATCH)
  END DO
ENDIF
!
IF (CISBA=='2-L' .OR. CISBA=='3-L') THEN
  ALLOCATE(XCGSAT (ILU))
  ALLOCATE(XC1SAT (ILU,NPATCH))
  ALLOCATE(XC2REF (ILU,NPATCH))
  ALLOCATE(XC3    (ILU,2,NPATCH))
  ALLOCATE(XC4B   (ILU))
  ALLOCATE(XACOEF (ILU))
  ALLOCATE(XPCOEF (ILU))
  ALLOCATE(XC4REF (ILU,NPATCH))
  XCGSAT(:)  = CGSAT_FUNC(XCLAY(:,1),XSAND(:,1))
  XC4B(:)    = C4B_FUNC(XCLAY(:,1))
  !
  XACOEF(:)  = ACOEF_FUNC(XCLAY(:,1))
  XPCOEF(:)  = PCOEF_FUNC(XCLAY(:,1))
  !
  DO JPATCH=1,NPATCH
    XC1SAT(:,JPATCH) = C1SAT_FUNC(XCLAY(:,1))
    XC2REF(:,JPATCH) = C2REF_FUNC(XCLAY(:,1))          
    XC4REF(:,JPATCH) = C4REF_FUNC(XCLAY(:,1),XSAND(:,1),       &
                                  XDG(:,2,            JPATCH), &
                                  XDG(:,NGROUND_LAYER,JPATCH)  )
    XC3     (:,1,JPATCH) = C3_FUNC(XCLAY(:,1))
    XC3     (:,2,JPATCH) = C3_FUNC(XCLAY(:,2))

  END DO
  !
  ALLOCATE(XCONDSAT_EXP(0,0,0))
  ALLOCATE(XEXP_DIF    (0,0,0))
  !
ELSE IF (CISBA=='DIF') THEN
  !
  ALLOCATE(XCGSAT (0))
  ALLOCATE(XC1SAT (0,0))
  ALLOCATE(XC2REF (0,0))
  ALLOCATE(XC3    (0,0,0))
  ALLOCATE(XC4B   (0))
  ALLOCATE(XC4REF (0,0))
  ALLOCATE(XACOEF (0))
  ALLOCATE(XPCOEF (0))
  !
  ALLOCATE(XCONDSAT_EXP (ILU,NGROUND_LAYER,NPATCH))
  ALLOCATE(XEXP_DIF     (ILU,NGROUND_LAYER,NPATCH))
  XCONDSAT_EXP(:,:,:) = XCONDSAT(:,:,:)
  XEXP_DIF    (:,:,:) = 1.0
  !
END IF
!
!*       5.2     Soil thermal characteristics:
!               --------------------------------
!
IF (CSCOND=='PL98'.OR.CISBA=='DIF') THEN
  ALLOCATE(XHCAPSOIL(ILU,NGROUND_LAYER))
  ! 
  CALL HEATCAPZ(XSAND,XWSAT,XHCAPSOIL)
ELSE
  ALLOCATE(XHCAPSOIL(0,0))
END IF

IF (CSCOND=='PL98') THEN
  ALLOCATE(XCONDDRY (ILU,NGROUND_LAYER))
  ALLOCATE(XCONDSLD (ILU,NGROUND_LAYER))
  ! 
  CALL THRMCONDZ(XSAND,XWSAT,XCONDDRY,XCONDSLD)
ELSE
  ALLOCATE(XCONDDRY (0,0))
  ALLOCATE(XCONDSLD (0,0))
END IF
!
!-------------------------------------------------------------------------------
ALLOCATE(XPCPS (ILU,NPATCH))
ALLOCATE(XPLVTT(ILU,NPATCH))
ALLOCATE(XPLSTT(ILU,NPATCH))
XPCPS (:,:) = XCPD
XPLVTT(:,:) = XLVTT
XPLSTT(:,:) = XLSTT
IF(CCPSURF=='DRY'.AND.LCPL_ARP) THEN
  CALL ABOR1_SFX('CCPSURF=DRY must not be used with LCPL_ARP')
ENDIF
!
!*       6.1    Initialize of the SGH scheme:'
!               ------------------------------
!
!Rainfall spatial distribution
!  
IF(CRAIN=='SGH')THEN
  ALLOCATE(XMUF(ILU))
  XMUF(:)=0.0
ELSE
  ALLOCATE(XMUF(0))
ENDIF
!
!Horton (also used by the flooding sheme)
! 
ALLOCATE(XKSAT_ICE(ILU,NPATCH))
!
IF(CISBA/='DIF')XD_ICE(:,:)=MIN(XDG(:,2,:),XD_ICE(:,:))
!
XD_ICE   (:,:) = MAX(XICE_DEPH_MAX,XD_ICE(:,:))
XKSAT_ICE(:,:) = XCONDSAT(:,1,:)
!
IF(CISBA=='DIF'.AND.CHORT=='SGH ')THEN
  CALL ABOR1_SFX('HORTON RUNOFF (CHORT) NOT IMPLEMENTED IN DIFFUSION SCHEME')
ENDIF
!  
!Topmodel
!  
ALLOCATE(ZM(ILU))
ALLOCATE(ZF(ILU,NPATCH))
!
ZM (:)   = XUNDEF
ZF (:,:) = XUNDEF
!
IF(CRUNOFF=='SGH ') THEN 
!
  IF(CISBA=='DIF')THEN
     CALL ABOR1_SFX('TOPMODEL (CRUNOFF=SGH) NOT IMPLEMENTED IN DIFFUSION SCHEME')
  ENDIF
!
  ALLOCATE(XFSAT(ILU))  
  XFSAT(:) = 0.0
!
  ALLOCATE(XTAB_FSAT(ILU,NDIMTAB))
  ALLOCATE(XTAB_WTOP(ILU,NDIMTAB))
!
  XTAB_FSAT(:,:) = 0.0
  XTAB_WTOP(:,:) = 0.0
!
  IF(HINIT/='PRE')THEN
!
    WHERE(XCLAY(:,1)==XUNDEF.AND.XTI_MEAN(:)/=XUNDEF) XTI_MEAN(:)=XUNDEF
!
    IF(CTOPREG/='DEF')THEN
       WRITE(ILUOUT,*)'!'
       WRITE(ILUOUT,*)'  YOU USE TOPMODEL WITHOUT THE REGRESSION    ' 
       WRITE(ILUOUT,*)' OF WOLOCK AND MCCABE (2000) (OPTION TOPREG) '
       WRITE(ILUOUT,*)'!'
    ENDIF
!      
    CALL INIT_TOP (CISBA, CTOPREG, ILUOUT, XPATCH, XDG, XWWILT(:,1),     &
                     XWSAT(:,1), XTI_MIN, XTI_MAX, XTI_MEAN, XTI_STD,      &
                     XTI_SKEW, XTAB_FSAT, XTAB_WTOP, ZM                    )  
!
  ENDIF
! 
ELSE                  
!
  ALLOCATE(XFSAT(0))
!  
  ALLOCATE(XTAB_FSAT(0,0))
  ALLOCATE(XTAB_WTOP(0,0))
!                  
ENDIF  
!
!Exponential decay
!  
IF(HINIT/='PRE'.AND.CKSAT=='SGH')THEN 
!
  IF(CRUNOFF=='SGH')THEN
!   Exponential decay factor calculate using soil properties (eq. 11, Decharme et al., J. Hydrometeor, 2006)
    DO JPATCH=1,NPATCH
       WHERE(ZM(:)/=XUNDEF)
             ZF(:,JPATCH)=(XWSAT(:,1)-XWWILT(:,1))/ZM(:)
       ELSEWHERE
             ZF(:,JPATCH)=4.0/XDG(:,2,JPATCH)
       ENDWHERE
       ZF(:,JPATCH)=MIN(ZF(:,JPATCH),XF_DECAY)
    ENDDO
!
  ELSE
!   Exponential decay factor calculate using soil depth (equivalent to Decharme et al.)
    IF(CISBA/='DIF')THEN
      ZF(:,:)=4.0/XDG(:,2,:)
      ZF(:,:)=MIN(ZF(:,:),XF_DECAY)
    ENDIF
  ENDIF
!
  IF(CISBA=='2-L' .OR. CISBA=='3-L') THEN
!
    DO JPATCH=1,NPATCH
       CALL EXP_DECAY_SOIL_FR(CISBA, ZF(:,JPATCH),XC1SAT(:,JPATCH),XC2REF(:,JPATCH), &
                                XDG(:,:,JPATCH),XD_ICE(:,JPATCH),XC4REF(:,JPATCH),   &
                                XC3(:,:,JPATCH),XCONDSAT(:,:,JPATCH),                &
                                XKSAT_ICE(:,JPATCH))  
    ENDDO
!
  ELSE
!     
    DO JPATCH=1,NPATCH
       CALL EXP_DECAY_SOIL_DIF(XDG(:,:,JPATCH),XROOTFRAC(:,:,JPATCH),ZF(:,JPATCH), &
                               XCONDSAT_EXP(:,:,JPATCH),XEXP_DIF(:,:,JPATCH)       )    
    ENDDO
!                
  ENDIF        
!                                    
ENDIF
!
DEALLOCATE(ZM)
DEALLOCATE(ZF)
!
!-------------------------------------------------------------------------------
!
!*       6.2    Initialize of TRIP or ESM coupling:'
!               ------------------------------------
!
IF(LCPL_ESM)THEN
   LTRIP=.TRUE.
   IF(.NOT.LGLACIER)THEN
     CALL ABOR1_SFX('LGLACIER MUST BE ACTIVATED WITH EARTH SYSTEM MODEL')
   ENDIF
ENDIF
!
IF(LGLACIER)THEN
   ALLOCATE(XCPL_ICEFLUX(ILU))
   XCPL_ICEFLUX(:) = 0.0
ELSE
   ALLOCATE(XCPL_ICEFLUX(0))
ENDIF
!
IF(LTRIP)THEN
!        
  ALLOCATE(XCPL_DRAIN (ILU))
  ALLOCATE(XCPL_RUNOFF(ILU))
  XCPL_DRAIN  = 0.0
  XCPL_RUNOFF = 0.0
!
  IF(LFLOOD)THEN
    !         
    IF(CISBA=='DIF')THEN
       CALL ABOR1_SFX('FLOOD SCHEME (LFLOOD) NOT IMPLEMENTED WITH DIFFUSION SCHEME')
    ENDIF
    !
    ALLOCATE(XFFLOOD      (ILU))
    ALLOCATE(XPIFLOOD     (ILU))
    ALLOCATE(XCPL_EFLOOD  (ILU))
    ALLOCATE(XCPL_PFLOOD  (ILU))
    ALLOCATE(XCPL_IFLOOD  (ILU))
    ALLOCATE(XFF          (ILU,NPATCH))
    ALLOCATE(XFFG         (ILU,NPATCH))
    ALLOCATE(XFFV         (ILU,NPATCH))  
    ALLOCATE(XFFROZEN     (ILU,NPATCH))  
    ALLOCATE(XALBF        (ILU,NPATCH))  
    ALLOCATE(XEMISF       (ILU,NPATCH))  
    XFFLOOD       = 0.0
    XPIFLOOD      = 0.0
    XCPL_EFLOOD   = 0.0
    XCPL_PFLOOD   = 0.0
    XCPL_IFLOOD   = 0.0
    XFF           = 0.0
    XFFG          = 0.0
    XFFV          = 0.0
    XFFROZEN      = 0.0
    XALBF         = 0.0
    XEMISF        = 0.0
  ELSE
    ALLOCATE(XFFLOOD      (0))
    ALLOCATE(XPIFLOOD     (0))
    ALLOCATE(XCPL_EFLOOD  (0))
    ALLOCATE(XCPL_PFLOOD  (0))
    ALLOCATE(XCPL_IFLOOD  (0))
    ALLOCATE(XFF        (0,0))
    ALLOCATE(XFFG       (0,0))
    ALLOCATE(XFFV       (0,0))
    ALLOCATE(XFFROZEN   (0,0))
    ALLOCATE(XALBF      (0,0))  
    ALLOCATE(XEMISF     (0,0))      
  ENDIF
  !
ELSE
!        
  ALLOCATE(XCPL_DRAIN (0))
  ALLOCATE(XCPL_RUNOFF(0))   
!  
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      7.     ISBA time-varying deep force-restore temperature initialization
!              ---------------------------------------------------------------
!
CALL SOILTEMP_ARP_PAR(HPROGRAM,LTEMP_ARP,NTEMPLAYER_ARP)
!
!-------------------------------------------------------------------------------
!
!*       8.     Physiographic Radiative fields:  
!               ------------------------------
!
!
!* dry and wet bare soil albedos
!
ALLOCATE(XALBNIR_DRY  (ILU))
ALLOCATE(XALBVIS_DRY  (ILU))
ALLOCATE(XALBUV_DRY   (ILU))
ALLOCATE(XALBNIR_WET  (ILU))
ALLOCATE(XALBVIS_WET  (ILU))
ALLOCATE(XALBUV_WET   (ILU))
!
CALL DRY_WET_SOIL_ALBEDOS(XSAND(:,1),XCLAY(:,1),                 &
                            XVEGTYPE,                              &
                            XALBNIR_DRY,XALBVIS_DRY,XALBUV_DRY,    &
                            XALBNIR_WET,XALBVIS_WET,XALBUV_WET     )  
!
!
!-------------------------------------------------------------------------------
!
!*       9.     Prints of cover parameters in a tex file
!               ----------------------------------------
!
CALL WRITE_COVER_TEX_ISBA    (NPATCH,NGROUND_LAYER,CISBA)
CALL WRITE_COVER_TEX_ISBA_PAR(NPATCH,NGROUND_LAYER,CISBA,CPHOTO)
!
!* if only physiographic fields are to be initialized, stop here.
!
IF (HINIT/='ALL') THEN
  CALL END_IO_SURF_n(HPROGRAM)
  IF (LHOOK) CALL DR_HOOK('COMPUTE_ISBA_PARAMETERS',1,ZHOOK_HANDLE)
  RETURN
END IF
!
!-------------------------------------------------------------------------------
!
!*      10.     Prognostic and semi-prognostic fields
!               -------------------------------------
!
CALL READ_ISBA_n(HPROGRAM)
!
!* z0 and vegetation fraction estimated from LAI
IF (CPHOTO=='LAI' .OR. CPHOTO=='LST' .OR. CPHOTO=='NIT' .OR. CPHOTO=='NCB') THEN
  DO JPATCH=1,NPATCH
     DO JILU=1,ILU    
        IF(XLAI(JILU,JPATCH)/=XUNDEF) THEN
           IF (XLAI(JILU,JPATCH).LT.XLAIMIN(JILU,JPATCH)) THEN
              XLAI(JILU,JPATCH)=XLAIMIN(JILU,JPATCH)
           ENDIF
           XZ0  (JILU,JPATCH) = Z0V_FROM_LAI(XLAI(JILU,JPATCH),XH_TREE(JILU,JPATCH),XVEGTYPE_PATCH(JILU,:,JPATCH))
           XVEG (JILU,JPATCH) = VEG_FROM_LAI(XLAI(JILU,JPATCH),XVEGTYPE_PATCH(JILU,:,JPATCH))
           XEMIS(JILU,JPATCH) = EMIS_FROM_VEG(XVEG(JILU,JPATCH),XVEGTYPE_PATCH(JILU,:,JPATCH))
        END IF  
     END DO
  END DO
END IF
!
!-------------------------------------------------------------------------------
!
!*      10bis.  Extrapolation of the prognostic and semi-prognostic fields
!                           LAND USE case 
!               -------------------------------------
!
IF (OLAND_USE) THEN
!
  DO JLAYER=1,SIZE(XTG(1,:,1))
    CALL INI_VAR_FROM_VEGTYPE_DATA(HPROGRAM,ILUOUT,'TEMP GRO', XTG(:,JLAYER,:))
  END DO
  DO JLAYER=1,NGROUND_LAYER
    CALL INI_VAR_FROM_VEGTYPE_DATA(HPROGRAM,ILUOUT,'WG      ', XWG(:,JLAYER,:))
    CALL INI_VAR_FROM_VEGTYPE_DATA(HPROGRAM,ILUOUT,'WGI     ', XWGI(:,JLAYER,:))
  END DO
  CALL INI_VAR_FROM_VEGTYPE_DATA(HPROGRAM,ILUOUT,'WR      ', XWR(:,:))
  CALL INI_VAR_FROM_VEGTYPE_DATA(HPROGRAM,ILUOUT,'RESA    ', XRESA(:,:))   
  CALL INI_VAR_FROM_VEGTYPE_DATA(HPROGRAM,ILUOUT,'ICE_STO ', XICE_STO(:,:))
  IF (CPHOTO=='LAI' .OR. CPHOTO=='LST' .OR. CPHOTO=='NIT' .OR. CPHOTO=='NCB') THEN
    CALL INI_VAR_FROM_VEGTYPE_DATA(HPROGRAM,ILUOUT,'AN      ', XAN(:,:))
    CALL INI_VAR_FROM_VEGTYPE_DATA(HPROGRAM,ILUOUT,'ANDAY   ', XANDAY(:,:))   
    CALL INI_VAR_FROM_VEGTYPE_DATA(HPROGRAM,ILUOUT,'ANFM    ', XANFM(:,:))
  END IF
!  
END IF
!
!-------------------------------------------------------------------------------
!
!*      11.     Canopy air fields:
!               -----------------
!
CALL READ_ISBA_CANOPY_n(HPROGRAM)
!
!-------------------------------------------------------------------------------
!
!*      12.     Roughness length option
!               -----------------------
!
CALL SET_ROUGH(LCANOPY,CROUGH)
!
!-------------------------------------------------------------------------------
!
!*      13.     Radiative fields and snow/flood fracion initialization:
!               -------------------------------------------------------
!
!* snow long-wave properties (not initialized in read_gr_snow)
!
CALL INIT_SNOW_LW(XEMISSN,TSNOW)
!
!* albedo per tile and averaged albedo, emissivity and radiative temperature
!
ALLOCATE(XALBNIR_SOIL(ILU,NPATCH))
ALLOCATE(XALBVIS_SOIL(ILU,NPATCH))
ALLOCATE(XALBUV_SOIL (ILU,NPATCH))
ALLOCATE(XALBNIR     (ILU,NPATCH))
ALLOCATE(XALBVIS     (ILU,NPATCH))
ALLOCATE(XALBUV      (ILU,NPATCH))
XALBNIR_SOIL(:,:) = XUNDEF
XALBVIS_SOIL(:,:) = XUNDEF
XALBUV_SOIL (:,:) = XUNDEF
XALBNIR     (:,:) = XUNDEF
XALBVIS     (:,:) = XUNDEF
XALBUV      (:,:) = XUNDEF
!
LSURF_DIAG_ALBEDO = .TRUE.
!
ALLOCATE(XEMIS_NAT   (ILU))
XEMIS_NAT (:) = XUNDEF
!
ALLOCATE(ZWG1(ILU,NPATCH))
ALLOCATE(ZTG1(ILU,NPATCH))
DO JPATCH=1,NPATCH
  ZWG1(:,JPATCH) = XWG(:,1,JPATCH)
  ZTG1(:,JPATCH) = XTG(:,1,JPATCH)
END DO
!
CALL CONVERT_PATCH_ISBA(CISBA,IDECADE,XCOVER,CPHOTO,'NAT',&
                          PWG1 = ZWG1, &
                          PALBNIR_SOIL=XALBNIR_SOIL, &
                          PALBVIS_SOIL=XALBVIS_SOIL, &
                          PALBUV_SOIL=XALBUV_SOIL )
!
!* Initialization of total albedo, emissivity and snow/flood fractions
!
ALLOCATE(XPSN (ILU,NPATCH))
ALLOCATE(XPSNG(ILU,NPATCH))
ALLOCATE(XPSNV(ILU,NPATCH))
XPSN  = 0.0
XPSNG = 0.0
XPSNV = 0.0
!
ALLOCATE(XDIR_ALB_WITH_SNOW(ILU,KSW,NPATCH))
ALLOCATE(XSCA_ALB_WITH_SNOW(ILU,KSW,NPATCH))
XDIR_ALB_WITH_SNOW = 0.0
XSCA_ALB_WITH_SNOW = 0.0
!
IF(TSNOW%SCHEME=='EBA')THEN
   ALLOCATE(XPSNV_A(ILU,NPATCH))
   XPSNV_A = 0.0
ELSE
   ALLOCATE(XPSNV_A(0,0))
ENDIF
!
PDIR_ALB = XUNDEF
PSCA_ALB = XUNDEF
PEMIS    = XUNDEF
PTSRAD   = XUNDEF
!
CALL AVERAGED_ALBEDO_EMIS_ISBA(LFLOOD, CALBEDO, PZENITH,                 &
                                 XVEG,XZ0,XLAI,ZTG1,                     &
                                 XPATCH,                                 &
                                 PSW_BANDS,                              &
                                 XALBNIR_VEG,XALBVIS_VEG,XALBUV_VEG,     &
                                 XALBNIR_SOIL,XALBVIS_SOIL,XALBUV_SOIL,  &
                                 XEMIS,                                  &
                                 TSNOW,                                  &
                                 XALBNIR,XALBVIS,XALBUV,                 &
                                 PDIR_ALB, PSCA_ALB,                     &
                                 XEMIS_NAT,ZTSRAD_NAT                    )  
!
PEMIS  = XEMIS_NAT
PTSRAD = ZTSRAD_NAT
!
DEALLOCATE(ZWG1)
DEALLOCATE(ZTG1)
!
!-------------------------------------------------------------------------------
!
!*      14.     ISBA diagnostics initialization
!               -------------------------------
!
IF(NPATCH<=1)LPATCH_BUDGET=.FALSE.
!
CALL DIAG_ISBA_INIT_n(HPROGRAM,ILU,ISWB)
!
!-------------------------------------------------------------------------------
!
!         End of IO
!
CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('COMPUTE_ISBA_PARAMETERS',1,ZHOOK_HANDLE)
!
END SUBROUTINE COMPUTE_ISBA_PARAMETERS
