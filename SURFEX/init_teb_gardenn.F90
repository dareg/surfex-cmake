!#############################################################
SUBROUTINE INIT_TEB_GARDEN_n(HPROGRAM,HINIT,                            &
                               KI,KSV,KSW,                                &
                               HSV,PCO2,PRHOA,                            &
                               PSW_BANDS,PDIR_ALB,PSCA_ALB,               &
                               PEMIS,PTSRAD,                              &
                               HATMFILE,HATMFILETYPE,                     &
                               HTEST                                      )  
!#############################################################
!
!!****  *INIT_TEB_GARDEN_n* - routine to initialize ISBA
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
!!	A. Lemonsu  *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    09/2009
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODD_DATA_TEB_GARDEN_n, ONLY : XDATA_VEGTYPE
USE MODD_READ_NAMELIST,     ONLY : LNAM_READ
!
USE MODD_TYPE_DATE_SURF
USE MODD_TYPE_SNOW
!
USE MODD_SNOW_PAR,        ONLY: XEMISSN
USE MODD_ISBA_PAR,        ONLY: XTAU_ICE
USE MODD_TEB_n,           ONLY: TTIME, XCOVER, XGARDEN
USE MODD_TEB_GARDEN_n,    ONLY: CROUGH,CISBA,CDIF,CPEDOTF,CPHOTO,CRUNOFF,CALBEDO,     &
                                  CSCOND, CC1DRY, CSOILFRZ, CDIFSFCOND, CSNOWRES,     &
                                  CRESPSL,NNBIOMASS, NNLITTER, NNLITTLEVS, NNSOILCARB,&
                                  XCLAY, XSAND, XWWILT, XWFC, XWSAT,                  &
                                  XVEG, XLAI, XRSMIN, XGAMMA, XRGL, XCV,              &
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
                                  XAOSIP,XAOSIM,XAOSJP,XAOSJM,                        &
                                  XHO2IP,XHO2IM,XHO2JP,XHO2JM,                        &
                                  XZ0EFFIP,XZ0EFFIM,XZ0EFFJP,XZ0EFFJM, XZ0REL,        &
                                  XVEGTYPE_PATCH,XROOTFRAC,XRUNOFFD,                  &
                                  XCGSAT, XC1SAT, XC2REF, XC3, XC4B, XACOEF, XPCOEF,  &
                                  XTAUICE, XACOEF, XPCOEF, XTAUICE, XBCOEF, XCONDSAT, &
                                  XHCAPSOIL, XCONDDRY, XCONDSLD, XC4REF, XMPOTSAT,    &
                                  XTDEEP, XGAMMAT, NGROUND_LAYER,                     &
                                  XCE_NITRO, XCF_NITRO,                               &
                                  XCNA_NITRO, XBSLAI_NITRO, CCPSURF, TSEED,           &
                                  TREAP, XWATSUP, XIRRIG, XCGMAX,                     &
                                  CKSAT, CTOPREG, CHORT,                              &
                                  XD_ICE, XKSAT_ICE, XPCPS, XPLVTT, XPLSTT, LCANOPY,  &
                                  XPSN, XPSNG, XPSNV, XPSNV_A, XCONDSAT_EXP, XEXP_DIF,&
                                  LPAR_GARDEN      
USE MODD_CH_TEB_n,        ONLY: CSV, CCH_NAMES, NBEQ, NSV_CHSBEG, NSV_CHSEND,         &
                                  CCHEM_SURF_FILE, NDSTEQ, NSV_DSTBEG, NSV_DSTEND,    &
                                  NSV_AERBEG, NSV_AEREND, NAEREQ, CDSTNAMES,          &
                                  CAER_NAMES, NSLTEQ, NSV_SLTBEG,                     &
                                  NSV_SLTEND, CSLTNAMES, CCH_DRY_DEP, LCH_BIO_FLUX  
USE MODD_CHS_AEROSOL,     ONLY: LVARSIGI, LVARSIGJ
USE MODD_DST_SURF,        ONLY: LVARSIG, CDSTYN, NDSTMDE,  NDST_MDEBEG, LRGFIX_DST
USE MODD_SLT_SURF,        ONLY: LVARSIG_SLT, CSLTYN, NSLTMDE,  NSLT_MDEBEG, LRGFIX_SLT
USE MODD_DIAG_MISC_TEB_n, ONLY: LSURF_DIAG_ALBEDO

USE MODD_DST_n
USE MODD_SLT_n
USE MODD_DATA_COVER_PAR,  ONLY: NVEGTYPE
USE MODD_SURF_PAR,        ONLY: XUNDEF, NUNDEF

USE MODD_ASSIM_GARDEN,    ONLY: LASSIM, CASSIM
USE MODD_AGRI_GARDEN,     ONLY: LAGRIP, XTHRESHOLD
USE MODD_AGRI_GARDEN_n,   ONLY: NIRRINUM, XTHRESHOLDSPT, LIRRIDAY, LIRRIGATE
USE MODD_SGH_PAR,         ONLY: NDIMTAB, XICE_DEPH_MAX, XF_DECAY
USE MODD_DEEPSOIL_GARDEN, ONLY: LPHYSDOMC, LDEEPSOIL, XTDEEP_CLI, XGAMMAT_CLI
USE MODD_CSTS,            ONLY: XCPD, XLVTT, XLSTT
USE MODD_SURF_ATM,        ONLY: LCPL_ARP, LCPL_ESM
USE MODD_SGH_PAR,         ONLY: NDIMTAB
!
USE MODN_TEB_GARDEN_n,    ONLY : XTSTEP
USE MODI_DEFAULT_ISBA
USE MODI_DEFAULT_CH_DEP
USE MODI_DEFAULT_CH_BIO_FLUX
USE MODI_READ_DEFAULT_TEB_GARDEN_n
USE MODI_ALLOCATE_TEB_GARDEN
USE MODI_ALLOC_DIAG_TEB_GARDEN
USE MODI_READ_TEB_GARDEN_CONF_n
USE MODI_READ_TEB_GARDEN_n
USE MODI_READ_PGD_TEB_GARDEN_n
USE MODI_CONVERT_PATCH_ISBA
USE MODI_INIT_FROM_DATA_GRDN_n
USE MODI_CH_INIT_NAMES
USE MODI_SUBSCALE_Z0EFF
USE MODI_DRY_WET_SOIL_ALBEDOS
USE MODI_SOIL_ALBEDO
USE MODI_INIT_SNOW_LW
USE MODI_AVG_ALBEDO_EMIS_GARDEN
USE MODI_THRMCONDZ
USE MODI_HEATCAPZ
USE MODI_EMIS_FROM_VEG
USE MODI_VEG_FROM_LAI
USE MODI_Z0V_FROM_LAI
USE MODI_SURF_PATCH
USE MODI_DIAG_TEB_GARDEN_INIT_n
USE MODI_READ_PREP_GARDEN_SNOW
USE MODI_CH_INIT_DEP_ISBA_n
USE MODI_OPEN_NAMELIST
USE MODI_CLOSE_NAMELIST
USE MODI_DST_INIT_NAMES
USE MODI_DST_INIT_MODES
USE MODI_SLT_INIT_NAMES
USE MODI_SLT_INIT_MODES
USE MODI_GET_LUOUT
USE MODI_EXP_DECAY_SOIL_FR
USE MODI_EXP_DECAY_SOIL_DIF
USE MODI_SET_ROUGH
!
USE MODE_SOIL

!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
USE MODI_GET_1D_MASK
!
USE MODI_INIT_DST_n
!
USE MODI_INIT_SLT_n
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=6),                   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
CHARACTER(LEN=3),                   INTENT(IN)  :: HINIT     ! choice of fields to initialize
INTEGER,                            INTENT(IN)  :: KI        ! number of points
INTEGER,                            INTENT(IN)  :: KSV       ! number of scalars
INTEGER,                            INTENT(IN)  :: KSW       ! number of short-wave spectral bands
CHARACTER(LEN=6), DIMENSION(KSV),   INTENT(IN)  :: HSV       ! name of all scalar variables
REAL,             DIMENSION(KI),    INTENT(IN)  :: PCO2      ! CO2 concentration (kg/m3)
REAL,             DIMENSION(KI),    INTENT(IN)  :: PRHOA     ! air density
REAL,             DIMENSION(KSW),   INTENT(IN)  :: PSW_BANDS ! middle wavelength of each band
REAL,             DIMENSION(KI,KSW),INTENT(OUT) :: PDIR_ALB  ! direct albedo for each band
REAL,             DIMENSION(KI,KSW),INTENT(OUT) :: PSCA_ALB  ! diffuse albedo for each band
REAL,             DIMENSION(KI),    INTENT(OUT) :: PEMIS     ! emissivity
REAL,             DIMENSION(KI),    INTENT(OUT) :: PTSRAD    ! radiative temperature
!
CHARACTER(LEN=28),                  INTENT(IN)  :: HATMFILE    ! atmospheric file name
CHARACTER(LEN=6),                   INTENT(IN)  :: HATMFILETYPE! atmospheric file type
CHARACTER(LEN=2),                   INTENT(IN)  :: HTEST       ! must be equal to 'OK'
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: JILU     ! loop increment
INTEGER           :: ILUOUT   ! unit of output listing file
INTEGER           :: ICH      ! unit of input chemistry file
!
!
INTEGER           :: IDECADE  ! decade of simulation
!
INTEGER :: JLAYER  ! loop counter on layers
INTEGER :: JPATCH  ! loop counter on tiles
INTEGER :: JVEGTYPE  ! loop counter on vegtypes
INTEGER :: IRESP   ! return code
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZWG1 ! work array for surface water content
REAL, DIMENSION(:,:), ALLOCATABLE :: ZTG1 ! work array for surface temperature
REAL, DIMENSION(SIZE(PCO2))       :: ZCO2  ! CO2 concentration  (kg/kg)
REAL, DIMENSION(:)  , ALLOCATABLE :: ZM
REAL, DIMENSION(:,:), ALLOCATABLE :: ZF
REAL                              :: ZOUT_TSTEP
!
CHARACTER(LEN=3)                  :: YRAIN 
LOGICAL                           :: GCANOPY_DRAG
LOGICAL                           :: GGLACIER
LOGICAL                           :: GTRIP
LOGICAL                           :: GFLOOD
LOGICAL                           :: GVEGUPD
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!               Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('INIT_TEB_GARDEN_N',0,ZHOOK_HANDLE)
CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('INIT_TEB_GARDEN_N: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF
!
!               Other little things
!
PDIR_ALB = XUNDEF
PSCA_ALB = XUNDEF
PEMIS    = XUNDEF
PTSRAD   = XUNDEF
!
LSURF_DIAG_ALBEDO = .FALSE.
!
IF (LNAM_READ) THEN
 !
 !*       0.     Defaults
 !               --------
 !
 !        0.1. Hard defaults
 !      
 !       Definition of default options for ISBA (in MODD_TEB_GARDEN_n)
 !       REM - TSTEP, OUT_TSTEP, CANOPY_DRAG are defined as local variables
 !             because they are already in init_teb.f90 (these options are 
 !             forced to the same values for TEB and urban green areas)
 !
 CALL DEFAULT_ISBA(XTSTEP, ZOUT_TSTEP,                           &
                     CROUGH,CRUNOFF,CALBEDO,CSCOND,              &
                     CC1DRY, CSOILFRZ, CDIFSFCOND, CSNOWRES,     &
                     CCPSURF, XCGMAX, CKSAT, CTOPREG,            &
                     YRAIN, CHORT, GFLOOD, GTRIP, GGLACIER,      &
                     GCANOPY_DRAG, GVEGUPD                       )  
 !
 CALL DEFAULT_CH_DEP(CCH_DRY_DEP)
 CALL DEFAULT_CH_BIO_FLUX(LCH_BIO_FLUX)
 !
ENDIF
!        0.2. Defaults from file header
!    
CALL READ_DEFAULT_TEB_GARDEN_n(HPROGRAM)
!
CALL READ_TEB_GARDEN_CONF_n(HPROGRAM)
!
!*       1.     Reading of configuration:
!               -------------------------
!
!* initialization of snow scheme (TSNOW defined in MODD_TEB_GARDEN_n)
!
IF (HINIT=='PRE') CALL READ_PREP_GARDEN_SNOW(HPROGRAM,TSNOW%SCHEME,TSNOW%NLAYER) 
!
!* other general options
!-------------------------------------------------------------------------------
!
!*       2.     Physiographic fields
!               --------------------
!
!
!*       2.1    Cover, soil and orographic fields:
!               ---------------------------------
!
CALL READ_PGD_TEB_GARDEN_n(HPROGRAM)
IF ( CPHOTO/='NON' .AND. NPATCH/=12) THEN
  CALL ABOR1_SFX('INIT_TEB_GARDEN_N: INCONSISTENCY BETWEEN CPHOTO AND NPATCH')
ENDIF   
IF ( CPHOTO/='LAI' .AND. CPHOTO/='LST' .AND. CPHOTO/='NIT' .AND. LAGRIP) THEN
  CALL ABOR1_SFX('INIT_TEB_GARDEN_N: INCONSISTENCY BETWEEN CPHOTO AND LAGRIP')
ENDIF
IF (HINIT=='PRE' .AND. TSNOW%SCHEME.NE.'3-L' .AND. TSNOW%SCHEME.NE.'CRO' .AND. CISBA=='DIF') THEN
    CALL ABOR1_SFX("INIT_TEB_GARDEN_N: WITH CISBA = DIF, CSNOW MUST BE 3-L OR CRO")
ENDIF
!
!
!* allocation of urban green area variables
!
CALL ALLOCATE_TEB_GARDEN(HINIT, KI, NVEGTYPE, NGROUND_LAYER,      &
                           NPATCH, KSW, NDIMTAB)  
!
!
!
!*       2.3    Physiographic data fields from land cover:
!               -----------------------------------------
!
IF (TTIME%TDATE%MONTH /= NUNDEF) THEN
  IDECADE = 3 * ( TTIME%TDATE%MONTH - 1 ) + MIN(TTIME%TDATE%DAY-1,29) / 10 + 1
ELSE
  IDECADE = 1
END IF
!
!
IF (.NOT. LPAR_GARDEN) THEN
  CALL CONVERT_PATCH_ISBA(CISBA,IDECADE,XCOVER,CPHOTO,                         &
                        'GRD',PVEG=XVEG,PLAI=XLAI,                             &
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
ELSE
 CALL INIT_FROM_DATA_GRDN_n(IDECADE,CPHOTO,                     &
                            XVEG,                               &
                            XLAI,XRSMIN,XGAMMA,XWRMAX_CF,       &
                            XRGL,XCV,XDG,XD_ICE,XZ0,XZ0_O_Z0H,  &
                            XALBNIR_VEG,XALBVIS_VEG,            &
                            XALBUV_VEG,XEMIS,                   &
                            XVEGTYPE,XROOTFRAC,                 &
                            XGMES,XBSLAI,XLAIMIN,XSEFOLD,XGC,   &
                            XDMAX, XF2I, LSTRESS, XH_TREE,XRE25,&
                            XCE_NITRO,XCF_NITRO,XCNA_NITRO      )  
END IF
!
WHERE (XGARDEN(:)==0.)
  XVEG(:,1)=0.
  XLAI(:,1)=0.
  XRSMIN(:,1)=40.
  XGAMMA(:,1)=0.
  XWRMAX_CF(:,1)=0.2
  XRGL(:,1)=100.
  XCV(:,1)=2.E-5
  XZ0(:,1)=0.013
  XZ0_O_Z0H(:,1)=10.
  XALBNIR_VEG(:,1)=0.30
  XALBVIS_VEG(:,1)=0.30
  XALBUV_VEG(:,1)=0.06
  XEMIS(:,1)=0.94
  XGMES(:,1)=0.020
  XBSLAI(:,1)=0.36
  XLAIMIN(:,1)=0.3
  XSEFOLD(:,1)=90*86400.
  XGC(:,1)=0.00025
  XDMAX(:,1)=0.1
  XF2I(:,1)=0.3
  XH_TREE(:,1)=0.
  XRE25(:,1)=3.6E-7
  XCE_NITRO(:,1)=7.68
  XCF_NITRO(:,1)=-4.33
  XCNA_NITRO(:,1)=1.3
END WHERE
DO JLAYER=1,NGROUND_LAYER
  WHERE (XGARDEN(:)==0.)
    XDG(:,JLAYER,1)=0.2*JLAYER
    XROOTFRAC(:,JLAYER,1)=0.2
  END WHERE
ENDDO
WHERE (XGARDEN(:)==0.) XD_ICE(:,1)=0.8*XDG(:,2,1)
DO JVEGTYPE=1,NVEGTYPE
  WHERE (XGARDEN(:)==0.)
    XVEGTYPE(:,JVEGTYPE)=0.
    XVEGTYPE(:,1)=1.
  END WHERE
ENDDO
!
!*       2.4    Fraction of each tile
!               ---------------------
!
CALL SURF_PATCH(NPATCH,XVEGTYPE,XPATCH,XVEGTYPE_PATCH)
!
!*       2.5    Masks for tiles (even if only ONE PATCH for the time being)
!               ---------------
!
DO JPATCH=1,NPATCH
  NSIZE_NATURE_P(JPATCH) = COUNT(XPATCH(:,JPATCH) > 0.0)
ENDDO
NR_NATURE_P(:,:) = 0
DO JPATCH=1,NPATCH
  CALL GET_1D_MASK( NSIZE_NATURE_P(JPATCH), KI, XPATCH(:,JPATCH), NR_NATURE_P(:NSIZE_NATURE_P(JPATCH),JPATCH))
ENDDO
!
!
!*       2.6    Miscellaneous fields for ISBA:
!               -----------------------------
!
!* default value for:
! lateral water flux, deep soil temperature climatology and its relaxation time-scale
!
XTDEEP (:) = XUNDEF
XGAMMAT(:) = XUNDEF
!
IF (LDEEPSOIL) THEN
   DO JILU = 1, KI
      XTDEEP (JILU) = XTDEEP_CLI (TTIME%TDATE%MONTH)
      XGAMMAT(JILU) = 1. / XGAMMAT_CLI(TTIME%TDATE%MONTH)
   END DO
ENDIF
!
!
!*       2.7    Irrigation
!               ----------
!
IF (LAGRIP) THEN
   !
   NIRRINUM (:,:) = 1
   LIRRIDAY (:,:) = .FALSE.                          
   LIRRIGATE(:,:) = .FALSE.                          
   !
   DO JILU = 1, KI
      DO JPATCH = 1, NPATCH
         XTHRESHOLDSPT(JILU,JPATCH) = XTHRESHOLD(NIRRINUM(JILU,JPATCH))
      END DO
   END DO
END IF
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
    CALL CH_INIT_DEP_ISBA_n(ICH,ILUOUT,CSV,KI)
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
     ALLOCATE (XSFDST(KI,NDSTEQ,NPATCH))  !Output array
     ALLOCATE (XSFDSTM(KI,NDSTEQ,NPATCH))  !Output array
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
     ALLOCATE (XSFSLT(KI,NSLTEQ,NPATCH))  !Output array
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
! Note that if ISBA/=DIF, alway CDIF = 'BC' and CPEDOTF = 'CH78'
  DO JLAYER=1,NGROUND_LAYER
     XBCOEF  (:,JLAYER) = BCOEF_FUNC     (XCLAY(:,JLAYER),XSAND(:,JLAYER),CPEDOTF)
     XMPOTSAT(:,JLAYER) = MATPOTSAT_FUNC (XCLAY(:,JLAYER),XSAND(:,JLAYER),CPEDOTF)
  END DO
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
! Van genuchten not yet implementd for garden
!  
ENDIF
!
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

IF (CISBA=='2-L' .OR. CISBA=='3-L') THEN
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

ELSE IF (CISBA=='DIF') THEN
  XCONDSAT_EXP(:,:,:) = XCONDSAT(:,:,:)
  XEXP_DIF    (:,:,:) = 1.0
END IF
!
!*       5.2     Soil thermal characteristics:
!               --------------------------------
!
IF (CSCOND=='PL98'.OR.CISBA=='DIF') THEN
  ! 
  CALL HEATCAPZ(XSAND,XWSAT,XHCAPSOIL)
ELSE
END IF

IF (CSCOND=='PL98') THEN
  CALL THRMCONDZ(XSAND,XWSAT,XCONDDRY,XCONDSLD)
ELSE
END IF
!
!-------------------------------------------------------------------------------
XPCPS(:,:) = XCPD
XPLVTT(:,:) = XLVTT
XPLSTT(:,:) = XLSTT
IF( CCPSURF=='DRY' .AND. LCPL_ARP ) THEN
  CALL ABOR1_SFX('CCPSURF=DRY must not be used with LCPL_ARP')
ENDIF
!
!*       6.1    Initialize of the SGH scheme:'
!               ------------------------------
!
!Horton (also used by the flooding sheme)
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
!Exponential decay
! 
ALLOCATE(ZM(KI))
ALLOCATE(ZF(KI,NPATCH))
!
ZM (:)   = XUNDEF
ZF (:,:) = XUNDEF
!
IF(HINIT/='PGD'.AND.HINIT/='PRE'.AND.CKSAT=='SGH')THEN 
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
                                XDG(:,:,JPATCH),XD_ICE(:,JPATCH),XC4REF(:,JPATCH),     &
                                XC3(:,:,JPATCH),XCONDSAT(:,:,JPATCH),      &
                                XKSAT_ICE(:,JPATCH))  
    ENDDO
!
  ELSEIF(CISBA=='DIF')THEN
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
!*       8.     Physiographic Radiative fields:  
!               ------------------------------
!
!
!* dry and wet bare soil albedos
!
CALL DRY_WET_SOIL_ALBEDOS(XSAND(:,1),XCLAY(:,1),                 &
                            XVEGTYPE,                              &
                            XALBNIR_DRY,XALBVIS_DRY,XALBUV_DRY,    &
                            XALBNIR_WET,XALBVIS_WET,XALBUV_WET     )  
!
!
!-------------------------------------------------------------------------------
!
IF (HINIT/='ALL') THEN
  IF (LHOOK) CALL DR_HOOK('INIT_TEB_GARDEN_N',1,ZHOOK_HANDLE)      
  RETURN
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      10.     Prognostic and semi-prognostic fields
!               -------------------------------------
!
!* allocation of urban green area variables
!
!
CALL READ_TEB_GARDEN_n(HPROGRAM)
!
CALL ALLOC_DIAG_TEB_GARDEN(KI,NGROUND_LAYER,KSW)
!
!* z0 and vegetation fraction estimated from LAI
IF (CPHOTO=='LAI' .OR. CPHOTO=='LST' .OR. CPHOTO=='NIT' .OR. CPHOTO=='NCB') THEN
  DO JPATCH=1,NPATCH
     DO JILU=1,KI
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
!*      12.     Roughness length option
!               -----------------------
!
CALL SET_ROUGH(LCANOPY,CROUGH)
!
!-------------------------------------------------------------------------------
!
!*      13.     Radiative fields:
!               ----------------
!
!* snow long-wave properties (not initialized in read_gr_snow)
!
CALL INIT_SNOW_LW(XEMISSN,TSNOW)
!
!* albedo per tile and averaged albedo, emissivity and radiative temperature
!
XALBNIR_SOIL(:,:) = XUNDEF
XALBVIS_SOIL(:,:) = XUNDEF
XALBUV_SOIL (:,:) = XUNDEF
XALBNIR     (:,:) = XUNDEF
XALBVIS     (:,:) = XUNDEF
XALBUV      (:,:) = XUNDEF
!
LSURF_DIAG_ALBEDO = .TRUE.
!
ALLOCATE(ZWG1(KI,NPATCH))
ALLOCATE(ZTG1(KI,NPATCH))
DO JPATCH=1,NPATCH
  ZWG1(:,JPATCH) = XWG(:,1,JPATCH)
  ZTG1(:,JPATCH) = XTG(:,1,JPATCH)
END DO
!
IF (.NOT. LPAR_GARDEN) THEN
  CALL SOIL_ALBEDO(CALBEDO,                               &
                     XWSAT(:,1),ZWG1,                       &
                     XALBVIS_DRY,XALBNIR_DRY,XALBUV_DRY,    &
                     XALBVIS_WET,XALBNIR_WET,XALBUV_WET,    &
                     XALBVIS_SOIL,XALBNIR_SOIL,XALBUV_SOIL  )  
ELSE
  CALL INIT_FROM_DATA_GRDN_n(IDECADE,CPHOTO,              &
                               PALBNIR_SOIL=XALBNIR_SOIL,   &
                               PALBVIS_SOIL=XALBVIS_SOIL,   &
                               PALBUV_SOIL=XALBUV_SOIL      )  
END IF
!
!
!* Initialization of total albedo, emissivity and snow/flood fractions
!
XPSN  = 0.0
XPSNG = 0.0
XPSNV = 0.0
!
!XDIR_ALB_WITH_SNOW = 0.0
!XSCA_ALB_WITH_SNOW = 0.0
!
IF(TSNOW%SCHEME=='EBA')THEN
   XPSNV_A = 0.0
ENDIF
! 
CALL AVG_ALBEDO_EMIS_GARDEN(CALBEDO,                                   &
                                 XVEG,XZ0,XLAI,ZTG1,                     &
                                 XPATCH,                                 &
                                 PSW_BANDS,                              &
                                 XALBNIR_VEG,XALBVIS_VEG,XALBUV_VEG,     &
                                 XALBNIR_SOIL,XALBVIS_SOIL,XALBUV_SOIL,  &
                                 XEMIS,                                  &
                                 TSNOW,                                  &
                                 XALBNIR,XALBVIS,XALBUV,                 &
                                 PDIR_ALB, PSCA_ALB,                     &
                                 PEMIS,PTSRAD                            )  
!
!
DEALLOCATE(ZWG1)
DEALLOCATE(ZTG1)
!
!-------------------------------------------------------------------------------
!
!*      14.     ISBA diagnostics initialization
!               -------------------------------
!
CALL DIAG_TEB_GARDEN_INIT_n(HPROGRAM,KI,KSW)
IF (LHOOK) CALL DR_HOOK('INIT_TEB_GARDEN_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE INIT_TEB_GARDEN_n
