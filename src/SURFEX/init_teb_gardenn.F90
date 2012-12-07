!#############################################################
SUBROUTINE INIT_TEB_GARDEN_n(HPROGRAM,HINIT,                          &
                             KI,KSV,KSW,                              &
                             HSV,PCO2,PRHOA,                          &
                             PSW_BANDS,PDIR_ALB,PSCA_ALB,             &
                             PEMIS,PTSRAD,                            &
                             HATMFILE,HATMFILETYPE,                   &
                             HTEST                                    )  
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
!!      B. Decharme 07/2011 : read pgd+prep
!!       
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_TEB_n,           ONLY: TTIME, XTSTEP, XCOVER, XGARDEN
USE MODD_TEB_GARDEN_n,    ONLY: CROUGH, CISBA, CPEDOTF, CPHOTO, CRUNOFF, CALBEDO,   &
                                CSCOND, CC1DRY, CSOILFRZ, CDIFSFCOND, CSNOWRES,     &
                                CRESPSL, LTR_ML, NNBIOMASS, NNLITTER, NNLITTLEVS,   &
                                NNSOILCARB, XCLAY, XSAND, XWWILT, XWFC, XW33, XWSAT,&
                                XVEG, XLAI, XRSMIN, XGAMMA, XRGL, XCV,              &
                                XDG, NWG_LAYER, XDROOT, XDG2, XDZG, XDZDIF,         &
                                XSOILGRID, XZ0, XZ0_O_Z0H,                          &
                                XALBNIR_VEG, XALBVIS_VEG, XALBUV_VEG, XQDGMES,      &
                                XEMIS, XVEGTYPE, XGMES, XRE25, XBSLAI, XLAIMIN, XGC,&
                                XDMAX, LSTRESS, XF2I, XAH, XAMAX, XANMAX, XBH,      &
                                XEPSO, XFZERO, XGAMM, XINCREASE, XQDAMAX, XQDGAMM,  &
                                XT1AMAX, XT1GMES, XT2AMAX, XT2GMES, XTAU_WOOD,      &
                                XSEFOLD, XH_TREE, XPATCH, NPATCH, XWRMAX_CF,        &
                                NR_NATURE_P, NSIZE_NATURE_P, XTURNOVER,             &
                                XALBNIR_DRY, XALBVIS_DRY, XALBUV_DRY,               &
                                XALBNIR_WET, XALBVIS_WET, XALBUV_WET,               &
                                XALBNIR_SOIL, XALBVIS_SOIL, XALBUV_SOIL,            &
                                XWG, XTG, TSNOW, XALBNIR, XALBVIS, XALBUV,          &
                                XAOSIP,XAOSIM,XAOSJP,XAOSJM,                        &
                                XHO2IP,XHO2IM,XHO2JP,XHO2JM,                        &
                                XZ0EFFIP,XZ0EFFIM,XZ0EFFJP,XZ0EFFJM, XZ0REL,        &
                                XVEGTYPE_PATCH,XROOTFRAC,XRUNOFFD, XSOILWGHT,       &
                                XCGSAT, XC1SAT, XC2REF, XC3, XC4B, XACOEF, XPCOEF,  &
                                XTAUICE, XBCOEF, XCONDSAT,                          &
                                XHCAPSOIL, XCONDDRY, XCONDSLD, XC4REF, XMPOTSAT,    &
                                XTDEEP, XGAMMAT, NGROUND_LAYER, XABC, XPOI, XMUS,   &
                                XCE_NITRO, XCF_NITRO, XFAPARC, XFAPIRC, XLAI_EFFC,  &
                                XCNA_NITRO, XBSLAI_NITRO, CCPSURF, TSEED,           &
                                TREAP, XWATSUP, XIRRIG, XCGMAX, XCDRAG,             &
                                CKSAT, CSOC, CTOPREG, CHORT,                        &
                                XD_ICE, XKSAT_ICE, XPCPS, XPLVTT, XPLSTT, LCANOPY,  &
                                XPSN, XPSNG, XPSNV, XPSNV_A, LPAR_GARDEN,           &
                                NLAYER_HORT, NLAYER_DUN, LSPINUPCARBS,              &
                                LSPINUPCARBW, NNBYEARSOLD, NSPINS, NSPINW
!
USE MODD_CH_TEB_n,        ONLY: CSV, CCH_NAMES, NBEQ, NSV_CHSBEG, NSV_CHSEND,       &
                                CCHEM_SURF_FILE, NDSTEQ, NSV_DSTBEG, NSV_DSTEND,    &
                                NSV_AERBEG, NSV_AEREND, NAEREQ, CDSTNAMES,          &
                                CAER_NAMES, NSLTEQ, NSV_SLTBEG,                     &
                                NSV_SLTEND, CSLTNAMES, CCH_DRY_DEP, LCH_BIO_FLUX 
!
USE MODD_DST_n
USE MODD_SLT_n
!
USE MODD_DEEPSOIL_GARDEN, ONLY: LPHYSDOMC, LDEEPSOIL, XTDEEP_CLI, XGAMMAT_CLI
USE MODD_AGRI_GARDEN,     ONLY: LAGRIP, XTHRESHOLD
USE MODD_AGRI_GARDEN_n,   ONLY: NIRRINUM, XTHRESHOLDSPT, LIRRIDAY, LIRRIGATE
!
USE MODD_DIAG_MISC_TEB_n, ONLY: LSURF_DIAG_ALBEDO
!
USE MODD_SGH_PAR,         ONLY: NDIMTAB, XICE_DEPH_MAX, XF_DECAY
!
USE MODD_DATA_COVER_PAR,  ONLY: NVEGTYPE
USE MODD_SURF_PAR,        ONLY: XUNDEF, NUNDEF
USE MODD_SNOW_PAR,        ONLY: XEMISSN
!
USE MODD_READ_NAMELIST,     ONLY : LNAM_READ
!
USE MODI_GET_LUOUT
USE MODI_ABOR1_SFX
USE MODI_INIT_IO_SURF_n
USE MODI_DEFAULT_ISBA
USE MODI_DEFAULT_CH_DEP
USE MODI_DEFAULT_CH_BIO_FLUX
USE MODI_READ_DEFAULT_TEB_GARDEN_n
USE MODI_READ_TEB_GARDEN_CONF_n
USE MODI_READ_PREP_GARDEN_SNOW
USE MODI_SET_SURFEX_FILEIN
USE MODI_READ_PGD_TEB_GARDEN_n
USE MODI_END_IO_SURF_n
!
USE MODI_ALLOCATE_PHYSIO
USE MODI_CONVERT_PATCH_ISBA
USE MODI_INIT_FROM_DATA_GRDN_n
USE MODI_COMMON_PARTS
USE MODI_EXP_DECAY_SOIL_DIF
USE MODI_EXP_DECAY_SOIL_FR
!
USE MODI_READ_TEB_GARDEN_n
USE MODI_COMMON_PARTS2
USE MODI_SOIL_ALBEDO
USE MODI_AVG_ALBEDO_EMIS_GARDEN
USE MODI_ALLOC_DIAG_TEB_GARDEN
USE MODI_DIAG_TEB_GARDEN_INIT_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
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
!
INTEGER           :: IDECADE  ! decade of simulation
!
INTEGER :: JLAYER   ! loop counter on layers
INTEGER :: JPATCH   ! loop counter on tiles
INTEGER :: JVEGTYPE
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZWG1 ! work array for surface water content
REAL, DIMENSION(:,:), ALLOCATABLE :: ZTG1 ! work array for surface temperature
!
REAL, DIMENSION(:), ALLOCATABLE   :: ZWORK
REAL, DIMENSION(:,:), ALLOCATABLE :: ZF
!
REAL                              :: ZOUT_TSTEP
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
                     CCPSURF, XCGMAX, XCDRAG, CKSAT, CSOC,       &
                     CTOPREG, YRAIN, CHORT, GFLOOD, GTRIP,       &
                     GGLACIER, GCANOPY_DRAG, GVEGUPD,            &
                     LSPINUPCARBS, LSPINUPCARBW                  )
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
!*       2.     Physiographic fields
!               --------------------
!
!
!*       2.1    Cover, soil and orographic fields:
!               ---------------------------------
!
CALL READ_PGD_TEB_GARDEN_n(HPROGRAM)
IF ( CPHOTO/='LAI' .AND. CPHOTO/='LST' .AND. CPHOTO/='NIT' .AND. CPHOTO/='NCB' .AND. LAGRIP) THEN
  CALL ABOR1_SFX('INIT_TEB_GARDEN_N: INCONSISTENCY BETWEEN CPHOTO AND LAGRIP')
ENDIF
IF (HINIT=='PRE' .AND. TSNOW%SCHEME.NE.'3-L' .AND. TSNOW%SCHEME.NE.'CRO' .AND. CISBA=='DIF') THEN
    CALL ABOR1_SFX("INIT_TEB_GARDEN_N: WITH CISBA = DIF, CSNOW MUST BE 3-L OR CRO")
ENDIF
!
!         End of IO
!
CALL END_IO_SURF_n(HPROGRAM)
IF (HINIT=='PGD') THEN
  IF (LHOOK) CALL DR_HOOK('INIT_TEB_GARDEN_N',1,ZHOOK_HANDLE)      
  RETURN
ENDIF
CALL SET_SURFEX_FILEIN(HPROGRAM,'PREP') ! restore input file name
!
!-----------------------------------------------------------------------------------------------------
! END READ PGD FILE
!-----------------------------------------------------------------------------------------------------
!  
!* allocation of urban green area variables
!
CALL ALLOCATE_PHYSIO(CPHOTO, CISBA, KI, NVEGTYPE, NGROUND_LAYER, NPATCH, &
                     XVEGTYPE, XLAI, XVEG, XZ0, XEMIS, XDG, XD_ICE,      &
                     XRSMIN, XGAMMA, XWRMAX_CF, XRGL, XCV,               &
                     XZ0_O_Z0H, XALBNIR_VEG, XALBVIS_VEG, XALBUV_VEG,    &
                     XH_TREE, XRE25, XLAIMIN, XBSLAI, XSEFOLD,           &
                     XGMES, XGC, XF2I, XDMAX, LSTRESS,                   &
                     XCE_NITRO, XCF_NITRO, XCNA_NITRO,                   &
                     TSEED, TREAP, XWATSUP, XIRRIG,                      &
                     XROOTFRAC, NWG_LAYER, XDROOT, XDG2                  )
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
!-----------------------------------------------------------------------------------------------------
!
IF (.NOT. LPAR_GARDEN) THEN        
  CALL CONVERT_PATCH_ISBA(CISBA,IDECADE,IDECADE,XCOVER,CPHOTO,LAGRIP,          &
                        'GRD',PVEG=XVEG,PLAI=XLAI,                             &
                        PRSMIN=XRSMIN,PGAMMA=XGAMMA,PWRMAX_CF=XWRMAX_CF,       &
                        PRGL=XRGL,PCV=XCV,PSOILGRID=XSOILGRID,                 &
                        PDG=XDG,KWG_LAYER=NWG_LAYER,PDROOT=XDROOT,PDG2=XDG2,   &
                        PZ0=XZ0,PZ0_O_Z0H=XZ0_O_Z0H,                           &
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
!
!-----------------------------------------------------------------------------------------------------
!

  IF (CISBA=='DIF') THEN
    DO JPATCH=1,NPATCH
      WHERE(XGARDEN(:)/=0.)
        NWG_LAYER(:,JPATCH)=NGROUND_LAYER 
        XDG2  (:,JPATCH)=0.0
        XDROOT(:,JPATCH)=0.0
      ENDWHERE
      DO JLAYER=NGROUND_LAYER,1,-1
        DO JILU=1,KI
          IF(XGARDEN(JILU)/=0..AND.XROOTFRAC(JILU,JLAYER,JPATCH)>=1.0)THEN
            XDG2  (JILU,JPATCH)=XDG(JILU,JLAYER,JPATCH)
            XDROOT(JILU,JPATCH)=XDG(JILU,JLAYER,JPATCH)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDIF
!
END IF
!
WHERE (XGARDEN(:)==0.)
  XVEG       (:,1)=0.
  XLAI       (:,1)=0.
  XRSMIN     (:,1)=40.
  XGAMMA     (:,1)=0.
  XWRMAX_CF  (:,1)=0.2
  XRGL       (:,1)=100.
  XCV        (:,1)=2.E-5
  XZ0        (:,1)=0.013
  XZ0_O_Z0H  (:,1)=10.
  XALBNIR_VEG(:,1)=0.30
  XALBVIS_VEG(:,1)=0.30
  XALBUV_VEG (:,1)=0.06
  XEMIS      (:,1)=0.94
ENDWHERE
IF (CPHOTO/='NON') THEN
  WHERE (XGARDEN(:)==0.)
    XH_TREE(:,1)=0.    
    XRE25  (:,1)=3.6E-7
    XLAIMIN(:,1)=0.3
    XBSLAI (:,1)=0.36
    XSEFOLD(:,1)=90*86400.
    XGMES  (:,1)=0.020
    XGC    (:,1)=0.00025
  END WHERE
  IF (CPHOTO/='AGS' .AND. CPHOTO/='LAI') THEN
    WHERE (XGARDEN(:)==0.)
      XF2I (:,1)=0.3          
      XDMAX(:,1)=0.1
    END WHERE
    IF (CPHOTO=='NIT' .OR. CPHOTO=='NCB') THEN
      WHERE (XGARDEN(:)==0.)
        XCE_NITRO (:,1)=7.68
        XCF_NITRO (:,1)=-4.33
        XCNA_NITRO(:,1)=1.3
      END WHERE
    ENDIF
  ENDIF
ENDIF
IF(CISBA/='DIF')THEN
  DO JLAYER=1,NGROUND_LAYER
    WHERE (XGARDEN(:)==0.)
      XDG(:,JLAYER,1)=0.2*JLAYER
    END WHERE
  ENDDO
ELSE
  WHERE (XGARDEN(:)==0.)
      XDG(:,1,1)=0.01
      XDG(:,2,1)=0.04
      XROOTFRAC(:,1,1)=0.
      XROOTFRAC(:,2,1)=0.
  END WHERE        
  DO JLAYER=3,NGROUND_LAYER
    WHERE (XGARDEN(:)==0.)
      XDG(:,JLAYER,1)=0.1*(JLAYER-2)
      XROOTFRAC(:,JLAYER,1)=0.
    END WHERE
  ENDDO        
  WHERE (XGARDEN(:)==0.) 
       NWG_LAYER(:,1)=NGROUND_LAYER
       XDROOT   (:,1)=0.0
       XDG2     (:,1)=XDG(:,NGROUND_LAYER-1,1)
  ENDWHERE
ENDIF  
WHERE (XGARDEN(:)==0.) 
   XD_ICE   (:,1)=0.8*XDG(:,2,1)
END WHERE
DO JVEGTYPE=1,NVEGTYPE
  WHERE (XGARDEN(:)==0.)
    XVEGTYPE(:,JVEGTYPE)=0.
    XVEGTYPE(:,1)=1.
  END WHERE
ENDDO
!
!* initialization of carbon scheme
!
NNBYEARSOLD = 0
NSPINS      = 1
NSPINW      = 1
!
LSPINUPCARBW = .FALSE.
LSPINUPCARBS = .FALSE.
!-------------------------------------------------------------------------------
!
CALL COMMON_PARTS(HPROGRAM, ILUOUT, KI, NPATCH, NGROUND_LAYER, TTIME%TDATE%MONTH,   &
                  XVEGTYPE, XPATCH, XVEGTYPE_PATCH, NSIZE_NATURE_P, NR_NATURE_P,    &
                  LDEEPSOIL, LPHYSDOMC, XTDEEP_CLI, XGAMMAT_CLI, XTDEEP, XGAMMAT,   &
                  LAGRIP, XTHRESHOLD, NIRRINUM, LIRRIDAY, LIRRIGATE, XTHRESHOLDSPT, &
                  CPHOTO, HINIT, LTR_ML, NNBIOMASS, PCO2, PRHOA, XABC, XPOI,  &
                  XGMES, XGC, XDMAX, XANMAX, XFZERO, XEPSO, XGAMM, XQDGAMM,   &
                  XQDGMES, XT1GMES, XT2GMES, XAMAX, XQDAMAX, XT1AMAX, XT2AMAX,&
                  XAH, XBH, XTAU_WOOD, XINCREASE, XTURNOVER,                  &
                  KSV, HSV, NBEQ, CSV, NAEREQ, NSV_CHSBEG, NSV_CHSEND,        &
                  NSV_AERBEG, NSV_AEREND, CCH_NAMES, CAER_NAMES, NDSTEQ,      &
                  NSV_DSTBEG, NSV_DSTEND, NSLTEQ, NSV_SLTBEG, NSV_SLTEND,     &
                  CDSTNAMES, CSLTNAMES, CCHEM_SURF_FILE,                      &
                  XSFDST, XSFDSTM, XSFSLT,                                    &
                  XAOSIP, XAOSIM, XAOSJP, XAOSJM, XHO2IP, XHO2IM, XHO2JP,     &
                  XHO2JM, XZ0, XZ0EFFIP, XZ0EFFIM, XZ0EFFJP, XZ0EFFJM, XZ0REL,&
                  XCLAY, XSAND, CPEDOTF,                                      &
                  XCONDSAT, XMPOTSAT, XBCOEF, XWWILT, XWFC, XW33, XWSAT,      &
                  XTAUICE, XCGSAT, XC1SAT, XC2REF, XC3, XC4B, XACOEF, XPCOEF, &
                  XC4REF, XPCPS, XPLVTT, XPLSTT,                              &
                  CSCOND, CISBA, XHCAPSOIL, XCONDDRY, XCONDSLD, CCPSURF,      &
                  XDG, XDROOT, XDG2, XROOTFRAC, XRUNOFFD, XDZG, XDZDIF,       &
                  XSOILWGHT, NWG_LAYER, NLAYER_HORT, NLAYER_DUN, XD_ICE,      &
                  XKSAT_ICE, XALBNIR_DRY, XALBVIS_DRY, XALBUV_DRY,            &
                  XALBNIR_WET, XALBVIS_WET, XALBUV_WET, XBSLAI_NITRO,         &
                  XCE_NITRO, XCNA_NITRO, XCF_NITRO                            )  
!
!------------------------------------------------------------------------------- 
!
IF(CISBA=='DIF'.AND.CSOC=='SGH')THEN
  CALL ABOR1_SFX('SUBGRID Soil organic matter effect (CSOC) NOT YET IMPLEMENTED FOR GARDEN')
ENDIF
!
IF(CKSAT=='SGH' .AND. HINIT/='PRE')THEN 
  !
  ALLOCATE(ZF(KI,NPATCH))
  ZF (:,:) = XUNDEF
  !  
  !Soil organic carbon effect and/or Exponential decay for DIF option 
  IF(CISBA=='DIF') THEN
    ALLOCATE(ZWORK(KI))
    ZWORK(:) = XUNDEF
    ZF(:,:) = 4.0/MERGE(XDROOT(:,:),XDG2(:,:),XDROOT(:,:)>0.0) 
  ELSE
    WHERE (ZF(:,:)==XUNDEF) ZF(:,:) =  4.0/XDG(:,2,:)
  ENDIF
  ZF(:,:)=MIN(ZF(:,:),XF_DECAY)
  !
  IF(CISBA=='DIF') THEN
    !   
    DO JPATCH=1,NPATCH    
      IF (NSIZE_NATURE_P(JPATCH) == 0 ) CYCLE 
      ZWORK(:) = MERGE(XDROOT(:,JPATCH),XDG2(:,JPATCH),XDROOT(:,JPATCH)>0.0)      
      CALL EXP_DECAY_SOIL_DIF(ZF(:,JPATCH),XDG(:,:,JPATCH),NWG_LAYER(:,JPATCH),ZWORK(:),&
                              XCONDSAT(:,:,JPATCH))   
    ENDDO  
    DEALLOCATE(ZWORK)

  !Exponential decay for ISBA-FR option
  ELSE
!
    DO JPATCH=1,NPATCH
       IF (NSIZE_NATURE_P(JPATCH) == 0 ) CYCLE
       CALL EXP_DECAY_SOIL_FR(CISBA, ZF(:,JPATCH),XC1SAT(:,JPATCH),XC2REF(:,JPATCH),    &
                                XDG(:,:,JPATCH),XD_ICE(:,JPATCH),XC4REF(:,JPATCH),      &
                                XC3(:,:,JPATCH),XCONDSAT(:,:,JPATCH),XKSAT_ICE(:,JPATCH))  
    ENDDO                       
! 
  ENDIF
  !  
  DEALLOCATE(ZF)
  !
ENDIF
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
CALL INIT_IO_SURF_n(HPROGRAM,'TOWN  ','TEB   ','READ ')
!
!*      10.     Prognostic and semi-prognostic fields
!               -------------------------------------
!
!* allocation of urban green area variables
!
!
CALL READ_TEB_GARDEN_n(HPROGRAM)
!
!-------------------------------------------------------------------------------
!
CALL COMMON_PARTS2(NPATCH, KI, LCANOPY, CROUGH, TSNOW, &
                   CPHOTO, XLAIMIN, XH_TREE, XVEGTYPE_PATCH, XLAI, XZ0, XVEG, XEMIS, &
                   LTR_ML, XFAPARC, XFAPIRC, XLAI_EFFC, XMUS, &
                   XALBNIR_SOIL, XALBVIS_SOIL, XALBUV_SOIL, XALBNIR, XALBVIS, XALBUV, &
                   LSURF_DIAG_ALBEDO, XPSN, XPSNG, XPSNV, XPSNV_A, &
                   PDIR_ALB, PSCA_ALB, PEMIS, PTSRAD ) 
!
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
DEALLOCATE(ZWG1)
!
CALL AVG_ALBEDO_EMIS_GARDEN(CALBEDO,                                &
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
DEALLOCATE(ZTG1)
!-------------------------------------------------------------------------------
!
!*      14.     ISBA diagnostics initialization
!               -------------------------------
!
!
CALL ALLOC_DIAG_TEB_GARDEN(KI,NGROUND_LAYER,KSW)
CALL DIAG_TEB_GARDEN_INIT_n(HPROGRAM,KI,KSW)
!
!-------------------------------------------------------------------------------
!
!         End of IO
!
CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('INIT_TEB_GARDEN_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE INIT_TEB_GARDEN_n
