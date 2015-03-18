!#############################################################
SUBROUTINE COMPUTE_ISBA_PARAMETERS(HPROGRAM,HINIT,OLAND_USE,            &
                             KI,KSV,KSW,                                &
                             HSV,PCO2,PRHOA,                            &
                             PZENITH,PSW_BANDS,PDIR_ALB,PSCA_ALB,       &
                             PEMIS,PTSRAD,PTSURF,                       &
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
!!      Modified by B. Decharme  (09/2012): Bug in exponential profile calculation with DIF
!!      F. Bouttier    08/13 : apply random perturbation patterns for ensembles
!!      B. Vincendon   03/14 : bug correction for CISBA=3L and CKSAT=EXP (TOPD coupling)
!!      Modified by B. Decharme  (04/2013): Subsurface runoff if SGH (DIF option only)
!!                                          Delete CTOPREG (never used)
!!                                          Delete NWG_LAYER_TOT, NWG_SIZE
!!                                          water table / Surface coupling
!!      P. Samuelsson  02/14 : MEB
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SFX_OASIS,  ONLY : LCPL_LAND, LCPL_FLOOD, LCPL_GW, LCPL_CALVING
!
USE MODD_ISBA_n,   ONLY : CROUGH, CISBA, CPEDOTF, CPHOTO, CRUNOFF, CALBEDO,   &
                          CSCOND, CRESPSL, LTR_ML, NNBIOMASS, NNLITTER,       &
                          NNLITTLEVS, NNSOILCARB, XCLAY, XSAND, XSOC,         &
                          XWWILT, XWFC, XWSAT, XWD0, XRM_PATCH, LPERM, LCOVER,&
                          XCOVER, XVEG, XLAI, XRSMIN, XGAMMA, XRGL, XCV,      &
                          XDG, NWG_LAYER, XDROOT, XDG2, XDZG, XDZDIF,         &
                          XZ0, XZ0_O_Z0H, XABC, XPOI, XANMAX, XFZERO, XEPSO,  &
                          XGAMM, XQDGAMM, XQDGMES, XT1GMES, XT2GMES, XAMAX,   &
                          XQDAMAX, XT1AMAX, XT2AMAX, XAH, XBH, XTAU_WOOD,     &
                          XINCREASE, XTURNOVER, XALBNIR_VEG, XALBVIS_VEG,     &
                          XALBUV_VEG, XEMIS, XVEGTYPE, XGMES, XRE25, XBSLAI,  &
                          XLAIMIN, XGC,XDMAX, LSTRESS, XF2I,                  &
                          XSEFOLD, XH_TREE, XPATCH, NPATCH, XWRMAX_CF,        &
                          NR_NATURE_P, NSIZE_NATURE_P,                        &
                          XALBNIR_DRY, XALBVIS_DRY, XALBUV_DRY,               &
                          XALBNIR_WET, XALBVIS_WET, XALBUV_WET,               &
                          XALBNIR_SOIL, XALBVIS_SOIL, XALBUV_SOIL,            &
                          XWG, XTG, TSNOW, XALBNIR, XALBVIS, XALBUV,          &
                          XEMIS_NAT, XFAPARC, XFAPIRC, XLAI_EFFC, XMUS,       &
                          XAOSIP,XAOSIM,XAOSJP,XAOSJM,                        &
                          XHO2IP,XHO2IM,XHO2JP,XHO2JM,                        &
                          XZ0EFFIP,XZ0EFFIM,XZ0EFFJP,XZ0EFFJM, XZ0REL,        &
                          XVEGTYPE_PATCH,XROOTFRAC,XRUNOFFD,XSOILWGHT,        &
                          XCGSAT, XC1SAT, XC2REF, XC3, XC4B, XACOEF, XPCOEF,  &
                          XTAUICE, XBCOEF, XCONDSAT,                          &
                          XHCAPSOIL, XCONDDRY, XCONDSLD, XC4REF, XMPOTSAT,    &
                          XTDEEP, XGAMMAT, NGROUND_LAYER, XSOILGRID, TTIME,   &
                          XCE_NITRO, XCF_NITRO,                               &
                          XCNA_NITRO, XBSLAI_NITRO, CCPSURF, TSEED,           &
                          TREAP, XWATSUP, XIRRIG, XCGMAX,                     &
                          CKSAT, CRAIN, LSOCP, LSOC, XFRACSOC,                &
                          XTI_MIN, XTI_MAX, XTI_MEAN, XTI_STD, XTI_SKEW,      &
                          XTAB_FSAT, XTAB_WTOP, XTAB_QTOP, XD_ICE, XKSAT_ICE, &
                          XFSAT, XMUF, LFLOOD, XFFLOOD, XFFROZEN,             &
                          XPIFLOOD, XCPL_EFLOOD, XCPL_PFLOOD, XCPL_IFLOOD,    &
                          XCPL_DRAIN, XCPL_RUNOFF, XCPL_RECHARGE, LGLACIER,   &
                          LTEMP_ARP, NTEMPLAYER_ARP, XPSN, XPSNG, XPSNV,      &
                          XPSNV_A, XFF, XFFG, XFFV, XPCPS, XPLVTT, XPLSTT,    &
                          LCANOPY, LCANOPY_DRAG, XDIR_ALB_WITH_SNOW,          &
                          XSCA_ALB_WITH_SNOW, XALBF, XEMISF, XCPL_ICEFLUX,    &
                          NLAYER_HORT, NLAYER_DUN, XF_PARAM,  XKANISO, XTOPQS,&
                          LGW, LWTD, XFWTD, XWTD, LCPL_RRM, LAGRI_TO_GRASS,   &
                          LPERTSURF,XPERTVEG,XPERTLAI,XPERTCV,XPERTALB,       &
                          XPERTZ0,                                            &
                          LMEB_PATCH,                                         &
                          XGNDLITTER, XZF_TALLVEG , XRGLGV,                   &
                          XGAMMAGV, XRSMINGV, XROOTFRACGV,                    &
                          XWRMAX_CFGV, XLAIGV, XZ0LITTER, XH_VEG, XTV
!
#ifdef TOPD
USE MODD_DUMMY_EXP_PROFILE,ONLY : XC_DEPTH_RATIO
#endif
!
USE MODD_CH_ISBA_n, ONLY : CSV, CCH_NAMES, NBEQ, NSV_CHSBEG, NSV_CHSEND,         &
                           CCHEM_SURF_FILE, NDSTEQ, NSV_DSTBEG, NSV_DSTEND,      &
                           NSV_AERBEG, NSV_AEREND, NAEREQ, CDSTNAMES, CAER_NAMES,&
                           NSLTEQ, NSV_SLTBEG,  NSV_SLTEND, CSLTNAMES,           &
                           LCH_BIO_FLUX, CCH_DRY_DEP  
!
USE MODD_DEEPSOIL,       ONLY : LPHYSDOMC, LDEEPSOIL, XTDEEP_CLI, XGAMMAT_CLI
USE MODD_AGRI,           ONLY : LAGRIP, XTHRESHOLD
USE MODD_AGRI_n,         ONLY : NIRRINUM, XTHRESHOLDSPT, LIRRIDAY, LIRRIGATE
!
USE MODD_DIAG_ISBA_n,      ONLY : LPATCH_BUDGET
USE MODD_DIAG_MISC_ISBA_n, ONLY : LSURF_DIAG_ALBEDO
!
USE MODD_SGH_PAR,        ONLY : NDIMTAB, XICE_DEPH_MAX, XF_DECAY
!
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE
USE MODD_SURF_PAR,       ONLY : XUNDEF, NUNDEF
USE MODD_SNOW_PAR,       ONLY : XEMISSN
!
USE MODD_TOPD_PAR, ONLY : NUNIT
USE MODD_TOPODYN, ONLY : NNCAT, NMESHT
USE MODD_SURF_ATM_n, ONLY : NR_NATURE, NDIM_FULL
!
USE MODD_DST_n
USE MODD_SLT_n
!
USE MODI_GET_LUOUT
USE MODI_ABOR1_SFX
USE MODI_INIT_IO_SURF_n
USE MODI_ALLOCATE_PHYSIO
USE MODI_INIT_ISBA_MIXPAR
USE MODI_CONVERT_PATCH_ISBA
USE MODI_INIT_VEG_PGD_n
USE MODI_INIT_TOP
USE MODI_EXP_DECAY_SOIL_FR
USE MODI_CARBON_INIT
USE MODI_SOILTEMP_ARP_PAR
USE MODI_END_IO_SURF_n
!
USE MODI_READ_SURF
USE MODI_READ_ISBA_n
USE MODI_INIT_ISBA_LANDUSE
USE MODI_READ_ISBA_CANOPY_n
USE MODI_INIT_VEG_n
USE MODI_AVERAGED_ALBEDO_EMIS_ISBA
USE MODI_DIAG_ISBA_INIT_n
USE MODI_INIT_SURF_TOPD
USE MODI_ISBA_SOC_PARAMETERS
!
USE MODI_READ_AND_SEND_MPI
USE MODI_ISBA_TO_TOPD
USE MODI_OPEN_FILE
USE MODI_CLOSE_FILE
USE MODI_FIX_MEB_VEG
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
REAL,             DIMENSION(KSW), INTENT(IN)  :: PSW_BANDS ! middle wavelength of each band
REAL,             DIMENSION(KI,KSW),INTENT(OUT) :: PDIR_ALB  ! direct albedo for each band
REAL,             DIMENSION(KI,KSW),INTENT(OUT) :: PSCA_ALB  ! diffuse albedo for each band
REAL,             DIMENSION(KI),  INTENT(OUT) :: PEMIS     ! emissivity
REAL,             DIMENSION(KI),  INTENT(OUT) :: PTSRAD    ! radiative temperature
REAL,             DIMENSION(KI),  INTENT(OUT) :: PTSURF    ! surface effective temperature         (K)
!
 CHARACTER(LEN=2),                 INTENT(IN)  :: HTEST       ! must be equal to 'OK'
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL, DIMENSION(NDIM_FULL)   :: ZF_PARAM, ZC_DEPTH_RATIO
!
REAL, DIMENSION(KI)     :: ZTSRAD_NAT !radiative temperature
REAL, DIMENSION(KI)     :: ZTSURF_NAT !effective temperature
!
REAL, DIMENSION(KI,NPATCH)  :: ZWG1 ! work array for surface water content
REAL, DIMENSION(KI,NPATCH)  :: ZTG1 ! work array for surface temperature
!
REAL, DIMENSION(KI)         :: ZM, ZWORK
REAL, DIMENSION(KI,NPATCH)  :: ZF, ZPERTBUF
!
INTEGER :: IDIM_FULL, JL
INTEGER           :: JILU     ! loop increment
INTEGER           :: ILUOUT   ! unit of output listing file
INTEGER           :: IRESP   ! return code
INTEGER           :: IDECADE, IDECADE2  ! decade of simulation
INTEGER           :: JPATCH  ! loop counter on tiles
INTEGER           :: INFOMPI
INTEGER           :: ISIZE_LMEB_PATCH  ! Number of patches with MEB=true
!
LOGICAL                           :: LWORK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!               Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('COMPUTE_ISBA_PARAMETERS',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('COMPUTE_ISBA_PARAMETERS: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF
!
ZF   (:,:) = XUNDEF
ZM   (:)   = XUNDEF
ZWORK(:)   = XUNDEF
!
!*       2.3    Physiographic data fields from land cover:
!               -----------------------------------------
!
 CALL ALLOCATE_PHYSIO(CPHOTO, CISBA, KI, NVEGTYPE, NGROUND_LAYER, NPATCH, &
                     XVEGTYPE, XLAI, XVEG, XZ0, XEMIS, XDG, XD_ICE,      &
                     XRSMIN, XGAMMA, XWRMAX_CF, XRGL, XCV,               &
                     XZ0_O_Z0H, XALBNIR_VEG, XALBVIS_VEG, XALBUV_VEG,    &
                     XH_TREE, XRE25, XLAIMIN, XBSLAI, XSEFOLD,           &
                     XGMES, XGC, XF2I, XDMAX, LSTRESS,                   &
                     XCE_NITRO, XCF_NITRO, XCNA_NITRO,                   &
                     TSEED, TREAP, XWATSUP, XIRRIG,                      &
                     XROOTFRAC, NWG_LAYER, XDROOT, XDG2,                 &
                     XGNDLITTER,XZF_TALLVEG,XRGLGV,XGAMMAGV,XRSMINGV,    &
                     XROOTFRACGV,XWRMAX_CFGV,XLAIGV,XZ0LITTER,XH_VEG     )
!
IF (TTIME%TDATE%MONTH /= NUNDEF) THEN
  IDECADE = 3 * ( TTIME%TDATE%MONTH - 1 ) + MIN(TTIME%TDATE%DAY-1,29) / 10 + 1
ELSE
  IDECADE = 1
END IF
!
IDECADE2 = IDECADE
!
 CALL INIT_ISBA_MIXPAR(CISBA,IDECADE,IDECADE2,XCOVER,LCOVER,CPHOTO,'NAT')
!
ISIZE_LMEB_PATCH=COUNT(LMEB_PATCH(:))
IF (ISIZE_LMEB_PATCH>0)  THEN
  CALL FIX_MEB_VEG(NPATCH)
ENDIF
!
 CALL CONVERT_PATCH_ISBA(CISBA,IDECADE,IDECADE2,XCOVER,LCOVER,CPHOTO,LAGRIP,   &
                        LPERM,LTR_ML,'NAT',PVEG=XVEG,PLAI=XLAI,                &
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
                        TPREAP=TREAP,PWATSUP=XWATSUP,PIRRIG=XIRRIG,            &
                        PGNDLITTER=XGNDLITTER,PZF_TALLVEG=XZF_TALLVEG,         &
                        PRGLGV=XRGLGV,   &
                        PGAMMAGV=XGAMMAGV,PRSMINGV=XRSMINGV,                   &
                        PROOTFRACGV=XROOTFRACGV,PWRMAX_CFGV=XWRMAX_CFGV,       &
                        PLAIGV=XLAIGV,PZ0LITTER=XZ0LITTER,PH_VEG=XH_VEG        )
!
!-------------------------------------------------------------------------------
!
CALL INIT_VEG_PGD_n(HPROGRAM, 'NATURE', ILUOUT, KI, NPATCH, NGROUND_LAYER,      &
                    TTIME%TDATE%MONTH,                                          &
                    XVEGTYPE, XPATCH, XVEGTYPE_PATCH, NSIZE_NATURE_P,           &
                    NR_NATURE_P, XRM_PATCH,                                     &
                    LDEEPSOIL, LPHYSDOMC, XTDEEP_CLI, XGAMMAT_CLI, XTDEEP,      &
                    XGAMMAT, LAGRIP, XTHRESHOLD, NIRRINUM, LIRRIDAY, LIRRIGATE, &
                    XTHRESHOLDSPT,                                              &
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
                    XCONDSAT, XMPOTSAT, XBCOEF, XWWILT, XWFC, XWSAT, XWD0,      &
                    XKANISO, CRUNOFF,                                           &
                    XTAUICE, XCGSAT, XC1SAT, XC2REF, XC3, XC4B, XACOEF, XPCOEF, &
                    XC4REF, XPCPS, XPLVTT, XPLSTT,                              &
                    CSCOND, CISBA, XHCAPSOIL, XCONDDRY, XCONDSLD, CCPSURF,      &
                    XDG, XDROOT, XDG2, XROOTFRAC, XRUNOFFD, XDZG, XDZDIF,       &
                    XSOILWGHT, NWG_LAYER, NLAYER_HORT, NLAYER_DUN, XD_ICE,      &
                    XKSAT_ICE, XALBNIR_DRY, XALBVIS_DRY, XALBUV_DRY,            &
                    XALBNIR_WET, XALBVIS_WET, XALBUV_WET, XBSLAI_NITRO,         &
                    XCE_NITRO, XCNA_NITRO, XCF_NITRO, XFWTD, XWTD               )  
!
!-------------------------------------------------------------------------------
!
!DIF option :
!    Anisotropy coeficient for hydraulic conductivity for topmodel drainage (Fan et al. 2006)
!    Soil organic matter effect and/or Exponential decay for DIF option
!    Must be call before INIT_TOP
!
!
IF(CISBA=='DIF') THEN
  !
  IF( CKSAT=='SGH' )THEN
    WRITE(ILUOUT,*)'THE KSAT EXP PROFILE WITH ISBA-DF IS NOT PHYSIC AND HAS BEEN REMOVED FOR NOW' 
    WRITE(ILUOUT,*)'A NEW PHYSICAL APPROACH WILL BE DEVELLOPED ACCOUNTING FOR COMPACTION IN ALL ' 
    WRITE(ILUOUT,*)'HYDRODYNAMIC PARAMETERS (WSAT, PSISAT, KSAT, B) AND NOT ONLY IN KSAT        ' 
    CALL ABOR1_SFX('CKSAT=SGH is not physic with ISBA-DF and has been removed for now')    
  ENDIF
  !  
  IF(LSOC)THEN   
    IF(.NOT.LSOCP)THEN
      WRITE(ILUOUT,*)'LSOC = T can be activated only if SOC data given in PGD fields'
      CALL ABOR1_SFX('LSOC = T can be activated only if SOC data given in PGD fields')
    ENDIF
    ALLOCATE(XFRACSOC(KI,NGROUND_LAYER))
    XFRACSOC(:,:)=0.0
    CALL ISBA_SOC_PARAMETERS(CRUNOFF,XPATCH,XDG,XSOC,XBCOEF,XMPOTSAT,   &
                             XCONDSAT,XWSAT,XHCAPSOIL,XCONDDRY,         &
                             XCONDSLD,XWFC,XWWILT,XWD0,XKANISO,XFRACSOC )
  ELSE
    ALLOCATE(XFRACSOC(0,0))
  ENDIF
!
ELSE
  ALLOCATE(XFRACSOC(0,0))
ENDIF
!
!Topmodel
!
!CRUNOFF used in hydro_sgh and isba_sgh_update
IF( CRUNOFF=='SGH ') THEN 
!
  ALLOCATE(XTAB_FSAT(KI,NDIMTAB))
  ALLOCATE(XTAB_WTOP(KI,NDIMTAB))
  ALLOCATE(XTAB_QTOP(KI,NDIMTAB))
!
  XTAB_FSAT(:,:) = 0.0
  XTAB_WTOP(:,:) = 0.0
  XTAB_QTOP(:,:) = 0.0
!
  IF(HINIT/='PRE')THEN
!
    WHERE(XCLAY(:,1)==XUNDEF.AND.XTI_MEAN(:)/=XUNDEF) XTI_MEAN(:)=XUNDEF
!
    CALL INIT_TOP (CISBA, ILUOUT, XPATCH, XRUNOFFD,          &
                   XWD0, XWSAT, XTI_MIN,                     &
                   XTI_MAX, XTI_MEAN, XTI_STD, XTI_SKEW,     &
                   XSOILWGHT, XTAB_FSAT, XTAB_WTOP,          &
                   XTAB_QTOP, ZM                             )  
!
!
    IF (CKSAT=='SGH' .AND. CISBA/='DIF') THEN
!     Exponential decay factor calculate using soil properties 
!     (eq. 11, Decharme et al., J. Hydrometeor, 2006)
      DO JILU=1,KI
        IF (ZM(JILU)/=XUNDEF) ZF(JILU,:) = (XWSAT(JILU,1)-XWD0(JILU,1))/ZM(JILU)
      ENDDO
!       
    ENDIF
!
  ENDIF
!
! Subsurface flow by layer (m/s)
  IF(CISBA=='DIF') THEN
    ALLOCATE(XTOPQS(KI,NGROUND_LAYER,NPATCH))
    XTOPQS(:,:,:)=0.0
  ELSE
    ALLOCATE(XTOPQS(0,0,0))
  ENDIF
!
ELSE                  
!  
  ALLOCATE(XTAB_FSAT(0,0))
  ALLOCATE(XTAB_WTOP(0,0))
  ALLOCATE(XTAB_QTOP(0,0))
  ALLOCATE(XTOPQS(0,0,0))  
!                  
ENDIF  
!
!Exponential decay for ISBA-FR option
!CKSAT used in hydro_soil.F90 and soil.F90
IF(HINIT/='PRE'.AND.CISBA/='DIF')THEN 
  !
  IF(CKSAT=='SGH') THEN
    !
    WHERE(ZF(:,:)==XUNDEF.AND.XDG(:,2,:)/=XUNDEF) 
      ZF(:,:) = 4.0/XDG(:,2,:)
    ENDWHERE
    ZF(:,:) = MIN(ZF(:,:),XF_DECAY)
    !
    ALLOCATE(XF_PARAM (KI))
    XF_PARAM(:) = ZF(:,1)
    !
    DO JPATCH=1,NPATCH
      IF (NSIZE_NATURE_P(JPATCH) == 0 ) CYCLE
        CALL EXP_DECAY_SOIL_FR(CISBA, ZF(:,JPATCH),XC1SAT(:,JPATCH),XC2REF(:,JPATCH),   &
                                XDG(:,:,JPATCH),XD_ICE(:,JPATCH),XC4REF(:,JPATCH),      &
                                XC3(:,:,JPATCH),XCONDSAT(:,:,JPATCH),XKSAT_ICE(:,JPATCH))  
    ENDDO                       
    !
  ELSEIF ( CKSAT=='EXP' .AND. CISBA=='3-L' ) THEN
    !
    ALLOCATE(XF_PARAM (KI))
    XF_PARAM(:) = XUNDEF
    !
    IF (HPROGRAM/='AROME ' .AND. HPROGRAM/='MESONH ') THEN
      !
      CALL OPEN_FILE('ASCII ',NUNIT,HFILE='carte_f_dc.txt',HFORM='FORMATTED',HACTION='READ ')
      DO JILU=1,NDIM_FULL
        READ(NUNIT,*) ZF_PARAM(JILU), ZC_DEPTH_RATIO(JILU)
      ENDDO
      CALL CLOSE_FILE('ASCII ',NUNIT)
      CALL READ_AND_SEND_MPI(ZF_PARAM,XF_PARAM,NR_NATURE)
#ifdef TOPD
IF (.NOT.ALLOCATED(XC_DEPTH_RATIO))    ALLOCATE(XC_DEPTH_RATIO (KI))
    XC_DEPTH_RATIO(:) = XUNDEF
      CALL READ_AND_SEND_MPI(ZC_DEPTH_RATIO,XC_DEPTH_RATIO,NR_NATURE)
#endif
      !
    ELSE
      WRITE(ILUOUT,*) "COMPUTE_ISBA_PARAMETERS: WITH CKSAT=EXP, IN NOT OFFLINE "//&
                      "MODE, TOPMODEL FILE FOR F_PARAM IS NOT READ "
    ENDIF
    !
    DO JPATCH=1,NPATCH
      WHERE (XF_PARAM(:)==XUNDEF.AND.XDG(:,2,JPATCH)/=XUNDEF)
        ZF(:,JPATCH) = 4.0/XDG(:,2,JPATCH)
      ELSEWHERE
        ZF(:,JPATCH) = XF_PARAM(:)
      ENDWHERE
    ENDDO
     ZF(:,:) = MIN(ZF(:,:),XF_DECAY)
    !
    DO JPATCH=1,NPATCH
      CALL EXP_DECAY_SOIL_FR(CISBA, ZF(:,JPATCH),XC1SAT(:,JPATCH),XC2REF(:,JPATCH), &
                             XDG(:,:,JPATCH),XD_ICE(:,JPATCH),XC4REF(:,JPATCH),   &
                             XC3(:,:,JPATCH),XCONDSAT(:,:,JPATCH),                &
                             XKSAT_ICE(:,JPATCH))  
    ENDDO    
    !    
  ENDIF
  ! 
ENDIF
!
!
!*       2.10   Soil carbon
!               -----------                        
!
IF (HINIT == 'ALL' .AND. CRESPSL=='CNT' .AND. CPHOTO == 'NCB') THEN
  CALL CARBON_INIT(NNBIOMASS, NNLITTER, NNLITTLEVS, NNSOILCARB)
ENDIF
!
!Rainfall spatial distribution
!CRAIN used in HYDRO_VEG and HYDRO_SGH and ISBA_SGH_UPDATE
IF(CRAIN=='SGH')THEN
  ALLOCATE(XMUF(KI))
  XMUF(:)=0.0
ELSE
  ALLOCATE(XMUF(0))
ENDIF
!
ALLOCATE(XFSAT(KI))  
XFSAT(:) = 0.0
!
!-------------------------------------------------------------------------------
!
!*       6.2    Initialize of SFX - RRM coupling:
!               ---------------------------------
!
! * Check some key :
!
IF(LCPL_CALVING)THEN
   IF(.NOT.LGLACIER)THEN
     CALL ABOR1_SFX('COMPUTE_ISBA_PARAMETERS: LGLACIER MUST BE ACTIVATED IF LCPL_CALVING')
   ENDIF
ENDIF
!
! * Initialize required coupling fields :
!
LCPL_RRM = .FALSE.
LFLOOD   = .FALSE.
LWTD     = .FALSE.
!
IF(LCPL_LAND)THEN
!    
  LCPL_RRM = .TRUE.
!
  ALLOCATE(XCPL_DRAIN (KI))
  ALLOCATE(XCPL_RUNOFF(KI))
  XCPL_DRAIN (:) = 0.0
  XCPL_RUNOFF(:) = 0.0
!
  IF(LGLACIER)THEN
     ALLOCATE(XCPL_ICEFLUX(KI))
     XCPL_ICEFLUX(:) = 0.0
  ELSE
     ALLOCATE(XCPL_ICEFLUX(0))
  ENDIF
!
  IF(LCPL_GW)THEN
    LWTD = .TRUE.
    ALLOCATE(XCPL_RECHARGE(KI))
    XCPL_RECHARGE(:) = 0.0
  ELSE
    ALLOCATE(XCPL_RECHARGE(0))
  ENDIF
!
  IF(LCPL_FLOOD)THEN
     LFLOOD = .TRUE.
     ALLOCATE(XCPL_EFLOOD(KI))
     ALLOCATE(XCPL_PFLOOD(KI))
     ALLOCATE(XCPL_IFLOOD(KI))
     XCPL_EFLOOD(:)= 0.0
     XCPL_PFLOOD(:)= 0.0
     XCPL_IFLOOD(:)= 0.0    
  ELSE
    ALLOCATE(XCPL_EFLOOD(0))
    ALLOCATE(XCPL_PFLOOD(0))
    ALLOCATE(XCPL_IFLOOD(0))     
  ENDIF     
!
ELSE
!
  ALLOCATE(XCPL_RUNOFF  (0))
  ALLOCATE(XCPL_DRAIN   (0))
  ALLOCATE(XCPL_ICEFLUX (0))
  ALLOCATE(XCPL_RECHARGE(0))
  ALLOCATE(XCPL_EFLOOD  (0))
  ALLOCATE(XCPL_PFLOOD  (0))
  ALLOCATE(XCPL_IFLOOD  (0))
!
ENDIF
!
IF(LWTD.AND..NOT.LGW)THEN
  WRITE(ILUOUT,*)'COMPUTE_ISBA_PARAMETERS: Groundwater map is required by SFX - Groundwater coupling '
  WRITE(ILUOUT,*)'COMPUTE_ISBA_PARAMETERS: Please check your pgd namelist where this map must be     '
  WRITE(ILUOUT,*)'COMPUTE_ISBA_PARAMETERS: specified (YGW and YGWFILETYPE, or XUNIF_GW, or LIMP_GW)  '
  CALL ABOR1_SFX('COMPUTE_ISBA_PARAMETERS: Groundwater map is required by SFX - Groundwater coupling')
ENDIF
!
! * Initialize flood scheme :
!
IF(LFLOOD)THEN
  ALLOCATE(XFFLOOD (KI))
  ALLOCATE(XPIFLOOD(KI))
  ALLOCATE(XFF     (KI,NPATCH))
  ALLOCATE(XFFG    (KI,NPATCH))
  ALLOCATE(XFFV    (KI,NPATCH))
  ALLOCATE(XFFROZEN(KI,NPATCH))
  ALLOCATE(XALBF   (KI,NPATCH))
  ALLOCATE(XEMISF  (KI,NPATCH)) 
  XFFLOOD       = 0.0
  XPIFLOOD      = 0.0
  XFF           = 0.0
  XFFG          = 0.0
  XFFV          = 0.0
  XFFROZEN      = 0.0
  XALBF         = 0.0
  XEMISF        = 0.0  
ELSE
  ALLOCATE(XFFLOOD   (0))
  ALLOCATE(XPIFLOOD  (0))
  ALLOCATE(XFF     (0,0))
  ALLOCATE(XFFG    (0,0))
  ALLOCATE(XFFV    (0,0))
  ALLOCATE(XFFROZEN(0,0))
  ALLOCATE(XALBF   (0,0))
  ALLOCATE(XEMISF  (0,0))
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
!*       9.     Prints of cover parameters in a tex file
!               ----------------------------------------
!
!* if only physiographic fields are to be initialized, stop here.
!
IF (HINIT/='ALL' .AND. HINIT/='SOD') THEN
  IF (LHOOK) CALL DR_HOOK('COMPUTE_ISBA_PARAMETERS',1,ZHOOK_HANDLE)
  RETURN
END IF
!
!-------------------------------------------------------------------------------
!
CALL INIT_IO_SURF_n(HPROGRAM,'NATURE','ISBA  ','READ ')
!
!*      10.     Prognostic and semi-prognostic fields
!               -------------------------------------
!
 CALL READ_ISBA_n(HPROGRAM)
!
IF (HINIT/='ALL') THEN
  CALL END_IO_SURF_n(HPROGRAM)
  IF (LHOOK) CALL DR_HOOK('COMPUTE_ISBA_PARAMETERS',1,ZHOOK_HANDLE)
  RETURN
END IF
!
IF (HINIT=='PRE' .AND. TSNOW%SCHEME.NE.'3-L' .AND. TSNOW%SCHEME.NE.'CRO' .AND. CISBA=='DIF') THEN
    CALL ABOR1_SFX("INIT_ISBAN: WITH CISBA = DIF, CSNOW MUST BE 3-L OR CRO")
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      11.  Extrapolation of the prognostic and semi-prognostic fields
!                           LAND USE case 
!               -------------------------------------
!
IF (OLAND_USE) THEN
   CALL INIT_ISBA_LANDUSE(HPROGRAM)  
END IF
!
!-------------------------------------------------------------------------------
!
!*      12.     Canopy air fields:
!               -----------------
!
 CALL READ_ISBA_CANOPY_n(HPROGRAM)
!
!-------------------------------------------------------------------------------
!
!*      13.     initialize radiative and physical properties
!               --------------------------------------------
!
ALLOCATE(XDIR_ALB_WITH_SNOW(KI,KSW,NPATCH))
ALLOCATE(XSCA_ALB_WITH_SNOW(KI,KSW,NPATCH))
XDIR_ALB_WITH_SNOW = 0.0
XSCA_ALB_WITH_SNOW = 0.0
!
 CALL INIT_VEG_n(NPATCH, KI, LCANOPY, CROUGH, LAGRI_TO_GRASS, TSNOW, &
                 CPHOTO, XLAIMIN, XH_TREE, XVEGTYPE_PATCH, XLAI, XZ0, XVEG, XEMIS, &
                 LTR_ML, XFAPARC, XFAPIRC, XLAI_EFFC, XMUS, &
                 XALBNIR_SOIL, XALBVIS_SOIL, XALBUV_SOIL, XALBNIR, XALBVIS, XALBUV, &
                 LSURF_DIAG_ALBEDO, XPSN, XPSNG, XPSNV, XPSNV_A, &
                 PDIR_ALB, PSCA_ALB, PEMIS, PTSRAD )
!
DO JPATCH=1,NPATCH
  ZWG1(:,JPATCH) = XWG(:,1,JPATCH)
  ZTG1(:,JPATCH) = XTG(:,1,JPATCH)
END DO
!
 CALL CONVERT_PATCH_ISBA(CISBA,IDECADE,IDECADE2,XCOVER,LCOVER,&
                          CPHOTO,LAGRIP,LPERM,LTR_ML,'NAT',   &
                          PWG1 = ZWG1,               &
                          PALBNIR_SOIL=XALBNIR_SOIL, &
                          PALBVIS_SOIL=XALBVIS_SOIL, &
                          PALBUV_SOIL=XALBUV_SOIL )
!
! Load randomly perturbed fields. Perturbation ratios are saved in case fields are reset later.
IF(LPERTSURF) THEN
!
  CALL READ_SURF(HPROGRAM,'VEG',XVEG(:,:),IRESP)
  ALLOCATE(XPERTVEG(KI))
  XPERTVEG(:)=XVEG(:,1)
!
  CALL READ_SURF(HPROGRAM,'LAI',XLAI(:,:),IRESP)
  ALLOCATE(XPERTLAI(KI))
  XPERTLAI(:)=XLAI(:,1)
!
  CALL READ_SURF(HPROGRAM,'CV',XCV(:,:),IRESP)
  ALLOCATE(XPERTCV(KI))
  XPERTCV(:)=XCV(:,1)
!
  CALL READ_SURF(HPROGRAM,'PERTALB',ZPERTBUF(:,:),IRESP)
  ALLOCATE(XPERTALB(KI))
  XPERTALB(:)=ZPERTBUF(:,1)
  WHERE(XALBNIR_VEG(:,1)/=XUNDEF)  XALBNIR_VEG(:,1) = XALBNIR_VEG(:,1) *( 1.+ XPERTALB(:) )
  WHERE(XALBVIS_VEG(:,1)/=XUNDEF)  XALBVIS_VEG(:,1) = XALBVIS_VEG(:,1) *( 1.+ XPERTALB(:) )
  WHERE(XALBUV_VEG(:,1)/=XUNDEF)   XALBUV_VEG(:,1)  = XALBUV_VEG(:,1)  *( 1.+ XPERTALB(:) )
  WHERE(XALBNIR_SOIL(:,1)/=XUNDEF) XALBNIR_SOIL(:,1)= XALBNIR_SOIL(:,1)*( 1.+ XPERTALB(:) )
  WHERE(XALBVIS_SOIL(:,1)/=XUNDEF) XALBVIS_SOIL(:,1)= XALBVIS_SOIL(:,1)*( 1.+ XPERTALB(:) )
  WHERE(XALBUV_SOIL(:,1)/=XUNDEF)  XALBUV_SOIL(:,1) = XALBUV_SOIL(:,1) *( 1.+ XPERTALB(:) )
!
  CALL READ_SURF(HPROGRAM,'PERTZ0LAND',ZPERTBUF(:,:),IRESP)
  ALLOCATE(XPERTZ0(KI))
  XPERTZ0(:)=ZPERTBUF(:,1)
  WHERE(XZ0(:,1)/=XUNDEF)      XZ0(:,1)     =XZ0(:,1)     *( 1.+ XPERTZ0(:) )
  WHERE(XZ0EFFIP(:,1)/=XUNDEF) XZ0EFFIP(:,1)=XZ0EFFIP(:,1)*( 1.+ XPERTZ0(:) )
  WHERE(XZ0EFFIM(:,1)/=XUNDEF) XZ0EFFIM(:,1)=XZ0EFFIM(:,1)*( 1.+ XPERTZ0(:) )
  WHERE(XZ0EFFJP(:,1)/=XUNDEF) XZ0EFFJP(:,1)=XZ0EFFJP(:,1)*( 1.+ XPERTZ0(:) )
  WHERE(XZ0EFFJM(:,1)/=XUNDEF) XZ0EFFJM(:,1)=XZ0EFFJM(:,1)*( 1.+ XPERTZ0(:) )
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       14.    Output radiative fields
!               -----------------------
!
ALLOCATE(XEMIS_NAT   (KI))
XEMIS_NAT (:) = XUNDEF
!
 CALL AVERAGED_ALBEDO_EMIS_ISBA(LFLOOD, CALBEDO, PZENITH,                &
                                 XVEG,XZ0,XLAI,                          &
                                 LMEB_PATCH,XGNDLITTER,XZ0LITTER,XLAIGV, &
                                 XZF_TALLVEG, XH_VEG, XTV,               &
                                 ZTG1,                                   &
                                 XPATCH,                                 &
                                 PSW_BANDS,                              &
                                 XALBNIR_VEG,XALBVIS_VEG,XALBUV_VEG,     &
                                 XALBNIR_SOIL,XALBVIS_SOIL,XALBUV_SOIL,  &
                                 XEMIS,                                  &
                                 TSNOW,                                  &
                                 XALBNIR,XALBVIS,XALBUV,                 &
                                 PDIR_ALB, PSCA_ALB,                     &
                                 XEMIS_NAT,ZTSRAD_NAT,ZTSURF_NAT         )  
!
PEMIS  = XEMIS_NAT
PTSRAD = ZTSRAD_NAT
PTSURF = ZTSURF_NAT
!
!-------------------------------------------------------------------------------
!
!*      15.     ISBA diagnostics initialization
!               -------------------------------
!
IF(NPATCH<=1) LPATCH_BUDGET=.FALSE.
!
 CALL DIAG_ISBA_INIT_n(HPROGRAM,KI,KSW)
!
!-------------------------------------------------------------------------------
!
 CALL INIT_SURF_TOPD(HPROGRAM,NDIM_FULL)
!
!-------------------------------------------------------------------------------
!
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('COMPUTE_ISBA_PARAMETERS',1,ZHOOK_HANDLE)
!
END SUBROUTINE COMPUTE_ISBA_PARAMETERS


