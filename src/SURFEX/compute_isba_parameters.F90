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
USE MODD_ISBA_CANOPY_n, ONLY : ICP => ISBA_CANOPY
!
USE MODD_SURF_ATM_GRID_n, ONLY : UG => SURF_ATM_GRID
!
USE MODD_IO_BUFF_n, ONLY : IOB => IO_BUFF
!
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
!
USE MODD_DATA_ISBA_n, ONLY : DTI => DATA_ISBA
USE MODD_ISBA_GRID_n, ONLY : IG => ISBA_GRID
!
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DGEI => DIAG_EVAP_ISBA
USE MODD_DIAG_SURF_ATM_n, ONLY : DGU => DIAG_SURF_ATM
USE MODD_GR_BIOG_n, ONLY : GB => GR_BIOG
!
USE MODD_SFX_OASIS,  ONLY : LCPL_LAND, LCPL_FLOOD, LCPL_GW, LCPL_CALVING
!
USE MODD_ISBA_n, ONLY : I => ISBA
!
#ifdef TOPD
USE MODD_DUMMY_EXP_PROFILE,ONLY : XC_DEPTH_RATIO
#endif
!
USE MODD_CH_ISBA_n, ONLY : CHI => CH_ISBA
!
USE MODD_DEEPSOIL,       ONLY : LPHYSDOMC, LDEEPSOIL, XTDEEP_CLI, XGAMMAT_CLI
USE MODD_AGRI,           ONLY : LAGRIP, XTHRESHOLD
USE MODD_AGRI_n, ONLY : AG => AGRI
!
USE MODD_DIAG_ISBA_n, ONLY : DGI => DIAG_ISBA
USE MODD_DIAG_MISC_ISBA_n, ONLY : DGMI => DIAG_MISC_ISBA
!
USE MODD_SGH_PAR,        ONLY : NDIMTAB, XICE_DEPH_MAX, XF_DECAY
!
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE
USE MODD_SURF_PAR,       ONLY : XUNDEF, NUNDEF
USE MODD_SNOW_PAR,       ONLY : XEMISSN
!
USE MODD_TOPD_PAR, ONLY : NUNIT
USE MODD_TOPODYN, ONLY : NNCAT, NMESHT
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
USE MODD_DST_n, ONLY : DST => DST
USE MODD_SLT_n, ONLY : SLT => SLT
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
REAL, DIMENSION(U%NDIM_FULL)   :: ZF_PARAM, ZC_DEPTH_RATIO
!
REAL, DIMENSION(KI)     :: ZTSRAD_NAT !radiative temperature
REAL, DIMENSION(KI)     :: ZTSURF_NAT !effective temperature
!
REAL, DIMENSION(KI,I%NPATCH)  :: ZWG1 ! work array for surface water content
REAL, DIMENSION(KI,I%NPATCH)  :: ZTG1 ! work array for surface temperature
!
REAL, DIMENSION(KI)         :: ZM, ZWORK
REAL, DIMENSION(KI,I%NPATCH)  :: ZF, ZPERTBUF
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
 CALL ALLOCATE_PHYSIO(I, &
                      I%CPHOTO, I%CISBA, KI, NVEGTYPE, I%NGROUND_LAYER, I%NPATCH, &
                     I%XVEGTYPE, I%XLAI, I%XVEG, I%XZ0, I%XEMIS, I%XDG, I%XD_ICE,      &
                     I%XRSMIN, I%XGAMMA, I%XWRMAX_CF, I%XRGL, I%XCV,               &
                     I%XZ0_O_Z0H, I%XALBNIR_VEG, I%XALBVIS_VEG, I%XALBUV_VEG,    &
                     I%XH_TREE, I%XRE25, I%XLAIMIN, I%XBSLAI, I%XSEFOLD,           &
                     I%XGMES, I%XGC, I%XF2I, I%XDMAX, I%LSTRESS,                   &
                     I%XCE_NITRO, I%XCF_NITRO, I%XCNA_NITRO,                   &
                     I%TSEED, I%TREAP, I%XWATSUP, I%XIRRIG,                      &
                     I%XROOTFRAC, I%NWG_LAYER, I%XDROOT, I%XDG2,                 &
                     I%XGNDLITTER,I%XZF_TALLVEG,I%XRGLGV,I%XGAMMAGV,I%XRSMINGV,    &
                     I%XROOTFRACGV,I%XWRMAX_CFGV,I%XLAIGV,I%XZ0LITTER,I%XH_VEG     )
!
IF (I%TTIME%TDATE%MONTH /= NUNDEF) THEN
  IDECADE = 3 * ( I%TTIME%TDATE%MONTH - 1 ) + MIN(I%TTIME%TDATE%DAY-1,29) / 10 + 1
ELSE
  IDECADE = 1
END IF
!
IDECADE2 = IDECADE
!
 CALL INIT_ISBA_MIXPAR(DTCO, DTI, IG, I, &
                       I%CISBA,IDECADE,IDECADE2,I%XCOVER,I%LCOVER,I%CPHOTO,'NAT')
!
ISIZE_LMEB_PATCH=COUNT(I%LMEB_PATCH(:))
IF (ISIZE_LMEB_PATCH>0)  THEN
  CALL FIX_MEB_VEG(DTI, IG, I, &
                   I%NPATCH)
ENDIF
!
 CALL CONVERT_PATCH_ISBA(DTCO, DTI, I, &
                         I%CISBA,IDECADE,IDECADE2,I%XCOVER,I%LCOVER,I%CPHOTO,LAGRIP,   &
                        I%LPERM,I%LTR_ML,'NAT',PVEG=I%XVEG,PLAI=I%XLAI,                &
                        PRSMIN=I%XRSMIN,PGAMMA=I%XGAMMA,PWRMAX_CF=I%XWRMAX_CF,       &
                        PRGL=I%XRGL,PCV=I%XCV,PSOILGRID=I%XSOILGRID,                 &
                        PDG=I%XDG,KWG_LAYER=I%NWG_LAYER,PDROOT=I%XDROOT,PDG2=I%XDG2,   &
                        PZ0=I%XZ0,PZ0_O_Z0H=I%XZ0_O_Z0H,                           &
                        PALBNIR_VEG=I%XALBNIR_VEG,PALBVIS_VEG=I%XALBVIS_VEG,       &
                        PALBUV_VEG=I%XALBUV_VEG,PEMIS_ECO=I%XEMIS,                 &
                        PVEGTYPE=I%XVEGTYPE,PROOTFRAC=I%XROOTFRAC,                 &
                        PGMES=I%XGMES,PBSLAI=I%XBSLAI,PLAIMIN=I%XLAIMIN,             &
                        PSEFOLD=I%XSEFOLD,PGC=I%XGC,                               &
                        PDMAX=I%XDMAX,PF2I=I%XF2I,OSTRESS=I%LSTRESS,PH_TREE=I%XH_TREE, &
                        PRE25=I%XRE25,PCE_NITRO=I%XCE_NITRO,PCF_NITRO=I%XCF_NITRO,   &
                        PCNA_NITRO=I%XCNA_NITRO,PD_ICE=I%XD_ICE,TPSEED=I%TSEED,      &
                        TPREAP=I%TREAP,PWATSUP=I%XWATSUP,PIRRIG=I%XIRRIG,            &
                        PGNDLITTER=I%XGNDLITTER,PZF_TALLVEG=I%XZF_TALLVEG,         &
                        PRGLGV=I%XRGLGV,   &
                        PGAMMAGV=I%XGAMMAGV,PRSMINGV=I%XRSMINGV,                   &
                        PROOTFRACGV=I%XROOTFRACGV,PWRMAX_CFGV=I%XWRMAX_CFGV,       &
                        PLAIGV=I%XLAIGV,PZ0LITTER=I%XZ0LITTER,PH_VEG=I%XH_VEG        )
!
!-------------------------------------------------------------------------------
!
CALL INIT_VEG_PGD_n(CHI, DTCO, DST, I, SLT, U, &
                    HPROGRAM, 'NATURE', ILUOUT, KI, I%NPATCH, I%NGROUND_LAYER,      &
                    I%TTIME%TDATE%MONTH,                                          &
                    I%XVEGTYPE, I%XPATCH, I%XVEGTYPE_PATCH, I%NSIZE_NATURE_P,           &
                    I%NR_NATURE_P, I%XRM_PATCH,                                     &
                    LDEEPSOIL, LPHYSDOMC, XTDEEP_CLI, XGAMMAT_CLI, I%XTDEEP,      &
                    I%XGAMMAT, LAGRIP, XTHRESHOLD, AG%NIRRINUM, AG%LIRRIDAY, AG%LIRRIGATE, &
                    AG%XTHRESHOLDSPT,                                              &
                    I%CPHOTO, HINIT, I%LTR_ML, I%NNBIOMASS, PCO2, PRHOA, I%XABC, I%XPOI,  &
                    I%XGMES, I%XGC, I%XDMAX, I%XANMAX, I%XFZERO, I%XEPSO, I%XGAMM, I%XQDGAMM,   & 
                    I%XQDGMES, I%XT1GMES, I%XT2GMES, I%XAMAX, I%XQDAMAX, I%XT1AMAX, I%XT2AMAX,&
                    I%XAH, I%XBH, I%XTAU_WOOD, I%XINCREASE, I%XTURNOVER,                  &
                    KSV, HSV, CHI%NBEQ, CHI%CSV, CHI%NAEREQ, CHI%NSV_CHSBEG, CHI%NSV_CHSEND,        &
                    CHI%NSV_AERBEG, CHI%NSV_AEREND, CHI%CCH_NAMES, CHI%CAER_NAMES, CHI%NDSTEQ,      &
                    CHI%NSV_DSTBEG, CHI%NSV_DSTEND, CHI%NSLTEQ, CHI%NSV_SLTBEG, CHI%NSV_SLTEND,     &
                    CHI%CDSTNAMES, CHI%CSLTNAMES, CHI%CCHEM_SURF_FILE,                      &
                    DST%XSFDST, DST%XSFDSTM, SLT%XSFSLT,                                    &
                    I%XAOSIP, I%XAOSIM, I%XAOSJP, I%XAOSJM, I%XHO2IP, I%XHO2IM, I%XHO2JP,     &
                    I%XHO2JM, I%XZ0, I%XZ0EFFIP, I%XZ0EFFIM, I%XZ0EFFJP, I%XZ0EFFJM, I%XZ0REL,&
                    I%XCLAY, I%XSAND, I%CPEDOTF,                                      &
                    I%XCONDSAT, I%XMPOTSAT, I%XBCOEF, I%XWWILT, I%XWFC, I%XWSAT, I%XWD0,      &
                    I%XKANISO, I%CRUNOFF,                                           &
                    I%XTAUICE, I%XCGSAT, I%XC1SAT, I%XC2REF, I%XC3, I%XC4B, I%XACOEF, I%XPCOEF, &
                    I%XC4REF, I%XPCPS, I%XPLVTT, I%XPLSTT,                              &
                    I%CSCOND, I%CISBA, I%XHCAPSOIL, I%XCONDDRY, I%XCONDSLD, I%CCPSURF,      &
                    I%XDG, I%XDROOT, I%XDG2, I%XROOTFRAC, I%XRUNOFFD, I%XDZG, I%XDZDIF,       &
                    I%XSOILWGHT, I%NWG_LAYER, I%NLAYER_HORT, I%NLAYER_DUN, I%XD_ICE,      &
                    I%XKSAT_ICE, I%XALBNIR_DRY, I%XALBVIS_DRY, I%XALBUV_DRY,            &
                    I%XALBNIR_WET, I%XALBVIS_WET, I%XALBUV_WET, I%XBSLAI_NITRO,         &
                    I%XCE_NITRO, I%XCNA_NITRO, I%XCF_NITRO, I%XFWTD, I%XWTD               )  
!
!-------------------------------------------------------------------------------
!
!DIF option :
!    Anisotropy coeficient for hydraulic conductivity for topmodel drainage (Fan et al. 2006)
!    Soil organic matter effect and/or Exponential decay for DIF option
!    Must be call before INIT_TOP
!
!
IF(I%CISBA=='DIF') THEN
  !
  IF( I%CKSAT=='SGH' )THEN
    WRITE(ILUOUT,*)'THE KSAT EXP PROFILE WITH ISBA-DF IS NOT PHYSIC AND HAS BEEN REMOVED FOR NOW' 
    WRITE(ILUOUT,*)'A NEW PHYSICAL APPROACH WILL BE DEVELLOPED ACCOUNTING FOR COMPACTION IN ALL ' 
    WRITE(ILUOUT,*)'HYDRODYNAMIC PARAMETERS (WSAT, PSISAT, KSAT, B) AND NOT ONLY IN KSAT        ' 
    CALL ABOR1_SFX('CKSAT=SGH is not physic with ISBA-DF and has been removed for now')    
  ENDIF
  !  
  IF(I%LSOC)THEN   
    IF(.NOT.I%LSOCP)THEN
      WRITE(ILUOUT,*)'LSOC = T can be activated only if SOC data given in PGD fields'
      CALL ABOR1_SFX('LSOC = T can be activated only if SOC data given in PGD fields')
    ENDIF
    ALLOCATE(I%XFRACSOC(KI,I%NGROUND_LAYER))
    I%XFRACSOC(:,:)=0.0
    CALL ISBA_SOC_PARAMETERS(I%CRUNOFF,I%XPATCH,I%XDG,I%XSOC,I%XBCOEF,I%XMPOTSAT,   &
                             I%XCONDSAT,I%XWSAT,I%XHCAPSOIL,I%XCONDDRY,         &
                             I%XCONDSLD,I%XWFC,I%XWWILT,I%XWD0,I%XKANISO,I%XFRACSOC )
  ELSE
    ALLOCATE(I%XFRACSOC(0,0))
  ENDIF
!
ELSE
  ALLOCATE(I%XFRACSOC(0,0))
ENDIF
!
!Topmodel
!
!CRUNOFF used in hydro_sgh and isba_sgh_update
IF( I%CRUNOFF=='SGH ') THEN 
!
  ALLOCATE(I%XTAB_FSAT(KI,NDIMTAB))
  ALLOCATE(I%XTAB_WTOP(KI,NDIMTAB))
  ALLOCATE(I%XTAB_QTOP(KI,NDIMTAB))
!
  I%XTAB_FSAT(:,:) = 0.0
  I%XTAB_WTOP(:,:) = 0.0
  I%XTAB_QTOP(:,:) = 0.0
!
  IF(HINIT/='PRE')THEN
!
    WHERE(I%XCLAY(:,1)==XUNDEF.AND.I%XTI_MEAN(:)/=XUNDEF) I%XTI_MEAN(:)=XUNDEF
!
    CALL INIT_TOP(I, &
                   I%CISBA, ILUOUT, I%XPATCH, I%XRUNOFFD,          &
                   I%XWD0, I%XWSAT, I%XTI_MIN,                     &
                   I%XTI_MAX, I%XTI_MEAN, I%XTI_STD, I%XTI_SKEW,     &
                   I%XSOILWGHT, I%XTAB_FSAT, I%XTAB_WTOP,          &
                   I%XTAB_QTOP, ZM                             )  
!
!
    IF (I%CKSAT=='SGH' .AND. I%CISBA/='DIF') THEN
!     Exponential decay factor calculate using soil properties 
!     (eq. 11, Decharme et al., J. Hydrometeor, 2006)
      DO JILU=1,KI
        IF (ZM(JILU)/=XUNDEF) ZF(JILU,:) = (I%XWSAT(JILU,1)-I%XWD0(JILU,1))/ZM(JILU)
      ENDDO
!       
    ENDIF
!
  ENDIF
!
! Subsurface flow by layer (m/s)
  IF(I%CISBA=='DIF') THEN
    ALLOCATE(I%XTOPQS(KI,I%NGROUND_LAYER,I%NPATCH))
    I%XTOPQS(:,:,:)=0.0
  ELSE
    ALLOCATE(I%XTOPQS(0,0,0))
  ENDIF
!
ELSE                  
!  
  ALLOCATE(I%XTAB_FSAT(0,0))
  ALLOCATE(I%XTAB_WTOP(0,0))
  ALLOCATE(I%XTAB_QTOP(0,0))
  ALLOCATE(I%XTOPQS(0,0,0))  
!                  
ENDIF  
!
!Exponential decay for ISBA-FR option
!CKSAT used in hydro_soil.F90 and soil.F90
IF(HINIT/='PRE'.AND.I%CISBA/='DIF')THEN 
  !
  IF(I%CKSAT=='SGH') THEN
    !
    WHERE(ZF(:,:)==XUNDEF.AND.I%XDG(:,2,:)/=XUNDEF) 
      ZF(:,:) = 4.0/I%XDG(:,2,:)
    ENDWHERE
    ZF(:,:) = MIN(ZF(:,:),XF_DECAY)
    !
    ALLOCATE(I%XF_PARAM (KI))
    I%XF_PARAM(:) = ZF(:,1)
    !
    DO JPATCH=1,I%NPATCH
      IF (I%NSIZE_NATURE_P(JPATCH) == 0 ) CYCLE
        CALL EXP_DECAY_SOIL_FR(I%CISBA, ZF(:,JPATCH),I%XC1SAT(:,JPATCH),I%XC2REF(:,JPATCH),   &
                                I%XDG(:,:,JPATCH),I%XD_ICE(:,JPATCH),I%XC4REF(:,JPATCH),      &
                                I%XC3(:,:,JPATCH),I%XCONDSAT(:,:,JPATCH),I%XKSAT_ICE(:,JPATCH))  
    ENDDO                       
    !
  ELSEIF ( I%CKSAT=='EXP' .AND. I%CISBA=='3-L' ) THEN
    !
    ALLOCATE(I%XF_PARAM (KI))
    I%XF_PARAM(:) = XUNDEF
    !
    IF (HPROGRAM/='AROME ' .AND. HPROGRAM/='MESONH ') THEN
      !
      CALL OPEN_FILE('ASCII ',NUNIT,HFILE='carte_f_dc.txt',HFORM='FORMATTED',HACTION='READ ')
      DO JILU=1,U%NDIM_FULL
        READ(NUNIT,*) ZF_PARAM(JILU), ZC_DEPTH_RATIO(JILU)
      ENDDO
      CALL CLOSE_FILE('ASCII ',NUNIT)
      CALL READ_AND_SEND_MPI(ZF_PARAM,I%XF_PARAM,U%NR_NATURE)
#ifdef TOPD
IF (.NOT.ALLOCATED(XC_DEPTH_RATIO))    ALLOCATE(XC_DEPTH_RATIO (KI))
    XC_DEPTH_RATIO(:) = XUNDEF
      CALL READ_AND_SEND_MPI(ZC_DEPTH_RATIO,XC_DEPTH_RATIO,U%NR_NATURE)
#endif
      !
    ELSE
      WRITE(ILUOUT,*) "COMPUTE_ISBA_PARAMETERS: WITH CKSAT=EXP, IN NOT OFFLINE "//&
                      "MODE, TOPMODEL FILE FOR F_PARAM IS NOT READ "
    ENDIF
    !
    DO JPATCH=1,I%NPATCH
      WHERE (I%XF_PARAM(:)==XUNDEF.AND.I%XDG(:,2,JPATCH)/=XUNDEF)
        ZF(:,JPATCH) = 4.0/I%XDG(:,2,JPATCH)
      ELSEWHERE
        ZF(:,JPATCH) = I%XF_PARAM(:)
      ENDWHERE
    ENDDO
     ZF(:,:) = MIN(ZF(:,:),XF_DECAY)
    !
    DO JPATCH=1,I%NPATCH
      CALL EXP_DECAY_SOIL_FR(I%CISBA, ZF(:,JPATCH),I%XC1SAT(:,JPATCH),I%XC2REF(:,JPATCH), &
                             I%XDG(:,:,JPATCH),I%XD_ICE(:,JPATCH),I%XC4REF(:,JPATCH),   &
                             I%XC3(:,:,JPATCH),I%XCONDSAT(:,:,JPATCH),                &
                             I%XKSAT_ICE(:,JPATCH))  
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
IF (HINIT == 'ALL' .AND. I%CRESPSL=='CNT' .AND. I%CPHOTO == 'NCB') THEN
  CALL CARBON_INIT(I%NNBIOMASS, I%NNLITTER, I%NNLITTLEVS, I%NNSOILCARB)
ENDIF
!
!Rainfall spatial distribution
!CRAIN used in HYDRO_VEG and HYDRO_SGH and ISBA_SGH_UPDATE
IF(I%CRAIN=='SGH')THEN
  ALLOCATE(I%XMUF(KI))
  I%XMUF(:)=0.0
ELSE
  ALLOCATE(I%XMUF(0))
ENDIF
!
ALLOCATE(I%XFSAT(KI))  
I%XFSAT(:) = 0.0
!
!-------------------------------------------------------------------------------
!
!*       6.2    Initialize of SFX - RRM coupling:
!               ---------------------------------
!
! * Check some key :
!
IF(LCPL_CALVING)THEN
   IF(.NOT.I%LGLACIER)THEN
     CALL ABOR1_SFX('COMPUTE_ISBA_PARAMETERS: LGLACIER MUST BE ACTIVATED IF LCPL_CALVING')
   ENDIF
ENDIF
!
! * Initialize required coupling fields :
!
I%LCPL_RRM = .FALSE.
I%LFLOOD   = .FALSE.
I%LWTD     = .FALSE.
!
IF(LCPL_LAND)THEN
!    
  I%LCPL_RRM = .TRUE.
!
  ALLOCATE(I%XCPL_DRAIN (KI))
  ALLOCATE(I%XCPL_RUNOFF(KI))
  I%XCPL_DRAIN (:) = 0.0
  I%XCPL_RUNOFF(:) = 0.0
!
  IF(I%LGLACIER)THEN
     ALLOCATE(I%XCPL_ICEFLUX(KI))
     I%XCPL_ICEFLUX(:) = 0.0
  ELSE
     ALLOCATE(I%XCPL_ICEFLUX(0))
  ENDIF
!
  IF(LCPL_GW)THEN
    I%LWTD = .TRUE.
    ALLOCATE(I%XCPL_RECHARGE(KI))
    I%XCPL_RECHARGE(:) = 0.0
  ELSE
    ALLOCATE(I%XCPL_RECHARGE(0))
  ENDIF
!
  IF(LCPL_FLOOD)THEN
     I%LFLOOD = .TRUE.
     ALLOCATE(I%XCPL_EFLOOD(KI))
     ALLOCATE(I%XCPL_PFLOOD(KI))
     ALLOCATE(I%XCPL_IFLOOD(KI))
     I%XCPL_EFLOOD(:)= 0.0
     I%XCPL_PFLOOD(:)= 0.0
     I%XCPL_IFLOOD(:)= 0.0    
  ELSE
    ALLOCATE(I%XCPL_EFLOOD(0))
    ALLOCATE(I%XCPL_PFLOOD(0))
    ALLOCATE(I%XCPL_IFLOOD(0))     
  ENDIF     
!
ELSE
!
  ALLOCATE(I%XCPL_RUNOFF  (0))
  ALLOCATE(I%XCPL_DRAIN   (0))
  ALLOCATE(I%XCPL_ICEFLUX (0))
  ALLOCATE(I%XCPL_RECHARGE(0))
  ALLOCATE(I%XCPL_EFLOOD  (0))
  ALLOCATE(I%XCPL_PFLOOD  (0))
  ALLOCATE(I%XCPL_IFLOOD  (0))
!
ENDIF
!
IF(I%LWTD.AND..NOT.I%LGW)THEN
  WRITE(ILUOUT,*)'COMPUTE_ISBA_PARAMETERS: Groundwater map is required by SFX - Groundwater coupling '
  WRITE(ILUOUT,*)'COMPUTE_ISBA_PARAMETERS: Please check your pgd namelist where this map must be     '
  WRITE(ILUOUT,*)'COMPUTE_ISBA_PARAMETERS: specified (YGW and YGWFILETYPE, or XUNIF_GW, or LIMP_GW)  '
  CALL ABOR1_SFX('COMPUTE_ISBA_PARAMETERS: Groundwater map is required by SFX - Groundwater coupling')
ENDIF
!
! * Initialize flood scheme :
!
IF(I%LFLOOD)THEN
  ALLOCATE(I%XFFLOOD (KI))
  ALLOCATE(I%XPIFLOOD(KI))
  ALLOCATE(I%XFF     (KI,I%NPATCH))
  ALLOCATE(I%XFFG    (KI,I%NPATCH))
  ALLOCATE(I%XFFV    (KI,I%NPATCH))
  ALLOCATE(I%XFFROZEN(KI,I%NPATCH))
  ALLOCATE(I%XALBF   (KI,I%NPATCH))
  ALLOCATE(I%XEMISF  (KI,I%NPATCH)) 
  I%XFFLOOD       = 0.0
  I%XPIFLOOD      = 0.0
  I%XFF           = 0.0
  I%XFFG          = 0.0
  I%XFFV          = 0.0
  I%XFFROZEN      = 0.0
  I%XALBF         = 0.0
  I%XEMISF        = 0.0  
ELSE
  ALLOCATE(I%XFFLOOD   (0))
  ALLOCATE(I%XPIFLOOD  (0))
  ALLOCATE(I%XFF     (0,0))
  ALLOCATE(I%XFFG    (0,0))
  ALLOCATE(I%XFFV    (0,0))
  ALLOCATE(I%XFFROZEN(0,0))
  ALLOCATE(I%XALBF   (0,0))
  ALLOCATE(I%XEMISF  (0,0))
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      7.     ISBA time-varying deep force-restore temperature initialization
!              ---------------------------------------------------------------
!
 CALL SOILTEMP_ARP_PAR(I, &
                       HPROGRAM,I%LTEMP_ARP,I%NTEMPLAYER_ARP)
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
IF (HINIT=='PRE' .AND. I%TSNOW%SCHEME.NE.'3-L' .AND. I%TSNOW%SCHEME.NE.'CRO' .AND. I%CISBA=='DIF') THEN
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
 CALL READ_ISBA_CANOPY_n(DTCO, IOB, ICP, I, U, &
                         HPROGRAM)
!
!-------------------------------------------------------------------------------
!
!*      13.     initialize radiative and physical properties
!               --------------------------------------------
!
ALLOCATE(I%XDIR_ALB_WITH_SNOW(KI,KSW,I%NPATCH))
ALLOCATE(I%XSCA_ALB_WITH_SNOW(KI,KSW,I%NPATCH))
I%XDIR_ALB_WITH_SNOW = 0.0
I%XSCA_ALB_WITH_SNOW = 0.0
!
 CALL INIT_VEG_n(I%NPATCH, KI, I%LCANOPY, I%CROUGH, I%LAGRI_TO_GRASS, I%TSNOW, &
                 I%CPHOTO, I%XLAIMIN, I%XH_TREE, I%XVEGTYPE_PATCH, I%XLAI, I%XZ0, I%XVEG, I%XEMIS, &
                 I%LTR_ML, I%XFAPARC, I%XFAPIRC, I%XLAI_EFFC, I%XMUS, &
                 I%XALBNIR_SOIL, I%XALBVIS_SOIL, I%XALBUV_SOIL, I%XALBNIR, I%XALBVIS, I%XALBUV, &
                 DGMI%LSURF_DIAG_ALBEDO, I%XPSN, I%XPSNG, I%XPSNV, I%XPSNV_A, &
                 PDIR_ALB, PSCA_ALB, PEMIS, PTSRAD )
!
DO JPATCH=1,I%NPATCH
  ZWG1(:,JPATCH) = I%XWG(:,1,JPATCH)
  ZTG1(:,JPATCH) = I%XTG(:,1,JPATCH)
END DO
!
 CALL CONVERT_PATCH_ISBA(DTCO, DTI, I, &
                         I%CISBA,IDECADE,IDECADE2,I%XCOVER,I%LCOVER,&
                          I%CPHOTO,LAGRIP,I%LPERM,I%LTR_ML,'NAT',   &
                          PWG1 = ZWG1,               &
                          PALBNIR_SOIL=I%XALBNIR_SOIL, &
                          PALBVIS_SOIL=I%XALBVIS_SOIL, &
                          PALBUV_SOIL=I%XALBUV_SOIL )
!
! Load randomly perturbed fields. Perturbation ratios are saved in case fields are reset later.
IF(I%LPERTSURF) THEN
!
  CALL READ_SURF(IOB, &
                 HPROGRAM,'VEG',I%XVEG(:,:),IRESP)
  ALLOCATE(I%XPERTVEG(KI))
  I%XPERTVEG(:)=I%XVEG(:,1)
!
  CALL READ_SURF(IOB, &
                 HPROGRAM,'LAI',I%XLAI(:,:),IRESP)
  ALLOCATE(I%XPERTLAI(KI))
  I%XPERTLAI(:)=I%XLAI(:,1)
!
  CALL READ_SURF(IOB, &
                 HPROGRAM,'CV',I%XCV(:,:),IRESP)
  ALLOCATE(I%XPERTCV(KI))
  I%XPERTCV(:)=I%XCV(:,1)
!
  CALL READ_SURF(IOB, &
                 HPROGRAM,'PERTALB',ZPERTBUF(:,:),IRESP)
  ALLOCATE(I%XPERTALB(KI))
  I%XPERTALB(:)=ZPERTBUF(:,1)
  WHERE(I%XALBNIR_VEG(:,1)/=XUNDEF)  I%XALBNIR_VEG(:,1) = I%XALBNIR_VEG(:,1) *( 1.+ I%XPERTALB(:) )
  WHERE(I%XALBVIS_VEG(:,1)/=XUNDEF)  I%XALBVIS_VEG(:,1) = I%XALBVIS_VEG(:,1) *( 1.+ I%XPERTALB(:) )
  WHERE(I%XALBUV_VEG(:,1)/=XUNDEF)   I%XALBUV_VEG(:,1)  = I%XALBUV_VEG(:,1)  *( 1.+ I%XPERTALB(:) )
  WHERE(I%XALBNIR_SOIL(:,1)/=XUNDEF) I%XALBNIR_SOIL(:,1)= I%XALBNIR_SOIL(:,1)*( 1.+ I%XPERTALB(:) )
  WHERE(I%XALBVIS_SOIL(:,1)/=XUNDEF) I%XALBVIS_SOIL(:,1)= I%XALBVIS_SOIL(:,1)*( 1.+ I%XPERTALB(:) )
  WHERE(I%XALBUV_SOIL(:,1)/=XUNDEF)  I%XALBUV_SOIL(:,1) = I%XALBUV_SOIL(:,1) *( 1.+ I%XPERTALB(:) )
!
  CALL READ_SURF(IOB, &
                 HPROGRAM,'PERTZ0LAND',ZPERTBUF(:,:),IRESP)
  ALLOCATE(I%XPERTZ0(KI))
  I%XPERTZ0(:)=ZPERTBUF(:,1)
  WHERE(I%XZ0(:,1)/=XUNDEF)      I%XZ0(:,1)     =I%XZ0(:,1)     *( 1.+ I%XPERTZ0(:) )
  WHERE(I%XZ0EFFIP(:,1)/=XUNDEF) I%XZ0EFFIP(:,1)=I%XZ0EFFIP(:,1)*( 1.+ I%XPERTZ0(:) )
  WHERE(I%XZ0EFFIM(:,1)/=XUNDEF) I%XZ0EFFIM(:,1)=I%XZ0EFFIM(:,1)*( 1.+ I%XPERTZ0(:) )
  WHERE(I%XZ0EFFJP(:,1)/=XUNDEF) I%XZ0EFFJP(:,1)=I%XZ0EFFJP(:,1)*( 1.+ I%XPERTZ0(:) )
  WHERE(I%XZ0EFFJM(:,1)/=XUNDEF) I%XZ0EFFJM(:,1)=I%XZ0EFFJM(:,1)*( 1.+ I%XPERTZ0(:) )
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       14.    Output radiative fields
!               -----------------------
!
ALLOCATE(I%XEMIS_NAT   (KI))
I%XEMIS_NAT (:) = XUNDEF
!
 CALL AVERAGED_ALBEDO_EMIS_ISBA(I, &
                                I%LFLOOD, I%CALBEDO, PZENITH,                &
                                 I%XVEG,I%XZ0,I%XLAI,                          &
                                 I%LMEB_PATCH,I%XGNDLITTER,I%XZ0LITTER,I%XLAIGV, &
                                 I%XZF_TALLVEG, I%XH_VEG, I%XTV,               &
                                 ZTG1,                                   &
                                 I%XPATCH,                                 &
                                 PSW_BANDS,                              &
                                 I%XALBNIR_VEG,I%XALBVIS_VEG,I%XALBUV_VEG,     &
                                 I%XALBNIR_SOIL,I%XALBVIS_SOIL,I%XALBUV_SOIL,  &
                                 I%XEMIS,                                  &
                                 I%TSNOW,                                  &
                                 I%XALBNIR,I%XALBVIS,I%XALBUV,                 &
                                 PDIR_ALB, PSCA_ALB,                     &
                                 I%XEMIS_NAT,ZTSRAD_NAT,ZTSURF_NAT         )  
!
PEMIS  = I%XEMIS_NAT
PTSRAD = ZTSRAD_NAT
PTSURF = ZTSURF_NAT
!
!-------------------------------------------------------------------------------
!
!*      15.     ISBA diagnostics initialization
!               -------------------------------
!
IF(I%NPATCH<=1) DGI%LPATCH_BUDGET=.FALSE.
!
 CALL DIAG_ISBA_INIT_n(IOB, &
                       CHI, DGEI, DGI, DGMI, DGU, GB, I, &
                       HPROGRAM,KI,KSW)
!
!-------------------------------------------------------------------------------
!
 CALL INIT_SURF_TOPD(DGEI, I, UG, U, &
                     HPROGRAM,U%NDIM_FULL)
!
!-------------------------------------------------------------------------------
!
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('COMPUTE_ISBA_PARAMETERS',1,ZHOOK_HANDLE)
!
END SUBROUTINE COMPUTE_ISBA_PARAMETERS


