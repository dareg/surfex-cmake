!#############################################################
SUBROUTINE COMPUTE_ISBA_PARAMETERS (DTCO, DGU, UG, U, IM, DST, SLT, SV, &
                                    HPROGRAM,HINIT,OLAND_USE,            &
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
USE MODD_SURFEX_n, ONLY : ISBA_MODEL_t
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DIAG_SURF_ATM_n, ONLY : DIAG_SURF_ATM_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_DST_n, ONLY : DST_t
USE MODD_SLT_n, ONLY : SLT_t
USE MODD_SV_n, ONLY : SV_t
!
USE MODD_SFX_OASIS,  ONLY : LCPL_LAND, LCPL_FLOOD, LCPL_GW, LCPL_CALVING
!
!
#ifdef TOPD
USE MODD_DUMMY_EXP_PROFILE,ONLY : XC_DEPTH_RATIO
#endif
!
USE MODD_ASSIM, ONLY : CASSIM_ISBA, LASSIM
!
USE MODD_DEEPSOIL,       ONLY : LPHYSDOMC, LDEEPSOIL, XTDEEP_CLI, XGAMMAT_CLI
USE MODD_AGRI,           ONLY : LAGRIP, XTHRESHOLD
!
!
USE MODD_SGH_PAR,        ONLY : NDIMTAB, XICE_DEPH_MAX, XF_DECAY
!
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE
USE MODD_SURF_PAR,       ONLY : XUNDEF, NUNDEF
USE MODD_SNOW_PAR,       ONLY : XEMISSN
!
USE MODD_TOPD_PAR, ONLY : NUNIT
USE MODD_TOPODYN, ONLY : NNCAT, NMESHT
!
!
USE MODE_RANDOM
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
!
TYPE(ISBA_MODEL_t), INTENT(INOUT) :: IM
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DIAG_SURF_ATM_t), INTENT(INOUT) :: DGU
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(DST_t), INTENT(INOUT) :: DST
TYPE(SLT_t), INTENT(INOUT) :: SLT
TYPE(SV_t), INTENT(INOUT) :: SV
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
REAL, DIMENSION(KI,IM%I%O%NPATCH)  :: ZWG1 ! work array for surface water content
REAL, DIMENSION(KI,IM%I%O%NPATCH)  :: ZTG1 ! work array for surface temperature
!
REAL, DIMENSION(KI)         :: ZM, ZWORK
REAL, DIMENSION(KI,IM%I%O%NPATCH)  :: ZF, ZPERTBUF
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
 CALL ALLOCATE_PHYSIO(IM%I, &
                      IM%I%O%CPHOTO, IM%I%O%CISBA, KI, NVEGTYPE, IM%I%O%NGROUND_LAYER, IM%I%O%NPATCH, &
                     IM%I%M%X%XVEGTYPE, IM%I%M%T%XLAI, IM%I%M%T%XVEG, IM%I%M%T%XZ0, IM%I%M%T%XEMIS, &
                     IM%I%M%X%XDG, IM%I%M%X%XD_ICE,      &
                     IM%I%M%T%XRSMIN, IM%I%M%T%XGAMMA, IM%I%M%T%XWRMAX_CF, IM%I%M%T%XRGL, IM%I%M%T%XCV,               &
                     IM%I%M%X%XZ0_O_Z0H, IM%I%M%T%XALBNIR_VEG, IM%I%M%T%XALBVIS_VEG, IM%I%M%T%XALBUV_VEG,    &
                     IM%I%M%X%XH_TREE, IM%I%M%X%XRE25, IM%I%M%T%XLAIMIN, IM%I%M%T%XBSLAI, IM%I%M%T%XSEFOLD,           &
                     IM%I%M%T%XGMES, IM%I%M%T%XGC, IM%I%M%T%XF2I, IM%I%M%X%XDMAX, IM%I%M%T%LSTRESS,                   &
                     IM%I%M%T%XCE_NITRO, IM%I%M%T%XCF_NITRO, IM%I%M%T%XCNA_NITRO,                   &
                     IM%I%M%I%TSEED, IM%I%M%I%TREAP, IM%I%M%I%XWATSUP, IM%I%M%I%XIRRIG,                      &
                     IM%I%M%X%XROOTFRAC, IM%I%M%X%NWG_LAYER, IM%I%M%X%XDROOT, IM%I%M%X%XDG2,                 &
                     IM%I%M%M%XGNDLITTER,IM%I%M%M%XRGLGV,IM%I%M%M%XGAMMAGV,IM%I%M%M%XRSMINGV,    &
                     IM%I%M%M%XROOTFRACGV,IM%I%M%M%XWRMAX_CFGV,IM%I%M%M%XLAIGV,IM%I%M%M%XZ0LITTER,IM%I%M%M%XH_VEG     )
!
IF (IM%I%I%TTIME%TDATE%MONTH /= NUNDEF) THEN
  IDECADE = 3 * ( IM%I%I%TTIME%TDATE%MONTH - 1 ) + MIN(IM%I%I%TTIME%TDATE%DAY-1,29) / 10 + 1
ELSE
  IDECADE = 1
END IF
!
IDECADE2 = IDECADE
!
 CALL INIT_ISBA_MIXPAR(DTCO, IM%DTI, IM%IG, IM%I%O, &
                       IM%I%O%CISBA,IDECADE,IDECADE2,IM%I%P%XCOVER,IM%I%P%LCOVER,IM%I%O%CPHOTO,'NAT')
!
ISIZE_LMEB_PATCH=COUNT(IM%I%O%LMEB_PATCH(:))
IF (ISIZE_LMEB_PATCH>0)  THEN
  CALL FIX_MEB_VEG(IM%DTI, IM%IG, IM%I%O%LMEB_PATCH, &
                   IM%I%O%NPATCH)
ENDIF
!
 CALL CONVERT_PATCH_ISBA(DTCO, IM%DTI, IM%I%O, &
                         IM%I%O%CISBA,IDECADE,IDECADE2,IM%I%P%XCOVER,IM%I%P%LCOVER,IM%I%O%CPHOTO,LAGRIP,   &
                        IM%I%O%LPERM,IM%I%O%LTR_ML,'NAT',PVEG=IM%I%M%T%XVEG,PLAI=IM%I%M%T%XLAI,                &
                        PRSMIN=IM%I%M%T%XRSMIN,PGAMMA=IM%I%M%T%XGAMMA,PWRMAX_CF=IM%I%M%T%XWRMAX_CF,       &
                        PRGL=IM%I%M%T%XRGL,PCV=IM%I%M%T%XCV,PSOILGRID=IM%I%O%XSOILGRID,                 &
                        PDG=IM%I%M%X%XDG,KWG_LAYER=IM%I%M%X%NWG_LAYER,PDROOT=IM%I%M%X%XDROOT,PDG2=IM%I%M%X%XDG2,   &
                        PZ0=IM%I%M%T%XZ0,PZ0_O_Z0H=IM%I%M%X%XZ0_O_Z0H,                           &
                        PALBNIR_VEG=IM%I%M%T%XALBNIR_VEG,PALBVIS_VEG=IM%I%M%T%XALBVIS_VEG,       &
                        PALBUV_VEG=IM%I%M%T%XALBUV_VEG,PEMIS_ECO=IM%I%M%T%XEMIS,                 &
                        PVEGTYPE=IM%I%M%X%XVEGTYPE,PROOTFRAC=IM%I%M%X%XROOTFRAC,                 &
                        PGMES=IM%I%M%T%XGMES,PBSLAI=IM%I%M%T%XBSLAI,PLAIMIN=IM%I%M%T%XLAIMIN,             &
                        PSEFOLD=IM%I%M%T%XSEFOLD,PGC=IM%I%M%T%XGC, PPERM=IM%I%P%XPERM,               &
                        PDMAX=IM%I%M%X%XDMAX,PF2I=IM%I%M%T%XF2I,OSTRESS=IM%I%M%T%LSTRESS,PH_TREE=IM%I%M%X%XH_TREE, &
                        PRE25=IM%I%M%X%XRE25,PCE_NITRO=IM%I%M%T%XCE_NITRO,PCF_NITRO=IM%I%M%T%XCF_NITRO,   &
                        PCNA_NITRO=IM%I%M%T%XCNA_NITRO,PD_ICE=IM%I%M%X%XD_ICE,TPSEED=IM%I%M%I%TSEED,      &
                        TPREAP=IM%I%M%I%TREAP,PWATSUP=IM%I%M%I%XWATSUP,PIRRIG=IM%I%M%I%XIRRIG,            &
                        PGNDLITTER=IM%I%M%M%XGNDLITTER,                                           &
                        PRGLGV=IM%I%M%M%XRGLGV,   &
                        PGAMMAGV=IM%I%M%M%XGAMMAGV,PRSMINGV=IM%I%M%M%XRSMINGV,                   &
                        PROOTFRACGV=IM%I%M%M%XROOTFRACGV,PWRMAX_CFGV=IM%I%M%M%XWRMAX_CFGV,       &
                        PLAIGV=IM%I%M%M%XLAIGV,PZ0LITTER=IM%I%M%M%XZ0LITTER,PH_VEG=IM%I%M%M%XH_VEG        )
!
!-------------------------------------------------------------------------------
!
CALL INIT_VEG_PGD_n(IM%CHI, DTCO, DST, SLT, U, &
                    IM%I%O%LAGRI_TO_GRASS, IM%I%P%LCOVER, IM%I%P%XCOVER, & 
                    HPROGRAM, 'NATURE', ILUOUT, KI, IM%I%O%NPATCH, IM%I%O%NGROUND_LAYER,      &
                    IM%I%I%TTIME%TDATE%MONTH,                                          &
                    IM%I%M%X%XVEGTYPE, IM%I%IP%XPATCH, IM%I%IP%XVEGTYPE_PATCH, IM%I%IP%NSIZE_NATURE_P,           &
                    IM%I%IP%NR_NATURE_P, IM%I%O%XRM_PATCH,                                     &
                    LDEEPSOIL, LPHYSDOMC, XTDEEP_CLI, XGAMMAT_CLI, IM%I%IP%XTDEEP,      &
                    IM%I%IP%XGAMMAT, LAGRIP, XTHRESHOLD, IM%AG%NIRRINUM, IM%AG%LIRRIDAY, IM%AG%LIRRIGATE, &
                    IM%AG%XTHRESHOLDSPT,                                              &
                    IM%I%O%CPHOTO, HINIT, IM%I%O%LTR_ML, IM%I%O%NNBIOMASS, PCO2, PRHOA, IM%I%IP%XABC, IM%I%IP%XPOI,  &
                    IM%I%M%T%XGMES, IM%I%M%T%XGC, IM%I%M%X%XDMAX, IM%I%IP%XANMAX, IM%I%IP%XFZERO, &
                    IM%I%IP%XEPSO, IM%I%IP%XGAMM, IM%I%IP%XQDGAMM,   & 
                    IM%I%IP%XQDGMES, IM%I%IP%XT1GMES, IM%I%IP%XT2GMES, IM%I%IP%XAMAX, &
                    IM%I%IP%XQDAMAX, IM%I%IP%XT1AMAX, IM%I%IP%XT2AMAX,&
                    IM%I%IP%XAH, IM%I%IP%XBH, IM%I%IP%XTAU_WOOD, IM%I%IP%XINCREASE, IM%I%IP%XTURNOVER,                  &
                    KSV, HSV, IM%CHI%SVI, IM%CHI%CCH_NAMES, IM%CHI%CAER_NAMES, &
                    IM%CHI%CDSTNAMES, IM%CHI%CSLTNAMES, IM%CHI%CCHEM_SURF_FILE,                      &
                    DST%XSFDST, DST%XSFDSTM, SLT%XSFSLT,                                    &
                    IM%I%P%XAOSIP, IM%I%P%XAOSIM, IM%I%P%XAOSJP, IM%I%P%XAOSJM, IM%I%P%XHO2IP, &
                    IM%I%P%XHO2IM, IM%I%P%XHO2JP,     &
                    IM%I%P%XHO2JM, IM%I%M%T%XZ0, IM%I%IP%XZ0EFFIP, IM%I%IP%XZ0EFFIM, IM%I%IP%XZ0EFFJP, &
                    IM%I%IP%XZ0EFFJM, IM%I%IP%XZ0REL,&
                    IM%I%P%XCLAY, IM%I%P%XSAND, IM%I%O%CPEDOTF,                                      &
                    IM%I%IP%XCONDSAT, IM%I%IP%XMPOTSAT, IM%I%IP%XBCOEF, IM%I%IP%XWWILT, IM%I%IP%XWFC, &
                    IM%I%IP%XWSAT, IM%I%IP%XWD0,      &
                    IM%I%IP%XKANISO, IM%I%O%CRUNOFF,                                           &
                    IM%I%IP%XTAUICE, IM%I%IP%XCGSAT, IM%I%IP%XC1SAT, IM%I%IP%XC2REF, IM%I%IP%XC3, &
                    IM%I%IP%XC4B, IM%I%IP%XACOEF, IM%I%IP%XPCOEF, &
                    IM%I%IP%XC4REF, IM%I%IP%XPCPS, IM%I%IP%XPLVTT, IM%I%IP%XPLSTT,                              &
                    IM%I%O%CSCOND, IM%I%O%CISBA, IM%I%IP%XHCAPSOIL, IM%I%IP%XCONDDRY, IM%I%IP%XCONDSLD, &
                    IM%I%O%CCPSURF,      &
                    IM%I%M%X%XDG, IM%I%M%X%XDROOT, IM%I%M%X%XDG2, IM%I%M%X%XROOTFRAC, IM%I%IP%XRUNOFFD, &
                    IM%I%IP%XDZG, IM%I%IP%XDZDIF,       &
                    IM%I%IP%XSOILWGHT, IM%I%M%X%NWG_LAYER, IM%I%O%NLAYER_HORT, IM%I%O%NLAYER_DUN, IM%I%M%X%XD_ICE,      &
                    IM%I%IP%XKSAT_ICE, IM%I%IP%XALBNIR_DRY, IM%I%IP%XALBVIS_DRY, IM%I%IP%XALBUV_DRY,            &
                    IM%I%IP%XALBNIR_WET, IM%I%IP%XALBVIS_WET, IM%I%IP%XALBUV_WET, IM%I%IP%XBSLAI_NITRO,         &
                    IM%I%M%T%XCE_NITRO, IM%I%M%T%XCNA_NITRO, IM%I%M%T%XCF_NITRO, IM%I%IP%XFWTD, IM%I%IP%XWTD               )  
!
!-------------------------------------------------------------------------------
!
!DIF option :
!    Anisotropy coeficient for hydraulic conductivity for topmodel drainage (Fan et al. 2006)
!    Soil organic matter effect and/or Exponential decay for DIF option
!    Must be call before INIT_TOP
!
!
IF(IM%I%O%CISBA=='DIF') THEN
  !
  IF( IM%I%O%CKSAT=='SGH' )THEN
    WRITE(ILUOUT,*)'THE KSAT EXP PROFILE WITH ISBA-DF IS NOT PHYSIC AND HAS BEEN REMOVED FOR NOW' 
    WRITE(ILUOUT,*)'A NEW PHYSICAL APPROACH WILL BE DEVELLOPED ACCOUNTING FOR COMPACTION IN ALL ' 
    WRITE(ILUOUT,*)'HYDRODYNAMIC PARAMETERS (WSAT, PSISAT, KSAT, B) AND NOT ONLY IN KSAT        ' 
    CALL ABOR1_SFX('CKSAT=SGH is not physic with ISBA-DF and has been removed for now')    
  ENDIF
  !  
  IF(IM%I%O%LSOC)THEN   
    IF(.NOT.IM%I%O%LSOCP)THEN
      WRITE(ILUOUT,*)'LSOC = T can be activated only if SOC data given in PGD fields'
      CALL ABOR1_SFX('LSOC = T can be activated only if SOC data given in PGD fields')
    ENDIF
    ALLOCATE(IM%I%I%XFRACSOC(KI,IM%I%O%NGROUND_LAYER))
    IM%I%I%XFRACSOC(:,:)=0.0
    CALL ISBA_SOC_PARAMETERS(IM%I%O%CRUNOFF,IM%I%IP%XPATCH,IM%I%M%X%XDG,IM%I%P%XSOC,IM%I%IP%XBCOEF,IM%I%IP%XMPOTSAT,   &
                             IM%I%IP%XCONDSAT,IM%I%IP%XWSAT,IM%I%IP%XHCAPSOIL,IM%I%IP%XCONDDRY,         &
                             IM%I%IP%XCONDSLD,IM%I%IP%XWFC,IM%I%IP%XWWILT,IM%I%IP%XWD0,IM%I%IP%XKANISO,IM%I%I%XFRACSOC )
  ELSE
    ALLOCATE(IM%I%I%XFRACSOC(0,0))
  ENDIF
!
ELSE
  ALLOCATE(IM%I%I%XFRACSOC(0,0))
ENDIF
!
!Topmodel
!
!CRUNOFF used in hydro_sgh and isba_sgh_update
IF( IM%I%O%CRUNOFF=='SGH ') THEN 
!
  ALLOCATE(IM%I%I%XTAB_FSAT(KI,NDIMTAB))
  ALLOCATE(IM%I%I%XTAB_WTOP(KI,NDIMTAB))
  ALLOCATE(IM%I%I%XTAB_QTOP(KI,NDIMTAB))
!
  IM%I%I%XTAB_FSAT(:,:) = 0.0
  IM%I%I%XTAB_WTOP(:,:) = 0.0
  IM%I%I%XTAB_QTOP(:,:) = 0.0
!
  IF(HINIT/='PRE' .AND. .NOT.LASSIM)THEN
!
    WHERE(IM%I%P%XCLAY(:,1)==XUNDEF.AND.IM%I%P%XTI_MEAN(:)/=XUNDEF) IM%I%P%XTI_MEAN(:)=XUNDEF
!
    CALL INIT_TOP(IM%I%O%NPATCH, IM%I%IP%NSIZE_NATURE_P, &
                   IM%I%O%CISBA, ILUOUT, IM%I%IP%XPATCH, IM%I%IP%XRUNOFFD,          &
                   IM%I%IP%XWD0, IM%I%IP%XWSAT, IM%I%P%XTI_MIN,                     &
                   IM%I%P%XTI_MAX, IM%I%P%XTI_MEAN, IM%I%P%XTI_STD, IM%I%P%XTI_SKEW,     &
                   IM%I%IP%XSOILWGHT, IM%I%I%XTAB_FSAT, IM%I%I%XTAB_WTOP,          &
                   IM%I%I%XTAB_QTOP, ZM                             )  
!
!
    IF (IM%I%O%CKSAT=='SGH' .AND. IM%I%O%CISBA/='DIF') THEN
!     Exponential decay factor calculate using soil properties 
!     (eq. 11, Decharme et al., J. Hydrometeor, 2006)
      DO JILU=1,KI
        IF (ZM(JILU)/=XUNDEF) ZF(JILU,:) = (IM%I%IP%XWSAT(JILU,1)-IM%I%IP%XWD0(JILU,1))/ZM(JILU)
      ENDDO
!       
    ENDIF
!
  ENDIF
!
! Subsurface flow by layer (m/s)
  IF(IM%I%O%CISBA=='DIF') THEN
    ALLOCATE(IM%I%I%XTOPQS(KI,IM%I%O%NGROUND_LAYER,IM%I%O%NPATCH))
    IM%I%I%XTOPQS(:,:,:)=0.0
  ELSE
    ALLOCATE(IM%I%I%XTOPQS(0,0,0))
  ENDIF
!
ELSE                  
!  
  ALLOCATE(IM%I%I%XTAB_FSAT(0,0))
  ALLOCATE(IM%I%I%XTAB_WTOP(0,0))
  ALLOCATE(IM%I%I%XTAB_QTOP(0,0))
  ALLOCATE(IM%I%I%XTOPQS(0,0,0))  
!                  
ENDIF  
!
!Exponential decay for ISBA-FR option
!CKSAT used in hydro_soil.F90 and soil.F90
IF(HINIT/='PRE'.AND.IM%I%O%CISBA/='DIF')THEN 
  !
  IF(IM%I%O%CKSAT=='SGH') THEN
    !
    WHERE(ZF(:,:)==XUNDEF.AND.IM%I%M%X%XDG(:,2,:)/=XUNDEF) 
      ZF(:,:) = 4.0/IM%I%M%X%XDG(:,2,:)
    ENDWHERE
    ZF(:,:) = MIN(ZF(:,:),XF_DECAY)
    !
    ALLOCATE(IM%I%I%XF_PARAM (KI))
    IM%I%I%XF_PARAM(:) = ZF(:,1)
    !
    DO JPATCH=1,IM%I%O%NPATCH
      IF (IM%I%IP%NSIZE_NATURE_P(JPATCH) == 0 ) CYCLE
        CALL EXP_DECAY_SOIL_FR(IM%I%O%CISBA, ZF(:,JPATCH),IM%I%IP%XC1SAT(:,JPATCH),IM%I%IP%XC2REF(:,JPATCH),   &
                                IM%I%M%X%XDG(:,:,JPATCH),IM%I%M%X%XD_ICE(:,JPATCH),IM%I%IP%XC4REF(:,JPATCH),      &
                                IM%I%IP%XC3(:,:,JPATCH),IM%I%IP%XCONDSAT(:,:,JPATCH),IM%I%IP%XKSAT_ICE(:,JPATCH))  
    ENDDO                       
    !
  ELSEIF ( IM%I%O%CKSAT=='EXP' .AND. IM%I%O%CISBA=='3-L' ) THEN
    !
    ALLOCATE(IM%I%I%XF_PARAM (KI))
    IM%I%I%XF_PARAM(:) = XUNDEF
    !
    IF (HPROGRAM/='AROME ' .AND. HPROGRAM/='MESONH ') THEN
      !
      CALL OPEN_FILE('ASCII ',NUNIT,HFILE='carte_f_dc.txt',HFORM='FORMATTED',HACTION='READ ')
      DO JILU=1,U%NDIM_FULL
        READ(NUNIT,*) ZF_PARAM(JILU), ZC_DEPTH_RATIO(JILU)
      ENDDO
      CALL CLOSE_FILE('ASCII ',NUNIT)
      CALL READ_AND_SEND_MPI(ZF_PARAM,IM%I%I%XF_PARAM,U%NR_NATURE)
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
    DO JPATCH=1,IM%I%O%NPATCH
      WHERE (IM%I%I%XF_PARAM(:)==XUNDEF.AND.IM%I%M%X%XDG(:,2,JPATCH)/=XUNDEF)
        ZF(:,JPATCH) = 4.0/IM%I%M%X%XDG(:,2,JPATCH)
      ELSEWHERE
        ZF(:,JPATCH) = IM%I%I%XF_PARAM(:)
      ENDWHERE
    ENDDO
     ZF(:,:) = MIN(ZF(:,:),XF_DECAY)
    !
    DO JPATCH=1,IM%I%O%NPATCH
      CALL EXP_DECAY_SOIL_FR(IM%I%O%CISBA, ZF(:,JPATCH),IM%I%IP%XC1SAT(:,JPATCH),IM%I%IP%XC2REF(:,JPATCH), &
                             IM%I%M%X%XDG(:,:,JPATCH),IM%I%M%X%XD_ICE(:,JPATCH),IM%I%IP%XC4REF(:,JPATCH),   &
                             IM%I%IP%XC3(:,:,JPATCH),IM%I%IP%XCONDSAT(:,:,JPATCH),                &
                             IM%I%IP%XKSAT_ICE(:,JPATCH))  
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
IF (HINIT == 'ALL' .AND. IM%I%O%CRESPSL=='CNT' .AND. IM%I%O%CPHOTO == 'NCB') THEN
  CALL CARBON_INIT(IM%I%O%NNBIOMASS, IM%I%O%NNLITTER, IM%I%O%NNLITTLEVS, IM%I%O%NNSOILCARB)
ENDIF
!
!Rainfall spatial distribution
!CRAIN used in HYDRO_VEG and HYDRO_SGH and ISBA_SGH_UPDATE
IF(IM%I%O%CRAIN=='SGH')THEN
  ALLOCATE(IM%I%I%XMUF(KI))
  IM%I%I%XMUF(:)=0.0
ELSE
  ALLOCATE(IM%I%I%XMUF(0))
ENDIF
!
ALLOCATE(IM%I%I%XFSAT(KI))  
IM%I%I%XFSAT(:) = 0.0
!
!-------------------------------------------------------------------------------
!
!*       6.2    Initialize of SFX - RRM coupling:
!               ---------------------------------
!
! * Check some key :
!
IF(LCPL_CALVING)THEN
   IF(.NOT.IM%I%O%LGLACIER)THEN
     CALL ABOR1_SFX('COMPUTE_ISBA_PARAMETERS: LGLACIER MUST BE ACTIVATED IF LCPL_CALVING')
   ENDIF
ENDIF
!
! * Initialize required coupling fields :
!
IM%I%O%LCPL_RRM = .FALSE.
IM%I%O%LFLOOD   = .FALSE.
IM%I%O%LWTD     = .FALSE.
!
IF(LCPL_LAND)THEN
!    
  IM%I%O%LCPL_RRM = .TRUE.
!
  ALLOCATE(IM%I%I%XCPL_DRAIN (KI))
  ALLOCATE(IM%I%I%XCPL_RUNOFF(KI))
  IM%I%I%XCPL_DRAIN (:) = 0.0
  IM%I%I%XCPL_RUNOFF(:) = 0.0
!
  IF(IM%I%O%LGLACIER)THEN
     ALLOCATE(IM%I%I%XCPL_ICEFLUX(KI))
     IM%I%I%XCPL_ICEFLUX(:) = 0.0
  ELSE
     ALLOCATE(IM%I%I%XCPL_ICEFLUX(0))
  ENDIF
!
  IF(LCPL_GW)THEN
    IM%I%O%LWTD = .TRUE.
    ALLOCATE(IM%I%I%XCPL_RECHARGE(KI))
    IM%I%I%XCPL_RECHARGE(:) = 0.0
  ELSE
    ALLOCATE(IM%I%I%XCPL_RECHARGE(0))
  ENDIF
!
  IF(LCPL_FLOOD)THEN
     IM%I%O%LFLOOD = .TRUE.
     ALLOCATE(IM%I%I%XCPL_EFLOOD(KI))
     ALLOCATE(IM%I%I%XCPL_PFLOOD(KI))
     ALLOCATE(IM%I%I%XCPL_IFLOOD(KI))
     IM%I%I%XCPL_EFLOOD(:)= 0.0
     IM%I%I%XCPL_PFLOOD(:)= 0.0
     IM%I%I%XCPL_IFLOOD(:)= 0.0    
  ELSE
    ALLOCATE(IM%I%I%XCPL_EFLOOD(0))
    ALLOCATE(IM%I%I%XCPL_PFLOOD(0))
    ALLOCATE(IM%I%I%XCPL_IFLOOD(0))     
  ENDIF     
!
ELSE
!
  ALLOCATE(IM%I%I%XCPL_RUNOFF  (0))
  ALLOCATE(IM%I%I%XCPL_DRAIN   (0))
  ALLOCATE(IM%I%I%XCPL_ICEFLUX (0))
  ALLOCATE(IM%I%I%XCPL_RECHARGE(0))
  ALLOCATE(IM%I%I%XCPL_EFLOOD  (0))
  ALLOCATE(IM%I%I%XCPL_PFLOOD  (0))
  ALLOCATE(IM%I%I%XCPL_IFLOOD  (0))
!
ENDIF
!
IF(IM%I%O%LWTD.AND..NOT.IM%I%O%LGW)THEN
  WRITE(ILUOUT,*)'COMPUTE_ISBA_PARAMETERS: Groundwater map is required by SFX - Groundwater coupling '
  WRITE(ILUOUT,*)'COMPUTE_ISBA_PARAMETERS: Please check your pgd namelist where this map must be     '
  WRITE(ILUOUT,*)'COMPUTE_ISBA_PARAMETERS: specified (YGW and YGWFILETYPE, or XUNIF_GW, or LIMP_GW)  '
  CALL ABOR1_SFX('COMPUTE_ISBA_PARAMETERS: Groundwater map is required by SFX - Groundwater coupling')
ENDIF
!
! * Initialize flood scheme :
!
IF(IM%I%O%LFLOOD)THEN
  ALLOCATE(IM%I%I%XFFLOOD (KI))
  ALLOCATE(IM%I%I%XPIFLOOD(KI))
  ALLOCATE(IM%I%I%XFF     (KI,IM%I%O%NPATCH))
  ALLOCATE(IM%I%I%XFFG    (KI,IM%I%O%NPATCH))
  ALLOCATE(IM%I%I%XFFV    (KI,IM%I%O%NPATCH))
  ALLOCATE(IM%I%I%XFFROZEN(KI,IM%I%O%NPATCH))
  ALLOCATE(IM%I%I%XALBF   (KI,IM%I%O%NPATCH))
  ALLOCATE(IM%I%I%XEMISF  (KI,IM%I%O%NPATCH)) 
  IM%I%I%XFFLOOD       = 0.0
  IM%I%I%XPIFLOOD      = 0.0
  IM%I%I%XFF           = 0.0
  IM%I%I%XFFG          = 0.0
  IM%I%I%XFFV          = 0.0
  IM%I%I%XFFROZEN      = 0.0
  IM%I%I%XALBF         = 0.0
  IM%I%I%XEMISF        = 0.0  
ELSE
  ALLOCATE(IM%I%I%XFFLOOD   (0))
  ALLOCATE(IM%I%I%XPIFLOOD  (0))
  ALLOCATE(IM%I%I%XFF     (0,0))
  ALLOCATE(IM%I%I%XFFG    (0,0))
  ALLOCATE(IM%I%I%XFFV    (0,0))
  ALLOCATE(IM%I%I%XFFROZEN(0,0))
  ALLOCATE(IM%I%I%XALBF   (0,0))
  ALLOCATE(IM%I%I%XEMISF  (0,0))
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      7.     ISBA time-varying deep force-restore temperature initialization
!              ---------------------------------------------------------------
!
 CALL SOILTEMP_ARP_PAR(IM%I%O%CISBA, IM%I%O%XSODELX, &
                       HPROGRAM,IM%I%O%LTEMP_ARP,IM%I%O%NTEMPLAYER_ARP)
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
IF (CASSIM_ISBA=="ENKF ") THEN
  !
  CALL INIT_RANDOM_SEED()
  !
ENDIF
!
CALL INIT_IO_SURF_n(DTCO, DGU, U, &
                    HPROGRAM,'NATURE','ISBA  ','READ ')
!
!*      10.     Prognostic and semi-prognostic fields
!               -------------------------------------
!
 CALL READ_ISBA_n(DTCO, IM%I, U, &
                  HPROGRAM)
!
IF (HINIT/='ALL') THEN
  CALL END_IO_SURF_n(HPROGRAM)
  IF (LHOOK) CALL DR_HOOK('COMPUTE_ISBA_PARAMETERS',1,ZHOOK_HANDLE)
  RETURN
END IF
!
IF (HINIT=='PRE' .AND. IM%I%R%TSNOW%SCHEME.NE.'3-L' .AND. &
        IM%I%R%TSNOW%SCHEME.NE.'CRO' .AND. IM%I%O%CISBA=='DIF') THEN
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
   CALL INIT_ISBA_LANDUSE(DTCO, IM%IG, IM%I, UG, U, &
                          HPROGRAM)  
END IF
!
!-------------------------------------------------------------------------------
!
!*      12.     Canopy air fields:
!               -----------------
!
 CALL READ_ISBA_CANOPY_n(DTCO, IM%ICP, IM%I%O%LCANOPY, U, &
                         HPROGRAM)
!
!-------------------------------------------------------------------------------
!
!*      13.     initialize radiative and physical properties
!               --------------------------------------------
!
ALLOCATE(IM%I%I%XDIR_ALB_WITH_SNOW(KI,KSW,IM%I%O%NPATCH))
ALLOCATE(IM%I%I%XSCA_ALB_WITH_SNOW(KI,KSW,IM%I%O%NPATCH))
IM%I%I%XDIR_ALB_WITH_SNOW = 0.0
IM%I%I%XSCA_ALB_WITH_SNOW = 0.0
!
 CALL INIT_VEG_n(IM%I%O%NPATCH, KI, IM%I%O%LCANOPY, IM%I%O%CROUGH, IM%I%O%LAGRI_TO_GRASS, IM%I%R%TSNOW, &
                 IM%I%O%CPHOTO, IM%I%M%T%XLAIMIN, IM%I%M%X%XH_TREE, IM%I%IP%XVEGTYPE_PATCH, IM%I%M%T%XLAI, &
                 IM%I%M%T%XZ0, IM%I%M%T%XVEG, IM%I%M%T%XEMIS, &
                 IM%I%O%LTR_ML, IM%I%R%XFAPARC, IM%I%R%XFAPIRC, IM%I%R%XLAI_EFFC, IM%I%R%XMUS, &
                 IM%I%M%A%XALBNIR_SOIL, IM%I%M%A%XALBVIS_SOIL, IM%I%M%A%XALBUV_SOIL, IM%I%M%T%XALBNIR, &
                 IM%I%M%T%XALBVIS, IM%I%M%T%XALBUV, &
                 IM%DGMI%LSURF_DIAG_ALBEDO, IM%I%R%XPSN, IM%I%R%XPSNG, IM%I%R%XPSNV, IM%I%R%XPSNV_A, &
                 PDIR_ALB, PSCA_ALB, PEMIS, PTSRAD )
!
DO JPATCH=1,IM%I%O%NPATCH
  ZWG1(:,JPATCH) = IM%I%R%XWG(:,1,JPATCH)
  ZTG1(:,JPATCH) = IM%I%R%XTG(:,1,JPATCH)
END DO
!
 CALL CONVERT_PATCH_ISBA(DTCO, IM%DTI, IM%I%O, &
                         IM%I%O%CISBA,IDECADE,IDECADE2,IM%I%P%XCOVER,IM%I%P%LCOVER,&
                          IM%I%O%CPHOTO,LAGRIP,IM%I%O%LPERM,IM%I%O%LTR_ML,'NAT',   &
                          PWG1=ZWG1, PWSAT=IM%I%IP%XWSAT,        &
                          PALBNIR_SOIL=IM%I%M%A%XALBNIR_SOIL, PALBVIS_SOIL=IM%I%M%A%XALBVIS_SOIL, &
                          PALBUV_SOIL=IM%I%M%A%XALBUV_SOIL ,&
                          PALBVIS_DRY=IM%I%IP%XALBVIS_DRY, PALBNIR_DRY=IM%I%IP%XALBNIR_DRY, &
                          PALBUV_DRY=IM%I%IP%XALBUV_DRY,  PALBVIS_WET=IM%I%IP%XALBVIS_WET, &
                          PALBNIR_WET=IM%I%IP%XALBNIR_WET, PALBUV_WET=IM%I%IP%XALBUV_WET)
!
! Load randomly perturbed fields. Perturbation ratios are saved in case fields are reset later.
IF(IM%I%O%LPERTSURF) THEN
!
  CALL READ_SURF(&
                 HPROGRAM,'VEG',IM%I%M%T%XVEG(:,:),IRESP)
  ALLOCATE(IM%I%I%XPERTVEG(KI))
  IM%I%I%XPERTVEG(:)=IM%I%M%T%XVEG(:,1)
!
  CALL READ_SURF(&
                 HPROGRAM,'LAI',IM%I%M%T%XLAI(:,:),IRESP)
  ALLOCATE(IM%I%I%XPERTLAI(KI))
  IM%I%I%XPERTLAI(:)=IM%I%M%T%XLAI(:,1)
!
  CALL READ_SURF(&
                 HPROGRAM,'CV',IM%I%M%T%XCV(:,:),IRESP)
  ALLOCATE(IM%I%I%XPERTCV(KI))
  IM%I%I%XPERTCV(:)=IM%I%M%T%XCV(:,1)
!
  CALL READ_SURF(&
                 HPROGRAM,'PERTALB',ZPERTBUF(:,:),IRESP)
  ALLOCATE(IM%I%I%XPERTALB(KI))
  IM%I%I%XPERTALB(:)=ZPERTBUF(:,1)
  WHERE(IM%I%M%T%XALBNIR_VEG(:,1)/=XUNDEF)  IM%I%M%T%XALBNIR_VEG(:,1) = IM%I%M%T%XALBNIR_VEG(:,1) *( 1.+ IM%I%I%XPERTALB(:) )
  WHERE(IM%I%M%T%XALBVIS_VEG(:,1)/=XUNDEF)  IM%I%M%T%XALBVIS_VEG(:,1) = IM%I%M%T%XALBVIS_VEG(:,1) *( 1.+ IM%I%I%XPERTALB(:) )
  WHERE(IM%I%M%T%XALBUV_VEG(:,1)/=XUNDEF)   IM%I%M%T%XALBUV_VEG(:,1)  = IM%I%M%T%XALBUV_VEG(:,1)  *( 1.+ IM%I%I%XPERTALB(:) )
  WHERE(IM%I%M%A%XALBNIR_SOIL(:,1)/=XUNDEF) IM%I%M%A%XALBNIR_SOIL(:,1)= IM%I%M%A%XALBNIR_SOIL(:,1)*( 1.+ IM%I%I%XPERTALB(:) )
  WHERE(IM%I%M%A%XALBVIS_SOIL(:,1)/=XUNDEF) IM%I%M%A%XALBVIS_SOIL(:,1)= IM%I%M%A%XALBVIS_SOIL(:,1)*( 1.+ IM%I%I%XPERTALB(:) )
  WHERE(IM%I%M%A%XALBUV_SOIL(:,1)/=XUNDEF)  IM%I%M%A%XALBUV_SOIL(:,1) = IM%I%M%A%XALBUV_SOIL(:,1) *( 1.+ IM%I%I%XPERTALB(:) )
!
  CALL READ_SURF(&
                 HPROGRAM,'PERTZ0LAND',ZPERTBUF(:,:),IRESP)
  ALLOCATE(IM%I%I%XPERTZ0(KI))
  IM%I%I%XPERTZ0(:)=ZPERTBUF(:,1)
  WHERE(IM%I%M%T%XZ0(:,1)/=XUNDEF)      IM%I%M%T%XZ0(:,1)     =IM%I%M%T%XZ0(:,1)     *( 1.+ IM%I%I%XPERTZ0(:) )
  WHERE(IM%I%IP%XZ0EFFIP(:,1)/=XUNDEF) IM%I%IP%XZ0EFFIP(:,1)=IM%I%IP%XZ0EFFIP(:,1)*( 1.+ IM%I%I%XPERTZ0(:) )
  WHERE(IM%I%IP%XZ0EFFIM(:,1)/=XUNDEF) IM%I%IP%XZ0EFFIM(:,1)=IM%I%IP%XZ0EFFIM(:,1)*( 1.+ IM%I%I%XPERTZ0(:) )
  WHERE(IM%I%IP%XZ0EFFJP(:,1)/=XUNDEF) IM%I%IP%XZ0EFFJP(:,1)=IM%I%IP%XZ0EFFJP(:,1)*( 1.+ IM%I%I%XPERTZ0(:) )
  WHERE(IM%I%IP%XZ0EFFJM(:,1)/=XUNDEF) IM%I%IP%XZ0EFFJM(:,1)=IM%I%IP%XZ0EFFJM(:,1)*( 1.+ IM%I%I%XPERTZ0(:) )
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       14.    Output radiative fields
!               -----------------------
!
ALLOCATE(IM%I%I%XEMIS_NAT   (KI))
IM%I%I%XEMIS_NAT (:) = XUNDEF
!
 CALL AVERAGED_ALBEDO_EMIS_ISBA(IM%I, &
                                IM%I%O%LFLOOD, IM%I%O%CALBEDO, PZENITH,                &
                                 IM%I%M%T%XVEG,IM%I%M%T%XZ0,IM%I%M%T%XLAI,                          &
                                 IM%I%O%LMEB_PATCH,IM%I%M%M%XGNDLITTER,&
                                 IM%I%M%M%XZ0LITTER,IM%I%M%M%XLAIGV, &
                                 IM%I%M%M%XH_VEG, IM%I%R%XTV,               &
                                 ZTG1,                                   &
                                 IM%I%IP%XPATCH,                                 &
                                 PSW_BANDS,                              &
                                 IM%I%M%T%XALBNIR_VEG,IM%I%M%T%XALBVIS_VEG,IM%I%M%T%XALBUV_VEG,     &
                                 IM%I%M%A%XALBNIR_SOIL,IM%I%M%A%XALBVIS_SOIL,IM%I%M%A%XALBUV_SOIL,  &
                                 IM%I%M%T%XEMIS,                                  &
                                 IM%I%R%TSNOW,                                  &
                                 IM%I%M%T%XALBNIR,IM%I%M%T%XALBVIS,IM%I%M%T%XALBUV,                 &
                                 PDIR_ALB, PSCA_ALB,                     &
                                 IM%I%I%XEMIS_NAT,ZTSRAD_NAT,ZTSURF_NAT         )  
!
PEMIS  = IM%I%I%XEMIS_NAT
PTSRAD = ZTSRAD_NAT
PTSURF = ZTSURF_NAT
!
!-------------------------------------------------------------------------------
!
!*      15.     ISBA diagnostics initialization
!               -------------------------------
!
IF(IM%I%O%NPATCH<=1) IM%DGI%LPATCH_BUDGET=.FALSE.
!
 CALL DIAG_ISBA_INIT_n(&
                       IM%CHI, IM%DGEI, IM%DGI, IM%DGMI, DGU, IM%GB, IM%I%O, &
                       IM%I%R%TSNOW%SCHEME, IM%I%R%TSNOW%NLAYER, SIZE(IM%I%IP%XABC), &
                       HPROGRAM,KI,KSW)
!
!-------------------------------------------------------------------------------
!
 CALL INIT_SURF_TOPD(IM%DGEI, IM%I, UG, U, &
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


