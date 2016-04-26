!#############################################################
SUBROUTINE COMPUTE_ISBA_PARAMETERS (DTCO, OREAD_BUDGETC, UG, U, IM, NDST, SLT, SV,   &
                                    HPROGRAM,HINIT,OLAND_USE,                  &
                                    KI,KSV,KSW,HSV,PCO2,PRHOA,                 &
                                    PZENITH,PSW_BANDS,PDIR_ALB,PSCA_ALB,       &
                                    PEMIS,PTSRAD,PTSURF, HTEST             )  
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
USE MODD_ISBA_n, ONLY : ISBA_P_t, ISBA_PE_t, ISBA_K_t
USE MODD_AGRI_n, ONLY : AGRI_t
USE MODD_SSO_n, ONLY : SSO_t
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_DST_n, ONLY : DST_NP_t, DST_t
USE MODD_SLT_n, ONLY : SLT_t
USE MODD_SV_n, ONLY : SV_t
USE MODD_GRID_n, ONLY : GRID_t
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
USE MODI_GET_Z0REL
USE MODI_SET_ROUGH
USE MODI_GET_LUOUT
USE MODI_ABOR1_SFX
USE MODI_INIT_IO_SURF_n
USE MODI_ALLOCATE_PHYSIO
USE MODI_INIT_ISBA_MIXPAR
USE MODI_CONVERT_PATCH_ISBA
USE MODI_CONVERT_PATCH_ISBA_1P
USE MODI_INIT_VEG_PGD_n
USE MODI_INIT_TOP
USE MODI_EXP_DECAY_SOIL_FR
USE MODI_CARBON_INIT
USE MODI_SOILTEMP_ARP_PAR
USE MODI_END_IO_SURF_n
!
USE MODI_MAKE_CHOICE_ARRAY
USE MODI_READ_SURF
USE MODI_READ_ISBA_n
USE MODI_INIT_ISBA_LANDUSE
USE MODI_READ_SBL_n
USE MODI_INIT_VEG_n
USE MODI_INIT_CHEMICAL_n
USE MODI_OPEN_NAMELIST
USE MODI_CH_INIT_DEP_ISBA_n
USE MODI_CLOSE_NAMELIST
USE MODI_INIT_DST
USE MODI_INIT_SLT
USE MODI_AVERAGED_ALBEDO_EMIS_ISBA
USE MODI_DIAG_ISBA_INIT_n
USE MODI_INIT_SURF_TOPD
USE MODI_ISBA_SOC_PARAMETERS
USE MODI_PACK_SAME_RANK
!
USE MODI_READ_AND_SEND_MPI
USE MODI_ISBA_TO_TOPD
USE MODI_OPEN_FILE
USE MODI_CLOSE_FILE
USE MODI_FIX_MEB_VEG
USE MODI_AV_PGD
USE MODI_SURF_PATCH
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
TYPE(ISBA_MODEL_t), INTENT(INOUT) :: IM
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
LOGICAL, INTENT(IN) :: OREAD_BUDGETC
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(DST_NP_t), INTENT(INOUT) :: NDST
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
TYPE(GRID_t), POINTER :: GK
TYPE(ISBA_P_t), POINTER :: PK
TYPE(ISBA_K_t), POINTER :: KK
TYPE(ISBA_PE_t), POINTER :: PEK
TYPE(AGRI_t), POINTER :: AGK
TYPE(SSO_t), POINTER :: ISSK
TYPE(DST_t), POINTER :: DSTK
!
REAL, DIMENSION(U%NDIM_FULL)   :: ZF_PARAM, ZC_DEPTH_RATIO
!
REAL, DIMENSION(KI)     :: ZTSRAD_NAT !radiative temperature
REAL, DIMENSION(KI)     :: ZTSURF_NAT !effective temperature
!
REAL, DIMENSION(U%NSIZE_NATURE)  :: ZWG1 ! work array for surface water content
REAL, DIMENSION(U%NSIZE_NATURE,IM%O%NPATCH)  :: ZTG1 ! work array for surface temperature
!
REAL, DIMENSION(KI)         :: ZM
REAL, DIMENSION(:,:), ALLOCATABLE :: ZWORK
REAL, DIMENSION(KI,IM%O%NPATCH)  :: ZF
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZDG_SOIL, ZDG_SOIL_P
REAL, DIMENSION(:), ALLOCATABLE :: ZSUM_PATCH
!
INTEGER :: ICH     ! unit of input chemistry file
INTEGER           :: JI, JL     ! loop increment
INTEGER           :: ILUOUT   ! unit of output listing file
INTEGER           :: IRESP   ! return code
INTEGER           :: IDECADE, IDECADE2  ! decade of simulation
INTEGER           :: JP  ! loop counter on tiles
INTEGER           :: ISIZE_LMEB_PATCH  ! Number of patches with MEB=true
!
LOGICAL :: GDIM
INTEGER :: JVEG, IVERSION, IBUGFIX, IMASK, ISIZE, JMAXLOC
!
 CHARACTER(LEN=4)  :: YLVL
 CHARACTER(LEN=12) :: YRECFM
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
ALLOCATE(IM%S%XVEGTYPE(KI,NVEGTYPE))
IF (IM%DTI%LDATA_VEGTYPE) THEN
  IM%S%XVEGTYPE = IM%DTI%XPAR_VEGTYPE
ELSE
  !classical ecoclimap case
  DO JVEG=1,NVEGTYPE
    CALL AV_PGD(DTCO, IM%S%XVEGTYPE(:,JVEG),IM%S%XCOVER ,DTCO%XDATA_VEGTYPE(:,JVEG),'NAT','ARI',IM%S%LCOVER)
  END DO
ENDIF
!
ALLOCATE(IM%S%XPATCH(KI,IM%O%NPATCH))
ALLOCATE(IM%S%XVEGTYPE_PATCH(KI,NVEGTYPE,IM%O%NPATCH))
 CALL SURF_PATCH(IM%O%NPATCH,IM%S%XVEGTYPE,IM%S%XPATCH,IM%S%XVEGTYPE_PATCH)
DO JP = 1,IM%O%NPATCH
  IM%NP%AL(JP)%NSIZE_P = COUNT(IM%S%XPATCH(:,JP) > 0.0)
ENDDO
!
IF (IM%O%XRM_PATCH/=0.) THEN
  !
  WRITE(ILUOUT,*) " REMOVE PATCH below 5 % add to dominant patch " 
  ! remove small fraction of PATCHES and add to MAIN PATCH
  DO JI = 1,KI
    !1) find most present patch maximum value 
    JMAXLOC = MAXVAL(MAXLOC(IM%S%XPATCH(JI,:)))
    !2) FIND small value of cover 
    DO JP = 1,IM%O%NPATCH
      IF ( IM%S%XPATCH(JI,JP)<IM%O%XRM_PATCH ) THEN
        IM%S%XPATCH(JI,JMAXLOC) = IM%S%XPATCH(JI,JMAXLOC) + IM%S%XPATCH(JI,JP)
        IM%S%XPATCH(JI,JP) = 0.0
       ENDIF
    ENDDO
  ENDDO
  !
ENDIF
!
!*       2.3    Physiographic data fields from land cover:
!               -----------------------------------------
!
IF (IM%S%TTIME%TDATE%MONTH /= NUNDEF) THEN
  IDECADE = 3 * ( IM%S%TTIME%TDATE%MONTH - 1 ) + MIN(IM%S%TTIME%TDATE%DAY-1,29) / 10 + 1
ELSE
  IDECADE = 1
END IF
!
IDECADE2 = IDECADE
!
 CALL INIT_ISBA_MIXPAR(DTCO, IM%DTI, IM%G%NDIM, IM%O, IDECADE,IDECADE2,IM%S%XCOVER,IM%S%LCOVER,'NAT')
!
ISIZE_LMEB_PATCH=COUNT(IM%O%LMEB_PATCH(:))
IF (ISIZE_LMEB_PATCH>0)  THEN
  CALL FIX_MEB_VEG(IM%DTI, IM%G%NDIM, IM%O%LMEB_PATCH, IM%O%NPATCH)
ENDIF
!
!
DO JP = 1, IM%O%NPATCH
  !
  ISIZE = IM%NP%AL(JP)%NSIZE_P
  PK => IM%NP%AL(JP)
  KK => IM%NK%AL(JP)
  PEK => IM%NPE%AL(JP)
  AGK => IM%NAG%AL(JP)
  ISSK => IM%NISS%AL(JP)
  !
  ALLOCATE(PK%NR_P    (ISIZE))
  CALL GET_1D_MASK(PK%NSIZE_P,KI,IM%S%XPATCH(:,JP),PK%NR_P) 
  !
  ALLOCATE(KK%XVEGTYPE(ISIZE,NVEGTYPE))
  DO JI = 1,PK%NSIZE_P
    KK%XVEGTYPE(JI,:) = IM%S%XVEGTYPE(PK%NR_P(JI),:)
  ENDDO
  !
  ALLOCATE(PK%XPATCH(ISIZE))
  ALLOCATE(PK%XVEGTYPE_PATCH (ISIZE,NVEGTYPE))
  !
  CALL SURF_PATCH(JP,IM%O%NPATCH,KK%XVEGTYPE,PK%XPATCH,PK%XVEGTYPE_PATCH)
  !
  !*       2.5    Masks for tiles
  !               ---------------
  !
  CALL ALLOCATE_PHYSIO(IM%O, KK, PK, PEK, NVEGTYPE  )
  !
  IF (IM%O%LPERM) THEN
    ALLOCATE(KK%XPERM(ISIZE))
    CALL PACK_SAME_RANK(PK%NR_P, IM%K%XPERM, KK%XPERM)
  ELSE
    ALLOCATE(KK%XPERM(0))
  ENDIF
  !
  CALL CONVERT_PATCH_ISBA_1P(DTCO, IM%DTI, IM%O, IDECADE, IDECADE2, IM%S%XCOVER, IM%S%LCOVER, &
                         LAGRIP, 'NAT', JP, KK, PK, PEK, .TRUE., .TRUE., .TRUE., .TRUE., .FALSE., .FALSE., &
                         PSOILGRID=IM%O%XSOILGRID, PPERM=KK%XPERM  )
  !
  !-------------------------------------------------------------------------------
  !
  ALLOCATE(KK%XSAND(ISIZE,IM%O%NGROUND_LAYER))
  ALLOCATE(KK%XCLAY(ISIZE,IM%O%NGROUND_LAYER))
  !
  ALLOCATE(ISSK%XAOSIP(ISIZE))
  ALLOCATE(ISSK%XAOSIM(ISIZE))
  ALLOCATE(ISSK%XAOSJP(ISIZE))
  ALLOCATE(ISSK%XAOSJM(ISIZE))
  ALLOCATE(ISSK%XHO2IP(ISIZE))
  ALLOCATE(ISSK%XHO2IM(ISIZE))
  ALLOCATE(ISSK%XHO2JP(ISIZE))
  ALLOCATE(ISSK%XHO2JM(ISIZE))
  !
  CALL PACK_SAME_RANK(PK%NR_P, IM%K%XSAND, KK%XSAND)
  CALL PACK_SAME_RANK(PK%NR_P, IM%K%XCLAY, KK%XCLAY)
  !
  CALL PACK_SAME_RANK(PK%NR_P,IM%ISS%XAOSIP,ISSK%XAOSIP)
  CALL PACK_SAME_RANK(PK%NR_P,IM%ISS%XAOSIM,ISSK%XAOSIM)
  CALL PACK_SAME_RANK(PK%NR_P,IM%ISS%XAOSJP,ISSK%XAOSJP)
  CALL PACK_SAME_RANK(PK%NR_P,IM%ISS%XAOSJM,ISSK%XAOSJM)
  CALL PACK_SAME_RANK(PK%NR_P,IM%ISS%XHO2IP,ISSK%XHO2IP)
  CALL PACK_SAME_RANK(PK%NR_P,IM%ISS%XHO2IM,ISSK%XHO2IM)
  CALL PACK_SAME_RANK(PK%NR_P,IM%ISS%XHO2JP,ISSK%XHO2JP)
  CALL PACK_SAME_RANK(PK%NR_P,IM%ISS%XHO2JM,ISSK%XHO2JM)
  !
  CALL INIT_VEG_PGD_n(ISSK, IM%O, IM%S, IM%K, KK, PK, PEK, AGK, KI,  &
                      HPROGRAM, 'NATURE', ILUOUT, ISIZE, IM%S%TTIME%TDATE%MONTH, &
                      LDEEPSOIL, LPHYSDOMC, XTDEEP_CLI, XGAMMAT_CLI,   & 
                      LAGRIP, XTHRESHOLD, HINIT, PCO2, PRHOA        )  
  !
  !Rainfall spatial distribution
  !CRAIN used in HYDRO_VEG and HYDRO_SGH and ISBA_SGH_UPDATE
  IF(IM%O%CRAIN=='SGH')THEN
    ALLOCATE(KK%XMUF(ISIZE))
    KK%XMUF(:)=0.0
  ELSE
    ALLOCATE(KK%XMUF(0))
  ENDIF
  !
  ALLOCATE(KK%XFSAT(ISIZE))  
  KK%XFSAT(:) = 0.0
  !
  ! * Initialize flood scheme :
  !
  ALLOCATE(KK%XFFLOOD (ISIZE))
  ALLOCATE(KK%XPIFLOOD(ISIZE))
  ALLOCATE(KK%XFF     (ISIZE))
  ALLOCATE(KK%XFFG    (ISIZE))
  ALLOCATE(KK%XFFV    (ISIZE))
  ALLOCATE(KK%XFFROZEN(ISIZE))
  ALLOCATE(KK%XALBF   (ISIZE))
  ALLOCATE(KK%XEMISF  (ISIZE)) 
  KK%XFFLOOD       = 0.0
  KK%XPIFLOOD      = 0.0
  KK%XFF           = 0.0
  KK%XFFG          = 0.0
  KK%XFFV          = 0.0
  KK%XFFROZEN      = 0.0
  KK%XALBF         = 0.0
  KK%XEMISF        = 0.0  
  !
ENDDO
!
IF(IM%O%CRAIN=='SGH')THEN
  ALLOCATE(IM%K%XMUF(ISIZE))
  IM%K%XMUF(:)=0.0
ENDIF

ALLOCATE(IM%ISS%XZ0REL(KI))
 CALL GET_Z0REL(IM%ISS)
!
IM%ISS%XAOSIP => NULL()
IM%ISS%XAOSIM => NULL()
IM%ISS%XAOSJP => NULL()
IM%ISS%XAOSJM => NULL()
IM%ISS%XHO2IP => NULL()
IM%ISS%XHO2IM => NULL()
IM%ISS%XHO2JP => NULL()
IM%ISS%XHO2JM => NULL()
!
!-------------------------------------------------------------------------------
!
!        3.  Initialize Chemical Deposition
!            ------------------------------
!
!        3.1 Chemical gazes
!            --------------
!
    !* for the time being, chemistry on vegetation works only for
    ! ISBA on nature tile (not for gardens), because subroutine INIT_CHEMICAL_n
    ! contains explicitely modules from ISBAn. It should be cleaned in a future
    ! version.
 CALL INIT_CHEMICAL_n(ILUOUT, KSV, HSV, IM%CHI%SVI, IM%CHI%CCH_NAMES, IM%CHI%CAER_NAMES,  &
                     HDSTNAMES=IM%CHI%CDSTNAMES, HSLTNAMES=IM%CHI%CSLTNAMES        )
!
IF (KSV /= 0) THEN
  !
  IF (IM%CHI%SVI%NBEQ > 0) THEN
    !* for the time being, chemistry deposition on vegetation works only for
    ! ISBA on nature tile (not for gardens), because subroutine CH_INIT_DEP_ISBA_n
    ! contains explicitely modules from ISBAn. It should be cleaned in a future
    ! version.
    CALL OPEN_NAMELIST(HPROGRAM, ICH, HFILE=IM%CHI%CCHEM_SURF_FILE)
    CALL CH_INIT_DEP_ISBA_n(IM%CHI, IM%NCHI, IM%NP, DTCO, IM%O%NPATCH, &
                            IM%S%LCOVER, IM%S%XCOVER, ICH, ILUOUT, KI)
    CALL CLOSE_NAMELIST(HPROGRAM, ICH)
  END IF
  !
  DO JP = 1,IM%O%NPATCH
    !
    DSTK => NDST%AL(JP)
    !
    IF (IM%CHI%SVI%NDSTEQ >=1) THEN
      !
      ALLOCATE (DSTK%XSFDST (ISIZE, IM%CHI%SVI%NDSTEQ))  !Output array
      ALLOCATE (DSTK%XSFDSTM(ISIZE, IM%CHI%SVI%NDSTEQ))  !Output array
      DSTK%XSFDST(:,:)  = 0.
      DSTK%XSFDSTM(:,:) = 0.     
      CALL INIT_DST(DSTK, U, HPROGRAM, PK%NSIZE_P, PK%NR_P, PK%XVEGTYPE_PATCH)    
    ELSE
      ALLOCATE(DSTK%XSFDST (0,0))
      ALLOCATE(DSTK%XSFDSTM(0,0))
    END IF
    !
  ENDDO
  !
  IF (IM%CHI%SVI%NSLTEQ >=1) THEN
    CALL INIT_SLT(SLT, HPROGRAM)
  END IF
  !
ENDIF
!
!-------------------------------------------------------------------------------
!
!DIF option :
!    Anisotropy coeficient for hydraulic conductivity for topmodel drainage (Fan et al. 2006)
!    Soil organic matter effect and/or Exponential decay for DIF option
!    Must be call before INIT_TOP
!
!
IF(IM%O%CISBA=='DIF') THEN
  !
  IF( IM%O%CKSAT=='SGH' )THEN
    WRITE(ILUOUT,*)'THE KSAT EXP PROFILE WITH ISBA-DF IS NOT PHYSIC AND HAS BEEN REMOVED FOR NOW' 
    WRITE(ILUOUT,*)'A NEW PHYSICAL APPROACH WILL BE DEVELLOPED ACCOUNTING FOR COMPACTION IN ALL ' 
    WRITE(ILUOUT,*)'HYDRODYNAMIC PARAMETERS (WSAT, PSISAT, KSAT, B) AND NOT ONLY IN KSAT        ' 
    CALL ABOR1_SFX('CKSAT=SGH is not physic with ISBA-DF and has been removed for now')    
  ENDIF
  !  
  IF(IM%O%LSOC)THEN   
    IF(.NOT.IM%O%LSOCP)THEN
      WRITE(ILUOUT,*)'LSOC = T can be activated only if SOC data given in PGD fields'
      CALL ABOR1_SFX('LSOC = T can be activated only if SOC data given in PGD fields')
    ENDIF
    !
    ALLOCATE(IM%S%XFRACSOC(KI,IM%O%NGROUND_LAYER))
      CALL ISBA_SOC_PARAMETERS(IM%O%CRUNOFF, IM%S%XSOC, IM%K, IM%NP, IM%S%XFRACSOC, &
                               IM%K%XWSAT, IM%K%XWFC, IM%K%XWWILT, IM%O%NPATCH )
    !
  ELSE
    ALLOCATE(IM%S%XFRACSOC(0,0))
  ENDIF
  !
ELSE
  ALLOCATE(IM%S%XFRACSOC(0,0))
ENDIF
!
!Topmodel
!
!
ZF   (:,:) = XUNDEF
ZM   (:)   = XUNDEF
!
!CRUNOFF used in hydro_sgh and isba_sgh_update
IF( IM%O%CRUNOFF=='SGH ') THEN 
!
  ALLOCATE(IM%S%XTAB_FSAT(KI,NDIMTAB))
  ALLOCATE(IM%S%XTAB_WTOP(KI,NDIMTAB))
  ALLOCATE(IM%S%XTAB_QTOP(KI,NDIMTAB))
!
  IM%S%XTAB_FSAT(:,:) = 0.0
  IM%S%XTAB_WTOP(:,:) = 0.0
  IM%S%XTAB_QTOP(:,:) = 0.0
!
  IF(HINIT/='PRE' .AND. .NOT.LASSIM)THEN
!
    WHERE(IM%K%XCLAY(:,1)==XUNDEF.AND.IM%S%XTI_MEAN(:)/=XUNDEF) IM%S%XTI_MEAN(:)=XUNDEF
!
    CALL INIT_TOP(IM%O, IM%S, IM%K, IM%NK, IM%NP, ILUOUT, ZM  )  
!
    IF (IM%O%CKSAT=='SGH' .AND. IM%O%CISBA/='DIF') THEN
!     Exponential decay factor calculate using soil properties 
!     (eq. 11, Decharme et al., J. Hydrometeor, 2006)
      DO JP = 1,IM%O%NPATCH
        !
        KK => IM%NK%AL(JP)
        PK => IM%NP%AL(JP)
        !
        DO JI = 1,PK%NSIZE_P
          IMASK = PK%NR_P(JI)
          !
          IF (ZM(IMASK)/=XUNDEF) ZF(JI,JP) = (IM%K%XWSAT(IMASK,1)-IM%K%XWD0(IMASK,1))/ZM(IMASK)
        ENDDO
        !
      ENDDO
!       
    ENDIF
!
  ENDIF
!
! Subsurface flow by layer (m/s)
  DO JP = 1,IM%O%NPATCH
    PK => IM%NP%AL(JP)
    IF(IM%O%CISBA=='DIF') THEN
      ALLOCATE(PK%XTOPQS(PK%NSIZE_P,IM%O%NGROUND_LAYER))
      PK%XTOPQS(:,:)=0.0
    ELSE
      ALLOCATE(PK%XTOPQS(0,0))
    ENDIF
  ENDDO
!
ELSE                  
  !  
  ALLOCATE(IM%S%XTAB_FSAT(0,0))
  ALLOCATE(IM%S%XTAB_WTOP(0,0))
  ALLOCATE(IM%S%XTAB_QTOP(0,0))
  DO JP = 1,IM%O%NPATCH
    PK => IM%NP%AL(JP)
    ALLOCATE(PK%XTOPQS(0,0))
  ENDDO  
  !                  
ENDIF  
!
!Exponential decay for ISBA-FR option
!CKSAT used in hydro_soil.F90 and soil.F90
IF(HINIT/='PRE'.AND.IM%O%CISBA/='DIF')THEN 
  !
  IF(IM%O%CKSAT=='SGH') THEN
    !
    DO JP = 1,IM%O%NPATCH
      PK => IM%NP%AL(JP)

      WHERE(ZF(:,JP)==XUNDEF.AND.PK%XDG(:,2)/=XUNDEF) 
        ZF(:,JP) = 4.0/PK%XDG(:,2)
      ENDWHERE
      ZF(:,JP) = MIN(ZF(:,JP),XF_DECAY)
      !
      CALL EXP_DECAY_SOIL_FR(IM%O%CISBA, ZF(1:PK%NSIZE_P,JP), PK )
      !
    ENDDO
    !
    ALLOCATE(IM%S%XF_PARAM (KI))
    IM%S%XF_PARAM(:) = XUNDEF
    DO JI = 1,IM%NP%AL(1)%NSIZE_P
      IM%S%XF_PARAM(IM%NP%AL(1)%NR_P(JI)) = ZF(JI,1)
    ENDDO       
    !
  ELSEIF ( IM%O%CKSAT=='EXP' .AND. IM%O%CISBA=='3-L' ) THEN
    !
    ALLOCATE(IM%S%XF_PARAM (KI))
    IM%S%XF_PARAM(:) = XUNDEF
    !
    IF (HPROGRAM/='AROME ' .AND. HPROGRAM/='MESONH ') THEN
      !
      CALL OPEN_FILE('ASCII ',NUNIT,HFILE='carte_f_dc.txt',HFORM='FORMATTED',HACTION='READ ')
      DO JI=1,U%NDIM_FULL
        READ(NUNIT,*) ZF_PARAM(JI), ZC_DEPTH_RATIO(JI)
      ENDDO
      CALL CLOSE_FILE('ASCII ',NUNIT)
      CALL READ_AND_SEND_MPI(ZF_PARAM,IM%S%XF_PARAM,U%NR_NATURE)
#ifdef TOPD
      IF (.NOT.ALLOCATED(XC_DEPTH_RATIO)) ALLOCATE(XC_DEPTH_RATIO (KI))
      XC_DEPTH_RATIO(:) = XUNDEF
      CALL READ_AND_SEND_MPI(ZC_DEPTH_RATIO,XC_DEPTH_RATIO,U%NR_NATURE)
#endif
      !
    ELSE
      WRITE(ILUOUT,*) "COMPUTE_ISBA_PARAMETERS: WITH CKSAT=EXP, IN NOT OFFLINE "//&
                      "MODE, TOPMODEL FILE FOR F_PARAM IS NOT READ "
    ENDIF
    !
    DO JP=1,IM%O%NPATCH
      !
      PK => IM%NP%AL(JP)
      !
      DO JI = 1,PK%NSIZE_P
        IMASK = PK%NR_P(JI)
!
        IF (IM%S%XF_PARAM(IMASK)==XUNDEF.AND.PK%XDG(JI,2)/=XUNDEF) THEN
          ZF(JI,JP) = 4.0/PK%XDG(JI,2)
        ELSE
          ZF(JI,JP) = IM%S%XF_PARAM(IMASK)
        ENDIF
      ENDDO
      !
      ZF(:,JP) = MIN(ZF(:,JP),XF_DECAY)
      !
      IF (ALLOCATED(XC_DEPTH_RATIO)) &
        CALL PACK_SAME_RANK(PK%NR_P,XC_DEPTH_RATIO,ZC_DEPTH_RATIO(1:PK%NSIZE_P))
      CALL EXP_DECAY_SOIL_FR(IM%O%CISBA, ZF(1:PK%NSIZE_P,JP), PK, ZC_DEPTH_RATIO(1:PK%NSIZE_P))  
      !
    ENDDO
    !    
  ENDIF
!  ! 
ENDIF
!
!
!*       2.10   Soil carbon
!               -----------                        
!
IF (HINIT == 'ALL' .AND. IM%O%CRESPSL=='CNT' .AND. IM%O%CPHOTO == 'NCB') THEN
  CALL CARBON_INIT(IM%O)
ENDIF
!
!
!-------------------------------------------------------------------------------
!
!*       6.2    Initialize of SFX - RRM coupling:
!               ---------------------------------
!
! * Check some key :
!
IF(LCPL_CALVING)THEN
   IF(.NOT.IM%O%LGLACIER)THEN
     CALL ABOR1_SFX('COMPUTE_ISBA_PARAMETERS: LGLACIER MUST BE ACTIVATED IF LCPL_CALVING')
   ENDIF
ENDIF
!
! * Initialize required coupling fields :
!
IM%O%LCPL_RRM = .FALSE.
IM%O%LFLOOD   = .FALSE.
IM%O%LWTD     = .FALSE.
!
IF(LCPL_LAND)THEN
!    
  IM%O%LCPL_RRM = .TRUE.
!
  ALLOCATE(IM%S%XCPL_DRAIN (KI))
  ALLOCATE(IM%S%XCPL_RUNOFF(KI))
  IM%S%XCPL_DRAIN (:) = 0.0
  IM%S%XCPL_RUNOFF(:) = 0.0
!
  IF(IM%O%LGLACIER)THEN
     ALLOCATE(IM%S%XCPL_ICEFLUX(KI))
     IM%S%XCPL_ICEFLUX(:) = 0.0
  ELSE
     ALLOCATE(IM%S%XCPL_ICEFLUX(0))
  ENDIF
!
  IF(LCPL_GW)THEN
    IM%O%LWTD = .TRUE.
    ALLOCATE(IM%S%XCPL_RECHARGE(KI))
    IM%S%XCPL_RECHARGE(:) = 0.0
  ELSE
    ALLOCATE(IM%S%XCPL_RECHARGE(0))
  ENDIF
!
  IF(LCPL_FLOOD)THEN
     IM%O%LFLOOD = .TRUE.
     ALLOCATE(IM%S%XCPL_EFLOOD(KI))
     ALLOCATE(IM%S%XCPL_PFLOOD(KI))
     ALLOCATE(IM%S%XCPL_IFLOOD(KI))
     IM%S%XCPL_EFLOOD(:)= 0.0
     IM%S%XCPL_PFLOOD(:)= 0.0
     IM%S%XCPL_IFLOOD(:)= 0.0    
  ELSE
    ALLOCATE(IM%S%XCPL_EFLOOD(0))
    ALLOCATE(IM%S%XCPL_PFLOOD(0))
    ALLOCATE(IM%S%XCPL_IFLOOD(0))     
  ENDIF     
!
ELSE
!
  ALLOCATE(IM%S%XCPL_RUNOFF  (0))
  ALLOCATE(IM%S%XCPL_DRAIN   (0))
  ALLOCATE(IM%S%XCPL_ICEFLUX (0))
  ALLOCATE(IM%S%XCPL_RECHARGE(0))
  ALLOCATE(IM%S%XCPL_EFLOOD  (0))
  ALLOCATE(IM%S%XCPL_PFLOOD  (0))
  ALLOCATE(IM%S%XCPL_IFLOOD  (0))
!
ENDIF
!
IF(IM%O%LWTD.AND..NOT.IM%O%LGW)THEN
  WRITE(ILUOUT,*)'COMPUTE_ISBA_PARAMETERS: Groundwater map is required by SFX - Groundwater coupling '
  WRITE(ILUOUT,*)'COMPUTE_ISBA_PARAMETERS: Please check your pgd namelist where this map must be     '
  WRITE(ILUOUT,*)'COMPUTE_ISBA_PARAMETERS: specified (YGW and YGWFILETYPE, or XUNIF_GW, or LIMP_GW)  '
  CALL ABOR1_SFX('COMPUTE_ISBA_PARAMETERS: Groundwater map is required by SFX - Groundwater coupling')
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      7.     ISBA time-varying deep force-restore temperature initialization
!              ---------------------------------------------------------------
!
 CALL SOILTEMP_ARP_PAR(IM%O, HPROGRAM)
!
!-------------------------------------------------------------------------------
!
DO JP = 1,IM%O%NPATCH
  !
  KK => IM%NK%AL(JP)
  PK => IM%NP%AL(JP)
  ISSK => IM%NISS%AL(JP)
  GK => IM%NG%AL(JP)
  !
  ALLOCATE(KK%XMPOTSAT(PK%NSIZE_P,IM%O%NGROUND_LAYER))
  ALLOCATE(KK%XBCOEF  (PK%NSIZE_P,IM%O%NGROUND_LAYER))
  ALLOCATE(KK%XWWILT  (PK%NSIZE_P,IM%O%NGROUND_LAYER))
  ALLOCATE(KK%XWFC    (PK%NSIZE_P,IM%O%NGROUND_LAYER))
  ALLOCATE(KK%XWSAT   (PK%NSIZE_P,IM%O%NGROUND_LAYER))
  ALLOCATE(KK%XWD0    (PK%NSIZE_P,IM%O%NGROUND_LAYER))
  ALLOCATE(KK%XKANISO (PK%NSIZE_P,IM%O%NGROUND_LAYER))
  !
  CALL PACK_SAME_RANK(PK%NR_P,IM%K%XMPOTSAT,KK%XMPOTSAT)
  CALL PACK_SAME_RANK(PK%NR_P,IM%K%XBCOEF,KK%XBCOEF)
  !
  CALL PACK_SAME_RANK(PK%NR_P,IM%K%XWWILT,KK%XWWILT)
  CALL PACK_SAME_RANK(PK%NR_P,IM%K%XWFC,KK%XWFC)
  CALL PACK_SAME_RANK(PK%NR_P,IM%K%XWSAT,KK%XWSAT)
  !   
  IF (IM%O%CRUNOFF=='SGH') THEN
    CALL PACK_SAME_RANK(PK%NR_P,IM%K%XWD0,KK%XWD0)
    CALL PACK_SAME_RANK(PK%NR_P,IM%K%XKANISO,KK%XKANISO)
  ENDIF
  !  
  IF (IM%O%CISBA=='2-L' .OR. IM%O%CISBA=='3-L') THEN
    ALLOCATE(KK%XC4B  (PK%NSIZE_P))
    ALLOCATE(KK%XACOEF(PK%NSIZE_P))
    ALLOCATE(KK%XPCOEF(PK%NSIZE_P))
    ALLOCATE(KK%XCGSAT(PK%NSIZE_P))
    CALL PACK_SAME_RANK(PK%NR_P,IM%K%XC4B,  KK%XC4B)
    CALL PACK_SAME_RANK(PK%NR_P,IM%K%XACOEF,KK%XACOEF)
    CALL PACK_SAME_RANK(PK%NR_P,IM%K%XPCOEF,KK%XPCOEF)
    CALL PACK_SAME_RANK(PK%NR_P,IM%K%XCGSAT,KK%XCGSAT)
  ENDIF
  !
  IF (IM%O%CSCOND=='PL98'.OR.IM%O%CISBA=='DIF') THEN
    ALLOCATE(KK%XHCAPSOIL(PK%NSIZE_P,IM%O%NGROUND_LAYER))
    ALLOCATE(KK%XCONDDRY (PK%NSIZE_P,IM%O%NGROUND_LAYER))
    ALLOCATE(KK%XCONDSLD (PK%NSIZE_P,IM%O%NGROUND_LAYER))
    CALL PACK_SAME_RANK(PK%NR_P,IM%K%XHCAPSOIL,KK%XHCAPSOIL)
    CALL PACK_SAME_RANK(PK%NR_P,IM%K%XCONDDRY ,KK%XCONDDRY)
    CALL PACK_SAME_RANK(PK%NR_P,IM%K%XCONDSLD ,KK%XCONDSLD)
  ENDIF
  !
  ALLOCATE(KK%XWDRAIN (PK%NSIZE_P))
  ALLOCATE(KK%XRUNOFFB(PK%NSIZE_P))
  CALL PACK_SAME_RANK(PK%NR_P,IM%K%XWDRAIN,KK%XWDRAIN)
  CALL PACK_SAME_RANK(PK%NR_P,IM%K%XRUNOFFB,KK%XRUNOFFB)
  !
  ALLOCATE(ISSK%XZ0REL(PK%NSIZE_P))
  ALLOCATE(ISSK%XSSO_SLOPE(PK%NSIZE_P))
  !
  CALL PACK_SAME_RANK(PK%NR_P,IM%ISS%XZ0REL,ISSK%XZ0REL)
  CALL PACK_SAME_RANK(PK%NR_P,IM%ISS%XSSO_SLOPE,ISSK%XSSO_SLOPE)
  !
  ALLOCATE(GK%XLAT(PK%NSIZE_P))
  ALLOCATE(GK%XLON(PK%NSIZE_P))
  !
  CALL PACK_SAME_RANK(PK%NR_P,IM%G%XLAT,GK%XLAT)
  CALL PACK_SAME_RANK(PK%NR_P,IM%G%XLON,GK%XLON)
  !
ENDDO
!
IM%K%XKANISO => NULL()
!
IM%K%XMPOTSAT => NULL()
IM%K%XBCOEF => NULL()
!
IM%K%XC4B => NULL()
IM%K%XACOEF => NULL()
IM%K%XPCOEF => NULL()  
IM%K%XCGSAT => NULL()
!
IM%K%XHCAPSOIL => NULL()
IM%K%XCONDDRY => NULL()
IM%K%XCONDSLD => NULL()
!
IM%K%XWDRAIN => NULL()
IM%K%XRUNOFFB => NULL()
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
CALL INIT_IO_SURF_n(DTCO, U, HPROGRAM,'NATURE','ISBA  ','READ ')
!
!*      10.     Prognostic and semi-prognostic fields
!               -------------------------------------
!
 CALL READ_ISBA_n(DTCO, IM%O, IM%S, IM%NP, IM%NPE, IM%K%XCLAY, U, HPROGRAM)
!
IF (HINIT/='ALL') THEN
  CALL END_IO_SURF_n(HPROGRAM)
  IF (LHOOK) CALL DR_HOOK('COMPUTE_ISBA_PARAMETERS',1,ZHOOK_HANDLE)
  RETURN
END IF
!
IF (HINIT=='PRE' .AND. IM%NPE%AL(1)%TSNOW%SCHEME.NE.'3-L' .AND. &
        IM%NPE%AL(1)%TSNOW%SCHEME.NE.'CRO' .AND. IM%O%CISBA=='DIF') THEN
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
  !
  CALL READ_SURF(HPROGRAM,'VERSION',IVERSION,IRESP)
  CALL READ_SURF(HPROGRAM,'BUG',IBUGFIX,IRESP)
  GDIM = (IVERSION>8 .OR. IVERSION==8 .AND. IBUGFIX>0)
  !  
  ALLOCATE(ZWORK(KI,IM%O%NPATCH))
  !
  !* read old patch fraction
  !
  DO JP = 1,IM%O%NPATCH
    ALLOCATE(IM%NP%AL(JP)%XPATCH_OLD(IM%NP%AL(JP)%NSIZE_P))
  ENDDO
  !
  CALL MAKE_CHOICE_ARRAY(HPROGRAM, IM%O%NPATCH, GDIM, 'PATCH', ZWORK)
  DO JP = 1,IM%O%NPATCH
    CALL PACK_SAME_RANK(IM%NP%AL(JP)%NR_P,ZWORK(:,JP),IM%NP%AL(JP)%XPATCH_OLD(:))
  ENDDO
  !
  !* read old soil layer thicknesses (m)
  !
  DO JP = 1,IM%O%NPATCH
    ALLOCATE(IM%NP%AL(JP)%XDG_OLD(IM%NP%AL(JP)%NSIZE_P,IM%O%NGROUND_LAYER))
  ENDDO
  !
  DO JL=1,IM%O%NGROUND_LAYER
    WRITE(YLVL,'(I4)') JL
    YRECFM='OLD_DG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    CALL MAKE_CHOICE_ARRAY(HPROGRAM, IM%O%NPATCH, GDIM, YRECFM, ZWORK)
    DO JP = 1,IM%O%NPATCH
      CALL PACK_SAME_RANK(IM%NP%AL(JP)%NR_P,ZWORK(:,JP),IM%NP%AL(JP)%XDG_OLD(:,JL))
    ENDDO
  END DO
  DEALLOCATE(ZWORK)
  !
   CALL INIT_ISBA_LANDUSE(DTCO, UG, U, IM%O, IM%NK, IM%NP, IM%NPE, IM%G%XMESH_SIZE, &
                          HPROGRAM)  
END IF
!
!-------------------------------------------------------------------------------
!
!*      12.     Canopy air fields:
!               -----------------
!
 CALL READ_SBL_n(DTCO, U, IM%SB, IM%O%LCANOPY, HPROGRAM, "NATURE")
!
!-------------------------------------------------------------------------------
!
!*      13.     initialize radiative and physical properties
!               --------------------------------------------
!
DO JP=1,IM%O%NPATCH
  PK => IM%NP%AL(JP)
  KK => IM%NK%AL(JP)
  ALLOCATE(KK%XDIR_ALB_WITH_SNOW(PK%NSIZE_P,KSW))
  ALLOCATE(KK%XSCA_ALB_WITH_SNOW(PK%NSIZE_P,KSW))
  KK%XDIR_ALB_WITH_SNOW = 0.0
  KK%XSCA_ALB_WITH_SNOW = 0.0
ENDDO
!
 CALL SET_ROUGH(IM%O%LCANOPY,IM%O%CROUGH)
!
DO JP = 1,IM%O%NPATCH
  KK => IM%NK%AL(JP)
  PK => IM%NP%AL(JP)
  PEK => IM%NPE%AL(JP)

  CALL INIT_VEG_n(IM%O, KK, PK, PEK, IM%ID%DM%LSURF_DIAG_ALBEDO, PDIR_ALB, PSCA_ALB, PEMIS, PTSRAD )
  !
  ZWG1(1:PK%NSIZE_P) = PEK%XWG(:,1)
  ZTG1(1:PK%NSIZE_P,JP) = PEK%XTG(:,1)
  !
  CALL CONVERT_PATCH_ISBA_1P(DTCO, IM%DTI, IM%O, IDECADE, IDECADE2, IM%S%XCOVER, IM%S%LCOVER,&
                         LAGRIP, 'NAT', JP, KK, PK, PEK, &
                         .FALSE., .FALSE., .FALSE., .FALSE., .TRUE., .FALSE., &
                         PWG1=ZWG1(1:PK%NSIZE_P), PWSAT=KK%XWSAT)
  !
ENDDO
!
! Load randomly perturbed fields. Perturbation ratios are saved in case fields are reset later.
IF(IM%O%LPERTSURF) THEN
  !
  CALL READ_SURF(HPROGRAM,'VERSION',IVERSION,IRESP)
  CALL READ_SURF(HPROGRAM,'BUG',IBUGFIX,IRESP)
  GDIM = (IVERSION>8 .OR. IVERSION==8 .AND. IBUGFIX>0)
  !
  ALLOCATE(ZWORK(KI,IM%O%NPATCH))
  !
  CALL MAKE_CHOICE_ARRAY(HPROGRAM, IM%O%NPATCH, GDIM, 'VEG', ZWORK)
  ALLOCATE(IM%S%XPERTVEG(KI))
  IM%S%XPERTVEG(:)=ZWORK(:,1)
!
  CALL MAKE_CHOICE_ARRAY(HPROGRAM, IM%O%NPATCH, GDIM, 'LAI', ZWORK)
  ALLOCATE(IM%S%XPERTLAI(KI))
  IM%S%XPERTLAI(:)=ZWORK(:,1)
!
  CALL MAKE_CHOICE_ARRAY(HPROGRAM, IM%O%NPATCH, GDIM, 'CV', ZWORK)
  ALLOCATE(IM%S%XPERTCV(KI))
  IM%S%XPERTCV(:)=ZWORK(:,1)
!
  CALL MAKE_CHOICE_ARRAY(HPROGRAM, IM%O%NPATCH, GDIM, 'PERTALB', ZWORK)
  ALLOCATE(IM%S%XPERTALB(KI))
  IM%S%XPERTALB(:)=ZWORK(:,1)

  PEK => IM%NPE%AL(1)
  ISSK => IM%NISS%AL(1)

  WHERE(PEK%XALBNIR_VEG(:)/=XUNDEF)  PEK%XALBNIR_VEG(:) = PEK%XALBNIR_VEG(:) *( 1.+ IM%S%XPERTALB(:) )
  WHERE(PEK%XALBVIS_VEG(:)/=XUNDEF)  PEK%XALBVIS_VEG(:) = PEK%XALBVIS_VEG(:) *( 1.+ IM%S%XPERTALB(:) )
  WHERE(PEK%XALBUV_VEG (:)/=XUNDEF)  PEK%XALBUV_VEG (:)  = PEK%XALBUV_VEG(:) *( 1.+ IM%S%XPERTALB(:) )
  WHERE(PEK%XALBNIR_SOIL(:)/=XUNDEF) PEK%XALBNIR_SOIL(:)= PEK%XALBNIR_SOIL(:) *( 1.+ IM%S%XPERTALB(:) )
  WHERE(PEK%XALBVIS_SOIL(:)/=XUNDEF) PEK%XALBVIS_SOIL(:)= PEK%XALBVIS_SOIL(:) *( 1.+ IM%S%XPERTALB(:) )
  WHERE(PEK%XALBUV_SOIL (:)/=XUNDEF) PEK%XALBUV_SOIL (:) = PEK%XALBUV_SOIL(:) *( 1.+ IM%S%XPERTALB(:) )
!
  CALL MAKE_CHOICE_ARRAY(HPROGRAM, IM%O%NPATCH, GDIM, 'PERTZ0LAND', ZWORK)
  ALLOCATE(IM%S%XPERTZ0(KI))
  IM%S%XPERTZ0(:)=ZWORK(:,1)
  WHERE(PEK%XZ0(:)/=XUNDEF)      PEK%XZ0(:)       =PEK%XZ0(:)      *( 1.+ IM%S%XPERTZ0(:) )
  WHERE(ISSK%XZ0EFFIP(:)/=XUNDEF) ISSK%XZ0EFFIP(:)=ISSK%XZ0EFFIP(:)*( 1.+ IM%S%XPERTZ0(:) )
  WHERE(ISSK%XZ0EFFIM(:)/=XUNDEF) ISSK%XZ0EFFIM(:)=ISSK%XZ0EFFIM(:)*( 1.+ IM%S%XPERTZ0(:) )
  WHERE(ISSK%XZ0EFFJP(:)/=XUNDEF) ISSK%XZ0EFFJP(:)=ISSK%XZ0EFFJP(:)*( 1.+ IM%S%XPERTZ0(:) )
  WHERE(ISSK%XZ0EFFJM(:)/=XUNDEF) ISSK%XZ0EFFJM(:)=ISSK%XZ0EFFJM(:)*( 1.+ IM%S%XPERTZ0(:) )
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       14.    Output radiative fields
!               -----------------------
!
ALLOCATE(IM%S%XEMIS_NAT   (KI))
IM%S%XEMIS_NAT (:) = XUNDEF
!
 CALL AVERAGED_ALBEDO_EMIS_ISBA(IM%O, IM%S, IM%NK, IM%NP, IM%NPE, &
                                PZENITH, ZTG1, PSW_BANDS, PDIR_ALB, PSCA_ALB, &
                                IM%S%XEMIS_NAT,ZTSRAD_NAT,ZTSURF_NAT         )  
!
PEMIS  = IM%S%XEMIS_NAT
PTSRAD = ZTSRAD_NAT
PTSURF = ZTSURF_NAT
!
!-------------------------------------------------------------------------------
!
!*      15.     ISBA diagnostics initialization
!               -------------------------------
!
IF(IM%O%NPATCH<=1) IM%ID%O%LPATCH_BUDGET=.FALSE.
!
 CALL DIAG_ISBA_INIT_n(IM%CHI, IM%ID%DE, IM%ID%DEC, IM%ID%DEP, IM%ID%DEPC, IM%ID%O, &
                       IM%ID%D, IM%ID%DC, IM%ID%DP, IM%ID%DPC, IM%ID%DM, IM%ID%DMP, &
                       OREAD_BUDGETC, IM%NGB, IM%GB, IM%O, IM%NP, IM%NPE%AL(1)%TSNOW%SCHEME, &
                       IM%NPE%AL(1)%TSNOW%NLAYER, SIZE(IM%S%XABC), HPROGRAM,KI,KSW)
!
!-------------------------------------------------------------------------------
!
 CALL INIT_SURF_TOPD(IM%ID%DEC, IM%O, IM%S, IM%K, IM%NP, IM%NPE, UG, U, HPROGRAM, U%NDIM_FULL)
!
!-------------------------------------------------------------------------------
!
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
!
IF (LHOOK) CALL DR_HOOK('COMPUTE_ISBA_PARAMETERS',1,ZHOOK_HANDLE)
!
END SUBROUTINE COMPUTE_ISBA_PARAMETERS


