!#############################################################
SUBROUTINE COMPUTE_ISBA_PARAMETERS (DTCO, OREAD_BUDGETC, UG, U, AG, CHI, DTI,  &
                                    DGI, GB, ICP, ISS, IG, I, DST, SLT, SV,    &
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
USE MODD_SURFEX_n, ONLY : ISBA_DIAG_t
!
USE MODD_AGRI_n, ONLY : AGRI_t
USE MODD_CH_ISBA_n, ONLY : CH_ISBA_t
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_GR_BIOG_n, ONLY : GR_BIOG_t
USE MODD_ISBA_n, ONLY : ISBA_t
!
USE MODD_GRID_n, ONLY : GRID_t
USE MODD_CANOPY_n, ONLY : CANOPY_t
USE MODD_SSO_n, ONLY : SSO_t
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
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
TYPE(AGRI_t), INTENT(INOUT) :: AG
TYPE(CH_ISBA_t), INTENT(INOUT) :: CHI
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTI
TYPE(ISBA_DIAG_t), INTENT(INOUT) :: DGI
TYPE(GR_BIOG_t), INTENT(INOUT) :: GB
TYPE(CANOPY_t), INTENT(INOUT) :: ICP
TYPE(SSO_t), INTENT(INOUT) :: ISS 
TYPE(GRID_t), INTENT(INOUT) :: IG
TYPE(ISBA_t), INTENT(INOUT) :: I
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
LOGICAL, INTENT(IN) :: OREAD_BUDGETC
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
REAL, DIMENSION(KI,I%O%NPATCH)  :: ZWG1 ! work array for surface water content
REAL, DIMENSION(KI,I%O%NPATCH)  :: ZTG1 ! work array for surface temperature
!
REAL, DIMENSION(KI)         :: ZM, ZWORK
REAL, DIMENSION(KI,I%O%NPATCH)  :: ZF, ZPERTBUF
!
INTEGER :: ICH     ! unit of input chemistry file
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
 CALL ALLOCATE_PHYSIO(I%O, I%M, KI, NVEGTYPE  )
!
IF (I%I%TTIME%TDATE%MONTH /= NUNDEF) THEN
  IDECADE = 3 * ( I%I%TTIME%TDATE%MONTH - 1 ) + MIN(I%I%TTIME%TDATE%DAY-1,29) / 10 + 1
ELSE
  IDECADE = 1
END IF
!
IDECADE2 = IDECADE
!
 CALL INIT_ISBA_MIXPAR(DTCO, DTI, IG%NDIM, I%O, IDECADE,IDECADE2,I%P%XCOVER,I%P%LCOVER,'NAT')
!
ISIZE_LMEB_PATCH=COUNT(I%O%LMEB_PATCH(:))
IF (ISIZE_LMEB_PATCH>0)  THEN
  CALL FIX_MEB_VEG(DTI, IG%NDIM, I%O%LMEB_PATCH, I%O%NPATCH)
ENDIF
!
 CALL CONVERT_PATCH_ISBA(DTCO, DTI, I%O, IDECADE, IDECADE2, I%P%XCOVER, I%P%LCOVER, &
                        LAGRIP, 'NAT', .FALSE., IMX=I%M%X, IMT=I%M%T, IMM=I%M%M, IMI=I%M%I, &
                        PSOILGRID=I%O%XSOILGRID, PPERM=I%P%XPERM  )
!
!-------------------------------------------------------------------------------
!
CALL INIT_VEG_PGD_n(ISS, DTCO, U, I%O, I%IP, I%P, I%M%X, I%M%T, AG, &
                    HPROGRAM, 'NATURE', ILUOUT, KI, I%I%TTIME%TDATE%MONTH, &
                    LDEEPSOIL, LPHYSDOMC, XTDEEP_CLI, XGAMMAT_CLI,   & 
                    LAGRIP, XTHRESHOLD, HINIT, PCO2, PRHOA        )  
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
 CALL INIT_CHEMICAL_n(ILUOUT, KSV, HSV, CHI%SVI, CHI%CCH_NAMES, CHI%CAER_NAMES,  &
                     HDSTNAMES=CHI%CDSTNAMES, HSLTNAMES=CHI%CSLTNAMES        )
!
IF (KSV /= 0) THEN
  !
  IF (CHI%SVI%NBEQ > 0) THEN
    !* for the time being, chemistry deposition on vegetation works only for
    ! ISBA on nature tile (not for gardens), because subroutine CH_INIT_DEP_ISBA_n
    ! contains explicitely modules from ISBAn. It should be cleaned in a future
    ! version.
    CALL OPEN_NAMELIST(HPROGRAM, ICH, HFILE=CHI%CCHEM_SURF_FILE)
    CALL CH_INIT_DEP_ISBA_n(CHI, DTCO, I%O%NPATCH, I%P%LCOVER, I%P%XCOVER, ICH, ILUOUT, KI)
    CALL CLOSE_NAMELIST(HPROGRAM, ICH)
  END IF
  !
  IF (CHI%SVI%NDSTEQ >=1) THEN
    ALLOCATE (DST%XSFDST (KI, CHI%SVI%NDSTEQ, I%O%NPATCH))  !Output array
    ALLOCATE (DST%XSFDSTM(KI, CHI%SVI%NDSTEQ, I%O%NPATCH))  !Output array
    DST%XSFDST(:,:,:)  = 0.
    DST%XSFDSTM(:,:,:) = 0.     
    CALL INIT_DST(DST, U, HPROGRAM,I%IP%NSIZE_NATURE_P,I%IP%NR_NATURE_P, I%O%NPATCH,I%IP%XVEGTYPE_PATCH)    
  ELSE
    ALLOCATE(DST%XSFDST (0,0,0))
    ALLOCATE(DST%XSFDSTM(0,0,0))
  END IF
  !
  IF (CHI%SVI%NSLTEQ >=1) THEN
    ALLOCATE (SLT%XSFSLT(KI,CHI%SVI%NSLTEQ,I%O%NPATCH))  !Output array
    CALL INIT_SLT(SLT, HPROGRAM)   
  ELSE
    ALLOCATE(SLT%XSFSLT(0,0,0))
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
IF(I%O%CISBA=='DIF') THEN
  !
  IF( I%O%CKSAT=='SGH' )THEN
    WRITE(ILUOUT,*)'THE KSAT EXP PROFILE WITH ISBA-DF IS NOT PHYSIC AND HAS BEEN REMOVED FOR NOW' 
    WRITE(ILUOUT,*)'A NEW PHYSICAL APPROACH WILL BE DEVELLOPED ACCOUNTING FOR COMPACTION IN ALL ' 
    WRITE(ILUOUT,*)'HYDRODYNAMIC PARAMETERS (WSAT, PSISAT, KSAT, B) AND NOT ONLY IN KSAT        ' 
    CALL ABOR1_SFX('CKSAT=SGH is not physic with ISBA-DF and has been removed for now')    
  ENDIF
  !  
  IF(I%O%LSOC)THEN   
    IF(.NOT.I%O%LSOCP)THEN
      WRITE(ILUOUT,*)'LSOC = T can be activated only if SOC data given in PGD fields'
      CALL ABOR1_SFX('LSOC = T can be activated only if SOC data given in PGD fields')
    ENDIF
    ALLOCATE(I%I%XFRACSOC(KI,I%O%NGROUND_LAYER))
    I%I%XFRACSOC(:,:)=0.0
    CALL ISBA_SOC_PARAMETERS(I%O%CRUNOFF,I%M%X%XDG,I%P%XSOC,I%IP,I%I%XFRACSOC )
  ELSE
    ALLOCATE(I%I%XFRACSOC(0,0))
  ENDIF
!
ELSE
  ALLOCATE(I%I%XFRACSOC(0,0))
ENDIF
!
!Topmodel
!
!CRUNOFF used in hydro_sgh and isba_sgh_update
IF( I%O%CRUNOFF=='SGH ') THEN 
!
  ALLOCATE(I%I%XTAB_FSAT(KI,NDIMTAB))
  ALLOCATE(I%I%XTAB_WTOP(KI,NDIMTAB))
  ALLOCATE(I%I%XTAB_QTOP(KI,NDIMTAB))
!
  I%I%XTAB_FSAT(:,:) = 0.0
  I%I%XTAB_WTOP(:,:) = 0.0
  I%I%XTAB_QTOP(:,:) = 0.0
!
  IF(HINIT/='PRE' .AND. .NOT.LASSIM)THEN
!
    WHERE(I%P%XCLAY(:,1)==XUNDEF.AND.I%P%XTI_MEAN(:)/=XUNDEF) I%P%XTI_MEAN(:)=XUNDEF
!
    CALL INIT_TOP(I%O, I%IP, I%P, I%I, ILUOUT, ZM           )  
!
!
    IF (I%O%CKSAT=='SGH' .AND. I%O%CISBA/='DIF') THEN
!     Exponential decay factor calculate using soil properties 
!     (eq. 11, Decharme et al., J. Hydrometeor, 2006)
      DO JILU=1,KI
        IF (ZM(JILU)/=XUNDEF) ZF(JILU,:) = (I%IP%XWSAT(JILU,1)-I%IP%XWD0(JILU,1))/ZM(JILU)
      ENDDO
!       
    ENDIF
!
  ENDIF
!
! Subsurface flow by layer (m/s)
  IF(I%O%CISBA=='DIF') THEN
    ALLOCATE(I%I%XTOPQS(KI,I%O%NGROUND_LAYER,I%O%NPATCH))
    I%I%XTOPQS(:,:,:)=0.0
  ELSE
    ALLOCATE(I%I%XTOPQS(0,0,0))
  ENDIF
!
ELSE                  
!  
  ALLOCATE(I%I%XTAB_FSAT(0,0))
  ALLOCATE(I%I%XTAB_WTOP(0,0))
  ALLOCATE(I%I%XTAB_QTOP(0,0))
  ALLOCATE(I%I%XTOPQS(0,0,0))  
!                  
ENDIF  
!
!Exponential decay for ISBA-FR option
!CKSAT used in hydro_soil.F90 and soil.F90
IF(HINIT/='PRE'.AND.I%O%CISBA/='DIF')THEN 
  !
  IF(I%O%CKSAT=='SGH') THEN
    !
    WHERE(ZF(:,:)==XUNDEF.AND.I%M%X%XDG(:,2,:)/=XUNDEF) 
      ZF(:,:) = 4.0/I%M%X%XDG(:,2,:)
    ENDWHERE
    ZF(:,:) = MIN(ZF(:,:),XF_DECAY)
    !
    ALLOCATE(I%I%XF_PARAM (KI))
    I%I%XF_PARAM(:) = ZF(:,1)
    !
    CALL EXP_DECAY_SOIL_FR(I%O%CISBA, ZF,I%IP, I%M%X )                 
    !
  ELSEIF ( I%O%CKSAT=='EXP' .AND. I%O%CISBA=='3-L' ) THEN
    !
    ALLOCATE(I%I%XF_PARAM (KI))
    I%I%XF_PARAM(:) = XUNDEF
    !
    IF (HPROGRAM/='AROME ' .AND. HPROGRAM/='MESONH ') THEN
      !
      CALL OPEN_FILE('ASCII ',NUNIT,HFILE='carte_f_dc.txt',HFORM='FORMATTED',HACTION='READ ')
      DO JILU=1,U%NDIM_FULL
        READ(NUNIT,*) ZF_PARAM(JILU), ZC_DEPTH_RATIO(JILU)
      ENDDO
      CALL CLOSE_FILE('ASCII ',NUNIT)
      CALL READ_AND_SEND_MPI(ZF_PARAM,I%I%XF_PARAM,U%NR_NATURE)
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
    DO JPATCH=1,I%O%NPATCH
      WHERE (I%I%XF_PARAM(:)==XUNDEF.AND.I%M%X%XDG(:,2,JPATCH)/=XUNDEF)
        ZF(:,JPATCH) = 4.0/I%M%X%XDG(:,2,JPATCH)
      ELSEWHERE
        ZF(:,JPATCH) = I%I%XF_PARAM(:)
      ENDWHERE
    ENDDO
     ZF(:,:) = MIN(ZF(:,:),XF_DECAY)
    !
    CALL EXP_DECAY_SOIL_FR(I%O%CISBA, ZF,I%IP, I%M%X)  
    !    
  ENDIF
  ! 
ENDIF
!
!
!*       2.10   Soil carbon
!               -----------                        
!
IF (HINIT == 'ALL' .AND. I%O%CRESPSL=='CNT' .AND. I%O%CPHOTO == 'NCB') THEN
  CALL CARBON_INIT(I%O)
ENDIF
!
!Rainfall spatial distribution
!CRAIN used in HYDRO_VEG and HYDRO_SGH and ISBA_SGH_UPDATE
IF(I%O%CRAIN=='SGH')THEN
  ALLOCATE(I%I%XMUF(KI))
  I%I%XMUF(:)=0.0
ELSE
  ALLOCATE(I%I%XMUF(0))
ENDIF
!
ALLOCATE(I%I%XFSAT(KI))  
I%I%XFSAT(:) = 0.0
!
!-------------------------------------------------------------------------------
!
!*       6.2    Initialize of SFX - RRM coupling:
!               ---------------------------------
!
! * Check some key :
!
IF(LCPL_CALVING)THEN
   IF(.NOT.I%O%LGLACIER)THEN
     CALL ABOR1_SFX('COMPUTE_ISBA_PARAMETERS: LGLACIER MUST BE ACTIVATED IF LCPL_CALVING')
   ENDIF
ENDIF
!
! * Initialize required coupling fields :
!
I%O%LCPL_RRM = .FALSE.
I%O%LFLOOD   = .FALSE.
I%O%LWTD     = .FALSE.
!
IF(LCPL_LAND)THEN
!    
  I%O%LCPL_RRM = .TRUE.
!
  ALLOCATE(I%I%XCPL_DRAIN (KI))
  ALLOCATE(I%I%XCPL_RUNOFF(KI))
  I%I%XCPL_DRAIN (:) = 0.0
  I%I%XCPL_RUNOFF(:) = 0.0
!
  IF(I%O%LGLACIER)THEN
     ALLOCATE(I%I%XCPL_ICEFLUX(KI))
     I%I%XCPL_ICEFLUX(:) = 0.0
  ELSE
     ALLOCATE(I%I%XCPL_ICEFLUX(0))
  ENDIF
!
  IF(LCPL_GW)THEN
    I%O%LWTD = .TRUE.
    ALLOCATE(I%I%XCPL_RECHARGE(KI))
    I%I%XCPL_RECHARGE(:) = 0.0
  ELSE
    ALLOCATE(I%I%XCPL_RECHARGE(0))
  ENDIF
!
  IF(LCPL_FLOOD)THEN
     I%O%LFLOOD = .TRUE.
     ALLOCATE(I%I%XCPL_EFLOOD(KI))
     ALLOCATE(I%I%XCPL_PFLOOD(KI))
     ALLOCATE(I%I%XCPL_IFLOOD(KI))
     I%I%XCPL_EFLOOD(:)= 0.0
     I%I%XCPL_PFLOOD(:)= 0.0
     I%I%XCPL_IFLOOD(:)= 0.0    
  ELSE
    ALLOCATE(I%I%XCPL_EFLOOD(0))
    ALLOCATE(I%I%XCPL_PFLOOD(0))
    ALLOCATE(I%I%XCPL_IFLOOD(0))     
  ENDIF     
!
ELSE
!
  ALLOCATE(I%I%XCPL_RUNOFF  (0))
  ALLOCATE(I%I%XCPL_DRAIN   (0))
  ALLOCATE(I%I%XCPL_ICEFLUX (0))
  ALLOCATE(I%I%XCPL_RECHARGE(0))
  ALLOCATE(I%I%XCPL_EFLOOD  (0))
  ALLOCATE(I%I%XCPL_PFLOOD  (0))
  ALLOCATE(I%I%XCPL_IFLOOD  (0))
!
ENDIF
!
IF(I%O%LWTD.AND..NOT.I%O%LGW)THEN
  WRITE(ILUOUT,*)'COMPUTE_ISBA_PARAMETERS: Groundwater map is required by SFX - Groundwater coupling '
  WRITE(ILUOUT,*)'COMPUTE_ISBA_PARAMETERS: Please check your pgd namelist where this map must be     '
  WRITE(ILUOUT,*)'COMPUTE_ISBA_PARAMETERS: specified (YGW and YGWFILETYPE, or XUNIF_GW, or LIMP_GW)  '
  CALL ABOR1_SFX('COMPUTE_ISBA_PARAMETERS: Groundwater map is required by SFX - Groundwater coupling')
ENDIF
!
! * Initialize flood scheme :
!
IF(I%O%LFLOOD)THEN
  ALLOCATE(I%I%XFFLOOD (KI))
  ALLOCATE(I%I%XPIFLOOD(KI))
  ALLOCATE(I%I%XFF     (KI,I%O%NPATCH))
  ALLOCATE(I%I%XFFG    (KI,I%O%NPATCH))
  ALLOCATE(I%I%XFFV    (KI,I%O%NPATCH))
  ALLOCATE(I%I%XFFROZEN(KI,I%O%NPATCH))
  ALLOCATE(I%I%XALBF   (KI,I%O%NPATCH))
  ALLOCATE(I%I%XEMISF  (KI,I%O%NPATCH)) 
  I%I%XFFLOOD       = 0.0
  I%I%XPIFLOOD      = 0.0
  I%I%XFF           = 0.0
  I%I%XFFG          = 0.0
  I%I%XFFV          = 0.0
  I%I%XFFROZEN      = 0.0
  I%I%XALBF         = 0.0
  I%I%XEMISF        = 0.0  
ELSE
  ALLOCATE(I%I%XFFLOOD   (0))
  ALLOCATE(I%I%XPIFLOOD  (0))
  ALLOCATE(I%I%XFF     (0,0))
  ALLOCATE(I%I%XFFG    (0,0))
  ALLOCATE(I%I%XFFV    (0,0))
  ALLOCATE(I%I%XFFROZEN(0,0))
  ALLOCATE(I%I%XALBF   (0,0))
  ALLOCATE(I%I%XEMISF  (0,0))
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      7.     ISBA time-varying deep force-restore temperature initialization
!              ---------------------------------------------------------------
!
 CALL SOILTEMP_ARP_PAR(I%O, HPROGRAM)
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
 CALL READ_ISBA_n(DTCO, I%O, I%R, I%P%XCLAY, I%M%T%XLAI, I%M%T%XBSLAI, U, HPROGRAM)
!
IF (HINIT/='ALL') THEN
  CALL END_IO_SURF_n(HPROGRAM)
  IF (LHOOK) CALL DR_HOOK('COMPUTE_ISBA_PARAMETERS',1,ZHOOK_HANDLE)
  RETURN
END IF
!
IF (HINIT=='PRE' .AND. I%R%TSNOW%SCHEME.NE.'3-L' .AND. &
        I%R%TSNOW%SCHEME.NE.'CRO' .AND. I%O%CISBA=='DIF') THEN
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
   CALL INIT_ISBA_LANDUSE(DTCO, UG, U, I%O, I%R, IG%XMESH_SIZE, &
                          I%M%X%XDG, I%M%X%XDG_OLD, I%M%T%XLAI, &
                          I%IP%XPATCH, I%IP%XPATCH_OLD, HPROGRAM)  
END IF
!
!-------------------------------------------------------------------------------
!
!*      12.     Canopy air fields:
!               -----------------
!
 CALL READ_ISBA_CANOPY_n(DTCO, ICP, I%O%LCANOPY, U, HPROGRAM)
!
!-------------------------------------------------------------------------------
!
!*      13.     initialize radiative and physical properties
!               --------------------------------------------
!
ALLOCATE(I%I%XDIR_ALB_WITH_SNOW(KI,KSW,I%O%NPATCH))
ALLOCATE(I%I%XSCA_ALB_WITH_SNOW(KI,KSW,I%O%NPATCH))
I%I%XDIR_ALB_WITH_SNOW = 0.0
I%I%XSCA_ALB_WITH_SNOW = 0.0
!
 CALL INIT_VEG_n(I%O, I%O%LCANOPY, KI, I%M%T, I%M%A, I%M%X%XH_TREE, I%R, I%IP%XVEGTYPE_PATCH, &
                 DGI%DM%LSURF_DIAG_ALBEDO,  &
                 PDIR_ALB, PSCA_ALB, PEMIS, PTSRAD )
!
DO JPATCH=1,I%O%NPATCH
  ZWG1(:,JPATCH) = I%R%XWG(:,1,JPATCH)
  ZTG1(:,JPATCH) = I%R%XTG(:,1,JPATCH)
END DO
!
 CALL CONVERT_PATCH_ISBA(DTCO, DTI, I%O, IDECADE,IDECADE2, I%P%XCOVER, I%P%LCOVER,&
                         LAGRIP, 'NAT', .FALSE., IMA=I%M%A, IP=I%IP, PWG1=ZWG1, PWSAT=I%IP%XWSAT)
!
! Load randomly perturbed fields. Perturbation ratios are saved in case fields are reset later.
IF(I%O%LPERTSURF) THEN
!
  CALL READ_SURF(HPROGRAM,'VEG',I%M%T%XVEG(:,:),IRESP)
  ALLOCATE(I%I%XPERTVEG(KI))
  I%I%XPERTVEG(:)=I%M%T%XVEG(:,1)
!
  CALL READ_SURF(HPROGRAM,'LAI',I%M%T%XLAI(:,:),IRESP)
  ALLOCATE(I%I%XPERTLAI(KI))
  I%I%XPERTLAI(:)=I%M%T%XLAI(:,1)
!
  CALL READ_SURF(HPROGRAM,'CV',I%M%T%XCV(:,:),IRESP)
  ALLOCATE(I%I%XPERTCV(KI))
  I%I%XPERTCV(:)=I%M%T%XCV(:,1)
!
  CALL READ_SURF(HPROGRAM,'PERTALB',ZPERTBUF(:,:),IRESP)
  ALLOCATE(I%I%XPERTALB(KI))
  I%I%XPERTALB(:)=ZPERTBUF(:,1)
  WHERE(I%M%T%XALBNIR_VEG(:,1)/=XUNDEF)  I%M%T%XALBNIR_VEG(:,1) = I%M%T%XALBNIR_VEG(:,1) *( 1.+ I%I%XPERTALB(:) )
  WHERE(I%M%T%XALBVIS_VEG(:,1)/=XUNDEF)  I%M%T%XALBVIS_VEG(:,1) = I%M%T%XALBVIS_VEG(:,1) *( 1.+ I%I%XPERTALB(:) )
  WHERE(I%M%T%XALBUV_VEG(:,1)/=XUNDEF)   I%M%T%XALBUV_VEG(:,1)  = I%M%T%XALBUV_VEG(:,1)  *( 1.+ I%I%XPERTALB(:) )
  WHERE(I%M%A%XALBNIR_SOIL(:,1)/=XUNDEF) I%M%A%XALBNIR_SOIL(:,1)= I%M%A%XALBNIR_SOIL(:,1)*( 1.+ I%I%XPERTALB(:) )
  WHERE(I%M%A%XALBVIS_SOIL(:,1)/=XUNDEF) I%M%A%XALBVIS_SOIL(:,1)= I%M%A%XALBVIS_SOIL(:,1)*( 1.+ I%I%XPERTALB(:) )
  WHERE(I%M%A%XALBUV_SOIL(:,1)/=XUNDEF)  I%M%A%XALBUV_SOIL(:,1) = I%M%A%XALBUV_SOIL(:,1) *( 1.+ I%I%XPERTALB(:) )
!
  CALL READ_SURF(HPROGRAM,'PERTZ0LAND',ZPERTBUF(:,:),IRESP)
  ALLOCATE(I%I%XPERTZ0(KI))
  I%I%XPERTZ0(:)=ZPERTBUF(:,1)
  WHERE(I%M%T%XZ0(:,1)/=XUNDEF)      I%M%T%XZ0(:,1)     =I%M%T%XZ0(:,1)     *( 1.+ I%I%XPERTZ0(:) )
  WHERE(ISS%XZ0EFFIP(:,1)/=XUNDEF) ISS%XZ0EFFIP(:,1)=ISS%XZ0EFFIP(:,1)*( 1.+ I%I%XPERTZ0(:) )
  WHERE(ISS%XZ0EFFIM(:,1)/=XUNDEF) ISS%XZ0EFFIM(:,1)=ISS%XZ0EFFIM(:,1)*( 1.+ I%I%XPERTZ0(:) )
  WHERE(ISS%XZ0EFFJP(:,1)/=XUNDEF) ISS%XZ0EFFJP(:,1)=ISS%XZ0EFFJP(:,1)*( 1.+ I%I%XPERTZ0(:) )
  WHERE(ISS%XZ0EFFJM(:,1)/=XUNDEF) ISS%XZ0EFFJM(:,1)=ISS%XZ0EFFJM(:,1)*( 1.+ I%I%XPERTZ0(:) )
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       14.    Output radiative fields
!               -----------------------
!
ALLOCATE(I%I%XEMIS_NAT   (KI))
I%I%XEMIS_NAT (:) = XUNDEF
!
 CALL AVERAGED_ALBEDO_EMIS_ISBA(I%O, I%M%T, I%M%M, I%M%A, I%IP, I%I, I%R, I%M%X%XVEGTYPE, &
                                PZENITH, ZTG1, PSW_BANDS, PDIR_ALB, PSCA_ALB, &
                                I%I%XEMIS_NAT,ZTSRAD_NAT,ZTSURF_NAT         )  
!
PEMIS  = I%I%XEMIS_NAT
PTSRAD = ZTSRAD_NAT
PTSURF = ZTSURF_NAT
!
!-------------------------------------------------------------------------------
!
!*      15.     ISBA diagnostics initialization
!               -------------------------------
!
IF(I%O%NPATCH<=1) DGI%O%LPATCH_BUDGET=.FALSE.
!
 CALL DIAG_ISBA_INIT_n(CHI, DGI%DE, DGI%DEC, DGI%DEP, DGI%DEPC, DGI%O, &
                       DGI%D, DGI%DC, DGI%DP, DGI%DPC, DGI%DM, DGI%DMP, &
                       OREAD_BUDGETC, GB, I%O, I%R%TSNOW%SCHEME, I%R%TSNOW%NLAYER, &
                       SIZE(I%IP%XABC), HPROGRAM,KI,KSW)
!
!-------------------------------------------------------------------------------
!
 CALL INIT_SURF_TOPD(DGI%DEC, I, UG, U, HPROGRAM,U%NDIM_FULL)
!
!-------------------------------------------------------------------------------
!
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('COMPUTE_ISBA_PARAMETERS',1,ZHOOK_HANDLE)
!
END SUBROUTINE COMPUTE_ISBA_PARAMETERS


