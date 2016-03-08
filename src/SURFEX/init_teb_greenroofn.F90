!#############################################################
SUBROUTINE INIT_TEB_GREENROOF_n (DTCO, U, DGMTO, TOP, GDO, &
                                 GRO, DTGR, GRIP, GRI, GRR, GRMX, GRMT, GRMA, &
                                 GRDG, &
                                 HPROGRAM,HINIT,KI,KSW,PSW_BANDS,KPATCH)
!#############################################################
!
!!****  *INIT_TEB_GREENROOF_n* - routine to initialize ISBA
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
!!      A. Lemonsu  *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    09/2009
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_DIAG_MISC_TEB_OPTION_n, ONLY : DIAG_MISC_TEB_OPTIONS_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t, ISBA_INIT_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_FIX_t, ISBA_PARAM_TIME_t, ISBA_PARAM_ALB_t
USE MODD_SURFEX_n, ONLY : TEB_VEG_DIAG_t
!
USE MODD_TYPE_DATE_SURF
USE MODD_TYPE_SNOW
!

!
USE MODD_DATA_COVER_PAR,       ONLY: NVEGTYPE
USE MODD_SURF_PAR,             ONLY: XUNDEF, NUNDEF
!
USE MODD_SURF_ATM,             ONLY: LCPL_ARP
!
USE MODI_GET_LUOUT
USE MODI_READ_PREP_GREENROOF_SNOW
USE MODI_ALLOCATE_TEB_GREENROOF
USE MODI_ABOR1_SFX
USE MODI_READ_TEB_GREENROOF_n
USE MODI_INIT_VEG_n
USE MODI_SOIL_ALBEDO
USE MODI_INIT_FROM_DATA_GREENROOF_n
USE MODI_AVG_ALBEDO_EMIS_GARDEN
USE MODI_DIAG_TEB_GREENROOF_INIT_n
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
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(DIAG_MISC_TEB_OPTIONS_t), INTENT(INOUT) :: DGMTO
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: GDO
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: GRO
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTGR
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: GRIP
TYPE(ISBA_INIT_t), INTENT(INOUT) :: GRI
TYPE(ISBA_PROG_t), INTENT(INOUT) :: GRR
TYPE(ISBA_PARAM_FIX_t), INTENT(INOUT) :: GRMX
TYPE(ISBA_PARAM_TIME_t), INTENT(INOUT) :: GRMT
TYPE(ISBA_PARAM_ALB_t), INTENT(INOUT) :: GRMA
TYPE(TEB_VEG_DIAG_t), INTENT(INOUT) :: GRDG
!
 CHARACTER(LEN=6),                   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=3),                   INTENT(IN)  :: HINIT     ! choice of fields to initialize
INTEGER,                            INTENT(IN)  :: KI        ! number of points
INTEGER,                            INTENT(IN)  :: KSW       ! number of short-wave spectral bands
REAL,             DIMENSION(KSW),   INTENT(IN)  :: PSW_BANDS ! middle wavelength of each band
INTEGER,                            INTENT(IN)  :: KPATCH
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: ILUOUT   ! unit of output listing file
!
INTEGER           :: IDECADE  ! decade of simulation
!
 CHARACTER(LEN=3) :: YPATCH ! patch identificator
!
REAL, DIMENSION(KI,1)               :: ZWG1 ! work array for surface water content
REAL, DIMENSION(KI,1)               :: ZTG1 ! work array for surface temperature
REAL, DIMENSION(KI,KSW)           :: ZDIR_ALB  ! direct albedo for each band
REAL, DIMENSION(KI,KSW)           :: ZSCA_ALB  ! diffuse albedo for each band
REAL, DIMENSION(KI)               :: ZEMIS     ! emissivity
REAL, DIMENSION(KI)               :: ZTSRAD    ! radiative temperature
!
REAL, DIMENSION(KI,NVEGTYPE,1) :: ZVEGTYPE_PATCH
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!               Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('INIT_TEB_GREENROOF_N',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
GRO%CC1DRY = GDO%CC1DRY
GRO%CSNOWRES = GDO%CSNOWRES
GRO%CCPSURF = GDO%CCPSURF
GRO%CSOILFRZ = GDO%CSOILFRZ
GRO%CDIFSFCOND = GDO%CDIFSFCOND
GRO%XCGMAX = GDO%XCGMAX
GRO%CALBEDO = GDO%CALBEDO
GRO%LNITRO_DILU = GDO%LNITRO_DILU
GRO%CROUGH = GDO%CROUGH
GRO%CRAIN = GDO%CRAIN
GRO%XCDRAG = GDO%XCDRAG
GRO%LCANOPY_DRAG = GDO%LCANOPY_DRAG
GRO%LVEGUPD = GDO%LVEGUPD
GRO%CRAIN = GDO%CRAIN
!
!-------------------------------------------------------------------------------
!
!*       1.     Reading of snow configuration:
!               ------------------------------
!
!* initialization of snow scheme (TSNOW defined in MODD_TEB_GREENROOF_n)
!
IF (HINIT=='PRE') THEN
   CALL READ_PREP_GREENROOF_SNOW(HPROGRAM,GRR%TSNOW%SCHEME,GRR%TSNOW%NLAYER)
!
   IF (GRR%TSNOW%SCHEME.NE.'3-L' .AND. GRR%TSNOW%SCHEME.NE.'CRO' &
           .AND. GRO%CISBA=='DIF') THEN
    CALL ABOR1_SFX("INIT_TEB_GREENROOF_n: WITH CISBA_GR = DIF, CSNOW MUST BE 3-L OR CRO")
  ENDIF
  IF (LHOOK) CALL DR_HOOK('INIT_TEB_GREENROOF_N',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!-------------------------------------------------------------------------------
!
 CALL ALLOCATE_TEB_GREENROOF(GRR, KI, GRO%NGROUND_LAYER, GRO%NNBIOMASS)  
!
!-------------------------------------------------------------------------------
!
IF( GRO%CCPSURF=='DRY' .AND. LCPL_ARP ) THEN
  CALL ABOR1_SFX('CCPSURF=DRY must not be used with LCPL_ARP')
ENDIF
!
!-------------------------------------------------------------------------------
!
IF (HINIT/='ALL') THEN
  IF (LHOOK) CALL DR_HOOK('INIT_TEB_GREENROOF_N',1,ZHOOK_HANDLE)      
  RETURN
ENDIF
!
!-------------------------------------------------------------------------------
! Variables needed to run isba
!
ALLOCATE(GRI%XFFLOOD (KI))
ALLOCATE(GRI%XFF     (KI,1))
ALLOCATE(GRI%XFFG    (KI,1))
ALLOCATE(GRI%XFFV    (KI,1))
ALLOCATE(GRI%XFFROZEN(KI,1))
ALLOCATE(GRI%XALBF   (KI,1))
ALLOCATE(GRI%XEMISF  (KI,1))
GRI%XFFLOOD = 0.0
GRI%XFF = 0.0
GRI%XFFG = 0.0
GRI%XFFV = 0.0
GRI%XFFROZEN = 0.0
GRI%XALBF = 0.0
GRI%XEMISF = 0.0
!
ALLOCATE(GRI%XFSAT(KI))  
GRI%XFSAT(:) = 0.0
!
!-------------------------------------------------------------------------------
!
!*       2.     Prognostic and semi-prognostic fields
!               -------------------------------------
!
!* allocation of urban green area variables
!
!
  YPATCH='   '
  IF (TOP%NTEB_PATCH>1) WRITE(YPATCH,FMT='(A,I1,A)') 'T',KPATCH,'_'
!
  CALL READ_TEB_GREENROOF_n(DTCO, U, GRO, GRR, GRMT, HPROGRAM,YPATCH)
!
 ZVEGTYPE_PATCH(:,:,1) = GRMX%XVEGTYPE(:,:)
 !
  CALL INIT_VEG_n(GRO, TOP%LCANOPY, KI, GRMT, GRMA, &
                 GRMX%XH_TREE, GRR, ZVEGTYPE_PATCH, &
                 DGMTO%LSURF_DIAG_ALBEDO, ZDIR_ALB, ZSCA_ALB, ZEMIS, ZTSRAD )
!
ZWG1(:,1) = GRR%XWG(:,1,1)
ZTG1(:,1) = GRR%XTG(:,1,1)
!
IF (.NOT. GRO%LPAR) THEN
  CALL SOIL_ALBEDO(GRO%CALBEDO, GRIP%XWSAT(:,1),ZWG1, GRIP, GRMA, "ALL" )
ELSE
  IF (TOP%TTIME%TDATE%MONTH /= NUNDEF) THEN
    IDECADE = 3 * ( TOP%TTIME%TDATE%MONTH - 1 ) + MIN(TOP%TTIME%TDATE%DAY-1,29) / 10 + 1
  ELSE
    IDECADE = 1
  END IF
  CALL INIT_FROM_DATA_GREENROOF_n(DTGR, IDECADE,GRO%CPHOTO, .FALSE., VMA=GRMA, VMIP = GRIP)
END IF
!
!
!
 CALL AVG_ALBEDO_EMIS_GARDEN(GRR, GRO%CALBEDO, GRMT,GRMA, &
                             ZTG1(:,1), PSW_BANDS,  ZDIR_ALB, ZSCA_ALB, ZEMIS,ZTSRAD )  
!
 CALL DIAG_TEB_GREENROOF_INIT_n(GRDG%D, GRDG%DP, GRDG%E, GRDG%EP, GRDG%M, &
                                KI, GRR%TSNOW%NLAYER)
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('INIT_TEB_GREENROOF_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!

END SUBROUTINE INIT_TEB_GREENROOF_n
