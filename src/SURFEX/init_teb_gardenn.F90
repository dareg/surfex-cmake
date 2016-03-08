!#############################################################
SUBROUTINE INIT_TEB_GARDEN_n (DTCO, UG, U, DGMTO, TOP, &
                              GDO, DTGD, GDIP, GDI, GDR, GDMX, GDMT, GDMA, &
                              GDDG, &
                              HPROGRAM, HINIT, KI, KSW, PSW_BANDS, KPATCH)
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
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DIAG_n, ONLY : DIAG_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
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

USE MODD_DATA_COVER_PAR,  ONLY: NVEGTYPE
USE MODD_SURF_PAR,        ONLY: XUNDEF, NUNDEF

USE MODD_SURF_ATM,        ONLY: LCPL_ARP
!
USE MODI_GET_LUOUT
USE MODI_READ_PREP_GARDEN_SNOW
USE MODI_ALLOCATE_TEB_GARDEN
USE MODI_ABOR1_SFX
USE MODI_READ_TEB_GARDEN_n
USE MODI_INIT_VEG_n
USE MODI_SOIL_ALBEDO
USE MODI_INIT_FROM_DATA_GRDN_n
USE MODI_AVG_ALBEDO_EMIS_GARDEN
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
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(DIAG_MISC_TEB_OPTIONS_t), INTENT(INOUT) :: DGMTO
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: GDO
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTGD
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: GDIP
TYPE(ISBA_INIT_t), INTENT(INOUT) :: GDI
TYPE(ISBA_PROG_t), INTENT(INOUT) :: GDR
TYPE(ISBA_PARAM_FIX_t), INTENT(INOUT) :: GDMX
TYPE(ISBA_PARAM_TIME_t), INTENT(INOUT) :: GDMT
TYPE(ISBA_PARAM_ALB_t), INTENT(INOUT) :: GDMA
TYPE(TEB_VEG_DIAG_t), INTENT(INOUT) :: GDDG
!
 CHARACTER(LEN=6),                   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=3),                   INTENT(IN)  :: HINIT     ! choice of fields to initialize
INTEGER,                            INTENT(IN)  :: KI        ! number of points
INTEGER,                            INTENT(IN)  :: KSW       ! number of short-wave spectral bands
REAL,             DIMENSION(KSW),   INTENT(IN)  :: PSW_BANDS ! middle wavelength of each band
INTEGER,                            INTENT(IN)  :: KPATCH
!
!
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
IF (LHOOK) CALL DR_HOOK('INIT_TEB_GARDEN_N',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
GDO%CRAIN = "DEF"
!
!*       1.     Reading of snow configuration:
!               ------------------------------
!
!* initialization of snow scheme (TSNOW defined in MODD_TEB_GARDEN_n)
!
IF (HINIT=='PRE') THEN
  CALL READ_PREP_GARDEN_SNOW(HPROGRAM,GDR%TSNOW%SCHEME,GDR%TSNOW%NLAYER)
!
  IF (GDR%TSNOW%SCHEME.NE.'3-L' .AND. &
                GDR%TSNOW%SCHEME.NE.'CRO' .AND. GDO%CISBA=='DIF') THEN
    CALL ABOR1_SFX("INIT_TEB_GARDEN_n: WITH CISBA = DIF, CSNOW MUST BE 3-L OR CRO")
  ENDIF
  IF (LHOOK) CALL DR_HOOK('INIT_TEB_GARDEN_N',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!-------------------------------------------------------------------------------
!
 CALL ALLOCATE_TEB_GARDEN(GDR, KI, GDO%NGROUND_LAYER, GDO%NNBIOMASS)  
!
!-------------------------------------------------------------------------------
!
IF( GDO%CCPSURF=='DRY' .AND. LCPL_ARP ) THEN
  CALL ABOR1_SFX('CCPSURF=DRY must not be used with LCPL_ARP')
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
! Variables needed to run isba
!
ALLOCATE(GDI%XFFLOOD (KI))
ALLOCATE(GDI%XFF     (KI,1))
ALLOCATE(GDI%XFFG    (KI,1))
ALLOCATE(GDI%XFFV    (KI,1))
ALLOCATE(GDI%XFFROZEN(KI,1))
ALLOCATE(GDI%XALBF   (KI,1))
ALLOCATE(GDI%XEMISF  (KI,1))
GDI%XFFLOOD = 0.0
GDI%XFF = 0.0
GDI%XFFG = 0.0
GDI%XFFV = 0.0
GDI%XFFROZEN = 0.0
GDI%XALBF = 0.0
GDI%XEMISF = 0.0
!
ALLOCATE(GDI%XFSAT(KI))  
GDI%XFSAT(:) = 0.0
!
!-------------------------------------------------------------------------------
!
!*      10.     Prognostic and semi-prognostic fields
!               -------------------------------------
!
!* allocation of urban green area variables
!
!
  YPATCH='   '
  IF (TOP%NTEB_PATCH>1) WRITE(YPATCH,FMT='(A,I1,A)') 'T',KPATCH,'_'
!
  CALL READ_TEB_GARDEN_n(DTCO, U, GDO, GDR, GDMT, HPROGRAM,YPATCH)
!
!
ZVEGTYPE_PATCH(:,:,1) = GDMX%XVEGTYPE(:,:)
!
 CALL INIT_VEG_n(GDO, TOP%LCANOPY, KI, GDMT, GDMA, &
                 GDMX%XH_TREE, GDR, ZVEGTYPE_PATCH, &
                 DGMTO%LSURF_DIAG_ALBEDO, ZDIR_ALB, ZSCA_ALB, ZEMIS, ZTSRAD )
!
ZWG1(:,1) = GDR%XWG(:,1,1)
ZTG1(:,1) = GDR%XTG(:,1,1)
!
IF (.NOT. GDO%LPAR) THEN
  CALL SOIL_ALBEDO(GDO%CALBEDO, GDIP%XWSAT(:,1),ZWG1, GDIP, GDMA, "ALL" )  
ELSE
  IF (TOP%TTIME%TDATE%MONTH /= NUNDEF) THEN
    IDECADE = 3 * ( TOP%TTIME%TDATE%MONTH - 1 ) + MIN(TOP%TTIME%TDATE%DAY-1,29) / 10 + 1
  ELSE
    IDECADE = 1
  END IF
  CALL INIT_FROM_DATA_GRDN_n(DTGD, IDECADE,GDO%CPHOTO, .FALSE.,VMA=GDMA, VMIP=GDIP )  
END IF
!
 CALL AVG_ALBEDO_EMIS_GARDEN(GDR, GDO%CALBEDO,  GDMT,  GDMA, &
                             ZTG1(:,1), PSW_BANDS, ZDIR_ALB, ZSCA_ALB, ZEMIS,ZTSRAD      )  
!
 CALL DIAG_TEB_GARDEN_INIT_n(GDDG%D, GDDG%DP, GDDG%E, GDDG%EP, GDDG%M, KI, GDR%TSNOW%NLAYER)
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('INIT_TEB_GARDEN_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE INIT_TEB_GARDEN_n
