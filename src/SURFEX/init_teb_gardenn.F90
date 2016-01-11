!#############################################################
SUBROUTINE INIT_TEB_GARDEN_n (DTCO, DGU, UG, U, DGMTO, TOP, GDM, &
                              HPROGRAM,HINIT,KI,KSW,PSW_BANDS,KPATCH)
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
USE MODD_SURFEX_n, ONLY : TEB_GARDEN_MODEL_t
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
TYPE(DIAG_t), INTENT(INOUT) :: DGU
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(DIAG_MISC_TEB_OPTIONS_t), INTENT(INOUT) :: DGMTO
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_GARDEN_MODEL_t), INTENT(INOUT) :: GDM
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
!*       1.     Reading of snow configuration:
!               ------------------------------
!
!* initialization of snow scheme (TSNOW defined in MODD_TEB_GARDEN_n)
!
IF (HINIT=='PRE') THEN
  CALL READ_PREP_GARDEN_SNOW(HPROGRAM,GDM%TV%R%CUR%TSNOW%SCHEME,GDM%TV%R%CUR%TSNOW%NLAYER)
!
  IF (GDM%TV%R%CUR%TSNOW%SCHEME.NE.'3-L' .AND. &
                GDM%TV%R%CUR%TSNOW%SCHEME.NE.'CRO' .AND. GDM%TV%O%CISBA=='DIF') THEN
    CALL ABOR1_SFX("INIT_TEB_GARDEN_n: WITH CISBA = DIF, CSNOW MUST BE 3-L OR CRO")
  ENDIF
  IF (LHOOK) CALL DR_HOOK('INIT_TEB_GARDEN_N',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!-------------------------------------------------------------------------------
!
 CALL ALLOCATE_TEB_GARDEN(GDM%TV%R,  &
                          KI, GDM%TV%O%NGROUND_LAYER, GDM%TV%O%NNBIOMASS)  
!
!-------------------------------------------------------------------------------
!
IF( GDM%TV%O%CCPSURF=='DRY' .AND. LCPL_ARP ) THEN
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
  CALL READ_TEB_GARDEN_n(DTCO, DGU, U, GDM, &
                         HPROGRAM,YPATCH)
!
!
ZVEGTYPE_PATCH(:,:,1) = GDM%TV%M%X%XVEGTYPE(:,:)
!
 CALL INIT_VEG_n(1, KI, TOP%LCANOPY, GDM%TV%O%CROUGH, .FALSE., GDM%TV%R%CUR%TSNOW, &
                   GDM%TV%O%CPHOTO, GDM%TV%M%T%CUR%XLAIMIN, GDM%TV%M%X%XH_TREE, &
                   ZVEGTYPE_PATCH, &
                   GDM%TV%M%T%CUR%XLAI, GDM%TV%M%T%CUR%XZ0, GDM%TV%M%T%CUR%XVEG, &
                   GDM%TV%M%T%CUR%XEMIS, &
                   GDM%TV%O%LTR_ML, GDM%TV%R%CUR%XFAPARC, GDM%TV%R%CUR%XFAPIRC, &
                   GDM%TV%R%CUR%XLAI_EFFC, &
                   GDM%TV%R%CUR%XMUS, GDM%TV%M%A%XALBNIR_SOIL, GDM%TV%M%A%XALBVIS_SOIL, &
                   GDM%TV%M%A%XALBUV_SOIL,GDM%TV%M%T%CUR%XALBNIR, GDM%TV%M%T%CUR%XALBVIS, &
                   GDM%TV%M%T%CUR%XALBUV, DGMTO%LSURF_DIAG_ALBEDO, GDM%TV%R%CUR%XPSN, &
                   GDM%TV%R%CUR%XPSNG, GDM%TV%R%CUR%XPSNV, GDM%TV%R%CUR%XPSNV_A,&
                   ZDIR_ALB, ZSCA_ALB, ZEMIS, ZTSRAD )
!
ZWG1(:,1) = GDM%TV%R%CUR%XWG(:,1,1)
ZTG1(:,1) = GDM%TV%R%CUR%XTG(:,1,1)
!
IF (.NOT. GDM%TV%O%LPAR) THEN
  CALL SOIL_ALBEDO(GDM%TV%O%CALBEDO,                               &
                     GDM%TV%IP%XWSAT(:,1),ZWG1,                       &
                     GDM%TV%IP%XALBVIS_DRY,GDM%TV%IP%XALBNIR_DRY,GDM%TV%IP%XALBUV_DRY,    &
                     GDM%TV%IP%XALBVIS_WET,GDM%TV%IP%XALBNIR_WET,GDM%TV%IP%XALBUV_WET,    &
                     GDM%TV%M%A%XALBVIS_SOIL,GDM%TV%M%A%XALBNIR_SOIL,GDM%TV%M%A%XALBUV_SOIL  )  
ELSE
  IF (TOP%TTIME%TDATE%MONTH /= NUNDEF) THEN
    IDECADE = 3 * ( TOP%TTIME%TDATE%MONTH - 1 ) + MIN(TOP%TTIME%TDATE%DAY-1,29) / 10 + 1
  ELSE
    IDECADE = 1
  END IF
  CALL INIT_FROM_DATA_GRDN_n(GDM%DTI, &
                             IDECADE,GDM%TV%O%CPHOTO,              &
                               PALBNIR_SOIL=GDM%TV%M%A%XALBNIR_SOIL,   &
                               PALBVIS_SOIL=GDM%TV%M%A%XALBVIS_SOIL,   &
                               PALBUV_SOIL=GDM%TV%M%A%XALBUV_SOIL      )  
END IF
!
 CALL AVG_ALBEDO_EMIS_GARDEN(GDM%TV%R, GDM%TV%O%CALBEDO,                  &
                             GDM%TV%M%T%CUR%XVEG(:,1),GDM%TV%M%T%CUR%XZ0(:,1),GDM%TV%M%T%CUR%XLAI(:,1),ZTG1(:,1),      &
                             PSW_BANDS,                             &
                             GDM%TV%M%T%CUR%XALBNIR_VEG(:,1),GDM%TV%M%T%CUR%XALBVIS_VEG(:,1),&
                             GDM%TV%M%T%CUR%XALBUV_VEG(:,1),     &
                             GDM%TV%M%A%XALBNIR_SOIL(:,1),GDM%TV%M%A%XALBVIS_SOIL(:,1),GDM%TV%M%A%XALBUV_SOIL(:,1),  &
                             GDM%TV%M%T%CUR%XEMIS(:,1), GDM%TV%R%CUR%TSNOW,                            &
                             GDM%TV%M%T%CUR%XALBNIR(:,1),GDM%TV%M%T%CUR%XALBVIS(:,1),GDM%TV%M%T%CUR%XALBUV(:,1),  &
                                 ZDIR_ALB, ZSCA_ALB,                     &
                                 ZEMIS,ZTSRAD                            )  
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('INIT_TEB_GARDEN_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE INIT_TEB_GARDEN_n
