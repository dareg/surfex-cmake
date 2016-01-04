!#############################################################
SUBROUTINE INIT_TEB_GREENROOF_n (DTCO, U, DGMTO, TOP, GDO, GRM, &
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
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_SURFEX_n, ONLY : TEB_GREENROOF_MODEL_t
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
TYPE(TEB_GREENROOF_MODEL_t), INTENT(INOUT) :: GRM
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
GRM%TV%O%CC1DRY = GDO%CC1DRY
GRM%TV%O%CSNOWRES = GDO%CSNOWRES
GRM%TV%O%CCPSURF = GDO%CCPSURF
GRM%TV%O%CSOILFRZ = GDO%CSOILFRZ
GRM%TV%O%CDIFSFCOND = GDO%CDIFSFCOND
GRM%TV%O%XCGMAX = GDO%XCGMAX
GRM%TV%O%CALBEDO = GDO%CALBEDO
GRM%TV%O%LNITRO_DILU = GDO%LNITRO_DILU
GRM%TV%O%CROUGH = GDO%CROUGH
GRM%TV%O%CRAIN = GDO%CRAIN
GRM%TV%O%XCDRAG = GDO%XCDRAG
GRM%TV%O%LCANOPY_DRAG = GDO%LCANOPY_DRAG
GRM%TV%O%LVEGUPD = GDO%LVEGUPD
!
!-------------------------------------------------------------------------------
!
!*       1.     Reading of snow configuration:
!               ------------------------------
!
!* initialization of snow scheme (TSNOW defined in MODD_TEB_GREENROOF_n)
!
IF (HINIT=='PRE') THEN
   CALL READ_PREP_GREENROOF_SNOW(HPROGRAM,GRM%TV%R%CUR%TSNOW%SCHEME,GRM%TV%R%CUR%TSNOW%NLAYER)
!
   IF (GRM%TV%R%CUR%TSNOW%SCHEME.NE.'3-L' .AND. GRM%TV%R%CUR%TSNOW%SCHEME.NE.'CRO' &
           .AND. GRM%TV%O%CISBA=='DIF') THEN
    CALL ABOR1_SFX("INIT_TEB_GREENROOF_n: WITH CISBA_GR = DIF, CSNOW MUST BE 3-L OR CRO")
  ENDIF
  IF (LHOOK) CALL DR_HOOK('INIT_TEB_GREENROOF_N',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!-------------------------------------------------------------------------------
!
 CALL ALLOCATE_TEB_GREENROOF(GRM%TV%R,  &
                             KI, GRM%TV%O%NGROUND_LAYER, GRM%TV%O%NNBIOMASS)  
!
!-------------------------------------------------------------------------------
!
IF( GRM%TV%O%CCPSURF=='DRY' .AND. LCPL_ARP ) THEN
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
  CALL READ_TEB_GREENROOF_n(DTCO, U, GRM, &
                            HPROGRAM,YPATCH)
!
 ZVEGTYPE_PATCH(:,:,1) = GRM%TV%M%X%XVEGTYPE(:,:)
 CALL INIT_VEG_n(1, KI, TOP%LCANOPY, GRM%TV%O%CROUGH, .FALSE., GRM%TV%R%CUR%TSNOW, &
                   GRM%TV%O%CPHOTO, GRM%TV%M%T%CUR%XLAIMIN, GRM%TV%M%X%XH_TREE, &
                   ZVEGTYPE_PATCH, &
                   GRM%TV%M%T%CUR%XLAI, GRM%TV%M%T%CUR%XZ0, GRM%TV%M%T%CUR%XVEG, &
                   GRM%TV%M%T%CUR%XEMIS, &
                   GRM%TV%O%LTR_ML, GRM%TV%R%CUR%XFAPARC, GRM%TV%R%CUR%XFAPIRC, &
                   GRM%TV%R%CUR%XLAI_EFFC, &
                   GRM%TV%R%CUR%XMUS, GRM%TV%M%A%XALBNIR_SOIL, GRM%TV%M%A%XALBVIS_SOIL, &
                   GRM%TV%M%A%XALBUV_SOIL,GRM%TV%M%T%CUR%XALBNIR, GRM%TV%M%T%CUR%XALBVIS, &
                   GRM%TV%M%T%CUR%XALBUV, DGMTO%LSURF_DIAG_ALBEDO, GRM%TV%R%CUR%XPSN, &
                   GRM%TV%R%CUR%XPSNG, GRM%TV%R%CUR%XPSNV, GRM%TV%R%CUR%XPSNV_A,&
                   ZDIR_ALB, ZSCA_ALB, ZEMIS, ZTSRAD )
!
ZWG1(:,1) = GRM%TV%R%CUR%XWG(:,1,1)
ZTG1(:,1) = GRM%TV%R%CUR%XTG(:,1,1)
!
IF (.NOT. GRM%TV%O%LPAR) THEN
  CALL SOIL_ALBEDO(GRM%TV%O%CALBEDO,                               &
                     GRM%TV%IP%XWSAT(:,1),ZWG1,                       &
                     GRM%TV%IP%XALBVIS_DRY,GRM%TV%IP%XALBNIR_DRY,GRM%TV%IP%XALBUV_DRY,    &
                     GRM%TV%IP%XALBVIS_WET,GRM%TV%IP%XALBNIR_WET,GRM%TV%IP%XALBUV_WET,    &
                     GRM%TV%M%A%XALBVIS_SOIL,GRM%TV%M%A%XALBNIR_SOIL,GRM%TV%M%A%XALBUV_SOIL  )  
ELSE
  IF (TOP%TTIME%TDATE%MONTH /= NUNDEF) THEN
    IDECADE = 3 * ( TOP%TTIME%TDATE%MONTH - 1 ) + MIN(TOP%TTIME%TDATE%DAY-1,29) / 10 + 1
  ELSE
    IDECADE = 1
  END IF
  CALL INIT_FROM_DATA_GREENROOF_n(GRM%DTI, &
                                  IDECADE,GRM%TV%O%CPHOTO,              &
                                  PALBNIR_SOIL=GRM%TV%M%A%XALBNIR_SOIL,   &
                                  PALBVIS_SOIL=GRM%TV%M%A%XALBVIS_SOIL,   &
                                  PALBUV_SOIL=GRM%TV%M%A%XALBUV_SOIL      )  
END IF
!
!
 CALL AVG_ALBEDO_EMIS_GARDEN(GRM%TV%R, GRM%TV%O%CALBEDO,                  &
                             GRM%TV%M%T%CUR%XVEG(:,1),GRM%TV%M%T%CUR%XZ0(:,1),GRM%TV%M%T%CUR%XLAI(:,1),ZTG1(:,1),      &
                             PSW_BANDS,                             &
                             GRM%TV%M%T%CUR%XALBNIR_VEG(:,1),GRM%TV%M%T%CUR%XALBVIS_VEG(:,1),&
                             GRM%TV%M%T%CUR%XALBUV_VEG(:,1),     &
                             GRM%TV%M%A%XALBNIR_SOIL(:,1),GRM%TV%M%A%XALBVIS_SOIL(:,1),GRM%TV%M%A%XALBUV_SOIL(:,1),  &
                             GRM%TV%M%T%CUR%XEMIS(:,1), GRM%TV%R%CUR%TSNOW,                            &
                             GRM%TV%M%T%CUR%XALBNIR(:,1),GRM%TV%M%T%CUR%XALBVIS(:,1),GRM%TV%M%T%CUR%XALBUV(:,1),  &
                                 ZDIR_ALB, ZSCA_ALB,                     &
                                 ZEMIS,ZTSRAD                            )  

!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('INIT_TEB_GREENROOF_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!

END SUBROUTINE INIT_TEB_GREENROOF_n
