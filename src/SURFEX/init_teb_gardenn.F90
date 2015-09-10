!#############################################################
SUBROUTINE INIT_TEB_GARDEN_n (DTCO, DTGD, DGMTO, DGU, IOB, UG, &
                              U, TGD, TGDO, TGDPE, TGDP, TOP, TVG, &
                              HPROGRAM,HINIT,KI,KSW,PSW_BANDS)
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
USE MODD_DATA_TEB_GARDEN_n, ONLY : DATA_TEB_GARDEN_t
USE MODD_DIAG_MISC_TEB_OPTION_n, ONLY : DIAG_MISC_TEB_OPTIONS_t
USE MODD_DIAG_SURF_ATM_n, ONLY : DIAG_SURF_ATM_t
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_TEB_GARDEN_n, ONLY : TEB_GARDEN_t
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TEB_GARDEN_OPTIONS_t
USE MODD_TEB_GARDEN_PGD_EVOL_n, ONLY : TEB_GARDEN_PGD_EVOL_t
USE MODD_TEB_GARDEN_PGD_n, ONLY : TEB_GARDEN_PGD_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_OPTIONS_t
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
USE MODI_GET_CURRENT_TEB_PATCH
USE MODI_READ_TEB_GARDEN_n
USE MODI_INIT_VEG_GARDEN_n
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
TYPE(DATA_TEB_GARDEN_t), INTENT(INOUT) :: DTGD
TYPE(DIAG_MISC_TEB_OPTIONS_t), INTENT(INOUT) :: DGMTO
TYPE(DIAG_SURF_ATM_t), INTENT(INOUT) :: DGU
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(TEB_GARDEN_t), INTENT(INOUT) :: TGD
TYPE(TEB_GARDEN_OPTIONS_t), INTENT(INOUT) :: TGDO
TYPE(TEB_GARDEN_PGD_EVOL_t), INTENT(INOUT) :: TGDPE
TYPE(TEB_GARDEN_PGD_t), INTENT(INOUT) :: TGDP
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_VEG_OPTIONS_t), INTENT(INOUT) :: TVG
!
 CHARACTER(LEN=6),                   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=3),                   INTENT(IN)  :: HINIT     ! choice of fields to initialize
INTEGER,                            INTENT(IN)  :: KI        ! number of points
INTEGER,                            INTENT(IN)  :: KSW       ! number of short-wave spectral bands
REAL,             DIMENSION(KSW),   INTENT(IN)  :: PSW_BANDS ! middle wavelength of each band
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
INTEGER :: JTEB_PATCH  ! loop counter on TEB patches
 CHARACTER(LEN=3) :: YPATCH ! patch identificator
!
REAL, DIMENSION(KI)               :: ZWG1 ! work array for surface water content
REAL, DIMENSION(KI)               :: ZTG1 ! work array for surface temperature
REAL, DIMENSION(KI,KSW)           :: ZDIR_ALB  ! direct albedo for each band
REAL, DIMENSION(KI,KSW)           :: ZSCA_ALB  ! diffuse albedo for each band
REAL, DIMENSION(KI)               :: ZEMIS     ! emissivity
REAL, DIMENSION(KI)               :: ZTSRAD    ! radiative temperature
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
  CALL READ_PREP_GARDEN_SNOW(HPROGRAM,TGD%CUR%TSNOW%SCHEME,TGD%CUR%TSNOW%NLAYER)
!
  IF (TGD%CUR%TSNOW%SCHEME.NE.'3-L' .AND. TGD%CUR%TSNOW%SCHEME.NE.'CRO' .AND. TVG%CISBA=='DIF') THEN
    CALL ABOR1_SFX("INIT_TEB_GARDEN_n: WITH CISBA = DIF, CSNOW MUST BE 3-L OR CRO")
  ENDIF
  IF (LHOOK) CALL DR_HOOK('INIT_TEB_GARDEN_N',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!-------------------------------------------------------------------------------
!
 CALL ALLOCATE_TEB_GARDEN(TGD, TVG, &
                          KI, TGDO%NGROUND_LAYER)  
!
!-------------------------------------------------------------------------------
!
IF( TVG%CCPSURF=='DRY' .AND. LCPL_ARP ) THEN
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
  CALL GET_CURRENT_TEB_PATCH(JTEB_PATCH)
  IF (TOP%NTEB_PATCH>1) WRITE(YPATCH,FMT='(A,I1,A)') 'T',JTEB_PATCH,'_'
!
  CALL READ_TEB_GARDEN_n(DTCO, DGU, IOB, U, TGD, TGDO, TGDPE, TGDP, TVG, &
                         HPROGRAM,YPATCH)
!
!
 CALL INIT_VEG_GARDEN_n(KI, TOP%LCANOPY, TVG%CROUGH, TGD%CUR%TSNOW, &
                   TVG%CPHOTO, TGDP%XLAIMIN, TGDP%XH_TREE, TGDP%XVEGTYPE, TGDPE%CUR%XLAI, &
                   TGDPE%CUR%XZ0, TGDPE%CUR%XVEG, TGDPE%CUR%XEMIS, &
                   TVG%LTR_ML, TGD%CUR%XFAPARC, TGD%CUR%XFAPIRC, TGD%CUR%XLAI_EFFC, TGD%CUR%XMUS, &
                   TGDP%XALBNIR_SOIL, TGDP%XALBVIS_SOIL, TGDP%XALBUV_SOIL, TGDPE%CUR%XALBNIR, &
                   TGDPE%CUR%XALBVIS, TGDPE%CUR%XALBUV, &
                   DGMTO%LSURF_DIAG_ALBEDO, TGD%CUR%XPSN, TGD%CUR%XPSNG, TGD%CUR%XPSNV, TGD%CUR%XPSNV_A, &
                   ZDIR_ALB, ZSCA_ALB, ZEMIS, ZTSRAD )
!
ZWG1(:) = TGD%CUR%XWG(:,1)
ZTG1(:) = TGD%CUR%XTG(:,1)
!
IF (.NOT. TGDO%LPAR_GARDEN) THEN
  CALL SOIL_ALBEDO(TVG%CALBEDO,                               &
                     TGDP%XWSAT(:,1),ZWG1,                       &
                     TGDP%XALBVIS_DRY,TGDP%XALBNIR_DRY,TGDP%XALBUV_DRY,    &
                     TGDP%XALBVIS_WET,TGDP%XALBNIR_WET,TGDP%XALBUV_WET,    &
                     TGDP%XALBVIS_SOIL,TGDP%XALBNIR_SOIL,TGDP%XALBUV_SOIL  )  
ELSE
  IF (TOP%TTIME%TDATE%MONTH /= NUNDEF) THEN
    IDECADE = 3 * ( TOP%TTIME%TDATE%MONTH - 1 ) + MIN(TOP%TTIME%TDATE%DAY-1,29) / 10 + 1
  ELSE
    IDECADE = 1
  END IF
  CALL INIT_FROM_DATA_GRDN_n(DTGD, &
                             IDECADE,TVG%CPHOTO,              &
                               PALBNIR_SOIL=TGDP%XALBNIR_SOIL,   &
                               PALBVIS_SOIL=TGDP%XALBVIS_SOIL,   &
                               PALBUV_SOIL=TGDP%XALBUV_SOIL      )  
END IF
!
 CALL AVG_ALBEDO_EMIS_GARDEN(TGD, &
                             TVG%CALBEDO,                                   &
                                 TGDPE%CUR%XVEG,TGDPE%CUR%XZ0,TGDPE%CUR%XLAI,ZTG1,                     &
                                 PSW_BANDS,                              &
                                 TGDP%XALBNIR_VEG,TGDP%XALBVIS_VEG,TGDP%XALBUV_VEG,     &
                                 TGDP%XALBNIR_SOIL,TGDP%XALBVIS_SOIL,TGDP%XALBUV_SOIL,  &
                                 TGDPE%CUR%XEMIS,                                  &
                                 TGD%CUR%TSNOW,                                  &
                                 TGDPE%CUR%XALBNIR,TGDPE%CUR%XALBVIS,TGDPE%CUR%XALBUV,                 &
                                 ZDIR_ALB, ZSCA_ALB,                     &
                                 ZEMIS,ZTSRAD                            )  
!
!
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('INIT_TEB_GARDEN_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE INIT_TEB_GARDEN_n
