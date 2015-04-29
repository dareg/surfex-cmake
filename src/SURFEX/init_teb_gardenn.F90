!#############################################################
SUBROUTINE INIT_TEB_GARDEN_n(HPROGRAM,HINIT,KI,KSW,PSW_BANDS)
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
USE MODD_DATA_TEB_GARDEN_n, ONLY : DTGD => DATA_TEB_GARDEN
!
USE MODD_TYPE_DATE_SURF
USE MODD_TYPE_SNOW
!
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS
USE MODD_TEB_VEG_n, ONLY : TVG => TEB_VEG_OPTIONS
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TGDO => TEB_GARDEN_OPTIONS
USE MODD_TEB_GARDEN_PGD_EVOL_n, ONLY : TGDPE => TEB_GARDEN_PGD_EVOL
USE MODD_TEB_GARDEN_PGD_n, ONLY : TGDP => TEB_GARDEN_PGD
USE MODD_TEB_GARDEN_n, ONLY : TGD => TEB_GARDEN
USE MODD_DIAG_MISC_TEB_OPTION_n, ONLY : DGMTO => DIAG_MISC_TEB_OPTIONS

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
  CALL READ_PREP_GARDEN_SNOW(HPROGRAM,TGD%TSNOW%SCHEME,TGD%TSNOW%NLAYER)
!
  IF (TGD%TSNOW%SCHEME.NE.'3-L' .AND. TGD%TSNOW%SCHEME.NE.'CRO' .AND. TVG%CISBA=='DIF') THEN
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
  CALL READ_TEB_GARDEN_n(HPROGRAM,YPATCH)
!
!
 CALL INIT_VEG_GARDEN_n(KI, TOP%LCANOPY, TVG%CROUGH, TGD%TSNOW, &
                   TVG%CPHOTO, TGDP%XLAIMIN, TGDP%XH_TREE, TGDP%XVEGTYPE, TGDPE%XLAI, TGDPE%XZ0, TGDPE%XVEG, TGDPE%XEMIS, &
                   TVG%LTR_ML, TGD%XFAPARC, TGD%XFAPIRC, TGD%XLAI_EFFC, TGD%XMUS, &
                   TGDP%XALBNIR_SOIL, TGDP%XALBVIS_SOIL, TGDP%XALBUV_SOIL, TGDPE%XALBNIR, TGDPE%XALBVIS, TGDPE%XALBUV, &
                   DGMTO%LSURF_DIAG_ALBEDO, TGD%XPSN, TGD%XPSNG, TGD%XPSNV, TGD%XPSNV_A, &
                   ZDIR_ALB, ZSCA_ALB, ZEMIS, ZTSRAD )
!
ZWG1(:) = TGD%XWG(:,1)
ZTG1(:) = TGD%XTG(:,1)
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
                                 TGDPE%XVEG,TGDPE%XZ0,TGDPE%XLAI,ZTG1,                     &
                                 PSW_BANDS,                              &
                                 TGDP%XALBNIR_VEG,TGDP%XALBVIS_VEG,TGDP%XALBUV_VEG,     &
                                 TGDP%XALBNIR_SOIL,TGDP%XALBVIS_SOIL,TGDP%XALBUV_SOIL,  &
                                 TGDPE%XEMIS,                                  &
                                 TGD%TSNOW,                                  &
                                 TGDPE%XALBNIR,TGDPE%XALBVIS,TGDPE%XALBUV,                 &
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
