!#############################################################
SUBROUTINE INIT_TEB_GREENROOF_n(HPROGRAM,HINIT,KI,KSW,PSW_BANDS)
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
USE MODD_DATA_TEB_GREENROOF_n, ONLY : DTGR => DATA_TEB_GREENROOF
!
USE MODD_TYPE_DATE_SURF
USE MODD_TYPE_SNOW
!
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS
USE MODD_TEB_VEG_n, ONLY : TVG => TEB_VEG_OPTIONS
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TGRO => TEB_GREENROOF_OPTIONS
USE MODD_TEB_GREENROOF_PGD_n, ONLY : TGRP => TEB_GREENROOF_PGD

USE MODD_TEB_GREENROOF_PGD_EVOL_n, ONLY : TGRPE => TEB_GREENROOF_PGD_EVOL
USE MODD_TEB_GREENROOF_n, ONLY : TGR => TEB_GREENROOF
USE MODD_DIAG_MISC_TEB_OPTION_n, ONLY : DGMTO => DIAG_MISC_TEB_OPTIONS
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
USE MODI_GET_CURRENT_TEB_PATCH
USE MODI_READ_TEB_GREENROOF_n
USE MODI_INIT_VEG_GARDEN_n
USE MODI_SOIL_ALBEDO
USE MODI_INIT_FROM_DATA_GREENROOF_n
USE MODI_AVG_ALBEDO_EMIS_GREENROOF
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
IF (LHOOK) CALL DR_HOOK('INIT_TEB_GREENROOF_N',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
!-------------------------------------------------------------------------------
!
!*       1.     Reading of snow configuration:
!               ------------------------------
!
!* initialization of snow scheme (TSNOW defined in MODD_TEB_GREENROOF_n)
!
IF (HINIT=='PRE') THEN
   CALL READ_PREP_GREENROOF_SNOW(HPROGRAM,TGR%CUR%TSNOW%SCHEME,TGR%CUR%TSNOW%NLAYER)
!
   IF (TGR%CUR%TSNOW%SCHEME.NE.'3-L' .AND. TGR%CUR%TSNOW%SCHEME.NE.'CRO' .AND. TGRO%CISBA_GR=='DIF') THEN
    CALL ABOR1_SFX("INIT_TEB_GREENROOF_n: WITH CISBA_GR = DIF, CSNOW MUST BE 3-L OR CRO")
  ENDIF
  IF (LHOOK) CALL DR_HOOK('INIT_TEB_GREENROOF_N',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!-------------------------------------------------------------------------------
!
 CALL ALLOCATE_TEB_GREENROOF(TGR, TVG, &
                             KI, TGRO%NLAYER_GR)  
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
  CALL GET_CURRENT_TEB_PATCH(JTEB_PATCH)
  IF (TOP%NTEB_PATCH>1) WRITE(YPATCH,FMT='(A,I1,A)') 'T',JTEB_PATCH,'_'
!
  CALL READ_TEB_GREENROOF_n(HPROGRAM,YPATCH)
!
!
 CALL INIT_VEG_GARDEN_n(KI, TOP%LCANOPY, TVG%CROUGH, TGR%CUR%TSNOW, &
                   TVG%CPHOTO, TGRP%XLAIMIN, TGRP%XH_TREE, TGRP%XVEGTYPE, &
                   TGRPE%CUR%XLAI, TGRPE%CUR%XZ0, TGRPE%CUR%XVEG, TGRPE%CUR%XEMIS, &
                   TGRO%LTR_ML_GR, TGR%CUR%XFAPARC, TGR%CUR%XFAPIRC, TGR%CUR%XLAI_EFFC, TGR%CUR%XMUS, &
                   TGRP%XALBNIR_SOIL, TGRP%XALBVIS_SOIL, TGRP%XALBUV_SOIL, &
                   TGRPE%CUR%XALBNIR, TGRPE%CUR%XALBVIS, TGRPE%CUR%XALBUV, &
                   DGMTO%LSURF_DIAG_ALBEDO, TGR%CUR%XPSN, TGR%CUR%XPSNG, TGR%CUR%XPSNV, TGR%CUR%XPSNV_A, &
                   ZDIR_ALB, ZSCA_ALB, ZEMIS, ZTSRAD )
!
ZWG1(:) = TGR%CUR%XWG(:,1)
ZTG1(:) = TGR%CUR%XTG(:,1)
!
IF (.NOT. TGRO%LPAR_GREENROOF) THEN
  CALL SOIL_ALBEDO(TVG%CALBEDO,                               &
                     TGRP%XWSAT(:,1),ZWG1,                       &
                     TGRP%XALBVIS_DRY,TGRP%XALBNIR_DRY,TGRP%XALBUV_DRY,    &
                     TGRP%XALBVIS_WET,TGRP%XALBNIR_WET,TGRP%XALBUV_WET,    &
                     TGRP%XALBVIS_SOIL,TGRP%XALBNIR_SOIL,TGRP%XALBUV_SOIL  )  
ELSE
  IF (TOP%TTIME%TDATE%MONTH /= NUNDEF) THEN
    IDECADE = 3 * ( TOP%TTIME%TDATE%MONTH - 1 ) + MIN(TOP%TTIME%TDATE%DAY-1,29) / 10 + 1
  ELSE
    IDECADE = 1
  END IF
  CALL INIT_FROM_DATA_GREENROOF_n(DTGR, TGRO, &
                                  IDECADE,TVG%CPHOTO,              &
                                  PALBNIR_SOIL=TGRP%XALBNIR_SOIL,   &
                                  PALBVIS_SOIL=TGRP%XALBVIS_SOIL,   &
                                  PALBUV_SOIL=TGRP%XALBUV_SOIL      )  
END IF
!
! 
 CALL AVG_ALBEDO_EMIS_GREENROOF(TGR, &
                                TVG%CALBEDO,                                &
                               TGRPE%CUR%XVEG,TGRPE%CUR%XZ0,TGRPE%CUR%XLAI,ZTG1,                     &
                               PSW_BANDS,                              &
                               TGRP%XALBNIR_VEG,TGRP%XALBVIS_VEG,TGRP%XALBUV_VEG,     &
                               TGRP%XALBNIR_SOIL,TGRP%XALBVIS_SOIL,TGRP%XALBUV_SOIL,  &
                               TGRPE%CUR%XEMIS,                                  &
                               TGR%CUR%TSNOW,                                  &
                               TGRPE%CUR%XALBNIR,TGRPE%CUR%XALBVIS,TGRPE%CUR%XALBUV,                 &
                               ZDIR_ALB, ZSCA_ALB,                     &
                               ZEMIS,ZTSRAD                            )  
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('INIT_TEB_GREENROOF_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE INIT_TEB_GREENROOF_n
