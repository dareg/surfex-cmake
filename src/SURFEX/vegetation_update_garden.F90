!     #########
    SUBROUTINE VEGETATION_UPDATE_GARDEN (DTCO, DTI, IG, I, T, TOP, DTGR, GDM, &
                                         TPTIME,PTSTEP,KLU)  
!   ##########################################################################
!
!!****  *GARDEN*  
!!
!!    PURPOSE
!!    -------
!
!     
!!**  METHOD
!     ------
!
!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!      A. Lemonsu          * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!    Original    05/2009
!!    P. Samuelsson  10/2014  Introduced MEB dummy variables in call to VEGETATION_UPDATE
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_GRID_n, ONLY : GRID_t
USE MODD_ISBA_n, ONLY : ISBA_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_SURFEX_n, ONLY : TEB_GARDEN_MODEL_t
!
USE MODD_TYPE_DATE_SURF,    ONLY: DATE_TIME


  
!
USE MODI_VEGETATION_UPDATE
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    Declarations of arguments
!
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTI
TYPE(GRID_t), INTENT(INOUT) :: IG
TYPE(ISBA_t), INTENT(INOUT) :: I
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTGR
TYPE(TEB_GARDEN_MODEL_t), INTENT(INOUT) :: GDM
!
TYPE(DATE_TIME)     , INTENT(IN)    :: TPTIME             ! current date and time from teb
REAL                , INTENT(IN)    :: PTSTEP             ! time step
INTEGER,              INTENT(IN)    :: KLU                ! number of points
!
!
!*      0.2    Declarations of local variables
!
REAL, DIMENSION(KLU,1) :: ZZ0EFFIP
REAL, DIMENSION(KLU,1) :: ZZ0EFFIM
REAL, DIMENSION(KLU,1) :: ZZ0EFFJP
REAL, DIMENSION(KLU,1) :: ZZ0EFFJM
REAL, DIMENSION(KLU) :: ZAOSIP
REAL, DIMENSION(KLU) :: ZAOSIM
REAL, DIMENSION(KLU) :: ZAOSJP
REAL, DIMENSION(KLU) :: ZAOSJM
REAL, DIMENSION(KLU) :: ZHO2IP
REAL, DIMENSION(KLU) :: ZHO2IM
REAL, DIMENSION(KLU) :: ZHO2JP
REAL, DIMENSION(KLU) :: ZHO2JM
REAL, DIMENSION(KLU,1) :: ZLAI
REAL, DIMENSION(KLU,1) :: ZVEG
REAL, DIMENSION(KLU,1) :: ZZ0
REAL, DIMENSION(KLU,1) :: ZALBNIR
REAL, DIMENSION(KLU,1) :: ZALBVIS
REAL, DIMENSION(KLU,1) :: ZALBUV
REAL, DIMENSION(KLU,1) :: ZEMIS
REAL, DIMENSION(KLU,1) :: ZRSMIN
REAL, DIMENSION(KLU,1) :: ZGAMMA
REAL, DIMENSION(KLU,1) :: ZWRMAX_CF
REAL, DIMENSION(KLU,1) :: ZRGL
REAL, DIMENSION(KLU,1) :: ZCV
REAL, DIMENSION(KLU,1) :: ZGMES
REAL, DIMENSION(KLU,1) :: ZBSLAI
REAL, DIMENSION(KLU,1) :: ZLAIMIN
REAL, DIMENSION(KLU,1) :: ZSEFOLD
REAL, DIMENSION(KLU,1) :: ZGC
REAL, DIMENSION(KLU,1) :: ZF2I
LOGICAL, DIMENSION(KLU,1) :: GSTRESS
REAL, DIMENSION(KLU,1) :: ZALBNIR_VEG
REAL, DIMENSION(KLU,1) :: ZALBVIS_VEG
REAL, DIMENSION(KLU,1) :: ZALBUV_VEG
REAL, DIMENSION(KLU,1) :: ZALBNIR_SOIL
REAL, DIMENSION(KLU,1) :: ZALBVIS_SOIL
REAL, DIMENSION(KLU,1) :: ZALBUV_SOIL
REAL, DIMENSION(KLU,1) :: ZCE_NITRO
REAL, DIMENSION(KLU,1) :: ZCF_NITRO
REAL, DIMENSION(KLU,1) :: ZCNA_NITRO
!
! MEB stuff
REAL, DIMENSION(KLU,1) :: ZGNDLITTER
REAL, DIMENSION(KLU,1) :: ZRGLGV
REAL, DIMENSION(KLU,1) :: ZGAMMAGV
REAL, DIMENSION(KLU,1) :: ZRSMINGV
REAL, DIMENSION(KLU,1) :: ZWRMAX_CFGV
REAL, DIMENSION(KLU,1) :: ZH_VEG
REAL, DIMENSION(KLU,1) :: ZLAIGV
REAL, DIMENSION(KLU,1) :: ZZ0LITTER
!
TYPE (DATE_TIME),  DIMENSION(KLU,1) :: TZSEED
TYPE (DATE_TIME), DIMENSION(KLU,1) :: TZREAP
REAL, DIMENSION(KLU,1) :: ZWATSUP
REAL, DIMENSION(KLU,1) :: ZIRRIG
LOGICAL :: GUPDATED              ! T if VEGETATION_UPDATE has reset fields
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*      1.     various initialisations
!              -----------------------
!
IF (LHOOK) CALL DR_HOOK('VEGETATION_UPDATE_GARDEN',0,ZHOOK_HANDLE)
!
!* orographic roughness not used
!
ZAOSIP = 0.
ZAOSIM = 0.
ZAOSJP = 0.
ZAOSJM = 0.
ZHO2IP = 0.
ZHO2IM = 0.
ZHO2JP = 0.
ZHO2JM = 0.
!
!* vegetation parameters to update
!
ZVEG = GDM%TV%M%T%CUR%XVEG
ZZ0 = GDM%TV%M%T%CUR%XZ0
ZALBNIR = GDM%TV%M%T%CUR%XALBNIR
ZALBVIS = GDM%TV%M%T%CUR%XALBVIS
ZALBUV = GDM%TV%M%T%CUR%XALBUV
ZEMIS = GDM%TV%M%T%CUR%XEMIS
ZRSMIN = GDM%TV%M%T%CUR%XRSMIN
ZGAMMA = GDM%TV%M%T%CUR%XGAMMA
ZWRMAX_CF = GDM%TV%M%T%CUR%XWRMAX_CF
ZRGL = GDM%TV%M%T%CUR%XRGL
ZCV = GDM%TV%M%T%CUR%XCV
ZGMES = GDM%TV%M%T%CUR%XGMES
ZBSLAI = GDM%TV%M%T%CUR%XBSLAI
ZLAIMIN = GDM%TV%M%T%CUR%XLAIMIN
ZSEFOLD = GDM%TV%M%T%CUR%XSEFOLD
ZGC = GDM%TV%M%T%CUR%XGC
ZF2I = GDM%TV%M%T%CUR%XF2I
GSTRESS = GDM%TV%M%T%CUR%LSTRESS
ZALBNIR_VEG = GDM%TV%M%T%CUR%XALBNIR_VEG
ZALBVIS_VEG = GDM%TV%M%T%CUR%XALBVIS_VEG
ZALBUV_VEG = GDM%TV%M%T%CUR%XALBUV_VEG
ZALBNIR_SOIL = GDM%TV%M%A%XALBNIR_SOIL
ZALBVIS_SOIL = GDM%TV%M%A%XALBVIS_SOIL
ZALBUV_SOIL = GDM%TV%M%A%XALBUV_SOIL
ZCE_NITRO = GDM%TV%M%T%CUR%XCE_NITRO
ZCF_NITRO = GDM%TV%M%T%CUR%XCF_NITRO
ZCNA_NITRO = GDM%TV%M%T%CUR%XCNA_NITRO
! --------------------------------------------------------------------------------------
! Vegetation update (in case of non-interactive vegetation):
! --------------------------------------------------------------------------------------
!
GUPDATED=.FALSE.
IF (GDM%TV%O%CPHOTO=='NON' .OR. GDM%TV%O%CPHOTO=='AGS' .OR. GDM%TV%O%CPHOTO=='AST') THEN
     CALL VEGETATION_UPDATE(DTCO, GDM%DTI, IG, I%O,  &
                            PTSTEP,TPTIME,TOP%XCOVER,TOP%LCOVER,                 &
                         GDM%TV%O%CISBA,(.NOT. GDM%TV%O%LPAR), &
                         GDM%TV%O%CPHOTO, .FALSE.,     &
                         GDM%TV%O%LTR_ML, 'GRD',                                  &
                         ZLAI,ZVEG,ZZ0,                                  &
                         ZALBNIR,ZALBVIS,ZALBUV,ZEMIS,                   &
                         ZRSMIN,ZGAMMA,ZWRMAX_CF,                        &
                         ZRGL,ZCV,                                       &
                         ZGMES,ZBSLAI,ZLAIMIN,ZSEFOLD,ZGC,         &
                         ZF2I, GSTRESS,                                  &
                         ZAOSIP,ZAOSIM,ZAOSJP,ZAOSJM,                    &
                         ZHO2IP,ZHO2IM,ZHO2JP,ZHO2JM,                    &
                         ZZ0EFFIP,ZZ0EFFIM,ZZ0EFFJP,ZZ0EFFJM,            &
                         GDM%TV%O%CALBEDO, ZALBNIR_VEG, ZALBVIS_VEG, ZALBUV_VEG,  &
                         ZALBNIR_SOIL, ZALBVIS_SOIL, ZALBUV_SOIL,        &
                         ZCE_NITRO, ZCF_NITRO, ZCNA_NITRO,               &
                         TZSEED, TZREAP, ZWATSUP, ZIRRIG,                &
                         ZGNDLITTER,ZRGLGV,ZGAMMAGV,                     &
                         ZRSMINGV, ZWRMAX_CFGV,                          &
                         ZH_VEG, ZLAIGV, ZZ0LITTER,                      &
                         GUPDATED, OABSENT=(T%CUR%XGARDEN==0.)                 )
END IF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
GDM%TV%M%T%CUR%XVEG = ZVEG
GDM%TV%M%T%CUR%XZ0 = ZZ0
GDM%TV%M%T%CUR%XALBNIR = ZALBNIR
GDM%TV%M%T%CUR%XALBVIS = ZALBVIS
GDM%TV%M%T%CUR%XALBUV = ZALBUV
GDM%TV%M%T%CUR%XEMIS = ZEMIS
GDM%TV%M%T%CUR%XRSMIN = ZRSMIN
GDM%TV%M%T%CUR%XGAMMA = ZGAMMA
GDM%TV%M%T%CUR%XWRMAX_CF = ZWRMAX_CF
GDM%TV%M%T%CUR%XRGL = ZRGL
GDM%TV%M%T%CUR%XCV = ZCV
GDM%TV%M%T%CUR%XGMES = ZGMES
GDM%TV%M%T%CUR%XBSLAI = ZBSLAI
GDM%TV%M%T%CUR%XLAIMIN = ZLAIMIN
GDM%TV%M%T%CUR%XSEFOLD = ZSEFOLD
GDM%TV%M%T%CUR%XGC = ZGC
GDM%TV%M%T%CUR%XF2I = ZF2I
GDM%TV%M%T%CUR%LSTRESS = GSTRESS
GDM%TV%M%T%CUR%XALBNIR_VEG = ZALBNIR_VEG
GDM%TV%M%T%CUR%XALBVIS_VEG = ZALBVIS_VEG
GDM%TV%M%T%CUR%XALBUV_VEG = ZALBUV_VEG
GDM%TV%M%A%XALBNIR_SOIL = ZALBNIR_SOIL
GDM%TV%M%A%XALBVIS_SOIL = ZALBVIS_SOIL
GDM%TV%M%A%XALBUV_SOIL = ZALBUV_SOIL
GDM%TV%M%T%CUR%XCE_NITRO = ZCE_NITRO
GDM%TV%M%T%CUR%XCF_NITRO = ZCF_NITRO
GDM%TV%M%T%CUR%XCNA_NITRO = ZCNA_NITRO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
IF (LHOOK) CALL DR_HOOK('VEGETATION_UPDATE_GARDEN',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE VEGETATION_UPDATE_GARDEN
