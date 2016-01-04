!#############################################################
SUBROUTINE CONVERT_PATCH_GARDEN (DTCO, DTI, TGDO, TVM, TOP, PPERM, &
                                 KLU,KDECADE)
!#############################################################
!
!!****  *CONVERT_PATCH_GARDEN* - routine to initialize GARDEN parameters 
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
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_TEB_VEG_PARAM_n, ONLY : TEB_VEG_PARAM_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
!
USE MODD_TYPE_DATE_SURF,    ONLY: DATE_TIME
!
USE MODI_CONVERT_PATCH_ISBA
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
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTI
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: TGDO
TYPE(TEB_VEG_PARAM_t), INTENT(INOUT) :: TVM
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
!
REAL, DIMENSION(:), INTENT(IN) :: PPERM
!
INTEGER, INTENT(IN) :: KLU     ! number of points
INTEGER, INTENT(IN) :: KDECADE ! number of decades
!
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL,             DIMENSION(KLU,1)               :: ZLAI
REAL,             DIMENSION(KLU,1)               :: ZVEG
REAL,             DIMENSION(KLU,1)               :: ZZ0
REAL,             DIMENSION(KLU,1)               :: ZEMIS
REAL,             DIMENSION(KLU,1)               :: ZRSMIN
REAL,             DIMENSION(KLU,1)               :: ZGAMMA
REAL,             DIMENSION(KLU,1)               :: ZWRMAX_CF
REAL,             DIMENSION(KLU,1)               :: ZRGL
REAL,             DIMENSION(KLU,1)               :: ZCV
REAL,             DIMENSION(KLU,1)               :: ZGMES
REAL,             DIMENSION(KLU,1)               :: ZBSLAI
REAL,             DIMENSION(KLU,1)               :: ZLAIMIN
REAL,             DIMENSION(KLU,1)               :: ZSEFOLD
REAL,             DIMENSION(KLU,1)               :: ZGC
REAL,             DIMENSION(KLU,1)               :: ZDMAX
REAL,             DIMENSION(KLU,1)               :: ZF2I
REAL,             DIMENSION(KLU,1)               :: ZALBNIR_VEG
REAL,             DIMENSION(KLU,1)               :: ZALBVIS_VEG
REAL,             DIMENSION(KLU,1)               :: ZALBUV_VEG
REAL,             DIMENSION(KLU,1)               :: ZCE_NITRO
REAL,             DIMENSION(KLU,1)               :: ZCF_NITRO
REAL,             DIMENSION(KLU,1)               :: ZCNA_NITRO
REAL,             DIMENSION(KLU,1)               :: ZRE25
REAL,             DIMENSION(KLU,1)               :: ZH_TREE
REAL,             DIMENSION(KLU,1)               :: ZZ0_O_Z0H
REAL,             DIMENSION(KLU,1)               :: ZD_ICE
REAL,             DIMENSION(KLU,TGDO%NGROUND_LAYER,1) :: ZROOTFRAC
REAL,             DIMENSION(KLU,TGDO%NGROUND_LAYER,1) :: ZDG
REAL,             DIMENSION(KLU,1)               :: ZDROOT
REAL,             DIMENSION(KLU,1)               :: ZDG2
INTEGER,          DIMENSION(KLU,1)               :: IWG_LAYER
LOGICAL,          DIMENSION(KLU,1)               :: GSTRESS
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!               Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('CONVERT_PATCH_GARDEN',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!
  CALL CONVERT_PATCH_ISBA(DTCO, DTI, TGDO, &
                          TGDO%CISBA,KDECADE,KDECADE,TOP%XCOVER,TOP%LCOVER,TGDO%CPHOTO,.FALSE.,  &
                        .FALSE.,TGDO%LTR_ML,'GRD',PVEG=ZVEG,PLAI=ZLAI,              &
                        PRSMIN=ZRSMIN,PGAMMA=ZGAMMA,PWRMAX_CF=ZWRMAX_CF,       &
                        PRGL=ZRGL,PCV=ZCV,PSOILGRID=TGDO%XSOILGRID,                 &
                        PDG=ZDG,KWG_LAYER=IWG_LAYER,PDROOT=ZDROOT,PDG2=ZDG2,   &
                        PZ0=ZZ0,PZ0_O_Z0H=ZZ0_O_Z0H,PPERM=PPERM,         &
                        PALBNIR_VEG=ZALBNIR_VEG,PALBVIS_VEG=ZALBVIS_VEG,       &
                        PALBUV_VEG=ZALBUV_VEG,PEMIS_ECO=ZEMIS,                 &
                        PVEGTYPE=TVM%X%XVEGTYPE,PROOTFRAC=ZROOTFRAC,                 &
                        PGMES=ZGMES,PBSLAI=ZBSLAI,PLAIMIN=ZLAIMIN,             &
                        PSEFOLD=ZSEFOLD,PGC=ZGC,                               &
                        PDMAX=ZDMAX,PF2I=ZF2I,OSTRESS=GSTRESS,PH_TREE=ZH_TREE, &
                        PRE25=ZRE25,PCE_NITRO=ZCE_NITRO,PCF_NITRO=ZCF_NITRO,   &
                        PCNA_NITRO=ZCNA_NITRO,PD_ICE=ZD_ICE                    )
!
TVM%T%CUR%XVEG         = ZVEG
TVM%T%CUR%XLAI         = ZLAI
TVM%T%CUR%XZ0          = ZZ0
TVM%T%CUR%XEMIS        = ZEMIS
TVM%T%CUR%XRSMIN       = ZRSMIN
TVM%T%CUR%XGAMMA       = ZGAMMA
TVM%T%CUR%XWRMAX_CF    = ZWRMAX_CF
TVM%T%CUR%XRGL         = ZRGL
TVM%T%CUR%XCV          = ZCV
TVM%T%CUR%XGMES        = ZGMES
TVM%T%CUR%XBSLAI       = ZBSLAI
TVM%T%CUR%XLAIMIN      = ZLAIMIN
TVM%T%CUR%XSEFOLD      = ZSEFOLD
TVM%T%CUR%XGC          = ZGC
TVM%X%XDMAX        = ZDMAX
TVM%T%CUR%XF2I         = ZF2I
TVM%T%CUR%LSTRESS      = GSTRESS
TVM%T%CUR%XALBNIR_VEG  = ZALBNIR_VEG
TVM%T%CUR%XALBVIS_VEG  = ZALBVIS_VEG
TVM%T%CUR%XALBUV_VEG   = ZALBUV_VEG
TVM%T%CUR%XCE_NITRO    = ZCE_NITRO
TVM%T%CUR%XCF_NITRO    = ZCF_NITRO
TVM%T%CUR%XCNA_NITRO   = ZCNA_NITRO
TVM%X%XH_TREE      = ZH_TREE
TVM%X%XRE25        = ZRE25
TVM%X%XROOTFRAC    = ZROOTFRAC
TVM%X%XDG          = ZDG
TVM%X%XDROOT       = ZDROOT
TVM%X%XDG2         = ZDG2
TVM%X%NWG_LAYER    = IWG_LAYER
TVM%X%XZ0_O_Z0H    = ZZ0_O_Z0H
TVM%X%XD_ICE       = ZD_ICE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('CONVERT_PATCH_GARDEN',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE CONVERT_PATCH_GARDEN
