!     #########
    SUBROUTINE VEGETATION_UPDATE_GREENROOF (TGRO, TGRPE, TGRP, T, TOP, TVG, &
                                            TPTIME,PTSTEP,KLU)  
!   ##########################################################################
!
!!****  *GREENROOF*  
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
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
USE MODD_DATA_ISBA_n, ONLY : DTI => DATA_ISBA
USE MODD_DATA_TEB_GARDEN_n, ONLY : DTGD => DATA_TEB_GARDEN
USE MODD_DATA_TEB_GREENROOF_n, ONLY : DTGR => DATA_TEB_GREENROOF
USE MODD_ISBA_GRID_n, ONLY : IG => ISBA_GRID
USE MODD_ISBA_n, ONLY : I => ISBA
!
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TEB_GREENROOF_OPTIONS_t
USE MODD_TEB_GREENROOF_PGD_EVOL_n, ONLY : TEB_GREENROOF_PGD_EVOL_t
USE MODD_TEB_GREENROOF_PGD_n, ONLY : TEB_GREENROOF_PGD_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_OPTIONS_t
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
TYPE(TEB_GREENROOF_OPTIONS_t), INTENT(INOUT) :: TGRO
TYPE(TEB_GREENROOF_PGD_EVOL_t), INTENT(INOUT) :: TGRPE
TYPE(TEB_GREENROOF_PGD_t), INTENT(INOUT) :: TGRP
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_VEG_OPTIONS_t), INTENT(INOUT) :: TVG
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
REAL, DIMENSION(KLU)   :: ZAOSIP
REAL, DIMENSION(KLU)   :: ZAOSIM
REAL, DIMENSION(KLU)   :: ZAOSJP
REAL, DIMENSION(KLU)   :: ZAOSJM
REAL, DIMENSION(KLU)   :: ZHO2IP
REAL, DIMENSION(KLU)   :: ZHO2IM
REAL, DIMENSION(KLU)   :: ZHO2JP
REAL, DIMENSION(KLU)   :: ZHO2JM
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
REAL, DIMENSION(KLU,1) :: ZDMAX
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
! MEB stuff
REAL, DIMENSION(KLU,1) :: ZGNDLITTER
REAL, DIMENSION(KLU,1) :: ZZF_TALLVEG
REAL, DIMENSION(KLU,1) :: ZRGLGV
REAL, DIMENSION(KLU,1) :: ZGAMMAGV
REAL, DIMENSION(KLU,1) :: ZRSMINGV
REAL, DIMENSION(KLU,1) :: ZWRMAX_CFGV
REAL, DIMENSION(KLU,1) :: ZH_VEG
REAL, DIMENSION(KLU,1) :: ZLAIGV
REAL, DIMENSION(KLU,1) :: ZZ0LITTER
!
TYPE (DATE_TIME), DIMENSION(KLU,1) :: TZSEED
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
IF (LHOOK) CALL DR_HOOK('VEGETATION_UPDATE_GREENROF',0,ZHOOK_HANDLE)
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
ZVEG(:,1)         = TGRPE%CUR%XVEG
ZZ0(:,1)          = TGRPE%CUR%XZ0
ZALBNIR(:,1)      = TGRPE%CUR%XALBNIR
ZALBVIS(:,1)      = TGRPE%CUR%XALBVIS
ZALBUV(:,1)       = TGRPE%CUR%XALBUV
ZEMIS(:,1)        = TGRPE%CUR%XEMIS
ZRSMIN(:,1)       = TGRP%XRSMIN
ZGAMMA(:,1)       = TGRP%XGAMMA
ZWRMAX_CF(:,1)    = TGRP%XWRMAX_CF
ZRGL(:,1)         = TGRP%XRGL
ZCV(:,1)          = TGRP%XCV
ZGMES(:,1)        = TGRP%XGMES
ZBSLAI(:,1)       = TGRP%XBSLAI
ZLAIMIN(:,1)      = TGRP%XLAIMIN
ZSEFOLD(:,1)      = TGRP%XSEFOLD
ZGC(:,1)          = TGRP%XGC
ZDMAX(:,1)        = TGRP%XDMAX
ZF2I(:,1)         = TGRP%XF2I
GSTRESS(:,1)      = TGRP%LSTRESS
ZALBNIR_VEG(:,1)  = TGRP%XALBNIR_VEG
ZALBVIS_VEG(:,1)  = TGRP%XALBVIS_VEG
ZALBUV_VEG(:,1)   = TGRP%XALBUV_VEG
ZALBNIR_SOIL(:,1) = TGRP%XALBNIR_SOIL
ZALBVIS_SOIL(:,1) = TGRP%XALBVIS_SOIL
ZALBUV_SOIL(:,1)  = TGRP%XALBUV_SOIL
ZCE_NITRO(:,1)    = TGRP%XCE_NITRO
ZCF_NITRO(:,1)    = TGRP%XCF_NITRO
ZCNA_NITRO(:,1)   = TGRP%XCNA_NITRO
! --------------------------------------------------------------------------------------
! Vegetation update (in case of non-interactive vegetation):
! --------------------------------------------------------------------------------------
!
GUPDATED=.FALSE.
IF (TVG%CPHOTO=='NON' .OR. TVG%CPHOTO=='AGS' .OR. TVG%CPHOTO=='AST') THEN
     CALL VEGETATION_UPDATE(DTCO, DTI, DTGD, DTGR, IG, I, TGRO, &
                            PTSTEP,TPTIME,TOP%XCOVER,TOP%LCOVER,                 &
                         TVG%CISBA,(.NOT. TGRO%LPAR_GREENROOF), TVG%CPHOTO, .FALSE.,  &
                         TVG%LTR_ML, 'GR ',                                  &
                         ZLAI,ZVEG,ZZ0,                                  &
                         ZALBNIR,ZALBVIS,ZALBUV,ZEMIS,                   &
                         ZRSMIN,ZGAMMA,ZWRMAX_CF,                        &
                         ZRGL,ZCV,                                       &
                         ZGMES,ZBSLAI,ZLAIMIN,ZSEFOLD,ZGC,ZDMAX,         &
                         ZF2I, GSTRESS,                                  &
                         ZAOSIP,ZAOSIM,ZAOSJP,ZAOSJM,                    &
                         ZHO2IP,ZHO2IM,ZHO2JP,ZHO2JM,                    &
                         ZZ0EFFIP,ZZ0EFFIM,ZZ0EFFJP,ZZ0EFFJM,            &
                         TVG%CALBEDO, ZALBNIR_VEG, ZALBVIS_VEG, ZALBUV_VEG,  &
                         ZALBNIR_SOIL, ZALBVIS_SOIL, ZALBUV_SOIL,        &
                         ZCE_NITRO, ZCF_NITRO, ZCNA_NITRO,               &
                         TZSEED, TZREAP, ZWATSUP, ZIRRIG,                &
                         ZGNDLITTER,ZZF_TALLVEG, ZRGLGV,ZGAMMAGV,        &
                         ZRSMINGV, ZWRMAX_CFGV,                          &
                         ZH_VEG, ZLAIGV, ZZ0LITTER,                      &
                         GUPDATED, OABSENT=(T%CUR%XGREENROOF==0.)              )
END IF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
TGRPE%CUR%XVEG         = ZVEG(:,1)
TGRPE%CUR%XZ0          = ZZ0(:,1)
TGRPE%CUR%XALBNIR      = ZALBNIR(:,1)
TGRPE%CUR%XALBVIS      = ZALBVIS(:,1)
TGRPE%CUR%XALBUV       = ZALBUV(:,1)
TGRPE%CUR%XEMIS        = ZEMIS(:,1)
TGRP%XRSMIN       = ZRSMIN(:,1)
TGRP%XGAMMA       = ZGAMMA(:,1)
TGRP%XWRMAX_CF    = ZWRMAX_CF(:,1)
TGRP%XRGL         = ZRGL(:,1)
TGRP%XCV          = ZCV(:,1)
TGRP%XGMES        = ZGMES(:,1)
TGRP%XBSLAI       = ZBSLAI(:,1)
TGRP%XLAIMIN      = ZLAIMIN(:,1)
TGRP%XSEFOLD      = ZSEFOLD(:,1)
TGRP%XGC          = ZGC(:,1)
TGRP%XDMAX        = ZDMAX(:,1)
TGRP%XF2I         = ZF2I(:,1)
TGRP%LSTRESS      = GSTRESS(:,1)
TGRP%XALBNIR_VEG  = ZALBNIR_VEG(:,1)
TGRP%XALBVIS_VEG  = ZALBVIS_VEG(:,1)
TGRP%XALBUV_VEG   = ZALBUV_VEG(:,1)
TGRP%XALBNIR_SOIL = ZALBNIR_SOIL(:,1)
TGRP%XALBVIS_SOIL = ZALBVIS_SOIL(:,1)
TGRP%XALBUV_SOIL  = ZALBUV_SOIL(:,1)
TGRP%XCE_NITRO    = ZCE_NITRO(:,1)
TGRP%XCF_NITRO    = ZCF_NITRO(:,1)
TGRP%XCNA_NITRO   = ZCNA_NITRO(:,1)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
IF (LHOOK) CALL DR_HOOK('VEGETATION_UPDATE_GREENROOF',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE VEGETATION_UPDATE_GREENROOF
