!     #########
      SUBROUTINE GARDEN_PROPERTIES(PDIR_SW, PSCA_SW, PSW_BANDS, KSW,      &
                                   PTS, PEMIS, PALB                       )  
!     ##########################################################################
!
!!****  *GARDEN_PROPERTIES*  
!!
!!    PURPOSE
!!    -------
!
!     Calculates grid-averaged albedo and emissivity (according to snow scheme)
!         
!!    EXTERNAL
!!    --------
!!
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      
!!    AUTHOR
!!    ------
!!
!!	S. Belair           * Meteo-France *
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_SURF_PAR, ONLY : XUNDEF
!
USE MODD_TEB_GARDEN_n,      ONLY : CISBA, TSNOW, XALBNIR, XALBVIS, XALBUV,    &
                                   XALBNIR_VEG, XALBVIS_VEG, XALBUV_VEG,      &
                                   XALBNIR_SOIL, XALBVIS_SOIL, XALBUV_SOIL,   &
                                   XVEG, XLAI, XZ0, XEMIS, XTG, LTR_ML,       &
                                   XPSN, XPSNV, XPSNG, XPSNV_A, XALBNIR_TVEG, &
                                   XALBVIS_TVEG, XALBNIR_TSOIL, XALBVIS_TSOIL  
USE MODD_DIAG_TEB_GARDEN_n, ONLY : XSNOWFREE_ALB_VEG, XSNOWFREE_ALB_SOIL,     &
                                   XSNOWFREE_ALB                                   
!
USE MODI_ISBA_PROPERTIES
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:), INTENT(IN)   :: PDIR_SW            ! direct incoming solar radiation
REAL, DIMENSION(:,:), INTENT(IN)   :: PSCA_SW            ! diffus incoming solar radiation
REAL, DIMENSION(:)  , INTENT(IN)   :: PSW_BANDS          ! mean wavelength of each shortwave band (m)
INTEGER,              INTENT(IN)   :: KSW                ! number of short-wave spectral bands
!
REAL, DIMENSION(:)  , INTENT(OUT)  :: PTS                ! radiative surface temperature
REAL, DIMENSION(:)  , INTENT(OUT)  :: PEMIS              ! green areas emissivity
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALB               ! green areas albedo
!
!-------------------------------------------------------------------------------
!
!*      0.2    Local variables
!              ---------------
!
INTEGER                        :: JLAYER
INTEGER                        :: JSWB
!
REAL, DIMENSION(SIZE(PALB))    :: ZTSNOSNOW ! surf. temp. on snow free part
REAL, DIMENSION(SIZE(PALB))    :: ZTSSNOW   ! surf. temp. on snow covered part
REAL, DIMENSION(SIZE(PALB))    :: ZANOSNOW  ! snow-free surface albedo
REAL, DIMENSION(SIZE(PALB))    :: ZASNOW    ! snow albedo
REAL, DIMENSION(SIZE(PALB))    :: ZENOSNOW  ! snow-free surface emissivity
REAL, DIMENSION(SIZE(PALB))    :: ZESNOW    ! snow emissivity
!
INTEGER  :: JP ! patch number
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!* only one patch for gardens
IF (LHOOK) CALL DR_HOOK('GARDEN_PROPERTIES',0,ZHOOK_HANDLE)
JP = 1
!
CALL ISBA_PROPERTIES(CISBA, LTR_ML, TSNOW, 1,                                     &
                     PDIR_SW, PSCA_SW, PSW_BANDS, KSW,                            &
                     XALBNIR(:,JP), XALBVIS(:,JP), XALBUV(:,JP),                  &
                     XALBNIR_VEG(:,JP), XALBVIS_VEG(:,JP), XALBUV_VEG(:,JP),      &
                     XALBNIR_SOIL(:,JP), XALBVIS_SOIL(:,JP), XALBUV_SOIL(:,JP),   &
                     XVEG(:,JP), XLAI(:,JP), XZ0(:,JP), XEMIS(:,JP),XTG(:,JP,JP), &
                     ZASNOW, ZANOSNOW, ZESNOW, ZENOSNOW, ZTSSNOW, ZTSNOSNOW,      &
                     XSNOWFREE_ALB_VEG, XSNOWFREE_ALB_SOIL,                       &
                     XALBNIR_TVEG, XALBVIS_TVEG, XALBNIR_TSOIL, XALBVIS_TSOIL,    &
                     XPSN(:,JP), XPSNV_A(:,JP), XPSNG(:,JP), XPSNV(:,JP)          )  
!
XSNOWFREE_ALB = ZANOSNOW
!
!* averaged albedo
PALB =  XPSN(:,JP) * ZASNOW              + (1.-XPSN(:,JP)) * ZANOSNOW
!* averaged emissivity
PEMIS=  XPSN(:,JP) * ZESNOW              + (1.-XPSN(:,JP)) * ZENOSNOW
!* averaged surface radiative temperature
!  (recomputed from emitted long wave)
PTS  =((XPSN(:,JP) * ZESNOW * ZTSSNOW**4 + (1.-XPSN(:,JP)) * ZENOSNOW * ZTSNOSNOW**4) / PEMIS)**0.25
IF (LHOOK) CALL DR_HOOK('GARDEN_PROPERTIES',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE GARDEN_PROPERTIES
