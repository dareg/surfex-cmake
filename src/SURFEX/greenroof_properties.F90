!     #########
      SUBROUTINE GREENROOF_PROPERTIES(PDIR_SW, PSCA_SW, PSW_BANDS, KSW,&
                                      PTS, PEMIS, PALB, PTA,           &  
                                      PALBNIR_TVEG, PALBVIS_TVEG,      &
                                      PALBNIR_TSOIL, PALBVIS_TSOIL     )  
!     ##########################################################################
!
!!****  *GREENROOF_PROPERTIES*  
!!
!!    PURPOSE
!!    -------
!
!     Based on garden_properties
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
!!      S. Belair           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original                      ?     
!!      C. de Munck and A. Lemonsu    09/2011  
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_TEB_n, ONLY : T => TEB
USE MODD_TEB_VEG_n, ONLY : TVG => TEB_VEG_OPTIONS
!
USE MODD_SURF_PAR,            ONLY : XUNDEF
!
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TGRO => TEB_GREENROOF_OPTIONS
!
USE MODD_TEB_GREENROOF_PGD_n, ONLY : TGRP => TEB_GREENROOF_PGD
USE MODD_TEB_GREENROOF_PGD_EVOL_n, ONLY : TGRPE => TEB_GREENROOF_PGD_EVOL
USE MODD_TEB_GREENROOF_n, ONLY : TGR => TEB_GREENROOF
!
USE MODI_ISBA_PROPERTIES
USE MODI_FLAG_TEB_GREENROOF_n
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
REAL, DIMENSION(:)  , INTENT(IN), OPTIONAL :: PTA        ! Air temperature (K)
!
REAL, DIMENSION(:)  , INTENT(OUT), OPTIONAL  :: PALBNIR_TVEG       ! nearIR  veg tot albedo
REAL, DIMENSION(:)  , INTENT(OUT), OPTIONAL  :: PALBVIS_TVEG       ! visible veg tot albedo
REAL, DIMENSION(:)  , INTENT(OUT), OPTIONAL  :: PALBNIR_TSOIL      ! nearIR  soil tot albedo
REAL, DIMENSION(:)  , INTENT(OUT), OPTIONAL  :: PALBVIS_TSOIL      ! visible soil tot albedo
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
REAL, DIMENSION(SIZE(PALB))    :: ZALBNIR_TVEG       ! nearIR  veg tot albedo
REAL, DIMENSION(SIZE(PALB))    :: ZALBVIS_TVEG       ! visible veg tot albedo
REAL, DIMENSION(SIZE(PALB))    :: ZALBNIR_TSOIL      ! nearIR  soil tot albedo
REAL, DIMENSION(SIZE(PALB))    :: ZALBVIS_TSOIL      ! visible soil tot albedo
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!* only one patch for green roofs
IF (LHOOK) CALL DR_HOOK('GREENROOF_PROPERTIES',0,ZHOOK_HANDLE)
!
!*      1.     Set physical values for points where there is no green roof
!              -----------------------------------------------------------
!
! This way, ISBA can run without problem for these points
!
 CALL FLAG_TEB_GREENROOF_n(TGR, TGRO, TGRPE, T, TVG, &
                           1)
!
!
!*      2.     Computes several properties of green roofs
!              ------------------------------------------
!
!
 CALL ISBA_PROPERTIES(TGRO%CISBA_GR, TGRO%LTR_ML_GR, TGR%CUR%TSNOW, 1,             &
                      PDIR_SW, PSCA_SW, PSW_BANDS, KSW,          &
                      TGRPE%CUR%XALBNIR, TGRPE%CUR%XALBVIS, TGRPE%CUR%XALBUV,                  &
                      TGRP%XALBNIR_VEG, TGRP%XALBVIS_VEG, TGRP%XALBUV_VEG,      &
                      TGRP%XALBNIR_SOIL, TGRP%XALBVIS_SOIL,                &
                      TGRP%XALBUV_SOIL,                               &
                      TGRPE%CUR%XVEG, TGRPE%CUR%XLAI, TGRPE%CUR%XZ0, TGRPE%CUR%XEMIS, TGR%CUR%XTG(:,1),          &
                      ZASNOW,ZANOSNOW,                           &
                      ZESNOW,ZENOSNOW,                           &
                      ZTSSNOW,ZTSNOSNOW,                         &
                      TGR%CUR%XSNOWFREE_ALB_VEG, TGR%CUR%XSNOWFREE_ALB_SOIL,     &
                      ZALBNIR_TVEG, ZALBVIS_TVEG, ZALBNIR_TSOIL, &
                      ZALBVIS_TSOIL,                             &
                      TGR%CUR%XPSN, TGR%CUR%XPSNV_A, TGR%CUR%XPSNG, TGR%CUR%XPSNV                )
!
TGR%CUR%XSNOWFREE_ALB = ZANOSNOW
!
!* averaged albedo
PALB =  TGR%CUR%XPSN(:)    * ZASNOW              + (1.-TGR%CUR%XPSN(:)) * ZANOSNOW
!* averaged emissivity
PEMIS=  TGR%CUR%XPSN(:)    * ZESNOW              + (1.-TGR%CUR%XPSN(:)) * ZENOSNOW
!* averaged surface radiative temperature
!  (recomputed from emitted long wave)
PTS  =((TGR%CUR%XPSN(:)    * ZESNOW * ZTSSNOW**4 + (1.-TGR%CUR%XPSN(:)) * ZENOSNOW * ZTSNOSNOW**4) / PEMIS)**0.25
!
IF(PRESENT(PALBNIR_TVEG))PALBNIR_TVEG(:)=ZALBNIR_TVEG(:)
IF(PRESENT(PALBVIS_TVEG))PALBVIS_TVEG(:)=ZALBVIS_TVEG(:)
IF(PRESENT(PALBNIR_TSOIL))PALBNIR_TSOIL(:)=ZALBNIR_TSOIL(:)
IF(PRESENT(PALBVIS_TSOIL))PALBVIS_TSOIL(:)=ZALBVIS_TSOIL(:)
!
IF (LHOOK) CALL DR_HOOK('GREENROOF_PROPERTIES',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE GREENROOF_PROPERTIES
