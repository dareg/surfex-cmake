!     #########
      SUBROUTINE ISBA_PROPERTIES(VO, VR, VMT, VMA, KPATCH,   &
                                 PDIR_SW, PSCA_SW, PSW_BANDS, KSW,     &
                                 PASNOW, PANOSNOW, PESNOW, PENOSNOW,      &
                                 PTSSNOW, PTSNOSNOW,                      &
                                 PALBNIR_TVEG, PALBVIS_TVEG,              &
                                 PALBNIR_TSOIL, PALBVIS_TSOIL         )  
!     ##########################################################################
!
!!****  *ISBA_PROPERTIES*  
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
!!
!!    MODIFICATIONS
!!    ------------- 
!!      
!!      P. Samuelsson  02/2012  MEB
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_TIME_t, ISBA_PARAM_ALB_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
!
USE MODD_TYPE_SNOW
USE MODD_SNOW_PAR   , ONLY : XEMISSN, XEMCRIN, XSNOWDMIN, &
                               XRHOSMAX_ES, XRHOSMIN_ES  
USE MODD_WATER_PAR  , ONLY : XEMISWAT
!
USE MODI_ISBA_SNOW_FRAC
USE MODI_ISBA_ALBEDO
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: VO
TYPE(ISBA_PROG_t), INTENT(INOUT) :: VR
TYPE(ISBA_PARAM_TIME_t), INTENT(INOUT) :: VMT
TYPE(ISBA_PARAM_ALB_t), INTENT(INOUT) :: VMA
!
INTEGER,              INTENT(IN)   :: KPATCH     ! patch being treated
!
REAL, DIMENSION(:,:), INTENT(IN)   :: PDIR_SW            ! direct incoming solar radiation
REAL, DIMENSION(:,:), INTENT(IN)   :: PSCA_SW            ! diffus incoming solar radiation
REAL, DIMENSION(:)  , INTENT(IN)   :: PSW_BANDS          ! mean wavelength of each shortwave band (m)
INTEGER,              INTENT(IN)   :: KSW                ! number of short-wave spectral bands            
!
REAL, DIMENSION(:)  , INTENT(OUT)  :: PASNOW    ! = snow albedo
REAL, DIMENSION(:)  , INTENT(OUT)  :: PANOSNOW  ! = snow free albedo 
REAL, DIMENSION(:)  , INTENT(OUT)  :: PESNOW    ! = snow emissivity
REAL, DIMENSION(:)  , INTENT(OUT)  :: PENOSNOW  ! = snow free emissivity
REAL, DIMENSION(:)  , INTENT(OUT)  :: PTSSNOW   ! = snow radiative temperature
REAL, DIMENSION(:)  , INTENT(OUT)  :: PTSNOSNOW ! = snow free radiative temperature
!
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALBNIR_TVEG       ! nearIR  veg tot albedo
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALBVIS_TVEG       ! visible veg tot albedo
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALBNIR_TSOIL      ! nearIR  soil tot albedo
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALBVIS_TSOIL      ! visible soil tot albedo
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PDIR_SW,1)) :: ZGLOBAL_SW                 ! global incoming SW rad.
REAL, DIMENSION(SIZE(VMT%XALBNIR))   :: ZALBF
REAL, DIMENSION(SIZE(VMT%XALBNIR))   :: ZFFV
REAL, DIMENSION(SIZE(VMT%XALBNIR))   :: ZFFG
!
LOGICAL, PARAMETER :: GMEB=.FALSE.
REAL, DIMENSION(SIZE(PDIR_SW,1))   :: ZP_MEB_SCA_SW, ZALBNIR_TSNOW, ZALBVIS_TSNOW
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ISBA_PROPERTIES',0,ZHOOK_HANDLE)
!
 CALL ISBA_SNOW_FRAC(VR%TSNOW%SCHEME, VR%TSNOW%WSNOW(:,:,KPATCH),  &
                     VR%TSNOW%RHO(:,:,KPATCH), VR%TSNOW%ALB(:,KPATCH), &
                     VMT%XVEG(:,1), VMT%XLAI(:,1), VMT%XZ0(:,1),       &
                     VR%XPSN(:,1), VR%XPSNV_A(:,1), VR%XPSNG(:,1), VR%XPSNV(:,1)  )  
!
!-------------------------------------------------------------------------------
!*      2.     Compute snow-free albedo
!              ------------------------
!
!* Snow-free surface albedo for each wavelength
!
ZALBF         = 0.
ZFFV          = 0.
ZFFG          = 0.
!
 CALL ISBA_ALBEDO(VR%TSNOW%SCHEME, VO%LTR_ML, GMEB, PDIR_SW, PSCA_SW, PSW_BANDS, KSW, &
                  VMT, VMA, VR, ZALBF, ZFFV, ZFFG, ZGLOBAL_SW, ZP_MEB_SCA_SW, &
                  PALBNIR_TVEG, PALBVIS_TVEG, PALBNIR_TSOIL, PALBVIS_TSOIL)

PANOSNOW(:) = VR%XSNOWFREE_ALB(:,1)
!-------------------------------------------------------------------------------
!
!*      3.     Compute aggeragted albedo and emissivity
!              ----------------------------------------
!
IF(VR%TSNOW%SCHEME == '3-L' .OR. VR%TSNOW%SCHEME == 'CRO' .OR. VO%CISBA == 'DIF')THEN
!
! NON-SNOW covered Grid averaged albedo and emissivity for explicit snow scheme:
!
  PASNOW(:) = VR%TSNOW%ALB(:,KPATCH)
  PESNOW(:) = VR%TSNOW%EMIS(:,KPATCH)
  PENOSNOW(:) = VMT%XEMIS(:,1)

  PTSSNOW(:)   = VR%TSNOW%TS(:,KPATCH)
  PTSNOSNOW(:) = VR%XTG(:,1,1)

ELSE
!
! Grid averaged albedo and emissivity for composite snow scheme:
!
  IF(VR%TSNOW%SCHEME =='EBA') THEN
!
    PASNOW(:) = VR%TSNOW%ALB(:,KPATCH)
    PESNOW(:) = XEMCRIN
    PENOSNOW(:) = VMT%XEMIS(:,1)

    PTSSNOW(:)   = VR%XTG(:,1,1)
    PTSNOSNOW(:) = VR%XTG(:,1,1)

  ELSE

    PASNOW(:) = VR%TSNOW%ALB(:,KPATCH)
    PESNOW(:) = XEMISSN
    PENOSNOW(:) = VMT%XEMIS(:,1)

    PTSSNOW(:)   = VR%XTG(:,1,1)
    PTSNOSNOW(:) = VR%XTG(:,1,1)

  ENDIF
!
ENDIF
IF (LHOOK) CALL DR_HOOK('ISBA_PROPERTIES',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ISBA_PROPERTIES
