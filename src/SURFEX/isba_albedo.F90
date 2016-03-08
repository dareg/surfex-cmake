!     #########
      SUBROUTINE ISBA_ALBEDO(HSNOW, OTR_ML, OMEB, PDIR_SW, PSCA_SW, PSW_BANDS, KSW, &
                             IMT, IMA, IR, PFALB, PFFV, PFFG, PGLOBAL_SW,           &
                             PMEB_SCA_SW, PALBNIR_TVEG, PALBVIS_TVEG,               &
                             PALBNIR_TSOIL, PALBVIS_TSOIL               )
!     ##########################################################################
!
!!****  *ISBA_ALBEDO*  
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
!!      Original    
!!      P. Samuelsson  02/2012  MEB
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_TIME_t, ISBA_PARAM_ALB_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
!
USE MODD_SURF_PAR,     ONLY : XUNDEF
!
USE MODI_ALBEDO_FROM_NIR_VIS
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
 CHARACTER(LEN=*)    , INTENT(IN)   :: HSNOW      ! ISBA snow scheme
LOGICAL,              INTENT(IN)   :: OTR_ML
LOGICAL,              INTENT(IN)   :: OMEB        ! True = patch with multi-energy balance 
!                                                 ! False = patch with classical ISBA
!
REAL, DIMENSION(:,:), INTENT(IN)   :: PDIR_SW            ! direct incoming solar radiation
REAL, DIMENSION(:,:), INTENT(IN)   :: PSCA_SW            ! diffus incoming solar radiation
REAL, DIMENSION(:)  , INTENT(IN)   :: PSW_BANDS          ! mean wavelength of each shortwave band (m)
INTEGER,              INTENT(IN)   :: KSW                ! number of short-wave spectral bands
!
TYPE(ISBA_PARAM_TIME_t), INTENT(INOUT) :: IMT
TYPE(ISBA_PARAM_ALB_t), INTENT(INOUT) :: IMA
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
!
REAL, DIMENSION(:)  , INTENT(IN)   :: PFALB              ! Floodplain albedo
REAL, DIMENSION(:)  , INTENT(IN)   :: PFFV               ! Floodplain fraction over vegetation
REAL, DIMENSION(:)  , INTENT(IN)   :: PFFG               ! Floodplain fraction over the ground
!
REAL, DIMENSION(:)  , INTENT(OUT)  :: PGLOBAL_SW         ! global incoming SW rad.
REAL, DIMENSION(:)  , INTENT(OUT)  :: PMEB_SCA_SW        ! diffuse incoming SW rad.
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALBNIR_TVEG       ! nearIR  veg tot albedo
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALBVIS_TVEG       ! visible veg tot albedo
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALBNIR_TSOIL      ! nearIR  soil tot albedo
REAL, DIMENSION(:)  , INTENT(OUT)  :: PALBVIS_TSOIL      ! visible soil tot albedo
!
!-------------------------------------------------------------------------------
!
!*      0.     Local variables
!              ---------------
!
INTEGER                          :: JLAYER
INTEGER                          :: JSWB
REAL, DIMENSION(SIZE(IMT%XALBNIR))      :: ZSW_UP
REAL, DIMENSION(SIZE(IMT%XALBNIR),KSW)  :: ZDIR_ALB_WITHOUT_SNOW
REAL, DIMENSION(SIZE(IMT%XALBNIR),KSW)  :: ZSCA_ALB_WITHOUT_SNOW
REAL, DIMENSION(SIZE(IMT%XALBNIR),KSW)  :: ZDIR_ALB_VEG_WITHOUT_SNOW
REAL, DIMENSION(SIZE(IMT%XALBNIR),KSW)  :: ZSCA_ALB_VEG_WITHOUT_SNOW
REAL, DIMENSION(SIZE(IMT%XALBNIR),KSW)  :: ZDIR_ALB_SOIL_WITHOUT_SNOW
REAL, DIMENSION(SIZE(IMT%XALBNIR),KSW)  :: ZSCA_ALB_SOIL_WITHOUT_SNOW
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*      2.     Compute snow-free albedo
!              ------------------------
!
!* Snow-free surface albedo for each wavelength
!
IF (LHOOK) CALL DR_HOOK('ISBA_ALBEDO',0,ZHOOK_HANDLE)
!
IF (OTR_ML )THEN
   IF (OMEB) THEN
      PALBNIR_TVEG (:) =               IMT%XALBNIR_VEG(:,1)
      PALBNIR_TSOIL(:) = ( 1.-PFFG(:))*IMA%XALBNIR_SOIL(:,1) + PFFG(:)*PFALB(:)   
      PALBVIS_TVEG (:) =               IMT%XALBVIS_VEG(:,1)
      PALBVIS_TSOIL(:) = ( 1.-PFFG(:))*IMA%XALBVIS_SOIL(:,1) + PFFG(:)*PFALB(:)
   ELSE
      PALBNIR_TVEG (:) = IMT%XALBNIR_VEG(:,1)
      PALBNIR_TSOIL(:) = IMA%XALBNIR_SOIL(:,1) 
      PALBVIS_TVEG (:) = IMT%XALBVIS_VEG(:,1)
      PALBVIS_TSOIL(:) = IMA%XALBVIS_SOIL(:,1) 
   ENDIF
ELSE
  PALBNIR_TVEG (:) = XUNDEF
  PALBNIR_TSOIL(:) = XUNDEF
  PALBVIS_TVEG (:) = XUNDEF
  PALBVIS_TSOIL(:) = XUNDEF
ENDIF
!
 CALL ALBEDO_FROM_NIR_VIS(PSW_BANDS, IMT%XALBNIR(:,1), IMT%XALBVIS(:,1), IMT%XALBUV(:,1),  &
                           ZDIR_ALB_WITHOUT_SNOW, ZSCA_ALB_WITHOUT_SNOW )  
!
!* total shortwave incoming radiation
!
  PGLOBAL_SW(:) = 0.
  PMEB_SCA_SW(:) = 0.
  DO JSWB=1,KSW
    PGLOBAL_SW(:) = PGLOBAL_SW(:) + (PDIR_SW(:,JSWB) + PSCA_SW(:,JSWB))
    PMEB_SCA_SW(:) = PMEB_SCA_SW(:) + (PSCA_SW(:,JSWB))
  END DO
!
!* snow-free global albedo (needed by ISBA)
!
  ZSW_UP(:) = 0. 
  DO JSWB=1,KSW
    ZSW_UP(:) =  ZSW_UP(:)                                       &
                 + ZDIR_ALB_WITHOUT_SNOW(:,JSWB) * PDIR_SW(:,JSWB) &
                 + ZSCA_ALB_WITHOUT_SNOW(:,JSWB) * PSCA_SW(:,JSWB)  
  END DO
  IR%XSNOWFREE_ALB(:,1) = XUNDEF
  WHERE(PGLOBAL_SW(:)>0.)  
       IR%XSNOWFREE_ALB(:,1) = ZSW_UP(:) / PGLOBAL_SW(:)
  ELSEWHERE
       IR%XSNOWFREE_ALB(:,1) = ZDIR_ALB_WITHOUT_SNOW(:,1)
  END WHERE
!
  IF(HSNOW == 'EBA') THEN
     CALL ALBEDO_FROM_NIR_VIS(PSW_BANDS,            &
               IMT%XALBNIR_VEG(:,1), IMT%XALBVIS_VEG(:,1), IMT%XALBUV_VEG(:,1), &
               ZDIR_ALB_VEG_WITHOUT_SNOW,            &
               ZSCA_ALB_VEG_WITHOUT_SNOW             )  
     ZSW_UP(:) = 0.
     DO JSWB=1,KSW
        ZSW_UP(:) =  ZSW_UP(:)                                           &
                     + ZDIR_ALB_VEG_WITHOUT_SNOW(:,JSWB) * PDIR_SW(:,JSWB) &
                     + ZSCA_ALB_VEG_WITHOUT_SNOW(:,JSWB) * PSCA_SW(:,JSWB)  
     END DO
     IR%XSNOWFREE_ALB_VEG(:,1) = XUNDEF
     WHERE(PGLOBAL_SW(:)>0.)  IR%XSNOWFREE_ALB_VEG(:,1) = ZSW_UP(:) / PGLOBAL_SW(:)
!
     CALL ALBEDO_FROM_NIR_VIS(PSW_BANDS,               &
               IMA%XALBNIR_SOIL(:,1), IMA%XALBVIS_SOIL(:,1), IMA%XALBUV_SOIL(:,1), &
               ZDIR_ALB_SOIL_WITHOUT_SNOW,              &
               ZSCA_ALB_SOIL_WITHOUT_SNOW               )  
     ZSW_UP(:) = 0.
     DO JSWB=1,KSW
        ZSW_UP(:) =  ZSW_UP(:)                                            &
                     + ZDIR_ALB_SOIL_WITHOUT_SNOW(:,JSWB) * PDIR_SW(:,JSWB) &
                     + ZSCA_ALB_SOIL_WITHOUT_SNOW(:,JSWB) * PSCA_SW(:,JSWB)  
     END DO
     IR%XSNOWFREE_ALB_SOIL(:,1) = XUNDEF
     WHERE(PGLOBAL_SW(:)>0.)  IR%XSNOWFREE_ALB_SOIL(:,1) = ZSW_UP(:) / PGLOBAL_SW(:)             
  ENDIF
IF (LHOOK) CALL DR_HOOK('ISBA_ALBEDO',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ISBA_ALBEDO
