!     ####################################################################
      SUBROUTINE ALBEDO(HALBEDO, IMT, IMA, PVEGTYPE, OMASK    )  
!     ####################################################################
!
!!****  *ALBEDO*  
!!
!!    PURPOSE
!!    -------
!  computes the albedo of for different types (patches) 
! of natural continental parts, from
! vegetation albedo and soil albedo.
! Soil albedo is estimated from sand fraction.
! A correction due to the soil humidity is used.
!
!     
!!**  METHOD
!!    ------
!
!!    EXTERNAL
!!    --------
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
!!      F.Solmon  /  V. Masson          
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     
!!                  01/2004  Externalization (V. Masson)
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_TIME_t, ISBA_PARAM_ALB_t
!
USE MODD_DATA_COVER_PAR, ONLY : NVT_SNOW
USE MODD_SNOW_PAR,       ONLY : XANSMAX
USE MODD_SURF_PAR,       ONLY : XUNDEF
!
USE MODI_VEGTYPE_TO_PATCH
USE MODI_SURF_PATCH
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=*),       INTENT(IN)   :: HALBEDO
! Albedo dependance wxith surface soil water content
!   "EVOL" = albedo evolves with soil wetness
!   "DRY " = constant albedo value for dry soil
!   "WET " = constant albedo value for wet soil
!   "MEAN" = constant albedo value for medium soil wetness
!
TYPE(ISBA_PARAM_TIME_t), INTENT(INOUT) :: IMT
TYPE(ISBA_PARAM_ALB_t), INTENT(INOUT) :: IMA
!
REAL,    DIMENSION(:,:), INTENT(IN), OPTIONAL :: PVEGTYPE ! vegetation type
LOGICAL, DIMENSION(:),   INTENT(IN), OPTIONAL :: OMASK    ! mask where computations are done
!
!*      0.2    declarations of local variables
!              -------------------------------
!
LOGICAL, DIMENSION(SIZE(IMT%XVEG,1)) :: GMASK
!
REAL, DIMENSION(SIZE(IMT%XVEG,1),SIZE(IMT%XVEG,2))  ::ZPATCH, ZSNOWPATCH 
INTEGER :: ISNOWPATCH !patch index for snow 
INTEGER :: IPATCH     ! number of patches
INTEGER :: JPATCH     !loop index for patches
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ALBEDO',0,ZHOOK_HANDLE)
IF (HALBEDO=='USER' .AND. LHOOK) CALL DR_HOOK('ALBEDO',1,ZHOOK_HANDLE)
IF (HALBEDO=='USER') RETURN
!
GMASK(:) = .TRUE.
IF (PRESENT(OMASK)) GMASK(:) = OMASK(:)
!
IPATCH = SIZE(IMT%XVEG,2)

DO JPATCH=1,IPATCH
  WHERE (GMASK(:))
    IMT%XALBVIS (:,JPATCH) = XUNDEF
    IMT%XALBNIR (:,JPATCH) = XUNDEF
    IMT%XALBUV  (:,JPATCH) = XUNDEF
  END WHERE
END DO
!
!
!
ZSNOWPATCH (:,:) =0.
!
IF (PRESENT(PVEGTYPE)) THEN
  ! calculation of patch surfaces  (weights for average)
  CALL SURF_PATCH(IPATCH,PVEGTYPE,ZPATCH)
  ! permanent snow fraction in the corresponding patch
  ISNOWPATCH= VEGTYPE_TO_PATCH (NVT_SNOW,IPATCH)
  WHERE(GMASK(:) .AND. ZPATCH(:,ISNOWPATCH)>0.)
    ZSNOWPATCH (:,ISNOWPATCH) = PVEGTYPE(:,NVT_SNOW)/ZPATCH(:,ISNOWPATCH)
  END WHERE
END IF
!
DO JPATCH=1,IPATCH
  WHERE (GMASK(:) .AND. IMT%XVEG(:,JPATCH)/=XUNDEF)

    IMT%XALBVIS(:,JPATCH) = ( (1.-IMT%XVEG(:,JPATCH)) * IMA%XALBVIS_SOIL(:,JPATCH)  &
                                + IMT%XVEG(:,JPATCH)  * IMT%XALBVIS_VEG (:,JPATCH)) &
           * (1-ZSNOWPATCH(:,JPATCH)) + XANSMAX  * ZSNOWPATCH(:,JPATCH)   
    !
    IMT%XALBNIR(:,JPATCH) = ( (1.-IMT%XVEG(:,JPATCH)) * IMA%XALBNIR_SOIL(:,JPATCH)  &
                                + IMT%XVEG(:,JPATCH)  * IMT%XALBNIR_VEG (:,JPATCH)) &
           * (1-ZSNOWPATCH(:,JPATCH)) + XANSMAX  * ZSNOWPATCH(:,JPATCH)   
    !
    IMT%XALBUV (:,JPATCH) = ( (1.-IMT%XVEG(:,JPATCH)) * IMA%XALBUV_SOIL (:,JPATCH)  &
                                + IMT%XVEG(:,JPATCH)  * IMT%XALBUV_VEG  (:,JPATCH)) &
           * (1-ZSNOWPATCH(:,JPATCH)) + XANSMAX  * ZSNOWPATCH(:,JPATCH)   
  END WHERE
END DO
IF (LHOOK) CALL DR_HOOK('ALBEDO',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE ALBEDO
