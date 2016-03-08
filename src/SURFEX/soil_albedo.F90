!     ####################################################################
      SUBROUTINE SOIL_ALBEDO(HALBEDO, PWSAT, PWG1, IP, IMA, HBAND)  
!     ####################################################################
!
!!****  *SOIL_ALBEDO*  
!!
!!    PURPOSE
!!    -------
!  computes the SOIL_ALBEDO of for different types (patches) 
! of natural continental parts.
!
! Soil SOIL_ALBEDO is estimated from sand fraction.
! A correction due to the soil humidity can be used.
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
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_ALB_t
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
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
! SOIL_ALBEDO dependance wxith surface soil water content
!   "EVOL" = SOIL_ALBEDO evolves with soil wetness
!   "DRY " = constant SOIL_ALBEDO value for dry soil
!   "WET " = constant SOIL_ALBEDO value for wet soil
!   "MEAN" = constant SOIL_ALBEDO value for medium soil wetness
!
REAL, DIMENSION(:), INTENT(IN)    :: PWSAT       ! saturation water content
REAL, DIMENSION(:,:), INTENT(IN)  :: PWG1        ! surface water content
!
!
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
TYPE(ISBA_PARAM_ALB_t), INTENT(INOUT) :: IMA
!
 CHARACTER(LEN=*), INTENT(IN) :: HBAND
!
!*      0.2    declarations of local variables
!              -------------------------------
!
REAL,    DIMENSION(SIZE(PWSAT)) :: ZX
!
INTEGER :: IPATCH     ! number of patches
INTEGER :: JPATCH     !loop index for patches
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SOIL_ALBEDO',0,ZHOOK_HANDLE)
IF (HALBEDO=='USER' .AND. LHOOK) CALL DR_HOOK('SOIL_ALBEDO',1,ZHOOK_HANDLE)
IF (HALBEDO=='USER') RETURN
!
IPATCH = SIZE(PWG1,2)

IF (TRIM(HBAND)=="VIS".OR.TRIM(HBAND)=="ALL") IMA%XALBVIS_SOIL = XUNDEF
IF (TRIM(HBAND)=="NIR".OR.TRIM(HBAND)=="ALL") IMA%XALBNIR_SOIL = XUNDEF
IF (TRIM(HBAND)=="UV" .OR.TRIM(HBAND)=="ALL") IMA%XALBUV_SOIL  = XUNDEF
!
SELECT CASE ( HALBEDO )
CASE ('EVOL')

  DO JPATCH=1,IPATCH
    ZX = MIN( PWG1(:,JPATCH)/PWSAT(:) , 1. )

    IF (TRIM(HBAND)=="VIS".OR.TRIM(HBAND)=="ALL") &
      WHERE (PWG1(:,JPATCH)/=XUNDEF) &
        IMA%XALBVIS_SOIL(:,JPATCH) = IP%XALBVIS_WET(:) + &
           (0.25*IP%XALBVIS_DRY(:)-IP%XALBVIS_WET(:))           &
             * (1. - ZX(:)) * ( ZX(:) + (IP%XALBVIS_DRY(:)-IP%XALBVIS_WET(:)) &
             /(0.25*IP%XALBVIS_DRY(:)-IP%XALBVIS_WET(:)) )  
    IF (TRIM(HBAND)=="NIR".OR.TRIM(HBAND)=="ALL") &
      WHERE (PWG1(:,JPATCH)/=XUNDEF) &      
        IMA%XALBNIR_SOIL(:,JPATCH) = IP%XALBNIR_WET(:) + &
           (0.25*IP%XALBNIR_DRY(:)-IP%XALBNIR_WET(:))           &
             * (1. - ZX(:)) * ( ZX(:) + (IP%XALBNIR_DRY(:)-IP%XALBNIR_WET(:)) &
             /(0.25*IP%XALBNIR_DRY(:)-IP%XALBNIR_WET(:)) )  
    IF (TRIM(HBAND)=="UV".OR.TRIM(HBAND)=="ALL") &
      WHERE (PWG1(:,JPATCH)/=XUNDEF) &      
        IMA%XALBUV_SOIL (:,JPATCH) = IP%XALBUV_WET (:) + &
             (0.25*IP%XALBUV_DRY (:)-IP%XALBUV_WET (:))           &
             * (1. - ZX(:)) * ( ZX(:) + (IP%XALBUV_DRY (:)-IP%XALBUV_WET (:)) &
             /(0.25*IP%XALBUV_DRY (:)-IP%XALBUV_WET (:)) )  

    !END WHERE
  END DO

CASE ('DRY ')
  IF (TRIM(HBAND)=="VIS".OR.TRIM(HBAND)=="ALL") &
          IMA%XALBVIS_SOIL(:,:) = SPREAD(IP%XALBVIS_DRY(:),2,IPATCH)
  IF (TRIM(HBAND)=="NIR".OR.TRIM(HBAND)=="ALL") &
          IMA%XALBNIR_SOIL(:,:) = SPREAD(IP%XALBNIR_DRY(:),2,IPATCH)
  IF (TRIM(HBAND)=="UV".OR.TRIM(HBAND)=="ALL") &
          IMA%XALBUV_SOIL (:,:) = SPREAD(IP%XALBUV_DRY (:),2,IPATCH)

CASE ('WET ')
  IF (TRIM(HBAND)=="VIS".OR.TRIM(HBAND)=="ALL") &
          IMA%XALBVIS_SOIL(:,:) = SPREAD(IP%XALBVIS_WET(:),2,IPATCH)
  IF (TRIM(HBAND)=="NIR".OR.TRIM(HBAND)=="ALL") &
          IMA%XALBNIR_SOIL(:,:) = SPREAD(IP%XALBNIR_WET(:),2,IPATCH)
  IF (TRIM(HBAND)=="UV".OR.TRIM(HBAND)=="ALL") &
          IMA%XALBUV_SOIL (:,:) = SPREAD(IP%XALBUV_WET (:),2,IPATCH)

CASE ('MEAN')
  IF (TRIM(HBAND)=="VIS".OR.TRIM(HBAND)=="ALL") &
          IMA%XALBVIS_SOIL(:,:) = 0.5 * ( SPREAD(IP%XALBVIS_DRY(:),2,IPATCH) + &
                                  SPREAD(IP%XALBVIS_WET(:),2,IPATCH) )
  IF (TRIM(HBAND)=="NIR".OR.TRIM(HBAND)=="ALL") &
          IMA%XALBNIR_SOIL(:,:) = 0.5 * ( SPREAD(IP%XALBNIR_DRY(:),2,IPATCH) + &
                                  SPREAD(IP%XALBNIR_WET(:),2,IPATCH) )
  IF (TRIM(HBAND)=="UV".OR.TRIM(HBAND)=="ALL") &
          IMA%XALBUV_SOIL (:,:) = 0.5 * ( SPREAD(IP%XALBUV_DRY (:),2,IPATCH) + &
                                  SPREAD(IP%XALBUV_WET (:),2,IPATCH) )

END SELECT
IF (LHOOK) CALL DR_HOOK('SOIL_ALBEDO',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SOIL_ALBEDO
