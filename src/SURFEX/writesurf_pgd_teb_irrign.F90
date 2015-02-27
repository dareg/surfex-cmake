!     #########
      SUBROUTINE WRITESURF_PGD_TEB_IRRIG_n(HPROGRAM)
!     ################################################
!
!!****  *WRITESURF_PGD_TEB_IRRIG_n* - writes TEB irrigation physiographic fields
!!                        
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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/2013
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_TEB_IRRIG_n, ONLY : LPAR_GD_IRRIG, LPAR_GR_IRRIG, LPAR_RD_IRRIG,                                 &
                             XGD_START_MONTH, XGD_END_MONTH, XGD_START_HOUR, XGD_END_HOUR, XGD_24H_IRRIG, &
                             XRD_START_MONTH, XRD_END_MONTH, XRD_START_HOUR, XRD_END_HOUR, XRD_24H_IRRIG, &
                             XGR_START_MONTH, XGR_END_MONTH, XGR_START_HOUR, XGR_END_HOUR, XGR_24H_IRRIG
!
USE MODI_WRITE_SURF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
 CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
 CHARACTER(LEN=100):: YCOMMENT       ! Comment string
INTEGER           :: JLAYER         ! loop index
INTEGER           :: JTIME          ! loop index
REAL, DIMENSION(:), ALLOCATABLE :: ZWORK
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_PGD_TEB_IRRIG_n',0,ZHOOK_HANDLE)
!
! Flag for irrigation of gardens
YRECFM='L_PAR_GD_IRR'
YCOMMENT='FLAG FOR SPECIFIED GARDEN IRRIGATION PARAMETERS'
 CALL WRITE_SURF(HPROGRAM,YRECFM,LPAR_GD_IRRIG,IRESP,HCOMMENT=YCOMMENT)
!
! Parameters describing irrigation
IF (LPAR_GD_IRRIG) THEN
!
  YRECFM='D_GD_SM_IRR'
  YCOMMENT='Start Month for Gardens Irrigation'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XGD_START_MONTH(:),IRESP,HCOMMENT=YCOMMENT)
!
  YRECFM='D_GD_EM_IRR'
  YCOMMENT='End   Month for Gardens Irrigation'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XGD_END_MONTH  (:),IRESP,HCOMMENT=YCOMMENT)
!
  YRECFM='D_GD_SH_IRR'
  YCOMMENT='Start Hour  for Gardens Irrigation'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XGD_START_HOUR (:),IRESP,HCOMMENT=YCOMMENT)
!
  YRECFM='D_GD_EH_IRR'
  YCOMMENT='End   Hour  for Gardens Irrigation'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XGD_END_HOUR   (:),IRESP,HCOMMENT=YCOMMENT)
!
  YRECFM='D_GD_IRRIG'
  YCOMMENT='24h mean Irrigation rate for Gardens Irrigation'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XGD_24H_IRRIG  (:),IRESP,HCOMMENT=YCOMMENT)
!
ENDIF
!
! Flag for irrigation of greenroofs
YRECFM='L_PAR_GR_IRR'
YCOMMENT='FLAG FOR SPECIFIED GREENROOFS IRRIGATION PARAMETERS'
 CALL WRITE_SURF(HPROGRAM,YRECFM,LPAR_GR_IRRIG,IRESP,HCOMMENT=YCOMMENT)
!
! Parameters describing irrigation
IF (LPAR_GR_IRRIG) THEN
!
  YRECFM='D_GR_SM_IRR'
  YCOMMENT='Start Month for Greenroofs Irrigation'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XGR_START_MONTH(:),IRESP,HCOMMENT=YCOMMENT)
!
  YRECFM='D_GR_EM_IRR'
  YCOMMENT='End   Month for Greenroofs Irrigation'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XGR_END_MONTH  (:),IRESP,HCOMMENT=YCOMMENT)
!
  YRECFM='D_GR_SH_IRR'
  YCOMMENT='Start Hour  for Greenroofs Irrigation'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XGR_START_HOUR (:),IRESP,HCOMMENT=YCOMMENT)
!
  YRECFM='D_GR_EH_IRR'
  YCOMMENT='End   Hour  for Greenroofs Irrigation'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XGR_END_HOUR   (:),IRESP,HCOMMENT=YCOMMENT)
!
  YRECFM='D_GR_IRRIG'
  YCOMMENT='24h mean Irrigation rate for Greenroofs Irrigation'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XGR_24H_IRRIG  (:),IRESP,HCOMMENT=YCOMMENT)
!
ENDIF
!
! Flag for watering of greenroofs
YRECFM='L_PAR_RD_IRR'
YCOMMENT='FLAG FOR SPECIFIED ROAD IRRIGATION PARAMETERS'
 CALL WRITE_SURF(HPROGRAM,YRECFM,LPAR_RD_IRRIG,IRESP,HCOMMENT=YCOMMENT)
!
! Parameters describing watering
IF (LPAR_RD_IRRIG) THEN
!
  YRECFM='D_RD_SM_IRR'
  YCOMMENT='Start Month for Roads Irrigation'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XRD_START_MONTH(:),IRESP,HCOMMENT=YCOMMENT)
!
  YRECFM='D_RD_EM_IRR'
  YCOMMENT='End   Month for Roads Irrigation'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XRD_END_MONTH  (:),IRESP,HCOMMENT=YCOMMENT)
!
  YRECFM='D_RD_SH_IRR'
  YCOMMENT='Start Hour  for Roads Irrigation'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XRD_START_HOUR (:),IRESP,HCOMMENT=YCOMMENT)
!
  YRECFM='D_RD_EH_IRR'
  YCOMMENT='End   Hour  for Roads Irrigation'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XRD_END_HOUR   (:),IRESP,HCOMMENT=YCOMMENT)
!
  YRECFM='D_RD_IRRIG'
  YCOMMENT='24h mean Irrigation rate for Roads Irrigation'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XRD_24H_IRRIG  (:),IRESP,HCOMMENT=YCOMMENT)
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_PGD_TEB_IRRIG_n',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITESURF_PGD_TEB_IRRIG_n
