!     #########
 SUBROUTINE DIAG_INLINE_ISBA_n (DIO, INI, DGIP, OCANOPY, PTA, PQA, PPA, PPS, PRHOA, PZONA, PMERA, &
                                  PHT, PHW, PSFTH, PSFTQ, PSFZON, PSFMER, PDIR_SW, PSCA_SW, PLW )  
!     ###############################################################################
!
!!****  *DIAG_INLINE_ISBA_n * - computes diagnostics during ISBA time-step
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     V. Masson 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      B. Decharme 08/2009 caculate cumulated diag LSURF_BUDGETC
!!      S. Riette   06/2009 CLS_2M becomes CLS_TQ, CLS_TQ and CLS_WIND have one
!!                          more argument (height of diagnostic)
!!------------------------------------------------------------------
!
!
!
!
USE MODD_DIAG_n, ONLY : DIAG_t, DIAG_OPTIONS_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_t
!
USE MODD_SURF_PAR,         ONLY : XUNDEF
!
USE MODI_PARAM_CLS
USE MODI_CLS_TQ
USE MODI_CLS_WIND
USE MODI_DIAG_SURF_BUDGET_ISBA
! 
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(DIAG_OPTIONS_t), INTENT(IN) :: DIO
TYPE(ISBA_INIT_t), INTENT(INOUT) :: INI
TYPE(DIAG_t), INTENT(INOUT) :: DGIP
!
LOGICAL, INTENT(IN) :: OCANOPY
!
REAL, DIMENSION(:), INTENT(IN)       :: PTA      ! atmospheric temperature
REAL, DIMENSION(:), INTENT(IN)       :: PQA      ! atmospheric specific humidity
REAL, DIMENSION(:), INTENT(IN)       :: PPA      ! atmospheric level pressure
REAL, DIMENSION(:), INTENT(IN)       :: PPS      ! surface pressure
REAL, DIMENSION(:), INTENT(IN)       :: PRHOA    ! air density
REAL, DIMENSION(:), INTENT(IN)       :: PZONA    ! zonal wind
REAL, DIMENSION(:), INTENT(IN)       :: PMERA    ! meridian wind
REAL, DIMENSION(:), INTENT(IN)       :: PHT      ! atmospheric level height
REAL, DIMENSION(:), INTENT(IN)       :: PHW      ! atmospheric level height for wind
REAL, DIMENSION(:,:), INTENT(IN)     :: PDIR_SW  ! direct  solar radiation (on horizontal surf.)
REAL, DIMENSION(:,:), INTENT(IN)     :: PSCA_SW  ! diffuse solar radiation (on horizontal surf.)
REAL, DIMENSION(:), INTENT(IN)       :: PLW      ! longwave radiation (on horizontal surf.)
!
REAL, DIMENSION(:), INTENT(IN)       :: PSFZON   ! zonal friction
REAL, DIMENSION(:), INTENT(IN)       :: PSFMER   ! meridian friction
REAL, DIMENSION(:), INTENT(IN)       :: PSFTH    ! heat flux (W/m2)
REAL, DIMENSION(:), INTENT(IN)       :: PSFTQ    ! water flux (kg/m2/s)
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PTA)) :: ZH
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DIAG_INLINE_ISBA_N',0,ZHOOK_HANDLE)
!
! * Near surface atmospheric variables
!
IF (.NOT. OCANOPY) THEN
!        
  IF (DIO%N2M==1) THEN
    CALL PARAM_CLS(DGIP, PTA, DGIP%XTS, PQA, PPA, PRHOA, PZONA, PMERA, &
                   PHT, PHW, PSFTH, PSFTQ, PSFZON, PSFMER    )  
  ELSE IF (DIO%N2M==2) THEN
    ZH(:)=2.          
    CALL CLS_TQ(PTA, PQA, PPA, PPS, PHT, DGIP%XCD, DGIP%XCH, DGIP%XRI, &
                DGIP%XTS, DGIP%XHU, DGIP%XZ0H, ZH, DGIP%XT2M, DGIP%XQ2M, DGIP%XHU2M )  
    ZH(:)=10.                
    CALL CLS_WIND(PZONA, PMERA, PHW, DGIP%XCD, DGIP%XCDN, DGIP%XRI, ZH, &
                 DGIP%XZON10M, DGIP%XMER10M  )  
  END IF
!
ELSE
  !        
  IF (DIO%N2M>=1) THEN
    DGIP%XT2M    = XUNDEF
    DGIP%XQ2M    = XUNDEF
    DGIP%XHU2M   = XUNDEF
    DGIP%XZON10M = XUNDEF
    DGIP%XMER10M = XUNDEF
  ENDIF
  !        
ENDIF
!
! * Surface energy budget
!
IF (DIO%LSURF_BUDGET.OR.DIO%LSURF_BUDGETC) THEN
   !
   CALL DIAG_SURF_BUDGET_ISBA(PDIR_SW, PSCA_SW, PLW, INI, DGIP)          
   !
   DGIP%XFMU = PSFZON
   DGIP%XFMV = PSFMER
   !
END IF
!
IF (DIO%LCOEF) THEN
  !
  !* Transfer coefficient
  !
  DGIP%XCD = DGIP%XCD
  DGIP%XCH = DGIP%XCH
  DGIP%XCE = DGIP%XCH
  !
  !* Roughness lengths
  !
  DGIP%XZ0    = DGIP%XZ0
  DGIP%XZ0H   = DGIP%XZ0H
  DGIP%XZ0EFF = DGIP%XZ0EFF
  !
ENDIF
!
IF (DIO%LSURF_VARS) THEN
  !
  !* Humidity at surface
  !
  DGIP%XQS = DGIP%XQS
  !
ENDIF
IF (LHOOK) CALL DR_HOOK('DIAG_INLINE_ISBA_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_INLINE_ISBA_n
