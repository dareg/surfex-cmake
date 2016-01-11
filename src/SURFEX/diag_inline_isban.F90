!     #########
 SUBROUTINE DIAG_INLINE_ISBA_n (DGI, PKD, OSURF_BUDGETC, OCANOPY, &
                                 PTA, PTS, PQA, PPA, PPS, PRHOA, PZONA, PMERA,  &
                                  PHT, PHW, PCD, PCDN, PCH, PRI, PHU, PZ0, PZ0H, &
                                  PZ0EFF, PSFTH, PSFTQ, PSFZON, PSFMER, PQS,     &
                                  PDIR_ALB, PSCA_ALB, PDIR_SW, PSCA_SW, PLW, PRN )  
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
USE MODD_DIAG_n, ONLY : DIAG_t
USE MODD_PACK_DIAG_ISBA, ONLY : PACK_DIAG_ISBA_t
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
TYPE(DIAG_t), INTENT(INOUT) :: DGI
TYPE(PACK_DIAG_ISBA_t), INTENT(INOUT) :: PKD
!
LOGICAL, INTENT(IN) :: OSURF_BUDGETC
LOGICAL, INTENT(IN) :: OCANOPY
!
REAL, DIMENSION(:), INTENT(IN)       :: PTA      ! atmospheric temperature
REAL, DIMENSION(:), INTENT(IN)       :: PTS      ! surface temperature
REAL, DIMENSION(:), INTENT(IN)       :: PQA      ! atmospheric specific humidity
REAL, DIMENSION(:), INTENT(IN)       :: PPA      ! atmospheric level pressure
REAL, DIMENSION(:), INTENT(IN)       :: PPS      ! surface pressure
REAL, DIMENSION(:), INTENT(IN)       :: PRHOA    ! air density
REAL, DIMENSION(:), INTENT(IN)       :: PZONA    ! zonal wind
REAL, DIMENSION(:), INTENT(IN)       :: PMERA    ! meridian wind
REAL, DIMENSION(:), INTENT(IN)       :: PHT      ! atmospheric level height
REAL, DIMENSION(:), INTENT(IN)       :: PHW      ! atmospheric level height for wind
REAL, DIMENSION(:), INTENT(IN)       :: PCD      ! drag coefficient for momentum
REAL, DIMENSION(:), INTENT(IN)       :: PCDN     ! neutral drag coefficient
REAL, DIMENSION(:), INTENT(IN)       :: PCH      ! drag coefficient for heat
REAL, DIMENSION(:), INTENT(IN)       :: PRI      ! Richardson number
REAL, DIMENSION(:), INTENT(IN)       :: PHU      ! near-surface humidity
REAL, DIMENSION(:), INTENT(IN)       :: PZ0      ! roughness length for momentum
REAL, DIMENSION(:), INTENT(IN)       :: PZ0H     ! roughness length for heat
REAL, DIMENSION(:), INTENT(IN)       :: PZ0EFF   ! effective roughness length (z0+z0rel)
REAL, DIMENSION(:), INTENT(IN)       :: PQS      ! humidity at surface 
REAL, DIMENSION(:,:), INTENT(IN)     :: PDIR_ALB ! direct albedo for each spectral band
REAL, DIMENSION(:,:), INTENT(IN)     :: PSCA_ALB ! diffuse albedo for each spectral band (-)
REAL, DIMENSION(:,:), INTENT(IN)     :: PDIR_SW  ! direct  solar radiation (on horizontal surf.)
REAL, DIMENSION(:,:), INTENT(IN)     :: PSCA_SW  ! diffuse solar radiation (on horizontal surf.)
REAL, DIMENSION(:), INTENT(IN)       :: PLW      ! longwave radiation (on horizontal surf.)
REAL, DIMENSION(:), INTENT(IN)       :: PRN      ! Surface net radiation
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
! * Richardson number
!
IF (DGI%N2M>=1) THEN
    PKD%DGIP%XRI     = PRI        
ENDIF
!
! * Near surface atmospheric variables
!
IF (.NOT. OCANOPY) THEN
!        
  IF (DGI%N2M==1) THEN
    CALL PARAM_CLS(PTA, PTS, PQA, PPA, PRHOA, PZONA, PMERA, PHT, PHW,  &
                     PSFTH, PSFTQ, PSFZON, PSFMER,                     &
                     PKD%DGIP%XT2M, PKD%DGIP%XQ2M, PKD%DGIP%XHU2M, PKD%DGIP%XZON10M, PKD%DGIP%XMER10M     )  
  ELSE IF (DGI%N2M==2) THEN
    ZH(:)=2.          
    CALL CLS_TQ(PTA, PQA, PPA, PPS, PHT,           &
                  PCD, PCH, PRI,                   &
                  PTS, PHU, PZ0H, ZH,              &
                  PKD%DGIP%XT2M, PKD%DGIP%XQ2M, PKD%DGIP%XHU2M          )  
    ZH(:)=10.                
    CALL CLS_WIND(PZONA, PMERA, PHW,               &
                    PCD, PCDN, PRI, ZH,            &
                    PKD%DGIP%XZON10M, PKD%DGIP%XMER10M           )  
  END IF
!
ELSE
  !        
  IF (DGI%N2M>=1) THEN
    PKD%DGIP%XT2M    = XUNDEF
    PKD%DGIP%XQ2M    = XUNDEF
    PKD%DGIP%XHU2M   = XUNDEF
    PKD%DGIP%XZON10M = XUNDEF
    PKD%DGIP%XMER10M = XUNDEF
  ENDIF
  !        
ENDIF
!
! * Surface energy budget
!
IF (DGI%LSURF_BUDGET.OR.OSURF_BUDGETC) THEN
   !
   CALL DIAG_SURF_BUDGET_ISBA(PDIR_SW, PSCA_SW, PDIR_ALB, PSCA_ALB,  &
                                PLW, PRN,                              &
                                PKD%DGIP%XSWD, PKD%DGIP%XSWU, PKD%DGIP%XSWBD, PKD%DGIP%XSWBU,      &
                                PKD%DGIP%XLWD, PKD%DGIP%XLWU   )          
   !
   PKD%DGIP%XFMU = PSFZON
   PKD%DGIP%XFMV = PSFMER
   !
END IF
!
IF (DGI%LCOEF) THEN
  !
  !* Transfer coefficient
  !
  PKD%DGIP%XCD = PCD
  PKD%DGIP%XCH = PCH
  PKD%DGIP%XCE = PCH
  !
  !* Roughness lengths
  !
  PKD%DGIP%XZ0  = PZ0
  PKD%DGIP%XZ0H = PZ0H
  PKD%DGIP%XZ0EFF         = PZ0EFF
  !
ENDIF
!
IF (DGI%LSURF_VARS) THEN
  !
  !* Humidity at surface
  !
  PKD%DGIP%XQS = PQS
  !
ENDIF
IF (LHOOK) CALL DR_HOOK('DIAG_INLINE_ISBA_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_INLINE_ISBA_n
