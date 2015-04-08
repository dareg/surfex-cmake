!     #########
       SUBROUTINE DIAG_INLINE_WATFLUX_n (PTSTEP, PTA, PTS, PQA, PPA, PPS, PRHOA, PZONA,  &
                                           PMERA, PHT, PHW, PCD, PCDN, PCH, PRI, PHU, PZ0, &
                                           PZ0H, PQSAT, PSFTH, PSFTQ, PSFZON, PSFMER,      &
                                           PDIR_SW, PSCA_SW, PLW, PDIR_ALB, PSCA_ALB,      &
                                           PEMIS, PTRAD, PRAIN, PSNOW, PTICE, PSFTH_ICE,   &
                                           PSFTQ_ICE                                       )  
!     ###############################################################################
!
!!****  *DIAG_INLINE_WATFLUX_n * - computes diagnostics during WATFLUX time-step
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
!!      B. Decharme 08/2009 : Diag for Earth System Model Coupling
!!      S. Riette   06/2009 CLS_2M becomes CLS_TQ, CLS_TQ and CLS_WIND have one
!!                          more argument (height of diagnostic)
!       B. decharme 04/2013 : Add EVAP and SUBL diag
!!------------------------------------------------------------------
!

!
!
USE MODD_CSTS,           ONLY : XTT
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_SFX_OASIS,      ONLY : LCPL_SEA
USE MODD_WATFLUX_n, ONLY : W => WATFLUX
USE MODD_DIAG_WATFLUX_n, ONLY : DGW => DIAG_WATFLUX
!
USE MODI_PARAM_CLS
USE MODI_CLS_TQ
USE MODI_CLS_WIND
USE MODI_DIAG_SURF_BUDGET_WATER
USE MODI_DIAG_SURF_BUDGETC_WATER
USE MODI_DIAG_CPL_ESM_WATER
! 
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,               INTENT(IN) :: PTSTEP ! atmospheric time-step                 (s)
REAL, DIMENSION(:), INTENT(IN) :: PTA    ! atmospheric temperature
REAL, DIMENSION(:), INTENT(IN) :: PTS    ! surface temperature
REAL, DIMENSION(:), INTENT(IN) :: PQA    ! atmospheric specific humidity
REAL, DIMENSION(:), INTENT(IN) :: PPA    ! atmospheric level pressure
REAL, DIMENSION(:), INTENT(IN) :: PPS    ! surface pressure
REAL, DIMENSION(:), INTENT(IN) :: PRHOA  ! air density
REAL, DIMENSION(:), INTENT(IN) :: PZONA  ! zonal wind
REAL, DIMENSION(:), INTENT(IN) :: PMERA  ! meridian wind
REAL, DIMENSION(:), INTENT(IN) :: PHT    ! atmospheric level height
REAL, DIMENSION(:), INTENT(IN) :: PHW    ! atmospheric level height for wind
REAL, DIMENSION(:), INTENT(IN) :: PCD    ! drag coefficient for momentum
REAL, DIMENSION(:), INTENT(IN) :: PCDN   ! neutral drag coefficient
REAL, DIMENSION(:), INTENT(IN) :: PCH    ! drag coefficient for heat
REAL, DIMENSION(:), INTENT(IN) :: PRI    ! Richardson number
REAL, DIMENSION(:), INTENT(IN) :: PHU    ! near-surface humidity
REAL, DIMENSION(:), INTENT(IN) :: PZ0    ! roughness length for momentum
REAL, DIMENSION(:), INTENT(IN) :: PZ0H   ! roughness length for heat
REAL, DIMENSION(:), INTENT(IN) :: PQSAT  ! humidity at saturation
REAL, DIMENSION(:), INTENT(IN) :: PSFZON ! zonal friction
REAL, DIMENSION(:), INTENT(IN) :: PSFMER ! meridian friction
REAL, DIMENSION(:), INTENT(IN) :: PSFTH  ! heat flux  (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PSFTQ  ! water flux (kg/m2/s)
REAL, DIMENSION(:,:),INTENT(IN):: PDIR_SW   ! direct  solar radiation (on horizontal surf.)
!                                           !                                       (W/m2)
REAL, DIMENSION(:,:),INTENT(IN):: PSCA_SW   ! diffuse solar radiation (on horizontal surf.)
!                                           !                                       (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PLW       ! longwave radiation (on horizontal surf.)
REAL, DIMENSION(:), INTENT(IN) :: PTRAD     ! radiative temperature                 (K)
REAL, DIMENSION(:,:),INTENT(IN):: PDIR_ALB  ! direct albedo for each spectral band  (-)
REAL, DIMENSION(:,:),INTENT(IN):: PSCA_ALB  ! diffuse albedo for each spectral band (-)
REAL, DIMENSION(:), INTENT(IN) :: PEMIS     ! emissivity                            (-)
!
REAL, DIMENSION(:), INTENT(IN) :: PRAIN     ! Rainfall (kg/m2/s)
REAL, DIMENSION(:), INTENT(IN) :: PSNOW     ! Snowfall (kg/m2/s)
REAL, DIMENSION(:), INTENT(IN) :: PSFTH_ICE ! heat flux  (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PSFTQ_ICE ! water flux (kg/m2/s)
REAL, DIMENSION(:), INTENT(IN) :: PTICE     ! Ice Surface Temperature
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PTA)) :: ZH
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DIAG_INLINE_WATFLUX_N',0,ZHOOK_HANDLE)
!
! * Mean surface temperature need to couple with AGCM
!
DGW%XDIAG_TS(:) = PTS(:)
!
IF (.NOT. W%LSBL) THEN
!
  IF (DGW%N2M==1) THEN
    CALL PARAM_CLS(PTA, PTS, PQA, PPA, PRHOA, PZONA, PMERA, PHT, PHW, &
                     PSFTH, PSFTQ, PSFZON, PSFMER,                    &
                     DGW%XT2M, DGW%XQ2M, DGW%XHU2M, DGW%XZON10M, DGW%XMER10M              )  
  ELSE IF (DGW%N2M==2) THEN
    ZH(:)=2.          
    CALL CLS_TQ(PTA, PQA, PPA, PPS, PHT,         &
                  PCD, PCH, PRI,                   &
                  PTS, PHU, PZ0H, ZH,              &
                  DGW%XT2M, DGW%XQ2M, DGW%XHU2M                )  
    ZH(:)=10.                
    CALL CLS_WIND(PZONA, PMERA, PHW,             &
                    PCD, PCDN, PRI, ZH,            &
                    DGW%XZON10M, DGW%XMER10M               )  
  END IF
!
  IF (DGW%N2M>=1) THEN
    !
    DGW%XT2M_MIN(:) = MIN(DGW%XT2M_MIN(:),DGW%XT2M(:))
    DGW%XT2M_MAX(:) = MAX(DGW%XT2M_MAX(:),DGW%XT2M(:))
    !
    DGW%XHU2M_MIN(:) = MIN(DGW%XHU2M_MIN(:),DGW%XHU2M(:))
    DGW%XHU2M_MAX(:) = MAX(DGW%XHU2M_MAX(:),DGW%XHU2M(:))
    !
    DGW%XWIND10M  (:) = SQRT(DGW%XZON10M(:)**2+DGW%XMER10M(:)**2)
    DGW%XWIND10M_MAX(:) = MAX(DGW%XWIND10M_MAX(:),DGW%XWIND10M(:))
    !
    !* Richardson number
    DGW%XRI = PRI
    !
  ENDIF
!
ELSE
  IF (DGW%N2M>=1) THEN
    DGW%XT2M    = XUNDEF
    DGW%XQ2M    = XUNDEF
    DGW%XHU2M   = XUNDEF
    DGW%XZON10M = XUNDEF
    DGW%XMER10M = XUNDEF
    DGW%XRI     = PRI
  ENDIF
ENDIF
!
IF (DGW%LSURF_BUDGET.OR.DGW%LSURF_BUDGETC) THEN
  !
  CALL  DIAG_SURF_BUDGET_WATER (XTT, PTS, PRHOA, PSFTH, PSFTQ,          &
                                  PDIR_SW, PSCA_SW, PLW,                &
                                  PDIR_ALB, PSCA_ALB, PEMIS, PTRAD,     &
                                  PSFZON, PSFMER,                       &
                                  DGW%XRN, DGW%XH, DGW%XLE, DGW%XLEI, DGW%XGFLUX,           &
                                  DGW%XSWD, DGW%XSWU, DGW%XSWBD, DGW%XSWBU, DGW%XLWD, DGW%XLWU, &
                                  DGW%XFMU, DGW%XFMV, DGW%XEVAP, DGW%XSUBL )  
  !
END IF
!
IF(DGW%LSURF_BUDGETC)THEN
  CALL DIAG_SURF_BUDGETC_WATER(PTSTEP, DGW%XRN, DGW%XH, DGW%XLE, DGW%XLEI, DGW%XGFLUX,  &
                                 DGW%XSWD, DGW%XSWU, DGW%XLWD, DGW%XLWU, DGW%XFMU, DGW%XFMV,&
                                 DGW%XEVAP, DGW%XSUBL                       )  
ENDIF
!
IF (DGW%LCOEF) THEN
  !
  !* Transfer coefficients
  !
  DGW%XCD = PCD
  DGW%XCH = PCH
  DGW%XCE = PCH
  !
  !* Roughness lengths
  !
  DGW%XZ0  = PZ0
  DGW%XZ0H = PZ0H
  !
ENDIF
!
IF (DGW%LSURF_VARS) THEN
  !
  !* Humidity at saturation
  !
  DGW%XQS = PQSAT
  !
ENDIF
!
! Diag for Earth System Model coupling
!
IF (LCPL_SEA) THEN
!
  CALL DIAG_CPL_ESM_WATER(PTSTEP,DGW%XZON10M,DGW%XMER10M,DGW%XFMU,DGW%XFMV,  &
                            DGW%XSWD,DGW%XSWU,DGW%XGFLUX,PSFTQ,PRAIN,      &
                            PSNOW,PLW,PTICE,PSFTH_ICE,         &
                            PSFTQ_ICE,PDIR_SW,PSCA_SW          )  
! 
ENDIF
IF (LHOOK) CALL DR_HOOK('DIAG_INLINE_WATFLUX_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_INLINE_WATFLUX_n
