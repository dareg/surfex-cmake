!     #########
       SUBROUTINE DIAG_SURF_BUDGETC_SEA (DGSC, DGSIC, &
                                          PTSTEP, PRN, PH, PLE, PLE_ICE, PGFLUX,&
                                         PSWD, PSWU, PLWD, PLWU, PFMU, PFMV,   &  
                                         PEVAP, PSUBL, OHANDLE_SIC,            &
                                         PRN_ICE, PH_ICE, PGFLUX_ICE,          &
                                         PSWU_ICE, PLWU_ICE, PFMU_ICE, PFMV_ICE)  
!     ########################################################################
!
!!****  *DIAG_SURF_BUDGETC_SEA * - Computes cumulated diagnostics over sea
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
!!     B. Decharme 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    08/2009
!!      S.Senesi    01/2014  Add fluxes on seaice
!!------------------------------------------------------------------
! 
!
!
!
USE MODD_DIAG_n, ONLY : DIAG_t
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(DIAG_t), INTENT(INOUT) :: DGSC
TYPE(DIAG_t), INTENT(INOUT) :: DGSIC
!
REAL,               INTENT(IN) :: PTSTEP    
REAL, DIMENSION(:), INTENT(IN) :: PRN      ! net radiation                         (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PH       ! sensible heat flux                    (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PLE      ! total latent heat flux                (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PLE_ICE  ! sublimation latent heat flux          (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PGFLUX   ! storage flux                          (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PEVAP    ! total evaporation                     (kg/m2/s)
REAL, DIMENSION(:), INTENT(IN) :: PSUBL    ! sublimation                           (kg/m2/s)
REAL, DIMENSION(:), INTENT(IN) :: PSWD     ! total incoming short wave radiation (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PSWU     ! total upward short wave radiation (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PLWD     ! Downward long wave radiation (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PLWU     ! upward long wave radiation (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PFMU     ! zonal wind stress
REAL, DIMENSION(:), INTENT(IN) :: PFMV     ! meridian wind stress
!
LOGICAL, INTENT(IN)         :: OHANDLE_SIC  ! Do we weight seaice and open sea fluxes
!
REAL, DIMENSION(:), INTENT(IN) :: PRN_ICE  ! net radiation                         (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PH_ICE   ! sensible heat flux                    (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PGFLUX_ICE!storage flux                          (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PSWU_ICE ! total upward short wave radiation (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PLWU_ICE ! upward long wave radiation (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PFMU_ICE ! zonal wind stress
REAL, DIMENSION(:), INTENT(IN) :: PFMV_ICE ! meridian wind stress
!
!
!*      0.2    declarations of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_BUDGETC_SEA',0,ZHOOK_HANDLE)
!
!* total incoming and outgoing SW
!
DGSC%XSWD(:) = DGSC%XSWD(:) + PSWD(:) * PTSTEP
DGSC%XSWU(:) = DGSC%XSWU(:) + PSWU(:) * PTSTEP
!
!*incoming outgoing LW
!
DGSC%XLWD(:) = DGSC%XLWD(:) + PLWD(:) * PTSTEP
DGSC%XLWU(:) = DGSC%XLWU(:) + PLWU(:) * PTSTEP
!
!* net radiation
!
DGSC%XRN(:) = DGSC%XRN(:) + PRN(:) * PTSTEP
!
!* sensible heat flux
!
DGSC%XH(:) = DGSC%XH(:) + PH(:) * PTSTEP 
!
!* latent heat flux (J/m2)
!
DGSC%XLE    (:) = DGSC%XLE    (:) + PLE    (:) * PTSTEP 
DGSC%XLEI   (:) = DGSC%XLEI   (:) + PLE_ICE(:) * PTSTEP 
IF (OHANDLE_SIC) DGSIC%XLE(:) = DGSC%XLEI (:)
!
!* evaporation and sublimation (kg/m2)
!
DGSC%XEVAP(:) = DGSC%XEVAP(:) + PEVAP(:) * PTSTEP
DGSC%XSUBL(:) = DGSC%XSUBL(:) + PSUBL(:) * PTSTEP
!
!* storage flux
!
DGSC%XGFLUX(:) = DGSC%XGFLUX(:) + PGFLUX(:) * PTSTEP
!
!* wind stress
!
DGSC%XFMU(:) = DGSC%XFMU(:) + PFMU(:) * PTSTEP 
DGSC%XFMV(:) = DGSC%XFMV(:) + PFMV(:) * PTSTEP
!
IF (OHANDLE_SIC) THEN
!
!* total incoming and outgoing SW
!
   DGSIC%XSWU(:) = DGSIC%XSWU(:) + PSWU_ICE(:) * PTSTEP
!
!*incoming outgoing LW
!
   DGSIC%XLWU(:) = DGSIC%XLWU(:) + PLWU_ICE(:) * PTSTEP
!
!* net radiation
!
   DGSIC%XRN(:) = DGSIC%XRN(:) + PRN_ICE(:) * PTSTEP
!
!* sensible heat flux
!
   DGSIC%XH(:) = DGSIC%XH(:) + PH_ICE(:) * PTSTEP 
!
!* storage flux
!
   DGSIC%XGFLUX(:) = DGSIC%XGFLUX(:) + PGFLUX_ICE(:) * PTSTEP 
!
!* wind stress
!
   DGSIC%XFMU(:) = DGSIC%XFMU(:) + PFMU_ICE(:) * PTSTEP 
   DGSIC%XFMV(:) = DGSIC%XFMV(:) + PFMV_ICE(:) * PTSTEP
!        
ENDIF
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_BUDGETC_SEA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_SURF_BUDGETC_SEA
