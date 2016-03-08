!     #########
       SUBROUTINE DIAG_SURF_BUDGETC_SEA (DGS, DGSI, DGSC, DGSIC, PTSTEP, OHANDLE_SIC)  
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
TYPE(DIAG_t), INTENT(INOUT) :: DGS
TYPE(DIAG_t), INTENT(INOUT) :: DGSI
TYPE(DIAG_t), INTENT(INOUT) :: DGSC
TYPE(DIAG_t), INTENT(INOUT) :: DGSIC
!
REAL,               INTENT(IN) :: PTSTEP    
!
LOGICAL, INTENT(IN)         :: OHANDLE_SIC  ! Do we weight seaice and open sea fluxes
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
DGSC%XSWD(:) = DGSC%XSWD(:) + DGS%XSWD(:) * PTSTEP
DGSC%XSWU(:) = DGSC%XSWU(:) + DGS%XSWU(:) * PTSTEP
!
!*incoming outgoing LW
!
DGSC%XLWD(:) = DGSC%XLWD(:) + DGS%XLWD(:) * PTSTEP
DGSC%XLWU(:) = DGSC%XLWU(:) + DGS%XLWU(:) * PTSTEP
!
!* net radiation
!
DGSC%XRN(:) = DGSC%XRN(:) + DGS%XRN(:) * PTSTEP
!
!* sensible heat flux
!
DGSC%XH(:) = DGSC%XH(:) + DGS%XH(:) * PTSTEP 
!
!* latent heat flux (J/m2)
!
DGSC%XLE    (:) = DGSC%XLE    (:) + DGS%XLE    (:) * PTSTEP 
DGSC%XLEI   (:) = DGSC%XLEI   (:) + DGS%XLEI(:) * PTSTEP 
IF (OHANDLE_SIC) DGSIC%XLE(:) = DGSC%XLEI (:)
!
!* evaporation and sublimation (kg/m2)
!
DGSC%XEVAP(:) = DGSC%XEVAP(:) + DGS%XEVAP(:) * PTSTEP
DGSC%XSUBL(:) = DGSC%XSUBL(:) + DGS%XSUBL(:) * PTSTEP
!
!* storage flux
!
DGSC%XGFLUX(:) = DGSC%XGFLUX(:) + DGS%XGFLUX(:) * PTSTEP
!
!* wind stress
!
DGSC%XFMU(:) = DGSC%XFMU(:) + DGS%XFMU(:) * PTSTEP 
DGSC%XFMV(:) = DGSC%XFMV(:) + DGS%XFMV(:) * PTSTEP
!
IF (OHANDLE_SIC) THEN
!
!* total incoming and outgoing SW
!
   DGSIC%XSWU(:) = DGSIC%XSWU(:) + DGSI%XSWU(:) * PTSTEP
!
!*incoming outgoing LW
!
   DGSIC%XLWU(:) = DGSIC%XLWU(:) + DGSI%XLWU(:) * PTSTEP
!
!* net radiation
!
   DGSIC%XRN(:) = DGSIC%XRN(:) + DGSI%XRN(:) * PTSTEP
!
!* sensible heat flux
!
   DGSIC%XH(:) = DGSIC%XH(:) + DGSI%XH(:) * PTSTEP 
!
!* storage flux
!
   DGSIC%XGFLUX(:) = DGSIC%XGFLUX(:) + DGSI%XGFLUX(:) * PTSTEP 
!
!* wind stress
!
   DGSIC%XFMU(:) = DGSIC%XFMU(:) + DGSI%XFMU(:) * PTSTEP 
   DGSIC%XFMV(:) = DGSIC%XFMV(:) + DGSI%XFMV(:) * PTSTEP
!        
ENDIF
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_BUDGETC_SEA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_SURF_BUDGETC_SEA
