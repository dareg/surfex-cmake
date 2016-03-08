!     #########
       SUBROUTINE DIAG_SURF_BUDGETC_WATER (DGW, DGWC, PTSTEP )  
!     #########################################################################
!
!!****  *DIAG_SURF_BUDGETC_WATER * - Computes cumulated diagnostics over water
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
TYPE(DIAG_t), INTENT(INOUT) :: DGW
TYPE(DIAG_t), INTENT(INOUT) :: DGWC
!
REAL,               INTENT(IN) :: PTSTEP    
!
!*      0.2    declarations of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
!* total incoming and outgoing SW
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_BUDGETC_WATER',0,ZHOOK_HANDLE)
DGWC%XSWD(:) = DGWC%XSWD(:) + DGW%XSWD(:) * PTSTEP
DGWC%XSWU(:) = DGWC%XSWU(:) + DGW%XSWU(:) * PTSTEP
!
!*incoming outgoing LW
!
DGWC%XLWD(:) = DGWC%XLWD(:) + DGW%XLWD(:) * PTSTEP
DGWC%XLWU(:) = DGWC%XLWU(:) + DGW%XLWU(:) * PTSTEP
!
!* net radiation
!
DGWC%XRN(:) = DGWC%XRN(:) + DGW%XRN(:) * PTSTEP
!
!* sensible heat flux
!
DGWC%XH(:) = DGWC%XH(:) + DGW%XH(:) * PTSTEP 
!
!* latent heat flux (J/m2)
!
DGWC%XLE (:) = DGWC%XLE (:) + DGW%XLE (:) * PTSTEP 
DGWC%XLEI(:) = DGWC%XLEI(:) + DGW%XLEI(:) * PTSTEP 
!
!* evaporation and sublimation (kg/m2)
!
DGWC%XEVAP(:) = DGWC%XEVAP(:) + DGW%XEVAP(:) * PTSTEP
DGWC%XSUBL(:) = DGWC%XSUBL(:) + DGW%XSUBL(:) * PTSTEP
!
!* storage flux
!
DGWC%XGFLUX(:) = DGWC%XGFLUX(:) + DGW%XGFLUX(:) * PTSTEP 
!
!* wind stress
!
DGWC%XFMU(:) = DGWC%XFMU(:) + DGW%XFMU(:) * PTSTEP 
DGWC%XFMV(:) = DGWC%XFMV(:) + DGW%XFMV(:) * PTSTEP
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_BUDGETC_WATER',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_SURF_BUDGETC_WATER
