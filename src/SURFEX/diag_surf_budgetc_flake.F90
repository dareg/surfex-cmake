!     #########
       SUBROUTINE DIAG_SURF_BUDGETC_FLAKE (DGF, DGFC, PTSTEP  )  
!     #########################################################################
!
!!****  *DIAG_SURF_BUDGETC_FLAKE * - Computes cumulated diagnostics over water
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
TYPE(DIAG_t), INTENT(INOUT) :: DGF
TYPE(DIAG_t), INTENT(INOUT) :: DGFC
!
REAL,               INTENT(IN) :: PTSTEP   
!
!*      0.2    declarations of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
!* total incoming and outgoing SW
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_BUDGETC_FLAKE',0,ZHOOK_HANDLE)
DGFC%XSWD(:) = DGFC%XSWD(:) + DGF%XSWD(:) * PTSTEP
DGFC%XSWU(:) = DGFC%XSWU(:) + DGF%XSWU(:) * PTSTEP
!
!*incoming outgoing LW
!
DGFC%XLWD(:) = DGFC%XLWD(:) + DGF%XLWD(:) * PTSTEP
DGFC%XLWU(:) = DGFC%XLWU(:) + DGF%XLWU(:) * PTSTEP
!
!* net radiation
!
DGFC%XRN(:) = DGFC%XRN(:) + DGF%XRN(:) * PTSTEP
!
!* sensible heat flux
!
DGFC%XH(:) = DGFC%XH(:) + DGF%XH(:) * PTSTEP 
!
!* latent heat flux
!
DGFC%XLE (:) = DGFC%XLE (:) + DGF%XLE (:) * PTSTEP 
DGFC%XLEI(:) = DGFC%XLEI(:) + DGF%XLEI(:) * PTSTEP 
!
!* evaporation and sublimation (kg/m2)
!
DGFC%XEVAP(:) = DGFC%XEVAP(:) + DGF%XEVAP(:) * PTSTEP
DGFC%XSUBL(:) = DGFC%XSUBL(:) + DGF%XSUBL(:) * PTSTEP
!
!* storage flux
!
DGFC%XGFLUX(:) = DGFC%XGFLUX(:) + DGF%XGFLUX(:) * PTSTEP 
!
!* wind stress
!
DGFC%XFMU(:) = DGFC%XFMU(:) + DGF%XFMU(:) * PTSTEP 
DGFC%XFMV(:) = DGFC%XFMV(:) + DGF%XFMV(:) * PTSTEP
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_BUDGETC_FLAKE',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_SURF_BUDGETC_FLAKE
