!     #########
       SUBROUTINE DIAG_CPL_ESM_FLAKE (PTSTEP,PRAIN,PSNOW,PSFTQ)  
!     #####################################################################
!
!!****  *DIAG_CPL_ESM_FLAKE * - Computes diagnostics over sea for 
!!                                Earth system model coupling
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
USE MODD_FLAKE_n, ONLY : F => FLAKE
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,               INTENT(IN) :: PTSTEP    ! atmospheric time-step
REAL, DIMENSION(:), INTENT(IN) :: PRAIN     ! Rainfall
REAL, DIMENSION(:), INTENT(IN) :: PSNOW     ! Snowfall
REAL, DIMENSION(:), INTENT(IN) :: PSFTQ     ! water flux
!
!*      0.2    declarations of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('DIAG_CPL_ESM_FLAKE',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
! Total flux
!-------------------------------------------------------------------------------------
!
!* Evaporation (kg/m2)
!
F%XCPL_FLAKE_EVAP(:) = F%XCPL_FLAKE_EVAP(:) + PTSTEP * PSFTQ(:)
!
!* Precip (kg/m2)
! 
F%XCPL_FLAKE_RAIN(:) = F%XCPL_FLAKE_RAIN(:) + PTSTEP * PRAIN(:) 
F%XCPL_FLAKE_SNOW(:) = F%XCPL_FLAKE_SNOW(:) + PTSTEP * PSNOW(:)
!
IF (LHOOK) CALL DR_HOOK('DIAG_CPL_ESM_FLAKE',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_CPL_ESM_FLAKE
