!     #########
       SUBROUTINE DIAG_SURF_BUDGETC_IDEAL(DGL, DGLC, PTSTEP, PRN, PH,          &
                                          PLE, PLEI, PGFLUX, PSWD, PSWU, PLWD, &
                                          PLWU, PFMU, PFMV, PEVAP, PSUBL       )  
!     #########################################################################
!
!!****  *DIAG_SURF_BUDGETC_IDEAL * - Computes cumulated diagnostics in IDEAL case
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
!!     P. Le Moigne 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    03/2015
!!------------------------------------------------------------------
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
TYPE(DIAG_t) :: DGL
TYPE(DIAG_t) :: DGLC
!
REAL,               INTENT(IN) :: PTSTEP    
REAL, DIMENSION(:), INTENT(IN) :: PRN      ! net radiation                         (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PH       ! sensible heat flux                    (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PLE      ! total latent heat flux                (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PLEI     ! sublimation latent heat flux          (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PGFLUX   ! storage flux                          (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PEVAP    ! total evaporation                     (kg/m2/s)
REAL, DIMENSION(:), INTENT(IN) :: PSUBL    ! sublimation                           (kg/m2/s)
REAL, DIMENSION(:), INTENT(IN) :: PSWD     ! total incoming short wave radiation (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PSWU     ! total upward short wave radiation (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PLWD     ! Downward long wave radiation (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PLWU     ! upward long wave radiation (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PFMU     ! zonal wind stress
REAL, DIMENSION(:), INTENT(IN) :: PFMV     ! meridian wind stress
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
!
!* total incoming and outgoing SW
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_BUDGETC_IDEAL',0,ZHOOK_HANDLE)
DGLC%XSWD(:) = DGLC%XSWD(:) + PSWD(:) * PTSTEP
DGLC%XSWU(:) = DGLC%XSWU(:) + PSWU(:) * PTSTEP
!
!*incoming outgoing LW
!
DGLC%XLWD(:) = DGLC%XLWD(:) + PLWD(:) * PTSTEP
DGLC%XLWU(:) = DGLC%XLWU(:) + PLWU(:) * PTSTEP
!
!* net radiation
!
DGLC%XRN(:) = DGLC%XRN(:) + PRN(:) * PTSTEP
!
!* sensible heat flux
!
DGLC%XH(:) = DGLC%XH(:) + PH(:) * PTSTEP 
!
!* latent heat flux
!
DGLC%XLE(:) = DGLC%XLE(:) + PLE (:) * PTSTEP 
DGLC%XLEI(:) = DGLC%XLEI(:) + PLEI(:) * PTSTEP 
!
!* evaporation and sublimation (kg/m2)
!
DGLC%XEVAP(:) = DGLC%XEVAP(:) + PEVAP(:) * PTSTEP
DGLC%XSUBL(:) = DGLC%XSUBL(:) + PSUBL(:) * PTSTEP
!
!* storage flux
!
DGLC%XGFLUX(:) = DGLC%XGFLUX(:) + PGFLUX(:) * PTSTEP 
!
!* wind stress
!
DGLC%XFMU(:) = DGLC%XFMU(:) + PFMU(:) * PTSTEP 
DGLC%XFMV(:) = DGLC%XFMV(:) + PFMV(:) * PTSTEP
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_BUDGETC_IDEAL',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_SURF_BUDGETC_IDEAL
