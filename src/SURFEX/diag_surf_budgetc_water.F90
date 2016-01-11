!     #########
       SUBROUTINE DIAG_SURF_BUDGETC_WATER (DGWC,  &
                                           PTSTEP, PRN, PH, PLE, PLEI, PGFLUX,  &
                                            PSWD, PSWU, PLWD, PLWU, PFMU, PFMV,&
                                            PEVAP, PSUBL                       )  
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
TYPE(DIAG_t), INTENT(INOUT) :: DGWC
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
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_BUDGETC_WATER',0,ZHOOK_HANDLE)
DGWC%XSWD(:) = DGWC%XSWD(:) + PSWD(:) * PTSTEP
DGWC%XSWU(:) = DGWC%XSWU(:) + PSWU(:) * PTSTEP
!
!*incoming outgoing LW
!
DGWC%XLWD(:) = DGWC%XLWD(:) + PLWD(:) * PTSTEP
DGWC%XLWU(:) = DGWC%XLWU(:) + PLWU(:) * PTSTEP
!
!* net radiation
!
DGWC%XRN(:) = DGWC%XRN(:) + PRN(:) * PTSTEP
!
!* sensible heat flux
!
DGWC%XH(:) = DGWC%XH(:) + PH(:) * PTSTEP 
!
!* latent heat flux (J/m2)
!
DGWC%XLE (:) = DGWC%XLE (:) + PLE (:) * PTSTEP 
DGWC%XLEI(:) = DGWC%XLEI(:) + PLEI(:) * PTSTEP 
!
!* evaporation and sublimation (kg/m2)
!
DGWC%XEVAP(:) = DGWC%XEVAP(:) + PEVAP(:) * PTSTEP
DGWC%XSUBL(:) = DGWC%XSUBL(:) + PSUBL(:) * PTSTEP
!
!* storage flux
!
DGWC%XGFLUX(:) = DGWC%XGFLUX(:) + PGFLUX(:) * PTSTEP 
!
!* wind stress
!
DGWC%XFMU(:) = DGWC%XFMU(:) + PFMU(:) * PTSTEP 
DGWC%XFMV(:) = DGWC%XFMV(:) + PFMV(:) * PTSTEP
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_BUDGETC_WATER',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_SURF_BUDGETC_WATER
