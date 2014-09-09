!     #########
       SUBROUTINE DIAG_SURF_BUDGETC_SEA (PTSTEP, PRN, PH, PLE, PLEI, PGFLUX,&
                                         PSWD, PSWU, PLWD, PLWU, PFMU, PFMV,&  
                                         PEVAP, PSUBL,                      &
                                         PRN_ICE, PH_ICE, PGFLUX_ICE,       &
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
USE MODD_DIAG_SEAFLUX_n, ONLY : XRNC, XHC, XLEC, XLEIC, XGFLUXC, XSWDC,  &
                                  XSWUC, XLWDC, XLWUC, XFMUC, XFMVC,     &
                                  XEVAPC, XSUBLC,                        &
                                  XRN_ICEC, XH_ICEC, XGFLUX_ICEC,        &
                                  XSWU_ICEC, XLWU_ICEC, XFMU_ICEC,       &
                                  XFMV_ICEC
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
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
XSWDC(:) = XSWDC(:) + PSWD(:) * PTSTEP
XSWUC(:) = XSWUC(:) + PSWU(:) * PTSTEP
!
XSWU_ICEC(:) = XSWU_ICEC(:) + PSWU_ICE(:) * PTSTEP
!
!*incoming outgoing LW
!
XLWDC(:) = XLWDC(:) + PLWD(:) * PTSTEP
XLWUC(:) = XLWUC(:) + PLWU(:) * PTSTEP
!
XLWU_ICEC(:) = XLWU_ICEC(:) + PLWu_ICE(:) * PTSTEP
!
!* net radiation
!
XRNC(:) = XRNC(:) + PRN(:) * PTSTEP
!
XRN_ICEC(:) = XRN_ICEC(:) + PRN_ICE(:) * PTSTEP
!
!* sensible heat flux
!
XHC(:) = XHC(:) + PH(:) * PTSTEP 
!
XH_ICEC(:) = XH_ICEC(:) + PH_ICE(:) * PTSTEP 
!
!* latent heat flux (J/m2)
!
XLEC (:) = XLEC (:) + PLE (:) * PTSTEP 
XLEIC(:) = XLEIC(:) + PLEI(:) * PTSTEP 
!
!* evaporation and sublimation (kg/m2)
!
XEVAPC(:) = XEVAPC(:) + PEVAP(:) * PTSTEP
XSUBLC(:) = XSUBLC(:) + PSUBL(:) * PTSTEP
!
!* storage flux
!
XGFLUXC(:) = XGFLUXC(:) + PGFLUX(:) * PTSTEP 
!
XGFLUX_ICEC(:) = XGFLUX_ICEC(:) + PGFLUX_ICE(:) * PTSTEP 
!
!* wind stress
!
XFMUC(:) = XFMUC(:) + PFMU(:) * PTSTEP 
XFMVC(:) = XFMVC(:) + PFMV(:) * PTSTEP
!
XFMU_ICEC(:) = XFMU_ICEC(:) + PFMU_ICE(:) * PTSTEP 
XFMV_ICEC(:) = XFMV_ICEC(:) + PFMV_ICE(:) * PTSTEP
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_BUDGETC_SEA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_SURF_BUDGETC_SEA
