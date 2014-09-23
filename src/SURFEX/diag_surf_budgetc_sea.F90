!     #########
       SUBROUTINE DIAG_SURF_BUDGETC_SEA (PTSTEP, PRN, PH, PLE, PLE_ICE, PGFLUX,&
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
USE MODD_DIAG_SEAFLUX_n, ONLY : XRNC, XHC, XLEC, XLEC_ICE, XGFLUXC, XSWDC,  &
                                  XSWUC, XLWDC, XLWUC, XFMUC, XFMVC,     &
                                  XEVAPC, XSUBLC,                        &
                                  XRNC_ICE, XHC_ICE, XGFLUXC_ICE,        &
                                  XSWUC_ICE, XLWUC_ICE, XFMUC_ICE,       &
                                  XFMVC_ICE
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
XSWDC(:) = XSWDC(:) + PSWD(:) * PTSTEP
XSWUC(:) = XSWUC(:) + PSWU(:) * PTSTEP
!
!*incoming outgoing LW
!
XLWDC(:) = XLWDC(:) + PLWD(:) * PTSTEP
XLWUC(:) = XLWUC(:) + PLWU(:) * PTSTEP
!
!* net radiation
!
XRNC(:) = XRNC(:) + PRN(:) * PTSTEP
!
!* sensible heat flux
!
XHC(:) = XHC(:) + PH(:) * PTSTEP 
!
!* latent heat flux (J/m2)
!
XLEC    (:) = XLEC    (:) + PLE    (:) * PTSTEP 
XLEC_ICE(:) = XLEC_ICE(:) + PLE_ICE(:) * PTSTEP 
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
!* wind stress
!
XFMUC(:) = XFMUC(:) + PFMU(:) * PTSTEP 
XFMVC(:) = XFMVC(:) + PFMV(:) * PTSTEP
!
IF (OHANDLE_SIC) THEN
!
!* total incoming and outgoing SW
!
   XSWUC_ICE(:) = XSWUC_ICE(:) + PSWU_ICE(:) * PTSTEP
!
!*incoming outgoing LW
!
   XLWUC_ICE(:) = XLWUC_ICE(:) + PLWU_ICE(:) * PTSTEP
!
!* net radiation
!
   XRNC_ICE(:) = XRNC_ICE(:) + PRN_ICE(:) * PTSTEP
!
!* sensible heat flux
!
   XHC_ICE(:) = XHC_ICE(:) + PH_ICE(:) * PTSTEP 
!
!* storage flux
!
   XGFLUXC(:) = XGFLUXC(:) + PGFLUX(:) * PTSTEP 
!
!* wind stress
!
   XFMUC(:) = XFMUC(:) + PFMU(:) * PTSTEP 
   XFMVC(:) = XFMVC(:) + PFMV(:) * PTSTEP
!        
ENDIF
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_BUDGETC_SEA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_SURF_BUDGETC_SEA
