!     #########
SUBROUTINE DIAG_ISBA_n (DGEI, DGEIC, DGI, DGIC, &
                        HPROGRAM,                                               &
                         PRN, PH, PLE, PLEI, PGFLUX, PRI, PCD, PCH, PCE, PQS,    &
                         PZ0, PZ0H, PT2M, PTS, PQ2M, PHU2M, PZON10M, PMER10M,    &
                         PSWD, PSWU, PLWD, PLWU, PSWBD, PSWBU, PFMU, PFMV,       &
                         PRNC, PHC, PLEC, PGFLUXC, PSWDC, PSWUC, PLWDC,          &
                         PLWUC, PFMUC, PFMVC, PT2M_MIN, PT2M_MAX, PLEIC,         &
                         PHU2M_MIN, PHU2M_MAX, PWIND10M, PWIND10M_MAX,           &
                            PEVAP, PEVAPC, PSUBL, PSUBLC                         )
!     ###############################################################################
!
!!****  *DIAG_ISBA_n * - Stores ISBA diagnostics
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
!!      Modified    01/2006 : sea flux parameterization.
!!      Modified    08/2009 : new diag
!       B. decharme 04/2013 : Add EVAP and SUBL diag
!!------------------------------------------------------------------
!
!
!
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_n, ONLY : DIAG_t
!
USE MODD_SURF_PAR,    ONLY : XUNDEF
!                                  
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEI
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEIC
TYPE(DIAG_t), INTENT(INOUT) :: DGI
TYPE(DIAG_t), INTENT(INOUT) :: DGIC
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM ! program calling surf. schemes
!
REAL, DIMENSION(:), INTENT(OUT) :: PRN      ! Net radiation       (W/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PH       ! Sensible heat flux  (W/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PLE      ! Total latent heat flux    (W/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PLEI     ! Sublimation latent heat flux    (W/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PGFLUX   ! Storage flux        (W/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PEVAP    ! Total evapotranspiration  (kg/m2/s)
REAL, DIMENSION(:), INTENT(OUT) :: PSUBL    ! Sublimation (kg/m2/s)
REAL, DIMENSION(:), INTENT(OUT) :: PRI      ! Richardson number   (-)
REAL, DIMENSION(:), INTENT(OUT) :: PCD      ! drag coefficient    (W/s2)
REAL, DIMENSION(:), INTENT(OUT) :: PCH      ! transf. coef heat   (W/s)
REAL, DIMENSION(:), INTENT(OUT) :: PCE      ! transf. coef vapor  (W/s/K)
REAL, DIMENSION(:), INTENT(OUT) :: PQS
REAL, DIMENSION(:), INTENT(OUT) :: PZ0      ! rough. length wind  (m)
REAL, DIMENSION(:), INTENT(OUT) :: PZ0H     ! rough. length heat  (m)
REAL, DIMENSION(:), INTENT(OUT) :: PTS      ! surface temperature (K)
REAL, DIMENSION(:), INTENT(OUT) :: PT2M     ! temperature at 2m   (K)
REAL, DIMENSION(:), INTENT(OUT) :: PQ2M     ! humidity at 2m      (kg/kg)
REAL, DIMENSION(:), INTENT(OUT) :: PHU2M    ! relative humidity at 2m (-)
REAL, DIMENSION(:), INTENT(OUT) :: PZON10M  ! zonal wind at 10m   (m/s)
REAL, DIMENSION(:), INTENT(OUT) :: PMER10M  ! meridian wind at 10m(m/s)
REAL, DIMENSION(:), INTENT(OUT) :: PSWD     ! incoming short-wave radiation (W/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PSWU     ! upward short-wave radiation (W/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PLWD     ! incoming long-wave radiation (W/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PLWU     ! upward long-wave radiation (W/m2)
REAL, DIMENSION(:,:), INTENT(OUT) :: PSWBD  ! incoming short-wave radiation by spectral band (W/m2)
REAL, DIMENSION(:,:), INTENT(OUT) :: PSWBU  ! upward short-wave radiation by spectral band (W/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PFMU     ! zonal momentum flux (Pa)
REAL, DIMENSION(:), INTENT(OUT) :: PFMV     ! meridian momentum flux (Pa)
REAL, DIMENSION(:), INTENT(OUT) :: PRNC     ! Net radiation       (J/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PHC      ! Sensible heat flux  (J/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PLEC     ! Total latent heat flux    (J/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PLEIC    ! Sublimation latent heat flux    (J/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PGFLUXC  ! Storage flux        (J/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PEVAPC   ! Total evapotranspiration  (kg/m2/s)
REAL, DIMENSION(:), INTENT(OUT) :: PSUBLC   ! Sublimation (kg/m2/s)
REAL, DIMENSION(:), INTENT(OUT) :: PSWDC    ! incoming short wave radiation (J/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PSWUC    ! outgoing short wave radiation (J/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PLWDC    ! incoming long wave radiation (J/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PLWUC    ! outgoing long wave radiation (J/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PFMUC    ! zonal friction
REAL, DIMENSION(:), INTENT(OUT) :: PFMVC    ! meridian friction
REAL, DIMENSION(:), INTENT(OUT) :: PT2M_MIN ! Minimum temperature at 2m   (K)
REAL, DIMENSION(:), INTENT(OUT) :: PT2M_MAX ! Maximum temperature at 2m   (K)
REAL, DIMENSION(:), INTENT(OUT) :: PHU2M_MIN! Minimum relative humidity at 2m (-)
REAL, DIMENSION(:), INTENT(OUT) :: PHU2M_MAX! Maximum relative humidity at 2m (-)
REAL, DIMENSION(:), INTENT(OUT) :: PWIND10M ! wind at 10m (m/s)
REAL, DIMENSION(:), INTENT(OUT) :: PWIND10M_MAX! Maximum wind at 10m (m/s)
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DIAG_ISBA_N',0,ZHOOK_HANDLE)
IF (DGI%LSURF_BUDGET) THEN
  PRN      = DGI%XRN
  PH       = DGI%XH
  PLE      = DGI%XLE
  PLEI     = DGI%XLEI
  PGFLUX   = DGI%XGFLUX
  PSWD     = DGI%XSWD
  PSWU     = DGI%XSWU
  PLWD     = DGI%XLWD
  PLWU     = DGI%XLWU
  PSWBD    = DGI%XSWBD
  PSWBU    = DGI%XSWBU
  PFMU     = DGI%XFMU
  PFMV     = DGI%XFMV
END IF
!
IF (DGEI%LSURF_EVAP_BUDGET) THEN
  PEVAP    = DGI%XEVAP
  PSUBL    = DGI%XSUBL
ENDIF
!
IF (DGEI%LSURF_BUDGETC) THEN
  PRNC      = DGIC%XRN
  PHC       = DGIC%XH
  PLEC      = DGIC%XLE
  PLEIC     = DGIC%XLEI
  PGFLUXC   = DGIC%XGFLUX
  PEVAPC    = DGIC%XEVAP
  PSUBLC    = DGIC%XSUBL
  PSWDC     = DGIC%XSWD
  PSWUC     = DGIC%XSWU
  PLWDC     = DGIC%XLWD
  PLWUC     = DGIC%XLWU
  PFMUC     = DGIC%XFMU
  PFMVC     = DGIC%XFMV
END IF
!
IF (DGI%N2M>=1 .OR. DGI%LSURF_BUDGET .OR. DGEI%LSURF_BUDGETC) PTS = DGI%XTS
!
IF (DGI%N2M>=1) THEN
  PRI      = DGI%XRI
  PT2M     = DGI%XT2M
  PT2M_MIN = DGI%XT2M_MIN
  PT2M_MAX = DGI%XT2M_MAX
  PQ2M     = DGI%XQ2M
  PHU2M    = DGI%XHU2M
  PHU2M_MIN= DGI%XHU2M_MIN
  PHU2M_MAX= DGI%XHU2M_MAX
  PZON10M  = DGI%XZON10M
  PMER10M  = DGI%XMER10M
  PWIND10M = DGI%XWIND10M
  PWIND10M_MAX = DGI%XWIND10M_MAX
END IF
!
IF (DGI%LCOEF) THEN
  PCD      = DGI%XCD
  PCH      = DGI%XCH
  PCE      = DGI%XCE
  PZ0      = DGI%XZ0
  PZ0H     = DGI%XZ0H
END IF
!
IF (DGI%LSURF_VARS) THEN
  PQS = DGI%XQS
ENDIF
!
IF (LHOOK) CALL DR_HOOK('DIAG_ISBA_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_ISBA_n
