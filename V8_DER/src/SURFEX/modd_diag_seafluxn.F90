!     ######################
      MODULE MODD_DIAG_SEAFLUX_n
!     ######################
!
!!****  *MODD_DIAG_SEAFLUX - declaration of diagnostics for SEAFLUX scheme
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2004
!!      Modified    01/2006 : sea flux parameterization.
!!      S.Senesi    01/2014 : add diags on seaice 
!
!*       0.   DECLARATIONS
!             ------------
!
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE DIAG_SEAFLUX_t
!------------------------------------------------------------------------------
!
  REAL    :: XDIAG_TSTEP  ! time step for diagnostics writing
!
  INTEGER :: N2M          ! flag for 2 meters (and 10 meters) quantities
  LOGICAL :: L2M_MIN_ZS   ! flag for 2 meters quantities evaluated on
!                         ! the minimum orographyy of the grid      
  LOGICAL :: LSURF_BUDGET ! flag for surface energy budget
  LOGICAL :: LRAD_BUDGET  ! flag for radiative energy budget
  LOGICAL :: LCOEF        ! flag for transfer coefficients
  LOGICAL :: LSURF_VARS   ! flag for surface variables
  LOGICAL :: LSURF_BUDGETC       ! flag for surface cumulated energy budget
  LOGICAL :: LRESET_BUDGETC      ! flag for surface cumulated energy budget  
!
!* averaged variables
!
  REAL, POINTER, DIMENSION(:)   :: XRI      ! Bulk-Richardson number           (-)
  REAL, POINTER, DIMENSION(:)   :: XRI_ICE  ! Seaice Bulk-Richardson number    (-)
  REAL, POINTER, DIMENSION(:)   :: XCD      ! drag coefficient for wind        (W/s2)
  REAL, POINTER, DIMENSION(:)   :: XCD_ICE  ! Seaice drag coefficient for wind (W/s2)
  REAL, POINTER, DIMENSION(:)   :: XCH      ! drag coefficient for heat        (W/s)
  REAL, POINTER, DIMENSION(:)   :: XCH_ICE  ! Seaice drag coefficient for heat (W/s)
  REAL, POINTER, DIMENSION(:)   :: XCE      ! drag coefficient for vapor       (W/s/K)
  REAL, POINTER, DIMENSION(:)   :: XZ0      ! roughness length for momentum    (m)
  REAL, POINTER, DIMENSION(:)   :: XZ0_ICE  ! Seaice roughness length for momentum (m)
  REAL, POINTER, DIMENSION(:)   :: XZ0H     ! roughness length for heat        (m)
  REAL, POINTER, DIMENSION(:)   :: XZ0H_ICE ! Seaice roughness length for heat (m)
  REAL, POINTER, DIMENSION(:)   :: XRN      ! net radiation at surface         (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XRN_ICE  ! Seaice net radiation at surface  (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XH       ! sensible heat flux               (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XH_ICE   ! Seaice  sensible heat flux       (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XLE      ! total latent heat flux           (W/m2) 
  REAL, POINTER, DIMENSION(:)   :: XLE_ICE     ! sublimation latent heat flux     (W/m2) 
  REAL, POINTER, DIMENSION(:)   :: XGFLUX   ! net soil-vegetation flux         (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XGFLUX_ICE ! net soil-vegetation flux (seaice) (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XEVAP    ! total evaporation                (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XSUBL    ! sublimation                      (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XT2M     ! air temperature at 2 meters      (K)
  REAL, POINTER, DIMENSION(:)   :: XT2M_ICE ! Seaice air temperature at 2 meters (K)
  REAL, POINTER, DIMENSION(:)   :: XT2M_MIN ! Minimum air temperature at 2 meters (K)
  REAL, POINTER, DIMENSION(:)   :: XT2M_MAX ! Maximum air temperature at 2 meters (K)
  REAL, POINTER, DIMENSION(:)   :: XQ2M     ! air humidity at 2 meters         (kg/kg)
  REAL, POINTER, DIMENSION(:)   :: XQ2M_ICE ! Seaice air humidity at 2 meters  (kg/kg)
  REAL, POINTER, DIMENSION(:)   :: XHU2M    ! air relative humidity at 2 meters(-)
  REAL, POINTER, DIMENSION(:)   :: XHU2M_ICE! Seaice air relative humidity at 2 meters(-)
  REAL, POINTER, DIMENSION(:)   :: XHU2M_MIN! Minimum relative humidity at 2 meters (-)
  REAL, POINTER, DIMENSION(:)   :: XHU2M_MAX! Maximum relative humidity at 2 meters (-)
  REAL, POINTER, DIMENSION(:)   :: XQS      ! air humidity at surface          (kg/kg)
  REAL, POINTER, DIMENSION(:)   :: XQS_ICE  ! Seaice air humidity at surface   (kg/kg)
  REAL, POINTER, DIMENSION(:)   :: XZON10M  ! zonal wind at 10 meters          (m/s)
  REAL, POINTER, DIMENSION(:)   :: XZON10M_ICE ! Seaice zonal wind at 10 meters(m/s)
  REAL, POINTER, DIMENSION(:)   :: XMER10M  ! meridian wind at 10 meters       (m/s)
  REAL, POINTER, DIMENSION(:)   :: XMER10M_ICE ! Seaice meridian wind at 10 meters (m/s)
  REAL, POINTER, DIMENSION(:)   :: XWIND10M ! wind at 10 meters                (m/s)
  REAL, POINTER, DIMENSION(:)   :: XWIND10M_ICE ! Seaice wind at 10 meters     (m/s)
  REAL, POINTER, DIMENSION(:)   :: XWIND10M_MAX! Maximum wind at 10 meters     (m/s)
  REAL, POINTER, DIMENSION(:)   :: XLWD     ! downward long wave radiation     (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XLWU     ! upward long wave radiation       (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XLWU_ICE ! Seaice upward long wave radiation (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XSWD     ! downward short wave radiation    (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XSWU     ! upward short wave radiation      (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XSWU_ICE ! Seaice upward short wave radiation (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XSWBD    ! downward short wave radiation by spectral band (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XSWBU    ! upward short wave radiation by spectral band (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XSWBU_ICE! Seaice upward short wave radiation by spectral band (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XFMU     ! horizontal momentum flux zonal    (kg/ms2)
  REAL, POINTER, DIMENSION(:)   :: XFMU_ICE ! Seaice horizontal momentum flux zonal (kg/ms2)
  REAL, POINTER, DIMENSION(:)   :: XFMV     ! horizontal momentum flux meridian (kg/ms2)
  REAL, POINTER, DIMENSION(:)   :: XFMV_ICE ! Seaice horizontal momentum flux meridian (kg/ms2)
!
  REAL, POINTER, DIMENSION(:)   :: XTS     ! surface temperature              (K)
  REAL, POINTER, DIMENSION(:)   :: XTSRAD  ! surface radiative temperature    (K)
  REAL, POINTER, DIMENSION(:)   :: XALBT   ! Total Albedo  
!
!* cumulated averaged variables
!
  REAL, POINTER, DIMENSION(:)   :: XRNC     ! net radiation at surface         (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XRNC_ICE ! Seaice net radiation at surface  (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XHC      ! sensible heat flux               (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XHC_ICE  ! Seaice sensible heat flux        (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XLEC     ! total latent heat flux           (J/m2) 
  REAL, POINTER, DIMENSION(:)   :: XLEC_ICE    ! sublimation latent heat flux     (J/m2) 
  REAL, POINTER, DIMENSION(:)   :: XGFLUXC  ! net soil-vegetation flux         (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XGFLUXC_ICE !Seaice net soil-vegetation flux(J/m2)
  REAL, POINTER, DIMENSION(:)   :: XEVAPC   ! total evaporation                (kg/m2)
  REAL, POINTER, DIMENSION(:)   :: XSUBLC   ! sublimation                      (kg/m2)
  REAL, POINTER, DIMENSION(:)   :: XLWDC    ! downward long wave radiation     (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XLWUC    ! upward long wave radiation       (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XLWUC_ICE! Seaice upward long wave radiation(J/m2)
  REAL, POINTER, DIMENSION(:)   :: XSWDC    ! downward short wave radiation    (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XSWUC    ! upward short wave radiation      (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XSWUC_ICE! Seaice upward short wave radiation(J/m2)
  REAL, POINTER, DIMENSION(:)   :: XFMUC    ! horizontal momentum flux zonal    (kg/ms)
  REAL, POINTER, DIMENSION(:)   :: XFMUC_ICE! Seaice horizontal momentum flux zonal (kg/ms)
  REAL, POINTER, DIMENSION(:)   :: XFMVC    ! horizontal momentum flux meridian (kg/ms)
  REAL, POINTER, DIMENSION(:)   :: XFMVC_ICE! Seaice horizontal momentum flux meridian (kg/ms)
!
!------------------------------------------------------------------------------
!
END TYPE DIAG_SEAFLUX_t

TYPE (DIAG_SEAFLUX_t), ALLOCATABLE, TARGET, SAVE :: DIAG_SEAFLUX_MODEL(:)

TYPE(DIAG_SEAFLUX_t), POINTER :: DIAG_SEAFLUX => NULL()
!$OMP THREADPRIVATE(DIAG_SEAFLUX)

CONTAINS

SUBROUTINE DIAG_SEAFLUX_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_DIAG_SEAFLUX_N:DIAG_SEAFLUX_GOTO_MODEL',0,ZHOOK_HANDLE)

DIAG_SEAFLUX => DIAG_SEAFLUX_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_DIAG_SEAFLUX_N:DIAG_SEAFLUX_GOTO_MODEL',1,ZHOOK_HANDLE)


END SUBROUTINE DIAG_SEAFLUX_GOTO_MODEL

SUBROUTINE DIAG_SEAFLUX_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_SEAFLUX_N:DIAG_SEAFLUX_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(DIAG_SEAFLUX_MODEL(KMODEL))
DIAG_SEAFLUX => DIAG_SEAFLUX_MODEL(KMODEL)
DO J=1,KMODEL
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XRI)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XRI_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XCD)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XCD_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XCH)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XCH_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XCE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XZ0)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XZ0_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XZ0H)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XZ0H_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XRN)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XRN_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XH)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XH_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XLE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XLE_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XGFLUX)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XGFLUX_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XEVAP)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XSUBL)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XT2M)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XT2M_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XT2M_MIN)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XT2M_MAX)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XQ2M)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XQ2M_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XHU2M)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XHU2M_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XHU2M_MIN)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XHU2M_MAX)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XQS)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XQS_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XZON10M)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XZON10M_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XMER10M)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XMER10M_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XWIND10M)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XWIND10M_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XWIND10M_MAX)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XLWD)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XLWU)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XLWU_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XSWD)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XSWU)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XSWU_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XSWBD)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XSWBU)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XSWBU_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XFMU)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XFMU_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XFMV)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XFMV_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XTS)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XTSRAD)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XALBT)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XRNC)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XRNC_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XHC)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XHC_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XLEC)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XLEC_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XGFLUXC)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XGFLUXC_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XEVAPC)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XSUBLC)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XLWDC)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XLWUC)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XLWUC_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XSWDC)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XSWUC)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XSWUC_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XFMUC)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XFMUC_ICE)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XFMVC)
  NULLIFY(DIAG_SEAFLUX_MODEL(J)%XFMVC_ICE)
ENDDO
DIAG_SEAFLUX_MODEL(:)%XDIAG_TSTEP=0.
DIAG_SEAFLUX_MODEL(:)%N2M=0
DIAG_SEAFLUX_MODEL(:)%L2M_MIN_ZS=.FALSE.
DIAG_SEAFLUX_MODEL(:)%LSURF_BUDGET=.FALSE.
DIAG_SEAFLUX_MODEL(:)%LRAD_BUDGET=.FALSE.
DIAG_SEAFLUX_MODEL(:)%LCOEF=.FALSE.
DIAG_SEAFLUX_MODEL(:)%LSURF_VARS=.FALSE.
DIAG_SEAFLUX_MODEL(:)%LSURF_BUDGETC=.FALSE.
DIAG_SEAFLUX_MODEL(:)%LRESET_BUDGETC=.FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_SEAFLUX_N:DIAG_SEAFLUX_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_SEAFLUX_ALLOC

SUBROUTINE DIAG_SEAFLUX_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_SEAFLUX_N:DIAG_SEAFLUX_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(DIAG_SEAFLUX_MODEL)) DEALLOCATE(DIAG_SEAFLUX_MODEL)
IF (ASSOCIATED(DIAG_SEAFLUX)) NULLIFY(DIAG_SEAFLUX)
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_SEAFLUX_N:DIAG_SEAFLUX_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_SEAFLUX_DEALLO

END MODULE MODD_DIAG_SEAFLUX_n
