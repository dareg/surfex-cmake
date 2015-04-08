!     ######################
      MODULE MODD_DIAG_WATFLUX_n
!     ######################
!
!!****  *MODD_DIAG_WATFLUX - declaration of diagnostics for WATFLUX scheme
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

TYPE DIAG_WATFLUX_t
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
  REAL, POINTER, DIMENSION(:)   :: XCD      ! drag coefficient for wind        (W/s2)
  REAL, POINTER, DIMENSION(:)   :: XCH      ! drag coefficient for heat        (W/s)
  REAL, POINTER, DIMENSION(:)   :: XCE      ! drag coefficient for vapor       (W/s/K)
  REAL, POINTER, DIMENSION(:)   :: XZ0      ! roughness length for momentum    (m)
  REAL, POINTER, DIMENSION(:)   :: XZ0H     ! roughness length for heat        (m)      
  REAL, POINTER, DIMENSION(:)   :: XRN      ! net radiation at surface         (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XH       ! sensible heat flux               (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XLE      ! total latent heat flux           (W/m2) 
  REAL, POINTER, DIMENSION(:)   :: XLEI     ! sublimation latent heat flux     (W/m2) 
  REAL, POINTER, DIMENSION(:)   :: XGFLUX   ! net soil-vegetation flux         (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XEVAP    ! total evaporation                (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XSUBL    ! sublimation                      (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XT2M     ! air temperature at 2 meters      (K)
  REAL, POINTER, DIMENSION(:)   :: XT2M_MIN ! Minimum air temperature at 2 meters      (K)
  REAL, POINTER, DIMENSION(:)   :: XT2M_MAX ! Maximum air temperature at 2 meters      (K)
  REAL, POINTER, DIMENSION(:)   :: XQ2M     ! air humidity at 2 meters         (kg/kg)
  REAL, POINTER, DIMENSION(:)   :: XHU2M    ! air relative humidity at 2 meters(-)
  REAL, POINTER, DIMENSION(:)   :: XHU2M_MIN! Minimum relative humidity at 2 meters    (-)
  REAL, POINTER, DIMENSION(:)   :: XHU2M_MAX! Maximum relative humidity at 2 meters    (-)
  REAL, POINTER, DIMENSION(:)   :: XQS      ! air humidity at surface          (kg/kg)
  REAL, POINTER, DIMENSION(:)   :: XZON10M  ! zonal wind at 10 meters          (m/s)
  REAL, POINTER, DIMENSION(:)   :: XMER10M  ! meridian wind at 10 meters       (m/s)
  REAL, POINTER, DIMENSION(:)   :: XWIND10M ! wind at 10 meters                (m/s)
  REAL, POINTER, DIMENSION(:)   :: XWIND10M_MAX! Maximum wind at 10 meters     (m/s)
  REAL, POINTER, DIMENSION(:)   :: XLWD     ! downward long wave radiation     (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XLWU     ! upward long wave radiation       (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XSWD     ! downward short wave radiation    (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XSWU     ! upward short wave radiation      (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XSWBD    ! downward short wave radiation by spectral band   (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XSWBU    ! upward short wave radiation by spectral band (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XFMU     ! horizontal momentum flux zonal   (kg/ms2)
  REAL, POINTER, DIMENSION(:)   :: XFMV     ! horizontal momentum flux meridian (kg/ms2)
  REAL, POINTER, DIMENSION(:)   :: XDIAG_TS ! water surface temperature (K)
!
!* cumulated averaged variables
!
  REAL, POINTER, DIMENSION(:)   :: XRNC     ! net radiation at surface         (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XHC      ! sensible heat flux               (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XLEC     ! total latent heat flux           (J/m2) 
  REAL, POINTER, DIMENSION(:)   :: XLEIC    ! sublimation latent heat flux     (J/m2) 
  REAL, POINTER, DIMENSION(:)   :: XGFLUXC  ! net soil-vegetation flux         (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XEVAPC   ! total evaporation                (kg/m2)
  REAL, POINTER, DIMENSION(:)   :: XSUBLC   ! sublimation                      (kg/m2)
  REAL, POINTER, DIMENSION(:)   :: XLWDC    ! downward long wave radiation     (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XLWUC    ! upward long wave radiation       (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XSWDC    ! downward short wave radiation    (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XSWUC    ! upward short wave radiation      (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XFMUC    ! horizontal momentum flux zonal    (kg/ms)
  REAL, POINTER, DIMENSION(:)   :: XFMVC    ! horizontal momentum flux meridian (kg/ms)
!  
!------------------------------------------------------------------------------
!

END TYPE DIAG_WATFLUX_t

TYPE(DIAG_WATFLUX_t), ALLOCATABLE, TARGET, SAVE :: DIAG_WATFLUX_MODEL(:)

TYPE(DIAG_WATFLUX_t), POINTER :: DIAG_WATFLUX => NULL()
!$OMP THREADPRIVATE(DIAG_WATFLUX)

CONTAINS

SUBROUTINE DIAG_WATFLUX_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_DIAG_WATFLUX_N:DIAG_WATFLUX_GOTO_MODEL',0,ZHOOK_HANDLE)

DIAG_WATFLUX => DIAG_WATFLUX_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_DIAG_WATFLUX_N:DIAG_WATFLUX_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE DIAG_WATFLUX_GOTO_MODEL

SUBROUTINE DIAG_WATFLUX_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_WATFLUX_N:DIAG_WATFLUX_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(DIAG_WATFLUX_MODEL(KMODEL))
DO J=1,KMODEL
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XRI)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XCD)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XCH)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XCE)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XZ0)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XZ0H)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XRN)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XH)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XLE)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XLEI)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XGFLUX)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XEVAP)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XSUBL)  
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XT2M)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XT2M_MIN)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XT2M_MAX)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XQ2M)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XHU2M)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XHU2M_MIN)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XHU2M_MAX)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XQS)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XZON10M)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XMER10M)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XWIND10M)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XWIND10M_MAX)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XLWD)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XLWU)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XSWD)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XSWU)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XSWBD)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XSWBU)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XFMU)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XFMV)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XDIAG_TS)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XRNC)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XHC)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XLEC)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XLEIC)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XGFLUXC)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XEVAPC)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XSUBLC)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XLWDC)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XLWUC)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XSWDC)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XSWUC)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XFMUC)
  NULLIFY(DIAG_WATFLUX_MODEL(J)%XFMVC)
ENDDO
DIAG_WATFLUX_MODEL(:)%XDIAG_TSTEP=0.
DIAG_WATFLUX_MODEL(:)%N2M=0
DIAG_WATFLUX_MODEL(:)%L2M_MIN_ZS=.FALSE.
DIAG_WATFLUX_MODEL(:)%LSURF_BUDGET=.FALSE.
DIAG_WATFLUX_MODEL(:)%LRAD_BUDGET=.FALSE.
DIAG_WATFLUX_MODEL(:)%LCOEF=.FALSE.
DIAG_WATFLUX_MODEL(:)%LSURF_VARS=.FALSE.
DIAG_WATFLUX_MODEL(:)%LSURF_BUDGETC=.FALSE.
DIAG_WATFLUX_MODEL(:)%LRESET_BUDGETC=.FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_WATFLUX_N:DIAG_WATFLUX_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_WATFLUX_ALLOC

SUBROUTINE DIAG_WATFLUX_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_WATFLUX_N:DIAG_WATFLUX_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(DIAG_WATFLUX_MODEL)) DEALLOCATE(DIAG_WATFLUX_MODEL)
IF (ASSOCIATED(DIAG_WATFLUX)) NULLIFY(DIAG_WATFLUX)
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_WATFLUX_N:DIAG_WATFLUX_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_WATFLUX_DEALLO

END MODULE MODD_DIAG_WATFLUX_n
