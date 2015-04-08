!     ######################
      MODULE MODD_DIAG_TEB_n
!     ######################
!
!!****  *MODD_DIAG_TEB - declaration of diagnostics for TEB scheme
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

TYPE DIAG_TEB_t
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
!
  LOGICAL :: LPGD         ! flag for writing of PGD files
  LOGICAL :: LPGD_FIX     ! flag for writing of PGD files for time
!                           invariant field
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
  REAL, POINTER, DIMENSION(:)   :: XLE      ! latent heat flux                 (W/m2) 
  REAL, POINTER, DIMENSION(:)   :: XGFLUX   ! net soil-vegetation flux         (W/m2)
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
  REAL, POINTER, DIMENSION(:)   :: XSFCO2   ! CO2 flux                         (m/s*kg_CO2/kg_air)
  REAL, POINTER, DIMENSION(:)   :: XLWD     ! downward long wave radiation     (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XLWU     ! upward long wave radiation       (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XSWD     ! downward short wave radiation    (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XSWU     ! upward short wave radiation      (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XSWBD    ! downward short wave radiation by spectral band   (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XSWBU    ! upward short wave radiation by spectral band (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XFMU     ! horizontal momentum flux zonal   (m2/s2)
  REAL, POINTER, DIMENSION(:)   :: XFMV     ! horizontal momentum flux meridian (m2/s2)             
  REAL, POINTER, DIMENSION(:)   :: XDIAG_TS ! arithmetic mean of surface temperature (K)
!------------------------------------------------------------------------------
!

END TYPE DIAG_TEB_t

TYPE(DIAG_TEB_t), ALLOCATABLE, TARGET, SAVE :: DIAG_TEB_MODEL(:)

TYPE(DIAG_TEB_t), POINTER :: DIAG_TEB => NULL()
!$OMP THREADPRIVATE(DIAG_TEB)

CONTAINS

SUBROUTINE DIAG_TEB_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_DIAG_TEB_N:DIAG_TEB_GOTO_MODEL',0,ZHOOK_HANDLE)

DIAG_TEB => DIAG_TEB_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_DIAG_TEB_N:DIAG_TEB_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE DIAG_TEB_GOTO_MODEL

SUBROUTINE DIAG_TEB_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_TEB_N:DIAG_TEB_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(DIAG_TEB_MODEL(KMODEL))
DO J=1,KMODEL
  NULLIFY(DIAG_TEB_MODEL(J)%XRI)
  NULLIFY(DIAG_TEB_MODEL(J)%XCD)
  NULLIFY(DIAG_TEB_MODEL(J)%XCH)
  NULLIFY(DIAG_TEB_MODEL(J)%XCE)
  NULLIFY(DIAG_TEB_MODEL(J)%XZ0)
  NULLIFY(DIAG_TEB_MODEL(J)%XZ0H)
  NULLIFY(DIAG_TEB_MODEL(J)%XRN)
  NULLIFY(DIAG_TEB_MODEL(J)%XH)
  NULLIFY(DIAG_TEB_MODEL(J)%XLE)
  NULLIFY(DIAG_TEB_MODEL(J)%XGFLUX)
  NULLIFY(DIAG_TEB_MODEL(J)%XT2M)
  NULLIFY(DIAG_TEB_MODEL(J)%XT2M_MIN)
  NULLIFY(DIAG_TEB_MODEL(J)%XT2M_MAX)
  NULLIFY(DIAG_TEB_MODEL(J)%XQ2M)
  NULLIFY(DIAG_TEB_MODEL(J)%XHU2M)
  NULLIFY(DIAG_TEB_MODEL(J)%XHU2M_MIN)
  NULLIFY(DIAG_TEB_MODEL(J)%XHU2M_MAX)
  NULLIFY(DIAG_TEB_MODEL(J)%XQS)
  NULLIFY(DIAG_TEB_MODEL(J)%XZON10M)
  NULLIFY(DIAG_TEB_MODEL(J)%XMER10M)
  NULLIFY(DIAG_TEB_MODEL(J)%XWIND10M)
  NULLIFY(DIAG_TEB_MODEL(J)%XWIND10M_MAX)
  NULLIFY(DIAG_TEB_MODEL(J)%XSFCO2)
  NULLIFY(DIAG_TEB_MODEL(J)%XLWD)
  NULLIFY(DIAG_TEB_MODEL(J)%XLWU)
  NULLIFY(DIAG_TEB_MODEL(J)%XSWD)
  NULLIFY(DIAG_TEB_MODEL(J)%XSWU)
  NULLIFY(DIAG_TEB_MODEL(J)%XSWBD)
  NULLIFY(DIAG_TEB_MODEL(J)%XSWBU)
  NULLIFY(DIAG_TEB_MODEL(J)%XFMU)
  NULLIFY(DIAG_TEB_MODEL(J)%XFMV)
  NULLIFY(DIAG_TEB_MODEL(J)%XDIAG_TS)
ENDDO
DIAG_TEB_MODEL(:)%XDIAG_TSTEP=0.
DIAG_TEB_MODEL(:)%N2M=0
DIAG_TEB_MODEL(:)%L2M_MIN_ZS=.FALSE.
DIAG_TEB_MODEL(:)%LSURF_BUDGET=.FALSE.
DIAG_TEB_MODEL(:)%LRAD_BUDGET=.FALSE.
DIAG_TEB_MODEL(:)%LCOEF=.FALSE.
DIAG_TEB_MODEL(:)%LSURF_VARS=.FALSE.
DIAG_TEB_MODEL(:)%LPGD=.FALSE.
DIAG_TEB_MODEL(:)%LPGD_FIX=.FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_TEB_N:DIAG_TEB_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_TEB_ALLOC

SUBROUTINE DIAG_TEB_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_TEB_N:DIAG_TEB_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(DIAG_TEB_MODEL)) DEALLOCATE(DIAG_TEB_MODEL)
IF (ASSOCIATED(DIAG_TEB)) NULLIFY(DIAG_TEB)
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_TEB_N:DIAG_TEB_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_TEB_DEALLO

END MODULE MODD_DIAG_TEB_n
