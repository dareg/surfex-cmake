!     ######################
      MODULE MODD_DIAG_IDEAL_n
!     ######################
!
!!****  *MODD_DIAG_IDEAL - declaration of diagnostics for IDEAL scheme
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
!!      P. Le Moigne   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       04/2009
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

TYPE DIAG_IDEAL_t
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
  REAL, POINTER, DIMENSION(:)   :: XQ2M     ! air humidity at 2 meters         (kg/kg)
  REAL, POINTER, DIMENSION(:)   :: XHU2M    ! relative humidity at 2 meters    (-)
  REAL, POINTER, DIMENSION(:)   :: XQS      ! humidity at surface              (kg/kg)
  REAL, POINTER, DIMENSION(:)   :: XZON10M  ! zonal wind at 10 meters          (m/s)
  REAL, POINTER, DIMENSION(:)   :: XMER10M  ! meridian wind at 10 meters       (m/s)
  REAL, POINTER, DIMENSION(:)   :: XLWD     ! downward long wave radiation     (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XLWU     ! upward long wave radiation       (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XSWD     ! downward short wave radiation    (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XSWU     ! upward short wave radiation      (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XSWBD    ! downward short wave radiation by spectral band   (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XSWBU    ! upward short wave radiation by spectral band (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XFMU     ! horizontal momentum flux zonal   (m2/s2)
  REAL, POINTER, DIMENSION(:)   :: XFMV     ! horizontal momentum flux meridian (m2/s2)             
!------------------------------------------------------------------------------
!

END TYPE DIAG_IDEAL_t

TYPE(DIAG_IDEAL_t), ALLOCATABLE, TARGET, SAVE :: DIAG_IDEAL_MODEL(:)

TYPE(DIAG_IDEAL_t), POINTER :: DIAG_IDEAL => NULL()
!$OMP THREADPRIVATE(DIAG_IDEAL)

CONTAINS

SUBROUTINE DIAG_IDEAL_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_DIAG_IDEAL_N:DIAG_IDEAL_GOTO_MODEL',0,ZHOOK_HANDLE)

DIAG_IDEAL => DIAG_IDEAL_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_DIAG_IDEAL_N:DIAG_IDEAL_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE DIAG_IDEAL_GOTO_MODEL

SUBROUTINE DIAG_IDEAL_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_IDEAL_N:DIAG_IDEAL_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(DIAG_IDEAL_MODEL(KMODEL))
DIAG_IDEAL => DIAG_IDEAL_MODEL(KMODEL)
DO J=1,KMODEL
  NULLIFY(DIAG_IDEAL_MODEL(J)%XRI)
  NULLIFY(DIAG_IDEAL_MODEL(J)%XCD)
  NULLIFY(DIAG_IDEAL_MODEL(J)%XCH)
  NULLIFY(DIAG_IDEAL_MODEL(J)%XCE)
  NULLIFY(DIAG_IDEAL_MODEL(J)%XZ0)
  NULLIFY(DIAG_IDEAL_MODEL(J)%XZ0H)
  NULLIFY(DIAG_IDEAL_MODEL(J)%XRN)
  NULLIFY(DIAG_IDEAL_MODEL(J)%XH)
  NULLIFY(DIAG_IDEAL_MODEL(J)%XLE)
  NULLIFY(DIAG_IDEAL_MODEL(J)%XGFLUX)
  NULLIFY(DIAG_IDEAL_MODEL(J)%XT2M)
  NULLIFY(DIAG_IDEAL_MODEL(J)%XQ2M)
  NULLIFY(DIAG_IDEAL_MODEL(J)%XHU2M)
  NULLIFY(DIAG_IDEAL_MODEL(J)%XQS)
  NULLIFY(DIAG_IDEAL_MODEL(J)%XZON10M)
  NULLIFY(DIAG_IDEAL_MODEL(J)%XMER10M)
  NULLIFY(DIAG_IDEAL_MODEL(J)%XLWD)
  NULLIFY(DIAG_IDEAL_MODEL(J)%XLWU)
  NULLIFY(DIAG_IDEAL_MODEL(J)%XSWD)
  NULLIFY(DIAG_IDEAL_MODEL(J)%XSWU)
  NULLIFY(DIAG_IDEAL_MODEL(J)%XSWBD)
  NULLIFY(DIAG_IDEAL_MODEL(J)%XSWBU)
  NULLIFY(DIAG_IDEAL_MODEL(J)%XFMU)
  NULLIFY(DIAG_IDEAL_MODEL(J)%XFMV)
ENDDO
DIAG_IDEAL_MODEL(:)%XDIAG_TSTEP=0.
DIAG_IDEAL_MODEL(:)%N2M=0
DIAG_IDEAL_MODEL(:)%L2M_MIN_ZS=.FALSE.
DIAG_IDEAL_MODEL(:)%LSURF_BUDGET=.FALSE.
DIAG_IDEAL_MODEL(:)%LRAD_BUDGET=.FALSE.
DIAG_IDEAL_MODEL(:)%LCOEF=.FALSE.
DIAG_IDEAL_MODEL(:)%LSURF_VARS=.FALSE.
DIAG_IDEAL_MODEL(:)%LSURF_BUDGETC=.FALSE.
DIAG_IDEAL_MODEL(:)%LRESET_BUDGETC=.FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_IDEAL_N:DIAG_IDEAL_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_IDEAL_ALLOC

SUBROUTINE DIAG_IDEAL_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_IDEAL_N:DIAG_IDEAL_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(DIAG_IDEAL_MODEL)) DEALLOCATE(DIAG_IDEAL_MODEL)
IF (ASSOCIATED(DIAG_IDEAL)) NULLIFY(DIAG_IDEAL)
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_IDEAL_N:DIAG_IDEAL_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_IDEAL_DEALLO

END MODULE MODD_DIAG_IDEAL_n
