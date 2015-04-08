!######################
MODULE MODD_DIAG_ISBA_n
!######################
!
!!****  *MODD_DIAG_ISBA - declaration of diagnostics for ISBA scheme
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
!!      P. Samuelsson 10/2014 : added min max for XT2M
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

TYPE DIAG_ISBA_t
!------------------------------------------------------------------------------
!
  REAL    :: XDIAG_TSTEP  ! time step for diagnostics writing
!
  INTEGER :: N2M          ! flag for 2 meters (and 10 meters) quantities
  LOGICAL :: L2M_MIN_ZS   ! flag for 2 meters quantities evaluated on
!                         ! the minimum orographyy of the grid      
  LOGICAL :: LSURF_BUDGET   ! flag for surface energy budget
  LOGICAL :: LRAD_BUDGET    ! flag for radiative energy budget
!
  LOGICAL :: LCOEF        ! flag for transfer coefficients
  LOGICAL :: LSURF_VARS   ! flag for surface variables
!
  LOGICAL :: LPGD          ! flag for writing of PGD files
  LOGICAL :: LPATCH_BUDGET ! flag for patch output
!
!* variables for each patch
!
  REAL, POINTER, DIMENSION(:,:) :: XRI     ! Bulk-Richardson number           (-)
  REAL, POINTER, DIMENSION(:,:) :: XCD     ! drag coefficient for wind        (W/s2)
  REAL, POINTER, DIMENSION(:,:) :: XCH     ! drag coefficient for heat        (W/s)
  REAL, POINTER, DIMENSION(:,:) :: XCE     ! drag coefficient for vapor       (W/s/K)
  REAL, POINTER, DIMENSION(:,:) :: XRN     ! net radiation at surface         (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XH      ! sensible heat flux               (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XLE     ! total latent heat flux           (W/m2) 
  REAL, POINTER, DIMENSION(:,:) :: XLEI    ! sublimation latent heat flux     (W/m2) 
  REAL, POINTER, DIMENSION(:,:) :: XGFLUX  ! net soil-vegetation flux         (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XTS     ! surface temperature              (K)
  REAL, POINTER, DIMENSION(:,:) :: XTSRAD  ! surface radiative temperature    (K)  
  REAL, POINTER, DIMENSION(:,:) :: XT2M    ! temperature at 2 meters          (K)
  REAL, POINTER, DIMENSION(:,:) :: XT2M_MIN! min temperature at 2 meters      (K)
  REAL, POINTER, DIMENSION(:,:) :: XT2M_MAX! max temperature at 2 meters      (K)
  REAL, POINTER, DIMENSION(:,:) :: XQ2M    ! humidity    at 2 meters          (kg/kg)
  REAL, POINTER, DIMENSION(:,:) :: XHU2M   ! relative humidity at 2 meters    (-)
  REAL, POINTER, DIMENSION(:,:) :: XQS     ! humidity at surface              (kg/kg)
  REAL, POINTER, DIMENSION(:,:) :: XZON10M ! zonal wind at 10 meters          (m/s)
  REAL, POINTER, DIMENSION(:,:) :: XMER10M ! meridian wind at 10 meters       (m/s)
  REAL, POINTER, DIMENSION(:,:) :: XWIND10M! wind at 10 meters                (m/s)
  REAL, POINTER, DIMENSION(:,:) :: XLWD    ! downward long wave radiation     (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XLWU    ! upward long wave radiation       (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XSWD    ! downward short wave radiation    (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XSWU    ! upward short wave radiation      (W/m2)
  REAL, POINTER, DIMENSION(:,:,:) :: XSWBD ! downward short wave radiation by spectral band   (W/m2)
  REAL, POINTER, DIMENSION(:,:,:) :: XSWBU ! upward short wave radiation by spectral band (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XFMU    ! horizontal momentum flux zonal    (Pa)
  REAL, POINTER, DIMENSION(:,:) :: XFMV    ! horizontal momentum flux meridian (Pa)             
  ! 
  REAL, POINTER, DIMENSION(:,:) :: XSWDC   ! downward short wave radiation     (J/m2)
  REAL, POINTER, DIMENSION(:,:) :: XSWUC   ! upward short wave radiation       (J/m2)
  REAL, POINTER, DIMENSION(:,:) :: XLWDC   ! downward long wave radiation      (J/m2)
  REAL, POINTER, DIMENSION(:,:) :: XLWUC   ! upward long wave radiation        (J/m2)
  REAL, POINTER, DIMENSION(:,:) :: XFMUC   ! horizontal momentum flux zonal    (Pa.s)
  REAL, POINTER, DIMENSION(:,:) :: XFMVC   ! horizontal momentum flux meridian (Pa.s)
  !
  REAL, POINTER, DIMENSION(:,:) :: XZ0_WITH_SNOW  ! roughness length for momentum
                                                  ! for vegetation and snow    (m)
  REAL, POINTER, DIMENSION(:,:) :: XZ0H_WITH_SNOW ! roughness length for heat
                                                  ! for vegetation and snow    (m)
  REAL, POINTER, DIMENSION(:,:) :: XZ0EFF         ! effective roughness length for heat
                                                  ! for vegetation and snow    (m)
!
!* averaged variables
!
  REAL, POINTER, DIMENSION(:)   :: XAVG_RI       ! Bulk-Richardson number           (-)
  REAL, POINTER, DIMENSION(:)   :: XAVG_CD       ! drag coefficient for wind        (W/s2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_CH       ! drag coefficient for heat        (W/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_CE       ! drag coefficient for vapor       (W/s/K)
  REAL, POINTER, DIMENSION(:)   :: XAVG_RN       ! net radiation at surface         (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_H        ! sensible heat flux               (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LE       ! total latent heat flux           (W/m2) 
  REAL, POINTER, DIMENSION(:)   :: XAVG_LEI      ! sublimation latent heat flux     (W/m2) 
  REAL, POINTER, DIMENSION(:)   :: XAVG_GFLUX    ! net soil-vegetation flux         (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_TS       ! surface temperature              (K)
  REAL, POINTER, DIMENSION(:)   :: XAVG_T2M      ! temperature at 2 meters          (K)
  REAL, POINTER, DIMENSION(:)   :: XAVG_T2M_MIN  ! Minimum temperature at 2 meters          (K)
  REAL, POINTER, DIMENSION(:)   :: XAVG_T2M_MAX  ! Maximum temperature at 2 meters          (K)
  REAL, POINTER, DIMENSION(:)   :: XAVG_Q2M      ! humidity    at 2 meters          (kg/kg)
  REAL, POINTER, DIMENSION(:)   :: XAVG_HU2M     ! relative humidity at 2 meters    (-)
  REAL, POINTER, DIMENSION(:)   :: XAVG_HU2M_MIN ! Minimum relative humidity at 2 meters    (-)
  REAL, POINTER, DIMENSION(:)   :: XAVG_HU2M_MAX ! Maximum relative humidity at 2 meters    (-)
  REAL, POINTER, DIMENSION(:)   :: XAVG_QS       ! humidity at surface              (kg/kg)
  REAL, POINTER, DIMENSION(:)   :: XAVG_ZON10M   ! zonal wind at 10 meters          (m/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_MER10M   ! meridian wind at 10 meters       (m/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_WIND10M  ! wind at 10 meters                (m/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_WIND10M_MAX  ! Maximum wind at 10 meters    (m/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_SFCO2    ! CO2 flux                         (m/s*kg_CO2/kg_air)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LWD      ! downward long wave radiation     (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LWU      ! upward long wave radiation       (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_SWD      ! downward short wave radiation    (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_SWU      ! upward short wave radiation      (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XAVG_SWBD     ! downward short wave radiation by spectral band   (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XAVG_SWBU     ! upward short wave radiation by spectral band (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_FMU      ! horizontal momentum flux zonal   (m2/s2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_FMV      ! horizontal momentum flux meridian (m2/s2) 
  REAL, POINTER, DIMENSION(:)   :: XAVG_LWDC     ! downward long wave radiation     (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LWUC     ! upward long wave radiation       (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_SWDC     ! downward short wave radiation    (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_SWUC     ! upward short wave radiation      (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_FMUC     ! horizontal momentum flux zonal   (Pa.s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_FMVC     ! horizontal momentum flux meridian (Pa.s)             
  
  !                                                
  REAL, POINTER, DIMENSION(:)   :: XAVG_Z0       ! roughness length for momentum
                                                 ! for vegetation and snow    (m)
  REAL, POINTER, DIMENSION(:)   :: XAVG_Z0H      ! roughness length for heat
                                                 ! for vegetation and snow    (m)
  REAL, POINTER, DIMENSION(:)   :: XAVG_Z0EFF    ! effective roughness length for heat
                                                 ! for vegetation and snow    (m)
!------------------------------------------------------------------------------
!

END TYPE DIAG_ISBA_t

TYPE(DIAG_ISBA_t), ALLOCATABLE, TARGET, SAVE :: DIAG_ISBA_MODEL(:)

TYPE(DIAG_ISBA_t), POINTER :: DIAG_ISBA => NULL()
!$OMP THREADPRIVATE(DIAG_ISBA)

CONTAINS

SUBROUTINE DIAG_ISBA_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_DIAG_ISBA_N:DIAG_ISBA_GOTO_MODEL',0,ZHOOK_HANDLE)

DIAG_ISBA => DIAG_ISBA_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_DIAG_ISBA_N:DIAG_ISBA_GOTO_MODEL',1,ZHOOK_HANDLE)


END SUBROUTINE DIAG_ISBA_GOTO_MODEL

SUBROUTINE DIAG_ISBA_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_ISBA_N:DIAG_ISBA_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(DIAG_ISBA_MODEL(KMODEL))
DO J=1,KMODEL
  NULLIFY(DIAG_ISBA_MODEL(J)%XRI)
  NULLIFY(DIAG_ISBA_MODEL(J)%XCD)
  NULLIFY(DIAG_ISBA_MODEL(J)%XCH)
  NULLIFY(DIAG_ISBA_MODEL(J)%XCE)
  NULLIFY(DIAG_ISBA_MODEL(J)%XRN)
  NULLIFY(DIAG_ISBA_MODEL(J)%XH)
  NULLIFY(DIAG_ISBA_MODEL(J)%XLE)
  NULLIFY(DIAG_ISBA_MODEL(J)%XLEI)
  NULLIFY(DIAG_ISBA_MODEL(J)%XGFLUX)
  NULLIFY(DIAG_ISBA_MODEL(J)%XTS)
  NULLIFY(DIAG_ISBA_MODEL(J)%XTSRAD)
  NULLIFY(DIAG_ISBA_MODEL(J)%XT2M)
  NULLIFY(DIAG_ISBA_MODEL(J)%XT2M_MIN)
  NULLIFY(DIAG_ISBA_MODEL(J)%XT2M_MAX)
  NULLIFY(DIAG_ISBA_MODEL(J)%XQ2M)
  NULLIFY(DIAG_ISBA_MODEL(J)%XHU2M)
  NULLIFY(DIAG_ISBA_MODEL(J)%XQS)
  NULLIFY(DIAG_ISBA_MODEL(J)%XZON10M)
  NULLIFY(DIAG_ISBA_MODEL(J)%XMER10M)
  NULLIFY(DIAG_ISBA_MODEL(J)%XWIND10M)
  NULLIFY(DIAG_ISBA_MODEL(J)%XLWD)
  NULLIFY(DIAG_ISBA_MODEL(J)%XLWU)
  NULLIFY(DIAG_ISBA_MODEL(J)%XSWD)
  NULLIFY(DIAG_ISBA_MODEL(J)%XSWU)
  NULLIFY(DIAG_ISBA_MODEL(J)%XSWBD)
  NULLIFY(DIAG_ISBA_MODEL(J)%XSWBU)
  NULLIFY(DIAG_ISBA_MODEL(J)%XFMU)
  NULLIFY(DIAG_ISBA_MODEL(J)%XFMV)
  NULLIFY(DIAG_ISBA_MODEL(J)%XSWDC)
  NULLIFY(DIAG_ISBA_MODEL(J)%XSWUC)
  NULLIFY(DIAG_ISBA_MODEL(J)%XLWDC)
  NULLIFY(DIAG_ISBA_MODEL(J)%XLWUC)
  NULLIFY(DIAG_ISBA_MODEL(J)%XFMUC)
  NULLIFY(DIAG_ISBA_MODEL(J)%XFMVC)
  NULLIFY(DIAG_ISBA_MODEL(J)%XZ0_WITH_SNOW)
  NULLIFY(DIAG_ISBA_MODEL(J)%XZ0H_WITH_SNOW)
  NULLIFY(DIAG_ISBA_MODEL(J)%XZ0EFF)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_RI)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_CD)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_CH)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_CE)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_RN)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_H)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_LE)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_LEI)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_GFLUX)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_TS)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_T2M)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_T2M_MIN)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_T2M_MAX)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_Q2M)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_HU2M)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_HU2M_MIN)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_HU2M_MAX)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_QS)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_ZON10M)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_MER10M)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_WIND10M)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_WIND10M_MAX)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_SFCO2)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_LWD)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_LWU)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_SWD)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_SWU)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_SWBD)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_SWBU)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_FMU)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_FMV)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_LWDC)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_LWUC)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_SWDC)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_SWUC)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_FMUC)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_FMVC)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_Z0)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_Z0H)
  NULLIFY(DIAG_ISBA_MODEL(J)%XAVG_Z0EFF)
ENDDO
DIAG_ISBA_MODEL(:)%XDIAG_TSTEP=0.
DIAG_ISBA_MODEL(:)%N2M=0
DIAG_ISBA_MODEL(:)%L2M_MIN_ZS=.FALSE.
DIAG_ISBA_MODEL(:)%LSURF_BUDGET=.FALSE.
DIAG_ISBA_MODEL(:)%LRAD_BUDGET=.FALSE.
DIAG_ISBA_MODEL(:)%LCOEF=.FALSE.
DIAG_ISBA_MODEL(:)%LSURF_VARS=.FALSE.
DIAG_ISBA_MODEL(:)%LPGD=.FALSE.
DIAG_ISBA_MODEL(:)%LPATCH_BUDGET=.FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_ISBA_N:DIAG_ISBA_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_ISBA_ALLOC

SUBROUTINE DIAG_ISBA_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_ISBA_N:DIAG_ISBA_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(DIAG_ISBA_MODEL)) DEALLOCATE(DIAG_ISBA_MODEL)
IF (ASSOCIATED(DIAG_ISBA)) NULLIFY(DIAG_ISBA)
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_ISBA_N:DIAG_ISBA_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_ISBA_DEALLO

END MODULE MODD_DIAG_ISBA_n
