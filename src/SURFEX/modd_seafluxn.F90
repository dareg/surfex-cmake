!     #################
      MODULE MODD_SEAFLUX_n
!     #################
!
!!****  *MODD_SEAFLUX_n - declaration of surface parameters for an inland water surface
!!
!!    PURPOSE
!!    -------
!     Declaration of surface parameters
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
!!      Original      01/2004
!!      S. Senesi     01/2014  adapt to fractional seaice, and to seaice scheme
!!      S. Belamari   03/2014  Include NZ0
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_TYPE_DATE_SURF
!
USE MODD_TYPES_GLT,   ONLY : T_GLT
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE SEAFLUX_t
!
! General surface: 
!
  REAL, POINTER, DIMENSION(:)   :: XZS     ! orography
  REAL, POINTER, DIMENSION(:,:) :: XCOVER  ! fraction of each ecosystem       (-)
  LOGICAL, POINTER, DIMENSION(:):: LCOVER  ! GCOVER(i)=T --> ith cover field is not 0.
  LOGICAL                       :: LSBL    ! T: SBL scheme between sea and atm. forcing level
!                                          ! F: no atmospheric layers below forcing level      
  LOGICAL                       :: LHANDLE_SIC ! T: we do weight seaice and open sea fluxes
  CHARACTER(LEN=6)              :: CSEAICE_SCHEME! Name of the seaice scheme 
  REAL, POINTER, DIMENSION(:)   :: XSEABATHY   ! bathymetry
!
  LOGICAL                       :: LINTERPOL_SST ! Quadratic interpolation of monthly SST
  CHARACTER(LEN=6)              :: CINTERPOL_SST ! Quadratic interpolation of monthly SST
  LOGICAL                       :: LINTERPOL_SSS ! Quadratic interpolation of monthly SSS
  CHARACTER(LEN=6)              :: CINTERPOL_SSS ! Quadratic interpolation of monthly SSS
  LOGICAL                       :: LINTERPOL_SIC ! Quadratic interpolation of monthly SIC
  CHARACTER(LEN=6)              :: CINTERPOL_SIC ! Quadratic interpolation of monthly SIC
  LOGICAL                       :: LINTERPOL_SIT ! Quadratic interpolation of monthly SIT
  CHARACTER(LEN=6)              :: CINTERPOL_SIT ! Quadratic interpolation of monthly SIT
  REAL                          :: XFREEZING_SST ! Value marking frozen sea in SST data
  REAL                          :: XSIC_EFOLDING_TIME ! For damping of SIC (days)
  REAL                          :: XSIT_EFOLDING_TIME ! For damping of SIT (days)
  REAL                          :: XSEAICE_TSTEP ! Sea ice model time step
  REAL                          :: XCD_ICE_CST   ! Turbulent exchange coefficient for seaice
  REAL                          :: XSI_FLX_DRV   ! Derivative of fluxes on seaice w.r.t to the temperature (W m-2 K-1)
  
!
! Type of formulation for the fluxes
!
  CHARACTER(LEN=6)                  :: CSEA_FLUX   ! type of flux computation
  CHARACTER(LEN=4)                  :: CSEA_ALB    ! type of albedo
  LOGICAL                           :: LPWG        ! flag for gust
  LOGICAL                           :: LPRECIP     ! flag for precip correction
  LOGICAL                           :: LPWEBB      ! flag for Webb correction
  INTEGER                           :: NZ0         ! set to 0,1 or 2 according to Z0 formulation
                                                   ! 0= ARPEGE / 1= Smith (1988) / 2= Direct
  INTEGER                           :: NGRVWAVES   ! set to 0,1 or 2 according to the 
                                                   ! gravity waves model used in coare30_flux
  REAL                              :: XICHCE      ! CE coef calculation for ECUME
  LOGICAL                           :: LPERTFLUX   ! flag for stochastic flux perturbation
!
! Sea/Ocean:
!
  REAL, POINTER, DIMENSION(:) :: XSST    ! sea surface temperature
  REAL, POINTER, DIMENSION(:) :: XSSS    ! sea surface salinity
  REAL, POINTER, DIMENSION(:) :: XTICE   ! sea ice temperature
  REAL, POINTER, DIMENSION(:) :: XSIC    ! sea ice concentration ( constraint for seaice scheme )
  REAL, POINTER, DIMENSION(:) :: XSST_INI! initial sea surface temperature
  REAL, POINTER, DIMENSION(:) :: XZ0     ! roughness length
  REAL, POINTER, DIMENSION(:) :: XZ0H    ! roughness length for heat
  REAL, POINTER, DIMENSION(:) :: XEMIS   ! emissivity
  REAL, POINTER, DIMENSION(:) :: XDIR_ALB! direct albedo
  REAL, POINTER, DIMENSION(:) :: XSCA_ALB! diffuse albedo
  REAL, POINTER, DIMENSION(:) :: XICE_ALB! sea-ice albedo from seaice model (ESM or embedded)
  REAL, POINTER, DIMENSION(:) :: XUMER   ! U component of sea current (for ESM coupling)
  REAL, POINTER, DIMENSION(:) :: XVMER   ! V component of sea current (for ESM coupling)
!
  REAL, POINTER, DIMENSION(:,:) :: XSST_MTH! monthly sea surface temperature (precedent, current and next)
  REAL, POINTER, DIMENSION(:,:) :: XSSS_MTH! monthly sea surface salinity    (precedent, current and next)
  REAL, POINTER, DIMENSION(:,:) :: XSIC_MTH! monthly sea ice cover           (precedent, current and next)
  REAL, POINTER, DIMENSION(:,:) :: XSIT_MTH! monthly sea ice thickness       (precedent, current and next)
  REAL, POINTER, DIMENSION(:)   :: XFSIC   ! nudging (or forcing) sea ice cover
  REAL, POINTER, DIMENSION(:)   :: XFSIT   ! nudging sea ice thickness
!
  REAL, POINTER, DIMENSION(:) :: XCPL_SEA_WIND ! 10m wind speed for ESM coupling
  REAL, POINTER, DIMENSION(:) :: XCPL_SEA_FWSU ! zonal wind stress for ESM coupling
  REAL, POINTER, DIMENSION(:) :: XCPL_SEA_FWSV ! meridian wind stress for ESM coupling
  REAL, POINTER, DIMENSION(:) :: XCPL_SEA_SNET ! Solar net heat flux
  REAL, POINTER, DIMENSION(:) :: XCPL_SEA_HEAT ! Non solar net heat flux
  REAL, POINTER, DIMENSION(:) :: XCPL_SEA_EVAP ! Evaporation for ESM coupling
  REAL, POINTER, DIMENSION(:) :: XCPL_SEA_RAIN ! Rainfall for ESM coupling
  REAL, POINTER, DIMENSION(:) :: XCPL_SEA_SNOW ! Snowfall for ESM coupling
  REAL, POINTER, DIMENSION(:) :: XCPL_SEA_FWSM ! wind stress for ESM coupling
!  
  REAL, POINTER, DIMENSION(:) :: XCPL_SEAICE_SNET ! Solar net heat flux for ESM coupling
  REAL, POINTER, DIMENSION(:) :: XCPL_SEAICE_HEAT ! Non solar net heat flux
  REAL, POINTER, DIMENSION(:) :: XCPL_SEAICE_EVAP ! Ice sublimation for ESM coupling
!
  REAL, POINTER, DIMENSION(:) :: XPERTFLUX     ! Stochastic flux perturbation pattern
!
! Sea-ice :
!
  TYPE(T_GLT)                        :: TGLT ! Sea-ice state , diagnostics and auxilliaries
                                             ! for the case of embedded Gelato Seaice model
!
! Date:
!
  TYPE (DATE_TIME)                      :: TTIME   ! current date and time
  TYPE (DATE_TIME)                      :: TZTIME  
  LOGICAL                               :: LTZTIME_DONE
  INTEGER                               :: JSX
!
! Time-step:
!
  REAL                                  :: XTSTEP  ! time step
!
  REAL                                  :: XOUT_TSTEP  ! output writing time step
!
!
!
END TYPE SEAFLUX_t

TYPE(SEAFLUX_t), ALLOCATABLE, TARGET, SAVE :: SEAFLUX_MODEL(:)

TYPE(SEAFLUX_t), POINTER :: SEAFLUX => NULL()
!$OMP THREADPRIVATE(SEAFLUX)

CONTAINS

SUBROUTINE SEAFLUX_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_SEAFLUX_N:SEAFLUX_GOTO_MODEL',0,ZHOOK_HANDLE)

SEAFLUX => SEAFLUX_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_SEAFLUX_N:SEAFLUX_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE SEAFLUX_GOTO_MODEL

SUBROUTINE SEAFLUX_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_SEAFLUX_N:SEAFLUX_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(SEAFLUX_MODEL(KMODEL))
DO J=1,KMODEL
  NULLIFY(SEAFLUX_MODEL(J)%XZS)
  NULLIFY(SEAFLUX_MODEL(J)%XCOVER)
  NULLIFY(SEAFLUX_MODEL(J)%LCOVER)
  NULLIFY(SEAFLUX_MODEL(J)%XSEABATHY)
  NULLIFY(SEAFLUX_MODEL(J)%XSST)
  NULLIFY(SEAFLUX_MODEL(J)%XSSS)
  NULLIFY(SEAFLUX_MODEL(J)%XSIC)
  NULLIFY(SEAFLUX_MODEL(J)%XTICE)
  NULLIFY(SEAFLUX_MODEL(J)%XSST_INI)
  NULLIFY(SEAFLUX_MODEL(J)%XZ0)
  NULLIFY(SEAFLUX_MODEL(J)%XZ0H)
  NULLIFY(SEAFLUX_MODEL(J)%XEMIS)
  NULLIFY(SEAFLUX_MODEL(J)%XDIR_ALB)
  NULLIFY(SEAFLUX_MODEL(J)%XSCA_ALB)
  NULLIFY(SEAFLUX_MODEL(J)%XICE_ALB)
  NULLIFY(SEAFLUX_MODEL(J)%XUMER)
  NULLIFY(SEAFLUX_MODEL(J)%XVMER)
  NULLIFY(SEAFLUX_MODEL(J)%XSST_MTH)
  NULLIFY(SEAFLUX_MODEL(J)%XSSS_MTH)
  NULLIFY(SEAFLUX_MODEL(J)%XSIC_MTH)
  NULLIFY(SEAFLUX_MODEL(J)%XSIT_MTH)
  NULLIFY(SEAFLUX_MODEL(J)%XFSIC)
  NULLIFY(SEAFLUX_MODEL(J)%XFSIT)
  NULLIFY(SEAFLUX_MODEL(J)%XCPL_SEA_WIND)
  NULLIFY(SEAFLUX_MODEL(J)%XCPL_SEA_FWSU)
  NULLIFY(SEAFLUX_MODEL(J)%XCPL_SEA_FWSV)
  NULLIFY(SEAFLUX_MODEL(J)%XCPL_SEA_SNET)
  NULLIFY(SEAFLUX_MODEL(J)%XCPL_SEA_HEAT)
  NULLIFY(SEAFLUX_MODEL(J)%XCPL_SEA_EVAP)
  NULLIFY(SEAFLUX_MODEL(J)%XCPL_SEA_RAIN)
  NULLIFY(SEAFLUX_MODEL(J)%XCPL_SEA_SNOW)
  NULLIFY(SEAFLUX_MODEL(J)%XCPL_SEA_FWSM)
  NULLIFY(SEAFLUX_MODEL(J)%XCPL_SEAICE_SNET)
  NULLIFY(SEAFLUX_MODEL(J)%XCPL_SEAICE_HEAT)
  NULLIFY(SEAFLUX_MODEL(J)%XCPL_SEAICE_EVAP)
  NULLIFY(SEAFLUX_MODEL(J)%XPERTFLUX)
ENDDO
SEAFLUX_MODEL(:)%LSBL=.FALSE.
SEAFLUX_MODEL(:)%LHANDLE_SIC=.FALSE.
SEAFLUX_MODEL(:)%CSEAICE_SCHEME='NONE  '
SEAFLUX_MODEL(:)%LINTERPOL_SST=.FALSE.
SEAFLUX_MODEL(:)%CINTERPOL_SST=' '
SEAFLUX_MODEL(:)%LINTERPOL_SSS=.FALSE.
SEAFLUX_MODEL(:)%CINTERPOL_SSS=' '
SEAFLUX_MODEL(:)%LINTERPOL_SIC=.FALSE.
SEAFLUX_MODEL(:)%CINTERPOL_SIC=' '
SEAFLUX_MODEL(:)%LINTERPOL_SIT=.FALSE.
SEAFLUX_MODEL(:)%CINTERPOL_SIT=' '
SEAFLUX_MODEL(:)%XFREEZING_SST=-1.8
SEAFLUX_MODEL(:)%XSIC_EFOLDING_TIME=0. ! means : no damping
SEAFLUX_MODEL(:)%XSIT_EFOLDING_TIME=0. ! means : no damping
SEAFLUX_MODEL(:)%XSEAICE_TSTEP=XUNDEF 
SEAFLUX_MODEL(:)%XCD_ICE_CST=0.
SEAFLUX_MODEL(:)%XSI_FLX_DRV=-20. 
SEAFLUX_MODEL(:)%CSEA_FLUX=' '
SEAFLUX_MODEL(:)%CSEA_ALB=' '
SEAFLUX_MODEL(:)%LPWG=.FALSE.
SEAFLUX_MODEL(:)%LPRECIP=.FALSE.
SEAFLUX_MODEL(:)%LPWEBB=.FALSE.
SEAFLUX_MODEL(:)%NZ0=0
SEAFLUX_MODEL(:)%NGRVWAVES=0
SEAFLUX_MODEL(:)%XICHCE=0.
SEAFLUX_MODEL(:)%LPERTFLUX=.FALSE.
SEAFLUX_MODEL(:)%JSX=0
SEAFLUX_MODEL(:)%LTZTIME_DONE = .FALSE.
SEAFLUX_MODEL(:)%XTSTEP=0.
SEAFLUX_MODEL(:)%XOUT_TSTEP=0.
IF (LHOOK) CALL DR_HOOK("MODD_SEAFLUX_N:SEAFLUX_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE SEAFLUX_ALLOC

SUBROUTINE SEAFLUX_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_SEAFLUX_N:SEAFLUX_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(SEAFLUX_MODEL)) DEALLOCATE(SEAFLUX_MODEL)
IF (ASSOCIATED(SEAFLUX)) NULLIFY(SEAFLUX)
IF (LHOOK) CALL DR_HOOK("MODD_SEAFLUX_N:SEAFLUX_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE SEAFLUX_DEALLO

END MODULE MODD_SEAFLUX_n
