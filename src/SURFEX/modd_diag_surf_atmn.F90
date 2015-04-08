!     ######################
      MODULE MODD_DIAG_SURF_ATM_n
!     ######################
!
!!****  *MODD_DIAG_SURF_ATM - declaration of diagnostics for the surface
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
!!      Modified    04/2009 : precip for/from restart file.
!!      Modified    08/2009 : BUDGETC for all tiles
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TYPE_DATE_SURF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE DIAG_SURF_ATM_t
!------------------------------------------------------------------------------
!
  REAL    :: XDIAG_TSTEP  ! time step for diagnostics writing
!
  INTEGER :: N2M          ! flag for 2 meters (and 10 meters) quantities
  LOGICAL :: LT2MMW       ! flag to perform modified weighting of 2m temperature
  LOGICAL :: L2M_MIN_ZS   ! flag for 2 meters quantities evaluated on
!                         ! the minimum orographyy of the grid 
  LOGICAL :: LSURF_BUDGET ! flag for surface energy budget
  LOGICAL :: LRAD_BUDGET  ! flag for radiative energy budget      
  LOGICAL :: LCOEF        ! flag for transfer coefficients
  LOGICAL :: LSURF_VARS   ! flag for surface variables
  LOGICAL :: LFRAC        ! flag for writing fractions of each four tiles
  LOGICAL :: LDIAG_GRID   ! flag for mean grid diag
  LOGICAL :: LSURF_BUDGETC       ! flag for surface cumulated energy budget
  LOGICAL :: LRESET_BUDGETC      ! flag for surface cumulated energy budget
  LOGICAL :: LREAD_BUDGETC       ! flag for surface cumulated energy budget
  LOGICAL :: LPROVAR_TO_DIAG     ! switch to write (or not) prognostic variable
                                 ! and allows puting field in diagnostics 
  LOGICAL    :: LSELECT          ! switch to control which fields are written
                                 ! (only those whose naem appears in in text array)
!  
  TYPE(DATE_TIME):: TIME_BUDGETC
!                                  
  CHARACTER(LEN=12), POINTER, DIMENSION(:) :: CSELECT  ! Name of ouput fields if LSELECT=true
!
!* variables for each tile
!
  REAL, POINTER, DIMENSION(:,:) :: XRI_TILE     ! Bulk-Richardson number           (-)
  REAL, POINTER, DIMENSION(:,:) :: XCD_TILE     ! drag coefficient for wind        (W/s2)
  REAL, POINTER, DIMENSION(:,:) :: XCH_TILE     ! drag coefficient for heat        (W/s)
  REAL, POINTER, DIMENSION(:,:) :: XCE_TILE     ! drag coefficient for vapor       (W/s/K)
  REAL, POINTER, DIMENSION(:,:) :: XRN_TILE     ! net radiation at surface         (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XH_TILE      ! sensible heat flux               (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XLE_TILE     ! total latent heat flux           (W/m2) 
  REAL, POINTER, DIMENSION(:,:) :: XLEI_TILE    ! sublimation latent heat flux     (W/m2) 
  REAL, POINTER, DIMENSION(:,:) :: XGFLUX_TILE  ! net soil-vegetation flux         (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XEVAP_TILE   ! total evapotranspiration         (kg/m2/s) 
  REAL, POINTER, DIMENSION(:,:) :: XSUBL_TILE   ! sublimation                      (kg/m2/s) 
  REAL, POINTER, DIMENSION(:,:) :: XTS_TILE     ! surface temperature              (K)
  REAL, POINTER, DIMENSION(:,:) :: XT2M_TILE    ! air temperature at 2 meters      (K)
  REAL, POINTER, DIMENSION(:,:) :: XT2M_MIN_TILE! Minimum air temperature at 2 meters (K)
  REAL, POINTER, DIMENSION(:,:) :: XT2M_MAX_TILE! Maximum air temperature at 2 meters (K)
  REAL, POINTER, DIMENSION(:,:) :: XQ2M_TILE    ! air humidity at 2 meters         (kg/kg)
  REAL, POINTER, DIMENSION(:,:) :: XHU2M_TILE   ! air relative humidity at 2 meters(-)
  REAL, POINTER, DIMENSION(:,:) :: XHU2M_MIN_TILE! Minimum air relative humidity at 2 meters(-)
  REAL, POINTER, DIMENSION(:,:) :: XHU2M_MAX_TILE! Maximum air relative humidity at 2 meters(-)
  REAL, POINTER, DIMENSION(:,:) :: XZON10M_TILE ! zonal wind at 10 meters          (m/s)
  REAL, POINTER, DIMENSION(:,:) :: XMER10M_TILE ! meridian wind at 10 meters       (m/s)
  REAL, POINTER, DIMENSION(:,:) :: XWIND10M_TILE! wind at 10 meters                (m/s)
  REAL, POINTER, DIMENSION(:,:) :: XWIND10M_MAX_TILE ! Maximumwind at 10 meters    (m/s)
  REAL, POINTER, DIMENSION(:,:) :: XQS_TILE
  REAL, POINTER, DIMENSION(:,:) :: XZ0_TILE     ! roughness length for momentum    (m)
  REAL, POINTER, DIMENSION(:,:) :: XZ0H_TILE    ! roughness length for heat        (m)
  REAL, POINTER, DIMENSION(:,:) :: XSWD_TILE    ! short wave downward radiation (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XSWU_TILE    ! short wave upward radiation (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XLWD_TILE    ! longt wave downward radiation (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XLWU_TILE    ! longt wave upward radiation (W/m2)
  REAL, POINTER, DIMENSION(:,:,:) :: XSWBD_TILE ! short wave downward radiation by spectral band(W/m2)
  REAL, POINTER, DIMENSION(:,:,:) :: XSWBU_TILE ! short wave upward radiation by spectral band(W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XFMU_TILE    ! zonal friction
  REAL, POINTER, DIMENSION(:,:) :: XFMV_TILE    ! meridian friction
!
!* Cumulated variables for each tile
!
  REAL, POINTER, DIMENSION(:,:) :: XRNC_TILE     ! net radiation at surface         (J/m2)
  REAL, POINTER, DIMENSION(:,:) :: XHC_TILE      ! sensible heat flux               (J/m2)
  REAL, POINTER, DIMENSION(:,:) :: XLEC_TILE     ! total latent heat flux           (J/m2) 
  REAL, POINTER, DIMENSION(:,:) :: XLEIC_TILE    ! sublimation latent heat flux     (J/m2) 
  REAL, POINTER, DIMENSION(:,:) :: XGFLUXC_TILE  ! net soil-vegetation flux         (J/m2)
  REAL, POINTER, DIMENSION(:,:) :: XEVAPC_TILE   ! total evapotranspiration         (kg/m2)
  REAL, POINTER, DIMENSION(:,:) :: XSUBLC_TILE   ! sublimation                      (kg/m2) 
  REAL, POINTER, DIMENSION(:,:) :: XSWDC_TILE    ! short wave downward radiation    (J/m2)
  REAL, POINTER, DIMENSION(:,:) :: XSWUC_TILE    ! short wave upward radiation      (J/m2)
  REAL, POINTER, DIMENSION(:,:) :: XLWDC_TILE    ! longt wave downward radiation    (J/m2)
  REAL, POINTER, DIMENSION(:,:) :: XLWUC_TILE    ! longt wave upward radiation      (J/m2)
  REAL, POINTER, DIMENSION(:,:) :: XFMUC_TILE    ! zonal friction
  REAL, POINTER, DIMENSION(:,:) :: XFMVC_TILE    ! meridian friction
!
!* averaged variables
!
  REAL, POINTER, DIMENSION(:)   :: XAVG_RI      ! Bulk-Richardson number           (-)
  REAL, POINTER, DIMENSION(:)   :: XAVG_CD       ! drag coefficient for wind        (W/s2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_CH       ! drag coefficient for heat        (W/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_CE       ! drag coefficient for vapor       (W/s/K)
  REAL, POINTER, DIMENSION(:)   :: XAVG_RN      ! net radiation at surface         (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_H       ! sensible heat flux               (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LE      ! total latent heat flux           (W/m2) 
  REAL, POINTER, DIMENSION(:)   :: XAVG_LEI     ! sublimation latent heat flux     (W/m2) 
  REAL, POINTER, DIMENSION(:)   :: XAVG_GFLUX   ! net soil-vegetation flux         (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_EVAP    ! total evapotranspiration         (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_SUBL    ! sublimation                      (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_TS      ! surface temperature              (K)
  REAL, POINTER, DIMENSION(:)   :: XAVG_T2M     ! air temperature at 2 meters      (K)
  REAL, POINTER, DIMENSION(:)   :: XAVG_Q2M     ! air humidity at 2 meters         (kg/kg)
  REAL, POINTER, DIMENSION(:)   :: XAVG_HU2M    ! air relative humidity at 2 meters(-)
  REAL, POINTER, DIMENSION(:)   :: XAVG_ZON10M  ! zonal wind at 10 meters          (m/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_MER10M  ! meridian wind at 10 meters       (m/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_SFCO2   ! CO2 flux                         (m/s*kg_CO2/kg_air)  
  REAL, POINTER, DIMENSION(:)   :: XAVG_T2M_MIN_ZS ! air temperature at 2 meters   (K)
  REAL, POINTER, DIMENSION(:)   :: XAVG_Q2M_MIN_ZS ! air humidity at 2 meters      (kg/kg)
  REAL, POINTER, DIMENSION(:)   :: XAVG_HU2M_MIN_ZS! air relative humidity at 2 m  (-)
  REAL, POINTER, DIMENSION(:)   :: XPS          ! air pressure at the surface      (Pa)
  REAL, POINTER, DIMENSION(:)   :: XRHOA        ! air density  at the surface      (kg/m3)
  REAL, POINTER, DIMENSION(:)   :: XAVG_QS
  REAL, POINTER, DIMENSION(:)   :: XAVG_Z0      ! roughness length for momentum    (m)
  REAL, POINTER, DIMENSION(:)   :: XAVG_Z0H     ! roughness length for heat        (m)
  REAL, POINTER, DIMENSION(:)   :: XAVG_SWD     ! short wave downward radiation (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_SWU     ! short wave upward radiation (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LWD     ! longt wave downward radiation (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LWU     ! longt wave upward radiation (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XAVG_SWBD    ! short wave downward radiation by spectral band(W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XAVG_SWBU    ! short wave upward radiation by spectral band(W/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_FMU     ! zonal friction
  REAL, POINTER, DIMENSION(:)   :: XAVG_FMV     ! meridian friction
  REAL, POINTER, DIMENSION(:)   :: XSSO_FMU     ! zonal friction    (with SSO)     (Pa)
  REAL, POINTER, DIMENSION(:)   :: XSSO_FMV     ! meridian friction (with SSO)     (Pa)
!
  REAL, POINTER, DIMENSION(:)   :: XDIAG_UREF   ! reference height for momentum    (m)
  REAL, POINTER, DIMENSION(:)   :: XDIAG_ZREF   ! reference height for heat        (m)
  REAL, POINTER, DIMENSION(:)   :: XDIAG_TRAD   ! radiative temperature at t       (K)
  REAL, POINTER, DIMENSION(:)   :: XDIAG_EMIS   ! surface emissivity at t          (-)
!
  REAL, POINTER, DIMENSION(:)   :: XAVG_T2M_MIN ! Minimun air temperature at 2 meters      (K)
  REAL, POINTER, DIMENSION(:)   :: XAVG_T2M_MAX ! Maximum air temperature at 2 meters      (K)
  REAL, POINTER, DIMENSION(:)   :: XAVG_HU2M_MIN! Minimun air relative humidity at 2 meters(-)
  REAL, POINTER, DIMENSION(:)   :: XAVG_HU2M_MAX! Maximum air relative humidity at 2 meters(-)
  REAL, POINTER, DIMENSION(:)   :: XAVG_WIND10M ! wind at 10 meters                      (m/s)
  REAL, POINTER, DIMENSION(:)   :: XAVG_WIND10M_MAX ! Maximum wind at 10 meters          (m/s)
!
!* cumulated averaged variables
!
  REAL, POINTER, DIMENSION(:)   :: XAVG_RNC      ! net radiation at surface         (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_HC       ! sensible heat flux               (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LEC      ! total latent heat flux           (J/m2) 
  REAL, POINTER, DIMENSION(:)   :: XAVG_LEIC     ! sublimation latent heat flux     (J/m2) 
  REAL, POINTER, DIMENSION(:)   :: XAVG_GFLUXC   ! net soil-vegetation flux         (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_EVAPC    ! total evapotranspiration         (kg/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_SUBLC    ! sublimation                      (kg/m2)  
  REAL, POINTER, DIMENSION(:)   :: XAVG_SWDC     ! short wave downward radiation    (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_SWUC     ! short wave upward radiation      (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LWDC     ! longt wave downward radiation    (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_LWUC     ! longt wave upward radiation      (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XAVG_FMUC     ! zonal friction
  REAL, POINTER, DIMENSION(:)   :: XAVG_FMVC     ! meridian friction
!
!------------------------------------------------------------------------------
!

END TYPE DIAG_SURF_ATM_t

TYPE(DIAG_SURF_ATM_t), ALLOCATABLE, TARGET, SAVE :: DIAG_SURF_ATM_MODEL(:)

TYPE(DIAG_SURF_ATM_t), POINTER :: DIAG_SURF_ATM => NULL()
!$OMP THREADPRIVATE(DIAG_SURF_ATM)

CONTAINS

SUBROUTINE DIAG_SURF_ATM_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_DIAG_SURF_ATM_N:DIAG_SURF_ATM_GOTO_MODEL',0,ZHOOK_HANDLE)

DIAG_SURF_ATM => DIAG_SURF_ATM_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_DIAG_SURF_ATM_N:DIAG_SURF_ATM_GOTO_MODEL',1,ZHOOK_HANDLE)
!
END SUBROUTINE DIAG_SURF_ATM_GOTO_MODEL

SUBROUTINE DIAG_SURF_ATM_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_SURF_ATM_N:DIAG_SURF_ATM_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(DIAG_SURF_ATM_MODEL(KMODEL))
DO J=1,KMODEL
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XRI_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XCD_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XCH_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XCE_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XRN_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XH_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XLE_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XLEI_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XGFLUX_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XEVAP_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XSUBL_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XTS_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XT2M_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XT2M_MIN_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XT2M_MAX_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XQ2M_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XHU2M_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XHU2M_MIN_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XHU2M_MAX_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XZON10M_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XMER10M_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XWIND10M_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XWIND10M_MAX_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XQS_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XZ0_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XZ0H_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XSWD_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XSWU_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XLWD_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XLWU_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XSWBD_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XSWBU_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XFMU_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XFMV_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XRNC_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XHC_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XLEC_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XLEIC_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XGFLUXC_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XEVAPC_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XSUBLC_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XSWDC_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XSWUC_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XLWDC_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XLWUC_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XFMUC_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XFMVC_TILE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_RI)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_CD)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_CH)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_CE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_RN)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_H)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_LE)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_LEI)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_GFLUX)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_EVAP)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_SUBL)  
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_TS)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_T2M)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_Q2M)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_HU2M)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_ZON10M)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_MER10M)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_SFCO2)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_T2M_MIN_ZS)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_Q2M_MIN_ZS)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_HU2M_MIN_ZS)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XPS)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XRHOA)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_QS)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_Z0)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_Z0H)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_SWD)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_SWU)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_LWD)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_LWU)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_SWBD)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_SWBU)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_FMU)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_FMV)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XSSO_FMU)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XSSO_FMV)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XDIAG_UREF)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XDIAG_ZREF)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XDIAG_TRAD)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XDIAG_EMIS)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_T2M_MIN)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_T2M_MAX)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_HU2M_MIN)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_HU2M_MAX)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_WIND10M)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_WIND10M_MAX)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_RNC)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_HC)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_LEC)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_LEIC)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_GFLUXC)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_EVAPC)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_SUBLC)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_SWDC)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_SWUC)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_LWDC)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_LWUC)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_FMUC)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%XAVG_FMVC)
  NULLIFY(DIAG_SURF_ATM_MODEL(J)%CSELECT)
ENDDO
DIAG_SURF_ATM_MODEL(:)%XDIAG_TSTEP=0.
DIAG_SURF_ATM_MODEL(:)%N2M=0
DIAG_SURF_ATM_MODEL(:)%LT2MMW=.FALSE.
DIAG_SURF_ATM_MODEL(:)%L2M_MIN_ZS=.FALSE.
DIAG_SURF_ATM_MODEL(:)%LSURF_BUDGET=.FALSE.
DIAG_SURF_ATM_MODEL(:)%LRAD_BUDGET=.FALSE.
DIAG_SURF_ATM_MODEL(:)%LCOEF=.FALSE.
DIAG_SURF_ATM_MODEL(:)%LSURF_VARS=.FALSE.
DIAG_SURF_ATM_MODEL(:)%LFRAC=.FALSE.
DIAG_SURF_ATM_MODEL(:)%LDIAG_GRID=.FALSE.
DIAG_SURF_ATM_MODEL(:)%LSURF_BUDGETC=.FALSE.
DIAG_SURF_ATM_MODEL(:)%LRESET_BUDGETC=.FALSE.
DIAG_SURF_ATM_MODEL(:)%LREAD_BUDGETC=.FALSE.
DIAG_SURF_ATM_MODEL(:)%LPROVAR_TO_DIAG=.FALSE.
DIAG_SURF_ATM_MODEL(:)%LSELECT=.FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_SURF_ATM_N:DIAG_SURF_ATM_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_SURF_ATM_ALLOC

SUBROUTINE DIAG_SURF_ATM_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_SURF_ATM_N:DIAG_SURF_ATM_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(DIAG_SURF_ATM_MODEL)) DEALLOCATE(DIAG_SURF_ATM_MODEL)
IF (ASSOCIATED(DIAG_SURF_ATM)) NULLIFY(DIAG_SURF_ATM)
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_SURF_ATM_N:DIAG_SURF_ATM_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_SURF_ATM_DEALLO

END MODULE MODD_DIAG_SURF_ATM_n
