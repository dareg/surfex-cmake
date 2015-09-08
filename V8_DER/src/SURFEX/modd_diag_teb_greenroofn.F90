!######################
MODULE MODD_DIAG_TEB_GREENROOF_n
!######################
!
!!****  *MODD_DIAG_TEB_GREENROOF - declaration of diagnostics for ISBA scheme
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
!!    Based on modd_diag_teb_gardenn
!!
!!    AUTHOR
!!    ------
!!      C. de Munck & A. Lemonsu   *Meteo France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       07/2011
!
!*       0.   DECLARATIONS
!             ------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE DIAG_TEB_GREENROOF_t
!------------------------------------------------------------------------------
!
!* variables for one patch
!
  REAL, POINTER, DIMENSION(:)   :: XRI     ! Bulk-Richardson number           (-)
  REAL, POINTER, DIMENSION(:)   :: XCD     ! drag coefficient for wind        (W/s2)
  REAL, POINTER, DIMENSION(:)   :: XCH     ! drag coefficient for heat        (W/s)
  REAL, POINTER, DIMENSION(:)   :: XCE     ! drag coefficient for vapor       (W/s/K)
  REAL, POINTER, DIMENSION(:)   :: XRN     ! net radiation at surface         (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XH      ! sensible heat flux               (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XGFLUX  ! net soil-vegetation flux         (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XTS     ! surface temperature              (K)
  REAL, POINTER, DIMENSION(:)   :: XTSRAD  ! radiative surface temperature    (K)
  REAL, POINTER, DIMENSION(:)   :: XQS     ! humidity at surface              (kg/kg)
  REAL, POINTER, DIMENSION(:)   :: XLWD    ! downward long wave radiation     (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XLWU    ! upward long wave radiation       (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XSWD    ! downward short wave radiation    (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XSWU    ! upward short wave radiation      (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XSWBD ! downward short wave radiation by spectral band   (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XSWBU ! upward short wave radiation by spectral band (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XFMU    ! horizontal momentum flux zonal   (m2/s2)
  REAL, POINTER, DIMENSION(:)   :: XFMV    ! horizontal momentum flux meridian (m2/s2)  
  !    
  REAL, POINTER, DIMENSION(:)   :: XZ0_WITH_SNOW  ! roughness length for momentum
                                                  ! for vegetation and snow    (m)
  REAL, POINTER, DIMENSION(:)   :: XZ0H_WITH_SNOW ! roughness length for heat
                                                  ! for vegetation and snow    (m)
  REAL, POINTER, DIMENSION(:)   :: XZ0EFF         ! effective roughness length for heat
                                                  ! for vegetation and snow    (m)
!
  REAL, POINTER, DIMENSION(:)   :: XLEI          ! sublimation latent heat flux     (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XLEG          ! latent heat of evaporation over the ground   (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XLEGI         ! surface soil ice sublimation                 (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XLEV          ! latent heat of evaporation over vegetation   (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XLES          ! latent heat of evaporation over the snow     (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XLER          ! evaporation from canopy water interception   (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XLETR         ! evapotranspiration of the vegetation         (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XEVAP         ! evapotranspiration                           (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XDRAIN        ! soil drainage flux                           (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XRUNOFF       ! sub-grid and supersaturation runoff          (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XHORT         ! sub-grid Horton runoff from the SGH scheme   (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XRRVEG        !  precipitation intercepted by the vegetation (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XMELT         ! snow melt                                    (kg/m2/s)       
  REAL, POINTER, DIMENSION(:)   :: XDRIP         ! dripping from the vegetation reservoir       (kg/m2/s)  
!
!* pack diag
!
  REAL, POINTER, DIMENSION(:)   :: XCG           ! heat capacity of the ground
  REAL, POINTER, DIMENSION(:)   :: XC1           ! coefficients for the moisure
  REAL, POINTER, DIMENSION(:)   :: XC2           ! equation.
  REAL, POINTER, DIMENSION(:)   :: XWGEQ         ! equilibrium volumetric water content
  REAL, POINTER, DIMENSION(:)   :: XCT           ! area-averaged heat capacity
  REAL, POINTER, DIMENSION(:)   :: XRS           ! stomatal resistance                            (s/m)
  REAL, POINTER, DIMENSION(:)   :: XCDN          ! neutral drag coefficient                      (-)
  REAL, POINTER, DIMENSION(:)   :: XHU           ! area averaged surface humidity coefficient    (-)
  REAL, POINTER, DIMENSION(:)   :: XHUG          ! baresoil surface humidity coefficient         (-)
  REAL, POINTER, DIMENSION(:)   :: XRESTORE      ! surface energy budget restore term            (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XUSTAR        ! friction velocity                             (m/s)
  REAL, POINTER, DIMENSION(:,:) :: XIACAN        ! PAR in the canopy at different gauss level    (micmolphot/m2/s)
!
! for ISBA-ES:3-L
  REAL, POINTER, DIMENSION(:,:) :: XSNOWTEMP     ! snow temperature profile (ISBA-ES:3-L)        (K)
  REAL, POINTER, DIMENSION(:,:) :: XSNOWLIQ      ! snow liquid water profile (ISBA-ES:3-L)       (m)
  REAL, POINTER, DIMENSION(:,:) :: XSNOWDZ       ! snow layer thicknesses                        (m)
  REAL, POINTER, DIMENSION(:)   :: XSNOWHMASS    ! heat content change due to mass changes in snowpack (J/m2)
  REAL, POINTER, DIMENSION(:)   :: XMELTADV      ! advective energy from snow melt water!
!
!* budget summation variables for one patch
!
!
  REAL, POINTER, DIMENSION(:)   :: XHV           ! Halstead coefficient
!      
  REAL, POINTER, DIMENSION(:,:) :: XSWI        ! Soil wetness index
  REAL, POINTER, DIMENSION(:,:) :: XTSWI       ! Total soil wetness index
!     
  REAL, POINTER, DIMENSION(:)   :: XTWSNOW       ! Total snow reservoir
  REAL, POINTER, DIMENSION(:)   :: XTDSNOW       ! Total snow height
!
  REAL, POINTER, DIMENSION(:)   :: XALBT             ! Total Albedo
  REAL, POINTER, DIMENSION(:)   :: XEMIST            ! averaged emissivity                     (-)
!
  REAL, POINTER, DIMENSION(:)   :: XSEUIL        ! Irrigation threshold
!
  REAL, POINTER, DIMENSION(:)   :: XGPP          ! Gross Primary Production
  REAL, POINTER, DIMENSION(:)   :: XRESP_AUTO    ! Autotrophic respiration
  REAL, POINTER, DIMENSION(:)   :: XRESP_ECO     ! Ecosystem respiration
!
!------------------------------------------------------------------------------
!

END TYPE DIAG_TEB_GREENROOF_t

TYPE(DIAG_TEB_GREENROOF_t), ALLOCATABLE, TARGET, SAVE :: DIAG_TEB_GREENROOF_MODEL(:)

TYPE(DIAG_TEB_GREENROOF_t), POINTER :: DIAG_TEB_GREENROOF => NULL()
!$OMP THREADPRIVATE(DIAG_TEB_GREENROOF)

CONTAINS

SUBROUTINE DIAG_TEB_GREENROOF_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_DIAG_TEB_GREENROOF_N:DIAG_TEB_GREENROOF_GOTO_MODEL',0,ZHOOK_HANDLE)

DIAG_TEB_GREENROOF => DIAG_TEB_GREENROOF_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_DIAG_TEB_GREENROOF_N:DIAG_TEB_GREENROOF_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE DIAG_TEB_GREENROOF_GOTO_MODEL

SUBROUTINE DIAG_TEB_GREENROOF_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_TEB_GREENROOF_N:DIAG_TEB_GREENROOF_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(DIAG_TEB_GREENROOF_MODEL(KMODEL))
DIAG_TEB_GREENROOF => DIAG_TEB_GREENROOF_MODEL(KMODEL)
DO J=1,KMODEL
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XRI)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XCD)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XCH)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XCE)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XRN)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XH)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XGFLUX)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XTS)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XTSRAD)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XQS)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XLWD)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XLWU)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XSWD)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XSWU)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XSWBD)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XSWBU)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XFMU)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XFMV)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XZ0_WITH_SNOW)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XZ0H_WITH_SNOW)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XZ0EFF)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XLEI)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XLEG)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XLEGI)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XLEV)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XLES)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XLER)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XLETR)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XEVAP)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XDRAIN)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XRUNOFF)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XHORT)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XRRVEG)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XMELT)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XDRIP)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XCG)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XC1)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XC2)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XWGEQ)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XCT)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XRS)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XCDN)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XHU)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XHUG)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XRESTORE)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XUSTAR)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XIACAN)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XSNOWTEMP)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XSNOWLIQ)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XSNOWDZ)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XSNOWHMASS)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XMELTADV)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XHV)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XSWI)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XTSWI)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XTWSNOW)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XTDSNOW)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XALBT)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XEMIST)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XSEUIL)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XGPP)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XRESP_AUTO)
  NULLIFY(DIAG_TEB_GREENROOF_MODEL(J)%XRESP_ECO)
ENDDO
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_TEB_GREENROOF_N:DIAG_TEB_GREENROOF_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_TEB_GREENROOF_ALLOC

SUBROUTINE DIAG_TEB_GREENROOF_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_TEB_GREENROOF_N:DIAG_TEB_GREENROOF_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(DIAG_TEB_GREENROOF_MODEL)) DEALLOCATE(DIAG_TEB_GREENROOF_MODEL)
IF (ASSOCIATED(DIAG_TEB_GREENROOF)) NULLIFY(DIAG_TEB_GREENROOF)
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_TEB_GREENROOF_N:DIAG_TEB_GREENROOF_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_TEB_GREENROOF_DEALLO

END MODULE MODD_DIAG_TEB_GREENROOF_n
