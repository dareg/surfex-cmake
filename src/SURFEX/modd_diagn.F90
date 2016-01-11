!######################
MODULE MODD_DIAG_n
!######################
!
!!****  *MODD_DIAG - declaration of diagnostics for ISBA scheme
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
USE MODD_TYPE_DATE_SURF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE DIAG_t
!------------------------------------------------------------------------------
!
  REAL    :: XDIAG_TSTEP  ! time step for diagnostics writing
!
  INTEGER :: N2M          ! flag for 2 meters (and 10 meters) quantities
  LOGICAL :: LT2MMW       ! flag to perform modified weighting of 2m temperature  
  LOGICAL :: L2M_MIN_ZS   ! flag for 2 meters quantities evaluated on
!                         ! the minimum orographyy of the grid      
  LOGICAL :: LSURF_BUDGET   ! flag for surface energy budget
  LOGICAL :: LRAD_BUDGET    ! flag for radiative energy budget
!
  LOGICAL :: LCOEF        ! flag for transfer coefficients
  LOGICAL :: LSURF_VARS   ! flag for surface variables
  LOGICAL :: LFRAC        ! flag for writing fractions of each four tiles
  LOGICAL :: LDIAG_GRID   ! flag for mean grid diag
!  
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
  LOGICAL :: LPGD          ! flag for writing of PGD files
  LOGICAL :: LPATCH_BUDGET ! flag for patch output
!
!* variables for each patch
!
!* averaged variables
!
  REAL, POINTER, DIMENSION(:)   :: XRI       ! Bulk-Richardson number           (-)
  REAL, POINTER, DIMENSION(:)   :: XCD       ! drag coefficient for wind        (W/s2)
  REAL, POINTER, DIMENSION(:)   :: XCDN       ! neutral drag coefficient                      (-)  
  REAL, POINTER, DIMENSION(:)   :: XCH       ! drag coefficient for heat        (W/s)
  REAL, POINTER, DIMENSION(:)   :: XCE       ! drag coefficient for vapor       (W/s/K)
!
  REAL, POINTER, DIMENSION(:)   :: XHU        ! area averaged surface humidity coefficient    (-)
  REAL, POINTER, DIMENSION(:)   :: XHUG       ! baresoil surface humidity coefficient         (-)
  REAL, POINTER, DIMENSION(:)   :: XHV        ! Halstead coefficient                          (-)  
!
  REAL, POINTER, DIMENSION(:)   :: XRN       ! net radiation at surface         (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XH        ! sensible heat flux               (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XLE       ! total latent heat flux           (W/m2) 
  REAL, POINTER, DIMENSION(:)   :: XLEI      ! sublimation latent heat flux     (W/m2) 
  REAL, POINTER, DIMENSION(:)   :: XGFLUX    ! net soil-vegetation flux         (W/m2)
!
  REAL, POINTER, DIMENSION(:)   :: XEVAP    ! total evaporation                (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XSUBL    ! sublimation                      (kg/m2/s)
!
  REAL, POINTER, DIMENSION(:)   :: XTS       ! surface temperature              (K)
  REAL, POINTER, DIMENSION(:)   :: XTSRAD    ! surface temperature              (K)
  REAL, POINTER, DIMENSION(:)   :: XALBT          ! Total Albedo  
  REAL, POINTER, DIMENSION(:)   :: XSWE     ! snow water equivalent (kg/m2)
!  
  REAL, POINTER, DIMENSION(:)   :: XT2M      ! temperature at 2 meters          (K)
  REAL, POINTER, DIMENSION(:)   :: XT2M_MIN  ! Minimum temperature at 2 meters          (K)
  REAL, POINTER, DIMENSION(:)   :: XT2M_MAX  ! Maximum temperature at 2 meters          (K)
  REAL, POINTER, DIMENSION(:)   :: XQ2M      ! humidity    at 2 meters          (kg/kg)
  REAL, POINTER, DIMENSION(:)   :: XHU2M     ! relative humidity at 2 meters    (-)
  REAL, POINTER, DIMENSION(:)   :: XHU2M_MIN ! Minimum relative humidity at 2 meters    (-)
  REAL, POINTER, DIMENSION(:)   :: XHU2M_MAX ! Maximum relative humidity at 2 meters    (-)
  REAL, POINTER, DIMENSION(:)   :: XQS       ! humidity at surface              (kg/kg)
  REAL, POINTER, DIMENSION(:)   :: XZON10M   ! zonal wind at 10 meters          (m/s)
  REAL, POINTER, DIMENSION(:)   :: XMER10M   ! meridian wind at 10 meters       (m/s)
  REAL, POINTER, DIMENSION(:)   :: XWIND10M  ! wind at 10 meters                (m/s)
  REAL, POINTER, DIMENSION(:)   :: XWIND10M_MAX  ! Maximum wind at 10 meters    (m/s)
!
  REAL, POINTER, DIMENSION(:)   :: XSFCO2    ! CO2 flux                         (m/s*kg_CO2/kg_air)
!
  REAL, POINTER, DIMENSION(:,:) :: XSWBD     ! downward short wave radiation by spectral band   (W/m2)
  REAL, POINTER, DIMENSION(:,:) :: XSWBU     ! upward short wave radiation by spectral band (W/m2)
!
  REAL, POINTER, DIMENSION(:)   :: XLWD      ! downward long wave radiation     (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XLWU      ! upward long wave radiation       (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XSWD      ! downward short wave radiation    (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XSWU      ! upward short wave radiation      (W/m2)
!
  REAL, POINTER, DIMENSION(:)   :: XFMU      ! horizontal momentum flux zonal   (m2/s2)
  REAL, POINTER, DIMENSION(:)   :: XFMV      ! horizontal momentum flux meridian (m2/s2) 
  !                                                
  REAL, POINTER, DIMENSION(:)   :: XZ0       ! roughness length for momentum
                                                 ! for vegetation and snow    (m)
  REAL, POINTER, DIMENSION(:)   :: XZ0H      ! roughness length for heat
                                                 ! for vegetation and snow    (m)
  REAL, POINTER, DIMENSION(:)   :: XZ0EFF    ! effective roughness length for heat
                                                 ! for vegetation and snow    (m)
!
  REAL, POINTER, DIMENSION(:)   :: XT2M_MIN_ZS ! air temperature at 2 meters   (K)
  REAL, POINTER, DIMENSION(:)   :: XQ2M_MIN_ZS ! air humidity at 2 meters      (kg/kg)
  REAL, POINTER, DIMENSION(:)   :: XHU2M_MIN_ZS! air relative humidity at 2 m  (-)



  REAL, POINTER, DIMENSION(:)   :: XPS          ! air pressure at the surface      (Pa)
  REAL, POINTER, DIMENSION(:)   :: XRHOA        ! air density  at the surface      (kg/m3)

  REAL, POINTER, DIMENSION(:)   :: XSSO_FMU     ! zonal friction    (with SSO)     (Pa)
  REAL, POINTER, DIMENSION(:)   :: XSSO_FMV     ! meridian friction (with SSO)     (Pa)
!

  REAL, POINTER, DIMENSION(:)   :: XUREF   ! reference height for momentum    (m)
  REAL, POINTER, DIMENSION(:)   :: XZREF   ! reference height for heat        (m)
  REAL, POINTER, DIMENSION(:)   :: XTRAD   ! radiative temperature at t       (K)
  REAL, POINTER, DIMENSION(:)   :: XEMIS   ! surface emissivity at t          (-)

!------------------------------------------------------------------------------
!
END TYPE DIAG_t
!
TYPE DIAG_PATCH_t
!
TYPE(DIAG_t), ALLOCATABLE :: AL(:) 
!
END TYPE DIAG_PATCH_t
!
CONTAINS
!
SUBROUTINE DIAG_PATCH_INIT(YDIAG_PATCH,KPATCH)
TYPE(DIAG_PATCH_t), INTENT(INOUT) :: YDIAG_PATCH 
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_N:DIAG_PATCH_INIT",0,ZHOOK_HANDLE)
 ALLOCATE(YDIAG_PATCH%AL(KPATCH))
DO JP=1,KPATCH
  CALL DIAG_INIT(YDIAG_PATCH%AL(JP))
ENDDO         
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_N:DIAG_PATCH_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_PATCH_INIT
!
!
SUBROUTINE DIAG_INIT(YDIAG)
TYPE(DIAG_t), INTENT(INOUT) :: YDIAG
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_N:DIAG_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YDIAG%XRI)
  NULLIFY(YDIAG%XCD)
  NULLIFY(YDIAG%XCDN)
  NULLIFY(YDIAG%XCH)
  NULLIFY(YDIAG%XCE)
  NULLIFY(YDIAG%XHU)
  NULLIFY(YDIAG%XHUG)
  NULLIFY(YDIAG%XHV)  
  NULLIFY(YDIAG%XRN)
  NULLIFY(YDIAG%XH)
  NULLIFY(YDIAG%XLE)
  NULLIFY(YDIAG%XLEI)
  NULLIFY(YDIAG%XGFLUX)
  NULLIFY(YDIAG%XEVAP)
  NULLIFY(YDIAG%XSUBL)  
  NULLIFY(YDIAG%XTS)
  NULLIFY(YDIAG%XTSRAD)
  NULLIFY(YDIAG%XALBT)  
  NULLIFY(YDIAG%XSWE)  
  NULLIFY(YDIAG%XT2M)
  NULLIFY(YDIAG%XT2M_MIN)
  NULLIFY(YDIAG%XT2M_MAX)
  NULLIFY(YDIAG%XQ2M)
  NULLIFY(YDIAG%XHU2M)
  NULLIFY(YDIAG%XQS)
  NULLIFY(YDIAG%XZON10M)
  NULLIFY(YDIAG%XMER10M)
  NULLIFY(YDIAG%XWIND10M)
  NULLIFY(YDIAG%XWIND10M_MAX)  
  NULLIFY(YDIAG%XLWD)
  NULLIFY(YDIAG%XLWU)
  NULLIFY(YDIAG%XSWD)
  NULLIFY(YDIAG%XSWU)
  NULLIFY(YDIAG%XSWBD)
  NULLIFY(YDIAG%XSWBU)
  NULLIFY(YDIAG%XFMU)
  NULLIFY(YDIAG%XFMV)
  NULLIFY(YDIAG%XZ0)
  NULLIFY(YDIAG%XZ0H)
  NULLIFY(YDIAG%XZ0EFF)
  NULLIFY(YDIAG%CSELECT)
  NULLIFY(YDIAG%XT2M_MIN_ZS)
  NULLIFY(YDIAG%XQ2M_MIN_ZS)
  NULLIFY(YDIAG%XHU2M_MIN_ZS)
  NULLIFY(YDIAG%XPS)
  NULLIFY(YDIAG%XRHOA)
  NULLIFY(YDIAG%XSSO_FMU)
  NULLIFY(YDIAG%XSSO_FMV)
  NULLIFY(YDIAG%XUREF)
  NULLIFY(YDIAG%XZREF)
  NULLIFY(YDIAG%XTRAD)
  NULLIFY(YDIAG%XEMIS)
YDIAG%XDIAG_TSTEP=0.
YDIAG%N2M=0
YDIAG%LT2MMW=.FALSE.
YDIAG%L2M_MIN_ZS=.FALSE.
YDIAG%LSURF_BUDGET=.FALSE.
YDIAG%LRAD_BUDGET=.FALSE.
YDIAG%LCOEF=.FALSE.
YDIAG%LSURF_VARS=.FALSE.
YDIAG%LFRAC=.FALSE.
YDIAG%LDIAG_GRID=.FALSE.
YDIAG%LPGD=.FALSE.
YDIAG%LPATCH_BUDGET=.FALSE.
YDIAG%LSURF_BUDGETC=.FALSE.
YDIAG%LRESET_BUDGETC=.FALSE.
YDIAG%LREAD_BUDGETC=.FALSE.
YDIAG%LPROVAR_TO_DIAG=.FALSE.
YDIAG%LSELECT=.FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_N:DIAG_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_INIT


END MODULE MODD_DIAG_n
