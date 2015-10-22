!     ############################
      MODULE MODD_DIAG_MISC_ISBA_n
!     ############################
!
!!****  *MODD_DIAG_MISC_ISBA - declaration of packed surface parameters for ISBA scheme
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
!!      Original       07/10/04
!!      A.L. Gibelin 04/2009 : Add respiration diagnostics
!!      A.L. Gibelin 05/2009 : Add carbon spinup
!!      A.L. Gibelin 07/2009 : Suppress RDK and transform GPP as a diagnostic
!!       B. Decharme 05/2012 : Carbon fluxes in diag_evap
!!       B. Decharme 07/2012 : New diag for DIF under LSURF_MISC_DIF key
!!                               F2 stress
!!                               Root zone swi, wg and wgi
!!                               swi, wg and wgi comparable to ISBA-FR-DG2 and DG3 layers
!!       B. Decharme 10/2012 : New diag for DIF 
!!                               active layer thickness over permafrost area
!!                               frozen layer thickness over non-permafrost area
!!       M. Lafaysse 09/2015 : new Crocus-MEPRA outputs
!-------------------------------------------------------------------------------
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

TYPE DIAG_MISC_ISBA_t
!------------------------------------------------------------------------------
!
  LOGICAL :: LSURF_MISC_BUDGET   ! flag for miscellaneous terms of isba scheme
  LOGICAL :: LSURF_DIAG_ALBEDO   ! flag to write out diagnostic albedo
  LOGICAL :: LSURF_MISC_DIF      ! flag for miscellaneous terms of isba-dif scheme
  LOGICAL :: LPROSNOW            ! flag for Crocus-MEPRA outputs
!
!* variables for each patch
!
  REAL, POINTER, DIMENSION(:,:) :: XHV           ! Halstead coefficient
!      
  REAL, POINTER, DIMENSION(:,:,:) :: XSWI        ! Soil wetness index
  REAL, POINTER, DIMENSION(:,:,:) :: XTSWI       ! Total soil wetness index
!     
  REAL, POINTER, DIMENSION(:,:)   :: XALT        ! Active layer thickness in permafrost area
  REAL, POINTER, DIMENSION(:,:)   :: XFLT        ! Frozen layer thickness in non-permmafrost area
!
  REAL, POINTER, DIMENSION(:,:,:) :: XSNOWLIQ    ! snow liquid water profile (ISBA-ES:3-L)
  REAL, POINTER, DIMENSION(:,:,:) :: XSNOWTEMP   ! snow temperature profile  (ISBA-ES:3-L)

  REAL, POINTER, DIMENSION(:,:,:) :: XSNOWDZ     ! snow layer thickness (Crocus)
  REAL, POINTER, DIMENSION(:,:,:) :: XSNOWDEND   ! dendricity (Crocus)
  REAL, POINTER, DIMENSION(:,:,:) :: XSNOWSPHER  ! sphericity (Crocus)
  REAL, POINTER, DIMENSION(:,:,:) :: XSNOWSIZE    ! grain size (Crocus)
  REAL, POINTER, DIMENSION(:,:,:) :: XSNOWSSA    ! snow specific surface area (Crocus)
  REAL, POINTER, DIMENSION(:,:,:) :: XSNOWTYPEMEPRA ! snow grain type (Crocus-MEPRA)
  REAL, POINTER, DIMENSION(:,:,:) :: XSNOWRAM    ! snow ram resistance (Crocus-MEPRA)
  REAL, POINTER, DIMENSION(:,:,:) :: XSNOWSHEAR  ! snow shear resistance (Crocus-MEPRA)
!     
  REAL, POINTER, DIMENSION(:,:) :: XTWSNOW       ! Total snow reservoir
  REAL, POINTER, DIMENSION(:,:) :: XTDSNOW       ! Total snow height
  REAL, POINTER, DIMENSION(:,:) :: XTTSNOW       ! Total snow temperature
  
  REAL, POINTER, DIMENSION(:,:) :: XSNOWDEPTH_1DAYS,XSNOWDEPTH_3DAYS,XSNOWDEPTH_5DAYS,XSNOWDEPTH_7DAYS ! fresh snow depth in 1, 3, 5, 7 days
  REAL, POINTER, DIMENSION(:,:) :: XSNOWSWE_1DAYS,XSNOWSWE_3DAYS,XSNOWSWE_5DAYS,XSNOWSWE_7DAYS! fresh snow water equivalent in 1, 3, 5, 7 days
  REAL, POINTER, DIMENSION(:,:) :: XSNOWRAM_SONDE ! penetration depth of the ram resistance sensor (2 DaN)
  REAL, POINTER, DIMENSION(:,:) :: XSNOW_WETTHICKNESS ! Thickness of wet snow at the top of the snowpack
  REAL, POINTER, DIMENSION(:,:) :: XSNOW_REFROZENTHICKNESS  ! Thickness of refrozen snow at the top of the snowpack
!
  REAL, POINTER, DIMENSION(:,:) :: XSYTMASS    ! Eroded/accumulated snow mass (Crocus/SYTRON) (kg/m2/s)
  REAL, POINTER, DIMENSION(:,:) :: XSYTMASSC   ! Eroded/accumulated snow mass (Crocus/SYTRON) (kg/m2)
  
  REAL, POINTER, DIMENSION(:,:) :: XDPSNG         ! Snow fraction over ground, diag at time t
  REAL, POINTER, DIMENSION(:,:) :: XDPSNV         ! Snow fraction over vegetation, diag at time t
  REAL, POINTER, DIMENSION(:,:) :: XDPSN          ! Total Snow fraction, diag at time t
  REAL, POINTER, DIMENSION(:,:) :: XALBT          ! Total Albedo
!    
  REAL, POINTER, DIMENSION(:,:) :: XDFSAT         ! Topmodel/dt92 saturated fraction
!
  REAL, POINTER, DIMENSION(:,:) :: XDFFG          ! Flood fraction over ground, diag at time t
  REAL, POINTER, DIMENSION(:,:) :: XDFFV          ! Flood fraction over vegetation, diag at time t
  REAL, POINTER, DIMENSION(:,:) :: XDFF           ! Total Flood fraction, diag at time t
!
  REAL, POINTER, DIMENSION(:,:) :: XSEUIL        ! Irrigation threshold
!
  REAL, POINTER, DIMENSION(:,:) :: XFAPAR        ! Fapar of vegetation
  REAL, POINTER, DIMENSION(:,:) :: XFAPIR        ! Fapir of vegetation
  REAL, POINTER, DIMENSION(:,:) :: XDFAPARC      ! Fapar of vegetation (cumul)
  REAL, POINTER, DIMENSION(:,:) :: XDFAPIRC      ! Fapir of vegetation (cumul)
  REAL, POINTER, DIMENSION(:,:) :: XFAPAR_BS     ! Fapar of bare soil
  REAL, POINTER, DIMENSION(:,:) :: XFAPIR_BS     ! Fapir of bare soil
  REAL, POINTER, DIMENSION(:,:) :: XDLAI_EFFC    ! Effective LAI (cumul)
!  
!* average variables
!
  REAL, POINTER, DIMENSION(:)   :: XAVG_HV       ! Halstead coefficient
  REAL, POINTER, DIMENSION(:)   :: XAVG_LAI      ! leaf average index
!     
  REAL, POINTER, DIMENSION(:,:)   :: XAVG_SWI      ! Liquid Soil wetness index by layer
  REAL, POINTER, DIMENSION(:,:)   :: XAVG_TSWI     ! Total soil wetness index by layer
  REAL, POINTER, DIMENSION(:)     :: XSOIL_SWI     ! Soil wetness index
  REAL, POINTER, DIMENSION(:)     :: XSOIL_TSWI    ! Total Soil wetness index
  REAL, POINTER, DIMENSION(:)     :: XSOIL_TWG     ! Soil water content (liquid+ice) (kg.m-2)
  REAL, POINTER, DIMENSION(:)     :: XSOIL_TWGI    ! Soil ice content (kg.m-2)
  REAL, POINTER, DIMENSION(:)     :: XSOIL_WG     ! Soil water content (liquid+ice) (m3.m-3)
  REAL, POINTER, DIMENSION(:)     :: XSOIL_WGI    ! Soil ice content (m3.m-3)  
!     
  REAL, POINTER, DIMENSION(:)   :: XAVG_ALT      ! Active layer thickness in permafrost area
  REAL, POINTER, DIMENSION(:)   :: XAVG_FLT      ! Frozen layer thickness in non-permmafrost area
!
  REAL, POINTER, DIMENSION(:)   :: XAVG_TWSNOW   ! Snow total reservoir
  REAL, POINTER, DIMENSION(:)   :: XAVG_TDSNOW   ! Snow total height
  REAL, POINTER, DIMENSION(:)   :: XAVG_TTSNOW   ! Snow total temperature
!
  REAL, POINTER, DIMENSION(:)   :: XAVG_PSNG     ! Snow fraction over ground
  REAL, POINTER, DIMENSION(:)   :: XAVG_PSNV     ! Snow fraction over vegetation
  REAL, POINTER, DIMENSION(:)   :: XAVG_PSN      ! Total Snow fraction
  REAL, POINTER, DIMENSION(:)   :: XAVG_ALBT     ! Total Albedo
!    
  REAL, POINTER, DIMENSION(:) :: XAVG_FSAT       ! Topmodel/dt92 saturated fraction
!      
  REAL, POINTER, DIMENSION(:) :: XAVG_FFG        ! Flood fraction over ground
  REAL, POINTER, DIMENSION(:) :: XAVG_FFV        ! Flood fraction over vegetation
  REAL, POINTER, DIMENSION(:) :: XAVG_FF         ! Total Flood fraction  
!
  REAL, POINTER, DIMENSION(:) :: XAVG_SNOWDEPTH_1DAYS,XAVG_SNOWDEPTH_3DAYS,XAVG_SNOWDEPTH_5DAYS,XAVG_SNOWDEPTH_7DAYS ! fresh snow depth in 1, 3, 5, 7 days
  REAL, POINTER, DIMENSION(:) :: XAVG_SNOWSWE_1DAYS,XAVG_SNOWSWE_3DAYS,XAVG_SNOWSWE_5DAYS,XAVG_SNOWSWE_7DAYS! fresh snow water equivalent in 1, 3, 5, 7 days
  REAL, POINTER, DIMENSION(:) :: XAVG_SNOWRAM_SONDE ! penetration depth of the ram resistance sensor (2 DaN)
  REAL, POINTER, DIMENSION(:) :: XAVG_SNOW_WETTHICKNESS ! Thickness of wet snow at the top of the snowpack
  REAL, POINTER, DIMENSION(:) :: XAVG_SNOW_REFROZENTHICKNESS  ! Thickness of refrozen snow at the top of the snowpack
!  
  REAL, POINTER, DIMENSION(:) :: XAVG_SYTMASS    ! Eroded/accumulated snow mass (Crocus/SYTRON) (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XAVG_SYTMASSC   ! Eroded/accumulated snow mass (Crocus/SYTRON) (kg/m2)
!
  REAL, POINTER, DIMENSION(:) :: XFRD2_TSWI      ! ISBA-FR-DG2 comparable soil wetness index (DIF option)
  REAL, POINTER, DIMENSION(:) :: XFRD2_TWG       ! ISBA-FR-DG2 comparable soil water content (liquid+ice) (DIF option)
  REAL, POINTER, DIMENSION(:) :: XFRD2_TWGI      ! ISBA-FR-DG2 comparable soil ice content (DIF option)  
  REAL, POINTER, DIMENSION(:) :: XFRD3_TSWI      ! ISBA-FR-Deep comparable soil wetness index (DIF option)
  REAL, POINTER, DIMENSION(:) :: XFRD3_TWG       ! ISBA-FR-Deep comparable soil water content (liquid+ice) (DIF option)
  REAL, POINTER, DIMENSION(:) :: XFRD3_TWGI      ! ISBA-FR-Deep comparable soil ice content (DIF option)
!
!------------------------------------------------------------------------------
!

END TYPE DIAG_MISC_ISBA_t



CONTAINS

!
!


!


SUBROUTINE DIAG_MISC_ISBA_INIT(YDIAG_MISC_ISBA)
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: YDIAG_MISC_ISBA
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_MISC_ISBA_N:DIAG_MISC_ISBA_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YDIAG_MISC_ISBA%XHV)
  NULLIFY(YDIAG_MISC_ISBA%XSWI)
  NULLIFY(YDIAG_MISC_ISBA%XTSWI)
  NULLIFY(YDIAG_MISC_ISBA%XALT)
  NULLIFY(YDIAG_MISC_ISBA%XFLT)
  NULLIFY(YDIAG_MISC_ISBA%XSNOWLIQ)
  NULLIFY(YDIAG_MISC_ISBA%XSNOWTEMP)
  NULLIFY(YDIAG_MISC_ISBA%XSYTMASS)
  NULLIFY(YDIAG_MISC_ISBA%XSYTMASSC)
  NULLIFY(YDIAG_MISC_ISBA%XTWSNOW)
  NULLIFY(YDIAG_MISC_ISBA%XTDSNOW)
  NULLIFY(YDIAG_MISC_ISBA%XTTSNOW)
  NULLIFY(YDIAG_MISC_ISBA%XSNOWDZ)
  NULLIFY(YDIAG_MISC_ISBA%XSNOWDEND)
  NULLIFY(YDIAG_MISC_ISBA%XSNOWSPHER)
  NULLIFY(YDIAG_MISC_ISBA%XSNOWSIZE)
  NULLIFY(YDIAG_MISC_ISBA%XSNOWTYPEMEPRA)
  NULLIFY(YDIAG_MISC_ISBA%XSNOWSSA)
  NULLIFY(YDIAG_MISC_ISBA%XSNOWRAM)
  NULLIFY(YDIAG_MISC_ISBA%XSNOWSHEAR)
  NULLIFY(YDIAG_MISC_ISBA%XSNOWDEPTH_1DAYS)
  NULLIFY(YDIAG_MISC_ISBA%XSNOWDEPTH_3DAYS)
  NULLIFY(YDIAG_MISC_ISBA%XSNOWDEPTH_5DAYS)
  NULLIFY(YDIAG_MISC_ISBA%XSNOWDEPTH_7DAYS) 
  NULLIFY(YDIAG_MISC_ISBA%XSNOWSWE_1DAYS)
  NULLIFY(YDIAG_MISC_ISBA%XSNOWSWE_3DAYS)
  NULLIFY(YDIAG_MISC_ISBA%XSNOWSWE_5DAYS)
  NULLIFY(YDIAG_MISC_ISBA%XSNOWSWE_7DAYS)
  NULLIFY(YDIAG_MISC_ISBA%XSNOWRAM_SONDE)
  NULLIFY(YDIAG_MISC_ISBA%XSNOW_REFROZENTHICKNESS)
  NULLIFY(YDIAG_MISC_ISBA%XSNOW_WETTHICKNESS)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_SNOWDEPTH_1DAYS)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_SNOWDEPTH_3DAYS)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_SNOWDEPTH_5DAYS)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_SNOWDEPTH_7DAYS) 
  NULLIFY(YDIAG_MISC_ISBA%XAVG_SNOWSWE_1DAYS)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_SNOWSWE_3DAYS)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_SNOWSWE_5DAYS)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_SNOWSWE_7DAYS)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_SNOWRAM_SONDE)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_SNOW_REFROZENTHICKNESS)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_SNOW_WETTHICKNESS)  
  NULLIFY(YDIAG_MISC_ISBA%XDPSNG)
  NULLIFY(YDIAG_MISC_ISBA%XDPSNV)
  NULLIFY(YDIAG_MISC_ISBA%XDPSN)
  NULLIFY(YDIAG_MISC_ISBA%XALBT)
  NULLIFY(YDIAG_MISC_ISBA%XDFFG)
  NULLIFY(YDIAG_MISC_ISBA%XDFFV)
  NULLIFY(YDIAG_MISC_ISBA%XDFF)
  NULLIFY(YDIAG_MISC_ISBA%XSEUIL)
  NULLIFY(YDIAG_MISC_ISBA%XFAPAR)
  NULLIFY(YDIAG_MISC_ISBA%XFAPIR)
  NULLIFY(YDIAG_MISC_ISBA%XDFAPARC)
  NULLIFY(YDIAG_MISC_ISBA%XDFAPIRC)
  NULLIFY(YDIAG_MISC_ISBA%XFAPAR_BS)
  NULLIFY(YDIAG_MISC_ISBA%XFAPIR_BS)
  NULLIFY(YDIAG_MISC_ISBA%XDLAI_EFFC)  
  NULLIFY(YDIAG_MISC_ISBA%XAVG_HV)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_LAI)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_SWI)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_TSWI)
  NULLIFY(YDIAG_MISC_ISBA%XSOIL_SWI)
  NULLIFY(YDIAG_MISC_ISBA%XSOIL_TSWI)
  NULLIFY(YDIAG_MISC_ISBA%XSOIL_TWG)
  NULLIFY(YDIAG_MISC_ISBA%XSOIL_TWGI)
  NULLIFY(YDIAG_MISC_ISBA%XSOIL_WG)
  NULLIFY(YDIAG_MISC_ISBA%XSOIL_WGI)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_ALT)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_FLT)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_TWSNOW)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_TDSNOW)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_TTSNOW)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_PSNG)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_PSNV)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_PSN)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_ALBT)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_FFG)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_FFV)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_FF)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_SYTMASS)
  NULLIFY(YDIAG_MISC_ISBA%XAVG_SYTMASSC)  
  NULLIFY(YDIAG_MISC_ISBA%XFRD2_TSWI)
  NULLIFY(YDIAG_MISC_ISBA%XFRD2_TWG)
  NULLIFY(YDIAG_MISC_ISBA%XFRD2_TWGI)
  NULLIFY(YDIAG_MISC_ISBA%XFRD3_TSWI)
  NULLIFY(YDIAG_MISC_ISBA%XFRD3_TWG)
  NULLIFY(YDIAG_MISC_ISBA%XFRD3_TWGI)  
YDIAG_MISC_ISBA%LSURF_MISC_BUDGET=.FALSE.
YDIAG_MISC_ISBA%LSURF_DIAG_ALBEDO=.FALSE.
YDIAG_MISC_ISBA%LSURF_MISC_DIF=.FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_MISC_ISBA_N:DIAG_MISC_ISBA_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_MISC_ISBA_INIT


END MODULE MODD_DIAG_MISC_ISBA_n
