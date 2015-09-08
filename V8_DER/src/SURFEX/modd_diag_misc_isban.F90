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
!!
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
!     
  REAL, POINTER, DIMENSION(:,:) :: XTWSNOW       ! Total snow reservoir
  REAL, POINTER, DIMENSION(:,:) :: XTDSNOW       ! Total snow height
  REAL, POINTER, DIMENSION(:,:) :: XTTSNOW       ! Total snow temperature
!
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

TYPE(DIAG_MISC_ISBA_t), ALLOCATABLE, TARGET, SAVE :: DIAG_MISC_ISBA_MODEL(:)

TYPE(DIAG_MISC_ISBA_t), POINTER :: DIAG_MISC_ISBA => NULL()
!$OMP THREADPRIVATE(DIAG_MISC_ISBA)

CONTAINS

SUBROUTINE DIAG_MISC_ISBA_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_DIAG_MISC_ISBA_N:DIAG_MISC_ISBA_GOTO_MODEL',0,ZHOOK_HANDLE)

DIAG_MISC_ISBA => DIAG_MISC_ISBA_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_DIAG_MISC_ISBA_N:DIAG_MISC_ISBA_GOTO_MODEL',1,ZHOOK_HANDLE)
!

END SUBROUTINE DIAG_MISC_ISBA_GOTO_MODEL

SUBROUTINE DIAG_MISC_ISBA_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_MISC_ISBA_N:DIAG_MISC_ISBA_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(DIAG_MISC_ISBA_MODEL(KMODEL))
DIAG_MISC_ISBA => DIAG_MISC_ISBA_MODEL(KMODEL)
DO J=1,KMODEL
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XHV)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XSWI)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XTSWI)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XALT)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XFLT)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XSNOWLIQ)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XSNOWTEMP)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XTWSNOW)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XTDSNOW)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XTTSNOW)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XDPSNG)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XDPSNV)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XDPSN)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XALBT)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XDFFG)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XDFFV)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XDFF)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XSEUIL)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XFAPAR)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XFAPIR)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XDFAPARC)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XDFAPIRC)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XFAPAR_BS)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XFAPIR_BS)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XDLAI_EFFC)  
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XAVG_HV)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XAVG_LAI)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XAVG_SWI)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XAVG_TSWI)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XSOIL_SWI)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XSOIL_TSWI)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XSOIL_TWG)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XSOIL_TWGI)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XSOIL_WG)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XSOIL_WGI)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XAVG_ALT)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XAVG_FLT)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XAVG_TWSNOW)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XAVG_TDSNOW)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XAVG_TTSNOW)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XAVG_PSNG)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XAVG_PSNV)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XAVG_PSN)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XAVG_ALBT)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XAVG_FFG)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XAVG_FFV)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XAVG_FF)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XFRD2_TSWI)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XFRD2_TWG)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XFRD2_TWGI)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XFRD3_TSWI)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XFRD3_TWG)
  NULLIFY(DIAG_MISC_ISBA_MODEL(J)%XFRD3_TWGI)  
ENDDO
DIAG_MISC_ISBA_MODEL(:)%LSURF_MISC_BUDGET=.FALSE.
DIAG_MISC_ISBA_MODEL(:)%LSURF_DIAG_ALBEDO=.FALSE.
DIAG_MISC_ISBA_MODEL(:)%LSURF_MISC_DIF=.FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_MISC_ISBA_N:DIAG_MISC_ISBA_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_MISC_ISBA_ALLOC

SUBROUTINE DIAG_MISC_ISBA_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_MISC_ISBA_N:DIAG_MISC_ISBA_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(DIAG_MISC_ISBA_MODEL)) DEALLOCATE(DIAG_MISC_ISBA_MODEL)
IF (ASSOCIATED(DIAG_MISC_ISBA)) NULLIFY(DIAG_MISC_ISBA)
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_MISC_ISBA_N:DIAG_MISC_ISBA_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_MISC_ISBA_DEALLO

END MODULE MODD_DIAG_MISC_ISBA_n
