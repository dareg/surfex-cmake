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
  REAL, POINTER, DIMENSION(:) :: XHV       ! Halstead coefficient
  REAL, POINTER, DIMENSION(:) :: XLAI      ! leaf average index  
!      
  REAL, POINTER, DIMENSION(:,:) :: XSWI        ! Soil wetness index
  REAL, POINTER, DIMENSION(:,:) :: XTSWI       ! Total soil wetness index
  REAL, POINTER, DIMENSION(:)   :: XSOIL_SWI     ! Soil wetness index
  REAL, POINTER, DIMENSION(:)   :: XSOIL_TSWI    ! Total Soil wetness index
  REAL, POINTER, DIMENSION(:)   :: XSOIL_TWG     ! Soil water content (liquid+ice) (kg.m-2)
  REAL, POINTER, DIMENSION(:)   :: XSOIL_TWGI    ! Soil ice content (kg.m-2)
  REAL, POINTER, DIMENSION(:)   :: XSOIL_WG     ! Soil water content (liquid+ice) (m3.m-3)
  REAL, POINTER, DIMENSION(:)   :: XSOIL_WGI    ! Soil ice content (m3.m-3)  
!     
  REAL, POINTER, DIMENSION(:) :: XFRD2_TSWI      ! ISBA-FR-DG2 comparable soil wetness index (DIF option)
  REAL, POINTER, DIMENSION(:) :: XFRD2_TWG       ! ISBA-FR-DG2 comparable soil water content (liquid+ice) (DIF option)
  REAL, POINTER, DIMENSION(:) :: XFRD2_TWGI      ! ISBA-FR-DG2 comparable soil ice content (DIF option)  
  REAL, POINTER, DIMENSION(:) :: XFRD3_TSWI      ! ISBA-FR-Deep comparable soil wetness index (DIF option)
  REAL, POINTER, DIMENSION(:) :: XFRD3_TWG       ! ISBA-FR-Deep comparable soil water content (liquid+ice) (DIF option)
  REAL, POINTER, DIMENSION(:) :: XFRD3_TWGI      ! ISBA-FR-Deep comparable soil ice content (DIF option)
!
  REAL, POINTER, DIMENSION(:)   :: XALT        ! Active layer thickness in permafrost area
  REAL, POINTER, DIMENSION(:)   :: XFLT        ! Frozen layer thickness in non-permmafrost area
!
  REAL, POINTER, DIMENSION(:) :: XRNSNOW    ! net radiative flux from snow (ISBA-ES:3-L)    (W/m2)
  REAL, POINTER, DIMENSION(:) :: XHSNOW     ! sensible heat flux from snow (ISBA-ES:3-L)    (W/m2)
  REAL, POINTER, DIMENSION(:) :: XHPSNOW    ! heat release from rainfall (ISBA-ES:3-L)      (W/m2)
  REAL, POINTER, DIMENSION(:) :: XGFLUXSNOW ! net surface energy flux into snowpack      
!                                               ! (ISBA-ES:3-L)                                 (W/m2)
  REAL, POINTER, DIMENSION(:) :: XUSTARSNOW ! friction velocity  over snow 
!                                               ! (ISBA-ES:3-L)                                 (m/s)
  REAL, POINTER, DIMENSION(:) :: XGRNDFLUX  ! soil/snow interface heat flux (ISBA-ES:3-L)   (W/m2)
  REAL, POINTER, DIMENSION(:) :: XSRSFC     ! snowfall over snowpack (ISBA-ES:3-L)          (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XRRSFC     ! rainfall over snowpack (ISBA-ES:3-L)          (kg/m2/s)
  REAL, POINTER, DIMENSION(:) :: XCDSNOW    ! snow drag coefficient (ISBA-ES:3-L)           (-)
  REAL, POINTER, DIMENSION(:) :: XCHSNOW    ! heat turbulent transfer coefficient 
!                                               ! (ISBA-ES:3-L)                                 (-)
  REAL, POINTER, DIMENSION(:,:)::XSNOWDZ    ! snow layer thicknesses                        (m)
  REAL, POINTER, DIMENSION(:) :: XSNOWHMASS ! heat content change due to mass
!                                               ! changes in snowpack: for budget

!
  REAL, POINTER, DIMENSION(:,:) :: XSNOWLIQ    ! snow liquid water profile (ISBA-ES:3-L)
  REAL, POINTER, DIMENSION(:,:) :: XSNOWTEMP   ! snow temperature profile  (ISBA-ES:3-L)
!     
  REAL, POINTER, DIMENSION(:) :: XTWSNOW       ! Total snow reservoir
  REAL, POINTER, DIMENSION(:) :: XTDSNOW       ! Total snow height
  REAL, POINTER, DIMENSION(:) :: XTTSNOW       ! Total snow temperature
!
  REAL, POINTER, DIMENSION(:) :: XPSNG         ! Snow fraction over ground, diag at time t
  REAL, POINTER, DIMENSION(:) :: XPSNV         ! Snow fraction over vegetation, diag at time t
  REAL, POINTER, DIMENSION(:) :: XPSN          ! Total Snow fraction, diag at time t
!    
  REAL, POINTER, DIMENSION(:) :: XFSAT         ! Topmodel/dt92 saturated fraction
!
  REAL, POINTER, DIMENSION(:) :: XFFG          ! Flood fraction over ground, diag at time t
  REAL, POINTER, DIMENSION(:) :: XFFV          ! Flood fraction over vegetation, diag at time t
  REAL, POINTER, DIMENSION(:) :: XFF           ! Total Flood fraction, diag at time t
!
  REAL, POINTER, DIMENSION(:) :: XSEUIL        ! Irrigation threshold
!
  REAL, POINTER, DIMENSION(:) :: XFAPAR        ! Fapar of vegetation
  REAL, POINTER, DIMENSION(:) :: XFAPIR        ! Fapir of vegetation
  REAL, POINTER, DIMENSION(:) :: XDFAPARC      ! Fapar of vegetation (cumul)
  REAL, POINTER, DIMENSION(:) :: XDFAPIRC      ! Fapir of vegetation (cumul)
  REAL, POINTER, DIMENSION(:) :: XFAPAR_BS     ! Fapar of bare soil
  REAL, POINTER, DIMENSION(:) :: XFAPIR_BS     ! Fapir of bare soil
  REAL, POINTER, DIMENSION(:) :: XDLAI_EFFC    ! Effective LAI (cumul)
!
  REAL, POINTER, DIMENSION(:) :: XCG        ! heat capacity of the ground
  REAL, POINTER, DIMENSION(:) :: XC1        ! coefficients for the moisure
  REAL, POINTER, DIMENSION(:) :: XC2        ! equation.
  REAL, POINTER, DIMENSION(:) :: XWGEQ      ! equilibrium volumetric water
!                                               ! content
  REAL, POINTER, DIMENSION(:) :: XCT        ! area-averaged heat capacity
  REAL, POINTER, DIMENSION(:) :: XRS        ! stomatal resistance                            (s/m)
!
!------------------------------------------------------------------------------
!
END TYPE DIAG_MISC_ISBA_t
!
TYPE DIAG_MISC_ISBA_PATCH_t
!
TYPE(DIAG_MISC_ISBA_t), ALLOCATABLE :: AL(:) 
!
END TYPE DIAG_MISC_ISBA_PATCH_t
!
CONTAINS
!
SUBROUTINE DIAG_MISC_ISBA_PATCH_INIT(YDIAG_MISC_ISBA_PATCH,KPATCH)
TYPE(DIAG_MISC_ISBA_PATCH_t), INTENT(INOUT) :: YDIAG_MISC_ISBA_PATCH 
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_MISC_ISBA_N:DIAG_MISC_ISBA_PATCH_INIT",0,ZHOOK_HANDLE)
 ALLOCATE(YDIAG_MISC_ISBA_PATCH%AL(KPATCH))
DO JP=1,KPATCH
  CALL DIAG_MISC_ISBA_INIT(YDIAG_MISC_ISBA_PATCH%AL(JP))
ENDDO         
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_MISC_ISBA_N:DIAG_MISC_ISBA_PATCH_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_MISC_ISBA_PATCH_INIT
!
SUBROUTINE DIAG_MISC_ISBA_INIT(YDIAG_MISC_ISBA)
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: YDIAG_MISC_ISBA
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_MISC_ISBA_N:DIAG_MISC_ISBA_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YDIAG_MISC_ISBA%XHV)
  NULLIFY(YDIAG_MISC_ISBA%XLAI)  
  NULLIFY(YDIAG_MISC_ISBA%XSWI)
  NULLIFY(YDIAG_MISC_ISBA%XTSWI)
  NULLIFY(YDIAG_MISC_ISBA%XSOIL_SWI)
  NULLIFY(YDIAG_MISC_ISBA%XSOIL_TSWI)
  NULLIFY(YDIAG_MISC_ISBA%XSOIL_TWG)
  NULLIFY(YDIAG_MISC_ISBA%XSOIL_TWGI)
  NULLIFY(YDIAG_MISC_ISBA%XSOIL_WG)
  NULLIFY(YDIAG_MISC_ISBA%XSOIL_WGI)
  NULLIFY(YDIAG_MISC_ISBA%XFRD2_TWG)
  NULLIFY(YDIAG_MISC_ISBA%XFRD2_TWGI)
  NULLIFY(YDIAG_MISC_ISBA%XFRD3_TSWI)
  NULLIFY(YDIAG_MISC_ISBA%XFRD3_TWG)
  NULLIFY(YDIAG_MISC_ISBA%XFRD3_TWGI)    
  NULLIFY(YDIAG_MISC_ISBA%XALT)
  NULLIFY(YDIAG_MISC_ISBA%XFLT)
  NULLIFY(YDIAG_MISC_ISBA%XRNSNOW)
  NULLIFY(YDIAG_MISC_ISBA%XHSNOW)
  NULLIFY(YDIAG_MISC_ISBA%XHPSNOW)
  NULLIFY(YDIAG_MISC_ISBA%XGFLUXSNOW)
  NULLIFY(YDIAG_MISC_ISBA%XUSTARSNOW)
  NULLIFY(YDIAG_MISC_ISBA%XGRNDFLUX)
  NULLIFY(YDIAG_MISC_ISBA%XSRSFC)
  NULLIFY(YDIAG_MISC_ISBA%XRRSFC)
  NULLIFY(YDIAG_MISC_ISBA%XCDSNOW)
  NULLIFY(YDIAG_MISC_ISBA%XCHSNOW)
  NULLIFY(YDIAG_MISC_ISBA%XSNOWDZ)
  NULLIFY(YDIAG_MISC_ISBA%XSNOWHMASS)  
  NULLIFY(YDIAG_MISC_ISBA%XSNOWLIQ)
  NULLIFY(YDIAG_MISC_ISBA%XSNOWTEMP)
  NULLIFY(YDIAG_MISC_ISBA%XTWSNOW)
  NULLIFY(YDIAG_MISC_ISBA%XTDSNOW)
  NULLIFY(YDIAG_MISC_ISBA%XTTSNOW)
  NULLIFY(YDIAG_MISC_ISBA%XPSNG)
  NULLIFY(YDIAG_MISC_ISBA%XPSNV)
  NULLIFY(YDIAG_MISC_ISBA%XPSN)
  NULLIFY(YDIAG_MISC_ISBA%XFFG)
  NULLIFY(YDIAG_MISC_ISBA%XFFV)
  NULLIFY(YDIAG_MISC_ISBA%XFF)
  NULLIFY(YDIAG_MISC_ISBA%XSEUIL)
  NULLIFY(YDIAG_MISC_ISBA%XFAPAR)
  NULLIFY(YDIAG_MISC_ISBA%XFAPIR)
  NULLIFY(YDIAG_MISC_ISBA%XDFAPARC)
  NULLIFY(YDIAG_MISC_ISBA%XDFAPIRC)
  NULLIFY(YDIAG_MISC_ISBA%XFAPAR_BS)
  NULLIFY(YDIAG_MISC_ISBA%XFAPIR_BS)
  NULLIFY(YDIAG_MISC_ISBA%XDLAI_EFFC) 
  NULLIFY(YDIAG_MISC_ISBA%XCG)
  NULLIFY(YDIAG_MISC_ISBA%XC1)
  NULLIFY(YDIAG_MISC_ISBA%XC2)
  NULLIFY(YDIAG_MISC_ISBA%XWGEQ)
  NULLIFY(YDIAG_MISC_ISBA%XCT)
  NULLIFY(YDIAG_MISC_ISBA%XRS)  
YDIAG_MISC_ISBA%LSURF_MISC_BUDGET=.FALSE.
YDIAG_MISC_ISBA%LSURF_DIAG_ALBEDO=.FALSE.
YDIAG_MISC_ISBA%LSURF_MISC_DIF=.FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_MISC_ISBA_N:DIAG_MISC_ISBA_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_MISC_ISBA_INIT


END MODULE MODD_DIAG_MISC_ISBA_n
