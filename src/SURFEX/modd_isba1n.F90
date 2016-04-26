!##################
MODULE MODD_ISBA1_n
!##################
!
!!****  *MODD_ISBA - declaration of packed surface parameters for ISBA scheme
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
!!      A. Boone   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       20/09/02
!!      A.L. Gibelin    04/2009 : BIOMASS and RESP_BIOMASS arrays 
!!      A.L. Gibelin    04/2009 : TAU_WOOD for NCB option 
!!      A.L. Gibelin    05/2009 : Add carbon spinup
!!      A.L. Gibelin    06/2009 : Soil carbon variables for CNT option
!!      A.L. Gibelin    07/2009 : Suppress RDK and transform GPP as a diagnostic
!!      A.L. Gibelin    07/2009 : Suppress PPST and PPSTF as outputs
!!      P. Samuelsson   02/2012 : MEB
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_ISBA_PGD_n, ONLY : ISBA_PGD_NP_t, ISBA_PGD_NP_INIT
USE MODD_ISBA_PARAM1_n
USE MODD_ISBA_INIT1_n
USE MODD_GRID_n, ONLY : GRID_NP_t, GRID_NP_INIT
USE MODD_AGRI1_n
USE MODD_SSO1_n
!
USE MODD_TYPE_SNOW1
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
TYPE ISBA_PROG_1P_t
!
! Prognostic variables:
!
! - Soil and vegetation heat and water:
!
REAL, POINTER, DIMENSION(:,:) :: XWG           ! soil volumetric water content profile   (m3/m3)
REAL, POINTER, DIMENSION(:,:) :: XWGI          ! soil liquid water equivalent volumetric 
!                                                ! ice content profile                     (m3/m3)
REAL, POINTER, DIMENSION(:)   :: XWR           ! liquid water retained on the
!                                                ! foliage of the vegetation
!                                                ! canopy                                  (kg/m2)
REAL, POINTER, DIMENSION(:,:) :: XTG           ! surface and sub-surface soil 
!                                                ! temperature profile                     (K)
!
! - Snow Cover:
!
TYPE(SURF_SNOW1) :: TSNOW                         ! snow state: 
!                                                ! scheme type/option                      (-)
!                                                ! number of layers                        (-)
!                                                ! snow (& liq. water) content             (kg/m2)
!                                                ! heat content                            (J/m2)
!                                                ! temperature                             (K)
!                                                ! density                                 (kg m-3)
!
REAL, POINTER, DIMENSION(:) :: XICE_STO        ! Glacier ice storage reservoir
!
! - For multi-energy balance:
!
REAL, POINTER, DIMENSION(:) :: XWRL            ! liquid water retained on litter          (kg/m2)
REAL, POINTER, DIMENSION(:) :: XWRLI           ! ice retained on litter          (kg/m2)
REAL, POINTER, DIMENSION(:) :: XWRVN           ! snow retained on the foliage
!                                                ! of the canopy vegetation                  (kg/m2)
REAL, POINTER, DIMENSION(:) :: XTV             ! canopy vegetation temperature             (K)
REAL, POINTER, DIMENSION(:) :: XTL             ! litter temperature             (K)
REAL, POINTER, DIMENSION(:) :: XTC             ! canopy air temperature                    (K)
REAL, POINTER, DIMENSION(:) :: XQC             ! canopy air specific humidity              (kg/kg)
!
! * Half prognostic fields
!
REAL, POINTER, DIMENSION(:)     :: XRESA         ! aerodynamic resistance                  (s/m)
!
! - Vegetation: Ags Prognostic (YPHOTO = 'AGS', 'LAI', 'AST', 'LST', 'NIT', 'NCB')
!
REAL, POINTER, DIMENSION(:) :: XAN           ! net CO2 assimilation                    (mg/m2/s)
REAL, POINTER, DIMENSION(:) :: XANDAY        ! daily net CO2 assimilation              (mg/m2)
REAL, POINTER, DIMENSION(:) :: XANFM         ! maximum leaf assimilation               (mg/m2/s)
REAL, POINTER, DIMENSION(:) :: XLE           ! evapotranspiration                      (W/m2)
!
REAL, POINTER, DIMENSION(:) :: XFAPARC       ! Fapar of vegetation (cumul)
REAL, POINTER, DIMENSION(:) :: XFAPIRC       ! Fapir of vegetation (cumul)
REAL, POINTER, DIMENSION(:) :: XLAI_EFFC     ! Effective LAI (cumul)
REAL, POINTER, DIMENSION(:) :: XMUS          ! cos zenithal angle (cumul)
!
REAL, POINTER, DIMENSION(:,:) :: XRESP_BIOMASS    ! daily cumulated respiration of 
!                                                   ! biomass                              (kg/m2/s)
REAL, POINTER, DIMENSION(:,:) :: XBIOMASS         ! biomass of previous day              (kg/m2) 
!
! - Soil carbon (ISBA-CC, YRESPSL = 'CNT')
!
REAL, POINTER, DIMENSION(:,:,:) :: XLITTER          ! litter pools                         (gC/m2)
REAL, POINTER, DIMENSION(:,:)   :: XSOILCARB        ! soil carbon pools                    (gC/m2) 
REAL, POINTER, DIMENSION(:,:)   :: XLIGNIN_STRUC    ! ratio Lignin/Carbon in structural
!                                                       litter                               (gC/m2)
!
REAL, POINTER, DIMENSION(:) :: XPSNG         ! Snow fraction over ground
REAL, POINTER, DIMENSION(:) :: XPSNV         ! Snow fraction over vegetation
REAL, POINTER, DIMENSION(:) :: XPSNV_A       ! Snow fraction over vegetation
REAL, POINTER, DIMENSION(:) :: XPSN
!
REAL, POINTER, DIMENSION(:)   :: XSNOWFREE_ALB     ! snow free albedo                        (-)
REAL, POINTER, DIMENSION(:)   :: XSNOWFREE_ALB_VEG ! snow free albedo for vegetation         (-)
REAL, POINTER, DIMENSION(:)   :: XSNOWFREE_ALB_SOIL! snow free albedo for soil
!
END TYPE ISBA_PROG_1P_t
!
!
TYPE ISBA_PROG_NP_t
!
TYPE(ISBA_PROG_1P_t), ALLOCATABLE :: AL(:)
!
END TYPE ISBA_PROG_NP_t
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
TYPE ISBA_NP_t
!
TYPE(ISBA_PGD_NP_t) :: P
TYPE(ISBA_PARAM_NP_t) :: M
TYPE(ISBA_INIT_NP_t) :: I
TYPE(ISBA_INIT_PGD_NP_t) :: IP
TYPE(ISBA_PROG_NP_t) :: R
!
TYPE(GRID_NP_t) :: G
TYPE(AGRI_NP_t) :: AG
TYPE(SSO_NP_t) :: ISS
!
END TYPE ISBA_NP_t
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
CONTAINS
!
SUBROUTINE ISBA_PROG_1P_INIT(R)
TYPE(ISBA_PROG_1P_t), INTENT(INOUT) :: R
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA1_N:ISBA_PROG_1P_INIT",0,ZHOOK_HANDLE)
!
  NULLIFY(R%XWG)
  NULLIFY(R%XWGI)
  NULLIFY(R%XWR)
  NULLIFY(R%XTG)
  NULLIFY(R%XICE_STO)
  NULLIFY(R%XWRL)
  NULLIFY(R%XWRLI)
  NULLIFY(R%XWRVN)
  NULLIFY(R%XTV)
  NULLIFY(R%XTL)
  NULLIFY(R%XTC)
  NULLIFY(R%XQC)
  NULLIFY(R%XRESA)
  NULLIFY(R%XAN)
  NULLIFY(R%XANDAY)
  NULLIFY(R%XANFM)
  NULLIFY(R%XLE)
  NULLIFY(R%XFAPARC)
  NULLIFY(R%XFAPIRC)
  NULLIFY(R%XLAI_EFFC)  
  NULLIFY(R%XMUS)   
  NULLIFY(R%XRESP_BIOMASS)
  NULLIFY(R%XBIOMASS)
  NULLIFY(R%XLITTER)
  NULLIFY(R%XSOILCARB)
  NULLIFY(R%XLIGNIN_STRUC)
  NULLIFY(R%XPSNG)
  NULLIFY(R%XPSNV)
  NULLIFY(R%XPSNV_A)
  NULLIFY(R%XPSN)
  NULLIFY(R%XSNOWFREE_ALB)
  NULLIFY(R%XSNOWFREE_ALB_VEG)
  NULLIFY(R%XSNOWFREE_ALB_SOIL)  
!
IF (LHOOK) CALL DR_HOOK("MODD_ISBA1_N:ISBA_PROG_1P_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_PROG_1P_INIT
!
SUBROUTINE ISBA_PROG_NP_INIT(R,KPATCH)
TYPE(ISBA_PROG_NP_t), INTENT(INOUT) :: R
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA1_N:ISBA_PROG_NP_INIT",0,ZHOOK_HANDLE)
!
ALLOCATE(R%AL(KPATCH))
DO JP=1,KPATCH
  CALL ISBA_PROG_1P_INIT(R%AL(JP))
ENDDO
!
IF (LHOOK) CALL DR_HOOK("MODD_ISBA1_N:ISBA_PROG_NP_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_PROG_NP_INIT

 !
SUBROUTINE ISBA_NP_INIT(IK,KPATCH)
TYPE(ISBA_NP_t), INTENT(INOUT) :: IK
INTEGER, INTENT(IN) :: KPATCH
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA1_N:ISBA_NP_INIT",0,ZHOOK_HANDLE)
!
CALL ISBA_PROG_NP_INIT(IK%R,KPATCH)
CALL ISBA_INIT_NP_INIT(IK%I,KPATCH)
CALL ISBA_INIT_PGD_NP_INIT(IK%IP,KPATCH)
CALL ISBA_PARAM_NP_INIT(IK%M%X,IK%M%M,IK%M%A,IK%M%I,KPATCH)
CALL ISBA_PARAM_TIME_NP_INIT(IK%M%T,KPATCH)
CALL ISBA_PGD_NP_INIT(IK%P,KPATCH)
!
CALL AGRI_NP_INIT(IK%AG,KPATCH)
CALL GRID_NP_INIT(IK%G,KPATCH)
CALL SSO_NP_INIT(IK%ISS,KPATCH)
!
IF (LHOOK) CALL DR_HOOK("MODD_ISBA1_N:ISBA_NP_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_NP_INIT
!

END MODULE MODD_ISBA1_n
