!##################
MODULE MODD_ISBA_n
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
USE MODD_ISBA_OPTIONS_n
USE MODD_ISBA_PGD_n
USE MODD_ISBA_PARAM_n
USE MODD_ISBA_INIT_n
!
USE MODD_TEB_VEG_PARAM_n
!
USE MODD_TYPE_SNOW
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
TYPE ISBA_PROG_t
!
! Prognostic variables:
!
! - Soil and vegetation heat and water:
!
REAL, POINTER, DIMENSION(:,:,:) :: XWG           ! soil volumetric water content profile   (m3/m3)
REAL, POINTER, DIMENSION(:,:,:) :: XWGI          ! soil liquid water equivalent volumetric 
!                                                ! ice content profile                     (m3/m3)
REAL, POINTER, DIMENSION(:,:)   :: XWR           ! liquid water retained on the
!                                                ! foliage of the vegetation
!                                                ! canopy                                  (kg/m2)
REAL, POINTER, DIMENSION(:,:,:) :: XTG           ! surface and sub-surface soil 
!                                                ! temperature profile                     (K)
!
! - Snow Cover:
!
TYPE(SURF_SNOW) :: TSNOW                         ! snow state: 
!                                                ! scheme type/option                      (-)
!                                                ! number of layers                        (-)
!                                                ! snow (& liq. water) content             (kg/m2)
!                                                ! heat content                            (J/m2)
!                                                ! temperature                             (K)
!                                                ! density                                 (kg m-3)
!
REAL, POINTER, DIMENSION(:,:) :: XICE_STO        ! Glacier ice storage reservoir
!
! - For multi-energy balance:
!
REAL, POINTER, DIMENSION(:,:) :: XWRL            ! liquid water retained on litter          (kg/m2)
REAL, POINTER, DIMENSION(:,:) :: XWRLI           ! ice retained on litter          (kg/m2)
REAL, POINTER, DIMENSION(:,:) :: XWRVN           ! snow retained on the foliage
!                                                ! of the canopy vegetation                  (kg/m2)
REAL, POINTER, DIMENSION(:,:) :: XTV             ! canopy vegetation temperature             (K)
REAL, POINTER, DIMENSION(:,:) :: XTL             ! litter temperature             (K)
REAL, POINTER, DIMENSION(:,:) :: XTC             ! canopy air temperature                    (K)
REAL, POINTER, DIMENSION(:,:) :: XQC             ! canopy air specific humidity              (kg/kg)
!
! * Half prognostic fields
!
REAL, POINTER, DIMENSION(:,:)     :: XRESA         ! aerodynamic resistance                  (s/m)
!
REAL, POINTER, DIMENSION(:)   :: XTSRAD_NAT        ! patch averaged radiative temperature    (K)
!
! - Vegetation: Ags Prognostic (YPHOTO = 'AGS', 'LAI', 'AST', 'LST', 'NIT', 'NCB')
!
REAL, POINTER, DIMENSION(:,:) :: XAN           ! net CO2 assimilation                    (mg/m2/s)
REAL, POINTER, DIMENSION(:,:) :: XANDAY        ! daily net CO2 assimilation              (mg/m2)
REAL, POINTER, DIMENSION(:,:) :: XANFM         ! maximum leaf assimilation               (mg/m2/s)
REAL, POINTER, DIMENSION(:,:) :: XLE           ! evapotranspiration                      (W/m2)
!
REAL, POINTER, DIMENSION(:,:) :: XFAPARC       ! Fapar of vegetation (cumul)
REAL, POINTER, DIMENSION(:,:) :: XFAPIRC       ! Fapir of vegetation (cumul)
REAL, POINTER, DIMENSION(:,:) :: XLAI_EFFC     ! Effective LAI (cumul)
REAL, POINTER, DIMENSION(:,:) :: XMUS          ! cos zenithal angle (cumul)
!
REAL, POINTER, DIMENSION(:,:,:) :: XRESP_BIOMASS    ! daily cumulated respiration of 
!                                                   ! biomass                              (kg/m2/s)
REAL, POINTER, DIMENSION(:,:,:) :: XBIOMASS         ! biomass of previous day              (kg/m2) 
!
! - Soil carbon (ISBA-CC, YRESPSL = 'CNT')
!
REAL, POINTER, DIMENSION(:,:,:,:) :: XLITTER          ! litter pools                         (gC/m2)
REAL, POINTER, DIMENSION(:,:,:)   :: XSOILCARB        ! soil carbon pools                    (gC/m2) 
REAL, POINTER, DIMENSION(:,:,:)   :: XLIGNIN_STRUC    ! ratio Lignin/Carbon in structural
!                                                       litter                               (gC/m2)
!
REAL, POINTER, DIMENSION(:,:) :: XPSNG         ! Snow fraction over ground
REAL, POINTER, DIMENSION(:,:) :: XPSNV         ! Snow fraction over vegetation
REAL, POINTER, DIMENSION(:,:) :: XPSNV_A       ! Snow fraction over vegetation
REAL, POINTER, DIMENSION(:,:) :: XPSN
!
REAL, POINTER, DIMENSION(:,:)   :: XSNOWFREE_ALB     ! snow free albedo                        (-)
REAL, POINTER, DIMENSION(:,:)   :: XSNOWFREE_ALB_VEG ! snow free albedo for vegetation         (-)
REAL, POINTER, DIMENSION(:,:)   :: XSNOWFREE_ALB_SOIL! snow free albedo for soil
!
!  - Assimilation: ENKF
!
REAL, POINTER, DIMENSION(:,:,:)   :: XRED_NOISE
REAL, POINTER, DIMENSION(:,:)     :: XINCR
!
END TYPE ISBA_PROG_t
!
!
TYPE ISBA_PROG_PATCH_t
!
TYPE(ISBA_PROG_t), ALLOCATABLE :: AL(:)
!
END TYPE ISBA_PROG_PATCH_t
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
TYPE ISBA_t
!
TYPE(ISBA_OPTIONS_t) :: O
TYPE(ISBA_PGD_t) :: P
TYPE(ISBA_PARAM_t) :: M
TYPE(ISBA_INIT_t) :: I
TYPE(ISBA_INIT_PGD_t) :: IP
TYPE(ISBA_PROG_t) :: R
!
END TYPE ISBA_t
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
CONTAINS
!
SUBROUTINE ISBA_PROG_INIT(YISBA_PROG)
TYPE(ISBA_PROG_t), INTENT(INOUT) :: YISBA_PROG
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_N:ISBA_PROG_INIT",0,ZHOOK_HANDLE)
!
NULLIFY(YISBA_PROG%XWG)
NULLIFY(YISBA_PROG%XWGI)
NULLIFY(YISBA_PROG%XWR)
NULLIFY(YISBA_PROG%XTG)
NULLIFY(YISBA_PROG%XICE_STO)
NULLIFY(YISBA_PROG%XWRL)
NULLIFY(YISBA_PROG%XWRLI)
NULLIFY(YISBA_PROG%XWRVN)
NULLIFY(YISBA_PROG%XTV)
NULLIFY(YISBA_PROG%XTL)
NULLIFY(YISBA_PROG%XTC)
NULLIFY(YISBA_PROG%XQC)
NULLIFY(YISBA_PROG%XRESA)
NULLIFY(YISBA_PROG%XTSRAD_NAT)
NULLIFY(YISBA_PROG%XAN)
NULLIFY(YISBA_PROG%XANDAY)
NULLIFY(YISBA_PROG%XANFM)
NULLIFY(YISBA_PROG%XLE)
NULLIFY(YISBA_PROG%XFAPARC)
NULLIFY(YISBA_PROG%XFAPIRC)
NULLIFY(YISBA_PROG%XLAI_EFFC)  
NULLIFY(YISBA_PROG%XMUS)   
NULLIFY(YISBA_PROG%XRESP_BIOMASS)
NULLIFY(YISBA_PROG%XBIOMASS)
NULLIFY(YISBA_PROG%XLITTER)
NULLIFY(YISBA_PROG%XSOILCARB)
NULLIFY(YISBA_PROG%XLIGNIN_STRUC)
NULLIFY(YISBA_PROG%XPSNG)
NULLIFY(YISBA_PROG%XPSNV)
NULLIFY(YISBA_PROG%XPSNV_A)
NULLIFY(YISBA_PROG%XSNOWFREE_ALB)
NULLIFY(YISBA_PROG%XSNOWFREE_ALB_VEG)
NULLIFY(YISBA_PROG%XSNOWFREE_ALB_SOIL)
NULLIFY(YISBA_PROG%XPSN)
NULLIFY(YISBA_PROG%XRED_NOISE)
NULLIFY(YISBA_PROG%XINCR)
!
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_N:ISBA_PROG_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_PROG_INIT

 !
SUBROUTINE ISBA_INIT(YISBA)
TYPE(ISBA_t), INTENT(INOUT) :: YISBA
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_N:ISBA_INIT",0,ZHOOK_HANDLE)
!
CALL ISBA_PROG_INIT(YISBA%R)
CALL ISBA_INIT_INIT(YISBA%I)
CALL ISBA_INIT_PGD_INIT(YISBA%IP)
CALL ISBA_PARAM_INIT(YISBA%M%X,YISBA%M%M,YISBA%M%A,YISBA%M%I)
CALL ISBA_PARAM_TIME_INIT(YISBA%M%T)
CALL ISBA_PGD_INIT(YISBA%P)
CALL ISBA_OPTIONS_INIT(YISBA%O)
!
IF (LHOOK) CALL DR_HOOK("MODD_ISBA_N:ISBA_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE ISBA_INIT
!

END MODULE MODD_ISBA_n
