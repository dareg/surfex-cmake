!     #########
      SUBROUTINE AVG_ALBEDO_EMIS_GARDEN (GDR, HALBEDO, GDMT, GDMA, &
                                 PTG1, PSW_BANDS, PDIR_ALB,PSCA_ALB, PEMIS, PTSRAD         )  
!     ###################################################
!
!!**** ** computes radiative fields used in GARDEN
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    01/2004
!!     A. Bogatchev 09/2005 EBA snow option
!!     B. Decharme  2008    The fraction of vegetation covered by snow must be
!                            <= to XPSNG
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
!
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_TIME_t, ISBA_PARAM_ALB_t
!
USE MODD_SURF_PAR,  ONLY : XUNDEF
!
USE MODD_TYPE_SNOW
!
USE MODD_SNOW_PAR,   ONLY : XEMISSN
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
!
USE MODI_ALBEDO
USE MODI_ALBEDO_FROM_NIR_VIS
USE MODI_ISBA_SNOW_FRAC
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
!
TYPE(ISBA_PROG_t), INTENT(INOUT) :: GDR
!
 CHARACTER(LEN=4),       INTENT(IN)   :: HALBEDO     ! albedo type
! Albedo dependance with surface soil water content
!   "EVOL" = albedo evolves with soil wetness
!   "DRY " = constant albedo value for dry soil
!   "WET " = constant albedo value for wet soil
!   "MEAN" = constant albedo value for medium soil wetness
!
TYPE(ISBA_PARAM_TIME_t), INTENT(INOUT) :: GDMT
TYPE(ISBA_PARAM_ALB_t), INTENT(INOUT) :: GDMA
!
REAL, DIMENSION(:),   INTENT(IN)   :: PTG1        ! soil surface temperature
REAL, DIMENSION(:),     INTENT(IN)   :: PSW_BANDS   ! middle wavelength of each band 
!
REAL, DIMENSION(:,:),   INTENT(OUT)  :: PDIR_ALB    ! averaged direct albedo  (per wavelength)
REAL, DIMENSION(:,:),   INTENT(OUT)  :: PSCA_ALB    ! averaged diffuse albedo (per wavelength)
REAL, DIMENSION(:),     INTENT(OUT)  :: PEMIS       ! averaged emissivity
REAL, DIMENSION(:),     INTENT(OUT)  :: PTSRAD      ! averaged radiaitve temp.
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
!
REAL, DIMENSION(SIZE(GDMT%XALBNIR_VEG(:,1))) :: ZALBNIR ! near-infra-red albedo with snow
REAL, DIMENSION(SIZE(GDMT%XALBVIS_VEG(:,1))) :: ZALBVIS ! visible albedo with snow
REAL, DIMENSION(SIZE(GDMT%XALBUV_VEG(:,1) )) :: ZALBUV  ! UV albedo with snow
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!
!*    1.      averaged albedo on natural continental surfaces (except prognostic snow)
!             -----------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('AVG_ALBEDO_EMIS_GARDEN',0,ZHOOK_HANDLE)
 CALL ALBEDO(HALBEDO, GDMT, GDMA )  

!
!*    2.      averaged albedo and emis. on natural continental surfaces (with prognostic snow)
!             ---------------------------------------------------------
!
ZALBNIR(:)=0.
ZALBVIS(:)=0.
ZALBUV (:)=0.
!
PDIR_ALB(:,:)=0.
PSCA_ALB(:,:)=0.
PEMIS   (:)  =0.
PTSRAD  (:)  =0.
!   
!
  CALL ISBA_SNOW_FRAC(GDR%TSNOW%SCHEME,           &
         GDR%TSNOW%WSNOW(:,:,1), GDR%TSNOW%RHO(:,:,1), GDR%TSNOW%ALB  (:,1),  &
         GDMT%XVEG(:,1), GDMT%XLAI(:,1), GDMT%XZ0(:,1),            &
         GDR%XPSN(:,1), GDR%XPSNV_A(:,1), GDR%XPSNG(:,1), GDR%XPSNV(:,1)      )
!
 WHERE (GDMT%XVEG(:,1)/=XUNDEF)
!
! albedo on this tile
!
    ZALBNIR(:) = (1.-GDR%XPSN(:,1))*GDMT%XALBNIR(:,1) + GDR%XPSN(:,1) *GDR%TSNOW%ALB (:,1)  
      
    ZALBVIS(:) = (1.-GDR%XPSN(:,1))*GDMT%XALBVIS(:,1) + GDR%XPSN(:,1) *GDR%TSNOW%ALB (:,1)  
      
    ZALBUV(:)  = (1.-GDR%XPSN(:,1))*GDMT%XALBUV(:,1) + GDR%XPSN(:,1) *GDR%TSNOW%ALB (:,1)  
  END WHERE
!
!* albedo for each wavelength
!
  CALL ALBEDO_FROM_NIR_VIS(PSW_BANDS,ZALBNIR, ZALBVIS, ZALBUV,  &
                           PDIR_ALB(:,:), PSCA_ALB(:,:)         )  
!
! emissivity
!
  WHERE (GDMT%XEMIS(:,1)/=XUNDEF)
    PEMIS(:)   = (1.-GDR%XPSN(:,1))*GDMT%XEMIS(:,1) + GDR%XPSN(:,1) *XEMISSN  
  END WHERE
!
!* radiative surface temperature
!
  IF (GDR%TSNOW%SCHEME=='D95' .OR. GDR%TSNOW%SCHEME=='EBA') THEN
    PTSRAD(:) = PTG1(:)
  ELSE IF (GDR%TSNOW%SCHEME=='3-L' .OR. GDR%TSNOW%SCHEME=='CRO') THEN
    WHERE (GDMT%XEMIS(:,1)/=XUNDEF)
    PTSRAD(:) =( ( (1.-GDR%XPSN(:,1))*PEMIS      (:)       *PTG1     (:)**4         &
                  +    GDR%XPSN(:,1) *GDR%TSNOW%EMIS(:,1)*GDR%TSNOW%TS(:,1)**4 ) )**0.25  &
                             / PEMIS(:)**0.25  
    END WHERE
  END IF
!
IF (LHOOK) CALL DR_HOOK('AVG_ALBEDO_EMIS_GARDEN',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE AVG_ALBEDO_EMIS_GARDEN
