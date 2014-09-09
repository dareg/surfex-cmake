!     #########
      SUBROUTINE AVERAGED_ALBEDO_EMIS_ISBA(OFLOOD, HALBEDO, &
                                 PZENITH,PVEG,PZ0,PLAI,PTG1,&
                                 PPATCH,                    &
                                 PSW_BANDS,                 &
                                 PALBNIR_VEG,PALBVIS_VEG,   &
                                 PALBUV_VEG,                &
                                 PALBNIR_SOIL,PALBVIS_SOIL, &
                                 PALBUV_SOIL,               &
                                 PEMIS_ECO,                 &
                                 TPSNOW,                    &
                                 PALBNIR_ECO,PALBVIS_ECO,   &
                                 PALBUV_ECO,                &
                                 PDIR_ALB,PSCA_ALB,         &
                                 PEMIS,PTSRAD,PTSURF        )  
!     ###################################################
!
!!**** ** computes radiative fields used in ISBA
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
!                            <= to ZSNG
!!     B. Decharme  2013    new coupling variable and optimization    
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_SURF_PAR,  ONLY : XUNDEF
!
USE MODD_TYPE_SNOW
!
USE MODD_ISBA_n,    ONLY : XPSN,XFF,XEMISF
!
USE MODI_ALBEDO
USE MODI_AVERAGE_RAD
USE MODI_UPDATE_RAD_ISBA_n
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
LOGICAL,                INTENT(IN)   :: OFLOOD
 CHARACTER(LEN=4),       INTENT(IN)   :: HALBEDO     ! albedo type
! Albedo dependance with surface soil water content
!   "EVOL" = albedo evolves with soil wetness
!   "DRY " = constant albedo value for dry soil
!   "WET " = constant albedo value for wet soil
!   "MEAN" = constant albedo value for medium soil wetness
!
REAL, DIMENSION(:,:),   INTENT(IN)   :: PVEG        ! vegetation fraction
REAL, DIMENSION(:,:),   INTENT(IN)   :: PZ0         ! roughness length
REAL, DIMENSION(:,:),   INTENT(IN)   :: PLAI        ! leaf area index
REAL, DIMENSION(:,:),   INTENT(IN)   :: PTG1        ! soil surface temperature
REAL, DIMENSION(:,:),   INTENT(IN)   :: PPATCH      ! tile fraction
REAL, DIMENSION(:),     INTENT(IN)   :: PSW_BANDS   ! middle wavelength of each band
REAL, DIMENSION(:),     INTENT(IN)   :: PZENITH     

REAL, DIMENSION(:,:),   INTENT(IN)   :: PALBNIR_VEG ! near-infra-red albedo of vegetation
REAL, DIMENSION(:,:),   INTENT(IN)   :: PALBVIS_VEG ! visible albedo of vegetation
REAL, DIMENSION(:,:),   INTENT(IN)   :: PALBUV_VEG  ! UV albedo of vegetation
REAL, DIMENSION(:,:),   INTENT(IN)   :: PALBNIR_SOIL! near-infra-red albedo of soil
REAL, DIMENSION(:,:),   INTENT(IN)   :: PALBVIS_SOIL! visible albedo of soil
REAL, DIMENSION(:,:),   INTENT(IN)   :: PALBUV_SOIL ! UV albedo of soil
REAL, DIMENSION(:,:),   INTENT(IN)   :: PEMIS_ECO   ! emissivity (soil+vegetation)
TYPE(SURF_SNOW),        INTENT(IN)   :: TPSNOW      ! prognostic snow cover
!
REAL, DIMENSION(:,:),   INTENT(OUT)  :: PALBNIR_ECO ! near-infra-red albedo (soil+vegetation)
REAL, DIMENSION(:,:),   INTENT(OUT)  :: PALBVIS_ECO ! visible albedo (soil+vegetation)
REAL, DIMENSION(:,:),   INTENT(OUT)  :: PALBUV_ECO  ! UV albedo (soil+vegetation)
!
REAL, DIMENSION(:,:),   INTENT(OUT)  :: PDIR_ALB    ! averaged direct albedo  (per wavelength)
REAL, DIMENSION(:,:),   INTENT(OUT)  :: PSCA_ALB    ! averaged diffuse albedo (per wavelength)
REAL, DIMENSION(:),     INTENT(OUT)  :: PEMIS       ! averaged emissivity
REAL, DIMENSION(:),     INTENT(OUT)  :: PTSRAD      ! averaged radiaitve temp.
REAL, DIMENSION(:),     INTENT(OUT)  :: PTSURF      ! surface effective temperature         (K)
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
!
REAL, DIMENSION(SIZE(PALBNIR_VEG,1),SIZE(PSW_BANDS),SIZE(PALBVIS_VEG,2)) :: ZDIR_ALB_PATCH 
!                                                     ! direct albedo
REAL, DIMENSION(SIZE(PALBNIR_VEG,1),SIZE(PSW_BANDS),SIZE(PALBVIS_VEG,2)) :: ZSCA_ALB_PATCH 
!                                                     ! diffuse albedo
REAL, DIMENSION(SIZE(PEMIS_ECO,  1),SIZE(PALBVIS_VEG,2)) :: ZEMIS_PATCH   ! emissivity with snow-flood
REAL, DIMENSION(SIZE(PEMIS_ECO,  1),SIZE(PALBVIS_VEG,2)) :: ZTSRAD_PATCH  ! Tsrad
REAL, DIMENSION(SIZE(PEMIS_ECO,  1),SIZE(PALBVIS_VEG,2)) :: ZTSURF_PATCH  ! Tsurf
REAL, DIMENSION(SIZE(PEMIS_ECO,  1),SIZE(PALBVIS_VEG,2)) :: ZEMIS         ! emissivity with flood
!
LOGICAL :: LEXPLICIT_SNOW ! snow scheme key
!
INTEGER :: INP, INI
INTEGER :: JP, JI ! loop on patches
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!*    0.      Init
!             ----
!
IF (LHOOK) CALL DR_HOOK('AVERAGED_ALBEDO_EMIS_ISBA',0,ZHOOK_HANDLE)
!
INI=SIZE(PPATCH,1)
INP=SIZE(PPATCH,2)
!
PDIR_ALB(:,:)=0.
PSCA_ALB(:,:)=0.
PEMIS   (:)  =0.
PTSRAD  (:)  =0.
PTSURF  (:)  =0.
!
ZDIR_ALB_PATCH(:,:,:)=0.
ZSCA_ALB_PATCH(:,:,:)=0.
ZEMIS_PATCH   (:,:  )=0.
!
LEXPLICIT_SNOW = (TPSNOW%SCHEME=='3-L'.OR.TPSNOW%SCHEME=='CRO')
!
ZTSRAD_PATCH (:,:) = PTG1(:,:)
ZTSURF_PATCH (:,:) = PTG1(:,:)
!
!
!*    1.      averaged albedo on natural continental surfaces (except prognostic snow)
!             -----------------------------------------------
!
 CALL ALBEDO(HALBEDO,                                    &
              PALBVIS_VEG,PALBNIR_VEG,PALBUV_VEG,PVEG,    &
              PALBVIS_SOIL,PALBNIR_SOIL,PALBUV_SOIL,      &
              PALBVIS_ECO,PALBNIR_ECO,PALBUV_ECO          )  

!
!*    2.      averaged albedo and emis. on natural continental surfaces (with prognostic snow)
!             ---------------------------------------------------------
!* Initialization of albedo for each wavelength, emissivity and snow/flood fractions
!
 CALL UPDATE_RAD_ISBA_n(OFLOOD, TPSNOW%SCHEME,PZENITH,PSW_BANDS,PVEG,PLAI, &
                         PZ0,PALBNIR_ECO,PALBVIS_ECO,PALBUV_ECO,PEMIS_ECO,&
                         ZDIR_ALB_PATCH,ZSCA_ALB_PATCH,ZEMIS_PATCH        )  
!
!* radiative surface temperature
!
ZEMIS(:,:)=PEMIS_ECO(:,:)
!
IF(LEXPLICIT_SNOW.AND.OFLOOD)THEN
  WHERE(XPSN(:,:)<1.0.AND.PEMIS_ECO(:,:)/=XUNDEF)
       ZEMIS(:,:) = ((1.-XFF(:,:)-XPSN(:,:))*PEMIS_ECO(:,:) + XFF(:,:)*XEMISF(:,:)) / (1.-XPSN(:,:))
  ENDWHERE
ENDIF
!
IF(LEXPLICIT_SNOW)THEN
  WHERE(PEMIS_ECO(:,:)/=XUNDEF.AND.ZEMIS_PATCH(:,:)/=0.)
       ZTSRAD_PATCH(:,:) = ( ( (1.-XPSN(:,:))*ZEMIS      (:,:)*PTG1     (:,:)**4            &
                             +     XPSN(:,:) *TPSNOW%EMIS(:,:)*TPSNOW%TS(:,:)**4 )   &
                           / ZEMIS_PATCH(:,:) )**0.25        
  ENDWHERE
ENDIF
!
!* averaged radiative fields
!
 CALL AVERAGE_RAD(PPATCH,                                                     &
                   ZDIR_ALB_PATCH, ZSCA_ALB_PATCH, ZEMIS_PATCH, ZTSRAD_PATCH, &
                   PDIR_ALB,       PSCA_ALB,       PEMIS,       PTSRAD        )  
!
!* averaged effective temperature
!
IF(LEXPLICIT_SNOW)THEN
  ZTSURF_PATCH(:,:) = PTG1(:,:)*(1.-XPSN(:,:)) + TPSNOW%TS(:,:)*XPSN(:,:)
ENDIF
!
DO JP=1,INP
  DO JI=1,INI
     PTSURF(JI) = PTSURF(JI) + PPATCH(JI,JP) * ZTSURF_PATCH(JI,JP)
  ENDDO
ENDDO
!
IF (LHOOK) CALL DR_HOOK('AVERAGED_ALBEDO_EMIS_ISBA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE AVERAGED_ALBEDO_EMIS_ISBA
