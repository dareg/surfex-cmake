!     #########
      SUBROUTINE AVERAGED_ALBEDO_EMIS_ISBA (IO, MT, MM, MA, IP, I, IR, PVEGTYPE, &
                                 PZENITH, PTG1, PSW_BANDS, PDIR_ALB, PSCA_ALB,     &
                                 PEMIS, PTSRAD, PTSURF, PDIR_SW, PSCA_SW           )
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
!!     P. Samuelsson 10/2014 MEB
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_TIME_t, ISBA_PARAM_MEB_t, ISBA_PARAM_ALB_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t, ISBA_INIT_t
USE MODD_ISBA_n, ONLY: ISBA_PROG_t
!
USE MODD_SURF_PAR,  ONLY : XUNDEF
!
USE MODD_TYPE_SNOW
!
USE MODD_CSTS,      ONLY : XSTEFAN
USE MODE_MEB,       ONLY : MEBPALPHAN
!
USE MODI_ALBEDO
USE MODI_AVERAGE_RAD
USE MODI_UPDATE_RAD_ISBA_n
USE MODI_ISBA_LWNET_MEB
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
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_PARAM_TIME_t), INTENT(INOUT) :: MT
TYPE(ISBA_PARAM_MEB_t), INTENT(INOUT) :: MM
TYPE(ISBA_PARAM_ALB_t), INTENT(INOUT) :: MA
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
TYPE(ISBA_INIT_t), INTENT(INOUT) :: I
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
!
REAL, DIMENSION(:,:), INTENT(IN) :: PVEGTYPE
REAL, DIMENSION(:,:),   INTENT(IN)   :: PTG1        ! soil surface temperature
REAL, DIMENSION(:),     INTENT(IN)   :: PZENITH     
REAL, DIMENSION(:),     INTENT(IN)   :: PSW_BANDS   ! middle wavelength of each band
!
REAL, DIMENSION(:,:),   INTENT(OUT)  :: PDIR_ALB    ! averaged direct albedo  (per wavelength)
REAL, DIMENSION(:,:),   INTENT(OUT)  :: PSCA_ALB    ! averaged diffuse albedo (per wavelength)
REAL, DIMENSION(:),     INTENT(OUT)  :: PEMIS       ! averaged emissivity
REAL, DIMENSION(:),     INTENT(OUT)  :: PTSRAD      ! averaged radiaitve temp.
REAL, DIMENSION(:),     INTENT(OUT)  :: PTSURF      ! surface effective temperature         (K)
!
REAL, DIMENSION(:,:),   INTENT(IN), OPTIONAL   :: PDIR_SW ! Downwelling direct SW radiation
REAL, DIMENSION(:,:),   INTENT(IN), OPTIONAL   :: PSCA_SW ! Downwelling diffuse SW radiation
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
!
REAL, DIMENSION(SIZE(MT%XALBNIR_VEG,1),SIZE(PSW_BANDS),SIZE(MT%XALBVIS_VEG,2)) :: ZDIR_ALB_PATCH 
!                                                     ! direct albedo
REAL, DIMENSION(SIZE(MT%XALBNIR_VEG,1),SIZE(PSW_BANDS),SIZE(MT%XALBVIS_VEG,2)) :: ZSCA_ALB_PATCH 
!                                                     ! diffuse albedo
REAL, DIMENSION(SIZE(MT%XEMIS,  1),SIZE(MT%XALBVIS_VEG,2)) :: ZEMIS_PATCH   ! emissivity with snow-flood
REAL, DIMENSION(SIZE(MT%XEMIS,  1),SIZE(MT%XALBVIS_VEG,2)) :: ZTSRAD_PATCH  ! Tsrad
REAL, DIMENSION(SIZE(MT%XEMIS,  1),SIZE(MT%XALBVIS_VEG,2)) :: ZTSURF_PATCH  ! Tsurf
REAL, DIMENSION(SIZE(MT%XEMIS,  1)) :: ZEMIS         ! emissivity with flood
!
REAL, DIMENSION(SIZE(MT%XEMIS,  1)) :: ZSNOWDEPTH    ! Total snow depth
REAL, DIMENSION(SIZE(MT%XEMIS,  1)) :: ZPALPHAN      ! Snow/canopy ratio factor 
REAL, DIMENSION(SIZE(MT%XEMIS,  1)) :: ZLW_RAD       ! Fake downwelling LW rad
REAL, DIMENSION(SIZE(MT%XEMIS,  1)) :: ZLW_UP        ! Upwelling LW rad
REAL, DIMENSION(SIZE(MT%XEMIS,  1)) :: ZLWNET_N      ! LW net for snow surface
REAL, DIMENSION(SIZE(MT%XEMIS,  1)) :: ZLWNET_V      ! LW net for canopy veg
REAL, DIMENSION(SIZE(MT%XEMIS,  1)) :: ZLWNET_G      ! LW net for ground
REAL, DIMENSION(SIZE(MT%XEMIS,  1)) :: ZDUMMY
REAL, DIMENSION(SIZE(MT%XEMIS,  1)) :: ZEMISF
REAL, DIMENSION(SIZE(MT%XEMIS,  1)) :: ZFF
!
LOGICAL :: LEXPLICIT_SNOW ! snow scheme key
!
INTEGER :: INP, INI
INTEGER :: JP, JI ! loop on patches
INTEGER :: JPATCH ! loop on patches
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!*    0.      Init
!             ----
!
IF (LHOOK) CALL DR_HOOK('AVERAGED_ALBEDO_EMIS_ISBA',0,ZHOOK_HANDLE)
!
INI=SIZE(IP%XPATCH,1)
INP=SIZE(IP%XPATCH,2)
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
LEXPLICIT_SNOW = (IR%TSNOW%SCHEME=='3-L'.OR.IR%TSNOW%SCHEME=='CRO')
!
ZTSRAD_PATCH (:,:) = PTG1(:,:)
ZTSURF_PATCH (:,:) = PTG1(:,:)
!
!
!*    1.      averaged albedo on natural continental surfaces (except prognostic snow)
!             -----------------------------------------------
!
 CALL ALBEDO(IO%CALBEDO, MT, MA )
!
!*    2.      averaged albedo and emis. on natural continental surfaces (with prognostic snow)
!             ---------------------------------------------------------
!
! A dummy downwelling LW radiation can be used for calculation of radiative surface temp 
!
ZLW_RAD(:) = 300.0
!    
!* Initialization of albedo for each wavelength, emissivity and snow/flood fractions
!
IF(PRESENT(PDIR_SW))THEN
!
! For the case when MEB patch albedo is requested downweeling SW is needed
!
  CALL UPDATE_RAD_ISBA_n(IO, MT, MM, MA, IP, I, IR, PVEGTYPE, &
                         PZENITH,PSW_BANDS,ZDIR_ALB_PATCH,ZSCA_ALB_PATCH,ZEMIS_PATCH,        &
                         PDIR_SW, PSCA_SW    )
ELSE
!
! For cases when MEB patch albedo is not requested no downweeling SW is needed
!
  CALL UPDATE_RAD_ISBA_n(IO, MT, MM, MA, IP, I, IR, PVEGTYPE, &
                         PZENITH,PSW_BANDS, ZDIR_ALB_PATCH,ZSCA_ALB_PATCH,ZEMIS_PATCH  )
ENDIF
!
!
!* radiative surface temperature
!
DO JPATCH=1,SIZE(MT%XALBVIS_VEG,2)
!
  IF(IO%LMEB_PATCH(JPATCH))THEN  ! MEB patches
!
!   ZPALPHAN is needed as input to ISBA_LWNET_MEB
!
    ZSNOWDEPTH(:) = SUM(IR%TSNOW%WSNOW(:,:,JPATCH)/IR%TSNOW%RHO(:,:,JPATCH),2)
    ZPALPHAN(:)   = MEBPALPHAN(ZSNOWDEPTH,MM%XH_VEG(:,JPATCH))
!
!   ZLWNET_N,ZLWNET_V,ZLWNET_G are needed for ZLW_UP and ZTSRAD_PATCH
!
    IF(IO%LFLOOD)THEN
      ZEMISF(:) = I%XEMISF(:,JPATCH)
      ZFF   (:) = I%XFF   (:,JPATCH)
    ELSE
      ZEMISF(:) = XUNDEF
      ZFF   (:) = 0.0
    ENDIF
!
    CALL ISBA_LWNET_MEB(MT%XLAI(:,JPATCH),IR%XPSN(:,JPATCH),ZPALPHAN,   &
                     IR%TSNOW%EMIS(:,JPATCH),ZEMISF(:),ZFF(:),          &
                     IR%XTV(:,JPATCH),PTG1(:,JPATCH),IR%TSNOW%TS(:,JPATCH),  &
                     ZLW_RAD,ZLWNET_N,ZLWNET_V,ZLWNET_G,                     &
                     ZDUMMY,ZDUMMY,ZDUMMY, ZDUMMY,ZDUMMY,ZDUMMY,             &
                     ZDUMMY,ZDUMMY,ZDUMMY, ZDUMMY,ZDUMMY,ZDUMMY            )
!
    ZLW_UP(:)   = ZLW_RAD(:) - (ZLWNET_V(:) + ZLWNET_G(:) + ZLWNET_N(:))
!
!   MEB patch radiative temperature
!
    WHERE (ZEMIS_PATCH(:,JPATCH)/=0.)
      ZTSRAD_PATCH(:,JPATCH) = ((ZLW_UP(:) - ZLW_RAD(:)*(1.0-ZEMIS_PATCH(:,JPATCH)))/ &
                              (XSTEFAN*ZEMIS_PATCH(:,JPATCH)))**0.25
    END WHERE
!
  ELSE   ! Non-MEB patches

    ZEMIS(:) = MT%XEMIS(:,JPATCH)
!
    IF(IO%LFLOOD.AND.LEXPLICIT_SNOW)THEN
      WHERE(IR%XPSN(:,JPATCH)<1.0.AND.MT%XEMIS(:,JPATCH)/=XUNDEF)          
        ZEMIS(:) = ((1.-I%XFF(:,JPATCH)-IR%XPSN(:,JPATCH))*MT%XEMIS(:,JPATCH) + I%XFF(:,JPATCH)*I%XEMISF(:,JPATCH)) &
                   /(1.-IR%XPSN(:,JPATCH))
      ENDWHERE   
    ENDIF
!
    IF (.NOT.LEXPLICIT_SNOW) THEN
      ZTSRAD_PATCH(:,JPATCH) = PTG1(:,JPATCH)
    ELSE IF (LEXPLICIT_SNOW) THEN
      WHERE (MT%XEMIS(:,JPATCH)/=XUNDEF .AND. ZEMIS_PATCH(:,JPATCH)/=0.)
        ZTSRAD_PATCH(:,JPATCH) =( ( (1.-IR%XPSN(:,JPATCH))*ZEMIS      (:)       *PTG1     (:,JPATCH)**4            &
                                    +    IR%XPSN(:,JPATCH) *IR%TSNOW%EMIS(:,JPATCH)*IR%TSNOW%TS(:,JPATCH)**4 ) )**0.25  &
                                 / ZEMIS_PATCH(:,JPATCH)**0.25  
      END WHERE
    END IF
  ENDIF
END DO
!
!* averaged radiative fields
!
 CALL AVERAGE_RAD(IP%XPATCH,                                                     &
                   ZDIR_ALB_PATCH, ZSCA_ALB_PATCH, ZEMIS_PATCH, ZTSRAD_PATCH, &
                   PDIR_ALB,       PSCA_ALB,       PEMIS,       PTSRAD        )  
!
!* averaged effective temperature
!
IF(LEXPLICIT_SNOW)THEN
  ZTSURF_PATCH(:,:) = PTG1(:,:)*(1.-IR%XPSN(:,:)) + IR%TSNOW%TS(:,:)*IR%XPSN(:,:)
ENDIF
!
DO JP=1,INP
  DO JI=1,INI
     PTSURF(JI) = PTSURF(JI) + IP%XPATCH(JI,JP) * ZTSURF_PATCH(JI,JP)
  ENDDO
ENDDO
!
IF (LHOOK) CALL DR_HOOK('AVERAGED_ALBEDO_EMIS_ISBA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE AVERAGED_ALBEDO_EMIS_ISBA
