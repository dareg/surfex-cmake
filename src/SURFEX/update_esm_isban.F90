!     ################################################################
      SUBROUTINE UPDATE_ESM_ISBA_n(KI,KSW,PZENITH,PSW_BANDS,PDIR_ALB,& 
                                   PSCA_ALB,PEMIS,PTSRAD,PTSURF      )
!     ################################################################
!
!!****  *UPDATE_ESM_ISBA_n* - update ISBA radiative and physical properties in Earth System Model 
!!                            after the call to OASIS coupler in order 
!!                            to close the energy budget between radiative scheme and surfex
!!
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!     B. Decharme 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    09/2009
!!      B. Decharme 06/2013 new coupling variables
!!      P. Samuelsson 10/2014 MEB
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_TYPE_SNOW
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_ISBA_n,   ONLY : NPATCH,XTG,TSNOW,XPSN,XVEG,XLAI,XZ0,    &
                          XALBNIR,XALBVIS,XALBUV,XEMIS,XPATCH,    &
                          LFLOOD,XFF,XEMISF,XEMIS_NAT,XTSRAD_NAT, &
                          LMEB_PATCH,XLAIGV,XVEGGV,XZ0GV, XH_VEG
!
USE MODI_AVERAGE_RAD
USE MODI_AVERAGE_TSURF
USE MODI_UPDATE_RAD_ISBA_n
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER,                            INTENT(IN)  :: KI        ! number of points
INTEGER,                            INTENT(IN)  :: KSW       ! number of short-wave spectral bands
!
REAL,             DIMENSION(KI),    INTENT(IN)  :: PZENITH   ! solar zenithal angle
REAL,             DIMENSION(KSW),   INTENT(IN)  :: PSW_BANDS ! short-wave spectral bands
!
REAL,             DIMENSION(KI,KSW),INTENT(OUT) :: PDIR_ALB  ! direct albedo for each band
REAL,             DIMENSION(KI,KSW),INTENT(OUT) :: PSCA_ALB  ! diffuse albedo for each band
REAL,             DIMENSION(KI),    INTENT(OUT) :: PEMIS     ! emissivity
REAL,             DIMENSION(KI),    INTENT(OUT) :: PTSRAD    ! radiative temperature
REAL,             DIMENSION(KI),    INTENT(OUT) :: PTSURF    ! surface effective temperature         (K)
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL, DIMENSION(KI,KSW,NPATCH) :: ZDIR_ALB_PATCH
REAL, DIMENSION(KI,KSW,NPATCH) :: ZSCA_ALB_PATCH
REAL, DIMENSION(KI,NPATCH)     :: ZEMIS_PATCH
REAL, DIMENSION(KI,NPATCH)     :: ZTSRAD_PATCH
REAL, DIMENSION(KI,NPATCH)     :: ZTSURF_PATCH
REAL, DIMENSION(KI,NPATCH)     :: ZEMIS          ! emissivity with flood
!
LOGICAL :: LEXPLICIT_SNOW ! snow scheme key
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*       1.     Defaults
!               --------
!
IF (LHOOK) CALL DR_HOOK('UPDATE_ESM_ISBA_N',0,ZHOOK_HANDLE)
!
ZDIR_ALB_PATCH(:,:,:) = 0.0
ZSCA_ALB_PATCH(:,:,:) = 0.0
ZEMIS_PATCH   (:,:  ) = 0.0
ZEMIS         (:,:  ) = XEMIS(:,:)
!
LEXPLICIT_SNOW = (TSNOW%SCHEME=='3-L'.OR.TSNOW%SCHEME=='CRO')
!
ZTSRAD_PATCH (:,:) = XTG(:,1,:)
ZTSURF_PATCH (:,:) = XTG(:,1,:)
!
!
!*       2.     Update nature albedo and emissivity
!               -----------------------------------
!
 CALL UPDATE_RAD_ISBA_n(LFLOOD,TSNOW%SCHEME,PZENITH,PSW_BANDS,XVEG,XLAI,XZ0, &
                         LMEB_PATCH,XLAIGV,XVEGGV,XZ0GV,XH_VEG,              &
                         XALBNIR,XALBVIS,XALBUV,XEMIS,                       &
                         ZDIR_ALB_PATCH,ZSCA_ALB_PATCH,ZEMIS_PATCH           )
!
!*       3.     radiative surface temperature
!               -----------------------------
!
IF(LEXPLICIT_SNOW.AND.LFLOOD)THEN
  WHERE(XPSN(:,:)<1.0.AND.XEMIS(:,:)/=XUNDEF)
       ZEMIS(:,:) = ((1.-XFF(:,:)-XPSN(:,:))*XEMIS(:,:) + XFF(:,:)*XEMISF(:,:)) / (1.-XPSN(:,:))
  ENDWHERE
ENDIF
!
IF(LEXPLICIT_SNOW)THEN
  WHERE(XEMIS(:,:)/=XUNDEF.AND.ZEMIS_PATCH(:,:)/=0.)
       ZTSRAD_PATCH(:,:) = ( ( (1.-XPSN(:,:))*ZEMIS     (:,:)*XTG   (:,1,:)**4     &
                             +     XPSN(:,:) *TSNOW%EMIS(:,:)*TSNOW%TS(:,:)**4 )   &
                           / ZEMIS_PATCH(:,:) )**0.25         
  ENDWHERE
ENDIF        
!
!
!*       4.     averaged fields
!               ---------------
!
 CALL AVERAGE_RAD(XPATCH,                                                     &
                   ZDIR_ALB_PATCH, ZSCA_ALB_PATCH, ZEMIS_PATCH, ZTSRAD_PATCH, &
                   PDIR_ALB,       PSCA_ALB,       XEMIS_NAT,   XTSRAD_NAT    )  
!
PEMIS = XEMIS_NAT
PTSRAD = XTSRAD_NAT
!
!* averaged effective temperature
!
IF(LEXPLICIT_SNOW)THEN
  ZTSURF_PATCH(:,:) = XTG(:,1,:)*(1.-XPSN(:,:)) + TSNOW%TS(:,:)*XPSN(:,:)
ENDIF
!
 CALL AVERAGE_TSURF(XPATCH, ZTSURF_PATCH, PTSURF)
!
IF (LHOOK) CALL DR_HOOK('UPDATE_ESM_ISBA_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE UPDATE_ESM_ISBA_n
