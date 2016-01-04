!     ################################################################
      SUBROUTINE UPDATE_ESM_ISBA_n (I, &
                                    KI,KSW,PZENITH,PSW_BANDS,PDIR_ALB,& 
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
!
USE MODD_ISBA_n, ONLY : ISBA_t
!
USE MODD_TYPE_SNOW
USE MODD_SURF_PAR, ONLY : XUNDEF
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
!
TYPE(ISBA_t), INTENT(INOUT) :: I
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
REAL, DIMENSION(KI,KSW,I%O%NPATCH) :: ZDIR_ALB_PATCH
REAL, DIMENSION(KI,KSW,I%O%NPATCH) :: ZSCA_ALB_PATCH
REAL, DIMENSION(KI,I%O%NPATCH)     :: ZEMIS_PATCH
REAL, DIMENSION(KI,I%O%NPATCH)     :: ZTSRAD_PATCH
REAL, DIMENSION(KI,I%O%NPATCH)     :: ZTSURF_PATCH
REAL, DIMENSION(KI,I%O%NPATCH)     :: ZEMIS          ! emissivity with flood
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
ZEMIS         (:,:  ) = I%M%T%XEMIS(:,:)
!
LEXPLICIT_SNOW = (I%R%TSNOW%SCHEME=='3-L'.OR.I%R%TSNOW%SCHEME=='CRO')
!
ZTSRAD_PATCH (:,:) = I%R%XTG(:,1,:)
ZTSURF_PATCH (:,:) = I%R%XTG(:,1,:)
!
!
!*       2.     Update nature albedo and emissivity
!               -----------------------------------
!
 CALL UPDATE_RAD_ISBA_n(I, &
                        I%O%LFLOOD,I%R%TSNOW%SCHEME,PZENITH,PSW_BANDS,I%M%T%XVEG,I%M%T%XLAI,I%M%T%XZ0, &
                         I%O%LMEB_PATCH,I%M%M%XLAIGV,I%M%M%XGNDLITTER,I%M%M%XZ0LITTER,I%M%M%XH_VEG,      &
                         I%M%T%XALBNIR,I%M%T%XALBVIS,I%M%T%XALBUV,I%M%T%XEMIS,                       &
                         ZDIR_ALB_PATCH,ZSCA_ALB_PATCH,ZEMIS_PATCH           )
!
!*       3.     radiative surface temperature
!               -----------------------------
!
IF(LEXPLICIT_SNOW.AND.I%O%LFLOOD)THEN
  WHERE(I%R%XPSN(:,:)<1.0.AND.I%M%T%XEMIS(:,:)/=XUNDEF)
       ZEMIS(:,:) = ((1.-I%I%XFF(:,:)-I%R%XPSN(:,:))*I%M%T%XEMIS(:,:) + &
                I%I%XFF(:,:)*I%I%XEMISF(:,:)) / (1.-I%R%XPSN(:,:))
  ENDWHERE
ENDIF
!
IF(LEXPLICIT_SNOW)THEN
  WHERE(I%M%T%XEMIS(:,:)/=XUNDEF.AND.ZEMIS_PATCH(:,:)/=0.)
       ZTSRAD_PATCH(:,:) = ( ( (1.-I%R%XPSN(:,:))*ZEMIS     (:,:)*I%R%XTG   (:,1,:)**4     &
                             +     I%R%XPSN(:,:) *I%R%TSNOW%EMIS(:,:)*I%R%TSNOW%TS(:,:)**4 )   &
                           / ZEMIS_PATCH(:,:) )**0.25         
  ENDWHERE
ENDIF        
!
!
!*       4.     averaged fields
!               ---------------
!
 CALL AVERAGE_RAD(I%IP%XPATCH,                                                     &
                   ZDIR_ALB_PATCH, ZSCA_ALB_PATCH, ZEMIS_PATCH, ZTSRAD_PATCH, &
                   PDIR_ALB,       PSCA_ALB,       I%I%XEMIS_NAT,   I%R%XTSRAD_NAT    )  
!
PEMIS = I%I%XEMIS_NAT
PTSRAD = I%R%XTSRAD_NAT
!
!* averaged effective temperature
!
IF(LEXPLICIT_SNOW)THEN
  ZTSURF_PATCH(:,:) = I%R%XTG(:,1,:)*(1.-I%R%XPSN(:,:)) + I%R%TSNOW%TS(:,:)*I%R%XPSN(:,:)
ENDIF
!
 CALL AVERAGE_TSURF(I%IP%XPATCH, ZTSURF_PATCH, PTSURF)
!
IF (LHOOK) CALL DR_HOOK('UPDATE_ESM_ISBA_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE UPDATE_ESM_ISBA_n
