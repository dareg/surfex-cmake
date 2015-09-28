!#########
SUBROUTINE SFX_OASIS_RECV(HPROGRAM,KI,KSW,PTIMEC,OOASIS_PUT,     &
                          PLAND_WTD,PLAND_FWTD,PLAND_FFLOOD,     &
                          PLAND_PIFLOOD,                         &
                          PSEA_SST,PSEA_UCU,PSEA_VCU,            &
                          PSEAICE_SIT,PSEAICE_CVR,PSEAICE_ALB    )
!########################################
!
!!****  *SFX_OASIS_RECV* - Receive coupling fields from oasis
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
!!      B. Decharme   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2013
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODN_SFX_OASIS,  ONLY : XTSTEP_CPL_LAND
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE MODD_SFX_OASIS
!
USE MODI_GET_LUOUT
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
#ifdef SFXOASIS
USE MOD_OASIS
#endif
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=*),       INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
INTEGER,                INTENT(IN)  :: KI        ! number of points on this proc
INTEGER,                INTENT(IN)  :: KSW       ! number of short-wave spectral bands
REAL,                   INTENT(IN)  :: PTIMEC    ! Cumulated run time step (s)
!
LOGICAL,                INTENT(OUT) :: OOASIS_PUT    ! To initialise all module
!
REAL, DIMENSION(KI),    INTENT(OUT) :: PLAND_WTD     ! Land water table depth (m)
REAL, DIMENSION(KI),    INTENT(OUT) :: PLAND_FWTD    ! Land grid-cell fraction of water table rise (-)
REAL, DIMENSION(KI),    INTENT(OUT) :: PLAND_FFLOOD  ! Land Floodplains fraction (-)
REAL, DIMENSION(KI),    INTENT(OUT) :: PLAND_PIFLOOD ! Land Potential flood infiltration (kg/m2)
!
REAL, DIMENSION(KI),    INTENT(OUT) :: PSEA_SST ! Sea surface temperature (K)
REAL, DIMENSION(KI),    INTENT(OUT) :: PSEA_UCU ! Sea u-current stress (Pa)
REAL, DIMENSION(KI),    INTENT(OUT) :: PSEA_VCU ! Sea v-current stress (Pa)
!
REAL, DIMENSION(KI),    INTENT(OUT) :: PSEAICE_SIT ! Sea-ice Temperature (K)
REAL, DIMENSION(KI),    INTENT(OUT) :: PSEAICE_CVR ! Sea-ice cover (-)
REAL, DIMENSION(KI),    INTENT(OUT) :: PSEAICE_ALB ! Sea-ice albedo (-)
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL, DIMENSION(KI,1) :: ZREAD
!
INTEGER               :: IDATE  ! current coupling time step (s)
INTEGER               :: IERR   ! Error info
INTEGER               :: ILUOUT
CHARACTER(LEN=50)     :: YCOMMENT
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
#ifdef SFXOASIS
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SFX_OASIS_RECV',0,ZHOOK_HANDLE)
!
!*       1.     Initialize :
!               ------------
!
CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
IDATE = INT(PTIMEC)
!
OOASIS_PUT = .FALSE.
!
!-------------------------------------------------------------------------------
!
!*       2.     Get Land surface variable :
!               ------------------------------------
!
IF(LCPL_LAND)THEN
!
! * Init river input fields
!
  ZREAD(:,:) = XUNDEF
!
  PLAND_WTD    (:) = XUNDEF
  PLAND_FWTD   (:) = XUNDEF
  PLAND_FFLOOD (:) = XUNDEF
  PLAND_PIFLOOD(:) = XUNDEF
!
! * Receive river input fields
!
  IF(LCPL_GW)THEN
!
    YCOMMENT='water table depth'
    CALL OASIS_GET(NWTD_ID,IDATE,ZREAD(:,:),IERR)
    CALL CHECK_RECV(ILUOUT,IERR,YCOMMENT,OOASIS_PUT)
    PLAND_WTD(:)=ZREAD(:,1)
!
    YCOMMENT='fraction of water table rise'
    CALL OASIS_GET(NFWTD_ID,IDATE,ZREAD(:,:),IERR)
    CALL CHECK_RECV(ILUOUT,IERR,YCOMMENT,OOASIS_PUT)
    PLAND_FWTD(:)=ZREAD(:,1)
!
  ENDIF
!
  IF(LCPL_FLOOD)THEN
!
    YCOMMENT='Flood fraction'
    CALL OASIS_GET(NFFLOOD_ID,IDATE,ZREAD(:,:),IERR)
    CALL CHECK_RECV(ILUOUT,IERR,YCOMMENT,OOASIS_PUT)
    PLAND_FFLOOD(:)=ZREAD(:,1)
!
    YCOMMENT='Potential flood infiltration'
    CALL OASIS_GET(NPIFLOOD_ID,IDATE,ZREAD(:,:),IERR)
    CALL CHECK_RECV(ILUOUT,IERR,YCOMMENT,OOASIS_PUT)
    PLAND_PIFLOOD(:)=ZREAD(:,1)*XTSTEP_CPL_LAND
!
  ENDIF
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       3.     Get Sea variables :
!               -----------------------------
!
!
IF(LCPL_SEA)THEN
!
! * Init ocean input fields
!
  ZREAD(:,:) = XUNDEF
!
  PSEA_SST (:) = XUNDEF
  PSEA_UCU (:) = XUNDEF
  PSEA_VCU (:) = XUNDEF
!
  PSEAICE_SIT (:) = XUNDEF
  PSEAICE_CVR (:) = XUNDEF
  PSEAICE_ALB (:) = XUNDEF
!
! * Receive ocean input fields
!
  YCOMMENT='Sea surface temperature'
  CALL OASIS_GET(NSEA_SST_ID,IDATE,ZREAD(:,:),IERR)
  CALL CHECK_RECV(ILUOUT,IERR,YCOMMENT,OOASIS_PUT)
  PSEA_SST(:)=ZREAD(:,1)
!
  YCOMMENT='Sea u-current stress'
  CALL OASIS_GET(NSEA_UCU_ID,IDATE,ZREAD(:,:),IERR)
  CALL CHECK_RECV(ILUOUT,IERR,YCOMMENT,OOASIS_PUT)
  PSEA_UCU(:)=ZREAD(:,1)
!
  YCOMMENT='Sea v-current stress'
  CALL OASIS_GET(NSEA_VCU_ID,IDATE,ZREAD(:,:),IERR)
  CALL CHECK_RECV(ILUOUT,IERR,YCOMMENT,OOASIS_PUT)
  PSEA_VCU(:)=ZREAD(:,1)
!
  IF(LCPL_SEAICE)THEN
!
    YCOMMENT='Sea-ice Temperature'
    CALL OASIS_GET(NSEAICE_SIT_ID,IDATE,ZREAD(:,:),IERR)
    CALL CHECK_RECV(ILUOUT,IERR,YCOMMENT,OOASIS_PUT)
    PSEAICE_SIT(:)=ZREAD(:,1)
!
    YCOMMENT='Sea-ice cover'
    CALL OASIS_GET(NSEAICE_CVR_ID,IDATE,ZREAD(:,:),IERR)
    CALL CHECK_RECV(ILUOUT,IERR,YCOMMENT,OOASIS_PUT)
    PSEAICE_CVR(:)=ZREAD(:,1)
!
    YCOMMENT='Sea-ice albedo'
    CALL OASIS_GET(NSEAICE_ALB_ID,IDATE,ZREAD(:,:),IERR)
    CALL CHECK_RECV(ILUOUT,IERR,YCOMMENT,OOASIS_PUT)
    PSEAICE_ALB(:)=ZREAD(:,1)
!
  ENDIF
!
ENDIF
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SFX_OASIS_RECV',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!
SUBROUTINE CHECK_RECV(KLUOUT,KERR,HCOMMENT,OPUT)
!
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
INTEGER,          INTENT(IN)  :: KLUOUT
INTEGER,          INTENT(IN)  :: KERR
CHARACTER(LEN=*), INTENT(IN)  :: HCOMMENT
LOGICAL,          INTENT(OUT) :: OPUT
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SFX_OASIS_RECV:CHECK_RECV',0,ZHOOK_HANDLE)
!
IF (KERR/=OASIS_OK.AND.KERR<OASIS_RECVD) THEN
   WRITE(KLUOUT,'(A,I4)')'Return OASIS code receiving '//TRIM(HCOMMENT)//' : ',KERR
   CALL ABOR1_SFX('SFX_OASIS_RECV: problem receiving '//TRIM(HCOMMENT)//' from OASIS')
ENDIF
!
IF (KERR>=OASIS_RECVD) THEN
   OPUT = .TRUE.
ENDIF 
!
IF (LHOOK) CALL DR_HOOK('SFX_OASIS_RECV:CHECK_RECV',1,ZHOOK_HANDLE)
!
END SUBROUTINE CHECK_RECV
!
!-------------------------------------------------------------------------------
#endif
!-------------------------------------------------------------------------------
!
END SUBROUTINE SFX_OASIS_RECV
