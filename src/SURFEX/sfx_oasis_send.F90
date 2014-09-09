!#########
SUBROUTINE SFX_OASIS_SEND(KLUOUT,KI,KDATE,OSEND_LAND,OSEND_LAKE,OSEND_SEA,      &
                          PLAND_RUNOFF,PLAND_DRAIN,PLAND_CALVING,PLAND_RECHARGE,&
                          PLAND_PFLOOD,PLAND_EFLOOD,PLAND_IFLOOD,               &
                          PLAKE_EVAP,PLAKE_RAIN,PLAKE_SNOW,                     &
                          PSEA_FWSU,PSEA_FWSV,PSEA_HEAT,PSEA_SNET,PSEA_WIND,    &
                          PSEA_FWSM,PSEA_EVAP,PSEA_RAIN,PSEA_SNOW,PSEA_WATF,    &
                          PSEAICE_HEAT,PSEAICE_SNET,PSEAICE_EVAP                )
!###########################################
!
!!****  *SFX_OASIS_SEND* - Send coupling fields
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
!!	B. Decharme   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2013
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_PAR,   ONLY : XUNDEF, NUNDEF
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
INTEGER,             INTENT(IN) :: KLUOUT
INTEGER,             INTENT(IN) :: KI            ! number of points
INTEGER,             INTENT(IN) :: KDATE  ! current coupling time step (s)
LOGICAL,             INTENT(IN) :: OSEND_LAND
LOGICAL,             INTENT(IN) :: OSEND_LAKE
LOGICAL,             INTENT(IN) :: OSEND_SEA
!
REAL, DIMENSION(KI), INTENT(IN) :: PLAND_RUNOFF    ! Cumulated Surface runoff             (kg/m2)
REAL, DIMENSION(KI), INTENT(IN) :: PLAND_DRAIN     ! Cumulated Deep drainage              (kg/m2)
REAL, DIMENSION(KI), INTENT(IN) :: PLAND_CALVING   ! Cumulated Calving flux               (kg/m2)
REAL, DIMENSION(KI), INTENT(IN) :: PLAND_RECHARGE  ! Cumulated Recharge to groundwater    (kg/m2)
REAL, DIMENSION(KI), INTENT(IN) :: PLAND_PFLOOD    ! Cumulated flood precip interception  (kg/m2)
REAL, DIMENSION(KI), INTENT(IN) :: PLAND_EFLOOD    ! Cumulated floodplains evaporation    (kg/m2)
REAL, DIMENSION(KI), INTENT(IN) :: PLAND_IFLOOD    ! Cumulated floodplains infiltration   (kg/m2)
!
REAL, DIMENSION(KI), INTENT(IN) :: PLAKE_EVAP  ! Cumulated Evaporation             (kg/m2)
REAL, DIMENSION(KI), INTENT(IN) :: PLAKE_RAIN  ! Cumulated Rainfall rate           (kg/m2)
REAL, DIMENSION(KI), INTENT(IN) :: PLAKE_SNOW  ! Cumulated Snowfall rate           (kg/m2)
!
REAL, DIMENSION(KI), INTENT(IN) :: PSEA_FWSU  ! Cumulated zonal wind stress       (Pa.s)
REAL, DIMENSION(KI), INTENT(IN) :: PSEA_FWSV  ! Cumulated meridian wind stress    (Pa.s)
REAL, DIMENSION(KI), INTENT(IN) :: PSEA_HEAT  ! Cumulated Non solar net heat flux (J/m2)
REAL, DIMENSION(KI), INTENT(IN) :: PSEA_SNET  ! Cumulated Solar net heat flux     (J/m2)
REAL, DIMENSION(KI), INTENT(IN) :: PSEA_WIND  ! Cumulated 10m wind speed          (m)
REAL, DIMENSION(KI), INTENT(IN) :: PSEA_FWSM  ! Cumulated wind stress             (Pa.s)
REAL, DIMENSION(KI), INTENT(IN) :: PSEA_EVAP  ! Cumulated Evaporation             (kg/m2)
REAL, DIMENSION(KI), INTENT(IN) :: PSEA_RAIN  ! Cumulated Rainfall rate           (kg/m2)
REAL, DIMENSION(KI), INTENT(IN) :: PSEA_SNOW  ! Cumulated Snowfall rate           (kg/m2)
REAL, DIMENSION(KI), INTENT(IN) :: PSEA_WATF  ! Cumulated outgoing freshwater flux (kg/m2)
!
REAL, DIMENSION(KI), INTENT(IN) :: PSEAICE_HEAT ! Cumulated Sea-ice non solar net heat flux (J/m2)
REAL, DIMENSION(KI), INTENT(IN) :: PSEAICE_SNET ! Cumulated Sea-ice solar net heat flux     (J/m2)
REAL, DIMENSION(KI), INTENT(IN) :: PSEAICE_EVAP ! Cumulated Sea-ice sublimation             (kg/m2)
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL, DIMENSION(KI,1) :: ZWRITE
!
CHARACTER(LEN=50)     :: YCOMMENT
INTEGER               :: IERR   ! Error info
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
#ifdef SFXOASIS
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SFX_OASIS_SEND',0,ZHOOK_HANDLE)
!
!*       1.     Initialize :
!               ------------
!
ZWRITE(:,:) = XUNDEF
!
!-------------------------------------------------------------------------------
!
!*       2.     Send land fields to OASIS:
!               --------------------------
!
IF(OSEND_LAND)THEN
!
! * Send river output fields
!
  YCOMMENT='Surface runoff over land'
  ZWRITE(:,1) = PLAND_RUNOFF(:)
  CALL OASIS_PUT(NRUNOFF_ID,KDATE,ZWRITE(:,:),IERR)
  CALL CHECK_SFX_SEND(KLUOUT,IERR,YCOMMENT)
!
  YCOMMENT='Deep drainage over land'
  ZWRITE(:,1) = PLAND_DRAIN(:)
  CALL OASIS_PUT(NDRAIN_ID,KDATE,ZWRITE(:,:),IERR)
  CALL CHECK_SFX_SEND(KLUOUT,IERR,YCOMMENT)
!
  IF(LCPL_CALVING)THEN
    YCOMMENT='calving flux over land'
    ZWRITE(:,1) = PLAND_CALVING(:)
    CALL OASIS_PUT(NCALVING_ID,KDATE,ZWRITE(:,:),IERR)
    CALL CHECK_SFX_SEND(KLUOUT,IERR,YCOMMENT)
  ENDIF
!
  IF(LCPL_GW)THEN
    YCOMMENT='groundwater recharge over land'
    ZWRITE(:,1) = PLAND_RECHARGE(:)
    CALL OASIS_PUT(NRECHARGE_ID,KDATE,ZWRITE(:,:),IERR)
    CALL CHECK_SFX_SEND(KLUOUT,IERR,YCOMMENT)
  ENDIF
!
  IF(LCPL_FLOOD)THEN
!          
    YCOMMENT='flood precip interception over land'
    ZWRITE(:,1) = PLAND_PFLOOD(:)
    CALL OASIS_PUT(NPFLOOD_ID,KDATE,ZWRITE(:,:),IERR)
    CALL CHECK_SFX_SEND(KLUOUT,IERR,YCOMMENT)
!          
    YCOMMENT='flood evaporation over land'
    ZWRITE(:,1) = PLAND_EFLOOD(:)
    CALL OASIS_PUT(NEFLOOD_ID,KDATE,ZWRITE(:,:),IERR)
    CALL CHECK_SFX_SEND(KLUOUT,IERR,YCOMMENT)
!          
    YCOMMENT='flood infiltration over land'
    ZWRITE(:,1) = PLAND_IFLOOD(:)
    CALL OASIS_PUT(NIFLOOD_ID,KDATE,ZWRITE(:,:),IERR)
    CALL CHECK_SFX_SEND(KLUOUT,IERR,YCOMMENT)
!    
  ENDIF
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       3.     Send lake fields to OASIS :
!               --------------------------
IF(OSEND_LAKE)THEN
!
! * Send output fields
!
  YCOMMENT='Evaporation over lake'
  ZWRITE(:,1) = PLAKE_EVAP(:)
  CALL OASIS_PUT(NLAKE_EVAP_ID,KDATE,ZWRITE(:,:),IERR)
  CALL CHECK_SFX_SEND(KLUOUT,IERR,YCOMMENT)
!
  YCOMMENT='Rainfall rate over lake'
  ZWRITE(:,1) = PLAKE_RAIN(:)
  CALL OASIS_PUT(NLAKE_RAIN_ID,KDATE,ZWRITE(:,:),IERR)
  CALL CHECK_SFX_SEND(KLUOUT,IERR,YCOMMENT)
!
  YCOMMENT='Snowfall rate over lake'
  ZWRITE(:,1) = PLAKE_SNOW(:)
  CALL OASIS_PUT(NLAKE_SNOW_ID,KDATE,ZWRITE(:,:),IERR)
  CALL CHECK_SFX_SEND(KLUOUT,IERR,YCOMMENT)
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       4.     Send sea fields to OASIS :
!               --------------------------
!
IF(OSEND_SEA)THEN
!
! * Send sea output fields
!
  IF(NSEA_FWSU_ID/=NUNDEF)THEN
    YCOMMENT='zonal wind stress over sea'
    ZWRITE(:,1) = PSEA_FWSU(:)
    CALL OASIS_PUT(NSEA_FWSU_ID,KDATE,ZWRITE(:,:),IERR)
    CALL CHECK_SFX_SEND(KLUOUT,IERR,YCOMMENT)
  ENDIF
!
  IF(NSEA_FWSV_ID/=NUNDEF)THEN
    YCOMMENT='meridian wind stress over sea'
    ZWRITE(:,1) = PSEA_FWSV(:)
    CALL OASIS_PUT(NSEA_FWSV_ID,KDATE,ZWRITE(:,:),IERR)
    CALL CHECK_SFX_SEND(KLUOUT,IERR,YCOMMENT)
  ENDIF
!
  IF(NSEA_HEAT_ID/=NUNDEF)THEN
    YCOMMENT='Non solar net heat flux over sea'
    ZWRITE(:,1) = PSEA_HEAT(:)
    CALL OASIS_PUT(NSEA_HEAT_ID,KDATE,ZWRITE(:,:),IERR)
    CALL CHECK_SFX_SEND(KLUOUT,IERR,YCOMMENT)
  ENDIF
!
  IF(NSEA_SNET_ID/=NUNDEF)THEN
    YCOMMENT='Solar net heat flux over sea'
    ZWRITE(:,1) = PSEA_SNET(:)
    CALL OASIS_PUT(NSEA_SNET_ID,KDATE,ZWRITE(:,:),IERR)
    CALL CHECK_SFX_SEND(KLUOUT,IERR,YCOMMENT)
  ENDIF
!
  IF(NSEA_WIND_ID/=NUNDEF)THEN
    YCOMMENT='10m wind speed over sea'
    ZWRITE(:,1) = PSEA_WIND(:)
    CALL OASIS_PUT(NSEA_WIND_ID,KDATE,ZWRITE(:,:),IERR)
    CALL CHECK_SFX_SEND(KLUOUT,IERR,YCOMMENT)
  ENDIF
!
  IF(NSEA_FWSM_ID/=NUNDEF)THEN
    YCOMMENT='wind stress over sea'
    ZWRITE(:,1) = PSEA_FWSM(:)
    CALL OASIS_PUT(NSEA_FWSM_ID,KDATE,ZWRITE(:,:),IERR)
    CALL CHECK_SFX_SEND(KLUOUT,IERR,YCOMMENT)
  ENDIF
!
  IF(NSEA_EVAP_ID/=NUNDEF)THEN
    YCOMMENT='Evaporation over sea'
    ZWRITE(:,1) = PSEA_EVAP(:)
    CALL OASIS_PUT(NSEA_EVAP_ID,KDATE,ZWRITE(:,:),IERR)
    CALL CHECK_SFX_SEND(KLUOUT,IERR,YCOMMENT)
  ENDIF
!
  IF(NSEA_RAIN_ID/=NUNDEF)THEN
    YCOMMENT='Rainfall rate over sea'
    ZWRITE(:,1) = PSEA_RAIN(:)
    CALL OASIS_PUT(NSEA_RAIN_ID,KDATE,ZWRITE(:,:),IERR)
    CALL CHECK_SFX_SEND(KLUOUT,IERR,YCOMMENT)
  ENDIF
!
  IF(NSEA_SNOW_ID/=NUNDEF)THEN
    YCOMMENT='Snowfall rate over sea'
    ZWRITE(:,1) = PSEA_SNOW(:)
    CALL OASIS_PUT(NSEA_SNOW_ID,KDATE,ZWRITE(:,:),IERR)
    CALL CHECK_SFX_SEND(KLUOUT,IERR,YCOMMENT)
  ENDIF
!
  IF(NSEA_WATF_ID/=NUNDEF)THEN
    YCOMMENT='Freshwater flux over sea'
    ZWRITE(:,1) = PSEA_WATF(:)
    CALL OASIS_PUT(NSEA_WATF_ID,KDATE,ZWRITE(:,:),IERR)
    CALL CHECK_SFX_SEND(KLUOUT,IERR,YCOMMENT)
  ENDIF

! * Sea-ice output fields
!
  IF(LCPL_SEAICE)THEN
!
    IF(NSEAICE_HEAT_ID/=NUNDEF)THEN
      YCOMMENT='Sea-ice non solar net heat flux over sea-ice'
      ZWRITE(:,1) = PSEAICE_HEAT(:)
      CALL OASIS_PUT(NSEAICE_HEAT_ID,KDATE,ZWRITE(:,:),IERR)
      CALL CHECK_SFX_SEND(KLUOUT,IERR,YCOMMENT)
    ENDIF
!
    IF(NSEAICE_SNET_ID/=NUNDEF)THEN
      YCOMMENT='Sea-ice solar net heat flux over sea-ice'
      ZWRITE(:,1) = PSEAICE_SNET(:)
      CALL OASIS_PUT(NSEAICE_SNET_ID,KDATE,ZWRITE(:,:),IERR)
      CALL CHECK_SFX_SEND(KLUOUT,IERR,YCOMMENT)
    ENDIF
!
    IF(NSEAICE_EVAP_ID/=NUNDEF)THEN
      YCOMMENT='Sea-ice sublimation over sea-ice'
      ZWRITE(:,1) = PSEAICE_EVAP(:)
      CALL OASIS_PUT(NSEAICE_EVAP_ID,KDATE,ZWRITE(:,:),IERR)
      CALL CHECK_SFX_SEND(KLUOUT,IERR,YCOMMENT)
    ENDIF
!
  ENDIF
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('SFX_OASIS_SEND',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!
SUBROUTINE CHECK_SFX_SEND(KLUOUT,KERR,HCOMMENT)
!
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
INTEGER,          INTENT(IN) :: KLUOUT
INTEGER,          INTENT(IN) :: KERR
CHARACTER(LEN=*), INTENT(IN) :: HCOMMENT
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('SFX_OASIS_SEND:CHECK_SFX_SEND',0,ZHOOK_HANDLE)
!
IF (KERR/=OASIS_OK.AND.KERR<OASIS_SENT) THEN
   WRITE(KLUOUT,'(A,I4)')'Return OASIS code from sending '//TRIM(HCOMMENT)//' : ',KERR
   CALL ABOR1_SFX('SFX_OASIS_SEND: problem sending '//TRIM(HCOMMENT))
ENDIF 
!
IF (LHOOK) CALL DR_HOOK('SFX_OASIS_SEND:CHECK_SFX_SEND',1,ZHOOK_HANDLE)
!
END SUBROUTINE CHECK_SFX_SEND
!
!-------------------------------------------------------------------------------
#endif
!-------------------------------------------------------------------------------
!
END SUBROUTINE SFX_OASIS_SEND
