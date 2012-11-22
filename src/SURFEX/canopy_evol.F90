!     #########################################
      SUBROUTINE CANOPY_EVOL(KI,KLVL,PTSTEP, KIMPL,                                     &
                              PZZ,PWIND,PTA,PQA,PPA,PRHOA,                              &
                              PSFLUX_U,PSFLUX_T,PSFLUX_Q,                               &
                              PFORC_U,PDFORC_UDU,PFORC_E,PDFORC_EDE,PFORC_T,PDFORC_TDT, &
                              PFORC_Q,PDFORC_QDQ,                                       &
                              PZ,PZF,PDZ,PDZF,PU,PTKE,PT,PQ,PLMO,PLM,PLEPS,PP,PUSTAR,   &
                              PALFAU,PBETAU,PALFATH,PBETATH,PALFAQ,PBETAQ, ONEUTRAL     ) 
!     #########################################
!
!!****  *CANOPY_EVOL* - evolution of canopy
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
!!	V. Masson   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/2006 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CSTS,        ONLY : XG, XRD, XCPD, XP00
USE MODD_SURF_PAR,    ONLY : XUNDEF
USE MODD_CANOPY_TURB, ONLY : XCMFS, XCSHF
!
USE MODI_RMC01_SURF
USE MODI_CANOPY_EVOL_WIND
USE MODI_CANOPY_EVOL_TKE
USE MODI_CANOPY_EVOL_TEMP
!
USE MODE_SBLS
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER,                  INTENT(IN)    :: KI        ! number of horizontal points
INTEGER,                  INTENT(IN)    :: KLVL      ! number of levels in canopy
REAL,                     INTENT(IN)    :: PTSTEP    ! atmospheric time-step                 (s)
INTEGER,                  INTENT(IN)    :: KIMPL     ! implicitation code: 
!                                                    ! 1 : computes only alfa and beta coupling
!                                                    !     coefficients for all variables
!                                                    ! 2 : computes temporal evolution of the
!                                                    !     variables
REAL, DIMENSION(:,:), INTENT(IN)    :: PZZ          ! Mixing length generic profile at mid levels (-)
REAL, DIMENSION(:),      INTENT(IN)    :: PWIND     ! wind speed                            (m/s)
REAL, DIMENSION(:),      INTENT(IN)    :: PTA       ! Air temperature                       (K)
REAL, DIMENSION(:),      INTENT(IN)    :: PQA       ! Air humidity                          (kg/m3)
REAL, DIMENSION(:),      INTENT(IN)    :: PPA       ! Pressure at forcing level             (Pa)
REAL, DIMENSION(:),      INTENT(IN)    :: PRHOA     ! Air density at forcing level          (kg/m3)
REAL, DIMENSION(:),      INTENT(IN)    :: PSFLUX_U  ! surface flux u'w'                     (m2/s2)
REAL, DIMENSION(:),      INTENT(IN)    :: PSFLUX_T  ! surface flux w'T'                     (Km/s)
REAL, DIMENSION(:),      INTENT(IN)    :: PSFLUX_Q  ! surface flux w'q'                     (kg/m2/s)
REAL, DIMENSION(:,:), INTENT(IN)    :: PFORC_U      ! tendency of wind due to canopy drag   (m/s2)
REAL, DIMENSION(:,:), INTENT(IN)    :: PDFORC_UDU   ! formal derivative of the tendency of
!                                                   ! wind due to canopy drag               (1/s)
REAL, DIMENSION(:,:), INTENT(IN)    :: PFORC_E      ! tendency of TKE  due to canopy drag   (m2/s3)
REAL, DIMENSION(:,:), INTENT(IN)    :: PDFORC_EDE   ! formal derivative of the tendency of
!                                                   ! TKE  due to canopy drag               (1/s)
REAL, DIMENSION(:,:), INTENT(IN)    :: PFORC_T      ! tendency of Temp due to canopy drag   (T/s)
REAL, DIMENSION(:,:), INTENT(IN)    :: PDFORC_TDT   ! formal derivative of the tendency of
!                                                   ! Temp due to canopy drag               (1/s)
REAL, DIMENSION(:,:), INTENT(IN)    :: PFORC_Q      ! tendency of Hum. due to canopy drag   (kg/m3/s)
REAL, DIMENSION(:,:), INTENT(IN)    :: PDFORC_QDQ   ! formal derivative of the tendency of
!                                                   ! Hum. due to canopy drag               (1/s)
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PZ        ! heights of canopy levels              (m)
REAL, DIMENSION(:,:), INTENT(IN)    :: PZF       ! heights of bottom of canopy levels    (m)
REAL, DIMENSION(:,:), INTENT(IN)    :: PDZ       ! depth   of canopy levels              (m)
REAL, DIMENSION(:,:), INTENT(IN)    :: PDZF      ! depth between canopy levels           (m)
REAL, DIMENSION(:,:), INTENT(INOUT) :: PU        ! wind speed at canopy levels           (m/s)
REAL, DIMENSION(:,:), INTENT(INOUT) :: PTKE      ! TKE at canopy levels                  (m2/s2)
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT        ! Temperature at canopy levels          (K)
REAL, DIMENSION(:,:), INTENT(INOUT) :: PQ        ! Humidity    at canopy levels          (kg/m3)
REAL, DIMENSION(:,:), INTENT(INOUT) :: PLMO      ! Monin-Obhukov length                  (m)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PLM       ! mixing length                         (m)
REAL, DIMENSION(:,:), INTENT(OUT)   :: PLEPS     ! dissipative length                    (m)
REAL, DIMENSION(:,:), INTENT(INOUT) :: PP        ! Pressure    at canopy levels          (Pa)
REAL, DIMENSION(:),      INTENT(OUT)   :: PUSTAR    ! friction velocity at forcing level    (m/s)
!
REAL, DIMENSION(:),      INTENT(OUT)   :: PALFAU   ! V+(1) = alfa u'w'(1) + beta
REAL, DIMENSION(:),      INTENT(OUT)   :: PBETAU   ! V+(1) = alfa u'w'(1) + beta
REAL, DIMENSION(:),      INTENT(OUT)   :: PALFATH  ! Th+(1) = alfa w'th'(1) + beta
REAL, DIMENSION(:),      INTENT(OUT)   :: PBETATH  ! Th+(1) = alfa w'th'(1) + beta
REAL, DIMENSION(:),      INTENT(OUT)   :: PALFAQ   ! Q+(1) = alfa w'q'(1) + beta
REAL, DIMENSION(:),      INTENT(OUT)   :: PBETAQ   ! Q+(1) = alfa w'q'(1) + beta
!
LOGICAL, OPTIONAL, INTENT(IN) :: ONEUTRAL
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)) :: ZK                 ! mixing coefficient
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)) :: ZDKDDVDZ           ! formal derivative of mixing coefficient according to variable vertical gradient
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)) :: ZTH                ! potential temperature at full levels
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)) :: ZEXN               ! Exner function        at full levels
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)) :: ZUW                ! friction at mid levels
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)) :: ZFORC_TH           ! tendency of Temp due to canopy drag   (T/s)
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)) :: ZDFORC_THDTH       ! formal derivative of the tendency of
!                                                              ! Temp due to canopy drag               (1/s)
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)) :: ZWTH               ! w'Th' at mid levels
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)) :: ZWQ                ! w'q'  at mid levels
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)) :: ZSFTH              ! heat flux at atmospheric forcing level
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)) :: ZSFRV              ! vapor flux at atmospheric forcing level
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)) :: ZRHOA              ! air density profile
!
REAL, DIMENSION(SIZE(PZZ,1))   :: ZTHA               ! potential temperature at forcing level
REAL, DIMENSION(SIZE(PZZ,1))   :: ZSFLUX_TH          ! Surface flux w'Th' (mK/s)
!
REAL                 :: ZZ0                ! a value of z0 just for first time-step init.
!
LOGICAL :: GNEUTRAL
!
INTEGER :: JLAYER                              ! loop counter on layers
INTEGER :: JI                                  ! loop counter
!        
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('CANOPY_EVOL',0,ZHOOK_HANDLE)
!
GNEUTRAL = .FALSE. 
!
IF (PRESENT(ONEUTRAL)) GNEUTRAL = ONEUTRAL
!
!*    1. First time step initialization
!        ------------------------------
!
!* first time step (i.e. wind profile is initially zero) : 
!  set a neutral wind profile in relation with forcing wind
DO JI=1,KI
  IF(PWIND(JI)>0. .AND. PU(JI,KLVL-1)==0.) THEN
    ZZ0 = PZ(JI,1)/2.
    PU(JI,:) = PWIND(JI) * LOG (PZ(JI,:)/ZZ0) / LOG (PZ(JI,KLVL)/ZZ0)
  END IF
END DO
!-------------------------------------------------------------------------------
!
!*    5. mixing and dissipative lengths (at full levels)
!        ------------------------------
!
CALL RMC01_SURF(PZZ,PLMO,PLM,PLEPS,GNEUTRAL)
!
!-------------------------------------------------------------------------------
!
!*    6. time evolution of wind in canopy
!        --------------------------------
!
!*    6.1 mixing coefficient for wind at mid layers (at half levels)
!         ---------------------------
!
ZK = -999.
DO JLAYER=2,KLVL
  ZK(:,JLAYER) = 0.5 * XCMFS * PLM(:,JLAYER)   * SQRT(PTKE(:,JLAYER)  ) &
               + 0.5 * XCMFS * PLM(:,JLAYER-1) * SQRT(PTKE(:,JLAYER-1))  
END DO
!
!
!*    6.2 mixing coefficient vertical derivative at mid layers (at half levels)
!         --------------------------------------
!
!* there is no formal dependency on wind
ZDKDDVDZ = 0.
!
!
!*    6.3  time evolution of wind in canopy
!          --------------------------------
!
CALL CANOPY_EVOL_WIND(KI,KLVL,PTSTEP,KIMPL,PWIND,ZK,ZDKDDVDZ,PSFLUX_U,PFORC_U,PDFORC_UDU,PDZ,PDZF,PU,ZUW,PALFAU,PBETAU)
!
!*    6.4  Friction velocity at top of SBL layers
!          --------------------------------------
!
PUSTAR = SQRT(ABS(ZUW(:,KLVL)))
!
!-------------------------------------------------------------------------------
!
IF (GNEUTRAL) THEN
  !
  ZTH  = 300.  
  ZWTH = 0.
  ZWQ  = 0.
  !
ELSE
  !
  !*    7. time evolution of temperature in canopy
  !        ---------------------------------------
  !
  !*    7.3 convertion into potential temperature (at half levels)
  !         -------------------------------------
  !
  DO JLAYER=1,KLVL
    PP(:,JLAYER) = PPA(:) + XG * PRHOA(:) * (PZ(:,KLVL) - PZ(:,JLAYER))
  END DO
  ZEXN = (PP/XP00)**(XRD/XCPD)
  !
  ZTH  = XUNDEF
  WHERE(PT/=XUNDEF) ZTH(:,:) = PT(:,:) / ZEXN(:,:)
  !
  ZTHA(:) = PTA(:) / ZEXN(:,KLVL)
  !
  !*    7.1 mixing coefficient for wind at mid layers (at half levels)
  !         ---------------------------
  !
  ZK = -999.
  DO JLAYER=2,KLVL
    ZK(:,JLAYER) = 0.5 * XCSHF * PLM(:,JLAYER)   * SQRT(PTKE(:,JLAYER)  ) &
                 + 0.5 * XCSHF * PLM(:,JLAYER-1) * SQRT(PTKE(:,JLAYER-1))  
  END DO
  !
  !*    7.2 mixing coefficient vertical derivative at mid layers (at half levels)
  !         --------------------------------------
  !
  !* there is no formal dependency on temperature
  ZDKDDVDZ = 0.
  !
  !*    7.4  conversion of canopy tendency into potential temperature tendency
  !          -----------------------------------------------------------------
  !
  ZSFLUX_TH    = PSFLUX_T / ZEXN(:,1)
  ZFORC_TH     = PFORC_T  / ZEXN
  ZDFORC_THDTH = PDFORC_TDT
  !
  !
  !*    7.5  time evolution of temperature in canopy
  !          ---------------------------------------
  !
  CALL CANOPY_EVOL_TEMP(KI,KLVL,PTSTEP,KIMPL,ZTHA,ZK,ZDKDDVDZ,ZSFLUX_TH,ZFORC_TH,ZDFORC_THDTH,PDZ,PDZF,ZTH,ZWTH,PALFATH,PBETATH)
  !
  !*    7.6  convertion into absolute temperature
  !          ------------------------------------
  !
  WHERE(PT/=XUNDEF) PT(:,:) = ZTH(:,:) * ZEXN(:,:)
  !
  !-------------------------------------------------------------------------------
  !
  !*    8. time evolution of Humidity in canopy
  !        ------------------------------------
  !
  CALL CANOPY_EVOL_TEMP(KI,KLVL,PTSTEP,KIMPL,PQA,ZK,ZDKDDVDZ,PSFLUX_Q,PFORC_Q,PDFORC_QDQ,PDZ,PDZF,PQ,ZWQ,PALFAQ,PBETAQ)
  !
  !-------------------------------------------------------------------------------
  IF (KIMPL==1 .AND. LHOOK) CALL DR_HOOK('CANOPY_EVOL',1,ZHOOK_HANDLE)
  IF (KIMPL==1) RETURN
  !-------------------------------------------------------------------------------
  !
ENDIF
!
!*    9. time evolution of TKE in canopy
!        -------------------------------
!
CALL CANOPY_EVOL_TKE(KI,KLVL,PTSTEP,PRHOA,PZ,PZF,PDZ,PDZF,PFORC_E,PDFORC_EDE, &
                      PU,ZTH,ZUW,ZWTH,ZWQ,PLEPS,PTKE                          ) 
!
!-------------------------------------------------------------------------------
!
IF (.NOT.GNEUTRAL) THEN
  !
  !*   10. Monin-Obuhkov length
  !        --------------------
  !
  !* MO length is estimated using the heat and vapor turbulent fluxes at atmospheric level
  !  (it avoids the problems of vertical variation of the fluxes in the canopy)
  !  However, friction flux MUST be taken as the maximum flux on the
  !  profile, in order to avoid unrealistically small MO length when using
  !  small time-steps
  !
  ZRHOA(:,:) = SPREAD(PRHOA(:),2,KLVL)
  !
  ZSFTH(:,:)  = ZWTH(:,:)
  ZSFRV(:,:)  = ZWQ (:,:) / ZRHOA(:,:)
  !
  PLMO(:,:) = LMO(SQRT(ABS(ZUW)),PT,PQ,ZSFTH,ZSFRV)
  !
  DO JLAYER=1,KLVL
    WHERE (PLMO(:,JLAYER)>0.) PLMO(:,JLAYER) = MAX(PLMO(:,JLAYER),PZ(:,KLVL))
    WHERE (PLMO(:,JLAYER)<0.) PLMO(:,JLAYER) = MIN(PLMO(:,JLAYER),-PZ(:,KLVL))     
  ENDDO
  !
  !-------------------------------------------------------------------------------
  !
  !*   11. Security at atmospheric forcing level
  !        -------------------------------------
  !
  PT(:,KLVL) = PTA(:)
  !
  PQ(:,KLVL) = PQA(:)
  !
ENDIF
!
PU(:,KLVL) = PWIND(:)
!
IF (LHOOK) CALL DR_HOOK('CANOPY_EVOL',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE CANOPY_EVOL

