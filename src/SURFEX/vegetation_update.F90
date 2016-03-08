!     #########
    SUBROUTINE VEGETATION_UPDATE (DTCO, DTI, KDIM, IO, IMT, IMM, IMI, IMA, &
                                  PTSTEP,TTIME,PCOVER,OCOVER,       &
                                  OAGRIP, HSFTYPE, OALB, MSS,       &
                                  ODUPDATED, OABSENT              )
!   ###############################################################
!!****  *VEGETATION EVOL*
!!
!!    PURPOSE
!!    -------
!
!     performs the time evolution of vegetation parameters
!       at UTC midnight for prescribed parameters, with effective change each ten days
!              
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      
!!    none
!!
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!      V. Masson          * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/03/03 
!!
!!      P Le Moigne 09/2005 AGS modifs of L. Jarlan
!!      P Samuelsson 10/2014 MEB
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!
!
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_TIME_t, ISBA_PARAM_MEB_t, &
                              ISBA_PARAM_ALB_t, ISBA_PARAM_IRRIG_t
!
USE MODD_SSO_n, ONLY : SSO_t
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
!
USE MODD_TYPE_DATE_SURF
!
USE MODI_INIT_ISBA_MIXPAR
USE MODI_CONVERT_PATCH_ISBA
USE MODI_INIT_FROM_DATA_GRDN_n
USE MODI_INIT_FROM_DATA_GREENROOF_n
USE MODI_SUBSCALE_Z0EFF
USE MODI_ALBEDO
USE MODI_UPDATE_DATA_COVER
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTI
INTEGER, INTENT(IN) :: KDIM
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
!
TYPE(ISBA_PARAM_TIME_t), INTENT(INOUT) :: IMT
TYPE(ISBA_PARAM_MEB_t), INTENT(INOUT) :: IMM
TYPE(ISBA_PARAM_ALB_t), INTENT(INOUT) :: IMA
TYPE(ISBA_PARAM_IRRIG_t), INTENT(INOUT) :: IMI
!
REAL,                 INTENT(IN)    :: PTSTEP  ! time step
TYPE(DATE_TIME),      INTENT(IN)    :: TTIME   ! UTC time
REAL,   DIMENSION(:,:), INTENT(IN)  :: PCOVER  ! cover types
LOGICAL, DIMENSION(:), INTENT(IN)   :: OCOVER
LOGICAL,              INTENT(IN)    :: OAGRIP
CHARACTER(LEN=*),     INTENT(IN)    :: HSFTYPE ! nature / garden
!
LOGICAL, INTENT(IN) :: OALB
!
TYPE(SSO_t), INTENT(INOUT) :: MSS
!
LOGICAL,              INTENT(OUT)   :: ODUPDATED  ! T if parameters are being reset
LOGICAL,DIMENSION(:), INTENT(IN), OPTIONAL :: OABSENT ! T where field is not defined
!
!*      0.2    declarations of local variables
!
INTEGER                                  :: IDECADE, IDECADE2  ! decade of simulation
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-----------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('VEGETATION_UPDATE',0,ZHOOK_HANDLE)
!
!*      2.     Non-interactive vegetation
!              --------------------------
!
!*      2.1    Decade
!              ------
!
IDECADE = 3 * ( TTIME%TDATE%MONTH - 1 ) + MIN(TTIME%TDATE%DAY-1,29) / 10 + 1
IDECADE2 = IDECADE
ODUPDATED=.FALSE.
!
!*      2.2    From ecoclimap
!              --------------
!
!* new decade?
  IF ( MOD(MIN(TTIME%TDATE%DAY,30),10)==1 .AND. TTIME%TIME - PTSTEP < 0.) THEN
    ODUPDATED=.TRUE.
!* time varying parameters
    IF (IO%LECOCLIMAP .OR. HSFTYPE=='NAT') THEN
!* new year ? --> recomputes data LAI and derivated parameters (usefull in case of ecoclimap2)
      CALL UPDATE_DATA_COVER(DTCO, DTI, KDIM, IO%NPATCH, IO%LMEB_PATCH, TTIME%TDATE%YEAR)  
      IF (HSFTYPE=='NAT') THEN
        CALL INIT_ISBA_MIXPAR(DTCO, DTI, KDIM, IO, IDECADE,IDECADE2,PCOVER,OCOVER,HSFTYPE)
        CALL CONVERT_PATCH_ISBA(DTCO, DTI, IO, IDECADE, IDECADE2, PCOVER, OCOVER,&
                              OAGRIP, HSFTYPE, OALB, IMT=IMT, IMM=IMM, IMI=IMI   )
      ELSE
        CALL CONVERT_PATCH_ISBA(DTCO, DTI, IO, IDECADE, IDECADE2, PCOVER, OCOVER,&
                              OAGRIP, HSFTYPE, OALB, IMT=IMT   )
      ENDIF
      IF ( IO%CALBEDO=='CM13') THEN
        CALL CONVERT_PATCH_ISBA(DTCO, DTI, IO, IDECADE, IDECADE2, PCOVER, OCOVER,&
                                OAGRIP,HSFTYPE, .FALSE., IMA=IMA )
      ENDIF
    ELSEIF (.NOT.OALB) THEN

      IF (HSFTYPE=='GRD') THEN
        CALL INIT_FROM_DATA_GRDN_n(DTI, IDECADE, IO%CPHOTO, .TRUE.,  VMT=IMT  )  
      ELSEIF (HSFTYPE=='GNR') THEN
        CALL INIT_FROM_DATA_GREENROOF_n(DTI,  IDECADE,IO%CPHOTO, .TRUE.,  VMT=IMT )
      ENDIF

    ENDIF
!
!* default values to avoid problems in physical routines
!  for points where there is no vegetation or soil to be simulated by ISBA.
    IF (PRESENT(OABSENT) .AND. .NOT.OALB) THEN
        WHERE (OABSENT(:))
          IMT%XVEG       (:,1) = 0.
          IMT%XLAI       (:,1) = 0.
          IMT%XRSMIN     (:,1) = 40.
          IMT%XGAMMA     (:,1) = 0.
          IMT%XWRMAX_CF  (:,1) = 0.2
          IMT%XRGL       (:,1) = 100.
          IMT%XCV        (:,1) = 2.E-5
          IMT%XZ0        (:,1) = 0.013
          IMT%XALBNIR_VEG(:,1) = 0.30
          IMT%XALBVIS_VEG(:,1) = 0.30
          IMT%XALBUV_VEG (:,1) = 0.06
          IMT%XEMIS      (:,1) = 0.94                
        END WHERE
        IF (IO%CPHOTO/='NON') THEN
          WHERE (OABSENT(:))
            IMT%XGMES      (:,1) = 0.020
            IMT%XBSLAI     (:,1) = 0.36
            IMT%XLAIMIN    (:,1) = 0.3
            IMT%XSEFOLD    (:,1) = 90*86400.
            IMT%XGC        (:,1) = 0.00025                  
          END WHERE
          IF (IO%CPHOTO/='AGS' .AND. IO%CPHOTO/='LAI') THEN
            WHERE (OABSENT(:)) IMT%XF2I       (:,1) = 0.3
            IF (IO%CPHOTO=='NIT' .OR. IO%CPHOTO=='NCB') THEN
              WHERE (OABSENT(:))
                IMT%XCE_NITRO  (:,1) = 7.68
                IMT%XCF_NITRO  (:,1) = -4.33
                IMT%XCNA_NITRO (:,1) = 1.3                      
              END WHERE
            ENDIF
          ENDIF
        ENDIF
    ENDIF

    IF (HSFTYPE=='NAT') THEN
!* albedo
      CALL ALBEDO(IO%CALBEDO, IMT,IMA )  
!
!* effective roughness length
      IF (.NOT.OALB) CALL SUBSCALE_Z0EFF(MSS,IMT%XZ0,.FALSE.  )  
    ENDIF

  END IF
IF (LHOOK) CALL DR_HOOK('VEGETATION_UPDATE',1,ZHOOK_HANDLE)
!
!*      2.3    Prescribed vegetation
!              ---------------------
!
!-----------------------------------------------------------------
!
END SUBROUTINE VEGETATION_UPDATE
