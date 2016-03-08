!     #########
SUBROUTINE ISBA_BUDGET (IO, IR, DGEIP, OWATER_BUDGET, HSNOW_ISBA, PTSTEP,&
                       PDG, PDZG, PWG_INI, PWGI_INI, PWR_INI, PSWE_INI, PRAIN, PSNOW, PEVAP  )
!     ###############################################################################
!
!!****  *ISBA_BUDGET * - water and energy budget for ISBA
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
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
!!      Original    07/2012
!!
!!------------------------------------------------------------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_CSTS,       ONLY : XRHOLW
!     
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEIP
!
LOGICAL, INTENT(IN) :: OWATER_BUDGET
!
 CHARACTER(LEN=*),     INTENT(IN)  :: HSNOW_ISBA ! 'DEF' = Default F-R snow scheme
!                                               !         (Douville et al. 1995)
!                                               ! '3-L' = 3-L snow scheme (option)
!                                               !         (Boone and Etchevers 2000)
!                                               ! 'CRO' = Crocus snow scheme
REAL,                  INTENT(IN) :: PTSTEP     ! timestep of the integration    (s)
!
REAL, DIMENSION(:,:),  INTENT(IN) :: PDG        ! soil layer depth               (m)
REAL, DIMENSION(:,:),  INTENT(IN) :: PDZG       ! soil layer thickness           (m)
!
REAL, DIMENSION(:),    INTENT(IN) :: PWG_INI    ! total wg at t-1                (kg m-2)
REAL, DIMENSION(:),    INTENT(IN) :: PWGI_INI   ! total wgi at t-1               (kg m-2)
REAL, DIMENSION(:),    INTENT(IN) :: PWR_INI    ! total wr at t-1                (kg m-2)
REAL, DIMENSION(:),    INTENT(IN) :: PSWE_INI   ! total swe at t-1               (kg m-2)
!
REAL, DIMENSION(:),    INTENT(IN)  :: PRAIN     ! Rain rate                      (kg/m2/s)
REAL, DIMENSION(:),    INTENT(IN)  :: PSNOW     ! Snow rate                      (kg/m2/s)
REAL, DIMENSION(:),    INTENT(IN)  :: PEVAP     ! total evaporative flux         (kg/m2/s)
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(IR%XWR,1)) :: ZINPUT
REAL, DIMENSION(SIZE(IR%XWR,1)) :: ZOUTPUT
REAL, DIMENSION(SIZE(IR%XWR,1)) :: ZTENDENCY
REAL, DIMENSION(SIZE(IR%XWR,1)) :: ZICEFLUX
REAL, DIMENSION(SIZE(IR%XWR,1)) :: ZSWE_T
REAL, DIMENSION(SIZE(IR%XWR,1)) :: ZWG_T
REAL, DIMENSION(SIZE(IR%XWR,1)) :: ZWGI_T
REAL, DIMENSION(SIZE(IR%XWR,1)) :: ZSNDRIFT
!
INTEGER :: INI, INL, INLS
INTEGER :: JI, JL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ISBA_BUDGET',0,ZHOOK_HANDLE)
!
!*      1.0    Init
!       -----------
!
INI =SIZE(IR%XWG,1)
INL =SIZE(IR%XWG,2)
INLS=SIZE(IR%TSNOW%WSNOW,2)
!
DGEIP%XDWG (:) = XUNDEF
DGEIP%XDWGI(:) = XUNDEF
DGEIP%XDWR (:) = XUNDEF
DGEIP%XDSWE(:) = XUNDEF
!
IF (HSNOW_ISBA=='3-L'.OR.HSNOW_ISBA=='CRO') THEN
  ZSNDRIFT(:) = DGEIP%XSNDRIFT(:)
ELSE
  ZSNDRIFT(:) = 0.0
ENDIF
!
!*      2.0    Comptut isba water budget in kg/m2/s
!       -------------------------------------------
!
IF(OWATER_BUDGET)THEN
!
! total swe at t in kg/m2
  ZSWE_T(:)=0.0
  DO JL=1,INLS
    DO JI=1,INI
      ZSWE_T(JI) = ZSWE_T(JI)+IR%TSNOW%WSNOW(JI,JL,1)
    ENDDO
  ENDDO
!
! total wg and wgi at t in kg/m2
  ZWG_T (:)= 0.0
  ZWGI_T(:)= 0.0
  IF(IO%CISBA=='DIF')THEN
    DO JL=1,INL
      DO JI=1,INI
        IF(IR%XWG(JI,JL,1)/=XUNDEF)THEN
          ZWG_T (JI) = ZWG_T (JI)+IR%XWG(JI,JL,1) *PDZG(JI,JL)*XRHOLW
          ZWGI_T(JI) = ZWGI_T(JI)+IR%XWGI(JI,JL,1)*PDZG(JI,JL)*XRHOLW
        ENDIF
      ENDDO
    ENDDO
  ELSE
    ZWG_T (:) = IR%XWG (:,2,1)*PDG(:,2)*XRHOLW
    ZWGI_T(:) = IR%XWGI(:,2,1)*PDG(:,2)*XRHOLW
    IF(IO%CISBA=='3-L')THEN
      ZWG_T(:)=ZWG_T(:)+IR%XWG(:,3,1)*(PDG(:,3)-PDG(:,2))*XRHOLW
    ENDIF
  ENDIF
!
! Comptut reservoir time tendencies in kg/m2/s
  DGEIP%XDWG (:) = (ZWG_T   (:)-PWG_INI (:))/PTSTEP
  DGEIP%XDWGI(:) = (ZWGI_T  (:)-PWGI_INI(:))/PTSTEP
  DGEIP%XDWR (:) = (IR%XWR(:,1)-PWR_INI (:))/PTSTEP
  DGEIP%XDSWE(:) = (ZSWE_T  (:)-PSWE_INI(:))/PTSTEP
!
! ice calving flux if used
  IF(IO%LGLACIER)THEN
    ZICEFLUX(:)=DGEIP%XICEFLUX(:)
  ELSE
    ZICEFLUX(:)=0.0
  ENDIF
!
! total input water in the system at t
  ZINPUT(:)=PRAIN(:)+PSNOW(:)+DGEIP%XIFLOOD(:)+DGEIP%XIRRIG_FLUX(:)
!
! total output water in the system at t
  ZOUTPUT(:) = PEVAP  (:)+DGEIP%XDRAIN  (:)+DGEIP%XRUNOFF (:) &
             + DGEIP%XPFLOOD(:)+ZICEFLUX(:)+ZSNDRIFT(:)
!
! total reservoir time tendencies at "t - (t-1)"
  ZTENDENCY(:) = DGEIP%XDWG(:)+DGEIP%XDWGI(:)+DGEIP%XDWR(:)+DGEIP%XDSWE(:)
!
! isba water budget (dw/dt=in-out) in kg/m2/s
  DGEIP%XWATBUD(:)=ZTENDENCY(:)-(ZINPUT(:)-ZOUTPUT(:))
!
ENDIF
!
!*      3.0    Comptut isba energy budget in W/m2
!       -----------------------------------------
!
! not yet implemented
!
IF (LHOOK) CALL DR_HOOK('ISBA_BUDGET',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------------
!
END SUBROUTINE ISBA_BUDGET
