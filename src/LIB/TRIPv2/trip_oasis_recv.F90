!#########
SUBROUTINE TRIP_OASIS_RECV(KLISTING,KLON,KLAT,PTIMEC,PRUNOFF,  & 
                           PDRAIN,PCALVING,PRECHARGE,PSRC_FLOOD)
!#############################################################################
!
!!****  *TRIP_OASIS_RECV* - Receive coupling fields
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
USE MODD_TRIP, ONLY : TP => TRIP
!
USE MODD_TRIP_PAR,   ONLY : XUNDEF
!
USE MODN_TRIP_OASIS, ONLY : XTSTEP_CPL_LAND
USE MODD_TRIP_OASIS
!
USE MODD_TRIP_GRID, ONLY : TPG => TRIP_GRID
!
USE MODI_FLOOD_REDISTRIB
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
#ifdef TRIPOASIS
USE MOD_OASIS
#endif
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER, INTENT(IN)               :: KLISTING
INTEGER, INTENT(IN)               :: KLON
INTEGER, INTENT(IN)               :: KLAT
REAL,    INTENT(IN)               :: PTIMEC        ! Cumulated run time step (s)
!
REAL, DIMENSION(:,:), INTENT(OUT) :: PRUNOFF       ! Surface runoff                  (kg/s)
REAL, DIMENSION(:,:), INTENT(OUT) :: PDRAIN        ! Deep drainage                   (kg/s)
REAL, DIMENSION(:,:), INTENT(OUT) :: PCALVING      ! Calving flux                    (kg/s)
REAL, DIMENSION(:,:), INTENT(OUT) :: PRECHARGE     ! Groundwater recharge            (kg/s)
REAL, DIMENSION(:,:), INTENT(OUT) :: PSRC_FLOOD    ! Input P-E-I flood source term   (kg/s)
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL, DIMENSION(KLON,KLAT) :: ZREAD
REAL, DIMENSION(KLON,KLAT) :: ZPFLOOD
REAL, DIMENSION(KLON,KLAT) :: ZEFLOOD
REAL, DIMENSION(KLON,KLAT) :: ZIFLOOD
REAL, DIMENSION(KLON,KLAT) :: ZSRC_FLOOD
REAL, DIMENSION(KLON,KLAT) :: ZRESIDU
!
CHARACTER(LEN=50)          :: YCOMMENT
INTEGER                    :: IDATE  ! current coupling time step (s)
INTEGER                    :: IERR   ! Error info
INTEGER                    :: IERR1  ! Error info
INTEGER                    :: IERR2  ! Error info
INTEGER                    :: IERR3  ! Error info
LOGICAL                    :: LFLDOK
INTEGER                    :: JVAR
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
#ifdef TRIPOASIS
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TRIP_OASIS_RECV',0,ZHOOK_HANDLE)
!
!*       1.     Initialize :
!               ------------
!
IDATE = INT(PTIMEC)
!
PDRAIN    (:,:) = 0.0
PRUNOFF   (:,:) = 0.0
PCALVING  (:,:) = 0.0
PRECHARGE (:,:) = 0.0
PSRC_FLOOD(:,:) = 0.0
!
!-------------------------------------------------------------------------------
!
!*       2.     Get coupling fields :
!               ---------------------
!
IF(LCPL_LAND)THEN
!
! * Receive river input fields
!
  ZREAD(:,:) = XUNDEF  
  YCOMMENT='Surface runoff'
  CALL OASIS_GET(NRUNOFF_ID,IDATE,ZREAD(:,:),IERR)
  CALL CHECK_TRIP_RECV(IERR,YCOMMENT)
  CALL KGM2_TO_KGS(IERR,ZREAD,PRUNOFF)
!
  ZREAD(:,:) = XUNDEF  
  YCOMMENT='Deep drainage'
  CALL OASIS_GET(NDRAIN_ID,IDATE,ZREAD(:,:),IERR)
  CALL CHECK_TRIP_RECV(IERR,YCOMMENT)
  CALL KGM2_TO_KGS(IERR,ZREAD,PDRAIN)
!
  IF(LCPL_CALVING)THEN
    ZREAD(:,:) = XUNDEF  
    YCOMMENT='calving flux'
    CALL OASIS_GET(NCALVING_ID,IDATE,ZREAD(:,:),IERR)
    CALL CHECK_TRIP_RECV(IERR,YCOMMENT)
    CALL KGM2_TO_KGS(IERR,ZREAD,PCALVING)
  ENDIF
!
  IF(LCPL_GW)THEN
    ZREAD(:,:) = XUNDEF  
    YCOMMENT='groundwater recharge'
    CALL OASIS_GET(NRECHARGE_ID,IDATE,ZREAD(:,:),IERR)
    CALL CHECK_TRIP_RECV(IERR,YCOMMENT)    
    CALL KGM2_TO_KGS(IERR,ZREAD,PRECHARGE)
  ENDIF 
!
  IF(LCPL_FLOOD)THEN
!
    ZPFLOOD(:,:)=XUNDEF
    ZEFLOOD(:,:)=XUNDEF
    ZIFLOOD(:,:)=XUNDEF
    ZRESIDU(:,:)=XUNDEF
!          
    YCOMMENT='flood precip interception'
    CALL OASIS_GET(NPFLOOD_ID,IDATE,ZPFLOOD(:,:),IERR1)
    CALL CHECK_TRIP_RECV(IERR1,YCOMMENT)
!
    YCOMMENT='flood evaporation'
    CALL OASIS_GET(NEFLOOD_ID,IDATE,ZEFLOOD(:,:),IERR2)
    CALL CHECK_TRIP_RECV(IERR2,YCOMMENT)
!
    YCOMMENT='flood infiltration'
    CALL OASIS_GET(NIFLOOD_ID,IDATE,ZIFLOOD(:,:),IERR3)
    CALL CHECK_TRIP_RECV(IERR3,YCOMMENT)
!
    LFLDOK=(IERR1>=OASIS_RECVD.AND.IERR1==IERR2.AND.IERR1==IERR3)
!
    IF(LFLDOK)THEN     
       IERR = IERR1
       CALL FLOOD_REDISTRIB(TP, TPG, &
                            KLON,KLAT,ZPFLOOD,ZEFLOOD,ZIFLOOD,ZREAD)
       ZSRC_FLOOD(:,:) = ZPFLOOD(:,:)-ZEFLOOD(:,:)-ZIFLOOD(:,:)
       CALL KGM2_TO_KGS(IERR,ZSRC_FLOOD,PSRC_FLOOD)
       WHERE(TPG%GMASK(:,:)) 
         PRUNOFF(:,:)=PRUNOFF(:,:)+ZRESIDU(:,:)
       ENDWHERE
    ENDIF
!
  ENDIF
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('TRIP_OASIS_RECV',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!
SUBROUTINE CHECK_TRIP_RECV(KERR,HCOMMENT)
!
USE MODI_ABORT_TRIP
!
IMPLICIT NONE
!
INTEGER,          INTENT(IN) :: KERR
CHARACTER(LEN=*), INTENT(IN) :: HCOMMENT
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('TRIP_OASIS_RECV:CHECK_TRIP_RECV',0,ZHOOK_HANDLE)
!
! Check receiving field 
!
IF (KERR/=OASIS_OK.AND.KERR<OASIS_RECVD) THEN
   WRITE(KLISTING,'(A,I4)')'Return code from receiving '//TRIM(HCOMMENT)//' : ',KERR
   CALL ABORT_TRIP('TRIP_OASIS_RECV: problem receiving '//TRIM(HCOMMENT))
ENDIF
!
IF (LHOOK) CALL DR_HOOK('TRIP_OASIS_RECV:CHECK_TRIP_RECV',1,ZHOOK_HANDLE)
!
END SUBROUTINE CHECK_TRIP_RECV
!
!-------------------------------------------------------------------------------
!
SUBROUTINE KGM2_TO_KGS(KERR,PIN,POUT)
!
IMPLICIT NONE
!
INTEGER,              INTENT(IN ) :: KERR
REAL, DIMENSION(:,:), INTENT(IN ) :: PIN
REAL, DIMENSION(:,:), INTENT(OUT) :: POUT
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('TRIP_OASIS_RECV:KGM2_TO_KGS',0,ZHOOK_HANDLE)
!
! kg/m2 -> kg/s
!
IF(KERR>=OASIS_RECVD)THEN
  WHERE(TPG%GMASK(:,:)) 
        POUT(:,:) = PIN(:,:) * TPG%XAREA(:,:) / XTSTEP_CPL_LAND
  ELSEWHERE
        POUT(:,:) = XUNDEF
  ENDWHERE
ELSE
  POUT(:,:) = 0.0
ENDIF
!
IF (LHOOK) CALL DR_HOOK('TRIP_OASIS_RECV:KGM2_TO_KGS',1,ZHOOK_HANDLE)
!
END SUBROUTINE KGM2_TO_KGS
!
!-------------------------------------------------------------------------------
#endif
!-------------------------------------------------------------------------------
!
END SUBROUTINE TRIP_OASIS_RECV
