!#########
SUBROUTINE TRIP_OASIS_SEND(KLISTING,KLON,KLAT,PTIMEC)
!############################################
!
!!****  *TRIP_OASIS_SEND* - Send coupling fields
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
USE MODD_TRIP_PAR,   ONLY : XUNDEF
!
USE MODN_TRIP_OASIS, ONLY : XTSTEP_CPL_SEA, XTSTEP_CPL_LAND
USE MODD_TRIP_OASIS
!
USE MODD_TRIP, ONLY : TP => TRIP
USE MODD_TRIP_GRID, ONLY : TPG => TRIP_GRID
!
USE MODI_ABORT_TRIP
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
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER                    :: IDATE   ! current coupling time step (s)
INTEGER                    :: IERR    ! Error info
INTEGER                    :: JVAR
CHARACTER(LEN=50)          :: YCOMMENT
!
REAL,    DIMENSION(KLON,KLAT) :: ZWRITE
LOGICAL, DIMENSION(KLON,KLAT) :: LMASK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
#ifdef TRIPOASIS
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TRIP_OASIS_SEND',0,ZHOOK_HANDLE)
!
!*       1.     Define current coupling time step in second :
!               ---------------------------------------------
!
IDATE = INT(PTIMEC)
!
!-------------------------------------------------------------------------------
!
!*       2.     Send coupling fields to land surface model :
!               -------------------------------------------
!
!
IF(LCPL_LAND.AND.MOD(PTIMEC,XTSTEP_CPL_LAND)==0.0)THEN
!
  IF(LCPL_GW)THEN
!          
    YCOMMENT='Water table depth'
    CALL OASIS_PUT(NWTD_ID,IDATE,TP%XCPL_WTD(:,:),IERR)
    CALL CHECK_TRIP_SEND(YCOMMENT)
!          
    YCOMMENT='Grid-cell fraction of WTD to rise'
    CALL OASIS_PUT(NFWTD_ID,IDATE,TP%XCPL_FWTD(:,:),IERR)
    CALL CHECK_TRIP_SEND(YCOMMENT)
!
  ENDIF
!
  IF(LCPL_FLOOD)THEN
!
    LMASK(:,:) = TPG%GMASK(:,:)
!
    YCOMMENT='Flood fraction'
    CALL MASK_TRIP(TP%XCPL_FFLOOD(:,:),ZWRITE(:,:),LMASK(:,:))
    CALL OASIS_PUT(NFFLOOD_ID,IDATE,ZWRITE(:,:),IERR)
    CALL CHECK_TRIP_SEND(YCOMMENT)
!          
    YCOMMENT='Flood potential infiltration'
    CALL MASK_TRIP(TP%XCPL_PIFLOOD(:,:),ZWRITE(:,:),LMASK(:,:))
    CALL OASIS_PUT(NPIFLOOD_ID,IDATE,ZWRITE(:,:),IERR)
    CALL CHECK_TRIP_SEND(YCOMMENT)
!
  ENDIF  
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       3.     Send coupling fields to ocean :
!               -------------------------------
!
!
IF(LCPL_SEA.AND.MOD(PTIMEC,XTSTEP_CPL_SEA)==0.0)THEN
!
! * Send sea output fields
!
  LMASK(:,:) = (TPG%NGRCN(:,:)==9.OR.TPG%NGRCN(:,:)==12)
!
  YCOMMENT='Discharge to ocean'
  CALL MASK_TRIP(TP%XCPL_RIVDIS(:,:),ZWRITE(:,:),LMASK(:,:))
  CALL OASIS_PUT(NRIVDIS_ID,IDATE,ZWRITE(:,:),IERR)
  CALL CHECK_TRIP_SEND(YCOMMENT)
!
! * Calving output fields
!
  IF(LCPL_CALVSEA)THEN
!
    LMASK(:,:) = TPG%GMASK_GRE(:,:)
!
    YCOMMENT='Calving flux over greenland'
    CALL MASK_TRIP(TP%XCPL_CALVGRE(:,:),ZWRITE(:,:),LMASK(:,:))
    CALL OASIS_PUT(NCALVGRE_ID,IDATE,ZWRITE(:,:),IERR)
    CALL CHECK_TRIP_SEND(YCOMMENT)
!
    LMASK(:,:) = TPG%GMASK_ANT(:,:)
!
    YCOMMENT='Calving flux over antarctica'
    CALL MASK_TRIP(TP%XCPL_CALVANT(:,:),ZWRITE(:,:),LMASK(:,:))
    CALL OASIS_PUT(NCALVANT_ID,IDATE,ZWRITE(:,:),IERR)
    CALL CHECK_TRIP_SEND(YCOMMENT)
!
  ENDIF
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('TRIP_OASIS_SEND',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!
SUBROUTINE CHECK_TRIP_SEND(HCOMMENT)
!
USE MODI_ABORT_TRIP
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(IN) :: HCOMMENT
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('TRIP_OASIS_SEND:CHECK_TRIP_SEND',0,ZHOOK_HANDLE)
!
! Check receiving field 
!
IF (IERR/=OASIS_OK.AND.IERR<OASIS_SENT) THEN
   WRITE(KLISTING,'(A,I4)')'Return code from sending '//TRIM(HCOMMENT)//' : ',IERR
   CALL ABORT_TRIP('TRIP_OASIS_SEND: problem sending '//TRIM(HCOMMENT))
ENDIF
!
IF (LHOOK) CALL DR_HOOK('TRIP_OASIS_SEND:CHECK_TRIP_SEND',1,ZHOOK_HANDLE)
!
END SUBROUTINE CHECK_TRIP_SEND
!
!-------------------------------------------------------------------------------
!
SUBROUTINE MASK_TRIP(PIN,POUT,OMASK)
!
IMPLICIT NONE
!
REAL,    DIMENSION(:,:), INTENT(INOUT) :: PIN
REAL,    DIMENSION(:,:), INTENT(OUT  ) :: POUT
LOGICAL, DIMENSION(:,:), INTENT(IN   ) :: OMASK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('TRIP_OASIS_SEND:MASK_TRIP',0,ZHOOK_HANDLE)
!
WHERE(OMASK(:,:)) 
      POUT(:,:) = PIN(:,:)
      PIN (:,:) = 0.0     
ELSEWHERE
      POUT(:,:) = XUNDEF
      PIN (:,:) = XUNDEF    
ENDWHERE
!
IF (LHOOK) CALL DR_HOOK('TRIP_OASIS_SEND:MASK_TRIP',1,ZHOOK_HANDLE)
!
END SUBROUTINE MASK_TRIP
!
!
!-------------------------------------------------------------------------------
#endif
!-------------------------------------------------------------------------------
!
END SUBROUTINE TRIP_OASIS_SEND
