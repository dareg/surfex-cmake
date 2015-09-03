!     #########
      SUBROUTINE FLOOD_REDISTRIB (TP, TPG, &
                                  KLON,KLAT,PREAD,PSRC_FLOOD,PRESIDU)  
!     #####################################################################
!
!!****  *FLOOD_REDISTRIB*  
!!
!!    PURPOSE
!!    -------
!
!     
!!**  METHOD
!!    ------
!
!     Direct calculation
!
!!    EXTERNAL
!!    --------
!
!     None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!      
!!    REFERENCE
!!    ---------
!!      
!!    AUTHOR
!!    ------
!!      B. Decharme
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/12/13 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!
!
USE MODD_TRIP, ONLY : TRIP_t
USE MODD_TRIP_GRID, ONLY : TRIP_GRID_t
!
USE MODD_TRIP_PAR, ONLY : XUNDEF, XRHOLW
!
!
USE MODI_ABORT_TRIP
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(TRIP_t), INTENT(INOUT) :: TP
TYPE(TRIP_GRID_t), INTENT(INOUT) :: TPG
!
INTEGER, INTENT(IN)               :: KLON
INTEGER, INTENT(IN)               :: KLAT
!
REAL, DIMENSION(:,:), INTENT(IN ) :: PREAD      ![kg/m2/s]
REAL, DIMENSION(:,:), INTENT(OUT) :: PSRC_FLOOD ![kg/m2/s]
REAL, DIMENSION(:,:), INTENT(OUT) :: PRESIDU
!
!
!*      0.2    declarations of local variables
!
REAL,    DIMENSION(TPG%NBASMAX) :: ZBAS_IN
REAL,    DIMENSION(TPG%NBASMAX) :: ZBAS_OUT
REAL,    DIMENSION(TPG%NBASMAX) :: ZBAS_FACTOR
INTEGER, DIMENSION(TPG%NBASMAX) :: IBAS_NCELL
!
REAL,    DIMENSION(KLON,KLAT) :: ZOUT
INTEGER, DIMENSION(KLON,KLAT) :: ILOC_FLOOD_IN
INTEGER, DIMENSION(KLON,KLAT) :: ILOC_FLUXE_IN
LOGICAL, DIMENSION(KLON,KLAT) :: GCOLLOCATED
!
REAL :: ZFLOOD_AREA
REAL :: ZTOT_AREA
REAL :: ZRESIDU
REAL :: ZREAD_IN,ZREAD_OUT
!
LOGICAL :: GRETURN
!
INTEGER :: JBAS, JLON, JLAT, IFLOOD_COUNT
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('FLOOD_REDISTRIB',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
! * Init
!-------------------------------------------------------------------------------
!
PRESIDU   (:,:) = 0.0
PSRC_FLOOD(:,:) = 0.0
!
GRETURN=.FALSE.
!
!-------------------------------------------------------------------------------
! * Redistribution or not ?
!-------------------------------------------------------------------------------
!
WHERE(TPG%GMASK(:,:).AND.TP%XFFLOOD(:,:)>0.0)
  ILOC_FLOOD_IN(:,:)=1
ELSEWHERE
  ILOC_FLOOD_IN(:,:)=0
ENDWHERE
!
WHERE(TPG%GMASK(:,:).AND.PREAD(:,:)/=0.0)
  ILOC_FLUXE_IN(:,:)=1
ELSEWHERE
  ILOC_FLUXE_IN(:,:)=0
ENDWHERE
!
!
!-------------------------------------------------------------------------------
! * If floodplains at t and fluxes at t-1 are well co-localized : return
!-------------------------------------------------------------------------------
!
IF(GRETURN)THEN
  IF (LHOOK) CALL DR_HOOK('TRIP_OASIS_RECV:FLOOD_REDISTRIB',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!
IF(ALL(ILOC_FLUXE_IN(:,:)==0))THEN
!
   GRETURN=.TRUE.
!
ELSE
!        
   DO JLAT=1,KLAT
      DO JLON=1,KLON
         IF(ILOC_FLUXE_IN(JLON,JLAT)==ILOC_FLOOD_IN(JLON,JLAT))THEN
           GCOLLOCATED(JLON,JLAT)=.TRUE.
         ELSE
           GCOLLOCATED(JLON,JLAT)=.FALSE.
         ENDIF
      ENDDO
   ENDDO
!
   IF(ALL(GCOLLOCATED(:,:)))THEN
     GRETURN=.TRUE.
   ENDIF
!
ENDIF
!
!-------------------------------------------------------------------------------
! * Comput cumulated flooded areas
!-------------------------------------------------------------------------------
!
ZFLOOD_AREA = 0.0
!
DO JLAT=1,KLAT
   DO JLON=1,KLON
      IF(TPG%GMASK(JLON,JLAT).AND.TP%XFFLOOD(JLON,JLAT)>0.0)THEN
        ZFLOOD_AREA = ZFLOOD_AREA + TPG%XAREA(JLON,JLAT)
      ENDIF   
   ENDDO
ENDDO
!
!-------------------------------------------------------------------------------
! * If no flooded areas, redistribute the redidue over the entire domain
!-------------------------------------------------------------------------------
!
IF(ZFLOOD_AREA==0.0)THEN
!
  ZTOT_AREA=0.0
  ZRESIDU  =0.0
!  
  DO JLAT=1,KLAT
     DO JLON=1,KLON
        IF(TPG%GMASK(JLON,JLAT))THEN
          ZTOT_AREA = ZTOT_AREA +                  TPG%XAREA(JLON,JLAT)
          ZRESIDU   = ZRESIDU   + PREAD(JLON,JLAT)*TPG%XAREA(JLON,JLAT) ![kg/s]
        ENDIF
     ENDDO
  ENDDO
!
  PRESIDU(:,:) = ZRESIDU * TPG%XAREA(:,:) / ZTOT_AREA ![kg/s]
  PSRC_FLOOD(:,:) = 0.0  
!
  IF (LHOOK) CALL DR_HOOK('TRIP_OASIS_RECV:FLOOD_REDISTRIB',1,ZHOOK_HANDLE)
  RETURN
!
ENDIF
!
!-------------------------------------------------------------------------------
! * If flooded areas at time t, redistribute the redidue over each basin
!-------------------------------------------------------------------------------
!
ZOUT(:,:) = 0.0
!
WHERE(TPG%GMASK(:,:).AND.TP%XFFLOOD(:,:)>0.0)
  ZOUT(:,:) = PREAD(:,:)
ENDWHERE
!
! Basin redistribution
!
ZBAS_IN   (:) = 0.0
ZBAS_OUT  (:) = 0.0
IBAS_NCELL(:) = 0
!
DO JLAT=1,KLAT
   DO JLON=1,KLON
      JBAS=TPG%NBASID(JLON,JLAT)
      IF(JBAS>0)THEN
        ZBAS_IN(JBAS)=ZBAS_IN(JBAS)+PREAD(JLON,JLAT)*TPG%XAREA(JLON,JLAT)
      ENDIF
      IF(TP%XFFLOOD(JLON,JLAT)>0.0)THEN
        IBAS_NCELL(JBAS)=IBAS_NCELL(JBAS)+1.0
        ZBAS_OUT  (JBAS)=ZBAS_OUT  (JBAS)+PREAD(JLON,JLAT)*TPG%XAREA(JLON,JLAT)
      ENDIF
   ENDDO
ENDDO
!
!
ZBAS_FACTOR(:)=1.0
!
WHERE(IBAS_NCELL(:)>0.AND.ZBAS_OUT(:)/=0.0)
  ZBAS_FACTOR(:)=ZBAS_FACTOR(:)+(ZBAS_IN(:)-ZBAS_OUT(:))/ZBAS_OUT(:)
ENDWHERE
!
DO JLAT=1,KLAT
   DO JLON=1,KLON
      JBAS=TPG%NBASID(JLON,JLAT)
      IF(JBAS>0)THEN
        ZOUT(JLON,JLAT) = ZOUT(JLON,JLAT) * ZBAS_FACTOR(JBAS)
      ENDIF
     ENDDO
ENDDO
!
!-------------------------------------------------------------------------------
! * Ensure final global conservation
!-------------------------------------------------------------------------------
!
ZREAD_IN = 0.0
ZREAD_OUT= 0.0
!
DO JLAT=1,KLAT
   DO JLON=1,KLON
      IF(TPG%GMASK(JLON,JLAT))THEN
        ZREAD_IN  = ZREAD_IN  + PREAD(JLON,JLAT)*TPG%XAREA(JLON,JLAT) ![kg/s]
        ZREAD_OUT = ZREAD_OUT + ZOUT (JLON,JLAT)*TPG%XAREA(JLON,JLAT) ![kg/s]              
      ENDIF
   ENDDO
ENDDO
!
WHERE(TPG%GMASK(:,:).AND.TP%XFFLOOD(:,:)>0.0)
  PSRC_FLOOD(:,:) = ZOUT(:,:) + (ZREAD_IN-ZREAD_OUT)/ZFLOOD_AREA ![kg/m2/s]
ELSEWHERE
  PSRC_FLOOD(:,:)=0.0        
ENDWHERE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('FLOOD_REDISTRIB',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE FLOOD_REDISTRIB
