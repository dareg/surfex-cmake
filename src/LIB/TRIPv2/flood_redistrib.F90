!     #########
      SUBROUTINE FLOOD_REDISTRIB(KLON,KLAT,PPFLOOD,PEFLOOD,PIFLOOD,PRESIDU)  
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
USE MODD_TRIP_PAR, ONLY : XUNDEF, XRHOLW
!
USE MODD_TRIP,       ONLY : TTRIP
USE MODD_TRIP_GRID,  ONLY : TGRID
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
INTEGER, INTENT(IN)               :: KLON
INTEGER, INTENT(IN)               :: KLAT
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PPFLOOD ![kg/m2]
REAL, DIMENSION(:,:), INTENT(INOUT) :: PEFLOOD ![kg/m2]
REAL, DIMENSION(:,:), INTENT(INOUT) :: PIFLOOD ![kg/m2]
REAL, DIMENSION(:,:), INTENT(OUT) :: PRESIDU
!
!
!*      0.2    declarations of local variables
!
REAL,    DIMENSION(TGRID%NBASMAX) :: ZBAS_P_IN
REAL,    DIMENSION(TGRID%NBASMAX) :: ZBAS_E_IN
REAL,    DIMENSION(TGRID%NBASMAX) :: ZBAS_I_IN
REAL,    DIMENSION(TGRID%NBASMAX) :: ZBAS_P_OUT
REAL,    DIMENSION(TGRID%NBASMAX) :: ZBAS_E_OUT
REAL,    DIMENSION(TGRID%NBASMAX) :: ZBAS_I_OUT
REAL,    DIMENSION(TGRID%NBASMAX) :: ZBAS_P_FACTOR
REAL,    DIMENSION(TGRID%NBASMAX) :: ZBAS_E_FACTOR
REAL,    DIMENSION(TGRID%NBASMAX) :: ZBAS_I_FACTOR
INTEGER, DIMENSION(TGRID%NBASMAX) :: IBAS_NCELL
!
REAL,    DIMENSION(KLON,KLAT) :: ZPFLD_OUT
REAL,    DIMENSION(KLON,KLAT) :: ZEFLD_OUT
REAL,    DIMENSION(KLON,KLAT) :: ZIFLD_OUT
INTEGER, DIMENSION(KLON,KLAT) :: ILOC_FLOOD_IN
INTEGER, DIMENSION(KLON,KLAT) :: ILOC_FLUXE_IN
!
REAL :: ZFLOOD_AREA
REAL :: ZTOT_AREA
REAL :: ZRESIDU
REAL :: ZP_IN,ZP_OUT
REAL :: ZE_IN,ZE_OUT
REAL :: ZI_IN,ZI_OUT
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
PRESIDU(:,:)=0.0
!
!-------------------------------------------------------------------------------
! * Redistribution or not ?
!-------------------------------------------------------------------------------
!
WHERE(TGRID%GMASK(:,:).AND.TTRIP%XFFLOOD(:,:)>0.0)
  ILOC_FLOOD_IN(:,:)=1
ELSEWHERE
  ILOC_FLOOD_IN(:,:)=0
ENDWHERE
!
WHERE(TGRID%GMASK(:,:).AND.(PPFLOOD(:,:)/=0.0.OR.PEFLOOD(:,:)/=0.0.OR.PIFLOOD(:,:)/=0.0))
  ILOC_FLUXE_IN(:,:)=1
ELSEWHERE
  ILOC_FLUXE_IN(:,:)=0
ENDWHERE
!
GRETURN=.TRUE.
DO JLAT=1,KLAT
   DO JLON=1,KLON
      IF(ILOC_FLUXE_IN(JLON,JLAT)==1.AND.ILOC_FLOOD_IN(JLON,JLAT)==0)THEN
        GRETURN=.FALSE.
      ENDIF
   ENDDO
ENDDO
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
!-------------------------------------------------------------------------------
! * Comput cumulated flooded areas
!-------------------------------------------------------------------------------
!
ZFLOOD_AREA = 0.0
!
DO JLAT=1,KLAT
   DO JLON=1,KLON
      IF(TGRID%GMASK(JLON,JLAT).AND.TTRIP%XFFLOOD(JLON,JLAT)>0.0)THEN
        ZFLOOD_AREA = ZFLOOD_AREA + TGRID%XAREA(JLON,JLAT)
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
        IF(TGRID%GMASK(JLON,JLAT))THEN
          ZTOT_AREA = ZTOT_AREA +                    TGRID%XAREA(JLON,JLAT)
          ZRESIDU   = ZRESIDU   + PPFLOOD(JLON,JLAT)*TGRID%XAREA(JLON,JLAT) &
                                - PEFLOOD(JLON,JLAT)*TGRID%XAREA(JLON,JLAT) &
                                - PIFLOOD(JLON,JLAT)*TGRID%XAREA(JLON,JLAT) ![kg]
        ENDIF
     ENDDO
  ENDDO
!
  PRESIDU(:,:) = ZRESIDU / ZTOT_AREA ![kg/m2]
  PPFLOOD(:,:) = 0.0
  PEFLOOD(:,:) = 0.0
  PIFLOOD(:,:) = 0.0
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
ZPFLD_OUT(:,:) = 0.0
ZEFLD_OUT(:,:) = 0.0
ZIFLD_OUT(:,:) = 0.0
!
WHERE(TGRID%GMASK(:,:).AND.TTRIP%XFFLOOD(:,:)>0.0)
   ZPFLD_OUT(:,:) = PPFLOOD(:,:)
   ZEFLD_OUT(:,:) = PEFLOOD(:,:)
   ZIFLD_OUT(:,:) = PIFLOOD(:,:)
ENDWHERE
!
! Basin redistribution
!
ZBAS_P_IN(:)=0.0
ZBAS_E_IN(:)=0.0
ZBAS_I_IN(:)=0.0
!
ZBAS_P_OUT(:)=0.0
ZBAS_E_OUT(:)=0.0
ZBAS_I_OUT(:)=0.0
!
IBAS_NCELL(:)=0
!
DO JBAS=TGRID%NBASMIN,TGRID%NBASMAX
   DO JLAT=1,KLAT
      DO JLON=1,KLON
         IF(TGRID%NBASID(JLON,JLAT)==JBAS)THEN
            IF(TTRIP%XFFLOOD(JLON,JLAT)>0.0)THEN
              IBAS_NCELL(JBAS)=IBAS_NCELL(JBAS)+1.0
              ZBAS_P_OUT(JBAS)=ZBAS_P_OUT(JBAS)+PPFLOOD(JLON,JLAT)*TGRID%XAREA(JLON,JLAT)
              ZBAS_E_OUT(JBAS)=ZBAS_E_OUT(JBAS)+PEFLOOD(JLON,JLAT)*TGRID%XAREA(JLON,JLAT)
              ZBAS_I_OUT(JBAS)=ZBAS_I_OUT(JBAS)+PIFLOOD(JLON,JLAT)*TGRID%XAREA(JLON,JLAT)
            ENDIF
            ZBAS_P_IN (JBAS)=ZBAS_P_IN(JBAS)+PPFLOOD(JLON,JLAT)*TGRID%XAREA(JLON,JLAT)
            ZBAS_E_IN (JBAS)=ZBAS_E_IN(JBAS)+PEFLOOD(JLON,JLAT)*TGRID%XAREA(JLON,JLAT)
            ZBAS_I_IN (JBAS)=ZBAS_I_IN(JBAS)+PIFLOOD(JLON,JLAT)*TGRID%XAREA(JLON,JLAT)
         ENDIF
      ENDDO
   ENDDO
ENDDO
!
ZBAS_P_FACTOR(:)=1.0
ZBAS_E_FACTOR(:)=1.0
ZBAS_I_FACTOR(:)=1.0
!
WHERE(IBAS_NCELL(:)>0.AND.ZBAS_P_OUT(:)/=0.0)
  ZBAS_P_FACTOR(:)=ZBAS_P_FACTOR(:)+(ZBAS_P_IN(:)-ZBAS_P_OUT(:))/ZBAS_P_OUT(:)
ENDWHERE
WHERE(IBAS_NCELL(:)>0.AND.ZBAS_E_OUT(:)/=0.0)
  ZBAS_E_FACTOR(:)=ZBAS_E_FACTOR(:)+(ZBAS_E_IN(:)-ZBAS_E_OUT(:))/ZBAS_E_OUT(:)
ENDWHERE
WHERE(IBAS_NCELL(:)>0.AND.ZBAS_I_OUT(:)/=0.0)
  ZBAS_I_FACTOR(:)=ZBAS_I_FACTOR(:)+(ZBAS_I_IN(:)-ZBAS_I_OUT(:))/ZBAS_I_OUT(:)
ENDWHERE
!
DO JBAS=TGRID%NBASMIN,TGRID%NBASMAX
     DO JLAT=1,KLAT
        DO JLON=1,KLON
           IF(TGRID%NBASID(JLON,JLAT)==JBAS)THEN
             ZPFLD_OUT(JLON,JLAT) = ZPFLD_OUT(JLON,JLAT) * ZBAS_P_FACTOR(JBAS)
             ZEFLD_OUT(JLON,JLAT) = ZEFLD_OUT(JLON,JLAT) * ZBAS_E_FACTOR(JBAS)
             ZIFLD_OUT(JLON,JLAT) = ZIFLD_OUT(JLON,JLAT) * ZBAS_I_FACTOR(JBAS)
           ENDIF
        ENDDO
     ENDDO
ENDDO
!
!-------------------------------------------------------------------------------
! * Ensure final global conservation
!-------------------------------------------------------------------------------
!
ZP_IN = 0.0
ZE_IN = 0.0
ZI_IN = 0.0
!
ZP_OUT = 0.0
ZE_OUT = 0.0
ZI_OUT = 0.0
!
DO JLAT=1,KLAT
   DO JLON=1,KLON
      IF(TGRID%GMASK(JLON,JLAT))THEN
        ZP_IN  = ZP_IN  + PPFLOOD  (JLON,JLAT)*TGRID%XAREA(JLON,JLAT) ![kg]
        ZE_IN  = ZE_IN  + PEFLOOD  (JLON,JLAT)*TGRID%XAREA(JLON,JLAT) ![kg]
        ZI_IN  = ZI_IN  + PIFLOOD  (JLON,JLAT)*TGRID%XAREA(JLON,JLAT) ![kg]
        ZP_OUT = ZP_OUT + ZPFLD_OUT(JLON,JLAT)*TGRID%XAREA(JLON,JLAT) ![kg]
        ZE_OUT = ZE_OUT + ZEFLD_OUT(JLON,JLAT)*TGRID%XAREA(JLON,JLAT) ![kg]
        ZI_OUT = ZI_OUT + ZIFLD_OUT(JLON,JLAT)*TGRID%XAREA(JLON,JLAT) ![kg]
      ENDIF
   ENDDO
ENDDO
!
WHERE(TGRID%GMASK(:,:).AND.TTRIP%XFFLOOD(:,:)>0.0)
  PPFLOOD(:,:)=ZPFLD_OUT(:,:) + (ZP_IN-ZP_OUT)/ZFLOOD_AREA ![kg/m2]
  PEFLOOD(:,:)=ZEFLD_OUT(:,:) + (ZE_IN-ZE_OUT)/ZFLOOD_AREA ![kg/m2]
  PIFLOOD(:,:)=ZIFLD_OUT(:,:) + (ZI_IN-ZI_OUT)/ZFLOOD_AREA ![kg/m2]
ELSEWHERE
  PPFLOOD(:,:)=0.0
  PEFLOOD(:,:)=0.0
  PIFLOOD(:,:)=0.0
ENDWHERE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('FLOOD_REDISTRIB',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE FLOOD_REDISTRIB
