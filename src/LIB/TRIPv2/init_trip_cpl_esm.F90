SUBROUTINE INIT_TRIP_CPL_ESM (TP, TPG, &
                              KLON,KLAT)  
!     ##################################
!
!!****  *INIT_TRIP_CPL_ESM*  
!!
!!    PURPOSE
!!    -------
!
!     TRIP cpl diag compuation
!     
!!      
!!    AUTHOR
!!    ------
!!      B. Decharme
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    12/12/13 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!
USE MODD_TRIP, ONLY : TRIP_t
USE MODD_TRIP_GRID, ONLY : TRIP_GRID_t
!
USE MODD_TRIP_OASIS, ONLY : LCPL_SEA, LCPL_LAND, LCPL_GW,    &
                            LCPL_FLOOD, LCPL_CALVSEA
!
USE MODD_TRIP_PAR,  ONLY : XUNDEF
!
!
USE MODI_GWF_CPL_UPDATE
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
INTEGER, INTENT(IN) :: KLON
INTEGER, INTENT(IN) :: KLAT
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(KLON,KLAT) :: ZHG_OLD    !Water table elevation at t-1    [m]
REAL, DIMENSION(KLON,KLAT) :: ZWTD       !Water table depth               [m]
REAL, DIMENSION(KLON,KLAT) :: ZFWTD      !Fraction of Water table to rise
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('INIT_TRIP_CPL_ESM',0,ZHOOK_HANDLE)
!
ZHG_OLD (:,:) = XUNDEF
ZWTD    (:,:) = XUNDEF
ZFWTD   (:,:) = XUNDEF
!
!-------------------------------------------------------------------------------
! * Allocate coupling variables
!-------------------------------------------------------------------------------
!
IF(LCPL_SEA)THEN
  ALLOCATE(TP%XCPL_RIVDIS(KLON,KLAT))
  TP%XCPL_RIVDIS(:,:) = XUNDEF
ELSE
  ALLOCATE(TP%XCPL_RIVDIS(0,0))
ENDIF
!
IF(LCPL_CALVSEA)THEN
  ALLOCATE(TP%XCPL_CALVGRE(KLON,KLAT))
  ALLOCATE(TP%XCPL_CALVANT(KLON,KLAT))
  TP%XCPL_CALVGRE(:,:) = XUNDEF
  TP%XCPL_CALVANT(:,:) = XUNDEF
ELSE
  ALLOCATE(TP%XCPL_CALVGRE(0,0))
  ALLOCATE(TP%XCPL_CALVANT(0,0))
ENDIF
!
IF(LCPL_LAND)THEN
  IF(LCPL_GW)THEN
    ALLOCATE(TP%XCPL_FWTD(KLON,KLAT))
    ALLOCATE(TP%XCPL_WTD (KLON,KLAT))
    TP%XCPL_FWTD(:,:) = XUNDEF
    TP%XCPL_WTD (:,:) = XUNDEF
  ELSE
    ALLOCATE(TP%XCPL_FWTD(0,0))
    ALLOCATE(TP%XCPL_WTD (0,0))
  ENDIF
  IF(LCPL_FLOOD)THEN
    ALLOCATE(TP%XCPL_FFLOOD (KLON,KLAT))
    ALLOCATE(TP%XCPL_PIFLOOD(KLON,KLAT))
    TP%XCPL_FFLOOD (:,:) = XUNDEF
    TP%XCPL_PIFLOOD(:,:) = XUNDEF
  ELSE
    ALLOCATE(TP%XCPL_FFLOOD (0,0))
    ALLOCATE(TP%XCPL_PIFLOOD(0,0))
  ENDIF
ELSE
  ALLOCATE(TP%XCPL_FFLOOD (0,0))
  ALLOCATE(TP%XCPL_PIFLOOD(0,0))
  ALLOCATE(TP%XCPL_FWTD(0,0))
  ALLOCATE(TP%XCPL_WTD (0,0))
ENDIF
!
!-------------------------------------------------------------------------------
! * Actualisation of coupling diagnostic:
!-------------------------------------------------------------------------------
!
IF(LCPL_SEA)THEN
  WHERE(TPG%NGRCN(:,:)==9.OR.TPG%NGRCN(:,:)==12)
    TP%XCPL_RIVDIS(:,:) = 0.0
  ENDWHERE 
ENDIF
!
IF(LCPL_CALVSEA)THEN
  WHERE(TPG%GMASK_GRE(:,:))
    TP%XCPL_CALVGRE(:,:) = 0.0
  ENDWHERE 
  WHERE(TPG%GMASK_ANT(:,:))
    TP%XCPL_CALVANT(:,:) = 0.0
  ENDWHERE 
ENDIF
!
IF(LCPL_LAND)THEN
!
! Water table depth and fraction of water table to rise
!
  IF(LCPL_GW)THEN
!
    CALL GWF_CPL_UPDATE(TP%XTABGW_H,TP%XTABGW_F,TPG%GMASK_GW,&
                        TP%XTOPO_RIV,TP%XHGROUND,ZHG_OLD,    &
                        ZWTD,ZFWTD                           )
!
    WHERE(TPG%GMASK_GW(:,:))
          TP%XCPL_WTD (:,:) = ZWTD (:,:)
          TP%XCPL_FWTD(:,:) = ZFWTD(:,:)
    ENDWHERE
!    
  ENDIF
!
! Flood fraction [-] and potential infiltration [kg/m2]
!       
  IF(LCPL_FLOOD)THEN
    TP%XCPL_FFLOOD (:,:) = TP%XFFLOOD    (:,:)
    TP%XCPL_PIFLOOD(:,:) = TP%XFLOOD_STO (:,:) / TPG%XAREA(:,:)
  ENDIF
!  
ENDIF
!
IF (LHOOK) CALL DR_HOOK('INIT_TRIP_CPL_ESM',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE INIT_TRIP_CPL_ESM
