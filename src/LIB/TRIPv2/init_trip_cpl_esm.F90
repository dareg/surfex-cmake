SUBROUTINE INIT_TRIP_CPL_ESM(KLON,KLAT)  
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
!!	B. Decharme     
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    12/12/13 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_TRIP_OASIS, ONLY : LCPL_SEA, LCPL_LAND, LCPL_GW,    &
                            LCPL_FLOOD, LCPL_CALVSEA
!
USE MODD_TRIP_PAR,  ONLY : XUNDEF
!
USE MODD_TRIP_GRID, ONLY : TGRID
USE MODD_TRIP,      ONLY : TTRIP
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
INTEGER, INTENT(IN) :: KLON
INTEGER, INTENT(IN) :: KLAT
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(KLON,KLAT) :: ZHG_OLD    !Water table elevation at t-1    [m]
REAL, DIMENSION(KLON,KLAT) :: ZWTD       !Water table depth               [m]
REAL, DIMENSION(KLON,KLAT) :: ZFWTD      !Fraction of Water table to rise
REAL, DIMENSION(KLON,KLAT) :: ZWTDELEV   !Water table depth / elevation   [m]
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
ZWTDELEV(:,:) = XUNDEF
!
!-------------------------------------------------------------------------------
! * Allocate coupling variables
!-------------------------------------------------------------------------------
!
IF(LCPL_SEA)THEN
  ALLOCATE(TTRIP%XCPL_RIVDIS(KLON,KLAT))
  TTRIP%XCPL_RIVDIS(:,:) = XUNDEF
ELSE
  ALLOCATE(TTRIP%XCPL_RIVDIS(0,0))
ENDIF
!
IF(LCPL_CALVSEA)THEN
  ALLOCATE(TTRIP%XCPL_CALVGRE(KLON,KLAT))
  ALLOCATE(TTRIP%XCPL_CALVANT(KLON,KLAT))
  TTRIP%XCPL_CALVGRE(:,:) = XUNDEF
  TTRIP%XCPL_CALVANT(:,:) = XUNDEF
ELSE
  ALLOCATE(TTRIP%XCPL_CALVGRE(0,0))
  ALLOCATE(TTRIP%XCPL_CALVANT(0,0))
ENDIF
!
IF(LCPL_LAND)THEN
  IF(LCPL_GW)THEN
    ALLOCATE(TTRIP%XCPL_FWTD(KLON,KLAT))
    ALLOCATE(TTRIP%XCPL_WTD (KLON,KLAT))
    TTRIP%XCPL_FWTD(:,:) = XUNDEF
    TTRIP%XCPL_WTD (:,:) = XUNDEF
  ELSE
    ALLOCATE(TTRIP%XCPL_FWTD(0,0))
    ALLOCATE(TTRIP%XCPL_WTD (0,0))
  ENDIF
  IF(LCPL_FLOOD)THEN
    ALLOCATE(TTRIP%XCPL_FFLOOD (KLON,KLAT))
    ALLOCATE(TTRIP%XCPL_PIFLOOD(KLON,KLAT))
    TTRIP%XCPL_FFLOOD (:,:) = XUNDEF
    TTRIP%XCPL_PIFLOOD(:,:) = XUNDEF
  ELSE
    ALLOCATE(TTRIP%XCPL_FFLOOD (0,0))
    ALLOCATE(TTRIP%XCPL_PIFLOOD(0,0))
  ENDIF
ELSE
  ALLOCATE(TTRIP%XCPL_FFLOOD (0,0))
  ALLOCATE(TTRIP%XCPL_PIFLOOD(0,0))
  ALLOCATE(TTRIP%XCPL_FWTD(0,0))
  ALLOCATE(TTRIP%XCPL_WTD (0,0))
ENDIF
!
!-------------------------------------------------------------------------------
! * Actualisation of coupling diagnostic:
!-------------------------------------------------------------------------------
!
IF(LCPL_SEA)THEN
  WHERE(TGRID%NGRCN(:,:)==9.OR.TGRID%NGRCN(:,:)==12)
    TTRIP%XCPL_RIVDIS(:,:) = 0.0
  ENDWHERE 
ENDIF
!
IF(LCPL_LAND)THEN
!
! Water table depth and fraction of water table to rise
!
  IF(LCPL_GW)THEN
!
    CALL GWF_CPL_UPDATE(TTRIP%XTABGW_H,TTRIP%XTABGW_F,TGRID%GMASK_GW, &
                        TTRIP%XTOPO_RIV,TTRIP%XELEV,TTRIP%XHGROUND,   &
                        ZHG_OLD,ZWTD,ZFWTD,ZWTDELEV                   )
!
    WHERE(TGRID%GMASK_GW(:,:))
          TTRIP%XCPL_WTD (:,:) = MAX(-1.0*ZWTD(:,:),0.0)
          TTRIP%XCPL_FWTD(:,:) = ZFWTD(:,:)
    ENDWHERE
!    
  ENDIF
!
! Flood fraction [-] and potential infiltration [kg/m2]
!       
  IF(LCPL_FLOOD)THEN
    TTRIP%XCPL_FFLOOD (:,:) = TTRIP%XFFLOOD    (:,:)
    TTRIP%XCPL_PIFLOOD(:,:) = TTRIP%XFLOOD_STO (:,:) / TGRID%XAREA(:,:)
  ENDIF
!  
ENDIF
!
IF (LHOOK) CALL DR_HOOK('INIT_TRIP_CPL_ESM',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE INIT_TRIP_CPL_ESM
