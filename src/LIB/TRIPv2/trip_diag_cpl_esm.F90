SUBROUTINE TRIP_DIAG_CPL_ESM(PTSTEP_RUN,PDISCHARGE,PCALVING,PWTD,PFWTD)  
!     #################################################################
!
!!****  *TRIP_DIAG_CPL_ESM*  
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
USE MODD_TRIP_PAR,   ONLY : XRHOLW
!
USE MODD_TRIP_OASIS, ONLY : LCPL_SEA, LCPL_LAND, LCPL_GW,    &
                            LCPL_FLOOD, LCPL_CALVSEA
!
USE MODD_TRIP_GRID, ONLY : TGRID
USE MODD_TRIP,      ONLY : TTRIP
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,                 INTENT(IN) :: PTSTEP_RUN !Run  timestep                   [s]
REAL, DIMENSION(:,:), INTENT(IN) :: PDISCHARGE !Cumulated river discharges      [kg]
REAL, DIMENSION(:,:), INTENT(IN) :: PCALVING   !Input claving flux from glacier [kg/s]
REAL, DIMENSION(:,:), INTENT(IN) :: PWTD       !Water table depth               [m]
REAL, DIMENSION(:,:), INTENT(IN) :: PFWTD      !Fraction of Water table to rise
!
!*      0.2    declarations of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TRIP_DIAG_CPL_ESM',0,ZHOOK_HANDLE)
!
!*       1.      Actualisation of sea coupling diagnostic:
!               ------------------------------------------
!
!
! River discharges to ocean [kg/m2]
!
IF(LCPL_SEA)THEN
  WHERE(TGRID%NGRCN(:,:)==9.OR.TGRID%NGRCN(:,:)==12)
    TTRIP%XCPL_RIVDIS(:,:) = TTRIP%XCPL_RIVDIS(:,:) + PDISCHARGE(:,:) / TGRID%XAREA(:,:)
  ENDWHERE
ENDIF
!
! Calving flux over greenland and antarctica [kg/m2]
!
IF(LCPL_CALVSEA)THEN
  WHERE(TGRID%GMASK_GRE(:,:))
    TTRIP%XCPL_CALVGRE(:,:) = TTRIP%XCPL_CALVGRE(:,:) + PCALVING(:,:) * PTSTEP_RUN / TGRID%XAREA(:,:)
  ENDWHERE
  WHERE(TGRID%GMASK_ANT(:,:))
    TTRIP%XCPL_CALVANT(:,:) = TTRIP%XCPL_CALVANT(:,:) + PCALVING(:,:) * PTSTEP_RUN / TGRID%XAREA(:,:)
  ENDWHERE
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       2.      Actualisation of land coupling diagnostic:
!               -------------------------------------------
!
IF(LCPL_LAND)THEN
!
! Water table depth and fraction of water table to rise
!
  IF(LCPL_GW)THEN
    WHERE(TGRID%GMASK_GW(:,:))
          TTRIP%XCPL_WTD (:,:) = MAX(PWTD(:,:),0.0)
          TTRIP%XCPL_FWTD(:,:) = PFWTD(:,:)
    ENDWHERE
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
IF (LHOOK) CALL DR_HOOK('TRIP_DIAG_CPL_ESM',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE TRIP_DIAG_CPL_ESM
