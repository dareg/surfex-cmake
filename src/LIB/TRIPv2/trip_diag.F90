SUBROUTINE TRIP_DIAG(TPDG, TP, TPG, &
                     PTSTEP,PSOUT,PSIN,PVEL,PHS,PGOUT,PGNEG,    &
                     PWTD,PFWTD,PQGCELL,PHGHS, &
                     PQFR,PQRF,PVFIN,PVFOUT,PHSF,PSRC_FLOOD,    &
                     PDRAIN,PRUNOFF,PDISCHARGE  )
!     #####################################################
!
!!****  *TRIP_DIAG*  
!!
!!    PURPOSE
!!    -------
!
!     TRIP diag compuation
!     
!!      
!!    AUTHOR
!!    ------
!!      B. Decharme
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    12/12/13 
!!      09/16   B. Decharme  limit wtd to -1000m
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!
USE MODD_TRIP_DIAG, ONLY : TRIP_DIAG_t
USE MODD_TRIP,      ONLY : TRIP_t
USE MODD_TRIP_GRID, ONLY : TRIP_GRID_t
!
USE MODN_TRIP_RUN,   ONLY : LDIAG_MISC
USE MODD_TRIP_OASIS, ONLY : LCPL_LAND, LCPL_FLOOD
!
USE MODN_TRIP,       ONLY : CGROUNDW, CVIT, LFLOOD
!
USE MODD_TRIP_PAR,   ONLY : XUNDEF, XRHOLW, XGWDZMAX
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
TYPE(TRIP_DIAG_t), INTENT(INOUT) :: TPDG
TYPE(TRIP_t),      INTENT(INOUT) :: TP
TYPE(TRIP_GRID_t), INTENT(INOUT) :: TPG
!
REAL,                 INTENT(IN)  :: PTSTEP     !Time step                     [s]
!
REAL, DIMENSION(:,:), INTENT(IN)  :: PSOUT      !streamflow                    [kg/s]
REAL, DIMENSION(:,:), INTENT(IN)  :: PSIN       !grid-cell input streamflow    [kg/s]
REAL, DIMENSION(:,:), INTENT(IN)  :: PVEL       !river velocity                [m/s]
REAL, DIMENSION(:,:), INTENT(IN)  :: PHS        !River heigh                   [m]
REAL, DIMENSION(:,:), INTENT(IN)  :: PGOUT      !Groundwater outflow           [kg/s]
REAL, DIMENSION(:,:), INTENT(IN)  :: PGNEG      !Groundwater intflow (neg)     [kg/s]
REAL, DIMENSION(:,:), INTENT(IN)  :: PWTD       !Water table depth for coupling[m]
REAL, DIMENSION(:,:), INTENT(IN)  :: PFWTD      !Fraction of water table to rise
REAL, DIMENSION(:,:), INTENT(IN)  :: PQGCELL    !lateral groundwater exchanges [kg/s]
REAL, DIMENSION(:,:), INTENT(IN)  :: PHGHS      !groundwater minus river heigh [m]
REAL, DIMENSION(:,:), INTENT(IN)  :: PQFR       !floodplains to river exchange [kg/s]
REAL, DIMENSION(:,:), INTENT(IN)  :: PQRF       !river to floodplains exchange [kg/s]
REAL, DIMENSION(:,:), INTENT(IN)  :: PVFIN      !QRF velocity                  [m/s]
REAL, DIMENSION(:,:), INTENT(IN)  :: PVFOUT     !QFR velocity                  [m/s]
REAL, DIMENSION(:,:), INTENT(IN)  :: PHSF       !river minus flodd heigh       [m]
REAL, DIMENSION(:,:), INTENT(IN)  :: PSRC_FLOOD !P-E-I flood source term       [kg/s]
REAL, DIMENSION(:,:), INTENT(IN)  :: PDRAIN     !Input drainage or recharge    [kg/s]
REAL, DIMENSION(:,:), INTENT(IN)  :: PRUNOFF    !Input surface runoff          [kg/s]
REAL, DIMENSION(:,:), INTENT(OUT) :: PDISCHARGE !Cumulated river discharges    [kg]
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PSOUT,1),SIZE(PSOUT,2)) :: ZWORK
!
REAL, DIMENSION(SIZE(PSOUT,1),SIZE(PSOUT,2)) :: ZQGWR       ![kg/s]
REAL, DIMENSION(SIZE(PSOUT,1),SIZE(PSOUT,2)) :: ZSURF_STO   ![kg.s/m2]
REAL, DIMENSION(SIZE(PSOUT,1),SIZE(PSOUT,2)) :: ZGROUND_STO ![kg.s/m2]
REAL, DIMENSION(SIZE(PSOUT,1),SIZE(PSOUT,2)) :: ZFLOOD_STO  ![kg.s/m2]
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TRIP_DIAG',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
! * Init
!-------------------------------------------------------------------------------
!
ZSURF_STO  (:,:) = 0.0
ZGROUND_STO(:,:) = 0.0
ZQGWR      (:,:) = 0.0
ZFLOOD_STO (:,:) = 0.0
!
IF(CGROUNDW/='DEF')THEN
  WHERE(TPG%GMASK_GW(:,:))
        ZQGWR(:,:) = PDRAIN(:,:)
  ENDWHERE
ENDIF
!
!-------------------------------------------------------------------------------
! * River mass, fluxes, and velocity
!-------------------------------------------------------------------------------
!
ZSURF_STO           (:,:) = TP%XSURF_STO(:,:) / TPG%XAREA(:,:)
TPDG%TDIAG%XSURF_STO(:,:) = TPDG%TDIAG%XSURF_STO(:,:) + ZSURF_STO(:,:) * PTSTEP
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
PDISCHARGE          (:,:) = PDISCHARGE      (:,:) + PSOUT(:,:) * PTSTEP          ![kg]
TPDG%TDIAG%XQDIS    (:,:) = TPDG%TDIAG%XQDIS(:,:) + PSOUT(:,:) * (PTSTEP/XRHOLW) ![m3]
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!For CMIP, QIN include all input flux at the surface (runoff and drainage) 
ZWORK          (:,:) = PSIN(:,:)+PRUNOFF(:,:)+PDRAIN(:,:)-ZQGWR(:,:)
TPDG%TDIAG%XQIN(:,:) = TPDG%TDIAG%XQIN(:,:) + ZWORK(:,:) * (PTSTEP/XRHOLW) ![m3]
!TPDG%TDIAG%XQIN(:,:) = TPDG%TDIAG%XQIN(:,:) + PSIN(:,:) * (PTSTEP/XRHOLW) ![m3]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
IF(CVIT=='VAR')THEN
!
   TPDG%TDIAG%XVEL(:,:) = TPDG%TDIAG%XVEL(:,:) + PVEL(:,:) * PTSTEP
!
   WHERE(TPG%GMASK_VEL(:,:))
         TPDG%TDIAG%XHS(:,:) = TPDG%TDIAG%XHS(:,:) + PHS (:,:) * PTSTEP
  ENDWHERE
!
ENDIF
!
!-------------------------------------------------------------------------------
! * Groundwater mass, fluxes, and depth
!-------------------------------------------------------------------------------
!
IF(CGROUNDW=='CST')THEN  
  WHERE(TPG%GMASK_GW(:,:))
        ZGROUND_STO(:,:) =  TP%XGROUND_STO(:,:) / TPG%XAREA(:,:)
  ENDWHERE
ELSEIF(CGROUNDW=='DIF')THEN
!
  WHERE(TPG%GMASK_GW(:,:))
        ZGROUND_STO        (:,:) = (XGWDZMAX+PWTD(:,:)) * TP%XWEFF(:,:) * XRHOLW
        TPDG%TDIAG%XHGROUND(:,:) = TPDG%TDIAG%XHGROUND(:,:) + TP%XHGROUND(:,:) * PTSTEP
        TPDG%TDIAG%XWTD    (:,:) = TPDG%TDIAG%XWTD    (:,:) - PWTD       (:,:) * PTSTEP !positive downward
  ENDWHERE
!
ENDIF
!
IF(CGROUNDW/='DEF')THEN
!
  ZWORK(:,:) =0.0
  WHERE(TPG%GMASK_GW(:,:))
        ZWORK(:,:) = PGOUT(:,:)+PGNEG(:,:)
  ENDWHERE
!
  WHERE(TPG%GMASK_GW(:,:))
        TPDG%TDIAG%XGROUND_STO(:,:) = TPDG%TDIAG%XGROUND_STO(:,:) + ZGROUND_STO(:,:) * PTSTEP
        TPDG%TDIAG%XQGF       (:,:) = TPDG%TDIAG%XQGF       (:,:) + ZWORK      (:,:) * PTSTEP / TPG%XAREA(:,:)
  ENDWHERE
! 
ENDIF
!
!-------------------------------------------------------------------------------
! * Floodplains
!-------------------------------------------------------------------------------
!
IF(LFLOOD)THEN    
!
  WHERE(TPG%GMASK_FLD(:,:))
!
        ZFLOOD_STO           (:,:) =  TP%XFLOOD_STO(:,:) / TPG%XAREA(:,:)
!
        TPDG%TDIAG%XFLOOD_STO(:,:) = TPDG%TDIAG%XFLOOD_STO(:,:) + ZFLOOD_STO(:,:)
        TPDG%TDIAG%XFF       (:,:) = TPDG%TDIAG%XFF       (:,:) + TP%XFFLOOD(:,:) * PTSTEP
        TPDG%TDIAG%XHF       (:,:) = TPDG%TDIAG%XHF       (:,:) + TP%XHFLOOD(:,:) * PTSTEP
!
  ENDWHERE
!
ENDIF
!
!-------------------------------------------------------------------------------
! * Forcing variables can be used to force TRIP offline
!-------------------------------------------------------------------------------
!
IF(LCPL_LAND)THEN
  TPDG%TDIAG%XRUNOFF(:,:) = TPDG%TDIAG%XRUNOFF(:,:) + PRUNOFF(:,:) * PTSTEP / TPG%XAREA(:,:)
  TPDG%TDIAG%XDRAIN (:,:) = TPDG%TDIAG%XDRAIN (:,:) + PDRAIN (:,:) * PTSTEP / TPG%XAREA(:,:)
ENDIF
!
IF(LCPL_FLOOD)THEN
  TPDG%TDIAG%XSOURCE(:,:) = TPDG%TDIAG%XSOURCE(:,:) + PSRC_FLOOD(:,:) * PTSTEP / TPG%XAREA(:,:)
ENDIF
!
!-------------------------------------------------------------------------------
! * MISC fields
!-------------------------------------------------------------------------------
!
IF(LDIAG_MISC)THEN
!  
  IF(CGROUNDW=='DIF')THEN
    WHERE(TPG%GMASK_GW(:,:))
          TPDG%TDIAG%XFWTD  (:,:) = TPDG%TDIAG%XFWTD  (:,:) + PFWTD  (:,:) * PTSTEP
          TPDG%TDIAG%XQGCELL(:,:) = TPDG%TDIAG%XQGCELL(:,:) + PQGCELL(:,:) * PTSTEP / TPG%XAREA(:,:)
          TPDG%TDIAG%XHGHS  (:,:) = TPDG%TDIAG%XHGHS  (:,:) + PHGHS  (:,:) * PTSTEP
    ENDWHERE
  ENDIF
!
  IF(LFLOOD)THEN    
    WHERE(TPG%GMASK_FLD(:,:))
          TPDG%TDIAG%XWF    (:,:) = TPDG%TDIAG%XWF    (:,:) + TP%XWFLOOD   (:,:) * PTSTEP
          TPDG%TDIAG%XLF    (:,:) = TPDG%TDIAG%XLF    (:,:) + TP%XFLOOD_LEN(:,:) * PTSTEP
          TPDG%TDIAG%XHSF   (:,:) = TPDG%TDIAG%XHSF   (:,:) + PHSF         (:,:) * PTSTEP
    ENDWHERE
  ENDIF
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('TRIP_DIAG',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE TRIP_DIAG
