SUBROUTINE TRIP_DIAG(PTSTEP,PSOUT,PSIN,PVEL,PHS,PGOUT,PGNEG,    &
                     PWTD,PFWTD,PWTDRIV,PWTDELEV,PQGCELL,PHGHS, &
                     PQFR,PQRF,PVFIN,PVFOUT,PHSF,PSRC_FLOOD,    &
                     PDRAIN,PRUNOFF,PDISCHARGE                  )  
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
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODN_TRIP_RUN,   ONLY : LDIAG_MISC
USE MODD_TRIP_OASIS, ONLY : LCPL_LAND
!
USE MODN_TRIP,       ONLY : CGROUNDW, CVIT, LFLOOD
!
USE MODD_TRIP_PAR,   ONLY : XRHOLW
!
USE MODD_TRIP_GRID, ONLY : TPG => TRIP_GRID
USE MODD_TRIP, ONLY : TP => TRIP
USE MODD_TRIP_DIAG, ONLY : TPDG => TRIP_DIAG
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
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
REAL, DIMENSION(:,:), INTENT(IN)  :: PWTDRIV    !Water table depth / topo riv  [m]
REAL, DIMENSION(:,:), INTENT(IN)  :: PWTDELEV   !Water table depth / mean elev [m]
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
!
REAL, DIMENSION(:,:), INTENT(OUT) :: PDISCHARGE !Cumulated river discharges    [kg]
!
!*      0.2    declarations of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TRIP_DIAG',0,ZHOOK_HANDLE)
!
! * Actualisation of river discharge diags
!       
PDISCHARGE (:,:) = PDISCHARGE (:,:) + PSOUT(:,:) * PTSTEP          ![kg]
TPDG%TDIAG%XQDIS(:,:) = TPDG%TDIAG%XQDIS(:,:) + PSOUT(:,:) * PTSTEP / XRHOLW ![m3]
!
! * Actualisation of input total flux in the river   
!   
IF(LDIAG_MISC)THEN
  TPDG%TDIAG%XQIN(:,:) = TPDG%TDIAG%XQIN (:,:) + PSIN (:,:) * PTSTEP / XRHOLW
ENDIF
!
! * Actualisation of input surface runoff and drainage (or recharge)
!
IF(LCPL_LAND.AND.LDIAG_MISC)THEN
  TPDG%TDIAG%XRUNOFF(:,:) = TPDG%TDIAG%XRUNOFF(:,:) + PRUNOFF(:,:) * PTSTEP / TPG%XAREA(:,:)
  TPDG%TDIAG%XDRAIN (:,:) = TPDG%TDIAG%XDRAIN (:,:) + PDRAIN (:,:) * PTSTEP / TPG%XAREA(:,:)
ENDIF
!
! * Actualisation of stream reservoir
!
TPDG%TDIAG%XSURF_STO(:,:) = TPDG%TDIAG%XSURF_STO(:,:) + TP%XSURF_STO(:,:) * PTSTEP / TPG%XAREA(:,:)
!
! * Actualisation of variable velocity diagnostic variables   
!
IF(CVIT=='VAR')THEN
   TPDG%TDIAG%XVEL(:,:) = TPDG%TDIAG%XVEL(:,:) + PVEL(:,:) * PTSTEP
   TPDG%TDIAG%XHS (:,:) = TPDG%TDIAG%XHS (:,:) + PHS (:,:) * PTSTEP
ENDIF
!
! * Actualisation of groundwater diagnostic variables   
!   
IF(CGROUNDW/='DEF')THEN
  TPDG%TDIAG%XQGF(:,:) = TPDG%TDIAG%XQGF(:,:) + (PGOUT(:,:)+PGNEG(:,:)) * PTSTEP / XRHOLW
ENDIF
!  
IF(CGROUNDW=='CST')THEN        
  TPDG%TDIAG%XGROUND_STO(:,:) = TPDG%TDIAG%XGROUND_STO(:,:) + TP%XGROUND_STO(:,:) * PTSTEP / TPG%XAREA(:,:)
ELSEIF(CGROUNDW=='DIF')THEN
  TPDG%TDIAG%XHGROUND(:,:) = TPDG%TDIAG%XHGROUND(:,:) + TP%XHGROUND (:,:) * PTSTEP
  TPDG%TDIAG%XWTD    (:,:) = TPDG%TDIAG%XWTD    (:,:) + PWTD           (:,:) * PTSTEP
  TPDG%TDIAG%XFWTD   (:,:) = TPDG%TDIAG%XFWTD   (:,:) + PFWTD          (:,:) * PTSTEP
  IF(LDIAG_MISC)THEN
    TPDG%TDIAG%XWTDRIV (:,:) = TPDG%TDIAG%XWTDRIV (:,:) + PWTDRIV (:,:) * PTSTEP
    TPDG%TDIAG%XWTDELEV(:,:) = TPDG%TDIAG%XWTDELEV(:,:) + PWTDELEV(:,:) * PTSTEP
    TPDG%TDIAG%XQGCELL (:,:) = TPDG%TDIAG%XQGCELL (:,:) + PQGCELL (:,:) * PTSTEP / XRHOLW
    TPDG%TDIAG%XHGHS   (:,:) = TPDG%TDIAG%XHGHS   (:,:) + PHGHS   (:,:) * PTSTEP
  ENDIF 
ENDIF
!
! * Actualisation of flooding scheme diagnostic variables   
!
IF(LFLOOD)THEN          
   TPDG%TDIAG%XFLOOD_STO(:,:) = TPDG%TDIAG%XFLOOD_STO(:,:) + TP%XFLOOD_STO(:,:) * PTSTEP / TPG%XAREA(:,:)
   TPDG%TDIAG%XFF       (:,:) = TPDG%TDIAG%XFF       (:,:) + TP%XFFLOOD   (:,:) * PTSTEP
   TPDG%TDIAG%XHF       (:,:) = TPDG%TDIAG%XHF       (:,:) + TP%XHFLOOD   (:,:) * PTSTEP
   IF(LDIAG_MISC)THEN
     TPDG%TDIAG%XQFR   (:,:) = TPDG%TDIAG%XQFR   (:,:) + PQFR            (:,:) * PTSTEP / XRHOLW
     TPDG%TDIAG%XQRF   (:,:) = TPDG%TDIAG%XQRF   (:,:) + PQRF            (:,:) * PTSTEP / XRHOLW
     TPDG%TDIAG%XVFIN  (:,:) = TPDG%TDIAG%XVFIN  (:,:) + PVFIN           (:,:) * PTSTEP
     TPDG%TDIAG%XVFOUT (:,:) = TPDG%TDIAG%XVFOUT (:,:) + PVFOUT          (:,:) * PTSTEP
     TPDG%TDIAG%XWF    (:,:) = TPDG%TDIAG%XWF    (:,:) + TP%XWFLOOD   (:,:) * PTSTEP
     TPDG%TDIAG%XLF    (:,:) = TPDG%TDIAG%XLF    (:,:) + TP%XFLOOD_LEN(:,:) * PTSTEP
     TPDG%TDIAG%XHSF   (:,:) = TPDG%TDIAG%XHSF   (:,:) + PHSF            (:,:) * PTSTEP
     TPDG%TDIAG%XSOURCE(:,:) = TPDG%TDIAG%XSOURCE(:,:) + PSRC_FLOOD      (:,:) * PTSTEP / TPG%XAREA(:,:)
   ENDIF  
ENDIF
!
IF (LHOOK) CALL DR_HOOK('TRIP_DIAG',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE TRIP_DIAG
