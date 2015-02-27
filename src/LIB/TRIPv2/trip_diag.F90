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
USE MODD_TRIP_GRID,  ONLY : TGRID
USE MODD_TRIP,       ONLY : TTRIP
USE MODD_TRIP_DIAG,  ONLY : TDIAG
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
TDIAG%XQDIS(:,:) = TDIAG%XQDIS(:,:) + PSOUT(:,:) * PTSTEP / XRHOLW ![m3]
!
! * Actualisation of input total flux in the river   
!   
IF(LDIAG_MISC)THEN
  TDIAG%XQIN(:,:) = TDIAG%XQIN (:,:) + PSIN (:,:) * PTSTEP / XRHOLW
ENDIF
!
! * Actualisation of input surface runoff and drainage (or recharge)
!
IF(LCPL_LAND.AND.LDIAG_MISC)THEN
  TDIAG%XRUNOFF(:,:) = TDIAG%XRUNOFF(:,:) + PRUNOFF(:,:) * PTSTEP / TGRID%XAREA(:,:)
  TDIAG%XDRAIN (:,:) = TDIAG%XDRAIN (:,:) + PDRAIN (:,:) * PTSTEP / TGRID%XAREA(:,:)
ENDIF
!
! * Actualisation of stream reservoir
!
TDIAG%XSURF_STO(:,:) = TDIAG%XSURF_STO(:,:) + TTRIP%XSURF_STO(:,:) * PTSTEP / TGRID%XAREA(:,:)
!
! * Actualisation of variable velocity diagnostic variables   
!
IF(CVIT=='VAR')THEN
   TDIAG%XVEL(:,:) = TDIAG%XVEL(:,:) + PVEL(:,:) * PTSTEP
   TDIAG%XHS (:,:) = TDIAG%XHS (:,:) + PHS (:,:) * PTSTEP
ENDIF
!
! * Actualisation of groundwater diagnostic variables   
!   
IF(CGROUNDW/='DEF')THEN
  TDIAG%XQGF(:,:) = TDIAG%XQGF(:,:) + (PGOUT(:,:)+PGNEG(:,:)) * PTSTEP / XRHOLW
ENDIF
!  
IF(CGROUNDW=='CST')THEN        
  TDIAG%XGROUND_STO(:,:) = TDIAG%XGROUND_STO(:,:) + TTRIP%XGROUND_STO(:,:) * PTSTEP / TGRID%XAREA(:,:)
ELSEIF(CGROUNDW=='DIF')THEN
  TDIAG%XHGROUND(:,:) = TDIAG%XHGROUND(:,:) + TTRIP%XHGROUND (:,:) * PTSTEP
  TDIAG%XWTD    (:,:) = TDIAG%XWTD    (:,:) + PWTD           (:,:) * PTSTEP
  TDIAG%XFWTD   (:,:) = TDIAG%XFWTD   (:,:) + PFWTD          (:,:) * PTSTEP
  IF(LDIAG_MISC)THEN
    TDIAG%XWTDRIV (:,:) = TDIAG%XWTDRIV (:,:) + PWTDRIV (:,:) * PTSTEP
    TDIAG%XWTDELEV(:,:) = TDIAG%XWTDELEV(:,:) + PWTDELEV(:,:) * PTSTEP
    TDIAG%XQGCELL (:,:) = TDIAG%XQGCELL (:,:) + PQGCELL (:,:) * PTSTEP / XRHOLW
    TDIAG%XHGHS   (:,:) = TDIAG%XHGHS   (:,:) + PHGHS   (:,:) * PTSTEP
  ENDIF 
ENDIF
!
! * Actualisation of flooding scheme diagnostic variables   
!
IF(LFLOOD)THEN          
   TDIAG%XFLOOD_STO(:,:) = TDIAG%XFLOOD_STO(:,:) + TTRIP%XFLOOD_STO(:,:) * PTSTEP / TGRID%XAREA(:,:)
   TDIAG%XFF       (:,:) = TDIAG%XFF       (:,:) + TTRIP%XFFLOOD   (:,:) * PTSTEP
   TDIAG%XHF       (:,:) = TDIAG%XHF       (:,:) + TTRIP%XHFLOOD   (:,:) * PTSTEP
   IF(LDIAG_MISC)THEN
     TDIAG%XQFR   (:,:) = TDIAG%XQFR   (:,:) + PQFR            (:,:) * PTSTEP / XRHOLW
     TDIAG%XQRF   (:,:) = TDIAG%XQRF   (:,:) + PQRF            (:,:) * PTSTEP / XRHOLW
     TDIAG%XVFIN  (:,:) = TDIAG%XVFIN  (:,:) + PVFIN           (:,:) * PTSTEP
     TDIAG%XVFOUT (:,:) = TDIAG%XVFOUT (:,:) + PVFOUT          (:,:) * PTSTEP
     TDIAG%XWF    (:,:) = TDIAG%XWF    (:,:) + TTRIP%XWFLOOD   (:,:) * PTSTEP
     TDIAG%XLF    (:,:) = TDIAG%XLF    (:,:) + TTRIP%XFLOOD_LEN(:,:) * PTSTEP
     TDIAG%XHSF   (:,:) = TDIAG%XHSF   (:,:) + PHSF            (:,:) * PTSTEP
     TDIAG%XSOURCE(:,:) = TDIAG%XSOURCE(:,:) + PSRC_FLOOD      (:,:) * PTSTEP / TGRID%XAREA(:,:)
   ENDIF  
ENDIF
!
IF (LHOOK) CALL DR_HOOK('TRIP_DIAG',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE TRIP_DIAG
