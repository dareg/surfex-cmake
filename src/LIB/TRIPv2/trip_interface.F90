!#################################################################
SUBROUTINE TRIP_INTERFACE (TPDG, TP, TPG, &
                            KLISTING,KLON,KLAT,PTIME,OPRINT, &
                           KNB_TSTEP_RUN,KNB_TSTEP_DIAG,    &
                           PTSTEP_RUN,PTSTEP_DIAG,PRUNOFF,  &
                           PDRAIN,PCALVING,PRECHARGE,       &
                           PSRC_FLOOD                       ) 
!#################################################################
!
!!****  *TRIP*  
!!
!!    PURPOSE
!!    -------
!
!     Driver for the TRIP river routing.
!     Here, we call the physical and the diag routines     
!
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
!!      Original    01/02/05 
!!      Modif.      28/05/08 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!
!
USE MODD_TRIP_DIAG, ONLY : TRIP_DIAG_t
USE MODD_TRIP, ONLY : TRIP_t
USE MODD_TRIP_GRID, ONLY : TRIP_GRID_t
!
USE MODN_TRIP_RUN,  ONLY : LDIAG_MISC
!
USE MODN_TRIP,      ONLY : CGROUNDW, LFLOOD, XTSTEP
!
USE MODD_TRIP_PAR,  ONLY : XRHOLW, XUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_TRIP
USE MODI_TRIP_DIAG_INIT
USE MODI_TRIP_DIAG
USE MODI_TRIP_DIAG_WRITE
USE MODI_TRIP_DIAG_CPL_ESM
!
USE MODI_GWF
USE MODI_GWF_CPL_UPDATE
!
USE MODI_ABORT_TRIP
!
USE MODI_GWF
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
!
TYPE(TRIP_DIAG_t), INTENT(INOUT) :: TPDG
TYPE(TRIP_t), INTENT(INOUT) :: TP
TYPE(TRIP_GRID_t), INTENT(INOUT) :: TPG
!
INTEGER,              INTENT(IN)    :: KLISTING       !Output file id
INTEGER,              INTENT(IN)    :: KLON           !Number of longitude
INTEGER,              INTENT(IN)    :: KLAT           !Number of latittude
REAL,                 INTENT(INOUT) :: PTIME          !Current time          (s)
LOGICAL,              INTENT(IN)    :: OPRINT         !print option          [-]
INTEGER,              INTENT(IN)    :: KNB_TSTEP_RUN  !TSTEP_RUN counter     [-]
REAL,                 INTENT(IN)    :: PTSTEP_RUN     !Run  timestep         [s]
REAL,                 INTENT(IN)    :: PTSTEP_DIAG    !Diag timestep         [s]
INTEGER,              INTENT(INOUT) :: KNB_TSTEP_DIAG !DIAG call counter     [-]
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PRUNOFF   !Input surface runoff            [kg/s]
REAL, DIMENSION(:,:), INTENT(IN)    :: PDRAIN    !Input free drainage             [kg/s]
REAL, DIMENSION(:,:), INTENT(IN)    :: PCALVING  !Input claving flux from glacier [kg/s]
REAL, DIMENSION(:,:), INTENT(IN)    :: PRECHARGE !Input goundwater recharge       [kg/s]
REAL, DIMENSION(:,:), INTENT(IN)    :: PSRC_FLOOD! Input P-E-I flood source term  [kg/s]
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(KLON,KLAT) :: ZRUNOFF    !Input surface runoff          [kg/s]
REAL, DIMENSION(KLON,KLAT) :: ZDRAIN     !Input drainage + recharge     [kg/s]
REAL, DIMENSION(KLON,KLAT) :: ZSRC_FLOOD !Input P-E-I flood source term [kg/s]
!
REAL, DIMENSION(KLON,KLAT) :: ZSOUT      !streamflow                    [kg/s]
REAL, DIMENSION(KLON,KLAT) :: ZSIN       !grid-cell input streamflow    [kg/s]
REAL, DIMENSION(KLON,KLAT) :: ZVEL       !river velocity                [m/s]
REAL, DIMENSION(KLON,KLAT) :: ZHS        !River heigh                   [m]
REAL, DIMENSION(KLON,KLAT) :: ZGOUT      !Groundwater outflow           [kg/s]
REAL, DIMENSION(KLON,KLAT) :: ZGNEG      !Groundwater intflow (neg)     [kg/s]
REAL, DIMENSION(KLON,KLAT) :: ZWTD       !Water table depth for coupling[m]
REAL, DIMENSION(KLON,KLAT) :: ZFWTD      !fraction of water table to rise
REAL, DIMENSION(KLON,KLAT) :: ZWTDRIV    !Water table depth / topo riv  [m]
REAL, DIMENSION(KLON,KLAT) :: ZWTDELEV   !Water table depth / elevation [m]
REAL, DIMENSION(KLON,KLAT) :: ZQGCELL    !lateral groundwater exchanges [kg/s]
REAL, DIMENSION(KLON,KLAT) :: ZHGHS      !groundwater minus river heigh [m]
REAL, DIMENSION(KLON,KLAT) :: ZQFR       !floodplains to river exchange [kg/s]
REAL, DIMENSION(KLON,KLAT) :: ZQRF       !river to floodplains exchange [kg/s]
REAL, DIMENSION(KLON,KLAT) :: ZVFIN      !QRF velocity                  [m/s]
REAL, DIMENSION(KLON,KLAT) :: ZVFOUT     !QFR velocity                  [m/s]
REAL, DIMENSION(KLON,KLAT) :: ZHSF       !river minus flodd heigh       [m]
REAL, DIMENSION(KLON,KLAT) :: ZDISCHARGE !river discharges              [kg]
REAL, DIMENSION(KLON,KLAT) :: ZHG_OLD    !Water table depth at t-1      [m]
!
REAL                       :: ZGSTO_ALL  !Global groundwater storage at t    [kg]
REAL                       :: ZGSTO2_ALL !Global groundwater storage at t-1  [kg]
REAL                       :: ZGIN_ALL   !Global gw recharge + lateral input [kg/m2/s]
REAL                       :: ZGOUT_ALL  !Global gw outflow                  [kg/m2/s]
!
INTEGER :: JTSTEP, ITSTEP
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TRIP_INTERFACE',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!*       1.     Initialize :
!               ------------
!
!Surface runoff treatment
!
ZRUNOFF(:,:) = PRUNOFF(:,:)
!
!Drainage and Calving treatment 
!calving over greenland and antarctica directly to ocean
!
WHERE(TPG%GMASK(:,:).AND..NOT.TPG%GMASK_GRE(:,:).AND..NOT.TPG%GMASK_ANT(:,:))
  ZDRAIN(:,:) = PDRAIN(:,:)+PRECHARGE(:,:)+PCALVING(:,:)
ELSEWHERE
  ZDRAIN(:,:) = PDRAIN(:,:)+PRECHARGE(:,:)
ENDWHERE
!
! Flood treatment
!
IF(LFLOOD)THEN 
  ZSRC_FLOOD(:,:) = PSRC_FLOOD(:,:)
  WHERE(TP%XFFLOOD(:,:)==1.0.AND.ZSRC_FLOOD(:,:)>0.0)
        ZRUNOFF   (:,:) = ZRUNOFF(:,:) + PSRC_FLOOD(:,:)
        ZSRC_FLOOD(:,:) = 0.0
  ENDWHERE
ELSE
  ZSRC_FLOOD(:,:) = XUNDEF  
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       2.     Initialize local diag :
!               -----------------------
!
CALL TRIP_DIAG_INIT(ZSOUT,ZSIN,ZVEL,ZHS,ZGOUT,ZGNEG,ZHG_OLD,   &
                    ZWTD,ZFWTD,ZWTDRIV,ZWTDELEV,ZQGCELL,ZHGHS, &
                    ZQFR,ZQRF,ZVFIN,ZVFOUT,ZHSF,ZSRC_FLOOD,    &
                    ZDISCHARGE,                                &
                    ZGSTO_ALL,ZGSTO2_ALL,ZGIN_ALL,ZGOUT_ALL    )
!
!-------------------------------------------------------------------------------
!
!*       3.     Call Groundwater dynamic :
!               --------------------------
!
IF(CGROUNDW=='DIF')THEN
  CALL GWF(TPG, &
           KLON,KLAT,OPRINT,PTSTEP_RUN,XTSTEP,          &
           TPG%GMASK_GW,TP%XNUM_AQUI,ZDRAIN,       &
           TPG%XLEN,TP%XWIDTH,TP%XHC_BED,       &
           TP%XTOPO_RIV,TP%XTAUG,TPG%XAREA,     &
           TP%XELEV,TP%XTRANS,TP%XWEFF,        &
           TP%XTABGW_F,TP%XTABGW_H,               &
           TP%XSURF_STO,TP%XHGROUND,ZHG_OLD,      &
           ZQGCELL,ZWTD,ZFWTD,ZWTDRIV,ZWTDELEV,         &
           ZHGHS,ZGOUT,ZGNEG,                           &
           ZGSTO_ALL,ZGSTO2_ALL,ZGIN_ALL,ZGOUT_ALL      )
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       4.     Call Trip river routines and actualisation of diagnostic :
!               ----------------------------------------------------------
!
ITSTEP = INT(PTSTEP_RUN/XTSTEP)
!
DO JTSTEP=1,ITSTEP !TRIP time step loop
!
  CALL TRIP(KLISTING,CGROUNDW,LFLOOD,OPRINT,XTSTEP,               &
            TPG%NGRCN,TPG%NSEQ,TPG%NNEXTX,TPG%NNEXTY,     &
            TPG%NSEQMAX,TPG%XAREA,TPG%XLEN,                 &
            TPG%GMASK_GW,TPG%GMASK_VEL,TPG%GMASK_FLD,       &
            TP%XTAUG,TP%XFLOOD_LEN,TP%XSLOPEBED,         &
            TP%XWIDTH,TP%XN,TP%XN_FLOOD,TP%XHC_BED,   &
            TP%XWFLOOD,TP%XTAB_F,TP%XTAB_H,              &
            TP%XTAB_VF,ZDRAIN,ZRUNOFF,ZSRC_FLOOD,              &
            TP%XGROUND_STO,TP%XSURF_STO,TP%XFLOOD_STO,   &
            ZSOUT,ZGOUT,ZHS,TP%XHFLOOD,ZVEL,TP%XFFLOOD,     &
            ZQFR,ZQRF,ZVFIN,ZVFOUT,ZHSF,ZSIN,KNB_TSTEP_RUN,       &
            JTSTEP,ITSTEP,ZGSTO_ALL,ZGSTO2_ALL,ZGIN_ALL,ZGOUT_ALL,&
            TP%XHGROUND,TP%XWEFF)
!
!  * Actualisation of diagnostic  
!
   IF(CGROUNDW=='DIF')THEN
      CALL GWF_CPL_UPDATE(TP%XTABGW_H,TP%XTABGW_F,TPG%GMASK_GW, &
                          TP%XTOPO_RIV,TP%XELEV,TP%XHGROUND,   &
                          ZHG_OLD,ZWTD,ZFWTD,ZWTDELEV )          
   ENDIF
!
   CALL TRIP_DIAG(TPDG, TP, TPG, &
                  XTSTEP,ZSOUT,ZSIN,ZVEL,ZHS,ZGOUT,ZGNEG,    &
                  ZWTD,ZFWTD,ZWTDRIV,ZWTDELEV,ZQGCELL,ZHGHS, &
                  ZQFR,ZQRF,ZVFIN,ZVFOUT,ZHSF,ZSRC_FLOOD,    &
                  ZDRAIN,ZRUNOFF,ZDISCHARGE                  )
!
!  * Time actualization  
!
   PTIME = PTIME + XTSTEP
!
!  * Write diagnostic  
!
   IF (MOD(PTIME,PTSTEP_DIAG) == 0. ) THEN
      KNB_TSTEP_DIAG = KNB_TSTEP_DIAG + 1
      CALL TRIP_DIAG_WRITE(TPDG, TPG, &
                           KLISTING,KLON,KLAT,KNB_TSTEP_DIAG,PTSTEP_DIAG)
   ENDIF
!
!  * end 
!
   IF(OPRINT)WRITE(KLISTING,*)' '
!
ENDDO ! * End TRIP time step loop
!
!
!-------------------------------------------------------------------------------
!
!*       5.      Actualisation of coupling diagnostic:
!               --------------------------------------
!
CALL TRIP_DIAG_CPL_ESM(TP, TPG, &
                       PTSTEP_RUN,ZDISCHARGE,PCALVING,ZWTD,ZFWTD)
!
IF (LHOOK) CALL DR_HOOK('TRIP_INTERFACE',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
END SUBROUTINE TRIP_INTERFACE
