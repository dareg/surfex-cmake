!#############################
SUBROUTINE AVERAGE_DIAG_EVAP_ISBA_n
!#############################
!
!
!!****  *AVERAGE_DIAG_EVAP_ISBA_n*  
!!
!!    PURPOSE
!!    -------
!      Average the cumulated diagnostics from all ISBA tiles
!     
!!**  METHOD
!!    ------
!
!!    EXTERNAL
!!    --------
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
!!	P. Le Moigne           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    11/03
!!      B. Decharme 2008     New diag for the water budget
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_n,           ONLY : XPATCH, LGLACIER
USE MODD_DIAG_EVAP_ISBA_n, ONLY : XRNC, XAVG_RNC, XHC, XAVG_HC,                  &
                                  XLEC, XAVG_LEC, XGFLUXC, XAVG_GFLUXC,          &
                                  XLEIC, XAVG_LEIC,                              &
                                  XLEG, XLEGC, XAVG_LEG, XAVG_LEGC,              &
                                  XLEGI, XLEGIC, XAVG_LEGI, XAVG_LEGIC,          &
                                  XLEV, XLEVC, XAVG_LEV, XAVG_LEVC,              &
                                  XLES, XLESC, XAVG_LES, XAVG_LESC,              &
                                  XLER, XLERC, XAVG_LER, XAVG_LERC,              &
                                  XLETR, XLETRC, XAVG_LETR, XAVG_LETRC,          &
                                  XEVAP, XEVAPC, XAVG_EVAP, XAVG_EVAPC,          &
                                  XDRAIN, XDRAINC, XAVG_DRAIN, XAVG_DRAINC,      &
                                  XRUNOFF, XRUNOFFC, XAVG_RUNOFF, XAVG_RUNOFFC,  &
                                  XMELT, XMELTC, XAVG_MELT, XAVG_MELTC,          &
                                  LSURF_EVAP_BUDGET, LSURF_BUDGETC,              &
                                  XHORT, XHORTC, XAVG_HORT, XAVG_HORTC,          &
                                  XDRIP, XDRIPC, XAVG_DRIP, XAVG_DRIPC,          &
                                  XIFLOOD, XIFLOODC, XAVG_IFLOOD, XAVG_IFLOODC,  &
                                  XPFLOOD, XPFLOODC, XAVG_PFLOOD, XAVG_PFLOODC,  &
                                  XLE_FLOOD, XLE_FLOODC, XAVG_LE_FLOOD,          &
                                  XAVG_LE_FLOODC, XLEI_FLOOD, XLEI_FLOODC,       &
                                  XAVG_LEI_FLOOD, XAVG_LEI_FLOODC,               &
                                  XICEFLUXC, XAVG_ICEFLUXC,                      &
                                  XRRVEG, XRRVEGC, XAVG_RRVEG, XAVG_RRVEGC  
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
!
!*      0.2    declarations of local variables
!
INTEGER :: JPATCH ! tile loop counter
INTEGER :: JJ
REAL, DIMENSION(SIZE(XPATCH,1)) :: ZSUMPATCH
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!       0.     Initialization
!              --------------
!
!
IF (LHOOK) CALL DR_HOOK('AVERAGE_DIAG_EVAP_ISBA_N',0,ZHOOK_HANDLE)
ZSUMPATCH(:) = 0.
DO JPATCH=1,SIZE(XPATCH,2)
  ZSUMPATCH(:) = ZSUMPATCH(:) + XPATCH(:,JPATCH)
END DO
!
!       1.     Surface Energy fluxes
!              -----------------------
!
IF (LSURF_EVAP_BUDGET) THEN
   XAVG_LEG(:)   = 0.
   XAVG_LEGI(:)  = 0.
   XAVG_LEV(:)   = 0.
   XAVG_LES(:)   = 0.
   XAVG_LER(:)   = 0.
   XAVG_LETR(:)  = 0.
   XAVG_EVAP(:)  = 0.
   XAVG_DRAIN(:) = 0.
   XAVG_RUNOFF(:)= 0.
   XAVG_HORT(:)  = 0.
   XAVG_DRIP(:)  = 0.
   XAVG_RRVEG(:) = 0.
   XAVG_MELT(:)  = 0.
   XAVG_IFLOOD(:)= 0.
   XAVG_PFLOOD(:)= 0.
   XAVG_LE_FLOOD(:)= 0.
   XAVG_LEI_FLOOD(:)= 0.

  DO JPATCH=1,SIZE(XPATCH,2)
!cdir nodep
    DO JJ=1,SIZE(ZSUMPATCH)
      IF (ZSUMPATCH(JJ) > 0.) THEN
!
! Latent heat of evaporation over the ground
!
        XAVG_LEG(JJ)  = XAVG_LEG(JJ) + XPATCH(JJ,JPATCH) * XLEG(JJ,JPATCH)
!
! Surface soil ice sublimation
!
        XAVG_LEGI(JJ) = XAVG_LEGI(JJ) + XPATCH(JJ,JPATCH) * XLEGI(JJ,JPATCH)
!
! Latent heat of evaporation over vegetation
!
        XAVG_LEV(JJ)  = XAVG_LEV(JJ) + XPATCH(JJ,JPATCH) * XLEV(JJ,JPATCH)
!
! Latent heat of evaporation over snow
!
        XAVG_LES(JJ)  = XAVG_LES(JJ) + XPATCH(JJ,JPATCH) * XLES(JJ,JPATCH)
!
! Evaporation from canopy water interception
!
        XAVG_LER(JJ)  = XAVG_LER(JJ) + XPATCH(JJ,JPATCH) * XLER(JJ,JPATCH)
!
! Evapotranspiration of the vegetation
!
        XAVG_LETR(JJ)  = XAVG_LETR(JJ) + XPATCH(JJ,JPATCH) * XLETR(JJ,JPATCH)
!
! Evapotranspiration
!
        XAVG_EVAP(JJ)  = XAVG_EVAP(JJ) + XPATCH(JJ,JPATCH) * XEVAP(JJ,JPATCH)
!
! Soil drainage flux
!
        XAVG_DRAIN(JJ)  = XAVG_DRAIN(JJ) + XPATCH(JJ,JPATCH) * XDRAIN(JJ,JPATCH)
!
! Supersaturation runoff
!
        XAVG_RUNOFF(JJ) = XAVG_RUNOFF(JJ) + XPATCH(JJ,JPATCH) * XRUNOFF(JJ,JPATCH)
!
! Horton runoff
!
        XAVG_HORT(JJ)  = XAVG_HORT(JJ) + XPATCH(JJ,JPATCH) * XHORT(JJ,JPATCH)
!
! Vegetation dripping
!
        XAVG_DRIP(JJ)  = XAVG_DRIP(JJ) + XPATCH(JJ,JPATCH) * XDRIP(JJ,JPATCH)
!
! Precipitation intercepted by the	vegetation
!
        XAVG_RRVEG(JJ)  = XAVG_RRVEG(JJ) + XPATCH(JJ,JPATCH) * XRRVEG(JJ,JPATCH)
!      
! Snow melt
!
        XAVG_MELT(JJ)  = XAVG_MELT(JJ) + XPATCH(JJ,JPATCH) * XMELT(JJ,JPATCH)
!      
! Flood infiltartion
!
        XAVG_IFLOOD(JJ) = XAVG_IFLOOD(JJ) + XPATCH(JJ,JPATCH) * XIFLOOD(JJ,JPATCH)
!      
! Precipitation intercepted by the floodplains
!     
        XAVG_PFLOOD(JJ) = XAVG_PFLOOD(JJ) + XPATCH(JJ,JPATCH) * XPFLOOD(JJ,JPATCH)
!      
! Floodplains evaporation
!     
        XAVG_LE_FLOOD (JJ) = XAVG_LE_FLOOD (JJ) + XPATCH(JJ,JPATCH) * XLE_FLOOD (JJ,JPATCH)
        XAVG_LEI_FLOOD(JJ) = XAVG_LEI_FLOOD(JJ) + XPATCH(JJ,JPATCH) * XLEI_FLOOD(JJ,JPATCH)
!
      ENDIF
    END DO
  ENDDO
END IF
!
!
!       2.     Surface Cumulated Energy fluxes
!              -------------------------------
!
IF (LSURF_BUDGETC) THEN
   XAVG_RNC(:)     = 0.
   XAVG_HC(:)      = 0.
   XAVG_LEC(:)     = 0.
   XAVG_GFLUXC(:)  = 0.
   XAVG_LEIC(:)    = 0.
   XAVG_LEGC(:)    = 0.
   XAVG_LEGIC(:)   = 0.
   XAVG_LEVC(:)    = 0.
   XAVG_LESC(:)    = 0.
   XAVG_LERC(:)    = 0.
   XAVG_LETRC(:)   = 0.
   XAVG_EVAPC(:)   = 0.
   XAVG_DRAINC(:)  = 0.
   XAVG_RUNOFFC(:) = 0.
   XAVG_HORTC(:)   = 0.
   XAVG_DRIPC(:)   = 0.
   XAVG_RRVEGC(:)  = 0.
   XAVG_MELTC(:)   = 0.  
   XAVG_IFLOODC(:) = 0.
   XAVG_PFLOODC(:) = 0.
   XAVG_LE_FLOODC(:)= 0.
   XAVG_LEI_FLOODC(:)= 0.

  DO JPATCH=1,SIZE(XPATCH,2)
!cdir nodep
    DO JJ=1,SIZE(ZSUMPATCH)
      IF (ZSUMPATCH(JJ) > 0.) THEN
!
! Net radiation
!
        XAVG_RNC(JJ)  = XAVG_RNC(JJ) + XPATCH(JJ,JPATCH) * XRNC(JJ,JPATCH)
!
! Sensible heat flux
!
        XAVG_HC(JJ)  = XAVG_HC(JJ) + XPATCH(JJ,JPATCH) * XHC(JJ,JPATCH)
!
! Total latent heat flux
!
        XAVG_LEC(JJ)  = XAVG_LEC(JJ) + XPATCH(JJ,JPATCH) * XLEC(JJ,JPATCH)
!
! Storage flux
!
        XAVG_GFLUXC(JJ)  = XAVG_GFLUXC(JJ) + XPATCH(JJ,JPATCH) * XGFLUXC(JJ,JPATCH)
!
! Total surface sublimation
!
        XAVG_LEIC(JJ)  = XAVG_LEIC(JJ) + XPATCH(JJ,JPATCH) * XLEIC(JJ,JPATCH)
!
! Latent heat of evaporation over the ground
!
        XAVG_LEGC(JJ)  = XAVG_LEGC(JJ) + XPATCH(JJ,JPATCH) * XLEGC(JJ,JPATCH)
!
! Surface soil ice sublimation
!
        XAVG_LEGIC(JJ)  = XAVG_LEGIC(JJ) + XPATCH(JJ,JPATCH) * XLEGIC(JJ,JPATCH)
!
! Latent heat of evaporation over vegetation
!
        XAVG_LEVC(JJ)  = XAVG_LEVC(JJ) + XPATCH(JJ,JPATCH) * XLEVC(JJ,JPATCH)
!
! Latent heat of evaporation over snow
!
        XAVG_LESC(JJ)  = XAVG_LESC(JJ) + XPATCH(JJ,JPATCH) * XLESC(JJ,JPATCH)
!
! Evaporation from canopy water interception
!
        XAVG_LERC(JJ)  = XAVG_LERC(JJ) + XPATCH(JJ,JPATCH) * XLERC(JJ,JPATCH)
!
! Evapotranspiration of the vegetation
!
        XAVG_LETRC(JJ)  = XAVG_LETRC(JJ) + XPATCH(JJ,JPATCH) * XLETRC(JJ,JPATCH)
!
! Evapotranspiration
!
        XAVG_EVAPC(JJ)  = XAVG_EVAPC(JJ) + XPATCH(JJ,JPATCH) * XEVAPC(JJ,JPATCH)
!
! Soil drainage flux
!
        XAVG_DRAINC(JJ)  = XAVG_DRAINC(JJ) + XPATCH(JJ,JPATCH) * XDRAINC(JJ,JPATCH)
!
! Supersaturation runoff
!
        XAVG_RUNOFFC(JJ)  = XAVG_RUNOFFC(JJ) + XPATCH(JJ,JPATCH) * XRUNOFFC(JJ,JPATCH)
!
! Horton runoff
!
        XAVG_HORTC(JJ)  = XAVG_HORTC(JJ) + XPATCH(JJ,JPATCH) * XHORTC(JJ,JPATCH)
!
! Vegetation dripping
!
        XAVG_DRIPC(JJ)  = XAVG_DRIPC(JJ) + XPATCH(JJ,JPATCH) * XDRIPC(JJ,JPATCH)
!
! precipitation intercepted by the	vegetation
!
        XAVG_RRVEGC(JJ)  = XAVG_RRVEGC(JJ) + XPATCH(JJ,JPATCH) * XRRVEGC(JJ,JPATCH)
!      
! Snow melt
!
        XAVG_MELTC(JJ)  = XAVG_MELTC(JJ) + XPATCH(JJ,JPATCH) * XMELTC(JJ,JPATCH)
!      
! Flood infiltartion
!
        XAVG_IFLOODC(JJ) = XAVG_IFLOODC(JJ) + XPATCH(JJ,JPATCH) * XIFLOODC(JJ,JPATCH)
!      
! Precipitation intercepted by the floodplains
!     
        XAVG_PFLOODC(JJ) = XAVG_PFLOODC(JJ) + XPATCH(JJ,JPATCH) * XPFLOODC(JJ,JPATCH)
!      
! Floodplains evaporation
!     
        XAVG_LE_FLOODC (JJ) = XAVG_LE_FLOODC (JJ) + XPATCH(JJ,JPATCH) * XLE_FLOODC (JJ,JPATCH)
        XAVG_LEI_FLOODC(JJ) = XAVG_LEI_FLOODC(JJ) + XPATCH(JJ,JPATCH) * XLEI_FLOODC(JJ,JPATCH)
!      
      ENDIF
    ENDDO
  END DO
END IF
!
IF(LSURF_BUDGETC.AND.LGLACIER)THEN 
  XAVG_ICEFLUXC(:)= 0.
  DO JPATCH=1,SIZE(XPATCH,2)
    WHERE (ZSUMPATCH(:) > 0.)
      XAVG_ICEFLUXC(:) = XAVG_ICEFLUXC(:) + XPATCH(:,JPATCH) * XICEFLUXC(:,JPATCH)      
    END WHERE
  END DO
END IF
!        
!
IF (LHOOK) CALL DR_HOOK('AVERAGE_DIAG_EVAP_ISBA_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE AVERAGE_DIAG_EVAP_ISBA_n
