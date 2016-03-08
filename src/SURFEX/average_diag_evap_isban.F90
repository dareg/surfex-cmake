!#############################
SUBROUTINE AVERAGE_DIAG_EVAP_ISBA_n (OSURF_BUDGETC, DGEI, DGEIC, DGEIP, DGEIPC, &
                                     OGLACIER, OMEB_PATCH, PPATCH, PTSTEP, PRAIN, PSNOW)
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
!!      P. Le Moigne           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    11/03
!!      B. Decharme 2008     New diag for the water budget
!!      B. Decharme 2012     New diag for snow 
!!                                        carbon
!!                                        isab water budget
!!                  2013                  Sublimation
!!                                        Subsurface runoff if SGH (DIF option only)
!!      P. Samuelsson 10/2014: MEB
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!
!
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t, DIAG_EVAP_ISBA_PATCH_t
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
LOGICAL, INTENT(IN) :: OSURF_BUDGETC
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEI
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEIC 
TYPE(DIAG_EVAP_ISBA_PATCH_t), INTENT(INOUT) :: DGEIP
TYPE(DIAG_EVAP_ISBA_PATCH_t), INTENT(INOUT) :: DGEIPC

!
LOGICAL, INTENT(IN) :: OGLACIER
LOGICAL, DIMENSION(:), INTENT(IN) :: OMEB_PATCH
REAL, DIMENSION(:,:), INTENT(IN) :: PPATCH
!
REAL,                  INTENT(IN) :: PTSTEP        ! time step (s)
REAL,    DIMENSION(:), INTENT(IN) :: PRAIN         ! rainfall rate
REAL,    DIMENSION(:), INTENT(IN) :: PSNOW         ! snowfall rate
!
!
!*      0.2    declarations of local variables
!
INTEGER :: JPATCH ! tile loop counter
INTEGER :: JJ
INTEGER           :: ISIZE_LMEB_PATCH   ! Number of patches where multi-energy balance should be applied
REAL, DIMENSION(SIZE(PPATCH,1)) :: ZSUMPATCH
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!       0.     Initialization
!              --------------
!
!
IF (LHOOK) CALL DR_HOOK('AVERAGE_DIAG_EVAP_ISBA_N',0,ZHOOK_HANDLE)
!
ISIZE_LMEB_PATCH=COUNT(OMEB_PATCH(:))
!
ZSUMPATCH(:) = 0.
DO JPATCH=1,SIZE(PPATCH,2)
   DO JJ=1,SIZE(PPATCH,1)
      ZSUMPATCH(JJ) = ZSUMPATCH(JJ) + PPATCH(JJ,JPATCH)
  ENDDO
ENDDO
!
!       1.     Surface Energy fluxes
!              -----------------------
!
IF (DGEI%LSURF_EVAP_BUDGET) THEN
!        
   DGEI%XLEG        (:) = 0.
   DGEI%XLEGI       (:) = 0.
   DGEI%XLEV        (:) = 0.
   DGEI%XLES        (:) = 0.
   DGEI%XLESL       (:) = 0.
   DGEI%XLER        (:) = 0.
   DGEI%XLETR       (:) = 0.
   DGEI%XSNDRIFT    (:) = 0.
   DGEI%XDRAIN      (:) = 0.
   DGEI%XQSB        (:) = 0.
   DGEI%XRUNOFF     (:) = 0.
   DGEI%XHORT       (:) = 0.
   DGEI%XDRIP       (:) = 0.
   DGEI%XRRVEG      (:) = 0.
   DGEI%XMELT       (:) = 0.
   DGEI%XIFLOOD     (:) = 0.
   DGEI%XPFLOOD     (:) = 0.
   DGEI%XLE_FLOOD   (:) = 0.
   DGEI%XLEI_FLOOD  (:) = 0.
   DGEI%XIRRIG_FLUX (:) = 0.
   DGEI%XGPP        (:) = 0.
   DGEI%XRESP_AUTO  (:) = 0.
   DGEI%XRESP_ECO   (:) = 0.
!
   IF (ISIZE_LMEB_PATCH>0) THEN
     DGEI%XLEVCV         (:) = 0.
     DGEI%XLESC          (:) = 0.
     DGEI%XLETRGV        (:) = 0.
     DGEI%XLETRCV        (:) = 0.
     DGEI%XLERGV         (:) = 0.
     DGEI%XLELITTER      (:) = 0.
     DGEI%XLELITTERI     (:) = 0.
     DGEI%XDRIPLIT       (:) = 0.
     DGEI%XRRLIT         (:) = 0.
     DGEI%XLERCV         (:) = 0.
     DGEI%XLE_C_A        (:) = 0.
     DGEI%XLE_V_C        (:) = 0.
     DGEI%XLE_G_C        (:) = 0.
     DGEI%XLE_N_C        (:) = 0.
     !
     DGEI%XSWNET_V       (:) = 0.
     DGEI%XSWNET_G       (:) = 0.
     DGEI%XSWNET_N       (:) = 0.
     DGEI%XSWNET_NS      (:) = 0.
     DGEI%XLWNET_V       (:) = 0.
     DGEI%XLWNET_G       (:) = 0.
     DGEI%XLWNET_N       (:) = 0.
     DGEI%XSWDOWN_GN     (:) = 0.
     DGEI%XLWDOWN_GN     (:) = 0.
     DGEI%XH_V_C         (:) = 0.
     DGEI%XH_G_C         (:) = 0.
     DGEI%XH_C_A         (:) = 0.
     DGEI%XH_N_C         (:) = 0.
     DGEI%XSR_GN         (:) = 0.
     DGEI%XMELTCV        (:) = 0.
     DGEI%XFRZCV         (:) = 0.
   ENDIF
!
  DO JPATCH=1,SIZE(PPATCH,2)
!cdir nodep
    DO JJ=1,SIZE(ZSUMPATCH)
      IF (ZSUMPATCH(JJ) > 0.) THEN
!
! Latent heat of evaporation over the ground
!
        DGEI%XLEG(JJ)  = DGEI%XLEG(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XLEG(JJ)
!
! Surface soil ice sublimation
!
        DGEI%XLEGI(JJ) = DGEI%XLEGI(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XLEGI(JJ)
!
! Latent heat of evaporation over vegetation
!
        DGEI%XLEV(JJ)  = DGEI%XLEV(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XLEV(JJ)
!
! Latent heat of sublimation over snow
!
        DGEI%XLES(JJ)  = DGEI%XLES(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XLES(JJ)
!
! Latent heat of evaporation of liquid water over snow
!
        DGEI%XLESL(JJ)  = DGEI%XLESL(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XLESL(JJ)
!
! Evaporation from canopy water interception
!
        DGEI%XLER(JJ)  = DGEI%XLER(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XLER(JJ)
!
! Evapotranspiration of the vegetation
!
        DGEI%XLETR(JJ)  = DGEI%XLETR(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XLETR(JJ)
!
! Blowing snow sublimation (ES or Crocus)
!
        DGEI%XSNDRIFT(JJ)  = DGEI%XSNDRIFT(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XSNDRIFT(JJ)
!
! Soil drainage flux
!
        DGEI%XDRAIN(JJ)  = DGEI%XDRAIN(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XDRAIN(JJ)
!
! Soil lateral subsurface flux
!
        DGEI%XQSB(JJ)  = DGEI%XQSB(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XQSB(JJ)        
!
! Supersaturation runoff
!
        DGEI%XRUNOFF(JJ) = DGEI%XRUNOFF(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XRUNOFF(JJ)
!
! Horton runoff
!
        DGEI%XHORT(JJ)  = DGEI%XHORT(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XHORT(JJ)
!
! Vegetation dripping
!
        DGEI%XDRIP(JJ)  = DGEI%XDRIP(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XDRIP(JJ)
!
! Precipitation intercepted by the vegetation
!
        DGEI%XRRVEG(JJ)  = DGEI%XRRVEG(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XRRVEG(JJ)
!      
! Snow melt
!
        DGEI%XMELT(JJ)  = DGEI%XMELT(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XMELT(JJ)
!      
! Flood infiltartion
!
        DGEI%XIFLOOD(JJ) = DGEI%XIFLOOD(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XIFLOOD(JJ)
!      
! Precipitation intercepted by the floodplains
!     
        DGEI%XPFLOOD(JJ) = DGEI%XPFLOOD(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XPFLOOD(JJ)
!      
! Floodplains evaporation
!     
        DGEI%XLE_FLOOD (JJ) = DGEI%XLE_FLOOD (JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XLE_FLOOD (JJ)
        DGEI%XLEI_FLOOD(JJ) = DGEI%XLEI_FLOOD(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XLEI_FLOOD(JJ)
!      
! irrigation rate (as soil input)
!
        DGEI%XIRRIG_FLUX(JJ)  = DGEI%XIRRIG_FLUX(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XIRRIG_FLUX(JJ)
!
! Gross primary production
!
        DGEI%XGPP(JJ) = DGEI%XGPP(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XGPP(JJ)
!
! Autotrophic respiration
!   
        DGEI%XRESP_AUTO(JJ) = DGEI%XRESP_AUTO(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XRESP_AUTO(JJ)
!
! Ecosystem respiration
!
        DGEI%XRESP_ECO(JJ) = DGEI%XRESP_ECO(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XRESP_ECO(JJ)  
!        
        IF (ISIZE_LMEB_PATCH>0) THEN
          DGEI%XLEVCV(JJ) = DGEI%XLEVCV(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XLEVCV(JJ)
          DGEI%XLESC(JJ) = DGEI%XLESC(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XLESC(JJ)
          DGEI%XLETRCV(JJ) = DGEI%XLETRCV(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XLETRCV(JJ)
          DGEI%XLELITTER(JJ) = DGEI%XLELITTER(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XLELITTER(JJ)
          DGEI%XLELITTERI(JJ) = DGEI%XLELITTERI(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XLELITTERI(JJ)
          DGEI%XDRIPLIT(JJ) = DGEI%XDRIPLIT(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XDRIPLIT(JJ)
          DGEI%XRRLIT(JJ) = DGEI%XRRLIT(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XRRLIT(JJ)
          DGEI%XLERCV(JJ) = DGEI%XLERCV(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XLERCV(JJ)
          DGEI%XLE_C_A(JJ) = DGEI%XLE_C_A(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XLE_C_A(JJ)
          DGEI%XLE_V_C(JJ) = DGEI%XLE_V_C(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XLE_V_C(JJ)
          DGEI%XLE_G_C(JJ) = DGEI%XLE_G_C(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XLE_G_C(JJ)
          DGEI%XLE_N_C(JJ) = DGEI%XLE_N_C(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XLE_N_C(JJ)
          DGEI%XSWNET_V(JJ) = DGEI%XSWNET_V(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XSWNET_V(JJ)
          DGEI%XSWNET_G(JJ) = DGEI%XSWNET_G(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XSWNET_G(JJ)
          DGEI%XSWNET_N(JJ) = DGEI%XSWNET_N(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XSWNET_N(JJ)
          DGEI%XSWNET_NS(JJ) = DGEI%XSWNET_NS(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XSWNET_NS(JJ)
          DGEI%XLWNET_V(JJ) = DGEI%XLWNET_V(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XLWNET_V(JJ)
          DGEI%XLWNET_G(JJ) = DGEI%XLWNET_G(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XLWNET_G(JJ)
          DGEI%XLWNET_N(JJ) = DGEI%XLWNET_N(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XLWNET_N(JJ)
          DGEI%XSWDOWN_GN(JJ) = DGEI%XSWDOWN_GN(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XSWDOWN_GN(JJ)
          DGEI%XLWDOWN_GN(JJ) = DGEI%XLWDOWN_GN(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XLWDOWN_GN(JJ)
          DGEI%XH_V_C(JJ) = DGEI%XH_V_C(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XH_V_C(JJ)
          DGEI%XH_G_C(JJ) = DGEI%XH_G_C(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XH_G_C(JJ)
          DGEI%XH_C_A(JJ) = DGEI%XH_C_A(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XH_C_A(JJ)
          DGEI%XH_N_C(JJ) = DGEI%XH_N_C(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XH_N_C(JJ)
          DGEI%XSR_GN(JJ) = DGEI%XSR_GN(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XSR_GN(JJ)
          DGEI%XMELTCV(JJ) = DGEI%XMELTCV(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XMELTCV(JJ)
          DGEI%XFRZCV(JJ) = DGEI%XFRZCV(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XFRZCV(JJ)
        ENDIF
        !
      ENDIF
    END DO
  ENDDO
!
! Isba water budget and reservoir time tendencies
!
  IF(DGEI%LWATER_BUDGET)THEN
!  
    DGEI%XRAINFALL  (:) = PRAIN(:) * PTSTEP
    DGEI%XSNOWFALL  (:) = PSNOW(:) * PTSTEP
    DGEI%XDWG   (:) = 0.0
    DGEI%XDWGI  (:) = 0.0
    DGEI%XDWR   (:) = 0.0
    DGEI%XDSWE  (:) = 0.0
    DGEI%XWATBUD(:) = 0.0
!
    DO JPATCH=1,SIZE(PPATCH,2)
!     cdir nodep
      DO JJ=1,SIZE(ZSUMPATCH)
        IF (ZSUMPATCH(JJ) > 0.) THEN
!
           DGEI%XDWG   (JJ) = DGEI%XDWG   (JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XDWG   (JJ)
           DGEI%XDWGI  (JJ) = DGEI%XDWGI  (JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XDWGI  (JJ)
           DGEI%XDWR   (JJ) = DGEI%XDWR   (JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XDWR   (JJ)
           DGEI%XDSWE  (JJ) = DGEI%XDSWE  (JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XDSWE  (JJ)
           DGEI%XWATBUD(JJ) = DGEI%XWATBUD(JJ) + PPATCH(JJ,JPATCH) * DGEIP%AL(JPATCH)%XWATBUD(JJ)
!
        ENDIF
      ENDDO
    ENDDO
!
  ENDIF
!
END IF
!
!
!       2.     Surface Cumulated Energy fluxes
!              -------------------------------
!
IF (OSURF_BUDGETC) THEN
   DGEIC%XLEG       (:) = 0.
   DGEIC%XLEGI      (:) = 0.
   DGEIC%XLEV       (:) = 0.
   DGEIC%XLES       (:) = 0.
   DGEIC%XLESL      (:) = 0.
   DGEIC%XLER       (:) = 0.
   DGEIC%XLETR      (:) = 0.
   DGEIC%XSNDRIFT   (:) = 0.
   DGEIC%XDRAIN     (:) = 0.
   DGEIC%XQSB       (:) = 0.
   DGEIC%XRUNOFF    (:) = 0.
   DGEIC%XHORT      (:) = 0.
   DGEIC%XDRIP      (:) = 0.
   DGEIC%XRRVEG     (:) = 0.
   DGEIC%XMELT      (:) = 0.
   DGEIC%XIFLOOD    (:) = 0.
   DGEIC%XPFLOOD    (:) = 0.
   DGEIC%XLE_FLOOD  (:) = 0.
   DGEIC%XLEI_FLOOD (:) = 0.
   DGEIC%XIRRIG_FLUX(:) = 0.
   DGEIC%XGPP       (:) = 0.
   DGEIC%XRESP_AUTO (:) = 0.
   DGEIC%XRESP_ECO  (:) = 0.
!
   IF (ISIZE_LMEB_PATCH>0) THEN
        DGEIC%XLEVCV    (:) = 0.
        DGEIC%XLESC     (:) = 0.
        DGEIC%XLETRGV   (:) = 0.
        DGEIC%XLETRCV   (:) = 0.
        DGEIC%XLERGV    (:) = 0.
        DGEIC%XLERCV    (:) = 0.
        DGEIC%XLE_C_A   (:) = 0.
        DGEIC%XLE_V_C   (:) = 0.
        DGEIC%XLE_G_C   (:) = 0.
        DGEIC%XLE_N_C   (:) = 0.
        DGEIC%XSWNET_V     (:) = 0.
        DGEIC%XSWNET_G     (:) = 0.
        DGEIC%XSWNET_N     (:) = 0.
        DGEIC%XSWNET_NS    (:) = 0.
        DGEIC%XLWNET_V     (:) = 0.
        DGEIC%XLWNET_G     (:) = 0.
        DGEIC%XLWNET_N     (:) = 0.
        DGEIC%XSWDOWN_GN   (:) = 0.
        DGEIC%XLWDOWN_GN   (:) = 0.
        DGEIC%XH_V_C       (:) = 0.
        DGEIC%XH_G_C       (:) = 0.
        DGEIC%XH_C_A       (:) = 0.
        DGEIC%XH_N_C       (:) = 0.
        DGEIC%XSR_GN       (:) = 0.
        DGEIC%XMELTCV      (:) = 0.
        DGEIC%XFRZCV       (:) = 0.
   ENDIF
!
  DO JPATCH=1,SIZE(PPATCH,2)
!cdir nodep
    DO JJ=1,SIZE(ZSUMPATCH)
      IF (ZSUMPATCH(JJ) > 0.) THEN
!
!
! Latent heat of evaporation over the ground
!
        DGEIC%XLEG(JJ)  = DGEIC%XLEG(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XLEG(JJ)
!
! Surface soil ice sublimation
!
        DGEIC%XLEGI(JJ)  = DGEIC%XLEGI(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XLEGI(JJ)
!
! Latent heat of evaporation over vegetation
!
        DGEIC%XLEV(JJ)  = DGEIC%XLEV(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XLEV(JJ)
!
! Latent heat of sublimation over snow
!
        DGEIC%XLES(JJ)  = DGEIC%XLES(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XLES(JJ)
!
! Latent heat of evaporation of liquid water over snow
!
        DGEIC%XLESL(JJ)  = DGEIC%XLESL(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XLESL(JJ)
!
! Evaporation from canopy water interception
!
        DGEIC%XLER(JJ)  = DGEIC%XLER(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XLER(JJ)
!
! Evapotranspiration of the vegetation
!
        DGEIC%XLETR(JJ)  = DGEIC%XLETR(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XLETR(JJ)
!
! Blowing snow sublimation (ES or Crocus)
!
        DGEIC%XSNDRIFT(JJ)  = DGEIC%XSNDRIFT(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XSNDRIFT(JJ)
!
! Soil drainage flux
!
        DGEIC%XDRAIN(JJ)  = DGEIC%XDRAIN(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XDRAIN(JJ)
!
! Soil lateral subsurface flux
!
        DGEIC%XQSB(JJ)  = DGEIC%XQSB(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XQSB(JJ)        
!
! Supersaturation runoff
!
        DGEIC%XRUNOFF(JJ)  = DGEIC%XRUNOFF(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XRUNOFF(JJ)
!
! Horton runoff
!
        DGEIC%XHORT(JJ)  = DGEIC%XHORT(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XHORT(JJ)
!
! Vegetation dripping
!
        DGEIC%XDRIP(JJ)  = DGEIC%XDRIP(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XDRIP(JJ)
!
! precipitation intercepted by the vegetation
!
        DGEIC%XRRVEG(JJ)  = DGEIC%XRRVEG(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XRRVEG(JJ)
!      
! Snow melt
!
        DGEIC%XMELT(JJ)  = DGEIC%XMELT(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XMELT(JJ)
!      
! Flood infiltartion
!
        DGEIC%XIFLOOD(JJ) = DGEIC%XIFLOOD(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XIFLOOD(JJ)
!      
! Precipitation intercepted by the floodplains
!     
        DGEIC%XPFLOOD(JJ) = DGEIC%XPFLOOD(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XPFLOOD(JJ)
!      
! Floodplains evaporation
!     
        DGEIC%XLE_FLOOD(JJ) = DGEIC%XLE_FLOOD(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XLE_FLOOD(JJ)
        DGEIC%XLEI_FLOOD(JJ) = DGEIC%XLEI_FLOOD(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XLEI_FLOOD(JJ)
!      
! irrigation rate (as soil input)
!
        DGEIC%XIRRIG_FLUX(JJ)  = DGEIC%XIRRIG_FLUX(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XIRRIG_FLUX(JJ)
!
! Gross primary production
!
        DGEIC%XGPP(JJ) = DGEIC%XGPP(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XGPP(JJ)
!
! Autotrophic respiration
!   
        DGEIC%XRESP_AUTO(JJ) = DGEIC%XRESP_AUTO(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XRESP_AUTO(JJ)
!
! Ecosystem respiration
!
        DGEIC%XRESP_ECO(JJ) = DGEIC%XRESP_ECO(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XRESP_ECO(JJ)
!      
        IF (ISIZE_LMEB_PATCH>0) THEN
          DGEIC%XLEVCV(JJ) = DGEIC%XLEVCV(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XLEVCV(JJ)
          DGEIC%XLESC(JJ) = DGEIC%XLESC(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XLESC(JJ)
!         DGEIC%XLETRGV(JJ) = DGEIC%XLETRGV(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XLETRGV(JJ)
          DGEIC%XLETRCV(JJ) = DGEIC%XLETRCV(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XLETRCV(JJ)
!         DGEIC%XLERGV(JJ) = DGEIC%XLERGV(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XLERGV(JJ)
          DGEIC%XLERCV(JJ) = DGEIC%XLERCV(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XLERCV(JJ)
          DGEIC%XLE_C_A(JJ) = DGEIC%XLE_C_A(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XLE_C_A(JJ)
          DGEIC%XLE_V_C(JJ) = DGEIC%XLE_V_C(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XLE_V_C(JJ)
          DGEIC%XLE_G_C(JJ) = DGEIC%XLE_G_C(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XLE_G_C(JJ)
          DGEIC%XLE_N_C(JJ) = DGEIC%XLE_N_C(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XLE_N_C(JJ)
          DGEIC%XSWNET_V(JJ) = DGEIC%XSWNET_V(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XSWNET_V(JJ)
          DGEIC%XSWNET_G(JJ) = DGEIC%XSWNET_G(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XSWNET_G(JJ)
          DGEIC%XSWNET_N(JJ) = DGEIC%XSWNET_N(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XSWNET_N(JJ)
          DGEIC%XSWNET_NS(JJ) = DGEIC%XSWNET_NS(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XSWNET_NS(JJ)
          DGEIC%XLWNET_V(JJ) = DGEIC%XLWNET_V(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XLWNET_V(JJ)
          DGEIC%XLWNET_G(JJ) = DGEIC%XLWNET_G(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XLWNET_G(JJ)
          DGEIC%XLWNET_N(JJ) = DGEIC%XLWNET_N(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XLWNET_N(JJ)
          DGEIC%XSWDOWN_GN(JJ) = DGEIC%XSWDOWN_GN(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XSWDOWN_GN(JJ)
          DGEIC%XLWDOWN_GN(JJ) = DGEIC%XLWDOWN_GN(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XLWDOWN_GN(JJ)
          DGEIC%XH_V_C(JJ) = DGEIC%XH_V_C(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XH_V_C(JJ)
          DGEIC%XH_G_C(JJ) = DGEIC%XH_G_C(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XH_G_C(JJ)
          DGEIC%XH_C_A(JJ) = DGEIC%XH_C_A(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XH_C_A(JJ)
          DGEIC%XH_N_C(JJ) = DGEIC%XH_N_C(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XH_N_C(JJ)
          DGEIC%XSR_GN(JJ) = DGEIC%XSR_GN(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XSR_GN(JJ)
          DGEIC%XMELTCV(JJ) = DGEIC%XMELTCV(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XMELTCV(JJ)
          DGEIC%XFRZCV(JJ) = DGEIC%XFRZCV(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XFRZCV(JJ)
        ENDIF
        !
      ENDIF
    ENDDO
  END DO
!
! Isba water budget and reservoir time tendencies
!
  IF(DGEI%LWATER_BUDGET)THEN
!  
    DGEIC%XRAINFALL  (:) = DGEIC%XRAINFALL (:) + PRAIN(:) * PTSTEP
    DGEIC%XSNOWFALL  (:) = DGEIC%XSNOWFALL (:) + PSNOW(:) * PTSTEP
    DGEIC%XDWG   (:) = 0.0
    DGEIC%XDWGI  (:) = 0.0
    DGEIC%XDWR   (:) = 0.0
    DGEIC%XDSWE  (:) = 0.0
    DGEIC%XWATBUD(:) = 0.0
!
    DO JPATCH=1,SIZE(PPATCH,2)
!     cdir nodep
      DO JJ=1,SIZE(ZSUMPATCH)
        IF (ZSUMPATCH(JJ) > 0.) THEN
!
           DGEIC%XDWG   (JJ) = DGEIC%XDWG   (JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XDWG   (JJ)
           DGEIC%XDWGI  (JJ) = DGEIC%XDWGI  (JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XDWGI  (JJ)
           DGEIC%XDWR   (JJ) = DGEIC%XDWR   (JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XDWR   (JJ)
           DGEIC%XDSWE  (JJ) = DGEIC%XDSWE  (JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XDSWE  (JJ)
           DGEIC%XWATBUD(JJ) = DGEIC%XWATBUD(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XWATBUD(JJ)
!
        ENDIF
      ENDDO
    ENDDO
!
  ENDIF
!
! Ice calving flux
!  
  IF(OGLACIER)THEN 
    DGEIC%XICEFLUX(:)= 0.
    DO JPATCH=1,SIZE(PPATCH,2)
!     cdir nodep  
      DO JJ=1,SIZE(ZSUMPATCH)
         IF(ZSUMPATCH(JJ) > 0.)THEN
            DGEIC%XICEFLUX(JJ) = DGEIC%XICEFLUX(JJ) + PPATCH(JJ,JPATCH) * DGEIPC%AL(JPATCH)%XICEFLUX(JJ)      
         ENDIF
      END DO
    END DO
  END IF
!  
END IF
!
IF (LHOOK) CALL DR_HOOK('AVERAGE_DIAG_EVAP_ISBA_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE AVERAGE_DIAG_EVAP_ISBA_n
