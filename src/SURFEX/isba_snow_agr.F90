!     #########
      SUBROUTINE ISBA_SNOW_AGR(IP, INI, IR, DGI, DGEI, DGMI, DGIP, DGEIP, &
                               HSNOW_ISBA, OMEB, PEXNS, PEXNA, PTA, PQA,  &
                               PZREF, PUREF, PDIRCOSZW, PVMOD, PRR, PSR,  &
                                PEMIS, PALB, PUSTAR, PLES3L, PLEL3L,      &
                                PEVAP3L, PQS3L, PALB3L, PGSFCSNOW,        &
                                PZGRNDFLUX, PFLSN_COR, PEMIST, PPALPHAN )
!     ##########################################################################
!
!
!!****  *ISBA_SNOW_AGR* aggregates snow free and snow fluxes
!!
!!    PURPOSE
!!    -------
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
!!    Noilhan and Planton (1989)
!!      
!!    AUTHOR
!!    ------
!!      V. Masson           * Meteo-France *
!!      (following A. Boone)
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/03/95 
!!      B. Decharme 01/2009  Floodplains 
!!      B. Decharme 01/2010  Effective surface temperature (for diag)
!!      B. Decharme 09/2012  Bug total sublimation flux: no DGEIP%XLESL
!!      B. Decharme 04/2013  Bug wrong radiative temperature
!!                           Sublimation diag flux
!!                           Qs for 3l or crocus (needed for coupling with atm)
!!      A. Boone    11/2014  MEB
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t, ISBA_INIT_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
USE MODD_DIAG_n, ONLY : DIAG_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!              -------------------------
!
!
!* general variables
!  -----------------
!
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
TYPE(ISBA_INIT_t), INTENT(INOUT) :: INI
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
TYPE(DIAG_t), INTENT(INOUT) :: DGI
TYPE(DIAG_t), INTENT(INOUT) :: DGIP
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEI
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEIP
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DGMI
!
 CHARACTER(LEN=*),     INTENT(IN)  :: HSNOW_ISBA ! 'DEF' = Default F-R snow scheme
!                                               !         (Douville et al. 1995)
!                                               ! '3-L' = 3-L snow scheme (option)
!                                               !         (Boone and Etchevers 2000)
!
LOGICAL,              INTENT(IN)  :: OMEB       ! True = patch with multi-energy balance 
!                                               ! False = patch with classical ISBA 
!
!* surface and atmospheric parameters
!  ----------------------------------
!
REAL, DIMENSION(:), INTENT(IN)  :: PEXNS     ! Exner function at the surface
REAL, DIMENSION(:), INTENT(IN)  :: PEXNA     ! Exner function
REAL, DIMENSION(:), INTENT(IN)  :: PTA       ! air temperature
REAL, DIMENSION(:), INTENT(IN)  :: PQA       ! air specific humidity
REAL, DIMENSION(:), INTENT(IN)  :: PZREF     ! reference height of the first atmospheric level
REAL, DIMENSION(:), INTENT(IN)  :: PUREF     ! reference height of the wind
REAL, DIMENSION(:), INTENT(IN)  :: PDIRCOSZW ! Cosinus of the angle between the normal to the surface and the vertical
REAL, DIMENSION(:), INTENT(IN)  :: PVMOD     ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)  :: PRR       ! Rain rate (in kg/m2/s)
REAL, DIMENSION(:), INTENT(IN)  :: PSR       ! Snow rate (in kg/m2/s)
!
!* surface parameters
!  ------------------
!
REAL, DIMENSION(:), INTENT(IN)  :: PALB       ! albedo 
REAL, DIMENSION(:), INTENT(IN)  :: PEMIS      ! emissivity
!  'D95'     : represents aggregated (snow + flood + snow-flood-free) albedo and emissivity
!  '3-L'     : represents                    flood + snow-flood-free  albedo and emissivity
!  'MEB+3-L' : represents aggregated (snow + flood + snow-flood-free) albedo and emissivity
!
!
!* snow fractions
!  --------------
!
REAL, DIMENSION(:), INTENT(IN)  :: PPALPHAN   ! fraction of the the explicit veg.
!                                             ! canopy buried by snow
!
!* ISBA-SNOW3L variables/parameters:
!  ---------------------------------
!
! Prognostic variables:
!
REAL, DIMENSION(:),   INTENT(IN) :: PALB3L      ! Snow albedo
REAL, DIMENSION(:),   INTENT(IN) :: PQS3L       ! Surface humidity
! 
! Diagnostics:
!
REAL, DIMENSION(:), INTENT(IN)    :: PZGRNDFLUX ! snow/soil-biomass interface flux (W/m2)
REAL, DIMENSION(:), INTENT(IN)    :: PFLSN_COR  ! snow/soil-biomass correction flux (W/m2)
!
REAL, DIMENSION(:), INTENT(INOUT) :: PGSFCSNOW  ! heat flux from snow sfc to sub sfc layers (W/m2)
REAL, DIMENSION(:), INTENT(IN)    :: PLES3L     ! sublimation from ISBA-ES(3L)
REAL, DIMENSION(:), INTENT(IN)    :: PLEL3L     ! evaporation heat flux of water in the snow (W/m2)
REAL, DIMENSION(:), INTENT(INOUT) :: PEVAP3L    ! evaporation flux over snow from ISBA-ES (kg/m2/s)
!
!* diagnostic variables
!  --------------------
!
REAL, DIMENSION(:), INTENT(INOUT) :: PEMIST   ! total surface emissivity
!
!* surface fluxes
!  --------------
!
REAL, DIMENSION(:), INTENT(INOUT) :: PUSTAR   ! friction velocity
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PTA))       :: ZWORK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ISBA_SNOW_AGR',0,ZHOOK_HANDLE)
!
ZWORK(:) = 0.
!
IF(OMEB)THEN
!
! MEB uses only an explicit scheme option.
! Fluxes enter here as "snow relative", here we
! transform back to "patch or grid box relative" (by incorporating
! the notion of fractional coverage)
!
   DGEIP%XLES     (:) = IR%XPSN(:,1) * PLES3L(:)
   DGEIP%XLESL    (:) = IR%XPSN(:,1) * PLEL3L(:)
   DGEIP%XLWNET_N (:) = IR%XPSN(:,1) * DGEIP%XLWNET_N (:)
   DGEIP%XSWNET_N (:) = IR%XPSN(:,1) * DGEIP%XSWNET_N (:)
   DGEIP%XSWNET_NS(:) = IR%XPSN(:,1) * DGEIP%XSWNET_NS(:)  
   DGEIP%XMELTADV (:) = IR%XPSN(:,1) * DGEIP%XMELTADV (:) 
   DGMI%XRNSNOW   (:) = IR%XPSN(:,1) * DGMI%XRNSNOW   (:)
   DGMI%XHSNOW    (:) = IR%XPSN(:,1) * DGMI%XHSNOW    (:)
   DGMI%XGFLUXSNOW(:) = IR%XPSN(:,1) * DGMI%XGFLUXSNOW(:)
   DGMI%XSNOWHMASS(:) = IR%XPSN(:,1) * DGMI%XSNOWHMASS(:)  
   DGMI%XHPSNOW   (:) = IR%XPSN(:,1) * DGMI%XHPSNOW   (:)  
   DGMI%XGRNDFLUX (:) = IR%XPSN(:,1) * (PZGRNDFLUX(:)+PFLSN_COR(:)) 
   PEVAP3L  (:) = IR%XPSN(:,1) * PEVAP3L(:)   
   PGSFCSNOW(:) = IR%XPSN(:,1) * PGSFCSNOW(:)

! Snow free (ground-based snow) diagnostics: canopy and ground blended (W m-2):
! NOTE that the effects of snow cover *fraction* are implicitly *included* in these fluxes 
! so do NOT multiply by snow fraction.
!
   DGEI%XLEG  (:) = DGEIP%XLEG   (:)
   DGEI%XLEGI (:) = DGEIP%XLEGI  (:)
   DGEI%XLEV  (:) = DGEIP%XLEVCV (:)
   DGEI%XLETR (:) = DGEIP%XLETRCV(:) 
! NOTE for now, this is same as total Ustar (includes snow)   
   DGEI%XUSTAR(:) = PUSTAR       (:)        
! LER does not include intercepted snow sublimation
   DGEI%XLER  (:) = DGEIP%XLEVCV(:) - DGEIP%XLETRCV(:) 

   DGI%XRN   (:) = DGEIP%XSWNET_V(:) + DGEIP%XSWNET_G(:) + DGEIP%XLWNET_V(:) + DGEIP%XLWNET_G(:)
   DGI%XH    (:) = DGEIP%XH_V_C(:) + DGEIP%XH_G_C(:)
   DGI%XLEI  (:) = DGEIP%XLEGI(:) + DGEIP%XLEI_FLOOD(:) + DGEIP%XLES(:) + DGEIP%XLESC(:)
  ! LE includes intercepted snow sublimation
   DGI%XLE   (:) = DGEI%XLEG(:) + DGEI%XLEGI(:) + DGEI%XLEV(:) + DGEIP%XLESC(:) + &
                   DGEIP%XLE_FLOOD(:) + DGEIP%XLEI_FLOOD(:)
   DGI%XGFLUX(:) = DGI%XRN(:) - DGI%XH(:) - DGI%XLE(:)
   DGI%XLEI  (:) = DGEIP%XLEGI(:) + DGEIP%XLEI_FLOOD(:) 
!
   PEMIST(:) = PEMIS(:)
!
! Effective surface temperature (for diag): for MEB:

   ZWORK   (:) =  PPALPHAN(:)*IR%XPSN(:,1)
   DGIP%XTS(:) = (1.0 - ZWORK(:))*IR%XTC(:,1) + ZWORK(:)*DGMI%XSNOWTEMP(:,1)
!
! Total heat FLUX into snow/soil/vegetation surface:
!
   DGIP%XGFLUX(:) = DGIP%XRN(:) - DGIP%XH(:) - IR%XLE(:,1) + DGMI%XHPSNOW(:) 
!
ELSE
!
! * 2. Using an explicit snow scheme option with composite soil/veg ISBA:
!      ------------------------------------------------------------------
!
   IF(HSNOW_ISBA == '3-L' .OR. HSNOW_ISBA == 'CRO')THEN
!
!     Save fluxes from Force-Restore snow/explicit snow-free
!     portion of grid box (vegetation/soil):
!
      DGEI%XLEG  (:) = DGEIP%XLEG (:)
      DGEI%XLEGI (:) = DGEIP%XLEGI(:)
      DGEI%XLEV  (:) = DGEIP%XLEV (:)
      DGEI%XLETR (:) = DGEIP%XLETR(:)
      DGEI%XLER  (:) = DGEIP%XLER (:)
      DGI%XRN    (:) = DGIP%XRN   (:)
      DGI%XH     (:) = DGIP%XH    (:)
      DGEI%XUSTAR(:) = PUSTAR     (:)      

      DGI%XLE   (:) = IR%XLE(:,1)
      DGI%XGFLUX(:) = DGIP%XGFLUX(:)
!  
      DGI%XLEI  (:)= DGEIP%XLEGI(:) + DGEIP%XLEI_FLOOD(:) + DGEIP%XLES(:)
!
!     Effective surface temperature (for diag):
!
      DGIP%XTS(:) = (1.-IR%XPSN(:,1))*IR%XTG(:,1,1)+IR%XPSN(:,1)*DGMI%XSNOWTEMP(:,1)
!
!     Effective surface radiating temperature:
!
      DGIP%XALBT (:) = PALB (:)*(1.-IR%XPSN(:,1)) + IR%XPSN(:,1)*PALB3L(:)
      PEMIST     (:) = PEMIS(:)*(1.-IR%XPSN(:,1)) + IR%XPSN(:,1)*IR%TSNOW%EMIS(:,1)
!  
      DGIP%XTSRAD(:) = ( ((1.-IR%XPSN(:,1))*PEMIS(:)*IR%XTG(:,1,1)**4 +   &
                               IR%XPSN(:,1) *IR%TSNOW%EMIS(:,1)*DGMI%XSNOWTEMP(:,1)**4     &
                            )/PEMIST(:) )**(0.25)  
!
!     Calculate actual fluxes from snow-free natural
!     portion of surface: NET flux from surface is the sum of
!     fluxes from snow free and snow covered portions
!     of natural portion of grid box when *ISBA-ES* in force.
!     when NOT in use, then these fluxes equal those above.
!
      DGIP%XRN   (:) = (1.-IR%XPSN(:,1)) * DGIP%XRN(:) + IR%XPSN(:,1) * DGMI%XRNSNOW(:)
      DGIP%XH    (:) = (1.-IR%XPSN(:,1)) * DGIP%XH (:) + IR%XPSN(:,1) * DGMI%XHSNOW(:)
!  
      DGEIP%XLEG (:) = (1.-IR%XPSNG(:,1)-INI%XFFG(:,1)) * DGEIP%XLEG(:)  
      DGEIP%XLEGI(:) = (1.-IR%XPSNG(:,1)-INI%XFFG(:,1)) * DGEIP%XLEGI(:)  
      DGEIP%XLEV (:) = (1.-IR%XPSNV(:,1)-INI%XFFV(:,1)) * DGEIP%XLEV(:)   
      DGEIP%XLETR(:) = (1.-IR%XPSNV(:,1)-INI%XFFV(:,1)) * DGEIP%XLETR(:)  
      DGEIP%XLER (:) = (1.-IR%XPSNV(:,1)-INI%XFFV(:,1)) * DGEIP%XLER(:)  
!
!     Total evapotranspiration flux (kg/m2/s):
!
      DGIP%XEVAP(:) = (DGEIP%XLEV(:) + DGEIP%XLEG(:))/IP%XPLVTT(:,1) + DGEIP%XLEGI(:)/IP%XPLSTT(:,1) + &
                                   DGEIP%XLE_FLOOD(:)/IP%XPLVTT(:,1) + DGEIP%XLEI_FLOOD(:)/IP%XPLSTT(:,1) + &
                                   IR%XPSN(:,1) * PEVAP3L(:)
!
!     ISBA-ES/SNOW3L fluxes:
!
      DGEIP%XLES     (:) = IR%XPSN(:,1) * PLES3L(:)
      DGEIP%XLESL    (:) = IR%XPSN(:,1) * PLEL3L(:)
      DGEIP%XSWNET_N (:) = IR%XPSN(:,1) * DGEIP%XSWNET_N (:)
      DGEIP%XSWNET_NS(:) = IR%XPSN(:,1) * DGEIP%XSWNET_NS(:)     
      DGMI%XRNSNOW   (:) = IR%XPSN(:,1) * DGMI%XRNSNOW   (:)
      DGMI%XHSNOW    (:) = IR%XPSN(:,1) * DGMI%XHSNOW    (:)
      DGMI%XGFLUXSNOW(:) = IR%XPSN(:,1) * DGMI%XGFLUXSNOW(:)
      DGMI%XSNOWHMASS(:) = IR%XPSN(:,1) * DGMI%XSNOWHMASS(:)  ! (J m-2)
      DGMI%XHPSNOW   (:) = IR%XPSN(:,1) * DGMI%XHPSNOW   (:)
      PGSFCSNOW      (:) = IR%XPSN(:,1) * PGSFCSNOW(:)
      PEVAP3L        (:) = IR%XPSN(:,1) * PEVAP3L  (:)
!
!     Total heat flux between snow and soil
!
      DGMI%XGRNDFLUX(:) = IR%XPSN(:,1) * (PZGRNDFLUX(:)+PFLSN_COR(:))
      DGEIP%XMELTADV(:) = IR%XPSN(:,1) * DGEIP%XMELTADV(:)
!
!     Total evaporative flux (W/m2) :
!
      IR%XLE(:,1) = DGEIP%XLEG(:) + DGEIP%XLEV(:) + DGEIP%XLES(:) + DGEIP%XLESL(:) + &
                   DGEIP%XLEGI(:) + DGEIP%XLE_FLOOD(:) + DGEIP%XLEI_FLOOD(:)
!
!     Total sublimation flux (W/m2) :
!
      DGIP%XLEI  (:) = DGEIP%XLES(:) + DGEIP%XLEGI(:) + DGEIP%XLEI_FLOOD(:)
!
!     Total sublimation flux (kg/m2/s):
!
      DGIP%XSUBL (:) = DGIP%XLEI(:)/IP%XPLSTT(:,1)
!
!     Total FLUX into snow/soil/vegetation surface:
!
      DGIP%XGFLUX(:) = DGIP%XRN(:) - DGIP%XH(:) - IR%XLE(:,1) + DGMI%XHPSNOW(:)  
!
!     surface humidity:
!
      DGIP%XQS   (:) = (1.-IR%XPSN(:,1)) * DGIP%XQS(:) + IR%XPSN(:,1) * PQS3L(:)
!
!     near-surface humidity :
!  
      DGIP%XHU   (:) = (1.-IR%XPSN(:,1)) * DGIP%XHU(:) + IR%XPSN(:,1)
!
!     Momentum fluxes:
!
      PUSTAR     (:) = SQRT( (1.-IR%XPSN(:,1))  * PUSTAR(:)**2  + &
                                 IR%XPSN(:,1) * DGMI%XUSTARSNOW(:)**2 )
!
!     Richardson number and Drag coeff:
!
      CALL COMPUT_RI_DRAG
!
   ELSE
!
      DGIP%XTS    (:)  = IR%XTG(:,1,1)
      DGIP%XTSRAD (:)  = IR%XTG(:,1,1)
      DGIP%XALBT  (:)  = PALB (:)
      PEMIST      (:)  = PEMIS(:)
!  
!     Total sublimation flux (W/m2) :
      DGIP%XLEI   (:)  = DGEIP%XLES(:) + DGEIP%XLEGI(:) + DGEIP%XLEI_FLOOD(:)
!
!     Total sublimation flux (kg/m2/s):
      DGIP%XSUBL  (:)  = DGIP%XLEI(:)/IP%XPLSTT(:,1)
!
   ENDIF
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('ISBA_SNOW_AGR',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!
SUBROUTINE COMPUT_RI_DRAG
!
USE MODD_SURF_ATM, ONLY : LDRAG_COEF_ARP, LRRGUST_ARP,   &
                          XRRSCALE, XRRGAMMA, XUTILGUST  
!
USE MODI_SURFACE_RI
USE MODI_SURFACE_AERO_COND
USE MODI_SURFACE_CD
USE MODI_SURFACE_CDCH_1DARP
USE MODI_WIND_THRESHOLD
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PTA)) :: ZFP, ZRRCOR, ZVMOD, ZAC, ZRA
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ISBA_SNOW_AGR:COMPUT_RI_DRAG',0,ZHOOK_HANDLE)
!
! * Richardson number
!
CALL SURFACE_RI(DGIP%XTS, DGIP%XQS, PEXNS, PEXNA, PTA, PQA,  &
                PZREF, PUREF, PDIRCOSZW, PVMOD, DGIP%XRI)  
!
! * Wind check
!
ZVMOD = WIND_THRESHOLD(PVMOD,PUREF)
!
! * Drag coefficient for heat and momentum
!
IF (LDRAG_COEF_ARP) THEN
   CALL SURFACE_CDCH_1DARP(PZREF, DGIP%XZ0EFF, DGIP%XZ0H, ZVMOD, PTA, IR%XTG(:,1,1), &
                             PQA, DGIP%XQS, DGIP%XCD, DGIP%XCDN, DGIP%XCH              )
ELSE
   CALL SURFACE_AERO_COND(DGIP%XRI, PZREF, PUREF, ZVMOD, DGIP%XZ0, DGIP%XZ0H, ZAC, ZRA, DGIP%XCH)
   CALL SURFACE_CD(DGIP%XRI, PZREF, PUREF, DGIP%XZ0EFF, DGIP%XZ0H, DGIP%XCD, DGIP%XCDN)
ENDIF
!
IF (LRRGUST_ARP) THEN
   ZFP(:)=MAX(0.0,PRR(:)+PSR(:))
   ZRRCOR(:)=SQRT(1.0+((((ZFP(:)/(ZFP(:)+XRRSCALE))**XRRGAMMA)*XUTILGUST)**2) &
       /(DGIP%XCD(:)*ZVMOD(:)**2))  
   DGIP%XCD  = DGIP%XCD  * ZRRCOR
   DGIP%XCH  = DGIP%XCH  * ZRRCOR
   DGIP%XCDN = DGIP%XCDN * ZRRCOR
ENDIF
!
IF (LHOOK) CALL DR_HOOK('ISBA_SNOW_AGR:COMPUT_RI_DRAG',1,ZHOOK_HANDLE)
!
END SUBROUTINE COMPUT_RI_DRAG
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ISBA_SNOW_AGR
