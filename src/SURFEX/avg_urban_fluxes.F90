!     #########
    SUBROUTINE AVG_URBAN_FLUXES(                                              &
                       PTS_TOWN, PEMIS_TOWN,                                    &
                       PT_CANYON, PQ_CANYON,                                    &
                       PT_LOWCAN, PQ_LOWCAN,                                    &
                       PTS_ROOF,PTS_ROAD,PTS_WALL,PTS_GARDEN,                   &
                       PTA, PQA, PRHOA, PPS,                                    &
                       PH_TRAFFIC, PLE_TRAFFIC, PH_INDUSTRY, PLE_INDUSTRY,      &
                       PBLD, PROAD, PGARDEN, PWALL_O_HOR, PWALL_O_GRND,         &
                       PEMIS_ROOF,                                              &
                       PESNOW_ROOF,                                             &
                       PLW_RAD, PLW_S_TO_W, PLW_S_TO_R, PLW_S_TO_G,             &
                       PLW_S_TO_NR,                                              &
                       PABS_LW_ROOF, PABS_LW_WALL, PABS_LW_ROAD, PABS_LW_GARDEN,&
                       PABS_LW_SNOW_ROOF, PABS_LW_SNOW_ROAD,                    &
                       PAC_ROOF, PAC_ROOF_WAT,                                  &
                       PAC_WALL, PAC_ROAD, PAC_ROAD_WAT, PAC_TOP,               &
                       PAC_GARDEN,                                              &
                       PQSAT_GARDEN, PAC_AGG_GARDEN, PHU_AGG_GARDEN,            &
                       PQSAT_ROOF, PQSAT_ROAD,                                  &
                       PDELT_ROOF, PDELT_ROAD,                                  &
                       PROOF_FRAC, PWALL_FRAC, PROAD_FRAC, PGARDEN_FRAC,        &
                       PTOTS_O_HORS,                                            &
                       PDF_ROOF, PDN_ROOF, PDF_ROAD, PDN_ROAD,                  &
                       PRN_ROOF, PH_ROOF, PLE_ROOF, PGFLUX_ROOF,                &
                       PRN_ROAD, PH_ROAD, PLE_ROAD, PGFLUX_ROAD,                &
                       PRN_GARDEN, PH_GARDEN, PLE_GARDEN, PGFLUX_GARDEN,        &
                       PRN_WALL, PH_WALL, PLE_WALL, PGFLUX_WALL,                &
                       PLEW_ROOF, PLESNOW_ROOF,                                 &
                       PLEW_ROAD, PLESNOW_ROAD, PHSNOW_ROAD,                    &
                       PEVAP_GARDEN,                                            &
                       PRN_GRND, PH_GRND, PLE_GRND, PGFLUX_GRND,                &
                       PRN_TOWN, PH_TOWN, PLE_TOWN, PGFLUX_TOWN, PEVAP_TOWN     )  
!   ##########################################################################
!
!!****  *AVG_URBAN_FLUXES* computes fluxes on urbanized surfaces  
!!
!!    PURPOSE
!!    -------
!         
!     
!!**  METHOD
!     ------
!
!
!
!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    MODD_CST
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!	V. Masson           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    23/01/98 
!!                     12/02 (A. Lemonsu) modifications of emissivity and Tstown
!!                     07/07 (P. LeMoigne) expression of latent heat fluxes as 
!!                           functions of w'theta' instead of w'T'
!!                     17/10 (G. Pigeon)  computation of anthropogenic heat due
!!                            to domestic heating
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,ONLY : XCPD, XLVTT, XLSTT, XSTEFAN
!
USE MODE_THERMOS
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
REAL, DIMENSION(:), INTENT(OUT)   :: PTS_TOWN          ! town surface temperature
REAL, DIMENSION(:), INTENT(OUT)   :: PEMIS_TOWN        ! town equivalent emissivity
REAL, DIMENSION(:), INTENT(INOUT) :: PT_CANYON         ! canyon air temperature
REAL, DIMENSION(:), INTENT(INOUT) :: PQ_CANYON         ! canyon air specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PT_LOWCAN         ! low canyon air temperature
REAL, DIMENSION(:), INTENT(IN)    :: PQ_LOWCAN         ! low canyon air specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PTS_ROOF          ! roof surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_ROAD          ! road surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_WALL          ! wall surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_GARDEN        ! green area surface temperature

REAL, DIMENSION(:), INTENT(IN)    :: PTA               ! temperature at roof level
REAL, DIMENSION(:), INTENT(IN)    :: PQA               ! specific humidity
                                                       ! at roof level
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA             ! air density
                                                       ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PPS               ! surface pressure
REAL, DIMENSION(:), INTENT(IN)    :: PH_TRAFFIC        ! anthropogenic sensible
!                                                      ! heat fluxes due to traffic
REAL, DIMENSION(:), INTENT(IN)    :: PLE_TRAFFIC       ! anthropogenic latent
!                                                      ! heat fluxes due to traffic
REAL, DIMENSION(:), INTENT(IN)    :: PH_INDUSTRY       ! anthropogenic sensible
!                                                      ! heat fluxes due to factories
REAL, DIMENSION(:), INTENT(IN)    :: PLE_INDUSTRY      ! anthropogenic latent
!                                                      ! heat fluxes due to factories
REAL, DIMENSION(:), INTENT(IN)    :: PBLD              ! fraction of buildings
REAL, DIMENSION(:), INTENT(IN)    :: PROAD             ! fraction of roads
REAL, DIMENSION(:), INTENT(IN)    :: PGARDEN           ! fraction of green areas
REAL, DIMENSION(:), INTENT(IN)    :: PWALL_O_HOR       ! wall Surf. / (bld+road+green) Surf.
REAL, DIMENSION(:), INTENT(IN)    :: PWALL_O_GRND      ! wall Surf. / ground (road+green) Surf.
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_ROOF        ! roof emissivity
! 
REAL, DIMENSION(:), INTENT(IN)    :: PESNOW_ROOF       ! snow roof emissivity
!
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD           ! incoming longwave rad.
!
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_W        ! LW contrib sky -> wall
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_R        ! LW contrib sky -> road
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_G        ! LW contrib sky -> green
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_NR       ! LW contrib sky -> road (snow)
!
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_ROOF      ! absorbed LW rad. by roof
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_WALL      ! absorbed LW rad. by wall
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_ROAD      ! absorbed LW rad. by road
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_GARDEN    ! absorbed LW rad. by green areas
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_SNOW_ROOF ! abs. LW rad. by snow
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_SNOW_ROAD ! abs. LW rad. by snow
!
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROOF          ! surface conductance
!                                                      ! for heat transfers
!                                                      ! above roofs
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROOF_WAT      ! surface conductance
!                                                      ! for heat transfers
!                                                      ! above roofs (for water)
REAL, DIMENSION(:), INTENT(IN)    :: PAC_WALL          ! surface conductance
!                                                      ! for heat transfer
!                                                      ! between wall and canyon
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROAD          ! surface conductance
!                                                      ! for heat transfers
!                                                      ! between road and canyon
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROAD_WAT      ! surface conductance
!                                                      ! for heat transfers
!                                                      ! inside canyon (for water)
REAL, DIMENSION(:), INTENT(IN)    :: PAC_TOP           ! aerodynamical conductance
!                                                      ! between atmosphere and
!                                                      ! canyon top
REAL, DIMENSION(:), INTENT(IN)    :: PAC_GARDEN        ! aerodynamical conductance
!                                                      ! between atmosphere and
!                                                      ! green areas
REAL, DIMENSION(:), INTENT(IN)    :: PQSAT_GARDEN      ! q_sat(Ts)
REAL, DIMENSION(:), INTENT(IN)    :: PAC_AGG_GARDEN    ! aggregated aerodyn resistance for green areas
REAL, DIMENSION(:), INTENT(IN)    :: PHU_AGG_GARDEN    ! aggregated relative humidity for green areas
!
REAL, DIMENSION(:), INTENT(IN)    :: PQSAT_ROOF        ! q_sat(Ts)
REAL, DIMENSION(:), INTENT(IN)    :: PQSAT_ROAD        ! q_sat(Ts)
REAL, DIMENSION(:), INTENT(IN)    :: PDELT_ROOF        ! water fraction on snow-free
REAL, DIMENSION(:), INTENT(IN)    :: PDELT_ROAD        ! roof and roads
REAL, DIMENSION(:), INTENT(IN)    :: PROOF_FRAC        ! roof, wall,
REAL, DIMENSION(:), INTENT(IN)    :: PWALL_FRAC        ! road, and green area
REAL, DIMENSION(:), INTENT(IN)    :: PROAD_FRAC        ! fractions
REAL, DIMENSION(:), INTENT(IN)    :: PGARDEN_FRAC      ! of exchange surf.
REAL, DIMENSION(:), INTENT(IN)    :: PTOTS_O_HORS      ! total canyon+roof surface
!                                                      ! over horizontal surface
REAL, DIMENSION(:), INTENT(IN)    :: PDF_ROOF          ! snow-free    roof fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROOF          ! snow-covered roof fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDF_ROAD          ! snow-free    road fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROAD          ! snow-covered road fraction
!
REAL, DIMENSION(:), INTENT(IN)    :: PRN_ROOF          ! net radiation over roof
REAL, DIMENSION(:), INTENT(IN)    :: PH_ROOF           ! sensible heat flux over roof
REAL, DIMENSION(:), INTENT(IN)    :: PLE_ROOF          ! latent heat flux over roof
REAL, DIMENSION(:), INTENT(IN)    :: PGFLUX_ROOF       ! flux through the roof
REAL, DIMENSION(:), INTENT(IN)    :: PRN_ROAD          ! net radiation over road
REAL, DIMENSION(:), INTENT(IN)    :: PH_ROAD           ! sensible heat flux over road
REAL, DIMENSION(:), INTENT(IN)    :: PLE_ROAD          ! latent heat flux over road
REAL, DIMENSION(:), INTENT(IN)    :: PGFLUX_ROAD       ! flux through the road
REAL, DIMENSION(:), INTENT(IN)    :: PRN_GARDEN        ! net radiation over green areas
REAL, DIMENSION(:), INTENT(IN)    :: PH_GARDEN         ! sensible heat flux over green areas
REAL, DIMENSION(:), INTENT(IN)    :: PLE_GARDEN        ! latent heat flux over green areas
REAL, DIMENSION(:), INTENT(IN)    :: PGFLUX_GARDEN     ! flux through the green areas
REAL, DIMENSION(:), INTENT(IN)    :: PRN_WALL          ! net radiation over wall
REAL, DIMENSION(:), INTENT(IN)    :: PH_WALL           ! sensible heat flux over wall
REAL, DIMENSION(:), INTENT(IN)    :: PLE_WALL          ! latent heat flux over wall
REAL, DIMENSION(:), INTENT(IN)    :: PGFLUX_WALL       ! flux through the wall
!
REAL, DIMENSION(:), INTENT(IN)    :: PLEW_ROOF         ! latent heat flux of snowfree roof
REAL, DIMENSION(:), INTENT(IN)    :: PLESNOW_ROOF      ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PLEW_ROAD         ! latent heat flux of snowfree road
REAL, DIMENSION(:), INTENT(IN)    :: PLESNOW_ROAD      ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PHSNOW_ROAD       ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PEVAP_GARDEN      ! evaporation over gardens
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_GRND          ! net radiation over ground
REAL, DIMENSION(:), INTENT(OUT)   :: PH_GRND           ! sensible heat flux over ground
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_GRND          ! latent heat flux over ground
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_GRND       ! flux through the ground
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_TOWN          ! net radiation over town
REAL, DIMENSION(:), INTENT(OUT)   :: PH_TOWN           ! sensible heat flux over town
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_TOWN          ! latent heat flux over town
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_TOWN       ! flux through the ground for town
REAL, DIMENSION(:), INTENT(OUT)   :: PEVAP_TOWN        ! evaporation (kg/m2/s)
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PLW_RAD))    :: ZLW_UP            ! upwards radiations
REAL, DIMENSION(SIZE(PROAD)) :: ZROAD, ZGARDEN
!
INTEGER :: JJ
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('AVG_URBAN_FLUXES',0,ZHOOK_HANDLE)
!
ZROAD(:)=0.
ZGARDEN(:)=0.
!
DO JJ=1,SIZE(PROAD)
!  
  IF (PROAD(JJ)+PGARDEN(JJ).NE.0.) THEN
    ZROAD(JJ)  = PROAD(JJ)  / (PROAD(JJ)+PGARDEN(JJ))
    ZGARDEN(JJ) =  PGARDEN(JJ) / (PROAD(JJ)+PGARDEN(JJ))
  ELSE
    ZROAD(JJ)=0.
    ZGARDEN(JJ)=0.
  ENDIF
!
!*      1.     Averaged fluxes for ground (green areas + road)
!              -----------------------------------------------
!
  PRN_GRND(JJ)    = ZROAD  (JJ) * PRN_ROAD  (JJ)  +  ZGARDEN(JJ) * PRN_GARDEN(JJ)
!
  PH_GRND (JJ)    = ZROAD  (JJ) * PH_ROAD   (JJ)  +  ZGARDEN(JJ) * PH_GARDEN (JJ)       &
                 +  PH_TRAFFIC(JJ)  
!
  PLE_GRND(JJ)    = ZROAD  (JJ) * PLE_ROAD  (JJ)  +  ZGARDEN(JJ) * PLE_GARDEN(JJ)       &
                 +  PLE_TRAFFIC(JJ)  
!
  PGFLUX_GRND(JJ) = ZROAD  (JJ) * PGFLUX_ROAD  (JJ)  +  ZGARDEN(JJ) * PGFLUX_GARDEN(JJ)
!
!-------------------------------------------------------------------------------
!
!*      2.     Averaged fluxes JJ built + green areas
!              -------------------------------------
!
  PRN_TOWN(JJ)    = PTOTS_O_HORS(JJ) * ( &
                   PROOF_FRAC  (JJ) * PRN_ROOF    (JJ)                  &
                 + PROAD_FRAC  (JJ) * PRN_ROAD    (JJ)                  &
                 + PGARDEN_FRAC(JJ) * PRN_GARDEN  (JJ)                  &
                 + PWALL_FRAC  (JJ) * PRN_WALL    (JJ))  
!
  PH_TOWN (JJ)    = PTOTS_O_HORS(JJ) * ( &
                   PROOF_FRAC  (JJ) * PH_ROOF    (JJ)                   &
                 + PROAD_FRAC  (JJ) * PH_ROAD    (JJ)                   &
                 + PGARDEN_FRAC(JJ) * PH_GARDEN  (JJ)                   &
                 + PWALL_FRAC  (JJ) * PH_WALL    (JJ))                  &
                 + PH_TRAFFIC(JJ)                                                      &
                 + PH_INDUSTRY(JJ)  
!
  PLE_TOWN(JJ)    = PTOTS_O_HORS(JJ) * ( &
                   PROOF_FRAC  (JJ) * PLE_ROOF  (JJ)                     &
                 + PROAD_FRAC  (JJ) * PLE_ROAD  (JJ)                     &
                 + PGARDEN_FRAC(JJ) * PLE_GARDEN(JJ)                     &
                 + PWALL_FRAC  (JJ) * PLE_WALL  (JJ))                    &
                 + PLE_TRAFFIC (JJ)                                                      &
                 + PLE_INDUSTRY(JJ)  
!
  PGFLUX_TOWN(JJ)= PTOTS_O_HORS(JJ) * ( &
                  PROOF_FRAC  (JJ) * PGFLUX_ROOF  (JJ)                   &
                + PROAD_FRAC  (JJ) * PGFLUX_ROAD  (JJ)                   &
                + PGARDEN_FRAC(JJ) * PGFLUX_GARDEN(JJ)                   &
                + PWALL_FRAC  (JJ) * PGFLUX_WALL  (JJ))   
!
!-------------------------------------------------------------------------------
!
!*      3.     Infra-red Radiative properties
!              ------------------------------
!
!*      3.1    Upward IR radiation for town
!              ----------------------------
!
  ZLW_UP(JJ) = PLW_RAD      (JJ)                                    &
            - ( PROOF_FRAC (JJ)*PDF_ROOF (JJ)*PABS_LW_ROOF      (JJ) &
               +PROOF_FRAC (JJ)*PDN_ROOF (JJ)*PABS_LW_SNOW_ROOF (JJ) &
               +PROAD_FRAC (JJ)*PDF_ROAD (JJ)*PABS_LW_ROAD      (JJ) &
               +PROAD_FRAC (JJ)*PDN_ROAD (JJ)*PABS_LW_SNOW_ROAD (JJ) &
               +PGARDEN_FRAC(JJ)*             PABS_LW_GARDEN     (JJ) &
               +PWALL_FRAC (JJ)             *PABS_LW_WALL      (JJ) &
              )*PTOTS_O_HORS(JJ)  
!
!*      3.2    Town emissivity
!              ---------------
!
  PEMIS_TOWN(JJ) =   PBLD      (JJ) * PDF_ROOF (JJ)*PEMIS_ROOF(JJ)                    &
                  + PBLD      (JJ) * PDN_ROOF (JJ)*PESNOW_ROOF(JJ)                   &
                  + PROAD     (JJ) * PDF_ROAD (JJ)                * PLW_S_TO_R (JJ)  &
                  + PROAD     (JJ) * PDN_ROAD (JJ)                * PLW_S_TO_NR(JJ)  &
                  + PGARDEN    (JJ)                               * PLW_S_TO_G (JJ)  &
                  + PWALL_FRAC(JJ) * PTOTS_O_HORS(JJ)             * PLW_S_TO_W (JJ)  
!
!*      3.3    Town radiative surface temperature
!              ----------------------------------
!
  PTS_TOWN(JJ)   = ((ZLW_UP(JJ) - PLW_RAD(JJ)*(1.-PEMIS_TOWN(JJ))) /PEMIS_TOWN(JJ)/XSTEFAN)**0.25
!
!-------------------------------------------------------------------------------
!
!*      4.     Averaged evaporative flux (kg/m2/s)
!              -----------------------------------
!
  PEVAP_TOWN(JJ) = PTOTS_O_HORS(JJ) * (  &
                  PROOF_FRAC (JJ) * (PDF_ROOF(JJ) *PLEW_ROOF(JJ) /XLVTT+PDN_ROOF(JJ) *PLESNOW_ROOF(JJ) /XLSTT)  &
                + PROAD_FRAC (JJ) * (PDF_ROAD(JJ) *PLEW_ROAD(JJ) /XLVTT+PDN_ROAD(JJ) *PLESNOW_ROAD(JJ) /XLSTT)  &
                + PGARDEN_FRAC(JJ) * PEVAP_GARDEN(JJ)                                               &
                + PWALL_FRAC (JJ) *  (PLE_WALL (JJ)/XLVTT))                                     &
                + (PLE_TRAFFIC(JJ) /XLVTT)                                                      &
                + (PLE_INDUSTRY(JJ)/XLVTT)  
!
!-------------------------------------------------------------------------------
!
!*      5.    Air canyon temperature at time t+dt
!             -----------------------------------
!
  PT_CANYON(JJ) = (  PTS_ROAD    (JJ) * PAC_ROAD (JJ) * PDF_ROAD (JJ) * ZROAD (JJ)  &
                  + PTS_GARDEN   (JJ) * PAC_GARDEN(JJ)                * ZGARDEN(JJ)   &
                  + PTS_WALL    (JJ) * PAC_WALL (JJ)                * PWALL_O_GRND(JJ)                  &
                  + PTA         (JJ) * PAC_TOP  (JJ)                                                   &
                  + PH_TRAFFIC  (JJ) / (1.-PBLD (JJ)) / PRHOA(JJ) / XCPD                                &
                  + PHSNOW_ROAD (JJ) * PDN_ROAD (JJ)  / PRHOA(JJ) / XCPD  )                             &
               / (                    PAC_ROAD (JJ) * PDF_ROAD (JJ) * ZROAD (JJ)  &
                  +                   PAC_GARDEN(JJ)                * ZGARDEN(JJ) &
                  +                   PAC_WALL (JJ) * PWALL_O_GRND(JJ)                                 &
                  +                   PAC_TOP  (JJ)                     )  
!
!-------------------------------------------------------------------------------
!
!*      6.     Air canyon specific humidity
!              ----------------------------
!
  PQ_CANYON(JJ) = (  PQSAT_ROAD   (JJ) * PAC_ROAD_WAT (JJ) * PDF_ROAD(JJ)      * ZROAD (JJ) * PDELT_ROAD(JJ)  &
                  + PQSAT_GARDEN  (JJ) * PAC_AGG_GARDEN(JJ) * PHU_AGG_GARDEN(JJ) * ZGARDEN(JJ)                  &
                  + PQA          (JJ) * PAC_TOP(JJ)                                                                            &
                  + PLE_TRAFFIC  (JJ) / (1.-PBLD(JJ)) / PRHOA(JJ) / XLVTT                                                       &
                  + PLESNOW_ROAD (JJ) * PDN_ROAD(JJ)  / PRHOA(JJ) / XLVTT     * ZROAD (JJ)                ) &
               / (  PAC_ROAD_WAT (JJ) * PDF_ROAD(JJ) * PDELT_ROAD(JJ)         * ZROAD (JJ)                  &
                  + PAC_AGG_GARDEN(JJ) * PHU_AGG_GARDEN(JJ)                    * ZGARDEN(JJ)                  &
                  + PAC_TOP(JJ)                                                                                             )  
!
ENDDO
IF (LHOOK) CALL DR_HOOK('AVG_URBAN_FLUXES',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE AVG_URBAN_FLUXES
