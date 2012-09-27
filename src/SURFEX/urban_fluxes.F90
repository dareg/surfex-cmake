!     #########
    SUBROUTINE URBAN_FLUXES(HIMPLICIT_WIND,                                     &
                       PT_CANYON, PQ_CANYON,                                    &
                       PT_LOWCAN, PQ_LOWCAN,                                    &
                       PTS_ROOF,PTS_ROAD,PTS_WALL,                              &
                       PPEW_A_COEF, PPEW_B_COEF,                                &
                       PEXNS, PEXNA, PTA, PQA, PRHOA,                           &
                       PVMOD,                                                   &
                       PH_TRAFFIC, PLE_TRAFFIC, PH_INDUSTRY, PLE_INDUSTRY,      &
                       PBLD, PROAD, PWALL_O_HOR,                                &
                       PEMIS_ROOF,                                              &
                       PABS_SW_ROOF, PABS_LW_ROOF, PABS_SW_WALL, PABS_LW_WALL,  &
                       PABS_SW_ROAD, PABS_LW_ROAD,                              &
                       PAC_ROOF, PAC_ROOF_WAT,                                  &
                       PAC_WALL, PAC_ROAD, PAC_ROAD_WAT,                        &
                       PCD,                                                     &
                       PQSAT_ROOF, PQSAT_ROAD,                                  &
                       PDELT_ROOF, PDELT_ROAD,                                  &
                       PROOF_FRAC, PWALL_FRAC, PROAD_FRAC,                      &
                       PTOTS_O_HORS,                                            &
                       PDF_ROOF, PDN_ROOF, PDF_ROAD, PDN_ROAD,                  &
                       PRN_ROOF, PH_ROOF, PLE_ROOF, PLEW_ROOF, PGFLUX_ROOF,     &
                       PRN_ROAD, PH_ROAD, PLE_ROAD, PLEW_ROAD, PGFLUX_ROAD,     &
                       PRN_WALL, PH_WALL, PLE_WALL, PGFLUX_WALL,                &
                       PRNSNOW_ROOF, PHSNOW_ROOF, PLESNOW_ROOF, PGSNOW_ROOF,    &
                       PRNSNOW_ROAD, PHSNOW_ROAD, PLESNOW_ROAD, PGSNOW_ROAD,    &
                       PRN_BLT,  PH_BLT,  PLE_BLT,  PGFLUX_BLT,                 &
                       PUSTAR_TOWN, PDQS_WALL, PDQS_ROAD, PDQS_ROOF, PQF_BLD,   &
                       PQF_BLDWFR, PDQS_TOWN, PQF_TOWN, PQF_WALL, PQF_ROOF,     &
                       PFLX_BLD_WALL, PFLX_BLD_ROOF, PFLX_BLD, PMELT_ROOF,      &
                       PDQS_SNOW_ROOF, PMELT_ROAD, PMELT_BLT                    )  
!   ##########################################################################
!
!!****  *URBAN_FLUXES* computes fluxes on urbanized surfaces  
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
!!      Modified    09/2012 : B. Decharme New wind implicitation
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,ONLY : XCPD, XLVTT, XLSTT, XSTEFAN, XLMTT
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
CHARACTER(LEN=*),     INTENT(IN)  :: HIMPLICIT_WIND   ! wind implicitation option
!                                                     ! 'OLD' = direct
!                                                     ! 'NEW' = Taylor serie, order 1
!
REAL, DIMENSION(:), INTENT(IN)    :: PT_CANYON    ! canyon air temperature
REAL, DIMENSION(:), INTENT(IN)    :: PQ_CANYON    ! canyon air specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PT_LOWCAN    ! low canyon air temperature
REAL, DIMENSION(:), INTENT(IN)    :: PQ_LOWCAN    ! low canyon air specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PTS_ROOF     ! roof surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_ROAD     ! road surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_WALL     ! wall surface temperature

REAL, DIMENSION(:), INTENT(IN)    :: PPEW_A_COEF  ! implicit coefficients
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_B_COEF  ! for wind coupling
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS        ! surface exner function
REAL, DIMENSION(:), INTENT(IN)    :: PEXNA        ! exner function

REAL, DIMENSION(:), INTENT(IN)    :: PTA          ! temperature at roof level
REAL, DIMENSION(:), INTENT(IN)    :: PQA          ! specific humidity
                                                  ! at roof level
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD        ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA        ! air density
                                                  ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PH_TRAFFIC   ! anthropogenic sensible
!                                                 ! heat fluxes due to traffic
REAL, DIMENSION(:), INTENT(IN)    :: PLE_TRAFFIC  ! anthropogenic latent
!                                                 ! heat fluxes due to traffic
REAL, DIMENSION(:), INTENT(IN)    :: PH_INDUSTRY  ! anthropogenic sensible
!                                                 ! heat fluxes due to factories
REAL, DIMENSION(:), INTENT(IN)    :: PLE_INDUSTRY ! anthropogenic latent
!                                                 ! heat fluxes due to factories
REAL, DIMENSION(:), INTENT(IN)    :: PBLD         ! fraction of buildings
REAL, DIMENSION(:), INTENT(IN)    :: PROAD        ! fraction of roads
REAL, DIMENSION(:), INTENT(IN)    :: PWALL_O_HOR  ! wall Surf. / (bld+road+green) Surf.
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_ROOF   ! roof emissivity
! 
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_ROOF ! absorbed SW rad. by roof
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_ROOF ! absorbed LW rad. by roof
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_WALL ! absorbed SW rad. by wall
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_WALL ! absorbed LW rad. by wall
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_ROAD ! absorbed SW rad. by road
REAL, DIMENSION(:), INTENT(IN)    :: PABS_LW_ROAD ! absorbed LW rad. by road
!
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROOF     ! surface conductance
!                                                 ! for heat transfers
!                                                 ! above roofs
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROOF_WAT ! surface conductance
!                                                 ! for heat transfers
!                                                 ! above roofs (for water)
REAL, DIMENSION(:), INTENT(IN)    :: PAC_WALL     ! surface conductance
!                                                 ! for heat transfer
!                                                 ! between wall and canyon
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROAD     ! surface conductance
!                                                 ! for heat transfers
!                                                 ! between road and canyon
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROAD_WAT ! surface conductance
!                                                 ! for heat transfers
!                                                 ! inside canyon (for water)
REAL, DIMENSION(:), INTENT(IN)    :: PCD          ! drag coefficient
REAL, DIMENSION(:), INTENT(IN)    :: PQSAT_ROOF   ! q_sat(Ts)
REAL, DIMENSION(:), INTENT(IN)    :: PQSAT_ROAD   ! q_sat(Ts)
REAL, DIMENSION(:), INTENT(IN)    :: PDELT_ROOF   ! water fraction on snow-free
REAL, DIMENSION(:), INTENT(IN)    :: PDELT_ROAD   ! roof and roads
REAL, DIMENSION(:), INTENT(IN)    :: PROOF_FRAC   ! roof, wall,
REAL, DIMENSION(:), INTENT(IN)    :: PWALL_FRAC   ! road, and green area
REAL, DIMENSION(:), INTENT(IN)    :: PROAD_FRAC   ! fractions
REAL, DIMENSION(:), INTENT(IN)    :: PTOTS_O_HORS ! total canyon+roof surface
!                                                 ! over horizontal surface
REAL, DIMENSION(:), INTENT(IN)    :: PDF_ROOF     ! snow-free    roof fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROOF     ! snow-covered roof fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDF_ROAD     ! snow-free    road fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROAD     ! snow-covered road fraction
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_ROOF     ! net radiation over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PH_ROOF      ! sensible heat flux over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_ROOF     ! latent heat flux over roof
REAL, DIMENSION(:), INTENT(OUT)   :: PLEW_ROOF    ! latent heat flux over snow-free roof
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_ROOF  ! flux through the roof
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_ROAD     ! net radiation over road
REAL, DIMENSION(:), INTENT(OUT)   :: PH_ROAD      ! sensible heat flux over road
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_ROAD     ! latent heat flux over road
REAL, DIMENSION(:), INTENT(OUT)   :: PLEW_ROAD    ! latent heat flux of snow-free road
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_ROAD  ! flux through the road
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_WALL     ! net radiation over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PH_WALL      ! sensible heat flux over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_WALL     ! latent heat flux over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_WALL  ! flux through the wall
!
REAL, DIMENSION(:), INTENT(IN)    :: PRNSNOW_ROOF ! net radiation over snow
REAL, DIMENSION(:), INTENT(IN)    :: PHSNOW_ROOF  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PLESNOW_ROOF ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PGSNOW_ROOF  ! flux under the snow
REAL, DIMENSION(:), INTENT(IN)    :: PRNSNOW_ROAD ! net radiation over snow
REAL, DIMENSION(:), INTENT(IN)    :: PHSNOW_ROAD  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PLESNOW_ROAD ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PGSNOW_ROAD  ! flux under the snow
!
REAL, DIMENSION(:), INTENT(OUT)   :: PRN_BLT      ! net radiation over built covers
REAL, DIMENSION(:), INTENT(OUT)   :: PH_BLT       ! sensible heat flux over built covers
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_BLT      ! latent heat flux over built covers
REAL, DIMENSION(:), INTENT(OUT)   :: PGFLUX_BLT   ! flux through the built covers
REAL, DIMENSION(:), INTENT(OUT)   :: PUSTAR_TOWN  ! friction velocity over town
REAL, DIMENSION(:), INTENT(IN)    :: PDQS_ROOF    ! storage inside roofs
REAL, DIMENSION(:), INTENT(IN)    :: PDQS_WALL    ! storage inside walls
REAL, DIMENSION(:), INTENT(IN)    :: PDQS_ROAD    ! storage inside roads
REAL, DIMENSION(:), INTENT(IN)    :: PQF_WALL     ! heating through wall 
REAL, DIMENSION(:), INTENT(IN)    :: PQF_ROOF     ! heating through roof
REAL, DIMENSION(:), INTENT(IN)    :: PFLX_BLD_ROOF! heat flx from bld to roof
REAL, DIMENSION(:), INTENT(IN)    :: PFLX_BLD_WALL! heat flx from bld to wall
REAL, DIMENSION(:), INTENT(OUT)   :: PQF_BLD      ! anthropogenic heat due to domestic heating
REAL, DIMENSION(:), INTENT(OUT)   :: PQF_BLDWFR   ! 2nd version of Qf
REAL, DIMENSION(:), INTENT(OUT)   :: PFLX_BLD     ! heat flx from bld to its structure
REAL, DIMENSION(:), INTENT(OUT)   :: PDQS_TOWN    ! storage inside town materials
REAL, DIMENSION(:), INTENT(OUT)   :: PQF_TOWN     ! total anthropogenic heat
!
REAL, DIMENSION(:), INTENT(IN)    :: PMELT_ROOF   ! snow melting on roof
REAL, DIMENSION(:), INTENT(IN)    :: PMELT_ROAD   ! snow melting on road
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT_BLT    ! snow melting for town
!
REAL, DIMENSION(:), INTENT(IN)    :: PDQS_SNOW_ROOF ! heat storage in roof snowpack
!
!*      0.2    declarations of local variables
!
! local variables for built area seb
!
REAL, DIMENSION(SIZE(PTA))  :: ZTOTS_O_BLT
REAL, DIMENSION(SIZE(PTA))  :: ZROOF_BLT
REAL, DIMENSION(SIZE(PTA))  :: ZROAD_BLT
REAL, DIMENSION(SIZE(PTA))  :: ZWALL_BLT
!
REAL, DIMENSION(SIZE(PTA))  :: ZH_ROOF_SNOWFREE
REAL, DIMENSION(SIZE(PTA))  :: ZRN_ROOF_SNOWFREE
!
REAL, DIMENSION(SIZE(PTA))  :: ZUSTAR2 ! square of friction velocity (m2/s2)
REAL, DIMENSION(SIZE(PTA))  :: ZVMOD   ! Wind
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*      1.     Fluxes at snow-free roofs
!              -------------------------
!
!                                            net radiation
!
IF (LHOOK) CALL DR_HOOK('URBAN_FLUXES',0,ZHOOK_HANDLE)
ZRN_ROOF_SNOWFREE(:) = PABS_SW_ROOF(:) + PABS_LW_ROOF(:)
!
!                                            sensible heat flux
!
ZH_ROOF_SNOWFREE(:) =(PTS_ROOF(:) - PTA(:) * PEXNS(:) / PEXNA(:)) &
                        * PAC_ROOF(:) * PRHOA(:) * XCPD / PEXNS(:)  
!
!                                            latent heat of evaporation from
!                                            the roof
!
PLEW_ROOF(:) =(PQSAT_ROOF(:) - PQA(:))                     &
               * PAC_ROOF_WAT(:) * PDELT_ROOF(:) * XLVTT * PRHOA(:)  
!
!-------------------------------------------------------------------------------
!
!*      2.     Fluxes at snow-free roads
!              -------------------------
!
!                                            net radiation
!
PRN_ROAD(:) = PABS_SW_ROAD(:) + PABS_LW_ROAD(:)
!
!                                            sensible heat flux
!
PH_ROAD(:) =  (PTS_ROAD(:) - PT_LOWCAN(:) * PEXNS(:) / PEXNA(:))     &
                 * PAC_ROAD(:)  * XCPD * PRHOA(:) / PEXNS(:)  
!
!                                            latent heat of evaporation from
!                                            the road
!
PLEW_ROAD(:) = ( PQSAT_ROAD(:) - PQ_LOWCAN(:) )                   &
               *   PAC_ROAD_WAT(:) * PDELT_ROAD(:) * XLVTT * PRHOA(:)  
!
!-------------------------------------------------------------------------------
!
!*      3.     Fluxes at walls
!              ---------------
!
!                                            net radiation
!
PRN_WALL(:) = PABS_SW_WALL(:) + PABS_LW_WALL(:)
!
!                                            sensible heat flux
!
PH_WALL(:) =  (PTS_WALL(:) - PT_CANYON(:) * PEXNS(:) / PEXNA(:))     &
                 * PAC_WALL(:) * XCPD * PRHOA(:) / PEXNS(:)  
!
!                                            latent heat of evaporation from
!                                            the ground
!
PLE_WALL(:) = 0.
!
!                                            heat flux into the ground
!
PGFLUX_WALL(:) = PRN_WALL(:) - PH_WALL(:) - PLE_WALL(:)
!
!-------------------------------------------------------------------------------
!
!*      4.     Snow-free and snow-covered surfaces averaging
!              ---------------------------------------------
!
!*      4.1    Roads
!              -----
!
!                                            heat flux into the ground
!
!
PGFLUX_ROAD (:) =  PDF_ROAD(:) * (PRN_ROAD(:) - PH_ROAD (:) - PLEW_ROAD(:) ) &
                   + PDN_ROAD(:) *  PGSNOW_ROAD(:)  
!
!
!                                            net radiation
!
PRN_ROAD(:) = PRN_ROAD(:) * PDF_ROAD(:) + PRNSNOW_ROAD(:) * PDN_ROAD(:)
!
!                                            sensible heat flux
!
PH_ROAD (:) = PH_ROAD (:) * PDF_ROAD(:) + PHSNOW_ROAD(:)  * PDN_ROAD(:)
!
!                                            total latent heat of evaporation from
!                                            the road (snow free + snow)
!
PLE_ROAD (:) = PLEW_ROAD    (:) * PDF_ROAD(:) +  PLESNOW_ROAD(:) * PDN_ROAD(:)
!
!*      4.2    Roofs
!              -----
!
!                                            heat flux into the ground
!
!
PGFLUX_ROOF (:) =  PDF_ROOF(:) * (ZRN_ROOF_SNOWFREE(:) - ZH_ROOF_SNOWFREE (:) - PLEW_ROOF(:) ) &
                   + PDN_ROOF(:) *  PGSNOW_ROOF(:)  
!
!
!                                            net radiation
!
PRN_ROOF(:) = ZRN_ROOF_SNOWFREE(:) * PDF_ROOF(:) + PRNSNOW_ROOF(:) * PDN_ROOF(:)
!
!                                            sensible heat flux
!
PH_ROOF (:) = ZH_ROOF_SNOWFREE (:) * PDF_ROOF(:) + PHSNOW_ROOF(:)  * PDN_ROOF(:)
!
!                                            total latent heat of evaporation from
!                                            the roof (snow free + snow)
!
PLE_ROOF (:) = PLEW_ROOF   (:) * PDF_ROOF(:)  + PLESNOW_ROOF(:) * PDN_ROOF(:)
!
!-------------------------------------------------------------------------------
!
!*      5.     Momentum fluxes
!              ---------------
!
ZUSTAR2(:) = 0.0
ZVMOD  (:) = PVMOD(:)
!
IF(HIMPLICIT_WIND=='OLD')THEN 
! old implicitation
  ZUSTAR2(:) = (PCD(:)*PVMOD(:)*PPEW_B_COEF(:))/    &
               (1.0-PRHOA(:)*PCD(:)*PVMOD(:)*PPEW_A_COEF(:))
ELSE
! new implicitation
  ZUSTAR2(:) = (PCD(:)*PVMOD(:)*(2.*PPEW_B_COEF(:)-PVMOD(:)))/ &
               (1.0-2.0*PRHOA(:)*PCD(:)*PVMOD(:)*PPEW_A_COEF(:))
!                   
  ZVMOD(:) = PRHOA(:)*PPEW_A_COEF(:)*ZUSTAR2(:) + PPEW_B_COEF(:)
  ZVMOD(:) = MAX(ZVMOD(:),0.)
!
  WHERE(PPEW_A_COEF(:)/= 0.)
        ZUSTAR2(:) = MAX( ( ZVMOD(:) - PPEW_B_COEF(:) ) / (PRHOA(:)*PPEW_A_COEF(:)), 0.)
  ENDWHERE
!               
ENDIF
!
PUSTAR_TOWN(:) = SQRT(ZUSTAR2(:))
!
!-------------------------------------------------------------------------------
!
!*      6.     Averaged fluxes
!              ---------------
!
!*      6.1    Built fraction
!              --------------
!
!
WHERE (PBLD(:) .GT. 0. .OR. PROAD(:) .GT. 0.)
 ZTOTS_O_BLT(:) = PWALL_O_HOR (:) / (PBLD(:)+PROAD(:)) + 1.
 ZROOF_BLT(:)   = PBLD        (:) / (PBLD(:)+PROAD(:)) / ZTOTS_O_BLT(:)
 ZROAD_BLT(:)   = PROAD       (:) / (PBLD(:)+PROAD(:)) / ZTOTS_O_BLT(:)
 ZWALL_BLT(:)   = PWALL_O_HOR (:) / (PBLD(:)+PROAD(:)) / ZTOTS_O_BLT(:)
ELSEWHERE
 ZTOTS_O_BLT(:) = 0.
 ZROOF_BLT(:)   = 0.
 ZROAD_BLT(:)   = 0.
 ZWALL_BLT(:)   = 0.
ENDWHERE
!
PRN_BLT (:)    =  ZROOF_BLT  (:) * PRN_ROOF    (:) * ZTOTS_O_BLT (:)  &
                 +  ZROAD_BLT  (:) * PRN_ROAD    (:) * ZTOTS_O_BLT (:)  &
                 +  ZWALL_BLT  (:) * PRN_WALL    (:) * ZTOTS_O_BLT (:)  
!
PH_BLT  (:)    =  ZROOF_BLT  (:) * PH_ROOF    (:) * ZTOTS_O_BLT (:)   &
                 +  ZROAD_BLT  (:) * PH_ROAD    (:) * ZTOTS_O_BLT (:)   &
                 +  ZWALL_BLT  (:) * PH_WALL    (:) * ZTOTS_O_BLT (:)   &
                 +  PH_TRAFFIC(:)                                       &
                 +  PH_INDUSTRY(:)  
!
PLE_BLT (:)    =  ZROOF_BLT  (:) * PLE_ROOF (:) * ZTOTS_O_BLT (:) &
                 +  ZROAD_BLT  (:) * PLE_ROAD (:) * ZTOTS_O_BLT (:) &
                 +  ZWALL_BLT  (:) * PLE_WALL (:) * ZTOTS_O_BLT (:) &
                 +  PLE_TRAFFIC(:)                                  &
                 +  PLE_INDUSTRY(:)  
!
PGFLUX_BLT (:) =  ZROOF_BLT  (:) * PGFLUX_ROOF (:) * ZTOTS_O_BLT (:)  &
                 +  ZROAD_BLT  (:) * PGFLUX_ROAD (:) * ZTOTS_O_BLT (:)  &
                 +  ZWALL_BLT  (:) * PGFLUX_WALL (:) * ZTOTS_O_BLT (:)  
!
PMELT_BLT (:) = PROOF_FRAC(:) * PMELT_ROOF(:) * PDN_ROOF(:) * PTOTS_O_HORS(:) &
                + PROAD_FRAC(:) * PMELT_ROAD(:) * PDN_ROAD(:) * PTOTS_O_HORS(:)  
!
PDQS_TOWN(:)  = PROOF_FRAC (:) * PDQS_ROOF (:) * PTOTS_O_HORS(:)                    &
                + PROAD_FRAC (:) * PDQS_ROAD (:) * PTOTS_O_HORS(:)                    &
                + PWALL_FRAC (:) * PDQS_WALL (:) * PTOTS_O_HORS(:)  
!
!PQF_BLD(:)    = PROOF_FRAC(:) * ( PH_ROOF(:)+PLE_ROOF(:)+PDQS_ROOF(:)-PRN_ROOF(:) ) &
!              * PTOTS_O_HORS(:)                                                     &
!              + PWALL_FRAC(:) * ( PH_WALL(:)+PLE_WALL(:)+PDQS_WALL(:)-PRN_WALL(:) ) &
!              * PTOTS_O_HORS(:)
!
PQF_BLD(:)= PROOF_FRAC(:) * PTOTS_O_HORS(:)                                                        &
              * ( (ZH_ROOF_SNOWFREE(:)+PLEW_ROOF(:)+PDQS_ROOF(:)-ZRN_ROOF_SNOWFREE(:))*PDF_ROOF(:)   &
                + (PDQS_ROOF(:)-PGSNOW_ROOF(:))*PDN_ROOF(:) )                                        &
            + PWALL_FRAC(:) * PTOTS_O_HORS(:)                                                        &
              * (PH_WALL(:)+PLE_WALL(:)+PDQS_WALL(:)-PRN_WALL(:))  
!
PQF_BLDWFR(:)=(PROOF_FRAC(:)*PQF_ROOF(:)+PWALL_FRAC(:)*PQF_WALL(:))*PTOTS_O_HORS(:)
!
PQF_TOWN(:)=PQF_BLD(:)+PH_TRAFFIC(:)+PH_INDUSTRY(:)
!
PFLX_BLD(:)=(PROOF_FRAC(:)*PFLX_BLD_ROOF(:)+ &
               PWALL_FRAC(:)*PFLX_BLD_WALL(:))*PTOTS_O_HORS(:)                
IF (LHOOK) CALL DR_HOOK('URBAN_FLUXES',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE URBAN_FLUXES
