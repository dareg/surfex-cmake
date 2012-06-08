!   ##########################################################################
    SUBROUTINE ROAD_LAYER_E_BUDGET(PT_ROAD, PTS_WALL, PQSAT_ROAD, &
                                      PT_LOWCAN, PQ_LOWCAN,       &
                                      PTS_GARDEN,                 &
                                      PTA, PQA, PPS,              &
                                      PLW_RAD,  PTSTEP,           &
                                      PHC_ROAD,PTC_ROAD,PD_ROAD,  &
                                      PDELT_ROAD, PDN_ROAD,       &
                                      PTSNOW_ROAD, PGSNOW_ROAD,   &
                                      PRHOA,                      &
                                      PAC_ROAD, PAC_ROAD_WAT,     &
                                      PABS_SW_ROAD, PABS_LW_ROAD, &
                                      PLW_W_TO_R, PLW_R_TO_R,     &
                                      PLW_G_TO_R, PLW_S_TO_R,     &
                                      PLW_NR_TO_R, PDQS_ROAD      )
!   ##########################################################################
!
!!****  *ROAD_LAYER_E_BUDGET*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the evoultion of roads surface temperatures
!         
!     
!!**  METHOD
!     ------
!
!    6 : equations for evolution of Ts_road 
!        **********************************
!
!
!     dTr_1(t) / dt = 1/(dr_1*Cr_1) * (  Rn_r - H_r - LE_r 
!                                      - 2*Kr_1*(Tr_1-Tr_2)/(dr_1 +dr_2)       )
!
!     dTr_k(t) / dt = 1/(dr_k*Cr_k) * (- 2*Kr_k-1*(Tr_k-Tr_k-1)/(dr_k-1 +dr_k) 
!                                      - 2*Kr_k  *(Tr_k-Tr_k+1)/(dr_k+1 +dr_k) )
!
!       with
!
!   K*_k  = (d*_k+ d*_k+1)/(d*_k/k*_k+ d*_k+1/k*_k+1)
!
!   Rn_r = abs_Rg_r
!  - sigma * emis_r                                                   * Ts_r**4 (t+dt)
!  +         emis_r                       *    SVF_r                  * LWR
!  + sigma * emis_r * emis_w              * (1-SVF_r)                 * Ts_w**4 (t+dt)
!  +         emis_r            (1-emis_w) * (1-SVF_r)   *      SVF_w  * LWR
!  + sigma * emis_r * emis_w * (1-emis_w) * (1-SVF_r)   * (1-2*SVF_w) * Ts_w**4 (t+dt)
!  + sigma * emis_r * emis_r * (1-emis_w) * (1-SVF_r)   *      SVF_w  * Ts_r**4 (t+dt)
!
!  H_r  = rho Cp CH V ( Ts_r (t+dt) - Ta_canyon )
!
!  LE_r = rho Lv CH V ( qs_r (t+dt) - qa_canyon )
!
!
! The system is implicited (or semi-implicited).
!
! ZIMPL=1    ---> implicit system
! ZIMPL=0.5  ---> semi-implicit system
! ZIMPL=0    ---> explicit system
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
!!                  21/11/01 (V. Masson and A. Lemonsu) bug of latent flux
!!                           for very strong evaporation (all reservoir emptied
!!                           in one time-step)
!!                     02/11 (V. Masson) split of the routine for roads and walls separately
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,ONLY : XCPD, XLVTT
!
USE MODE_THERMOS
!
USE MODI_TRIDIAG_GROUND
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_ROAD    ! road layers temperatures
REAL, DIMENSION(:), INTENT(IN)    :: PTS_WALL     ! wall surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PQSAT_ROAD   ! q_sat(Ts)
REAL, DIMENSION(:), INTENT(IN)    :: PT_LOWCAN    ! low canyon air temperature
REAL, DIMENSION(:), INTENT(IN)    :: PQ_LOWCAN    ! and specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PTS_GARDEN   ! green area surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTA          ! atmospheric air temperature
REAL, DIMENSION(:), INTENT(IN)    :: PQA          ! and specific humidity at roof level
REAL, DIMENSION(:), INTENT(IN)    :: PPS          ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD      ! atmospheric infrared radiation
REAL,               INTENT(IN)    :: PTSTEP     ! time step
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_ROAD     ! heat capacity for road layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_ROAD     ! thermal conductivity for road layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_ROAD      ! depth of road layers
REAL, DIMENSION(:), INTENT(IN)    :: PDELT_ROAD   ! fraction of water
REAL, DIMENSION(:), INTENT(IN)    :: PDN_ROAD     ! road snow fraction
REAL, DIMENSION(:), INTENT(IN)    :: PTSNOW_ROAD  ! road snow temperature
REAL, DIMENSION(:), INTENT(IN)    :: PGSNOW_ROAD  ! road snow conduction
!                                                 ! heat fluxes at mantel
!                                                 ! base
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA        ! rho
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROAD     ! aerodynamical conductance
!                                                 ! between road and canyon
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROAD_WAT ! aerodynamical conductance
!                                                 ! between road and canyon
!                                                 ! (for water)
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_ROAD ! absorbed solar radiation
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_ROAD ! absorbed infrared rad.
!
REAL, DIMENSION(:), INTENT(IN)    :: PLW_W_TO_R    ! LW interactions wall  -> road 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_R_TO_R    ! LW interactions road  -> road 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_G_TO_R    ! LW interactions green -> road 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_R    ! LW interactions sky   -> road 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_NR_TO_R   ! LW interactions road(snow) -> road 
!
REAL, DIMENSION(:), INTENT(INOUT) :: PDQS_ROAD     !heat storage inside the road
!
!*      0.2    declarations of local variables
!
!
REAL :: ZIMPL=0.5      ! implicit coefficient
REAL :: ZEXPL=0.5      ! explicit coefficient
!
REAL, DIMENSION(SIZE(PTA),SIZE(PT_ROAD,2)) ::  ZA,& ! lower diag.
                                               ZB,& ! main  diag.
                                               ZC,& ! upper diag.
                                               ZY,& ! r.h.s.
                                               ZX   ! solution

!
REAL, DIMENSION(SIZE(PTA)) :: ZDN_ROAD    ! snow-covered surface fraction on road
REAL, DIMENSION(SIZE(PTA)) :: ZDF_ROAD    ! snow-free surface fraction on road
!
REAL, DIMENSION(SIZE(PTA)) :: ZDQSAT_ROAD ! dq_sat/dTs
REAL, DIMENSION(SIZE(PTA)) :: ZRHO_ACF_R  ! rho * conductance
!                                         !     * snow-free f.
REAL, DIMENSION(SIZE(PTA)) :: ZRHO_ACF_R_WAT ! rho * conductance for water
!                                         !     * snow-free f.
REAL, DIMENSION(SIZE(PTA),SIZE(PT_ROAD,2)) :: ZMTC_O_D_ROAD
! mean thermal conductivity over distance between 2 layers
REAL, DIMENSION(SIZE(PTA),SIZE(PT_ROAD,2)) :: ZHC_D_ROAD
! thermal capacity times layer depth
REAL, DIMENSION(SIZE(PTA)) :: ZTS_ROAD    ! road surface temperature
REAL, DIMENSION(SIZE(PT_ROAD,1)) :: ZEI_ROAD  ! internal energy of roads
REAL, DIMENSION(SIZE(PT_ROAD,1)) :: ZPEI_ROAD ! internal energy of roads at the
                                              ! following time step
!
INTEGER :: IROAD_LAYER           ! number of road layers
INTEGER :: JLAYER, JJ            ! loop counter
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ROAD_LAYER_E_BUDGET',0,ZHOOK_HANDLE)
!
!*      1.     Layer thermal properties
!              ------------------------
!
!*      1.1    Roads
!              -----
!
IROAD_LAYER = SIZE(PT_ROAD,2)
ZMTC_O_D_ROAD(:,:) = 0.
!
DO JLAYER=1,IROAD_LAYER-1
!
  DO JJ=1,SIZE(PT_ROAD,1)
!
    ZMTC_O_D_ROAD(JJ,JLAYER) = 2./(  PD_ROAD(JJ,JLAYER  )/PTC_ROAD(JJ,JLAYER  ) &
                                + PD_ROAD(JJ,JLAYER+1)/PTC_ROAD(JJ,JLAYER+1) )
    ZHC_D_ROAD   (JJ,JLAYER) = PHC_ROAD(JJ,JLAYER) * PD_ROAD (JJ,JLAYER)
  ENDDO
!
END DO
!
DO JJ=1,SIZE(PT_ROAD,1)
  ZMTC_O_D_ROAD(JJ,IROAD_LAYER) = 2. * PTC_ROAD(JJ,IROAD_LAYER) &
                                      / PD_ROAD (JJ,IROAD_LAYER)
  ZHC_D_ROAD   (JJ,IROAD_LAYER) = PHC_ROAD(JJ,IROAD_LAYER) &
                               * PD_ROAD (JJ,IROAD_LAYER)
ENDDO
!
!-------------------------------------------------------------------------------
!
!*      2.    Preliminaries
!             -------------
!
!
!*      2.1    snow-free surface fraction
!              --------------------------
!
ZDN_ROAD(:)=PDN_ROAD(:)
ZDF_ROAD(:)=1.-ZDN_ROAD(:)
!
!*      2.3    Surface temperatures
!              --------------------
!
ZTS_ROAD(:) = PT_ROAD(:,1)
!
!
!*      2.2    flux properties
!              ---------------
!
ZRHO_ACF_R(:) = PRHOA(:) * PAC_ROAD(:) * ZDF_ROAD(:)
ZRHO_ACF_R_WAT(:) = PRHOA(:) * PAC_ROAD_WAT(:) * ZDF_ROAD(:)
!
!*      2.4    qsat, dqsat/dTs, and humidity for roads
!              ---------------------------------------
!
ZDQSAT_ROAD(:) = DQSAT(ZTS_ROAD(:),PPS(:),PQSAT_ROAD(:))
!
!
ZEI_ROAD(:)=PHC_ROAD(:,1)*PD_ROAD(:,1)*PT_ROAD(:,1)
DO JLAYER=2,IROAD_LAYER
  DO JJ=1,SIZE(PT_ROAD,1)
     ZEI_ROAD(JJ)=ZEI_ROAD(JJ)+(PHC_ROAD(JJ,JLAYER)*PD_ROAD(JJ,JLAYER)* &
                        PT_ROAD(JJ,JLAYER))
  ENDDO
END DO
!
!-------------------------------------------------------------------------------
!
!*      3.     First road layers coefficients
!              ------------------------------
!
DO JJ=1,SIZE(PT_ROAD,1)

  ZA(JJ,1) =   0.

  ZB(JJ,1) =   ZHC_D_ROAD(JJ,1)/PTSTEP                             &
     + ZIMPL * ( - ZDF_ROAD(JJ) * 4.*ZTS_ROAD(JJ)**3 * PLW_R_TO_R(JJ)            &
                 + ZRHO_ACF_R(JJ) * XCPD                                       &
                 + ZRHO_ACF_R_WAT(JJ) * XLVTT * PDELT_ROAD(JJ) * ZDQSAT_ROAD(JJ) &
                 + ZMTC_O_D_ROAD(JJ,1)                                         &
               )

  ZC(JJ,1) =                                                  &
       ZIMPL * ( - ZMTC_O_D_ROAD(JJ,1)                                    &
                 )

  ZY(JJ,1) =   ZHC_D_ROAD(JJ,1)/PTSTEP*PT_ROAD(JJ,1)           &
                 + ZDF_ROAD(JJ) * PABS_SW_ROAD(JJ)                         &
                 + ZDF_ROAD(JJ) * PLW_RAD     (JJ)    * PLW_S_TO_R (JJ)     &
                 + ZDF_ROAD(JJ) * PTS_WALL    (JJ)**4 * PLW_W_TO_R (JJ)     &
                 + ZDF_ROAD(JJ) * ZTS_ROAD    (JJ)**4 * PLW_R_TO_R (JJ)     &
                 + ZDF_ROAD(JJ) * PTSNOW_ROAD (JJ)**4 * PLW_NR_TO_R(JJ)     &
                 + ZDF_ROAD(JJ) * PTS_GARDEN  (JJ)**4 * PLW_G_TO_R (JJ)     &
                 + ZRHO_ACF_R(JJ) * XCPD * PT_LOWCAN(JJ)                   &
                 + ZRHO_ACF_R_WAT(JJ) * XLVTT * PDELT_ROAD(JJ)             &
                         * (  PQ_LOWCAN(JJ) - PQSAT_ROAD(JJ)  )            &
                 + ZDN_ROAD(JJ) * PGSNOW_ROAD(JJ)                          &
     + ZIMPL * ( - ZDF_ROAD(JJ) * 4.*ZTS_ROAD(JJ)**4 * PLW_R_TO_R(JJ)       &
                 + ZRHO_ACF_R_WAT(JJ) * XLVTT * PDELT_ROAD(JJ)             &
                               * ZDQSAT_ROAD(JJ) * ZTS_ROAD(JJ)            &
               )                                                         &
     + ZEXPL * ( - ZRHO_ACF_R(JJ) * XCPD * ZTS_ROAD(JJ)                    &
                 - ZMTC_O_D_ROAD(JJ,1) * PT_ROAD(JJ,1)                     &
                 + ZMTC_O_D_ROAD(JJ,1) * PT_ROAD(JJ,2)                     &
               )
!-------------------------------------------------------------------------------
!
!*      8.     Last road layers coefficients
!              ------------------------------
!
  ZA(JJ,IROAD_LAYER) =                                     &
       ZIMPL * ( - ZMTC_O_D_ROAD(JJ,IROAD_LAYER-1)                     &
                 )

  ZB(JJ,IROAD_LAYER) =   ZHC_D_ROAD(JJ,IROAD_LAYER) / PTSTEP  &
     + ZIMPL * (   ZMTC_O_D_ROAD(JJ,IROAD_LAYER-1)                     &
               )

  ZC(JJ,IROAD_LAYER) =   0.
!
  ZY(JJ,IROAD_LAYER) =   ZHC_D_ROAD(JJ,IROAD_LAYER) / PTSTEP   &
                                             * PT_ROAD(JJ,IROAD_LAYER) &
     + ZEXPL * (   ZMTC_O_D_ROAD(JJ,IROAD_LAYER-1)                     &
                       * PT_ROAD(JJ,IROAD_LAYER-1)                     &
                 - ZMTC_O_D_ROAD(JJ,IROAD_LAYER-1)                     &
                       * PT_ROAD(JJ,IROAD_LAYER  )                     &
               )
     
ENDDO
!
!-------------------------------------------------------------------------------
!
!*      8.     Other road layers coefficients
!              ------------------------------
!
DO JLAYER=2,IROAD_LAYER-1

  DO JJ=1,SIZE(ZA,1)

    ZA(JJ,JLAYER) =                                                     &
         ZIMPL * ( - ZMTC_O_D_ROAD(JJ,JLAYER-1)                         &
                 )

    ZB(JJ,JLAYER) =   ZHC_D_ROAD(JJ,JLAYER)/PTSTEP                       &
       + ZIMPL * (   ZMTC_O_D_ROAD(JJ,JLAYER-1)                         &
                   + ZMTC_O_D_ROAD(JJ,JLAYER  )                         &
                 )

    ZC(JJ,JLAYER) =                                                     &
         ZIMPL * ( - ZMTC_O_D_ROAD(JJ,JLAYER  )                         &
                 )
!
    ZY(JJ,JLAYER) =   ZHC_D_ROAD(JJ,JLAYER)/PTSTEP*PT_ROAD(JJ,JLAYER)     &
       + ZEXPL * (   ZMTC_O_D_ROAD(JJ,JLAYER-1) * PT_ROAD(JJ,JLAYER-1)   &
                   - ZMTC_O_D_ROAD(JJ,JLAYER-1) * PT_ROAD(JJ,JLAYER  )   &
                   - ZMTC_O_D_ROAD(JJ,JLAYER  ) * PT_ROAD(JJ,JLAYER  )   &
                   + ZMTC_O_D_ROAD(JJ,JLAYER  ) * PT_ROAD(JJ,JLAYER+1)   &
                 )
  ENDDO
END DO
!
!-------------------------------------------------------------------------------
!
!*     10.     Tri-diagonal system resolution
!              ------------------------------
!
CALL TRIDIAG_GROUND(ZA,ZB,ZC,ZY,ZX)
!
DO JLAYER=1,IROAD_LAYER
  PT_ROAD(:,JLAYER) = ZX(:,JLAYER)
END DO  
! 
!-------------------------------------------------------------------------------
!
!*     12.    Road and wall absorbed infra-red radiation on snow-free surfaces
!             ----------------------------------------------------------------
!
!* radiative surface temperature used during the energy balance (linearized at t+)
!
DO JJ=1,SIZE(PTS_WALL)
!
  ZTS_ROAD(JJ) = (ZTS_ROAD(JJ)**4 + 4.*ZIMPL*ZTS_ROAD(JJ)**3 * (PT_ROAD(JJ,1) - ZTS_ROAD(JJ)))**0.25
!
!* absorbed LW
  PABS_LW_ROAD       (JJ) = PLW_S_TO_R (JJ)*PLW_RAD       (JJ)        &
                       + PLW_W_TO_R (JJ)*PTS_WALL      (JJ)**4       &
                       + PLW_R_TO_R (JJ)*ZTS_ROAD      (JJ)**4       &
                       + PLW_G_TO_R (JJ)*PTS_GARDEN    (JJ)**4       &
                       + PLW_NR_TO_R(JJ)*PTSNOW_ROAD   (JJ)**4
!
ENDDO
!
!-------------------------------------------------------------------------------
!
!
!*      13     heat storage inside roads
!              -------------------------
!
!       13.1   internal energy of the roads at the next time step
!              --------------------------------------------------
!
ZPEI_ROAD(:)=PHC_ROAD(:,1)*PD_ROAD(:,1)*PT_ROAD(:,1)
DO JLAYER=2,IROAD_LAYER
  ZPEI_ROAD(:)=ZPEI_ROAD(:)+(PHC_ROAD(:,JLAYER)*PD_ROAD(:,JLAYER)* &
                        PT_ROAD(:,JLAYER))
END DO
!
!       13.2   heat storage flux inside roads 
!              ------------------------------
!
PDQS_ROAD(:)=(ZPEI_ROAD(:)-ZEI_ROAD(:))/PTSTEP
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ROAD_LAYER_E_BUDGET',1,ZHOOK_HANDLE)
!
END SUBROUTINE ROAD_LAYER_E_BUDGET
