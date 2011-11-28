!   ##########################################################################
    SUBROUTINE WALL_LAYER_E_BUDGET(PTS_ROAD, PT_WALL,                  &
                                   PT_CANYON, PQ_CANYON,               &
                                   PTS_GARDEN,                         &
                                   PTA, PQA, PPS,                      &
                                   PLW_RAD,  PTSTEP,                   &
                                   PHC_WALL,PTC_WALL,PD_WALL,          &
                                   PTI_BLD, PTI_BLDWFR, PAC_BLD,       &
                                   PTSNOW_ROAD, PRHOA, PAC_WALL,       &
                                   PABS_SW_WALL, PABS_LW_WALL,         &
                                   PLW_W_TO_W, PLW_R_TO_W, PLW_G_TO_W, &
                                   PLW_S_TO_W, PLW_NR_TO_W,            &
                                   PDQS_WALL, PQF_WALL, PFLX_BLD_WALL  )
!   ##########################################################################
!
!!****  *ROAD_WALL_LAYER_E_BUDGET*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the evoultion of roads and walls surface temperatures
!         
!     
!!**  METHOD
!     ------
!
!    6 : equations for evolution of Ts_road and Ts_wall simultaneously
!        *************************************************************
!
!     dTw_k(t) / dt = 1/(dw_k*Cw_k) * (- 2*Kw_k-1*(Tw_k-Tw_k-1)/(dw_k-1 +dw_k) 
!                                      - 2*Kw_k  *(Tw_k-Tw_k+1)/(dw_k+1 +dw_k) )
!
!     dTw_1(t) / dt = 1/(dw_1*Cw_1) * (  Rn_w - H_w - LE_w 
!                                      - 2*Kw_1*(Tw_1-Tw_2)/(dw_1 +dw_2)       )
!
!
!       with
!
!   K*_k  = (d*_k+ d*_k+1)/(d*_k/k*_k+ d*_k+1/k*_k+1)
!
!   Rn_w = abs_Rg_w 
!  - sigma * emis_w                                                   * Ts_w**4 (t+dt)
!  +         emis_w                       *      SVF_w                * LWR
!  + sigma * emis_w * emis_r              *      SVF_w                * Ts_r**4 (t+dt)
!  + sigma * emis_w * emis_w              * (1-2*SVF_w)               * Ts_w**4 (t+dt)
!  +         emis_w            (1-emis_r) *      SVF_r  *      SVF_w  * LWR
!  +         emis_w            (1-emis_w) *      SVF_w  * (1-2*SVF_w) * LWR
!  + sigma * emis_w * emis_w * (1-emis_r) *      SVF_w  * (1-  SVF_r) * Ts_w**4 (t+dt)
!  + sigma * emis_w * emis_w * (1-emis_w) * (1-2*SVF_w) * (1-2*SVF_w) * Ts_w**4 (t+dt)
!  + sigma * emis_w * emis_r * (1-emis_w) *      SVF_w  * (1-2*SVF_w) * Ts_r**4 (t+dt)
!
!  H_w  = rho Cp CH V ( Ts_w (t+dt) - Ta_canyon )
!
!  LE_w = rho Lv CH V ( qs_w (t+dt) - qa_canyon )
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
!!                     02/11 (V. Masson) splits the routine for road and walls separately
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,ONLY : XCPD
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
REAL, DIMENSION(:), INTENT(IN)    :: PTS_ROAD     ! road surface temperature
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_WALL    ! wall layers temperatures
REAL, DIMENSION(:), INTENT(IN)    :: PT_CANYON    ! air canyon temperature
REAL, DIMENSION(:), INTENT(IN)    :: PQ_CANYON    ! and specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PTS_GARDEN   ! green area surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTA          ! atmospheric air temperature
REAL, DIMENSION(:), INTENT(IN)    :: PQA          ! and specific humidity at roof level
REAL, DIMENSION(:), INTENT(IN)    :: PPS          ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD      ! atmospheric infrared radiation
REAL,               INTENT(IN)    :: PTSTEP     ! time step
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_WALL     ! heat capacity for wall layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_WALL     ! thermal conductivity for wall layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_WALL      ! depth of wall layers
REAL, DIMENSION(:), INTENT(IN)    :: PTI_BLD      ! inside building temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTI_BLDWFR   ! inside building temperature without heating
REAL, DIMENSION(:), INTENT(IN)    :: PAC_BLD      ! aerodynamical conductance
                                                  ! inside the building itself
REAL, DIMENSION(:), INTENT(IN)    :: PTSNOW_ROAD  ! road snow temperature
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA        ! rho
REAL, DIMENSION(:), INTENT(IN)    :: PAC_WALL     ! aerodynamical conductance
!                                                 ! between wall and canyon
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_WALL ! absorbed solar radiation
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_WALL ! absorbed infrared rad.
!
REAL, DIMENSION(:), INTENT(IN)    :: PLW_W_TO_W    ! LW interactions wall  -> wall 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_R_TO_W    ! LW interactions road  -> wall 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_G_TO_W    ! LW interactions green -> wall 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_W    ! LW interactions sky   -> wall 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_NR_TO_W   ! LW interactions road(snow) -> wall 
!
REAL, DIMENSION(:), INTENT(INOUT) :: PDQS_WALL     !heat storage inside the wall 
REAL, DIMENSION(:), INTENT(INOUT) :: PQF_WALL      !flux from bld to wall due to space heating
REAL, DIMENSION(:), INTENT(INOUT) :: PFLX_BLD_WALL !flux from bld to wall
!
!*      0.2    declarations of local variables
!
!
REAL :: ZIMPL=0.5      ! implicit coefficient
REAL :: ZEXPL=0.5      ! explicit coefficient
!
REAL, DIMENSION(SIZE(PTA),SIZE(PT_WALL,2)) ::  ZA,& ! lower diag.
                                               ZB,& ! main  diag.
                                               ZC,& ! upper diag.
                                               ZY,& ! r.h.s.
                                               ZX   ! solution

!
REAL, DIMENSION(SIZE(PTA)) :: ZRHO_AC_W   ! rho * conductance (for walls)
REAL, DIMENSION(SIZE(PTA),SIZE(PT_WALL,2)) :: ZMTC_O_D_WALL
! mean thermal conductivity over distance between 2 layers
REAL, DIMENSION(SIZE(PTA),SIZE(PT_WALL,2)) :: ZHC_D_WALL
! thermal capacity times layer depth
REAL, DIMENSION(SIZE(PTA)) :: ZTS_WALL    ! wall surface temperature
REAL, DIMENSION(SIZE(PT_WALL,1)) :: ZEI_WALL  ! internal energy of walls
REAL, DIMENSION(SIZE(PT_WALL,1)) :: ZPEI_WALL ! internal energy of walls at time t+
REAL, DIMENSION(SIZE(PT_WALL,1)) :: ZTI_WALL    !temperature of internal wall layer at time t
REAL, DIMENSION(SIZE(PT_WALL,1)) :: ZFLX_BLD_WALL_EQ !heat flux between inside of
                                                   !the building and the wall without heating
!
INTEGER :: IWALL_LAYER           ! number of wall layers
INTEGER :: ILAYER                ! current layer
INTEGER :: JLAYER, JJ            ! loop counter
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WALL_LAYER_E_BUDGET',0,ZHOOK_HANDLE)
!
!*      1.     Layer thermal properties
!              ------------------------
!
!*      1.1    Walls
!              -----
!
IWALL_LAYER = SIZE(PT_WALL,2)
ZMTC_O_D_WALL(:,:) = 0.
!
DO JLAYER=1,IWALL_LAYER-1
!
  DO JJ=1,SIZE(PT_WALL,1)
!
    ZMTC_O_D_WALL(JJ,JLAYER) = 2./(  PD_WALL(JJ,JLAYER  )/PTC_WALL(JJ,JLAYER  ) &
                                + PD_WALL(JJ,JLAYER+1)/PTC_WALL(JJ,JLAYER+1) )
    ZHC_D_WALL   (JJ,JLAYER) = PHC_WALL(JJ,JLAYER) * PD_WALL (JJ,JLAYER)
  ENDDO
!
END DO
!
DO JJ=1,SIZE(PT_WALL,1)
  ZMTC_O_D_WALL(JJ,IWALL_LAYER) = 2. * PTC_WALL(JJ,IWALL_LAYER) &
                                  / PD_WALL (JJ,IWALL_LAYER)
  ZMTC_O_D_WALL(JJ,IWALL_LAYER) = 1./(  1./ZMTC_O_D_WALL(JJ,IWALL_LAYER)    &
                                   + 1./(XCPD*PRHOA(JJ)*PAC_BLD(JJ))      )
!
  ZHC_D_WALL   (JJ,IWALL_LAYER) = PHC_WALL(JJ,IWALL_LAYER) &
                             * PD_WALL (JJ,IWALL_LAYER)
ENDDO
!
!-------------------------------------------------------------------------------
!
!*      2.    Preliminaries
!             -------------
!
!*      storage current temperature of internal wall layer
ZTI_WALL(:)=PT_WALL(:,IWALL_LAYER)
!
!
!*      2.3    Surface temperatures
!              --------------------
!
ZTS_WALL(:) = PT_WALL(:,1)
!
DO JJ=1,SIZE(PT_WALL,1)
!*      2.2    flux properties
  ZRHO_AC_W (JJ) = PRHOA(JJ) * PAC_WALL(JJ)
!*      2.5    internal energy of walls/roads at the current time step
  ZEI_WALL(JJ)=PHC_WALL(JJ,1)*PD_WALL(JJ,1)*PT_WALL(JJ,1)
!
ENDDO
!
DO JLAYER=2,IWALL_LAYER
  DO JJ=1,SIZE(PT_WALL,1)
      ZEI_WALL(JJ)=ZEI_WALL(JJ)+(PHC_WALL(JJ,JLAYER)*PD_WALL(JJ,JLAYER)* &
                          PT_WALL(JJ,JLAYER))
   ENDDO
END DO                                
!
!-------------------------------------------------------------------------------
!
!*      3.    Inside wall layer coefficients
!             ------------------------------
!
!
DO JJ=1,SIZE(PT_WALL,1)
!
  ZA(JJ,1) =   0.

  ZB(JJ,1) =   ZHC_D_WALL(JJ,IWALL_LAYER) / PTSTEP                     &
     + ZIMPL * (   ZMTC_O_D_WALL(JJ,IWALL_LAYER  )                     &
                 + ZMTC_O_D_WALL(JJ,IWALL_LAYER-1)                     &
               )

  ZC(JJ,1) =                                                           &
       ZIMPL * ( - ZMTC_O_D_WALL(JJ,IWALL_LAYER-1)                     &
               )
!
  ZY(JJ,1) =   ZHC_D_WALL(JJ,IWALL_LAYER) / PTSTEP                          &
                                             * PT_WALL(JJ,IWALL_LAYER)      &
                 + ZMTC_O_D_WALL(JJ,IWALL_LAYER) * PTI_BLDWFR(JJ)           &
     + ZEXPL * ( - ZMTC_O_D_WALL(JJ,IWALL_LAYER  )                          &
                       * PT_WALL(JJ,IWALL_LAYER  )                          &
                 - ZMTC_O_D_WALL(JJ,IWALL_LAYER-1)                          &
                       * PT_WALL(JJ,IWALL_LAYER  )                          &
                 + ZMTC_O_D_WALL(JJ,IWALL_LAYER-1)                          &
                       * PT_WALL(JJ,IWALL_LAYER-1)                          &
               )
!-------------------------------------------------------------------------------
!
!*      3.    Outer wall layer coefficients
!             -----------------------------
!
!
  ZA(JJ,IWALL_LAYER) =                                                      &
       ZIMPL * ( - ZMTC_O_D_WALL(JJ,1)                                      &
               )

  ZB(JJ,IWALL_LAYER) =   ZHC_D_WALL(JJ,1)/PTSTEP                            &
     + ZIMPL * ( -          4. *ZTS_WALL(JJ)**3 * PLW_W_TO_W(JJ)            &
                 + ZRHO_AC_W(JJ) * XCPD                                     &
                 + ZMTC_O_D_WALL(JJ,1)                                      &
               )

  ZC(JJ,IWALL_LAYER) = 0.

  ZY(JJ,IWALL_LAYER) =   ZHC_D_WALL(JJ,1)/PTSTEP*ZTS_WALL(JJ)                &
                 + PABS_SW_WALL(JJ)                                          &
                 + PLW_RAD     (JJ)    * PLW_S_TO_W (JJ)                     &
                 + ZTS_WALL    (JJ)**4 * PLW_W_TO_W (JJ)                     &
                 + PTS_ROAD    (JJ)**4 * PLW_R_TO_W (JJ)                     &
                 + PTSNOW_ROAD (JJ)**4 * PLW_NR_TO_W(JJ)                     &
                 + PTS_GARDEN  (JJ)**4 * PLW_G_TO_W (JJ)                     &
                 + ZRHO_AC_W(JJ) * XCPD * PT_CANYON(JJ)                      &
     + ZIMPL * (                                                             &
                 - 4.*ZTS_WALL(JJ)**4 * PLW_W_TO_W(JJ)                       &
               )                                                             &
     + ZEXPL * (   ZMTC_O_D_WALL(JJ,1) * PT_WALL(JJ,2)                       &
                 - ZMTC_O_D_WALL(JJ,1) * PT_WALL(JJ,1)                       &
                 - ZRHO_AC_W(JJ) * XCPD *ZTS_WALL(JJ)                        &
               )

END DO
!
!-------------------------------------------------------------------------------
!
!*      5.     Other wall layers coefficients
!              ------------------------------
!
DO JLAYER=2,IWALL_LAYER-1

  ILAYER=IWALL_LAYER-JLAYER+1

  DO JJ=1,SIZE(PT_WALL,1)

    ZA(JJ,ILAYER) =                                                       &
         ZIMPL * ( - ZMTC_O_D_WALL(JJ,JLAYER  )                           &
                 )

    ZB(JJ,ILAYER) =   ZHC_D_WALL(JJ,JLAYER)/PTSTEP                        &
       + ZIMPL * (   ZMTC_O_D_WALL(JJ,JLAYER  )                           &
                   + ZMTC_O_D_WALL(JJ,JLAYER-1)                           &
                 )

    ZC(JJ,ILAYER) =                                                       &
         ZIMPL * ( - ZMTC_O_D_WALL(JJ,JLAYER-1)                           &
                 )
!
    ZY(JJ,ILAYER) =   ZHC_D_WALL(JJ,JLAYER)/PTSTEP * PT_WALL(JJ,JLAYER)  &
        + ZEXPL * (    ZMTC_O_D_WALL(JJ,JLAYER  ) * PT_WALL(JJ,JLAYER+1) &
                     - ZMTC_O_D_WALL(JJ,JLAYER  ) * PT_WALL(JJ,JLAYER  ) &
                     - ZMTC_O_D_WALL(JJ,JLAYER-1) * PT_WALL(JJ,JLAYER  ) &
                     + ZMTC_O_D_WALL(JJ,JLAYER-1) * PT_WALL(JJ,JLAYER-1) &
                  )
  ENDDO
END DO
!
!-------------------------------------------------------------------------------
!
!*      6.     Tri-diagonal system resolution
!              ------------------------------
!
CALL TRIDIAG_GROUND(ZA,ZB,ZC,ZY,ZX)
!
DO JLAYER=1,IWALL_LAYER
  ILAYER=IWALL_LAYER-JLAYER+1
  PT_WALL(:,JLAYER) = ZX(:,ILAYER)
END DO  
!
!-------------------------------------------------------------------------------
!
!*     7.  computation of flux between bld and wall without heating
!          --------------------------------------------------------
!
DO JJ=1,SIZE(PT_WALL,1)

  ! computation of flux between bld and wall without heating
  ZFLX_BLD_WALL_EQ(JJ)=ZMTC_O_D_WALL(JJ,IWALL_LAYER)*(PTI_BLDWFR(JJ) &
                                  -ZIMPL*ZTI_WALL(JJ) &
                                  -ZEXPL*PT_WALL(JJ,IWALL_LAYER) &
                                         )
!
!*    11. computation of heat fluxes with heating
!
  ZY(JJ,1) =   ZHC_D_WALL(JJ,IWALL_LAYER) / PTSTEP                          &
                                             * PT_WALL(JJ,IWALL_LAYER)      &
                 + ZMTC_O_D_WALL(JJ,IWALL_LAYER) * PTI_BLD(JJ)              &
     + ZEXPL * ( - ZMTC_O_D_WALL(JJ,IWALL_LAYER  )                          &
                       * PT_WALL(JJ,IWALL_LAYER  )                          &
                 - ZMTC_O_D_WALL(JJ,IWALL_LAYER-1)                          &
                       * PT_WALL(JJ,IWALL_LAYER  )                          &
                 + ZMTC_O_D_WALL(JJ,IWALL_LAYER-1)                          &
                       * PT_WALL(JJ,IWALL_LAYER-1)                          &
               )
ENDDO
!               
CALL TRIDIAG_GROUND(ZA,ZB,ZC,ZY,ZX)
!
DO JLAYER=1,IWALL_LAYER
  ILAYER=IWALL_LAYER-JLAYER+1
  PT_WALL(:,JLAYER) = ZX(:,ILAYER)
END DO 
!
DO JJ=1,SIZE(PT_WALL,1)
  PFLX_BLD_WALL(JJ)=ZMTC_O_D_WALL(JJ,IWALL_LAYER)*(PTI_BLD(JJ)       &
                                  -ZIMPL*ZTI_WALL(JJ)            &
                                  -ZEXPL*PT_WALL(JJ,IWALL_LAYER) &
                                         )
!
  PQF_WALL(JJ)=PFLX_BLD_WALL(JJ)-ZFLX_BLD_WALL_EQ(JJ)
!
ENDDO
!-------------------------------------------------------------------------------
!
!*     12.    Wall absorbed infra-red radiation on snow-free surfaces
!             -------------------------------------------------------
!
!* radiative surface temperature used during the enrgy balance (linearized at t+)
!
DO JJ=1,SIZE(ZTS_WALL)
!
  ZTS_WALL(JJ) = (ZTS_WALL(JJ)**4 + 4.*ZIMPL*ZTS_WALL(JJ)**3 * (PT_WALL(JJ,1) - ZTS_WALL(JJ)))**0.25
!
!* absorbed LW

  PABS_LW_WALL       (JJ) = PLW_S_TO_W (JJ)*PLW_RAD      (JJ)         &
                       + PLW_W_TO_W (JJ)*ZTS_WALL     (JJ)**4        &
                       + PLW_R_TO_W (JJ)*PTS_ROAD     (JJ)**4        &
                       + PLW_G_TO_W (JJ)*PTS_GARDEN   (JJ)**4        &
                       + PLW_NR_TO_W(JJ)*PTSNOW_ROAD  (JJ)**4
!
ENDDO
!-------------------------------------------------------------------------------
!
!
!*      13     heat storage inside walls
!              -------------------------
!
!       13.1   internal energy of the walls at the next time step
!              --------------------------------------------------
!
ZPEI_WALL(:)=PHC_WALL(:,1)*PD_WALL(:,1)*PT_WALL(:,1)
DO JLAYER=2,IWALL_LAYER
    ZPEI_WALL(:)=ZPEI_WALL(:)+(PHC_WALL(:,JLAYER)*PD_WALL(:,JLAYER)* &
                        PT_WALL(:,JLAYER))
END DO
!
!       13.2   heat storage flux inside walls 
!              ------------------------------
!
PDQS_WALL(:)=(ZPEI_WALL(:)-ZEI_WALL(:))/PTSTEP
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('WALL_LAYER_E_BUDGET',1,ZHOOK_HANDLE)
!
END SUBROUTINE WALL_LAYER_E_BUDGET
