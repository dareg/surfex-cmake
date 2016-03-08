!   ##########################################################################
    SUBROUTINE ROAD_LAYER_E_BUDGET(T, B, PTSTEP, PDN_RD, PRHOA, PAC_RD, PAC_RD_WAT, &
                                   PLW_RAD, PPS, PQSAT_RD, PDELT_RD, PEXNS,           &
                                   PABS_SW_RD, PGSNOW_RD, PQ_LOWCAN, PT_LOWCAN,       &
                                   PTS_WALL_A, PTS_WALL_B, PTSNOW_RD, PTS_GARDEN,       &
                                   PLW_WA_TO_R, PLW_WB_TO_R, PLW_S_TO_R, PLW_WIN_TO_R,    &
                                   PEMIT_LW_RD, PDQS_RD, PABS_LW_RD, PHFREE_RD,   &
                                   PLEFREE_RD, PIMB_RD, PRR )
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
!!      V. Masson           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    23/01/98 
!!                  21/11/01 (V. Masson and A. Lemonsu) bug of latent flux
!!                           for very strong evaporation (all reservoir emptied
!!                           in one time-step)
!!                     02/11 (V. Masson) split of the routine for roads and walls separately
!!      G. Pigeon      09/2012: add heating/cooling of rain from air temperature
!!                             to surface road temp. for the road energy budget 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_TEB_n, ONLY : TEB_1P_t
USE MODD_BEM_n, ONLY : BEM_1P_t
!
USE MODD_CSTS,ONLY : XCPD, XLVTT, XSTEFAN, XCL
!
USE MODE_THERMOS
!
USE MODI_LAYER_E_BUDGET
USE MODI_LAYER_E_BUDGET_GET_COEF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(TEB_1P_t), INTENT(INOUT) :: T
TYPE(BEM_1P_t), INTENT(INOUT) :: B
!
REAL,               INTENT(IN)    :: PTSTEP       ! time step
REAL, DIMENSION(:), INTENT(IN)    :: PDN_RD     ! road snow fraction
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA        ! rho
REAL, DIMENSION(:), INTENT(IN)    :: PAC_RD     ! aerodynamical conductance
!                                                 ! between road and canyon
REAL, DIMENSION(:), INTENT(IN)    :: PAC_RD_WAT ! aerodynamical conductance
!                                                 ! between road and canyon
!                                                 ! (for water)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD      ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PPS          ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PQSAT_RD   ! q_sat(Ts)
REAL, DIMENSION(:), INTENT(IN)    :: PDELT_RD   ! fraction of water
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS        ! surface Exner function
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_RD ! absorbed solar radiation
REAL, DIMENSION(:), INTENT(IN)    :: PGSNOW_RD  ! road snow conduction
!                                                 ! heat fluxes at mantel
!                                                 ! base
REAL, DIMENSION(:), INTENT(IN)    :: PQ_LOWCAN    ! and specific humidity
REAL, DIMENSION(:), INTENT(IN)    :: PT_LOWCAN    ! low canyon air temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_WALL_A   ! wall surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_WALL_B   ! wall surface temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTSNOW_RD  ! road snow temperature
REAL, DIMENSION(:), INTENT(IN)    :: PTS_GARDEN   ! green area surface temperature
!
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WA_TO_R   ! LW interactions wall  -> road 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WB_TO_R   ! LW interactions wall  -> road 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_S_TO_R    ! LW interactions sky   -> road 
REAL, DIMENSION(:), INTENT(IN)    :: PLW_WIN_TO_R ! LW interactions window -> road 
!
REAL, DIMENSION(:), INTENT(OUT)   :: PEMIT_LW_RD! LW flux emitted by the road (W/m2 of road)
REAL, DIMENSION(:), INTENT(OUT)   :: PDQS_RD    !heat storage inside the road
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_RD ! absorbed infrared rad.
REAL, DIMENSION(:), INTENT(OUT)   :: PHFREE_RD  ! sensible heat flux on the
                                                  ! snow free part of the road [W m-2]
REAL, DIMENSION(:), INTENT(OUT)   :: PLEFREE_RD ! latent heat flux on the
                                                  ! snow free part of the road [W m-2]
REAL, DIMENSION(:), INTENT(OUT)   :: PIMB_RD    ! road residual energy imbalance 
                                                  ! for verification [W m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PRR          ! rain rate [kg m-2 s-1]

!
!*      0.2    declarations of local variables
!
REAL :: ZIMPL=1.0      ! implicit coefficient
REAL :: ZEXPL=0.0      ! explicit coefficient
!
REAL, DIMENSION(SIZE(T%XT_ROAD,1),SIZE(T%XT_ROAD,2)) :: ZA,& ! lower diag.
                                                    ZB,& ! main  diag.
                                                    ZC,& ! upper diag.
                                                    ZY   ! r.h.s.
!
REAL, DIMENSION(SIZE(PPS)) :: ZDN_RD    ! snow-covered surface fraction on road
REAL, DIMENSION(SIZE(PPS)) :: ZDF_RD    ! snow-free surface fraction on road
!
REAL, DIMENSION(SIZE(PPS)) :: ZDQSAT_RD ! dq_sat/dTs
REAL, DIMENSION(SIZE(PPS)) :: ZRHO_ACF_R  ! rho * conductance
!                                         !     * snow-free f.
REAL, DIMENSION(SIZE(PPS)) :: ZRHO_ACF_R_WAT ! rho * conductance for water
!                                         !     * snow-free f.
! thermal capacity times layer depth
REAL, DIMENSION(SIZE(PPS)) :: ZTS_RD    ! road surface temperature
REAL, DIMENSION(SIZE(PPS)) :: ZHEAT_RR    ! heat used too cool/heat the rain from the roof
REAL, DIMENSION(SIZE(PPS)) :: ZT_SKY      ! road surface temperature
!
INTEGER :: IRD_LAYER           ! number of road layers
INTEGER :: JJ            ! loop counter
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ROAD_LAYER_E_BUDGET',0,ZHOOK_HANDLE)
!
 CALL LAYER_E_BUDGET_GET_COEF( T%XT_ROAD, PTSTEP, ZIMPL, T%XHC_ROAD, T%XTC_ROAD, T%XD_ROAD, &
                              ZA, ZB, ZC, ZY )
!
!*      1.     Layer thermal properties
!              ------------------------
!
IRD_LAYER = SIZE(T%XT_ROAD,2)
!
DO JJ=1, SIZE(PDN_RD) 
  !
  ZDN_RD(JJ) = PDN_RD (JJ)
  ZDF_RD(JJ) = 1. - ZDN_RD (JJ)
  !
  !*      2.3    Surface temperatures
  !              --------------------
  !
  ZTS_RD(JJ) = T%XT_ROAD(JJ,1)
  !
  !*      2.2    flux properties
  !              ---------------
  !
  ZRHO_ACF_R    (JJ) = PRHOA(JJ) * PAC_RD(JJ)     * ZDF_RD(JJ)
  ZRHO_ACF_R_WAT(JJ) = PRHOA(JJ) * PAC_RD_WAT(JJ) * ZDF_RD(JJ)
  !
  !*     2.4   Sky temperature
  !            ---------------
  !
  ZT_SKY(JJ) = (PLW_RAD(JJ)/XSTEFAN)**0.25
  !
ENDDO
!
!*      2.4    qsat, dqsat/dTs, and humidity for roads
!              ---------------------------------------
! 
ZDQSAT_RD(:) = DQSAT(ZTS_RD(:),PPS(:),PQSAT_RD(:))
!
!-------------------------------------------------------------------------------
!
!*      3.     First road layers coefficients (in contact with outdoor env.)
!              -------------------------------------------------------------
!
!print*,'ZY ',ZY(32,1)
!print*,'EXNS ',PEXNS(32)
!print*,'RHO_ACF_R ',ZRHO_ACF_R(32)
!print*,'T_LOWCAN ',PT_LOWCAN(32)
!print*,'TS_RD ',ZTS_RD(32)
!print*,'DF_RD ',ZDF_RD(32)
!print*,'ABS_SW_RD ',PABS_SW_RD(32)
!print*,'DN_RD ',ZDN_RD(32)
!print*,'GSNOW_RD ',PGSNOW_RD(32)
!print*,'RHO_ACF_R_WAT ',ZRHO_ACF_R_WAT(32)
!print*,'DELT_RD ',PDELT_RD(32)
!print*,'Q_LOWCAN ',PQ_LOWCAN(32)
!print*,'QSAT_RD ',PQSAT_RD(32)
!print*,'DQSAT_RD ',ZDQSAT_RD(32)
!print*,'LW_S_TO_R ',PLW_S_TO_R(32)
!print*,'T_SKY ',ZT_SKY(32)
!print*,'PLW_WIN_TO_R ',PLW_WIN_TO_R(32)
!print*,'T_WIN1 ',B%XT_WIN1(32)
!print*,'LW_WA_TO_R ',PLW_WA_TO_R(32)
!print*,'TS_WALL_A ',PTS_WALL_A(32)
!print*,'LW_WB_TO_R ',PLW_WB_TO_R(32)
!print*,'TS_WALL_B ',PTS_WALL_B(32)
!print*,'RR ',PRR(32)
!
DO JJ=1,SIZE(T%XT_ROAD,1)
  !
  ZB(JJ,1) = ZB(JJ,1) + ZIMPL * XCPD/PEXNS(JJ) * ZRHO_ACF_R(JJ) &
                      + ZIMPL * XLVTT * ZRHO_ACF_R_WAT(JJ) * PDELT_RD(JJ) * ZDQSAT_RD(JJ)
  !
  ZY(JJ,1) = ZY(JJ,1)  &
             + XCPD/PEXNS(JJ) * ZRHO_ACF_R(JJ) * ( PT_LOWCAN(JJ) - ZEXPL * ZTS_RD(JJ) ) &
             + ZDF_RD(JJ)*PABS_SW_RD(JJ) + ZDN_RD(JJ)*PGSNOW_RD(JJ)               &
             + XLVTT * ZRHO_ACF_R_WAT(JJ) * PDELT_RD(JJ)                                &
               * ( PQ_LOWCAN(JJ) - PQSAT_RD(JJ) + ZIMPL * ZDQSAT_RD(JJ) * ZTS_RD(JJ) )     
  !
  ZB(JJ,1) = ZB(JJ,1) &
             + ZIMPL * ZDF_RD(JJ) * ( PLW_S_TO_R(JJ) + PLW_WA_TO_R(JJ) + &
                                        PLW_WB_TO_R(JJ) + PLW_WIN_TO_R(JJ) + &
                                        PRR(JJ) * XCL ) ! heat/cool rain
  !
  ZY(JJ,1) = ZY(JJ,1) &
             + ZDF_RD(JJ) * (                                             &
               PLW_S_TO_R  (JJ) * (ZT_SKY    (JJ) - ZEXPL * ZTS_RD(JJ))   &
             + PLW_WIN_TO_R(JJ) * (B%XT_WIN1 (JJ) - ZEXPL * ZTS_RD(JJ))   &
             + PLW_WA_TO_R (JJ) * (PTS_WALL_A(JJ) - ZEXPL * ZTS_RD(JJ))   & 
             + PLW_WB_TO_R (JJ) * (PTS_WALL_B(JJ) - ZEXPL * ZTS_RD(JJ))   &
             + PRR(JJ) * XCL *    (PT_LOWCAN (JJ) - ZEXPL * ZTS_RD(JJ) ))   !heat/cool rain     
  !     
ENDDO
!
!
!print*,'A ',ZA(32,:)
!print*,'B ',ZB(32,:)
!print*,'C ',ZC(32,:)
!print*,'Y ',ZY(32,1)
 CALL LAYER_E_BUDGET( T%XT_ROAD, PTSTEP, ZIMPL, T%XHC_ROAD, T%XTC_ROAD, T%XD_ROAD, &
                     ZA, ZB, ZC, ZY, PDQS_RD )
!
!-------------------------------------------------------------------------------
!
!*     12.    Road and wall absorbed infra-red radiation on snow-free surfaces
!             ----------------------------------------------------------------
!
!* absorbed LW
DO JJ=1,SIZE(T%XT_ROAD,1)
  !
  ! surface temperature used in energy balance
  ZTS_RD(JJ) = ZEXPL *  ZTS_RD(JJ) + ZIMPL * T%XT_ROAD(JJ,1)
  PABS_LW_RD(JJ) = PLW_S_TO_R  (JJ) * (ZT_SKY(JJ)     - ZTS_RD(JJ)) + &
                     PLW_WA_TO_R (JJ) * (PTS_WALL_A(JJ) - ZTS_RD(JJ)) + &
                     PLW_WB_TO_R (JJ) * (PTS_WALL_B(JJ) - ZTS_RD(JJ)) + &
                     PLW_WIN_TO_R(JJ) * (B%XT_WIN1(JJ)  - ZTS_RD(JJ))
  !
  !*     9.    Road emitted LW radiation on snow-free surfaces
  !            -----------------------------------------------
  PEMIT_LW_RD(JJ) = XSTEFAN * T%XT_ROAD(JJ,1)**4 + &
                      (1 - T%XEMIS_ROAD(JJ))/T%XEMIS_ROAD(JJ) * PABS_LW_RD(JJ)
  !
  !*      10.     road and wall sensible heat flux
  !              --------------------------------
  !
  PHFREE_RD(JJ) = ZRHO_ACF_R(JJ) * XCPD/PEXNS(JJ) * &
                   ( ZIMPL*T%XT_ROAD(JJ,1) + ZEXPL*ZTS_RD(JJ) - PT_LOWCAN(JJ) )
  !
  !*      11     road latent heat flux
  !              ---------------------
  !
  PLEFREE_RD(JJ) = ZRHO_ACF_R_WAT(JJ) * XLVTT * PDELT_RD(JJ) * &
                    ( PQSAT_RD(JJ) - PQ_LOWCAN(JJ) +             &
                     ZIMPL * ZDQSAT_RD(JJ) * (T%XT_ROAD(JJ,1) - ZTS_RD(JJ)) )
  ZHEAT_RR(JJ) = PRR(JJ) * XCL * (ZTS_RD(JJ) - PT_LOWCAN(JJ))
  !
  !*      12     heat storage inside roads
  !              -------------------------
  !
  !*      13     road energy residual imbalance for verification
  !              -----------------------------------------------
  !
  PIMB_RD(JJ) = PABS_SW_RD(JJ) + PABS_LW_RD(JJ) - PDQS_RD(JJ) &
               - ZDF_RD(JJ) * ( PHFREE_RD(JJ) + PLEFREE_RD(JJ)) &
               - ZDN_RD(JJ) *   PGSNOW_RD(JJ)
  !
ENDDO
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ROAD_LAYER_E_BUDGET',1,ZHOOK_HANDLE)
!
END SUBROUTINE ROAD_LAYER_E_BUDGET

