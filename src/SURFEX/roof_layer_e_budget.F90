!     #########
    SUBROUTINE ROOF_LAYER_E_BUDGET(TOP, T, B, PQSAT_RF, PAC_BLD, PTSTEP, PDN_RF,   &
                                   PRHOA, PAC_RF, PAC_RF_WAT, PLW_RAD, PPS,        &
                                   PDELT_RF, PTA, PQA, PEXNA, PEXNS, PABS_SW_RF,   &
                                   PGSNOW_RF, PFLX_BLD_RF, PDQS_RF, PABS_LW_RF,    &
                                   PHFREE_RF, PLEFREE_RF, PIMB_RF,                 &
                                   PG_GRF_RF, PRADHT_IN, PTS_FLOOR, PTI_WL,    &
                                   PRAD_RF_WL, PRAD_RF_WIN, PRAD_RF_FLOOR,         &
                                   PRAD_RF_MASS, PCONV_RF_BLD, PRR, PLOAD_IN_RF )
!   ##########################################################################
!
!!****  *ROOF_LAYER_E_BUDGET*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the evoultion of surface temperature of roofs
!         
!     
!!**  METHOD
!     ------
!
!
!
!
!    5 : equation for evolution of Ts_roof
!        *********************************
!
!     dTt_1(t) / dt = 1/(dt_1*Ct_1) * (  Rn - H - LE
!                                      - 2*Kt_1*(Tt_1-Tt_2)/(dt_1 +dt_2)       )
!
!     dTt_k(t) / dt = 1/(dt_k*Ct_k) * (- 2*Kt_k-1*(Tt_k-Tt_k-1)/(dt_k-1 +dt_k) 
!                                      - 2*Kt_k  *(Tt_k-Tt_k+1)/(dt_k+1 +dt_k) )
!
!       with
!
!       K*_k  = (d*_k+ d*_k+1)/(d*_k/k*_k+ d*_k+1/k*_k+1)
!
!       Rn = (dir_Rg + sca_Rg) (1-a) + emis * ( Rat - sigma Ts**4 (t+dt) )
!
!       H  = rho Cp CH V ( Ts (t+dt) - Tas )
!
!       LE = rho Lv CH V ( qs (t+dt) - qas )
!
!      where the as subscript denotes atmospheric values at ground level
!      (and not at first half level)
!
!      The tridiagonal systel is linearized with
!
!       using      Ts**4 (t+dt) = Ts**4 (t) + 4*Ts**3 (t) * ( Ts(t+dt) - Ts(t) )
!
!       and  qs (t+dt) = Hu(t) * qsat(t) + Hu(t) dqsat/dT * ( Ts(t+dt) - Ts(t) )
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
!!                  17/10/05 (G. Pigeon) computation of storage inside the roofs
!!                  26/04/12 (G. Pigeon) add term of heating of rain (new arg PRR+XCL)
!!                     09/12 (G. Pigeon) modif of indoor conv. coef and implicit calculation
!!                     10/12 (G. Pigeon) add indoor solar heat load
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_n, ONLY : TEB_1P_t
USE MODD_BEM_n, ONLY : BEM_1P_t
!
USE MODD_SURF_PAR,  ONLY : XUNDEF
USE MODD_CSTS,ONLY : XCPD, XLVTT, XSTEFAN, XCL
!
USE MODE_THERMOS
!
USE MODI_LAYER_E_BUDGET
USE MODI_LAYER_E_BUDGET_GET_COEF
USE MODE_CONV_DOE
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_1P_t), INTENT(INOUT) :: T
TYPE(BEM_1P_t), INTENT(INOUT) :: B
!
REAL, DIMENSION(:), INTENT(INOUT) :: PQSAT_RF     ! q_sat(Ts)
REAL, DIMENSION(:), INTENT(IN)    :: PAC_BLD        ! aerodynamical resistance
                                                    ! inside building itself
REAL,               INTENT(IN)    :: PTSTEP         ! time step
REAL, DIMENSION(:), INTENT(IN)    :: PDN_RF       ! roof snow fraction
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA          ! air density
REAL, DIMENSION(:), INTENT(IN)    :: PAC_RF       ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(IN)    :: PAC_RF_WAT   ! aerodynamical conductance (for water)
REAL, DIMENSION(:), INTENT(IN)    :: PLW_RAD        ! atmospheric infrared radiation
REAL, DIMENSION(:), INTENT(IN)    :: PPS            ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PDELT_RF     ! fraction of water
REAL, DIMENSION(:), INTENT(IN)    :: PTA            ! air temperature at roof level
REAL, DIMENSION(:), INTENT(IN)    :: PQA            ! air specific humidity
                                                    ! at roof level
REAL, DIMENSION(:), INTENT(IN)    :: PEXNA          ! exner function
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS          ! surface exner function
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_RF   ! absorbed solar radiation
REAL, DIMENSION(:), INTENT(IN)    :: PGSNOW_RF    ! roof snow conduction
!                                                   ! heat fluxes at mantel
!                                                   ! base
REAL, DIMENSION(:), INTENT(IN)    :: PG_GRF_RF ! heat conduction flux
!                                                        between greenroof and
!                                                        structural roof
REAL, DIMENSION(:), INTENT(OUT)   :: PFLX_BLD_RF  ! flux from bld to roof
REAL, DIMENSION(:), INTENT(OUT)   :: PDQS_RF      ! heat storage inside the roofs
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_RF   ! absorbed infra-red rad.
REAL, DIMENSION(:), INTENT(OUT)   :: PHFREE_RF    ! sensible heat flux of the
                                                    ! snow free part of the roof
REAL, DIMENSION(:), INTENT(OUT)   :: PLEFREE_RF   ! latent heat flux of the
                                                    ! snow free part of the roof
REAL, DIMENSION(:), INTENT(OUT)   :: PIMB_RF      ! residual energy imbalance
                                                    ! of the roof for
REAL, DIMENSION(:), INTENT(IN)    :: PRADHT_IN      ! Indoor radiant heat transfer coefficient
                                                    ! [W K-1 m-2]
REAL, DIMENSION(:), INTENT(IN)    :: PTS_FLOOR      ! surf. floor temp. (contact with bld air)
REAL, DIMENSION(:), INTENT(IN)    :: PTI_WL       ! indoor wall temp.
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_RF_WL ! rad. fluxes from roof to wall [W m-2(roof)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_RF_WIN  ! rad. fluxes from roof to win [W m-2(roof)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_RF_FLOOR! rad. fluxes from roof to floor [W m-2(roof)]
REAL, DIMENSION(:), INTENT(OUT)   :: PRAD_RF_MASS ! rad. fluxes from roof to mass [W m-2(roof)]
REAL, DIMENSION(:), INTENT(OUT)   :: PCONV_RF_BLD ! conv. fluxes from roof to bld [W m-2(roof)]
REAL, DIMENSION(:), INTENT(IN)    :: PRR ! rain rate [kg m-2 s-1]
REAL, DIMENSION(:), INTENT(IN)    :: PLOAD_IN_RF ! solar + int heat gain on roof W/m2 [roof]
!
!*      0.2    declarations of local variables
!
REAL :: ZIMPL = 1.0        ! implicit coefficient
REAL :: ZEXPL = 0.0        ! explicit coefficient
!
REAL, DIMENSION(SIZE(PTA)) :: ZDF_RF ! snow-free fraction
REAL, DIMENSION(SIZE(PTA),SIZE(T%XT_ROOF,2)) :: ZA,& ! lower diag.
                                              ZB,& ! main  diag.
                                              ZC,& ! upper diag.
                                              ZY   ! r.h.s.
!
REAL, DIMENSION(SIZE(PTA)) :: ZDQSAT_RF      ! dq_sat/dTs
REAL, DIMENSION(SIZE(PTA)) :: ZRHO_ACF_RF    ! conductance * rho
REAL, DIMENSION(SIZE(PTA)) :: ZRHO_ACF_RF_WAT! conductance * rho (for water)
REAL, DIMENSION(SIZE(PTA)) :: ZMTC_O_D_RF_IN ! thermal capacity times layer depth
REAL, DIMENSION(SIZE(PTA)) :: ZTS_RF         ! roof surface temperature at previous time step
REAL, DIMENSION(SIZE(PTA)) :: ZTRAD_RF       ! roof radiative surface temperature at intermediate time step
REAL, DIMENSION(SIZE(PTA)) :: ZTAER_RF       ! roof aerodyn. surface temperature at intermediate time step
REAL, DIMENSION(SIZE(PTA)) :: ZHEAT_RR         ! heat used too cool/heat the rain from the roof
REAL, DIMENSION(SIZE(PTA)) :: ZTI_RF         ! temperature of internal roof layer used for radiative exchanges
REAL, DIMENSION(SIZE(PTA)) :: ZTI_RF_CONV    ! temperature of internal roof layer used for convective exchanges
REAL, DIMENSION(SIZE(PTA)) :: ZCHTC_IN_RF      ! Indoor roof convec heat transfer coefficient
                                               ! [W K-1 m-2(bld)]
!
INTEGER :: JJ
INTEGER :: IRF_LAYER           ! number of roof layers
INTEGER :: JLAYER                ! loop counter
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ROOF_LAYER_E_BUDGET',0,ZHOOK_HANDLE)
!
PRAD_RF_WL(:) = XUNDEF
PRAD_RF_WIN(:)  = XUNDEF
PRAD_RF_FLOOR(:)= XUNDEF
PRAD_RF_MASS(:) = XUNDEF
PCONV_RF_BLD(:) = XUNDEF
!
! *Convection heat transfer coefficients [W m-2 K-1] from EP Engineering Reference
!
IRF_LAYER = SIZE(T%XT_ROOF,2)
!
ZCHTC_IN_RF(:) = CHTC_DOWN_DOE(T%XT_ROOF(:,IRF_LAYER), B%XTI_BLD(:))
DO JJ=1,SIZE(ZCHTC_IN_RF)
   ZCHTC_IN_RF(JJ) = MAX(1., ZCHTC_IN_RF(JJ))
ENDDO
!
 CALL LAYER_E_BUDGET_GET_COEF(T%XT_ROOF, PTSTEP, ZIMPL, T%XHC_ROOF, T%XTC_ROOF, T%XD_ROOF, &
                              ZA, ZB, ZC, ZY )
!
!
DO JJ=1,SIZE(PDN_RF)
  !
  ZDF_RF(JJ) = 1. - PDN_RF(JJ)
  !
  ZTS_RF(JJ) = T%XT_ROOF(JJ,1)
  ZTI_RF(JJ) = T%XT_ROOF(JJ,IRF_LAYER)
  !
  !*      2.     Roof Ts coefficients
  !              --------------------
  !
  ZRHO_ACF_RF    (JJ) = PRHOA(JJ) * PAC_RF    (JJ)
  ZRHO_ACF_RF_WAT(JJ) = PRHOA(JJ) * PAC_RF_WAT(JJ)
  !
  IF (TOP%CBEM .EQ. 'DEF') THEN
    ZMTC_O_D_RF_IN(JJ) = 2. * T%XTC_ROOF(JJ,IRF_LAYER) / T%XD_ROOF (JJ,IRF_LAYER)
    ZMTC_O_D_RF_IN(JJ) = 1./(  1./ZMTC_O_D_RF_IN(JJ) + 1./(XCPD*PRHOA(JJ)*PAC_BLD(JJ)) ) 
  ENDIF
  !
ENDDO
!
!*      2.1    dqsat/dTs, and humidity for roofs
!              ---------------------------------
!
ZDQSAT_RF(:) = DQSAT(ZTS_RF(:),PPS(:),PQSAT_RF(:))
!
!*      2.2    coefficients
!              ------------
! 
DO JJ=1,SIZE(T%XT_ROOF,1)
  !
  ZB(JJ,1) = ZB(JJ,1) + ZDF_RF(JJ) * (1.-T%XGREENROOF(JJ)) * (                           &
                        ZIMPL * ( XCPD/PEXNS(JJ) * ZRHO_ACF_RF(JJ)                       &
                        + XLVTT * ZRHO_ACF_RF_WAT(JJ) * PDELT_RF(JJ) * ZDQSAT_RF(JJ)     &
                        + XSTEFAN * T%XEMIS_ROOF(JJ) * 4.*ZTS_RF(JJ)**3                  &
                        + PRR(JJ) * XCL )) !! heating/cooling of rain 
  !
  ZY(JJ,1) = ZY(JJ,1) + (1.-T%XGREENROOF(JJ))                                                        &
                      * (PDN_RF(JJ)*PGSNOW_RF(JJ) + ZDF_RF(JJ) * ( PABS_SW_RF(JJ)                    &
                        + XCPD * ZRHO_ACF_RF(JJ) * ( PTA(JJ)/PEXNA(JJ) - ZEXPL*ZTS_RF(JJ)/PEXNS(JJ)) &
                        + T%XEMIS_ROOF(JJ)*PLW_RAD(JJ)                                               &                 
                        + XLVTT * ZRHO_ACF_RF_WAT(JJ) * PDELT_RF(JJ)                                 &
                          * ( PQA(JJ) - PQSAT_RF(JJ) + ZIMPL * ZDQSAT_RF(JJ) * ZTS_RF(JJ) )          &
                        + XSTEFAN * T%XEMIS_ROOF(JJ) * ZTS_RF(JJ)**4 * ( 3.*ZIMPL-ZEXPL )            &
                        + PRR(JJ) * XCL * (PTA(JJ) - ZEXPL * ZTS_RF(JJ)) ) ) & !! heating/cooling of rain
                        + T%XGREENROOF(JJ)*PG_GRF_RF(JJ)
  !
  IF (TOP%CBEM=="DEF") THEN
    !
    ZB(JJ,IRF_LAYER) = ZB(JJ,IRF_LAYER) + ZIMPL * ZMTC_O_D_RF_IN(JJ)
    !
    ZY(JJ,IRF_LAYER) = ZY(JJ,IRF_LAYER) &
                         + ZMTC_O_D_RF_IN(JJ) * B%XTI_BLD(JJ) &
                         - ZEXPL * ZMTC_O_D_RF_IN(JJ) * T%XT_ROOF(JJ,IRF_LAYER)
    !
  ELSEIF (TOP%CBEM=="BEM") THEN
    !
    ZB(JJ, IRF_LAYER) = ZB(JJ,IRF_LAYER) + ZIMPL * &
                         (ZCHTC_IN_RF(JJ) * 4./3. + PRADHT_IN(JJ) * &
                         (B%XF_FLOOR_MASS(JJ) + B%XF_FLOOR_WIN(JJ) + &
                          B%XF_FLOOR_WALL(JJ) + B%XF_FLOOR_ROOF(JJ) ))

    ZY(JJ,IRF_LAYER) = ZY(JJ,IRF_LAYER) + &
       ZCHTC_IN_RF(JJ) * (B%XTI_BLD(JJ) - 1./3. *  T%XT_ROOF(JJ, IRF_LAYER)*(4*ZEXPL - 1.)) + &
       PRADHT_IN(JJ) * ( &
          B%XF_FLOOR_MASS (JJ) * (B%XT_MASS(JJ,1) - ZEXPL * T%XT_ROOF(JJ,IRF_LAYER)) + &
          B%XF_FLOOR_WIN  (JJ) * (B%XT_WIN2  (JJ) - ZEXPL * T%XT_ROOF(JJ,IRF_LAYER)) + &
          B%XF_FLOOR_WALL (JJ) * (PTI_WL     (JJ) - ZEXPL * T%XT_ROOF(JJ,IRF_LAYER)) + &
          B%XF_FLOOR_ROOF (JJ) * (PTS_FLOOR  (JJ) - ZEXPL * T%XT_ROOF(JJ,IRF_LAYER))) + &
          PLOAD_IN_RF(JJ)
    !
  ENDIF
  !
ENDDO
!print*,'ZY ',ZY(1,IRF_LAYER)
!print*,'CHTC_IN_RF ',ZCHTC_IN_RF(1)
!print*,'TI_BLD ',B%XTI_BLD(1)
!print*,'RADHT_IN ',PRADHT_IN(1)
!print*,'F_FLOOR_MASS ',B%XF_FLOOR_MASS(1)
!print*,'T_MASS ',B%XT_MASS(1,1)
!print*,'T_ROOF ',T%XT_ROOF(1,IRF_LAYER)
!print*,'F_FLOOR_WIN ',B%XF_FLOOR_WIN(1)
!print*,'T_WIN2 ',B%XT_WIN2(1)
!print*,'F_FLOOR_WALL ',B%XF_FLOOR_WALL(1)
!print*,'F_FLOOR_ROOF ',B%XF_FLOOR_ROOF(1)
!print*,'TI_WL ',PTI_WL(1)
!print*,'TS_FLOOR ',PTS_FLOOR(1)
!print*,'LOAD_IN_RF ',PLOAD_IN_RF(1)
!
!print*,'A ',ZA(1,:)
!print*,'B ',ZB(1,:)
!print*,'C ',ZC(1,:)
!print*,'Y ',ZY(1,:)
 CALL LAYER_E_BUDGET( T%XT_ROOF, PTSTEP, ZIMPL, T%XHC_ROOF, T%XTC_ROOF, T%XD_ROOF, &
                     ZA, ZB, ZC, ZY, PDQS_RF )
!
!-------------------------------------------------------------------------------
!
!*     diagnostic: computation of flux between bld and internal roof layer
DO JJ=1,SIZE(T%XT_ROOF,1)
  !
  ZTI_RF_CONV(JJ) = 4./3. * ZIMPL * T%XT_ROOF(JJ, IRF_LAYER) + 1./3. * ZTI_RF(JJ) * (4*ZEXPL -1.)
  ZTI_RF(JJ) = ZEXPL * ZTI_RF(JJ) + ZIMPL * T%XT_ROOF(JJ, IRF_LAYER) 
  SELECT CASE(TOP%CBEM)
  CASE("DEF")
     PFLX_BLD_RF(JJ) = ZMTC_O_D_RF_IN(JJ) * (B%XTI_BLD(JJ) - ZTI_RF(JJ))
  CASE("BEM")
     PRAD_RF_WL(JJ)   = PRADHT_IN  (JJ) * (ZTI_RF(JJ) - PTI_WL(JJ))
     PRAD_RF_WIN(JJ)  = PRADHT_IN  (JJ) * (ZTI_RF(JJ) - B%XT_WIN2(JJ))
     PRAD_RF_FLOOR(JJ)= PRADHT_IN  (JJ) * (ZTI_RF(JJ) - PTS_FLOOR(JJ))
     PRAD_RF_MASS(JJ) = PRADHT_IN  (JJ) * (ZTI_RF(JJ) - B%XT_MASS(JJ,1))
     PCONV_RF_BLD(JJ) = ZCHTC_IN_RF(JJ) * (ZTI_RF_CONV(JJ) - B%XTI_BLD (JJ))
     PFLX_BLD_RF(JJ)  = -(PRAD_RF_WL(JJ) + PRAD_RF_WIN(JJ) + PRAD_RF_FLOOR(JJ) + &
                            PRAD_RF_MASS(JJ) + PCONV_RF_BLD(JJ))
  ENDSELECT
  
  !
  !*      8.     Infra-red radiation absorbed by roofs
  !              -------------------------------------
  !
  !print*,JJ,ZTS_RF(JJ),T%XT_ROOF(JJ,1)
  !* radiative surface temperature at intermediate time step
  ZTRAD_RF(JJ) = ( ZTS_RF(JJ)**4 + 4.*ZIMPL*ZTS_RF(JJ)**3 * (T%XT_ROOF(JJ,1) - ZTS_RF(JJ)) )**0.25
  !
  !* absorbed LW
  PABS_LW_RF(JJ) = T%XEMIS_ROOF(JJ) * (PLW_RAD(JJ) - XSTEFAN * ZTRAD_RF(JJ)** 4)
  !
  !*      9.     Sensible heat flux between snow free roof and air
  !              -------------------------------------------------
  !
  !* aerodynamic surface temperature at the intermediate time step
  ZTAER_RF(JJ) = ZEXPL * ZTS_RF(JJ) + ZIMPL * T%XT_ROOF(JJ,1)
  PHFREE_RF(JJ) = ZRHO_ACF_RF(JJ) * XCPD * ( ZTAER_RF(JJ)/PEXNS(JJ) - PTA(JJ)/PEXNA(JJ) )
  !
  ZHEAT_RR(JJ) = PRR(JJ) * XCL * (ZTAER_RF(JJ) - PTA(JJ))
  !
  !*      10.     Latent heat flux between snow free roof and air
  !              -------------------------------------------------
  !
  PLEFREE_RF(JJ) = ZRHO_ACF_RF_WAT(JJ) * XLVTT * PDELT_RF(JJ) * &
                   ( PQSAT_RF(JJ) - PQA(JJ) +                     &
                       ZIMPL * ZDQSAT_RF(JJ) * (T%XT_ROOF(JJ,1) - ZTS_RF(JJ)) ) 
  !
  !      13.     Energy imbalance for verification
  !              ---------------------------------
  PIMB_RF(JJ) = PABS_SW_RF(JJ) + PABS_LW_RF(JJ) - PDQS_RF(JJ) &
               - ZDF_RF(JJ) * ( PHFREE_RF(JJ) + PLEFREE_RF(JJ)) &
               - PDN_RF(JJ) * PGSNOW_RF(JJ) + PFLX_BLD_RF(JJ)
  !
ENDDO
!
!*      11.     New saturated specified humidity near the roof surface
!              ------------------------------------------------------
!
PQSAT_RF(:) =  QSAT(T%XT_ROOF(:,1),PPS(:))
!
!-------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ROOF_LAYER_E_BUDGET',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------
!
END SUBROUTINE ROOF_LAYER_E_BUDGET
