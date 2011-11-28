!     #########
    SUBROUTINE ROOF_LAYER_E_BUDGET(PT_ROOF, PQSAT_ROOF,                        &
                               PTA, PQA, PPS,                                    &
                               PLW, PTSTEP,                                      &
                               PEMIS_ROOF, PHC_ROOF, PTC_ROOF, PD_ROOF,          &
                               PTI_BLD, PTI_BLDWFR, PAC_BLD, PDELT_ROOF,         &
                               PDELTSNOW_ROOF, PGSNOW_ROOF,                      &
                               PRHOA, PAC_ROOF, PAC_ROOF_WAT,                    &
                               PABS_SW_ROOF, PABS_LW_ROOF,                       &
                               PDQS_ROOF, PQF_ROOF, PFLX_BLD_ROOF                )  
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
!!	V. Masson           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    23/01/98 
!!                  17/10/05 (G. Pigeon) computation of storage inside the roofs
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,ONLY : XCPD, XLVTT, XSTEFAN
!
USE MODE_THERMOS
!
USE MODI_TRIDIAG_GROUND
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
REAL, DIMENSION(:,:), INTENT(INOUT) :: PT_ROOF      ! roof layers temperatures
REAL, DIMENSION(:), INTENT(INOUT) :: PQSAT_ROOF     ! q_sat(Ts)
REAL, DIMENSION(:), INTENT(IN)    :: PTA            ! air temperature at roof level
REAL, DIMENSION(:), INTENT(IN)    :: PQA            ! air specific humidity
                                                    ! at roof level
REAL, DIMENSION(:), INTENT(IN)    :: PPS            ! pressure at the surface
REAL, DIMENSION(:), INTENT(IN)    :: PLW            ! atmospheric infrared radiation
REAL,                 INTENT(IN)    :: PTSTEP       ! time step
REAL, DIMENSION(:), INTENT(IN)    :: PEMIS_ROOF     ! roof emissivity
REAL, DIMENSION(:,:), INTENT(IN)  :: PHC_ROOF       ! heat capacity for roof layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PTC_ROOF       ! thermal conductivity for roof layers
REAL, DIMENSION(:,:), INTENT(IN)  :: PD_ROOF        ! depth of roof layers
REAL, DIMENSION(:), INTENT(IN)    :: PTI_BLD        ! inside building temp.
REAL, DIMENSION(:), INTENT(IN)    :: PTI_BLDWFR     ! inside building temp. computed
                                                    ! with its own evolution equation
REAL, DIMENSION(:), INTENT(IN)    :: PAC_BLD        ! aerodynamical resistance
                                                    ! inside building itself
REAL, DIMENSION(:), INTENT(IN)    :: PDELT_ROOF     ! fraction of water
REAL, DIMENSION(:), INTENT(IN)    :: PDELTSNOW_ROOF ! roof snow fraction
REAL, DIMENSION(:), INTENT(IN)    :: PGSNOW_ROOF    ! roof snow conduction
!                                                   ! heat fluxes at mantel
!                                                   ! base
REAL, DIMENSION(:), INTENT(IN)    :: PRHOA          ! air density
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROOF       ! aerodynamical conductance
REAL, DIMENSION(:), INTENT(IN)    :: PAC_ROOF_WAT   ! aerodynamical conductance (for water)
REAL, DIMENSION(:), INTENT(IN)    :: PABS_SW_ROOF   ! absorbed solar radiation
REAL, DIMENSION(:), INTENT(OUT)   :: PABS_LW_ROOF   ! absorbed infra-red rad.
REAL, DIMENSION(:), INTENT(OUT)   :: PDQS_ROOF      ! heat storage inside the roofs
REAL, DIMENSION(:), INTENT(OUT)   :: PQF_ROOF       ! flux from bld to roof due to space heating
REAL, DIMENSION(:), INTENT(OUT)   :: PFLX_BLD_ROOF  ! flux from bld to roof
!
!*      0.2    declarations of local variables
!
REAL :: ZIMPL = 0.5        ! implicit coefficient
REAL :: ZEXPL = 0.5        ! explicit coefficient
!
REAL, DIMENSION(SIZE(PTA)) :: ZDELTFREE_ROOF ! snow-free fraction
REAL, DIMENSION(SIZE(PTA),SIZE(PT_ROOF,2)) :: ZA ! lower diag.
REAL, DIMENSION(SIZE(PTA),SIZE(PT_ROOF,2)) :: ZB ! main  diag.
REAL, DIMENSION(SIZE(PTA),SIZE(PT_ROOF,2)) :: ZC ! upper diag.
REAL, DIMENSION(SIZE(PTA),SIZE(PT_ROOF,2)) :: ZY ! r.h.s.
REAL, DIMENSION(SIZE(PTA),SIZE(PT_ROOF,2)) :: ZX ! solution
!
REAL, DIMENSION(SIZE(PTA)) :: ZRHO_AC_ROOF ! conductance * rho
REAL, DIMENSION(SIZE(PTA)) :: ZRHO_AC_ROOF_WAT ! ! conductance * rho (for water)
REAL, DIMENSION(SIZE(PTA)) :: ZDQSAT_ROOF  ! dq_sat/dTs
REAL, DIMENSION(SIZE(PTA),SIZE(PT_ROOF,2)) :: ZMTC_O_D
! mean thermal conductivity over distance between 2 layers
REAL, DIMENSION(SIZE(PTA),SIZE(PT_ROOF,2)) :: ZHC_D_ROOF
! thermal capacity times layer depth
REAL, DIMENSION(SIZE(PTA)) :: ZTS_ROOF         ! roof surface temperature
REAL, DIMENSION(SIZE(PT_ROOF,1)) :: ZEI_ROOF   ! internal energy of roofs
REAL, DIMENSION(SIZE(PT_ROOF,1)) :: ZPEI_ROOF  ! internal energy of roofs at t+
REAL, DIMENSION(SIZE(PTA)) :: ZTI_ROOF         ! temperature of internal roof layer at t
REAL, DIMENSION(SIZE(PTA)) :: ZFLX_BLD_ROOF_EQ ! heat flux between inside of the
                                               ! building and the roof with PTI_BLDWFR
!
INTEGER :: IROOF_LAYER           ! number of roof layers
INTEGER :: JLAYER                ! loop counter
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ROOF_LAYER_E_BUDGET',0,ZHOOK_HANDLE)
ZDELTFREE_ROOF=0.
!
ZDELTFREE_ROOF=1.-PDELTSNOW_ROOF
!
!*      1.     Layer thermal properties
!              ------------------------
!
IROOF_LAYER = SIZE(PT_ROOF,2)
ZMTC_O_D(:,:) = 0.
ZHC_D_ROOF(:,:) = 0.
!
DO JLAYER=1,IROOF_LAYER-1
  ZMTC_O_D(:,JLAYER) = 2./(   PD_ROOF(:,JLAYER  )/PTC_ROOF(:,JLAYER  ) &
                              + PD_ROOF(:,JLAYER+1)/PTC_ROOF(:,JLAYER+1) )  
  ZHC_D_ROOF(:,JLAYER) = PHC_ROOF(:,JLAYER) * PD_ROOF(:,JLAYER)
END DO
!
ZMTC_O_D(:,IROOF_LAYER) = 2. * PTC_ROOF(:,IROOF_LAYER) &
                               / PD_ROOF (:,IROOF_LAYER)  
ZMTC_O_D(:,IROOF_LAYER) = 1./(  1./ZMTC_O_D(:,IROOF_LAYER)      &
                                + 1./(XCPD*PRHOA(:)*PAC_BLD(:))   )  
!
ZHC_D_ROOF(:,IROOF_LAYER) = PHC_ROOF(:,IROOF_LAYER) &
                            * PD_ROOF (:,IROOF_LAYER)  
!
!-------------------------------------------------------------------------------
!
!*       Preliminaries
!*       -------------
!
!*       computation of internal energy at the current time step
ZEI_ROOF(:)=PHC_ROOF(:,1)*PD_ROOF(:,1)*PT_ROOF(:,1)
DO JLAYER=2,IROOF_LAYER
        ZEI_ROOF(:)=ZEI_ROOF(:)+(PHC_ROOF(:,JLAYER)*PD_ROOF(:,JLAYER)* &
                          PT_ROOF(:,JLAYER))  
END DO
!*       storage current temperature of internal roof layer (modif G. Pigeon)
ZTI_ROOF(:)=PT_ROOF(:,IROOF_LAYER)
!
!*      2.     Roof Ts coefficients
!              --------------------
!
!
ZRHO_AC_ROOF    (:) = PRHOA(:) * PAC_ROOF    (:)
ZRHO_AC_ROOF_WAT(:) = PRHOA(:) * PAC_ROOF_WAT(:)
!
ZTS_ROOF(:) = PT_ROOF(:,1)
!
!*      2.1    dqsat/dTs, and humidity for roofs
!              ---------------------------------
!
ZDQSAT_ROOF(:) = DQSAT(ZTS_ROOF(:),PPS(:),PQSAT_ROOF(:))
!
!
!*      2.2    coefficients
!              ------------
!
ZA(:,1) =   0.

ZB(:,1) =   ZHC_D_ROOF(:,1)/PTSTEP                                           &
    + ZIMPL * ( ZDELTFREE_ROOF(:) * (                                          &
                 4. * PEMIS_ROOF(:) * XSTEFAN * ZTS_ROOF(:)**3                 &
               + ZRHO_AC_ROOF(:) * XCPD                                        &
               + ZRHO_AC_ROOF_WAT(:) * XLVTT * PDELT_ROOF(:) * ZDQSAT_ROOF(:)  &
                                    )                                          &
               + ZMTC_O_D(:,1)                                                 &
              )  

ZC(:,1) =                                                                    &
      ZIMPL * ( - ZMTC_O_D(:,1)                                                &
              )  
!
ZY(:,1) =   ZHC_D_ROOF(:,1)/PTSTEP * PT_ROOF(:,1)                            &
            + ZDELTFREE_ROOF(:) * (                                            &
            + PABS_SW_ROOF(:) + PEMIS_ROOF(:)*PLW(:)                           &
            + ZRHO_AC_ROOF(:) * XCPD * PTA(:)                                  &
            - ZRHO_AC_ROOF_WAT(:) * XLVTT * PDELT_ROOF(:) * PQSAT_ROOF(:)      &
            + ZRHO_AC_ROOF_WAT(:) * XLVTT * PDELT_ROOF(:) * PQA(:)             &
                                  )                                            &
            + PDELTSNOW_ROOF(:) * PGSNOW_ROOF(:)                               &
    + ZIMPL * (   ZDELTFREE_ROOF(:) * (                                        &
                    3. * PEMIS_ROOF(:) * XSTEFAN * ZTS_ROOF(:)** 4             &
                  + ZRHO_AC_ROOF_WAT(:) * XLVTT * PDELT_ROOF(:)                &
                                    * ZDQSAT_ROOF(:) * ZTS_ROOF(:)             &
                                      )                                        &
              )                                                                &
    + ZEXPL * (   ZDELTFREE_ROOF(:) * (                                        &
                - PEMIS_ROOF(:) * XSTEFAN * ZTS_ROOF(:)** 4                    &
                - ZRHO_AC_ROOF(:) * XCPD * ZTS_ROOF(:)                         &
                                      )                                        &
                - ZMTC_O_D(:,1) * PT_ROOF(:,1)                                 &
                + ZMTC_O_D(:,1) * PT_ROOF(:,2)                                 &
              )  
!
!-------------------------------------------------------------------------------
!
!*      3.     Other layers coefficients
!              -------------------------
!
DO JLAYER=2,IROOF_LAYER-1
  ZA(:,JLAYER) =                                                       &
           ZIMPL * ( - ZMTC_O_D(:,JLAYER-1)                              &
                   )  

  ZB(:,JLAYER) =   ZHC_D_ROOF(:,JLAYER)/PTSTEP                         &
         + ZIMPL * (   ZMTC_O_D(:,JLAYER-1)                              &
                     + ZMTC_O_D(:,JLAYER  )                              &
                   )  

  ZC(:,JLAYER) =                                                       &
           ZIMPL * ( - ZMTC_O_D(:,JLAYER  )                              &
                     )  
!
  ZY(:,JLAYER) =   ZHC_D_ROOF(:,JLAYER)/PTSTEP * PT_ROOF(:,JLAYER)     &
         + ZEXPL * (   ZMTC_O_D(:,JLAYER-1) * PT_ROOF(:,JLAYER-1)        &
                     - ZMTC_O_D(:,JLAYER-1) * PT_ROOF(:,JLAYER  )        &
                     - ZMTC_O_D(:,JLAYER  ) * PT_ROOF(:,JLAYER  )        &
                     + ZMTC_O_D(:,JLAYER  ) * PT_ROOF(:,JLAYER+1)        &
                   )  
END DO
!
!-------------------------------------------------------------------------------
!
!*      4.     Inside layer coefficients (with PTI_BLDWFR)
!              -------------------------------------------
!
ZA(:,IROOF_LAYER) =                                                        &
              ZIMPL * ( - ZMTC_O_D(:,IROOF_LAYER-1)                          &
                      )  

ZB(:,IROOF_LAYER) =   ZHC_D_ROOF(:,IROOF_LAYER) / PTSTEP                   &
            + ZIMPL * (   ZMTC_O_D(:,IROOF_LAYER-1)                          &
                        + ZMTC_O_D(:,IROOF_LAYER  )                          &
                      )  

ZC(:,IROOF_LAYER) =   0.
!
ZY(:,IROOF_LAYER) =     ZHC_D_ROOF(:,IROOF_LAYER)/PTSTEP                   &
                                                  * PT_ROOF(:,IROOF_LAYER)   &
                        + ZMTC_O_D(:,IROOF_LAYER) * PTI_BLDWFR(:)               &
            + ZEXPL * (   ZMTC_O_D(:,IROOF_LAYER-1)                          &
                         * PT_ROOF(:,IROOF_LAYER-1)                          &
                        - ZMTC_O_D(:,IROOF_LAYER-1)                          &
                         * PT_ROOF(:,IROOF_LAYER  )                          &
                        - ZMTC_O_D(:,IROOF_LAYER  )                          &
                         * PT_ROOF(:,IROOF_LAYER  )                          &
                      )  
!
!-------------------------------------------------------------------------------
!
!*      5.     Tri-diagonal system resolution (with PTI_BLDWFR)
!              ------------------------------------------------
!
CALL TRIDIAG_GROUND(ZA,ZB,ZC,ZY,ZX)
!
PT_ROOF(:,:) = ZX(:,:) ! calculé avec PTI_BLDWFR
!* diagnostic: computation of flux between bld and internal roof layer in the
!  case of no domestic heating i-e internal temperature not set to its min value
ZFLX_BLD_ROOF_EQ(:)=ZMTC_O_D(:,IROOF_LAYER)*(PTI_BLDWFR(:) &
                                          -ZIMPL*ZTI_ROOF(:) &
                                          -ZEXPL*PT_ROOF(:,IROOF_LAYER))  
!
!
!-------------------------------------------------------------------------------
!
!*      6.     Inside layer coefficients (with PTI_BLD)
!              ----------------------------------------
!
ZA(:,IROOF_LAYER) =                                                        &
              ZIMPL * ( - ZMTC_O_D(:,IROOF_LAYER-1)                          &
                      )  

ZB(:,IROOF_LAYER) =   ZHC_D_ROOF(:,IROOF_LAYER) / PTSTEP                   &
            + ZIMPL * (   ZMTC_O_D(:,IROOF_LAYER-1)                          &
                        + ZMTC_O_D(:,IROOF_LAYER  )                          &
                      )  

ZC(:,IROOF_LAYER) =   0.
!
ZY(:,IROOF_LAYER) =     ZHC_D_ROOF(:,IROOF_LAYER)/PTSTEP                   &
                                                  * PT_ROOF(:,IROOF_LAYER)   &
                        + ZMTC_O_D(:,IROOF_LAYER) * PTI_BLD(:)               &
            + ZEXPL * (   ZMTC_O_D(:,IROOF_LAYER-1)                          &
                         * PT_ROOF(:,IROOF_LAYER-1)                          &
                        - ZMTC_O_D(:,IROOF_LAYER-1)                          &
                         * PT_ROOF(:,IROOF_LAYER  )                          &
                        - ZMTC_O_D(:,IROOF_LAYER  )                          &
                         * PT_ROOF(:,IROOF_LAYER  )                          &
                      )  
!
!
!-------------------------------------------------------------------------------
!
!*      7.     Tri-diagonal system resolution
!              ------------------------------
!
CALL TRIDIAG_GROUND(ZA,ZB,ZC,ZY,ZX)
!
PT_ROOF(:,:) = ZX(:,:)
!
!*     diagnostic: computation of total flux between bld and internal roof layer
PFLX_BLD_ROOF(:)=ZMTC_O_D(:,IROOF_LAYER)*(PTI_BLD(:) &
                                    -ZIMPL*ZTI_ROOF(:) &
                                    -ZEXPL*PT_ROOF(:,IROOF_LAYER))  
!*     diagnostic: computation of contribution of heating parametrisation
PQF_ROOF(:)=PFLX_BLD_ROOF(:)-ZFLX_BLD_ROOF_EQ(:)
!
!-------------------------------------------------------------------------------
!
!*      8.     Infra-red radiation absorbed by roofs
!              -------------------------------------
!
!* radiative surface temperature used during the enrgy balance (linearized at t+)
ZTS_ROOF = (ZTS_ROOF(:)**4 + 4.*ZIMPL*ZTS_ROOF**3 * (PT_ROOF(:,1) - ZTS_ROOF))**0.25
!
!* absorbed LW
PABS_LW_ROOF(:) = PEMIS_ROOF(:) * (PLW(:) - XSTEFAN * ZTS_ROOF(:)** 4)
!
!-------------------------------------------------------------------------------
!
!*      9.     New saturated specified humidity near the roof surface
!              ------------------------------------------------------
!
!
PQSAT_ROOF(:) =  QSAT(PT_ROOF(:,1),PPS(:))
!
!-------------------------------------------------------------------------------
!
!*     10.     heat storage inside roofs
!              -------------------------
!
!      10.1    internal energy of the roofs at the next time step
!              --------------------------------------------------
!
ZPEI_ROOF(:)=PHC_ROOF(:,1)*PD_ROOF(:,1)*PT_ROOF(:,1)
DO JLAYER=2,IROOF_LAYER
        ZPEI_ROOF(:)=ZPEI_ROOF(:)+(PHC_ROOF(:,JLAYER)*PD_ROOF(:,JLAYER)* &
                          PT_ROOF(:,JLAYER))  
END DO
!
!      10.2    heat storage flux inside roofs 
!              ------------------------------
!
PDQS_ROOF(:)=(ZPEI_ROOF(:)-ZEI_ROOF(:))/PTSTEP
IF (LHOOK) CALL DR_HOOK('ROOF_LAYER_E_BUDGET',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------
!
END SUBROUTINE ROOF_LAYER_E_BUDGET
