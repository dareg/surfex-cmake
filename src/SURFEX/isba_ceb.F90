!#########
SUBROUTINE ISBA_CEB(IO, IP, INI, IMX, IMT, IR, DGIP, DGEIP, DGMI,      &
                    HSNOW_ISBA, HIMPLICIT_WIND, PTSTEP, PPEW_A_COEF,   &
                    PPEW_B_COEF, PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF,&
                    PPEQ_B_COEF, PSW_RAD, PLW_RAD, PEXNS, PEXNA, PTA,  &
                    PVMOD, PQA, PRR, PSR, PPS, PZREF, PUREF, PDIRCOSZW,&
                    PF5, PFFG_NOSNOW, PFFV_NOSNOW, PRHOA, PCS,         &
                    PSOILCONDZ, PSOILHCAPZ, PFROZEN1, PTDEEP_A,        &
                    PGRNDFLUX, PFLSN_COR, PSNOW_THRUFAL, PDELTA, PHUGI,&
                    PALBT, PEMIST, PDEEP_FLUX, PUSTAR, PAC_AGG, PHU_AGG )
!     ##########################################################################
!
!
!!****  *ISBA_CEB*  
!!
!!    PURPOSE
!!    -------
!       Monitor for the calculation of the surface composit energy budget
!!      call of drag, e_budget, and isba_fluxes
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
!!      B. Decharme           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/03/95 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_FIX_t, ISBA_PARAM_TIME_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t, ISBA_INIT_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
USE MODD_DIAG_n, ONLY : DIAG_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE MODD_SURF_ATM,   ONLY : LCPL_ARP
!
USE MODI_DRAG
USE MODI_E_BUDGET
USE MODI_ISBA_FLUXES
!
USE MODE_THERMOS
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!              -------------------------
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_PARAM_FIX_t), INTENT(INOUT) :: IMX
TYPE(ISBA_PARAM_TIME_t), INTENT(INOUT) :: IMT
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
TYPE(ISBA_INIT_t), INTENT(INOUT) :: INI
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
TYPE(DIAG_t), INTENT(INOUT) :: DGIP
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEIP
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DGMI
!
CHARACTER(LEN=*),     INTENT(IN) :: HSNOW_ISBA ! 'DEF' = Default F-R snow scheme
!                                              !         (Douville et al. 1995)
!                                              ! '3-L' = 3-L snow scheme (option)
!                                              !         (Boone and Etchevers 2000)
 CHARACTER(LEN=*),    INTENT(IN) :: HIMPLICIT_WIND   ! wind implicitation option
!                                                    ! 'OLD' = direct
!                                                    ! 'NEW' = Taylor serie, order 1
REAL,                 INTENT(IN) :: PTSTEP    ! timestep of the integration
REAL, DIMENSION(:),  INTENT(IN)  :: PPEW_A_COEF, PPEW_B_COEF,                   &
                                    PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF,      &
                                    PPEQ_B_COEF  
!                                  PPEW_A_COEF = A-wind coefficient (m2s/kg)
!                                  PPEW_B_COEF = B-wind coefficient (m/s)
!                                  PPET_A_COEF = A-air temperature coefficient
!                                  PPET_B_COEF = B-air temperature coefficient
!                                  PPEQ_A_COEF = A-air specific humidity coefficient
!                                  PPEQ_B_COEF = B-air specific humidity coefficient
!
REAL, DIMENSION(:), INTENT(IN)   :: PSW_RAD, PLW_RAD
!                                     PSW_RAD = incoming solar radiation
!                                     PLW_RAD = atmospheric infrared radiation

REAL, DIMENSION(:), INTENT(IN)   :: PEXNA, PEXNS, PTA, PVMOD, PQA, PRR, PSR, PPS
!                                     PEXNA= Exner function near surface atmospheric variables
!                                     PEXNS   = Exner function at the surface
!                                     PTA = 2m temperature
!                                     PVMOD = module of the horizontal wind
!                                             NOTE it should be input as >= 1. (m)
!                                     PQA = specific humidity
!                                     PPS = surface pressure
!                                     PRR = rain rate    
!                                     PSR = snow rate             
!
REAL, DIMENSION(:), INTENT(IN)   :: PZREF, PUREF
!                                     PZREF = reference height of the first
!                                             atmospheric level 
!                                     PUREF = reference height of the wind
!                                             NOTE this is different from ZZREF
!                                             ONLY in stand-alone/forced mode,
!                                             NOT when coupled to a model (MesoNH)
!
REAL, DIMENSION(:), INTENT(IN)    :: PF5
!                                     PF5 = water stress function for Hv
REAL, DIMENSION(:), INTENT(IN)    :: PDIRCOSZW 
!                                     PDIRCOSZW = Cosinus of the angle between the normal to the surface and the vertical
!
REAL, DIMENSION(:), INTENT(IN)   :: PFFG_NOSNOW, PFFV_NOSNOW
!
REAL, DIMENSION(:), INTENT(IN)   :: PRHOA
!                                   PRHOA = near-ground air density
!
REAL, DIMENSION(:), INTENT(IN)   :: PCS
!                                     PEMIS = emissivity
!                                     PCS    = heat capacity of the snow (K m2 J-1)
REAL, DIMENSION(:), INTENT(IN)   :: PGRNDFLUX, PFLSN_COR, PSNOW_THRUFAL 
!                                     PGRNDFLUX = soil/snow interface flux (W/m2) using
!                                                 ISBA-SNOW3L option
!                                     PFLSN_COR = soil/snow interface correction flux to conserve energy (W/m2)
!                                     PSNOW_THRUFAL  = snow runoff/melt leaving pack and available
!                                                  at the surface for runoff or infiltration
!                                                  [kg/(m2 s)]
REAL, DIMENSION(:,:), INTENT(IN)  :: PSOILCONDZ, PSOILHCAPZ
!                                     PSOILCONDZ= ISBA-DF Soil conductivity profile  [W/(m K)]
!                                     PSOILHCAPZ=ISBA-DF Soil heat capacity profile [J/(m3 K)]
!
REAL, DIMENSION(:), INTENT(IN)    :: PFROZEN1
!                                     PFROZEN1 = ice fraction in supurficial soil
!
REAL, DIMENSION(:), INTENT(IN)     :: PTDEEP_A
!                                      PTDEEP_A = Deep soil temperature
!                                                 coefficient depending on flux
!inout
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PDELTA
!                                      PDELTA = fraction of the foliage covered
!                                              by intercepted water
!out
!
REAL, DIMENSION(:), INTENT(OUT)  :: PHUGI
!                                     PHUGI = ground (ice) relative humidity
!
REAL, DIMENSION(:), INTENT(OUT)  :: PALBT, PEMIST
!                                     PALBT  = averaged albedo
!                                     PEMIST = averaged emissivity
!
REAL, DIMENSION(:), INTENT(OUT)   :: PDEEP_FLUX ! Heat flux at bottom of ISBA (W/m2)
!
REAL, DIMENSION(:), INTENT(OUT)     :: PUSTAR
!                                      PUSTAR  = friction velocity
!
REAL, DIMENSION(:),  INTENT(OUT) :: PAC_AGG  ! aggregated aerodynamic conductance
                                     ! for evaporative flux calculations
REAL, DIMENSION(:),  INTENT(OUT) :: PHU_AGG  ! aggregated relative humidity
                                    ! for evaporative flux calculations
!
!
!*      0.2    declarations of local parameters
!
REAL, PARAMETER            :: ZDEPTH_COR = 0.6
!                             ZDEPTH_COR = depth over which the correction flux is applied
!
REAL, PARAMETER            :: ZDTG1_COR = 10.0 
!                             ZDTG1_COR = Delta temperature limit to comput the correction flux (K)
!
!
!*      0.3    declarations of local variables
!
REAL, DIMENSION(SIZE(PTA)) :: ZQSAT     ! expression for the saturation specific humidity 
!
REAL, DIMENSION(SIZE(PTA)) :: ZDQSAT    ! expression for the saturation specific humidity derivative
!
REAL, DIMENSION(SIZE(PTA)) :: ZTA_IC, ZQA_IC, ZUSTAR2_IC ! TA, QA and friction updated values
!                                                        ! if implicit coupling with atmosphere used.
!
REAL, DIMENSION(SIZE(PTA)) :: ZLEG_DELTA  ! soil evaporation delta fn
REAL, DIMENSION(SIZE(PTA)) :: ZLEGI_DELTA ! soil sublimation delta fn
!
REAL, DIMENSION(SIZE(PTA)) :: ZT2M     ! restore temperature before time integration (K)
REAL, DIMENSION(SIZE(PTA)) :: ZTSM     ! surface temperature before time integration (K)
!
REAL, DIMENSION(SIZE(IR%XTG(:,:,1),1),SIZE(IR%XTG(:,:,1),2)) :: ZFLUX_COR, ZLAYERHCAP
!                                           ZFLUX_COR = correction flux by layer to conserve energy (W/m2)
!
REAL, DIMENSION(SIZE(PTA)) :: ZGRNDFLUX, ZTOTALHCAP, ZWORK
!
INTEGER                    :: INJ, INL, JI, JL
LOGICAL                    :: LEXPLICIT_SNOW
!
!*      0.4    declarations of local time spliting variables
!
! Working arrays for flux averaging over time split
!
REAL, DIMENSION(SIZE(PTA)) :: ZDEEP_FLUX, ZLE_FLOOD, ZLEI_FLOOD, &
                              ZRN, ZH, ZLE, ZLEG, ZLEV,  &
                              ZLES, ZLER, ZLETR, ZEVAP,      &
                              ZGFLUX, ZMELTADV, ZMELT,           &
                              ZRESTORE, ZLEGI, ZUSTAR2,          &
                              ZAC_AGG, ZHU_AGG
!
REAL, DIMENSION(SIZE(PTA)) :: ZDEEP_FLUX_SUM, ZLE_FLOOD_SUM, ZLEI_FLOOD_SUM, &
                              ZRN_SUM, ZH_SUM, ZLE_SUM, ZLEG_SUM, ZLEV_SUM,  &
                              ZLES_SUM, ZLER_SUM, ZLETR_SUM, ZEVAP_SUM,      &
                              ZGFLUX_SUM, ZMELTADV_SUM, ZMELT_SUM,           &
                              ZRESTORE_SUM, ZLEGI_SUM, ZUSTAR2_SUM,          &
                              ZAC_AGG_SUM, ZHU_AGG_SUM
!
REAL, PARAMETER            :: ZTSPLIT  = 300. ! s Minimum time tstep required to time-split energy budget
INTEGER                    :: ITSPLIT, JSPLIT
REAL                       :: ZTSTEP, ZNSPLIT
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ISBA_CEB',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!*      1.0    Preliminaries
!              -------------
!
INJ=SIZE(IR%XTG(:,:,1),1)
INL=SIZE(IR%XTG(:,:,1),2)
!
!local init
!
ZQSAT      (:) = XUNDEF
ZTA_IC     (:) = XUNDEF
ZQA_IC     (:) = XUNDEF
ZUSTAR2_IC (:) = XUNDEF
ZLEG_DELTA (:) = XUNDEF
ZLEGI_DELTA(:) = XUNDEF
ZDQSAT     (:) = XUNDEF
!
ZDEEP_FLUX(:) = XUNDEF
ZLE_FLOOD (:) = XUNDEF
ZLEI_FLOOD(:) = XUNDEF
ZRN       (:) = XUNDEF
ZH        (:) = XUNDEF
ZLE       (:) = XUNDEF
ZLEG      (:) = XUNDEF
ZLEV      (:) = XUNDEF
ZLES      (:) = XUNDEF
ZLER      (:) = XUNDEF
ZLETR     (:) = XUNDEF
ZEVAP     (:) = XUNDEF
ZGFLUX    (:) = XUNDEF
ZMELTADV  (:) = XUNDEF
ZMELT     (:) = XUNDEF
ZRESTORE  (:) = XUNDEF
ZLEGI     (:) = XUNDEF
ZAC_AGG   (:) = XUNDEF
ZHU_AGG   (:) = XUNDEF
!
ZUSTAR2_SUM   (:) = 0.0
ZEVAP_SUM     (:) = 0.0
!
ZRN_SUM       (:) = 0.0
ZH_SUM        (:) = 0.0
ZGFLUX_SUM    (:) = 0.0
ZLE_SUM       (:) = 0.0
!
ZLEG_SUM      (:) = 0.0
ZLEGI_SUM     (:) = 0.0
ZLEV_SUM      (:) = 0.0
ZLES_SUM      (:) = 0.0
ZLER_SUM      (:) = 0.0
ZLETR_SUM     (:) = 0.0
ZLE_FLOOD_SUM (:) = 0.0
ZLEI_FLOOD_SUM(:) = 0.0
!
ZDEEP_FLUX_SUM(:) = 0.0
ZMELTADV_SUM  (:) = 0.0
ZMELT_SUM     (:) = 0.0
ZRESTORE_SUM  (:) = 0.0
ZAC_AGG_SUM   (:) = 0.0
ZHU_AGG_SUM   (:) = 0.0
!
!-------------------------------------------------------------------------------
!
ZGRNDFLUX(:  ) = PGRNDFLUX(:)
!
ZFLUX_COR(:,:) = 0.0
!
!-------------------------------------------------------------------------------
!
!*      2.0    Correction flux to conserv energy budget
!              ----------------------------------------
!
LEXPLICIT_SNOW=(HSNOW_ISBA == '3-L' .OR. HSNOW_ISBA == 'CRO')
!
IF(LEXPLICIT_SNOW.AND.IO%CISBA/='DIF')THEN
!
  ZFLUX_COR(:,1)=IR%XPSN(:,1)*PFLSN_COR(:)
!
ELSEIF(LEXPLICIT_SNOW.AND.IO%CISBA=='DIF')THEN
!
  ZLAYERHCAP  (:,:) = 0.0
  ZTOTALHCAP    (:) = 0.0
!
! To conserv energy, the correction flux is distributed at least
! over the first layers of the soil, ZDEPTH_COR. This method prevent 
! numerical oscillations especially when explicit snow vanishes
!
  ZWORK(:) = MIN(IMX%XDG(:,INL,1),ZDEPTH_COR)
!
  ZLAYERHCAP(:,1)= 1.0/DGMI%XCT(:)
  ZTOTALHCAP(:  )= 1.0/DGMI%XCT(:)
  DO JL=2,INL
    DO JI=1,INJ
      ZLAYERHCAP(JI,JL) = PSOILHCAPZ(JI,JL) * MIN( IP%XDZG(JI,JL,1), &
                MAX(0.0,ZWORK(JI)-IMX%XDG(JI,JL,1)+IP%XDZG(JI,JL,1)) )
      ZTOTALHCAP(JI   ) = ZTOTALHCAP(JI) + ZLAYERHCAP(JI,JL)
    ENDDO
  ENDDO
!
  DO JL=1,INL
    DO JI=1,INJ
      IF(ZTOTALHCAP(JI)>0.0)THEN
        ZFLUX_COR(JI,JL) = IR%XPSN(JI,1)*PFLSN_COR(JI)*ZLAYERHCAP(JI,JL)/ZTOTALHCAP(JI)
      ENDIF
    ENDDO
  ENDDO
!
! The second correction is computed if the delta temperature
! due to snow/soil ground flux is superior to ZDTG1_COR (K)
! Especially relevant when PPSN ~ 1 over vegetated area
!
  ZWORK(:)=PTSTEP*DGMI%XCT(:)*IR%XPSN(:,1)*ABS(PGRNDFLUX(:))
!
  WHERE(ZTOTALHCAP(:)>0.0.AND.ZWORK(:)>=ZDTG1_COR)
    ZGRNDFLUX(:) = PGRNDFLUX(:)*ZLAYERHCAP(:,1)/ZTOTALHCAP(:)
  ENDWHERE
!
  DO JL=2,INL
    DO JI=1,INJ
      IF(ZTOTALHCAP(JI)>0.0.AND.ZWORK(JI)>=ZDTG1_COR)THEN
        ZFLUX_COR(JI,JL)=ZFLUX_COR(JI,JL)+IR%XPSN(JI,1)*PGRNDFLUX(JI) &
                                         *ZLAYERHCAP(JI,JL)/ZTOTALHCAP(JI)
      ENDIF
    ENDDO
  ENDDO
!
ENDIF
!
!-------------------------------------------------------------------------------
!
IF(LCPL_ARP)THEN
  ITSPLIT = 1
ELSE
  ITSPLIT  = MAX(1,NINT(PTSTEP/ZTSPLIT))  ! number of split-time steps
ENDIF
!
ZNSPLIT  = REAL(ITSPLIT)
!
ZTSTEP  = PTSTEP/ZNSPLIT             ! split time step
!
!-------------------------------------------------------------------------------
!
DO JSPLIT=1,ITSPLIT
!
!  Save surface and sub-surface temperature values at beginning of time step for 
!  budget and flux calculations:
!
   ZTSM(:) = IR%XTG(:,1,1)
   ZT2M(:) = IR%XTG(:,2,1)
!
!
!*      3.0    Aerodynamic drag and heat transfer coefficients
!              -----------------------------------------------
!
  !print*,JSPLIT,'tgc1 ',IR%XTG(1,1,1)
   CALL DRAG(IO%CISBA, HSNOW_ISBA, IO%CCPSURF, ZTSTEP, IR%XTG(:,1,1), IR%XWG(:,1,1), &
             IR%XWGI(:,1,1), PEXNS, PEXNA, PTA, PVMOD, PQA, PRR, PSR, PPS, DGMI%XRS, &
             IMT%XVEG(:,1), DGIP%XZ0, DGIP%XZ0EFF, DGIP%XZ0H, IP%XWFC(:,1),          &
             IP%XWSAT(:,1), IR%XPSNG(:,1), IR%XPSNV(:,1), PZREF, PUREF, PDIRCOSZW,   &
             PDELTA, PF5, IR%XRESA(:,1), DGIP%XCH, DGIP%XCD, DGIP%XCDN, DGIP%XRI,    &
             DGIP%XHUG, PHUGI, DGMI%XHV, DGIP%XHU, IP%XPCPS(:,1), DGIP%XQS,          &
             INI%XFFG(:,1), INI%XFFV(:,1), INI%XFF(:,1), PFFG_NOSNOW, PFFV_NOSNOW,   &
             ZLEG_DELTA, ZLEGI_DELTA, IR%XWR(:,1), PRHOA, IP%XPLVTT(:,1), PQSAT=ZQSAT )  
!
!
  !print*,JSPLIT,'tgc2 ',IR%XTG(1,1,1)
!
!*      4.0    Resolution of the surface and soil energy budget
!              ------------------------------------------------
!
   CALL E_BUDGET(IO, IP, INI, IMT, IR, DGIP, DGMI, HSNOW_ISBA, HIMPLICIT_WIND,  &
                 ZTSTEP, PUREF, PPEW_A_COEF, PPEW_B_COEF, PPET_A_COEF,          &
                 PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF, PVMOD, ZTSM, ZT2M,      &
                 PSW_RAD, PLW_RAD, PTA, PQA, PPS, PRHOA, PEXNS, PEXNA,          &
                 PHUGI, ZLEG_DELTA, ZLEGI_DELTA,  ZGRNDFLUX, ZFLUX_COR,         &
                 PSOILCONDZ, PSOILHCAPZ, PALBT, PEMIST, ZQSAT, ZDQSAT,          &
                 PFROZEN1, PTDEEP_A, ZTA_IC, ZQA_IC, ZUSTAR2_IC, ZDEEP_FLUX,    &
                 ZRESTORE          )
!
 !print*,JSPLIT,'tgc3 ',IR%XTG(1,1,1)
!
!*      5.0    Energy and momentum fluxes
!              --------------------------
!
!*******************************************************************************
! WARNING: at this stage, ZALBT and ZEMIST have two different meanings according
!          to the ISBA snow-scheme option:
!  'D95' : they represent aggregated (snow + flood + snow-flood-free) albedo and emissivity
!  '3-L' : they represent                    flood + snow-flood-free  albedo and emissivity
!*******************************************************************************
!
   CALL ISBA_FLUXES(IO, IP, INI, IMT, IR, DGMI, HSNOW_ISBA, ZTSTEP,           &
                    PSW_RAD, PLW_RAD, ZTA_IC, ZQA_IC, PRHOA, PEXNS, PEXNA,    &
                    DGIP%XHUG, PHUGI, ZLEG_DELTA, ZLEGI_DELTA, PDELTA, PF5,   &
                    PCS, ZTSM, ZT2M, PFROZEN1, PALBT, PEMIST, ZQSAT, ZDQSAT,  &
                    PSNOW_THRUFAL, ZRN, ZH, ZLE, ZLEG, ZLEGI, ZLEV, ZLES,     &
                    ZLER, ZLETR, ZEVAP, ZGFLUX, ZMELTADV, ZMELT, PSOILCONDZ,  &
                    ZLE_FLOOD, ZLEI_FLOOD ) 
!
 !print*,JSPLIT,'tgc4 ',IR%XTG(1,1,1)
!
!*      6.0    Aggregated coefficients
!              -----------------------
!
!  Compute aggregated coefficients for evaporation
!  Sum(LEV+LEG+LEGI+LES) = ACagg * Lv * RHOA * (HUagg.Qsat - Qa)
!
   ZAC_AGG(:) =   1. / IR%XRESA(:,1) / IP%XPLVTT(:,1)     &
             * ( IP%XPLVTT(:,1) *    IMT%XVEG(:,1) * (1.-IR%XPSNV(:,1)) * DGMI%XHV(:)   &
               + IP%XPLVTT(:,1) *(1.-IMT%XVEG(:,1))* (1.-IR%XPSNG(:,1)) * (1.-PFROZEN1(:))  &
               + IP%XPLSTT(:,1) *(1.-IMT%XVEG(:,1))* (1.-IR%XPSNG(:,1)) * PFROZEN1(:)  &
               + IP%XPLSTT(:,1) *                         IR%XPSN(:,1)   )  
!
   WHERE(ZAC_AGG(:)>0.0)
     ZHU_AGG(:) =   1. / (IR%XRESA(:,1) * ZAC_AGG(:)) / IP%XPLVTT(:,1)         &
                  * ( IP%XPLVTT(:,1)*    IMT%XVEG(:,1) *(1.-IR%XPSNV(:,1))                 *DGMI%XHV(:)   &
                    + IP%XPLVTT(:,1)*(1.-IMT%XVEG(:,1))*(1.-IR%XPSNG(:,1))*(1.-PFROZEN1(:))*DGIP%XHUG(:)  &
                    + IP%XPLSTT(:,1)*(1.-IMT%XVEG(:,1))*(1.-IR%XPSNG(:,1))*    PFROZEN1(:) *PHUGI(:) &
                    + IP%XPLSTT(:,1)*                 IR%XPSN(:,1)                            )  
   ENDWHERE
!
   ZUSTAR2_SUM   (:) = ZUSTAR2_SUM   (:) + ZUSTAR2_IC(:)
!
   ZEVAP_SUM     (:) = ZEVAP_SUM     (:) + ZEVAP     (:)
!
   ZRN_SUM       (:) = ZRN_SUM       (:) + ZRN       (:)
   ZH_SUM        (:) = ZH_SUM        (:) + ZH        (:)
   ZGFLUX_SUM    (:) = ZGFLUX_SUM    (:) + ZGFLUX    (:)
   ZLE_SUM       (:) = ZLE_SUM       (:) + ZLE       (:)
!
   ZLEG_SUM      (:) = ZLEG_SUM      (:) + ZLEG      (:)
   ZLEGI_SUM     (:) = ZLEGI_SUM     (:) + ZLEGI     (:)   
   ZLEV_SUM      (:) = ZLEV_SUM      (:) + ZLEV      (:)
   ZLES_SUM      (:) = ZLES_SUM      (:) + ZLES      (:)
   ZLER_SUM      (:) = ZLER_SUM      (:) + ZLER      (:)
   ZLETR_SUM     (:) = ZLETR_SUM     (:) + ZLETR     (:)
   ZLE_FLOOD_SUM (:) = ZLE_FLOOD_SUM (:) + ZLE_FLOOD (:)
   ZLEI_FLOOD_SUM(:) = ZLEI_FLOOD_SUM(:) + ZLEI_FLOOD(:)   
!
   ZDEEP_FLUX_SUM(:) = ZDEEP_FLUX_SUM(:) + ZDEEP_FLUX(:)
   ZMELTADV_SUM  (:) = ZMELTADV_SUM  (:) + ZMELTADV  (:)
   ZMELT_SUM     (:) = ZMELT_SUM     (:) + ZMELT     (:)
   ZRESTORE_SUM  (:) = ZRESTORE_SUM  (:) + ZRESTORE  (:)
   ZAC_AGG_SUM   (:) = ZAC_AGG_SUM   (:) + ZAC_AGG   (:)
   ZHU_AGG_SUM   (:) = ZHU_AGG_SUM   (:) + ZHU_AGG   (:)
!
!-------------------------------------------------------------------------------
!
ENDDO
!
PUSTAR    (:) = SQRT(ZUSTAR2_SUM(:)/ZNSPLIT)
!
DGIP%XEVAP     (:) = ZEVAP_SUM     (:) / ZNSPLIT
!
DGIP%XRN       (:) = ZRN_SUM       (:) / ZNSPLIT
DGIP%XH        (:) = ZH_SUM        (:) / ZNSPLIT
DGIP%XGFLUX    (:) = ZGFLUX_SUM    (:) / ZNSPLIT
IR%XLE       (:,1) = ZLE_SUM       (:) / ZNSPLIT
!
DGEIP%XLEG      (:) = ZLEG_SUM      (:) / ZNSPLIT
DGEIP%XLEGI     (:) = ZLEGI_SUM     (:) / ZNSPLIT
DGEIP%XLEV      (:) = ZLEV_SUM      (:) / ZNSPLIT
DGEIP%XLES      (:) = ZLES_SUM      (:) / ZNSPLIT
DGEIP%XLER      (:) = ZLER_SUM      (:) / ZNSPLIT
DGEIP%XLETR     (:) = ZLETR_SUM     (:) / ZNSPLIT
DGEIP%XLE_FLOOD (:) = ZLE_FLOOD_SUM (:) / ZNSPLIT
DGEIP%XLEI_FLOOD(:) = ZLEI_FLOOD_SUM(:) / ZNSPLIT
!
PDEEP_FLUX(:) = ZDEEP_FLUX_SUM(:) / ZNSPLIT
DGEIP%XMELTADV  (:) = ZMELTADV_SUM  (:) / ZNSPLIT
DGEIP%XMELT     (:) = ZMELT_SUM     (:) / ZNSPLIT
DGEIP%XRESTORE  (:) = ZRESTORE_SUM  (:) / ZNSPLIT
PAC_AGG   (:) = ZAC_AGG_SUM   (:) / ZNSPLIT
PHU_AGG   (:) = ZHU_AGG_SUM   (:) / ZNSPLIT
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ISBA_CEB',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ISBA_CEB
