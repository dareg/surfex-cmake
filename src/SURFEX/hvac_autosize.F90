!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #############################################################
SUBROUTINE HVAC_AUTOSIZE (B, BOP, G, T, TOP, KI,KLUOUT)
!     #############################################################
!!
!!    PURPOSE
!!    -------
!!
!!    Calibrates HVAC systems for TEB-BEM
!!      
!!    AUTHOR
!!    ------
!!
!!      G. Pigeon           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!    Original    05/2011
!!    modified    08/2013 add solar panels (V. Masson)
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_BEM_n, ONLY : BEM_t
USE MODD_BEM_OPTION_n, ONLY : BEM_OPTIONS_t
USE MODD_GRID_n, ONLY : GRID_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_IRRIG_n, ONLY : TEB_IRRIG_t, TEB_IRRIG_INIT
USE MODD_DIAG_MISC_TEB_n, ONLY : DIAG_MISC_TEB_t, DIAG_MISC_TEB_INIT
USE MODD_DIAG_MISC_TEB_OPTION_n, ONLY : DIAG_MISC_TEB_OPTIONS_t, DIAG_MISC_TEB_OPTIONS_INIT
USE MODD_DIAG_n, ONLY : DIAG_t, DIAG_INIT
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t, DIAG_EVAP_ISBA_INIT
!
USE MODD_CSTS,  ONLY : XCPD, XPI, XP00, XRD, XSTEFAN
!
USE MODI_TEB
USE MODI_SUNPOS
USE MODI_SW_DAYCYCLE
USE MODI_URBAN_LW_COEF
USE MODI_URBAN_SOLAR_ABS
USE MODI_GET_SIZES_PARALLEL
USE MODI_DIAG_MISC_TEB_INIT_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
#ifdef AIX64 
USE OMP_LIB
#endif
!
IMPLICIT NONE
!
#ifndef AIX64
  INCLUDE 'omp_lib.h'
#endif
!
!*    0.1    Declaration of arguments
!            ------------------------
!
TYPE(BEM_t), INTENT(INOUT) :: B
TYPE(BEM_OPTIONS_t), INTENT(INOUT) :: BOP
TYPE(GRID_t), INTENT(INOUT) :: G
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
!
INTEGER,       INTENT(IN)    :: KI     ! number of points
INTEGER,       INTENT(IN)    :: KLUOUT ! output listing logical unit
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
!-------------------------------------------------------------------------
!
TYPE(DIAG_t) :: YGRD
TYPE(DIAG_EVAP_ISBA_t) :: YGRDE
TYPE(DIAG_MISC_TEB_t) :: YDMTC
TYPE(DIAG_MISC_TEB_t) :: YDMT
TYPE(DIAG_MISC_TEB_OPTIONS_t) :: YDMTO
!
! Options TOP
!
LOGICAL :: GCANOPY
 CHARACTER(LEN=3) :: YBEM
 CHARACTER(LEN=6) :: YZ0H 
 CHARACTER(LEN=5) :: YCH_BEM  
!
! Teb Fields T
!
REAL, DIMENSION(KI) :: ZGARDEN
!
! Bem Options BOP
!
 CHARACTER(LEN=6)    :: YCOOL_COIL
 CHARACTER(LEN=6)    :: YHEAT_COIL
!
! Bem!
 CHARACTER(LEN=4),DIMENSION(KI) :: YNATVENT
REAL, DIMENSION(KI) :: ZF_WATER_COND
REAL, DIMENSION(KI) :: ZAUX_MAX
REAL, DIMENSION(KI) :: ZT_WIN2
LOGICAL, DIMENSION(KI) :: GNATVENT_NIGHT
!
! Road irrigation (not used here)
!
TYPE(TEB_IRRIG_t) :: YIR
!
!-------------------------------------------------------------------------
!
REAL, PARAMETER     :: ZTSTEP = 300.
INTEGER, PARAMETER  :: JPYEAR = 2004   ! Current year (UTC)
!
!local variable
!
INTEGER             :: IMONTH = 7  ! Current month (UTC)
INTEGER             :: IDAY   = 12 ! Current day (UTC)
REAL                :: ZTIME  = 0.  ! Time at start of the run (s)
!
 CHARACTER(LEN=3)    :: YIMPLICIT_WIND = 'NEW'
!
! Loop control indexes
!
INTEGER             :: JJ
INTEGER             :: JFORC_STEP
INTEGER             :: INB_STEP_ATM
!
! Intermediate variables
!
REAL, DIMENSION(KI) :: ZU_ROOF
REAL, DIMENSION(KI) :: ZU_WALL
!
REAL,DIMENSION(KI) :: ZT1    ! intermediate variable
REAL,DIMENSION(KI) :: ZTN    ! intermediate variable
!
REAL, DIMENSION(KI) :: ZT_SKY
!
REAL, DIMENSION(KI) :: ZTOT_SW
REAL, DIMENSION(KI) :: ZTOUT_EQ
!
! Arguments to TEB
!
REAL, DIMENSION(KI) :: ZT_CANYON
REAL, DIMENSION(KI) :: ZQ_CANYON
REAL, DIMENSION(KI) :: ZU_CANYON
REAL, DIMENSION(KI) :: ZZ_LOWCAN
!
REAL, DIMENSION(KI) :: ZPS
REAL, DIMENSION(KI) :: ZPA
REAL, DIMENSION(KI) :: ZEXNS
REAL, DIMENSION(KI) :: ZEXNA
REAL, DIMENSION(KI) :: ZTA
REAL, DIMENSION(KI) :: ZQA
REAL, DIMENSION(KI) :: ZRHOA
REAL, DIMENSION(KI) :: ZLW_RAD
REAL, DIMENSION(KI) :: ZRR
REAL, DIMENSION(KI) :: ZSR
REAL, DIMENSION(KI) :: ZZREF
!
REAL, DIMENSION(KI) :: ZPEW_A_COEF
REAL, DIMENSION(KI) :: ZPEW_B_COEF
!
REAL, DIMENSION(KI) :: ZTS_GARDEN
!
REAL, DIMENSION(KI) :: ZDF_ROOF
REAL, DIMENSION(KI) :: ZDN_ROOF
REAL, DIMENSION(KI) :: ZDN_ROAD
REAL, DIMENSION(KI) :: ZDF_ROAD
!
REAL, DIMENSION(KI) :: ZQSAT_ROOF
REAL, DIMENSION(KI) :: ZQSAT_ROAD
REAL, DIMENSION(KI) :: ZDELT_ROOF
REAL, DIMENSION(KI) :: ZDELT_ROAD
!
 CHARACTER(LEN=4)    :: YSNOW_ROOF
 CHARACTER(LEN=4)    :: YSNOW_ROAD
!
REAL, DIMENSION(KI) :: ZLW_WA_TO_WB   ! longwave exchange coefficients
REAL, DIMENSION(KI) :: ZLW_WA_TO_R
REAL, DIMENSION(KI) :: ZLW_WB_TO_R
REAL, DIMENSION(KI) :: ZLW_WA_TO_NR
REAL, DIMENSION(KI) :: ZLW_WB_TO_NR
REAL, DIMENSION(KI) :: ZLW_R_TO_WA
REAL, DIMENSION(KI) :: ZLW_R_TO_WB
REAL, DIMENSION(KI) :: ZLW_G_TO_WA
REAL, DIMENSION(KI) :: ZLW_G_TO_WB
REAL, DIMENSION(KI) :: ZLW_S_TO_WA
REAL, DIMENSION(KI) :: ZLW_S_TO_WB
REAL, DIMENSION(KI) :: ZLW_S_TO_R
REAL, DIMENSION(KI) :: ZLW_S_TO_NR
REAL, DIMENSION(KI) :: ZLW_NR_TO_WA
REAL, DIMENSION(KI) :: ZLW_NR_TO_WB
REAL, DIMENSION(KI) :: ZLW_NR_TO_WIN
REAL, DIMENSION(KI) :: ZLW_WA_TO_WIN
REAL, DIMENSION(KI) :: ZLW_WB_TO_WIN
REAL, DIMENSION(KI) :: ZLW_G_TO_WIN
REAL, DIMENSION(KI) :: ZLW_R_TO_WIN
REAL, DIMENSION(KI) :: ZLW_S_TO_WIN
REAL, DIMENSION(KI) :: ZLW_WIN_TO_WA
REAL, DIMENSION(KI) :: ZLW_WIN_TO_WB
REAL, DIMENSION(KI) :: ZLW_WIN_TO_R
REAL, DIMENSION(KI) :: ZLW_WIN_TO_NR
!
REAL, DIMENSION(KI) :: ZLEW_ROOF     ! latent heat flux over roof (snow)
!
REAL, DIMENSION(KI) :: ZLEW_ROAD     ! latent heat flux over road (snow)
!
REAL, DIMENSION(KI) :: ZLE_WALL_A    ! latent heat flux over wall
REAL, DIMENSION(KI) :: ZLE_WALL_B    ! latent heat flux over wall
 
!
REAL, DIMENSION(KI) :: ZRNSNOW_ROOF  ! net radiation over snow
REAL, DIMENSION(KI) :: ZHSNOW_ROOF   ! sensible heat flux over snow
REAL, DIMENSION(KI) :: ZLESNOW_ROOF  ! latent heat flux over snow
REAL, DIMENSION(KI) :: ZGSNOW_ROOF   ! flux under the snow
REAL, DIMENSION(KI) :: ZMELT_ROOF    ! snow melt
!
REAL, DIMENSION(KI) :: ZRNSNOW_ROAD  ! net radiation over snow
REAL, DIMENSION(KI) :: ZHSNOW_ROAD   ! sensible heat flux over snow
REAL, DIMENSION(KI) :: ZLESNOW_ROAD  ! latent heat flux over snow
REAL, DIMENSION(KI) :: ZGSNOW_ROAD   ! flux under the snow
REAL, DIMENSION(KI) :: ZMELT_ROAD    ! snow melt
!
REAL, DIMENSION(KI) :: ZUW_ROAD      ! Momentum flux for roads
REAL, DIMENSION(KI) :: ZUW_ROOF      ! Momentum flux for roofs
REAL, DIMENSION(KI) :: ZDUWDU_ROAD   !
REAL, DIMENSION(KI) :: ZDUWDU_ROOF   !
REAL, DIMENSION(KI) :: ZUSTAR_TOWN   ! friction velocity over town
!
REAL, DIMENSION(KI) :: ZCD           ! town averaged drag coefficient
REAL, DIMENSION(KI) :: ZCDN          ! town averaged neutral drag coefficient
REAL, DIMENSION(KI) :: ZCH_TOWN      ! town averaged heat transfer coefficient
REAL, DIMENSION(KI) :: ZRI_TOWN      ! town averaged Richardson number
REAL, DIMENSION(KI) :: ZRESA_TOWN    ! town aerodynamical resistance
REAL, DIMENSION(KI) :: ZDQS_TOWN     ! heat storage inside town
REAL, DIMENSION(KI) :: ZQF_TOWN      ! total anthropogenic heat
REAL, DIMENSION(KI) :: ZQF_BLD       ! anthropogenic heat flux of domestic heating  
REAL, DIMENSION(KI) :: ZFLX_BLD      ! heat flux between inside of the bld
REAL, DIMENSION(KI) :: ZAC_ROOF      ! roof conductance
REAL, DIMENSION(KI) :: ZAC_ROAD       ! road conductance
REAL, DIMENSION(KI) :: ZAC_WALL      ! wall conductance
REAL, DIMENSION(KI) :: ZAC_TOP       ! top conductance
REAL, DIMENSION(KI) :: ZAC_GARDEN     ! garden conductance
REAL, DIMENSION(KI) :: ZAC_ROOF_WAT  ! roof water conductance
REAL, DIMENSION(KI) :: ZAC_ROAD_WAT  ! roof water conductance 
REAL, DIMENSION(KI) :: ZEMIT_LW_FAC
REAL, DIMENSION(KI) :: ZEMIT_LW_ROAD
REAL, DIMENSION(KI) :: ZT_RAD_IND
REAL, DIMENSION(KI) :: ZHU_BLD
REAL, DIMENSION(KI) :: ZTSUN
!
! Arguments to urban_solar_abs
!
REAL, DIMENSION(KI) :: ZZENITH
REAL, DIMENSION(KI) :: ZAZIM
!
REAL, DIMENSION(KI) :: ZDIR_SW
REAL, DIMENSION(KI) :: ZSCA_SW
!
REAL, DIMENSION(KI) :: ZFRAC_PANEL
REAL, DIMENSION(KI) :: ZALB_PANEL
REAL, DIMENSION(KI) :: ZSVF_GARDEN
!
REAL, DIMENSION(KI) :: ZREC_SW_ROAD
REAL, DIMENSION(KI) :: ZREC_SW_SNOW_ROAD
REAL, DIMENSION(KI) :: ZREC_SW_WALL_A
REAL, DIMENSION(KI) :: ZREC_SW_WALL_B
REAL, DIMENSION(KI) :: ZREC_SW_GARDEN
REAL, DIMENSION(KI) :: ZREC_SW_ROOF
!
REAL, DIMENSION(KI) :: ZDIR_ALB_TOWN
REAL, DIMENSION(KI) :: ZSCA_ALB_TOWN
REAL, DIMENSION(KI) :: ZSW_RAD_GARDEN
REAL, DIMENSION(KI) :: ZREC_SW_WIN
REAL, DIMENSION(KI) :: ZREF_SW_GRND
REAL, DIMENSION(KI) :: ZREF_SW_FAC
!
REAL, DIMENSION(KI) :: ZALB_GREENROOF
REAL, DIMENSION(KI) :: ZALB_GARDEN
!
! Arguments to urban_lw_coef
!
REAL, DIMENSION(KI) :: ZEMIS_GARDEN
!
REAL, DIMENSION(KI) :: ZLW_WA_TO_G
REAL, DIMENSION(KI) :: ZLW_WB_TO_G
REAL, DIMENSION(KI) :: ZLW_S_TO_G
REAL, DIMENSION(KI) :: ZLW_WIN_TO_G
!
!! GREGOIRE 13/03
!
!
!new for shading
REAL, DIMENSION(KI) :: ZE_SHADING
LOGICAL, DIMENSION(KI) :: GSHAD_DAY
LOGICAL, DIMENSION(KI) :: GSHADE
!
!
! Case greenroof
REAL, DIMENSION(KI) :: ZUW_GREENROOF
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*      0.     Initialization
!              --------------
!
IF (LHOOK) CALL DR_HOOK('HVAC_AUTOSIZE',0,ZHOOK_HANDLE)
!
!    Design parameters
!
ZPS   = 101325.
ZPA = ZPS
ZEXNS = (ZPS/XP00)**(XRD/XCPD)
ZEXNA = (ZPA/XP00)**(XRD/XCPD)
ZQA   = 0.011
ZLW_RAD= 300.
ZRR = 0.0 
ZSR = 0.0 
ZZREF =  50.
!
ZPEW_A_COEF  = 0.5  
ZPEW_B_COEF  = 0.5 
!
! initial value for air temperature and outdoor wall/roof/window/road temperature
ZT_CANYON(:) = 10.7/2 * SIN(2*XPI/(24*3600) * (ZTIME+16*3600)) + (B%XT_SIZE_MAX(:)-10.7/2)
!
ZQ_CANYON = 0.011
ZU_CANYON = 2.5 
ZZ_LOWCAN = ZZREF
!
!
ZTS_GARDEN = 300.
!
ZDN_ROOF = 0.0
ZDF_ROOF = 1.0
ZDN_ROAD = 0.0
ZDF_ROAD = 1.0
!
ZQSAT_ROOF = 0.015
ZQSAT_ROAD = 0.015
ZDELT_ROOF = 0.0 
ZDELT_ROAD = 0.0 
!
!
!
!*      A.     Autosize of the heating system
!              ---------------------------------
!
ZU_ROOF(:) = 0.0
DO JJ=1,TOP%NROOF_LAYER
  ZU_ROOF(:) = ZU_ROOF(:) + T%XD_ROOF(:,JJ)/T%XTC_ROOF(:,JJ) 
END DO
ZU_ROOF(:) = ZU_ROOF(:) + 1./10. + 1./25.         
ZU_ROOF(:) = 1. / ZU_ROOF(:)
!
ZU_WALL(:) = 0.0
DO JJ=1,TOP%NWALL_LAYER
  ZU_WALL(:) = ZU_WALL(:) + T%XD_WALL(:,JJ)/T%XTC_WALL(:,JJ)
END DO
ZU_WALL(:) = ZU_WALL(:) + 1./10. + 1./25.         
ZU_WALL(:) = 1. / ZU_WALL(:)
!
ZRHOA = 1.30
ZT_SKY = 253.15
ZTOUT_EQ(:) = (B%XT_SIZE_MIN(:) + ZT_SKY(:))/2.
!
!   Heating Coil Capacity [W m-2(bld)]
B%XCAP_SYS_HEAT(:) = ZU_WALL(:) * T%XWALL_O_BLD(:) * (B%XTHEAT_TARGET(:) - ZTOUT_EQ(:)) &
                 + B%XU_WIN(:)  * B%XGLAZ_O_BLD(:) * (B%XTHEAT_TARGET(:) - ZTOUT_EQ(:)) &
                 + ZU_ROOF(:)                      * (B%XTHEAT_TARGET(:) - ZTOUT_EQ(:)) &
                 - B%XQIN(:) * T%XBLD_HEIGHT(:) / B%XFLOOR_HEIGHT(:)*          &
                   (1 - B%XQIN_FLAT(:))                                        &
                 + B%XINF(:)    * T%XBLD_HEIGHT(:) / 3600* ZRHOA(:) * XCPD *   &
                   (B%XTHEAT_TARGET(:) - B%XT_SIZE_MIN(:))                     &
                 + B%XV_VENT(:) * T%XBLD_HEIGHT(:) / 3600* ZRHOA(:) * XCPD *   &
                   (B%XTHEAT_TARGET(:) - B%XT_SIZE_MIN(:))
!
!   Rated air flow rate [kg s-1 m-2(bld)]
B%XM_SYS_RAT(:) = B%XCAP_SYS_HEAT(:)/XCPD/(323.15 - B%XTHEAT_TARGET(:))
!
!    Initialization
B%XCAP_SYS_RAT  = 0.
!
!*      B.     Autosize of the cooling system
!              -----------------------------------
!
ZRHOA = 1.15
!
!-------------------------------------------------------
!
!    Initial values
!
! Options TOP
!
YBEM = TOP%CBEM
TOP%CBEM = "BEM" 
!
YZ0H = TOP%CZ0H
TOP%CZ0H = 'KAND07'
!
GCANOPY = TOP%LCANOPY
TOP%LCANOPY = .FALSE.
!
YCH_BEM = TOP%CCH_BEM
TOP%CCH_BEM = 'DOE-2'
!
!
! Teb Fields T
!
ALLOCATE(T%XWS_ROOF(KI))
ALLOCATE(T%XWS_ROAD(KI))
T%XWS_ROOF = 0.0
T%XWS_ROAD = 0.0
!
ALLOCATE(T%TSNOW_ROOF%WSNOW(KI,1))
ALLOCATE(T%TSNOW_ROOF%T(KI,1))
ALLOCATE(T%TSNOW_ROOF%ALB(KI))
ALLOCATE(T%TSNOW_ROOF%RHO(KI,1))
ALLOCATE(T%TSNOW_ROOF%TS(KI))
ALLOCATE(T%TSNOW_ROOF%EMIS(KI))
!
YSNOW_ROOF = T%TSNOW_ROOF%SCHEME 
T%TSNOW_ROOF%SCHEME = 'NONE'
T%TSNOW_ROOF%WSNOW  = 0.0
T%TSNOW_ROOF%T      = 273.0
T%TSNOW_ROOF%RHO    = 0.0
T%TSNOW_ROOF%ALB    = 0.8
T%TSNOW_ROOF%TS     = 273.0
T%TSNOW_ROOF%EMIS   = 0.0
!
ALLOCATE(T%TSNOW_ROAD%WSNOW(KI,1))
ALLOCATE(T%TSNOW_ROAD%T(KI,1))
ALLOCATE(T%TSNOW_ROAD%ALB(KI))
ALLOCATE(T%TSNOW_ROAD%RHO(KI,1))
ALLOCATE(T%TSNOW_ROAD%TS(KI))
ALLOCATE(T%TSNOW_ROAD%EMIS(KI))
!
YSNOW_ROAD = T%TSNOW_ROAD%SCHEME 
T%TSNOW_ROAD%SCHEME = 'NONE'
T%TSNOW_ROAD%WSNOW  = 0.0
T%TSNOW_ROAD%T      = 273.0
T%TSNOW_ROAD%RHO    = 0.0
T%TSNOW_ROAD%ALB    = 0.8
T%TSNOW_ROAD%TS     = 273.0
T%TSNOW_ROAD%EMIS   = 1.0
!
ZGARDEN = T%XGARDEN
T%XGARDEN = 0.
!
!------------------------------------------------
!BEM 
!
YCOOL_COIL = BOP%CCOOL_COIL
BOP%CCOOL_COIL = 'IDEAL '
YHEAT_COIL = BOP%CHEAT_COIL
BOP%CHEAT_COIL ='IDEAL '
!
!!TI_BLD
ALLOCATE(B%XTI_BLD(KI))
!! !ZTI_BLD   = 297.16 ! indoor air temperature
DO JJ=1,KI
    B%XTI_BLD(JJ) = MAX(B%XTHEAT_TARGET(JJ),ZT_CANYON(JJ)) ! indoor air temperature
ENDDO
!
!T_FLOOR, T_MASS
ALLOCATE(B%XT_FLOOR(KI,BOP%NFLOOR_LAYER))
ALLOCATE(B%XT_MASS(KI,BOP%NFLOOR_LAYER))
!
DO JJ=1,BOP%NFLOOR_LAYER
   B%XT_FLOOR(:,JJ)  = B%XTI_BLD(:) ! building floor temperature
   B%XT_MASS(:,JJ)   = B%XTI_BLD(:) ! building mass temperature
ENDDO
!
!OUTDOOR WINDOW TEMPERATURE
ALLOCATE(B%XT_WIN1(KI))
B%XT_WIN1(:) = ZT_CANYON(:)
!! 
ALLOCATE(B%XT_WIN2(KI))
B%XT_WIN2(:) = B%XTI_BLD(:)
!
YNATVENT(:) = B%CNATVENT(:)
B%CNATVENT(:) = 'NONE'
ZF_WATER_COND(:) = B%XF_WATER_COND(:)
B%XF_WATER_COND(:) = 0.
ZAUX_MAX = B%XAUX_MAX
B%XAUX_MAX = 0.
!
ALLOCATE(B%XQI_BLD(KI))
B%XQI_BLD = 0.011
!
GNATVENT_NIGHT(:) = B%LNATVENT_NIGHT(:)
B%LNATVENT_NIGHT(:) = .FALSE.
!
!
! TEB
!
ALLOCATE(T%XT_ROOF(KI,TOP%NROOF_LAYER))
ALLOCATE(T%XT_ROAD(KI,TOP%NROAD_LAYER))
ALLOCATE(T%XT_WALL_A(KI,TOP%NWALL_LAYER))
ALLOCATE(T%XT_WALL_B(KI,TOP%NWALL_LAYER))
!
!ROOF
T%XT_ROOF  (:,TOP%NROOF_LAYER)   = B%XTI_BLD(:)   ! roof layers temperatures
T%XT_ROOF(:,1) = ZT_CANYON(:)
ZT1(:) = T%XT_ROOF(:,1)
ZTN(:) = T%XT_ROOF(:,TOP%NROOF_LAYER)
IF (TOP%NROOF_LAYER .GT. 2) CALL INTERP_PROFTWALL(ZT1, ZTN, T%XD_ROOF, T%XT_ROOF)
!
!ROAD
DO JJ=1,TOP%NROAD_LAYER
   T%XT_ROAD(:,JJ) = ZT_CANYON(:)
ENDDO
!
!WALL_A
T%XT_WALL_A(:,TOP%NWALL_LAYER)   = B%XTI_BLD(:)   ! wall layers temperatures
T%XT_WALL_A(:,1) = ZT_CANYON(:)
ZT1(:)=T%XT_WALL_A(:,1)
ZTN(:)=T%XT_WALL_A(:,TOP%NWALL_LAYER)
IF (TOP%NWALL_LAYER .GT. 2) CALL INTERP_PROFTWALL(ZT1, ZTN, T%XD_WALL, T%XT_WALL_A)
!
!WALL_B
T%XT_WALL_B = T%XT_WALL_A
!
!
!* road watering (not used)
YIR%LPAR_RD_IRRIG   = .FALSE.
!
ALLOCATE(YIR%XRD_START_MONTH(KI))
ALLOCATE(YIR%XRD_END_MONTH(KI))
ALLOCATE(YIR%XRD_START_HOUR(KI))
ALLOCATE(YIR%XRD_END_HOUR(KI))
ALLOCATE(YIR%XRD_24H_IRRIG(KI))
!
YIR%XRD_START_MONTH = 1.
YIR%XRD_END_MONTH   = 1.
YIR%XRD_START_HOUR  = 0.
YIR%XRD_END_HOUR    = 24.
YIR%XRD_24H_IRRIG   = 0.
!
!----------------------------------------------
!
 CALL DIAG_INIT(YGRD)
 CALL DIAG_EVAP_ISBA_INIT(YGRDE)
!
! greenroofs   are not taken into account in the building's HVAC equipment sizing process
YGRD%XRN     (:) = 0.
YGRD%XH      (:) = 0.
YGRD%XLE     (:) = 0.
YGRD%XGFLUX  (:) = 0.
YGRDE%XRUNOFF(:) = 0.
YGRDE%XDRAIN (:) = 0.
!
 CALL DIAG_MISC_TEB_INIT(YDMT)
!
ZUW_GREENROOF    (:) = 0.
!
ZAC_GARDEN    (:) = 0.
!
ZFRAC_PANEL  = 0.
ZALB_PANEL   = 0.1
ZSVF_GARDEN   (:) = 0.
ZALB_GARDEN   (:) = 0.
ZALB_GREENROOF(:) = 0.
GSHADE        (:) = .FALSE.
GSHAD_DAY     (:) = .FALSE.
ZE_SHADING(:) = 0.
!
 CALL DIAG_MISC_TEB_OPTIONS_INIT(YDMTO)
 CALL DIAG_MISC_TEB_INIT(YDMTC)
!
 YDMTO%LSURF_MISC_BUDGET = .TRUE.
 CALL DIAG_MISC_TEB_INIT_n(YDMTC, YDMT, YDMTO, TOP, KI, 0)
!
!* one supposes zero conduction heat flux between the greenroof and the nnroof.
YDMT%XG_GREENROOF_ROOF(:) = 0.
!
YDMT%XQIN(:) = B%XQIN(:)
!
YDMT%XT_SYS(:) = B%XTI_BLD(:)
YDMT%XQ_SYS(:) = B%XQI_BLD(:)
!
YDMT%XTHEAT_TARGET(:) = B%XTHEAT_TARGET(:)
YDMT%XTCOOL_TARGET(:) = B%XTCOOL_TARGET(:)
!
INB_STEP_ATM = 3600*24*4/ZTSTEP
DO JFORC_STEP= 1,INB_STEP_ATM
!
! Daily outdoor air temperature cycle
  ZT_CANYON(:) = 10.7/2 * SIN(2*XPI/(24*3600) * (ZTIME+16*3600))  &
                 + (B%XT_SIZE_MAX(:)-10.7/2)
  ZTA(:) = ZT_CANYON(:)
!
!
!*      B.1     Solar radiation
!               ---------------
!
  CALL SUNPOS(JPYEAR, IMONTH, IDAY, ZTIME, G%XLON, G%XLAT, ZTSUN, ZZENITH, ZAZIM)
!
  CALL SW_DAYCYCLE(KI, ZZENITH, ZTOT_SW)
!
  ZDIR_SW(:) = 0.88 * ZTOT_SW(:) * 0.85 ! manual adjustment
  ZSCA_SW(:) = 0.12 * ZTOT_SW(:) * 0.85 ! manual adjustment
  WHERE (ZDIR_SW < 0.0) 
    ZDIR_SW = 0.0
  END WHERE
  WHERE (ZSCA_SW < 0.0) 
    ZSCA_SW = 0.0
  END WHERE

!
! solar panels are not taken into account in the building's HVAC equipment sizing process
  CALL URBAN_SOLAR_ABS(TOP, T, B, YDMT, ZDIR_SW, ZSCA_SW, ZZENITH, ZAZIM,&
                       ZFRAC_PANEL, ZALB_PANEL, ZALB_GARDEN, ZSVF_GARDEN,&
                       ZALB_GREENROOF, ZDN_ROOF, ZDF_ROOF, ZDN_ROAD,     &
                       ZDF_ROAD, ZREC_SW_ROAD, ZREC_SW_SNOW_ROAD,        &
                       ZREC_SW_WALL_A, ZREC_SW_WALL_B, ZREC_SW_GARDEN,   &
                       ZREC_SW_ROOF, ZDIR_ALB_TOWN, ZSCA_ALB_TOWN,       &
                       ZSW_RAD_GARDEN, ZREC_SW_WIN, ZREF_SW_GRND,        &
                       ZREF_SW_FAC, ZE_SHADING, GSHAD_DAY, GSHADE  )
!
!
!*      B.2     LW properties
!               -------------
!
  ZEMIS_GARDEN = 1.0
  CALL URBAN_LW_COEF(B, T, ZLW_RAD, ZEMIS_GARDEN,                       &
                     T%TSNOW_ROAD%TS, ZTS_GARDEN,                       &  
                     ZLW_WA_TO_WB, ZLW_WA_TO_R, ZLW_WB_TO_R,            &
                     ZLW_WA_TO_NR,ZLW_WB_TO_NR,                         &
                     ZLW_WA_TO_G, ZLW_WB_TO_G,                          &
                     ZLW_WA_TO_WIN, ZLW_WB_TO_WIN,                      &
                     ZLW_R_TO_WA, ZLW_R_TO_WB, ZLW_R_TO_WIN,            &
                     ZLW_G_TO_WA, ZLW_G_TO_WB, ZLW_G_TO_WIN,            &
                     ZLW_S_TO_WA, ZLW_S_TO_WB, ZLW_S_TO_R,              &
                     ZLW_S_TO_NR, ZLW_S_TO_G, ZLW_S_TO_WIN,             &
                     ZLW_WIN_TO_WA, ZLW_WIN_TO_WB, ZLW_WIN_TO_R,        &
                     ZLW_WIN_TO_NR, ZLW_WIN_TO_G,                       &
                     ZLW_NR_TO_WA, ZLW_NR_TO_WB, ZLW_NR_TO_WIN          )
!
!*      B.3     TEB simulation
!               -------------
!
  CALL TEB  (TOP, T, BOP, B, YIR, YDMT, YGRD, YGRDE, YIMPLICIT_WIND, ZTSUN,  &
             ZT_CANYON, ZQ_CANYON, ZU_CANYON, ZT_CANYON, ZQ_CANYON,    &
             ZU_CANYON, ZZ_LOWCAN, ZPEW_A_COEF, ZPEW_B_COEF,           &
             ZPEW_A_COEF, ZPEW_B_COEF, ZPS, ZPA, ZEXNS, ZEXNA, ZTA,    &
             ZQA, ZRHOA, ZLW_RAD, ZRR, ZSR, ZZREF, ZZREF, ZU_CANYON,   &
             T%XH_TRAFFIC, T%XLE_TRAFFIC, ZTSTEP, ZDF_ROOF, ZDN_ROOF,  &
             ZDF_ROAD, ZDN_ROAD, ZQSAT_ROOF, ZQSAT_ROAD,               &
             ZDELT_ROOF, ZDELT_ROAD, ZTS_GARDEN, ZLEW_ROOF,            &
             ZUW_GREENROOF, ZLEW_ROAD, ZLE_WALL_A, ZLE_WALL_B,         &
             ZRNSNOW_ROOF, ZHSNOW_ROOF, ZLESNOW_ROOF, ZGSNOW_ROOF,     &
             ZMELT_ROOF, ZRNSNOW_ROAD, ZHSNOW_ROAD, ZLESNOW_ROAD,      &
             ZGSNOW_ROAD, ZMELT_ROAD, ZUW_ROAD, ZUW_ROOF, ZDUWDU_ROAD, &
             ZDUWDU_ROOF, ZUSTAR_TOWN, ZCD, ZCDN, ZCH_TOWN, ZRI_TOWN,  &
             ZRESA_TOWN, ZDQS_TOWN, ZQF_TOWN,                          &
             ZQF_BLD, ZFLX_BLD, ZAC_ROOF, ZAC_ROAD, ZAC_WALL, ZAC_TOP, &
             ZAC_GARDEN, ZAC_ROOF_WAT, ZAC_ROAD_WAT, ZLW_WA_TO_WB,     &
             ZLW_WA_TO_R, ZLW_WB_TO_R, ZLW_WA_TO_NR, ZLW_WB_TO_NR,     &
             ZLW_R_TO_WA, ZLW_R_TO_WB, ZLW_G_TO_WA, ZLW_G_TO_WB,       &
             ZLW_S_TO_WA, ZLW_S_TO_WB, ZLW_S_TO_R, ZLW_S_TO_NR,        &
             ZLW_NR_TO_WA, ZLW_NR_TO_WB, ZLW_NR_TO_WIN, ZLW_WA_TO_WIN, &
             ZLW_WB_TO_WIN, ZLW_G_TO_WIN, ZLW_R_TO_WIN, ZLW_S_TO_WIN,  &
             ZLW_WIN_TO_WA, ZLW_WIN_TO_WB, ZLW_WIN_TO_R, ZLW_WIN_TO_NR,&
             IDAY, ZEMIT_LW_FAC, ZEMIT_LW_ROAD,    &
             ZT_RAD_IND, ZHU_BLD, ZTIME, ZE_SHADING  )
! 
!   Time update
  ZTIME = ZTIME + ZTSTEP
  IF (ZTIME >= 86400) THEN
    ZTIME = 0.0
    IDAY = IDAY + 1
  END IF
!
ENDDO
!
B%XQIN(:) = YDMT%XQIN(:)
B%XTHEAT_TARGET(:) = YDMT%XTHEAT_TARGET(:)
B%XTCOOL_TARGET(:) = YDMT%XTCOOL_TARGET(:)
!
 CALL TEB_IRRIG_INIT(YIR)
 CALL DIAG_MISC_TEB_INIT(YDMT)
 CALL DIAG_MISC_TEB_INIT(YDMTC)
!----------------------------------------------------
!
! Options
!
TOP%CZ0H = YZ0H
TOP%CBEM = YBEM
TOP%LCANOPY = GCANOPY
TOP%CCH_BEM = YCH_BEM
!
! Teb Fields T
!
DEALLOCATE(T%XWS_ROOF,T%XWS_ROAD)
!
T%TSNOW_ROOF%SCHEME = YSNOW_ROOF
DEALLOCATE(T%TSNOW_ROOF%WSNOW)
DEALLOCATE(T%TSNOW_ROOF%T)
DEALLOCATE(T%TSNOW_ROOF%ALB)
DEALLOCATE(T%TSNOW_ROOF%RHO)
DEALLOCATE(T%TSNOW_ROOF%TS)
DEALLOCATE(T%TSNOW_ROOF%EMIS)
!
T%TSNOW_ROAD%SCHEME = YSNOW_ROAD
DEALLOCATE(T%TSNOW_ROAD%WSNOW)
DEALLOCATE(T%TSNOW_ROAD%T)
DEALLOCATE(T%TSNOW_ROAD%ALB)
DEALLOCATE(T%TSNOW_ROAD%RHO)
DEALLOCATE(T%TSNOW_ROAD%TS)
DEALLOCATE(T%TSNOW_ROAD%EMIS)
!
T%XGARDEN = ZGARDEN 
!
!------------------------------------------------
!BEM 
!
BOP%CCOOL_COIL = YCOOL_COIL
BOP%CHEAT_COIL = YHEAT_COIL
!
!!TI_BLD
DEALLOCATE(B%XTI_BLD)
!
!T_FLOOR, T_MASS
DEALLOCATE(B%XT_FLOOR,B%XT_MASS)
!
!OUTDOOR WINDOW TEMPERATURE
DEALLOCATE(B%XT_WIN1)
DEALLOCATE(B%XT_WIN2)
!
B%CNATVENT(:) = YNATVENT(:) 
B%XF_WATER_COND(:) = ZF_WATER_COND(:)
B%XAUX_MAX = ZAUX_MAX
!
DEALLOCATE(B%XQI_BLD)
!
B%LNATVENT_NIGHT(:) = GNATVENT_NIGHT(:)
!
! TEB
!
DEALLOCATE(T%XT_ROOF,T%XT_ROAD,T%XT_WALL_A,T%XT_WALL_B)
!
!--------------------------------------------------------
!
! -----------------------------------------------------------
! Print autosize results
! -----------------------------------------------------------
!
WRITE(KLUOUT,*) ' '
WRITE(KLUOUT,*) '    --------------------------------'
WRITE(KLUOUT,*) '      HVAC AUTOSIZE CALCULATIONS '
WRITE(KLUOUT,*) ' '
WRITE(KLUOUT,*) '    Rated mass flow rate:'
WRITE(KLUOUT,*) '    ',MAXVAL(B%XM_SYS_RAT), 'kg s-1 m-2(bld)'
WRITE(KLUOUT,*) '    ',MINVAL(B%XM_SYS_RAT), 'kg s-1 m-2(bld)'
WRITE(KLUOUT,*) '    Rated cooling system capacity:'
WRITE(KLUOUT,*) '    ',MAXVAL(B%XCAP_SYS_RAT), 'W m-2(bld)'
WRITE(KLUOUT,*) '    ',MINVAL(B%XCAP_SYS_RAT), 'W m-2(bld)'
WRITE(KLUOUT,*) '    Rated heating sysem capacity:'
WRITE(KLUOUT,*) '    ',MAXVAL(B%XCAP_SYS_HEAT), 'W m-2(bld)'
WRITE(KLUOUT,*) '    ',MINVAL(B%XCAP_SYS_HEAT), 'W m-2(bld)'
WRITE(KLUOUT,*) '    --------------------------------'
WRITE(KLUOUT,*) ' '
IF (LHOOK) CALL DR_HOOK('HVAC_AUTOSIZE',1,ZHOOK_HANDLE)
!
CONTAINS 
!
SUBROUTINE INTERP_PROFTWALL(PT1, PTN, PD, PT)
!interpolation of vertical profile for 'wall' : roof/wall
!arguments
REAL, DIMENSION(:), INTENT(IN)    :: PT1 !temperature layer 1
REAL, DIMENSION(:), INTENT(IN)    :: PTN !temperature layer N
REAL, DIMENSION(:,:), INTENT(IN)  :: PD  !depth of all layers
REAL, DIMENSION(:,:), INTENT(OUT) :: PT  !temperature of all layers
!local variables
INTEGER :: ILAYER ! number of layers
REAL, DIMENSION(SIZE(PT1)) :: ZDN ! total depth from mid layer 1 to mid layer n
REAL, DIMENSION(SIZE(PT1)) :: ZD  ! sequential depth in the calculation
INTEGER :: JJ, JI

ILAYER=SIZE(PD,2)
DO JI=1,KI
   ZDN(JI) = 0.5 * PD(JI,1)
   DO JJ=2,ILAYER-1
      ZDN(JI) = ZDN(JI) + PD(JI,JJ)
   ENDDO
   ZDN(JI) = ZDN(JI) + 0.5 * PD(JI,ILAYER)
ENDDO
DO JI=1,KI
   ZD(JI) = 0.5*PD(JI,1)
   DO JJ=2,ILAYER-1
      ZD(JI) = ZD(JI) + 0.5*PD(JI,JJ)
      PT(JI,JJ) = PT1(JI) + (PTN(JI)-PT1(JI)) / ZDN(JI) * ZD(JI) 
      ZD(JI) = ZD(JI) + 0.5 * PD(JI,JJ)
   ENDDO
   PT(JI,1) = PT1(JI)
   PT(JI,ILAYER) = PTN(JI)
ENDDO
END SUBROUTINE INTERP_PROFTWALL

END SUBROUTINE HVAC_AUTOSIZE
