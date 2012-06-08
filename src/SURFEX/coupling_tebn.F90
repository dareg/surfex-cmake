!     ###############################################################################
SUBROUTINE COUPLING_TEB_n(HPROGRAM, HCOUPLING,                                             &
               PTSTEP, KYEAR, KMONTH, KDAY, PTIME, KI, KSV, KSW, PTSUN, PZENITH, PAZIM,    &
               PZREF, PUREF, PZS, PU, PV, PQA, PTA, PRHOA, PSV, PCO2, HSV,                 &
               PRAIN, PSNOW, PLW, PDIR_SW, PSCA_SW, PSW_BANDS, PPS, PPA,                   &
               PSFTQ, PSFTH, PSFTS, PSFCO2, PSFU, PSFV,                                    &
               PTRAD, PDIR_ALB, PSCA_ALB, PEMIS,                                           &
               PPEW_A_COEF, PPEW_B_COEF,                                                   &
               PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,                         &
               HTEST                                                                       )
!     ###############################################################################
!
!!****  *COUPLING_TEB_n * - Driver for TEB 
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     V. Masson 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!                  10/2005 (G.Pigeon) transfer of domestic heating
!!      S. Riette   06/2009 Initialisation of XT, XQ, XU and XTKE on canopy levels
!!      S. Riette   01/2010 Use of interpol_sbl to compute 10m wind diagnostic
!!---------------------------------------------------------------
!
!
USE MODD_CSTS,         ONLY : XRD, XCPD, XP00, XLVTT, XPI
USE MODD_SURF_PAR,     ONLY : XUNDEF
!
USE MODD_TEB_n,        ONLY : TTIME,LCANOPY,CZ0H,                                      &
                              XT_CANYON, XQ_CANYON, XTI_BLD,                           &
                              XT_ROOF, XT_ROAD, XT_WALL, XWS_ROOF, XWS_ROAD,           &
                              TSNOW_ROOF, TSNOW_ROAD,                                  &
                              XH_TRAFFIC, XLE_TRAFFIC, XH_INDUSTRY, XLE_INDUSTRY,      &
                              XZ0_TOWN, XBLD, XGARDEN, XROAD,                          &
                              XBLD_HEIGHT, XWALL_O_HOR, XCAN_HW_RATIO,                 &
                              XALB_ROOF, XEMIS_ROOF, XHC_ROOF,XTC_ROOF, XD_ROOF,       &
                              XALB_ROAD, XEMIS_ROAD, XHC_ROAD,XTC_ROAD, XD_ROAD,       &
                              XALB_WALL, XEMIS_WALL, XHC_WALL,XTC_WALL, XD_WALL,       &
                              XSVF_ROAD, XSVF_WALL,                                    &
                              XSVF_GARDEN,                                             &
                              XQSAT_ROOF, XQSAT_ROAD, XDELT_ROOF, XDELT_ROAD                       
USE MODD_CH_TEB_n,     ONLY : CSV, CCH_DRY_DEP, XDEP, NBEQ, NSV_CHSBEG, NSV_CHSEND,    &
                              NSV_DSTBEG, NSV_DSTEND, NAEREQ, NDSTEQ, NSLTEQ,          &
                              NSV_AERBEG, NSV_AEREND, NSV_SLTBEG, NSV_SLTEND
USE MODD_TEB_CANOPY_n, ONLY : XZ, XU, NLVL, XTKE, XT, XQ, XLMO, XLM, XLEPS, XZF, XDZ, XDZF, XP
USE MODD_DIAG_TEB_n,   ONLY : N2M, XT2M, XQ2M, XHU2M, XZON10M, XMER10M
USE MODD_DST_n,        ONLY : XEMISRADIUS_DST, XEMISSIG_DST
USE MODD_SLT_n,        ONLY : XEMISRADIUS_SLT, XEMISSIG_SLT
USE MODD_DST_SURF
USE MODD_SLT_SURF
!
USE MODE_DSLT_SURF
USE MODE_THERMOS
USE MODE_SBLS
USE MODE_COUPLING_CANOPY
!
USE MODI_SM10
USE MODI_ADD_FORECAST_TO_DATE_SURF
USE MODI_DIAG_INLINE_TEB_n
USE MODI_DIAG_MISC_TEB_n
USE MODI_CH_AER_DEP
USE MODI_CH_DEP_TOWN
USE MODI_DSLT_DEP
USE MODI_TEB_GARDEN
USE MODI_TEB_CANOPY
! 
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
USE MODI_CANOPY_EVOL
USE MODI_CANOPY_GRID_UPDATE
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=6),    INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
CHARACTER(LEN=1),    INTENT(IN)  :: HCOUPLING ! type of coupling
                                              ! 'E' : explicit
                                              ! 'I' : implicit
INTEGER,             INTENT(IN)  :: KYEAR     ! current year (UTC)
INTEGER,             INTENT(IN)  :: KMONTH    ! current month (UTC)
INTEGER,             INTENT(IN)  :: KDAY      ! current day (UTC)
REAL,                INTENT(IN)  :: PTIME     ! current time since midnight (UTC, s)
INTEGER,             INTENT(IN)  :: KI        ! number of points
INTEGER,             INTENT(IN)  :: KSV       ! number of scalars
INTEGER,             INTENT(IN)  :: KSW       ! number of short-wave spectral bands
REAL, DIMENSION(KI), INTENT(IN)  :: PTSUN     ! solar time                    (s from midnight)
REAL,                INTENT(IN)  :: PTSTEP    ! atmospheric time-step                 (s)
REAL, DIMENSION(KI), INTENT(IN)  :: PZREF     ! height of T,q forcing                 (m)
REAL, DIMENSION(KI), INTENT(IN)  :: PUREF     ! height of wind forcing                (m)
!
REAL, DIMENSION(KI), INTENT(IN)  :: PTA       ! air temperature forcing               (K)
REAL, DIMENSION(KI), INTENT(IN)  :: PQA       ! air humidity forcing                  (kg/m3)
REAL, DIMENSION(KI), INTENT(IN)  :: PRHOA     ! air density                           (kg/m3)
REAL, DIMENSION(KI,KSV),INTENT(IN) :: PSV     ! scalar variables
!                                             ! chemistry:   first char. in HSV: '#'  (molecule/m3)
!                                             !
CHARACTER(LEN=6), DIMENSION(KSV),INTENT(IN):: HSV  ! name of all scalar variables
REAL, DIMENSION(KI), INTENT(IN)  :: PU        ! zonal wind                            (m/s)
REAL, DIMENSION(KI), INTENT(IN)  :: PV        ! meridian wind                         (m/s)
REAL, DIMENSION(KI,KSW),INTENT(IN) :: PDIR_SW ! direct  solar radiation (on horizontal surf.)
!                                             !                                       (W/m2)
REAL, DIMENSION(KI,KSW),INTENT(IN) :: PSCA_SW ! diffuse solar radiation (on horizontal surf.)
!                                             !                                       (W/m2)
REAL, DIMENSION(KSW),INTENT(IN)  :: PSW_BANDS ! mean wavelength of each shortwave band (m)
REAL, DIMENSION(KI), INTENT(IN)  :: PZENITH   ! zenithal angle       (radian from the vertical)
REAL, DIMENSION(KI), INTENT(IN)  :: PAZIM     ! azimuthal angle      (radian from North, clockwise)
REAL, DIMENSION(KI), INTENT(IN)  :: PLW       ! longwave radiation (on horizontal surf.)
!                                             !                                       (W/m2)
REAL, DIMENSION(KI), INTENT(IN)  :: PPS       ! pressure at atmospheric model surface (Pa)
REAL, DIMENSION(KI), INTENT(IN)  :: PPA       ! pressure at forcing level             (Pa)
REAL, DIMENSION(KI), INTENT(IN)  :: PZS       ! atmospheric model orography           (m)
REAL, DIMENSION(KI), INTENT(IN)  :: PCO2      ! CO2 concentration in the air          (kg/m3)
REAL, DIMENSION(KI), INTENT(IN)  :: PSNOW     ! snow precipitation                    (kg/m2/s)
REAL, DIMENSION(KI), INTENT(IN)  :: PRAIN     ! liquid precipitation                  (kg/m2/s)
!
!
REAL, DIMENSION(KI), INTENT(OUT) :: PSFTH     ! flux of heat                          (W/m2)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFTQ     ! flux of water vapor                   (kg/m2/s)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFU      ! zonal momentum flux                   (Pa)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFV      ! meridian momentum flux                (Pa)
REAL, DIMENSION(KI), INTENT(OUT) :: PSFCO2    ! flux of CO2                           (kg/m2/s)
REAL, DIMENSION(KI,KSV),INTENT(OUT):: PSFTS   ! flux of scalar var.                   (kg/m2/s)
!
REAL, DIMENSION(KI), INTENT(OUT) :: PTRAD     ! radiative temperature                 (K)
REAL, DIMENSION(KI,KSW),INTENT(OUT):: PDIR_ALB! direct albedo for each spectral band  (-)
REAL, DIMENSION(KI,KSW),INTENT(OUT):: PSCA_ALB! diffuse albedo for each spectral band (-)
REAL, DIMENSION(KI), INTENT(OUT) :: PEMIS     ! emissivity                            (-)
!
REAL, DIMENSION(KI), INTENT(IN) :: PPEW_A_COEF! implicit coefficients
REAL, DIMENSION(KI), INTENT(IN) :: PPEW_B_COEF! needed if HCOUPLING='I'
REAL, DIMENSION(KI), INTENT(IN) :: PPET_A_COEF
REAL, DIMENSION(KI), INTENT(IN) :: PPEQ_A_COEF
REAL, DIMENSION(KI), INTENT(IN) :: PPET_B_COEF
REAL, DIMENSION(KI), INTENT(IN) :: PPEQ_B_COEF
CHARACTER(LEN=2),    INTENT(IN) :: HTEST ! must be equal to 'OK'
!
!
!*      0.2    declarations of local variables
!
INTEGER                     :: JSWB        ! loop counter on shortwave spectral bands
!         
REAL, DIMENSION(SIZE(PTA))  :: ZQA         ! specific humidity                 (kg/kg)
REAL, DIMENSION(SIZE(PTA))  :: ZEXNA       ! Exner function at forcing level
REAL, DIMENSION(SIZE(PTA))  :: ZEXNS       ! Exner function at surface level
REAL, DIMENSION(SIZE(PTA))  :: ZWIND       ! wind
!
! Ouput Diagnostics:
!
REAL, DIMENSION(SIZE(PTA))  :: ZU_CANYON   ! wind in canyon
REAL, DIMENSION(SIZE(PTA))  :: ZT_CANYON   ! temperature in canyon
REAL, DIMENSION(SIZE(PTA))  :: ZQ_CANYON   ! specific humidity in canyon
!
REAL, DIMENSION(SIZE(PTA))  :: ZRN_ROOF    ! net radiation on roof
REAL, DIMENSION(SIZE(PTA))  :: ZH_ROOF     ! sensible heat flux on roof
REAL, DIMENSION(SIZE(PTA))  :: ZLE_ROOF    ! latent heat flux on roof
REAL, DIMENSION(SIZE(PTA))  :: ZLEW_ROOF   ! latent heat flux on snowfree roof
REAL, DIMENSION(SIZE(PTA))  :: ZGFLUX_ROOF ! storage flux in roof
REAL, DIMENSION(SIZE(PTA))  :: ZRUNOFF_ROOF! water runoff from roof
REAL, DIMENSION(SIZE(PTA))  :: ZRN_ROAD    ! net radiation on road
REAL, DIMENSION(SIZE(PTA))  :: ZH_ROAD     ! sensible heat flux on road
REAL, DIMENSION(SIZE(PTA))  :: ZLE_ROAD    ! latent heat flux on road
REAL, DIMENSION(SIZE(PTA))  :: ZLEW_ROAD   ! latent heat flux on snowfree road
REAL, DIMENSION(SIZE(PTA))  :: ZGFLUX_ROAD ! storage flux in road
REAL, DIMENSION(SIZE(PTA))  :: ZRUNOFF_ROAD! water runoff from road
REAL, DIMENSION(SIZE(PTA))  :: ZRN_WALL    ! net radiation on walls
REAL, DIMENSION(SIZE(PTA))  :: ZH_WALL     ! sensible heat flux on walls
REAL, DIMENSION(SIZE(PTA))  :: ZLE_WALL    ! latent heat flux on walls
REAL, DIMENSION(SIZE(PTA))  :: ZGFLUX_WALL ! storage flux in walls
REAL, DIMENSION(SIZE(PTA))  :: ZRN_GARDEN  ! net radiation on green areas
REAL, DIMENSION(SIZE(PTA))  :: ZH_GARDEN   ! sensible heat flux on green areas
REAL, DIMENSION(SIZE(PTA))  :: ZLE_GARDEN  ! latent heat flux on green areas
REAL, DIMENSION(SIZE(PTA))  :: ZGFLUX_GARDEN!storage flux in green areas
REAL, DIMENSION(SIZE(PTA))  :: ZRN_BLT     ! net radiation on built surf 
REAL, DIMENSION(SIZE(PTA))  :: ZH_BLT      ! sensible heat flux on built surf 
REAL, DIMENSION(SIZE(PTA))  :: ZLE_BLT     ! latent heat flux on built surf 
REAL, DIMENSION(SIZE(PTA))  :: ZGFLUX_BLT  ! storage flux in built surf 
REAL, DIMENSION(SIZE(PTA))  :: ZRN_GRND    ! net radiation on ground built surf
REAL, DIMENSION(SIZE(PTA))  :: ZH_GRND     ! sensible heat flux on ground built surf
REAL, DIMENSION(SIZE(PTA))  :: ZLE_GRND    ! latent heat flux on ground built surf
REAL, DIMENSION(SIZE(PTA))  :: ZGFLUX_GRND ! storage flux in ground built surf
REAL, DIMENSION(SIZE(PTA))  :: ZRNSNOW_ROOF  ! net radiation over snow
REAL, DIMENSION(SIZE(PTA))  :: ZHSNOW_ROOF   ! sensible heat flux over snow
REAL, DIMENSION(SIZE(PTA))  :: ZLESNOW_ROOF  ! latent heat flux over snow
REAL, DIMENSION(SIZE(PTA))  :: ZGSNOW_ROOF   ! flux under the snow
REAL, DIMENSION(SIZE(PTA))  :: ZMELT_ROOF    ! snow melt
REAL, DIMENSION(SIZE(PTA))  :: ZRNSNOW_ROAD  ! net radiation over snow
REAL, DIMENSION(SIZE(PTA))  :: ZHSNOW_ROAD   ! sensible heat flux over snow
REAL, DIMENSION(SIZE(PTA))  :: ZLESNOW_ROAD  ! latent heat flux over snow
REAL, DIMENSION(SIZE(PTA))  :: ZGSNOW_ROAD   ! flux under the snow
REAL, DIMENSION(SIZE(PTA))  :: ZMELT_ROAD    ! snow melt
!
REAL, DIMENSION(SIZE(PTA))  :: ZRN           ! net radiation over town
REAL, DIMENSION(SIZE(PTA))  :: ZH            ! sensible heat flux over town
REAL, DIMENSION(SIZE(PTA))  :: ZLE           ! latent heat flux over town
REAL, DIMENSION(SIZE(PTA))  :: ZGFLUX        ! flux through the ground
REAL, DIMENSION(SIZE(PTA))  :: ZQF_BLD       ! domestic heating
REAL, DIMENSION(SIZE(PTA))  :: ZQF_BLDWFR    ! domestic heating
REAL, DIMENSION(SIZE(PTA))  :: ZFLX_BLD      ! flux from bld
REAL, DIMENSION(SIZE(PTA))  :: ZTI_BLD_EQ    ! internal temperature witout heating
REAL, DIMENSION(SIZE(PTA))  :: ZTI_BLDWFR    ! internal temperature witout F/R
REAL, DIMENSION(SIZE(PTA))  :: ZDQS_TOWN     ! storage inside town materials
REAL, DIMENSION(SIZE(PTA))  :: ZQF_TOWN      ! total anthropogenic heat
REAL, DIMENSION(SIZE(PTA))  :: ZEVAP         ! evaporation (km/m2/s)
REAL, DIMENSION(SIZE(PTA))  :: ZRUNOFF       ! runoff over the ground
REAL, DIMENSION(SIZE(PTA))  :: ZCD           ! drag coefficient
REAL, DIMENSION(SIZE(PTA))  :: ZCDN          ! neutral drag coefficient
REAL, DIMENSION(SIZE(PTA))  :: ZCH           ! heat drag
REAL, DIMENSION(SIZE(PTA))  :: ZRI           ! Richardson number
REAL, DIMENSION(SIZE(PTA))  :: ZUW_GRND      ! momentum flux for ground built surf
REAL, DIMENSION(SIZE(PTA))  :: ZUW_ROOF      ! momentum flux for roofs
REAL, DIMENSION(SIZE(PTA))  :: ZDUWDU_GRND   !
REAL, DIMENSION(SIZE(PTA))  :: ZDUWDU_ROOF   !
REAL, DIMENSION(SIZE(PTA))  :: ZUSTAR        ! friction velocity
!
REAL, DIMENSION(SIZE(PTA))  :: ZDIR_ALB      ! direct albedo of town
REAL, DIMENSION(SIZE(PTA))  :: ZSCA_ALB      ! diffuse albedo of town
!
REAL, DIMENSION(SIZE(PTA))  :: ZH_TRAFFIC    ! anthropogenic sensible
!                                            ! heat fluxes due to traffic
REAL, DIMENSION(SIZE(PTA))  :: ZLE_TRAFFIC   ! anthropogenic latent
!                                            ! heat fluxes due to traffic
REAL, DIMENSION(SIZE(PTA))  :: ZRESA_TOWN    ! aerodynamical resistance
REAL, DIMENSION(SIZE(PTA))  :: ZAC_ROAD      ! road aerodynamical conductance
REAL, DIMENSION(SIZE(PTA))  :: ZAC_GARDEN    ! green area aerodynamical conductance
REAL, DIMENSION(SIZE(PTA))  :: ZAC_GRND      ! ground built surf aerodynamical conductance
REAL, DIMENSION(SIZE(PTA))  :: ZAC_ROAD_WAT  ! road water aerodynamical conductance
REAL, DIMENSION(SIZE(PTA))  :: ZAC_GARDEN_WAT! green area water aerodynamical conductance
REAL, DIMENSION(SIZE(PTA))  :: ZAC_GRND_WAT  ! ground built surf water aerodynamical conductance
!
REAL                        :: ZBEGIN_TRAFFIC_TIME ! start traffic time (solar time, s)
REAL                        :: ZEND_TRAFFIC_TIME   ! end traffic time   (solar time, s)
REAL, DIMENSION(SIZE(PTA))  :: ZDIR_SW       ! total direct SW
REAL, DIMENSION(SIZE(PTA))  :: ZSCA_SW       ! total diffuse SW
REAL, DIMENSION(SIZE(PTA))  :: ZPEW_A_COEF   ! implicit coefficients
REAL, DIMENSION(SIZE(PTA))  :: ZPEW_B_COEF   ! needed if HCOUPLING='I'
!
!***** CANOPY  *****
REAL, DIMENSION(SIZE(PTA))        :: ZSFLUX_U  ! Surface flux u'w' (m2/s2)
REAL, DIMENSION(SIZE(PTA))        :: ZSFLUX_T  ! Surface flux w'T' (mK/s)
REAL, DIMENSION(SIZE(PTA))        :: ZSFLUX_Q  ! Surface flux w'q' (kgm2/s)
REAL, DIMENSION(SIZE(PTA),NLVL)   :: ZFORC_U   ! tendency due to drag force for wind
REAL, DIMENSION(SIZE(PTA),NLVL)   :: ZDFORC_UDU! formal derivative of
!                                              ! tendency due to drag force for wind
REAL, DIMENSION(SIZE(PTA),NLVL)   :: ZFORC_E   ! tendency due to drag force for TKE
REAL, DIMENSION(SIZE(PTA),NLVL)   :: ZDFORC_EDE! formal derivative of
!                                              ! tendency due to drag force for TKE
REAL, DIMENSION(SIZE(PTA),NLVL)   :: ZFORC_T   ! tendency due to drag force for Temp
REAL, DIMENSION(SIZE(PTA),NLVL)   :: ZDFORC_TDT! formal derivative of
!                                              ! tendency due to drag force for Temp
REAL, DIMENSION(SIZE(PTA),NLVL)   :: ZFORC_Q   ! tendency due to drag force for hum
REAL, DIMENSION(SIZE(PTA),NLVL)   :: ZDFORC_QDQ! formal derivative of
!                                              ! tendency due to drag force for hum.
REAL, DIMENSION(SIZE(PTA))        :: ZE_ROOF
REAL, DIMENSION(SIZE(PTA))        :: ZE_GRND
REAL, DIMENSION(SIZE(PTA))        :: ZT_LOWCAN  ! temperature at lowest canyon level (K)
REAL, DIMENSION(SIZE(PTA))        :: ZQ_LOWCAN  ! humidity    at lowest canyon level (kg/kg)
REAL, DIMENSION(SIZE(PTA))        :: ZU_LOWCAN  ! wind        at lowest canyon level (m/s)
REAL, DIMENSION(SIZE(PTA))        :: ZZ_LOWCAN  ! height      of lowest canyon level (m)
REAL, DIMENSION(SIZE(PTA))        :: ZPEW_A_COEF_LOWCAN   ! implicit coefficients for wind coupling
REAL, DIMENSION(SIZE(PTA))        :: ZPEW_B_COEF_LOWCAN   ! between first canopy level and road
REAL, DIMENSION(SIZE(PTA))        :: ZTA        ! temperature at canyon level just above roof (K)
REAL, DIMENSION(SIZE(PTA))        :: ZPA        ! pressure    at canyon level just above roof (K)
REAL, DIMENSION(SIZE(PTA))        :: ZUA        ! wind        at canyon level just above roof (m/s)
REAL, DIMENSION(SIZE(PTA))        :: ZUREF      ! height      of canyon level just above roof (m)
REAL, DIMENSION(SIZE(PTA))        :: ZZREF      ! height      of canyon level just above roof (m)
REAL, DIMENSION(SIZE(PTA))        :: ZLAMBDA_F  ! frontal density (-)
REAL, DIMENSION(SIZE(PTA))        :: ZLMO       ! Monin-Obukhov length at canopy height (m)
REAL, DIMENSION(SIZE(PTA),NLVL)   :: ZL         ! Mixing length generic profile at mid levels
REAL, DIMENSION(SIZE(PTA))        :: ZCOEF
!
REAL, DIMENSION(SIZE(PTA))        :: ZALFAU   ! V+(1) = alfa u'w'(1) + beta
REAL, DIMENSION(SIZE(PTA))        :: ZBETAU   ! V+(1) = alfa u'w'(1) + beta
REAL, DIMENSION(SIZE(PTA))        :: ZALFAT   ! Th+(1) = alfa w'th'(1) + beta
REAL, DIMENSION(SIZE(PTA))        :: ZBETAT   ! Th+(1) = alfa w'th'(1) + beta
REAL, DIMENSION(SIZE(PTA))        :: ZALFAQ   ! Q+(1) = alfa w'q'(1) + beta
REAL, DIMENSION(SIZE(PTA))        :: ZBETAQ   ! Q+(1) = alfa w'q'(1) + beta
!***** CANOPY  *****
REAL, DIMENSION(SIZE(PTA))        :: ZWAKE      ! reduction of average wind speed
                                                ! in canyon due to direction average.
!
! absorbed solar and infra-red radiation by road, wall and roof
!                                                      
REAL, DIMENSION(SIZE(PTA)) :: ZABS_SW_ROAD
REAL, DIMENSION(SIZE(PTA)) :: ZABS_SW_WALL
REAL, DIMENSION(SIZE(PTA)) :: ZABS_SW_ROOF
REAL, DIMENSION(SIZE(PTA)) :: ZABS_SW_GARDEN
REAL, DIMENSION(SIZE(PTA)) :: ZABS_SW_SNOW_ROAD
REAL, DIMENSION(SIZE(PTA)) :: ZABS_SW_SNOW_ROOF
REAL, DIMENSION(SIZE(PTA)) :: ZABS_LW_SNOW_ROAD
REAL, DIMENSION(SIZE(PTA)) :: ZABS_LW_SNOW_ROOF
REAL, DIMENSION(SIZE(PTA)) :: ZABS_LW_ROAD
REAL, DIMENSION(SIZE(PTA)) :: ZABS_LW_WALL
REAL, DIMENSION(SIZE(PTA)) :: ZABS_LW_ROOF
REAL, DIMENSION(SIZE(PTA)) :: ZABS_LW_GARDEN
!
REAL, DIMENSION(SIZE(PTA))  :: ZWIND10M
REAL, DIMENSION(SIZE(PTA))  :: ZWIND10M_MAX
REAL, DIMENSION(SIZE(PTA))  :: ZT2M_MIN
REAL, DIMENSION(SIZE(PTA))  :: ZT2M_MAX
REAL, DIMENSION(SIZE(PTA))  :: ZHU2M_MIN
REAL, DIMENSION(SIZE(PTA))  :: ZHU2M_MAX
!
REAL                       :: ZCONVERTFACM0_SLT, ZCONVERTFACM0_DST
REAL                       :: ZCONVERTFACM3_SLT, ZCONVERTFACM3_DST
REAL                       :: ZCONVERTFACM6_SLT, ZCONVERTFACM6_DST
!
INTEGER                           :: JI
INTEGER                           :: JJ
INTEGER                           :: JLAYER
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
! Preliminaries:
!-------------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('COUPLING_TEB_N',0,ZHOOK_HANDLE)
IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('COUPLING_TEBN: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF

!-------------------------------------------------------------------------------------
!
! scalar fluxes
!
PSFTS(:,:) = 0.
!
! broadband radiative fluxes
!
ZDIR_SW(:) = 0.
ZSCA_SW(:) = 0.
DO JSWB=1,KSW
  DO JJ=1,SIZE(PDIR_SW,1)
    ZDIR_SW(JJ) = ZDIR_SW(JJ) + PDIR_SW(JJ,JSWB)
    ZSCA_SW(JJ) = ZSCA_SW(JJ) + PSCA_SW(JJ,JSWB)
  ENDDO
END DO
!
DO JJ=1,SIZE(PTA)
! specific humidity (conversion from kg/m3 to kg/kg)
!
  ZQA(JJ) = PQA(JJ) / PRHOA(JJ)
!
! wind
!
  ZWIND(JJ) = SQRT(PU(JJ)**2+PV(JJ)**2)
!
ENDDO
! method of wind coupling
!
IF (HCOUPLING=='I') THEN
  ZPEW_A_COEF = PPEW_A_COEF
  ZPEW_B_COEF = PPEW_B_COEF
ELSE
  ZPEW_A_COEF =  0.
  ZPEW_B_COEF =  ZWIND
END IF
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Time evolution
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
TTIME%TIME = TTIME%TIME + PTSTEP
CALL ADD_FORECAST_TO_DATE_SURF(TTIME%TDATE%YEAR,TTIME%TDATE%MONTH,TTIME%TDATE%DAY,TTIME%TIME)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Anthropogenic fluxes (except building heating)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
ZBEGIN_TRAFFIC_TIME = 21600.
ZEND_TRAFFIC_TIME   = 64800.
!
WHERE(       PTSUN>ZBEGIN_TRAFFIC_TIME   &
      .AND.  PTSUN<ZEND_TRAFFIC_TIME     )
  ZH_TRAFFIC  (:) = XH_TRAFFIC   (:)
  ZLE_TRAFFIC (:) = XLE_TRAFFIC  (:)
ELSEWHERE
  ZH_TRAFFIC  (:) = 0.
  ZLE_TRAFFIC (:) = 0.   
END WHERE
!
!--------------------------------------------------------------------------------------
!  Canyon forcing for TEB
!--------------------------------------------------------------------------------------
!
IF (LCANOPY) THEN
!-------------------------------------------------------------------------------------
! Updates canopy vertical grid as a function of forcing height
!-------------------------------------------------------------------------------------
!
!* determines where is the forcing level and modifies the upper levels of the canopy grid
!
  CALL CANOPY_GRID_UPDATE(KI,NLVL,XBLD_HEIGHT,XBLD_HEIGHT+PUREF,XZ,XZF,XDZ,XDZF)
!
!* Initialisations of T, Q, TKE and wind at first time step
!

  IF(ANY(XT(:,:) == XUNDEF)) THEN
    DO JLAYER=1,NLVL
      XT(:,JLAYER) = XT_ROAD(:,1)        
    END  DO
    XQ(:,:) = 0.
    XTKE(:,:) = 1.
    XU(:,:) = 0.
  ENDIF
!
!* default forcing above roof: forcing level
ZUREF(:)     = PUREF(:)
ZZREF(:)     = PZREF(:)
ZUA(:)       = XU(:,NLVL)
ZTA(:)       = XT(:,NLVL)
ZQA(:)       = XQ(:,NLVL)/PRHOA(:)
ZPA(:)       = XP(:,NLVL)
!* for the time being, only one value is kept for wall in-canyon forcing, in the middle of the canyon
ZU_CANYON(:) = ZUA(:)
XT_CANYON(:) = ZTA(:)
XQ_CANYON(:) = ZQA(:)
  DO JLAYER=1,NLVL-1
    DO JI=1,KI
      !* finds middle canyon layer
      IF (XZ(JI,JLAYER)<XBLD_HEIGHT(JI)/2. .AND. XZ(JI,JLAYER+1)>=XBLD_HEIGHT(JI)/2.) THEN
        ZCOEF(JI) = (XBLD_HEIGHT(JI)/2.-XZ(JI,JLAYER))/(XZ(JI,JLAYER+1)-XZ(JI,JLAYER))
        ZU_CANYON(JI) = XU(JI,JLAYER) + ZCOEF(JI) * (XU(JI,JLAYER+1)-XU(JI,JLAYER))
        XT_CANYON(JI) = XT(JI,JLAYER) + ZCOEF(JI) * (XT(JI,JLAYER+1)-XT(JI,JLAYER))
        XQ_CANYON(JI) =(XQ(JI,JLAYER) + ZCOEF(JI) * (XQ(JI,JLAYER+1)-XQ(JI,JLAYER)))/PRHOA(JI)
      END IF
      !* finds layer just above roof (at least 1m above roof)
      IF (XZ(JI,JLAYER)<XBLD_HEIGHT(JI)+1. .AND. XZ(JI,JLAYER+1)>=XBLD_HEIGHT(JI)+1.) THEN
        ZUREF(JI) = XZ(JI,JLAYER+1) - XBLD_HEIGHT(JI)
        ZZREF(JI) = XZ(JI,JLAYER+1) - XBLD_HEIGHT(JI)
        ZTA  (JI) = XT(JI,JLAYER+1)
        ZQA  (JI) = XQ(JI,JLAYER+1)/PRHOA(JI)
        !ZUA  (JI) = XU(JI,JLAYER+1)
        ZUA  (JI) = MAX(XU(JI,JLAYER+1) - 2.*SQRT(XTKE(JI,JLAYER+1)) , XU(JI,JLAYER+1)/3.)
        ZPA  (JI) = XP(JI,JLAYER+1)
        ZLMO (JI) = XLMO(JI,JLAYER+1)
      END IF
    END DO
  END DO
  ZU_CANYON= MAX(ZU_CANYON,0.2)
  ZU_LOWCAN=XU(:,1)
  ZT_LOWCAN=XT(:,1)
  ZQ_LOWCAN=XQ(:,1) / PRHOA(:)
  ZZ_LOWCAN=XZ(:,1)
  WHERE(ZPA==XUNDEF) ZPA = PPA   ! security for first time step
!
!-------------------------------------------------------------------------------------
! determine the vertical profile for mixing and dissipative lengths (at full levels)
!-------------------------------------------------------------------------------------
!
! frontal density
  ZLAMBDA_F(:) = XCAN_HW_RATIO*XBLD / (0.5*XPI)
!
  CALL SM10(XZ,XBLD_HEIGHT,ZLAMBDA_F,ZL)
!
!-------------------------------------------------------------------------------------
! computes coefficients for implicitation
!-------------------------------------------------------------------------------------
!
  ZUW_GRND(:)      = 0.
  ZDUWDU_GRND(:)   = 0.
  ZUW_ROOF(:)      = 0.
  ZDUWDU_ROOF(:)   = 0.
  ZH_GRND(:)       = 0.
  ZH_WALL(:)       = 0.
  ZH_ROOF(:)       = 0.
  ZE_GRND(:)       = 0.
  ZE_ROOF(:)       = 0.
  ZAC_GRND(:)      = 0.
  ZAC_GRND_WAT(:)  = 0.
  ZSFLUX_U(:)      = 0.
  ZSFLUX_T(:)      = 0.
  ZSFLUX_Q(:)      = 0.
!
  DO JLAYER=1,NLVL-1
      !* Monin-Obuhkov theory not used inside the urban canopy
      ! => neutral mixing  if layer is below : (roof level +1 meter)
      WHERE (XZ(:,JLAYER)<=XBLD_HEIGHT(:)+1.) XLMO(:,JLAYER) = XUNDEF
  ENDDO
!
!
!* computes tendencies on wind and Tke due to canopy
  CALL  TEB_CANOPY(KI,NLVL,XZ,XZF,XDZ,XDZF,XBLD,XBLD_HEIGHT,XWALL_O_HOR,PPA,PRHOA,XU, &
                    ZDUWDU_GRND, ZUW_ROOF, ZDUWDU_ROOF,                               &
                    ZH_WALL,ZH_ROOF,ZE_ROOF,ZAC_GRND,ZAC_GRND_WAT,                    &
                    ZFORC_U,ZDFORC_UDU,ZFORC_E,ZDFORC_EDE,ZFORC_T,ZDFORC_TDT,ZFORC_Q,ZDFORC_QDQ)
!
!* computes coefficients for implicitation
  CALL CANOPY_EVOL(KI,NLVL,PTSTEP,1,                         &
                     ZL,ZWIND,PTA,PQA,PPA,PRHOA,             &
                     ZSFLUX_U,ZSFLUX_T,ZSFLUX_Q,             &
                     ZFORC_U,ZDFORC_UDU,ZFORC_E,ZDFORC_EDE,  &
                     ZFORC_T,ZDFORC_TDT,ZFORC_Q,ZDFORC_QDQ,  &
                     XZ,XZF,XDZ,XDZF,XU,XTKE,XT,XQ,XLMO,     &
                     XLM,XLEPS,XP,ZUSTAR,                    &
                     ZALFAU,ZBETAU,ZALFAT,ZBETAT,ZALFAQ,ZBETAQ)
!
  ZPEW_A_COEF_LOWCAN = - ZALFAU / PRHOA
  ZPEW_B_COEF_LOWCAN = ZBETAU  
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ELSE              ! no canopy case
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  DO JI=1,KI
!* skimming flow for h/w>1 (maximum effect of direction on wind in the canyon);
!* isolated flow for h/w<0.5 (wind is the same in large streets for all dir.)
!* wake flow between.
!
    ZWAKE(JI)= 1. + (2./XPI-1.) * 2. * (XCAN_HW_RATIO(JI)-0.5)
    ZWAKE(JI)= MAX(MIN(ZWAKE(JI),1.),2./XPI)
!
!* Estimation of canyon wind speed from wind just above roof level
!  (at 1.33h). Wind at 1.33h is estimated using the log law.
!
   IF (XBLD_HEIGHT(JI) .GT. 0.) THEN
    ZU_CANYON(JI) = ZWAKE(JI) * EXP(-XCAN_HW_RATIO(JI)/4.) * ZWIND(JI)     &
              * LOG( (           2.* XBLD_HEIGHT(JI)/3.) / XZ0_TOWN(JI))   &
              / LOG( (PUREF(JI)+ 2.* XBLD_HEIGHT(JI)/3.) / XZ0_TOWN(JI))
    ZZ_LOWCAN(JI) = XBLD_HEIGHT(JI) / 2.
   ELSE
    ZU_CANYON(JI) = ZWIND(JI)
    ZZ_LOWCAN(JI) = PZREF(JI)
   ENDIF
 END DO
!
!* Without SBL scheme, canyon air is assumed at mid height
  ZU_LOWCAN=ZU_CANYON
  ZT_LOWCAN=XT_CANYON
  ZQ_LOWCAN=XQ_CANYON
  ZUREF    =PUREF
  ZZREF    =PZREF
  ZTA      =PTA
  ZUA      =ZWIND
  ZPA      =PPA
  ZPEW_A_COEF_LOWCAN =  0.
  ZPEW_B_COEF_LOWCAN =  ZU_CANYON
END IF
!
! Exner functions
!
ZEXNS     (:) = (PPS(:)/XP00)**(XRD/XCPD)
ZEXNA     (:) = (ZPA(:)/XP00)**(XRD/XCPD)

!--------------------------------------------------------------------------------------
! Over Urban surfaces/towns:
!--------------------------------------------------------------------------------------
!
ZT_CANYON(:) = XT_CANYON(:)
ZQ_CANYON(:) = XQ_CANYON(:)
!
ZLESNOW_ROOF(:) = 0.
ZLESNOW_ROAD(:) = 0.
CALL TEB_GARDEN      (CZ0H, HCOUPLING, TTIME, ZT_CANYON, ZQ_CANYON, ZU_CANYON,         &
                      ZT_LOWCAN, ZQ_LOWCAN, ZU_LOWCAN, ZZ_LOWCAN,                      &
                      XTI_BLD,                                                         &
                      XT_ROOF, XT_ROAD, XT_WALL, XWS_ROOF,XWS_ROAD,                    &
                      TSNOW_ROOF%SCHEME,                                               &
                      TSNOW_ROOF%WSNOW(:,:,1), TSNOW_ROOF%T(:,:,1),                    &
                      TSNOW_ROOF%RHO(:,:,1), TSNOW_ROOF%ALB(:,1),                      &
                      TSNOW_ROOF%TS(:,1), TSNOW_ROOF%EMIS(:,1),                        &
                      TSNOW_ROAD%SCHEME,                                               &
                      TSNOW_ROAD%WSNOW(:,:,1), TSNOW_ROAD%T(:,:,1),                    &
                      TSNOW_ROAD%RHO(:,:,1), TSNOW_ROAD%ALB(:,1),                      &
                      TSNOW_ROAD%TS(:,1), TSNOW_ROAD%EMIS(:,1),                        &
                      ZPEW_A_COEF, ZPEW_B_COEF,                                        &
                      ZPEW_A_COEF_LOWCAN, ZPEW_B_COEF_LOWCAN,                          &
                      PPS, ZPA, ZEXNS, ZEXNA, ZTA, ZQA, PRHOA, PCO2,                   &
                      PLW, PDIR_SW, PSCA_SW, PSW_BANDS, KSW, PZENITH,                  &
                      PRAIN, PSNOW, ZZREF, ZUREF, ZUA,                                 &
                      ZH_TRAFFIC, ZLE_TRAFFIC, XH_INDUSTRY, XLE_INDUSTRY,              &
                      PTSTEP,                                                          &
                      XZ0_TOWN,                                                        &
                      XBLD, XGARDEN, XROAD,                                            &
                      XBLD_HEIGHT, XWALL_O_HOR, XCAN_HW_RATIO,                         &
                      XALB_ROOF, XEMIS_ROOF,                                           &
                      XHC_ROOF,XTC_ROOF,XD_ROOF,                                       &
                      XALB_ROAD, XEMIS_ROAD, XSVF_ROAD,                                &
                      XHC_ROAD,XTC_ROAD,XD_ROAD,                                       &
                      XALB_WALL, XEMIS_WALL, XSVF_WALL,                                &
                      XSVF_GARDEN,                                                     &
                      XHC_WALL,XTC_WALL,XD_WALL,                                       &
                      ZRN_ROOF, ZH_ROOF, ZLE_ROOF, ZLEW_ROOF, ZGFLUX_ROOF,             &
                      ZRUNOFF_ROOF,                                                    &
                      ZRN_ROAD, ZH_ROAD, ZLE_ROAD, ZLEW_ROAD, ZGFLUX_ROAD,             &
                      ZRUNOFF_ROAD,                                                    &
                      ZRN_WALL, ZH_WALL, ZLE_WALL, ZGFLUX_WALL,                        &
                      ZRN_GARDEN,ZH_GARDEN,ZLE_GARDEN,ZGFLUX_GARDEN,                   &
                      ZRN_BLT,ZH_BLT,ZLE_BLT,ZGFLUX_BLT,                               &
                      ZRNSNOW_ROOF, ZHSNOW_ROOF, ZLESNOW_ROOF, ZGSNOW_ROOF,            &
                      ZMELT_ROOF,                                                      &
                      ZRNSNOW_ROAD, ZHSNOW_ROAD, ZLESNOW_ROAD, ZGSNOW_ROAD,            &
                      ZMELT_ROAD,                                                      &
                      ZRN_GRND, ZH_GRND, ZLE_GRND, ZGFLUX_GRND,                        &
                      ZRN, ZH, ZLE, ZGFLUX, ZEVAP, ZRUNOFF, PSFCO2,                    &
                      ZUW_GRND, ZUW_ROOF, ZDUWDU_GRND, ZDUWDU_ROOF,                    &
                      ZUSTAR, ZCD, ZCDN, ZCH, ZRI,                                     &
                      PTRAD, PEMIS, ZDIR_ALB, ZSCA_ALB, ZRESA_TOWN, ZDQS_TOWN,         &
                      ZQF_TOWN, ZQF_BLD, ZQF_BLDWFR, ZTI_BLD_EQ, ZTI_BLDWFR,           &
                      ZFLX_BLD, ZAC_ROAD, ZAC_GARDEN,                                  &
                      ZAC_ROAD_WAT, ZAC_GARDEN_WAT,                                    &
                      ZABS_SW_ROOF,ZABS_LW_ROOF,                                       &
                      ZABS_SW_SNOW_ROOF,ZABS_LW_SNOW_ROOF,                             &
                      ZABS_SW_ROAD,ZABS_LW_ROAD,                                       &
                      ZABS_SW_SNOW_ROAD,ZABS_LW_SNOW_ROAD,                             &
                      ZABS_SW_WALL,ZABS_LW_WALL,                                       &
                      ZABS_SW_GARDEN,ZABS_LW_GARDEN                                    )
!
IF (.NOT. LCANOPY) THEN
   XT_CANYON(:) = ZT_CANYON(:)
   XQ_CANYON(:) = ZQ_CANYON(:)
ENDIF
!
!-------------------------------------------------------------------------------------
! Use of the canopy version of TEB
!-------------------------------------------------------------------------------------
!
IF (LCANOPY) THEN
!-------------------------------------------------------------------------------------
! diagnostic of evaporation over roofs and ground built surfaces
!-------------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------------
! diagnostic of evaporation over roofs and ground built surfaces
!-------------------------------------------------------------------------------------
!
ZE_ROOF = ZLE_ROOF / XLVTT
ZE_GRND = ZLE_GRND / XLVTT
!
!-------------------------------------------------------------------------------------
! Computes the impact of canopy and surfaces on air
!-------------------------------------------------------------------------------------
!
WHERE (XROAD(:)+XGARDEN(:).NE.0.)
        ZAC_GRND(:) = (XROAD(:)*ZAC_ROAD(:) + XGARDEN(:)*ZAC_GARDEN(:)) / (XROAD(:)+XGARDEN(:))
        ZAC_GRND_WAT(:) = (XROAD(:)*ZAC_ROAD_WAT(:) + XGARDEN(:)*ZAC_GARDEN_WAT(:)) / (XROAD(:)+XGARDEN(:))
ELSEWHERE
        ZAC_GRND(:) = 0.
        ZAC_GRND_WAT(:) = 0.
ENDWHERE
!
ZSFLUX_U (:) = ZUW_GRND(:) * (1.-XBLD(:))
ZSFLUX_T (:) = ZH_GRND (:) * (1.-XBLD(:)) / XCPD / PRHOA(:) 
ZSFLUX_Q (:) = ZE_GRND (:) * (1.-XBLD(:))
!
!
CALL TEB_CANOPY(KI,NLVL,XZ,XZF,XDZ,XDZF,XBLD,XBLD_HEIGHT,XWALL_O_HOR,PPA,PRHOA,XU,         &
                ZDUWDU_GRND, ZUW_ROOF, ZDUWDU_ROOF,                                        &
                ZH_WALL,ZH_ROOF,ZE_ROOF,ZAC_GRND,ZAC_GRND_WAT,                             &
                ZFORC_U,ZDFORC_UDU,ZFORC_E,ZDFORC_EDE,ZFORC_T,ZDFORC_TDT,ZFORC_Q,ZDFORC_QDQ)
!
!-------------------------------------------------------------------------------------
!* Evolution of canopy air due to these impacts
!-------------------------------------------------------------------------------------
!
CALL CANOPY_EVOL(KI,NLVL,PTSTEP,2,                                            &
                 ZL,ZWIND,PTA,PQA,PPA,PRHOA,                                  &
                 ZSFLUX_U,ZSFLUX_T,ZSFLUX_Q,                                  &
                 ZFORC_U,ZDFORC_UDU,ZFORC_E,ZDFORC_EDE,                       &
                 ZFORC_T,ZDFORC_TDT,ZFORC_Q,ZDFORC_QDQ,                       &
                 XZ,XZF,XDZ,XDZF,XU,XTKE,XT,XQ,XLMO,XLM,XLEPS,XP,             &
                 ZUSTAR,                                                      &
                 ZALFAU,ZBETAU,ZALFAT,ZBETAT,ZALFAQ,ZBETAQ                    )
!
!
END IF
!
!-------------------------------------------------------------------------------------
! Outputs:
!-------------------------------------------------------------------------------------
!
! Momentum fluxes
!
PSFU = 0.
PSFV = 0.
DO JJ=1,SIZE(PU)
  IF (ZWIND(JJ)>0.) THEN
    ZCOEF(JJ) = - PRHOA(JJ) * ZUSTAR(JJ)**2 / ZWIND(JJ)
    PSFU(JJ) = ZCOEF(JJ) * PU(JJ)
    PSFV(JJ) = ZCOEF(JJ) * PV(JJ)
  ENDIF
ENDDO
!
! Heat and CO2 fluxes
!
PSFTH(:)        = ZH(:)
PSFTQ(:)        = ZEVAP(:)
!
DO JSWB=1,SIZE(PSW_BANDS)
  DO JJ=1,SIZE(ZDIR_ALB)
    PDIR_ALB(JJ,JSWB) = ZDIR_ALB(JJ)
    PSCA_ALB(JJ,JSWB) = ZSCA_ALB(JJ)
  ENDDO
END DO
!
!-------------------------------------------------------------------------------------
! Scalar fluxes:
!-------------------------------------------------------------------------------------
!
IF (NBEQ>0) THEN
  IF (CCH_DRY_DEP == "WES89") THEN
    CALL CH_DEP_TOWN(ZRESA_TOWN,  ZUSTAR, PTA, PTRAD, XWALL_O_HOR,&
                     PSV(:,NSV_CHSBEG:NSV_CHSEND),        &
                     CSV(NSV_CHSBEG:NSV_CHSEND),             &
                     XDEP(:,1:NBEQ)  )
   
    DO JI=NSV_CHSBEG,NSV_CHSEND
!cdir nodep
      DO JJ=1,SIZE(PSFTS,1)
        PSFTS(JJ,JI) = - PSV(JJ,JI) * XDEP(JJ,JI-NSV_CHSBEG+1)
      ENDDO
    ENDDO

    IF (NAEREQ > 0 ) THEN
      CALL CH_AER_DEP(PSV(:,NSV_AERBEG:NSV_AEREND),&
                         PSFTS(:,NSV_AERBEG:NSV_AEREND),&
                         ZUSTAR,ZRESA_TOWN,PTA,PRHOA)   
    END IF

  ELSE
    DO JI=NSV_CHSBEG,NSV_CHSEND
      PSFTS(:,JI) =0.
    ENDDO
    IF(NSV_AERBEG.LT.NSV_AEREND) THEN
      DO JI=NSV_AERBEG,NSV_AEREND
        PSFTS(:,JI) =0.
      ENDDO
    ENDIF
  ENDIF
ENDIF

IF (NDSTEQ>0) THEN
  CALL DSLT_DEP(PSV(:,NSV_DSTBEG:NSV_DSTEND), PSFTS(:,NSV_DSTBEG:NSV_DSTEND),   &
                ZUSTAR, ZRESA_TOWN, PTA, PRHOA, XEMISSIG_DST, XEMISRADIUS_DST,  &
                JPMODE_DST, XDENSITY_DST, XMOLARWEIGHT_DST, ZCONVERTFACM0_DST,  &
                ZCONVERTFACM6_DST, ZCONVERTFACM3_DST, LVARSIG_DST, LRGFIX_DST,  &
                CVERMOD  )  

  CALL MASSFLUX2MOMENTFLUX(         &
    PSFTS(:,NSV_DSTBEG:NSV_DSTEND), & !I/O ![kg/m2/sec] In: flux of only mass, out: flux of moments
    PRHOA,                          & !I [kg/m3] air density
    XEMISRADIUS_DST,                &!I [um] emitted radius for the modes (max 3)
    XEMISSIG_DST,                   &!I [-] emitted sigma for the different modes (max 3)
    NDSTMDE,                        &
    ZCONVERTFACM0_DST,              &
    ZCONVERTFACM6_DST,              &
    ZCONVERTFACM3_DST,              &
    LVARSIG_DST, LRGFIX_DST         )  
ENDIF
IF (NSLTEQ>0) THEN
  CALL DSLT_DEP(PSV(:,NSV_SLTBEG:NSV_SLTEND), PSFTS(:,NSV_SLTBEG:NSV_SLTEND),   &
                ZUSTAR, ZRESA_TOWN, PTA, PRHOA, XEMISSIG_SLT, XEMISRADIUS_SLT,  &
                JPMODE_SLT, XDENSITY_SLT, XMOLARWEIGHT_SLT, ZCONVERTFACM0_SLT,  &
                ZCONVERTFACM6_SLT, ZCONVERTFACM3_SLT, LVARSIG_SLT, LRGFIX_SLT,  &
                CVERMOD  )  

  CALL MASSFLUX2MOMENTFLUX(         &
    PSFTS(:,NSV_SLTBEG:NSV_SLTEND), & !I/O ![kg/m2/sec] In: flux of only mass, out: flux of moments
    PRHOA,                          & !I [kg/m3] air density
    XEMISRADIUS_SLT,                &!I [um] emitted radius for the modes (max 3)
    XEMISSIG_SLT,                   &!I [-] emitted sigma for the different modes (max 3)
    NSLTMDE,                        &
    ZCONVERTFACM0_SLT,              &
    ZCONVERTFACM6_SLT,              &
    ZCONVERTFACM3_SLT,              &
    LVARSIG_SLT, LRGFIX_SLT         ) 
ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Inline diagnostics
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
CALL DIAG_INLINE_TEB_n(PTA, PTRAD, ZQA, PPA, PPS, PRHOA, PU, PV, PZREF, PUREF, &
                       ZCD, ZCDN, ZRI, ZCH, XZ0_TOWN,                          &
                       PTRAD, PEMIS, PDIR_ALB, PSCA_ALB,                       &
                       PLW, PDIR_SW, PSCA_SW,                                  &
                       PSFTH, PSFTQ, PSFU, PSFV,                               &
                       ZRN, ZH, ZLE, ZGFLUX                                    )
!
!
!-------------------------------------------------------------------------------------
!
CALL DIAG_MISC_TEB_n(PTSTEP, ZDQS_TOWN, ZQF_BLD, ZQF_BLDWFR, ZQF_TOWN,         &
                     ZFLX_BLD, ZTI_BLD_EQ, ZTI_BLDWFR,                         &
                     ZRN_ROAD, ZH_ROAD, ZLE_ROAD, ZGFLUX_ROAD,                 &
                     ZRN_WALL, ZH_WALL, ZGFLUX_WALL,                           &
                     ZRN_ROOF, ZH_ROOF, ZLE_ROOF, ZGFLUX_ROOF, ZRUNOFF,        &
                     ZRN_GARDEN,ZH_GARDEN,ZLE_GARDEN,ZGFLUX_GARDEN,            &
                     ZRN_BLT,ZH_BLT,ZLE_BLT,ZGFLUX_BLT,                        &
                     ZABS_SW_ROOF,ZABS_LW_ROOF,                                &
                     ZABS_SW_SNOW_ROOF,ZABS_LW_SNOW_ROOF,                      &
                     ZABS_SW_ROAD,ZABS_LW_ROAD,                                &
                     ZABS_SW_SNOW_ROAD,ZABS_LW_SNOW_ROAD,                      &
                     ZABS_SW_WALL,ZABS_LW_WALL,                                &
                     ZABS_SW_GARDEN,ZABS_LW_GARDEN                             )
!
!-------------------------------------------------------------------------------------
!          
IF (LCANOPY) THEN
  IF (N2M>0) CALL INIT_2M_10M( XP(:,2), XT(:,2), XQ(:,2),  XU, XZ, &
                               PU, PV, ZWIND, PRHOA,               &
                               XT2M, XQ2M, XHU2M, XZON10M, XMER10M,&
                               ZWIND10M, ZWIND10M_MAX, ZT2M_MIN,   &
                               ZT2M_MAX, ZHU2M_MIN, ZHU2M_MAX )
ENDIF
!
IF (LHOOK) CALL DR_HOOK('COUPLING_TEB_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE COUPLING_TEB_n


