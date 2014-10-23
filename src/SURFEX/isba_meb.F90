!     #########
      SUBROUTINE ISBA_MEB(OMEB, OFORC_MEASURE, TPTIME, HISBA, KWG_LAYER,       &  
        PTSTEP, PPS, PZENITH, PSCA_SW, PSW_RAD, PVMOD, PRR, PSR, PRHOA,        &
        PTA, PQA, PH_VEG, PDIRCOSZW,                                           &
        PEXNS, PEXNA, PPET_A_COEF, PPET_B_COEF, PPEQ_A_COEF, PPEQ_B_COEF,      &
        PPEW_A_COEF, PPEW_B_COEF,                                              &
        PZREF, PUREF, PCH, PCD, PCDN, PRI, PRESA, PHUG, PHV, PHU, PQS,         &
        PZ0G_WITHOUT_SNOW,                                                     &
        PZ0_MEBV, PZ0H_MEBV, PZ0EFF_MEBV,                                      &
        PZ0_MEBN, PZ0H_MEBN, PZ0EFF_MEBN,                                      &
        PZ0_WITH_SNOW, PZ0H_WITH_SNOW, PZ0EFF,                                 &
        PTV, PTG, PTC, PQC, PWR, PWRVN, PWG, PWGI,                             &
        PWRMAX_CF, PRGL, PRSMIN, PGAMMA, PRS,                                  &
        PALBNIR_TVEG, PALBVIS_TVEG,PALBNIR_TSOIL, PALBVIS_TSOIL, PFALB,        &
        PSNOWALB, PSNOWALBVIS, PSNOWALBNIR, PSNOWALBFIR,                       &
        PVEG, PFF, PPSN, PPALPHAN, PZF_TALLVEG, PLAI, PROOTFRAC,               &
        PWSAT, PWFC, PWWILT,                                                   &
        PSNOWRHO, PSNOWSWE, PSNOWHEAT, PSNOWTEMP, PSNOWDZ, PEMISNOW, PFEMIS,   &
        PSWUP, PSWNET_N, PSWNET_V, PSWNET_G, PSWNET_NS, PALBT, PSWDOWN_GN,     &
        PLW_RAD, PLWNET_N, PLWNET_V, PLWNET_G, PLWDOWN_GN                      )
!     ##########################################################################
!
!                             
!!****  *ISBA_MEB*  
!!
!!    PURPOSE
!!    -------
!       Monitor for the calculation of the surface fluxes and of the
!     prognostic variables of the surface over natural areas 
!     with an explicit vegetation layer
!
!     NOTE...currently MEB can be coupled with 
!     HISBA='DIF' or '3-L' soil options
!     HSNOW='3-L' snow scheme
!     Soon, HSNOW=CRO and HPHOTO/=NON (i.e. Ags will be added)
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
!!	A. Boone           * Meteo-France *
!!      P. Samuelsson      * SMHI *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2014
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_CSTS,           ONLY : XCPD 
!
USE MODD_TYPE_DATE_SURF, ONLY : DATE_TIME
!
USE MODE_THERMOS
USE MODE_MEB,            ONLY : SNOW_INTERCEPT_EFF
!
USE MODI_SOILSTRESS
USE MODI_WET_LEAVES_FRAC
USE MODI_VEG
USE MODI_SNOW_LEAVES_FRAC_MEB
USE MODI_PREPS_FOR_MEB_EBUD_RAD
USE MODI_ISBA_SWNET_MEB
USE MODI_ISBA_LWNET_MEB
USE MODI_DRAG_MEB
USE MODI_E_BUDGET_MEB
USE MODI_ISBA_FLUXES_MEB
USE MODI_SNOW_LOAD_MEB
USE MODI_HYDRO_VEG
USE MODI_SNOW3L_ISBA

!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!              -------------------------
!
!
!* general variables
!  -----------------
!
TYPE(DATE_TIME),      INTENT(IN)    :: TPTIME        ! current date and time
!
LOGICAL,              INTENT(IN)    :: OMEB          ! True = patch with multi-energy balance 
!                                                    ! False = patch with classical ISBA 
LOGICAL,              INTENT(IN)    :: OFORC_MEASURE ! 
!
CHARACTER(LEN=*),     INTENT(IN)    :: HISBA         ! type of ISBA version:
!                                                    ! '2-L' (default)
!                                                    ! '3-L'
!                                                    ! 'DIF'
!
INTEGER, DIMENSION(:),INTENT(IN)    :: KWG_LAYER     ! Number of soil moisture layers (DIF option)
!
REAL,                 INTENT(IN)    :: PTSTEP        ! Model time step (s)
REAL, DIMENSION(:),   INTENT(IN)    :: PPS           ! Pressure [Pa]
REAL, DIMENSION(:),   INTENT(IN)    :: PZENITH       ! solar zenith angle
REAL, DIMENSION(:),   INTENT(IN)    :: PSW_RAD       ! solar (shortwave) incoming radiation [W/m2]
REAL, DIMENSION(:),   INTENT(IN)    :: PLW_RAD       ! thermal (longwave) incoming radiation [W/m2]
REAL, DIMENSION(:),   INTENT(IN)    :: PSCA_SW       ! solar diffuse incoming radiation [W/m2]
REAL, DIMENSION(:),   INTENT(IN)    :: PEXNA         ! Exner function: forcing level (-)
REAL, DIMENSION(:),   INTENT(IN)    :: PEXNS         ! Exner function: surface (-)
REAL, DIMENSION(:),   INTENT(IN)    :: PRR           ! Rain rate (kg/m2/s)
REAL, DIMENSION(:),   INTENT(IN)    :: PSR           ! Snow rate (kg/m2/s)
REAL, DIMENSION(:),   INTENT(IN)    :: PRHOA         ! air density (kg/m3)
REAL, DIMENSION(:),   INTENT(IN)    :: PVMOD         ! modulus of the wind
!                                                    ! parallel to the orography (m/s)
REAL, DIMENSION(:),   INTENT(IN)    :: PTA           ! Temperature of atmosphere (K)
REAL, DIMENSION(:),   INTENT(IN)    :: PQA           ! specific humidity of atmosphere (kg/kg)
REAL, DIMENSION(:),   INTENT(IN)    :: PH_VEG        ! height of vegetation
REAL, DIMENSION(:),   INTENT(IN)    :: PZREF         ! normal distance of the first
!                                                    ! atmospheric level to the
!                                                    ! orography (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PUREF         ! reference height of the wind (m)
!                                                    ! NOTE this is different from ZZREF
!                                                    ! ONLY in stand-alone/forced mode,
!                                                    ! NOT when coupled to a model (MesoNH)
REAL, DIMENSION(:),   INTENT(IN)    :: PDIRCOSZW     ! Director Cosinus along the z
!                                                    ! direction at the surface w-point
REAL, DIMENSION(:),   INTENT(IN)    :: PZF_TALLVEG   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLAI          ! vegetation Leaf Area Index (m2/m2)
REAL, DIMENSION(:),   INTENT(IN)    :: PVEG          ! fraction of litter on the surface
REAL, DIMENSION(:),   INTENT(IN)    :: PRGL          ! maximum solar radiation
!                                                    ! usable in photosynthesis
REAL, DIMENSION(:),   INTENT(IN)    :: PRSMIN        ! minimum stomatal resistance (s/m)
REAL, DIMENSION(:),   INTENT(IN)    :: PGAMMA        ! coefficient for the calculation
!                                                    ! of the surface stomatal resistance
REAL, DIMENSION(:),   INTENT(IN)    :: PFF           ! Floodplain fraction at the surface
REAL, DIMENSION(:),   INTENT(IN)    :: PPSN          ! fraction of the grid covered
!                                                    ! by snow
REAL, DIMENSION(:),   INTENT(IN)    :: PPALPHAN      ! snow/canopy transition coefficient
REAL, DIMENSION(:),   INTENT(IN)    :: PFALB         ! Floodplain albedo
REAL, DIMENSION(:),   INTENT(IN)    :: PALBNIR_TVEG  ! albedo of vegetation in NIR 
!                                                    ! (needed for LM_TR or MEB)
REAL, DIMENSION(:),   INTENT(IN)    :: PALBVIS_TVEG  ! albedo of vegetation in VIS 
!                                                    ! (needed for LM_TR or MEB)
REAL, DIMENSION(:),   INTENT(IN)    :: PALBNIR_TSOIL ! albedo of bare soil in NIR 
!                                                    ! (needed for LM_TR or MEB)
REAL, DIMENSION(:),   INTENT(IN)    :: PALBVIS_TSOIL ! albedo of bare soil in VIS 
REAL, DIMENSION(:),   INTENT(IN)    :: PWRMAX_CF     ! maximum vegetation interception storage (kg/m2) 
REAL, DIMENSION(:),   INTENT(IN)    :: PFEMIS        ! Floodplain emissivity (-)
REAL, DIMENSION(:,:), INTENT(IN)    :: PROOTFRAC     ! vegetation root fraction (-)
REAL, DIMENSION(:,:), INTENT(IN)    :: PWFC          ! field capacity profile               (m3/m3)
REAL, DIMENSION(:,:), INTENT(IN)    :: PWWILT        ! wilting point profile                (m3/m3)
REAL, DIMENSION(:,:), INTENT(IN)    :: PWSAT         ! porosity profile                     (m3/m3)
REAL, DIMENSION(:,:), INTENT(IN)    :: PWG, PWGI     ! PWG  = soil liquid volumetric water content (m3/m3)
!                                                    ! PWGI = soil frozen volumetric water content (m3/m3)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0G_WITHOUT_SNOW ! roughness length for momentum at snow-free canopy floor (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0_MEBV      ! roughness length for momentum over MEB vegetation part of patch (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0H_MEBV     ! roughness length for heat over MEB vegetation part of path (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0EFF_MEBV   ! roughness length for momentum over MEB vegetation part of patch (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0_MEBN      ! roughness length for momentum over MEB snow part of patch (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0H_MEBN     ! roughness length for heat over MEB snow part of path (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0EFF_MEBN   ! roughness length for momentum over MEB snow part of patch (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0_WITH_SNOW ! roughness length for momentum
!                                                    ! (with snow taken into account) (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0H_WITH_SNOW ! roughness length for heat
!                                                    ! (with snow taken into account) (m)
REAL, DIMENSION(:),   INTENT(IN)    :: PZ0EFF        ! roughness length for momentum (m)
!
! implicit atmospheric coupling coefficients:
!
REAL, DIMENSION(:),   INTENT(IN)    :: PPET_A_COEF, PPET_B_COEF, &
                                       PPEQ_A_COEF, PPEQ_B_COEF, &
                                       PPEW_A_COEF, PPEW_B_COEF  
!                                                    ! PPEW_A_COEF  A-wind coefficient
!                                                    ! PPEW_B_COEF  B-wind coefficient
!                                                    ! PPET_A_COEF  A-air temperature coefficient
!                                                    ! PPET_B_COEF  B-air temperature coefficient
!                                                    ! PPEQ_A_COEF  A-air specific humidity coefficient
!                                                    ! PPEQ_B_COEF  B-air specific humidity coefficient
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PSNOWALB      ! Snow albedo
REAL, DIMENSION(:),   INTENT(INOUT) :: PSNOWALBVIS   ! Snow VIS albedo
REAL, DIMENSION(:),   INTENT(INOUT) :: PSNOWALBNIR   ! Snow NIR albedo
REAL, DIMENSION(:),   INTENT(INOUT) :: PSNOWALBFIR   ! Snow FIR albedo
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWSWE      ! Snow model layer liquid water equivalent or 
!                                                    ! SWE (kg m-2)  
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWHEAT     ! Snow layer heat content (J/m3) 
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWRHO      ! Snow layer average density (kg/m3)
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWTEMP     ! Snow layer average temperature (K)
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSNOWDZ       ! Snow layer thickness (m)
REAL, DIMENSION(:,:), INTENT(INOUT) :: PTG           ! Soil layer average temperature (K)
REAL, DIMENSION(:),   INTENT(INOUT) :: PTV           ! Canopy vegetation temperature (K)
REAL, DIMENSION(:),   INTENT(INOUT) :: PTC           ! Canopy air temperature [K]
REAL, DIMENSION(:),   INTENT(INOUT) :: PQC           ! Canopy air specific humidity [kg/kg]
REAL, DIMENSION(:),   INTENT(INOUT) :: PWR           ! liquid water retained on the foliage
!                                                    ! of the canopy vegetation [kg/m2]
REAL, DIMENSION(:),   INTENT(INOUT) :: PWRVN         ! liquid water equiv of snow retained on the foliage
!                                                    ! of the canopy vegetation [kg/m2]
REAL, DIMENSION(:),   INTENT(INOUT) :: PRESA         ! aerodynamic resistance (s/m)
!
REAL, DIMENSION(:),   INTENT(OUT)   :: PEMISNOW      ! Snow surface emissivity (-)
REAL, DIMENSION(:),   INTENT(OUT)   :: PSWNET_N      ! net snow shortwave radiation [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PSWNET_NS     ! net snow shortwave radiation for 
!                                                    ! the *surface* snow layer 
REAL, DIMENSION(:),   INTENT(OUT)   :: PSWNET_V      ! net vegetation canopy shortwave radiation 
!                                                    ! [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PSWNET_G      ! net surface (ground) shortwave radiation [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PSWUP         ! net upwelling shortwave radiation [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PALBT         ! total surface albedo
REAL, DIMENSION(:),   INTENT(OUT)   :: PSWDOWN_GN    ! total shortwave radiation transmitted through 
                                                     ! the vegetation canopy
REAL, DIMENSION(:),   INTENT(OUT)   :: PLWNET_V      ! net vegetation canopy longwave radiation [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PLWNET_G      ! net ground longwave radiation [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PLWNET_N      ! net snow longwave radiation [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PLWDOWN_GN    ! total shortwave radiation transmitted through and emitted by the canopy
!                                                    ! reaching the snowpack/ground understory (explicit part) [W/m2]
REAL, DIMENSION(:),   INTENT(OUT)   :: PRS           ! surface stomatal resistance (s/m)
REAL, DIMENSION(:),   INTENT(OUT)   :: PCH           ! drag coefficient for heat
REAL, DIMENSION(:),   INTENT(OUT)   :: PCD           ! drag coefficient for momentum
REAL, DIMENSION(:),   INTENT(OUT)   :: PCDN          ! neutral drag coefficient for momentum
REAL, DIMENSION(:),   INTENT(OUT)   :: PRI           ! Richardson number
REAL, DIMENSION(:),   INTENT(OUT)   :: PHV           ! Halstead coefficient
REAL, DIMENSION(:),   INTENT(OUT)   :: PHU           ! grid-area humidity of the soil
REAL, DIMENSION(:),   INTENT(OUT)   :: PHUG          ! ground relative humidity
REAL, DIMENSION(:),   INTENT(OUT)   :: PQS           ! surface humidity (kg/kg)
!
!
!*      0.2    declarations of local variables
!
!
REAL, PARAMETER                                    :: ZTSTEP_EB     = 300. ! s Minimum time tstep required 
!                                                                          !   to time-split MEB energy budget
INTEGER                                            :: JTSPLIT_EB           ! number of time splits
INTEGER                                            :: JDT                  ! time split loop index
!
REAL                                               :: ZTSTEP               ! Local time split timestep (s)
REAL, DIMENSION(SIZE(PPS))                         :: ZWORK                ! Working variable [*]
REAL, DIMENSION(SIZE(PSNOWSWE,1),SIZE(PSNOWSWE,2)) :: ZSNOWCOND            ! snow thermal conductivity  [W/(m K)] 
REAL, DIMENSION(SIZE(PSNOWSWE,1),SIZE(PSNOWSWE,2)) :: ZSNOWHCAP            ! snow heat capacity [J/(m3 K)]
REAL, DIMENSION(SIZE(PPS))                         :: ZCHIP                ! 
REAL, DIMENSION(SIZE(PPS))                         :: ZTAU_N               ! 
REAL, DIMENSION(SIZE(PPS))                         :: ZALBG                ! Effective ground albedo
REAL, DIMENSION(SIZE(PPS))                         :: ZSIGMA_F             ! LW transmission factor
REAL, DIMENSION(SIZE(PPS))                         :: ZSIGMA_FN            ! LW transmission factor - including buried (snow) 
!                                                                          ! vegetation effect
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_V_DTV        ! LW Jacobian: flux derrivative d LWnet_v/dTv [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_V_DTG        ! LW Jacobian: flux derrivative d LWnet_v/dTg [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_V_DTN        ! LW Jacobian: flux derrivative d LWnet_v/dTn [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_G_DTV        ! LW Jacobian: flux derrivative d LWnet_g/dTv [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_G_DTG        ! LW Jacobian: flux derrivative d LWnet_g/dTg [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_G_DTN        ! LW Jacobian: flux derrivative d LWnet_g/dTn [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_N_DTV        ! LW Jacobian: flux derrivative d LWnet_n/dTv [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_N_DTG        ! LW Jacobian: flux derrivative d LWnet_n/dTg [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZDLWNET_N_DTN        ! LW Jacobian: flux derrivative d LWnet_n/dTn [W/(m K2)]
REAL, DIMENSION(SIZE(PPS))                         :: ZWRMAX               ! [kg/m2]
REAL, DIMENSION(SIZE(PPS))                         :: ZF2                  ! water stress coefficient for vegetation canopy (-)  
REAL, DIMENSION(SIZE(PWG,1),SIZE(PWG,2))           :: ZF2WGHT              ! water stress factor for vegetation canopy (-)  
REAL, DIMENSION(SIZE(PPS))                         :: ZF5                  ! water stress coefficient for vegetation canopy to force 
REAL, DIMENSION(SIZE(PPS))                         :: ZLAIN                ! reduced (by burying by ground-based snow) LAI (m2/m2)
REAL, DIMENSION(SIZE(PPS))                         :: ZRS                  ! stomatal resistance (s/m)
REAL, DIMENSION(SIZE(PPS))                         :: ZRSN                 ! stomatal resistance of non-snow-buried canopy (s/m)
!                                                                          ! Etv=>0 as F2=>0 (-)  
REAL, DIMENSION(SIZE(PPS))                         :: ZWRVNMAX             ! maximum snow water equivalent interception capacity (kg/m2)
REAL, DIMENSION(SIZE(PPS))                         :: ZPSNCV               ! intercepted canopy snow fraction (-) NOTE! Not the same as the
!                                                                          ! ground-based snowpack
REAL, DIMENSION(SIZE(PPS))                         :: ZMELTVN              ! intercepted canopy snow net freeze/melt rate (kg m-2 s-1)
!                                                                          ! (if it is < 0, this signifies freezing)
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMA_TA            ! linear transform energy budget coefficient for Ta
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMB_TA            ! linear transform energy budget coefficient for Ta
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMA_TC            ! linear transform energy budget coefficient for Tc
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMB_TC            ! linear transform energy budget coefficient for Tc
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMA_TN            ! linear transform energy budget coefficient for Tn
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMB_TN            ! linear transform energy budget coefficient for Tn
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMA_TG            ! linear transform energy budget coefficient for Tg
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMB_TG            ! linear transform energy budget coefficient for Tg
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMA_TV            ! linear transform energy budget coefficient for Tv
REAL, DIMENSION(SIZE(PPS))                         :: ZTHRMB_TV            ! linear transform energy budget coefficient for Tv
REAL, DIMENSION(SIZE(PPS))                         :: ZPET_A_COEF          ! atmospheric coupling coefficient: Ta
REAL, DIMENSION(SIZE(PPS))                         :: ZPET_B_COEF          ! atmospheric coupling coefficient: Ta
REAL, DIMENSION(SIZE(PPS))                         :: ZKVN                 ! snow interception efficiency
REAL, DIMENSION(SIZE(PPS))                         :: ZVELC                ! wind speed at the top of the canopy (m/s)
REAL, DIMENSION(SIZE(PPS))                         :: ZDELTA               ! fraction of the foliage
!                                                                          ! covered with intercepted water (-)
REAL, DIMENSION(SIZE(PPS))                         :: ZHUGI                ! humidity over frozen bare ground (-)
REAL, DIMENSION(SIZE(PPS))                         :: ZHVN                 ! Halstead coefficient vegetation canopy above snow (-) 
REAL, DIMENSION(SIZE(PPS))                         :: ZLEG_DELTA           ! soil evaporation delta fn (-) 
REAL, DIMENSION(SIZE(PPS))                         :: ZLEGI_DELTA          ! soil sublimation delta fn (-) 
REAL, DIMENSION(SIZE(PPS))                         :: ZHSGL                ! surface halstead cofficient for bare soil (-)
REAL, DIMENSION(SIZE(PPS))                         :: ZHSGF                ! surface halstead cofficient for bare soil ice  (-)
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_C_A            ! turb transfer coef between vegetation canopy air and atmosphere (kg/m2/s) 
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_N_A            ! ...between the snow on the ground and atmosphere    (kg/m2/s) 
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_G_C            ! ...between snow-free ground and canopy air     (kg/m2/s)    
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_N_C            ! ...between snow on the ground and canopy air   (kg/m2/s)     
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_VG_C           ! ...between vegetation canopy over snow-free ground and canopy air   (kg/m2/s) 
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_VN_C           ! ...between vegetation canopy over the snow on the ground and canopy air  (kg/m2/s)  
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_V_C            ! ...between vegetation canopy and canopy air  (kg/m2/s)               
REAL, DIMENSION(SIZE(PPS))                         :: ZFLXC_MOM            ! Effective drag coefficient for momentum [kg/(m2 s)]    
REAL, DIMENSION(SIZE(PPS))                         :: ZQSATG               ! saturation specific humidity for PTG (ground surface: kg kg-1)    
REAL, DIMENSION(SIZE(PPS))                         :: ZQSATV               ! saturation specific humidity for PTV (vegetation canopy: kg kg-1) 
REAL, DIMENSION(SIZE(PPS))                         :: ZQSATC               ! saturation specific humidity for PTC (canopy air: kg kg-1)      
REAL, DIMENSION(SIZE(PPS))                         :: ZQSATN               ! saturation specific humidity for PSNOWTEMP (snow surface: kg kg-1) 
REAL, DIMENSION(SIZE(PPS))                         :: ZDELTACVK            ! canopy interception capacity fraction including K-factor (-)  


REAL, DIMENSION(SIZE(PPS))                         :: ZRS_GV, ZLAI, ZDELTAX, ZHVG, ZWR !aabtmptest to remove

!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*      1.0    Preliminaries
!              -------------
!
IF (LHOOK) CALL DR_HOOK('ISBA_MEB',0,ZHOOK_HANDLE)
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      1.0    Preliminaries for energy and radiation budget
!              ---------------------------------------------
!
CALL PREPS_FOR_MEB_EBUD_RAD(PPS,                                     &
        PLAI,PSNOWRHO,PSNOWSWE,PSNOWHEAT,                            &
        PSNOWTEMP,PSNOWDZ,ZSNOWCOND,ZSNOWHCAP,PEMISNOW,ZTAU_N,       &
        ZSIGMA_F,ZCHIP)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      2.0    Shortwave radiative transfer
!              ----------------------------  
!
!
! NOTE, currently MEB only uses 2 of 3 potential snow albedo spectral bands
!
CALL ISBA_SWNET_MEB(PLAI,PZF_TALLVEG,                                        &
        PSNOWALBVIS,PSNOWALBNIR,                                             &
        PALBNIR_TVEG, PALBVIS_TVEG,PALBNIR_TSOIL, PALBVIS_TSOIL, PFALB,      &
        PFF,PPSN,PPALPHAN,ZTAU_N,PZENITH,PSCA_SW,PSW_RAD,                    &
        PSWUP,PSWNET_N,PSWNET_V,PSWNET_G,PSWNET_NS,                          &
        PALBT,ZALBG,PSWDOWN_GN                                               )

!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      3.0    Longwave radiative transfer
!              ---------------------------  
!
CALL ISBA_LWNET_MEB(PLAI,PPSN,PPALPHAN,                                 &
        PEMISNOW,PFEMIS,PFF,                                            &
        PTV,PTG(:,1),PSNOWTEMP(:,1),                                    &
        PLW_RAD,PLWNET_N,PLWNET_V,PLWNET_G,                             &
        ZDLWNET_V_DTV,ZDLWNET_V_DTG,ZDLWNET_V_DTN,                      &
        ZDLWNET_G_DTV,ZDLWNET_G_DTG,ZDLWNET_G_DTN,                      &
        ZDLWNET_N_DTV,ZDLWNET_N_DTG,ZDLWNET_N_DTN,                      &
        ZSIGMA_F,ZSIGMA_FN,PLWDOWN_GN                                   )
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      4.0    Fraction of leaves occupied by intercepted water
!              ------------------------------------------------
!
! Vegetation canopy:
!
! First, compute an effective veg fraction: it can only be < unity if vegetation is buried by snowpack...
!
ZWORK(:) = PLAI(:)*(1.0 - PPSN(:) + PPSN(:)*(1.0 - PPALPHAN(:))) 
! 
CALL WET_LEAVES_FRAC(PWR, ZWORK, PWRMAX_CF, PZ0_MEBV, PLAI, ZWRMAX, ZDELTA) 
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      5.0    Plant stress due to soil water deficit
!              --------------------------------------
!
! Vegetation canopy:
!
CALL SOILSTRESS(HISBA, ZF2,                      &
        PROOTFRAC, PWSAT, PWFC, PWWILT,          &
        PWG, PWGI, KWG_LAYER, ZF2WGHT, ZF5 )  
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      6.0    Plant stress, stomatal resistance and, possibly, CO2 assimilation
!              --------------------------------------------------------------------
!
! MEB-NOTE here assumed HPHOTO=='DEF' for now, Ags options to be added later 
!
!
! Canopy vegetation (no snow, or snow below the main part of the canopy):
!
CALL VEG(PSW_RAD, PTC, PQC, PPS, PRGL, PLAI, PRSMIN, PGAMMA, ZF2, ZRS)
!
! Possibly snow-buried canopy vegetation:
!
ZLAIN(:) = PLAI(:)*(1.0-PPALPHAN(:)) ! LAI of possibly buried vegetation canopy
!
CALL VEG(PSW_RAD, PTC, PQC, PPS, PRGL, ZLAIN, PRSMIN, PGAMMA, ZF2, ZRSN)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      6.0    Canopy snow (intercepted) needed diagnostics:
!              ---------------------------------------------
!
CALL SNOW_LEAVES_FRAC_MEB(PPSN,PPALPHAN,PWRVN,PTV,ZCHIP,PLAI,        &
                             ZWRVNMAX,ZPSNCV,ZMELTVN)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!*      7.0    Aerodynamic drag and heat/mass transfer/fluxes 
!              and energy budget solution
!              ----------------------------------------------
!
! NOTE, this assumes thermodynamic variable herein is potential T

ZPET_A_COEF(:)  = -PPET_A_COEF(:)*XCPD 
ZPET_B_COEF(:)  =  PPET_B_COEF(:)*XCPD
ZTHRMA_TA(:)    =  XCPD/PEXNA(:)
ZTHRMB_TA(:)    =  0.0
ZTHRMA_TC(:)    =  XCPD/PEXNS(:)
ZTHRMB_TC(:)    =  0.0
ZTHRMA_TN(:)    =  XCPD/PEXNS(:)
ZTHRMB_TN(:)    =  0.0
ZTHRMA_TG(:)    =  XCPD/PEXNS(:)
ZTHRMB_TG(:)    =  0.0
ZTHRMA_TV(:)    =  XCPD/PEXNS(:)
ZTHRMB_TV(:)    =  0.0

! understory canopy resistance already computed, set here

ZRS_GV(:)       = 5000. !aabtmptest to remove
ZLAI(:)         = 0.01  !aabtmptest to remove
ZDELTAX(:)      = 0.0   !aabtmptest to remove
ZHVG(:)         = 0.0   !aabtmptest to remove
ZWR(:)          = 0.0   !aabtmptest to remove

! Possibly split time step if large: 
! Although the energy budget is fully implicit, a very small canopy heat capacity 
! (and neglect of canopy air space heat capacity) can possibly lead to
! numerical shocks, especially during transition periods between stable and unstable 
! regimes. Thus, for relatively large steps, a simple time split scheme is activated.
! Note that soil moisture is held constant, while turbulent exchange coefficients are updated during the split.
! Also, experience shows that splitting at least once for moderately sized time steps 
! is quite effective in removing any lingering small but possible oscillations.
! Finally, for *very* small time steps (such as those for high res runs), no split performed.

JTSPLIT_EB      = 1 + INT(PTSTEP/ZTSTEP_EB)  ! number of split-time steps
ZTSTEP          = PTSTEP/JTSPLIT_EB          ! split time step...for relatively small time steps, no split

LOOP_TIME_SPLIT_EB: DO JDT=1,JTSPLIT_EB

!*      7a.    Aerodynamic drag and heat transfer coefficients
!              -----------------------------------------------

!   CALL DRAG_MEB(OFORC_MEASURE,                                         &
!              PTG(:,1), PTC, PTV, PSNOWTEMP(:,1), PTA, PQC, PQA, PVMOD, &
!              PWG(:,1), PWGI(:,1), PWSAT(:,1), PWFC(:,1),               &
!              PEXNS, PEXNA, PPS,                                        &
!              PRR, PSR, PRHOA, PZ0G_WITHOUT_SNOW,                       &
!              PZ0_MEBV, PZ0H_MEBV, PZ0EFF_MEBV,                         &
!              PZ0_MEBN, PZ0H_MEBN, PZ0EFF_MEBN,                         &
!              PZ0_WITH_SNOW, PZ0H_WITH_SNOW, PZ0EFF,                    &
!              PSNOWSWE(:,1),                                            &
!              PWR, ZWR, ZCHIP, ZTSTEP, ZRS_GV, ZRS, ZRSN,               &
!              PPSN, PPALPHAN, PZREF, PUREF, PH_VEG, PDIRCOSZW,          &
!              ZPSNCV, ZDELTA, ZDELTAX, ZLAI, PLAI,                      &
!              PCH, PCD, PCDN, PRI, PRESA, ZVELC,                        &
!              PHUG, ZHUGI, PHV, ZHVN, ZHVG, PHU, PQS, PRS,              &
!              ZLEG_DELTA, ZLEGI_DELTA, ZHSGL, ZHSGF,                    &
!              ZFLXC_C_A, ZFLXC_N_A, ZFLXC_G_C, ZFLXC_N_C,               &    
!              ZFLXC_VG_C, ZFLXC_VN_C, ZFLXC_MOM,                        &
!              ZQSATG, ZQSATV, ZQSATC, ZQSATN, ZDELTACVK                 )
!
   ZKVN(:) = SNOW_INTERCEPT_EFF(ZCHIP,ZVELC,ZWRVNMAX)


ENDDO LOOP_TIME_SPLIT_EB
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
IF (LHOOK) CALL DR_HOOK('ISBA_MEB',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE ISBA_MEB
