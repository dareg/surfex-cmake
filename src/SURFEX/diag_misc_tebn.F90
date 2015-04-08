!     #########
       SUBROUTINE DIAG_MISC_TEB_n(PTSTEP, PDQS_TOWN,PQF_BLD,PQF_TOWN, PFLX_BLD,             &
                                    PRUNOFF_TOWN,                                           &
                                    PRN_ROAD, PH_ROAD, PLE_ROAD, PGFLUX_ROAD,               &
                                    PRUNOFF_ROAD, PIRRIG_ROAD,                              &
                                    PRN_WALL_A, PH_WALL_A, PGFLUX_WALL_A,                   &
                                    PRN_WALL_B, PH_WALL_B, PGFLUX_WALL_B,                   &
                                    PRN_ROOF, PH_ROOF, PLE_ROOF, PGFLUX_ROOF,               &
                                    PRUNOFF_ROOF,                                           &
                                    PRN_STRLROOF, PH_STRLROOF,                              &
                                    PLE_STRLROOF, PGFLUX_STRLROOF,                          &
                                    PRUNOFF_STRLROOF,                                       &
                                    PRN_GREENROOF, PH_GREENROOF,                            &
                                    PLE_GREENROOF, PGFLUX_GREENROOF, PG_GREENROOF_ROOF,     &
                                    PRUNOFF_GREENROOF, PDRAIN_GREENROOF, PIRRIG_GREENROOF,  &
                                    PRN_GARDEN,PH_GARDEN,PLE_GARDEN,PGFLUX_GARDEN,          &
                                    PRUNOFF_GARDEN, PDRAIN_GARDEN, PIRRIG_GARDEN,           &
                                    PRN_BLT,PH_BLT,PLE_BLT,PGFLUX_BLT,                      &
                                    PABS_SW_ROOF,PABS_LW_ROOF,                              &
                                    PABS_SW_SNOW_ROOF,PABS_LW_SNOW_ROOF,                    &
                                    PABS_SW_ROAD,PABS_LW_ROAD,                              &
                                    PABS_SW_SNOW_ROAD,PABS_LW_SNOW_ROAD,                    &
                                    PABS_SW_WALL_A, PABS_LW_WALL_A,                         &
                                    PABS_SW_WALL_B, PABS_LW_WALL_B,                         &
                                    PABS_SW_GARDEN,PABS_LW_GARDEN,                          &  
                                    PABS_SW_GREENROOF,PABS_LW_GREENROOF,                    &  
                                    PH_BLD_COOL, PT_BLD_COOL,                               &     
                                    PH_BLD_HEAT, PLE_BLD_COOL, PLE_BLD_HEAT,                &
                                    PH_WASTE, PLE_WASTE, PHVAC_COOL,                        &
                                    PHVAC_HEAT, PCAP_SYS, PM_SYS, PCOP,                     &
                                    PQ_SYS, PT_SYS, PTR_SW_WIN, PFAN_POWER,                 &
                                    PABS_SW_WIN, PABS_LW_WIN,                               &
                                    PTCOOL_TARGET, PTHEAT_TARGET, PQIN,                     &
                                    PABS_SW_PANEL, PABS_LW_PANEL, PRN_PANEL,                &
                                    PH_PANEL, PTHER_PROD_PANEL, PPHOT_PROD_PANEL,           &
                                    PPROD_PANEL, PTHER_PROD_BLD, PPHOT_PROD_BLD             )  
!     ###############################################################################
!
!!****  *DIAG_MISC-TEB_n * - additional diagnostics for TEB
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
!!     P. Le Moigne 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2005
!!------------------------------------------------------------------
!
!
!
USE MODD_DIAG_MISC_TEB_OPTION_n, ONLY : DGMTO => DIAG_MISC_TEB_OPTIONS
USE MODD_DIAG_MISC_TEB_n, ONLY : DGMT => DIAG_MISC_TEB
!
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS 
USE MODI_CUMUL_DIAG_TEB_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
       REAL,               INTENT(IN) :: PTSTEP            ! time step
       REAL, DIMENSION(:), INTENT(IN) :: PQF_BLD           ! domestic heating
       REAL, DIMENSION(:), INTENT(IN) :: PFLX_BLD          ! heat flux from bld
       REAL, DIMENSION(:), INTENT(IN) :: PQF_TOWN          ! total anthropogenic heat
       REAL, DIMENSION(:), INTENT(IN) :: PDQS_TOWN         ! storage inside town mat.
       REAL, DIMENSION(:), INTENT(IN) :: PRN_ROAD          ! net radiation for roads
       REAL, DIMENSION(:), INTENT(IN) :: PH_ROAD           ! sensible heat flux for roads
       REAL, DIMENSION(:), INTENT(IN) :: PLE_ROAD          ! latent heat flux for roads
       REAL, DIMENSION(:), INTENT(IN) :: PGFLUX_ROAD       ! storage flux for roads
       REAL, DIMENSION(:), INTENT(IN) :: PRUNOFF_ROAD      ! runoff for roads       
       REAL, DIMENSION(:), INTENT(IN) :: PIRRIG_ROAD       ! water supply for watering of roads       
       REAL, DIMENSION(:), INTENT(IN) :: PRN_WALL_A        ! net radiation for wall
       REAL, DIMENSION(:), INTENT(IN) :: PH_WALL_A         ! sensible heat flux for walls
       REAL, DIMENSION(:), INTENT(IN) :: PGFLUX_WALL_A     ! storage flux for walls
       REAL, DIMENSION(:), INTENT(IN) :: PRN_WALL_B        ! net radiation for wall
       REAL, DIMENSION(:), INTENT(IN) :: PH_WALL_B         ! sensible heat flux for walls
       REAL, DIMENSION(:), INTENT(IN) :: PGFLUX_WALL_B     ! storage flux for walls
       REAL, DIMENSION(:), INTENT(IN) :: PRN_ROOF          ! net radiation for roofs
       REAL, DIMENSION(:), INTENT(IN) :: PH_ROOF           ! sensible heat flux for roofs
       REAL, DIMENSION(:), INTENT(IN) :: PLE_ROOF          ! latent heat flux for roofs
       REAL, DIMENSION(:), INTENT(IN) :: PGFLUX_ROOF       ! storage flux for roofs       
       REAL, DIMENSION(:), INTENT(IN) :: PRUNOFF_ROOF      ! aggregated runoff for roof
       REAL, DIMENSION(:), INTENT(IN) :: PRUNOFF_TOWN      ! aggregated runoff for town
       REAL, DIMENSION(:), INTENT(IN) :: PRN_STRLROOF      ! net radiation for structural roofs
       REAL, DIMENSION(:), INTENT(IN) :: PH_STRLROOF       ! sensible heat flux for structural roofs
       REAL, DIMENSION(:), INTENT(IN) :: PLE_STRLROOF      ! latent heat flux for structural roofs
       REAL, DIMENSION(:), INTENT(IN) :: PGFLUX_STRLROOF   ! storage flux for structural roofs       
       REAL, DIMENSION(:), INTENT(IN) :: PRUNOFF_STRLROOF  ! runoff for structural roofs       
       REAL, DIMENSION(:), INTENT(IN) :: PRN_GREENROOF     ! net radiation for green roofs
       REAL, DIMENSION(:), INTENT(IN) :: PH_GREENROOF      ! sensible heat flux for green roofs
       REAL, DIMENSION(:), INTENT(IN) :: PLE_GREENROOF     ! latent heat flux for green roofs
       REAL, DIMENSION(:), INTENT(IN) :: PGFLUX_GREENROOF  ! storage flux for green roofs
       REAL, DIMENSION(:), INTENT(IN) :: PG_GREENROOF_ROOF ! heat flux between green/structural roofs
       REAL, DIMENSION(:), INTENT(IN) :: PRUNOFF_GREENROOF ! runoff for green roofs       
       REAL, DIMENSION(:), INTENT(IN) :: PDRAIN_GREENROOF  ! total vertical drainage for green roofs       
       REAL, DIMENSION(:), INTENT(IN) :: PIRRIG_GREENROOF  ! water supply from green roof ground irrigation
       REAL, DIMENSION(:), INTENT(IN) :: PRN_GARDEN        ! net radiation for GARDEN areas
       REAL, DIMENSION(:), INTENT(IN) :: PH_GARDEN         ! sensible heat flux for GARDEN areas
       REAL, DIMENSION(:), INTENT(IN) :: PLE_GARDEN        ! latent heat flux for GARDEN areas
       REAL, DIMENSION(:), INTENT(IN) :: PGFLUX_GARDEN     ! storage flux for GARDEN areas
       REAL, DIMENSION(:), INTENT(IN) :: PRUNOFF_GARDEN    ! runoff for green areas
       REAL, DIMENSION(:), INTENT(IN) :: PDRAIN_GARDEN     ! drainage for green areas
       REAL, DIMENSION(:), INTENT(IN) :: PIRRIG_GARDEN     ! water supply for irrigation for green areas
       REAL, DIMENSION(:), INTENT(IN) :: PRN_BLT           ! net radiation for built surf
       REAL, DIMENSION(:), INTENT(IN) :: PH_BLT            ! sensible heat flux for built surf
       REAL, DIMENSION(:), INTENT(IN) :: PLE_BLT           ! latent heat flux for built surf
       REAL, DIMENSION(:), INTENT(IN) :: PGFLUX_BLT        ! storage flux for built surf
!
       REAL, DIMENSION(:), INTENT(IN) :: PABS_SW_ROOF      ! Sdown absorbed by roofs
       REAL, DIMENSION(:), INTENT(IN) :: PABS_SW_SNOW_ROOF ! Sdown absorbed by snow on roofs
       REAL, DIMENSION(:), INTENT(IN) :: PABS_LW_ROOF      ! Ldown absorbed by roofs
       REAL, DIMENSION(:), INTENT(IN) :: PABS_LW_SNOW_ROOF ! Ldown absorbed by snow on roofs
       REAL, DIMENSION(:), INTENT(IN) :: PABS_SW_ROAD      ! Sdown absorbed by roads
       REAL, DIMENSION(:), INTENT(IN) :: PABS_SW_SNOW_ROAD ! Sdown absorbed by snow on roads
       REAL, DIMENSION(:), INTENT(IN) :: PABS_LW_ROAD      ! Ldown absorbed by roads
       REAL, DIMENSION(:), INTENT(IN) :: PABS_LW_SNOW_ROAD ! Ldown absorbed by snow on roads
       REAL, DIMENSION(:), INTENT(IN) :: PABS_SW_WALL_A    ! Sdown absorbed by walls
       REAL, DIMENSION(:), INTENT(IN) :: PABS_LW_WALL_A    ! Ldown absorbed by walls
       REAL, DIMENSION(:), INTENT(IN) :: PABS_SW_WALL_B    ! Sdown absorbed by walls
       REAL, DIMENSION(:), INTENT(IN) :: PABS_LW_WALL_B    ! Ldown absorbed by walls
       REAL, DIMENSION(:), INTENT(IN) :: PABS_SW_GARDEN    ! Sdown absorbed by GARDEN areas
       REAL, DIMENSION(:), INTENT(IN) :: PABS_LW_GARDEN    ! Ldown absorbed by GARDEN areas
       REAL, DIMENSION(:), INTENT(IN) :: PABS_SW_GREENROOF ! Sdown absorbed by green roofs
       REAL, DIMENSION(:), INTENT(IN) :: PABS_LW_GREENROOF ! Ldown absorbed by green roofs
!  new arguments after BEM
       REAL, DIMENSION(:), INTENT(IN) :: PH_BLD_COOL       ! Sensible cooling energy demand  
                                                           ! of the building [W m-2(bld)]
       REAL, DIMENSION(:), INTENT(IN) :: PT_BLD_COOL       ! Total cooling energy demand  
                                                           ! of the building [W m-2(bld)]
       REAL, DIMENSION(:), INTENT(IN) :: PH_BLD_HEAT       ! Heating energy demand       
                                                           ! of the building [W m-2(bld)]
       REAL, DIMENSION(:), INTENT(IN) :: PLE_BLD_COOL      ! Latent cooling energy demand 
                                                           ! of the building [W m-2(bld)]
       REAL, DIMENSION(:), INTENT(IN) :: PLE_BLD_HEAT      ! Latent heating energy demand 
                                                           ! of the building [W m-2(bld)]
       REAL, DIMENSION(:), INTENT(IN) :: PH_WASTE          ! Sensible waste heat from HVAC system
                                                           ! [W m-2(bld)]
       REAL, DIMENSION(:), INTENT(IN) :: PLE_WASTE         ! Latent waste heat from HVAC system
                                                           ! [W m-2(bld)]
       REAL, DIMENSION(:), INTENT(IN) :: PHVAC_COOL        ! Energy consumption of the cooling system
                                                           ! [W m-2(bld)]
       REAL, DIMENSION(:), INTENT(IN) :: PHVAC_HEAT        ! Energy consumption of the heating system
                                                           ! [W m-2(bld)]
       REAL, DIMENSION(:), INTENT(IN) :: PCAP_SYS          ! Actual capacity of the cooling system
                                                           ! [W m-2(bld)] 
       REAL, DIMENSION(:), INTENT(IN) :: PM_SYS            ! Actual HVAC mass flow rate 
                                                           ! [kg s-1 m-2(bld)]
       REAL, DIMENSION(:), INTENT(IN) :: PCOP              ! COP of the cooling system
       REAL, DIMENSION(:), INTENT(IN) :: PQ_SYS            ! Supply air specific humidity [kg kg-1]
       REAL, DIMENSION(:), INTENT(IN) :: PT_SYS            ! Supply air temperature [K]
       REAL, DIMENSION(:), INTENT(IN) :: PTR_SW_WIN        ! Solar radiation transmitted throught
                                                           ! windows [W m-2(bld)]
       REAL, DIMENSION(:), INTENT(IN) :: PFAN_POWER        ! HVAC fan power
       REAL, DIMENSION(:), INTENT(IN) :: PABS_SW_WIN       ! window absorbed shortwave radiation [W m-2] 
       REAL, DIMENSION(:), INTENT(IN) :: PABS_LW_WIN       ! absorbed infrared rad. [W m-2]
       !
       REAL, DIMENSION(:), INTENT(IN) :: PTCOOL_TARGET     ! Cooling system set point modulated by bld_occ_calendar [K]
       REAL, DIMENSION(:), INTENT(IN) :: PTHEAT_TARGET     ! Heating system set point modulated by bld_occ_calendar [K] 
       REAL, DIMENSION(:), INTENT(IN) :: PQIN              ! Building interal heat load modulated by bld_occ_calendar [W m-2(floor)]
!
       REAL, DIMENSION(:), INTENT(IN) :: PABS_SW_PANEL     ! absorbed solar    energy by solar panels [W m-2(panel)]
       REAL, DIMENSION(:), INTENT(IN) :: PABS_LW_PANEL     ! absorbed longwave energy by solar panels [W m-2(panel)]
       REAL, DIMENSION(:), INTENT(IN) :: PRN_PANEL         ! net radiation            of solar panels [W m-2(panel)]
       REAL, DIMENSION(:), INTENT(IN) :: PH_PANEL          ! sensible heat flux       of solar panels [W m-2(panel)]
       REAL, DIMENSION(:), INTENT(IN) :: PTHER_PROD_PANEL  ! thermal      production  of solar panels [W m-2(panel)]
       REAL, DIMENSION(:), INTENT(IN) :: PPHOT_PROD_PANEL  ! photovoltaic production  of solar panels [W m-2(panel)]
       REAL, DIMENSION(:), INTENT(IN) :: PPROD_PANEL       ! averaged     production  of solar panels [W m-2(panel)]
       REAL, DIMENSION(:), INTENT(IN) :: PTHER_PROD_BLD    ! thermal      production  of solar panels [W m-2(bld)]
       REAL, DIMENSION(:), INTENT(IN) :: PPHOT_PROD_BLD    ! photovoltaic production  of solar panels [W m-2(bld)]
!
!*      0.2    declarations of local variables
!
       REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DIAG_MISC_TEB_N',0,ZHOOK_HANDLE)
IF (DGMTO%LSURF_MISC_BUDGET) THEN
   DGMT%XQF_BLD            =  PQF_BLD
   DGMT%XFLX_BLD           =  PFLX_BLD
   DGMT%XQF_TOWN           =  PQF_TOWN
   DGMT%XDQS_TOWN          =  PDQS_TOWN
   DGMT%XRN_ROAD           = PRN_ROAD
   DGMT%XH_ROAD            = PH_ROAD
   DGMT%XLE_ROAD           = PLE_ROAD
   DGMT%XGFLUX_ROAD        = PGFLUX_ROAD
   DGMT%XRN_WALL_A         = PRN_WALL_A
   DGMT%XH_WALL_A          = PH_WALL_A
   DGMT%XGFLUX_WALL_A      = PGFLUX_WALL_A
   DGMT%XRN_WALL_B         = PRN_WALL_B
   DGMT%XH_WALL_B          = PH_WALL_B
   DGMT%XGFLUX_WALL_B      = PGFLUX_WALL_B
   DGMT%XRN_ROOF           = PRN_ROOF
   DGMT%XH_ROOF            = PH_ROOF
   DGMT%XLE_ROOF           = PLE_ROOF
   DGMT%XGFLUX_ROOF        = PGFLUX_ROOF   
   DGMT%XRN_STRLROOF       = PRN_STRLROOF
   DGMT%XH_STRLROOF        = PH_STRLROOF
   DGMT%XLE_STRLROOF       = PLE_STRLROOF
   DGMT%XGFLUX_STRLROOF    = PGFLUX_STRLROOF
   DGMT%XRN_GREENROOF      = PRN_GREENROOF
   DGMT%XH_GREENROOF       = PH_GREENROOF
   DGMT%XLE_GREENROOF      = PLE_GREENROOF
   DGMT%XGFLUX_GREENROOF   = PGFLUX_GREENROOF
   DGMT%XG_GREENROOF_ROOF  = PG_GREENROOF_ROOF
   DGMT%XRUNOFF_TOWN       = PRUNOFF_TOWN
   DGMT%XRUNOFF_GARDEN     = PRUNOFF_GARDEN
   DGMT%XRUNOFF_ROAD       = PRUNOFF_ROAD
   DGMT%XRUNOFF_ROOF       = PRUNOFF_ROOF
   DGMT%XRUNOFF_STRLROOF   = PRUNOFF_STRLROOF
   DGMT%XRUNOFF_GREENROOF  = PRUNOFF_GREENROOF
   DGMT%XDRAIN_GARDEN      = PDRAIN_GARDEN
   DGMT%XDRAIN_GREENROOF   = PDRAIN_GREENROOF
   DGMT%XIRRIG_ROAD        = PIRRIG_ROAD
   DGMT%XIRRIG_GARDEN      = PIRRIG_GARDEN
   DGMT%XIRRIG_GREENROOF   = PIRRIG_GREENROOF
   DGMT%XRN_GARDEN         = PRN_GARDEN
   DGMT%XH_GARDEN          = PH_GARDEN
   DGMT%XLE_GARDEN         = PLE_GARDEN
   DGMT%XGFLUX_GARDEN      = PGFLUX_GARDEN  
   DGMT%XRN_BLT            = PRN_BLT  
   DGMT%XH_BLT             = PH_BLT  
   DGMT%XLE_BLT            = PLE_BLT  
   DGMT%XGFLUX_BLT         = PGFLUX_BLT    
!
   DGMT%XABS_SW_ROOF       = PABS_SW_ROOF
   DGMT%XABS_LW_ROOF       = PABS_LW_ROOF
   DGMT%XABS_SW_SNOW_ROOF  = PABS_SW_SNOW_ROOF
   DGMT%XABS_LW_SNOW_ROOF  = PABS_LW_SNOW_ROOF
   DGMT%XABS_SW_ROAD       = PABS_SW_ROAD
   DGMT%XABS_LW_ROAD       = PABS_LW_ROAD
   DGMT%XABS_SW_SNOW_ROAD  = PABS_SW_SNOW_ROAD
   DGMT%XABS_LW_SNOW_ROAD  = PABS_LW_SNOW_ROAD
   DGMT%XABS_SW_WALL_A     = PABS_SW_WALL_A
   DGMT%XABS_LW_WALL_A     = PABS_LW_WALL_A
   DGMT%XABS_SW_WALL_B     = PABS_SW_WALL_B
   DGMT%XABS_LW_WALL_B     = PABS_LW_WALL_B
   DGMT%XABS_SW_GARDEN     = PABS_SW_GARDEN
   DGMT%XABS_LW_GARDEN     = PABS_LW_GARDEN
   DGMT%XABS_SW_GREENROOF  = PABS_SW_GREENROOF
   DGMT%XABS_LW_GREENROOF  = PABS_LW_GREENROOF
   !
   IF (TOP%CBEM=='BEM') THEN
     DGMT%XH_BLD_COOL = PH_BLD_COOL 
     DGMT%XT_BLD_COOL = PT_BLD_COOL  
     DGMT%XH_BLD_HEAT = PH_BLD_HEAT  
     DGMT%XLE_BLD_COOL= PLE_BLD_COOL  
     DGMT%XLE_BLD_HEAT= PLE_BLD_HEAT 
     DGMT%XH_WASTE    = PH_WASTE      
     DGMT%XLE_WASTE   = PLE_WASTE     
     DGMT%XHVAC_COOL  = PHVAC_COOL    
     DGMT%XHVAC_HEAT  = PHVAC_HEAT     
     DGMT%XCAP_SYS    = PCAP_SYS        
     DGMT%XM_SYS      = PM_SYS         
     DGMT%XCOP        = PCOP          
     DGMT%XQ_SYS      = PQ_SYS     
     DGMT%XT_SYS      = PT_SYS  
     DGMT%XTR_SW_WIN  = PTR_SW_WIN
     DGMT%XFAN_POWER  = PFAN_POWER 
     !
     DGMT%XABS_SW_WIN = PABS_SW_WIN 
     DGMT%XABS_LW_WIN = PABS_LW_WIN
     !
     DGMT%XTCOOL_CUR_TARGET  = PTCOOL_TARGET    
     DGMT%XTHEAT_CUR_TARGET  = PTHEAT_TARGET    
     DGMT%XCUR_QIN           = PQIN    
   ENDIF
   !
   IF (TOP%LSOLAR_PANEL) THEN
     DGMT%XABS_SW_PANEL    = PABS_SW_PANEL
     DGMT%XABS_LW_PANEL    = PABS_LW_PANEL
     DGMT%XRN_PANEL        = PRN_PANEL
     DGMT%XH_PANEL         = PH_PANEL
     DGMT%XTHER_PROD_PANEL = PTHER_PROD_PANEL
     DGMT%XPHOT_PROD_PANEL = PPHOT_PROD_PANEL
     DGMT%XPROD_PANEL      = PPROD_PANEL
     DGMT%XTHER_PROD_BLD   = PTHER_PROD_BLD
     DGMT%XPHOT_PROD_BLD   = PPHOT_PROD_BLD
   END IF
   !
   ! cumulated diagnostics 
   ! ---------------------
   !
   CALL CUMUL_DIAG_TEB_n(PTSTEP)
   !
END IF
!
IF (LHOOK) CALL DR_HOOK('DIAG_MISC_TEB_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_MISC_TEB_n
