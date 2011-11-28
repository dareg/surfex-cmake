!     #########
       SUBROUTINE DIAG_MISC_TEB_n(PTSTEP, PDQS_TOWN,PQF_BLD,PQF_BLDWFR,PQF_TOWN,          &
                                    PFLX_BLD,PTI_BLD_EQ,PTI_BLDWFR,                         &
                                    PRN_ROAD, PH_ROAD, PLE_ROAD, PGFLUX_ROAD,               &
                                    PRN_WALL, PH_WALL, PGFLUX_WALL,                         &
                                    PRN_ROOF, PH_ROOF, PLE_ROOF, PGFLUX_ROOF,               &
                                    PRUNOFF,                                                &
                                    PRN_GARDEN,PH_GARDEN,PLE_GARDEN,PGFLUX_GARDEN,          &
                                    PRN_BLT,PH_BLT,PLE_BLT,PGFLUX_BLT,                      &
                                    PABS_SW_ROOF,PABS_LW_ROOF,                              &
                                    PABS_SW_SNOW_ROOF,PABS_LW_SNOW_ROOF,                    &
                                    PABS_SW_ROAD,PABS_LW_ROAD,                              &
                                    PABS_SW_SNOW_ROAD,PABS_LW_SNOW_ROAD,                    &
                                    PABS_SW_WALL,PABS_LW_WALL,                              &
                                    PABS_SW_GARDEN,PABS_LW_GARDEN                           )  
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
USE MODD_DIAG_MISC_TEB_n,    ONLY : XQF_BLD, XQF_TOWN, XDQS_TOWN,                 &
                                      XFLX_BLD,XTI_BLD_EQ,XQF_BLDWFR,               &
                                      LSURF_MISC_BUDGET,LSURF_BUDGETC,              &
                                      XTI_BLDWFR,XRN_ROAD, XH_ROAD, XLE_ROAD,       &
                                      XGFLUX_ROAD, XRN_WALL, XH_WALL,               &
                                      XGFLUX_WALL, XRN_ROOF, XH_ROOF, XLE_ROOF,     &
                                      XGFLUX_ROOF, XRUNOFF, XRUNOFFC,               &
                                      XRN_GARDEN,XH_GARDEN,XLE_GARDEN,XGFLUX_GARDEN,&
                                      XRN_BLT,XH_BLT,XLE_BLT,XGFLUX_BLT,            &
                                      XABS_SW_ROOF,XABS_LW_ROOF,                    &
                                      XABS_SW_SNOW_ROOF,XABS_LW_SNOW_ROOF,          &
                                      XABS_SW_ROAD,XABS_LW_ROAD,                    &
                                      XABS_SW_SNOW_ROAD,XABS_LW_SNOW_ROAD,          &
                                      XABS_SW_WALL,XABS_LW_WALL,                    &
                                      XABS_SW_GARDEN,XABS_LW_GARDEN  
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
       REAL,               INTENT(IN) :: PTSTEP        ! time step
       REAL, DIMENSION(:), INTENT(IN) :: PQF_BLD       ! domestic heating
       REAL, DIMENSION(:), INTENT(IN) :: PQF_BLDWFR    ! domestic heating
       REAL, DIMENSION(:), INTENT(IN) :: PFLX_BLD      ! heat flux from bld
       REAL, DIMENSION(:), INTENT(IN) :: PTI_BLD_EQ    ! internal T without heating
       REAL, DIMENSION(:), INTENT(IN) :: PTI_BLDWFR    ! internal T without F/R
       REAL, DIMENSION(:), INTENT(IN) :: PQF_TOWN      ! total anthropogenic heat
       REAL, DIMENSION(:), INTENT(IN) :: PDQS_TOWN     ! storage inside town mat.
       REAL, DIMENSION(:), INTENT(IN) :: PRN_ROAD      ! net radiation for roads
       REAL, DIMENSION(:), INTENT(IN) :: PH_ROAD       ! sensible heat flux for roads
       REAL, DIMENSION(:), INTENT(IN) :: PLE_ROAD      ! latent heat flux for roads
       REAL, DIMENSION(:), INTENT(IN) :: PGFLUX_ROAD   ! storage flux for roads
       REAL, DIMENSION(:), INTENT(IN) :: PRN_WALL      ! net radiation for wall
       REAL, DIMENSION(:), INTENT(IN) :: PH_WALL       ! sensible heat flux for walls
       REAL, DIMENSION(:), INTENT(IN) :: PGFLUX_WALL   ! storage flux for walls
       REAL, DIMENSION(:), INTENT(IN) :: PRN_ROOF      ! net radiation for roofs
       REAL, DIMENSION(:), INTENT(IN) :: PH_ROOF       ! sensible heat flux for roofs
       REAL, DIMENSION(:), INTENT(IN) :: PLE_ROOF      ! latent heat flux for roofs
       REAL, DIMENSION(:), INTENT(IN) :: PGFLUX_ROOF   ! storage flux for roofs       
       REAL, DIMENSION(:), INTENT(IN) :: PRUNOFF       ! runoff for town    
       REAL, DIMENSION(:), INTENT(IN) :: PRN_GARDEN    ! net radiation for GARDEN areas
       REAL, DIMENSION(:), INTENT(IN) :: PH_GARDEN     ! sensible heat flux for GARDEN areas
       REAL, DIMENSION(:), INTENT(IN) :: PLE_GARDEN    ! latent heat flux for GARDEN areas
       REAL, DIMENSION(:), INTENT(IN) :: PGFLUX_GARDEN ! storage flux for GARDEN areas
       REAL, DIMENSION(:), INTENT(IN) :: PRN_BLT       ! net radiation for built surf
       REAL, DIMENSION(:), INTENT(IN) :: PH_BLT        ! sensible heat flux for built surf
       REAL, DIMENSION(:), INTENT(IN) :: PLE_BLT       ! latent heat flux for built surf
       REAL, DIMENSION(:), INTENT(IN) :: PGFLUX_BLT    ! storage flux for built surf
!
       REAL, DIMENSION(:), INTENT(IN) :: PABS_SW_ROOF      ! Sdown absorbed by roofs
       REAL, DIMENSION(:), INTENT(IN) :: PABS_SW_SNOW_ROOF ! Sdown absorbed by snow on roofs
       REAL, DIMENSION(:), INTENT(IN) :: PABS_LW_ROOF      ! Ldown absorbed by roofs
       REAL, DIMENSION(:), INTENT(IN) :: PABS_LW_SNOW_ROOF ! Ldown absorbed by snow on roofs
       REAL, DIMENSION(:), INTENT(IN) :: PABS_SW_ROAD      ! Sdown absorbed by roads
       REAL, DIMENSION(:), INTENT(IN) :: PABS_SW_SNOW_ROAD ! Sdown absorbed by snow on roads
       REAL, DIMENSION(:), INTENT(IN) :: PABS_LW_ROAD      ! Ldown absorbed by roads
       REAL, DIMENSION(:), INTENT(IN) :: PABS_LW_SNOW_ROAD ! Ldown absorbed by snow on roads
       REAL, DIMENSION(:), INTENT(IN) :: PABS_SW_WALL      ! Sdown absorbed by walls
       REAL, DIMENSION(:), INTENT(IN) :: PABS_LW_WALL      ! Ldown absorbed by walls
       REAL, DIMENSION(:), INTENT(IN) :: PABS_SW_GARDEN    ! Sdown absorbed by GARDEN areas
       REAL, DIMENSION(:), INTENT(IN) :: PABS_LW_GARDEN    ! Ldown absorbed by GARDEN areas
       REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*      0.2    declarations of local variables
!
!
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DIAG_MISC_TEB_N',0,ZHOOK_HANDLE)
IF (LSURF_MISC_BUDGET) THEN
   XQF_BLD            =  PQF_BLD
   XQF_BLDWFR         =  PQF_BLDWFR
   XFLX_BLD           =  PFLX_BLD
   XTI_BLD_EQ         =  PTI_BLD_EQ
   XTI_BLDWFR         =  PTI_BLDWFR
   XQF_TOWN           =  PQF_TOWN
   XDQS_TOWN          =  PDQS_TOWN
   XRN_ROAD           = PRN_ROAD
   XH_ROAD            = PH_ROAD
   XLE_ROAD           = PLE_ROAD
   XGFLUX_ROAD        = PGFLUX_ROAD
   XRN_WALL           = PRN_WALL
   XH_WALL            = PH_WALL
   XGFLUX_WALL        = PGFLUX_WALL
   XRN_ROOF           = PRN_ROOF
   XH_ROOF            = PH_ROOF
   XLE_ROOF           = PLE_ROOF
   XGFLUX_ROOF        = PGFLUX_ROOF   
   XRUNOFF            = PRUNOFF
   XRN_GARDEN         = PRN_GARDEN
   XH_GARDEN          = PH_GARDEN
   XLE_GARDEN         = PLE_GARDEN
   XGFLUX_GARDEN      = PGFLUX_GARDEN  
   XRN_BLT            = PRN_BLT  
   XH_BLT             = PH_BLT  
   XLE_BLT            = PLE_BLT  
   XGFLUX_BLT         = PGFLUX_BLT    
!
   XABS_SW_ROOF       = PABS_SW_ROOF
   XABS_LW_ROOF       = PABS_LW_ROOF
   XABS_SW_SNOW_ROOF  = PABS_SW_SNOW_ROOF
   XABS_LW_SNOW_ROOF  = PABS_LW_SNOW_ROOF
   XABS_SW_ROAD       = PABS_SW_ROAD
   XABS_LW_ROAD       = PABS_LW_ROAD
   XABS_SW_SNOW_ROAD  = PABS_SW_SNOW_ROAD
   XABS_LW_SNOW_ROAD  = PABS_LW_SNOW_ROAD
   XABS_SW_WALL       = PABS_SW_WALL
   XABS_LW_WALL       = PABS_LW_WALL
   XABS_SW_GARDEN     = PABS_SW_GARDEN
   XABS_LW_GARDEN     = PABS_LW_GARDEN
END IF
!
IF (LSURF_BUDGETC) THEN
   XRUNOFFC           = XRUNOFFC + PRUNOFF * PTSTEP
ENDIF
IF (LHOOK) CALL DR_HOOK('DIAG_MISC_TEB_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_MISC_TEB_n
