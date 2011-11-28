!     #########
      SUBROUTINE AVERAGED_ALBEDO_TEB(PZENITH,                      &
                       PBLD, PGARDEN, PROAD,                         &
                       PWALL_O_HOR, PCAN_HW_RATIO,                   &
                       PALB_ROOF,                                    &
                       PALB_ROAD, PSVF_ROAD,                         &
                       PALB_WALL, PSVF_WALL,                         &
                       PALB_GARDEN, PSVF_GARDEN,                     &
                       TSNOW_ROOF, TSNOW_ROAD,                       &
                       PDIR_ALB_TOWN, PSCA_ALB_TOWN                  )  
!     ###################################################
!
!!**** *AVERAGED_ALBEDO_TEB* computes averaged albedo for TEB scheme
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    01/2004
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_TYPE_SNOW
!
USE MODI_URBAN_SOLAR_ABS
USE MODE_SURF_SNOW_FRAC
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
REAL, DIMENSION(:), INTENT(IN) :: PZENITH      ! zenithal solar angle
!
REAL, DIMENSION(:), INTENT(IN) :: PBLD         ! building fraction
REAL, DIMENSION(:), INTENT(IN) :: PGARDEN      ! green area fraction
REAL, DIMENSION(:), INTENT(IN) :: PROAD        ! road fraction
REAL, DIMENSION(:), INTENT(IN) :: PWALL_O_HOR  ! vertical surf. / horizontal surf.
REAL, DIMENSION(:), INTENT(IN) :: PSVF_ROAD    ! sky-view-factor from roads
REAL, DIMENSION(:), INTENT(IN) :: PSVF_WALL    ! sky-view-factor from walls
REAL, DIMENSION(:), INTENT(IN) :: PSVF_GARDEN  ! sky-view-factor from green areas
REAL, DIMENSION(:), INTENT(IN) :: PCAN_HW_RATIO! canyon height/width ratio
!
REAL, DIMENSION(:), INTENT(IN) :: PALB_ROOF    ! roof albedo
REAL, DIMENSION(:), INTENT(IN) :: PALB_ROAD    ! road albedo
REAL, DIMENSION(:), INTENT(IN) :: PALB_WALL    ! wall albedo
REAL, DIMENSION(:), INTENT(IN) :: PALB_GARDEN  ! green areas albedo
TYPE(SURF_SNOW),    INTENT(IN) :: TSNOW_ROOF   ! snow on roofs
TYPE(SURF_SNOW),    INTENT(IN) :: TSNOW_ROAD   ! snow on roads
!
REAL, DIMENSION(:), INTENT(OUT):: PDIR_ALB_TOWN ! direct albedo
REAL, DIMENSION(:), INTENT(OUT):: PSCA_ALB_TOWN ! diffuse albedo
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
REAL, DIMENSION(SIZE(PBLD)) :: ZDN_ROOF       ! snow fraction 
REAL, DIMENSION(SIZE(PBLD)) :: ZDN_ROAD       ! on the surface
REAL, DIMENSION(SIZE(PBLD)) :: ZDF_ROOF       ! free-snow fraction 
REAL, DIMENSION(SIZE(PBLD)) :: ZDF_ROAD       ! on the surface
LOGICAL, DIMENSION(SIZE(PBLD)) :: GMASK       ! .false. (= no snow precip.)
!
!
REAL, DIMENSION(SIZE(PBLD)) :: ZDIR_SW        ! direct and diffuse shortwave radiation
REAL, DIMENSION(SIZE(PBLD)) :: ZSCA_SW        ! to mimic radiation behaviour of town
!
REAL, DIMENSION(SIZE(PBLD)) :: ZABS_SW_WALL       ! shortwave absorbed by walls
REAL, DIMENSION(SIZE(PBLD)) :: ZABS_SW_ROAD       ! shortwave absorbed by roads
REAL, DIMENSION(SIZE(PBLD)) :: ZABS_SW_ROOF       ! shortwave absorbed by roofs
REAL, DIMENSION(SIZE(PBLD)) :: ZABS_SW_GARDEN     ! shortwave absorbed by green areas
REAL, DIMENSION(SIZE(PBLD)) :: ZABS_SW_SNOW_ROAD  ! shortwave absorbed by snow
REAL, DIMENSION(SIZE(PBLD)) :: ZABS_SW_SNOW_ROOF  ! on roads, roofs,
!
REAL, DIMENSION(SIZE(PBLD)) :: ZREC_SW_ROAD       ! shortwave received by roads
REAL, DIMENSION(SIZE(PBLD)) :: ZREC_SW_WALL       ! shortwave received by walls
REAL, DIMENSION(SIZE(PBLD)) :: ZREC_SW_GARDEN     ! shortwave received by green areas
REAL, DIMENSION(SIZE(PBLD)) :: ZREC_SW_SNOW_ROAD  ! shortwave received by snow on roads
!
REAL, DIMENSION(SIZE(PBLD)) :: ZSW_RAD_GARDEN ! total solar radiation reaching green areas
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!* snow fractions
!  --------------
!
IF (LHOOK) CALL DR_HOOK('AVERAGED_ALBEDO_TEB',0,ZHOOK_HANDLE)
GMASK(:) = .FALSE.
CALL SNOW_FRAC_ROAD(TSNOW_ROAD%WSNOW(:,1,1),GMASK,ZDN_ROAD,ZDF_ROAD)
CALL SNOW_FRAC_ROOF(TSNOW_ROOF%WSNOW(:,1,1),GMASK,ZDN_ROOF,ZDF_ROOF)
!
!
!* town  direct and diffuse albedo
!  -------------------------------
!
ZDIR_SW=1.
ZSCA_SW=1.
!
CALL URBAN_SOLAR_ABS(ZDIR_SW, ZSCA_SW, PZENITH,                    &
                       PBLD, PGARDEN, PROAD,                         &
                       PWALL_O_HOR, PCAN_HW_RATIO,                   &
                       PALB_ROOF,                                    &
                       PALB_ROAD, PSVF_ROAD, PALB_WALL, PSVF_WALL,   &
                       PALB_GARDEN, PSVF_GARDEN,                     &
                       TSNOW_ROOF%ALB(:,1), TSNOW_ROAD%ALB(:,1),     &
                       ZDN_ROOF, ZDF_ROOF, ZDN_ROAD, ZDF_ROAD,       &
                       ZABS_SW_ROOF, ZABS_SW_ROAD, ZABS_SW_WALL,     &
                       ZABS_SW_GARDEN,                               &
                       ZABS_SW_SNOW_ROOF, ZABS_SW_SNOW_ROAD,         &
                       ZREC_SW_ROAD,  ZREC_SW_SNOW_ROAD,             &
                       ZREC_SW_WALL,                                 &
                       ZREC_SW_GARDEN,                               &
                       PDIR_ALB_TOWN, PSCA_ALB_TOWN,                 &
                       ZSW_RAD_GARDEN                                )  
IF (LHOOK) CALL DR_HOOK('AVERAGED_ALBEDO_TEB',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE AVERAGED_ALBEDO_TEB
