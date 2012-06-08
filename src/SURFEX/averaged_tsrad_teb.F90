!     #########
      SUBROUTINE AVERAGED_TSRAD_TEB(PEMIS_ROOF,PTS_ROOF,     &
                                      PEMIS_ROAD,PTS_ROAD,     &
                                      PEMIS_WALL,PTS_WALL,     &
                                      PEMIS_GARDEN,PTS_GARDEN, &
                                      TSNOW_ROOF,TSNOW_ROAD,   &
                                      PROAD, PGARDEN,          &
                                      PBLD,PWALL_O_HOR,        &
                                      PSVF_ROAD,PSVF_WALL,     &
                                      PSVF_GARDEN,             &
                                      PEMIS,PTSRAD             )  
!     ###################################################
!
!!**** *AVERAGED_TSRAD_TEB* computes averaged emissivity and radiative surface
!!                          temperature for TEB scheme
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
USE MODD_TYPE_SNOW
!
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_CSTS,     ONLY : XSTEFAN
!
USE MODI_URBAN_LW_COEF
!
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
REAL, DIMENSION(:), INTENT(IN) :: PEMIS_ROOF ! roof emissivity
REAL, DIMENSION(:), INTENT(IN) :: PTS_ROOF   ! roof surface temperature
REAL, DIMENSION(:), INTENT(IN) :: PEMIS_ROAD ! road emissivity
REAL, DIMENSION(:), INTENT(IN) :: PTS_ROAD   ! road surface temperature
REAL, DIMENSION(:), INTENT(IN) :: PEMIS_WALL ! wall emissivity
REAL, DIMENSION(:), INTENT(IN) :: PTS_WALL   ! wall surface temperature
REAL, DIMENSION(:), INTENT(IN) :: PEMIS_GARDEN!green area emissivity (snowfree)
REAL, DIMENSION(:), INTENT(IN) :: PTS_GARDEN ! green area surf. temp.
TYPE(SURF_SNOW),    INTENT(IN) :: TSNOW_ROOF ! snow on roofs
TYPE(SURF_SNOW),    INTENT(IN) :: TSNOW_ROAD ! snow on roads
REAL, DIMENSION(:), INTENT(IN) :: PROAD      ! road fraction
REAL, DIMENSION(:), INTENT(IN) :: PGARDEN    ! green area fraction
REAL, DIMENSION(:), INTENT(IN) :: PBLD       ! building fraction
REAL, DIMENSION(:), INTENT(IN) :: PWALL_O_HOR! vertical surf. / horizontal surf.
REAL, DIMENSION(:), INTENT(IN) :: PSVF_ROAD  ! sky-view-factor from roads
REAL, DIMENSION(:), INTENT(IN) :: PSVF_WALL  ! sky-view-factor from walls
REAL, DIMENSION(:), INTENT(IN) :: PSVF_GARDEN! sky-view-factor from green areas
REAL, DIMENSION(:), INTENT(OUT):: PEMIS      ! averaged emissivity (all tiles)
REAL, DIMENSION(:), INTENT(OUT):: PTSRAD     ! averaged radiaitve temp. (all tiles)
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZDN_ROOF       ! snow fraction 
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZDN_ROAD       ! on the surface
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZDF_ROOF       ! free-snow fraction 
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZDF_ROAD       ! on the surface
LOGICAL, DIMENSION(SIZE(PEMIS_ROOF)) :: GMASK       ! .false. (= no snow precip.)
!
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZLW_W_TO_W     ! longwave exchange coefficients
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZLW_W_TO_R     
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZLW_W_TO_G     
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZLW_W_TO_NR
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZLW_NR_TO_W
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZLW_NR_TO_R
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZLW_NR_TO_G
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZLW_NR_TO_NR
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZLW_S_TO_W
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZLW_S_TO_R
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZLW_S_TO_G
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZLW_S_TO_NR
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZLW_R_TO_W
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZLW_R_TO_R
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZLW_R_TO_G
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZLW_R_TO_NR
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZLW_G_TO_W
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZLW_G_TO_R
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZLW_G_TO_G
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZLW_G_TO_NR
!
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZLW_RAD     ! incoming LW to mimic
!                                                ! radiation behaviour of town
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZABS_LW_WALL! longwave absorbed by walls
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZABS_LW_ROAD! longwave absorbed by roads
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZABS_LW_ROOF! longwave absorbed by roofs
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZABS_LW_SNOW_ROAD! longwave absorbed by snow
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZABS_LW_SNOW_ROOF! on roads and roofs
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZABS_LW_GARDEN   ! longwave absorbed by gardens
REAL, DIMENSION(SIZE(PEMIS_ROOF)) :: ZLW_UP      ! outgoing longwave
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!* snow fractions
!  --------------
!
IF (LHOOK) CALL DR_HOOK('AVERAGED_TSRAD_TEB',0,ZHOOK_HANDLE)
GMASK(:) = .FALSE.
CALL SNOW_FRAC_ROAD(TSNOW_ROAD%WSNOW(:,1,1),GMASK,ZDN_ROAD,ZDF_ROAD)
CALL SNOW_FRAC_ROOF(TSNOW_ROOF%WSNOW(:,1,1),GMASK,ZDN_ROOF,ZDF_ROOF)
!
!* long-wave trapping coefficients
!  -------------------------------
!
CALL URBAN_LW_COEF(TSNOW_ROAD%EMIS(:,1),  PSVF_ROAD,               &
                     PEMIS_WALL,            PSVF_WALL,               &
                     PEMIS_GARDEN,          PSVF_GARDEN,             &
                     PROAD, PGARDEN,                                 &
                     ZDF_ROAD,  ZDN_ROAD,  PEMIS_ROAD,               &
                     ZLW_W_TO_W,  ZLW_W_TO_NR,  ZLW_W_TO_G,          &
                     ZLW_NR_TO_W, ZLW_NR_TO_NR, ZLW_NR_TO_G,         &
                     ZLW_G_TO_W,  ZLW_G_TO_NR,  ZLW_G_TO_G,          &
                     ZLW_S_TO_W,  ZLW_S_TO_NR,  ZLW_S_TO_G,          &
                     ZLW_R_TO_W,  ZLW_R_TO_NR,  ZLW_R_TO_G           )  

CALL URBAN_LW_COEF(PEMIS_ROAD,            PSVF_ROAD,               &
                     PEMIS_WALL,            PSVF_WALL,               &
                     PEMIS_GARDEN,          PSVF_GARDEN,             &
                     PROAD, PGARDEN,                                 &
                     ZDN_ROAD, ZDF_ROAD, TSNOW_ROAD%EMIS(:,1),       &
                     ZLW_W_TO_W,  ZLW_W_TO_R,  ZLW_W_TO_G,           &
                     ZLW_R_TO_W,  ZLW_R_TO_R,  ZLW_R_TO_G,           &
                     ZLW_G_TO_W,  ZLW_G_TO_R,  ZLW_G_TO_G,           &
                     ZLW_S_TO_W,  ZLW_S_TO_R,  ZLW_S_TO_G,           &
                     ZLW_NR_TO_W, ZLW_NR_TO_R, ZLW_NR_TO_G           )  
!
!
!
!* town averaged emissivity
!  ------------------------
!
!
PEMIS(:) =   PBLD (:) * ZDF_ROOF(:) * PEMIS_ROOF(:)                        &
           + PBLD (:) * ZDN_ROOF(:) * TSNOW_ROOF%EMIS(:,1)                 &
           + PROAD(:) * ZDF_ROAD(:)                       * ZLW_S_TO_R (:) &
           + PROAD(:) * ZDN_ROAD(:)                       * ZLW_S_TO_NR(:) &
           + PWALL_O_HOR(:)                               * ZLW_S_TO_W (:) &
           + PGARDEN(:)                                   * ZLW_S_TO_G (:)  

!
!* town radiative surface temperature
!  ----------------------------------
!
! fixed incoming LW (W/m2)
ZLW_RAD(:)= XSTEFAN * (PTS_ROAD(:) ** 4)

!
! LW absorbed by roofs
ZABS_LW_ROOF(:) = PEMIS_ROOF(:) * (ZLW_RAD(:) - XSTEFAN * PTS_ROOF(:)** 4)
!
! LW absorbed by roads
ZABS_LW_ROAD(:) =  ZLW_S_TO_R (:)*ZLW_RAD(:)                  &
                   + ZLW_R_TO_R (:)*PTS_ROAD     (:)  **4       &
                   + ZLW_W_TO_R (:)*PTS_WALL     (:)  **4       &
                   + ZLW_NR_TO_R(:)*TSNOW_ROAD%TS(:,1)**4       &
                   + ZLW_G_TO_R (:)*PTS_GARDEN   (:)  **4  
!
! LW absorbed by walls
ZABS_LW_WALL(:) =  ZLW_S_TO_W (:)*ZLW_RAD(:)                  &
                   + ZLW_W_TO_W (:)*PTS_WALL     (:)  **4       &
                   + ZLW_R_TO_W (:)*PTS_ROAD     (:)  **4       &
                   + ZLW_NR_TO_W(:)*TSNOW_ROAD%TS(:,1)**4       &
                   + ZLW_G_TO_W (:)*PTS_GARDEN   (:)  **4  
!
!* LW absorbed by snow on roof
ZABS_LW_SNOW_ROOF(:) = TSNOW_ROOF%EMIS(:,1)*ZLW_RAD(:)                 &
                       - TSNOW_ROOF%EMIS(:,1)*XSTEFAN*TSNOW_ROOF%TS(:,1)**4   
!
!* LW absorbed by snow on road
ZABS_LW_SNOW_ROAD(:) =  ZLW_S_TO_NR (:)*ZLW_RAD(:)                   &
                        + ZLW_NR_TO_NR(:)*TSNOW_ROAD%TS(:,1)**4        &
                        + ZLW_W_TO_NR (:)*PTS_WALL     (:)  **4        &
                        + ZLW_R_TO_NR (:)*PTS_ROAD     (:)  **4        &
                        + ZLW_G_TO_NR (:)*PTS_GARDEN   (:)  **4  
!
!* LW absorbed by gardens
ZABS_LW_GARDEN(:) =  ZLW_S_TO_G (:)*ZLW_RAD(:)                  &
                    + ZLW_NR_TO_G(:)*TSNOW_ROAD%TS(:,1)**4        &
                    + ZLW_W_TO_G (:)*PTS_WALL     (:)  **4        &
                    + ZLW_R_TO_G (:)*PTS_ROAD     (:)  **4        &
                    + ZLW_G_TO_G (:)*PTS_GARDEN   (:)  **4  
!
!* outgoing longwave radiation
ZLW_UP(:) = ZLW_RAD(:)                                             &
            - (     PBLD(:)    *ZDF_ROOF(:)*ZABS_LW_ROOF     (:)   &
               +    PBLD(:)    *ZDN_ROOF(:)*ZABS_LW_SNOW_ROOF(:)   &
               +    PROAD(:)   *ZDF_ROAD(:)*ZABS_LW_ROAD     (:)   &
               +    PROAD(:)   *ZDN_ROAD(:)*ZABS_LW_SNOW_ROAD(:)   &
               +PWALL_O_HOR(:)             *ZABS_LW_WALL     (:)   &
               +PGARDEN(:)                 *ZABS_LW_GARDEN   (:))  
!
!* town radiative surface temperature
!
PTSRAD(:)   = ((ZLW_UP(:) - ZLW_RAD(:)*(1.-PEMIS(:))) /PEMIS(:)/XSTEFAN)**0.25
IF (LHOOK) CALL DR_HOOK('AVERAGED_TSRAD_TEB',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE AVERAGED_TSRAD_TEB
