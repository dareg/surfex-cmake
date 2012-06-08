!     #########
      SUBROUTINE CONVERT_COVER_TEB   (PCOVER,                             &
                                        PZ0_TOWN,                           &
                                        PALB_ROOF,                          &
                                        PEMIS_ROOF,PHC_ROOF,PTC_ROOF,       &
                                        PD_ROOF,                            &
                                        PALB_ROAD,                          &
                                        PEMIS_ROAD,PHC_ROAD,PTC_ROAD,       &
                                        PD_ROAD,                            &
                                        PALB_WALL,                          &
                                        PEMIS_WALL,PHC_WALL,PTC_WALL,       &
                                        PD_WALL,                            &
                                        PBLD_HEIGHT,                        &
                                        PWALL_O_HOR,PBLD,                   &
                                        PH_TRAFFIC,  PLE_TRAFFIC,           &
                                        PH_INDUSTRY, PLE_INDUSTRY           )  
!     ##############################################################
!
!!**** *CONVERT_COVER* convert surface cover classes into secondary
!!                     physiographic variables for TEB
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
!     
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_DATA_COVER, ONLY : XDATA_Z0_TOWN, XDATA_ALB_ROOF,                    &
                              XDATA_EMIS_ROOF, XDATA_HC_ROOF, XDATA_TC_ROOF,  &
                              XDATA_D_ROOF, XDATA_ALB_ROAD, XDATA_EMIS_ROAD,  &
                              XDATA_HC_ROAD, XDATA_TC_ROAD, XDATA_D_ROAD,     &
                              XDATA_ALB_WALL, XDATA_EMIS_WALL, XDATA_HC_WALL, &
                              XDATA_TC_WALL, XDATA_D_WALL, XDATA_BLD_HEIGHT,  &
                              XDATA_H_TRAFFIC, XDATA_LE_TRAFFIC,              &
                              XDATA_H_INDUSTRY, XDATA_LE_INDUSTRY  
!
USE MODD_DATA_COVER_n, ONLY : XDATA_BLD, XDATA_WALL_O_HOR
!
USE MODI_AV_PGD
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PCOVER
!
REAL, DIMENSION(:),   INTENT(OUT), OPTIONAL   :: PZ0_TOWN
REAL, DIMENSION(:),   INTENT(OUT), OPTIONAL   :: PALB_ROOF
REAL, DIMENSION(:),   INTENT(OUT), OPTIONAL   :: PEMIS_ROOF
REAL, DIMENSION(:,:), INTENT(OUT), OPTIONAL   :: PHC_ROOF
REAL, DIMENSION(:,:), INTENT(OUT), OPTIONAL   :: PTC_ROOF
REAL, DIMENSION(:,:), INTENT(OUT), OPTIONAL   :: PD_ROOF
REAL, DIMENSION(:),   INTENT(OUT), OPTIONAL   :: PALB_ROAD
REAL, DIMENSION(:),   INTENT(OUT), OPTIONAL   :: PEMIS_ROAD
REAL, DIMENSION(:,:), INTENT(OUT), OPTIONAL   :: PHC_ROAD
REAL, DIMENSION(:,:), INTENT(OUT), OPTIONAL   :: PTC_ROAD
REAL, DIMENSION(:,:), INTENT(OUT), OPTIONAL   :: PD_ROAD
REAL, DIMENSION(:),   INTENT(OUT), OPTIONAL   :: PALB_WALL
REAL, DIMENSION(:),   INTENT(OUT), OPTIONAL   :: PEMIS_WALL
REAL, DIMENSION(:,:), INTENT(OUT), OPTIONAL   :: PHC_WALL
REAL, DIMENSION(:,:), INTENT(OUT), OPTIONAL   :: PTC_WALL
REAL, DIMENSION(:,:), INTENT(OUT), OPTIONAL   :: PD_WALL
REAL, DIMENSION(:),   INTENT(OUT), OPTIONAL   :: PBLD_HEIGHT
REAL, DIMENSION(:),   INTENT(OUT), OPTIONAL   :: PWALL_O_HOR
REAL, DIMENSION(:),   INTENT(OUT), OPTIONAL   :: PBLD
REAL, DIMENSION(:),   INTENT(OUT), OPTIONAL   :: PH_TRAFFIC
REAL, DIMENSION(:),   INTENT(OUT), OPTIONAL   :: PLE_TRAFFIC
REAL, DIMENSION(:),   INTENT(OUT), OPTIONAL   :: PH_INDUSTRY
REAL, DIMENSION(:),   INTENT(OUT), OPTIONAL   :: PLE_INDUSTRY
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER :: JLAYER ! loop counter on surface layers
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*    2.      SECONDARY VARIABLES
!             -------------------
!
!*    2.2     fields on artificial surfaces only
!             ----------------------------------
!
IF (LHOOK) CALL DR_HOOK('CONVERT_COVER_TEB',0,ZHOOK_HANDLE)
IF (PRESENT(PZ0_TOWN)) &
  CALL AV_PGD (PZ0_TOWN    ,PCOVER ,XDATA_Z0_TOWN    (:),'TWN','CDN')  
IF (PRESENT(PBLD)) &
  CALL AV_PGD (PBLD        ,PCOVER ,XDATA_BLD        (:),'TWN','ARI')  
IF (PRESENT(PALB_ROOF)) &
  CALL AV_PGD (PALB_ROOF   ,PCOVER ,XDATA_ALB_ROOF   (:),'BLD','ARI')  
IF (PRESENT(PEMIS_ROOF)) &
  CALL AV_PGD (PEMIS_ROOF  ,PCOVER ,XDATA_EMIS_ROOF  (:),'BLD','ARI')  
!
IF (PRESENT(PHC_ROOF)) THEN
  IF ( SIZE(PHC_ROOF,2) > SIZE(XDATA_HC_ROOF,2) ) CALL ABOR1_SFX('CONVERT_COVER_TEB: (1)')
  DO JLAYER=1,SIZE(PHC_ROOF,2)
    CALL AV_PGD (PHC_ROOF(:,JLAYER)    ,PCOVER ,XDATA_HC_ROOF    (:,JLAYER),'BLD','INV')
  END DO
ENDIF
IF (PRESENT(PTC_ROOF)) THEN
  IF ( SIZE(PTC_ROOF,2) > SIZE(XDATA_TC_ROOF,2) ) CALL ABOR1_SFX('CONVERT_COVER_TEB: (2)')
  DO JLAYER=1,SIZE(PTC_ROOF,2)
    CALL AV_PGD (PTC_ROOF(:,JLAYER)    ,PCOVER ,XDATA_TC_ROOF    (:,JLAYER),'BLD','ARI')
  END DO
ENDIF
IF (PRESENT(PD_ROOF)) THEN
  DO JLAYER=1,SIZE(PD_ROOF,2)
    IF ( SIZE(PD_ROOF,2) > SIZE(XDATA_D_ROOF,2) ) CALL ABOR1_SFX('CONVERT_COVER_TEB: (3)')
    CALL AV_PGD (PD_ROOF (:,JLAYER)    ,PCOVER ,XDATA_D_ROOF     (:,JLAYER),'BLD','ARI')
  END DO
ENDIF
!
IF (PRESENT(PALB_ROAD)) &
  CALL AV_PGD (PALB_ROAD   ,PCOVER ,XDATA_ALB_ROAD   (:),'STR','ARI')  
IF (PRESENT(PEMIS_ROAD)) &
  CALL AV_PGD (PEMIS_ROAD  ,PCOVER ,XDATA_EMIS_ROAD  (:),'STR','ARI')  
!
IF (PRESENT(PHC_ROAD)) THEN
  IF ( SIZE(PHC_ROAD,2) > SIZE(XDATA_HC_ROAD,2) ) CALL ABOR1_SFX('CONVERT_COVER_TEB: (4)')
  DO JLAYER=1,SIZE(PHC_ROAD,2)
    CALL AV_PGD (PHC_ROAD(:,JLAYER)    ,PCOVER ,XDATA_HC_ROAD    (:,JLAYER),'STR','INV')
  END DO
ENDIF
IF (PRESENT(PTC_ROAD)) THEN
  IF ( SIZE(PTC_ROAD,2) > SIZE(XDATA_TC_ROAD,2) ) CALL ABOR1_SFX('CONVERT_COVER_TEB: (5)')
  DO JLAYER=1,SIZE(PTC_ROAD,2)
    CALL AV_PGD (PTC_ROAD(:,JLAYER)    ,PCOVER ,XDATA_TC_ROAD    (:,JLAYER),'STR','ARI')
  END DO
ENDIF
IF (PRESENT(PD_ROAD)) THEN
  IF ( SIZE(PD_ROAD,2) > SIZE(XDATA_D_ROAD,2) ) CALL ABOR1_SFX('CONVERT_COVER_TEB: (6)')
  DO JLAYER=1,SIZE(PD_ROAD,2)
    CALL AV_PGD (PD_ROAD (:,JLAYER)    ,PCOVER ,XDATA_D_ROAD     (:,JLAYER),'STR','ARI')
  END DO
ENDIF
!
IF (PRESENT(PALB_WALL)) &
  CALL AV_PGD (PALB_WALL   ,PCOVER ,XDATA_ALB_WALL   (:),'BLD','ARI')  
IF (PRESENT(PEMIS_WALL)) &
  CALL AV_PGD (PEMIS_WALL  ,PCOVER ,XDATA_EMIS_WALL  (:),'BLD','ARI')  
!
IF (PRESENT(PHC_WALL)) THEN
  IF ( SIZE(PHC_WALL,2) > SIZE(XDATA_HC_WALL,2) ) CALL ABOR1_SFX('CONVERT_COVER_TEB: (7)')
  DO JLAYER=1,SIZE(PHC_WALL,2)
    CALL AV_PGD (PHC_WALL(:,JLAYER)    ,PCOVER ,XDATA_HC_WALL    (:,JLAYER),'BLD','INV')
  END DO
ENDIF
IF (PRESENT(PTC_WALL)) THEN
  IF ( SIZE(PTC_WALL,2) > SIZE(XDATA_TC_WALL,2) ) CALL ABOR1_SFX('CONVERT_COVER_TEB: (8)')
  DO JLAYER=1,SIZE(PTC_WALL,2)
    CALL AV_PGD (PTC_WALL(:,JLAYER)    ,PCOVER ,XDATA_TC_WALL    (:,JLAYER),'BLD','ARI')
  END DO
ENDIF
IF (PRESENT(PD_WALL)) THEN
  IF ( SIZE(PD_WALL,2) > SIZE(XDATA_D_WALL,2) ) CALL ABOR1_SFX('CONVERT_COVER_TEB: (9)')
  DO JLAYER=1,SIZE(PD_WALL,2)
    CALL AV_PGD (PD_WALL (:,JLAYER)    ,PCOVER ,XDATA_D_WALL     (:,JLAYER),'BLD','ARI')
  END DO
ENDIF
!
IF (PRESENT(PBLD_HEIGHT)) &
  CALL AV_PGD (PBLD_HEIGHT ,PCOVER ,XDATA_BLD_HEIGHT (:),'BLD','ARI')  
IF (PRESENT(PWALL_O_HOR)) &
  CALL AV_PGD (PWALL_O_HOR     ,PCOVER ,XDATA_WALL_O_HOR   (:),'BLD','ARI')  
IF (PRESENT(PH_TRAFFIC)) &
  CALL AV_PGD (PH_TRAFFIC      ,PCOVER ,XDATA_H_TRAFFIC    (:),'TWN','ARI')  
IF (PRESENT(PLE_TRAFFIC)) &
  CALL AV_PGD (PLE_TRAFFIC     ,PCOVER ,XDATA_LE_TRAFFIC   (:),'TWN','ARI')  
IF (PRESENT(PH_INDUSTRY)) &
  CALL AV_PGD (PH_INDUSTRY     ,PCOVER ,XDATA_H_INDUSTRY   (:),'TWN','ARI')  
IF (PRESENT(PLE_INDUSTRY)) &
  CALL AV_PGD (PLE_INDUSTRY    ,PCOVER ,XDATA_LE_INDUSTRY  (:),'TWN','ARI')  
IF (LHOOK) CALL DR_HOOK('CONVERT_COVER_TEB',1,ZHOOK_HANDLE)
!

!-------------------------------------------------------------------------------
!
END SUBROUTINE CONVERT_COVER_TEB
