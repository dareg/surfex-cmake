!     #########
    SUBROUTINE URBAN_LW_COEF(PEMIS_ROAD, PSVF_ROAD,                          &
                               PEMIS_WALL, PSVF_WALL,                          &
                               PEMIS_GARDEN, PSVF_GARDEN,                      &
                               PROAD, PGARDEN,                                 &
                               PDN_ROAD, PDF_ROAD, PESNOW_ROAD,                &
                               PLW_W_TO_W,  PLW_W_TO_R,  PLW_W_TO_G,           &
                               PLW_R_TO_W,  PLW_R_TO_R,  PLW_R_TO_G,           &
                               PLW_G_TO_W,  PLW_G_TO_R,  PLW_G_TO_G,           &
                               PLW_S_TO_W,  PLW_S_TO_R,  PLW_S_TO_G,           &
                               PLW_NR_TO_W, PLW_NR_TO_R, PLW_NR_TO_G           )  
!   ##########################################################################
!
!!****  *URBAN_LW_COEF*  
!!
!!    PURPOSE
!!    -------
!
!     Computes the coefficients before each of the temperatures in the
!     radiative budgets
!         
!     
!!**  METHOD
!     ------
!
! without snow, the radiative budgets read:
!
!   Rn_w = abs_Rg_w 
!  - sigma * emis_w                                                   * Ts_w**4 (t+dt)
!  +         emis_w                       *      SVF_w                * LWR
!  + sigma * emis_w * emis_r              *      SVF_w                * Ts_r**4 (t+dt)
!  + sigma * emis_w * emis_w              * (1-2*SVF_w)               * Ts_w**4 (t+dt)
!  +         emis_w            (1-emis_r) *      SVF_r  *      SVF_w  * LWR
!  +         emis_w            (1-emis_w) *      SVF_w  * (1-2*SVF_w) * LWR
!  + sigma * emis_w * emis_w * (1-emis_r) *      SVF_w  * (1-  SVF_r) * Ts_w**4 (t+dt)
!  + sigma * emis_w * emis_w * (1-emis_w) * (1-2*SVF_w) * (1-2*SVF_w) * Ts_w**4 (t+dt)
!  + sigma * emis_w * emis_r * (1-emis_w) *      SVF_w  * (1-2*SVF_w) * Ts_r**4 (t+dt)
!
!   Rn_r = abs_Rg_r
!  - sigma * emis_r                                                   * Ts_r**4 (t+dt)
!  +         emis_r                       *    SVF_r                  * LWR
!  + sigma * emis_r * emis_w              * (1-SVF_r)                 * Ts_w**4 (t+dt)
!  +         emis_r            (1-emis_w) * (1-SVF_r)   *      SVF_w  * LWR
!  + sigma * emis_r * emis_w * (1-emis_w) * (1-SVF_r)   * (1-2*SVF_w) * Ts_w**4 (t+dt)
!  + sigma * emis_r * emis_r * (1-emis_w) * (1-SVF_r)   *      SVF_w  * Ts_r**4 (t+dt)
!
!
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
!!      Original    08/09/98 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,ONLY : XSTEFAN
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
!
REAL, DIMENSION(:), INTENT(IN)  :: PEMIS_ROAD  ! road emissivity
REAL, DIMENSION(:), INTENT(IN)  :: PSVF_ROAD   ! road sky view factor
REAL, DIMENSION(:), INTENT(IN)  :: PEMIS_WALL  ! wall emissivity
REAL, DIMENSION(:), INTENT(IN)  :: PSVF_WALL   ! wall sky view factor
REAL, DIMENSION(:), INTENT(IN)  :: PEMIS_GARDEN! GARDEN area emissivity
REAL, DIMENSION(:), INTENT(IN)  :: PSVF_GARDEN ! GARDEN area sky view factor
REAL, DIMENSION(:), INTENT(IN)  :: PROAD       ! road fraction
REAL, DIMENSION(:), INTENT(IN)  :: PGARDEN     ! GARDEN area fraction
REAL, DIMENSION(:), INTENT(IN)  :: PDN_ROAD    ! road snow-covered surface fraction
REAL, DIMENSION(:), INTENT(IN)  :: PDF_ROAD    ! road snow-free surface fraction
REAL, DIMENSION(:), INTENT(IN)  :: PESNOW_ROAD ! road snow emissivity
!
REAL, DIMENSION(:), INTENT(OUT) :: PLW_W_TO_W  ! L.W. interactions wall->wall (from first Temp. on second Temp.)
REAL, DIMENSION(:), INTENT(OUT) :: PLW_W_TO_R  ! L.W. interactions wall->road
REAL, DIMENSION(:), INTENT(OUT) :: PLW_W_TO_G  ! L.W. interactions wall->GARDEN areas
REAL, DIMENSION(:), INTENT(OUT) :: PLW_R_TO_W  ! L.W. interactions road->wall
REAL, DIMENSION(:), INTENT(OUT) :: PLW_R_TO_R  ! L.W. interactions road->road
REAL, DIMENSION(:), INTENT(OUT) :: PLW_R_TO_G  ! L.W. interactions road->GARDEN areas
REAL, DIMENSION(:), INTENT(OUT) :: PLW_G_TO_W  ! L.W. interactions GARDEN areas->wall
REAL, DIMENSION(:), INTENT(OUT) :: PLW_G_TO_R  ! L.W. interactions GARDEN areas->road
REAL, DIMENSION(:), INTENT(OUT) :: PLW_G_TO_G  ! L.W. interactions GARDEN areas->GARDEN areas
REAL, DIMENSION(:), INTENT(OUT) :: PLW_S_TO_W  ! L.W. interactions sky->wall
REAL, DIMENSION(:), INTENT(OUT) :: PLW_S_TO_R  ! L.W. interactions sky->road
REAL, DIMENSION(:), INTENT(OUT) :: PLW_S_TO_G  ! L.W. interactions sky->GARDEN areas
REAL, DIMENSION(:), INTENT(OUT) :: PLW_NR_TO_W ! L.W. interactions snow(road)->wall
REAL, DIMENSION(:), INTENT(OUT) :: PLW_NR_TO_R ! L.W. interactions snow(road)->road
REAL, DIMENSION(:), INTENT(OUT) :: PLW_NR_TO_G ! L.W. interactions snow(road)->GARDEN areas
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PDN_ROAD))  :: ZAVEMIS_R   ! averaged road emissivity
!
REAL, DIMENSION(SIZE(PDN_ROAD))  :: ZROAD       !
REAL, DIMENSION(SIZE(PDN_ROAD))  :: ZGARDEN     !
!
INTEGER :: JJ
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('URBAN_LW_COEF',0,ZHOOK_HANDLE)
!
ZROAD(:)=0.
ZGARDEN(:)=0.
!
DO JJ=1,SIZE(PROAD)
!
  IF (PROAD(JJ)+PGARDEN(JJ).NE.0.) THEN
    ZROAD(JJ)  = PROAD(JJ)  / (PROAD(JJ)+PGARDEN(JJ))
    ZGARDEN(JJ) =  PGARDEN(JJ) / (PROAD(JJ)+PGARDEN(JJ))
  ELSE
    ZROAD(JJ)=0.
    ZGARDEN(JJ)=0.
  ENDIF
!
!
  ZAVEMIS_R(JJ) =  PDF_ROAD (JJ) * PEMIS_ROAD  (JJ) &
                + PDN_ROAD (JJ) * PESNOW_ROAD (JJ)  
!
!
  PLW_W_TO_W(JJ)  =     XSTEFAN * PEMIS_WALL(JJ)                                                               &
                     * ( - 1.                                                                                &
                         + PEMIS_WALL(JJ) * (1.-2.*PSVF_WALL(JJ))                                              &
                         + PEMIS_WALL(JJ) * (1.-PEMIS_WALL(JJ)) * (1.-2.*PSVF_WALL(JJ))**2                      &
                         + PEMIS_WALL(JJ) * (1.- ZAVEMIS_R(JJ)) * PSVF_WALL(JJ) * (1.-PSVF_ROAD(JJ)) *ZROAD(JJ)   &
                         + PEMIS_WALL(JJ) * (1.- PEMIS_GARDEN(JJ))*PSVF_WALL(JJ) * (1.-PSVF_GARDEN(JJ))*ZGARDEN(JJ)  )  
!
  PLW_W_TO_R(JJ)  =     XSTEFAN * PEMIS_WALL(JJ) * PEMIS_ROAD(JJ) * (1.-PSVF_ROAD(JJ))                           &
                     * (1. + (1.-PEMIS_WALL(JJ)) * (1.-2.*PSVF_WALL(JJ))                                       )  
!
  PLW_W_TO_G(JJ)  =     XSTEFAN * PEMIS_WALL(JJ) * PEMIS_GARDEN(JJ) * (1.-PSVF_GARDEN(JJ))                       &
                     * (1. + (1.-PEMIS_WALL(JJ)) * (1.-2.*PSVF_WALL(JJ))                                       )  
!
  PLW_R_TO_W(JJ)  =     XSTEFAN * PEMIS_WALL(JJ) * PEMIS_ROAD(JJ) * PSVF_WALL(JJ) * ZROAD(JJ) * PDF_ROAD(JJ)       &
                     * ( 1. + (1.-PEMIS_WALL(JJ)) * (1.-2.*PSVF_WALL(JJ))                                      )  
!
  PLW_R_TO_R(JJ)  =     XSTEFAN * PEMIS_ROAD(JJ)                                                               &
                     * (-1. + PEMIS_ROAD(JJ) * (1.-PEMIS_WALL(JJ))                                             &
                              * (1.-PSVF_ROAD(JJ)) * PSVF_WALL(JJ) * ZROAD(JJ) * PDF_ROAD(JJ)                    )  
!
  PLW_R_TO_G(JJ)  =     XSTEFAN * PEMIS_ROAD(JJ) * PEMIS_GARDEN(JJ) * (1.-PEMIS_WALL(JJ))                        &
                     * (1.-PSVF_GARDEN(JJ)) * PSVF_WALL(JJ) * ZROAD(JJ) * PDF_ROAD(JJ)  
!
  PLW_G_TO_W(JJ)  =     XSTEFAN * PEMIS_WALL(JJ) * PEMIS_GARDEN(JJ) * PSVF_WALL(JJ) * ZGARDEN(JJ)                 &
                     * (1. + (1.-PEMIS_WALL(JJ)) * (1.-2.*PSVF_WALL(JJ))                                       )  
!
  PLW_G_TO_R(JJ)  =     XSTEFAN * PEMIS_ROAD(JJ) * PEMIS_GARDEN(JJ) * (1.-PEMIS_WALL(JJ))                        &
                     * (1.-PSVF_ROAD(JJ)) * PSVF_WALL(JJ) * ZGARDEN(JJ)  
!
  PLW_G_TO_G(JJ)  =     XSTEFAN * PEMIS_GARDEN(JJ)                                                             &
                     * (-1. + PEMIS_GARDEN(JJ) * (1.-PEMIS_WALL(JJ))                                           &
                              * (1.-PSVF_GARDEN(JJ)) * PSVF_WALL(JJ) * ZGARDEN(JJ)                              )  
!
  PLW_S_TO_W(JJ)  =     PEMIS_WALL(JJ) * PSVF_WALL(JJ)                                                          &
                     * (1. + (1.-ZAVEMIS_R(JJ))  * PSVF_ROAD(JJ) * ZROAD(JJ)                                       &
                           + (1.-PEMIS_GARDEN(JJ))* PSVF_GARDEN(JJ) * ZGARDEN(JJ)                                  &
                           + (1.-PEMIS_WALL(JJ)) * (1.-2.*PSVF_WALL(JJ))                                       )  
!
  PLW_S_TO_R(JJ)  =     PEMIS_ROAD(JJ) * PSVF_ROAD(JJ)                                                          &
                     + PEMIS_ROAD(JJ) * (1. - PEMIS_WALL(JJ)) * PSVF_WALL(JJ) * (1.-PSVF_ROAD(JJ))  
!
  PLW_S_TO_G(JJ)  =     PEMIS_GARDEN(JJ) * PSVF_GARDEN(JJ)                                                      &
                     + PEMIS_GARDEN(JJ) * (1. - PEMIS_WALL(JJ)) * PSVF_WALL(JJ) * (1.-PSVF_GARDEN(JJ))  
!
  PLW_NR_TO_W(JJ) =     XSTEFAN * PEMIS_WALL(JJ) * PESNOW_ROAD(JJ) * PSVF_WALL(JJ) * ZROAD(JJ) * PDN_ROAD(JJ)      &
                     * ( 1. + (1.-PEMIS_WALL(JJ)) * (1.-2.*PSVF_WALL(JJ))                                      )  
!
  PLW_NR_TO_R(JJ) =     XSTEFAN * PEMIS_ROAD(JJ) * (1.-PEMIS_WALL(JJ)) * PESNOW_ROAD(JJ)                         &
                     * (1.-PSVF_ROAD(JJ)) * PSVF_WALL(JJ) * ZROAD(JJ) * PDN_ROAD(JJ)  
!
  PLW_NR_TO_G(JJ) =     XSTEFAN * PEMIS_GARDEN(JJ) * (1.-PEMIS_WALL(JJ)) * PESNOW_ROAD(JJ)                        &
                     * (1.-PSVF_GARDEN(JJ)) * PSVF_WALL(JJ) * ZROAD(JJ) * PDN_ROAD(JJ)  
!
ENDDO
IF (LHOOK) CALL DR_HOOK('URBAN_LW_COEF',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE URBAN_LW_COEF
