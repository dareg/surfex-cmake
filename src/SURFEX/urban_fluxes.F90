!     #########
    SUBROUTINE URBAN_FLUXES(TOP, T, B, DGMT, HIMPLICIT_WIND, PT_CANYON, PPEW_A_COEF, PPEW_B_COEF, &
                            PEXNS, PRHOA, PVMOD, PH_TRAFFIC, PLE_TRAFFIC, PAC_WL, PCD, PDF_RF, &
                            PDN_RF, PDF_RD, PDN_RD, PRNSN_RF, PHSN_RF, PLESN_RF, PGSN_RF, &
                            PRNSN_RD, PHSN_RD, PLESN_RD, PGSN_RD, PMELT_RF, PDQS_RF, PMELT_RD, &
                            PDQS_RD, PDQS_WL_A, PDQS_WL_B, PFLX_BLD_RF, PFLX_BLD_WL_A, &
                            PFLX_BLD_WL_B, PFLX_BLD_FL, PFLX_BLD_MA, PE_SHADING, PLEW_RF, PLEW_RD, &
                            PLE_WL_A, PLE_WL_B, PMELT_BLT, PUSTAR_TWN     )
!   ##########################################################################
!
!!****  *URBAN_FLUXES* computes fluxes on urbanized surfaces  
!!
!!    PURPOSE
!!    -------
!         
!     
!!**  METHOD
!     ------
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
!!      V. Masson           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    23/01/98 
!!                     12/02 (A. Lemonsu) modifications of emissivity and Tstown
!!                     07/07 (P. LeMoigne) expression of latent heat fluxes as 
!!                           functions of w'theta' instead of w'T'
!!                     17/10 (G. Pigeon)  computation of anthropogenic heat due
!!                            to domestic heating
!!      Modified    09/2012 : B. Decharme New wind implicitation
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_n, ONLY : TEB_1P_t
USE MODD_BEM_n, ONLY : BEM_1P_t
USE MODD_DIAG_MISC_TEB_n, ONLY : DIAG_MISC_TEB_1P_t
!
USE MODD_SURF_PAR, ONLY : XUNDEF
USE MODD_CSTS,     ONLY : XCPD
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_1P_t), INTENT(INOUT) :: T
TYPE(BEM_1P_t), INTENT(INOUT) :: B
TYPE(DIAG_MISC_TEB_1P_t), INTENT(INOUT) :: DGMT
!
 CHARACTER(LEN=*),     INTENT(IN)  :: HIMPLICIT_WIND   ! wind implicitation option
!                                                     ! 'OLD' = direct
!                                                     ! 'NEW' = Taylor serie, order 1
REAL, DIMENSION(:), INTENT(IN)    :: PT_CANYON    ! canyon air temperature

REAL, DIMENSION(:), INTENT(IN)    :: PPEW_A_COEF  ! implicit coefficients
REAL, DIMENSION(:), INTENT(IN)    :: PPEW_B_COEF  ! for wind coupling
REAL, DIMENSION(:), INTENT(IN)    :: PEXNS        ! surface exner function

REAL, DIMENSION(:), INTENT(IN)    :: PRHOA        ! air density
                                                  ! at the lowest level
REAL, DIMENSION(:), INTENT(IN)    :: PVMOD        ! module of the horizontal wind
REAL, DIMENSION(:), INTENT(IN)    :: PH_TRAFFIC    ! anthropogenic sensible
!                                                  ! heat fluxes due to traffic
REAL, DIMENSION(:), INTENT(IN)    :: PLE_TRAFFIC   ! anthropogenic latent
!                                                  ! heat fluxes due to traffic
REAL, DIMENSION(:), INTENT(IN)    :: PAC_WL     ! surface conductance
!                                                 ! for heat transfer
!                                                 ! between wall and canyon
REAL, DIMENSION(:), INTENT(IN)    :: PCD          ! drag coefficient
REAL, DIMENSION(:), INTENT(IN)    :: PDF_RF     ! snow-free    roof fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDN_RF     ! snow-covered roof fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDF_RD     ! snow-free    road fraction
REAL, DIMENSION(:), INTENT(IN)    :: PDN_RD     ! snow-covered road fraction
!
REAL, DIMENSION(:), INTENT(IN)    :: PRNSN_RF ! net radiation over snow
REAL, DIMENSION(:), INTENT(IN)    :: PHSN_RF  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PLESN_RF ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PGSN_RF  ! flux under the snow
REAL, DIMENSION(:), INTENT(IN)    :: PRNSN_RD ! net radiation over snow
REAL, DIMENSION(:), INTENT(IN)    :: PHSN_RD  ! sensible heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PLESN_RD ! latent heat flux over snow
REAL, DIMENSION(:), INTENT(IN)    :: PGSN_RD  ! flux under the snow
!
REAL, DIMENSION(:), INTENT(IN)    :: PMELT_RF   ! snow melting on roof
REAL, DIMENSION(:), INTENT(IN)    :: PDQS_RF    ! storage inside roofs
REAL, DIMENSION(:), INTENT(IN)    :: PMELT_RD   ! snow melting on road
REAL, DIMENSION(:), INTENT(IN)    :: PDQS_RD    ! storage inside roads
REAL, DIMENSION(:), INTENT(IN)    :: PDQS_WL_A  ! storage inside walls
REAL, DIMENSION(:), INTENT(IN)    :: PDQS_WL_B  ! storage inside walls
REAL, DIMENSION(:), INTENT(IN)    :: PFLX_BLD_RF! heat flx from bld to roof
REAL, DIMENSION(:), INTENT(IN)    :: PFLX_BLD_WL_A! heat flx from bld to wall A
REAL, DIMENSION(:), INTENT(IN)    :: PFLX_BLD_WL_B! heat flx from bld to wall B
REAL, DIMENSION(:), INTENT(IN)    :: PFLX_BLD_FL! heat flx from bld to floor
REAL, DIMENSION(:), INTENT(IN)    :: PFLX_BLD_MA! heat flx from bld to mass
REAL, DIMENSION(:), INTENT(IN)    :: PE_SHADING   ! energy not ref., nor absorbed, nor
                                                  !trans. by glazing [W
                                                  !m-2(win)]
!
REAL, DIMENSION(:), INTENT(IN)   :: PLEW_RF     ! latent heat flux over snow-free roof
REAL, DIMENSION(:), INTENT(IN)   :: PLEW_RD     ! latent heat flux of snow-free road
!
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_WL_A   ! latent heat flux over wall
REAL, DIMENSION(:), INTENT(OUT)   :: PLE_WL_B   ! latent heat flux over wall
!
REAL, DIMENSION(:), INTENT(OUT)   :: PMELT_BLT    ! snow melting for town
!
REAL, DIMENSION(:), INTENT(OUT)   :: PUSTAR_TWN  ! friction velocity over town
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PRHOA)) :: ZH_RF_SNFREE
REAL, DIMENSION(SIZE(PRHOA)) :: ZRN_RF_SNFREE
REAL, DIMENSION(SIZE(PRHOA)) :: ZUSTAR2 ! square of friction velocity (m2/s2)
REAL, DIMENSION(SIZE(PRHOA)) :: ZVMOD   ! Wind
REAL, DIMENSIOn(SIZE(PRHOA)) :: ZINTER
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('URBAN_FLUXES',0,ZHOOK_HANDLE)
!
!*      1.     Fluxes at snow-free roofs
!              -------------------------
!
!                                            net radiation
!
!
ZRN_RF_SNFREE(:) = DGMT%XABS_SW_ROOF(:) + DGMT%XABS_LW_ROOF(:)
!
!                                            sensible heat flux
!
ZH_RF_SNFREE(:) = DGMT%XH_ROOF(:)
!
!-------------------------------------------------------------------------------
!
!*      2.     Fluxes at snow-free roads
!              -------------------------
!
!                                            net radiation
!
DGMT%XRN_ROAD(:) = DGMT%XABS_SW_ROAD(:) + DGMT%XABS_LW_ROAD(:)
!
!-------------------------------------------------------------------------------
!
!*      3.     Fluxes at walls
!              ---------------
!
!                                            net radiation
!
DGMT%XRN_WALL_A(:) = DGMT%XABS_SW_WALL_A(:) + DGMT%XABS_LW_WALL_A(:)
DGMT%XRN_WALL_B(:) = DGMT%XABS_SW_WALL_B(:) + DGMT%XABS_LW_WALL_B(:)
!
IF (TOP%CBEM=="BEM") THEN
  ZINTER(:) = DGMT%XABS_SW_WIN (:) + DGMT%XABS_LW_WIN (:)
  DGMT%XRN_WALL_A(:) = DGMT%XRN_WALL_A(:) * (1.-B%XGR(:))  + ZINTER(:) * B%XGR(:)
  DGMT%XRN_WALL_B(:) = DGMT%XRN_WALL_B(:) * (1.-B%XGR(:))  + ZINTER(:) * B%XGR(:)
ENDIF
!
!                                            heat flux into the ground
!
DGMT%XGFLUX_WALL_A(:) = DGMT%XRN_WALL_A(:) - DGMT%XH_WALL_A(:)            
DGMT%XGFLUX_WALL_B(:) = DGMT%XRN_WALL_B(:) - DGMT%XH_WALL_B(:)            
!
!                                            sensible heat flux
!
!before -> PH_WL in [W.m-2(wall)]
ZINTER(:) = PAC_WL(:)*XCPD*PRHOA(:)/PEXNS(:) * (B%XT_WIN1(:)-PT_CANYON(:)) + PE_SHADING(:)
!
DGMT%XH_WALL_A(:) = (1. - B%XGR(:)) * DGMT%XH_WALL_A(:) + B%XGR(:) * ZINTER(:)
DGMT%XH_WALL_B(:) = (1. - B%XGR(:)) * DGMT%XH_WALL_B(:) + B%XGR(:) * ZINTER(:)
!
IF (TOP%CBEM=="BEM") THEN
  ZINTER(:) = B%XF_WASTE_CAN(:) * DGMT%XH_WASTE(:) / T%XWALL_O_HOR(:)
  DGMT%XH_WALL_A(:) = DGMT%XH_WALL_A(:) + ZINTER(:)
  DGMT%XH_WALL_B(:) = DGMT%XH_WALL_B(:) + ZINTER(:)
ENDIF
!
IF (TOP%CBEM=="BEM") THEN
  !after PH_WALL in [W.m-2(facade=wall + win)]
  ZINTER(:) = B%XF_WASTE_CAN(:) * DGMT%XLE_WASTE(:) / T%XWALL_O_HOR(:)
ELSE
  ZINTER(:) = 0.
ENDIF
PLE_WL_A(:) = ZINTER(:)
PLE_WL_B(:) = ZINTER(:)
!
!-------------------------------------------------------------------------------
!
!*      4.     Snow-free and snow-covered surfaces averaging
!              ---------------------------------------------
!
!*      4.1    Roads
!              -----
!
!                                            heat flux into the ground
!
DGMT%XGFLUX_ROAD (:) =  PDF_RD(:) * (DGMT%XRN_ROAD(:) - DGMT%XH_ROAD (:) - PLEW_RD(:) )&
                       + PDN_RD(:) * PGSN_RD(:)  
!
!                                            net radiation
!
DGMT%XRN_ROAD(:) = DGMT%XRN_ROAD(:) * PDF_RD(:) + PRNSN_RD(:) * PDN_RD(:)
!
!                                            sensible heat flux
!                                            total latent heat of evaporation from
!                                            the road (snow free + snow)
!
! sensible heat flux
!
DGMT%XH_ROAD  (:) = DGMT%XH_ROAD (:) * PDF_RD(:) + PHSN_RD(:) * PDN_RD(:)
! total latent heat of evaporation from  the road (snow free + snow)
!
DGMT%XLE_ROAD (:) = PLEW_RD(:) * PDF_RD(:) + PLESN_RD(:) * PDN_RD(:)
!
!*      4.2    Roofs
!              -----
!
!                                            heat flux into the ground
!
DGMT%XGFLUX_STRLROOF(:) =  PDF_RF(:) * (ZRN_RF_SNFREE(:) - ZH_RF_SNFREE(:) - PLEW_RF(:)) + PDN_RF(:) *  PGSN_RF(:)  
DGMT%XGFLUX_ROOF    (:) = (1.-T%XGREENROOF(:)) * DGMT%XGFLUX_STRLROOF(:) + T%XGREENROOF(:) * DGMT%XGFLUX_GREENROOF(:)
!
!
!                                            net radiation
! 
DGMT%XRN_STRLROOF   (:) = ZRN_RF_SNFREE(:) * PDF_RF(:) + PRNSN_RF(:) * PDN_RF(:)
DGMT%XRN_ROOF       (:) = (1.-T%XGREENROOF(:)) * DGMT%XRN_STRLROOF(:) + T%XGREENROOF(:) * DGMT%XRN_GREENROOF(:)
!
!                                            sensible heat flux
!                                            total latent heat of evaporation from
!                                            the roof (snow free + snow)
!
! sensible heat flux
!
DGMT%XH_STRLROOF    (:) = DGMT%XH_ROOF(:) * PDF_RF(:) + PHSN_RF(:)  * PDN_RF(:)
DGMT%XH_ROOF        (:) = (1.-T%XGREENROOF(:)) * DGMT%XH_STRLROOF(:) + T%XGREENROOF(:) * DGMT%XH_GREENROOF(:)
!
!
! total latent heat of evaporation from the roof (snow free + snow)
!
DGMT%XLE_STRLROOF   (:) = PLEW_RF(:) * PDF_RF(:)  + PLESN_RF(:) * PDN_RF(:)
DGMT%XLE_ROOF       (:) = (1.-T%XGREENROOF(:)) * DGMT%XLE_STRLROOF(:) + T%XGREENROOF(:) * DGMT%XLE_GREENROOF(:) 
!
IF (TOP%CBEM=="BEM") THEN
  DGMT%XH_ROOF (:) = DGMT%XH_ROOF (:) + (1 - B%XF_WASTE_CAN(:)) * DGMT%XH_WASTE (:)/T%XBLD(:)
  DGMT%XLE_ROOF(:) = DGMT%XLE_ROOF(:) + (1 - B%XF_WASTE_CAN(:)) * DGMT%XLE_WASTE(:)/T%XBLD(:)
ENDIF

!-------------------------------------------------------------------------------
!
!*      5.     Momentum fluxes
!              ---------------
!
PUSTAR_TWN(:) = XUNDEF
!
IF (.NOT. TOP%LCANOPY) THEN
  !
  ZUSTAR2(:) = 0.0
  ZVMOD  (:) = PVMOD(:)
  !
  IF(HIMPLICIT_WIND=='OLD')THEN 
    ! old implicitation
    ZUSTAR2(:) = (PCD(:)*PVMOD(:)*PPEW_B_COEF(:))/ (1.0-PRHOA(:)*PCD(:)*PVMOD(:)*PPEW_A_COEF(:))
  ELSE
    ! new implicitation
    ZUSTAR2(:) = (PCD(:)*PVMOD(:)*(2.*PPEW_B_COEF(:)-PVMOD(:)))/ &
                 (1.0-2.0*PRHOA(:)*PCD(:)*PVMOD(:)*PPEW_A_COEF(:))
    !                   
    ZVMOD(:) = PRHOA(:)*PPEW_A_COEF(:)*ZUSTAR2(:) + PPEW_B_COEF(:)
    ZVMOD(:) = MAX(ZVMOD(:),0.)
    !
    WHERE(PPEW_A_COEF(:)/= 0.)
      ZUSTAR2(:) = MAX( ( ZVMOD(:) - PPEW_B_COEF(:) ) / (PRHOA(:)*PPEW_A_COEF(:)), 0.)
    ENDWHERE
    !               
  ENDIF
  !
  PUSTAR_TWN(:) = SQRT(ZUSTAR2(:))     
  !  
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      6.     Averaged fluxes
!              ---------------
!
!*      6.1    Built fraction
!              --------------
!
DGMT%XRN_BLT (:)    = ( T%XBLD(:)        * DGMT%XRN_ROOF(:)        &
                  +     T%XROAD(:)       * DGMT%XRN_ROAD(:)        &
                  + 0.5*T%XWALL_O_HOR(:) * DGMT%XRN_WALL_A(:)      &
                  + 0.5*T%XWALL_O_HOR(:) * DGMT%XRN_WALL_B(:))     &
                  / (T%XROAD(:) + T%XBLD(:))
!
DGMT%XH_BLT  (:)    = ( T%XBLD(:)        * DGMT%XH_ROOF(:)         &
                  +     T%XROAD(:)       * DGMT%XH_ROAD(:)         &
                  + 0.5*T%XWALL_O_HOR(:) * DGMT%XH_WALL_A(:)       &   
                  + 0.5*T%XWALL_O_HOR(:) * DGMT%XH_WALL_B(:))      &   
                  / (T%XROAD(:) + T%XBLD(:))
!
DGMT%XLE_BLT (:)    = ( T%XBLD(:)        * DGMT%XLE_ROOF (:)       &
                  +     T%XROAD(:)       * DGMT%XLE_ROAD (:)       &
                  + 0.5*T%XWALL_O_HOR(:) * PLE_WL_A (:)     & 
                  + 0.5*T%XWALL_O_HOR(:) * PLE_WL_B (:))    & 
                  / (T%XROAD(:) + T%XBLD(:))
!
DGMT%XGFLUX_BLT (:) = ( T%XBLD(:)        * DGMT%XGFLUX_ROOF (:)    &
                  +     T%XROAD(:)       * DGMT%XGFLUX_ROAD (:)    &
                  + 0.5*T%XWALL_O_HOR(:) * DGMT%XGFLUX_WALL_A (:)  &
                  + 0.5*T%XWALL_O_HOR(:) * DGMT%XGFLUX_WALL_B (:)) &
                  / (T%XROAD(:) + T%XBLD(:))
!
PMELT_BLT  (:) = (     T%XBLD(:)         * PMELT_RF(:) * PDN_RF(:)  &
                  +    T%XROAD(:)        * PMELT_RD(:) * PDN_RD(:)) &
                  / (T%XROAD(:) + T%XBLD(:))
!
DGMT%XDQS_TOWN  (:) = (  T%XBLD(:)         * PDQS_RF (:)      &
                  +      T%XROAD(:)      * PDQS_RD (:)      &
                  + 0.5*T%XWALL_O_HOR(:) * PDQS_WL_A (:)    &
                  + 0.5*T%XWALL_O_HOR(:) * PDQS_WL_B (:) )  &
                / (T%XROAD(:) + T%XBLD(:))
!
IF (TOP%CBEM == "DEF") THEN
  DGMT%XQF_BLD(:) = ( ZH_RF_SNFREE(:) + PLEW_RF(:) +          & 
                 PDQS_RF(:) - ZRN_RF_SNFREE(:)  ) * PDF_RF(:) &
             + ( PDQS_RF(:) - PGSN_RF      (:)  ) * PDN_RF(:) &
             + 0.5*T%XWALL_O_HOR(:)/T%XBLD(:) * ( DGMT%XH_WALL_A(:) + PLE_WL_A(:) + PDQS_WL_A(:) - DGMT%XRN_WALL_A(:) ) &
             + 0.5*T%XWALL_O_HOR(:)/T%XBLD(:) * ( DGMT%XH_WALL_B(:) + PLE_WL_B(:) + PDQS_WL_B(:) - DGMT%XRN_WALL_B(:) )
  DGMT%XFLX_BLD(:)= XUNDEF
ELSEIF (TOP%CBEM == "BEM") THEN
  DGMT%XQF_BLD(:) = DGMT%XQIN(:)*B%XN_FLOOR(:) + DGMT%XHVAC_COOL(:) + DGMT%XHVAC_HEAT(:)
  DGMT%XFLX_BLD(:)=  PFLX_BLD_RF(:) + 0.5*T%XWALL_O_HOR(:)/T%XBLD(:)*PFLX_BLD_WL_A(:) &
                 +                    0.5*T%XWALL_O_HOR(:)/T%XBLD(:)*PFLX_BLD_WL_B(:) &
                 + PFLX_BLD_FL(:) + PFLX_BLD_MA(:)  
ENDIF
!
DGMT%XQF_TOWN(:)= T%XBLD(:)*DGMT%XQF_BLD(:) + PH_TRAFFIC(:) + T%XH_INDUSTRY(:) + PLE_TRAFFIC(:) + T%XLE_INDUSTRY(:)
!
!Flux from the building to its structure -> we need to add the component to the
!floor, the mass and the window
!
!
IF (LHOOK) CALL DR_HOOK('URBAN_FLUXES',1,ZHOOK_HANDLE)
!
END SUBROUTINE URBAN_FLUXES
