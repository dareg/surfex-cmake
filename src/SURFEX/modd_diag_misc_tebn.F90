!     ############################
      MODULE MODD_DIAG_MISC_TEB_n
!     ############################
!
!!****  *MODD_DIAG_MISC_TEB - declaration of packed surface parameters for TEB scheme
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      P. Le Moigne   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       07/10/04
!!      C de Munck        02/13  adding runoff contributions for teb garden  
!!      V. Masson      06/2013  splits module in two
!
!
!*       0.   DECLARATIONS
!             ------------
!
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
TYPE DIAG_MISC_TEB_1P_t
!
!* miscellaneous variables
!
  REAL, POINTER, DIMENSION(:)   :: XZ0_TOWN  ! town roughness length
  REAL, POINTER, DIMENSION(:)   :: XQF_BLD   ! domestic heating
  REAL, POINTER, DIMENSION(:)   :: XFLX_BLD ! heat flux from bld
  REAL, POINTER, DIMENSION(:)   :: XQF_TOWN  ! total anthropogenic heat
  REAL, POINTER, DIMENSION(:)   :: XDQS_TOWN ! storage inside building
!
  REAL, POINTER, DIMENSION(:)   :: XH_WALL_A   ! wall sensible heat flux          (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XH_WALL_B   ! wall sensible heat flux          (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XH_ROOF     ! roof sensible heat flux          (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XH_ROAD     ! road sensible heat flux          (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XH_STRLROOF ! structural roof sens. heat flux  (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XH_BLT      ! built surf sensible heat flux    (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XRN_WALL_A  ! net radiation at wall            (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XRN_WALL_B  ! net radiation at wall            (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XRN_ROOF    ! net radiation at roof            (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XRN_ROAD    ! net radiation at road            (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XRN_STRLROOF !net radiation at structural roofs(W/m2)
  REAL, POINTER, DIMENSION(:)   :: XRN_BLT     ! net radiation at built surf      (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XGFLUX_WALL_A !net wall conduction flux        (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XGFLUX_WALL_B !net wall conduction flux        (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XGFLUX_ROOF ! net roof conduction flux         (W/m2)                                         
  REAL, POINTER, DIMENSION(:)   :: XGFLUX_ROAD ! net road conduction flux         (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XGFLUX_STRLROOF !net structural roof cond flux (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XGFLUX_BLT  ! net built surf conduction flux   (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XLE_ROOF    ! roof latent heat flux            (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XLE_ROAD    ! road latent heat flux            (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XLE_STRLROOF !structural roof latent heat flux (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XLE_BLT     ! built surf latent heat flux      (W/m2)
!
  REAL, POINTER, DIMENSION(:)   :: XRUNOFF_TOWN      ! aggregated water runoff for town      (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XIRRIG_GARDEN     ! summer ground irrigation rate         (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XRUNOFF_ROAD      ! water runoff for roads                (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XIRRIG_ROAD       ! road man-made watering rate           (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XRUNOFF_ROOF      ! aggregated water runoff for roofs     (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XRUNOFF_STRLROOF  ! water runoff for structural roofs     (kg/m2/s)
  REAL, POINTER, DIMENSION(:)   :: XIRRIG_GREENROOF  ! summer ground irrigation rate         (kg/m2/s)
!
  REAL, POINTER, DIMENSION(:)   :: XABS_SW_ROOF      ! absorbed shortwave radiation over roofs
  REAL, POINTER, DIMENSION(:)   :: XABS_SW_SNOW_ROOF ! absorbed longwave radiation over roofs
  REAL, POINTER, DIMENSION(:)   :: XABS_LW_ROOF      ! absorbed shortwave radiation over roofs
  REAL, POINTER, DIMENSION(:)   :: XABS_LW_SNOW_ROOF ! absorbed longwave radiation over roofs
  REAL, POINTER, DIMENSION(:)   :: XABS_SW_ROAD      ! absorbed shortwave radiation over roads
  REAL, POINTER, DIMENSION(:)   :: XABS_SW_SNOW_ROAD ! absorbed longwave radiation over roads
  REAL, POINTER, DIMENSION(:)   :: XABS_LW_ROAD      ! absorbed shortwave radiation over roads
  REAL, POINTER, DIMENSION(:)   :: XABS_LW_SNOW_ROAD ! absorbed longwave radiation over roads
  REAL, POINTER, DIMENSION(:)   :: XABS_SW_WALL_A    ! absorbed shortwave radiation over walls
  REAL, POINTER, DIMENSION(:)   :: XABS_SW_WALL_B    ! absorbed shortwave radiation over walls
  REAL, POINTER, DIMENSION(:)   :: XABS_LW_WALL_A    ! absorbed shortwave radiation over walls
  REAL, POINTER, DIMENSION(:)   :: XABS_LW_WALL_B    ! absorbed shortwave radiation over walls
  REAL, POINTER, DIMENSION(:)   :: XABS_SW_GARDEN    ! absorbed shortwave radiation over green areas
  REAL, POINTER, DIMENSION(:)   :: XABS_LW_GARDEN    ! absorbed shortwave radiation over green areas
  REAL, POINTER, DIMENSION(:)   :: XABS_SW_GREENROOF ! absorbed shortwave radiation over green roofs
  REAL, POINTER, DIMENSION(:)   :: XABS_LW_GREENROOF ! absorbed shortwave radiation over green roofs
  REAL, POINTER, DIMENSION(:)   :: XG_GREENROOF_ROOF ! Heat flux between green roof and structural roof
  REAL, POINTER, DIMENSION(:)   :: XABS_SW_PANEL     ! absorbed shortwave radiation over solar panels
  REAL, POINTER, DIMENSION(:)   :: XABS_LW_PANEL     ! absorbed longwave  radiation over solar panels
!
  REAL, POINTER, DIMENSION(:)   :: XRN_PANEL         ! net radiation           over solar panels (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XH_PANEL          ! sensible heat flux      over solar panels (W/m2)
  REAL, POINTER, DIMENSION(:)   :: XTHER_PROD_PANEL  ! thermal      production of   solar panels (W/m2 thermal panel)
  REAL, POINTER, DIMENSION(:)   :: XPHOT_PROD_PANEL  ! photovoltaic production of   solar panels (W/m2 photovoltaic panel)
  REAL, POINTER, DIMENSION(:)   :: XPROD_PANEL       !              production of   solar panels (W/m2 panel)
  REAL, POINTER, DIMENSION(:)   :: XTHER_PROD_BLD    ! thermal      production of   solar panels (W/m2 bld)
  REAL, POINTER, DIMENSION(:)   :: XPHOT_PROD_BLD    ! photovoltaic production of   solar panels (W/m2 bld)

  REAL, POINTER, DIMENSION(:)   :: XH_BLD_COOL       ! Sensible cooling energy demand  
                                                     ! of the building [W m-2(bld)]
  REAL, POINTER, DIMENSION(:)   :: XT_BLD_COOL       ! Total cooling energy demand  
                                                     ! of the building [W m-2(bld)]
  REAL, POINTER, DIMENSION(:)   :: XH_BLD_HEAT       ! Heating energy demand       
                                                     ! of the building [W m-2(bld)]
  REAL, POINTER, DIMENSION(:)   :: XLE_BLD_COOL      ! Latent cooling energy demand 
                                                     ! of the building [W m-2(bld)]
  REAL, POINTER, DIMENSION(:)   :: XLE_BLD_HEAT      ! Latent heating energy demand 
                                                     ! of the building [W m-2(bld)]
  REAL, POINTER, DIMENSION(:)   :: XH_WASTE          ! Sensible waste heat from HVAC system
                                                     ! [W m-2(tot)]
  REAL, POINTER, DIMENSION(:)   :: XLE_WASTE         ! Latent waste heat from HVAC system
                                                     ! [W m-2(tot)]
  REAL, POINTER, DIMENSION(:)   :: XHVAC_COOL        ! Energy consumption of the cooling system
                                                     ! [W m-2(bld)]
  REAL, POINTER, DIMENSION(:)   :: XHVAC_HEAT        ! Energy consumption of the heating system
                                                     ! [W m-2(bld)]
  REAL, POINTER, DIMENSION(:)   :: XCAP_SYS          ! Actual capacity of the cooling system
                                                     ! [W m-2(bld)] 
  REAL, POINTER, DIMENSION(:)   :: XM_SYS            ! Actual HVAC mass flow rate 
                                                     ! [kg s-1 m-2(bld)]
  REAL, POINTER, DIMENSION(:)   :: XCOP              ! COP of the cooling system
  REAL, POINTER, DIMENSION(:)   :: XQ_SYS            ! Supply air specific humidity [kg kg-1]
  REAL, POINTER, DIMENSION(:)   :: XT_SYS            ! Supply air temperature [K]
  REAL, POINTER, DIMENSION(:)   :: XTR_SW_WIN        ! Solar radiation transmitted throught
                                                     ! windows [W m-2(bld)]
  REAL, POINTER, DIMENSION(:)   :: XFAN_POWER        ! HVAC fan power
  REAL, POINTER, DIMENSION(:)   :: XABS_SW_WIN       ! window absorbed shortwave radiation [W m-2] 
  REAL, POINTER, DIMENSION(:)   :: XABS_LW_WIN       ! absorbed infrared rad. [W m-2]

  REAL, POINTER, DIMENSION(:)   :: XEMIT_LW_FAC      ! LW flux emitted by the facade (W/m2 facade)
  REAL, POINTER, DIMENSION(:)   :: XEMIT_LW_GRND     ! LW flux emitted by the ground (W/m2 ground = road + garden)
  REAL, POINTER, DIMENSION(:)   :: XT_RAD_IND     !Indoor mean radiant temperature [K]
  REAL, POINTER, DIMENSION(:)   :: XREF_SW_GRND ! total solar rad reflected by ground
  REAL, POINTER, DIMENSION(:)   :: XREF_SW_FAC ! total solar rad reflected by facade
  REAL, POINTER, DIMENSION(:)   :: XHU_BLD !Indoor relative humidity
!
  REAL, POINTER, DIMENSION(:)   :: XTCOOL_TARGET ! current cooling setpoint of indoor air
  REAL, POINTER, DIMENSION(:)   :: XTHEAT_TARGET ! current heating setpoint of indoor air
  REAL, POINTER, DIMENSION(:)   :: XQIN          ! current internal heat gains [W m-2(floor)]

!------------------------------------------------------------------------------
!

END TYPE DIAG_MISC_TEB_1P_t

TYPE DIAG_MISC_TEB_t
  !
  TYPE(DIAG_MISC_TEB_1P_t), POINTER :: ALP(:) => NULL()
  TYPE(DIAG_MISC_TEB_1P_t), POINTER :: CUR => NULL()
  !
END TYPE DIAG_MISC_TEB_t
!
!
!
CONTAINS


!




SUBROUTINE DIAG_MISC_TEB_GOTO_PATCH(YDIAG_MISC_TEB,KTO_PATCH)
TYPE(DIAG_MISC_TEB_t), INTENT(INOUT) :: YDIAG_MISC_TEB
INTEGER, INTENT(IN) :: KTO_PATCH
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('MODD_DIAG_MISC_TEB_N:DIAG_MISC_TEB_GOTO_PATCH',0,ZHOOK_HANDLE)

YDIAG_MISC_TEB%CUR => YDIAG_MISC_TEB%ALP(KTO_PATCH)

IF (LHOOK) CALL DR_HOOK('MODD_DIAG_MISC_TEB_N:DIAG_MISC_TEB_GOTO_PATCH',1,ZHOOK_HANDLE)

END SUBROUTINE DIAG_MISC_TEB_GOTO_PATCH

SUBROUTINE DIAG_MISC_TEB_1P_INIT(YDIAG_MISC_TEB_1P)
TYPE(DIAG_MISC_TEB_1P_t), INTENT(INOUT) :: YDIAG_MISC_TEB_1P
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_MISC_TEB_N:DIAG_MISC_TEB_INIT_1P",0,ZHOOK_HANDLE)
  NULLIFY(YDIAG_MISC_TEB_1P%XZ0_TOWN)
  NULLIFY(YDIAG_MISC_TEB_1P%XQF_BLD)
  NULLIFY(YDIAG_MISC_TEB_1P%XFLX_BLD)
  NULLIFY(YDIAG_MISC_TEB_1P%XQF_TOWN)
  NULLIFY(YDIAG_MISC_TEB_1P%XDQS_TOWN)
  NULLIFY(YDIAG_MISC_TEB_1P%XH_WALL_A)
  NULLIFY(YDIAG_MISC_TEB_1P%XH_WALL_B)
  NULLIFY(YDIAG_MISC_TEB_1P%XH_ROOF)
  NULLIFY(YDIAG_MISC_TEB_1P%XH_ROAD)
  NULLIFY(YDIAG_MISC_TEB_1P%XH_STRLROOF)
  NULLIFY(YDIAG_MISC_TEB_1P%XH_BLT)
  NULLIFY(YDIAG_MISC_TEB_1P%XRN_WALL_A)
  NULLIFY(YDIAG_MISC_TEB_1P%XRN_WALL_B)
  NULLIFY(YDIAG_MISC_TEB_1P%XRN_ROOF)
  NULLIFY(YDIAG_MISC_TEB_1P%XRN_ROAD)
  NULLIFY(YDIAG_MISC_TEB_1P%XRN_STRLROOF)
  NULLIFY(YDIAG_MISC_TEB_1P%XRN_BLT)
  NULLIFY(YDIAG_MISC_TEB_1P%XGFLUX_WALL_A)
  NULLIFY(YDIAG_MISC_TEB_1P%XGFLUX_WALL_B)
  NULLIFY(YDIAG_MISC_TEB_1P%XGFLUX_ROOF)
  NULLIFY(YDIAG_MISC_TEB_1P%XGFLUX_ROAD)
  NULLIFY(YDIAG_MISC_TEB_1P%XGFLUX_STRLROOF)
  NULLIFY(YDIAG_MISC_TEB_1P%XGFLUX_BLT)
  NULLIFY(YDIAG_MISC_TEB_1P%XLE_ROOF)
  NULLIFY(YDIAG_MISC_TEB_1P%XLE_ROAD)
  NULLIFY(YDIAG_MISC_TEB_1P%XLE_STRLROOF)
  NULLIFY(YDIAG_MISC_TEB_1P%XLE_BLT)
  NULLIFY(YDIAG_MISC_TEB_1P%XRUNOFF_TOWN)
  NULLIFY(YDIAG_MISC_TEB_1P%XIRRIG_GARDEN)
  NULLIFY(YDIAG_MISC_TEB_1P%XRUNOFF_ROAD)
  NULLIFY(YDIAG_MISC_TEB_1P%XIRRIG_ROAD)
  NULLIFY(YDIAG_MISC_TEB_1P%XRUNOFF_ROOF)
  NULLIFY(YDIAG_MISC_TEB_1P%XRUNOFF_STRLROOF)
  NULLIFY(YDIAG_MISC_TEB_1P%XIRRIG_GREENROOF)
  NULLIFY(YDIAG_MISC_TEB_1P%XABS_LW_GREENROOF)
  NULLIFY(YDIAG_MISC_TEB_1P%XG_GREENROOF_ROOF)
  NULLIFY(YDIAG_MISC_TEB_1P%XABS_SW_PANEL)
  NULLIFY(YDIAG_MISC_TEB_1P%XABS_LW_PANEL)
  NULLIFY(YDIAG_MISC_TEB_1P%XRN_PANEL)
  NULLIFY(YDIAG_MISC_TEB_1P%XH_PANEL)
  NULLIFY(YDIAG_MISC_TEB_1P%XTHER_PROD_PANEL)
  NULLIFY(YDIAG_MISC_TEB_1P%XPHOT_PROD_PANEL)
  NULLIFY(YDIAG_MISC_TEB_1P%XPROD_PANEL)
  NULLIFY(YDIAG_MISC_TEB_1P%XTHER_PROD_BLD)
  NULLIFY(YDIAG_MISC_TEB_1P%XPHOT_PROD_BLD)
  NULLIFY(YDIAG_MISC_TEB_1P%XT_BLD_COOL)
  NULLIFY(YDIAG_MISC_TEB_1P%XH_BLD_COOL)
  NULLIFY(YDIAG_MISC_TEB_1P%XH_BLD_HEAT)
  NULLIFY(YDIAG_MISC_TEB_1P%XLE_BLD_COOL)
  NULLIFY(YDIAG_MISC_TEB_1P%XLE_BLD_HEAT)
  NULLIFY(YDIAG_MISC_TEB_1P%XHVAC_COOL)
  NULLIFY(YDIAG_MISC_TEB_1P%XHVAC_HEAT)
  NULLIFY(YDIAG_MISC_TEB_1P%XH_WASTE)
  NULLIFY(YDIAG_MISC_TEB_1P%XLE_WASTE)
  NULLIFY(YDIAG_MISC_TEB_1P%XCAP_SYS)
  NULLIFY(YDIAG_MISC_TEB_1P%XM_SYS)
  NULLIFY(YDIAG_MISC_TEB_1P%XCOP)
  NULLIFY(YDIAG_MISC_TEB_1P%XQ_SYS)
  NULLIFY(YDIAG_MISC_TEB_1P%XT_SYS)
  NULLIFY(YDIAG_MISC_TEB_1P%XTR_SW_WIN)
  NULLIFY(YDIAG_MISC_TEB_1P%XFAN_POWER)
  NULLIFY(YDIAG_MISC_TEB_1P%XABS_SW_WIN)
  NULLIFY(YDIAG_MISC_TEB_1P%XABS_LW_WIN)
  NULLIFY(YDIAG_MISC_TEB_1P%XEMIT_LW_GRND)
  NULLIFY(YDIAG_MISC_TEB_1P%XEMIT_LW_FAC)
  NULLIFY(YDIAG_MISC_TEB_1P%XT_RAD_IND)
  NULLIFY(YDIAG_MISC_TEB_1P%XREF_SW_GRND)
  NULLIFY(YDIAG_MISC_TEB_1P%XREF_SW_FAC)
  NULLIFY(YDIAG_MISC_TEB_1P%XHU_BLD)
  NULLIFY(YDIAG_MISC_TEB_1P%XTCOOL_TARGET)
  NULLIFY(YDIAG_MISC_TEB_1P%XTHEAT_TARGET)
  NULLIFY(YDIAG_MISC_TEB_1P%XQIN)
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_MISC_TEB_N:DIAG_MISC_TEB_INIT_1P",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_MISC_TEB_1P_INIT
!
SUBROUTINE DIAG_MISC_TEB_INIT(YDIAG_MISC_TEB,KPATCH)
TYPE(DIAG_MISC_TEB_t), INTENT(INOUT) :: YDIAG_MISC_TEB
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_MISC_TEB_N:DIAG_MISC_TEB_INIT",0,ZHOOK_HANDLE)
IF (.NOT.ASSOCIATED(YDIAG_MISC_TEB%ALP)) THEN
  ALLOCATE(YDIAG_MISC_TEB%ALP(KPATCH))
  YDIAG_MISC_TEB%CUR => YDIAG_MISC_TEB%ALP(1)
  DO JP=1,KPATCH
    CALL DIAG_MISC_TEB_1P_INIT(YDIAG_MISC_TEB%ALP(JP))
  ENDDO
ELSE
  DO JP=1,KPATCH
    CALL DIAG_MISC_TEB_1P_INIT(YDIAG_MISC_TEB%ALP(JP))
  ENDDO
  YDIAG_MISC_TEB%CUR => NULL()
  DEALLOCATE(YDIAG_MISC_TEB%ALP)
ENDIF  
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_MISC_TEB_N:DIAG_MISC_TEB_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_MISC_TEB_INIT



END MODULE MODD_DIAG_MISC_TEB_n
