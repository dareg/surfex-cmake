!     ################
      MODULE MODD_BEM_n
!     ################
!
!!****  *MODD_BEM_n - declaration of parameters and option for BEM
!!
!!    PURPOSE
!!    -------
!     Declaration of surface parameters
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
!!      B. Bueno   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       10/2010
!!      G. Pigeon      06/2011 add LSHAD_DAY
!!      G. Pigeon      07/2011 add LNATVENT_NIGHT
!!      G. Pigeon      08/2011 change from MODD_BLD -> MODD_BEM
!!      G. Pigeon      10/2011 add indoor relative surf. and view factors
!!      G. Pigeon      09/2012 add TRAN_WIN
!!      G. Pigeon      10/2012 add XF_WIN_WIN
!!      V. Masson      06/2013 splits module in two
!
!*       0.   DECLARATIONS
!             ------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!--------------------------------------------------------------------------
!
TYPE BEM_1P_t
!
! Floor parameters
!
  REAL, POINTER, DIMENSION(:,:) :: XHC_FLOOR     ! floor layers heat capacity        (J/K/m3)
  REAL, POINTER, DIMENSION(:,:) :: XTC_FLOOR     ! floor layers thermal conductivity (W/K/m)
  REAL, POINTER, DIMENSION(:,:) :: XD_FLOOR      ! depth of floor layers             (m)
!
! HVAC parameters
!
  REAL, POINTER, DIMENSION(:)   :: XTCOOL_TARGET ! cooling setpoint of indoor air
  REAL, POINTER, DIMENSION(:)   :: XTHEAT_TARGET ! heating setpoint of indoor air
  REAL, POINTER, DIMENSION(:)   :: XF_WASTE_CAN  ! fraction of waste heat released into the canyon
  REAL, POINTER, DIMENSION(:)   :: XEFF_HEAT     ! efficiency of the heating system
!
! Indoor parameters
!
  REAL, POINTER, DIMENSION(:)   :: XTI_BLD       ! building interior temperature    (K)
  REAL, POINTER, DIMENSION(:,:) :: XT_FLOOR      ! floor layer temperatures         (K)
  REAL, POINTER, DIMENSION(:,:) :: XT_MASS       ! Air cooled building internal th. mass temperature (K)
!
  REAL, POINTER, DIMENSION(:)   :: XQIN          ! internal heat gains [W m-2(floor)]
  REAL, POINTER, DIMENSION(:)   :: XQIN_FRAD     ! radiant fraction of internal heat gains
  REAL, POINTER, DIMENSION(:)   :: XSHGC         ! solar heat gain coef. of windows
  REAL, POINTER, DIMENSION(:)   :: XSHGC_SH      ! solar heat gain coef. of windows + shading
  REAL, POINTER, DIMENSION(:)   :: XU_WIN        ! window U-factor [K m W-2]
  REAL, POINTER, DIMENSION(:)   :: XTRAN_WIN     ! window transmittance (-)
  REAL, POINTER, DIMENSION(:)   :: XGR           ! glazing ratio
  REAL, POINTER, DIMENSION(:)   :: XFLOOR_HEIGHT ! building floor height [m]
  REAL, POINTER, DIMENSION(:)   :: XINF          ! infiltration/ventilation flow rate [AC/H]
!
! New parameters
!
  REAL, POINTER, DIMENSION(:)   :: XF_WATER_COND  ! fraction of evaporation for condensers (cooling system)
  REAL, POINTER, DIMENSION(:)   :: XAUX_MAX      ! Auxiliar variable for autosize calcs
  REAL, POINTER, DIMENSION(:)   :: XQIN_FLAT     ! Latent franction of internal heat gains
  REAL, POINTER, DIMENSION(:)   :: XHR_TARGET    ! Relative humidity setpoint
  REAL, POINTER, DIMENSION(:)   :: XT_WIN2       ! Indoor window temperature [K]
  REAL, POINTER, DIMENSION(:)   :: XQI_BLD       ! Indoor air specific humidity [kg kg-1]
  REAL, POINTER, DIMENSION(:)   :: XV_VENT       ! Ventilation flow rate [AC/H]
  REAL, POINTER, DIMENSION(:)   :: XCAP_SYS_HEAT ! Capacity of the heating system 
                                                 ! [W m-2(bld)]
  REAL, POINTER, DIMENSION(:)   :: XCAP_SYS_RAT  ! Rated capacity of the cooling system
                                                 ! [W m-2(bld)]
  REAL, POINTER, DIMENSION(:)   :: XT_ADP        ! Apparatus dewpoint temperature of the
                                                 ! cooling coil [K]
  REAL, POINTER, DIMENSION(:)   :: XM_SYS_RAT    ! Rated HVAC mass flow rate 
                                                 ! [kg s-1 m-2(bld)]
  REAL, POINTER, DIMENSION(:)   :: XCOP_RAT      ! Rated COP of the cooling system
  REAL, POINTER, DIMENSION(:)   :: XT_WIN1       ! outdoor window temperature [K]
  REAL, POINTER, DIMENSION(:)   :: XALB_WIN      ! window albedo
  REAL, POINTER, DIMENSION(:)   :: XABS_WIN      ! window absortance
  REAL, POINTER, DIMENSION(:)   :: XT_SIZE_MAX   ! Maximum outdoor air temperature for
                                                 ! HVAC sizing [K]
  REAL, POINTER, DIMENSION(:)   :: XT_SIZE_MIN   ! Minimum outdoor air temperature for
                                                 ! HVAC sizing [K]
  REAL, POINTER, DIMENSION(:)   :: XUGG_WIN      ! Window glass-to-glass U-factor [K m W-2]
  LOGICAL, POINTER, DIMENSION(:):: LSHADE        ! flag to activate shading devices -> LOGICAL in the code
  REAL,    POINTER, DIMENSION(:):: XSHADE        ! flag to activate shading devices -> REAL for i/o 0. or 1.
  CHARACTER(LEN=4), POINTER, DIMENSION(:) :: CNATVENT ! flag to activate natural ventilation 'NONE', 'MANU', 'AUTO'
  REAL,    POINTER, DIMENSION(:):: XNATVENT      ! flag to describe surventilation system for i/o 
                                                 ! 0 for NONE, 1 for MANU and 2 for AUTO
  LOGICAL, POINTER, DIMENSION(:):: LSHAD_DAY     !Has shading been necessary this day ?
  LOGICAL, POINTER, DIMENSION(:):: LNATVENT_NIGHT !Has nocturnal surventilation been necessary and possible this night ?
  !
  !indoor relative surfaces and view factors
  REAL, POINTER, DIMENSION(:) :: XN_FLOOR        ! Number of floors     
  REAL, POINTER, DIMENSION(:) :: XGLAZ_O_BLD    ! Window area [m2_win/m2_bld]
  REAL, POINTER, DIMENSION(:) :: XMASS_O_BLD    ! Mass area [m2_mass/m2_bld]
  REAL, POINTER, DIMENSION(:) :: XFLOOR_HW_RATIO ! H/W ratio of 1 floor level
  REAL, POINTER, DIMENSION(:) :: XF_FLOOR_MASS   ! View factor floor-mass
  REAL, POINTER, DIMENSION(:) :: XF_FLOOR_WALL   ! View factor floor-wall
  REAL, POINTER, DIMENSION(:) :: XF_FLOOR_WIN    ! View factor floor-window
  REAL, POINTER, DIMENSION(:) :: XF_FLOOR_ROOF   ! View factor floor-roof
  REAL, POINTER, DIMENSION(:) :: XF_WALL_FLOOR   ! View factor wall-floor
  REAL, POINTER, DIMENSION(:) :: XF_WALL_MASS    ! View factor wall-mass
  REAL, POINTER, DIMENSION(:) :: XF_WALL_WIN     ! View factor wall-win
  REAL, POINTER, DIMENSION(:) :: XF_WIN_FLOOR    ! View factor win-floor
  REAL, POINTER, DIMENSION(:) :: XF_WIN_MASS     ! View factor win-mass
  REAL, POINTER, DIMENSION(:) :: XF_WIN_WALL     ! View factor win-wall
  REAL, POINTER, DIMENSION(:) :: XF_WIN_WIN      ! indoor View factor win-win
  REAL, POINTER, DIMENSION(:) :: XF_MASS_FLOOR   ! View factor mass-floor
  REAL, POINTER, DIMENSION(:) :: XF_MASS_WALL    ! View factor mass-wall
  REAL, POINTER, DIMENSION(:) :: XF_MASS_WIN     ! View factor mass-window


! 
END TYPE BEM_1P_t
!
TYPE BEM_t
  !
  TYPE(BEM_1P_t), POINTER :: ALP(:) => NULL()
  TYPE(BEM_1P_t), POINTER :: CUR => NULL()
  !
END TYPE BEM_t
!
TYPE(BEM_t), ALLOCATABLE, TARGET, SAVE :: BEM_MODEL(:)

TYPE(BEM_t), POINTER :: BEM => NULL()
!$OMP THREADPRIVATE(BEM)

CONTAINS

!----------------------------------------------------------------------------

SUBROUTINE BEM_GOTO_MODEL(KFROM, KTO, LKFROM)
INTEGER, INTENT(IN) :: KFROM, KTO
LOGICAL, INTENT(IN) :: LKFROM
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_BEM_N:BEM_GOTO_MODEL',0,ZHOOK_HANDLE)

BEM => BEM_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_BEM_N:BEM_GOTO_MODEL',1,ZHOOK_HANDLE)
!
END SUBROUTINE BEM_GOTO_MODEL
!
SUBROUTINE BEM_GOTO_PATCH(KTO_PATCH)
INTEGER, INTENT(IN) :: KTO_PATCH
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current patch is set to patch KTO_PATCH
IF (LHOOK) CALL DR_HOOK('MODD_BEM_N:BEM_GOTO_PATCH',0,ZHOOK_HANDLE)

BEM%CUR => BEM%ALP(KTO_PATCH)

IF (LHOOK) CALL DR_HOOK('MODD_BEM_N:BEM_GOTO_PATCH',1,ZHOOK_HANDLE)

END SUBROUTINE BEM_GOTO_PATCH
!
!
SUBROUTINE BEM_ALLOC(KMODEL,KPATCH)
INTEGER, INTENT(IN) :: KMODEL
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: J,JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_BEM_N:BEM_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(BEM_MODEL(KMODEL))
DO J=1,KMODEL
 ALLOCATE(BEM_MODEL(J)%ALP(KPATCH))
 DO JP=1,KPATCH
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XF_WATER_COND)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XHC_FLOOR)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XTC_FLOOR)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XD_FLOOR)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XTCOOL_TARGET)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XTHEAT_TARGET)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XTI_BLD)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XT_FLOOR)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XT_MASS)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XQIN)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XQIN_FRAD)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XSHGC)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XSHGC_SH)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XU_WIN)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XTRAN_WIN)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XGR)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XFLOOR_HEIGHT)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XEFF_HEAT)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XINF)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XF_WASTE_CAN)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XAUX_MAX)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XQIN_FLAT)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XHR_TARGET)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XT_WIN2)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XQI_BLD)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XV_VENT)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XCAP_SYS_HEAT)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XCAP_SYS_RAT)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XT_ADP)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XM_SYS_RAT)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XCOP_RAT)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XT_WIN1)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XALB_WIN)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XABS_WIN)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XT_SIZE_MAX)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XT_SIZE_MIN)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XUGG_WIN)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%LSHAD_DAY)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%LNATVENT_NIGHT)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%LSHADE)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XSHADE)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%CNATVENT)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XNATVENT)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XN_FLOOR)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XGLAZ_O_BLD)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XMASS_O_BLD)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XFLOOR_HW_RATIO)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XF_FLOOR_MASS)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XF_FLOOR_WALL)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XF_FLOOR_WIN)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XF_FLOOR_ROOF)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XF_WALL_FLOOR)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XF_WALL_MASS)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XF_WALL_WIN)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XF_WIN_FLOOR)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XF_WIN_MASS)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XF_WIN_WALL)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XF_WIN_WIN)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XF_MASS_FLOOR)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XF_MASS_WALL)
  NULLIFY(BEM_MODEL(J)%ALP(JP)%XF_MASS_WIN)
 ENDDO
ENDDO
IF (LHOOK) CALL DR_HOOK("MODD_BEM_N:BEM_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE BEM_ALLOC
!
SUBROUTINE BEM_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_BEM_N:BEM_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(BEM_MODEL)) DEALLOCATE(BEM_MODEL)
IF (ASSOCIATED(BEM)) NULLIFY(BEM)
IF (LHOOK) CALL DR_HOOK("MODD_BEM_N:BEM_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE BEM_DEALLO
!
END MODULE MODD_BEM_n
