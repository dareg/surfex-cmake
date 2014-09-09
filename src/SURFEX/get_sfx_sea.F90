!     #########
      SUBROUTINE GET_SFX_SEA(OCPL_SEAICE,OWATER,                       &
                              PSEA_FWSU,PSEA_FWSV,PSEA_HEAT,PSEA_SNET, &
                              PSEA_WIND,PSEA_FWSM,PSEA_EVAP,PSEA_RAIN, &
                              PSEA_SNOW,PSEA_WATF,                     &
                              PSEAICE_HEAT,PSEAICE_SNET,PSEAICE_EVAP   )  
!     ############################################################################
!
!!****  *GET_SFX_SEA* - routine to get some variables from surfex to
!                        a oceanic general circulation model
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	B. Decharme      *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2013
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE MODD_SURF_ATM_n, ONLY : NR_SEA, NSIZE_SEA, XSEA, &
                            NR_WATER, NSIZE_WATER,   &
                            XWATER, CWATER
!
USE MODD_SEAFLUX_n, ONLY : XCPL_SEA_WIND,XCPL_SEA_EVAP,XCPL_SEA_HEAT, &
                           XCPL_SEA_SNET,XCPL_SEA_FWSU,XCPL_SEA_FWSV, &
                           XCPL_SEA_RAIN,XCPL_SEA_SNOW,XCPL_SEA_FWSM, &
                           XCPL_SEAICE_EVAP,XCPL_SEAICE_HEAT,         &
                           XCPL_SEAICE_SNET  
!
USE MODD_WATFLUX_n, ONLY : XCPL_WATER_WIND,XCPL_WATER_EVAP,XCPL_WATER_HEAT, &
                           XCPL_WATER_SNET,XCPL_WATER_FWSU,XCPL_WATER_FWSV, &
                           XCPL_WATER_RAIN,XCPL_WATER_SNOW,XCPL_WATER_FWSM, &
                           XCPL_WATERICE_EVAP,XCPL_WATERICE_HEAT,           &
                           XCPL_WATERICE_SNET  
!
USE MODI_UNPACK_SAME_RANK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
LOGICAL,            INTENT(IN)  :: OCPL_SEAICE ! sea-ice / ocean key
LOGICAL,            INTENT(IN)  :: OWATER      ! water included in sea smask
!
REAL, DIMENSION(:), INTENT(OUT) :: PSEA_FWSU  ! Cumulated zonal wind stress       (Pa.s)
REAL, DIMENSION(:), INTENT(OUT) :: PSEA_FWSV  ! Cumulated meridian wind stress    (Pa.s)
REAL, DIMENSION(:), INTENT(OUT) :: PSEA_HEAT  ! Cumulated Non solar net heat flux (J/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PSEA_SNET  ! Cumulated Solar net heat flux     (J/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PSEA_WIND  ! Cumulated 10m wind speed          (m)
REAL, DIMENSION(:), INTENT(OUT) :: PSEA_FWSM  ! Cumulated wind stress             (Pa.s)
REAL, DIMENSION(:), INTENT(OUT) :: PSEA_EVAP  ! Cumulated Evaporation             (kg/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PSEA_RAIN  ! Cumulated Rainfall rate           (kg/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PSEA_SNOW  ! Cumulated Snowfall rate           (kg/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PSEA_WATF  ! Cumulated Net outgoing water flux (kg/m2)
!
REAL, DIMENSION(:), INTENT(OUT) :: PSEAICE_HEAT ! Cumulated Sea-ice non solar net heat flux (J/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PSEAICE_SNET ! Cumulated Sea-ice solar net heat flux     (J/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PSEAICE_EVAP ! Cumulated Sea-ice sublimation             (kg/m2)
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL, DIMENSION(SIZE(PSEA_HEAT))  :: ZWIND
REAL, DIMENSION(SIZE(PSEA_HEAT))  :: ZFWSU
REAL, DIMENSION(SIZE(PSEA_HEAT))  :: ZFWSV
REAL, DIMENSION(SIZE(PSEA_HEAT))  :: ZSNET
REAL, DIMENSION(SIZE(PSEA_HEAT))  :: ZHEAT
REAL, DIMENSION(SIZE(PSEA_HEAT))  :: ZEVAP
REAL, DIMENSION(SIZE(PSEA_HEAT))  :: ZRAIN
REAL, DIMENSION(SIZE(PSEA_HEAT))  :: ZSNOW
REAL, DIMENSION(SIZE(PSEA_HEAT))  :: ZFWSM
!
REAL, DIMENSION(SIZE(PSEA_HEAT))  :: ZSNET_ICE
REAL, DIMENSION(SIZE(PSEA_HEAT))  :: ZHEAT_ICE
REAL, DIMENSION(SIZE(PSEA_HEAT))  :: ZEVAP_ICE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GET_SFX_SEA',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!*       1.0   Initialization
!              --------------
!
PSEA_FWSU (:) = XUNDEF
PSEA_FWSV (:) = XUNDEF
PSEA_HEAT (:) = XUNDEF
PSEA_SNET (:) = XUNDEF
PSEA_WIND (:) = XUNDEF
PSEA_FWSM (:) = XUNDEF
PSEA_EVAP (:) = XUNDEF
PSEA_RAIN (:) = XUNDEF
PSEA_SNOW (:) = XUNDEF
PSEA_WATF (:) = XUNDEF
!
PSEAICE_HEAT (:) = XUNDEF
PSEAICE_SNET (:) = XUNDEF
PSEAICE_EVAP (:) = XUNDEF
!
ZFWSU (:) = XUNDEF
ZFWSV (:) = XUNDEF
ZHEAT (:) = XUNDEF
ZSNET (:) = XUNDEF
ZWIND (:) = XUNDEF
ZFWSM (:) = XUNDEF
ZEVAP (:) = XUNDEF
ZRAIN (:) = XUNDEF
ZSNOW (:) = XUNDEF
!
ZHEAT_ICE (:) = XUNDEF
ZSNET_ICE (:) = XUNDEF
ZEVAP_ICE (:) = XUNDEF
!
!*       2.0   Get variable over sea
!              ---------------------
!
IF(NSIZE_SEA>0)THEN
!
  CALL UNPACK_SAME_RANK(NR_SEA(:),XCPL_SEA_WIND(:),PSEA_WIND(:),XUNDEF)
  CALL UNPACK_SAME_RANK(NR_SEA(:),XCPL_SEA_FWSU(:),PSEA_FWSU(:),XUNDEF)
  CALL UNPACK_SAME_RANK(NR_SEA(:),XCPL_SEA_FWSV(:),PSEA_FWSV(:),XUNDEF)
  CALL UNPACK_SAME_RANK(NR_SEA(:),XCPL_SEA_SNET(:),PSEA_SNET(:),XUNDEF)
  CALL UNPACK_SAME_RANK(NR_SEA(:),XCPL_SEA_HEAT(:),PSEA_HEAT(:),XUNDEF)
  CALL UNPACK_SAME_RANK(NR_SEA(:),XCPL_SEA_EVAP(:),PSEA_EVAP(:),XUNDEF)
  CALL UNPACK_SAME_RANK(NR_SEA(:),XCPL_SEA_RAIN(:),PSEA_RAIN(:),XUNDEF)
  CALL UNPACK_SAME_RANK(NR_SEA(:),XCPL_SEA_SNOW(:),PSEA_SNOW(:),XUNDEF)
  CALL UNPACK_SAME_RANK(NR_SEA(:),XCPL_SEA_FWSM(:),PSEA_FWSM(:),XUNDEF)
  XCPL_SEA_WIND(:) = 0.0
  XCPL_SEA_EVAP(:) = 0.0
  XCPL_SEA_HEAT(:) = 0.0
  XCPL_SEA_SNET(:) = 0.0
  XCPL_SEA_FWSU(:) = 0.0
  XCPL_SEA_FWSV(:) = 0.0
  XCPL_SEA_RAIN(:) = 0.0
  XCPL_SEA_SNOW(:) = 0.0
  XCPL_SEA_FWSM(:) = 0.0
!
  IF (OCPL_SEAICE) THEN
    CALL UNPACK_SAME_RANK(NR_SEA(:),XCPL_SEAICE_SNET(:),PSEAICE_SNET(:),XUNDEF)
    CALL UNPACK_SAME_RANK(NR_SEA(:),XCPL_SEAICE_HEAT(:),PSEAICE_HEAT(:),XUNDEF)
    CALL UNPACK_SAME_RANK(NR_SEA(:),XCPL_SEAICE_EVAP(:),PSEAICE_EVAP(:),XUNDEF)
    XCPL_SEAICE_SNET(:) = 0.0
    XCPL_SEAICE_EVAP(:) = 0.0
    XCPL_SEAICE_HEAT(:) = 0.0  
  ENDIF
!
ENDIF
!
!*       3.0   Get variable over water without Flake
!              -------------------------------------
!
IF (OWATER.AND.NSIZE_WATER>0) THEN
!
  CALL UNPACK_SAME_RANK(NR_WATER(:),XCPL_WATER_WIND(:),ZWIND(:),XUNDEF)
  CALL UNPACK_SAME_RANK(NR_WATER(:),XCPL_WATER_FWSU(:),ZFWSU(:),XUNDEF)
  CALL UNPACK_SAME_RANK(NR_WATER(:),XCPL_WATER_FWSV(:),ZFWSV(:),XUNDEF)
  CALL UNPACK_SAME_RANK(NR_WATER(:),XCPL_WATER_SNET(:),ZSNET(:),XUNDEF)
  CALL UNPACK_SAME_RANK(NR_WATER(:),XCPL_WATER_HEAT(:),ZHEAT(:),XUNDEF)
  CALL UNPACK_SAME_RANK(NR_WATER(:),XCPL_WATER_EVAP(:),ZEVAP(:),XUNDEF)
  CALL UNPACK_SAME_RANK(NR_WATER(:),XCPL_WATER_RAIN(:),ZRAIN(:),XUNDEF)
  CALL UNPACK_SAME_RANK(NR_WATER(:),XCPL_WATER_SNOW(:),ZSNOW(:),XUNDEF)
  CALL UNPACK_SAME_RANK(NR_WATER(:),XCPL_WATER_FWSM(:),ZFWSM(:),XUNDEF)
!
  WHERE(XWATER(:)>0.0) 
    PSEA_WIND(:) = (XSEA(:)*PSEA_WIND(:)+XWATER(:)*ZWIND(:))/(XSEA(:)+XWATER(:))
    PSEA_FWSU(:) = (XSEA(:)*PSEA_FWSU(:)+XWATER(:)*ZFWSU(:))/(XSEA(:)+XWATER(:))
    PSEA_FWSV(:) = (XSEA(:)*PSEA_FWSV(:)+XWATER(:)*ZFWSV(:))/(XSEA(:)+XWATER(:))
    PSEA_SNET(:) = (XSEA(:)*PSEA_SNET(:)+XWATER(:)*ZSNET(:))/(XSEA(:)+XWATER(:))
    PSEA_HEAT(:) = (XSEA(:)*PSEA_HEAT(:)+XWATER(:)*ZHEAT(:))/(XSEA(:)+XWATER(:))
    PSEA_EVAP(:) = (XSEA(:)*PSEA_EVAP(:)+XWATER(:)*ZEVAP(:))/(XSEA(:)+XWATER(:))
    PSEA_RAIN(:) = (XSEA(:)*PSEA_RAIN(:)+XWATER(:)*ZRAIN(:))/(XSEA(:)+XWATER(:))
    PSEA_SNOW(:) = (XSEA(:)*PSEA_SNOW(:)+XWATER(:)*ZSNOW(:))/(XSEA(:)+XWATER(:))
    PSEA_FWSM(:) = (XSEA(:)*PSEA_FWSM(:)+XWATER(:)*ZFWSM(:))/(XSEA(:)+XWATER(:))
  ENDWHERE 
!
  XCPL_WATER_WIND(:) = 0.0
  XCPL_WATER_EVAP(:) = 0.0
  XCPL_WATER_HEAT(:) = 0.0
  XCPL_WATER_SNET(:) = 0.0
  XCPL_WATER_FWSU(:) = 0.0
  XCPL_WATER_FWSV(:) = 0.0
  XCPL_WATER_RAIN(:) = 0.0
  XCPL_WATER_SNOW(:) = 0.0
  XCPL_WATER_FWSM(:) = 0.0
!
  IF (OCPL_SEAICE) THEN
    CALL UNPACK_SAME_RANK(NR_WATER(:),XCPL_WATERICE_SNET(:),ZSNET_ICE(:),XUNDEF)
    CALL UNPACK_SAME_RANK(NR_WATER(:),XCPL_WATERICE_HEAT(:),ZHEAT_ICE(:),XUNDEF)
    CALL UNPACK_SAME_RANK(NR_WATER(:),XCPL_WATERICE_EVAP(:),ZEVAP_ICE(:),XUNDEF)
    WHERE(XWATER(:)>0.0)     
      PSEAICE_SNET(:) = (XSEA(:)*PSEAICE_SNET(:)+XWATER(:)*ZSNET_ICE(:))/(XSEA(:)+XWATER(:))
      PSEAICE_HEAT(:) = (XSEA(:)*PSEAICE_HEAT(:)+XWATER(:)*ZHEAT_ICE(:))/(XSEA(:)+XWATER(:))
      PSEAICE_EVAP(:) = (XSEA(:)*PSEAICE_EVAP(:)+XWATER(:)*ZEVAP_ICE(:))/(XSEA(:)+XWATER(:))
    ENDWHERE  
    XCPL_WATERICE_SNET(:) = 0.0
    XCPL_WATERICE_EVAP(:) = 0.0
    XCPL_WATERICE_HEAT(:) = 0.0
  ENDIF  
! 
ENDIF
!
!*       4.0   Net outgoing water flux
!              -----------------------
!
IF(NSIZE_SEA>0)THEN
!
  PSEA_WATF(:) = PSEA_EVAP(:) - PSEA_RAIN(:) - PSEA_SNOW(:)
!
ENDIF
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GET_SFX_SEA',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE GET_SFX_SEA
