!#########
SUBROUTINE SFX_OASIS_SEND_OL(HPROGRAM,KI,PTIMEC,PSTEP_SURF,KSIZE_OMP)
!###########################################
!
!!****  *SFX_OASIS_SEND_OL* - Offline driver to send coupling fields
!!
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
!!	B. Decharme   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2013
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURFEX_OMP, ONLY :  NINDX1SFX, NINDX2SFX, NBLOCK, NBLOCKTOT, &
                             INIT_DIM, RESET_DIM
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE MODN_SFX_OASIS,  ONLY : XTSTEP_CPL_LAND, &
                            XTSTEP_CPL_LAKE, &
                            XTSTEP_CPL_SEA , &
                            LWATER
!
USE MODD_SFX_OASIS,  ONLY : LCPL_LAND,LCPL_GW,       &
                            LCPL_FLOOD,LCPL_CALVING, &
                            LCPL_LAKE,               &
                            LCPL_SEA,LCPL_SEAICE
!
USE MODI_GOTO_SURFEX
USE MODI_GET_SFX_LAND
USE MODI_GET_SFX_LAKE
USE MODI_GET_SFX_SEA
!
USE MODI_GET_LUOUT
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
#ifdef AIX64
!$ USE OMP_LIB
#endif
!
IMPLICIT NONE
!
#ifndef AIX64
!$ INCLUDE 'omp_lib.h'
#endif
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=*),      INTENT(IN) :: HPROGRAM
INTEGER,               INTENT(IN) :: KI            ! number of points
REAL,                  INTENT(IN) :: PTIMEC        ! Cumulated run time step (s)
REAL,                  INTENT(IN) :: PSTEP_SURF    ! Model time step (s)
INTEGER, DIMENSION(:), INTENT(IN) :: KSIZE_OMP
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL, DIMENSION(KI)   :: ZLAND_RUNOFF    ! Cumulated Surface runoff             (kg/m2)
REAL, DIMENSION(KI)   :: ZLAND_DRAIN     ! Cumulated Deep drainage              (kg/m2)
REAL, DIMENSION(KI)   :: ZLAND_CALVING   ! Cumulated Calving flux               (kg/m2)
REAL, DIMENSION(KI)   :: ZLAND_RECHARGE  ! Cumulated Recharge to groundwater    (kg/m2)
REAL, DIMENSION(KI)   :: ZLAND_PFLOOD    ! Cumulated flood precip interception  (kg/m2)
REAL, DIMENSION(KI)   :: ZLAND_EFLOOD    ! Cumulated floodplains evaporation    (kg/m2)
REAL, DIMENSION(KI)   :: ZLAND_IFLOOD    ! Cumulated floodplains infiltration   (kg/m2)
!
REAL, DIMENSION(KI)   :: ZLAKE_EVAP  ! Cumulated Evaporation             (kg/m2)
REAL, DIMENSION(KI)   :: ZLAKE_RAIN  ! Cumulated Rainfall rate           (kg/m2)
REAL, DIMENSION(KI)   :: ZLAKE_SNOW  ! Cumulated Snowfall rate           (kg/m2)
!
REAL, DIMENSION(KI)   :: ZSEA_FWSU  ! Cumulated zonal wind stress       (Pa.s)
REAL, DIMENSION(KI)   :: ZSEA_FWSV  ! Cumulated meridian wind stress    (Pa.s)
REAL, DIMENSION(KI)   :: ZSEA_HEAT  ! Cumulated Non solar net heat flux (J/m2)
REAL, DIMENSION(KI)   :: ZSEA_SNET  ! Cumulated Solar net heat flux     (J/m2)
REAL, DIMENSION(KI)   :: ZSEA_WIND  ! Cumulated 10m wind speed          (m)
REAL, DIMENSION(KI)   :: ZSEA_FWSM  ! Cumulated wind stress             (Pa.s)
REAL, DIMENSION(KI)   :: ZSEA_EVAP  ! Cumulated Evaporation             (kg/m2)
REAL, DIMENSION(KI)   :: ZSEA_RAIN  ! Cumulated Rainfall rate           (kg/m2)
REAL, DIMENSION(KI)   :: ZSEA_SNOW  ! Cumulated Snowfall rate           (kg/m2)
REAL, DIMENSION(KI)   :: ZSEA_WATF  ! Cumulated net freshwater rate     (kg/m2)
!
REAL, DIMENSION(KI)   :: ZSEAICE_HEAT ! Cumulated Sea-ice non solar net heat flux (J/m2)
REAL, DIMENSION(KI)   :: ZSEAICE_SNET ! Cumulated Sea-ice solar net heat flux     (J/m2)
REAL, DIMENSION(KI)   :: ZSEAICE_EVAP ! Cumulated Sea-ice sublimation             (kg/m2)
!
INTEGER               :: IDATE  ! current coupling time step (s)
INTEGER               :: ILUOUT
INTEGER               :: INKPROMA
!
LOGICAL               :: GSEND_LAND
LOGICAL               :: GSEND_LAKE
LOGICAL               :: GSEND_SEA
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SFX_OASIS_SEND_OL',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!*       1.     Initialize proc by proc :
!               -------------------------
!
CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
IDATE = INT(PTIMEC-PSTEP_SURF)
!
GSEND_LAND=(LCPL_LAND.AND.MOD(PTIMEC,XTSTEP_CPL_LAND)==0.0)
GSEND_LAKE=(LCPL_LAKE.AND.MOD(PTIMEC,XTSTEP_CPL_LAKE)==0.0)
GSEND_SEA =(LCPL_SEA .AND.MOD(PTIMEC,XTSTEP_CPL_SEA )==0.0)
!
IF(GSEND_LAND)THEN
  ZLAND_RUNOFF  (:) = XUNDEF
  ZLAND_DRAIN   (:) = XUNDEF
  ZLAND_CALVING (:) = XUNDEF
  ZLAND_RECHARGE(:) = XUNDEF
  ZLAND_PFLOOD  (:) = XUNDEF
  ZLAND_EFLOOD  (:) = XUNDEF
  ZLAND_IFLOOD  (:) = XUNDEF
ENDIF
!
IF(GSEND_LAKE)THEN
  ZLAKE_EVAP (:) = XUNDEF
  ZLAKE_RAIN (:) = XUNDEF
  ZLAKE_SNOW (:) = XUNDEF
ENDIF
!
IF(GSEND_SEA)THEN
  ZSEA_FWSU (:) = XUNDEF
  ZSEA_FWSV (:) = XUNDEF
  ZSEA_HEAT (:) = XUNDEF
  ZSEA_SNET (:) = XUNDEF
  ZSEA_WIND (:) = XUNDEF
  ZSEA_FWSM (:) = XUNDEF
  ZSEA_EVAP (:) = XUNDEF
  ZSEA_RAIN (:) = XUNDEF
  ZSEA_SNOW (:) = XUNDEF
  ZSEA_WATF (:) = XUNDEF
  !
  ZSEAICE_HEAT (:) = XUNDEF
  ZSEAICE_SNET (:) = XUNDEF
  ZSEAICE_EVAP (:) = XUNDEF
ENDIF
!
!-------------------------------------------------------------------------------
!
!$OMP PARALLEL PRIVATE(INKPROMA)
!
!$ NBLOCK = OMP_GET_THREAD_NUM()
!
IF (NBLOCK==NBLOCKTOT) THEN
   CALL INIT_DIM(KSIZE_OMP,0,INKPROMA,NINDX1SFX,NINDX2SFX)
ELSE
   CALL INIT_DIM(KSIZE_OMP,NBLOCK,INKPROMA,NINDX1SFX,NINDX2SFX)
ENDIF
!
IF (NBLOCK==0) THEN
   CALL GOTO_SURFEX(NBLOCKTOT,.TRUE.)
ELSE
   CALL GOTO_SURFEX(NBLOCK,.TRUE.)
ENDIF
!
!*       2.     get local fields :
!               ------------------
!
IF(GSEND_LAND)THEN
!
! * Get river output fields
!
  CALL GET_SFX_LAND(LCPL_GW,LCPL_FLOOD,LCPL_CALVING,                           &
                    ZLAND_RUNOFF (NINDX1SFX:NINDX2SFX),ZLAND_DRAIN   (NINDX1SFX:NINDX2SFX),&
                    ZLAND_CALVING(NINDX1SFX:NINDX2SFX),ZLAND_RECHARGE(NINDX1SFX:NINDX2SFX),&
                    ZLAND_PFLOOD (NINDX1SFX:NINDX2SFX),ZLAND_EFLOOD  (NINDX1SFX:NINDX2SFX),&
                    ZLAND_IFLOOD (NINDX1SFX:NINDX2SFX)                               )
!
ENDIF
!
IF(GSEND_LAKE)THEN
!
! * Get output fields
!
  CALL GET_SFX_LAKE(ZLAKE_EVAP(NINDX1SFX:NINDX2SFX),ZLAKE_RAIN(NINDX1SFX:NINDX2SFX), &
                    ZLAKE_SNOW(NINDX1SFX:NINDX2SFX)                            )
!
ENDIF
!
IF(GSEND_SEA)THEN
!
! * Get sea output fields
!
  CALL GET_SFX_SEA(LCPL_SEAICE,LWATER,                                                                                    &
                   ZSEA_FWSU   (NINDX1SFX:NINDX2SFX),ZSEA_FWSV   (NINDX1SFX:NINDX2SFX),ZSEA_HEAT   (NINDX1SFX:NINDX2SFX),&
                   ZSEA_SNET   (NINDX1SFX:NINDX2SFX),ZSEA_WIND   (NINDX1SFX:NINDX2SFX),ZSEA_FWSM   (NINDX1SFX:NINDX2SFX),&
                   ZSEA_EVAP   (NINDX1SFX:NINDX2SFX),ZSEA_RAIN   (NINDX1SFX:NINDX2SFX),ZSEA_SNOW   (NINDX1SFX:NINDX2SFX),&
                   ZSEA_WATF   (NINDX1SFX:NINDX2SFX),                                                                    &
                   ZSEAICE_HEAT(NINDX1SFX:NINDX2SFX),ZSEAICE_SNET(NINDX1SFX:NINDX2SFX),ZSEAICE_EVAP(NINDX1SFX:NINDX2SFX) )
!
ENDIF
!
CALL RESET_DIM(KI,INKPROMA,NINDX1SFX,NINDX2SFX)
!
!$OMP END PARALLEL
!    
!-------------------------------------------------------------------------------
!
!*       3.     Send fields to OASIS proc by proc:
!               ----------------------------------
!
IF(GSEND_LAND.OR.GSEND_LAKE.OR.GSEND_SEA)THEN
!
  CALL SFX_OASIS_SEND(ILUOUT,KI,IDATE,GSEND_LAND,GSEND_LAKE,GSEND_SEA,      &
                      ZLAND_RUNOFF,ZLAND_DRAIN,ZLAND_CALVING,ZLAND_RECHARGE,&
                      ZLAND_PFLOOD,ZLAND_EFLOOD,ZLAND_IFLOOD,               &
                      ZLAKE_EVAP,ZLAKE_RAIN,ZLAKE_SNOW,                     &
                      ZSEA_FWSU,ZSEA_FWSV,ZSEA_HEAT,ZSEA_SNET,ZSEA_WIND,    &
                      ZSEA_FWSM,ZSEA_EVAP,ZSEA_RAIN,ZSEA_SNOW,ZSEA_WATF,    &
                      ZSEAICE_HEAT,ZSEAICE_SNET,ZSEAICE_EVAP                )
                      
!
ENDIF
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SFX_OASIS_SEND_OL',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SFX_OASIS_SEND_OL
