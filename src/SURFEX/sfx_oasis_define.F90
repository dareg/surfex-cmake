!#########
SUBROUTINE SFX_OASIS_DEFINE(HPROGRAM,KNPTS,KPARAL)
!###################################################
!
!!****  *SFX_OASIS_DEFINE* - Definitions for exchange of coupling fields
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
!!      B. Decharme   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2013
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_PAR,  ONLY : NUNDEF
!
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
USE MODD_ISBA_n, ONLY : I => ISBA
!
USE MODN_SFX_OASIS
USE MODD_SFX_OASIS
!
USE MODI_GET_LUOUT
USE MODI_ABOR1_SFX
USE MODI_SFX_OASIS_CHECK
!
#ifdef SFXOASIS
USE MOD_OASIS
#endif
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=6),        INTENT(IN) :: HPROGRAM    ! program calling surf. schemes
INTEGER,                 INTENT(IN) :: KNPTS  ! Number of grid point on this proc
INTEGER, DIMENSION(:),   INTENT(IN) :: KPARAL
!
!
!*       0.2   Declarations of local parameter
!              -------------------------------
!
INTEGER, DIMENSION(2), PARAMETER  :: IVAR_NODIMS  = (/1,1/) ! rank and number of bundles in coupling field
!
!
!*       0.3   Declarations of local variables
!              -------------------------------
!
INTEGER, DIMENSION(2)          :: IVAR_SHAPE  ! indexes for the coupling field local dimension
!
INTEGER                        :: IPART_ID ! Local partition ID
INTEGER                        :: IERR     ! Error info
!
INTEGER                        :: ILUOUT, IFLAG
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SFX_OASIS_DEFINE',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
#ifdef SFXOASIS
!-------------------------------------------------------------------------------
!
!
!*       0.     Initialize :
!               ------------
!
CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
CALL SFX_OASIS_CHECK(ILUOUT)
!
!-------------------------------------------------------------------------------
!
!*       1.     Define parallel partitions:
!               ---------------------------
!
CALL OASIS_DEF_PARTITION(IPART_ID,KPARAL(:),IERR)
!
IF(IERR/=OASIS_OK)THEN
   WRITE(ILUOUT,*)'SFX_OASIS_DEFINE: OASIS def partition problem, err = ',IERR
   CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def partition problem')
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       2.     Coupling fields shape :
!               -----------------------
!
IVAR_SHAPE(1)= 1
IVAR_SHAPE(2)= KNPTS
!
!-------------------------------------------------------------------------------
!
!*       3.     Sea variables for Surfex - Oasis coupling :
!               -------------------------------------------
!
IF(LCPL_SEA)THEN
!
! Sea output fields
!
  IF(LEN_TRIM(CSEA_FWSU)/=0)THEN
    CALL OASIS_DEF_VAR(NSEA_FWSU_ID,CSEA_FWSU,IPART_ID,IVAR_NODIMS,OASIS_OUT,IVAR_SHAPE,OASIS_DOUBLE,IERR)
    IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for sea zonal wind stress')
  ELSE
    NSEA_FWSU_ID=NUNDEF
  ENDIF
!
  IF(LEN_TRIM(CSEA_FWSV)/=0)THEN
    CALL OASIS_DEF_VAR(NSEA_FWSV_ID,CSEA_FWSV,IPART_ID,IVAR_NODIMS,OASIS_OUT,IVAR_SHAPE,OASIS_DOUBLE,IERR)
    IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for sea meridian wind stress')
  ELSE
    NSEA_FWSV_ID=NUNDEF
  ENDIF
!
  IF(LEN_TRIM(CSEA_HEAT)/=0)THEN
    CALL OASIS_DEF_VAR(NSEA_HEAT_ID,CSEA_HEAT,IPART_ID,IVAR_NODIMS,OASIS_OUT,IVAR_SHAPE,OASIS_DOUBLE,IERR)
    IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for sea Non solar net heat flux')
  ELSE
    NSEA_HEAT_ID=NUNDEF
  ENDIF
!
  IF(LEN_TRIM(CSEA_SNET)/=0)THEN
    CALL OASIS_DEF_VAR(NSEA_SNET_ID,CSEA_SNET,IPART_ID,IVAR_NODIMS,OASIS_OUT,IVAR_SHAPE,OASIS_DOUBLE,IERR)
    IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for sea Solar net heat')
  ELSE
    NSEA_SNET_ID=NUNDEF
  ENDIF
!
  IF(LEN_TRIM(CSEA_WIND)/=0)THEN
    CALL OASIS_DEF_VAR(NSEA_WIND_ID,CSEA_WIND,IPART_ID,IVAR_NODIMS,OASIS_OUT,IVAR_SHAPE,OASIS_DOUBLE,IERR)
    IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for sea 10m wind speed')
  ELSE
    NSEA_WIND_ID=NUNDEF
  ENDIF
!
  IF(LEN_TRIM(CSEA_FWSM)/=0)THEN
    CALL OASIS_DEF_VAR(NSEA_FWSM_ID,CSEA_FWSM,IPART_ID,IVAR_NODIMS,OASIS_OUT,IVAR_SHAPE,OASIS_DOUBLE,IERR)
    IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for sea wind stress')
  ELSE
    NSEA_FWSM_ID=NUNDEF
  ENDIF
!
  IF(LEN_TRIM(CSEA_EVAP)/=0)THEN
    CALL OASIS_DEF_VAR(NSEA_EVAP_ID,CSEA_EVAP,IPART_ID,IVAR_NODIMS,OASIS_OUT,IVAR_SHAPE,OASIS_DOUBLE,IERR)
    IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for sea Evaporation')
  ELSE
    NSEA_EVAP_ID=NUNDEF
  ENDIF
!
  IF(LEN_TRIM(CSEA_RAIN)/=0)THEN
    CALL OASIS_DEF_VAR(NSEA_RAIN_ID,CSEA_RAIN,IPART_ID,IVAR_NODIMS,OASIS_OUT,IVAR_SHAPE,OASIS_DOUBLE,IERR)
    IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for sea Rainfall rate')
  ELSE
    NSEA_RAIN_ID=NUNDEF
  ENDIF
!
  IF(LEN_TRIM(CSEA_SNOW)/=0)THEN
    CALL OASIS_DEF_VAR(NSEA_SNOW_ID,CSEA_SNOW,IPART_ID,IVAR_NODIMS,OASIS_OUT,IVAR_SHAPE,OASIS_DOUBLE,IERR)
    IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for sea Snowfall rate')
  ELSE
    NSEA_SNOW_ID=NUNDEF
  ENDIF
!
  IF(LEN_TRIM(CSEA_WATF)/=0)THEN
    CALL OASIS_DEF_VAR(NSEA_WATF_ID,CSEA_WATF,IPART_ID,IVAR_NODIMS,OASIS_OUT,IVAR_SHAPE,OASIS_DOUBLE,IERR)
    IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for sea freshwater rate')
  ELSE
    NSEA_WATF_ID=NUNDEF
  ENDIF
!
! Sea intput fields
!
  CALL OASIS_DEF_VAR(NSEA_SST_ID,CSEA_SST,IPART_ID,IVAR_NODIMS,OASIS_IN,IVAR_SHAPE,OASIS_DOUBLE,IERR)
  IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for Sea surface temperature')
!
  CALL OASIS_DEF_VAR(NSEA_UCU_ID,CSEA_UCU,IPART_ID,IVAR_NODIMS,OASIS_IN,IVAR_SHAPE,OASIS_DOUBLE,IERR)
  IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for Sea u-current stress')
!
  CALL OASIS_DEF_VAR(NSEA_VCU_ID,CSEA_VCU,IPART_ID,IVAR_NODIMS,OASIS_IN,IVAR_SHAPE,OASIS_DOUBLE,IERR)
  IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for Sea v-current stress')
!
! Particular case due to Sea-ice
!
  IF(LCPL_SEAICE)THEN
!
!   Output fields
!
    IF(LEN_TRIM(CSEAICE_HEAT)/=0)THEN
      CALL OASIS_DEF_VAR(NSEAICE_HEAT_ID,CSEAICE_HEAT,IPART_ID,IVAR_NODIMS,OASIS_OUT,IVAR_SHAPE,OASIS_DOUBLE,IERR)
      IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for Sea-ice non solar net heat')
    ELSE
      NSEAICE_HEAT_ID=NUNDEF
    ENDIF
!
    IF(LEN_TRIM(CSEAICE_SNET)/=0)THEN
      CALL OASIS_DEF_VAR(NSEAICE_SNET_ID,CSEAICE_SNET,IPART_ID,IVAR_NODIMS,OASIS_OUT,IVAR_SHAPE,OASIS_DOUBLE,IERR)
      IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for Sea-ice solar net heat flux')
    ELSE
      NSEAICE_SNET_ID=NUNDEF
    ENDIF
!
    IF(LEN_TRIM(CSEAICE_EVAP)/=0)THEN
      CALL OASIS_DEF_VAR(NSEAICE_EVAP_ID,CSEAICE_EVAP,IPART_ID,IVAR_NODIMS,OASIS_OUT,IVAR_SHAPE,OASIS_DOUBLE,IERR)
      IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for Sea-ice sublimation')
    ELSE
      NSEAICE_EVAP_ID=NUNDEF
    ENDIF
!
!   Intput fields
!
    CALL OASIS_DEF_VAR(NSEAICE_SIT_ID,CSEAICE_SIT,IPART_ID,IVAR_NODIMS,OASIS_IN,IVAR_SHAPE,OASIS_DOUBLE,IERR)
    IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for Sea-ice non solar net heat')
!
    CALL OASIS_DEF_VAR(NSEAICE_CVR_ID,CSEAICE_CVR,IPART_ID,IVAR_NODIMS,OASIS_IN,IVAR_SHAPE,OASIS_DOUBLE,IERR)
    IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for Sea-ice non solar net heat')
!
    CALL OASIS_DEF_VAR(NSEAICE_ALB_ID,CSEAICE_ALB,IPART_ID,IVAR_NODIMS,OASIS_IN,IVAR_SHAPE,OASIS_DOUBLE,IERR)
    IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for Sea-ice non solar net heat')
!
  ENDIF
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       4.     Lake variables for Surfex - Oasis coupling :
!               -------------------------------------------
!
IF(LCPL_LAKE)THEN
!
! Output fields
!
  CALL OASIS_DEF_VAR(NLAKE_EVAP_ID,CLAKE_EVAP,IPART_ID,IVAR_NODIMS,OASIS_OUT,IVAR_SHAPE,OASIS_DOUBLE,IERR)
  IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for lake Evaporation')
!
  CALL OASIS_DEF_VAR(NLAKE_RAIN_ID,CLAKE_RAIN,IPART_ID,IVAR_NODIMS,OASIS_OUT,IVAR_SHAPE,OASIS_DOUBLE,IERR)
  IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for lake Rainfall rate')
!
  CALL OASIS_DEF_VAR(NLAKE_SNOW_ID,CLAKE_SNOW,IPART_ID,IVAR_NODIMS,OASIS_OUT,IVAR_SHAPE,OASIS_DOUBLE,IERR)
  IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for lake Snowfall rate')
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       5.     Land surface variables for Surfex - Oasis coupling :
!               ----------------------------------------------------
!
IF(LCPL_LAND)THEN
!
! Output Surface runoff
!
  CALL OASIS_DEF_VAR(NRUNOFF_ID,CRUNOFF,IPART_ID,IVAR_NODIMS,OASIS_OUT,IVAR_SHAPE,OASIS_DOUBLE,IERR)
  IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for land Surface runoff')
!
! Output Calving flux
!
  IF(LCPL_CALVING)THEN
!
!     Output Calving flux
      CALL OASIS_DEF_VAR(NCALVING_ID,CCALVING,IPART_ID,IVAR_NODIMS,OASIS_OUT,IVAR_SHAPE,OASIS_DOUBLE,IERR)  
      IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for land Calving flux')
!
  ENDIF  
!
! Output Deep drainage
!
  CALL OASIS_DEF_VAR(NDRAIN_ID,CDRAIN,IPART_ID,IVAR_NODIMS,OASIS_OUT,IVAR_SHAPE,OASIS_DOUBLE,IERR)  
  IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for land Deep drainage')
!
! Particular case due to water table depth / surface coupling
!
  IF(LCPL_GW)THEN
!
!     Output groundwater recharge
      CALL OASIS_DEF_VAR(NRECHARGE_ID,CRECHARGE,IPART_ID,IVAR_NODIMS,OASIS_OUT,IVAR_SHAPE,OASIS_DOUBLE,IERR)  
      IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for land groundwater recharge')
!      
!     Input Water table depth
      CALL OASIS_DEF_VAR(NWTD_ID,CWTD,IPART_ID,IVAR_NODIMS,OASIS_IN,IVAR_SHAPE,OASIS_DOUBLE,IERR)  
      IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for land Water table depth')
!          
!     Input grid-cell fraction of WTD to rise
      CALL OASIS_DEF_VAR(NFWTD_ID,CFWTD,IPART_ID,IVAR_NODIMS,OASIS_IN,IVAR_SHAPE,OASIS_DOUBLE,IERR)  
      IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for land fraction of WTD to rise')
!   
  ENDIF      
!
! Particular case due to floodplains coupling
!
  IF(LCPL_FLOOD)THEN
!
!   Output Flood precip interception
    CALL OASIS_DEF_VAR(NPFLOOD_ID,CPFLOOD,IPART_ID,IVAR_NODIMS,OASIS_OUT,IVAR_SHAPE,OASIS_DOUBLE,IERR)  
    IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for land Flood precip interception')
!
!   Output Flood evaporation
    CALL OASIS_DEF_VAR(NEFLOOD_ID,CEFLOOD,IPART_ID,IVAR_NODIMS,OASIS_OUT,IVAR_SHAPE,OASIS_DOUBLE,IERR)  
    IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for land Flood evaporation')
!
!   Output Flood infiltration
    CALL OASIS_DEF_VAR(NIFLOOD_ID,CIFLOOD,IPART_ID,IVAR_NODIMS,OASIS_OUT,IVAR_SHAPE,OASIS_DOUBLE,IERR)  
    IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for land Flood infiltration')
!          
!   Input floodplains fraction
    CALL OASIS_DEF_VAR(NFFLOOD_ID,CFFLOOD,IPART_ID,IVAR_NODIMS,OASIS_IN,IVAR_SHAPE,OASIS_DOUBLE,IERR)  
    IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for land Floodplains fraction')
!          
!   Input floodplains potential infiltration
    CALL OASIS_DEF_VAR(NPIFLOOD_ID,CPIFLOOD,IPART_ID,IVAR_NODIMS,OASIS_IN,IVAR_SHAPE,OASIS_DOUBLE,IERR)  
    IF(IERR/=OASIS_OK) CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS def var problem for land Floodplains potential infiltration')
!    
  ENDIF
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       6.     End of declaration phase:
!               --------------
!
CALL OASIS_ENDDEF(IERR)
!
IF(IERR/=OASIS_OK)THEN
   WRITE(ILUOUT,*)'SFX_OASIS_DEFINE: OASIS enddef problem, err = ',IERR
   CALL ABOR1_SFX('SFX_OASIS_DEFINE: OASIS enddef problem')
ENDIF
!
!-------------------------------------------------------------------------------
#endif
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SFX_OASIS_DEFINE',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SFX_OASIS_DEFINE
