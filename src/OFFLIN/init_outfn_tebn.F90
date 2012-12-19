!     #########
       SUBROUTINE INIT_OUTFN_TEB_n(HPROGRAM,KLUOUT)
!     ###############################
!
!
!!****  *INIT_OUTFN_TEB_n* -  create output files and defines variables
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
!!      P. LeMoigne   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04-04  P. LeMoigne
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_OL_FILEID,       ONLY : XVAR_TO_FILEOUT, XID, XOUT
USE MODD_TEB_n,           ONLY: NROOF_LAYER, NROAD_LAYER, NWALL_LAYER, &
                                XZ0_TOWN, TTIME
USE MODD_TEB_CANOPY_n,    ONLY: NLVL
USE MODD_DIAG_TEB_n
USE MODD_DIAG_MISC_TEB_n ,ONLY : LSURF_MISC_BUDGET, XQF_BLD,      &
                                 XFLX_BLD, XQF_TOWN, XDQS_TOWN,   &
                                 XH_WALL_A, XH_WALL_B,XH_ROOF,XH_ROAD
!
USE MODN_IO_OFFLINE,      ONLY : XTSTEP_OUTPUT
!
USE MODI_GET_DIM_FULL_n
USE MODI_OL_DEFINE_DIM
USE MODI_GET_DATE_OL
USE MODI_CREATE_FILE
USE MODI_DEF_VAR_NETCDF
USE MODI_OL_WRITE_COORD
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
include 'netcdf.inc'

CHARACTER(LEN=6),  INTENT(IN) :: HPROGRAM
INTEGER,           INTENT(IN) :: KLUOUT
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER                          :: INI
INTEGER                          :: IDIM1
CHARACTER(LEN=13),DIMENSION(1)   :: YUNIT1, YUNIT2
REAL,DIMENSION(:), POINTER       :: ZX, ZY
INTEGER, DIMENSION(:), POINTER   :: IDIMS, IDDIM
CHARACTER(LEN=100), DIMENSION(:), POINTER :: YNAME_DIM
!
CHARACTER(LEN=40),DIMENSION(1)   :: YDATE
INTEGER                          :: IFILE_ID
CHARACTER(LEN=50)                :: YFILE
CHARACTER(LEN=100), DIMENSION(1) :: YATT_TITLE, YATT
CHARACTER(LEN=3)                 :: YPAS, YLVL
INTEGER                          :: JLAYER
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------

! 1. Compute output lenght dimension
!-----------------------------------
IF (LHOOK) CALL DR_HOOK('INIT_OUTFN_TEB_N',0,ZHOOK_HANDLE)
CALL GET_DIM_FULL_n(INI)
CALL OL_DEFINE_DIM(HPROGRAM, KLUOUT, INI, IDIM1, YUNIT1, YUNIT2, &
                   ZX, ZY, IDIMS, IDDIM, YNAME_DIM)
CALL GET_DATE_OL(TTIME,XTSTEP_OUTPUT,YDATE(1))
!
! 4. Create output file for fluxes values
!----------------------------------------------------------

IF (ALLOCATED(XVAR_TO_FILEOUT)) DEALLOCATE(XVAR_TO_FILEOUT)
IF (ALLOCATED(XID)) DEALLOCATE(XID)
ALLOCATE(XVAR_TO_FILEOUT(0))
ALLOCATE(XID(0))
XOUT=0
!
YATT_TITLE(1)='units'
!
YFILE='TEB_DIAGNOSTICS.OUT.nc'
CALL CREATE_FILE(YFILE,IDIMS,YNAME_DIM,IFILE_ID,IDDIM)
YATT='-'
!
IF (LSURF_MISC_BUDGET) THEN
  CALL DEF_VAR_NETCDF(IFILE_ID,'Z0_TOWN'  ,'Town_Rougness_Length',IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'XQF_BLD'  ,'Domestic_heating',IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'XQF_BLDWFR'  ,'Domestic_heating',IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'XFLX_BLD'  ,'Heat flux from bld',IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'XQF_TOWN' ,'Anthropogenic heat',IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'XDQS_TOWN'  ,'Storage',IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'XH_WALL_A'  ,'Wall sensible heat flux',IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'XH_WALL_B'  ,'Wall sensible heat flux',IDDIM,YATT_TITLE,YATT)  
  CALL DEF_VAR_NETCDF(IFILE_ID,'XH_ROOF'  ,'Roof sensible heat flux',IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'XH_ROAD'  ,'Road sensible heat flux',IDDIM,YATT_TITLE,YATT)
ENDIF
!
IF (LCOEF)  THEN
  YATT='W/s2'
  CALL DEF_VAR_NETCDF(IFILE_ID,'CD_TEB'   ,'Drag_Coefficient_For_Momentum   '   ,IDDIM,YATT_TITLE,YATT)
  YATT='W/s'
  CALL DEF_VAR_NETCDF(IFILE_ID,'CH_TEB'   ,'Drag_Coefficient_For_Heat       '   ,IDDIM,YATT_TITLE,YATT)
  YATT='W/s/K'
  CALL DEF_VAR_NETCDF(IFILE_ID,'CE_TEB'   ,'Drag_Coefficient_For_Evaporation'   ,IDDIM,YATT_TITLE,YATT)
  YATT='m'
  CALL DEF_VAR_NETCDF(IFILE_ID,'Z0_TEB'   ,'Roughness_Length_For_Momentum'   ,IDDIM,YATT_TITLE,YATT)
  YATT='m'
  CALL DEF_VAR_NETCDF(IFILE_ID,'Z0H_TEB'  ,'Roughness_Length_For_Heat'   ,IDDIM,YATT_TITLE,YATT)
ENDIF
!
IF (LSURF_VARS) THEN
  YATT='kg/kg'
  CALL DEF_VAR_NETCDF(IFILE_ID,'QS_TEB'   ,'Surface_Humidity   '   ,IDDIM,YATT_TITLE,YATT)
ENDIF
!
IF (N2M>0) THEN
  YATT='K'
  CALL DEF_VAR_NETCDF(IFILE_ID,'T2M_TEB' ,'2m_Temperature         '   ,IDDIM,YATT_TITLE,YATT)
  YATT='kg/kg'
  CALL DEF_VAR_NETCDF(IFILE_ID,'Q2M_TEB' ,'2m_Specific_Humidity   '   ,IDDIM,YATT_TITLE,YATT)
  YATT='(-)'
  CALL DEF_VAR_NETCDF(IFILE_ID,'HU2M_TEB','2m_Relative_Humidity   '   ,IDDIM,YATT_TITLE,YATT)
  YATT='m/s'
  CALL DEF_VAR_NETCDF(IFILE_ID,'ZON10M_TEB','10m_Zonal_wind       '   ,IDDIM,YATT_TITLE,YATT)
  YATT='m/s'
  CALL DEF_VAR_NETCDF(IFILE_ID,'MER10M_TEB','2m_Meridian_Wind     '   ,IDDIM,YATT_TITLE,YATT)
ENDIF

IF (N2M>0) THEN
  YATT='-'
  CALL DEF_VAR_NETCDF(IFILE_ID,'RI_TEB'   ,'Averaged_Richardson_Number'                ,IDDIM,YATT_TITLE,YATT)
ENDIF
!
YATT='W/m2'
CALL DEF_VAR_NETCDF(IFILE_ID,'RN_TEB'   ,'Averaged_Net_Radiation'                    ,IDDIM,YATT_TITLE,YATT)
CALL DEF_VAR_NETCDF(IFILE_ID,'H_TEB'    ,'Averaged_Sensible_Heat_Flux'               ,IDDIM,YATT_TITLE,YATT)
CALL DEF_VAR_NETCDF(IFILE_ID,'LE_TEB'   ,'Averaged_Latent_Heat_Flux  '               ,IDDIM,YATT_TITLE,YATT)
CALL DEF_VAR_NETCDF(IFILE_ID,'GFLUX_TEB','Averaged_Ground_Heat_Flux  '               ,IDDIM,YATT_TITLE,YATT)
CALL DEF_VAR_NETCDF(IFILE_ID,'SWD_TEB'  ,'Averaged_Downward_SW       '               ,IDDIM,YATT_TITLE,YATT)
CALL DEF_VAR_NETCDF(IFILE_ID,'SWU_TEB'  ,'Averaged_Upward_SW         '               ,IDDIM,YATT_TITLE,YATT)
CALL DEF_VAR_NETCDF(IFILE_ID,'LWD_TEB'  ,'Averaged_Downward_LW       '               ,IDDIM,YATT_TITLE,YATT)
CALL DEF_VAR_NETCDF(IFILE_ID,'LWU_TEB'  ,'Averaged_Upward_LW         '               ,IDDIM,YATT_TITLE,YATT)
CALL DEF_VAR_NETCDF(IFILE_ID,'FMU_TEB'  ,'Averaged_Zonal_Wind_Stress '               ,IDDIM,YATT_TITLE,YATT)
CALL DEF_VAR_NETCDF(IFILE_ID,'FMV_TEB'  ,'Averaged_Merid_Wind_Stress '               ,IDDIM,YATT_TITLE,YATT)
!
CALL OL_WRITE_COORD(YFILE,IFILE_ID,IDDIM,YATT_TITLE,YNAME_DIM,YUNIT1,YUNIT2,IDIM1,YDATE,ZX,ZY)
!
!
!YFILE='TEB_DIAG_CUMUL.OUT.nc'
!CALL CREATE_FILE(YFILE,IDIMS,YNAME_DIM,IFILE_ID,IDDIM)
!YATT='dimensionless'
!
!CALL DEF_VAR_NETCDF(IFILE_ID,'RNC_TEB'  ,'Averaged_Cumulated_Net_Radiation'        ,IDDIM,YATT_TITLE,YATT)
!CALL DEF_VAR_NETCDF(IFILE_ID,'HC_TEB'   ,'Averaged_Cumulated_Sensible_Heat_Flux'   ,IDDIM,YATT_TITLE,YATT)
!CALL DEF_VAR_NETCDF(IFILE_ID,'LEC_TEB'  ,'Averaged_Cumulated_Latent_Heat_Flux'     ,IDDIM,YATT_TITLE,YATT)
!CALL DEF_VAR_NETCDF(IFILE_ID,'GFLUXC_TEB','Averaged_Cumulated_Ground_Heat_Flux'    ,IDDIM,YATT_TITLE,YATT)
!
!CALL OL_WRITE_COORD(YFILE,IFILE_ID,IDDIM,YATT_TITLE,YNAME_DIM,YUNIT1,YUNIT2,IDIM1,YDATE,ZX,ZY)
!
!
! 5. Create output file for prognostic and semi-pronostic variable in teb
!------------------------------------------------------------------------

YFILE='TEB_PROGNOSTIC.OUT.nc'
CALL CREATE_FILE(YFILE,IDIMS,YNAME_DIM,IFILE_ID,IDDIM)
!
! 5.1 Temperatures
YATT='K'
!
! Roof temperatures
DO JLAYER=1,NROOF_LAYER
  WRITE(YPAS,'(I3)') JLAYER ; YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
  CALL DEF_VAR_NETCDF(IFILE_ID,'T_ROOF'//YLVL,'Roof_Temperature_Layer_'//YLVL ,IDDIM,YATT_TITLE,YATT)
END DO
! Road temperatures
DO JLAYER=1,NROAD_LAYER
  WRITE(YPAS,'(I3)') JLAYER ; YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
  CALL DEF_VAR_NETCDF(IFILE_ID,'T_ROAD'//YLVL,'Road_Temperature_Layer_'//YLVL ,IDDIM,YATT_TITLE,YATT)
END DO
! Wall temperatures
DO JLAYER=1,NWALL_LAYER
  WRITE(YPAS,'(I3)') JLAYER ; YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
  CALL DEF_VAR_NETCDF(IFILE_ID,'T_WALL_A'//YLVL,'Wall_Temperature_Layer_'//YLVL ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'T_WALL_B'//YLVL,'Wall_Temperature_Layer_'//YLVL ,IDDIM,YATT_TITLE,YATT)  
END DO
! Internal building temperature
CALL DEF_VAR_NETCDF(IFILE_ID,'TI_BLD','Internal_Building_Temperature',IDDIM,YATT_TITLE,YATT)
!
! Deep road temperature
CALL DEF_VAR_NETCDF(IFILE_ID,'TI_ROAD','Deep_Road_Temperature',IDDIM,YATT_TITLE,YATT)
!
! 5.2 Water contents
!
YATT = 'kg/m2'
! Roof water content
CALL DEF_VAR_NETCDF(IFILE_ID,'WS_ROOF','Roof_Water_Content_Layer' ,IDDIM,YATT_TITLE,YATT)
!
! Road water content
CALL DEF_VAR_NETCDF(IFILE_ID,'WS_ROAD','Road_Water_Content_Layer' ,IDDIM,YATT_TITLE,YATT)
!
! 5.3 semi pronostic variables
!
! temperature of canyon air
YATT = 'K'
CALL DEF_VAR_NETCDF(IFILE_ID,'T_CANYON','Canyon_Air_Temperature',IDDIM,YATT_TITLE,YATT)
!
! humidity of canyon air
YATT = 'kg/kg'
CALL DEF_VAR_NETCDF(IFILE_ID,'Q_CANYON','Canyon_Air_Humidity',IDDIM,YATT_TITLE,YATT)
!
CALL OL_WRITE_COORD(YFILE,IFILE_ID,IDDIM,YATT_TITLE,YNAME_DIM,YUNIT1,YUNIT2,IDIM1,YDATE,ZX,ZY)
!
! 6. Create output file for prognostic variables in teb canopy
!------------------------------------------------------------------------
!
YFILE='TEB_CANOPY.OUT.nc'
CALL CREATE_FILE(YFILE,IDIMS,YNAME_DIM,IFILE_ID,IDDIM)
!
! 6.1 Heights
YATT = 'm'
!
DO JLAYER=1,NLVL
  WRITE(YPAS,'(I2.2,1X)') JLAYER ; YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
  CALL DEF_VAR_NETCDF(IFILE_ID,'TEB_CAN_Z'//YLVL,'Height_of_canopy_Layer_'//YLVL ,IDDIM,YATT_TITLE,YATT)
END DO
!
! 6.2 Wind
YATT = 'm/s'
DO JLAYER=1,NLVL
  WRITE(YPAS,'(I2.2,1X)') JLAYER ; YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
  CALL DEF_VAR_NETCDF(IFILE_ID,'TEB_CAN_U'//YLVL,'Wind_of_canopy_Layer_'//YLVL ,IDDIM,YATT_TITLE,YATT)
END DO
!
! 6.3 Temperature
YATT = 'K'
DO JLAYER=1,NLVL
  WRITE(YPAS,'(I2.2,1X)') JLAYER ; YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
  CALL DEF_VAR_NETCDF(IFILE_ID,'TEB_CAN_T'//YLVL,'Temperature_of_canopy_Layer_'//YLVL ,IDDIM,YATT_TITLE,YATT)
END DO
!
! 6.4 Temperature
YATT = 'kg/m3'
DO JLAYER=1,NLVL
  WRITE(YPAS,'(I2.2,1X)') JLAYER ; YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
  CALL DEF_VAR_NETCDF(IFILE_ID,'TEB_CAN_Q'//YLVL,'Humidity_of_canopy_Layer_'//YLVL ,IDDIM,YATT_TITLE,YATT)
END DO
!
! 6.5 Turbulence
YATT = 'm2/s2'
DO JLAYER=1,NLVL
  WRITE(YPAS,'(I2.2,1X)') JLAYER ; YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
  CALL DEF_VAR_NETCDF(IFILE_ID,'TEB_CAN_E'//YLVL,'TKE_of_canopy_Layer_'//YLVL ,IDDIM,YATT_TITLE,YATT)
END DO
!
CALL OL_WRITE_COORD(YFILE,IFILE_ID,IDDIM,YATT_TITLE,YNAME_DIM,YUNIT1,YUNIT2,IDIM1,YDATE,ZX,ZY)
!
IF (LHOOK) CALL DR_HOOK('INIT_OUTFN_TEB_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE INIT_OUTFN_TEB_n
