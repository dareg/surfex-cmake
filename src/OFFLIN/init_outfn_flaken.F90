!     #########
       SUBROUTINE INIT_OUTFN_FLAKE_n(HPROGRAM,KLUOUT)
!     ###############################
!
!!****  *INIT_OUTFN_FLAKE_n* -  create output files and defines variables
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
!!      Modified 06-10 by S. Faroux
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------

USE MODD_OL_FILEID,    ONLY : XVAR_TO_FILEOUT, XID, XOUT
USE MODD_DIAG_FLAKE_n
USE MODD_FLAKE_n,      ONLY : LSBL, TTIME
USE MODD_FLAKE_SBL_n,  ONLY : NLVL
USE MODN_IO_OFFLINE,   ONLY : XTSTEP_OUTPUT
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
INTEGER, INTENT(IN)           :: KLUOUT
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER                          :: INI, ICANLVL
INTEGER                          :: IDIM1
 CHARACTER(LEN=13),DIMENSION(1)   :: YUNIT1, YUNIT2
REAL,DIMENSION(:), POINTER       :: ZX, ZY
INTEGER, DIMENSION(:), POINTER   :: IDIMS, IDDIM
 CHARACTER(LEN=100), DIMENSION(:), POINTER :: YNAME_DIM
!
 CHARACTER(LEN=40),DIMENSION(1)   :: YDATE
INTEGER                          :: IFILE_ID, IVAR_ID
 CHARACTER(LEN=50)                :: YFILE
 CHARACTER(LEN=100), DIMENSION(1) :: YATT_TITLE, YATT
 CHARACTER(LEN=2)                 :: YLVLV
INTEGER                          :: JLAYER
INTEGER                          :: JRET
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------

! 1. Compute output lenght dimension
!-----------------------------------
IF (LHOOK) CALL DR_HOOK('INIT_OUTFN_FLAKE_N',0,ZHOOK_HANDLE)
 CALL GET_DIM_FULL_n(INI)
!
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
YFILE='WATFLUX_DIAGNOSTICS.OUT.nc'
 CALL CREATE_FILE(YFILE,IDIMS,YNAME_DIM,IFILE_ID,IDDIM)
YATT='dimensionless'
!
IF (LCOEF) THEN
  YATT='W/s2'
  CALL DEF_VAR_NETCDF(IFILE_ID,'CD_WAT'   ,'Drag_Coefficient_For_Momentum   '   ,IDDIM,YATT_TITLE,YATT)
  YATT='W/s'
  CALL DEF_VAR_NETCDF(IFILE_ID,'CH_WAT'   ,'Drag_Coefficient_For_Heat       '   ,IDDIM,YATT_TITLE,YATT)
  YATT='W/s/K'
  CALL DEF_VAR_NETCDF(IFILE_ID,'CE_WAT'   ,'Drag_Coefficient_For_Evaporation'   ,IDDIM,YATT_TITLE,YATT)
  YATT='m'
  CALL DEF_VAR_NETCDF(IFILE_ID,'Z0_WAT'   ,'Roughness_Length_For_Momentum'   ,IDDIM,YATT_TITLE,YATT)
  YATT='m'
  CALL DEF_VAR_NETCDF(IFILE_ID,'Z0H_WAT'  ,'Roughness_Length_For_Heat'   ,IDDIM,YATT_TITLE,YATT)
ENDIF
!
IF (LSURF_VARS) THEN
  YATT='kg/kg'
  CALL DEF_VAR_NETCDF(IFILE_ID,'QS_WAT'   ,'Surface_Humidity   '   ,IDDIM,YATT_TITLE,YATT)
ENDIF
!
IF (N2M>0) THEN
  YATT='K'
  CALL DEF_VAR_NETCDF(IFILE_ID,'T2M_WAT' ,'2m_Temperature         '   ,IDDIM,YATT_TITLE,YATT)
  YATT='kg/kg'
  CALL DEF_VAR_NETCDF(IFILE_ID,'Q2M_WAT' ,'2m_Specific_Humidity   '   ,IDDIM,YATT_TITLE,YATT)
  YATT='(-)'
  CALL DEF_VAR_NETCDF(IFILE_ID,'HU2M_WAT','2m_Relative_Humidity   '   ,IDDIM,YATT_TITLE,YATT)
  YATT='m/s'
  CALL DEF_VAR_NETCDF(IFILE_ID,'ZON10M_WAT','10m_Zonal_wind       '   ,IDDIM,YATT_TITLE,YATT)
  YATT='m/s'
  CALL DEF_VAR_NETCDF(IFILE_ID,'MER10M_WAT','2m_Meridian_Wind     '   ,IDDIM,YATT_TITLE,YATT)
ENDIF
!
IF (N2M>0) THEN
  YATT='-'
   CALL DEF_VAR_NETCDF(IFILE_ID,'RI_WAT'      ,'Averaged_Richardson_Number'                ,IDDIM,YATT_TITLE,YATT)
ENDIF
!
IF(LSURF_BUDGET) THEN
  YATT='W/m2'
  CALL DEF_VAR_NETCDF(IFILE_ID,'RN_WAT'      ,'Averaged_Net_Radiation'                    ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'H_WAT'       ,'Averaged_Sensible_Heat_Flux'               ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LE_WAT'      ,'Averaged_Latent_Heat_Flux  '               ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'GFLUX_WAT'   ,'Averaged_Ground_Heat_Flux  '               ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'SWD_WAT'     ,'Averaged_Downward_SW       '               ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'SWU_WAT'     ,'Averaged_Upward_SW         '               ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LWD_WAT'     ,'Averaged_Downward_LW       '               ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LWU_WAT'     ,'Averaged_Upward_LW         '               ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'FMU_WAT'     ,'Averaged_Zonal_Wind_Stress '               ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'FMV_WAT'     ,'Averaged_Merid_Wind_Stress '               ,IDDIM,YATT_TITLE,YATT)
ENDIF
!
 CALL OL_WRITE_COORD(YFILE,IFILE_ID,IDDIM,YATT_TITLE,YNAME_DIM,YUNIT1,YUNIT2,IDIM1,YDATE,ZX,ZY)
!
!YFILE='WATFLUX_DIAG_CUMUL.OUT.nc'
!CALL CREATE_FILE(YFILE,IDIMS,YNAME_DIM,IFILE_ID,IDDIM)
!YATT      (1)='dimensionless'
!CALL DEF_VAR_NETCDF(IFILE_ID,'RNC_WAT'     ,'Averaged_Cumulated_Net_Radiation'        ,IDDIM,YATT_TITLE,YATT)
!CALL DEF_VAR_NETCDF(IFILE_ID,'HC_WAT'      ,'Averaged_Cumulated_Sensible_Heat_Flux'   ,IDDIM,YATT_TITLE,YATT)
!CALL DEF_VAR_NETCDF(IFILE_ID,'LEC_WAT'     ,'Averaged_Cumulated_Latent_Heat_Flux'     ,IDDIM,YATT_TITLE,YATT)
!CALL DEF_VAR_NETCDF(IFILE_ID,'GFLUXC_WAT'   ,'Averaged_Cumulated_Ground_Heat_Flux'    ,IDDIM,YATT_TITLE,YATT)
!CALL OL_WRITE_COORD(YFILE,IFILE_ID,IDDIM,YATT_TITLE,YNAME_DIM,YUNIT1,YUNIT2,IDIM1,YDATE,ZX,ZY)
!
!
YFILE='WATFLUX_PROGNOSTIC.OUT.nc'
 CALL CREATE_FILE(YFILE,IDIMS,YNAME_DIM,IFILE_ID,IDDIM)
YATT='dimensionless'
!
YATT='K'
 CALL DEF_VAR_NETCDF(IFILE_ID,'TS_WATER'   , 'Averaged_Water_S_Temperature'                 ,IDDIM,YATT_TITLE,YATT)
 CALL DEF_VAR_NETCDF(IFILE_ID,'T_MNW'      , 'Averaged_Water_Temperature  '                 ,IDDIM,YATT_TITLE,YATT)
 CALL DEF_VAR_NETCDF(IFILE_ID,'T_WML'      , 'Mixed_layer_wat_temperature '                 ,IDDIM,YATT_TITLE,YATT)
 CALL DEF_VAR_NETCDF(IFILE_ID,'T_BOT'      , 'Bottom_Water_Temperature    '                 ,IDDIM,YATT_TITLE,YATT)
YATT='m' 
 CALL DEF_VAR_NETCDF(IFILE_ID,'H_ML'       , 'Mixed_Layer_Depth           '                 ,IDDIM,YATT_TITLE,YATT)
YATT=' '
 CALL DEF_VAR_NETCDF(IFILE_ID,'CT'         , 'Termocline_Shape_Factor        '              ,IDDIM,YATT_TITLE,YATT)
YATT='K'
 CALL DEF_VAR_NETCDF(IFILE_ID,'T_SNOW'     , 'Temperature at the air-snow interface'        ,IDDIM,YATT_TITLE,YATT)
 CALL DEF_VAR_NETCDF(IFILE_ID,'T_ICE'      , 'Ice_surface_Temperature  '                    ,IDDIM,YATT_TITLE,YATT)
 CALL DEF_VAR_NETCDF(IFILE_ID,'T_B1'       , 'Temperature of the upper layer of sediments ' ,IDDIM,YATT_TITLE,YATT)
YATT='m'
 CALL DEF_VAR_NETCDF(IFILE_ID,'H_SNOW'     , 'Snow thickness'              ,IDDIM,YATT_TITLE,YATT)
 CALL DEF_VAR_NETCDF(IFILE_ID,'H_ICE'      , 'Ice thickness'               ,IDDIM,YATT_TITLE,YATT)
 CALL DEF_VAR_NETCDF(IFILE_ID,'H_B1'       , 'Thickness of the upper layer of sediments'     ,IDDIM,YATT_TITLE,YATT)
!
IF (LSBL) THEN
  DO JLAYER=1,ICANLVL
    WRITE(YLVLV,'(i2.2)') JLAYER
    CALL DEF_VAR_NETCDF(IFILE_ID,'WAT_SBL_Z'//YLVLV,'Canopy height',   IDDIM,YATT_TITLE,(/'m'/))
    CALL DEF_VAR_NETCDF(IFILE_ID,'WAT_SBL_U'//YLVLV,'Canopy wind',     IDDIM,YATT_TITLE,(/'m/s'/))
    CALL DEF_VAR_NETCDF(IFILE_ID,'WAT_SBL_T'//YLVLV,'Canopy temp',     IDDIM,YATT_TITLE,(/'K'/))
    CALL DEF_VAR_NETCDF(IFILE_ID,'WAT_SBL_Q'//YLVLV,'Canopy humidity', IDDIM,YATT_TITLE,(/'kg/m3'/))
    CALL DEF_VAR_NETCDF(IFILE_ID,'WAT_SBL_E'//YLVLV,'Canopy TKE',      IDDIM,YATT_TITLE,(/'m2/s2'/))
  END DO
ENDIF
!
 CALL OL_WRITE_COORD(YFILE,IFILE_ID,IDDIM,YATT_TITLE,YNAME_DIM,YUNIT1,YUNIT2,IDIM1,YDATE,ZX,ZY)
!
IF (LHOOK) CALL DR_HOOK('INIT_OUTFN_FLAKE_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE INIT_OUTFN_FLAKE_n
