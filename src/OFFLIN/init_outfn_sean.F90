!     #########
       SUBROUTINE INIT_OUTFN_SEA_n(HPROGRAM,KLUOUT)
!     ###############################
!
!!****  *INIT_OUTFN_SEA_n* -  create output files and defines variables
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
USE MODD_DIAG_SURF_ATM_n,ONLY : LPROVAR_TO_DIAG
!
USE MODD_OL_FILEID,      ONLY : XVAR_TO_FILEOUT, XID, XOUT
USE MODD_SEAFLUX_n,      ONLY : LSBL, TTIME
USE MODD_SEAFLUX_SBL_n,  ONLY : NLVL
USE MODD_DIAG_SEAFLUX_n
USE MODD_OCEAN_n
USE MODD_OCEAN_GRID_n, ONLY: NOCKMIN, NOCKMAX
USE MODD_DIAG_OCEAN_n, ONLY: LDIAG_OCEAN, XTOCMOY, XSOCMOY, XUOCMOY, XVOCMOY, XDOCMOY
!
USE MODN_IO_OFFLINE,    ONLY: XTSTEP_OUTPUT
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
INTEGER                          :: INI
INTEGER                          :: IDIM1, INDIMS
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
IF (LHOOK) CALL DR_HOOK('INIT_OUTFN_SEA_N',0,ZHOOK_HANDLE)
 CALL GET_DIM_FULL_n(INI)
 CALL OL_DEFINE_DIM(HPROGRAM, KLUOUT, INI, IDIM1, YUNIT1, YUNIT2, &
                   ZX, ZY, IDIMS, IDDIM, YNAME_DIM)
 CALL GET_DATE_OL(TTIME,XTSTEP_OUTPUT,YDATE(1))
!
INDIMS = SIZE(IDDIM)
!
! 4. Create output file for fluxes values
!----------------------------------------------------------
!
IF (ALLOCATED(XVAR_TO_FILEOUT)) DEALLOCATE(XVAR_TO_FILEOUT)
IF (ALLOCATED(XID)) DEALLOCATE(XID)
ALLOCATE(XVAR_TO_FILEOUT(0))
ALLOCATE(XID(0))
XOUT=0
!
YATT_TITLE(1)='units'
!
YFILE='SEAFLUX_DIAGNOSTICS.OUT.nc'
 CALL CREATE_FILE(YFILE,IDIMS,YNAME_DIM,IFILE_ID,IDDIM)
YATT='dimensionless'
!
IF(LPROVAR_TO_DIAG) THEN
  YATT='K'
  CALL DEF_VAR_NETCDF(IFILE_ID,'SST'   ,'Sea_Surface_Temperature'               ,IDDIM,YATT_TITLE,YATT)
ENDIF
!
IF (LCOEF) THEN
  YATT='W/s2'
  CALL DEF_VAR_NETCDF(IFILE_ID,'CD_SEA'   ,'Drag_Coefficient_For_Momentum   ',IDDIM,YATT_TITLE,YATT)
  YATT='W/s'
  CALL DEF_VAR_NETCDF(IFILE_ID,'CH_SEA'   ,'Drag_Coefficient_For_Heat       ',IDDIM,YATT_TITLE,YATT)
  YATT='W/s/K'
  CALL DEF_VAR_NETCDF(IFILE_ID,'CE_SEA'   ,'Drag_Coefficient_For_Evaporation',IDDIM,YATT_TITLE,YATT)
  YATT='m'
  CALL DEF_VAR_NETCDF(IFILE_ID,'Z0_SEA'   ,'Roughness_Length_For_Momentum   ',IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'Z0H_SEA'  ,'Roughness_Length_For_Heat       ',IDDIM,YATT_TITLE,YATT)
ENDIF
!
IF (LSURF_VARS) THEN
   YATT='kg/kg'
   CALL DEF_VAR_NETCDF(IFILE_ID,'QS_SEA'   ,'Surface_Humidity   '   ,IDDIM,YATT_TITLE,YATT)
ENDIF
!
IF (N2M>0) THEN
   YATT='K'
   CALL DEF_VAR_NETCDF(IFILE_ID,'T2M_SEA' ,'2m_Temperature         '   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'T2MMIN_SEA' ,'Minimum_2m_Temperature   '   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'T2MMAX_SEA' ,'Maximum_2m_Temperature   '   ,IDDIM,YATT_TITLE,YATT)
   YATT='kg/kg'
   CALL DEF_VAR_NETCDF(IFILE_ID,'Q2M_SEA' ,'2m_Specific_Humidity   '   ,IDDIM,YATT_TITLE,YATT)
   YATT='(-)'
   CALL DEF_VAR_NETCDF(IFILE_ID,'HU2M_SEA','2m_Relative_Humidity   '   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'HU2MMIN_SEA','Minimum_2m_Relative_Humidity   ' ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'HU2MMAX_SEA','Maximum_2m_Relative_Humidity   ' ,IDDIM,YATT_TITLE,YATT)
   YATT='m/s'
   CALL DEF_VAR_NETCDF(IFILE_ID,'ZON10M_SEA','10m_Zonal_wind       '   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'MER10M_SEA','10m_Meridian_Wind     '   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'W10M_SEA','10m_Wind     '   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'W10MMAX_SEA','Maximum_10m_Wind  '   ,IDDIM,YATT_TITLE,YATT)
ENDIF
!
IF (N2M>0) THEN
   YATT='-'
   CALL DEF_VAR_NETCDF(IFILE_ID,'RI_SEA'   ,'Averaged_Richardson_Number'                ,IDDIM,YATT_TITLE,YATT)
ENDIF
!
IF (LSURF_BUDGET) THEN
   YATT='W/m2'
   CALL DEF_VAR_NETCDF(IFILE_ID,'RN_SEA'   ,'Averaged_Net_Radiation'                    ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'H_SEA'    ,'Averaged_Sensible_Heat_Flux'               ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'LE_SEA'   ,'Averaged_Total_Latent_Heat_Flux  '         ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'LEI_SEA'  ,'Averaged_SublimationLatent_Heat_Flux  '    ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'GFLUX_SEA','Averaged_Ground_Heat_Flux  '               ,IDDIM,YATT_TITLE,YATT)
   IF(LRAD_BUDGET)THEN
      CALL DEF_VAR_NETCDF(IFILE_ID,'SWD_SEA'  ,'Averaged_Downward_SW       '             ,IDDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(IFILE_ID,'SWU_SEA'  ,'Averaged_Upward_SW         '             ,IDDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(IFILE_ID,'LWD_SEA'  ,'Averaged_Downward_LW       '             ,IDDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(IFILE_ID,'LWU_SEA'  ,'Averaged_Upward_LW         '             ,IDDIM,YATT_TITLE,YATT)
   ENDIF
   YATT='kg/ms2'
   CALL DEF_VAR_NETCDF(IFILE_ID,'FMU_SEA'  ,'Averaged_Zonal_Wind_Stress '               ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'FMV_SEA'  ,'Averaged_Merid_Wind_Stress '               ,IDDIM,YATT_TITLE,YATT)
ENDIF
!
IF (LDIAG_OCEAN) THEN
     ! Mean cmo temperature
     YATT='K'
     CALL DEF_VAR_NETCDF(IFILE_ID,'TOML','Mean cmo Tempe ',IDDIM,YATT_TITLE,YATT)
     ! Mean cmo salinity
     YATT='psu'
     CALL DEF_VAR_NETCDF(IFILE_ID,'SOML','Mean cmo Salin ',IDDIM,YATT_TITLE,YATT)
     ! Mean cmo U-current
     YATT='m/s'
     CALL DEF_VAR_NETCDF(IFILE_ID,'UOML','Mean cmo U-cur ',IDDIM,YATT_TITLE,YATT)
     ! Mean cmo V-current
     YATT='m/s'
     CALL DEF_VAR_NETCDF(IFILE_ID,'VOML','Mean cmo V-cur ',IDDIM,YATT_TITLE,YATT)
     ! Mean cmo density
     YATT='Kg/m3'
     CALL DEF_VAR_NETCDF(IFILE_ID,'DOML','Mean cmo Densi ',IDDIM,YATT_TITLE,YATT)

ENDIF
!
 CALL OL_WRITE_COORD(YFILE,IFILE_ID,IDDIM,YATT_TITLE,YNAME_DIM,YUNIT1,YUNIT2,IDIM1,YDATE,ZX,ZY)
!
!
IF (LSURF_BUDGETC) THEN
  !        
  YFILE='SEAFLUX_DIAG_CUMUL.OUT.nc'
  CALL CREATE_FILE(YFILE,IDIMS,YNAME_DIM,IFILE_ID,IDDIM)
  !
  YATT='J/m2'
  CALL DEF_VAR_NETCDF(IFILE_ID,'RNC_SEA'  ,'Cumulated_Averaged_Net_Radiation'        ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'HC_SEA'   ,'Cumulated_Averaged_Sensible_Heat_Flux'   ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEC_SEA'  ,'Cumulated_Averaged_Total_Latent_Heat_Flux',IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEIC_SEA' ,'Cumulated_Averaged_Sublimation_Latent_Heat_Flux',IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'GFLUXC_SEA','Cumulated_Averaged_Ground_Heat_Flux'    ,IDDIM,YATT_TITLE,YATT)
  IF(LRAD_BUDGET)THEN
    CALL DEF_VAR_NETCDF(IFILE_ID,'SWDC_SEA'  ,'Cumulated_Averaged_Downward_SW  '    ,IDDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(IFILE_ID,'SWUC_SEA'  ,'Cumulated_Averaged_Upward_SW    '    ,IDDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(IFILE_ID,'LWDC_SEA'  ,'Cumulated_Averaged_Downward_LW  '    ,IDDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(IFILE_ID,'LWUC_SEA'  ,'Cumulated_Averaged_Upward_LW     '   ,IDDIM,YATT_TITLE,YATT)
  ENDIF
  YATT='kg/ms'
  CALL DEF_VAR_NETCDF(IFILE_ID,'FMUC_SEA'  ,'Cumulated_Averaged_Zonal_Wind_Stress '  ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'FMVC_SEA'  ,'Cumulated_Averaged_Merid_Wind_Stress '  ,IDDIM,YATT_TITLE,YATT)
  !   
  CALL OL_WRITE_COORD(YFILE,IFILE_ID,IDDIM,YATT_TITLE,YNAME_DIM,YUNIT1,YUNIT2,IDIM1,YDATE,ZX,ZY)
  !
ENDIF
!
!
IF(.NOT.LPROVAR_TO_DIAG) THEN
  !
  YFILE='SEAFLUX_PROGNOSTIC.OUT.nc'
  CALL CREATE_FILE(YFILE,IDIMS,YNAME_DIM,IFILE_ID,IDDIM)
  !
  YATT='K'
  CALL DEF_VAR_NETCDF(IFILE_ID,'SST'   ,'Sea_Surface_Temperature'                    ,IDDIM,YATT_TITLE,YATT)
  !
  YATT='m'
  CALL DEF_VAR_NETCDF(IFILE_ID,'Z0SEA' ,'Roughness_Length'                           ,IDDIM,YATT_TITLE,YATT)
  !
  DO JLAYER=NOCKMIN+1,NOCKMAX
     WRITE(YPAS,'(I1.1,1X)') JLAYER
     IF (JLAYER>=10) WRITE(YPAS,'(I2.2,1X)') JLAYER
     YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
     ! Ocean temperature
     YATT='K'
     CALL DEF_VAR_NETCDF(IFILE_ID,'TEMP_OC'//YLVL,'Ocean Tempe '//YLVL ,IDDIM,YATT_TITLE,YATT)
     ! Ocean salinity
     YATT='psu'
     CALL DEF_VAR_NETCDF(IFILE_ID,'SALT_OC'//YLVL,'Ocean Salinity '//YLVL ,IDDIM,YATT_TITLE,YATT)
     ! Ocean U-current
     YATT='m/s'
     CALL DEF_VAR_NETCDF(IFILE_ID,'UCUR_OC'//YLVL,'Ocean U-cur '//YLVL ,IDDIM,YATT_TITLE,YATT)
     ! Ocean V-current
     YATT='m/s'
     CALL DEF_VAR_NETCDF(IFILE_ID,'VCUR_OC'//YLVL,'Ocean V-cur '//YLVL ,IDDIM,YATT_TITLE,YATT)
     ! Ocean TKE
     YATT='J'
     CALL DEF_VAR_NETCDF(IFILE_ID,'TKE_OC'//YLVL,'Ocean TKE '//YLVL ,IDDIM,YATT_TITLE,YATT)
     ! Ocean mixing coeff
     YATT='m2/s'
     CALL DEF_VAR_NETCDF(IFILE_ID,'KMEL_OC'//YLVL,'Ocean KMEL '//YLVL ,IDDIM,YATT_TITLE,YATT)
     ! Relaxation temperature
     YATT='K/s'
     CALL DEF_VAR_NETCDF(IFILE_ID,'T_OC_REL'//YLVL,'Ocean T rel'//YLVL ,IDDIM,YATT_TITLE,YATT)
     ! Relaxation salinity
     YATT='psu/s'
     CALL DEF_VAR_NETCDF(IFILE_ID,'S_OC_REL'//YLVL,'Ocean S rel'//YLVL ,IDDIM,YATT_TITLE,YATT)
     ! Relaxation U
     YATT='m/s2'
     CALL DEF_VAR_NETCDF(IFILE_ID,'U_OC_REL'//YLVL,'Ocean U rel'//YLVL ,IDDIM,YATT_TITLE,YATT)
     ! Relaxation V
     YATT='m/s2'
     CALL DEF_VAR_NETCDF(IFILE_ID,'V_OC_REL'//YLVL,'Ocean V rel'//YLVL ,IDDIM,YATT_TITLE,YATT)
  END DO  
!
  IF (LSBL) THEN
    ! 6.1 Heights
    YATT = 'm'
    DO JLAYER=1,NLVL
      WRITE(YPAS,'(I2.2,1X)') JLAYER ; YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
      CALL DEF_VAR_NETCDF(IFILE_ID,'SEA_SBL_Z'//YLVL,'Height_of_canopy_Layer_'//YLVL ,IDDIM,YATT_TITLE,YATT)
    END DO
    ! 6.2 Wind
    YATT = 'm/s'
    DO JLAYER=1,NLVL
      WRITE(YPAS,'(I2.2,1X)') JLAYER ; YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
      CALL DEF_VAR_NETCDF(IFILE_ID,'SEA_SBL_U'//YLVL,'Wind_of_canopy_Layer_'//YLVL ,IDDIM,YATT_TITLE,YATT)
    END DO
    ! 6.3 Temperature
    YATT = 'K'
    DO JLAYER=1,NLVL
      WRITE(YPAS,'(I2.2,1X)') JLAYER ; YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
      CALL DEF_VAR_NETCDF(IFILE_ID,'SEA_SBL_T'//YLVL,'Temperature_of_canopy_Layer_'//YLVL ,IDDIM,YATT_TITLE,YATT)
    END DO
    ! 6.4 Temperature
    YATT = 'kg/m3'
    DO JLAYER=1,NLVL
      WRITE(YPAS,'(I2.2,1X)') JLAYER ; YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
      CALL DEF_VAR_NETCDF(IFILE_ID,'SEA_SBL_Q'//YLVL,'Humidity_of_canopy_Layer_'//YLVL ,IDDIM,YATT_TITLE,YATT)
    END DO
    ! 6.5 Turbulence
    YATT = 'm2/s2'
    DO JLAYER=1,NLVL
      WRITE(YPAS,'(I2.2,1X)') JLAYER ; YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
      CALL DEF_VAR_NETCDF(IFILE_ID,'SEA_SBL_E'//YLVL,'TKE_of_canopy_Layer_'//YLVL ,IDDIM,YATT_TITLE,YATT)
    END DO
  ENDIF
  !
  CALL OL_WRITE_COORD(YFILE,IFILE_ID,IDDIM,YATT_TITLE,YNAME_DIM,YUNIT1,YUNIT2,IDIM1,YDATE,ZX,ZY)
  !
ENDIF ! End .NOT.LPROVAR_TO_DIAG
IF (LHOOK) CALL DR_HOOK('INIT_OUTFN_SEA_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE INIT_OUTFN_SEA_n
