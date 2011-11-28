!     #########
       SUBROUTINE INIT_OUTFN_ISBA_n(HPROGRAM,KLUOUT)
!     ###############################
!
!
!!****  *INIT_OUTFN_ISBA_n* -  create output files and defines variables
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
!!      F. Habets   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07-03
!!      modified    11-03, by P. Le Moigne   *Meteo France*
!!      modified    05-04, by P. Le Moigne : surf_atm diagnostics moved at the
!!                                           right place
!!      modified    10-04, by P. Le Moigne : add new diagnostics
!!      modified    10-04, by P. Le Moigne : add Halstead coefficient
!!      modified     2008, by B. Decharme  : limit the number of diag
!!                                           Add floodplains diag
!!      modified    04-09, by A.L. Gibelin : Add respiration diagnostics
!!      modified    05-09, by A.L. Gibelin : Add carbon spinup
!!      modified    07-09, by A.L. Gibelin : Add carbon prognostic variables
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIAG_SURF_ATM_n,ONLY : LPROVAR_TO_DIAG, LSELECT, CSELECT
!
USE MODD_OL_FILEID,        ONLY : XVAR_TO_FILEOUT, XID, XOUT
USE MODD_IO_SURF_OL,       ONLY : NSTEP_OUTPUT
USE MODD_ISBA_n,           ONLY : CPHOTO, CRUNOFF, CRAIN, LCANOPY, TSNOW, &
                                    LFLOOD, LGLACIER, LTEMP_ARP, CISBA,    &
                                    NTEMPLAYER_ARP, CRESPSL, TTIME 
USE MODD_DATA_COVER_PAR,   ONLY : NVEGTYPE
USE MODD_DIAG_ISBA_n
USE MODD_DIAG_EVAP_ISBA_n
USE MODD_DIAG_MISC_ISBA_n ,ONLY : LSURF_MISC_BUDGET, LWOOD_SPIN, LSOILCARB_SPIN
USE MODD_ASSIM ,           ONLY : LASSIM, CASSIM
USE MODD_AGRI  ,           ONLY : LAGRIP
USE MODD_DIAG_ISBA_n,      ONLY : LPGD
!
USE MODI_GET_DIM_FOR_2D
USE MODI_GET_COORD_FOR_2D
USE MODI_GET_DATE_OL
USE MODI_FILL_ID_OL
USE MODI_CREATE_FILE
USE MODI_DEF_VAR_NETCDF
USE MODI_GET_COORD_n
!
USE MODN_IO_OFFLINE
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_GET_DIM_FULL_n
!
USE MODI_GET_ISBA_CONF_n
!
USE MODI_GET_OFFLINE_CONF
!
IMPLICIT NONE
include 'netcdf.inc'

CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM
INTEGER, INTENT(IN) :: KLUOUT
!
!*       0.2   Declarations of local variables
!              -------------------------------
INTEGER, DIMENSION(:),ALLOCATABLE   :: JDIM, IDDIM, IDIMS
CHARACTER(LEN=100), DIMENSION(:), ALLOCATABLE :: YNAME_DIM
CHARACTER(LEN=100), DIMENSION(1) :: YATT_TITLE, YATT
CHARACTER(LEN=50)                :: YFILE
CHARACTER(LEN=3)                 :: YPAS, YLVL
CHARACTER(LEN=40),DIMENSION(1)   :: HDATE
CHARACTER(LEN=2)                 :: YLVLV
REAL                             :: ZTSTEP_OUTPUT, ZTSTEP_FORC
INTEGER                          :: IOUTPUT, IRET, INB_FORC, INI
INTEGER                          :: INPATCH, INLVLD, INLVLS, IL, INBIOMASS, &
                                      INLITTER, INLITTLEVS, INSOILCARB 
INTEGER                          :: JLAYER, JVEG, JNBIOMASS,            &
                                      JNLITTER, JNLITTLEVS, JNSOILCARB  
INTEGER                          :: IFILE_ID, IVAR_ID
INTEGER                          :: ICANLVL
INTEGER                          :: IDIM1, IDIM2, I0
CHARACTER(LEN=13),DIMENSION(1)   :: YUNIT1, YUNIT2
CHARACTER(LEN=3)                 :: HTYPE
INTEGER                          :: IRESP          ! IRESP  : return-code if a problem appears
CHARACTER(LEN=100)               :: YCOMMENT       ! Comment string
REAL,DIMENSION(:), ALLOCATABLE   :: ZX, ZY
INTEGER                          :: JRET
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!
!-------------------------------------------------------------------------------

! 1. Compute output lenght dimension
!-----------------------------------

IF (LHOOK) CALL DR_HOOK('INIT_OUTFN_ISBA_N',0,ZHOOK_HANDLE)
CALL GET_DIM_FULL_n(INI)
CALL GET_ISBA_CONF_n(INPATCH,INLVLD,INLVLS, INBIOMASS, &
                       INLITTER, INLITTLEVS, INSOILCARB)  
CALL GET_OFFLINE_CONF(ZTSTEP_OUTPUT)
!
!
ICANLVL = 6
!
CALL GET_DIM_FOR_2D(IDIM1,IDIM2,HTYPE)
IF (IDIM1.NE.0 .AND. .NOT.LWRITE_COORD) THEN
        ALLOCATE(ZX(IDIM1))
        ALLOCATE(ZY(IDIM2))
        ALLOCATE(IDIMS(4))
        ALLOCATE(IDDIM(4))
        ALLOCATE(JDIM(3))
        ALLOCATE(YNAME_DIM(4))
        CALL GET_COORD_FOR_2D(KLUOUT,INI,ZX,ZY)
        ! 2. define dimension length
        IDIMS(1)=IDIM1
        IDIMS(2)=IDIM2
        ! 3. define dimension name
        IF (HTYPE.EQ.'LON') THEN
          YNAME_DIM(1)='lon'
          YNAME_DIM(2)='lat'
          YUNIT1(1)='degrees_east'
          YUNIT2(1)='degrees_north'
        ELSE
          YNAME_DIM(1)='xx'
          YNAME_DIM(2)='yy'
          YUNIT1(1)='meters'
          YUNIT2(1)='meters'          
        ENDIF
ELSE
        IF (LWRITE_COORD) THEN
          ALLOCATE(ZX(INI))
          ALLOCATE(ZY(INI))
          CALL GET_COORD_n(HPROGRAM,INI,ZX,ZY)
        ENDIF
        ALLOCATE(IDIMS(3))
        ALLOCATE(IDDIM(3))
        ALLOCATE(JDIM(2))
        ALLOCATE(YNAME_DIM(3))
        IDIMS(1)=INI
        YNAME_DIM(1)='Number_of_points'
ENDIF
I0=SIZE(IDDIM)
IDIMS(I0-1)=INPATCH
IDIMS(I0)=NF_UNLIMITED
YNAME_DIM(I0-1)='Number_of_Tile'
YNAME_DIM(I0)='time'
!
CALL GET_DATE_OL(TTIME,XTSTEP_OUTPUT,HDATE(1))
!
! 4. Create output file for prognostic variables
!----------------------------------------------------------
!
IF (ALLOCATED(XVAR_TO_FILEOUT)) DEALLOCATE(XVAR_TO_FILEOUT)
IF (ALLOCATED(XID)) DEALLOCATE(XID)
ALLOCATE(XVAR_TO_FILEOUT(0))
ALLOCATE(XID(0))
XOUT=0
!
IF(.NOT.LPROVAR_TO_DIAG)THEN
!
  YFILE='ISBA_PROGNOSTIC.OUT.nc'
  CALL CREATE_FILE(YFILE,IDIMS,YNAME_DIM,IFILE_ID,IDDIM)
!
  YATT_TITLE(1)='units'
!
  IF(LTEMP_ARP)THEN
     IL=NTEMPLAYER_ARP
  ELSE
     IL=INLVLD
  ENDIF
  !
  IF(LWRITE_COORD)THEN
     CALL DEF_VAR_NETCDF(IFILE_ID,'LON','longitudes',IDDIM(1:1),YATT_TITLE,(/'deg'/))
     CALL DEF_VAR_NETCDF(IFILE_ID,'LAT','latitudes ',IDDIM(1:1),YATT_TITLE,(/'deg'/))
  ENDIF  
  !
  IF (IDIM1.NE.0 .AND. .NOT.LWRITE_COORD) THEN
        CALL DEF_VAR_NETCDF(IFILE_ID,TRIM(YNAME_DIM(1)),'',IDDIM(1:1),YATT_TITLE,YUNIT1)
        CALL DEF_VAR_NETCDF(IFILE_ID,TRIM(YNAME_DIM(2)),'',IDDIM(2:2),YATT_TITLE,YUNIT2)
  ENDIF
  CALL DEF_VAR_NETCDF(IFILE_ID,'time','',IDDIM(I0:I0),YATT_TITLE,HDATE)
  !
  DO JLAYER=1,IL
     WRITE(YPAS,'(i3)') JLAYER ; YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
     CALL DEF_VAR_NETCDF(IFILE_ID,'TG'//YLVL,'Soil_temp_layer_'//YLVL ,IDDIM,YATT_TITLE,(/'Kelvin'/))
  ENDDO
!
  DO JLAYER=1,INLVLD
     WRITE(YPAS,'(i3)') JLAYER ; YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
     CALL DEF_VAR_NETCDF(IFILE_ID,'WG'//YLVL,'Soil_liquid_layer_'//YLVL ,IDDIM,YATT_TITLE,(/'m3/m3'/))
     CALL DEF_VAR_NETCDF(IFILE_ID,'WGI'//YLVL,'Soil_ice'//YLVL ,IDDIM,YATT_TITLE,(/'m3/m3'/))
  ENDDO
  !
  CALL DEF_VAR_NETCDF(IFILE_ID,'WR','Interception_reservoir',IDDIM,YATT_TITLE,(/'mm'/))
  CALL DEF_VAR_NETCDF(IFILE_ID,'RESA','Aerodynamic_resistance',IDDIM,YATT_TITLE,(/'s/m'/))
  !
  IF(LGLACIER)THEN
     CALL DEF_VAR_NETCDF(IFILE_ID,'ICE_STO','Glacier_reservoir',IDDIM,YATT_TITLE,(/'Kg/m2'/))
  ENDIF
  !
  DO JLAYER=1,INLVLS
     WRITE(YPAS,'(i3)') JLAYER ; YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
     CALL DEF_VAR_NETCDF(IFILE_ID,'WSNOW_VEG'//YLVL,'Snow_Water_Equivalent_layer_'//YLVL ,IDDIM,YATT_TITLE,(/'Kg/m2'/))
     CALL DEF_VAR_NETCDF(IFILE_ID,'RSNOW_VEG'//YLVL,'Snow_density_layer_'//YLVL ,IDDIM,YATT_TITLE,(/'Kg/m3'/))
    IF (TSNOW%SCHEME=='3-L' .OR. TSNOW%SCHEME=='CRO') THEN   
      CALL DEF_VAR_NETCDF(IFILE_ID,'HSNOW_VEG'//YLVL,'Snow_heat_layer'//YLVL ,IDDIM,YATT_TITLE,(/'J/m2'/))
    ELSE
      CALL DEF_VAR_NETCDF(IFILE_ID,'TSNOW_VEG'//YLVL,'Snow_temp_layer'//YLVL ,IDDIM,YATT_TITLE,(/'K'/))
    ENDIF
    IF (TSNOW%SCHEME=='CRO') THEN   
       CALL DEF_VAR_NETCDF(IFILE_ID,'SGRAN1_VEG'//YLVL,'Snow_grain_parameter1_layer_'//YLVL ,IDDIM,YATT_TITLE,(/'-'/))
       CALL DEF_VAR_NETCDF(IFILE_ID,'SGRAN2_VEG'//YLVL,'Snow_grain_parameter2_layer_'//YLVL ,IDDIM,YATT_TITLE,(/'-'/))
       CALL DEF_VAR_NETCDF(IFILE_ID,'SHIST_VEG'//YLVL,'Snow_historical_param_layer_'//YLVL ,IDDIM,YATT_TITLE,(/'-'/))
       CALL DEF_VAR_NETCDF(IFILE_ID,'SAGE_VEG'//YLVL,'Snow_age_param_layer_'//YLVL ,IDDIM,YATT_TITLE,(/'-'/))
    ENDIF
  ENDDO
  !
  CALL DEF_VAR_NETCDF(IFILE_ID,'ASNOW_VEG','Snow_albedo',IDDIM,YATT_TITLE,(/'-'/))
  !
  IF (CPHOTO /= 'NON') THEN
     CALL DEF_VAR_NETCDF(IFILE_ID,'AN','Net CO2 Assimilation',IDDIM,YATT_TITLE,(/'kgCO2/kgair m/s'/))
     CALL DEF_VAR_NETCDF(IFILE_ID,'ANFM','Leaf CO2 Assimilation',IDDIM,YATT_TITLE,(/'kgCO2/kgair m/s'/))
     CALL DEF_VAR_NETCDF(IFILE_ID,'ANDAY','Daily Net CO2 Assimilation',IDDIM,YATT_TITLE,(/'kgCO2/m2/day'/))
  ENDIF
  !
  IF (CPHOTO == 'NIT' .OR. CPHOTO == 'NCB') THEN
    DO JNBIOMASS=1,INBIOMASS
      WRITE(YPAS,'(I3)') JNBIOMASS ; YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
      CALL DEF_VAR_NETCDF(IFILE_ID,'BIOMASS'//YLVL,'Plant biomass'//YLVL,IDDIM,YATT_TITLE,(/'kg/m2'/))
    END DO
  ENDIF
  !
  IF (CRESPSL=='CNT') THEN
  ! 
    DO JNLITTER=1,INLITTER
      DO JNLITTLEVS=1,INLITTLEVS
        WRITE(YPAS,'(I1,A1,I1)') JNLITTER,'_',JNLITTLEVS
        YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
        CALL DEF_VAR_NETCDF(IFILE_ID,'LITTER'//YLVL,'Litter pool'//YLVL,IDDIM,YATT_TITLE,(/'gC/m2'/))
      END DO
    END DO  
  !
    DO JNSOILCARB=1,INSOILCARB
      WRITE(YPAS,'(I3)') JNSOILCARB ; YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
      CALL DEF_VAR_NETCDF(IFILE_ID,'SOILCARB'//YLVL,'Soil carbon pool'//YLVL,IDDIM,YATT_TITLE,(/'gC/m2'/))
    END DO
  !
    DO JNLITTLEVS=1,INLITTLEVS
      WRITE(YPAS,'(I3)') JNLITTLEVS ; YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
      CALL DEF_VAR_NETCDF(IFILE_ID,'LIGNIN_STR'//YLVL,'Ratio Lignin/Carbon in structural litter'//YLVL,IDDIM,YATT_TITLE,(/'gC/m2'/))
    END DO
  !
  ENDIF


  IF (LCANOPY) THEN
    IF (IDIM1.NE.0 .AND. .NOT.LWRITE_COORD) THEN
      JDIM(1)=IDDIM(1)
      JDIM(2)=IDDIM(2)
      JDIM(3)=IDDIM(4)
    ELSE
      JDIM(1)=IDDIM(1)
      JDIM(2)=IDDIM(3)
    ENDIF
  DO JLAYER=1,ICANLVL
    WRITE(YLVLV,'(i2.2)') JLAYER
    CALL DEF_VAR_NETCDF(IFILE_ID,'ISBA_CAN_Z'//YLVLV,'Canopy height',JDIM,YATT_TITLE,(/'m'/))
    CALL DEF_VAR_NETCDF(IFILE_ID,'ISBA_CAN_U'//YLVLV,'Canopy wind',JDIM,YATT_TITLE,(/'m/s'/))
    CALL DEF_VAR_NETCDF(IFILE_ID,'ISBA_CAN_T'//YLVLV,'Canopy temp',JDIM,YATT_TITLE,(/'K'/))
    CALL DEF_VAR_NETCDF(IFILE_ID,'ISBA_CAN_Q'//YLVLV,'Canopy humidity',JDIM,YATT_TITLE,(/'kg/m3'/))
    CALL DEF_VAR_NETCDF(IFILE_ID,'ISBA_CAN_E'//YLVLV,'Canopy TKE',JDIM,YATT_TITLE,(/'m2/s2'/))
  END DO
  ENDIF
!
  JRET=NF_OPEN(YFILE,NF_WRITE,IFILE_ID)
  CALL FILL_ID_OL(IFILE_ID)
  IF (IDIM1.NE.0 .AND. .NOT.LWRITE_COORD) THEN
    JRET=NF_INQ_VARID   (IFILE_ID,YNAME_DIM(1),IVAR_ID)
    JRET=NF_PUT_VAR_DOUBLE(IFILE_ID,IVAR_ID,ZX)
    JRET=NF_INQ_VARID   (IFILE_ID,YNAME_DIM(2),IVAR_ID)
    JRET=NF_PUT_VAR_DOUBLE(IFILE_ID,IVAR_ID,ZY)
  ELSEIF (LWRITE_COORD) THEN
    JRET=NF_INQ_VARID   (IFILE_ID,'LON',IVAR_ID)
    JRET=NF_PUT_VAR_DOUBLE(IFILE_ID,IVAR_ID,ZX)
    JRET=NF_INQ_VARID   (IFILE_ID,'LAT',IVAR_ID)
    JRET=NF_PUT_VAR_DOUBLE(IFILE_ID,IVAR_ID,ZY)
  ENDIF
ENDIF
!
! 4. Create output file for fluxes values
!----------------------------------------------------------

YFILE='ISBA_DIAGNOSTICS.OUT.nc'
CALL CREATE_FILE(YFILE,IDIMS,YNAME_DIM,IFILE_ID,IDDIM)
YATT_TITLE(1)='units'
YATT      (1)='dimensionless'

IF (IDIM1.NE.0 .AND. .NOT.LWRITE_COORD) THEN
  JDIM(1)=IDDIM(1)
  JDIM(2)=IDDIM(2)  
  JDIM(3)=IDDIM(4)
ELSE
  JDIM(1)=IDDIM(1)
  JDIM(2)=IDDIM(3)
ENDIF

IF (IDIM1.NE.0 .AND. .NOT.LWRITE_COORD) THEN
      CALL DEF_VAR_NETCDF(IFILE_ID,TRIM(YNAME_DIM(1)),'',IDDIM(1:1),YATT_TITLE,YUNIT1)
      CALL DEF_VAR_NETCDF(IFILE_ID,TRIM(YNAME_DIM(2)),'',IDDIM(2:2),YATT_TITLE,YUNIT2)
ENDIF
CALL DEF_VAR_NETCDF(IFILE_ID,'time','',IDDIM(I0:I0),YATT_TITLE,HDATE)

IF(LWRITE_COORD)THEN
   CALL DEF_VAR_NETCDF(IFILE_ID,'LON','longitudes',JDIM(1:1),YATT_TITLE,(/'deg'/))
   CALL DEF_VAR_NETCDF(IFILE_ID,'LAT','latitudes ',JDIM(1:1),YATT_TITLE,(/'deg'/))
ENDIF

IF (LCOEF) THEN
   YATT='W/s2'
   CALL DEF_VAR_NETCDF(IFILE_ID,'CD_ISBA'   ,'Drag_Coefficient_For_Momentum   '   ,JDIM,YATT_TITLE,YATT)
   YATT='W/s'
   CALL DEF_VAR_NETCDF(IFILE_ID,'CH_ISBA'   ,'Drag_Coefficient_For_Heat       '   ,JDIM,YATT_TITLE,YATT)
   YATT='W/s/K'
   CALL DEF_VAR_NETCDF(IFILE_ID,'CE_ISBA'   ,'Drag_Coefficient_For_Evaporation'   ,JDIM,YATT_TITLE,YATT)
   YATT='m'
   CALL DEF_VAR_NETCDF(IFILE_ID,'Z0_ISBA'   ,'Roughness_Length_For_Momentum'   ,JDIM,YATT_TITLE,YATT)
   YATT='m'
   CALL DEF_VAR_NETCDF(IFILE_ID,'Z0H_ISBA'  ,'Roughness_Length_For_Heat'   ,JDIM,YATT_TITLE,YATT)
ENDIF

IF (LSURF_VARS) THEN
   YATT='kg/kg'
   CALL DEF_VAR_NETCDF(IFILE_ID,'QS_ISBA'   ,'Surface_Humidity   '   ,JDIM,YATT_TITLE,YATT)
ENDIF

IF (N2M>0) THEN
   YATT='dimensionless'
   CALL DEF_VAR_NETCDF(IFILE_ID,'RI_ISBA'   ,'Averaged_Richardson_Number'     ,JDIM,YATT_TITLE,YATT) 
ENDIF
!
IF (N2M>0) THEN
   YATT='K'
   CALL DEF_VAR_NETCDF(IFILE_ID,'T2M_ISBA' ,'2m_Temperature         '   ,JDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'T2MMIN_ISBA' ,'Minimum_2m_Temperature         '   ,JDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'T2MMAX_ISBA' ,'Maximum_2m_Temperature         '   ,JDIM,YATT_TITLE,YATT)
   YATT='kg/kg'
   CALL DEF_VAR_NETCDF(IFILE_ID,'Q2M_ISBA' ,'2m_Specific_Humidity   '   ,JDIM,YATT_TITLE,YATT)
   YATT='(-)'
   CALL DEF_VAR_NETCDF(IFILE_ID,'HU2M_ISBA','2m_Relative_Humidity   '   ,JDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'HU2MMIN_ISBA','Minimum_2m_Relative_Humidity   '   ,JDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'HU2MMAX_ISBA','Maximum_2m_Relative_Humidity   '   ,JDIM,YATT_TITLE,YATT)
   YATT='m/s'
   CALL DEF_VAR_NETCDF(IFILE_ID,'ZON10M_ISBA','10m_Zonal_wind       '   ,JDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'MER10M_ISBA','10m_Meridian_Wind     '   ,JDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'W10M_ISBA','10m_Wind     '   ,JDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'W10MMAX_ISBA','Maximum_10m_Wind     '   ,JDIM,YATT_TITLE,YATT)
ENDIF

IF (LSURF_BUDGET)  THEN
  YATT='W/m2'
  CALL DEF_VAR_NETCDF(IFILE_ID,'RN_ISBA'     ,'Averaged_Net_Radiation'                                 ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'H_ISBA'      ,'Averaged_Sensible_Heat_Flux'                            ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LE_ISBA'     ,'Averaged_Total_Latent_Heat_Flux  '                      ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEI_ISBA'    ,'Averaged_Sublimation_Latent_Heat_Flux  '                ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'GFLUX_ISBA'  ,'Averaged_Ground_Heat_Flux  '                            ,JDIM,YATT_TITLE,YATT)
  IF(LRAD_BUDGET)THEN
    CALL DEF_VAR_NETCDF(IFILE_ID,'SWD_ISBA'    ,'Averaged_Downward_SW       '                          ,JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(IFILE_ID,'SWU_ISBA'    ,'Averaged_Upward_SW         '                          ,JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(IFILE_ID,'LWD_ISBA'    ,'Averaged_Downward_LW       '                          ,JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(IFILE_ID,'LWU_ISBA'    ,'Averaged_Upward_LW         '                          ,JDIM,YATT_TITLE,YATT)
  ENDIF
  YATT='kg/ms2'
  CALL DEF_VAR_NETCDF(IFILE_ID,'FMU_ISBA'    ,'Averaged_Zonal_Wind_Stress '                        ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'FMV_ISBA'    ,'Averaged_Merid_Wind_Stress '                        ,JDIM,YATT_TITLE,YATT)
ENDIF

IF (LPATCH_BUDGET.AND.LAGRIP .AND. (CPHOTO=='NIT' .OR. CPHOTO=='LAI' .OR. CPHOTO=='LST' .OR. CPHOTO=='NCB') ) THEN
  CALL DEF_VAR_NETCDF(IFILE_ID,'IRRISEUIL' ,'Irrigation_Threshold' ,IDDIM,YATT_TITLE,YATT)
  !DO JLAYER=1,INLVLD
  !   WRITE(YPAS,'(I3)') JLAYER ; YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
  !   CALL DEF_VAR_NETCDF(IFILE_ID,'SWI'//YLVL,'Soil_Wetness_Index'//YLVL ,IDDIM,YATT_TITLE,(/'-'/))
  !   CALL DEF_VAR_NETCDF(IFILE_ID,'TSWI'//YLVL,'Total_SWI_(solid+liquid)'//YLVL ,IDDIM,YATT_TITLE,(/'-'/))
  !ENDDO
ENDIF

IF (LSURF_EVAP_BUDGET) THEN

  IF(LPATCH_BUDGET)THEN      
  YATT='W/m2'
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEG'         ,'Ground_Evaporation_Heat_Flux'                             ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEGI'        ,'Soil_Ice_Sublimation'                                     ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEV'         ,'Vegetation_Evaporation_Heat_Flux'                         ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LES'         ,'Snow_Evaporation_Heat_Flux'                               ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LER'         ,'Canopy_Water_Interception_Evaporation'                    ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LETR'        ,'Vegetation_Evapotranspiration'                            ,IDDIM,YATT_TITLE,YATT)
  YATT='kg/m2s'
  CALL DEF_VAR_NETCDF(IFILE_ID,'EVAP'        ,'Evapotranspiration'                                       ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'DRAIN'       ,'Soil_Drainage_Flux'                                       ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'RUNOFF'      ,'Supersaturation_Runoff'                                   ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'HORTON'      ,'Horton_Surface_Runoff'                                    ,IDDIM,YATT_TITLE,YATT)  

  CALL DEF_VAR_NETCDF(IFILE_ID,'DRIVEG'      ,'Dripping_from_the_vegetation_reservoir'                  ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'SNOMLT'      ,'Snow_melt_flux'                                          ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'RRVEG'       ,'Precipitation_Intercepted_by_Vegetation'                 ,IDDIM,YATT_TITLE,YATT)
  IF(LFLOOD)THEN
  CALL DEF_VAR_NETCDF(IFILE_ID,'IFLOOD'      ,'Floodplains_infiltration'                                ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'PFLOOD'      ,'Precipitation_intercepted_by_the_floodplains'            ,IDDIM,YATT_TITLE,YATT)
  YATT='W/m2'
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEF'         ,'Floodplains_evaporation_Heat_Flux'                       ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEIF'        ,'Floodplains_Frozen_evaporation_Heat_Flux'                ,IDDIM,YATT_TITLE,YATT)
  ENDIF
  ENDIF
  YATT='W/m2'
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEG_ISBA'    ,'Averaged_Ground_Evaporation_Heat_Flux'                    ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEGI_ISBA'   ,'Averaged_Soil_Ice_Sublimation'                            ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEV_ISBA'    ,'Averaged_Vegetation_Evaporation_Heat_Flux'                ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LES_ISBA'    ,'Averaged_Snow_Evaporation_Heat_Flux'                      ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LER_ISBA'    ,'Averaged_Canopy_Water_Interception_Evaporation'           ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LETR_ISBA'   ,'Averaged_Vegetation_Evapotranspiration'                   ,JDIM,YATT_TITLE,YATT)
  YATT='kg/m2s'
  CALL DEF_VAR_NETCDF(IFILE_ID,'EVAP_ISBA'   ,'Averaged_Evapotranspiration'                              ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'DRAIN_ISBA'  ,'Averaged_Soil_Drainage_Flux'                              ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'RUNOFF_ISBA' ,'Averaged_Supersaturation_Runoff'                          ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'HORTON_ISBA' ,'Averaged_Horton_Surface_Runoff'                           ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'DRIVEG_ISBA'  ,'Averaged_Dripping_from_the_vegetation_reservoir'        ,JDIM,YATT_TITLE,YATT)


  CALL DEF_VAR_NETCDF(IFILE_ID,'RRVEG_ISBA'   ,'Averaged_Precipitation_Intercepted_by_Vegetation'       ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'SNOMLT_ISBA'  ,'Averaged_Snow_melt_flux'                                ,JDIM,YATT_TITLE,YATT)
  IF(LFLOOD)THEN
  CALL DEF_VAR_NETCDF(IFILE_ID,'IFLOOD_ISBA'  ,'Averaged_Floodplains_infiltration'                      ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'PFLOOD_ISBA'  ,'Averaged_Precipitation_intercepted_by_the floodplains'  ,JDIM,YATT_TITLE,YATT)
  YATT='W/m2'
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEF_ISBA'     ,'Averaged_Floodplains_evaporation_Heat_Flux'             ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEIF_ISBA'    ,'Averaged_Floodplains_Frozen_evaporation_Heat_Flux'      ,JDIM,YATT_TITLE,YATT)
  ENDIF
ENDIF
!
!
IF (TSNOW%SCHEME=='3-L' .OR. TSNOW%SCHEME=='CRO') THEN   
   YATT='K'
   IF (LSURF_MISC_BUDGET) THEN
     CALL DEF_VAR_NETCDF(IFILE_ID,'TS_ISBA' ,'Surface_Temperature_(isba+snow3l)    '   ,JDIM,YATT_TITLE,YATT)
   ENDIF
   IF (LPATCH_BUDGET)THEN
      CALL DEF_VAR_NETCDF(IFILE_ID,'TS_PATCH' ,'Surface_Temperature_(isba+snow3l)'   ,IDDIM,YATT_TITLE,YATT)
   ENDIF
ENDIF
!
IF (LPATCH_BUDGET.AND.LSURF_BUDGET)  THEN
  YATT='W/m2'
  CALL DEF_VAR_NETCDF(IFILE_ID,'H_PATCH'         ,'Sensible_Heat_Flux'                          ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LE_PATCH'        ,'Total_Latent_Heat_Flux'                      ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEI_PATCH'       ,'Sublimatiob_Latent_Heat_Flux'                ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'GFLUX_PATCH'     ,'Ground_Heat_Flux'                            ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'RN_PATCH'        ,'Net_Radiation'                               ,IDDIM,YATT_TITLE,YATT)
  IF(LRAD_BUDGET)THEN
    CALL DEF_VAR_NETCDF(IFILE_ID,'SWD_PATCH'    ,'Downward_SW       '                            ,IDDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(IFILE_ID,'SWU_PATCH'    ,'Upward_SW         '                            ,IDDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(IFILE_ID,'LWD_PATCH'    ,'Downward_LW       '                            ,IDDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(IFILE_ID,'LWU_PATCH'    ,'Upward_LW         '                            ,IDDIM,YATT_TITLE,YATT)
  ENDIF
  YATT='kg/ms2'
  CALL DEF_VAR_NETCDF(IFILE_ID,'FMU_PATCH'      ,'Zonal_Wind_Stress '                            ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'FMV_PATCH'      ,'Merid_Wind_Stress '                            ,IDDIM,YATT_TITLE,YATT)
ENDIF

IF (LSURF_MISC_BUDGET) THEN

  YATT      (1)='-'
  CALL DEF_VAR_NETCDF(IFILE_ID,'HV_ISBA'      ,'Halstead_coefficient                  '       ,JDIM,YATT_TITLE,YATT)
  YATT      (1)='m'
  CALL DEF_VAR_NETCDF(IFILE_ID,'Z0REL'        ,'Output_Z0REL'                                ,JDIM(1:2),YATT_TITLE,YATT)
  DO JLAYER=1,INLVLS
    WRITE(YPAS,'(i3)') JLAYER ; YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
    IF (TSNOW%SCHEME=='3-L') THEN   
       CALL DEF_VAR_NETCDF(IFILE_ID,'SNOWTEMP'//YLVL,'Snow_Temp_layer'//YLVL ,IDDIM,YATT_TITLE,(/'K'/))
       CALL DEF_VAR_NETCDF(IFILE_ID,'SNOWLIQ'//YLVL,'Snow_liquid_water_layer_'//YLVL ,IDDIM,YATT_TITLE,(/'m'/))
    ELSE
       CALL DEF_VAR_NETCDF(IFILE_ID,'TSNOW_VEG'//YLVL,'Snow_Temp_layer'//YLVL ,IDDIM,YATT_TITLE,(/'K'/))
    ENDIF
  ENDDO

  DO JLAYER=1,INLVLD
     WRITE(YPAS,'(I3)') JLAYER ; YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
     CALL DEF_VAR_NETCDF(IFILE_ID,'SWI'//YLVL(:LEN_TRIM(YLVL))//'_ISBA', &
                                    'Soil_Wetness_Index'//YLVL ,JDIM,YATT_TITLE,(/'-'/))  
     CALL DEF_VAR_NETCDF(IFILE_ID,'TSWI'//YLVL(:LEN_TRIM(YLVL))//'_ISBA', &
                                    'Total_SWI_(liquid+solid)'//YLVL ,JDIM,YATT_TITLE,(/'-'/))  
     IF(LPATCH_BUDGET)THEN      
     CALL DEF_VAR_NETCDF(IFILE_ID,'SWI'//YLVL,'Soil_Wetness_Index'//YLVL ,IDDIM,YATT_TITLE,(/'-'/))
     CALL DEF_VAR_NETCDF(IFILE_ID,'TSWI'//YLVL,'Total_SWI_(liquid+solid)'//YLVL ,IDDIM,YATT_TITLE,(/'-'/))
     ENDIF
  ENDDO

  YATT      (1)='-'
  CALL DEF_VAR_NETCDF(IFILE_ID,'TSWI_T_ISBA'    ,'total_swi_over_entire_soil '                  ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'PSNG_ISBA'      ,'Snow_frac_over_ground      '                  ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'PSNV_ISBA'      ,'Snow_frac_over_veg         '                  ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'PSN_ISBA '      ,'Snow_fraction              '                  ,JDIM,YATT_TITLE,YATT)
!  
  CALL DEF_VAR_NETCDF(IFILE_ID,'TALB_ISBA'      ,'Surface total albedo       '                  ,JDIM,YATT_TITLE,YATT)
!
  YATT      (1)='kg/m2'
  CALL DEF_VAR_NETCDF(IFILE_ID,'WGTOT_T_ISBA'  ,'Total_soil_water_reservoir_(liquid+solid)'      ,JDIM,YATT_TITLE,YATT)
  YATT      (1)='kg/m2'
  CALL DEF_VAR_NETCDF(IFILE_ID,'WGI_T_ISBA'    ,'Total_soil_ice_reservoir'                      ,JDIM,YATT_TITLE,YATT)
  YATT      (1)='kg/m2'
  CALL DEF_VAR_NETCDF(IFILE_ID,'WSNOW_T_ISBA'    ,'Total_snow_reservoir '              ,JDIM,YATT_TITLE,YATT)
  YATT      (1)='m'
  CALL DEF_VAR_NETCDF(IFILE_ID,'DSNOW_T_ISBA'    ,'Total_snow_depth '                  ,JDIM,YATT_TITLE,YATT)
  YATT      (1)='K'
  CALL DEF_VAR_NETCDF(IFILE_ID,'TSNOW_T_ISBA'    ,'Total_snow_temperature '            ,JDIM,YATT_TITLE,YATT)  
! 
  IF(CRUNOFF=='SGH ')THEN
    YATT (1)='-'
    CALL DEF_VAR_NETCDF(IFILE_ID,'FSAT_ISBA'      ,'Topmodel_saturated_grid-cell_fraction_(SGH)'   ,JDIM,YATT_TITLE,YATT) 
  ENDIF 
  IF(CRAIN=='SGH ')THEN
    YATT (1)='-'
    CALL DEF_VAR_NETCDF(IFILE_ID,'MUF_ISBA'       ,'Fraction_of_rainfall_reaching_the ground_(SGH)',JDIM,YATT_TITLE,YATT)  
  ENDIF  

  IF (CPHOTO /= 'NON'.AND.LPATCH_BUDGET) THEN
    YATT (1)='(kgCO2/m2/s)'
    CALL DEF_VAR_NETCDF(IFILE_ID,'GPP'        ,'Gross Primary Production',IDDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(IFILE_ID,'RESP_AUTO'  ,'Autotrophic Respiration ',IDDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(IFILE_ID,'RESP_ECO'   ,'Ecosystem Respiration   ',IDDIM,YATT_TITLE,YATT)          
  ENDIF

  IF(LFLOOD)THEN
    YATT      (1)='-'
    CALL DEF_VAR_NETCDF(IFILE_ID,'FFG_ISBA'      ,'flood_frac_over_ground      '                  ,JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(IFILE_ID,'FFV_ISBA'      ,'flood_frac_over_veg         '                  ,JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(IFILE_ID,'FF_ISBA '      ,'flood_fraction              '                  ,JDIM,YATT_TITLE,YATT)
    YATT (1)='-'
    CALL DEF_VAR_NETCDF(IFILE_ID,'FFLOOD_ISBA'  ,'Potential_floodplain_grid-cell_fraction'   ,JDIM,YATT_TITLE,YATT) 
    YATT (1)='m'
    CALL DEF_VAR_NETCDF(IFILE_ID,'HF_ISBA'      ,'Floodplain_height'               ,JDIM,YATT_TITLE,YATT)
    YATT (1)='kg/m2/s'
    CALL DEF_VAR_NETCDF(IFILE_ID,'IPF_ISBA'     ,'Potential_floodplain_infiltration',JDIM,YATT_TITLE,YATT) 
  ENDIF
  
ENDIF

IF (CPHOTO == 'NCB' .AND. LWOOD_SPIN) THEN

  DO JNBIOMASS=1,INBIOMASS
    WRITE(YPAS,'(I3)') JNBIOMASS ; YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
    CALL DEF_VAR_NETCDF(IFILE_ID,'INCREASE'//YLVL,'Biomass_increase_'//YLVL,IDDIM,YATT_TITLE,(/'kg/m2/day'/))
  ENDDO

ENDIF

IF (CRESPSL == 'CNT' .AND. LSOILCARB_SPIN) THEN

  DO JNBIOMASS=1,INBIOMASS
    WRITE(YPAS,'(I3)') JNBIOMASS ; YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
    CALL DEF_VAR_NETCDF(IFILE_ID,'TURNOVER'//YLVL,'Biomass_turnover_'//YLVL,IDDIM,YATT_TITLE,(/'gC/m2/s'/))
  ENDDO

ENDIF
!
IF(LPROVAR_TO_DIAG)THEN
!
  YATT_TITLE(1)='units'
!
  IF(LTEMP_ARP)THEN
     IL=NTEMPLAYER_ARP
  ELSEIF(CISBA/='DIFF')THEN
     IL=INLVLD-1
  ELSE
     IL=INLVLD
  ENDIF
!
  YATT='K'
!  
  DO JLAYER=1,IL
     WRITE(YPAS,'(i3)') JLAYER
     YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
     CALL DEF_VAR_NETCDF(IFILE_ID,'TG'//YLVL(:LEN_TRIM(YLVL))//'_ISBA','Soil_temp_layer_'//YLVL ,JDIM,YATT_TITLE,YATT)
  ENDDO
!
  YATT='kg/m2'
!
  DO JLAYER=1,INLVLD
     WRITE(YPAS,'(i3)') JLAYER
     YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
     CALL DEF_VAR_NETCDF(IFILE_ID,'WG'//YLVL(:LEN_TRIM(YLVL))//'_ISBA','Soil_liquid_layer_'//YLVL ,JDIM,YATT_TITLE,YATT)
  ENDDO
!  
  IF(CISBA/='DIFF')THEN
     IL=INLVLD-1
  ELSE
     IL=INLVLD
  ENDIF
  DO JLAYER=1,IL
     WRITE(YPAS,'(i3)') JLAYER
     YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
     CALL DEF_VAR_NETCDF(IFILE_ID,'WGI'//YLVL(:LEN_TRIM(YLVL))//'_ISBA','Soil_ice_layer_'//YLVL ,JDIM,YATT_TITLE,YATT)
  ENDDO
!
  CALL DEF_VAR_NETCDF(IFILE_ID,'WR_ISBA','Interception_reservoir',JDIM,YATT_TITLE,YATT) 
!
  IF(LGLACIER)THEN
     CALL DEF_VAR_NETCDF(IFILE_ID,'ICE_STO_ISBA','Glacier_reservoir',JDIM,YATT_TITLE,YATT)
  ENDIF
!
  YATT='-'
  CALL DEF_VAR_NETCDF(IFILE_ID,'ASNOW_ISBA','Snow_Albedo',JDIM,YATT_TITLE,YATT)
!
  IF(TSNOW%SCHEME=='3-L'  .OR. TSNOW%SCHEME=='CRO')THEN
!
    DO JLAYER=1,INLVLS
!    
       WRITE(YPAS,'(i3)') JLAYER
       YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
!
       YATT='kg/m2'
!     
       CALL DEF_VAR_NETCDF(IFILE_ID,'WSNOW_'//YLVL(:LEN_TRIM(YLVL))//'_ISBA', &
                             'Snow_Water_Equivalent_layer_'//YLVL,JDIM,YATT_TITLE,YATT)  
!
       YATT='m'
!                         
       CALL DEF_VAR_NETCDF(IFILE_ID,'DSNOW_'//YLVL(:LEN_TRIM(YLVL))//'_ISBA', &
                             'Snow_Depth_layer_'//YLVL,JDIM,YATT_TITLE,YATT)  
!
        YATT='K'                        
        CALL DEF_VAR_NETCDF(IFILE_ID,'TSNOW_'//YLVL(:LEN_TRIM(YLVL))//'_ISBA', &
                              'Snow_Temperature_layer_'//YLVL,JDIM,YATT_TITLE,YATT)  
!                            
    ENDDO
!    
  ENDIF
!
ENDIF    
!
JRET=NF_OPEN(YFILE,NF_WRITE,IFILE_ID)
CALL FILL_ID_OL(IFILE_ID)
IF (IDIM1.NE.0 .AND. .NOT.LWRITE_COORD) THEN
  JRET=NF_INQ_VARID   (IFILE_ID,YNAME_DIM(1),IVAR_ID)
  JRET=NF_PUT_VAR_DOUBLE(IFILE_ID,IVAR_ID,ZX)
  JRET=NF_INQ_VARID   (IFILE_ID,YNAME_DIM(2),IVAR_ID)
  JRET=NF_PUT_VAR_DOUBLE(IFILE_ID,IVAR_ID,ZY)
ELSEIF (LWRITE_COORD) THEN
  JRET=NF_INQ_VARID   (IFILE_ID,'LON',IVAR_ID)
  JRET=NF_PUT_VAR_DOUBLE(IFILE_ID,IVAR_ID,ZX)
  JRET=NF_INQ_VARID   (IFILE_ID,'LAT',IVAR_ID)
  JRET=NF_PUT_VAR_DOUBLE(IFILE_ID,IVAR_ID,ZY)
ENDIF
!
!
IF (IDIM1.NE.0 .AND. LWRITE_COORD) THEN
  DEALLOCATE(ZX)
  DEALLOCATE(ZY)
  ALLOCATE(ZX(IDIM1))
  ALLOCATE(ZY(IDIM2))
  CALL GET_COORD_FOR_2D(KLUOUT,INI,ZX,ZY)
  DEALLOCATE(IDIMS)
  DEALLOCATE(IDDIM)
  DEALLOCATE(JDIM)
  DEALLOCATE(YNAME_DIM)
  ALLOCATE(IDIMS(4))
  ALLOCATE(IDDIM(4))
  ALLOCATE(JDIM(3))
  ALLOCATE(YNAME_DIM(4))
  ! 2. define dimension length
  IDIMS(1)=IDIM1
  IDIMS(2)=IDIM2
  ! 3. define dimension name
  IF (HTYPE.EQ.'LON') THEN
    YNAME_DIM(1)='lon'
    YNAME_DIM(2)='lat'
    YUNIT1(1)='degrees_east'
    YUNIT2(1)='degrees_north'
  ELSE
    YNAME_DIM(1)='xx'
    YNAME_DIM(2)='yy'
    YUNIT1(1)='meters'
    YUNIT2(1)='meters'          
  ENDIF
  IDIMS(3)=INPATCH
  IDIMS(4)=NF_UNLIMITED
  YNAME_DIM(3)='Number_of_Tile'
  YNAME_DIM(4)='time'
  I0=4
ENDIF
!
!
IF (LSURF_BUDGETC) THEN

  YFILE='ISBA_DIAG_CUMUL.OUT.nc'
  CALL CREATE_FILE(YFILE,IDIMS,YNAME_DIM,IFILE_ID,IDDIM)
  YATT_TITLE(1)='units'
  YATT      (1)='dimensionless'

  IF (IDIM1.NE.0) THEN
    JDIM(1)=IDDIM(1)
    JDIM(2)=IDDIM(2)  
    JDIM(3)=IDDIM(4)
  ELSE
    JDIM(1)=IDDIM(1)
    JDIM(2)=IDDIM(3)
  ENDIF
!
  IF (IDIM1.NE.0) THEN
      CALL DEF_VAR_NETCDF(IFILE_ID,TRIM(YNAME_DIM(1)),'',IDDIM(1:1),YATT_TITLE,YUNIT1)
      CALL DEF_VAR_NETCDF(IFILE_ID,TRIM(YNAME_DIM(2)),'',IDDIM(2:2),YATT_TITLE,YUNIT2)
  ENDIF
  CALL DEF_VAR_NETCDF(IFILE_ID,'time','',IDDIM(I0:I0),YATT_TITLE,HDATE)
!
  IF(LPATCH_BUDGET)THEN      
  YATT='J/m2'
  CALL DEF_VAR_NETCDF(IFILE_ID,'RNC'         ,'Cumulated_Net_Radiation'                                 ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'HC'          ,'Cumulated_Sensible_Heat_Flux'                            ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEC'         ,'Cumulated_Total_Latent_Heat_Flux'                        ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEIC'        ,'Cumulated_Sublimation_Latent_Heat_Flux'                  ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'GFLUXC'      ,'Cumulated_Ground_Heat_Flux'                              ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEGC'        ,'Cumulated_Ground_Evaporation_Heat_Flux'                  ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEGIC'       ,'Cumulated_Soil_Ice_Sublimation'                          ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEVC'        ,'Cumulated_Vegetation_Evaporation_Heat_Flux'              ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LESC'        ,'Cumulated_Snow_Evaporation_Heat_Flux'                    ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LERC'        ,'Cumulated_Canopy_Water_Interception_Evaporation'         ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LETRC'       ,'Cumulated_Vegetation_Evapotranspiration'                 ,IDDIM,YATT_TITLE,YATT)
  IF(LRAD_BUDGET)THEN
    CALL DEF_VAR_NETCDF(IFILE_ID,'SWDC'      ,'Cumulated_Downward_SW       '                            ,IDDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(IFILE_ID,'SWUC'      ,'Cumulated_Upward_SW         '                            ,IDDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(IFILE_ID,'LWDC'      ,'Cumulated_Downward_LW       '                            ,IDDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(IFILE_ID,'LWUC'      ,'Cumulated_Upward_LW         '                            ,IDDIM,YATT_TITLE,YATT)
  ENDIF

  YATT='kg/ms2'
  CALL DEF_VAR_NETCDF(IFILE_ID,'FMUC'        ,'Cumulated_Zonal_Wind_Stress '                            ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'FMVC'        ,'Cumulated_Merid_Wind_Stress '                            ,IDDIM,YATT_TITLE,YATT)  
  YATT='kg/m2'  
  CALL DEF_VAR_NETCDF(IFILE_ID,'EVAPC'       ,'Cumulated_Evapotranspiration'                            ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'DRAINC'      ,'Cumulated_Soil_Drainage_Flux'                            ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'RUNOFFC'     ,'Cumulated_Supersaturation_Runoff'                        ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'HORTONC'     ,'Cumulated_Horton_Runoff'                                 ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'DRIVEGC'     ,'Cumulated_Dripping_from_the_vegetation_reservoir'        ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'SNOMLTC'     ,'Cumulated_Snow_melt_flux'                                ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'RRVEGC'      ,'Cumulated_Precipitation_Intercepted_by_Vegetation'       ,IDDIM,YATT_TITLE,YATT)
 
  IF(LGLACIER)THEN
  CALL DEF_VAR_NETCDF(IFILE_ID,'ICE_FC'      ,'Cumulated_Glacier_ice_flux'                              ,IDDIM,YATT_TITLE,YATT)
  ENDIF
  IF(LFLOOD)THEN
  CALL DEF_VAR_NETCDF(IFILE_ID,'IFLOODC'     ,'Cumulated_Floodplains_infiltration'                      ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'PFLOODC'     ,'Cumulated_Precipitation_intercepted_by_the_floodplains'  ,IDDIM,YATT_TITLE,YATT)
  YATT='J/m2'
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEFC'        ,'Cumulated_Floodplains_evaporation_Heat_Flux'             ,IDDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEIFC'       ,'Cumulated_Floodplains_Frozen_evaporation_Heat_Flux'      ,IDDIM,YATT_TITLE,YATT)
  ENDIF
  ENDIF
!  
  YATT='J/m2'
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEGC_ISBA'   ,'Averaged_Cumulated_Ground_Evaporation_Heat_Flux'          ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEGIC_ISBA'  ,'Averaged_Cumulated_Soil_Ice_Sublimation'                  ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEVC_ISBA'   ,'Averaged_Cumulated_Vegetation_Evaporation_Heat_Flux'      ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LESC_ISBA'   ,'Averaged_Cumulated_Snow_Evaporation_Heat_Flux'            ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LERC_ISBA'   ,'Averaged_Cumulated_Canopy_Water_Interception_Evaporation' ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LETRC_ISBA'  ,'Averaged_Cumulated_Vegetation_Evapotranspiration'         ,JDIM,YATT_TITLE,YATT)
  YATT='kg/m2'  
  CALL DEF_VAR_NETCDF(IFILE_ID,'EVAPC_ISBA'  ,'Averaged_Cumulated_Evapotranspiration'                    ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'DRAINC_ISBA' ,'Averaged_Cumulated_Soil_Drainage_Flux'                    ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'RUNOFFC_ISBA','Averaged_Cumulated_Supersaturation_Runoff'                ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'HORTONC_ISBA','Averaged_Cumulated_Horton_Surface_Runoff'                 ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'DRIVEGC_ISBA','Averaged_Dripping_from_the_vegetation_reservoir'          ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'SNOMLTC_ISBA','Averaged_Cumulated_Snow_melt_flux'                        ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'RRVEGC_ISBA','Averaged_Cumulated_Precipitation_Intercepted_by_Vegetation',&
         JDIM,YATT_TITLE,YATT)     

  IF(LGLACIER)THEN
  CALL DEF_VAR_NETCDF(IFILE_ID,'ICE_FC_ISBA' ,'Averaged_Cumulated_Glacier_ice_flux'                      ,JDIM,YATT_TITLE,YATT)
  ENDIF
  IF(LFLOOD)THEN
  CALL DEF_VAR_NETCDF(IFILE_ID,'IFLOODC_ISBA','Averaged_Cumulated_Floodplains_infiltration'              ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'PFLOODC_ISBA','Averaged_Cumulated_Precip_intercepted_by_the_floodplains' ,JDIM,YATT_TITLE,YATT)
  YATT='J/m2'
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEFC_ISBA'   ,'Averaged_Cumulated_Flood_evaporation_Heat_Flux'           ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEIFC_ISBA'  ,'Averaged_Cumulated_Flood_Frozen_evaporation_Heat_Flux'    ,JDIM,YATT_TITLE,YATT)
  ENDIF 

  YATT='J/m2'
  CALL DEF_VAR_NETCDF(IFILE_ID,'RNC_ISBA'    ,'Averaged_Cumulated_Net_Radiation'                         ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'HC_ISBA'     ,'Averaged_Cumulated_Sensible_Heat_Flux'                    ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEC_ISBA'    ,'Averaged_Cumulated_Total_Latent_Heat_Flux'                ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'LEIC_ISBA'   ,'Averaged_Cumulated_Sublimation_Latent_Heat_Flux'          ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'GFLUXC_ISBA' ,'Averaged_Cumulated_Ground_Heat_Flux'                      ,JDIM,YATT_TITLE,YATT)
  IF(LRAD_BUDGET)THEN
    CALL DEF_VAR_NETCDF(IFILE_ID,'SWDC_ISBA' ,'Averaged_Cumulated_Downward_SW       '                    ,JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(IFILE_ID,'SWUC_ISBA' ,'Averaged_Cumulated_Upward_SW         '                    ,JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(IFILE_ID,'LWDC_ISBA' ,'Averaged_Cumulated_Downward_LW       '                    ,JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(IFILE_ID,'LWUC_ISBA' ,'Averaged_Cumulated_Upward_LW         '                    ,JDIM,YATT_TITLE,YATT)
  ENDIF
  YATT='kg/ms2'
  CALL DEF_VAR_NETCDF(IFILE_ID,'FMUC_ISBA'   ,'Averaged_Cumulated_Zonal_Wind_Stress '                ,JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(IFILE_ID,'FMVC_ISBA'   ,'Averaged_Cumulated_Merid_Wind_Stress '                ,JDIM,YATT_TITLE,YATT)
!  
  JRET=NF_OPEN(YFILE,NF_WRITE,IFILE_ID)
  CALL FILL_ID_OL(IFILE_ID)
  IF (IDIM1.NE.0) THEN
    JRET=NF_INQ_VARID   (IFILE_ID,YNAME_DIM(1),IVAR_ID)
    JRET=NF_PUT_VAR_DOUBLE(IFILE_ID,IVAR_ID,ZX)
    JRET=NF_INQ_VARID   (IFILE_ID,YNAME_DIM(2),IVAR_ID)
    JRET=NF_PUT_VAR_DOUBLE(IFILE_ID,IVAR_ID,ZY)
  ENDIF
ENDIF


! 6. Create file for vegetation parameter values
!----------------------------------------------------------

IF(LASSIM) THEN
   IF(CASSIM=='PLUS ') THEN
      YFILE='ISBA_VEG_EVOLUTION_P.OUT.nc'
      CALL CREATE_FILE(YFILE,IDIMS,YNAME_DIM,IFILE_ID,IDDIM)
      YATT_TITLE(1)='units'
      YATT      (1)='dimensionless'
      IF (IDIM1.NE.0) THEN
        CALL DEF_VAR_NETCDF(IFILE_ID,TRIM(YNAME_DIM(1)),'',IDDIM(1:1),YATT_TITLE,YUNIT1)
        CALL DEF_VAR_NETCDF(IFILE_ID,TRIM(YNAME_DIM(2)),'',IDDIM(2:2),YATT_TITLE,YUNIT2)
      ENDIF
      CALL DEF_VAR_NETCDF(IFILE_ID,'time','',IDDIM(I0:I0),YATT_TITLE,HDATE)
      CALL DEF_VAR_NETCDF(IFILE_ID,'LAIp'   ,'Output_LAI_ISBA' ,IDDIM,YATT_TITLE,YATT)
      JRET=NF_OPEN(YFILE,NF_WRITE,IFILE_ID)
      CALL FILL_ID_OL(IFILE_ID)
      IF (IDIM1.NE.0) THEN
        JRET=NF_INQ_VARID   (IFILE_ID,YNAME_DIM(1),IVAR_ID)
        JRET=NF_PUT_VAR_DOUBLE(IFILE_ID,IVAR_ID,ZX)
        JRET=NF_INQ_VARID   (IFILE_ID,YNAME_DIM(2),IVAR_ID)
        JRET=NF_PUT_VAR_DOUBLE(IFILE_ID,IVAR_ID,ZY)
      ENDIF
   ELSEIF(CASSIM=='AVERA') THEN
      YFILE='ISBA_VEG_EVOLUTION_A.OUT.nc'
      CALL CREATE_FILE(YFILE,IDIMS,YNAME_DIM,IFILE_ID,IDDIM)
      YATT_TITLE(1)='units'
      YATT      (1)='dimensionless'
      IF (IDIM1.NE.0) THEN
        CALL DEF_VAR_NETCDF(IFILE_ID,TRIM(YNAME_DIM(1)),'',IDDIM(1:1),YATT_TITLE,YUNIT1)
        CALL DEF_VAR_NETCDF(IFILE_ID,TRIM(YNAME_DIM(2)),'',IDDIM(2:2),YATT_TITLE,YUNIT2)
      ENDIF      
      CALL DEF_VAR_NETCDF(IFILE_ID,'time','',IDDIM(I0:I0),YATT_TITLE,HDATE)
      CALL DEF_VAR_NETCDF(IFILE_ID,'LAIa'   ,'Output_LAI_ISBA' ,IDDIM,YATT_TITLE,YATT)
      JRET=NF_OPEN(YFILE,NF_WRITE,IFILE_ID)
      CALL FILL_ID_OL(IFILE_ID)
      IF (IDIM1.NE.0) THEN
        JRET=NF_INQ_VARID   (IFILE_ID,YNAME_DIM(1),IVAR_ID)
        JRET=NF_PUT_VAR_DOUBLE(IFILE_ID,IVAR_ID,ZX)
        JRET=NF_INQ_VARID   (IFILE_ID,YNAME_DIM(2),IVAR_ID)
        JRET=NF_PUT_VAR_DOUBLE(IFILE_ID,IVAR_ID,ZY)
      ENDIF
   ELSEIF(CASSIM=='2DVAR') THEN
      YFILE='ISBA_VEG_EVOLUTION.OUT.nc'
      CALL CREATE_FILE(YFILE,IDIMS,YNAME_DIM,IFILE_ID,IDDIM)
      YATT_TITLE(1)='units'
      YATT      (1)='dimensionless'
      IF (IDIM1.NE.0) THEN
        CALL DEF_VAR_NETCDF(IFILE_ID,TRIM(YNAME_DIM(1)),'',IDDIM(1:1),YATT_TITLE,YUNIT1)
        CALL DEF_VAR_NETCDF(IFILE_ID,TRIM(YNAME_DIM(2)),'',IDDIM(2:2),YATT_TITLE,YUNIT2)
      ENDIF      
      CALL DEF_VAR_NETCDF(IFILE_ID,'time','',IDDIM(I0:I0),YATT_TITLE,HDATE)
      CALL DEF_VAR_NETCDF(IFILE_ID,'LAI'   ,'Output_LAI_ISBA' ,IDDIM,YATT_TITLE,YATT)
      JRET=NF_OPEN(YFILE,NF_WRITE,IFILE_ID)
      CALL FILL_ID_OL(IFILE_ID)
      IF (IDIM1.NE.0) THEN
        JRET=NF_INQ_VARID   (IFILE_ID,YNAME_DIM(1),IVAR_ID)
        JRET=NF_PUT_VAR_DOUBLE(IFILE_ID,IVAR_ID,ZX)
        JRET=NF_INQ_VARID   (IFILE_ID,YNAME_DIM(2),IVAR_ID)
        JRET=NF_PUT_VAR_DOUBLE(IFILE_ID,IVAR_ID,ZY)
      ENDIF
   ENDIF
ELSEIF(LPGD)THEN
   YFILE='ISBA_VEG_EVOLUTION.OUT.nc'
   CALL CREATE_FILE(YFILE,IDIMS,YNAME_DIM,IFILE_ID,IDDIM)  
   YATT_TITLE(1)='units'
   YATT      (1)='dimensionless'

   IF (IDIM1.NE.0) THEN
      CALL DEF_VAR_NETCDF(IFILE_ID,TRIM(YNAME_DIM(1)),'',IDDIM(1:1),YATT_TITLE,YUNIT1)
      CALL DEF_VAR_NETCDF(IFILE_ID,TRIM(YNAME_DIM(2)),'',IDDIM(2:2),YATT_TITLE,YUNIT2)
   ENDIF
   CALL DEF_VAR_NETCDF(IFILE_ID,'time','',IDDIM(I0:I0),YATT_TITLE,HDATE)
   CALL DEF_VAR_NETCDF(IFILE_ID,'VEG'         ,'Output_vegetation_fraction'         ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'LAI'         ,'Output_LAI_ISBA'                    ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'Z0VEG'       ,'Roughness_Length_Vegetation'        ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'PATCH'       ,'Fraction_Of_Patch'                  ,IDDIM(1:I0),YATT_TITLE,YATT)
  
   DO JVEG=1,NVEGTYPE
      WRITE(YPAS,'(i2)') JVEG 
      YLVLV=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
      CALL DEF_VAR_NETCDF(IFILE_ID,'VEGTYPE_P'//YLVLV,'fraction_of_vegetation_type_'//YLVLV ,IDDIM(1:I0),YATT_TITLE,(/'-'/))
   ENDDO   

   CALL DEF_VAR_NETCDF(IFILE_ID,'EMIS_ISBA'   ,'Emissivity_Of_Vegetation'           ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'RSMIN'       ,'Minimal_Stomatal_Resistance'        ,IDDIM(1:I0),YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'GAMMA'       ,'Coefficient_Computation_Rsmin'      ,IDDIM(1:I0),YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'CV'          ,'Vegetal_Thermal_Inertia'            ,IDDIM(1:I0),YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'RGL'         ,'Max_Solar_Radiation_Photosynthesis' ,IDDIM(1:I0),YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'WRMAX_CF'    ,'Coefficient_Max_Water_Interception' ,IDDIM(1:I0),YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'ALBNIR_SOIL' ,'Output_ALBNIR_SOIL'                 ,IDDIM(1:I0),YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'ALBVIS_SOIL' ,'Output_ALBVIS_SOIL'                 ,IDDIM(1:I0),YATT_TITLE,YATT)

   IF (LAGRIP .AND. (CPHOTO=='NIT' .OR. CPHOTO=='LAI' .OR. CPHOTO=='LST' .OR. CPHOTO=='NCB') ) THEN
      CALL DEF_VAR_NETCDF(IFILE_ID,'WATSUP' ,'Water_Supply_Irrigation' ,IDDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(IFILE_ID,'IRRIG'  ,'Fraction_Of_Irrigated_Vegetation' ,IDDIM,YATT_TITLE,YATT)
      !CALL DEF_VAR_NETCDF(IFILE_ID,'SEED'  ,'Seeding_Date' ,IDDIM,YATT_TITLE,YATT)
      !CALL DEF_VAR_NETCDF(IFILE_ID,'REAP'  ,'Reaping_Date' ,IDDIM,YATT_TITLE,YATT)
   END IF

   JRET=NF_OPEN(YFILE,NF_WRITE,IFILE_ID)
   CALL FILL_ID_OL(IFILE_ID)
   IF (IDIM1.NE.0) THEN
     JRET=NF_INQ_VARID   (IFILE_ID,YNAME_DIM(1),IVAR_ID)
     JRET=NF_PUT_VAR_DOUBLE(IFILE_ID,IVAR_ID,ZX)
     JRET=NF_INQ_VARID   (IFILE_ID,YNAME_DIM(2),IVAR_ID)
     JRET=NF_PUT_VAR_DOUBLE(IFILE_ID,IVAR_ID,ZY)
    ENDIF
ENDIF
IF (LHOOK) CALL DR_HOOK('INIT_OUTFN_ISBA_N',1,ZHOOK_HANDLE)

END SUBROUTINE INIT_OUTFN_ISBA_n
