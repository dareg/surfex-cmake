!----------------------------
PROGRAM PRE_INPUT_EXPERIMENT
!----------------------------
!!
!!    PURPOSE
!!    -------
!!   This program prepares the input files for offline run:
!!   ie a file named PARAMS.nc containing information for the run
!!   and FORCING.nc which contains meteorological forcing      
!!
!!    METHOD
!!    ------
!!   
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    P. LeMoigne                  Meteo-France
!!    S. Lafont    05/2009 correct size (1:JNPTS,:) of output arrays for writing netcdf
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original     06/06
!!
!----------------------------------------------------------------------------
!
USE MODD_SURF_CONF, ONLY: CPROGNAME
USE MODD_CSTS
USE MODD_TYPE_DATE_SURF
USE MODD_DIAG_SURF_ATM_n,ONLY : N2M, L2M_MIN_ZS, LSURF_BUDGET,     &
                                LRAD_BUDGET, LCOEF, XDIAG_TSTEP,   &
                                LFRAC, LSURF_VARS, LDIAG_GRID,     &
                                LSURF_BUDGETC, LRESET_BUDGETC,     &
                                LPROVAR_TO_DIAG, LSELECT, CSELECT
USE MODE_POS_SURF
!
USE MODI_ALLOC_SURFEX
USE MODI_GOTO_SURFEX
USE MODI_DEALLOC_SURFEX
USE MODI_CREATE_FILE
USE MODI_SUNPOS
USE MODI_WRITE_SURF
USE MODI_WRITE_NETCDF
USE MODI_DEF_VAR_NETCDF
USE MODI_MY_FORC
USE MODI_OPEN_NAMELIST
USE MODI_DEFAULT_DIAG_SURF_ATM
USE MODI_GET_DATE_OL
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
!----------------------------------------------------------------------------
!
!*    0.     Declaration of parameters
!            -------------------------
!      
IMPLICIT NONE
!
include 'netcdf.inc'
!----------------------------------------------------------------------------
!
INTEGER :: INI                  ! number of forcing points
INTEGER :: INPTS                ! number of forcing time-step in input file
INTEGER :: JNPTS                ! number of forcing time-step written in netcdf file!
REAL    :: ZTSTEPFRC            ! time step of atmospheric forcing (s)
!
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZCO2      ! CO2 concentration (kg/m3) 
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZDIR_SW   ! Solar direct   radiation (W/m2)
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZSCA_SW   ! Solar diffused radiation (W/m2)
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZLW       ! Longwave radiation (W/m2)
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZWINDSPEED! Wind speed (m/s)
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZWINDDIR  ! Wind dir. (deg. from N, clockwise)
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZRAIN     ! rain rate (kg/m2/s)
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZSNOW     ! snow rate (kg/m2/s)
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZTA       ! temperature (K)
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZQA       ! humidity (kg/m3)
REAL, DIMENSION(:,:),  ALLOCATABLE :: ZPS       ! pressure (Pa)
REAL, DIMENSION(:),    ALLOCATABLE :: ZZREF     ! height of temperature forcing (m)
REAL, DIMENSION(:),    ALLOCATABLE :: ZUREF     ! height of wind forcing (m)
REAL, DIMENSION(:),    ALLOCATABLE :: ZZS       ! orography (m)
REAL, DIMENSION(:),    ALLOCATABLE :: ZLON      ! longitude (degrees)
REAL, DIMENSION(:),    ALLOCATABLE :: ZLAT      ! latitude  (degrees)
!
REAL, DIMENSION(:),    ALLOCATABLE :: ZTIME
!
!----------------------------------------------------------------------------
!      
!*    1.     Declaration of variables
!            ------------------------
!
CHARACTER(LEN=9), PARAMETER       :: YFILE_PARA    = 'PARAMS.nc'
CHARACTER(LEN=10), PARAMETER      :: YFILE_FORCING_OUT = 'FORCING.nc'
!
CHARACTER(LEN=6)                  :: YPROG = 'OFFLIN'
!
CHARACTER(LEN=100)                :: YCOMMENT = ' '
!
TYPE (DATE_TIME)                  :: TDTCUR
!
INTEGER                           :: IFILE_ID
INTEGER                           :: IVAR_ID
INTEGER                           :: IRES
INTEGER                           :: IRET
INTEGER,DIMENSION(2)              :: IDDIM,IDIMS
CHARACTER(LEN=100) ,DIMENSION(2)  :: YNAME_DIM
CHARACTER(LEN=100) ,DIMENSION(2)  :: YATT_TITLE,YATT
!
!
REAL, DIMENSION(:),    ALLOCATABLE :: ZZENITH, ZAZIM, ZTSUN
!
INTEGER                                 :: JT ! loop counter on times
REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: ZF ! field to write
CHARACTER(LEN=12) ::  YEXPER    ! experiment name
INTEGER :: ILUNAM, ILUOUT
CHARACTER(LEN=28)  :: YLUOUT    ='LISTING_FORCING'   ! name of the listing
LOGICAL :: GFOUND
REAL(KIND=JPRB)   :: ZHOOK_HANDLE
!
CHARACTER(LEN=12) :: YEXPERIMENT_NAME
INTEGER           :: NUMBER_GRID_CELLS
INTEGER           :: NUMBER_OF_TIME_STEPS_INPUT
INTEGER           :: NUMBER_OF_TIME_STEPS_FINAL
REAL              :: ZATM_FORC_STEP, ZDEN
CHARACTER(LEN=6)  :: YFORCING_FILETYPE       ! output file type:'ASCII ', 'BINARY', 'NETCDF'
NAMELIST/NAM_MY_PARAM/YEXPERIMENT_NAME,NUMBER_GRID_CELLS,NUMBER_OF_TIME_STEPS_INPUT, &
NUMBER_OF_TIME_STEPS_FINAL,ZATM_FORC_STEP,YFORCING_FILETYPE
!----------------------------------------------------------------------------
!
!*    1.     Initialization of parameter variables of the simulation
!            -------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PRE_INPUT_EXPERIMENT',0,ZHOOK_HANDLE)
CALL ALLOC_SURFEX(1)
CALL GOTO_SURFEX(1,.TRUE.)
!
CALL GET_LUOUT('ASCII ',ILUOUT)
OPEN(UNIT=ILUOUT,FILE=ADJUSTL(ADJUSTR(YLUOUT)//'.txt'),FORM='FORMATTED',ACTION='WRITE')
!
YFORCING_FILETYPE='NETCDF'
!
ILUOUT=0
CALL OPEN_NAMELIST('ASCII ',ILUNAM,'MY_PARAM.nam                ')
CALL POSNAM(ILUNAM,'NAM_MY_PARAM',GFOUND,ILUOUT)
IF (GFOUND) READ (UNIT=ILUNAM,NML=NAM_MY_PARAM)
CPROGNAME = 'ASCII '
!
YEXPER    = YEXPERIMENT_NAME
INI       = NUMBER_GRID_CELLS
INPTS     = NUMBER_OF_TIME_STEPS_INPUT
JNPTS     = NUMBER_OF_TIME_STEPS_FINAL
ZTSTEPFRC = ZATM_FORC_STEP
!
print*,' > ==================================================================='
print*,' > PRE_INPUT_EXPERIMENT: YEXPER             = ',YEXPER
print*,' > PRE_INPUT_EXPERIMENT: INI                = ',INI   
print*,' > PRE_INPUT_EXPERIMENT: INPTS              = ',INPTS 
print*,' > PRE_INPUT_EXPERIMENT: JNPTS              = ',JNPTS 
print*,' > PRE_INPUT_EXPERIMENT: ZTSTEPFRC          = ',ZTSTEPFRC 
print*,' > PRE_INPUT_EXPERIMENT: YFORCING_FILETYPE  = ',YFORCING_FILETYPE
print*,' > ==================================================================='
!
ALLOCATE(ZZREF(INI))
ALLOCATE(ZUREF(INI))
ALLOCATE(ZZS(INI))
ALLOCATE(ZLON(INI))
ALLOCATE(ZLAT(INI))
ALLOCATE(ZTA(INPTS,INI))
ALLOCATE(ZQA(INPTS,INI))
ALLOCATE(ZPS(INPTS,INI))
ALLOCATE(ZWINDSPEED(INPTS,INI))
ALLOCATE(ZWINDDIR(INPTS,INI))
ALLOCATE(ZDIR_SW(INPTS,INI))
ALLOCATE(ZSCA_SW(INPTS,INI))
ALLOCATE(ZLW(INPTS,INI))
ALLOCATE(ZRAIN(INPTS,INI))
ALLOCATE(ZSNOW(INPTS,INI))
ALLOCATE(ZCO2(INPTS,INI))
!
ALLOCATE(CSELECT(200))
CALL DEFAULT_DIAG_SURF_ATM(N2M, LSURF_BUDGET, L2M_MIN_ZS, LRAD_BUDGET, &
                           LCOEF, LSURF_VARS, LSURF_BUDGETC,           &
                           LRESET_BUDGETC, LSELECT, LPROVAR_TO_DIAG,   &
                           LDIAG_GRID, LFRAC, XDIAG_TSTEP, CSELECT     )

!
! 
!*    2.     Initialization of forcing variables
!            -----------------------------------
!
CALL MY_FORC(YEXPER,INI,INPTS,ZTSTEPFRC,      &
             TDTCUR%TDATE%YEAR,TDTCUR%TDATE%MONTH,      &
             TDTCUR%TDATE%DAY,TDTCUR%TIME,              &
             ZLON, ZLAT, ZZS, ZZREF, ZUREF,             &
             ZTA, ZQA, ZPS, ZWINDSPEED, ZWINDDIR,       &
             ZDIR_SW, ZSCA_SW, ZLW, ZRAIN, ZSNOW, ZCO2  )
!
!----------------------------------------------------------------------------
!      
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!---- DONT MODIFY BELOW HERE !!! --------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!     
!        3.5    solar time computation
!               ----------------------
CALL INI_CSTS
! 
ALLOCATE(ZTSUN(INI))
ALLOCATE(ZZENITH(INI))
ALLOCATE(ZAZIM(INI))
!
CALL SUNPOS ( TDTCUR%TDATE%YEAR, TDTCUR%TDATE%MONTH, TDTCUR%TDATE%DAY, TDTCUR%TIME, ZLON, ZLAT, ZTSUN, ZZENITH, ZAZIM)
!
!----------------------------------------------------------------------------
!
IF (YFORCING_FILETYPE == 'BINARY') THEN
!      
!*    4.     Writing in binary files
!            -----------------------
!
ALLOCATE(ZF(INI))
!
CALL OPEN_CLOSE_BIN_ASC_FORC('CONF ','BINARY',INI,'W')
CALL OPEN_CLOSE_BIN_ASC_FORC('OPEN ','BINARY',INI,'W')
!
!* configuration file is written for both BINARY and ASCII forcing file types
WRITE(21,*) 'N' !for the forcing swap: set Y to swap, N not to swap
WRITE(21,*) INI
WRITE(21,*) JNPTS
WRITE(21,*) ZTSTEPFRC
WRITE(21,*) TDTCUR%TDATE%YEAR
WRITE(21,*) TDTCUR%TDATE%MONTH
WRITE(21,*) TDTCUR%TDATE%DAY
WRITE(21,*) TDTCUR%TIME
WRITE(21,FMT='(50(F15.8))') ZLON
WRITE(21,FMT='(50(F15.8))') ZLAT
WRITE(21,FMT='(50(F15.8))') ZZS
WRITE(21,FMT='(50(F15.8))') ZZREF
WRITE(21,FMT='(50(F15.8))') ZUREF
!
DO JT=1,INPTS
  ZF = ZTA(JT,:)
  WRITE(22,REC=JT) ZF(:)
  ZF = ZQA(JT,:)
  WRITE(23,REC=JT) ZF(:)
  ZF = ZWINDSPEED(JT,:)
  WRITE(24,REC=JT) ZF(:)
  ZF = ZLW(JT,:)
  WRITE(25,REC=JT) ZF(:)
  ZF = ZDIR_SW(JT,:)
  WRITE(26,REC=JT) ZF(:)
  ZF = ZSCA_SW(JT,:)
  WRITE(27,REC=JT) ZF(:)
  ZF = ZRAIN(JT,:)
  WRITE(28,REC=JT) ZF(:)
  ZF = ZSNOW(JT,:)
  WRITE(29,REC=JT) ZF(:)
  ZF = ZPS(JT,:)
  WRITE(30,REC=JT) ZF(:)
  ZF = ZWINDDIR(JT,:)
  WRITE(31,REC=JT) ZF(:)
  ZF = ZCO2(JT,:)
  WRITE(32,REC=JT) ZF(:)
END DO
CALL OPEN_CLOSE_BIN_ASC_FORC('CLOSE','BINARY',INI,'W')
!
DEALLOCATE(ZF)
!
!----------------------------------------------------------------------------
!
ELSE IF (YFORCING_FILETYPE == 'ASCII ') THEN
!      
!*    5.     Writing in ASCII files
!            ----------------------
!
CALL OPEN_CLOSE_BIN_ASC_FORC('CONF ','ASCII ',INI,'W')
CALL OPEN_CLOSE_BIN_ASC_FORC('OPEN ','ASCII ',INI,'W')
!
!* configuration file is written for both BINARY and ASCII forcing file types
WRITE(21,*) INI
WRITE(21,*) JNPTS
WRITE(21,*) ZTSTEPFRC
WRITE(21,*) TDTCUR%TDATE%YEAR
WRITE(21,*) TDTCUR%TDATE%MONTH
WRITE(21,*) TDTCUR%TDATE%DAY
WRITE(21,*) TDTCUR%TIME
WRITE(21,FMT='(50(F15.8))') ZLON
WRITE(21,FMT='(50(F15.8))') ZLAT
WRITE(21,FMT='(50(F15.8))') ZZS
WRITE(21,FMT='(50(F15.8))') ZZREF
WRITE(21,FMT='(50(F15.8))') ZUREF
!
!* writes forcing in ASCII files
DO JT=1,INPTS
  WRITE(UNIT=22,FMT='(50(F20.5))') ZTA(JT,:)
  WRITE(UNIT=23,FMT='(50(F20.5))') ZQA(JT,:)
  WRITE(UNIT=24,FMT='(50(F20.5))') ZWINDSPEED(JT,:)
  WRITE(UNIT=25,FMT='(50(F20.5))') ZLW(JT,:)
  WRITE(UNIT=26,FMT='(50(F20.5))') ZDIR_SW(JT,:)
  WRITE(UNIT=27,FMT='(50(F20.5))') ZSCA_SW(JT,:)
  WRITE(UNIT=28,FMT='(50(F20.5))') ZRAIN(JT,:)
  WRITE(UNIT=29,FMT='(50(F20.5))') ZSNOW(JT,:)
  WRITE(UNIT=30,FMT='(50(F20.5))') ZPS(JT,:)
  WRITE(UNIT=31,FMT='(50(F20.5))') ZWINDDIR(JT,:)
  WRITE(UNIT=32,FMT='(50(F20.5))') ZCO2(JT,:)
END DO
CALL OPEN_CLOSE_BIN_ASC_FORC('CLOSE','ASCII ',INI,'W')
!
!----------------------------------------------------------------------------
!
ELSE IF (YFORCING_FILETYPE == 'NETCDF') THEN
!      
!*    4.     Writing of PARAMS.nc file
!            -------------------------
!
!----------------------------------------------------------------------------
!      
!        4.1    define dimensions
!               -----------------
!
IDIMS(1) = INI
!
!----------------------------------------------------------------------------
!      
!        4.2    define dimension names
!               ----------------------
!
YNAME_DIM(1) = 'Number_of_points'
!
!----------------------------------------------------------------------------
!      
!*       5.     Writing of FORCING.nc file
!               --------------------------
!
!----------------------------------------------------------------------------
!      
!        5.1    define dimensions
!               -----------------
!
IDIMS(1) = INI    ! space dimension
IDIMS(2) = JNPTS  ! time dimension   
!
!----------------------------------------------------------------------------
!      
!        5.2    define dimension names
!               ----------------------
!
YNAME_DIM(1) = 'Number_of_points'
YNAME_DIM(2) = 'time'
!
!----------------------------------------------------------------------------
!    
!        5.3    create file
!               -----------
!
CALL CREATE_FILE(YFILE_FORCING_OUT,IDIMS,YNAME_DIM(1:2),IFILE_ID,IDDIM)
!
!----------------------------------------------------------------------------
!
ALLOCATE(ZTIME(INPTS))
!
IF (ZTSTEPFRC == FLOOR(ZTSTEPFRC/86400.)*86400) THEN 
  ZDEN = 86400.
ELSEIF (ZTSTEPFRC == FLOOR(ZTSTEPFRC/3600.)*3600) THEN
  ZDEN = 3600.
ELSEIF (ZTSTEPFRC == FLOOR(ZTSTEPFRC/60.)*60) THEN
  ZDEN = 60.
ELSE
  ZDEN = 1.
ENDIF

DO JT = 1, INPTS
  ZTIME(JT) = ZTSTEPFRC/ZDEN * (JT-1)
ENDDO
!
YATT_TITLE(1) = 'units'
TDTCUR%TIME = TDTCUR%TIME + (ZTSTEPFRC-FLOOR(ZTSTEPFRC/86400)*86400)
CALL GET_DATE_OL(TDTCUR,ZTSTEPFRC,YATT(1))
!
CALL WRITE_NETCDF(IFILE_ID,'FORC_TIME_STEP','Forcing_Time_Step',ZTSTEPFRC)
!
CALL WRITE_NETCDF(IFILE_ID,'time','',ZTIME(:),IDDIM(2),YATT_TITLE(1:1),YATT(1:1))
!
CALL WRITE_NETCDF(IFILE_ID,'LON','Longitude', ZLON, IDDIM(1))
CALL WRITE_NETCDF(IFILE_ID,'LAT','Latitude', ZLAT, IDDIM(1))
!
CALL WRITE_NETCDF(IFILE_ID,'ZS','Surface_Orography',ZZS, IDDIM(1))
!
YATT_TITLE(1)='units'
YATT(1) = 'm'
CALL WRITE_NETCDF(IFILE_ID,'ZREF','Reference_Height',ZZREF,IDDIM(1),YATT_TITLE(1:1),YATT(1:1))
CALL WRITE_NETCDF(IFILE_ID,'UREF','Reference_Height_for_Wind',ZUREF,IDDIM(1),YATT_TITLE(1:1),YATT(1:1))
!
! 2D VARIABLES WITH 2 COMMENTS
!
YATT_TITLE(1) = 'measurement_height' 
YATT      (1) = '2m'
YATT_TITLE(2) = 'units'
!
YATT      (2) = 'K'
CALL WRITE_NETCDF(IFILE_ID,'Tair','Near_Surface_Air_Temperature',TRANSPOSE(ZTA(1:JNPTS,:)),&
        IDDIM(1),IDDIM(2),YATT_TITLE(1:2),YATT(1:2))
!
YATT      (2) = 'Kg/Kg'
CALL WRITE_NETCDF(IFILE_ID,'Qair','Near_Surface_Specific_Humidity',TRANSPOSE(ZQA(1:JNPTS,:)),&
        IDDIM(1),IDDIM(2),YATT_TITLE(1:2),YATT(1:2))
!
YATT      (2) = 'm/s'
CALL WRITE_NETCDF(IFILE_ID,'Wind','Wind_Speed',TRANSPOSE(ZWINDSPEED(1:JNPTS,:)),&
        IDDIM(1),IDDIM(2),YATT_TITLE(1:2),YATT(1:2))
!
!
! 2D VARIABLES WITH 1 COMMENT
!
YATT_TITLE(1) = 'units' 
!
YATT(1) = 'W/m2'
CALL WRITE_NETCDF(IFILE_ID,'DIR_SWdown','Surface_Indicent_Direct_Shortwave_Radiation' ,TRANSPOSE(ZDIR_SW(1:JNPTS,:)),&
        IDDIM(1),IDDIM(2),YATT_TITLE(1:1),YATT(1:1))
!
YATT(1) = 'W/m2'
CALL WRITE_NETCDF(IFILE_ID,'SCA_SWdown','Surface_Incident_Diffuse_Shortwave_Radiation' ,TRANSPOSE(ZSCA_SW(1:JNPTS,:)),&
        IDDIM(1),IDDIM(2),YATT_TITLE(1:1),YATT(1:1))
!
YATT(1) = 'W/m2'
CALL WRITE_NETCDF(IFILE_ID,'LWdown','Surface_Incident_Longwave_Radiation' ,TRANSPOSE(ZLW(1:JNPTS,:)),&
        IDDIM(1),IDDIM(2),YATT_TITLE(1:1),YATT(1:1))
!
YATT(1) = 'Pa'
CALL WRITE_NETCDF(IFILE_ID,'PSurf','Surface_Pressure',TRANSPOSE(ZPS(1:JNPTS,:)),&
        IDDIM(1),IDDIM(2),YATT_TITLE(1:1),YATT(1:1))
!
YATT(1) = 'Kg/m2/s'
CALL WRITE_NETCDF(IFILE_ID,'Rainf','Rainfall_Rate',TRANSPOSE(ZRAIN(1:JNPTS,:)),&
        IDDIM(1),IDDIM(2),YATT_TITLE(1:1),YATT(1:1))
!
YATT(1) = 'Kg/m2/s'
CALL WRITE_NETCDF(IFILE_ID,'Snowf','Snowfall_Rate',TRANSPOSE(ZSNOW(1:JNPTS,:)),&
        IDDIM(1),IDDIM(2),YATT_TITLE(1:1),YATT(1:1))
!
YATT(1) = 'Kg/m3'
CALL WRITE_NETCDF(IFILE_ID,'CO2air','Near_Surface_CO2_Concentration',TRANSPOSE(ZCO2(1:JNPTS,:)),&
        IDDIM(1),IDDIM(2),YATT_TITLE(1:1),YATT(1:1))
!
YATT(1) = 'deg'
CALL WRITE_NETCDF(IFILE_ID,'Wind_DIR','Wind_Direction',TRANSPOSE(ZWINDDIR(1:JNPTS,:)),&
        IDDIM(1),IDDIM(2),YATT_TITLE(1:1),YATT(1:1))
!
!              2.4 closing file
!     ------------
!
IRET=NF_CLOSE(IFILE_ID)

ELSE
        PRINT*,' ABORT: CHECK YFORCING_FILETYPE '
        IF (LHOOK) CALL DR_HOOK('PRE_INPUT_EXPERIMENT',1,ZHOOK_HANDLE)
        STOP
ENDIF
!
CLOSE(ILUOUT)
CALL DEALLOC_SURFEX
IF (LHOOK) CALL DR_HOOK('PRE_INPUT_EXPERIMENT',1,ZHOOK_HANDLE)

END
