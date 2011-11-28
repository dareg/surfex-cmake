!     #########
       SUBROUTINE INIT_OUTFN_SURF_ATM_n(HPROGRAM,KLUOUT)
!     ###############################
!
!
!!****  *INIT_OUTFN_SURF_ATM_n* -  create output files and defines variables
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
!!      Original    05-04
!!      Modified 06-10 by S. Faroux
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------

USE MODD_DIAG_SURF_ATM_n, ONLY : LFRAC, LCOEF, N2M, LSURF_BUDGET, LSURF_BUDGETC, LSELECT
USE MODD_OL_FILEID,       ONLY : XVAR_TO_FILEOUT, XID, XOUT
USE MODN_IO_OFFLINE,      ONLY : XTSTEP_OUTPUT
USE MODD_SURF_ATM_n,      ONLY : TTIME
!
USE MODI_GET_DIM_FOR_2D
USE MODI_GET_COORD_FOR_2D
USE MODI_GET_DATE_OL
USE MODI_FILL_ID_OL
USE MODI_CREATE_FILE
USE MODI_DEF_VAR_NETCDF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_GET_DIM_FULL_n
!
IMPLICIT NONE
include 'netcdf.inc'

CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM
INTEGER, INTENT(IN) :: KLUOUT
!
!*       0.2   Declarations of local variables
!              -------------------------------
INTEGER, DIMENSION(:), ALLOCATABLE  :: IDDIM, IDIMS
INTEGER, DIMENSION(:), ALLOCATABLE  :: IDDIM1, IDIMS1
CHARACTER(LEN=100), DIMENSION(:),ALLOCATABLE :: YNAME_DIM1
CHARACTER(LEN=100), DIMENSION(:),ALLOCATABLE :: YNAME_DIM
CHARACTER(LEN=100), DIMENSION(1) :: YATT_TITLE, YATT
CHARACTER(LEN=50)                :: YFILE
INTEGER                          :: IRET, INB_FORC, INI
INTEGER                          :: IFILE_ID, IVAR_ID
CHARACTER(LEN=40),DIMENSION(1)   :: HDATE
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
IF (LHOOK) CALL DR_HOOK('INIT_OUTFN_SURF_ATM_N',0,ZHOOK_HANDLE)
CALL GET_DIM_FULL_n(INI)

CALL GET_DIM_FOR_2D(IDIM1,IDIM2,HTYPE)
IF (IDIM1.NE.0) THEN
        ALLOCATE(ZX(IDIM1))
        ALLOCATE(ZY(IDIM2))
        ALLOCATE(IDIMS(3))
        ALLOCATE(IDDIM(3))
        ALLOCATE(YNAME_DIM(3))
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
        ALLOCATE(IDIMS(2))
        ALLOCATE(IDDIM(2))
        ALLOCATE(YNAME_DIM(2))
        IDIMS(1)=INI
        YNAME_DIM(1)='Number_of_points'
ENDIF
I0=SIZE(IDDIM)
IDIMS(I0)=NF_UNLIMITED
YNAME_DIM(I0)='time'
ALLOCATE(IDIMS1(I0-1))
ALLOCATE(IDDIM1(I0-1))
ALLOCATE(YNAME_DIM1(I0-1))

! 2. define dimension length
!---------------------------
IDIMS1=IDIMS(1:I0-1)

! 3. define dimension name
!---------------------------
YNAME_DIM1=YNAME_DIM(1:I0-1)

CALL GET_DATE_OL(TTIME,XTSTEP_OUTPUT,HDATE(1))

! 4. Create output file for fraction of tiles
!--------------------------------------------

IF (ALLOCATED(XVAR_TO_FILEOUT)) DEALLOCATE(XVAR_TO_FILEOUT)
IF (ALLOCATED(XID)) DEALLOCATE(XID)
ALLOCATE(XVAR_TO_FILEOUT(0))
ALLOCATE(XID(0))
XOUT=0

YFILE='SURF_ATM.OUT.nc'
CALL CREATE_FILE(YFILE,IDIMS1,YNAME_DIM1,IFILE_ID,IDDIM1)
!PRINT*,'  CREATING : ',YFILE

IF (IDIM1.NE.0) THEN
      CALL DEF_VAR_NETCDF(IFILE_ID,TRIM(YNAME_DIM(1)),'',IDDIM(1:1),YATT_TITLE,YUNIT1)
      CALL DEF_VAR_NETCDF(IFILE_ID,TRIM(YNAME_DIM(2)),'',IDDIM(2:2),YATT_TITLE,YUNIT2)
ENDIF
IF (LFRAC) THEN
   YATT_TITLE(1)='units'
   YATT(1)='%'
   CALL DEF_VAR_NETCDF(IFILE_ID,'FRAC_SEA   ','Fraction_of_sea   ',IDDIM1,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'FRAC_WATER ','Fraction_of_water ',IDDIM1,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'FRAC_TOWN  ','Fraction_of_town  ',IDDIM1,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'FRAC_NATURE','Fraction_of_nature',IDDIM1,YATT_TITLE,YATT)
ENDIF
!
JRET=NF_OPEN(YFILE,NF_WRITE,IFILE_ID)
CALL FILL_ID_OL(IFILE_ID)
IF (IDIM1.NE.0) THEN
  JRET=NF_INQ_VARID   (IFILE_ID,YNAME_DIM(1),IVAR_ID)
  JRET=NF_PUT_VAR_DOUBLE(IFILE_ID,IVAR_ID,ZX)
  JRET=NF_INQ_VARID   (IFILE_ID,YNAME_DIM(2),IVAR_ID)
  JRET=NF_PUT_VAR_DOUBLE(IFILE_ID,IVAR_ID,ZY)
ENDIF

! 5. Create output file for diagnostic variables
!-----------------------------------------------
YFILE='SURF_ATM_DIAGNOSTICS.OUT.nc'
CALL CREATE_FILE(YFILE,IDIMS,YNAME_DIM,IFILE_ID,IDDIM)
!PRINT*,'  CREATING : ',YFILE

YATT_TITLE='units'

IF (IDIM1.NE.0) THEN
      CALL DEF_VAR_NETCDF(IFILE_ID,TRIM(YNAME_DIM(1)),'',IDDIM(1:1),YATT_TITLE,YUNIT1)
      CALL DEF_VAR_NETCDF(IFILE_ID,TRIM(YNAME_DIM(2)),'',IDDIM(2:2),YATT_TITLE,YUNIT2)
ENDIF
CALL DEF_VAR_NETCDF(IFILE_ID,'time','',IDDIM(I0:I0),YATT_TITLE,HDATE)

IF (N2M>0) THEN
   YATT='(-)'
   CALL DEF_VAR_NETCDF(IFILE_ID,'RI'   ,'Averaged_Richardson_Number      '   ,IDDIM,YATT_TITLE,YATT)
ENDIF
!
IF (LCOEF) THEN
   YATT='W/s2'
   CALL DEF_VAR_NETCDF(IFILE_ID,'CD'   ,'Drag_Coefficient_For_Momentum   '   ,IDDIM,YATT_TITLE,YATT)
   YATT='W/s'
   CALL DEF_VAR_NETCDF(IFILE_ID,'CH'   ,'Drag_Coefficient_For_Heat       '   ,IDDIM,YATT_TITLE,YATT)
   YATT='W/s/K'
   CALL DEF_VAR_NETCDF(IFILE_ID,'CE'   ,'Drag_Coefficient_For_Evaporation'   ,IDDIM,YATT_TITLE,YATT)
   YATT='m'
   CALL DEF_VAR_NETCDF(IFILE_ID,'Z0'   ,'Roughness_Length_For_Momentum'   ,IDDIM,YATT_TITLE,YATT)
   YATT='m'
   CALL DEF_VAR_NETCDF(IFILE_ID,'Z0H'  ,'Roughness_Length_For_Heat'   ,IDDIM,YATT_TITLE,YATT)
   YATT='m'
   CALL DEF_VAR_NETCDF(IFILE_ID,'UREF'  ,'Reference_Height_For_Momentum'   ,IDDIM,YATT_TITLE,YATT)
   YATT='m'
   CALL DEF_VAR_NETCDF(IFILE_ID,'ZREF'  ,'Reference_Height_For_Heat'   ,IDDIM,YATT_TITLE,YATT)
ENDIF
!
IF (LSURF_BUDGET)  THEN
   YATT='W/m2'
   CALL DEF_VAR_NETCDF(IFILE_ID,'RN'   ,'Averaged_Net_Radiation     '   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'H'    ,'Averaged_Sensible_Heat_Flux'   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'LE'   ,'Averaged_Total_Latent_Heat_Flux  ',IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'LEI'  ,'Averaged_Sublimation_Latent_Heat_Flux  ',IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'GFLUX','Averaged_Ground_Heat_Flux  '   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'SWD'  ,'Averaged_Downward_SW       '   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'SWU'  ,'Averaged_Upward_SW         '   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'LWD'  ,'Averaged_Downward_LW       '   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'LWU'  ,'Averaged_Upward_LW         '   ,IDDIM,YATT_TITLE,YATT)
   YATT='kg/ms2'
   CALL DEF_VAR_NETCDF(IFILE_ID,'FMU'  ,'Averaged_Zonal_Wind_Stress '   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'FMV'  ,'Averaged_Merid_Wind_Stress '   ,IDDIM,YATT_TITLE,YATT)
ENDIF
!
IF (LSURF_BUDGETC) THEN
   YATT='J/m2'
   CALL DEF_VAR_NETCDF(IFILE_ID,'RNC'   ,'Cumulated_Averaged_Net_Radiation     '   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'HC'    ,'Cumulated_Averaged_Sensible_Heat_Flux'   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'LEC'   ,'Cumulated_Averaged_Total_Latent_Heat_Flux  ',IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'LEIC'  ,'Cumulated_Averaged_Sublimation_Latent_Heat_Flux  ',IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'GFLUXC','Cumulated_Averaged_Ground_Heat_Flux  '   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'SWDC'  ,'Cumulated_Averaged_Downward_SW       '   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'SWUC'  ,'Cumulated_Averaged_Upward_SW         '   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'LWDC'  ,'Cumulated_Averaged_Downward_LW       '   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'LWUC'  ,'Cumulated_Averaged_Upward_LW         '   ,IDDIM,YATT_TITLE,YATT)
   YATT='kg/ms'
   CALL DEF_VAR_NETCDF(IFILE_ID,'FMUC'  ,'Cumulated_Averaged_Zonal_Wind_Stress '   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'FMVC'  ,'Cumulated_Averaged_Merid_Wind_Stress '   ,IDDIM,YATT_TITLE,YATT)
ENDIF
!
IF (N2M>=1.OR.LSURF_BUDGET.OR.LSURF_BUDGETC) THEN
   YATT='K'
   CALL DEF_VAR_NETCDF(IFILE_ID,'TS' ,'Surface_Temperature         '   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'TSRAD' ,'Surface_Skin,Temperature '   ,IDDIM,YATT_TITLE,YATT)
   YATT='-'
   CALL DEF_VAR_NETCDF(IFILE_ID,'EMIS' ,'Surface_emissivity '   ,IDDIM,YATT_TITLE,YATT)
ENDIF
!
IF (N2M>0) THEN
   YATT='K'
   CALL DEF_VAR_NETCDF(IFILE_ID,'T2M_MIN_ZS' ,'2m_Temperature         '   ,IDDIM,YATT_TITLE,YATT)
   YATT='kg/kg'
   CALL DEF_VAR_NETCDF(IFILE_ID,'Q2M_MIN_ZS' ,'2m_Specific_Humidity   '   ,IDDIM,YATT_TITLE,YATT)
   YATT='(-)'
   CALL DEF_VAR_NETCDF(IFILE_ID,'HU2M_MIN_ZS','2m_Relative_Humidity   '   ,IDDIM,YATT_TITLE,YATT)
   YATT='K'
   CALL DEF_VAR_NETCDF(IFILE_ID,'T2M' ,'2m_Temperature         '   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'T2MMIN' ,'Minimum_2m_Temperature         '   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'T2MMAX' ,'Maximum_2m_Temperature         '   ,IDDIM,YATT_TITLE,YATT)
   YATT='kg/kg'
   CALL DEF_VAR_NETCDF(IFILE_ID,'Q2M' ,'2m_Specific_Humidity   '   ,IDDIM,YATT_TITLE,YATT)
   YATT='(-)'
   CALL DEF_VAR_NETCDF(IFILE_ID,'HU2M','2m_Relative_Humidity   '   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'HU2MMIN','Minimum_2m_Relative_Humidity   '   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'HU2MMAX','Maximum_2m_Relative_Humidity   '   ,IDDIM,YATT_TITLE,YATT)
   YATT='m/s'
   CALL DEF_VAR_NETCDF(IFILE_ID,'ZON10M','10m_Zonal_wind       '   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'MER10M','10m_Meridian_Wind     '   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'W10M','10m_Wind     '   ,IDDIM,YATT_TITLE,YATT)
   CALL DEF_VAR_NETCDF(IFILE_ID,'W10MMAX','Maximum_10m_Wind     ' ,IDDIM,YATT_TITLE,YATT)
ENDIF
!
JRET=NF_OPEN(YFILE,NF_WRITE,IFILE_ID)
CALL FILL_ID_OL(IFILE_ID)
IF (IDIM1.NE.0) THEN
  JRET=NF_INQ_VARID   (IFILE_ID,YNAME_DIM(1),IVAR_ID)
  JRET=NF_PUT_VAR_DOUBLE(IFILE_ID,IVAR_ID,ZX)
  JRET=NF_INQ_VARID   (IFILE_ID,YNAME_DIM(2),IVAR_ID)
  JRET=NF_PUT_VAR_DOUBLE(IFILE_ID,IVAR_ID,ZY)
ENDIF
IF (LHOOK) CALL DR_HOOK('INIT_OUTFN_SURF_ATM_N',1,ZHOOK_HANDLE)

END SUBROUTINE INIT_OUTFN_SURF_ATM_n
