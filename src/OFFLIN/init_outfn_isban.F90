!     #########
       SUBROUTINE INIT_OUTFN_ISBA_n (CHI, DE, DGO, DMI, GB, SB, IO, & 
                                     TPSNOW, TPTIME, UG, U, HSELECT,    &
                                     OPROVAR_TO_DIAG, HPROGRAM, KLUOUT)
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
!!      modified    11-03,by P. Le Moigne   *Meteo France*
!!      modified    05-04,by P. Le Moigne : surf_atm diagnostics moved at the
!!                                           right place
!!      modified    10-04,by P. Le Moigne : add new diagnostics
!!      modified    10-04,by P. Le Moigne : add Halstead coefficient
!!      modified     2008,by B. Decharme  : limit the number of diag
!!                                           Add floodplains diag
!!      modified    04-09,by A.L. Gibelin : Add respiration diagnostics
!!      modified    05-09,by A.L. Gibelin : Add carbon spinup
!!      modified    07-09,by A.L. Gibelin : Add carbon prognostic variables
!!  
!!      modified    09-12,by B. Decharme  : delete LPROVAR_TO_DIAG for prognostic variables
!!                                           delete NWG_LAYER
!!                                           Erroneous description in diag comments
!!      modified    06-13,by B. Decharme  : good dimension for Tg,Wg,et Wgi
!!                                           bug : TSN_VEG if Snowlayer = 1 ; 
!!                                           bug : TSRAD_P and not TTSRAD_P
!!                                           add diag (Qsb,Subl) and Snow noted SN
!!      modified    10-14,by P. Samuelsson: Added MEB output
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_TYPE_SNOW
USE MODD_TYPE_DATE_SURF
!
!
USE MODD_CH_ISBA_n,ONLY : CH_ISBA_t
USE MODD_DIAG_EVAP_ISBA_n,ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_n,ONLY : DIAG_OPTIONS_t
USE MODD_DIAG_MISC_ISBA_n,ONLY : DIAG_MISC_ISBA_t
USE MODD_GR_BIOG_n,ONLY : GR_BIOG_t
USE MODD_CANOPY_n,ONLY : CANOPY_t
USE MODD_ISBA_OPTIONS_n,ONLY : ISBA_OPTIONS_t
USE MODD_SURF_ATM_GRID_n,ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n,ONLY : SURF_ATM_t
!
USE MODD_SURF_PAR,        ONLY : NUNDEF
USE MODD_DATA_COVER_PAR,  ONLY : NVEGTYPE
USE MODD_OL_FILEID,       ONLY : XVAR_TO_FILEOUT,XID,XOUT
USE MODD_ASSIM,          ONLY : LASSIM, CASSIM, CASSIM_ISBA, NVAR
USE MODD_AGRI,          ONLY : LAGRIP
!
!
USE MODN_IO_OFFLINE,      ONLY : XTSTEP_OUTPUT
!
USE MODI_GET_DIM_FULL_n
USE MODI_GET_ISBA_CONF_n
USE MODI_OL_DEFINE_DIM
USE MODI_GET_DATE_OL
USE MODI_CREATE_FILE
USE MODI_DEF_VAR_NETCDF
USE MODI_OL_WRITE_COORD
!
USE YOMHOOK ,ONLY : LHOOK,  DR_HOOK
USE PARKIND1,ONLY : JPRB
!
IMPLICIT NONE
include'netcdf.inc'

!
TYPE(CH_ISBA_t),INTENT(INOUT) :: CHI
TYPE(DIAG_EVAP_ISBA_t),INTENT(INOUT) :: DE
TYPE(DIAG_OPTIONS_t),INTENT(INOUT) :: DGO
TYPE(DIAG_MISC_ISBA_t),INTENT(INOUT) :: DMI
TYPE(GR_BIOG_t),INTENT(INOUT) :: GB
TYPE(CANOPY_t),INTENT(INOUT) :: SB
TYPE(ISBA_OPTIONS_t),INTENT(INOUT) :: IO
TYPE(SURF_SNOW),INTENT(INOUT) :: TPSNOW
TYPE(DATE_TIME),INTENT(INOUT) :: TPTIME
TYPE(SURF_ATM_GRID_t),INTENT(INOUT) :: UG
TYPE(SURF_ATM_t),INTENT(INOUT) :: U
!
 CHARACTER(LEN=*),DIMENSION(:),INTENT(IN) :: HSELECT
 LOGICAL,INTENT(IN) :: OPROVAR_TO_DIAG
!
 CHARACTER(LEN=6),INTENT(IN) :: HPROGRAM
INTEGER,         INTENT(IN) :: KLUOUT
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
 CHARACTER(LEN=100),DIMENSION(:),POINTER :: YNAME_DIM
 CHARACTER(LEN=100),DIMENSION(1) :: YATT_TITLE,YATT
 CHARACTER(LEN=40),DIMENSION(1)   :: YDATE
 CHARACTER(LEN=13),DIMENSION(1)   :: YUNIT1,YUNIT2
 CHARACTER(LEN=100)               :: YCOMMENT  
 CHARACTER(LEN=50)                :: YFILE
 CHARACTER(LEN=12)                :: YRECFM
 CHARACTER(LEN=3)                 :: YPAS,YPAT
 CHARACTER(LEN=6) :: YLVL
 CHARACTER(LEN=3)                 :: YISBA
 CHARACTER(LEN=2)                 :: YLVLV
 CHARACTER(LEN=1) :: YVAR 
! 
REAL,DIMENSION(:),POINTER       :: ZX,ZY
!
INTEGER,DIMENSION(:),POINTER   :: IDIMS,JDIM
INTEGER                          :: INI,INPATCH,INLVLD,INLVLS,INBIOMASS,&
                                    INLITTER,INLITTLEVS,INSOILCARB,JP
INTEGER                          :: JL,JVEG,JNBIOMASS,JNLITTER,JNLITTLEVS,JNSOILCARB
INTEGER                          :: IDIM1,INDIMS
INTEGER                          :: IFILE_ID,IDIMID,JSV
INTEGER                          :: IL,JRET, JVAR
INTEGER                          :: ISIZE_LMEB_PATCH   ! Number of patches where multi-energy balance should be applied
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------

! 1. Compute output lenght dimension
!-----------------------------------

IF (LHOOK) CALL DR_HOOK('INIT_OUTFN_ISBA_N',0,ZHOOK_HANDLE)
!
ISIZE_LMEB_PATCH=COUNT(IO%LMEB_PATCH(:))
!
 CALL GET_DIM_FULL_n(U%NDIM_FULL,INI)
 CALL GET_ISBA_CONF_n(IO,TPSNOW%NLAYER,YISBA,INPATCH,INLVLD,INLVLS,INBIOMASS,&
                       INLITTER,INLITTLEVS,INSOILCARB)  
!
 CALL OL_DEFINE_DIM(UG,U%NSIZE_FULL,HPROGRAM,KLUOUT,INI,IDIM1,YUNIT1,YUNIT2,&
                   ZX,ZY,IDIMS,JDIM,YNAME_DIM)
 CALL GET_DATE_OL(TPTIME,XTSTEP_OUTPUT,YDATE(1))
!
INDIMS = SIZE(JDIM)
!
! 4. Create output file for prognostic variables
!----------------------------------------------------------
!
IF (ALLOCATED(XVAR_TO_FILEOUT)) DEALLOCATE(XVAR_TO_FILEOUT)
IF (ALLOCATED(XID)) DEALLOCATE(XID)
ALLOCATE(XID(0))
XOUT=0
!
!
YFILE='ISBA_PROGNOSTIC.OUT.nc'
 CALL CREATE_FILE(YFILE,IDIMS,YNAME_DIM,IFILE_ID,JDIM)
JRET=NF_REDEF(IFILE_ID) 
!
!
DO JP = 1,IO%NPATCH
!
IF (JP>=10) THEN
  WRITE(YPAT,'(I2.2)') JP
ELSE
  WRITE(YPAT,'(I1.1)') JP
ENDIF
YPAT = "P"//ADJUSTL(YPAT(:LEN_TRIM(YPAT)))
!
IF(IO%LTEMP_ARP)THEN
  IL=IO%NTEMPLAYER_ARP
ELSEIF(IO%CISBA=='DIF')THEN
  IL=INLVLD
ELSE
  IL=2
ENDIF
!
YATT_TITLE(1)='units'
!
DO JL=1,IL
  WRITE(YPAS,'(I3)') JL ; YLVL = TRIM(ADJUSTL(YPAS))//YPAT
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'TG'//YLVL,'Soil_temp_layer_'//YLVL,JDIM,YATT_TITLE,(/'Kelvin'/))
ENDDO
!
IL=INLVLD
!
DO JL=1,IL
  WRITE(YPAS,'(I3)') JL ; YLVL=TRIM(ADJUSTL(YPAS))//YPAT
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WG'//YLVL,'Soil_liquid_layer_'//YLVL,JDIM,YATT_TITLE,(/'m3/m3'/))
ENDDO
!  
IF(IO%CISBA/='DIF')THEN
   IL=2
ENDIF
DO JL=1,IL
  WRITE(YPAS,'(I3)') JL ; YLVL=TRIM(ADJUSTL(YPAS))//YPAT
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WGI'//YLVL,'Soil_ice_layer_'//YLVL,JDIM,YATT_TITLE,(/'m3/m3'/))
ENDDO  
!
 CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WR'//YPAT,'Interception_reservoir'//YPAT,JDIM,YATT_TITLE,(/'mm'/))
 CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'RESA'//YPAT,'Aerodynamic_resistance'//YPAT,JDIM,YATT_TITLE,(/'s/m'/))
!
IF (ISIZE_LMEB_PATCH>0) THEN
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WRV'//YPAT,'MEB: water intercepted on canopy vegetation leaves'//YPAT,&
          JDIM,YATT_TITLE,(/'mm'/))
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WRVN'//YPAT,'MEB: snow intercepted on canopy vegetation leaves'//YPAT,&
          JDIM,YATT_TITLE,(/'mm'/))
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'TV'//YPAT,'MEB: canopy vegetation temperature'//YPAT,JDIM,YATT_TITLE,(/'K'/))
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'TC'//YPAT,'MEB: vegetation canopy air temperature'//YPAT,JDIM,YATT_TITLE,(/'K'/))
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'QC'//YPAT,'MEB: vegetation canopy specifc humidity'//YPAT,JDIM,YATT_TITLE,(/'kg/kg'/))
ENDIF
!
IF(IO%LGLACIER)THEN
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'ICE_STO'//YPAT,'Glacier_reservoir'//YPAT,JDIM,YATT_TITLE,(/'Kg/m2'/))
ENDIF
!
DO JL=1,INLVLS
  WRITE(YPAS,'(I3)') JL
  YLVL = TRIM(ADJUSTL(YPAS))//YPAT
  !
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WSN_VEG'//YLVL,'Snow_Water_Equivalent_layer_'//YLVL,JDIM,YATT_TITLE,(/'Kg/m2'/))
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'RSN_VEG'//YLVL,'Snow_density_layer_'//YLVL,JDIM,YATT_TITLE,(/'Kg/m3'/))
  IF (TPSNOW%SCHEME=='3-L' .OR. TPSNOW%SCHEME=='CRO') THEN   
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'HSN_VEG'//YLVL,'Snow_heat_layer'//YLVL,JDIM,YATT_TITLE,(/'J/m2'/))
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SAG_VEG'//YLVL,'Snow_age_param_layer_'//YLVL,JDIM,YATT_TITLE,(/'days_since_snowfall'/))
  ELSEIF(TPSNOW%SCHEME=='1-L') THEN
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'TSN_VEG'//YLVL,'Snow_temp_layer'//YLVL,JDIM,YATT_TITLE,(/'K'/))
  ENDIF
  IF (TPSNOW%SCHEME=='CRO') THEN   
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SG1_VEG'//YLVL,'Snow_grain_parameter1_layer_'//YLVL,JDIM,YATT_TITLE,(/'-'/))
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SG2_VEG'//YLVL,'Snow_grain_parameter2_layer_'//YLVL,JDIM,YATT_TITLE,(/'-'/))
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SHI_VEG'//YLVL,'Snow_historical_param_layer_'//YLVL,JDIM,YATT_TITLE,(/'-'/))
  ENDIF
ENDDO
!
 CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'ASN_VEG'//YPAT,'Snow_albedo'//YPAT,JDIM,YATT_TITLE,(/'-'/))
!
IF (IO%CPHOTO /='NON') THEN
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'AN'//YPAT,'Net CO2 Assimilation'//YPAT,JDIM,YATT_TITLE,(/'kgCO2/kgair m/s'/))
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'ANFM'//YPAT,'Leaf CO2 Assimilation'//YPAT,JDIM,YATT_TITLE,(/'kgCO2/kgair m/s'/))
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'ANDAY'//YPAT,'Daily Net CO2 Assimilation'//YPAT,JDIM,YATT_TITLE,(/'kgCO2/m2/day'/))
ENDIF
!
IF (IO%CPHOTO =='NIT' .OR. IO%CPHOTO =='NCB') THEN
  DO JNBIOMASS=1,INBIOMASS
    WRITE(YPAS,'(I3)') JNBIOMASS
    YLVL = TRIM(ADJUSTL(YPAS))//YPAT
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'BIOMA'//YLVL,'Plant biomass'//YLVL,JDIM,YATT_TITLE,(/'kgDM/m2'/))
  END DO
ENDIF
!
IF (IO%CRESPSL=='CNT') THEN
  DO JNLITTER=1,INLITTER
    DO JNLITTLEVS=1,INLITTLEVS
      WRITE(YPAS,'(I1,A1,I1)') JNLITTER,'_',JNLITTLEVS
      YLVL = TRIM(ADJUSTL(YPAS))//YPAT
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LITTER'//YLVL,'Litter pool'//YLVL,JDIM,YATT_TITLE,(/'gC/m2'/))
    END DO
  END DO  
  DO JNSOILCARB=1,INSOILCARB
    WRITE(YPAS,'(I3)') JNSOILCARB
    YLVL=TRIM(ADJUSTL(YPAS))//YPAT
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SOILCARB'//YLVL,'Soil carbon pool'//YLVL,JDIM,YATT_TITLE,(/'gC/m2'/))
  END DO
  DO JNLITTLEVS=1,INLITTLEVS
    WRITE(YPAS,'(I3)') JNLITTLEVS
    YLVL=TRIM(ADJUSTL(YPAS))//YPAT
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LIGNIN_STR'//YLVL,&
            'Ratio Lignin/Carbon in structural litter'//YLVL,JDIM,YATT_TITLE,(/'gC/m2'/))
  END DO
ENDIF
!
ENDDO
!
IF (IO%LCANOPY) THEN
  DO JL=1,SB%NLVL
    WRITE(YLVLV,'(i2.2)') JL
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'ISBA_CAN_Z'//YLVLV,'Canopy height',JDIM,YATT_TITLE,(/'m'/))
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'ISBA_CAN_U'//YLVLV,'Canopy wind',JDIM,YATT_TITLE,(/'m/s'/))
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'ISBA_CAN_T'//YLVLV,'Canopy temp',JDIM,YATT_TITLE,(/'K'/))
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'ISBA_CAN_Q'//YLVLV,'Canopy humidity',JDIM,YATT_TITLE,(/'kg/m3'/))
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'ISBA_CAN_E'//YLVLV,'Canopy TKE',JDIM,YATT_TITLE,(/'m2/s2'/))
  END DO
ENDIF
!
 CALL OL_WRITE_COORD(HSELECT,YFILE,IFILE_ID,JDIM,YATT_TITLE,YNAME_DIM,YUNIT1,YUNIT2,IDIM1,YDATE,ZX,ZY)
!
!
! 4. Create output file for fluxes values
!----------------------------------------------------------
!
YFILE='ISBA_DIAGNOSTICS.OUT.nc'
 CALL CREATE_FILE(YFILE,IDIMS,YNAME_DIM,IFILE_ID,JDIM)
JRET=NF_REDEF(IFILE_ID) 
YATT ='dimensionless'
!
IF (DGO%LCOEF) THEN
  YATT ='W/s2'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'CD_ISBA','Drag_Coefficient_For_Momentum',JDIM,YATT_TITLE,YATT)
  YATT ='W/s'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'CH_ISBA','Drag_Coefficient_For_Heat',JDIM,YATT_TITLE,YATT)
  YATT ='W/s/K'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'CE_ISBA','Drag_Coefficient_For_Evaporation',JDIM,YATT_TITLE,YATT)
  YATT ='m'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'Z0_ISBA','Roughness_Length_For_Momentum',JDIM,YATT_TITLE,YATT)
  YATT ='m'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'Z0H_ISBA','Roughness_Length_For_Heat',JDIM,YATT_TITLE,YATT)
ENDIF
!
IF (DGO%LSURF_VARS) THEN
  YATT ='kg/kg'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'QS_ISBA','Surface_Humidity',JDIM,YATT_TITLE,YATT)
ENDIF
!
IF (DGO%N2M>0) THEN
  !
  YATT ='dimensionless'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'RI_ISBA','Averaged_Richardson_Number',JDIM,YATT_TITLE,YATT) 
  YATT ='K'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'T2M_ISBA','2m_Temperature',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'T2MMIN_ISBA','Minimum_2m_Temperature',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'T2MMAX_ISBA','Maximum_2m_Temperature',JDIM,YATT_TITLE,YATT)
  YATT ='kg/kg'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'Q2M_ISBA','2m_Specific_Humidity',JDIM,YATT_TITLE,YATT)
  YATT ='(-)'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'HU2M_ISBA','2m_Relative_Humidity',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'HU2MMIN_ISBA','Minimum_2m_Relative_Humidity',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'HU2MMAX_ISBA','Maximum_2m_Relative_Humidity',JDIM,YATT_TITLE,YATT)
  YATT ='m/s'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'ZON10M_ISBA','10m_Zonal_wind',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'MER10M_ISBA','10m_Meridian_Wind',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'W10M_ISBA','10m_Wind',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'W10MMAX_ISBA','Maximum_10m_Wind',JDIM,YATT_TITLE,YATT)
  !
  IF(DGO%LPATCH_BUDGET) THEN
    DO JP = 1,IO%NPATCH
      IF (JP>=10) THEN
        WRITE(YPAT,'(I2.2)') JP
      ELSE
        WRITE(YPAT,'(I1.1)') JP
      ENDIF
      YPAT = "P"//ADJUSTL(YPAT(:LEN_TRIM(YPAT)))

      YATT='K'
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'T2M_'//YPAT,'2m_Temperature'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'T2MMIN_'//YPAT,'Minimum_2m_Temperature'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'T2MMAX_'//YPAT,'Maximum_2m_Temperature'//YPAT,JDIM,YATT_TITLE,YATT)
      YATT='kg/kg'
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'Q2M_'//YPAT,'2m_Specific_Humidity'//YPAT,JDIM,YATT_TITLE,YATT)
      YATT='(-)'
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'HU2M_'//YPAT,'2m_Relative_Humidity'//YPAT,JDIM,YATT_TITLE,YATT)
      YATT='m/s'
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'ZON10M_'//YPAT,'10m_Zonal_wind'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'MER10M_'//YPAT,'10m_Meridian_Wind'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'W10M_'//YPAT,'10m_Wind'//YPAT,JDIM,YATT_TITLE,YATT)
    ENDDO
  ENDIF
  !
ENDIF
!
IF (DGO%LSURF_BUDGET)  THEN
  !
  YATT ='W/m2'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'RN_ISBA','Averaged_Net_Radiation',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'H_ISBA','Averaged_Sensible_Heat_Flux',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LE_ISBA','Averaged_Total_Latent_Heat_Flux',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEI_ISBA','Averaged_Sublimation_Latent_Heat_Flux',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'GFLUX_ISBA','Averaged_Ground_Heat_Flux',JDIM,YATT_TITLE,YATT)
  !
  IF(DGO%LRAD_BUDGET)THEN
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SWD_ISBA','Averaged_Downward_SW',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SWU_ISBA','Averaged_Upward_SW',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LWD_ISBA','Averaged_Downward_LW',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LWU_ISBA','Averaged_Upward_LW',JDIM,YATT_TITLE,YATT)
  ENDIF
  !
  YATT ='Pa'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FMU_ISBA','Averaged_Zonal_Wind_Stress',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FMV_ISBA','Averaged_Merid_Wind_Stress',JDIM,YATT_TITLE,YATT)
  !
  IF (DGO%LPATCH_BUDGET) THEN
    !
    DO JP = 1,IO%NPATCH
      IF (JP>=10) THEN
        WRITE(YPAT,'(I2.2)') JP
      ELSE
        WRITE(YPAT,'(I1.1)') JP
      ENDIF
      YPAT = "P"//ADJUSTL(YPAT(:LEN_TRIM(YPAT)))

      YATT ='W/m2'
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'RN_'//YPAT,'Net_Radiation'//YPAT,JDIM,YATT_TITLE,YATT) 
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'H_'//YPAT,'Sensible_Heat_Flux'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LE_'//YPAT,'Total_Latent_Heat_Flux'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEI_'//YPAT,'Sublimatiob_Latent_Heat_Flux'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'GFLUX_'//YPAT,'Ground_Heat_Flux'//YPAT,JDIM,YATT_TITLE,YATT)
      !
      IF(DGO%LRAD_BUDGET) THEN
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SWD_'//YPAT,'Downward_SW'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SWU_'//YPAT,'Upward_SW'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LWD_'//YPAT,'Downward_LW'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LWU_'//YPAT,'Upward_LW'//YPAT,JDIM,YATT_TITLE,YATT)
      ENDIF
      !
      YATT ='Pa'
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FMU_'//YPAT,'Zonal_Wind_Stress'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FMV_'//YPAT,'Merid_Wind_Stress'//YPAT,JDIM,YATT_TITLE,YATT)
      !
    ENDDO
  ENDIF
  !
ENDIF
!
!
IF (DGO%LPATCH_BUDGET.AND.LAGRIP .AND. (IO%CPHOTO=='NIT' .OR. IO%CPHOTO=='LAI' &
        .OR. IO%CPHOTO=='LST' .OR. IO%CPHOTO=='NCB')) THEN
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'IRRISEUIL','Irrigation_Threshold',JDIM,YATT_TITLE,YATT)
ENDIF
!
IF (DE%LSURF_EVAP_BUDGET) THEN
  !
  YATT ='W/m2'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEG_ISBA','Averaged_Ground_Evaporation_Heat_Flux',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEGI_ISBA','Averaged_Soil_Ice_Sublimation',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEV_ISBA','Averaged_Vegetation_Evaporation_Heat_Flux',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LES_ISBA','Averaged_Snow_Sublimation_Heat_Flux',JDIM,YATT_TITLE,YATT)
  IF (TPSNOW%SCHEME=='3-L' .OR. TPSNOW%SCHEME=='CRO') THEN
     CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LESL_ISBA','Averaged_Snow_Evaporation_Heat_Flux',JDIM,YATT_TITLE,YATT)
  ENDIF
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LER_ISBA','Averaged_Canopy_Direct_Evaporation_Heat_Flux',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LETR_ISBA','Averaged_Vegetation_Transpiration_Heat_Flux',JDIM,YATT_TITLE,YATT)
  YATT ='kg/m2s'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'EVAP_ISBA','Averaged_Evapotranspiration',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SUBL_ISBA','Averaged_Sublimation_of_ice/snow',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DRAIN_ISBA','Averaged_Soil_Drainage_Flux',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'RUNOFF_ISBA','Averaged_Supersaturation_Runoff',JDIM,YATT_TITLE,YATT)
  IF (TPSNOW%SCHEME=='3-L' .OR. TPSNOW%SCHEME=='CRO') THEN
     CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SNDRIF_ISBA','Averaged_blowing_snow_sublimation',JDIM,YATT_TITLE,YATT)
  ENDIF
  IF(IO%CRUNOFF=='SGH'.AND.IO%CISBA=='DIF')THEN  
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'QSB_ISBA','Averaged_lateral_subsurface_flow',JDIM,YATT_TITLE,YATT)
  ENDIF
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'HORTON_ISBA','Averaged_Horton_Surface_Runoff',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DRIVEG_ISBA','Averaged_Dripping_from_the_vegetation_reservoir',JDIM,YATT_TITLE,YATT)
  !
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'RRVEG_ISBA','Averaged_Precipitation_Intercepted_by_Vegetation',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SNOMLT_ISBA','Averaged_Snow_melt_flux',JDIM,YATT_TITLE,YATT)
  IF(LAGRIP) CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'IRRIG_ISBA','Averaged_irrigation_rate',JDIM,YATT_TITLE,YATT)
  !
  IF (ISIZE_LMEB_PATCH>0) THEN
    YATT ='W/m2'
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEVCV_ISBA',&
            'MEB: total evapotranspiration from vegetation canopy overstory',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LESC_ISBA',&
            'MEB: total snow sublimation from vegetation canopy overstory',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LETRGV_ISBA','MEB: transpiration from understory vegetation',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LETRCV_ISBA','MEB: transpiration from overstory canopy vegetation',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LERGV_ISBA',&
            'MEB: interception evaporation from understory vegetation',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LERCV_ISBA',&
            'MEB: interception evaporation from overstory canopy vegetation',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LE_V_C_ISBA',&
            'MEB: latent heat flux from vegetation canopy overstory',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LE_G_C_ISBA','MEB: latent heat flux from understory',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LE_C_A_ISBA',&
            'MEB: latent heat flux from canopy air space to the atmosphere',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LE_N_C_ISBA','MEB: latent heat flux from the snow on the ground',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SWNET_V_ISBA','MEB: net vegetation canopy shortwave radiation',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SWNET_G_ISBA','MEB: net ground shortwave radiation',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SWNET_N_ISBA','MEB: net snow shortwave radiation',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SWNET_NS_ISBA',&
            'MEB: net snow shortwave radiation for surface layer',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LWNET_V_ISBA','MEB: net vegetation canopy longwave radiation',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LWNET_G_ISBA','MEB: net ground longwave radiation',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LWNET_N_ISBA','MEB: net snow longwave radiation',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'H_V_C_ISBA',&
            'MEB: sensible heat flux from vegetation canopy overstory',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'H_G_C_ISBA','MEB: sensible heat flux from understory',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'H_C_A_ISBA',&
            'MEB: sensible heat flux from canopy air space to the atmosphere',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'H_N_C_ISBA','MEB: sensible heat flux from the snow on the ground',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SWDOWN_GN_ISBA','MEB: SW reaching the snowpack/ground understory',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LWDOWN_GN_ISBA','MEB: LW reaching the snowpack/ground understory',JDIM,YATT_TITLE,YATT)
    YATT ='kg/m2s'
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'EVAP_N_C_ISBA',&
            'MEB: Total evap from snow on the ground to canopy air space',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'EVAP_G_C_ISBA','MEB: Total evap from ground to canopy air space',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SR_GN_ISBA','MEB: total snow reaching the ground snow',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'MELTCV_ISBA',&
            'MEB: snow melt rate from the overstory snow reservoir',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FRZCV_ISBA',&
            'MEB: snow refreeze rate from the overstory snow reservoir',JDIM,YATT_TITLE,YATT)
  ENDIF
  !
  IF(IO%LFLOOD) THEN
    YATT ='kg/m2s'
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'IFLOOD_ISBA','Averaged_Floodplains_infiltration',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'PFLOOD_ISBA','Averaged_Precipitation_intercepted_by_the floodplains',JDIM,YATT_TITLE,YATT)
    YATT ='W/m2'
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEF_ISBA','Averaged_Floodplains_evaporation_Heat_Flux',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEIF_ISBA','Averaged_Floodplains_Frozen_evaporation_Heat_Flux',JDIM,YATT_TITLE,YATT)
  ENDIF
  IF(IO%CPHOTO/='NON')THEN
    YATT ='kgCO2/m2/s'
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'GPP_ISBA','Averaged_gross_primary_production',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'R_AUTO_ISBA','Averaged_autotrophic_respiration',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'R_ECO_ISBA','Averaged_ecosystem_respiration',JDIM,YATT_TITLE,YATT)
  ENDIF
  IF(DE%LWATER_BUDGET)THEN 
    YATT ='kg/m2s'
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'RAINF_ISBA','Averaged_input_rainfall_rate',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SNOWF_ISBA','Averaged_input_snowfall_rate',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DWG_ISBA','Averaged_change_in_liquid_soil_moisture',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DWGI_ISBA','Averaged_change_in_solid_soil_moisture',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DWR_ISBA','Averaged_change_in_canopy_water',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DSWE_ISBA','Averaged_change_in_snow_water_equivalent',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WATBUD_ISBA','Averaged_isba_water_budget_as_residue',JDIM,YATT_TITLE,YATT)
  ENDIF
  !
  IF(DGO%LPATCH_BUDGET) THEN
    DO JP = 1,IO%NPATCH
      IF (JP>=10) THEN
        WRITE(YPAT,'(I2.2)') JP
      ELSE
        WRITE(YPAT,'(I1.1)') JP
      ENDIF
      YPAT = "P"//ADJUSTL(YPAT(:LEN_TRIM(YPAT)))
          
      YATT ='W/m2'
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEG_'//YPAT,'Ground_Evaporation_Heat_Flux'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEGI_'//YPAT,'Soil_Ice_Sublimation'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEV_'//YPAT,'Vegetation_Evaporation_Heat_Flux'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LES_'//YPAT,'Snow_Sublimation_Heat_Flux'//YPAT,JDIM,YATT_TITLE,YATT)
      IF (TPSNOW%SCHEME=='3-L' .OR. TPSNOW%SCHEME=='CRO') THEN
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LESL_'//YPAT,'Snow_Evaporation_Heat_Flux'//YPAT,JDIM,YATT_TITLE,YATT)
      ENDIF
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LER_'//YPAT,'Canopy_Direct_Evaporation_Heat_Flux'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LETR_'//YPAT,'Vegetation_Transpiration_Heat_Flux'//YPAT,JDIM,YATT_TITLE,YATT)
      YATT ='kg/m2s'
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'EVAP_'//YPAT,'Evapotranspiration'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SUBL_'//YPAT,'Sublimation_of_ice/snow'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DRAIN_'//YPAT,'Soil_Drainage_Flux'//YPAT,JDIM,YATT_TITLE,YATT)
      IF (TPSNOW%SCHEME=='3-L' .OR. TPSNOW%SCHEME=='CRO') THEN
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SNDRIF_'//YPAT,'blowing_snow_sublimation'//YPAT,JDIM,YATT_TITLE,YATT)
      ENDIF
      IF(IO%CRUNOFF=='SGH'.AND.IO%CISBA=='DIF')THEN  
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'QSB_'//YPAT,'lateral_subsurface_flow'//YPAT,JDIM,YATT_TITLE,YATT)
      ENDIF
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'RUNOFF_'//YPAT,'Supersaturation_Runoff'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'HORTON_'//YPAT,'Horton_Surface_Runoff'//YPAT,JDIM,YATT_TITLE,YATT)  
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DRIVEG_'//YPAT,'Dripping_from_the_vegetation_reservoir'//YPAT,JDIM,YATT_TITLE,YATT)
      !
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'RRVEG_'//YPAT,'Precipitation_Intercepted_by_Vegetation'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SNOMLT_'//YPAT,'Snow_melt_flux'//YPAT,JDIM,YATT_TITLE,YATT)
      IF(LAGRIP) CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'IRRIG_'//YPAT,'Irrigation_rate'//YPAT,JDIM,YATT_TITLE,YATT)
      !
      IF (ISIZE_LMEB_PATCH>0) THEN
        YATT ='W/m2'
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEVCV_'//YPAT,&
                'MEB: total evapotranspiration from vegetation canopy overstory'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LESC_'//YPAT,&
                'MEB: total snow sublimation from vegetation canopy overstory'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LETRGV_'//YPAT,'MEB: transpiration from understory vegetation'//YPAT,&
                JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LETRCV_'//YPAT,'MEB: transpiration from overstory canopy vegetation'//YPAT,&
                JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LERGV_'//YPAT,&
                'MEB: interception evaporation from understory vegetation'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LERCV_'//YPAT,&
                'MEB: interception evaporation from overstory canopy vegetation'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LE_V_C_'//YPAT,&
                'MEB: latent heat flux from vegetation canopy overstory'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LE_G_C_'//YPAT,'MEB: latent heat flux from understory'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LE_C_A_'//YPAT,&
                'MEB: latent heat flux from canopy air space to the atmosphere'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LE_N_C_'//YPAT,'MEB: latent heat flux from the snow on the ground'//YPAT,&
                JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SWNET_V_'//YPAT,'MEB: net vegetation canopy shortwave radiation'//YPAT,&
                JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SWNET_G_'//YPAT,'MEB: net ground shortwave radiation'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SWNET_N_'//YPAT,'MEB: net snow shortwave radiation'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SWNET_NS_'//YPAT,'MEB: net snow shortwave radiation for surface layer'//YPAT,&
                JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LWNET_V_'//YPAT,'MEB: net vegetation canopy longwave radiation'//YPAT,&
                JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LWNET_G_'//YPAT,'MEB: net ground longwave radiation'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LWNET_N_'//YPAT,'MEB: net snow longwave radiation'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'H_V_C_'//YPAT,&
                'MEB: sensible heat flux from vegetation canopy overstory'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'H_G_C_'//YPAT,'MEB: sensible heat flux from understory'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'H_C_A_'//YPAT,&
                'MEB: sensible heat flux from canopy air space to the atmosphere'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'H_N_C_'//YPAT,'MEB: sensible heat flux from the snow on the ground'//YPAT,&
                JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SWDOWN_GN_'//YPAT,'MEB: SW reaching the snowpack/ground understory'//YPAT,&
                JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LWDOWN_GN_'//YPAT,'MEB: LW reaching the snowpack/ground understory'//YPAT,&
                JDIM,YATT_TITLE,YATT)
        YATT ='kg/m2s'
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'EVAP_N_C_'//YPAT,&
                'MEB: Total evap from snow on the ground to canopy air space'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'EVAP_G_C_'//YPAT,'MEB: Total evap from ground to canopy air space'//YPAT,&
                JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SR_GN_'//YPAT,'MEB: total snow reaching the ground snow'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'MELTCV_'//YPAT,'MEB: snow melt rate from the overstory snow reservoir'//YPAT,&
                JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FRZCV_'//YPAT,&
                'MEB: snow refreeze rate from the overstory snow reservoir'//YPAT,JDIM,YATT_TITLE,YATT)
      ENDIF
      !
      IF(IO%LFLOOD)THEN
        YATT ='kg/m2s'
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'IFLOOD_'//YPAT,'Floodplains_infiltration'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'PFLOOD_'//YPAT,'Precipitation_intercepted_by_the_floodplains'//YPAT,&
                JDIM,YATT_TITLE,YATT)
        YATT ='W/m2'
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEF_'//YPAT,'Floodplains_evaporation_Heat_Flux'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEIF_'//YPAT,'Floodplains_Frozen_evaporation_Heat_Flux'//YPAT,JDIM,YATT_TITLE,YATT)
      ENDIF 
      IF(IO%CPHOTO/='NON')THEN
        YATT ='kgCO2/m2/s'
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'GPP_'//YPAT,'gross_primary_production'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'R_AUTO_'//YPAT,'autotrophic_respiration'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'R_ECO_'//YPAT,'ecosystem_respiration'//YPAT,JDIM,YATT_TITLE,YATT)
      ENDIF
      IF(DE%LWATER_BUDGET)THEN 
        YATT ='kg/m2s'
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DWG_'//YPAT,'change_in_liquid_soil_moisture'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DWGI_'//YPAT,'change_in_solid_soil_moisture'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DWR_'//YPAT,'change_in_water_on_canopy'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DSWE_'//YPAT,'change_in_snow_water_equivalent'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WATBUD_'//YPAT,'isba_water_budget_as_residue'//YPAT,JDIM,YATT_TITLE,YATT)
      ENDIF
      !
    ENDDO
    !
  ENDIF
  !
ENDIF
!
IF (DMI%LSURF_MISC_BUDGET) THEN
  !
  IL=INLVLD
  !
  DO JL=1,IL
     WRITE(YPAS,'(I3)') JL
     YLVL = TRIM(ADJUSTL(YPAS))
     CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SWI'//YLVL(:LEN_TRIM(YLVL))//'_ISBA','Soil_Wetness_Index'//YLVL,JDIM,YATT_TITLE,(/'-'/))  
     CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'TSWI'//YLVL(:LEN_TRIM(YLVL))//'_ISBA',&
             'Total_SWI_(liquid+solid)'//YLVL,JDIM,YATT_TITLE,(/'-'/))  
     IF(DGO%LPATCH_BUDGET)THEN
       DO JP = 1,IO%NPATCH
         IF (JP>=10) THEN
           WRITE(YPAT,'(I2.2)') JP
         ELSE
           WRITE(YPAT,'(I1.1)') JP
         ENDIF
         YPAT = "P"//ADJUSTL(YPAT(:LEN_TRIM(YPAT)))
         YLVL = TRIM(ADJUSTL(YPAS))//YPAT
         CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SWI'//YLVL,'Soil_Wetness_Index'//YLVL,JDIM,YATT_TITLE,(/'-'/))
         CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'TSWI'//YLVL,'Total_SWI_(liquid+solid)'//YLVL,JDIM,YATT_TITLE,(/'-'/))
       ENDDO
     ENDIF
  ENDDO
  !  
  YATT ='-'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SWI_T_ISBA','SWI_over_entire_soil',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'TSWI_T_ISBA','Total_SWI_over_entire_soil',JDIM,YATT_TITLE,YATT)
  IF(IO%CISBA=='DIF'.AND.DMI%LSURF_MISC_DIF)THEN
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'TSWI_D2_ISBA','Total_SWI_over_comparable_FR-DG2_reservoir',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'TSWI_D3_ISBA','Total_SWI_over_comparable_FR-DG3_reservoir',JDIM,YATT_TITLE,YATT)
  ENDIF  
  !
  YATT ='kg/m2'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WGTOT_T_ISBA','Total_soil_water_reservoir_(liquid+solid)',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WGI_T_ISBA','Total_soil_ice_reservoir',JDIM,YATT_TITLE,YATT)
  YATT ='m3/m3'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WGTOT_ISBA','Total_volumetric_soil_water_content_(liquid+solid)',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WGI_ISBA','Total_volumetric_soil_ice_content',JDIM,YATT_TITLE,YATT)  
  IF(IO%CISBA=='DIF'.AND.DMI%LSURF_MISC_DIF)THEN
    YATT ='m3/m3'
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WG_D2_ISBA','soil_liquid_water_over_comparable_FR-DG2_reservoir',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WGI_D2_ISBA','soil_ice_over_comparable_FR-DG2_reservoir',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WG_D3_ISBA','soil_liquid_water_comparable_FR-DG3_reservoir',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WGI_D3_ISBA','soil_ice_over_comparable_FR-DG3_reservoir',JDIM,YATT_TITLE,YATT)
  ENDIF  
  IF(IO%CISBA=='DIF')THEN
    YATT ='m'
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'ALT_ISBA','permafrost_active_layer_thickness',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FLT_ISBA','non-permafrost_frozen_layer_thickness',JDIM,YATT_TITLE,YATT)
  ENDIF
  !
  IF(IO%LGW)THEN
    YATT ='-'
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FWTD_ISBA','grid-cell_fraction_of_water_table_to_rise',JDIM,YATT_TITLE,YATT)
    YATT ='m'
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WTD_ISBA','water_table_depth_from_RRM_model_or_observation',JDIM,YATT_TITLE,YATT)
  ENDIF  
  !
  IF (TPSNOW%SCHEME=='3-L' .OR. TPSNOW%SCHEME=='CRO') THEN
    YATT ='K'
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'TS_ISBA','Surface_Temperature_(isba+snow3l)',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'TSRAD_ISBA','Surface_Radiative_Temperature_(isba+snow3l)',JDIM,YATT_TITLE,YATT)
    IF (DGO%LPATCH_BUDGET) THEN
      DO JP = 1,IO%NPATCH
        IF (JP>=10) THEN
          WRITE(YPAT,'(I2.2)') JP
        ELSE
          WRITE(YPAT,'(I1.1)') JP
        ENDIF
        YPAT = "P"//ADJUSTL(YPAT(:LEN_TRIM(YPAT)))

        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'TS_'//YPAT,'Surface_Temperature_(isba+snow3l)'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'TSRAD_'//YPAT,'total_radiative_surface_Temperature_(isba+snow3l)'//YPAT,&
                JDIM,YATT_TITLE,YATT)
      ENDDO
    ENDIF
  ENDIF
  !
  IF (TPSNOW%SCHEME=='3-L' .OR. TPSNOW%SCHEME=='CRO') THEN  
    DO JL=1,INLVLS
      WRITE(YPAS,'(I3)') JL  
      DO JP = 1,IO%NPATCH
        IF (JP>=10) THEN
          WRITE(YPAT,'(I2.2)') JP
        ELSE
          WRITE(YPAT,'(I1.1)') JP
        ENDIF
        YPAT = "P"//ADJUSTL(YPAT(:LEN_TRIM(YPAT)))
        YLVL = TRIM(ADJUSTL(YPAS))//YPAT
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SNOWTEMP'//YLVL,'Snow_Temp_layer'//YLVL,JDIM,YATT_TITLE,(/'K'/))
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SNOWLIQ'//YLVL,'Snow_liquid_water_layer_'//YLVL,JDIM,YATT_TITLE,(/'m'/))
      ENDDO
    ENDDO
  ENDIF
  ! 
  IF(IO%CRAIN=='SGH')THEN
    YATT ='-'
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'MUF_ISBA','Fraction_of_rainfall_reaching_the ground_(SGH)',JDIM,YATT_TITLE,YATT)  
  ENDIF  
  !
  YATT ='-'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'PSNG_ISBA','Snow_frac_over_ground',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'PSNV_ISBA','Snow_frac_over_veg',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'PSN_ISBA','Snow_fraction',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'TALB_ISBA','Surface total albedo',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'HV_ISBA','Halstead_coefficient',JDIM,YATT_TITLE,YATT)
  IF(IO%CPHOTO/='NON')THEN
    YATT ='kg/kg'
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LAI_ISBA','leaf_area_index',JDIM,YATT_TITLE,YATT) 
  ENDIF 
  YATT ='kg/m2'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WSN_T_ISBA','Total_snow_reservoir',JDIM,YATT_TITLE,YATT)
  YATT ='m'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DSN_T_ISBA','Total_snow_depth',JDIM,YATT_TITLE,YATT)
  YATT ='K'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'TSN_T_ISBA','Total_snow_temperature',JDIM,YATT_TITLE,YATT)  
  !
  IF(IO%CRUNOFF=='SGH'.OR.IO%CRUNOFF=='DT92')THEN
    YATT ='-'
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FSAT_ISBA','Soil_saturated_grid-cell_fraction',JDIM,YATT_TITLE,YATT)
  ENDIF   
  !
  IF(IO%LFLOOD)THEN
    YATT ='-'
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FFG_ISBA','flood_frac_over_ground',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FFV_ISBA','flood_frac_over_veg',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FF_ISBA','flood_fraction',JDIM,YATT_TITLE,YATT)
    YATT ='-'
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FFLOOD_ISBA','Potential_floodplain_grid-cell_fraction',JDIM,YATT_TITLE,YATT)
    YATT (1)='kg/m2'
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'PIFLOOD_ISBA','Potential_floodplain_infiltration',JDIM,YATT_TITLE,YATT)
  ENDIF
  ! 
  IF(DGO%LPATCH_BUDGET)THEN
    ! 
    DO JP = 1,IO%NPATCH
      IF (JP>=10) THEN
        WRITE(YPAT,'(I2.2)') JP
      ELSE
        WRITE(YPAT,'(I1.1)') JP
      ENDIF
      YPAT = "P"//ADJUSTL(YPAT(:LEN_TRIM(YPAT)))

      YATT (1)='-'
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'PSNG_'//YPAT,'snow_fraction_per_patch_over_ground'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'PSNV_'//YPAT,'snow_fraction_per_patch_over_vegetation'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'PSN_'//YPAT,'total_snow_fraction_per_patch'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'TALB_'//YPAT,'total_albedo_per_patch'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'HV_'//YPAT,'Halstead_coefficient_per_patch'//YPAT,JDIM,YATT_TITLE,YATT)
      YATT(1)='kg/m2'
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WSN_T_'//YPAT,'Total_snow_reservoir_per_patch'//YPAT,JDIM,YATT_TITLE,YATT)
      YATT(1)='m'
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DSN_T_'//YPAT,'Total_snow_depth_per_patch'//YPAT,JDIM,YATT_TITLE,YATT)
      YATT(1)='K'
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'TSN_T_'//YPAT,'Total_snow_temperature_per_patch'//YPAT,JDIM,YATT_TITLE,YATT)
      ! 
      IF(IO%CRUNOFF=='SGH'.OR.IO%CRUNOFF=='DT92')THEN
        YATT(1)='-'
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FSAT_'//YPAT,'Soil_saturated_fraction_per_patch'//YPAT,JDIM,YATT_TITLE,YATT)
      ENDIF
      !
      IF(IO%CISBA=='DIF')THEN
        YATT ='m'
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'ALT_'//YPAT,'permafrost_active_layer_thickness_per_patch'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FLT_'//YPAT,'non-permafrost_frozen_layer_thickness_per_patch'//YPAT,&
                JDIM,YATT_TITLE,YATT)
      ENDIF
      !
      IF(IO%LFLOOD)THEN
        YATT(1)='-'
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FFG_'//YPAT,'flood_frac_per_patch_over_ground'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FFV_'//YPAT,'flood_frac_per_patch_over_veg'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FF_'//YPAT,'total_flood_fraction_per_patch'//YPAT,JDIM,YATT_TITLE,YATT)
      ENDIF
      ! 
      IF (IO%LTR_ML) THEN
        YATT (1)='(-)'
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FAPAR'//YPAT,'Fapar of vegetation'//YPAT,JDIM,YATT_TITLE,YATT)
        YATT (1)='(-)'
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FAPIR'//YPAT,'Fapir of vegetation'//YPAT,JDIM,YATT_TITLE,YATT)
        YATT (1)='(-)'
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DFAPARC'//YPAT,'Fapar of vegetation (daily cumul)'//YPAT,JDIM,YATT_TITLE,YATT)
        YATT (1)='(-)'
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DFAPIRC'//YPAT,'Fapir of vegetation (daily cumul)'//YPAT,JDIM,YATT_TITLE,YATT)
        YATT (1)='(-)'
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FAPAR_BS'//YPAT,'Fapar of bare soil'//YPAT,JDIM,YATT_TITLE,YATT)
        YATT (1)='(-)'
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FAPIR_BS'//YPAT,'Fapir of bare soil'//YPAT,JDIM,YATT_TITLE,YATT)
        YATT (1)='m2/m2'
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DLAI_EFFC'//YPAT,'Effective LAI (daily cumul)'//YPAT,JDIM,YATT_TITLE,YATT)
      ENDIF
      !
      IF (LASSIM .AND. CASSIM_ISBA=="EKF  ") THEN
        YATT (1)='-'
        DO JVAR = 1,NVAR
          WRITE(YVAR,FMT='(I1.1)') JVAR
          CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'ANAL_INCR'//YVAR//YPAT,'Analysis increment'//YPAT,JDIM,YATT_TITLE,YATT)
        ENDDO
      ENDIF
      !      
    ENDDO
    !
  ENDIF
  !  
ENDIF
!
IF (CHI%SVI%NBEQ>0 .AND. CHI%CCH_DRY_DEP=="WES89 ") THEN
  !
  YATT="(m/s)"
  !
  DO JSV = 1,SIZE(CHI%CCH_NAMES,1)
    !
    YRECFM ='DV_NAT_'//TRIM(CHI%CCH_NAMES(JSV))
    WRITE(YCOMMENT,'(A7,I3.3)')'DV_NAT_',JSV
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,YRECFM,YCOMMENT,JDIM,YATT_TITLE,YATT)
    !
  ENDDO
  !
END IF
!
IF (CHI%SVI%NBEQ>0 .AND. CHI%LCH_BIO_FLUX) THEN
  !
  IF (ASSOCIATED(GB%XFISO)) THEN
    YRECFM='FISO'
    WRITE(YCOMMENT,'(A21)')'FISO (molecules/m2/s)'
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,YRECFM,YCOMMENT,JDIM,YATT_TITLE,YATT)  
  END IF
  !
  IF (ASSOCIATED(GB%XFISO)) THEN
    YRECFM='FMONO'
    WRITE(YCOMMENT,'(A22)')'FMONO (molecules/m2/s)'
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,YRECFM,YCOMMENT,JDIM,YATT_TITLE,YATT)  
  END IF
  !
ENDIF
!
IF (CHI%LCH_NO_FLUX) THEN
  !
  IF (ASSOCIATED(GB%XNOFLUX)) THEN
    YRECFM='NOFLUX'
    WRITE(YCOMMENT,'(A21)')'NOFLUX (molecules/m2/s)'
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,YRECFM,YCOMMENT,JDIM,YATT_TITLE,YATT)
  END IF
  !
END IF
!
IF(OPROVAR_TO_DIAG)THEN
  !
  IF(IO%LTEMP_ARP)THEN
    IL=IO%NTEMPLAYER_ARP
  ELSEIF(IO%CISBA=='DIF')THEN
     IL=INLVLD    
  ELSE
    IL=2
  ENDIF
  !
  YATT ='K'
  DO JL=1,IL
     WRITE(YPAS,'(I3)') JL
     YLVL=TRIM(ADJUSTL(YPAS))
     CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'TG'//YLVL(:LEN_TRIM(YLVL))//'_ISBA','Soil_temp_layer_'//YLVL,JDIM,YATT_TITLE,YATT)
  ENDDO
  !
  IL=INLVLD
  !
  YATT ='m3/m3'
  DO JL=1,IL
     WRITE(YPAS,'(I3)') JL
     YLVL=TRIM(ADJUSTL(YPAS))
     CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WG'//YLVL(:LEN_TRIM(YLVL))//'_ISBA','Soil_liquid_layer_'//YLVL,JDIM,YATT_TITLE,YATT)
  ENDDO
  YATT ='kg/m2'
  DO JL=1,IL
     WRITE(YPAS,'(I3)') JL
     YLVL=TRIM(ADJUSTL(YPAS))
     CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SOILM'//YLVL(:LEN_TRIM(YLVL))//'_ISBA',&
             'Soil_moisture_(liquid)_layer_'//YLVL,JDIM,YATT_TITLE,YATT)
  ENDDO  
  !  
  IF(IO%CISBA/='DIF')THEN
    IL=2
  ENDIF
  !
  YATT ='m3/m3'
  DO JL=1,IL
    WRITE(YPAS,'(I3)') JL
    YLVL=TRIM(ADJUSTL(YPAS))
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WGI'//YLVL(:LEN_TRIM(YLVL))//'_ISBA','Soil_ice_layer_'//YLVL,JDIM,YATT_TITLE,YATT)
  ENDDO
  YATT ='kg/m2'
  DO JL=1,IL
     WRITE(YPAS,'(I3)') JL
     YLVL=TRIM(ADJUSTL(YPAS))
     CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SOILI'//YLVL(:LEN_TRIM(YLVL))//'_ISBA','Soil_ice_mass_layer_'//YLVL,JDIM,YATT_TITLE,YATT)
  ENDDO  
  !
  YATT ='kg/m2'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WR_ISBA','Interception_reservoir',JDIM,YATT_TITLE,YATT) 
  !  
  IF(IO%LGLACIER)THEN
     YATT ='kg/m2'
     CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'ICE_STO_ISBA','Glacier_reservoir',JDIM,YATT_TITLE,YATT)
  ENDIF
  !
  YATT='-'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'ASN_ISBA','Snow_Albedo',JDIM,YATT_TITLE,YATT)
  !
  IF(TPSNOW%SCHEME=='3-L'  .OR. TPSNOW%SCHEME=='CRO')THEN
    DO JL=1,INLVLS   
      WRITE(YPAS,'(I3)') JL
      YLVL=TRIM(ADJUSTL(YPAS))
      YATT ='kg/m2'
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WSN_'//YLVL(:LEN_TRIM(YLVL))//'_ISBA',&
              'Snow_Water_Equivalent_layer_'//YLVL,JDIM,YATT_TITLE,YATT)  
      YATT ='m'
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DSN_'//YLVL(:LEN_TRIM(YLVL))//'_ISBA','Snow_Depth_layer_'//YLVL,JDIM,YATT_TITLE,YATT)  
      YATT ='K'
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'TSN_'//YLVL(:LEN_TRIM(YLVL))//'_ISBA',&
              'Snow_Temperature_layer_'//YLVL,JDIM,YATT_TITLE,YATT)
      YATT ='day_since_snowfall'
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'AGSN_'//YLVL(:LEN_TRIM(YLVL))//'_ISBA','Snow_age_layer_'//YLVL,JDIM,YATT_TITLE,YATT)
    ENDDO 
  ENDIF
  !
  IF(IO%CPHOTO=='NIT'.OR.IO%CPHOTO=='NCB')THEN
    YATT ='kgDM/m2'
    DO JNBIOMASS=1,INBIOMASS
      WRITE(YPAS,'(I3)') JNBIOMASS
      YLVL=TRIM(ADJUSTL(YPAS))
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'BIOM'//YLVL(:LEN_TRIM(YLVL))//'_ISBA','Biomass_reservoir_'//YLVL,JDIM,YATT_TITLE,YATT)
    ENDDO
  ENDIF
  !
  IF(IO%CRESPSL=='CNT')THEN
    YATT ='gC/m2'
    DO JNLITTER=1,INLITTER
      DO JNLITTLEVS=1,INLITTLEVS
        WRITE(YPAS,'(I1,A1,I1)') JNLITTER,'_',JNLITTLEVS
        YLVL = TRIM(ADJUSTL(YPAS))
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LIT'//YLVL(:LEN_TRIM(YLVL))//'_ISBA','Litter_pool'//YLVL,JDIM,YATT_TITLE,YATT)
      END DO
    END DO  
    DO JNSOILCARB=1,INSOILCARB
      WRITE(YPAS,'(I3)') JNSOILCARB
      YLVL=TRIM(ADJUSTL(YPAS))
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SCARB'//YLVL(:LEN_TRIM(YLVL))//'_ISBA','Soil_carbon_pool'//YLVL,JDIM,YATT_TITLE,YATT)
    END DO
    YATT ='-'
    DO JNLITTLEVS=1,INLITTLEVS
      WRITE(YPAS,'(I3)') JNLITTLEVS
      YLVL=TRIM(ADJUSTL(YPAS))
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LIGSTR'//YLVL(:LEN_TRIM(YLVL))//'_ISBA',&
              'Ratio_Lignin/Carbon_in_structural_litter'//YLVL,JDIM,YATT_TITLE,YATT)
    END DO
  ENDIF
ENDIF    
!
 CALL OL_WRITE_COORD(HSELECT,YFILE,IFILE_ID,JDIM,YATT_TITLE,YNAME_DIM,YUNIT1,YUNIT2,IDIM1,YDATE,ZX,ZY)
!
IF (DGO%LSURF_BUDGETC) THEN
  !
  YFILE='ISBA_DIAG_CUMUL.OUT.nc'
  CALL CREATE_FILE(YFILE,IDIMS,YNAME_DIM,IFILE_ID,JDIM)
  JRET=NF_REDEF(IFILE_ID)
  YATT(1)='dimensionless'
  !
  IF(DGO%LPATCH_BUDGET)THEN
    DO JP = 1,IO%NPATCH
      IF (JP>=10) THEN
        WRITE(YPAT,'(I2.2)') JP
      ELSE
        WRITE(YPAT,'(I1.1)') JP
      ENDIF
      YPAT = "P"//ADJUSTL(YPAT(:LEN_TRIM(YPAT)))

      YATT='J/m2'
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'RNC_'//YPAT,'Cumulated_Net_Radiation'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'HC_'//YPAT,'Cumulated_Sensible_Heat_Flux'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEC_'//YPAT,'Cumulated_Total_Latent_Heat_Flux'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEIC_'//YPAT,'Cumulated_Sublimation_Latent_Heat_Flux'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'GFLUXC_'//YPAT,'Cumulated_Ground_Heat_Flux'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEGC_'//YPAT,'Cumulated_Ground_Evaporation_Heat_Flux'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEGIC_'//YPAT,'Cumulated_Soil_Ice_Sublimation'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEVC_'//YPAT,'Cumulated_Vegetation_Evaporation_Heat_Flux'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LESC_'//YPAT,'Cumulated_Snow_Sublimation_Heat_Flux'//YPAT,JDIM,YATT_TITLE,YATT)
      IF (TPSNOW%SCHEME=='3-L' .OR. TPSNOW%SCHEME=='CRO') THEN
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LESLC_'//YPAT,'Cumulated_Snow_Evaporation_Heat_Flux'//YPAT,JDIM,YATT_TITLE,YATT)
      ENDIF
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LERC_'//YPAT,'Cumulated_Canopy_Direct_Evaporation_Heat_Flux'//YPAT,&
              JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LETRC_'//YPAT,'Cumulated_Vegetation_Transpiration_Heat_Flux'//YPAT,&
              JDIM,YATT_TITLE,YATT)
      IF(DGO%LRAD_BUDGET)THEN
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SWDC_'//YPAT,'Cumulated_Downward_SW'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SWUC_'//YPAT,'Cumulated_Upward_SW'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LWDC_'//YPAT,'Cumulated_Downward_LW'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LWUC_'//YPAT,'Cumulated_Upward_LW'//YPAT,JDIM,YATT_TITLE,YATT)
      ENDIF
      !
      YATT='Pa.s'
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FMUC_'//YPAT,'Cumulated_Zonal_Wind_Stress'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FMVC_'//YPAT,'Cumulated_Merid_Wind_Stress'//YPAT,JDIM,YATT_TITLE,YATT)  
      YATT='kg/m2'  
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'EVAPC_'//YPAT,'Cumulated_Evapotranspiration'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SUBLC_'//YPAT,'Cumulated_Sublimation_of_ice/snow'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DRAINC_'//YPAT,'Cumulated_Soil_Drainage_Flux'//YPAT,JDIM,YATT_TITLE,YATT)
      IF (TPSNOW%SCHEME=='3-L' .OR. TPSNOW%SCHEME=='CRO') THEN
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SNDRIFC_'//YPAT,'Cumulated_blowing_snow_sublimation'//YPAT,JDIM,YATT_TITLE,YATT)
      ENDIF
      IF(IO%CRUNOFF=='SGH'.AND.IO%CISBA=='DIF')THEN  
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'QSBC_'//YPAT,'Cumulated_lateral_subsurface_flow'//YPAT,JDIM,YATT_TITLE,YATT)
      ENDIF    
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'RUNOFFC_'//YPAT,'Cumulated_Supersaturation_Runoff'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'HORTONC_'//YPAT,'Cumulated_Horton_Runoff'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DRIVEGC_'//YPAT,'Cumulated_Dripping_from_the_vegetation_reservoir'//YPAT,&
              JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SNOMLTC_'//YPAT,'Cumulated_Snow_melt_flux'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'RRVEGC_'//YPAT,'Cumulated_Precipitation_Intercepted_by_Vegetation'//YPAT,&
              JDIM,YATT_TITLE,YATT)
      IF(LAGRIP) CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'IRRIGC_'//YPAT,'Cumulated_irrigation_rate'//YPAT,JDIM,YATT_TITLE,YATT)
      !
      IF(IO%LGLACIER) THEN
        YATT='kg/m2'  
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'ICE_FC_'//YPAT,'Cumulated_Glacier_ice_flux'//YPAT,JDIM,YATT_TITLE,YATT)
      ENDIF
      IF(IO%LFLOOD) THEN
        YATT='kg/m2'  
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'IFLOODC_'//YPAT,'Cumulated_Floodplains_infiltration'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'PFLOODC_'//YPAT,&
                'Cumulated_Precipitation_intercepted_by_the_floodplains'//YPAT,JDIM,YATT_TITLE,YATT)
        YATT='J/m2'
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEFC_'//YPAT,'Cumulated_Floodplains_evaporation_Heat_Flux'//YPAT,&
                JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEIFC_'//YPAT,'Cumulated_Floodplains_Frozen_evaporation_Heat_Flux'//YPAT,&
                JDIM,YATT_TITLE,YATT)
      ENDIF
      IF(IO%CPHOTO/='NON')THEN
        YATT ='kgCO2/m2'
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'GPPC_'//YPAT,'Cumulated_gross_primary_production'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'RC_AUTO_'//YPAT,'Cumulated_autotrophic_respiration'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'RC_ECO_'//YPAT,'Cumulated_ecosystem_respiration'//YPAT,JDIM,YATT_TITLE,YATT)
      ENDIF
      IF(DE%LWATER_BUDGET)THEN 
        YATT ='kg/m2'
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DWGC_'//YPAT,'Cumulated_change_in_liquid_soil_moisture'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DWGIC_'//YPAT,'Cumulated_change_in_solid_soil_moisture'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DWRC_'//YPAT,'Cumulated_change_in_canopy_water'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DSWEC_'//YPAT,'Cumulated_change_in_snow_water_equivalent'//YPAT,JDIM,YATT_TITLE,YATT)
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WATBUDC_'//YPAT,'Cumulated_isba_water_budget_as_residue'//YPAT,JDIM,YATT_TITLE,YATT)
      ENDIF
    ENDDO
  ENDIF
  !  
  YATT='J/m2'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEGC_ISBA','Averaged_Cumulated_Ground_Evaporation_Heat_Flux',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEGIC_ISBA','Averaged_Cumulated_Soil_Ice_Sublimation',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEVC_ISBA','Averaged_Cumulated_Vegetation_Evaporation_Heat_Flux',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LESC_ISBA','Averaged_Cumulated_Snow_Sublimation_Heat_Flux',JDIM,YATT_TITLE,YATT)
  IF (TPSNOW%SCHEME=='3-L' .OR. TPSNOW%SCHEME=='CRO') THEN
     CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LESLC_ISBA','Averaged_Cumulated_Snow_Evaporation_Heat_Flux',JDIM,YATT_TITLE,YATT)
  ENDIF
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LERC_ISBA','Averaged_Cumulated_Canopy_Direct_Evaporation_Heat_Flux',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LETRC_ISBA','Averaged_Cumulated_Vegetation_Transpiration_Heat_Flux',JDIM,YATT_TITLE,YATT)
  YATT='kg/m2'  
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'EVAPC_ISBA','Averaged_Cumulated_Evapotranspiration',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SUBLC_ISBA','Averaged_Cumulated_Sublimation_of_ice/snow',JDIM,YATT_TITLE,YATT)
  IF (TPSNOW%SCHEME=='3-L' .OR. TPSNOW%SCHEME=='CRO') THEN
     CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SNDRIFC_ISBA','Averaged_Cumulated_blowing_snow_sublimation',JDIM,YATT_TITLE,YATT)
  ENDIF
  IF(IO%CRUNOFF=='SGH'.AND.IO%CISBA=='DIF')THEN  
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'QSBC_ISBA','Averaged_Cumulated_lateral_subsurface_flow',JDIM,YATT_TITLE,YATT)
  ENDIF
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DRAINC_ISBA','Averaged_Cumulated_Soil_Drainage_Flux',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'RUNOFFC_ISBA','Averaged_Cumulated_Supersaturation_Runoff',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'HORTONC_ISBA','Averaged_Cumulated_Horton_Surface_Runoff',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DRIVEGC_ISBA','Averaged_Dripping_from_the_vegetation_reservoir',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SNOMLTC_ISBA','Averaged_Cumulated_Snow_melt_flux',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'RRVEGC_ISBA',&
          'Averaged_Cumulated_Precipitation_Intercepted_by_Vegetation',JDIM,YATT_TITLE,YATT)
  IF(LAGRIP) CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'IRRIGC_ISBA','Averaged_Cumulated_irrigation_rate',JDIM,YATT_TITLE,YATT)
  !
  IF(IO%LGLACIER)THEN
  YATT='kg/m2'  
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'ICE_FC_ISBA','Averaged_Cumulated_Glacier_ice_flux',JDIM,YATT_TITLE,YATT)
  ENDIF
  IF(IO%LFLOOD)THEN
  YATT='kg/m2'  
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'IFLOODC_ISBA','Averaged_Cumulated_Floodplains_infiltration',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'PFLOODC_ISBA',&
          'Averaged_Cumulated_Precip_intercepted_by_the_floodplains',JDIM,YATT_TITLE,YATT)
  YATT='J/m2'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEFC_ISBA','Averaged_Cumulated_Flood_evaporation_Heat_Flux',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEIFC_ISBA','Averaged_Cumulated_Flood_Frozen_evaporation_Heat_Flux',JDIM,YATT_TITLE,YATT)
  ENDIF 
  IF(IO%CPHOTO/='NON')THEN
    YATT ='kgCO2/m2'
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'GPPC_ISBA','Averaged_Cumulated_gross_primary_production',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'RC_AUTO_ISBA','Averaged_Cumulated_autotrophic_respiration',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'RC_ECO_ISBA','Averaged_Cumulated_ecosystem_respiration',JDIM,YATT_TITLE,YATT)
  ENDIF
  IF(DE%LWATER_BUDGET)THEN 
    YATT ='kg/m2'
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'RAINFC_ISBA','Averaged_Cumulated_input_rainfall_rate',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SNOWFC_ISBA','Averaged_Cumulated_input_snowfall_rate',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DWGC_ISBA','Averaged_Cumulated_change_in_liquid_soil_moisture',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DWGIC_ISBA','Averaged_Cumulated_change_in_solid_soil_moisture',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DWRC_ISBA','Averaged_Cumulated_change_in_canopy_water',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DSWEC_ISBA','Averaged_Cumulated_change_in_snow_water_equivalent',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WATBUDC_ISBA','Averaged_Cumulated_isba_water_budget_as_residue',JDIM,YATT_TITLE,YATT)
  ENDIF
  !
  YATT='J/m2'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'RNC_ISBA','Averaged_Cumulated_Net_Radiation',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'HC_ISBA','Averaged_Cumulated_Sensible_Heat_Flux',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEC_ISBA','Averaged_Cumulated_Total_Latent_Heat_Flux',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LEIC_ISBA','Averaged_Cumulated_Sublimation_Latent_Heat_Flux',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'GFLUXC_ISBA','Averaged_Cumulated_Ground_Heat_Flux',JDIM,YATT_TITLE,YATT)
  IF(DGO%LRAD_BUDGET)THEN
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SWDC_ISBA','Averaged_Cumulated_Downward_SW',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'SWUC_ISBA','Averaged_Cumulated_Upward_SW',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LWDC_ISBA','Averaged_Cumulated_Downward_LW',JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LWUC_ISBA','Averaged_Cumulated_Upward_LW',JDIM,YATT_TITLE,YATT)
  ENDIF
  YATT='Pa.s'
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FMUC_ISBA','Averaged_Cumulated_Zonal_Wind_Stress',JDIM,YATT_TITLE,YATT)
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FMVC_ISBA','Averaged_Cumulated_Merid_Wind_Stress',JDIM,YATT_TITLE,YATT)
  !  
  CALL OL_WRITE_COORD(HSELECT,YFILE,IFILE_ID,JDIM,YATT_TITLE,YNAME_DIM,YUNIT1,YUNIT2,IDIM1,YDATE,ZX,ZY)
  !
ENDIF


! 6. Create file for vegetation parameter values
!----------------------------------------------------------

IF(LASSIM) THEN
  IF(CASSIM=='PLUS') THEN
    YFILE='ISBA_VEG_EVOLUTION_P.OUT.nc'
    CALL CREATE_FILE(YFILE,IDIMS,YNAME_DIM,IFILE_ID,JDIM)
    JRET=NF_REDEF(IFILE_ID)
    YATT='dimensionless'
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LAIp','Output_LAI_ISBA',JDIM,YATT_TITLE,YATT)
    CALL OL_WRITE_COORD(HSELECT,YFILE,IFILE_ID,JDIM,YATT_TITLE,YNAME_DIM,YUNIT1,YUNIT2,IDIM1,YDATE,ZX,ZY)
  ELSEIF(CASSIM=='AVERA') THEN
    YFILE='ISBA_VEG_EVOLUTION_A.OUT.nc'
    CALL CREATE_FILE(YFILE,IDIMS,YNAME_DIM,IFILE_ID,JDIM)
    JRET=NF_REDEF(IFILE_ID)
    YATT ='dimensionless'
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LAIa','Output_LAI_ISBA',JDIM,YATT_TITLE,YATT)
    CALL OL_WRITE_COORD(HSELECT,YFILE,IFILE_ID,JDIM,YATT_TITLE,YNAME_DIM,YUNIT1,YUNIT2,IDIM1,YDATE,ZX,ZY)
  ELSEIF(CASSIM=='2DVAR') THEN
    YFILE='ISBA_VEG_EVOLUTION.OUT.nc'
    CALL CREATE_FILE(YFILE,IDIMS,YNAME_DIM,IFILE_ID,JDIM)
    JRET=NF_REDEF(IFILE_ID)
    YATT='dimensionless'
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LAI','Output_LAI_ISBA',JDIM,YATT_TITLE,YATT)
    CALL OL_WRITE_COORD(HSELECT,YFILE,IFILE_ID,JDIM,YATT_TITLE,YNAME_DIM,YUNIT1,YUNIT2,IDIM1,YDATE,ZX,ZY)
  ENDIF
ELSEIF(DGO%LPGD)THEN
  YFILE='ISBA_VEG_EVOLUTION.OUT.nc'
  CALL CREATE_FILE(YFILE,IDIMS,YNAME_DIM,IFILE_ID,JDIM)  
  JRET=NF_REDEF(IFILE_ID)
  !
  DO JP = 1,IO%NPATCH
    IF (JP>=10) THEN
      WRITE(YPAT,'(I2.2)') JP
    ELSE
      WRITE(YPAT,'(I1.1)') JP
    ENDIF
    YPAT = "P"//ADJUSTL(YPAT(:LEN_TRIM(YPAT)))  
    !
    YATT ='dimensionless'
    !
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'VEG'//YPAT,'Output_vegetation_fraction'//YPAT,JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'LAI'//YPAT,'Output_LAI_per_patch'//YPAT,JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'Z0VEG'//YPAT,'Roughness_Length_Vegetation'//YPAT,JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'PATCH'//YPAT,'Fraction_Of_Patch'//YPAT,JDIM(1:INDIMS-1),YATT_TITLE,YATT)
    !
    DO JL=1,INLVLD
      WRITE(YPAS,'(I3)') JL ; YLVL=TRIM(ADJUSTL(YPAS))//YPAT
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DG'//YLVL,'soil_depth_layer_'//YLVL,JDIM(1:INDIMS-1),YATT_TITLE,(/'m'/))
    ENDDO
    
  ENDDO
  !
  IF(INPATCH>1)THEN
    DO JL=1,INLVLD
      WRITE(YPAS,'(I3)') JL ; YLVL=TRIM(ADJUSTL(YPAS))
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DG'//YLVL(:LEN_TRIM(YLVL))//'_ISBA',&
              'averaged_soil_depth_layer_'//YLVL,JDIM(1:1),YATT_TITLE,(/'m'/))
    ENDDO
  ENDIF
  !
  CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'Z0REL','orography_roughness_length',JDIM(1:1),YATT_TITLE,(/'m'/))
  !
  !
  DO JL=1,INLVLD
    WRITE(YPAS,'(I3)') JL ; YLVL=TRIM(ADJUSTL(YPAS))
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WSAT'//YLVL,'soil_porosity_layer_'//YLVL,JDIM(1:1),YATT_TITLE,(/'m3/m3'/))
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WFC'//YLVL,'field_capacity_layer_'//YLVL,JDIM(1:1),YATT_TITLE,(/'m3/m3'/)) 
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WWILT'//YLVL,'wilting_point_layer_'//YLVL,JDIM(1:1),YATT_TITLE,(/'m3/m3'/))
  ENDDO
  !
  IF(IO%CISBA=='DIF')THEN
    !
    DO JP = 1,IO%NPATCH
      IF (JP>=10) THEN
        WRITE(YPAT,'(I2.2)') JP
      ELSE
        WRITE(YPAT,'(I1.1)') JP
      ENDIF
      YPAT = "P"//ADJUSTL(YPAT(:LEN_TRIM(YPAT)))
      !
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DROOT_DIF'//YPAT,'Root_depth_in_ISBA-DIF'//YPAT,JDIM(1:INDIMS-1),YATT_TITLE,(/'m'/))
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DG2_DIF'//YPAT,'DG2_depth_in_ISBA-DIF'//YPAT,JDIM(1:INDIMS-1),YATT_TITLE,(/'m'/))
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'RUNOFFD'//YPAT,'Runoff_depth_in_ISBA-DIF'//YPAT,JDIM(1:INDIMS-1),YATT_TITLE,(/'m'/))
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DTOT_DIF'//YPAT,&
              'Total_soil_depth_for_moisture_in_ISBA-DIF'//YPAT,JDIM(1:INDIMS-1),YATT_TITLE,(/'m'/))
      !
      DO JL=1,INLVLD
        WRITE(YPAS,'(I3)') JL ; YLVL=TRIM(ADJUSTL(YPAS))//YPAT
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'ROOTFRAC'//YLVL,'root_fraction_layer_'//YLVL,JDIM(1:INDIMS-1),YATT_TITLE,(/'-'/))
      ENDDO
      !    
    ENDDO
    !
    IF(INPATCH>1)THEN
     CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DG2_DIF_ISBA','averaged_DG2_depth_in_ISBA-DIF',JDIM(1:1),YATT_TITLE,(/'m'/))
     CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DTOTDF_ISBA',&
             'averaged_Total_soil_depth_for_moisture_in_ISBA-DIF',JDIM(1:1),YATT_TITLE,(/'m'/))
    ENDIF
    !
    IF(IO%LSOC)THEN
      DO JL=1,INLVLD
        WRITE(YPAS,'(I3)') JL ; YLVL=TRIM(ADJUSTL(YPAS))
        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'FRACSOC'//YLVL,'SOC_fraction_layer_'//YLVL,JDIM(1:1),YATT_TITLE,(/'-'/))
      ENDDO
    ENDIF
    !
  ENDIF
  !
  IF(IO%CHORT=='SGH')THEN
    DO JP = 1,IO%NPATCH
      IF (JP>=10) THEN
        WRITE(YPAT,'(I2.2)') JP
      ELSE
        WRITE(YPAT,'(I1.1)') JP
      ENDIF
      YPAT = "P"//ADJUSTL(YPAT(:LEN_TRIM(YPAT)))
      !          
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'DICE'//YPAT,'soil_ice_depth_for_runoff'//YPAT,JDIM(1:INDIMS-1),YATT_TITLE,(/'m'/))
    ENDDO
  ENDIF   
  !    
  DO JVEG=1,NVEGTYPE
     WRITE(YPAS,'(i2)') JVEG 
     YLVLV=TRIM(ADJUSTL(YPAS))
     CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'VEGTYPE'//YLVLV,'fraction_of_vegtype_in_the_grid_cell',JDIM(1:1),YATT_TITLE,(/'-'/))
  ENDDO
  !    
  IF(INPATCH>1.AND.NVEGTYPE/=INPATCH)THEN
    DO JP = 1,IO%NPATCH
      IF (JP>=10) THEN
        WRITE(YPAT,'(I2.2)') JP
      ELSE
        WRITE(YPAT,'(I1.1)') JP
      ENDIF
      YPAT = "P"//ADJUSTL(YPAT(:LEN_TRIM(YPAT)))
      DO JVEG=1,NVEGTYPE
        WRITE(YPAS,'(i2)') JVEG 
        YLVLV=TRIM(ADJUSTL(YPAS))//YPAT

        CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'VEGTY_'//YLVLV,&
                'fraction_of_vegtype_in_each_patch'//YLVLV,JDIM(1:INDIMS-1),YATT_TITLE,(/'-'/))
      ENDDO
    ENDDO
  ENDIF
  !
  DO JP = 1,IO%NPATCH
    IF (JP>=10) THEN
      WRITE(YPAT,'(I2.2)') JP
    ELSE
      WRITE(YPAT,'(I1.1)') JP
    ENDIF
    YPAT = "P"//ADJUSTL(YPAT(:LEN_TRIM(YPAT)))
  
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'EMIS'//YPAT,'Emissivity_Of_Vegetation'//YPAT,JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'RSMIN'//YPAT,'Minimal_Stomatal_Resistance'//YPAT,JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'GAMMA'//YPAT,'Coefficient_Computation_Rsmin'//YPAT,JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'CV'//YPAT,'Vegetal_Thermal_Inertia'//YPAT,JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'RGL'//YPAT,'Max_Solar_Radiation_Photosynthesis'//YPAT,JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WRMAX_CF'//YPAT,'Coefficient_Max_Water_Interception'//YPAT,JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'ALBNIR_SOIL'//YPAT,'Output_ALBNIR_SOIL'//YPAT,JDIM(1:INDIMS-1),YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'ALBVIS_SOIL'//YPAT,'Output_ALBVIS_SOIL'//YPAT,JDIM(1:INDIMS-1),YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'ALBUV_SOIL'//YPAT,'soil_UV_albedo'//YPAT,JDIM(1:INDIMS-1),YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'ALBNIR_ISBA'//YPAT,'total_near-infra-red albedo'//YPAT,JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'ALBVIS_ISBA'//YPAT,'total_visible_albedo'//YPAT,JDIM,YATT_TITLE,YATT)
    CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'ALBUV_ISBA'//YPAT,'total_UV_albedo'//YPAT,JDIM,YATT_TITLE,YATT)
    !
    IF (ISIZE_LMEB_PATCH>0) THEN
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'H_VEG'//YPAT,'MEB: height_of_vegetation'//YPAT,JDIM,YATT_TITLE,YATT)
    ENDIF
    !  
    IF (LAGRIP .AND. (IO%CPHOTO=='NIT' .OR. IO%CPHOTO=='LAI' .OR. IO%CPHOTO=='LST' .OR. IO%CPHOTO=='NCB') ) THEN
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'WATSUP'//YPAT,'Water_Supply_Irrigation'//YPAT,JDIM,YATT_TITLE,YATT)
      CALL DEF_VAR_NETCDF(HSELECT,IFILE_ID,'IRRIG'//YPAT,'Fraction_Of_Irrigated_Vegetation'//YPAT,JDIM,YATT_TITLE,YATT)
    END IF
    !
  ENDDO
  !
  CALL OL_WRITE_COORD(HSELECT,YFILE,IFILE_ID,JDIM,YATT_TITLE,YNAME_DIM,YUNIT1,YUNIT2,IDIM1,YDATE,ZX,ZY)
  !
ENDIF
!
DEALLOCATE(JDIM)
!
IF (LHOOK) CALL DR_HOOK('INIT_OUTFN_ISBA_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE INIT_OUTFN_ISBA_n
