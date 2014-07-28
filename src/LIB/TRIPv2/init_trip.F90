!####################################################################################
SUBROUTINE INIT_TRIP (KYEAR,KMONTH,KDAY,PTIME,KLON,KLAT,PTSTEP_RUN,PTSTEP_DIAG,ORESTART)
!####################################################################################
!
!!****  *INIT_TRIP*  
!!
!!    PURPOSE
!!    -------
!
!     Initialize TRIP variables and parameters.
!     
!!**  METHOD
!!    ------
!
!     Direct calculation
!
!!    EXTERNAL
!!    --------
!
!     None
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
!!	B. Decharme     
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/02/05 
!!      For surfex  21/05/08 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODN_TRIP, ONLY : CGROUNDW, CVIT, LFLOOD,  &
                      XCVEL, XRATMED, XTSTEP
!
USE MODD_TRIP_PAR
USE MODD_TRIP_LISTING, ONLY : NLISTING
!
USE MODD_TRIP_GRID, ONLY : TGRID
USE MODD_TRIP,      ONLY : TTRIP
!
USE MODE_TRIP_GRID
USE MODE_TRIP_INIT
USE MODE_TRIP_FUNCTION
USE MODE_RW_TRIP
!
USE MODI_ABORT_TRIP
!USE MODI_READ_TRIP_GRID
USE MODI_INIT_TRIP_DIAG
USE MODI_INIT_RESTART_TRIP
USE MODI_GET_LONLAT_TRIP
USE MODI_INIT_TRIP_CPL_ESM
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!-------------------------------------------------------------------------------
!
!*      0.1    declarations of arguments
!
INTEGER,          INTENT(OUT) :: KYEAR   !date UTC
INTEGER,          INTENT(OUT) :: KMONTH  !date UTC
INTEGER,          INTENT(OUT) :: KDAY    !date UTC
REAL,             INTENT(OUT) :: PTIME   !date UTC
!
INTEGER,          INTENT(OUT) :: KLON      ! number of points in longitude
INTEGER,          INTENT(OUT) :: KLAT      ! number of points in latitude
!
REAL,             INTENT(IN) :: PTSTEP_RUN
REAL,             INTENT(IN) :: PTSTEP_DIAG
!
LOGICAL,          INTENT(IN) :: ORESTART
!
!-------------------------------------------------------------------------------
!
!*      0.2    declarations of local variables
!
CHARACTER(LEN=13), PARAMETER         :: YFILE_PARAM  ='TRIP_PARAM.nc'
CHARACTER(LEN=12), PARAMETER         :: YFILE_INIT   ='TRIP_PREP.nc'
CHARACTER(LEN=15), PARAMETER         :: YFILE_RESTART='TRIP_RESTART.nc'
CHARACTER(LEN=19), PARAMETER         :: YDIAG        ='TRIP_DIAG.nc'
CHARACTER(LEN=18), PARAMETER         :: YRUN         ='TRIP_DIAG_RUN.nc'
! 
CHARACTER(LEN=6)                     :: YTIME
CHARACTER(LEN=50)                    :: YFILE
CHARACTER(LEN=20)                    :: YVAR 
!
REAL,DIMENSION(4)                    :: ZDATE
!
REAL,DIMENSION(:,:),ALLOCATABLE      :: ZREAD
REAL,DIMENSION(:,:),ALLOCATABLE      :: ZHSTREAM
REAL,DIMENSION(:,:),ALLOCATABLE      :: ZWORK
REAL,DIMENSION(:,:),ALLOCATABLE      :: ZHG_OLD
REAL,DIMENSION(:,:),ALLOCATABLE      :: ZWTD
REAL,DIMENSION(:,:),ALLOCATABLE      :: ZFWTD
REAL,DIMENSION(:,:),ALLOCATABLE      :: ZWTDELEV
!
REAL, DIMENSION(:),ALLOCATABLE       :: ZLON
REAL, DIMENSION(:),ALLOCATABLE       :: ZLAT
!
REAL    :: ZLONMIN   ! minimum longitude (degrees)
REAL    :: ZLONMAX   ! maximum longitude (degrees)
REAL    :: ZLATMIN   ! minimum latitude  (degrees)
REAL    :: ZLATMAX   ! maximum latitude  (degrees)
REAL    :: ZGRID_RES ! 1° or 0.5° resolution
!
INTEGER :: IWORK, IFLOOD, INI, JLON, JLAT
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
! * Output attribut for netcdf diag file
!-------------------------------------------------------------------------------
!
CHARACTER(LEN=50) :: YTITLE, YUNITTIME
!                    YTITLE    = Title of each output file
!                    YUNITTIME = Time unit in each output file if present
!
!-------------------------------------------------------------------------------
! * Initilyse TRIP
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('INIT_TRIP',0,ZHOOK_HANDLE)
!
WRITE(NLISTING,*)''
WRITE(NLISTING,*)''
WRITE(NLISTING,*)'!'
WRITE(NLISTING,*)'! TRIP RUN !!!!!!!!!!!!!'
WRITE(NLISTING,*)'!'
!
WRITE(NLISTING,*)''
WRITE(NLISTING,*)'        INITIALYSE TRIP            '
WRITE(NLISTING,*)''
!
!-------------------------------------------------------------------------------
! * Read date
!-------------------------------------------------------------------------------
!
YVAR='date'
CALL READ_TRIP(NLISTING,YFILE_INIT,YVAR,ZDATE)
!
KYEAR  = INT(ZDATE(1))
KMONTH = INT(ZDATE(2))
KDAY   = INT(ZDATE(3))
PTIME  = ZDATE(4)
!
!-------------------------------------------------------------------------------
!
! * Get TRIP grid configuration
!
CALL GET_TRIP_GRID(TGRID%XTRIP_GRID,ZLONMIN,ZLONMAX,ZLATMIN,ZLATMAX,ZGRID_RES,KLON,KLAT)
!
ALLOCATE(ZLON(KLON))
ALLOCATE(ZLAT(KLAT))
ZLON(:)=XUNDEF
ZLAT(:)=XUNDEF
!
CALL GET_LONLAT_TRIP(KLON,KLAT,ZLON,ZLAT)
!
!-------------------------------------------------------------------------------
! * Check options
!-------------------------------------------------------------------------------
!
IF(CVIT == 'VAR')THEN
   IF(XCVEL/=0.5)THEN
      XCVEL=0.5
      WRITE(NLISTING,*)'!!! ATTENTION : You use the velocity scheme and XCVEL is not 0.5 m/s !!!'
      WRITE(NLISTING,*)'!!! ATTENTION : XCVEL put at 0.5 m/s !!!'
  ENDIF
ELSE
   IF(XCVEL<0.1)THEN
      WRITE(NLISTING,*)'!!! ATTENTION : XCVEL < 0.1 m/s !!! Not good !!!'
      CALL ABORT_TRIP('INIT_TRIP: ATTENTION : XCVEL < 0.1 m/s !!! Not good !!!')
  ENDIF     
ENDIF      
!
IF(CGROUNDW=='DIF')THEN
  IF(CVIT /= 'VAR')THEN
    WRITE(NLISTING,*)'! You cannot use the groundwater scheme without the variable velocity scheme !!!'
    CALL ABORT_TRIP('INIT_TRIP: You cannot use the groundwater scheme without the variable velocity scheme !!!')
  ENDIF   
  IF(ZGRID_RES>0.5)THEN
    WRITE(NLISTING,*)'! You cannot use the groundwater scheme with another resolution than 0.5 or 1/12 !!!'
    CALL ABORT_TRIP('INIT_TRIP: You cannot use the groundwater scheme with another resolution than 0.5 or 1/12 !!!')
  ENDIF
ENDIF
!
IF(LFLOOD)THEN
  IF(CGROUNDW=='DEF')THEN
    WRITE(NLISTING,*)'! ATTENTION : You use the flooding scheme without groundwater scheme !!!'
  ENDIF
  IF(CVIT /= 'VAR')THEN
    WRITE(NLISTING,*)'! You cannot use the flooding scheme without the variable velocity scheme !!!'
    CALL ABORT_TRIP('INIT_TRIP: You cannot use the flooding scheme without the variable velocity scheme !!!')
  ENDIF
  IF(XTSTEP>1800.)THEN
    WRITE(NLISTING,*)'!'
    WRITE(NLISTING,*)'! For flooding, the TRIP time step is too big      !!!'
    WRITE(NLISTING,*)'! XTSTEP must be equal or inferior to 1800s   !!!'
    WRITE(NLISTING,*)'!'
    CALL ABORT_TRIP('INIT_TRIP: For flooding, the TRIP time step is too big      !!!')
  ENDIF
ENDIF
!
IF(MOD(PTSTEP_RUN,XTSTEP)*MOD(XTSTEP,PTSTEP_RUN)/=0.)THEN
  WRITE(NLISTING,*)'! XTSTEP_RUN and XTSTEP are not good !!!'     
  WRITE(NLISTING,*)'! XTSTEP_RUN =', PTSTEP_RUN    
  WRITE(NLISTING,*)'! XTSTEP     =', XTSTEP    
  CALL ABORT_TRIP('INIT_TRIP: PTSTEP_RUN and XTSTEP are not good !!!')
ENDIF
!
IF(MOD(PTSTEP_DIAG,XTSTEP)*MOD(XTSTEP,PTSTEP_DIAG)/=0.)THEN
  WRITE(NLISTING,*)'! PTSTEP_DIAG and XTSTEP are not good !!!'  
  WRITE(NLISTING,*)'! XTSTEP_DIAG =', PTSTEP_DIAG    
  WRITE(NLISTING,*)'! XTSTEP      =', XTSTEP  
  CALL ABORT_TRIP('INIT_TRIP: PTSTEP_DIAG and XTSTEP are not good !!!')
ENDIF
!
WRITE(NLISTING,*)'! ',ZGRID_RES,'° TRIP run !!!'  
!
IF(ZGRID_RES<0.5.AND.XRATMED==1.4)THEN
     WRITE(NLISTING,*)'! meandering ratio is 1.4 at 0.5° or 1° resolution !!!'   
     WRITE(NLISTING,*)'! for other resolution change XRATMED in namelist  !!!' 
     CALL ABORT_TRIP('INIT_TRIP: meandering ratio is 1.4 at 0.5° or 1° resolution !!!')
ENDIF
!
!-------------------------------------------------------------------------------
! * Allocate trip grid arguments
!-------------------------------------------------------------------------------
!
ALLOCATE(TGRID%XAREA    (KLON,KLAT))
ALLOCATE(TGRID%XLEN     (KLON,KLAT))
ALLOCATE(TGRID%NGRCN    (KLON,KLAT))       
ALLOCATE(TGRID%NSEQ     (KLON,KLAT))   
ALLOCATE(TGRID%NNEXTX   (KLON,KLAT))    
ALLOCATE(TGRID%NNEXTY   (KLON,KLAT))   
ALLOCATE(TGRID%NBASID   (KLON,KLAT))
!
ALLOCATE(TGRID%GMASK    (KLON,KLAT))
ALLOCATE(TGRID%GMASK_VEL(KLON,KLAT))
ALLOCATE(TGRID%GMASK_FLD(KLON,KLAT))
ALLOCATE(TGRID%GMASK_GW (KLON,KLAT))
ALLOCATE(TGRID%GMASK_GRE(KLON,KLAT))
ALLOCATE(TGRID%GMASK_ANT(KLON,KLAT))
!
TGRID%GMASK    (:,:) = .FALSE.
TGRID%GMASK_VEL(:,:) = .FALSE.
TGRID%GMASK_FLD(:,:) = .FALSE.
TGRID%GMASK_GW (:,:) = .FALSE.
TGRID%GMASK_GRE(:,:) = .FALSE.
TGRID%GMASK_ANT(:,:) = .FALSE.
!
!-------------------------------------------------------------------------------
! * Allocate trip arguments
!-------------------------------------------------------------------------------
!
ALLOCATE(TTRIP%XSURF_STO(KLON,KLAT))
!
IF(CGROUNDW/='DEF')THEN
  ALLOCATE(TTRIP%XTAUG(KLON,KLAT))
ELSE
  ALLOCATE(TTRIP%XTAUG(0,0))
ENDIF
!
IF(CGROUNDW=='CST')THEN
  ALLOCATE(TTRIP%XGROUND_STO(KLON,KLAT))      
ELSE
  ALLOCATE(TTRIP%XGROUND_STO(0,0))      
ENDIF
!
IF(CGROUNDW=='DIF')THEN
  ALLOCATE(TTRIP%XHGROUND (KLON,KLAT))
  ALLOCATE(TTRIP%XWEFF    (KLON,KLAT))
  ALLOCATE(TTRIP%XTRANS   (KLON,KLAT))
  ALLOCATE(TTRIP%XNUM_AQUI(KLON,KLAT))
  ALLOCATE(TTRIP%XELEV    (KLON,KLAT))
  ALLOCATE(TTRIP%XTOPO_RIV(KLON,KLAT))  
ELSE
  ALLOCATE(TTRIP%XHGROUND (0,0))
  ALLOCATE(TTRIP%XWEFF    (0,0))
  ALLOCATE(TTRIP%XTRANS   (0,0))
  ALLOCATE(TTRIP%XNUM_AQUI(0,0))
  ALLOCATE(TTRIP%XELEV    (0,0))
  ALLOCATE(TTRIP%XTOPO_RIV(0,0))  
ENDIF
!
IF(LFLOOD.OR.CGROUNDW=='DIF')THEN
  ALLOCATE(TTRIP%XHC_BED(KLON,KLAT))
ELSE
  ALLOCATE(TTRIP%XHC_BED(0,0))
ENDIF
!
IF(LFLOOD)THEN
  ALLOCATE(TTRIP%XN_FLOOD  (KLON,KLAT))      
  ALLOCATE(TTRIP%XFLOOD_STO(KLON,KLAT))
  ALLOCATE(TTRIP%XWFLOOD   (KLON,KLAT))
  ALLOCATE(TTRIP%XFLOOD_LEN(KLON,KLAT))
  ALLOCATE(TTRIP%XHFLOOD   (KLON,KLAT))
  ALLOCATE(TTRIP%XFFLOOD   (KLON,KLAT))
ELSE
  ALLOCATE(TTRIP%XN_FLOOD  (0,0))      
  ALLOCATE(TTRIP%XFLOOD_STO(0,0))
  ALLOCATE(TTRIP%XWFLOOD   (0,0))
  ALLOCATE(TTRIP%XFLOOD_LEN(0,0))
  ALLOCATE(TTRIP%XHFLOOD   (0,0))
  ALLOCATE(TTRIP%XFFLOOD   (0,0))
ENDIF
!
IF(CVIT=='VAR')THEN
  ALLOCATE(TTRIP%XSLOPEBED(KLON,KLAT))      
  ALLOCATE(TTRIP%XWIDTH   (KLON,KLAT))      
  ALLOCATE(TTRIP%XN       (KLON,KLAT))
ELSE
  ALLOCATE(TTRIP%XSLOPEBED(0,0))      
  ALLOCATE(TTRIP%XWIDTH   (0,0))      
  ALLOCATE(TTRIP%XN       (0,0))
ENDIF
!
!-------------------------------------------------------------------------------
! * Compute and Read TRIP parameter
!-------------------------------------------------------------------------------
!
ALLOCATE(ZREAD(KLON,KLAT))
!
! * Flow direction 
!
YVAR ='FLOWDIR'
CALL READ_TRIP(NLISTING,YFILE_PARAM,YVAR,ZREAD)
WHERE(ZREAD==XUNDEF)ZREAD=0.0
TGRID%NGRCN(:,:)=INT(ZREAD(:,:))
WHERE(TGRID%NGRCN(:,:)>0)TGRID%GMASK(:,:)=.TRUE.
!
! * Rriver sequence
!
YVAR ='RIVSEQ'
CALL READ_TRIP(NLISTING,YFILE_PARAM,YVAR,ZREAD)
WHERE(ZREAD==XUNDEF)ZREAD=0.0
TGRID%NSEQ(:,:)=INT(ZREAD(:,:))
!
! * Maximum river sequence value
!
TGRID%NSEQMAX = MAXVAL(TGRID%NSEQ(:,:))
!
! * Basin number id
!
YVAR ='NUM_BAS'
CALL READ_TRIP(NLISTING,YFILE_PARAM,YVAR,ZREAD)
WHERE(ZREAD==XUNDEF)ZREAD=0.0
TGRID%NBASID(:,:)=INT(ZREAD(:,:))
!
TGRID%NBASMIN = MINVAL(TGRID%NBASID(:,:),TGRID%NBASID(:,:)>0)
TGRID%NBASMAX = MAXVAL(TGRID%NBASID(:,:),TGRID%NBASID(:,:)>0)
!
! * Set down stream
!
CALL SETNEXT(KLON,KLAT,TGRID%NGRCN,TGRID%NNEXTX,TGRID%NNEXTY)
!
! * Set area size
!
CALL SETAREA(KLAT,ZLATMIN,ZGRID_RES,TGRID%XAREA)
!
! * Distance between grids with the meandering ratio
!
YVAR ='RIVLEN'
CALL READ_TRIP(NLISTING,YFILE_PARAM,YVAR,TGRID%XLEN)
WHERE(.NOT.TGRID%GMASK(:,:))TGRID%XLEN(:,:)=XUNDEF
!
! * Land mask for Greenland and Antarctica
!
YVAR ='GREEN_ANT'
CALL READ_TRIP(NLISTING,YFILE_PARAM,YVAR,ZREAD)
DO JLAT=1,KLAT
   DO JLON=1,KLON
      IF(ZREAD(JLON,JLAT)>0.0.AND.(ZLAT(JLAT)<0.0))THEN
          TGRID%GMASK_ANT(JLON,JLAT)=.TRUE.
      ELSEIF(ZREAD(JLON,JLAT)>0.0.AND.(ZLAT(JLAT)>=0.0))THEN
          TGRID%GMASK_GRE(JLON,JLAT)=.TRUE.
      ENDIF
   ENDDO
ENDDO
!
DEALLOCATE(ZREAD)
DEALLOCATE(ZLON)
DEALLOCATE(ZLAT)
!
! * Variable velocity schemes variables
!
IF(CVIT == 'VAR')THEN
!
  YVAR ='N_RIV'
  CALL READ_TRIP(NLISTING,YFILE_PARAM,YVAR,TTRIP%XN)
!
  YVAR ='WIDTHRIV'
  CALL READ_TRIP(NLISTING,YFILE_PARAM,YVAR,TTRIP%XWIDTH)
!
  YVAR ='SLOPERIV'
  CALL READ_TRIP(NLISTING,YFILE_PARAM,YVAR,TTRIP%XSLOPEBED)
!  
  TGRID%GMASK_VEL(:,:)=TGRID%GMASK(:,:)
  WHERE(TGRID%NSEQ(:,:)==0.OR.TTRIP%XWIDTH(:,:)>=XUNDEF-1.0)
        TTRIP%XSLOPEBED(:,:)= 0.0      
        TTRIP%XWIDTH   (:,:)= 0.0      
        TTRIP%XN       (:,:)= 0.0
        TGRID%GMASK_VEL(:,:)=.FALSE.
  ENDWHERE
!
ENDIF 
!
! * Groundwater parameters
!
IF(CGROUNDW/='DEF')THEN
!
! * Calculate the groundwater transfert time
!
  YVAR ='TAUG'
  CALL READ_TRIP(NLISTING,YFILE_PARAM,YVAR,TTRIP%XTAUG)
!
  WHERE(TGRID%NSEQ(:,:)==0)
        TTRIP%XTAUG(:,:)=XUNDEF
  ENDWHERE
!
  TGRID%GMASK_GW(:,:)=TGRID%GMASK(:,:)
!
  WHERE(TTRIP%XTAUG(:,:)/=XUNDEF)
        TTRIP%XTAUG(:,:)=TTRIP%XTAUG(:,:)*XDAY
        TGRID%GMASK_GW(:,:)=.TRUE.
  ELSEWHERE
        TGRID%GMASK_GW(:,:)=.FALSE.
  ENDWHERE
!
ENDIF
!
IF(CGROUNDW == 'DIF')THEN
!
! * Diffusive groundwater scheme variables
!
!
  YVAR = 'NUM_AQUI'
  CALL READ_TRIP(NLISTING,YFILE_PARAM,YVAR,TTRIP%XNUM_AQUI)
!
  WHERE(TTRIP%XNUM_AQUI(:,:)==XUNDEF)
        TGRID%GMASK_GW(:,:)=.FALSE.
  ENDWHERE
!
  YVAR = 'WEFF'
  CALL READ_TRIP(NLISTING,YFILE_PARAM,YVAR,TTRIP%XWEFF)
!
  YVAR = 'TRANS'
  CALL READ_TRIP(NLISTING,YFILE_PARAM,YVAR,TTRIP%XTRANS)
!
  YVAR = 'ELEV'
  CALL READ_TRIP(NLISTING,YFILE_PARAM,YVAR,TTRIP%XELEV)
!
  YVAR = 'TOPO_RIV'
  CALL READ_TRIP(NLISTING,YFILE_PARAM,YVAR,TTRIP%XTOPO_RIV)
!
  ALLOCATE(TTRIP%XTABGW_F (KLON,KLAT,NDIMTAB))      
  ALLOCATE(TTRIP%XTABGW_H (KLON,KLAT,NDIMTAB))      
!
  YVAR ='TABGW_F'
  CALL READ_TRIP(NLISTING,YFILE_PARAM,YVAR,TTRIP%XTABGW_F)
!
  YVAR ='TABGW_H'
  CALL READ_TRIP(NLISTING,YFILE_PARAM,YVAR,TTRIP%XTABGW_H)
!
ELSE
!
  ALLOCATE(TTRIP%XTABGW_F (0,0,0))      
  ALLOCATE(TTRIP%XTABGW_H (0,0,0))      
!
ENDIF
!
IF (LFLOOD.OR.CGROUNDW == 'DIF') THEN
!
  YVAR ='RIVDEPTH'     
  CALL READ_TRIP(NLISTING,YFILE_PARAM,YVAR,TTRIP%XHC_BED)  
!
ENDIF  
!
! * Calculate floodplains parameters
!
IF(LFLOOD)THEN
!
  YVAR ='NFLOOD'
  CALL READ_TRIP(NLISTING,YFILE_PARAM,YVAR,TTRIP%XN_FLOOD)
!  
  TGRID%GMASK_FLD(:,:)=TGRID%GMASK_VEL(:,:)
  WHERE(TTRIP%XN_FLOOD(:,:)==XUNDEF)
        TGRID%GMASK_FLD(:,:)=.FALSE.
  ENDWHERE
!
  ALLOCATE(TTRIP%XTAB_F (KLON,KLAT,NDIMTAB))      
  ALLOCATE(TTRIP%XTAB_H (KLON,KLAT,NDIMTAB))      
  ALLOCATE(TTRIP%XTAB_VF(KLON,KLAT,NDIMTAB))      
!
  YVAR ='TABF'
  CALL READ_TRIP(NLISTING,YFILE_PARAM,YVAR,TTRIP%XTAB_F)
!
  YVAR ='TABH'
  CALL READ_TRIP(NLISTING,YFILE_PARAM,YVAR,TTRIP%XTAB_H)
!
  YVAR ='TABVF'
  CALL READ_TRIP(NLISTING,YFILE_PARAM,YVAR,TTRIP%XTAB_VF)
!
ELSE
!
  ALLOCATE(TTRIP%XTAB_F (0,0,0))      
  ALLOCATE(TTRIP%XTAB_H (0,0,0))      
  ALLOCATE(TTRIP%XTAB_VF(0,0,0))      
!
ENDIF
!
!-------------------------------------------------------------------------------
! * Read initial and historical variables
!-------------------------------------------------------------------------------
!
YVAR ='SURF_STO'
CALL READ_TRIP(NLISTING,YFILE_INIT,YVAR,TTRIP%XSURF_STO)
!
IF(CGROUNDW=='CST')THEN
!        
  YVAR ='GROUND_STO'
  CALL READ_TRIP(NLISTING,YFILE_INIT,YVAR,TTRIP%XGROUND_STO)
!  
ELSEIF(CGROUNDW=='DIF')THEN
!
  YVAR ='HGROUND'
  CALL READ_TRIP(NLISTING,YFILE_INIT,YVAR,TTRIP%XHGROUND)
!
ENDIF
!
IF(LFLOOD)THEN
!
  YVAR ='FLOOD_STO'
  CALL READ_TRIP(NLISTING,YFILE_INIT,YVAR,TTRIP%XFLOOD_STO)
!
  YVAR ='FFLOOD'
  CALL READ_TRIP(NLISTING,YFILE_INIT,YVAR,TTRIP%XFFLOOD)
!
  YVAR ='HFLOOD'
  CALL READ_TRIP(NLISTING,YFILE_INIT,YVAR,TTRIP%XHFLOOD)
!  
ENDIF
!
!-------------------------------------------------------------------------------
! * Initial Conditions
!-------------------------------------------------------------------------------
!
ALLOCATE(ZHSTREAM(KLON,KLAT))
!
IF(CVIT == 'VAR')THEN       
!
  WHERE(TTRIP%XWIDTH(:,:)>0.0)
        ZHSTREAM(:,:)=MAX(0.0,TTRIP%XSURF_STO(:,:)/(XRHOLW*TGRID%XLEN(:,:)*TTRIP%XWIDTH(:,:)))
  ELSEWHERE
        ZHSTREAM(:,:)=0.0
  ENDWHERE
!
ENDIF
!
! * Fraction, width, water depth of floodplains
! 
IWORK=0
!
IF(LFLOOD)THEN 
!        
  WHERE(TGRID%GMASK(:,:).AND.TGRID%GMASK_FLD(:,:))
        TTRIP%XFLOOD_LEN(:,:) = XRATMED*SQRT(TTRIP%XFFLOOD(:,:)*TGRID%XAREA(:,:))
  ELSEWHERE
        TTRIP%XFLOOD_LEN(:,:) = 0.0
        TTRIP%XWFLOOD   (:,:) = 0.0
        TTRIP%XFFLOOD   (:,:) = 0.0
        TTRIP%XHFLOOD   (:,:) = 0.0
        TTRIP%XFLOOD_STO(:,:) = 0.0
  ENDWHERE
!
!
  WHERE(TTRIP%XFFLOOD(:,:)>0.0.AND.TGRID%GMASK_FLD(:,:))
        TTRIP%XWFLOOD(:,:) = TGRID%XAREA(:,:) * TTRIP%XFFLOOD(:,:) / TTRIP%XFLOOD_LEN(:,:)
  ELSEWHERE
        TTRIP%XWFLOOD(:,:) = 0.0
  ENDWHERE
!  
  INI   =COUNT(TGRID%GMASK)
  IWORK =COUNT(TGRID%GMASK_FLD)
  IFLOOD=COUNT(TTRIP%XFFLOOD>0.0)
!
ENDIF
!
!-------------------------------------------------------------------------------
! * Initialize coupling variables
!-------------------------------------------------------------------------------
!
CALL INIT_TRIP_CPL_ESM(KLON,KLAT)
!
!-------------------------------------------------------------------------------
!
WRITE(NLISTING,*)'Coupling_time_step :           ',PTSTEP_RUN
WRITE(NLISTING,*)'Output_time_step :             ',PTSTEP_DIAG
WRITE(NLISTING,*)'TRIP_time_step :               ',XTSTEP
WRITE(NLISTING,*)''
WRITE(NLISTING,*)'Sequence max :                 ',TGRID%NSEQMAX
WRITE(NLISTING,*)''
WRITE(NLISTING,*)'MEANDERING RATIO FIXED TO      ',XRATMED
WRITE(NLISTING,*)'CELL LENGTH MIN, MAX (km):     ',MINVAL(TGRID%XLEN/1.E3, TGRID%GMASK),  &
                                                   MAXVAL(TGRID%XLEN/1.E3, TGRID%GMASK)                                 
WRITE(NLISTING,*)'CELL AREA MIN, MAX (km²):      ',MINVAL(TGRID%XAREA/1.E6,TGRID%GMASK),  &
                                                   MAXVAL(TGRID%XAREA/1.E6,TGRID%GMASK)  
WRITE(NLISTING,*)''
IF(CGROUNDW=='CST')THEN
  WRITE(NLISTING,*)''
  WRITE(NLISTING,*)'Ground transf. time MIN, MAX : ',MINVAL(TTRIP%XTAUG/XDAY,TGRID%GMASK_GW),  &
                                                     MAXVAL(TTRIP%XTAUG/XDAY,TGRID%GMASK_GW)  
ELSEIF(CGROUNDW=='DIF')THEN
  WRITE(NLISTING,*)''
  WRITE(NLISTING,*)'Ground transf. time MIN, MAX : ',MINVAL(TTRIP%XTAUG/XDAY,TGRID%GMASK_GW),  &
                                                     MAXVAL(TTRIP%XTAUG/XDAY,TGRID%GMASK_GW)  
  WRITE(NLISTING,*)'Effective porosity MIN, MAX  : ',MINVAL(TTRIP%XWEFF,     TGRID%GMASK_GW),  &
                                                     MAXVAL(TTRIP%XWEFF,     TGRID%GMASK_GW)  
  WRITE(NLISTING,*)'Transmissivity      MIN, MAX : ',MINVAL(TTRIP%XTRANS,    TGRID%GMASK_GW),  &
                                                     MAXVAL(TTRIP%XTRANS,    TGRID%GMASK_GW)
ENDIF
IF(CVIT == 'VAR')THEN
  WRITE(NLISTING,*)'WIDTH_RIVER MIN, MAX :         ',MINVAL(TTRIP%XWIDTH,   TGRID%GMASK_VEL),  &
                                                     MAXVAL(TTRIP%XWIDTH,   TGRID%GMASK_VEL)
  WRITE(NLISTING,*)'N MANNING COEF MIN, MAX :      ',MINVAL(TTRIP%XN,       TGRID%GMASK_VEL),  &
                                                     MAXVAL(TTRIP%XN,       TGRID%GMASK_VEL)
  WRITE(NLISTING,*)'RIVER SLOPE MIN, MAX :         ',MINVAL(TTRIP%XSLOPEBED,TGRID%GMASK_VEL),  &
                                                     MAXVAL(TTRIP%XSLOPEBED,TGRID%GMASK_VEL)  
  WRITE(NLISTING,*)''
  WRITE(NLISTING,*)'Initial river height         : ',MINVAL(ZHSTREAM,  TGRID%GMASK_VEL),  &
                                                     MAXVAL(ZHSTREAM,  TGRID%GMASK_VEL)  
ENDIF
WRITE(NLISTING,*)''
WRITE(NLISTING,*)'Initial river storage        : ',MINVAL(TTRIP%XSURF_STO, TGRID%GMASK),  &
                                                   MAXVAL(TTRIP%XSURF_STO, TGRID%GMASK)  
IF(CGROUNDW=='CST')THEN                                                 
  WRITE(NLISTING,*)''
  WRITE(NLISTING,*)'Initial ground storage       : ',MINVAL(TTRIP%XGROUND_STO,TGRID%GMASK_GW),  &
                                                     MAXVAL(TTRIP%XGROUND_STO,TGRID%GMASK_GW)  
ELSEIF(CGROUNDW=='DIF')THEN
  WRITE(NLISTING,*)''
  WRITE(NLISTING,*)'Initial groundwater elevation : ',MINVAL(TTRIP%XHGROUND,TGRID%GMASK_GW),  &
                                                      MAXVAL(TTRIP%XHGROUND,TGRID%GMASK_GW)  
  WHERE(TGRID%GMASK_GW(:,:))
        ZHSTREAM(:,:)=TTRIP%XTOPO_RIV(:,:)-TTRIP%XHC_BED(:,:)+ZHSTREAM(:,:)
  ENDWHERE
  WRITE(NLISTING,*)'Initial river elevation       : ',MINVAL(ZHSTREAM,TGRID%GMASK_GW),  &
                                                      MAXVAL(ZHSTREAM,TGRID%GMASK_GW)
ENDIF
WRITE(NLISTING,*)''

IF(LFLOOD)THEN 
  WRITE(NLISTING,*)'N FLOOD FIXED TO               ',MINVAL(TTRIP%XN_FLOOD,TGRID%GMASK_FLD)
  WRITE(NLISTING,*)'RIVER DEPTH MIN, MAX :         ',MINVAL(TTRIP%XHC_BED, TGRID%GMASK_FLD),  &
                                                     MAXVAL(TTRIP%XHC_BED, TGRID%GMASK_FLD)  
  WRITE(NLISTING,*)''
  WRITE(NLISTING,*)'Number of potential flood cell : ',IWORK,'on',INI
  WRITE(NLISTING,*)'          %                      ',100.0*(FLOAT(IWORK)/FLOAT(INI))
  WRITE(NLISTING,*)'Number of actual flood cell :    ',IFLOOD,'on',IWORK
  WRITE(NLISTING,*)'          %                      ',100.0*(FLOAT(IFLOOD)/FLOAT(IWORK))
  WRITE(NLISTING,*)'% of flooded area in the domain :',100.0*SUM(TTRIP%XFFLOOD*TGRID%XAREA)/SUM(TGRID%XAREA)
  WRITE(NLISTING,*)''
  WRITE(NLISTING,*)'Initial flood depth (m) :      ',MINVAL(TTRIP%XHFLOOD,                 TGRID%GMASK_FLD), &
                                                     MAXVAL(TTRIP%XHFLOOD,                 TGRID%GMASK_FLD)  
  WRITE(NLISTING,*)'Initial flood fraction :       ',MINVAL(TTRIP%XFFLOOD,                 TGRID%GMASK_FLD), &
                                                     MAXVAL(TTRIP%XFFLOOD,                 TGRID%GMASK_FLD)  
  WRITE(NLISTING,*)'Initial flood volume m3/1E9 :  ',MINVAL(TTRIP%XFLOOD_STO/(XRHOLW*1.E9),TGRID%GMASK_FLD), &
                                                     MAXVAL(TTRIP%XFLOOD_STO/(XRHOLW*1.E9),TGRID%GMASK_FLD)   
  WRITE(NLISTING,*)'Initial flood length (km):     ',MINVAL(TTRIP%XFLOOD_LEN/1.E3,         TGRID%GMASK_FLD), &
                                                     MAXVAL(TTRIP%XFLOOD_LEN/1.E3,         TGRID%GMASK_FLD)                           
  WRITE(NLISTING,*)'Initial flood WIDTH (km) :     ',MINVAL(TTRIP%XWFLOOD/1.E3,            TGRID%GMASK_FLD), &
                                                     MAXVAL(TTRIP%XWFLOOD/1.E3,            TGRID%GMASK_FLD)      
  WRITE(NLISTING,*)''
ENDIF
!
DEALLOCATE(ZHSTREAM)
!
!-------------------------------------------------------------------------------
! * Create high frequency diag file
!-------------------------------------------------------------------------------
!
YFILE  = YDIAG  
YTITLE = 'TRIP high frequency outputs'  
CALL OUTPUT_DATE(YUNITTIME,XTIME_DIAG)
CALL INIT_TRIP_DIAG(NLISTING,YFILE,KLON,KLAT,YTITLE,YUNITTIME,.TRUE.) 
!
!-------------------------------------------------------------------------------
! * Create run mean diag file
!-------------------------------------------------------------------------------
!
YFILE     = YRUN
YTITLE    = 'TRIP run mean outputs'
WRITE(YTIME,'(i4.4,i2.2)') KYEAR, KMONTH
YUNITTIME = 'months since '//YTIME(1:4)//'-'//YTIME(5:LEN_TRIM(YTIME))//'-15'
CALL INIT_TRIP_DIAG(NLISTING,YFILE,KLON,KLAT,YTITLE,YUNITTIME,.TRUE.)
!
!
!-------------------------------------------------------------------------------
! * Create restart file
!-------------------------------------------------------------------------------
!
IF(ORESTART)THEN        
  YTITLE   ='TRIP restart variables'
  YUNITTIME='-'
  CALL INIT_RESTART_TRIP(NLISTING,YFILE_RESTART,KLON,KLAT,YTITLE,YUNITTIME,.FALSE.)
ENDIF
IF (LHOOK) CALL DR_HOOK('INIT_TRIP',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
! * CONTAINS
!-------------------------------------------------------------------------------
!
CONTAINS
!
SUBROUTINE OUTPUT_DATE(HTIMEUNIT,PTIME_DIAG)
!
IMPLICIT NONE
!
CHARACTER(LEN=*), INTENT(OUT) :: HTIMEUNIT
REAL,             INTENT(OUT) :: PTIME_DIAG
!
INTEGER, DIMENSION(3)         :: ITIME, IDATE
INTEGER                       :: INDAYS
INTEGER                       :: ISEC
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('INIT_TRIP:OUTPUT_DATE',0,ZHOOK_HANDLE)
!
INDAYS = FLOOR(PTSTEP_DIAG/86400.)
ISEC=MAX(0,NINT(PTIME-(PTSTEP_DIAG-INDAYS*86400)))
ITIME(1)=FLOOR(ISEC/3600.)
ITIME(2)=FLOOR((ISEC-ITIME(1)*3600)/60.)
ITIME(3)=ISEC-ITIME(1)*3600-ITIME(2)*60 
!
IF (PTSTEP_DIAG == FLOOR(PTSTEP_DIAG/86400.)*86400) THEN 
  HTIMEUNIT ='days since '
  PTIME_DIAG=86400.
ELSEIF (PTSTEP_DIAG == FLOOR(PTSTEP_DIAG/3600.)*3600) THEN
  HTIMEUNIT ='hours since '
  PTIME_DIAG=3600.
ELSEIF (PTSTEP_DIAG == FLOOR(PTSTEP_DIAG/60.)*60) THEN
  HTIMEUNIT ='minutes since '
  PTIME_DIAG=60.
ELSE
  HTIMEUNIT ='seconds since '
  PTIME_DIAG=1.
ENDIF
!
IDATE(1) = KYEAR
IDATE(2) = KMONTH
IDATE(3) = KDAY
!
CALL WRITE_TIME(IDATE(1),1,"-",HTIMEUNIT)
CALL WRITE_TIME(IDATE(2),0,"-",HTIMEUNIT)
CALL WRITE_TIME(IDATE(3),0,"",HTIMEUNIT)
CALL WRITE_TIME(ITIME(1),1,":",HTIMEUNIT)
CALL WRITE_TIME(ITIME(2),0,":",HTIMEUNIT)
CALL WRITE_TIME(ITIME(3),0,"",HTIMEUNIT)
!
IF (LHOOK) CALL DR_HOOK('INIT_TRIP:OUTPUT_DATE',1,ZHOOK_HANDLE)
!
END SUBROUTINE OUTPUT_DATE
!
SUBROUTINE WRITE_TIME(ITIME,ISPACE,HSEP,HTDATE)
!
INTEGER, INTENT(IN)             :: ITIME
INTEGER, INTENT(IN)             :: ISPACE
 CHARACTER(LEN=*), INTENT(IN)    :: HSEP
 CHARACTER(LEN=*), INTENT(INOUT) :: HTDATE
 CHARACTER(LEN=10)               :: YPAS
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!      
IF (LHOOK) CALL DR_HOOK('INIT_TRIP:WRITE_TIME',0,ZHOOK_HANDLE)
IF (ITIME.LT.10) THEN
  WRITE(YPAS,'(i1)') ITIME
  IF (ISPACE==1) THEN
    HTDATE=trim(HTDATE)//" 0"//trim(YPAS)//HSEP
  ELSE
    HTDATE=trim(HTDATE)//"0"//trim(YPAS)//HSEP
  ENDIF
ELSE
  IF (ITIME.LT.100) THEN
    WRITE(YPAS,'(i2)') ITIME
  ELSE
    WRITE(YPAS,'(i4)') ITIME
  ENDIF
  IF (ISPACE==1) THEN
    HTDATE=trim(HTDATE)//" "//trim(YPAS)//HSEP
  ELSE
    HTDATE=trim(HTDATE)//trim(YPAS)//HSEP
  ENDIF  
ENDIF
IF (LHOOK) CALL DR_HOOK('INIT_TRIP:WRITE_TIME',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_TIME
!
!-------------------------------------------------------------------------------
! * END
!-------------------------------------------------------------------------------
END SUBROUTINE INIT_TRIP





