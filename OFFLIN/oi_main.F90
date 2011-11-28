PROGRAM OI_MAIN
! ------------------------------------------------------------------------------------------
!  *****************************************************************************************
!
!  Program to perform within SURFEX 
!  a soil analysis for water content and temperature 
!  using the Meteo-France optimum interpolation technique of Giard and Bazile (2000)
!
!  Derived from CANARI subroutines externalized by Lora Taseva (Dec. 2007)
!
!  Author : Jean-Francois Mahfouf (01/2008)
!
!  Modifications : 
!   (05/2008)  : The I/O of this version follow the newly available LFI format in SURFEX  
!   (01/2009)  : Read directly atmospheric FA files using XRD library instead of using "edf"
!   (06/2009)  : Modifications to allow the assimilation of ASCAT superficial soil moisture
!   (09/2010)  : More parameters to goto_surfex
!   (03/2011)  : Initialization of ZEVAPTR (F.Bouyssel)
! ******************************************************************************************
! ------------------------------------------------------------------------------------------
 USE MODD_TYPE_DATE_SURF
 USE MODI_INI_CSTS
 USE MODI_INI_ASSIM
 USE MODI_READ_SURF
 USE MODI_WRITE_SURF
 USE MODD_OL_FILEID
 USE MODI_INIT_IO_SURF_n
 USE MODI_END_IO_SURF_n
 USE MODI_OPEN_NAMELIST
 USE MODI_CLOSE_NAMELIST
 USE MODE_POS_SURF
  USE MODN_IO_OFFLINE, ONLY : CSURF_FILETYPE
 USE MODD_IO_SURF_LFI,ONLY : CFILEIN_LFI, CFILEOUT_LFI
 USE MODD_IO_SURF_FA, ONLY : CFILEIN_FA, CDNOMC
 USE MODD_CSTS,       ONLY : XDAY, XPI, XG, XRHOLW, XLVTT, NDAYSEC
 USE MODD_SURF_PAR,   ONLY : XUNDEF
 USE MODI_CONVERT_COVER_FRAC
 USE MODI_GET_SIZE_FULL_n
 USE MODI_READ_COVER_n
 USE MODD_SURF_ATM_n,      ONLY :  CSEA,      CWATER,      CTOWN,      CNATURE, &
                                     XSEA,      XWATER,      XTOWN,      XNATURE,&
                                     NSIZE_SEA, NSIZE_WATER, NSIZE_TOWN, NSIZE_NATURE, &
                                     NR_SEA,    NR_WATER,    NR_TOWN,&
                                     NR_NATURE, XCOVER, NDIM_FULL, NSIZE_FULL, &
                                     NDIM_NATURE, NDIM_SEA, NDIM_WATER, NDIM_TOWN  
 USE MODD_ASSIM
 USE MODD_IO_SURF_FA,      ONLY :  NDGUX,  NDLUX,  PERPK,  PELON0, PELAT0, &
                                     PEDELX, PEDELY, PELON1, PELAT1, PEBETA  
!
! -------------------------------------------------------------------------------------
!
!
 USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
 USE PARKIND1  ,ONLY : JPRB
!
 USE MODI_GET_1D_MASK
 USE MODI_GOTO_SURFEX
 USE MODI_INI_DATA_COVER
 USE MODI_IO_BUFF_CLEAN_n
 USE MODI_OI_BC_SOIL_MOISTURE
 USE MODI_OI_CACSTS
 USE MODI_OI_HOR_EXTRAPOL_SURF
 USE MODI_OI_LATLON_CONF_PROJ
 USE MODI_TRANS_CHAINE
 USE MODI_READ_ALL_NAMELISTS
 USE MODI_ALLOC_SURFEX
 USE MODI_DEALLOC_SURFEX
 USE MODI_FLAG_UPDATE
!
 IMPLICIT NONE
 INTEGER                                  :: IDAT
 CHARACTER(LEN=28)                        :: YNAMELIST = 'OPTIONS.nam                 '
!
!    Declarations of local variables
!
 CHARACTER(LEN=3), PARAMETER              :: YINIT     = 'ALL'
 CHARACTER(LEN=6)                         :: YPROGRAM  = 'LFI   '
 CHARACTER(LEN=6)                         :: YPROGRAM2 = 'FA    '
 CHARACTER(LEN=2)                         :: CMONTH
 INTEGER                                  :: IYEAR                      ! current year (UTC)
 INTEGER                                  :: IMONTH                     ! current month (UTC)
 INTEGER                                  :: IDAY                       ! current day (UTC)
 INTEGER                                  :: IHOUR                      ! current day (UTC)
 INTEGER                                  :: NSSSSS                     ! current time since start of the run (s)
 INTEGER                                  :: IRESP                      ! return code
 TYPE (DATE_TIME)                         :: TTIME                      ! Current date and time  
!
! Arrays for soil OI analysis
!
 REAL, DIMENSION (:), ALLOCATABLE         :: PWS, PWP, PTS, PTP, PTL, PSST, PSNS, PLAI, PVEG, PRSMIN, PD2, PSAB, PARG,  &
                                               PLAT, PLON, PTSUN, PZENITH, PAZIMSOL, PTCLS, PHCLS, PRAIN, PSNOW,          &
                                               PWIND, PSWD, PSWS, PUCLS, PVCLS, PEVAP, PEVAPTR, PT2M_O, PHU2M_O,          &
                                               PTS_O, ZT2INC, ZH2INC, ZWS, ZWP, ZTL, ZTS, ZTP, ZTCLS, ZHCLS, ZUCLS,       &
                                               ZVCLS, PSSTC, PWPINC1, PWPINC2, PWPINC3, PT2MBIAS, PH2MBIAS, PRRCN, PRRCL, &
                                               PRRSN, PRRSL, PATMNEB, PITM, PALBF, PEMISF, PZ0F, PIVEG, PZ0H, PTSC,       &
                                               PTPC, PWSC, PWPC, PSNC, ZEVAP, ZEVAPTR, PGELAT, PGELAM, PGEMU,             &
                                               ZWSINC, ZWPINC,   ZTSINC, ZTPINC, ZTLINC, ZSNINC, ZSNS, ZPX, ZPY,            &
                                               PSM_O,  PSIG_SMO, PLSM_O, PWS_O,  ZWGINC , &
                                               PLST, PTRD3, ZSST, ZLST, ZALT  
!
 INTEGER                                  :: I,J
 CHARACTER(LEN=10)                        :: YVAR    ! Name of the prognostic variable (in LFI file)
 CHARACTER(LEN=100)                       :: YPREFIX ! Prefix of the prognostic variable  (in LFI file)
 INTEGER                                  :: ILUOUT  ! ascii output unit number
 INTEGER                                  :: ILUNAM  ! namelist unit number
 INTEGER                                  :: INOBS   ! number of observations
 LOGICAL                                  :: GFOUND
 LOGICAL                                  :: LALADSURF
 REAL                                     :: TPRT
!
 REAL                                     :: PLAT0,PLON0,PRPK,PLATOR,PLONOR,DELX,DELY,PBETA,ZTHRES
!
 LOGICAL, DIMENSION(:), ALLOCATABLE       :: OINTERP_LST, OINTERP_SST
 REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
 NAMELIST/NAM_NACVEG/NECHGU,RCLIMCA,RCLISST,SIGH2MO,SIGT2MO,SIGWGO,SIGWGB,SIGW2B, &
                       LOBSWG,LOBS2M,LIMVEG,SPRECIP2,RTHR_QC,SIGWGO_MAX,RSCAL_JAC, &
                       LPRINT,LAROME  
! ----------------------------------------------------------------------------------
!
 IF (LHOOK) CALL DR_HOOK('OI_MAIN',0,ZHOOK_HANDLE)
 CALL ALLOC_SURFEX(1)
 PRINT *,'--------------------------------------------------------------------------'
 PRINT *,'|                                                                        |'
 PRINT *,'|                             ENTER OI_ASSIM                             |'
 PRINT *,'|                                                                        |'
 PRINT *,'--------------------------------------------------------------------------'
!
 CALL READ_ALL_NAMELISTS(CSURF_FILETYPE,'ALL',.FALSE.)
!
 CALL GOTO_SURFEX(1,.TRUE.)
!
!   Initializations
!
 LALADSURF = .TRUE.
!
 CALL INI_CSTS
 CALL INI_ASSIM
!
!   Read namelist (to modify default values in NACVEG set by INI_ASSIM)
!
 ILUOUT = 111
!
 CALL OPEN_NAMELIST(CSURF_FILETYPE,ILUNAM,YNAMELIST)
 CALL POSNAM(ILUNAM,'NAM_NACVEG',GFOUND,ILUOUT)
 IF (GFOUND) READ (UNIT=ILUNAM,NML=NAM_NACVEG)
 CALL CLOSE_NAMELIST(CSURF_FILETYPE,ILUNAM)
!
!   Update some constants dependant from NACVEG
!
!  scaling of soil moisture increments when assimilation window is different
!  from 6 hours
 RSCALDW = REAL(NECHGU)/6.0
!  half assimilation window in sec
  ITRAD   = NECHGU*1800
!
 CALL INI_DATA_COVER
!
!   File handling definition
!
 CFILEIN_LFI = 'PREP'        ! input PREP file (surface fields)
!
!   Read grid dimension for allocation
!
 CALL INIT_IO_SURF_n(YPROGRAM,'FULL  ','SURF  ','READ ')
!
!   Find current time
!
 CALL READ_SURF(YPROGRAM,'DTCUR',TTIME,IRESP)
!
!   Time initializations 
!
 IYEAR  = TTIME%TDATE%YEAR
 IMONTH = TTIME%TDATE%MONTH
 IDAY   = TTIME%TDATE%DAY
 NSSSSS = TTIME%TIME
 IF (NSSSSS > NDAYSEC) NSSSSS = NSSSSS - NDAYSEC
 CALL TRANS_CHAINE(CMONTH,IMONTH,2)
!
!   Reading grid characteristics to perform nature mask
!
 CALL READ_SURF(YPROGRAM,'SEA   ',CSEA   ,IRESP)
 CALL READ_SURF(YPROGRAM,'WATER ',CWATER ,IRESP)
 CALL READ_SURF(YPROGRAM,'NATURE',CNATURE,IRESP)
 CALL READ_SURF(YPROGRAM,'TOWN  ',CTOWN  ,IRESP)
!
 CALL READ_SURF(YPROGRAM,'DIM_FULL  ',NDIM_FULL,  IRESP)
 CALL READ_SURF(YPROGRAM,'DIM_SEA   ',NDIM_SEA,   IRESP)
 CALL READ_SURF(YPROGRAM,'DIM_NATURE',NDIM_NATURE,IRESP)
 CALL READ_SURF(YPROGRAM,'DIM_WATER ',NDIM_WATER, IRESP)
 CALL READ_SURF(YPROGRAM,'DIM_TOWN  ',NDIM_TOWN,  IRESP)
!
!   Get total dimension of domain (excluding extension zone)
!
 CALL GET_SIZE_FULL_n(YPROGRAM,NDIM_FULL,NSIZE_FULL)
 CALL READ_COVER_n(YPROGRAM)
!
!   Perform masks (only nature used)
!
 ALLOCATE(XSEA   (NDIM_FULL))
 ALLOCATE(XNATURE(NDIM_FULL))
 ALLOCATE(XWATER (NDIM_FULL))
 ALLOCATE(XTOWN  (NDIM_FULL))
!
 CALL CONVERT_COVER_FRAC(XCOVER,XSEA,XNATURE,XTOWN,XWATER)
!
 NSIZE_NATURE = COUNT(XNATURE(:) > 0.0)
 NSIZE_TOWN   = COUNT(XTOWN(:)   > 0.0)
 NSIZE_WATER  = COUNT(XWATER(:)  > 0.0)
 NSIZE_SEA    = COUNT(XSEA(:)    > 0.0)
!
 ALLOCATE(NR_NATURE (NSIZE_NATURE))
 ALLOCATE(NR_TOWN   (NSIZE_TOWN  ))
 ALLOCATE(NR_WATER  (NSIZE_WATER ))
 ALLOCATE(NR_SEA    (NSIZE_SEA   ))
!
 CALL GET_1D_MASK( NSIZE_SEA,    NDIM_FULL, XSEA   , NR_SEA   )
 CALL GET_1D_MASK( NSIZE_WATER,  NDIM_FULL, XWATER , NR_WATER )
 CALL GET_1D_MASK( NSIZE_TOWN,   NDIM_FULL, XTOWN  , NR_TOWN  )
 CALL GET_1D_MASK( NSIZE_NATURE, NDIM_FULL, XNATURE, NR_NATURE)
!
! Allocate arrays 
! 
 ALLOCATE (PWS(NDIM_FULL))
 ALLOCATE (PWP(NDIM_FULL))
 ALLOCATE (PTS(NDIM_FULL))
 ALLOCATE (PTP(NDIM_FULL))
 ALLOCATE (PTL(NDIM_FULL))
 ALLOCATE (PSST(NDIM_FULL))
 ALLOCATE (PSNS(NDIM_FULL))
 ALLOCATE (PLAI(NDIM_FULL))
 ALLOCATE (PVEG(NDIM_FULL))
 ALLOCATE (PRSMIN(NDIM_FULL))
 ALLOCATE (PD2(NDIM_FULL))
 ALLOCATE (PSAB(NDIM_FULL))
 ALLOCATE (PARG(NDIM_FULL))
 ALLOCATE (PTCLS(NDIM_FULL))
 ALLOCATE (PHCLS(NDIM_FULL))
 ALLOCATE (PUCLS(NDIM_FULL))
 ALLOCATE (PVCLS(NDIM_FULL))
 ALLOCATE (PEVAP(NDIM_FULL))
 ALLOCATE (PLST(NDIM_FULL))
 ALLOCATE (PTRD3(NDIM_FULL))
!
 ALLOCATE (OINTERP_LST(NDIM_FULL))
 ALLOCATE (OINTERP_SST(NDIM_FULL))
 ALLOCATE (ZLST(NDIM_FULL))
 ALLOCATE (ZSST(NDIM_FULL))
 ALLOCATE (ZALT(NDIM_FULL))
!
!  Read prognostic variables
!
 CALL READ_SURF(YPROGRAM,'WG1',       PWS,   IRESP)
 CALL READ_SURF(YPROGRAM,'WG2',       PWP,   IRESP)
 CALL READ_SURF(YPROGRAM,'TG1',       PTS,   IRESP)
 CALL READ_SURF(YPROGRAM,'TG2',       PTP,   IRESP)
 CALL READ_SURF(YPROGRAM,'WGI2',      PTL,   IRESP)
 CALL READ_SURF(YPROGRAM,'WSNOW_VEG1',PSNS,  IRESP)
 CALL READ_SURF(YPROGRAM,'SST',       PSST,  IRESP)
 CALL READ_SURF(YPROGRAM,'TS_WATER',  PLST,  IRESP)
 IF (NSIZE_TOWN > 0 .AND. LAROME) THEN         
   CALL READ_SURF(YPROGRAM,'T_ROAD3',   PTRD3,  IRESP)
 ELSE
   PTRD3(:) = XUNDEF
 ENDIF
!
 CALL READ_SURF(YPROGRAM,'T2M',       PTCLS, IRESP)
 CALL READ_SURF(YPROGRAM,'HU2M',      PHCLS, IRESP)
 CALL READ_SURF(YPROGRAM,'ZON10M',    PUCLS, IRESP)
 CALL READ_SURF(YPROGRAM,'MER10M',    PVCLS, IRESP)
!
! Read constant surface fields
!
 CALL READ_SURF(YPROGRAM,'RSMIN',     PRSMIN,IRESP)
 CALL READ_SURF(YPROGRAM,'DG2',       PD2,   IRESP)
 CALL READ_SURF(YPROGRAM,'SAND',      PSAB,  IRESP)
 CALL READ_SURF(YPROGRAM,'CLAY',      PARG,  IRESP)
 CALL READ_SURF(YPROGRAM,'LAI',       PLAI,  IRESP)
 CALL READ_SURF(YPROGRAM,'VEG',       PVEG,  IRESP)
 CALL READ_SURF(YPROGRAM,'ZS',        ZALT,  IRESP) 
!
! PRINT 
!
 IF (LPRINT) THEN
  J = NR_NATURE(1)
  PRINT *,'value in PREP file => WG1       ',PWS(j)
  PRINT *,'value in PREP file => WG2       ',PWP(j)
  PRINT *,'value in PREP file => TG1       ',PTS(j)
  PRINT *,'value in PREP file => TG2       ',PTP(j)
  PRINT *,'value in PREP file => WGI2      ',PTL(j)
  PRINT *,'value in PREP file => WSNOW_VEG1',PSNS(j)
  PRINT *,'value in PREP file => SST       ',PSST(j)
  PRINT *,'value in PREP file => LAI       ',PLAI(j)
  PRINT *,'value in PREP file => VEG       ',PVEG(j)
  PRINT *,'value in PREP file => RSMIN     ',PRSMIN(j)
  PRINT *,'value in PREP file => DATA_DG2  ',PD2(j)
  PRINT *,'value in PREP file => SAND      ',PSAB(j)
  PRINT *,'value in PREP file => CLAY      ',PARG(j)
  PRINT *,'value in PREP file => ZS        ',ZALT(j)  
 ENDIF
!
 CALL END_IO_SURF_n(YPROGRAM)
 CALL IO_BUFF_CLEAN_n
!
!  Interface (allocate arrays)
!
 ALLOCATE (PLAT(NDIM_FULL))
 ALLOCATE (PLON(NDIM_FULL))
 ALLOCATE (ZPX(NDIM_FULL))
 ALLOCATE (ZPY(NDIM_FULL))
 ALLOCATE (PEVAPTR(NDIM_FULL))
 ALLOCATE (ZWP(NDIM_FULL))
 ALLOCATE (ZWS(NDIM_FULL))
 ALLOCATE (ZTL(NDIM_FULL))
 ALLOCATE (ZTS(NDIM_FULL))
 ALLOCATE (ZTP(NDIM_FULL))
 ALLOCATE (ZSNS(NDIM_FULL))
 ALLOCATE (ZTCLS(NDIM_FULL))
 ALLOCATE (ZHCLS(NDIM_FULL))
 ALLOCATE (ZUCLS(NDIM_FULL))
 ALLOCATE (ZVCLS(NDIM_FULL))
 ALLOCATE (PSSTC(NDIM_FULL))
 ALLOCATE (PWPINC1(NDIM_FULL))
 ALLOCATE (PWPINC2(NDIM_FULL))
 ALLOCATE (PWPINC3(NDIM_FULL))
 ALLOCATE (PT2MBIAS(NDIM_FULL))
 ALLOCATE (PH2MBIAS(NDIM_FULL))
 ALLOCATE (PRRCN(NDIM_FULL))
 ALLOCATE (PRRCL(NDIM_FULL))
 ALLOCATE (PRRSN(NDIM_FULL))
 ALLOCATE (PRRSL(NDIM_FULL))
 ALLOCATE (PATMNEB(NDIM_FULL))
 ALLOCATE (PITM(NDIM_FULL))
 ALLOCATE (PALBF(NDIM_FULL))
 ALLOCATE (PEMISF(NDIM_FULL))
 ALLOCATE (PZ0F(NDIM_FULL))
 ALLOCATE (PIVEG(NDIM_FULL))
 ALLOCATE (PZ0H(NDIM_FULL))
 ALLOCATE (PTSC(NDIM_FULL))
 ALLOCATE (PTPC(NDIM_FULL))
 ALLOCATE (PWSC(NDIM_FULL))
 ALLOCATE (PWPC(NDIM_FULL)) 
 ALLOCATE (PSNC(NDIM_FULL))
 ALLOCATE (ZEVAP(NDIM_FULL))
 ALLOCATE (ZEVAPTR(NDIM_FULL))
 ALLOCATE (PGELAT(NDIM_FULL))
 ALLOCATE (PGELAM(NDIM_FULL))
 ALLOCATE (PGEMU(NDIM_FULL))
 ALLOCATE (PT2M_O(NDIM_FULL))
 ALLOCATE (PHU2M_O(NDIM_FULL))
 ALLOCATE (PTS_O(NDIM_FULL))
 ALLOCATE (PSM_O(NDIM_FULL))
 ALLOCATE (PSIG_SMO(NDIM_FULL))
 ALLOCATE (PLSM_O(NDIM_FULL))
 ALLOCATE (PWS_O(NDIM_FULL))
 ALLOCATE (ZWGINC(NDIM_FULL))
!
!  Read atmospheric forecast fields from FA files 
!
 CFILEIN_FA = 'FG_OI_MAIN'        ! input forecast
!
!  Open FA file (LAM version with extension zone)
!
 CALL INIT_IO_SURF_n(YPROGRAM2,'EXTZON','SURF  ','READ ')
!
!  Read model forecast quantities
!
 IF (LAROME) THEN
   CALL READ_SURF(YPROGRAM2,'SURFACCPLUIE',    PRRSL  ,IRESP)
   CALL READ_SURF(YPROGRAM2,'SURFACCNEIGE',    PRRSN  ,IRESP)
   CALL READ_SURF(YPROGRAM2,'SURFACCGRAUPEL',  PRRCN  ,IRESP)
   PRRCL(:) = 0.0
!   CALL READ_SURF(YPROGRAM2,'SURFIND.VEG.DOMI',PIVEG  ,IRESP) 
   PIVEG(:) = 0.0
 ELSE
   CALL READ_SURF(YPROGRAM2,'SURFPREC.EAU.CON',PRRCL  ,IRESP)
   CALL READ_SURF(YPROGRAM2,'SURFPREC.EAU.GEC',PRRSL  ,IRESP)
   CALL READ_SURF(YPROGRAM2,'SURFPREC.NEI.CON',PRRCN  ,IRESP)
   CALL READ_SURF(YPROGRAM2,'SURFPREC.NEI.GEC',PRRSN  ,IRESP)
   PIVEG(:) = 0.0
 ENDIF
 CALL READ_SURF(YPROGRAM2,'ATMONEBUL.BASSE ',PATMNEB,IRESP)
 CALL READ_SURF(YPROGRAM2,'SURFIND.TERREMER',PITM   ,IRESP) 
 CALL READ_SURF(YPROGRAM2,'SURFFLU.LAT.MEVA',PEVAP  ,IRESP) ! accumulated fluxes (not available in LFI)
 IF (.NOT.LALADSURF) THEN
   CALL READ_SURF(YPROGRAM2,'SURFXEVAPOTRANSP',PEVAPTR,IRESP) ! not in ALADIN SURFEX
 ELSE
   PEVAPTR(:) = 0.0
 ENDIF
!
!  Close FA file
! 
 CALL END_IO_SURF_n(YPROGRAM2)
 CALL IO_BUFF_CLEAN_n
 PRINT *,'READ FG_OI_MAIN OK'
!
!  Define FA file name for CANARI analysis
!
 CFILEIN_FA = 'CANARI'        ! input CANARI analysis
!
!  Open FA file 
!
 CALL INIT_IO_SURF_n(YPROGRAM2,'EXTZON','SURF  ','READ ')
!
!  Read CANARI analysis
!
 CALL READ_SURF(YPROGRAM2,'CLSTEMPERATURE  ',PT2M_O ,IRESP)
 CALL READ_SURF(YPROGRAM2,'CLSHUMI.RELATIVE',PHU2M_O,IRESP)
 CALL READ_SURF(YPROGRAM2,'SURFTEMPERATURE ',PTS_O  ,IRESP)
!
!  Close CANARI file
!
 CALL END_IO_SURF_n(YPROGRAM2)
 CALL IO_BUFF_CLEAN_n
 PRINT *,'READ CANARI OK'
!
!  Read ASCAT SM observations (in percent)
!
 INOBS = 0
 IF (LOBSWG) THEN
   OPEN(UNIT=111,FILE='ASCAT_SM.DAT')
   DO I=1,NDIM_FULL
     READ(111,*) PSM_O(I),PSIG_SMO(I),PLSM_O(I)
     IF (PLSM_O(I) < 1.0)          PSM_O(I) = 999.0 ! data rejection if not on land
     IF (PSIG_SMO(I) > SIGWGO_MAX) PSM_O(I) = 999.0 ! data rejection of error too large
     IF (PSM_O(I) /= 999.0) INOBS = INOBS + 1
   ENDDO
   CLOSE(UNIT=111)
   PRINT *,'READ ASCAT SM OK'
 ELSE
   PSM_O(:)    = 999.0
   PSIG_SMO(:) = 999.0
   PLSM_O(:)   = 0.0
 ENDIF
 PRINT *,' NUMBER OF ASCAT OBSERVATIONS AFTER INITIAL CHECKS  :: ',INOBS
 INOBS = 0
!
! Perform bias correction of SM observations
!
 CALL OI_BC_SOIL_MOISTURE(NDIM_FULL,PSM_O,PSAB,PWS_O)
!
!  Define FA file name for surface climatology
!
 CFILEIN_FA = 'clim_isba'               ! input climatology
 CDNOMC     = 'climat'                  ! new frame name
!
!  Open FA file 
!
 CALL INIT_IO_SURF_n(YPROGRAM2,'EXTZON','SURF  ','READ ')
!
!  Read climatology file (snow water equivalent)
!
 CALL READ_SURF(YPROGRAM2,'SURFRESERV.NEIGE',PSNC  ,IRESP)
!
!  Close climatology file
!
 CALL END_IO_SURF_n(YPROGRAM2)
 CALL IO_BUFF_CLEAN_n
 PRINT *,'READ CLIMATOLOGY OK'
!
 PLAT0  = PELAT0 
 PLON0  = PELON0 
 PLATOR = PELAT1 
 PLONOR = PELON1 
 PRPK   = PERPK  
 PBETA  = PEBETA 
 DELX   = PEDELX 
 DELY   = PEDELY 
!
 IF (PLONOR > 180.0) PLONOR = PLONOR - 360.0
 IF (PLON0  > 180.0) PLON0  = PLON0  - 360.0
!
 DO J=1,NDGUX
  DO I=1,NDLUX
    ZPX((J-1)*NDLUX + I) = DELX*REAL(I-1)
    ZPY((J-1)*NDLUX + I) = DELY*REAL(J-1)
  ENDDO
 ENDDO
!
 CALL OI_LATLON_CONF_PROJ(NDIM_FULL,PLAT0,PLON0,PRPK,PBETA,PLATOR,PLONOR,ZPX,ZPY,PLAT,PLON)
!
!  Allocate arrays to produce analysis increments  
!
 ALLOCATE (ZT2INC(NDIM_FULL))
 ALLOCATE (ZH2INC(NDIM_FULL))
 ALLOCATE (ZWSINC(NDIM_FULL))
 ALLOCATE (ZWPINC(NDIM_FULL))
 ALLOCATE (ZTLINC(NDIM_FULL))
 ALLOCATE (ZTSINC(NDIM_FULL))
 ALLOCATE (ZTPINC(NDIM_FULL))
 ALLOCATE (ZSNINC(NDIM_FULL))
!
! Screen-level innovations
!
 ZT2INC(:) = PT2M_O(:) - PTCLS(:)
 ZH2INC(:) = PHU2M_O(:) - PHCLS(:)
!
! Threshold for background check
!
 ZTHRES=RTHR_QC*SQRT(SIGWGO**2 + SIGWGB**2)
!
! Superficial soil moisture innovations in (m3/m3)
!
 DO I=1,NDIM_FULL
   IF (PWS_O(I) /= 999.0) THEN
     ZWGINC(I) = PWS_O(I) - PWS(I)
     IF (ABS(ZWGINC(I)) > ZTHRES) THEN 
       ZWGINC(I) = 0.0 ! background check
     ELSE
       INOBS = INOBS + 1
     ENDIF
   ELSE
     ZWGINC(I) = 0.0
   ENDIF
 ENDDO
 PRINT *,' NUMBER OF ASCAT OBSERVATIONS AFTER BACKGROUND CHECK  :: ',INOBS
!
 PRINT *,'           '
 PRINT *,'Mean T2m increments  ',SUM(ZT2INC)/SIZE(ZT2INC)
 PRINT *,'Mean HU2m increments ',SUM(ZH2INC)/SIZE(ZH2INC)
 PRINT *,'           '
!
! Interface (define arrays and perform unit conversions)
!
 PARG(:)     = PARG(:)*100.0
 PSAB(:)     = PSAB(:)*100.0
!
 ZWS(:) = XUNDEF
 ZWP(:) = XUNDEF
 ZTL(:) = XUNDEF 
 WHERE (PWS(:)/=XUNDEF) 
   ZWS(:)      = PWS(:)*RD1*XRHOLW     ! conversion of m3/m3 -> mm
   ZWP(:)      = PWP(:)*PD2(:)*XRHOLW  ! conversion of m3/m3 -> mm
   ZTL(:)      = PTL(:)*PD2(:)*XRHOLW  ! conversion of m3/m3 -> mm
 END WHERE
 ZTCLS(:)    = PTCLS(:)
 ZHCLS(:)    = PHCLS(:)
 ZUCLS(:)    = PUCLS(:)
 ZVCLS(:)    = PVCLS(:)
 PSSTC(:)    = PTS_O(:)
 PWPINC1(:)  = XUNDEF
 PWPINC2(:)  = XUNDEF
 PWPINC3(:)  = XUNDEF
 PT2MBIAS(:) = XUNDEF
 PH2MBIAS(:) = XUNDEF
!
! Sea-ice surface properties
!
 PALBF(:)    = XUNDEF
 PEMISF(:)   = XUNDEF
 PZ0F(:)     = XUNDEF
 PZ0H(:)     = XUNDEF
!
! Climatological arrays set to missing values
!
 PSNC(:)     =  PSNS(:) ! need to read the snow climatology 
 PWSC(:)     =  XUNDEF
 PWPC(:)     =  XUNDEF
 PTSC(:)     =  XUNDEF
 PTPC(:)     =  XUNDEF
! 
 DO I=1,NDIM_FULL
   PGELAT(I)   = PLAT(I) 
   PGELAM(I)   = PLON(I) 
   PGEMU(I)    = SIN(PLAT(I)*XPI/180.)
 ENDDO
!
 ZEVAP(:)   =  (PEVAP(:)/XLVTT*XDAY)/(NECHGU*3600.) ! conversion W/m2 -> mm/day
 ZEVAPTR(:) =  PEVAPTR(:)*XDAY 
 ZSNS(:)    =  PSNS(:)
!
 DO I=1,NDIM_FULL
   ZTS(I) = PTS(I)
   ZTP(I) = PTP(I)
 ENDDO
!
 IDAT = IYEAR*10000. + IMONTH*100. + IDAY
!
!  Soil analysis based on optimal interpolation
!
 CALL OI_CACSTS(NDIM_FULL,ZT2INC,ZH2INC,ZWGINC,PWS_O,                  &
                  IDAT,NSSSSS,                                           &
                  ZTP,ZWP,ZTL,ZSNS,ZTS,ZWS,                              &
                  ZTCLS,ZHCLS,ZUCLS,ZVCLS,PSSTC,PWPINC1,PWPINC2,PWPINC3, &
                  PT2MBIAS,PH2MBIAS,                                     &
                  PRRCL,PRRSL,PRRCN,PRRSN,PATMNEB,ZEVAP,ZEVAPTR,         &
                  PITM,PVEG,PALBF,PEMISF,PZ0F,                           &
                  PIVEG,PARG,PD2,PSAB,PLAI,PRSMIN,PZ0H,                  &
                  PTSC,PTPC,PWSC,PWPC,PSNC,                              &
                  PGELAT,PGELAM,PGEMU)  
!
 PRINT *,'after OI_CACSTS'
!
!  Store increments
!
 ZWSINC(:) = 0.0
 ZWPINC(:) = 0.0
 ZTLINC(:) = 0.0
 ZSNINC(:) = 0.0
!
 WHERE (PWS(:)/=XUNDEF)
   ZWSINC(:) = ZWS(:) - PWS(:)*(RD1*XRHOLW)    
   ZWPINC(:) = ZWP(:) - PWP(:)*(PD2(:)*XRHOLW) 
   ZTLINC(:) = ZTL(:) - PTL(:)*(PD2(:)*XRHOLW) 
   ZSNINC(:) = ZSNS(:) - PSNS(:)
 END WHERE
!
!  Define soil moiture analyses over NATURE points
!
 WHERE (PWS(:)/=XUNDEF)
   PWS(:)  = ZWS(:)/(RD1*XRHOLW)
   PWP(:)  = ZWP(:)/(PD2(:)*XRHOLW)
   PTL(:)  = ZTL(:)/(PD2(:)*XRHOLW)
   PSNS(:) = ZSNS(:)
 END WHERE
!
!  Perform temperature analysis according to surface types
!
 OINTERP_LST(:) = .FALSE.
 OINTERP_SST(:) = .FALSE.
!
 ZTSINC(:) = 0.0
 ZTPINC(:) = 0.0
!
! a) Temperature analysis of NATURE points
!
 WHERE (PTS(:)/=XUNDEF)
   ZTSINC(:) = ZTS(:) - PTS(:)
   ZTPINC(:) = ZTP(:) - PTP(:)
   PTS(:)    = ZTS(:)
   PTP(:)    = ZTP(:)
 END WHERE
!
! b) Temperature analysis of SEA and LAKE points
!
 DO I=1,NDIM_FULL
   IF (PITM(I) < 0.5) THEN
     IF (PSST(I)/=XUNDEF) THEN
       ZTSINC(I) = PTS_O(I) - PSST(I)
       PSST(I) = PTS_O(I)   ! canari
     ELSEIF (PLST(I)/=XUNDEF) THEN
       PLST(I) = PTS_O(I)   ! canari
     ENDIF
   ELSE
     IF (PSST(I)/=XUNDEF) THEN
       PSST(I) = XUNDEF
       OINTERP_SST(I) = .TRUE.
     ELSEIF (PLST(I)/=XUNDEF) THEN
       PLST(I) = XUNDEF
       OINTERP_LST(I) = .TRUE.
     ENDIF
   ENDIF
 ENDDO
!
! c) Temperature analysis of TOWN points
!
 WHERE (PTRD3(:)/=XUNDEF)
  PTRD3(:) = PTRD3(:) + ZT2INC(:)/(2.0*XPI)
 END WHERE
!
! Search for the nearest grid point values for lake and sea points
! at locations where the water fraction is less than 50 % 
! and therefore no useful information is given from the SST analysis
! A standard temperature gradient is applied to account for the atitude differences
!
 ZLST(:) = PLST(:)
 CALL OI_HOR_EXTRAPOL_SURF(NDIM_FULL,PLAT,PLON,ZLST,PLAT,PLON,PLST,ZALT,OINTERP_LST)
 ZSST(:) = PSST(:)
 CALL OI_HOR_EXTRAPOL_SURF(NDIM_FULL,PLAT,PLON,ZSST,PLAT,PLON,PSST,ZALT,OINTERP_SST)
!
! PRINT values produced by OI_HO_EXTRAPOL_SURF
!
 IF (LPRINT) THEN
   DO I=1,NDIM_FULL
    IF (OINTERP_LST(I)) THEN
      PRINT *,'Lake surface temperature set to ',PLST(I),'from nearest neighbour at I=',I
    ENDIF
    IF (OINTERP_SST(I)) THEN
      PRINT *,'Sea surface temperature set to ',PSST(I),'from nearest neighbour at I=',I
    ENDIF
   ENDDO
 ENDIF
!
! PRINT statistics of the soil analysis
!
 PRINT *,'---------------------------------------------------------------'
 PRINT *,'Mean WS increments over NATURE ',SUM(ZWSINC,XNATURE > 0.)/NDIM_NATURE
 PRINT *,'Mean WP increments over NATURE ',SUM(ZWPINC,XNATURE > 0.)/NDIM_NATURE
 PRINT *,'Mean TS increments over NATURE ',SUM(ZTSINC,XNATURE > 0.)/NDIM_NATURE
 PRINT *,'Mean TP increments over NATURE ',SUM(ZTPINC,XNATURE > 0.)/NDIM_NATURE
 PRINT *,'Mean TL increments over NATURE ',SUM(ZTLINC,XNATURE > 0.)/NDIM_NATURE
 PRINT *,'Mean SN increments over NATURE ',SUM(ZSNINC,XNATURE > 0.)/NDIM_NATURE
 PRINT *,'---------------------------------------------------------------'
 PRINT *,'Mean WS increments over SEA    ',SUM(ZWSINC,XSEA > 0.)/NDIM_SEA
 PRINT *,'Mean WP increments over SEA    ',SUM(ZWPINC,XSEA > 0.)/NDIM_SEA
 PRINT *,'Mean TS increments over SEA    ',SUM(ZTSINC,XSEA > 0.)/NDIM_SEA
 PRINT *,'Mean TP increments over SEA    ',SUM(ZTPINC,XSEA > 0.)/NDIM_SEA
 PRINT *,'Mean TL increments over SEA    ',SUM(ZTLINC,XSEA > 0.)/NDIM_SEA
 PRINT *,'Mean SN increments over SEA    ',SUM(ZSNINC,XSEA > 0.)/NDIM_SEA
 PRINT *,'---------------------------------------------------------------'
!
!   Write analysis in LFI file PREP
!
 CFILEOUT_LFI='PREP'
 CALL FLAG_UPDATE(.FALSE.,.FALSE.,.TRUE.,.FALSE.,.FALSE.)
 CALL INIT_IO_SURF_n(YPROGRAM,'FULL  ','SURF  ','WRITE')
!
 YVAR='WG1'
 YPREFIX='X_Y_WG1 (m3/m3)                                   '
 CALL WRITE_SURF(YPROGRAM,YVAR,PWS,IRESP,HCOMMENT=YPREFIX)
 YVAR='WG2'
 YPREFIX='X_Y_WG2 (m3/m3)                                   '
 CALL WRITE_SURF(YPROGRAM,YVAR,PWP,IRESP,HCOMMENT=YPREFIX)
 YVAR='WGI2'
 YPREFIX='X_Y_WGI2 (m3/m3)                                  '
 CALL WRITE_SURF(YPROGRAM,YVAR,PTL,IRESP,HCOMMENT=YPREFIX)
 YVAR='TG1'
 YPREFIX='X_Y_TG1 (K)                                       '
 CALL WRITE_SURF(YPROGRAM,YVAR,PTS,IRESP,HCOMMENT=YPREFIX)
 YVAR='TG2'
 YPREFIX='X_Y_TG2 (K)                                       '
 CALL WRITE_SURF(YPROGRAM,YVAR,PTP,IRESP,HCOMMENT=YPREFIX)
 YVAR='SST'
 YPREFIX='X_Y_SST (K)                                       '
 CALL WRITE_SURF(YPROGRAM,YVAR,PSST,IRESP,HCOMMENT=YPREFIX)
 YVAR='TS_WATER'
 YPREFIX='X_Y_TS_WATER (K)                                  '
 CALL WRITE_SURF(YPROGRAM,YVAR,PLST,IRESP,HCOMMENT=YPREFIX)
 IF (NSIZE_TOWN > 0 .AND. LAROME) THEN 
   YVAR='T_ROAD3'
   YPREFIX='X_Y_T_ROAD3 (K)                                   '
   CALL WRITE_SURF(YPROGRAM,YVAR,PTRD3,IRESP,HCOMMENT=YPREFIX)
 ENDIF
 YVAR='WSNOW_VEG1'
 YPREFIX='X_Y_WSNOW_VEG1 (kg/m2)                            '
 CALL WRITE_SURF(YPROGRAM,YVAR,PSNS,IRESP,HCOMMENT=YPREFIX)
!
 CALL END_IO_SURF_n(YPROGRAM)
 CALL IO_BUFF_CLEAN_n
 PRINT *,'after write in PREP file'
!
CALL DEALLOC_SURFEX
! -------------------------------------------------------------------------------------
 IF (LHOOK) CALL DR_HOOK('OI_MAIN',1,ZHOOK_HANDLE)
 END PROGRAM OI_MAIN
