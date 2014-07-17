! *****************************************************************************************
PROGRAM SODA
!
! ------------------------------------------------------------------------------------------
!!
!!    SODA: SURFEX Offline Data Assimilation
!!
!!    PURPOSE
!!    -------
!!    Program to perform surface data within SURFEX 
!!
!!
!!    METHOD
!!    ------
!!    Different methods for different tiles
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    T. Aspelien                  met.no
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original         04/2012
!!
!! 03/2014 E. Martin change indices names in OMP module according to GMAP changes
!----------------------------------------------------------------------------
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO, WLOG_MPI, PREP_LOG_MPI, NPROC, NCOMM,   &
                            NINDEX, NSIZE_TASK, END_LOG_MPI, NSIZE
USE MODD_SURFEX_OMP, ONLY : NINDX2SFX, NWORK, XWORK, XWORK2, XWORK3, &
                            NWORK_FULL, XWORK_FULL, XWORK2_FULL, NWORK2, &
                            NWORK2_FULL
!
USE MODD_MASK, ONLY: NMASK_FULL
!
USE MODD_TYPE_DATE_SURF, ONLY : DATE_TIME
!
USE MODD_SURF_CONF, ONLY : CSOFTWARE
USE MODD_SURF_PAR,  ONLY : XUNDEF,NUNDEF
!
USE MODD_SURF_ATM_GRID_n, ONLY : XLON, XLAT
!
USE MODD_SURF_ATM_n, ONLY : NSIZE_NATURE, NSIZE_FULL
USE MODD_ISBA_n,     ONLY : NPATCH, XWG, XTG, XLAI, XBIOMASS, XRESP_BIOMASS
!
USE MODD_ASSIM, ONLY : LASSIM, LAROME, LALADSURF, CASSIM_ISBA, NVAR, XF, XF_PATCH, &
                       NOBSTYPE, XAT2M_ISBA, XAHU2M_ISBA, CVAR, COBS, NECHGU, XI,  &
                       NBOUTPUT, XLAI_PASS, XBIO_PASS, CBIO, NIVAR, NIPERT, XYO,   &
                       NOBS, LOBSFILE, NPRINTLEV
!
USE MODD_FORC_ATM,       ONLY : CSV, XDIR_ALB, XSCA_ALB, XEMIS, XTSRAD, XTSUN, XZS, &
                                XZREF, XUREF, XTA, XQA, XSV, XU, XV, XSW_BANDS,     &
                                XZENITH, XAZIM, XCO2, XRHOA 
!
#ifdef ARO 
USE MODD_IO_SURF_ARO,ONLY : NGPTOT, NGPTOT_CAP, NPROMA, NINDX1, NINDX2, NBLOCK, NKPROMA
#endif
!
#ifdef NC
USE MODD_IO_SURF_NC,   ONLY : CFILEIN_NC, CFILEIN_NC_SAVE, CFILEPGD_NC, CFILEOUT_NC, LDEF, &
                              CLUOUT_NC
#endif
#ifdef ASC
USE MODD_IO_SURF_ASC,  ONLY : CFILEIN, CFILEIN_SAVE, CFILEPGD, CFILEOUT
#endif
#ifdef FA
USE MODD_IO_SURF_FA,   ONLY : CFILEIN_FA, CFILEIN_FA_SAVE, CFILEPGD_FA, CDNOMC, CFILEOUT_FA, &
                              NUNIT_FA, IVERBFA, LFANOCOMPACT
#endif
#ifdef LFI
USE MODD_IO_SURF_LFI,    ONLY : CFILEIN_LFI, CFILEIN_LFI_SAVE, &
                                CFILEPGD_LFI, CFILE_LFI, CLUOUT_LFI, CFILEOUT_LFI 
#endif
!
USE MODN_IO_OFFLINE,     ONLY : NAM_IO_OFFLINE, CNAMELIST, CPGDFILE, CPREPFILE, CSURFFILE, &
                                CSURF_FILETYPE, CTIMESERIES_FILETYPE, LLAND_USE,           &
                                LDIAG_FA_NOCOMPACT, LOUT_TIMENAME
!
USE MODE_POS_SURF,  ONLY : POSNAM
!
USE MODI_ABOR1_SFX
USE MODI_ALLOC_SURFEX
USE MODI_DEALLOC_SURFEX
USE MODI_GOTO_SURFEX
USE MODI_GET_LUOUT
USE MODI_OPEN_NAMELIST
USE MODI_CLOSE_NAMELIST
USE MODI_READ_ALL_NAMELISTS
USE MODI_INIT_IO_SURF_n
USE MODI_END_IO_SURF_n
USE MODI_READ_SURF
USE MODI_IO_BUFF_CLEAN_n
USE MODI_GET_SIZE_FULL_n
USE MODI_INIT_SURF_ATM_n
USE MODI_ASSIM_SURF_ATM_n
USE MODI_WRITE_SURF_ATM_n
USE MODI_WRITE_DIAG_SURF_ATM_n
USE MODI_ASSIM_SET_SST
USE MODI_ADD_FORECAST_TO_DATE_SURF
USE MODI_FLAG_UPDATE
USE MODI_FLAG_DIAG_UPDATE
!
USE MODE_EKF, ONLY : GET_FILE_NAME
!
#ifdef OFF
USE MODI_FANDAR
#endif
!
USE YOMHOOK,             ONLY : LHOOK,DR_HOOK
USE PARKIND1,            ONLY : JPRB
!
IMPLICIT NONE
!
#ifndef NOMPI
INCLUDE 'mpif.h'
#endif
!
!*    0.     Declaration of local variables
!            ------------------------------
!
 CHARACTER(LEN=200) :: YMFILE     ! Name of the observation, perturbed or reference file!
 CHARACTER(LEN=3)  :: YINIT
 CHARACTER(LEN=2), PARAMETER  :: YTEST        = 'OK'          ! must be equal to 'OK'
 CHARACTER(LEN=28)            :: YATMFILE  ='   '  ! name of the Atmospheric file
 CHARACTER(LEN=6)             :: YATMFILETYPE ='      '                     ! type of the Atmospheric file
 CHARACTER(LEN=28)            :: YLUOUT    ='LISTING_SODA                '  ! name of listing
 CHARACTER(LEN=6)             :: YPROGRAM2 = 'FA    '
 CHARACTER(LEN=28)            :: YFILEIN
 CHARACTER(LEN=1)             :: YVAR
!
REAL,ALLOCATABLE, DIMENSION(:)   :: ZLSM                ! Land-Sea mask
REAL,ALLOCATABLE, DIMENSION(:)   :: ZCON_RAIN           ! Amount of convective liquid precipitation
REAL,ALLOCATABLE, DIMENSION(:)   :: ZSTRAT_RAIN         ! Amount of stratiform liquid precipitation
REAL,ALLOCATABLE, DIMENSION(:)   :: ZCON_SNOW           ! Amount of convective solid precipitation
REAL,ALLOCATABLE, DIMENSION(:)   :: ZSTRAT_SNOW         ! Amount of stratiform solid precipitation
REAL,ALLOCATABLE, DIMENSION(:)   :: ZCON_GRAUPEL        ! Amount of convective graupel pecipitation (AROME)
REAL,ALLOCATABLE, DIMENSION(:)   :: ZCLOUDS             ! Cloudcover
REAL,ALLOCATABLE, DIMENSION(:)   :: ZEVAPTR             ! Evaporation
REAL,ALLOCATABLE, DIMENSION(:)   :: ZEVAP               ! Evaporation
REAL,ALLOCATABLE, DIMENSION(:)   :: ZTSC                ! Climatological surface temperature
REAL,ALLOCATABLE, DIMENSION(:)   :: ZTS                 ! Surface temperature
REAL,ALLOCATABLE, DIMENSION(:)   :: ZT2M                ! Screen level temperature
REAL,ALLOCATABLE, DIMENSION(:)   :: ZHU2M               ! Screen level relative humidity
REAL,ALLOCATABLE, DIMENSION(:)   :: ZSNC
REAL,ALLOCATABLE, DIMENSION(:)   :: ZSWE                ! Snow water equvivalent (amount of snow on the ground)
REAL,ALLOCATABLE, DIMENSION(:)   :: ZSWEC               ! Climatological snow water equvivalent (amount of snow on the ground)
REAL,ALLOCATABLE, DIMENSION(:)   :: ZUCLS
REAL,ALLOCATABLE, DIMENSION(:)   :: ZVCLS
REAL,ALLOCATABLE, DIMENSION(:)   :: ZSST                ! SST from external file
REAL,ALLOCATABLE, DIMENSION(:)   :: ZSIC                ! SIC from external file
REAL,ALLOCATABLE, DIMENSION(:)   :: ZLAT
REAL,ALLOCATABLE, DIMENSION(:)   :: ZLON
!
REAL    :: ZTIME
REAL    :: ZTIME_OUT           ! output time since start of the run (s)
REAL(KIND=JPRB)                  :: ZHOOK_HANDLE
!
LOGICAL, ALLOCATABLE, DIMENSION(:) :: GD_MASKEXT
LOGICAL :: GLKEEPEXTZONE
LOGICAL :: GFOUND
!
TYPE (DATE_TIME)                 :: TTIME               ! Current date and time  
!
INTEGER, DIMENSION(11)  :: IDATEF
INTEGER :: INI          ! grid dimension
INTEGER :: ISV                 ! Number of scalar species
INTEGER :: ISW                 ! Number of radiative bands 
INTEGER :: IYEAR, IMONTH, IDAY, IHOUR
INTEGER :: IYEAR_OUT, IMONTH_OUT, IDAY_OUT
INTEGER :: L,I,J,INBPERT
INTEGER :: INW, JNW
INTEGER :: ISTEP
INTEGER :: IOBS, IIOBS
INTEGER :: IGPCOMP
INTEGER :: ILUOUT
INTEGER :: ILUNAM
INTEGER :: IRET, INB
INTEGER :: IRESP, ISTAT               ! Response value
INTEGER :: INFOMPI, ILEVEL
!
! ******************************************************************************************
!
#ifndef NOMPI
 CALL MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,ILEVEL,INFOMPI)
#endif
!
IF (LHOOK) CALL DR_HOOK('SODA',0,ZHOOK_HANDLE)
!
!
#ifndef NOMPI
NCOMM = MPI_COMM_WORLD
 CALL MPI_COMM_SIZE(NCOMM,NPROC,INFOMPI)
 CALL MPI_COMM_RANK(NCOMM,NRANK,INFOMPI)
#endif
!
!
IF (LHOOK) CALL DR_HOOK('SODA',0,ZHOOK_HANDLE)
!
 CALL PREP_LOG_MPI
!
!--------------------------------------
!
WRITE(*,*)
WRITE(*,*) '   ------------------------------------'
WRITE(*,*) '   |               SODA               |'
WRITE(*,*) '   | SURFEX OFFLINE DATA ASSIMILATION |'
WRITE(*,*) '   ------------------------------------'
WRITE(*,*)

CSOFTWARE = 'SODA'

! Read offline specific things
 CALL OPEN_NAMELIST('ASCII ',ILUNAM,CNAMELIST)
 CALL POSNAM(ILUNAM,'NAM_IO_OFFLINE',GFOUND)
IF (GFOUND) READ (UNIT=ILUNAM,NML=NAM_IO_OFFLINE)
 CALL CLOSE_NAMELIST('ASCII ',ILUNAM)

! Open ascii outputfile for writing
#ifdef LFI
CLUOUT_LFI =  ADJUSTL(ADJUSTR(YLUOUT)//'.txt')
#endif
#ifdef NC
CLUOUT_NC = ADJUSTL(ADJUSTR(YLUOUT)//'.txt')
#endif
 CALL GET_LUOUT(CSURF_FILETYPE,ILUOUT)
OPEN(UNIT=ILUOUT,FILE=ADJUSTL(ADJUSTR(YLUOUT)//'.txt'),FORM='FORMATTED',ACTION='WRITE')


! Reading all namelist (also assimilation)
CALL READ_ALL_NAMELISTS(CSURF_FILETYPE,'ALL',.FALSE.)

! Go to SURFEX
CALL ALLOC_SURFEX(1)
CALL GOTO_SURFEX(1,.TRUE.)
!
! Setting input files read from namelist
IF ( CSURF_FILETYPE == "LFI   " ) THEN
#ifdef LFI
  CFILEIN_LFI      = CPREPFILE
  CFILE_LFI        = CPREPFILE
  CFILEIN_LFI_SAVE = CPREPFILE
  CFILEPGD_LFI     = CPGDFILE
  CFILEOUT_LFI     = CSURFFILE
#endif
ELSEIF ( CSURF_FILETYPE == "FA    " ) THEN
#ifdef FA
  CFILEIN_FA      = ADJUSTL(ADJUSTR(CPREPFILE)//'.fa')
  CFILEIN_FA_SAVE = ADJUSTL(ADJUSTR(CPREPFILE)//'.fa')
  CFILEPGD_FA     = ADJUSTL(ADJUSTR(CPGDFILE)//'.fa')
  CFILEOUT_FA  = ADJUSTL(ADJUSTR(CSURFFILE)//'.fa')
#endif
ELSEIF ( CSURF_FILETYPE == "ASCII " ) THEN
#ifdef ASC
  CFILEIN      = ADJUSTL(ADJUSTR(CPREPFILE)//'.txt')
  CFILEIN_SAVE = ADJUSTL(ADJUSTR(CPREPFILE)//'.txt')
  CFILEPGD     = ADJUSTL(ADJUSTR(CPGDFILE)//'.txt')
  CFILEOUT = ADJUSTL(ADJUSTR(CSURFFILE)//'.txt')
#endif
ELSEIF ( CSURF_FILETYPE == "NC    " ) THEN
#ifdef ASC
  CFILEIN_NC      = ADJUSTL(ADJUSTR(CPREPFILE)//'.nc')
  CFILEIN_NC_SAVE = ADJUSTL(ADJUSTR(CPREPFILE)//'.nc')
  CFILEPGD_NC     = ADJUSTL(ADJUSTR(CPGDFILE)//'.nc')
  CFILEOUT_NC = ADJUSTL(ADJUSTR(CSURFFILE)//'.nc')
#endif
ELSE
  CALL ABOR1_SFX(TRIM(CSURF_FILETYPE)//" is not implemented!")
ENDIF
!
! Initialize time information
IYEAR    = NUNDEF
IMONTH   = NUNDEF
IDAY     = NUNDEF
ZTIME    = XUNDEF
 CALL INIT_IO_SURF_n(CSURF_FILETYPE,'FULL  ','SURF  ','READ ')
 CALL READ_SURF(CSURF_FILETYPE,'DIM_FULL  ',INI,  IRESP)
 CALL READ_SURF(CSURF_FILETYPE,'DTCUR     ',TTIME,  IRESP)
 CALL END_IO_SURF_n(CSURF_FILETYPE)
!
NINDX2SFX = INI
NSIZE = INI
IF (ALLOCATED(NMASK_FULL)) DEALLOCATE(NMASK_FULL)
ALLOCATE(NWORK(INI))
ALLOCATE(XWORK(INI))
ALLOCATE(NWORK2(INI,10))
ALLOCATE(XWORK2(INI,10))
ALLOCATE(XWORK3(INI,10,10))
IF (NRANK==NPIO) THEN
  ALLOCATE(NWORK_FULL(INI))
  ALLOCATE(XWORK_FULL(INI))
  ALLOCATE(NWORK2_FULL(INI,10))  
  ALLOCATE(XWORK2_FULL(INI,10))
ELSE
  ALLOCATE(NWORK_FULL(0))
  ALLOCATE(XWORK_FULL(0))
  ALLOCATE(NWORK2_FULL(0,0))  
  ALLOCATE(XWORK2_FULL(0,0))
ENDIF
!
ISW = 0
ISV = 0
ALLOCATE(CSV(ISV))
ALLOCATE(XCO2(INI))
ALLOCATE(XRHOA(INI))
ALLOCATE(XZENITH(INI))
ALLOCATE(XAZIM(INI))
ALLOCATE(XSW_BANDS(ISW))
ALLOCATE(XDIR_ALB(INI,ISW))
ALLOCATE(XSCA_ALB(INI,ISW))
ALLOCATE(XEMIS(INI))
ALLOCATE(XTSRAD(INI))
!
! Indicate that zenith and azimuth angles are not initialized
XZENITH = XUNDEF
XAZIM   = XUNDEF
!
IF ( CASSIM_ISBA == 'EKF  ' ) THEN
 ! Has to do initialization for all the perturbations + 
 ! control + the real run at last
 INBPERT = NVAR + 2
ELSE
 INBPERT = 1
ENDIF

WRITE(*,*) "INITIALIZING SURFEX..."
!
YINIT = 'ALL'
!
IYEAR  = TTIME%TDATE%YEAR
IMONTH = TTIME%TDATE%MONTH
IDAY   = TTIME%TDATE%DAY
ZTIME  = TTIME%TIME
!
NOBS = 0
!
IHOUR = 0
ZTIME = FLOAT(NECHGU) * 3600.
! BEGINNING OF TIME LOOP
DO ISTEP = 1,NBOUTPUT
  !
  ! Update date
  CALL ADD_FORECAST_TO_DATE_SURF(IYEAR, IMONTH, IDAY, ZTIME)
  ZTIME = ZTIME + FLOAT(NECHGU) * 3600.
  IHOUR = IHOUR + NECHGU
  !  
  DO NIPERT = INBPERT,1,-1
    !
    ! If we have more than one initialization to do
    ! For last initialization, we must re-do the first.
    ! Could be avoided by introducing knowlegde of LASSIM on this level
    IF ( CASSIM_ISBA == 'EKF  ' ) THEN
      !    
      IF ( NIPERT<INBPERT ) THEN
        YMFILE = "PREP_"
        CALL GET_FILE_NAME(IYEAR,IMONTH,IDAY,IHOUR,YMFILE)              
        WRITE(YVAR,'(I1.1)') NIPERT-1
        YFILEIN = TRIM(YMFILE)//"_EKF_PERT"//YVAR
      ELSE
        YFILEIN = "PREP_INIT"
      ENDIF
      !
      IF ( CSURF_FILETYPE == "LFI   " ) THEN
        YFILEIN = TRIM(YFILEIN)
#ifdef LFI
        CFILEIN_LFI      = YFILEIN 
        CFILE_LFI        = YFILEIN
        CFILEIN_LFI_SAVE = YFILEIN
#endif    
      ELSEIF ( CSURF_FILETYPE == "FA    " ) THEN
        YFILEIN = TRIM(YFILEIN)//'.fa'
#ifdef FA
        CFILEIN_FA      = YFILEIN
        CFILEIN_FA_SAVE = YFILEIN
#endif    
      ELSEIF ( CSURF_FILETYPE == "ASCII " ) THEN
        YFILEIN = TRIM(YFILEIN)//'.txt'
#ifdef ASC
        CFILEIN      = YFILEIN
        CFILEIN_SAVE = YFILEIN
#endif    
      ELSEIF ( CSURF_FILETYPE == "NC    " ) THEN
        YFILEIN = TRIM(YFILEIN)//'.nc'
#ifdef NC
        CFILEIN_NC      = YFILEIN
        CFILEIN_NC_SAVE = YFILEIN
#endif    
      ELSE
        CALL ABOR1_SFX(TRIM(CSURF_FILETYPE)//" is not implemented!")
      ENDIF
      !
    ENDIF 
    !
    ! Initialize the SURFEX interface
    CALL IO_BUFF_CLEAN_n
    CALL INIT_SURF_ATM_n(CSURF_FILETYPE,YINIT, LLAND_USE, INI, ISV, ISW,  &
                         CSV, XCO2, XRHOA, XZENITH, XAZIM, XSW_BANDS,     &
                         XDIR_ALB, XSCA_ALB, XEMIS, XTSRAD,               &
                         IYEAR, IMONTH, IDAY, ZTIME,                      &
                          YATMFILE, YATMFILETYPE, YTEST              )
    !
    IF ( CASSIM_ISBA=='EKF  ' ) THEN
      !
      IF ( ISTEP==1 .AND. NIPERT==INBPERT ) THEN
        ALLOCATE(XLAI_PASS(NSIZE_NATURE,NPATCH)) 
        ALLOCATE(XBIO_PASS(NSIZE_NATURE,NPATCH))     
        ALLOCATE(XI       (NSIZE_NATURE,NPATCH,NVAR    ))
        ALLOCATE(XF       (NSIZE_NATURE,NPATCH,NVAR+1,NVAR    ))
        ALLOCATE(XF_PATCH (NSIZE_NATURE,NPATCH,NVAR+1,NOBSTYPE*NBOUTPUT))
      ENDIF
      !
      IF ( NIPERT<INBPERT ) THEN
        !
        ! Set the global state values for this control value
        DO IOBS = 1,NOBSTYPE
          IIOBS = (ISTEP-1)*NOBSTYPE + IOBS
          SELECT CASE (TRIM(COBS(IOBS)))
            CASE("T2M")
              XF_PATCH(:,:,NIPERT,IIOBS) = XAT2M_ISBA(:,:)
            CASE("HU2M")
              XF_PATCH(:,:,NIPERT,IIOBS) = XAHU2M_ISBA(:,:)
            CASE("WG1")
              XF_PATCH(:,:,NIPERT,IIOBS) = XWG(:,1,:)
            CASE("LAI")
              XF_PATCH(:,:,NIPERT,IIOBS) = XLAI(:,:)
            CASE DEFAULT
              CALL ABOR1_SFX("Mapping of "//COBS(IOBS)//" is not defined in SODA_CONTROL!")
          END SELECT
        ENDDO
        !
        ! Prognostic fields for assimilation (Control vector)
        DO L = 1,NVAR
          SELECT CASE (TRIM(CVAR(L)))
            CASE("TG1")
              XF(:,:,NIPERT,L) = XTG(:,1,:)
            CASE("TG2")
              XF(:,:,NIPERT,L) = XTG(:,2,:)
            CASE("WG1")
              XF(:,:,NIPERT,L) = XWG(:,1,:)
            CASE("WG2")
              XF(:,:,NIPERT,L) = XWG(:,2,:)
            CASE("LAI")
              XF(:,:,NIPERT,L) = XLAI(:,:)
            CASE DEFAULT
              CALL ABOR1_SFX("Mapping of "//TRIM(CVAR(L))//" is not defined in SODA_CONTROL!")
          END SELECT
        ENDDO
        !
        IF ( NIPERT==1 ) THEN
          !
          IF ( NPATCH==1 .AND. TRIM(CBIO)/="LAI" ) THEN
            CALL ABOR1_SFX("Mapping of "//CBIO//" is not defined in EKF with NPATCH=1!")
          ENDIF
          SELECT CASE (TRIM(CBIO))
            CASE("BIOMA1","BIOMASS1")
              XBIO_PASS(:,:) = XBIOMASS(:,1,:)
            CASE("BIOMA2","BIOMASS2")
              XBIO_PASS(:,:) = XBIOMASS(:,2,:)
            CASE("RESPI1","RESP_BIOM1")
              XBIO_PASS(:,:) = XRESP_BIOMASS(:,1,:)
            CASE("RESPI2","RESP_BIOM2")
              XBIO_PASS(:,:) = XRESP_BIOMASS(:,2,:)
            CASE("LAI")
              XBIO_PASS(:,:) = XLAI(:,:)
            CASE DEFAULT
              CALL ABOR1_SFX("Mapping of "//CBIO//" is not defined in EKF!")
          END SELECT
          !
        ENDIF
        !
      ELSE
        !
        DO L = 1,NVAR
          SELECT CASE (TRIM(CVAR(L)))
            CASE("TG1")
              XI(:,:,L) = XTG(:,1,:)
            CASE("TG2")
              XI(:,:,L) = XTG(:,2,:)
            CASE("WG1")
              XI(:,:,L) = XWG(:,1,:)
            CASE("WG2")
              XI(:,:,L) = XWG(:,2,:)
            CASE("LAI")
              XI(:,:,L) = XLAI(:,:)
            CASE DEFAULT
              CALL ABOR1_SFX("Mapping of "//TRIM(CVAR(L))//" is not defined in SODA_CONTROL!")
          END SELECT
        ENDDO        
        XLAI_PASS(:,:) = XLAI(:,:)
        !
      ENDIF
      !
    ENDIF
    !
  ENDDO
  !
  IF ( CASSIM_ISBA=="EKF  " ) THEN
    !
    IF (ISTEP==1) ALLOCATE(XYO(NSIZE_NATURE,NOBSTYPE*NBOUTPUT))
    !
    IF ( LOBSFILE ) THEN
      !
      YMFILE = 'CANARI_NATURE_'
      CALL GET_FILE_NAME(IYEAR,IMONTH,IDAY,IHOUR,YMFILE)
      OPEN(UNIT=55,FILE=TRIM(YMFILE)//".DAT",FORM='FORMATTED',STATUS='OLD',IOSTAT=ISTAT) 
      IF ( ISTAT==0 ) THEN
        !   If it exists, read observations
        DO I = 1,NSIZE_NATURE
          READ (55,*)  (XYO(I,NOBS+J),J=1,NOBSTYPE)
        ENDDO
        NOBS = NOBS + NOBSTYPE      
        IF ( NPRINTLEV > 0 ) WRITE(*,*) 'read in obs: ', XYO(1,:), NOBS
        CLOSE(55)
      ENDIF
      !
    ELSE
      !
      NOBS = NOBSTYPE
      !
    ENDIF
    !
  ENDIF
  !
ENDDO
!
IF ( NOBS==0 .AND. CASSIM_ISBA=="EKF  " ) THEN
  IF ( NPRINTLEV > 0 ) WRITE(*,*) 'No observations read for LAI in OBS file - stop'
  CALL ABOR1_SFX("SODA_CONTROL: No observations read for LAI in OBS file - stop")
ENDIF
!
 CALL GET_SIZE_FULL_n(CSURF_FILETYPE,INI,NSIZE_FULL)
!
WRITE(*,*) "READING input files..."
! Normal reading of needed FA fields
ALLOCATE(ZLSM        (INI))
ALLOCATE(ZCON_RAIN   (INI))
ALLOCATE(ZSTRAT_RAIN (INI))
ALLOCATE(ZCON_SNOW   (INI))
ALLOCATE(ZSTRAT_SNOW (INI))
ALLOCATE(ZCON_GRAUPEL(INI))
ALLOCATE(ZCLOUDS     (INI))
ALLOCATE(ZEVAPTR     (INI))
ALLOCATE(ZEVAP       (INI))
ALLOCATE(ZTSC        (INI))
ALLOCATE(ZSWE        (INI))
ALLOCATE(ZSWEC       (INI))
ALLOCATE(ZTS         (INI))
ALLOCATE(ZT2M        (INI))
ALLOCATE(ZHU2M       (INI))
ALLOCATE(ZUCLS       (INI))
ALLOCATE(ZVCLS       (INI))
ALLOCATE(ZSST        (INI))
ALLOCATE(ZSIC        (INI))

!  Read atmospheric forecast fields from FA files 
#ifdef FA
CFILEIN_FA = 'FG_OI_MAIN'
CDNOMC     = 'oimain'                  ! new frame name
#endif
!  Open FA file (LAM version with extension zone)
 CALL INIT_IO_SURF_n(YPROGRAM2,'EXTZON','SURF  ','READ ')
!
!  Read model forecast quantities
IF (LAROME) THEN  
  CALL READ_SURF(YPROGRAM2,'SURFACCPLUIE',    ZCON_RAIN    ,IRESP)
  CALL READ_SURF(YPROGRAM2,'SURFACCNEIGE',    ZCON_SNOW    ,IRESP)
  CALL READ_SURF(YPROGRAM2,'SURFACCGRAUPEL',  ZCON_GRAUPEL ,IRESP)
  ! So far graupel has not been used
  !ZCON_SNOW=ZCON_SNOW+ZCON_GRAUPEL
  ZSTRAT_RAIN(:) = 0.0
  ZSTRAT_SNOW(:) = 0.0  
ELSE    
  CALL READ_SURF(YPROGRAM2,'SURFPREC.EAU.CON',ZCON_RAIN    ,IRESP)
  CALL READ_SURF(YPROGRAM2,'SURFPREC.EAU.GEC',ZSTRAT_RAIN  ,IRESP)
  CALL READ_SURF(YPROGRAM2,'SURFPREC.NEI.CON',ZCON_SNOW    ,IRESP)
  CALL READ_SURF(YPROGRAM2,'SURFPREC.NEI.GEC',ZSTRAT_SNOW  ,IRESP)
  ZCON_GRAUPEL(:) = 0.0
ENDIF
!
 CALL READ_SURF(YPROGRAM2,'ATMONEBUL.BASSE ',ZCLOUDS,IRESP)
 CALL READ_SURF(YPROGRAM2,'SURFIND.TERREMER',ZLSM   ,IRESP)
 CALL READ_SURF(YPROGRAM2,'SURFFLU.LAT.MEVA',ZEVAP  ,IRESP) ! accumulated fluxes (not available in LFI)
!
IF (.NOT.LALADSURF) THEN    
  CALL READ_SURF(YPROGRAM2,'SURFXEVAPOTRANSP',ZEVAPTR,IRESP) ! not in ALADIN SURFEX
ELSE
  ZEVAPTR(:) = 0.0
ENDIF

!  Close FA file
 CALL END_IO_SURF_n(YPROGRAM2)
 CALL IO_BUFF_CLEAN_n
WRITE(*,*)'READ FG_OI_MAIN OK'

!  Define FA file name for CANARI analysis
#ifdef FA
CFILEIN_FA = 'CANARI'        ! input CANARI analysis
CDNOMC     = 'canari'                  ! new frame name
#endif
!  Open FA file 
 CALL INIT_IO_SURF_n(YPROGRAM2,'EXTZON','SURF  ','READ ')

!  Read CANARI analysis
 CALL READ_SURF(YPROGRAM2,'CLSTEMPERATURE  ',ZT2M ,IRESP)
 CALL READ_SURF(YPROGRAM2,'CLSHUMI.RELATIVE',ZHU2M,IRESP)
 CALL READ_SURF(YPROGRAM2,'SURFTEMPERATURE ',ZTS  ,IRESP)
 CALL READ_SURF(YPROGRAM2,'SURFRESERV.NEIGE',ZSWE ,IRESP)
 CALL READ_SURF(YPROGRAM2,'CLSVENT.ZONAL   ',ZUCLS,IRESP)
 CALL READ_SURF(YPROGRAM2,'CLSVENT.MERIDIEN',ZVCLS,IRESP)  

!  Close CANARI file
 CALL END_IO_SURF_n(YPROGRAM2)
 CALL IO_BUFF_CLEAN_n
 WRITE(*,*) 'READ CANARI OK'

!  Define FA file name for surface climatology
#ifdef FA
CFILEIN_FA = 'clim_isba'               ! input climatology
CDNOMC     = 'climat'                  ! new frame name
#endif
!  Open FA file 
 CALL INIT_IO_SURF_n(YPROGRAM2,'EXTZON','SURF  ','READ ')

!  Read climatology file
 CALL READ_SURF(YPROGRAM2,'SURFRESERV.NEIGE',ZSWEC,IRESP)
 CALL READ_SURF(YPROGRAM2,'SURFTEMPERATURE' ,ZTSC ,IRESP)

!  Close climatology file
 CALL END_IO_SURF_n(YPROGRAM2)
 CALL IO_BUFF_CLEAN_n
WRITE(*,*) 'READ CLIMATOLOGY OK'

 CALL ASSIM_SET_SST(NSIZE_FULL,ZLSM,ZSST,ZSIC,YTEST)

IF ( .NOT. LASSIM ) CALL ABOR1_SFX("YOU CAN'T RUN SODA WITHOUT SETTING LASSIM=.TRUE. IN THE ASSIM NAMELIST")

ALLOCATE(GD_MASKEXT(INI))
GD_MASKEXT(:) = .FALSE.
!
ALLOCATE(ZLON(INI))
ALLOCATE(ZLAT(INI))
ZLON(:) = XLON(:)
ZLAT(:) = XLAT(:)        
!
GLKEEPEXTZONE = .TRUE.
!
WRITE(*,*) 'PERFORMIMG OFFLINE SURFEX DATA ASSIMILATION...'
 CALL ASSIM_SURF_ATM_n(CSURF_FILETYPE,INI,   &
                      ZCON_RAIN,  ZSTRAT_RAIN, ZCON_SNOW, ZSTRAT_SNOW, &
                      ZCLOUDS,    ZLSM,        ZEVAPTR,   ZEVAP,       &
                      ZSWEC,      ZTSC,       &
                      ZTS,        ZT2M,        ZHU2M,     ZSWE,        &
                      ZSST,       ZSIC,       &
                      YTEST, GD_MASKEXT, ZLON, ZLAT, GLKEEPEXTZONE )

DEALLOCATE(ZCON_RAIN)
DEALLOCATE(ZSTRAT_RAIN)
DEALLOCATE(ZCON_SNOW)
DEALLOCATE(ZSTRAT_SNOW)
DEALLOCATE(ZCON_GRAUPEL)
DEALLOCATE(ZCLOUDS)
DEALLOCATE(ZLSM)
DEALLOCATE(ZEVAPTR)
DEALLOCATE(ZEVAP)
DEALLOCATE(ZTSC)
DEALLOCATE(ZSWEC)
DEALLOCATE(ZTS)
DEALLOCATE(ZT2M)
DEALLOCATE(ZHU2M)
DEALLOCATE(ZSWE)
DEALLOCATE(ZUCLS)
DEALLOCATE(ZVCLS)
DEALLOCATE(ZSST)
DEALLOCATE(ZSIC)
!
!       
INW = 1
IF (CSURF_FILETYPE=="NC    ") INW = 2
!
ZTIME_OUT  = ZTIME
IDAY_OUT   = IDAY
IMONTH_OUT = IMONTH
IYEAR_OUT  = IYEAR
!
IF(LOUT_TIMENAME)THEN
  ! if true, change the name of output file at the end of a day
  ! (ex: 19860502_00h00 -> 19860501_24h00)                     
  IF(ZTIME==0.0)THEN
    ZTIME_OUT = 86400.
    IDAY_OUT   = IDAY-1
    IF(IDAY_OUT==0)THEN
      IMONTH_OUT = IMONTH - 1
      IF(IMONTH_OUT==0)THEN
        IMONTH_OUT=12
        IYEAR_OUT = IYEAR - 1
      ENDIF
      SELECT CASE (IMONTH_OUT)
        CASE(4,6,9,11)
          IDAY_OUT=30
        CASE(1,3,5,7:8,10,12)
          IDAY_OUT=31
        CASE(2)
          IF( ((MOD(IYEAR_OUT,4)==0).AND.(MOD(IYEAR_OUT,100)/=0)) .OR. (MOD(IYEAR_OUT,400)==0))THEN 
            IDAY_OUT=29
          ELSE
           IDAY_OUT=28
         ENDIF
      END SELECT
    ENDIF
  ENDIF
  !
ENDIF
!
IF (CSURF_FILETYPE=='FA    ') THEN
  CDNOMC = 'header'
  LFANOCOMPACT = LDIAG_FA_NOCOMPACT
  IDATEF(1)= IYEAR_OUT
  IDATEF(2)= IMONTH_OUT
  IDATEF(3)= IDAY_OUT
  IDATEF(4)= FLOOR(ZTIME_OUT/3600.)
  IDATEF(5)= FLOOR(ZTIME_OUT/60.) - IDATEF(4) * 60 
  IDATEF(6)= NINT(ZTIME_OUT) - IDATEF(4) * 3600 - IDATEF(5) * 60
  IDATEF(7:11) = 0
  CALL FAITOU(IRET,NUNIT_FA,.TRUE.,CFILEOUT_FA,'UNKNOWN',.TRUE.,.FALSE.,IVERBFA,0,INB,CDNOMC)
  CALL FANDAR(IRET,NUNIT_FA,IDATEF)
END IF
!
LDEF = .TRUE.
DO JNW = 1,INW
  CALL FLAG_UPDATE(.FALSE.,.TRUE.,.FALSE.,.FALSE.)
  CALL FLAG_DIAG_UPDATE(.FALSE.,.TRUE.,0,.FALSE.,.FALSE.,.FALSE.,&
                        .FALSE.,0,0,.FALSE.,.FALSE.,&
                        .FALSE.,.FALSE.,.FALSE.,.FALSE.,&
                        .FALSE.,.FALSE.,.FALSE.,&
                        .FALSE.,.FALSE.)  
  ! Store results from assimilation
  CALL WRITE_SURF_ATM_n(CSURF_FILETYPE,'ALL',LLAND_USE)
  CALL WRITE_DIAG_SURF_ATM_n(CSURF_FILETYPE,'ALL')
  LDEF = .FALSE.
  CALL IO_BUFF_CLEAN_n
  !
ENDDO  
!
IF (CSURF_FILETYPE=='FA    ') THEN
  CALL FAIRME(IRET,NUNIT_FA,'UNKNOWN')
END IF
!
!*    3.     Close parallelized I/O
!            ----------------------
!
WRITE(ILUOUT,*) ' '
WRITE(ILUOUT,*) '    -----------------------'
WRITE(ILUOUT,*) '    | SODA ENDS CORRECTLY |'
WRITE(ILUOUT,*) '    -----------------------'
!
WRITE(*,*) ' '
WRITE(*,*) '    -----------------------'
WRITE(*,*) '    | SODA ENDS CORRECTLY |'
WRITE(*,*) '    -----------------------'
!
CLOSE(ILUOUT)
!
DEALLOCATE(NWORK)
DEALLOCATE(XWORK)
DEALLOCATE(XWORK2)
DEALLOCATE(XWORK3)
DEALLOCATE(NWORK_FULL)
DEALLOCATE(XWORK_FULL)
DEALLOCATE(XWORK2_FULL)
!
CALL DEALLOC_SURFEX
!
IF (LHOOK) CALL DR_HOOK('SODA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END PROGRAM SODA
