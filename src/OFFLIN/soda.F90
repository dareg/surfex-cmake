!
! *****************************************************************************************
PROGRAM SODA
! ******************************************************************************************
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
!----------------------------------------------------------------------------
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO
USE MODD_SURFEX_OMP, ONLY : NINDX2, NWORK, XWORK, XWORK2, XWORK3, &
                            NWORK_FULL, XWORK_FULL, XWORK2_FULL
!
USE MODD_TYPE_DATE_SURF, ONLY : DATE_TIME
USE MODD_SURF_PAR,       ONLY : XUNDEF,NUNDEF
USE MODD_ASSIM,          ONLY : LAROME,LALADSURF,LPRT,LSIM,LBEV,CASSIM_ISBA
!
USE MODD_FORC_ATM,       ONLY : CSV         ,&! name of all scalar variables
                                XDIR_ALB    ,&! direct albedo for each band
                                XSCA_ALB    ,&! diffuse albedo for each band
                                XEMIS       ,&! emissivity
                                XTSRAD      ,&! radiative temperature
                                XTSUN       ,&! solar time                    (s from midnight)
                                XZS         ,&! orography                             (m)
                                XZREF       ,&! height of T,q forcing                 (m)
                                XUREF       ,&! height of wind forcing                (m)
                                XTA         ,&! air temperature forcing               (K)
                                XQA         ,&! air humidity forcing                  (kg/m3)
                                XSV         ,&! scalar variables
                                XU          ,&! zonal wind                            (m/s)
                                XV          ,&! meridian wind                         (m/s)
                                XSW_BANDS   ,&! mean wavelength of each shortwave band (m)
                                XZENITH     ,&! zenithal angle       (radian from the vertical)
                                XAZIM       ,&! azimuthal angle      (radian from North, clockwise)
                                XCO2        ,&! CO2 concentration in the air          (kg/m3)
                                XRHOA         ! air density  at the surface           (kg/m3)
!
USE MODD_SURF_ATM_n,     ONLY : NDIM_NATURE
!
#ifdef ASC
USE MODD_IO_SURF_ASC,    ONLY : CFILEIN,CFILEIN_SAVE,&
                                CFILEPGD
#endif
#ifdef FA
USE MODD_IO_SURF_FA,     ONLY : CFILEIN_FA,CFILEIN_FA_SAVE,&
                                CFILEPGD_FA,CDNOMC  
#endif
#ifdef LFI
USE MODD_IO_SURF_LFI,    ONLY : CFILEIN_LFI,CFILEIN_LFI_SAVE,&
                                CFILEPGD_LFI,CFILE_LFI,CLUOUT_LFI 
#endif
!
USE MODN_IO_OFFLINE,     ONLY : NAM_IO_OFFLINE,CNAMELIST,CPGDFILE,CPREPFILE,&
                                CSURF_FILETYPE,LLAND_USE,LRESTART
!
USE MODE_POS_SURF,       ONLY : POSNAM
!
USE YOMHOOK,             ONLY : LHOOK,DR_HOOK
USE PARKIND1,            ONLY : JPRB
!
USE MODI_ALLOC_SURFEX
USE MODI_GET_LUOUT
USE MODI_OPEN_NAMELIST
USE MODI_CLOSE_NAMELIST
USE MODI_ABOR1_SFX
USE MODI_INIT_IO_SURF_n
USE MODI_READ_SURF
USE MODI_END_IO_SURF_n
USE MODI_READ_ALL_NAMELISTS
USE MODI_GOTO_SURFEX
USE MODI_GOTO_TRIP
USE MODI_INIT_SURF_ATM_n
USE MODI_INIT_SURF_TRIP_n
USE MODI_ASSIM_NATURE_ISBA_EKF
USE MODI_IO_BUFF_CLEAN_n
USE MODI_ASSIM_SURF_ATM_n
USE MODI_DEALLOC_SURFEX
!
IMPLICIT NONE
!
!*    0.     Declaration of local variables
!            ------------------------------
!
 CHARACTER(LEN=3), PARAMETER      :: YINIT        = 'ALL'
 CHARACTER(LEN=2), PARAMETER      :: YTEST        = 'OK'          ! must be equal to 'OK'
INTEGER                          :: ILUOUT
INTEGER                          :: ILUNAM
INTEGER                          :: IYEAR, IMONTH, IDAY,IHOUR,NSSSSS,ndaysec
REAL                             :: ZTIME
LOGICAL                          :: GFOUND
 CHARACTER(LEN=28)                :: YATMFILE  ='                            '  ! name of the Atmospheric file
 CHARACTER(LEN=6)                 :: YATMFILETYPE ='      '                     ! type of the Atmospheric file
 CHARACTER(LEN=28)                :: YLUOUT    ='LISTING_SODA                '  ! name of listing
INTEGER                          :: IRET, INB
REAL(KIND=JPRB)                  :: ZHOOK_HANDLE
REAL                             :: ZTIMEC              ! current duration since start of the run (s)
REAL,ALLOCATABLE, DIMENSION(:)   :: PCON_RAIN           ! Amount of convective liquid precipitation
REAL,ALLOCATABLE, DIMENSION(:)   :: PSTRAT_RAIN         ! Amount of stratiform liquid precipitation
REAL,ALLOCATABLE, DIMENSION(:)   :: PCON_SNOW           ! Amount of convective solid precipitation
REAL,ALLOCATABLE, DIMENSION(:)   :: PSTRAT_SNOW         ! Amount of stratiform solid precipitation
REAL,ALLOCATABLE, DIMENSION(:)   :: PCON_GRAUPEL        ! Amount of convective graupel pecipitation (AROME)
REAL,ALLOCATABLE, DIMENSION(:)   :: PCLOUDS             ! Cloudcover
REAL,ALLOCATABLE, DIMENSION(:)   :: PLSM                ! Land-Sea mask
REAL,ALLOCATABLE, DIMENSION(:)   :: PEVAPTR             ! Evaporation
REAL,ALLOCATABLE, DIMENSION(:)   :: PEVAP               ! Evaporation
REAL,ALLOCATABLE, DIMENSION(:)   :: PTSC                ! Climatological surface temperature
REAL,ALLOCATABLE, DIMENSION(:)   :: PSWEC               ! Climatological snow water equvivalent (amount of snow on the ground)
REAL,ALLOCATABLE, DIMENSION(:)   :: PTS                 ! Surface temperature
REAL,ALLOCATABLE, DIMENSION(:)   :: PT2M                ! Screen level temperature
REAL,ALLOCATABLE, DIMENSION(:)   :: PHU2M               ! Screen level relative humidity
REAL,ALLOCATABLE, DIMENSION(:)   :: PSWE                ! Snow water equvivalent (amount of snow on the ground)
TYPE (DATE_TIME)                 :: TTIME               ! Current date and time  
 CHARACTER(LEN=6)                 :: YPROGRAM2 = 'FA    '
INTEGER                          :: INI                 ! grid dimension
INTEGER                          :: KSV                 ! Number of scalar species
INTEGER                          :: KSW                 ! Number of radiative bands 
INTEGER                          :: ITRIP_COUNT         ! day counter
INTEGER                          :: IRESP               ! Response value
INTEGER                          :: ITRIP_MONTH         ! mont counter for TRIP
REAL                             :: ZDURATION           ! duration of run                     (s)
INTEGER                          :: I

IF (LHOOK) CALL DR_HOOK('SODA',0,ZHOOK_HANDLE)

WRITE(*,*)
WRITE(*,*) '   ------------------------------------'
WRITE(*,*) '   |               SODA               |'
WRITE(*,*) '   | SURFEX OFFLINE DATA ASSIMILATION |'
WRITE(*,*) '   ------------------------------------'
WRITE(*,*)

! Allocate SURFEX
 CALL ALLOC_SURFEX(1)

! Open ascii outputfile for writing
#ifdef LFI
CLUOUT_LFI =  ADJUSTL(ADJUSTR(YLUOUT)//'.txt')
#endif
 CALL GET_LUOUT('LFI   ',ILUOUT)
OPEN(UNIT=ILUOUT,FILE=ADJUSTL(ADJUSTR(YLUOUT)//'.txt'),FORM='FORMATTED',ACTION='WRITE')

! Read offline specific things
 CALL OPEN_NAMELIST('ASCII ',ILUNAM,CNAMELIST)
 CALL POSNAM(ILUNAM,'NAM_IO_OFFLINE',GFOUND)
IF (GFOUND) READ (UNIT=ILUNAM,NML=NAM_IO_OFFLINE)
 CALL CLOSE_NAMELIST('ASCII ',ILUNAM)

! Setting input files read from namelist
if ( CSURF_FILETYPE == "LFI   " ) then
#ifdef LFI
  CFILEIN_LFI      = CPREPFILE
  CFILE_LFI        = CPREPFILE
  CFILEIN_LFI_SAVE = CPREPFILE
  CFILEPGD_LFI     = CPGDFILE
#endif
elseif ( CSURF_FILETYPE == "FA    " ) then
#ifdef FA
  CFILEIN_FA      = ADJUSTL(ADJUSTR(CPREPFILE)//'.fa')
  CFILEIN_FA_SAVE = ADJUSTL(ADJUSTR(CPREPFILE)//'.fa')
  CFILEPGD_FA     = ADJUSTL(ADJUSTR(CPGDFILE)//'.fa')
#endif
elseif ( CSURF_FILETYPE == "ASCII " ) then
#ifdef ASC
  CFILEIN      = ADJUSTL(ADJUSTR(CPREPFILE)//'.txt')
  CFILEIN_SAVE = ADJUSTL(ADJUSTR(CPREPFILE)//'.txt')
  CFILEPGD     = ADJUSTL(ADJUSTR(CPGDFILE)//'.txt')
#endif
else
  CALL ABOR1_SFX(TRIM(CSURF_FILETYPE)//" is not implemented!")
endif


! Reading all namelist (also assimilation)
 CALL READ_ALL_NAMELISTS(CSURF_FILETYPE,'ALL',.FALSE.)

! Go to SURFEX
 CALL GOTO_SURFEX(1,.TRUE.)
 CALL GOTO_TRIP(1,.TRUE.)

! Initialize time information
IYEAR    = NUNDEF
IMONTH   = NUNDEF
IDAY     = NUNDEF
ZTIME    = XUNDEF
 CALL INIT_IO_SURF_n(CSURF_FILETYPE,'FULL  ','SURF  ','READ ')
 CALL READ_SURF(CSURF_FILETYPE,'DIM_FULL  ',INI,  IRESP)
 CALL READ_SURF(CSURF_FILETYPE,'DTCUR     ',TTIME,  IRESP)
 CALL END_IO_SURF_n(CSURF_FILETYPE)

NINDX2 = INI
ALLOCATE(NWORK(INI))
ALLOCATE(XWORK(INI))
ALLOCATE(XWORK2(INI,10))
ALLOCATE(XWORK3(INI,10,10))
IF (NRANK==NPIO) THEN
  ALLOCATE(NWORK_FULL(INI))
  ALLOCATE(XWORK_FULL(INI))
  ALLOCATE(XWORK2_FULL(INI,10))
ELSE
  ALLOCATE(NWORK_FULL(0))
  ALLOCATE(XWORK_FULL(0))
  ALLOCATE(XWORK2_FULL(0,0))
ENDIF

IYEAR  = TTIME%TDATE%YEAR
IMONTH = TTIME%TDATE%MONTH
IDAY   = TTIME%TDATE%DAY
ZTIME  = TTIME%TIME
IF (ZTIME > NDAYSEC) ZTIME = ZTIME - NDAYSEC

KSW=0
KSV=0
ALLOCATE(CSV(KSV))
ALLOCATE(XCO2(INI))
ALLOCATE(XRHOA(INI))
ALLOCATE(XZENITH(INI))
ALLOCATE(XAZIM(INI))
ALLOCATE(XSW_BANDS(KSW))
ALLOCATE(XDIR_ALB(INI,KSW))
ALLOCATE(XSCA_ALB(INI,KSW))
ALLOCATE(XEMIS(INI))
ALLOCATE(XTSRAD(INI))

WRITE(*,*) "INITIALIZING SURFEX..."
! Initialize the SURFEX interface
 CALL INIT_SURF_ATM_n(CSURF_FILETYPE,YINIT, LLAND_USE,             &
                       INI, kSV, KSW,                             &
                       CSV,XCO2,XRHOA,                            &
                       XZENITH,XAZIM,XSW_BANDS,XDIR_ALB,XSCA_ALB, &
                       XEMIS,XTSRAD,                              &
                       IYEAR, IMONTH, IDAY, ZTIME,                &
                       YATMFILE, YATMFILETYPE, YTEST              )

! Initialyse the SURFACE-TRIP interface
 CALL INIT_SURF_TRIP_n(CSURF_FILETYPE,INI,KSW,LRESTART,IYEAR,IMONTH,&
                       ZDURATION,ITRIP_MONTH,ITRIP_COUNT,XZENITH,  &
                       XSW_BANDS,XEMIS,XTSRAD,XDIR_ALB,XSCA_ALB    )

! For EKF we only need CANARI FA fields for the analysis
! To save time we can do LPRT LSIM and LBEV without needing CANARI results
IF (( CASSIM_ISBA == 'EKF  ' ) .AND. ( ( LPRT ) .OR. ( LSIM ) .OR. ( LBEV ) )) THEN
  ALLOCATE(PT2M(NDIM_NATURE))
  ALLOCATE(PHU2M(NDIM_NATURE))
  PT2M=999.0
  PHU2M=999.0
  CALL ASSIM_NATURE_ISBA_EKF(CSURF_FILETYPE, NDIM_NATURE,&
                             PT2M,           PHU2M,      &
                             YTEST )
  DEALLOCATE(PT2M)
  DEALLOCATE(PHU2M)
ELSE
  WRITE(*,*) "READING input files..."
  ! Normal reading of needed FA fields
  ALLOCATE(PCON_RAIN(INI))
  ALLOCATE(PSTRAT_RAIN(INI))
  ALLOCATE(PCON_SNOW(INI))
  ALLOCATE(PSTRAT_SNOW(INI))
  ALLOCATE(PCON_GRAUPEL(INI))
  ALLOCATE(PCLOUDS(INI))
  ALLOCATE(PLSM(INI))
  ALLOCATE(PEVAPTR(INI))
  ALLOCATE(PEVAP(INI))
  ALLOCATE(PTSC(INI))
  ALLOCATE(PSWEC(INI))
  ALLOCATE(PTS(INI))
  ALLOCATE(PT2M(INI))
  ALLOCATE(PHU2M(INI))
  ALLOCATE(PSWE(INI))

  !  Read atmospheric forecast fields from FA files 
#ifdef FA
  CFILEIN_FA = 'FG_OI_MAIN'
#endif
  !  Open FA file (LAM version with extension zone)
  CALL INIT_IO_SURF_n(YPROGRAM2,'EXTZON','SURF  ','READ ')

  !  Read model forecast quantities
  IF (LAROME) THEN
    CALL READ_SURF(YPROGRAM2,'SURFACCPLUIE',    PCON_RAIN    ,IRESP)
    CALL READ_SURF(YPROGRAM2,'SURFACCNEIGE',    PCON_SNOW    ,IRESP)
    CALL READ_SURF(YPROGRAM2,'SURFACCGRAUPEL',  PCON_GRAUPEL ,IRESP)
  ELSE
    CALL READ_SURF(YPROGRAM2,'SURFPREC.EAU.CON',PCON_RAIN  ,IRESP)
    CALL READ_SURF(YPROGRAM2,'SURFPREC.EAU.GEC',PSTRAT_RAIN  ,IRESP)
    CALL READ_SURF(YPROGRAM2,'SURFPREC.NEI.CON',PCON_SNOW  ,IRESP)
    CALL READ_SURF(YPROGRAM2,'SURFPREC.NEI.GEC',PSTRAT_SNOW  ,IRESP)
  ENDIF
  CALL READ_SURF(YPROGRAM2,'ATMONEBUL.BASSE ',PCLOUDS,IRESP)
  CALL READ_SURF(YPROGRAM2,'SURFIND.TERREMER',PLSM   ,IRESP)
  CALL READ_SURF(YPROGRAM2,'SURFFLU.LAT.MEVA',PEVAP  ,IRESP) ! accumulated fluxes (not available in LFI)
  IF (.NOT.LALADSURF) THEN
    CALL READ_SURF(YPROGRAM2,'SURFXEVAPOTRANSP',PEVAPTR,IRESP) ! not in ALADIN SURFEX
  ELSE
    PEVAPTR(:) = 0.0
  ENDIF

  !  Close FA file
  CALL END_IO_SURF_n(YPROGRAM2)
  CALL IO_BUFF_CLEAN_n
  WRITE(*,*)'READ FG_OI_MAIN OK'

  !  Define FA file name for CANARI analysis
#ifdef FA
  CFILEIN_FA = 'CANARI'        ! input CANARI analysis
#endif

  !  Open FA file 
  CALL INIT_IO_SURF_n(YPROGRAM2,'EXTZON','SURF  ','READ ')

  !  Read CANARI analysis
  CALL READ_SURF(YPROGRAM2,'CLSTEMPERATURE  ',PT2M ,IRESP)
  CALL READ_SURF(YPROGRAM2,'CLSHUMI.RELATIVE',PHU2M,IRESP)
  CALL READ_SURF(YPROGRAM2,'SURFTEMPERATURE ',PTS  ,IRESP)
  CALL READ_SURF(YPROGRAM2,'SURFRESERV.NEIGE',PSWE ,IRESP)


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
  CALL READ_SURF(YPROGRAM2,'SURFRESERV.NEIGE',PSWEC  ,IRESP)
  CALL READ_SURF(YPROGRAM2,'SURFTEMPERATURE',PTSC ,IRESP)
  
  !  Close climatology file
  CALL END_IO_SURF_n(YPROGRAM2)
  CALL IO_BUFF_CLEAN_n
  WRITE(*,*) 'READ CLIMATOLOGY OK'

  WRITE(*,*) 'PERFORMIMG OFFLINE SURFEX DATA ASSIMILATION...'
  CALL ASSIM_SURF_ATM_n(CSURF_FILETYPE,INI,                              &
                        PCON_RAIN,  PSTRAT_RAIN, PCON_SNOW, PSTRAT_SNOW, &
                        PCLOUDS,    PLSM,        PEVAPTR,   PEVAP,       &
                        PSWEC,      PTSC,                                &
                        PTS,        PT2M,        PHU2M,     PSWE,        &
                        YTEST )


  DEALLOCATE(PCON_RAIN)
  DEALLOCATE(PSTRAT_RAIN)
  DEALLOCATE(PCON_SNOW)
  DEALLOCATE(PSTRAT_SNOW)
  DEALLOCATE(PCON_GRAUPEL)
  DEALLOCATE(PCLOUDS)
  DEALLOCATE(PLSM)
  DEALLOCATE(PEVAPTR)
  DEALLOCATE(PEVAP)
  DEALLOCATE(PTSC)
  DEALLOCATE(PSWEC)
  DEALLOCATE(PTS)
  DEALLOCATE(PT2M)
  DEALLOCATE(PHU2M)
  DEALLOCATE(PSWE)

ENDIF
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
IF (LHOOK) CALL DR_HOOK('SODA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END PROGRAM SODA
