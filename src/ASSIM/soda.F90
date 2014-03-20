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
!! 03/2014 E. Martin change indices names in OMP module according to GMAP changes
!----------------------------------------------------------------------------
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO
USE MODD_SURFEX_OMP, ONLY : NINDX2SFX, NWORK, XWORK, XWORK2, XWORK3, &
                            NWORK_FULL, XWORK_FULL, XWORK2_FULL
!
USE MODD_TYPE_DATE_SURF, ONLY : DATE_TIME
USE MODD_SURF_PAR,       ONLY : XUNDEF,NUNDEF
USE MODD_ASSIM,          ONLY : LASSIM,LAROME,LALADSURF,CASSIM_ISBA,LREAD_SST_FROM_FILE,&
                              & NVAR,XF,YF_PATCH,NOBSTYPE,XAT2M_ISBA,&
                              & XAHU2M_ISBA,XVAR,XOBS
USE MODD_ISBA_n,     ONLY : NPATCH,XWG,XTG
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
USE MODD_SURF_ATM_n,     ONLY : NSIZE_NATURE,NSIZE_FULL
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
                                CFILEPGD_LFI,CFILE_LFI,CLUOUT_LFI,CFILEOUT_LFI 
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
USE MODI_WRITE_SURF_ATM_n
USE MODI_WRITE_DIAG_SURF_ATM_n
USE MODI_ASSIM_READ_SST_FROM_FILE


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
CHARACTER(LEN=28)                :: CDSURFEX_FNAME_SAVE
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
REAL,ALLOCATABLE, DIMENSION(:)   :: PSST                ! SST from external file
REAL,ALLOCATABLE, DIMENSION(:)   :: PSIC                ! SIC from external file
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
INTEGER                          :: IIVAR,NTIMES
CHARACTER                        :: CVAR
INTEGER                          :: IVAR_COUNT
INTEGER                          :: IOBS

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
  CFILEOUT_LFI     = CPREPFILE
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

NINDX2SFX = INI
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

! Indicate that zenith and azimuth angles are not initialized

XZENITH = XUNDEF
XAZIM   = XUNDEF

IF ( CASSIM_ISBA == 'EKF  ' ) THEN
 ! Has to do initialization for all the perturbations + 
 ! control + the real run at last
 NTIMES=NVAR+2
ELSE
 NTIMES=1
ENDIF

WRITE(*,*) "INITIALIZING SURFEX..."

CDSURFEX_FNAME_SAVE=CFILEIN_LFI
DO IVAR_COUNT=1,NTIMES

  ! If we have more than one initialization to do
  IF ( NTIMES > 1 ) THEN
    ! For last initialization, we must re-do the first.
    ! Could be avoided by introducing knowlegde of LASSIM on this level
    IF ( IVAR_COUNT /= NTIMES ) THEN
      WRITE(CVAR,'(I1.1)') IVAR_COUNT-1
#ifdef LFI
      CFILEIN_LFI="PREP_EKF_PERT"//CVAR
#endif
    ELSE
#ifdef LFI
      CFILEIN_LFI=CDSURFEX_FNAME_SAVE
#endif
    ENDIF
  ELSE
#ifdef LFI
    CFILEIN_LFI=CDSURFEX_FNAME_SAVE
#endif
  ENDIF
#ifdef LFI
  CFILE_LFI        = CFILEIN_LFI
  CFILEIN_LFI_SAVE = CFILEIN_LFI
#endif

  ! Initialize the SURFEX interface
  CALL INIT_SURF_ATM_n(CSURF_FILETYPE,YINIT, LLAND_USE,             &
                         INI, kSV, KSW,                             &
                         CSV,XCO2,XRHOA,                            &
                         XZENITH,XAZIM,XSW_BANDS,XDIR_ALB,XSCA_ALB, &
                         XEMIS,XTSRAD,                              &
                         IYEAR, IMONTH, IDAY, ZTIME,                &
                         YATMFILE, YATMFILETYPE, YTEST              )

  IF ( CASSIM_ISBA == 'EKF  ' ) THEN
    IF ( IVAR_COUNT == 1 ) THEN
      ALLOCATE(XF(NSIZE_NATURE,NPATCH,NVAR+1,NVAR))
      ALLOCATE(YF_PATCH(NSIZE_NATURE,NPATCH,NVAR+1,NOBSTYPE))
    ENDIF
    IF (( IVAR_COUNT > 0 ) .AND. (IVAR_COUNT < NTIMES ))  THEN
      ! Set the global state values for this control value
      DO IOBS = 1,NOBSTYPE
        SELECT CASE (TRIM(XOBS(IOBS)))
          CASE("T2M")
            YF_PATCH(:,:,IVAR_COUNT,IOBS) = XAT2M_ISBA(:,:)
          CASE("HU2M")
            YF_PATCH(:,:,IVAR_COUNT,IOBS) = XAHU2M_ISBA(:,:)
          CASE("WG1")
            YF_PATCH(:,:,IVAR_COUNT,IOBS) = XWG(:,1,:)
          CASE DEFAULT
            CALL ABOR1_SFX("Mapping of "//XOBS(IOBS)//" is not defined in AROINI_SURFC!")
        END SELECT
      ENDDO

      ! Prognostic fields for assimilation (Control vector)
      DO IIVAR = 1,NVAR
        SELECT CASE (TRIM(XVAR(IIVAR)))
          CASE("TG1")
            XF(:,:,IVAR_COUNT,IIVAR) = XTG(:,1,:)
          CASE("TG2")
            XF(:,:,IVAR_COUNT,IIVAR) = XTG(:,2,:)
          CASE("WG1")
            XF(:,:,IVAR_COUNT,IIVAR) = XWG(:,1,:)
          CASE("WG2")
            XF(:,:,IVAR_COUNT,IIVAR) = XWG(:,2,:)
          CASE DEFAULT
            CALL ABOR1_SFX("Mapping of "//TRIM(XVAR(IIVAR))//" is not defined in AROINI_SURFC!")
        END SELECT
      ENDDO
    ENDIF
  ENDIF
ENDDO
! Initialyse the SURFACE-TRIP interface
CALL INIT_SURF_TRIP_n(CSURF_FILETYPE,INI,KSW,LRESTART,IYEAR,IMONTH,&
                       ZDURATION,ITRIP_MONTH,ITRIP_COUNT,XZENITH,  &
                       XSW_BANDS,XEMIS,XTSRAD,XDIR_ALB,XSCA_ALB    )

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
ALLOCATE(PSST(INI))
ALLOCATE(PSIC(INI))

!  Read atmospheric forecast fields from FA files 
#ifdef FA
CFILEIN_FA = 'FG_OI_MAIN'
#endif
!  Open FA file (LAM version with extension zone)
CALL INIT_IO_SURF_n(YPROGRAM2,'EXTZON','SURF  ','READ ')

!  Read model forecast quantities
IF (LAROME) THEN
  CALL READ_SURF(YPROGRAM2,'SURFACCPLUIE',    PCON_RAIN    ,IRESP)
  PSTRAT_RAIN(:) = 0.0
  CALL READ_SURF(YPROGRAM2,'SURFACCNEIGE',    PCON_SNOW    ,IRESP)
  PSTRAT_SNOW(:) = 0.0
  CALL READ_SURF(YPROGRAM2,'SURFACCGRAUPEL',  PCON_GRAUPEL ,IRESP)
  ! So far graupel has not been used
  !PCON_SNOW=PCON_SNOW+PCON_GRAUPEL
ELSE
  CALL READ_SURF(YPROGRAM2,'SURFPREC.EAU.CON',PCON_RAIN    ,IRESP)
  CALL READ_SURF(YPROGRAM2,'SURFPREC.EAU.GEC',PSTRAT_RAIN  ,IRESP)
  CALL READ_SURF(YPROGRAM2,'SURFPREC.NEI.CON',PCON_SNOW    ,IRESP)
  CALL READ_SURF(YPROGRAM2,'SURFPREC.NEI.GEC',PSTRAT_SNOW  ,IRESP)
  PCON_GRAUPEL(:) = 0.0
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

IF ( LREAD_SST_FROM_FILE ) CALL ASSIM_READ_SST_FROM_FILE(CSURF_FILETYPE,NSIZE_FULL,PLSM,PSST,PSIC,YTEST)

IF ( .NOT. LASSIM ) CALL ABOR1_SFX("YOU CAN'T RUN SODA WITHOUT SETTING LASSIM=.TRUE. IN THE ASSIM NAMELIST")

WRITE(*,*) 'PERFORMIMG OFFLINE SURFEX DATA ASSIMILATION...'
CALL ASSIM_SURF_ATM_n(CSURF_FILETYPE,INI,                              &
                      PCON_RAIN,  PSTRAT_RAIN, PCON_SNOW, PSTRAT_SNOW, &
                      PCLOUDS,    PLSM,        PEVAPTR,   PEVAP,       &
                      PSWEC,      PTSC,                                &
                      PTS,        PT2M,        PHU2M,     PSWE,        &
                      PSST,       PSIC,                                &
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
DEALLOCATE(PSST)
DEALLOCATE(PSIC)

! Store results from assimilation
CALL WRITE_SURF_ATM_n(CSURF_FILETYPE,'ALL',LLAND_USE)
CALL WRITE_DIAG_SURF_ATM_n(CSURF_FILETYPE,'ALL')

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
