! -------------------------------------------------
PROGRAM OFFLINE
!
! -------------------------------------------------
! Driver structure
! ----------------
! 1. Initializations
! 2. Temporal loops
!   2.a Read forcing
!   2.b Interpolate forcing in time
!   2.c Run surface
!   2.d Write prognostics and diagnostics variables
!
! -------------------------------------------------
USE MODD_FORC_ATM,  ONLY: CSV         ,&! name of all scalar variables
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
                            XDIR_SW     ,&! direct  solar radiation (on horizontal surf.)
                            XSCA_SW     ,&! diffuse solar radiation (on horizontal surf.)
                            XSW_BANDS   ,&! mean wavelength of each shortwave band (m)
                            XZENITH     ,&! zenithal angle       (radian from the vertical)
                            XZENITH2    ,&! zenithal angle       (radian from the vertical)
                            XAZIM       ,&! azimuthal angle      (radian from North, clockwise)
                            XLW         ,&! longwave radiation (on horizontal surf.)
                            XPS         ,&! pressure at atmospheric model surface (Pa)
                            XPA         ,&! pressure at forcing level             (Pa)
                            XRHOA       ,&! density at forcing level              (kg/m3)
                            XCO2        ,&! CO2 concentration in the air          (kg/m3)
                            XSNOW       ,&! snow precipitation                    (kg/m2/s)
                            XRAIN       ,&! liquid precipitation                  (kg/m2/s)
                            XSFTH       ,&! flux of heat                          (W/m2)
                            XSFTQ       ,&! flux of water vapor                   (kg/m2/s)
                            XSFU        ,&! zonal momentum flux                   (m/s)
                            XSFV        ,&! meridian momentum flux                (m/s)
                            XSFCO2      ,&! flux of CO2                           (kg/m2/s)
                            XSFTS       ,&! flux of scalar var.                   (kg/m2/s)
                            XPEW_A_COEF ,&! implicit coefficients
                            XPEW_B_COEF ,&! needed if HCOUPLING='I'
                            XPET_A_COEF ,&
                            XPEQ_A_COEF ,&
                            XPET_B_COEF ,&
                            XPEQ_B_COEF  
!
USE MODI_OL_READ_ATM_CONF
USE MODI_COMPARE_OROGRAPHY
USE MODI_OL_READ_ATM
USE MODI_OL_ALLOC_ATM
USE MODI_OL_TIME_INTERP_ATM
USE MODI_SUNPOS
USE MODI_OPEN_NAMELIST
USE MODI_CLOSE_NAMELIST
USE MODI_READ_ALL_NAMELISTS
USE MODI_TEST_NAM_VAR_SURF
USE MODI_GET_SURF_VAR_n
USE MODI_FLAG_UPDATE
USE MODI_FLAG_DIAG_UPDATE
USE MODI_CLOSE_FILEOUT_OL
USE MODI_OPEN_FILEIN_OL
USE MODI_CLOSE_FILEIN_OL
USE MODI_ADD_FORECAST_TO_DATE_SURF
USE MODI_ALLOC_SURFEX
USE MODI_DEALLOC_SURFEX
USE MODI_DIAG_SURF_ATM_n
USE MODI_GET_LUOUT
!
USE MODD_SURF_CONF,  ONLY : CPROGNAME
USE MODD_CSTS,       ONLY : XPI, XDAY, XRV, XRD, XG
USE MODD_IO_SURF_ASC,ONLY : CFILEIN,CFILEIN_SAVE,CFILEOUT,CFILEPGD
USE MODD_SURF_PAR
USE MODD_IO_SURF_FA, ONLY : CFILEIN_FA, CFILEIN_FA_SAVE,       &
                            CFILEOUT_FA, NUNIT_FA, CDNOMC,     &
                            IVERBFA, LFANOCOMPACT, CFILEPGD_FA  
USE MODD_IO_SURF_LFI,ONLY : CFILEIN_LFI, CFILEIN_LFI_SAVE, CLUOUT_LFI, CFILEOUT_LFI, &
                            LMNH_COMPATIBLE, CFILEPGD_LFI  
USE MODD_IO_SURF_OL, ONLY : XSTART, XCOUNT, XSTRIDE,           &
                              LDEFINED_NATURE, LDEFINED_SEA,    &
                              LDEFINED_WATER,  LDEFINED_TOWN,   &
                              LDEFINED_SURF_ATM, LPARTW,        &
                              XSTARTW, XCOUNTW, LTIME_WRITTEN,  &
                              NSTEP_OUTPUT  
!
USE MODD_WRITE_BIN,  ONLY : NWRITE
!
USE MODD_SURF_ATM,   ONLY : LCPL_ESM
!
USE MODE_POS_SURF
!
USE MODN_IO_OFFLINE
!
! --------------------------------------------------------------------------------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
USE MODI_COUPLING_SURF_ATM_n
USE MODI_COUPLING_SURF_TRIP_n
USE MODI_GOTO_SURFEX
USE MODI_GOTO_TRIP
USE MODI_INIT_SURF_ATM_n
USE MODI_INIT_SURF_TRIP_n
USE MODI_IO_BUFF_CLEAN_n
USE MODI_OPEN_CLOSE_BIN_ASC_FORC
USE MODI_WRITE_DIAG_SURF_ATM_n
USE MODI_WRITE_HEADER_MNH
USE MODI_WRITE_SURF_ATM_n
USE MODI_INIT_SURF_LANDUSE_n
!
IMPLICIT NONE
!
!*      0.    declarations of local variables
!
CHARACTER(LEN=3), PARAMETER       :: YINIT     = 'ALL'
!
CHARACTER(LEN=28)                 :: YLUOUT    = 'LISTING_OFFLINE             '
!
INTEGER                           :: IYEAR               ! current year (UTC)
INTEGER                           :: IMONTH              ! current month (UTC)
INTEGER                           :: IDAY                ! current day (UTC)
REAL                              :: ZTIME               ! current time since start of the run (s)
REAL                              :: ZTIMEC              ! current duration since start of the run (s)
!
INTEGER                           :: IYEAR_OUT           ! output year name
INTEGER                           :: IMONTH_OUT          ! output month name
INTEGER                           :: IDAY_OUT            ! output day name
REAL                              :: ZTIME_OUT           ! output time since start of the run (s)
!
INTEGER, DIMENSION(11)  :: IDATEF
!
CHARACTER(LEN=28), PARAMETER      :: YATMFILE     = '                            '
CHARACTER(LEN=6),  PARAMETER      :: YATMFILETYPE = '      '
CHARACTER(LEN=2),  PARAMETER      :: YTEST        = 'OK'          ! must be equal to 'OK'
!
REAL, DIMENSION(:), POINTER       :: ZLAT                ! latitude                         (rad)
REAL, DIMENSION(:), POINTER       :: ZLON                ! longitude                        (rad)
REAL, DIMENSION(:), POINTER       :: ZZS_FORC            ! orography                        (m)  
REAL, DIMENSION(:), POINTER       :: ZZREF               ! Forcing level for T
REAL, DIMENSION(:), POINTER       :: ZUREF               ! Forcing level for U
!
REAL                              :: ZTSTEP              ! atmospheric time-step            (s)
!
INTEGER                           :: INI                 ! grid dimension
INTEGER                           :: JLOOP               ! loop counter
INTEGER                           :: ISCAL               ! Number of scalar species
INTEGER                           :: IBANDS              ! Number of radiative bands 
INTEGER                           :: INB_STEP_ATM        ! Number of atmospheric time-steps
INTEGER                           :: INB_ATM             ! Number of Isba time-steps 
                                                         ! within a forcing time-step
INTEGER                           :: ID_FORC             ! indice of forcing in the file
INTEGER                           :: INB_LINES           ! nb of lines to read in the forcing file
INTEGER                           :: IDMAX               ! nb of lines to read in the forcing file at last 
INTEGER                           :: JFORC_STEP          ! atmospheric loop index
INTEGER                           :: JSURF_STEP          ! isba loop index
INTEGER                           :: ICOUNT              ! day counter 
INTEGER                           :: ITRIP_COUNT         ! day counter
INTEGER                           :: ITRIP_MONTH
REAL                              :: ZDURATION           ! duration of run                     (s)
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZTA                 ! air temperature forcing               (K)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZQA                 ! air humidity forcing                  (kg/m3)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZWIND               ! wind speed                            (m/s)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZSCA_SW             ! diffuse solar radiation (on horizontal surf.)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZDIR_SW             ! direct  solar radiation (on horizontal surf.)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZLW                 ! longwave radiation (on horizontal surf.)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZSNOW               ! snow precipitation                    (kg/m2/s)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZRAIN               ! liquid precipitation                  (kg/m2/s)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZPS                 ! pressure at forcing level             (Pa)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCO2                ! CO2 concentration in the air          (kg/m3)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZDIR                ! wind direction
INTEGER                           :: ILUOUT              ! ascii output unit number
INTEGER                           :: ILUNAM              ! namelist unit number
INTEGER                           :: IRET                ! error return code
INTEGER                           :: INB 
CHARACTER(LEN=14)                 :: YTAG                
LOGICAL                           :: GFOUND              ! return logical when reading namelist
REAL, DIMENSION(:),   ALLOCATABLE :: ZSW                 ! total solar radiation (on horizontal surf.)
!
! Inquiry mode arrays:
!
REAL, DIMENSION(:), ALLOCATABLE   :: ZSEA, ZWATER, ZNATURE, ZTOWN
REAL, DIMENSION(:), ALLOCATABLE   :: ZT2M, ZQ2M
REAL, DIMENSION(:), ALLOCATABLE   :: ZZ0, ZZ0H, ZQS
REAL, DIMENSION(:), ALLOCATABLE   :: ZQS_SEA, ZQS_WATER, ZQS_NATURE, ZQS_TOWN
REAL, DIMENSION(:), ALLOCATABLE   :: ZPSNG, ZPSNV
REAL, DIMENSION(:), ALLOCATABLE   :: ZZ0EFF
REAL, DIMENSION(:), ALLOCATABLE   :: ZZS
INTEGER :: ISERIES
REAL(KIND=JPRB) :: ZHOOK_HANDLE
! --------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('OFFLINE',0,ZHOOK_HANDLE)
CALL ALLOC_SURFEX(1)
!
!*      0.    Open ascii file for writing
!
CLUOUT_LFI =  ADJUSTL(ADJUSTR(YLUOUT)//'.txt')
CALL GET_LUOUT('ASCII ',ILUOUT)
OPEN(UNIT=ILUOUT,FILE=ADJUSTL(ADJUSTR(YLUOUT)//'.txt'),FORM='FORMATTED',ACTION='WRITE')
!
!*      0.1   Open namelist
!
CALL OPEN_NAMELIST('ASCII ',ILUNAM,CNAMELIST)
!
CALL POSNAM(ILUNAM,'NAM_IO_OFFLINE',GFOUND,ILUOUT)
IF (GFOUND) READ (UNIT=ILUNAM,NML=NAM_IO_OFFLINE)
!
CALL TEST_NAM_VAR_SURF(ILUOUT,'CSURF_FILETYPE',CSURF_FILETYPE,'ASCII ','LFI   ','FA    ')
CALL TEST_NAM_VAR_SURF(ILUOUT,'CTIMESERIES_FILETYPE',CTIMESERIES_FILETYPE,'NETCDF','TEXTE ','BINARY',&
                                                                            'ASCII ','LFI   ','FA    ',&
                                                                            'NONE  ','OFFLIN')  
CALL TEST_NAM_VAR_SURF(ILUOUT,'CFORCING_FILETYPE',CFORCING_FILETYPE,'NETCDF','ASCII ','BINARY')
IF (CTIMESERIES_FILETYPE=='NETCDF') CTIMESERIES_FILETYPE='OFFLIN'
!
CFILEPGD     = ADJUSTL(ADJUSTR(CPGDFILE)//'.txt')
CFILEPGD_LFI = CPGDFILE
CFILEPGD_FA  = ADJUSTL(ADJUSTR(CPGDFILE)//'.fa')
!
CFILEIN     = ADJUSTL(ADJUSTR(CPREPFILE)//'.txt')
CFILEIN_LFI = CPREPFILE
CFILEIN_FA  = ADJUSTL(ADJUSTR(CPREPFILE)//'.fa')
!
CFILEIN_SAVE     = CFILEIN
CFILEIN_LFI_SAVE = CFILEIN_LFI
CFILEIN_FA_SAVE  = CFILEIN_FA
!
CALL CLOSE_NAMELIST('ASCII ',ILUNAM)
!
CALL READ_ALL_NAMELISTS(CSURF_FILETYPE,'ALL',.FALSE.)
!
CALL GOTO_SURFEX(1,.TRUE.)
CALL GOTO_TRIP(1,.TRUE.)
!
!*      0.2   Assume FA filetype consistency 
!
CPROGNAME = CSURF_FILETYPE
!
! --------------------------------------------------------------------------------------
!
!*      1.    Initializations
!
LCPL_ESM = .FALSE.
!
!       netcdf file handling
!
XSTART            = NUNDEF
XSTRIDE           = NUNDEF
XCOUNT            = NUNDEF
XSTARTW           = 0
XCOUNTW           = 1
LPARTW            = .TRUE.
LDEFINED_SURF_ATM = .FALSE.
LDEFINED_NATURE   = .FALSE.
LDEFINED_TOWN     = .FALSE.
LDEFINED_WATER    = .FALSE.
LDEFINED_SEA      = .FALSE.
!
!       ascii or FA file handling
!
!
IF (CFORCING_FILETYPE=='ASCII ' .OR. CFORCING_FILETYPE=='BINARY') CALL OPEN_CLOSE_BIN_ASC_FORC('CONF ',CFORCING_FILETYPE,0,'R')
IF (CFORCING_FILETYPE=='NETCDF') CALL OPEN_FILEIN_OL
!
!       configuration of run
!
CALL OL_READ_ATM_CONF(CSURF_FILETYPE, CFORCING_FILETYPE,            &
                        ZDURATION, ZTSTEP, INI, IYEAR, IMONTH, IDAY,  &
                        ZTIME, ZLAT, ZLON, ZZS_FORC, ZZREF, ZUREF     )  
!
!*     time steps coherence check 
!
IF ( (MOD(XTSTEP_OUTPUT,ZTSTEP)*MOD(ZTSTEP,XTSTEP_OUTPUT) /= 0) .OR. (MOD(ZTSTEP,XTSTEP_SURF) /= 0) ) THEN
   WRITE(ILUOUT,*)' FORCING  AND OUTPUT/SURFACE TIME STEP SHOULD BE MULTIPLE', &
     NINT(ZTSTEP),NINT(XTSTEP_OUTPUT),NINT(XTSTEP_SURF)    
   CALL ABOR1_SFX('OFFLINE: FORCING  AND OUTPUT/SURFACE TIME STEP SHOULD BE MULTIPLE')
ENDIF
!
IF ( ZTIME /= 0. .AND. MOD(ZTIME,XTSTEP_SURF) /= 0  ) THEN
   WRITE(ILUOUT,*)' INITIAL AND SURFACE TIME STEP SHOULD BE MULTIPLE', &
   NINT(ZTIME),NINT(XTSTEP_SURF)  
   CALL ABOR1_SFX('OFFLINE: INITIAL AND SURFACE TIME STEP SHOULD BE MULTIPLE')
ENDIF
!
INB_STEP_ATM  = INT(ZDURATION / ZTSTEP)
INB_ATM       = INT(ZTSTEP / XTSTEP_SURF)
NSTEP_OUTPUT  = INT(ZDURATION / XTSTEP_OUTPUT)
!
!* opens forcing files (if ASCII or BINARY)
!
IF (CFORCING_FILETYPE=='ASCII ' .OR. CFORCING_FILETYPE=='BINARY') CALL OPEN_CLOSE_BIN_ASC_FORC('OPEN ',CFORCING_FILETYPE,INI,'R')
!
!
!       allocation of variables
!
IBANDS = 1
ISCAL  = 1
!
CALL OL_ALLOC_ATM(INI,IBANDS,ISCAL)
!
XZS   = ZZS_FORC
XZREF = ZZREF
XUREF = ZUREF
!
!       compare orography
!
CALL COMPARE_OROGRAPHY (CSURF_FILETYPE, LSET_FORC_ZS, 200.)
!
!       miscellaneous initialization
!
ICOUNT = 0
ZTIMEC = 0.
!
CALL SUNPOS(IYEAR, IMONTH, IDAY, ZTIME, ZLON, ZLAT, XTSUN, XZENITH, XAZIM)
!
!number of lines read in forcing files
INB_LINES=1
IF (NB_READ_FORC.EQ.1) THEN
  INB_LINES=INB_STEP_ATM
ELSEIF (NB_READ_FORC.NE.0) THEN
  !to be sure the number of readings will be NB_READ_FORC as a maximum
  INB_LINES=CEILING(1.*(INB_STEP_ATM+1)/NB_READ_FORC)
ENDIF
!number of lines to be read effectively
IDMAX=INB_LINES+1
!effective number of readings of the forcing files
NB_READ_FORC=CEILING(1.*(INB_STEP_ATM+1)/INB_LINES)
!
!       allocate local atmospheric variables
!
IF (.NOT.ALLOCATED(ZTA)) ALLOCATE(ZTA    (INB_LINES+1,INI)) 
IF (.NOT.ALLOCATED(ZQA))ALLOCATE(ZQA    (INB_LINES+1,INI))
IF (.NOT.ALLOCATED(ZWIND))ALLOCATE(ZWIND  (INB_LINES+1,INI))
IF (.NOT.ALLOCATED(ZDIR_SW))ALLOCATE(ZDIR_SW(INB_LINES+1,INI))
IF (.NOT.ALLOCATED(ZSCA_SW))ALLOCATE(ZSCA_SW(INB_LINES+1,INI))
IF (.NOT.ALLOCATED(ZLW))ALLOCATE(ZLW    (INB_LINES+1,INI))
IF (.NOT.ALLOCATED(ZSNOW))ALLOCATE(ZSNOW  (INB_LINES+1,INI))
IF (.NOT.ALLOCATED(ZRAIN))ALLOCATE(ZRAIN  (INB_LINES+1,INI))
IF (.NOT.ALLOCATED(ZPS))ALLOCATE(ZPS    (INB_LINES+1,INI))
IF (.NOT.ALLOCATED(ZCO2))ALLOCATE(ZCO2   (INB_LINES+1,INI))
IF (.NOT.ALLOCATED(ZDIR))ALLOCATE(ZDIR   (INB_LINES+1,INI))
!
IF (.NOT.ALLOCATED(ZSW))ALLOCATE(ZSW    (INI))
!
!      computes initial air co2 concentration and  density
!
CALL OL_READ_ATM(CSURF_FILETYPE, CFORCING_FILETYPE, 1,             &
                   ZTA,ZQA,ZWIND,ZDIR_SW,ZSCA_SW,ZLW,ZSNOW,ZRAIN,ZPS,&
                   ZCO2,ZDIR,LLIMIT_QAIR                           )  
!     
XCO2(:)  = ZCO2(1,:)
XRHOA (:) = ZPS(1,:) / (XRD * ZTA(1,:) * ( 1.+((XRV/XRD)-1.)*ZQA(1,:) ) + XG * XZREF )
!                 
!       surface Initialisation     
!
CALL IO_BUFF_CLEAN_n
CALL INIT_SURF_ATM_n(CSURF_FILETYPE,YINIT, LLAND_USE,             &
                       INI, ISCAL, IBANDS,                        &
                       CSV,XCO2,XRHOA,                            &
                       XZENITH,XAZIM,XSW_BANDS,XDIR_ALB,XSCA_ALB, &
                       XEMIS,XTSRAD,                              &
                       IYEAR, IMONTH, IDAY, ZTIME,                &
                       YATMFILE, YATMFILETYPE, YTEST              )
!
!   Land use or/and vegetation dynamic
!                  
CALL INIT_SURF_LANDUSE_n(CSURF_FILETYPE,YINIT,LLAND_USE,          &
                       INI, ISCAL, IBANDS,                        &
                       CSV,XCO2,XRHOA,                            &
                       XZENITH,XAZIM,XSW_BANDS,XDIR_ALB,XSCA_ALB, &
                       XEMIS,XTSRAD,                              &
                       IYEAR, IMONTH, IDAY, ZTIME,                &
                       YATMFILE, YATMFILETYPE, YTEST              )
!
! Initialyse the SURFACE-TRIP interface
!
CALL INIT_SURF_TRIP_n(CSURF_FILETYPE,INI,IBANDS,LRESTART,IYEAR,IMONTH,&
                        ZDURATION,ITRIP_MONTH,ITRIP_COUNT,XZENITH,      &
                        XSW_BANDS,XEMIS,XTSRAD,XDIR_ALB,XSCA_ALB        )  
!
!
! --------------------------------------------------------------------------------------
NWRITE = 0
!
!*      2.    Temporal loops
!
DO JFORC_STEP=1,INB_STEP_ATM
   !
   ! read Forcing
   !
   !indice of forcing line in forcing arrays
   ID_FORC=JFORC_STEP-INT(JFORC_STEP/INB_LINES)*INB_LINES
   IF (ID_FORC==0) ID_FORC=INB_LINES
   !new forcings to read
   IF (ID_FORC==1 .AND. JFORC_STEP.NE.1) THEN
     !if last part of forcing, the last point has to be adjusted on the end of
     !files
     IF (JFORC_STEP/INB_LINES==NB_READ_FORC-1) THEN 
       IDMAX=INB_STEP_ATM-JFORC_STEP+1+1
       !for ascii and binary forcing files
       ZTA(IDMAX,:)=ZTA(SIZE(ZTA,1),:)
       ZQA(IDMAX,:)=ZQA(SIZE(ZTA,1),:)
       ZWIND(IDMAX,:)=ZWIND(SIZE(ZTA,1),:)
       ZDIR_SW(IDMAX,:)=ZDIR_SW(SIZE(ZTA,1),:)
       ZSCA_SW(IDMAX,:)=ZSCA_SW(SIZE(ZTA,1),:)
       ZLW(IDMAX,:)=ZLW(SIZE(ZTA,1),:)
       ZSNOW(IDMAX,:)=ZSNOW(SIZE(ZTA,1),:)
       ZRAIN(IDMAX,:)=ZRAIN(SIZE(ZTA,1),:)
       ZPS(IDMAX,:)=ZPS(SIZE(ZTA,1),:)
       ZCO2(IDMAX,:)=ZCO2(SIZE(ZTA,1),:)
       ZDIR(IDMAX,:)=ZDIR(SIZE(ZTA,1),:)
     ENDIF
     CALL OL_READ_ATM(CSURF_FILETYPE, CFORCING_FILETYPE, JFORC_STEP,    &
                      ZTA(1:IDMAX,:),ZQA(1:IDMAX,:),ZWIND(1:IDMAX,:), &
                      ZDIR_SW(1:IDMAX,:),ZSCA_SW(1:IDMAX,:),ZLW(1:IDMAX,:), &
                      ZSNOW(1:IDMAX,:),ZRAIN(1:IDMAX,:),ZPS(1:IDMAX,:),&
                      ZCO2(1:IDMAX,:),ZDIR(1:IDMAX,:),LLIMIT_QAIR         )
   ENDIF  
   !  
   DO JSURF_STEP=1,INB_ATM
       !
       ! time interpolation of the forcing
       !
       CALL SUNPOS(IYEAR, IMONTH, IDAY, ZTIME, ZLON, ZLAT, XTSUN, XZENITH, XAZIM)
       CALL SUNPOS(IYEAR, IMONTH, IDAY, ZTIME+XTSTEP_SURF, ZLON, ZLAT, XTSUN, XZENITH2, XAZIM)
       !interpolation between beginning and end of current forcing time step
       CALL OL_TIME_INTERP_ATM(JSURF_STEP,INB_ATM,                               &
                               ZTA(ID_FORC:ID_FORC+1,:),ZQA(ID_FORC:ID_FORC+1,:),&
                               ZWIND(ID_FORC:ID_FORC+1,:),ZDIR_SW(ID_FORC:ID_FORC+1,:),&
                               ZSCA_SW(ID_FORC:ID_FORC+1,:),ZLW(ID_FORC:ID_FORC+1,:),&
                               ZSNOW(ID_FORC:ID_FORC+1,:),ZRAIN(ID_FORC:ID_FORC+1,:),&
                               ZPS(ID_FORC:ID_FORC+1,:),ZCO2(ID_FORC:ID_FORC+1,:),&
                               ZDIR(ID_FORC:ID_FORC+1,:)                      )  
       !
       ! coherence between solar zenithal angle and radiation
       !
       ZSW(:) = 0.
       DO JLOOP=1,SIZE(XDIR_SW,2)
         ZSW(:) = ZSW(:) + XDIR_SW(:,JLOOP) + XSCA_SW(:,JLOOP)
       END DO
       WHERE (ZSW(:)>0.)
         XZENITH  = MIN (XZENITH ,XPI/2.-0.01)
         XZENITH2 = MIN (XZENITH2,XPI/2.-0.01)
       ELSEWHERE
         XZENITH  = MAX (XZENITH ,XPI/2.)
         XZENITH2 = MAX (XZENITH2,XPI/2.)
       END WHERE
       !
       ! updates time
       ZTIMEC= ZTIMEC+ XTSTEP_SURF
       !
       ! run Surface
       !
       CALL IO_BUFF_CLEAN_n
       CALL COUPLING_SURF_ATM_n(CSURF_FILETYPE, 'E', ZTIMEC,                   &
                XTSTEP_SURF, IYEAR, IMONTH, IDAY, ZTIME, INI, ISCAL, IBANDS,   &
                XTSUN, XZENITH, XZENITH2, XAZIM,                               &
                XZREF, XUREF, XZS, XU, XV, XQA, XTA, XRHOA, XSV, XCO2, CSV,    &
                XRAIN, XSNOW, XLW, XDIR_SW, XSCA_SW, XSW_BANDS, XPS, XPA,      &
                XSFTQ, XSFTH, XSFTS, XSFCO2, XSFU, XSFV,                       &
                XTSRAD, XDIR_ALB, XSCA_ALB, XEMIS,                             &
                XPEW_A_COEF, XPEW_B_COEF,                                      &
                XPET_A_COEF, XPEQ_A_COEF, XPET_B_COEF, XPEQ_B_COEF, YTEST      )  
       !
       CALL COUPLING_SURF_TRIP_n(CSURF_FILETYPE,INI,IBANDS,LRESTART,IYEAR,  &
                                   ITRIP_MONTH,ITRIP_COUNT,ZTIME+XTSTEP_SURF, &
                                   ZDURATION,XZENITH,XSW_BANDS,XEMIS,XTSRAD,  &
                                   XDIR_ALB,XSCA_ALB                          )  
       !
       ZTIME = ZTIME + XTSTEP_SURF
       CALL ADD_FORECAST_TO_DATE_SURF(IYEAR, IMONTH, IDAY, ZTIME)
       !
       ! ecrit Surface
       !
       !
       IF (MOD(ZTIMEC,XTSTEP_OUTPUT) == 0. .AND. CTIMESERIES_FILETYPE/='NONE  ') THEN
          !* name of the file
          IF (CTIMESERIES_FILETYPE=="ASCII " .OR. &
                CTIMESERIES_FILETYPE=="LFI   " .OR. &
                CTIMESERIES_FILETYPE=="FA    "      ) THEN  
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
            WRITE(YTAG,FMT='(I4.4,I2.2,I2.2,A1,I2.2,A1,I2.2)') IYEAR_OUT,IMONTH_OUT,IDAY_OUT,&
                   '_',INT(ZTIME_OUT/3600.),'h',NINT(ZTIME_OUT)/60-60*INT(ZTIME_OUT/3600.)  
            CFILEOUT    = ADJUSTL(ADJUSTR(CSURFFILE)//'.'//YTAG//'.txt')
            CFILEOUT_LFI= ADJUSTL(ADJUSTR(CSURFFILE)//'.'//YTAG)
            CFILEOUT_FA = ADJUSTL(ADJUSTR(CSURFFILE)//'.'//YTAG//'.fa')
            IF (CTIMESERIES_FILETYPE=='FA    ') THEN 
              IDATEF(1)= IYEAR
              IDATEF(2)= IMONTH
              IDATEF(3)= IDAY
              IDATEF(4)= NINT(ZTIME/3600.) 
              IDATEF(5)= NINT(ZTIME/60.) - IDATEF(4) * 60
              IDATEF(6)= NINT(ZTIME) - IDATEF(4) * 3600 - IDATEF(5) * 60
              IDATEF(7:11) = 0                      
              CALL FAITOU(IRET,NUNIT_FA,.TRUE.,CFILEOUT_FA,'NEW',.TRUE.,.FALSE.,IVERBFA,0,INB,CDNOMC) 
              CALL FANDAR(IRET,NUNIT_FA,IDATEF)
            END IF
          END IF
          !
          XSTARTW = XSTARTW + 1
          NWRITE  = NWRITE  + 1
          LTIME_WRITTEN(:)=.FALSE.          
          CALL IO_BUFF_CLEAN_n
          CALL WRITE_SURF_ATM_n(CTIMESERIES_FILETYPE,'ALL',LLAND_USE)
          CALL DIAG_SURF_ATM_n(CTIMESERIES_FILETYPE)
          CALL WRITE_DIAG_SURF_ATM_n(CTIMESERIES_FILETYPE,'ALL')
          !
          IF (CTIMESERIES_FILETYPE=='FA    ') THEN
            CALL FAIRME(IRET,NUNIT_FA,'UNKNOWN')
          END IF
          !* add informations in the file
          IF (CTIMESERIES_FILETYPE=='LFI   ' .AND. LMNH_COMPATIBLE) CALL WRITE_HEADER_MNH
       ENDIF
       !
   END DO
   !
   IF (LPRINT) THEN
      IF (MOD(ZTIMEC,XDAY) == 0.) THEN
         ICOUNT = ICOUNT + 1
         WRITE(*,'(A7,I4,A2,I4)')'  DAY :',ICOUNT,' /',INT(ZDURATION/XDAY)
      ENDIF
   ENDIF
   !
END DO
!
! --------------------------------------------------------------------------------------
!
!*    3.     write restart file
!            ------------------
!
IF ( LRESTART ) THEN
  !* name of the file
  CFILEOUT    = ADJUSTL(ADJUSTR(CSURFFILE)//'.txt')
  CFILEOUT_LFI= CSURFFILE
  CFILEOUT_FA = ADJUSTL(ADJUSTR(CSURFFILE)//'.fa')

  !* opens the file
  IF (CSURF_FILETYPE=='FA    ') THEN
    LFANOCOMPACT = .TRUE.
    IDATEF(1)=IYEAR
    IDATEF(2)=IMONTH
    IDATEF(3)=IDAY
    IDATEF(4)=NINT(ZTIME/3600.) 
    IDATEF(5)=NINT(ZTIME/60.) - IDATEF(4) * 60 
    IDATEF(6)=NINT(ZTIME) - IDATEF(4) * 3600 - IDATEF(5) * 6 
    IDATEF(7:11)=0  
    CALL FAITOU(IRET,NUNIT_FA,.TRUE.,CFILEOUT_FA,'NEW',.TRUE.,.FALSE.,IVERBFA,0,INB,CDNOMC)
    CALL FANDAR(IRET,NUNIT_FA,IDATEF)
  END IF
  !
  CALL FLAG_UPDATE(.FALSE.,.TRUE.,.FALSE.,.FALSE.)
  !
  !* writes into the file
  CALL IO_BUFF_CLEAN_n
  CALL WRITE_SURF_ATM_n(CSURF_FILETYPE,'ALL',LLAND_USE)
  IF(CSURF_FILETYPE/='FA    ')THEN
     CALL FLAG_DIAG_UPDATE(.FALSE.,.TRUE.,0,.FALSE.,.FALSE.,.FALSE.,&
                           .FALSE.,0,0,.FALSE.,.FALSE.,.FALSE.,.FALSE.,&
                           .FALSE.,.FALSE.,.FALSE.,.FALSE.)
     CALL WRITE_DIAG_SURF_ATM_n(CSURF_FILETYPE,'ALL')
  ENDIF
  !* closes the file
  IF (CSURF_FILETYPE=='FA    ') THEN
    CALL FAIRME(IRET,NUNIT_FA,'UNKNOWN')
  END IF
  !* add informations in the file
  IF (CSURF_FILETYPE=='LFI   ' .AND. LMNH_COMPATIBLE) CALL WRITE_HEADER_MNH

END IF
!
! --------------------------------------------------------------------------------------
!
!*    4.     inquiry mode
!            ------------
!
IF ( LINQUIRE ) THEN
!
   ALLOCATE( ZSEA       ( INI ) )
   ALLOCATE( ZWATER     ( INI ) )
   ALLOCATE( ZNATURE    ( INI ) )
   ALLOCATE( ZTOWN      ( INI ) )
   ALLOCATE( ZT2M       ( INI ) )
   ALLOCATE( ZQ2M       ( INI ) )
   ALLOCATE( ZZ0        ( INI ) )
   ALLOCATE( ZZ0H       ( INI ) )
   ALLOCATE( ZQS_SEA    ( INI ) )
   ALLOCATE( ZQS_WATER  ( INI ) )
   ALLOCATE( ZQS_NATURE ( INI ) )
   ALLOCATE( ZQS_TOWN   ( INI ) )
   ALLOCATE( ZQS        ( INI ) )
   ALLOCATE( ZPSNG      ( INI ) )
   ALLOCATE( ZPSNV      ( INI ) )
   ALLOCATE( ZZ0EFF     ( INI ) )
   ALLOCATE( ZZS        ( INI ) )
   !
   ISERIES = 0
   CALL GET_SURF_VAR_n(CSURF_FILETYPE,INI,ISERIES,PSEA=ZSEA,PWATER=ZWATER,PNATURE=ZNATURE,PTOWN=ZTOWN, &
                         PT2M=ZT2M,PQ2M=ZQ2M,PQS=ZQS,PZ0=ZZ0,PZ0H=ZZ0H,PZ0EFF=ZZ0EFF,PQS_SEA=ZQS_SEA,    &
                         PQS_WATER=ZQS_WATER,PQS_NATURE=ZQS_NATURE,PQS_TOWN=ZQS_TOWN,            &
                         PPSNG=ZPSNG,PPSNV=ZPSNV,PZS=ZZS                                         )  
   !   
   WRITE(*,'(A32,I4,A3,I4)') ' GRID BOXES CONTAINING SEA    : ',COUNT( ZSEA    (:) > 0. ),' / ',INI
   WRITE(*,'(A32,I4,A3,I4)') ' GRID BOXES CONTAINING WATER  : ',COUNT( ZWATER  (:) > 0. ),' / ',INI
   WRITE(*,'(A32,I4,A3,I4)') ' GRID BOXES CONTAINING NATURE : ',COUNT( ZNATURE (:) > 0. ),' / ',INI
   WRITE(*,'(A32,I4,A3,I4)') ' GRID BOXES CONTAINING TOWN   : ',COUNT( ZTOWN   (:) > 0. ),' / ',INI
   WRITE(*,*)'ZZ0    = ',ZZ0
   WRITE(*,*)'ZZ0EFF = ',ZZ0EFF
   WRITE(*,*)'ZZS = ',ZZS
   WRITE(*,*)'MINVAL(ZZS) = ',MINVAL(ZZS),' MAXVAL(ZZS) = ',MAXVAL(ZZS)
   !
ENDIF   
!
! --------------------------------------------------------------------------------------
IF (CFORCING_FILETYPE=='ASCII ' .OR. CFORCING_FILETYPE=='BINARY') CALL OPEN_CLOSE_BIN_ASC_FORC('CLOSE',CFORCING_FILETYPE,0,'R')

IF (CFORCING_FILETYPE=='NETCDF') CALL CLOSE_FILEIN_OL
IF (CTIMESERIES_FILETYPE=='OFFLIN') CALL CLOSE_FILEOUT_OL
!
!
!*    5.     Close parallelized I/O
!            ----------------------
!
CLOSE(ILUOUT)
!
WRITE(*,*) ' '
WRITE(*,*) '    --------------------------'
WRITE(*,*) '    | OFFLINE ENDS CORRECTLY |'
WRITE(*,*) '    --------------------------'
WRITE(*,*) ' '
!
CALL DEALLOC_SURFEX
IF (LHOOK) CALL DR_HOOK('OFFLINE',1,ZHOOK_HANDLE)
! --------------------------------------------------------------------------------------
!
END PROGRAM OFFLINE
