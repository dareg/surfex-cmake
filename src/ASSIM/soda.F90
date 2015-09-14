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
!  05/2013 B. Decharme New coupling variables XTSURF (for AGCM)
!----------------------------------------------------------------------------
!
USE MODD_OFF_SURFEX_n
USE MODE_MODELN_SURFEX_HANDLER
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
USE MODD_SURF_CONF, ONLY : CPROGNAME, CSOFTWARE
USE MODD_SURF_PAR,  ONLY : XUNDEF,NUNDEF
!
USE MODD_ASSIM, ONLY : LASSIM, LAROME, LALADSURF, CASSIM_ISBA, NVAR, XF, XF_PATCH, &
                       NOBSTYPE, XAT2M_ISBA, XAHU2M_ISBA, CVAR, COBS, NECHGU, XI,  &
                       NBOUTPUT, XLAI_PASS, XBIO_PASS, CBIO, NIVAR, XYO, NIFIC, &
                       NOBS, LOBSFILE, NPRINTLEV, LREAD_ALL, NENS,LBIAS_CORRECTION
!
USE MODD_FORC_ATM,       ONLY : CSV, XDIR_ALB, XSCA_ALB, XEMIS, XTSRAD, XTSUN, XZS, &
                                XZREF, XUREF, XTA, XQA, XSV, XU, XV, XSW_BANDS,     &
                                XZENITH, XAZIM, XCO2, XRHOA, XTSURF
!
USE MODD_WRITE_BIN,  ONLY : NWRITE
!
#ifdef SFX_ARO
USE MODD_IO_SURF_ARO,ONLY : NGPTOT, NGPTOT_CAP, NPROMA, NINDX1, NINDX2, NBLOCK, NKPROMA
#endif
!
#ifdef SFX_OL
USE MODD_IO_SURF_OL, ONLY : XSTART, XCOUNT, XSTRIDE, XSTARTW, XCOUNTW, LTIME_WRITTEN,  &
                            LDEFINED_NATURE, LDEFINED_SEA, LDEFINED_WATER, &
                            LDEFINED_TOWN, LDEFINED_SURF_ATM, LPARTW
                            
#endif
!
#ifdef SFX_NC
USE MODD_IO_SURF_NC,   ONLY : CFILEIN_NC, CFILEIN_NC_SAVE, CFILEPGD_NC, CFILEOUT_NC, LDEF, &
                              CLUOUT_NC
#endif
#ifdef SFX_ASC
USE MODD_IO_SURF_ASC,  ONLY : CFILEIN, CFILEIN_SAVE, CFILEPGD, CFILEOUT, LCREATED
#endif
#ifdef SFX_FA
USE MODD_IO_SURF_FA,   ONLY : CFILEIN_FA, CFILEIN_FA_SAVE, CFILEPGD_FA, CDNOMC, CFILEOUT_FA, &
                              NUNIT_FA, IVERBFA, LFANOCOMPACT
#endif
#ifdef SFX_LFI
USE MODD_IO_SURF_LFI,    ONLY : CFILEIN_LFI, CFILEIN_LFI_SAVE, &
                                CFILEPGD_LFI, CFILE_LFI, CLUOUT_LFI, CFILEOUT_LFI 
#endif
!
USE MODN_IO_OFFLINE,     ONLY : NAM_IO_OFFLINE, CNAMELIST, CPGDFILE, CPREPFILE, CSURFFILE, &
                                CSURF_FILETYPE, CTIMESERIES_FILETYPE, LLAND_USE, YALG_MPI, &
                                LDIAG_FA_NOCOMPACT, LOUT_TIMENAME, XIO_FRAC, LRESTART_2M
!
USE MODE_POS_SURF,  ONLY : POSNAM
!
USE MODI_SET_SURFEX_FILEIN
USE MODI_INIT_INDEX_MPI
USE MODI_GATHER_AND_WRITE_MPI
USE MODI_READ_AND_SEND_MPI
USE MODI_ABOR1_SFX
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
USE MODI_TEST_NAM_VAR_SURF
USE MODI_INIT_OUTPUT_NC_n
!
USE MODE_EKF, ONLY : GET_FILE_NAME, SET_FILEIN
!
USE YOMHOOK,             ONLY : LHOOK,DR_HOOK
USE PARKIND1,            ONLY : JPRB
!
IMPLICIT NONE
!
#ifdef SFX_MPI
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
 CHARACTER(LEN=3)             :: YVAR
!
 CHARACTER(LEN=100) :: YNAME
 CHARACTER(LEN=10)  :: YRANK
 CHARACTER(LEN=3) :: YENS
!
REAL,ALLOCATABLE, DIMENSION(:,:) :: ZYO_NAT, ZYO
REAL,ALLOCATABLE, DIMENSION(:)   :: ZNATURE
REAL,ALLOCATABLE, DIMENSION(:)   :: ZLSM                ! Land-Sea mask
REAL,ALLOCATABLE, DIMENSION(:)   :: ZCON_RAIN           ! Amount of convective liquid precipitation
REAL,ALLOCATABLE, DIMENSION(:)   :: ZSTRAT_RAIN         ! Amount of stratiform liquid precipitation
REAL,ALLOCATABLE, DIMENSION(:)   :: ZCON_SNOW           ! Amount of convective solid precipitation
REAL,ALLOCATABLE, DIMENSION(:)   :: ZSTRAT_SNOW         ! Amount of stratiform solid precipitation
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
!
LOGICAL, ALLOCATABLE, DIMENSION(:) :: GD_MASKEXT
LOGICAL :: GLKEEPEXTZONE
LOGICAL :: GFOUND
!
TYPE (DATE_TIME)                 :: TTIME               ! Current date and time  
!
 CHARACTER(LEN=14)                :: YTAG
!
INTEGER, DIMENSION(11)  :: IDATEF
INTEGER :: IDIM_FULL
INTEGER :: INI          ! grid dimension
INTEGER :: ISV                 ! Number of scalar species
INTEGER :: ISW                 ! Number of radiative bands 
INTEGER :: IYEAR, IMONTH, IDAY, IHOUR
INTEGER :: IYEAR_OUT, IMONTH_OUT, IDAY_OUT
INTEGER :: JL,JI,JJ,INB,ICPT
INTEGER :: INW, JNW
INTEGER :: ISTEP
INTEGER :: IOBS, IIOBS
INTEGER :: IGPCOMP
INTEGER :: ILUOUT
INTEGER :: ILUNAM
INTEGER :: IRET, INBFA
INTEGER :: IRESP, ISTAT               ! Response value
INTEGER :: INFOMPI, ILEVEL
INTEGER :: ISIZE, IENS
!
! Flag diag :
!
INTEGER :: I2M, IBEQ, IDSTEQ
LOGICAL :: GFRAC, GDIAG_GRID, GSURF_BUDGET, GRAD_BUDGET, GCOEF,    &
           GSURF_VARS, GDIAG_OCEAN, GDIAG_SEAICE, GWATER_PROFILE, &
           GSURF_EVAP_BUDGET, GFLOOD,  GPGD_ISBA, GCH_NO_FLUX_ISBA,&
           GSURF_MISC_BUDGET_ISBA, GPGD_TEB, GSURF_MISC_BUDGET_TEB
!
REAL(KIND=JPRB)                  :: ZHOOK_HANDLE
! ******************************************************************************************
!
INFOMPI=1
!
#ifdef SFX_MPI
 CALL MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,ILEVEL,INFOMPI)
#endif
!
IF (LHOOK) CALL DR_HOOK('SODA',0,ZHOOK_HANDLE)
!
 CSOFTWARE = 'SODA'
!
#ifdef SFX_MPI
NCOMM = MPI_COMM_WORLD
 CALL MPI_COMM_SIZE(NCOMM,NPROC,INFOMPI)
 CALL MPI_COMM_RANK(NCOMM,NRANK,INFOMPI)
#endif
!
 CALL PREP_LOG_MPI
!
!--------------------------------------
!
IF (NRANK==NPIO) THEN
  WRITE(*,*)
  WRITE(*,*) '   ------------------------------------'
  WRITE(*,*) '   |               SODA               |'
  WRITE(*,*) '   | SURFEX OFFLINE DATA ASSIMILATION |'
  WRITE(*,*) '   ------------------------------------'
  WRITE(*,*)
ENDIF
!
WRITE(YRANK,FMT='(I10)') NRANK
YNAME=TRIM(YLUOUT)//ADJUSTL(YRANK)
!
! Open ascii outputfile for writing
#ifdef SFX_LFI
CLUOUT_LFI =  ADJUSTL(ADJUSTR(YLUOUT)//'.txt')
#endif
#ifdef SFX_NC
CLUOUT_NC = ADJUSTL(ADJUSTR(YLUOUT)//'.txt')
#endif
 CALL GET_LUOUT('ASCII ',ILUOUT)
OPEN(UNIT=ILUOUT,FILE=ADJUSTL(ADJUSTR(YNAME)//'.txt'),FORM='FORMATTED',ACTION='WRITE')
!
! Read offline specific things
 CALL OPEN_NAMELIST('ASCII ',ILUNAM,CNAMELIST)
 CALL POSNAM(ILUNAM,'NAM_IO_OFFLINE',GFOUND)
IF (GFOUND) READ (UNIT=ILUNAM,NML=NAM_IO_OFFLINE)
 CALL CLOSE_NAMELIST('ASCII ',ILUNAM)
!
IF (NPROC==1) THEN 
  XIO_FRAC=1.
ELSE
  XIO_FRAC = MAX(MIN(XIO_FRAC,1.),0.)
ENDIF
!
 CALL TEST_NAM_VAR_SURF(ILUOUT,'CSURF_FILETYPE',CSURF_FILETYPE,'ASCII ','LFI   ','FA    ','NC    ')
 CALL TEST_NAM_VAR_SURF(ILUOUT,'CTIMESERIES_FILETYPE',CTIMESERIES_FILETYPE,'NETCDF','TEXTE ','BINARY',&
                                                                            'ASCII ','LFI   ','FA    ',&
                                                                            'NONE  ','OFFLIN','NC    ')  
!
IF (CTIMESERIES_FILETYPE=='NETCDF') CTIMESERIES_FILETYPE='OFFLIN'
!
! Setting input files read from namelist
IF ( CSURF_FILETYPE == "LFI   " ) THEN
#ifdef SFX_LFI
  CFILEIN_LFI      = CPREPFILE
  CFILE_LFI        = CPREPFILE
  CFILEIN_LFI_SAVE = CPREPFILE
  CFILEPGD_LFI     = CPGDFILE
#endif
ELSEIF ( CSURF_FILETYPE == "FA    " ) THEN
#ifdef SFX_FA
  CFILEIN_FA      = ADJUSTL(ADJUSTR(CPREPFILE)//'.fa')
  CFILEIN_FA_SAVE = ADJUSTL(ADJUSTR(CPREPFILE)//'.fa')
  CFILEPGD_FA     = ADJUSTL(ADJUSTR(CPGDFILE)//'.fa')
#endif
ELSEIF ( CSURF_FILETYPE == "ASCII " ) THEN
#ifdef SFX_ASC
  CFILEIN      = ADJUSTL(ADJUSTR(CPREPFILE)//'.txt')
  CFILEIN_SAVE = ADJUSTL(ADJUSTR(CPREPFILE)//'.txt')
  CFILEPGD     = ADJUSTL(ADJUSTR(CPGDFILE)//'.txt')
#endif
ELSEIF ( CSURF_FILETYPE == "NC    " ) THEN
#ifdef SFX_NC
  CFILEIN_NC      = ADJUSTL(ADJUSTR(CPREPFILE)//'.nc')
  CFILEIN_NC_SAVE = ADJUSTL(ADJUSTR(CPREPFILE)//'.nc')
  CFILEPGD_NC     = ADJUSTL(ADJUSTR(CPGDFILE)//'.nc')
#endif
ELSE
  CALL ABOR1_SFX(TRIM(CSURF_FILETYPE)//" is not implemented!")
ENDIF
!
! Allocation of Surfex Types
 CALL SURFEX_ALLOC(1)
 YSURF_CUR => YSURF_LIST(1)
!
! Reading all namelist (also assimilation)
CALL READ_ALL_NAMELISTS(YSURF_CUR, &
                        CSURF_FILETYPE,'ALL',.FALSE.)
!
!*     0.2.    Goto model of Surfex Types 
!
ICURRENT_MODEL = 1
!
CPROGNAME = CSURF_FILETYPE
!
IF (NRANK==NPIO) THEN
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
ENDIF
!
LREAD_ALL = .TRUE.
!
CALL INIT_INDEX_MPI(YSURF_CUR, &
                    CSURF_FILETYPE, YALG_MPI, XIO_FRAC, .FALSE.)
!
!
! Initialize time information
IYEAR    = NUNDEF
IMONTH   = NUNDEF
IDAY     = NUNDEF
ZTIME    = XUNDEF
 CALL SET_SURFEX_FILEIN(CSURF_FILETYPE,'PREP')
CALL INIT_IO_SURF_n(YSURF_CUR%DTCO, YSURF_CUR%DGU, YSURF_CUR%IOB, YSURF_CUR%U, &
                        CSURF_FILETYPE,'FULL  ','SURF  ','READ ') 
 CALL READ_SURF(YSURF_CUR%IOB, &
                CSURF_FILETYPE,'DIM_FULL  ',IDIM_FULL,  IRESP)
 CALL READ_SURF(YSURF_CUR%IOB, &
                CSURF_FILETYPE,'DTCUR     ',TTIME,  IRESP)
 CALL END_IO_SURF_n(CSURF_FILETYPE)
!
 CALL GET_SIZE_FULL_n(YSURF_CUR%U,CSURF_FILETYPE,IDIM_FULL,INI)
!
IF (ALLOCATED(NMASK_FULL)) DEALLOCATE(NMASK_FULL)
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
ALLOCATE(XTSURF(INI))
!
! Indicate that zenith and azimuth angles are not initialized
XZENITH = XUNDEF
XAZIM   = XUNDEF
XCO2    = 0.
XRHOA   = 1.
!
IF ( CASSIM_ISBA == 'EKF  ' ) THEN
  ! Has to do initialization for all the perturbations + 
  ! control + the real run at last
  INB = NVAR + 2
  ISIZE = NVAR
ELSEIF ( CASSIM_ISBA == 'OI   ' ) THEN
  INB = 1
ELSEIF ( CASSIM_ISBA == 'ENKF ' ) THEN
  INB = NENS
  IF (LBIAS_CORRECTION) INB = INB + 1
  ISIZE = NENS
ENDIF
!
IF (NRANK==NPIO) WRITE(*,*) "INITIALIZING SURFEX..."
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
  LREAD_ALL = .FALSE.  
  !
  DO NIFIC = INB,1,-1
    !
    ! If we have more than one initialization to do
    ! For last initialization, we must re-do the first.
    ! Could be avoided by introducing knowlegde of LASSIM on this level
    IF ( CASSIM_ISBA == 'EKF  ' .OR. CASSIM_ISBA == 'ENKF ' ) THEN
      !
      IF (CASSIM_ISBA == 'EKF  ') THEN
        IF ( NIFIC<INB ) THEN
          YMFILE = "PREP_"
          CALL GET_FILE_NAME(IYEAR,IMONTH,IDAY,IHOUR,YMFILE)              
          WRITE(YVAR,'(I1.1)') NIFIC-1
          YFILEIN = TRIM(YMFILE)//"_EKF_PERT"//ADJUSTL(YVAR)
        ELSE
          YFILEIN = "PREP_INIT"
        ENDIF
      ELSEIF (CASSIM_ISBA == 'ENKF ') THEN
        YMFILE = "PREP_"
        CALL GET_FILE_NAME(IYEAR,IMONTH,IDAY,IHOUR,YMFILE)              
        WRITE(YVAR,'(I1.1)') NIFIC
        YFILEIN = TRIM(YMFILE)//"_EKF_ENS"//ADJUSTL(YVAR)
      ENDIF
      !
      CALL SET_FILEIN(YFILEIN)
      !
    ENDIF 
    !
    CALL DEALLOC_SURF_ATM_n(YSURF_CUR)
    !
    IF ( CASSIM_ISBA=='EKF  ' .AND. NIFIC==1 ) LREAD_ALL = .TRUE.
    !    
    ! Initialize the SURFEX interface
    CALL IO_BUFF_CLEAN_n(YSURF_CUR%IOB)
    CALL INIT_SURF_ATM_n(YSURF_CUR, &
                            CSURF_FILETYPE,YINIT, LLAND_USE, INI, ISV, ISW,  &
                         CSV, XCO2, XRHOA, XZENITH, XAZIM, XSW_BANDS,     &
                         XDIR_ALB, XSCA_ALB, XEMIS, XTSRAD, XTSURF,       &
                         IYEAR, IMONTH, IDAY, ZTIME,                      &
                          YATMFILE, YATMFILETYPE, YTEST              )
    !
    IF ( CASSIM_ISBA=='EKF  ' .OR. CASSIM_ISBA=='ENKF ' ) THEN
      !
      IF ( ISTEP==1 .AND. NIFIC==INB ) THEN
        ALLOCATE(XLAI_PASS(YSURF_CUR%U%NSIZE_NATURE,YSURF_CUR%IM%I%NPATCH)) 
        ALLOCATE(XBIO_PASS(YSURF_CUR%U%NSIZE_NATURE,YSURF_CUR%IM%I%NPATCH))     
        IF (CASSIM_ISBA=='EKF  ') ALLOCATE(XI       (YSURF_CUR%U%NSIZE_NATURE,YSURF_CUR%IM%I%NPATCH,ISIZE   ))
        ALLOCATE(XF       (YSURF_CUR%U%NSIZE_NATURE,YSURF_CUR%IM%I%NPATCH,ISIZE+1,NVAR    ))
        ALLOCATE(XF_PATCH (YSURF_CUR%U%NSIZE_NATURE,YSURF_CUR%IM%I%NPATCH,ISIZE+1,NOBSTYPE*NBOUTPUT))
      ENDIF
      !
      IF ( CASSIM_ISBA=='EKF  ' .AND. NIFIC<INB .OR. CASSIM_ISBA=='ENKF ') THEN
        !
        ! Set the global state values for this control value
        DO IOBS = 1,NOBSTYPE
          IIOBS = (ISTEP-1)*NOBSTYPE + IOBS
          DO JI=1,YSURF_CUR%IM%I%NPATCH
            SELECT CASE (TRIM(COBS(IOBS)))
              CASE("T2M")
                XF_PATCH(:,JI,NIFIC,IIOBS) = XAT2M_ISBA(:,1)
              CASE("HU2M")
                XF_PATCH(:,JI,NIFIC,IIOBS) = XAHU2M_ISBA(:,1)
              CASE("WG1")
                XF_PATCH(:,JI,NIFIC,IIOBS) = YSURF_CUR%IM%I%XWG(:,1,JI)
              CASE("LAI")
                XF_PATCH(:,JI,NIFIC,IIOBS) = YSURF_CUR%IM%I%XLAI(:,JI)
              CASE DEFAULT
                CALL ABOR1_SFX("Mapping of "//COBS(IOBS)//" is not defined in SODA!")
            END SELECT
          ENDDO
        ENDDO
        !
        ! Prognostic fields for assimilation (Control vector)
        DO JL = 1,NVAR
          SELECT CASE (TRIM(CVAR(JL)))
            CASE("TG1")
              XF(:,:,NIFIC,JL) = YSURF_CUR%IM%I%XTG(:,1,:)
            CASE("TG2")
              XF(:,:,NIFIC,JL) = YSURF_CUR%IM%I%XTG(:,2,:)
            CASE("WG1")
              XF(:,:,NIFIC,JL) = YSURF_CUR%IM%I%XWG(:,1,:)
            CASE("WG2")
              XF(:,:,NIFIC,JL) = YSURF_CUR%IM%I%XWG(:,2,:)
            CASE("WG3")
              XF(:,:,NIFIC,JL) = YSURF_CUR%IM%I%XWG(:,3,:)              
            CASE("LAI")
              XF(:,:,NIFIC,JL) = YSURF_CUR%IM%I%XLAI(:,:)
            CASE DEFAULT
              CALL ABOR1_SFX("Mapping of "//TRIM(CVAR(JL))//" is not defined in SODA!")
          END SELECT
        ENDDO
        !
        IF ( NIFIC==1 ) THEN
          !
          DO JL = 1,NVAR
            IF (TRIM(CVAR(JL))=="LAI") THEN
              IF ( YSURF_CUR%IM%I%NPATCH==1 .AND. TRIM(CBIO)/="LAI" ) THEN
                CALL ABOR1_SFX("Mapping of "//CBIO//" is not defined in EKF with NPATCH=1!")
              ENDIF
              SELECT CASE (TRIM(CBIO))
                CASE("BIOMA1","BIOMASS1")
                  XBIO_PASS(:,:) = YSURF_CUR%IM%I%XBIOMASS(:,1,:)
                CASE("BIOMA2","BIOMASS2")
                  XBIO_PASS(:,:) = YSURF_CUR%IM%I%XBIOMASS(:,2,:)
                CASE("RESPI1","RESP_BIOM1")
                  XBIO_PASS(:,:) = YSURF_CUR%IM%I%XRESP_BIOMASS(:,1,:)
                CASE("RESPI2","RESP_BIOM2")
                  XBIO_PASS(:,:) = YSURF_CUR%IM%I%XRESP_BIOMASS(:,2,:)
                CASE("LAI")
                  XBIO_PASS(:,:) = YSURF_CUR%IM%I%XLAI(:,:)
                CASE DEFAULT
                  CALL ABOR1_SFX("Mapping of "//CBIO//" is not defined in EKF!")
              END SELECT
              !
              XLAI_PASS(:,:) = YSURF_CUR%IM%I%XLAI(:,:)          
              !
            ENDIF
            !
          ENDDO
        ENDIF
        !
      ELSE
        !
        DO JL = 1,NVAR
          SELECT CASE (TRIM(CVAR(JL)))
            CASE("TG1")
              XI(:,:,JL) = YSURF_CUR%IM%I%XTG(:,1,:)
            CASE("TG2")
              XI(:,:,JL) = YSURF_CUR%IM%I%XTG(:,2,:)
            CASE("WG1")
              XI(:,:,JL) = YSURF_CUR%IM%I%XWG(:,1,:)
            CASE("WG2")
              XI(:,:,JL) = YSURF_CUR%IM%I%XWG(:,2,:)
            CASE("WG3")
              XI(:,:,JL) = YSURF_CUR%IM%I%XWG(:,3,:)               
            CASE("LAI")
              XI(:,:,JL) = YSURF_CUR%IM%I%XLAI(:,:)
            CASE DEFAULT
              CALL ABOR1_SFX("Mapping of "//TRIM(CVAR(JL))//" is not defined in SODA!")
          END SELECT
        ENDDO        
        !
      ENDIF
      !
    ENDIF
    !
  ENDDO
  !
  IF ( CASSIM_ISBA=="EKF  " .OR. CASSIM_ISBA=="ENKF " ) THEN
    !
    IF (ISTEP==1) ALLOCATE(XYO(YSURF_CUR%U%NSIZE_NATURE,NOBSTYPE*NBOUTPUT))
    !
    IF ( LOBSFILE ) THEN
      !
      IF (ISTEP==1 .AND. NRANK==NPIO) THEN
        ALLOCATE(ZYO_NAT(YSURF_CUR%U%NDIM_NATURE,NOBSTYPE))
        ALLOCATE(ZYO(YSURF_CUR%U%NDIM_FULL,NOBSTYPE))
        ALLOCATE(ZNATURE(YSURF_CUR%U%NDIM_FULL))
      ENDIF
      !      
      YMFILE = 'CANARI_NATURE_'
      CALL GET_FILE_NAME(IYEAR,IMONTH,IDAY,IHOUR,YMFILE)
      ISTAT = 0
      IF (NRANK==NPIO) OPEN(UNIT=55,FILE=TRIM(YMFILE)//".DAT",FORM='FORMATTED',STATUS='OLD',IOSTAT=ISTAT) 
      IF ( ISTAT==0 ) THEN
        IF (NPROC>1) CALL GATHER_AND_WRITE_MPI(YSURF_CUR%U%XNATURE,ZNATURE)
        IF (NRANK==NPIO) THEN               
          !   If it exists, read observations
          DO JI = 1,YSURF_CUR%U%NSIZE_NATURE
            READ (55,*)  (ZYO_NAT(JI,NOBS+JJ),JJ=1,NOBSTYPE)
          ENDDO
          IF (NPROC>1) THEN       
            ZYO(:,:) = XUNDEF
            ICPT = 0
            DO JI = 1,YSURF_CUR%U%NDIM_FULL
              IF (ZNATURE(JI)>0.) THEN
                ICPT = ICPT + 1
                ZYO(JI,:) = ZYO_NAT(ICPT,:)
              ENDIF
            ENDDO
          ENDIF
        ENDIF
        IF (NPROC>1) THEN
          CALL READ_AND_SEND_MPI(ZYO,XYO,YSURF_CUR%U%NR_NATURE)
        ELSE
          XYO(:,:) = ZYO_NAT(:,:)
        ENDIF          
        NOBS = NOBS + NOBSTYPE      
        IF ( NPRINTLEV > 0 ) WRITE(*,*) 'read in obs: ', XYO(1,:), NOBS
        CLOSE(55)
      ENDIF
      !
      IF (ISTEP==NBOUTPUT .AND. NRANK==NPIO) DEALLOCATE(ZYO_NAT,ZYO,ZNATURE)
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
  CALL ABOR1_SFX("SODA: No observations read for LAI in OBS file - stop")
ENDIF
!
!
IF (NRANK==NPIO) WRITE(*,*) "READING input files..."
!
! Normal reading of needed FA fields
ALLOCATE(ZLSM        (INI))
ALLOCATE(ZCON_RAIN   (INI))
ALLOCATE(ZSTRAT_RAIN (INI))
ALLOCATE(ZCON_SNOW   (INI))
ALLOCATE(ZSTRAT_SNOW (INI))
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
!
ZTS(:) = XUNDEF
!
IF (CASSIM_ISBA=="OI   " .OR. .NOT.LOBSFILE) THEN
  !
  !  Read atmospheric forecast fields from FA files 
#ifdef SFX_FA
  CFILEIN_FA = 'FG_OI_MAIN'
  CDNOMC     = 'oimain'                  ! new frame name
#endif
  !  Open FA file (LAM version with extension zone)
CALL INIT_IO_SURF_n(YSURF_CUR%DTCO, YSURF_CUR%DGU, YSURF_CUR%IOB, YSURF_CUR%U, &
                        YPROGRAM2,'EXTZON','SURF  ','READ ') 
  !
  !  Read model forecast quantities
  IF (LAROME) THEN  
    CALL READ_SURF(YSURF_CUR%IOB, &
                YPROGRAM2,'SURFACCPLUIE',    ZSTRAT_RAIN    ,IRESP)
    CALL READ_SURF(YSURF_CUR%IOB, &
                YPROGRAM2,'SURFACCNEIGE',  ZSTRAT_SNOW  ,IRESP)
    CALL READ_SURF(YSURF_CUR%IOB, &
                YPROGRAM2,'SURFACCGRAUPEL',  ZCON_SNOW   ,IRESP)
    ! So far graupel has not been used
    !ZCON_SNOW=ZCON_SNOW+ZCON_GRAUPEL
    ZCON_RAIN(:) = 0.0
  ELSE    
    CALL READ_SURF(YSURF_CUR%IOB, &
                YPROGRAM2,'SURFPREC.EAU.CON',ZCON_RAIN    ,IRESP)
    CALL READ_SURF(YSURF_CUR%IOB, &
                YPROGRAM2,'SURFPREC.EAU.GEC',ZSTRAT_RAIN  ,IRESP)
    CALL READ_SURF(YSURF_CUR%IOB, &
                YPROGRAM2,'SURFPREC.NEI.CON',ZCON_SNOW    ,IRESP)
    CALL READ_SURF(YSURF_CUR%IOB, &
                YPROGRAM2,'SURFPREC.NEI.GEC',ZSTRAT_SNOW  ,IRESP)
  ENDIF
  !
  CALL READ_SURF(YSURF_CUR%IOB, &
                YPROGRAM2,'ATMONEBUL.BASSE ',ZCLOUDS,IRESP)
  CALL READ_SURF(YSURF_CUR%IOB, &
                YPROGRAM2,'SURFIND.TERREMER',ZLSM   ,IRESP)
  CALL READ_SURF(YSURF_CUR%IOB, &
                YPROGRAM2,'SURFFLU.LAT.MEVA',ZEVAP  ,IRESP) ! accumulated fluxes (not available in LFI)
  !
  IF (.NOT.LALADSURF) THEN    
    CALL READ_SURF(YSURF_CUR%IOB, &
                YPROGRAM2,'SURFXEVAPOTRANSP',ZEVAPTR,IRESP) ! not in ALADIN SURFEX
  ELSE
    ZEVAPTR(:) = 0.0
  ENDIF
  !
  !  Close FA file
  CALL END_IO_SURF_n(YPROGRAM2)
  CALL IO_BUFF_CLEAN_n(YSURF_CUR%IOB)
  WRITE(*,*)'READ FG_OI_MAIN OK'
  !
  !
  !  Define FA file name for CANARI analysis
#ifdef SFX_FA
  CFILEIN_FA = 'CANARI'        ! input CANARI analysis
  CDNOMC     = 'canari'                  ! new frame name
#endif
!  Open FA file 
CALL INIT_IO_SURF_n(YSURF_CUR%DTCO, YSURF_CUR%DGU, YSURF_CUR%IOB, YSURF_CUR%U, &
                        YPROGRAM2,'EXTZON','SURF  ','READ ')
!
!  Read CANARI analysis
  CALL READ_SURF(YSURF_CUR%IOB, &
                YPROGRAM2,'CLSTEMPERATURE  ',ZT2M ,IRESP)
  CALL READ_SURF(YSURF_CUR%IOB, &
                YPROGRAM2,'CLSHUMI.RELATIVE',ZHU2M,IRESP)
  CALL READ_SURF(YSURF_CUR%IOB, &
                YPROGRAM2,'SURFTEMPERATURE ',ZTS  ,IRESP)
  CALL READ_SURF(YSURF_CUR%IOB, &
                YPROGRAM2,'SURFRESERV.NEIGE',ZSWE ,IRESP)
  CALL READ_SURF(YSURF_CUR%IOB, &
                 YPROGRAM2,'CLSVENT.ZONAL   ',ZUCLS,IRESP)
  CALL READ_SURF(YSURF_CUR%IOB, &
                YPROGRAM2,'CLSVENT.MERIDIEN',ZVCLS,IRESP)  

!  Close CANARI file
  CALL END_IO_SURF_n(YPROGRAM2)
  CALL IO_BUFF_CLEAN_n(YSURF_CUR%IOB)
  WRITE(*,*) 'READ CANARI OK'
  !
ENDIF
!
IF (CASSIM_ISBA=="OI   ") THEN
  !
  !  Define FA file name for surface climatology
#ifdef SFX_FA
  CFILEIN_FA = 'clim_isba'               ! input climatology
  CDNOMC     = 'climat'                  ! new frame name
#endif
!  Open FA file 
CALL INIT_IO_SURF_n(YSURF_CUR%DTCO, YSURF_CUR%DGU, YSURF_CUR%IOB, YSURF_CUR%U, &
                        YPROGRAM2,'EXTZON','SURF  ','READ ')
!
!  Read climatology file
  CALL READ_SURF(YSURF_CUR%IOB, &
                YPROGRAM2,'SURFRESERV.NEIGE',ZSWEC,IRESP)
  CALL READ_SURF(YSURF_CUR%IOB, &
                YPROGRAM2,'SURFTEMPERATURE' ,ZTSC ,IRESP)

!  Close climatology file
  CALL END_IO_SURF_n(YPROGRAM2)
  CALL IO_BUFF_CLEAN_n(YSURF_CUR%IOB)
  WRITE(*,*) 'READ CLIMATOLOGY OK'
  !
ENDIF
!
 CALL ASSIM_SET_SST(YSURF_CUR%DTCO, YSURF_CUR%DGU, YSURF_CUR%IOB, YSURF_CUR%SM%S, YSURF_CUR%U, &
                    INI,ZLSM,ZSST,ZSIC,YTEST)

IF ( .NOT. LASSIM ) CALL ABOR1_SFX("YOU CAN'T RUN SODA WITHOUT SETTING LASSIM=.TRUE. IN THE ASSIM NAMELIST")
!
ALLOCATE(GD_MASKEXT(INI))
GD_MASKEXT(:) = .FALSE.
!
ALLOCATE(ZLON(INI))
ALLOCATE(ZLAT(INI))
ZLON(:) = YSURF_CUR%UG%XLON(:)
ZLAT(:) = YSURF_CUR%UG%XLAT(:)        
!
GLKEEPEXTZONE = .TRUE.
!
IF (NRANK==NPIO) WRITE(*,*) 'PERFORMIMG OFFLINE SURFEX DATA ASSIMILATION...'
 CALL ASSIM_SURF_ATM_n(YSURF_CUR%IM%DGMI, YSURF_CUR%IM%IG, YSURF_CUR%IM%I, YSURF_CUR%SM%S, &
                       YSURF_CUR%U, YSURF_CUR%TM%T, YSURF_CUR%TM%TOP, YSURF_CUR%WM%W, &
                       CSURF_FILETYPE,INI,   &
                      ZCON_RAIN,  ZSTRAT_RAIN, ZCON_SNOW, ZSTRAT_SNOW, &
                      ZCLOUDS,    ZLSM,        ZEVAPTR,   ZEVAP,       &
                      ZSWEC,      ZTSC,       &
                      ZTS,        ZT2M,        ZHU2M,     ZSWE,        &
                      ZSST,       ZSIC,  ZUCLS,     ZVCLS,       &
                      YTEST, GD_MASKEXT, ZLON, ZLAT, GLKEEPEXTZONE )
!
DEALLOCATE(ZCON_RAIN)
DEALLOCATE(ZSTRAT_RAIN)
DEALLOCATE(ZCON_SNOW)
DEALLOCATE(ZSTRAT_SNOW)
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
IF (CTIMESERIES_FILETYPE=="OFFLIN") THEN
  CALL INIT_OUTPUT_OL_n (YSURF_CUR)
ENDIF
!
ZTIME_OUT  = ZTIME
IDAY_OUT   = IDAY
IMONTH_OUT = IMONTH
IYEAR_OUT  = IYEAR
!
IF (NRANK==NPIO) THEN
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
  CFILEOUT_NC = ADJUSTL(ADJUSTR(CSURFFILE)//'.'//YTAG//'.nc')  
  !
  IF (CSURF_FILETYPE=='FA    ') THEN
#ifdef SFX_FA
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
#endif
  END IF
  !
ENDIF
!
ISIZE = 1
IF (CASSIM_ISBA=="ENKF ") THEN
  ISIZE = NENS
  IF (LBIAS_CORRECTION) ISIZE = ISIZE + 1
ENDIF
!
NWRITE = 1
XSTARTW = 1
LTIME_WRITTEN(:)=.FALSE.
!
DO IENS = 1,ISIZE
  !
  IF (CASSIM_ISBA=="ENKF ") THEN
    !
    YMFILE = "PREP_"
    CALL GET_FILE_NAME(IYEAR,IMONTH,IDAY,IHOUR,YMFILE)              
    WRITE(YVAR,'(I3)') IENS
    YFILEIN = TRIM(YMFILE)//"_EKF_ENS"//ADJUSTL(YVAR)
    CALL SET_FILEIN(YFILEIN)
    !
    LREAD_ALL = .TRUE.
    !
    CALL DEALLOC_SURF_ATM_n(YSURF_CUR)    
    !
    ! Initialize the SURFEX interface
    CALL IO_BUFF_CLEAN_n(YSURF_CUR%IOB)
    CALL INIT_SURF_ATM_n(YSURF_CUR, &
                            CSURF_FILETYPE,YINIT, LLAND_USE, INI, ISV, ISW,  &
                         CSV, XCO2, XRHOA, XZENITH, XAZIM, XSW_BANDS,     &
                         XDIR_ALB, XSCA_ALB, XEMIS, XTSRAD, XTSURF,       &
                         IYEAR, IMONTH, IDAY, ZTIME,                      &
                          YATMFILE, YATMFILETYPE, YTEST              )
    !
    !
    DO JL=1,NVAR
      !
      ! Update the modified values
      SELECT CASE (TRIM(CVAR(JL)))
        CASE("TG1")
          YSURF_CUR%IM%I%XTG(:,1,:) = XF(:,:,IENS,JL)
        CASE("TG2")
          YSURF_CUR%IM%I%XTG(:,2,:) = XF(:,:,IENS,JL)
        CASE("WG1")
          YSURF_CUR%IM%I%XWG(:,1,:) = XF(:,:,IENS,JL)
        CASE("WG2")
          YSURF_CUR%IM%I%XWG(:,2,:) = XF(:,:,IENS,JL)
        CASE("WG3")
          YSURF_CUR%IM%I%XWG(:,3,:) = XF(:,:,IENS,JL)          
        CASE("LAI") 
          YSURF_CUR%IM%I%XLAI(:,:) = XF(:,:,IENS,JL)
        CASE DEFAULT
          CALL ABOR1_SFX("Mapping of "//TRIM(CVAR(JL))//" is not defined in EKF!")
      END SELECT
    ENDDO
    !
  ENDIF
  !
#ifdef SFX_NC  
  LDEF = .TRUE.
#endif
  !
  IF (CTIMESERIES_FILETYPE=="NC    ") THEN
    CALL INIT_OUTPUT_NC_n (YSURF_CUR%TM%BDD, YSURF_CUR%CHE, YSURF_CUR%CHN, YSURF_CUR%CHU, &
                               YSURF_CUR%SM%DTS, YSURF_CUR%TM%DTT, YSURF_CUR%DTZ, YSURF_CUR%IM%I, &
                               YSURF_CUR%UG, YSURF_CUR%U, YSURF_CUR%DGU)                   
  ENDIF
  !
  INW = 1
  IF (CTIMESERIES_FILETYPE=="NC    ") INW = 2
  !
  DO JNW = 1,INW
    CALL IO_BUFF_CLEAN_n(YSURF_CUR%IOB)
    CALL WRITE_SURF_ATM_n(YSURF_CUR, &
                        CTIMESERIES_FILETYPE,'ALL',LLAND_USE)
    CALL WRITE_DIAG_SURF_ATM_n(YSURF_CUR, &
                             CTIMESERIES_FILETYPE,'ALL')
#ifdef SFX_NC
      LDEF = .FALSE.
#endif
  ENDDO
  !
  CALL FLAG_UPDATE(YSURF_CUR%IM%DGI, YSURF_CUR%DGU,.FALSE.,.TRUE.,.FALSE.,.FALSE.)
  !
  IF (LRESTART_2M) THEN
    I2M       = 1
    GPGD_ISBA = .TRUE.
  ELSE
    I2M       = 0
    GPGD_ISBA = .FALSE.
  ENDIF  
  GFRAC                  = .TRUE.  
  GDIAG_GRID             = .TRUE.
  GSURF_BUDGET           = .FALSE.
  GRAD_BUDGET            = .FALSE.
  GCOEF                  = .FALSE.
  GSURF_VARS             = .FALSE.
  IBEQ                   = 0
  IDSTEQ                 = 0
  GDIAG_OCEAN            = .FALSE.
  GDIAG_SEAICE           = .FALSE.
  GWATER_PROFILE         = .FALSE.
  GSURF_EVAP_BUDGET      = .FALSE.
  GFLOOD                 = .FALSE.
  GPGD_ISBA              = .FALSE.  
  GCH_NO_FLUX_ISBA       = .FALSE.
  GSURF_MISC_BUDGET_ISBA = .FALSE.
  GPGD_TEB               = .FALSE.
  GSURF_MISC_BUDGET_TEB  = .FALSE.  
  !
    CALL FLAG_DIAG_UPDATE(YSURF_CUR%FM%CHF, YSURF_CUR%IM%CHI, YSURF_CUR%SM%CHS, YSURF_CUR%TM%CHT, &
                          YSURF_CUR%WM%CHW, YSURF_CUR%IM%DGEI, YSURF_CUR%FM%DGF, YSURF_CUR%IM%DGI, &
                          YSURF_CUR%FM%DGMF, YSURF_CUR%IM%DGMI, YSURF_CUR%TM%DGMTO, YSURF_CUR%SM%DGO, &
                          YSURF_CUR%SM%DGS, YSURF_CUR%SM%DGSI, YSURF_CUR%DGU, YSURF_CUR%TM%DGT, &
                          YSURF_CUR%WM%DGW, YSURF_CUR%IM%I, YSURF_CUR%U, &          
                        GFRAC, GDIAG_GRID, I2M, GSURF_BUDGET, GRAD_BUDGET, GCOEF,  &
                        GSURF_VARS, IBEQ, IDSTEQ, GDIAG_OCEAN, GDIAG_SEAICE,       &
                        GWATER_PROFILE,                                            &
                        GSURF_EVAP_BUDGET, GFLOOD,  GPGD_ISBA, GCH_NO_FLUX_ISBA,   &
                        GSURF_MISC_BUDGET_ISBA, GPGD_TEB, GSURF_MISC_BUDGET_TEB    )
  ! 
  YENS = '   '
  IF (ISIZE>1) WRITE(YENS,'(I3)') IENS
  !
  IF ( CSURF_FILETYPE == "LFI   " ) THEN
#ifdef SFX_LFI
    CFILEOUT_LFI     = TRIM(TRIM(CSURFFILE)//ADJUSTL(YENS))
#endif
  ELSEIF ( CSURF_FILETYPE == "FA    " ) THEN
#ifdef SFX_FA
    CFILEOUT_FA  = ADJUSTL(TRIM(ADJUSTR(CSURFFILE)//ADJUSTL(YENS))//'.fa')
#endif
  ELSEIF ( CSURF_FILETYPE == "ASCII " ) THEN
#ifdef SFX_ASC
    CFILEOUT = ADJUSTL(TRIM(ADJUSTR(CSURFFILE)//ADJUSTL(YENS))//'.txt')
    LCREATED = .FALSE.
#endif
  ELSEIF ( CSURF_FILETYPE == "NC    " ) THEN
#ifdef SFX_NC
    CFILEOUT_NC = ADJUSTL(TRIM(ADJUSTR(CSURFFILE)//ADJUSTL(YENS))//'.nc')
#endif
  ELSE
    CALL ABOR1_SFX(TRIM(CSURF_FILETYPE)//" is not implemented!")
  ENDIF
  !
#ifdef SFX_NC
  LDEF = .TRUE.
#endif
  !  
  IF (CSURF_FILETYPE=="NC    ") THEN
    CALL INIT_OUTPUT_NC_n (YSURF_CUR%TM%BDD, YSURF_CUR%CHE, YSURF_CUR%CHN, YSURF_CUR%CHU, &
                         YSURF_CUR%SM%DTS, YSURF_CUR%TM%DTT, YSURF_CUR%DTZ, YSURF_CUR%IM%I, &
                         YSURF_CUR%UG, YSURF_CUR%U, YSURF_CUR%DGU)
  ENDIF
  !  
  INW = 1
  IF (CSURF_FILETYPE=="NC    ") INW = 2
  !  
  DO JNW = 1,INW
    !
    CALL IO_BUFF_CLEAN_n(YSURF_CUR%IOB)
    !  
    ! Store results from assimilation
    CALL WRITE_SURF_ATM_n(YSURF_CUR, &
                        CSURF_FILETYPE,'ALL',LLAND_USE)
    IF (YSURF_CUR%DGU%LREAD_BUDGETC.AND..NOT.YSURF_CUR%IM%DGEI%LRESET_BUDGETC) THEN
      CALL WRITE_DIAG_SURF_ATM_n(YSURF_CUR, &
                             CSURF_FILETYPE,'ALL')
    ENDIF
#ifdef SFX_NC
    LDEF = .FALSE.
#endif
    !
  ENDDO  
  !
ENDDO
!
IF (NRANK==NPIO .AND. CSURF_FILETYPE=='FA    ') THEN
#ifdef SFX_FA
  CALL FAIRME(IRET,NUNIT_FA,'UNKNOWN')
#endif
END IF
!
!*    3.     Close parallelized I/O
!            ----------------------
!
IF (NRANK==NPIO) THEN
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
ENDIF
!
CLOSE(ILUOUT)
!
 CALL SURFEX_DEALLO
!
IF (ALLOCATED(NINDEX)) DEALLOCATE(NINDEX)
IF (ALLOCATED(NSIZE_TASK)) DEALLOCATE(NSIZE_TASK)
!
IF (ASSOCIATED(NWORK)) DEALLOCATE(NWORK)
IF (ASSOCIATED(XWORK)) DEALLOCATE(XWORK)
IF (ASSOCIATED(NWORK2)) DEALLOCATE(NWORK2)
IF (ASSOCIATED(XWORK2)) DEALLOCATE(XWORK2)
IF (ASSOCIATED(XWORK3)) DEALLOCATE(XWORK3)
IF (ASSOCIATED(NWORK_FULL)) DEALLOCATE(NWORK_FULL)
IF (ASSOCIATED(XWORK_FULL)) DEALLOCATE(XWORK_FULL)
IF (ASSOCIATED(NWORK2_FULL)) DEALLOCATE(NWORK2_FULL)
IF (ASSOCIATED(XWORK2_FULL)) DEALLOCATE(XWORK2_FULL)
!
 CALL END_LOG_MPI
!
IF (LHOOK) CALL DR_HOOK('SODA',1,ZHOOK_HANDLE)
!
#ifdef SFX_MPI
 CALL MPI_FINALIZE(INFOMPI)
#endif
!
!-------------------------------------------------------------------------------
!
END PROGRAM SODA
