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
USE MODD_BEM_n, ONLY : B => BEM
USE MODD_BEM_OPTION_n, ONLY : BOP => BEM_OPTIONS
USE MODD_CH_EMIS_FIELD_n, ONLY : CHE => CH_EMIS_FIELD
USE MODD_CH_SNAP_n, ONLY : CHN => CH_EMIS_SNAP
USE MODD_DIAG_CUMUL_TEB_n, ONLY : DGCT => DIAG_CUMUL_TEB
USE MODD_DIAG_MISC_TEB_n, ONLY : DGMT => DIAG_MISC_TEB
USE MODD_DST_n, ONLY : DST => DST
USE MODD_GR_BIOG_n, ONLY : GB => GR_BIOG
USE MODD_SV_n, ONLY : SV => SV
USE MODD_TEB_GARDEN_PGD_EVOL_n, ONLY : TGDPE => TEB_GARDEN_PGD_EVOL
USE MODD_TEB_GARDEN_PGD_n, ONLY : TGDP => TEB_GARDEN_PGD
USE MODD_TEB_PANEL_n, ONLY : TPN => TEB_PANEL
!
USE MODD_CH_SURF_n, ONLY : CHU => CH_SURF
USE MODD_DIAG_IDEAL_n, ONLY : DGL => DIAG_IDEAL
USE MODD_DIAG_UTCI_TEB_n, ONLY : DGUT => DIAG_UTCI_TEB
USE MODD_FLAKE_n, ONLY : F => FLAKE
USE MODD_OCEAN_n, ONLY : O => OCEAN
USE MODD_SURF_ATM_SSO_n, ONLY : USS => SURF_ATM_SSO
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TGRO => TEB_GREENROOF_OPTIONS
USE MODD_TEB_n, ONLY : T => TEB
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS
USE MODD_TEB_VEG_n, ONLY : TVG => TEB_VEG_OPTIONS
USE MODD_WATFLUX_n, ONLY : W => WATFLUX
!
USE MODD_SEAFLUX_n, ONLY : S => SEAFLUX
!
USE MODD_IO_BUFF_n, ONLY : IOB => IO_BUFF
!
USE MODD_CH_FLAKE_n, ONLY : CHF => CH_FLAKE
USE MODD_CH_ISBA_n, ONLY : CHI => CH_ISBA
USE MODD_CH_SEAFLUX_n, ONLY : CHS => CH_SEAFLUX
USE MODD_CH_TEB_n, ONLY : CHT => CH_TEB
USE MODD_CH_WATFLUX_n, ONLY : CHW => CH_WATFLUX
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DGEI => DIAG_EVAP_ISBA
USE MODD_DIAG_FLAKE_n, ONLY : DGF => DIAG_FLAKE
USE MODD_DIAG_ISBA_n, ONLY : DGI => DIAG_ISBA
USE MODD_DIAG_MISC_FLAKE_n, ONLY : DGMF => DIAG_MISC_FLAKE
USE MODD_DIAG_MISC_ISBA_n, ONLY : DGMI => DIAG_MISC_ISBA
USE MODD_DIAG_MISC_TEB_OPTION_n, ONLY : DGMTO => DIAG_MISC_TEB_OPTIONS
USE MODD_DIAG_OCEAN_n, ONLY : DGO => DIAG_OCEAN
USE MODD_DIAG_SEAFLUX_n, ONLY : DGS => DIAG_SEAFLUX
USE MODD_DIAG_SEAICE_n, ONLY : DGSI => DIAG_SEAICE
USE MODD_DIAG_SURF_ATM_n, ONLY : DGU => DIAG_SURF_ATM
USE MODD_DIAG_TEB_n, ONLY : DGT => DIAG_TEB
USE MODD_DIAG_WATFLUX_n, ONLY : DGW => DIAG_WATFLUX
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
USE MODD_SURF_ATM_GRID_n, ONLY : UG => SURF_ATM_GRID
!
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
USE MODD_ISBA_n, ONLY : I => ISBA
!
USE MODD_ASSIM, ONLY : LASSIM, LAROME, LALADSURF, CASSIM_ISBA, NVAR, XF, XF_PATCH, &
                       NOBSTYPE, XAT2M_ISBA, XAHU2M_ISBA, CVAR, COBS, NECHGU, XI,  &
                       NBOUTPUT, XLAI_PASS, XBIO_PASS, CBIO, NIVAR, NIPERT, XYO,   &
                       NOBS, LOBSFILE, NPRINTLEV
!
USE MODD_FORC_ATM,       ONLY : CSV, XDIR_ALB, XSCA_ALB, XEMIS, XTSRAD, XTSUN, XZS, &
                                XZREF, XUREF, XTA, XQA, XSV, XU, XV, XSW_BANDS,     &
                                XZENITH, XAZIM, XCO2, XRHOA, XTSURF
!
#ifdef SFX_ARO
USE MODD_IO_SURF_ARO,ONLY : NGPTOT, NGPTOT_CAP, NPROMA, NINDX1, NINDX2, NBLOCK, NKPROMA
#endif
!
#ifdef SFX_NC
USE MODD_IO_SURF_NC,   ONLY : CFILEIN_NC, CFILEIN_NC_SAVE, CFILEPGD_NC, CFILEOUT_NC, LDEF, &
                              CLUOUT_NC
#endif
#ifdef SFX_ASC
USE MODD_IO_SURF_ASC,  ONLY : CFILEIN, CFILEIN_SAVE, CFILEPGD, CFILEOUT
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
 CHARACTER(LEN=1)             :: YVAR
!
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
INTEGER :: JL,JI,JJ,INBPERT
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
! Flag diag :
!
INTEGER :: I2M, IBEQ, IDSTEQ
LOGICAL :: GFRAC, GDIAG_GRID, GSURF_BUDGET, GRAD_BUDGET, GCOEF,    &
           GSURF_VARS, GDIAG_OCEAN, GDIAG_SEAICE, GWATER_PROFILE, &
           GSURF_EVAP_BUDGET, GFLOOD,  GPGD_ISBA, GCH_NO_FLUX_ISBA,&
           GSURF_MISC_BUDGET_ISBA, GPGD_TEB, GSURF_MISC_BUDGET_TEB
!
! ******************************************************************************************
!
#ifdef SFX_MPI
 CALL MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,ILEVEL,INFOMPI)
#endif
!
IF (LHOOK) CALL DR_HOOK('SODA',0,ZHOOK_HANDLE)
!
!
#ifdef SFX_MPI
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
#ifdef SFX_LFI
CLUOUT_LFI =  ADJUSTL(ADJUSTR(YLUOUT)//'.txt')
#endif
#ifdef SFX_NC
CLUOUT_NC = ADJUSTL(ADJUSTR(YLUOUT)//'.txt')
#endif
 CALL GET_LUOUT(CSURF_FILETYPE,ILUOUT)
OPEN(UNIT=ILUOUT,FILE=ADJUSTL(ADJUSTR(YLUOUT)//'.txt'),FORM='FORMATTED',ACTION='WRITE')

! Allocation of Surfex Types
CALL ALLOC_SURFEX(1)

! Reading all namelist (also assimilation)
CALL READ_ALL_NAMELISTS(CHF, CHI, CHS, CHU, CHT, CHW, &
                        DGEI, DGF, DGL, DGI, DGMF, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                        F, I, O, S, USS, TGRO, T, TOP, TVG, W, &
                        CSURF_FILETYPE,'ALL',.FALSE.)

! Go to SURFEX
CALL GOTO_SURFEX(1,.TRUE.)
!
! Setting input files read from namelist
IF ( CSURF_FILETYPE == "LFI   " ) THEN
#ifdef SFX_LFI
  CFILEIN_LFI      = CPREPFILE
  CFILE_LFI        = CPREPFILE
  CFILEIN_LFI_SAVE = CPREPFILE
  CFILEPGD_LFI     = CPGDFILE
  CFILEOUT_LFI     = CSURFFILE
#endif
ELSEIF ( CSURF_FILETYPE == "FA    " ) THEN
#ifdef SFX_FA
  CFILEIN_FA      = ADJUSTL(ADJUSTR(CPREPFILE)//'.fa')
  CFILEIN_FA_SAVE = ADJUSTL(ADJUSTR(CPREPFILE)//'.fa')
  CFILEPGD_FA     = ADJUSTL(ADJUSTR(CPGDFILE)//'.fa')
  CFILEOUT_FA  = ADJUSTL(ADJUSTR(CSURFFILE)//'.fa')
#endif
ELSEIF ( CSURF_FILETYPE == "ASCII " ) THEN
#ifdef SFX_ASC
  CFILEIN      = ADJUSTL(ADJUSTR(CPREPFILE)//'.txt')
  CFILEIN_SAVE = ADJUSTL(ADJUSTR(CPREPFILE)//'.txt')
  CFILEPGD     = ADJUSTL(ADJUSTR(CPGDFILE)//'.txt')
  CFILEOUT = ADJUSTL(ADJUSTR(CSURFFILE)//'.txt')
#endif
ELSEIF ( CSURF_FILETYPE == "NC    " ) THEN
#ifdef SFX_NC
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
 CALL READ_SURF(IOB, &
                CSURF_FILETYPE,'DIM_FULL  ',INI,  IRESP)
 CALL READ_SURF(IOB, &
                CSURF_FILETYPE,'DTCUR     ',TTIME,  IRESP)
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
ALLOCATE(XTSURF(INI))
!
! Indicate that zenith and azimuth angles are not initialized
XZENITH = XUNDEF
XAZIM   = XUNDEF
XRHOA   = 1.
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
#ifdef SFX_LFI
        CFILEIN_LFI      = YFILEIN 
        CFILE_LFI        = YFILEIN
        CFILEIN_LFI_SAVE = YFILEIN
#endif
      ELSEIF ( CSURF_FILETYPE == "FA    " ) THEN
        YFILEIN = TRIM(YFILEIN)//'.fa'
#ifdef SFX_FA
        CFILEIN_FA      = YFILEIN
        CFILEIN_FA_SAVE = YFILEIN
#endif
      ELSEIF ( CSURF_FILETYPE == "ASCII " ) THEN
        YFILEIN = TRIM(YFILEIN)//'.txt'
#ifdef SFX_ASC
        CFILEIN      = YFILEIN
        CFILEIN_SAVE = YFILEIN
#endif
      ELSEIF ( CSURF_FILETYPE == "NC    " ) THEN
        YFILEIN = TRIM(YFILEIN)//'.nc'
#ifdef SFX_NC
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
    CALL IO_BUFF_CLEAN_n(IOB)
    CALL INIT_SURF_ATM_n(CSURF_FILETYPE,YINIT, LLAND_USE, INI, ISV, ISW,  &
                         CSV, XCO2, XRHOA, XZENITH, XAZIM, XSW_BANDS,     &
                         XDIR_ALB, XSCA_ALB, XEMIS, XTSRAD, XTSURF,       &
                         IYEAR, IMONTH, IDAY, ZTIME,                      &
                          YATMFILE, YATMFILETYPE, YTEST              )
    !
    IF ( CASSIM_ISBA=='EKF  ' ) THEN
      !
      IF ( ISTEP==1 .AND. NIPERT==INBPERT ) THEN
        ALLOCATE(XLAI_PASS(U%NSIZE_NATURE,I%NPATCH)) 
        ALLOCATE(XBIO_PASS(U%NSIZE_NATURE,I%NPATCH))     
        ALLOCATE(XI       (U%NSIZE_NATURE,I%NPATCH,NVAR    ))
        ALLOCATE(XF       (U%NSIZE_NATURE,I%NPATCH,NVAR+1,NVAR    ))
        ALLOCATE(XF_PATCH (U%NSIZE_NATURE,I%NPATCH,NVAR+1,NOBSTYPE*NBOUTPUT))
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
              XF_PATCH(:,:,NIPERT,IIOBS) = I%XWG(:,1,:)
            CASE("LAI")
              XF_PATCH(:,:,NIPERT,IIOBS) = I%XLAI(:,:)
            CASE DEFAULT
              CALL ABOR1_SFX("Mapping of "//COBS(IOBS)//" is not defined in SODA!")
          END SELECT
        ENDDO
        !
        ! Prognostic fields for assimilation (Control vector)
        DO JL = 1,NVAR
          SELECT CASE (TRIM(CVAR(JL)))
            CASE("TG1")
              XF(:,:,NIPERT,JL) = I%XTG(:,1,:)
            CASE("TG2")
              XF(:,:,NIPERT,JL) = I%XTG(:,2,:)
            CASE("WG1")
              XF(:,:,NIPERT,JL) = I%XWG(:,1,:)
            CASE("WG2")
              XF(:,:,NIPERT,JL) = I%XWG(:,2,:)
            CASE("LAI")
              XF(:,:,NIPERT,JL) = I%XLAI(:,:)
            CASE DEFAULT
              CALL ABOR1_SFX("Mapping of "//TRIM(CVAR(JL))//" is not defined in SODA_!")
          END SELECT
        ENDDO
        !
        IF ( NIPERT==1 ) THEN
          !
          IF ( I%NPATCH==1 .AND. TRIM(CBIO)/="LAI" ) THEN
            CALL ABOR1_SFX("Mapping of "//CBIO//" is not defined in EKF with NPATCH=1!")
          ENDIF
          SELECT CASE (TRIM(CBIO))
            CASE("BIOMA1","BIOMASS1")
              XBIO_PASS(:,:) = I%XBIOMASS(:,1,:)
            CASE("BIOMA2","BIOMASS2")
              XBIO_PASS(:,:) = I%XBIOMASS(:,2,:)
            CASE("RESPI1","RESP_BIOM1")
              XBIO_PASS(:,:) = I%XRESP_BIOMASS(:,1,:)
            CASE("RESPI2","RESP_BIOM2")
              XBIO_PASS(:,:) = I%XRESP_BIOMASS(:,2,:)
            CASE("LAI")
              XBIO_PASS(:,:) = I%XLAI(:,:)
            CASE DEFAULT
              CALL ABOR1_SFX("Mapping of "//CBIO//" is not defined in EKF!")
          END SELECT
          !
          XLAI_PASS(:,:) = I%XLAI(:,:)          
          !
        ENDIF
        !
      ELSE
        !
        DO JL = 1,NVAR
          SELECT CASE (TRIM(CVAR(JL)))
            CASE("TG1")
              XI(:,:,JL) = I%XTG(:,1,:)
            CASE("TG2")
              XI(:,:,JL) = I%XTG(:,2,:)
            CASE("WG1")
              XI(:,:,JL) = I%XWG(:,1,:)
            CASE("WG2")
              XI(:,:,JL) = I%XWG(:,2,:)
            CASE("LAI")
              XI(:,:,JL) = I%XLAI(:,:)
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
  IF ( CASSIM_ISBA=="EKF  " ) THEN
    !
    IF (ISTEP==1) ALLOCATE(XYO(U%NSIZE_NATURE,NOBSTYPE*NBOUTPUT))
    !
    IF ( LOBSFILE ) THEN
      !
      YMFILE = 'CANARI_NATURE_'
      CALL GET_FILE_NAME(IYEAR,IMONTH,IDAY,IHOUR,YMFILE)
      OPEN(UNIT=55,FILE=TRIM(YMFILE)//".DAT",FORM='FORMATTED',STATUS='OLD',IOSTAT=ISTAT) 
      IF ( ISTAT==0 ) THEN
        !   If it exists, read observations
        DO JI = 1,U%NSIZE_NATURE
          READ (55,*)  (XYO(JI,NOBS+JJ),JJ=1,NOBSTYPE)
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
  CALL ABOR1_SFX("SODA: No observations read for LAI in OBS file - stop")
ENDIF
!
 CALL GET_SIZE_FULL_n(U, &
                      CSURF_FILETYPE,INI,U%NSIZE_FULL)
!
WRITE(*,*) "READING input files..."
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
IF (CASSIM_ISBA=="OI   " .OR. .NOT.LOBSFILE) THEN
  !
  !  Read atmospheric forecast fields from FA files 
#ifdef SFX_FA
  CFILEIN_FA = 'FG_OI_MAIN'
  CDNOMC     = 'oimain'                  ! new frame name
#endif
  !  Open FA file (LAM version with extension zone)
  CALL INIT_IO_SURF_n(YPROGRAM2,'EXTZON','SURF  ','READ ')
  !
  !  Read model forecast quantities
  IF (LAROME) THEN  
   CALL READ_SURF(IOB, &
                YPROGRAM2,'SURFACCPLUIE',    ZCON_RAIN    ,IRESP)
     CALL READ_SURF(IOB, &
                YPROGRAM2,'SURFACCNEIGE',  ZSTRAT_SNOW  ,IRESP)
    CALL READ_SURF(IOB, &
                YPROGRAM2,'SURFACCGRAUPEL',  ZCON_SNOW   ,IRESP)
    ! So far graupel has not been used
    !ZCON_SNOW=ZCON_SNOW+ZCON_GRAUPEL
    ZCON_RAIN(:) = 0.0
  ELSE    
    CALL READ_SURF(IOB, &
                YPROGRAM2,'SURFPREC.EAU.CON',ZCON_RAIN    ,IRESP)
    CALL READ_SURF(IOB, &
                YPROGRAM2,'SURFPREC.EAU.GEC',ZSTRAT_RAIN  ,IRESP)
    CALL READ_SURF(IOB, &
                YPROGRAM2,'SURFPREC.NEI.CON',ZCON_SNOW    ,IRESP)
    CALL READ_SURF(IOB, &
                YPROGRAM2,'SURFPREC.NEI.GEC',ZSTRAT_SNOW  ,IRESP)
  ENDIF
  !
  CALL READ_SURF(IOB, &
                YPROGRAM2,'ATMONEBUL.BASSE ',ZCLOUDS,IRESP)
  CALL READ_SURF(IOB, &
                YPROGRAM2,'SURFIND.TERREMER',ZLSM   ,IRESP)
  CALL READ_SURF(IOB, &
                YPROGRAM2,'SURFFLU.LAT.MEVA',ZEVAP  ,IRESP) ! accumulated fluxes (not available in LFI)
  !
  IF (.NOT.LALADSURF) THEN    
    CALL READ_SURF(IOB, &
                YPROGRAM2,'SURFXEVAPOTRANSP',ZEVAPTR,IRESP) ! not in ALADIN SURFEX
  ELSE
    ZEVAPTR(:) = 0.0
  ENDIF
  !
  !  Close FA file
  CALL END_IO_SURF_n(YPROGRAM2)
  CALL IO_BUFF_CLEAN_n(IOB)
  WRITE(*,*)'READ FG_OI_MAIN OK'
  !
ENDIF
!
!  Define FA file name for CANARI analysis
#ifdef SFX_FA
CFILEIN_FA = 'CANARI'        ! input CANARI analysis
CDNOMC     = 'canari'                  ! new frame name
#endif
!  Open FA file 
 CALL INIT_IO_SURF_n(YPROGRAM2,'EXTZON','SURF  ','READ ')

!  Read CANARI analysis
 CALL READ_SURF(IOB, &
                YPROGRAM2,'CLSTEMPERATURE  ',ZT2M ,IRESP)
 CALL READ_SURF(IOB, &
                YPROGRAM2,'CLSHUMI.RELATIVE',ZHU2M,IRESP)
 CALL READ_SURF(IOB, &
                YPROGRAM2,'SURFTEMPERATURE ',ZTS  ,IRESP)
 CALL READ_SURF(IOB, &
                YPROGRAM2,'SURFRESERV.NEIGE',ZSWE ,IRESP)
 CALL READ_SURF(IOB, &
                YPROGRAM2,'CLSVENT.ZONAL   ',ZUCLS,IRESP)
 CALL READ_SURF(IOB, &
                YPROGRAM2,'CLSVENT.MERIDIEN',ZVCLS,IRESP)  

!  Close CANARI file
 CALL END_IO_SURF_n(YPROGRAM2)
 CALL IO_BUFF_CLEAN_n(IOB)
 WRITE(*,*) 'READ CANARI OK'

!  Define FA file name for surface climatology
#ifdef SFX_FA
CFILEIN_FA = 'clim_isba'               ! input climatology
CDNOMC     = 'climat'                  ! new frame name
#endif
!  Open FA file 
 CALL INIT_IO_SURF_n(YPROGRAM2,'EXTZON','SURF  ','READ ')

!  Read climatology file
 CALL READ_SURF(IOB, &
                YPROGRAM2,'SURFRESERV.NEIGE',ZSWEC,IRESP)
 CALL READ_SURF(IOB, &
                YPROGRAM2,'SURFTEMPERATURE' ,ZTSC ,IRESP)

!  Close climatology file
 CALL END_IO_SURF_n(YPROGRAM2)
 CALL IO_BUFF_CLEAN_n(IOB)
WRITE(*,*) 'READ CLIMATOLOGY OK'

 CALL ASSIM_SET_SST(IOB, S, U, &
                    U%NSIZE_FULL,ZLSM,ZSST,ZSIC,YTEST)

IF ( .NOT. LASSIM ) CALL ABOR1_SFX("YOU CAN'T RUN SODA WITHOUT SETTING LASSIM=.TRUE. IN THE ASSIM NAMELIST")

ALLOCATE(GD_MASKEXT(INI))
GD_MASKEXT(:) = .FALSE.
!
ALLOCATE(ZLON(INI))
ALLOCATE(ZLAT(INI))
ZLON(:) = UG%XLON(:)
ZLAT(:) = UG%XLAT(:)        
!
GLKEEPEXTZONE = .TRUE.
!
WRITE(*,*) 'PERFORMIMG OFFLINE SURFEX DATA ASSIMILATION...'
 CALL ASSIM_SURF_ATM_n(CSURF_FILETYPE,INI,   &
                      ZCON_RAIN,  ZSTRAT_RAIN, ZCON_SNOW, ZSTRAT_SNOW, &
                      ZCLOUDS,    ZLSM,        ZEVAPTR,   ZEVAP,       &
                      ZSWEC,      ZTSC,       &
                      ZTS,        ZT2M,        ZHU2M,     ZSWE,        &
                      ZSST,       ZSIC,  ZUCLS,     ZVCLS,       &
                      YTEST, GD_MASKEXT, ZLON, ZLAT, GLKEEPEXTZONE )

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
INW = 1
IF (CTIMESERIES_FILETYPE=="NC    ") INW = 2
#ifdef SFX_NC
LDEF = .TRUE.
#endif
DO JNW = 1,INW
  CALL IO_BUFF_CLEAN_n(IOB)
  CALL WRITE_SURF_ATM_n(CTIMESERIES_FILETYPE,'ALL',LLAND_USE)
  CALL WRITE_DIAG_SURF_ATM_n(B, BOP, CHE, CHF, CHI, CHS, CHN, CHU, CHT, CHW, &
                             DGCT, DGEI, DGF, DGI, DGMF, DGMI, DGMT, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, DST, &
                             F, GB, I, S, UG, U, SV, TGDPE, TGDP, T, TOP, TPN, TVG, W, &
                             CTIMESERIES_FILETYPE,'ALL')
#ifdef SFX_NC
  LDEF = .FALSE.
#endif
ENDDO
!
INW = 1
IF (CSURF_FILETYPE=="NC    ") INW = 2
#ifdef SFX_NC
LDEF = .TRUE.
#endif
DO JNW = 1,INW
  !
  CALL FLAG_UPDATE(DGI, DGU, &
                   .FALSE.,.TRUE.,.FALSE.,.FALSE.)
  !
  GFRAC                  = .FALSE.
  GDIAG_GRID             = .TRUE.
  I2M                    = 0
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
  CALL FLAG_DIAG_UPDATE(CHF, CHI, CHS, CHT, CHW, DGEI, DGF, DGI, DGMF, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGW, I, U, &
                        GFRAC, GDIAG_GRID, I2M, GSURF_BUDGET, GRAD_BUDGET, GCOEF,  &
                        GSURF_VARS, IBEQ, IDSTEQ, GDIAG_OCEAN, GDIAG_SEAICE,       &
                        GWATER_PROFILE,                                            &
                        GSURF_EVAP_BUDGET, GFLOOD,  GPGD_ISBA, GCH_NO_FLUX_ISBA,   &
                        GSURF_MISC_BUDGET_ISBA, GPGD_TEB, GSURF_MISC_BUDGET_TEB    )
  ! 
  ! Store results from assimilation
  CALL WRITE_SURF_ATM_n(CSURF_FILETYPE,'ALL',LLAND_USE)
  CALL WRITE_DIAG_SURF_ATM_n(B, BOP, CHE, CHF, CHI, CHS, CHN, CHU, CHT, CHW, &
                             DGCT, DGEI, DGF, DGI, DGMF, DGMI, DGMT, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, DST, &
                             F, GB, I, S, UG, U, SV, TGDPE, TGDP, T, TOP, TPN, TVG, W, &
                             CSURF_FILETYPE,'ALL')
#ifdef SFX_NC
  LDEF = .FALSE.
#endif
  CALL IO_BUFF_CLEAN_n(IOB)
  !
ENDDO  
!
IF (CSURF_FILETYPE=='FA    ') THEN
#ifdef SFX_FA
  CALL FAIRME(IRET,NUNIT_FA,'UNKNOWN')
#endif
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
