!     ######spl
        PROGRAM NCPOST
!
!!    MODIFICATIONS
!!    -------------
!!      B. Decharme : partition pgd/prep (grid attributes are only in the PGD file)
!!
!-------------------------------------------------------------------------------
!
USE MODD_BEM_OPTION_n, ONLY : BOP => BEM_OPTIONS
USE MODD_BLD_DESCRIPTION_n, ONLY : BDD => BLD_DESC
USE MODD_CH_EMIS_FIELD_n, ONLY : CHE => CH_EMIS_FIELD
USE MODD_CH_ISBA_n, ONLY : CHI => CH_ISBA
USE MODD_CH_SEAFLUX_n, ONLY : CHS => CH_SEAFLUX
USE MODD_CH_SNAP_n, ONLY : CHN => CH_EMIS_SNAP
USE MODD_CH_SURF_n, ONLY : CHU => CH_SURF
USE MODD_CH_TEB_n, ONLY : CHT => CH_TEB
USE MODD_CH_WATFLUX_n, ONLY : CHW => CH_WATFLUX
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
USE MODD_DATA_SEAFLUX_n, ONLY : DTS => DATA_SEAFLUX
USE MODD_DATA_TEB_n, ONLY : DTT => DATA_TEB
USE MODD_DATA_TSZ0_n, ONLY : DTZ => DATA_TSZ0
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DGEI => DIAG_EVAP_ISBA
USE MODD_DIAG_FLAKE_n, ONLY : DGF => DIAG_FLAKE
USE MODD_DIAG_ISBA_n, ONLY : DGI => DIAG_ISBA
USE MODD_DIAG_MISC_ISBA_n, ONLY : DGMI => DIAG_MISC_ISBA
USE MODD_DIAG_MISC_TEB_OPTION_n, ONLY : DGMTO => DIAG_MISC_TEB_OPTIONS
USE MODD_DIAG_OCEAN_n, ONLY : DGO => DIAG_OCEAN
USE MODD_DIAG_SEAFLUX_n, ONLY : DGS => DIAG_SEAFLUX
USE MODD_DIAG_SEAICE_n, ONLY : DGSI => DIAG_SEAICE
USE MODD_DIAG_SURF_ATM_n, ONLY : DGU => DIAG_SURF_ATM
USE MODD_DIAG_TEB_n, ONLY : DGT => DIAG_TEB
USE MODD_DIAG_UTCI_TEB_n, ONLY : DGUT => DIAG_UTCI_TEB
USE MODD_DIAG_WATFLUX_n, ONLY : DGW => DIAG_WATFLUX
USE MODD_FLAKE_n, ONLY : F => FLAKE
USE MODD_FLAKE_SBL_n, ONLY : FSB => FLAKE_SBL
USE MODD_GR_BIOG_n, ONLY : GB => GR_BIOG
USE MODD_ISBA_CANOPY_n, ONLY : ICP => ISBA_CANOPY
USE MODD_ISBA_n, ONLY : I => ISBA
USE MODD_OCEAN_n, ONLY : O => OCEAN
USE MODD_SEAFLUX_n, ONLY : S => SEAFLUX
USE MODD_SEAFLUX_SBL_n, ONLY : SSB => SEAFLUX_SBL
USE MODD_SURF_ATM_GRID_n, ONLY : UG => SURF_ATM_GRID
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
USE MODD_SV_n, ONLY : SV => SV
USE MODD_TEB_CANOPY_n, ONLY : TCP => TEB_CANOPY
USE MODD_TEB_GARDEN_n, ONLY : TGD => TEB_GARDEN
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TGDO => TEB_GARDEN_OPTIONS
USE MODD_TEB_GREENROOF_n, ONLY : TGR => TEB_GREENROOF
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TGRO => TEB_GREENROOF_OPTIONS
USE MODD_TEB_n, ONLY : T => TEB
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS
USE MODD_TEB_VEG_n, ONLY : TVG => TEB_VEG_OPTIONS
USE MODD_WATFLUX_n, ONLY : W => WATFLUX
USE MODD_WATFLUX_SBL_n, ONLY : WSB => WATFLUX_SBL
!
USE MODD_IO_BUFF_n, ONLY : IOB => IO_BUFF
!
        USE MODD_IO_SURF_ASC
        USE MODD_SURF_PAR
        USE MODI_OPEN_FILEIN_OL
        USE MODI_CLOSE_FILEIN_OL
        USE MODI_READ_SURF
        USE MODE_POS_SURF
!
        USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
        USE PARKIND1  ,ONLY : JPRB
!
        USE MODI_GOTO_SURFEX
!
        USE MODI_END_IO_SURF_n
        USE MODI_INIT_IO_SURF_n
        IMPLICIT NONE        

        REAL, ALLOCATABLE, DIMENSION(:)   ::   ZLOC
        REAL, ALLOCATABLE, DIMENSION(:)   ::   ZWRK
        REAL, ALLOCATABLE, DIMENSION(:)   ::   XLON
        REAL, ALLOCATABLE, DIMENSION(:)   ::   XLAT
        INTEGER, ALLOCATABLE, DIMENSION(:)::   IWRK2
        CHARACTER(LEN=50)                 ::   YCOMMENT
        CHARACTER(LEN=50)                 ::   NOM_ARTICLE
        CHARACTER(LEN=12)                 ::   HREC
        CHARACTER(LEN=1)                  ::   PATCHFLAG
        CHARACTER(LEN=2)                  ::   YPAS,YLVL
        CHARACTER(LEN=10)                 ::   CGRID_TYPE
        LOGICAL                           ::   GFOUND
        LOGICAL                           ::   LINITS ! true if PGD has been run
        LOGICAL                           ::   LSXNAM ! true if NCPOST.nam present
        LOGICAL                           ::   LCOORD ! true if LONLAT.dat present
        LOGICAL                           ::   LGEO=.TRUE.  !

        INTEGER    ::   IRET
        INTEGER    ::   INI
        INTEGER    ::   INJ
        INTEGER    ::   IF, IC, IP
        INTEGER    ::   IFIELD, IWFIELD
        INTEGER    ::   IPATCH, JPATCH
        INTEGER    ::   IBEG, IEND

        

        !plm
        !=====================================================================
        real, allocatable, dimension(:,:,:)   ::   zfield3d
        real, allocatable, dimension(:,:)     ::   zfield2d
        character (len=40) :: cfile
        character (len=56) :: comlink
        integer    ::   inb_forc
        integer    ::   ji
        REAL(KIND=JPRB) :: ZHOOK_HANDLE
        !=====================================================================

        IF (LHOOK) CALL DR_HOOK('NCPOST',0,ZHOOK_HANDLE)
        CALL GOTO_SURFEX(1,.TRUE.)

        !=====================================================================
        !*
        !** get domain size and read latitudes and longitudes 
        !*
        !=====================================================================

        INQUIRE(FILE='LONLAT.dat',EXIST=LCOORD)
        IF (.NOT.LCOORD) THEN

           INQUIRE(FILE='PGD.txt', EXIST=LINITS)
      
           IF (.NOT. LINITS) THEN
              WRITE(*,*)' Now grid attributes are only in the PGD file'
              WRITE(*,*)' NO INPUT FILE FOUND FOR NCPOST'
              WRITE(*,*)' YOU SHOULD AT LEAST RUN PGD! '
              STOP
           ELSE
              CFILEIN='PGD.txt'
           ENDIF

           CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                               'ASCII ','FULL  ','SURF  ','READ ')

           CALL READ_SURF(IOB, &
                          'ASCII ','DIM_FULL', INI, IRET)
           CALL READ_SURF(IOB, &
                          'ASCII ','GRID_TYPE', CGRID_TYPE, IRET)

        
           ALLOCATE(XLON(INI))
           ALLOCATE(XLAT(INI))

           IF (CGRID_TYPE=='GAUSS     ') THEN
              CALL POSNAM(NUNIT,'FULL  '//' '//'LONGAUSS',GFOUND,NLUOUT)
           ELSE
              CALL POSNAM(NUNIT,'FULL  '//' '//'XLON',GFOUND,NLUOUT)
           ENDIF

           READ(NUNIT,FMT=*)
           READ(NUNIT,FMT='(A50)') YCOMMENT
           READ(NUNIT,FMT=*,ERR=100) XLON(:)

           IF (CGRID_TYPE=='GAUSS     ') THEN
              CALL POSNAM(NUNIT,'FULL  '//' '//'LATGAUSS',GFOUND,NLUOUT)
           ELSE
              CALL POSNAM(NUNIT,'FULL  '//' '//'XLAT',GFOUND,NLUOUT)
           ENDIF

           READ(NUNIT,FMT=*)
           READ(NUNIT,FMT='(A50)') YCOMMENT
           READ(NUNIT,FMT=*,ERR=100) XLAT(:)

           OPEN(UNIT=30,FILE='LONLAT.dat',FORM='FORMATTED')
           DO IP=1,INI
              WRITE(30,*)XLON(IP),XLAT(IP)
           ENDDO

           CALL END_IO_SURF_n('ASCII ')
                   
        ENDIF

        !=====================================================================
        !*
        !** read fields from netcdf output file
        !*
        !=====================================================================

        INQUIRE(FILE='NCPOST.nam',EXIST=LSXNAM)
        IF (.NOT.LSXNAM) THEN
           WRITE(*,*)' > NCPOST.nam does not exist'
           STOP
        ENDIF
        OPEN(UNIT=46,FILE='NCPOST.nam',FORM='FORMATTED')
        READ(46,'(A1,1X,A6,1X,A16,1X,A40)')PATCHFLAG,CMASK,HREC,CFILE

        CALL OPEN_FILEIN_OL
        CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                               'OFFLIN','FULL  ','SURF  ','READ ')
        CALL READ_SURF(IOB, &
                       'OFFLIN','DIM_FULL', INI, IRET)
        ALLOCATE(XLON(INI))
        ALLOCATE(XLAT(INI))
        OPEN(UNIT=30,FILE='LONLAT.dat',FORM='FORMATTED')
        DO IP=1,INI
           READ(30,*)XLON(IP),XLAT(IP)
        ENDDO

        CALL READ_SURF(IOB, &
                       'OFFLIN','NB_TIMESTP', INB_FORC, IRET)
        CALL READ_SURF(IOB, &
                          'OFFLIN','PATCH_NUMBER', IPATCH, IRET)
        CALL system('rm SXPOST.nc')
        comlink='ln -s '//CFILE//' SXPOST.nc'
        CALL system(comlink)

        IF (CMASK == 'FORC') THEN
           allocate(zfield2d(inb_forc-1,ini))
           CALL READ_SURF(IOB, &
                          'OFFLIN',HREC,zfield2d(:,:), IRET)
           do ji=1,ini
              write(50,*)xlon(ji),xlat(ji),zfield2d(1,ji)
           enddo
        ELSEIF (CMASK == 'SIMU') THEN
           IF (PATCHFLAG == '+') THEN
              allocate(zfield3d(ini,ipatch,inb_forc-1))
              CALL READ_SURF(IOB, &
                             'OFFLIN',HREC,zfield3d(:,:,:), IRET)
              do ji=1,ini
                 write(50,*)xlon(ji),xlat(ji),zfield3d(ji,1,1)
              enddo
           ELSE IF (PATCHFLAG == '-') THEN
              allocate(zfield2d(ini,inb_forc-1))
              CALL READ_SURF(IOB, &
                             'OFFLIN',HREC,zfield2d(:,:), IRET)
              do ji=1,ini
                 write(50,*)xlon(ji),xlat(ji),zfield2d(ji,1)
              enddo
           ENDIF
        ELSE
           write(*,*)' > ',CMASK,'NOT ALLOWED (only FORC|SIMU)'
           write(*,*)' > Update NCPOST.nam'
           STOP
        ENDIF

        CALL CLOSE_FILEIN_OL

        STOP
 100    CONTINUE
        WRITE(NLUOUT,*) ' '
        WRITE(NLUOUT,*) ' ERROR WHEN READING ARTICLE',HREC
        WRITE(NLUOUT,*) ' '
        

        IF (LHOOK) CALL DR_HOOK('NCPOST',1,ZHOOK_HANDLE)
        CONTAINS

        SUBROUTINE ERR_STOP(HREC,CFILEIN,NLUOUT)
        CHARACTER(LEN=12)   ::   HREC
        CHARACTER(LEN=*)    ::   CFILEIN
        INTEGER             ::   NLUOUT
        REAL(KIND=JPRB) :: ZHOOK_HANDLE
        IF (LHOOK) CALL DR_HOOK('ERR_STOP',0,ZHOOK_HANDLE)
        WRITE(NLUOUT,*) ' '
        WRITE(NLUOUT,*) ' ARTICLE ',TRIM(HREC),' NOT FOUND IN FILE ', CFILEIN
        WRITE(NLUOUT,*) ' '
        WRITE(*,*) ' '
        WRITE(*,*) ' ARTICLE ',TRIM(HREC),' NOT FOUND IN FILE ', CFILEIN
        WRITE(*,*) ' '
        STOP
        IF (LHOOK) CALL DR_HOOK('ERR_STOP',1,ZHOOK_HANDLE)
        END SUBROUTINE ERR_STOP

        !=====================================================================


        END PROGRAM NCPOST
