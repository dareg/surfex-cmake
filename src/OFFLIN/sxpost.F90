!     ######spl
        PROGRAM SXPOST
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
        USE MODI_READ_SURF
        USE MODI_GET_LUOUT
        USE MODE_POS_SURF
        USE MODD_IO_SURF_OL, ONLY : XSTART,XCOUNT,XSTRIDE,LPARTR
        
!
        USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
        USE PARKIND1  ,ONLY : JPRB
!
        USE MODI_ABOR1_SFX
!
        USE MODI_GOTO_SURFEX
        USE MODI_GET_LUOUT
!
        USE MODI_ALLOC_SURFEX
        USE MODI_DEALLOC_SURFEX
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
        CHARACTER(LEN=6)                  ::   CMASK_SAVE
        LOGICAL                           ::   GFOUND
        LOGICAL                           ::   LINITS ! true if PGD has been run
        LOGICAL                           ::   LINITP ! true if PREP has been run
        LOGICAL                           ::   LSXNAM ! true if SXPOST.nam present
        LOGICAL                           ::   LCOORD ! true if LONLAT.dat present
        LOGICAL                           ::   LGEO=.TRUE.  !

        CHARACTER(LEN=28)  :: YLUOUT ='LISTING_SXPOST                '  ! name of listing

        INTEGER    ::   IRET
        INTEGER    ::   INI, INI_N
        INTEGER    ::   INJ
        INTEGER    ::   IF, IC, IP
        INTEGER    ::   IFIELD, IWFIELD
        INTEGER    ::   IPATCH, JPATCH
        INTEGER    ::   IBEG, IEND
        REAL(KIND=JPRB) :: ZHOOK_HANDLE

        
        !=====================================================================
        !*
        !** check if file exists
        !*
        !=====================================================================
        IF (LHOOK) CALL DR_HOOK('SXPOST',0,ZHOOK_HANDLE)
        CALL ALLOC_SURFEX(1)

        CALL GET_LUOUT('ASCII ',NLUOUT)
        OPEN(UNIT=NLUOUT,FILE=ADJUSTL(ADJUSTR(YLUOUT)//'.txt'),&
                        FORM='FORMATTED',ACTION='WRITE')

        INQUIRE(FILE='SXPOST.nam', EXIST=LSXNAM)
        IF (.NOT.LSXNAM) THEN
           WRITE(*,*)' > SXPOST.nam missing'
           CALL ABOR1_SFX('SXPOST: NAMELIST SXPOST.nam MISSING')
        ENDIF

        INQUIRE(FILE='PGD.txt',    EXIST=LINITS)
        INQUIRE(FILE='PREP.txt',   EXIST=LINITP)

        IF (.NOT. LINITP .AND. .NOT. LINITS) THEN
           WRITE(*,*)' NO INPUT FILE FOUND FOR SXPOST'
           WRITE(*,*)' YOU SHOULD AT LEAST RUN PGD! '
           CALL ABOR1_SFX('SXPOST: NO INPUT FILE')
        ENDIF

        CFILEIN = 'PGD.txt'

        !=====================================================================
        !*
        !** get number of patches
        !*
        !=====================================================================

        CALL GOTO_SURFEX(1,.TRUE.)

        CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                            'ASCII ','FULL  ','SURF  ','READ ')

        CALL READ_SURF(IOB, &
                       'ASCII ','DIM_FULL', INI, IRET)
        CALL READ_SURF(IOB, &
                       'ASCII ','GRID_TYPE', CGRID_TYPE, IRET)
        CALL READ_SURF(IOB, &
                       'ASCII ','DIM_NATURE', INI_N, IRET)        
           
        CALL END_IO_SURF_n('ASCII ')

        NFULL = INI

        IF (INI_N.NE.0) THEN
          CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                            'ASCII ','NATURE','SURF  ','READ ')

          CALL READ_SURF(IOB, &
                       'ASCII ','PATCH_NUMBER', IPATCH, IRET)

          CALL END_IO_SURF_n('ASCII ')
        ENDIF

        !=====================================================================
        !*
        !** get domain size and read latitudes and longitudes 
        !*
        !=====================================================================
        CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                            'ASCII ','FULL  ','SURF  ','READ ')

        OPEN(UNIT=45,FILE='SXPOST.nam',FORM='FORMATTED')
        READ(45,*)IFIELD

        INQUIRE(FILE='LONLAT.dat',EXIST=LCOORD)
        ALLOCATE(XLON(INI))
        ALLOCATE(XLAT(INI))
        OPEN(UNIT=30,FILE='LONLAT.dat',FORM='FORMATTED')

        IF (LCOORD) THEN
          DO IP=1,INI
            READ(30,*)XLON(IP),XLAT(IP)
          ENDDO
        ELSE
          WRITE(*,*) 'LONLAT.DAT file missing !'
          IF (CGRID_TYPE=='GAUSS     ') THEN
            CALL POSNAM(NUNIT,'FULL  '//' '//'LONGAUSS',GFOUND,NLUOUT)
          ELSE
            CALL POSNAM(NUNIT,'FULL  '//' '//'XLON',GFOUND,NLUOUT)
          ENDIF

         IF (.NOT.GFOUND) THEN
           CALL ERR_STOP('XLON            ',CFILEIN,NLUOUT)
         ELSE
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

           DO IP=1,INI
             WRITE(30,*)XLON(IP),XLAT(IP)
           ENDDO
         ENDIF
       ENDIF

       CALL END_IO_SURF_n('ASCII ')

       IF (IFIELD==0) STOP

       !=====================================================================
       !*
       !** read 2d fields from PGD.txt or PREP.txt if exists
       !*
       !=====================================================================

       ! Search var first in PREP file
       IC=0

       DO IF=1,IFIELD

         READ(45,'(A1,1X,A6,1X,A16)') PATCHFLAG,CMASK,HREC
         CMASK_SAVE = CMASK

         IF (PATCHFLAG == '+') THEN
           INJ = INI * IPATCH
         ELSE IF (PATCHFLAG == '-') THEN
           INJ = INI
         ELSE
           PRINT*,' '
           PRINT*,' WRONG PATCHFLAG IN SXPOST.nam '
           PRINT*,' USE + FOR PATCHED VARIABLES '
           PRINT*,' USE - FOR UNPATCHED VARIABLES '
           PRINT*,' '
           PRINT*,' SYNTAX OF SXPOST.nam SHOULD LOOK: '
           PRINT*,' '
           PRINT*,'2 '
           PRINT*,'- FULL   ZS  '
           PRINT*,'+ NATURE TG1 '
           PRINT*,' '
           CALL ABOR1_SFX('SXPOST: WRONG PATCHFLAG')
         ENDIF

         ALLOCATE(ZWRK(INJ))
         IC=IC+1
         
         IF (LINITP) CFILEIN = 'PREP.txt'
         CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                            'ASCII ',CMASK_SAVE,'SURF  ','READ ')
         CALL POSNAM(NUNIT,CMASK//' '//HREC,GFOUND,NLUOUT)
         IF (.NOT.GFOUND .AND. LINITP)THEN
           ! Search now in PGD file
           CALL END_IO_SURF_n('ASCII ')
           CFILEIN = 'PGD.txt'
           CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                            'ASCII ',CMASK_SAVE,'SURF  ','READ ')
           CALL POSNAM(NUNIT,CMASK//' '//HREC,GFOUND,NLUOUT)
         ENDIF
         IF (.NOT.GFOUND) CALL ERR_STOP(HREC,CFILEIN,NLUOUT)

         READ(NUNIT,FMT='(A50)') NOM_ARTICLE
         READ(NUNIT,FMT='(A50)') YCOMMENT
         READ(NUNIT,FMT=*,ERR=100) ZWRK
         ALLOCATE(ZLOC(INJ))
         ZLOC(:)=ZWRK(:)
         WHERE(ZLOC(:)==999.) ZLOC(:)=-999.

         PRINT*,CMASK,' ',HREC,' ','MINVAL = ',MINVAL(ZWRK(:)),&
                     ' MAXVAL = ',MAXVAL(ZLOC(:))
         DEALLOCATE(ZLOC)

         IWFIELD=1
         IF (PATCHFLAG == '+') THEN
           DO JPATCH=1,IPATCH
             WRITE(YPAS,'(I2)') JPATCH
             YLVL=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
             OPEN(UNIT=30,FILE=TRIM(HREC)//'_p'//TRIM(YLVL)//'.dat',FORM='FORMATTED')
             IBEG=INI*(JPATCH-1)+1
             IEND=INI*JPATCH
             IF (LGEO) THEN
               DO IP=1,INI
                 !IF (ZWRK(IP)/=XUNDEF) THEN
                 WRITE(30,*)XLON(IP),XLAT(IP),ZWRK(INI*(JPATCH-1)+IP)
                 !ENDIF
               ENDDO
             ELSE
               WRITE(30,*)NOM_ARTICLE
               WRITE(30,*)YCOMMENT
               DO IP=1,INI
                 WRITE(30,*)ZWRK(IP)
               ENDDO
               !WRITE(30,*)INI
               !WRITE(30,'(60F16.8)')ZWRK(IBEG:IEND)
             ENDIF
             CLOSE(30)
           ENDDO
         ELSE
           OPEN(UNIT=30,FILE=TRIM(HREC)//'.dat',FORM='FORMATTED')
           IF (LGEO) THEN
             DO IP=1,INI
               !IF (ZWRK(IP)/=XUNDEF) THEN
               WRITE(30,*)XLAT(IP),XLON(IP),ZWRK(IP)
               !ENDIF
             ENDDO
           ELSE
             WRITE(30,*)NOM_ARTICLE
             WRITE(30,*)YCOMMENT
             DO IP=1,INI
               WRITE(30,*)ZWRK(IP)
             ENDDO
             !WRITE(30,*)INI
             !WRITE(30,'(60F16.8)')ZWRK(:)
           ENDIF
           CLOSE(30)
         ENDIF

         DEALLOCATE(ZWRK)
         CALL END_IO_SURF_n('ASCII ')

       ENDDO

       CLOSE(NLUOUT)
       CALL DEALLOC_SURFEX
       IF (LHOOK) CALL DR_HOOK('SXPOST',1,ZHOOK_HANDLE)


       STOP
 100   CONTINUE
       WRITE(NLUOUT,*) ' '
       WRITE(NLUOUT,*) ' ERROR WHEN READING ARTICLE',HREC
       WRITE(NLUOUT,*) ' '
       CLOSE(NLUOUT)
       CALL DEALLOC_SURFEX
       IF (LHOOK) CALL DR_HOOK('SXPOST',1,ZHOOK_HANDLE)

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
       CALL ABOR1_SFX('SXPOST: ARTICLE '//HREC//' NOT FOUND')
       IF (LHOOK) CALL DR_HOOK('ERR_STOP',1,ZHOOK_HANDLE)
       END SUBROUTINE ERR_STOP

       END PROGRAM SXPOST
