!     #########
      SUBROUTINE READWRITE_EMIS_FIELD_n (U, &
                                         HPROGRAM)
!     #######################################################################
!
!-----------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
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
USE MODD_DIAG_SURF_ATM_n, ONLY : DGU => DIAG_SURF_ATM
!
USE MODD_IO_BUFF_n, ONLY : IOB => IO_BUFF
!
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODI_GET_LUOUT
USE MODI_INIT_IO_SURF_n
USE MODI_END_IO_SURF_n
USE MODI_READ_SURF
USE MODI_WRITE_SURF
!
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
!
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
 CHARACTER(LEN=6) :: HPROGRAM
!
!*       0.2   declarations of local variables
!
INTEGER             :: IRESP  ! I/O error code
 CHARACTER (LEN=16)  :: YRECFM ! article name
 CHARACTER (LEN=100) :: YCOMMENT ! comment
INTEGER             :: ILUOUT   ! Unit number for prints
INTEGER             :: JSPEC    ! Loop index for emission species
INTEGER             :: IEMISPEC_NBR    ! number of emitted chemical species
 CHARACTER(LEN=40)   :: YEMISPEC_NAME   ! species name
INTEGER             :: IEMISPEC_NTIMES ! number of emission times
 CHARACTER(LEN=3)    :: YSURF ! surface type
INTEGER,DIMENSION(:),ALLOCATABLE :: ITIMES ! emission times for a species
REAL, DIMENSION(:,:),ALLOCATABLE :: ZWORK  ! work array read in the file
!
INTEGER           :: IVERSION       ! version of surfex file being read
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('READWRITE_EMIS_FIELD_N',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
!-------------------------------------------------------------------------------
!
 CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                     HPROGRAM,'FULL  ','SURF  ','READ ')
!* ascendant compatibility
YRECFM='VERSION'
 CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,IVERSION,IRESP)
!
YRECFM='EMISFILE_NBR'
IF (IVERSION<4) YRECFM='EMISFILE_GR_NBR'
 CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,IEMISPEC_NBR,IRESP,YCOMMENT)
 CALL END_IO_SURF_n(HPROGRAM)
!
IF (IRESP/=0) THEN
  CALL ABOR1_SFX('READWRITE_EMIS_FIELDN: PROBLEM READING NUMBER OF 2D CHEMICAL EMISSION FIELDS')
END IF
!
 CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                     HPROGRAM,'FULL  ','SURF  ','WRITE')
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,IEMISPEC_NBR,IRESP,YCOMMENT)
 CALL END_IO_SURF_n(HPROGRAM)
!
!-------------------------------------------------------------------------------
!
 CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                     HPROGRAM,'FULL  ','SURF  ','READ ')
YRECFM='EMISPEC_NBR'
IF (IVERSION<4) YRECFM='EMISPEC_GR_NBR'
 CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,IEMISPEC_NBR,IRESP,YCOMMENT)
 CALL END_IO_SURF_n(HPROGRAM)
!
IF (IRESP/=0) THEN
  CALL ABOR1_SFX('READWRITE_EMIS_FIELDN: PROBLEM READING NUMBER OF EMITTED CHEMICAL SPECIES')
END IF
!
 CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                     HPROGRAM,'FULL  ','SURF  ','WRITE')
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,IEMISPEC_NBR,IRESP,YCOMMENT)
 CALL END_IO_SURF_n(HPROGRAM)
!
!-------------------------------------------------------------------------------
!
DO JSPEC=1,IEMISPEC_NBR
  CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                     HPROGRAM,'FULL  ','SURF  ','READ ')
  WRITE(YRECFM,'("EMISNAME",I3.3)') JSPEC
  CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,YEMISPEC_NAME,IRESP,YCOMMENT)
  CALL END_IO_SURF_n(HPROGRAM)
!
  IF (IRESP/=0) THEN
    CALL ABOR1_SFX('READWRITE_EMIS_FIELDN: PROBLEM WHEN READING THE NAME OF EMITTED CHEMICAL SPECIES '//YRECFM)
  END IF
  READ(YCOMMENT,'(A3,24x,I5)') YSURF, IEMISPEC_NTIMES
  !
  CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                     HPROGRAM,'FULL  ','SURF  ','WRITE')
  CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,YEMISPEC_NAME,IRESP,YCOMMENT)
  CALL END_IO_SURF_n(HPROGRAM)
!  
!-------------------------------------------------------------------------------
!
  ALLOCATE(ITIMES(IEMISPEC_NTIMES))
  ALLOCATE(ZWORK(U%NSIZE_FULL,IEMISPEC_NTIMES))
!
!-------------------------------------------------------------------------------
!
  CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                     HPROGRAM,'FULL  ','SURF  ','READ ')
  YRECFM='E_'//TRIM(YEMISPEC_NAME)
  CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,ZWORK,IRESP,YCOMMENT)
  CALL END_IO_SURF_n(HPROGRAM)
  !
  IF (IRESP/=0) THEN
    CALL ABOR1_SFX('READWRITE_EMIS_FIELDN: PROBLEM WHEN READING THE EMISSION DATA '//YRECFM)
  END IF
  !
  CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                     HPROGRAM,'FULL  ','SURF  ','WRITE')
  CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,ZWORK,IRESP,YCOMMENT)
  CALL END_IO_SURF_n(HPROGRAM)
!
!-------------------------------------------------------------------------------
!
  CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                     HPROGRAM,'FULL  ','SURF  ','READ ')
  WRITE(YRECFM,'("EMISTIMES",I3.3)') JSPEC
  CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,ITIMES,IRESP,YCOMMENT,'-')
  CALL END_IO_SURF_n(HPROGRAM)

  IF (IRESP/=0) THEN
    CALL ABOR1_SFX('READWRITE_EMIS_FIELDN: PROBLEM WHEN READING THE EMISSION TIMES '//YRECFM)
  END IF

  CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                     HPROGRAM,'FULL  ','SURF  ','WRITE')
  CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,YRECFM,ITIMES,IRESP,YCOMMENT,'-')
  CALL END_IO_SURF_n(HPROGRAM)
!
!-------------------------------------------------------------------------------
!
  DEALLOCATE(ITIMES)
  DEALLOCATE(ZWORK)
!
!-------------------------------------------------------------------------------
END DO
IF (LHOOK) CALL DR_HOOK('READWRITE_EMIS_FIELD_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE READWRITE_EMIS_FIELD_n
