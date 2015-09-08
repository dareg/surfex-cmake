!     ###########
      PROGRAM PGD
!     ###########
!!
!!    PURPOSE
!!    -------
!!   This program prepares the physiographic data fields.
!!
!!    METHOD
!!    ------
!!   
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    F. Mereyde                  Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original     21/07/95
!!    Modification 26/07/95       Treatment of orography and subgrid-scale
!!                                orography roughness length (V. Masson)
!!    Modification 22/05/96       Variable CSTORAGE_TYPE (V. Masson)
!!    Modification 25/05/96       Modification of splines, correction on z0rel
!!                                and set limits for some surface varaibles
!!    Modification 12/06/96       Treatment of a rare case for ZPGDZ0EFF (Masson)
!!    Modification 22/11/96       removes the filtering. It will have to be 
!!                                performed in ADVANCED_PREP_PGD (Masson)
!!    Modification 15/03/99       **** MAJOR MODIFICATION **** (Masson)
!!                                PGD fields are now defined from the cover
!!                                type fractions in the grid meshes
!!                                User can still include its own data, and
!!                                even additional (dummy) fields
!!    Modificatio 06/00           patch approach, for vegetation related variable (Solmon/Masson)
!                                  averaging is performed on subclass(=patch) of nature
!!                08/03/01        add chemical emission treatment (D.Gazen)
!!    Modification 15/10/01       allow namelists in different orders (I.Mallet)
!!    Modification    07/11       new routine write_pgd_surf_atmn.F90 for writing PGD field (B.Decharme)
!!                                flag_update now in write_pgd_surf_atmn.F90 (B.Decharme)
!!
!!
!!                   ################################
!!    13/10/03       EXTERNALIZED VERSION (V. Masson)
!!                   ################################
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_DATA_BEM_n, ONLY : DTB => DATA_BEM
USE MODD_DATA_TEB_GARDEN_n, ONLY : DTGD => DATA_TEB_GARDEN
USE MODD_DATA_TEB_GREENROOF_n, ONLY : DTGR => DATA_TEB_GREENROOF
USE MODD_TEB_GARDEN_PGD_n, ONLY : TGDP => TEB_GARDEN_PGD
USE MODD_TEB_GREENROOF_PGD_n, ONLY : TGRP => TEB_GREENROOF_PGD
USE MODD_TEB_GRID_n, ONLY : TG => TEB_GRID
USE MODD_TEB_IRRIG_n, ONLY : TIR => TEB_IRRIG
!
USE MODD_DATA_ISBA_n, ONLY : DTI => DATA_ISBA
USE MODD_DUMMY_SURF_FIELDS_n, ONLY : DUU => DUMMY_SURF_FIELDS
USE MODD_FLAKE_GRID_n, ONLY : FG => FLAKE_GRID
USE MODD_ISBA_GRID_n, ONLY : IG => ISBA_GRID
USE MODD_SEAFLUX_GRID_n, ONLY : SG => SEAFLUX_GRID
USE MODD_SURF_ATM_SSO_n, ONLY : USS => SURF_ATM_SSO
USE MODD_WATFLUX_GRID_n, ONLY : WG => WATFLUX_GRID
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
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
!
USE MODD_SURF_ATM_GRID_n, ONLY : UG => SURF_ATM_GRID
!
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
USE MODD_IO_BUFF_n, ONLY : IOB => IO_BUFF
!
USE MODD_SURFEX_OMP, ONLY : NWORK, NWORK2, XWORK, XWORK2, XWORK3, &
                            NWORK_FULL, NWORK2_FULL, XWORK_FULL, XWORK2_FULL
!
USE MODD_IO_SURF_ASC
USE MODD_IO_SURF_FA
USE MODD_IO_SURF_LFI
USE MODD_IO_SURF_NC
USE MODD_SURF_CONF
USE MODI_OPEN_NAMELIST
USE MODI_CLOSE_NAMELIST
!      
USE MODI_GET_LONLAT_n
!
USE MODI_GOTO_SURFEX
USE MODI_IO_BUFF_CLEAN_n
USE MODI_PGD_OROG_FILTER
USE MODI_PGD_SURF_ATM
USE MODI_PGD_GRID_SURF_ATM
USE MODI_SPLIT_GRID
USE MODI_WRITE_HEADER_FA
USE MODI_WRITE_HEADER_MNH
USE MODI_WRITE_PGD_SURF_ATM_n
!
USE MODE_POS_SURF
!
USE MODN_IO_OFFLINE
USE MODN_WRITE_SURF_ATM
!
USE MODD_IO_SURF_FA, ONLY : LFANOCOMPACT
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ALLOC_SURFEX
USE MODI_DEALLOC_SURFEX
USE MODI_GET_LUOUT
!
IMPLICIT NONE
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER            :: ILUOUT
INTEGER            :: ILUNAM
LOGICAL            :: GFOUND
!
 CHARACTER(LEN=28)  :: YLUOUT    ='LISTING_PGD'   ! name of the listing
!
INTEGER            :: INW, JNW
INTEGER            :: IRET      
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PGD',0,ZHOOK_HANDLE)
!
 CALL ALLOC_SURFEX(1)
CSOFTWARE='PGD    '
 CALL GOTO_SURFEX(1,.TRUE.)
!
!*    1.      Set default names and parallelized I/O
!             --------------------------------------
!
 CALL GET_LUOUT('ASCII ',ILUOUT)
CLUOUT_LFI =  ADJUSTL(ADJUSTR(YLUOUT)//'.txt')
OPEN(UNIT=ILUOUT,FILE=ADJUSTL(ADJUSTR(YLUOUT)//'.txt'),FORM='FORMATTED',ACTION='WRITE')
!
!     1.3     output file name read in namelist
!             ---------------------------------
 CALL OPEN_NAMELIST('ASCII ',ILUNAM,CNAMELIST)
 CALL POSNAM(ILUNAM,'NAM_IO_OFFLINE',GFOUND)
IF (GFOUND) READ (UNIT=ILUNAM,NML=NAM_IO_OFFLINE)
 CALL POSNAM(ILUNAM,'NAM_WRITE_SURF_ATM',GFOUND)
IF (GFOUND) READ (UNIT=ILUNAM,NML=NAM_WRITE_SURF_ATM)
 CALL CLOSE_NAMELIST('ASCII ',ILUNAM)
!
CFILEOUT     = ADJUSTL(ADJUSTR(CPGDFILE)//'.txt')      ! output of PGD program
CFILEOUT_FA  = ADJUSTL(ADJUSTR(CPGDFILE)//'.fa')
CFILEOUT_LFI = CPGDFILE
CFILEOUT_NC  = ADJUSTL(ADJUSTR(CPGDFILE)//'.nc')
!
!*    2.      Preparation of surface physiographic fields
!             -------------------------------------------
!
 CALL PGD_GRID_SURF_ATM(IOB, &
                        UG, U, &
                        CSURF_FILETYPE,'                            ','      ',.FALSE.)
!
 CALL SPLIT_GRID(UG, U, &
                 'OFFLIN')
!
 CALL PGD_SURF_ATM(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTB, DTCO, &
                               DTI, DTS, DTGD, DTGR, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, &
                               DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, DUU, FG, F, FSB, &
                               GB, IOB, ICP, IG, I, O, SG, S, SSB, UG, U, &
                               USS, SV, TCP, TGD, TGDO, TGDP, TGR, TGRO, TGRP, TG, TIR, &
                               T, TOP, TVG, WG, W, WSB, &
                   CSURF_FILETYPE,'                            ','      ',.FALSE.)
!
 CALL PGD_OROG_FILTER(U, &
                      CSURF_FILETYPE)
!
!*    3.      writing of surface physiographic fields
!             ---------------------------------------
!
!* building of the header for the opening of the file in case of Arpege file
IF (CSURF_FILETYPE=='FA    ') THEN
  LFANOCOMPACT = .TRUE.
  CALL WRITE_HEADER_FA(IOB, UG, &
                       CSURF_FILETYPE,'PGD') 
END IF
!
LDEF = .TRUE.
!
INW = 1
IF (CSURF_FILETYPE=="NC    ") INW = 2
!
DO JNW = 1,INW
  !
  IF (LWRITE_COORD) CALL GET_LONLAT_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTS, DTT, &
                                DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, &
                                DGUT, DGW, F, FSB, GB, ICP, I, O, S, SSB, SV, &
                                TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                                      DTCO, IOB, UG, U, &
                                      CSURF_FILETYPE)
  !
  !* writing of the fields
 CALL IO_BUFF_CLEAN_n(IOB)
  
  ! FLAG_UPDATE now in WRITE_PGD_SURF_ATM_n
  CALL WRITE_PGD_SURF_ATM_n(DTB, DTGD, DTGR, TGDP, TGRP, TG, TIR, &
                            BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTCO, DTI, &
                                       DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, &
                                       DGU, DGT, DGUT, DGW, DUU, FG, F, FSB, GB, IOB, ICP, &
                                       IG, I, O, SG, S, SSB, UG, U, USS, SV, TCP, &
                                       TGD, TGDO, TGR, TGRO, T, TOP, TVG, WG, W, WSB, &
                            CSURF_FILETYPE)
  !
  LDEF = .FALSE.
  CALL IO_BUFF_CLEAN_n(IOB)  
  !
ENDDO
!
!* closes the file
IF (CSURF_FILETYPE=='FA    ') THEN
  CALL FAIRME(IRET,NUNIT_FA,'UNKNOWN')
END IF
!
!* add informations in the file
IF (CSURF_FILETYPE=='LFI   ' .AND. LMNH_COMPATIBLE) CALL WRITE_HEADER_MNH
!
!*    3.     Close parallelized I/O
!            ----------------------
!
WRITE(ILUOUT,*) ' '
WRITE(ILUOUT,*) '    ----------------------'
WRITE(ILUOUT,*) '    | PGD ENDS CORRECTLY |'
WRITE(ILUOUT,*) '    ----------------------'
!
WRITE(*,*) ' '
WRITE(*,*) '    ----------------------'
WRITE(*,*) '    | PGD ENDS CORRECTLY |'
WRITE(*,*) '    ----------------------'
      !
CLOSE(ILUOUT)
 CALL DEALLOC_SURFEX
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
IF (LHOOK) CALL DR_HOOK('PGD',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END PROGRAM PGD
