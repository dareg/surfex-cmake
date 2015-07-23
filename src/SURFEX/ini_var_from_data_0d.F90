!     #########
      SUBROUTINE INI_VAR_FROM_DATA_0D (DTCO, IOB, UG, U, USS, &
                                       HPROGRAM, HATYPE,  HNAME, HTYPE, HFNAM, &
                                        HFTYP, PUNIF, PFIELD, OPRESENT)
!     ##############################################################
!
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!
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
!!    S. Faroux        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    16/11/10
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
!
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
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SURF_ATM_SSO_n, ONLY : SURF_ATM_SSO_t
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
!
USE MODD_PGDWORK,       ONLY : CATYPE
!
USE MODI_PGD_FIELD
USE MODI_READ_FROM_SURFEX_FILE
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SURF_ATM_SSO_t), INTENT(INOUT) :: USS
!
 CHARACTER(LEN=6), INTENT(IN) :: HPROGRAM
 CHARACTER(LEN=3), INTENT(IN) :: HATYPE
 CHARACTER(LEN=*), INTENT(IN) :: HNAME
 CHARACTER(LEN=3), INTENT(IN) :: HTYPE
 CHARACTER(LEN=28), INTENT(IN) :: HFNAM
 CHARACTER(LEN=6), INTENT(IN) :: HFTYP
REAL, INTENT(IN) :: PUNIF
REAL, DIMENSION(:), INTENT(OUT) :: PFIELD
LOGICAL, INTENT(OUT) :: OPRESENT
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!

!-------------------------------------------------------------------------------
!
!*    1.      Initializations
!             ---------------
!
IF (LHOOK) CALL DR_HOOK('INI_VAR_FROM_DATA_0D',0,ZHOOK_HANDLE)
!
PFIELD(:)=XUNDEF
OPRESENT=.FALSE.
!
IF (HFTYP.EQ.'FA    ' .OR. HFTYP.EQ.'ASCII ' .OR. HFTYP.EQ.'LFI   ') THEN

  OPRESENT=.TRUE.
  SELECT CASE (HTYPE)
    CASE ('LAN')
      CALL READ_FROM_SURFEX_FILE(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTCO, DTS, &
                                   DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, &
                                   DGT, DGUT, DGW, F, FSB, GB, ICP, I, O, S, SSB, &
                                   UG, U, SV, TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, &
                                   W, WSB, &
                                 IOB, &
                                 HFTYP,HFNAM,'SURF  ','      ',PFIELD)
    CASE ('TWN')
      CALL READ_FROM_SURFEX_FILE(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTCO, DTS, &
                                   DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, &
                                   DGT, DGUT, DGW, F, FSB, GB, ICP, I, O, S, SSB, &
                                   UG, U, SV, TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, &
                                   W, WSB, &
                                 IOB, &
                                 HFTYP,HFNAM,'TOWN  ','      ',PFIELD)              
    CASE ('NAT')
      CALL READ_FROM_SURFEX_FILE(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTCO, DTS, &
                                   DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, &
                                   DGT, DGUT, DGW, F, FSB, GB, ICP, I, O, S, SSB, &
                                   UG, U, SV, TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, &
                                   W, WSB, &
                                 IOB, &
                                 HFTYP,HFNAM,'NATURE','      ',PFIELD)              
    CASE ('SEA')
      CALL READ_FROM_SURFEX_FILE(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTCO, DTS, &
                                   DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, &
                                   DGT, DGUT, DGW, F, FSB, GB, ICP, I, O, S, SSB, &
                                   UG, U, SV, TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, &
                                   W, WSB, &
                                 IOB, &
                                 HFTYP,HFNAM,'SEA   ','      ',PFIELD)              
    CASE ('WAT')
      CALL READ_FROM_SURFEX_FILE(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTCO, DTS, &
                                   DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, &
                                   DGT, DGUT, DGW, F, FSB, GB, ICP, I, O, S, SSB, &
                                   UG, U, SV, TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, &
                                   W, WSB, &
                                 IOB, &
                                 HFTYP,HFNAM,'WATER ','      ',PFIELD)              
   END SELECT

ELSE

  CATYPE = HATYPE
  CALL PGD_FIELD(DTCO, UG, U, USS, &
                 HPROGRAM,HNAME,HTYPE,HFNAM,HFTYP,PUNIF,PFIELD(:),OPRESENT=OPRESENT)

ENDIF
!
IF (LHOOK) CALL DR_HOOK('INI_VAR_FROM_DATA_0D',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INI_VAR_FROM_DATA_0D

