!     #################################################################################
SUBROUTINE INIT_PGD_SURF_ATM(HPROGRAM,HINIT,HATMFILE,HATMFILETYPE, &
                               KYEAR, KMONTH, KDAY, PTIME            )  
!     #################################################################################
!
!!****  *INIT_PGD_SURF_ATM* - Call surface initialization for PGD fields only
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     V. Masson 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      B. Decharme  04/2013 new coupling variables
!!------------------------------------------------------------------
!
!
!
!
USE MODD_AGRI_n, ONLY : AG => AGRI
USE MODD_BEM_n, ONLY : B => BEM
USE MODD_BEM_OPTION_n, ONLY : BOP => BEM_OPTIONS
USE MODD_BLD_DESCRIPTION_n, ONLY : BDD => BLD_DESC
USE MODD_CH_EMIS_FIELD_n, ONLY : CHE => CH_EMIS_FIELD
USE MODD_CH_FLAKE_n, ONLY : CHF => CH_FLAKE
USE MODD_CH_ISBA_n, ONLY : CHI => CH_ISBA
USE MODD_CH_SEAFLUX_n, ONLY : CHS => CH_SEAFLUX
USE MODD_CH_SNAP_n, ONLY : CHN => CH_EMIS_SNAP
USE MODD_CH_SURF_n, ONLY : CHU => CH_SURF
USE MODD_CH_TEB_n, ONLY : CHT => CH_TEB
USE MODD_CH_WATFLUX_n, ONLY : CHW => CH_WATFLUX
USE MODD_DATA_BEM_n, ONLY : DTB => DATA_BEM
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
USE MODD_DATA_ISBA_n, ONLY : DTI => DATA_ISBA
USE MODD_DATA_SEAFLUX_n, ONLY : DTS => DATA_SEAFLUX
USE MODD_DATA_TEB_GARDEN_n, ONLY : DTGD => DATA_TEB_GARDEN
USE MODD_DATA_TEB_GREENROOF_n, ONLY : DTGR => DATA_TEB_GREENROOF
USE MODD_DATA_TEB_n, ONLY : DTT => DATA_TEB
USE MODD_DATA_TSZ0_n, ONLY : DTZ => DATA_TSZ0
USE MODD_DIAG_CUMUL_TEB_n, ONLY : DGCT => DIAG_CUMUL_TEB
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DGEI => DIAG_EVAP_ISBA
USE MODD_DIAG_FLAKE_n, ONLY : DGF => DIAG_FLAKE
USE MODD_DIAG_IDEAL_n, ONLY : DGL => DIAG_IDEAL
USE MODD_DIAG_ISBA_n, ONLY : DGI => DIAG_ISBA
USE MODD_DIAG_MISC_FLAKE_n, ONLY : DGMF => DIAG_MISC_FLAKE
USE MODD_DIAG_MISC_ISBA_n, ONLY : DGMI => DIAG_MISC_ISBA
USE MODD_DIAG_MISC_TEB_n, ONLY : DGMT => DIAG_MISC_TEB
USE MODD_DIAG_MISC_TEB_OPTION_n, ONLY : DGMTO => DIAG_MISC_TEB_OPTIONS
USE MODD_DIAG_OCEAN_n, ONLY : DGO => DIAG_OCEAN
USE MODD_DIAG_SEAFLUX_n, ONLY : DGS => DIAG_SEAFLUX
USE MODD_DIAG_SEAICE_n, ONLY : DGSI => DIAG_SEAICE
USE MODD_DIAG_SURF_ATM_n, ONLY : DGU => DIAG_SURF_ATM
USE MODD_DIAG_TEB_n, ONLY : DGT => DIAG_TEB
USE MODD_DIAG_UTCI_TEB_n, ONLY : DGUT => DIAG_UTCI_TEB
USE MODD_DIAG_WATFLUX_n, ONLY : DGW => DIAG_WATFLUX
USE MODD_DST_n, ONLY : DST => DST
USE MODD_DUMMY_SURF_FIELDS_n, ONLY : DUU => DUMMY_SURF_FIELDS
USE MODD_FLAKE_GRID_n, ONLY : FG => FLAKE_GRID
USE MODD_FLAKE_n, ONLY : F => FLAKE
USE MODD_FLAKE_SBL_n, ONLY : FSB => FLAKE_SBL
USE MODD_GR_BIOG_GARDEN_n, ONLY : GBGD => GR_BIOG_GARDEN
USE MODD_GR_BIOG_GREENROOF_n, ONLY : GBGR => GR_BIOG_GREENROOF
USE MODD_GR_BIOG_n, ONLY : GB => GR_BIOG
USE MODD_IO_BUFF_n, ONLY : IOB => IO_BUFF
USE MODD_ISBA_CANOPY_n, ONLY : ICP => ISBA_CANOPY
USE MODD_ISBA_GRID_n, ONLY : IG => ISBA_GRID
USE MODD_ISBA_n, ONLY : I => ISBA
USE MODD_OCEAN_n, ONLY : O => OCEAN
USE MODD_OCEAN_REL_n, ONLY : OR => OCEAN_REL
USE MODD_SEAFLUX_GRID_n, ONLY : SG => SEAFLUX_GRID
USE MODD_SEAFLUX_n, ONLY : S => SEAFLUX
USE MODD_SEAFLUX_SBL_n, ONLY : SSB => SEAFLUX_SBL
USE MODD_SLT_n, ONLY : SLT => SLT
USE MODD_SSO_CANOPY_n, ONLY : SSCP => SSO_CANOPY
USE MODD_SURF_ATM_GRID_n, ONLY : UG => SURF_ATM_GRID
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
USE MODD_SURF_ATM_SSO_n, ONLY : USS => SURF_ATM_SSO
USE MODD_SV_n, ONLY : SV => SV
USE MODD_TEB_CANOPY_n, ONLY : TCP => TEB_CANOPY
USE MODD_TEB_GARDEN_n, ONLY : TGD => TEB_GARDEN
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TGDO => TEB_GARDEN_OPTIONS
USE MODD_TEB_GARDEN_PGD_EVOL_n, ONLY : TGDPE => TEB_GARDEN_PGD_EVOL
USE MODD_TEB_GARDEN_PGD_n, ONLY : TGDP => TEB_GARDEN_PGD
USE MODD_TEB_GREENROOF_n, ONLY : TGR => TEB_GREENROOF
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TGRO => TEB_GREENROOF_OPTIONS
USE MODD_TEB_GREENROOF_PGD_EVOL_n, ONLY : TGRPE => TEB_GREENROOF_PGD_EVOL
USE MODD_TEB_GREENROOF_PGD_n, ONLY : TGRP => TEB_GREENROOF_PGD
USE MODD_TEB_GRID_n, ONLY : TG => TEB_GRID
USE MODD_TEB_IRRIG_n, ONLY : TIR => TEB_IRRIG
USE MODD_TEB_n, ONLY : T => TEB
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS
USE MODD_TEB_PANEL_n, ONLY : TPN => TEB_PANEL
USE MODD_TEB_VEG_n, ONLY : TVG => TEB_VEG_OPTIONS
USE MODD_WATFLUX_GRID_n, ONLY : WG => WATFLUX_GRID
USE MODD_WATFLUX_n, ONLY : W => WATFLUX
USE MODD_WATFLUX_SBL_n, ONLY : WSB => WATFLUX_SBL
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_INIT_SURF_ATM_n
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=3),   INTENT(IN)  :: HINIT     ! fields to initialize 'ALL', 'PRE', 'PGD'
 CHARACTER(LEN=28), INTENT(IN)   :: HATMFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),  INTENT(IN)   :: HATMFILETYPE! type of the Atmospheric file
INTEGER,           INTENT(IN)   :: KYEAR       ! year
INTEGER,           INTENT(IN)   :: KMONTH      ! month
INTEGER,           INTENT(IN)   :: KDAY        ! day
REAL,              INTENT(IN)   :: PTIME       ! time
!
!
!*      0.2    declarations of local variables
!
 CHARACTER(LEN=6), DIMENSION(0)  :: YSV       ! name of all scalar variables
REAL,             DIMENSION(0)  :: ZCO2      ! CO2 concentration (kg/m3)
REAL,             DIMENSION(0)  :: ZRHOA     ! air density (kg/m3)
REAL,             DIMENSION(0)  :: ZZENITH   ! solar zenithal angle
REAL,             DIMENSION(0)  :: ZAZIM     ! solar azimuthal angle (rad from N, clock)
REAL,             DIMENSION(1)  :: ZSW_BANDS ! middle wavelength of each band
REAL,             DIMENSION(0,1):: ZDIR_ALB  ! direct albedo for each band
REAL,             DIMENSION(0,1):: ZSCA_ALB  ! diffuse albedo for each band
REAL,             DIMENSION(0)  :: ZEMIS     ! emissivity
REAL,             DIMENSION(0)  :: ZTSRAD    ! radiative temperature
REAL,             DIMENSION(0)  :: ZTSURF    ! radiative temperature
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
!* initialization of PGD fields of output domain
!
IF (LHOOK) CALL DR_HOOK('INIT_PGD_SURF_ATM',0,ZHOOK_HANDLE)
 CALL INIT_SURF_ATM_n(AG, B, BOP, BDD, CHE, CHF, CHI, CHS, CHN, CHU, CHT, &
                            CHW, DTB, DTCO, DTI, DTS, DTGD, DTGR, DTT, DTZ, DGCT, DGEI, &
                            DGF, DGL, DGI, DGMF, DGMI, DGMT, DGMTO, DGO, DGS, DGSI, DGU, &
                            DGT, DGUT, DGW, DST, DUU, FG, F, FSB, GBGD, GBGR, GB, &
                            IOB, ICP, IG, I, O, OR, SG, S, SSB, SLT, SSCP, &
                            UG, U, USS, SV, TCP, TGD, TGDO, TGDPE, TGDP, TGR, TGRO, &
                            TGRPE, TGRP, TG, TIR, T, TOP, TPN, TVG, WG, W, WSB, &
                            HPROGRAM,HINIT,.FALSE.,                     &
                      0,0,1,                                     &
                      YSV,ZCO2,ZRHOA,                            &
                      ZZENITH,ZAZIM,ZSW_BANDS,ZDIR_ALB,ZSCA_ALB, &
                      ZEMIS,ZTSRAD,ZTSURF,                       &
                      KYEAR, KMONTH, KDAY, PTIME,                &
                      HATMFILE,HATMFILETYPE, 'OK'                )  
IF (LHOOK) CALL DR_HOOK('INIT_PGD_SURF_ATM',1,ZHOOK_HANDLE)

!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE INIT_PGD_SURF_ATM
