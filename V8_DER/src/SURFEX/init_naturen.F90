!     #############################################################
      SUBROUTINE INIT_NATURE_n (AG, BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTCO, &
                                  DTI, DTS, DTT, DTZ, DGEI, DGF, DGL, DGI, DGMI, DGMTO, DGO, &
                                  DGS, DGSI, DGU, DGT, DGUT, DGW, DST, F, FSB, GB, IOB, &
                                  ICP, IG, I, O, S, SSB, SLT, UG, U, SV, TCP, &
                                  TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                                  HPROGRAM,HINIT,OLAND_USE,                    &
                                   KI,KSV,KSW,                                &
                                   HSV,PCO2,PRHOA,                            &
                                   PZENITH,PAZIM,PSW_BANDS,PDIR_ALB,PSCA_ALB, &
                                   PEMIS,PTSRAD,PTSURF,                       &
                                   KYEAR, KMONTH,KDAY, PTIME,                 &
                                   HATMFILE,HATMFILETYPE,                     &
                                   HTEST                                      )  
!     #############################################################
!
!!****  *INIT_NATURE_n* - routine to choose initialization of vegetation scheme
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    27/09/96 
!!       V.Masson   18/08/97 call to fmread directly with dates and strings
!!       V.Masson   15/03/99 new PGD treatment with COVER types
!        F.Solmon  06/00   adaptation for patch approach
!!      B. Decharme  04/2013 new coupling variables
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
!
!
USE MODD_AGRI_n, ONLY : AGRI_t
USE MODD_BEM_OPTION_n, ONLY : BEM_OPTIONS_t
USE MODD_BLD_DESCRIPTION_n, ONLY : BLD_DESC_t
USE MODD_CH_EMIS_FIELD_n, ONLY : CH_EMIS_FIELD_t
USE MODD_CH_ISBA_n, ONLY : CH_ISBA_t
USE MODD_CH_SEAFLUX_n, ONLY : CH_SEAFLUX_t
USE MODD_CH_SNAP_n, ONLY : CH_EMIS_SNAP_t
USE MODD_CH_SURF_n, ONLY : CH_SURF_t
USE MODD_CH_TEB_n, ONLY : CH_TEB_t
USE MODD_CH_WATFLUX_n, ONLY : CH_WATFLUX_t
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_DATA_SEAFLUX_n, ONLY : DATA_SEAFLUX_t
USE MODD_DATA_TEB_n, ONLY : DATA_TEB_t
USE MODD_DATA_TSZ0_n, ONLY : DATA_TSZ0_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_FLAKE_n, ONLY : DIAG_FLAKE_t
USE MODD_DIAG_IDEAL_n, ONLY : DIAG_IDEAL_t
USE MODD_DIAG_ISBA_n, ONLY : DIAG_ISBA_t
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t
USE MODD_DIAG_MISC_TEB_OPTION_n, ONLY : DIAG_MISC_TEB_OPTIONS_t
USE MODD_DIAG_OCEAN_n, ONLY : DIAG_OCEAN_t
USE MODD_DIAG_SEAFLUX_n, ONLY : DIAG_SEAFLUX_t
USE MODD_DIAG_SEAICE_n, ONLY : DIAG_SEAICE_t
USE MODD_DIAG_SURF_ATM_n, ONLY : DIAG_SURF_ATM_t
USE MODD_DIAG_TEB_n, ONLY : DIAG_TEB_t
USE MODD_DIAG_UTCI_TEB_n, ONLY : DIAG_UTCI_TEB_t
USE MODD_DIAG_WATFLUX_n, ONLY : DIAG_WATFLUX_t
USE MODD_DST_n, ONLY : DST_t
USE MODD_FLAKE_n, ONLY : FLAKE_t
USE MODD_FLAKE_SBL_n, ONLY : FLAKE_SBL_t
USE MODD_GR_BIOG_n, ONLY : GR_BIOG_t
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
USE MODD_ISBA_CANOPY_n, ONLY : ISBA_CANOPY_t
USE MODD_ISBA_GRID_n, ONLY : ISBA_GRID_t
USE MODD_ISBA_n, ONLY : ISBA_t
USE MODD_OCEAN_n, ONLY : OCEAN_t
USE MODD_SEAFLUX_n, ONLY : SEAFLUX_t
USE MODD_SEAFLUX_SBL_n, ONLY : SEAFLUX_SBL_t
USE MODD_SLT_n, ONLY : SLT_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SV_n, ONLY : SV_t
USE MODD_TEB_CANOPY_n, ONLY : TEB_CANOPY_t
USE MODD_TEB_GARDEN_n, ONLY : TEB_GARDEN_t
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TEB_GARDEN_OPTIONS_t
USE MODD_TEB_GREENROOF_n, ONLY : TEB_GREENROOF_t
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TEB_GREENROOF_OPTIONS_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_OPTIONS_t
USE MODD_WATFLUX_n, ONLY : WATFLUX_t
USE MODD_WATFLUX_SBL_n, ONLY : WATFLUX_SBL_t
!
USE MODD_CSTS,       ONLY : XTT
!
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_INIT_IDEAL_FLUX
!
USE MODI_INIT_ISBA_n
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
!
TYPE(AGRI_t), INTENT(INOUT) :: AG
TYPE(BEM_OPTIONS_t), INTENT(INOUT) :: BOP
TYPE(BLD_DESC_t), INTENT(INOUT) :: BDD
TYPE(CH_EMIS_FIELD_t), INTENT(INOUT) :: CHE
TYPE(CH_ISBA_t), INTENT(INOUT) :: CHI
TYPE(CH_SEAFLUX_t), INTENT(INOUT) :: CHS
TYPE(CH_EMIS_SNAP_t), INTENT(INOUT) :: CHN
TYPE(CH_SURF_t), INTENT(INOUT) :: CHU
TYPE(CH_TEB_t), INTENT(INOUT) :: CHT
TYPE(CH_WATFLUX_t), INTENT(INOUT) :: CHW
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTI
TYPE(DATA_SEAFLUX_t), INTENT(INOUT) :: DTS
TYPE(DATA_TEB_t), INTENT(INOUT) :: DTT
TYPE(DATA_TSZ0_t), INTENT(INOUT) :: DTZ
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEI
TYPE(DIAG_FLAKE_t), INTENT(INOUT) :: DGF
TYPE(DIAG_IDEAL_t), INTENT(INOUT) :: DGL
TYPE(DIAG_ISBA_t), INTENT(INOUT) :: DGI
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DGMI
TYPE(DIAG_MISC_TEB_OPTIONS_t), INTENT(INOUT) :: DGMTO
TYPE(DIAG_OCEAN_t), INTENT(INOUT) :: DGO
TYPE(DIAG_SEAFLUX_t), INTENT(INOUT) :: DGS
TYPE(DIAG_SEAICE_t), INTENT(INOUT) :: DGSI
TYPE(DIAG_SURF_ATM_t), INTENT(INOUT) :: DGU
TYPE(DIAG_TEB_t), INTENT(INOUT) :: DGT
TYPE(DIAG_UTCI_TEB_t), INTENT(INOUT) :: DGUT
TYPE(DIAG_WATFLUX_t), INTENT(INOUT) :: DGW
TYPE(DST_t), INTENT(INOUT) :: DST
TYPE(FLAKE_t), INTENT(INOUT) :: F
TYPE(FLAKE_SBL_t), INTENT(INOUT) :: FSB
TYPE(GR_BIOG_t), INTENT(INOUT) :: GB
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
TYPE(ISBA_CANOPY_t), INTENT(INOUT) :: ICP
TYPE(ISBA_GRID_t), INTENT(INOUT) :: IG
TYPE(ISBA_t), INTENT(INOUT) :: I
TYPE(OCEAN_t), INTENT(INOUT) :: O
TYPE(SEAFLUX_t), INTENT(INOUT) :: S
TYPE(SEAFLUX_SBL_t), INTENT(INOUT) :: SSB
TYPE(SLT_t), INTENT(INOUT) :: SLT
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SV_t), INTENT(INOUT) :: SV
TYPE(TEB_CANOPY_t), INTENT(INOUT) :: TCP
TYPE(TEB_GARDEN_t), INTENT(INOUT) :: TGD
TYPE(TEB_GARDEN_OPTIONS_t), INTENT(INOUT) :: TGDO
TYPE(TEB_GREENROOF_t), INTENT(INOUT) :: TGR
TYPE(TEB_GREENROOF_OPTIONS_t), INTENT(INOUT) :: TGRO
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_VEG_OPTIONS_t), INTENT(INOUT) :: TVG
TYPE(WATFLUX_t), INTENT(INOUT) :: W
TYPE(WATFLUX_SBL_t), INTENT(INOUT) :: WSB
!
 CHARACTER(LEN=6),                 INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=3),                 INTENT(IN)  :: HINIT     ! choice of fields to initialize
LOGICAL,                          INTENT(IN)  :: OLAND_USE ! 
INTEGER,                          INTENT(IN)  :: KI        ! number of points
INTEGER,                          INTENT(IN)  :: KSV       ! number of scalars
INTEGER,                          INTENT(IN)  :: KSW       ! number of short-wave spectral bands
 CHARACTER(LEN=6), DIMENSION(KSV), INTENT(IN)  :: HSV       ! name of all scalar variables
REAL,             DIMENSION(KI),  INTENT(IN)  :: PCO2      ! CO2 concentration (kg/m3)
REAL,             DIMENSION(KI),  INTENT(IN)  :: PRHOA     ! air density
REAL,             DIMENSION(KI),  INTENT(IN)  :: PZENITH   ! solar zenithal angle
REAL,             DIMENSION(KI),  INTENT(IN)  :: PAZIM     ! solar azimuthal angle (rad from N, clock)
REAL,             DIMENSION(KSW), INTENT(IN)  :: PSW_BANDS ! middle wavelength of each band
REAL,             DIMENSION(KI,KSW),INTENT(OUT) :: PDIR_ALB  ! direct albedo for each band
REAL,             DIMENSION(KI,KSW),INTENT(OUT) :: PSCA_ALB  ! diffuse albedo for each band
REAL,             DIMENSION(KI),  INTENT(OUT) :: PEMIS     ! emissivity
REAL,             DIMENSION(KI),  INTENT(OUT) :: PTSRAD    ! radiative temperature
REAL,             DIMENSION(KI),  INTENT(OUT) :: PTSURF    ! surface effective temperature         (K)
INTEGER,                          INTENT(IN)  :: KYEAR     ! current year (UTC)
INTEGER,                          INTENT(IN)  :: KMONTH    ! current month (UTC)
INTEGER,                          INTENT(IN)  :: KDAY      ! current day (UTC)
REAL,                             INTENT(IN)  :: PTIME     ! current time since
                                                          !  midnight (UTC, s)
!
 CHARACTER(LEN=28),                INTENT(IN)  :: HATMFILE    ! atmospheric file name
 CHARACTER(LEN=6),                 INTENT(IN)  :: HATMFILETYPE! atmospheric file type
 CHARACTER(LEN=2),                 INTENT(IN)  :: HTEST       ! must be equal to 'OK'
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
!*       2.     Selection of surface scheme
!               ---------------------------
!
IF (LHOOK) CALL DR_HOOK('INIT_NATURE_N',0,ZHOOK_HANDLE)
IF (U%CNATURE=='NONE  ') THEN
  PDIR_ALB=0.
  PSCA_ALB=0.
  PEMIS   =1.
  PTSRAD  =XTT
  PTSURF  =XTT
ELSE IF (U%CNATURE=='FLUX  ') THEN
  CALL INIT_IDEAL_FLUX(DGL, &
                       HPROGRAM,HINIT,KI,KSV,KSW,HSV,PCO2,PRHOA,     &
                           PZENITH,PAZIM,PSW_BANDS,PDIR_ALB,PSCA_ALB,  &
                           PEMIS,PTSRAD,PTSURF,'OK'                    )  
ELSE IF (U%CNATURE=='ISBA  ' .OR. U%CNATURE=='TSZ0') THEN
  CALL INIT_ISBA_n(AG, BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTCO, &
                            DTI, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, &
                            DGSI, DGU, DGT, DGUT, DGW, DST, F, FSB, GB, IOB, ICP, &
                            IG, I, O, S, SSB, SLT, UG, U, SV, TCP, TGD, &
                            TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                   HPROGRAM,HINIT,OLAND_USE,KI,KSV,KSW,HSV,PCO2,PRHOA,     &
                     PZENITH,PAZIM,PSW_BANDS,PDIR_ALB,PSCA_ALB,    &
                     PEMIS,PTSRAD,PTSURF,                          &
                     KYEAR,KMONTH,KDAY,PTIME,HATMFILE,HATMFILETYPE,&
                     'OK'                                          )  
END IF
IF (LHOOK) CALL DR_HOOK('INIT_NATURE_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INIT_NATURE_n
