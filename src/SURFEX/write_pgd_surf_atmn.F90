!     ####################################
      SUBROUTINE WRITE_PGD_SURF_ATM_n (BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTCO, DTI, &
                                       DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, &
                                       DGU, DGT, DGUT, DGW, DUU, FG, F, FSB, GB, IOB, ICP, &
                                       IG, I, O, SG, S, SSB, UG, U, USS, SV, TCP, &
                                       TGD, TGDO, TGR, TGRO, T, TOP, TVG, WG, W, WSB, &
                                       HPROGRAM)
!     ####################################
!
!!****  *WRITE_PGD_SURF_ATM_n* - routine to write pgd surface variables 
!!                               in their respective files or in file
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
!!      B. Decharme   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/2011 according to previous write_surf_atmn.f90
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
!
!
!
!
!
!
!
!
!
USE MODD_DATA_BEM_n, ONLY : DTB => DATA_BEM
USE MODD_DATA_TEB_GARDEN_n, ONLY : DTGD => DATA_TEB_GARDEN
USE MODD_DATA_TEB_GREENROOF_n, ONLY : DTGR => DATA_TEB_GREENROOF
USE MODD_TEB_GARDEN_PGD_n, ONLY : TGDP => TEB_GARDEN_PGD
USE MODD_TEB_GREENROOF_PGD_n, ONLY : TGRP => TEB_GREENROOF_PGD
USE MODD_TEB_GRID_n, ONLY : TG => TEB_GRID
USE MODD_TEB_IRRIG_n, ONLY : TIR => TEB_IRRIG
!
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
USE MODD_DUMMY_SURF_FIELDS_n, ONLY : DUMMY_SURF_FIELDS_t
USE MODD_FLAKE_GRID_n, ONLY : FLAKE_GRID_t
USE MODD_FLAKE_n, ONLY : FLAKE_t
USE MODD_FLAKE_SBL_n, ONLY : FLAKE_SBL_t
USE MODD_GR_BIOG_n, ONLY : GR_BIOG_t
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
USE MODD_ISBA_CANOPY_n, ONLY : ISBA_CANOPY_t
USE MODD_ISBA_GRID_n, ONLY : ISBA_GRID_t
USE MODD_ISBA_n, ONLY : ISBA_t
USE MODD_OCEAN_n, ONLY : OCEAN_t
USE MODD_SEAFLUX_GRID_n, ONLY : SEAFLUX_GRID_t
USE MODD_SEAFLUX_n, ONLY : SEAFLUX_t
USE MODD_SEAFLUX_SBL_n, ONLY : SEAFLUX_SBL_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SURF_ATM_SSO_n, ONLY : SURF_ATM_SSO_t
USE MODD_SV_n, ONLY : SV_t
USE MODD_TEB_CANOPY_n, ONLY : TEB_CANOPY_t
USE MODD_TEB_GARDEN_n, ONLY : TEB_GARDEN_t
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TEB_GARDEN_OPTIONS_t
USE MODD_TEB_GREENROOF_n, ONLY : TEB_GREENROOF_t
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TEB_GREENROOF_OPTIONS_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_OPTIONS_t
USE MODD_WATFLUX_GRID_n, ONLY : WATFLUX_GRID_t
USE MODD_WATFLUX_n, ONLY : WATFLUX_t
USE MODD_WATFLUX_SBL_n, ONLY : WATFLUX_SBL_t
!
USE MODD_SURF_CONF,       ONLY : CPROGNAME
USE MODD_SURF_PAR,        ONLY : NVERSION, NBUGFIX
USE MODD_IO_SURF_FA,      ONLY : LFANOCOMPACT
!
USE MODI_INIT_IO_SURF_n
USE MODI_WRITE_SURF
USE MODI_WRITE_PGD_SEA_n
USE MODI_WRITE_PGD_INLAND_WATER_n
USE MODI_WRITE_PGD_NATURE_n
USE MODI_WRITE_PGD_TOWN_n
USE MODI_END_IO_SURF_n
!
USE MODI_FLAG_UPDATE
!
USE MODI_WRITESURF_COVER_n
USE MODI_WRITESURF_SSO_n
USE MODI_WRITESURF_DUMMY_n
USE MODI_WRITESURF_SNAP_n
USE MODI_WRITESURF_CH_EMIS_n
USE MODI_WRITE_GRID
!
USE MODI_WRITE_ECOCLIMAP2_DATA
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
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
TYPE(DUMMY_SURF_FIELDS_t), INTENT(INOUT) :: DUU
TYPE(FLAKE_GRID_t), INTENT(INOUT) :: FG
TYPE(FLAKE_t), INTENT(INOUT) :: F
TYPE(FLAKE_SBL_t), INTENT(INOUT) :: FSB
TYPE(GR_BIOG_t), INTENT(INOUT) :: GB
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
TYPE(ISBA_CANOPY_t), INTENT(INOUT) :: ICP
TYPE(ISBA_GRID_t), INTENT(INOUT) :: IG
TYPE(ISBA_t), INTENT(INOUT) :: I
TYPE(OCEAN_t), INTENT(INOUT) :: O
TYPE(SEAFLUX_GRID_t), INTENT(INOUT) :: SG
TYPE(SEAFLUX_t), INTENT(INOUT) :: S
TYPE(SEAFLUX_SBL_t), INTENT(INOUT) :: SSB
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SURF_ATM_SSO_t), INTENT(INOUT) :: USS
TYPE(SV_t), INTENT(INOUT) :: SV
TYPE(TEB_CANOPY_t), INTENT(INOUT) :: TCP
TYPE(TEB_GARDEN_t), INTENT(INOUT) :: TGD
TYPE(TEB_GARDEN_OPTIONS_t), INTENT(INOUT) :: TGDO
TYPE(TEB_GREENROOF_t), INTENT(INOUT) :: TGR
TYPE(TEB_GREENROOF_OPTIONS_t), INTENT(INOUT) :: TGRO
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_VEG_OPTIONS_t), INTENT(INOUT) :: TVG
TYPE(WATFLUX_GRID_t), INTENT(INOUT) :: WG
TYPE(WATFLUX_t), INTENT(INOUT) :: W
TYPE(WATFLUX_SBL_t), INTENT(INOUT) :: WSB
!
 CHARACTER(LEN=6),    INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
 CHARACTER(LEN=3)   :: YWRITE
 CHARACTER(LEN=100) :: YCOMMENT
INTEGER            :: IRESP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITE_PGD_SURF_ATM_N',0,ZHOOK_HANDLE)
!
!*       0.     Initialize some options:
!               ------------------------
!
CPROGNAME = HPROGRAM
!
 CALL FLAG_UPDATE(DGI, DGU, &
                  .FALSE.,.TRUE.,.FALSE.,.FALSE.)
!
!*       1.     Configuration and cover fields:
!               ------------------------------
!
!
!         Initialisation for IO
!
 CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                     HPROGRAM,'FULL  ','SURF  ','WRITE')
!
YWRITE='PGD'
YCOMMENT='(-)'
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,'VERSION',NVERSION,IRESP,YCOMMENT)
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,'BUG    ',NBUGFIX ,IRESP,YCOMMENT)
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,'STORAGETYPE',YWRITE,IRESP,YCOMMENT)
!
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,'SEA   ',U%CSEA   ,IRESP,YCOMMENT)
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,'WATER ',U%CWATER ,IRESP,YCOMMENT)
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,'NATURE',U%CNATURE,IRESP,YCOMMENT)
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,'TOWN  ',U%CTOWN  ,IRESP,YCOMMENT)
!
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,'DIM_FULL  ',U%NDIM_FULL,  IRESP,HCOMMENT=YCOMMENT)
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,'DIM_SEA   ',U%NDIM_SEA,   IRESP,HCOMMENT=YCOMMENT)
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,'DIM_NATURE',U%NDIM_NATURE,IRESP,HCOMMENT=YCOMMENT)
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,'DIM_WATER ',U%NDIM_WATER, IRESP,HCOMMENT=YCOMMENT)
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,'DIM_TOWN  ',U%NDIM_TOWN,  IRESP,HCOMMENT=YCOMMENT)
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,'ECOCLIMAP ',U%LECOCLIMAP ,IRESP,YCOMMENT)
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,'WATER_TO_NAT',U%LWATER_TO_NATURE,IRESP,YCOMMENT)
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,'TOWN_TO_ROCK',U%LTOWN_TO_ROCK,IRESP,YCOMMENT)
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,'GARDEN',U%LGARDEN,IRESP,YCOMMENT)
IF (HPROGRAM.NE.'BINARY' .AND. HPROGRAM.NE.'TEXTE ') THEN
   CALL WRITE_ECOCLIMAP2_DATA(DGU, IOB, U, &
                              HPROGRAM)
ENDIF
!
 CALL WRITE_GRID(DGU, IOB, U, &
                 HPROGRAM,UG%CGRID,UG%XGRID_PAR,UG%XLAT,UG%XLON,UG%XMESH_SIZE,IRESP,USS%XZ0EFFJPDIR)
!
 CALL WRITESURF_COVER_n(DGU, IOB, &
                        U, &
                        HPROGRAM)
 CALL WRITESURF_SSO_n(DGU, IOB, U, &
                      USS, &
                      HPROGRAM)
 CALL WRITESURF_DUMMY_n(DGU, IOB, U, &
                        DUU, &
                        HPROGRAM)
!
YCOMMENT='CH_EMIS'
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,'CH_EMIS',CHU%LCH_EMIS,IRESP,HCOMMENT=YCOMMENT)
!
IF (CHU%LCH_EMIS) THEN
  YCOMMENT='CH_EMIS_OPT'
  CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,'CH_EMIS_OPT',CHU%CCH_EMIS,IRESP,HCOMMENT=YCOMMENT)
END IF
!
IF (CHU%LCH_EMIS) THEN
  IF (CHU%CCH_EMIS=='AGGR') THEN
    CALL WRITESURF_CH_EMIS_n(DGU, IOB, U, &
                             CHE, &
                             HPROGRAM)
  ELSE IF (CHU%CCH_EMIS=='SNAP') THEN
    CALL WRITESURF_SNAP_n(DGU, IOB, U, &
                          CHN, &
                          HPROGRAM)
  ENDIF
ENDIF
!
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
!
!
!*       2.     Sea
!               ---
!
IF (U%NDIM_SEA>0) CALL WRITE_PGD_SEA_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTCO, DTT, &
                                   DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, &
                                   DGUT, DGW, F, FSB, GB, IOB, ICP, I, O, SSB, UG, &
                                   SV, TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                                   DTS, SG, S, &
                                       U, &
                                       HPROGRAM)
!
!
!*       3.     Inland water
!               ------------
!
IF (U%NDIM_WATER>0) CALL WRITE_PGD_INLAND_WATER_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTCO, DTS, &
                                            DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, &
                                            DGT, DGUT, DGW, FSB, GB, IOB, ICP, I, O, S, SSB, &
                                            UG, SV, TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, WSB, &
                                            FG, F, WG, W, &
                                                  U, &
                                                  HPROGRAM)
!
!
!*       4.     Vegetation scheme
!               -----------------
!
IF (U%NDIM_NATURE>0) CALL WRITE_PGD_NATURE_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTCO, DTS, &
                                      DTT, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, &
                                      DGUT, DGW, F, FSB, GB, IOB, ICP, O, S, SSB, UG, &
                                      SV, TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                                      DTI, DTZ, IG, I, &
                                             U, &
                                             HPROGRAM)
!
!
!*       5.     Urban scheme
!               ------------
!
IF (U%NDIM_TOWN>0) CALL WRITE_PGD_TOWN_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTB, DTCO, &
                                    DTS, DTGD, DTGR, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, &
                                    DGS, DGSI, DGU, DGT, DGUT, DGW, F, FSB, GB, IOB, ICP, &
                                    I, O, S, SSB, UG, SV, TCP, TGD, TGDO, TGDP, TGR, &
                                    TGRO, TGRP, TG, TIR, T, TOP, TVG, W, WSB, &
                                         U, &
                                         HPROGRAM)
!
!
IF (LHOOK) CALL DR_HOOK('WRITE_PGD_SURF_ATM_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_PGD_SURF_ATM_n
