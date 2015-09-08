!     #################################################################################
SUBROUTINE WRITE_DIAG_SURF_ATM_n (BDD, DTCO, DTS, DTT, DTZ, FSB, IOB, ICP, O, SSB, TCP, &
                                   TGD, TGDO, TGR, TGRO, WSB, &
                                   B, BOP, CHE, CHF, CHI, CHS, CHN, CHU, CHT, CHW, &
                                  DGCT, DGEI, DGF, DGI, DGMF, DGMI, DGMT, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, DST, &
                                  F, GB, I, S, UG, U, SV, TGDPE, TGDP, T, TOP, TPN, TVG, W, &
                                  HPROGRAM,HWRITE)
!     #################################################################################
!
!!****  *WRITE_DIAG_SURF_ATM_n * - Chooses the surface schemes for diagnostics writing
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
!!------------------------------------------------------------------
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
!
USE MODD_BLD_DESCRIPTION_n, ONLY : BLD_DESC_t
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DATA_SEAFLUX_n, ONLY : DATA_SEAFLUX_t
USE MODD_DATA_TEB_n, ONLY : DATA_TEB_t
USE MODD_DATA_TSZ0_n, ONLY : DATA_TSZ0_t
USE MODD_FLAKE_SBL_n, ONLY : FLAKE_SBL_t
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
USE MODD_ISBA_CANOPY_n, ONLY : ISBA_CANOPY_t
USE MODD_OCEAN_n, ONLY : OCEAN_t
USE MODD_SEAFLUX_SBL_n, ONLY : SEAFLUX_SBL_t
USE MODD_TEB_CANOPY_n, ONLY : TEB_CANOPY_t
USE MODD_TEB_GARDEN_n, ONLY : TEB_GARDEN_t
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TEB_GARDEN_OPTIONS_t
USE MODD_TEB_GREENROOF_n, ONLY : TEB_GREENROOF_t
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TEB_GREENROOF_OPTIONS_t
USE MODD_WATFLUX_SBL_n, ONLY : WATFLUX_SBL_t
!
USE MODD_BEM_n, ONLY : BEM_t
USE MODD_BEM_OPTION_n, ONLY : BEM_OPTIONS_t
USE MODD_CH_EMIS_FIELD_n, ONLY : CH_EMIS_FIELD_t
USE MODD_CH_FLAKE_n, ONLY : CH_FLAKE_t
USE MODD_CH_ISBA_n, ONLY : CH_ISBA_t
USE MODD_CH_SEAFLUX_n, ONLY : CH_SEAFLUX_t
USE MODD_CH_SNAP_n, ONLY : CH_EMIS_SNAP_t
USE MODD_CH_SURF_n, ONLY : CH_SURF_t
USE MODD_CH_TEB_n, ONLY : CH_TEB_t
USE MODD_CH_WATFLUX_n, ONLY : CH_WATFLUX_t
USE MODD_DIAG_CUMUL_TEB_n, ONLY : DIAG_CUMUL_TEB_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_FLAKE_n, ONLY : DIAG_FLAKE_t
USE MODD_DIAG_ISBA_n, ONLY : DIAG_ISBA_t
USE MODD_DIAG_MISC_FLAKE_n, ONLY : DIAG_MISC_FLAKE_t
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t
USE MODD_DIAG_MISC_TEB_n, ONLY : DIAG_MISC_TEB_t
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
USE MODD_GR_BIOG_n, ONLY : GR_BIOG_t
USE MODD_ISBA_n, ONLY : ISBA_t
USE MODD_SEAFLUX_n, ONLY : SEAFLUX_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SV_n, ONLY : SV_t
USE MODD_TEB_GARDEN_PGD_EVOL_n, ONLY : TEB_GARDEN_PGD_EVOL_t
USE MODD_TEB_GARDEN_PGD_n, ONLY : TEB_GARDEN_PGD_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_PANEL_n, ONLY : TEB_PANEL_t
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_OPTIONS_t
USE MODD_WATFLUX_n, ONLY : WATFLUX_t
!
USE MODD_SURF_CONF,      ONLY : CPROGNAME
!
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
!
USE MODI_INIT_IO_SURF_n
USE MODI_WRITE_SURF
USE MODI_END_IO_SURF_n
!
USE MODI_WRITE_DIAG_NATURE_n 
USE MODI_WRITE_DIAG_SEA_n 
USE MODI_WRITE_DIAG_INLAND_WATER_n 
USE MODI_WRITE_DIAG_TOWN_n 
!
USE MODI_WRITE_DIAG_SEB_SURF_ATM_n
!
USE MODI_WRITE_DIAG_CH_AGGR_n
USE MODI_WRITE_DIAG_CH_SNAP_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
!
TYPE(BLD_DESC_t), INTENT(INOUT) :: BDD
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DATA_SEAFLUX_t), INTENT(INOUT) :: DTS
TYPE(DATA_TEB_t), INTENT(INOUT) :: DTT
TYPE(DATA_TSZ0_t), INTENT(INOUT) :: DTZ
TYPE(FLAKE_SBL_t), INTENT(INOUT) :: FSB
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
TYPE(ISBA_CANOPY_t), INTENT(INOUT) :: ICP
TYPE(OCEAN_t), INTENT(INOUT) :: O
TYPE(SEAFLUX_SBL_t), INTENT(INOUT) :: SSB
TYPE(TEB_CANOPY_t), INTENT(INOUT) :: TCP
TYPE(TEB_GARDEN_t), INTENT(INOUT) :: TGD
TYPE(TEB_GARDEN_OPTIONS_t), INTENT(INOUT) :: TGDO
TYPE(TEB_GREENROOF_t), INTENT(INOUT) :: TGR
TYPE(TEB_GREENROOF_OPTIONS_t), INTENT(INOUT) :: TGRO
TYPE(WATFLUX_SBL_t), INTENT(INOUT) :: WSB
!
TYPE(BEM_t), INTENT(INOUT) :: B
TYPE(BEM_OPTIONS_t), INTENT(INOUT) :: BOP
TYPE(CH_EMIS_FIELD_t), INTENT(INOUT) :: CHE
TYPE(CH_FLAKE_t), INTENT(INOUT) :: CHF
TYPE(CH_ISBA_t), INTENT(INOUT) :: CHI
TYPE(CH_SEAFLUX_t), INTENT(INOUT) :: CHS
TYPE(CH_EMIS_SNAP_t), INTENT(INOUT) :: CHN
TYPE(CH_SURF_t), INTENT(INOUT) :: CHU
TYPE(CH_TEB_t), INTENT(INOUT) :: CHT
TYPE(CH_WATFLUX_t), INTENT(INOUT) :: CHW
TYPE(DIAG_CUMUL_TEB_t), INTENT(INOUT) :: DGCT
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEI
TYPE(DIAG_FLAKE_t), INTENT(INOUT) :: DGF
TYPE(DIAG_ISBA_t), INTENT(INOUT) :: DGI
TYPE(DIAG_MISC_FLAKE_t), INTENT(INOUT) :: DGMF
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DGMI
TYPE(DIAG_MISC_TEB_t), INTENT(INOUT) :: DGMT
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
TYPE(GR_BIOG_t), INTENT(INOUT) :: GB
TYPE(ISBA_t), INTENT(INOUT) :: I
TYPE(SEAFLUX_t), INTENT(INOUT) :: S
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SV_t), INTENT(INOUT) :: SV
TYPE(TEB_GARDEN_PGD_EVOL_t), INTENT(INOUT) :: TGDPE
TYPE(TEB_GARDEN_PGD_t), INTENT(INOUT) :: TGDP
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_PANEL_t), INTENT(INOUT) :: TPN
TYPE(TEB_VEG_OPTIONS_t), INTENT(INOUT) :: TVG
TYPE(WATFLUX_t), INTENT(INOUT) :: W
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=3),   INTENT(IN)  :: HWRITE    ! 'PGD' : only physiographic fields are written
!                                            ! 'ALL' : all fields are written
!
!
!*      0.2    declarations of local variables
!
 CHARACTER(LEN=100) :: YCOMMENT
INTEGER            :: IRESP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_SURF_ATM_N',0,ZHOOK_HANDLE)
CPROGNAME = HPROGRAM
!
IF (U%NDIM_SEA    >0) CALL WRITE_DIAG_SEA_n(BOP, BDD, CHE, CHI, CHN, CHU, CHT, CHW, DTCO, DTS, DTT, &
                              DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGT, DGUT, DGW, F, FSB, &
                              GB, IOB, ICP, I, O, SSB, UG, SV, TCP, TGD, TGDO, &
                              TGR, TGRO, T, TOP, TVG, W, WSB, &
                                            CHS, DGO, DGS, DGSI, DGU, S, U, &
                                                     HPROGRAM,HWRITE)
IF (U%NDIM_WATER  >0) CALL WRITE_DIAG_INLAND_WATER_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, DTCO, DTS, DTT, &
                                       DTZ, DGEI, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGT, DGUT, FSB, &
                                       GB, IOB, ICP, I, O, S, SSB, UG, SV, TCP, TGD, &
                                       TGDO, TGR, TGRO, T, TOP, TVG, WSB, &
                                                     CHF, CHW, DGF, DGMF, DGU, DGW, F, U, W, &
                                                     HPROGRAM,HWRITE)
IF (U%NDIM_NATURE >0) CALL WRITE_DIAG_NATURE_n(BOP, BDD, CHE, CHS, CHN, CHU, CHT, CHW, DTCO, DTS, DTT, &
                                 DTZ, DGF, DGMTO, DGO, DGS, DGSI, DGT, DGUT, DGW, F, FSB, &
                                 IOB, ICP, O, S, SSB, UG, SV, TCP, TGD, TGDO, TGR, &
                                 TGRO, T, TOP, TVG, W, WSB, &
                                               CHI, DGEI, DGI, DGMI, DGU, DST, GB, I, U, &
                                                     HPROGRAM,HWRITE)
IF (U%NDIM_TOWN   >0) CALL WRITE_DIAG_TOWN_n(BDD, CHE, CHI, CHS, CHN, CHU, CHW, DTCO, DTS, DTT, DTZ, &
                               DGEI, DGF, DGI, DGMI, DGO, DGS, DGSI, DGU, DGW, F, FSB, &
                               GB, IOB, ICP, I, O, S, SSB, UG, SV, TCP, TGD, &
                               TGDO, TGR, TGRO, W, WSB, &
                                             B, BOP, CHT, DGCT, DGMT, DGMTO, DGT, DGUT, U, TGDPE, TGDP, T, TOP, TPN, TVG, &
                                                     HPROGRAM,HWRITE)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! Writing
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
!
IF (DGU%XDIAG_TSTEP==XUNDEF .OR. ABS(NINT(U%TTIME%TIME/DGU%XDIAG_TSTEP)*DGU%XDIAG_TSTEP-U%TTIME%TIME)<1.E-3 ) THEN
  !
  IF (DGU%LFRAC) THEN
    CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                        HPROGRAM,'FULL  ','SURF  ','WRITE')
    YCOMMENT = '(-)'
    CALL WRITE_SURF(DGU, IOB, U, &
                    HPROGRAM,'FRAC_SEA   ',U%XSEA,   IRESP,HCOMMENT=YCOMMENT)
    CALL WRITE_SURF(DGU, IOB, U, &
                    HPROGRAM,'FRAC_NATURE',U%XNATURE,IRESP,HCOMMENT=YCOMMENT)
    CALL WRITE_SURF(DGU, IOB, U, &
                    HPROGRAM,'FRAC_WATER ',U%XWATER, IRESP,HCOMMENT=YCOMMENT)
    CALL WRITE_SURF(DGU, IOB, U, &
                    HPROGRAM,'FRAC_TOWN  ',U%XTOWN,  IRESP,HCOMMENT=YCOMMENT)
    CALL END_IO_SURF_n(HPROGRAM)
  END IF
  !
  IF (HWRITE/='PGD'.AND.DGU%LDIAG_GRID) CALL WRITE_DIAG_SEB_SURF_ATM_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTCO, DTS, &
                                             DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGT, &
                                             DGUT, DGW, F, FSB, GB, IOB, ICP, I, O, S, SSB, &
                                             U, SV, TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, &
                                             WSB, &
                                                                       DGU, UG, &
                                                                       HPROGRAM)
  !
  IF (CHU%LCH_EMIS .AND. SV%NBEQ>0 .AND. CHU%LCH_SURF_EMIS) THEN
    IF (CHU%CCH_EMIS=='AGGR') THEN 
      CALL WRITE_DIAG_CH_AGGR_n(BOP, BDD, CHI, CHS, CHN, CHU, CHT, CHW, DTCO, DTS, DTT, &
                                        DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, &
                                        DGUT, DGW, F, FSB, GB, IOB, ICP, I, O, S, SSB, &
                                        UG, U, SV, TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, &
                                        W, WSB, &
                                CHE, &
                                HPROGRAM)
    ELSE IF (CHU%CCH_EMIS=='SNAP') THEN
      CALL WRITE_DIAG_CH_SNAP_n(BOP, BDD, CHE, CHI, CHS, CHU, CHT, CHW, DTCO, DTS, DTT, &
                                        DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, &
                                        DGUT, DGW, F, FSB, GB, IOB, ICP, I, O, S, SSB, &
                                        UG, U, SV, TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, &
                                        W, WSB, &
                                CHN, &
                                HPROGRAM)
    END IF
  END IF
  !  
END IF
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_SURF_ATM_N',1,ZHOOK_HANDLE)
!
!--------------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_DIAG_SURF_ATM_n
