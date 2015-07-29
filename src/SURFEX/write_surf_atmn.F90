!     ####################################
      SUBROUTINE WRITE_SURF_ATM_n (B, BOP, BDD, CHE, CHF, CHI, CHS, CHN, CHU, CHT, CHW, &
                                   DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMF, DGMI, DGMTO, DGO, &
                                   DGS, DGSI, DGU, DGT, DGUT, DGW, DST, F, FSB, GB, IOB, &
                                   ICP, I, O, OR, S, SSB, SSCP, UG, U, USS, SV, &
                                   TCP, TGD, TGDO, TGDPE, TGR, TGRO, TGRPE, T, TOP, TPN, TVG, &
                                   W, WSB, &
                                   HPROGRAM,HWRITE,OLAND_USE)
!     ####################################
!
!!****  *WRITE_SURF_ATM_n* - routine to write surface variables 
!!                           in their respective files or in file
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
!!      Original    01/2003
!!      Modified    06/2007, P.LeMoigne: do not write pgd fields in
!!                                       historical files
!!      Modified    03/2009, B.Decharme: keys for arrange cover
!!      Modified    04/2009, B.Decharme: write precipitation forcing into the restart file for ARPEGE/ALADIN run
!       Modified    06/2009, B.Decharme: flag to desactivate writing of horizontal grid 
!       Modified    08/2009, B.Decharme: BUDGETC for all tiles
!       Modified    07/2011, B.Decharme: delete write pgd fields
!       Modified    07/2011, B.Decharme: land_use key for writing semi-prognostic variables
!       Modified    05/2012, B.Decharme: supress LPROVAR_TO_DIAG to write prognostic fields if user want
!       Modified    05/2013, B.Decharme: WRITESURF_PRECIP becomes WRITESURF_CPL_GCM
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
USE MODD_BEM_n, ONLY : BEM_t
USE MODD_BEM_OPTION_n, ONLY : BEM_OPTIONS_t
USE MODD_BLD_DESCRIPTION_n, ONLY : BLD_DESC_t
USE MODD_CH_EMIS_FIELD_n, ONLY : CH_EMIS_FIELD_t
USE MODD_CH_FLAKE_n, ONLY : CH_FLAKE_t
USE MODD_CH_ISBA_n, ONLY : CH_ISBA_t
USE MODD_CH_SEAFLUX_n, ONLY : CH_SEAFLUX_t
USE MODD_CH_SNAP_n, ONLY : CH_EMIS_SNAP_t
USE MODD_CH_SURF_n, ONLY : CH_SURF_t
USE MODD_CH_TEB_n, ONLY : CH_TEB_t
USE MODD_CH_WATFLUX_n, ONLY : CH_WATFLUX_t
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DATA_SEAFLUX_n, ONLY : DATA_SEAFLUX_t
USE MODD_DATA_TEB_n, ONLY : DATA_TEB_t
USE MODD_DATA_TSZ0_n, ONLY : DATA_TSZ0_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_FLAKE_n, ONLY : DIAG_FLAKE_t
USE MODD_DIAG_ISBA_n, ONLY : DIAG_ISBA_t
USE MODD_DIAG_MISC_FLAKE_n, ONLY : DIAG_MISC_FLAKE_t
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
USE MODD_ISBA_n, ONLY : ISBA_t
USE MODD_OCEAN_n, ONLY : OCEAN_t
USE MODD_OCEAN_REL_n, ONLY : OCEAN_REL_t
USE MODD_SEAFLUX_n, ONLY : SEAFLUX_t
USE MODD_SEAFLUX_SBL_n, ONLY : SEAFLUX_SBL_t
USE MODD_SSO_CANOPY_n, ONLY : SSO_CANOPY_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SURF_ATM_SSO_n, ONLY : SURF_ATM_SSO_t
USE MODD_SV_n, ONLY : SV_t
USE MODD_TEB_CANOPY_n, ONLY : TEB_CANOPY_t
USE MODD_TEB_GARDEN_n, ONLY : TEB_GARDEN_t
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TEB_GARDEN_OPTIONS_t
USE MODD_TEB_GARDEN_PGD_EVOL_n, ONLY : TEB_GARDEN_PGD_EVOL_t
USE MODD_TEB_GREENROOF_n, ONLY : TEB_GREENROOF_t
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TEB_GREENROOF_OPTIONS_t
USE MODD_TEB_GREENROOF_PGD_EVOL_n, ONLY : TEB_GREENROOF_PGD_EVOL_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_PANEL_n, ONLY : TEB_PANEL_t
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_OPTIONS_t
USE MODD_WATFLUX_n, ONLY : WATFLUX_t
USE MODD_WATFLUX_SBL_n, ONLY : WATFLUX_SBL_t
!
USE MODD_SURF_CONF,       ONLY : CPROGNAME
USE MODD_SURF_PAR,        ONLY : NVERSION, NBUGFIX
USE MODD_WRITE_SURF_ATM,  ONLY : LNOWRITE_CANOPY
!
USE MODI_INIT_IO_SURF_n
USE MODI_WRITE_SURF
USE MODI_WRITE_SEA_n
USE MODI_WRITE_INLAND_WATER_n
USE MODI_WRITE_NATURE_n
USE MODI_WRITE_TOWN_n
USE MODI_END_IO_SURF_n
USE MODI_WRITE_GRID
!
USE MODI_WRITESURF_ATM_CONF_n
USE MODI_WRITESURF_SSO_CANOPY_n
USE MODI_WRITESURF_CPL_GCM_n
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
TYPE(BEM_t), INTENT(INOUT) :: B
TYPE(BEM_OPTIONS_t), INTENT(INOUT) :: BOP
TYPE(BLD_DESC_t), INTENT(INOUT) :: BDD
TYPE(CH_EMIS_FIELD_t), INTENT(INOUT) :: CHE
TYPE(CH_FLAKE_t), INTENT(INOUT) :: CHF
TYPE(CH_ISBA_t), INTENT(INOUT) :: CHI
TYPE(CH_SEAFLUX_t), INTENT(INOUT) :: CHS
TYPE(CH_EMIS_SNAP_t), INTENT(INOUT) :: CHN
TYPE(CH_SURF_t), INTENT(INOUT) :: CHU
TYPE(CH_TEB_t), INTENT(INOUT) :: CHT
TYPE(CH_WATFLUX_t), INTENT(INOUT) :: CHW
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DATA_SEAFLUX_t), INTENT(INOUT) :: DTS
TYPE(DATA_TEB_t), INTENT(INOUT) :: DTT
TYPE(DATA_TSZ0_t), INTENT(INOUT) :: DTZ
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEI
TYPE(DIAG_FLAKE_t), INTENT(INOUT) :: DGF
TYPE(DIAG_ISBA_t), INTENT(INOUT) :: DGI
TYPE(DIAG_MISC_FLAKE_t), INTENT(INOUT) :: DGMF
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
TYPE(ISBA_t), INTENT(INOUT) :: I
TYPE(OCEAN_t), INTENT(INOUT) :: O
TYPE(OCEAN_REL_t), INTENT(INOUT) :: OR
TYPE(SEAFLUX_t), INTENT(INOUT) :: S
TYPE(SEAFLUX_SBL_t), INTENT(INOUT) :: SSB
TYPE(SSO_CANOPY_t), INTENT(INOUT) :: SSCP
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SURF_ATM_SSO_t), INTENT(INOUT) :: USS
TYPE(SV_t), INTENT(INOUT) :: SV
TYPE(TEB_CANOPY_t), INTENT(INOUT) :: TCP
TYPE(TEB_GARDEN_t), INTENT(INOUT) :: TGD
TYPE(TEB_GARDEN_OPTIONS_t), INTENT(INOUT) :: TGDO
TYPE(TEB_GARDEN_PGD_EVOL_t), INTENT(INOUT) :: TGDPE
TYPE(TEB_GREENROOF_t), INTENT(INOUT) :: TGR
TYPE(TEB_GREENROOF_OPTIONS_t), INTENT(INOUT) :: TGRO
TYPE(TEB_GREENROOF_PGD_EVOL_t), INTENT(INOUT) :: TGRPE
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_PANEL_t), INTENT(INOUT) :: TPN
TYPE(TEB_VEG_OPTIONS_t), INTENT(INOUT) :: TVG
TYPE(WATFLUX_t), INTENT(INOUT) :: W
TYPE(WATFLUX_SBL_t), INTENT(INOUT) :: WSB
!
 CHARACTER(LEN=6),    INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=3),    INTENT(IN)  :: HWRITE    ! 'PREP' : does not write SBL XUNDEF fields
!                                             ! 'ALL' : all fields are written
LOGICAL,             INTENT(IN)  :: OLAND_USE !
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
 CHARACTER(LEN=100) :: YCOMMENT
INTEGER            :: IRESP
LOGICAL            :: LSAVE_SELECT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITE_SURF_ATM_N',0,ZHOOK_HANDLE)
CPROGNAME = HPROGRAM
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
LSAVE_SELECT=DGU%LSELECT
DGU%LSELECT     =.FALSE.
!
YCOMMENT='(-)'
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,'VERSION',NVERSION,IRESP,YCOMMENT)
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,'BUG    ',NBUGFIX ,IRESP,YCOMMENT)
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,'STORAGETYPE',HWRITE,IRESP,YCOMMENT)
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,'DIM_FULL  ',U%NDIM_FULL,IRESP,HCOMMENT=YCOMMENT)
!
YCOMMENT='s'
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,'DTCUR',U%TTIME,IRESP,YCOMMENT)
!
DGU%LSELECT=LSAVE_SELECT
!
 CALL WRITE_GRID(DGU, IOB, U, &
                 HPROGRAM,UG%CGRID,UG%XGRID_PAR,UG%XLAT,UG%XLON,UG%XMESH_SIZE,IRESP)
!
 CALL WRITESURF_ATM_CONF_n(CHU, DGU, USS, &
                           HPROGRAM)
!
IF (HWRITE/='PRE') CALL WRITESURF_SSO_CANOPY_n(DGU, IOB, U, &
                                               SSCP, &
                                                         HPROGRAM,(USS%CROUGH=='BE04' .AND. .NOT. LNOWRITE_CANOPY))
!
 CALL WRITESURF_CPL_GCM_n(DGU, IOB, &
                          U, &
                          HPROGRAM)
!
YCOMMENT='flag for accumulated variables'
 CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,'BUDC',DGU%LSURF_BUDGETC,IRESP,HCOMMENT=YCOMMENT)
!
IF (DGU%LSURF_BUDGETC) THEN
   YCOMMENT='time of beginning of accumulation'
   CALL WRITE_SURF(DGU, IOB, U, &
                 HPROGRAM,'TBUDC',DGU%TIME_BUDGETC,IRESP,HCOMMENT=YCOMMENT)   
END IF
!  
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
!
!*       3.     Sea
!               ---
!
IF (U%NDIM_SEA>0) CALL WRITE_SEA_n(BOP, BDD, CHE, CHI, CHN, CHU, CHT, CHW, DTCO, DTS, DTT, &
                               DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGS, DGT, DGUT, DGW, F, &
                               FSB, GB, IOB, ICP, I, UG, SV, TCP, TGD, TGDO, TGR, &
                               TGRO, T, TOP, TVG, W, WSB, &
                                   CHS, DGO, DGSI, DGU, O, OR, S, SSB, U, &
                                   HPROGRAM,HWRITE)
!
!
!*       4.     Inland water
!               ------------
!
IF (U%NDIM_WATER>0) CALL WRITE_INLAND_WATER_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, DTCO, DTS, DTT, &
                                        DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGT, DGUT, &
                                        DGW, GB, IOB, ICP, I, O, S, SSB, UG, SV, TCP, &
                                        TGD, TGDO, TGR, TGRO, T, TOP, TVG, &
                                              CHF, CHW, DGMF, DGU, F, FSB, U, W, WSB, &
                                              HPROGRAM,HWRITE)
!
!
!*       5.     Vegetation scheme
!               -----------------
!
IF (U%NDIM_NATURE>0) CALL WRITE_NATURE_n(BOP, BDD, CHE, CHS, CHN, CHU, CHT, CHW, DTCO, DTS, DTT, &
                                  DTZ, DGF, DGMTO, DGO, DGS, DGSI, DGT, DGUT, DGW, F, FSB, &
                                  GB, IOB, O, S, SSB, UG, SV, TCP, TGD, TGDO, TGR, &
                                  TGRO, T, TOP, TVG, W, WSB, &
                                         CHI, DGEI, DGI, DGMI, DGU, DST, ICP, I, U, &
                                         HPROGRAM,HWRITE,OLAND_USE)
!
!
!*       6.     Urban scheme
!               ------------
!
IF (U%NDIM_TOWN>0) CALL WRITE_TOWN_n(BDD, CHE, CHI, CHS, CHN, CHU, CHW, DTCO, DTS, DTZ, DGEI, &
                                DGF, DGI, DGMI, DGO, DGS, DGSI, DGW, F, FSB, GB, IOB, &
                                ICP, I, O, S, SSB, UG, SV, W, WSB, &
                                     B, BOP, CHT, DTT, DGMTO, DGU, DGT, DGUT, U, &
                                     TCP, TGD, TGDO, TGDPE, TGR, TGRO, TGRPE, T, TOP, TPN, TVG, &
                                     HPROGRAM,HWRITE)
IF (LHOOK) CALL DR_HOOK('WRITE_SURF_ATM_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_SURF_ATM_n
