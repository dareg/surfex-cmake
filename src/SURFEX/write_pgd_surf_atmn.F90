!     ####################################
      SUBROUTINE WRITE_PGD_SURF_ATM_n(HPROGRAM)
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
USE MODD_BEM_OPTION_n, ONLY : BOP => BEM_OPTIONS
USE MODD_BLD_DESCRIPTION_n, ONLY : BDD => BLD_DESC
USE MODD_CH_ISBA_n, ONLY : CHI => CH_ISBA
USE MODD_CH_SEAFLUX_n, ONLY : CHS => CH_SEAFLUX
USE MODD_CH_TEB_n, ONLY : CHT => CH_TEB
USE MODD_CH_WATFLUX_n, ONLY : CHW => CH_WATFLUX
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
USE MODD_DATA_TEB_n, ONLY : DTT => DATA_TEB
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DGEI => DIAG_EVAP_ISBA
USE MODD_DIAG_FLAKE_n, ONLY : DGF => DIAG_FLAKE
USE MODD_DIAG_MISC_ISBA_n, ONLY : DGMI => DIAG_MISC_ISBA
USE MODD_DIAG_MISC_TEB_OPTION_n, ONLY : DGMTO => DIAG_MISC_TEB_OPTIONS
USE MODD_DIAG_OCEAN_n, ONLY : DGO => DIAG_OCEAN
USE MODD_DIAG_SEAFLUX_n, ONLY : DGS => DIAG_SEAFLUX
USE MODD_DIAG_SEAICE_n, ONLY : DGSI => DIAG_SEAICE
USE MODD_DIAG_TEB_n, ONLY : DGT => DIAG_TEB
USE MODD_DIAG_UTCI_TEB_n, ONLY : DGUT => DIAG_UTCI_TEB
USE MODD_DIAG_WATFLUX_n, ONLY : DGW => DIAG_WATFLUX
USE MODD_FLAKE_SBL_n, ONLY : FSB => FLAKE_SBL
USE MODD_GR_BIOG_n, ONLY : GB => GR_BIOG
USE MODD_ISBA_CANOPY_n, ONLY : ICP => ISBA_CANOPY
USE MODD_OCEAN_n, ONLY : O => OCEAN
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
USE MODD_WATFLUX_SBL_n, ONLY : WSB => WATFLUX_SBL
!
USE MODD_IO_BUFF_n, ONLY : IOB => IO_BUFF
!
USE MODD_DATA_SEAFLUX_n, ONLY : DTS => DATA_SEAFLUX
USE MODD_SEAFLUX_GRID_n, ONLY : SG => SEAFLUX_GRID
USE MODD_SEAFLUX_n, ONLY : S => SEAFLUX
!
USE MODD_DATA_ISBA_n, ONLY : DTI => DATA_ISBA
USE MODD_DATA_TSZ0_n, ONLY : DTZ => DATA_TSZ0
USE MODD_ISBA_GRID_n, ONLY : IG => ISBA_GRID
USE MODD_ISBA_n, ONLY : I => ISBA
!
USE MODD_FLAKE_GRID_n, ONLY : FG => FLAKE_GRID
USE MODD_FLAKE_n, ONLY : F => FLAKE
USE MODD_WATFLUX_GRID_n, ONLY : WG => WATFLUX_GRID
USE MODD_WATFLUX_n, ONLY : W => WATFLUX
!
USE MODD_CH_SNAP_n, ONLY : CHN => CH_EMIS_SNAP
!
USE MODD_DUMMY_SURF_FIELDS_n, ONLY : DUU => DUMMY_SURF_FIELDS
!
USE MODD_CH_EMIS_FIELD_n, ONLY : CHE => CH_EMIS_FIELD
!
USE MODD_DIAG_ISBA_n, ONLY : DGI => DIAG_ISBA
USE MODD_DIAG_SURF_ATM_n, ONLY : DGU => DIAG_SURF_ATM
!
USE MODD_SURF_CONF,       ONLY : CPROGNAME
USE MODD_SURF_PAR,        ONLY : NVERSION, NBUGFIX
USE MODD_SURF_ATM_GRID_n, ONLY : UG => SURF_ATM_GRID
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
USE MODD_SURF_ATM_SSO_n, ONLY : USS => SURF_ATM_SSO
USE MODD_CH_SURF_n, ONLY : CHU => CH_SURF
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
IF (U%NDIM_SEA>0) CALL WRITE_PGD_SEA_n(DTS, SG, S, &
                                       U, &
                                       HPROGRAM)
!
!
!*       3.     Inland water
!               ------------
!
IF (U%NDIM_WATER>0) CALL WRITE_PGD_INLAND_WATER_n(FG, F, WG, W, &
                                                  U, &
                                                  HPROGRAM)
!
!
!*       4.     Vegetation scheme
!               -----------------
!
IF (U%NDIM_NATURE>0) CALL WRITE_PGD_NATURE_n(DTI, DTZ, IG, I, &
                                             U, &
                                             HPROGRAM)
!
!
!*       5.     Urban scheme
!               ------------
!
IF (U%NDIM_TOWN>0) CALL WRITE_PGD_TOWN_n(U, &
                                         HPROGRAM)
!
!
IF (LHOOK) CALL DR_HOOK('WRITE_PGD_SURF_ATM_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_PGD_SURF_ATM_n
