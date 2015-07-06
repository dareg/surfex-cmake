!     ####################################
      SUBROUTINE WRITE_SURF_ATM_n(HPROGRAM,HWRITE,OLAND_USE)
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
USE MODD_BLD_DESCRIPTION_n, ONLY : BDD => BLD_DESC
USE MODD_CH_EMIS_FIELD_n, ONLY : CHE => CH_EMIS_FIELD
USE MODD_CH_SNAP_n, ONLY : CHN => CH_EMIS_SNAP
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
USE MODD_DATA_SEAFLUX_n, ONLY : DTS => DATA_SEAFLUX
USE MODD_DATA_TSZ0_n, ONLY : DTZ => DATA_TSZ0
USE MODD_DIAG_FLAKE_n, ONLY : DGF => DIAG_FLAKE
USE MODD_DIAG_SEAFLUX_n, ONLY : DGS => DIAG_SEAFLUX
USE MODD_DIAG_WATFLUX_n, ONLY : DGW => DIAG_WATFLUX
USE MODD_GR_BIOG_n, ONLY : GB => GR_BIOG
USE MODD_SV_n, ONLY : SV => SV
!
USE MODD_BEM_n, ONLY : B => BEM
USE MODD_BEM_OPTION_n, ONLY : BOP => BEM_OPTIONS
USE MODD_CH_TEB_n, ONLY : CHT => CH_TEB
USE MODD_DATA_TEB_n, ONLY : DTT => DATA_TEB
USE MODD_DIAG_MISC_TEB_OPTION_n, ONLY : DGMTO => DIAG_MISC_TEB_OPTIONS
USE MODD_DIAG_TEB_n, ONLY : DGT => DIAG_TEB
USE MODD_DIAG_UTCI_TEB_n, ONLY : DGUT => DIAG_UTCI_TEB
USE MODD_TEB_CANOPY_n, ONLY : TCP => TEB_CANOPY
USE MODD_TEB_GARDEN_n, ONLY : TGD => TEB_GARDEN
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TGDO => TEB_GARDEN_OPTIONS
USE MODD_TEB_GARDEN_PGD_EVOL_n, ONLY : TGDPE => TEB_GARDEN_PGD_EVOL
USE MODD_TEB_GREENROOF_n, ONLY : TGR => TEB_GREENROOF
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TGRO => TEB_GREENROOF_OPTIONS
USE MODD_TEB_GREENROOF_PGD_EVOL_n, ONLY : TGRPE => TEB_GREENROOF_PGD_EVOL
USE MODD_TEB_n, ONLY : T => TEB
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS
USE MODD_TEB_PANEL_n, ONLY : TPN => TEB_PANEL
USE MODD_TEB_VEG_n, ONLY : TVG => TEB_VEG_OPTIONS
!
USE MODD_IO_BUFF_n, ONLY : IOB => IO_BUFF
!
USE MODD_CH_SEAFLUX_n, ONLY : CHS => CH_SEAFLUX
USE MODD_DIAG_OCEAN_n, ONLY : DGO => DIAG_OCEAN
USE MODD_DIAG_SEAICE_n, ONLY : DGSI => DIAG_SEAICE
USE MODD_OCEAN_n, ONLY : O => OCEAN
USE MODD_OCEAN_REL_n, ONLY : OR => OCEAN_REL
USE MODD_SEAFLUX_n, ONLY : S => SEAFLUX
USE MODD_SEAFLUX_SBL_n, ONLY : SSB => SEAFLUX_SBL
!
USE MODD_CH_ISBA_n, ONLY : CHI => CH_ISBA
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DGEI => DIAG_EVAP_ISBA
USE MODD_DIAG_ISBA_n, ONLY : DGI => DIAG_ISBA
USE MODD_DIAG_MISC_ISBA_n, ONLY : DGMI => DIAG_MISC_ISBA
USE MODD_DST_n, ONLY : DST => DST
USE MODD_ISBA_CANOPY_n, ONLY : ICP => ISBA_CANOPY
USE MODD_ISBA_n, ONLY : I => ISBA
!
USE MODD_CH_FLAKE_n, ONLY : CHF => CH_FLAKE
USE MODD_CH_WATFLUX_n, ONLY : CHW => CH_WATFLUX
USE MODD_DIAG_MISC_FLAKE_n, ONLY : DGMF => DIAG_MISC_FLAKE
USE MODD_FLAKE_n, ONLY : F => FLAKE
USE MODD_FLAKE_SBL_n, ONLY : FSB => FLAKE_SBL
USE MODD_WATFLUX_n, ONLY : W => WATFLUX
USE MODD_WATFLUX_SBL_n, ONLY : WSB => WATFLUX_SBL
!
USE MODD_CH_SURF_n, ONLY : CHU => CH_SURF
!
USE MODD_SSO_CANOPY_n, ONLY : SSCP => SSO_CANOPY
!
USE MODD_SURF_CONF,       ONLY : CPROGNAME
USE MODD_SURF_PAR,        ONLY : NVERSION, NBUGFIX
USE MODD_WRITE_SURF_ATM,  ONLY : LNOWRITE_CANOPY
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
USE MODD_SURF_ATM_SSO_n, ONLY : USS => SURF_ATM_SSO
USE MODD_DIAG_SURF_ATM_n, ONLY : DGU => DIAG_SURF_ATM  
USE MODD_SURF_ATM_GRID_n, ONLY : UG => SURF_ATM_GRID
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
IF (U%NDIM_SEA>0) CALL WRITE_SEA_n(CHS, DGO, DGSI, DGU, O, OR, S, SSB, U, &
                                   HPROGRAM,HWRITE)
!
!
!*       4.     Inland water
!               ------------
!
IF (U%NDIM_WATER>0) CALL WRITE_INLAND_WATER_n(CHF, CHW, DGMF, DGU, F, FSB, U, W, WSB, &
                                              HPROGRAM,HWRITE)
!
!
!*       5.     Vegetation scheme
!               -----------------
!
IF (U%NDIM_NATURE>0) CALL WRITE_NATURE_n(CHI, DGEI, DGI, DGMI, DGU, DST, ICP, I, U, &
                                         HPROGRAM,HWRITE,OLAND_USE)
!
!
!*       6.     Urban scheme
!               ------------
!
IF (U%NDIM_TOWN>0) CALL WRITE_TOWN_n(B, BOP, CHT, DTT, DGMTO, DGU, DGT, DGUT, U, &
                                     TCP, TGD, TGDO, TGDPE, TGR, TGRO, TGRPE, T, TOP, TPN, TVG, &
                                     HPROGRAM,HWRITE)
IF (LHOOK) CALL DR_HOOK('WRITE_SURF_ATM_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_SURF_ATM_n
