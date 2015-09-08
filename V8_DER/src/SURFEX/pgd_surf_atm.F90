!     ###########################################################
      SUBROUTINE PGD_SURF_ATM (BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTB, DTCO, &
                               DTI, DTS, DTGD, DTGR, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, &
                               DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, DUU, FG, F, FSB, &
                               GB, IOB, ICP, IG, I, O, SG, S, SSB, UG, U, &
                               USS, SV, TCP, TGD, TGDO, TGDP, TGR, TGRO, TGRP, TG, TIR, &
                               T, TOP, TVG, WG, W, WSB, &
                               HPROGRAM,HFILE,HFILETYPE,OZS)
!     ###########################################################
!!
!!    PURPOSE
!!    -------
!!   This program prepares the physiographic data fields.
!!
!!    METHOD
!!    ------
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
!!    V. Masson                   Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original     13/10/03
!!      A. Lemonsu      05/2009         Ajout de la clef LGARDEN pour TEB
!!      J. Escobar      11/2013         Add USE MODI_READ_NAM_PGD_CHEMISTRY
!!      B. Decharme     02/2014         Add LRM_RIVER
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
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
USE MODD_BEM_OPTION_n, ONLY : BEM_OPTIONS_t
USE MODD_BLD_DESCRIPTION_n, ONLY : BLD_DESC_t
USE MODD_CH_EMIS_FIELD_n, ONLY : CH_EMIS_FIELD_t
USE MODD_CH_ISBA_n, ONLY : CH_ISBA_t
USE MODD_CH_SEAFLUX_n, ONLY : CH_SEAFLUX_t
USE MODD_CH_SNAP_n, ONLY : CH_EMIS_SNAP_t
USE MODD_CH_SURF_n, ONLY : CH_SURF_t
USE MODD_CH_TEB_n, ONLY : CH_TEB_t
USE MODD_CH_WATFLUX_n, ONLY : CH_WATFLUX_t
USE MODD_DATA_BEM_n, ONLY : DATA_BEM_t
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_DATA_SEAFLUX_n, ONLY : DATA_SEAFLUX_t
USE MODD_DATA_TEB_GARDEN_n, ONLY : DATA_TEB_GARDEN_t
USE MODD_DATA_TEB_GREENROOF_n, ONLY : DATA_TEB_GREENROOF_t
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
USE MODD_TEB_GARDEN_PGD_n, ONLY : TEB_GARDEN_PGD_t
USE MODD_TEB_GREENROOF_n, ONLY : TEB_GREENROOF_t
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TEB_GREENROOF_OPTIONS_t
USE MODD_TEB_GREENROOF_PGD_n, ONLY : TEB_GREENROOF_PGD_t
USE MODD_TEB_GRID_n, ONLY : TEB_GRID_t
USE MODD_TEB_IRRIG_n, ONLY : TEB_IRRIG_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_OPTIONS_t
USE MODD_WATFLUX_GRID_n, ONLY : WATFLUX_GRID_t
USE MODD_WATFLUX_n, ONLY : WATFLUX_t
USE MODD_WATFLUX_SBL_n, ONLY : WATFLUX_SBL_t
!
USE MODD_SURF_CONF,       ONLY : CPROGNAME
USE MODD_PGD_GRID,        ONLY : LLATLONMASK
!
USE MODI_GET_LUOUT
USE MODI_READ_PGD_ARRANGE_COVER
USE MODI_READ_PGD_COVER_GARDEN
USE MODI_INI_DATA_COVER
USE MODI_READ_PGD_SCHEMES
USE MODI_READ_NAM_PGD_CHEMISTRY
USE MODI_READ_NAM_WRITE_COVER_TEX
USE MODI_WRITE_COVER_TEX_START
USE MODI_WRITE_COVER_TEX_COVER
USE MODI_LATLON_GRID
USE MODI_PUT_PGD_GRID
USE MODI_LATLONMASK
USE MODI_PGD_FRAC
USE MODI_PGD_COVER
USE MODI_PGD_OROGRAPHY
USE MODI_PGD_NATURE
USE MODI_PGD_TOWN
USE MODI_PGD_INLAND_WATER
USE MODI_PGD_SEA
USE MODI_PGD_DUMMY
USE MODI_PGD_CHEMISTRY
USE MODI_PGD_CHEMISTRY_SNAP
USE MODI_WRITE_COVER_TEX_END
USE MODI_INIT_READ_DATA_COVER
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
!
IMPLICIT NONE
!
!*    0.1    Declaration of dummy arguments
!            ------------------------------
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
TYPE(DATA_BEM_t), INTENT(INOUT) :: DTB
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTI
TYPE(DATA_SEAFLUX_t), INTENT(INOUT) :: DTS
TYPE(DATA_TEB_GARDEN_t), INTENT(INOUT) :: DTGD
TYPE(DATA_TEB_GREENROOF_t), INTENT(INOUT) :: DTGR
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
TYPE(TEB_GARDEN_PGD_t), INTENT(INOUT) :: TGDP
TYPE(TEB_GREENROOF_t), INTENT(INOUT) :: TGR
TYPE(TEB_GREENROOF_OPTIONS_t), INTENT(INOUT) :: TGRO
TYPE(TEB_GREENROOF_PGD_t), INTENT(INOUT) :: TGRP
TYPE(TEB_GRID_t), INTENT(INOUT) :: TG
TYPE(TEB_IRRIG_t), INTENT(INOUT) :: TIR
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_VEG_OPTIONS_t), INTENT(INOUT) :: TVG
TYPE(WATFLUX_GRID_t), INTENT(INOUT) :: WG
TYPE(WATFLUX_t), INTENT(INOUT) :: W
TYPE(WATFLUX_SBL_t), INTENT(INOUT) :: WSB
!
CHARACTER(LEN=6),     INTENT(IN)  :: HPROGRAM ! program calling
CHARACTER(LEN=28),    INTENT(IN)  :: HFILE    ! atmospheric file name
CHARACTER(LEN=6),     INTENT(IN)  :: HFILETYPE! atmospheric file type
LOGICAL,              INTENT(IN)  :: OZS      ! .true. if orography is imposed by atm. model
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
LOGICAL :: LRM_RIVER   !delete inland river coverage. Default is false
!
INTEGER :: ILUOUT ! logical unit of output listing file
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('PGD_SURF_ATM',0,ZHOOK_HANDLE)
!
LRM_RIVER = .FALSE.
!
CPROGNAME=HPROGRAM
!
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
!*    1.      Set default constant values 
!             ---------------------------
!
 CALL READ_PGD_ARRANGE_COVER(HPROGRAM,U%LWATER_TO_NATURE,U%LTOWN_TO_ROCK)
!
 CALL READ_PGD_COVER_GARDEN(HPROGRAM,U%LGARDEN)
!
 CALL INIT_READ_DATA_COVER(HPROGRAM)
!
 CALL INI_DATA_COVER(DTCO, U)
!
!*    1.2     surface schemes
 CALL READ_PGD_SCHEMES(HPROGRAM,U%CNATURE,U%CSEA,U%CTOWN,U%CWATER)
!
!*    1.3     prints all parameters in a Latex file
 CALL READ_NAM_WRITE_COVER_TEX(HPROGRAM)
!
 CALL WRITE_COVER_TEX_START(HPROGRAM)
 CALL WRITE_COVER_TEX_COVER
!-------------------------------------------------------------------------------
!
!*    2.      Grid
!             ----
!
ALLOCATE(UG%XLAT(U%NSIZE_FULL))
ALLOCATE(UG%XLON(U%NSIZE_FULL))
ALLOCATE(UG%XMESH_SIZE(U%NSIZE_FULL))
ALLOCATE(UG%XJPDIR(U%NSIZE_FULL))
 CALL LATLON_GRID(UG%CGRID,UG%NGRID_PAR,U%NSIZE_FULL,ILUOUT,UG%XGRID_PAR,UG%XLAT,UG%XLON,UG%XMESH_SIZE,UG%XJPDIR)
!
!
!*    2.3     Stores the grid in the module MODD_PGD_GRID
!
 CALL PUT_PGD_GRID(UG%CGRID,U%NSIZE_FULL,UG%NGRID_PAR,UG%XGRID_PAR)
!
!*    2.4     mask to limit the number of input data to read
 CALL LATLONMASK      (UG%CGRID,UG%NGRID_PAR,UG%XGRID_PAR,LLATLONMASK)
!
!-------------------------------------------------------------------------------
!
!*    3.      surface cover
!             -------------
!
 CALL PGD_FRAC(DTCO, UG, U, USS, &
               HPROGRAM,U%LECOCLIMAP)
IF (U%LECOCLIMAP) CALL PGD_COVER(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTS, DTT, &
                             DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, &
                             DGUT, DGW, F, FSB, GB, IOB, ICP, I, O, S, SSB, &
                             SV, TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                                 DTCO, UG, U, USS, &
                                 HPROGRAM,LRM_RIVER)
!
!-------------------------------------------------------------------------------
!
!*    4.      Orography
!             ---------
!
 CALL PGD_OROGRAPHY(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTCO, DTS, &
                                DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, &
                                DGT, DGUT, DGW, F, FSB, GB, IOB, ICP, I, O, S, &
                                SSB, UG, U, USS, SV, TCP, TGD, TGDO, TGR, TGRO, T, &
                                TOP, TVG, W, WSB, &
                    HPROGRAM,U%XSEA,U%XWATER,HFILE,HFILETYPE,OZS)
!
!_______________________________________________________________________________
!
!*    5.      Additionnal fields for nature scheme
!             ------------------------------------
!
IF (U%NDIM_NATURE>0) CALL PGD_NATURE(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTCO, DTI, &
                             DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, &
                             DGU, DGT, DGUT, DGW, F, FSB, GB, IOB, ICP, IG, I, &
                             O, S, SSB, UG, U, USS, SV, TCP, TGD, TGDO, TGR, &
                             TGRO, T, TOP, TVG, W, WSB, &
                                     HPROGRAM,U%LECOCLIMAP)  
!_______________________________________________________________________________
!
!*    6.      Additionnal fields for town scheme
!             ----------------------------------
!
IF (U%NDIM_TOWN>0) CALL PGD_TOWN(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTB, DTCO, &
                           DTI, DTS, DTGD, DTGR, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, &
                           DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, F, FSB, GB, IOB, &
                           ICP, I, O, S, SSB, UG, U, USS, SV, TCP, TGD, &
                           TGDO, TGDP, TGR, TGRO, TGRP, TG, TIR, T, TOP, TVG, W, &
                           WSB, &
                                 HPROGRAM,U%LECOCLIMAP,U%LGARDEN)  
!_______________________________________________________________________________
!
!*    7.      Additionnal fields for inland water scheme
!             ------------------------------------------
!
IF (U%NDIM_WATER>0) CALL PGD_INLAND_WATER(DTCO, FG, F, UG, U, USS, WG, W, &
                                          HPROGRAM,U%LECOCLIMAP,LRM_RIVER)   
!_______________________________________________________________________________
!
!*    8.      Additionnal fields for sea scheme
!             ---------------------------------
!
IF (U%NDIM_SEA>0) CALL PGD_SEA(DTCO, DTS, SG, S, UG, U, USS, &
                               HPROGRAM)  
!
!_______________________________________________________________________________
!
!*    9.      Dummy fields
!             ------------
!
 CALL PGD_DUMMY(DTCO, DUU, UG, U, USS, &
                HPROGRAM)
!_______________________________________________________________________________
!
!*   10.      Chemical Emission fields
!             ------------------------
!
 CALL READ_NAM_PGD_CHEMISTRY(HPROGRAM,CHU%CCH_EMIS)
IF (CHU%CCH_EMIS=='SNAP') THEN
  CALL PGD_CHEMISTRY_SNAP(CHN, DTCO, UG, U, USS, &
                          HPROGRAM,CHU%LCH_EMIS)
ELSE IF (CHU%CCH_EMIS=='AGGR') THEN
  CALL PGD_CHEMISTRY(CHE, DTCO, UG, U, USS, &
                     HPROGRAM,CHU%LCH_EMIS)
ENDIF
!_______________________________________________________________________________
!
!*   11.     Writing in cover latex file
!            ---------------------------
!
 CALL WRITE_COVER_TEX_END(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('PGD_SURF_ATM',1,ZHOOK_HANDLE)
!_______________________________________________________________________________
!
END SUBROUTINE PGD_SURF_ATM
