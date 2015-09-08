!     #########
      SUBROUTINE WRITE_DIAG_PGD_GRDN_n (BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, DTCO, DTS, &
                                         DTT, DTZ, DGEI, DGF, DGI, DGMI, DGO, DGS, DGSI, DGU, DGT, &
                                         DGUT, DGW, F, FSB, GB, IOB, ICP, I, O, S, SSB, &
                                         UG, U, SV, TCP, TGD, TGDO, TGR, TGRO, T, TOP, W, &
                                         WSB, &
                                         DGMTO, TGDPE, TGDP, TVG, &
                                        HPROGRAM)
!     #########################################
!
!!****  *WRITE_DIAG_PGD_TEB_GARDEN_n* - writes the ISBA physiographic diagnostic fields
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
!!      Original    01/2004 
!!      Modified    10/2004 by P. Le Moigne: add XZ0REL, XVEGTYPE_PATCH
!!      Modified    11/2005 by P. Le Moigne: limit length of VEGTYPE_PATCH field names
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
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
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DATA_SEAFLUX_n, ONLY : DATA_SEAFLUX_t
USE MODD_DATA_TEB_n, ONLY : DATA_TEB_t
USE MODD_DATA_TSZ0_n, ONLY : DATA_TSZ0_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_FLAKE_n, ONLY : DIAG_FLAKE_t
USE MODD_DIAG_ISBA_n, ONLY : DIAG_ISBA_t
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t
USE MODD_DIAG_OCEAN_n, ONLY : DIAG_OCEAN_t
USE MODD_DIAG_SEAFLUX_n, ONLY : DIAG_SEAFLUX_t
USE MODD_DIAG_SEAICE_n, ONLY : DIAG_SEAICE_t
USE MODD_DIAG_SURF_ATM_n, ONLY : DIAG_SURF_ATM_t
USE MODD_DIAG_TEB_n, ONLY : DIAG_TEB_t
USE MODD_DIAG_UTCI_TEB_n, ONLY : DIAG_UTCI_TEB_t
USE MODD_DIAG_WATFLUX_n, ONLY : DIAG_WATFLUX_t
USE MODD_FLAKE_n, ONLY : FLAKE_t
USE MODD_FLAKE_SBL_n, ONLY : FLAKE_SBL_t
USE MODD_GR_BIOG_n, ONLY : GR_BIOG_t
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
USE MODD_ISBA_CANOPY_n, ONLY : ISBA_CANOPY_t
USE MODD_ISBA_n, ONLY : ISBA_t
USE MODD_OCEAN_n, ONLY : OCEAN_t
USE MODD_SEAFLUX_n, ONLY : SEAFLUX_t
USE MODD_SEAFLUX_SBL_n, ONLY : SEAFLUX_SBL_t
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
USE MODD_WATFLUX_n, ONLY : WATFLUX_t
USE MODD_WATFLUX_SBL_n, ONLY : WATFLUX_SBL_t
!
USE MODD_DIAG_MISC_TEB_OPTION_n, ONLY : DIAG_MISC_TEB_OPTIONS_t
USE MODD_TEB_GARDEN_PGD_EVOL_n, ONLY : TEB_GARDEN_PGD_EVOL_t
USE MODD_TEB_GARDEN_PGD_n, ONLY : TEB_GARDEN_PGD_t
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_OPTIONS_t
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
!
!
USE MODD_IO_SURF_FA, ONLY : LFANOCOMPACT, LPREP
!
USE MODI_INIT_IO_SURF_n
USE MODI_WRITE_SURF
USE MODI_END_IO_SURF_n
!
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
TYPE(DATA_SEAFLUX_t), INTENT(INOUT) :: DTS
TYPE(DATA_TEB_t), INTENT(INOUT) :: DTT
TYPE(DATA_TSZ0_t), INTENT(INOUT) :: DTZ
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEI
TYPE(DIAG_FLAKE_t), INTENT(INOUT) :: DGF
TYPE(DIAG_ISBA_t), INTENT(INOUT) :: DGI
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DGMI
TYPE(DIAG_OCEAN_t), INTENT(INOUT) :: DGO
TYPE(DIAG_SEAFLUX_t), INTENT(INOUT) :: DGS
TYPE(DIAG_SEAICE_t), INTENT(INOUT) :: DGSI
TYPE(DIAG_SURF_ATM_t), INTENT(INOUT) :: DGU
TYPE(DIAG_TEB_t), INTENT(INOUT) :: DGT
TYPE(DIAG_UTCI_TEB_t), INTENT(INOUT) :: DGUT
TYPE(DIAG_WATFLUX_t), INTENT(INOUT) :: DGW
TYPE(FLAKE_t), INTENT(INOUT) :: F
TYPE(FLAKE_SBL_t), INTENT(INOUT) :: FSB
TYPE(GR_BIOG_t), INTENT(INOUT) :: GB
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
TYPE(ISBA_CANOPY_t), INTENT(INOUT) :: ICP
TYPE(ISBA_t), INTENT(INOUT) :: I
TYPE(OCEAN_t), INTENT(INOUT) :: O
TYPE(SEAFLUX_t), INTENT(INOUT) :: S
TYPE(SEAFLUX_SBL_t), INTENT(INOUT) :: SSB
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
TYPE(WATFLUX_t), INTENT(INOUT) :: W
TYPE(WATFLUX_SBL_t), INTENT(INOUT) :: WSB
!
TYPE(DIAG_MISC_TEB_OPTIONS_t), INTENT(INOUT) :: DGMTO
TYPE(TEB_GARDEN_PGD_EVOL_t), INTENT(INOUT) :: TGDPE
TYPE(TEB_GARDEN_PGD_t), INTENT(INOUT) :: TGDP
TYPE(TEB_VEG_OPTIONS_t), INTENT(INOUT) :: TVG
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
 CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
 CHARACTER(LEN=100):: YCOMMENT       ! Comment string
 CHARACTER(LEN=2)  :: YLVLV, YPAS
!
INTEGER           :: JL, JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!         Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_PGD_GRDN_N',0,ZHOOK_HANDLE)
 CALL INIT_IO_SURF_n(BOP, BDD, CHE, CHI, CHS, CHN, CHU, CHT, CHW, &
                      DTCO, DTS, DTT, DTZ, DGEI, DGF, DGI, DGMI, DGMTO, DGO, DGS, DGSI, DGU, DGT, DGUT, DGW, &
                      F, FSB, GB, IOB, ICP, I, O, S, SSB, UG, U, SV, &
                      TCP, TGD, TGDO, TGR, TGRO, T, TOP, TVG, W, WSB, &
                     HPROGRAM,'TOWN  ','TEB   ','WRITE')
!
!* Leaf Area Index
!
IF (TVG%CPHOTO=='NON' .OR. TVG%CPHOTO=='AGS' .OR. TVG%CPHOTO=='AST') THEN
  !
  YRECFM='GD_LAI'
  YCOMMENT='leaf area index (-)'
  !
  CALL WRITE_SURF(DGU, IOB, U, &
                  HPROGRAM,YRECFM,TGDPE%CUR%XLAI(:),IRESP,HCOMMENT=YCOMMENT)
  !
ENDIF
!
!-------------------------------------------------------------------------------
!
!* Vegetation fraction
!
YRECFM='GD_VEG'
YCOMMENT='vegetation fraction (-)'
!
 CALL WRITE_SURF(DGU, IOB, U, &
                  HPROGRAM,YRECFM,TGDPE%CUR%XVEG(:),IRESP,HCOMMENT=YCOMMENT)
!
!* Surface roughness length (without snow)
!
YRECFM='GD_Z0VEG'
YCOMMENT='surface roughness length (without snow) (M)'
!
 CALL WRITE_SURF(DGU, IOB, U, &
                  HPROGRAM,YRECFM,TGDPE%CUR%XZ0(:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
!* Soil depth for each patch
!
DO JL=1,SIZE(TGDP%XDG,2)
  WRITE(YRECFM,FMT='(A5,I1)') 'GD_DG',JL
  YCOMMENT='soil depth'//' (M)'
  CALL WRITE_SURF(DGU, IOB, U, &
                  HPROGRAM,YRECFM,TGDP%XDG(:,JL),IRESP,HCOMMENT=YCOMMENT)
END DO
!
!-------------------------------------------------------------------------------
! For Earth System Model
IF(LFANOCOMPACT.AND..NOT.LPREP)THEN
   CALL END_IO_SURF_n(HPROGRAM)
   IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_PGD_GRDN_N',1,ZHOOK_HANDLE)
   RETURN
ENDIF
!
!-------------------------------------------------------------------------------
!
!* Runoff soil ice depth for each patch
!
IF(TVG%CHORT=='SGH')THEN
  YRECFM='GD_DICE'
  YCOMMENT='soil ice depth for runoff (m)'
  CALL WRITE_SURF(DGU, IOB, U, &
                  HPROGRAM,YRECFM,TGDP%XD_ICE(:),IRESP,HCOMMENT=YCOMMENT)
ENDIF
!
!-------------------------------------------------------------------------------
!
!* Fraction of each vegetation type for each patch
!
DO JL=1,SIZE(TGDP%XVEGTYPE,2)
  WRITE(YPAS,'(I2)') JL 
  YLVLV=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
  WRITE(YRECFM,FMT='(A12)') 'GD_VEGTY_P'//YLVLV
  YCOMMENT='fraction of each vegetation type '//' (-)'
  CALL WRITE_SURF(DGU, IOB, U, &
                  HPROGRAM,YRECFM,TGDP%XVEGTYPE(:,JL),IRESP,HCOMMENT=YCOMMENT)
END DO
!-------------------------------------------------------------------------------
!
!* other surface parameters
!
YRECFM='GD_RSMIN'
YCOMMENT='minimum stomatal resistance (SM-1)'
 CALL WRITE_SURF(DGU, IOB, U, &
                  HPROGRAM,YRECFM,TGDP%XRSMIN(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='GD_GAMMA'
YCOMMENT='coefficient for RSMIN calculation (-)'
 CALL WRITE_SURF(DGU, IOB, U, &
                  HPROGRAM,YRECFM,TGDP%XGAMMA(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='GD_CV'
YCOMMENT='vegetation thermal inertia coefficient (-)'
 CALL WRITE_SURF(DGU, IOB, U, &
                  HPROGRAM,YRECFM,TGDP%XCV(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='GD_RGL'
YCOMMENT='maximum solar radiation usable in photosynthesis (-)'
 CALL WRITE_SURF(DGU, IOB, U, &
                  HPROGRAM,YRECFM,TGDP%XRGL(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='GD_EMIS_ISBA'
YCOMMENT='surface emissivity (-)'
 CALL WRITE_SURF(DGU, IOB, U, &
                  HPROGRAM,YRECFM,TGDPE%CUR%XEMIS(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='GD_WRMAX_CF'
YCOMMENT='coefficient for maximum water interception (-)'
 CALL WRITE_SURF(DGU, IOB, U, &
                  HPROGRAM,YRECFM,TGDP%XWRMAX_CF(:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
IF (DGMTO%LSURF_DIAG_ALBEDO) THEN
!
!* Soil albedos
!
!
   YRECFM='GD_ALBNIR_S'
   YCOMMENT='soil near-infra-red albedo (-)'
   CALL WRITE_SURF(DGU, IOB, U, &
                  HPROGRAM,YRECFM,TGDP%XALBNIR_SOIL(:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
   YRECFM='GD_ALBVIS_S'
   YCOMMENT='soil visible albedo (-)'
   CALL WRITE_SURF(DGU, IOB, U, &
                  HPROGRAM,YRECFM,TGDP%XALBVIS_SOIL(:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
   YRECFM='GD_ALBUV_S'
   YCOMMENT='soil UV albedo (-)'
   CALL WRITE_SURF(DGU, IOB, U, &
                  HPROGRAM,YRECFM,TGDP%XALBUV_SOIL(:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
!* albedos
!
   YRECFM='GD_ALBNIR_T'
   YCOMMENT='total near-infra-red albedo (-)'
   CALL WRITE_SURF(DGU, IOB, U, &
                  HPROGRAM,YRECFM,TGDPE%CUR%XALBNIR(:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
   YRECFM='GD_ALBVIS_T'
   YCOMMENT='total visible albedo (-)'
   CALL WRITE_SURF(DGU, IOB, U, &
                  HPROGRAM,YRECFM,TGDPE%CUR%XALBVIS(:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
   YRECFM='GD_ALBUV_T'
   YCOMMENT='total UV albedo (-)'
   CALL WRITE_SURF(DGU, IOB, U, &
                  HPROGRAM,YRECFM,TGDPE%CUR%XALBUV(:),IRESP,HCOMMENT=YCOMMENT)
!
END IF
!
!-------------------------------------------------------------------------------
!
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_PGD_GRDN_N',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE WRITE_DIAG_PGD_GRDN_n
