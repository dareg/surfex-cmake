!     #########
      SUBROUTINE PGD_TEB (BDD, BOP, DTB, DTCO, DTI, DTGD, DTGR, DTT, DGU, IOB, &
                          UG, U, USS, TGDO, TGDP, TGRO, TGRP, TG, TIR, T, TOP, TVG, &
                          HPROGRAM,OECOCLIMAP,OGARDEN)
!     ##############################################################
!
!!**** *PGD_TEB* monitor for averaging and interpolations of TEB physiographic fields
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
!!    V. Masson        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    10/12/97
!!    A. Lemonsu      05/2009         Key for garden option
!!    G. Pigeon     /09/12: WALL, ROOF, FLOOR, MASS LAYER default to 5
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
!
USE MODD_BLD_DESCRIPTION_n, ONLY : BLD_DESC_t
USE MODD_BEM_OPTION_n, ONLY : BEM_OPTIONS_t
USE MODD_DATA_BEM_n, ONLY : DATA_BEM_t
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_DATA_TEB_GARDEN_n, ONLY : DATA_TEB_GARDEN_t
USE MODD_DATA_TEB_GREENROOF_n, ONLY : DATA_TEB_GREENROOF_t
USE MODD_DATA_TEB_n, ONLY : DATA_TEB_t
USE MODD_DIAG_SURF_ATM_n, ONLY : DIAG_SURF_ATM_t
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SURF_ATM_SSO_n, ONLY : SURF_ATM_SSO_t
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TEB_GARDEN_OPTIONS_t
USE MODD_TEB_GARDEN_PGD_n, ONLY : TEB_GARDEN_PGD_t
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TEB_GREENROOF_OPTIONS_t
USE MODD_TEB_GREENROOF_PGD_n, ONLY : TEB_GREENROOF_PGD_t
USE MODD_TEB_GRID_n, ONLY : TEB_GRID_t
USE MODD_TEB_IRRIG_n, ONLY : TEB_IRRIG_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_OPTIONS_t
!
USE MODD_DATA_COVER_PAR, ONLY : JPCOVER
!
USE MODI_GET_SURF_SIZE_n
USE MODI_PACK_PGD
USE MODI_PGD_TEB_PAR
USE MODI_PGD_TEB_VEG
USE MODI_GET_LUOUT
USE MODI_READ_NAM_PGD_TEB
USE MODI_TEST_NAM_VAR_SURF
USE MODI_PGD_BEM_PAR
USE MODI_ABOR1_SFX
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_WRITE_COVER_TEX_TEB
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
!
TYPE(BLD_DESC_t), INTENT(INOUT) :: BDD
TYPE(BEM_OPTIONS_t), INTENT(INOUT) :: BOP
TYPE(DATA_BEM_t), INTENT(INOUT) :: DTB
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTI
TYPE(DATA_TEB_GARDEN_t), INTENT(INOUT) :: DTGD
TYPE(DATA_TEB_GREENROOF_t), INTENT(INOUT) :: DTGR
TYPE(DATA_TEB_t), INTENT(INOUT) :: DTT
TYPE(DIAG_SURF_ATM_t), INTENT(INOUT) :: DGU
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SURF_ATM_SSO_t), INTENT(INOUT) :: USS
TYPE(TEB_GARDEN_OPTIONS_t), INTENT(INOUT) :: TGDO
TYPE(TEB_GARDEN_PGD_t), INTENT(INOUT) :: TGDP
TYPE(TEB_GREENROOF_OPTIONS_t), INTENT(INOUT) :: TGRO
TYPE(TEB_GREENROOF_PGD_t), INTENT(INOUT) :: TGRP
TYPE(TEB_GRID_t), INTENT(INOUT) :: TG
TYPE(TEB_IRRIG_t), INTENT(INOUT) :: TIR
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_VEG_OPTIONS_t), INTENT(INOUT) :: TVG
!
 CHARACTER(LEN=6), INTENT(IN)  :: HPROGRAM   ! Type of program
LOGICAL,          INTENT(IN)  :: OECOCLIMAP ! T if parameters are computed with ecoclimap
!                                           ! F if all parameters must be specified
LOGICAL,          INTENT(IN)  :: OGARDEN    ! T if urban green areas
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER         :: ILUOUT    ! output listing logical unit
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*    1.      Initializations of defaults
!             ---------------------------
!
IF (LHOOK) CALL DR_HOOK('PGD_TEB',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)

TOP%NROOF_LAYER  = 5
TOP%NROAD_LAYER  = 5
TOP%NWALL_LAYER  = 5
BOP%NFLOOR_LAYER = 5
!
!-------------------------------------------------------------------------------
!
!*    2.      Reading of namelist
!             -------------------
!
 CALL READ_NAM_PGD_TEB(HPROGRAM,TOP%NTEB_PATCH,TOP%CBEM,BOP%CCOOL_COIL,BOP%CHEAT_COIL,BOP%LAUTOSIZE,&
                      TOP%NROAD_LAYER,TOP%NROOF_LAYER,TOP%NWALL_LAYER,BOP%NFLOOR_LAYER,        &
                      TOP%LGREENROOF,TOP%LHYDRO,TOP%LSOLAR_PANEL                           )
!
!-------------------------------------------------------------------------------
!
!*    3.      Coherence of options
!             --------------------
!
 CALL TEST_NAM_VAR_SURF(ILUOUT,'CBLD',TOP%CBEM,'DEF','BEM ')
 CALL TEST_NAM_VAR_SURF(ILUOUT,'CCOOL_COIL',BOP%CCOOL_COIL,'IDEAL ','DXCOIL')
 CALL TEST_NAM_VAR_SURF(ILUOUT,'CHEAT_COIL',BOP%CHEAT_COIL,'IDEAL ','FINCAP')
!
IF (.NOT. OGARDEN) THEN
  IF (TOP%LGREENROOF) CALL ABOR1_SFX('ERROR: You cannot activate LGREENROOF if LGARDEN is FALSE')
  IF (TOP%LHYDRO    ) CALL ABOR1_SFX('ERROR: You cannot activate LHYDRO     if LGARDEN is FALSE')
ENDIF
!
!-------------------------------------------------------------------------------
!
!*    4.      Number of points and packing
!             ----------------------------
!
 CALL GET_SURF_SIZE_n(DTCO, U, &
                      'TOWN  ',TG%NDIM)
!
ALLOCATE(TOP%LCOVER     (JPCOVER))
ALLOCATE(TOP%XCOVER     (TG%NDIM,JPCOVER))
ALLOCATE(TOP%XZS        (TG%NDIM))
ALLOCATE(TG%XLAT       (TG%NDIM))
ALLOCATE(TG%XLON       (TG%NDIM))
ALLOCATE(TG%XMESH_SIZE (TG%NDIM))
!
 CALL PACK_PGD(DTCO, U, &
               HPROGRAM, 'TOWN  ',                    &
                TG%CGRID,  TG%XGRID_PAR,                   &
                TOP%LCOVER, TOP%XCOVER, TOP%XZS,                 &
                TG%XLAT, TG%XLON, TG%XMESH_SIZE               )  
!
!-------------------------------------------------------------------------------
!
!*    5.      TEB specific fields
!             -------------------
!
TOP%LECOCLIMAP = OECOCLIMAP
 CALL PGD_TEB_PAR(DTCO, DGU, IOB, UG, U, USS, BDD, DTT, DTI, TG, &
                  HPROGRAM,OGARDEN,TOP%LGREENROOF,TOP%CBLD_ATYPE)
!
!-------------------------------------------------------------------------------
!
!*    6.      Prints of cover parameters in a tex file
!             ----------------------------------------
!
IF (OECOCLIMAP) CALL WRITE_COVER_TEX_TEB
!
!
!-------------------------------------------------------------------------------
!
!*    7.      Case of urban green areas (and hydrology)
!             -----------------------------------------
!
TOP%LGARDEN       = OGARDEN
!
IF (TOP%LGARDEN) CALL PGD_TEB_VEG(DTCO, DTGD, DTGR, UG, U, USS, TGDO, TGDP, TGRO, TGRP, TG, &
                              TIR, TOP, TVG, &
                                  HPROGRAM)
!
!-------------------------------------------------------------------------------
!
!*    8.      Case of Building Energy Model
!             -----------------------------
!
IF (TOP%CBEM .EQ. 'BEM') CALL PGD_BEM_PAR(DTCO, DGU, IOB, UG, U, USS, DTB, DTI, TG, &
                                          HPROGRAM,BOP%LAUTOSIZE)
!
IF (LHOOK) CALL DR_HOOK('PGD_TEB',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE PGD_TEB
