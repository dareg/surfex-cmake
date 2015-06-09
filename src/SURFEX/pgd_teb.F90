!     #########
      SUBROUTINE PGD_TEB(HPROGRAM,OECOCLIMAP,OGARDEN)
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
USE MODD_BLD_DESCRIPTION_n, ONLY : BDD => BLD_DESC
USE MODD_DATA_TEB_n, ONLY : DTT => DATA_TEB
!
USE MODD_DATA_BEM_n, ONLY : DTB => DATA_BEM
USE MODD_DATA_ISBA_n, ONLY : DTI => DATA_ISBA
!
USE MODD_DATA_COVER_PAR, ONLY : JPCOVER
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS
USE MODD_BEM_OPTION_n, ONLY : BOP => BEM_OPTIONS
USE MODD_TEB_GRID_n, ONLY : TG => TEB_GRID
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
 CALL GET_SURF_SIZE_n('TOWN  ',TG%NDIM)
!
ALLOCATE(TOP%LCOVER     (JPCOVER))
ALLOCATE(TOP%XCOVER     (TG%NDIM,JPCOVER))
ALLOCATE(TOP%XZS        (TG%NDIM))
ALLOCATE(TG%XLAT       (TG%NDIM))
ALLOCATE(TG%XLON       (TG%NDIM))
ALLOCATE(TG%XMESH_SIZE (TG%NDIM))
!
 CALL PACK_PGD(HPROGRAM, 'TOWN  ',                    &
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
 CALL PGD_TEB_PAR(BDD, DTI, DTT, TG, &
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
IF (TOP%LGARDEN) CALL PGD_TEB_VEG(HPROGRAM)
!
!-------------------------------------------------------------------------------
!
!*    8.      Case of Building Energy Model
!             -----------------------------
!
IF (TOP%CBEM .EQ. 'BEM') CALL PGD_BEM_PAR(DTB, DTI, TG, &
                                          HPROGRAM,BOP%LAUTOSIZE)
!
IF (LHOOK) CALL DR_HOOK('PGD_TEB',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE PGD_TEB
