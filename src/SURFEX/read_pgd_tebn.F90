!     #########
      SUBROUTINE READ_PGD_TEB_n (BOP, BDD, DTB, DTCO, DTT, IOB, U, TG, TOP, &
                                 HPROGRAM)
!     #########################################
!
!!****  *READ_PGD_TEB_n* - reads TEB physiographic fields
!!                       
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
USE MODD_DATA_BEM_n, ONLY : DATA_BEM_t
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DATA_TEB_n, ONLY : DATA_TEB_t
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_TEB_GRID_n, ONLY : TEB_GRID_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
!
USE MODD_TYPE_DATE_SURF
!
USE MODD_DATA_COVER_PAR, ONLY : JPCOVER
!
!
USE MODE_READ_SURF_COV, ONLY : READ_SURF_COV
!
USE MODI_READ_SURF
USE MODI_READ_GRID
USE MODI_READ_LCOVER
USE MODI_READ_PGD_TEB_PAR_n
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_GET_TYPE_DIM_n
!
USE MODI_READ_LECOCLIMAP
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
TYPE(BEM_OPTIONS_t), INTENT(INOUT) :: BOP
TYPE(BLD_DESC_t), INTENT(INOUT) :: BDD
TYPE(DATA_BEM_t), INTENT(INOUT) :: DTB
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DATA_TEB_t), INTENT(INOUT) :: DTT
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(TEB_GRID_t), INTENT(INOUT) :: TG
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! calling program
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! Error code after redding
!
 CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
INTEGER           :: IVERSION
INTEGER           :: IBUGFIX
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!* 1D physical dimension
!
IF (LHOOK) CALL DR_HOOK('READ_PGD_TEB_N',0,ZHOOK_HANDLE)
YRECFM='SIZE_TOWN'
 CALL GET_TYPE_DIM_n(DTCO, U, &
                     'TOWN  ',TG%NDIM)
!
!*       2.     Other dimension initializations:
!               --------------------------------
!
 CALL READ_SURF(IOB, &
                HPROGRAM,'VERSION',IVERSION,IRESP)
 CALL READ_SURF(IOB, &
                HPROGRAM,'BUG',IBUGFIX,IRESP)
!
!* number of TEB patches
!
IF (IVERSION<7 .OR. IVERSION==7 .AND. IBUGFIX<=2) THEN
  TOP%NTEB_PATCH=1
ELSE
  YRECFM='TEB_PATCH'
  CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,TOP%NTEB_PATCH,IRESP)
END IF
!
!* number of road and roof layers
!
YRECFM='ROAD_LAYER'
 CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,TOP%NROAD_LAYER,IRESP)

YRECFM='ROOF_LAYER'
 CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,TOP%NROOF_LAYER,IRESP)

YRECFM='WALL_LAYER'
 CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,TOP%NWALL_LAYER,IRESP)
!
!
!* type of averaging for Buildings (to allow ascendant compatibility)
!* type of Building Energy Model
!
IF (IVERSION<7 .OR.( IVERSION==7 .AND. IBUGFIX<=2)) THEN
  TOP%CBLD_ATYPE='ARI'
  TOP%CBEM = 'DEF'
ELSE
  YRECFM='BLD_ATYPE'
  CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,TOP%CBLD_ATYPE,IRESP)
  YRECFM='BEM'
  CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,TOP%CBEM,IRESP)
END IF
!
IF (TOP%CBEM=="BEM") THEN
  YRECFM='FLOOR_LAYER'
  CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,BOP%NFLOOR_LAYER,IRESP)
  YRECFM='COOL_COIL'
  CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,BOP%CCOOL_COIL,IRESP)
  YRECFM='HEAT_COIL'
  CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,BOP%CHEAT_COIL,IRESP)
  YRECFM='AUTOSIZE'
  CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,BOP%LAUTOSIZE,IRESP)
ENDIF
!
!* Case of urban green roofs
!
IF (TOP%LGARDEN) THEN
  IF (IVERSION<7 .OR.( IVERSION==7 .AND. IBUGFIX<=2)) THEN
    TOP%LGREENROOF = .FALSE.
  ELSE
    YRECFM='LGREENROOF'
    CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,TOP%LGREENROOF,IRESP)
  END IF
!
!* Case of urban hydrology
!
  IF (IVERSION<7 .OR.( IVERSION==7 .AND. IBUGFIX<=3)) THEN
    TOP%LHYDRO = .FALSE.
  ELSE
    YRECFM='LURBAN_HYDRO'
    CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,TOP%LHYDRO,IRESP)
  END IF
ENDIF
!
!* Solar panels
!
IF (IVERSION<7 .OR.( IVERSION==7 .AND. IBUGFIX<=3)) THEN
  TOP%LSOLAR_PANEL = .FALSE.
ELSE
  YRECFM='SOLAR_PANEL'
  CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,TOP%LSOLAR_PANEL,IRESP)
END IF
!
!
!
!*       3.     Physiographic data fields:
!               -------------------------
!
!* cover classes
!
ALLOCATE(TOP%LCOVER(JPCOVER))
 CALL READ_LCOVER(IOB, &
                  HPROGRAM,TOP%LCOVER)
!
ALLOCATE(TOP%XCOVER(TG%NDIM,COUNT(TOP%LCOVER)))
 CALL READ_SURF_COV(IOB, &
                    HPROGRAM,'COVER',TOP%XCOVER(:,:),TOP%LCOVER,IRESP)
!
!* orography
!
ALLOCATE(TOP%XZS(TG%NDIM))
YRECFM='ZS'
 CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,TOP%XZS(:),IRESP)
!
!
!* latitude, longitude 
!
ALLOCATE(TG%XLAT      (TG%NDIM))
ALLOCATE(TG%XLON      (TG%NDIM))
ALLOCATE(TG%XMESH_SIZE(TG%NDIM))
 CALL READ_GRID(IOB, &
                HPROGRAM,TG%CGRID,TG%XGRID_PAR,TG%XLAT,TG%XLON,TG%XMESH_SIZE,IRESP)
!
!
!-------------------------------------------------------------------------------
!
!*       4.     Physiographic data fields not to be computed by ecoclimap
!               ---------------------------------------------------------
!
 CALL READ_LECOCLIMAP(IOB, &
                      HPROGRAM,TOP%LECOCLIMAP)
!
 CALL READ_PGD_TEB_PAR_n(DTCO, IOB, U, &
                         BDD, DTB, DTT, TG, TOP, &
                         HPROGRAM,TG%NDIM,'-')
IF (LHOOK) CALL DR_HOOK('READ_PGD_TEB_N',1,ZHOOK_HANDLE)
!
!
!------------------------------------------------------------------------------
!
END SUBROUTINE READ_PGD_TEB_n
