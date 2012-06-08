!     #########
      SUBROUTINE READ_PGD_TEB_n(HPROGRAM)
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
!!	V. Masson   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2003 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_TYPE_DATE_SURF
!
USE MODD_DATA_COVER_PAR, ONLY : JPCOVER
!
USE MODD_TEB_n,          ONLY : XCOVER, XZS,                           &
                                  NROOF_LAYER, NROAD_LAYER, NWALL_LAYER, &
                                  TTIME, LCOVER, LECOCLIMAP  
USE MODD_TEB_GRID_n,     ONLY : XLAT, XLON, XMESH_SIZE, CGRID, XGRID_PAR, NDIM
!
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
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! calling program
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! Error code after redding
!
CHARACTER(LEN=16) :: YRECFM         ! Name of the article to be read
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!* 1D physical dimension
!
IF (LHOOK) CALL DR_HOOK('READ_PGD_TEB_N',0,ZHOOK_HANDLE)
YRECFM='SIZE_TOWN'
CALL GET_TYPE_DIM_n('TOWN  ',NDIM)
!
!*       2.     Other dimension initializations:
!               --------------------------------
!
!* number of road and roof layers
!
YRECFM='ROAD_LAYER'
CALL READ_SURF(HPROGRAM,YRECFM,NROAD_LAYER,IRESP)

YRECFM='ROOF_LAYER'
CALL READ_SURF(HPROGRAM,YRECFM,NROOF_LAYER,IRESP)

YRECFM='WALL_LAYER'
CALL READ_SURF(HPROGRAM,YRECFM,NWALL_LAYER,IRESP)
!
!
!
!
!*       3.     Physiographic data fields:
!               -------------------------
!
!* cover classes
!
ALLOCATE(LCOVER(JPCOVER))
CALL READ_LCOVER(HPROGRAM,LCOVER)
!
ALLOCATE(XCOVER(NDIM,JPCOVER))
CALL READ_SURF(HPROGRAM,'COVER',XCOVER(:,:),LCOVER,IRESP)
!
!* orography
!
ALLOCATE(XZS(NDIM))
YRECFM='ZS'
CALL READ_SURF(HPROGRAM,YRECFM,XZS(:),IRESP)
!
!
!* latitude, longitude 
!
ALLOCATE(XLAT      (NDIM))
ALLOCATE(XLON      (NDIM))
ALLOCATE(XMESH_SIZE(NDIM))
CALL READ_GRID(HPROGRAM,CGRID,XGRID_PAR,XLAT,XLON,XMESH_SIZE,IRESP)
!
!
!-------------------------------------------------------------------------------
!
!*       4.     Physiographic data fields not to be computed by ecoclimap
!               ---------------------------------------------------------
!
CALL READ_LECOCLIMAP(HPROGRAM,LECOCLIMAP)
!
CALL READ_PGD_TEB_PAR_n(HPROGRAM,NDIM,'-')
IF (LHOOK) CALL DR_HOOK('READ_PGD_TEB_N',1,ZHOOK_HANDLE)
!
!
!------------------------------------------------------------------------------
!
END SUBROUTINE READ_PGD_TEB_n
