!     #########
      SUBROUTINE WRITESURF_PGD_SEAFLUX_n(HPROGRAM)
!     ###################################################
!
!!****  *WRITE_SEAFLUX_n* - writes SEAFLUX fields
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
!!      B. Decharme 07/2011 : delete argument HWRITE
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SEAFLUX_n, ONLY : S => SEAFLUX
USE MODD_SEAFLUX_GRID_n, ONLY : SG => SEAFLUX_GRID
USE MODD_DATA_SEAFLUX_n, ONLY : DTS => DATA_SEAFLUX
!
USE MODD_DATA_COVER_PAR, ONLY : JPCOVER
!
USE MODE_WRITE_SURF_COV, ONLY : WRITE_SURF_COV
!
USE MODI_WRITE_SURF
USE MODI_WRITE_GRID
USE MODI_WRITESURF_PGD_SEAF_PAR_n
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
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
 CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
 CHARACTER(LEN=100):: YCOMMENT       ! Comment string
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!
!
!*       2.     Physiographic data fields:
!               -------------------------
!
!* cover classes
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_PGD_SEAFLUX_N',0,ZHOOK_HANDLE)
YRECFM='COVER_LIST'
YCOMMENT='(LOGICAL LIST)'
 CALL WRITE_SURF(HPROGRAM,YRECFM,S%LCOVER(:),IRESP,HCOMMENT=YCOMMENT,HDIR='-')
!
YCOMMENT='COVER FIELDS'
 CALL WRITE_SURF_COV(HPROGRAM,'COVER',S%XCOVER(:,:),S%LCOVER,IRESP,HCOMMENT=YCOMMENT)
!
!
!* orography
!
YRECFM='ZS'
YCOMMENT='ZS'
 CALL WRITE_SURF(HPROGRAM,YRECFM,S%XZS(:),IRESP,HCOMMENT=YCOMMENT)
!
!* bathymetry
!
YRECFM='BATHY'
YCOMMENT='BATHY'
 CALL WRITE_SURF(HPROGRAM,YRECFM,S%XSEABATHY(:),IRESP,HCOMMENT=YCOMMENT)
!
!* latitude, longitude
!
 CALL WRITE_GRID(HPROGRAM,SG%CGRID,SG%XGRID_PAR,SG%XLAT,SG%XLON,SG%XMESH_SIZE,IRESP)
!
!* sst
!
YRECFM='SST_DATA'
YCOMMENT='(LOGICAL)'
 CALL WRITE_SURF(HPROGRAM,YRECFM,DTS%LSST_DATA,IRESP,HCOMMENT=YCOMMENT)
!
IF (DTS%LSST_DATA) CALL WRITESURF_PGD_SEAF_PAR_n(DTS, &
                                                 HPROGRAM)
IF (LHOOK) CALL DR_HOOK('WRITESURF_PGD_SEAFLUX_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITESURF_PGD_SEAFLUX_n
