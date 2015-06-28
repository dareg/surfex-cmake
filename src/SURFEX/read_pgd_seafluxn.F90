!     #########
      SUBROUTINE READ_PGD_SEAFLUX_n(HPROGRAM)
!     #########################################
!
!!****  *READ_PGD_SEAFLUX_n* - routine to read SEAFLUX physiographic fields
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
USE MODD_IO_BUFF_n, ONLY : IOB => IO_BUFF
!
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
USE MODD_TYPE_DATE_SURF
!
USE MODD_DATA_COVER_PAR, ONLY : JPCOVER
USE MODD_SEAFLUX_n, ONLY : S => SEAFLUX
USE MODD_SEAFLUX_GRID_n, ONLY : SG => SEAFLUX_GRID
USE MODD_DATA_SEAFLUX_n, ONLY : DTS => DATA_SEAFLUX
!
USE MODE_READ_SURF_COV, ONLY : READ_SURF_COV
!
USE MODI_READ_SURF
USE MODI_READ_GRID
USE MODI_READ_LCOVER
USE MODI_READ_PGD_SEAFLUX_PAR_n
!
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_GET_TYPE_DIM_n
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
 CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
!
INTEGER           :: IVERSION   ! surface version
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!* 1D physical dimension
!
IF (LHOOK) CALL DR_HOOK('READ_PGD_SEAFLUX_N',0,ZHOOK_HANDLE)
YRECFM='SIZE_SEA'
 CALL GET_TYPE_DIM_n(DTCO, U, &
                     'SEA   ',SG%NDIM)
!
!
!*       2.     Physiographic data fields:
!               -------------------------
!
!* cover classes
!
ALLOCATE(S%LCOVER(JPCOVER))
 CALL READ_LCOVER(HPROGRAM,S%LCOVER)
!
ALLOCATE(S%XCOVER(SG%NDIM,JPCOVER))
 CALL READ_SURF_COV(IOB, &
                    HPROGRAM,'COVER',S%XCOVER(:,:),S%LCOVER,IRESP)
!
!* orography
!
ALLOCATE(S%XZS(SG%NDIM))
S%XZS(:) = 0.
!
YRECFM='VERSION'
 CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,IVERSION,IRESP)
!
!* bathymetry
!
ALLOCATE(S%XSEABATHY(SG%NDIM))
IF (IVERSION<=3) THEN
  S%XSEABATHY(:) = -300.
ELSE
  YRECFM='BATHY'
  CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,S%XSEABATHY(:),IRESP)
END IF
!
!* latitude, longitude 
!
ALLOCATE(SG%XLAT      (SG%NDIM))
ALLOCATE(SG%XLON      (SG%NDIM))
ALLOCATE(SG%XMESH_SIZE(SG%NDIM))
 CALL READ_GRID(IOB, &
                HPROGRAM,SG%CGRID,SG%XGRID_PAR,SG%XLAT,SG%XLON,SG%XMESH_SIZE,IRESP)
!
!
!* sst
!
!
IF (IVERSION<3) THEN
  DTS%LSST_DATA = .FALSE.
ELSE
  YRECFM='SST_DATA'
  CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,DTS%LSST_DATA,IRESP)
END IF
!
IF (DTS%LSST_DATA) CALL READ_PGD_SEAFLUX_PAR_n(IOB, &
                                               DTS, SG, &
                                               HPROGRAM,SG%NDIM)
IF (LHOOK) CALL DR_HOOK('READ_PGD_SEAFLUX_N',1,ZHOOK_HANDLE)
!
!------------------------------------------------------------------------------!
END SUBROUTINE READ_PGD_SEAFLUX_n
