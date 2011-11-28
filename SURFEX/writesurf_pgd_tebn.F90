!     #########
      SUBROUTINE WRITESURF_PGD_TEB_n(HPROGRAM,HWRITE)
!     ###############################################
!
!!****  *WRITE_PGD_TEB_n* - writes TEB fields
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
!
USE MODD_TEB_n,          ONLY : NROOF_LAYER, NROAD_LAYER, NWALL_LAYER, &
                                  XZS,XCOVER, LCOVER, LECOCLIMAP, LGARDEN  
USE MODD_TEB_GRID_n,     ONLY : XLAT, XLON, XMESH_SIZE, CGRID, XGRID_PAR
USE MODD_DATA_COVER_PAR, ONLY : JPCOVER
USE MODD_TEB_GARDEN_n,   ONLY : CISBA, CDIF, CPEDOTF, CPHOTO,      &
                                NGROUND_LAYER, NNBIOMASS, &
                                XCLAY, XSAND, XRUNOFFB, XWDRAIN  
!
USE MODI_WRITE_SURF
USE MODI_WRITE_GRID
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_WRITESURF_PGD_TEB_PAR_n
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling
CHARACTER(LEN=3),  INTENT(IN)  :: HWRITE   ! 'PGD' : only physiographic fields are written
!                                          ! 'ALL' : all fields are written

!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
CHARACTER(LEN=16) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=100):: YCOMMENT       ! Comment string
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*       1.     Dimension initializations:
!               -------------------------
!
!
!* number of roof layers
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_PGD_TEB_N',0,ZHOOK_HANDLE)
YRECFM='ROOF_LAYER'
YCOMMENT=YRECFM
CALL WRITE_SURF(HPROGRAM,YRECFM,NROOF_LAYER,IRESP,HCOMMENT=YCOMMENT)
!
!* number of road layers
!
YRECFM='ROAD_LAYER'
YCOMMENT=YRECFM
CALL WRITE_SURF(HPROGRAM,YRECFM,NROAD_LAYER,IRESP,HCOMMENT=YCOMMENT)
!
!* number of wall layers
!
YRECFM='WALL_LAYER'
YCOMMENT=YRECFM
CALL WRITE_SURF(HPROGRAM,YRECFM,NWALL_LAYER,IRESP,HCOMMENT=YCOMMENT)
!
!* flag indicating if fields are computed from ecoclimap or not
!
YRECFM='ECOCLIMAP'
YCOMMENT=YRECFM
CALL WRITE_SURF(HPROGRAM,YRECFM,LECOCLIMAP,IRESP,HCOMMENT=YCOMMENT)
!
!
!------------------------------------------------------------------------------
!
! * ISBA fields for urban green areas
! 
IF (LGARDEN) THEN
!        
!* soil scheme option
!
YRECFM='TWN_ISBA'
YCOMMENT=YRECFM
CALL WRITE_SURF(HPROGRAM,YRECFM,CISBA,IRESP,HCOMMENT=YCOMMENT)
!
!* Theta(psi) function
!
YRECFM='TWN_THETAPSI'
YCOMMENT=YRECFM
CALL WRITE_SURF(HPROGRAM,YRECFM,CDIF,IRESP,HCOMMENT=YCOMMENT)
!
!* Pedo-transfert function
!
YRECFM='TWN_PEDOTF'
YCOMMENT=YRECFM
CALL WRITE_SURF(HPROGRAM,YRECFM,CPEDOTF,IRESP,HCOMMENT=YCOMMENT)
!
!* type of photosynthesis
!
YRECFM='TWN_PHOTO'
YCOMMENT=YRECFM
CALL WRITE_SURF(HPROGRAM,YRECFM,CPHOTO,IRESP,HCOMMENT=YCOMMENT)
!
!* number of soil layers
!
YRECFM='TWN_LAYER'
YCOMMENT=YRECFM
CALL WRITE_SURF(HPROGRAM,YRECFM,NGROUND_LAYER,IRESP,HCOMMENT=YCOMMENT)
!
!* number of biomass pools
!
YRECFM='TWN_NBIOMASS'
YCOMMENT=YRECFM
CALL WRITE_SURF(HPROGRAM,YRECFM,NNBIOMASS,IRESP,HCOMMENT=YCOMMENT)
!
! * clay fraction
!
YRECFM='TWN_CLAY'
YCOMMENT='X_Y_TWN_CLAY'
CALL WRITE_SURF(HPROGRAM,YRECFM,XCLAY(:,1),IRESP,HCOMMENT=YCOMMENT)
!        
! * sand fraction
!
YRECFM='TWN_SAND'
YCOMMENT='X_Y_TWN_SAND'
CALL WRITE_SURF(HPROGRAM,YRECFM,XSAND(:,1),IRESP,HCOMMENT=YCOMMENT)
!        
! * orographic runoff coefficient
!
YRECFM='TWN_RUNOFFB'
YCOMMENT='X_Y_TWN_RUNOFFB'
CALL WRITE_SURF(HPROGRAM,YRECFM,XRUNOFFB,IRESP,HCOMMENT=YCOMMENT)
!        
! * subgrid drainage coefficient
!
YRECFM='TWN_WDRAIN'
YCOMMENT='X_Y_TWN_WDRAIN'
CALL WRITE_SURF(HPROGRAM,YRECFM,XWDRAIN,IRESP,HCOMMENT=YCOMMENT)
!
ENDIF
!
!------------------------------------------------------------------------------
!
!*       2.     Physiographic data fields:
!               -------------------------
!
!* cover classes
!
YRECFM='COVER_LIST'
YCOMMENT='(LOGICAL LIST)'
CALL WRITE_SURF(HPROGRAM,YRECFM,LCOVER(:),IRESP,HCOMMENT=YCOMMENT,HDIR='-')
!
YCOMMENT='COVER FIELDS'
CALL WRITE_SURF(HPROGRAM,'COVER',XCOVER(:,:),LCOVER,IRESP,HCOMMENT=YCOMMENT)
!
!* orography
!
YRECFM='ZS'
YCOMMENT='ZS'
CALL WRITE_SURF(HPROGRAM,YRECFM,XZS(:),IRESP,HCOMMENT=YCOMMENT)
!
!* latitude, longitude
!
CALL WRITE_GRID(HPROGRAM,CGRID,XGRID_PAR,XLAT,XLON,XMESH_SIZE,IRESP)
!
!-------------------------------------------------------------------------------
CALL WRITESURF_PGD_TEB_PAR_n(HPROGRAM)
!
IF (LHOOK) CALL DR_HOOK('WRITESURF_PGD_TEB_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE WRITESURF_PGD_TEB_n
