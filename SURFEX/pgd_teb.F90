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
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_DATA_COVER_PAR, ONLY : JPCOVER
USE MODD_TEB_n,          ONLY : XCOVER, LCOVER, XZS,                  &
                                  NROAD_LAYER, NWALL_LAYER, NROOF_LAYER, &
                                  LECOCLIMAP, LGARDEN  
USE MODD_TEB_GRID_n,     ONLY : CGRID, XGRID_PAR, XLAT, XLON, XMESH_SIZE, NDIM
!
USE MODI_GET_SURF_SIZE_n
USE MODI_PACK_PGD
USE MODI_PGD_TEB_PAR
USE MODI_PGD_TEB_GARDEN
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
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
!
!*    0.3    Declaration of namelists
!            ------------------------
!
!-------------------------------------------------------------------------------
!
!*    1.      Initializations of defaults
!             ---------------------------
!
IF (LHOOK) CALL DR_HOOK('PGD_TEB',0,ZHOOK_HANDLE)

NROOF_LAYER = 3
NROAD_LAYER = 3
NWALL_LAYER = 3
!
!-------------------------------------------------------------------------------
!
!*    2.      Reading of namelist
!             -------------------
!
!-------------------------------------------------------------------------------
!
!*    3.      Coherence of options
!             --------------------
!
!-------------------------------------------------------------------------------
!
!*    4.      Number of points and packing
!             ----------------------------
!
CALL GET_SURF_SIZE_n('TOWN  ',NDIM)
!
ALLOCATE(LCOVER     (JPCOVER))
ALLOCATE(XCOVER     (NDIM,JPCOVER))
ALLOCATE(XZS        (NDIM))
ALLOCATE(XLAT       (NDIM))
ALLOCATE(XLON       (NDIM))
ALLOCATE(XMESH_SIZE (NDIM))
!
CALL PACK_PGD(HPROGRAM, 'TOWN  ',                    &
                CGRID,  XGRID_PAR,                     &
                LCOVER, XCOVER, XZS,                   &
                XLAT, XLON, XMESH_SIZE                 )  
!
!-------------------------------------------------------------------------------
!
!*    9.      TEB specific fields
!             -------------------
!
LECOCLIMAP = OECOCLIMAP
CALL PGD_TEB_PAR(HPROGRAM,NROOF_LAYER,NROAD_LAYER,NWALL_LAYER,OGARDEN)
!
!-------------------------------------------------------------------------------
!
!*   10.     Prints of cover parameters in a tex file
!            ----------------------------------------
!
IF (OECOCLIMAP) CALL WRITE_COVER_TEX_TEB
!
!-------------------------------------------------------------------------------
!
!*   11.     Case of urban green areas
!            -------------------------
!
LGARDEN       = OGARDEN
!
IF (OGARDEN) CALL PGD_TEB_GARDEN(HPROGRAM,OECOCLIMAP)
IF (LHOOK) CALL DR_HOOK('PGD_TEB',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE PGD_TEB
