!     ########################################
      SUBROUTINE GET_BLD_CONF_n(ODATA_BLDTYPE, ODATA_BLD_AGE, ODATA_USETYPE, &
                KDESC_ROOF_LAYER, KDESC_ROAD_LAYER, KDESC_WALL_LAYER, &
                KDESC_FLOOR_LAYER, KDESC_CODE, KDESC_USE, KDESC_AGE, KDESC_BLD)  
!     ########################################
!
!!****  *GET_BLD_CONF_n* - routine to get some ISBA fields
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
!!      Original    04/2008
!!      A.L. Gibelin 07/2009 : Dimensions for carbon options
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DATA_TEB_n, ONLY : LDATA_BLDTYPE, LDATA_BLD_AGE, LDATA_USETYPE
USE MODD_BLD_DESCRIPTION_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
LOGICAL, INTENT(OUT) :: ODATA_BLDTYPE
LOGICAL, INTENT(OUT) :: ODATA_BLD_AGE
LOGICAL, INTENT(OUT) :: ODATA_USETYPE
INTEGER, INTENT(OUT) :: KDESC_ROOF_LAYER
INTEGER, INTENT(OUT) :: KDESC_ROAD_LAYER
INTEGER, INTENT(OUT) :: KDESC_WALL_LAYER
INTEGER, INTENT(OUT) :: KDESC_FLOOR_LAYER
INTEGER, INTENT(OUT) :: KDESC_CODE
INTEGER, INTENT(OUT) :: KDESC_USE
INTEGER, INTENT(OUT) :: KDESC_AGE
INTEGER, INTENT(OUT) :: KDESC_BLD
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GET_BLD_CONF_N',0,ZHOOK_HANDLE)
ODATA_BLDTYPE = LDATA_BLDTYPE
ODATA_BLD_AGE = LDATA_BLD_AGE
ODATA_USETYPE = LDATA_USETYPE
KDESC_ROOF_LAYER = NDESC_ROOF_LAYER
KDESC_ROAD_LAYER = NDESC_ROAD_LAYER
KDESC_WALL_LAYER = NDESC_WALL_LAYER
KDESC_FLOOR_LAYER = NDESC_FLOOR_LAYER
KDESC_CODE = NDESC_CODE
KDESC_USE = NDESC_USE
KDESC_AGE = NDESC_AGE
KDESC_BLD = NDESC_BLD
IF (LHOOK) CALL DR_HOOK('GET_BLD_CONF_N',1,ZHOOK_HANDLE)
!
!==============================================================================
!
END SUBROUTINE GET_BLD_CONF_n
