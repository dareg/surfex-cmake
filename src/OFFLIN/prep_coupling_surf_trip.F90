!     #########
      SUBROUTINE PREP_COUPLING_SURF_TRIP(HPROGRAM,KNI,OFLOOD_ISBA,HGRID)
!     ##################################################################
!
!!****  *PREP_COUPLING_SURF_TRIP* - routine to prepare the SURFACE-TRIP coupling
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
!!	B. Decharme   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/2008
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODI_GET_GRID_CONF_ISBA_n
USE MODI_READ_NAM_GRID_TRIP
USE MODI_GET_GRID_CONF_TRIP_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
USE MODI_GET_LUOUT
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=6),                 INTENT(IN)  :: HPROGRAM    ! program calling surf. schemes
 CHARACTER(LEN=10),                INTENT(IN)  :: HGRID       
LOGICAL,                          INTENT(IN)  :: OFLOOD_ISBA
INTEGER,                          INTENT(IN)  :: KNI         ! Surfex grid dimension
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
LOGICAL :: LFLOOD_TRIP
!
REAL    :: ZLONMIN_ISBA
REAL    :: ZLONMIN_TRIP
REAL    :: ZLONMAX_ISBA
REAL    :: ZLONMAX_TRIP
REAL    :: ZLATMIN_ISBA
REAL    :: ZLATMIN_TRIP
REAL    :: ZLATMAX_ISBA
REAL    :: ZLATMAX_TRIP
REAL    :: ZRES_ISBA
REAL    :: ZRES_TRIP
!
INTEGER :: ICOMP_TRIP
INTEGER :: ICOMP_ISBA
!
INTEGER :: ILON_ISBA
INTEGER :: ILON_TRIP
INTEGER :: ILAT_ISBA
INTEGER :: ILAT_TRIP
!
INTEGER :: ILUOUT, JL, K, I, J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PREP_COUPLING_SURF_TRIP',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
!*       1.0   Test ISBA - TRIP grid compatibility
!              -------------------------------------
!
IF(HGRID/="LONLAT REG")THEN
  CALL ABOR1_SFX('PREP_COUPLING_SURF_TRIPN: ISBA-TRIP REQUIRE LONLAT REG GRID')
ENDIF
!
 CALL READ_NAM_GRID_TRIP(HPROGRAM)
 CALL GET_GRID_CONF_TRIP_n(ZLONMIN_TRIP,ZLONMAX_TRIP,ZLATMIN_TRIP,ZLATMAX_TRIP,ZRES_TRIP,ILON_TRIP,ILAT_TRIP)
!
 CALL GET_GRID_CONF_ISBA_n(ZLONMIN_ISBA,ZLONMAX_ISBA,ZLATMIN_ISBA,ZLATMAX_ISBA,ILON_ISBA,ILAT_ISBA,JL)
ZRES_ISBA = (ZLONMAX_ISBA - ZLONMIN_ISBA) / ILON_ISBA
!
IF(ZRES_TRIP/=ZRES_ISBA)THEN
  WRITE(ILUOUT,*)'PREP_COUPLING_SURF_TRIPN: TRIP RESOLUTION = ',ZRES_TRIP
  WRITE(ILUOUT,*)'PREP_COUPLING_SURF_TRIPN: ISBA RESOLUTION = ',ZRES_ISBA
  CALL ABOR1_SFX('PREP_COUPLING_SURF_TRIPN: ISBA AND TRIP REQUIRE SAME RESOLUTION GRID')
ENDIF
!
IF(ZLONMIN_ISBA/=ZLONMIN_TRIP.OR.ZLONMAX_ISBA/=ZLONMAX_TRIP)THEN
  CALL ABOR1_SFX('PREP_COUPLING_SURF_TRIPN: WRONG CONFIGURATION FOR LONGITUDE')
ENDIF
!
IF(ZLATMIN_ISBA/=ZLATMIN_TRIP.OR.ZLATMAX_ISBA/=ZLATMAX_TRIP)THEN
  CALL ABOR1_SFX('PREP_COUPLING_SURF_TRIPN: WRONG CONFIGURATION FOR LATITUDE')
ENDIF        
!
IF (LHOOK) CALL DR_HOOK('PREP_COUPLING_SURF_TRIP',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE PREP_COUPLING_SURF_TRIP      
