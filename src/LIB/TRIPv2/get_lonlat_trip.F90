!###############################################
 SUBROUTINE GET_LONLAT_TRIP(KLON,KLAT,PLON,PLAT)
!###############################################
!
!!****  *GET_LONLAT_TRIP* - routine to get the TRIP longitude and latitude
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
!!      Original    01/2013 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_TRIP_GRID
USE MODD_TRIP_GRID, ONLY : TGRID
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
INTEGER,               INTENT(IN ) :: KLON
INTEGER,               INTENT(IN ) :: KLAT
REAL, DIMENSION(KLON), INTENT(OUT) :: PLON
REAL, DIMENSION(KLAT), INTENT(OUT) :: PLAT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*       0.2   Declarations of local variables
!              -------------------------------
!---------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GET_LONLAT_TRIP',0,ZHOOK_HANDLE)
CALL GET_TRIP_GRID(TGRID%XTRIP_GRID,PLON=PLON,PLAT=PLAT)
IF (LHOOK) CALL DR_HOOK('GET_LONLAT_TRIP',1,ZHOOK_HANDLE)
!    
!---------------------------------------------------------------------------
!
END SUBROUTINE GET_LONLAT_TRIP
