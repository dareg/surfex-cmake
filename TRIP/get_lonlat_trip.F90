!###############################################
 SUBROUTINE GET_LONLAT_TRIP (TPG, &
                             KLON,KLAT,PLON,PLAT)
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
!!      B. Decharme   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2013 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_TRIP_GRID, ONLY : TRIP_GRID_t
!
USE MODE_TRIP_GRID
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
TYPE(TRIP_GRID_t), INTENT(INOUT) :: TPG
!
INTEGER,               INTENT(IN ) :: KLON
INTEGER,               INTENT(IN ) :: KLAT
REAL, DIMENSION(KLON), INTENT(OUT) :: PLON
REAL, DIMENSION(KLAT), INTENT(OUT) :: PLAT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!*       0.2   Declarations of local variables
!              -------------------------------
!---------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GET_LONLAT_TRIP',0,ZHOOK_HANDLE)
 CALL GET_TRIP_GRID(TPG%XTRIP_GRID,PLON=PLON,PLAT=PLAT)
IF (LHOOK) CALL DR_HOOK('GET_LONLAT_TRIP',1,ZHOOK_HANDLE)
!    
!---------------------------------------------------------------------------
!
END SUBROUTINE GET_LONLAT_TRIP
