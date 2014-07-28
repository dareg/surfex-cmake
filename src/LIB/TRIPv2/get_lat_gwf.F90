!###############################################
 SUBROUTINE GET_LAT_GWF(KLAT,PRES,PLAT)
!###############################################
!
!!****  *GET_LAT_GWF* - routine to get the TRIP longitude and latitude
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
INTEGER,               INTENT(IN ) :: KLAT
REAL,                  INTENT(OUT) :: PRES
REAL, DIMENSION(KLAT), INTENT(OUT) :: PLAT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*       0.2   Declarations of local variables
!              -------------------------------
!---------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GET_LAT_GWF',0,ZHOOK_HANDLE)
CALL GET_TRIP_GRID(TGRID%XTRIP_GRID,PRES=PRES,PLAT=PLAT)
IF (LHOOK) CALL DR_HOOK('GET_LAT_GWF',1,ZHOOK_HANDLE)
!    
!---------------------------------------------------------------------------
!
END SUBROUTINE GET_LAT_GWF
