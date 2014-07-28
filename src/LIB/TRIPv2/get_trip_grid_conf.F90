!     #########
      SUBROUTINE GET_TRIP_GRID_CONF(PLONMIN,PLONMAX,PLATMIN,PLATMAX,PRES,KLON,KLAT)
!     #########################################
!
!!****  *GET_TRIP_GRID_CONF* - routine to get the TRIP grid configuration
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
REAL,               INTENT(OUT) :: PLONMIN
REAL,               INTENT(OUT) :: PLONMAX
REAL,               INTENT(OUT) :: PLATMIN
REAL,               INTENT(OUT) :: PLATMAX
REAL,               INTENT(OUT) :: PRES
INTEGER,            INTENT(OUT) :: KLON
INTEGER,            INTENT(OUT) :: KLAT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*       0.2   Declarations of local variables
!              -------------------------------
!---------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GET_TRIP_GRID_CONF',0,ZHOOK_HANDLE)
!
CALL GET_TRIP_GRID(TGRID%XTRIP_GRID,PLONMIN,PLONMAX,PLATMIN,PLATMAX,PRES,KLON,KLAT)
!
IF (LHOOK) CALL DR_HOOK('GET_TRIP_GRID_CONF',1,ZHOOK_HANDLE)
!    
!---------------------------------------------------------------------------
!
END SUBROUTINE GET_TRIP_GRID_CONF
