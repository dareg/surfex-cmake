!     #########
      SUBROUTINE GET_GRID_CONF_ISBA_n(PLONMIN,PLONMAX,PLATMIN,PLATMAX,KX,KY,KL)
!     #########################################
!
!!****  *GET_GRID_CONF_ISBA_n* - routine to get the ISBA grid configuration
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
!!      Original    01/2004 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_GRIDTYPE_LONLAT_REG
USE MODD_ISBA_GRID_n,    ONLY : CGRID, XGRID_PAR 
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
INTEGER,            INTENT(OUT) :: KX
INTEGER,            INTENT(OUT) :: KY
INTEGER,            INTENT(OUT) :: KL
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
 CHARACTER(LEN=100) :: YCOMMENT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!---------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GET_GRID_CONF_ISBA_N',0,ZHOOK_HANDLE)
SELECT CASE (CGRID)
!     
!  CASE("CONF PROJ ")
!
!  CASE("CARTESIAN ")

  CASE("LONLAT REG")
    CALL GET_GRIDTYPE_LONLAT_REG(XGRID_PAR,PLONMIN,PLONMAX,PLATMIN,PLATMAX,KX,KY,KL)
!
!  CASE("GAUSS     ")
!
!  CASE("NONE      ", "IGN       ")

END SELECT
IF (LHOOK) CALL DR_HOOK('GET_GRID_CONF_ISBA_N',1,ZHOOK_HANDLE)
!
!---------------------------------------------------------------------------
!
END SUBROUTINE GET_GRID_CONF_ISBA_n
