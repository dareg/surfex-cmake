!     #########
      SUBROUTINE DEFAULT_FLAKE(PTSTEP,POUT_TSTEP,OSEDIMENTS,HFLAKE_SNOW, &
        OWATFLX)  
!     ########################################################################
!
!!****  *DEFAULT_FLAKE* - routine to set default values for the configuration for FLAKE scheme
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
USE MODD_SURF_PAR,   ONLY : XUNDEF
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
!
REAL,              INTENT(OUT) :: PTSTEP     ! time step for run
REAL,              INTENT(OUT) :: POUT_TSTEP ! time step for writing
!
LOGICAL,          INTENT(OUT) :: OSEDIMENTS 
LOGICAL,          INTENT(OUT) :: OWATFLX 
CHARACTER(LEN=6), INTENT(OUT) :: HFLAKE_SNOW 
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DEFAULT_FLAKE',0,ZHOOK_HANDLE)
PTSTEP     = XUNDEF
POUT_TSTEP = XUNDEF
!
OSEDIMENTS  = .TRUE.
HFLAKE_SNOW = "NON"
OWATFLX = .FALSE.
IF (LHOOK) CALL DR_HOOK('DEFAULT_FLAKE',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE DEFAULT_FLAKE
