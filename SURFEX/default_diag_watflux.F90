!     #########
      SUBROUTINE DEFAULT_DIAG_WATFLUX(K2M,OSURF_BUDGET,ORAD_BUDGET,PDIAG_TSTEP, &
                                        OSURF_BUDGETC,ORESET_BUDGETC              )  
!     ########################################################################
!
!!****  *DEFAULT_DIAG_WATFLUX* - routine to set default values for the choice of diagnostics
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
INTEGER,  INTENT(OUT) :: K2M           ! flag for operational 2m quantities
LOGICAL,  INTENT(OUT) :: OSURF_BUDGET  ! flag for surface budget
LOGICAL,  INTENT(OUT) :: ORAD_BUDGET   ! flag for radiative budget
REAL,     INTENT(OUT) :: PDIAG_TSTEP   ! time-step for writing
LOGICAL,  INTENT(OUT) :: OSURF_BUDGETC ! flag for cumulated surface budget
LOGICAL,  INTENT(OUT) :: ORESET_BUDGETC! flag for cumulated surface budget
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DEFAULT_DIAG_WATFLUX',0,ZHOOK_HANDLE)
K2M = 0
OSURF_BUDGET = .FALSE.
ORAD_BUDGET  = .FALSE.
PDIAG_TSTEP  = XUNDEF
OSURF_BUDGETC= .FALSE.
ORESET_BUDGETC= .FALSE.
IF (LHOOK) CALL DR_HOOK('DEFAULT_DIAG_WATFLUX',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DEFAULT_DIAG_WATFLUX
