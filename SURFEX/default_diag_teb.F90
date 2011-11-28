!     #########
      SUBROUTINE DEFAULT_DIAG_TEB (K2M,OSURF_BUDGET, &
        OSURF_MISC_BUDGET,OSURF_BUDGETC,OPGD,OPGD_FIX, &
        ORAD_BUDGET,PDIAG_TSTEP,ORESET_BUDGETC)  
!     #################################################################################################################
!
!!****  *DEFAULT_DIAG_TEB * - routine to set default values for the choice of diagnostics
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
!!      Modified by P. Le Moigne, 11/2004: add budget switch 
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
INTEGER,  INTENT(OUT) :: K2M                ! flag for operational 2m quantities
LOGICAL,  INTENT(OUT) :: OSURF_BUDGET       ! flag for surface budget
LOGICAL,  INTENT(OUT) :: OSURF_MISC_BUDGET  ! flag for surface miscellaneous budget
LOGICAL,  INTENT(OUT) :: OSURF_BUDGETC      ! flag for cumulated surface budget
LOGICAL,  INTENT(OUT) :: ORESET_BUDGETC     ! flag for cumulated surface budget
LOGICAL,  INTENT(OUT) :: OPGD               ! flag for PGD fields
LOGICAL,  INTENT(OUT) :: OPGD_FIX           ! flag for PGD fields
LOGICAL,  INTENT(OUT) :: ORAD_BUDGET        ! flag for radiative budget
REAL,     INTENT(OUT) :: PDIAG_TSTEP        ! time-step for writing
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DEFAULT_DIAG_TEB',0,ZHOOK_HANDLE)
K2M               = 0
OSURF_BUDGET      = .FALSE.
OSURF_MISC_BUDGET = .FALSE.
OSURF_BUDGETC     = .FALSE.
ORESET_BUDGETC    = .FALSE.
OPGD              = .FALSE.
OPGD_FIX          = .TRUE.
ORAD_BUDGET       = .FALSE.
PDIAG_TSTEP       = XUNDEF
IF (LHOOK) CALL DR_HOOK('DEFAULT_DIAG_TEB',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DEFAULT_DIAG_TEB
