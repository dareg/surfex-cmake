!     #########
      SUBROUTINE DEFAULT_DIAG_ISBA(K2M,OSURF_BUDGET,O2M_MIN_ZS,ORAD_BUDGET, &
                                   OCOEF,OSURF_VARS,OSURF_EVAP_BUDGET,      &
                                   OSURF_MISC_BUDGET,OSURF_BUDGETC,         &
                                   OPATCH_BUDGET,OWOOD_SPIN,OSOILCARB_SPIN, &
                                   OPGD,ORESET_BUDGETC,PDIAG_TSTEP          )  
!     #################################################################################################################
!
!!****  *DEFAULT_DIAG_ISBA* - routine to set default values for the choice of diagnostics
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
!!      Modified by B. Decharme , 06/2009: add patch budget switch 
!!      Modified by A.L. Gibelin, 04/2009: add carbon spinup
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
LOGICAL,  INTENT(OUT) :: O2M_MIN_ZS
LOGICAL,  INTENT(OUT) :: ORAD_BUDGET        ! flag for radiative budget
LOGICAL,  INTENT(OUT) :: OCOEF
LOGICAL,  INTENT(OUT) :: OSURF_VARS
LOGICAL,  INTENT(OUT) :: OSURF_EVAP_BUDGET  ! flag for surface evaporation budget
LOGICAL,  INTENT(OUT) :: OSURF_MISC_BUDGET  ! flag for surface miscellaneous budget
LOGICAL,  INTENT(OUT) :: OSURF_BUDGETC      ! flag for cumulated surface budget
LOGICAL,  INTENT(OUT) :: OPATCH_BUDGET      ! flag for patch output
LOGICAL,  INTENT(OUT) :: OWOOD_SPIN         ! flag for wood spinup
LOGICAL,  INTENT(OUT) :: OSOILCARB_SPIN     ! flag for soil carbon spinup
LOGICAL,  INTENT(OUT) :: OPGD               ! flag for PGD fields
LOGICAL,  INTENT(OUT) :: ORESET_BUDGETC     ! flag for cumulated surface budget
REAL,     INTENT(OUT) :: PDIAG_TSTEP        ! time-step for writing
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DEFAULT_DIAG_ISBA',0,ZHOOK_HANDLE)
K2M               = 0
OSURF_BUDGET      = .FALSE.
!
O2M_MIN_ZS        = .FALSE.
ORAD_BUDGET       = .FALSE.
!
OCOEF             = .FALSE.
OSURF_VARS        = .FALSE.
!
OSURF_EVAP_BUDGET = .FALSE.
OSURF_MISC_BUDGET = .FALSE.
!
OSURF_BUDGETC     = .FALSE.
!
OPATCH_BUDGET     = .TRUE.
!
OWOOD_SPIN        = .FALSE.
OSOILCARB_SPIN    = .FALSE.
!
OPGD              = .FALSE.
ORESET_BUDGETC    = .FALSE.
!
PDIAG_TSTEP       = XUNDEF
IF (LHOOK) CALL DR_HOOK('DEFAULT_DIAG_ISBA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DEFAULT_DIAG_ISBA
