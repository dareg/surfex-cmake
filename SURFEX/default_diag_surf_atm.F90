!     #########
      SUBROUTINE DEFAULT_DIAG_SURF_ATM(K2M,O2M_MIN_ZS,OSURF_BUDGET,ORAD_BUDGET, &
                                         OCOEF,OSURF_VARS,OFRAC,PDIAG_TSTEP,      &
                                         ODIAG_GRID,OSURF_BUDGETC,ORESET_BUDGETC, &
                                         OPROVAR_TO_DIAG, OSELECT, CSELECT        )                                         
!     ########################################################################
!
!!****  *DEFAULT_DIAG_SURF_ATM* - routine to set default values for the choice of diagnostics
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
!!      Modified    01/2006 : sea flux parameterization.
!!      B. Decharme   2008    flag for mean grid diag
!!      B. Decharme   2009    flag for cumulative budget and to write selected diags
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
LOGICAL,  INTENT(OUT) :: O2M_MIN_ZS    ! flag for 2m quantities on min.  orography
LOGICAL,  INTENT(OUT) :: OSURF_BUDGET  ! flag for surface budget
LOGICAL,  INTENT(OUT) :: ORAD_BUDGET   ! flag for radiative budget
LOGICAL,  INTENT(OUT) :: OCOEF         ! flag for transfer coefficients
LOGICAL,  INTENT(OUT) :: OSURF_VARS    ! flag for surface variables
LOGICAL,  INTENT(OUT) :: OFRAC         ! flag for fractions of tiles
LOGICAL,  INTENT(OUT) :: ODIAG_GRID    ! flag for mean grid diag
REAL,     INTENT(OUT) :: PDIAG_TSTEP   ! time-step for writing
LOGICAL,  INTENT(OUT) :: OSURF_BUDGETC      ! flag for cumulated surface budget
LOGICAL,  INTENT(OUT) :: ORESET_BUDGETC     ! flag for cumulated surface budget
LOGICAL,  INTENT(OUT) :: OPROVAR_TO_DIAG    ! switch to write (or not) prognostic variable
LOGICAL,  INTENT(OUT) :: OSELECT       ! switch to control which fields are written
CHARACTER(LEN=12), DIMENSION(200), INTENT(OUT), OPTIONAL :: CSELECT  
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DEFAULT_DIAG_SURF_ATM',0,ZHOOK_HANDLE)
K2M          = 0
O2M_MIN_ZS   = .FALSE.
!
OSURF_BUDGET = .FALSE.
ORAD_BUDGET  = .FALSE.
OCOEF        = .FALSE.
OSURF_VARS   = .FALSE.
OFRAC        = .FALSE.
!
OSURF_BUDGETC     = .FALSE.
ORESET_BUDGETC    = .FALSE.
!
ODIAG_GRID   = .TRUE.
!
PDIAG_TSTEP  = XUNDEF
!
OPROVAR_TO_DIAG    = .FALSE.
OSELECT            = .FALSE.
IF (PRESENT(CSELECT)) CSELECT(:) = '            '
IF (LHOOK) CALL DR_HOOK('DEFAULT_DIAG_SURF_ATM',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DEFAULT_DIAG_SURF_ATM
