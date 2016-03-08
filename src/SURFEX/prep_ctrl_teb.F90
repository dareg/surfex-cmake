!     #########
      SUBROUTINE PREP_CTRL_TEB (DTO, OSURF_EVAP_BUDGET,OSURF_MISC_BUDGET,OUTCI,KLUOUT)  
!     #################################################################################################################
!
!!****  *PREP_CTRL_TEB * - routine to check that diagnostics are switched off
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
!!      P. Le Moigne   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2007 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIAG_n, ONLY : DIAG_OPTIONS_t
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DTO
!
LOGICAL,  INTENT(INOUT) :: OSURF_EVAP_BUDGET  ! flag for surface evaporation budget
LOGICAL,  INTENT(INOUT) :: OSURF_MISC_BUDGET  ! flag for surface miscellaneous budget
LOGICAL,  INTENT(INOUT) :: OUTCI              ! flag for UTCI fields
INTEGER,  INTENT(IN)    :: KLUOUT             ! unit number
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PREP_CTRL_TEB',0,ZHOOK_HANDLE)
DTO%N2M = 0
!
DTO%LSURF_BUDGET  = .FALSE.
DTO%L2M_MIN_ZS    = .FALSE.
DTO%LRAD_BUDGET   = .FALSE.
DTO%LCOEF         = .FALSE.
DTO%LSURF_VARS    = .FALSE.
!
OSURF_EVAP_BUDGET = .FALSE.
OSURF_MISC_BUDGET = .FALSE.
OUTCI             = .FALSE.
!
WRITE(KLUOUT,*)'TEB  DIAGNOSTICS DESACTIVATED'
IF (LHOOK) CALL DR_HOOK('PREP_CTRL_TEB',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE PREP_CTRL_TEB 
