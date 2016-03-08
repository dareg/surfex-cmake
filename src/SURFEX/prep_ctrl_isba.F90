!     #########
SUBROUTINE PREP_CTRL_ISBA(DIO,OSURF_EVAP_BUDGET,OSURF_MISC_BUDGET,OSURF_MISC_DIF,KLUOUT )  
!     #################################################################################################################
!
!!****  *PREP_CTRL_ISBA* - routine to check that diagnostics are switched off
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
!!      Modified by A.L. Gibelin, 04/2009: add carbon spinup
!!
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
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DIO
!
LOGICAL,  INTENT(INOUT) :: OSURF_EVAP_BUDGET  ! flag for surface evaporation budget
LOGICAL,  INTENT(INOUT) :: OSURF_MISC_BUDGET  ! flag for surface miscellaneous budget
LOGICAL,  INTENT(INOUT) :: OSURF_MISC_DIF     ! flag for surface miscellaneous dif variables
INTEGER,  INTENT(IN)    :: KLUOUT             ! unit number
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PREP_CTRL_ISBA',0,ZHOOK_HANDLE)
DIO%N2M = 0
!
DIO%LSURF_BUDGET  = .FALSE.
DIO%L2M_MIN_ZS    = .FALSE.
DIO%LRAD_BUDGET   = .FALSE.
DIO%LCOEF         = .FALSE.
DIO%LSURF_VARS    = .FALSE.
!
DIO%LSURF_BUDGETC     = .FALSE.
!
DIO%LPATCH_BUDGET     = .FALSE.
!
OSURF_EVAP_BUDGET = .FALSE.
OSURF_MISC_BUDGET = .FALSE.
OSURF_MISC_DIF    = .FALSE.
!
WRITE(KLUOUT,*)'ISBA DIAGNOSTICS DESACTIVATED'
IF (LHOOK) CALL DR_HOOK('PREP_CTRL_ISBA',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE PREP_CTRL_ISBA
