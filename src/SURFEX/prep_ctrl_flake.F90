!     #########
SUBROUTINE PREP_CTRL_FLAKE(DFO, KLUOUT,OWATER_PROFILE)  
!     #################################################################################################################
!
!!****  *PREP_CTRL_FLAKE* - routine to check that diagnostics are switched off
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
!
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DFO
!
INTEGER,  INTENT(IN)    :: KLUOUT             ! unit number
LOGICAL,  INTENT(INOUT) :: OWATER_PROFILE     !
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PREP_CTRL_FLAKE',0,ZHOOK_HANDLE)
!
DFO%N2M = 0
!
DFO%LSURF_BUDGET  = .FALSE.
DFO%L2M_MIN_ZS    = .FALSE.
DFO%LRAD_BUDGET   = .FALSE.
DFO%LCOEF         = .FALSE.
DFO%LSURF_VARS    = .FALSE.
DFO%LSURF_BUDGETC = .FALSE.
!
OWATER_PROFILE = .FALSE.
!
WRITE(KLUOUT,*)'FLAKE DIAGNOSTICS DESACTIVATED'
!
IF (LHOOK) CALL DR_HOOK('PREP_CTRL_FLAKE',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE PREP_CTRL_FLAKE
