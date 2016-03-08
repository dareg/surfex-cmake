!     #########
      SUBROUTINE PREP_CTRL_IDEAL(DLO,KLUOUT)  
!     #################################################################################################################
!
!!****  *PREP_CTRL_IDEAL* - routine to check that diagnostics are switched off
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
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DLO
!
INTEGER,  INTENT(IN)    :: KLUOUT             ! unit number
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PREP_CTRL_IDEAL',0,ZHOOK_HANDLE)
DLO%N2M = 0
!
DLO%LSURF_BUDGET  = .FALSE.
DLO%L2M_MIN_ZS    = .FALSE.
DLO%LRAD_BUDGET   = .FALSE.
DLO%LCOEF         = .FALSE.
DLO%LSURF_VARS    = .FALSE.
DLO%LSURF_BUDGETC = .FALSE.
!
WRITE(KLUOUT,*)'IDEAL DIAGNOSTICS DESACTIVATED'
IF (LHOOK) CALL DR_HOOK('PREP_CTRL_IDEAL',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE PREP_CTRL_IDEAL
