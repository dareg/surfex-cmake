!     #########
      SUBROUTINE PREP_CTRL_WATFLUX(DWO,KLUOUT)  
!     #################################################################################################################
!
!!****  *PREP_CTRL_WATFLUX* - routine to check that diagnostics are switched off
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
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DWO
INTEGER,  INTENT(IN)    :: KLUOUT             ! unit number
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PREP_CTRL_WATFLUX',0,ZHOOK_HANDLE)
DWO%N2M = 0
!
DWO%LSURF_BUDGET  = .FALSE.
DWO%L2M_MIN_ZS    = .FALSE.
DWO%LRAD_BUDGET   = .FALSE.
DWO%LCOEF         = .FALSE.
DWO%LSURF_VARS    = .FALSE.
DWO%LSURF_BUDGETC = .FALSE.
!
WRITE(KLUOUT,*)'WATFLUX DIAGNOSTICS DESACTIVATED'
IF (LHOOK) CALL DR_HOOK('PREP_CTRL_WATFLUX',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE PREP_CTRL_WATFLUX
