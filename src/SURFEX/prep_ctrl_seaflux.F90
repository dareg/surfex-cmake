!     #########
      SUBROUTINE PREP_CTRL_SEAFLUX(DSO,ODIAG_OCEAN,ODIAG_MISC_SEAICE,KLUOUT)  
!     #################################################################################################################
!
!!****  *PREP_CTRL_SEAFLUX* - routine to check that diagnostics are switched off
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
!!      Modified    09/2013 : S. Senesi : manage ODIAG_SEAICE
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
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DSO
!
LOGICAL, INTENT(INOUT) :: ODIAG_OCEAN
LOGICAL,  INTENT(INOUT) :: ODIAG_MISC_SEAICE       ! flag for seaice variables
INTEGER,  INTENT(IN)    :: KLUOUT             ! unit number
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PREP_CTRL_SEAFLUX',0,ZHOOK_HANDLE)
DSO%N2M = 0
!
DSO%LSURF_BUDGET  = .FALSE.
DSO%L2M_MIN_ZS    = .FALSE.
DSO%LRAD_BUDGET   = .FALSE.
DSO%LCOEF         = .FALSE.
DSO%LSURF_VARS    = .FALSE.
!
ODIAG_OCEAN   = .FALSE.
!
ODIAG_MISC_SEAICE  = .FALSE.
!
DSO%LSURF_BUDGETC = .FALSE.
!
WRITE(KLUOUT,*)'SEAFLUX DIAGNOSTICS DESACTIVATED'
IF (LHOOK) CALL DR_HOOK('PREP_CTRL_SEAFLUX',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE PREP_CTRL_SEAFLUX
