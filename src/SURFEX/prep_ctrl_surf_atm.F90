!     #########
      SUBROUTINE PREP_CTRL_SURF_ATM(DUO,ONOWRITE_TEXFILE,KLUOUT)  
!     ########################################################################
!
!!****  *PREP_CTRL_SURF_ATM* - routine to check that diagnostics are switched off
!!                             
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
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DUO
LOGICAL,  INTENT(INOUT) :: ONOWRITE_TEXFILE    ! flag for surface variables
INTEGER,  INTENT(IN)    :: KLUOUT        ! unit number
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PREP_CTRL_SURF_ATM',0,ZHOOK_HANDLE)
DUO%N2M = 0
!
DUO%LSURF_BUDGET  = .FALSE.
DUO%L2M_MIN_ZS    = .FALSE.
DUO%LRAD_BUDGET   = .FALSE.
DUO%LCOEF         = .FALSE.
DUO%LSURF_VARS    = .FALSE.
!
DUO%LSURF_BUDGETC     = .FALSE.
DUO%LRESET_BUDGETC    = .FALSE.
!
ONOWRITE_TEXFILE  = .TRUE.
DUO%LSELECT           = .FALSE.
DUO%LPROVAR_TO_DIAG   = .FALSE.
!
WRITE(KLUOUT,*)'SURF_ATM DIAGNOSTICS DESACTIVATED'
IF (LHOOK) CALL DR_HOOK('PREP_CTRL_SURF_ATM',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE PREP_CTRL_SURF_ATM
