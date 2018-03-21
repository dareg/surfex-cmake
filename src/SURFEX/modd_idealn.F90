!     #################
      MODULE MODD_IDEAL_n
!     #################
!
!!****  *MODD_IDEAL_n - declaration of surface parameters for an inland water surface
!!
!!    PURPOSE
!!    -------
!     Declaration of surface parameters
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2004
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TYPE_DATE_SURF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE IDEAL_t
!
! Time-step:
!
  REAL                                  :: XTSTEP  ! time step
!
  REAL                                  :: XOUT_TSTEP  ! output writing time step
!
!
!
END TYPE IDEAL_t



CONTAINS

!




SUBROUTINE IDEAL_INIT(YIDEAL)
TYPE(IDEAL_t), INTENT(INOUT) :: YIDEAL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_IDEAL_N:IDEAL_INIT",0,ZHOOK_HANDLE)
YIDEAL%XTSTEP=0.
YIDEAL%XOUT_TSTEP=0.
IF (LHOOK) CALL DR_HOOK("MODD_IDEAL_N:IDEAL_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE IDEAL_INIT


END MODULE MODD_IDEAL_n
