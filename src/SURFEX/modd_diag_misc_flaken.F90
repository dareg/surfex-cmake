!     ############################
      MODULE MODD_DIAG_MISC_FLAKE_n
!     ############################
!
!!****  *MODD_DIAG_MISC_FLAKE - declaration of diagnostic variables for FLAKE scheme
!!
!!    PURPOSE
!!    -------
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
!!      P. Le Moigne   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       07/10/04
!
!
!*       0.   DECLARATIONS
!             ------------
!
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE DIAG_MISC_FLAKE_t
!------------------------------------------------------------------------------
!
  LOGICAL :: LWATER_PROFILE   ! flag for miscellaneous terms of FLake scheme
!
!* miscellaneous variables
!
  REAL, POINTER, DIMENSION(:) :: XZWAT_PROFILE ! depth of output levels (m) in namelist
  REAL, POINTER, DIMENSION(:)   :: XZW_PROFILE ! depth of output levels (m)
  REAL, POINTER, DIMENSION(:,:) :: XTW_PROFILE ! Water temperature in output levels (K)
!
!
!------------------------------------------------------------------------------
!

END TYPE DIAG_MISC_FLAKE_t



CONTAINS

!




SUBROUTINE DIAG_MISC_FLAKE_INIT(YDIAG_MISC_FLAKE)
TYPE(DIAG_MISC_FLAKE_t), INTENT(INOUT) :: YDIAG_MISC_FLAKE
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_MISC_FLAKE_N:DIAG_MISC_FLAKE_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YDIAG_MISC_FLAKE%XZWAT_PROFILE)
  NULLIFY(YDIAG_MISC_FLAKE%XZW_PROFILE)
  NULLIFY(YDIAG_MISC_FLAKE%XTW_PROFILE)
YDIAG_MISC_FLAKE%LWATER_PROFILE=.FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_MISC_FLAKE_N:DIAG_MISC_FLAKE_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_MISC_FLAKE_INIT


END MODULE MODD_DIAG_MISC_FLAKE_n
