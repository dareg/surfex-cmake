!     ######################
      MODULE MODD_DIAG_SEAICE_n
!     ######################
!
!!****  *MODD_DIAG_SEAICE - declaration of diagnostics for sea ice model 
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
!!      S.Sénési   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       08/2013
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

TYPE DIAG_SEAICE_t
!------------------------------------------------------------------------------
!
  LOGICAL :: LDIAG_SEAICE                    ! flag for SeaIce variables
!
!* averaged variables
!
  REAL, POINTER, DIMENSION(:)   :: XSIT  ! Sea ice thickness 
  REAL, POINTER, DIMENSION(:)   :: XSND  ! Sea ice snow depth 
  REAL, POINTER, DIMENSION(:)   :: XMLT  ! Sea mixed layer temp. seen by Gelato 
!------------------------------------------------------------------------------
!

END TYPE DIAG_SEAICE_t



CONTAINS

!





SUBROUTINE DIAG_SEAICE_INIT(YDIAG_SEAICE)
TYPE(DIAG_SEAICE_t), INTENT(INOUT) :: YDIAG_SEAICE
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_SEAICE_N:DIAG_SEAICE_INIT",0,ZHOOK_HANDLE)
  NULLIFY(YDIAG_SEAICE%XSIT)
  NULLIFY(YDIAG_SEAICE%XSND)
  NULLIFY(YDIAG_SEAICE%XMLT)
YDIAG_SEAICE%LDIAG_SEAICE=.FALSE.
IF (LHOOK) CALL DR_HOOK("MODD_DIAG_SEAICE_N:DIAG_SEAICE_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE DIAG_SEAICE_INIT


END MODULE MODD_DIAG_SEAICE_n
