!     #########
      SUBROUTINE DEFAULT_ASSIM(OASSIM,HASSIM)
!     ########################################################################
!
!!****  *DEFAULT_ISBA* - routine to set default values for the configuration for ISBA assimilation scheme
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
!!	L. Jarlan  *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    02/2005
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
LOGICAL,           INTENT(OUT) :: OASSIM     ! assimilation or not
CHARACTER(LEN=5),  INTENT(OUT) :: HASSIM   ! type of corrections
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!                                          ! 'PLUS '
!                                          ! '2DVAR'
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('DEFAULT_ASSIM',0,ZHOOK_HANDLE)
OASSIM = .FALSE.
HASSIM = "PLUS "
IF (LHOOK) CALL DR_HOOK('DEFAULT_ASSIM',1,ZHOOK_HANDLE)
!
END SUBROUTINE DEFAULT_ASSIM
