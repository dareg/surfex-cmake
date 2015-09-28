!     ########################################
      SUBROUTINE GET_OFFLINE_CONF(PTSTEP_OUTPUT)
!     ########################################
!
!
!!****  *GET_OFFLINE_CONF* - routine to get some ISBA fields
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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2008
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODN_IO_OFFLINE,     ONLY : XTSTEP_OUTPUT
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
REAL, INTENT(OUT) :: PTSTEP_OUTPUT ! time step of output time series
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GET_OFFLINE_CONF',0,ZHOOK_HANDLE)
PTSTEP_OUTPUT = XTSTEP_OUTPUT
IF (LHOOK) CALL DR_HOOK('GET_OFFLINE_CONF',1,ZHOOK_HANDLE)
!
!==============================================================================
!
END SUBROUTINE GET_OFFLINE_CONF
