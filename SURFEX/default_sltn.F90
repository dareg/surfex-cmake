!     #########
      SUBROUTINE DEFAULT_SLT_n
!     ########################################################################
!
!!****  *DEFAULT_SLT_n* - routine to set default values for the configuration for SLT
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
!!	Alf Grini CNRM
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    03/2005 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_SLT_SURF,   ONLY : CEMISPARAM_SLT
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!-------------------------------------------------------------------------------
!
! Set initial values of variables. These are modified by namelist

REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('DEFAULT_SLT_N',0,ZHOOK_HANDLE)
CEMISPARAM_SLT='Vig01'
IF (LHOOK) CALL DR_HOOK('DEFAULT_SLT_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE DEFAULT_SLT_n
