!     #########
      SUBROUTINE DEFAULT_DST_n
!     ########################################################################
!
!!****  *DEFAULT_DST_n* - routine to set default values for the configuration for DST
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
USE MODD_DST_SURF,   ONLY : CEMISPARAM, CVERMOD, XFLX_MSS_FDG_FCT
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

IF (LHOOK) CALL DR_HOOK('DEFAULT_DST_N',0,ZHOOK_HANDLE)
XFLX_MSS_FDG_FCT = 12.0e-4 
CEMISPARAM='AMMA'
CVERMOD='NONE  '
IF (LHOOK) CALL DR_HOOK('DEFAULT_DST_N',1,ZHOOK_HANDLE)

!-------------------------------------------------------------------------------
!
END SUBROUTINE DEFAULT_DST_n
