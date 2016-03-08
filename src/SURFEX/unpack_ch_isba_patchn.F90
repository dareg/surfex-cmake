!     #########
SUBROUTINE UNPACK_CH_ISBA_PATCH_n (CHI, PKCI, KMASK, KSIZE, KPATCH)
!##############################################
!
!!****  *UNPACK_CH_ISBA_PATCH_n* - unpacks ISBA prognostic variables
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     A. Boone
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!------------------------------------------------------------------
!
!
!
USE MODD_CH_ISBA_n, ONLY : CH_ISBA_t
USE MODD_PACK_CH_ISBA, ONLY : PACK_CH_ISBA_t
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!
TYPE(CH_ISBA_t), INTENT(INOUT) :: CHI
TYPE(PACK_CH_ISBA_t), INTENT(INOUT) :: PKCI
!
INTEGER, INTENT(IN)               :: KSIZE, KPATCH
!
INTEGER, DIMENSION(:), INTENT(IN) :: KMASK
!
INTEGER :: JSV
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
! Only save values for patches which are in use:
!
IF (LHOOK) CALL DR_HOOK('UNPACK_CH_ISBA_PATCH_N',0,ZHOOK_HANDLE)
CHI%XDEP(:,:,KPATCH) = XUNDEF
!
DO JSV=1,SIZE(CHI%XDEP,2)
  CHI%XDEP(:,JSV,KPATCH) = PKCI%XDEP(KMASK(:),JSV) 
END DO
!
PKCI%XSOILRC_SO2 => NULL()
PKCI%XSOILRC_O3  => NULL()
!
DEALLOCATE(PKCI%XBLOCK_SIMPLE)
DEALLOCATE(PKCI%XDEP)
!
IF (LHOOK) CALL DR_HOOK('UNPACK_CH_ISBA_PATCH_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE UNPACK_CH_ISBA_PATCH_n
