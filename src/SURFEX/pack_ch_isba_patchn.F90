!     #########
SUBROUTINE PACK_CH_ISBA_PATCH_n (CHI, PKCI, KMASK, KSIZE, KPATCH)
!##############################################
!
!
!!****  *PACK_CH_ISBA_PATCH_n * - packs chemistry variables
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
!!     V. Masson
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
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!------------------------------------------------------------------------
!
! Packed surface module variables:
!
IF (LHOOK) CALL DR_HOOK('PACK_CH_ISBA_PATCH_N',0,ZHOOK_HANDLE)
!
ALLOCATE(PKCI%XBLOCK_SIMPLE(KSIZE,2))
!
PKCI%XSOILRC_SO2 => PKCI%XBLOCK_SIMPLE(:,1)
PKCI%XSOILRC_O3 => PKCI%XBLOCK_SIMPLE(:,2)
!
ALLOCATE(PKCI%XDEP(KSIZE,CHI%SVI%NBEQ))
!
!------------------------------------------------------------------------
!
PKCI%XSOILRC_SO2   (:)    =    CHI%XSOILRC_SO2   (KMASK(:), KPATCH)
PKCI%XSOILRC_O3    (:)    =    CHI%XSOILRC_O3    (KMASK(:), KPATCH)
!
IF (LHOOK) CALL DR_HOOK('PACK_CH_ISBA_PATCH_N',1,ZHOOK_HANDLE)
!
!------------------------------------------------------------------------
!
END SUBROUTINE PACK_CH_ISBA_PATCH_n
