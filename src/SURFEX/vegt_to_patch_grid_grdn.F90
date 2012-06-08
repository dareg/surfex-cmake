!     #########
SUBROUTINE VEGT_TO_PATCH_GRID_GRDN(PFIELDOUT,PW)
!        ################################################
!!
!!****  *VEGTYPE_GRID_TO_PATCH_GRID* averages fields from all (12) vegtypes 
!!                                   on only a few patches
!!    PURPOSE
!!    -------
!
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
!!
!!	V. Masson          * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!
!-------------------------------------------------------------------------------

!
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_TEB_GARDEN_n,   ONLY : NPATCH, XVEGTYPE_PATCH, XPATCH
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE
!
USE MODI_VEGTYPE_TO_PATCH
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:,:), INTENT(IN)                          :: PFIELDOUT
REAL, DIMENSION(:,:,:), INTENT(OUT)                         :: PW
!
!
!*      0.2    declarations of local variables
!
INTEGER                       :: JPATCH    ! loop on patches
INTEGER                       :: JVEGTYPE  ! loop on vegtypes
INTEGER                       :: JLAYER    ! loop on layers
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!---------------------------------------------------------------------------
!
!* averages from vegtypes to chosen number of patches
IF (LHOOK) CALL DR_HOOK('VEGT_TO_PATCH_GRID_GRDN',0,ZHOOK_HANDLE)
PW(:,:,:) = 0.
DO JVEGTYPE=1,NVEGTYPE
  JPATCH = VEGTYPE_TO_PATCH(JVEGTYPE,NPATCH)
  DO JLAYER=1,SIZE(PW,2)
    PW(:,JLAYER,JPATCH) = PW(:,JLAYER,JPATCH)                                            &
                          + XVEGTYPE_PATCH(:,JVEGTYPE,JPATCH) * PFIELDOUT(:,JLAYER,JVEGTYPE)  
  END DO
END DO
!
!* insures undefined value when patch is not present
DO JPATCH=1,NPATCH
  DO JLAYER=1,SIZE(PW,2)
    WHERE(XPATCH(:,JPATCH)==0.) PW(:,JLAYER,JPATCH) = XUNDEF
  END DO
END DO
!
WHERE( ABS(PW-XUNDEF)/XUNDEF < 1.E-6 ) PW = XUNDEF
IF (LHOOK) CALL DR_HOOK('VEGT_TO_PATCH_GRID_GRDN',1,ZHOOK_HANDLE)
!
!---------------------------------------------------------------------------
!
END SUBROUTINE VEGT_TO_PATCH_GRID_GRDN
