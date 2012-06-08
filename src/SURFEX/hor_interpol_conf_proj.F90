!     #########
SUBROUTINE HOR_INTERPOL_CONF_PROJ(KLUOUT,PFIELDIN,PFIELDOUT)
!     #################################################################################
!
!
USE MODD_PREP,           ONLY : XLAT_OUT, XLON_OUT, LINTERP
USE MODD_GRID_CONF_PROJ, ONLY : XX, XY, NX, NY, XLAT0, XLON0, XLATORI, &
                                  XLONORI, XRPK, XBETA  
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE MODE_GRIDTYPE_CONF_PROJ
USE MODI_BILIN
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
INTEGER,            INTENT(IN)  :: KLUOUT    ! logical unit of output listing
REAL, DIMENSION(:,:), INTENT(IN)  :: PFIELDIN  ! field to interpolate horizontally
REAL, DIMENSION(:,:), INTENT(OUT) :: PFIELDOUT ! interpolated field
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(:),   ALLOCATABLE :: ZX       ! X coordinate
REAL, DIMENSION(:),   ALLOCATABLE :: ZY       ! Y coordinate
INTEGER                           :: INO      ! output number of points
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZFIELDIN ! input field
!
INTEGER                           :: JI       ! loop index
INTEGER                           :: JJ       ! loop index
INTEGER                           :: JL       ! loop index
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------------------
!
!*      1.    Allocations
!
IF (LHOOK) CALL DR_HOOK('HOR_INTERPOL_CONF_PROJ',0,ZHOOK_HANDLE)
INO = SIZE(XLAT_OUT)
!
ALLOCATE(ZX      (INO))
ALLOCATE(ZY      (INO))
!
!*      2.    Transformation of latitudes/longitudes into metric coordinates of output grid
!
CALL XY_CONF_PROJ(XLAT0,XLON0,XRPK,XBETA,XLATORI,XLONORI, &
                    ZX,ZY,XLAT_OUT,XLON_OUT          )  
!
!*      3.    Put input field on its squared grid
!
ALLOCATE(ZFIELDIN(NX,NY,SIZE(PFIELDIN,2)))
!
DO JJ=1,NY
  DO JI=1,NX
    ZFIELDIN(JI,JJ,:) = PFIELDIN(JI+NX*(JJ-1),:)
  END DO
END DO
!
!*      4.    Interpolation with bilinear
!
DO JL=1,SIZE(PFIELDIN,2)
  CALL BILIN(KLUOUT,XX,XY,ZFIELDIN(:,:,JL),ZX,ZY,PFIELDOUT(:,JL),LINTERP)
END DO
!
!
!
!*      5.    Deallocations
!
!
DEALLOCATE(ZX)
DEALLOCATE(ZY)
DEALLOCATE(ZFIELDIN)
IF (LHOOK) CALL DR_HOOK('HOR_INTERPOL_CONF_PROJ',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
END SUBROUTINE HOR_INTERPOL_CONF_PROJ
