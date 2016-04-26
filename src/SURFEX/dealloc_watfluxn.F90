!     #################################################################################
SUBROUTINE DEALLOC_WATFLUX_n (CHW, G, W)
!     #################################################################################
!
!!****  *DEALLOC_WATFLUX_n * - Deallocate all arrays
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
!
USE MODD_CH_WATFLUX_n, ONLY : CH_WATFLUX_t
USE MODD_GRID_n, ONLY : GRID_t
USE MODD_WATFLUX_n, ONLY : WATFLUX_t
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
!

!
TYPE(CH_WATFLUX_t), INTENT(INOUT) :: CHW
TYPE(GRID_t), INTENT(INOUT) :: G
TYPE(WATFLUX_t), INTENT(INOUT) :: W
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('DEALLOC_WATFLUX_N',0,ZHOOK_HANDLE)
IF (ASSOCIATED(W%LCOVER ))  DEALLOCATE(W%LCOVER )
IF (ASSOCIATED(W%XCOVER ))  DEALLOCATE(W%XCOVER )
IF (ASSOCIATED(W%XZS    ))  DEALLOCATE(W%XZS    )
IF (ASSOCIATED(W%XTS    ))  DEALLOCATE(W%XTS    )
IF (ASSOCIATED(W%XZ0    ))  DEALLOCATE(W%XZ0    )
IF (ASSOCIATED(W%XEMIS  ))  DEALLOCATE(W%XEMIS  )
!
IF (ASSOCIATED(W%XDIR_ALB))  DEALLOCATE(W%XDIR_ALB)
IF (ASSOCIATED(W%XSCA_ALB))  DEALLOCATE(W%XSCA_ALB)
!
!-------------------------------------------------------------------------------------
!
IF (ASSOCIATED(G%XGRID_PAR )) DEALLOCATE(G%XGRID_PAR )
IF (ASSOCIATED(G%XLAT      )) DEALLOCATE(G%XLAT      )
IF (ASSOCIATED(G%XLON      )) DEALLOCATE(G%XLON      )
IF (ASSOCIATED(G%XMESH_SIZE)) DEALLOCATE(G%XMESH_SIZE)
!
!-------------------------------------------------------------------------------------
!
IF(ASSOCIATED(CHW%XDEP))      DEALLOCATE(CHW%XDEP)
IF(ASSOCIATED(CHW%CCH_NAMES)) DEALLOCATE(CHW%CCH_NAMES)
IF(ASSOCIATED(CHW%SVW%CSV))       DEALLOCATE(CHW%SVW%CSV)
!
!-------------------------------------------------------------------------------------
!
IF(ASSOCIATED(W%XCPL_WATER_WIND))      DEALLOCATE(W%XCPL_WATER_WIND)
IF(ASSOCIATED(W%XCPL_WATER_FWSU))      DEALLOCATE(W%XCPL_WATER_FWSU)
IF(ASSOCIATED(W%XCPL_WATER_FWSV))      DEALLOCATE(W%XCPL_WATER_FWSV)
IF(ASSOCIATED(W%XCPL_WATER_SNET))      DEALLOCATE(W%XCPL_WATER_SNET)
IF(ASSOCIATED(W%XCPL_WATER_HEAT))      DEALLOCATE(W%XCPL_WATER_HEAT)
IF(ASSOCIATED(W%XCPL_WATER_EVAP))      DEALLOCATE(W%XCPL_WATER_EVAP)
IF(ASSOCIATED(W%XCPL_WATER_RAIN))      DEALLOCATE(W%XCPL_WATER_RAIN)
IF(ASSOCIATED(W%XCPL_WATER_SNOW))      DEALLOCATE(W%XCPL_WATER_SNOW)
IF (LHOOK) CALL DR_HOOK('DEALLOC_WATFLUX_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DEALLOC_WATFLUX_n


