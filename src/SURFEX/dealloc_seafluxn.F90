!     #################################################################################
SUBROUTINE DEALLOC_SEAFLUX_n
!     #################################################################################
!
!!****  *DEALLOC_SEAFLUX_n * - Deallocate all arrays
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
!!      S. Belamari 03/2014   other _SEA_ variables
!!      S. Senesi   09/2013   introduce sea-ice-cover ans sea-surface salinity
!!------------------------------------------------------------------
!
!
USE MODD_SEAFLUX_n,      ONLY : LCOVER, XCOVER, XZS, XSST, XSSS, XSIC, XZ0, XZ0H,  &
                                XSEABATHY, XEMIS, XDIR_ALB, XSCA_ALB, TGLT, XFSIC, &
                                XFSIT, XCPL_SEA_WIND, XCPL_SEA_FWSU, XCPL_SEA_FWSV,&
                                XCPL_SEA_SNET, XCPL_SEA_HEAT, XCPL_SEA_EVAP,       &
                                XCPL_SEA_RAIN, XCPL_SEA_SNOW  
USE MODD_SEAFLUX_GRID_n, ONLY : XGRID_PAR, XLAT, XLON, XMESH_SIZE
USE MODD_CH_SEAFLUX_n,   ONLY : XDEP, CCH_NAMES, CSV
!
USE MODI_GLTOOLS_DEALLOC
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
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DEALLOC_SEAFLUX_N',0,ZHOOK_HANDLE)
!
IF (ASSOCIATED(LCOVER ))   DEALLOCATE(LCOVER )
IF (ASSOCIATED(XCOVER ))   DEALLOCATE(XCOVER )
IF (ASSOCIATED(XZS    ))   DEALLOCATE(XZS    )
IF (ASSOCIATED(XSST   ))   DEALLOCATE(XSST   )
IF (ASSOCIATED(XSSS   ))   DEALLOCATE(XSSS   )
IF (ASSOCIATED(XSIC   ))   DEALLOCATE(XSIC   )
IF (ASSOCIATED(XFSIC  ))   DEALLOCATE(XFSIC  )
IF (ASSOCIATED(XFSIT  ))   DEALLOCATE(XFSIT  )
IF (ASSOCIATED(XZ0    ))   DEALLOCATE(XZ0    )
IF (ASSOCIATED(XZ0H   ))   DEALLOCATE(XZ0H   )
IF (ASSOCIATED(XSEABATHY)) DEALLOCATE(XSEABATHY)
IF (ASSOCIATED(XEMIS  ))   DEALLOCATE(XEMIS  )
IF (ASSOCIATED(XDIR_ALB))  DEALLOCATE(XDIR_ALB)
IF (ASSOCIATED(XSCA_ALB))  DEALLOCATE(XSCA_ALB)
!
!-------------------------------------------------------------------------------------
!
IF (ASSOCIATED(XGRID_PAR )) DEALLOCATE(XGRID_PAR )
IF (ASSOCIATED(XLAT      )) DEALLOCATE(XLAT      )
IF (ASSOCIATED(XLON      )) DEALLOCATE(XLON      )
IF (ASSOCIATED(XMESH_SIZE)) DEALLOCATE(XMESH_SIZE)
!
!-------------------------------------------------------------------------------------
!
IF(ASSOCIATED(XDEP))      DEALLOCATE(XDEP)
IF(ASSOCIATED(CCH_NAMES)) DEALLOCATE(CCH_NAMES)
IF(ASSOCIATED(CSV))       DEALLOCATE(CSV)
!
!-------------------------------------------------------------------------------------
!
IF(ASSOCIATED(XCPL_SEA_WIND))      DEALLOCATE(XCPL_SEA_WIND)
IF(ASSOCIATED(XCPL_SEA_FWSU))      DEALLOCATE(XCPL_SEA_FWSU)
IF(ASSOCIATED(XCPL_SEA_FWSV))      DEALLOCATE(XCPL_SEA_FWSV)
IF(ASSOCIATED(XCPL_SEA_SNET))      DEALLOCATE(XCPL_SEA_SNET)
IF(ASSOCIATED(XCPL_SEA_HEAT))      DEALLOCATE(XCPL_SEA_HEAT)
IF(ASSOCIATED(XCPL_SEA_EVAP))      DEALLOCATE(XCPL_SEA_EVAP)
IF(ASSOCIATED(XCPL_SEA_RAIN))      DEALLOCATE(XCPL_SEA_RAIN)
IF(ASSOCIATED(XCPL_SEA_SNOW))      DEALLOCATE(XCPL_SEA_SNOW)
!
!-------------------------------------------------------------------------------------
!
IF (ASSOCIATED(TGLT%bat)) CALL GLTOOLS_DEALLOC(TGLT)
!
IF (LHOOK) CALL DR_HOOK('DEALLOC_SEAFLUX_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DEALLOC_SEAFLUX_n


