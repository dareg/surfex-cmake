!###############
MODULE MODD_TRIP
!###############
!
!!****  *MODD_TRIP - declaration of surface variable for TRIP RRM
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	B. Decharme   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       21/05/08
!
!*       0.   DECLARATIONS
!             ------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE TRIP_t
!
!-------------------------------------------------------------------------------
!
! Input river geometry Parameters :
!
REAL, POINTER, DIMENSION(:,:) :: XTAUG          ! ground water transfer time   [s]
REAL, POINTER, DIMENSION(:,:) :: XSLOPEBED      ! river bed slopes             [m/m]
REAL, POINTER, DIMENSION(:,:) :: XWIDTH         ! river widths                 [m]
REAL, POINTER, DIMENSION(:,:) :: XN             ! Manning roughness coeficient [-] (0.03 to 0.065)
REAL, POINTER, DIMENSION(:,:) :: XN_FLOOD       ! Manning coeficient over floodplains  [-] (currently 0.1)
REAL, POINTER, DIMENSION(:,:) :: XHC_BED        ! River bed depth              [m]
REAL, POINTER, DIMENSION(:,:) :: XWEFF          ! Porosité efficace
REAL, POINTER, DIMENSION(:,:) :: XTRANS         ! Transmissivité
REAL, POINTER, DIMENSION(:,:) :: XNUM_AQUI      ! Numéro aquifère
REAL, POINTER, DIMENSION(:,:) :: XELEV          ! Elevation                    [m]
REAL, POINTER, DIMENSION(:,:) :: XTOPO_RIV      ! River elevation              [m]
!
!-------------------------------------------------------------------------------
!
! Time varing variables :
!
REAL, POINTER, DIMENSION(:,:) :: XSURF_STO        ! river channel storage        [kg]
REAL, POINTER, DIMENSION(:,:) :: XGROUND_STO      ! groundwater storage          [kg]
REAL, POINTER, DIMENSION(:,:) :: XFLOOD_STO       ! Floodplain water storage     [kg]
REAL, POINTER, DIMENSION(:,:) :: XHGROUND         ! Groudwater height            [m]
REAL, POINTER, DIMENSION(:,:) :: XHFLOOD          ! Floodplain water depth       [m]
REAL, POINTER, DIMENSION(:,:) :: XFFLOOD          ! Floodplain grid-cell fraction [-]
REAL, POINTER, DIMENSION(:,:) :: XWFLOOD          ! Floodplain width             [m]
REAL, POINTER, DIMENSION(:,:) :: XFLOOD_LEN       ! Floodplain lenght            [m]
!
!-------------------------------------------------------------------------------
!
! Floodplain fonctions :
!
REAL, POINTER, DIMENSION(:,:,:) :: XTAB_F         ! Flood fraction array
REAL, POINTER, DIMENSION(:,:,:) :: XTAB_H         ! Topo height array
REAL, POINTER, DIMENSION(:,:,:) :: XTAB_VF        ! Flood volume array
!
!-------------------------------------------------------------------------------
!
! Groundwater fonctions :
!
REAL, POINTER, DIMENSION(:,:,:) :: XTABGW_F       ! Groundwater fraction array
REAL, POINTER, DIMENSION(:,:,:) :: XTABGW_H       ! Topo height array
!
!-------------------------------------------------------------------------------
!
! Coupling variable with other models :
!
REAL, POINTER, DIMENSION(:,:) :: XCPL_FWTD         ! grid-cell fraction of water table to rise
REAL, POINTER, DIMENSION(:,:) :: XCPL_WTD          ! Water table depth            [m]
!
REAL, POINTER, DIMENSION(:,:) :: XCPL_RIVDIS       ! River discharges  [kg/m2]
REAL, POINTER, DIMENSION(:,:) :: XCPL_FFLOOD       ! Flood fraction    [-]
REAL, POINTER, DIMENSION(:,:) :: XCPL_PIFLOOD      ! Floodplains potential infiltration  [kg/m2]
!
REAL, POINTER, DIMENSION(:,:) :: XCPL_CALVGRE      ! Calving flux over greenland
REAL, POINTER, DIMENSION(:,:) :: XCPL_CALVANT      ! Calving flux over antarctica
!
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!
END TYPE TRIP_t
!
TYPE(TRIP_t), ALLOCATABLE, TARGET, SAVE :: TRIP_MODEL(:)
!
TYPE(TRIP_t), POINTER  :: TTRIP=>NULL()
!$OMP THREADPRIVATE(TTRIP)
!
CONTAINS
!
SUBROUTINE TRIP_GOTO_MODEL(KTO)
INTEGER, INTENT(IN) :: KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_TRIP:TRIP_GOTO_MODEL',0,ZHOOK_HANDLE)
!
TTRIP=>TRIP_MODEL(KTO)
!
IF (LHOOK) CALL DR_HOOK('MODD_TRIP:TRIP_GOTO_MODEL',1,ZHOOK_HANDLE)
!
END SUBROUTINE TRIP_GOTO_MODEL

SUBROUTINE TRIP_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TRIP:TRIP_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(TRIP_MODEL(KMODEL))
DO J=1,KMODEL
!
  NULLIFY(TRIP_MODEL(J)%XTAUG)
  NULLIFY(TRIP_MODEL(J)%XSLOPEBED)
  NULLIFY(TRIP_MODEL(J)%XWIDTH)
  NULLIFY(TRIP_MODEL(J)%XN)
  NULLIFY(TRIP_MODEL(J)%XN_FLOOD)
  NULLIFY(TRIP_MODEL(J)%XHC_BED)
  NULLIFY(TRIP_MODEL(J)%XSURF_STO)
  NULLIFY(TRIP_MODEL(J)%XGROUND_STO)
  NULLIFY(TRIP_MODEL(J)%XFLOOD_STO)
  NULLIFY(TRIP_MODEL(J)%XHFLOOD)
  NULLIFY(TRIP_MODEL(J)%XFFLOOD)
  NULLIFY(TRIP_MODEL(J)%XWFLOOD)
  NULLIFY(TRIP_MODEL(J)%XFLOOD_LEN)
  NULLIFY(TRIP_MODEL(J)%XWEFF)
  NULLIFY(TRIP_MODEL(J)%XTRANS)
  NULLIFY(TRIP_MODEL(J)%XNUM_AQUI)
  NULLIFY(TRIP_MODEL(J)%XELEV)
  NULLIFY(TRIP_MODEL(J)%XTOPO_RIV)
  NULLIFY(TRIP_MODEL(J)%XTABGW_H)
  NULLIFY(TRIP_MODEL(J)%XHGROUND)
  NULLIFY(TRIP_MODEL(J)%XTAB_F)
  NULLIFY(TRIP_MODEL(J)%XTAB_H)
  NULLIFY(TRIP_MODEL(J)%XTAB_VF)
  NULLIFY(TRIP_MODEL(J)%XTABGW_F)
  NULLIFY(TRIP_MODEL(J)%XTABGW_H)
!
  NULLIFY(TRIP_MODEL(J)%XCPL_WTD)
  NULLIFY(TRIP_MODEL(J)%XCPL_FWTD)
  NULLIFY(TRIP_MODEL(J)%XCPL_RIVDIS)
  NULLIFY(TRIP_MODEL(J)%XCPL_FFLOOD)
  NULLIFY(TRIP_MODEL(J)%XCPL_PIFLOOD)
  NULLIFY(TRIP_MODEL(J)%XCPL_CALVGRE)
  NULLIFY(TRIP_MODEL(J)%XCPL_CALVANT)
!  
ENDDO
IF (LHOOK) CALL DR_HOOK("MODD_TRIP:TRIP_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE TRIP_ALLOC

SUBROUTINE TRIP_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TRIP:TRIP_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(TRIP_MODEL)) DEALLOCATE(TRIP_MODEL)
IF (LHOOK) CALL DR_HOOK("MODD_TRIP:TRIP_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE TRIP_DEALLO
!
END MODULE MODD_TRIP
