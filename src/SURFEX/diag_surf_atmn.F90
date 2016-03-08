!     ########
SUBROUTINE DIAG_SURF_ATM_n (DGI, DGS, TD, DFO, DGF, DGFC, DWO, DGW, DGWC, &
                            DLO, DGL, DGLC, DUO, DGU, DGUP, DGUC, DGUPC,  U, USS, HPROGRAM)
!     #################################################################################
!
!!****  *DIAG_SURF_ATM_n * - Chooses the surface schemes for diagnostics
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
!!      Modified    01/2006 : sea flux parameterization.
!!      Modified    08/2008 : cumulated fluxes
!       B. decharme 04/2013 : Add EVAP and SUBL diag
!!------------------------------------------------------------------
!
USE MODD_SURFEX_n, ONLY : ISBA_DIAG_t, SEAFLUX_DIAG_t, TEB_DIAG_t
!
USE MODD_DIAG_n, ONLY : DIAG_t, DIAG_PATCH_t, DIAG_OPTIONS_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SSO_n, ONLY : SSO_t
!
USE MODD_SURF_CONF,      ONLY : CPROGNAME
USE MODD_DATA_COVER_PAR, ONLY : NTILESFC
!
USE MODI_DIAG_NATURE_n 
USE MODI_DIAG_SEA_n 
USE MODI_DIAG_INLAND_WATER_n 
USE MODI_DIAG_TOWN_n 
USE MODI_AVERAGE_DIAG
!
USE MODI_MINZS_VERT_SHIFT
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(ISBA_DIAG_t), INTENT(INOUT) :: DGI
TYPE(SEAFLUX_DIAG_t), INTENT(INOUT) :: DGS
TYPE(TEB_DIAG_t), INTENT(INOUT) :: TD
!
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DFO
TYPE(DIAG_t), INTENT(INOUT) :: DGF
TYPE(DIAG_t), INTENT(INOUT) :: DGFC
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DWO
TYPE(DIAG_t), INTENT(INOUT) :: DGW
TYPE(DIAG_t), INTENT(INOUT) :: DGWC
!
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DLO
TYPE(DIAG_t), INTENT(INOUT) :: DGL
TYPE(DIAG_t), INTENT(INOUT) :: DGLC
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DUO
TYPE(DIAG_t), INTENT(INOUT) :: DGU
TYPE(DIAG_PATCH_t), INTENT(INOUT) :: DGUP
TYPE(DIAG_t), INTENT(INOUT) :: DGUC
TYPE(DIAG_PATCH_t), INTENT(INOUT) :: DGUPC
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SSO_t), INTENT(INOUT) :: USS
!
CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
!
!
!*      0.2    declarations of local variables
!
INTEGER :: JTILE                        ! loop on type of surface
LOGICAL :: GNATURE, GTOWN, GWATER, GSEA ! .T. if the corresponding surface is represented
INTEGER :: JSW                          ! number of spectral whort wave bands
!
REAL, DIMENSION(SIZE(U%XSEA),NTILESFC) :: ZFRAC_TILE! fraction of each tile
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
! Preliminaries: Tile related operations
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_ATM_N',0,ZHOOK_HANDLE)
CPROGNAME = HPROGRAM
!
! FLAGS for the various surfaces:
!
GSEA      = U%NDIM_SEA    >0
GWATER    = U%NDIM_WATER  >0
GTOWN     = U%NDIM_TOWN   >0
GNATURE   = U%NDIM_NATURE >0
!
! Tile counter:
!
JTILE     = 0 
!
! Fractions for each tile:
!
ZFRAC_TILE(:,:)    = 0.0
!
! Number of spectral short wave bands for detailed radiation budget
JSW = SIZE(DGUP%AL(1)%XSWBD,2)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! SEA Tile calculations:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
! first, pack vector...then call ALMA routine
!
JTILE               = JTILE + 1
!
IF(GSEA)THEN
! 
  ZFRAC_TILE(:,JTILE) = U%XSEA(:)
!
  CALL DIAG_SEA_n(DLO, DGL, DGLC, DGS%O, DGS%D, DGS%DC, U, &
                  HPROGRAM, DGUP%AL(1), DGUPC%AL(1), U%NR_SEA)
!
ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! INLAND WATER Tile calculations:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
JTILE               = JTILE + 1
!
IF(GWATER)THEN
!
  ZFRAC_TILE(:,JTILE) = U%XWATER(:)
!
  CALL DIAG_INLAND_WATER_n(DFO, DGF, DGFC, DLO, DGL, DGLC, DWO, DGW, DGWC, &
                           U, HPROGRAM, DGUP%AL(2), DGUPC%AL(2), U%NR_WATER)
!
ENDIF 
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! NATURAL SURFACE Tile calculations:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
JTILE               = JTILE + 1
!
IF(GNATURE)THEN
!
    ZFRAC_TILE(:,JTILE) = U%XNATURE(:)
!
  CALL DIAG_NATURE_n(DGI%DE, DLO, DGL, DGLC, DGI%O, DGI%D, DGI%DC, U, &
                     HPROGRAM, DGUP%AL(3), DGUPC%AL(3), U%NR_NATURE)   
!
ENDIF 
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! URBAN Tile calculations:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
JTILE               = JTILE + 1
!
IF(GTOWN)THEN
!
    ZFRAC_TILE(:,JTILE) = U%XTOWN(:)
!
  CALL DIAG_TOWN_n(DLO, DGL, DGLC, TD%O, TD%D, U, HPROGRAM, &
                   DGUP%AL(4), DGUPC%AL(4), U%NR_TOWN)  
!
ENDIF 
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! Grid box average fluxes/properties:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
CALL AVERAGE_DIAG(ZFRAC_TILE, DUO, DGU, DGUP, DGUC, DGUPC)              
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! Quantities at 2 meters above the minimum orography of the grid mesh
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
IF (DUO%L2M_MIN_ZS) CALL GET_2M
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_ATM_N',1,ZHOOK_HANDLE)
CONTAINS
!=======================================================================================
SUBROUTINE GET_2M
!
REAL, DIMENSION(SIZE(U%XSEA)) :: ZPS         ! surface air pressure
REAL, DIMENSION(SIZE(U%XSEA)) :: ZRHOA       ! surface air density
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_ATM_n:GET_2M',0,ZHOOK_HANDLE)
!
CALL MINZS_VERT_SHIFT(DGU, U%XZS, USS%XMIN_ZS, ZPS, ZRHOA)  
DGU%XHU2M_MIN_ZS = DGU%XHU2M
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_ATM_n:GET_2M',1,ZHOOK_HANDLE)
!
END SUBROUTINE GET_2M
!
!=======================================================================================
END SUBROUTINE DIAG_SURF_ATM_n
