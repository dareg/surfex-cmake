!     ########
SUBROUTINE DIAG_SURF_ATM_n (ID, DS, TD, DFO, DF, DFC, DWO, DW, DWC, &
                            DLO, DL, DLC, DUO, DU, DUP, DUC, DUPC,  U, USS, HPROGRAM)
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
TYPE(ISBA_DIAG_t), INTENT(INOUT) :: ID
TYPE(SEAFLUX_DIAG_t), INTENT(INOUT) :: DS
TYPE(TEB_DIAG_t), INTENT(INOUT) :: TD
!
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DFO
TYPE(DIAG_t), INTENT(INOUT) :: DF
TYPE(DIAG_t), INTENT(INOUT) :: DFC
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DWO
TYPE(DIAG_t), INTENT(INOUT) :: DW
TYPE(DIAG_t), INTENT(INOUT) :: DWC
!
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DLO
TYPE(DIAG_t), INTENT(INOUT) :: DL
TYPE(DIAG_t), INTENT(INOUT) :: DLC
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DUO
TYPE(DIAG_t), INTENT(INOUT) :: DU
TYPE(DIAG_PATCH_t), INTENT(INOUT) :: DUP
TYPE(DIAG_t), INTENT(INOUT) :: DUC
TYPE(DIAG_PATCH_t), INTENT(INOUT) :: DUPC
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
JSW = SIZE(DUP%AL(1)%XSWBD,2)
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
  CALL DIAG_SEA_n(DLO, DL, DLC, DS%O, DS%D, DS%DC, U, &
                  HPROGRAM, DUP%AL(1), DUPC%AL(1), U%NR_SEA)
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
  CALL DIAG_INLAND_WATER_n(DFO, DF, DFC, DLO, DL, DLC, DWO, DW, DWC, &
                           U, HPROGRAM, DUP%AL(2), DUPC%AL(2), U%NR_WATER)
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
  CALL DIAG_NATURE_n(ID%DE, DLO, DL, DLC, ID%O, ID%D, ID%DC, U, &
                     HPROGRAM, DUP%AL(3), DUPC%AL(3), U%NR_NATURE)   
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
  CALL DIAG_TOWN_n(DLO, DL, DLC, TD%O, TD%D, U, HPROGRAM, &
                   DUP%AL(4), DUPC%AL(4), U%NR_TOWN)  
!
ENDIF 
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! Grid box average fluxes/properties:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
CALL AVERAGE_DIAG(ZFRAC_TILE, DUO, DU, DUP, DUC, DUPC)              
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
CALL MINZS_VERT_SHIFT(DU, U%XZS, USS%XMIN_ZS, ZPS, ZRHOA)  
DU%XHU2M_MIN_ZS = DU%XHU2M
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_ATM_n:GET_2M',1,ZHOOK_HANDLE)
!
END SUBROUTINE GET_2M
!
!=======================================================================================
END SUBROUTINE DIAG_SURF_ATM_n
