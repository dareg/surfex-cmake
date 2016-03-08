!     #########
SUBROUTINE DIAG_INLAND_WATER_n (DFO, DGF, DGFC, DLO, DGL, DGLC, DWO, DGW, DGWC, U, &
                                HPROGRAM, DGUP, DGUPC, KMASK       )
!     ###############################################################################
!
!!****  *DIAG_INLAND_WATER_n * - Chooses the surface schemes for lakes diagnostics
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
!!      Modified    08/2009 : cumulated diag & t2m min/max
!!       V.Masson   10/2013 Adds min and max 2m parameters
!       B. decharme 04/2013 : Add EVAP and SUBL diag
!!------------------------------------------------------------------
!
USE MODE_DIAG
!
USE MODD_DIAG_n, ONLY : DIAG_t, DIAG_OPTIONS_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DFO
TYPE(DIAG_t), INTENT(INOUT) :: DGF
TYPE(DIAG_t), INTENT(INOUT) :: DGFC
TYPE(DIAG_t), INTENT(INOUT) :: DGL
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DLO
TYPE(DIAG_t), INTENT(INOUT) :: DGLC
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DWO
TYPE(DIAG_t), INTENT(INOUT) :: DGW
TYPE(DIAG_t), INTENT(INOUT) :: DGWC
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM ! program calling surf. schemes
!
TYPE(DIAG_t), INTENT(INOUT) :: DGUP
TYPE(DIAG_t), INTENT(INOUT) :: DGUPC
!
INTEGER, DIMENSION(:), INTENT(IN) :: KMASK
!
!*      0.2    declarations of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DIAG_INLAND_WATER_N',0,ZHOOK_HANDLE)
IF (U%CWATER=='WATFLX') THEN
  CALL DIAG_EVAP(DWO, DGW, DGWC, HPROGRAM, DGUP, DGUPC, KMASK)
ELSE IF (U%CWATER=='FLAKE ') THEN
  CALL DIAG_EVAP(DFO, DGF, DGFC, HPROGRAM, DGUP, DGUPC, KMASK)       
ELSE IF (U%CWATER=='FLUX  ') THEN
  CALL DIAG_EVAP(DLO, DGL, DGLC, HPROGRAM, DGUP, DGUPC, KMASK)                     
ELSE IF (U%CWATER=='NONE  ') THEN
  CALL INIT_BUD(DWO,DGUP,DGUPC,XUNDEF)         
END IF
IF (LHOOK) CALL DR_HOOK('DIAG_INLAND_WATER_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_INLAND_WATER_n
