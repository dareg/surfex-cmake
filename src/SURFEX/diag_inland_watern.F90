!     #########
SUBROUTINE DIAG_INLAND_WATER_n (DFO, DF, DFC, DLO, DL, DLC, DWO, DW, DWC, U, &
                                HPROGRAM, DUP, DUPC, KMASK       )
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
TYPE(DIAG_t), INTENT(INOUT) :: DF
TYPE(DIAG_t), INTENT(INOUT) :: DFC
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DLO
TYPE(DIAG_t), INTENT(INOUT) :: DL
TYPE(DIAG_t), INTENT(INOUT) :: DLC
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DWO
TYPE(DIAG_t), INTENT(INOUT) :: DW
TYPE(DIAG_t), INTENT(INOUT) :: DWC
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM ! program calling surf. schemes
!
TYPE(DIAG_t), INTENT(INOUT) :: DUP
TYPE(DIAG_t), INTENT(INOUT) :: DUPC
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
  CALL DIAG_EVAP(DWO, DW, DWC, HPROGRAM, DUP, DUPC, KMASK)
ELSE IF (U%CWATER=='FLAKE ') THEN
  CALL DIAG_EVAP(DFO, DF, DFC, HPROGRAM, DUP, DUPC, KMASK)       
ELSE IF (U%CWATER=='FLUX  ') THEN
  CALL DIAG_EVAP(DLO, DL, DLC, HPROGRAM, DUP, DUPC, KMASK)                     
ELSE IF (U%CWATER=='NONE  ') THEN
  CALL INIT_BUD(DWO,DUP,DUPC,XUNDEF)         
END IF
IF (LHOOK) CALL DR_HOOK('DIAG_INLAND_WATER_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_INLAND_WATER_n
