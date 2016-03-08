!     #########
SUBROUTINE DIAG_SEA_n (DLO, DGL, DGLC, DSO, DGS, DGSC, U, HPROGRAM, DGUP, DGUPC, KMASK )
!     #####################################################################
!
!!****  *DIAG_SEA_n * - Chooses the surface schemes for sea diagnostics
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
!!      Modified    08/2009 : new diag
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
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DLO
TYPE(DIAG_t), INTENT(INOUT) :: DGL
TYPE(DIAG_t), INTENT(INOUT) :: DGLC
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DSO
TYPE(DIAG_t), INTENT(INOUT) :: DGS
TYPE(DIAG_t), INTENT(INOUT) :: DGSC
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM ! program calling surf. schemes
!
TYPE(DIAG_t), INTENT(INOUT) :: DGUP
TYPE(DIAG_t), INTENT(INOUT) :: DGUPC
!
INTEGER, DIMENSION(:), INTENT(IN) :: KMASK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DIAG_SEA_N',0,ZHOOK_HANDLE)
IF (U%CSEA=='SEAFLX') THEN
  CALL DIAG_EVAP(DSO, DGS, DGSC, HPROGRAM, DGUP, DGUPC, KMASK)
ELSEIF (U%CSEA=='FLUX') THEN
  CALL DIAG_EVAP(DLO, DGL, DGLC, HPROGRAM, DGUP, DGUPC, KMASK)            
ELSE IF (U%CSEA=='NONE  ') THEN
  CALL INIT_BUD(DSO,DGUP,DGUPC,XUNDEF)
END IF
IF (LHOOK) CALL DR_HOOK('DIAG_SEA_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_SEA_n
