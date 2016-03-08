!     #########
SUBROUTINE DIAG_TOWN_n (DLO, DGL, DGLC, DTO, DGT, U, HPROGRAM, DGUP, DGUPC, KMASK )
!     ######################################################################
!
!!****  *DIAG_TOWN_n * - Chooses the surface schemes for town diagnostics
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
!!      Modified    09/2012 : new PLEI diag required by atmospheric model
!       B. decharme 04/2013 : Add EVAP and SUBL diag
!!------------------------------------------------------------------
!
USE MODE_DIAG
!
USE MODD_DIAG_n, ONLY : DIAG_t, DIAG_OPTIONS_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_CSTS,       ONLY : XTT, XLSTT, XLVTT
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DLO
TYPE(DIAG_t), INTENT(INOUT) :: DGL
TYPE(DIAG_t), INTENT(INOUT) :: DGLC
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DTO
TYPE(DIAG_t), INTENT(INOUT) :: DGT
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
REAL, DIMENSION(SIZE(DGUP%XRN)) :: ZDELTA
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DIAG_TOWN_N',0,ZHOOK_HANDLE)
IF (U%CTOWN=='TEB   ') THEN

  CALL DIAG(DTO, DGT, HPROGRAM, DGUP, KMASK)
!
!!!!! important, diagd should be computed in teb !!!!!!
!
! diag not yet inplemeted for TEB (these diag are required for the climate model)
!
! Ok with atmospheric model but LEI (latent heat of sublimation w/m2), EVAP (total evapotranspiration kg/m2/s),
! and SUBL (sublimation kg/m2/s) must by implemented in TEB as well as theirs cumulative values
! Not good if LCPL_ARP = TRUE in ISBA (ALARO)
!
  IF (SIZE(DGUP%XLEI)>0) THEN
    DGUP%XLEI (:) = XUNDEF
    DGUP%XEVAP(:) = XUNDEF
    DGUP%XSUBL(:) = XUNDEF
    WHERE(DGUP%XLE(:)/=XUNDEF)
      ZDELTA(:) = MAX(0.0,SIGN(1.0,XTT-DGUP%XTS(:)))
      DGUP%XEVAP (:) = (DGUP%XLE(:) * ZDELTA(:))/XLSTT + (DGUP%XLE(:) * (1.0-ZDELTA(:)))/XLVTT
      DGUP%XLEI  (:) = DGUP%XLE(:) * ZDELTA(:)
      DGUP%XSUBL (:) = DGUP%XLEI(:)/XLSTT
    ENDWHERE
  ENDIF
!
  IF (DTO%LSURF_BUDGETC) CALL INIT_SURF_BUD(DGUPC,XUNDEF)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      
ELSE IF (U%CTOWN=='FLUX  ') THEN
  CALL DIAG_EVAP(DLO, DGL, DGLC, HPROGRAM, DGUP, DGUPC, KMASK)          
ELSE IF (U%CTOWN=='NONE  ') THEN
  CALL INIT_BUD(DTO,DGUP,DGUPC,XUNDEF)         
END IF
IF (LHOOK) CALL DR_HOOK('DIAG_TOWN_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_TOWN_n
