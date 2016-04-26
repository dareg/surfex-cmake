!     #########
SUBROUTINE DIAG_NATURE_n (DE, DLO, DL, DLC, DIO, DI, DIC, U, &
                          HPROGRAM, DUP, DUPC, KMASK )
!     ###############################################################################
!
!!****  *DIAG_NATURE_n * - Chooses the surface schemes for diagnostics over
!!    natural continental parts
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
!       P. Le Moigne 03/2015 : add diagnostics IDEAL case
!!------------------------------------------------------------------
!
USE MODE_DIAG
!
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
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
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DE
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DLO
TYPE(DIAG_t), INTENT(INOUT) :: DL
TYPE(DIAG_t), INTENT(INOUT) :: DLC
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DIO
TYPE(DIAG_t), INTENT(INOUT) :: DI
TYPE(DIAG_t), INTENT(INOUT) :: DIC
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
INTEGER :: ISIZE, JJ
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DIAG_NATURE_N',0,ZHOOK_HANDLE)
IF (U%CNATURE=='ISBA  ' .OR. U%CNATURE=='TSZ0  ' ) THEN
  !  
  CALL DIAG_CUMUL(DIO, DI, DIC, HPROGRAM, DUP, DUPC, KMASK)
  !
  ISIZE = SIZE(KMASK)
  !
  IF (DE%LSURF_EVAP_BUDGET) THEN
    DO JJ=1,ISIZE        
      DUP%XEVAP    (KMASK(JJ))  = DI%XEVAP     (JJ)
      DUP%XSUBL    (KMASK(JJ))  = DI%XSUBL     (JJ)
    ENDDO
  ENDIF
  !                   
ELSE IF (U%CNATURE=='FLUX  ') THEN
  CALL DIAG_EVAP(DLO, DL, DLC, HPROGRAM, DUP, DUPC, KMASK)
ELSE IF (U%CNATURE=='NONE  ') THEN
  CALL INIT_BUD(DIO,DUP,DUPC,XUNDEF)
END IF
IF (LHOOK) CALL DR_HOOK('DIAG_NATURE_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_NATURE_n
