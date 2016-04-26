!##################
MODULE MODD_TEB_VEG_PARAM_n
!##################
!
!!****  *MODD_ISBA - declaration of packed surface parameters for ISBA scheme
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
!!      A. Boone   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       20/09/02
!!      A.L. Gibelin    04/2009 : BIOMASS and RESP_BIOMASS arrays 
!!      A.L. Gibelin    04/2009 : TAU_WOOD for NCB option 
!!      A.L. Gibelin    05/2009 : Add carbon spinup
!!      A.L. Gibelin    06/2009 : Soil carbon variables for CNT option
!!      A.L. Gibelin    07/2009 : Suppress RDK and transform GPP as a diagnostic
!!      A.L. Gibelin    07/2009 : Suppress PPST and PPSTF as outputs
!!      P. Samuelsson   02/2012 : MEB
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_ISBA_PARAM1_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
TYPE TEB_VEG_PARAM_TIME_t
  !
  TYPE(ISBA_PARAM_TIME_1P_t), POINTER :: ALP(:) => NULL()
  TYPE(ISBA_PARAM_TIME_1P_t), POINTER :: CUR => NULL()
  !
END TYPE TEB_VEG_PARAM_TIME_t        
!
TYPE TEB_VEG_PARAM_t
!
TYPE(TEB_VEG_PARAM_TIME_t) :: T
!
TYPE(ISBA_PARAM_FIX_1P_t) :: X
TYPE(ISBA_PARAM_MEB_1P_t) :: M
TYPE(ISBA_PARAM_ALB_1P_t) :: A
TYPE(ISBA_PARAM_IRRIG_1P_t) :: I
!
END TYPE TEB_VEG_PARAM_t
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
!
SUBROUTINE TEB_VEG_PARAM_TIME_GOTO_PATCH(YTEB_VEG_PARAM_TIME,KTO_PATCH)
TYPE(TEB_VEG_PARAM_TIME_t), INTENT(INOUT) :: YTEB_VEG_PARAM_TIME
INTEGER, INTENT(IN) :: KTO_PATCH
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current patch is set to patch KTO_PATCH
IF (LHOOK) CALL DR_HOOK('MODD_TEB_VEG_PARAM_N:TEB_VEG_PARAM_TIME_GOTO_PATCH',0,ZHOOK_HANDLE)
YTEB_VEG_PARAM_TIME%CUR => YTEB_VEG_PARAM_TIME%ALP(KTO_PATCH)

IF (LHOOK) CALL DR_HOOK('MODD_TEB_VEG_PARAM_N:TEB_VEG_PARAM_TIME_GOTO_PATCH',1,ZHOOK_HANDLE)
!
END SUBROUTINE TEB_VEG_PARAM_TIME_GOTO_PATCH
!
SUBROUTINE TEB_VEG_PARAM_TIME_INIT(YTEB_VEG_PARAM_TIME,KPATCH)
TYPE(TEB_VEG_PARAM_TIME_t), INTENT(INOUT) :: YTEB_VEG_PARAM_TIME
INTEGER, INTENT(IN) :: KPATCH
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_TEB_VEG_PARAM_N:TEB_VEG_PARAM_TIME_INIT",0,ZHOOK_HANDLE)
 ALLOCATE(YTEB_VEG_PARAM_TIME%ALP(KPATCH))
 YTEB_VEG_PARAM_TIME%CUR => YTEB_VEG_PARAM_TIME%ALP(1)
DO JP=1,KPATCH
  CALL ISBA_PARAM_TIME_1P_INIT(YTEB_VEG_PARAM_TIME%ALP(JP))
ENDDO 
IF (LHOOK) CALL DR_HOOK("MODD_TEB_VEG_PARAM_N:TEB_VEG_PARAM_TIME_INIT",1,ZHOOK_HANDLE)
END SUBROUTINE TEB_VEG_PARAM_TIME_INIT
!
END MODULE MODD_TEB_VEG_PARAM_n
