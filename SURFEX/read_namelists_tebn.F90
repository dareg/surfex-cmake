!     #########
SUBROUTINE READ_NAMELISTS_TEB_n(HPROGRAM, HINIT)
!     #######################################################
!
!---------------------------------------------------------------------------
!
USE MODN_TEB_n                          
!
USE MODI_DEFAULT_TEB
USE MODI_DEFAULT_CH_DEP
USE MODI_DEFAULT_DIAG_TEB
USE MODI_READ_DEFAULT_TEB_n
USe MODI_READ_TEB_CONF_n
!
USE MODI_READ_NAM_PREP_TEB_n
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_READ_TEB_CONF_n
IMPLICIT NONE
!
CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
CHARACTER(LEN=3),   INTENT(IN)  :: HINIT     ! choice of fields to initialize
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('READ_NAMELISTS_TEB_N',0,ZHOOK_HANDLE)
CALL DEFAULT_TEB(CZ0H,XTSTEP,XOUT_TSTEP)
!
CALL DEFAULT_CH_DEP(CCH_DRY_DEP)
!
CALL DEFAULT_DIAG_TEB(N2M,LSURF_BUDGET,LSURF_MISC_BUDGET,LSURF_BUDGETC,LPGD,     &
                         LPGD_FIX,LRAD_BUDGET,XDIAG_TSTEP,LRESET_BUDGETC)   
!               
CALL READ_DEFAULT_TEB_n(HPROGRAM)
!
CALL READ_TEB_CONF_n(HPROGRAM) 
!      
!
IF (HINIT=='PRE') CALL READ_NAM_PREP_TEB_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('READ_NAMELISTS_TEB_N',1,ZHOOK_HANDLE)
!
!------------------------------------
!
END SUBROUTINE READ_NAMELISTS_TEB_n
