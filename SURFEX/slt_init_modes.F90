!     #########
SUBROUTINE SLT_INIT_MODES (KSLTEQ,KSV_SLTBEG,KSV_SLTEND,OVARSIG,ORGFIX, KSLT_MDEBEG,KSLTMDE)
!!    ###########################################
!!
!!*** *SLT_INIT_MODES*
!!
!!    PURPOSE
!!    -------
!!    Find the number of sea salt modes to be transported
!!    Each mode needs 3 moments to be described, so logically, the number of modes is
!!    The number of sea salt tracers divided by 3
!!     
!!
!!    REFERENCE
!!    ---------
!!    Modified slt_init_names (march 2005)    
!!
!!    AUTHOR
!!    ------
!!    Alf Grini / P. Tulet
!!
!!    MODIFICATIONS
!!    -------------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!

INTEGER,                         INTENT(IN) :: KSLTEQ       ! number of sea salt variables
INTEGER,                         INTENT(IN) :: KSV_SLTBEG   ! First number of sea salt tracer
INTEGER,                         INTENT(IN) :: KSV_SLTEND   ! Last number of sea salt tracer
LOGICAL,                         INTENT(IN) :: OVARSIG      ! type of standard deviation (fixed or variable)
LOGICAL,                         INTENT(IN) :: ORGFIX       !type of mean radius
INTEGER,                         INTENT(OUT) :: KSLT_MDEBEG  ! Place in scalar list of sea saltmass in first mode
INTEGER,                         INTENT(OUT) :: KSLTMDE     ! Number of sea salt modes
REAL(KIND=JPRB) :: ZHOOK_HANDLE


!Check if you have a multiple of 3 sea salt related variables, and 
!Set the number of modes to the number of sea salt related variables
!divided by 3

IF (LHOOK) CALL DR_HOOK('SLT_INIT_MODES',0,ZHOOK_HANDLE)
IF (OVARSIG) THEN
  IF(mod((KSV_SLTEND - KSV_SLTBEG + 1),3).ne.0.)THEN
     CALL ABOR1_SFX('SLT_INIT_MODES: (1) WRONG SEA SALT RELATE VARIABLES')
  ELSE
     KSLT_MDEBEG=KSV_SLTBEG
     KSLTMDE=(KSV_SLTEND - KSV_SLTBEG + 1)/3
  ENDIF
ELSE IF (ORGFIX) THEN
   KSLT_MDEBEG=KSV_SLTBEG
   KSLTMDE=KSV_SLTEND - KSV_SLTBEG + 1
ELSE
  IF(mod((KSV_SLTEND - KSV_SLTBEG + 1),2).ne.0.)THEN
     CALL ABOR1_SFX('SLT_INIT_MODES: (2) WRONG SEA SALT RELATE VARIABLES')
  ELSE
     KSLT_MDEBEG=KSV_SLTBEG
     KSLTMDE=(KSV_SLTEND - KSV_SLTBEG + 1)/2
  END IF
END IF
IF (LHOOK) CALL DR_HOOK('SLT_INIT_MODES',1,ZHOOK_HANDLE)

END SUBROUTINE SLT_INIT_MODES
