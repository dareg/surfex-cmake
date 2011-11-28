!     #########
SUBROUTINE DST_INIT_MODES (KDSTEQ,KSV_DSTBEG,KSV_DSTEND,OVARSIG,ORGFIX,KDST_MDEBEG,KDSTMDE)
!!    ###########################################
!!
!!*** *DST_INIT_MODES*
!!
!!    PURPOSE
!!    -------
!!    Find the number of dust modes to be transported
!!    Each mode needs 3 moments to be described, so logically, the number of modes is
!!    The number of dust tracers divided by 3
!!     
!!
!!    REFERENCE
!!    ---------
!!    Modified dst_init_names (march 2005)    
!!
!!    AUTHOR
!!    ------
!!    Alf Grini <alf.grini@cnrm.meteo.fr>
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

INTEGER,                         INTENT(IN) :: KDSTEQ       ! number of dust variables
INTEGER,                         INTENT(IN) :: KSV_DSTBEG   ! First number of dust tracer
INTEGER,                         INTENT(IN) :: KSV_DSTEND   ! Last number of dust tracer
LOGICAL,                         INTENT(IN) :: OVARSIG      ! type of standard deviation (fixed or variable)
LOGICAL,                         INTENT(IN) :: ORGFIX       !type of mean radius
INTEGER,                         INTENT(OUT) :: KDST_MDEBEG  ! Place in scalar list of dustmass in first mode
INTEGER,                         INTENT(OUT) :: KDSTMDE     ! Number of dust modes
REAL(KIND=JPRB) :: ZHOOK_HANDLE


!Check if you have a multiple of 3 dust related variables, and 
!Set the number of modes to the number of dust related variables
!divided by 3
IF (LHOOK) CALL DR_HOOK('DST_INIT_MODES',0,ZHOOK_HANDLE)
IF (OVARSIG) THEN !case three moments by modes
  IF(mod((KSV_DSTEND - KSV_DSTBEG + 1),3).ne.0.)THEN
    CALL ABOR1_SFX('DST_INIT_MODES: (1) WRONG NUMBER OF DUST VARIABLES')
  ELSE
   KDST_MDEBEG=KSV_DSTBEG
   KDSTMDE=(KSV_DSTEND - KSV_DSTBEG + 1)/3
  ENDIF
ELSE IF (ORGFIX) THEN ! case one moment by modes
   KDST_MDEBEG=KSV_DSTBEG
   KDSTMDE=KSV_DSTEND - KSV_DSTBEG + 1
ELSE  ! case two moments by modes
  IF(mod((KSV_DSTEND - KSV_DSTBEG + 1),2).ne.0.)THEN
   CALL ABOR1_SFX('DST_INIT_MODES: (1) WRONG NUMBER OF DUST VARIABLES')
  ELSE
   KDST_MDEBEG=KSV_DSTBEG
   KDSTMDE=(KSV_DSTEND - KSV_DSTBEG + 1)/2
END IF
END IF
IF (LHOOK) CALL DR_HOOK('DST_INIT_MODES',1,ZHOOK_HANDLE)

END SUBROUTINE DST_INIT_MODES
