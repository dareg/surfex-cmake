!     #########
SUBROUTINE TOWN_PRESENCE (HFILETYPE,OTEB,HDIR)
!     #################################################################################
!
!
!
!
USE MODI_READ_SURF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HFILETYPE ! type of input file
LOGICAL,            INTENT(OUT) :: OTEB      ! TRUE if TEB data exist in the file
 CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: HDIR
!
!*      0.2    declarations of local variables
!
 CHARACTER(LEN=1) :: YDIR
INTEGER           :: IRESP     ! reading return code
 CHARACTER(LEN=6)  :: YTOWN     ! scheme for towns in input file
INTEGER           :: IDIM_TOWN ! number of TEB points in input file
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
!*      1.     reads if TEB fields exist in the input file
!              -------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TOWN_PRESENCE',0,ZHOOK_HANDLE)
!
YDIR = 'H'
IF (PRESENT(HDIR)) YDIR = HDIR
!
YTOWN = ''
IDIM_TOWN = 0
 CALL READ_SURF(HFILETYPE,'TOWN',YTOWN,IRESP,YDIR)
 CALL READ_SURF(HFILETYPE,'DIM_TOWN',IDIM_TOWN,IRESP,YDIR)
!
OTEB = (YTOWN=='TEB   ') .AND. (IDIM_TOWN > 0)
IF (LHOOK) CALL DR_HOOK('TOWN_PRESENCE',1,ZHOOK_HANDLE)
!
!---------------------------------------------------------------------------------------
!
END SUBROUTINE TOWN_PRESENCE
