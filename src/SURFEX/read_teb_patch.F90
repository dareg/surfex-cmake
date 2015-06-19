!     #######################
      SUBROUTINE READ_TEB_PATCH(HFILEPGD,HFILEPGDTYPE,KTEB_PATCH)
!     #######################
!
USE MODD_IO_BUFF_n, ONLY : IOB => IO_BUFF
!
USE MODI_READ_SURF
USE MODI_TOWN_PRESENCE
!
USE MODI_OPEN_AUX_IO_SURF
USE MODI_CLOSE_AUX_IO_SURF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!* dummy arguments
!  ---------------
!
CHARACTER(LEN=28),  INTENT(IN)  :: HFILEPGD     ! name of file
CHARACTER(LEN=6),   INTENT(IN)  :: HFILEPGDTYPE ! type of file
INTEGER,            INTENT(OUT) :: KTEB_PATCH! number of TEB patches
!
!
!* local variables
!  ---------------
!
 CHARACTER(LEN=12) :: YRECFM     ! Name of the article to be read
INTEGER           :: IRESP      ! reading return code
!
INTEGER           :: IVERSION   ! surface version
INTEGER           :: IBUGFIX    ! surface bugfix
LOGICAL           :: GTOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('READ_TEB_PATCH',0,ZHOOK_HANDLE)
!
CALL OPEN_AUX_IO_SURF(IOB, &
                      HFILEPGD,HFILEPGDTYPE,'FULL  ')
YRECFM='VERSION'
CALL READ_SURF(IOB, &
               HFILEPGDTYPE,YRECFM,IVERSION,IRESP)
YRECFM='BUG'
CALL READ_SURF(IOB, &
               HFILEPGDTYPE,YRECFM,IBUGFIX,IRESP)
!
CALL TOWN_PRESENCE(HFILEPGDTYPE,GTOWN)
CALL CLOSE_AUX_IO_SURF(HFILEPGD,HFILEPGDTYPE)
!
IF (IVERSION<7 .OR. (IVERSION==7 .AND. IBUGFIX<=2).OR..NOT.GTOWN) THEN
  KTEB_PATCH = 1
ELSE
  YRECFM='TEB_PATCH'
  CALL OPEN_AUX_IO_SURF(IOB, &
                      HFILEPGD,HFILEPGDTYPE,'TOWN  ')
  CALL READ_SURF(IOB, &
               HFILEPGDTYPE,YRECFM,KTEB_PATCH,IRESP)
  CALL CLOSE_AUX_IO_SURF(HFILEPGD,HFILEPGDTYPE)
END IF
!
IF (LHOOK) CALL DR_HOOK('READ_TEB_PATCH',1,ZHOOK_HANDLE)
!
!------------------------------------------------------------------------------
!
END SUBROUTINE READ_TEB_PATCH
