!     #######################
      SUBROUTINE READ_LECOCLIMAP (IOB, &
                                  HPROGRAM,OECOCLIMAP)
!     #######################
!
!
!
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
!
USE MODI_READ_SURF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!* dummy arguments
!  ---------------
!
!
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
!
 CHARACTER(LEN=6),     INTENT(IN)    :: HPROGRAM  ! program calling surf. schemes
LOGICAL,              INTENT(OUT)   :: OECOCLIMAP! flag for ecoclimap
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
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('READ_LECOCLIMAP',0,ZHOOK_HANDLE)
YRECFM='VERSION'
 CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,IVERSION,IRESP)
YRECFM='BUG'
 CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,IBUGFIX,IRESP)
!
IF (IVERSION<1 .OR. (IVERSION==1 .AND. IBUGFIX==0)) THEN
  OECOCLIMAP = .TRUE.
ELSE
  YRECFM='ECOCLIMAP'
  CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,OECOCLIMAP,IRESP)
END IF
IF (LHOOK) CALL DR_HOOK('READ_LECOCLIMAP',1,ZHOOK_HANDLE)
!
!------------------------------------------------------------------------------
!
END SUBROUTINE READ_LECOCLIMAP
