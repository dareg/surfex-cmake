!     #######################
      SUBROUTINE READ_LECOCLIMAP (HPROGRAM,OECOCLIMAP,HDIR)
!     #######################
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
 CHARACTER(LEN=6),     INTENT(IN)    :: HPROGRAM  ! program calling surf. schemes
LOGICAL,              INTENT(OUT)   :: OECOCLIMAP! flag for ecoclimap
 CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: HDIR
!
!
!* local variables
!  ---------------
!
 CHARACTER(LEN=1) :: YDIR
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
!
YDIR = 'H'
IF (PRESENT(HDIR)) YDIR = HDIR
!
YRECFM='VERSION'
 CALL READ_SURF(&
                HPROGRAM,YRECFM,IVERSION,IRESP,YDIR)
YRECFM='BUG'
 CALL READ_SURF(&
                HPROGRAM,YRECFM,IBUGFIX,IRESP,YDIR)
!
IF (IVERSION<1 .OR. (IVERSION==1 .AND. IBUGFIX==0)) THEN
  OECOCLIMAP = .TRUE.
ELSE
  YRECFM='ECOCLIMAP'
  CALL READ_SURF(&
                HPROGRAM,YRECFM,OECOCLIMAP,IRESP,YDIR)
END IF
IF (LHOOK) CALL DR_HOOK('READ_LECOCLIMAP',1,ZHOOK_HANDLE)
!
!------------------------------------------------------------------------------
!
END SUBROUTINE READ_LECOCLIMAP
