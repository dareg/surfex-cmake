SUBROUTINE TEST_RECORD_LEN(HPROGRAM,HREC,ONOWRITE)
!
USE MODI_GET_LUOUT
USE MODD_DIAG_SURF_ATM_n,  ONLY : LSELECT, CSELECT
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM ! calling program
CHARACTER(LEN=16),  INTENT(IN)  :: HREC     ! name of the article to be written
LOGICAL,            INTENT(OUT) :: ONOWRITE ! flag for article to be written
!
INTEGER :: IFIELD,JFIELD
INTEGER :: ILUOUT  ! listing logical unit
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODI_WRITE_SURF:TEST_RECORD_LEN',0,ZHOOK_HANDLE)
IF (LEN_TRIM(HREC)>16) THEN
  CALL GET_LUOUT(HPROGRAM,ILUOUT)
  WRITE(ILUOUT,*) '----------------------------------------------'
  WRITE(ILUOUT,*) 'Error occured when writing a field            '
  WRITE(ILUOUT,*) 'The name of the field is too long             '
  WRITE(ILUOUT,*) 'The name must not be longer than 16 characters'
  WRITE(ILUOUT,*) 'Please shorten the name of your field         '
  WRITE(ILUOUT,FMT='(A32,A16,A1)') ' The field name currently is : "',HREC,'"'
  WRITE(ILUOUT,*) '----------------------------------------------'
  CALL ABOR1_SFX('TEST_RECORD_LEN: FIELD NAME TOO LONG --> '//HREC)
END IF
!
! if output fields selection is active, test if this field is to be written
IF (LSELECT)  THEN
   IFIELD=COUNT(CSELECT /= '            ')
   ONOWRITE=.TRUE.
   DO JFIELD=1,IFIELD
      IF ( TRIM(CSELECT(JFIELD))==TRIM(HREC) ) THEN
         ONOWRITE=.FALSE.
      ENDIF
   ENDDO
ELSE
   ONOWRITE=.FALSE.
ENDIF
IF (LHOOK) CALL DR_HOOK('MODI_WRITE_SURF:TEST_RECORD_LEN',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE TEST_RECORD_LEN
