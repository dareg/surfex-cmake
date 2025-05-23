!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
SUBROUTINE ASSIM_GATHER_WRITE_INCREMENTS(HPROGRAM,HVAR,PFIELD_OLD,PFIELD_NEW,PUNDEF,LSTAT)
  USE PARKIND1,        ONLY : JPIM
  USE YOMHOOK ,        ONLY : LHOOK, DR_HOOK, JPHOOK
  USE MODD_SURF_PAR,   ONLY : XUNDEF
  USE MODD_SURFEX_MPI, ONLY : NINDEX
  USE MODD_ASSIM,      ONLY : LPIO
  USE MODI_GET_LUOUT
  USE MODI_ASSIM_GATHER_FIELD
  USE MODI_ABOR1_SFX
  USE MODD_SURFEX_HOST
  IMPLICIT NONE
  CHARACTER(LEN=6),  INTENT(IN) :: HPROGRAM
  CHARACTER(LEN=*),  INTENT(IN) :: HVAR
  REAL,DIMENSION(:), INTENT(IN) :: PFIELD_OLD
  REAL,DIMENSION(:), INTENT(IN) :: PFIELD_NEW
  REAL,              INTENT(IN) :: PUNDEF
  LOGICAL,OPTIONAL,  INTENT(IN) :: LSTAT
  INTEGER                       :: IOBS,ISIZE,JI,ILUOUT
  REAL,DIMENSION(:),ALLOCATABLE :: ZWORK,ZWORK_NEW,ZWORK_OLD
  CHARACTER(LEN=70)             :: YOUTFORMAT1, YOUTFORMAT2, YOUTFORMAT3

  LOGICAL,ALLOCATABLE,DIMENSION(:) :: LMASK
  LOGICAL                          :: OSTAT

  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('ASSIM_GATHER_WRITE_INCREMENTS',0,ZHOOK_HANDLE)

  OSTAT=.FALSE.
  IF (PRESENT(LSTAT)) OSTAT=LSTAT

  CALL GET_LUOUT(HPROGRAM,ILUOUT)
  ! Get global size
  IF(HPROGRAM == 'AROME ') THEN
    ISIZE = YRSURFEX_HOST%ASSIM_GET_SIZE_GLOBAL(PFIELD_OLD)
  ELSE
    ISIZE=SIZE(NINDEX)
  ENDIF

  IF ( LPIO ) THEN
    ALLOCATE(ZWORK_OLD(ISIZE))
    ALLOCATE(ZWORK_NEW(ISIZE))
  ELSE
    ALLOCATE(ZWORK_OLD(0))
    ALLOCATE(ZWORK_NEW(0))
  ENDIF

  CALL ASSIM_GATHER_FIELD(HPROGRAM,PFIELD_OLD,ZWORK_OLD)
  CALL ASSIM_GATHER_FIELD(HPROGRAM,PFIELD_NEW,ZWORK_NEW)

  IF ( LPIO ) THEN

    ALLOCATE(ZWORK(ISIZE))
    ALLOCATE(LMASK(ISIZE))

    LMASK(:)=.FALSE. 
    LMASK(:)=(ZWORK_OLD(:) /= PUNDEF .AND. ZWORK_OLD(:) /= XUNDEF .AND. ZWORK_NEW(:) /= PUNDEF .AND. ZWORK_NEW(:) /= XUNDEF)
    IOBS=COUNT(LMASK)
    WHERE ( LMASK(:) )
      ZWORK(:) = ZWORK_NEW(:) - ZWORK_OLD(:)
    ELSEWHERE
      ZWORK(:) = 0.
    ENDWHERE
    IF (IOBS > 0 ) THEN
      YOUTFORMAT1 = "('Mean increment for ',A,': ',E14.4,'   Defined values:',I9)"
      WRITE(*,YOUTFORMAT1) TRIM(HVAR), SUM(ZWORK(:))/REAL(IOBS), IOBS
    ELSE
      WRITE(ILUOUT,*) 'No defined values found for variable '//TRIM(HVAR)
    ENDIF
    IF ( OSTAT ) THEN
      IF ( COUNT(LMASK) > 0 ) THEN
        WRITE(ILUOUT,'(A)') 'Variable: '//TRIM(HVAR)
        YOUTFORMAT2 = "('New field, min,mean,max: ',3F12.4)"
        WRITE(ILUOUT,YOUTFORMAT2) MINVAL(ZWORK_NEW,MASK=LMASK),&
          SUM(ZWORK_NEW,MASK=LMASK)/COUNT(LMASK),MAXVAL(ZWORK_NEW,MASK=LMASK)
        YOUTFORMAT3 = "('Old field, min,mean,max: ',3F12.4)"
        WRITE(ILUOUT,YOUTFORMAT3) MINVAL(ZWORK_OLD,MASK=LMASK),&
          SUM(ZWORK_OLD,MASK=LMASK)/COUNT(LMASK),MAXVAL(ZWORK_OLD,MASK=LMASK)
      ELSE
         WRITE(ILUOUT,*) 'All values masked out'
      ENDIF
    ENDIF
    DEALLOCATE(LMASK)
    DEALLOCATE(ZWORK)
  ENDIF
  DEALLOCATE(ZWORK_NEW)
  DEALLOCATE(ZWORK_OLD)

  IF (LHOOK) CALL DR_HOOK('ASSIM_GATHER_WRITE_INCREMENTS',1,ZHOOK_HANDLE)
END SUBROUTINE ASSIM_GATHER_WRITE_INCREMENTS

