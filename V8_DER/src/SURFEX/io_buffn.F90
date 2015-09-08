!     #######################################################
      SUBROUTINE IO_BUFF_n (IOB, &
                            HREC,HACTION,OKNOWN)
!     #######################################################
!
!!****  *IO_BUFF_n* - function to check if the field has already been read/written
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    08/2007 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
!
!
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
!
 CHARACTER(LEN=12),  INTENT(IN) :: HREC     ! field to read or write
 CHARACTER(LEN=1),   INTENT(IN) :: HACTION  ! 'R' : file being read
                                           ! 'W' : file being written
!
LOGICAL,            INTENT(OUT):: OKNOWN   ! T : field has already been read/written
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JLOOP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('IO_BUFF_N',0,ZHOOK_HANDLE)
IF (HACTION=='R' .OR. HACTION=='W') THEN
  OKNOWN=.FALSE.
  DO JLOOP=1,IOB%NREC
    OKNOWN=(HREC==IOB%CREC(JLOOP))
    IF (OKNOWN .AND. LHOOK) CALL DR_HOOK('IO_BUFF_N',1,ZHOOK_HANDLE)
    IF (OKNOWN) RETURN
  END DO
  IOB%NREC=IOB%NREC+1
  IOB%CREC(IOB%NREC)=HREC
END IF
IF (LHOOK) CALL DR_HOOK('IO_BUFF_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE IO_BUFF_n
