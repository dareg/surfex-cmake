MODULE MODE_WRITE_SURF_ASC
!
INTERFACE WRITE_SURF0_ASC
        MODULE PROCEDURE WRITE_SURFX0_ASC
        MODULE PROCEDURE WRITE_SURFN0_ASC
        MODULE PROCEDURE WRITE_SURFL0_ASC
        MODULE PROCEDURE WRITE_SURFC0_ASC
END INTERFACE
INTERFACE WRITE_SURFN_ASC
        MODULE PROCEDURE WRITE_SURFX1_ASC
        MODULE PROCEDURE WRITE_SURFN1_ASC
        MODULE PROCEDURE WRITE_SURFL1_ASC
        MODULE PROCEDURE WRITE_SURFX2_ASC
END INTERFACE
INTERFACE WRITE_SURFT_ASC
        MODULE PROCEDURE WRITE_SURFT0_ASC
        MODULE PROCEDURE WRITE_SURFT1_ASC
        MODULE PROCEDURE WRITE_SURFT2_ASC
END INTERFACE
!
CONTAINS
!
!     #############################################################
      SUBROUTINE WRITE_SURFX0_ASC(HREC,PFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  * - routine to write a real scalar
!
USE MODI_ERROR_WRITE_SURF_ASC
!
USE MODD_IO_SURF_ASC,        ONLY : NUNIT, CMASK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_IO_BUFF_n
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=16),  INTENT(IN) :: HREC     ! name of the article to be read
REAL,               INTENT(IN) :: PFIELD   ! the real scalar to be read
INTEGER,            INTENT(OUT):: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(IN) :: HCOMMENT ! comment string
!
!*      0.2   Declarations of local variables
!
LOGICAL         :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFX0_ASC',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',GKNOWN)
IF (GKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFX0_ASC',1,ZHOOK_HANDLE)
IF (GKNOWN) RETURN
!
WRITE(NUNIT,FMT=*,ERR=100) '&'//CMASK//' '//HREC
WRITE(NUNIT,FMT='(A50)',ERR=100) HCOMMENT(1:50)
WRITE(NUNIT,FMT=*,ERR=100) PFIELD
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFX0_ASC',1,ZHOOK_HANDLE)
RETURN
!
100 CONTINUE
CALL ERROR_WRITE_SURF_ASC(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFX0_ASC',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFX0_ASC
!
!     #############################################################
      SUBROUTINE WRITE_SURFN0_ASC(HREC,KFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  * - routine to write an integer
!
USE MODI_ERROR_WRITE_SURF_ASC
!
USE MODD_IO_SURF_ASC,        ONLY : NUNIT, NMASK, CMASK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_IO_BUFF_n
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=16),  INTENT(IN) :: HREC     ! name of the article to be read
INTEGER,            INTENT(IN) :: KFIELD   ! the integer to be read
INTEGER,            INTENT(OUT):: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(IN) :: HCOMMENT ! comment string
!
!*      0.2   Declarations of local variables
!
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFN0_ASC',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',GKNOWN)
IF (GKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFN0_ASC',1,ZHOOK_HANDLE)
IF (GKNOWN) RETURN
!
WRITE(NUNIT,FMT=*,ERR=100) '&'//CMASK//' '//HREC
WRITE(NUNIT,FMT='(A50)',ERR=100) HCOMMENT(1:50)
WRITE(NUNIT,FMT=*,ERR=100) KFIELD
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFN0_ASC',1,ZHOOK_HANDLE)
RETURN
100 CONTINUE
CALL ERROR_WRITE_SURF_ASC(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFN0_ASC',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFN0_ASC
!
!     #############################################################
      SUBROUTINE WRITE_SURFL0_ASC(HREC,OFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  * - routine to write a logical
!
USE MODI_ERROR_WRITE_SURF_ASC
!
USE MODD_IO_SURF_ASC,        ONLY : NUNIT, CMASK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_IO_BUFF_n
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=16),  INTENT(IN) :: HREC     ! name of the article to be read
LOGICAL,            INTENT(IN) :: OFIELD   ! array containing the data field
INTEGER,            INTENT(OUT):: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(IN) :: HCOMMENT ! comment string
!
!*      0.2   Declarations of local variables
!
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFL0_ASC',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',GKNOWN)
IF (GKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFL0_ASC',1,ZHOOK_HANDLE)
IF (GKNOWN) RETURN
!
WRITE(NUNIT,FMT=*,ERR=100) '&'//CMASK//' '//HREC
WRITE(NUNIT,FMT='(A50)',ERR=100) HCOMMENT(1:50)
WRITE(NUNIT,FMT=*,ERR=100) OFIELD
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFL0_ASC',1,ZHOOK_HANDLE)
RETURN
100 CONTINUE
CALL ERROR_WRITE_SURF_ASC(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFL0_ASC',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFL0_ASC
!
!     #############################################################
      SUBROUTINE WRITE_SURFC0_ASC(HREC,HFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  * - routine to write a character
!
USE MODI_ERROR_WRITE_SURF_ASC
!
USE MODD_IO_SURF_ASC,        ONLY : NUNIT, CMASK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_IO_BUFF_n
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=16),  INTENT(IN)  :: HREC      ! name of the article to be read
CHARACTER(LEN=40),  INTENT(IN)  :: HFIELD    ! the integer to be read
INTEGER,            INTENT(OUT) :: KRESP     ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(IN)  :: HCOMMENT  ! comment string
!
!*      0.2   Declarations of local variables
!
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFC0_ASC',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',GKNOWN)
IF (GKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFC0_ASC',1,ZHOOK_HANDLE)
IF (GKNOWN) RETURN
!
WRITE(NUNIT,FMT=*,ERR=100) '&'//CMASK//' '//HREC
WRITE(NUNIT,FMT='(A50)',ERR=100) HCOMMENT(1:50)
WRITE(NUNIT,FMT='(A40)',ERR=100) HFIELD

IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFC0_ASC',1,ZHOOK_HANDLE)
RETURN
100 CONTINUE
CALL ERROR_WRITE_SURF_ASC(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFC0_ASC',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFC0_ASC
!
!     #############################################################
      SUBROUTINE WRITE_SURFX1_ASC(HREC,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  * - routine to fill a write 1D array for the externalised surface 
!
USE MODI_ERROR_WRITE_SURF_ASC
USE MODI_UNPACK_SAME_RANK
!
USE MODD_IO_SURF_ASC,        ONLY : NUNIT, NMASK, NFULL, CMASK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_IO_BUFF_n
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=16),   INTENT(IN) :: HREC     ! name of the article to be read
REAL, DIMENSION(:), INTENT(IN) :: PFIELD   ! array containing the data field
INTEGER,             INTENT(OUT):: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),  INTENT(IN) :: HCOMMENT ! comment string
CHARACTER(LEN=1),    INTENT(IN) :: HDIR     ! type of field :
                                            ! 'H' : field with
                                            !       horizontal spatial dim.
                                            ! '-' : no horizontal dim.
!
!*      0.2   Declarations of local variables
!
REAL, DIMENSION(NFULL) :: ZWORK   ! work array read in the file
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFX1_ASC',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',GKNOWN)
IF (GKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFX1_ASC',1,ZHOOK_HANDLE)
IF (GKNOWN) RETURN
!
WRITE(NUNIT,FMT=*,ERR=100) '&'//CMASK//' '//HREC
WRITE(NUNIT,FMT='(A50)',ERR=100) HCOMMENT(1:50)
!
IF (HDIR=='-') THEN
  WRITE(NUNIT,FMT='(50D20.8)',ERR=100) PFIELD
ELSE
  CALL UNPACK_SAME_RANK(NMASK,PFIELD,ZWORK(:))
  WRITE(NUNIT,FMT='(50D20.8)',ERR=100) ZWORK
END IF
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFX1_ASC',1,ZHOOK_HANDLE)
RETURN
!
100 CONTINUE
CALL ERROR_WRITE_SURF_ASC(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFX1_ASC',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFX1_ASC
!
!     #############################################################
      SUBROUTINE WRITE_SURFN1_ASC(HREC,KFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  * - routine to write an integer array
!
USE MODI_ERROR_WRITE_SURF_ASC
USE MODI_UNPACK_SAME_RANK
!
USE MODD_IO_SURF_ASC,        ONLY : NUNIT, NMASK, NFULL, CMASK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_IO_BUFF_n
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=16),      INTENT(IN) :: HREC     ! name of the article to be read
INTEGER, DIMENSION(:),  INTENT(IN) :: KFIELD   ! the integer to be read
INTEGER,                INTENT(OUT):: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),     INTENT(IN) :: HCOMMENT ! comment string
CHARACTER(LEN=1),       INTENT(IN) :: HDIR     ! type of field :
                                               ! 'H' : field with
                                               !       horizontal spatial dim.
                                               ! '-' : no horizontal dim.
!
!*      0.2   Declarations of local variables
!
INTEGER, DIMENSION(NFULL) :: IWORK  ! work array read in the file
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!---------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFN1_ASC',0,ZHOOK_HANDLE)
KRESP = 0
!
CALL IO_BUFF_n(HREC,'W',GKNOWN)
IF (GKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFN1_ASC',1,ZHOOK_HANDLE)
IF (GKNOWN) RETURN
!
WRITE(NUNIT,FMT=*,ERR=100) '&'//CMASK//' '//HREC
WRITE(NUNIT,FMT='(A50)',ERR=100) HCOMMENT(1:50)
!
IF (HDIR=='H' .AND. HREC/="-") THEN
  CALL UNPACK_SAME_RANK(NMASK,KFIELD,IWORK(:))
  WRITE(NUNIT,FMT='(100I8)',ERR=100) IWORK
ELSE
  WRITE(NUNIT,FMT='(100I8)',ERR=100) KFIELD
END IF
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFN1_ASC',1,ZHOOK_HANDLE)
RETURN
!
100 CONTINUE
CALL ERROR_WRITE_SURF_ASC(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFN1_ASC',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFN1_ASC
!
!     #############################################################
      SUBROUTINE WRITE_SURFL1_ASC(HREC,OFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  * - routine to write a logical array
!
USE MODI_ERROR_WRITE_SURF_ASC
!
USE MODD_IO_SURF_ASC,        ONLY : NUNIT, CMASK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_IO_BUFF_n
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=16),      INTENT(IN) :: HREC     ! name of the article to be read
LOGICAL, DIMENSION(:),  INTENT(IN) :: OFIELD   ! array containing the data field
INTEGER,                INTENT(OUT):: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),     INTENT(IN) :: HCOMMENT ! comment string
CHARACTER(LEN=1),       INTENT(IN) :: HDIR     ! type of field :
                                               ! 'H' : field with
                                               !       horizontal spatial dim.
                                               ! '-' : no horizontal dim.
!IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFL1_ASC',1,ZHOOK_HANDLE)

!*      0.2   Declarations of local variables
!
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFL1_ASC',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',GKNOWN)
IF (GKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFL1_ASC',1,ZHOOK_HANDLE)
IF (GKNOWN) RETURN
!
WRITE(NUNIT,FMT=*,ERR=100) '&'//CMASK//' '//HREC
WRITE(NUNIT,FMT='(A50)',ERR=100) HCOMMENT(1:50)
WRITE(NUNIT,FMT=*,ERR=100) OFIELD
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFL1_ASC',1,ZHOOK_HANDLE)
RETURN
!
100 CONTINUE
CALL ERROR_WRITE_SURF_ASC(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFL1_ASC',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFL1_ASC
!
!     #############################################################
      SUBROUTINE WRITE_SURFX2_ASC(HREC,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  * - routine to fill a write 2D array for the externalised surface 
!
USE MODI_ERROR_WRITE_SURF_ASC
USE MODI_UNPACK_SAME_RANK
!
USE MODD_IO_SURF_ASC,        ONLY : NUNIT, NMASK, NFULL, CMASK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_IO_BUFF_n
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=16),        INTENT(IN) :: HREC     ! name of the article to be read
REAL, DIMENSION(:,:),     INTENT(IN) :: PFIELD   ! array containing the data field
INTEGER,                  INTENT(OUT):: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),       INTENT(IN) :: HCOMMENT ! comment string
CHARACTER(LEN=1),         INTENT(IN) :: HDIR     ! type of field :
                                                 ! 'H' : field with
                                                 !       horizontal spatial dim.
                                                 ! '-' : no horizontal dim.
!
!*      0.2   Declarations of local variables
! 
REAL, DIMENSION(NFULL,SIZE(PFIELD,2)) :: ZWORK   ! work array read in the file
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFX2_ASC',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',GKNOWN)
IF (GKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFX2_ASC',1,ZHOOK_HANDLE)
IF (GKNOWN) RETURN
!
WRITE(NUNIT,FMT=*,ERR=100) '&'//CMASK//' '//HREC
WRITE(NUNIT,FMT='(A50)',ERR=100) HCOMMENT(1:50)
!
IF (HDIR=='-') THEN
  WRITE(NUNIT,FMT='(50D20.8)',ERR=100) PFIELD
ELSE
  CALL UNPACK_SAME_RANK(NMASK,PFIELD,ZWORK(:,:))
  WRITE(NUNIT,FMT='(50D20.8)',ERR=100) ZWORK
END IF
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFX2_ASC',1,ZHOOK_HANDLE)
RETURN
!
100 CONTINUE
CALL ERROR_WRITE_SURF_ASC(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFX2_ASC',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFX2_ASC
!
!     #############################################################
      SUBROUTINE WRITE_SURFT0_ASC(HREC,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
!     #############################################################
!
!!****  * - routine to write a date
!
USE MODI_ERROR_WRITE_SURF_ASC
!
USE MODD_IO_SURF_ASC,        ONLY : NUNIT, CMASK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_IO_BUFF_n
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=16),  INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER,            INTENT(IN)  :: KYEAR    ! year
INTEGER,            INTENT(IN)  :: KMONTH   ! month
INTEGER,            INTENT(IN)  :: KDAY     ! day
REAL,               INTENT(IN)  :: PTIME    ! time
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(IN)  :: HCOMMENT ! comment string

!*      0.2   Declarations of local variables
!
INTEGER, DIMENSION(3) :: ITDATE
LOGICAL               :: GKNOWN
REAL(KIND=JPRB)       :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFT0_ASC',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',GKNOWN)
IF (GKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFT0_ASC',1,ZHOOK_HANDLE)
IF (GKNOWN) RETURN
!
ITDATE(1) = KYEAR
ITDATE(2) = KMONTH
ITDATE(3) = KDAY
!
WRITE(NUNIT,FMT=*,ERR=100) '&'//CMASK//' '//TRIM(HREC)//'%TDATE'
WRITE(NUNIT,FMT='(A50)',ERR=100) HCOMMENT(1:50)
WRITE(NUNIT,FMT=*,ERR=100) ITDATE(:)
!
WRITE(NUNIT,FMT=*,ERR=100) '&'//CMASK//' '//TRIM(HREC)//'%TIME'
WRITE(NUNIT,FMT='(A50)',ERR=100) HCOMMENT(1:50)
WRITE(NUNIT,FMT=*,ERR=100) PTIME
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFT0_ASC',1,ZHOOK_HANDLE)
RETURN
!
100 CONTINUE
CALL ERROR_WRITE_SURF_ASC(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFT0_ASC',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFT0_ASC
!
!     #############################################################
      SUBROUTINE WRITE_SURFT1_ASC(HREC,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
!     #############################################################
!
!!****  * - routine to write a date
!
USE MODI_ERROR_WRITE_SURF_ASC
!
USE MODD_IO_SURF_ASC,        ONLY : NUNIT, CMASK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_IO_BUFF_n
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=16),     INTENT(IN) :: HREC     ! name of the article to be read
INTEGER, DIMENSION(:), INTENT(IN) :: KYEAR    ! year
INTEGER, DIMENSION(:), INTENT(IN) :: KMONTH   ! month
INTEGER, DIMENSION(:), INTENT(IN) :: KDAY     ! day
REAL,    DIMENSION(:), INTENT(IN) :: PTIME    ! time
INTEGER,               INTENT(OUT):: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),    INTENT(IN) :: HCOMMENT ! comment string

!*      0.2   Declarations of local variables
!
INTEGER, DIMENSION(3,SIZE(KYEAR)) :: ITDATE
LOGICAL                           :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFT1_ASC',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',GKNOWN)
IF (GKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFT1_ASC',1,ZHOOK_HANDLE)
IF (GKNOWN) RETURN
!
ITDATE(1,:) = KYEAR  (:)
ITDATE(2,:) = KMONTH (:)
ITDATE(3,:) = KDAY   (:)

WRITE(NUNIT,FMT=*,ERR=100) '&'//CMASK//' '//TRIM(HREC)//'%TDATE'
WRITE(NUNIT,FMT='(A50)',ERR=100) HCOMMENT(1:50)
WRITE(NUNIT,FMT=*,ERR=100) ITDATE(:,:)
!
WRITE(NUNIT,FMT=*,ERR=100) '&'//CMASK//' '//TRIM(HREC)//'%TIME'
WRITE(NUNIT,FMT='(A50)',ERR=100) HCOMMENT(1:50)
WRITE(NUNIT,FMT=*,ERR=100) PTIME
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFT1_ASC',1,ZHOOK_HANDLE)
RETURN
!
100 CONTINUE
CALL ERROR_WRITE_SURF_ASC(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFT1_ASC',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFT1_ASC
!
!     #############################################################
      SUBROUTINE WRITE_SURFT2_ASC(HREC,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
!     #############################################################
!
!!****  * - routine to write a date
!
USE MODI_ERROR_WRITE_SURF_ASC
!
USE MODD_IO_SURF_ASC,        ONLY : NUNIT, CMASK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_IO_BUFF_n
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=16),       INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER, DIMENSION(:,:), INTENT(IN)  :: KYEAR    ! year
INTEGER, DIMENSION(:,:), INTENT(IN)  :: KMONTH   ! month
INTEGER, DIMENSION(:,:), INTENT(IN)  :: KDAY     ! day
REAL,    DIMENSION(:,:), INTENT(IN)  :: PTIME    ! time
INTEGER,                 INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),      INTENT(IN)  :: HCOMMENT ! comment string

!*      0.2   Declarations of local variables
!
INTEGER, DIMENSION(3,SIZE(KYEAR,1),SIZE(KYEAR,2)) :: ITDATE
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFT2_ASC',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',GKNOWN)
IF (GKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFT2_ASC',1,ZHOOK_HANDLE)
IF (GKNOWN) RETURN
!
ITDATE(1,:,:) = KYEAR  (:,:)
ITDATE(2,:,:) = KMONTH (:,:)
ITDATE(3,:,:) = KDAY   (:,:)

WRITE(NUNIT,FMT=*,ERR=100) '&'//CMASK//' '//TRIM(HREC)//'%TDATE'
WRITE(NUNIT,FMT='(A50)',ERR=100) HCOMMENT(1:50)
WRITE(NUNIT,FMT=*,ERR=100) ITDATE(:,:,:)
!
WRITE(NUNIT,FMT=*,ERR=100) '&'//CMASK//' '//TRIM(HREC)//'%TIME'
WRITE(NUNIT,FMT='(A50)',ERR=100) HCOMMENT(1:50)
WRITE(NUNIT,FMT=*,ERR=100) PTIME
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFT2_ASC',1,ZHOOK_HANDLE)
RETURN
!
100 CONTINUE
CALL ERROR_WRITE_SURF_ASC(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_ASC:WRITE_SURFT2_ASC',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_SURFT2_ASC
!
END MODULE MODE_WRITE_SURF_ASC
