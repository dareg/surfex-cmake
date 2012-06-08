MODULE MODE_READ_SURF_ASC
!
INTERFACE READ_SURF0_ASC
        MODULE PROCEDURE READ_SURFX0_ASC
        MODULE PROCEDURE READ_SURFN0_ASC
        MODULE PROCEDURE READ_SURFL0_ASC
        MODULE PROCEDURE READ_SURFC0_ASC
END INTERFACE
INTERFACE READ_SURFN_ASC
        MODULE PROCEDURE READ_SURFX1_ASC
        MODULE PROCEDURE READ_SURFN1_ASC
        MODULE PROCEDURE READ_SURFL1_ASC
        MODULE PROCEDURE READ_SURFX2_ASC
END INTERFACE
INTERFACE READ_SURFT_ASC
        MODULE PROCEDURE READ_SURFT0_ASC
        MODULE PROCEDURE READ_SURFT1_ASC
        MODULE PROCEDURE READ_SURFT2_ASC
END INTERFACE
!
CONTAINS
!
!     #############################################################
      SUBROUTINE READ_SURFX0_ASC(HREC,PFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READX0* - routine to read a real scalar
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READX0 is
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      S.Malardel      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     01/08/03
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODE_POS_SURF
USE MODI_ERROR_READ_SURF_ASC
!
USE MODD_IO_SURF_ASC,        ONLY : NUNIT, NLUOUT, CMASK
USE MODD_SURF_PAR
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
CHARACTER(LEN=16), INTENT(IN)  :: HREC     ! name of the article to be read
REAL,              INTENT(OUT) :: PFIELD   ! the real scalar to be read
INTEGER,           INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),INTENT(OUT) :: HCOMMENT ! comment
!
!*      0.2   Declarations of local variables
!
LOGICAL :: GFOUND
CHARACTER(LEN=50):: YCOMMENT
CHARACTER(LEN=6) :: YMASK
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFX0_ASC',0,ZHOOK_HANDLE)
KRESP=0
!
YMASK=CMASK
CALL IO_BUFF_n(HREC,'R',GKNOWN)
IF (GKNOWN) YMASK='FULL  '
!
CALL POSNAM(NUNIT,YMASK//' '//HREC,GFOUND,NLUOUT)
IF (.NOT. GFOUND) CALL POSNAM(NUNIT,'FULL  '//' '//HREC,GFOUND,NLUOUT) ! used for auxilliary files

READ(NUNIT,FMT=*,END=100)
READ(NUNIT,FMT='(A50)') YCOMMENT
READ(NUNIT,FMT=*,ERR=100) PFIELD

HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFX0_ASC',1,ZHOOK_HANDLE)
RETURN
100 CONTINUE
CALL ERROR_READ_SURF_ASC(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFX0_ASC',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFX0_ASC
!
!     #############################################################
      SUBROUTINE READ_SURFX1_ASC(HREC,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READX1* - routine to fill a real 1D array for the externalised surface 
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READ_SURFX1 is
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      S.Malardel      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     01/08/03
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODE_POS_SURF
USE MODI_ERROR_READ_SURF_ASC
USE MODI_PACK_SAME_RANK
!
USE MODD_IO_SURF_ASC,        ONLY : NUNIT, NLUOUT, NMASK, NFULL, CMASK
USE MODD_SURF_PAR
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
CHARACTER(LEN=16),   INTENT(IN)  :: HREC     ! name of the article to be read
REAL, DIMENSION(:), INTENT(OUT) :: PFIELD   ! array containing the data field
INTEGER,             INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),  INTENT(OUT) :: HCOMMENT ! comment
CHARACTER(LEN=1),    INTENT(IN)  :: HDIR     ! type of field :
                                             ! 'H' : field with
                                             !       horizontal spatial dim.
                                             ! '-' : no horizontal dim.

!
!*      0.2   Declarations of local variables
!
LOGICAL                           :: GFOUND
CHARACTER(LEN=50):: YCOMMENT
REAL, DIMENSION(:),   ALLOCATABLE :: ZWORK   ! work array read in the file

CHARACTER(LEN=6) :: YMASK
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFX1_ASC',0,ZHOOK_HANDLE)
KRESP=0
!
!
IF (HDIR=='A') THEN
  CALL POSNAM(NUNIT,CMASK//' '//HREC,GFOUND,NLUOUT)
  IF (.NOT. GFOUND) CALL POSNAM(NUNIT,'FULL  '//' '//HREC,GFOUND,NLUOUT)
ELSE
  YMASK=CMASK
  CALL IO_BUFF_n(HREC,'R',GKNOWN)
  IF (GKNOWN) YMASK='FULL  '
!
  CALL POSNAM(NUNIT,YMASK//' '//HREC,GFOUND,NLUOUT)
END IF

IF (HDIR=='-' .OR. HDIR=='A') THEN
  ALLOCATE(ZWORK(SIZE(PFIELD)))
ELSE
  ALLOCATE(ZWORK(NFULL))
END IF

READ(NUNIT,FMT=*,END=100)
READ(NUNIT,FMT='(A50)') YCOMMENT
READ(NUNIT,FMT=*,ERR=100) ZWORK

HCOMMENT = YCOMMENT
!
IF (HDIR=='-' .OR. HDIR=='A') THEN
  PFIELD = ZWORK
ELSE
  CALL PACK_SAME_RANK(NMASK,ZWORK(:),PFIELD)
ENDIF
!  
DEALLOCATE(ZWORK)
!  
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFX1_ASC',1,ZHOOK_HANDLE)
RETURN
100 CONTINUE
CALL ERROR_READ_SURF_ASC(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFX1_ASC',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFX1_ASC
!
!     #############################################################
      SUBROUTINE READ_SURFX2_ASC(HREC,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READX2* - routine to fill a real 2D array for the externalised surface 
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READ_SURFX2 is
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      S.Malardel      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     01/08/03
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODE_POS_SURF
USE MODI_ERROR_READ_SURF_ASC
USE MODI_PACK_SAME_RANK
!
USE MODD_IO_SURF_ASC,        ONLY : NUNIT, NLUOUT, NMASK, NFULL, CMASK
USE MODD_SURF_PAR
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
CHARACTER(LEN=16),        INTENT(IN)  :: HREC     ! name of the article to be read
REAL, DIMENSION(:,:), INTENT(OUT) :: PFIELD   ! array containing the data field
INTEGER,                  INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),       INTENT(OUT) :: HCOMMENT ! comment
CHARACTER(LEN=1),         INTENT(IN)  :: HDIR     ! type of field :
                                                  ! 'H' : field with
                                                  !       horizontal spatial dim.
                                                  ! '-' : no horizontal dim.
!
!*      0.2   Declarations of local variables
! 
LOGICAL                           :: GFOUND
CHARACTER(LEN=50):: YCOMMENT
REAL, DIMENSION(:,:), ALLOCATABLE :: ZWORK   ! work array read in the file
CHARACTER(LEN=6) :: YMASK
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFX2_ASC',0,ZHOOK_HANDLE)
KRESP=0
!
IF (HDIR=='A') THEN
  CALL POSNAM(NUNIT,CMASK//' '//HREC,GFOUND,NLUOUT)
  IF (.NOT. GFOUND) CALL POSNAM(NUNIT,'FULL  '//' '//HREC,GFOUND,NLUOUT)
ELSE
  YMASK=CMASK
  CALL IO_BUFF_n(HREC,'R',GKNOWN)
  IF (GKNOWN) YMASK='FULL  '
  !
  CALL POSNAM(NUNIT,YMASK//' '//HREC,GFOUND,NLUOUT)
END IF

IF (HDIR=='-' .OR. HDIR=='A') THEN
  ALLOCATE(ZWORK(SIZE(PFIELD,1),SIZE(PFIELD,2)))
ELSE
  ALLOCATE(ZWORK(NFULL,SIZE(PFIELD,2)))
END IF

READ(NUNIT,FMT=*,END=100)
READ(NUNIT,FMT='(A50)') YCOMMENT
READ(NUNIT,FMT=*,ERR=100) ZWORK

HCOMMENT = YCOMMENT
!
IF (HDIR=='-' .OR. HDIR=='A') THEN
  PFIELD = ZWORK
ELSE
  CALL PACK_SAME_RANK(NMASK,ZWORK(:,:),PFIELD)
END IF

!  
DEALLOCATE(ZWORK)
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFX2_ASC',1,ZHOOK_HANDLE)
RETURN
100 CONTINUE
CALL ERROR_READ_SURF_ASC(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFX2_ASC',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFX2_ASC
!
!     #############################################################
      SUBROUTINE READ_SURFN0_ASC(HREC,KFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READN0* - routine to read an integer
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READN0 is
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      S.Malardel      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     01/08/03
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODE_POS_SURF
USE MODI_ERROR_READ_SURF_ASC
!
USE MODD_IO_SURF_ASC,        ONLY : NUNIT, NLUOUT, NMASK, CMASK
USE MODD_SURF_PAR
!
!
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
INTEGER,            INTENT(OUT) :: KFIELD   ! the integer to be read
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT ! comment
!
!*      0.2   Declarations of local variables
!
LOGICAL :: GFOUND
CHARACTER(LEN=50):: YCOMMENT
CHARACTER(LEN=6) :: YMASK
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFN0_ASC',0,ZHOOK_HANDLE)
KRESP=0
!
YMASK=CMASK
CALL IO_BUFF_n(HREC,'R',GKNOWN)
IF (GKNOWN) YMASK='FULL  '
!
CALL POSNAM(NUNIT,YMASK//' '//HREC,GFOUND,NLUOUT)
IF (.NOT. GFOUND) CALL POSNAM(NUNIT,'FULL  '//' '//HREC,GFOUND,NLUOUT) ! used for auxilliary files

READ(NUNIT,FMT=*,END=100)
READ(NUNIT,FMT='(A50)') YCOMMENT
READ(NUNIT,FMT=*,ERR=100) KFIELD

HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFN0_ASC',1,ZHOOK_HANDLE)
RETURN
100 CONTINUE
CALL ERROR_READ_SURF_ASC(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFN0_ASC',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFN0_ASC
!
!     #############################################################
      SUBROUTINE READ_SURFN1_ASC(HREC,KFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READN0* - routine to read an integer
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READN0 is
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      S.Malardel      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     01/08/03
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODE_POS_SURF
USE MODI_ERROR_READ_SURF_ASC
USE MODI_PACK_SAME_RANK
!
USE MODD_IO_SURF_ASC,        ONLY : NUNIT, NLUOUT, NMASK, NFULL, CMASK
USE MODD_SURF_PAR
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
CHARACTER(LEN=16),      INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER, DIMENSION(:), INTENT(OUT) :: KFIELD   ! the integer to be read
INTEGER,                INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),     INTENT(OUT) :: HCOMMENT ! comment
CHARACTER(LEN=1),       INTENT(IN)  :: HDIR     ! type of field :
                                                ! 'H' : field with
                                                !       horizontal spatial dim.
                                                ! '-' : no horizontal dim.
!
!*      0.2   Declarations of local variables
!
LOGICAL                            :: GFOUND
CHARACTER(LEN=50):: YCOMMENT
INTEGER, DIMENSION(:), ALLOCATABLE :: IWORK  ! work array read in the file
CHARACTER(LEN=6) :: YMASK
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!---------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFN1_ASC',0,ZHOOK_HANDLE)
KRESP = 0
!
YMASK=CMASK
CALL IO_BUFF_n(HREC,'R',GKNOWN)
IF (GKNOWN) YMASK='FULL  '
!
IF (HDIR=="-" .OR. HDIR=='A') THEN
  CALL POSNAM(NUNIT,YMASK//' '//HREC,GFOUND,NLUOUT)
!  IF (.NOT. GFOUND) CALL POSNAM(NUNIT,'FULL  '//' '//HREC,GFOUND,NLUOUT) ! used for auxilliary files
  READ(NUNIT,FMT=*,END=100)
  READ(NUNIT,FMT='(A50)') YCOMMENT
  READ(NUNIT,FMT=*,ERR=100) KFIELD

ELSE
  CALL POSNAM(NUNIT,YMASK//' '//HREC,GFOUND,NLUOUT)        

  ALLOCATE(IWORK(NFULL))

  READ(NUNIT,FMT=*,END=100)
  READ(NUNIT,FMT='(A50)') YCOMMENT
  READ(NUNIT,FMT=*,ERR=100) IWORK

  !
  CALL PACK_SAME_RANK(NMASK,IWORK(:),KFIELD)
  !  
  DEALLOCATE(IWORK)
ENDIF

HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFN1_ASC',1,ZHOOK_HANDLE)
RETURN
100 CONTINUE
CALL ERROR_READ_SURF_ASC(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFN1_ASC',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFN1_ASC
!
!     #############################################################
      SUBROUTINE READ_SURFC0_ASC(HREC,HFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READC0* - routine to read a character
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READC0 is
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      S.Malardel      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     01/08/03
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODE_POS_SURF
USE MODI_ERROR_READ_SURF_ASC
!
USE MODD_IO_SURF_ASC,        ONLY : NUNIT, NLUOUT, CMASK
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
CHARACTER(LEN=40),  INTENT(OUT) :: HFIELD    ! the integer to be read
INTEGER,            INTENT(OUT) :: KRESP     ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT  ! comment
!
!*      0.2   Declarations of local variables
!
LOGICAL           :: GFOUND
CHARACTER(LEN=50):: YCOMMENT
CHARACTER(LEN=6) :: YMASK
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFC0_ASC',0,ZHOOK_HANDLE)
KRESP=0
!
YMASK=CMASK
CALL IO_BUFF_n(HREC,'R',GKNOWN)
IF (GKNOWN) YMASK='FULL  '
!
CALL POSNAM(NUNIT,YMASK//' '//HREC,GFOUND,NLUOUT)
IF (.NOT. GFOUND) CALL POSNAM(NUNIT,'FULL  '//' '//HREC,GFOUND,NLUOUT) ! used for auxilliary files
!
READ(NUNIT,FMT=*,END=100)
READ(NUNIT,FMT='(A50)') YCOMMENT
READ(NUNIT,FMT='(A40)',ERR=100) HFIELD

HCOMMENT = YCOMMENT
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFC0_ASC',1,ZHOOK_HANDLE)
RETURN
100 CONTINUE
CALL ERROR_READ_SURF_ASC(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFC0_ASC',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFC0_ASC
!
!     #############################################################
      SUBROUTINE READ_SURFL1_ASC(HREC,OFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READL1* - routine to read a logical array
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READL1 is
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      S.Malardel      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     01/08/03
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODE_POS_SURF
USE MODI_ERROR_READ_SURF_ASC
!
USE MODD_IO_SURF_ASC,        ONLY : NUNIT, NLUOUT, CMASK
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
CHARACTER(LEN=16),      INTENT(IN)  :: HREC     ! name of the article to be read
LOGICAL, DIMENSION(:), INTENT(OUT) :: OFIELD   ! array containing the data field
INTEGER,                INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),     INTENT(OUT) :: HCOMMENT ! comment
CHARACTER(LEN=1),       INTENT(IN)  :: HDIR     ! type of field :
                                                ! 'H' : field with
                                                !       horizontal spatial dim.
                                                ! '-' : no horizontal dim.
!
!*      0.2   Declarations of local variables
!
LOGICAL           :: GFOUND
CHARACTER(LEN=50):: YCOMMENT
CHARACTER(LEN=6) :: YMASK
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFL1_ASC',0,ZHOOK_HANDLE)
KRESP=0
!
YMASK=CMASK
CALL IO_BUFF_n(HREC,'R',GKNOWN)
IF (GKNOWN) YMASK='FULL  '
!
CALL POSNAM(NUNIT,YMASK//' '//HREC,GFOUND,NLUOUT)
IF (.NOT. GFOUND) CALL POSNAM(NUNIT,'FULL  '//' '//HREC,GFOUND,NLUOUT) ! used for auxilliary files
!
READ(NUNIT,FMT=*,END=100)
READ(NUNIT,FMT='(A50)') YCOMMENT
READ(NUNIT,FMT=*,ERR=100) OFIELD

HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFL1_ASC',1,ZHOOK_HANDLE)
RETURN
100 CONTINUE
CALL ERROR_READ_SURF_ASC(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFL1_ASC',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFL1_ASC
!
!
!     #############################################################
      SUBROUTINE READ_SURFL0_ASC(HREC,OFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READL0* - routine to read a logical
!!
!!    PURPOSE
!!    -------
!
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      S.Malardel      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     01/08/03
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODE_POS_SURF
USE MODI_ERROR_READ_SURF_ASC
!
USE MODD_IO_SURF_ASC,        ONLY : NUNIT, NLUOUT, CMASK
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
LOGICAL,            INTENT(OUT) :: OFIELD   ! array containing the data field
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT ! comment
!
!*      0.2   Declarations of local variables
!
LOGICAL           :: GFOUND
CHARACTER(LEN=50):: YCOMMENT
CHARACTER(LEN=6) :: YMASK
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFL0_ASC',0,ZHOOK_HANDLE)
KRESP=0
!
YMASK=CMASK
CALL IO_BUFF_n(HREC,'R',GKNOWN)
IF (GKNOWN) YMASK='FULL  '
!
CALL POSNAM(NUNIT,YMASK//' '//HREC,GFOUND,NLUOUT)
IF (.NOT. GFOUND) CALL POSNAM(NUNIT,'FULL  '//' '//HREC,GFOUND,NLUOUT) ! used for auxilliary files
!
READ(NUNIT,FMT=*,END=100)
READ(NUNIT,FMT='(A50)') YCOMMENT
READ(NUNIT,FMT=*,ERR=100) OFIELD

HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFL0_ASC',1,ZHOOK_HANDLE)
RETURN
100 CONTINUE
CALL ERROR_READ_SURF_ASC(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFL0_ASC',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFL0_ASC
!
!     #############################################################
      SUBROUTINE READ_SURFT0_ASC(HREC,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READT0* - routine to read a date
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READT0 is
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      V. MASSON      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     18/08/97
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODE_POS_SURF
USE MODI_ERROR_READ_SURF_ASC
!
USE MODD_IO_SURF_ASC,        ONLY : NUNIT, NLUOUT, CMASK
USE MODD_SURF_PAR
!
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
INTEGER,            INTENT(OUT) :: KYEAR    ! year
INTEGER,            INTENT(OUT) :: KMONTH   ! month
INTEGER,            INTENT(OUT) :: KDAY     ! day
REAL,               INTENT(OUT) :: PTIME    ! year
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT ! comment

!*      0.2   Declarations of local variables
!
LOGICAL               :: GFOUND
CHARACTER(LEN=50):: YCOMMENT
INTEGER, DIMENSION(3) :: ITDATE
CHARACTER(LEN=6) :: YMASK
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFT0_ASC',0,ZHOOK_HANDLE)
KRESP=0
!
YMASK=CMASK
CALL IO_BUFF_n(HREC,'R',GKNOWN)
IF (GKNOWN) YMASK='FULL  '
!
!-------------------------------------------------------------------------------
!
CALL POSNAM(NUNIT,YMASK//' '//TRIM(HREC)//'%TDATE',GFOUND,NLUOUT)
!IF (.NOT. GFOUND) CALL POSNAM(NUNIT,'FULL  '//' '//HREC,GFOUND,NLUOUT) ! used for auxilliary files

READ(NUNIT,FMT=*,END=100)
READ(NUNIT,FMT='(A50)') YCOMMENT
READ(NUNIT,FMT=*,ERR=100) ITDATE(:)

KYEAR  = ITDATE(1)
KMONTH = ITDATE(2)
KDAY   = ITDATE(3)

!-------------------------------------------------------------------------------
!
CALL POSNAM(NUNIT,YMASK//' '//TRIM(HREC)//'%TIME',GFOUND,NLUOUT)
IF (.NOT. GFOUND) CALL POSNAM(NUNIT,'FULL  '//' '//HREC,GFOUND,NLUOUT) ! used for auxilliary files

READ(NUNIT,FMT=*,END=100)
READ(NUNIT,FMT='(A50)') YCOMMENT
READ(NUNIT,FMT=*,ERR=100) PTIME

!-------------------------------------------------------------------------------
!
HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFT0_ASC',1,ZHOOK_HANDLE)
RETURN
100 CONTINUE
CALL ERROR_READ_SURF_ASC(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFT0_ASC',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_SURFT0_ASC
!
!     #############################################################
      SUBROUTINE READ_SURFT1_ASC(HREC,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READT2* - routine to read a date
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READT2 is
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      V. MASSON      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     18/08/97
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODE_POS_SURF
USE MODI_ERROR_READ_SURF_ASC
!
USE MODD_IO_SURF_ASC,        ONLY : NUNIT, NLUOUT, CMASK
USE MODD_SURF_PAR
!
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
INTEGER, DIMENSION(:), INTENT(OUT) :: KYEAR    ! year
INTEGER, DIMENSION(:), INTENT(OUT) :: KMONTH   ! month
INTEGER, DIMENSION(:), INTENT(OUT) :: KDAY     ! day
REAL,    DIMENSION(:), INTENT(OUT) :: PTIME    ! year
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT ! comment

!*      0.2   Declarations of local variables
!
LOGICAL               :: GFOUND
CHARACTER(LEN=50):: YCOMMENT
INTEGER, DIMENSION(3,SIZE(KYEAR)) :: ITDATE
CHARACTER(LEN=6) :: YMASK
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFT1_ASC',0,ZHOOK_HANDLE)
KRESP=0
!
YMASK=CMASK
CALL IO_BUFF_n(HREC,'R',GKNOWN)
IF (GKNOWN) YMASK='FULL  '
!
!-------------------------------------------------------------------------------
!
CALL POSNAM(NUNIT,YMASK//' '//TRIM(HREC)//'%TDATE',GFOUND,NLUOUT)
!-------------------------------------------------------------------------------

READ(NUNIT,FMT=*,END=100)
READ(NUNIT,FMT='(A50)') YCOMMENT
READ(NUNIT,FMT=*,ERR=100) ITDATE(:,:)

KYEAR  (:) = ITDATE(1,:)
KMONTH (:) = ITDATE(2,:)
KDAY   (:) = ITDATE(3,:)

!-------------------------------------------------------------------------------
!
CALL POSNAM(NUNIT,CMASK//' '//TRIM(HREC)//'%TIME',GFOUND,NLUOUT)

READ(NUNIT,FMT=*,END=100)
READ(NUNIT,FMT='(A50)') YCOMMENT
READ(NUNIT,FMT=*,ERR=100) PTIME

!-------------------------------------------------------------------------------
!
HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFT1_ASC',1,ZHOOK_HANDLE)
RETURN
100 CONTINUE
CALL ERROR_READ_SURF_ASC(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFT1_ASC',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_SURFT1_ASC
!
!
!     #############################################################
      SUBROUTINE READ_SURFT2_ASC(HREC,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READT2* - routine to read a date
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READT2 is
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!
!!      V. MASSON      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     18/08/97
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODE_POS_SURF
USE MODI_ERROR_READ_SURF_ASC
!
USE MODD_IO_SURF_ASC,        ONLY : NUNIT, NLUOUT, CMASK
USE MODD_SURF_PAR
!
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
INTEGER, DIMENSION(:,:), INTENT(OUT) :: KYEAR    ! year
INTEGER, DIMENSION(:,:), INTENT(OUT) :: KMONTH   ! month
INTEGER, DIMENSION(:,:), INTENT(OUT) :: KDAY     ! day
REAL,    DIMENSION(:,:), INTENT(OUT) :: PTIME    ! year
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT ! comment

!*      0.2   Declarations of local variables
!
LOGICAL               :: GFOUND
CHARACTER(LEN=50):: YCOMMENT
INTEGER, DIMENSION(3,SIZE(KYEAR,1),SIZE(KYEAR,2)) :: ITDATE
CHARACTER(LEN=6) :: YMASK
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFT2_ASC',0,ZHOOK_HANDLE)
KRESP=0
!
YMASK=CMASK
CALL IO_BUFF_n(HREC,'R',GKNOWN)
IF (GKNOWN) YMASK='FULL  '
!
!-------------------------------------------------------------------------------
!
CALL POSNAM(NUNIT,YMASK//' '//TRIM(HREC)//'%TDATE',GFOUND,NLUOUT)

READ(NUNIT,FMT=*,END=100)
READ(NUNIT,FMT='(A50)') YCOMMENT
READ(NUNIT,FMT=*,ERR=100) ITDATE(:,:,:)

KYEAR  (:,:) = ITDATE(1,:,:)
KMONTH (:,:) = ITDATE(2,:,:)
KDAY   (:,:) = ITDATE(3,:,:)

!-------------------------------------------------------------------------------
!
CALL POSNAM(NUNIT,YMASK//' '//TRIM(HREC)//'%TIME',GFOUND,NLUOUT)

READ(NUNIT,FMT=*,END=100)
READ(NUNIT,FMT='(A50)') YCOMMENT
READ(NUNIT,FMT=*,ERR=100) PTIME

!-------------------------------------------------------------------------------
!
HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFT2_ASC',1,ZHOOK_HANDLE)
RETURN
100 CONTINUE
CALL ERROR_READ_SURF_ASC(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_ASC:READ_SURFT2_ASC',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_SURFT2_ASC
!
END MODULE MODE_READ_SURF_ASC
