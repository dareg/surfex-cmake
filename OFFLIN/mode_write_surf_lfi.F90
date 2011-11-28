!     ######spl
MODULE MODE_WRITE_SURF_LFI
!
USE MODI_GET_LUOUT
INTERFACE WRITE_SURF0_LFI
        MODULE PROCEDURE WRITE_SURFX0_LFI
        MODULE PROCEDURE WRITE_SURFN0_LFI
        MODULE PROCEDURE WRITE_SURFL0_LFI
        MODULE PROCEDURE WRITE_SURFC0_LFI
END INTERFACE
INTERFACE WRITE_SURFN_LFI
        MODULE PROCEDURE WRITE_SURFX1_LFI
        MODULE PROCEDURE WRITE_SURFN1_LFI
        MODULE PROCEDURE WRITE_SURFL1_LFI
        MODULE PROCEDURE WRITE_SURFX2_LFI
END INTERFACE
INTERFACE WRITE_SURFT_LFI
        MODULE PROCEDURE WRITE_SURFT0_LFI
END INTERFACE
!
CONTAINS
!
!     #############################################################
      SUBROUTINE WRITE_SURFX0_LFI(HREC,PFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  * - routine to write a real scalar
!
USE MODI_FMWRIT
USE MODI_ERROR_WRITE_SURF_LFI
!
USE MODD_IO_SURF_LFI,        ONLY : CFILEOUT_LFI, CLUOUT_LFI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_GET_SURF_UNDEF
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
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFX0_LFI',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',GKNOWN)
IF (GKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFX0_LFI',1,ZHOOK_HANDLE)
IF (GKNOWN) RETURN
!
CALL FMWRIT(CFILEOUT_LFI,HREC,CLUOUT_LFI,1,PFIELD,4,100,HCOMMENT,KRESP)
!
CALL ERROR_WRITE_SURF_LFI(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFX0_LFI',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFX0_LFI
!
!     #############################################################
      SUBROUTINE WRITE_SURFN0_LFI(HREC,KFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  * - routine to write an integer
!
USE MODI_FMWRIT
USE MODI_ERROR_WRITE_SURF_LFI
!
USE MODD_IO_SURF_LFI,        ONLY : CFILEOUT_LFI, CLUOUT_LFI, &
                                      LMNH_COMPATIBLE, NIU, NIB, NIE, NJU, NJB, NJE  
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_GET_SURF_UNDEF
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
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFN0_LFI',0,ZHOOK_HANDLE)
KRESP=0
!
IF (LMNH_COMPATIBLE .AND. HREC=='IMAX') THEN
  NIU = KFIELD+2
  NIB = 2
  NIE = KFIELD+1
END IF
IF (LMNH_COMPATIBLE .AND. HREC=='JMAX') THEN
  NJU = KFIELD+2
  NJB = 2
  NJE = KFIELD+1
END IF
!
CALL IO_BUFF_n(HREC,'W',GKNOWN)
IF (GKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFN0_LFI',1,ZHOOK_HANDLE)
IF (GKNOWN) RETURN
!
CALL FMWRIT(CFILEOUT_LFI,HREC,CLUOUT_LFI,1,KFIELD,4,100,HCOMMENT,KRESP)
!
CALL ERROR_WRITE_SURF_LFI(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFN0_LFI',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFN0_LFI
!
!     #############################################################
      SUBROUTINE WRITE_SURFL0_LFI(HREC,OFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  * - routine to write a logical
!
USE MODI_FMWRIT
USE MODI_ERROR_WRITE_SURF_LFI
!
USE MODD_IO_SURF_LFI,        ONLY : CFILEOUT_LFI, CLUOUT_LFI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_GET_SURF_UNDEF
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
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFL0_LFI',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',GKNOWN)
IF (GKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFL0_LFI',1,ZHOOK_HANDLE)
IF (GKNOWN) RETURN
!
CALL FMWRIT(CFILEOUT_LFI,HREC,CLUOUT_LFI,1,OFIELD,4,100,HCOMMENT,KRESP)
!
CALL ERROR_WRITE_SURF_LFI(HREC,KRESP)
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFL0_LFI',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFL0_LFI
!
!     #############################################################
      SUBROUTINE WRITE_SURFC0_LFI(HREC,HFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  * - routine to write a character
!
USE MODI_FMWRIT
USE MODI_ERROR_WRITE_SURF_LFI
!
USE MODD_IO_SURF_LFI,        ONLY : CFILEOUT_LFI, CLUOUT_LFI, LMNH_COMPATIBLE, LCARTESIAN
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_GET_SURF_UNDEF
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
LOGICAL          :: GCARTESIAN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFC0_LFI',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',GKNOWN)
IF (GKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFC0_LFI',1,ZHOOK_HANDLE)
IF (GKNOWN) RETURN
!
CALL FMWRIT(CFILEOUT_LFI,HREC,CLUOUT_LFI,1,HFIELD,4,100,HCOMMENT,KRESP)
!
IF (HREC=="GRID_TYPE") LMNH_COMPATIBLE = (HFIELD=="CARTESIAN " .OR. HFIELD=="CONF PROJ ")
IF (HREC=="GRID_TYPE" .AND. LMNH_COMPATIBLE) LCARTESIAN=(HFIELD=="CARTESIAN ")
!
CALL ERROR_WRITE_SURF_LFI(HREC,KRESP)
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFC0_LFI',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFC0_LFI
!
!     #############################################################
      SUBROUTINE WRITE_SURFX1_LFI(HREC,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  * - routine to fill a write 1D array for the externalised surface 
!
USE MODI_FMWRIT
USE MODI_ERROR_WRITE_SURF_LFI
USE MODI_UNPACK_SAME_RANK
!
USE MODD_IO_SURF_LFI,        ONLY : CFILEOUT_LFI, CLUOUT_LFI, NMASK, NFULL, &
                                    LMNH_COMPATIBLE, NIU, NIB, NIE, NJU, NJB, NJE
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_GET_SURF_UNDEF
USE MODI_IO_BUFF_n
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=16),   INTENT(IN) :: HREC     ! name of the article to be read
REAL, DIMENSION(:),  INTENT(IN) :: PFIELD   ! array containing the data field
INTEGER,             INTENT(OUT):: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),  INTENT(IN) :: HCOMMENT ! comment string
CHARACTER(LEN=1),    INTENT(IN) :: HDIR     ! type of field :
                                            ! 'H' : field with
                                            !       horizontal spatial dim.
                                            ! '-' : no horizontal dim.
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=20)        :: YREC
REAL, DIMENSION(NFULL)   :: ZWORK   ! work array read in the file
REAL, DIMENSION(NIU,NJU) :: ZWORK2D ! work array read in a MNH file
REAL                     :: ZUNDEF  ! default value
LOGICAL                  :: GKNOWN
INTEGER                  :: JI, JJ
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFX1_LFI',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',GKNOWN)
IF (GKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFX1_LFI',1,ZHOOK_HANDLE)
IF (GKNOWN) RETURN
!
IF (HDIR=='H') THEN
  !
  CALL UNPACK_SAME_RANK(NMASK,PFIELD,ZWORK(:))
  CALL GET_SURF_UNDEF(ZUNDEF)
  WHERE(ZWORK==ZUNDEF) ZWORK=999.
  !
  IF (.NOT. LMNH_COMPATIBLE) THEN
    CALL FMWRIT(CFILEOUT_LFI,HREC,CLUOUT_LFI,NFULL,ZWORK,4,100,HCOMMENT,KRESP)
  ELSE
    !
    ZWORK2D(:,:) = 999.
    DO JJ=1,NJE-NJB+1
      DO JI=1,NIE-NIB+1
        ZWORK2D(NIB+JI-1,NJB+JJ-1) = ZWORK(JI+(NIE-NIB+1)*(JJ-1))
      END DO
    END DO
    !
    IF     (HREC=='DX              ' .OR. HREC=='XX              ') THEN
      YREC = 'XHAT'
      CALL WRITE_IN_LFI_X1_FOR_MNH(HREC,YREC,ZWORK2D(NIB:NIE,NJB),KRESP,HCOMMENT,NIU,NIB,NIE)
    ELSEIF (HREC=='DY              ' .OR. HREC=='YY              ') THEN
      YREC = 'YHAT'
      CALL WRITE_IN_LFI_X1_FOR_MNH(HREC,YREC,ZWORK2D(NIB,NJB:NJE),KRESP,HCOMMENT,NJU,NJB,NJE)
    ELSE
      CALL FMWRIT(CFILEOUT_LFI,HREC,CLUOUT_LFI,SIZE(ZWORK2D),ZWORK2D,4,100,HCOMMENT,KRESP)
    ENDIF
    !
  END IF
  !
ELSE
  CALL FMWRIT(CFILEOUT_LFI,HREC,CLUOUT_LFI,SIZE(PFIELD),PFIELD,4,100,HCOMMENT,KRESP)
END IF
!
CALL ERROR_WRITE_SURF_LFI(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFX1_LFI',1,ZHOOK_HANDLE)
!
CONTAINS
!
!     #############################################################
      SUBROUTINE WRITE_IN_LFI_X1_FOR_MNH(HREC,HREC2,PFIELD,KRESP,HCOMMENT,KU,KB,KE)
!     #############################################################
!
!!****  * - routine to fill a write 2D array for the externalised surface 
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=16),        INTENT(IN) :: HREC     ! name of the article to be read
CHARACTER(LEN=20),        INTENT(IN) :: HREC2    ! name of the article to be read
REAL, DIMENSION(:),       INTENT(IN) :: PFIELD   ! array containing the data field
INTEGER,                  INTENT(OUT):: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),       INTENT(IN) :: HCOMMENT ! comment string
INTEGER,                  INTENT(IN) :: KU
INTEGER,                  INTENT(IN) :: KB
INTEGER,                  INTENT(IN) :: KE
!
!*      0.2   Declarations of local variables
! 
REAL, DIMENSION(KU)      :: ZWORK ! 1D work array read in the file
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFX1_LFI:WRITE_IN_LFI_X1_FOR_MNH',0,ZHOOK_HANDLE)
!
ZWORK(:) = 0.
!
SELECT CASE(HREC)
  !
  CASE('DX              ','DY              ')
    IF (KB/=KE) THEN
      IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFX1_LFI:WRITE_IN_LFI_X1_FOR_MNH',1,ZHOOK_HANDLE)
      RETURN
    ENDIF
    ZWORK(1) = - PFIELD(1)*0.5  ! 1D case
    ZWORK(2) =   PFIELD(1)*0.5
    ZWORK(3) =   PFIELD(1)*1.5
  !
  CASE('XX              ','YY              ')
    IF (KB==KE) THEN
      IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFX1_LFI:WRITE_IN_LFI_X1_FOR_MNH',1,ZHOOK_HANDLE)
      RETURN
    ENDIF          
    ZWORK(KB+1:KE)   = 0.5 * PFIELD(1:KE-2) + 0.5 * PFIELD(2:KE-1)
    ZWORK(KB)        = 1.5 * PFIELD(1)      - 0.5 * PFIELD(2)
    ZWORK(KB-1)      = 2. * ZWORK(KB) - ZWORK(KB+1)
    ZWORK(KE+1)      = 2. * ZWORK(KE) - ZWORK(KE-1)
  !  
END SELECT
!
CALL FMWRIT(CFILEOUT_LFI,HREC2,CLUOUT_LFI,KU,ZWORK,4,100,HCOMMENT,KRESP)
CALL ERROR_WRITE_SURF_LFI(HREC2,KRESP)
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFX1_LFI:WRITE_IN_LFI_X1_FOR_MNH',1,ZHOOK_HANDLE)
END SUBROUTINE WRITE_IN_LFI_X1_FOR_MNH
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFX1_LFI
!
!     #############################################################
      SUBROUTINE WRITE_SURFN1_LFI(HREC,KFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  * - routine to write an integer array
!
USE MODI_FMWRIT
USE MODI_ERROR_WRITE_SURF_LFI
USE MODI_UNPACK_SAME_RANK
!
USE MODD_IO_SURF_LFI,        ONLY : CFILEOUT_LFI, CLUOUT_LFI, NMASK, NFULL
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_GET_SURF_UNDEF
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
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFN1_LFI',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',GKNOWN)
IF (GKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFN1_LFI',1,ZHOOK_HANDLE)
IF (GKNOWN) RETURN
!
IF (HDIR=='H') THEN
  CALL UNPACK_SAME_RANK(NMASK,KFIELD,IWORK(:))
  CALL FMWRIT(CFILEOUT_LFI,HREC,CLUOUT_LFI,NFULL,IWORK,4,100,HCOMMENT,KRESP)
ELSE
  CALL FMWRIT(CFILEOUT_LFI,HREC,CLUOUT_LFI,SIZE(KFIELD),KFIELD,4,100,HCOMMENT,KRESP)
END IF
!
CALL ERROR_WRITE_SURF_LFI(HREC,KRESP)
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFN1_LFI',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFN1_LFI
!
!
!     #############################################################
      SUBROUTINE WRITE_SURFL1_LFI(HREC,OFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  * - routine to write a logical array
!
USE MODI_FMWRIT
USE MODI_ERROR_WRITE_SURF_LFI
!
USE MODD_IO_SURF_LFI,        ONLY : CFILEOUT_LFI, CLUOUT_LFI
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
USE MODI_GET_SURF_UNDEF
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
!
!*      0.2   Declarations of local variables
!
INTEGER         :: ILUOUT ! listing logical unit
LOGICAL         :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFL1_LFI',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',GKNOWN)
IF (GKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFL1_LFI',1,ZHOOK_HANDLE)
IF (GKNOWN) RETURN
!
IF (HDIR=='H') THEN
  CALL GET_LUOUT('LFI   ',ILUOUT)
  WRITE(ILUOUT,*) 'Error: 1D logical vector for writing on an horizontal grid:'
  WRITE(ILUOUT,*) 'this option is not coded in WRITE_SURFL1_LFI'
  CALL ABOR1_SFX('MODE_WRITE_SURF_LFI: 1D LOGICAL VECTOR FOR WRITING NOT CODED IN WRITE_SURFL1_LFI')
ELSE
  CALL FMWRIT(CFILEOUT_LFI,HREC,CLUOUT_LFI,SIZE(OFIELD),OFIELD,4,100,HCOMMENT,KRESP)
END IF
!
CALL ERROR_WRITE_SURF_LFI(HREC,KRESP)
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFL1_LFI',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFL1_LFI
!
!     #############################################################
      SUBROUTINE WRITE_SURFX2_LFI(HREC,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  * - routine to fill a write 2D array for the externalised surface 
!
USE MODI_FMWRIT
USE MODI_ERROR_WRITE_SURF_LFI
USE MODI_UNPACK_SAME_RANK
!
USE MODD_IO_SURF_LFI,        ONLY : CFILEOUT_LFI, CLUOUT_LFI, NMASK, NFULL, &
                                      LMNH_COMPATIBLE  
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_GET_SURF_UNDEF
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
REAL                                  :: ZUNDEF  ! default value
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFX2_LFI',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',GKNOWN)
IF (GKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFX2_LFI',1,ZHOOK_HANDLE)
IF (GKNOWN) RETURN
!
IF (HDIR=='H') THEN
  !
  CALL UNPACK_SAME_RANK(NMASK,PFIELD,ZWORK(:,:))
  CALL GET_SURF_UNDEF(ZUNDEF)
  WHERE(ZWORK==ZUNDEF) ZWORK=999.
  !
  IF (.NOT. LMNH_COMPATIBLE) THEN
    CALL FMWRIT(CFILEOUT_LFI,HREC,CLUOUT_LFI,SIZE(ZWORK),ZWORK,4,100,HCOMMENT,KRESP)
  ELSE
    CALL WRITE_IN_LFI_X2_FOR_MNH(HREC,ZWORK,KRESP,HCOMMENT)
  END IF
  !
ELSE
  CALL FMWRIT(CFILEOUT_LFI,HREC,CLUOUT_LFI,SIZE(PFIELD),PFIELD,4,100,HCOMMENT,KRESP)
END IF
!  
CALL ERROR_WRITE_SURF_LFI(HREC,KRESP)
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFX2_LFI',1,ZHOOK_HANDLE)
!
CONTAINS
!
!     #############################################################
      SUBROUTINE WRITE_IN_LFI_X2_FOR_MNH(HREC,PFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  * - routine to fill a write 2D array for the externalised surface 
!
USE MODI_FMWRIT
USE MODI_ERROR_WRITE_SURF_LFI
!
USE MODD_IO_SURF_LFI,        ONLY : CFILEOUT_LFI, CLUOUT_LFI, &
                                    NIU, NIB, NIE, NJU, NJB, NJE  
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=16),        INTENT(IN) :: HREC     ! name of the article to be read
REAL, DIMENSION(:,:),     INTENT(IN) :: PFIELD   ! array containing the data field
INTEGER,                  INTENT(OUT):: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),       INTENT(IN) :: HCOMMENT ! comment string
!
!*      0.2   Declarations of local variables
! 
REAL, DIMENSION(NIU,NJU,SIZE(PFIELD,2)) :: ZWORK3D ! work array read in a MNH file
INTEGER :: JI, JJ
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFX2_LFI:WRITE_IN_LFI_X2_FOR_MNH',0,ZHOOK_HANDLE)
!
ZWORK3D(:,:,:) = 999.
DO JJ=1,NJE-NJB+1
  DO JI=1,NIE-NIB+1
    ZWORK3D(NIB+JI-1,NJB+JJ-1,:) = PFIELD(JI+(NIE-NIB+1)*(JJ-1),:)
  END DO
END DO
!
CALL FMWRIT(CFILEOUT_LFI,HREC,CLUOUT_LFI,SIZE(ZWORK3D),ZWORK3D,4,100,HCOMMENT,KRESP)
!  
CALL ERROR_WRITE_SURF_LFI(HREC,KRESP)
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFX2_LFI:WRITE_IN_LFI_X2_FOR_MNH',1,ZHOOK_HANDLE)
END SUBROUTINE WRITE_IN_LFI_X2_FOR_MNH
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFX2_LFI
!
!     #############################################################
      SUBROUTINE WRITE_SURFT0_LFI(HREC,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
!     #############################################################
!
!!****  * - routine to write a date
!
USE MODI_FMWRIT
USE MODI_ERROR_WRITE_SURF_LFI
!
USE MODD_IO_SURF_LFI,        ONLY : CFILEOUT_LFI, CLUOUT_LFI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_GET_SURF_UNDEF
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
CHARACTER(LEN=16)     :: YREC     ! Name of the article to be written
INTEGER, DIMENSION(3) :: ITDATE
INTEGER               :: IRET
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFT0_LFI',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',GKNOWN)
IF (GKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFT0_LFI',1,ZHOOK_HANDLE)
IF (GKNOWN) RETURN
!
ITDATE(1) = KYEAR
ITDATE(2) = KMONTH
ITDATE(3) = KDAY
!
YREC=TRIM(HREC)//'%TDATE'
CALL FMWRIT(CFILEOUT_LFI,YREC,CLUOUT_LFI,3,ITDATE,4,100,HCOMMENT,KRESP)
CALL ERROR_WRITE_SURF_LFI(HREC,KRESP)
!
YREC=TRIM(HREC)//'%TIME'
CALL FMWRIT(CFILEOUT_LFI,YREC,CLUOUT_LFI,1,PTIME,4,100,HCOMMENT,KRESP)
CALL ERROR_WRITE_SURF_LFI(HREC,KRESP)
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_LFI:WRITE_SURFT0_LFI',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFT0_LFI
!
END MODULE MODE_WRITE_SURF_LFI
