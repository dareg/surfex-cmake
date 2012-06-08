MODULE MODE_WRITE_SURF_FA
!
INTERFACE WRITE_SURF0_FA
        MODULE PROCEDURE WRITE_SURFX0_FA
        MODULE PROCEDURE WRITE_SURFN0_FA
        MODULE PROCEDURE WRITE_SURFL0_FA
        MODULE PROCEDURE WRITE_SURFC0_FA
END INTERFACE
INTERFACE WRITE_SURFN_FA
        MODULE PROCEDURE WRITE_SURFX1_FA
        MODULE PROCEDURE WRITE_SURFN1_FA
        MODULE PROCEDURE WRITE_SURFL1_FA
        MODULE PROCEDURE WRITE_SURFX2_FA
END INTERFACE
INTERFACE WRITE_SURFT_FA
        MODULE PROCEDURE WRITE_SURFT0_FA
        MODULE PROCEDURE WRITE_SURFT2_FA
END INTERFACE
!
CONTAINS
!
!     #############################################################
      SUBROUTINE WRITE_SURFX0_FA(HREC,PFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  * - routine to write a real scalar
!
USE MODI_ERROR_WRITE_SURF_FA
!
USE MODD_IO_SURF_FA, ONLY : NUNIT_FA, CMASK, LFANOCOMPACT
!
USE MODE_FASURFEX
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
LOGICAL          :: LKNOWN
CHARACTER(LEN=16):: YNAME                  ! Field Name
INTEGER          :: INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL
REAL(KIND=JPRB)  :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:ERROR_WRITE_SURF_FA:WRITE_SURFX0_FA',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',LKNOWN)
IF (LKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFX0_FA',1,ZHOOK_HANDLE)
IF (LKNOWN) RETURN
!
IF(LFANOCOMPACT)THEN
  CALL FAVEUR(KRESP,NUNIT_FA,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
  ! -- Pour ecrire sans compactage
  CALL FAGOTE(KRESP,NUNIT_FA,-1,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
ENDIF
!
YNAME=TRIM(CMASK)//TRIM(HREC)
CALL  FAECR_R(KRESP,NUNIT_FA,YNAME,PFIELD)
IF (KRESP/=0) THEN
  CALL ERROR_WRITE_SURF_FA(HREC,KRESP)
ENDIF
!
IF(LFANOCOMPACT)THEN
  ! On remet la valeur par defaut 
  CALL FAGOTE(KRESP,NUNIT_FA,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFX0_FA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFX0_FA
!
!     #############################################################
      SUBROUTINE WRITE_SURFN0_FA(HREC,KFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  * - routine to write an integer
!
USE MODI_ERROR_WRITE_SURF_FA
!
USE MODE_FASURFEX
!
USE MODD_IO_SURF_FA, ONLY : NUNIT_FA, NMASK, CMASK, LFANOCOMPACT
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
LOGICAL          :: LKNOWN
CHARACTER(LEN=16):: YNAME                  ! Field Name
INTEGER          :: INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL
REAL(KIND=JPRB)  :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFN0_FA',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',LKNOWN)
IF (LKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFN0_FA',1,ZHOOK_HANDLE)
IF (LKNOWN) RETURN
!
IF(LFANOCOMPACT)THEN
  CALL FAVEUR(KRESP,NUNIT_FA,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
  ! -- Pour ecrire sans compactage
  CALL FAGOTE(KRESP,NUNIT_FA,-1,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
ENDIF
!
YNAME=TRIM(CMASK)//TRIM(HREC)
CALL  FAECR_I(KRESP,NUNIT_FA,YNAME,KFIELD)
IF (KRESP/=0) THEN
  CALL ERROR_WRITE_SURF_FA(HREC,KRESP)
ENDIF
!
IF(LFANOCOMPACT)THEN
  ! On remet la valeur par defaut 
  CALL FAGOTE(KRESP,NUNIT_FA,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFN0_FA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFN0_FA
!
!     #############################################################
      SUBROUTINE WRITE_SURFL0_FA(HREC,OFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  * - routine to write a logical
!
USE MODI_ERROR_WRITE_SURF_FA
!
USE MODD_IO_SURF_FA, ONLY : NUNIT_FA, CMASK, LFANOCOMPACT
!
USE MODE_FASURFEX
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
LOGICAL          :: LKNOWN
CHARACTER(LEN=16):: YNAME ! Field Name
INTEGER          :: INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL
REAL(KIND=JPRB)  :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFL0_FA',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',LKNOWN)
IF (LKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFL0_FA',1,ZHOOK_HANDLE)
IF (LKNOWN) RETURN
!
IF(LFANOCOMPACT)THEN
  CALL FAVEUR(KRESP,NUNIT_FA,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
  ! -- Pour ecrire sans compactage
  CALL FAGOTE(KRESP,NUNIT_FA,-1,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
ENDIF
!
YNAME=TRIM(CMASK)//TRIM(HREC)
CALL  FAECR_L(KRESP,NUNIT_FA,YNAME,OFIELD)
IF (KRESP/=0) THEN
  CALL ERROR_WRITE_SURF_FA(HREC,KRESP)
ENDIF
!
IF(LFANOCOMPACT)THEN
  ! On remet la valeur par defaut 
  CALL FAGOTE(KRESP,NUNIT_FA,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFL0_FA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFL0_FA
!
!     #############################################################
      SUBROUTINE WRITE_SURFC0_FA(HREC,HFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  * - routine to write a character
!
USE MODI_ERROR_WRITE_SURF_FA
!
USE MODD_IO_SURF_FA, ONLY : NUNIT_FA, CMASK, LFANOCOMPACT
!
USE MODE_FASURFEX
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
LOGICAL          :: LKNOWN
CHARACTER,DIMENSION(40)  :: YFIELD
CHARACTER(LEN=16)        :: YNAME ! Field Name
INTEGER                  :: INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL
REAL(KIND=JPRB)          :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFC0_FA',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',LKNOWN)
IF (LKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFC0_FA',1,ZHOOK_HANDLE)
IF (LKNOWN) RETURN
!
IF(LFANOCOMPACT)THEN
  CALL FAVEUR(KRESP,NUNIT_FA,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
  ! -- Pour ecrire sans compactage
  CALL FAGOTE(KRESP,NUNIT_FA,-1,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
ENDIF
!
READ(HFIELD,'(40A1)') YFIELD
YNAME=TRIM(CMASK)//TRIM(HREC)
CALL  FAECR_C(KRESP,NUNIT_FA,YNAME,40,YFIELD)
IF (KRESP/=0) THEN
  CALL ERROR_WRITE_SURF_FA(HREC,KRESP)
ENDIF
!
IF(LFANOCOMPACT)THEN
  ! On remet la valeur par defaut 
  CALL FAGOTE(KRESP,NUNIT_FA,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFC0_FA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFC0_FA
!
!     #############################################################
      SUBROUTINE WRITE_SURFX1_FA(HREC,KL,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  * - routine to fill a write 1D array for the externalised surface 
!
USE MODI_ERROR_WRITE_SURF_FA
USE MODI_UNPACK_SAME_RANK
!
USE MODD_SURF_PAR,    ONLY : XUNDEF
USE MODD_IO_SURF_FA,  ONLY : NUNIT_FA, NMASK, NFULL, CMASK, &
                             LFANOCOMPACT 
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
INTEGER,             INTENT(IN) :: KL       ! number of points
REAL, DIMENSION(KL), INTENT(IN) :: PFIELD   ! array containing the data field
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
REAL                   :: ZMEAN, ZCOUNT
LOGICAL                :: LKNOWN
INTEGER                :: I,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL
REAL(KIND=JPRB)        :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFX1_FA',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',LKNOWN)
IF (LKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFX1_FA',1,ZHOOK_HANDLE)
IF (LKNOWN) RETURN
!
IF(HDIR=='H')THEN
  CALL UNPACK_SAME_RANK(NMASK,PFIELD,ZWORK(:))
ELSE !no horizontal dim. case (not masked)
  ZWORK(1:KL)=PFIELD(1:KL)
  ZWORK(KL+1:NFULL)=SUM(PFIELD(1:KL))/REAL(KL)
ENDIF
!
IF(LFANOCOMPACT)THEN
  CALL FAVEUR(KRESP,NUNIT_FA,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
  ! -- Pour ecrire sans compactage
  CALL FAGOTE(KRESP,NUNIT_FA,-1,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
  CALL FAIENC(KRESP,NUNIT_FA,'S1D_',0,HREC,ZWORK,.FALSE.)
  IF (KRESP/=0) THEN
    CALL ERROR_WRITE_SURF_FA(HREC,KRESP)
  ENDIF  
  ! On remet la valeur par defaut 
  CALL FAGOTE(KRESP,NUNIT_FA,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
ELSE
  ZMEAN =0.0
  ZCOUNT=0.0
  DO I=1,NFULL
     IF(ZWORK(I)/=XUNDEF)THEN
        ZMEAN =ZMEAN+ZWORK(I)
        ZCOUNT=ZCOUNT+1.0
     ENDIF
  ENDDO
  IF (ZCOUNT.GT.0.0) ZMEAN=ZMEAN/ZCOUNT
  WHERE(ZWORK(:)==XUNDEF)ZWORK(:)=ZMEAN
  CALL FAIENC(KRESP,NUNIT_FA,'S1D_',0,HREC,ZWORK,.FALSE.)
  IF (KRESP/=0) THEN
    CALL ERROR_WRITE_SURF_FA(HREC,KRESP)
  ENDIF
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFX1_FA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFX1_FA
!
!     #############################################################
      SUBROUTINE WRITE_SURFN1_FA(HREC,KL,KFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  * - routine to write an integer array
!
USE MODI_ERROR_WRITE_SURF_FA
USE MODI_UNPACK_SAME_RANK
!
USE MODD_IO_SURF_FA, ONLY : NUNIT_FA, NMASK, NFULL, CMASK, LFANOCOMPACT
!
USE MODE_FASURFEX
!
USE MODE_GRIDTYPE_GAUSS
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
INTEGER,                INTENT(IN) :: KL       ! number of points
INTEGER, DIMENSION(KL), INTENT(IN) :: KFIELD   ! array containing the data field
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
LOGICAL                   :: LKNOWN
CHARACTER(LEN=16)         :: YNAME ! Field Name
INTEGER                   :: INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL
REAL(KIND=JPRB)           :: ZHOOK_HANDLE
!---------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFN1_FA',0,ZHOOK_HANDLE)
KRESP = 0
!
CALL IO_BUFF_n(HREC,'W',LKNOWN)
IF (LKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFN1_FA',1,ZHOOK_HANDLE)
IF (LKNOWN) RETURN
!
IF(LFANOCOMPACT)THEN
  CALL FAVEUR(KRESP,NUNIT_FA,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
  ! -- Pour ecrire sans compactage
  CALL FAGOTE(KRESP,NUNIT_FA,-1,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
ENDIF
!
YNAME=TRIM(CMASK)//TRIM(HREC)
IF (HREC=="-") THEN
  CALL  FAECR_I_D(KRESP,NUNIT_FA,YNAME,KL,KFIELD)
  IF (KRESP/=0) THEN
    CALL ERROR_WRITE_SURF_FA(HREC,KRESP)
  ENDIF
ELSE
  IF (HDIR=='H') THEN
    CALL UNPACK_SAME_RANK(NMASK,KFIELD,IWORK(:))
  ELSE
    IWORK = KFIELD
  END IF
  CALL  FAECR_I_D(KRESP,NUNIT_FA,YNAME,KL,IWORK)
  IF (KRESP/=0) THEN
    CALL ERROR_WRITE_SURF_FA(HREC,KRESP)
  ENDIF
ENDIF
!
IF(LFANOCOMPACT)THEN
  ! On remet la valeur par defaut 
  CALL FAGOTE(KRESP,NUNIT_FA,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFN1_FA',1,ZHOOK_HANDLE)
RETURN
!
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFN1_FA
!
!
!     #############################################################
      SUBROUTINE WRITE_SURFL1_FA(HREC,KL,OFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  * - routine to write a logical array
!
USE MODI_ERROR_WRITE_SURF_FA
!
USE MODD_IO_SURF_FA, ONLY : NUNIT_FA, CMASK, LFANOCOMPACT
!
USE MODE_FASURFEX
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
INTEGER,             INTENT(IN) :: KL       ! number of points
LOGICAL, DIMENSION(KL), INTENT(IN) :: OFIELD   ! array containing the data field
INTEGER,                INTENT(OUT):: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),     INTENT(IN) :: HCOMMENT ! comment string
CHARACTER(LEN=1),       INTENT(IN) :: HDIR     ! type of field :
                                               ! 'H' : field with
                                               !       horizontal spatial dim.
                                               ! '-' : no horizontal dim.
!
!*      0.2   Declarations of local variables
!
LOGICAL          :: LKNOWN
CHARACTER(LEN=16):: YNAME ! Field Name
INTEGER          :: INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL
REAL(KIND=JPRB)  :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFL1_FA',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',LKNOWN)
IF (LKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFL1_FA',1,ZHOOK_HANDLE)
IF (LKNOWN) RETURN
!
IF(LFANOCOMPACT)THEN
  CALL FAVEUR(KRESP,NUNIT_FA,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
  ! -- Pour ecrire sans compactage
  CALL FAGOTE(KRESP,NUNIT_FA,-1,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
ENDIF
!
YNAME=TRIM(CMASK)//TRIM(HREC)
CALL  FAECR_L_D(KRESP,NUNIT_FA,YNAME,KL,OFIELD)
IF (KRESP/=0) THEN
  CALL ERROR_WRITE_SURF_FA(HREC,KRESP)
ENDIF
!
IF(LFANOCOMPACT)THEN
  ! On remet la valeur par defaut 
  CALL FAGOTE(KRESP,NUNIT_FA,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFL1_FA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFL1_FA
!
!     #############################################################
      SUBROUTINE WRITE_SURFX2_FA(HREC,KL1,KL2,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  * - routine to fill a write 2D array for the externalised surface 
!
USE MODI_ERROR_WRITE_SURF_FA
USE MODI_UNPACK_SAME_RANK
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_IO_SURF_FA, ONLY : NUNIT_FA, NMASK, NFULL, &
                            CMASK, LFANOCOMPACT  
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
INTEGER,                  INTENT(IN) :: KL1      ! number of points
INTEGER,                  INTENT(IN) :: KL2      ! 2nd dimension
REAL, DIMENSION(KL1,KL2), INTENT(IN) :: PFIELD   ! array containing the data field
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
REAL, DIMENSION(SIZE(PFIELD,2))       :: ZMEAN, ZCOUNT
CHARACTER(LEN=4)  :: YSUFFIX
CHARACTER(LEN=2)  :: YPATCH
INTEGER           :: I, JL ! loop counter
LOGICAL           :: LKNOWN
INTEGER           :: INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL
REAL(KIND=JPRB)   :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFX2_FA',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',LKNOWN)
IF (LKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFX2_FA',1,ZHOOK_HANDLE)
IF (LKNOWN) RETURN
!
CALL UNPACK_SAME_RANK(NMASK,PFIELD,ZWORK(:,:))
!
IF(LFANOCOMPACT)THEN
  CALL FAVEUR(KRESP,NUNIT_FA,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
  ! -- Pour ecrire sans compactage
  CALL FAGOTE(KRESP,NUNIT_FA,-1,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
  DO JL=1,SIZE(ZWORK,2)
    WRITE(YPATCH,'(I2.2)')JL
    YSUFFIX='S'//YPATCH//'_'
    CALL FAIENC(KRESP,NUNIT_FA,YSUFFIX,0,HREC,ZWORK(:,JL),.FALSE.)
    IF (KRESP/=0) THEN
       CALL ERROR_WRITE_SURF_FA(HREC,KRESP)
    ENDIF
  END DO
  ! On remet la valeur par defaut 
  CALL FAGOTE(KRESP,NUNIT_FA,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
ELSE
  ZMEAN (:)=0.0
  ZCOUNT(:)=0.0
  DO I=1,NFULL
     WHERE(ZWORK(I,:)/=XUNDEF)
        ZMEAN (:)=ZMEAN(:)+ZWORK(I,:)
        ZCOUNT(:)=ZCOUNT(:)+1.0
     ENDWHERE
  ENDDO
  WHERE(ZCOUNT(:)>0.0)ZMEAN(:)=ZMEAN(:)/ZCOUNT(:)        
  DO JL=1,SIZE(ZWORK,2)
    WHERE(ZWORK(:,JL)==XUNDEF)ZWORK(:,JL)=ZMEAN(JL)
    WRITE(YPATCH,'(I2.2)')JL
    YSUFFIX='S'//YPATCH//'_'
    CALL FAIENC(KRESP,NUNIT_FA,YSUFFIX,0,HREC,ZWORK(:,JL),.FALSE.)
    IF (KRESP/=0) THEN
       CALL ERROR_WRITE_SURF_FA(HREC,KRESP)
    ENDIF
  END DO
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFX2_FA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFX2_FA
!
!     #############################################################
      SUBROUTINE WRITE_SURFT0_FA(HREC,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
!     #############################################################
!
!!****  * - routine to write a date
!
USE MODI_ERROR_WRITE_SURF_FA
!
USE MODD_IO_SURF_FA, ONLY : CMASK, NUNIT_FA, LFANOCOMPACT
!
USE MODE_FASURFEX
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
INTEGER               :: IRET
INTEGER               :: IHOUR, IMIN
LOGICAL               :: LKNOWN
CHARACTER(LEN=16)     :: YNAME ! Field Name
INTEGER               :: INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFT0_FA',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',LKNOWN)
IF (LKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFT0_FA',1,ZHOOK_HANDLE)
IF (LKNOWN) RETURN
!
IF (HREC=='DTCUR') THEN
  IHOUR = NINT(PTIME)/3600
  IMIN  = NINT(PTIME)/60 - IHOUR * 60
  CALL FANDAR(IRET,NUNIT_FA,(/ KYEAR, KMONTH, KDAY, IHOUR, IMIN, 1, 0, 0, 0, 0, 0 /))
END IF
!
ITDATE(1) = KYEAR
ITDATE(2) = KMONTH
ITDATE(3) = KDAY
!
IF(LFANOCOMPACT)THEN
  CALL FAVEUR(KRESP,NUNIT_FA,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
  ! -- Pour ecrire sans compactage
  CALL FAGOTE(KRESP,NUNIT_FA,-1,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
ENDIF
!
YNAME=TRIM(CMASK)//TRIM(HREC)//'%TDATE'
CALL  FAECR_I_D(KRESP,NUNIT_FA,YNAME,3,ITDATE)
IF (KRESP/=0) THEN
  CALL ERROR_WRITE_SURF_FA(HREC,KRESP)
ENDIF
!
YNAME=TRIM(CMASK)//TRIM(HREC)//'%TIME'
CALL  FAECR_R(KRESP,NUNIT_FA,YNAME,PTIME)
IF (KRESP/=0) THEN
  CALL ERROR_WRITE_SURF_FA(HREC,KRESP)
ENDIF
!
IF(LFANOCOMPACT)THEN
  ! On remet la valeur par defaut 
  CALL FAGOTE(KRESP,NUNIT_FA,INGRIB,INBPDG,INBCSP,ISTRON,IPUILA,IDMOPL)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFT0_FA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFT0_FA
!
!     #############################################################
      SUBROUTINE WRITE_SURFT2_FA(HREC,KL1,KL2,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
!     #############################################################
!
!!****  * - routine to write a date
!
USE MODI_ERROR_WRITE_SURF_FA
!
USE MODD_IO_SURF_FA, ONLY : NUNIT_FA, CMASK, NLUOUT
!
USE MODE_FASURFEX
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_IO_BUFF_n
!
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=16),  INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER,                      INTENT(IN) :: KL1      ! number of points
INTEGER,                      INTENT(IN) :: KL2      ! 2nd dimension
INTEGER, DIMENSION(KL1,KL2), INTENT(IN)  :: KYEAR    ! year
INTEGER, DIMENSION(KL1,KL2), INTENT(IN)  :: KMONTH   ! month
INTEGER, DIMENSION(KL1,KL2), INTENT(IN)  :: KDAY     ! day
REAL,    DIMENSION(KL1,KL2), INTENT(IN)  :: PTIME    ! time
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(IN)  :: HCOMMENT ! comment string

!*      0.2   Declarations of local variables
!
INTEGER, DIMENSION(3,SIZE(KYEAR,1),SIZE(KYEAR,2)) :: ITDATE
LOGICAL          :: LKNOWN
CHARACTER(LEN=16):: YNAME ! Field Name
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFT2_FA',0,ZHOOK_HANDLE)
KRESP=0
!
CALL IO_BUFF_n(HREC,'W',LKNOWN)
IF (LKNOWN .AND. LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFT2_FA',1,ZHOOK_HANDLE)
IF (LKNOWN) RETURN
!
ITDATE(1,:,:) = KYEAR  (:,:)
ITDATE(2,:,:) = KMONTH (:,:)
ITDATE(3,:,:) = KDAY   (:,:)
!
YNAME=TRIM(CMASK)//TRIM(HREC)
WRITE(NLUOUT,*) ' WRITE_SURFT2_FA : time in 2 dimensions not yet implemented : YNAME=',YNAME,'ITDATE=',ITDATE
CALL ABOR1_SFX('MODE_WRITE_SURF_FA:WRITE_SURFT2_FA: time in 2 dimensions not yet implemented')
!
IF (LHOOK) CALL DR_HOOK('MODE_WRITE_SURF_FA:WRITE_SURFT2_FA',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE WRITE_SURFT2_FA
!
END MODULE MODE_WRITE_SURF_FA
