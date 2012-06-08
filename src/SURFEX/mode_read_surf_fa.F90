MODULE MODE_READ_SURF_FA
!
INTERFACE READ_SURF0_FA
        MODULE PROCEDURE READ_SURFX0_FA
        MODULE PROCEDURE READ_SURFN0_FA
        MODULE PROCEDURE READ_SURFL0_FA
        MODULE PROCEDURE READ_SURFC0_FA
END INTERFACE
INTERFACE READ_SURFN_FA
        MODULE PROCEDURE READ_SURFX1_FA
        MODULE PROCEDURE READ_SURFN1_FA
        MODULE PROCEDURE READ_SURFL1_FA
        MODULE PROCEDURE READ_SURFX2_FA
END INTERFACE
INTERFACE READ_SURFT_FA
        MODULE PROCEDURE READ_SURFT0_FA
        MODULE PROCEDURE READ_SURFT2_FA
END INTERFACE
!
CONTAINS
!
!     #############################################################
      SUBROUTINE READ_SURFX0_FA(HREC,PFIELD,KRESP,HCOMMENT)
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
USE MODI_ERROR_READ_SURF_FA
!
USE MODD_IO_SURF_FA,        ONLY : NUNIT_FA, NLUOUT, CMASK
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
CHARACTER(LEN=16):: YNAME ! Field Name
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFX0_FA',0,ZHOOK_HANDLE)
KRESP=0
!
YMASK=CMASK
CALL IO_BUFF_n(HREC,'R',GKNOWN)
IF (GKNOWN) YMASK='FULL  '
!
YNAME=TRIM(YMASK)//TRIM(HREC)
CALL FALIT_R(KRESP,NUNIT_FA,YNAME,PFIELD)
IF (KRESP/=0) THEN
  CALL ERROR_READ_SURF_FA(HREC,KRESP)
ENDIF
!
YCOMMENT = TRIM(YNAME)
HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFX0_FA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFX0_FA
!
!     #############################################################
      SUBROUTINE READ_SURFX1_FA(HREC,KL,PFIELD,KRESP,HCOMMENT,HDIR)
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
USE MODI_ERROR_READ_SURF_FA
USE MODI_PACK_SAME_RANK
!
USE MODD_IO_SURF_FA,        ONLY : NUNIT_FA, NLUOUT, NMASK, NFULL, NFULL_EXT, &
                                     NDGL, NDLON, NDGUX, NDLUX  
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
CHARACTER(LEN=16),   INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER,             INTENT(IN)  :: KL       ! number of points
REAL, DIMENSION(KL), INTENT(OUT) :: PFIELD   ! array containing the data field
INTEGER,             INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),  INTENT(OUT) :: HCOMMENT ! comment
CHARACTER(LEN=1),    INTENT(IN)  :: HDIR     ! type of field :
                                             ! 'H' : field with
                                             !       horizontal spatial dim.
                                             ! '-' : no horizontal dim.

!
!*      0.2   Declarations of local variables
!
LOGICAL                    :: GFOUND
CHARACTER(LEN=50)          :: YCOMMENT
CHARACTER(LEN=3)           :: YPREF
CHARACTER(LEN=13)          :: YSUFF
REAL, DIMENSION(NFULL)     :: ZWORK   ! work array read in the file
REAL, DIMENSION(NFULL_EXT) :: ZWORK2
INTEGER ::  I,J
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------------
!
 IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFX1_FA',0,ZHOOK_HANDLE)
 KRESP=0
!
 YPREF=HREC(1:3)
 YSUFF=HREC(4:16) 
!
 IF (YPREF=='CLS' .OR. YPREF=='SUR' .OR. YPREF=='PRO' .OR. YPREF=='ATM') THEN
   CALL FACILE(KRESP,NUNIT_FA,HREC(1:4),0,HREC(5:16),ZWORK2,.FALSE.)
   IF (KRESP/=0) THEN
     CALL ERROR_READ_SURF_FA(HREC,KRESP)
   ENDIF
   DO J=1,NDGUX
     DO I=1,NDLUX
       ZWORK((J-1)*NDLUX + I) = ZWORK2((J-1)*NDLON + I)
     ENDDO
   ENDDO
   YCOMMENT = TRIM(HREC)
 ELSE
   CALL FACILE(KRESP,NUNIT_FA,'S1D_',0,HREC,ZWORK,.FALSE.)
   IF (KRESP/=0) THEN
     CALL ERROR_READ_SURF_FA(HREC,KRESP)
   ENDIF   
   YCOMMENT = 'S1D_'//TRIM(HREC)
 ENDIF
!
 HCOMMENT = YCOMMENT
!
IF(HDIR=='H')THEN
   CALL PACK_SAME_RANK(NMASK,ZWORK(:),PFIELD)
ELSE !no horizontal dim. case (not masked)
   PFIELD(:)=ZWORK(1:KL)
ENDIF
! 
 IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFX1_FA',1,ZHOOK_HANDLE)
! 
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFX1_FA
!
!     #############################################################
      SUBROUTINE READ_SURFX2_FA(HREC,KL1,KL2,PFIELD,KRESP,HCOMMENT,HDIR)
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
USE MODI_ERROR_READ_SURF_FA
USE MODI_PACK_SAME_RANK
!
USE MODD_IO_SURF_FA,        ONLY : NUNIT_FA, NLUOUT, NMASK, NFULL
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
CHARACTER(LEN=16),        INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER,                  INTENT(IN)  :: KL1      ! number of points
INTEGER,                  INTENT(IN)  :: KL2      ! 2nd dimension
REAL, DIMENSION(KL1,KL2), INTENT(OUT) :: PFIELD   ! array containing the data field
INTEGER,                  INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100),       INTENT(OUT) :: HCOMMENT ! comment
CHARACTER(LEN=1),         INTENT(IN)  :: HDIR     ! type of field :
                                                  ! 'H' : field with
                                                  !       horizontal spatial dim.
                                                  ! '-' : no horizontal dim.
!
!*      0.2   Declarations of local variables
! 
LOGICAL                               :: GFOUND
CHARACTER(LEN=50)                     :: YCOMMENT
REAL, DIMENSION(NFULL,SIZE(PFIELD,2)) :: ZWORK   ! work array read in the file
CHARACTER(LEN=4)                      :: YSUFFIX
CHARACTER(LEN=2)                      :: YPATCH

INTEGER :: JL ! loop counter
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFX2_FA',0,ZHOOK_HANDLE)
KRESP=0
!
DO JL=1,SIZE(ZWORK,2)
  WRITE(YPATCH,'(I2.2)')JL
  YSUFFIX='S'//YPATCH//'_'
  CALL FACILE(KRESP,NUNIT_FA,YSUFFIX,JL,HREC,ZWORK(:,JL),.FALSE.)
   IF (KRESP/=0) THEN
     YCOMMENT = YSUFFIX//TRIM(HREC)
     CALL ERROR_READ_SURF_FA(YCOMMENT,KRESP)
   ENDIF  
END DO
!
YCOMMENT = 'PATCH_'//TRIM(HREC)
HCOMMENT = YCOMMENT
!
CALL PACK_SAME_RANK(NMASK,ZWORK(:,:),PFIELD)
!  
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFX2_FA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFX2_FA
!
!     #############################################################
      SUBROUTINE READ_SURFN0_FA(HREC,KFIELD,KRESP,HCOMMENT)
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
USE MODI_ERROR_READ_SURF_FA
!
USE MODD_IO_SURF_FA,        ONLY : NUNIT_FA, NLUOUT, NMASK, CMASK
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
CHARACTER(LEN=16):: YNAME ! Field Name
REAL(KIND=JPRB)  :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFN0_FA',0,ZHOOK_HANDLE)
KRESP=0
!
!
YMASK=CMASK
CALL IO_BUFF_n(HREC,'R',GKNOWN)
IF (GKNOWN) YMASK='FULL  '
!
!
YNAME=TRIM(YMASK)//TRIM(HREC)
CALL FALIT_I(KRESP,NUNIT_FA,YNAME,KFIELD)
IF (KRESP/=0) THEN
  CALL ERROR_READ_SURF_FA(HREC,KRESP)
ENDIF
!
YCOMMENT = YNAME
HCOMMENT = YCOMMENT
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFN0_FA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFN0_FA
!
!     #############################################################
      SUBROUTINE READ_SURFN1_FA(HREC,KL,KFIELD,KRESP,HCOMMENT,HDIR)
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
USE MODI_ERROR_READ_SURF_FA
USE MODI_PACK_SAME_RANK
!
USE MODD_IO_SURF_FA,        ONLY : NUNIT_FA, NLUOUT, NMASK, NFULL, CMASK
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
CHARACTER(LEN=16),      INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER,                INTENT(IN)  :: KL       ! number of points
INTEGER, DIMENSION(KL), INTENT(OUT) :: KFIELD   ! the integer to be read
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
INTEGER, DIMENSION(NFULL) :: IWORK  ! work array read in the file
CHARACTER(LEN=6) :: YMASK
CHARACTER(LEN=16):: YNAME ! Field Name
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!---------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFN1_FA',0,ZHOOK_HANDLE)
KRESP = 0
!
YMASK=CMASK
CALL IO_BUFF_n(HREC,'R',GKNOWN)
IF (GKNOWN) YMASK='FULL  '
!
YNAME=TRIM(YMASK)//TRIM(HREC)
IF (HDIR=="-") THEN
  CALL FALIT_I_D(KRESP,NUNIT_FA,YNAME,KL,KFIELD)
  IF (KRESP/=0) THEN
    CALL ERROR_READ_SURF_FA(HREC,KRESP)
  ENDIF

ELSE
  CALL FALIT_I_D(KRESP,NUNIT_FA,YNAME,NFULL,IWORK)
  IF (KRESP/=0) THEN
    CALL ERROR_READ_SURF_FA(HREC,KRESP)
  ENDIF
  CALL PACK_SAME_RANK(NMASK,IWORK(:),KFIELD)
ENDIF
!
YCOMMENT = YNAME
HCOMMENT = YCOMMENT
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFN1_FA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFN1_FA
!
!     #############################################################
      SUBROUTINE READ_SURFC0_FA(HREC,HFIELD,KRESP,HCOMMENT)
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
USE MODI_ERROR_READ_SURF_FA
!
USE MODD_IO_SURF_FA,        ONLY : NUNIT_FA, NLUOUT, CMASK
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
CHARACTER(LEN=40),  INTENT(OUT) :: HFIELD    ! the integer to be read
INTEGER,            INTENT(OUT) :: KRESP     ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT  ! comment
!
!*      0.2   Declarations of local variables
!
LOGICAL           :: GFOUND
CHARACTER,DIMENSION(40) :: YFIELD
CHARACTER(LEN=50):: YCOMMENT
CHARACTER(LEN=6) :: YMASK
CHARACTER(LEN=16):: YNAME ! Field Name
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFC0_FA',0,ZHOOK_HANDLE)
KRESP=0
!
YMASK=CMASK
CALL IO_BUFF_n(HREC,'R',GKNOWN)
IF (GKNOWN) YMASK='FULL  '
!
YNAME=TRIM(YMASK)//TRIM(HREC)
CALL FALIT_C(KRESP,NUNIT_FA,YNAME,40,YFIELD)
IF (KRESP/=0) THEN
  CALL ERROR_READ_SURF_FA(HREC,KRESP)
ENDIF
WRITE(HFIELD,'(40A1)') YFIELD(:)
!
YCOMMENT = YNAME
HCOMMENT = YCOMMENT
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFC0_FA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFC0_FA
!
!     #############################################################
      SUBROUTINE READ_SURFL1_FA(HREC,KL,OFIELD,KRESP,HCOMMENT,HDIR)
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
USE MODI_ERROR_READ_SURF_FA
!
USE MODD_IO_SURF_FA,        ONLY : NUNIT_FA, NLUOUT, CMASK
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
CHARACTER(LEN=16),      INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER,                INTENT(IN)  :: KL       ! number of points
LOGICAL, DIMENSION(KL), INTENT(OUT) :: OFIELD   ! array containing the data field
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
CHARACTER(LEN=16):: YNAME ! Field Name
CHARACTER(LEN=6) :: YMASK
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFL1_FA',0,ZHOOK_HANDLE)
KRESP=0
!
YMASK=CMASK
CALL IO_BUFF_n(HREC,'R',GKNOWN)
IF (GKNOWN) YMASK='FULL  '
!
YNAME=TRIM(YMASK)//TRIM(HREC)
CALL FALIT_L_D(KRESP,NUNIT_FA,YNAME,KL,OFIELD)
IF (KRESP/=0) THEN
  CALL ERROR_READ_SURF_FA(HREC,KRESP)
ENDIF
!
YCOMMENT = YNAME
HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFL1_FA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFL1_FA
!
!
!     #############################################################
      SUBROUTINE READ_SURFL0_FA(HREC,OFIELD,KRESP,HCOMMENT)
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
USE MODI_ERROR_READ_SURF_FA
!
USE MODD_IO_SURF_FA,        ONLY : NUNIT_FA, NLUOUT, CMASK
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
LOGICAL,            INTENT(OUT) :: OFIELD   ! array containing the data field
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT ! comment
!
!*      0.2   Declarations of local variables
!
LOGICAL           :: GFOUND
CHARACTER(LEN=50):: YCOMMENT
CHARACTER(LEN=6) :: YMASK
CHARACTER(LEN=16):: YNAME ! Field Name
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFL0_FA',0,ZHOOK_HANDLE)
KRESP=0
!
YMASK=CMASK
CALL IO_BUFF_n(HREC,'R',GKNOWN)
IF (GKNOWN) YMASK='FULL  '
!
YNAME=TRIM(YMASK)//TRIM(HREC)
CALL FALIT_L(KRESP,NUNIT_FA,YNAME,OFIELD)
IF (KRESP/=0) THEN
  CALL ERROR_READ_SURF_FA(HREC,KRESP)
ENDIF
!
YCOMMENT = YNAME
HCOMMENT = YCOMMENT
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFL0_FA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFL0_FA
!
!     #############################################################
      SUBROUTINE READ_SURFT0_FA(HREC,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
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
USE MODI_ERROR_READ_SURF_FA
!
USE MODD_IO_SURF_FA,        ONLY : NUNIT_FA, NLUOUT, CMASK
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
!
CHARACTER(LEN=6) :: YMASK
CHARACTER(LEN=16):: YNAME ! Field Name
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFT0_FA',0,ZHOOK_HANDLE)
KRESP=0
!
YMASK=CMASK
CALL IO_BUFF_n(HREC,'R',GKNOWN)
IF (GKNOWN) YMASK='FULL  '
!
!
!-------------------------------------------------------------------------------
!
YNAME=TRIM(YMASK)//TRIM(HREC)//'%TDATE'
CALL FALIT_I_D(KRESP,NUNIT_FA,YNAME,3,ITDATE)
IF (KRESP/=0) THEN
  CALL ERROR_READ_SURF_FA(HREC,KRESP)
ENDIF
!
KYEAR  = ITDATE(1)
KMONTH = ITDATE(2)
KDAY   = ITDATE(3)
!
!-------------------------------------------------------------------------------
!
YNAME=TRIM(YMASK)//TRIM(HREC)//'%TIME'
CALL FALIT_R(KRESP,NUNIT_FA,YNAME,PTIME)
IF (KRESP/=0) THEN
  CALL ERROR_READ_SURF_FA(HREC,KRESP)
ENDIF
!
!-------------------------------------------------------------------------------
!
YCOMMENT = TRIM(HREC)
HCOMMENT = YCOMMENT
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFT0_FA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_SURFT0_FA
!
!     #############################################################
      SUBROUTINE READ_SURFT2_FA(HREC,KL1,KL2,KYEAR,KMONTH,KDAY,PTIME,KRESP,HCOMMENT)
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
USE MODI_ERROR_READ_SURF_FA
!
USE MODD_IO_SURF_FA,        ONLY : NUNIT_FA, NLUOUT, CMASK
!
USE MODE_FASURFEX
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
USE MODI_IO_BUFF_n
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=16),  INTENT(IN)  :: HREC     ! name of the article to be read
INTEGER                                  :: KL1, KL2
INTEGER, DIMENSION(KL1,KL2), INTENT(OUT) :: KYEAR    ! year
INTEGER, DIMENSION(KL1,KL2), INTENT(OUT) :: KMONTH   ! month
INTEGER, DIMENSION(KL1,KL2), INTENT(OUT) :: KDAY     ! day
REAL,    DIMENSION(KL1,KL2), INTENT(OUT) :: PTIME    ! year
INTEGER,            INTENT(OUT) :: KRESP    ! KRESP  : return-code if a problem appears
CHARACTER(LEN=100), INTENT(OUT) :: HCOMMENT ! comment

!*      0.2   Declarations of local variables
!
LOGICAL               :: GFOUND
CHARACTER(LEN=50):: YCOMMENT
INTEGER, DIMENSION(3,SIZE(KYEAR,1),SIZE(KYEAR,2)) :: ITDATE
CHARACTER(LEN=6) :: YMASK
CHARACTER(LEN=16):: YNAME ! Field Name
LOGICAL          :: GKNOWN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFT2_FA',0,ZHOOK_HANDLE)
KRESP=0
!
YMASK=CMASK
CALL IO_BUFF_n(HREC,'R',GKNOWN)
IF (GKNOWN) YMASK='FULL  '
!
!
!-------------------------------------------------------------------------------
!
YNAME=TRIM(CMASK)//TRIM(HREC)
WRITE(NLUOUT,*) ' READ_SURFT2_FA : time in 2 dimensions not yet implemented : YNAME=',YNAME
CALL ABOR1_SFX('MODE_READ_SURF_FA:READ_SURFT2_FA: time in 2 dimensions not yet implemented')
!
!-------------------------------------------------------------------------------
!
HCOMMENT = YCOMMENT
!
IF (LHOOK) CALL DR_HOOK('MODE_READ_SURF_FA:READ_SURFT2_FA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_SURFT2_FA
!
END MODULE MODE_READ_SURF_FA
