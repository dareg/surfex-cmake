!##################
MODULE MODI_READ_SURF
!##################
!
  INTERFACE READ_SURF
!
     SUBROUTINE READ_SURFX0(HPROGRAM,HREC,PFIELD,KRESP,HCOMMENT)
CHARACTER(LEN=6)          ,INTENT(IN) :: HPROGRAM ! calling program

CHARACTER(LEN=*)          ,INTENT(IN) :: HREC   ! name of the article to be read

REAL, &
                           INTENT(INOUT):: PFIELD ! real scalar to be read  

INTEGER                  ,INTENT(OUT) :: KRESP          ! KRESP  : return-code if a problem appears 
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read
!
END SUBROUTINE READ_SURFX0
!
     SUBROUTINE READ_SURFX1(HPROGRAM,HREC,PFIELD,KRESP,HCOMMENT,HDIR)
CHARACTER(LEN=6)         ,INTENT(IN)  ::HPROGRAM ! calling program

CHARACTER(LEN=*)          ,INTENT(IN) ::HREC   ! name of the article to be read

REAL, DIMENSION(:), &
                             INTENT(OUT)::PFIELD ! array containing the data field  

INTEGER                  ,INTENT(OUT) :: KRESP          ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read
CHARACTER(LEN=1),OPTIONAL,INTENT(IN)  :: HDIR       ! type of field :
!                                                   ! 'H' : field with
!                                                   !       horizontal spatial dim.
!                                                   ! '-' : no horizontal dim.
!
END SUBROUTINE READ_SURFX1
!
!
     SUBROUTINE READ_SURFX2(HPROGRAM,HREC,PFIELD,KRESP,HCOMMENT,HDIR)
CHARACTER(LEN=6)         ,INTENT(IN)  ::HPROGRAM ! calling program

CHARACTER(LEN=*)          ,INTENT(IN) ::HREC   ! name of the article to be read

REAL, DIMENSION(:,:), &
                             INTENT(OUT)::PFIELD ! array containing the data field  

INTEGER                  ,INTENT(OUT) :: KRESP          ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read
CHARACTER(LEN=1),OPTIONAL,INTENT(IN)  :: HDIR       ! type of field :
!                                                   ! 'H' : field with
!                                                   !       horizontal spatial dim.
!                                                   ! '-' : no horizontal dim.

!
END SUBROUTINE READ_SURFX2
!
      SUBROUTINE READ_SURFX2COV(HPROGRAM,HREC,PFIELD,OFLAG,KRESP,HCOMMENT,HDIR)
!
!
CHARACTER(LEN=6)       ,  INTENT(IN)  :: HPROGRAM ! calling program

CHARACTER(LEN=*)          ,INTENT(IN) ::HREC   ! name of the article to be read

REAL, DIMENSION(:,:), &
                           INTENT(OUT)::PFIELD ! array containing the data field
LOGICAL,DIMENSION(:),      INTENT(IN) ::OFLAG  ! mask for array filling

INTEGER                  ,INTENT(OUT) :: KRESP          ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read
CHARACTER(LEN=1),OPTIONAL,INTENT(IN)  :: HDIR       ! type of field :
!                                                   ! 'H' : field with
!                                                   !       horizontal spatial dim.
!                                                   ! '-' : no horizontal dim.
END SUBROUTINE READ_SURFX2COV
!
     SUBROUTINE READ_SURFX3(HPROGRAM,HREC,PFIELD,KRESP,HCOMMENT,HDIR)
CHARACTER(LEN=6)         ,INTENT(IN)  ::HPROGRAM ! calling program

CHARACTER(LEN=*)          ,INTENT(IN) ::HREC   ! name of the article to be read

REAL, DIMENSION(:,:,:), &
                             INTENT(OUT)::PFIELD ! array containing the data field  

INTEGER                  ,INTENT(OUT) :: KRESP          ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read
CHARACTER(LEN=1),OPTIONAL,INTENT(IN)  :: HDIR       ! type of field :
!                                                   ! 'H' : field with
!                                                   !       horizontal spatial dim.
!                                                   ! '-' : no horizontal dim.

!
END SUBROUTINE READ_SURFX3
!
     SUBROUTINE READ_SURFN0(HPROGRAM,HREC,KFIELD,KRESP,HCOMMENT)
CHARACTER(LEN=6)         ,INTENT(IN)  ::HPROGRAM ! calling program

CHARACTER(LEN=*)          ,INTENT(IN) ::HREC   ! name of the article to be read

INTEGER, &
                             INTENT(OUT)::KFIELD ! integer to be read  

INTEGER                  ,INTENT(OUT) :: KRESP          ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read
!
!
END SUBROUTINE READ_SURFN0
!
     SUBROUTINE READ_SURFN1(HPROGRAM,HREC,KFIELD,KRESP,HCOMMENT,HDIR)
CHARACTER(LEN=6)         ,INTENT(IN)  ::HPROGRAM ! calling program

CHARACTER(LEN=*)          ,INTENT(IN) ::HREC   ! name of the article to be read

INTEGER, DIMENSION(:), &
                             INTENT(OUT)::KFIELD ! integer to be read  

INTEGER                  ,INTENT(OUT) :: KRESP          ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read
CHARACTER(LEN=1),OPTIONAL,INTENT(IN)  :: HDIR       ! type of field :
!                                                   ! 'H' : field with
!                                                   !       horizontal spatial dim.
!                                                   ! '-' : no horizontal dim.
!
!
END SUBROUTINE READ_SURFN1
!

!
     SUBROUTINE READ_SURFC0(HPROGRAM,HREC,HFIELD,KRESP,HCOMMENT)
CHARACTER(LEN=6)         ,INTENT(IN)  :: HPROGRAM   ! calling program

CHARACTER(LEN=*)          ,INTENT(IN) :: HREC       ! name of the article to be read

CHARACTER(LEN=*), &
                             INTENT(OUT):: HFIELD     ! caracter to be read  

INTEGER                  ,INTENT(OUT) :: KRESP      ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read
!
END SUBROUTINE READ_SURFC0
!
      SUBROUTINE READ_SURFL1(HPROGRAM,HREC,OFIELD,KRESP,HCOMMENT,HDIR)
CHARACTER(LEN=6)       ,    INTENT(IN)::HPROGRAM ! calling program

CHARACTER(LEN=*)          ,INTENT(IN) ::HREC   ! name of the article to be read

LOGICAL, DIMENSION(:), &
                             INTENT(OUT)::OFIELD ! array containing the data field  

INTEGER                  ,INTENT(OUT) :: KRESP          ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read
CHARACTER(LEN=1),OPTIONAL,INTENT(IN)  :: HDIR       ! type of field :
!                                                   ! 'H' : field with
!                                                   !       horizontal spatial dim.
!                                                   ! '-' : no horizontal dim.

!
END SUBROUTINE READ_SURFL1
!
      SUBROUTINE READ_SURFL0(HPROGRAM,HREC,OFIELD,KRESP,HCOMMENT)
CHARACTER(LEN=6)       ,    INTENT(IN)::HPROGRAM ! calling program

CHARACTER(LEN=*)          ,INTENT(IN) ::HREC   ! name of the article to be read

LOGICAL,                  INTENT(OUT)::OFIELD ! array containing the data field

INTEGER                  ,INTENT(OUT) :: KRESP          ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read
!
END SUBROUTINE READ_SURFL0
!
      SUBROUTINE READ_SURFT0(HPROGRAM,HREC,TFIELD,KRESP,HCOMMENT)
!
USE MODD_TYPE_DATE_SURF
!
CHARACTER(LEN=6)       ,  INTENT(IN)  :: HPROGRAM ! calling program

CHARACTER(LEN=*)          ,INTENT(IN) ::HREC   ! name of the article to be read
TYPE (DATE_TIME),        INTENT(INOUT)::TFIELD ! array containing the data field  

INTEGER                  ,INTENT(OUT) :: KRESP          ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read
!
END SUBROUTINE READ_SURFT0
!
      SUBROUTINE READ_SURFT1(HPROGRAM,HREC,TFIELD,KRESP,HCOMMENT)
!
USE MODD_TYPE_DATE_SURF
!
CHARACTER(LEN=6)       ,  INTENT(IN)  :: HPROGRAM ! calling program

CHARACTER(LEN=*)          ,INTENT(IN) ::HREC   ! name of the article to be read
TYPE (DATE_TIME), DIMENSION(:), INTENT(INOUT)::TFIELD ! array containing the data field  

INTEGER                  ,INTENT(OUT) :: KRESP          ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read
!
END SUBROUTINE READ_SURFT1
!
      SUBROUTINE READ_SURFT2(HPROGRAM,HREC,TFIELD,KRESP,HCOMMENT)
!
USE MODD_TYPE_DATE_SURF
!
CHARACTER(LEN=6)       ,  INTENT(IN)  :: HPROGRAM ! calling program

CHARACTER(LEN=*)          ,INTENT(IN) ::HREC   ! name of the article to be read
TYPE (DATE_TIME), DIMENSION(:,:), INTENT(INOUT)::TFIELD ! array containing the data field  

INTEGER                  ,INTENT(OUT) :: KRESP          ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read
!
END SUBROUTINE READ_SURFT2
!
END INTERFACE
!
END MODULE MODI_READ_SURF
!
!     #############################################################
      SUBROUTINE READ_SURFX0(HPROGRAM,HREC,PFIELD,KRESP,HCOMMENT)
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
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
#ifdef OL
USE MODE_READ_SURF_OL, ONLY: READ_SURF0_OL
#endif
#ifdef ASC
USE MODE_READ_SURF_ASC, ONLY: READ_SURF0_ASC
#endif
#ifdef FA
USE MODE_READ_SURF_FA, ONLY: READ_SURF0_FA
#endif
#ifdef LFI
USE MODE_READ_SURF_LFI, ONLY: READ_SURF0_LFI
#endif
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=6)       ,  INTENT(IN)  :: HPROGRAM ! calling program

CHARACTER(LEN=*)          ,INTENT(IN) ::HREC   ! name of the article to be read

REAL, &
                             INTENT(OUT)::PFIELD ! the real scalar to be read  

INTEGER                  ,INTENT(OUT) :: KRESP          ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER(LEN=16)  :: YREC
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFX0',0,ZHOOK_HANDLE)
YREC = HREC
!
IF (HPROGRAM=='MESONH') THEN
#ifdef MNH
  CALL READ_SURFX0_MNH(YREC,PFIELD,KRESP,YCOMMENT)
#endif
ENDIF
!
IF (HPROGRAM=='OFFLIN') THEN
#ifdef OL   
  CALL READ_SURF0_OL(YREC,PFIELD,KRESP,YCOMMENT)
#endif
ENDIF
!
IF (HPROGRAM=='ASCII ') THEN
#ifdef ASC 
  CALL READ_SURF0_ASC(YREC,PFIELD,KRESP,YCOMMENT)
#endif
ENDIF
!
IF (HPROGRAM=='AROME ') THEN
#ifdef ARO 
  CALL READ_SURFX0_ARO(YREC,PFIELD,KRESP,YCOMMENT)
#endif
ENDIF
!
IF (HPROGRAM=='FA    ') THEN
#ifdef FA 
  CALL READ_SURF0_FA(YREC,PFIELD,KRESP,YCOMMENT)
#endif
ENDIF
!
IF (HPROGRAM=='LFI   ') THEN
#ifdef LFI
  CALL READ_SURF0_LFI(YREC,PFIELD,KRESP,YCOMMENT)
#endif
ENDIF
!-------------------------------------------------------------------------------
IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFX0',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFX0
!
!     #############################################################
      SUBROUTINE READ_SURFX1(HPROGRAM,HREC,PFIELD,KRESP,HCOMMENT,HDIR)
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
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
#ifdef OL
USE MODE_READ_SURF_OL, ONLY: READ_SURFN_OL
#endif
#ifdef ASC
USE MODE_READ_SURF_ASC, ONLY: READ_SURFN_ASC
#endif
#ifdef FA
USE MODE_READ_SURF_FA, ONLY: READ_SURFN_FA
#endif
#ifdef LFI
USE MODE_READ_SURF_LFI, ONLY: READ_SURFN_LFI
#endif
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=6)       ,  INTENT(IN)  :: HPROGRAM ! calling program

CHARACTER(LEN=*)          ,INTENT(IN) ::HREC   ! name of the article to be read

REAL, DIMENSION(:), &
                             INTENT(OUT)::PFIELD ! array containing the data field  

INTEGER                  ,INTENT(OUT) :: KRESP          ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read
CHARACTER(LEN=1),OPTIONAL,INTENT(IN)  :: HDIR       ! type of field :
!                                                   ! 'H' : field with
!                                                   !       horizontal spatial dim.
!                                                   ! '-' : no horizontal dim.
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER(LEN=16)  :: YREC
CHARACTER(LEN=1)   :: YDIR
INTEGER            :: IL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFX1',0,ZHOOK_HANDLE)
YREC = HREC
YDIR = 'H'
IF (PRESENT(HDIR)) YDIR = HDIR
!
IL = SIZE(PFIELD)
!-------------------------------------------------------------------------------
!
!
IF (HPROGRAM=='MESONH') THEN
#ifdef MNH 
  CALL READ_SURFX1_MNH(YREC,IL,PFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='OFFLIN') THEN
#ifdef OL 
  CALL READ_SURFN_OL(YREC,PFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='ASCII ') THEN
#ifdef ASC 
  CALL READ_SURFN_ASC(YREC,PFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='AROME ') THEN
#ifdef ARO
  CALL READ_SURFX1_ARO(YREC,IL,PFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='FA    ') THEN
#ifdef FA 
  CALL READ_SURFN_FA(YREC,IL,PFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='LFI   ') THEN
#ifdef LFI 
  CALL READ_SURFN_LFI(YREC,PFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!-------------------------------------------------------------------------------
IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFX1',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFX1
!
!     #############################################################
      SUBROUTINE READ_SURFX2(HPROGRAM,HREC,PFIELD,KRESP,HCOMMENT,HDIR)
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
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
#ifdef OL
USE MODE_READ_SURF_OL, ONLY: READ_SURFN_OL
#endif
#ifdef ASC
USE MODE_READ_SURF_ASC, ONLY: READ_SURFN_ASC
#endif
#ifdef FA
USE MODE_READ_SURF_FA, ONLY: READ_SURFN_FA
#endif
#ifdef LFI
USE MODE_READ_SURF_LFI, ONLY: READ_SURFN_LFI
#endif
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=6)       ,  INTENT(IN)  :: HPROGRAM ! calling program

CHARACTER(LEN=*)          ,INTENT(IN) ::HREC   ! name of the article to be read

REAL, DIMENSION(:,:), &
                             INTENT(OUT)::PFIELD ! array containing the data field  

INTEGER                  ,INTENT(OUT) :: KRESP          ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read
CHARACTER(LEN=1),OPTIONAL,INTENT(IN)  :: HDIR       ! type of field :
!                                                   ! 'H' : field with
!                                                   !       horizontal spatial dim.
!                                                   ! '-' : no horizontal dim.
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER(LEN=16)  :: YREC
CHARACTER(LEN=1)   :: YDIR
INTEGER            :: IL1, IL2
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFX2',0,ZHOOK_HANDLE)
YREC = HREC
YDIR = 'H'
IF (PRESENT(HDIR)) YDIR = HDIR
!
IL1 = SIZE(PFIELD,1)
IL2 = SIZE(PFIELD,2)
!-------------------------------------------------------------------------------
!
!
IF (HPROGRAM=='MESONH') THEN
#ifdef MNH 
  CALL READ_SURFX2_MNH(YREC,IL1,IL2,PFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='OFFLIN') THEN
#ifdef OL 
  CALL READ_SURFN_OL(YREC,PFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='ASCII ') THEN
#ifdef ASC
  CALL READ_SURFN_ASC(YREC,PFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='AROME ') THEN
#ifdef ARO 
  CALL READ_SURFX2_ARO(YREC,IL1,IL2,PFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='FA    ') THEN
#ifdef FA 
  CALL READ_SURFN_FA(YREC,IL1,IL2,PFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='LFI   ') THEN
#ifdef LFI 
  CALL READ_SURFN_LFI(YREC,PFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
!
!-------------------------------------------------------------------------------
IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFX2',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFX2
!
!     #############################################################
      SUBROUTINE READ_SURFX2COV(HPROGRAM,HREC,PFIELD,OFLAG,KRESP,HCOMMENT,HDIR)
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
#ifdef OL
USE MODE_READ_SURF_OL, ONLY: READ_SURFN_OL
#endif
#ifdef ASC
USE MODE_READ_SURF_ASC, ONLY: READ_SURFN_ASC
#endif
#ifdef FA
USE MODE_READ_SURF_FA, ONLY: READ_SURFN_FA
#endif
#ifdef LFI
USE MODE_READ_SURF_LFI, ONLY: READ_SURFN_LFI
#endif
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!*      0.    DECLARATIONS
!             ------------
!
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=6)       ,  INTENT(IN)  :: HPROGRAM ! calling program

CHARACTER(LEN=*)          ,INTENT(IN) ::HREC   ! name of the article to be read

REAL, DIMENSION(:,:), &
                           INTENT(OUT)::PFIELD ! array containing the data field
LOGICAL,DIMENSION(:),      INTENT(IN) ::OFLAG  ! mask for array filling

INTEGER                  ,INTENT(OUT) :: KRESP          ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read
CHARACTER(LEN=1),OPTIONAL,INTENT(IN)  :: HDIR       ! type of field :
!                                                   ! 'H' : field with
!                                                   !       horizontal spatial dim.
!                                                   ! '-' : no horizontal dim.
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER(LEN=16)  :: YREC
CHARACTER(LEN=1)   :: YDIR
INTEGER            :: JCOVER
INTEGER            :: IL1, IL2
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
YREC = HREC
YDIR = 'H'
IF (PRESENT(HDIR)) YDIR = HDIR
!
IL1 = SIZE(PFIELD,1)
IL2 = SIZE(PFIELD,2)
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFX2COV',0,ZHOOK_HANDLE)
!
IF (HPROGRAM=='MESONH') THEN
#ifdef MNH        
    CALL READ_SURFX2COV_MNH(YREC,IL1,IL2,PFIELD,OFLAG,KRESP,YCOMMENT,YDIR)
#endif
ELSE
  !
  DO JCOVER=1,IL2
    !
    WRITE(YREC,'(A5,I3.3)') 'COVER',JCOVER
    YCOMMENT='X_Y_'//YREC
    PFIELD(:,JCOVER)=0.
    IF (.NOT. OFLAG(JCOVER)) CYCLE
    !
    IF (HPROGRAM=='OFFLIN') THEN
#ifdef OL        
      CALL READ_SURFN_OL(YREC,PFIELD(:,JCOVER),KRESP,YCOMMENT,YDIR)
#endif
    ENDIF
    !
    IF (HPROGRAM=='ASCII ') THEN
#ifdef ASC        
      CALL READ_SURFN_ASC(YREC,PFIELD(:,JCOVER),KRESP,YCOMMENT,YDIR)
#endif
    ENDIF
    !
    IF (HPROGRAM=='AROME ') THEN
#ifdef ARO
      CALL READ_SURFX1_ARO(YREC,IL1,PFIELD(:,JCOVER),KRESP,YCOMMENT,YDIR)
#endif
    ENDIF
    !
    IF (HPROGRAM=='FA    ') THEN
#ifdef FA 
      CALL READ_SURFN_FA(YREC,IL1,PFIELD(:,JCOVER),KRESP,YCOMMENT,YDIR)
#endif
    ENDIF
    !
    IF (HPROGRAM=='LFI   ') THEN
#ifdef LFI 
      CALL READ_SURFN_LFI(YREC,PFIELD(:,JCOVER),KRESP,YCOMMENT,YDIR)
#endif
    ENDIF
  END DO
  !
ENDIF
!
!-------------------------------------------------------------------------------
IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFX2COV',1,ZHOOK_HANDLE)
!
END SUBROUTINE READ_SURFX2COV
!
!     #############################################################

!     #############################################################
      SUBROUTINE READ_SURFX3(HPROGRAM,HREC,PFIELD,KRESP,HCOMMENT,HDIR)
!     #############################################################
!
!!****  *READX3* - routine to fill a real 3D array for the externalised surface 
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READ_SURFX3 is
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
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
#ifdef OL
USE MODE_READ_SURF_OL, ONLY: READ_SURFN_OL
#endif
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=6)       ,  INTENT(IN)  :: HPROGRAM ! calling program

CHARACTER(LEN=*)          ,INTENT(IN) ::HREC   ! name of the article to be read

REAL, DIMENSION(:,:,:), &
                             INTENT(OUT)::PFIELD ! array containing the data field  

INTEGER                  ,INTENT(OUT) :: KRESP          ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read
CHARACTER(LEN=1),OPTIONAL,INTENT(IN)  :: HDIR       ! type of field :
!                                                   ! 'H' : field with
!                                                   !       horizontal spatial dim.
!                                                   ! '-' : no horizontal dim.
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER(LEN=16)  :: YREC
CHARACTER(LEN=1)   :: YDIR
INTEGER            :: IL1, IL2, IL3
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFX3',0,ZHOOK_HANDLE)
YREC = HREC
YDIR = 'H'
IF (PRESENT(HDIR)) YDIR = HDIR
!
IL1 = SIZE(PFIELD,1)
IL2 = SIZE(PFIELD,2)
IL3 = SIZE(PFIELD,3)
!-------------------------------------------------------------------------------
!
!
!plmIF (HPROGRAM=='MESONH') THEN
!plm  CALL READ_SURFX3_MNH(YREC,IL1,IL2,PFIELD,KRESP,YCOMMENT,YDIR)
!plmENDIF
!
IF (HPROGRAM=='OFFLIN') THEN
#ifdef OL 
  CALL READ_SURFN_OL(YREC,PFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
!plmIF (HPROGRAM=='ASCII ') THEN
!plm  CALL READ_SURFX3_ASC(YREC,IL1,IL2,PFIELD,KRESP,YCOMMENT,YDIR)
!plmENDIF
!
!plmIF (HPROGRAM=='AROME ') THEN
!plm  CALL READ_SURFX3_ARO(YREC,IL1,IL2,PFIELD,KRESP,YCOMMENT,YDIR)
!plmENDIF
!
!plmIF (HPROGRAM=='FA    ') THEN
!plm  CALL READ_SURFX3_FA(YREC,IL1,IL2,PFIELD,KRESP,YCOMMENT,YDIR)
!plmENDIF
!
!
!-------------------------------------------------------------------------------
IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFX3',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFX3
!
!     #############################################################
      SUBROUTINE READ_SURFN0(HPROGRAM,HREC,KFIELD,KRESP,HCOMMENT)
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
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
#ifdef OL
USE MODE_READ_SURF_OL, ONLY: READ_SURF0_OL
#endif
#ifdef ASC
USE MODE_READ_SURF_ASC, ONLY: READ_SURF0_ASC
#endif
#ifdef FA
USE MODE_READ_SURF_FA, ONLY: READ_SURF0_FA
#endif
#ifdef LFI
USE MODE_READ_SURF_LFI, ONLY: READ_SURF0_LFI
#endif
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=6)       ,  INTENT(IN)  :: HPROGRAM ! calling program

CHARACTER(LEN=*)          ,INTENT(IN) ::HREC   ! name of the article to be read

INTEGER, &
                             INTENT(OUT)::KFIELD ! the integer to be read  

INTEGER                  ,INTENT(OUT) :: KRESP          ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read
!
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER(LEN=16)  :: YREC
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFN0',0,ZHOOK_HANDLE)
YREC = HREC
!
!-------------------------------------------------------------------------------
!
IF (HPROGRAM=='MESONH') THEN
#ifdef MNH
  CALL READ_SURFN0_MNH(YREC,KFIELD,KRESP,YCOMMENT)
#endif
ENDIF
!
IF (HPROGRAM=='OFFLIN') THEN
#ifdef OL 
  CALL READ_SURF0_OL(YREC,KFIELD,KRESP,YCOMMENT)
#endif
ENDIF
!
IF (HPROGRAM=='ASCII ') THEN
#ifdef ASC 
  CALL READ_SURF0_ASC(YREC,KFIELD,KRESP,YCOMMENT)
#endif
ENDIF
!
IF (HPROGRAM=='AROME ') THEN
#ifdef ARO 
  CALL READ_SURFN0_ARO(YREC,KFIELD,KRESP,YCOMMENT)
#endif
ENDIF
!
IF (HPROGRAM=='FA    ') THEN
#ifdef FA 
  CALL READ_SURF0_FA(YREC,KFIELD,KRESP,YCOMMENT)
#endif
ENDIF
!
IF (HPROGRAM=='LFI   ') THEN
#ifdef LFI 
  CALL READ_SURF0_LFI(YREC,KFIELD,KRESP,YCOMMENT)
#endif
ENDIF
!
!-------------------------------------------------------------------------------
IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFN0',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFN0
!
!     #############################################################
      SUBROUTINE READ_SURFN1(HPROGRAM,HREC,KFIELD,KRESP,HCOMMENT,HDIR)
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
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
#ifdef OL
USE MODE_READ_SURF_OL, ONLY: READ_SURFN_OL
#endif
#ifdef ASC
USE MODE_READ_SURF_ASC, ONLY: READ_SURFN_ASC
#endif
#ifdef FA
USE MODE_READ_SURF_FA, ONLY: READ_SURFN_FA
#endif
#ifdef LFI
USE MODE_READ_SURF_LFI, ONLY: READ_SURFN_LFI
#endif
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=6)       ,  INTENT(IN)  :: HPROGRAM ! calling program

CHARACTER(LEN=*)          ,INTENT(IN) ::HREC   ! name of the article to be read

INTEGER, DIMENSION(:),&
                             INTENT(OUT)::KFIELD ! the integer to be read  

INTEGER                  ,INTENT(OUT) :: KRESP          ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read
CHARACTER(LEN=1),OPTIONAL,INTENT(IN)  :: HDIR       ! type of field :
!                                                   ! 'H' : field with
!                                                   !       horizontal spatial dim.
!                                                   ! '-' : no horizontal dim.
!
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER(LEN=16)  :: YREC
CHARACTER(LEN=1)   :: YDIR
INTEGER            :: IL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFN1',0,ZHOOK_HANDLE)
YREC = HREC
YDIR = 'H'
IF (PRESENT(HDIR)) YDIR = HDIR
!
IL = SIZE(KFIELD,1)
!-------------------------------------------------------------------------------
!
IF (HPROGRAM=='MESONH') THEN
#ifdef MNH 
  CALL READ_SURFN1_MNH(YREC,IL,KFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='OFFLIN') THEN
#ifdef OL 
  CALL READ_SURFN_OL(YREC,KFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='ASCII ') THEN
#ifdef ASC 
  CALL READ_SURFN_ASC(YREC,KFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='AROME ') THEN
#ifdef ARO 
  CALL READ_SURFN1_ARO(YREC,IL,KFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='FA    ') THEN
#ifdef FA
  CALL READ_SURFN_FA(YREC,IL,KFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
IF (HPROGRAM=='LFI   ') THEN
#ifdef LFI 
  CALL READ_SURFN_LFI(YREC,KFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
!-------------------------------------------------------------------------------
IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFN1',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFN1

!     #############################################################
      SUBROUTINE READ_SURFC0(HPROGRAM,HREC,HFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READC0* - routine to read an integer
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
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
#ifdef OL
USE MODE_READ_SURF_OL, ONLY: READ_SURF0_OL
#endif
#ifdef ASC
USE MODE_READ_SURF_ASC, ONLY: READ_SURF0_ASC
#endif
#ifdef FA
USE MODE_READ_SURF_FA, ONLY: READ_SURF0_FA
#endif
#ifdef LFI
USE MODE_READ_SURF_LFI, ONLY: READ_SURF0_LFI
#endif
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=6)       ,  INTENT(IN)  :: HPROGRAM ! calling program

CHARACTER(LEN=*)          ,INTENT(IN) ::HREC   ! name of the article to be read

CHARACTER(LEN=*), &
                             INTENT(OUT)::HFIELD ! the integer to be read  

INTEGER                  ,INTENT(OUT) :: KRESP          ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER(LEN=16)  :: YREC
CHARACTER(LEN=40)  :: YFIELD
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFC0',0,ZHOOK_HANDLE)
YREC = HREC
!
!-------------------------------------------------------------------------------
!
IF (HPROGRAM=='MESONH') THEN
#ifdef MNH 
  CALL READ_SURFC0_MNH(YREC,YFIELD,KRESP,YCOMMENT)
#endif
ELSE IF (HPROGRAM=='OFFLIN') THEN
#ifdef OL 
  CALL READ_SURF0_OL(YREC,YFIELD,KRESP,YCOMMENT)
#endif
ELSE IF (HPROGRAM=='ASCII ') THEN
#ifdef ASC 
  CALL READ_SURF0_ASC(YREC,YFIELD,KRESP,YCOMMENT)
#endif
ELSE IF (HPROGRAM=='AROME ') THEN
#ifdef ARO 
  CALL READ_SURFC0_ARO(YREC,YFIELD,KRESP,YCOMMENT)
#endif
ELSE IF (HPROGRAM=='FA    ') THEN
#ifdef FA
  CALL READ_SURF0_FA(YREC,YFIELD,KRESP,YCOMMENT)
#endif
ELSE IF (HPROGRAM=='LFI   ') THEN
#ifdef LFI 
  CALL READ_SURF0_LFI(YREC,YFIELD,KRESP,YCOMMENT)
#endif
ENDIF
!
HFIELD = YFIELD(1:LEN(HFIELD))
!
!-------------------------------------------------------------------------------
IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFC0',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFC0
!
!     #############################################################
      SUBROUTINE READ_SURFL1(HPROGRAM,HREC,OFIELD,KRESP,HCOMMENT,HDIR)
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
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
#ifdef OL
USE MODE_READ_SURF_OL, ONLY: READ_SURFN_OL
#endif
#ifdef ASC
USE MODE_READ_SURF_ASC, ONLY: READ_SURFN_ASC
#endif
#ifdef FA
USE MODE_READ_SURF_FA, ONLY: READ_SURFN_FA
#endif
#ifdef LFI
USE MODE_READ_SURF_LFI, ONLY: READ_SURFN_LFI
#endif
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=6)       ,  INTENT(IN)  :: HPROGRAM ! calling program

CHARACTER(LEN=*)          ,INTENT(IN) ::HREC   ! name of the article to be read

LOGICAL, DIMENSION(:), &
                             INTENT(OUT)::OFIELD ! array containing the data field  

INTEGER                  ,INTENT(OUT) :: KRESP          ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read
CHARACTER(LEN=1),OPTIONAL,INTENT(IN)  :: HDIR       ! type of field :
!                                                   ! 'H' : field with
!                                                   !       horizontal spatial dim.
!                                                   ! '-' : no horizontal dim.
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER(LEN=16)  :: YREC
CHARACTER(LEN=1)   :: YDIR
INTEGER            :: IL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFL1',0,ZHOOK_HANDLE)
YREC = HREC
YDIR = 'H'
IF (PRESENT(HDIR)) YDIR = HDIR
!
IL = SIZE(OFIELD)
!-------------------------------------------------------------------------------
!
IF (HPROGRAM=='MESONH') THEN
#ifdef MNH 
  CALL READ_SURFL1_MNH(YREC,IL,OFIELD,KRESP,YCOMMENT,YDIR)
#endif
ELSE IF (HPROGRAM=='OFFLIN') THEN
#ifdef OL 
  CALL READ_SURFN_OL(YREC,OFIELD,KRESP,YCOMMENT,YDIR)
#endif
ELSE IF (HPROGRAM=='ASCII ') THEN
#ifdef ASC 
  CALL READ_SURFN_ASC(YREC,OFIELD,KRESP,YCOMMENT,YDIR)
#endif
ELSE IF (HPROGRAM=='AROME ') THEN
#ifdef ARO 
  CALL READ_SURFL1_ARO(YREC,IL,OFIELD,KRESP,YCOMMENT,YDIR)
#endif
ELSE IF (HPROGRAM=='FA    ') THEN
#ifdef FA 
  CALL READ_SURFN_FA(YREC,IL,OFIELD,KRESP,YCOMMENT,YDIR)
#endif
ELSE IF (HPROGRAM=='LFI   ') THEN
#ifdef LFI 
  CALL READ_SURFN_LFI(YREC,OFIELD,KRESP,YCOMMENT,YDIR)
#endif
ENDIF
!
!-------------------------------------------------------------------------------
IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFL1',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFL1
!
!     #############################################################
      SUBROUTINE READ_SURFL0(HPROGRAM,HREC,OFIELD,KRESP,HCOMMENT)
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
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
#ifdef OL
USE MODE_READ_SURF_OL, ONLY: READ_SURF0_OL
#endif
#ifdef ASC
USE MODE_READ_SURF_ASC, ONLY: READ_SURF0_ASC
#endif
#ifdef FA
USE MODE_READ_SURF_FA, ONLY: READ_SURF0_FA
#endif
#ifdef LFI
USE MODE_READ_SURF_LFI, ONLY: READ_SURF0_LFI
#endif
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=6)       ,  INTENT(IN)  :: HPROGRAM ! calling program

CHARACTER(LEN=*)          ,INTENT(IN) ::HREC   ! name of the article to be read

LOGICAL,                  INTENT(OUT)::OFIELD ! array containing the data field

INTEGER                  ,INTENT(OUT) :: KRESP          ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER(LEN=16)  :: YREC
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFL0',0,ZHOOK_HANDLE)
YREC = HREC
!
!-------------------------------------------------------------------------------
!
IF (HPROGRAM=='MESONH') THEN
#ifdef MNH 
  CALL READ_SURFL0_MNH(YREC,OFIELD,KRESP,YCOMMENT)
#endif
ELSE IF (HPROGRAM=='OFFLIN') THEN
#ifdef OL 
  CALL READ_SURF0_OL(YREC,OFIELD,KRESP,YCOMMENT)
#endif
ELSE IF (HPROGRAM=='ASCII ') THEN
#ifdef ASC 
  CALL READ_SURF0_ASC(YREC,OFIELD,KRESP,YCOMMENT)
#endif
ELSE IF (HPROGRAM=='AROME ') THEN
#ifdef ARO 
  CALL READ_SURFL0_ARO(YREC,OFIELD,KRESP,YCOMMENT)
#endif
ELSE IF (HPROGRAM=='FA    ') THEN
#ifdef FA 
  CALL READ_SURF0_FA(YREC,OFIELD,KRESP,YCOMMENT)
#endif
ELSE IF (HPROGRAM=='LFI   ') THEN
#ifdef LFI 
  CALL READ_SURF0_LFI(YREC,OFIELD,KRESP,YCOMMENT)
#endif
ENDIF
!
!-------------------------------------------------------------------------------
IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFL0',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURFL0

!     #############################################################
      SUBROUTINE READ_SURFT0(HPROGRAM,HREC,TFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READT0* - routine to read a MESO-NH date_time scalar
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
USE MODI_GET_LUOUT
USE MODD_TYPE_DATE_SURF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
#ifdef OL
USE MODE_READ_SURF_OL, ONLY: READ_SURFT_OL
#endif
#ifdef ASC
USE MODE_READ_SURF_ASC, ONLY: READ_SURFT_ASC
#endif
#ifdef FA
USE MODE_READ_SURF_FA, ONLY: READ_SURFT_FA
#endif
#ifdef LFI
USE MODE_READ_SURF_LFI, ONLY: READ_SURFT_LFI
#endif
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=6)       ,  INTENT(IN)  :: HPROGRAM ! calling program

CHARACTER(LEN=*)          ,INTENT(IN) ::HREC   ! name of the article to be read
TYPE (DATE_TIME), &
                           INTENT(INOUT)::TFIELD ! array containing the data field  

INTEGER                  ,INTENT(OUT) :: KRESP          ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read

!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER(LEN=16)  :: YREC
INTEGER            :: ILUOUT
!
REAL    :: ZTIME
INTEGER :: IDAY
INTEGER :: IMONTH
INTEGER :: IYEAR
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFT0',0,ZHOOK_HANDLE)
YREC = HREC
!
!-------------------------------------------------------------------------------
!
IF (HPROGRAM=='MESONH') THEN
#ifdef MNH 
  CALL READ_SURFT0_MNH(YREC,IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)
#endif
ELSE IF (HPROGRAM=='OFFLIN') THEN
#ifdef OL 
  CALL READ_SURFT_OL(YREC,IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)
#endif
ELSE IF (HPROGRAM=='ASCII ') THEN
#ifdef ASC 
  CALL READ_SURFT_ASC(YREC,IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)
#endif
ELSE IF (HPROGRAM=='AROME ') THEN
#ifdef ARO 
  CALL READ_SURFT0_ARO(YREC,IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)
#endif
ELSE IF (HPROGRAM=='FA    ') THEN
#ifdef FA 
  CALL READ_SURFT_FA(YREC,IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)
#endif
ELSE IF (HPROGRAM=='LFI   ') THEN
#ifdef LFI 
  CALL READ_SURFT_LFI(YREC,IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)
#endif
ENDIF
!
IF (KRESP==-2) THEN
  CALL GET_LUOUT(HPROGRAM,ILUOUT)
  WRITE(ILUOUT,*) '-------'
  WRITE(ILUOUT,*) 'WARNING'
  WRITE(ILUOUT,*) '-------'
  WRITE(ILUOUT,*) ' '
  WRITE(ILUOUT,*) 'Date is not present file'
  WRITE(ILUOUT,*) 'Forcing value is kept'
  WRITE(ILUOUT,*) ' '
ELSE
  TFIELD%TDATE=DATE(IYEAR,IMONTH,IDAY)  
  TFIELD%TIME =ZTIME
END IF
!-------------------------------------------------------------------------------
IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFT0',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_SURFT0

!     #############################################################
      SUBROUTINE READ_SURFT1(HPROGRAM,HREC,TFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READT2* - routine to read a MESO-NH date_time 1d array
!!
!!    PURPOSE
!!    -------
!
!       The purpose of READT1 is
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
!!      P. Le Moigne      *METEO-FRANCE*
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      original                                                     09/2007
!!      G. TANGUY 03/2009 : add CALL READ_SURFT1_MNH
!!      E.BAZILE 10/2010 appel à read_surft1_aro
!----------------------------------------------------------------------------
!
!*      0.    DECLARATIONS
!             ------------
!
USE MODI_GET_LUOUT
USE MODD_TYPE_DATE_SURF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
#ifdef ASC
USE MODE_READ_SURF_ASC, ONLY: READ_SURFT_ASC
#endif
#ifdef LFI
USE MODE_READ_SURF_LFI, ONLY: READ_SURFT_LFI
#endif
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=6)       ,  INTENT(IN)  :: HPROGRAM ! calling program

CHARACTER(LEN=*)          ,INTENT(IN) ::HREC   ! name of the article to be read
TYPE (DATE_TIME), DIMENSION(:), &
                           INTENT(INOUT)::TFIELD ! array containing the data field  

INTEGER                  ,INTENT(OUT) :: KRESP          ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read

!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER(LEN=16)  :: YREC
INTEGER            :: ILUOUT
!
INTEGER :: IL1
REAL , DIMENSION(:), ALLOCATABLE   :: ZTIME
INTEGER, DIMENSION(:), ALLOCATABLE :: IDAY
INTEGER, DIMENSION(:), ALLOCATABLE :: IMONTH
INTEGER, DIMENSION(:), ALLOCATABLE :: IYEAR
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFT1',0,ZHOOK_HANDLE)
YREC = HREC
!
IL1=SIZE(TFIELD,1)
!
ALLOCATE(ZTIME (IL1))
ALLOCATE(IDAY  (IL1))
ALLOCATE(IMONTH(IL1))
ALLOCATE(IYEAR (IL1))
!-------------------------------------------------------------------------------
!
IF (HPROGRAM=='MESONH') THEN
#ifdef MNH 
  CALL READ_SURFT1_MNH(YREC,IL1,IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)
#endif
ELSE IF (HPROGRAM=='OFFLIN') THEN
  CALL ABOR1_SFX('READ_SURFT1: NOT AVAILABLE FOR OFFLIN')
ELSE IF (HPROGRAM=='ASCII ') THEN
#ifdef ASC 
  CALL READ_SURFT_ASC(YREC,IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)
#endif
ELSE IF (HPROGRAM=='AROME ') THEN
#ifdef ARO        
  CALL READ_SURFT1_ARO(YREC,IL1,IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)
#endif        
ELSE IF (HPROGRAM=='FA    ') THEN
  CALL ABOR1_SFX('READ_SURFT1: NOT AVAILABLE FOR FA')
ELSE IF (HPROGRAM=='LFI   ') THEN
#ifdef ASC 
  CALL READ_SURFT_LFI(YREC,IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)
#endif        
ENDIF
!
IF (KRESP==-2) THEN
  CALL GET_LUOUT(HPROGRAM,ILUOUT)
  WRITE(ILUOUT,*) '-------'
  WRITE(ILUOUT,*) 'WARNING'
  WRITE(ILUOUT,*) '-------'
  WRITE(ILUOUT,*) ' '
  WRITE(ILUOUT,*) 'Date is not present file'
  WRITE(ILUOUT,*) 'Forcing value is kept'
  WRITE(ILUOUT,*) ' '
ELSE
  TFIELD(:)%TDATE%YEAR  = IYEAR (:)
  TFIELD(:)%TDATE%MONTH = IMONTH(:)
  TFIELD(:)%TDATE%DAY   = IDAY  (:)  
  TFIELD(:)%TIME        = ZTIME (:)        
END IF
!
DEALLOCATE(ZTIME)
DEALLOCATE(IDAY)
DEALLOCATE(IMONTH)
DEALLOCATE(IYEAR)
!
!-------------------------------------------------------------------------------
IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFT1',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_SURFT1

!     #############################################################
      SUBROUTINE READ_SURFT2(HPROGRAM,HREC,TFIELD,KRESP,HCOMMENT)
!     #############################################################
!
!!****  *READT2* - routine to read a MESO-NH date_time 2d array
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
USE MODI_GET_LUOUT
USE MODD_TYPE_DATE_SURF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
#ifdef ASC
USE MODE_READ_SURF_ASC, ONLY: READ_SURFT_ASC
#endif
#ifdef FA
USE MODE_READ_SURF_FA, ONLY: READ_SURFT_FA
#endif
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=6)       ,  INTENT(IN)  :: HPROGRAM ! calling program

CHARACTER(LEN=*)          ,INTENT(IN) ::HREC   ! name of the article to be read
TYPE (DATE_TIME), DIMENSION(:,:), &
                           INTENT(INOUT)::TFIELD ! array containing the data field  

INTEGER                  ,INTENT(OUT) :: KRESP          ! KRESP  : return-code if a problem appears
CHARACTER(LEN=*),OPTIONAL,INTENT(OUT) :: HCOMMENT   ! name of the article to be read

!*      0.2   Declarations of local variables
!
CHARACTER(LEN=100) :: YCOMMENT
CHARACTER(LEN=16)  :: YREC
INTEGER            :: ILUOUT
!
INTEGER :: IL1, IL2
REAL , DIMENSION(:,:), ALLOCATABLE   :: ZTIME
INTEGER, DIMENSION(:,:), ALLOCATABLE :: IDAY
INTEGER, DIMENSION(:,:), ALLOCATABLE :: IMONTH
INTEGER, DIMENSION(:,:), ALLOCATABLE :: IYEAR
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFT2',0,ZHOOK_HANDLE)
YREC = HREC
!
IL1=SIZE(TFIELD,1)
IL2=SIZE(TFIELD,2)
!
ALLOCATE(ZTIME (IL1,IL2))
ALLOCATE(IDAY  (IL1,IL2))
ALLOCATE(IMONTH(IL1,IL2))
ALLOCATE(IYEAR (IL1,IL2))
!-------------------------------------------------------------------------------
!
IF (HPROGRAM=='MESONH') THEN
  CALL ABOR1_SFX('READ_SURFT2: NOT AVAILABLE FOR MESONH')
ELSE IF (HPROGRAM=='OFFLIN') THEN
  CALL ABOR1_SFX('READ_SURFT2: NOT AVAILABLE FOR OFFLIN')
ELSE IF (HPROGRAM=='ASCII ') THEN
#ifdef ASC 
  CALL READ_SURFT_ASC(YREC,IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)
#endif
ELSE IF (HPROGRAM=='AROME ') THEN
  CALL ABOR1_SFX('READ_SURFT2: NOT AVAILABLE FOR AROME')
ELSE IF (HPROGRAM=='FA    ') THEN
#ifdef FA 
  CALL READ_SURFT_FA(YREC,IL1,IL2,IYEAR,IMONTH,IDAY,ZTIME,KRESP,YCOMMENT)
#endif
ELSE IF (HPROGRAM=='LFI   ') THEN
  CALL ABOR1_SFX('READ_SURFT2: NOT AVAILABLE FOR LFI')
ENDIF

!
IF (KRESP==-2) THEN
  CALL GET_LUOUT(HPROGRAM,ILUOUT)
  WRITE(ILUOUT,*) '-------'
  WRITE(ILUOUT,*) 'WARNING'
  WRITE(ILUOUT,*) '-------'
  WRITE(ILUOUT,*) ' '
  WRITE(ILUOUT,*) 'Date is not present file'
  WRITE(ILUOUT,*) 'Forcing value is kept'
  WRITE(ILUOUT,*) ' '
ELSE
  TFIELD(:,:)%TDATE%YEAR  = IYEAR (:,:)
  TFIELD(:,:)%TDATE%MONTH = IMONTH(:,:)
  TFIELD(:,:)%TDATE%DAY   = IDAY  (:,:)  
  TFIELD(:,:)%TIME        = ZTIME (:,:)        
END IF
!
DEALLOCATE(ZTIME)
DEALLOCATE(IDAY)
DEALLOCATE(IMONTH)
DEALLOCATE(IYEAR)
!
!-------------------------------------------------------------------------------
IF (PRESENT(HCOMMENT)) HCOMMENT = YCOMMENT
IF (LHOOK) CALL DR_HOOK('MODI_READ_SURF:READ_SURFT2',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE READ_SURFT2
