!     ######spl
      SUBROUTINE FMREADT0(HFILEM,HRECFM,HFIPRI,KLENG,TFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!     #############################################################
!
!!****  *FMREADT0* - routine to read a date_time scalar into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMREADT0 is to call FM_READ without interface module
!      and to retrieve the date_time information
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_READ
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
USE MODD_FMDECLAR
USE MODD_TYPE_DATE_SURF, ONLY: DATE_TIME, DATE
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
TYPE (DATE_TIME), &
                           INTENT(OUT)::TFIELD ! array containing the data field
INTEGER,                   INTENT(OUT)::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(OUT)::KLENCH ! length of comment string

CHARACTER(LEN=*)          ,INTENT(OUT)::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
CHARACTER(LEN=12)              :: YRECFM    ! Name of the article to be read
CHARACTER(LEN=JPXKRK)          :: YCOMMENT 
INTEGER(KIND=8), DIMENSION(3)  :: ITDATE
REAL(KIND=8)                   :: ZFIELD
!-------------------------------------------------------------------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('FMREADT0',0,ZHOOK_HANDLE)
YRECFM=TRIM(HRECFM)//'%TDATE'
CALL FM_READ(HFILEM,YRECFM,HFIPRI,3,ITDATE,KGRID,KLENCH,YCOMMENT,KRESP)
TFIELD%TDATE=DATE(ITDATE(1),ITDATE(2),ITDATE(3))  
HCOMMENT = YCOMMENT
!
YRECFM=TRIM(HRECFM)//'%TIME'
CALL FM_READ(HFILEM,YRECFM,HFIPRI,1,ZFIELD,KGRID,KLENCH,YCOMMENT,KRESP)
TFIELD%TIME=ZFIELD
HCOMMENT = YCOMMENT
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('FMREADT0',1,ZHOOK_HANDLE)
END SUBROUTINE FMREADT0
