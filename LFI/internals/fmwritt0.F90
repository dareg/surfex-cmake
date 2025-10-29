!     ######spl
      SUBROUTINE FMWRITT0(HFILEM,HRECFM,HFIPRI,KLENG,TFIELD,KGRID,&
                           KLENCH,HCOMMENT,KRESP)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!     #############################################################
!
!!****  *FMWRITT0* - routine to write a date scalar into a "FM"-file
!!
!!    PURPOSE
!!    -------
!
!       The purpose of FMWRITT0 is to split a date_time scalar
!      and to call FM_WRIT without interface module
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!      FM_WRIT
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
USE MODD_TYPE_DATE_SURF, ONLY: DATE_TIME
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=*)          ,INTENT(IN) ::HFILEM   ! file name
CHARACTER(LEN=*)          ,INTENT(IN) ::HRECFM   ! name of the article to be written

CHARACTER(LEN=*)          ,INTENT(IN) ::HFIPRI   ! file for prints in FM

INTEGER,                   INTENT(IN) ::KLENG  ! length of the data field
TYPE (DATE_TIME), &
                           INTENT(IN) ::TFIELD ! array containing the data field
INTEGER,                   INTENT(IN) ::KGRID  ! C-grid indicator (u,v,w,T)
INTEGER,                   INTENT(IN) ::KLENCH ! length of comment string

CHARACTER(LEN=*)          ,INTENT(IN) ::HCOMMENT ! comment string

INTEGER,                   INTENT(OUT)::KRESP  ! return-code if problems araised
!
!*      0.2   Declarations of local variables
!
INTEGER(KIND=8), DIMENSION(3)  :: ITDATE    ! date array
CHARACTER(LEN=12)              :: YRECFM    ! Name of the article to be written
CHARACTER(LEN=JPXKRK)          :: YCOMMENT  ! Comment string
!
!-------------------------------------------------------------------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('FMWRITT0',0,ZHOOK_HANDLE)
YRECFM=TRIM(HRECFM)//'%TDATE'   ! array of rank 3 for date is written in file
YCOMMENT='YYYYMMDD'
ITDATE(1)=TFIELD%TDATE%YEAR
ITDATE(2)=TFIELD%TDATE%MONTH
ITDATE(3)=TFIELD%TDATE%DAY
CALL FM_WRIT(HFILEM,YRECFM,HFIPRI,3,ITDATE,0,8,YCOMMENT,KRESP)
!
YRECFM=TRIM(HRECFM)//'%TIME'
YCOMMENT='SECONDS'
CALL FM_WRIT(HFILEM,YRECFM,HFIPRI,1,TFIELD%TIME,0,7,YCOMMENT,KRESP)
!
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('FMWRITT0',1,ZHOOK_HANDLE)
END SUBROUTINE FMWRITT0
