!     #########
      SUBROUTINE INTERPOL_SST_MTH(KYEAR,KMONTH,KDAY,HFLAG,POUT)
!     #######################################################
!
!!****  *INTERPOL_SST_MTH* - Interpolation of monthly SST, SSS, SIT or SIC
!!
!!    PURPOSE
!!    -------
!!
!!
!!**  METHOD
!!    ------
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      
!     B.Decharme  Meteo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    28/01/10
!!      Modified    02/2014   S. Senesi : allow to work on SSS, SIT and SIC fields
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SEAFLUX_n, ONLY : S => SEAFLUX
!
USE MODI_INTERPOL_QUADRA
!
USE MODI_ABOR1_SFX
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!------------------------
! 
INTEGER,          INTENT(IN) :: KYEAR  ! year of date
INTEGER,          INTENT(IN) :: KMONTH ! month of date
INTEGER,          INTENT(IN) :: KDAY   ! day of date
CHARACTER(LEN=1), INTENT(IN) :: HFLAG  ! 'T' for SST, 'S' for SSS, 'H' for SIT, 'C' for SIC
!
REAL, DIMENSION(:), INTENT(OUT) :: POUT   ! Sea surface temperature or salinity, or SIC or SIT at time t 
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
REAL, DIMENSION(SIZE(POUT)) :: ZFIELD ! Field at time t 
!
REAL            :: ZDAT,ZNDAT
INTEGER         :: IMTH1,IMTH2,IMTH3
INTEGER         :: INDAYS ! number of days in KMONTH
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!*       1.    Number of days in a month
!              -------------------------
!
IF (LHOOK) CALL DR_HOOK('INTERPOL_SST_MTH',0,ZHOOK_HANDLE)
IF ( (HFLAG/='S') .AND. (HFLAG/='T') .AND. (HFLAG/='H') .AND. (HFLAG/='C') )THEN
   CALL ABOR1_SFX('FATAL ERROR in INTERPOL_SST_MTH : HFLAG not S nor T nor C nor H. !')
ENDIF
SELECT CASE (KMONTH)
    CASE(4,6,9,11)
      INDAYS=30
    CASE(1,3,5,7:8,10,12)
      INDAYS=31
    CASE(2)
      IF( ((MOD(KYEAR,4)==0).AND.(MOD(KYEAR,100)/=0)) .OR. (MOD(KYEAR,400)==0))THEN
        INDAYS=29
      ELSE
        INDAYS=28
      ENDIF
END SELECT
!
!
!-------------------------------------------------------------------------------
!
!*       2.    SST or SSS  Interpolation using previous, current and next month
!              --------------------------------------------------------
!
ZDAT = REAL(KDAY)
ZNDAT= REAL(INDAYS)
!
! The current month correspond to index 2 (or KMONTH+1 if ANNUAL)
!
IF (((HFLAG=='T').AND.(TRIM(S%CINTERPOL_SST)=='MONTH')) .OR. &
    ((HFLAG=='S').AND.(TRIM(S%CINTERPOL_SSS)=='MONTH')) .OR. &
    ((HFLAG=='H').AND.(TRIM(S%CINTERPOL_SIT)=='MONTH')) .OR. &
    ((HFLAG=='C').AND.(TRIM(S%CINTERPOL_SIC)=='MONTH'))) THEN
   IMTH1=1
   IMTH2=2
   IMTH3=3
ELSE
  IMTH1=KMONTH
  IMTH2=KMONTH+1
  IMTH3=KMONTH+2
ENDIF

IF (HFLAG .EQ. 'T' ) THEN 
   CALL INTERPOL_QUADRA(ZDAT,ZNDAT,S%XSST_MTH(:,IMTH1),S%XSST_MTH(:,IMTH2),S%XSST_MTH(:,IMTH3),ZFIELD)
   POUT(:) = ZFIELD(:)
ELSEIF (HFLAG .EQ. 'S' ) THEN 
   CALL INTERPOL_QUADRA(ZDAT,ZNDAT,S%XSSS_MTH(:,IMTH1),S%XSSS_MTH(:,IMTH2),S%XSSS_MTH(:,IMTH3),ZFIELD)
   POUT(:) = MAX(0.0,ZFIELD(:))
ELSEIF (HFLAG .EQ. 'H' ) THEN 
   CALL INTERPOL_QUADRA(ZDAT,ZNDAT,S%XSIT_MTH(:,IMTH1),S%XSIT_MTH(:,IMTH2),S%XSIT_MTH(:,IMTH3),ZFIELD)
   POUT(:) = MAX(0.0,ZFIELD(:))
ELSE ! SIC
   CALL INTERPOL_QUADRA(ZDAT,ZNDAT,S%XSIC_MTH(:,IMTH1),S%XSIC_MTH(:,IMTH2),S%XSIC_MTH(:,IMTH3),ZFIELD)
   POUT(:) = MAX(0.0,MIN(1.0,ZFIELD(:)))
ENDIF
IF (LHOOK) CALL DR_HOOK('INTERPOL_SST_MTH',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INTERPOL_SST_MTH
