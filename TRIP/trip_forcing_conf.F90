!######
SUBROUTINE TRIP_FORCING_CONF(KLISTING,KYEAR,KMONTH,KDAY,PTIME,  &
                             HFILE_FRC,HREADFRC,HDRAIN,HRUNOFF, &
                             KLON,KLAT,PTSTEP_RUN,KNB_TSTEP_RUN,&
                             PRUNTIME                          )
!#######################################################################
!
!!****  *TRIP_FORCING_CONF* - prepare the dimenssions (xt or xyt) of run
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
!!      B. decharme   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/2008 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_TRIP_PAR,    ONLY : XDAY
!
USE MODI_ABORT_TRIP
USE MODI_READ_DIMLEN
USE MODI_READ_FORCING_DATE
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER,                INTENT(IN)  :: KLISTING
INTEGER,                INTENT(IN)  :: KYEAR
INTEGER,                INTENT(IN)  :: KMONTH
INTEGER,                INTENT(IN)  :: KDAY
REAL,                   INTENT(IN)  :: PTIME
!
 CHARACTER(LEN=15),      INTENT(IN)  :: HFILE_FRC
 CHARACTER(LEN=6),       INTENT(IN)  :: HREADFRC 
 CHARACTER(LEN=8),       INTENT(IN)  :: HDRAIN 
 CHARACTER(LEN=8),       INTENT(IN)  :: HRUNOFF 
INTEGER,                INTENT(IN)  :: KLON
INTEGER,                INTENT(IN)  :: KLAT
REAL,                   INTENT(IN)  :: PTSTEP_RUN
!
INTEGER,                INTENT(OUT) :: KNB_TSTEP_RUN
REAL,                   INTENT(OUT) :: PRUNTIME
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER, DIMENSION (:), ALLOCATABLE :: IDIMLEN
!
INTEGER :: INDIM
INTEGER :: ILON
INTEGER :: ILAT
INTEGER :: INI
INTEGER :: IWORK1, IWORK2
!
INTEGER :: IYEAR
INTEGER :: IMONTH
INTEGER :: IDAY
REAL    :: ZTIME
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
! Read the configuration of the run
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TRIP_FORCING_CONF',0,ZHOOK_HANDLE)
!
IF(HREADFRC=='VECTOR')THEN
  INDIM=2
ELSE
  INDIM=3
ENDIF 
!
ALLOCATE(IDIMLEN(INDIM))
IDIMLEN=0
!
 CALL READ_DIMLEN(KLISTING,HFILE_FRC,HDRAIN,INDIM,IDIMLEN)
!
IF(HREADFRC=='VECTOR')THEN
!        
  INI    = IDIMLEN(1)
  IWORK1 = IDIMLEN(2)
!  
  IF(INI/=KLON*KLAT)THEN
     WRITE(KLISTING,*)'For the variable DRAIN : '
     WRITE(KLISTING,*)'TRIP_FORCING_CONF : The number of points in the forcing variables (',INI,')'
     WRITE(KLISTING,*)'are /= of the number of points in the TRIP domain'
     WRITE(KLISTING,*)'NLON = ',KLON,' and NLAT = ',KLAT,' so the number of points = ',KLON*KLAT
     CALL ABORT_TRIP('TRIP_FORCING_CONF: number of points in the forcing variables are /= than TRIP domain')
  ENDIF
!
ELSE
!        
  ILON   = IDIMLEN(1)
  ILAT   = IDIMLEN(2)
  IWORK1 = IDIMLEN(3)
  INI    = ILON*ILAT
!  
  IF(INI/=KLON*KLAT)THEN
     WRITE(KLISTING,*)'For the variable DRAIN : '
     WRITE(KLISTING,*)'TRIP_FORCING_CONF : The number of points in the forcing variables (',ILON*ILAT,')'
     WRITE(KLISTING,*)'are /= of the number of points in the TRIP domain'
     WRITE(KLISTING,*)'NLON = ',KLON,' and NLAT = ',KLAT,' while FRC_LON = ',ILON,' and FRC_LAT = ',ILAT
     CALL ABORT_TRIP('TRIP_FORCING_CONF: number of points in the forcing variables are /= than TRIP domain')
  ENDIF
!        
ENDIF 
!
IDIMLEN=0
!
 CALL READ_DIMLEN(KLISTING,HFILE_FRC,HRUNOFF,INDIM,IDIMLEN)
!
IF(HREADFRC=='VECTOR')THEN
!        
  INI    = IDIMLEN(1)
  IWORK2 = IDIMLEN(2)
!  
  IF(INI/=KLON*KLAT)THEN
     WRITE(KLISTING,*)'For the variable RUNOFF : '
     WRITE(KLISTING,*)'TRIP_FORCING_CONF : The number of points in the forcing variables (',INI,')'
     WRITE(KLISTING,*)'are /= of the number of points in the TRIP domain'
     WRITE(KLISTING,*)'NLON = ',KLON,' and NLAT = ',KLAT,' so the number of points = ',KLON*KLAT
     CALL ABORT_TRIP('TRIP_FORCING_CONF: number of points in the forcing variables are /= than TRIP domain')
  ENDIF
!
ELSE
!        
  ILON   = IDIMLEN(1)
  ILAT   = IDIMLEN(2)
  IWORK2 = IDIMLEN(3)
  INI    = ILON*ILAT
!  
  IF(INI/=KLON*KLAT)THEN
     WRITE(KLISTING,*)'For the variable RUNOFF : '
     WRITE(KLISTING,*)'TRIP_FORCING_CONF : The number of points in the forcing variables (',ILON*ILAT,')'
     WRITE(KLISTING,*)'are different than the number of points in the TRIP domain'
     WRITE(KLISTING,*)'NLON = ',KLON,' and NLAT = ',KLAT,' while FRC_LON = ',ILON,' and FRC_LAT = ',ILAT
     CALL ABORT_TRIP('TRIP_FORCING_CONF: number of points in the forcing variables are /= than TRIP domain')
  ENDIF
!        
ENDIF 
!
DEALLOCATE(IDIMLEN)
!
!-------------------------------------------------------------------------------
! Configuration of the run
!-------------------------------------------------------------------------------
!
! * Number of time step during the run
!
IF(IWORK1/=IWORK2)THEN
  WRITE(KLISTING,*)'TRIP_FORCING_CONF : The number of time step are different for each forcing variable !'
  WRITE(KLISTING,*)'NB_TSTEP for DRAIN = ',IWORK1,' while NB_TSTEP for RUNOFF = ',IWORK2
  CALL ABORT_TRIP('TRIP_FORCING_CONF: The number of time step are different for each forcing variable')
ELSE
  KNB_TSTEP_RUN = IWORK1
ENDIF
!
PRUNTIME = REAL(KNB_TSTEP_RUN) * PTSTEP_RUN
!
! * Date the run
!
 CALL READ_FORCING_DATE(KLISTING,HFILE_FRC,KNB_TSTEP_RUN, &
                       IYEAR,IMONTH,IDAY,ZTIME           )
!
IF ( (KYEAR /= IYEAR) .OR. (KMONTH /= IMONTH) .OR. (KDAY /= IDAY) .OR. (PTIME /= ZTIME)) THEN
  WRITE(KLISTING,*)' DATE INCONSISTENCY: RESTART FILE = ',KYEAR,KMONTH,KDAY,PTIME
  WRITE(KLISTING,*)' DATE INCONSISTENCY: FORCING FILE = ',IYEAR,IMONTH,IDAY,ZTIME
  CALL ABORT_TRIP('TRIP_FORCING_CONF: DATE INCONSISTENCY')
ENDIF
!
IF (LHOOK) CALL DR_HOOK('TRIP_FORCING_CONF',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE TRIP_FORCING_CONF
