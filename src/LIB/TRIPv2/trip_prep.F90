!######
PROGRAM TRIP_PREP
!###########
!
!!****  *TRIP_PREP* - driver for TRIP fields preparation
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     B. Decharme 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/2008
!!------------------------------------------------------------------
!
USE MODD_TRIP_GRID, ONLY : TPG => TRIP_GRID
!
USE MODD_TRIP_LISTING
!
USE MODI_GOTO_TRIP
USE MODI_INIT_TRIP_PAR
USE MODI_PREP_TRIP_RUN
!
USE MODI_OPEN_TRIP_NAMELIST
USE MODI_CLOSE_TRIP_NAMELIST
USE MODI_ABORT_TRIP
USE MODI_TRIP_POSNAM
USE MODI_READ_NAM_TRIP
USE MODI_READ_NAM_TRIP_GRID
!
USE MODI_ALLOC_TRIP
USE MODI_DEALLOC_TRIP
USE MODI_TRIP_OASIS_END
!
USE MODI_TRIP_OASIS_INIT
USE MODI_TRIP_OASIS_READ_NAM
USE MODI_TRIP_OASIS_PREP
USE MODI_TRIP_OASIS_END
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
#ifdef TRIPOASIS
INCLUDE 'mpif.h'
#endif
!
INTEGER  :: ILON
INTEGER  :: ILAT
!
INTEGER  :: NYEAR           ! current year  (UTC)
INTEGER  :: NMONTH          ! current month (UTC)
INTEGER  :: NDAY            ! current day   (UTC)
REAL     :: XTIME           ! current time    (s)
!
INTEGER  :: ILUNAM          ! namelist unit number
LOGICAL  :: GFOUND          ! return logical when reading namelist
INTEGER  :: ILOCAL_COMM     ! Local communicator
LOGICAL  :: GOASIS          ! OASIS used(default=.false.)
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
NAMELIST/NAM_START_DATE/NYEAR,NMONTH,NDAY,XTIME
!
! --------------------------------------------------------------------------------------
! * 0. MPI and OASIS must be initialized before any DR_HOOK call
! --------------------------------------------------------------------------------------
!
CALL TRIP_OASIS_INIT(GOASIS,ILOCAL_COMM)
!
IF (LHOOK) CALL DR_HOOK('TRIP_PREP',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
CALL ALLOC_TRIP(1)
CALL GOTO_TRIP (1)
!-------------------------------------------------------------------------------
!
CLISTING = CLISTING_PREP
!
OPEN (UNIT=NLISTING,FILE=CLISTING_PREP,FORM='FORMATTED',ACTION='WRITE')
!
! --------------------------------------------------------------------------------------
!* 1. Time configuration: start date of the run
! --------------------------------------------------------------------------------------
!
CALL OPEN_TRIP_NAMELIST(ILUNAM)
CALL TRIP_POSNAM(ILUNAM,'NAM_START_DATE',GFOUND,NLISTING)
IF (GFOUND) THEN
   READ (UNIT=ILUNAM,NML=NAM_START_DATE)
ELSE
  WRITE(NLISTING,*)'NAM_START_DATE not found in namelist'
  WRITE(NLISTING,*)'NYEAR, NMONTH, NDAY and XTIME must be initialized'
  WRITE(NLISTING,*)'as the date of the beginning of the run'
  CALL ABORT_TRIP('NAM_START_DATE not found in namelist')        
ENDIF
CALL CLOSE_TRIP_NAMELIST(ILUNAM)
!
! --------------------------------------------------------------------------------------
!* 2. Initializations
! --------------------------------------------------------------------------------------
!
CALL INIT_TRIP_PAR
!
! --------------------------------------------------------------------------------------
!* 3. Read grid and physical option namelists
! --------------------------------------------------------------------------------------
!
CALL READ_NAM_TRIP(NLISTING)
!
CALL READ_NAM_TRIP_GRID(TPG, &
                        NLISTING)
!
IF(GOASIS)THEN
  CALL TRIP_OASIS_READ_NAM(NLISTING)
ENDIF
!
! --------------------------------------------------------------------------------------
!* 4. TRIP parameters preparation
! --------------------------------------------------------------------------------------
!
CALL PREP_TRIP_RUN(NYEAR,NMONTH,NDAY,XTIME,ILON,ILAT)
!
IF(GOASIS)THEN
  CALL TRIP_OASIS_PREP(TPG, &
                       NLISTING,ILON,ILAT)
ENDIF
!
! --------------------------------------------------------------------------------------
!
WRITE(NLISTING,*) ' '
WRITE(NLISTING,*) '    ----------------------------'
WRITE(NLISTING,*) '    | TRIP PREP ENDS CORRECTLY |'
WRITE(NLISTING,*) '    ----------------------------'
!
WRITE(*,*) ' '
WRITE(*,*) '    ----------------------------'
WRITE(*,*) '    | TRIP PREP ENDS CORRECTLY |'
WRITE(*,*) '    ----------------------------'
!
CLOSE(NLISTING)
!
!-------------------------------------------------------------------------------
CALL DEALLOC_TRIP
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TRIP_PREP',1,ZHOOK_HANDLE)
!
! --------------------------------------------------------------------------------------
! * 3. MPI and OASIS must be finalized after the last DR_HOOK call
! --------------------------------------------------------------------------------------
!
IF(GOASIS)THEN
  CALL TRIP_OASIS_END
ENDIF
!
!-------------------------------------------------------------------------------------
!
END PROGRAM TRIP_PREP
