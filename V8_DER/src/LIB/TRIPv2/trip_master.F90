!###################################################################
PROGRAM TRIP_MASTER
!###################################################################
!
!!****  *TRIP_MASTER*  
!!
!!    PURPOSE
!!    -------
!!   
!!    Driver for TRIP from CNRM
!!      
!!    REFERENCE
!!    ---------
!!      
!!    AUTHOR
!!    ------
!!      B. Decharme
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/08 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_TRIP, ONLY : TP => TRIP
!
USE MODD_TRIP_DIAG, ONLY : TPDG => TRIP_DIAG
!
USE MODD_TRIP_GRID, ONLY : TPG => TRIP_GRID
!
USE MODD_TRIP_LISTING
!
USE MODN_TRIP_RUN, ONLY : LRESTART, LPRINT,   &
                          XTSTEP_RUN, XTSTEP_DIAG
!
USE MODD_TRIP_PAR, ONLY : XUNDEF, NUNDEF, XDAY
!
USE MODI_READ_NAM_TRIP_RUN
USE MODI_READ_NAM_TRIP
!
USE MODI_ABORT_TRIP
USE MODI_GET_TRIP_GRID_CONF
!
USE MODI_INIT_TRIP
USE MODI_INIT_TRIP_PAR
USE MODI_TRIP_RUN_CONF
USE MODI_TRIP_RESTART
USE MODI_TRIP_DIAG_RUN
USE MODI_TRIP_RUN
!
USE MODI_TRIP_OASIS_INIT
USE MODI_TRIP_OASIS_READ_NAM
USE MODI_TRIP_OASIS_DEFINE
USE MODI_TRIP_OASIS_END
!
USE MODI_ALLOC_TRIP
USE MODI_DEALLOC_TRIP
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
INTEGER                           :: IYEAR               ! current year         (UTC)
INTEGER                           :: IMONTH              ! current month        (UTC)
INTEGER                           :: IDAY                ! current day          (UTC)
REAL                              :: ZTIME               ! current time           (s)
REAL                              :: ZRUNTIME            ! total simulated time   (s)
!
INTEGER                           :: INB_TSTEP_RUN       ! number of time step in the run
INTEGER                           :: ILON                ! Number of longitude
INTEGER                           :: ILAT                ! Number of latittude
!
INTEGER                           :: INB_OL              ! number of time step if forcing offline
INTEGER                           :: ILON_OL             ! Number of longitude if forcing offline
INTEGER                           :: ILAT_OL             ! Number of latittude if forcing offline
!
INTEGER                           :: IERR                ! Error value
INTEGER                           :: INPROC              ! Number of processes
INTEGER                           :: IRANK               ! Local process number
INTEGER                           :: ILOCAL_COMM         ! Local communicator
LOGICAL                           :: GOASIS              ! OASIS used(default=.false.)
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! --------------------------------------------------------------------------------------
! * 0. MPI and OASIS must be initialized before any DR_HOOK call
! --------------------------------------------------------------------------------------
!
CALL TRIP_OASIS_INIT(GOASIS,ILOCAL_COMM,ZRUNTIME)
!
! --------------------------------------------------------------------------------------
! * 1. Alloc trip variables and open listing
! --------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TRIP_MASTER',0,ZHOOK_HANDLE)
!
CALL ALLOC_TRIP(1)
!
CALL INIT_TRIP_PAR
!
OPEN(UNIT=NLISTING,FILE=CLISTING,FORM='FORMATTED',ACTION='WRITE')
!
! --------------------------------------------------------------------------------------
! * 2. Check run attributes
! --------------------------------------------------------------------------------------
!
!Inquire if trip is parallel or not: TRIP is only a monoprocess model for now
!
INPROC = NUNDEF
IRANK  = NUNDEF
!
#ifdef TRIPOASIS
IF(.NOT.GOASIS)THEN
  ILOCAL_COMM=MPI_COMM_WORLD
ENDIF
CALL MPI_COMM_SIZE(ILOCAL_COMM,INPROC,IERR)
CALL MPI_COMM_RANK(ILOCAL_COMM,IRANK,IERR)  
IF(INPROC==NUNDEF.OR.IRANK==NUNDEF)THEN
  WRITE(NLISTING,*)'TRIP_MASTER: PROBLEM WITH MPI, INPROC = ',INPROC
  WRITE(NLISTING,*)'TRIP_MASTER: PROBLEM WITH MPI, IRANK  = ',IRANK
  CALL ABORT_TRIP('TRIP_MASTER: PROBLEM WITH MPI')
ENDIF
#endif
!
IF(INPROC>1.AND.INPROC<NUNDEF)THEN
  WRITE(NLISTING,*)'TRIP_MASTER: TRIP NOT YET PARALLELIZED, NPROC SHOULD BE 1'
  CALL ABORT_TRIP('TRIP_MASTER: TRIP NOT YET PARALLELIZED')
ENDIF
!
IF(GOASIS)THEN
  WRITE(NLISTING,*) '!!!!!!!!!!!!!!!!!!!!!!!'
  WRITE(NLISTING,*) '    OASIS is used      ' 
  WRITE(NLISTING,*) '                       '
  WRITE(NLISTING,*) 'Number of processes   :', INPROC
  WRITE(NLISTING,*) 'Local process number  :', IRANK
  WRITE(NLISTING,*) 'Local communicator    :', ILOCAL_COMM
  WRITE(NLISTING,*) '!!!!!!!!!!!!!!!!!!!!!!!'
  WRITE(NLISTING,*) '                       '
ELSE
  WRITE(NLISTING,*) '!!!!!!!!!!!!!!!!!!!!!!!'
  WRITE(NLISTING,*) '   TRIP run offline    ' 
  WRITE(NLISTING,*) '!!!!!!!!!!!!!!!!!!!!!!!'
  WRITE(NLISTING,*) '                       '
ENDIF
!
! --------------------------------------------------------------------------------------
! * 3. read namelists
! --------------------------------------------------------------------------------------
!
CALL READ_NAM_TRIP_RUN(NLISTING)
!
CALL READ_NAM_TRIP(NLISTING)
!
IF(GOASIS)THEN
  CALL TRIP_OASIS_READ_NAM(NLISTING,ZRUNTIME)
ENDIF
!
! --------------------------------------------------------------------------------------
! * 4. TRIP initializations
! --------------------------------------------------------------------------------------
!
CALL GOTO_TRIP (1)
!
CALL READ_NAM_TRIP_GRID(TPG, &
                        NLISTING)
!
CALL INIT_TRIP(TPDG, TP, TPG, &
               IYEAR,IMONTH,IDAY,ZTIME,ILON,ILAT,XTSTEP_RUN,XTSTEP_DIAG,LRESTART)
!
! --------------------------------------------------------------------------------------
! * 5. TRIP - OASIS  grid, partitions and local field definitions
! --------------------------------------------------------------------------------------
!
IF(GOASIS)THEN
  CALL TRIP_OASIS_DEFINE(NLISTING,ILON,ILAT)
ENDIF
!
! --------------------------------------------------------------------------------------
! * 6. Get run configuration
! --------------------------------------------------------------------------------------
!
CALL TRIP_RUN_CONF(NLISTING,GOASIS,IYEAR,IMONTH,IDAY,ZTIME, &
                   ILON,ILAT,INB_TSTEP_RUN,ZRUNTIME         )
!
IF(GOASIS)THEN
  INB_OL  = 0
  ILON_OL = 0
  ILAT_OL = 0
ELSE
  INB_OL  = INB_TSTEP_RUN
  ILON_OL = ILON
  ILAT_OL = ILAT
ENDIF
!
! --------------------------------------------------------------------------------------
! * 7. Read and prepare drainage and runoff if offline
! --------------------------------------------------------------------------------------
!
CALL TRIP_RUN(TPDG, TP, TPG, &
              GOASIS,                           &
              NLISTING,ILON,ILAT,INB_TSTEP_RUN, &
              ZRUNTIME,ILON_OL,ILAT_OL,INB_OL,  &
              IYEAR,IMONTH,IDAY,ZTIME           )
!
!-------------------------------------------------------------------------------
! * 9. Store run mean diagnostic and write restart
!-------------------------------------------------------------------------------
!
CALL TRIP_DIAG_RUN(TPDG, TPG, &
                   NLISTING,ILON,ILAT,ZRUNTIME)
!
IF(LRESTART)THEN
   CALL TRIP_RESTART(TP, TPG, &
                     NLISTING,IYEAR,IMONTH,IDAY,ZTIME,ILON,ILAT)
ENDIF
!
! --------------------------------------------------------------------------------------
! * 10. End of run
! --------------------------------------------------------------------------------------
!
CLOSE(NLISTING)
!
WRITE(*,*) ' '
WRITE(*,*) '    ------------------------------'
WRITE(*,*) '    | TRIP MASTER ENDS CORRECTLY |'
WRITE(*,*) '    ------------------------------'
WRITE(*,*) ' '
!
CALL DEALLOC_TRIP
!
IF (LHOOK) CALL DR_HOOK('TRIP_MASTER',1,ZHOOK_HANDLE)
!
! --------------------------------------------------------------------------------------
! * 11. MPI and OASIS must be finalized after the last DR_HOOK call
! --------------------------------------------------------------------------------------
!
IF(GOASIS)THEN
  CALL TRIP_OASIS_END
ENDIF
!
!-------------------------------------------------------------------------------
END PROGRAM TRIP_MASTER
