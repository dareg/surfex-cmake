!     #########
SUBROUTINE OL_READ_ATM_CONF_NETCDF(HSURF_FILETYPE,                &
                                     PDURATION, PTSTEP_FORC, KNI, &
                                     KYEAR, KMONTH, KDAY, PTIME,  &
                                     PLAT, PLON, PZS,             &
                                     PZREF, PUREF                 )  
!
!==================================================================
!!****  *OL_READ_ATM_CONF* - Initialization routine
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
!!	F. Habets   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      Modified by P. Le Moigne (04/2005): cleaning and checking
!!      Modified by P. Le Moigne (04/2006): init_io_surf for nature
!!                  with GTMSK to read dimensions.
!==================================================================
!
USE MODD_TYPE_DATE_SURF
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO, NCOMM, NPROC, XTIME_COMM_READ, XTIME_NPIO_READ
!
USE MODI_SET_SURFEX_FILEIN
USE MODI_GET_LUOUT
USE MODI_INIT_IO_SURF_n
USE MODI_READ_SURF
USE MODI_END_IO_SURF_n
USE MODI_GET_SIZE_FULL_n
USE MODI_ABOR1_SFX
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
#ifndef NOMPI
INCLUDE "mpif.h"
#endif
!
CHARACTER(LEN=6), INTENT(IN)  :: HSURF_FILETYPE
INTEGER,          INTENT(OUT) :: KNI
INTEGER,          INTENT(OUT) :: KYEAR, KMONTH, KDAY
REAL,             INTENT(OUT) :: PDURATION,PTSTEP_FORC
REAL,             INTENT(OUT) :: PTIME
REAL, DIMENSION(:),  POINTER  :: PLAT, PLON
REAL, DIMENSION(:),  POINTER  :: PZS 
REAL, DIMENSION(:),  POINTER  :: PZREF, PUREF
!
CHARACTER(LEN=100)            :: YUNITS, YPOUB, YFMT 
INTEGER                       :: IYEAR, IMONTH, IDAY
REAL                          :: ZTIME, ZHOUR, ZMIN, ZSEC
INTEGER                       :: IRET, INB_FORC
INTEGER                       :: INI, IDIM_FULL
INTEGER                       :: ILUOUT
TYPE (DATE_TIME)              :: TTIME
INTEGER                       :: INFOMPI
DOUBLE PRECISION   :: XTIME0
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!==================================================================
!
!*      0.    IO initialization
!
IF (LHOOK) CALL DR_HOOK('OL_READ_ATM_CONF_NETCDF',0,ZHOOK_HANDLE)
!
IF (NRANK==NPIO) THEN
  !
#ifndef NOMPI  
  XTIME0 = MPI_WTIME()
#endif
  !  
  CALL GET_LUOUT(HSURF_FILETYPE,ILUOUT)
  !
  !*      1.    Read parameters from netcdf forcing file
  !
!$OMP SINGLE
  !  
  YUNITS = ""
  CALL READ_SURF_DIM_OL(YUNITS, INB_FORC, INI, IRET)
  !
!$OMP END SINGLE COPYPRIVATE(YUNITS,INB_FORC,INI,IRET)
  !  
#ifndef NOMPI  
  XTIME_NPIO_READ = XTIME_NPIO_READ + (MPI_WTIME() - XTIME0)
#endif
  !  
ENDIF
!
CALL READ_SURF('OFFLIN','FORC_TIME_STEP'  ,PTSTEP_FORC   ,IRET)
!
IF (NPROC>1) THEN
#ifndef NOMPI        
  XTIME0 = MPI_WTIME()
!$OMP SINGLE  
  CALL MPI_BCAST(INB_FORC,KIND(INB_FORC)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)
!$OMP END SINGLE COPYPRIVATE(INB_FORC)  
  XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
#endif
ENDIF
!
PDURATION = ( INB_FORC - 1 ) * PTSTEP_FORC
!
!*      2.    Read full grid dimension and date
!
CALL SET_SURFEX_FILEIN(HSURF_FILETYPE,'PREP')
CALL INIT_IO_SURF_n(HSURF_FILETYPE,'FULL  ','SURF  ','READ ')
!
CALL READ_SURF(HSURF_FILETYPE,'DIM_FULL',IDIM_FULL,IRET)
CALL READ_SURF(HSURF_FILETYPE,'DTCUR',TTIME,IRET)
!
KYEAR  = TTIME%TDATE%YEAR
KMONTH = TTIME%TDATE%MONTH
KDAY   = TTIME%TDATE%DAY
PTIME  = TTIME%TIME
!
CALL END_IO_SURF_n(HSURF_FILETYPE)
!
!*      5.    Geographical initialization
!
CALL GET_SIZE_FULL_n('OFFLIN ',IDIM_FULL,KNI) 
!
ALLOCATE(PLON(KNI))
ALLOCATE(PLAT(KNI))
ALLOCATE(PZS(KNI))
ALLOCATE(PZREF(KNI))
ALLOCATE(PUREF(KNI))
!
CALL READ_SURF('OFFLIN','LAT',PLAT,IRET)
CALL READ_SURF('OFFLIN','LON',PLON,IRET)
CALL READ_SURF('OFFLIN','ZS',PZS,IRET)
CALL READ_SURF('OFFLIN','ZREF',PZREF,IRET)
CALL READ_SURF('OFFLIN','UREF',PUREF,IRET)
!
!*      6.    Check the consistency
!
IF (NRANK == NPIO) THEN
  !
  IF (IDIM_FULL /= INI) THEN
    WRITE(ILUOUT,*)' NUMBER OF GRID POINTS INCONSISTENCY: ',KNI,'/',INI
    CALL ABOR1_SFX('OL_READ_ATM_CONF_NETCDF: NUMBER OF GRID POINTS INCONSISTENCY')
  ENDIF
  !
!$OMP SINGLE
  !
  IF (YUNITS(1:1)=='d') THEN
    YFMT="(A11"
  ELSEIF (YUNITS(1:1)=='h') THEN
    YFMT="(A12"
  ELSEIF (YUNITS(1:1)=='s' .OR. YUNITS(1:1)=='m') THEN
    YFMT="(A14"
  ENDIF
  !
  YFMT = TRIM(YFMT)//",I4,2(A1,I2),3(A1,F2.0))"
  !
  READ(YUNITS,FMT=YFMT) YPOUB,IYEAR,YPOUB,IMONTH,YPOUB,IDAY,&
          YPOUB,ZHOUR,YPOUB,ZMIN,YPOUB,ZSEC
  !
  ZTIME = ZHOUR * 3600. + ZMIN * 60. + ZSEC
  !
  IF ( (KYEAR /= IYEAR) .OR. (KMONTH /= IMONTH) .OR. (KDAY /= IDAY) ) THEN
    WRITE(ILUOUT,*)' DATE INCONSISTANCY: ',KYEAR,KMONTH,KDAY,'/',IYEAR,IMONTH,IDAY
    CALL ABOR1_SFX('OL_READ_ATM_CONF_NETCDF: DATE INCONSISTENCY')
  ENDIF
  !
  IF ( PTIME /= ZTIME ) THEN
    WRITE(ILUOUT,*)' TIME INCONSISTANCY: ',PTIME,'/',ZTIME
    CALL ABOR1_SFX('OL_READ_ATM_CONF_NETCDF: TIME INCONSISTENCY')
  ENDIF
  !
!$OMP END SINGLE
  !
ENDIF
!
IF (LHOOK) CALL DR_HOOK('OL_READ_ATM_CONF_NETCDF',1,ZHOOK_HANDLE)
!==================================================================
CONTAINS
!
!     #############################################################
      SUBROUTINE READ_SURF_DIM_OL(HUNITS,KSIZE,KNI,KRESP)
!     #############################################################
!
USE MODI_OL_FIND_FILE_READ
!
IMPLICIT NONE
INCLUDE "netcdf.inc"
!
!*      0.1   Declarations of arguments
!
CHARACTER(LEN=100), INTENT(OUT) :: HUNITS
INTEGER,            INTENT(OUT) :: KSIZE
INTEGER,            INTENT(OUT) :: KNI
INTEGER,            INTENT(OUT) :: KRESP    
!
!*      0.2   Declarations of local variables
!
INTEGER :: IFILE_ID,IVAR_ID,INDIMS,JRET,JDIM
INTEGER,DIMENSION(2) :: IDIMIDS,IDIMLEN
INTEGER,DIMENSION(5) :: IRET
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('OL_READ_ATM_CONF_NETCDF:READ_SURF_DIM_OL',0,ZHOOK_HANDLE)
KRESP=0

! 0. find filename
! -----------------
CALL OL_FIND_FILE_READ('time',IFILE_ID)
 
IF (IFILE_ID.NE.0) THEN
    
  ! 1. Find id of the variable
  !----------------------------
  IRET(1)=NF_INQ_VARID   (IFILE_ID,'time',IVAR_ID)
  IRET(2)=NF_INQ_VARNDIMS(IFILE_ID,IVAR_ID,INDIMS)
  IRET(3)=NF_INQ_VARDIMID(IFILE_ID,IVAR_ID,IDIMIDS(1:1))
  IRET(4)=NF_INQ_DIMLEN(IFILE_ID,IDIMIDS(1),KSIZE)
  IRET(5)=NF_GET_ATT_TEXT(IFILE_ID,IVAR_ID,'units',HUNITS)

  
  ! 3. Check for errors
  !--------------------
  DO JRET=1,5
    IF (IFILE_ID==0.OR.IRET(JRET).NE.NF_NOERR) KRESP=1
  ENDDO

  IRET(1)=NF_INQ_VARID   (IFILE_ID,'LON',IVAR_ID)
  IRET(2)=NF_INQ_VARNDIMS(IFILE_ID,IVAR_ID,INDIMS)
  IRET(3)=NF_INQ_VARDIMID(IFILE_ID,IVAR_ID,IDIMIDS(:))
  IDIMLEN(:)=1.
  DO JDIM=1,INDIMS
    IRET(4)=NF_INQ_DIMLEN(IFILE_ID,IDIMIDS(JDIM),IDIMLEN(JDIM))
  ENDDO
  KNI = IDIMLEN(1) * IDIMLEN(2)

  DO JRET=1,4
    IF (IFILE_ID==0.OR.IRET(JRET).NE.NF_NOERR) KRESP=1
  ENDDO

ENDIF
!
IF (LHOOK) CALL DR_HOOK('OL_READ_ATM_CONF_NETCDF:READ_SURF_DIM_OL',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURF_DIM_OL
!
END SUBROUTINE OL_READ_ATM_CONF_NETCDF
