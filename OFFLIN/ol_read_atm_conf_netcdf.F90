!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
SUBROUTINE OL_READ_ATM_CONF_NETCDF (DTCO, U, HGRID, HSURF_FILETYPE, ODELAYEDSTART_NC, KDATESTOP,  &
                                    PDURATION, PTSTEP_FORC, KNI, KYEAR, KMONTH, KDAY, PTIME, &
                                    PLAT, PLON, PZS, PZREF, PUREF,KTIMESTARTINDEX         )  
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
!!      F. Habets   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      Modified by P. Le Moigne (04/2005): cleaning and checking
!!      Modified by P. Le Moigne (04/2006): init_io_surf for nature
!!                  with GTMSK to read dimensions.
!!      Modified by Matthieu Lafaysse 2012-11-12
!==================================================================
!
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_TYPE_DATE_SURF, ONLY: DATE_TIME
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO, NCOMM, NPROC, XTIME_COMM_READ, XTIME_NPIO_READ, LSFX_MPI
!
USE MODI_SET_SURFEX_FILEIN
USE MODI_GET_LUOUT
USE MODI_INIT_IO_SURF_n
USE MODI_READ_SURF
USE MODI_END_IO_SURF_n
USE MODI_GET_SIZE_FULL_n
USE MODI_ABOR1_SFX
!
USE MODE_DATES_NETCDF, ONLY : NETCDF2DATE
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
USE NETCDF
!
IMPLICIT NONE
!
#ifdef SFX_MPI
INCLUDE "mpif.h"
#endif
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
 CHARACTER(LEN=*), INTENT(IN)  :: HGRID
 CHARACTER(LEN=6), INTENT(IN)  :: HSURF_FILETYPE
LOGICAL, INTENT(IN)            :: ODELAYEDSTART_NC ! Allow the simulation to start from a different time step than the first record of a netcdf file
INTEGER,DIMENSION(4),INTENT(IN) :: KDATESTOP !Allow the simulation to end at a different time step than the last record of a netcdf file 
INTEGER,          INTENT(OUT) :: KNI
INTEGER,          INTENT(OUT) :: KYEAR, KMONTH, KDAY
REAL,             INTENT(OUT) :: PDURATION,PTSTEP_FORC
REAL,             INTENT(OUT) :: PTIME
REAL, DIMENSION(:),  POINTER  :: PLAT, PLON
REAL, DIMENSION(:),  POINTER  :: PZS 
REAL, DIMENSION(:),  POINTER  :: PZREF, PUREF
INTEGER,          INTENT(OUT) :: KTIMESTARTINDEX ! index from which we start reading FORCING.nc
!
 CHARACTER(LEN=100)           :: YUNITS
REAL, DIMENSION(:), POINTER   :: ZTIMEFILE
REAL,DIMENSION(:),ALLOCATABLE :: ZLAT1D, ZLON1D
REAL                          :: ZTIME
REAL                          :: ZFIRSTTIMEFILE
INTEGER                       :: IYEAR, IMONTH, IDAY
INTEGER                       :: IRET, INB_FORC
INTEGER                       :: INI, IDIM_FULL
INTEGER                       :: ILUOUT
INTEGER                       :: JINDEX
INTEGER                       :: INFOMPI
TYPE (DATE_TIME)              :: TTIME
TYPE (DATE_TIME), DIMENSION(1) :: TZ_FIRSTDATEFILE
TYPE (DATE_TIME), DIMENSION(:), ALLOCATABLE :: TZ_DATEFILE
INTEGER :: JPOINT, INLAT1D, INLON1D
LOGICAL :: GLONLAT1D
DOUBLE PRECISION   :: XTIME0
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!==================================================================
!
!*      0.    IO initialization
!
IF (LHOOK) CALL DR_HOOK('OL_READ_ATM_CONF_NETCDF',0,ZHOOK_HANDLE)
!
IF (NRANK==NPIO) THEN
  !
#ifdef SFX_MPI
  IF (LSFX_MPI) XTIME0 = MPI_WTIME()
#endif
  !
  CALL GET_LUOUT(HSURF_FILETYPE,ILUOUT)
  !
  !*      1.    Read parameters from netcdf forcing file
  !
  YUNITS = ""
  CALL READ_SURF_DIM_OL(YUNITS, INB_FORC, INI, ZFIRSTTIMEFILE, IRET,&
                        GLONLAT1D, INLON1D, INLAT1D, ZTIMEFILE)
  !
#ifdef SFX_MPI
  IF (LSFX_MPI) XTIME_NPIO_READ = XTIME_NPIO_READ + (MPI_WTIME() - XTIME0)
#endif
  !
ENDIF
!
IF (NPROC>1) THEN
#ifdef SFX_MPI
  IF (LSFX_MPI) THEN
    XTIME0 = MPI_WTIME()
    CALL MPI_BCAST(INB_FORC,KIND(INB_FORC)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)
    CALL MPI_BCAST(GLONLAT1D,1,MPI_LOGICAL,NPIO,NCOMM,INFOMPI)
    CALL MPI_BCAST(INLON1D,KIND(INLON1D)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)
    CALL MPI_BCAST(INLAT1D,KIND(INLAT1D)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)    
    XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
  ENDIF
#endif
ENDIF
!
 CALL READ_SURF(&
                'OFFLIN','FRC_TIME_STP'  ,PTSTEP_FORC   ,IRET)
!
PDURATION = ( INB_FORC - 1 ) * PTSTEP_FORC
!
!*      2.    Read full grid dimension and date
!
 CALL SET_SURFEX_FILEIN(HSURF_FILETYPE,'PREP')
CALL INIT_IO_SURF_n(DTCO, U, HSURF_FILETYPE,'FULL  ','SURF  ','READ ')  
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
 CALL GET_SIZE_FULL_n('OFFLIN ',IDIM_FULL,U%NSIZE_FULL,KNI) 
!
ALLOCATE(PLON(KNI))
ALLOCATE(PLAT(KNI))
ALLOCATE(PZS (KNI))
ALLOCATE(PZREF(KNI))
ALLOCATE(PUREF(KNI))
!

IF (GLONLAT1D) THEN
  ALLOCATE(ZLAT1D(INLAT1D))
  ALLOCATE(ZLON1D(INLON1D))
  ! we have to read PLAT and PLON due to parallel IO.
  ! however lat / lon have just INLAT1D and INLON1D dimension lengths in the netcdf file
  CALL READ_SURF('OFFLIN','LAT',PLAT,IRET)
  CALL READ_SURF('OFFLIN','LON',PLON,IRET)
  ZLAT1D=PLAT(1:INLAT1D)
  ZLON1D=PLON(1:INLON1D)
ELSE
  CALL READ_SURF('OFFLIN','LAT',PLAT,IRET)
  CALL READ_SURF('OFFLIN','LON',PLON,IRET)
ENDIF
!
 CALL READ_SURF('OFFLIN','ZS',PZS,IRET)
 CALL READ_SURF('OFFLIN','ZREF',PZREF,IRET)
 CALL READ_SURF('OFFLIN','UREF',PUREF,IRET)
!
! Modif Matthieu Lafaysse
! If LON LAT REGULAR grid and 1 dimension LON LAT variables
! Affect 1d lon lat to total tables
IF (GLONLAT1D) THEN
  !
  IF ( INLAT1D*INLON1D==IDIM_FULL ) THEN
    DO JPOINT = 1,IDIM_FULL
      PLAT(JPOINT) = ZLAT1D((JPOINT-1)/INLON1D+1)
      PLON(JPOINT) = ZLON1D(MOD(JPOINT-1,INLON1D)+1)
    END DO
    DEALLOCATE(ZLAT1D)
    DEALLOCATE(ZLON1D)
  ELSE
    WRITE(ILUOUT,*)' NUMBER OF GRID POINTS INCONSISTENCY: ',IDIM_FULL,'/',INLAT1D*INLON1D
    CALL ABOR1_SFX('OL_READ_ATM_CONF_NETCDF: NUMBER OF GRID POINTS INCONSISTENCY')
  END IF
END IF
!
!*      6.    Check the consistency
!
IF (NRANK == NPIO) THEN
  !
  IF (IDIM_FULL /= INI) THEN
    WRITE(ILUOUT,*)' NUMBER OF GRID POINTS INCONSISTENCY: ',IDIM_FULL,'/',INI
    CALL ABOR1_SFX('OL_READ_ATM_CONF_NETCDF: NUMBER OF GRID POINTS INCONSISTENCY')
  ENDIF
  !
  IF ( ODELAYEDSTART_NC .OR. KDATESTOP(1)/=0 ) THEN
    !
    ! Convert time variable in dates
    ALLOCATE(TZ_DATEFILE(INB_FORC))
    CALL NETCDF2DATE(ZTIMEFILE, YUNITS, TZ_DATEFILE)
    !
  END IF
  DEALLOCATE(ZTIMEFILE)
  !
  IF ( ODELAYEDSTART_NC ) THEN
    !
    KTIMESTARTINDEX = -1
    !
    ! Look for a date equal to the date prescribed in OPTIONS.NAM
    DO JINDEX = 1,INB_FORC
      IYEAR  = TZ_DATEFILE(JINDEX)%TDATE%YEAR
      IMONTH = TZ_DATEFILE(JINDEX)%TDATE%MONTH
      IDAY   = TZ_DATEFILE(JINDEX)%TDATE%DAY
      ZTIME  = TZ_DATEFILE(JINDEX)%TIME* 3600. 
      
      IF ( KYEAR==IYEAR .AND. KMONTH==IMONTH .AND. KDAY==IDAY .AND. PTIME==ZTIME ) THEN
        KTIMESTARTINDEX = JINDEX
        INB_FORC        = INB_FORC-KTIMESTARTINDEX+1
        PDURATION       = ( INB_FORC - 1 ) * PTSTEP_FORC
        EXIT
      END IF
      
    END DO
    
    ! Error if the initial date is not found
    IF ( KTIMESTARTINDEX==-1 ) THEN
      WRITE(ILUOUT,*)'IN THE FORCING FILE, WE CAN NOT FIND THIS INITIAL DATE :',KYEAR,KMONTH,KDAY,PTIME
      CALL ABOR1_SFX('OL_READ_ATM_CONF_NETCDF: DATE INCONSISTENCY')
    END IF
    
  ELSE

    KTIMESTARTINDEX = 1

    CALL NETCDF2DATE((/ZFIRSTTIMEFILE/),YUNITS,TZ_FIRSTDATEFILE)
    IYEAR  = TZ_FIRSTDATEFILE(1)%TDATE%YEAR
    IMONTH = TZ_FIRSTDATEFILE(1)%TDATE%MONTH
    IDAY   = TZ_FIRSTDATEFILE(1)%TDATE%DAY
    ZTIME  = TZ_FIRSTDATEFILE(1)%TIME* 3600.
    !
    IF ( (KYEAR /= IYEAR) .OR. (KMONTH /= IMONTH) .OR. (KDAY /= IDAY) ) THEN
      WRITE(ILUOUT,*)' DATE INCONSISTENCY: ',KYEAR,KMONTH,KDAY,'/',IYEAR,IMONTH,IDAY
      CALL ABOR1_SFX('OL_READ_ATM_CONF_NETCDF: DATE INCONSISTENCY')
    ENDIF
    !
    IF ( PTIME /= ZTIME ) THEN
      WRITE(ILUOUT,*)' TIME INCONSISTENCY: ',PTIME,'/',ZTIME
      CALL ABOR1_SFX('OL_READ_ATM_CONF_NETCDF: TIME INCONSISTENCY')
    ENDIF
    !
  ENDIF
  !
  ! Look for the prescribed final date to modify PDURATION
  ! If we don't find it, ignore.
  IF ( KDATESTOP(1)/=0 ) THEN
    !
    DO JINDEX = KTIMESTARTINDEX,INB_FORC+KTIMESTARTINDEX-1
      IYEAR  = TZ_DATEFILE(JINDEX)%TDATE%YEAR
      IMONTH = TZ_DATEFILE(JINDEX)%TDATE%MONTH
      IDAY   = TZ_DATEFILE(JINDEX)%TDATE%DAY
      ZTIME  = TZ_DATEFILE(JINDEX)%TIME* 3600.
      IF ( KDATESTOP(1)==IYEAR .AND. KDATESTOP(2)==IMONTH .AND. KDATESTOP(3)==IDAY &
                .AND. KDATESTOP(4)==INT(ZTIME) ) THEN
        INB_FORC  = JINDEX-KTIMESTARTINDEX+1
        PDURATION = ( INB_FORC - 1 ) * PTSTEP_FORC
        EXIT
      END IF
    END DO
    !
  END IF
  !
  IF ( ODELAYEDSTART_NC .OR. KDATESTOP(1)/=0 ) THEN
    DEALLOCATE(TZ_DATEFILE)
  END IF  
  !
ENDIF
!
IF (NPROC>1) THEN
#ifdef SFX_MPI       
  IF (LSFX_MPI) THEN
    XTIME0 = MPI_WTIME()
    CALL MPI_BCAST(KTIMESTARTINDEX,KIND(KTIMESTARTINDEX)/4,MPI_INTEGER,NPIO,NCOMM,INFOMPI)  
    CALL MPI_BCAST(PDURATION,KIND(PDURATION)/4,MPI_REAL,NPIO,NCOMM,INFOMPI)
    XTIME_COMM_READ = XTIME_COMM_READ + (MPI_WTIME() - XTIME0)
  ENDIF
#endif
ENDIF
!
IF (LHOOK) CALL DR_HOOK('OL_READ_ATM_CONF_NETCDF',1,ZHOOK_HANDLE)
!==================================================================
CONTAINS
!
!     #############################################################
      SUBROUTINE READ_SURF_DIM_OL(HUNITS,KSIZE,KNI,PFIRSTTIMEFILE,KRESP,&
                                  OLONLAT1D,KNLON1D,KNLAT1D,PTIMEFILE)
!     #############################################################
!
USE MODI_OL_FIND_FILE_READ
!
USE NETCDF
!
IMPLICIT NONE
!
!*      0.1   Declarations of arguments
!
 CHARACTER(LEN=100), INTENT(OUT) :: HUNITS
INTEGER,            INTENT(OUT) :: KSIZE
INTEGER,            INTENT(OUT) :: KNI
REAL,               INTENT(OUT) :: PFIRSTTIMEFILE
INTEGER,            INTENT(OUT) :: KRESP    
LOGICAL,            INTENT(OUT) :: OLONLAT1D
INTEGER,            INTENT(OUT) :: KNLON1D,KNLAT1D
REAL,DIMENSION(:), POINTER, INTENT(OUT) :: PTIMEFILE
!
!*      0.2   Declarations of local variables
!
INTEGER :: IFILE_ID,IVAR_ID,INDIMS,JRET,JDIM,ITYPE
INTEGER,DIMENSION(2) :: IDIMIDS,IDIMLEN
CHARACTER(20),DIMENSION(2) :: YDIMNAMES
INTEGER,DIMENSION(8) :: IRET

REAL*4, DIMENSION(:), ALLOCATABLE :: ZTIMEFILE4
REAL, DIMENSION(:), ALLOCATABLE   :: ZTIMEFILE
INTEGER, DIMENSION(:), ALLOCATABLE :: ITIMEFILE
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('OL_READ_ATM_CONF_NETCDF:READ_SURF_DIM_OL',0,ZHOOK_HANDLE)
KRESP=0

! 0. find filename
! -----------------
 CALL OL_FIND_FILE_READ('time',IFILE_ID)
 
IF (IFILE_ID.NE.0) THEN
    
  ! 1. Find id of the variable
  !----------------------------
  IRET(1)=NF90_INQ_VARID   (IFILE_ID,'time',IVAR_ID)
  IRET(2)=NF90_INQUIRE_VARIABLE(IFILE_ID,IVAR_ID,NDIMS=INDIMS)
  IRET(3)=NF90_INQUIRE_VARIABLE(IFILE_ID,IVAR_ID,DIMIDS=IDIMIDS(1:1))
  IRET(4)=NF90_INQUIRE_DIMENSION(IFILE_ID,IDIMIDS(1),LEN=KSIZE)
  IRET(5)=NF90_GET_ATT(IFILE_ID,IVAR_ID,'units',HUNITS)

  IRET(6)=NF90_INQUIRE_VARIABLE(IFILE_ID,IVAR_ID,XTYPE=ITYPE)

  ALLOCATE(PTIMEFILE(KSIZE))

  SELECT CASE (ITYPE)
    CASE (NF90_DOUBLE)
      ALLOCATE(ZTIMEFILE(KSIZE))
      IRET(7)=NF90_GET_VAR(IFILE_ID,IVAR_ID,ZTIMEFILE)
      PFIRSTTIMEFILE=ZTIMEFILE(1)
      PTIMEFILE(:)=ZTIMEFILE(:)
      DEALLOCATE(ZTIMEFILE)
    CASE (NF90_FLOAT)
      ALLOCATE(ZTIMEFILE4(KSIZE))
      IRET(7)=REAL(NF90_GET_VAR(IFILE_ID,IVAR_ID,ZTIMEFILE4))
      PFIRSTTIMEFILE=ZTIMEFILE4(1)
      PTIMEFILE(:)=ZTIMEFILE4(:)
      DEALLOCATE(ZTIMEFILE4)
    CASE (NF90_INT)
      ALLOCATE(ITIMEFILE(KSIZE))
      IRET(7)=REAL(NF90_GET_VAR(IFILE_ID,IVAR_ID,ITIMEFILE))
      PFIRSTTIMEFILE=ITIMEFILE(1)
      PTIMEFILE(:)=ITIMEFILE(:)
      DEALLOCATE(ITIMEFILE)
    CASE DEFAULT
      CALL ABOR1_SFX('OL_READ_ATM_CONF_NETCDF: TYPE OF TIME VARIABLE NOT KNOWN')
  END SELECT
  
  ! 3. Check for errors
  !--------------------
  DO JRET=1,7
    IF (IFILE_ID==0.OR.IRET(JRET).NE.NF90_NOERR) KRESP=1
  ENDDO

ENDIF

 CALL OL_FIND_FILE_READ('LON',IFILE_ID)
!
OLONLAT1D = .FALSE.
KNLAT1D = 0
KNLON1D = 0
!
IF (IFILE_ID.NE.0) THEN

  IRET(1)=NF90_INQ_VARID(IFILE_ID,'LON',IVAR_ID)
  IRET(2)=NF90_INQUIRE_VARIABLE(IFILE_ID,IVAR_ID,NDIMS=INDIMS)
  IRET(3)=NF90_INQUIRE_VARIABLE(IFILE_ID,IVAR_ID,DIMIDS=IDIMIDS(:))
  IDIMLEN(:)=1.
  DO JDIM=1,INDIMS
    IRET(4)=NF90_INQUIRE_DIMENSION(IFILE_ID,IDIMIDS(JDIM),NAME=YDIMNAMES(JDIM),LEN=IDIMLEN(JDIM))
  ENDDO
!
  IF ( HGRID.EQ.'LONLAT REG' .AND. INDIMS==1 ) THEN
    IF (TRIM(YDIMNAMES(1))/="Number_of_points") THEN
      ! Modif Matthieu Lafaysse
      ! If LON LAT REGULAR grid and 1 dimension LON LAT variables
      ! total dimension is nlon*nlat
    
      OLONLAT1D = .TRUE.
      CALL OL_FIND_FILE_READ('LAT',IFILE_ID)
      IRET(5) = NF90_INQ_VARID(IFILE_ID,'LAT',IVAR_ID)
      IRET(6) = NF90_INQUIRE_VARIABLE(IFILE_ID,IVAR_ID,NDIMS=INDIMS)
      IF ( INDIMS==1 ) THEN
        IRET(7) = NF90_INQUIRE_VARIABLE (IFILE_ID,IVAR_ID,DIMIDS=IDIMIDS(2:2))
        IRET(8) = NF90_INQUIRE_DIMENSION(IFILE_ID,IDIMIDS(2),LEN=IDIMLEN(2))
        KNLON1D = IDIMLEN(1)
        KNLAT1D = IDIMLEN(2)
      ELSE
        KRESP = 1
      END IF
    END IF
  END IF
  KNI = IDIMLEN(1) * IDIMLEN(2)

  DO JRET=1,4
    IF (IFILE_ID==0.OR.IRET(JRET).NE.NF90_NOERR) KRESP=1
  ENDDO

ENDIF
!
IF (LHOOK) CALL DR_HOOK('OL_READ_ATM_CONF_NETCDF:READ_SURF_DIM_OL',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
END SUBROUTINE READ_SURF_DIM_OL
!
END SUBROUTINE OL_READ_ATM_CONF_NETCDF
