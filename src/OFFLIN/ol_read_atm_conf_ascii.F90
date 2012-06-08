!     #########
SUBROUTINE OL_READ_ATM_CONF_ASCII (HSURF_FILETYPE, HFORCING_FILETYPE,  &
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
USE MODI_READ_SURF
USE MODI_GET_LUOUT
USE MODD_TYPE_DATE_SURF
USE MODD_ARCH, ONLY : LITTLE_ENDIAN_ARCH
USE MODD_IO_SURF_ASC, ONLY : NNI_FORC
!==================================================================
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
USE MODI_END_IO_SURF_n
USE MODI_INIT_IO_SURF_n
IMPLICIT NONE
!
CHARACTER(LEN=6), INTENT(IN)  :: HSURF_FILETYPE
CHARACTER(LEN=6), INTENT(IN)  :: HFORCING_FILETYPE
INTEGER,          INTENT(OUT) :: KNI
INTEGER,          INTENT(OUT) :: KYEAR, KMONTH, KDAY
REAL,             INTENT(OUT) :: PDURATION,PTSTEP_FORC
REAL,             INTENT(OUT) :: PTIME
REAL, DIMENSION(:),  POINTER  :: PLAT, PLON
REAL, DIMENSION(:),  POINTER  :: PZS 
REAL, DIMENSION(:),  POINTER  :: PZREF, PUREF
!
INTEGER                       :: IYEAR, IMONTH, IDAY
REAL                          :: ZTIME
INTEGER                       :: IRET, INB_FORC
CHARACTER(LEN=1)              :: YSWAP
INTEGER                       :: INI
INTEGER                       :: ILUOUT
TYPE (DATE_TIME)              :: TTIME
REAL                          :: ZWORK
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!==================================================================
!
IF (LHOOK) CALL DR_HOOK('OL_READ_ATM_CONF_ASCII',0,ZHOOK_HANDLE)
CALL GET_LUOUT(HSURF_FILETYPE,ILUOUT) 
!
!*      1.    Define configuration parameters
!
IF (HFORCING_FILETYPE == 'BINARY') READ(21,*) YSWAP
IF (YSWAP.EQ.'Y') THEN 
        LITTLE_ENDIAN_ARCH=.NOT.LITTLE_ENDIAN_ARCH
        WRITE(ILUOUT,*) '*******************************************************************'
        WRITE(ILUOUT,*) 'Architecture of the machine needs to swap LITTLE_ENDIAN_ARCH to ', &
                        LITTLE_ENDIAN_ARCH  
        WRITE(ILUOUT,*) '*******************************************************************'
ENDIF
!
READ(21,*) INI
READ(21,*) INB_FORC
READ(21,*) PTSTEP_FORC
!
NNI_FORC = INI
!
PDURATION = ( INB_FORC - 1 ) * PTSTEP_FORC
!
!*      2.    Read full grid dimension and date
!
CALL INIT_IO_SURF_n(HSURF_FILETYPE,'FULL  ','SURF  ','READ ')
!
CALL READ_SURF(HSURF_FILETYPE,'DIM_FULL',KNI,IRET)
CALL READ_SURF(HSURF_FILETYPE,'DTCUR',TTIME,IRET)
!
CALL END_IO_SURF_n(HSURF_FILETYPE)
!
READ(21,*) KYEAR
READ(21,*) KMONTH
READ(21,*) KDAY
READ(21,*) PTIME
!
!
!*      3.    Check the consistency
!
!
IF (KNI /= INI .AND. INI/=1) THEN
   WRITE(ILUOUT,*)' NUMBER OF GRID POINTS INCONSISTENCY: ',KNI,'/',INI
   CALL ABOR1_SFX('OL_READ_ATM_CONF_ASCII: NUMBER OF GRID POINTS INCONSISTENCY')
ENDIF
!
!*      4.    Geographical initialization
!
ALLOCATE(PLON(KNI))
ALLOCATE(PLAT(KNI))
ALLOCATE(PZS(KNI))
ALLOCATE(PZREF(KNI))
ALLOCATE(PUREF(KNI))
!
IF (INI==1) THEN
READ(UNIT=21,FMT='(F15.8)') ZWORK ; PLON(:) = ZWORK
READ(UNIT=21,FMT='(F15.8)') ZWORK ; PLAT(:) = ZWORK
READ(UNIT=21,FMT='(F15.8)') ZWORK ; PZS(:) = ZWORK
READ(UNIT=21,FMT='(F15.8)') ZWORK ; PZREF(:) = ZWORK
READ(UNIT=21,FMT='(F15.8)') ZWORK ; PUREF(:) = ZWORK
ELSE
READ(UNIT=21,FMT='(50(F15.8))') PLON(:)
READ(UNIT=21,FMT='(50(F15.8))') PLAT(:)
READ(UNIT=21,FMT='(50(F15.8))') PZS(:)
READ(UNIT=21,FMT='(50(F15.8))') PZREF(:)
READ(UNIT=21,FMT='(50(F15.8))') PUREF(:)
END IF
!
!

!
!  check date and time
!
IYEAR  = TTIME%TDATE%YEAR
IMONTH = TTIME%TDATE%MONTH
IDAY   = TTIME%TDATE%DAY
ZTIME  = TTIME%TIME
!
IF ( (KYEAR /= IYEAR) .OR. (KMONTH /= IMONTH) .OR. (KDAY /= IDAY) ) THEN
   WRITE(ILUOUT,*)' DATE INCONSISTANCY: ',KYEAR,KMONTH,KDAY,'/',IYEAR,IMONTH,IDAY
   CALL ABOR1_SFX('OL_READ_ATM_CONF_ASCII: DATE INCONSISTENCY')
ENDIF
!
IF ( PTIME /= ZTIME ) THEN
   WRITE(ILUOUT,*)' TIME INCONSISTANCY: ',PTIME,'/',ZTIME
   CALL ABOR1_SFX('OL_READ_ATM_CONF_ASCII: TIME INCONSISTENCY')
ENDIF
IF (LHOOK) CALL DR_HOOK('OL_READ_ATM_CONF_ASCII',1,ZHOOK_HANDLE)
!

!==================================================================
!
END SUBROUTINE OL_READ_ATM_CONF_ASCII
