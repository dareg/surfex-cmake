!------------------------
PROGRAM PREP
!------------------------
!!
!!    PURPOSE
!!    -------
!!   This program prepares the initial file for offline run
!!
!!    METHOD
!!    ------
!!   
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    P. LeMoigne                  Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original     22/04/04
!!    B. Decharme      2008 : FA type and TRIP prep coupling
!!
!----------------------------------------------------------------------------
!
USE MODE_POS_SURF
!
USE MODD_SURFEX_OMP, ONLY : NWORK, NWORK2, XWORK, XWORK2, XWORK3, NBLOCKTOT, &
                             NWORK_FULL, NWORK2_FULL, XWORK_FULL, XWORK2_FULL
!
USE MODN_IO_OFFLINE, ONLY : LLAND_USE
USE MODD_IO_SURF_ASC
USE MODD_IO_SURF_FA
USE MODD_IO_SURF_LFI
USE MODD_IO_SURF_NC
USE MODD_SURF_PAR
USE MODD_SURF_CONF, ONLY : CSOFTWARE
!
USE MODD_SURF_ATM, ONLY : LCPL_ESM
!
USE MODD_SURF_ATM_n, ONLY : TTIME
!
USE MODI_OPEN_NAMELIST
USE MODI_CLOSE_NAMELIST
USE MODI_READ_ALL_NAMELISTS
USE MODI_ALLOC_SURFEX
USE MODI_DEALLOC_SURFEX
USE MODI_GET_LUOUT
!
USE MODI_GOTO_SURFEX
USE MODI_GOTO_TRIP
USE MODI_INIT_PGD_SURF_ATM
USE MODI_IO_BUFF_CLEAN_n
USE MODI_PREP_SURF_ATM
USE MODI_PREP_SURF_TRIP
USE MODI_WRITE_DIAG_SURF_ATM_n
USE MODI_WRITE_HEADER_MNH
USE MODI_WRITE_SURF_ATM_n
!
#ifdef OFF
USE MODI_FANDAR
#endif
!
USE MODI_GET_LONLAT_n
USE MODI_FLAG_UPDATE
!
USE MODN_IO_OFFLINE
!------------------------------------------------------------------------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
#ifndef AIX64
INCLUDE 'omp_lib.h'
#endif
!
!*    0.     Declaration of local variables
!            ------------------------------
!
INTEGER            :: ILUOUT
INTEGER            :: ILUNAM
INTEGER            :: IYEAR, IMONTH, IDAY
REAL               :: ZTIME
LOGICAL            :: GFOUND

REAL, DIMENSION(0) :: ZZS
 CHARACTER(LEN=28)  :: YATMFILE  ='                            '  ! name of the Atmospheric file
 CHARACTER(LEN=6)   :: YATMFILETYPE ='      '                     ! type of the Atmospheric file
 CHARACTER(LEN=28)  :: YPGDFILE  ='                            '  ! name of the pgd file
 CHARACTER(LEN=6)   :: YPGDFILETYPE ='      '                     ! type of the pgd file
 CHARACTER(LEN=28)  :: YLUOUT    ='LISTING_PREP                '  ! name of listing
!
INTEGER, DIMENSION(11)  :: IDATEF
!
INTEGER :: IRET, INB
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!------------------------------------------------------------------------------
!
!
!*    1.      Set default names and parallelized I/O
!             --------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PREP',0,ZHOOK_HANDLE)
 CALL ALLOC_SURFEX(1)
CSOFTWARE='PREP'
!
!     1.1     initializations
!             ---------------
!
IYEAR    = NUNDEF
IMONTH   = NUNDEF
IDAY     = NUNDEF
ZTIME    = XUNDEF
!
LCPL_ESM = .FALSE.
LPREP    = .TRUE.
!
!     1.2     output listing
!             --------------
CLUOUT_LFI =  ADJUSTL(ADJUSTR(YLUOUT)//'.txt')
 CALL GET_LUOUT('ASCII ',ILUOUT)
OPEN (UNIT=ILUOUT,FILE=ADJUSTL(ADJUSTR(YLUOUT)//'.txt'),FORM='FORMATTED',ACTION='WRITE')
!
!     1.3     output file name read in namelist
!             ---------------------------------
 CALL OPEN_NAMELIST('ASCII ',ILUNAM,CNAMELIST)
!
 CALL POSNAM(ILUNAM,'NAM_IO_OFFLINE',GFOUND)
IF (GFOUND) READ (UNIT=ILUNAM,NML=NAM_IO_OFFLINE)
!
CFILEPGD     = ADJUSTL(ADJUSTR(CPGDFILE)//'.txt')
CFILEPGD_FA  = ADJUSTL(ADJUSTR(CPGDFILE)//'.fa')
CFILEPGD_LFI = CPGDFILE
CFILEPGD_NC  = ADJUSTL(ADJUSTR(CPGDFILE)//'.nc')
!
CFILEIN     = ADJUSTL(ADJUSTR(CPGDFILE)//'.txt')      ! output of PGD program
CFILEIN_FA  = ADJUSTL(ADJUSTR(CPGDFILE)//'.fa')
CFILEIN_LFI = CPGDFILE
CFILEIN_NC  = ADJUSTL(ADJUSTR(CPGDFILE)//'.nc')
!
CFILEIN_SAVE     = CFILEIN
CFILEIN_LFI_SAVE = CFILEIN_LFI
CFILEIN_FA_SAVE  = CFILEIN_FA
CFILEIN_NC_SAVE  = CFILEIN_NC
!
CFILEOUT    = ADJUSTL(ADJUSTR(CPREPFILE)//'.txt')
CFILEOUT_FA = ADJUSTL(ADJUSTR(CPREPFILE)//'.fa')
CFILEOUT_LFI= CPREPFILE
CFILEOUT_NC = ADJUSTL(ADJUSTR(CPREPFILE)//'.nc')
!
 CALL CLOSE_NAMELIST('ASCII ',ILUNAM)
!
 CALL READ_ALL_NAMELISTS(CSURF_FILETYPE,'PRE',.FALSE.)
!
 CALL GOTO_SURFEX(1,.TRUE.)
 CALL GOTO_TRIP(1,.TRUE.)
!
!*    2.      Preparation of surface physiographic fields
!             -------------------------------------------
!
!$OMP PARALLEL
!$ NBLOCKTOT = OMP_GET_NUM_THREADS()
!$OMP END PARALLEL
!
 CALL IO_BUFF_CLEAN_n
 CALL INIT_PGD_SURF_ATM(CSURF_FILETYPE,'PRE',YATMFILE,YATMFILETYPE, &
                         IYEAR, IMONTH, IDAY, ZTIME            ) 
!
 CALL IO_BUFF_CLEAN_n
 CALL PREP_SURF_ATM(CSURF_FILETYPE,YATMFILE,YATMFILETYPE,YPGDFILE,YPGDFILETYPE)
!
 CALL PREP_SURF_TRIP(CSURF_FILETYPE)
!
 CALL FLAG_UPDATE(.FALSE.,.TRUE.,.FALSE.,.FALSE.)
!
!* opens the file
IF (CSURF_FILETYPE=='FA    ') THEN
  LFANOCOMPACT = .TRUE.
  IDATEF(1)=TTIME%TDATE%YEAR
  IDATEF(2)=TTIME%TDATE%MONTH
  IDATEF(3)=TTIME%TDATE%DAY
  IDATEF(4)=NINT(TTIME%TIME/3600.) 
  IDATEF(5)=NINT(TTIME%TIME/60.) - IDATEF(4) * 60 
  IDATEF(6)=1 
  IDATEF(7:11)=0  
  CALL FAITOU(IRET,NUNIT_FA,.TRUE.,CFILEOUT_FA,'NEW',.TRUE.,.FALSE.,IVERBFA,0,INB,CDNOMC)
  CALL FANDAR(IRET,NUNIT_FA,IDATEF)
END IF
!
IF (LWRITE_COORD) CALL GET_LONLAT_n(CSURF_FILETYPE)
!
!* writes into the file
 CALL IO_BUFF_CLEAN_n
 CALL WRITE_SURF_ATM_n(CSURF_FILETYPE,'PRE',LLAND_USE) !no pgd field
 CALL WRITE_DIAG_SURF_ATM_n(CSURF_FILETYPE,'ALL')
!
!* closes the file
IF (CSURF_FILETYPE=='FA    ') THEN
  CALL FAIRME(IRET,NUNIT_FA,'UNKNOWN')
END IF
!
!* add informations in the file
IF (CSURF_FILETYPE=='LFI   ' .AND. LMNH_COMPATIBLE) CALL WRITE_HEADER_MNH
!
!
!*    3.     Close parallelized I/O
!            ----------------------
!
WRITE(ILUOUT,*) ' '
WRITE(ILUOUT,*) '    -----------------------'
WRITE(ILUOUT,*) '    | PREP ENDS CORRECTLY |'
WRITE(ILUOUT,*) '    -----------------------'
!
WRITE(*,*) ' '
WRITE(*,*) '    -----------------------'
WRITE(*,*) '    | PREP ENDS CORRECTLY |'
WRITE(*,*) '    -----------------------'
!
CLOSE(ILUOUT)
!
 CALL DEALLOC_SURFEX
!
IF (ASSOCIATED(NWORK)) DEALLOCATE(NWORK)
IF (ASSOCIATED(XWORK)) DEALLOCATE(XWORK)
IF (ASSOCIATED(NWORK2)) DEALLOCATE(NWORK2)
IF (ASSOCIATED(XWORK2)) DEALLOCATE(XWORK2)
IF (ASSOCIATED(XWORK3)) DEALLOCATE(XWORK3)
IF (ASSOCIATED(NWORK_FULL)) DEALLOCATE(NWORK_FULL)
IF (ASSOCIATED(XWORK_FULL)) DEALLOCATE(XWORK_FULL)
IF (ASSOCIATED(NWORK2_FULL)) DEALLOCATE(NWORK2_FULL)
IF (ASSOCIATED(XWORK2_FULL)) DEALLOCATE(XWORK2_FULL)
!
IF (LHOOK) CALL DR_HOOK('PREP',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END PROGRAM PREP
