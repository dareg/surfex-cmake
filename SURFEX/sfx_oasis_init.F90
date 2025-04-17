!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
SUBROUTINE SFX_OASIS_INIT(HNAMELIST,KLOCAL_COMM,HINIT,LDXIOS_ARPIFS) 
!!
!!
!!     PURPOSE
!!     --------
!!
!!     If needed, initialize XIOS I/O scheme and coupled mode communication
!!
!!
!!     METHOD
!!     ------
!!
!!     Read Surfex Oasis namelist to set Oasis activation flag LOASIS.
!!
!!     If LDXIOS_ARPIFS is TRUE, Xios has been previously initialized by 
!!     host model, together with Oasis (provided Xios xml setting 'using_oasis' is
!!     consistent) and there is only to read and set RUNTIME
!!
!!     Otherwise, if MODD:LXIOS is True, Oasis will be setup by setting
!!     up Xios; if it is False, we init Oasis by oasis_init_comp,
!!     and let it provide a local communicator (using oasis_get_localcomm)
!!
!!     If Oasis is used, final step will read Runtime in Oasis namelist
!!
!!     Note : OASIS-MCT interface must be initialized before any DR_HOOK call
!!
!!
!!     EXTERNAL
!!     --------
!!
!!
!!     REFERENCE
!!     ---------
!!
!!     S. Valcke et al., 2013: OASIS-MCT User Guide 
!!     CERFACS, Toulouse, France, 50 pp.
!!     https://verc.enes.org/oasis/oasis-dedicated-user-support-1/documentation/oasis3-mct-user-guide
!!
!!     XIOS Reference guide - Yann Meurdesoif - 10/10/2014 - 
!!     svn co -r 515 http://forge.ipsl.jussieu.fr/ioserver/svn/XIOS/branchs/xios-1.0 <dir> ; cd <dir>/doc ; ....
!!
!!
!!     AUTHOR
!!     ------
!!
!!     B. Decharme, CNRM
!!
!!     MODIFICATION
!!     --------------
!!
!!     Original    10/2013
!!     S.Sénési    08/2015 - handle XIOS
!!     B.Decharme  09/2016 - no CALL ABORT if no namelist in Arpege
!!     S.Sénési    04/2019 - introduce LDXIOS_ARPIFS and LDXIOS_OASIS
!!     D.StMartin  11/2019 - bugfix case (LXIOS_ARPIFS, LXIOS_SFX, LOASIS) = (T, F, F)
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SFX_OASIS, ONLY : LOASIS, CMODEL_NAME, XRUNTIME

USE MODD_XIOS, ONLY : LXIOS_SFX=>LXIOS
USE MODD_XIAS, ONLY : LXIAS

USE MODI_ABOR1_SFX
USE MODD_SURFEX_MPI, ONLY : LSFX_MPI
!
!
#ifdef WXIOS
USE XIOS, ONLY : XIOS_INITIALIZE
#endif
!
#ifdef CPLOASIS
USE MOD_OASIS
#endif
USE MODD_SURFEX_HOST
!
IMPLICIT NONE
!
#ifdef CPLOASIS
INCLUDE 'mpif.h'
#endif
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=28), INTENT(IN )           :: HNAMELIST
INTEGER,           INTENT(INOUT)         :: KLOCAL_COMM ! value of local communicator
CHARACTER(LEN=3),  INTENT(IN),  OPTIONAL :: HINIT       ! choice of fields to initialize
LOGICAL,           INTENT(IN),  OPTIONAL :: LDXIOS_ARPIFS ! Did host model initalize XIos
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
CHARACTER(LEN=9)   :: YWORD, YTIMERUN
CHARACTER(LEN=1000):: YLINE, YFOUND
INTEGER            :: IERR, IWORK, IRANK
INTEGER            :: ICOMP_ID
INTEGER            :: ITIMERUN
LOGICAL            :: GFOUND
CHARACTER(LEN=3)   :: YINIT
LOGICAL            :: LLXIOS_ARPIFS 
!
!
!*       0.3   Declarations of namelist variables
!              ----------------------------------
!
NAMELIST/NAM_OASIS/LOASIS,CMODEL_NAME
!
!-------------------------------------------------------------------------------
!
! ATTENTION : Do not introduce DR_HOOK in this routine
!
!*       0.     Initialization:
!               ---------------
!
LOASIS      = .FALSE.
CMODEL_NAME = 'surfex'
XRUNTIME    = 0.0
!
YINIT = 'ALL'
IF(PRESENT(HINIT))YINIT=HINIT
!
LLXIOS_ARPIFS=.FALSE.
IF(PRESENT(LDXIOS_ARPIFS)) LLXIOS_ARPIFS=LDXIOS_ARPIFS
!-------------------------------------------------------------------------------
!
!*       1.     Read namelist:
!               --------------
!
IF(LEN_TRIM(HNAMELIST)/=0)THEN
!
  OPEN(UNIT=11,FILE=HNAMELIST,ACTION='READ',FORM="FORMATTED",POSITION="REWIND",STATUS='OLD',IOSTAT=IERR)   
!
  IF (IERR /= 0) THEN
    WRITE(*,'(A)' )'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    WRITE(*,'(A)' )' WARNING WARNING WARNING WARNING WARNING     '
    WRITE(*,'(A)' )' ---------------------------------------     '
    WRITE(*,'(2A)')'SFX_OASIS_INIT: SFX NAMELIST FILE NOT FOUND: ',TRIM(HNAMELIST)
    WRITE(*,'(A)' )'-------------------------------------------  '     
    WRITE(*,'(A)' )'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    IF (ASSOCIATED (YRSURFEX_HOST)) THEN
      CALL YRSURFEX_HOST%ABORT ('SFX_OASIS_INIT: SFX NAMELIST FILE NOT FOUND')
    ENDIF
    CALL ABORT
  ELSE
    READ (UNIT=11,NML=NAM_OASIS,IOSTAT=IERR)
    CLOSE(UNIT=11)
  ENDIF
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       2.     Setup OASIS (possibly via XIOS)
!               -------------------------------
!
IF (LXIOS_SFX) THEN
!
#ifdef WXIOS
   IF (.NOT. LLXIOS_ARPIFS ) THEN
      ! NOTE : XIOS_INITIALIZE will call OASIS_INIT_COMP and 
      ! OASIS_GET_LOCALCOMM if its own config file calls for Oasis
!$OMP SINGLE
      CALL XIOS_INITIALIZE(CMODEL_NAME, return_comm=KLOCAL_COMM)
      LXIAS=.TRUE.
!$OMP END SINGLE
      ENDIF
#else
!
   WRITE(*,*) 'SFX_OASIS_INIT : BINARY WAS NOT COMPILED WITH XIOS SUPPORT '
   CALL ABOR1_SFX('SFX_OASIS_INIT : BINARY WAS NOT COMPILED WITH XIOS SUPPORT')
!
#endif
!
!
ELSE  ! (Surfex doesn't use XIOS)
!
!
#ifdef CPLOASIS
  IF (LOASIS ) THEN
      ! Maybe Arpege uses Xios, in which case Xios already did an Oasis init
      IF (.NOT. LLXIOS_ARPIFS ) THEN
         IRANK=0
         CALL OASIS_INIT_COMP(ICOMP_ID,CMODEL_NAME,IERR)  
         IF (IERR/=OASIS_OK) THEN
           WRITE(*,'(A)'   )'SFX : Error initializing OASIS'
           WRITE(*,'(A,I4)')'SFX : Return code from oasis_init_comp : ',IERR
           CALL OASIS_ABORT(ICOMP_ID,CMODEL_NAME,'SFX_OASIS_INIT: Error initializing OASIS')
           CALL ABORT
           STOP
         ENDIF
         CALL OASIS_GET_LOCALCOMM(KLOCAL_COMM,IERR) 
         IF (IERR/=OASIS_OK) THEN
            IF(IRANK==0)THEN
            WRITE(*,'(A)'   )'SFX : Error getting local communicator from OASIS'
            WRITE(*,'(A,I4)')'SFX : Return code from oasis_get_local_comm : ',IERR
           ENDIF
           CALL OASIS_ABORT(ICOMP_ID,CMODEL_NAME,'SFX_OASIS_INIT: Error getting local communicator')
           CALL ABORT
           STOP
         ENDIF
      ENDIF
!
  ELSE
    IF (.NOT. LLXIOS_ARPIFS ) THEN
      ! IF ARPEGE USES XIOS, LOCALCOMM IS ALREADY SET IN INIT_XIOS
      KLOCAL_COMM=0
      RETURN
    ENDIF
  ENDIF
#else
  IF (.NOT. LLXIOS_ARPIFS ) THEN
    ! IF ARPEGE USES XIOS, LOCALCOMM IS ALREADY SET IN INIT_XIOS
    KLOCAL_COMM=0
    RETURN
  ENDIF
#endif
!
ENDIF
!
#ifdef SFX_MPI
IF (LSFX_MPI) CALL MPI_COMM_RANK(KLOCAL_COMM,IRANK,IWORK)
#endif
IF(IRANK==0)THEN
   WRITE(*,'(A)')'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   IF (LOASIS) WRITE(*,'(A)')'OASIS used for model : '//TRIM(CMODEL_NAME)
   IF (LXIOS_SFX)  WRITE(*,'(A)')'XIOS  used for model : '//TRIM(CMODEL_NAME)
   IF (LLXIOS_ARPIFS)  WRITE(*,'(A)')'XIOS  used for model : ARPIFS'
   WRITE(*,'(A)')'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
ENDIF
!
IF(YINIT=='PRE')THEN
   RETURN
ENDIF

#ifdef CPLOASIS
IF (LOASIS) THEN
!
!-------------------------------------------------------------------------------
!
!*       5.     Read total simulated time in namcouple
!               --------------------------------------
!
 OPEN (UNIT=11,FILE ='namcouple',STATUS='OLD',FORM ='FORMATTED',POSITION="REWIND",IOSTAT=IERR)
 IF (IERR /= 0) THEN
    IF(IRANK==0)THEN
       WRITE(*,'(A)'   )'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
       WRITE(*,'(A)'   )'SFX : OASIS namcouple not found'
       WRITE(*,'(A)'   )'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    ENDIF
    CALL ABORT
    STOP
 ENDIF
!
 YTIMERUN=' $RUNTIME'
 ITIMERUN=-1
!
 DO WHILE (ITIMERUN==-1)
    READ (UNIT = 11,FMT = '(A9)',IOSTAT=IERR) YWORD
    IF(IERR/=0)THEN
       IF(IRANK==0)THEN
          WRITE(*,'(A)'   )'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
          WRITE(*,'(A)'   )'SFX : Problem $RUNTIME empty in namcouple'
          WRITE(*,'(A)'   )'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
       ENDIF
       CALL ABORT
       STOP           
    ENDIF
    IF (YWORD==YTIMERUN)THEN
       READ (UNIT = 11,FMT = '(A1000)',IOSTAT=IERR) YLINE
       IF(IERR/=0)THEN
          IF(IRANK==0)THEN
             WRITE(*,'(A)'   )'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
             WRITE(*,'(A)'   )'SFX : Problem looking for $RUNTIME in namcouple'
             WRITE(*,'(A)'   )'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
          ENDIF
          CALL ABORT
          STOP           
       ENDIF
       CALL FOUND_TIMERUN (YLINE, YFOUND, 1000, GFOUND)
       IF (GFOUND) THEN
          READ (YFOUND,FMT = '(I100)',IOSTAT=IERR) ITIMERUN
          IF(IERR/=0)THEN
             IF(IRANK==0)THEN
                WRITE(*,'(A)'   )'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                WRITE(*,'(A)'   )'SFX : Problem reading $RUNTIME in namcouple'
                WRITE(*,'(2A)'  )'$RUNTIME = ', TRIM(YFOUND)
                WRITE(*,'(A)'   )'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
             ENDIF
             CALL ABORT
             STOP
          ENDIF
       ENDIF
    ENDIF
 ENDDO
 CLOSE(11)
!
 XRUNTIME = REAL(ITIMERUN)
!
ENDIF
#endif
!
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!
SUBROUTINE FOUND_TIMERUN(HIN, HOUT, KLEN, OFOUND)
!
IMPLICIT NONE
!
INTEGER ,          INTENT (IN   ) :: KLEN
CHARACTER (LEN=*), INTENT (INOUT) :: HIN 
CHARACTER (LEN=*), INTENT (INOUT) :: HOUT
LOGICAL,           INTENT (OUT  ) :: OFOUND
!
!* ---------------------------- Local declarations -------------------
!
CHARACTER(LEN=1), PARAMETER :: YBLANK = ' '
CHARACTER(LEN=1), PARAMETER :: YNADA  = '#'

CHARACTER(LEN=KLEN) :: YLINE
CHARACTER(LEN=KLEN) :: YWORK
!
INTEGER             :: ILEN
INTEGER             :: IERR
!
!
!*    1. Skip line if it is a comment
!        ----------------------------
!
DO WHILE (HIN(1:1)==YNADA)
   READ (UNIT = 11, FMT = '(A20)',IOSTAT=IERR) YLINE 
   IF(IERR/=0)THEN
       IF(IRANK==0)THEN
         WRITE(*,'(A)'   )'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         WRITE(*,'(A)'   )'SFX : Problem looking for $RUNTINE line in namcouple'
         WRITE(*,'(A)'   )'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!' 
       ENDIF
       CALL ABORT
       STOP           
   ENDIF
   HIN(1:KLEN) = YLINE(1:KLEN)
ENDDO 
!
!* Fill HOUT with blanks
!
HOUT = YBLANK
!
!* Fill temporary string and remove leading blanks
!
YWORK = ADJUSTL(HIN)
!
IF(LEN_TRIM(YWORK)<=0)THEN
   OFOUND = .FALSE.
   RETURN
ENDIF
!
!* Find the length of this set of characters
!
ILEN = INDEX(YWORK,YBLANK) - 1
!
!* Copy to HOUT
!
HOUT(1:ILEN) = YWORK(1:ILEN)
!
OFOUND = .TRUE.
!
END SUBROUTINE FOUND_TIMERUN
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SFX_OASIS_INIT

