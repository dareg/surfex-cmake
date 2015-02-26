!     #########
      SUBROUTINE SET_VEGTYPES_FRACTIONS(HPROGRAM)
!     ##############################################################
!
!!**** *SET_VEGTYPES_FRACTIONS* monitor for averaging and interpolations of cover fractions
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!
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
!!    V. Masson        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    10/12/97
!!
!!       Modified 08/12/05, P. Le Moigne: user defined fields       
!!       Modified    07/11, R. Alkama   : 'netcdf' => 'offlin'       
!!                                        removes very small values due to computation precision
!!                   03/13, R. Alkama   : from 12 to 19 vegtypes
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_OL_FILEID, ONLY : XVAR_TO_FILEIN
!
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE
USE MODD_ISBA_n,         ONLY : NGROUND_LAYER
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_DATA_ISBA_n,    ONLY : XPAR_VEGTYPE, LDATA_VEGTYPE
!
#ifdef SFX_ASC
USE MODI_SET_SURFEX_FILE_NAME_ASC
#endif
#ifdef SFX_FA
USE MODI_SET_SURFEX_FILE_NAME_FA
#endif
#ifdef SFX_LFI
USE MODI_SET_SURFEX_FILE_NAME_LFI
#endif
#ifdef SFX_NC
USE MODI_SET_SURFEX_FILE_NAME_NC
#endif
!
USE MODI_GET_LUOUT
USE MODI_OPEN_NAMELIST
USE MODI_CLOSE_NAMELIST
!
USE MODI_OPEN_FILEIN_OL
USE MODI_CLOSE_FILEIN_OL
!
USE MODI_READ_FROM_SURFEX_FILE
!
USE MODE_POS_SURF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
USE MODI_EXTRAPOL_FIELDS
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
 CHARACTER(LEN=6),    INTENT(IN)    :: HPROGRAM     ! Type of program
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER               :: ILUOUT    ! output listing logical unit
INTEGER               :: ILUNAM    ! namelist file  logical unit
LOGICAL               :: GFOUND    ! true if namelist is found
!
INTEGER               :: JVEGTYPE    ! loop counter on patch
!
!*    0.3    Declaration of namelists
!            ------------------------
!
! name of files containing data
!
 CHARACTER(LEN=28)          :: CFNAM_VEGTYPE    ! fractions of each vegtypes
!
! types of file containing data
!
 CHARACTER(LEN=6)           :: CFTYP_VEGTYPE    ! fractions of each vegtypes
!
 CHARACTER(LEN=28)     :: HFILEIN
!
LOGICAL :: GOPEN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
NAMELIST/NAM_LAND_USE/CFNAM_VEGTYPE,CFTYP_VEGTYPE  

!-------------------------------------------------------------------------------
!
!*    1.      Initializations
!             ---------------
!
IF (LHOOK) CALL DR_HOOK('SET_VEGTYPES_FRACTIONS',0,ZHOOK_HANDLE)
CFNAM_VEGTYPE     = '                            '
CFTYP_VEGTYPE     = '      '
!
!-------------------------------------------------------------------------------
!
!*    2.      Input file for cover types
!             --------------------------
!
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
 CALL OPEN_NAMELIST(HPROGRAM,ILUNAM)
!
 CALL POSNAM(ILUNAM,'NAM_LAND_USE',GFOUND,ILUOUT)
IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_LAND_USE)
!
 CALL CLOSE_NAMELIST(HPROGRAM,ILUNAM)
!
!-------------------------------------------------------------------------------
!
!*    3.      Uniform fields are prescribed
!             -----------------------------
!
IF(CFTYP_VEGTYPE=='NETCDF') CFTYP_VEGTYPE='OFFLIN'
!
IF (CFTYP_VEGTYPE=='ASCII ') THEN
#ifdef SFX_ASC
  CALL SET_SURFEX_FILE_NAME_ASC(HNAME_OUT=HFILEIN)
#endif
ELSEIF (CFTYP_VEGTYPE=='FA    ') THEN
#ifdef SFX_FA
  CALL SET_SURFEX_FILE_NAME_FA(HNAME_OUT=HFILEIN)
#endif
ELSEIF (CFTYP_VEGTYPE=='LFI   ') THEN
#ifdef SFX_LFI
  CALL SET_SURFEX_FILE_NAME_LFI(HNAME_OUT=HFILEIN)
#endif
ELSEIF (CFTYP_VEGTYPE=='NC    ') THEN
#ifdef SFX_NC
  CALL SET_SURFEX_FILE_NAME_NC(HNAME_OUT=HFILEIN)
#endif
ENDIF
!
GOPEN = .FALSE.
IF(CFTYP_VEGTYPE=='OFFLIN' .AND. .NOT.ALLOCATED(XVAR_TO_FILEIN)) THEN
  GOPEN = .TRUE.
  CALL OPEN_FILEIN_OL
ENDIF
!
IF (CFTYP_VEGTYPE=='FA    '.OR.CFTYP_VEGTYPE=='ASCII '.OR.CFTYP_VEGTYPE=='LFI   ' &
        .OR.CFTYP_VEGTYPE=='OFFLIN' .OR.CFTYP_VEGTYPE=='NC    ') THEN
!        
  LDATA_VEGTYPE=.TRUE.
!
  CALL READ_FROM_SURFEX_FILE(CFTYP_VEGTYPE,CFNAM_VEGTYPE,'NATURE','      ',XPAR_VEGTYPE(:,1),HNAM='VEGTY_P1')
  CALL READ_FROM_SURFEX_FILE(CFTYP_VEGTYPE,CFNAM_VEGTYPE,'NATURE','      ',XPAR_VEGTYPE(:,2),HNAM='VEGTY_P2')
  CALL READ_FROM_SURFEX_FILE(CFTYP_VEGTYPE,CFNAM_VEGTYPE,'NATURE','      ',XPAR_VEGTYPE(:,3),HNAM='VEGTY_P3')
  CALL READ_FROM_SURFEX_FILE(CFTYP_VEGTYPE,CFNAM_VEGTYPE,'NATURE','      ',XPAR_VEGTYPE(:,4),HNAM='VEGTY_P4')
  CALL READ_FROM_SURFEX_FILE(CFTYP_VEGTYPE,CFNAM_VEGTYPE,'NATURE','      ',XPAR_VEGTYPE(:,5),HNAM='VEGTY_P5')
  CALL READ_FROM_SURFEX_FILE(CFTYP_VEGTYPE,CFNAM_VEGTYPE,'NATURE','      ',XPAR_VEGTYPE(:,6),HNAM='VEGTY_P6')
  CALL READ_FROM_SURFEX_FILE(CFTYP_VEGTYPE,CFNAM_VEGTYPE,'NATURE','      ',XPAR_VEGTYPE(:,7),HNAM='VEGTY_P7')
  CALL READ_FROM_SURFEX_FILE(CFTYP_VEGTYPE,CFNAM_VEGTYPE,'NATURE','      ',XPAR_VEGTYPE(:,8),HNAM='VEGTY_P8')
  CALL READ_FROM_SURFEX_FILE(CFTYP_VEGTYPE,CFNAM_VEGTYPE,'NATURE','      ',XPAR_VEGTYPE(:,9),HNAM='VEGTY_P9')
  CALL READ_FROM_SURFEX_FILE(CFTYP_VEGTYPE,CFNAM_VEGTYPE,'NATURE','      ',XPAR_VEGTYPE(:,10),HNAM='VEGTY_P10')
  CALL READ_FROM_SURFEX_FILE(CFTYP_VEGTYPE,CFNAM_VEGTYPE,'NATURE','      ',XPAR_VEGTYPE(:,11),HNAM='VEGTY_P11')
  CALL READ_FROM_SURFEX_FILE(CFTYP_VEGTYPE,CFNAM_VEGTYPE,'NATURE','      ',XPAR_VEGTYPE(:,12),HNAM='VEGTY_P12')
  CALL READ_FROM_SURFEX_FILE(CFTYP_VEGTYPE,CFNAM_VEGTYPE,'NATURE','      ',XPAR_VEGTYPE(:,13),HNAM='VEGTY_P13')
  CALL READ_FROM_SURFEX_FILE(CFTYP_VEGTYPE,CFNAM_VEGTYPE,'NATURE','      ',XPAR_VEGTYPE(:,14),HNAM='VEGTY_P14')
  CALL READ_FROM_SURFEX_FILE(CFTYP_VEGTYPE,CFNAM_VEGTYPE,'NATURE','      ',XPAR_VEGTYPE(:,18),HNAM='VEGTY_P15')
  CALL READ_FROM_SURFEX_FILE(CFTYP_VEGTYPE,CFNAM_VEGTYPE,'NATURE','      ',XPAR_VEGTYPE(:,15),HNAM='VEGTY_P16')
  CALL READ_FROM_SURFEX_FILE(CFTYP_VEGTYPE,CFNAM_VEGTYPE,'NATURE','      ',XPAR_VEGTYPE(:,16),HNAM='VEGTY_P17')
  CALL READ_FROM_SURFEX_FILE(CFTYP_VEGTYPE,CFNAM_VEGTYPE,'NATURE','      ',XPAR_VEGTYPE(:,17),HNAM='VEGTY_P18')
  CALL READ_FROM_SURFEX_FILE(CFTYP_VEGTYPE,CFNAM_VEGTYPE,'NATURE','      ',XPAR_VEGTYPE(:,19),HNAM='VEGTY_P19')
!
ENDIF
!
IF (GOPEN) CALL CLOSE_FILEIN_OL
!
! removes very small values due to computation precision
!
WHERE(XPAR_VEGTYPE < 1.E-8)XPAR_VEGTYPE(:,:)=0.0
!
IF (CFTYP_VEGTYPE=='ASCII ') THEN
#ifdef SFX_ASC
  CALL SET_SURFEX_FILE_NAME_ASC(HNAME_IN=HFILEIN)
#endif
ELSEIF (CFTYP_VEGTYPE=='FA    ') THEN
#ifdef SFX_FA
  CALL SET_SURFEX_FILE_NAME_FA(HNAME_IN=HFILEIN)
#endif
ELSEIF (CFTYP_VEGTYPE=='LFI   ') THEN
#ifdef SFX_LFI
  CALL SET_SURFEX_FILE_NAME_LFI(HNAME_IN=HFILEIN)
#endif
ELSEIF (CFTYP_VEGTYPE=='NC    ') THEN
#ifdef SFX_NC
  CALL SET_SURFEX_FILE_NAME_NC(HNAME_IN=HFILEIN)
#endif
ENDIF
!
IF (LDATA_VEGTYPE) THEN 
  IF (MAXVAL(ABS(SUM(XPAR_VEGTYPE,2)-1.))>1.E-6) THEN
    JVEGTYPE=COUNT(SUM(XPAR_VEGTYPE,2) .GT. 1.E19)
    WRITE(ILUOUT,*) ' '
    WRITE(ILUOUT,*) '******************************************************************************'
    WRITE(ILUOUT,*) '* Error in ISBA data field preparation                                       *'
    WRITE(ILUOUT,*) '* Sum of XPAR_VEGTYPE on all vegtypes is not equal to 1. for all grid point  *'
    WRITE(ILUOUT,*) '* nbr of indef VEGTYPE =',JVEGTYPE, ' /  total nbr =', SIZE(XPAR_VEGTYPE(:,1))    
    WRITE(ILUOUT,*) '* MAXVAL of SUM(XPAR_VEGTYPE,2) =', MAXVAL(SUM(XPAR_VEGTYPE,2))
    WRITE(ILUOUT,*) '* MAXLOC of SUM(XPAR_VEGTYPE,2) =', MAXLOC(SUM(XPAR_VEGTYPE,2))
    WRITE(ILUOUT,*) '******************************************************************************'
    WRITE(ILUOUT,*) ' '
    CALL ABOR1_SFX('SET_VEGTYPES_FRACTIONS: SUM OF ALL XPAR_VEGTYPE MUST BE 1.')
  ENDIF
ENDIF
!
IF (LDATA_VEGTYPE) CALL EXTRAPOL_FIELDS(HPROGRAM,ILUOUT)
!
IF (LHOOK) CALL DR_HOOK('SET_VEGTYPES_FRACTIONS',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SET_VEGTYPES_FRACTIONS
