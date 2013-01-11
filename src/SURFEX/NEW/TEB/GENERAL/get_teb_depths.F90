!     #########
      SUBROUTINE GET_TEB_DEPTHS(HFILEPGDTYPE, PD_ROOF, PD_ROAD, PD_WALL, PD_FLOOR)
!     ##############################################################
!
!!**** *CONVERT_COVER* 
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
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
!!    V. Masson        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    01/2004
!     
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_DATA_COVER,     ONLY : XDATA_D_ROOF, XDATA_D_ROAD, XDATA_D_WALL, XDATA_D_FLOOR
USE MODD_DATA_COVER_PAR, ONLY : JPCOVER
!
USE MODI_READ_SURF
USE MODI_AV_PGD
USE MODI_OLD_NAME
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
CHARACTER(LEN=6),   INTENT(IN)  :: HFILEPGDTYPE ! type of input file
!
REAL, DIMENSION(:,:), INTENT(OUT), OPTIONAL   :: PD_ROOF
REAL, DIMENSION(:,:), INTENT(OUT), OPTIONAL   :: PD_ROAD
REAL, DIMENSION(:,:), INTENT(OUT), OPTIONAL   :: PD_WALL
REAL, DIMENSION(:,:), INTENT(OUT), OPTIONAL   :: PD_FLOOR
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
LOGICAL, DIMENSION(JPCOVER)          :: GCOVER ! flag to read the covers
REAL,    DIMENSION(:,:), ALLOCATABLE :: ZCOVER ! cover fractions
REAL,    DIMENSION(:,:), ALLOCATABLE :: ZD     ! depth of surface layers

INTEGER           :: IVERSION       ! surface version
INTEGER           :: IBUGFIX        ! surface bugfix version
CHARACTER(LEN=16) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=16) :: YRECFM1        ! Name of the article to be read
CHARACTER(LEN=16) :: YRECFM2        ! Name of the article to be read
CHARACTER(LEN=3)  :: YAREA          ! Area where field is to be averaged
INTEGER           :: IRESP          ! reading return code
LOGICAL           :: GDATA          ! T if depth is to be read in the file
REAL, DIMENSION(SIZE(XDATA_D_ROOF,1),SIZE(XDATA_D_ROOF,2)) :: ZDATA
INTEGER :: ILAYER                   ! number on surface layers
INTEGER :: JLAYER                   ! loop counter on surface layers
INTEGER :: ILU                      ! number of points
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!*    2.      SECONDARY VARIABLES
!             -------------------
!
!*    2.2     fields on artificial surfaces only
!             ----------------------------------
!
IF (LHOOK) CALL DR_HOOK('GET_TEB_DEPTHS',0,ZHOOK_HANDLE)
!
YRECFM='VERSION'
CALL READ_SURF(HFILEPGDTYPE,YRECFM,IVERSION,IRESP)
YRECFM='BUG'
CALL READ_SURF(HFILEPGDTYPE,YRECFM,IBUGFIX,IRESP)
!
IF (PRESENT(PD_ROOF)) THEN
  ZDATA = XDATA_D_ROOF
  YRECFM1 = 'L_D_ROOF'
  YRECFM2 = 'D_D_ROOF'
  ILU     = SIZE(PD_ROOF,1)
  ILAYER  = SIZE(PD_ROOF,2)
  YAREA   = 'BLD'
END IF
IF (PRESENT(PD_WALL)) THEN
  ZDATA = XDATA_D_WALL
  YRECFM1 = 'L_D_WALL'
  YRECFM2 = 'D_D_WALL'
  ILU     = SIZE(PD_WALL,1)
  ILAYER  = SIZE(PD_WALL,2)
  YAREA   = 'BLD'
END IF
IF (PRESENT(PD_ROAD)) THEN
  ZDATA = XDATA_D_ROAD
  YRECFM1 = 'L_D_ROAD'
  YRECFM2 = 'D_D_ROAD'
  ILU     = SIZE(PD_ROAD,1)
  ILAYER  = SIZE(PD_ROAD,2)
  YAREA   = 'STR'
END IF
IF (PRESENT(PD_FLOOR)) THEN
  ZDATA = XDATA_D_FLOOR
  YRECFM1 = 'L_D_FLOOR'
  YRECFM2 = 'D_D_FLOOR'
  ILU     = SIZE(PD_FLOOR,1)
  ILAYER  = SIZE(PD_FLOOR,2)
  YAREA   = 'BLD'
END IF

ALLOCATE(ZD(ILU,ILAYER))
!
!* read if the depths description are written in the file
IF (IVERSION<7 .OR. (IVERSION==7 .AND. IBUGFIX<=2)) THEN
  GDATA = .FALSE.
ELSE
  CALL READ_SURF(HFILEPGDTYPE,YRECFM1,GDATA,IRESP)
END IF
!
!* depths are read in the file
IF (GDATA) THEN
  DO JLAYER=1,ILAYER
    WRITE(YRECFM,FMT='(A,I1)') TRIM(YRECFM2),JLAYER
    CALL READ_SURF(HFILEPGDTYPE,YRECFM,ZD(:,JLAYER),IRESP,HDIR='A')
  END DO
!
ELSE
!* depths are deduced from the cover types
  !* reading of the cover to obtain the thickness of layers
  CALL OLD_NAME(HFILEPGDTYPE,'COVER_LIST      ',YRECFM)
  CALL READ_SURF(HFILEPGDTYPE,YRECFM,GCOVER(:),IRESP,HDIR='-')
  !* reading of the cover fractions
  ALLOCATE(ZCOVER(ILU,JPCOVER))
  YRECFM='COVER'
  CALL READ_SURF(HFILEPGDTYPE,YRECFM,ZCOVER(:,:),GCOVER,IRESP,HDIR='A')
  !
  !* deduces the depths of each layer
  DO JLAYER=1,ILAYER
    CALL AV_PGD (ZD(:,JLAYER), ZCOVER, ZDATA(:,JLAYER),YAREA,'ARI')
  END DO
  DEALLOCATE(ZCOVER)
ENDIF
!
IF (PRESENT(PD_ROOF )) PD_ROOF  = ZD
IF (PRESENT(PD_WALL )) PD_WALL  = ZD
IF (PRESENT(PD_ROAD )) PD_ROAD  = ZD
IF (PRESENT(PD_FLOOR)) PD_FLOOR = ZD
!
DEALLOCATE(ZD)
!
IF (LHOOK) CALL DR_HOOK('GET_TEB_DEPTHS',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE GET_TEB_DEPTHS
