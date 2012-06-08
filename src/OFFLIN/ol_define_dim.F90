SUBROUTINE OL_DEFINE_DIM(HPROGRAM, KLUOUT, KNI, KDIM1, HUNIT1, HUNIT2, &
                         PX, PY, KDIMS, KDDIM, HNAME_DIM, KNPATCH)
!     #######################################################
!!****  *OL_DEFINE_DIM* - 
!!
!!    PURPOSE
!!    -------
!!
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
!!	S. Faroux   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    06/2010 
!!      07/2011     add specific computation for IGN grid (B. Decharme)
!-------------------------------------------------------------------------------                         
!
USE MODD_SURF_ATM_GRID_n, ONLY : CGRID, NGRID_PAR, XGRID_PAR
USE MODD_IO_SURF_OL, ONLY: NMASK_IGN
!
USE MODN_IO_OFFLINE, ONLY : LWRITE_COORD
!
USE MODI_GET_GRID_DIM
USE MODI_GET_GRID_COORD
USE MODI_GET_COORD_n
!
USE MODE_GRIDTYPE_IGN
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
include 'netcdf.inc'
!
CHARACTER(LEN=6),  INTENT(IN)    :: HPROGRAM
INTEGER, INTENT(IN)              :: KLUOUT
INTEGER, INTENT(IN)              :: KNI
INTEGER, INTENT(OUT)             :: KDIM1
CHARACTER(LEN=13) , DIMENSION(:), INTENT(OUT) :: HUNIT1, HUNIT2
REAL,DIMENSION(:), POINTER                :: PX, PY
INTEGER, DIMENSION(:), POINTER            :: KDIMS, KDDIM
CHARACTER(LEN=100), DIMENSION(:), POINTER :: HNAME_DIM
INTEGER, OPTIONAL, INTENT(IN)    :: KNPATCH
!
REAL, DIMENSION(KNI)             :: ZXX, ZYY
CHARACTER(LEN=3)                 :: YTYPE
INTEGER                          :: INDIMS, IDIM2
INTEGER                          :: I, J, K, L
LOGICAL                          :: GRECT     ! T if rectangular grid
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('OL_DEFINE_DIM',0,ZHOOK_HANDLE)
!
KDIM1=0
IDIM2=0
!
IF (.NOT.LWRITE_COORD) THEN
  !
  IF ( CGRID.EQ.'CONF PROJ ' .OR. CGRID.EQ.'CARTESIAN '&
  .OR. CGRID.EQ.'LONLAT REG' .OR. CGRID.EQ.'IGN' ) THEN
    YTYPE='XY '
    IF (CGRID.EQ.'LONLAT REG') YTYPE='LON'
    CALL GET_GRID_DIM(CGRID,NGRID_PAR,XGRID_PAR,GRECT,KDIM1,IDIM2)
  ENDIF
  !
ENDIF
!
INDIMS = 2
IF ( KDIM1.NE.0       ) INDIMS = 3
IF ( PRESENT(KNPATCH) ) INDIMS = INDIMS + 1
!
ALLOCATE(KDIMS(INDIMS))
ALLOCATE(KDDIM(INDIMS))
ALLOCATE(HNAME_DIM(INDIMS))
!
IF ( KDIM1.NE.0 ) THEN
  KDIMS(1) = KDIM1
  KDIMS(2) = IDIM2
  IF (YTYPE.EQ.'LON') THEN
    HNAME_DIM(1) = 'lon'
    HNAME_DIM(2) = 'lat'
    HUNIT1(1)    = 'degrees_east'
    HUNIT2(1)    = 'degrees_north'
  ELSE
    HNAME_DIM(1) = 'xx'
    HNAME_DIM(2) = 'yy'
    HUNIT1(1)    = 'meters'
    HUNIT2(1)    = 'meters'
  ENDIF
  ALLOCATE(PX(KDIM1))
  ALLOCATE(PY(IDIM2))
ELSE
  KDIMS(1) = KNI
  HNAME_DIM(1) = 'Number_of_points' 
  IF (LWRITE_COORD) THEN
    ALLOCATE(PX(KNI))
    ALLOCATE(PY(KNI))
  ENDIF
ENDIF
!
IF (LWRITE_COORD) THEN
  CALL GET_COORD_n(HPROGRAM,KNI,PX,PY)
ELSEIF ( CGRID.EQ.'CONF PROJ '.OR. CGRID.EQ.'CARTESIAN '.OR. &
         CGRID.EQ.'LONLAT REG' ) THEN
  !
  CALL GET_GRID_COORD(KLUOUT,PX=ZXX,PY=ZYY)
  !
  IF (ASSOCIATED(PX)) THEN
    DO J=1,SIZE(PX)
      PX(J)=ZXX(J)
    ENDDO
  ENDIF
  IF (ASSOCIATED(PY)) THEN
    DO J=1,SIZE(PY)
      PY(J)=ZYY((J-1)*(KNI/SIZE(PY))+1)
    ENDDO
  ENDIF
!
ELSEIF(CGRID.EQ.'IGN       ')THEN
  !
  CALL GET_GRIDTYPE_IGN(XGRID_PAR,PX=ZXX,PY=ZYY,PXALL=PX,PYALL=PY)
  !
  IF (.NOT.ALLOCATED(NMASK_IGN))THEN
    ALLOCATE(NMASK_IGN(KNI))
    L=0
    DO J=1,SIZE(PY)    
      DO I=1,SIZE(PX)
        L=L+1
        DO K=1,KNI
          IF((ZXX(K)==PX(I)).AND.(ZYY(K)==PY(J)))THEN
            NMASK_IGN(K)=L
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  !
ENDIF
!
!
IF (PRESENT(KNPATCH)) THEN
  KDIMS     (INDIMS-1) = KNPATCH
  HNAME_DIM (INDIMS-1) = 'Number_of_Tile'
ENDIF
!
KDIMS     (INDIMS) = NF_UNLIMITED
HNAME_DIM (INDIMS) = 'time'
!
IF (LHOOK) CALL DR_HOOK('OL_DEFINE_DIM',1,ZHOOK_HANDLE)
!
END SUBROUTINE OL_DEFINE_DIM
