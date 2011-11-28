!     ######spl
      SUBROUTINE ZOOM_PGD_SEAFLUX(HPROGRAM,HINIFILE,HINIFILETYPE,HFILE,HFILETYPE)
!     ##############################################################
!
!!**** *PGD_SEAFLUX* monitor for averaging and interpolations of SEAFLUX physiographic fields
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
!!    P. Le Moigne     Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    09/2008
!!    G. TANGUY   03/2009 : add reading and interpolation of XDATA_SST and 
!!                          TDATA_SST in the case LDATA_SST=T
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_PGD_GRID,       ONLY : NL
USE MODD_DATA_COVER_PAR,  ONLY : JPCOVER
USE MODD_SEAFLUX_n,       ONLY : XCOVER, LCOVER, XZS, XSEABATHY
USE MODD_SEAFLUX_GRID_n,  ONLY : CGRID, XGRID_PAR, XLAT, XLON, XMESH_SIZE, NDIM
USE MODD_TYPE_DATE_SURF

USE MODD_DATA_SEAFLUX_n,    ONLY : LSST_DATA,NTIME, XDATA_SST, TDATA_SST
!
USE MODD_PREP,             ONLY : CINGRID_TYPE, CINTERP_TYPE, LINTERP
!
USE MODI_READ_SURF
!
USE MODI_GET_SURF_SIZE_n
USE MODI_PACK_PGD
!
USE MODI_OPEN_AUX_IO_SURF
USE MODI_CLOSE_AUX_IO_SURF
USE MODI_GET_LUOUT
USE MODI_PREP_GRID_EXTERN
USE MODI_HOR_INTERPOL
USE MODI_PREP_OUTPUT_GRID
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_CLEAN_PREP_OUTPUT_GRID
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
CHARACTER(LEN=6),    INTENT(IN)    :: HPROGRAM     ! Type of program
CHARACTER(LEN=28),    INTENT(IN)  :: HINIFILE    ! input atmospheric file name
CHARACTER(LEN=6),     INTENT(IN)  :: HINIFILETYPE! input atmospheric file type
CHARACTER(LEN=28),    INTENT(IN)  :: HFILE       ! output file name
CHARACTER(LEN=6),     INTENT(IN)  :: HFILETYPE   ! output file type
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
!*    0.3    Declaration of namelists
!            ------------------------
!
CHARACTER(LEN=28)        :: YSEABATHY         ! file name for bathymetrie
CHARACTER(LEN=6)         :: YSEABATHYFILETYPE ! bathymetry data file type
CHARACTER(LEN=28)        :: YNCVARNAME        ! variable to read in netcdf
                                              ! file
REAL                     :: XUNIF_SEABATHY    ! uniform value of bathymetry
!
INTEGER :: IRESP
INTEGER :: ILUOUT
INTEGER :: INI

REAL, DIMENSION(:),   ALLOCATABLE   :: ZFIELD    ! field read
REAL, DIMENSION(:),   POINTER       :: ZSEABATHY
REAL, DIMENSION(:,:), POINTER       :: ZWORK1         
REAL, DIMENSION(:,:), ALLOCATABLE   :: ZWORK2 


CHARACTER(LEN=16) :: YRECFM         ! Name of the article to be read
INTEGER           :: JTIME          ! loop index
REAL, DIMENSION(:,:),   POINTER :: ZDATA_SST 

CHARACTER(LEN=100):: YCOMMENT       ! Comment string
REAL(KIND=JPRB) :: ZHOOK_HANDLE



!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ZOOM_PGD_SEAFLUX',0,ZHOOK_HANDLE)
CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
!*      1.     Preparation of IO for reading in the file
!              -----------------------------------------
!
!* Note that all points are read, even those without physical meaning.
!  These points will not be used during the horizontal interpolation step.
!  Their value must be defined as XUNDEF.
!
!
CALL OPEN_AUX_IO_SURF(HINIFILE,HINIFILETYPE,'FULL  ')
!

!-------------------------------------------------------------------------------
!
!*    1.      Initializations of defaults
!             ---------------------------
!
!-------------------------------------------------------------------------------
!
!*    2.      Reading of namelist
!             -------------------
!
!
!-------------------------------------------------------------------------------
!
!*    3.      Coherence of options
!             --------------------
!
!-------------------------------------------------------------------------------
!
!*    4.      Bathymetry
!             ----------
!
!-------------------------------------------------------------------------------
!
!*    5.      Number of points and packing
!             ----------------------------
!
CALL GET_SURF_SIZE_n('SEA   ',NDIM)
!
ALLOCATE(LCOVER     (JPCOVER))
ALLOCATE(XCOVER     (NDIM,JPCOVER))
ALLOCATE(XZS        (NDIM))
ALLOCATE(XLAT       (NDIM))
ALLOCATE(XLON       (NDIM))
ALLOCATE(XMESH_SIZE (NDIM))
!
CALL PACK_PGD(HPROGRAM, 'SEA   ',                    &
                CGRID,  XGRID_PAR, LCOVER,             &
                XCOVER, XZS,                           &
                XLAT, XLON, XMESH_SIZE                 )  
!
ALLOCATE(XSEABATHY (NDIM))

!------------------------------------------------------------------------------
!
!*      2.     Reading of grid
!              ---------------
!
CALL PREP_GRID_EXTERN(HINIFILETYPE,ILUOUT,CINGRID_TYPE,CINTERP_TYPE,INI)
!
CALL PREP_OUTPUT_GRID(ILUOUT,CGRID,XGRID_PAR,XLAT,XLON)
!
!------------------------------------------------------------------------------
!
!*      3.     Reading of fields
!              -----------------
!
!
ALLOCATE(ZFIELD(INI))
!
ALLOCATE(ZSEABATHY(INI))
CALL READ_SURF(HPROGRAM,'BATHY',ZFIELD,IRESP,HDIR='A')
ZSEABATHY(:) = ZFIELD(:)
!
DEALLOCATE(ZFIELD)
!


!============================================================
! G. TANGUY 03/2009
! reading of fields for SST_DATA
CALL READ_SURF(HPROGRAM,'SST_DATA',LSST_DATA,IRESP)
IF (LSST_DATA) THEN
  CALL READ_SURF(HPROGRAM,'NDATA_SEATIM',NTIME,IRESP)
  ALLOCATE(ZDATA_SST(INI,NTIME))


  DO JTIME=1,NTIME
    WRITE(YRECFM,FMT='(A9,I3.3)') 'DATA_SST_',JTIME
    YCOMMENT='(-)'
    CALL READ_SURF(HPROGRAM,YRECFM,ZDATA_SST(:,JTIME),IRESP,HDIR='A')
  END DO
   
 
    ALLOCATE(TDATA_SST(NTIME))   
  YRECFM='TDATA_SST'
  YCOMMENT='(-)'
  CALL READ_SURF(HPROGRAM,YRECFM,TDATA_SST,IRESP)
ENDIF
!============================================================


CALL CLOSE_AUX_IO_SURF(HINIFILE,HINIFILETYPE)
!------------------------------------------------------------------------------
!
!*      4.     Interpolations
!              --------------
!
!* mask where interpolations must be done
!
LINTERP(:) = .TRUE.
!
!* interpolations
!
ALLOCATE(ZWORK1(INI,1))
ALLOCATE(ZWORK2(NDIM,1))
ZWORK1(:,1)=ZSEABATHY(:)
CALL HOR_INTERPOL(ILUOUT,ZWORK1,ZWORK2)        
XSEABATHY(:)=ZWORK2(:,1)
DEALLOCATE(ZSEABATHY)
DEALLOCATE(ZWORK1,ZWORK2)
!

!============================================================
! G. TANGUY 03/2009
! interpolation of SST_DATA
IF (LSST_DATA) THEN
  ALLOCATE(XDATA_SST(NDIM,NTIME))
  DO JTIME=1,NTIME
    ALLOCATE(ZWORK1(INI,1))

    ALLOCATE(ZWORK2(NDIM,1))
    ZWORK1(:,1)=ZDATA_SST(:,JTIME)
    CALL HOR_INTERPOL(ILUOUT,ZWORK1,ZWORK2)
    XDATA_SST(:,JTIME)=ZWORK2(:,1)
    DEALLOCATE(ZWORK1)
    DEALLOCATE(ZWORK2)

  END DO
  DEALLOCATE(ZDATA_SST)

ENDIF

!============================================================


CALL CLEAN_PREP_OUTPUT_GRID
IF (LHOOK) CALL DR_HOOK('ZOOM_PGD_SEAFLUX',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE ZOOM_PGD_SEAFLUX
