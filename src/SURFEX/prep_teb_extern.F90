!     #########
SUBROUTINE PREP_TEB_EXTERN(HPROGRAM,HSURF,HFILE,HFILETYPE,HFILEPGD,HFILEPGDTYPE,KLUOUT,PFIELD)
!     #################################################################################
!
USE MODD_TYPE_DATE_SURF
!
USE MODI_PREP_GRID_EXTERN
USE MODI_READ_SURF
USE MODI_CONVERT_COVER_TEB
USE MODI_INTERP_GRID
USE MODI_OPEN_AUX_IO_SURF
USE MODI_CLOSE_AUX_IO_SURF
USE MODI_TOWN_PRESENCE
USE MODI_OLD_NAME
!
USE MODD_PREP,       ONLY : CINGRID_TYPE, CINTERP_TYPE
USE MODD_PREP_TEB,   ONLY : XGRID_ROAD, XGRID_WALL, XGRID_ROOF, &
                            XTI_BLD, XWS_ROOF, XWS_ROAD
USE MODD_DATA_COVER_PAR, ONLY : JPCOVER
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
CHARACTER(LEN=7),   INTENT(IN)  :: HSURF     ! type of field
CHARACTER(LEN=28),  INTENT(IN)  :: HFILE     ! name of file
CHARACTER(LEN=6),   INTENT(IN)  :: HFILETYPE ! type of input file
CHARACTER(LEN=28),  INTENT(IN)  :: HFILEPGD     ! name of file
CHARACTER(LEN=6),   INTENT(IN)  :: HFILEPGDTYPE ! type of input file
INTEGER,            INTENT(IN)  :: KLUOUT    ! logical unit of output listing
REAL,DIMENSION(:,:), POINTER    :: PFIELD    ! field to interpolate horizontally
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZFIELD         ! field read
REAL, DIMENSION(:,:), ALLOCATABLE :: ZDEPTH         ! depth of each layer
REAL, DIMENSION(:),   ALLOCATABLE :: ZDEPTH_TOT     ! total depth of surface
!
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZD  ! intermediate array
!
CHARACTER(LEN=16) :: YRECFM         ! Name of the article to be read
INTEGER           :: IRESP          ! reading return code
INTEGER           :: ILAYER         ! number of layers
INTEGER           :: JLAYER         ! loop counter
!
INTEGER           :: INI            ! total 1D dimension
!
LOGICAL, DIMENSION(JPCOVER)          :: GCOVER ! flag to read the covers
REAL,    DIMENSION(:,:), ALLOCATABLE :: ZCOVER ! cover fractions
LOGICAL                              :: GTEB      ! flag if TEB fields are present
INTEGER                              :: IPATCH    ! number of soil temperature patches
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
!*      1.     Preparation of IO for reading in the file
!              -----------------------------------------
!
!* Note that all points are read, even those without physical meaning.
!  These points will not be used during the horizontal interpolation step.
!  Their value must be defined as XUNDEF.
!
IF (LHOOK) CALL DR_HOOK('PREP_TEB_EXTERN',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
!*      2.     Reading of grid
!              ---------------
!
!* reads the grid
CALL OPEN_AUX_IO_SURF(HFILEPGD,HFILEPGDTYPE,'TOWN  ')
CALL PREP_GRID_EXTERN(HFILEPGDTYPE,KLUOUT,CINGRID_TYPE,CINTERP_TYPE,INI)
!
!
!* reads if TEB fields exist in the input file
CALL TOWN_PRESENCE(HFILEPGDTYPE,GTEB)
!
!---------------------------------------------------------------------------------------
!
!*     3.      Orography
!              ---------
!
IF (HSURF=='ZS     ') THEN
  ALLOCATE(PFIELD(INI,1))
  YRECFM='ZS'
  CALL READ_SURF(HFILEPGDTYPE,YRECFM,PFIELD(:,1),IRESP,HDIR='A')
  CALL CLOSE_AUX_IO_SURF(HFILEPGD,HFILEPGDTYPE)  
!
!---------------------------------------------------------------------------------------
ELSE
!---------------------------------------------------------------------------------------
!
!*     4.     TEB fields are read
!             -------------------
!
 IF (GTEB) THEN
!---------------------------------------------------------------------------------------
  SELECT CASE(HSURF)
!---------------------------------------------------------------------------------------
!
!*     4.1    Profile of temperatures in roads, roofs or walls
!             ------------------------------------------------
!
  CASE('T_ROAD','T_ROOF','T_WALL')
    !
    !* reading of number of layers
    IF (HSURF=='T_ROAD') YRECFM='ROAD_LAYER'
    IF (HSURF=='T_ROOF') YRECFM='ROOF_LAYER'
    IF (HSURF=='T_WALL') YRECFM='WALL_LAYER'
    CALL READ_SURF(HFILEPGDTYPE,YRECFM,ILAYER,IRESP)
    !
    !* reading of the cover to obtain the thickness of layers
    CALL OLD_NAME(HFILEPGDTYPE,'COVER_LIST      ',YRECFM)
    CALL READ_SURF(HFILEPGDTYPE,YRECFM,GCOVER(:),IRESP,HDIR='-')
    !
    ALLOCATE(ZCOVER(INI,JPCOVER))
    YRECFM='COVER'
    CALL READ_SURF(HFILEPGDTYPE,YRECFM,ZCOVER(:,:),GCOVER,IRESP,HDIR='A')
    !
    CALL CLOSE_AUX_IO_SURF(HFILEPGD,HFILEPGDTYPE)
    !
    !* reading of the profile
    ALLOCATE(ZFIELD(INI,ILAYER))
    CALL OPEN_AUX_IO_SURF(HFILE,HFILETYPE,'TOWN  ')
    DO JLAYER=1,ILAYER
      WRITE(YRECFM,'(A6,I1.1,A9)') HSURF,JLAYER,'         '
      CALL READ_SURF(HFILETYPE,YRECFM,ZFIELD(:,JLAYER),IRESP,HDIR='A')
    END DO
    CALL CLOSE_AUX_IO_SURF(HFILE,HFILETYPE)
    !
    ALLOCATE(ZD(INI,ILAYER))
    IF (HSURF=='T_ROAD') CALL CONVERT_COVER_TEB(ZCOVER,PD_ROAD=ZD)
    IF (HSURF=='T_ROOF') CALL CONVERT_COVER_TEB(ZCOVER,PD_ROOF=ZD)
    IF (HSURF=='T_WALL') CALL CONVERT_COVER_TEB(ZCOVER,PD_WALL=ZD)
    DEALLOCATE(ZCOVER)
    !
    !* recovers middle layer depth (from the surface)
    ALLOCATE(ZDEPTH    (INI,ILAYER))
    ALLOCATE(ZDEPTH_TOT(INI))
    ZDEPTH    (:,1)=ZD(:,1)/2.
    ZDEPTH_TOT(:)  =ZD(:,1)
    DO JLAYER=2,ILAYER
      ZDEPTH    (:,JLAYER) = ZDEPTH_TOT(:) + ZD(:,JLAYER)/2.
      ZDEPTH_TOT(:) = ZDEPTH_TOT(:) + ZD(:,JLAYER)
    END DO
    !
    !* in case of wall or roof, normalizes by total wall or roof thickness
    IF (HSURF=='T_ROOF' .OR. HSURF=='T_WALL') THEN
      DO JLAYER=1,ILAYER
        ZDEPTH(:,JLAYER) = ZDEPTH(:,JLAYER) / ZDEPTH_TOT(:)
      END DO
    END IF
    !
    !* interpolation on the fine vertical grid
    IF (HSURF=='T_ROAD') THEN
      ALLOCATE(PFIELD(SIZE(ZFIELD,1),SIZE(XGRID_ROAD)))
      CALL INTERP_GRID(ZDEPTH,ZFIELD,XGRID_ROAD,PFIELD)
    END IF
    IF (HSURF=='T_ROOF') THEN
      ALLOCATE(PFIELD(SIZE(ZFIELD,1),SIZE(XGRID_ROOF)))
      CALL INTERP_GRID(ZDEPTH,ZFIELD,XGRID_ROOF,PFIELD)
    END IF
    IF (HSURF=='T_WALL') THEN
      ALLOCATE(PFIELD(SIZE(ZFIELD,1),SIZE(XGRID_WALL)))
      CALL INTERP_GRID(ZDEPTH,ZFIELD,XGRID_WALL,PFIELD)
    END IF
    !
    !* end
    DEALLOCATE(ZD)
    DEALLOCATE(ZFIELD)
    DEALLOCATE(ZDEPTH)
    DEALLOCATE(ZDEPTH_TOT)
!
!---------------------------------------------------------------------------------------
!
!*      4.2    Other variables
!              ---------------
!
  CASE DEFAULT
    ALLOCATE(PFIELD(INI,1))
    YRECFM=HSURF
    IF (HSURF=='T_CAN  ') YRECFM='T_CANYON'
    IF (HSURF=='Q_CAN  ') YRECFM='Q_CANYON'
    CALL CLOSE_AUX_IO_SURF(HFILEPGD,HFILEPGDTYPE)
    CALL OPEN_AUX_IO_SURF(HFILE,HFILETYPE,'TOWN  ')    
    CALL READ_SURF(HFILETYPE,YRECFM,PFIELD(:,1),IRESP,HDIR='A')
    CALL CLOSE_AUX_IO_SURF(HFILE,HFILETYPE)    
!
!---------------------------------------------------------------------------------------
  END SELECT
!---------------------------------------------------------------------------------------
!
!*     5.     Subtitutes if TEB fields do not exist
!             -------------------------------------
!
 ELSE

  CALL CLOSE_AUX_IO_SURF(HFILEPGD,HFILEPGDTYPE)

  SELECT CASE(HSURF)

   !* temperature profiles
   CASE('T_ROAD','T_ROOF','T_WALL')
    !* reading of the soil surface temperature
    CALL OPEN_AUX_IO_SURF(HFILEPGD,HFILEPGDTYPE,'NATURE')    
    CALL READ_SURF(HFILEPGDTYPE,'PATCH_NUMBER',IPATCH,IRESP)
    CALL CLOSE_AUX_IO_SURF(HFILEPGD,HFILEPGDTYPE)
    ALLOCATE(ZFIELD(INI,IPATCH))
    CALL OPEN_AUX_IO_SURF(HFILE,HFILETYPE,'NATURE')
    CALL READ_SURF(HFILETYPE,'TG1',ZFIELD(:,:),IRESP,HDIR='A')
    CALL CLOSE_AUX_IO_SURF(HFILE,HFILETYPE)
    !* fills the whole temperature profile by this soil temperature
    IF (HSURF=='T_ROAD') ILAYER=SIZE(XGRID_ROAD)
    IF (HSURF=='T_ROOF') ILAYER=SIZE(XGRID_ROOF)
    IF (HSURF=='T_WALL') ILAYER=SIZE(XGRID_WALL)
    ALLOCATE(PFIELD(INI,ILAYER))
    DO JLAYER=1,ILAYER
      PFIELD(:,JLAYER) = ZFIELD(:,1)
    END DO
    DEALLOCATE(ZFIELD)

   !* water reservoirs
   CASE('WS_ROOF','WS_ROAD')
    ALLOCATE(PFIELD(INI,1))
    IF (HSURF=='WS_ROOF') PFIELD = XWS_ROOF
    IF (HSURF=='WS_ROAD') PFIELD = XWS_ROAD

   !* canyon or deep road temperature
   CASE('T_CAN','TI_ROAD')
    !* reading of the deep soil surface temperature
    CALL OPEN_AUX_IO_SURF(HFILEPGD,HFILEPGDTYPE,'NATURE')
    CALL READ_SURF(HFILEPGDTYPE,'PATCH_NUMBER',IPATCH,IRESP)
    CALL CLOSE_AUX_IO_SURF(HFILEPGD,HFILEPGDTYPE)
    ALLOCATE(ZFIELD(INI,IPATCH))
    CALL OPEN_AUX_IO_SURF(HFILE,HFILETYPE,'NATURE')
    CALL READ_SURF(HFILETYPE,'TG2',ZFIELD(:,:),IRESP,HDIR='A')
    CALL CLOSE_AUX_IO_SURF(HFILE,HFILETYPE)
    !* sets the temperature equal to this deep soil temperature
    ALLOCATE(PFIELD(INI,1))
    PFIELD(:,1) = ZFIELD(:,1)
    DEALLOCATE(ZFIELD)

   !* building temperature
   CASE('TI_BLD')
    ALLOCATE(PFIELD(INI,1))
    PFIELD = XTI_BLD

   !* other fields
   CASE DEFAULT
    ALLOCATE(PFIELD(INI,1))
    PFIELD = 0.

  END SELECT
 END IF
!-------------------------------------------------------------------------------------
END IF
!-------------------------------------------------------------------------------------
!
!*      6.     End of IO
!              ---------
!
IF (LHOOK) CALL DR_HOOK('PREP_TEB_EXTERN',1,ZHOOK_HANDLE)
!
!---------------------------------------------------------------------------------------
!
END SUBROUTINE PREP_TEB_EXTERN
