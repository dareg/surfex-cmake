!#########
SUBROUTINE SFX_OASIS_PREP(HPROGRAM)
!###################################################
!
!!****  *SFX_OASIS_PREP* - Prepare grid areas and mask file for SFX-OASIS coupling
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
!!	B. Decharme   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2013
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
!
USE MODD_SURF_ATM_GRID_n,ONLY : XLAT, XLON, XMESH_SIZE
!
USE MODD_SURF_ATM_n,     ONLY : XSEA, XWATER, XNATURE, XTOWN, &
                                CSEA, CWATER, NDIM_FULL,      &
                                NR_NATURE
!
USE MODD_ISBA_n,         ONLY : CISBA, LWTD, LGLACIER, LGW, XGW
!
USE MODN_SFX_OASIS
USE MODD_SFX_OASIS
!
USE MODI_GET_LUOUT
USE MODI_ABOR1_SFX
USE MODI_GET_MESH_CORNER
USE MODI_UNPACK_SAME_RANK
USE MODI_SFX_OASIS_CHECK
!
#ifdef SFXOASIS
USE MOD_OASIS
#endif
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=6),        INTENT(IN) :: HPROGRAM    ! program calling surf. schemes
!
!
!*       0.2   Declarations of local parameter
!              -------------------------------
!
INTEGER,           PARAMETER  :: INC = 4    ! Number of grid-cell corners
!
CHARACTER(LEN=4),  PARAMETER  :: YSFX_LAND = 'slan'
CHARACTER(LEN=4),  PARAMETER  :: YSFX_QSB  = 'sdra'
CHARACTER(LEN=4),  PARAMETER  :: YSFX_GW   = 'sgw '
CHARACTER(LEN=4),  PARAMETER  :: YSFX_SEA  = 'ssea'
CHARACTER(LEN=4),  PARAMETER  :: YSFX_LAKE = 'slak'
!
!*       0.3   Declarations of local variables
!              -------------------------------
!
REAL,    DIMENSION(NDIM_FULL)       :: ZGW        ! frac groundwater
REAL,    DIMENSION(NDIM_FULL)       :: ZMASK_LAND ! land-sea mask for rrm coupling
REAL,    DIMENSION(NDIM_FULL)       :: ZMASK_GW   ! groundwater mask for rrm coupling
REAL,    DIMENSION(NDIM_FULL)       :: ZMASK_LAKE ! lake mask for ogcm coupling
REAL,    DIMENSION(NDIM_FULL)       :: ZMASK_SEA  ! sea-land mask for ogcm coupling
!
REAL,    DIMENSION(NDIM_FULL,1)     :: ZLON
REAL,    DIMENSION(NDIM_FULL,1)     :: ZLAT
REAL,    DIMENSION(NDIM_FULL,1)     :: ZAREA
INTEGER, DIMENSION(NDIM_FULL,1)     :: IMASK
!
REAL,    DIMENSION(NDIM_FULL,1,INC) :: ZCORNER_LON
REAL,    DIMENSION(NDIM_FULL,1,INC) :: ZCORNER_LAT
!
INTEGER, DIMENSION(2)          :: IVAR_SHAPE  ! indexes for the coupling field local dimension
!
INTEGER                        :: IPART_ID ! Local partition ID
INTEGER                        :: IERR     ! Error info
!
INTEGER                        :: ILUOUT, IFLAG
!
INTEGER                        :: JI, JC
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SFX_OASIS_PREP',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
#ifdef SFXOASIS
!-------------------------------------------------------------------------------
!
!
!*       0.     Initialize :
!               ------------
!
CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
CALL SFX_OASIS_CHECK(ILUOUT)
!
!-------------------------------------------------------------------------------
!
!*       2.     Get grid definition :
!               ---------------------
!
CALL GET_MESH_CORNER(ILUOUT,ZCORNER_LAT(:,1,:),ZCORNER_LON(:,1,:))
!
ZLON(:,1)=XLON(:)
ZLAT(:,1)=XLAT(:)
!
IF(LGW)THEN
  CALL UNPACK_SAME_RANK(NR_NATURE(:),XGW(:),ZGW(:))
  WHERE(ZGW(:)==XUNDEF)ZGW(:)=0.0
ELSE
  ZGW(:) = 0.0
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       3.     Comput masks :
!               --------------
!
ZMASK_LAND(:) = XNATURE(:)+XTOWN(:)
ZMASK_SEA (:) = XSEA   (:)
ZMASK_GW  (:) = ZGW    (:)
IF(CWATER=='FLAKE ')THEN
  ZMASK_LAKE(:) = XWATER (:)
ELSE
  ZMASK_LAKE(:) = XUNDEF
ENDIF
IF(LCPL_SEA.AND.LWATER)THEN
  ZMASK_SEA (:) = XSEA (:)+XWATER(:)
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       5.     Write grid definition :
!               -----------------------
!
!
!
CALL OASIS_START_GRIDS_WRITING(IFLAG)
!
!*       1.1    Grid definition for Land surface :
!               ----------------------------------
!
IF(LCPL_LAND)THEN  
!
  ZAREA(:,1) = XMESH_SIZE(:) * ZMASK_LAND(:)
  !0 = not masked ; 1 = masked
  WHERE(ZAREA(:,1)>0.0)
        IMASK(:,1) = 0
  ELSEWHERE
        IMASK(:,1) = 1
  ENDWHERE
  CALL OASIS_WRITE_GRID  (YSFX_LAND,NDIM_FULL,1,ZLON(:,:),ZLAT(:,:))  
  CALL OASIS_WRITE_CORNER(YSFX_LAND,NDIM_FULL,1,INC,ZCORNER_LON(:,:,:),ZCORNER_LAT(:,:,:))
  CALL OASIS_WRITE_AREA  (YSFX_LAND,NDIM_FULL,1,ZAREA(:,:))
  CALL OASIS_WRITE_MASK  (YSFX_LAND,NDIM_FULL,1,IMASK(:,:))
!
  IF(LCPL_GW)THEN
    WHERE(ZMASK_LAND(:)>0.0)
          ZAREA(:,1) = XMESH_SIZE(:) * (1.0-ZMASK_GW(:))
    ELSEWHERE
          ZAREA(:,1) = 0.0
    ENDWHERE
  ELSE
    ZAREA(:,1) = XMESH_SIZE(:) * ZMASK_LAND(:)
  ENDIF
  !0 = not masked ; 1 = masked
  WHERE(ZAREA(:,1)>0.0)
        IMASK(:,1) = 0
  ELSEWHERE
        IMASK(:,1) = 1
  ENDWHERE
  CALL OASIS_WRITE_GRID  (YSFX_QSB,NDIM_FULL,1,ZLON(:,:),ZLAT(:,:))  
  CALL OASIS_WRITE_CORNER(YSFX_QSB,NDIM_FULL,1,INC,ZCORNER_LON(:,:,:),ZCORNER_LAT(:,:,:))
  CALL OASIS_WRITE_AREA  (YSFX_QSB,NDIM_FULL,1,ZAREA(:,:))
  CALL OASIS_WRITE_MASK  (YSFX_QSB,NDIM_FULL,1,IMASK(:,:))
!
ENDIF
!
! groundwater surface coupling case
!
IF(LCPL_GW)THEN       
  ZAREA(:,1) = XMESH_SIZE(:) * ZMASK_GW(:)
  !0 = not masked ; 1 = masked
  WHERE(ZAREA(:,1)>0.0)
        IMASK(:,1) = 0
  ELSEWHERE
        IMASK(:,1) = 1
  ENDWHERE
  CALL OASIS_WRITE_GRID  (YSFX_GW,NDIM_FULL,1,ZLON(:,:),ZLAT(:,:))  
  CALL OASIS_WRITE_CORNER(YSFX_GW,NDIM_FULL,1,INC,ZCORNER_LON(:,:,:),ZCORNER_LAT(:,:,:))
  CALL OASIS_WRITE_AREA  (YSFX_GW,NDIM_FULL,1,ZAREA(:,:))
  CALL OASIS_WRITE_MASK  (YSFX_GW,NDIM_FULL,1,IMASK(:,:))
ENDIF
!
!*       1.2    Grid definition for lake surface :
!               ----------------------------------
!
IF(LCPL_LAKE)THEN
  ZAREA(:,1) = XMESH_SIZE(:) * ZMASK_LAKE(:)
  !0 = not masked ; 1 = masked
  WHERE(ZAREA(:,1)>0.0)
        IMASK(:,1) = 0
  ELSEWHERE
        IMASK(:,1) = 1
  ENDWHERE
  CALL OASIS_WRITE_GRID  (YSFX_LAKE,NDIM_FULL,1,ZLON(:,:),ZLAT(:,:))  
  CALL OASIS_WRITE_CORNER(YSFX_LAKE,NDIM_FULL,1,INC,ZCORNER_LON(:,:,:),ZCORNER_LAT(:,:,:))
  CALL OASIS_WRITE_AREA  (YSFX_LAKE,NDIM_FULL,1,ZAREA(:,:))
  CALL OASIS_WRITE_MASK  (YSFX_LAKE,NDIM_FULL,1,IMASK(:,:))
ENDIF
!
!*       1.3    Grid definition for sea/water :
!               -------------------------------
!
IF(LCPL_SEA)THEN     
  ZAREA(:,1) = XMESH_SIZE(:) * ZMASK_SEA(:)
  !0 = not masked ; 1 = masked
  WHERE(ZAREA(:,1)>0.0)
        IMASK(:,1) = 0
  ELSEWHERE
        IMASK(:,1) = 1
  ENDWHERE
  CALL OASIS_WRITE_GRID  (YSFX_SEA,NDIM_FULL,1,ZLON(:,:),ZLAT(:,:))  
  CALL OASIS_WRITE_CORNER(YSFX_SEA,NDIM_FULL,1,INC,ZCORNER_LON(:,:,:),ZCORNER_LAT(:,:,:))
  CALL OASIS_WRITE_AREA  (YSFX_SEA,NDIM_FULL,1,ZAREA(:,:))
  CALL OASIS_WRITE_MASK  (YSFX_SEA,NDIM_FULL,1,IMASK(:,:))
ENDIF
!
CALL OASIS_TERMINATE_GRIDS_WRITING()
!
CALL OASIS_ENDDEF(IERR)
!
IF(IERR/=OASIS_OK)THEN
   WRITE(ILUOUT,*)'SFX_OASIS_PREP: OASIS enddef problem, err = ',IERR
   CALL ABOR1_SFX('SFX_OASIS_PREP: OASIS enddef problem')
ENDIF
!
!-------------------------------------------------------------------------------
#endif
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('SFX_OASIS_PREP',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SFX_OASIS_PREP
