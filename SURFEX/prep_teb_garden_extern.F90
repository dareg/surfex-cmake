!     #########
SUBROUTINE PREP_TEB_GARDEN_EXTERN(HPROGRAM,HSURF,HFILE,HFILETYPE,KLUOUT,PFIELD)
!     #################################################################################
!
!!****  *PREP_TEB_GARDEN_EXTERN* - initializes ISBA fields from operational GRIB
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     V. Masson 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!------------------------------------------------------------------
!

!
USE MODE_READ_EXTERN
!
USE MODD_TYPE_DATE_SURF
!
USE MODI_PREP_GRID_EXTERN
USE MODI_READ_SURF
USE MODI_INTERP_GRID
USE MODI_OPEN_AUX_IO_SURF
USE MODI_CLOSE_AUX_IO_SURF
!
USE MODD_PREP,           ONLY : CINGRID_TYPE, CINTERP_TYPE
USE MODD_PREP_TEB_GARDEN,ONLY : XGRID_SOIL, XWR_DEF
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE
USE MODD_SURF_PAR,       ONLY : XUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_PUT_ON_ALL_VEGTYPES
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
CHARACTER(LEN=7),   INTENT(IN)  :: HSURF     ! type of field
CHARACTER(LEN=28),  INTENT(IN)  :: HFILE     ! name of file
CHARACTER(LEN=6),   INTENT(IN)  :: HFILETYPE ! type of input file
INTEGER,            INTENT(IN)  :: KLUOUT    ! logical unit of output listing
REAL,DIMENSION(:,:,:), POINTER  :: PFIELD    ! field to interpolate horizontally (on final soil grid)
!
!*      0.2    declarations of local variables
!
CHARACTER(LEN=16) :: YRECFM         ! Name of the article to be read
INTEGER           :: IRESP          ! reading return code
INTEGER           :: INI            ! total 1D dimension
INTEGER           :: IPATCH         ! number of patch
!
REAL, DIMENSION(:,:,:), POINTER     :: ZFIELD         ! field read on initial MNH vertical soil grid, all patches
REAL, DIMENSION(:,:),   POINTER     :: ZFIELD1        ! field read on initial MNH vertical soil grid, one patch
REAL, DIMENSION(:,:,:), POINTER     :: ZD             ! depth of field in the soil
REAL, DIMENSION(:,:), POINTER     :: ZD1            ! depth of field in the soil, one patch
REAL, DIMENSION(:,:), ALLOCATABLE   :: ZOUT         !
INTEGER                             :: JPATCH         ! loop counter for patch
CHARACTER(LEN=7)                    :: YSURF     ! type of field
LOGICAL                         :: GGARDEN   ! T if gardens are present in the file
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!------------------------------------------------------------------------------
!
!*      1.     Preparation of IO for reading in the file
!              -----------------------------------------
!
!* Note that all points are read, even those without physical meaning.
!  These points will not be used during the horizontal interpolation step.
!  Their value must be defined as XUNDEF.
!
IF (LHOOK) CALL DR_HOOK('PREP_TEB_GARDEN_EXTERN',0,ZHOOK_HANDLE)
CALL OPEN_AUX_IO_SURF(HFILE,HFILETYPE,'TOWN  ')
!
!------------------------------------------------------------------------------
!
!*      2.     Reading of grid
!              ---------------
!
CALL PREP_GRID_EXTERN(HFILETYPE,KLUOUT,CINGRID_TYPE,CINTERP_TYPE,INI)
!
!
!---------------------------------------------------------------------------------------
!
!*      3.     Transformation into physical quantity to be interpolated
!              --------------------------------------------------------
!
SELECT CASE(HSURF)
!
!*     3.      Orography
!              ---------
!
  CASE('ZS     ')
    ALLOCATE(PFIELD(INI,1,1))
    YRECFM='ZS'
    CALL READ_SURF(HPROGRAM,YRECFM,PFIELD(:,1,1),IRESP,HDIR='A')
!
!--------------------------------------------------------------------------
!
!
!*      3.1    Profile of temperature, water or ice in the soil
!
  CASE('TG    ','WG    ','WGI   ')
!* choice if one reads garden fields (if present) or ISBA fields
    CALL READ_SURF(HPROGRAM,'GARDEN',GGARDEN,IRESP)
     IF (GGARDEN) THEN
       YSURF = 'TWN_'//HSURF(1:3)
     ELSE
       YSURF = HSURF
     END IF
!* reading of the profile and its depth definition
     CALL READ_EXTERN_ISBA(HPROGRAM,KLUOUT,INI,YSURF,ZFIELD,ZD)
! 
     ALLOCATE(ZFIELD1(SIZE(ZFIELD,1),SIZE(ZFIELD,2)))
     ALLOCATE(ZD1(SIZE(ZFIELD,1),SIZE(ZFIELD,2)))
     ALLOCATE(ZOUT(SIZE(ZFIELD,1),SIZE(XGRID_SOIL)))
     ALLOCATE(PFIELD(SIZE(ZFIELD,1),SIZE(XGRID_SOIL),SIZE(ZFIELD,3)))
!
     DO JPATCH=1,SIZE(ZFIELD,3)
        ZFIELD1(:,:)=ZFIELD(:,:,JPATCH)
        ZD1(:,:)=ZD(:,:,JPATCH)
        CALL INTERP_GRID(ZD1,ZFIELD1,XGRID_SOIL,ZOUT)
        PFIELD(:,:,JPATCH)=ZOUT(:,:)
     END DO
!
     DEALLOCATE(ZFIELD)
     DEALLOCATE(ZOUT)
     DEALLOCATE(ZFIELD1)
     DEALLOCATE(ZD)
!
!--------------------------------------------------------------------------
!
!*      3.4    Water content intercepted on leaves, LAI
!
  CASE('WR     ')
     ALLOCATE(PFIELD(INI,1,NVEGTYPE))
     !* choice if one reads garden fields (if present) or ISBA fields
     CALL READ_SURF(HPROGRAM,'GARDEN',GGARDEN,IRESP)
     IF (GGARDEN) THEN
       IPATCH = 1             
       YRECFM = 'TWN_WR'
     ELSE
       YRECFM = 'PATCH_NUMBER'
       CALL READ_SURF(HFILETYPE,YRECFM,IPATCH,IRESP)
       YRECFM = 'WR'
     END IF
     ALLOCATE(ZFIELD(INI,1,IPATCH))
     CALL READ_SURF(HPROGRAM,YRECFM,ZFIELD(:,1,:),IRESP,HDIR='A')
     CALL PUT_ON_ALL_VEGTYPES(INI,1,1,NVEGTYPE,ZFIELD,PFIELD)
     DEALLOCATE(ZFIELD)
!
  CASE('LAI    ')
     ALLOCATE(PFIELD(INI,1,NVEGTYPE))
     PFIELD(:,:,:) = XUNDEF
!
END SELECT
!
!
!---------------------------------------------------------------------------
!
!*      6.     End of IO
!              ---------
!
CALL CLOSE_AUX_IO_SURF(HFILE,HFILETYPE)
IF (LHOOK) CALL DR_HOOK('PREP_TEB_GARDEN_EXTERN',1,ZHOOK_HANDLE)
!
!---------------------------------------------------------------------------
!---------------------------------------------------------------------------
END SUBROUTINE PREP_TEB_GARDEN_EXTERN
