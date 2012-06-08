!     ###########################################################
      SUBROUTINE ZOOM_PGD_ISBA_FULL(HPROGRAM,HINIFILE,HINIFILETYPE)
!     ###########################################################

!!
!!    PURPOSE
!!    -------
!!   This program prepares the physiographic data fields.
!!
!!    METHOD
!!    ------
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
!!    V. Masson                   Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original     13/10/03
!!    B. Decharme      2008  XWDRAIN
!!    M.Tomasini   17/04/12 All COVER physiographic fields are now 
!!                          interpolated for spawning => 
!!                          ABOR1_SFX if (.NOT.OECOCLIMAP) in comment
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_ISBA_n,           ONLY : XCLAY, XSAND, XRUNOFFB, XWDRAIN, &
                                  NGROUND_LAYER  
USE MODD_ISBA_GRID_n,      ONLY : XLAT, XLON, CGRID, XGRID_PAR, NDIM
USE MODD_PREP,             ONLY : CINGRID_TYPE, CINTERP_TYPE, LINTERP
!
USE MODI_GET_LUOUT
USE MODI_OPEN_AUX_IO_SURF
USE MODI_PREP_GRID_EXTERN
USE MODI_PREP_OUTPUT_GRID
USE MODI_READ_SURF
USE MODI_CLOSE_AUX_IO_SURF
USE MODI_HOR_INTERPOL
USE MODI_GET_TYPE_DIM_n
USE MODI_READ_PGD_ISBA_PAR_n
USE MODI_CLEAN_PREP_OUTPUT_GRID
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of dummy arguments
!            ------------------------------
!
CHARACTER(LEN=6),       INTENT(IN)  :: HPROGRAM     ! program calling
CHARACTER(LEN=28),      INTENT(IN)  :: HINIFILE     ! input atmospheric file name
CHARACTER(LEN=6),       INTENT(IN)  :: HINIFILETYPE ! input atmospheric file type
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER :: IRESP
INTEGER :: ILUOUT
INTEGER :: INI     ! total 1D dimension (input grid)
INTEGER :: JLAYER  ! loop counter
REAL, DIMENSION(:),   ALLOCATABLE :: ZFIELD    ! field read
REAL, DIMENSION(:,:), POINTER     :: ZSAND   ! sand   on all surface points
REAL, DIMENSION(:,:), POINTER     :: ZCLAY   ! clay   on all surface points
REAL, DIMENSION(:,:), POINTER     :: ZRUNOFFB! runoff coef. on all surface points
REAL, DIMENSION(:,:), POINTER     :: ZWDRAIN ! drainage coef. on all surface points
REAL, DIMENSION(:,:), ALLOCATABLE :: ZOUTB   ! runoff coef. on all surface points
REAL, DIMENSION(:,:), ALLOCATABLE :: ZOUTW   ! drainage coef. on all surface points
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ZOOM_PGD_ISBA_FULL',0,ZHOOK_HANDLE)
CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
!*      1.     Preparation of IO for reading in the file
!              -----------------------------------------
!
!* Note that all points are read, even those without physical meaning.
!  These points will not be used during the horizontal interpolation step.
!  Their value must be defined as XUNDEF.
!
CALL OPEN_AUX_IO_SURF(HINIFILE,HINIFILETYPE,'FULL  ')
!
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
ALLOCATE(ZSAND(INI,NGROUND_LAYER))
CALL READ_SURF(HPROGRAM,'SAND',ZFIELD,IRESP,HDIR='A')
DO JLAYER=1,NGROUND_LAYER
  ZSAND(:,JLAYER) = ZFIELD(:)
END DO
!
ALLOCATE(ZCLAY(INI,NGROUND_LAYER))
CALL READ_SURF(HPROGRAM,'CLAY',ZFIELD,IRESP,HDIR='A')
DO JLAYER=1,NGROUND_LAYER
  ZCLAY(:,JLAYER) = ZFIELD(:)
END DO
!
ALLOCATE(ZRUNOFFB(INI,1))
CALL READ_SURF(HPROGRAM,'RUNOFFB',ZFIELD,IRESP,HDIR='A')
ZRUNOFFB(:,1) = ZFIELD(:)
!
ALLOCATE(ZWDRAIN(INI,1))
CALL READ_SURF(HPROGRAM,'WDRAIN',ZFIELD,IRESP,HDIR='A')
ZWDRAIN(:,1) = ZFIELD(:)
!
DEALLOCATE(ZFIELD)
!
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
CALL HOR_INTERPOL(ILUOUT,ZSAND,XSAND)
CALL HOR_INTERPOL(ILUOUT,ZCLAY,XCLAY)
ALLOCATE(ZOUTB(SIZE(XRUNOFFB),1))
CALL HOR_INTERPOL(ILUOUT,ZRUNOFFB,ZOUTB)
XRUNOFFB(:) = ZOUTB(:,1)
DEALLOCATE(ZOUTB)
ALLOCATE(ZOUTW(SIZE(XWDRAIN),1))
CALL HOR_INTERPOL(ILUOUT,ZWDRAIN,ZOUTW)
XWDRAIN(:) = ZOUTW(:,1)
DEALLOCATE(ZOUTW)
!
CALL GET_TYPE_DIM_n('NATURE',NDIM)
CALL READ_PGD_ISBA_PAR_n(HPROGRAM,INI,.FALSE.,HDIR='A')
!
CALL CLOSE_AUX_IO_SURF(HINIFILE,HINIFILETYPE)
!
CALL CLEAN_PREP_OUTPUT_GRID
IF (LHOOK) CALL DR_HOOK('ZOOM_PGD_ISBA_FULL',1,ZHOOK_HANDLE)
!------------------------------------------------------------------------------
!
END SUBROUTINE ZOOM_PGD_ISBA_FULL
