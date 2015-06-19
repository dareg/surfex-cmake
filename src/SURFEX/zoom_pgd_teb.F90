!     ###########################################################
      SUBROUTINE ZOOM_PGD_TEB(HPROGRAM,HINIFILE,HINIFILETYPE,OECOCLIMAP,OGARDEN)
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
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
!
USE MODD_IO_BUFF_n, ONLY : IOB => IO_BUFF
!
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
!
USE MODD_BLD_DESCRIPTION_n, ONLY : BDD => BLD_DESC
USE MODD_DATA_BEM_n, ONLY : DTB => DATA_BEM
USE MODD_DATA_TEB_n, ONLY : DTT => DATA_TEB
!
USE MODD_SURF_ATM_GRID_n, ONLY : UG => SURF_ATM_GRID
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
USE MODD_DATA_COVER_PAR,  ONLY : JPCOVER
USE MODD_TEB_GRID_n, ONLY : TG => TEB_GRID
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS
USE MODD_BEM_OPTION_n, ONLY : BOP => BEM_OPTIONS
!
USE MODD_PREP,            ONLY : CINGRID_TYPE, CINTERP_TYPE, LINTERP
!
USE MODI_GET_LUOUT
USE MODI_ABOR1_SFX
USE MODI_OPEN_AUX_IO_SURF
USE MODI_GET_SURF_SIZE_n
USE MODI_PACK_PGD
USE MODI_PREP_GRID_EXTERN
USE MODI_PREP_OUTPUT_GRID
USE MODI_READ_SURF
USE MODI_READ_PGD_TEB_PAR_n
USE MODI_CLOSE_AUX_IO_SURF
USE MODI_CLEAN_PREP_OUTPUT_GRID
USE MODI_GOTO_WRAPPER_TEB_PATCH
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
!
IMPLICIT NONE
!
!*    0.1    Declaration of dummy arguments
!            ------------------------------
!
 CHARACTER(LEN=6),     INTENT(IN)  :: HPROGRAM    ! program calling
 CHARACTER(LEN=28),    INTENT(IN)  :: HINIFILE    ! file to read
 CHARACTER(LEN=6),     INTENT(IN)  :: HINIFILETYPE! file type
LOGICAL,              INTENT(IN)  :: OECOCLIMAP  ! flag to use ecoclimap
LOGICAL,              INTENT(IN)  :: OGARDEN     ! flag to use garden
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER :: IRESP   ! error return code
INTEGER :: ILUOUT  ! output listing logical unit
INTEGER :: INI     ! total 1D dimension (input grid)
INTEGER :: JLAYER  ! loop counter
INTEGER :: ILU     ! total 1D dimension (output grid, TOWN points only)
INTEGER :: JPATCH  ! TEB patch
REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER           :: IVERSION
INTEGER           :: IBUGFIX
!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ZOOM_PGD_TEB',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
TOP%LECOCLIMAP = OECOCLIMAP
TOP%LGARDEN = OGARDEN
!
IF (.NOT. OECOCLIMAP) THEN
  WRITE(ILUOUT,*) 'ERROR'
  WRITE(ILUOUT,*) 'Ecoclimap is not used'
  WRITE(ILUOUT,*) 'Routine zoom_pgd_teb.f90 must be updated'
  WRITE(ILUOUT,*) 'to interpolate all TEB physiographic fields'
  CALL ABOR1_SFX('ZOOM_PGD_TEB: ECOCLIMAP NOT USED, ROUTINE MUST BE UPDATED')
END IF
!
!
!*      1.     Preparation of IO for reading in the file
!              -----------------------------------------
!
!* Note that all points are read, even those without physical meaning.
!  These points will not be used during the horizontal interpolation step.
!  Their value must be defined as XUNDEF.
!
 CALL OPEN_AUX_IO_SURF(IOB, &
                       HINIFILE,HINIFILETYPE,'FULL  ')
!
 CALL GOTO_WRAPPER_TEB_PATCH(1)
!-------------------------------------------------------------------------------
!
!*    2.      Number of points and packing of general fields
!             ----------------------------------------------
!
!
 CALL GET_SURF_SIZE_n(DTCO, U, &
                      'TOWN  ',ILU)
!
ALLOCATE(TOP%LCOVER     (JPCOVER))
ALLOCATE(TOP%XZS        (ILU))
ALLOCATE(TG%XLAT       (ILU))
ALLOCATE(TG%XLON       (ILU))
ALLOCATE(TG%XMESH_SIZE (ILU))
!
 CALL PACK_PGD(HPROGRAM, 'TOWN  ',                      &
                TG%CGRID,  TG%XGRID_PAR,                     &
                TOP%LCOVER, TOP%XCOVER, TOP%XZS,                   &
                TG%XLAT, TG%XLON, TG%XMESH_SIZE                 )  
!
TG%NDIM = ILU
!
!
 CALL READ_SURF(IOB, &
                HPROGRAM,'VERSION',IVERSION,IRESP)
 CALL READ_SURF(IOB, &
                HPROGRAM,'BUG',IBUGFIX,IRESP)
!------------------------------------------------------------------------------
!
!*      3.     Reading of grid
!              ---------------
!
 CALL PREP_GRID_EXTERN(HINIFILETYPE,ILUOUT,CINGRID_TYPE,CINTERP_TYPE,INI)
!
 CALL PREP_OUTPUT_GRID(UG, U, &
                       ILUOUT,TG%CGRID,TG%XGRID_PAR,TG%XLAT,TG%XLON)
!
!
!------------------------------------------------------------------------------
!
!*      4.     Reading & interpolation of fields
!              ---------------------------------
!
!
IF (IVERSION<7 .OR. IVERSION==7 .AND. IBUGFIX<=2) THEN
  TOP%NTEB_PATCH=1
ELSE
  CALL READ_SURF(IOB, &
                HPROGRAM,'TEB_PATCH',TOP%NTEB_PATCH,IRESP)
END IF

!
 CALL READ_SURF(IOB, &
                HPROGRAM,'ROOF_LAYER',TOP%NROOF_LAYER,IRESP)
 CALL READ_SURF(IOB, &
                HPROGRAM,'ROAD_LAYER',TOP%NROAD_LAYER,IRESP)
 CALL READ_SURF(IOB, &
                HPROGRAM,'WALL_LAYER',TOP%NWALL_LAYER,IRESP)
!
IF (IVERSION<7 .OR.( IVERSION==7 .AND. IBUGFIX<=2)) THEN
  TOP%CBLD_ATYPE='ARI'
  TOP%CBEM = 'DEF'
ELSE
  CALL READ_SURF(IOB, &
                HPROGRAM,'BLD_ATYPE' ,TOP%CBLD_ATYPE,IRESP)
  CALL READ_SURF(IOB, &
                HPROGRAM,'BEM'       ,TOP%CBEM      ,IRESP)
END IF
!
IF (TOP%CBEM/='DEF') THEN
  CALL READ_SURF(IOB, &
                HPROGRAM,'FLOOR_LAYER',BOP%NFLOOR_LAYER,IRESP)
END IF
!
DO JPATCH=1,TOP%NTEB_PATCH
  CALL GOTO_WRAPPER_TEB_PATCH(JPATCH)
  CALL READ_PGD_TEB_PAR_n(BDD, DTB, DTT, TG, TOP, &
                          HPROGRAM,INI,'A')
!
!------------------------------------------------------------------------------
!
!*      5.     Gardens
!              -------
!
  IF (TOP%LGARDEN) CALL ZOOM_PGD_TEB_GARDEN
END DO
!
 CALL CLOSE_AUX_IO_SURF(HINIFILE,HINIFILETYPE)
!
 CALL CLEAN_PREP_OUTPUT_GRID
!
!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ZOOM_PGD_TEB',1,ZHOOK_HANDLE)
!------------------------------------------------------------------------------
!
CONTAINS
!
SUBROUTINE ZOOM_PGD_TEB_GARDEN
!
USE MODI_HOR_INTERPOL
!
USE MODD_TEB_VEG_n, ONLY : TVG => TEB_VEG_OPTIONS
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TGDO => TEB_GARDEN_OPTIONS
USE MODD_TEB_GARDEN_PGD_n, ONLY : TGDP => TEB_GARDEN_PGD
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:), POINTER     :: ZIN     ! field  on all surface points
!
REAL, DIMENSION(INI)              :: ZFIELD  ! field read
REAL, DIMENSION(ILU,1)            :: ZOUT    ! final field
REAL(KIND=JPRB) :: ZHOOK_HANDLE
 CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
!
IF (LHOOK) CALL DR_HOOK('ZOOM_PGD_TEB:ZOOM_PGD_TEB_GARDEN',0,ZHOOK_HANDLE)
!
LINTERP(:) = .TRUE.
!
IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
  CALL READ_SURF(IOB, &
                HPROGRAM,'GD_LAYER',TGDO%NGROUND_LAYER,IRESP)
  CALL READ_SURF(IOB, &
                HPROGRAM,'GD_ISBA',TVG%CISBA,IRESP)
  CALL READ_SURF(IOB, &
                HPROGRAM,'GD_PHOTO',TVG%CPHOTO,IRESP)
  CALL READ_SURF(IOB, &
                HPROGRAM,'GD_PEDOTF',TVG%CPEDOTF,IRESP)
  TVG%NNBIOMASS=1
  IF (TVG%CPHOTO=='NIT') TVG%NNBIOMASS=3  
ELSE
  CALL READ_SURF(IOB, &
                HPROGRAM,'TWN_LAYER',TGDO%NGROUND_LAYER,IRESP)
  CALL READ_SURF(IOB, &
                HPROGRAM,'TWN_ISBA',TVG%CISBA,IRESP)
  CALL READ_SURF(IOB, &
                HPROGRAM,'TWN_PHOTO',TVG%CPHOTO,IRESP)
  CALL READ_SURF(IOB, &
                HPROGRAM,'TWN_PEDOTF',TVG%CPEDOTF,IRESP)
  CALL READ_SURF(IOB, &
                HPROGRAM,'TWN_NBIOMASS',TVG%NNBIOMASS,IRESP)
ENDIF
!
!* sand
!
ALLOCATE(ZIN(INI,TGDO%NGROUND_LAYER))
YRECFM='TWN_SAND'
IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) YRECFM='GD_SAND'
 CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,ZFIELD,IRESP,HDIR='A')
DO JLAYER=1,TGDO%NGROUND_LAYER
  ZIN(:,JLAYER) = ZFIELD(:)
END DO
ALLOCATE(TGDP%XSAND(ILU,TGDO%NGROUND_LAYER))
 CALL HOR_INTERPOL(ILUOUT,ZIN,TGDP%XSAND)
DEALLOCATE(ZIN)
!
!* clay
!
ALLOCATE(ZIN(INI,TGDO%NGROUND_LAYER))
YRECFM='TWN_CLAY'
IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) YRECFM='GD_CLAY'
 CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,ZFIELD,IRESP,HDIR='A')
DO JLAYER=1,TGDO%NGROUND_LAYER
  ZIN(:,JLAYER) = ZFIELD(:)
END DO
ALLOCATE(TGDP%XCLAY(ILU,TGDO%NGROUND_LAYER))
 CALL HOR_INTERPOL(ILUOUT,ZIN,TGDP%XCLAY)
DEALLOCATE(ZIN)
!
!* runoff & drainage
!
ALLOCATE(ZIN(INI,1))
YRECFM='TWN_RUNOFFB'
IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) YRECFM='GD_RUNOFFB'
CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,ZFIELD,IRESP,HDIR='A')
ZIN(:,1) = ZFIELD(:)
ALLOCATE(TGDP%XRUNOFFB(ILU))
 CALL HOR_INTERPOL(ILUOUT,ZIN,ZOUT)
TGDP%XRUNOFFB(:) = ZOUT(:,1)
!
IF (IVERSION<=3) THEN
  TGDP%XWDRAIN = 0.
ELSE
 YRECFM='TWN_WDRAIN'
 IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) YRECFM='GD_WDRAIN'
 CALL READ_SURF(IOB, &
                HPROGRAM,YRECFM,ZFIELD,IRESP,HDIR='A')
 ZIN(:,1) = ZFIELD(:)
 ALLOCATE(TGDP%XWDRAIN(ILU))
 CALL HOR_INTERPOL(ILUOUT,ZIN,ZOUT)
 TGDP%XWDRAIN(:) = ZOUT(:,1)
ENDIF
!
DEALLOCATE(ZIN)
!
!* other garden parameters
!
 CALL READ_SURF(IOB, &
                HPROGRAM,'PAR_GARDEN',TGDO%LPAR_GARDEN,IRESP)
!
!!
IF (TGDO%LPAR_GARDEN) THEN
  WRITE(ILUOUT,*) 'ERROR'
  WRITE(ILUOUT,*) 'Specific garden fields are prescribed'
  WRITE(ILUOUT,*) 'Routine zoom_pgd_teb.f90 must be updated'
  WRITE(ILUOUT,*) 'to interpolate all TEB physiographic garden fields'
  CALL ABOR1_SFX('ZOOM_PGD_TEB: GARDEN fields used, ROUTINE MUST BE UPDATED')
END IF
!
IF (LHOOK) CALL DR_HOOK('ZOOM_PGD_TEB:ZOOM_PGD_TEB_GARDEN',1,ZHOOK_HANDLE)
!
END SUBROUTINE ZOOM_PGD_TEB_GARDEN
!_______________________________________________________________________________
!
END SUBROUTINE ZOOM_PGD_TEB
