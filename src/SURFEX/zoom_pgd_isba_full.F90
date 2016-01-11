!     ###########################################################
      SUBROUTINE ZOOM_PGD_ISBA_FULL (CHI, DTCO, DTI, IG, I, UG, U, &
                                     HPROGRAM,HINIFILE,HINIFILETYPE)
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
!
!
!
!
!
USE MODD_CH_ISBA_n, ONLY : CH_ISBA_t
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_GRID_n, ONLY : GRID_t
USE MODD_ISBA_n, ONLY : ISBA_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
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
USE MODI_ABOR1_SFX
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of dummy arguments
!            ------------------------------
!
!
TYPE(CH_ISBA_t), INTENT(INOUT) :: CHI
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTI
TYPE(GRID_t), INTENT(INOUT) :: IG
TYPE(ISBA_t), INTENT(INOUT) :: I
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
!
 CHARACTER(LEN=6),       INTENT(IN)  :: HPROGRAM     ! program calling
 CHARACTER(LEN=28),      INTENT(IN)  :: HINIFILE     ! input atmospheric file name
 CHARACTER(LEN=6),       INTENT(IN)  :: HINIFILETYPE ! input atmospheric file type
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER :: IVERSION, IBUGFIX
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
 CALL OPEN_AUX_IO_SURF(&
                       HINIFILE,HINIFILETYPE,'FULL  ')
!
 CALL READ_SURF(&
                HINIFILETYPE,'VERSION',IVERSION,IRESP)
 CALL READ_SURF(&
                HINIFILETYPE,'BUG',IBUGFIX,IRESP) 
!
!------------------------------------------------------------------------------
!
!*      2.     Reading of grid
!              ---------------
!
 CALL PREP_GRID_EXTERN(&
                       HINIFILETYPE,ILUOUT,CINGRID_TYPE,CINTERP_TYPE,INI)
!
 CALL PREP_OUTPUT_GRID(UG, U, &
                       ILUOUT,IG%CGRID,IG%XGRID_PAR,IG%XLAT,IG%XLON)
!
!------------------------------------------------------------------------------
!
!*      3.     Reading of fields
!              -----------------
!
!
ALLOCATE(ZFIELD(INI))
!
ALLOCATE(ZSAND(INI,I%O%NGROUND_LAYER))
 CALL READ_SURF(&
                HPROGRAM,'SAND',ZFIELD,IRESP,HDIR='A')
DO JLAYER=1,I%O%NGROUND_LAYER
  ZSAND(:,JLAYER) = ZFIELD(:)
END DO
!
ALLOCATE(ZCLAY(INI,I%O%NGROUND_LAYER))
 CALL READ_SURF(&
                HPROGRAM,'CLAY',ZFIELD,IRESP,HDIR='A')
DO JLAYER=1,I%O%NGROUND_LAYER
  ZCLAY(:,JLAYER) = ZFIELD(:)
END DO
!
!* Soil organic carbon profile
!
IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
   CALL READ_SURF(&
                HPROGRAM,'SOCP',I%O%LSOCP,IRESP)
ELSE
   I%O%LSOCP=.FALSE.
ENDIF
!
IF(I%O%LSOCP)THEN
!  
  ALLOCATE(I%P%XSOC (INI,I%O%NGROUND_LAYER))
!
  CALL READ_SURF(&
                HPROGRAM,'SOC_TOP',I%P%XSOC(:,1),IRESP)
  CALL READ_SURF(&
                HPROGRAM,'SOC_SUB',I%P%XSOC(:,2),IRESP)
!
  DO JLAYER=2,I%O%NGROUND_LAYER
    I%P%XSOC (:,JLAYER)=I%P%XSOC (:,2)
  END DO
!
ELSE
!  
  ALLOCATE(I%P%XSOC (0,1))
!
ENDIF
!
!* permafrost distribution
!
IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
   CALL READ_SURF(&
                HPROGRAM,'PERMAFROST',I%O%LPERM,IRESP)
ELSE
   I%O%LPERM=.FALSE.
ENDIF
!
IF(I%O%LPERM)THEN
!  
  ALLOCATE(I%P%XPERM (INI))
  CALL READ_SURF(&
                HPROGRAM,'PERM',I%P%XPERM(:),IRESP)
!
ELSE
!  
  ALLOCATE(I%P%XPERM (0))
!
ENDIF
!
!* groundwater distribution
!
IF (IVERSION>=8) THEN
   CALL READ_SURF(&
                HPROGRAM,'GWKEY',I%O%LGW,IRESP)
ELSE
   I%O%LGW=.FALSE.
ENDIF
!
IF(I%O%LGW)THEN
!  
  ALLOCATE(I%P%XGW (INI))
  CALL READ_SURF(&
                HPROGRAM,'GWFRAC',I%P%XGW(:),IRESP)
!
ELSE
!  
  ALLOCATE(I%P%XGW (0))
!
ENDIF
!
IF (IVERSION>7 .OR. (IVERSION==7 .AND. IBUGFIX>=3)) THEN
   CALL READ_SURF(&
                HPROGRAM,'NO',I%O%LNOF,IRESP)
ELSE
   I%O%LNOF = .FALSE.
ENDIF
!
!SOILNOX
!
IF (CHI%LCH_NO_FLUX) THEN
  !
  IF (I%O%LNOF) THEN
    !
    ALLOCATE(I%P%XPH(INI))
    CALL READ_SURF(&
                HPROGRAM,'PH',I%P%XPH(:),IRESP)
    !
    ALLOCATE(I%P%XFERT(INI))
    CALL READ_SURF(&
                HPROGRAM,'FERT',I%P%XFERT(:),IRESP)
    !
  ELSE
    CALL ABOR1_SFX("READ_PGD_ISBAn: WITH LCH_NO_FLUX=T, PH AND FERT FIELDS ARE GIVEN AT PGD STEP")
  ENDIF
  !
ELSE
  ALLOCATE(I%P%XPH (0))
  ALLOCATE(I%P%XFERT(0))
END IF
!
ALLOCATE(ZRUNOFFB(INI,1))
 CALL READ_SURF(&
                HPROGRAM,'RUNOFFB',ZFIELD,IRESP,HDIR='A')
ZRUNOFFB(:,1) = ZFIELD(:)
!
ALLOCATE(ZWDRAIN(INI,1))
 CALL READ_SURF(&
                HPROGRAM,'WDRAIN',ZFIELD,IRESP,HDIR='A')
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
 CALL HOR_INTERPOL(DTCO, U, &
                   ILUOUT,ZSAND,I%P%XSAND)
 CALL HOR_INTERPOL(DTCO, U, &
                   ILUOUT,ZCLAY,I%P%XCLAY)
ALLOCATE(ZOUTB(SIZE(I%P%XRUNOFFB),1))
 CALL HOR_INTERPOL(DTCO, U, &
                   ILUOUT,ZRUNOFFB,ZOUTB)
I%P%XRUNOFFB(:) = ZOUTB(:,1)
DEALLOCATE(ZOUTB)
ALLOCATE(ZOUTW(SIZE(I%P%XWDRAIN),1))
 CALL HOR_INTERPOL(DTCO, U, &
                   ILUOUT,ZWDRAIN,ZOUTW)
I%P%XWDRAIN(:) = ZOUTW(:,1)
DEALLOCATE(ZOUTW)
!
 CALL GET_TYPE_DIM_n(DTCO, U, &
                     'NATURE',IG%NDIM)
 CALL READ_PGD_ISBA_PAR_n(DTCO, U, &
                          DTI, IG, I%O, &
                          HPROGRAM,INI,.FALSE.,HDIR='A')
!
 CALL CLOSE_AUX_IO_SURF(HINIFILE,HINIFILETYPE)
!
 CALL CLEAN_PREP_OUTPUT_GRID
IF (LHOOK) CALL DR_HOOK('ZOOM_PGD_ISBA_FULL',1,ZHOOK_HANDLE)
!------------------------------------------------------------------------------
!
END SUBROUTINE ZOOM_PGD_ISBA_FULL
