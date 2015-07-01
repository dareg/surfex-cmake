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
USE MODD_IO_BUFF_n, ONLY : IOB => IO_BUFF
!
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
!
USE MODD_DATA_ISBA_n, ONLY : DTI => DATA_ISBA
!
USE MODD_SURF_ATM_GRID_n, ONLY : UG => SURF_ATM_GRID
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
USE MODD_ISBA_n, ONLY : I => ISBA
USE MODD_CH_ISBA_n, ONLY : CHI => CH_ISBA                                  
USE MODD_ISBA_GRID_n, ONLY : IG => ISBA_GRID
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
 CALL OPEN_AUX_IO_SURF(IOB, &
                       HINIFILE,HINIFILETYPE,'FULL  ')
!
 CALL READ_SURF(IOB, &
                HINIFILETYPE,'VERSION',IVERSION,IRESP)
 CALL READ_SURF(IOB, &
                HINIFILETYPE,'BUG',IBUGFIX,IRESP) 
!
!------------------------------------------------------------------------------
!
!*      2.     Reading of grid
!              ---------------
!
 CALL PREP_GRID_EXTERN(IOB, &
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
ALLOCATE(ZSAND(INI,I%NGROUND_LAYER))
 CALL READ_SURF(IOB, &
                HPROGRAM,'SAND',ZFIELD,IRESP,HDIR='A')
DO JLAYER=1,I%NGROUND_LAYER
  ZSAND(:,JLAYER) = ZFIELD(:)
END DO
!
ALLOCATE(ZCLAY(INI,I%NGROUND_LAYER))
 CALL READ_SURF(IOB, &
                HPROGRAM,'CLAY',ZFIELD,IRESP,HDIR='A')
DO JLAYER=1,I%NGROUND_LAYER
  ZCLAY(:,JLAYER) = ZFIELD(:)
END DO
!
!* Soil organic carbon profile
!
IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
   CALL READ_SURF(IOB, &
                HPROGRAM,'SOCP',I%LSOCP,IRESP)
ELSE
   I%LSOCP=.FALSE.
ENDIF
!
IF(I%LSOCP)THEN
!  
  ALLOCATE(I%XSOC (INI,I%NGROUND_LAYER))
!
  CALL READ_SURF(IOB, &
                HPROGRAM,'SOC_TOP',I%XSOC(:,1),IRESP)
  CALL READ_SURF(IOB, &
                HPROGRAM,'SOC_SUB',I%XSOC(:,2),IRESP)
!
  DO JLAYER=2,I%NGROUND_LAYER
    I%XSOC (:,JLAYER)=I%XSOC (:,2)
  END DO
!
ELSE
!  
  ALLOCATE(I%XSOC (0,1))
!
ENDIF
!
!* permafrost distribution
!
IF (IVERSION>7 .OR. IVERSION==7 .AND. IBUGFIX>=3) THEN
   CALL READ_SURF(IOB, &
                HPROGRAM,'PERMAFROST',I%LPERM,IRESP)
ELSE
   I%LPERM=.FALSE.
ENDIF
!
IF(I%LPERM)THEN
!  
  ALLOCATE(I%XPERM (INI))
  CALL READ_SURF(IOB, &
                HPROGRAM,'PERM',I%XPERM(:),IRESP)
!
ELSE
!  
  ALLOCATE(I%XPERM (0))
!
ENDIF
!
!* groundwater distribution
!
IF (IVERSION>=8) THEN
   CALL READ_SURF(IOB, &
                HPROGRAM,'GWKEY',I%LGW,IRESP)
ELSE
   I%LGW=.FALSE.
ENDIF
!
IF(I%LGW)THEN
!  
  ALLOCATE(I%XGW (INI))
  CALL READ_SURF(IOB, &
                HPROGRAM,'GWFRAC',I%XGW(:),IRESP)
!
ELSE
!  
  ALLOCATE(I%XGW (0))
!
ENDIF
!
IF (IVERSION>7 .OR. (IVERSION==7 .AND. IBUGFIX>=3)) THEN
   CALL READ_SURF(IOB, &
                HPROGRAM,'NO',I%LNOF,IRESP)
ELSE
   I%LNOF = .FALSE.
ENDIF
!
!SOILNOX
!
IF (CHI%LCH_NO_FLUX) THEN
  !
  IF (I%LNOF) THEN
    !
    ALLOCATE(I%XPH(INI))
    CALL READ_SURF(IOB, &
                HPROGRAM,'PH',I%XPH(:),IRESP)
    !
    ALLOCATE(I%XFERT(INI))
    CALL READ_SURF(IOB, &
                HPROGRAM,'FERT',I%XFERT(:),IRESP)
    !
  ELSE
    CALL ABOR1_SFX("READ_PGD_ISBAn: WITH LCH_NO_FLUX=T, PH AND FERT FIELDS ARE GIVEN AT PGD STEP")
  ENDIF
  !
ELSE
  ALLOCATE(I%XPH (0))
  ALLOCATE(I%XFERT(0))
END IF
!
ALLOCATE(ZRUNOFFB(INI,1))
 CALL READ_SURF(IOB, &
                HPROGRAM,'RUNOFFB',ZFIELD,IRESP,HDIR='A')
ZRUNOFFB(:,1) = ZFIELD(:)
!
ALLOCATE(ZWDRAIN(INI,1))
 CALL READ_SURF(IOB, &
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
                   ILUOUT,ZSAND,I%XSAND)
 CALL HOR_INTERPOL(DTCO, U, &
                   ILUOUT,ZCLAY,I%XCLAY)
ALLOCATE(ZOUTB(SIZE(I%XRUNOFFB),1))
 CALL HOR_INTERPOL(DTCO, U, &
                   ILUOUT,ZRUNOFFB,ZOUTB)
I%XRUNOFFB(:) = ZOUTB(:,1)
DEALLOCATE(ZOUTB)
ALLOCATE(ZOUTW(SIZE(I%XWDRAIN),1))
 CALL HOR_INTERPOL(DTCO, U, &
                   ILUOUT,ZWDRAIN,ZOUTW)
I%XWDRAIN(:) = ZOUTW(:,1)
DEALLOCATE(ZOUTW)
!
 CALL GET_TYPE_DIM_n(DTCO, U, &
                     'NATURE',IG%NDIM)
 CALL READ_PGD_ISBA_PAR_n(DTI, IG, I, &
                          HPROGRAM,INI,.FALSE.,HDIR='A')
!
 CALL CLOSE_AUX_IO_SURF(HINIFILE,HINIFILETYPE)
!
 CALL CLEAN_PREP_OUTPUT_GRID
IF (LHOOK) CALL DR_HOOK('ZOOM_PGD_ISBA_FULL',1,ZHOOK_HANDLE)
!------------------------------------------------------------------------------
!
END SUBROUTINE ZOOM_PGD_ISBA_FULL
