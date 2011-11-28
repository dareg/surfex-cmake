!     #########
SUBROUTINE READ_SURF_ATM     (HPROGRAM, HFILE, PFIELD,                    &
                                KFORC_STEP, KNB, KRESP, KINIT               )  
!**************************************************************************
!
!!    PURPOSE
!!    -------
!         Read in the ascii file the atmospheric forcing for the actual time
!         step KFORC_STEP, and for the next one.
!         The two time step are needed for the time interpolation of the
!         forcing.
!         If the end of the file  is reached, set the two step to the last
!         values.
!         Return undef value if the variable is not present
!!
!!**  METHOD
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
!!
!!    AUTHOR
!!    ------
!!	A. Lemonsu  *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     03/2008
!          
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_IO_SURF_OL, ONLY : XSTART,XCOUNT,XSTRIDE,LPARTR
!
USE MODD_ARCH, ONLY : LITTLE_ENDIAN_ARCH
!
USE MODE_CHAR2REAL
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!

! global variables
REAL, DIMENSION(:,:),INTENT(INOUT) :: PFIELD
INTEGER,INTENT(IN)               :: KFORC_STEP
INTEGER,INTENT(IN)               :: KNB 
INTEGER,INTENT(IN)               :: KRESP
INTEGER,INTENT(IN)               :: KINIT
CHARACTER(LEN=6)    ,INTENT(IN)  :: HPROGRAM
CHARACTER(LEN=15)   ,INTENT(IN)  :: HFILE

! local variables
INTEGER                          :: ILUOUT
INTEGER                          :: ICOUNT
INTEGER                          :: I,INI,INB
INTEGER                          :: IRET,INB_FORC,INPATCH
CHARACTER(LEN=4), DIMENSION(:), ALLOCATABLE  :: ZF
REAL*4, DIMENSION(:), ALLOCATABLE :: ZFIELD
LOGICAL                          :: GSWAP              ! T: swap has been done
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('READ_SURF_ATM',0,ZHOOK_HANDLE)
INI = SIZE(PFIELD,2)
ALLOCATE(ZF(INI))
ALLOCATE(ZFIELD(INI))
!
IF (HPROGRAM == 'ASCII ') THEN
!
 IF (KFORC_STEP .EQ. 1) THEN
  REWIND(KINIT)
  DO I=1,KNB
   READ(UNIT=KINIT,FMT='(50(F20.5))') PFIELD(I,:)
  ENDDO
 ELSE
  PFIELD(1,:) = PFIELD(KNB,:)
  DO I=2,KNB   
    READ(UNIT=KINIT,FMT='(50(F20.5))') PFIELD(I,:)
  ENDDO
 ENDIF
!
ELSE IF (HPROGRAM == 'BINARY') THEN
!
 IF (KFORC_STEP .EQ. 1) THEN
  GSWAP = .FALSE.
  I=1
  DO I=1,KNB
    IF (LITTLE_ENDIAN_ARCH) THEN
      READ(UNIT=KINIT,REC=KFORC_STEP+I-1) ZF(:)
      PFIELD(I,:) = ZF(:)
    ELSE
      READ(UNIT=KINIT,REC=KFORC_STEP+I-1) ZFIELD(:)
      PFIELD(I,:)=ZFIELD(:)
    ENDIF
    IF (      ANY(ABS(PFIELD(I,:))>0. .AND. ABS(PFIELD(I,:))<1.E-30) &
        .OR. ANY(ABS(PFIELD(I,:))>1.E6)                       ) THEN  
      CALL ABOR1_SFX('READ_SURF_ATM: SWAP SET IN YOUR PARAMS_CONFIG FILE SEEMS '//&
        'INAPPROPRIATE - VERIFY ')  
    END IF
  ENDDO
 ELSE
  PFIELD(1,:) = PFIELD(KNB,:)
  DO I=2,KNB
    IF (LITTLE_ENDIAN_ARCH) THEN
      READ(UNIT=KINIT,REC=KFORC_STEP+I-1) ZF(:)
      PFIELD(I,:) = ZF(:)
    ELSE
      READ(UNIT=KINIT,REC=KFORC_STEP+I-1) ZFIELD(:)
      PFIELD(I,:)=ZFIELD(:)
    ENDIF
    IF (      ANY(ABS(PFIELD(I,:))>0. .AND. ABS(PFIELD(I,:))<1.E-30) &
        .OR. ANY(ABS(PFIELD(I,:))>1.E6)                       ) THEN  
      CALL ABOR1_SFX('READ_SURF_ATM: SWAP SET IN YOUR PARAMS_CONFIG FILE SEEMS '//&
        'INAPPROPRIATE - VERIFY  ')  
     END IF
   ENDDO
 ENDIF
!
ENDIF

DEALLOCATE(ZF)
DEALLOCATE(ZFIELD)

LPARTR=.FALSE.
IF (LHOOK) CALL DR_HOOK('READ_SURF_ATM',1,ZHOOK_HANDLE)

END SUBROUTINE READ_SURF_ATM
