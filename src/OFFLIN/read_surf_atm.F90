!     #########
SUBROUTINE READ_SURF_ATM     (HPROGRAM, PFIELD, KFORC_STEP, KNB, KINIT)  
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
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO, NINDEX, XTIME_COMM_READ, XTIME_NPIO_READ
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_IO_SURF_OL, ONLY : XSTART,XCOUNT,XSTRIDE,LPARTR
USE MODD_IO_SURF_ASC,ONLY : NNI_FORC
!
USE MODD_ARCH, ONLY : LITTLE_ENDIAN_ARCH
!
USE MODE_CHAR2REAL
!
USE MODI_ABOR1_SFX
USE MODI_READ_AND_SEND_MPI
USE MODI_GATHER_AND_WRITE_MPI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
#ifndef NOMPI
INCLUDE 'mpif.h'
#endif
!
! global variables
REAL, DIMENSION(:,:),INTENT(INOUT) :: PFIELD
INTEGER,INTENT(IN)               :: KFORC_STEP
INTEGER,INTENT(IN)               :: KNB 
INTEGER,INTENT(IN)               :: KINIT
CHARACTER(LEN=6)    ,INTENT(IN)  :: HPROGRAM

! local variables
INTEGER                          :: I,INI, J
CHARACTER(LEN=4), DIMENSION(:), ALLOCATABLE  :: YF
CHARACTER(LEN=4) :: YWORK
DOUBLE PRECISION :: XTIME0
REAL*4, DIMENSION(:), ALLOCATABLE :: ZFIELD4
REAL*4                            :: ZWORK4
REAL, DIMENSION(:,:), ALLOCATABLE :: ZFIELD
REAL                              :: ZWORK
LOGICAL                          :: GSWAP              ! T: swap has been done
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('READ_SURF_ATM',0,ZHOOK_HANDLE)
!
IF (NRANK==NPIO) THEN
  INI = SIZE(NINDEX)
  ALLOCATE(ZFIELD(INI,SIZE(PFIELD,2)))
  IF (HPROGRAM == 'BINARY') THEN
    ALLOCATE(ZFIELD4(INI))          
    ALLOCATE(YF(INI))
  ENDIF
ELSE
  ALLOCATE(ZFIELD(0,0)) 
  ALLOCATE(ZFIELD4(0))
  ALLOCATE(YF(0))
ENDIF
!
CALL GATHER_AND_WRITE_MPI(PFIELD,ZFIELD)
!
IF (NRANK==NPIO) THEN
  !
#ifndef NOMPI
  XTIME0 = MPI_WTIME() 
#endif 
  !  
  IF (HPROGRAM == 'ASCII ') THEN
    !
    IF (KFORC_STEP .EQ. 1) THEN
      REWIND(KINIT)
      DO I=1,KNB
        IF (NNI_FORC==1) THEN
          READ(UNIT=KINIT,FMT='(F20.5)') ZWORK
          ZFIELD(:,I) = ZWORK
        ELSE
          READ(UNIT=KINIT,FMT='(50(F20.5))') ZFIELD(:,I)
        END IF
      ENDDO
    ELSE
      ZFIELD(:,1) = ZFIELD(:,KNB)
      DO I=2,KNB
        IF (NNI_FORC==1) THEN
          READ(UNIT=KINIT,FMT='(F20.5)') ZWORK
          ZFIELD(:,I) = ZWORK
        ELSE
          READ(UNIT=KINIT,FMT='(50(F20.5))') ZFIELD(:,I)
        END IF
      ENDDO
    ENDIF
    !
  ELSE IF (HPROGRAM == 'BINARY') THEN
    !
    IF (KFORC_STEP .EQ. 1) THEN
      GSWAP = .FALSE.
      DO I=1,KNB
        IF (LITTLE_ENDIAN_ARCH) THEN
          IF (NNI_FORC==1) THEN
            READ(UNIT=KINIT,REC=KFORC_STEP+I-1) YWORK
            YF(:) = YWORK
          ELSE
            READ(UNIT=KINIT,REC=KFORC_STEP+I-1) YF(:)
          END IF
          ZFIELD(:,I) = YF(:)
        ELSE
          IF (NNI_FORC==1) THEN
            READ(UNIT=KINIT,REC=KFORC_STEP+I-1) ZWORK4
            ZFIELD4(:) = ZWORK4
          ELSE
            READ(UNIT=KINIT,REC=KFORC_STEP+I-1) ZFIELD4(:)
          END IF
          ZFIELD(:,I) = ZFIELD4(:)
        ENDIF
        IF (     ANY(ABS(ZFIELD(:,I))>0. .AND. ABS(ZFIELD(:,I))<1.E-30) &
            .OR. ANY(ABS(ZFIELD(:,I))>1.E6)                       ) THEN  
          CALL ABOR1_SFX('READ_SURF_ATM: SWAP SET IN YOUR PARAMS_CONFIG FILE SEEMS '//&
            'INAPPROPRIATE - VERIFY ')  
        END IF
      ENDDO
    ELSE
      ZFIELD(:,1) = ZFIELD(:,KNB)
      DO I=2,KNB
        IF (LITTLE_ENDIAN_ARCH) THEN
          IF (NNI_FORC==1) THEN
            READ(UNIT=KINIT,REC=KFORC_STEP+I-1) YWORK
            YF(:) = YWORK
          ELSE
            READ(UNIT=KINIT,REC=KFORC_STEP+I-1) YF(:)
          END IF
          ZFIELD(:,I) = YF(:)
        ELSE
          IF (NNI_FORC==1) THEN
            READ(UNIT=KINIT,REC=KFORC_STEP+I-1) ZWORK4
            ZFIELD4(:) = ZWORK4
          ELSE
            READ(UNIT=KINIT,REC=KFORC_STEP+I-1) ZFIELD4(:)
          END IF                
          ZFIELD(:,I) = ZFIELD4(:)
        ENDIF
        IF (     ANY(ABS(ZFIELD(:,I))>0. .AND. ABS(ZFIELD(:,I))<1.E-30) &
            .OR. ANY(ABS(ZFIELD(:,I))>1.E6)                       ) THEN  
          CALL ABOR1_SFX('READ_SURF_ATM: SWAP SET IN YOUR PARAMS_CONFIG FILE SEEMS '//&
            'INAPPROPRIATE - VERIFY  ')  
        END IF
      ENDDO
    ENDIF
    !
  ENDIF
  !
#ifndef NOMPI
  XTIME_NPIO_READ = XTIME_NPIO_READ + (MPI_WTIME() - XTIME0)
#endif
  !  
ENDIF
!
CALL READ_AND_SEND_MPI(ZFIELD,PFIELD)
!
DEALLOCATE(ZFIELD)
IF (HPROGRAM=='BINARY') THEN
  DEALLOCATE(YF)
  DEALLOCATE(ZFIELD4)
ENDIF
!
LPARTR=.FALSE.
!
IF (LHOOK) CALL DR_HOOK('READ_SURF_ATM',1,ZHOOK_HANDLE)

END SUBROUTINE READ_SURF_ATM
