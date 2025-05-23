!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
MODULE MODI_GATHER_AND_WRITE_MPI
!----------------------------------------------------
!!    MODIFICATIONS
!!    -------------
!!      Original       
!!      J.Escobar      10/06/2013: replace DOUBLE PRECISION by REAL to handle problem for promotion of real on IBM SP
!----------------------------------------------------
!
INTERFACE GATHER_AND_WRITE_MPI
!
SUBROUTINE GATHER_AND_WRITE_MPI_N1D(KWORK,KWORK2,KMASK,KMISS)
!
INTEGER, DIMENSION(:), INTENT(IN) :: KWORK
INTEGER, DIMENSION(:), INTENT(OUT) :: KWORK2
!
INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: KMASK
INTEGER, OPTIONAL, INTENT(IN) :: KMISS
!
END SUBROUTINE GATHER_AND_WRITE_MPI_N1D
!
SUBROUTINE GATHER_AND_WRITE_MPI_N2D(KWORK,KWORK2,KMASK,KMISS)
!
INTEGER, DIMENSION(:,:), INTENT(IN) :: KWORK
INTEGER, DIMENSION(:,:), INTENT(OUT) :: KWORK2
!
INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: KMASK
INTEGER, OPTIONAL, INTENT(IN) :: KMISS
!
END SUBROUTINE GATHER_AND_WRITE_MPI_N2D
!
SUBROUTINE GATHER_AND_WRITE_MPI_N3D(KWORK,KWORK2,KMASK,KMISS)
!
INTEGER, DIMENSION(:,:,:), INTENT(IN) :: KWORK
INTEGER, DIMENSION(:,:,:), INTENT(OUT) :: KWORK2
!
INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: KMASK
INTEGER, OPTIONAL, INTENT(IN) :: KMISS
!
END SUBROUTINE GATHER_AND_WRITE_MPI_N3D
!
SUBROUTINE GATHER_AND_WRITE_MPI_X1D(PWORK,PWORK2,KMASK,PMISS)
!
REAL, DIMENSION(:), INTENT(IN) :: PWORK
REAL(KIND=8), DIMENSION(:), INTENT(OUT) :: PWORK2
!
INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: KMASK
REAL, OPTIONAL, INTENT(IN) :: PMISS
!
END SUBROUTINE GATHER_AND_WRITE_MPI_X1D
!
SUBROUTINE GATHER_AND_WRITE_MPI_X2D(PWORK,PWORK2,KMASK,PMISS)
!
REAL, DIMENSION(:,:), INTENT(IN) :: PWORK
REAL(KIND=8), DIMENSION(:,:), INTENT(OUT) :: PWORK2
!
INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: KMASK
REAL, OPTIONAL, INTENT(IN) :: PMISS
!
END SUBROUTINE GATHER_AND_WRITE_MPI_X2D
!
SUBROUTINE GATHER_AND_WRITE_MPI_X3D(PWORK,PWORK2,KMASK,PMISS)
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PWORK
REAL(KIND=8), DIMENSION(:,:,:), INTENT(OUT) :: PWORK2
!
INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: KMASK
REAL, OPTIONAL, INTENT(IN) :: PMISS
!
END SUBROUTINE GATHER_AND_WRITE_MPI_X3D
!
SUBROUTINE GATHER_AND_WRITE_MPI_X1DK4(PWORK,PWORK2,KMASK,PMISS)
!
REAL, DIMENSION(:), INTENT(IN) :: PWORK
REAL(KIND=4), DIMENSION(:), INTENT(OUT) :: PWORK2
!
INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: KMASK
REAL, OPTIONAL, INTENT(IN) :: PMISS
!
END SUBROUTINE GATHER_AND_WRITE_MPI_X1DK4
!
SUBROUTINE GATHER_AND_WRITE_MPI_X2DK4(PWORK,PWORK2,KMASK,PMISS)
!
REAL, DIMENSION(:,:), INTENT(IN) :: PWORK
REAL(KIND=4), DIMENSION(:,:), INTENT(OUT) :: PWORK2
!
INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: KMASK
REAL, OPTIONAL, INTENT(IN) :: PMISS
!
END SUBROUTINE GATHER_AND_WRITE_MPI_X2DK4
!
SUBROUTINE GATHER_AND_WRITE_MPI_X3DK4(PWORK,PWORK2,KMASK,PMISS)
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PWORK
REAL(KIND=4), DIMENSION(:,:,:), INTENT(OUT) :: PWORK2
!
INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: KMASK
REAL, OPTIONAL, INTENT(IN) :: PMISS
!
END SUBROUTINE GATHER_AND_WRITE_MPI_X3DK4
!
END INTERFACE
!
END MODULE MODI_GATHER_AND_WRITE_MPI
!
SUBROUTINE GATHER_AND_WRITE_MPI_N1D(KWORK,KWORK2,KMASK,KMISS)
!
USE MODD_SURFEX_MPI, ONLY : NINDEX, NPROC, NRANK, NCOMM, NPIO, NSIZE, &
                            XTIME_CALC_WRITE, XTIME_COMM_WRITE, &
                            IDX_W, WLOG_MPI, LSFX_MPI
!
USE MODI_UNPACK_SAME_RANK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
#ifdef SFX_MPI
INCLUDE "mpif.h"
#endif
!
INTEGER, DIMENSION(:), INTENT(IN) :: KWORK
INTEGER, DIMENSION(:), INTENT(OUT) :: KWORK2
!
INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: KMASK
INTEGER, INTENT(IN), OPTIONAL :: KMISS
!
INTEGER, DIMENSION(NSIZE) :: IINTER
INTEGER, DIMENSION(NSIZE) :: IWORK
REAL   :: XTIME0
!
#ifdef SFX_MPI
INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ISTATUS
#endif
INTEGER :: ICPT
INTEGER :: I,J, IP1, IS1
INTEGER :: INFOMPI
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE, ZHOOK_HANDLE_OMP
!
IF (LHOOK) CALL DR_HOOK('GATHER_AND_WRITE_MPI_N1D',0,ZHOOK_HANDLE)
!
IWORK(:) = 0
!
#ifdef SFX_MPI
IF (LSFX_MPI) XTIME0 = MPI_WTIME()
#endif
!
IF (PRESENT(KMASK)) THEN
  IF (PRESENT(KMISS)) THEN
    CALL UNPACK_SAME_RANK(KMASK,KWORK,IWORK(:),KMISS)
  ELSE
    CALL UNPACK_SAME_RANK(KMASK,KWORK,IWORK(:))
  END IF
ELSE
  IWORK(1:SIZE(KWORK)) = KWORK(:)
ENDIF
!
#ifdef SFX_MPI
IF (LSFX_MPI) THEN
  XTIME_CALC_WRITE = XTIME_CALC_WRITE + (MPI_WTIME() - XTIME0)
  XTIME0 = MPI_WTIME()
ENDIF
#endif
!
IF (NRANK/=NPIO) THEN
  !
  IDX_W = IDX_W + 1
  !  
#ifdef SFX_MPI
  IF (LSFX_MPI) THEN
    XTIME0 = MPI_WTIME()
    CALL MPI_SEND(IWORK,SIZE(IWORK)*KIND(IWORK)/4,MPI_INTEGER,NPIO,IDX_W,NCOMM,INFOMPI)
    XTIME_COMM_WRITE = XTIME_COMM_WRITE + (MPI_WTIME() - XTIME0)
  ENDIF
#endif
  !
ELSE
  !
  KWORK2(:) = 0        
  !
  IDX_W = IDX_W + 1 
  !
  DO I=0,NPROC-1
    !
#ifdef SFX_MPI
    IF (LSFX_MPI) XTIME0 = MPI_WTIME()
#endif
    !
    IF (I/=NPIO) THEN
#ifdef SFX_MPI
      IF (LSFX_MPI) CALL MPI_RECV(IINTER,SIZE(IINTER)*KIND(IINTER)/4,MPI_INTEGER,I,IDX_W,NCOMM,ISTATUS,INFOMPI)
#endif
    ELSE
      IINTER(:) = IWORK(:)
    ENDIF
    !
#ifdef SFX_MPI 
    IF (LSFX_MPI) THEN
      XTIME_COMM_WRITE = XTIME_COMM_WRITE + (MPI_WTIME() - XTIME0)
      XTIME0 = MPI_WTIME()
    ENDIF
#endif
    !
    ICPT = 0
    !    
    DO J=1,SIZE(NINDEX)
      !
      IF ( NINDEX(J)==I ) THEN
        ICPT = ICPT + 1
        KWORK2(J) = IINTER(ICPT)
      ENDIF
      !
    ENDDO
    !
#ifdef SFX_MPI
    IF (LSFX_MPI) XTIME_CALC_WRITE = XTIME_CALC_WRITE + (MPI_WTIME() - XTIME0)
#endif
    !
  ENDDO
  !
ENDIF
!
IF (LHOOK) CALL DR_HOOK('GATHER_AND_WRITE_MPI_N1D',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE GATHER_AND_WRITE_MPI_N1D
!
!**************************************************************************
!
SUBROUTINE GATHER_AND_WRITE_MPI_N2D(KWORK,KWORK2,KMASK,KMISS)
!
USE MODD_SURFEX_MPI, ONLY : NINDEX, NPROC, NRANK, NCOMM, NPIO, NSIZE, &
                            XTIME_CALC_WRITE, XTIME_COMM_WRITE, &
                            IDX_W, WLOG_MPI, LSFX_MPI
!
USE MODI_UNPACK_SAME_RANK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
#ifdef SFX_MPI
INCLUDE "mpif.h"
#endif
!
INTEGER, DIMENSION(:,:), INTENT(IN) :: KWORK
INTEGER, DIMENSION(:,:), INTENT(OUT) :: KWORK2
!
INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: KMASK
INTEGER, INTENT(IN), OPTIONAL :: KMISS
!
INTEGER, DIMENSION(NSIZE,SIZE(KWORK2,2)) :: IINTER
INTEGER, DIMENSION(NSIZE,SIZE(KWORK2,2)) :: IWORK
REAL   :: XTIME0
!
#ifdef SFX_MPI
INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ISTATUS
#endif
INTEGER :: ICPT
INTEGER :: I,J
INTEGER :: INFOMPI
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('GATHER_AND_WRITE_MPI_N2D',0,ZHOOK_HANDLE)
!
IWORK(:,:) = 0
!
#ifdef SFX_MPI
IF (LSFX_MPI) XTIME0 = MPI_WTIME()
#endif
!
IF (SIZE(KWORK,1)>0) THEN
  IF (PRESENT(KMASK)) THEN
    IF (PRESENT(KMISS)) THEN
      CALL UNPACK_SAME_RANK(KMASK,KWORK,IWORK(:,:),KMISS)
    ELSE
      CALL UNPACK_SAME_RANK(KMASK,KWORK,IWORK(:,:))
    END IF
  ELSE
    IWORK(1:SIZE(KWORK,1),:) = KWORK(:,:)
  ENDIF
ENDIF
!
#ifdef SFX_MPI
IF (LSFX_MPI) THEN
  XTIME_CALC_WRITE = XTIME_CALC_WRITE + (MPI_WTIME() - XTIME0)
  XTIME0 = MPI_WTIME()
ENDIF
#endif
!
IF (NRANK/=NPIO) THEN
  !
  IDX_W = IDX_W + 1
  !  
#ifdef SFX_MPI
  IF (LSFX_MPI) THEN
    XTIME0 = MPI_WTIME()
    CALL MPI_SEND(IWORK(:,:),SIZE(IWORK)*KIND(IWORK)/4,MPI_INTEGER,NPIO,IDX_W,NCOMM,INFOMPI)
    XTIME_COMM_WRITE = XTIME_COMM_WRITE + (MPI_WTIME() - XTIME0)
  ENDIF
#endif
  !
ELSE
  !
  IDX_W = IDX_W + 1 
  !
  DO I=1,NPROC
    !
#ifdef SFX_MPI
    IF (LSFX_MPI) XTIME0 = MPI_WTIME()
#endif
    !
    IF (I<NPROC) THEN
#ifdef SFX_MPI
      IF (LSFX_MPI) CALL MPI_RECV(IINTER,SIZE(IINTER)*KIND(IINTER)/4,MPI_INTEGER,I,IDX_W,NCOMM,ISTATUS,INFOMPI)
#endif
    ELSE
      IINTER(:,:) = IWORK(:,:)
    ENDIF
    !
#ifdef SFX_MPI
    IF (LSFX_MPI) THEN
      XTIME_COMM_WRITE = XTIME_COMM_WRITE + (MPI_WTIME() - XTIME0)
      XTIME0 = MPI_WTIME()
    ENDIF
#endif
    !
    ICPT = 0
    !  
    DO J=1,SIZE(NINDEX)
      !
      IF ( NINDEX(J)==MOD(I,NPROC) ) THEN
        ICPT = ICPT + 1
        KWORK2(J,:) = IINTER(ICPT,:)
      ENDIF
      !
    ENDDO
    !
#ifdef SFX_MPI
    IF (LSFX_MPI) XTIME_CALC_WRITE = XTIME_CALC_WRITE + (MPI_WTIME() - XTIME0)
#endif
    !
  ENDDO
  !
ENDIF
!
IF (LHOOK) CALL DR_HOOK('GATHER_AND_WRITE_MPI_N2D',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE GATHER_AND_WRITE_MPI_N2D
!
!**************************************************************************
!
SUBROUTINE GATHER_AND_WRITE_MPI_N3D(KWORK,KWORK2,KMASK,KMISS)
!
USE MODD_SURFEX_MPI, ONLY : NINDEX, NPROC, NRANK, NCOMM, NPIO, NSIZE, &
                            XTIME_CALC_WRITE, XTIME_COMM_WRITE, &
                            IDX_W, WLOG_MPI, LSFX_MPI
!
USE MODI_UNPACK_SAME_RANK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
#ifdef SFX_MPI
INCLUDE "mpif.h"
#endif
!
INTEGER, DIMENSION(:,:,:), INTENT(IN) :: KWORK
INTEGER, DIMENSION(:,:,:), INTENT(OUT) :: KWORK2
!
INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: KMASK
INTEGER, INTENT(IN), OPTIONAL :: KMISS
!
INTEGER, DIMENSION(NSIZE,SIZE(KWORK2,2),SIZE(KWORK2,3)) :: IINTER
INTEGER, DIMENSION(NSIZE,SIZE(KWORK,2),SIZE(KWORK,3)) :: IWORK
!
DOUBLE PRECISION   :: XTIME0
!
#ifdef SFX_MPI
INTEGER, DIMENSION(MPI_STATUS_SIZE) :: ISTATUS
#endif
INTEGER :: ICPT
INTEGER :: I,J
INTEGER :: INFOMPI
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('GATHER_AND_WRITE_MPI_N3D',0,ZHOOK_HANDLE)
!
IWORK(:,:,:) = 0
!
#ifdef SFX_MPI
IF (LSFX_MPI) XTIME0 = MPI_WTIME()
#endif
!
IF (SIZE(KWORK,1)>0) THEN
  IF (PRESENT(KMASK)) THEN
    IF (PRESENT(KMISS)) THEN
      CALL UNPACK_SAME_RANK(KMASK,KWORK(:,:,:),IWORK(:,:,:),KMISS)
    ELSE
      CALL UNPACK_SAME_RANK(KMASK,KWORK(:,:,:),IWORK(:,:,:))
    END IF
  ELSE
    IWORK(1:SIZE(KWORK,1),:,:) = KWORK(:,:,:)
  ENDIF
ENDIF
!
#ifdef SFX_MPI
IF (LSFX_MPI) THEN
  XTIME_CALC_WRITE = XTIME_CALC_WRITE + (MPI_WTIME() - XTIME0)
  XTIME0 = MPI_WTIME()
ENDIF
#endif
!
IF (NRANK/=NPIO) THEN
  !
  IDX_W = IDX_W + 1
  !  
#ifdef SFX_MPI
  IF (LSFX_MPI) THEN
    XTIME0 = MPI_WTIME()
    CALL MPI_SEND(IWORK(:,:,:),SIZE(IWORK)*KIND(IWORK)/4,MPI_INTEGER,NPIO,IDX_W,NCOMM,INFOMPI)
    XTIME_COMM_WRITE = XTIME_COMM_WRITE + (MPI_WTIME() - XTIME0)
  ENDIF
#endif
  !
ELSE
  ! 
  KWORK2(:,:,:) = 0
  !
  IDX_W = IDX_W + 1 
  !
  DO I=0,NPROC-1
    !
#ifdef SFX_MPI    
    IF (LSFX_MPI) XTIME0 = MPI_WTIME()
#endif    
    !
    IF (I/=NPIO) THEN
#ifdef SFX_MPI
      CALL MPI_RECV(IINTER,SIZE(IINTER)*KIND(IINTER)/4,MPI_INTEGER,I,IDX_W,NCOMM,ISTATUS,INFOMPI)
#endif
    ELSE
      IINTER(:,:,:) = IWORK(:,:,:)
    ENDIF
    !
#ifdef SFX_MPI     
    IF (LSFX_MPI) THEN
      XTIME_COMM_WRITE = XTIME_COMM_WRITE + (MPI_WTIME() - XTIME0) 
      XTIME0 = MPI_WTIME()
    ENDIF
#endif     
    !
    ICPT = 0
    !  
    DO J=1,SIZE(NINDEX)
      !
      IF ( NINDEX(J)==I ) THEN
        ICPT = ICPT + 1
        KWORK2(J,:,:) = IINTER(ICPT,:,:)
      ENDIF
      !
    ENDDO
    !
#ifdef SFX_MPI     
    IF (LSFX_MPI) XTIME_CALC_WRITE = XTIME_CALC_WRITE + (MPI_WTIME() - XTIME0)
#endif     
    !
  ENDDO
  !
ENDIF
!
IF (LHOOK) CALL DR_HOOK('GATHER_AND_WRITE_MPI_N3D',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE GATHER_AND_WRITE_MPI_N3D
!
!**************************************************************************
!
SUBROUTINE GATHER_AND_WRITE_MPI_X1D(PWORK,PWORK2,KMASK,PMISS)
!
USE MODI_GATHER_AND_WRITE_MPI_K4
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO, LSFX_MPI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
REAL, DIMENSION(:), INTENT(IN) :: PWORK
REAL(KIND=KIND(PWORK)), DIMENSION(:), INTENT(OUT) :: PWORK2
INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: KMASK
REAL, INTENT(IN), OPTIONAL :: PMISS
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('GATHER_AND_WRITE_MPI_X1D',0,ZHOOK_HANDLE)
!
IF (PRESENT(KMASK)) THEN
  IF (PRESENT(PMISS)) THEN
    CALL GATHER_AND_WRITE_MPI_K4(PWORK,PWORK2,KMASK,PMISS)
  ELSE
    CALL GATHER_AND_WRITE_MPI_K4(PWORK,PWORK2,KMASK)
  END IF
ELSE
  CALL GATHER_AND_WRITE_MPI_K4(PWORK,PWORK2)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('GATHER_AND_WRITE_MPI_X1D',1,ZHOOK_HANDLE)
!
END SUBROUTINE GATHER_AND_WRITE_MPI_X1D
!
!**************************************************************************
!
SUBROUTINE GATHER_AND_WRITE_MPI_X2D(PWORK,PWORK2,KMASK,PMISS)
!
USE MODI_GATHER_AND_WRITE_MPI_K4
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO, LSFX_MPI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:), INTENT(IN) :: PWORK
REAL(KIND=KIND(PWORK)), DIMENSION(:,:), INTENT(OUT) :: PWORK2
INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: KMASK
REAL, INTENT(IN), OPTIONAL :: PMISS 
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('GATHER_AND_WRITE_MPI_X2D',0,ZHOOK_HANDLE)
!
IF (PRESENT(KMASK)) THEN
  IF(PRESENT(PMISS)) THEN
    CALL GATHER_AND_WRITE_MPI_K4(PWORK,PWORK2,KMASK,PMISS)
  ELSE
    CALL GATHER_AND_WRITE_MPI_K4(PWORK,PWORK2,KMASK)
  END IF
ELSE
  CALL GATHER_AND_WRITE_MPI_K4(PWORK,PWORK2)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('GATHER_AND_WRITE_MPI_X2D',1,ZHOOK_HANDLE)
!
END SUBROUTINE GATHER_AND_WRITE_MPI_X2D
!
!**************************************************************************
!
SUBROUTINE GATHER_AND_WRITE_MPI_X3D(PWORK,PWORK2,KMASK,PMISS)
!
USE MODI_GATHER_AND_WRITE_MPI_K4
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO, LSFX_MPI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PWORK
REAL(KIND=KIND(PWORK)), DIMENSION(:,:,:), INTENT(OUT) :: PWORK2
INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: KMASK
REAL, INTENT(IN), OPTIONAL :: PMISS 
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('GATHER_AND_WRITE_MPI_X3D',0,ZHOOK_HANDLE)
!
IF (PRESENT(KMASK)) THEN
  IF(PRESENT(PMISS)) THEN
    CALL GATHER_AND_WRITE_MPI_K4(PWORK,PWORK2,KMASK,PMISS)
  ELSE
    CALL GATHER_AND_WRITE_MPI_K4(PWORK,PWORK2,KMASK)
  END IF
ELSE
  CALL GATHER_AND_WRITE_MPI_K4(PWORK,PWORK2)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('GATHER_AND_WRITE_MPI_X3D',1,ZHOOK_HANDLE)
!
END SUBROUTINE GATHER_AND_WRITE_MPI_X3D
!
!**************************************************************************
!
SUBROUTINE GATHER_AND_WRITE_MPI_X1DK4(PWORK,PWORK2,KMASK,PMISS)
!
USE MODI_GATHER_AND_WRITE_MPI_K4
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO, LSFX_MPI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
REAL, DIMENSION(:), INTENT(IN) :: PWORK
REAL(KIND=4), DIMENSION(:), INTENT(OUT) :: PWORK2
INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: KMASK
REAL, INTENT(IN), OPTIONAL :: PMISS 
!
REAL, DIMENSION(:), ALLOCATABLE :: ZINTER
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('GATHER_AND_WRITE_MPI_X1DK4',0,ZHOOK_HANDLE)
!
ALLOCATE(ZINTER(SIZE(PWORK2)))
IF (PRESENT(KMASK)) THEN
  IF(PRESENT(PMISS)) THEN
    CALL GATHER_AND_WRITE_MPI_K4(PWORK,ZINTER,KMASK,PMISS)
  ELSE
    CALL GATHER_AND_WRITE_MPI_K4(PWORK,ZINTER,KMASK)
  END IF
ELSE
  CALL GATHER_AND_WRITE_MPI_K4(PWORK,ZINTER)
ENDIF
!
IF (NRANK==NPIO) THEN
  PWORK2(:) = ZINTER(:)
ENDIF
DEALLOCATE(ZINTER)
!
IF (LHOOK) CALL DR_HOOK('GATHER_AND_WRITE_MPI_X1DK4',1,ZHOOK_HANDLE)
!
END SUBROUTINE GATHER_AND_WRITE_MPI_X1DK4
!
!**************************************************************************
!
SUBROUTINE GATHER_AND_WRITE_MPI_X2DK4(PWORK,PWORK2,KMASK,PMISS)
!
USE MODI_GATHER_AND_WRITE_MPI_K4
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO, LSFX_MPI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:), INTENT(IN) :: PWORK
REAL(KIND=4), DIMENSION(:,:), INTENT(OUT) :: PWORK2
INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: KMASK
REAL, INTENT(IN), OPTIONAL :: PMISS
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZINTER
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('GATHER_AND_WRITE_MPI_X2DK4',0,ZHOOK_HANDLE)
!
ALLOCATE(ZINTER(SIZE(PWORK2,1),SIZE(PWORK2,2)))
IF (PRESENT(KMASK)) THEN
  IF (PRESENT(PMISS)) THEN
    CALL GATHER_AND_WRITE_MPI_K4(PWORK,ZINTER,KMASK,PMISS)
  ELSE
    CALL GATHER_AND_WRITE_MPI_K4(PWORK,ZINTER,KMASK)
  END IF
ELSE
  CALL GATHER_AND_WRITE_MPI_K4(PWORK,ZINTER)
ENDIF
!
IF (NRANK==NPIO) THEN
  PWORK2(:,:) = ZINTER(:,:)
ENDIF
DEALLOCATE(ZINTER)
!
IF (LHOOK) CALL DR_HOOK('GATHER_AND_WRITE_MPI_X2DK4',1,ZHOOK_HANDLE)
!
END SUBROUTINE GATHER_AND_WRITE_MPI_X2DK4

!**************************************************************************
!
SUBROUTINE GATHER_AND_WRITE_MPI_X3DK4(PWORK,PWORK2,KMASK,PMISS)
!
USE MODI_GATHER_AND_WRITE_MPI_K4
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO, LSFX_MPI
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PWORK
REAL(KIND=4), DIMENSION(:,:,:), INTENT(OUT) :: PWORK2
INTEGER, DIMENSION(:), INTENT(IN), OPTIONAL :: KMASK
REAL, INTENT(IN), OPTIONAL :: PMISS
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZINTER
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('GATHER_AND_WRITE_MPI_X3DK4',0,ZHOOK_HANDLE)
!
ALLOCATE(ZINTER(SIZE(PWORK2,1),SIZE(PWORK2,2),SIZE(PWORK2,3)))
IF (PRESENT(KMASK)) THEN
  IF (PRESENT(PMISS)) THEN
    CALL GATHER_AND_WRITE_MPI_K4(PWORK,ZINTER,KMASK,PMISS)
  ELSE
    CALL GATHER_AND_WRITE_MPI_K4(PWORK,ZINTER,KMASK)
  END IF
ELSE
  CALL GATHER_AND_WRITE_MPI_K4(PWORK,ZINTER)
ENDIF
!
IF (NRANK==NPIO) THEN
  PWORK2(:,:,:) = ZINTER(:,:,:)
ENDIF
DEALLOCATE(ZINTER)
!
IF (LHOOK) CALL DR_HOOK('GATHER_AND_WRITE_MPI_X3DK4',1,ZHOOK_HANDLE)
!
END SUBROUTINE GATHER_AND_WRITE_MPI_X3DK4
