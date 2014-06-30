! *****************************************************************************************
PROGRAM SODA
! ******************************************************************************************
!
USE MODD_SURFEX_MPI, ONLY : NCOMM, NPROC, NRANK
!
USE MODI_SODA_CONTROL
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
INTEGER :: ILEVEL, INFOMPI
REAL (KIND=JPRB) :: ZHOOK_HANDLE
!
#ifndef NOMPI
 CALL MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,ILEVEL,INFOMPI)
#endif
!
IF (LHOOK) CALL DR_HOOK('SODA',0,ZHOOK_HANDLE)
!
!
#ifndef NOMPI
NCOMM = MPI_COMM_WORLD
 CALL MPI_COMM_SIZE(NCOMM,NPROC,INFOMPI)
 CALL MPI_COMM_RANK(NCOMM,NRANK,INFOMPI)
#endif
!
 CALL SODA_CONTROL (ODINLINE = .FALSE.)
!
IF (LHOOK) CALL DR_HOOK('SODA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END PROGRAM SODA
