
!     #########
      SUBROUTINE MAKE_LCOVER(OCOVER)
!     ##############################################################
!
!!**** *PGD_COVER* monitor for averaging and interpolations of cover fractions
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    10/12/97
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO, NPROC, NCOMM
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
#ifndef NOMPI
INCLUDE "mpif.h"
#endif
!
!*    0.1    Declaration of arguments
!            ------------------------
!
LOGICAL, DIMENSION(:), INTENT(INOUT) :: OCOVER
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
INTEGER :: INFOMPI, JPROC, JCOVER
!
LOGICAL, DIMENSION(:,:), ALLOCATABLE :: GCOVER_ALL
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!---------------------------------------------------------------
!
!*    1.      Initializations
!             ---------------
!
IF (LHOOK) CALL DR_HOOK('MAKE_LCOVER',0,ZHOOK_HANDLE)
!
IF (NRANK==NPIO) THEN
  ALLOCATE(GCOVER_ALL(SIZE(OCOVER),0:NPROC-1))
ELSE
  ALLOCATE(GCOVER_ALL(0,0))
ENDIF

IF (NPROC>1) THEN
#ifndef NOMPI
  CALL MPI_GATHER(OCOVER,SIZE(OCOVER),MPI_LOGICAL,GCOVER_ALL,SIZE(OCOVER),&
                  MPI_LOGICAL,NPIO,NCOMM,INFOMPI)
#endif
ELSE
  GCOVER_ALL(:,0) = OCOVER(:)
ENDIF
!
IF (NRANK==NPIO) THEN
  OCOVER(:) = .FALSE.
  DO JPROC = 0,NPROC-1
    DO JCOVER=1,SIZE(OCOVER)
      IF (GCOVER_ALL(JCOVER,JPROC)) OCOVER(JCOVER) = .TRUE.
    ENDDO
  ENDDO
ENDIF
DEALLOCATE(GCOVER_ALL)
!
IF (NPROC>1) THEN
#ifndef NOMPI
  CALL MPI_BCAST(OCOVER,SIZE(OCOVER),MPI_LOGICAL,NPIO,NCOMM,INFOMPI)
#endif
ENDIF
!
IF (LHOOK) CALL DR_HOOK('MAKE_LCOVER',1,ZHOOK_HANDLE)
!
END SUBROUTINE MAKE_LCOVER
