subroutine dr_hook_procinfo(kmyproc, knproc)
USE PARKIND1  ,ONLY : JPIM     ,JPRB
#ifdef SFX_MPL
use mpl_data_module, only : MPL_RANK,MPL_NUMPROC
#endif
implicit none
#ifdef SFX_MPI
INCLUDE 'mpif.h'
#endif
INTEGER(KIND=JPIM),intent(out) :: kmyproc, knproc
INTEGER(KIND=JPIM) :: IRNK,ISZ,IERR
LOGICAL :: LLINIT
#if defined(SFX_MPL)
kmyproc = mpl_rank
knproc = mpl_numproc
#else
IRNK=0
ISZ=1
#ifdef SFX_MPI
CALL MPI_INITIALIZED(LLINIT, IERR)
IF(LLINIT) THEN
  CALL MPI_COMM_RANK(MPI_COMM_WORLD,IRNK,IERR)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,ISZ,IERR)
ENDIF
#endif
kmyproc = IRNK+1
knproc = ISZ
#endif
end subroutine dr_hook_procinfo
