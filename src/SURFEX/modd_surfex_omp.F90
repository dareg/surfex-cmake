!     ######################
      MODULE MODD_SURFEX_OMP
!     ######################
!
!!****  *MODD_SURFEX_OMP
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      S. Faroux   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       26/06/12
!!      Modified    11/2013 by J.Escobar :add !$ to inhibit completly omp
!!                                 dependency
!
!*       0.   DECLARATIONS
!             ------------
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
#ifdef AIX64 
 USE OMP_LIB
#endif
!
IMPLICIT NONE
!
#ifndef AIX64
  INCLUDE 'omp_lib.h'
#endif
!
!RJ: this broke non openmp version before
!RJ: OMP_GET_THREAD_NUM() returns 0 for first omp thread
!RJ: OMP_GET_NUM_THREADS() returns 1 for omp thread count
#ifdef RJ_OFIX
INTEGER :: NBLOCKTOT = 1
INTEGER :: NBLOCK = 0
#else
INTEGER :: NBLOCKTOT = 1
INTEGER :: NBLOCK = 1
#endif
!$OMP THREADPRIVATE(NBLOCK)
INTEGER :: IDC = 0
!
END MODULE MODD_SURFEX_OMP

