!     #######################################################
      SUBROUTINE END_IO_SURF_NC_n(HPROGRAM)
!     #######################################################
!
!!****  *END_IO_SURF_NC_n* - routine to close IO files
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
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
!!	S.Malardel   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    09/2003 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO
!
USE MODD_IO_SURF_NC, ONLY : LMASK, NID_NC
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INCLUDE "netcdf.inc"
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! main program
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: IRET
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('END_IO_SURF_NC_N',0,ZHOOK_HANDLE)
!
!$OMP BARRIER
!
IF (NRANK==NPIO) THEN
!$OMP SINGLE
  IRET = NF_CLOSE(NID_NC)
!$OMP END SINGLE
ENDIF
!
LMASK = .FALSE.
IF (LHOOK) CALL DR_HOOK('END_IO_SURF_NC_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE END_IO_SURF_NC_n
