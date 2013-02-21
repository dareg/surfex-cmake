!     #########
      SUBROUTINE CLOSE_NAMELIST_LFI(HPROGRAM,KLUNAM)
!     #######################################################
!
!!****  *CLOSE_NAMELIST_LFI* - closes namelists read by surface in MESOHN
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
!!	V. Masson   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2003 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_FMDECLAR,    ONLY : CNAMFI
USE MODD_IO_SURF_LFI, ONLY : CLUOUT_LFI

!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
 CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! main program
INTEGER,           INTENT(IN)  :: KLUNAM   ! logical unit of namelist
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: IRESP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!* closes the namelist
!  -------------------
!
IF (LHOOK) CALL DR_HOOK('CLOSE_NAMELIST_LFI',0,ZHOOK_HANDLE)
 CALL FMFREE(CNAMFI(KLUNAM),CLUOUT_LFI,IRESP)
CLOSE(KLUNAM)
IF (LHOOK) CALL DR_HOOK('CLOSE_NAMELIST_LFI',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE CLOSE_NAMELIST_LFI
