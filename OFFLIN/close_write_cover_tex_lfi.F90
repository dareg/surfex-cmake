!     ##################################
      SUBROUTINE CLOSE_WRITE_COVER_TEX_LFI(KTEX)
!     ##################################
!
!!****  *CLOSE_WRITE_COVER_TEX_LFI* - closes cover listing file
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
USE MODD_IO_SURF_LFI,  ONLY : CLUOUT_LFI
!
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER, INTENT(IN) :: KTEX  ! logical unit to close
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
CHARACTER(LEN=20) :: YTEX           ! name of tex file
INTEGER :: IRESP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!*       5.     Prints of cover parameters in a tex file
!               ----------------------------------------
!
IF (LHOOK) CALL DR_HOOK('CLOSE_WRITE_COVER_TEX_LFI',0,ZHOOK_HANDLE)
YTEX = 'class_cover_data.tex'
CALL FMFREE(YTEX,CLUOUT_LFI,IRESP)
CLOSE(KTEX)
IF (LHOOK) CALL DR_HOOK('CLOSE_WRITE_COVER_TEX_LFI',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE CLOSE_WRITE_COVER_TEX_LFI
