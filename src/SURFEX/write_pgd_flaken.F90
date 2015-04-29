!     #########
      SUBROUTINE WRITE_PGD_FLAKE_n(HPROGRAM)
!     ####################################
!
!!****  *WRITE_PGD_FLAKE_n* - routine to write pgd surface variables in their respective files
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
!!      B. Decharme   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/2011
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_FLAKE_GRID_n, ONLY : FG => FLAKE_GRID
USE MODD_FLAKE_n, ONLY : F => FLAKE
!
USE MODI_INIT_IO_SURF_n
USE MODI_WRITESURF_PGD_FLAKE_n
USE MODI_END_IO_SURF_n
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
 CHARACTER(LEN=6),    INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!                                             ! 'ALL' : all fields are written
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('WRITE_PGD_FLAKE_N',0,ZHOOK_HANDLE)
 CALL INIT_IO_SURF_n(HPROGRAM,'WATER ','FLAKE ','WRITE')
!
!*       1.     Selection of surface scheme
!               ---------------------------
!
 CALL WRITESURF_PGD_FLAKE_n(FG, F, &
                            HPROGRAM)
!
!-------------------------------------------------------------------------------
!
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('WRITE_PGD_FLAKE_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_PGD_FLAKE_n
