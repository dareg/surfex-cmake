!     #########
      SUBROUTINE WRITE_PGD_SEAFLUX_n (DTS, SG, S, &
                                      HPROGRAM)
!     ####################################
!
!!****  *WRITE_PGD_SEAFLUX_n* - routine to write pgd surface variables in their respective files
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
!
!
USE MODD_DATA_SEAFLUX_n, ONLY : DATA_SEAFLUX_t
USE MODD_SEAFLUX_GRID_n, ONLY : SEAFLUX_GRID_t
USE MODD_SEAFLUX_n, ONLY : SEAFLUX_t
!
USE MODI_INIT_IO_SURF_n
USE MODI_WRITESURF_PGD_SEAFLUX_n
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
!
TYPE(DATA_SEAFLUX_t), INTENT(INOUT) :: DTS
TYPE(SEAFLUX_GRID_t), INTENT(INOUT) :: SG
TYPE(SEAFLUX_t), INTENT(INOUT) :: S
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
!         Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('WRITE_PGD_SEAFLUX_N',0,ZHOOK_HANDLE)
 CALL INIT_IO_SURF_n(HPROGRAM,'SEA   ','SEAFLX','WRITE')
!
!*       1.     Selection of surface scheme
!               ---------------------------
!
 CALL WRITESURF_PGD_SEAFLUX_n(DTS, SG, S, &
                              HPROGRAM)
!
!-------------------------------------------------------------------------------
!
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('WRITE_PGD_SEAFLUX_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_PGD_SEAFLUX_n
