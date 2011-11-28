!     #########
      SUBROUTINE WRITE_FLAKE_n(HPROGRAM,HWRITE)
!     ####################################
!
!!****  *WRITE_FLAKE_n* - routine to write surface variables in their respective files
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
USE MODD_WRITE_SURF_ATM, ONLY : LNOWRITE_CANOPY, LNOWRITE_COVERS
USE MODD_DIAG_SURF_ATM_n, ONLY:LSELECT
USE MODI_INIT_IO_SURF_n
USE MODI_WRITESURF_FLAKE_n
USE MODI_WRITESURF_FLAKE_SBL_n
USE MODI_WRITESURF_PGD_FLAKE_n
USE MODI_WRITESURF_FLAKE_CONF_n
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
CHARACTER(LEN=3),    INTENT(IN)  :: HWRITE    ! 'PGD' : only physiographic fields are written
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!                                             ! 'ALL' : all fields are written
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!-------------------------------------------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('WRITE_FLAKE_N',0,ZHOOK_HANDLE)
CALL INIT_IO_SURF_n(HPROGRAM,'WATER ','FLAKE ','WRITE')
!
!*       1.     Selection of surface scheme
!               ---------------------------
!
IF (HWRITE/='PGD') THEN
   !
   CALL WRITESURF_FLAKE_CONF_n(HPROGRAM)
   CALL WRITESURF_FLAKE_n(HPROGRAM)
   !
   IF ((.NOT.LNOWRITE_CANOPY).OR.LSELECT) CALL WRITESURF_FLAKE_SBL_n(HPROGRAM)
   !
ENDIF
!
IF ((.NOT.LNOWRITE_COVERS).OR.LSELECT) CALL WRITESURF_PGD_FLAKE_n(HPROGRAM,HWRITE)
!
!-------------------------------------------------------------------------------
!
!         End of IO
!
CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('WRITE_FLAKE_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_FLAKE_n
