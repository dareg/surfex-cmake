!     #########
      SUBROUTINE WRITE_TEB_n(HPROGRAM,HWRITE)
!     ####################################
!
!!****  *WRITE_TEB_n* - routine to write surface variables in their respective files
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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2003
!!      B. Decharme 07/2011 : Suppress pgd output
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_TEB_CANOPY_n, ONLY : TCP => TEB_CANOPY
!
USE MODD_WRITE_SURF_ATM, ONLY : LNOWRITE_CANOPY
USE MODD_DIAG_SURF_ATM_n, ONLY : DGU => DIAG_SURF_ATM
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS
!
USE MODI_INIT_IO_SURF_n
USE MODI_WRITESURF_TEB_n
USE MODI_WRITESURF_TEB_CONF_n
USE MODI_END_IO_SURF_n
USE MODI_WRITESURF_TEB_CANOPY_n
USE MODI_GOTO_TEB
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
 CHARACTER(LEN=3),    INTENT(IN)  :: HWRITE    ! 'PREP' : does not write SBL XUNDEF fields
!                                             ! 'ALL' : all fields are written
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JPATCH
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!
!         Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('WRITE_TEB_N',0,ZHOOK_HANDLE)
 CALL INIT_IO_SURF_n(HPROGRAM,'TOWN  ','TEB   ','WRITE')
!
!*       1.     Selection of surface scheme
!               ---------------------------
!
 CALL WRITESURF_TEB_CONF_n(HPROGRAM)
!
DO JPATCH=1,TOP%NTEB_PATCH
  CALL GOTO_TEB(JPATCH)
  CALL WRITESURF_TEB_n(HPROGRAM,JPATCH,HWRITE)
END DO
!     
 CALL GOTO_TEB(1)
IF ((.NOT.LNOWRITE_CANOPY).OR.DGU%LSELECT) CALL WRITESURF_TEB_CANOPY_n(TCP, TOP, &
                                                                       HPROGRAM,HWRITE)
!
!-------------------------------------------------------------------------------
!
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('WRITE_TEB_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_TEB_n
