!     #########
      SUBROUTINE WRITE_FLAKE_n (DTCO, DGU, IOB, U, FM, &
                                HPROGRAM,HWRITE)
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
!
USE MODD_SURFEX_n, ONLY : FLAKE_MODEL_t
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DIAG_SURF_ATM_n, ONLY : DIAG_SURF_ATM_t
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_WRITE_SURF_ATM, ONLY : LNOWRITE_CANOPY
USE MODI_INIT_IO_SURF_n
USE MODI_WRITESURF_FLAKE_n
USE MODI_WRITESURF_FLAKE_SBL_n
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
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DIAG_SURF_ATM_t), INTENT(INOUT) :: DGU
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(FLAKE_MODEL_t), INTENT(INOUT) :: FM
!
 CHARACTER(LEN=6),    INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=3),    INTENT(IN)  :: HWRITE    ! 'PREP' : does not write SBL XUNDEF fields
!                                             ! 'ALL' : all fields are written
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('WRITE_FLAKE_N',0,ZHOOK_HANDLE)
CALL INIT_IO_SURF_n(DTCO, DGU, IOB, U, &
                     HPROGRAM,'WATER ','FLAKE ','WRITE')
!
!*       1.     Selection of surface scheme
!               ---------------------------
!
 CALL WRITESURF_FLAKE_CONF_n(FM%CHF, FM%DGMF, FM%F, &
                             HPROGRAM)
 CALL WRITESURF_FLAKE_n(DGU, IOB, U, &
                        FM%F, &
                        HPROGRAM)
!
IF ((.NOT.LNOWRITE_CANOPY).OR.DGU%LSELECT) CALL WRITESURF_FLAKE_SBL_n(DGU, IOB, U, &
                                                                      FM%F, FM%FSB, &
                                                                      HPROGRAM,HWRITE)
!
!
!-------------------------------------------------------------------------------
!
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('WRITE_FLAKE_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_FLAKE_n
