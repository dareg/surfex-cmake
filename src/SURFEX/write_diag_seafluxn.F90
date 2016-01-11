!     #########
SUBROUTINE WRITE_DIAG_SEAFLUX_n (DTCO, DGU, U, SM, &   
                                 HPROGRAM,HWRITE)
!     ###############################################################################
!
!!****  *WRITE_DIAG_SEAFLUX_n * - diagnostics for SEAFLUX
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     V. Masson 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      Modified    09/2013 : S. Senesi : call WRITE_DIAG_SEB_SEAICE_n
!!------------------------------------------------------------------
!
USE MODD_SFX_OASIS,      ONLY : LCPL_SEAICE
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DIAG_n, ONLY : DIAG_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
!
USE MODD_SURFEX_n, ONLY : SEAFLUX_MODEL_t
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
!
USE MODI_WRITE_DIAG_SEB_SEAFLUX_n
USE MODI_WRITE_DIAG_SEB_OCEAN_n
USE MODI_WRITE_DIAG_SEB_SEAICE_n
! 
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DIAG_t), INTENT(INOUT) :: DGU
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SEAFLUX_MODEL_t), INTENT(INOUT) :: SM
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM ! program calling surf. schemes
 CHARACTER(LEN=3),   INTENT(IN)  :: HWRITE   ! 'PGD' : only physiographic fields are written
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!                                           ! 'ALL' : all fields are written
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_SEAFLUX_N',0,ZHOOK_HANDLE)
IF (HWRITE/='PGD') THEN
!        
   IF (SM%DGS%XDIAG_TSTEP==XUNDEF .OR. &
           ABS(NINT(SM%S%TTIME%TIME/SM%DGS%XDIAG_TSTEP)*SM%DGS%XDIAG_TSTEP-SM%S%TTIME%TIME)<1.E-3 ) THEN
      CALL WRITE_DIAG_SEB_SEAFLUX_n(DTCO, DGU, U, SM%CHS, SM%DGS, SM%DGSC, SM%DGMSI, SM%S, &
                                    HPROGRAM)
      IF (SM%DGO%LDIAG_OCEAN)  CALL WRITE_DIAG_SEB_OCEAN_n(DTCO, DGU, U, SM%DGO, &
                                                        HPROGRAM)
      IF (SM%S%LHANDLE_SIC.OR.LCPL_SEAICE) CALL WRITE_DIAG_SEB_SEAICE_n(DTCO, DGU, U, SM%DGS, SM%DGSI, SM%DGSIC, SM%S, &
                                                          HPROGRAM)                                                
      IF (SM%DGMSI%LDIAG_MISC_SEAICE) CALL WRITE_DIAG_MISC_SEAICE_n(DTCO, DGU, U, SM%DGMSI, SM%S, &
                                                          HPROGRAM)
   END IF
!        
ENDIF
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_SEAFLUX_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_DIAG_SEAFLUX_n
