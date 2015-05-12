!     #########
SUBROUTINE WRITE_DIAG_SEAFLUX_n(HPROGRAM,HWRITE)
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

!
USE MODD_CH_SEAFLUX_n, ONLY : CHS => CH_SEAFLUX
USE MODD_DIAG_SURF_ATM_n, ONLY : DGU => DIAG_SURF_ATM
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_SEAFLUX_n, ONLY : S => SEAFLUX
USE MODD_DIAG_OCEAN_n, ONLY : DGO => DIAG_OCEAN
USE MODD_DIAG_SEAFLUX_n, ONLY : DGS => DIAG_SEAFLUX
USE MODD_DIAG_SEAICE_n, ONLY : DGSI => DIAG_SEAICE
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
   IF (DGS%XDIAG_TSTEP==XUNDEF .OR. ABS(NINT(S%TTIME%TIME/DGS%XDIAG_TSTEP)*DGS%XDIAG_TSTEP-S%TTIME%TIME)<1.E-3 ) THEN
      CALL WRITE_DIAG_SEB_SEAFLUX_n(CHS, DGS, DGSI, DGU, S, &
                                    HPROGRAM)
      IF (DGO%LDIAG_OCEAN)  CALL WRITE_DIAG_SEB_OCEAN_n(DGO, &
                                                        HPROGRAM)
      IF (DGSI%LDIAG_SEAICE) CALL WRITE_DIAG_SEB_SEAICE_n(DGS, DGSI, DGU, S, &
                                                          HPROGRAM)
   END IF
!        
ENDIF
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_SEAFLUX_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_DIAG_SEAFLUX_n
