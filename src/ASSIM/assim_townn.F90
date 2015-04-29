!     ###############################################################################
SUBROUTINE ASSIM_TOWN_n(HPROGRAM,KI,PT2M_O,HTEST)

!     ###############################################################################
!
!!****  *ASSIM_TOWN_n * - Chooses the surface schemes for TOWN parts  
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     T. Aspelien
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/2012
!!--------------------------------------------------------------------
!
USE MODD_TEB_n, ONLY : T => TEB
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS
!
USE MODD_CSTS,          ONLY : XPI
USE MODN_IO_OFFLINE,    ONLY : CSURF_FILETYPE
!
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
USE YOMHOOK,            ONLY : LHOOK,   DR_HOOK
USE PARKIND1,           ONLY : JPRB
!
USE MODI_ABOR1_SFX
USE MODI_ASSIM_TEB_n
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=6),   INTENT(IN) :: HPROGRAM  ! program calling surf. schemes
INTEGER,            INTENT(IN) :: KI
REAL,DIMENSION(kI), INTENT(IN) :: PT2M_O
CHARACTER(LEN=2),   INTENT(IN) :: HTEST ! must be equal to 'OK'
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
!
REAL(KIND=JPRB)                    :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('ASSIM_TOWN_N',0,ZHOOK_HANDLE)

IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('ASSIM_TOWN_n: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF

IF (U%CTOWN=='TEB   ') THEN
  CALL ASSIM_TEB_n(U, T, TOP, &
                   HPROGRAM,KI,PT2M_O,HTEST)
ELSE
  WRITE(*,*) 'No assimilation done for scheme: ',TRIM(U%CTOWN)
END IF

IF (LHOOK) CALL DR_HOOK('ASSIM_TOWN_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE ASSIM_TOWN_n
