!     #########
      SUBROUTINE GET_VAR_TOWN_n(HPROGRAM,KI,PQS,PZ0,PZ0H)
!     ###################################################
!
!!****  *GET_VAR_TOWN_n* - routine to get variables defined only over town
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
!!      P. Le Moigne *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    02/2006
!       M. Jidane   08/2008 Z0 and Z0H recovery from town tiles
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODI_GET_LUOUT
USE MODD_SURF_PAR,   ONLY   : XUNDEF
!
USE MODD_DIAG_TEB_n, ONLY : DGT => DIAG_TEB
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
 CHARACTER(LEN=6),     INTENT(IN)     :: HPROGRAM
INTEGER,              INTENT(IN)     :: KI      ! Number of points
REAL, DIMENSION(KI),  INTENT(OUT)    :: PQS     ! surface humidity
REAL, DIMENSION(KI),  INTENT(OUT)    :: PZ0     ! surface roughness length
REAL, DIMENSION(KI),  INTENT(OUT)    :: PZ0H    ! surface roughness length for heat
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: ILUOUT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GET_VAR_TOWN_N',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!-------------------------------------------------------------------------------
!
IF (DGT%LSURF_VARS) THEN 
        PQS      = DGT%XQS      
   ELSE 
        PQS      = XUNDEF      
ENDIF           
IF (DGT%LCOEF) THEN 
        PZ0      = DGT%XZ0      
        PZ0H     = DGT%XZ0H
   ELSE 
        PZ0      = XUNDEF      
        PZ0H     = XUNDEF
ENDIF           
IF (LHOOK) CALL DR_HOOK('GET_VAR_TOWN_N',1,ZHOOK_HANDLE)
!
!==============================================================================
!
END SUBROUTINE GET_VAR_TOWN_n
