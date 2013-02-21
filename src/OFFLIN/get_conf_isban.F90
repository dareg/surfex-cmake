!     #########
SUBROUTINE GET_CONF_ISBA_n(OTRIP,OFLOOD,HGRID,KDIMTAB)
!###############################################
!
!!**** *GET_CONF_ISBA_n* get the ISBA options configuration
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    B. Decharme         Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    05/2008
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_DIAG_EVAP_ISBA_n, ONLY : LSURF_BUDGETC
USE MODD_ISBA_n,           ONLY : LTRIP, LFLOOD
USE MODD_ISBA_GRID_n,      ONLY : CGRID
USE MODD_SGH_PAR,          ONLY : NDIMTAB
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
LOGICAL,           INTENT(OUT)           :: OTRIP
LOGICAL,           INTENT(OUT)           :: OFLOOD
 CHARACTER(LEN=10), INTENT(OUT), OPTIONAL :: HGRID
INTEGER,           INTENT(OUT), OPTIONAL :: KDIMTAB
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*    0.2    Declaration of other local variables
!            ------------------------------------
!
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GET_CONF_ISBA_N',0,ZHOOK_HANDLE)
OTRIP = LTRIP
!
OFLOOD = LFLOOD
!
IF(PRESENT(HGRID)) HGRID = CGRID
!
IF(PRESENT(KDIMTAB)) KDIMTAB = NDIMTAB
IF (LHOOK) CALL DR_HOOK('GET_CONF_ISBA_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE GET_CONF_ISBA_n
