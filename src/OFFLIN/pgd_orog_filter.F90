!     #########
      SUBROUTINE PGD_OROG_FILTER(HPROGRAM)
!     ##############################################################
!
!!**** *PGD_OROGRAPHY* monitor for averaging and interpolations of cover fractions
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    10/12/97
!!    12/2008 E. Martin : add case 'MAX' for choice of orography
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_PGD_GRID,       ONLY : CGRID, XGRID_PAR
USE MODD_SURF_ATM_n,     ONLY : XZS, XSEA
!
USE MODI_READ_NAM_PGD_OROG_FILTER
USE MODI_OROGRAPHY_FILTER
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
 CHARACTER(LEN=6),     INTENT(IN)  :: HPROGRAM ! program calling
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
!*    0.3    Declaration of namelists
!            ------------------------
!
INTEGER                  :: NZSFILTER   ! number of orographic spatial filter iterations
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*    1.      Initializations
!             ---------------
!
IF (LHOOK) CALL DR_HOOK('PGD_OROGRAPHY',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!*    2.      Reading of namelist
!             -------------------
!
 CALL READ_NAM_PGD_OROG_FILTER(HPROGRAM, NZSFILTER )  
!
!-------------------------------------------------------------------------------
!
!*   11.      Filtering of orography
!             ----------------------
!
 CALL OROGRAPHY_FILTER(CGRID, XGRID_PAR, XSEA, NZSFILTER, XZS)
!
IF (LHOOK) CALL DR_HOOK('PGD_OROG_FILTER',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE PGD_OROG_FILTER
