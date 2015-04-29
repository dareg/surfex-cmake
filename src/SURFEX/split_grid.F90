!     ###########################################################
      SUBROUTINE SPLIT_GRID(HPROGRAM)
!     ###########################################################
!!
!!    PURPOSE
!!    -------
!!   This program splits a PGD grid on several processors (according to host program)
!!
!!    METHOD
!!    ------
!!   
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson                   Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original     08/11
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_SURF_ATM_GRID_n, ONLY : UG => SURF_ATM_GRID
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
USE MODI_SPLIT_GRID_CONF_PROJ
USE MODI_SPLIT_GRID_CARTESIAN
USE MODI_GET_SIZE_FULL_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
!
IMPLICIT NONE
!
!*    0.1    Declaration of dummy arguments
!            ------------------------------
!
 CHARACTER(LEN=6),     INTENT(IN)  :: HPROGRAM ! program calling
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
 CHARACTER(LEN=100) :: YCOMMENT
INTEGER :: IRESP ! error return code
!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SPLIT_GRID',0,ZHOOK_HANDLE)
!
SELECT CASE(UG%CGRID)

  CASE('CONF PROJ ')
    CALL SPLIT_GRID_CONF_PROJ(HPROGRAM,U%NDIM_FULL,U%NSIZE_FULL,UG%NGRID_PAR,UG%XGRID_PAR)
  CASE('CARTESIAN ')
    CALL SPLIT_GRID_CARTESIAN(HPROGRAM,U%NDIM_FULL,U%NSIZE_FULL,UG%NGRID_PAR,UG%XGRID_PAR)
  CASE DEFAULT
    CALL GET_SIZE_FULL_n(U, &
                         HPROGRAM,U%NDIM_FULL,U%NSIZE_FULL)

END SELECT
!

IF (LHOOK) CALL DR_HOOK('SPLIT_GRID',1,ZHOOK_HANDLE)
!_______________________________________________________________________________
!
END SUBROUTINE SPLIT_GRID
