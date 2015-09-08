!     ##################
      MODULE MODD_SURF_ATM_GRID_n
!     ##################
!
!!****  *MODD_SURF_ATM_GRID - declaration of SURF_ATM grid
!!
!!    PURPOSE
!!    -------
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      V. Masson  *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2004
!
!*       0.   DECLARATIONS
!             ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

TYPE SURF_ATM_GRID_t
!-------------------------------------------------------------------------------
!
! Grid definition
!
  CHARACTER(LEN=10)               :: CGRID       ! grid type
!                                              ! "NONE        " : no grid computations
!                                              ! "CONF PROJ   " : conformal projection
!
  REAL, POINTER,     DIMENSION(:) :: XGRID_PAR   ! lits of parameters used to define the grid
!                                              ! (depends on value of CGRID)
  REAL, POINTER,     DIMENSION(:) :: XGRID_FULL_PAR   ! lits of parameters used to define the grid
!                                                     ! (depends on value of CGRID)
  INTEGER                         :: NGRID_PAR   ! size of XGRID_PAR
!
  INTEGER, POINTER, DIMENSION(:,:) :: NNEAR
!-------------------------------------------------------------------------------
!
! General surface parameters:
!
  REAL, POINTER, DIMENSION(:) :: XLAT        ! latitude (degrees +North)               (-)
  REAL, POINTER, DIMENSION(:) :: XLON        ! longitude (degrees +East)               (-)
  REAL, POINTER, DIMENSION(:) :: XMESH_SIZE  ! mesh size                               (m2)
  REAL, POINTER, DIMENSION(:) :: XJPDIR      ! heading of J direction (deg from N clockwise)
!-------------------------------------------------------------------------------
!

END TYPE SURF_ATM_GRID_t

TYPE(SURF_ATM_GRID_t), ALLOCATABLE, TARGET, SAVE :: SURF_ATM_GRID_MODEL(:)

TYPE(SURF_ATM_GRID_t), POINTER :: SURF_ATM_GRID => NULL()
!$OMP THREADPRIVATE(SURF_ATM_GRID)

CONTAINS

SUBROUTINE SURF_ATM_GRID_GOTO_MODEL(KFROM, KTO, LKFROM)
LOGICAL, INTENT(IN) :: LKFROM
INTEGER, INTENT(IN) :: KFROM, KTO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! Current model is set to model KTO
IF (LHOOK) CALL DR_HOOK('MODD_SURF_ATM_GRID_N:SURF_ATM_GRID_GOTO_MODEL',0,ZHOOK_HANDLE)

SURF_ATM_GRID => SURF_ATM_GRID_MODEL(KTO)

IF (LHOOK) CALL DR_HOOK('MODD_SURF_ATM_GRID_N:SURF_ATM_GRID_GOTO_MODEL',1,ZHOOK_HANDLE)

END SUBROUTINE SURF_ATM_GRID_GOTO_MODEL

SUBROUTINE SURF_ATM_GRID_ALLOC(KMODEL)
INTEGER, INTENT(IN) :: KMODEL
INTEGER :: J
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_SURF_ATM_GRID_N:SURF_ATM_GRID_ALLOC",0,ZHOOK_HANDLE)
ALLOCATE(SURF_ATM_GRID_MODEL(KMODEL))
SURF_ATM_GRID => SURF_ATM_GRID_MODEL(KMODEL)
DO J=1,KMODEL
  NULLIFY(SURF_ATM_GRID_MODEL(J)%XGRID_PAR)
  NULLIFY(SURF_ATM_GRID_MODEL(J)%NNEAR)
  NULLIFY(SURF_ATM_GRID_MODEL(J)%XGRID_FULL_PAR)
  NULLIFY(SURF_ATM_GRID_MODEL(J)%XLAT)
  NULLIFY(SURF_ATM_GRID_MODEL(J)%XLON)
  NULLIFY(SURF_ATM_GRID_MODEL(J)%XMESH_SIZE)
  NULLIFY(SURF_ATM_GRID_MODEL(J)%XJPDIR)
ENDDO
SURF_ATM_GRID_MODEL(:)%CGRID=' '
SURF_ATM_GRID_MODEL(:)%NGRID_PAR=0
IF (LHOOK) CALL DR_HOOK("MODD_SURF_ATM_GRID_N:SURF_ATM_GRID_ALLOC",1,ZHOOK_HANDLE)
END SUBROUTINE SURF_ATM_GRID_ALLOC

SUBROUTINE SURF_ATM_GRID_DEALLO
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK("MODD_SURF_ATM_GRID_N:SURF_ATM_GRID_DEALLO",0,ZHOOK_HANDLE)
IF (ALLOCATED(SURF_ATM_GRID_MODEL)) DEALLOCATE(SURF_ATM_GRID_MODEL)

IF (LHOOK) CALL DR_HOOK("MODD_SURF_ATM_GRID_N:SURF_ATM_GRID_DEALLO",1,ZHOOK_HANDLE)
END SUBROUTINE SURF_ATM_GRID_DEALLO

END MODULE MODD_SURF_ATM_GRID_n
