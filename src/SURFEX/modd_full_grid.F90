!     ##################
      MODULE MODD_FULL_GRID
!     ##################
!
!!****  *MODD_FULL_GRID - declaration of SURF_ATM grid
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
 CHARACTER(LEN=10)               :: CGRID       ! grid type
!                                              ! "NONE        " : no grid computations
!                                              ! "CONF PROJ   " : conformal projection
 REAL, POINTER,     DIMENSION(:) :: XGRID_FULL_PAR   ! lits of parameters used to define the grid
!                                                     ! (depends on value of CGRID)
 INTEGER                         :: NGRID_PAR   ! size of XGRID_PAR
!
!
END MODULE MODD_FULL_GRID
