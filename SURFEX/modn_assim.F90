!     ##################
      MODULE MODN_ASSIM
!     ##################
!
!!****  *MODN_ASSIM - declaration of keys for ISBA assimilation scheme (2DVAR, Bouyssel et al.)
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
!!	L. Jarlan   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       23/02/05
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_ASSIM, ONLY:  LASSIM, CASSIM     
!
IMPLICIT NONE
!
NAMELIST/NAM_ASSIM/LASSIM,CASSIM
!
END MODULE MODN_ASSIM
