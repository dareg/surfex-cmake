!     ##################
      MODULE MODN_ASSIM_GARDEN
!     ##################
!
!!****  *MODN_ASSIM_GARDEN - declaration of keys for ISBA assimilation scheme (2DVAR, Bouyssel et al.)
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
USE MODD_ASSIM_GARDEN, ONLY:  LASSIM, CASSIM     
!
IMPLICIT NONE
!
NAMELIST/NAM_ASSIM/LASSIM,CASSIM
!
END MODULE MODN_ASSIM_GARDEN
