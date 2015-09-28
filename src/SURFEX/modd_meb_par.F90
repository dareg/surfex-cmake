!     ######################
      MODULE MODD_MEB_PAR
!     ######################
!
!!****  *MODD_MEB_PAR* - declaration of parameters related
!!                          to the MEB parameterizations
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     parameters related to the surface parameterizations of MEB.
!
!!
!!      
!!
!!    AUTHOR
!!    ------
!!      P. Samuelsson
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       09/2013
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
! Extinction coefficient for view factor for long-wave radiation 
!
REAL, SAVE       :: XTAU_LW
!
! MEB resistance increase factor for canopy air sapce
!
REAL, SAVE       :: XRAGNC_FACTOR
!
! MEB maximum vegetation-intercepted water fraction
!
REAL, SAVE       :: XKDELTA_WR
!
! Constants
! ---------
!
! Cuurently MEB considers 2 spectral bands, with weights below.
! As more bands are eventually added, one must adjust these weights
! and add more so that the total = NSPB_MEB 
!
INTEGER, PARAMETER   :: NSPB_MEB     = 2      ! number of bands (below)
!
REAL,    PARAMETER   :: XSW_WGHT_VIS = 0.48
!
REAL,    PARAMETER   :: XSW_WGHT_NIR = 0.52
!
!-------------------------------------------------------------------------------
!
END MODULE MODD_MEB_PAR












