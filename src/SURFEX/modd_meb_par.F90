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
!!	P. Samuelsson
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
! Extinction coefficient for view factor for short-wave radiation 
!
REAL, SAVE       :: XTAU_V_CF
!
! MEB resistance increase factor for canopy air sapce
!
REAL, SAVE       :: XRAGNC_FACTOR
!
! MEB maximum vegetation-intercepted water fraction
!
REAL, SAVE       :: XKDELTA_WR
!
!-------------------------------------------------------------------------------
!
END MODULE MODD_MEB_PAR












