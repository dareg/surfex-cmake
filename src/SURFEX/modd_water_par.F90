!     ######################
      MODULE MODD_WATER_PAR
!     ######################
!
!!****  *MODD_WATER_PAR* - declaration of parameters related
!!                          to the water parameterization
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     parameters related to the surface parameterization of sea or
!     water.
!
!!
!!      
!!
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       01/2004                     
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
REAL, SAVE       :: XALBWAT 
!                   water global albedo (option "UNIF")
!
REAL, SAVE       :: XALBSCA_WAT
!                   water diffuse albedo
!
REAL, SAVE       :: XALBCOEF_TA96
!                   water direct albedo coefficient (TA96 computation)
!
REAL, SAVE       :: XEMISWAT
!                   water emissivity
!
REAL, SAVE       :: XALBSEAICE 
!                   sea ice global albedo
!
REAL, SAVE       :: XALBWATICE 
!                   water ice global albedo
!
REAL, SAVE       :: XALBWATSNOW 
!                   water snow global albedo (for lake)
!
REAL, SAVE       :: XEMISWATICE 
!                   sea ice emissivity
!-------------------------------------------------------------------------------
!
END MODULE MODD_WATER_PAR












