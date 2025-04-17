!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!###################
FUNCTION CHTC_ROUGH_DOE_0D(PCHTCN, PCHTCS, PROUGH) RESULT(PCHTC)
!#############################################
!
!!****  *CHTC_ROUGH_DOE* - 
!!
!!    PURPOSE
!!    -------
!     function to compute convective surface coefficient for
!     a rough surface from the natural convection coef, the
!     smooth surface convective coef and the roughness coef
!
!!**  METHOD
!!    ------
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      
!!    REFERENCE
!!    ---------
!!    EnergyPlus, Engineering Reference, DOE-2 model for convection on outside
!!    surfaces, p68
!!
!!    AUTHOR
!!    ------
!!      G. Pigeon       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/00/2012 
!
!-------------------------------------------------------------------------------
!*       0.    DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!
REAL, INTENT(IN)                :: PCHTCN  ! Convective heat transfer coefficient for natural conv. [W/(m2.K)]
REAL, INTENT(IN)                :: PCHTCS  ! Convective heat transfer coefficient for a smooth surface [W/(m2.K)]
REAL, INTENT(IN)                :: PROUGH  ! Convective heat transfer coefficient for a smooth surface [W/(m2.K)]
REAL                            :: PCHTC   ! Convective heat transfer coefficient [W/(m2.K)]
!
!*       0.2   Declarations of local variables
!
!----------------------------------------------------------------------
!
!*       1.    COMPUTE THE CHTC
!              ----------------
!
PCHTC = PCHTCN + PROUGH * (PCHTCS - PCHTCN)
!
END FUNCTION CHTC_ROUGH_DOE_0D
