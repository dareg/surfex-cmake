!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!###################
FUNCTION CHTC_SMOOTH_LEE_DOE_0D(PCHTCN, PVMOD) RESULT(PCHTC)
!
!!****  *CHTC_SMOOTH_LEE_DOE* - 
!!
!!    PURPOSE
!!    -------
!     function to compute convective surface coefficient for
!     a leeward smooth surface from the natural convection coef and the
!     wind speed
!
!!**  METHOD
!!    ------
!!
!!    from EnergyPlus Engineering Reference, average the leeward/windward coef 
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
REAL,INTENT(IN)                :: PCHTCN  ! Convective heat transfer coefficient for natural conv. [W/(m2.K)]
REAL,INTENT(IN)                :: PVMOD   ! wind speed (m/s)
REAL                           :: PCHTC   ! Convective heat transfer coefficient [W/(m2.K)]
!
!*       0.2   Declarations of local variables
!
REAL, PARAMETER :: ZALEE = 2.86 ! coef for leeward facade
REAL, PARAMETER :: ZBLEE = 0.617 ! coef for leeward facade
!----------------------------------------------------------------------
!
!*       2.    COMPUTE THE CHTC
!              ----------------
!
PCHTC = SQRT(PCHTCN**2+(ZALEE*PVMOD**ZBLEE)**2)
!
!-------------------------------------------------------------------------------
!
END FUNCTION CHTC_SMOOTH_LEE_DOE_0D
