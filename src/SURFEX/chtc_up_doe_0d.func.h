!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!###################
FUNCTION CHTC_UP_DOE_0D(PTS, PTA) RESULT(PCHTC)
!
!!****  *CHTC_UP_DOE* - 
!!
!!    PURPOSE
!!    -------
!     function to compute convective surface coefficient for
!     an upward surface from surface temperature and air
!     temperature
!
!!**  METHOD
!!    ------
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
REAL, INTENT(IN)                :: PTS     ! Surface temperature (Kelvin)
REAL, INTENT(IN)                :: PTA     ! Air temperature (Kelvin)
REAL                            :: PCHTC   ! Convective heat transfer coefficient [W/(m2.K)]
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE THE CHTC
!              ----------------
!
IF (PTA < PTS) THEN
   PCHTC  = 9.482 * (PTS - PTA)**(1./3.) / (7.283 - 1)
ELSE 
   PCHTC  = 1.810 * (PTA - PTS)**(1./3.) / (1.382 + 1)  
ENDIF
!
END FUNCTION CHTC_UP_DOE_0D
