!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!###################
MODULE MODE_CONV_DOE
!###################
!
!!****  *MODE_CONV_DOE* -
!!
!!    PURPOSE
!!    -------
!      
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!       NONE          
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!    G. Pigeon
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/08/12 
!
!--------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!--------------------------------------------------------------------------------
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
CONTAINS
!----------------------------
!#############################################
FUNCTION CHTC_VERT_DOE(PTS, PTA) RESULT(PCHTC)
!#############################################
!
!!****  *CHTC_VERT_DOE* - 
!!
!!    PURPOSE
!!    -------
!     function to compute convective surface coefficient for
!     a vertical surface from surface temperature and air
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
!
REAL, DIMENSION(:), INTENT(IN)                :: PTS     ! Surface temperature (Kelvin)
REAL, DIMENSION(:), INTENT(IN)                :: PTA     ! Air temperature (Kelvin)
REAL, DIMENSION(SIZE(PTS))                    :: PCHTC   ! Convective heat transfer coefficient [W/(m2.K)]
!
INTEGER :: J
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_VERT_DOE',0,ZHOOK_HANDLE)
!
!*       1.    COMPUTE THE CHTC
!              ----------------
!
DO J=1,SIZE(PCHTC)
  PCHTC(J) = CHTC_VERT_DOE_0D(PTS(J), PTA(J))
ENDDO
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_VERT_DOE',1,ZHOOK_HANDLE)
!
END FUNCTION CHTC_VERT_DOE
!#########################
!
!#############################################
FUNCTION CHTC_UP_DOE(PTS, PTA) RESULT(PCHTC)
!#############################################
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
REAL, DIMENSION(:), INTENT(IN)                :: PTS     ! Surface temperature (Kelvin)
REAL, DIMENSION(:), INTENT(IN)                :: PTA     ! Air temperature (Kelvin)
REAL, DIMENSION(SIZE(PTS))                    :: PCHTC   ! Convective heat transfer coefficient [W/(m2.K)]
!
INTEGER :: J
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_UP_DOE',0,ZHOOK_HANDLE)
!
!*       1.    COMPUTE THE CHTC
!              ----------------
!
DO J=1,SIZE(PCHTC)
  PCHTC(J) = CHTC_UP_DOE_0D(PTS(J), PTA(J))
ENDDO
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_UP_DOE',1,ZHOOK_HANDLE)
!
END FUNCTION CHTC_UP_DOE
!#######################
!
!#############################################
FUNCTION CHTC_DOWN_DOE(PTS, PTA) RESULT(PCHTC)
!#############################################
!
!!****  *CHTC_DOWN_DOE* - 
!!
!!    PURPOSE
!!    -------
!     function to compute convective surface coefficient for
!     an downward surface from surface temperature and air
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
REAL, DIMENSION(:), INTENT(IN)                :: PTS     ! Surface temperature (Kelvin)
REAL, DIMENSION(:), INTENT(IN)                :: PTA     ! Air temperature (Kelvin)
REAL, DIMENSION(SIZE(PTS))                    :: PCHTC   ! Convective heat transfer coefficient [W/(m2.K)]
!
INTEGER :: J
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_DOWN_DOE',0,ZHOOK_HANDLE)
!
!*       1.    COMPUTE THE CHTC
!              ----------------
!
DO J=1,SIZE(PCHTC)
  PCHTC(J) = CHTC_DOWN_DOE_0D(PTS(J), PTA(J))
ENDDO
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_DOWN_DOE',1,ZHOOK_HANDLE)
!
END FUNCTION CHTC_DOWN_DOE
!#########################
!
!#############################################
FUNCTION CHTC_SMOOTH_LEE_DOE(PCHTCN, PVMOD) RESULT(PCHTC)
!#############################################
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
REAL, DIMENSION(:), INTENT(IN)                :: PCHTCN  ! Convective heat transfer coefficient for natural conv. [W/(m2.K)]
REAL, DIMENSION(:), INTENT(IN)                :: PVMOD   ! wind speed (m/s)
REAL, DIMENSION(SIZE(PCHTCN))                 :: PCHTC   ! Convective heat transfer coefficient [W/(m2.K)]
!
!*       0.2   Declarations of local variables
!
INTEGER :: J
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_SMOOTH_LEE_DOE',0,ZHOOK_HANDLE)
!
!*       2.    COMPUTE THE CHTC
!              ----------------
!
DO J=1,SIZE(PCHTC)
  PCHTC(J) = CHTC_SMOOTH_LEE_DOE_0D(PCHTCN(J), PVMOD(J))
ENDDO
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_SMOOTH_LEE_DOE',1,ZHOOK_HANDLE)
!
END FUNCTION CHTC_SMOOTH_LEE_DOE
!#########################

!#############################################
FUNCTION CHTC_SMOOTH_WIND_DOE(PCHTCN, PVMOD) RESULT(PCHTC)
!#############################################
!
!!****  *CHTC_SMOOTH_WIND_DOE* - 
!!
!!    PURPOSE
!!    -------
!     function to compute convective surface coefficient for
!     a windward smooth surface from the natural convection coef and the
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
REAL, DIMENSION(:), INTENT(IN)                :: PCHTCN  ! Convective heat transfer coefficient for natural conv. [W/(m2.K)]
REAL, DIMENSION(:), INTENT(IN)                :: PVMOD   ! wind speed (m/s)
REAL, DIMENSION(SIZE(PCHTCN))                 :: PCHTC   ! Convective heat transfer coefficient [W/(m2.K)]
!
!*       0.2   Declarations of local variables
!
INTEGER :: J
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_SMOOTH_WIND_DOE',0,ZHOOK_HANDLE)
!
!*       2.    COMPUTE THE CHTC
!              ----------------
!
DO J=1,SIZE(PCHTC)
  PCHTC(J) = CHTC_SMOOTH_WIND_DOE_0D(PCHTCN(J), PVMOD(J))
ENDDO
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_SMOOTH_WIND_DOE',1,ZHOOK_HANDLE)
!
END FUNCTION CHTC_SMOOTH_WIND_DOE
!#########################
!#############################################
FUNCTION CHTC_ROUGH_DOE(PCHTCN, PCHTCS, PROUGH) RESULT(PCHTC)
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
REAL, DIMENSION(:), INTENT(IN)                :: PCHTCN  ! Convective heat transfer coefficient for natural conv. [W/(m2.K)]
REAL, DIMENSION(:), INTENT(IN)                :: PCHTCS  ! Convective heat transfer coefficient for a smooth surface [W/(m2.K)]
REAL, DIMENSION(:), INTENT(IN)                :: PROUGH  ! Convective heat transfer coefficient for a smooth surface [W/(m2.K)]
REAL, DIMENSION(SIZE(PCHTCN))                 :: PCHTC   ! Convective heat transfer coefficient [W/(m2.K)]
!
!*       0.2   Declarations of local variables
!
INTEGER :: J
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!----------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_ROUGH_DOE',0,ZHOOK_HANDLE)
!
!*       1.    COMPUTE THE CHTC
!              ----------------
!
DO J=1,SIZE(PCHTC)
  PCHTC(J) = CHTC_ROUGH_DOE_0D(PCHTCN(J), PCHTCS(J), PROUGH(J))
ENDDO
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('MODE_CONV_DOE:CHTC_ROUGH_DOE',1,ZHOOK_HANDLE)
!
END FUNCTION CHTC_ROUGH_DOE

INCLUDE "chtc_up_doe_0d.func.h"
INCLUDE "chtc_smooth_wind_doe_0d.func.h"
INCLUDE "chtc_smooth_lee_doe_0d.func.h"
INCLUDE "chtc_rough_doe_0d.func.h"
INCLUDE "chtc_vert_doe_0d.func.h"
INCLUDE "chtc_down_doe_0d.func.h"

END MODULE MODE_CONV_DOE
