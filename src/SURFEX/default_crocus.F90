!     #########
      SUBROUTINE DEFAULT_CROCUS(HSNOWDRIFT,OSNOWDRIFT_SUBLIM,OSNOW_ABS_ZENITH,&
                 HSNOWMETAMO,HSNOWRAD,OSNOWSYTRON,HSNOWFALL,HSNOWCOND,HSNOWHOLD,&
                 HSNOWCOMP,HSNOWZREF)  
!     ########################################################################
!
!!****  *DEFAULT_ISBA* - routine to set default values for the configuration for Crocus
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      M. Lafaysse   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/2012
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
! Snowdrift scheme and selection of properties of falling snow  
CHARACTER(*), INTENT(OUT)     :: HSNOWDRIFT
! Logicals to activate / disactivate blowing snow sublimation  
LOGICAL, INTENT(OUT)          :: OSNOWDRIFT_SUBLIM
LOGICAL, INTENT(OUT)          :: OSNOW_ABS_ZENITH
! Logical to activate / disactivate Sytron                                          
LOGICAL, INTENT(OUT)          :: OSNOWSYTRON
!
! Snow metamorphism scheme and radiative transfer scheme <bber
 CHARACTER(*), INTENT(OUT) :: HSNOWMETAMO,HSNOWRAD,HSNOWFALL, HSNOWCOND, HSNOWHOLD, HSNOWCOMP, HSNOWZREF!added HSNOWFALL HSNOWHOLD and HSNOWCOND bber>
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!                                          
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DEFAULT_CROCUS',0,ZHOOK_HANDLE)
!
HSNOWDRIFT        = 'DFLT'
OSNOWDRIFT_SUBLIM = .FALSE.
OSNOW_ABS_ZENITH = .FALSE.
OSNOWSYTRON=.FALSE.
!
HSNOWMETAMO = 'B92'
HSNOWRAD    = 'B92'
!New multiphysics options Cluzet et al 2016
HSNOWFALL   = 'V12'
HSNOWCOND   = 'Y81'
HSNOWHOLD   = 'B92'
HSNOWCOMP   = 'B92'
HSNOWZREF   = 'CST'
!
IF (LHOOK) CALL DR_HOOK('DEFAULT_CROCUS',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE DEFAULT_CROCUS
