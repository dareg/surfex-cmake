!     #########
      SUBROUTINE DEFAULT_CROCUS(OSNOWDRIFT,OSNOWDRIFT_SUBLIM,OSNOW_ABS_ZENITH,&
                 HSNOWMETAMO,HSNOWRAD,OSNOWSYTRON, &
		 OSNOWCOMPACT_BOOL, OSNOWMAK_BOOL, &
		 OPRODSNOWMAK, OSNOWMAK_PROP, OSNOWTILLER, OSELF_PROD)  
!
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
! Logicals to activate / disactivate snowdrift                                          
LOGICAL, INTENT(OUT)          :: OSNOWDRIFT
LOGICAL, INTENT(OUT)          :: OSNOWDRIFT_SUBLIM
LOGICAL, INTENT(OUT)          :: OSNOW_ABS_ZENITH
! Logical to activate / disactivate Sytron                                          
LOGICAL, INTENT(OUT)          :: OSNOWSYTRON
!
! Snow metamorphism scheme and radiative transfer scheme
 CHARACTER(*), INTENT(OUT) :: HSNOWMETAMO,HSNOWRAD
!
! Snowmaking option p.spandre 19/11/2013                                          
!
LOGICAL, INTENT(OUT)         	  :: OSNOWCOMPACT_BOOL
LOGICAL, DIMENSION(:), INTENT(OUT):: OPRODSNOWMAK
LOGICAL, INTENT(OUT)        	  :: OSNOWMAK_BOOL
LOGICAL, INTENT(OUT)     	  :: OSNOWMAK_PROP
LOGICAL, INTENT(OUT)       	  :: OSNOWTILLER
LOGICAL, INTENT(OUT)        	  :: OSELF_PROD
!
!*       0.2   Declarations of local variables
!              -------------------------------
!                                          
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DEFAULT_CROCUS',0,ZHOOK_HANDLE)
!
OSNOWDRIFT        = .TRUE.
OSNOWDRIFT_SUBLIM = .FALSE.
OSNOW_ABS_ZENITH = .FALSE.
OSNOWSYTRON=.FALSE.
!
! Snowmaking and Grooming default option and pressure default value | P.Spandre     20160211                                    
OSNOWCOMPACT_BOOL=.FALSE.                                    
OSNOWMAK_BOOL =.FALSE.
OPRODSNOWMAK = .FALSE.
OSNOWMAK_PROP =.FALSE.
OSNOWTILLER =.FALSE.
OSELF_PROD =.FALSE.
!
! 
HSNOWMETAMO = 'B92'
HSNOWRAD    = 'B92'
!
IF (LHOOK) CALL DR_HOOK('DEFAULT_CROCUS',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE DEFAULT_CROCUS
