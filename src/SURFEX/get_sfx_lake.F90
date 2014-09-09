!     #########
      SUBROUTINE GET_SFX_LAKE(PLAKE_EVAP,PLAKE_RAIN,PLAKE_SNOW)  
!     ############################################################################
!
!!****  *GET_SFX_LAKE* - routine to get some variables from surfex to
!                        a oceanic general circulation model
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
!!	B. Decharme      *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2013
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE MODD_SURF_ATM_n, ONLY : NR_WATER, NSIZE_WATER
!
USE MODD_FLAKE_n,    ONLY : XCPL_FLAKE_EVAP, &
                            XCPL_FLAKE_RAIN, &
                            XCPL_FLAKE_SNOW  
!
USE MODI_UNPACK_SAME_RANK
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
REAL, DIMENSION(:), INTENT(OUT) :: PLAKE_EVAP  ! Cumulated Evaporation             (kg/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PLAKE_RAIN  ! Cumulated Rainfall rate           (kg/m2)
REAL, DIMENSION(:), INTENT(OUT) :: PLAKE_SNOW  ! Cumulated Snowfall rate           (kg/m2)
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GET_SFX_LAKE',0,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
!*       1.0   Initialization
!              --------------
!
PLAKE_EVAP (:) = XUNDEF
PLAKE_RAIN (:) = XUNDEF
PLAKE_SNOW (:) = XUNDEF
!
!*       2.0   Get variable over lake
!              ----------------------
!
IF(NSIZE_WATER>0)THEN
!
  CALL UNPACK_SAME_RANK(NR_WATER(:),XCPL_FLAKE_EVAP(:),PLAKE_EVAP(:),XUNDEF)
  CALL UNPACK_SAME_RANK(NR_WATER(:),XCPL_FLAKE_RAIN(:),PLAKE_RAIN(:),XUNDEF)
  CALL UNPACK_SAME_RANK(NR_WATER(:),XCPL_FLAKE_SNOW(:),PLAKE_SNOW(:),XUNDEF)
  XCPL_FLAKE_EVAP(:) = 0.0
  XCPL_FLAKE_RAIN(:) = 0.0
  XCPL_FLAKE_SNOW(:) = 0.0
!
ENDIF
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('GET_SFX_LAKE',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE GET_SFX_LAKE
