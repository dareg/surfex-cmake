!   #############
FUNCTION CONTROL_MOIST_FUNC (PMOIST_IN,PWWILT,PWFC,PWSAT) RESULT (PMOISTFUNC_RESULT)

!   ###############################################################
!!**   CONTROL_MOIST_FUNC
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!     Moisture control factor for decomposition.
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      
!!    none
!!
!!    REFERENCE
!!    ---------
!!
!!      Krinner et al., Global Biochemical Cycles, 2005
!!      Modified for Wfc < W < Wsat following Probert et al., Agricultural Systems, 1998
!!      Gibelin et al. 2008, AFM
!!      
!!    AUTHOR
!!    ------
!!
!!	A.-L. Gibelin           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    23/06/09
!!      
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE

!*       0.1 input

! soil water content (m3/m3)
REAL, DIMENSION(:), INTENT(IN)                          :: PMOIST_IN
! wilting point (m3/m3)
REAL, DIMENSION(:), INTENT(IN)                          :: PWWILT
! field capacity (m3/m3)
REAL, DIMENSION(:), INTENT(IN)                          :: PWFC
! porosity (m3/m3)
REAL, DIMENSION(:), INTENT(IN)                          :: PWSAT

!*       0.2 result

! moisture control factor
REAL, DIMENSION(SIZE(PMOIST_IN))                        :: PMOISTFUNC_RESULT

!*       0.3 local

! relative humidity (-)
REAL, DIMENSION(SIZE(PMOIST_IN))                        :: ZMOIST_REL
! relative humidity (-)
REAL                                                    :: ZSEUIL
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!
!*       1 Initialisations
!
IF (LHOOK) CALL DR_HOOK('CONTROL_MOIST_FUNC',0,ZHOOK_HANDLE)
ZMOIST_REL(:) = 0.
ZSEUIL = 0.05

!
!*       2 Calculates moisture control factor
!
WHERE (PMOIST_IN(:) .LE. PWWILT(:))
  PMOISTFUNC_RESULT(:) = ZSEUIL
ELSEWHERE (PMOIST_IN(:) .GT. PWWILT(:) .AND. PMOIST_IN(:) .LE. PWFC(:))
  ZMOIST_REL(:) = (PMOIST_IN(:)-PWWILT(:)) / (PWFC(:)-PWWILT(:))
  ZMOIST_REL(:) = MIN(1.,MAX(0.,ZMOIST_REL(:)))

  PMOISTFUNC_RESULT(:) = -1.1 * ZMOIST_REL(:) * ZMOIST_REL(:) + 2.4 * ZMOIST_REL(:) - 0.29
  PMOISTFUNC_RESULT(:) = MAX( ZSEUIL, MIN( 1., PMOISTFUNC_RESULT(:) ) )
ELSEWHERE (PMOIST_IN(:) .GT. PWFC(:) .AND. PMOIST_IN(:) .LE. PWSAT(:))
  PMOISTFUNC_RESULT(:) = 1./(2.*(PWFC(:)-PWSAT(:)))*PMOIST_IN(:) + &
                           (PWFC(:)-2.*PWSAT(:))/(2.*(PWFC(:)-PWSAT(:)))  
ELSEWHERE
  PMOISTFUNC_RESULT(:) = 0.5
ENDWHERE
IF (LHOOK) CALL DR_HOOK('CONTROL_MOIST_FUNC',1,ZHOOK_HANDLE)

END FUNCTION CONTROL_MOIST_FUNC
