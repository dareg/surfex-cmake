!     #########
SUBROUTINE CARBON_SOIL (PTSTEP, PSAND,                                        &
                          PSOILCARBON_INPUT, PCONTROL_TEMP, PCONTROL_MOIST,     &
                          PSOILCARB, PRESP_HETERO_SOIL)  

!   ###############################################################
!!**  CARBON_SOIL 
!!
!!    PURPOSE
!!    -------
!!    Calculates soil carbon pools evolution.
!!
!!**  METHOD
!!    ------
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
!!      Parton et al., Biogeochemestry, 1988
!!      Krinner et al., Global Biochemical Cycles, 2005
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
USE MODD_CO2V_PAR,       ONLY : XTAU_SOILCARB
USE MODD_CSTS,           ONLY : XDAY
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE


!*       0.1 input

! time step in s
REAL, INTENT(IN)                                                  :: PTSTEP
! sand fraction (between 0 and 1)
REAL, DIMENSION(:), INTENT(IN)                                    :: PSAND
! quantity of carbon going into carbon pools from litter decomposition
!   (gC/m**2/day)
REAL, DIMENSION(:,:), INTENT(IN)                                  :: PSOILCARBON_INPUT
! temperature control of heterotrophic respiration
REAL, DIMENSION(:,:), INTENT(IN)                                  :: PCONTROL_TEMP
! moisture control of heterotrophic respiration
REAL, DIMENSION(:,:), INTENT(IN)                                  :: PCONTROL_MOIST

!*       0.2 modified fields

! carbon pool: active, slow, or passive (gC/m**2)
REAL, DIMENSION(:,:), INTENT(INOUT)                               :: PSOILCARB

!*       0.3 output

! soil heterotrophic respiration (in gC/day/m**2)
REAL, DIMENSION(:), INTENT(OUT)                                   :: PRESP_HETERO_SOIL

!*       0.4 local

! time step in days
REAL                                                              :: ZDT
! flux fractions within carbon pools
REAL, DIMENSION(SIZE(PSOILCARB,1),SIZE(PSOILCARB,2),SIZE(PSOILCARB,2)) :: ZFRAC_CARB
! fraction of carbon flux which goes into heterotrophic respiration
REAL, DIMENSION(SIZE(PSOILCARB,1),SIZE(PSOILCARB,2))                   :: ZFRAC_RESP
! total flux out of carbon pools (gC/m**2)
REAL, DIMENSION(SIZE(PSOILCARB,1),SIZE(PSOILCARB,2))                   :: ZFLUXTOT
! fluxes between carbon pools (gC/m**2)
REAL, DIMENSION(SIZE(PSOILCARB,1),SIZE(PSOILCARB,2),SIZE(PSOILCARB,2)) :: ZFLUX
! dimensions
INTEGER                                                           :: INSOILCARB
! indices
INTEGER                                                           :: K,KK
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! correspondence between array indices and litter levels
! LT_ABOVE = 1
! LT_BELOW = 2
! correspondence between array indices and soil carbon pools
! SL_ACTIVE = 1
! SL_SLOW = 2
! SL_PASSIVE = 3
!-------------------------------------------------------------------------------

!
!*       1 Initialisations
!

!
!*       1.1 dimensions
!
IF (LHOOK) CALL DR_HOOK('CARBON_SOIL',0,ZHOOK_HANDLE)
INSOILCARB = SIZE(PSOILCARB,2)

!
!*       1.2 get soil "constants"
!

!*       1.2.1 flux fractions between carbon pools: depend on soil texture, 
!       recalculated each time

!*       1.2.1.1 from active pool: depends on soil texture

ZFRAC_CARB(:,1,1) = 0.0
ZFRAC_CARB(:,1,3) = 0.004
ZFRAC_CARB(:,1,2) = 1. - ( .85 - .68 * (1.-PSAND(:)) ) - ZFRAC_CARB(:,1,3)

!*       1.2.1.2 from slow pool

ZFRAC_CARB(:,2,2) = .0
ZFRAC_CARB(:,2,1) = .42
ZFRAC_CARB(:,2,3) = .03

!*       1.2.1.3 from passive pool

ZFRAC_CARB(:,3,3) = .0
ZFRAC_CARB(:,3,1) = .45
ZFRAC_CARB(:,3,2) = .0


!
!*       1.3 set output to zero
!

PRESP_HETERO_SOIL(:) = 0.0

!
!*       2 input into carbon pools
!

ZDT = PTSTEP/XDAY

PSOILCARB(:,:) = PSOILCARB(:,:) + PSOILCARBON_INPUT(:,:) * ZDT

!
!*       3 fluxes within carbon reservoirs + respiration
!

!
!*       3.1 determine fraction of flux that is respiration
!     diagonal elements of frac_carb are zero

ZFRAC_RESP(:,:) = 1. - ZFRAC_CARB(:,:,1) - ZFRAC_CARB(:,:,2) - &
                   ZFRAC_CARB(:,:,3)   

!
!*       3.2 calculate fluxes
!


!*       3.2.1 flux out of pools

DO K = 1, INSOILCARB

  ! determine total flux out of pool

  ZFLUXTOT(:,K) = PTSTEP/XTAU_SOILCARB(K) * PSOILCARB(:,K) * &
                   PCONTROL_MOIST(:,2) * PCONTROL_TEMP(:,2)  

  IF ( K .EQ. 1 ) THEN
    ZFLUXTOT(:,K) = ZFLUXTOT(:,K) * ( 1. - .75 * (1.-PSAND(:)) )
  ENDIF

  ! decrease this carbon pool

  PSOILCARB(:,K) = PSOILCARB(:,K) - ZFLUXTOT(:,K)

  ! fluxes towards the other pools (k -> kk)

  DO KK = 1, INSOILCARB
    ZFLUX(:,K,KK) = ZFRAC_CARB(:,K,KK) * ZFLUXTOT(:,K)
  ENDDO

ENDDO

!*       3.2.2 respiration

PRESP_HETERO_SOIL(:) = &
       ( ZFRAC_RESP(:,1) * ZFLUXTOT(:,1) + &
         ZFRAC_RESP(:,2) * ZFLUXTOT(:,2) + &
         ZFRAC_RESP(:,3) * ZFLUXTOT(:,3)  ) / ZDT  

!*       3.2.3 add fluxes to active, slow, and passive pools

DO K = 1, INSOILCARB
  PSOILCARB(:,K) = PSOILCARB(:,K) + &
                    ZFLUX(:,1,K) + ZFLUX(:,3,K) + ZFLUX(:,2,K)  
ENDDO
IF (LHOOK) CALL DR_HOOK('CARBON_SOIL',1,ZHOOK_HANDLE)
!
END SUBROUTINE CARBON_SOIL
