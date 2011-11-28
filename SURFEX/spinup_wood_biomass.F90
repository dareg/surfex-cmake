!     #########
SUBROUTINE SPINUP_WOOD_BIOMASS( PVEGTYPE,PBIOMASS,PINCREASE,PTAU_WOOD )
!   ###############################################################
!!****  *SPINUP_WOOD_BIOMASS* - spinup of wood biomass
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!     Calculations for spinup of woody biomass.
!!     These calculations are derived from nitro_biomass_carbon.f90.
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
!!      Gibelin et al. 2008, AFM
!!        Modelling energy and CO2 fluxes with an interactive vegetation land surface model -
!!        Evaluation at high and middle latitudes.
!!      
!!    AUTHOR
!!    ------
!!
!!	A.-L. Gibelin           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    02/05/09
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE, NVT_NO, NVT_ROCK, NVT_SNOW, &
  NVT_TREE, NVT_CONI, NVT_EVER, NVT_C3, NVT_C4, NVT_IRR, NVT_GRAS, NVT_TROG, NVT_PARK  
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,   DIMENSION(:), INTENT(IN) :: PTAU_WOOD          ! residence time in wood (s)
!
REAL, DIMENSION(:,:), INTENT(IN) :: PVEGTYPE           ! fraction of each vegetation
!
REAL,   DIMENSION(:,:), INTENT(INOUT) :: PBIOMASS      ! biomass reservoirs
!
REAL,   DIMENSION(:,:), INTENT(IN)    :: PINCREASE     ! biomass increase
!
!*      0.2    declarations of local variables
!
REAL,    DIMENSION(SIZE(PBIOMASS,1),SIZE(PBIOMASS,2))  :: ZBIOMASS   
!                                   temporary biomass reservoirs
REAL,    DIMENSION(SIZE(PBIOMASS,1),SIZE(PBIOMASS,2))  :: ZDECLINE   
!                                   biomass decline (storage+mortality)
REAL, PARAMETER                                        :: ZDAY = 86400.  
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
! 1 - Initialisations
! -------------------
!
IF (LHOOK) CALL DR_HOOK('SPINUP_WOOD_BIOMASS',0,ZHOOK_HANDLE)
ZDECLINE(:,:) = 0.0
ZBIOMASS(:,:) = 0.0

ZBIOMASS(:,5) = PBIOMASS(:,5)
ZBIOMASS(:,6) = PBIOMASS(:,6)

!
! 2 - Wood biomass decline during one day
! ---------------------------------------
!
WHERE (PVEGTYPE(:,NVT_TREE)+PVEGTYPE(:,NVT_CONI)+PVEGTYPE(:,NVT_EVER) > 0.5)
!
  ZDECLINE(:,5) = ZBIOMASS(:,5)*(1.0-EXP(-1.0/(PTAU_WOOD(:)/ZDAY)))
  ZDECLINE(:,6) = ZBIOMASS(:,6)*(1.0-EXP(-1.0/(PTAU_WOOD(:)/ZDAY)))
!
ENDWHERE
           
!
! 3 - Evolution of woody reservoirs
! ---------------------------------
!
WHERE (PVEGTYPE(:,NVT_TREE)+PVEGTYPE(:,NVT_CONI)+PVEGTYPE(:,NVT_EVER) > 0.5)
!
  ZBIOMASS(:,5) = ZBIOMASS(:,5) + PINCREASE(:,5) - ZDECLINE(:,5)
  ZBIOMASS(:,6) = ZBIOMASS(:,6) + PINCREASE(:,6) - ZDECLINE(:,6)
!
END WHERE

!
! 4 - New values of woody biomass
! -------------------------------
!
  PBIOMASS(:,5) = ZBIOMASS(:,5)
  PBIOMASS(:,6) = ZBIOMASS(:,6)
IF (LHOOK) CALL DR_HOOK('SPINUP_WOOD_BIOMASS',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE SPINUP_WOOD_BIOMASS
