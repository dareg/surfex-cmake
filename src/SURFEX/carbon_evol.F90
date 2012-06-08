!     #########
    SUBROUTINE CARBON_EVOL(HRESPSL, HPHOTO, PTSTEP,                           &
                             PRHOA, PTG, PWG, PWFC, PWWILT, PWSAT, PSAND,       &
                             PRE25, PLAI, PRESP_BIOMASS_INST, PTURNOVER,        &
                             PLITTER, PLIGNIN_STRUC , PSOILCARB,                &
                             PRESP_AUTO, PRESP_ECO                              )  
!   ###############################################################
!!****  *CARBON EVOL*
!!
!!    PURPOSE
!!    -------
!!
!!    Diagnoses respiration carbon fluxes and performs the time evolution of 
!!    carbon pools in the case of 'CNT' option (ISBA-CC) 
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
!!      Gibelin et al. 2008, AFM
!!        Modelling energy and CO2 fluxes with an interactive vegetation land surface model -
!!        Evaluation at high and middle latitudes.
!!      
!!    AUTHOR
!!    ------
!!
!!	A.L. Gibelin       * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    22/06/09
!!      S.QUEGUINER 09/2011 Cas 'DEF'- condition si LAI=UNDEF->ZRESP_SOIL_TOT=0
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CO2V_PAR,       ONLY : XMC, XMCO2, XPCCO2
USE MODD_CSTS,           ONLY : XDAY, XTT
!USE MODD_DATA_COVER_PAR, ONLY : NVT_NO, NVT_ROCK, NVT_SNOW
USE MODD_SURF_PAR,       ONLY : XUNDEF
!
USE MODI_CONTROL_MOIST_FUNC
USE MODI_CONTROL_TEMP_FUNC
USE MODI_CARBON_LITTER
USE MODI_CARBON_SOIL

!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
CHARACTER(LEN=3), INTENT(IN) :: HRESPSL                ! Soil Respiration
!                                                      ! 'DEF' = Norman 1992
!                                                      ! 'PRM' = Rivalland PhD Thesis (2003)
!                                                      ! 'CNT' = CENTURY model (Gibelin 2008)
CHARACTER(LEN=3), INTENT(IN) :: HPHOTO                 ! type of photosynthesis
!
REAL, INTENT(IN)             :: PTSTEP                 ! time step
!
REAL, DIMENSION(:), INTENT(IN)        :: PRHOA         ! air density (kg/m3)
REAL, DIMENSION(:,:), INTENT(IN)      :: PTG           ! soil layer average temperatures (K)
REAL, DIMENSION(:,:), INTENT(IN)      :: PWG           ! soil liquid volumetric water content (m3/m3)
REAL, DIMENSION(:,:), INTENT(IN)      :: PWFC          ! field capacity profile (m3/m3)
REAL, DIMENSION(:,:), INTENT(IN)      :: PWWILT        ! wilting point profile (m3/m3)
REAL, DIMENSION(:,:), INTENT(IN)      :: PWSAT         ! porosity profile (m3/m3)
REAL, DIMENSION(:,:), INTENT(IN)      :: PSAND         ! profile of sand fraction
!
REAL,DIMENSION(:), INTENT(IN)         :: PRE25         ! Ecosystem respiration parameter (kg m-2 s-1)
REAL, DIMENSION(:), INTENT(IN)        :: PLAI          ! leaf area index
REAL, DIMENSION(:,:), INTENT(IN)      :: PRESP_BIOMASS_INST ! instantaneous respiration of biomass (kgCO2/kgair m/s)
REAL, DIMENSION(:,:), INTENT(IN)      :: PTURNOVER     ! biomass turnover going into litter
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PLITTER       ! litter pools
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PLIGNIN_STRUC ! L/C ratio in structural litter
REAL, DIMENSION(:,:), INTENT(INOUT)   :: PSOILCARB     ! soil carbon pools
REAL, DIMENSION(:), INTENT(OUT)       :: PRESP_AUTO    ! autotrophic respiration (kgCO2/kgair m/s)
REAL, DIMENSION(:), INTENT(OUT)       :: PRESP_ECO     ! total ecosystem respiration (kgCO2/kgair m/s)
!
!*      0.2    declarations of local variables
!
REAL,  DIMENSION(SIZE(PTG,1))    :: ZT2                ! soil temperature (C) 
!
REAL, DIMENSION(SIZE(PTG,1))     :: ZRESP_SOIL_TOT     ! total soil respiration (kgCO2/kgair m/s)
REAL, DIMENSION(SIZE(PTG,1))     :: ZRESP_AUTO_ABOVE   ! total above ground biomass respiration (kgCO2/kgair m/s)
REAL, DIMENSION(SIZE(PTG,1))     :: ZRESP_HETERO       ! total heterotrophic respiration (kgCO2/kgair m/s)
!
REAL, DIMENSION(SIZE(PSOILCARB,1),SIZE(PSOILCARB,2)) :: ZSOILCARBON_INPUT
!                  quantity of carbon going into carbon pools 
!                  from litter decomposition (gC/m2/day)
!
REAL, DIMENSION(SIZE(PSOILCARB,1)) :: ZRESP_HETERO_DAY_LITTER 
!                  litter heterotrophic respiration (gC/m2/day)
REAL, DIMENSION(SIZE(PSOILCARB,1)) :: ZRESP_HETERO_DAY_SOIL   
!                  soil heterotrophic respiration (gC/m2/day)
!
REAL, DIMENSION(SIZE(PLIGNIN_STRUC,1),SIZE(PLIGNIN_STRUC,2)) :: ZCONTROL_MOIST, &
                                                                  ZCONTROL_TEMP  
!                  ZCONTROL_MOIST = moisture control of heterotrophic respiration
!                  ZCONTROL_TEMP = temperature control of heterotrophic respiration
!
INTEGER                          :: JNBIOMASS
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-----------------------------------------------------------------
!
!*      1.     Preliminaries
!              -------------
!
IF (LHOOK) CALL DR_HOOK('CARBON_EVOL',0,ZHOOK_HANDLE)
ZRESP_SOIL_TOT(:)          = XUNDEF
ZRESP_AUTO_ABOVE(:)        = XUNDEF
ZRESP_HETERO(:)            = XUNDEF
ZSOILCARBON_INPUT(:,:)     = XUNDEF
ZRESP_HETERO_DAY_LITTER(:) = XUNDEF
ZRESP_HETERO_DAY_SOIL(:)   = XUNDEF
ZCONTROL_MOIST(:,:)        = XUNDEF
ZCONTROL_TEMP(:,:)         = XUNDEF
!
! convert soil temperature from K to C
!
ZT2(:)   = PTG(:,2) - XTT

!
!*      2.     Autotrophic respiration
!              -----------------------
!
ZRESP_AUTO_ABOVE(:)=0.

IF (HPHOTO=='NIT') THEN
  DO JNBIOMASS=1,2
    ZRESP_AUTO_ABOVE(:) = ZRESP_AUTO_ABOVE(:) + PRESP_BIOMASS_INST(:,JNBIOMASS)
  ENDDO
ELSE IF (HPHOTO=='NCB') THEN
  DO JNBIOMASS=1,3
    ZRESP_AUTO_ABOVE(:) = ZRESP_AUTO_ABOVE(:) + PRESP_BIOMASS_INST(:,JNBIOMASS)
  ENDDO
ELSE IF (HPHOTO=='AGS' .OR. HPHOTO=='AST' .OR. HPHOTO=='LAI' .OR. HPHOTO=='LST') THEN
  ZRESP_AUTO_ABOVE(:) = PRESP_BIOMASS_INST(:,1)
ENDIF

PRESP_AUTO(:)=SUM(PRESP_BIOMASS_INST,2)

!
!*      3.     Soil and ecosystem respiration
!              ------------------------------
!
!
IF (HRESPSL == 'DEF') THEN

  ! Soil respiration from Norman et al 1992 (kgCO2/kgair m/s)
  WHERE (PLAI == XUNDEF)
    ZRESP_SOIL_TOT(:) = 0.0
  ELSEWHERE
    ZRESP_SOIL_TOT(:) = 4.4E-8 / PRHOA * &
                       (13.5+5.4*PLAI(:))*PWG(:,1) * 2.**(0.1*(ZT2(:)-25.0))  
  ENDWHERE

  ! RESP_ECO is diagnosed from RESP_SOIL_TOT and RESP_AUTO_ABOVE

  PRESP_ECO(:) = ZRESP_SOIL_TOT(:) + ZRESP_AUTO_ABOVE(:)
  
ELSE IF (HRESPSL == 'PRM') THEN

  ! Ecosystem respiration : Q10 model following Albergel et al. 2009 for SMOSREX
  ! (kgCO2/kgair m/s)
  ! RESP_ECO is directly calculated by the parameterization

  PRESP_ECO(:) = PRE25(:)/PRHOA(:) * MIN(PWG(:,1)/PWFC(:,1),1.) * 2.**(0.1*(ZT2(:)-25.0))

ELSE IF (HRESPSL == 'CNT') THEN

  ! heterotrophic respiration following CENTURY model from Gibelin et al. 2008

! Calculations of ZCONTROL_TEMP and ZCONTROL_MOIST should depend on YISBA option.
! Following calculations are valid for '2-L' option.

  ZCONTROL_TEMP(:,1) = CONTROL_TEMP_FUNC (PTG(:,1))
  ZCONTROL_TEMP(:,2) = CONTROL_TEMP_FUNC (PTG(:,2))
  ZCONTROL_MOIST(:,1) = CONTROL_MOIST_FUNC (PWG(:,1),PWWILT(:,1),PWFC(:,1),PWSAT(:,1))
  ZCONTROL_MOIST(:,2) = CONTROL_MOIST_FUNC (PWG(:,2),PWWILT(:,2),PWFC(:,2),PWSAT(:,2))
  
  CALL CARBON_LITTER (PTSTEP,PTURNOVER,PLITTER,PLIGNIN_STRUC,                 &
                     ZCONTROL_TEMP,ZCONTROL_MOIST,                              &
                     ZRESP_HETERO_DAY_LITTER,ZSOILCARBON_INPUT)  

  CALL CARBON_SOIL (PTSTEP,PSAND(:,2),ZSOILCARBON_INPUT,                      &
                     ZCONTROL_TEMP,ZCONTROL_MOIST,PSOILCARB,ZRESP_HETERO_DAY_SOIL)  

  ! Transform units : gC/m2/day -> kgCO2/kgair m/s
  ZRESP_HETERO(:) = (ZRESP_HETERO_DAY_LITTER(:) + ZRESP_HETERO_DAY_SOIL(:)) &
                       * (XMCO2/XMC) / (1000. * PRHOA(:)* XDAY)  
  
  ! RESP_ECO is diagnosed from RESP_HETERO and RESP_AUTO

  PRESP_ECO(:) = ZRESP_HETERO(:) + PRESP_AUTO(:)
  
ENDIF
IF (LHOOK) CALL DR_HOOK('CARBON_EVOL',1,ZHOOK_HANDLE)

!
!-----------------------------------------------------------------
!
END SUBROUTINE CARBON_EVOL
