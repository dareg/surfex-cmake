!     #########
SUBROUTINE NITRO_CARBON_DECLINE(OMASK,PTSTEP,HRESPSL,PBSLAI_NITRO,PSEFOLD,PGMES,PANMAX,PANDAY,  &
       PLAT,PLAIMIN,PVEGTYPE,PTG,PANFM,PLAI,PBIOMASS,PRESP_BIOMASS,PBIOMASS_LEAF,  &
       PTAU_WOOD,PINCREASE,PTURNOVER )  
!
!   ###############################################################
!!**  NITRO_CARBON_DECLINE 
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!     Calvet and Soussana (2001) and Gibelin et al. (2006) for nitrogen dilution.
!!     Gibelin et al. (2008) : New biomass reservoirs, and new method for allocation, 
!!     mortality and respiration.
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
!! Calvet and Soussana (2001), "Modelling CO2-enrichment effects using an
!! interactive vegetation SVAT scheme", Agricultural and Forest Meteorology, Vol. 108
!! pp. 129-152
!! Gibelin et al. (2008), "Modelling energy and CO2 fluxes with an interactive vegetation 
!! land surface model - Evaluation at high and middle latitudes", 
!! Agricultural and Forest Meteorology, Vol. 148 , pp. 1611-1628
!!      
!!    AUTHOR
!!    ------
!!
!!	A.-L. Gibelin           * Meteo-France *
!!      (following Belair)
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    23/09/05
!!      A.L. Gibelin 04/2009 :  adaptation to SURFEX environment
!!      A.   Barbu   01/2011 :  modification of active biomass,leaf reservoir (see nitro_decline.f90)
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,           ONLY : XPI, XDAY, XTT
USE MODD_CO2V_PAR,       ONLY : XPCCO2, XCC_NIT, XCA_NIT, XRESPFACTOR_NIT, &
                                  XCOEFF_MAINT_RESP_ZERO, XSLOPE_MAINT_RESP, XMC, XMCO2, XPCCO2  
USE MODD_DATA_COVER_PAR, ONLY : NVT_TREE, NVT_EVER, NVT_CONI
!
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
LOGICAL,DIMENSION(:), INTENT(IN) :: OMASK            ! where solar time since midnight (s)
!
REAL,                 INTENT(IN) :: PTSTEP           ! time step
CHARACTER(LEN=3),     INTENT(IN) :: HRESPSL          ! Soil Respiration
!                                                    ! 'DEF' = Norman 1992
!                                                    ! 'PRM' = Rivalland PhD Thesis (2003)
!                                                    ! 'CNT' = CENTURY model (Gibelin 2008)
!
REAL,   DIMENSION(:), INTENT(IN) :: PBSLAI_NITRO     ! ratio of biomass to LAI
REAL,   DIMENSION(:), INTENT(IN) :: PSEFOLD          ! e-folding time for senescence (s)
REAL,   DIMENSION(:), INTENT(IN) :: PTAU_WOOD        ! residence time in wood (s)
REAL,   DIMENSION(:), INTENT(IN) :: PGMES            ! mesophyll conductance (m s-1)
REAL,   DIMENSION(:), INTENT(IN) :: PANMAX           ! maximum photosynthesis rate
REAL,   DIMENSION(:), INTENT(IN) :: PLAT             ! latitude of each grid point
REAL,   DIMENSION(:), INTENT(IN) :: PLAIMIN          ! minimum LAI
REAL, DIMENSION(:,:), INTENT(IN) :: PVEGTYPE         ! fraction of each vegetation
REAL, DIMENSION(:,:), INTENT(IN) :: PTG              ! soil layer average temperatures (K)
REAL,   DIMENSION(:), INTENT(IN) :: PANDAY           ! daily net CO2 accumulation
!
REAL,   DIMENSION(:), INTENT(INOUT) :: PANFM         ! maximum leaf assimilation
REAL,   DIMENSION(:), INTENT(INOUT) :: PLAI          ! leaf area index (LAI) 
REAL, DIMENSION(:,:), INTENT(INOUT) :: PBIOMASS      ! biomass reservoirs
REAL, DIMENSION(:,:), INTENT(INOUT) :: PRESP_BIOMASS ! cumulated daily biomass respiration
!
REAL,   DIMENSION(:), INTENT(OUT)   :: PBIOMASS_LEAF ! temporary leaf biomass
REAL, DIMENSION(:,:), INTENT(OUT)   :: PINCREASE     ! increment of biomass
REAL, DIMENSION(:,:), INTENT(OUT)   :: PTURNOVER     ! biomass turnover going into litter
!
!
!*      0.2    declarations of local variables
!
REAL                            :: ZCC_NITRO           ! c coefficient for nitrogen dilution
REAL                            :: ZBIOMASST_LIM       ! threshold value of ZBIOMASST 
!                                                        in nitrogen dilution theory
REAL                            :: ZBMCOEF
REAL,    DIMENSION(SIZE(PLAI))  :: ZXSEFOLD            ! e-folding time for senescence 
!                                                        corrected (days)
REAL,    DIMENSION(SIZE(PLAI))  :: ZLAIB_NITRO         ! LAI correction parameter used 
!                                                        in sefold calculation
!
REAL,    DIMENSION(SIZE(PLAI))  :: ZASSIM              ! assimilation
REAL,    DIMENSION(SIZE(PLAI))  :: ZBIOMASST           ! leaf + active structural biomass
!
REAL,    DIMENSION(SIZE(PLAI),SIZE(PBIOMASS,2))  :: ZBIOMASS      ! temporary biomass reservoirs
REAL,    DIMENSION(SIZE(PLAI),SIZE(PBIOMASS,2))  :: ZDECLINE      ! biomass decline (storage+mortality)
REAL,    DIMENSION(SIZE(PLAI),SIZE(PBIOMASS,2))  :: ZSTORAGE      ! storage (part of decline)
REAL,    DIMENSION(SIZE(PLAI))                   :: ZMORT_LEAF    ! leaf mortality
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! correspondence between array indices and biomass compartments
! LEAF = 1
! STRUCT_ACT = 2
! STRUCT_PAS = 3
! STRUCT_BELOW = 4
! WOOD_ABOVE = 5
! WOOD_BELOW = 6
!
!-------------------------------------------------------------------------------
!
! 1 - Initialisations
! -------------------
!
IF (LHOOK) CALL DR_HOOK('NITRO_CARBON_DECLINE',0,ZHOOK_HANDLE)
ZCC_NITRO           = 0.0
ZBIOMASST_LIM       = 0.0
ZXSEFOLD(:)         = 0.0
ZLAIB_NITRO(:)      = 0.0
ZBIOMASST(:)        = 0.0
ZASSIM(:)           = 0.0
ZBIOMASS(:,:)       = 0.0
ZDECLINE(:,:)       = 0.0
ZSTORAGE(:,:)       = 0.0
ZMORT_LEAF(:)       = 0.0
!---------------------------------------------------
!
ZBMCOEF     = XMC/(XMCO2*XPCCO2)
!
!-----------------------------------------------------------------
!avoid possible but unlikely negative values for biomass:        
!
PBIOMASS(:,1) = MAX(PBIOMASS(:,1),0.0)
!
! current leaf biomass value:
!
PBIOMASS_LEAF(:) = PBIOMASS(:,1)
!
!-------------------------------------------------------------------------------
!
! Once a day (at midnight),repartition of net assimilation and mortality 
! into different biomass compartments.
!
! 2 - Evolution of leaf biomass and senescence calculations
! ---------------------------------------------------------
!
! coef c for biomass in kg/m2
!
  ZCC_NITRO = XCC_NIT/10.**XCA_NIT
!
! threshold value for leaf biomass and total above ground biomass in nitrogen
! dilution theory
!
  ZBIOMASST_LIM = ZCC_NITRO**(1.0/XCA_NIT)
!
WHERE( OMASK(:) )
!
!
! LAI correction for shadow effect
!
  ZLAIB_NITRO(:) = MAX( 5.76-0.64*ATAN(ABS(PLAT(:))*XPI/180.),3.8 )
!
! leaf life expectancy
!
  ZXSEFOLD(:) = PSEFOLD(:)*MIN(1.0, PANFM(:)/PANMAX(:))*   &
         MAX(((PGMES(:)*1000.)**0.321)*PLAI(:)/ZLAIB_NITRO(:),1.)/XDAY  
!
! avoid possible but unlikely division by zero
!
  ZXSEFOLD(:) = MAX(1.0E-8,ZXSEFOLD(:))
!
! limitation of leaf life expectancy
!
! OLD   ZXSEFOLD(:) = MAX(5.,ZXSEFOLD(:))
! Bug Seb: Following Marita's work limitation of the senesence
! HERE
  ZXSEFOLD(:) = MAX(PSEFOLD(:)/XDAY/10.0,ZXSEFOLD(:))
!
! senesence of active biomass
!
  ZDECLINE(:,1) = MIN(PBIOMASS_LEAF(:)-PLAIMIN(:)*PBSLAI_NITRO(:), &
         PBIOMASS_LEAF(:)*(1.0-EXP(-1.0/ZXSEFOLD(:))))  
!
! avoid negative values due to computation precision
!
  ZDECLINE(:,1) = MAX(ZDECLINE(:,1),0.0)
!
! daily active biomass assimilation
!
  ZASSIM(:) = PANDAY(:)*ZBMCOEF
  PINCREASE(:,1) = ZASSIM(:)
!
! current leaf biomass with assimilation and senescence
!
  PBIOMASS_LEAF(:) = PBIOMASS_LEAF(:) - ZDECLINE(:,1)
!
END WHERE
!
!
!-------------------------------------------------------------------------------
!
! 3 - Evolution of active structural biomass
! ------------------------------------------
!
WHERE( OMASK(:))
  !
  WHERE (ZASSIM(:) >= ZDECLINE(:,1))
    !
    ! 3.1 - Growing phase : plant nitrogen decline theory
    !
    ! the growth allometric law is applied
    !
    ZBIOMASST(:) = MAX(PBIOMASS_LEAF(:), &
                        (PBIOMASS_LEAF(:)/ZCC_NITRO)**(1.0/(1.0-XCA_NIT)))  
    !
    ! active structural biomass increment and storage
    !
    ZBIOMASS(:,2) = ZBIOMASST(:)-PBIOMASS_LEAF(:)
    ZDECLINE(:,2) = ZBIOMASS(:,2)*(1.0-EXP(-1.0*XDAY/PSEFOLD(:)))
    ZSTORAGE(:,1) = ZBIOMASS(:,2)-PBIOMASS(:,2) +ZDECLINE(:,2)+PRESP_BIOMASS(:,2)
    PINCREASE(:,2) = ZSTORAGE(:,1)
    !
  ELSEWHERE
    !
    ! 3.2 - Senescence phase
    !
    ! the active structural biomass dies exponentially at the lowest rate
    !
    ZSTORAGE(:,1) = 0.0
    PINCREASE(:,2) = ZSTORAGE(:,1)
    ZDECLINE(:,2) = PBIOMASS(:,2)*(1.0-EXP(-1.0*XDAY/PSEFOLD(:)))
    ZBIOMASS(:,2) = PBIOMASS(:,2) - ZDECLINE(:,2)-PRESP_BIOMASS(:,2)
    !
    !  Avoid negative values of biomass
    !
    ZBIOMASS(:,2) = MAX(ZBIOMASS(:,2),0.0)
    ZBIOMASST(:) = PBIOMASS_LEAF(:)+ZBIOMASS(:,2)
    !
  END WHERE
  !
END WHERE
!
!
! 3.3 - Flow to the passive structural biomass: cut or growth after senescence
! Biomass is taken from active structural biomass, not from senescence of leaves
!
WHERE( OMASK(:) )
!
   PINCREASE(:,3) = -MIN(ZSTORAGE(:,1),0.0)
   ZSTORAGE(:,1) = MAX(0.0,ZSTORAGE(:,1))
!
END WHERE
!
! 3.4 - Mass conservation : leaf biomass sensecence must be >= structural storage
!
WHERE( ZSTORAGE(:,1) > ZDECLINE(:,1) .AND. OMASK(:) )
!
   ZDECLINE(:,2) = PBIOMASS(:,2)*(1.0 - EXP(-1.0*XDAY/PSEFOLD(:)))
   ZBIOMASST(:) = PBIOMASS(:,1) + PBIOMASS(:,2) &
                    - ZDECLINE(:,2) - PRESP_BIOMASS(:,2)  
   PBIOMASS_LEAF(:) = ZCC_NITRO*(ZBIOMASST(:)**(1.0-XCA_NIT))
   ZBIOMASS(:,2)  = ZBIOMASST(:) - PBIOMASS_LEAF(:)
   ZDECLINE(:,1) = PBIOMASS(:,1) - PBIOMASS_LEAF(:)
   ZSTORAGE(:,1) = ZBIOMASS(:,2) - PBIOMASS(:,2) &
                          + ZDECLINE(:,2) + PRESP_BIOMASS(:,2)  
   PINCREASE(:,2) = ZSTORAGE(:,1)
!
END WHERE
!
!-------------------------------------------------------------------------------
!
! 4 - Evolution of other biomass pools and final calculations
! -----------------------------------------------------------
!
! 4.1 - Mortality of leaf biomass
!
WHERE( OMASK(:) )
!
  ZMORT_LEAF(:) = MAX(0.0, ZDECLINE(:,1) - ZSTORAGE(:,1))
!
END WHERE
!
!
WHERE ( OMASK(:))
  !
  WHERE (ZASSIM(:) > ZDECLINE(:,1))
    !
    ! Growing phase, all leaf decline is used as storage.
    ! Remaining mortality is stored in roots.
    !
    ZSTORAGE(:,1) = ZSTORAGE(:,1) + ZMORT_LEAF(:)
    PINCREASE(:,4) = ZMORT_LEAF(:)
    ZMORT_LEAF(:) = 0.
    !
  ELSEWHERE
    !
    ! Senescence, a part of mortality is stored in roots, limited by assimilation
    ! rate.
    !
    ZSTORAGE(:,1) = ZSTORAGE(:,1) +  MIN(MAX(0.5*ZASSIM(:),0.),0.5*ZMORT_LEAF(:))
    PINCREASE(:,4) = MIN(MAX(0.5*ZASSIM(:),0.),0.5*ZMORT_LEAF(:))
    ZMORT_LEAF(:) = ZMORT_LEAF(:) - MIN(MAX(0.5*ZASSIM(:),0.),0.5*ZMORT_LEAF(:))
    !
  END WHERE
  !
END WHERE
!
!
! 4.2 - Evolution of the other reservoirs
!
! 4.2.1 - senesence, avoiding negative values of biomass
!
WHERE (OMASK(:)) 
  !
  ZDECLINE(:,3) = MIN(PBIOMASS(:,3)*(1.0-EXP(-1.0*XDAY/(PSEFOLD(:)/4.))), &
                      PBIOMASS(:,3)-PRESP_BIOMASS(:,3))
  !            
  ZDECLINE(:,4) = MIN(PBIOMASS(:,4)*(1.0-EXP(-1.0*XDAY/PSEFOLD(:))), &
                      PBIOMASS(:,4)-PRESP_BIOMASS(:,4))
  !
  WHERE (PVEGTYPE(:,NVT_TREE)+PVEGTYPE(:,NVT_CONI)+PVEGTYPE(:,NVT_EVER) > 0.5)
    !
    ! Woody
    !
    ZDECLINE(:,5) = PBIOMASS(:,5)*(1.0-EXP(-1.0*XDAY/PTAU_WOOD(:)))

    ZDECLINE(:,6) = PBIOMASS(:,6)*(1.0-EXP(-1.0*XDAY/PTAU_WOOD(:)))

  ELSEWHERE (PVEGTYPE(:,NVT_TREE)+PVEGTYPE(:,NVT_CONI)+PVEGTYPE(:,NVT_EVER) < 0.5)
    !
    ! Herbaceous
    !
    ZDECLINE(:,5) = 0.
    ZDECLINE(:,6) = 0.
    !
  END WHERE
  !
END WHERE
!
!
! 4.2.2 - storage (part of decline used as input for other reservoirs)
!
WHERE ( OMASK(:))
  !
  ! Woody
  !
  WHERE (PVEGTYPE(:,NVT_TREE)+PVEGTYPE(:,NVT_CONI)+PVEGTYPE(:,NVT_EVER) > 0.5)
    !
    WHERE (ZASSIM(:) >= ZDECLINE(:,1))
      !
      ! Growing phase, all decline is used as storage
      !
      ZSTORAGE(:,2) = ZDECLINE(:,2)
      ZSTORAGE(:,3) = ZDECLINE(:,3)
      ZSTORAGE(:,4) = ZDECLINE(:,4)
      ZSTORAGE(:,5) = 0.
      ZSTORAGE(:,6) = 0.
      !
      PINCREASE(:,4) = PINCREASE(:,4) + 0.3*ZSTORAGE(:,2) + 0.3*ZSTORAGE(:,3)
      PINCREASE(:,5) = 0.7*ZSTORAGE(:,2) + 0.7*ZSTORAGE(:,3)
      PINCREASE(:,6) = ZSTORAGE(:,4)
      !
    ELSEWHERE  
      !
      ! Senescence, only a part of decline is used as storage
      !
      ZSTORAGE(:,2) = 0.5*ZDECLINE(:,2)
      ZSTORAGE(:,3) = 0.5*ZDECLINE(:,3)
      ZSTORAGE(:,4) = 0.5*ZDECLINE(:,4)
      ZSTORAGE(:,5) = 0.
      ZSTORAGE(:,6) = 0.
      !
      PINCREASE(:,5) = ZSTORAGE(:,2) + ZSTORAGE(:,3)
      PINCREASE(:,6) = ZSTORAGE(:,4)
      !
    END WHERE
    !
  ELSEWHERE (PVEGTYPE(:,NVT_TREE)+PVEGTYPE(:,NVT_CONI)+PVEGTYPE(:,NVT_EVER) < 0.5)
    !
    ! Herbaceous
    !
    WHERE (ZASSIM(:) >= ZDECLINE(:,1))
      !
      ! Growing phase, all decline is used as storage
      !
      ZSTORAGE(:,2) = ZDECLINE(:,2)
      ZSTORAGE(:,3) = ZDECLINE(:,3)
      ZSTORAGE(:,4) = 0.
      ZSTORAGE(:,5) = 0.
      ZSTORAGE(:,6) = 0.
      !
      PINCREASE(:,4) = PINCREASE(:,4) + ZSTORAGE(:,2) + ZSTORAGE(:,3)
      !
    ELSEWHERE  
      !
      ! Senescence, no storage
      !
      ZSTORAGE(:,2) = 0.
      ZSTORAGE(:,3) = 0.
      ZSTORAGE(:,4) = 0.
      ZSTORAGE(:,5) = 0.
      ZSTORAGE(:,6) = 0.
      !
    END WHERE
    !
  END WHERE
  !
END WHERE
!
! 4.2.3 - mortality (senescence - storage) and turnover
!
IF (HRESPSL=='CNT') THEN
  WHERE( OMASK(:) )
!
    PTURNOVER(:,1) = ZMORT_LEAF(:)*1000.*XPCCO2/XDAY
    PTURNOVER(:,2) = (ZDECLINE(:,2) - ZSTORAGE(:,2))*1000.*XPCCO2/XDAY
    PTURNOVER(:,3) = (ZDECLINE(:,3) - ZSTORAGE(:,3))*1000.*XPCCO2/XDAY
    PTURNOVER(:,4) = (ZDECLINE(:,4) - ZSTORAGE(:,4))*1000.*XPCCO2/XDAY
    PTURNOVER(:,5) = (ZDECLINE(:,5) - ZSTORAGE(:,5))*1000.*XPCCO2/XDAY
    PTURNOVER(:,6) = (ZDECLINE(:,6) - ZSTORAGE(:,6))*1000.*XPCCO2/XDAY
!
  ENDWHERE
ENDIF
!
! 4.2.4 - evolution of reservoirs
!
WHERE (OMASK(:))
  !
  ZBIOMASS(:,3) = PBIOMASS(:,3) + PINCREASE(:,3) - ZDECLINE(:,3) - PRESP_BIOMASS(:,3)
  ZBIOMASS(:,4) = PBIOMASS(:,4) + PINCREASE(:,4) - ZDECLINE(:,4) - PRESP_BIOMASS(:,4)
  !
  WHERE (PVEGTYPE(:,NVT_TREE)+PVEGTYPE(:,NVT_CONI)+PVEGTYPE(:,NVT_EVER) > 0.5)
    !
    ! Woody
    !
    ZBIOMASS(:,5) = PBIOMASS(:,5) + PINCREASE(:,5) - ZDECLINE(:,5)
    ZBIOMASS(:,6) = PBIOMASS(:,6) + PINCREASE(:,6) - ZDECLINE(:,6)
    !
  ELSEWHERE (PVEGTYPE(:,NVT_TREE)+PVEGTYPE(:,NVT_CONI)+PVEGTYPE(:,NVT_EVER) < 0.5)  
    !
    ! Herbaceous
    !
    ZBIOMASS(:,5) = 0.
    ZBIOMASS(:,6) = 0.
    !
  END WHERE
  !
END WHERE
!
!
! 4.3 - Re-initialisations for next time step
!
WHERE( OMASK(:) )
!
! re-initialisation of biomass compartments values: X(day) <-- X(day-1)
!
! Add net accumulated CO2 assimilation 
!
  PBIOMASS_LEAF(:) = PBIOMASS_LEAF(:) + ZASSIM(:)
!
!debug seb
  PBIOMASS(:,1) = PBIOMASS_LEAF(:)
  PBIOMASS(:,2) = ZBIOMASS(:,2)
  PBIOMASS(:,3) = ZBIOMASS(:,3)
  PBIOMASS(:,4) = ZBIOMASS(:,4)
  PBIOMASS(:,5) = ZBIOMASS(:,5)
  PBIOMASS(:,6) = ZBIOMASS(:,6)
!
! re-initialisation of respiration and assimilation terms
!
  PRESP_BIOMASS(:,2) = 0.0
  PRESP_BIOMASS(:,3) = 0.0
  PRESP_BIOMASS(:,4) = 0.0
  PANFM(:) = 0.0
!
END WHERE
!
!-------------------------------------------------------------------------------
!
! Other calculations at every time step.
!
! 5 - Respiration of structural biomass pools (no respiration of woody pools)
! ---------------------------------------------------------------------------
!
PRESP_BIOMASS(:,2) = PRESP_BIOMASS(:,2) + &
        PBIOMASS(:,2)*XRESPFACTOR_NIT*2.0**((PTG(:,1)-XTT-25.)/10.0)*PTSTEP  
PRESP_BIOMASS(:,2) = MIN(PRESP_BIOMASS(:,2), PBIOMASS(:,2))
      
PRESP_BIOMASS(:,3) = PRESP_BIOMASS(:,3) + &
        PBIOMASS(:,3)*MAX( 0., XCOEFF_MAINT_RESP_ZERO * &
        (1.+XSLOPE_MAINT_RESP*(PTG(:,1)-XTT)))*PTSTEP  
PRESP_BIOMASS(:,3) = MIN(PRESP_BIOMASS(:,3), PBIOMASS(:,3))

PRESP_BIOMASS(:,4) = PRESP_BIOMASS(:,4) + &
        PBIOMASS(:,4)*MAX( 0., XCOEFF_MAINT_RESP_ZERO * &
        (1.+XSLOPE_MAINT_RESP*(PTG(:,2)-XTT)))*PTSTEP  
PRESP_BIOMASS(:,4) = MIN(PRESP_BIOMASS(:,4), PBIOMASS(:,4))
IF (LHOOK) CALL DR_HOOK('NITRO_CARBON_DECLINE',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE NITRO_CARBON_DECLINE
