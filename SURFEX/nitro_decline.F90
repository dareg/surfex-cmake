!     #########
SUBROUTINE NITRO_DECLINE(OMASK,PTSTEP,PBSLAI_NITRO,PSEFOLD,PGMES,PANMAX,PANDAY,       &
       PLAT,PLAIMIN,PTG,PANFM,PLAI,PBIOMASS,PRESP_BIOMASS,PBIOMASS_LEAF)  
!
!   ###############################################################
!!**  NITRO_DECLINE 
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!     Calvet and Soussana (2001) and Gibelin et al. (2006)
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
!!      
!!    AUTHOR
!!    ------
!!
!!	A. Boone           * Meteo-France *
!!      V. Rivalland
!!      (following Belair)
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    27/01/03 
!!
!!      P Le Moigne  09/2005 : AGS modifs of L. Jarlan
!!      A.L. Gibelin 04/2009 : BIOMASS and RESP_BIOMASS arrays
!!      A.L. Gibelin 04/2009 : Suppress unused arguments
!!      A.L. Gibelin 04/2009 : Suppress unused modules and add ONLY
!!      A. Barbu     01/2011 : modification of leaf reservoir
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CSTS,           ONLY : XPI, XDAY, XTT
USE MODD_CO2V_PAR,       ONLY : XCC_NIT, XCA_NIT, XRESPFACTOR_NIT, XMC, XMCO2, XPCCO2
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
!
REAL,   DIMENSION(:), INTENT(IN) :: PBSLAI_NITRO     ! ratio of biomass to LAI
REAL,   DIMENSION(:), INTENT(IN) :: PSEFOLD          ! e-folding time for senescence (s)
REAL,   DIMENSION(:), INTENT(IN) :: PGMES            ! mesophyll conductance (m s-1)
REAL,   DIMENSION(:), INTENT(IN) :: PANMAX           ! maximum photosynthesis rate
REAL,   DIMENSION(:), INTENT(IN) :: PLAT             !latitude of each grid point
REAL,   DIMENSION(:), INTENT(IN) :: PLAIMIN          ! minimum LAI
REAL, DIMENSION(:,:), INTENT(IN) :: PTG              ! soil layer average temperatures (K)
REAL,   DIMENSION(:), INTENT(IN) :: PANDAY           ! daily net CO2 accumulation
!
REAL,   DIMENSION(:), INTENT(INOUT) :: PANFM         ! maximum leaf assimilation
REAL,   DIMENSION(:), INTENT(INOUT) :: PLAI          ! leaf area index (LAI) 
REAL, DIMENSION(:,:), INTENT(INOUT) :: PBIOMASS      ! biomass reservoirs
REAL, DIMENSION(:,:), INTENT(INOUT) :: PRESP_BIOMASS ! cumulated daily biomass respiration
!
REAL,   DIMENSION(:), INTENT(OUT)   :: PBIOMASS_LEAF ! temporary leaf biomass
!
!
!*      0.2    declarations of local variables
!
REAL                            :: ZCC_NITRO         ! c coefficient for nitrogen dilution
REAL                            :: ZBIOMASST_LIM     ! threshold value of ZBIOMASST 
!                                                      in nitrogen dilution theory
!
REAL                            :: ZBMCOEF
!
REAL,    DIMENSION(SIZE(PLAI))  :: ZXSEFOLD          ! e-folding time for senescence corrected (days)
REAL,    DIMENSION(SIZE(PLAI))  :: ZLAIB_NITRO       ! LAI correction parameter used in sefold calculation
REAL,    DIMENSION(SIZE(PLAI))  :: ZASSIM            ! assimilation
REAL,    DIMENSION(SIZE(PLAI))  :: ZBIOMASST         ! leaf+structural dry above ground biomass
!
REAL,    DIMENSION(SIZE(PLAI),SIZE(PBIOMASS,2)) :: ZBIOMASS      ! temporary biomass reservoirs
REAL,    DIMENSION(SIZE(PLAI),SIZE(PBIOMASS,2)) :: ZDECLINE      ! biomass decline (storage+mortality)
REAL,    DIMENSION(SIZE(PLAI),SIZE(PBIOMASS,2)) :: ZSTORAGE      ! storage
REAL,    DIMENSION(SIZE(PLAI))                  :: ZMORT_LEAF    ! leaf mortality
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
! correspondance between indices of arrays and biomass compartments
! LEAF = 1
! STRUCT_ACT = 2
! STRUCT_PAS = 3
!
!-------------------------------------------------------------------------------
!
! 1 - Initialisations
! -------------------
!
IF (LHOOK) CALL DR_HOOK('NITRO_DECLINE',0,ZHOOK_HANDLE)
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
! current biomass value:
!
PBIOMASS_LEAF(:)         = PBIOMASS(:,1)
!
!-------------------------------------------------------------------------------
!
! Once a day (at midnight),repartition of net assimilation and mortality 
! into different biomass compartments.
!
! 2 - Evolution of leaf biomass
! -----------------------------

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
! LAI correction for shadow effect
!
  ZLAIB_NITRO(:) = MAX( 5.76-0.64*ATAN(ABS(PLAT(:))*XPI/180.),3.8 )
!
! leaf life expectancy
! PSEFOLD  is in second
! ZXSEFOLD is in Day
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
! Following Marita's work limitation of the senesence
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
!
! current leaf biomass with assimilation and senescence
!
  PBIOMASS_LEAF(:) = PBIOMASS_LEAF(:) - ZDECLINE(:,1)
!
! senesence of deep-structural biomass
!
  ZBIOMASS(:,3) = PBIOMASS(:,3)
  ZDECLINE(:,3) = ZBIOMASS(:,3)*(1.0-EXP(-1.0*XDAY/PSEFOLD(:)))
!
END WHERE
!
!
!-------------------------------------------------------------------------------
!
! 3 - Evolution of structural biomass
! -----------------------------------
!
WHERE( OMASK(:)) 
  ! 
  WHERE (ZASSIM(:) >= ZDECLINE(:,1))
    !
    ! 3.1 - Growing phase : plant nitrogen decline theory
    !
    ! the growth allometric law is applied
    ! repartition of total biomass
    !
    ZBIOMASST(:) = MAX(PBIOMASS_LEAF(:), &
                          (PBIOMASS_LEAF(:)/ZCC_NITRO)**(1.0/(1.0-XCA_NIT)))  
    !
    ! above-ground structure biomass increment and storage
    !
    ZBIOMASS(:,2) = ZBIOMASST(:)-PBIOMASS_LEAF(:)
    ZDECLINE(:,2) = ZBIOMASS(:,2)*(1.0-EXP(-1.0*XDAY/PSEFOLD(:)))
    ZSTORAGE(:,1) = ZBIOMASS(:,2)-PBIOMASS(:,2)+ZDECLINE(:,2)+PRESP_BIOMASS(:,2)
    !
  ELSE WHERE
    !
    ! 3.2 - Senescence phase
    !
    ! the structural biomass dies exponetially at the lowest rate
    !
    ZSTORAGE(:,1) = 0.0
    ZDECLINE(:,2) = PBIOMASS(:,2)*(1.0-EXP(-1.0*XDAY/PSEFOLD(:)))
    ZBIOMASS(:,2) = PBIOMASS(:,2)-ZDECLINE(:,2)-PRESP_BIOMASS(:,2)
    !
    !  Avoid negative values of biomass
    !  No test on ZDECLINE(:,2) as it is not used after, or recalculated
    !  No test on PRESP_BIOMASS(:,2) as it should be smaller than PBIOMASS(:,2)
    !  otherwise there are irrealistic values of temperature 
    !
    ZBIOMASS(:,2) = MAX(ZBIOMASS(:,2),0.0)
    !
    ZBIOMASST(:) = PBIOMASS_LEAF(:)+ZBIOMASS(:,2)
    !
  END WHERE
  !
END WHERE
!
!
! 3.3 - Flow to the deep structural biomass
!
WHERE( OMASK(:) )
!
   ZBIOMASS(:,3) = ZBIOMASS(:,3)-MIN(ZSTORAGE(:,1),0.0)
   ZSTORAGE(:,1) = MAX(0.0,ZSTORAGE(:,1))
!
END WHERE
!
!
! 3.4 - Mass conservation : leaf biomass mortality <= structural storage
!
WHERE( ZSTORAGE(:,1) > ZDECLINE(:,1) .AND. OMASK(:) )
!
   ZDECLINE(:,2) = PBIOMASS(:,2)*(1.0 - EXP(-1.0*XDAY/PSEFOLD(:)))
   ZBIOMASST(:) = PBIOMASS(:,1) + PBIOMASS(:,2) - ZDECLINE(:,2) - PRESP_BIOMASS(:,2)
   PBIOMASS_LEAF(:) = ZCC_NITRO*(ZBIOMASST(:)**(1.0-XCA_NIT))
   ZBIOMASS(:,2) = ZBIOMASST(:) - PBIOMASS_LEAF(:)
   ZDECLINE(:,1) = PBIOMASS(:,1) - PBIOMASS_LEAF(:)
   ZSTORAGE(:,1) = ZBIOMASS(:,2) - PBIOMASS(:,2) + ZDECLINE(:,2) + PRESP_BIOMASS(:,2)
!
END WHERE
!
!-------------------------------------------------------------------------------
!
! 4 - Evolution of other biomass pools and final calculations
! -----------------------------------------------------------
!
!
! mortality
!
WHERE( OMASK(:) )
!
  ZMORT_LEAF(:) = MAX(0.0, ZDECLINE(:,1) - ZSTORAGE(:,1))
!
END WHERE
!
! emergency deep structural biomass
!
WHERE( (ZBIOMASST(:) <= ZBIOMASST_LIM) .AND. (ZXSEFOLD(:) > 1.0) &
       .AND. OMASK(:) )  
!
  ZBIOMASS(:,3) = ZBIOMASS(:,3) + ZMORT_LEAF(:)
!
END WHERE
!
WHERE( OMASK(:) )
!
! Deep structural decline
!
  ZBIOMASS(:,3) = ZBIOMASS(:,3) - ZDECLINE(:,3) - PRESP_BIOMASS(:,3)
! 
! Avoid negative values of biomass
! No test on ZDECLINE(:,3) as it is not used after
! No test on PRESP_BIOMASS(:,3) as it should be smaller than PBIOMASS(:,3)
! otherwise there are irrealistic values of temperature 
!
  ZBIOMASS(:,3) = MAX(ZBIOMASS(:,3),0.0)
!
! Add net accumulated CO2 assimilation 
!
  PBIOMASS_LEAF(:) = PBIOMASS_LEAF(:) + ZASSIM(:)
!
! re-initialisation of biomass compartments values: X(day) <-- X(day-1)
!
  PBIOMASS(:,1) = PBIOMASS_LEAF(:)
  PBIOMASS(:,2) = ZBIOMASS(:,2)
  PBIOMASS(:,3) = ZBIOMASS(:,3)
!
! re-initialisation of respiration ,mortality and assimilation terms
!
  PRESP_BIOMASS(:,2) = 0.0
  PRESP_BIOMASS(:,3) = 0.0
  PANFM(:) = 0.0
!
END WHERE
!
!-------------------------------------------------------------------------------
!
! Other calculations at every time step.
!
! 5 - Respiration of structural biomass pools
! -------------------------------------------
!
PRESP_BIOMASS(:,2) = PRESP_BIOMASS(:,2)  + PBIOMASS(:,2) *XRESPFACTOR_NIT &
       *2.0**((PTG(:,1)-XTT-25.)/10.0)*PTSTEP  

PRESP_BIOMASS(:,3) = PRESP_BIOMASS(:,3) + PBIOMASS(:,3)*XRESPFACTOR_NIT &
       *2.0**((PTG(:,2)-XTT-25.)/10.0)*PTSTEP  
IF (LHOOK) CALL DR_HOOK('NITRO_DECLINE',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE NITRO_DECLINE
