!     #########
    SUBROUTINE VEGETATION_EVOL(HPHOTO, HRESPSL, HALBEDO, OAGRIP, OTR_ML,  &
                               PTSTEP, KMONTH, KDAY, PTIME, PLAT, PRHOA,  &
                               PTG, PALBNIR_VEG, PALBVIS_VEG, PALBUV_VEG, &
                               PALBNIR_SOIL, PALBVIS_SOIL, PALBUV_SOIL,   &
                               PVEGTYPE, PSEFOLD, PANMAX, PH_TREE, PBSLAI,& 
                               PLAIMIN, P_CO2, PCE_NITRO, PCF_NITRO,      &
                               PCNA_NITRO, PBSLAI_NITRO, PGMES, PTAU_WOOD,&
                               TPSEED, TPREAP, PAOSIP, PAOSIM, PAOSJP,    &
                               PAOSJM, PHO2IP, PHO2IM, PHO2JP, PHO2JM,    &
                               PZ0EFFIP, PZ0EFFIM, PZ0EFFJP, PZ0EFFJM,    &
                               PLAI, PVEG, PZ0, PALBNIR, PALBVIS, PALBUV, &
                               PEMIS, PANFM, PANDAY, PBIOMASS, PRESP_BIOMASS, &
                               PRESP_BIOMASS_INST, PINCREASE, PTURNOVER   )  
!   ###############################################################
!!****  *VEGETATION EVOL*
!!
!!    PURPOSE
!!    -------
!
!     performs the time evolution of vegetation parameters
!     at solar midnight in the case of interactive vegetation (ISBA-Ags)
!              
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
!!      
!!    AUTHOR
!!    ------
!!
!!	V. Masson          * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/03/03 
!!      P. Le Moigne 12/2004 : NIT version 
!!      P Le Moigne  09/2005 : AGS modifs of L. Jarlan
!!      A.L. Gibelin 04/2009 : BIOMASS and RESP_BIOMASS arrays
!!      A.L. Gibelin 04/2009 : Add NCB option 
!!      D. Carrer    01/2012 : epresentation of nitrogen dilution fct of CO2 (from Calvet et al. 2008)
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_CO2V_PAR,       ONLY : XMC, XMCO2, XPCCO2, XRESPFACTOR_NIT, &
                                XCOEFF_MAINT_RESP_ZERO, XSLOPE_MAINT_RESP, &
                                XCNAREF, XPARAM                   
USE MODD_CSTS,           ONLY : XDAY, XTT, XMD
!
USE MODI_ALBEDO
USE MODI_LAIGAIN
USE MODI_LAILOSS
USE MODI_NITRO_DECLINE
USE MODI_EMIS_FROM_VEG
USE MODI_VEG_FROM_LAI
USE MODI_Z0V_FROM_LAI
USE MODI_SUBSCALE_Z0EFF
USE MODD_TYPE_DATE_SURF
USE MODD_SURF_PAR
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
CHARACTER(LEN=3),     INTENT(IN)    :: HPHOTO  ! type of photosynthesis
!                                              ! 'NON'
!                                              ! 'AGS'
!                                              ! 'LAI'
CHARACTER(LEN=3),     INTENT(IN)    :: HRESPSL ! Soil Respiration
!                                              ! 'DEF' = Norman 1992
!                                              ! 'PRM' = Rivalland PhD Thesis (2003)
!                                              ! 'CNT' = CENTURY model (Gibelin 2008)
CHARACTER(LEN=4),     INTENT(IN)    :: HALBEDO ! albedo type
!                                              ! 'DRY ' 
!                                              ! 'EVOL' 
!                                              ! 'WET ' 
!                                              ! 'USER'
LOGICAL,              INTENT(IN)    :: OAGRIP  ! agricultural practices
LOGICAL,              INTENT(IN)    :: OTR_ML  ! new radiative transfert
!
REAL,                 INTENT(IN)    :: PTSTEP  ! time step
INTEGER,              INTENT(IN)    :: KMONTH  ! current month
INTEGER,              INTENT(IN)    :: KDAY    ! current day
REAL,                 INTENT(IN)    :: PTIME   ! current time since midnight
REAL,   DIMENSION(:), INTENT(IN)    :: PLAT    ! latitude of each grid point
REAL,   DIMENSION(:), INTENT(IN)    :: PRHOA   ! air density
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PTG     ! soil layer average temperatures (K)
REAL,   DIMENSION(:), INTENT(IN)    :: PALBVIS_VEG ! visible, near infra-red and UV
REAL,   DIMENSION(:), INTENT(IN)    :: PALBNIR_VEG ! albedo of the vegetation
REAL,   DIMENSION(:), INTENT(IN)    :: PALBUV_VEG  !
REAL,   DIMENSION(:), INTENT(IN)    :: PALBVIS_SOIL! visible, near infra-red and UV
REAL,   DIMENSION(:), INTENT(IN)    :: PALBNIR_SOIL! soil albedo
REAL,   DIMENSION(:), INTENT(IN)    :: PALBUV_SOIL !
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PVEGTYPE! fraction of each vegetation type
REAL,   DIMENSION(:), INTENT(IN)    :: PSEFOLD ! e-folding time for senescence (s)
REAL,   DIMENSION(:), INTENT(IN)    :: PANMAX  ! maximum photosynthesis rate
REAL,   DIMENSION(:), INTENT(IN)    :: PH_TREE ! height of trees
REAL,   DIMENSION(:), INTENT(IN)    :: PBSLAI  ! ratio of biomass to LAI
REAL,   DIMENSION(:), INTENT(IN)    :: PLAIMIN ! minimum LAI
!
REAL,   DIMENSION(:), INTENT(IN)    :: P_CO2 ! CO2 concentration [ppmm]
!
REAL,   DIMENSION(:), INTENT(IN)    :: PCE_NITRO    ! leaf aera ratio sensibility to nitrogen 
!                                                     concentration (10**2 m2 kg-1)
REAL,   DIMENSION(:), INTENT(IN)    :: PCF_NITRO    ! lethal minimum value of leaf aera ratio 
!                                                     (m2 kg-1)
REAL,   DIMENSION(:), INTENT(IN)    :: PCNA_NITRO   ! nitrogen concentration of active biomass (%)
REAL,   DIMENSION(:), INTENT(IN)    :: PBSLAI_NITRO ! ratio of biomass to LAI
!
REAL,   DIMENSION(:), INTENT(IN)    :: PGMES      ! mesophyll conductance (m s-1)
REAL,   DIMENSION(:), INTENT(IN)    :: PTAU_WOOD  ! residence time in wood (s)
!
!
TYPE (DATE_TIME),   DIMENSION(:), INTENT(IN) :: TPSEED ! seeding date
TYPE (DATE_TIME),   DIMENSION(:), INTENT(IN) :: TPREAP ! reaping date
!
REAL, DIMENSION(:), INTENT(IN)  :: PAOSIP  ! A/S for increasing x
REAL, DIMENSION(:), INTENT(IN)  :: PAOSIM  ! A/S for decreasing x
REAL, DIMENSION(:), INTENT(IN)  :: PAOSJP  ! A/S for increasing y
REAL, DIMENSION(:), INTENT(IN)  :: PAOSJM  ! A/S for decreasing y
REAL, DIMENSION(:), INTENT(IN)  :: PHO2IP  ! h/2 for increasing x
REAL, DIMENSION(:), INTENT(IN)  :: PHO2IM  ! h/2 for decreasing x
REAL, DIMENSION(:), INTENT(IN)  :: PHO2JP  ! h/2 for increasing y
REAL, DIMENSION(:), INTENT(IN)  :: PHO2JM  ! h/2 for decreasing y
!
REAL, DIMENSION(:), INTENT(INOUT) :: PZ0EFFIP! roughness length for increasing x
REAL, DIMENSION(:), INTENT(INOUT) :: PZ0EFFIM! roughness length for decreasing x
REAL, DIMENSION(:), INTENT(INOUT) :: PZ0EFFJP! roughness length for increasing y
REAL, DIMENSION(:), INTENT(INOUT) :: PZ0EFFJM! roughness length for decreasing y
!
REAL,   DIMENSION(:), INTENT(INOUT) :: PLAI    ! leaf area index (LAI) 
REAL,   DIMENSION(:), INTENT(INOUT) :: PVEG    ! vegetation fraction
REAL,   DIMENSION(:), INTENT(INOUT) :: PZ0     ! roughness length: momentum
REAL,   DIMENSION(:), INTENT(INOUT) :: PALBNIR ! snow-free near-infra-red albedo
REAL,   DIMENSION(:), INTENT(INOUT) :: PALBVIS ! snow-free visible albedo
REAL,   DIMENSION(:), INTENT(INOUT) :: PALBUV  ! snow-free UV albedo
REAL,   DIMENSION(:), INTENT(INOUT) :: PEMIS   ! snow-free emissivity
!
REAL,   DIMENSION(:), INTENT(INOUT) :: PANFM              ! maximum leaf assimilation
REAL,   DIMENSION(:), INTENT(INOUT) :: PANDAY             ! daily net CO2 assimilation
REAL, DIMENSION(:,:), INTENT(INOUT) :: PBIOMASS           ! biomass of day-1
REAL, DIMENSION(:,:), INTENT(INOUT) :: PRESP_BIOMASS      ! daily cumulated respiration of biomass (kg/m2/day)
REAL, DIMENSION(:,:), INTENT(INOUT) :: PRESP_BIOMASS_INST ! instantaneous respiration of biomass (kgCO2/kgair m/s)
!
REAL, DIMENSION(:,:), INTENT(OUT)   :: PINCREASE          ! increment of biomass
REAL, DIMENSION(:,:), INTENT(OUT)   :: PTURNOVER          ! biomass turnover going into litter
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PRESP_BIOMASS,1),SIZE(PRESP_BIOMASS,2)) :: ZRESP_BIOMASS_LAST ! biomass at t-1 (kg/m2/s)
REAL,    DIMENSION(SIZE(PLAI))    :: ZBIOMASS_LEAF   ! temporary leaf biomass 
REAL,    DIMENSION(SIZE(PLAI))    :: ZBSLAI_NITRO    ! (Calvet et al. 2008) ratio of biomass to LAI
                                                     ! with representation of nitrogen dilution
REAL,    DIMENSION(SIZE(PLAI)) :: ZCO2, ZCNA_NITRO   ! fct of CO2        
REAL,    DIMENSION(SIZE(PLAI)) ::  ZCNAREF, ZPARAM
INTEGER, DIMENSION(1)          :: IDMAX
!
LOGICAL, DIMENSION(SIZE(PLAI))    :: GMASK_AGRI
LOGICAL                           :: GMASK
INTEGER                           :: JWRK, JJ
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-----------------------------------------------------------------
!
!*      1.     Preliminaries
!              -------------
!
IF (LHOOK) CALL DR_HOOK('VEGETATION_EVOL',0,ZHOOK_HANDLE)
!
! Mask where vegetation evolution is performed (just before solar midnight)
GMASK = ( PTIME - PTSTEP < 0. ) .AND. ( PTIME >= 0. )
!
! Save RESP_BIOMASS at t-1
IF (GMASK) THEN
  PRESP_BIOMASS     (:,1) = 0.0
  ZRESP_BIOMASS_LAST(:,:) = 0.0
ELSE
  PRESP_BIOMASS     (:,1) = PRESP_BIOMASS(:,1) + PRESP_BIOMASS_INST(:,1) * (PTSTEP*PRHOA(:)*XMC)/(XPCCO2*XMCO2)
  ZRESP_BIOMASS_LAST(:,:) = PRESP_BIOMASS(:,:)
ENDIF
!
!*      2.     Interactive vegetation
!              ----------------------
!
!  LAI daily mortality and assimilation
!
ZBIOMASS_LEAF(:) = PBIOMASS(:,1)
!
IF (GMASK) THEN
  IF (HPHOTO=='LAI' .OR. HPHOTO=='LST') THEN
    CALL LAILOSS(PVEG, PSEFOLD, PANMAX, PANDAY, PANFM, ZBIOMASS_LEAF)  
    CALL LAIGAIN(PBSLAI, PLAIMIN, PVEG, ZBIOMASS_LEAF, PLAI, PANDAY)
    PBIOMASS(:,1) = ZBIOMASS_LEAF(:)
  ELSE IF (HPHOTO=='NIT' .OR. HPHOTO=='NCB') THEN
    ! ratio of biomass to LAI with representation of nitrogen dilution fct of CO2
    DO JJ = 1, SIZE(PLAI)
      IDMAX = MAXLOC(PVEGTYPE(JJ,:))
      ZCNAREF(JJ) = XCNAREF(IDMAX(1))
      ZPARAM (JJ) = XPARAM (IDMAX(1))
      !--------representation of nitrogen dilution fct of CO2 (Calvet et al. 2008)-------------
      IF ( (PCE_NITRO (JJ)*PCNA_NITRO(JJ)+PCF_NITRO (JJ))/=0. .AND. ZCNAREF(JJ).NE.0. ) THEN 
        ZCO2        (JJ) = P_CO2(JJ)/1.e-6/XMCO2*XMD  ! (ppmm ->  ppm)
        ZCNA_NITRO  (JJ) = ZCNAREF(JJ)*EXP(-0.048*EXP(ZPARAM(JJ)-ZCNAREF(JJ)/6.3)*ALOG(ZCO2(JJ)/371.))
        ZBSLAI_NITRO(JJ) = 1. / (PCE_NITRO(JJ)*PCNA_NITRO(JJ)+PCF_NITRO(JJ))
      ELSE
        ZBSLAI_NITRO(JJ) = PBSLAI_NITRO(JJ)
      ENDIF
    ENDDO      
    CALL NITRO_DECLINE(HPHOTO, HRESPSL, OTR_ML,                             &
                       ZBSLAI_NITRO, PSEFOLD, PGMES, PANMAX, PANDAY,        &
                       PLAT, PLAIMIN, PVEGTYPE, PTAU_WOOD,                  &
                       PANFM, PLAI, PBIOMASS, PRESP_BIOMASS, ZBIOMASS_LEAF, &
                       PINCREASE, PTURNOVER                               )
    CALL LAIGAIN(ZBSLAI_NITRO, PLAIMIN, PVEG, ZBIOMASS_LEAF, PLAI, PANDAY)
  ENDIF
  ! CASE CPHOTO=AST reinitialise  PANDAY and PANFM 
  PANDAY=0.0
  PANFM =0.0  
ENDIF
!
!
IF (HPHOTO == 'NIT' .OR. HPHOTO=='NCB') THEN
  !
  ! 5 - Respiration of structural biomass pools
  !
  PRESP_BIOMASS(:,2) = PRESP_BIOMASS(:,2)  + PBIOMASS(:,2) * XRESPFACTOR_NIT &
       * 2.0**((PTG(:,1)-XTT-25.)/10.0) * PTSTEP  

  IF (HPHOTO == 'NIT') THEN
    !
    PRESP_BIOMASS(:,3) = PRESP_BIOMASS(:,3) + PBIOMASS(:,3) * XRESPFACTOR_NIT &
       * 2.0**((PTG(:,2)-XTT-25.)/10.0) * PTSTEP  
    !
  ELSEIF (HPHOTO == 'NCB') THEN
    !
    PRESP_BIOMASS(:,2) = MIN(PRESP_BIOMASS(:,2), PBIOMASS(:,2))
    ! 
    PRESP_BIOMASS(:,3) = PRESP_BIOMASS(:,3) + PBIOMASS(:,3) * MAX( 0., &
        XCOEFF_MAINT_RESP_ZERO * (1. + XSLOPE_MAINT_RESP*(PTG(:,1)-XTT))) * PTSTEP  
    PRESP_BIOMASS(:,3) = MIN(PRESP_BIOMASS(:,3), PBIOMASS(:,3))
    ! 
    PRESP_BIOMASS(:,4) = PRESP_BIOMASS(:,4) + PBIOMASS(:,4) * MAX( 0., &
        XCOEFF_MAINT_RESP_ZERO * (1. + XSLOPE_MAINT_RESP*(PTG(:,2)-XTT))) * PTSTEP  
    PRESP_BIOMASS(:,4) = MIN(PRESP_BIOMASS(:,4), PBIOMASS(:,4))
    !
  ENDIF
  !
  ! Instantaneous respiration (kgCO2/kgair m/s)
  !
  DO JWRK=2,SIZE(PRESP_BIOMASS,2)
      PRESP_BIOMASS_INST(:,JWRK) = (PRESP_BIOMASS(:,JWRK) - ZRESP_BIOMASS_LAST(:,JWRK)) &
                                     * XPCCO2*XMCO2/(PTSTEP*PRHOA(:)*XMC)  
  ENDDO
ENDIF

!*      3.     Agricultural practices
!              ----------------------
!
IF (OAGRIP) THEN
  !
  GMASK_AGRI(:) = .FALSE.
  WHERE ( TPSEED(:)%TDATE%MONTH /= NUNDEF .AND. ( KMONTH < TPSEED(:)%TDATE%MONTH .OR. &
         (KMONTH == TPSEED(:)%TDATE%MONTH .AND. KDAY < TPSEED(:)%TDATE%DAY) ) )  GMASK_AGRI(:) = .TRUE.
  WHERE ( TPREAP(:)%TDATE%MONTH /= NUNDEF .AND. ( KMONTH > TPREAP(:)%TDATE%MONTH .OR. &
         (KMONTH == TPREAP(:)%TDATE%MONTH .AND. KDAY >= TPREAP(:)%TDATE%DAY) ) ) GMASK_AGRI(:) = .TRUE. 
  !
  WHERE (GMASK_AGRI(:))
    PLAI(:)             = PLAIMIN(:)
    ZBIOMASS_LEAF(:)    = PLAI(:) * ZBSLAI_NITRO(:)
  END WHERE

  IF (HPHOTO == 'NIT' .OR. HPHOTO == 'NCB') THEN
    !
    WHERE (GMASK_AGRI(:))
      PBIOMASS(:,1)       = 0.0
      PBIOMASS(:,2)       = 0.0
      PBIOMASS(:,3)       = 0.0
      PRESP_BIOMASS(:,2)  = 0.0
      PRESP_BIOMASS(:,3)  = 0.0
    END WHERE
    !
    IF (HPHOTO == 'NCB') THEN
      !
      WHERE (GMASK_AGRI(:)) 
        PBIOMASS(:,4)       = 0.0
        PBIOMASS(:,5)       = 0.0
        PBIOMASS(:,6)       = 0.0
        PRESP_BIOMASS(:,4)  = 0.0
      END WHERE
      !
    ENDIF
    !
  ENDIF
  !
ENDIF
!
!*      4.     Physical parameters depending on vegetation
!              -------------------------------------------
!
IF (GMASK) THEN
  !
  WHERE( PVEG(:) > 0. )
    ! Evolution of vegetation fraction and roughness length due to LAI change
    PZ0 (:) = Z0V_FROM_LAI(PLAI(:),PH_TREE(:),PVEGTYPE(:,:)) 
    PVEG(:) = VEG_FROM_LAI(PLAI(:),PVEGTYPE(:,:))
    !
    ! Evolution of radiative parameters due to vegetation fraction change
    PEMIS(:)= EMIS_FROM_VEG(PVEG(:),PVEGTYPE(:,:))
  END WHERE
  !
  CALL ALBEDO(HALBEDO,                                  &
              PALBVIS_VEG,PALBNIR_VEG,PALBUV_VEG,PVEG,  &
              PALBVIS_SOIL,PALBNIR_SOIL,PALBUV_SOIL,    &
              PALBVIS,PALBNIR,PALBUV                    )  
  !
  ! Evolution of effective roughness length due to new surface roughness length
  !
  CALL SUBSCALE_Z0EFF(PAOSIP,PAOSIM,PAOSJP,PAOSJM,         &
                      PHO2IP,PHO2IM,PHO2JP,PHO2JM,PZ0,     &
                      PZ0EFFIP,PZ0EFFIM,PZ0EFFJP,PZ0EFFJM  ) 
  !
ENDIF
!
IF (LHOOK) CALL DR_HOOK('VEGETATION_EVOL',1,ZHOOK_HANDLE)
!-----------------------------------------------------------------
!
END SUBROUTINE VEGETATION_EVOL
