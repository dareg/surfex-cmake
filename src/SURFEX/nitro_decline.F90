!     #########
SUBROUTINE NITRO_DECLINE(IO, IP, IMT, IR, OWOOD, PBSLAI_NITRO, PLAT, PBIOMASS_LEAF )  
!
!   ###############################################################
!!**  NITRO_DECLINE 
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
!!      A.-L. Gibelin           * Meteo-France *
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
!!      A.L. Gibelin 04/2009 : adaptation to SURFEX environment
!!      A.   Barbu   01/2011 : modification of active biomass,leaf reservoir (see nitro_decline.f90)
!!      C.   Delire  04/2012 : spinup wood carbon
!!      R.   Alkama  04/2012 : 19 vegtype rather than 12
!!      B.   Decharme 05/2012: Optimization
!!                              ZCC_NITRO and ZBIOMASST_LIM in modd_co2v_par.F90
!!      C.   Delire   01/2014 : sapwood respiration from IBIS

!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_TIME_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
!
USE MODD_CSTS,           ONLY : XPI, XDAY
USE MODD_CO2V_PAR,       ONLY : XPCCO2, XCC_NIT, XCA_NIT, XMC, &
                                XMCO2, XCC_NITRO, XBIOMASST_LIM 
USE MODD_DATA_COVER_PAR, ONLY : NVT_TEBD, NVT_BONE, NVT_TRBE, NVT_TRBD, &
                                NVT_TEBE, NVT_TENE, NVT_BOBD, NVT_BOND, &
                                NVT_SHRB
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
TYPE(ISBA_PARAM_TIME_t), INTENT(INOUT) :: IMT
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
!
!*      0.1    declarations of arguments
!
LOGICAL, DIMENSION(:), INTENT(IN) :: OWOOD
!
REAL,   DIMENSION(:), INTENT(IN) :: PBSLAI_NITRO     ! ratio of biomass to LAI
REAL,   DIMENSION(:), INTENT(IN) :: PLAT             ! latitude of each grid point
!
REAL,   DIMENSION(:), INTENT(OUT)   :: PBIOMASS_LEAF ! temporary leaf biomass
!
!*      0.2    declarations of local variables
!
REAL                            :: ZBMCOEF
REAL,    DIMENSION(SIZE(IMT%XLAI,1))  :: ZXSEFOLD        ! e-folding time for senescence corrected (days)
REAL,    DIMENSION(SIZE(IMT%XLAI,1))  :: ZLAIB_NITRO     ! LAI correction parameter used in sefold calculation
REAL,    DIMENSION(SIZE(IMT%XLAI,1))  :: ZASSIM          ! assimilation
REAL,    DIMENSION(SIZE(IMT%XLAI,1))  :: ZBIOMASST       ! leaf + active structural biomass
!
REAL, DIMENSION(SIZE(IMT%XLAI,1),SIZE(IR%XBIOMASS,2))  :: ZINCREASE
REAL, DIMENSION(SIZE(IMT%XLAI,1),SIZE(IR%XBIOMASS,2))  :: ZBIOMASS      ! temporary biomass reservoirs
REAL, DIMENSION(SIZE(IMT%XLAI,1),SIZE(IR%XBIOMASS,2))  :: ZDECLINE      ! biomass decline (storage+mortality) (kgDM m-2 day-1)
REAL, DIMENSION(SIZE(IMT%XLAI,1),SIZE(IR%XBIOMASS,2))  :: ZSTORAGE      ! storage (part of decline kgDM m-2 day-1)
REAL, DIMENSION(SIZE(IMT%XLAI,1))                   :: ZMORT_LEAF    ! leaf mortality
!
REAL, DIMENSION(SIZE(IMT%XLAI,1))                   :: ZWORK,ZRESP
LOGICAL, DIMENSION(SIZE(IMT%XLAI,1))                :: GMASK_ASSIM
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
INTEGER :: JSPIN, JI, INI
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
IF (LHOOK) CALL DR_HOOK('NITRO_DECLINE',0,ZHOOK_HANDLE)
!
INI = SIZE(IMT%XLAI,1)
!
ZXSEFOLD(:)         = 0.0
ZLAIB_NITRO(:)      = 0.0
ZBIOMASST(:)        = 0.0
ZASSIM(:)           = 0.0
ZBIOMASS(:,:)       = 0.0
ZDECLINE(:,:)       = 0.0
ZINCREASE(:,:)      = 0.0
ZSTORAGE(:,:)       = 0.0
ZMORT_LEAF(:)       = 0.0
!---------------------------------------------------
!
ZBMCOEF     = XMC/(XMCO2*XPCCO2)
!
!-----------------------------------------------------------------
!avoid possible but unlikely negative values for biomass:        
!
IR%XBIOMASS(:,1,1) = MAX(IR%XBIOMASS(:,1,1),0.0)
!
! current leaf biomass value:
!
PBIOMASS_LEAF(:) = IR%XBIOMASS(:,1,1)
!print*,'leaf1 ',PBIOMASS_LEAF(1)
!
!-------------------------------------------------------------------------------
!
! Once a day (at midnight),repartition of net assimilation and mortality 
! into different biomass compartments.
!
! 2 - Evolution of leaf biomass and senescence calculations
! ---------------------------------------------------------
!
! coef c for biomass in kg/m2 now in modd_co2v_par.F90 (XCC_NITRO)
!
! LAI correction for shadow effect
IF (IO%LTR_ML) THEN
  ZLAIB_NITRO(:) = 5.30
ELSE
  ZLAIB_NITRO(:) = MAX( 5.76-0.64*ATAN(ABS(PLAT(:))*XPI/180.),3.8 )
ENDIF
!
!
! leaf life expectancy
!
!print*,'LAT ',PLAT(1)
!print*,'GMES ',IMT%XGMES(1,1)
!print*,'LAI ',IMT%XLAI(1,1)
ZWORK(:) = 0.0
WHERE(IMT%XGMES(:,1)>0.0)
      ZWORK(:) = 0.321*LOG(IMT%XGMES(:,1)*1000.)
      ZWORK(:) = EXP(ZWORK(:))*IMT%XLAI(:,1)/ZLAIB_NITRO(:)
ENDWHERE
! before optimization
!ZXSEFOLD(:)= IMT%XSEFOLD(:,1) * MAX(((IMT%XGMES(:,1)*1000.)**0.321)*IMT%XLAI(:,1)/ZLAIB_NITRO(:), 1.) * ...
!print*,'ZWORK ',ZWORK(1)
!print*,'XSEFOLD ',IMT%XSEFOLD(1,1)
!print*,'ANFM ',IR%XANFM(1,1)
!print*,'ANMAX ',IP%XANMAX(1,1)
ZXSEFOLD(:) = IMT%XSEFOLD(:,1) * MAX(1.0,ZWORK(:)) * MIN(1.0,IR%XANFM(:,1)/IP%XANMAX(:,1)) / XDAY
!
! avoid possible but unlikely division by zero
!
ZXSEFOLD(:) = MAX(1.0E-8,ZXSEFOLD(:))
!
! limitation of leaf life expectancy
!
! OLD   ZXSEFOLD(:) = MAX(5.,ZXSEFOLD(:))
! Following Marita's work limitation of the senesence
ZXSEFOLD(:) = MAX(IMT%XSEFOLD(:,1)/XDAY/10.0,ZXSEFOLD(:))
!print*,'SEFOLD ',ZXSEFOLD(1)
!
! senesence of active biomass
!
ZDECLINE(:,1) = MIN(PBIOMASS_LEAF(:)-IMT%XLAIMIN(:,1)*PBSLAI_NITRO(:), &
                    PBIOMASS_LEAF(:)*(1.0-EXP(-1.0/ZXSEFOLD(:))))
!
! avoid negative values due to computation precision
!
ZDECLINE(:,1) = MAX(ZDECLINE(:,1),0.0)
!
! current leaf biomass with assimilation and senescence
!
PBIOMASS_LEAF(:) = PBIOMASS_LEAF(:) - ZDECLINE(:,1)
!print*,'leaf2 ',PBIOMASS_LEAF(1)
!
! daily active biomass assimilation
!
ZASSIM(:) = IR%XANDAY(:,1)*ZBMCOEF
!
!-------------------------------------------------------------------------------
!
! 3 - Evolution of active structural biomass
! ------------------------------------------
!
ZWORK(:) = 0.0
WHERE(PBIOMASS_LEAF(:)>0.0)
      ZWORK(:) = (1.0/(1.0-XCA_NIT))*LOG(PBIOMASS_LEAF(:)/XCC_NITRO)
      ZWORK(:) = EXP(ZWORK(:))
ENDWHERE
!
WHERE (ZASSIM(:) >= ZDECLINE(:,1))
  !
  ! 3.1 - Growing phase : plant nitrogen decline theory
  !
  ! the growth allometric law is applied
  ! repartition of total biomass    
  !
  !before optimization
  !ZBIOMASST(:)= MAX(PBIOMASS_LEAF(:), (PBIOMASS_LEAF(:)/XCC_NITRO)**(1.0/(1.0-XCA_NIT)))  
  ZBIOMASST(:) = MAX(PBIOMASS_LEAF(:), ZWORK(:))
  !
  ! active structural biomass increment and storage
  !
  ZBIOMASS(:,2)  = ZBIOMASST(:)  - PBIOMASS_LEAF(:)
  ZDECLINE(:,2)  = ZBIOMASS(:,2) * (1.0-EXP(-1.0*XDAY/IMT%XSEFOLD(:,1)))
  ZSTORAGE(:,1)  = ZBIOMASS(:,2) - IR%XBIOMASS(:,2,1) + ZDECLINE(:,2) + IR%XRESP_BIOMASS(:,2,1)
  !
ELSE WHERE
  !
  ! 3.2 - Senescence phase
  !
  ! the active structural biomass dies exponentially at the lowest rate
  !
  ZSTORAGE(:,1) = 0.0
  ZDECLINE(:,2) = IR%XBIOMASS(:,2,1) * (1.0-EXP(-1.0*XDAY/IMT%XSEFOLD(:,1)))
  ZBIOMASS(:,2) = IR%XBIOMASS(:,2,1) - ZDECLINE(:,2) - IR%XRESP_BIOMASS(:,2,1)
  !
  !  Avoid negative values of biomass
  !  No test on ZDECLINE(:,2) as it is not used after, or recalculated
  !  No test on IR%XRESP_BIOMASS(:,2,1) as it should be smaller than IR%XBIOMASS(:,2,1)
  !  otherwise there are irrealistic values of temperature     
  !
  ZBIOMASS(:,2) = MAX(ZBIOMASS(:,2),0.0)
  !
  ZBIOMASST(:) = PBIOMASS_LEAF(:) + ZBIOMASS(:,2)
  !
END WHERE
!
! 3.3 - Flow to the passive structural biomass: cut or growth after senescence
! Biomass is taken from active structural biomass, not from senescence of leaves
! 
ZINCREASE(:,1) = ZASSIM(:)
ZINCREASE(:,2) = ZSTORAGE(:,1)
ZINCREASE(:,3) = -MIN(ZSTORAGE(:,1),0.0)
!
ZSTORAGE (:,1) = MAX(0.0,ZSTORAGE(:,1))
!
! 3.4 - Mass conservation : leaf biomass sensecence must be >= structural storage
!
WHERE( ZSTORAGE(:,1) > ZDECLINE(:,1))
  ZDECLINE(:,2)    = IR%XBIOMASS(:,2,1) * (1.0 - EXP(-1.0*XDAY/IMT%XSEFOLD(:,1)))
  ZBIOMASST(:)     = IR%XBIOMASS(:,1,1) + IR%XBIOMASS(:,2,1) - ZDECLINE(:,2) - IR%XRESP_BIOMASS(:,2,1)  
END WHERE
!
ZWORK(:) = 0.0
WHERE( ZBIOMASST(:) > 0.0)
      ZWORK(:) = (1.0-XCA_NIT)*LOG(ZBIOMASST(:))
      ZWORK(:) = EXP(ZWORK(:))
ENDWHERE
!
WHERE( ZSTORAGE(:,1) > ZDECLINE(:,1))
  !   
  !before optimization
  !PBIOMASS_LEAF(:)= ZCC_NITRO * (ZBIOMASST(:)**(1.0-XCA_NIT))
  PBIOMASS_LEAF(:) = XCC_NITRO * ZWORK(:)
  ZBIOMASS(:,2)    = ZBIOMASST(:)  - PBIOMASS_LEAF(:)
  ZDECLINE(:,1)    = IR%XBIOMASS(:,1,1) - PBIOMASS_LEAF(:)
  ZSTORAGE(:,1)    = ZBIOMASS(:,2) - IR%XBIOMASS(:,2,1) + ZDECLINE(:,2) + IR%XRESP_BIOMASS(:,2,1)  
  !
  ZINCREASE(:,2) = ZSTORAGE(:,1)
  !
END WHERE
!print*,'leaf3 ',PBIOMASS_LEAF(1)
!
!-------------------------------------------------------------------------------
!
! 4 - Evolution of other biomass pools and final calculations
! -----------------------------------------------------------
!
! 4.1 - Mortality of leaf biomass
!
ZMORT_LEAF(:) = MAX(0.0, ZDECLINE(:,1) - ZSTORAGE(:,1))
!
ZBIOMASS(:,3) = IR%XBIOMASS(:,3,1)
!
IF (IO%CPHOTO=='NIT') THEN
  !
  ! senesence of deep-structural biomass
  !
  ZDECLINE(:,3) = ZBIOMASS(:,3)*(1.0-EXP(-1.0*XDAY/IMT%XSEFOLD(:,1)))          
  !
  ! threshold value for leaf biomass and total above ground biomass in nitrogen
  ! dilution theory now in modd_co2v_par.F90 (XBIOMASST_LIM)
  !
  ! emergency deep structural biomass
  WHERE((ZBIOMASST(:) <= XBIOMASST_LIM) .AND. (ZXSEFOLD(:) > 1.0))
    ZBIOMASS(:,3) = ZBIOMASS(:,3) + ZMORT_LEAF(:)
  END WHERE
  !
ELSEIF (IO%CPHOTO=='NCB') THEN
  !
  ! 4.2 - Evolution of the other reservoirs
  ! 4.2.1 - senesence, avoiding negative values of biomass
  !
  ZDECLINE(:,3) = MIN(IR%XBIOMASS(:,3,1)*(1.0-EXP(-1.0*XDAY/(IMT%XSEFOLD(:,1)/4.))), &
                      IR%XBIOMASS(:,3,1)-IR%XRESP_BIOMASS(:,3,1))            
  ZDECLINE(:,4) = MIN(IR%XBIOMASS(:,4,1)*(1.0-EXP(-1.0*XDAY/IMT%XSEFOLD(:,1))), &
                      IR%XBIOMASS(:,4,1)-IR%XRESP_BIOMASS(:,4,1))
  !
  WHERE (OWOOD(:))
    ! Woody
    ZDECLINE(:,5) = MIN(IR%XBIOMASS(:,5,1)*(1.0-EXP(-1.0*XDAY/IP%XTAU_WOOD(:,1))), &
                        IR%XBIOMASS(:,5,1)-IR%XRESP_BIOMASS(:,5,1))
    ZDECLINE(:,6) = IR%XBIOMASS(:,6,1)*(1.0-EXP(-1.0*XDAY/IP%XTAU_WOOD(:,1)))
  ELSEWHERE
    ! Herbaceous
    ZDECLINE(:,5) = 0.
    ZDECLINE(:,6) = 0.
  END WHERE
  !
  ! 4.2.2 - storage (part of decline used as input for other reservoirs)
  !
  GMASK_ASSIM (:)=(ZASSIM(:) >= ZDECLINE(:,1))
  !
  WHERE (GMASK_ASSIM(:))
    !
    ! Remaining mortality is stored in roots.
    ZINCREASE(:,4)   = ZMORT_LEAF(:)
    !      
    ! Growing phase, all leaf decline is used as storage.
    ZSTORAGE(:,1)    = ZSTORAGE(:,1) + ZINCREASE(:,4)
    ZMORT_LEAF(:)    = ZMORT_LEAF(:) - ZINCREASE(:,4)
    !      
    ZSTORAGE(:,2)    = ZDECLINE(:,2)
    ZSTORAGE(:,3)    = ZDECLINE(:,3)
    !   
  ELSEWHERE
    !
    ! Senescence, a part of mortality is stored in roots, limited by assimilation rate.
    ZINCREASE(:,4)   = MIN(MAX(0.5*ZASSIM(:),0.) , 0.5*ZMORT_LEAF(:))
    !
    ZSTORAGE(:,1)    = ZSTORAGE(:,1) + ZINCREASE(:,4)
    ZMORT_LEAF(:)    = ZMORT_LEAF(:) - ZINCREASE(:,4)
    !   
  END WHERE
  !
  WHERE(GMASK_ASSIM(:).AND.OWOOD(:))
      ! Woody
      ZSTORAGE(:,4)  = ZDECLINE(:,4)
      !
      ZINCREASE(:,4) = ZINCREASE(:,4) + 0.3* (ZSTORAGE(:,2) + ZSTORAGE(:,3))
      ZINCREASE(:,5) =                  0.7* (ZSTORAGE(:,2) + ZSTORAGE(:,3))
      ZINCREASE(:,6) = ZSTORAGE(:,4)
      !
  ELSEWHERE(GMASK_ASSIM(:).AND..NOT.OWOOD(:))
      ! Herbaceous
      ZSTORAGE(:,4)  = 0.
      !
      ZINCREASE(:,4) = ZINCREASE(:,4) + ZSTORAGE(:,2) + ZSTORAGE(:,3)
      !
  END WHERE
  !
  WHERE (.NOT.GMASK_ASSIM(:).AND.OWOOD(:))
      ! Woody
      ! Senescence, only a part of decline is used as storage
      ZSTORAGE(:,2)  = 0.5*ZDECLINE(:,2)
      ZSTORAGE(:,3)  = 0.5*ZDECLINE(:,3)
      ZSTORAGE(:,4)  = 0.5*ZDECLINE(:,4)
      !
      ZINCREASE(:,5) = ZSTORAGE(:,2) + ZSTORAGE(:,3)
      ZINCREASE(:,6) = ZSTORAGE(:,4)
      !
  ELSEWHERE(.NOT.GMASK_ASSIM(:).AND..NOT.OWOOD(:))
      !  Herbaceous
      ! Senescence, no storage
      ZSTORAGE(:,2)  = 0.
      ZSTORAGE(:,3)  = 0.
      ZSTORAGE(:,4)  = 0.
      !
  END WHERE
  !
  ZSTORAGE(:,5) = 0.
  ZSTORAGE(:,6) = 0.
  !
  ! 4.2.3 - mortality (senescence - storage) and turnover
  !
  IF (IO%CRESPSL=='CNT') THEN
    IP%XTURNOVER(:,1,1) = ZMORT_LEAF(:)*1000.*XPCCO2/XDAY
    IP%XTURNOVER(:,2,1) = (ZDECLINE(:,2) - ZSTORAGE(:,2))*1000.*XPCCO2/XDAY
    IP%XTURNOVER(:,3,1) = (ZDECLINE(:,3) - ZSTORAGE(:,3))*1000.*XPCCO2/XDAY
    IP%XTURNOVER(:,4,1) = (ZDECLINE(:,4) - ZSTORAGE(:,4))*1000.*XPCCO2/XDAY
    IP%XTURNOVER(:,5,1) = (ZDECLINE(:,5) - ZSTORAGE(:,5))*1000.*XPCCO2/XDAY
    IP%XTURNOVER(:,6,1) = (ZDECLINE(:,6) - ZSTORAGE(:,6))*1000.*XPCCO2/XDAY
  ENDIF
  !
ENDIF
!
!
! 4.3 - Re-initialisations for next time step
!
ZBIOMASS(:,3) = ZBIOMASS(:,3) + ZINCREASE(:,3) - ZDECLINE(:,3) - IR%XRESP_BIOMASS(:,3,1)
!
! Add net accumulated CO2 assimilation 
PBIOMASS_LEAF(:) = PBIOMASS_LEAF(:) + ZASSIM(:)
!print*,'leaf4 ',PBIOMASS_LEAF(1)
!
! re-initialisation of biomass compartments values: X(day) <-- X(day-1)
IR%XBIOMASS(:,1,1) = PBIOMASS_LEAF(:)
IR%XBIOMASS(:,2,1) = ZBIOMASS(:,2)
IR%XBIOMASS(:,3,1) = ZBIOMASS(:,3)
!
! re-initialisation of respiration and assimilation terms
IR%XRESP_BIOMASS(:,2,1) = 0.0
IR%XRESP_BIOMASS(:,3,1) = 0.0
IR%XANFM(:,1) = 0.0
!
!
! 4.2.4 - evolution of reservoirs
!
IF (IO%CPHOTO=='NIT') THEN
  !
  IR%XBIOMASS(:,3,1) = MAX(IR%XBIOMASS(:,3,1),0.0)
  !
ELSEIF (IO%CPHOTO=='NCB') THEN
  !
  ZBIOMASS(:,4) = IR%XBIOMASS(:,4,1) + ZINCREASE(:,4) - ZDECLINE(:,4) - IR%XRESP_BIOMASS(:,4,1)
  !
!
  ZBIOMASS(:,5) = IR%XBIOMASS(:,5,1)
  ZBIOMASS(:,6) = IR%XBIOMASS(:,6,1)
  ZRESP(:) = IR%XRESP_BIOMASS(:,5,1)
!
  DO JSPIN = 1, IO%NSPINW
    DO JI = 1,INI
       IF(OWOOD(JI))THEN
         !Woody
         ZBIOMASS(JI,5) = ZBIOMASS(JI,5) + ZINCREASE(JI,5) - ZDECLINE(JI,5) - ZRESP(JI)
         ZBIOMASS(JI,6) = ZBIOMASS(JI,6) + ZINCREASE(JI,6) - ZDECLINE(JI,6)
         ZDECLINE(JI,5) = ZBIOMASS(JI,5)*(1.0-EXP((-1.0*XDAY)/IP%XTAU_WOOD(JI,1)))
         ZDECLINE(JI,6) = ZBIOMASS(JI,6)*(1.0-EXP((-1.0*XDAY)/IP%XTAU_WOOD(JI,1)))
         IF (IR%XBIOMASS(JI,5,1).GT.0.0) &
                 ZRESP(JI) = IR%XRESP_BIOMASS(JI,5,1)/IR%XBIOMASS(JI,5,1) * ZBIOMASS(JI,5)  
       ELSE   
         !Herbaceous
         ZBIOMASS(JI,5) = 0.
         ZBIOMASS(JI,6) = 0.
       ENDIF
    ENDDO
  ENDDO
!
  IR%XBIOMASS(:,4,1) = ZBIOMASS(:,4)
  IR%XBIOMASS(:,5,1) = ZBIOMASS(:,5)
  IR%XBIOMASS(:,6,1) = ZBIOMASS(:,6)
  !
  IR%XRESP_BIOMASS(:,4,1) = 0.0
  IR%XRESP_BIOMASS(:,5,1) = 0.0
  !
  IP%XINCREASE(:,:,1) = ZINCREASE(:,:)
!  
ENDIF
!
IF (LHOOK) CALL DR_HOOK('NITRO_DECLINE',1,ZHOOK_HANDLE)
!
END SUBROUTINE NITRO_DECLINE
