!     #########
      SUBROUTINE HYDRO(IO, P, IP, INI, IMX, IMT, IMI, AG, IR, DGEIP, DGMI, &
                       HSNOW_ISBA, OMEB, PTSTEP, PVEG, PWRMAX, PSNOW_THRUFAL, &
                       PEVAPCOR, PSUBVCOR, PSOILHCAPZ, PF2WGHT, PF2, PPS,  &
                       PIRRIG_GR, PDELHEATG, PDELHEATG_SFC, PDELPHASEG,    &
                       PDELPHASEG_SFC                         )
!     #####################################################################
!
!!****  *HYDRO*  
!!
!!    PURPOSE
!!    -------
!
!     Calculates the evolution of the water variables, i.e., the superficial
!     and deep-soil volumetric water content (wg and w2), the equivalent
!     liquid water retained in the vegetation canopy (Wr), the equivalent
!     water of the snow canopy (Ws), and also of the albedo and density of
!     the snow (i.e., SNOWALB and SNOWRHO).  Also determine the runoff and drainage
!     into the soil.
!         
!     
!!**  METHOD
!!    ------
!
!!    EXTERNAL
!!    --------
!!
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!    Noilhan and Planton (1989)
!!    Belair (1995)
!!      
!!    AUTHOR
!!    ------
!!
!!      S. Belair           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original    14/03/95 
!!                  31/08/98 (V. Masson and F. Habets) add Dumenil et Todini
!!                           runoff scheme
!!                  31/08/98 (V. Masson and A. Boone) add the third soil-water
!!                           reservoir (WG3,D3)
!!                  19/07/05 (P. LeMoigne) bug in runoff computation if isba-2L
!!                  10/10/05 (P. LeMoigne) bug in hydro-soil calling sequence
!!                  25/05/08 (B. Decharme) Add floodplains
!!                  27/11/09 (A. Boone)    Add possibility to do time-splitting when
!!                                         calling hydro_soildif (DIF option only)
!!                                         for *very* large time steps (30min to 1h+).
!!                                         For *usual* sized time steps, time step
!!                                         NOT split.
!!                     08/11 (B. Decharme) DIF optimization
!!                     09/12 (B. Decharme) Bug in wg2 ice energy budget
!!                     10/12 (B. Decharme) EVAPCOR snow correction in DIF
!!                                         Add diag IRRIG_FLUX
!!                     04/13 (B. Decharme) Pass soil phase changes routines here
!!                                         Apply physical limits on wg in hydro_soil.F90
!!                                         Subsurface runoff if SGH (DIF option only)
!!                                         water table / surface coupling
!!                  02/2013  (C. de Munck) specified irrigation rate of ground added
!!                  10/2014  (A. Boone)    MEB added
!!                  07/15    (B. Decharme) Numerical adjustement for F2 soilstress function
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_FIX_t, ISBA_PARAM_TIME_t, ISBA_PARAM_IRRIG_t
USE MODD_ISBA_PGD_n, ONLY : ISBA_PGD_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t, ISBA_INIT_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
USE MODD_AGRI_n, ONLY : AGRI_t
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_MISC_ISBA_n, ONLY : DIAG_MISC_ISBA_t
!
USE MODD_CSTS,      ONLY : XRHOLW, XDAY, XTT, XLVTT, XLSTT, XLMTT
USE MODD_ISBA_PAR,  ONLY : XWGMIN, XDENOM_MIN
USE MODD_SURF_PAR,  ONLY : XUNDEF, NUNDEF
!
#ifdef TOPD
USE MODD_COUPLING_TOPD, ONLY : LCOUPL_TOPD, XAS_NATURE, XATOP, XRUNOFF_TOP, NMASKT_PATCH
#endif
!
USE MODI_HYDRO_VEG
USE MODI_HYDRO_SNOW
USE MODI_HYDRO_SOIL
USE MODI_HYDRO_SOILDIF                                          
USE MODI_HYDRO_SGH
USE MODI_ICE_SOILDIF              
USE MODI_ICE_SOILFR              
!
USE MODE_THERMOS
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_PARAM_FIX_t), INTENT(INOUT) :: IMX
TYPE(ISBA_PARAM_TIME_t), INTENT(INOUT) :: IMT
TYPE(ISBA_PARAM_IRRIG_t), INTENT(INOUT) :: IMI
TYPE(ISBA_PGD_t), INTENT(INOUT) :: P
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
TYPE(ISBA_INIT_t), INTENT(INOUT) :: INI
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
TYPE(AGRI_t), INTENT(INOUT) :: AG
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEIP
TYPE(DIAG_MISC_ISBA_t), INTENT(INOUT) :: DGMI
!
 CHARACTER(LEN=*),     INTENT(IN)  :: HSNOW_ISBA ! 'DEF' = Default F-R snow scheme
!                                               !         (Douville et al. 1995)
!                                               ! '3-L' = 3-L snow scheme (option)
!                                               !         (Boone and Etchevers 2001)
LOGICAL, INTENT(IN)                :: OMEB   ! True  = patch with multi-energy balance 
!                                            ! False = patch with classical (composite) ISBA
REAL, INTENT(IN)                    :: PTSTEP
!                                      timestep of the integration
!
REAL, DIMENSION(:), INTENT(IN)    :: PVEG, PWRMAX
!                                      PVEG = fraction of vegetation 
!                                      PWRMAX = maximum equivalent water content
!                                               in the vegetation canopy
!
REAL, DIMENSION(:), INTENT(IN)    :: PSNOW_THRUFAL, PEVAPCOR, PSUBVCOR
!                                    PSNOW_THRUFAL = rate that liquid water leaves snow pack: 
!                                               *ISBA-ES* [kg/(m2 s)]
!                                    PEVAPCOR = correction if evaporation from snow exceeds
!                                               actual amount on the surface [kg/(m2 s)]
!                                    PSUBVCOR = correction if sublimation from snow intercepted 
!                                               on the MEB canopy exceeds snow available as it 
!                                               disappears [kg/(m2 s)]
!
REAL, DIMENSION(:), INTENT(IN)    :: PPS, PF2                                       
!                                    PPS  = surface pressure (Pa)
!                                    PF2  = total water stress factor (-)
!
REAL, DIMENSION(:,:), INTENT(IN)  :: PF2WGHT
!                                    PF2WGHT   = water stress factor (profile) (-)
!
REAL, DIMENSION(:,:), INTENT(IN) :: PSOILHCAPZ
!                                   PSOILHCAPZ = ISBA-DF Soil heat capacity profile [J/(m3 K)]
!
REAL, DIMENSION(:), INTENT(INOUT) :: PDELHEATG, PDELHEATG_SFC
!                                      PDELHEATG_SFC = change in heat storage of the surface soil layer over the current time step (W m-2)
!                                      PDELHEATG     = change in heat storage of the entire soil column over the current time step (W m-2)
!
REAL, DIMENSION(:), INTENT(OUT)   :: PDELPHASEG, PDELPHASEG_SFC
!                                     PDELPHASEG     = latent heating due to soil freeze-thaw in the entire soil column  (W m-2)
!                                     PDELPHASEG_SFC = latent heating due to soil freeze-thaw in the surface soil layer  (W m-2)
!
REAL   ,DIMENSION(:),INTENT(IN)    :: PIRRIG_GR ! ground irrigation rate (kg/m2/s)
!
!*      0.2    declarations of local variables
!
!
INTEGER                         :: JJ, JL      ! loop control                                       
INTEGER                         :: INDT, JDT   ! Time splitting indicies
INTEGER                         :: INJ, INL, IDEPTH ! (ISBA-DF option)
!
REAL                            :: ZTSTEP      ! maximum time split time step (<= PTSTEP)
!                                              ! ONLY used for DIF option.
!
REAL, DIMENSION(SIZE(PVEG))     :: ZPG, ZPG_MELT, ZDUNNE,                            &
                                   ZLEV, ZLEG, ZLEGI, ZLETR, ZPSNV,                  &
                                   ZRR, ZDG3, ZWG3, ZWSAT_AVG, ZWWILT_AVG, ZWFC_AVG, &
                                   ZRUNOFF, ZDRAIN, ZHORTON, ZEVAPCOR, ZQSB 
!                                      Prognostic variables of ISBA at 't-dt'
!                                      ZPG = total water reaching the ground
!                                      ZPG_MELT = snowmelt reaching the ground 
!                                      ZDUNNE  = Dunne runoff
!                                 ZLEV, ZLEG, ZLEGI, ZLETR = Evapotranspiration amounts
!                                      from the non-explicit snow area *ISBA-ES*
!                                 ZPSNV = used to calculate interception of liquid
!                                      water by the vegetation in FR snow method:
!                                      For ES snow method, precipitation already modified
!                                      so set this to zero here for this option.
!                                 ZWSAT_AVG, ZWWILT_AVG, ZWFC_AVG = Average water and ice content
!                                      values over the soil depth D2 (for calculating surface runoff)
!                                 ZDG3, ZWG3, ZRUNOFF, ZDRAIN, ZQSB and ZHORTON are working variables only used for DIF option
!                                 ZEVAPCOR = correction if evaporation from snow exceeds
!                                               actual amount on the surface [m/s]
!
REAL, DIMENSION(SIZE(PVEG))     :: ZDWGI1, ZDWGI2, ZKSFC_IVEG
!                                      ZDWGI1 = surface layer liquid water equivalent 
!                                               volumetric ice content time tendency
!                                      ZDWGI2 = deep-soil layer liquid water equivalent 
!                                               volumetric ice content time tendency
!                                      ZKSFC_IVEG = non-dimensional vegetation insolation coefficient
!
REAL, DIMENSION(SIZE(PVEG))    :: ZWGI_EXCESS, ZF2
!                                 ZWGI_EXCESS = Soil ice excess water content
!                                 ZF2         = Soilstress function for transpiration
!
REAL, DIMENSION(SIZE(IR%XWG,1),SIZE(IR%XWG,2)) :: ZQSAT, ZQSATI, ZTI, ZPS
!                                           For specific humidity at saturation computation (ISBA-DIF)
!
REAL, DIMENSION(SIZE(IR%XWG,1),SIZE(IR%XWG,2)) :: ZWGI0
!                                      ZWGI0 = initial soil ice content (m3 m-3) before update
!                                              for budget diagnostics     
!
!*      0.3    declarations of local parameters
!
REAL, PARAMETER             :: ZINSOLFRZ_VEG = 0.20  ! (-)       Vegetation insolation coefficient
!
REAL, PARAMETER             :: ZINSOLFRZ_LAI = 30.0  ! (m2 m-2)  Vegetation insolation coefficient

REAL, PARAMETER             :: ZTIMEMAX      = 900.  ! s  Maximum timescale without time spliting
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!*       0.     Initialization:
!               ---------------
!
IF (LHOOK) CALL DR_HOOK('HYDRO',0,ZHOOK_HANDLE)
!
JDT    = 0
INDT   = 0
ZTSTEP = 0.0
!
ZPG(:)           = 0.0
ZPG_MELT(:)      = 0.0
ZDUNNE(:)        = 0.0
!
ZWSAT_AVG(:)     = 0.0
ZWWILT_AVG(:)    = 0.0
ZWFC_AVG(:)      = 0.0
!
ZRR(:)           = DGMI%XRRSFC(:)
!
ZDRAIN(:)        = 0.
ZHORTON(:)       = 0.
ZRUNOFF(:)       = 0.
ZWGI_EXCESS(:)   = 0.
ZEVAPCOR(:)      = 0.
ZQSB    (:)      = 0.
!
DGEIP%XDRAIN(:)        = 0.
DGEIP%XRUNOFF(:)       = 0.
DGEIP%XHORT(:)       = 0.
DGEIP%XQSB   (:)       = 0.
!
PDELPHASEG(:)    = 0.0
PDELPHASEG_SFC(:)= 0.0
ZWGI0(:,:)       = 0.0
!
ZF2(:) = MAX(XDENOM_MIN,PF2(:))
!
! Initialize evaporation components: variable definitions
! depend on snow or explicit canopy scheme:
!
IF(OMEB)THEN
!
! MEB uses explicit snow scheme by default, but fluxes already aggregated
! for snow and floods so no need to multiply by fractions here. 
!
   ZLEV(:)          = DGEIP%XLEV(:)
   ZLETR(:)         = DGEIP%XLETR(:)
   ZLEG(:)          = DGEIP%XLEG(:)
   ZLEGI(:)         = DGEIP%XLEGI(:)
   ZPSNV(:)         = 0.0
!
   ZEVAPCOR(:)      = PEVAPCOR(:) + PSUBVCOR(:)
!
ELSE
!
! Initialize evaporation components: variable definitions
! depend on snow scheme:
!
   IF(HSNOW_ISBA == '3-L' .OR. HSNOW_ISBA == 'CRO' .OR. IO%CISBA == 'DIF')THEN
      ZLEV(:)          = (1.0-IR%XPSNV(:,1)-INI%XFFV(:,1)) * DGEIP%XLEV(:)
      ZLETR(:)         = (1.0-IR%XPSNV(:,1)-INI%XFFV(:,1)) * DGEIP%XLETR(:)
      ZLEG(:)          = (1.0-IR%XPSNG(:,1)-INI%XFFG(:,1)) * DGEIP%XLEG(:)
      ZLEGI(:)         = (1.0-IR%XPSNG(:,1)-INI%XFFG(:,1)) * DGEIP%XLEGI(:)
      ZPSNV(:)         = 0.0
   ELSE
      ZLEV(:)          = DGEIP%XLEV(:)
      ZLETR(:)         = DGEIP%XLETR(:)
      ZLEG(:)          = DGEIP%XLEG(:)
      ZLEGI(:)         = DGEIP%XLEGI(:)
      ZPSNV(:)         = IR%XPSNV(:,1)+INI%XFFV(:,1)
   ENDIF
!
   ZEVAPCOR(:)         = PEVAPCOR(:) 

ENDIF
!
! Initialize average soil hydrological parameters
! over the entire soil column: if Isba Force-Restore
! is in use, then parameter profile is constant
! so simply use first element of this array: if
! the Diffusion option is in force, the relevant
! calculation is done later within this routine.
!
IF(IO%CISBA == '2-L' .OR. IO%CISBA == '3-L')THEN  
   ZWSAT_AVG(:)     = IP%XWSAT(:,1)
   ZWWILT_AVG(:)    = IP%XWWILT(:,1)
   ZWFC_AVG(:)      = IP%XWFC(:,1)
ENDIF
!
IF (IO%CISBA == '3-L') THEN                                   
   ZDG3(:) = IMX%XDG(:,3,1)
   ZWG3(:) = IR%XWG (:,3,1)
ELSE
   ZDG3(:) = XUNDEF
   ZWG3(:) = XUNDEF
END IF
!
!-------------------------------------------------------------------------------
!
!*       1.     EVOLUTION OF THE EQUIVALENT WATER CONTENT Wr
!               --------------------------------------------
!
!
!
IF(.NOT.OMEB)THEN ! Canopy Int & Irrig Already accounted for if MEB in use.
!
   DGEIP%XIRRIG_FLUX(:)=0.0
!
!* add irrigation over vegetation to liquid precipitation (rr)
!
!
   IF (ASSOCIATED(AG%LIRRIGATE)) THEN
     IF (SIZE(AG%LIRRIGATE,1)>0) THEN
       WHERE (AG%LIRRIGATE(:,1) .AND. IMI%XIRRIG(:,1)>0. .AND. &
             IMI%XIRRIG(:,1) /= XUNDEF .AND. (PF2(:)<AG%XTHRESHOLDSPT(:,1)) )
         DGEIP%XIRRIG_FLUX(:) = IMI%XWATSUP(:,1) / XDAY           
         ZRR     (:) = ZRR(:) + IMI%XWATSUP(:,1) / XDAY
         AG%LIRRIDAY(:,1)    = .TRUE.           
       END WHERE
     ENDIF
   ENDIF
!
!* interception reservoir and dripping computation
!
   CALL HYDRO_VEG(IO%CRAIN, PTSTEP, INI%XMUF, ZRR, ZLEV, ZLETR, PVEG, &
                  ZPSNV,  IR%XWR(:,1), PWRMAX, ZPG, DGEIP%XDRIP, DGEIP%XRRVEG  ) 
!
!
!
ELSE
!
! For MEB case, interception interactions already computed and DGMI%XRRSFC represents
! water falling (drip and not intercepted by vegetation) outside of snow covered
! areas. Part for snow covered areas (net outflow at base of snowpack) accounted
! for in PSNOW_THRUFAL.
!
   ZPG(:) = DGMI%XRRSFC(:)
!
ENDIF
!
!* add irrigation over ground to potential soil infiltration (pg)
!
DGEIP%XIRRIG_FLUX(:) = DGEIP%XIRRIG_FLUX(:) + PIRRIG_GR(:)
!
ZPG(:) = ZPG(:) + PIRRIG_GR(:)
!
!-------------------------------------------------------------------------------
!
!*       2.     EVOLUTION OF THE EQUIVALENT WATER CONTENT snowSWE 
!               -------------------------------------------------
!
!*       3.     EVOLUTION OF SNOW ALBEDO 
!               ------------------------
!
!*       4.     EVOLUTION OF SNOW DENSITY 
!               -------------------------
!
! Boone and Etchevers '3-L' snow option
IF(HSNOW_ISBA == '3-L' .OR. HSNOW_ISBA == 'CRO' .OR. IO%CISBA == 'DIF')THEN
!
  ZPG_MELT(:)   = ZPG_MELT(:)   + PSNOW_THRUFAL(:)          ! [kg/(m2 s)]
!
! Note that 'melt' is referred to as rain and meltwater
! running off from the snowpack in a timestep for ISBA-ES,
! not the actual amount of ice converted to liquid.
!
  DGEIP%XMELT(:) = DGEIP%XMELT(:) + PSNOW_THRUFAL(:)          ! [kg/(m2 s)]
!
ELSE
  !
  CALL HYDRO_SNOW(IO%LGLACIER, PTSTEP, IP%XVEGTYPE_PATCH(:,:,1), DGMI%XSRSFC, &
                  DGEIP%XLES, DGEIP%XMELT, IR%TSNOW, ZPG_MELT)
  !
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       5.     Sub Grid Hydrology
!               ------------------
!
! - Dunne runoff  : Dumenil et Todini (1992) or Topmodel
! - Horton runoff : Direct or exponential precipitation distribution
! - Floodplains interception and infiltration
!
 CALL HYDRO_SGH(IO, P, IP, INI, IMX, IR, DGEIP, DGMI, PTSTEP, ZPG, ZPG_MELT, ZDUNNE )         
!
!-------------------------------------------------------------------------------
!
!*       6.     EVOLUTION OF THE SOIL WATER CONTENT
!               -----------------------------------
!
!*       7.     EFFECT OF MELTING/FREEZING ON SOIL ICE AND LIQUID WATER CONTENTS
!               ----------------------------------------------------------------
!
!*       8.     DRAINAGE FROM THE DEEP SOIL
!               ---------------------------
!
!*      9.     RUN-OFF 
!               -------
!                                     when the soil water exceeds saturation, 
!                                     there is fast-time-response runoff
!
!
! -----------------------------------------------------------------
! Time splitting parameter for *very large time steps* since Richard
! and/or soil freezing equations are very non-linear 
! NOTE for NWP/GCM type applications, the time step is generally not split
! (usually just for offline applications with a time step on order of 
! 15 minutes to an hour for example)
! ------------------------------------------------------------------
!
INDT = 1
IF(PTSTEP>=ZTIMEMAX)THEN
  INDT = MAX(2,NINT(PTSTEP/ZTIMEMAX))
ENDIF
!
ZTSTEP  = PTSTEP/REAL(INDT)
!
! ------------------------------------------------------------------
! The values for the two coefficients (multiplied by VEG and LAI) 
! in the expression below are from 
! Giard and Bazile (2000), Mon. Wea. Rev.: they model the effect of insolation due to
! vegetation cover. This used by both 'DEF' (code blocks 3.-4.) and 'DIF' options.
! ------------------------------------------------------------------
!
WHERE(IMT%XLAI(:,1)/=XUNDEF .AND. PVEG(:)/=0.)
    ZKSFC_IVEG(:) = (1.0-ZINSOLFRZ_VEG*PVEG(:)) * (1.0-(IMT%XLAI(:,1)/ZINSOLFRZ_LAI))
ELSEWHERE
    ZKSFC_IVEG(:) = 1.0 ! No vegetation
ENDWHERE
!
!
ZWGI0 (:,:) = IR%XWGI(:,:,1) ! save initial ice content before phase changes and sublimation
!
IF (IO%CISBA=='DIF') THEN                
!
  INJ = SIZE(IMX%XDG(:,:,1),1)
  INL = MAXVAL(IMX%NWG_LAYER(:,1))
!
! Initialize some field
! ---------------------
!
  ZPS(:,:)=XUNDEF
  ZTI(:,:)=XUNDEF
  DO JL=1,INL
     DO JJ=1,INJ
        IDEPTH=IMX%NWG_LAYER(JJ,1)
        IF(JL<=IDEPTH)THEN
          ZPS(JJ,JL) = PPS(JJ)
          ZTI(JJ,JL) = MIN(XTT,IR%XTG(JJ,JL,1))
        ENDIF
     ENDDO
  ENDDO
!
! Compute specific humidity at saturation for the vapor conductivity
! ------------------------------------------------------------------
!
  ZQSAT (:,:) = QSAT (IR%XTG(:,:,1),ZPS(:,:),IMX%NWG_LAYER(:,1),INL)
  ZQSATI(:,:) = QSATI(ZTI(:,:),ZPS(:,:),IMX%NWG_LAYER(:,1),INL)
!
! Soil water sink terms: convert from (W m-2) and (kg m-2 s-1) to (m s-1)
! ------------------------------------------------------------------
!
  ZPG     (:) =  ZPG    (:)        / XRHOLW
  ZEVAPCOR(:) = PEVAPCOR(:)        / XRHOLW
  ZLEG    (:) =  ZLEG   (:)        /(XRHOLW*XLVTT)
  ZLETR   (:) = (ZLETR  (:)/ZF2(:))/(XRHOLW*XLVTT)
  ZLEGI   (:) = ZLEGI   (:)        /(XRHOLW*XLSTT)
!
  DO JDT = 1,INDT
!                      
    CALL HYDRO_SOILDIF(IO, IP, INI, IMX, IR, ZTSTEP, ZPG, ZLETR, ZLEG, ZEVAPCOR,  &
                       PF2WGHT, PPS, ZQSAT, ZQSATI, ZDRAIN, ZHORTON, INL, ZQSB )
!
    CALL ICE_SOILDIF(IP, IMX, IR, ZTSTEP, ZKSFC_IVEG, ZLEGI, PSOILHCAPZ, ZWGI_EXCESS )
!
    DGEIP%XDRAIN(:) = DGEIP%XDRAIN(:) + (ZDRAIN(:)+ZQSB(:)+ZWGI_EXCESS(:))/REAL(INDT)
    DGEIP%XQSB  (:) = DGEIP%XQSB  (:) + ZQSB   (:)/REAL(INDT)
    DGEIP%XHORT (:) = DGEIP%XHORT (:) + ZHORTON(:)/REAL(INDT)
!
! Output diagnostics:
! Compute latent heating from phase change only in surface layer and total soil column,
! then adjust surface and total soil heat content to maintain balance.
!
    PDELPHASEG_SFC(:) = (IR%XWGI(:,1,1)-ZWGI0(:,1))*(XLMTT*XRHOLW/PTSTEP)*IP%XDZG(:,1,1) + &
                         ZLEGI(:)*(XRHOLW*XLSTT)
    PDELPHASEG(:)     = PDELPHASEG_SFC(:)
    DO JL=2,INL
      DO JJ=1,INJ
        PDELPHASEG(JJ) = PDELPHASEG(JJ) + (IR%XWGI(JJ,JL,1)-ZWGI0(JJ,JL))*&
                           (XLMTT*XRHOLW/PTSTEP)*IP%XDZG(JJ,JL,1)
      ENDDO
    ENDDO
    PDELHEATG_SFC(:)     = PDELHEATG_SFC(:) + PDELPHASEG_SFC(:)
    PDELHEATG(:)         = PDELHEATG(:)     + PDELPHASEG(:)

  ENDDO
!
ELSE
!
  DO JDT = 1,INDT
!
!   Only layer 1 and 2 are used for soil freezing (ZWG3 not used)
    CALL ICE_SOILFR(IO, IP, INI, IMX, IR, DGMI, HSNOW_ISBA, ZTSTEP, &
                    ZKSFC_IVEG, ZDWGI1, ZDWGI2)
!
    CALL HYDRO_SOIL(IO, P, IP, IMX, IR, DGMI, ZTSTEP,    &
                    ZLETR, ZLEG, ZPG, ZEVAPCOR, ZDG3,    &
                    ZWSAT_AVG, ZWFC_AVG, ZDWGI1, ZDWGI2, &
                    ZLEGI, ZWG3, ZRUNOFF, ZDRAIN, ZWWILT_AVG )
!
    DGEIP%XDRAIN (:)  = DGEIP%XDRAIN (:) + ZDRAIN (:)/REAL(INDT)
    DGEIP%XRUNOFF(:)  = DGEIP%XRUNOFF(:) + ZRUNOFF(:)/REAL(INDT)
!    
  ENDDO
!
! Output diagnostics:
! Compute latent heating from phase change only in surface layer and total soil column,
! then adjust surface and total soil heat content to maintain balance.
!
  PDELPHASEG_SFC(:)    = (IR%XWGI(:,1,1)-ZWGI0(:,1))*(XLMTT*XRHOLW/PTSTEP)*IMX%XDG(:,1,1) + ZLEGI(:)
  PDELPHASEG(:)        = (IR%XWGI(:,2,1)-ZWGI0(:,2))*(XLMTT*XRHOLW/PTSTEP)*IMX%XDG(:,2,1)
  PDELHEATG_SFC(:)     = PDELHEATG_SFC(:) + PDELPHASEG_SFC(:)
  PDELHEATG(:)         = PDELHEATG(:)     + PDELPHASEG(:)
!
  IF (IO%CISBA == '3-L') IR%XWG(:,3,1) = ZWG3(:)
!
#ifdef TOPD
  IF (LCOUPL_TOPD) THEN
    !runoff topo cumule (kg/m²)
    DO JJ=1,SIZE(NMASKT_PATCH)
      IF  (NMASKT_PATCH(JJ)/=0) THEN
        IF ( XATOP(NMASKT_PATCH(JJ))/=XUNDEF) THEN
          XRUNOFF_TOP(NMASKT_PATCH(JJ)) = XRUNOFF_TOP(NMASKT_PATCH(JJ)) + &
                                          (DGEIP%XRUNOFF(JJ)+ DGEIP%XHORT(JJ))*XATOP(NMASKT_PATCH(JJ))*PTSTEP
          IF (IO%CRUNOFF=='TOPD') THEN
            XRUNOFF_TOP(NMASKT_PATCH(JJ)) = XRUNOFF_TOP(NMASKT_PATCH(JJ)) + ZDUNNE(JJ)*PTSTEP
          ELSE
            ! ZDUNNE contains only saturated pixels on mesh so only catchment
            XRUNOFF_TOP(NMASKT_PATCH(JJ)) = XRUNOFF_TOP(NMASKT_PATCH(JJ)) + &
                                            ZDUNNE(JJ)*XATOP(NMASKT_PATCH(JJ))*PTSTEP
          ENDIF  
        ENDIF
      ENDIF
      ! ZDUNNE concerns all the mesh so not only catchment =>*XATOP
    ENDDO
  ENDIF
#endif
  !
ENDIF
!
!-------------------------------------------------------------------------------
!
! Add sub-grid surface and subsurface runoff to saturation excess:
!
DGEIP%XRUNOFF(:) = DGEIP%XRUNOFF(:) + ZDUNNE(:) + DGEIP%XHORT(:)
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('HYDRO',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE HYDRO
