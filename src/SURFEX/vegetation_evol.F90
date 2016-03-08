!     #########
    SUBROUTINE VEGETATION_EVOL(IO, IP, IMX, IMT, IMA, IMI, IR, OAGRIP, PTSTEP, KMONTH, KDAY, PTIME, &
                               PLAT, PRHOA, P_CO2, MSS, PRESP_BIOMASS_INST, PSWDIR)  
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
!!      V. Masson          * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/03/03 
!!      P. Le Moigne 12/2004 : NIT version 
!!      P Le Moigne  09/2005 : AGS modifs of L. Jarlan
!!      A.L. Gibelin 04/2009 : BIOMASS and RESP_BIOMASS arrays
!!      A.L. Gibelin 04/2009 : Add NCB option 
!!      D. Carrer    01/2012 : representation of nitrogen dilution fct of CO2 (from Calvet et al. 2008)
!!      B. Decharme  05/2012 : Optimization and ISBA-DIF coupling
!!      C. Delire    01/2014 : IBIS respiration for tropical evergreen
!!      R. Seferian  05/2015 : expanding of Nitrogen dilution option to the complete formulation proposed by Yin et al. GCB 2002 
!!Seferian & Delire  06/2015 : accouting for living woody biomass respiration (expanding work of E Joetzjer to all woody PFTs) 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_FIX_t, ISBA_PARAM_TIME_t, ISBA_PARAM_ALB_t, &
                              ISBA_PARAM_IRRIG_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
!
USE MODD_SSO_n, ONLY : SSO_t
!
USE MODD_CO2V_PAR,       ONLY : XMC, XMCO2, XPCCO2, XRESPFACTOR_NIT,       &
                                XCOEFF_MAINT_RESP_ZERO, XSLOPE_MAINT_RESP, &
                                XPARAM, XPARCF, XDILUDEC
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
USE MODD_DATA_COVER_PAR, ONLY : NVT_TEBD, NVT_TRBE, NVT_BONE,   &
                                NVT_TRBD, NVT_TEBE, NVT_TENE,   &
                                NVT_BOBD, NVT_BOND, NVT_SHRB,   &
                                NVT_TRBE, NVT_C3, NVT_C4,       &
                                NVT_IRR, NVT_GRAS
!
USE MODD_SURF_PAR
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
TYPE(ISBA_PARAM_FIX_t), INTENT(INOUT) :: IMX
TYPE(ISBA_PARAM_TIME_t), INTENT(INOUT) :: IMT
TYPE(ISBA_PARAM_ALB_t), INTENT(INOUT) :: IMA
TYPE(ISBA_PARAM_IRRIG_t), INTENT(INOUT) :: IMI
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
!
LOGICAL,              INTENT(IN)    :: OAGRIP  ! agricultural practices
!
REAL,                 INTENT(IN)    :: PTSTEP  ! time step
INTEGER,              INTENT(IN)    :: KMONTH  ! current month
INTEGER,              INTENT(IN)    :: KDAY    ! current day
REAL,                 INTENT(IN)    :: PTIME   ! current time since midnight
REAL,   DIMENSION(:), INTENT(IN)    :: PLAT    ! latitude of each grid point
REAL,   DIMENSION(:), INTENT(IN)    :: PRHOA   ! air density
!
REAL,   DIMENSION(:), INTENT(IN)    :: P_CO2 ! CO2 concentration [ppmm]
!
TYPE(SSO_t), INTENT(INOUT) :: MSS
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PRESP_BIOMASS_INST ! instantaneous respiration of biomass (kgCO2/kgair m/s)
!
REAL, DIMENSION(:),   INTENT(IN),   OPTIONAL :: PSWDIR    ! Global incoming shortwave radiation (W m-2)
!
!*      0.2    declarations of local parameter
!
REAL, PARAMETER                   :: ZCOEF1 = 10.0
REAL, PARAMETER                   :: ZCOEF2 = 25.0
REAL, PARAMETER                   :: ZDEPTH = 1.0   !Temp depth m
!
REAL, PARAMETER                   :: ZWOOD_IBIS=0.0125
REAL, PARAMETER                   :: ZROOT_IBIS=1.25 
REAL, PARAMETER                   :: ZCIBIS1   =3500.
REAL, PARAMETER                   :: ZCIBIS2   =1./288.
REAL, PARAMETER                   :: ZNDAY     =365.
!
REAL, PARAMETER                   :: ZCDILU1 = -0.048
REAL, PARAMETER                   :: ZCDILU2 = 6.3
REAL, PARAMETER                   :: ZCDILU3 = 371.
! Required for Yin et al., nitrogen dilu param
REAL, PARAMETER                   :: ZPHOTON    = 2.010402e-3 ! conversion coef for W m-2 in photon m-2
REAL, PARAMETER                   :: ZDEPTH_VEG = 0.40        !Depth in meters for daily temperature
REAL, PARAMETER                   :: ZTEMP_VEG  = 23.         !Average temperature of the vegetation
REAL, PARAMETER                   :: ZDECIDUS   = 0.75        !Coef for decidus trees
!
!*      0.3    declarations of local variables
!
REAL, DIMENSION(SIZE(IR%XRESP_BIOMASS,1),SIZE(IR%XRESP_BIOMASS,2)) :: ZRESP_BIOMASS_LAST ! biomass at t-1 (kg_DM/m2/day)
REAL,    DIMENSION(SIZE(IMT%XLAI,1))    :: ZBIOMASS_LEAF   ! temporary leaf biomass 
REAL,    DIMENSION(SIZE(IMT%XLAI,1))    :: ZBSLAI_NITRO    ! (Calvet et al. 2008) ratio of biomass to LAI
                                                     ! with representation of nitrogen dilution
REAL,    DIMENSION(SIZE(IMT%XLAI,1)) :: ZCO2, ZCNA_NITRO   ! fct of CO2        
REAL,    DIMENSION(SIZE(IMT%XLAI,1)) :: ZPARAM
REAL,    DIMENSION(SIZE(IMT%XLAI,1)) :: ZHTREE, ZSAPFRAC   ! tree height & sap fraction used for estimation of 
                                                     ! sapwood fraction
!
REAL                              :: ZLOG2, ZWORK
!
REAL, DIMENSION(SIZE(IR%XTG,1))      :: ZTG_VEG      ! surface temperature   (C)
REAL, DIMENSION(SIZE(IR%XTG,1))      :: ZTG_SOIL     ! soil temperature   (C)
REAL, DIMENSION(SIZE(IR%XTG,1))      :: ZDG_SOIL     ! soil depth for DIF (m)
REAL                              :: ZWGHT_SOIL   ! Weight for DIF (m)
!
LOGICAL, DIMENSION(SIZE(IMT%XLAI,1))    :: GWOOD,GHERB
LOGICAL, DIMENSION(SIZE(IMT%XLAI,1))    :: GMASK_AGRI
LOGICAL                           :: GMASK
INTEGER                           :: INI, INL, JI, JL, IDEPTH, JTYPE
!
REAL,    DIMENSION(SIZE(IP%XVEGTYPE_PATCH,1),SIZE(IP%XVEGTYPE_PATCH,2)) :: ZPARAM_TYPE
!
! * Azote
REAL,    DIMENSION(SIZE(IMT%XLAI,1)) :: ZFERT
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-----------------------------------------------------------------
!
!*      1.     Preliminaries
!              -------------
!
IF (LHOOK) CALL DR_HOOK('VEGETATION_EVOL',0,ZHOOK_HANDLE)
!
INI=SIZE(IR%XTG,1)
INL=SIZE(IR%XTG,2)
!
ZLOG2 = LOG(2.0)
!
ZTG_SOIL(:) = 0.0
ZTG_VEG (:) = 0.0
!
! Define herbaceous and woody patches
GHERB(:) = ( IP%XVEGTYPE_PATCH(:,NVT_TEBD,1) + IP%XVEGTYPE_PATCH(:,NVT_TRBE,1) + IP%XVEGTYPE_PATCH(:,NVT_BONE,1)    &
&          + IP%XVEGTYPE_PATCH(:,NVT_TRBD,1) + IP%XVEGTYPE_PATCH(:,NVT_TEBE,1) + IP%XVEGTYPE_PATCH(:,NVT_TENE,1)    &
&          + IP%XVEGTYPE_PATCH(:,NVT_BOBD,1) + IP%XVEGTYPE_PATCH(:,NVT_BOND,1) + IP%XVEGTYPE_PATCH(:,NVT_SHRB,1)<0.5)
GWOOD(:) = (.NOT.GHERB (:))
!
! Mask where vegetation evolution is performed (just before solar midnight)
GMASK = ( PTIME - PTSTEP < 0. ) .AND. ( PTIME >= 0. )
!
! Save RESP_BIOMASS at t-1
IF (GMASK) THEN
  IR%XRESP_BIOMASS(:,1,1) = 0.0
  ZRESP_BIOMASS_LAST(:,:) = 0.0
ELSE
  IR%XRESP_BIOMASS(:,1,1) = IR%XRESP_BIOMASS(:,1,1) + PRESP_BIOMASS_INST(:,1) * (PTSTEP*PRHOA(:)*XMC)/(XPCCO2*XMCO2)
  ZRESP_BIOMASS_LAST(:,:) = IR%XRESP_BIOMASS(:,:,1)
ENDIF
!
!*      2.     Interactive vegetation
!              ----------------------
!
!  LAI daily mortality and assimilation
!
ZBIOMASS_LEAF(:) = IR%XBIOMASS(:,1,1)
!
IF (GMASK) THEN
!        
  IF (IO%CPHOTO=='LAI' .OR. IO%CPHOTO=='LST')THEN
!
    CALL LAILOSS(IMT, IP, IR, ZBIOMASS_LEAF)  
    CALL LAIGAIN(IMT%XBSLAI(:,1), IMT, IR, ZBIOMASS_LEAF)
    IR%XBIOMASS(:,1,1) = ZBIOMASS_LEAF(:)
!    
  ELSE IF (IO%CPHOTO=='NIT' .OR. IO%CPHOTO=='NCB') THEN
!    
    IP%XINCREASE (:,:,1) = 0.0
    IP%XTURNOVER(:,:,1) = 0.0
    ZBSLAI_NITRO(:  ) = IP%XBSLAI_NITRO(:,1) 
!
    IF(IO%LNITRO_DILU)THEN
!
!     * Compute Vegetation temperature
!       We use the temperature of the second layer of the soil (<40cm)
!       since the parametrization employs a daily temperature
!
      IF(IO%CISBA/='DIF')THEN        
        ZTG_VEG(:) = IR%XTG(:,2,1)
      ELSE 
        DO JI=1,INI
           IDEPTH=IMX%NWG_LAYER(JI,1)
           ZDG_SOIL(JI)=MIN(ZDEPTH_VEG,IMX%XDG(JI,IDEPTH,1))
        ENDDO  
        DO JL=1,INL
           DO JI=1,INI     
              ZWGHT_SOIL=MIN(IP%XDZG(JI,JL,1),MAX(0.0,ZDG_SOIL(JI)-IMX%XDG(JI,JL,1)+IP%XDZG(JI,JL,1)))        
              ZTG_VEG(JI)=ZTG_VEG(JI)+IR%XTG(JI,JL,1)*ZWGHT_SOIL/ZDG_SOIL(JI)
           ENDDO
        ENDDO 
      ENDIF
!
      ZPARAM(:) = 0.0
      ZFERT (:) = 0.0
      DO JTYPE=1,SIZE(IP%XVEGTYPE_PATCH,2)
        DO JI = 1,INI
            ZPARAM_TYPE(JI,JTYPE) = XDILUDEC(JTYPE) * (ZDECIDUS + 1.1 * ZPHOTON * XPARCF * PSWDIR(JI)       &
                                  + (ZTG_VEG(JI)-XTT)/ZTEMP_VEG - 0.33 * ZFERT(JI))                         &
                                  + (1 - XDILUDEC(JTYPE)) * (1.1 * ZPHOTON * XPARCF * PSWDIR(JI) &
                                  + (ZTG_VEG(JI)-XTT)/ZTEMP_VEG - 0.33 * ZFERT(JI))
            ZPARAM(JI) = ZPARAM(JI) + ZPARAM_TYPE(JI,JTYPE) * IP%XVEGTYPE_PATCH(JI,JTYPE,1)
        ENDDO 
      ENDDO  

      WHERE((IMT%XCE_NITRO(:,1)*IMT%XCNA_NITRO(:,1)+IMT%XCF_NITRO(:,1))/=0.0.AND.IMT%XCNA_NITRO(:,1)/=0.0)
            ZCO2        (:) = P_CO2(:)*(XMD/(1.E-6*XMCO2))  ! (ppmm ->  ppm)
            ZCNA_NITRO  (:) = IMT%XCNA_NITRO(:,1) * &
                    EXP(ZCDILU1*EXP(ZPARAM(:)-IMT%XCNA_NITRO(:,1)/ZCDILU2) * ALOG(MAX(1.,ZCO2(:)/ZCDILU3)))
            ZBSLAI_NITRO(:) = 1. / (IMT%XCE_NITRO(:,1)*ZCNA_NITRO(:)+IMT%XCF_NITRO(:,1))
      ENDWHERE
!
    ENDIF
!    
    IF(ANY(IMT%XLAI(:,1)/=XUNDEF))THEN
      CALL NITRO_DECLINE(IO, IP, IMT, IR, GWOOD, ZBSLAI_NITRO, PLAT, ZBIOMASS_LEAF)
       !print*,'BIOMASS ',ZBIOMASS_LEAF(1)
       !print*,'LAIMIN ',IMT%XLAIMIN(1,1)
       !print*,'BSLAI ',ZBSLAI_NITRO(1)      
      CALL LAIGAIN(ZBSLAI_NITRO, IMT, IR, ZBIOMASS_LEAF)
    ENDIF
!    
  ENDIF
!  
! CASE CPHOTO=AST reinitialise  IR%XANDAY(:,1) and IR%XANFM(:,1) 
  IR%XANDAY(:,1)=0.0
  IR%XANFM(:,1) =0.0
!
ENDIF
!
!
IF (IO%CPHOTO == 'NIT' .OR. IO%CPHOTO=='NCB') THEN
  !
  ! * soil temperature in K (over 1m depth for DIF)
  !
  ZTG_VEG(:) = IR%XTG(:,1,1)
  !
  IF(IO%CISBA/='DIF')THEN        
    ZTG_SOIL(:) = IR%XTG(:,2,1)
  ELSE       
    DO JI=1,INI
       IDEPTH=IMX%NWG_LAYER(JI,1)
       ZDG_SOIL(JI)=MIN(ZDEPTH,IMX%XDG(JI,IDEPTH,1))
    ENDDO  
    DO JL=1,INL
       DO JI=1,INI     
          ZWGHT_SOIL=MIN(IP%XDZG(JI,JL,1),MAX(0.0,ZDG_SOIL(JI)-IMX%XDG(JI,JL,1)+IP%XDZG(JI,JL,1)))        
          ZTG_SOIL(JI)=ZTG_SOIL(JI)+IR%XTG(JI,JL,1)*ZWGHT_SOIL/ZDG_SOIL(JI)
       ENDDO
    ENDDO 
  ENDIF
  !
  !
  ! * Respiration of structural biomass pools
  !
  WHERE(GWOOD(:))
  ! IBIS respiration with either respiration factor rwood=0.0125 - otherwise rroot=1.25 
  ! (Kucharik et al, 2000, eq 6-8) Soil temp in K         
    IR%XRESP_BIOMASS(:,2,1) = IR%XRESP_BIOMASS(:,2,1) + IR%XBIOMASS(:,2,1) * PTSTEP &
                              * MAX(0.,ZROOT_IBIS*EXP(ZCIBIS1*(ZCIBIS2-1./ZTG_VEG(:)))/(ZNDAY*XDAY)) 
  ELSEWHERE 
    IR%XRESP_BIOMASS(:,2,1) = IR%XRESP_BIOMASS(:,2,1) + IR%XBIOMASS(:,2,1) * XRESPFACTOR_NIT    &
                              * EXP((ZLOG2/ZCOEF1)*(ZTG_VEG(:)-XTT-ZCOEF2)) * PTSTEP  
  ! before optimization                   * 2.0**((IR%XTG(:,2,1)-XTT-ZCOEF2)/ZCOEF1) * PTSTEP               
  ENDWHERE
  !
  IF (IO%CPHOTO == 'NIT') THEN
    !
    IR%XRESP_BIOMASS(:,3,1) = IR%XRESP_BIOMASS(:,3,1) + IR%XBIOMASS(:,3,1) * XRESPFACTOR_NIT &
                              * EXP((ZLOG2/ZCOEF1)*(ZTG_SOIL(:)-XTT-ZCOEF2)) * PTSTEP  
    ! before optimization                   * 2.0**((IR%XTG(:,2,1)-XTT-ZCOEF2)/ZCOEF1) * PTSTEP               
    !
  ELSEIF (IO%CPHOTO == 'NCB') THEN
    !
    IR%XRESP_BIOMASS(:,2,1) = MIN(IR%XRESP_BIOMASS(:,2,1), IR%XBIOMASS(:,2,1))
    ! 
    IR%XRESP_BIOMASS(:,3,1) = IR%XRESP_BIOMASS(:,3,1) + IR%XBIOMASS(:,3,1) * &
            MAX( 0., XCOEFF_MAINT_RESP_ZERO * (1. + XSLOPE_MAINT_RESP*(ZTG_VEG(:)-XTT))) * PTSTEP  
    IR%XRESP_BIOMASS(:,3,1) = MIN(IR%XRESP_BIOMASS(:,3,1), IR%XBIOMASS(:,3,1))
    ! 
    WHERE(GWOOD(:))
    ! Resp IBIS (Soil temp in K)
      IR%XRESP_BIOMASS(:,4,1) = IR%XRESP_BIOMASS(:,4,1) + IR%XBIOMASS(:,4,1) * PTSTEP &
                        * MAX(0.,ZROOT_IBIS * EXP(ZCIBIS1*(ZCIBIS2-1./ZTG_SOIL(:)))/(ZNDAY*XDAY))
    ELSEWHERE 
    IR%XRESP_BIOMASS(:,4,1) = IR%XRESP_BIOMASS(:,4,1) + IR%XBIOMASS(:,4,1) * &
             MAX( 0., XCOEFF_MAINT_RESP_ZERO * (1. + XSLOPE_MAINT_RESP*(ZTG_SOIL(:)-XTT))) * PTSTEP  
    ENDWHERE
    !
    IR%XRESP_BIOMASS(:,4,1) = MIN(IR%XRESP_BIOMASS(:,4,1), IR%XBIOMASS(:,4,1))
    !
    WHERE( (GWOOD(:)).AND.(IR%XBIOMASS(:,5,1)>0.) )
    ! IBIS estimation of sapwood fraction based on the height of tree, sapspeed and 
    ! max transpiration rates. Conversion from DM to C. To be changed with DGVM.  (Soil temp in K)        
      ZHTREE(:) = 2.5*0.75*(IR%XBIOMASS(:,1,1)+IR%XBIOMASS(:,2,1)+IR%XBIOMASS(:,3,1)+&
                            IR%XBIOMASS(:,4,1)+IR%XBIOMASS(:,5,1)+IR%XBIOMASS(:,6,1))*0.4
      ZSAPFRAC(:) = MIN(0.5, MAX(0.05,0.0025/25.*ZHTREE(:)*0.75*400/(IR%XBIOMASS(:,5,1)*0.4)))
      !ZSAPFRAC(:) = 0.5
      
      IR%XRESP_BIOMASS(:,5,1) = IR%XRESP_BIOMASS(:,5,1) + IR%XBIOMASS(:,5,1) * ZSAPFRAC(:) * PTSTEP &
                                 * MAX(0.,ZWOOD_IBIS*EXP(ZCIBIS1*(ZCIBIS2-1./ZTG_VEG(:)))/(ZNDAY*XDAY))
      IR%XRESP_BIOMASS(:,5,1) = MIN(IR%XRESP_BIOMASS(:,5,1), IR%XBIOMASS(:,5,1))
    ELSEWHERE
      IR%XRESP_BIOMASS(:,5,1) = 0.0
    ENDWHERE

    !
  ENDIF
  !
  ! * Instantaneous respiration (kgCO2/kgair m/s)
  !
  DO JL=2,SIZE(IR%XRESP_BIOMASS(:,:,1),2)
      PRESP_BIOMASS_INST(:,JL) = (IR%XRESP_BIOMASS(:,JL,1) - ZRESP_BIOMASS_LAST(:,JL)) &
                                     * XPCCO2*XMCO2/(PTSTEP*PRHOA(:)*XMC)                              
  ENDDO
 !  
ENDIF

!*      3.     Agricultural practices
!              ----------------------
!
IF (OAGRIP) THEN
  !
  GMASK_AGRI(:) = .FALSE.
  WHERE ( IMI%TSEED(:,1)%TDATE%MONTH /= NUNDEF .AND. ( KMONTH < IMI%TSEED(:,1)%TDATE%MONTH .OR. &
         (KMONTH == IMI%TSEED(:,1)%TDATE%MONTH .AND. KDAY < IMI%TSEED(:,1)%TDATE%DAY) ) )  GMASK_AGRI(:) = .TRUE.
  WHERE ( IMI%TREAP(:,1)%TDATE%MONTH /= NUNDEF .AND. ( KMONTH > IMI%TREAP(:,1)%TDATE%MONTH .OR. &
         (KMONTH == IMI%TREAP(:,1)%TDATE%MONTH .AND. KDAY >= IMI%TREAP(:,1)%TDATE%DAY) ) ) GMASK_AGRI(:) = .TRUE. 
  !
  WHERE (GMASK_AGRI(:))
    IMT%XLAI(:,1)             = IMT%XLAIMIN(:,1)
    ZBIOMASS_LEAF(:)    = IMT%XLAI(:,1) * ZBSLAI_NITRO(:)
  END WHERE

  IF (IO%CPHOTO == 'NIT' .OR. IO%CPHOTO == 'NCB') THEN
    !
    WHERE (GMASK_AGRI(:))
      IR%XBIOMASS(:,1,1)       = 0.0
      IR%XBIOMASS(:,2,1)       = 0.0
      IR%XBIOMASS(:,3,1)       = 0.0
      IR%XRESP_BIOMASS(:,2,1)  = 0.0
      IR%XRESP_BIOMASS(:,3,1)  = 0.0
    END WHERE
    !
    IF (IO%CPHOTO == 'NCB') THEN
      !
      WHERE (GMASK_AGRI(:)) 
        IR%XBIOMASS(:,4,1)       = 0.0
        IR%XBIOMASS(:,5,1)       = 0.0
        IR%XBIOMASS(:,6,1)       = 0.0
        IR%XRESP_BIOMASS(:,4,1)  = 0.0
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
  !if (size(IMT%XLAI,1)>31) THEN
  !print*,'LAI ',IMT%XLAI(32,1)
  !print*,'VEGTYPE ',IMX%XVEGTYPE(32,:)
  !print*,'AGRI_TO_GRASS ',IO%LAGRI_TO_GRASS
  !endif
  WHERE( IMT%XVEG(:,1) > 0. )
    ! Evolution of vegetation fraction and roughness length due to LAI change
    IMT%XZ0(:,1) = Z0V_FROM_LAI(IMT%XLAI(:,1),IMX%XH_TREE(:,1),IP%XVEGTYPE_PATCH(:,:,1),IO%LAGRI_TO_GRASS) 
    IMT%XVEG(:,1) = VEG_FROM_LAI(IMT%XLAI(:,1),IP%XVEGTYPE_PATCH(:,:,1),IO%LAGRI_TO_GRASS)
    !
    ! Evolution of radiative parameters due to vegetation fraction change
    IMT%XEMIS(:,1)= EMIS_FROM_VEG(IMT%XVEG(:,1),IP%XVEGTYPE_PATCH(:,:,1))
  END WHERE
  !
  CALL ALBEDO(IO%CALBEDO, IMT, IMA )  
  !
  ! Evolution of effective roughness length due to new surface roughness length
  !
  IF (ASSOCIATED(MSS%XAOSIP)) THEN
    IF (SIZE(MSS%XAOSIP)>0) THEN
      CALL SUBSCALE_Z0EFF(MSS,SPREAD(IMT%XZ0(:,1),2,SIZE(MSS%XZ0EFFIP,2)),.FALSE. )
    ENDIF
  ENDIF
  !
ENDIF
!
!print*,'Z0EFFIP ',MSS%XZ0EFFIP(1,:)
!print*,'Z0EFFIM ',MSS%XZ0EFFIM(1,:)
!print*,'Z0EFFJP ',MSS%XZ0EFFJP(1,:)
!print*,'Z0EFFJM ',MSS%XZ0EFFJM(1,:)
!
!print*,'LAI ',IMT%XLAI(1,1)
!print*,'VEG ',IMT%XVEG(1,1)
!print*,'Z0 ',IMT%XZ0(1,1)
!print*,'ALBNIR ',IMT%XALBNIR(1,1)
!print*,'ALBVIS ',IMT%XALBVIS(1,1)
!print*,'ALBUV ',IMT%XALBUV(1,1)
!print*,'EMIS ',IMT%XEMIS(1,1)
!
!print*,'ANFM ',IR%XANFM(1,1)
!print*,'ANDAY ',IR%XANDAY(1,1)
!print*,'BIOMASS ',IR%XBIOMASS(1,:,1)
!print*,'RESP_BIOMASS ',IR%XRESP_BIOMASS(1,:,1)
!print*,'RESP_BIOMASS_INST ',PRESP_BIOMASS_INST(1,:)
!
!print*,'INCREASE ',IP%XINCREASE(1,:,1)
!print*,'TURNOVER ',IP%XTURNOVER(1,:,1)
!
IF (LHOOK) CALL DR_HOOK('VEGETATION_EVOL',1,ZHOOK_HANDLE)
!-----------------------------------------------------------------
!
END SUBROUTINE VEGETATION_EVOL
