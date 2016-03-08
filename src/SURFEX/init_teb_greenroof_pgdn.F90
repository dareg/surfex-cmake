!#############################################################
SUBROUTINE INIT_TEB_GREENROOF_PGD_n (DTCO, U, OCH_BIO_FLUX, TG, PGREENROOF, TOP, &
                                     GDO, VO, VMX, VMT, VP, VIP, DT, GB, &
                                     HPROGRAM,HINIT,OREAD_PGD, KI, KVERSION, PCO2, PRHOA)
!#############################################################
!
!!****  *INIT_TEB_GREENROOF_PGD_n* - routine to initialize ISBA
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      A. Lemonsu  *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    09/2009
!!                  11/2013 (B. Decharme) No exp profile with DIF
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SSO_n, ONLY : SSO_t, SSO_INIT
USE MODD_GRID_n, ONLY : GRID_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_FIX_t, ISBA_PARAM_TIME_t
USE MODD_ISBA_PGD_n, ONLY : ISBA_PGD_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t
!
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_GR_BIOG_n, ONLY : GR_BIOG_t
!
USE MODD_AGRI_n, ONLY : AGRI_t
!
USE MODD_TYPE_DATE_SURF
USE MODD_TYPE_SNOW
!


USE MODD_DATA_COVER_PAR,       ONLY: NVEGTYPE
USE MODD_SURF_PAR,             ONLY: XUNDEF, NUNDEF

USE MODD_SGH_PAR,              ONLY: XF_DECAY
!
USE MODI_GET_LUOUT
USE MODI_ALLOCATE_TEB_GREENROOF_PGD
USE MODI_READ_PGD_TEB_GREENROOF_n
USE MODI_CONVERT_PATCH_ISBA
USE MODI_INIT_FROM_DATA_GREENROOF_n
USE MODI_INIT_VEG_PGD_n
USE MODI_EXP_DECAY_SOIL_FR
USE MODI_ABOR1_SFX
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
LOGICAL, INTENT(IN) :: OCH_BIO_FLUX
TYPE(GRID_t), INTENT(INOUT) :: TG
REAL, DIMENSION(:), INTENT(IN) :: PGREENROOF
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: GDO
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: VO
TYPE(ISBA_PARAM_FIX_t), INTENT(INOUT) :: VMX
TYPE(ISBA_PARAM_TIME_t), INTENT(INOUT) :: VMT
TYPE(ISBA_PGD_t), INTENT(INOUT) :: VP
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: VIP
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DT
TYPE(GR_BIOG_t), INTENT(INOUT) :: GB
!
 CHARACTER(LEN=6),                   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=3),                   INTENT(IN)  :: HINIT     ! choice of fields to initialize
LOGICAL,                            INTENT(IN)  :: OREAD_PGD ! flag to read PGD fields in the file
INTEGER,                            INTENT(IN)  :: KI        ! number of points
INTEGER,                            INTENT(IN)  :: KVERSION  ! version number of the file being read
REAL,             DIMENSION(KI),    INTENT(IN)  :: PCO2        ! CO2 concentration (kg/m3)
REAL,             DIMENSION(KI),    INTENT(IN)  :: PRHOA       ! air density
!
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
TYPE(SSO_t) :: YGRSS
TYPE(AGRI_t) :: YAG
!
INTEGER           :: JILU     ! loop increment
INTEGER           :: ILUOUT   ! unit of output listing file
!
INTEGER           :: IDECADE  ! decade of simulation
!
INTEGER :: JVEGTYPE, JLAYER  ! loop counter on layers
!
REAL, DIMENSION(KI)               :: ZF
REAL, DIMENSION(KI)               :: ZWORK
!
!*       0.3   Soil parameter values for organic matter - from Lawrence and Slater (2008):
!              ----------------------------------------------------------------------------------
!
REAL, PARAMETER   :: ZWSAT_OM      = 0.9       ! Porosity of OM (m3/m3)
REAL, PARAMETER   :: ZCONDSAT_OM   = 2.8E-4    ! Saturated hydraulic conductivity for OM (m/s)
REAL, PARAMETER   :: ZMPOTSAT_OM   = -10.3E-3  ! Saturated matric potential for OM (m)
REAL, PARAMETER   :: ZBCOEF_OM     = 2.7       ! CH78 b-parameter for OM (-)
!
REAL, PARAMETER   :: ZCONDDRY_OM   = 0.05      ! Dry thermal conductivity for OM (W/m/K)
REAL, PARAMETER   :: ZCONDSLD_OM   = 0.25      ! Soil solids thermal conductivity for OM (W/m/K)
REAL, PARAMETER   :: ZHCAPSOIL_OM  = 2.5E+6    ! Soil heat capacity for OM
!
REAL, PARAMETER   :: ZMPOT_WWILT   = -150.     ! Matric potential at wilting point (m)
REAL, PARAMETER   :: ZHYDCOND_WFC  = 1.157E-9  ! Hydraulic conductivity at field capacity (m/s)
!
REAL, DIMENSION(0) :: ZTDEEP_CLI, ZGAMMAT_CLI, ZTHRESHOLD
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!               Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('INIT_TEB_GREENROOF_PGD_n',0,ZHOOK_HANDLE)
!
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
CALL SSO_INIT(YGRSS)
!-------------------------------------------------------------------------------
!
!*       2.     Physiographic fields
!               --------------------
!
VO%NPATCH = 1
VO%NNBIOMASS = GDO%NNBIOMASS
VO%CPHOTO = GDO%CPHOTO
VO%CPEDOTF = GDO%CPEDOTF
VO%CRESPSL = GDO%CRESPSL
!
!*       2.1    Cover, soil and orographic fields:
!               ---------------------------------
!
IF (OREAD_PGD) &
 CALL READ_PGD_TEB_GREENROOF_n(OCH_BIO_FLUX, DTCO, DT, GB, U, &
                               VO, VP, TG%NDIM, HPROGRAM,KVERSION)
!
!
!* allocation of green roofs variables
!
 CALL ALLOCATE_TEB_GREENROOF_PGD(VMX, VMT, VP, VIP, &
                                 OREAD_PGD, KI, NVEGTYPE, VO%NGROUND_LAYER)
!
!*       2.2    Physiographic data fields from land cover:
!               -----------------------------------------
!
IF (TOP%TTIME%TDATE%MONTH /= NUNDEF) THEN
  IDECADE = 3 * ( TOP%TTIME%TDATE%MONTH - 1 ) + MIN(TOP%TTIME%TDATE%DAY-1,29) / 10 + 1
ELSE
  IDECADE = 1
END IF
!
!
IF (.NOT. VO%LPAR) THEN
  CALL CONVERT_PATCH_ISBA(DTCO, DT, VO, IDECADE, IDECADE, TOP%XCOVER, TOP%LCOVER,&
                       .FALSE., 'GRD', .FALSE., VMX, VMT, PSOILGRID=VO%XSOILGRID    )   
ELSE
  CALL INIT_FROM_DATA_GREENROOF_n(DT, IDECADE,VO%CPHOTO, .FALSE., VMX=VMX, VMT=VMT  ) 
  ! 
  IF (VO%CISBA=='DIF') THEN
    WHERE(PGREENROOF(:)/=0.)
      VMX%NWG_LAYER(:,1)=VO%NGROUND_LAYER
      VMX%XDG2  (:,1)=0.0
      VMX%XDROOT(:,1)=0.0
    ENDWHERE
    DO JLAYER=VO%NGROUND_LAYER,1,-1
      DO JILU=1,KI
        IF(PGREENROOF(JILU)/=0..AND.VMX%XROOTFRAC(JILU,JLAYER,1)>=1.0)THEN
          VMX%XDG2  (JILU,1)=VMX%XDG(JILU,JLAYER,1)
          VMX%XDROOT(JILU,1)=VMX%XDG(JILU,JLAYER,1)
        ENDIF
      ENDDO
    ENDDO
  ENDIF
END IF
!
WHERE (PGREENROOF(:)==0.)
  ! GARDEN default values /may need changing for green roofs
  VP%XSOC    (:,1) = 0.5
  VP%XSOC    (:,2) = 0.5
  VP%XSAND   (:,1) = 0.33
  VP%XSAND   (:,2) = 0.33
  VP%XCLAY   (:,1) = 0.33
  VP%XCLAY   (:,2) = 0.33
  VMT%XVEG       (:  ,1) = 0.
  VMT%XLAI       (:  ,1) = 0.
  VMT%XRSMIN     (:  ,1) = 40.
  VMT%XGAMMA     (:  ,1) = 0.
  VMT%XWRMAX_CF  (:  ,1) = 0.2
  VMT%XRGL       (:  ,1) = 100.
  VMT%XCV        (:  ,1) = 2.E-5
  VMT%XZ0        (: ,1 ) = 0.013
  VMX%XZ0_O_Z0H  (: ,1 ) = 10.
  VMT%XALBNIR_VEG(: ,1 ) = 0.30
  VMT%XALBVIS_VEG(:  ,1) = 0.30
  VMT%XALBUV_VEG (:  ,1) = 0.06
  VMT%XEMIS      (: ,1 ) = 0.94
END WHERE
IF (VO%CPHOTO/='NON') THEN
  WHERE (PGREENROOF(:)==0.)
    VMT%XGMES      (:  ,1) = 0.020
    VMT%XBSLAI     (:  ,1) = 0.36
    VMT%XLAIMIN    (:  ,1) = 0.3
    VMT%XSEFOLD    (:  ,1) = 90*86400.
    VMX%XH_TREE    (: ,1 ) = 0.
    VMX%XRE25      (: ,1 ) = 3.6E-7    
    VMT%XGC        (: ,1 ) = 0.00025
  END WHERE
  IF (VO%CPHOTO/='AGS' .AND. VO%CPHOTO/='LAI') THEN
    WHERE (PGREENROOF(:)==0.)     
      VMX%XDMAX      (:  ,1) = 0.1
      VMT%XF2I       (: ,1 ) = 0.3
    END WHERE
    IF (VO%CPHOTO=='NIT' .OR. VO%CPHOTO=='NCB') THEN
      WHERE (PGREENROOF(:)==0.)          
        VMT%XCE_NITRO  (: ,1 ) = 7.68
        VMT%XCF_NITRO  (: ,1 ) = -4.33
        VMT%XCNA_NITRO (: ,1 ) = 1.3
      END WHERE
    ENDIF
  ENDIF
ENDIF  
IF(VO%CISBA/='DIF')THEN
  DO JLAYER=1,VO%NGROUND_LAYER
    WHERE (PGREENROOF(:)==0.)
      VMX%XDG(:,JLAYER,1)=0.2*JLAYER
    END WHERE
  ENDDO
ELSE
  WHERE (PGREENROOF(:)==0.) 
    VMX%XDG(:,1,1)=0.01
    VMX%XDG(:,2,1)=0.04
    VMX%XROOTFRAC(:,1,1)=0.
    VMX%XROOTFRAC(:,2,1)=0.
  END WHERE        
  DO JLAYER=3,VO%NGROUND_LAYER
    WHERE (PGREENROOF(:)==0.)
      VMX%XDG(:,JLAYER,1)=0.1*(JLAYER-2)
      VMX%XROOTFRAC(:,JLAYER,1)=0.
    END WHERE
  ENDDO               
  WHERE (PGREENROOF(:)==0.) 
    VMX%NWG_LAYER(:,1)=VO%NGROUND_LAYER
    VMX%XDROOT   (:,1)=0.0
    VMX%XDG2     (:,1)=VMX%XDG(:,VO%NGROUND_LAYER-1,1)
  ENDWHERE    
ENDIF  
WHERE (PGREENROOF(:)==0.) 
  VMX%XD_ICE(:,1)=0.8*VMX%XDG(:,2,1)
END WHERE  
DO JVEGTYPE=1,NVEGTYPE
  WHERE (PGREENROOF(:)==0.)
    VMX%XVEGTYPE(:,JVEGTYPE)=0.
    VMX%XVEGTYPE(:,1)=1.
  END WHERE
ENDDO
!
ALLOCATE(YGRSS%XAOSIP(0))
!
  CALL INIT_VEG_PGD_n(YGRSS, DTCO, U, VO, VIP, VP, VMX, VMT, YAG, &
                      HPROGRAM, 'TOWN  ',ILUOUT, KI, TOP%TTIME%TDATE%MONTH,    &
                      .FALSE., .FALSE., ZTDEEP_CLI, ZGAMMAT_CLI, &
                      .FALSE., ZTHRESHOLD, HINIT, PCO2, PRHOA  )
!
!-------------------------------------------------------------------------------
!
!*       5.1     Soil thermal characteristics for greenroofs:
!               ----------------------------------------------
!
! WARNING: must be done before soil hydraulic characteristics (because of WSAT)
! Estimation of WSAT_MI for use in HEATCAPZ and THRMCONDZ for mineral fraction
! and allow weighted combination with regard to OM & no-OM fractions:
!
IF (VO%CSCOND=='PL98' .OR. VO%CISBA=='DIF') THEN
  DO JLAYER=1,VO%NGROUND_LAYER
    VIP%XHCAPSOIL(:,JLAYER) = VP%XSOC(:,JLAYER) * ZHCAPSOIL_OM + &
                        (1-VP%XSOC(:,JLAYER)) * VIP%XHCAPSOIL(:,JLAYER)  
    VIP%XCONDDRY(:,JLAYER) = (ZCONDDRY_OM         * VIP%XCONDDRY(:,JLAYER))    &
                         /(  VP%XSOC(:,JLAYER)  * VIP%XCONDDRY(:,JLAYER) +   &
                          (1-VP%XSOC(:,JLAYER)) * ZCONDDRY_OM)
    VIP%XCONDSLD(:,JLAYER) = (ZCONDSLD_OM         * VIP%XCONDSLD(:,JLAYER))    &
                         /(  VP%XSOC(:,JLAYER)  * VIP%XCONDSLD(:,JLAYER) +   &
                          (1-VP%XSOC(:,JLAYER)) * ZCONDSLD_OM)
  ENDDO
END IF
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Validation case : experimental values for Nancy 2011 case
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Substrate layer
DO JLAYER=1,4
    VIP%XCONDDRY (:,JLAYER) = 0.15
    VIP%XHCAPSOIL(:,JLAYER) = 1342000.
ENDDO
! Drainage layer
DO JLAYER=5,6
    VIP%XCONDDRY (:,JLAYER) = 0.09
    VIP%XHCAPSOIL(:,JLAYER) = 331500.
ENDDO
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!*       5.2     Soil thermal characteristics:
!               --------------------------------
!
DO JLAYER=1,VO%NGROUND_LAYER
  VIP%XCONDSAT(:,JLAYER,1) =   VP%XSOC(:,JLAYER)* ZCONDSAT_OM   &
                        +(1-VP%XSOC(:,JLAYER))* VIP%XCONDSAT(:,JLAYER,1)
END DO
!
! Note that if ISBA/=DIF, always CDIF = 'BC' and CPEDOTF = 'CH78'
DO JLAYER=1,VO%NGROUND_LAYER
  VIP%XBCOEF  (:,JLAYER) =    VP%XSOC(:,JLAYER) * ZBCOEF_OM        &
                       +(1-VP%XSOC(:,JLAYER))* VIP%XBCOEF(:,JLAYER)
  VIP%XMPOTSAT(:,JLAYER) =    VP%XSOC(:,JLAYER) * ZMPOTSAT_OM      &
                       +(1-VP%XSOC(:,JLAYER))* VIP%XMPOTSAT(:,JLAYER)
END DO
!        
DO JLAYER=1,VO%NGROUND_LAYER
   VIP%XWSAT (:,JLAYER) =    VP%XSOC(:,JLAYER)* ZWSAT_OM            &
                     +(1-VP%XSOC(:,JLAYER))* VIP%XWSAT(:,JLAYER)
   VIP%XWWILT(:,JLAYER) = EXP(((LOG(-1*ZMPOT_WWILT)-LOG(-1*VIP%XMPOTSAT(:,JLAYER)))   &
                    / (-1*VIP%XBCOEF(:,JLAYER)))+LOG(VIP%XWSAT(:,JLAYER)))
   VIP%XWFC  (:,JLAYER) = EXP(((LOG(ZHYDCOND_WFC)-LOG(VIP%XCONDSAT(:,JLAYER,1)))        &
                    / (2*VIP%XBCOEF(:,JLAYER)+3))+LOG(VIP%XWSAT(:,JLAYER)))
END DO
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Validation case : experimental values for Nancy 2011 case
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Substrate layer
DO JLAYER=1,4
  VIP%XWSAT   (:,JLAYER) = 0.674     ! Value tested
  VIP%XCONDSAT(:,JLAYER,1) = 2.162E-3  ! Value tested
  VIP%XMPOTSAT(:,JLAYER) = -0.932    ! Value tested
  VIP%XBCOEF  (:,JLAYER) = 3.9       ! Value tested
  VIP%XWWILT  (:,JLAYER) = 0.15      ! from OBS-NANCY
  VIP%XWFC    (:,JLAYER) = 0.37      ! from OBS-NANCY
ENDDO
! Drainage layer
DO JLAYER=5,6
   VIP%XWSAT   (:,JLAYER) = 0.9       ! Value tested
   VIP%XCONDSAT(:,JLAYER,1) = 3.32E-3   ! Value tested
   VIP%XMPOTSAT(:,JLAYER) = -0.121    ! Value tested
   VIP%XBCOEF  (:,JLAYER) = 2.7       ! Value tested
   VIP%XWWILT  (:,JLAYER) = 0.15      ! sert à initialiser le WG ds la couche
   VIP%XWFC    (:,JLAYER) = 0.37      ! sert à initialiser le WG ds la couche
ENDDO
!-------------------------------------------------------------------------------
!
!*       6.1    Initialize of the SGH scheme:'
!               ------------------------------
!
IF(VO%CKSAT=='SGH' .AND. VO%CISBA/='DIF' .AND. HINIT/='PRE')THEN 
  ZF(:)=MIN(4.0/VMX%XDG(:,2,1),XF_DECAY)
  CALL EXP_DECAY_SOIL_FR(VO%CISBA, SPREAD(ZF,2,1),VIP, VMX)
ENDIF
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('INIT_TEB_GREENROOF_PGD_n',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE INIT_TEB_GREENROOF_PGD_n
