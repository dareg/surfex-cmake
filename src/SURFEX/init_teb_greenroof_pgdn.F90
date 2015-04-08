!#############################################################
SUBROUTINE INIT_TEB_GREENROOF_PGD_n(HPROGRAM,HINIT,OREAD_PGD, KI, KSV, HSV, KVERSION, PCO2, PRHOA)
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
USE MODD_TYPE_DATE_SURF
USE MODD_TYPE_SNOW
!
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS
USE MODD_TEB_n, ONLY : T => TEB
USE MODD_TEB_VEG_n, ONLY : TVG => TEB_VEG_OPTIONS
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TGRO => TEB_GREENROOF_OPTIONS

USE MODD_TEB_GREENROOF_PGD_EVOL_n, ONLY : TGRPE => TEB_GREENROOF_PGD_EVOL
USE MODD_TEB_GREENROOF_PGD_n, ONLY : TGRP => TEB_GREENROOF_PGD
USE MODD_CH_TEB_n, ONLY : CHT => CH_TEB

USE MODD_DATA_COVER_PAR,       ONLY: NVEGTYPE
USE MODD_SURF_PAR,             ONLY: XUNDEF, NUNDEF

USE MODD_SGH_PAR,              ONLY: NDIMTAB, XF_DECAY
!
USE MODI_GET_LUOUT
USE MODI_ALLOCATE_TEB_GREENROOF_PGD
USE MODI_READ_PGD_TEB_GREENROOF_n
USE MODI_CONVERT_PATCH_TEB_GREENROOF
USE MODI_INIT_FROM_DATA_GREENROOF_n
USE MODI_INIT_VEG_PGD_GARDEN_n
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
 CHARACTER(LEN=6),                   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=3),                   INTENT(IN)  :: HINIT     ! choice of fields to initialize
LOGICAL,                            INTENT(IN)  :: OREAD_PGD ! flag to read PGD fields in the file
INTEGER,                            INTENT(IN)  :: KI        ! number of points
INTEGER,                            INTENT(IN)  :: KSV       ! number of scalars
 CHARACTER(LEN=6), DIMENSION(KSV),   INTENT(IN)  :: HSV       ! name of all scalar variables
INTEGER,                            INTENT(IN)  :: KVERSION  ! version number of the file being read
REAL,             DIMENSION(KI),    INTENT(IN)  :: PCO2        ! CO2 concentration (kg/m3)
REAL,             DIMENSION(KI),    INTENT(IN)  :: PRHOA       ! air density
!
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
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
!                                              ! = 0.1 mm/day
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!               Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('INIT_TEB_GREENROOF_PGD_n',0,ZHOOK_HANDLE)
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
!-------------------------------------------------------------------------------
!
!*       2.     Physiographic fields
!               --------------------
!
!
!*       2.1    Cover, soil and orographic fields:
!               ---------------------------------
!
IF (OREAD_PGD) &
 CALL READ_PGD_TEB_GREENROOF_n(HPROGRAM,KVERSION)
!
!
!* allocation of green roofs variables
!
 CALL ALLOCATE_TEB_GREENROOF_PGD(OREAD_PGD, KI, NVEGTYPE, TGRO%NLAYER_GR, NDIMTAB)
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
IF (.NOT. TGRO%LPAR_GREENROOF) THEN
  CALL CONVERT_PATCH_TEB_GREENROOF(KI,IDECADE)
ELSE
 CALL INIT_FROM_DATA_GREENROOF_n(IDECADE,TVG%CPHOTO,                     &
                                 TGRP%XOM_GR,                             &
                                 TGRP%XSAND_GR, TGRP%XCLAY_GR, TGRPE%XVEG,           &
                                 TGRPE%XLAI,TGRP%XRSMIN,TGRP%XGAMMA,TGRP%XWRMAX_CF,       &
                                 TGRP%XRGL,TGRP%XCV,TGRP%XDG,TGRP%XD_ICE,TGRPE%XZ0,TGRP%XZ0_O_Z0H,  &
                                 TGRP%XALBNIR_VEG,TGRP%XALBVIS_VEG,            &
                                 TGRP%XALBUV_VEG,TGRPE%XEMIS,                   &
                                 TGRP%XVEGTYPE,TGRP%XROOTFRAC,                 &
                                 TGRP%XGMES,TGRP%XBSLAI,TGRP%XLAIMIN,TGRP%XSEFOLD,TGRP%XGC,   &
                                 TGRP%XDMAX, TGRP%XF2I, TGRP%LSTRESS, TGRP%XH_TREE,TGRP%XRE25,&
                                 TGRP%XCE_NITRO,TGRP%XCF_NITRO,TGRP%XCNA_NITRO      )  
  IF (TGRO%CISBA_GR=='DIF') THEN
    WHERE(T%XGREENROOF(:)/=0.)
      TGRP%NWG_LAYER(:)=TGRO%NLAYER_GR 
      TGRP%XDG2  (:)=0.0
      TGRP%XDROOT(:)=0.0
    ENDWHERE
    DO JLAYER=TGRO%NLAYER_GR,1,-1
      DO JILU=1,KI
        IF(T%XGREENROOF(JILU)/=0..AND.TGRP%XROOTFRAC(JILU,JLAYER)>=1.0)THEN
          TGRP%XDG2  (JILU)=TGRP%XDG(JILU,JLAYER)
          TGRP%XDROOT(JILU)=TGRP%XDG(JILU,JLAYER)
        ENDIF
      ENDDO
    ENDDO
  ENDIF
END IF
!
WHERE (T%XGREENROOF(:)==0.)
  ! GARDEN default values /may need changing for green roofs
  TGRP%XOM_GR     (:,1) = 0.5
  TGRP%XOM_GR     (:,2) = 0.5
  TGRP%XSAND_GR   (:,1) = 0.33
  TGRP%XSAND_GR   (:,2) = 0.33
  TGRP%XCLAY_GR   (:,1) = 0.33
  TGRP%XCLAY_GR   (:,2) = 0.33
  TGRPE%XVEG       (:  ) = 0.
  TGRPE%XLAI       (:  ) = 0.
  TGRP%XRSMIN     (:  ) = 40.
  TGRP%XGAMMA     (:  ) = 0.
  TGRP%XWRMAX_CF  (:  ) = 0.2
  TGRP%XRGL       (:  ) = 100.
  TGRP%XCV        (:  ) = 2.E-5
  TGRPE%XZ0        (:  ) = 0.013
  TGRP%XZ0_O_Z0H  (:  ) = 10.
  TGRP%XALBNIR_VEG(:  ) = 0.30
  TGRP%XALBVIS_VEG(:  ) = 0.30
  TGRP%XALBUV_VEG (:  ) = 0.06
  TGRPE%XEMIS      (:  ) = 0.94
END WHERE
IF (TVG%CPHOTO/='NON') THEN
  WHERE (T%XGREENROOF(:)==0.)
    TGRP%XGMES      (:  ) = 0.020
    TGRP%XBSLAI     (:  ) = 0.36
    TGRP%XLAIMIN    (:  ) = 0.3
    TGRP%XSEFOLD    (:  ) = 90*86400.
    TGRP%XH_TREE    (:  ) = 0.
    TGRP%XRE25      (:  ) = 3.6E-7    
    TGRP%XGC        (:  ) = 0.00025
  END WHERE
  IF (TVG%CPHOTO/='AGS' .AND. TVG%CPHOTO/='LAI') THEN
    WHERE (T%XGREENROOF(:)==0.)     
      TGRP%XDMAX      (:  ) = 0.1
      TGRP%XF2I       (:  ) = 0.3
    END WHERE
    IF (TVG%CPHOTO=='NIT' .OR. TVG%CPHOTO=='NCB') THEN
      WHERE (T%XGREENROOF(:)==0.)          
        TGRP%XCE_NITRO  (:  ) = 7.68
        TGRP%XCF_NITRO  (:  ) = -4.33
        TGRP%XCNA_NITRO (:  ) = 1.3
      END WHERE
    ENDIF
  ENDIF
ENDIF  
IF(TGRO%CISBA_GR/='DIF')THEN
  DO JLAYER=1,TGRO%NLAYER_GR
    WHERE (T%XGREENROOF(:)==0.)
      TGRP%XDG(:,JLAYER)=0.2*JLAYER
    END WHERE
  ENDDO
ELSE
  WHERE (T%XGREENROOF(:)==0.) 
    TGRP%XDG(:,1)=0.01
    TGRP%XDG(:,2)=0.04
    TGRP%XROOTFRAC(:,1)=0.
    TGRP%XROOTFRAC(:,2)=0.
  END WHERE        
  DO JLAYER=3,TGRO%NLAYER_GR
    WHERE (T%XGREENROOF(:)==0.)
      TGRP%XDG(:,JLAYER)=0.1*(JLAYER-2)
      TGRP%XROOTFRAC(:,JLAYER)=0.
    END WHERE
  ENDDO               
  WHERE (T%XGREENROOF(:)==0.) 
    TGRP%NWG_LAYER(:)=TGRO%NLAYER_GR
    TGRP%XDROOT   (:)=0.0
    TGRP%XDG2     (:)=TGRP%XDG(:,TGRO%NLAYER_GR-1)
  ENDWHERE    
ENDIF  
WHERE (T%XGREENROOF(:)==0.) 
  TGRP%XD_ICE(:)=0.8*TGRP%XDG(:,2)
END WHERE  
DO JVEGTYPE=1,NVEGTYPE
  WHERE (T%XGREENROOF(:)==0.)
    TGRP%XVEGTYPE(:,JVEGTYPE)=0.
    TGRP%XVEGTYPE(:,1)=1.
  END WHERE
ENDDO
!
 CALL INIT_VEG_PGD_GARDEN_n(HPROGRAM, ILUOUT, KI, TGRO%NLAYER_GR, TOP%TTIME%TDATE%MONTH,    &
                        TGRP%XVEGTYPE, TGRP%XTDEEP, TGRP%XGAMMAT, TVG%CPHOTO, HINIT, TGRO%LTR_ML_GR,        &
                        TGRO%CRUNOFF_GR,                                                 &
                        TVG%NNBIOMASS, PCO2, PRHOA, TGRP%XABC, TGRP%XPOI,                         &
                        TGRP%XGMES, TGRP%XGC, TGRP%XDMAX, TGRP%XANMAX, TGRP%XFZERO, TGRP%XEPSO, TGRP%XGAMM, TGRP%XQDGAMM,   &
                        TGRP%XQDGMES, TGRP%XT1GMES, TGRP%XT2GMES, TGRP%XAMAX, TGRP%XQDAMAX, TGRP%XT1AMAX, TGRP%XT2AMAX,&
                        TGRP%XAH, TGRP%XBH,                                                   &
                        KSV, HSV, CHT%NBEQ, CHT%CSV, CHT%NAEREQ, CHT%NSV_CHSBEG, CHT%NSV_CHSEND,        &
                        CHT%NSV_AERBEG, CHT%NSV_AEREND, CHT%CCH_NAMES, CHT%CAER_NAMES, CHT%NDSTEQ,      &
                        CHT%NSV_DSTBEG, CHT%NSV_DSTEND, CHT%NSLTEQ, CHT%NSV_SLTBEG, CHT%NSV_SLTEND,     &
                        CHT%CDSTNAMES, CHT%CSLTNAMES, CHT%CCHEM_SURF_FILE,                      &
                        TGRP%XCLAY_GR, TGRP%XSAND_GR, TVG%CPEDOTF,                                &
                        TGRP%XCONDSAT, TGRP%XMPOTSAT, TGRP%XBCOEF, TGRP%XWWILT, TGRP%XWFC, TGRP%XWSAT,            &
                        TGRP%XTAUICE, TGRP%XCGSAT, TGRP%XC1SAT, TGRP%XC2REF, TGRP%XC3, TGRP%XC4B, TGRP%XACOEF, TGRP%XPCOEF, &
                        TGRP%XC4REF, TGRP%XPCPS, TGRP%XPLVTT, TGRP%XPLSTT,                              &
                        TGRO%CSCOND_GR, TGRO%CISBA_GR, TGRP%XHCAPSOIL, TGRP%XCONDDRY, TGRP%XCONDSLD, TVG%CCPSURF,&
                        TGRP%XDG, TGRP%XDROOT, TGRP%XDG2, TGRP%XROOTFRAC, TGRP%XRUNOFFD, TGRP%XDZG, TGRP%XDZDIF,       &
                        TGRP%XSOILWGHT, TGRP%NWG_LAYER, TGRO%NLAYER_HORT_GR, TGRO%NLAYER_DUN_GR, TGRP%XD_ICE,&
                        TGRP%XKSAT_ICE, TGRP%XALBNIR_DRY, TGRP%XALBVIS_DRY, TGRP%XALBUV_DRY,            &
                        TGRP%XALBNIR_WET, TGRP%XALBVIS_WET, TGRP%XALBUV_WET, TGRP%XBSLAI_NITRO,         &
                        TGRP%XCE_NITRO, TGRP%XCNA_NITRO, TGRP%XCF_NITRO                            )
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
IF (TGRO%CSCOND_GR=='PL98' .OR. TGRO%CISBA_GR=='DIF') THEN
  DO JLAYER=1,TGRO%NLAYER_GR
     TGRP%XHCAPSOIL(:,JLAYER) =    TGRP%XOM_GR(:,JLAYER)  * ZHCAPSOIL_OM +      &
                           (1-TGRP%XOM_GR(:,JLAYER)) * TGRP%XHCAPSOIL(:,JLAYER)  
     TGRP%XCONDDRY(:,JLAYER) = (ZCONDDRY_OM         * TGRP%XCONDDRY(:,JLAYER))    &
                         /(  TGRP%XOM_GR(:,JLAYER)  * TGRP%XCONDDRY(:,JLAYER) +   &
                          (1-TGRP%XOM_GR(:,JLAYER)) * ZCONDDRY_OM)
     TGRP%XCONDSLD(:,JLAYER) = (ZCONDSLD_OM         * TGRP%XCONDSLD(:,JLAYER))    &
                         /(  TGRP%XOM_GR(:,JLAYER)  * TGRP%XCONDSLD(:,JLAYER) +   &
                          (1-TGRP%XOM_GR(:,JLAYER)) * ZCONDSLD_OM)
  ENDDO
END IF
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Validation case : experimental values for Nancy 2011 case
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Substrate layer
DO JLAYER=1,4
    TGRP%XCONDDRY (:,JLAYER) = 0.15
    TGRP%XHCAPSOIL(:,JLAYER) = 1342000.
ENDDO
! Drainage layer
DO JLAYER=5,6
    TGRP%XCONDDRY (:,JLAYER) = 0.09
    TGRP%XHCAPSOIL(:,JLAYER) = 331500.
ENDDO
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!*       5.2     Soil thermal characteristics:
!               --------------------------------
!
DO JLAYER=1,TGRO%NLAYER_GR
  TGRP%XCONDSAT(:,JLAYER) =   TGRP%XOM_GR(:,JLAYER)* ZCONDSAT_OM   &
                        +(1-TGRP%XOM_GR(:,JLAYER))* TGRP%XCONDSAT(:,JLAYER)
END DO
!
! Note that if ISBA/=DIF, always CDIF = 'BC' and CPEDOTF = 'CH78'
DO JLAYER=1,TGRO%NLAYER_GR
  TGRP%XBCOEF  (:,JLAYER) =    TGRP%XOM_GR(:,JLAYER) * ZBCOEF_OM        &
                       +(1-TGRP%XOM_GR(:,JLAYER))* TGRP%XBCOEF(:,JLAYER)
  TGRP%XMPOTSAT(:,JLAYER) =    TGRP%XOM_GR(:,JLAYER) * ZMPOTSAT_OM      &
                       +(1-TGRP%XOM_GR(:,JLAYER))* TGRP%XMPOTSAT(:,JLAYER)
END DO
!        
DO JLAYER=1,TGRO%NLAYER_GR
   TGRP%XWSAT (:,JLAYER) =    TGRP%XOM_GR(:,JLAYER)* ZWSAT_OM            &
                     +(1-TGRP%XOM_GR(:,JLAYER))* TGRP%XWSAT(:,JLAYER)
   TGRP%XWWILT(:,JLAYER) = EXP(((LOG(-1*ZMPOT_WWILT)-LOG(-1*TGRP%XMPOTSAT(:,JLAYER)))   &
                    / (-1*TGRP%XBCOEF(:,JLAYER)))+LOG(TGRP%XWSAT(:,JLAYER)))
   TGRP%XWFC  (:,JLAYER) = EXP(((LOG(ZHYDCOND_WFC)-LOG(TGRP%XCONDSAT(:,JLAYER)))        &
                    / (2*TGRP%XBCOEF(:,JLAYER)+3))+LOG(TGRP%XWSAT(:,JLAYER)))
END DO
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Validation case : experimental values for Nancy 2011 case
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Substrate layer
DO JLAYER=1,4
  TGRP%XWSAT   (:,JLAYER) = 0.674     ! Value tested
  TGRP%XCONDSAT(:,JLAYER) = 2.162E-3  ! Value tested
  TGRP%XMPOTSAT(:,JLAYER) = -0.932    ! Value tested
  TGRP%XBCOEF  (:,JLAYER) = 3.9       ! Value tested
  TGRP%XWWILT  (:,JLAYER) = 0.15      ! from OBS-NANCY
  TGRP%XWFC    (:,JLAYER) = 0.37      ! from OBS-NANCY
ENDDO
! Drainage layer
DO JLAYER=5,6
   TGRP%XWSAT   (:,JLAYER) = 0.9       ! Value tested
   TGRP%XCONDSAT(:,JLAYER) = 3.32E-3   ! Value tested
   TGRP%XMPOTSAT(:,JLAYER) = -0.121    ! Value tested
   TGRP%XBCOEF  (:,JLAYER) = 2.7       ! Value tested
   TGRP%XWWILT  (:,JLAYER) = 0.15      ! sert à initialiser le WG ds la couche
   TGRP%XWFC    (:,JLAYER) = 0.37      ! sert à initialiser le WG ds la couche
ENDDO
!-------------------------------------------------------------------------------
!
!*       6.1    Initialize of the SGH scheme:'
!               ------------------------------
!
IF(TGRO%CKSAT_GR=='SGH' .AND. TGRO%CISBA_GR/='DIF' .AND. HINIT/='PRE')THEN 
  ZF(:)=MIN(4.0/TGRP%XDG(:,2),XF_DECAY)
  CALL EXP_DECAY_SOIL_FR(TGRO%CISBA_GR, ZF(:),TGRP%XC1SAT(:),TGRP%XC2REF(:),TGRP%XDG(:,:),TGRP%XD_ICE(:),&
                         TGRP%XC4REF(:),TGRP%XC3(:,:),TGRP%XCONDSAT(:,:),TGRP%XKSAT_ICE(:))
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
