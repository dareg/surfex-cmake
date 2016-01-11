!#############################################################
SUBROUTINE INIT_TEB_GREENROOF_PGD_n (DTCO, U, CHI, DTI, I, DST, SLT, CHT, TG, T, TOP, GDO, GRM, &
                                     HPROGRAM,HINIT,OREAD_PGD, KI, KSV, HSV, KVERSION, PCO2, PRHOA)
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
USE MODD_CH_ISBA_n, ONLY : CH_ISBA_t
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_ISBA_n, ONLY : ISBA_t
USE MODD_DST_n, ONLY : DST_t
USE MODD_SLT_n, ONLY : SLT_t
USE MODD_CH_TEB_n, ONLY : CH_TEB_t
USE MODD_GRID_n, ONLY : GRID_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_SURFEX_n, ONLY : TEB_GREENROOF_MODEL_t
!
USE MODD_TYPE_DATE_SURF
USE MODD_TYPE_SNOW
!


USE MODD_DATA_COVER_PAR,       ONLY: NVEGTYPE
USE MODD_SURF_PAR,             ONLY: XUNDEF, NUNDEF

USE MODD_SGH_PAR,              ONLY: NDIMTAB, XF_DECAY
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
TYPE(CH_ISBA_t), INTENT(INOUT) :: CHI
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTI
TYPE(ISBA_t), INTENT(INOUT) :: I
TYPE(DST_t), INTENT(INOUT) :: DST
TYPE(SLT_t), INTENT(INOUT) :: SLT
TYPE(CH_TEB_t), INTENT(INOUT) :: CHT
TYPE(GRID_t), INTENT(INOUT) :: TG
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: GDO
TYPE(TEB_GREENROOF_MODEL_t), INTENT(INOUT) :: GRM
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
!
        REAL, DIMENSION(:,:), POINTER :: ZPATCH
        REAL, DIMENSION(:,:,:), POINTER :: ZVEGTYPE_PATCH
        INTEGER, DIMENSION(:), POINTER :: ISIZE_NATURE_P
        INTEGER, DIMENSION(:,:), POINTER :: IR_NATURE_P
        REAL, DIMENSION(0) :: ZTDEEP_CLI, ZGAMMAT_CLI, ZTHRESHOLD
        INTEGER, DIMENSION(:,:), POINTER :: IIRRINUM
        LOGICAL, DIMENSION(:,:), POINTER :: GIRRIDAY
        LOGICAL, DIMENSION(:,:), POINTER :: GIRRIGATE
        REAL, DIMENSION(:,:), POINTER :: ZTHRESHOLDSPT
        REAL, DIMENSION(:,:,:), POINTER :: ZINCREASE
        REAL, DIMENSION(:,:,:), POINTER :: ZTURNOVER
        REAL, DIMENSION(:,:,:), POINTER :: ZSFDST
        REAL, DIMENSION(:,:,:), POINTER :: ZSFDSTM
        REAL, DIMENSION(:,:,:), POINTER :: ZSFSLT
        REAL, DIMENSION(0) :: ZAOSIP, ZAOSIM, ZAOSJP, ZAOSJM
        REAL, DIMENSION(0) :: ZHO2IP, ZHO2IM, ZHO2JP, ZHO2JM
        REAL, DIMENSION(0,0) :: ZZ0
        REAL, DIMENSION(:,:), POINTER :: ZZ0EFFIP
        REAL, DIMENSION(:,:), POINTER :: ZZ0EFFIM
        REAL, DIMENSION(:,:), POINTER :: ZZ0EFFJP
        REAL, DIMENSION(:,:), POINTER :: ZZ0EFFJM
        REAL, DIMENSION(:), POINTER :: ZZ0REL
        REAL, DIMENSION(:,:), POINTER :: ZTAU_WOOD
        REAL, DIMENSION(:,:), POINTER :: ZKANISO
        REAL, DIMENSION(:,:), POINTER :: ZWD0
        REAL, DIMENSION(:), POINTER   :: ZFWTD ! grid-cell fraction of water table to rise
        REAL, DIMENSION(:), POINTER   :: ZWTD  ! water table depth from Obs, TRIP or MODCOU! = 0.1 mm/day
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
GRM%TV%O%NNBIOMASS = GDO%NNBIOMASS
GRM%TV%O%CPHOTO = GDO%CPHOTO
GRM%TV%O%CPEDOTF = GDO%CPEDOTF
GRM%TV%O%CRESPSL = GDO%CRESPSL
!
!*       2.1    Cover, soil and orographic fields:
!               ---------------------------------
!
IF (OREAD_PGD) &
 CALL READ_PGD_TEB_GREENROOF_n(CHT, DTCO, GRM%DTI, GRM%GB, U, GRM%TV%O, GRM%TV%P, GRM%TV%IP, TG, &
                               HPROGRAM,KVERSION)
!
!
!* allocation of green roofs variables
!
 CALL ALLOCATE_TEB_GREENROOF_PGD(GRM%TV%M, GRM%TV%P, GRM%TV%IP, &
                                 OREAD_PGD, KI, NVEGTYPE, GRM%TV%O%NGROUND_LAYER, NDIMTAB)
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
IF (.NOT. GRM%TV%O%LPAR) THEN
          CALL CONVERT_PATCH_ISBA(DTCO, DTI, GRM%TV%O, &
                        GRM%TV%O%CISBA,IDECADE,IDECADE,TOP%XCOVER,TOP%LCOVER,&
                        GRM%TV%O%CPHOTO,.FALSE.,  &
                        .FALSE.,GRM%TV%O%LTR_ML,'GRD',PVEG=GRM%TV%M%T%CUR%XVEG,PLAI=GRM%TV%M%T%CUR%XLAI,              &
                        PRSMIN=GRM%TV%M%T%CUR%XRSMIN,PGAMMA=GRM%TV%M%T%CUR%XGAMMA,PWRMAX_CF=GRM%TV%M%T%CUR%XWRMAX_CF,       &
                        PRGL=GRM%TV%M%T%CUR%XRGL,PCV=GRM%TV%M%T%CUR%XCV,PSOILGRID=GRM%TV%O%XSOILGRID,                 &
                        PDG=GRM%TV%M%X%XDG,KWG_LAYER=GRM%TV%M%X%NWG_LAYER,PDROOT=GRM%TV%M%X%XDROOT,PDG2=GRM%TV%M%X%XDG2,   &
                        PZ0=GRM%TV%M%T%CUR%XZ0,PZ0_O_Z0H=GRM%TV%M%X%XZ0_O_Z0H,PPERM=GRM%TV%P%XPERM,         &
                        PALBNIR_VEG=GRM%TV%M%T%CUR%XALBNIR_VEG,PALBVIS_VEG=GRM%TV%M%T%CUR%XALBVIS_VEG,       &
                        PALBUV_VEG=GRM%TV%M%T%CUR%XALBUV_VEG,PEMIS_ECO=GRM%TV%M%T%CUR%XEMIS,                 &
                        PVEGTYPE=GRM%TV%M%X%XVEGTYPE,PROOTFRAC=GRM%TV%M%X%XROOTFRAC,                 &
                        PGMES=GRM%TV%M%T%CUR%XGMES,PBSLAI=GRM%TV%M%T%CUR%XBSLAI,PLAIMIN=GRM%TV%M%T%CUR%XLAIMIN,             &
                        PSEFOLD=GRM%TV%M%T%CUR%XSEFOLD,PGC=GRM%TV%M%T%CUR%XGC,                               &
                        PDMAX=GRM%TV%M%X%XDMAX,PF2I=GRM%TV%M%T%CUR%XF2I,OSTRESS=GRM%TV%M%T%CUR%LSTRESS,PH_TREE=GRM%TV%M%X%XH_TREE, &
                        PRE25=GRM%TV%M%X%XRE25,PCE_NITRO=GRM%TV%M%T%CUR%XCE_NITRO,PCF_NITRO=GRM%TV%M%T%CUR%XCF_NITRO,   &
                        PCNA_NITRO=GRM%TV%M%T%CUR%XCNA_NITRO,PD_ICE=GRM%TV%M%X%XD_ICE                    )   
ELSE
 CALL INIT_FROM_DATA_GREENROOF_n(GRM%DTI, IDECADE,GRM%TV%O%CPHOTO,        &
                                 GRM%TV%M%T%CUR%XVEG, GRM%TV%M%T%CUR%XLAI,GRM%TV%M%T%CUR%XRSMIN, &
                                 GRM%TV%M%T%CUR%XGAMMA,GRM%TV%M%T%CUR%XWRMAX_CF,GRM%TV%M%T%CUR%XRGL,&
                                 GRM%TV%M%T%CUR%XCV,&
                                 GRM%TV%M%X%XDG,GRM%TV%M%X%XD_ICE,GRM%TV%M%T%CUR%XZ0,GRM%TV%M%X%XZ0_O_Z0H,  &
                                 GRM%TV%M%T%CUR%XALBNIR_VEG,GRM%TV%M%T%CUR%XALBVIS_VEG,GRM%TV%M%T%CUR%XALBUV_VEG,&
                                 GRM%TV%M%T%CUR%XEMIS,GRM%TV%M%X%XVEGTYPE,GRM%TV%M%X%XROOTFRAC,          &
                                 GRM%TV%M%T%CUR%XGMES,GRM%TV%M%T%CUR%XBSLAI,GRM%TV%M%T%CUR%XLAIMIN,&
                                 GRM%TV%M%T%CUR%XSEFOLD,GRM%TV%M%T%CUR%XGC,   &
                                 GRM%TV%M%X%XDMAX, GRM%TV%M%T%CUR%XF2I, GRM%TV%M%T%CUR%LSTRESS, &
                                 GRM%TV%M%X%XH_TREE,GRM%TV%M%X%XRE25,&
                                 GRM%TV%M%T%CUR%XCE_NITRO,GRM%TV%M%T%CUR%XCF_NITRO,GRM%TV%M%T%CUR%XCNA_NITRO      )  
  IF (GRM%TV%O%CISBA=='DIF') THEN
    WHERE(T%CUR%XGREENROOF(:)/=0.)
      GRM%TV%M%X%NWG_LAYER(:,1)=GRM%TV%O%NGROUND_LAYER
      GRM%TV%M%X%XDG2  (:,1)=0.0
      GRM%TV%M%X%XDROOT(:,1)=0.0
    ENDWHERE
    DO JLAYER=GRM%TV%O%NGROUND_LAYER,1,-1
      DO JILU=1,KI
        IF(T%CUR%XGREENROOF(JILU)/=0..AND.GRM%TV%M%X%XROOTFRAC(JILU,JLAYER,1)>=1.0)THEN
          GRM%TV%M%X%XDG2  (JILU,1)=GRM%TV%M%X%XDG(JILU,JLAYER,1)
          GRM%TV%M%X%XDROOT(JILU,1)=GRM%TV%M%X%XDG(JILU,JLAYER,1)
        ENDIF
      ENDDO
    ENDDO
  ENDIF
END IF
!
WHERE (T%CUR%XGREENROOF(:)==0.)
  ! GARDEN default values /may need changing for green roofs
  GRM%TV%P%XSOC    (:,1) = 0.5
  GRM%TV%P%XSOC    (:,2) = 0.5
  GRM%TV%P%XSAND   (:,1) = 0.33
  GRM%TV%P%XSAND   (:,2) = 0.33
  GRM%TV%P%XCLAY   (:,1) = 0.33
  GRM%TV%P%XCLAY   (:,2) = 0.33
  GRM%TV%M%T%CUR%XVEG       (:  ,1) = 0.
  GRM%TV%M%T%CUR%XLAI       (:  ,1) = 0.
  GRM%TV%M%T%CUR%XRSMIN     (:  ,1) = 40.
  GRM%TV%M%T%CUR%XGAMMA     (:  ,1) = 0.
  GRM%TV%M%T%CUR%XWRMAX_CF  (:  ,1) = 0.2
  GRM%TV%M%T%CUR%XRGL       (:  ,1) = 100.
  GRM%TV%M%T%CUR%XCV        (:  ,1) = 2.E-5
  GRM%TV%M%T%CUR%XZ0        (: ,1 ) = 0.013
  GRM%TV%M%X%XZ0_O_Z0H  (: ,1 ) = 10.
  GRM%TV%M%T%CUR%XALBNIR_VEG(: ,1 ) = 0.30
  GRM%TV%M%T%CUR%XALBVIS_VEG(:  ,1) = 0.30
  GRM%TV%M%T%CUR%XALBUV_VEG (:  ,1) = 0.06
  GRM%TV%M%T%CUR%XEMIS      (: ,1 ) = 0.94
END WHERE
IF (GRM%TV%O%CPHOTO/='NON') THEN
  WHERE (T%CUR%XGREENROOF(:)==0.)
    GRM%TV%M%T%CUR%XGMES      (:  ,1) = 0.020
    GRM%TV%M%T%CUR%XBSLAI     (:  ,1) = 0.36
    GRM%TV%M%T%CUR%XLAIMIN    (:  ,1) = 0.3
    GRM%TV%M%T%CUR%XSEFOLD    (:  ,1) = 90*86400.
    GRM%TV%M%X%XH_TREE    (: ,1 ) = 0.
    GRM%TV%M%X%XRE25      (: ,1 ) = 3.6E-7    
    GRM%TV%M%T%CUR%XGC        (: ,1 ) = 0.00025
  END WHERE
  IF (GRM%TV%O%CPHOTO/='AGS' .AND. GRM%TV%O%CPHOTO/='LAI') THEN
    WHERE (T%CUR%XGREENROOF(:)==0.)     
      GRM%TV%M%X%XDMAX      (:  ,1) = 0.1
      GRM%TV%M%T%CUR%XF2I       (: ,1 ) = 0.3
    END WHERE
    IF (GRM%TV%O%CPHOTO=='NIT' .OR. GRM%TV%O%CPHOTO=='NCB') THEN
      WHERE (T%CUR%XGREENROOF(:)==0.)          
        GRM%TV%M%T%CUR%XCE_NITRO  (: ,1 ) = 7.68
        GRM%TV%M%T%CUR%XCF_NITRO  (: ,1 ) = -4.33
        GRM%TV%M%T%CUR%XCNA_NITRO (: ,1 ) = 1.3
      END WHERE
    ENDIF
  ENDIF
ENDIF  
IF(GRM%TV%O%CISBA/='DIF')THEN
  DO JLAYER=1,GRM%TV%O%NGROUND_LAYER
    WHERE (T%CUR%XGREENROOF(:)==0.)
      GRM%TV%M%X%XDG(:,JLAYER,1)=0.2*JLAYER
    END WHERE
  ENDDO
ELSE
  WHERE (T%CUR%XGREENROOF(:)==0.) 
    GRM%TV%M%X%XDG(:,1,1)=0.01
    GRM%TV%M%X%XDG(:,2,1)=0.04
    GRM%TV%M%X%XROOTFRAC(:,1,1)=0.
    GRM%TV%M%X%XROOTFRAC(:,2,1)=0.
  END WHERE        
  DO JLAYER=3,GRM%TV%O%NGROUND_LAYER
    WHERE (T%CUR%XGREENROOF(:)==0.)
      GRM%TV%M%X%XDG(:,JLAYER,1)=0.1*(JLAYER-2)
      GRM%TV%M%X%XROOTFRAC(:,JLAYER,1)=0.
    END WHERE
  ENDDO               
  WHERE (T%CUR%XGREENROOF(:)==0.) 
    GRM%TV%M%X%NWG_LAYER(:,1)=GRM%TV%O%NGROUND_LAYER
    GRM%TV%M%X%XDROOT   (:,1)=0.0
    GRM%TV%M%X%XDG2     (:,1)=GRM%TV%M%X%XDG(:,GRM%TV%O%NGROUND_LAYER-1,1)
  ENDWHERE    
ENDIF  
WHERE (T%CUR%XGREENROOF(:)==0.) 
  GRM%TV%M%X%XD_ICE(:,1)=0.8*GRM%TV%M%X%XDG(:,2,1)
END WHERE  
DO JVEGTYPE=1,NVEGTYPE
  WHERE (T%CUR%XGREENROOF(:)==0.)
    GRM%TV%M%X%XVEGTYPE(:,JVEGTYPE)=0.
    GRM%TV%M%X%XVEGTYPE(:,1)=1.
  END WHERE
ENDDO
!
        NULLIFY(ZPATCH)
        NULLIFY(ZVEGTYPE_PATCH)
        NULLIFY(ISIZE_NATURE_P)
        NULLIFY(IR_NATURE_P)
        NULLIFY(IIRRINUM)
        NULLIFY(GIRRIDAY)
        NULLIFY(GIRRIGATE)
        NULLIFY(ZTHRESHOLDSPT)
        NULLIFY(ZINCREASE)
        NULLIFY(ZTURNOVER)
        NULLIFY(ZSFDST)
        NULLIFY(ZSFDSTM)
        NULLIFY(ZSFSLT)
        NULLIFY(ZZ0EFFIP)
        NULLIFY(ZZ0EFFIM)
        NULLIFY(ZZ0EFFJP)
        NULLIFY(ZZ0EFFJM)
        NULLIFY(ZZ0REL)
        NULLIFY(ZFWTD)
        NULLIFY(ZWTD)
        NULLIFY(ZWD0)
        NULLIFY(ZTAU_WOOD)
        NULLIFY(ZKANISO)

         CALL INIT_VEG_PGD_n(CHI, DTCO, DST, SLT, U, &
                     GRM%TV%O%LAGRI_TO_GRASS, TOP%LCOVER, TOP%XCOVER, &
                     HPROGRAM, 'TOWN  ',ILUOUT, KI, 1, GRM%TV%O%NGROUND_LAYER, TOP%TTIME%TDATE%MONTH,    &
                     GRM%TV%M%X%XVEGTYPE,ZPATCH, ZVEGTYPE_PATCH, ISIZE_NATURE_P, IR_NATURE_P,    &
                  0.0, .FALSE., .FALSE., ZTDEEP_CLI, ZGAMMAT_CLI, &
                  GRM%TV%IP%XTDEEP, GRM%TV%IP%XGAMMAT,  .FALSE., &
                  ZTHRESHOLD, IIRRINUM, GIRRIDAY, GIRRIGATE, ZTHRESHOLDSPT, &
                  GRM%TV%O%CPHOTO, HINIT,GRM%TV%O%LTR_ML, GRM%TV%O%NNBIOMASS, &
                  PCO2, PRHOA, GRM%TV%IP%XABC, GRM%TV%IP%XPOI,  &
                  GRM%TV%M%T%CUR%XGMES, GRM%TV%M%T%CUR%XGC, GRM%TV%M%X%XDMAX, &
                  GRM%TV%IP%XANMAX, GRM%TV%IP%XFZERO, GRM%TV%IP%XEPSO, GRM%TV%IP%XGAMM, GRM%TV%IP%XQDGAMM,   &
                  GRM%TV%IP%XQDGMES, GRM%TV%IP%XT1GMES, GRM%TV%IP%XT2GMES, GRM%TV%IP%XAMAX, GRM%TV%IP%XQDAMAX, &
                  GRM%TV%IP%XT1AMAX, GRM%TV%IP%XT2AMAX,GRM%TV%IP%XAH, GRM%TV%IP%XBH,&           
                  ZTAU_WOOD, ZINCREASE, ZTURNOVER,                  &
                  KSV, HSV, CHT%SVT, CHT%CCH_NAMES, CHT%CAER_NAMES,CHT%CDSTNAMES, CHT%CSLTNAMES, &
                  CHT%CCHEM_SURF_FILE,                     &
                  ZSFDST, ZSFDSTM, ZSFSLT,                                    &
                  ZAOSIP, ZAOSIM, ZAOSJP, ZAOSJM, ZHO2IP, ZHO2IM, ZHO2JP,     &
                  ZHO2JM, ZZ0, ZZ0EFFIP, ZZ0EFFIM, ZZ0EFFJP, ZZ0EFFJM, ZZ0REL,&
                  GRM%TV%P%XCLAY, GRM%TV%P%XSAND, GRM%TV%O%CPEDOTF,      &
                  GRM%TV%IP%XCONDSAT, GRM%TV%IP%XMPOTSAT, GRM%TV%IP%XBCOEF, GRM%TV%IP%XWWILT, &
                  GRM%TV%IP%XWFC, GRM%TV%IP%XWSAT, ZWD0, ZKANISO, GRM%TV%O%CRUNOFF,  &
                  GRM%TV%IP%XTAUICE, GRM%TV%IP%XCGSAT, GRM%TV%IP%XC1SAT, &
                  GRM%TV%IP%XC2REF, GRM%TV%IP%XC3, GRM%TV%IP%XC4B, GRM%TV%IP%XACOEF, GRM%TV%IP%XPCOEF, &
                  GRM%TV%IP%XC4REF, GRM%TV%IP%XPCPS, GRM%TV%IP%XPLVTT, GRM%TV%IP%XPLSTT,        &
                  GRM%TV%O%CSCOND, GRM%TV%O%CISBA, GRM%TV%IP%XHCAPSOIL, GRM%TV%IP%XCONDDRY, &
                  GRM%TV%IP%XCONDSLD, GRM%TV%O%CCPSURF, GRM%TV%M%X%XDG, GRM%TV%M%X%XDROOT, GRM%TV%M%X%XDG2, &
                  GRM%TV%M%X%XROOTFRAC, GRM%TV%IP%XRUNOFFD, GRM%TV%IP%XDZG, GRM%TV%IP%XDZDIF,       &
                  GRM%TV%IP%XSOILWGHT, GRM%TV%M%X%NWG_LAYER, GRM%TV%O%NLAYER_HORT, &
                  GRM%TV%O%NLAYER_DUN, GRM%TV%M%X%XD_ICE,  &
                  GRM%TV%IP%XKSAT_ICE, GRM%TV%IP%XALBNIR_DRY, GRM%TV%IP%XALBVIS_DRY, GRM%TV%IP%XALBUV_DRY,   &
                  GRM%TV%IP%XALBNIR_WET, GRM%TV%IP%XALBVIS_WET, GRM%TV%IP%XALBUV_WET, GRM%TV%IP%XBSLAI_NITRO, &
                  GRM%TV%M%T%CUR%XCE_NITRO, GRM%TV%M%T%CUR%XCNA_NITRO, GRM%TV%M%T%CUR%XCF_NITRO,          &
                  ZFWTD, ZWTD               )
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
IF (GRM%TV%O%CSCOND=='PL98' .OR. GRM%TV%O%CISBA=='DIF') THEN
  DO JLAYER=1,GRM%TV%O%NGROUND_LAYER
     GRM%TV%IP%XHCAPSOIL(:,JLAYER) =    GRM%TV%P%XSOC(:,JLAYER)  * ZHCAPSOIL_OM +      &
                           (1-GRM%TV%P%XSOC(:,JLAYER)) * GRM%TV%IP%XHCAPSOIL(:,JLAYER)  
     GRM%TV%IP%XCONDDRY(:,JLAYER) = (ZCONDDRY_OM         * GRM%TV%IP%XCONDDRY(:,JLAYER))    &
                         /(  GRM%TV%P%XSOC(:,JLAYER)  * GRM%TV%IP%XCONDDRY(:,JLAYER) +   &
                          (1-GRM%TV%P%XSOC(:,JLAYER)) * ZCONDDRY_OM)
     GRM%TV%IP%XCONDSLD(:,JLAYER) = (ZCONDSLD_OM         * GRM%TV%IP%XCONDSLD(:,JLAYER))    &
                         /(  GRM%TV%P%XSOC(:,JLAYER)  * GRM%TV%IP%XCONDSLD(:,JLAYER) +   &
                          (1-GRM%TV%P%XSOC(:,JLAYER)) * ZCONDSLD_OM)
  ENDDO
END IF
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Validation case : experimental values for Nancy 2011 case
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Substrate layer
DO JLAYER=1,4
    GRM%TV%IP%XCONDDRY (:,JLAYER) = 0.15
    GRM%TV%IP%XHCAPSOIL(:,JLAYER) = 1342000.
ENDDO
! Drainage layer
DO JLAYER=5,6
    GRM%TV%IP%XCONDDRY (:,JLAYER) = 0.09
    GRM%TV%IP%XHCAPSOIL(:,JLAYER) = 331500.
ENDDO
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!*       5.2     Soil thermal characteristics:
!               --------------------------------
!
DO JLAYER=1,GRM%TV%O%NGROUND_LAYER
  GRM%TV%IP%XCONDSAT(:,JLAYER,1) =   GRM%TV%P%XSOC(:,JLAYER)* ZCONDSAT_OM   &
                        +(1-GRM%TV%P%XSOC(:,JLAYER))* GRM%TV%IP%XCONDSAT(:,JLAYER,1)
END DO
!
! Note that if ISBA/=DIF, always CDIF = 'BC' and CPEDOTF = 'CH78'
DO JLAYER=1,GRM%TV%O%NGROUND_LAYER
  GRM%TV%IP%XBCOEF  (:,JLAYER) =    GRM%TV%P%XSOC(:,JLAYER) * ZBCOEF_OM        &
                       +(1-GRM%TV%P%XSOC(:,JLAYER))* GRM%TV%IP%XBCOEF(:,JLAYER)
  GRM%TV%IP%XMPOTSAT(:,JLAYER) =    GRM%TV%P%XSOC(:,JLAYER) * ZMPOTSAT_OM      &
                       +(1-GRM%TV%P%XSOC(:,JLAYER))* GRM%TV%IP%XMPOTSAT(:,JLAYER)
END DO
!        
DO JLAYER=1,GRM%TV%O%NGROUND_LAYER
   GRM%TV%IP%XWSAT (:,JLAYER) =    GRM%TV%P%XSOC(:,JLAYER)* ZWSAT_OM            &
                     +(1-GRM%TV%P%XSOC(:,JLAYER))* GRM%TV%IP%XWSAT(:,JLAYER)
   GRM%TV%IP%XWWILT(:,JLAYER) = EXP(((LOG(-1*ZMPOT_WWILT)-LOG(-1*GRM%TV%IP%XMPOTSAT(:,JLAYER)))   &
                    / (-1*GRM%TV%IP%XBCOEF(:,JLAYER)))+LOG(GRM%TV%IP%XWSAT(:,JLAYER)))
   GRM%TV%IP%XWFC  (:,JLAYER) = EXP(((LOG(ZHYDCOND_WFC)-LOG(GRM%TV%IP%XCONDSAT(:,JLAYER,1)))        &
                    / (2*GRM%TV%IP%XBCOEF(:,JLAYER)+3))+LOG(GRM%TV%IP%XWSAT(:,JLAYER)))
END DO
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Validation case : experimental values for Nancy 2011 case
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Substrate layer
DO JLAYER=1,4
  GRM%TV%IP%XWSAT   (:,JLAYER) = 0.674     ! Value tested
  GRM%TV%IP%XCONDSAT(:,JLAYER,1) = 2.162E-3  ! Value tested
  GRM%TV%IP%XMPOTSAT(:,JLAYER) = -0.932    ! Value tested
  GRM%TV%IP%XBCOEF  (:,JLAYER) = 3.9       ! Value tested
  GRM%TV%IP%XWWILT  (:,JLAYER) = 0.15      ! from OBS-NANCY
  GRM%TV%IP%XWFC    (:,JLAYER) = 0.37      ! from OBS-NANCY
ENDDO
! Drainage layer
DO JLAYER=5,6
   GRM%TV%IP%XWSAT   (:,JLAYER) = 0.9       ! Value tested
   GRM%TV%IP%XCONDSAT(:,JLAYER,1) = 3.32E-3   ! Value tested
   GRM%TV%IP%XMPOTSAT(:,JLAYER) = -0.121    ! Value tested
   GRM%TV%IP%XBCOEF  (:,JLAYER) = 2.7       ! Value tested
   GRM%TV%IP%XWWILT  (:,JLAYER) = 0.15      ! sert à initialiser le WG ds la couche
   GRM%TV%IP%XWFC    (:,JLAYER) = 0.37      ! sert à initialiser le WG ds la couche
ENDDO
!-------------------------------------------------------------------------------
!
!*       6.1    Initialize of the SGH scheme:'
!               ------------------------------
!
IF(GRM%TV%O%CKSAT=='SGH' .AND. GRM%TV%O%CISBA/='DIF' .AND. HINIT/='PRE')THEN 
  ZF(:)=MIN(4.0/GRM%TV%M%X%XDG(:,2,1),XF_DECAY)
  CALL EXP_DECAY_SOIL_FR(GRM%TV%O%CISBA, ZF(:),GRM%TV%IP%XC1SAT(:,1),GRM%TV%IP%XC2REF(:,1),&
                         GRM%TV%M%X%XDG(:,:,1),GRM%TV%M%X%XD_ICE(:,1),GRM%TV%IP%XC4REF(:,1),&
                         GRM%TV%IP%XC3(:,:,1),GRM%TV%IP%XCONDSAT(:,:,1),GRM%TV%IP%XKSAT_ICE(:,1))
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
