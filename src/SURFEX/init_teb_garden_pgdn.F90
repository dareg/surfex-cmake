!#############################################################
        SUBROUTINE INIT_TEB_GARDEN_PGD_n (DTCO, U, CHI, DTI, I, DST, SLT, CHT, TG, T, TOP, GDM, &
                                          HPROGRAM, HINIT, OREAD_PGD,KI, KSV, HSV, KVERSION, KBUGFIX, &
                                          PCO2, PRHOA)
        !#############################################################
        !
        !!****  *INIT_TEB_GARDEN_PGD_n* - routine to initialize ISBA
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
        USE MODD_SURFEX_n, ONLY : TEB_GARDEN_MODEL_t
        !
        USE MODD_TYPE_DATE_SURF
        USE MODD_TYPE_SNOW
        !
                                        
        USE MODD_DATA_COVER_PAR,  ONLY: NVEGTYPE
        USE MODD_SURF_PAR,        ONLY: XUNDEF, NUNDEF

        USE MODD_SGH_PAR,         ONLY: NDIMTAB, XF_DECAY
        !
        USE MODI_GET_LUOUT
        USE MODI_ALLOCATE_TEB_GARDEN_PGD
        USE MODI_READ_PGD_TEB_GARDEN_n
        USE MODI_CONVERT_PATCH_ISBA
        USE MODI_INIT_FROM_DATA_GRDN_n
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
        TYPE(TEB_GARDEN_MODEL_t), INTENT(INOUT) :: GDM
        !
         CHARACTER(LEN=6),                   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
         CHARACTER(LEN=3),                   INTENT(IN)  :: HINIT     ! choice of fields to initialize
        LOGICAL,                            INTENT(IN)  :: OREAD_PGD ! flag to read PGD fields in the file
        INTEGER,                            INTENT(IN)  :: KI        ! number of points
        INTEGER,                            INTENT(IN)  :: KSV       ! number of scalars
         CHARACTER(LEN=6), DIMENSION(KSV),   INTENT(IN)  :: HSV       ! name of all scalar variables
        INTEGER,                            INTENT(IN)  :: KVERSION  ! version number of the file being read
        INTEGER,                            INTENT(IN)  :: KBUGFIX
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
        INTEGER :: JVEGTYPE, JLAYER  ! loop counter on vegtypes
        !
        REAL, DIMENSION(KI)               :: ZF
        REAL, DIMENSION(KI)               :: ZWORK
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
        REAL, DIMENSION(:), POINTER   :: ZWTD  ! water table depth from Obs, TRIP or MODCOU
        !
        REAL(KIND=JPRB) :: ZHOOK_HANDLE
        !
        !-------------------------------------------------------------------------------
        !
        !               Initialisation for IO
        !
        IF (LHOOK) CALL DR_HOOK('INIT_TEB_GARDEN_PGD_n',0,ZHOOK_HANDLE)
         CALL GET_LUOUT(HPROGRAM,ILUOUT)
        !
        !-------------------------------------------------------------------------------
        !
        !*       2.     Physiographic fields
        !               --------------------
        !
        !* allocation of urban green area variables
        !
         CALL ALLOCATE_TEB_GARDEN_PGD(GDM%TV%M, GDM%TV%P, GDM%TV%IP, &
                                      OREAD_PGD, KI, NVEGTYPE, GDM%TV%O%NGROUND_LAYER, NDIMTAB)  
        !
        !
        !*       2.1    Cover, soil and orographic fields:
        !               ---------------------------------
        !
        IF (OREAD_PGD) &
         CALL READ_PGD_TEB_GARDEN_n(CHT, DTCO, GDM%DTI, GDM%GB, U, GDM%TV%O, GDM%TV%P, GDM%TV%IP, TG, TOP, &
                                    HPROGRAM,KVERSION,KBUGFIX)
        !
        !
        !*       2.3    Physiographic data fields from land cover:
        !               -----------------------------------------
        !
        IF (TOP%TTIME%TDATE%MONTH /= NUNDEF) THEN
          IDECADE = 3 * ( TOP%TTIME%TDATE%MONTH - 1 ) + MIN(TOP%TTIME%TDATE%DAY-1,29) / 10 + 1
        ELSE
          IDECADE = 1
        END IF
        !
        !
        IF (.NOT. GDM%TV%O%LPAR) THEN
          CALL CONVERT_PATCH_ISBA(DTCO, DTI, GDM%TV%O, &
                        GDM%TV%O%CISBA,IDECADE,IDECADE,TOP%XCOVER,TOP%LCOVER,&
                        GDM%TV%O%CPHOTO,.FALSE.,  &
                        .FALSE.,GDM%TV%O%LTR_ML,'GRD',PVEG=GDM%TV%M%T%CUR%XVEG,PLAI=GDM%TV%M%T%CUR%XLAI,              &
                        PRSMIN=GDM%TV%M%T%CUR%XRSMIN,PGAMMA=GDM%TV%M%T%CUR%XGAMMA,PWRMAX_CF=GDM%TV%M%T%CUR%XWRMAX_CF,       &
                        PRGL=GDM%TV%M%T%CUR%XRGL,PCV=GDM%TV%M%T%CUR%XCV,PSOILGRID=GDM%TV%O%XSOILGRID,                 &
                        PDG=GDM%TV%M%X%XDG,KWG_LAYER=GDM%TV%M%X%NWG_LAYER,PDROOT=GDM%TV%M%X%XDROOT,PDG2=GDM%TV%M%X%XDG2,   &
                        PZ0=GDM%TV%M%T%CUR%XZ0,PZ0_O_Z0H=GDM%TV%M%X%XZ0_O_Z0H,PPERM=GDM%TV%P%XPERM,         &
                        PALBNIR_VEG=GDM%TV%M%T%CUR%XALBNIR_VEG,PALBVIS_VEG=GDM%TV%M%T%CUR%XALBVIS_VEG,       &
                        PALBUV_VEG=GDM%TV%M%T%CUR%XALBUV_VEG,PEMIS_ECO=GDM%TV%M%T%CUR%XEMIS,                 &
                        PVEGTYPE=GDM%TV%M%X%XVEGTYPE,PROOTFRAC=GDM%TV%M%X%XROOTFRAC,                 &
                        PGMES=GDM%TV%M%T%CUR%XGMES,PBSLAI=GDM%TV%M%T%CUR%XBSLAI,PLAIMIN=GDM%TV%M%T%CUR%XLAIMIN,             &
                        PSEFOLD=GDM%TV%M%T%CUR%XSEFOLD,PGC=GDM%TV%M%T%CUR%XGC,                               &
                        PDMAX=GDM%TV%M%X%XDMAX,PF2I=GDM%TV%M%T%CUR%XF2I,OSTRESS=GDM%TV%M%T%CUR%LSTRESS,PH_TREE=GDM%TV%M%X%XH_TREE, &
                        PRE25=GDM%TV%M%X%XRE25,PCE_NITRO=GDM%TV%M%T%CUR%XCE_NITRO,PCF_NITRO=GDM%TV%M%T%CUR%XCF_NITRO,   &
                        PCNA_NITRO=GDM%TV%M%T%CUR%XCNA_NITRO,PD_ICE=GDM%TV%M%X%XD_ICE                    )   

        ELSE
         CALL INIT_FROM_DATA_GRDN_n(GDM%DTI, &
                                    IDECADE,GDM%TV%O%CPHOTO, GDM%TV%M%T%CUR%XVEG, &
                                    GDM%TV%M%T%CUR%XLAI,GDM%TV%M%T%CUR%XRSMIN,GDM%TV%M%T%CUR%XGAMMA,&
                                    GDM%TV%M%T%CUR%XWRMAX_CF, GDM%TV%M%T%CUR%XRGL,GDM%TV%M%T%CUR%XCV,GDM%TV%M%X%XDG,&
                                    GDM%TV%M%X%XD_ICE,GDM%TV%M%T%CUR%XZ0,GDM%TV%M%X%XZ0_O_Z0H,&
                                    GDM%TV%M%T%CUR%XALBNIR_VEG,GDM%TV%M%T%CUR%XALBVIS_VEG,     &
                                    GDM%TV%M%T%CUR%XALBUV_VEG,GDM%TV%M%T%CUR%XEMIS,      &
                                    GDM%TV%M%X%XVEGTYPE,GDM%TV%M%X%XROOTFRAC,GDM%TV%M%T%CUR%XGMES,&
                                    GDM%TV%M%T%CUR%XBSLAI,GDM%TV%M%T%CUR%XLAIMIN,GDM%TV%M%T%CUR%XSEFOLD,GDM%TV%M%T%CUR%XGC,   &
                                    GDM%TV%M%X%XDMAX, GDM%TV%M%T%CUR%XF2I, GDM%TV%M%T%CUR%LSTRESS, GDM%TV%M%X%XH_TREE,&
                                    GDM%TV%M%X%XRE25,GDM%TV%M%T%CUR%XCE_NITRO,GDM%TV%M%T%CUR%XCF_NITRO,&
                                    GDM%TV%M%T%CUR%XCNA_NITRO      )  

          IF (GDM%TV%O%CISBA=='DIF') THEN
            WHERE(T%CUR%XGARDEN(:)/=0.)
              GDM%TV%M%X%NWG_LAYER(:,1)=GDM%TV%O%NGROUND_LAYER 
              GDM%TV%M%X%XDG2  (:,1)=0.0
              GDM%TV%M%X%XDROOT(:,1)=0.0
            ENDWHERE
            DO JLAYER=GDM%TV%O%NGROUND_LAYER,1,-1
              DO JILU=1,KI
                IF(T%CUR%XGARDEN(JILU)/=0..AND.GDM%TV%M%X%XROOTFRAC(JILU,JLAYER,1)>=1.0)THEN
                  GDM%TV%M%X%XDG2  (JILU,:)=GDM%TV%M%X%XDG(JILU,JLAYER,:)
                  GDM%TV%M%X%XDROOT(JILU,:)=GDM%TV%M%X%XDG(JILU,JLAYER,:)
                ENDIF
              ENDDO
            ENDDO
          ENDIF

        END IF
        !
        WHERE (T%CUR%XGARDEN(:)==0.)
          GDM%TV%M%T%CUR%XVEG(:,1)=0.
          GDM%TV%M%T%CUR%XLAI(:,1)=0.
          GDM%TV%M%T%CUR%XRSMIN(:,1)=40.
          GDM%TV%M%T%CUR%XGAMMA(:,1)=0.
          GDM%TV%M%T%CUR%XWRMAX_CF(:,1)=0.2
          GDM%TV%M%T%CUR%XRGL(:,1)=100.
          GDM%TV%M%T%CUR%XCV(:,1)=2.E-5
          GDM%TV%M%T%CUR%XZ0(:,1)=0.013
          GDM%TV%M%X%XZ0_O_Z0H(:,1)=10.
          GDM%TV%M%T%CUR%XALBNIR_VEG(:,1)=0.30
          GDM%TV%M%T%CUR%XALBVIS_VEG(:,1)=0.30
          GDM%TV%M%T%CUR%XALBUV_VEG(:,1)=0.06
          GDM%TV%M%T%CUR%XEMIS(:,1)=0.94
        ENDWHERE  
        IF (GDM%TV%O%CPHOTO/='NON') THEN
          WHERE (T%CUR%XGARDEN(:)==0.)
            GDM%TV%M%T%CUR%XGMES(:,1)=0.020
            GDM%TV%M%T%CUR%XBSLAI(:,1)=0.36
            GDM%TV%M%T%CUR%XLAIMIN(:,1)=0.3
            GDM%TV%M%T%CUR%XSEFOLD(:,1)=90*86400.
            GDM%TV%M%X%XH_TREE(:,1)=0.
            GDM%TV%M%X%XRE25(:,1)=3.6E-7
            GDM%TV%M%T%CUR%XGC(:,1)=0.00025
          END WHERE
          IF (GDM%TV%O%CPHOTO/='AGS' .AND. GDM%TV%O%CPHOTO/='LAI') THEN
            WHERE (T%CUR%XGARDEN(:)==0.) 
              GDM%TV%M%X%XDMAX(:,1)=0.1
              GDM%TV%M%T%CUR%XF2I(:,1)=0.3
            END WHERE
            IF (GDM%TV%O%CPHOTO=='NIT' .OR. GDM%TV%O%CPHOTO=='NCB') THEN
              WHERE (T%CUR%XGARDEN(:)==0.)      
                GDM%TV%M%T%CUR%XCE_NITRO(:,1)=7.68
                GDM%TV%M%T%CUR%XCF_NITRO(:,1)=-4.33
                GDM%TV%M%T%CUR%XCNA_NITRO(:,1)=1.3
              END WHERE
            ENDIF
          ENDIF
        ENDIF
        IF(GDM%TV%O%CISBA/='DIF')THEN
          DO JLAYER=1,GDM%TV%O%NGROUND_LAYER
            WHERE (T%CUR%XGARDEN(:)==0.)
              GDM%TV%M%X%XDG(:,JLAYER,1)=0.2*JLAYER
            END WHERE
          ENDDO
        ELSE
          WHERE (T%CUR%XGARDEN(:)==0.) 
            GDM%TV%M%X%XDG(:,1,1)=0.01
            GDM%TV%M%X%XDG(:,2,1)=0.04
            GDM%TV%M%X%XROOTFRAC(:,1,1)=0.
            GDM%TV%M%X%XROOTFRAC(:,2,1)=0.
          END WHERE        
          DO JLAYER=3,GDM%TV%O%NGROUND_LAYER
            WHERE (T%CUR%XGARDEN(:)==0.)
              GDM%TV%M%X%XDG(:,JLAYER,1)=0.1*(JLAYER-2)
              GDM%TV%M%X%XROOTFRAC(:,JLAYER,1)=0.
            END WHERE
          ENDDO               
          WHERE (T%CUR%XGARDEN(:)==0.) 
            GDM%TV%M%X%NWG_LAYER(:,1)=GDM%TV%O%NGROUND_LAYER
            GDM%TV%M%X%XDROOT   (:,1)=0.0
            GDM%TV%M%X%XDG2     (:,1)=GDM%TV%M%X%XDG(:,GDM%TV%O%NGROUND_LAYER-1,1)
          ENDWHERE    
        ENDIF  
        WHERE (T%CUR%XGARDEN(:)==0.) 
          GDM%TV%M%X%XD_ICE(:,1)=0.8*GDM%TV%M%X%XDG(:,2,1)
        END WHERE  
        DO JVEGTYPE=1,NVEGTYPE
          WHERE (T%CUR%XGARDEN(:)==0.)
            GDM%TV%M%X%XVEGTYPE(:,JVEGTYPE)=0.
            GDM%TV%M%X%XVEGTYPE(:,1)=1.
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
        !
         CALL INIT_VEG_PGD_n(CHI, DTCO, DST, SLT, U, &
                     GDM%TV%O%LAGRI_TO_GRASS, TOP%LCOVER, TOP%XCOVER, &
                     HPROGRAM, 'TOWN  ',ILUOUT, KI, 1, GDM%TV%O%NGROUND_LAYER, TOP%TTIME%TDATE%MONTH,    &
                     GDM%TV%M%X%XVEGTYPE,ZPATCH, ZVEGTYPE_PATCH, ISIZE_NATURE_P, IR_NATURE_P,    &
                  0.0, .FALSE., .FALSE., ZTDEEP_CLI, ZGAMMAT_CLI, &
                  GDM%TV%IP%XTDEEP, GDM%TV%IP%XGAMMAT,  .FALSE., &
                  ZTHRESHOLD, IIRRINUM, GIRRIDAY, GIRRIGATE, ZTHRESHOLDSPT, &
                  GDM%TV%O%CPHOTO, HINIT,GDM%TV%O%LTR_ML, GDM%TV%O%NNBIOMASS, &
                  PCO2, PRHOA, GDM%TV%IP%XABC, GDM%TV%IP%XPOI,  &
                  GDM%TV%M%T%CUR%XGMES, GDM%TV%M%T%CUR%XGC, GDM%TV%M%X%XDMAX, &
                  GDM%TV%IP%XANMAX, GDM%TV%IP%XFZERO, GDM%TV%IP%XEPSO, GDM%TV%IP%XGAMM, GDM%TV%IP%XQDGAMM,   &
                  GDM%TV%IP%XQDGMES, GDM%TV%IP%XT1GMES, GDM%TV%IP%XT2GMES, GDM%TV%IP%XAMAX, GDM%TV%IP%XQDAMAX, &
                  GDM%TV%IP%XT1AMAX, GDM%TV%IP%XT2AMAX,GDM%TV%IP%XAH, GDM%TV%IP%XBH,&           
                  ZTAU_WOOD, ZINCREASE, ZTURNOVER,                  &
                  KSV, HSV, CHT%SVT, CHT%CCH_NAMES, CHT%CAER_NAMES,CHT%CDSTNAMES, CHT%CSLTNAMES, &
                  CHT%CCHEM_SURF_FILE,                     &
                  ZSFDST, ZSFDSTM, ZSFSLT,                                    &
                  ZAOSIP, ZAOSIM, ZAOSJP, ZAOSJM, ZHO2IP, ZHO2IM, ZHO2JP,     &
                  ZHO2JM, ZZ0, ZZ0EFFIP, ZZ0EFFIM, ZZ0EFFJP, ZZ0EFFJM, ZZ0REL,&
                  GDM%TV%P%XCLAY, GDM%TV%P%XSAND, GDM%TV%O%CPEDOTF,      &
                  GDM%TV%IP%XCONDSAT, GDM%TV%IP%XMPOTSAT, GDM%TV%IP%XBCOEF, GDM%TV%IP%XWWILT, &
                  GDM%TV%IP%XWFC, GDM%TV%IP%XWSAT, ZWD0, ZKANISO, GDM%TV%O%CRUNOFF,  &
                  GDM%TV%IP%XTAUICE, GDM%TV%IP%XCGSAT, GDM%TV%IP%XC1SAT, &
                  GDM%TV%IP%XC2REF, GDM%TV%IP%XC3, GDM%TV%IP%XC4B, GDM%TV%IP%XACOEF, GDM%TV%IP%XPCOEF, &
                  GDM%TV%IP%XC4REF, GDM%TV%IP%XPCPS, GDM%TV%IP%XPLVTT, GDM%TV%IP%XPLSTT,        &
                  GDM%TV%O%CSCOND, GDM%TV%O%CISBA, GDM%TV%IP%XHCAPSOIL, GDM%TV%IP%XCONDDRY, &
                  GDM%TV%IP%XCONDSLD, GDM%TV%O%CCPSURF, GDM%TV%M%X%XDG, GDM%TV%M%X%XDROOT, GDM%TV%M%X%XDG2, &
                  GDM%TV%M%X%XROOTFRAC, GDM%TV%IP%XRUNOFFD, GDM%TV%IP%XDZG, GDM%TV%IP%XDZDIF,       &
                  GDM%TV%IP%XSOILWGHT, GDM%TV%M%X%NWG_LAYER, GDM%TV%O%NLAYER_HORT, &
                  GDM%TV%O%NLAYER_DUN, GDM%TV%M%X%XD_ICE,  &
                  GDM%TV%IP%XKSAT_ICE, GDM%TV%IP%XALBNIR_DRY, GDM%TV%IP%XALBVIS_DRY, GDM%TV%IP%XALBUV_DRY,   &
                  GDM%TV%IP%XALBNIR_WET, GDM%TV%IP%XALBVIS_WET, GDM%TV%IP%XALBUV_WET, GDM%TV%IP%XBSLAI_NITRO, &
                  GDM%TV%M%T%CUR%XCE_NITRO, GDM%TV%M%T%CUR%XCNA_NITRO, GDM%TV%M%T%CUR%XCF_NITRO,          &
                  ZFWTD, ZWTD               )         
!
!-------------------------------------------------------------------------------
!
IF(GDM%TV%O%CISBA=='DIF'.AND.GDM%TV%O%LSOC)THEN
  CALL ABOR1_SFX('INIT_TEB_GARDEN_PGDn: SUBGRID Soil organic matter'//&
                 ' effect (LSOC) NOT YET IMPLEMENTED FOR GARDEN')
ELSEIF (GDM%TV%O%CISBA=='3-L'.AND.GDM%TV%O%CKSAT=='EXP') THEN 
  CALL ABOR1_SFX('INIT_TEB_GARDEN_PGDn: topmodel exponential decay not implemented for garden')
ENDIF
!
IF(GDM%TV%O%CKSAT=='SGH' .AND. GDM%TV%O%CISBA/='DIF' .AND. HINIT/='PRE')THEN 
  ZF(:)=MIN(4.0/GDM%TV%M%X%XDG(:,2,1),XF_DECAY)
  CALL EXP_DECAY_SOIL_FR(GDM%TV%O%CISBA, ZF(:),GDM%TV%IP%XC1SAT(:,1),GDM%TV%IP%XC2REF(:,1),&
                         GDM%TV%M%X%XDG(:,:,1),GDM%TV%M%X%XD_ICE(:,1),GDM%TV%IP%XC4REF(:,1),&
                         GDM%TV%IP%XC3(:,:,1),GDM%TV%IP%XCONDSAT(:,:,1),GDM%TV%IP%XKSAT_ICE(:,1))
ENDIF
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('INIT_TEB_GARDEN_PGD_n',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE INIT_TEB_GARDEN_PGD_n
