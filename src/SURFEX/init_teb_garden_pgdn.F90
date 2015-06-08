!#############################################################
SUBROUTINE INIT_TEB_GARDEN_PGD_n(HPROGRAM,HINIT, OREAD_PGD,KI, KSV, HSV, KVERSION, KBUGFIX, &
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
USE MODD_DATA_TEB_GARDEN_n, ONLY : DTGD => DATA_TEB_GARDEN
!
USE MODD_TYPE_DATE_SURF
USE MODD_TYPE_SNOW
!
USE MODD_TEB_n, ONLY : T => TEB
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS
USE MODD_TEB_VEG_n, ONLY : TVG => TEB_VEG_OPTIONS
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TGDO => TEB_GARDEN_OPTIONS
USE MODD_TEB_GARDEN_PGD_EVOL_n, ONLY : TGDPE => TEB_GARDEN_PGD_EVOL
USE MODD_TEB_GARDEN_PGD_n, ONLY : TGDP => TEB_GARDEN_PGD
USE MODD_CH_TEB_n, ONLY : CHT => CH_TEB
                                
USE MODD_DATA_COVER_PAR,  ONLY: NVEGTYPE
USE MODD_SURF_PAR,        ONLY: XUNDEF, NUNDEF

USE MODD_SGH_PAR,         ONLY: NDIMTAB, XF_DECAY
!
USE MODI_GET_LUOUT
USE MODI_ALLOCATE_TEB_GARDEN_PGD
USE MODI_READ_PGD_TEB_GARDEN_n
USE MODI_CONVERT_PATCH_GARDEN
USE MODI_INIT_FROM_DATA_GRDN_n
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
 CALL ALLOCATE_TEB_GARDEN_PGD(TGDPE, TGDP, &
                              OREAD_PGD, KI, NVEGTYPE, TGDO%NGROUND_LAYER, NDIMTAB)  
!
!
!*       2.1    Cover, soil and orographic fields:
!               ---------------------------------
!
IF (OREAD_PGD) &
 CALL READ_PGD_TEB_GARDEN_n(HPROGRAM,KVERSION,KBUGFIX)
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
IF (.NOT. TGDO%LPAR_GARDEN) THEN
  CALL CONVERT_PATCH_GARDEN(KI,IDECADE)
ELSE
 CALL INIT_FROM_DATA_GRDN_n(DTGD, &
                            IDECADE,TVG%CPHOTO,                     &
                            TGDPE%CUR%XVEG,                               &
                            TGDPE%CUR%XLAI,TGDP%XRSMIN,TGDP%XGAMMA,TGDP%XWRMAX_CF,       &
                            TGDP%XRGL,TGDP%XCV,TGDP%XDG,TGDP%XD_ICE,TGDPE%CUR%XZ0,TGDP%XZ0_O_Z0H,  &
                            TGDP%XALBNIR_VEG,TGDP%XALBVIS_VEG,            &
                            TGDP%XALBUV_VEG,TGDPE%CUR%XEMIS,                   &
                            TGDP%XVEGTYPE,TGDP%XROOTFRAC,                 &
                            TGDP%XGMES,TGDP%XBSLAI,TGDP%XLAIMIN,TGDP%XSEFOLD,TGDP%XGC,   &
                            TGDP%XDMAX, TGDP%XF2I, TGDP%LSTRESS, TGDP%XH_TREE,TGDP%XRE25,&
                            TGDP%XCE_NITRO,TGDP%XCF_NITRO,TGDP%XCNA_NITRO      )  

  IF (TVG%CISBA=='DIF') THEN
    WHERE(T%CUR%XGARDEN(:)/=0.)
      TGDP%NWG_LAYER(:)=TGDO%NGROUND_LAYER 
      TGDP%XDG2  (:)=0.0
      TGDP%XDROOT(:)=0.0
    ENDWHERE
    DO JLAYER=TGDO%NGROUND_LAYER,1,-1
      DO JILU=1,KI
        IF(T%CUR%XGARDEN(JILU)/=0..AND.TGDP%XROOTFRAC(JILU,JLAYER)>=1.0)THEN
          TGDP%XDG2  (JILU)=TGDP%XDG(JILU,JLAYER)
          TGDP%XDROOT(JILU)=TGDP%XDG(JILU,JLAYER)
        ENDIF
      ENDDO
    ENDDO
  ENDIF

END IF
!

WHERE (T%CUR%XGARDEN(:)==0.)
  TGDPE%CUR%XVEG(:)=0.
  TGDPE%CUR%XLAI(:)=0.
  TGDP%XRSMIN(:)=40.
  TGDP%XGAMMA(:)=0.
  TGDP%XWRMAX_CF(:)=0.2
  TGDP%XRGL(:)=100.
  TGDP%XCV(:)=2.E-5
  TGDPE%CUR%XZ0(:)=0.013
  TGDP%XZ0_O_Z0H(:)=10.
  TGDP%XALBNIR_VEG(:)=0.30
  TGDP%XALBVIS_VEG(:)=0.30
  TGDP%XALBUV_VEG(:)=0.06
  TGDPE%CUR%XEMIS(:)=0.94
ENDWHERE  
IF (TVG%CPHOTO/='NON') THEN
  WHERE (T%CUR%XGARDEN(:)==0.)
    TGDP%XGMES(:)=0.020
    TGDP%XBSLAI(:)=0.36
    TGDP%XLAIMIN(:)=0.3
    TGDP%XSEFOLD(:)=90*86400.
    TGDP%XH_TREE(:)=0.
    TGDP%XRE25(:)=3.6E-7
    TGDP%XGC(:)=0.00025
  END WHERE
  IF (TVG%CPHOTO/='AGS' .AND. TVG%CPHOTO/='LAI') THEN
    WHERE (T%CUR%XGARDEN(:)==0.) 
      TGDP%XDMAX(:)=0.1
      TGDP%XF2I(:)=0.3
    END WHERE
    IF (TVG%CPHOTO=='NIT' .OR. TVG%CPHOTO=='NCB') THEN
      WHERE (T%CUR%XGARDEN(:)==0.)      
        TGDP%XCE_NITRO(:)=7.68
        TGDP%XCF_NITRO(:)=-4.33
        TGDP%XCNA_NITRO(:)=1.3
      END WHERE
    ENDIF
  ENDIF
ENDIF
IF(TVG%CISBA/='DIF')THEN
  DO JLAYER=1,TGDO%NGROUND_LAYER
    WHERE (T%CUR%XGARDEN(:)==0.)
      TGDP%XDG(:,JLAYER)=0.2*JLAYER
    END WHERE
  ENDDO
ELSE
  WHERE (T%CUR%XGARDEN(:)==0.) 
    TGDP%XDG(:,1)=0.01
    TGDP%XDG(:,2)=0.04
    TGDP%XROOTFRAC(:,1)=0.
    TGDP%XROOTFRAC(:,2)=0.
  END WHERE        
  DO JLAYER=3,TGDO%NGROUND_LAYER
    WHERE (T%CUR%XGARDEN(:)==0.)
      TGDP%XDG(:,JLAYER)=0.1*(JLAYER-2)
      TGDP%XROOTFRAC(:,JLAYER)=0.
    END WHERE
  ENDDO               
  WHERE (T%CUR%XGARDEN(:)==0.) 
    TGDP%NWG_LAYER(:)=TGDO%NGROUND_LAYER
    TGDP%XDROOT   (:)=0.0
    TGDP%XDG2     (:)=TGDP%XDG(:,TGDO%NGROUND_LAYER-1)
  ENDWHERE    
ENDIF  
WHERE (T%CUR%XGARDEN(:)==0.) 
  TGDP%XD_ICE(:)=0.8*TGDP%XDG(:,2)
END WHERE  
DO JVEGTYPE=1,NVEGTYPE
  WHERE (T%CUR%XGARDEN(:)==0.)
    TGDP%XVEGTYPE(:,JVEGTYPE)=0.
    TGDP%XVEGTYPE(:,1)=1.
  END WHERE
ENDDO
!
 CALL INIT_VEG_PGD_GARDEN_n(HPROGRAM, ILUOUT, KI, TGDO%NGROUND_LAYER, TOP%TTIME%TDATE%MONTH,    &
                        TGDP%XVEGTYPE, TGDP%XTDEEP, TGDP%XGAMMAT, TVG%CPHOTO, HINIT, TVG%LTR_ML, TVG%CRUNOFF,  &
                        TVG%NNBIOMASS, PCO2, PRHOA, TGDP%XABC, TGDP%XPOI,                         &
                        TGDP%XGMES, TGDP%XGC, TGDP%XDMAX, TGDP%XANMAX, TGDP%XFZERO, TGDP%XEPSO, TGDP%XGAMM, TGDP%XQDGAMM,   &
                        TGDP%XQDGMES, TGDP%XT1GMES, TGDP%XT2GMES, TGDP%XAMAX, TGDP%XQDAMAX, TGDP%XT1AMAX, TGDP%XT2AMAX,&
                        TGDP%XAH, TGDP%XBH,                                                   &
                        KSV, HSV, CHT%NBEQ, CHT%CSV, CHT%NAEREQ, CHT%NSV_CHSBEG, CHT%NSV_CHSEND,        &
                        CHT%NSV_AERBEG, CHT%NSV_AEREND, CHT%CCH_NAMES, CHT%CAER_NAMES, CHT%NDSTEQ,      &
                        CHT%NSV_DSTBEG, CHT%NSV_DSTEND, CHT%NSLTEQ, CHT%NSV_SLTBEG, CHT%NSV_SLTEND,     &
                        CHT%CDSTNAMES, CHT%CSLTNAMES, CHT%CCHEM_SURF_FILE,                      &
                        TGDP%XCLAY, TGDP%XSAND, TVG%CPEDOTF,                                      &
                        TGDP%XCONDSAT, TGDP%XMPOTSAT, TGDP%XBCOEF, TGDP%XWWILT, TGDP%XWFC, TGDP%XWSAT,            &
                        TGDP%XTAUICE, TGDP%XCGSAT, TGDP%XC1SAT, TGDP%XC2REF, TGDP%XC3, TGDP%XC4B, TGDP%XACOEF, TGDP%XPCOEF, &
                        TGDP%XC4REF, TGDP%XPCPS, TGDP%XPLVTT, TGDP%XPLSTT,                              &
                        TVG%CSCOND, TVG%CISBA, TGDP%XHCAPSOIL, TGDP%XCONDDRY, TGDP%XCONDSLD, TVG%CCPSURF,      &
                        TGDP%XDG, TGDP%XDROOT, TGDP%XDG2, TGDP%XROOTFRAC, TGDP%XRUNOFFD, TGDP%XDZG, TGDP%XDZDIF,       &
                        TGDP%XSOILWGHT, TGDP%NWG_LAYER, TGDO%NLAYER_HORT, TGDO%NLAYER_DUN, TGDP%XD_ICE,      &
                        TGDP%XKSAT_ICE, TGDP%XALBNIR_DRY, TGDP%XALBVIS_DRY, TGDP%XALBUV_DRY,            &
                        TGDP%XALBNIR_WET, TGDP%XALBVIS_WET, TGDP%XALBUV_WET, TGDP%XBSLAI_NITRO,         &
                        TGDP%XCE_NITRO, TGDP%XCNA_NITRO, TGDP%XCF_NITRO                            )
!
!-------------------------------------------------------------------------------
!
IF(TVG%CISBA=='DIF'.AND.TVG%LSOC)THEN
  CALL ABOR1_SFX('INIT_TEB_GARDEN_PGDn: SUBGRID Soil organic matter'//&
                 ' effect (LSOC) NOT YET IMPLEMENTED FOR GARDEN')
ELSEIF (TVG%CISBA=='3-L'.AND.TVG%CKSAT=='EXP') THEN 
  CALL ABOR1_SFX('INIT_TEB_GARDEN_PGDn: topmodel exponential decay not implemented for garden')
ENDIF
!
IF(TVG%CKSAT=='SGH' .AND. TVG%CISBA/='DIF' .AND. HINIT/='PRE')THEN 
  ZF(:)=MIN(4.0/TGDP%XDG(:,2),XF_DECAY)
  CALL EXP_DECAY_SOIL_FR(TVG%CISBA, ZF(:),TGDP%XC1SAT(:),TGDP%XC2REF(:),TGDP%XDG(:,:),TGDP%XD_ICE(:),&
                         TGDP%XC4REF(:),TGDP%XC3(:,:),TGDP%XCONDSAT(:,:),TGDP%XKSAT_ICE(:))
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
