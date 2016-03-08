!#############################################################
SUBROUTINE INIT_TEB_GARDEN_PGD_n (DTCO, U, OCH_BIO_FLUX, TG, PGARDEN, TOP, &
                                 VO, VMX, VMT, VP, VIP, DT, GB,  &
                                 HPROGRAM, HINIT, OREAD_PGD,KI, &
                                 KVERSION, KBUGFIX, PCO2, PRHOA)
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
!&
!
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
!!  11/2013 (B. Decharme) No exp profile with DIF
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!      ------------
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

USE MODD_DATA_COVER_PAR,  ONLY: NVEGTYPE
USE MODD_SURF_PAR,ONLY: XUNDEF, NUNDEF

USE MODD_SGH_PAR, ONLY: XF_DECAY
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
!      -------------------------
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
LOGICAL, INTENT(IN) :: OCH_BIO_FLUX
TYPE(GRID_t), INTENT(INOUT) :: TG
REAL, DIMENSION(:), INTENT(IN) :: PGARDEN
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: VO
TYPE(ISBA_PARAM_FIX_t), INTENT(INOUT) :: VMX
TYPE(ISBA_PARAM_TIME_t), INTENT(INOUT) :: VMT
TYPE(ISBA_PGD_t), INTENT(INOUT) :: VP
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: VIP
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DT
TYPE(GR_BIOG_t), INTENT(INOUT) :: GB
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=3),   INTENT(IN)  :: HINIT     ! choice of fields to initialize
LOGICAL,    INTENT(IN)  :: OREAD_PGD ! flag to read PGD fields in the file
INTEGER,    INTENT(IN)  :: KI! number of points
INTEGER,    INTENT(IN)  :: KVERSION  ! version number of the file being read
INTEGER,    INTENT(IN)  :: KBUGFIX
REAL,     DIMENSION(KI),    INTENT(IN)  :: PCO2! CO2 concentration (kg/m3)
REAL,     DIMENSION(KI),    INTENT(IN)  :: PRHOA       ! air density
!
!
!
!*       0.2   Declarations of local variables
!      -------------------------------
!
TYPE(SSO_t) :: YGDSS
TYPE(AGRI_t) :: YAG
!
INTEGER   :: JILU     ! loop increment
INTEGER   :: ILUOUT   ! unit of output listing file
!
INTEGER   :: IDECADE  ! decade of simulation
!
INTEGER :: JVEGTYPE, JLAYER  ! loop counter on vegtypes
!
REAL, DIMENSION(KI)       :: ZF
REAL, DIMENSION(KI)       :: ZWORK
!
REAL, DIMENSION(0) :: ZTDEEP_CLI, ZGAMMAT_CLI, ZTHRESHOLD
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!       Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('INIT_TEB_GARDEN_PGD_n',0,ZHOOK_HANDLE)
!
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
CALL SSO_INIT(YGDSS)
!-------------------------------------------------------------------------------
!
!*       2.     Physiographic fields
!       --------------------
!
!* allocation of urban green area variables
!
 CALL ALLOCATE_TEB_GARDEN_PGD(VMX, VMT, VP, VIP, OREAD_PGD, KI, NVEGTYPE, VO%NGROUND_LAYER )  
!
!
!*       2.1    Cover, soil and orographic fields:
!       ---------------------------------
!
IF (OREAD_PGD) &
 CALL READ_PGD_TEB_GARDEN_n(OCH_BIO_FLUX, DTCO, DT, GB, U, &
                            VO, VP, TG%NDIM, TOP, HPROGRAM,KVERSION,KBUGFIX)
!
!
!*       2.3    Physiographic data fields from land cover:
!       -----------------------------------------
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
                        .FALSE.,'GRD', .FALSE., VMX, VMT, PSOILGRID=VO%XSOILGRID  )   

ELSE

  CALL INIT_FROM_DATA_GRDN_n(DT, IDECADE, VO%CPHOTO, .FALSE., VMX=VMX, VMT=VMT )

  IF (VO%CISBA=='DIF') THEN
    WHERE(PGARDEN(:)/=0.)
      VMX%NWG_LAYER(:,1)=VO%NGROUND_LAYER 
      VMX%XDG2  (:,1)=0.0
      VMX%XDROOT(:,1)=0.0
    ENDWHERE
    DO JLAYER=VO%NGROUND_LAYER,1,-1
      DO JILU=1,KI
        IF(PGARDEN(JILU)/=0..AND.VMX%XROOTFRAC(JILU,JLAYER,1)>=1.0)THEN
          VMX%XDG2  (JILU,:)=VMX%XDG(JILU,JLAYER,:)
          VMX%XDROOT(JILU,:)=VMX%XDG(JILU,JLAYER,:)
        ENDIF
      ENDDO
    ENDDO
  ENDIF

END IF
!

WHERE (PGARDEN(:)==0.)
  VMT%XVEG(:,1)=0.
  VMT%XLAI(:,1)=0.
  VMT%XRSMIN(:,1)=40.
  VMT%XGAMMA(:,1)=0.
  VMT%XWRMAX_CF(:,1)=0.2
  VMT%XRGL(:,1)=100.
  VMT%XCV(:,1)=2.E-5
  VMT%XZ0(:,1)=0.013
  VMX%XZ0_O_Z0H(:,1)=10.
  VMT%XALBNIR_VEG(:,1)=0.30
  VMT%XALBVIS_VEG(:,1)=0.30
  VMT%XALBUV_VEG(:,1)=0.06
  VMT%XEMIS(:,1)=0.94
ENDWHERE  
IF (VO%CPHOTO/='NON') THEN
  WHERE (PGARDEN(:)==0.)
    VMT%XGMES(:,1)=0.020
    VMT%XBSLAI(:,1)=0.36
    VMT%XLAIMIN(:,1)=0.3
    VMT%XSEFOLD(:,1)=90*86400.
    VMX%XH_TREE(:,1)=0.
    VMX%XRE25(:,1)=3.6E-7
    VMT%XGC(:,1)=0.00025
  END WHERE
  IF (VO%CPHOTO/='AGS' .AND. VO%CPHOTO/='LAI') THEN
    WHERE (PGARDEN(:)==0.) 
      VMX%XDMAX(:,1)=0.1
      VMT%XF2I(:,1)=0.3
    END WHERE
    IF (VO%CPHOTO=='NIT' .OR. VO%CPHOTO=='NCB') THEN
      WHERE (PGARDEN(:)==0.)      
        VMT%XCE_NITRO(:,1)=7.68
        VMT%XCF_NITRO(:,1)=-4.33
        VMT%XCNA_NITRO(:,1)=1.3
      END WHERE
    ENDIF
  ENDIF
ENDIF
IF(VO%CISBA/='DIF')THEN
  DO JLAYER=1,VO%NGROUND_LAYER
    WHERE (PGARDEN(:)==0.)
      VMX%XDG(:,JLAYER,1)=0.2*JLAYER
    END WHERE
  ENDDO
ELSE
  WHERE (PGARDEN(:)==0.) 
    VMX%XDG(:,1,1)=0.01
    VMX%XDG(:,2,1)=0.04
    VMX%XROOTFRAC(:,1,1)=0.
    VMX%XROOTFRAC(:,2,1)=0.
  END WHERE
  DO JLAYER=3,VO%NGROUND_LAYER
    WHERE (PGARDEN(:)==0.)
      VMX%XDG(:,JLAYER,1)=0.1*(JLAYER-2)
      VMX%XROOTFRAC(:,JLAYER,1)=0.
    END WHERE
  ENDDO       
  WHERE (PGARDEN(:)==0.) 
    VMX%NWG_LAYER(:,1)=VO%NGROUND_LAYER
    VMX%XDROOT   (:,1)=0.0
    VMX%XDG2     (:,1)=VMX%XDG(:,VO%NGROUND_LAYER-1,1)
  ENDWHERE    
ENDIF  
WHERE (PGARDEN(:)==0.) 
  VMX%XD_ICE(:,1)=0.8*VMX%XDG(:,2,1)
END WHERE  
DO JVEGTYPE=1,NVEGTYPE
  WHERE (PGARDEN(:)==0.)
    VMX%XVEGTYPE(:,JVEGTYPE)=0.
    VMX%XVEGTYPE(:,1)=1.
  END WHERE
ENDDO
!
ALLOCATE(YGDSS%XAOSIP(0))
!
 CALL INIT_VEG_PGD_n(YGDSS, DTCO, U, VO, VIP, VP, VMX, VMT, YAG, &
     HPROGRAM, 'TOWN  ',ILUOUT, KI, TOP%TTIME%TDATE%MONTH,    &
     .FALSE., .FALSE., ZTDEEP_CLI, ZGAMMAT_CLI, &
     .FALSE., ZTHRESHOLD, HINIT, PCO2, PRHOA ) 
!
!-------------------------------------------------------------------------------
!
IF(VO%CISBA=='DIF'.AND.VO%LSOC)THEN
  CALL ABOR1_SFX('INIT_TEB_GARDEN_PGDn: SUBGRID Soil organic matter'//&
 ' effect (LSOC) NOT YET IMPLEMENTED FOR GARDEN')
ELSEIF (VO%CISBA=='3-L'.AND.VO%CKSAT=='EXP') THEN 
  CALL ABOR1_SFX('INIT_TEB_GARDEN_PGDn: topmodel exponential decay not implemented for garden')
ENDIF
!
IF(VO%CKSAT=='SGH' .AND. VO%CISBA/='DIF' .AND. HINIT/='PRE')THEN 
  ZF(:)=MIN(4.0/VMX%XDG(:,2,1),XF_DECAY)
  CALL EXP_DECAY_SOIL_FR(VO%CISBA, SPREAD(ZF,2,1),VIP,VMX)
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
