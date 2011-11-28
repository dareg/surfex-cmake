!#############################################################
SUBROUTINE INIT_ISBA_n    (HPROGRAM,HINIT,OLAND_USE,                    &
                             KI,KSV,KSW,                                &
                             HSV,PCO2,PRHOA,                            &
                             PZENITH,PAZIM,PSW_BANDS,PDIR_ALB,PSCA_ALB, &
                             PEMIS,PTSRAD,                              &
                             KYEAR, KMONTH,KDAY, PTIME,                 &
                             HATMFILE,HATMFILETYPE,                     &
                             HTEST                                      )  
!#############################################################
!
!!****  *INIT_ISBA_n* - routine to initialize ISBA
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
!!	V. Masson   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      Modified by P. Le Moigne (11/2004): miscellaneous diagnostics
!!      Modified by P. Le Moigne (06/2006): seeding and irrigation    
!!      Modified by B. Decharme    (2008) : SGH and Flooding scheme
!!      Modified by B. Decharme  (01/2009): optional deep soil temperature as in Arpege
!!      Modified by R. Hamdi     (01/2009): Cp and L
!!      Modified by B. Decharme  (06/2009): read topographic index statistics
!!      Modified by P. Le Moigne (01/2009): Beljaars sso
!!      Modified by B. Decharme  (08/2009): Active Trip coupling variable if Earth System Model
!!      A.L. Gibelin   04/09 : change BSLAI_NITRO initialisation
!!      A.L. Gibelin   04/09 : modifications for CENTURY model 
!!      A.L. Gibelin   06/09 : soil carbon initialisation
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_TYPE_DATE_SURF
USE MODD_TYPE_SNOW
!
USE MODD_ISBA_n,   ONLY : CROUGH,CISBA,CPHOTO,CRUNOFF,CALBEDO,CSCOND,           &
                            CC1DRY, CSOILFRZ, CDIFSFCOND, CSNOWRES, CRESPSL,    &
                            NNLITTER, NNLITTLEVS, NNSOILCARB, NPATCH,           &
                            TSNOW, TTIME, XTSTEP, XOUT_TSTEP,                   &
                            LTRIP, LFLOOD, LGLACIER, LVEGUPD, LCANOPY_DRAG,     &
                            CCPSURF, CHORT, XCGMAX, CKSAT, CTOPREG, CRAIN
!
USE MODD_CH_ISBA_n,      ONLY : LCH_BIO_FLUX, CCH_DRY_DEP  

USE MODD_DIAG_ISBA_n,    ONLY : N2M, LSURF_BUDGET, LRAD_BUDGET,          &
                                  XDIAG_TSTEP, LPGD,  L2M_MIN_ZS, LCOEF, &
                                  LSURF_VARS, LPATCH_BUDGET  
USE MODD_DIAG_EVAP_ISBA_n, ONLY : LSURF_EVAP_BUDGET, LSURF_BUDGETC, LRESET_BUDGETC
USE MODD_DIAG_MISC_ISBA_n, ONLY : LSURF_MISC_BUDGET, LSURF_DIAG_ALBEDO,     &
                                    LWOOD_SPIN, LSOILCARB_SPIN  
USE MODD_SURF_PAR,       ONLY : XUNDEF, NUNDEF
USE MODD_AGRI,           ONLY : LAGRIP
!
USE MODD_READ_NAMELIST,  ONLY : LNAM_READ
!
USE MODI_INIT_IO_SURF_n
!
USE MODI_GET_LUOUT
USE MODI_ABOR1_SFX
USE MODI_DEFAULT_ISBA
USE MODI_DEFAULT_CH_DEP
USE MODI_DEFAULT_CH_BIO_FLUX
USE MODI_DEFAULT_DIAG_ISBA
USE MODI_READ_DEFAULT_ISBA_n
USE MODI_READ_ISBA_CONF_n
USE MODI_READ_PREP_ISBA_SNOW
USE MODI_READ_PREP_ISBA_CARBON
USE MODI_READ_SURF
USE MODI_PREP_CTRL_ISBA
USE MODI_READ_ISBA_DATE
USE MODI_READ_PGD_ISBA_n
USE MODI_COMPUTE_ISBA_PARAMETERS
!
USE MODI_END_IO_SURF_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_READ_NAM_PREP_ISBA_n
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=6),                 INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
CHARACTER(LEN=3),                 INTENT(IN)  :: HINIT     ! choice of fields to initialize
LOGICAL,                          INTENT(IN)  :: OLAND_USE !
INTEGER,                          INTENT(IN)  :: KI        ! number of points
INTEGER,                          INTENT(IN)  :: KSV       ! number of scalars
INTEGER,                          INTENT(IN)  :: KSW       ! number of short-wave spectral bands
CHARACTER(LEN=6), DIMENSION(KSV), INTENT(IN)  :: HSV       ! name of all scalar variables
REAL,             DIMENSION(KI),  INTENT(IN)  :: PCO2      ! CO2 concentration (kg/m3)
REAL,             DIMENSION(KI),  INTENT(IN)  :: PRHOA     ! air density
REAL,             DIMENSION(KI),  INTENT(IN)  :: PZENITH   ! solar zenithal angle
REAL,             DIMENSION(KI),  INTENT(IN)  :: PAZIM     ! solar azimuthal angle (rad from N, clock)
REAL,             DIMENSION(KSW), INTENT(IN)  :: PSW_BANDS ! middle wavelength of each band
REAL,             DIMENSION(KI,KSW),INTENT(OUT) :: PDIR_ALB  ! direct albedo for each band
REAL,             DIMENSION(KI,KSW),INTENT(OUT) :: PSCA_ALB  ! diffuse albedo for each band
REAL,             DIMENSION(KI),  INTENT(OUT) :: PEMIS     ! emissivity
REAL,             DIMENSION(KI),  INTENT(OUT) :: PTSRAD    ! radiative temperature
INTEGER,                          INTENT(IN)  :: KYEAR     ! current year (UTC)
INTEGER,                          INTENT(IN)  :: KMONTH    ! current month (UTC)
INTEGER,                          INTENT(IN)  :: KDAY      ! current day (UTC)
REAL,                             INTENT(IN)  :: PTIME     ! current time since
                                                          !  midnight (UTC, s)
!
CHARACTER(LEN=28),                INTENT(IN)  :: HATMFILE    ! atmospheric file name
CHARACTER(LEN=6),                 INTENT(IN)  :: HATMFILETYPE! atmospheric file type
CHARACTER(LEN=2),                 INTENT(IN)  :: HTEST       ! must be equal to 'OK'
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: ILU      ! sizes of ISBA arrays
INTEGER           :: JILU     ! loop increment
INTEGER           :: ILUOUT   ! unit of output listing file
INTEGER           :: ICH      ! unit of input chemistry file
!
INTEGER           :: IDECADE  ! decade of simulation
!
INTEGER           :: IVERSION       ! surface version
!
INTEGER :: JLAYER  ! loop counter on layers
INTEGER :: JPATCH  ! loop counter on tiles
INTEGER :: IRESP   ! return code
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZWG1 ! work array for surface water content
REAL, DIMENSION(:,:), ALLOCATABLE :: ZTG1 ! work array for surface temperature
REAL, DIMENSION(SIZE(PCO2))       :: ZCO2  ! CO2 concentration  (kg/kg)
REAL, DIMENSION(SIZE(PTSRAD))     :: ZTSRAD_NAT !radiative temperature
!
REAL, DIMENSION(:), ALLOCATABLE   :: ZM
REAL, DIMENSION(:,:), ALLOCATABLE :: ZF
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!               Initialisation for IO
!
!
IF (LHOOK) CALL DR_HOOK('INIT_ISBA_N',0,ZHOOK_HANDLE)
CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('INIT_ISBAN: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF
!
!               Other little things
!
!*       0.     Defaults
!               --------
!
LSURF_DIAG_ALBEDO = .FALSE.
!
IF (LNAM_READ) THEN
 !
 !        0.1. Hard defaults
 !      
 CALL DEFAULT_ISBA(XTSTEP, XOUT_TSTEP,                           &
                     CROUGH,CRUNOFF,CALBEDO,CSCOND,              &
                     CC1DRY, CSOILFRZ, CDIFSFCOND, CSNOWRES,     &
                     CCPSURF, XCGMAX, CKSAT, CTOPREG,            &
                     CRAIN, CHORT, LFLOOD, LTRIP , LGLACIER,     &
                     LCANOPY_DRAG, LVEGUPD                       )  
 !                  
 CALL DEFAULT_CH_DEP(CCH_DRY_DEP)
 CALL DEFAULT_CH_BIO_FLUX(LCH_BIO_FLUX)                  
 CALL DEFAULT_DIAG_ISBA(N2M,LSURF_BUDGET,LSURF_EVAP_BUDGET,        &
                          LSURF_MISC_BUDGET,LSURF_BUDGETC,LPGD,      &
                          LRAD_BUDGET, LWOOD_SPIN, LSOILCARB_SPIN,   &
                          XDIAG_TSTEP,LRESET_BUDGETC,LPATCH_BUDGET)  
 !
ENDIF
!
!        0.2. Defaults from file header
!    
CALL READ_DEFAULT_ISBA_n(HPROGRAM)
!
CALL READ_ISBA_CONF_n(HPROGRAM)
!
!
!*       1.     Reading of configuration:
!               -------------------------
!
!* initialization of snow scheme
!
IF (HINIT=='PRE') THEN 
  CALL READ_PREP_ISBA_SNOW(HPROGRAM,TSNOW%SCHEME,TSNOW%NLAYER) 
!
!* initialization of soil carbon scheme
!
  CALL READ_PREP_ISBA_CARBON(HPROGRAM,CRESPSL)
!
  IF (CRESPSL=='CNT') THEN
    NNLITTER = 2
    NNLITTLEVS = 2
    NNSOILCARB = 3
  ELSE
    NNLITTER = 0
    NNLITTLEVS = 0
    NNSOILCARB = 0
  ENDIF

ELSEIF (HINIT=='ALL') THEN

  CALL INIT_IO_SURF_n(HPROGRAM,'NATURE','ISBA  ','READ ')
  CALL READ_SURF(HPROGRAM,'VERSION',IVERSION,IRESP)
  IF (IVERSION<6) THEN
    CRESPSL='DEF'
  ELSE  
    CALL READ_SURF(HPROGRAM,'RESPSL',CRESPSL,IRESP)
    CALL READ_SURF(HPROGRAM,'NLITTER',NNLITTER,IRESP)
    CALL READ_SURF(HPROGRAM,'NLITTLEVS',NNLITTLEVS,IRESP)
    CALL READ_SURF(HPROGRAM,'NSOILCARB',NNSOILCARB,IRESP)
  ENDIF
  CALL END_IO_SURF_n(HPROGRAM)

ENDIF
!
!-------------------------------------------------------------------------------
!
!*       2.     Physiographic fields
!               --------------------
!
!
!* date
!
SELECT CASE (HINIT)
  CASE ('PGD')
    TTIME%TDATE%YEAR = NUNDEF
    TTIME%TDATE%MONTH= NUNDEF
    TTIME%TDATE%DAY  = NUNDEF
    TTIME%TIME       = XUNDEF

  CASE ('PRE')
    CALL PREP_CTRL_ISBA(N2M,LSURF_BUDGET,L2M_MIN_ZS,LRAD_BUDGET,LCOEF,LSURF_VARS,&
                          LSURF_EVAP_BUDGET,LSURF_MISC_BUDGET,LSURF_BUDGETC,     &
                          LPATCH_BUDGET,LWOOD_SPIN,LSOILCARB_SPIN,ILUOUT )    
    IF (LNAM_READ) CALL READ_NAM_PREP_ISBA_n(HPROGRAM)                        
    CALL READ_ISBA_DATE(HPROGRAM,HINIT,ILUOUT,HATMFILE,HATMFILETYPE,KYEAR,KMONTH,KDAY,PTIME,TTIME)

  CASE DEFAULT
    CALL INIT_IO_SURF_n(HPROGRAM,'NATURE','ISBA  ','READ ')
    CALL READ_SURF(HPROGRAM,'DTCUR',TTIME,IRESP)
    CALL END_IO_SURF_n(HPROGRAM)
END SELECT
!
! initialization for I/O
!
CALL INIT_IO_SURF_n(HPROGRAM,'NATURE','ISBA  ','READ ')
!
!
!*       2.1    Cover, soil and orographic fields:
!               ---------------------------------
!
CALL READ_PGD_ISBA_n(HPROGRAM,OLAND_USE)
IF ( CPHOTO/='NON' .AND. NPATCH/=12) THEN
  CALL ABOR1_SFX('INIT_ISBAN: INCONSISTENCY BETWEEN CPHOTO AND NPATCH')
ENDIF
IF ( CPHOTO/='LAI' .AND. CPHOTO/='LST' .AND. CPHOTO/='NIT' .AND. CPHOTO/='NCB' .AND. LAGRIP) THEN
  CALL ABOR1_SFX('INIT_ISBAN: INCONSISTENCY BETWEEN CPHOTO AND LAGRIP')
ENDIF
IF ( CPHOTO/='NCB' .AND. CRESPSL=='CNT') THEN
  CALL ABOR1_SFX('INIT_ISBAN: INCONSISTENCY BETWEEN CPHOTO AND CRESPSL')
ENDIF
IF (HINIT=='PRE' .AND. TSNOW%SCHEME.NE.'3-L' .AND. TSNOW%SCHEME.NE.'CRO' .AND. CISBA=='DIF') THEN
    CALL ABOR1_SFX("INIT_ISBAN: WITH CISBA = DIF, CSNOW MUST BE 3-L OR CRO")
ENDIF
!
CALL END_IO_SURF_n(HPROGRAM)
!
IF (OLAND_USE) THEN
  IF (LHOOK) CALL DR_HOOK('INIT_ISBA_N',1,ZHOOK_HANDLE)
  RETURN
END IF
!
CALL COMPUTE_ISBA_PARAMETERS(HPROGRAM,HINIT,OLAND_USE,                  &
                             KI,KSV,KSW,                                &
                             HSV,PCO2,PRHOA,                            &
                             PZENITH,PAZIM,PSW_BANDS,PDIR_ALB,PSCA_ALB, &
                             PEMIS,PTSRAD,                              &
                             KYEAR, KMONTH,KDAY, PTIME,                 &
                             HATMFILE,HATMFILETYPE,                     &
                             HTEST                                )
!  
IF (LHOOK) CALL DR_HOOK('INIT_ISBA_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE INIT_ISBA_n
