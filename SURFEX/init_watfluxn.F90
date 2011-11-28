!     #############################################################
      SUBROUTINE INIT_WATFLUX_n(HPROGRAM,HINIT,                            &
                                  KI,KSV,KSW,                                &
                                  HSV,PCO2,PRHOA,                            &
                                  PZENITH,PAZIM,PSW_BANDS,PDIR_ALB,PSCA_ALB, &
                                  PEMIS,PTSRAD,                              &
                                  KYEAR, KMONTH,KDAY, PTIME,                 &
                                  HATMFILE,HATMFILETYPE,                     &
                                  HTEST                                      )  
!     #############################################################
!
!!****  *INIT_WATFLUX_n* - routine to initialize WATFLUX
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
!!      Original    01/2003
!!      B. Decharme 08/2009 : specific treatment for water/ice in the Earth System Model 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_READ_NAMELIST,  ONLY : LNAM_READ
USE MODD_SURF_ATM,       ONLY : LCPL_ESM
!
USE MODD_CSTS,           ONLY : XTT
USE MODD_WATFLUX_n,      ONLY : XCOVER, XDIR_ALB, XSCA_ALB, XEMIS, XTS, XTSTEP,   &
                                  CWAT_ALB, XOUT_TSTEP, TTIME, LSBL, CINTERPOL_TS,  &
                                  LINTERPOL_TS, XTICE, XICE_ALB  
USE MODD_DIAG_WATFLUX_n, ONLY : N2M, LSURF_BUDGET, LRAD_BUDGET, XDIAG_TSTEP, L2M_MIN_ZS, &
                                  LCOEF, LSURF_VARS, LSURF_BUDGETC, LRESET_BUDGETC  
USE MODD_CH_WATFLUX_n,   ONLY : XDEP, CCH_DRY_DEP, CSV, CCH_NAMES, &
                                  NBEQ, NSV_CHSBEG, NSV_CHSEND,  &
                                  NAEREQ, NSV_AERBEG, NSV_AEREND, CAER_NAMES,&
                                  NSV_DSTBEG, NSV_DSTEND, NDSTEQ, CDSTNAMES, &
                                  NSV_SLTBEG, NSV_SLTEND, NSLTEQ, CSLTNAMES  
USE MODD_CHS_AEROSOL,    ONLY: LVARSIGI, LVARSIGJ
USE MODD_DST_SURF,       ONLY: LVARSIG, CDSTYN, NDSTMDE, NDST_MDEBEG, LRGFIX_DST
USE MODD_SLT_SURF,       ONLY: LVARSIG_SLT, CSLTYN, NSLTMDE, NSLT_MDEBEG, LRGFIX_SLT

USE MODD_SURF_PAR,       ONLY : XUNDEF, NUNDEF
!
USE MODI_INIT_IO_SURF_n
USE MODI_DEFAULT_CH_DEP
USE MODI_DEFAULT_WATFLUX
USE MODI_DEFAULT_DIAG_WATFLUX
USE MODI_READ_DEFAULT_WATFLUX_n
USE MODI_READ_WATFLUX_CONF_n
USE MODI_READ_WATFLUX_n
USE MODI_READ_PGD_WATFLUX_n
USE MODI_DIAG_WATFLUX_INIT_n
USE MODI_END_IO_SURF_n
USE MODI_GET_LUOUT
USE MODI_READ_SURF
USE MODI_READ_WATFLUX_DATE
USE MODI_READ_NAM_PREP_WATFLUX_n
USE MODI_CH_INIT_NAMES
USE MODI_DST_INIT_NAMES
USE MODI_DST_INIT_MODES
USE MODI_SLT_INIT_NAMES
USE MODI_SLT_INIT_MODES
USE MODI_PREP_CTRL_WATFLUX
USE MODI_UPDATE_RAD_SEAWAT
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
USE MODI_WRITE_COVER_TEX_WATER
!
USE MODI_READ_WATFLUX_SBL_n
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=6),                 INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
CHARACTER(LEN=3),                 INTENT(IN)  :: HINIT     ! choice of fields to initialize
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
INTEGER           :: ILU    ! sizes of WATFLUX arrays
INTEGER           :: ILUOUT ! unit of output listing file
INTEGER           :: IRESP  ! return code
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!         Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('INIT_WATFLUX_N',0,ZHOOK_HANDLE)
CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
IF (HTEST/='OK') THEN
  CALL ABOR1_SFX('INIT_WATFLUXN: FATAL ERROR DURING ARGUMENT TRANSFER')
END IF
!
!         Other little things
!
!
IF (LNAM_READ) THEN
 !
 !
 !*       0.     Defaults
 !               --------
 !
 !        0.1. Hard defaults
 !      
 CALL DEFAULT_WATFLUX(XTSTEP,XOUT_TSTEP,CWAT_ALB,CINTERPOL_TS)
 CALL DEFAULT_CH_DEP(CCH_DRY_DEP)
 CALL DEFAULT_DIAG_WATFLUX(N2M,LSURF_BUDGET,LRAD_BUDGET,XDIAG_TSTEP, &
                             LSURF_BUDGETC,LRESET_BUDGETC              )  
 !
ENDIF
!
!        0.2. Defaults from file header
!    
CALL READ_DEFAULT_WATFLUX_n(HPROGRAM)
!
!*       1.1    Reading of configuration:
!               -------------------------
!
!
CALL READ_WATFLUX_CONF_n(HPROGRAM)
!
LINTERPOL_TS=.FALSE.
IF(LCPL_ESM)THEN       
! No TS water interpolation in Earth System Model
  CINTERPOL_TS='NONE  '
  LINTERPOL_TS=.FALSE.
ELSEIF(CINTERPOL_TS/='NONE  ')THEN
  LINTERPOL_TS=.TRUE.
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       1.     Cover fields and grid:
!               ---------------------
!* date
!
SELECT CASE (HINIT)
  CASE ('PGD')
    TTIME%TDATE%YEAR = NUNDEF
    TTIME%TDATE%MONTH= NUNDEF
    TTIME%TDATE%DAY  = NUNDEF
    TTIME%TIME       = XUNDEF

  CASE ('PRE')
    CALL PREP_CTRL_WATFLUX(N2M,LSURF_BUDGET,L2M_MIN_ZS,LRAD_BUDGET,LCOEF,LSURF_VARS,&
                             ILUOUT,LSURF_BUDGETC )  
    IF (LNAM_READ) CALL READ_NAM_PREP_WATFLUX_n(HPROGRAM)                 
    CALL READ_WATFLUX_DATE(HPROGRAM,HINIT,ILUOUT,HATMFILE,HATMFILETYPE,KYEAR,KMONTH,KDAY,PTIME,TTIME)

  CASE DEFAULT
    CALL INIT_IO_SURF_n(HPROGRAM,'WATER ','WATFLX','READ ')
    CALL READ_SURF(HPROGRAM,'DTCUR',TTIME,IRESP)
    CALL END_IO_SURF_n(HPROGRAM)
END SELECT
!
!         Initialisation for IO
!
CALL INIT_IO_SURF_n(HPROGRAM,'WATER ','WATFLX','READ ')
!
!         Reading of the fields
!
CALL READ_PGD_WATFLUX_n(HPROGRAM)
!
!-------------------------------------------------------------------------------
!
!* if only physiographic fields are to be initialized, stop here.
!
CALL WRITE_COVER_TEX_WATER
!
IF (HINIT/='ALL') THEN
  CALL END_IO_SURF_n(HPROGRAM)
  IF (LHOOK) CALL DR_HOOK('INIT_WATFLUX_N',1,ZHOOK_HANDLE)
  RETURN
END IF
!
!-------------------------------------------------------------------------------
!
!*       2.     Prognostic and cover fields:
!               ---------------------------
!
CALL READ_WATFLUX_n(HPROGRAM)
!
ILU = SIZE(XCOVER,1)
!
!
!*       3.     Specific fields when using earth system model (Ice temperature)
!               ---------------------------------------------------------------
!
IF(LCPL_ESM)THEN
  ALLOCATE(XTICE   (ILU))
  ALLOCATE(XICE_ALB(ILU))
  XTICE   (:)=XUNDEF
  XICE_ALB(:)=XUNDEF
ELSE
  ALLOCATE(XTICE   (0))
  ALLOCATE(XICE_ALB(0))
ENDIF
!
!*       4.     Albedo, emissivity and output radiative fields on open water and ice
!               --------------------------------------------------------------------
!
ALLOCATE(XDIR_ALB (ILU))
ALLOCATE(XSCA_ALB (ILU))
ALLOCATE(XEMIS    (ILU))
XDIR_ALB = 0.0
XSCA_ALB = 0.0
XEMIS    = 0.0
!
CALL UPDATE_RAD_SEAWAT(CWAT_ALB,XTS,PZENITH,XTT,XEMIS,XDIR_ALB,&
                         XSCA_ALB,PDIR_ALB,PSCA_ALB,PEMIS,PTSRAD )  
!
!-------------------------------------------------------------------------------
!
!*       5.     SBL air fields:
!               --------------
!
CALL READ_WATFLUX_SBL_n(HPROGRAM)
!
!-------------------------------------------------------------------------------
!
!*       6.     Chemistry / dust
!               ----------------
!
IF (KSV /= 0) THEN
  ALLOCATE(CSV(KSV))
  CALL CH_INIT_NAMES(ILUOUT,HSV,NBEQ,NAEREQ,CSV,&
                       NSV_CHSBEG,NSV_CHSEND, &
                       NSV_AERBEG, NSV_AEREND,&
                       LVARSIGI, LVARSIGJ)  

  IF (NBEQ > 0 ) THEN
    ALLOCATE(CCH_NAMES(NBEQ))
    CCH_NAMES(:) = CSV(NSV_CHSBEG:NSV_CHSEND)
  ELSE
    ALLOCATE(CCH_NAMES(0))
  END IF

  IF (NAEREQ > 0 ) THEN
    ALLOCATE(CAER_NAMES(NAEREQ))
    CAER_NAMES(:) = CSV(NSV_AERBEG:NSV_AEREND)
  ELSE
    ALLOCATE(CAER_NAMES(0))
  END IF

   CALL DST_INIT_NAMES(         &
          ILUOUT,                   &!I [idx] index of writing unit
          HSV,                     &!I [char] list of scalar variables
          NDSTEQ,                  &!O [nbr] number of dust related tracers
          NSV_DSTBEG,              &!O [idx] first dust related scalar variable
          NSV_DSTEND,              &!O [idx] last dust related scalar variable
          LVARSIG,                 &!O type of standard deviation (fixed or variable)
          LRGFIX_DST,              &!O type of mean radius (fixed or variable)
          CDSTYN                  &!O [char] Y/N, are dust related scalars present?
          )  

    IF (NDSTEQ >=1) THEN
   CALL DST_INIT_MODES(         &
          NDSTEQ,                   &!I [nbr] number of dust related variables in scalar list
          NSV_DSTBEG,              &!I [idx] index of first dust related variable in scalar list
          NSV_DSTEND,              &!I [idx] index of last dust related variable in scalar list
          LVARSIG,                 &!I type of standard deviation (fixed or variable)
          LRGFIX_DST,              &!I type of mean radius (fixed or variable)
          NDST_MDEBEG,             &!O [idx] index of mass for first mode in scalar list
          NDSTMDE                 &!O [nbr] number of modes to be transported
          )  
     IF(.NOT. ASSOCIATED(CDSTNAMES)) ALLOCATE (CDSTNAMES(NDSTEQ))
     CDSTNAMES(:) = HSV(NSV_DSTBEG:NSV_DSTEND)
    END IF


   CALL SLT_INIT_NAMES(         &
          ILUOUT,                   &!I [idx] index of writing unit
          HSV,                     &!I [char] list of scalar variables
          NSLTEQ,                  &!O [nbr] number of sea salt related tracers
          NSV_SLTBEG,              &!O [idx] first sea salt related scalar variable
          NSV_SLTEND,              &!O [idx] last sea salt related scalar variable
          LVARSIG_SLT,             &!O type of standard deviation (fixed or variable)
          LRGFIX_SLT,              &!O type of mean radius (fixed or variable)
          CSLTYN                  &!O [char] Y/N, are sea salt related scalars present?
          )  

    IF (NSLTEQ >=1) THEN
   CALL SLT_INIT_MODES(         &
          NSLTEQ,                   &!I [nbr] number of sea salt related variables in scalar list
          NSV_SLTBEG,              &!I [idx] index of first sea salt related variable in scalar list
          NSV_SLTEND,              &!I [idx] index of last sea salt related variable in scalar list
          LVARSIG_SLT,             &!I type of standard deviation (fixed or variable)
          LRGFIX_SLT,              &!I type of mean radius (fixed or variable)
          NSLT_MDEBEG,             &!O [idx] index of mass for first mode in scalar list
          NSLTMDE                 &!O [nbr] number of modes to be transported
          )  
     IF(.NOT. ASSOCIATED(CSLTNAMES)) ALLOCATE (CSLTNAMES(NSLTEQ))
     CSLTNAMES(:) = HSV(NSV_SLTBEG:NSV_SLTEND)
    END IF

ELSE
  ALLOCATE(CSV      (0))
ENDIF

!* depositiion scheme
!

IF (NBEQ>0 .AND. CCH_DRY_DEP=='WES89') THEN
  ALLOCATE(XDEP(ILU,NBEQ))
ELSE
  ALLOCATE(XDEP(0,0))
END IF

!
!-------------------------------------------------------------------------------
!
!*       7.     diagnostics initialization
!               --------------------------
!
CALL DIAG_WATFLUX_INIT_n(HPROGRAM,ILU,KSW)
!
!-------------------------------------------------------------------------------
!
!         End of IO
!
CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('INIT_WATFLUX_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE INIT_WATFLUX_n
