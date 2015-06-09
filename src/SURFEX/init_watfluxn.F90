!     #############################################################
      SUBROUTINE INIT_WATFLUX_n(HPROGRAM,HINIT,                            &
                                  KI,KSV,KSW,                                &
                                  HSV,PCO2,PRHOA,                            &
                                  PZENITH,PAZIM,PSW_BANDS,PDIR_ALB,PSCA_ALB, &
                                  PEMIS,PTSRAD,PTSURF,                       &
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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2003
!!      B. Decharme 08/2009 : specific treatment for water/ice in the Earth System Model 
!!      B. Decharme 07/2011 : read pgd+prep 
!!       B.Decharme 04/2013 new coupling variables
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_IO_BUFF_n, ONLY : IOB => IO_BUFF
!
USE MODD_DIAG_SURF_ATM_n, ONLY : DGU => DIAG_SURF_ATM
!
USE MODD_READ_NAMELIST,  ONLY : LNAM_READ
USE MODD_SFX_OASIS,      ONLY : LCPL_SEA, LCPL_SEAICE
!
USE MODD_CSTS,           ONLY : XTT
USE MODD_WATFLUX_n, ONLY : W => WATFLUX
USE MODD_DIAG_WATFLUX_n, ONLY : DGW => DIAG_WATFLUX
USE MODD_CH_WATFLUX_n, ONLY : CHW => CH_WATFLUX
USE MODD_CHS_AEROSOL,    ONLY: LVARSIGI, LVARSIGJ
USE MODD_DST_SURF,       ONLY: LVARSIG_DST, NDSTMDE, NDST_MDEBEG, LRGFIX_DST
USE MODD_SLT_SURF,       ONLY: LVARSIG_SLT, NSLTMDE, NSLT_MDEBEG, LRGFIX_SLT
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
USE MODI_INIT_CHEMICAL_n
USE MODI_PREP_CTRL_WATFLUX
USE MODI_UPDATE_RAD_WATER
!
USE MODI_READ_WATFLUX_SBL_n
USE MODI_SET_SURFEX_FILEIN
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
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
REAL,             DIMENSION(KI),  INTENT(OUT) :: PTSURF    ! surface effective temperature         (K)
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
PDIR_ALB = XUNDEF
PSCA_ALB = XUNDEF
PEMIS    = XUNDEF
PTSRAD   = XUNDEF
PTSURF   = XUNDEF
!
IF (LNAM_READ) THEN
 !
 !
 !*       0.     Defaults
 !               --------
 !
 !        0.1. Hard defaults
 !      
 CALL DEFAULT_WATFLUX(W%XTSTEP,W%XOUT_TSTEP,W%CWAT_ALB,W%CINTERPOL_TS)
 CALL DEFAULT_CH_DEP(CHW%CCH_DRY_DEP)
 CALL DEFAULT_DIAG_WATFLUX(DGW%N2M,DGW%LSURF_BUDGET,DGW%L2M_MIN_ZS,DGW%LRAD_BUDGET,DGW%LCOEF,DGW%LSURF_VARS, &
                           DGW%LSURF_BUDGETC,DGW%LRESET_BUDGETC,DGW%XDIAG_TSTEP        )  
 !
ENDIF
!
!        0.2. Defaults from file header
!    
 CALL READ_DEFAULT_WATFLUX_n(CHW, DGW, W, &
                             HPROGRAM)
!
!*       1.1    Reading of configuration:
!               -------------------------
!
!
 CALL READ_WATFLUX_CONF_n(CHW, DGW, W, &
                          HPROGRAM)
!
W%LINTERPOL_TS=.FALSE.
IF(LCPL_SEA)THEN       
! No TS water interpolation in Earth System Model
  W%CINTERPOL_TS='NONE  '
  W%LINTERPOL_TS=.FALSE.
ELSEIF(W%CINTERPOL_TS/='NONE  ')THEN
  W%LINTERPOL_TS=.TRUE.
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
    W%TTIME%TDATE%YEAR = NUNDEF
    W%TTIME%TDATE%MONTH= NUNDEF
    W%TTIME%TDATE%DAY  = NUNDEF
    W%TTIME%TIME       = XUNDEF

  CASE ('PRE')
    CALL PREP_CTRL_WATFLUX(DGW%N2M,DGW%LSURF_BUDGET,DGW%L2M_MIN_ZS,DGW%LRAD_BUDGET,DGW%LCOEF,DGW%LSURF_VARS,&
                             ILUOUT,DGW%LSURF_BUDGETC )  
    IF (LNAM_READ) CALL READ_NAM_PREP_WATFLUX_n(HPROGRAM)                 
    CALL READ_WATFLUX_DATE(IOB, &
                           HPROGRAM,HINIT,ILUOUT,HATMFILE,HATMFILETYPE,KYEAR,KMONTH,KDAY,PTIME,W%TTIME)

  CASE DEFAULT
    CALL INIT_IO_SURF_n(HPROGRAM,'WATER ','WATFLX','READ ')
    CALL READ_SURF(HPROGRAM,'DTCUR',W%TTIME,IRESP)
    CALL END_IO_SURF_n(HPROGRAM)
END SELECT
!
!-----------------------------------------------------------------------------------------------------
! READ PGD FILE
!-----------------------------------------------------------------------------------------------------
!
!        1.3. Schemes used
!
!         Initialisation for IO
!
 CALL SET_SURFEX_FILEIN(HPROGRAM,'PGD ') ! change input file name to pgd name
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
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
 CALL SET_SURFEX_FILEIN(HPROGRAM,'PREP') ! restore input file name
!
!-----------------------------------------------------------------------------------------------------
! END READ PGD FILE
!-----------------------------------------------------------------------------------------------------
!
IF (HINIT/='ALL' .AND. HINIT/='SOD') THEN
  IF (LHOOK) CALL DR_HOOK('INIT_WATFLUX_N',1,ZHOOK_HANDLE)
  RETURN
END IF
!
!-------------------------------------------------------------------------------
!
!         Initialisation for IO
!
 CALL INIT_IO_SURF_n(HPROGRAM,'WATER ','WATFLX','READ ')
!
!
!*       2.     Prognostic and cover fields:
!               ---------------------------
!
 CALL READ_WATFLUX_n(HPROGRAM)
!
IF (HINIT/='ALL') THEN
  CALL END_IO_SURF_n(HPROGRAM)
  IF (LHOOK) CALL DR_HOOK('INIT_WATFLUX_N',1,ZHOOK_HANDLE)
  RETURN
END IF
!
ILU = SIZE(W%XCOVER,1)
!
!
!*       3.     Specific fields when using earth system model (Ice temperature)
!               ---------------------------------------------------------------
!
IF(LCPL_SEAICE)THEN
  ALLOCATE(W%XTICE   (ILU))
  ALLOCATE(W%XICE_ALB(ILU))
  W%XTICE   (:)=XUNDEF
  W%XICE_ALB(:)=XUNDEF
ELSE
  ALLOCATE(W%XTICE   (0))
  ALLOCATE(W%XICE_ALB(0))
ENDIF
!
!*       4.     Albedo, emissivity and temperature fields on open water and ice
!               ---------------------------------------------------------------
!
ALLOCATE(W%XDIR_ALB (ILU))
ALLOCATE(W%XSCA_ALB (ILU))
ALLOCATE(W%XEMIS    (ILU))
W%XDIR_ALB = 0.0
W%XSCA_ALB = 0.0
W%XEMIS    = 0.0
!
 CALL UPDATE_RAD_WATER(W%CWAT_ALB,W%XTS,PZENITH,XTT,W%XEMIS,W%XDIR_ALB,&
                       W%XSCA_ALB,PDIR_ALB,PSCA_ALB,PEMIS,PTSRAD )  
!
PTSURF(:) = W%XTS(:)
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
 CALL INIT_CHEMICAL_n(ILUOUT, KSV, HSV, CHW%NBEQ, CHW%CSV, CHW%NAEREQ,            &
                     CHW%NSV_CHSBEG, CHW%NSV_CHSEND, CHW%NSV_AERBEG, CHW%NSV_AEREND, &
                     CHW%CCH_NAMES, CHW%CAER_NAMES, CHW%NDSTEQ, CHW%NSV_DSTBEG,      &
                     CHW%NSV_DSTEND, CHW%NSLTEQ, CHW%NSV_SLTBEG, CHW%NSV_SLTEND,     &
                     HDSTNAMES=CHW%CDSTNAMES, HSLTNAMES=CHW%CSLTNAMES        )
!
!* depositiion scheme
!

IF (CHW%NBEQ>0 .AND. CHW%CCH_DRY_DEP=='WES89') THEN
  ALLOCATE(CHW%XDEP(ILU,CHW%NBEQ))
ELSE
  ALLOCATE(CHW%XDEP(0,0))
END IF
!
!-------------------------------------------------------------------------------
!
!*       7.     diagnostics initialization
!               --------------------------
!
 CALL DIAG_WATFLUX_INIT_n(DGU, DGW, W, &
                          HPROGRAM,ILU,KSW)
!
!-------------------------------------------------------------------------------
!
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('INIT_WATFLUX_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE INIT_WATFLUX_n
