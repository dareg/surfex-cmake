!     #########
      SUBROUTINE DIAG_FLAKE_INIT_n(HPROGRAM,KLU,KSW)
!     #####################
!
!!****  *DIAG_FLAKE_INIT_n* - routine to initialize FLAKE diagnostic variables
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
!!      Original    01/2004 
!!       V.Masson   10/2013 Adds min and max 2m parameters
!!      B. Decharme  04/2013 new diag
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_SFX_OASIS,      ONLY : LCPL_LAKE
USE MODD_DIAG_SURF_ATM_n,ONLY : LREAD_BUDGETC
!
USE MODD_DIAG_FLAKE_n,   ONLY : N2M, LSURF_BUDGET, LCOEF, LSURF_VARS,   &
                                  LSURF_BUDGETC, LRESET_BUDGETC,        &
                                  XRN, XH, XLE, XLEI, XGFLUX, XEVAP,    &
                                  XSUBL, XRI, XCD, XCH, XCE, XZ0, XZ0H, &
                                  XT2M, XQ2M, XHU2M, XT2M_MIN, XT2M_MAX,&
                                  XZON10M, XMER10M, XQS,                &
                                  XSWD, XSWU, XLWD, XLWU,               &
                                  XSWBD, XSWBU, XFMU, XFMV,             &
                                  XRNC, XHC, XLEC, XLEIC, XGFLUXC,      &
                                  XEVAPC, XSUBLC,                       &
                                  XSWDC, XSWUC, XLWDC, XLWUC, XFMUC,    &
                                  XFMVC, XDIAG_TS, XHU2M_MIN, XHU2M_MAX,&
                                  XWIND10M, XWIND10M_MAX, XALBT, XSWE  
!
USE MODD_DIAG_MISC_FLAKE_n,    ONLY : LWATER_PROFILE , XZWAT_PROFILE,     &
                                      XZW_PROFILE, XTW_PROFILE
!
USE MODD_FLAKE_n,          ONLY : XCPL_FLAKE_EVAP, &
                                  XCPL_FLAKE_RAIN, &
                                  XCPL_FLAKE_SNOW
!
USE MODI_READ_SURF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
INTEGER, INTENT(IN) :: KLU   ! size of arrays
INTEGER, INTENT(IN) :: KSW   ! number of SW spectral bands
CHARACTER(LEN=6), INTENT(IN):: HPROGRAM  ! program calling
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IVERSION
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
CHARACTER(LEN=12) :: YREC           ! Name of the article to be read
REAL(KIND=JPRB)   :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!* surface energy budget
!
IF (LHOOK) CALL DR_HOOK('DIAG_FLAKE_INIT_N',0,ZHOOK_HANDLE)
!
ALLOCATE(XDIAG_TS(KLU))
XDIAG_TS = XUNDEF
!
IF (LSURF_BUDGET.OR.LSURF_BUDGETC) THEN
  ALLOCATE(XRN     (KLU))
  ALLOCATE(XH      (KLU))
  ALLOCATE(XLE     (KLU))
  ALLOCATE(XLEI    (KLU))
  ALLOCATE(XGFLUX  (KLU))
  ALLOCATE(XEVAP   (KLU))
  ALLOCATE(XSUBL   (KLU))
  ALLOCATE(XSWD    (KLU))
  ALLOCATE(XSWU    (KLU))
  ALLOCATE(XLWD    (KLU))
  ALLOCATE(XLWU    (KLU))
  ALLOCATE(XSWBD   (KLU,KSW))
  ALLOCATE(XSWBU   (KLU,KSW))
  ALLOCATE(XFMU    (KLU))
  ALLOCATE(XFMV    (KLU))
  ALLOCATE(XALBT   (KLU))
  ALLOCATE(XSWE    (KLU))
  !
  XRN      = XUNDEF
  XH       = XUNDEF
  XLE      = XUNDEF
  XLEI     = XUNDEF
  XGFLUX   = XUNDEF
  XEVAP    = XUNDEF
  XSUBL    = XUNDEF  
  XSWD     = XUNDEF
  XSWU     = XUNDEF
  XLWD     = XUNDEF
  XLWU     = XUNDEF
  XSWBD    = XUNDEF
  XSWBU    = XUNDEF
  XFMU     = XUNDEF
  XFMV     = XUNDEF
  XALBT    = XUNDEF
  XSWE     = XUNDEF
ELSE
  ALLOCATE(XRN     (0))
  ALLOCATE(XH      (0))
  ALLOCATE(XLE     (0))
  ALLOCATE(XLEI    (0))
  ALLOCATE(XGFLUX  (0))
  ALLOCATE(XEVAP   (0))
  ALLOCATE(XSUBL   (0))  
  ALLOCATE(XSWD    (0))
  ALLOCATE(XSWU    (0))
  ALLOCATE(XLWD    (0))
  ALLOCATE(XLWU    (0))
  ALLOCATE(XSWBD   (0,0))
  ALLOCATE(XSWBU   (0,0))
  ALLOCATE(XFMU    (0))
  ALLOCATE(XFMV    (0))
  ALLOCATE(XALBT   (0))
  ALLOCATE(XSWE    (0))
END IF
!
!* cumulative surface energy budget
!
IF (LSURF_BUDGETC) THEN
!    
  ALLOCATE(XRNC    (KLU))
  ALLOCATE(XHC     (KLU))
  ALLOCATE(XLEC    (KLU))
  ALLOCATE(XLEIC   (KLU))
  ALLOCATE(XGFLUXC (KLU))
  ALLOCATE(XEVAPC  (KLU))
  ALLOCATE(XSUBLC  (KLU))  
  ALLOCATE(XSWDC   (KLU))
  ALLOCATE(XSWUC   (KLU))
  ALLOCATE(XLWDC   (KLU))
  ALLOCATE(XLWUC   (KLU))
  ALLOCATE(XFMUC   (KLU))
  ALLOCATE(XFMVC   (KLU))
!
  IF (.NOT. LREAD_BUDGETC) THEN        
     XRNC    = 0.0
     XHC     = 0.0
     XLEC    = 0.0
     XLEIC   = 0.0
     XGFLUXC = 0.0
     XEVAPC  = 0.0
     XSUBLC  = 0.0
     XSWDC   = 0.0
     XSWUC   = 0.0
     XLWDC   = 0.0
     XLWUC   = 0.0
     XFMUC   = 0.0
     XFMVC   = 0.0
  ELSEIF (LREAD_BUDGETC.AND.LRESET_BUDGETC) THEN
     XRNC    = 0.0
     XHC     = 0.0
     XLEC    = 0.0
     XLEIC   = 0.0
     XGFLUXC = 0.0
     XEVAPC  = 0.0
     XSUBLC  = 0.0     
     XSWDC   = 0.0
     XSWUC   = 0.0
     XLWDC   = 0.0
     XLWUC   = 0.0
     XFMUC   = 0.0
     XFMVC   = 0.0
  ELSE
     CALL READ_SURF(HPROGRAM,'VERSION',IVERSION,IRESP)
     IF (IVERSION<8)THEN
       XRNC    = 0.0
       XHC     = 0.0
       XLEC    = 0.0
       XLEIC   = 0.0
       XGFLUXC = 0.0
       XEVAPC  = 0.0
       XSUBLC  = 0.0     
       XSWDC   = 0.0
       XSWUC   = 0.0
       XLWDC   = 0.0
       XLWUC   = 0.0
       XFMUC   = 0.0
       XFMVC   = 0.0             
     ELSE
       YREC='RNC_WAT'
       CALL READ_SURF(HPROGRAM,YREC,XRNC,IRESP)
       YREC='HC_WAT'
       CALL READ_SURF(HPROGRAM,YREC,XHC ,IRESP)
       YREC='LEC_WAT'
       CALL READ_SURF(HPROGRAM,YREC,XLEC,IRESP)
       YREC='LEIC_WAT'
       CALL READ_SURF(HPROGRAM,YREC,XLEIC,IRESP)     
       YREC='GFLUXC_WAT'
       CALL READ_SURF(HPROGRAM,YREC,XGFLUXC,IRESP)
       YREC='SWDC_WAT'
       CALL READ_SURF(HPROGRAM,YREC,XSWDC,IRESP)
       YREC='SWUC_WAT'
       CALL READ_SURF(HPROGRAM,YREC,XSWUC,IRESP)
       YREC='LWDC_WAT'
       CALL READ_SURF(HPROGRAM,YREC,XLWDC,IRESP)
       YREC='LWUC_WAT'
       CALL READ_SURF(HPROGRAM,YREC,XLWUC,IRESP)
       YREC='FMUC_WAT'
       CALL READ_SURF(HPROGRAM,YREC,XFMUC,IRESP)
       YREC='FMVC_WAT'
       CALL READ_SURF(HPROGRAM,YREC,XFMVC,IRESP)
       YREC='EVAPC_WAT'
        CALL READ_SURF(HPROGRAM,YREC,XEVAPC,IRESP)
        YREC='SUBLC_WAT'
        CALL READ_SURF(HPROGRAM,YREC,XSUBLC,IRESP)              
      ENDIF
!
  ENDIF   
ELSE
  ALLOCATE(XRNC    (0))
  ALLOCATE(XHC     (0))
  ALLOCATE(XLEC    (0))
  ALLOCATE(XLEIC   (0))
  ALLOCATE(XGFLUXC (0))
  ALLOCATE(XEVAPC  (0))
  ALLOCATE(XSUBLC  (0))  
  ALLOCATE(XSWDC   (0))
  ALLOCATE(XSWUC   (0))
  ALLOCATE(XLWDC   (0))
  ALLOCATE(XLWUC   (0))
  ALLOCATE(XFMUC   (0))
  ALLOCATE(XFMVC   (0))  
ENDIF
!
!* parameters at 2m
!
IF (N2M>=1) THEN
  ALLOCATE(XRI     (KLU))
  ALLOCATE(XT2M    (KLU))
  ALLOCATE(XT2M_MIN(KLU))
  ALLOCATE(XT2M_MAX(KLU))
  ALLOCATE(XQ2M    (KLU))
  ALLOCATE(XHU2M   (KLU))
  ALLOCATE(XHU2M_MIN(KLU))
  ALLOCATE(XHU2M_MAX(KLU))
  ALLOCATE(XZON10M (KLU))
  ALLOCATE(XMER10M (KLU))
  ALLOCATE(XWIND10M (KLU))
  ALLOCATE(XWIND10M_MAX(KLU))
  !
  XRI      = XUNDEF
  XT2M     = XUNDEF
  XT2M_MIN = XUNDEF
  XT2M_MAX = 0.0
  XQ2M     = XUNDEF
  XHU2M    = XUNDEF
  XHU2M_MIN= XUNDEF
  XHU2M_MAX=-XUNDEF
  XZON10M  = XUNDEF
  XMER10M  = XUNDEF
  XWIND10M = XUNDEF
  XWIND10M_MAX = 0.0
ELSE
  ALLOCATE(XRI      (0))
  ALLOCATE(XT2M     (0))
  ALLOCATE(XT2M_MIN (0))
  ALLOCATE(XT2M_MAX (0))
  ALLOCATE(XQ2M     (0))
  ALLOCATE(XHU2M    (0))
  ALLOCATE(XHU2M_MIN(0))
  ALLOCATE(XHU2M_MAX(0))
  ALLOCATE(XZON10M  (0))
  ALLOCATE(XMER10M  (0))
  ALLOCATE(XWIND10M (0))
  ALLOCATE(XWIND10M_MAX(0))
END IF
!
!* transfer coefficients
!
IF (LCOEF) THEN
  ALLOCATE(XCD     (KLU))
  ALLOCATE(XCH     (KLU))
  ALLOCATE(XCE     (KLU))
  ALLOCATE(XZ0     (KLU))
  ALLOCATE(XZ0H    (KLU))
  !
  XCD      = XUNDEF
  XCH      = XUNDEF
  XCE      = XUNDEF
  XZ0      = XUNDEF
  XZ0H     = XUNDEF
ELSE
  ALLOCATE(XCD     (0))
  ALLOCATE(XCH     (0))
  ALLOCATE(XCE     (0))
  ALLOCATE(XZ0     (0))
  ALLOCATE(XZ0H    (0))
END IF
!
!* surface humidity
!
IF (LSURF_VARS) THEN
  ALLOCATE(XQS     (KLU))
  !
  XQS      = XUNDEF
ELSE
  ALLOCATE(XQS     (0))
END IF
!
!* Flake temperature profile
!
IF (LWATER_PROFILE) THEN
   ALLOCATE (XZW_PROFILE(COUNT(XZWAT_PROFILE/= XUNDEF))) 
   ALLOCATE (XTW_PROFILE(COUNT(XZWAT_PROFILE/= XUNDEF),KLU)) 
   XZW_PROFILE=XZWAT_PROFILE(:COUNT(XZWAT_PROFILE /= XUNDEF))
 ELSE
   ALLOCATE (XZW_PROFILE(0)) 
   ALLOCATE (XTW_PROFILE(0,0)) 
 END IF
!
!* Coupling field with earth systme model
!
!
IF(LCPL_LAKE)THEN
!
  ALLOCATE(XCPL_FLAKE_EVAP(KLU))
  ALLOCATE(XCPL_FLAKE_RAIN(KLU))
  ALLOCATE(XCPL_FLAKE_SNOW(KLU))
  XCPL_FLAKE_EVAP(:) = 0.0
  XCPL_FLAKE_RAIN(:) = 0.0
  XCPL_FLAKE_SNOW(:) = 0.0
!
ELSE
!
  ALLOCATE(XCPL_FLAKE_EVAP(0))
  ALLOCATE(XCPL_FLAKE_RAIN(0))
  ALLOCATE(XCPL_FLAKE_SNOW(0))
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('DIAG_FLAKE_INIT_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_FLAKE_INIT_n
