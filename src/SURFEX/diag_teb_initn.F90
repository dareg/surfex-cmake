!     #########
      SUBROUTINE DIAG_TEB_INIT_n(HPROGRAM,KLU,KSW)
!     #####################
!
!!****  *DIAG_TEB_INIT_n* - routine to initialize TEB diagnostic variables
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
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
USE MODD_TYPE_DATE_SURF
USE MODD_DIAG_TEB_n, ONLY : N2M, LSURF_BUDGET, LCOEF, LSURF_VARS, &
                              XRN, XH, XLE, XGFLUX, XRI,            &
                              XCD, XCH, XCE, XZ0, XZ0H,             &
                              XT2M, XQ2M, XHU2M,                    &
                              XZON10M, XMER10M, XQS,                &
                              XSWD, XSWU, XSWBD, XSWBU, XLWD, XLWU, &
                              XFMU, XFMV  
!
USE MODD_DIAG_SURF_ATM_n,   ONLY : LREAD_BUDGETC
USE MODD_DIAG_MISC_TEB_n,   ONLY : LSURF_MISC_BUDGET,LSURF_BUDGETC,         &
                                     LRESET_BUDGETC,                          &
                                     XQF_BLD, XQF_TOWN, XDQS_TOWN,            &
                                     XFLX_BLD,XTI_BLD_EQ,XQF_BLDWFR,          &
                                     XTI_BLDWFR,                              &
                                     XRN_ROAD, XH_ROAD, XLE_ROAD, XGFLUX_ROAD,&
                                     XRN_WALL, XH_WALL, XGFLUX_WALL,          &
                                     XRN_ROOF, XH_ROOF, XLE_ROOF, XGFLUX_ROOF,&
                                     XRUNOFF, XRUNOFFC,                       &
                                     XRN_GARDEN,XH_GARDEN,XLE_GARDEN,         &
                                     XGFLUX_GARDEN,                           &
                                     XRN_BLT,XH_BLT,XLE_BLT,XGFLUX_BLT,       &
                                     XABS_SW_ROOF ,XABS_SW_SNOW_ROOF,         &
                                     XABS_LW_ROOF ,XABS_LW_SNOW_ROOF,         &
                                     XABS_SW_ROAD ,XABS_SW_SNOW_ROAD,         &
                                     XABS_LW_ROAD ,XABS_LW_SNOW_ROAD,         &
                                     XABS_SW_WALL ,                           &
                                     XABS_LW_WALL ,                           &
                                     XABS_SW_GARDEN,XABS_LW_GARDEN  
!
USE MODI_READ_SURF
!
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
INTEGER, INTENT(IN) :: KSW   ! spectral bands
CHARACTER(LEN=6), INTENT(IN):: HPROGRAM  ! program calling
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
CHARACTER(LEN=16) :: YREC           ! Name of the article to be read
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!* surface energy budget
!
IF (LHOOK) CALL DR_HOOK('DIAG_TEB_INIT_N',0,ZHOOK_HANDLE)
IF (LSURF_BUDGET) THEN
  ALLOCATE(XRN     (KLU))
  ALLOCATE(XH      (KLU))
  ALLOCATE(XLE     (KLU))
  ALLOCATE(XGFLUX  (KLU))
  ALLOCATE(XSWD    (KLU))
  ALLOCATE(XSWU    (KLU))
  ALLOCATE(XSWBD   (KLU,KSW))
  ALLOCATE(XSWBU   (KLU,KSW))
  ALLOCATE(XLWD    (KLU))
  ALLOCATE(XLWU    (KLU))
  ALLOCATE(XFMU    (KLU))
  ALLOCATE(XFMV    (KLU))
  !
  XRN      = XUNDEF
  XH       = XUNDEF
  XLE      = XUNDEF
  XGFLUX   = XUNDEF
  XSWD     = XUNDEF
  XSWU     = XUNDEF
  XSWBD    = XUNDEF
  XSWBU    = XUNDEF
  XLWD     = XUNDEF
  XLWU     = XUNDEF
  XFMU     = XUNDEF
  XFMV     = XUNDEF
ELSE
  ALLOCATE(XRN     (0))
  ALLOCATE(XH      (0))
  ALLOCATE(XLE     (0))
  ALLOCATE(XGFLUX  (0))
  ALLOCATE(XSWD    (0))
  ALLOCATE(XSWU    (0))
  ALLOCATE(XLWD    (0))
  ALLOCATE(XLWU    (0))
  ALLOCATE(XSWBD   (0,0))
  ALLOCATE(XSWBU   (0,0))
  ALLOCATE(XFMU    (0))
  ALLOCATE(XFMV    (0))
END IF
!
!* parameters at 2m
!
IF (N2M>=1) THEN
  ALLOCATE(XRI     (KLU))
  ALLOCATE(XT2M    (KLU))
  ALLOCATE(XQ2M    (KLU))
  ALLOCATE(XHU2M   (KLU))
  ALLOCATE(XZON10M (KLU))
  ALLOCATE(XMER10M (KLU))
  !
  XRI      = XUNDEF
  XT2M     = XUNDEF
  XQ2M     = XUNDEF
  XHU2M    = XUNDEF
  XZON10M  = XUNDEF
  XMER10M  = XUNDEF
ELSE
  ALLOCATE(XRI      (0))
  ALLOCATE(XT2M     (0))
  ALLOCATE(XQ2M     (0))
  ALLOCATE(XHU2M    (0))
  ALLOCATE(XZON10M  (0))
  ALLOCATE(XMER10M  (0))
END IF
!
!* miscellaneous fields
!
IF (LSURF_MISC_BUDGET) THEN
  ALLOCATE(XQF_BLD    (KLU))
  ALLOCATE(XQF_BLDWFR (KLU))
  ALLOCATE(XFLX_BLD   (KLU))
  ALLOCATE(XTI_BLD_EQ (KLU))
  ALLOCATE(XTI_BLDWFR (KLU))
  ALLOCATE(XQF_TOWN   (KLU))
  ALLOCATE(XDQS_TOWN  (KLU))
  ALLOCATE(XRN_ROAD     (KLU))
  ALLOCATE(XH_ROAD      (KLU))
  ALLOCATE(XLE_ROAD     (KLU))
  ALLOCATE(XGFLUX_ROAD  (KLU))
  ALLOCATE(XRN_WALL     (KLU))
  ALLOCATE(XH_WALL      (KLU))
  ALLOCATE(XGFLUX_WALL  (KLU))
  ALLOCATE(XRN_ROOF     (KLU))
  ALLOCATE(XH_ROOF      (KLU))
  ALLOCATE(XLE_ROOF     (KLU))
  ALLOCATE(XGFLUX_ROOF  (KLU))
  ALLOCATE(XRUNOFF  (KLU))
  ALLOCATE(XRN_GARDEN   (KLU))
  ALLOCATE(XH_GARDEN    (KLU))
  ALLOCATE(XLE_GARDEN   (KLU))
  ALLOCATE(XGFLUX_GARDEN(KLU))
  ALLOCATE(XRN_BLT      (KLU))
  ALLOCATE(XH_BLT       (KLU))
  ALLOCATE(XLE_BLT      (KLU))
  ALLOCATE(XGFLUX_BLT   (KLU))
  !
  ALLOCATE(XABS_SW_ROOF      (KLU))
  ALLOCATE(XABS_SW_SNOW_ROOF (KLU))
  ALLOCATE(XABS_LW_ROOF      (KLU))
  ALLOCATE(XABS_LW_SNOW_ROOF (KLU))
  ALLOCATE(XABS_SW_ROAD      (KLU))
  ALLOCATE(XABS_SW_SNOW_ROAD (KLU))
  ALLOCATE(XABS_LW_ROAD      (KLU))
  ALLOCATE(XABS_LW_SNOW_ROAD (KLU))
  ALLOCATE(XABS_SW_WALL      (KLU))
  ALLOCATE(XABS_LW_WALL      (KLU))
  ALLOCATE(XABS_SW_GARDEN    (KLU))
  ALLOCATE(XABS_LW_GARDEN    (KLU))
  !
  XQF_BLD     = XUNDEF
  XQF_BLDWFR  = XUNDEF
  XFLX_BLD    = XUNDEF
  XTI_BLD_EQ  = XUNDEF
  XTI_BLDWFR  = XUNDEF
  XQF_TOWN    = XUNDEF
  XDQS_TOWN   = XUNDEF
  XRN_ROAD      = XUNDEF
  XH_ROAD       = XUNDEF
  XLE_ROAD      = XUNDEF
  XGFLUX_ROAD   = XUNDEF
  XRN_WALL      = XUNDEF
  XH_WALL       = XUNDEF
  XGFLUX_WALL   = XUNDEF
  XRN_ROOF      = XUNDEF
  XH_ROOF       = XUNDEF
  XLE_ROOF      = XUNDEF
  XGFLUX_ROOF   = XUNDEF 
  XRUNOFF   = XUNDEF 
  XRN_GARDEN    = XUNDEF
  XH_GARDEN     = XUNDEF
  XLE_GARDEN    = XUNDEF
  XGFLUX_GARDEN = XUNDEF  
  XRN_BLT       = XUNDEF
  XH_BLT        = XUNDEF
  XLE_BLT       = XUNDEF
  XGFLUX_BLT    = XUNDEF  
!
  XABS_SW_ROOF       = XUNDEF  
  XABS_SW_SNOW_ROOF  = XUNDEF  
  XABS_LW_ROOF       = XUNDEF  
  XABS_LW_SNOW_ROOF  = XUNDEF  
  XABS_SW_ROAD       = XUNDEF  
  XABS_SW_SNOW_ROAD  = XUNDEF  
  XABS_LW_ROAD       = XUNDEF  
  XABS_LW_SNOW_ROAD  = XUNDEF  
  XABS_SW_WALL       = XUNDEF  
  XABS_LW_WALL       = XUNDEF  
  XABS_SW_GARDEN     = XUNDEF  
  XABS_LW_GARDEN     = XUNDEF  
ELSE
  ALLOCATE(XQF_BLD    (0))
  ALLOCATE(XQF_BLDWFR (0))
  ALLOCATE(XFLX_BLD   (0))
  ALLOCATE(XTI_BLD_EQ (0))
  ALLOCATE(XTI_BLDWFR (0))
  ALLOCATE(XQF_TOWN   (0))
  ALLOCATE(XDQS_TOWN  (0))
  ALLOCATE(XRN_ROAD     (0))
  ALLOCATE(XH_ROAD      (0))
  ALLOCATE(XLE_ROAD     (0))
  ALLOCATE(XGFLUX_ROAD  (0))
  ALLOCATE(XRN_WALL     (0))
  ALLOCATE(XH_WALL      (0))
  ALLOCATE(XGFLUX_WALL  (0))
  ALLOCATE(XRN_ROOF     (0))
  ALLOCATE(XH_ROOF      (0))
  ALLOCATE(XLE_ROOF     (0))
  ALLOCATE(XGFLUX_ROOF  (0))
  ALLOCATE(XRUNOFF  (0))
  ALLOCATE(XRN_GARDEN   (0))
  ALLOCATE(XH_GARDEN    (0))
  ALLOCATE(XLE_GARDEN   (0))
  ALLOCATE(XGFLUX_GARDEN(0))
  ALLOCATE(XRN_BLT      (0))
  ALLOCATE(XH_BLT       (0))
  ALLOCATE(XLE_BLT      (0))
  ALLOCATE(XGFLUX_BLT   (0))
!
  ALLOCATE(XABS_SW_ROOF      (0))
  ALLOCATE(XABS_SW_SNOW_ROOF (0))
  ALLOCATE(XABS_LW_ROOF      (0))
  ALLOCATE(XABS_LW_SNOW_ROOF (0))
  ALLOCATE(XABS_SW_ROAD      (0))
  ALLOCATE(XABS_SW_SNOW_ROAD (0))
  ALLOCATE(XABS_LW_ROAD      (0))
  ALLOCATE(XABS_LW_SNOW_ROAD (0))
  ALLOCATE(XABS_SW_WALL      (0))
  ALLOCATE(XABS_LW_WALL      (0))
  ALLOCATE(XABS_SW_GARDEN    (0))
  ALLOCATE(XABS_LW_GARDEN    (0))
ENDIF
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
!* surface cumulated energy budget
!
IF (LSURF_BUDGETC) THEN
  ALLOCATE(XRUNOFFC     (KLU))
  !
  IF (.NOT. LREAD_BUDGETC) THEN
        XRUNOFFC = 0.
  ELSEIF (LREAD_BUDGETC.AND.LRESET_BUDGETC) THEN
        XRUNOFFC = 0.
  ELSE
        YREC='RUNOFFC_TWN'
        CALL READ_SURF(HPROGRAM,YREC,XRUNOFFC,IRESP)
  ENDIF
  ELSE
  ALLOCATE(XRUNOFFC     (0))
ENDIF    
IF (LHOOK) CALL DR_HOOK('DIAG_TEB_INIT_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_TEB_INIT_n
