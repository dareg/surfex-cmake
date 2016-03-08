!     #########
      SUBROUTINE INIT_FROM_DATA_GREENROOF_n (DTGR,  KDECADE, HPHOTO, OUPDATE, VMX, VMT, VMA, VMIP)  
!     ##############################################################
!
!!**** *CONVERT_COVER* convert surface cover classes into secondary 
!!                     physiographic variables for ISBA
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!    Based on init_from_data_grdnn
!!    
!!    AUTHOR
!!    ------
!!
!!    C. de Munck & A. Lemonsu        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original   08/2011
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!

!
!
!
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_FIX_t, ISBA_PARAM_TIME_t, ISBA_PARAM_ALB_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_SOIL_ALBEDO
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
!
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTGR
!
INTEGER,                INTENT(IN)    :: KDECADE
 CHARACTER(LEN=*),       INTENT(IN)    :: HPHOTO  ! type of photosynthesis
!
LOGICAL, INTENT(IN) :: OUPDATE
!
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT), OPTIONAL :: VMIP
TYPE(ISBA_PARAM_FIX_t), INTENT(INOUT), OPTIONAL :: VMX
TYPE(ISBA_PARAM_TIME_t), INTENT(INOUT), OPTIONAL :: VMT
TYPE(ISBA_PARAM_ALB_t), INTENT(INOUT), OPTIONAL :: VMA
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZWG1
REAL, DIMENSION(:), ALLOCATABLE   :: ZWGSAT
!
INTEGER :: ITIME
INTEGER :: ILUOUT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!*    1.      TIME INITIALIZATION
!             -------------------
!
! data every month
IF (LHOOK) CALL DR_HOOK('INIT_FROM_DATA_GREENROOF_N',0,ZHOOK_HANDLE)
IF (DTGR%NTIME==12) THEN
  ITIME = (KDECADE+2)/3      
ELSE
  ITIME = 1
END IF
!
!*    2.      SECONDARY VARIABLES
!             -------------------
!
!
!*    2.1     fields on natural surfaces only, taking into account patches/ 
!             -------------------------------
!
IF (PRESENT(VMX)) THEN
!
  VMX%XVEGTYPE = DTGR%XPAR_VEGTYPE
!
  VMX%XZ0_O_Z0H = DTGR%XPAR_Z0_O_Z0H
!
!---------------------------------------------------------------------------------
!
!* soil layers
!  -----------
!
  VMX%XDG = DTGR%XPAR_DG
!
!* cumulative root fraction
!
  VMX%XROOTFRAC = DTGR%XPAR_ROOTFRAC
!
!* soil ice for runoff
!
  VMX%XD_ICE = DTGR%XPAR_DICE
!
  IF (SIZE(VMX%XDMAX)>0) VMX%XDMAX = DTGR%XPAR_DMAX

  IF (SIZE(VMX%XRE25)>0) VMX%XRE25 = DTGR%XPAR_RE25
!
ENDIF
!
IF (PRESENT(VMT)) THEN
!
! vegetation fraction
! -------------------
!
  VMT%XVEG =  DTGR%XPAR_VEG (:,ITIME,:)
!
! Leaf Area Index
! ---------------
!
  VMT%XLAI = DTGR%XPAR_LAI (:,ITIME,:)
!
! roughness length
! ----------------
!
  VMT%XZ0 =  DTGR%XPAR_Z0 (:,ITIME,:)
!
!emis-eco
!--------
!
  VMT%XEMIS =  DTGR%XPAR_EMIS (:,ITIME,:)
! 
 IF (.NOT.OUPDATE) THEN
!---------------------------------------------------------------------------------
! 
!* 1/Rsmin
!
  IF (SIZE(VMT%XRSMIN)>0) VMT%XRSMIN = DTGR%XPAR_RSMIN
!
!* other vegetation parameters
!
  VMT%XGAMMA = DTGR%XPAR_GAMMA
  VMT%XWRMAX_CF = DTGR%XPAR_WRMAX_CF
!
!
  VMT%XRGL = DTGR%XPAR_RGL
  VMT%XCV = DTGR%XPAR_CV
!
!
!---------------------------------------------------------------------------------
  VMT%XALBNIR_VEG = DTGR%XPAR_ALBNIR_VEG
  VMT%XALBVIS_VEG = DTGR%XPAR_ALBVIS_VEG
  VMT%XALBUV_VEG  = DTGR%XPAR_ALBUV_VEG
!
  IF (SIZE(VMT%XGMES)>0) VMT%XGMES = DTGR%XPAR_GMES

  IF (SIZE(VMT%XBSLAI)>0) VMT%XBSLAI = DTGR%XPAR_BSLAI

  IF (SIZE(VMT%XSEFOLD)>0) VMT%XSEFOLD = DTGR%XPAR_SEFOLD

  IF (SIZE(VMT%XGC)>0) VMT%XGC = DTGR%XPAR_GC
!
  IF (SIZE(VMT%XLAIMIN)>0) VMT%XLAIMIN = DTGR%XPAR_LAIMIN

  IF (SIZE(VMT%XCE_NITRO)>0) VMT%XCE_NITRO = DTGR%XPAR_CE_NITRO

  IF (SIZE(VMT%XCF_NITRO)>0) VMT%XCF_NITRO = DTGR%XPAR_CF_NITRO

  IF (SIZE(VMT%XCNA_NITRO)>0) VMT%XCNA_NITRO = DTGR%XPAR_CNA_NITRO

  IF (SIZE(VMT%XF2I)>0) VMT%XF2I = DTGR%XPAR_F2I
!
  IF (SIZE(VMT%LSTRESS)>0) VMT%LSTRESS = DTGR%LPAR_STRESS
!
 ENDIF
!
ENDIF
!
IF (PRESENT(VMA)) THEN
  !
  ALLOCATE(ZWGSAT(SIZE(VMIP%XALBVIS_DRY)))
  ALLOCATE(ZWG1(SIZE(VMIP%XALBVIS_DRY),1))
  ZWGSAT(:) = 0.
  ZWG1(:,:) = 0.
  CALL SOIL_ALBEDO('DRY',ZWGSAT, ZWG1, VMIP, VMA, "ALL" ) 
  DEALLOCATE(ZWGSAT,ZWG1)
  !
ENDIF
!
IF (LHOOK) CALL DR_HOOK('INIT_FROM_DATA_GREENROOF_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INIT_FROM_DATA_GREENROOF_n
