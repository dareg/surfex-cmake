!     #########
      SUBROUTINE INIT_FROM_DATA_GRDN_n (DTGD, KDECADE, HPHOTO, OUPDATE, VMX, VMT, VMA, VMIP)    
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
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original   01/2004
!     
!!      P Le Moigne 09/2005 AGS modifs of L. Jarlan
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_DATA_ISBA_n, ONLY : DATA_ISBA_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_FIX_t, ISBA_PARAM_TIME_t, ISBA_PARAM_ALB_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t
!
USE MODI_SOIL_ALBEDO
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
!
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTGD
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
!-------------------------------------------------------------------------------
!
!*    1.      TIME INITIALIZATION
!             -------------------
!
! data every month
IF (LHOOK) CALL DR_HOOK('INIT_FROM_DATA_GRDN_N',0,ZHOOK_HANDLE)
IF (DTGD%NTIME==12) THEN
  ITIME = (KDECADE+2)/3    
ELSEIF (DTGD%NTIME==1) THEN
  ITIME = 1
ENDIF
!
!*    2.      SECONDARY VARIABLES
!             -------------------
!
!*    2.1     fields on natural surfaces only, taking into account patches/ 
!             -------------------------------
!
!
IF (PRESENT(VMX)) THEN
!
  IF (SIZE(VMX%XH_TREE)>0) VMX%XH_TREE = DTGD%XPAR_H_TREE
!
  VMX%XVEGTYPE = DTGD%XPAR_VEGTYPE
!
 VMX%XZ0_O_Z0H = DTGD%XPAR_Z0_O_Z0H
!
!---------------------------------------------------------------------------------
!
!* soil layers
!  -----------
!
 VMX%XDG = DTGD%XPAR_DG
!
!* cumulative root fraction
!
  IF (SIZE(VMX%XROOTFRAC)>0) VMX%XROOTFRAC = DTGD%XPAR_ROOTFRAC
!
!* soil ice for runoff
!
 VMX%XD_ICE = DTGD%XPAR_DICE
!
  IF (SIZE(VMX%XDMAX)>0) VMX%XDMAX = DTGD%XPAR_DMAX

  IF (SIZE(VMX%XRE25)>0) VMX%XRE25 = DTGD%XPAR_RE25

ENDIF
!
IF (PRESENT(VMT)) THEN
!
! vegetation fraction
! -------------------
!
 VMT%XVEG =  DTGD%XPAR_VEG (:,ITIME,:)
!
! Leaf Aera Index
! ---------------
!
 VMT%XLAI = DTGD%XPAR_LAI (:,ITIME,:)
!
! roughness length
! ----------------
!
  VMT%XZ0 =  DTGD%XPAR_Z0 (:,ITIME,:)
!
!emis-eco
!--------
!
  VMT%XEMIS =  DTGD%XPAR_EMIS (:,ITIME,:)
!
 IF (.NOT.OUPDATE) THEN
!---------------------------------------------------------------------------------
! 
!* 1/Rsmin
!
  VMT%XRSMIN = DTGD%XPAR_RSMIN
!
!* other vegetation parameters
!
  VMT%XGAMMA = DTGD%XPAR_GAMMA
  VMT%XWRMAX_CF = DTGD%XPAR_WRMAX_CF
!
!
  VMT%XRGL = DTGD%XPAR_RGL
  VMT%XCV = DTGD%XPAR_CV
!
!---------------------------------------------------------------------------------
  VMT%XALBNIR_VEG = DTGD%XPAR_ALBNIR_VEG
  VMT%XALBVIS_VEG = DTGD%XPAR_ALBVIS_VEG
  VMT%XALBUV_VEG = DTGD%XPAR_ALBUV_VEG
!
  IF (SIZE(VMT%XGMES)>0) VMT%XGMES = DTGD%XPAR_GMES

  IF (SIZE(VMT%XBSLAI)>0) VMT%XBSLAI = DTGD%XPAR_BSLAI

  IF (SIZE(VMT%XSEFOLD)>0) VMT%XSEFOLD = DTGD%XPAR_SEFOLD

  IF (SIZE(VMT%XGC)>0) VMT%XGC = DTGD%XPAR_GC
!
  IF (SIZE(VMT%XLAIMIN)>0) VMT%XLAIMIN = DTGD%XPAR_LAIMIN
!
  IF (SIZE(VMT%XCE_NITRO)>0) VMT%XCE_NITRO = DTGD%XPAR_CE_NITRO
!
  IF (SIZE(VMT%XCF_NITRO)>0) VMT%XCF_NITRO = DTGD%XPAR_CF_NITRO
!
  IF (SIZE(VMT%XCNA_NITRO)>0) VMT%XCNA_NITRO = DTGD%XPAR_CNA_NITRO
!
  IF (SIZE(VMT%XF2I)>0) VMT%XF2I = DTGD%XPAR_F2I
!
  IF (SIZE(VMT%LSTRESS)>0) VMT%LSTRESS = DTGD%LPAR_STRESS
!
 ENDIF
 !
ENDIF
!
IF (PRESENT(VMA).AND.PRESENT(VMIP)) THEN
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
IF (LHOOK) CALL DR_HOOK('INIT_FROM_DATA_GRDN_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INIT_FROM_DATA_GRDN_n
