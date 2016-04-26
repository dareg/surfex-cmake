!     #########
      SUBROUTINE INIT_FROM_DATA_TEB_VEG_1P_n (DTI, K, P, PEK, KDECADE, HPHOTO, OUPDATE, OFIX, OTIME, OALB)    
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
USE MODD_ISBA_n, ONLY : ISBA_K_t, ISBA_P_t, ISBA_PE_t
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
TYPE(DATA_ISBA_t), INTENT(INOUT) :: DTI
!
INTEGER,                INTENT(IN)    :: KDECADE
 CHARACTER(LEN=*),       INTENT(IN)    :: HPHOTO  ! type of photosynthesis
!
LOGICAL, INTENT(IN) :: OUPDATE
LOGICAL, INTENT(IN) :: OFIX
LOGICAL, INTENT(IN) :: OTIME
LOGICAL, INTENT(IN) :: OALB
!
TYPE(ISBA_K_t), INTENT(INOUT) :: K
TYPE(ISBA_P_t), INTENT(INOUT) :: P
TYPE(ISBA_PE_t), INTENT(INOUT) :: PEK
!
REAL, DIMENSION(:), ALLOCATABLE :: ZWG1
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
IF (LHOOK) CALL DR_HOOK('INIT_FROM_DATA_TEB_VEG_1P_N',0,ZHOOK_HANDLE)
IF (DTI%NTIME==12) THEN
  ITIME = (KDECADE+2)/3    
ELSEIF (DTI%NTIME==1) THEN
  ITIME = 1
ENDIF
!
!
IF (OFIX) THEN
!

  IF (SIZE(P%XH_TREE)>0) P%XH_TREE(:) = DTI%XPAR_H_TREE(:,1)
!
  P%XZ0_O_Z0H(:) = DTI%XPAR_Z0_O_Z0H(:,1)
!
!---------------------------------------------------------------------------------
!
!* soil layers
!  -----------
!
  P%XDG(:,:) = DTI%XPAR_DG(:,:,1)
!
!* cumulative root fraction
!
  IF (SIZE(P%XROOTFRAC)>0) P%XROOTFRAC(:,:) = DTI%XPAR_ROOTFRAC(:,:,1)
!
!* soil ice for runoff
!
 P%XD_ICE(:) = DTI%XPAR_DICE(:,1)
!
  IF (SIZE(P%XDMAX)>0) P%XDMAX(:) = DTI%XPAR_DMAX(:,1)

  IF (SIZE(P%XRE25)>0) P%XRE25(:) = DTI%XPAR_RE25(:,1)


ENDIF
!
!
!*    2.      SECONDARY VARIABLES
!             -------------------
!
!*    2.1     fields on natural surfaces only, taking into account patches/ 
!             -------------------------------
!
IF (OTIME) THEN
!
! vegetation fraction
! -------------------
!
 PEK%XVEG(:) =  DTI%XPAR_VEG (:,ITIME,1)
!
! Leaf Aera Index
! ---------------
!
 PEK%XLAI(:) = DTI%XPAR_LAI (:,ITIME,1)
!
! roughness length
! ----------------
!
  PEK%XZ0(:) =  DTI%XPAR_Z0 (:,ITIME,1)
!
!emis-eco
!--------
!
  PEK%XEMIS(:) =  DTI%XPAR_EMIS (:,ITIME,1)
!
 IF (.NOT.OUPDATE) THEN
!---------------------------------------------------------------------------------
! 
!* 1/Rsmin
!
  PEK%XRSMIN(:) = DTI%XPAR_RSMIN(:,1)
!
!* other vegetation parameters
!
  PEK%XGAMMA(:) = DTI%XPAR_GAMMA(:,1)
  PEK%XWRMAX_CF(:) = DTI%XPAR_WRMAX_CF(:,1)
!
!
  PEK%XRGL(:) = DTI%XPAR_RGL(:,1)
  PEK%XCV(:) = DTI%XPAR_CV(:,1)
!
!---------------------------------------------------------------------------------
  PEK%XALBNIR_VEG(:) = DTI%XPAR_ALBNIR_VEG(:,1)
  PEK%XALBVIS_VEG(:) = DTI%XPAR_ALBVIS_VEG(:,1)
  PEK%XALBUV_VEG(:) = DTI%XPAR_ALBUV_VEG(:,1)
!
  IF (SIZE(PEK%XGMES)>0) PEK%XGMES(:) = DTI%XPAR_GMES(:,1)

  IF (SIZE(PEK%XBSLAI)>0) PEK%XBSLAI(:) = DTI%XPAR_BSLAI(:,1)

  IF (SIZE(PEK%XSEFOLD)>0) PEK%XSEFOLD(:) = DTI%XPAR_SEFOLD(:,1)

  IF (SIZE(PEK%XGC)>0) PEK%XGC(:) = DTI%XPAR_GC(:,1)
!
  IF (SIZE(PEK%XLAIMIN)>0) PEK%XLAIMIN(:) = DTI%XPAR_LAIMIN(:,1)
!
  IF (SIZE(PEK%XCE_NITRO)>0) PEK%XCE_NITRO(:) = DTI%XPAR_CE_NITRO(:,1)
!
  IF (SIZE(PEK%XCF_NITRO)>0) PEK%XCF_NITRO(:) = DTI%XPAR_CF_NITRO(:,1)
!
  IF (SIZE(PEK%XCNA_NITRO)>0) PEK%XCNA_NITRO(:) = DTI%XPAR_CNA_NITRO(:,1)
!
  IF (SIZE(PEK%XF2I)>0) PEK%XF2I(:) = DTI%XPAR_F2I(:,1)
!
  IF (SIZE(PEK%LSTRESS)>0) PEK%LSTRESS(:) = DTI%LPAR_STRESS(:,1)
!
 ENDIF
ENDIF
!
IF (OALB) THEN
  !
  ALLOCATE(ZWGSAT(SIZE(K%XALBVIS_DRY)))
  ALLOCATE(ZWG1(SIZE(K%XALBVIS_DRY)))
  ZWGSAT(:) = 0.
  ZWG1  (:) = 0.
  CALL SOIL_ALBEDO('DRY',ZWGSAT, ZWG1, K, PEK, "ALL" ) 
  DEALLOCATE(ZWGSAT,ZWG1)
  !
ENDIF
!
IF (LHOOK) CALL DR_HOOK('INIT_FROM_DATA_TEB_VEG_1P_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INIT_FROM_DATA_TEB_VEG_1P_n
