!#############################################################
SUBROUTINE INIT_VEG_n(IO, OCANOPY, KI, MT, MA, PH_TREE, IR, PVEGTYPE_PATCH, &
                      OSURF_DIAG_ALBEDO, PDIR_ALB, PSCA_ALB, PEMIS_OUT, PTSRAD )  
!#############################################################
!
!!****  *INIT_VEG_n* - routine to initialize ISBA
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
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_TIME_t, ISBA_PARAM_ALB_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
!
USE MODD_TYPE_SNOW
USE MODD_SNOW_PAR,       ONLY : XEMISSN
USE MODD_SURF_PAR,       ONLY : XUNDEF, NUNDEF
!
USE MODI_SET_ROUGH
USE MODI_INIT_SNOW_LW
USE MODI_Z0V_FROM_LAI
USE MODI_VEG_FROM_LAI
USE MODI_EMIS_FROM_VEG
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
TYPE(ISBA_PARAM_TIME_t), INTENT(INOUT) :: MT
TYPE(ISBA_PARAM_ALB_t), INTENT(INOUT) :: MA
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
!
LOGICAL, INTENT(IN) :: OCANOPY
INTEGER, INTENT(IN) :: KI
REAL, DIMENSION(:,:), INTENT(IN) :: PH_TREE
REAL, DIMENSION(:,:,:), INTENT(IN) :: PVEGTYPE_PATCH
!
LOGICAL, INTENT(OUT) :: OSURF_DIAG_ALBEDO
!
REAL, DIMENSION(:,:), INTENT(OUT) :: PDIR_ALB
REAL, DIMENSION(:,:), INTENT(OUT) :: PSCA_ALB
REAL, DIMENSION(:), INTENT(OUT) :: PEMIS_OUT
REAL, DIMENSION(:), INTENT(OUT) :: PTSRAD
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JP  ! loop counter on tiles
INTEGER :: JI     ! loop increment
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!               Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('INIT_VEG_n',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!*      1.     Roughness length option
!              -----------------------
!
 CALL SET_ROUGH(OCANOPY,IO%CROUGH)
!
!-------------------------------------------------------------------------------
!
!*      2.     Radiative fields and snow/flood fracion initialization:
!              -------------------------------------------------------
!
!* snow long-wave properties (not initialized in read_gr_snow)
!
 CALL INIT_SNOW_LW(XEMISSN,IR%TSNOW)
!
!-------------------------------------------------------------------------------
!
!* z0 and vegetation fraction estimated from LAI
IF (IO%CPHOTO=='LAI' .OR. IO%CPHOTO=='LST' .OR. IO%CPHOTO=='NIT' .OR. IO%CPHOTO=='NCB') THEN
  DO JP=1,IO%NPATCH
    DO JI=1,KI    
      IF(MT%XLAI(JI,JP)/=XUNDEF) THEN
        IF (MT%XLAI(JI,JP).LT.MT%XLAIMIN(JI,JP)) THEN
          MT%XLAI(JI,JP) = MT%XLAIMIN(JI,JP)
        ENDIF
           MT%XZ0  (JI,JP) = Z0V_FROM_LAI(MT%XLAI(JI,JP),PH_TREE(JI,JP),PVEGTYPE_PATCH(JI,:,JP),IO%LAGRI_TO_GRASS)
           MT%XVEG (JI,JP) = VEG_FROM_LAI(MT%XLAI(JI,JP),PVEGTYPE_PATCH(JI,:,JP),IO%LAGRI_TO_GRASS)
           MT%XEMIS(JI,JP) = EMIS_FROM_VEG(MT%XVEG(JI,JP),PVEGTYPE_PATCH(JI,:,JP))
        END IF  
     END DO
  END DO
END IF
!
!-------------------------------------------------------------------------------
!
IF (IO%LTR_ML) THEN
  ALLOCATE(IR%XFAPARC   (KI, IO%NPATCH))
  ALLOCATE(IR%XFAPIRC   (KI, IO%NPATCH))
  ALLOCATE(IR%XLAI_EFFC (KI, IO%NPATCH))
  ALLOCATE(IR%XMUS      (KI, IO%NPATCH))
  IR%XFAPARC   (:,:) = 0.
  IR%XFAPIRC   (:,:) = 0.
  IR%XLAI_EFFC (:,:) = 0.
  IR%XMUS      (:,:) = 0.
ELSE
  ALLOCATE(IR%XFAPARC   (0,0))
  ALLOCATE(IR%XFAPIRC   (0,0))
  ALLOCATE(IR%XLAI_EFFC (0,0))
  ALLOCATE(IR%XMUS      (0,0))
ENDIF        
!
!-------------------------------------------------------------------------------
!
!* albedo per tile and averaged albedo, emissivity and radiative temperature
!
ALLOCATE(MA%XALBNIR_SOIL(KI,IO%NPATCH))
ALLOCATE(MA%XALBVIS_SOIL(KI,IO%NPATCH))
ALLOCATE(MA%XALBUV_SOIL (KI,IO%NPATCH))
ALLOCATE(MT%XALBNIR     (KI,IO%NPATCH))
ALLOCATE(MT%XALBVIS     (KI,IO%NPATCH))
ALLOCATE(MT%XALBUV      (KI,IO%NPATCH))
MA%XALBNIR_SOIL(:,:) = XUNDEF
MA%XALBVIS_SOIL(:,:) = XUNDEF
MA%XALBUV_SOIL (:,:) = XUNDEF
MT%XALBNIR     (:,:) = XUNDEF
MT%XALBVIS     (:,:) = XUNDEF
MT%XALBUV      (:,:) = XUNDEF
!
OSURF_DIAG_ALBEDO = .TRUE.
!
!* Initialization of total albedo, emissivity and snow/flood fractions
!
ALLOCATE(IR%XPSN (KI,IO%NPATCH))
ALLOCATE(IR%XPSNG(KI,IO%NPATCH))
ALLOCATE(IR%XPSNV(KI,IO%NPATCH))
IR%XPSN  = 0.0
IR%XPSNG = 0.0
IR%XPSNV = 0.0
!
IF(IR%TSNOW%SCHEME=='EBA')THEN
   ALLOCATE(IR%XPSNV_A(KI,IO%NPATCH))
   IR%XPSNV_A = 0.0
ELSE
   ALLOCATE(IR%XPSNV_A(0,0))
ENDIF
!
PDIR_ALB = XUNDEF
PSCA_ALB = XUNDEF
PEMIS_OUT= XUNDEF
PTSRAD   = XUNDEF
!
IF (LHOOK) CALL DR_HOOK('INIT_VEG_n',1,ZHOOK_HANDLE)
!
END SUBROUTINE INIT_VEG_n
