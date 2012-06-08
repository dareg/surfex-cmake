!     #########
      SUBROUTINE WRITE_DIAG_PGD_GRDN_n(HPROGRAM)
!     #########################################
!
!!****  *WRITE_DIAG_PGD_TEB_GARDEN_n* - writes the ISBA physiographic diagnostic fields
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
!!      Modified    10/2004 by P. Le Moigne: add XZ0REL, XVEGTYPE_PATCH
!!      Modified    11/2005 by P. Le Moigne: limit length of VEGTYPE_PATCH field names
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_TEB_n,          ONLY : XGARDEN
USE MODD_TEB_GARDEN_n,   ONLY : NPATCH, CPHOTO, CHORT,          &
                                  XLAI, XVEG, XZ0,XALBNIR_SOIL,XALBVIS_SOIL,XALBUV_SOIL,&
                                  XRSMIN, XGAMMA, XRGL, XCV, XEMIS, XDG, XWRMAX_CF,     &
                                  XZ0REL, XVEGTYPE_PATCH, XALBNIR, XALBVIS, XALBUV,     &
                                  XPATCH, XWATSUP, TSEED, TREAP, XIRRIG, XD_ICE  
!
USE MODD_DIAG_MISC_TEB_n,ONLY : LSURF_DIAG_ALBEDO
!
USE MODD_IO_SURF_FA, ONLY : LFANOCOMPACT, LPREP
!
USE MODI_INIT_IO_SURF_n
USE MODI_WRITE_SURF
USE MODI_END_IO_SURF_n
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
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
CHARACTER(LEN=16) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=100):: YCOMMENT       ! Comment string
CHARACTER(LEN=2)  :: YLVLV, YPAS
!
INTEGER           :: JL, JP
REAL, DIMENSION(:,:), ALLOCATABLE :: ZWORK
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!         Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_PGD_GRDN_N',0,ZHOOK_HANDLE)
CALL INIT_IO_SURF_n(HPROGRAM,'TOWN  ','TEB   ','WRITE')
!
YRECFM='FRAC_GARDEN'
CALL WRITE_SURF(HPROGRAM,YRECFM,XGARDEN(:),IRESP,HCOMMENT=YCOMMENT)
!
!* Leaf Area Index
!
IF (CPHOTO=='NON' .OR. CPHOTO=='AGS' .OR. CPHOTO=='AST') THEN
  !
  YRECFM='G_LAI'
  YCOMMENT='leaf area index (-)'
  !
  CALL WRITE_SURF(HPROGRAM,YRECFM,XLAI(:,:),IRESP,HCOMMENT=YCOMMENT)
  !
ENDIF
!
!-------------------------------------------------------------------------------
!
!* Vegetation fraction
!
YRECFM='G_VEG'
YCOMMENT='vegetation fraction (-)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XVEG(:,:),IRESP,HCOMMENT=YCOMMENT)
!
!* Surface roughness length (without snow)
!
YRECFM='G_Z0VEG'
YCOMMENT='surface roughness length (without snow) (M)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XZ0(:,:),IRESP,HCOMMENT=YCOMMENT)
!
!* Fraction for each patch
!
YRECFM='G_PATCH'
YCOMMENT='fraction for each patch (-)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XPATCH(:,:),IRESP,HCOMMENT=YCOMMENT)
!-------------------------------------------------------------------------------
!
!* Soil depth for each patch
!
ALLOCATE(ZWORK(SIZE(XDG,1),SIZE(XDG,3)))
DO JL=1,SIZE(XDG,2)
  WRITE(YRECFM,FMT='(A4,I1)') 'G_DG',JL
  YCOMMENT='soil depth'//' (M)'
  DO JP=1,SIZE(XDG,3)
    ZWORK(:,JP) = XDG(:,JL,JP)
  END DO
  CALL WRITE_SURF(HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
END DO
DEALLOCATE(ZWORK)
!
!-------------------------------------------------------------------------------
! For Earth System Model
IF(LFANOCOMPACT.AND..NOT.LPREP)THEN
   CALL END_IO_SURF_n(HPROGRAM)
   RETURN
ENDIF
!
!-------------------------------------------------------------------------------
!
!* Runoff soil ice depth for each patch
!
IF(CHORT=='SGH')THEN
  YRECFM='G_DICE'
  YCOMMENT='soil ice depth for runoff (m)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XD_ICE(:,:),IRESP,HCOMMENT=YCOMMENT)
ENDIF
!
!-------------------------------------------------------------------------------
!
YRECFM='G_Z0REL'
YCOMMENT='orography roughness length (M)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XZ0REL(:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
!* Fraction of each vegetation type for each patch
!
ALLOCATE(ZWORK(SIZE(XVEGTYPE_PATCH,1),SIZE(XVEGTYPE_PATCH,3)))
DO JL=1,SIZE(XVEGTYPE_PATCH,2)
  WRITE(YPAS,'(I2)') JL 
  YLVLV=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
  WRITE(YRECFM,FMT='(A11)') 'G_VEGTYPE_P'//YLVLV
  YCOMMENT='fraction of each vegetation type for each patch'//' (-)'
  DO JP=1,SIZE(XVEGTYPE_PATCH,3)
    ZWORK(:,JP) = XVEGTYPE_PATCH(:,JL,JP)
  END DO
  CALL WRITE_SURF(HPROGRAM,YRECFM,ZWORK,IRESP,HCOMMENT=YCOMMENT)
END DO
DEALLOCATE(ZWORK)
!-------------------------------------------------------------------------------
!
!* other surface parameters
!
YRECFM='G_RSMIN'
YCOMMENT='minimum stomatal resistance (SM-1)'
CALL WRITE_SURF(HPROGRAM,YRECFM,XRSMIN(:,:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='G_GAMMA'
YCOMMENT='coefficient for RSMIN calculation (-)'
CALL WRITE_SURF(HPROGRAM,YRECFM,XGAMMA(:,:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='G_CV'
YCOMMENT='vegetation thermal inertia coefficient (-)'
CALL WRITE_SURF(HPROGRAM,YRECFM,XCV(:,:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='G_RGL'
YCOMMENT='maximum solar radiation usable in photosynthesis (-)'
CALL WRITE_SURF(HPROGRAM,YRECFM,XRGL(:,:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='G_EMIS_ISBA'
YCOMMENT='surface emissivity (-)'
CALL WRITE_SURF(HPROGRAM,YRECFM,XEMIS(:,:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='G_WRMAX_CF'
YCOMMENT='coefficient for maximum water interception (-)'
CALL WRITE_SURF(HPROGRAM,YRECFM,XWRMAX_CF(:,:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
IF (LSURF_DIAG_ALBEDO) THEN
!
!* Soil albedos
!
!
   YRECFM='GALBNIR_SOIL'
   YCOMMENT='soil near-infra-red albedo (-)'
   CALL WRITE_SURF(HPROGRAM,YRECFM,XALBNIR_SOIL(:,:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
   YRECFM='GALBVIS_SOIL'
   YCOMMENT='soil visible albedo (-)'
   CALL WRITE_SURF(HPROGRAM,YRECFM,XALBVIS_SOIL(:,:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
   YRECFM='GALBUV_SOIL'
   YCOMMENT='soil UV albedo (-)'
   CALL WRITE_SURF(HPROGRAM,YRECFM,XALBUV_SOIL(:,:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
!* albedos
!
   YRECFM='GALBNIR_ISBA'
   YCOMMENT='total near-infra-red albedo (-)'
   CALL WRITE_SURF(HPROGRAM,YRECFM,XALBNIR(:,:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
   YRECFM='GALBVIS_ISBA'
   YCOMMENT='total visible albedo (-)'
   CALL WRITE_SURF(HPROGRAM,YRECFM,XALBVIS(:,:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
   YRECFM='GALBUV_ISBA'
   YCOMMENT='total UV albedo (-)'
   CALL WRITE_SURF(HPROGRAM,YRECFM,XALBUV(:,:),IRESP,HCOMMENT=YCOMMENT)
!
END IF
!
!-------------------------------------------------------------------------------
!
!         End of IO
!
CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_PGD_GRDN_N',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE WRITE_DIAG_PGD_GRDN_n
