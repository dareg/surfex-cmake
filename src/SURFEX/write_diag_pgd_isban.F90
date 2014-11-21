!     #########
      SUBROUTINE WRITE_DIAG_PGD_ISBA_n(HPROGRAM)
!     #########################################
!
!!****  *WRITE_DIAG_PGD_ISBA_n* - writes the ISBA physiographic diagnostic fields
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
!!      Modified    11/2013 by B. Decharme : XPATCH now in writesurf_isban.F90
!!      Modified    10/2014 by P. Samuelsson: MEB variables
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_PAR,   ONLY : XUNDEF, NUNDEF
USE MODD_ISBA_n,     ONLY : NPATCH, CPHOTO, CHORT, CISBA, XPATCH,                   &
                              XLAI, XVEG, XZ0,XALBNIR_SOIL,XALBVIS_SOIL,XALBUV_SOIL,&
                              XRSMIN, XGAMMA, XRGL, XCV, XEMIS, XDG, XWRMAX_CF,     &
                              XZ0REL, XVEGTYPE_PATCH, XALBNIR, XALBVIS, XALBUV,     &
                              XWATSUP, TSEED, TREAP, XIRRIG, XD_ICE,                &
                              XROOTFRAC, NWG_LAYER, XDROOT, XDG2,                   &
                              XWSAT, XWFC, XWWILT, XRUNOFFD, LSOC, XFRACSOC,        &
                              LMEB_PATCH,                                           &
                              XGNDLITTER,XZF_TALLVEG, XRGLGV,XGAMMAGV,              &
                              XRSMINGV, XWRMAX_CFGV,                                &
                              XH_VEG, XLAIGV, XZ0GV, XROOTFRACGV
USE MODD_AGRI,       ONLY : LAGRIP
!
USE MODD_DIAG_MISC_ISBA_n,ONLY : LSURF_DIAG_ALBEDO
!
USE MODD_IO_SURF_FA, ONLY : LFANOCOMPACT, LPREP
!
USE MODD_CH_ISBA_n,  ONLY : XSOILRC_SO2, XSOILRC_O3, CCH_DRY_DEP, NBEQ
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
REAL, DIMENSION(SIZE(XDG,1),SIZE(XDG,3)) :: ZWORK ! Work array
REAL, DIMENSION(SIZE(XDG,1),SIZE(XDG,2)) :: ZDG   ! Work array
REAL, DIMENSION(SIZE(XDG,1)            ) :: ZDG2
REAL, DIMENSION(SIZE(XDG,1)            ) :: ZDTOT
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=100):: YCOMMENT       ! Comment string
CHARACTER(LEN=2)  :: YLVLV, YPAS
CHARACTER(LEN=4)  :: YLVL
!
INTEGER         :: JJ, JL, JP, ILAYER
INTEGER           :: ISIZE_LMEB_PATCH   ! Number of patches where multi-energy balance should be applied
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!         Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_PGD_ISBA_N',0,ZHOOK_HANDLE)
!
ISIZE_LMEB_PATCH=COUNT(LMEB_PATCH(:))
!
 CALL INIT_IO_SURF_n(HPROGRAM,'NATURE','ISBA  ','WRITE')
!
!-------------------------------------------------------------------------------
!
!* Leaf Area Index
!
IF (CPHOTO=='NON' .OR. CPHOTO=='AGS' .OR. CPHOTO=='AST') THEN
  !
  YRECFM='LAI'
  YCOMMENT='leaf area index (-)'
  !
  CALL WRITE_SURF(HPROGRAM,YRECFM,XLAI(:,:),IRESP,HCOMMENT=YCOMMENT)
  !
  IF (ISIZE_LMEB_PATCH>0) THEN
    !
    YRECFM='LAIGV'
    YCOMMENT='MEB: understory leaf area index (-)'
    !
    CALL WRITE_SURF(HPROGRAM,YRECFM,XLAIGV(:,:),IRESP,HCOMMENT=YCOMMENT)
    !
  ENDIF
  !
ENDIF
!
!-------------------------------------------------------------------------------
!
!* Vegetation fraction
!
YRECFM='VEG'
YCOMMENT='vegetation fraction (-)'
!
 CALL WRITE_SURF(HPROGRAM,YRECFM,XVEG(:,:),IRESP,HCOMMENT=YCOMMENT)
!
!* Surface roughness length (without snow)
!
YRECFM='Z0VEG'
YCOMMENT='surface roughness length (without snow) (m)'
!
 CALL WRITE_SURF(HPROGRAM,YRECFM,XZ0(:,:),IRESP,HCOMMENT=YCOMMENT)
!
IF (ISIZE_LMEB_PATCH>0) THEN
  !
  YRECFM='GNDLITTER'
  YCOMMENT='MEB: ground litter fraction (-)'
  !
  CALL WRITE_SURF(HPROGRAM,YRECFM,XGNDLITTER(:,:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='Z0LITTER'
  YCOMMENT='MEB: ground litter roughness length (without snow) (m)'
  !
  CALL WRITE_SURF(HPROGRAM,YRECFM,XZ0GV(:,:),IRESP,HCOMMENT=YCOMMENT)
  !
ENDIF
!
!-------------------------------------------------------------------------------
!
!* Soil depth for each patch
!
DO JL=1,SIZE(XDG,2)
  IF (JL<10) THEN
    WRITE(YRECFM,FMT='(A2,I1)') 'DG',JL
  ELSE
    WRITE(YRECFM,FMT='(A2,I2)') 'DG',JL          
  ENDIF
  YCOMMENT='soil depth'//' (M)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XDG(:,JL,:),IRESP,HCOMMENT=YCOMMENT)
END DO
!
!* Averaged Soil depth
!
IF(NPATCH>1)THEN
!        
  ZDG(:,:)=0.0
  DO JP=1,NPATCH
     DO JL=1,SIZE(XDG,2)
        DO JJ=1,SIZE(XDG,1) 
           ZDG(JJ,JL)=ZDG(JJ,JL)+XPATCH(JJ,JP)*XDG(JJ,JL,JP)
        ENDDO
     ENDDO
  ENDDO
!
  DO JL=1,SIZE(XDG,2)
    WRITE(YLVL,'(I4)')JL
    YRECFM='DG'//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))
    YRECFM=YRECFM(:LEN_TRIM(YRECFM))//'_ISBA'
    YCOMMENT='averaged soil depth layer '//ADJUSTL(YLVL(:LEN_TRIM(YLVL)))//' (m)'
    CALL WRITE_SURF(HPROGRAM,YRECFM,ZDG(:,JL),IRESP,HCOMMENT=YCOMMENT)
  END DO
!        
ENDIF
!
!-------------------------------------------------------------------------------
!
IF(CISBA=='DIF')THEN
!
  ZDG2 (:)=0.0
  ZDTOT(:)=0.0
  ZWORK(:,:)=XUNDEF
  DO JP=1,SIZE(XDG,3)
     DO JJ=1,SIZE(XDG,1)
        ZDG2(JJ)=ZDG2(JJ)+XPATCH(JJ,JP)*XDG2(JJ,JP)
        JL=NWG_LAYER(JJ,JP)
        IF(JL/=NUNDEF)THEN
          ZWORK(JJ,JP)=XDG(JJ,JL,JP)
          ZDTOT(JJ)=ZDTOT(JJ)+XPATCH(JJ,JP)*XDG(JJ,JL,JP)
        ENDIF
     ENDDO
  ENDDO
!
!* Root depth
!
  YRECFM='DROOT_DIF'
  YCOMMENT='Root depth in ISBA-DIF'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XDROOT(:,:),IRESP,HCOMMENT=YCOMMENT)
!
  YRECFM='DG2_DIF'
  YCOMMENT='DG2 depth in ISBA-DIF'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XDG2(:,:),IRESP,HCOMMENT=YCOMMENT)
!  
  IF(NPATCH>1)THEN
    YRECFM='DG2_DIF_ISBA'
    YCOMMENT='Averaged DG2 depth in ISBA-DIF'
    CALL WRITE_SURF(HPROGRAM,YRECFM,ZDG2(:),IRESP,HCOMMENT=YCOMMENT)          
  ENDIF  
!
!* Runoff depth
!
  YRECFM='RUNOFFD'
  YCOMMENT='Runoff deph in ISBA-DIF'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XRUNOFFD(:,:),IRESP,HCOMMENT=YCOMMENT)
!
!* Total soil depth for mositure
!
  YRECFM='DTOT_DIF'
  YCOMMENT='Total soil depth for moisture in ISBA-DIF'
  CALL WRITE_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP,HCOMMENT=YCOMMENT)
!
  IF(NPATCH>1)THEN
    YRECFM='DTOTDF_ISBA'
    YCOMMENT='Averaged Total soil depth for moisture in ISBA-DIF'
    CALL WRITE_SURF(HPROGRAM,YRECFM,ZDTOT(:),IRESP,HCOMMENT=YCOMMENT)          
  ENDIF
!
!* Root fraction for each patch
!
  DO JL=1,SIZE(XROOTFRAC,2)
     IF (JL<10) THEN
       WRITE(YRECFM,FMT='(A8,I1)') 'ROOTFRAC',JL
     ELSE
       WRITE(YRECFM,FMT='(A8,I2)') 'ROOTFRAC',JL          
     ENDIF  
     YCOMMENT='root fraction by layer (-)'
     ZWORK(:,:)=XUNDEF
     DO JJ=1,SIZE(XDG,1)
        WHERE(JL<=NWG_LAYER(JJ,:).AND.NWG_LAYER(JJ,:)/=NUNDEF)
              ZWORK(JJ,:)=XROOTFRAC(JJ,JL,:)
        ENDWHERE
     ENDDO
     CALL WRITE_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP,HCOMMENT=YCOMMENT)
  END DO
  !
  IF (ISIZE_LMEB_PATCH>0) THEN
    DO JL=1,SIZE(XROOTFRACGV,2)
       IF (JL<10) THEN
         WRITE(YRECFM,FMT='(A10,I1)') 'ROOTFRACGV',JL
       ELSE
         WRITE(YRECFM,FMT='(A10,I2)') 'ROOTFRACGV',JL          
       ENDIF  
       YCOMMENT='MEB: understory root fraction by layer (-)'
       ZWORK(:,:)=XUNDEF
       DO JJ=1,SIZE(XDG,1)
          WHERE(JL<=NWG_LAYER(JJ,:).AND.NWG_LAYER(JJ,:)/=NUNDEF)
                ZWORK(JJ,:)=XROOTFRACGV(JJ,JL,:)
          ENDWHERE
       ENDDO
       CALL WRITE_SURF(HPROGRAM,YRECFM,ZWORK(:,:),IRESP,HCOMMENT=YCOMMENT)
    END DO
  ENDIF
!
!* SOC fraction for each layer
!
  IF(LSOC)THEN
    DO JL=1,SIZE(XDG,2)
     IF (JL<10) THEN
       WRITE(YRECFM,FMT='(A7,I1)') 'FRACSOC',JL
     ELSE
       WRITE(YRECFM,FMT='(A7,I2)') 'FRACSOC',JL          
     ENDIF  
     YCOMMENT='SOC fraction by layer (-)'
     CALL WRITE_SURF(HPROGRAM,YRECFM,XFRACSOC(:,JL),IRESP,HCOMMENT=YCOMMENT)
    END DO
  ENDIF
!
ENDIF        
!
!-------------------------------------------------------------------------------
!
DO JL=1,SIZE(XDG,2)
   IF (JL<10) THEN
     WRITE(YRECFM,FMT='(A4,I1)') 'WSAT',JL
   ELSE
     WRITE(YRECFM,FMT='(A4,I2)') 'WSAT',JL          
   ENDIF  
  YCOMMENT='soil porosity by layer (m3/m3)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XWSAT(:,JL),IRESP,HCOMMENT=YCOMMENT)
ENDDO
!
DO JL=1,SIZE(XDG,2)
   IF (JL<10) THEN
     WRITE(YRECFM,FMT='(A3,I1)') 'WFC',JL
   ELSE
     WRITE(YRECFM,FMT='(A3,I2)') 'WFC',JL          
   ENDIF  
  YCOMMENT='field capacity by layer (m3/m3)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XWFC(:,JL),IRESP,HCOMMENT=YCOMMENT)
ENDDO
!
DO JL=1,SIZE(XDG,2)
   IF (JL<10) THEN
     WRITE(YRECFM,FMT='(A5,I1)') 'WWILT',JL
   ELSE
     WRITE(YRECFM,FMT='(A5,I2)') 'WWILT',JL          
   ENDIF  
  YCOMMENT='wilting point by layer (m3/m3)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XWWILT(:,JL),IRESP,HCOMMENT=YCOMMENT)
ENDDO     
!
!-------------------------------------------------------------------------------
! For Earth System Model
IF(LFANOCOMPACT.AND..NOT.LPREP)THEN
  CALL END_IO_SURF_n(HPROGRAM)
  IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_PGD_ISBA_N',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!
!-------------------------------------------------------------------------------
!
YRECFM='Z0REL'
YCOMMENT='orography roughness length (M)'
!
 CALL WRITE_SURF(HPROGRAM,YRECFM,XZ0REL(:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
!* Runoff soil ice depth for each patch
!
IF(CHORT=='SGH'.AND.CISBA/='DIF')THEN
  YRECFM='DICE'
  YCOMMENT='soil ice depth for runoff (m)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XD_ICE(:,:),IRESP,HCOMMENT=YCOMMENT)
ENDIF
!
!-------------------------------------------------------------------------------
!
!* Fraction of each vegetation type for each patch
!
DO JL=1,SIZE(XVEGTYPE_PATCH,2)
  WRITE(YPAS,'(I2)') JL 
  YLVLV=ADJUSTL(YPAS(:LEN_TRIM(YPAS)))
  WRITE(YRECFM,FMT='(A9)') 'VEGTY_P'//YLVLV
  YCOMMENT='fraction of each vegetation type for each patch'//' (-)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XVEGTYPE_PATCH(:,JL,:),IRESP,HCOMMENT=YCOMMENT)
END DO
!-------------------------------------------------------------------------------
!
!* other surface parameters
!
YRECFM='RSMIN'
YCOMMENT='minimum stomatal resistance (sm-1)'
 CALL WRITE_SURF(HPROGRAM,YRECFM,XRSMIN(:,:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='GAMMA'
YCOMMENT='coefficient for RSMIN calculation (-)'
 CALL WRITE_SURF(HPROGRAM,YRECFM,XGAMMA(:,:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='CV'
YCOMMENT='vegetation thermal inertia coefficient (-)'
 CALL WRITE_SURF(HPROGRAM,YRECFM,XCV(:,:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='RGL'
YCOMMENT='maximum solar radiation usable in photosynthesis (-)'
 CALL WRITE_SURF(HPROGRAM,YRECFM,XRGL(:,:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='EMIS_ISBA'
YCOMMENT='surface emissivity (-)'
 CALL WRITE_SURF(HPROGRAM,YRECFM,XEMIS(:,:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='WRMAX_CF'
YCOMMENT='coefficient for maximum water interception (-)'
 CALL WRITE_SURF(HPROGRAM,YRECFM,XWRMAX_CF(:,:),IRESP,HCOMMENT=YCOMMENT)
!
IF (ISIZE_LMEB_PATCH>0) THEN
  !
  YRECFM='RSMINGV'
  YCOMMENT='MEB: understory minimum stomatal resistance (sm-1)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XRSMINGV(:,:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='GAMMAGV'
  YCOMMENT='MEB: understory coefficient for RSMIN calculation (-)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XGAMMAGV(:,:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='RGLGV'
  YCOMMENT='MEB: understory maximum solar radiation usable in photosynthesis (-)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XRGLGV(:,:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='WRMAX_CFGV'
  YCOMMENT='MEB: understory coefficient for maximum water interception (-)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XWRMAX_CFGV(:,:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='ZF_TALLVEG'
  YCOMMENT='MEB: identification variable for tall vegetation (-)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XZF_TALLVEG(:,:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='H_VEG'
  YCOMMENT='MEB: height of vegetation (m)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XH_VEG(:,:),IRESP,HCOMMENT=YCOMMENT)
  !
ENDIF
!
!-------------------------------------------------------------------------------
!
IF (LSURF_DIAG_ALBEDO) THEN
!
!* Soil albedos
!
!
   YRECFM='ALBNIR_S'
   YCOMMENT='soil near-infra-red albedo (-)'
   CALL WRITE_SURF(HPROGRAM,YRECFM,XALBNIR_SOIL(:,:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
   YRECFM='ALBVIS_S'
   YCOMMENT='soil visible albedo (-)'
   CALL WRITE_SURF(HPROGRAM,YRECFM,XALBVIS_SOIL(:,:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
   YRECFM='ALBUV_S'
   YCOMMENT='soil UV albedo (-)'
   CALL WRITE_SURF(HPROGRAM,YRECFM,XALBUV_SOIL(:,:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
!* albedos
!
   YRECFM='ALBNIR_ISBA'
   YCOMMENT='total near-infra-red albedo (-)'
   CALL WRITE_SURF(HPROGRAM,YRECFM,XALBNIR(:,:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
   YRECFM='ALBVIS_ISBA'
   YCOMMENT='total visible albedo (-)'
   CALL WRITE_SURF(HPROGRAM,YRECFM,XALBVIS(:,:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
   YRECFM='ALBUV_ISBA'
   YCOMMENT='total UV albedo (-)'
   CALL WRITE_SURF(HPROGRAM,YRECFM,XALBUV(:,:),IRESP,HCOMMENT=YCOMMENT)
!
END IF
!
!-------------------------------------------------------------------------------
!
!* chemical soil resistances
!
IF (CCH_DRY_DEP=='WES89' .AND. NBEQ>0) THEN
  YRECFM='SOILRC_SO2'
  YCOMMENT='bare soil resistance for SO2 (?)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XSOILRC_SO2(:,:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='SOILRC_O3'
  YCOMMENT='bare soil resistance for O3 (?)'
  CALL WRITE_SURF(HPROGRAM,YRECFM,XSOILRC_O3(:,:),IRESP,HCOMMENT=YCOMMENT)
END IF
!
!-------------------------------------------------------------------------------
!
IF (LAGRIP .AND. (CPHOTO=='LAI' .OR. CPHOTO=='LST' .OR. CPHOTO=='NIT' .OR. CPHOTO=='NCB') ) THEN
!
!* seeding and reaping
!
!
  YRECFM='TSEED'
  YCOMMENT='date of seeding (-)'
!
  CALL WRITE_SURF(HPROGRAM,YRECFM,TSEED(:,:),IRESP,HCOMMENT=YCOMMENT)
!
  YRECFM='TREAP'
  YCOMMENT='date of reaping (-)'
!
  CALL WRITE_SURF(HPROGRAM,YRECFM,TREAP(:,:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
!* irrigated fraction
!
  YRECFM='IRRIG'
  YCOMMENT='flag for irrigation (irrigation if >0.) (-)'
!
  CALL WRITE_SURF(HPROGRAM,YRECFM,XIRRIG(:,:),IRESP,HCOMMENT=YCOMMENT)
!
!-------------------------------------------------------------------------------
!
!* water supply for irrigation
!
  YRECFM='WATSUP'
  YCOMMENT='water supply during irrigation process (mm)'
!
  CALL WRITE_SURF(HPROGRAM,YRECFM,XWATSUP(:,:),IRESP,HCOMMENT=YCOMMENT)
!
ENDIF
!-------------------------------------------------------------------------------
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_PGD_ISBA_N',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE WRITE_DIAG_PGD_ISBA_n
