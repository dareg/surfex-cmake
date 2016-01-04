!     #########
SUBROUTINE PREP_ISBA (DTCO, ICP, IG, I, UG, U, USS, &
                      HPROGRAM,HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!     #################################################################################
!
!!****  *PREP_ISBA* - Prepares ISBA fields
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!     V. Masson 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      Modified by P. Le Moigne (11/2004): AGS fields
!!      Modified by B. Decharme   (2008)  : Floodplains
!!      Modified by B. Decharme  (01/2009): Consistency with Arpege deep soil
!!                                          temperature
!!      Modified by B. Decharme  (03/2009): Consistency with Arpege permanent
!!                                          snow/ice treatment
!!      A.L. Gibelin 04/2009 : BIOMASS and RESP_BIOMASS arrays 
!!      A.L. Gibelin 06/2009 : Soil carbon variables for CNT option
!!      Modified by S. Riette    (06/2009): PREP_ISBA_CANOPY has no more arg.
!!      Modified by S. Riette    (04/2010): ecmwf ice content is computed during
!!                                          grib reading (no longer here)
!!      B. Decharme  (10/2012): coherence between soil temp and liquid/solid water with DIF
!!                              bug in biomass prognostic fields calculation
!!      B. Decharme  (06/2013): XPSNV_A for EBA snow scheme not allocated
!!      M. Lafaysse (04/2014) : LSNOW_PREP_PERM
!!      B. Decharme  (04/2013): Good computation for coherence between soil temp and 
!!                              liquid/solid water with DIF (results don't change)
!!                              if lglacier in input file, do not initialize again
!!      P. Samuelsson            (10/2014): MEB
!!------------------------------------------------------------------
!
!
!
!
USE MODD_SURFEX_MPI, ONLY : NRANK
!
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_ISBA_CANOPY_n, ONLY : ISBA_CANOPY_t
USE MODD_ISBA_GRID_n, ONLY : ISBA_GRID_t
USE MODD_ISBA_n, ONLY : ISBA_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SURF_ATM_SSO_n, ONLY : SURF_ATM_SSO_t
!
USE MODI_PREP_HOR_ISBA_FIELD
USE MODI_PREP_VER_ISBA
USE MODI_PREP_OUTPUT_GRID
USE MODI_GET_LUOUT
USE MODI_PREP_ISBA_CANOPY
!
USE MODD_READ_NAMELIST,  ONLY : LNAM_READ
USE MODD_SURF_ATM,       ONLY : LVERTSHIFT
USE MODD_DATA_COVER_PAR, ONLY : NVT_SNOW
!
!                           
USE MODD_DEEPSOIL,    ONLY : LPHYSDOMC
USE MODD_CSTS,        ONLY : XTT, XG, XLMTT
USE MODD_SNOW_PAR,    ONLY : XEMISSN
USE MODD_ISBA_PAR,    ONLY : XWGMIN
!
USE MODD_CO2V_PAR,    ONLY : XANFMINIT
USE MODD_SURF_PAR,    ONLY : XUNDEF
USE MODD_PREP,        ONLY : XZS_LS

USE MODD_PREP_SNOW,   ONLY : LSNOW_PREP_PERM
!
USE MODN_PREP_ISBA
USE MODN_PREP_ISBA_SNOW, ONLY : LSWEMAX, XSWEMAX 
!
USE MODI_VEGTYPE_TO_PATCH
USE MODI_PREP_PERM_SNOW
USE MODI_INIT_SNOW_LW
USE MODI_AVERAGED_ALBEDO_EMIS_ISBA
USE MODI_PREP_HOR_ISBA_CC_FIELD
USE MODI_SOIL_ALBEDO
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_CLEAN_PREP_OUTPUT_GRID
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(ISBA_CANOPY_t), INTENT(INOUT) :: ICP
TYPE(ISBA_GRID_t), INTENT(INOUT) :: IG
TYPE(ISBA_t), INTENT(INOUT) :: I
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SURF_ATM_SSO_t), INTENT(INOUT) :: USS
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=28),  INTENT(IN)  :: HATMFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),   INTENT(IN)  :: HATMFILETYPE! type of the Atmospheric file
 CHARACTER(LEN=28),  INTENT(IN)  :: HPGDFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),   INTENT(IN)  :: HPGDFILETYPE! type of the Atmospheric file
!
!*      0.2    declarations of local variables
!
INTEGER :: ILUOUT, INI
INTEGER :: JP, JL, JJ
INTEGER :: ISNOW          ! patch number where permanent snow is
REAL    :: ZWORK, ZLOG, ZWTOT, ZMATPOT, ZWL
!
REAL,             DIMENSION(1)   :: ZSW_BANDS ! middle wavelength of each band
REAL,             DIMENSION(SIZE(I%M%T%XLAI,1),SIZE(I%M%T%XLAI,2)) :: ZDIR_ALB  ! direct albedo for each band
REAL,             DIMENSION(SIZE(I%M%T%XLAI,1),SIZE(I%M%T%XLAI,2)) :: ZSCA_ALB  ! diffuse albedo for each band
REAL,             DIMENSION(SIZE(I%M%T%XLAI,1))   :: ZEMIS     ! emissivity
REAL,             DIMENSION(SIZE(I%M%T%XLAI,1))   :: ZZENITH   ! solar zenithal angle
REAL,             DIMENSION(SIZE(I%M%T%XLAI,1))   :: ZTSURF     ! surface effective temperature
!
LOGICAL         :: GPERMSNOW
LOGICAL         :: GTEMP2WGI
LOGICAL         :: GWG
LOGICAL         :: GWGI
LOGICAL         :: GTG
!
REAL            :: SMAX
!
INTEGER         :: ISIZE_LMEB_PATCH
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PREP_ISBA',0,ZHOOK_HANDLE)
!
!*      1.     Default of configuration
!
GPERMSNOW = .TRUE.
GWG       = .TRUE.
GWGI      = .TRUE.
GTG       = .TRUE.
!
ISIZE_LMEB_PATCH=COUNT(I%O%LMEB_PATCH(:))
!
!*      1.1    Default
!
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
 CALL PREP_OUTPUT_GRID(UG, U, &
                       ILUOUT,IG%CGRID,IG%XGRID_PAR,IG%XLAT,IG%XLON)
!
!-------------------------------------------------------------------------------------
!
!*      2.     Reading and horizontal interpolations
!
!
!*      2.0    Large scale orography
!
 CALL PREP_HOR_ISBA_FIELD(DTCO, IG, I, UG, U, USS, &
                          HPROGRAM,'ZS     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!*      2.1    Soil Water reservoirs
!
 CALL PREP_HOR_ISBA_FIELD(DTCO, IG, I, UG, U, USS, &
                          HPROGRAM,'WG     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,GWG)
!
!*      2.2    Soil ice reservoirs
!
 CALL PREP_HOR_ISBA_FIELD(DTCO, IG, I, UG, U, USS, &
                          HPROGRAM,'WGI    ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,GWGI)
!
!*      2.3    Leaves interception water reservoir
!
 CALL PREP_HOR_ISBA_FIELD(DTCO, IG, I, UG, U, USS, &
                          HPROGRAM,'WR     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!*      2.4    Temperature profile
!
 CALL PREP_HOR_ISBA_FIELD(DTCO, IG, I, UG, U, USS, &
                          HPROGRAM,'TG     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,GTG)
!
!*      2.5    Snow variables
!
 CALL PREP_HOR_ISBA_FIELD(DTCO, IG, I, UG, U, USS, &
                          HPROGRAM,'SN_VEG ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,GPERMSNOW)
!
!*      2.6    LAI
!
 CALL PREP_HOR_ISBA_FIELD(DTCO, IG, I, UG, U, USS, &
                          HPROGRAM,'LAI    ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!*      2.7    GLACIER
!
IF(I%O%LGLACIER)THEN
  CALL PREP_HOR_ISBA_FIELD(DTCO, IG, I, UG, U, USS, &
                          HPROGRAM,'ICE_STO',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
ENDIF
!
!*      2.8    Canopy vegetation temperature and interception reservoirs and air variables
!
IF(ISIZE_LMEB_PATCH>0)THEN
  CALL PREP_HOR_ISBA_FIELD(DTCO, IG, I, UG, U, USS, &
                          HPROGRAM,'TV     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  CALL PREP_HOR_ISBA_FIELD(DTCO, IG, I, UG, U, USS, &
                        HPROGRAM,'TL     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  CALL PREP_HOR_ISBA_FIELD(DTCO, IG, I, UG, U, USS, &
                          HPROGRAM,'WRL    ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  CALL PREP_HOR_ISBA_FIELD(DTCO, IG, I, UG, U, USS, &
                          HPROGRAM,'WRLI   ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  CALL PREP_HOR_ISBA_FIELD(DTCO, IG, I, UG, U, USS, &
                          HPROGRAM,'WRVN   ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  CALL PREP_HOR_ISBA_FIELD(DTCO, IG, I, UG, U, USS, &
                          HPROGRAM,'TC     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  CALL PREP_HOR_ISBA_FIELD(DTCO, IG, I, UG, U, USS, &
                          HPROGRAM,'QC     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
ENDIF
!
!-------------------------------------------------------------------------------------
!
!*      3.    Physical limitation: 
!
! No ice for force restore third layer:
IF (I%O%CISBA == '3-L') THEN
   DO JP=1,I%O%NPATCH
      WHERE(I%R%XWG(:,3,JP) /= XUNDEF)
        I%R%XWG(:,3,JP)  = MIN(I%R%XWG(:,3,JP)+I%R%XWGI(:,3,JP),I%IP%XWSAT(:,3))
        I%R%XWGI(:,3,JP) = 0.
      END WHERE
   ENDDO
ENDIF
!
! Total water content should not exceed saturation:
DO JP=1,I%O%NPATCH
   WHERE(I%R%XWG(:,:,JP) /= XUNDEF .AND. (I%R%XWG(:,:,JP) + I%R%XWGI(:,:,JP)) > I%IP%XWSAT(:,:) )
      I%R%XWGI(:,:,JP) = I%IP%XWSAT(:,:) - I%R%XWG(:,:,JP)
   END WHERE
ENDDO
!
!-------------------------------------------------------------------------------------
!
!*      3.     Vertical interpolations of all variables
!
IF(LVERTSHIFT)THEN
  CALL PREP_VER_ISBA(I)
ENDIF
!
DEALLOCATE(XZS_LS)
!-------------------------------------------------------------------------------------
!
!*      4.     Treatment of permanent snow
!
IF (GPERMSNOW.AND.LSNOW_PREP_PERM) THEN
  ISNOW = VEGTYPE_TO_PATCH(NVT_SNOW,I%O%NPATCH)
  CALL PREP_PERM_SNOW(I, &
                      I%R%TSNOW,I%R%XTG(:,:,ISNOW),I%IP%XVEGTYPE_PATCH(:,:,ISNOW),ISNOW)
ENDIF
!
CALL INIT_SNOW_LW(XEMISSN,I%R%TSNOW)
!
IF (LPHYSDOMC) THEN
   I%R%TSNOW%WSNOW(:,:,:)=0.
ENDIF 
!------------------------------------------------------------------------------------- 
! 
!*      4.b     Possibility for setting an upper limit on the initial snow water equivalent field 
IF (LSWEMAX) THEN 
  SMAX = MAXVAL(I%R%TSNOW%WSNOW(:,:,:)) 
  WRITE(*,*) ' MAX(Snow content (kg/m2)): ', SMAX 
  WRITE(*,*) ' Set MAX to', XSWEMAX, '(kg/m2)' 
  I%R%TSNOW%WSNOW(:,:,:) = MIN(I%R%TSNOW%WSNOW(:,:,:),XSWEMAX) 
  SMAX = MAXVAL(I%R%TSNOW%WSNOW(:,:,:)) 
  WRITE(*,*) ' MAX(Snow content (kg/m2)): ', SMAX 
ENDIF
!
!-------------------------------------------------------------------------------------
!
!*      5.     coherence between soil temperature and liquid/solid water
!
GTEMP2WGI=(GWG.OR.GWGI.OR.GTG)
!
IF (I%O%CISBA == 'DIF'.AND.GTEMP2WGI) THEN
   INI=SIZE(I%IP%XWSAT,1)
   DO JP=1,I%O%NPATCH
      DO JL=1,I%O%NGROUND_LAYER
         DO JJ=1,INI
            IF(I%R%XWG(JJ,JL,JP)/=XUNDEF)THEN
!     
!             total soil moisture
              ZWTOT = I%R%XWG(JJ,JL,JP)+I%R%XWGI(JJ,JL,JP)
              ZWTOT = MIN(ZWTOT,I%IP%XWSAT(JJ,JL))
!              
!             total matric potential
!             psi=mpotsat*(w/wsat)**(-bcoef)
              ZWORK   = ZWTOT/I%IP%XWSAT(JJ,JL)
              ZLOG    = I%IP%XBCOEF(JJ,JL)*LOG(ZWORK)
              ZMATPOT = I%IP%XMPOTSAT(JJ,JL)*EXP(-ZLOG)
!
!             soil liquid water content computation
!             w=wsat*(psi/mpotsat)**(-1/bcoef)
              ZMATPOT       = MIN(I%IP%XMPOTSAT(JJ,JL),XLMTT*(I%R%XTG(JJ,JL,JP)-XTT)/(XG*I%R%XTG(JJ,JL,JP)))
              ZWORK         = MAX(1.0,ZMATPOT/I%IP%XMPOTSAT(JJ,JL))
              ZLOG          = LOG(ZWORK)
              ZWL           = I%IP%XWSAT(JJ,JL)*EXP(-ZLOG/I%IP%XBCOEF(JJ,JL))
              ZWL           = MAX(ZWL,XWGMIN)
              I%R%XWG(JJ,JL,JP) = MIN(ZWL,ZWTOT )
!        
!             soil ice computation    
              I%R%XWGI(JJ,JL,JP) = MAX(0.0,ZWTOT-I%R%XWG(JJ,JL,JP))
! 
!             supress numerical artefact
              IF(I%R%XTG(JJ,JL,JP)>=XTT)THEN
                I%R%XWG (JJ,JL,JP) = MIN(I%R%XWG(JJ,JL,JP)+I%R%XWGI(JJ,JL,JP),I%IP%XWSAT(JJ,JL))
                I%R%XWGI(JJ,JL,JP) = 0.0
              ENDIF
!
            ENDIF
        ENDDO        
      ENDDO        
   ENDDO
ENDIF
!
!-------------------------------------------------------------------------------------
!
!*      6.     Half prognostic fields
!              The only variable used from the AVERAGED_ALBEDO_EMIS_ISBA call
!              is XTSRAD_NAT. All other variables are treated as dummies.
!
ALLOCATE(I%R%XRESA(SIZE(I%M%T%XLAI,1),SIZE(I%M%T%XLAI,2)))
I%R%XRESA = 100.
!
ALLOCATE(I%R%XTSRAD_NAT(SIZE(I%M%T%XLAI,1)))
ZZENITH(:)=0.
ZSW_BANDS(:)=0.
!
ALLOCATE(I%M%T%XALBNIR(SIZE(I%M%T%XLAI,1),SIZE(I%M%T%XLAI,2)))
ALLOCATE(I%M%T%XALBVIS(SIZE(I%M%T%XLAI,1),SIZE(I%M%T%XLAI,2)))
ALLOCATE(I%M%T%XALBUV(SIZE(I%M%T%XLAI,1),SIZE(I%M%T%XLAI,2)))
I%M%T%XALBNIR = 0.0
I%M%T%XALBVIS = 0.0
I%M%T%XALBUV = 0.0
!
ALLOCATE(I%M%A%XALBNIR_SOIL(SIZE(I%M%T%XLAI,1),SIZE(I%M%T%XLAI,2)))
ALLOCATE(I%M%A%XALBVIS_SOIL(SIZE(I%M%T%XLAI,1),SIZE(I%M%T%XLAI,2)))
ALLOCATE(I%M%A%XALBUV_SOIL(SIZE(I%M%T%XLAI,1),SIZE(I%M%T%XLAI,2)))
 CALL SOIL_ALBEDO (I%O%CALBEDO, I%IP%XWSAT(:,1),I%R%XWG(:,1,:),     &
                    I%IP%XALBVIS_DRY,I%IP%XALBNIR_DRY,I%IP%XALBUV_DRY, &
                    I%IP%XALBVIS_WET,I%IP%XALBNIR_WET,I%IP%XALBUV_WET, &
                    I%M%A%XALBVIS_SOIL,I%M%A%XALBNIR_SOIL,I%M%A%XALBUV_SOIL )
!
ALLOCATE(I%R%XPSN   (SIZE(I%M%T%XLAI,1),SIZE(I%M%T%XLAI,2)))
ALLOCATE(I%R%XPSNG  (SIZE(I%M%T%XLAI,1),SIZE(I%M%T%XLAI,2)))
ALLOCATE(I%R%XPSNV  (SIZE(I%M%T%XLAI,1),SIZE(I%M%T%XLAI,2)))
ALLOCATE(I%R%XPSNV_A(SIZE(I%M%T%XLAI,1),SIZE(I%M%T%XLAI,2)))
I%R%XPSN    = 0.0
I%R%XPSNG   = 0.0
I%R%XPSNV   = 0.0
I%R%XPSNV_A = 0.0
ALLOCATE(I%I%XDIR_ALB_WITH_SNOW(SIZE(I%M%T%XLAI,1),1,SIZE(I%M%T%XLAI,2)))
ALLOCATE(I%I%XSCA_ALB_WITH_SNOW(SIZE(I%M%T%XLAI,1),1,SIZE(I%M%T%XLAI,2)))
I%I%XDIR_ALB_WITH_SNOW = 0.0
I%I%XSCA_ALB_WITH_SNOW = 0.0
CALL AVERAGED_ALBEDO_EMIS_ISBA(I, &
                               .FALSE., I%O%CALBEDO, ZZENITH,                &
                                 I%M%T%XVEG,I%M%T%XZ0,I%M%T%XLAI,                          &
                                 I%O%LMEB_PATCH,I%M%M%XGNDLITTER,I%M%M%XZ0LITTER,I%M%M%XLAIGV, &
                                 I%M%M%XH_VEG, I%R%XTV,               &
                                 I%R%XTG(:,1,:),I%IP%XPATCH, ZSW_BANDS,           &
                                 I%M%T%XALBNIR_VEG,I%M%T%XALBVIS_VEG,I%M%T%XALBUV_VEG,     &
                                 I%M%A%XALBNIR_SOIL,I%M%A%XALBVIS_SOIL,I%M%A%XALBUV_SOIL,  &
                                 I%M%T%XEMIS,                                  &
                                 I%R%TSNOW,                                  &
                                 I%M%T%XALBNIR,I%M%T%XALBVIS,I%M%T%XALBUV,                 &
                                 ZDIR_ALB, ZSCA_ALB,                     &
                                 ZEMIS,I%R%XTSRAD_NAT,ZTSURF                 )
DEALLOCATE(I%R%XPSN)
DEALLOCATE(I%R%XPSNG)
DEALLOCATE(I%R%XPSNV)
DEALLOCATE(I%R%XPSNV_A)
DEALLOCATE(I%I%XDIR_ALB_WITH_SNOW)
DEALLOCATE(I%I%XSCA_ALB_WITH_SNOW)
!
!-------------------------------------------------------------------------------------
!
!*      7.     Isba-Ags prognostic fields
!
IF (I%O%CPHOTO /= 'NON') THEN
!
   ALLOCATE(I%R%XAN(SIZE(I%M%T%XLAI,1),SIZE(I%M%T%XLAI,2)))
   I%R%XAN = 0.
!
   ALLOCATE(I%R%XANDAY(SIZE(I%M%T%XLAI,1),SIZE(I%M%T%XLAI,2)))
   I%R%XANDAY = 0.
!
   ALLOCATE(I%R%XANFM(SIZE(I%M%T%XLAI,1),SIZE(I%M%T%XLAI,2)))
   I%R%XANFM = XANFMINIT
!
   ALLOCATE(I%R%XLE(SIZE(I%M%T%XLAI,1),SIZE(I%M%T%XLAI,2)))
   I%R%XLE = 0.
!
   ALLOCATE(I%R%XRESP_BIOMASS(SIZE(I%M%T%XLAI,1),I%O%NNBIOMASS,SIZE(I%M%T%XLAI,2)))
   I%R%XRESP_BIOMASS(:,:,:) = 0.
!
ENDIF
!
IF (I%O%CPHOTO == 'AGS' .OR. I%O%CPHOTO == 'AST') THEN
!
   ALLOCATE(I%R%XBIOMASS(SIZE(I%M%T%XLAI,1),I%O%NNBIOMASS,SIZE(I%M%T%XLAI,2)))
   I%R%XBIOMASS(:,:,:) = 0.
!
ELSEIF (I%O%CPHOTO == 'LAI' .OR. I%O%CPHOTO == 'LST') THEN
!
   ALLOCATE(I%R%XBIOMASS(SIZE(I%M%T%XLAI,1),I%O%NNBIOMASS,SIZE(I%M%T%XLAI,2)))
   WHERE(I%M%T%XLAI(:,:)/=XUNDEF)
         I%R%XBIOMASS(:,1,:) = I%M%T%XLAI(:,:) * I%M%T%XBSLAI(:,:)
   ELSEWHERE
         I%R%XBIOMASS(:,1,:) = 0.0
   ENDWHERE
!
ELSEIF (I%O%CPHOTO == 'NIT' .OR. I%O%CPHOTO == 'NCB') THEN
!
   CALL PREP_HOR_ISBA_CC_FIELD(DTCO, U, &
                               IG, I, &
                               HPROGRAM,'BIOMASS ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)   
!
ENDIF
!
!-------------------------------------------------------------------------------------
!
!*      8.     Isba-CC prognostic fields
!
IF (I%O%CRESPSL == 'CNT') THEN
!
!*      8.1    Litter
!
 CALL PREP_HOR_ISBA_CC_FIELD(DTCO, U, &
                               IG, I, &
                               HPROGRAM,'LITTER  ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!*      8.2    Soil carbon
!
 CALL PREP_HOR_ISBA_CC_FIELD(DTCO, U, &
                               IG, I, &
                               HPROGRAM,'SOILCARB',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!*      8.2    lignin
!
 CALL PREP_HOR_ISBA_CC_FIELD(DTCO, U, &
                               IG, I, &
                               HPROGRAM,'LIGNIN  ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
ENDIF
!
!-------------------------------------------------------------------------------------
 CALL CLEAN_PREP_OUTPUT_GRID
!-------------------------------------------------------------------------------------
!
!*      10.     Preparation of canopy air variables
!
!
I%O%LCANOPY = LISBA_CANOPY
IF (I%O%LCANOPY) CALL PREP_ISBA_CANOPY(ICP, IG)
!
IF (LHOOK) CALL DR_HOOK('PREP_ISBA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE PREP_ISBA
