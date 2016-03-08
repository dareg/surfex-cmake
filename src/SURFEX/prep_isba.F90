!     #########
SUBROUTINE PREP_ISBA (DTCO, UG, U, USS, ICP, IG, IO, PZS, IP, I, IR, MT, MX, MM, MA, &
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
!
!
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SSO_n, ONLY : SSO_t
!
USE MODD_CANOPY_n, ONLY : CANOPY_t
USE MODD_GRID_n, ONLY : GRID_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_t, ISBA_INIT_PGD_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_TIME_t, ISBA_PARAM_FIX_t, ISBA_PARAM_ALB_t, &
                              ISBA_PARAM_MEB_t
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
!
USE MODI_PREP_HOR_ISBA_FIELD
USE MODI_PREP_VER_ISBA
USE MODI_PREP_OUTPUT_GRID
USE MODI_GET_LUOUT
USE MODI_PREP_ISBA_CANOPY
!
USE MODD_READ_NAMELIST,  ONLY : LNAM_READ
USE MODD_SURF_ATM,       ONLY : LVERTSHIFT
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
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SSO_t), INTENT(INOUT) :: USS
!
TYPE(CANOPY_t), INTENT(INOUT) :: ICP
TYPE(GRID_t), INTENT(INOUT) :: IG
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
REAL, DIMENSION(:), INTENT(IN) :: PZS
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
TYPE(ISBA_INIT_t), INTENT(INOUT) :: I
TYPE(ISBA_PROG_t), INTENT(INOUT) :: IR
TYPE(ISBA_PARAM_TIME_t), INTENT(INOUT) :: MT
TYPE(ISBA_PARAM_FIX_t), INTENT(INOUT) :: MX
TYPE(ISBA_PARAM_MEB_t), INTENT(INOUT) :: MM
TYPE(ISBA_PARAM_ALB_t), INTENT(INOUT) :: MA
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
REAL    :: ZWORK, ZLOG, ZWTOT, ZMATPOT, ZWL
!
REAL,             DIMENSION(1)   :: ZSW_BANDS ! middle wavelength of each band
REAL,             DIMENSION(SIZE(MT%XLAI,1),SIZE(MT%XLAI,2)) :: ZDIR_ALB  ! direct albedo for each band
REAL,             DIMENSION(SIZE(MT%XLAI,1),SIZE(MT%XLAI,2)) :: ZSCA_ALB  ! diffuse albedo for each band
REAL,             DIMENSION(SIZE(MT%XLAI,1))   :: ZEMIS     ! emissivity
REAL,             DIMENSION(SIZE(MT%XLAI,1))   :: ZZENITH   ! solar zenithal angle
REAL,             DIMENSION(SIZE(MT%XLAI,1))   :: ZTSURF     ! surface effective temperature
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
ISIZE_LMEB_PATCH=COUNT(IO%LMEB_PATCH(:))
!
!*      1.1    Default
!
 CALL GET_LUOUT(HPROGRAM,ILUOUT)
!
 CALL PREP_OUTPUT_GRID(UG%G, IG, U%NSIZE_FULL, ILUOUT)
!
!-------------------------------------------------------------------------------------
!
!*      2.     Reading and horizontal interpolations
!
!
!*      2.0    Large scale orography
!
 CALL PREP_HOR_ISBA_FIELD(DTCO, UG, U, USS, IG, IP, IO, IR, I%TTIME, MT%XLAI, MX, &
                          HPROGRAM,'ZS     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!*      2.1    Soil Water reservoirs
!
 CALL PREP_HOR_ISBA_FIELD(DTCO, UG, U, USS, IG, IP, IO, IR, I%TTIME, MT%XLAI, MX, &
                          HPROGRAM,'WG     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,GWG)
!
!*      2.2    Soil ice reservoirs
!
 CALL PREP_HOR_ISBA_FIELD(DTCO, UG, U, USS, IG, IP, IO, IR, I%TTIME, MT%XLAI, MX, &
                          HPROGRAM,'WGI    ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,GWGI)
!
!*      2.3    Leaves interception water reservoir
!
 CALL PREP_HOR_ISBA_FIELD(DTCO, UG, U, USS, IG, IP, IO, IR, I%TTIME, MT%XLAI, MX, &
                          HPROGRAM,'WR     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!*      2.4    Temperature profile
!
 CALL PREP_HOR_ISBA_FIELD(DTCO, UG, U, USS, IG, IP, IO, IR, I%TTIME, MT%XLAI, MX, &
                          HPROGRAM,'TG     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,GTG)
!
!*      2.5    Snow variables
!
 CALL PREP_HOR_ISBA_FIELD(DTCO, UG, U, USS, IG, IP, IO, IR, I%TTIME, MT%XLAI, MX, &
                          HPROGRAM,'SN_VEG ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,GPERMSNOW)
!
!*      2.6    LAI
!
 CALL PREP_HOR_ISBA_FIELD(DTCO, UG, U, USS, IG, IP, IO, IR, I%TTIME, MT%XLAI, MX, &
                          HPROGRAM,'LAI    ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!*      2.7    GLACIER
!
IF(IO%LGLACIER)THEN
  CALL PREP_HOR_ISBA_FIELD(DTCO, UG, U, USS, IG, IP, IO, IR, I%TTIME, MT%XLAI, MX, &
                          HPROGRAM,'ICE_STO',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
ENDIF
!
!*      2.8    Canopy vegetation temperature and interception reservoirs and air variables
!
IF(ISIZE_LMEB_PATCH>0)THEN
  CALL PREP_HOR_ISBA_FIELD(DTCO, UG, U, USS, IG, IP, IO, IR, I%TTIME, MT%XLAI, MX, &
                          HPROGRAM,'TV     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  CALL PREP_HOR_ISBA_FIELD(DTCO, UG, U, USS, IG, IP, IO, IR, I%TTIME, MT%XLAI, MX, &
                        HPROGRAM,'TL     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  CALL PREP_HOR_ISBA_FIELD(DTCO, UG, U, USS, IG, IP, IO, IR, I%TTIME, MT%XLAI, MX, &
                          HPROGRAM,'WRL    ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  CALL PREP_HOR_ISBA_FIELD(DTCO, UG, U, USS, IG, IP, IO, IR, I%TTIME, MT%XLAI, MX, &
                          HPROGRAM,'WRLI   ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  CALL PREP_HOR_ISBA_FIELD(DTCO, UG, U, USS, IG, IP, IO, IR, I%TTIME, MT%XLAI, MX, &
                          HPROGRAM,'WRVN   ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  CALL PREP_HOR_ISBA_FIELD(DTCO, UG, U, USS, IG, IP, IO, IR, I%TTIME, MT%XLAI, MX, &
                          HPROGRAM,'TC     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
  CALL PREP_HOR_ISBA_FIELD(DTCO, UG, U, USS, IG, IP, IO, IR, I%TTIME, MT%XLAI, MX, &
                          HPROGRAM,'QC     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
ENDIF
!
!-------------------------------------------------------------------------------------
!
!*      3.    Physical limitation: 
!
! No ice for force restore third layer:
IF (IO%CISBA == '3-L') THEN
   DO JP=1,IO%NPATCH
      WHERE(IR%XWG(:,3,JP) /= XUNDEF)
        IR%XWG(:,3,JP)  = MIN(IR%XWG(:,3,JP)+IR%XWGI(:,3,JP),IP%XWSAT(:,3))
        IR%XWGI(:,3,JP) = 0.
      END WHERE
   ENDDO
ENDIF
!
! Total water content should not exceed saturation:
DO JP=1,IO%NPATCH
   WHERE(IR%XWG(:,:,JP) /= XUNDEF .AND. (IR%XWG(:,:,JP) + IR%XWGI(:,:,JP)) > IP%XWSAT(:,:) )
      IR%XWGI(:,:,JP) = IP%XWSAT(:,:) - IR%XWG(:,:,JP)
   END WHERE
ENDDO
!
!-------------------------------------------------------------------------------------
!
!*      3.     Vertical interpolations of all variables
!
IF(LVERTSHIFT)THEN
  CALL PREP_VER_ISBA(IO, IR, PZS, MX%XDG)
ENDIF
!
DEALLOCATE(XZS_LS)
!-------------------------------------------------------------------------------------
!
!*      4.     Treatment of permanent snow
!
IF (GPERMSNOW.AND.LSNOW_PREP_PERM) CALL PREP_PERM_SNOW(IO, IP, IR)
!
CALL INIT_SNOW_LW(XEMISSN,IR%TSNOW)
!
IF (LPHYSDOMC) THEN
   IR%TSNOW%WSNOW(:,:,:)=0.
ENDIF 
!------------------------------------------------------------------------------------- 
! 
!*      4.b     Possibility for setting an upper limit on the initial snow water equivalent field 
IF (LSWEMAX) THEN 
  SMAX = MAXVAL(IR%TSNOW%WSNOW(:,:,:)) 
  WRITE(*,*) ' MAX(Snow content (kg/m2)): ', SMAX 
  WRITE(*,*) ' Set MAX to', XSWEMAX, '(kg/m2)' 
  IR%TSNOW%WSNOW(:,:,:) = MIN(IR%TSNOW%WSNOW(:,:,:),XSWEMAX) 
  SMAX = MAXVAL(IR%TSNOW%WSNOW(:,:,:)) 
  WRITE(*,*) ' MAX(Snow content (kg/m2)): ', SMAX 
ENDIF
!
!-------------------------------------------------------------------------------------
!
!*      5.     coherence between soil temperature and liquid/solid water
!
GTEMP2WGI=(GWG.OR.GWGI.OR.GTG)
!
IF (IO%CISBA == 'DIF'.AND.GTEMP2WGI) THEN
   INI=SIZE(IP%XWSAT,1)
   DO JP=1,IO%NPATCH
      DO JL=1,IO%NGROUND_LAYER
         DO JJ=1,INI
            IF(IR%XWG(JJ,JL,JP)/=XUNDEF)THEN
!     
!             total soil moisture
              ZWTOT = IR%XWG(JJ,JL,JP)+IR%XWGI(JJ,JL,JP)
              ZWTOT = MIN(ZWTOT,IP%XWSAT(JJ,JL))
!              
!             total matric potential
!             psi=mpotsat*(w/wsat)**(-bcoef)
              ZWORK   = ZWTOT/IP%XWSAT(JJ,JL)
              ZLOG    = IP%XBCOEF(JJ,JL)*LOG(ZWORK)
              ZMATPOT = IP%XMPOTSAT(JJ,JL)*EXP(-ZLOG)
!
!             soil liquid water content computation
!             w=wsat*(psi/mpotsat)**(-1/bcoef)
              ZMATPOT       = MIN(IP%XMPOTSAT(JJ,JL),XLMTT*(IR%XTG(JJ,JL,JP)-XTT)/(XG*IR%XTG(JJ,JL,JP)))
              ZWORK         = MAX(1.0,ZMATPOT/IP%XMPOTSAT(JJ,JL))
              ZLOG          = LOG(ZWORK)
              ZWL           = IP%XWSAT(JJ,JL)*EXP(-ZLOG/IP%XBCOEF(JJ,JL))
              ZWL           = MAX(ZWL,XWGMIN)
              IR%XWG(JJ,JL,JP) = MIN(ZWL,ZWTOT )
!        
!             soil ice computation    
              IR%XWGI(JJ,JL,JP) = MAX(0.0,ZWTOT-IR%XWG(JJ,JL,JP))
! 
!             supress numerical artefact
              IF(IR%XTG(JJ,JL,JP)>=XTT)THEN
                IR%XWG (JJ,JL,JP) = MIN(IR%XWG(JJ,JL,JP)+IR%XWGI(JJ,JL,JP),IP%XWSAT(JJ,JL))
                IR%XWGI(JJ,JL,JP) = 0.0
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
ALLOCATE(IR%XRESA(SIZE(MT%XLAI,1),SIZE(MT%XLAI,2)))
IR%XRESA = 100.
!
ALLOCATE(IR%XTSRAD_NAT(SIZE(MT%XLAI,1)))
ZZENITH(:)=0.
ZSW_BANDS(:)=0.
!
ALLOCATE(MT%XALBNIR(SIZE(MT%XLAI,1),SIZE(MT%XLAI,2)))
ALLOCATE(MT%XALBVIS(SIZE(MT%XLAI,1),SIZE(MT%XLAI,2)))
ALLOCATE(MT%XALBUV(SIZE(MT%XLAI,1),SIZE(MT%XLAI,2)))
MT%XALBNIR = 0.0
MT%XALBVIS = 0.0
MT%XALBUV = 0.0
!
ALLOCATE(MA%XALBNIR_SOIL(SIZE(MT%XLAI,1),SIZE(MT%XLAI,2)))
ALLOCATE(MA%XALBVIS_SOIL(SIZE(MT%XLAI,1),SIZE(MT%XLAI,2)))
ALLOCATE(MA%XALBUV_SOIL(SIZE(MT%XLAI,1),SIZE(MT%XLAI,2)))
 CALL SOIL_ALBEDO (IO%CALBEDO, IP%XWSAT(:,1),IR%XWG(:,1,:), IP, MA, "ALL" )
!
ALLOCATE(IR%XPSN   (SIZE(MT%XLAI,1),SIZE(MT%XLAI,2)))
ALLOCATE(IR%XPSNG  (SIZE(MT%XLAI,1),SIZE(MT%XLAI,2)))
ALLOCATE(IR%XPSNV  (SIZE(MT%XLAI,1),SIZE(MT%XLAI,2)))
ALLOCATE(IR%XPSNV_A(SIZE(MT%XLAI,1),SIZE(MT%XLAI,2)))
IR%XPSN    = 0.0
IR%XPSNG   = 0.0
IR%XPSNV   = 0.0
IR%XPSNV_A = 0.0
ALLOCATE(I%XDIR_ALB_WITH_SNOW(SIZE(MT%XLAI,1),1,SIZE(MT%XLAI,2)))
ALLOCATE(I%XSCA_ALB_WITH_SNOW(SIZE(MT%XLAI,1),1,SIZE(MT%XLAI,2)))
I%XDIR_ALB_WITH_SNOW = 0.0
I%XSCA_ALB_WITH_SNOW = 0.0
CALL AVERAGED_ALBEDO_EMIS_ISBA(IO, MT, MM, MA, IP, I, IR, MX%XVEGTYPE, &
                               ZZENITH, IR%XTG(:,1,:), ZSW_BANDS, ZDIR_ALB, ZSCA_ALB,   &
                               ZEMIS,IR%XTSRAD_NAT,ZTSURF                 )
DEALLOCATE(IR%XPSN)
DEALLOCATE(IR%XPSNG)
DEALLOCATE(IR%XPSNV)
DEALLOCATE(IR%XPSNV_A)
DEALLOCATE(I%XDIR_ALB_WITH_SNOW)
DEALLOCATE(I%XSCA_ALB_WITH_SNOW)
!
!-------------------------------------------------------------------------------------
!
!*      7.     Isba-Ags prognostic fields
!
IF (IO%CPHOTO /= 'NON') THEN
!
   ALLOCATE(IR%XAN(SIZE(MT%XLAI,1),SIZE(MT%XLAI,2)))
   IR%XAN = 0.
!
   ALLOCATE(IR%XANDAY(SIZE(MT%XLAI,1),SIZE(MT%XLAI,2)))
   IR%XANDAY = 0.
!
   ALLOCATE(IR%XANFM(SIZE(MT%XLAI,1),SIZE(MT%XLAI,2)))
   IR%XANFM = XANFMINIT
!
   ALLOCATE(IR%XLE(SIZE(MT%XLAI,1),SIZE(MT%XLAI,2)))
   IR%XLE = 0.
!
   ALLOCATE(IR%XRESP_BIOMASS(SIZE(MT%XLAI,1),IO%NNBIOMASS,SIZE(MT%XLAI,2)))
   IR%XRESP_BIOMASS(:,:,:) = 0.
!
ENDIF
!
IF (IO%CPHOTO == 'AGS' .OR. IO%CPHOTO == 'AST') THEN
!
   ALLOCATE(IR%XBIOMASS(SIZE(MT%XLAI,1),IO%NNBIOMASS,SIZE(MT%XLAI,2)))
   IR%XBIOMASS(:,:,:) = 0.
!
ELSEIF (IO%CPHOTO == 'LAI' .OR. IO%CPHOTO == 'LST') THEN
!
   ALLOCATE(IR%XBIOMASS(SIZE(MT%XLAI,1),IO%NNBIOMASS,SIZE(MT%XLAI,2)))
   WHERE(MT%XLAI(:,:)/=XUNDEF)
         IR%XBIOMASS(:,1,:) = MT%XLAI(:,:) * MT%XBSLAI(:,:)
   ELSEWHERE
         IR%XBIOMASS(:,1,:) = 0.0
   ENDWHERE
!
ELSEIF (IO%CPHOTO == 'NIT' .OR. IO%CPHOTO == 'NCB') THEN
!
   CALL PREP_HOR_ISBA_CC_FIELD(DTCO, U, SIZE(IG%XLAT), IP, IO, IR, MT%XLAI, MX%XVEGTYPE,  &
                               HPROGRAM,'BIOMASS ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)   
!
ENDIF
!
!-------------------------------------------------------------------------------------
!
!*      8.     Isba-CC prognostic fields
!
IF (IO%CRESPSL == 'CNT') THEN
!
!*      8.1    Litter
!
 CALL PREP_HOR_ISBA_CC_FIELD(DTCO, U, SIZE(IG%XLAT), IP, IO, IR, MT%XLAI, MX%XVEGTYPE,  &
                               HPROGRAM,'LITTER  ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!*      8.2    Soil carbon
!
 CALL PREP_HOR_ISBA_CC_FIELD(DTCO, U, SIZE(IG%XLAT), IP, IO, IR, MT%XLAI, MX%XVEGTYPE,  &
                               HPROGRAM,'SOILCARB',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!*      8.2    lignin
!
 CALL PREP_HOR_ISBA_CC_FIELD(DTCO, U, SIZE(IG%XLAT), IP, IO, IR, MT%XLAI, MX%XVEGTYPE,  &
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
IO%LCANOPY = LISBA_CANOPY
IF (IO%LCANOPY) CALL PREP_ISBA_CANOPY(ICP, IG%NDIM)
!
IF (LHOOK) CALL DR_HOOK('PREP_ISBA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE PREP_ISBA
