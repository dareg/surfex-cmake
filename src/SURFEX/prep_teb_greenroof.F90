!     #########
SUBROUTINE PREP_TEB_GREENROOF (DTCO, UG, U, USS, IG, TG, T, TOP, GRM, &
                               HPROGRAM,HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)
!     #################################################################################
!
!!****  *PREP_TEB_GREENROOF* - Prepares ISBA fields for greenroofs
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!    Based on "prep_teb_garden"
!!
!!    REFERENCE
!!    ---------
!!      
!!
!!    AUTHOR
!!    ------
!!      A. Lemonsu & C. de Munck 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/2011
!!------------------------------------------------------------------
!
!
!
!
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SURF_ATM_SSO_n, ONLY : SURF_ATM_SSO_t
USE MODD_ISBA_GRID_n, ONLY : ISBA_GRID_t
USE MODD_TEB_GRID_n, ONLY : TEB_GRID_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_SURFEX_n, ONLY : TEB_GREENROOF_MODEL_t
!
USE MODI_PREP_HOR_TEB_GREENROOF_FIELD
USE MODI_PREP_VER_TEB_GREENROOF
!
                                ! A FAIRE :
                                ! IL FAUT RAJOUTER TSNOW
                                ! ----------------------
USE MODD_SURF_ATM,       ONLY : LVERTSHIFT
USE MODD_CSTS,           ONLY : XTT
USE MODD_SNOW_PAR,       ONLY : XZ0SN
USE MODD_ISBA_PAR,       ONLY : XWGMIN
USE MODD_CO2V_PAR,       ONLY : XCC_NIT, XCA_NIT, XANFMINIT
USE MODD_SURF_PAR,       ONLY : XUNDEF
!
USE MODN_PREP_ISBA
USE MODE_POS_SURF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SURF_ATM_SSO_t), INTENT(INOUT) :: USS
TYPE(ISBA_GRID_t), INTENT(INOUT) :: IG
TYPE(TEB_GRID_t), INTENT(INOUT) :: TG
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_GREENROOF_MODEL_t), INTENT(INOUT) :: GRM
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM    ! program calling surf. schemes
 CHARACTER(LEN=28),  INTENT(IN)  :: HATMFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),   INTENT(IN)  :: HATMFILETYPE! type of the Atmospheric file
 CHARACTER(LEN=28),  INTENT(IN)  :: HPGDFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),   INTENT(IN)  :: HPGDFILETYPE! type of the Atmospheric file
!
INTEGER,            INTENT(IN)  :: KPATCH
!
!*      0.2    declarations of local variables
!
INTEGER :: JP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
!*      1.     Default of configuration
!
!*      1.1    Default
!
!
!-------------------------------------------------------------------------------------
!
!*      2.     Reading and horizontal interpolations
!
!
!*      2.1    Soil Water reservoirs
!
IF (LHOOK) CALL DR_HOOK('PREP_TEB_GREENROOF',0,ZHOOK_HANDLE)
!
 CALL PREP_HOR_TEB_GREENROOF_FIELD(DTCO, IG, UG, U, USS, GRM%TV%R, GRM%TV%O, GRM%TV%M, GRM%TV%IP, &
                                         TG, TOP, &
                                   HPROGRAM,'WG     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)
!
!*      2.2    Soil ice reservoirs
!
 CALL PREP_HOR_TEB_GREENROOF_FIELD(DTCO, IG, UG, U, USS, GRM%TV%R, GRM%TV%O, GRM%TV%M, GRM%TV%IP, &
                                         TG, TOP, &
                                   HPROGRAM,'WGI    ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)
!
!*      2.3    Leaves interception water reservoir
!
 CALL PREP_HOR_TEB_GREENROOF_FIELD(DTCO, IG, UG, U, USS, GRM%TV%R, GRM%TV%O, GRM%TV%M, GRM%TV%IP, &
                                         TG, TOP, &
                                   HPROGRAM,'WR     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)
!
!*      2.4    Temperature profile
!
 CALL PREP_HOR_TEB_GREENROOF_FIELD(DTCO, IG, UG, U, USS, GRM%TV%R, GRM%TV%O, GRM%TV%M, GRM%TV%IP, &
                                         TG, TOP, &
                                   HPROGRAM,'TG     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)
!
! Initializing deep GR temp. with that of the outer layer of the structural roof 
!
GRM%TV%IP%XTDEEP(:) = T%CUR%XT_ROOF(:,1)
!
!*      2.5    Snow variables
!
 CALL PREP_HOR_TEB_GREENROOF_FIELD(DTCO, IG, UG, U, USS, GRM%TV%R, GRM%TV%O, GRM%TV%M, GRM%TV%IP, &
                                         TG, TOP, &
                                   HPROGRAM,'SN_VEG ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)
!
!*      2.6    LAI
!
 CALL PREP_HOR_TEB_GREENROOF_FIELD(DTCO, IG, UG, U, USS, GRM%TV%R, GRM%TV%O, GRM%TV%M, GRM%TV%IP, &
                                         TG, TOP, &
                                   HPROGRAM,'LAI    ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)
!
!-------------------------------------------------------------------------------------
!
!*      3.    Physical limitations: 
!
! 3.1  If whole ice reservoir is empty (grib from ecmwf case) and surface temperature is
!      lower than -10C, then ice content is maximum and water content minimum
!
IF (ALL(GRM%TV%R%CUR%XWGI(:,:,1)==0.)) THEN
   WHERE(GRM%TV%R%CUR%XTG(:,1:SIZE(GRM%TV%R%CUR%XWG,2),1) < XTT-10.)
      GRM%TV%R%CUR%XWGI(:,:,1) = GRM%TV%IP%XWSAT(:,:)-XWGMIN
      GRM%TV%R%CUR%XWG (:,:,1) = XWGMIN
   END WHERE
ENDIF
!
!
! 3.2.  Total water content should not exceed saturation:
WHERE(GRM%TV%R%CUR%XWG(:,:,1) /= XUNDEF .AND. (GRM%TV%R%CUR%XWG(:,:,1) + GRM%TV%R%CUR%XWGI(:,:,1)) > GRM%TV%IP%XWSAT(:,:) )
   GRM%TV%R%CUR%XWGI(:,:,1) = GRM%TV%IP%XWSAT(:,:) - GRM%TV%R%CUR%XWG(:,:,1)
END WHERE
!
!-------------------------------------------------------------------------------------
!
!*      4.     Vertical interpolations of all variables
!
IF(LVERTSHIFT)THEN
  CALL PREP_VER_TEB_GREENROOF(GRM%TV%R, GRM%TV%O, GRM%TV%P, TOP, GRM%TV%M%X%XDG(:,:,1))
ENDIF
!
!
!-------------------------------------------------------------------------------------
!
!*      5.     Half prognostic fields
!
ALLOCATE(GRM%TV%R%CUR%XRESA(SIZE(GRM%TV%M%T%CUR%XLAI),1))
GRM%TV%R%CUR%XRESA(:,1) = 100.
!
!-------------------------------------------------------------------------------------
!
!*      6.     Isba-Ags prognostic fields
!
IF (GRM%TV%O%CPHOTO /= 'NON') THEN
!
   ALLOCATE(GRM%TV%R%CUR%XAN(SIZE(GRM%TV%M%T%CUR%XLAI),1))
   GRM%TV%R%CUR%XAN = 0.
!
   ALLOCATE(GRM%TV%R%CUR%XANDAY(SIZE(GRM%TV%M%T%CUR%XLAI),1))
   GRM%TV%R%CUR%XANDAY = 0.
!
   ALLOCATE(GRM%TV%R%CUR%XANFM(SIZE(GRM%TV%M%T%CUR%XLAI),1))
   GRM%TV%R%CUR%XANFM = XANFMINIT
!
   ALLOCATE(GRM%TV%R%CUR%XLE(SIZE(GRM%TV%M%T%CUR%XLAI),1))
   GRM%TV%R%CUR%XLE = 0.
!
ENDIF
!
IF (GRM%TV%O%CPHOTO == 'AGS' .OR. GRM%TV%O%CPHOTO == 'AST') THEN
!
   ALLOCATE(GRM%TV%R%CUR%XBIOMASS(SIZE(GRM%TV%M%T%CUR%XLAI),GRM%TV%O%NNBIOMASS,1))
   GRM%TV%R%CUR%XBIOMASS(:,1,1) = 0.
!
   ALLOCATE(GRM%TV%R%CUR%XRESP_BIOMASS(SIZE(GRM%TV%M%T%CUR%XLAI),GRM%TV%O%NNBIOMASS,1))
   GRM%TV%R%CUR%XRESP_BIOMASS(:,:,1) = 0.
!
ELSEIF (GRM%TV%O%CPHOTO == 'LAI' .OR. GRM%TV%O%CPHOTO == 'LST') THEN
!
   ALLOCATE(GRM%TV%R%CUR%XBIOMASS(SIZE(GRM%TV%M%T%CUR%XLAI),GRM%TV%O%NNBIOMASS,1))
   GRM%TV%R%CUR%XBIOMASS(:,1,1) = GRM%TV%M%T%CUR%XLAI(:,1) * GRM%TV%M%T%CUR%XBSLAI(:,1)
!
   ALLOCATE(GRM%TV%R%CUR%XRESP_BIOMASS(SIZE(GRM%TV%M%T%CUR%XLAI),GRM%TV%O%NNBIOMASS,1))
   GRM%TV%R%CUR%XRESP_BIOMASS(:,:,1) = 0.
!
ELSEIF (GRM%TV%O%CPHOTO == 'NIT' .OR. GRM%TV%O%CPHOTO == 'NCB') THEN
!
   ALLOCATE(GRM%TV%R%CUR%XBIOMASS(SIZE(GRM%TV%M%T%CUR%XLAI),GRM%TV%O%NNBIOMASS,1))
   GRM%TV%R%CUR%XBIOMASS(:,1,1) = GRM%TV%M%T%CUR%XLAI(:,1) * GRM%TV%IP%XBSLAI_NITRO(:,1)
   GRM%TV%R%CUR%XBIOMASS(:,2,1) = MAX( 0., (GRM%TV%R%CUR%XBIOMASS(:,1,1)/ (XCC_NIT/10.**XCA_NIT))  &
                              **(1.0/(1.0-XCA_NIT)) - GRM%TV%R%CUR%XBIOMASS(:,1,1) )  
   GRM%TV%R%CUR%XBIOMASS(:,3:GRM%TV%O%NNBIOMASS,1) = 0.
!
   ALLOCATE(GRM%TV%R%CUR%XRESP_BIOMASS(SIZE(GRM%TV%M%T%CUR%XLAI),GRM%TV%O%NNBIOMASS,1))
   GRM%TV%R%CUR%XRESP_BIOMASS(:,:,1) = 0.
!
ENDIF
!
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PREP_TEB_GREENROOF',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE PREP_TEB_GREENROOF
