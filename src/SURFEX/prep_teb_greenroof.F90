!     #########
SUBROUTINE PREP_TEB_GREENROOF (DTCO, UG, U, USS, TG, TOP, GRR, GRO, GRMX, GRMT, GRIP,  &
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
USE MODD_SSO_n, ONLY : SSO_t
!
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_FIX_t, ISBA_PARAM_TIME_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t
USE MODD_GRID_n, ONLY : GRID_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
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
TYPE(SSO_t), INTENT(INOUT) :: USS
!
!
TYPE(ISBA_PROG_t), INTENT(INOUT) :: GRR
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: GRO
TYPE(ISBA_PARAM_FIX_t), INTENT(INOUT) :: GRMX
TYPE(ISBA_PARAM_TIME_t), INTENT(INOUT) :: GRMT
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: GRIP
!
TYPE(GRID_t), INTENT(INOUT) :: TG
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
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
 CALL PREP_HOR_TEB_GREENROOF_FIELD(DTCO, UG, U, USS, GRR, GRO, GRMX, GRMT%XLAI, GRIP, TG, TOP, &
                                   HPROGRAM,'WG     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)
!
!*      2.2    Soil ice reservoirs
!
 CALL PREP_HOR_TEB_GREENROOF_FIELD(DTCO, UG, U, USS, GRR, GRO, GRMX, GRMT%XLAI, GRIP, TG, TOP, &
                                   HPROGRAM,'WGI    ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)
!
!*      2.3    Leaves interception water reservoir
!
 CALL PREP_HOR_TEB_GREENROOF_FIELD(DTCO, UG, U, USS, GRR, GRO, GRMX, GRMT%XLAI, GRIP, TG, TOP, &
                                   HPROGRAM,'WR     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)
!
!*      2.4    Temperature profile
!
 CALL PREP_HOR_TEB_GREENROOF_FIELD(DTCO, UG, U, USS, GRR, GRO, GRMX, GRMT%XLAI, GRIP, TG, TOP, &
                                   HPROGRAM,'TG     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)
!
!*      2.5    Snow variables
!
 CALL PREP_HOR_TEB_GREENROOF_FIELD(DTCO, UG, U, USS, GRR, GRO, GRMX, GRMT%XLAI, GRIP, TG, TOP, &
                                   HPROGRAM,'SN_VEG ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)
!
!*      2.6    LAI
!
 CALL PREP_HOR_TEB_GREENROOF_FIELD(DTCO, UG, U, USS, GRR, GRO, GRMX, GRMT%XLAI, GRIP, TG, TOP, &
                                   HPROGRAM,'LAI    ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)
!
!-------------------------------------------------------------------------------------
!
!*      3.    Physical limitations: 
!
! 3.1  If whole ice reservoir is empty (grib from ecmwf case) and surface temperature is
!      lower than -10C, then ice content is maximum and water content minimum
!
IF (ALL(GRR%XWGI(:,:,1)==0.)) THEN
   WHERE(GRR%XTG(:,1:SIZE(GRR%XWG,2),1) < XTT-10.)
      GRR%XWGI(:,:,1) = GRIP%XWSAT(:,:)-XWGMIN
      GRR%XWG (:,:,1) = XWGMIN
   END WHERE
ENDIF
!
!
! 3.2.  Total water content should not exceed saturation:
WHERE(GRR%XWG(:,:,1) /= XUNDEF .AND. (GRR%XWG(:,:,1) + GRR%XWGI(:,:,1)) > GRIP%XWSAT(:,:) )
   GRR%XWGI(:,:,1) = GRIP%XWSAT(:,:) - GRR%XWG(:,:,1)
END WHERE
!
!-------------------------------------------------------------------------------------
!
!*      4.     Vertical interpolations of all variables
!
IF(LVERTSHIFT)THEN
  CALL PREP_VER_TEB_GREENROOF(GRR, GRO, TOP%XZS, GRMX%XDG(:,:,1))
ENDIF
!
!
!-------------------------------------------------------------------------------------
!
!*      5.     Half prognostic fields
!
ALLOCATE(GRR%XRESA(SIZE(GRMT%XLAI),1))
GRR%XRESA(:,1) = 100.
!
!-------------------------------------------------------------------------------------
!
!*      6.     Isba-Ags prognostic fields
!
IF (GRO%CPHOTO /= 'NON') THEN
!
   ALLOCATE(GRR%XAN(SIZE(GRMT%XLAI),1))
   GRR%XAN = 0.
!
   ALLOCATE(GRR%XANDAY(SIZE(GRMT%XLAI),1))
   GRR%XANDAY = 0.
!
   ALLOCATE(GRR%XANFM(SIZE(GRMT%XLAI),1))
   GRR%XANFM = XANFMINIT
!
   ALLOCATE(GRR%XLE(SIZE(GRMT%XLAI),1))
   GRR%XLE = 0.
!
ENDIF
!
IF (GRO%CPHOTO == 'AGS' .OR. GRO%CPHOTO == 'AST') THEN
!
   ALLOCATE(GRR%XBIOMASS(SIZE(GRMT%XLAI),GRO%NNBIOMASS,1))
   GRR%XBIOMASS(:,1,1) = 0.
!
   ALLOCATE(GRR%XRESP_BIOMASS(SIZE(GRMT%XLAI),GRO%NNBIOMASS,1))
   GRR%XRESP_BIOMASS(:,:,1) = 0.
!
ELSEIF (GRO%CPHOTO == 'LAI' .OR. GRO%CPHOTO == 'LST') THEN
!
   ALLOCATE(GRR%XBIOMASS(SIZE(GRMT%XLAI),GRO%NNBIOMASS,1))
   GRR%XBIOMASS(:,1,1) = GRMT%XLAI(:,1) * GRMT%XBSLAI(:,1)
!
   ALLOCATE(GRR%XRESP_BIOMASS(SIZE(GRMT%XLAI),GRO%NNBIOMASS,1))
   GRR%XRESP_BIOMASS(:,:,1) = 0.
!
ELSEIF (GRO%CPHOTO == 'NIT' .OR. GRO%CPHOTO == 'NCB') THEN
!
   ALLOCATE(GRR%XBIOMASS(SIZE(GRMT%XLAI),GRO%NNBIOMASS,1))
   GRR%XBIOMASS(:,1,1) = GRMT%XLAI(:,1) * GRIP%XBSLAI_NITRO(:,1)
   GRR%XBIOMASS(:,2,1) = MAX( 0., (GRR%XBIOMASS(:,1,1)/ (XCC_NIT/10.**XCA_NIT))  &
                              **(1.0/(1.0-XCA_NIT)) - GRR%XBIOMASS(:,1,1) )  
   GRR%XBIOMASS(:,3:GRO%NNBIOMASS,1) = 0.
!
   ALLOCATE(GRR%XRESP_BIOMASS(SIZE(GRMT%XLAI),GRO%NNBIOMASS,1))
   GRR%XRESP_BIOMASS(:,:,1) = 0.
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
