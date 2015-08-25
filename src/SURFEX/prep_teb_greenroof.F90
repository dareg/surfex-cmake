!     #########
SUBROUTINE PREP_TEB_GREENROOF (DTCO, IOB, IG, I, UG, U, USS, TGR, TGRO, TGRPE, TGRP, &
                               TG, T, TOP, TVG, &
                               HPROGRAM,HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
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
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
USE MODD_ISBA_GRID_n, ONLY : ISBA_GRID_t
USE MODD_ISBA_n, ONLY : ISBA_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SURF_ATM_SSO_n, ONLY : SURF_ATM_SSO_t
USE MODD_TEB_GREENROOF_n, ONLY : TEB_GREENROOF_t
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TEB_GREENROOF_OPTIONS_t
USE MODD_TEB_GREENROOF_PGD_EVOL_n, ONLY : TEB_GREENROOF_PGD_EVOL_t
USE MODD_TEB_GREENROOF_PGD_n, ONLY : TEB_GREENROOF_PGD_t
USE MODD_TEB_GRID_n, ONLY : TEB_GRID_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_OPTIONS_t
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
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
TYPE(ISBA_GRID_t), INTENT(INOUT) :: IG
TYPE(ISBA_t), INTENT(INOUT) :: I
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SURF_ATM_SSO_t), INTENT(INOUT) :: USS
TYPE(TEB_GREENROOF_t), INTENT(INOUT) :: TGR
TYPE(TEB_GREENROOF_OPTIONS_t), INTENT(INOUT) :: TGRO
TYPE(TEB_GREENROOF_PGD_EVOL_t), INTENT(INOUT) :: TGRPE
TYPE(TEB_GREENROOF_PGD_t), INTENT(INOUT) :: TGRP
TYPE(TEB_GRID_t), INTENT(INOUT) :: TG
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_VEG_OPTIONS_t), INTENT(INOUT) :: TVG
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM    ! program calling surf. schemes
 CHARACTER(LEN=28),  INTENT(IN)  :: HATMFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),   INTENT(IN)  :: HATMFILETYPE! type of the Atmospheric file
 CHARACTER(LEN=28),  INTENT(IN)  :: HPGDFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),   INTENT(IN)  :: HPGDFILETYPE! type of the Atmospheric file
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
 CALL PREP_HOR_TEB_GREENROOF_FIELD(DTCO, IOB, IG, I, UG, U, USS, TGR, TGRO, TGRPE, TGRP, &
                                         TG, TOP, &
                                   HPROGRAM,'WG     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!*      2.2    Soil ice reservoirs
!
 CALL PREP_HOR_TEB_GREENROOF_FIELD(DTCO, IOB, IG, I, UG, U, USS, TGR, TGRO, TGRPE, TGRP, &
                                         TG, TOP, &
                                   HPROGRAM,'WGI    ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!*      2.3    Leaves interception water reservoir
!
 CALL PREP_HOR_TEB_GREENROOF_FIELD(DTCO, IOB, IG, I, UG, U, USS, TGR, TGRO, TGRPE, TGRP, &
                                         TG, TOP, &
                                   HPROGRAM,'WR     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!*      2.4    Temperature profile
!
 CALL PREP_HOR_TEB_GREENROOF_FIELD(DTCO, IOB, IG, I, UG, U, USS, TGR, TGRO, TGRPE, TGRP, &
                                         TG, TOP, &
                                   HPROGRAM,'TG     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
! Initializing deep GR temp. with that of the outer layer of the structural roof 
!
TGRP%XTDEEP(:) = T%CUR%XT_ROOF(:,1)
!
!*      2.5    Snow variables
!
 CALL PREP_HOR_TEB_GREENROOF_FIELD(DTCO, IOB, IG, I, UG, U, USS, TGR, TGRO, TGRPE, TGRP, &
                                         TG, TOP, &
                                   HPROGRAM,'SN_VEG ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!*      2.6    LAI
!
 CALL PREP_HOR_TEB_GREENROOF_FIELD(DTCO, IOB, IG, I, UG, U, USS, TGR, TGRO, TGRPE, TGRP, &
                                         TG, TOP, &
                                   HPROGRAM,'LAI    ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!-------------------------------------------------------------------------------------
!
!*      3.    Physical limitations: 
!
! 3.1  If whole ice reservoir is empty (grib from ecmwf case) and surface temperature is
!      lower than -10C, then ice content is maximum and water content minimum
!
IF (ALL(TGR%CUR%XWGI(:,:)==0.)) THEN
   WHERE(TGR%CUR%XTG(:,1:SIZE(TGR%CUR%XWG,2)) < XTT-10.)
      TGR%CUR%XWGI(:,:) = TGRP%XWSAT(:,:)-XWGMIN
      TGR%CUR%XWG (:,:) = XWGMIN
   END WHERE
ENDIF
!
!
! 3.2.  Total water content should not exceed saturation:
WHERE(TGR%CUR%XWG(:,:) /= XUNDEF .AND. (TGR%CUR%XWG(:,:) + TGR%CUR%XWGI(:,:)) > TGRP%XWSAT(:,:) )
   TGR%CUR%XWGI(:,:) = TGRP%XWSAT(:,:) - TGR%CUR%XWG(:,:)
END WHERE
!
!-------------------------------------------------------------------------------------
!
!*      4.     Vertical interpolations of all variables
!
IF(LVERTSHIFT)THEN
  CALL PREP_VER_TEB_GREENROOF(TGR, TGRO, TGRP, TOP)
ENDIF
!
!
!-------------------------------------------------------------------------------------
!
!*      5.     Half prognostic fields
!
ALLOCATE(TGR%CUR%XRESA(SIZE(TGRPE%CUR%XLAI)))
TGR%CUR%XRESA(:) = 100.
!
!-------------------------------------------------------------------------------------
!
!*      6.     Isba-Ags prognostic fields
!
IF (TVG%CPHOTO /= 'NON') THEN
!
   ALLOCATE(TGR%CUR%XAN(SIZE(TGRPE%CUR%XLAI)))
   TGR%CUR%XAN = 0.
!
   ALLOCATE(TGR%CUR%XANDAY(SIZE(TGRPE%CUR%XLAI)))
   TGR%CUR%XANDAY = 0.
!
   ALLOCATE(TGR%CUR%XANFM(SIZE(TGRPE%CUR%XLAI)))
   TGR%CUR%XANFM = XANFMINIT
!
   ALLOCATE(TGR%CUR%XLE(SIZE(TGRPE%CUR%XLAI)))
   TGR%CUR%XLE = 0.
!
ENDIF
!
IF (TVG%CPHOTO == 'AGS' .OR. TVG%CPHOTO == 'AST') THEN
!
   ALLOCATE(TGR%CUR%XBIOMASS(SIZE(TGRPE%CUR%XLAI),TVG%NNBIOMASS))
   TGR%CUR%XBIOMASS(:,1) = 0.
!
   ALLOCATE(TGR%CUR%XRESP_BIOMASS(SIZE(TGRPE%CUR%XLAI),TVG%NNBIOMASS))
   TGR%CUR%XRESP_BIOMASS(:,:) = 0.
!
ELSEIF (TVG%CPHOTO == 'LAI' .OR. TVG%CPHOTO == 'LST') THEN
!
   ALLOCATE(TGR%CUR%XBIOMASS(SIZE(TGRPE%CUR%XLAI),TVG%NNBIOMASS))
   TGR%CUR%XBIOMASS(:,1) = TGRPE%CUR%XLAI(:) * TGRP%XBSLAI(:)
!
   ALLOCATE(TGR%CUR%XRESP_BIOMASS(SIZE(TGRPE%CUR%XLAI),TVG%NNBIOMASS))
   TGR%CUR%XRESP_BIOMASS(:,:) = 0.
!
ELSEIF (TVG%CPHOTO == 'NIT' .OR. TVG%CPHOTO == 'NCB') THEN
!
   ALLOCATE(TGR%CUR%XBIOMASS(SIZE(TGRPE%CUR%XLAI),TVG%NNBIOMASS))
   TGR%CUR%XBIOMASS(:,1) = TGRPE%CUR%XLAI(:) * TGRP%XBSLAI_NITRO(:)
   TGR%CUR%XBIOMASS(:,2) = MAX( 0., (TGR%CUR%XBIOMASS(:,1)/ (XCC_NIT/10.**XCA_NIT))  &
                              **(1.0/(1.0-XCA_NIT)) - TGR%CUR%XBIOMASS(:,1) )  
   TGR%CUR%XBIOMASS(:,3:TVG%NNBIOMASS) = 0.
!
   ALLOCATE(TGR%CUR%XRESP_BIOMASS(SIZE(TGRPE%CUR%XLAI),TVG%NNBIOMASS))
   TGR%CUR%XRESP_BIOMASS(:,:) = 0.
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
