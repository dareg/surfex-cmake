!     #########
SUBROUTINE PREP_TEB_GARDEN (TGD, TGDO, TGDPE, TGDP, TG, TOP, TVG, &
                            HPROGRAM,HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!     #################################################################################
!
!!****  *PREP_TEB_GARDEN* - Prepares ISBA fields
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
!!------------------------------------------------------------------
!
!
!
!
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
USE MODD_IO_BUFF_n, ONLY : IOB => IO_BUFF
USE MODD_ISBA_GRID_n, ONLY : IG => ISBA_GRID
USE MODD_ISBA_n, ONLY : I => ISBA
USE MODD_SURF_ATM_GRID_n, ONLY : UG => SURF_ATM_GRID
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
USE MODD_SURF_ATM_SSO_n, ONLY : USS => SURF_ATM_SSO
!
USE MODD_TEB_GARDEN_n, ONLY : TEB_GARDEN_t
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TEB_GARDEN_OPTIONS_t
USE MODD_TEB_GARDEN_PGD_EVOL_n, ONLY : TEB_GARDEN_PGD_EVOL_t
USE MODD_TEB_GARDEN_PGD_n, ONLY : TEB_GARDEN_PGD_t
USE MODD_TEB_GRID_n, ONLY : TEB_GRID_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_OPTIONS_t
!
USE MODI_PREP_HOR_TEB_GARDEN_FIELD
USE MODI_PREP_VER_TEB_GARDEN
!
USE MODD_SURF_ATM,       ONLY : LVERTSHIFT
!
                                ! A FAIRE :
                                ! IL FAUT RAJOUTER TSNOW
                                ! ----------------------
USE MODD_CSTS,        ONLY : XTT
USE MODD_SNOW_PAR,    ONLY : XZ0SN
USE MODD_ISBA_PAR,    ONLY : XWGMIN
USE MODD_CO2V_PAR,    ONLY : XANFMINIT, XCA_NIT, XCC_NIT
USE MODD_SURF_PAR,    ONLY : XUNDEF
!
USE MODN_PREP_ISBA
USE MODE_POS_SURF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(TEB_GARDEN_t), INTENT(INOUT) :: TGD
TYPE(TEB_GARDEN_OPTIONS_t), INTENT(INOUT) :: TGDO
TYPE(TEB_GARDEN_PGD_EVOL_t), INTENT(INOUT) :: TGDPE
TYPE(TEB_GARDEN_PGD_t), INTENT(INOUT) :: TGDP
TYPE(TEB_GRID_t), INTENT(INOUT) :: TG
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_VEG_OPTIONS_t), INTENT(INOUT) :: TVG
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=28),  INTENT(IN)  :: HATMFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),   INTENT(IN)  :: HATMFILETYPE! type of the Atmospheric file
 CHARACTER(LEN=28),  INTENT(IN)  :: HPGDFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),   INTENT(IN)  :: HPGDFILETYPE! type of the Atmospheric file
!
!*      0.2    declarations of local variables
!
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
IF (LHOOK) CALL DR_HOOK('PREP_TEB_GARDEN',0,ZHOOK_HANDLE)
 CALL PREP_HOR_TEB_GARDEN_FIELD(DTCO, IOB, IG, I, UG, U, USS, &
                                TGD, TGDO, TGDPE, TGDP, TG, TOP, TVG, &
                                HPROGRAM,'WG     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!*      2.2    Soil ice reservoirs
!
 CALL PREP_HOR_TEB_GARDEN_FIELD(DTCO, IOB, IG, I, UG, U, USS, &
                                TGD, TGDO, TGDPE, TGDP, TG, TOP, TVG, &
                                HPROGRAM,'WGI    ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!*      2.3    Leaves interception water reservoir
!
 CALL PREP_HOR_TEB_GARDEN_FIELD(DTCO, IOB, IG, I, UG, U, USS, &
                                TGD, TGDO, TGDPE, TGDP, TG, TOP, TVG, &
                                HPROGRAM,'WR     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!*      2.4    Temperature profile
!
 CALL PREP_HOR_TEB_GARDEN_FIELD(DTCO, IOB, IG, I, UG, U, USS, &
                                TGD, TGDO, TGDPE, TGDP, TG, TOP, TVG, &
                                HPROGRAM,'TG     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!*      2.5    Snow variables
!
 CALL PREP_HOR_TEB_GARDEN_FIELD(DTCO, IOB, IG, I, UG, U, USS, &
                                TGD, TGDO, TGDPE, TGDP, TG, TOP, TVG, &
                                HPROGRAM,'SN_VEG ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)

!
!*      2.6    LAI
!
IF (TVG%CPHOTO/='NON' .AND. TVG%CPHOTO/='AGS' .AND. TVG%CPHOTO/='LST')  &
 CALL PREP_HOR_TEB_GARDEN_FIELD(DTCO, IOB, IG, I, UG, U, USS, &
                                TGD, TGDO, TGDPE, TGDP, TG, TOP, TVG, &
                                HPROGRAM,'LAI    ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!-------------------------------------------------------------------------------------
!
!*      3.    Physical limitation: 
!
! If whole ice reservoir is empty (grib from ecmwf case) and surface temperature is
! lower than -10C, then ice content is maximum and water content minimum
!
IF (ALL(TGD%CUR%XWGI(:,:)==0.)) THEN
   WHERE(TGD%CUR%XTG(:,1:SIZE(TGD%CUR%XWG,2)) < XTT-10.)
       TGD%CUR%XWGI(:,:) = TGDP%XWSAT(:,:)-XWGMIN
       TGD%CUR%XWG (:,:) = XWGMIN
   END WHERE
ENDIF
!
! No ice for force restore third layer:
IF (TVG%CISBA == '3-L') THEN
      WHERE(TGD%CUR%XWG(:,3)/=XUNDEF.AND.TGD%CUR%XWGI(:,3)/=XUNDEF)
        TGD%CUR%XWG(:,3)  = MIN(TGD%CUR%XWG(:,3)+TGD%CUR%XWGI(:,3),TGDP%XWSAT(:,3))
        TGD%CUR%XWGI(:,3) = 0.
      END WHERE
ENDIF
!
! Total water content should not exceed saturation:
WHERE(TGD%CUR%XWG(:,:) /= XUNDEF .AND. (TGD%CUR%XWG(:,:) + TGD%CUR%XWGI(:,:)) > TGDP%XWSAT(:,:) )
     TGD%CUR%XWGI(:,:) = TGDP%XWSAT(:,:) - TGD%CUR%XWG(:,:)
END WHERE
!
!-------------------------------------------------------------------------------------
!
!*      3.     Vertical interpolations of all variables
!
IF(LVERTSHIFT)THEN
  CALL PREP_VER_TEB_GARDEN(TGD, TGDO, TGDP, TOP, TVG)
ENDIF
!
!
!-------------------------------------------------------------------------------------
!
!*      5.     Half prognostic fields
!
ALLOCATE(TGD%CUR%XRESA(SIZE(TGDPE%CUR%XLAI,1)))
TGD%CUR%XRESA = 100.
!
!-------------------------------------------------------------------------------------
!
!*      6.     Isba-Ags prognostic fields
!
IF (TVG%CPHOTO /= 'NON') THEN
!
   ALLOCATE(TGD%CUR%XAN(SIZE(TGDPE%CUR%XLAI,1)))
   TGD%CUR%XAN = 0.
!
   ALLOCATE(TGD%CUR%XANDAY(SIZE(TGDPE%CUR%XLAI,1)))
   TGD%CUR%XANDAY = 0.
!
   ALLOCATE(TGD%CUR%XANFM(SIZE(TGDPE%CUR%XLAI,1)))
   TGD%CUR%XANFM = XANFMINIT
!
   ALLOCATE(TGD%CUR%XLE(SIZE(TGDPE%CUR%XLAI,1)))
   TGD%CUR%XLE = 0.
!
ENDIF
!
IF (TVG%CPHOTO == 'AGS' .OR. TVG%CPHOTO == 'AST') THEN
!
   ALLOCATE(TGD%CUR%XBIOMASS(SIZE(TGDPE%CUR%XLAI,1),TVG%NNBIOMASS))
   TGD%CUR%XBIOMASS(:,1) = 0.
!
   ALLOCATE(TGD%CUR%XRESP_BIOMASS(SIZE(TGDPE%CUR%XLAI,1),TVG%NNBIOMASS))
   TGD%CUR%XRESP_BIOMASS(:,:) = 0.
!
ELSEIF (TVG%CPHOTO == 'LAI' .OR. TVG%CPHOTO == 'LST') THEN
!
   ALLOCATE(TGD%CUR%XBIOMASS(SIZE(TGDPE%CUR%XLAI,1),TVG%NNBIOMASS))
   TGD%CUR%XBIOMASS(:,1) = TGDPE%CUR%XLAI(:) * TGDP%XBSLAI(:)
!
   ALLOCATE(TGD%CUR%XRESP_BIOMASS(SIZE(TGDPE%CUR%XLAI,1),TVG%NNBIOMASS))
   TGD%CUR%XRESP_BIOMASS(:,:) = 0.
!
ELSEIF (TVG%CPHOTO == 'NIT' .OR. TVG%CPHOTO == 'NCB') THEN
!
   ALLOCATE(TGD%CUR%XBIOMASS(SIZE(TGDPE%CUR%XLAI,1),TVG%NNBIOMASS))
   TGD%CUR%XBIOMASS(:,1) = TGDPE%CUR%XLAI(:) * TGDP%XBSLAI_NITRO(:)
   TGD%CUR%XBIOMASS(:,2) = MAX( 0., (TGD%CUR%XBIOMASS(:,1)/ (XCC_NIT/10.**XCA_NIT))  &
                              **(1.0/(1.0-XCA_NIT)) - TGD%CUR%XBIOMASS(:,1) )  
   TGD%CUR%XBIOMASS(:,3:TVG%NNBIOMASS) = 0.
!
   ALLOCATE(TGD%CUR%XRESP_BIOMASS(SIZE(TGDPE%CUR%XLAI,1),TVG%NNBIOMASS))
   TGD%CUR%XRESP_BIOMASS(:,:) = 0.
!
ENDIF
!
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PREP_TEB_GARDEN',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE PREP_TEB_GARDEN
