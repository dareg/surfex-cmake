!     #########
SUBROUTINE PREP_TEB_GARDEN(HPROGRAM,HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
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
USE MODI_PREP_HOR_TEB_GARDEN_FIELD
USE MODI_PREP_VER_TEB_GARDEN
!
USE MODD_SURF_ATM,       ONLY : LVERTSHIFT
!
USE MODD_TEB_VEG_n, ONLY : TVG => TEB_VEG_OPTIONS
USE MODD_TEB_GARDEN_PGD_EVOL_n, ONLY : TGDPE => TEB_GARDEN_PGD_EVOL
USE MODD_TEB_GARDEN_PGD_n, ONLY : TGDP => TEB_GARDEN_PGD
USE MODD_TEB_GARDEN_n, ONLY : TGD => TEB_GARDEN
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
 CALL PREP_HOR_TEB_GARDEN_FIELD(HPROGRAM,'WG     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!*      2.2    Soil ice reservoirs
!
 CALL PREP_HOR_TEB_GARDEN_FIELD(HPROGRAM,'WGI    ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!*      2.3    Leaves interception water reservoir
!
 CALL PREP_HOR_TEB_GARDEN_FIELD(HPROGRAM,'WR     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!*      2.4    Temperature profile
!
 CALL PREP_HOR_TEB_GARDEN_FIELD(HPROGRAM,'TG     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!*      2.5    Snow variables
!
 CALL PREP_HOR_TEB_GARDEN_FIELD(HPROGRAM,'SN_VEG ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)

!
!*      2.6    LAI
!
IF (TVG%CPHOTO/='NON' .AND. TVG%CPHOTO/='AGS' .AND. TVG%CPHOTO/='LST')  &
 CALL PREP_HOR_TEB_GARDEN_FIELD(HPROGRAM,'LAI    ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!
!-------------------------------------------------------------------------------------
!
!*      3.    Physical limitation: 
!
! If whole ice reservoir is empty (grib from ecmwf case) and surface temperature is
! lower than -10C, then ice content is maximum and water content minimum
!
IF (ALL(TGD%XWGI(:,:)==0.)) THEN
   WHERE(TGD%XTG(:,1:SIZE(TGD%XWG,2)) < XTT-10.)
       TGD%XWGI(:,:) = TGDP%XWSAT(:,:)-XWGMIN
       TGD%XWG (:,:) = XWGMIN
   END WHERE
ENDIF
!
! No ice for force restore third layer:
IF (TVG%CISBA == '3-L') THEN
      WHERE(TGD%XWG(:,3)/=XUNDEF.AND.TGD%XWGI(:,3)/=XUNDEF)
        TGD%XWG(:,3)  = MIN(TGD%XWG(:,3)+TGD%XWGI(:,3),TGDP%XWSAT(:,3))
        TGD%XWGI(:,3) = 0.
      END WHERE
ENDIF
!
! Total water content should not exceed saturation:
WHERE(TGD%XWG(:,:) /= XUNDEF .AND. (TGD%XWG(:,:) + TGD%XWGI(:,:)) > TGDP%XWSAT(:,:) )
     TGD%XWGI(:,:) = TGDP%XWSAT(:,:) - TGD%XWG(:,:)
END WHERE
!
!-------------------------------------------------------------------------------------
!
!*      3.     Vertical interpolations of all variables
!
IF(LVERTSHIFT)THEN
  CALL PREP_VER_TEB_GARDEN
ENDIF
!
!
!-------------------------------------------------------------------------------------
!
!*      5.     Half prognostic fields
!
ALLOCATE(TGD%XRESA(SIZE(TGDPE%XLAI,1)))
TGD%XRESA = 100.
!
!-------------------------------------------------------------------------------------
!
!*      6.     Isba-Ags prognostic fields
!
IF (TVG%CPHOTO /= 'NON') THEN
!
   ALLOCATE(TGD%XAN(SIZE(TGDPE%XLAI,1)))
   TGD%XAN = 0.
!
   ALLOCATE(TGD%XANDAY(SIZE(TGDPE%XLAI,1)))
   TGD%XANDAY = 0.
!
   ALLOCATE(TGD%XANFM(SIZE(TGDPE%XLAI,1)))
   TGD%XANFM = XANFMINIT
!
   ALLOCATE(TGD%XLE(SIZE(TGDPE%XLAI,1)))
   TGD%XLE = 0.
!
ENDIF
!
IF (TVG%CPHOTO == 'AGS' .OR. TVG%CPHOTO == 'AST') THEN
!
   ALLOCATE(TGD%XBIOMASS(SIZE(TGDPE%XLAI,1),TVG%NNBIOMASS))
   TGD%XBIOMASS(:,1) = 0.
!
   ALLOCATE(TGD%XRESP_BIOMASS(SIZE(TGDPE%XLAI,1),TVG%NNBIOMASS))
   TGD%XRESP_BIOMASS(:,:) = 0.
!
ELSEIF (TVG%CPHOTO == 'LAI' .OR. TVG%CPHOTO == 'LST') THEN
!
   ALLOCATE(TGD%XBIOMASS(SIZE(TGDPE%XLAI,1),TVG%NNBIOMASS))
   TGD%XBIOMASS(:,1) = TGDPE%XLAI(:) * TGDP%XBSLAI(:)
!
   ALLOCATE(TGD%XRESP_BIOMASS(SIZE(TGDPE%XLAI,1),TVG%NNBIOMASS))
   TGD%XRESP_BIOMASS(:,:) = 0.
!
ELSEIF (TVG%CPHOTO == 'NIT' .OR. TVG%CPHOTO == 'NCB') THEN
!
   ALLOCATE(TGD%XBIOMASS(SIZE(TGDPE%XLAI,1),TVG%NNBIOMASS))
   TGD%XBIOMASS(:,1) = TGDPE%XLAI(:) * TGDP%XBSLAI_NITRO(:)
   TGD%XBIOMASS(:,2) = MAX( 0., (TGD%XBIOMASS(:,1)/ (XCC_NIT/10.**XCA_NIT))  &
                              **(1.0/(1.0-XCA_NIT)) - TGD%XBIOMASS(:,1) )  
   TGD%XBIOMASS(:,3:TVG%NNBIOMASS) = 0.
!
   ALLOCATE(TGD%XRESP_BIOMASS(SIZE(TGDPE%XLAI,1),TVG%NNBIOMASS))
   TGD%XRESP_BIOMASS(:,:) = 0.
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
