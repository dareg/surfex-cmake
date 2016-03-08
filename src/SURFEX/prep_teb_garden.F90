!     #########
SUBROUTINE PREP_TEB_GARDEN (DTCO, UG, U, USS, TG, TOP, GDR, GDO, GDMX, GDMT, GDIP, &
                            HPROGRAM,HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)
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
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SSO_n, ONLY : SSO_t
USE MODD_GRID_n, ONLY : GRID_t
!
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
!
USE MODD_ISBA_n, ONLY : ISBA_PROG_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
USE MODD_ISBA_PARAM_n, ONLY : ISBA_PARAM_FIX_t, ISBA_PARAM_TIME_t
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t
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
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SSO_t), INTENT(INOUT) :: USS
TYPE(GRID_t), INTENT(INOUT) :: TG
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
!
TYPE(ISBA_PROG_t), INTENT(INOUT) :: GDR
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: GDO
TYPE(ISBA_PARAM_FIX_t), INTENT(INOUT) :: GDMX
TYPE(ISBA_PARAM_TIME_t), INTENT(INOUT) :: GDMT
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: GDIP
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=28),  INTENT(IN)  :: HATMFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),   INTENT(IN)  :: HATMFILETYPE! type of the Atmospheric file
 CHARACTER(LEN=28),  INTENT(IN)  :: HPGDFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),   INTENT(IN)  :: HPGDFILETYPE! type of the Atmospheric file
!
INTEGER,            INTENT(IN)  :: KPATCH
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
 CALL PREP_HOR_TEB_GARDEN_FIELD(DTCO, UG, U, USS, GDR, GDO, GDMX, GDMT%XLAI, GDIP, TG, TOP, &
                                HPROGRAM,'WG     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)
!
!*      2.2    Soil ice reservoirs
!
 CALL PREP_HOR_TEB_GARDEN_FIELD(DTCO, UG, U, USS, GDR, GDO, GDMX, GDMT%XLAI, GDIP, TG, TOP, &         
                                HPROGRAM,'WGI    ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)
!
!*      2.3    Leaves interception water reservoir
!
 CALL PREP_HOR_TEB_GARDEN_FIELD(DTCO, UG, U, USS, GDR, GDO, GDMX, GDMT%XLAI, GDIP, TG, TOP, &
                                HPROGRAM,'WR     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)
!
!*      2.4    Temperature profile
!
 CALL PREP_HOR_TEB_GARDEN_FIELD(DTCO, UG, U, USS, GDR, GDO, GDMX, GDMT%XLAI, GDIP, TG, TOP, &
                                HPROGRAM,'TG     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)
!
!*      2.5    Snow variables
!
 CALL PREP_HOR_TEB_GARDEN_FIELD(DTCO, UG, U, USS, GDR, GDO, GDMX, GDMT%XLAI, GDIP, TG, TOP, &
                                HPROGRAM,'SN_VEG ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)

!
!*      2.6    LAI
!
IF (GDO%CPHOTO/='NON' .AND. GDO%CPHOTO/='AGS' .AND. GDO%CPHOTO/='LST')  &
 CALL PREP_HOR_TEB_GARDEN_FIELD(DTCO, UG, U, USS, GDR, GDO, GDMX, GDMT%XLAI, GDIP, TG, TOP, &
                                HPROGRAM,'LAI    ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)
!
!-------------------------------------------------------------------------------------
!
!*      3.    Physical limitation: 
!
! If whole ice reservoir is empty (grib from ecmwf case) and surface temperature is
! lower than -10C, then ice content is maximum and water content minimum
!
IF (ALL(GDR%XWGI(:,:,1)==0.)) THEN
   WHERE(GDR%XTG(:,1:SIZE(GDR%XWG,2),1) < XTT-10.)
       GDR%XWGI(:,:,1) = GDIP%XWSAT(:,:)-XWGMIN
       GDR%XWG (:,:,1) = XWGMIN
   END WHERE
ENDIF
!
! No ice for force restore third layer:
IF (GDO%CISBA == '3-L') THEN
      WHERE(GDR%XWG(:,3,1)/=XUNDEF.AND.GDR%XWGI(:,3,1)/=XUNDEF)
        GDR%XWG(:,3,1)  = MIN(GDR%XWG(:,3,1)+GDR%XWGI(:,3,1),GDIP%XWSAT(:,3))
        GDR%XWGI(:,3,1) = 0.
      END WHERE
ENDIF
!
! Total water content should not exceed saturation:
WHERE(GDR%XWG(:,:,1) /= XUNDEF .AND. &
                  (GDR%XWG(:,:,1) + GDR%XWGI(:,:,1)) > GDIP%XWSAT(:,:) )
     GDR%XWGI(:,:,1) = GDIP%XWSAT(:,:) - GDR%XWG(:,:,1)
END WHERE
!
!-------------------------------------------------------------------------------------
!
!*      3.     Vertical interpolations of all variables
!
IF(LVERTSHIFT)THEN
  CALL PREP_VER_TEB_GARDEN(GDR, GDO, TOP%XZS, GDMX%XDG(:,:,1))
ENDIF
!
!
!-------------------------------------------------------------------------------------
!
!*      5.     Half prognostic fields
!
ALLOCATE(GDR%XRESA(SIZE(GDMT%XLAI,1),1))
GDR%XRESA = 100.
!
!-------------------------------------------------------------------------------------
!
!*      6.     Isba-Ags prognostic fields
!
IF (GDO%CPHOTO /= 'NON') THEN
!
   ALLOCATE(GDR%XAN(SIZE(GDMT%XLAI,1),1))
   GDR%XAN = 0.
!
   ALLOCATE(GDR%XANDAY(SIZE(GDMT%XLAI,1),1))
   GDR%XANDAY = 0.
!
   ALLOCATE(GDR%XANFM(SIZE(GDMT%XLAI,1),1))
   GDR%XANFM = XANFMINIT
!
   ALLOCATE(GDR%XLE(SIZE(GDMT%XLAI,1),1))
   GDR%XLE = 0.
!
ENDIF
!
IF (GDO%CPHOTO == 'AGS' .OR. GDO%CPHOTO == 'AST') THEN
!
   ALLOCATE(GDR%XBIOMASS(SIZE(GDMT%XLAI,1),GDO%NNBIOMASS,1))
   GDR%XBIOMASS(:,1,1) = 0.
!
   ALLOCATE(GDR%XRESP_BIOMASS(SIZE(GDMT%XLAI,1),GDO%NNBIOMASS,1))
   GDR%XRESP_BIOMASS(:,:,1) = 0.
!
ELSEIF (GDO%CPHOTO == 'LAI' .OR. GDO%CPHOTO == 'LST') THEN
!
   ALLOCATE(GDR%XBIOMASS(SIZE(GDMT%XLAI,1),GDO%NNBIOMASS,1))
   GDR%XBIOMASS(:,1,1) = GDMT%XLAI(:,1) * GDMT%XBSLAI(:,1)
!
   ALLOCATE(GDR%XRESP_BIOMASS(SIZE(GDMT%XLAI,1),GDO%NNBIOMASS,1))
   GDR%XRESP_BIOMASS(:,:,:) = 0.
!
ELSEIF (GDO%CPHOTO == 'NIT' .OR. GDO%CPHOTO == 'NCB') THEN
!
   ALLOCATE(GDR%XBIOMASS(SIZE(GDMT%XLAI,1),GDO%NNBIOMASS,1))
   GDR%XBIOMASS(:,1,1) = GDMT%XLAI(:,1) * GDIP%XBSLAI_NITRO(:,1)
   GDR%XBIOMASS(:,2,1) = MAX( 0., (GDR%XBIOMASS(:,1,1)/ (XCC_NIT/10.**XCA_NIT))  &
                              **(1.0/(1.0-XCA_NIT)) - GDR%XBIOMASS(:,1,1) )  
   GDR%XBIOMASS(:,3:GDO%NNBIOMASS,1) = 0.
!
   ALLOCATE(GDR%XRESP_BIOMASS(SIZE(GDMT%XLAI,1),GDO%NNBIOMASS,1))
   GDR%XRESP_BIOMASS(:,:,:) = 0.
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
