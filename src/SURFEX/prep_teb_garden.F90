!     #########
SUBROUTINE PREP_TEB_GARDEN (DTCO, UG, U, USS, TG, TOP, GDM, &
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
USE MODD_SURF_ATM_SSO_n, ONLY : SURF_ATM_SSO_t
USE MODD_GRID_n, ONLY : GRID_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_SURFEX_n, ONLY : TEB_GARDEN_MODEL_t
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
TYPE(SURF_ATM_SSO_t), INTENT(INOUT) :: USS
TYPE(GRID_t), INTENT(INOUT) :: TG
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_GARDEN_MODEL_t), INTENT(INOUT) :: GDM
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
 CALL PREP_HOR_TEB_GARDEN_FIELD(DTCO, UG, U, USS, &
                                GDM%TV%R, GDM%TV%O, GDM%TV%M, GDM%TV%IP, TG, TOP, &
                                HPROGRAM,'WG     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)
!
!*      2.2    Soil ice reservoirs
!
 CALL PREP_HOR_TEB_GARDEN_FIELD(DTCO, UG, U, USS, &
                                GDM%TV%R, GDM%TV%O, GDM%TV%M, GDM%TV%IP, TG, TOP, &         
                                HPROGRAM,'WGI    ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)
!
!*      2.3    Leaves interception water reservoir
!
 CALL PREP_HOR_TEB_GARDEN_FIELD(DTCO, UG, U, USS, &
                                GDM%TV%R, GDM%TV%O, GDM%TV%M, GDM%TV%IP, TG, TOP, &
                                HPROGRAM,'WR     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)
!
!*      2.4    Temperature profile
!
 CALL PREP_HOR_TEB_GARDEN_FIELD(DTCO, UG, U, USS, &
                                GDM%TV%R, GDM%TV%O, GDM%TV%M, GDM%TV%IP, TG, TOP, &
                                HPROGRAM,'TG     ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)
!
!*      2.5    Snow variables
!
 CALL PREP_HOR_TEB_GARDEN_FIELD(DTCO, UG, U, USS, &
                                GDM%TV%R, GDM%TV%O, GDM%TV%M, GDM%TV%IP, TG, TOP, &
                                HPROGRAM,'SN_VEG ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)

!
!*      2.6    LAI
!
IF (GDM%TV%O%CPHOTO/='NON' .AND. GDM%TV%O%CPHOTO/='AGS' .AND. GDM%TV%O%CPHOTO/='LST')  &
 CALL PREP_HOR_TEB_GARDEN_FIELD(DTCO, UG, U, USS, &
                                GDM%TV%R, GDM%TV%O, GDM%TV%M, GDM%TV%IP, TG, TOP, &
                                HPROGRAM,'LAI    ',HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE,KPATCH)
!
!-------------------------------------------------------------------------------------
!
!*      3.    Physical limitation: 
!
! If whole ice reservoir is empty (grib from ecmwf case) and surface temperature is
! lower than -10C, then ice content is maximum and water content minimum
!
IF (ALL(GDM%TV%R%CUR%XWGI(:,:,1)==0.)) THEN
   WHERE(GDM%TV%R%CUR%XTG(:,1:SIZE(GDM%TV%R%CUR%XWG,2),1) < XTT-10.)
       GDM%TV%R%CUR%XWGI(:,:,1) = GDM%TV%IP%XWSAT(:,:)-XWGMIN
       GDM%TV%R%CUR%XWG (:,:,1) = XWGMIN
   END WHERE
ENDIF
!
! No ice for force restore third layer:
IF (GDM%TV%O%CISBA == '3-L') THEN
      WHERE(GDM%TV%R%CUR%XWG(:,3,1)/=XUNDEF.AND.GDM%TV%R%CUR%XWGI(:,3,1)/=XUNDEF)
        GDM%TV%R%CUR%XWG(:,3,1)  = MIN(GDM%TV%R%CUR%XWG(:,3,1)+GDM%TV%R%CUR%XWGI(:,3,1),GDM%TV%IP%XWSAT(:,3))
        GDM%TV%R%CUR%XWGI(:,3,1) = 0.
      END WHERE
ENDIF
!
! Total water content should not exceed saturation:
WHERE(GDM%TV%R%CUR%XWG(:,:,1) /= XUNDEF .AND. &
                  (GDM%TV%R%CUR%XWG(:,:,1) + GDM%TV%R%CUR%XWGI(:,:,1)) > GDM%TV%IP%XWSAT(:,:) )
     GDM%TV%R%CUR%XWGI(:,:,1) = GDM%TV%IP%XWSAT(:,:) - GDM%TV%R%CUR%XWG(:,:,1)
END WHERE
!
!-------------------------------------------------------------------------------------
!
!*      3.     Vertical interpolations of all variables
!
IF(LVERTSHIFT)THEN
  CALL PREP_VER_TEB_GARDEN(GDM%TV%R, GDM%TV%O, GDM%TV%P, TOP, GDM%TV%M%X%XDG(:,:,1))
ENDIF
!
!
!-------------------------------------------------------------------------------------
!
!*      5.     Half prognostic fields
!
ALLOCATE(GDM%TV%R%CUR%XRESA(SIZE(GDM%TV%M%T%CUR%XLAI,1),1))
GDM%TV%R%CUR%XRESA = 100.
!
!-------------------------------------------------------------------------------------
!
!*      6.     Isba-Ags prognostic fields
!
IF (GDM%TV%O%CPHOTO /= 'NON') THEN
!
   ALLOCATE(GDM%TV%R%CUR%XAN(SIZE(GDM%TV%M%T%CUR%XLAI,1),1))
   GDM%TV%R%CUR%XAN = 0.
!
   ALLOCATE(GDM%TV%R%CUR%XANDAY(SIZE(GDM%TV%M%T%CUR%XLAI,1),1))
   GDM%TV%R%CUR%XANDAY = 0.
!
   ALLOCATE(GDM%TV%R%CUR%XANFM(SIZE(GDM%TV%M%T%CUR%XLAI,1),1))
   GDM%TV%R%CUR%XANFM = XANFMINIT
!
   ALLOCATE(GDM%TV%R%CUR%XLE(SIZE(GDM%TV%M%T%CUR%XLAI,1),1))
   GDM%TV%R%CUR%XLE = 0.
!
ENDIF
!
IF (GDM%TV%O%CPHOTO == 'AGS' .OR. GDM%TV%O%CPHOTO == 'AST') THEN
!
   ALLOCATE(GDM%TV%R%CUR%XBIOMASS(SIZE(GDM%TV%M%T%CUR%XLAI,1),GDM%TV%O%NNBIOMASS,1))
   GDM%TV%R%CUR%XBIOMASS(:,1,1) = 0.
!
   ALLOCATE(GDM%TV%R%CUR%XRESP_BIOMASS(SIZE(GDM%TV%M%T%CUR%XLAI,1),GDM%TV%O%NNBIOMASS,1))
   GDM%TV%R%CUR%XRESP_BIOMASS(:,:,1) = 0.
!
ELSEIF (GDM%TV%O%CPHOTO == 'LAI' .OR. GDM%TV%O%CPHOTO == 'LST') THEN
!
   ALLOCATE(GDM%TV%R%CUR%XBIOMASS(SIZE(GDM%TV%M%T%CUR%XLAI,1),GDM%TV%O%NNBIOMASS,1))
   GDM%TV%R%CUR%XBIOMASS(:,1,1) = GDM%TV%M%T%CUR%XLAI(:,1) * GDM%TV%M%T%CUR%XBSLAI(:,1)
!
   ALLOCATE(GDM%TV%R%CUR%XRESP_BIOMASS(SIZE(GDM%TV%M%T%CUR%XLAI,1),GDM%TV%O%NNBIOMASS,1))
   GDM%TV%R%CUR%XRESP_BIOMASS(:,:,:) = 0.
!
ELSEIF (GDM%TV%O%CPHOTO == 'NIT' .OR. GDM%TV%O%CPHOTO == 'NCB') THEN
!
   ALLOCATE(GDM%TV%R%CUR%XBIOMASS(SIZE(GDM%TV%M%T%CUR%XLAI,1),GDM%TV%O%NNBIOMASS,1))
   GDM%TV%R%CUR%XBIOMASS(:,1,1) = GDM%TV%M%T%CUR%XLAI(:,1) * GDM%TV%IP%XBSLAI_NITRO(:,1)
   GDM%TV%R%CUR%XBIOMASS(:,2,1) = MAX( 0., (GDM%TV%R%CUR%XBIOMASS(:,1,1)/ (XCC_NIT/10.**XCA_NIT))  &
                              **(1.0/(1.0-XCA_NIT)) - GDM%TV%R%CUR%XBIOMASS(:,1,1) )  
   GDM%TV%R%CUR%XBIOMASS(:,3:GDM%TV%O%NNBIOMASS,1) = 0.
!
   ALLOCATE(GDM%TV%R%CUR%XRESP_BIOMASS(SIZE(GDM%TV%M%T%CUR%XLAI,1),GDM%TV%O%NNBIOMASS,1))
   GDM%TV%R%CUR%XRESP_BIOMASS(:,:,:) = 0.
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
