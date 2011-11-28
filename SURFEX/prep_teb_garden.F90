!     #########
SUBROUTINE PREP_TEB_GARDEN(HPROGRAM,HATMFILE,HATMFILETYPE,OCANOPY)
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
USE MODD_READ_NAMELIST,  ONLY : LNAM_READ
USE MODD_SURF_ATM,       ONLY : LVERTSHIFT
USE MODD_DATA_COVER_PAR, ONLY : NVT_SNOW
!
USE MODD_TEB_GARDEN_n,   ONLY : TSNOW, XRESA, XLAI, CPHOTO, CRESPSL,          &
                                  XAN, XANFM, XANDAY, XLE,                      &
                                  NNBIOMASS, NNLITTER, NNLITTLEVS, NNSOILCARB,  &
                                  XBSLAI, XBSLAI_NITRO, XBIOMASS, XRESP_BIOMASS,&
                                  XLITTER, XSOILCARB, XLIGNIN_STRUC,                                            &
                                  NPATCH, XWSAT, XWG, XWGI, CISBA, XTG, LCANOPY,&
                                  XPATCH, XVEGTYPE_PATCH  
!
USE MODD_DEEPSOIL,    ONLY : LPHYSDOMC
USE MODD_CSTS,        ONLY : XTT
USE MODD_SNOW_PAR,    ONLY : XZ0SN
USE MODD_ISBA_PAR,    ONLY : XWGMIN
!
USE MODD_TEB_GRID_n,  ONLY : CGRID, XGRID_PAR, XLAT, XLON
USE MODD_CO2V_PAR,    ONLY : XANFMINIT, XCA_NIT, XCC_NIT
USE MODD_SURF_PAR,    ONLY : XUNDEF
USE MODD_PREP,        ONLY : XZS_LS
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
LOGICAL,            INTENT(IN)  :: OCANOPY     ! option canopy
!
!
!*      0.2    declarations of local variables
!
INTEGER :: ILUOUT
INTEGER :: JP
LOGICAL :: GFOUND         ! Return code when searching namelist
INTEGER :: ILUNAM         ! logical unit of namelist file
INTEGER :: ISNOW          ! patch number where permanent snow is
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
CALL PREP_HOR_TEB_GARDEN_FIELD(HPROGRAM,'WG     ',HATMFILE,HATMFILETYPE)
!
!*      2.2    Soil ice reservoirs
!
CALL PREP_HOR_TEB_GARDEN_FIELD(HPROGRAM,'WGI    ',HATMFILE,HATMFILETYPE)
!
!*      2.3    Leaves interception water reservoir
!
CALL PREP_HOR_TEB_GARDEN_FIELD(HPROGRAM,'WR     ',HATMFILE,HATMFILETYPE)
!
!*      2.4    Temperature profile
!
CALL PREP_HOR_TEB_GARDEN_FIELD(HPROGRAM,'TG     ',HATMFILE,HATMFILETYPE)
!
!*      2.5    Snow variables
!
CALL PREP_HOR_TEB_GARDEN_FIELD(HPROGRAM,'SN_VEG ',HATMFILE,HATMFILETYPE)

!
!*      2.6    LAI
!
CALL PREP_HOR_TEB_GARDEN_FIELD(HPROGRAM,'LAI    ',HATMFILE,HATMFILETYPE)
!
!-------------------------------------------------------------------------------------
!
!*      3.    Physical limitation: 
!
! If whole ice reservoir is empty (grib from ecmwf case) and surface temperature is
! lower than -10°C, then ice content is maximum and water content minimum
!
DO JP=1,NPATCH
   IF (ALL(XWGI(:,:,JP)==0.)) THEN
      WHERE(XTG(:,1:SIZE(XWG,2),JP) < XTT-10.)
         XWGI(:,:,JP) = XWSAT(:,:)-XWGMIN
         XWG (:,:,JP) = XWGMIN
      END WHERE
   ENDIF
ENDDO
!
! No ice for force restore third layer:
IF (CISBA == '3-L') THEN
   DO JP=1,NPATCH
      WHERE(XWG(:,3,JP) /= XUNDEF)
        XWG(:,3,JP)  = MIN(XWG(:,3,JP)+XWGI(:,3,JP),XWSAT(:,3))
        XWGI(:,3,JP) = 0.
      END WHERE
   ENDDO
ENDIF
!
! Total water content should not exceed saturation:
DO JP=1,NPATCH
   WHERE(XWG(:,:,JP) /= XUNDEF .AND. (XWG(:,:,JP) + XWGI(:,:,JP)) > XWSAT(:,:) )
      XWGI(:,:,JP) = XWSAT(:,:) - XWG(:,:,JP)
   END WHERE
ENDDO
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
ALLOCATE(XRESA(SIZE(XLAI,1),SIZE(XLAI,2)))
XRESA = 100.
!
!-------------------------------------------------------------------------------------
!
!*      6.     Isba-Ags prognostic fields
!
IF (CPHOTO /= 'NON') THEN
!
   ALLOCATE(XAN(SIZE(XLAI,1),SIZE(XLAI,2)))
   XAN = 0.
!
   ALLOCATE(XANDAY(SIZE(XLAI,1),SIZE(XLAI,2)))
   XANDAY = 0.
!
   ALLOCATE(XANFM(SIZE(XLAI,1),SIZE(XLAI,2)))
   XANFM = XANFMINIT
!
   ALLOCATE(XLE(SIZE(XLAI,1),SIZE(XLAI,2)))
   XLE = 0.
!
ENDIF
!
IF (CPHOTO == 'AGS' .OR. CPHOTO == 'AST') THEN
!
   ALLOCATE(XBIOMASS(SIZE(XLAI,1),NNBIOMASS,SIZE(XLAI,2)))
   XBIOMASS(:,1,:) = 0.
!
   ALLOCATE(XRESP_BIOMASS(SIZE(XLAI,1),NNBIOMASS,SIZE(XLAI,2)))
   XRESP_BIOMASS(:,:,:) = 0.
!
ELSEIF (CPHOTO == 'LAI' .OR. CPHOTO == 'LST') THEN
!
   ALLOCATE(XBIOMASS(SIZE(XLAI,1),NNBIOMASS,SIZE(XLAI,2)))
   XBIOMASS(:,1,:) = XLAI(:,:) * XBSLAI(:,:)
!
   ALLOCATE(XRESP_BIOMASS(SIZE(XLAI,1),NNBIOMASS,SIZE(XLAI,2)))
   XRESP_BIOMASS(:,:,:) = 0.
!
ELSEIF (CPHOTO == 'NIT' .OR. CPHOTO == 'NCB') THEN
!
   ALLOCATE(XBIOMASS(SIZE(XLAI,1),NNBIOMASS,SIZE(XLAI,2)))
   XBIOMASS(:,1,:) = XLAI(:,:) * XBSLAI_NITRO(:,:)
   XBIOMASS(:,2,:) = MAX( 0., (XBIOMASS(:,1,:)/ (XCC_NIT/10.**XCA_NIT))  &
                              **(1.0/(1.0-XCA_NIT)) - XBIOMASS(:,1,:) )  
   XBIOMASS(:,3:NNBIOMASS,:) = 0.
!
   ALLOCATE(XRESP_BIOMASS(SIZE(XLAI,1),NNBIOMASS,SIZE(XLAI,2)))
   XRESP_BIOMASS(:,:,:) = 0.
!
ENDIF
!
!-------------------------------------------------------------------------------------
!
!*      7.     Isba-CC prognostic fields
!
IF (CRESPSL == 'CNT') THEN
!
   ALLOCATE(XLITTER(SIZE(XLAI,1),NNLITTER,NNLITTLEVS,SIZE(XLAI,2)))
   XLITTER(:,:,:,:) = 0.
!
   ALLOCATE(XSOILCARB(SIZE(XLAI,1),NNSOILCARB,SIZE(XLAI,2)))
   XSOILCARB(:,:,:) = 0.
!
   ALLOCATE(XLIGNIN_STRUC(SIZE(XLAI,1),NNLITTLEVS,SIZE(XLAI,2)))
   XLIGNIN_STRUC(:,:,:) = 0.
!
ENDIF
IF (LHOOK) CALL DR_HOOK('PREP_TEB_GARDEN',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE PREP_TEB_GARDEN
