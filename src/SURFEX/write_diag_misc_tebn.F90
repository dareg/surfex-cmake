!     #########
      SUBROUTINE WRITE_DIAG_MISC_TEB_n(HPROGRAM)
!     #################################
!
!!****  *WRITE_DIAG_MISC_TEB* - writes the TEB diagnostic fields
!!
!!    PURPOSE
!!    -------
!!
!!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	P. Le Moigne   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2004
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODI_INIT_IO_SURF_n
USE MODI_WRITE_SURF
USE MODI_END_IO_SURF_n
USE MODD_TEB_n,ONLY : XZ0_TOWN, XTI_BLD
USE MODD_DIAG_MISC_TEB_n,ONLY : LSURF_MISC_BUDGET, LSURF_BUDGETC,              &
                                  XQF_BLD,XFLX_BLD,                              &
                                  XQF_TOWN, XDQS_TOWN,XTI_BLD_EQ,XQF_BLDWFR,     &
                                  XTI_BLDWFR,                                    &
                                  XRN_ROAD, XH_ROAD, XLE_ROAD, XGFLUX_ROAD,      &
                                  XRN_WALL, XH_WALL, XGFLUX_WALL,                &
                                  XRN_ROOF, XH_ROOF, XLE_ROOF, XGFLUX_ROOF,      &
                                  XRUNOFF, XRUNOFFC,                             &
                                  XRN_GARDEN,XH_GARDEN,XLE_GARDEN,XGFLUX_GARDEN, &
                                  XRN_BLT,XH_BLT,XLE_BLT,XGFLUX_BLT,             &
                                  XABS_SW_ROOF ,XABS_SW_SNOW_ROOF,               &
                                  XABS_LW_ROOF ,XABS_LW_SNOW_ROOF,               &
                                  XABS_SW_ROAD ,XABS_SW_SNOW_ROAD,               &
                                  XABS_LW_ROAD ,XABS_LW_SNOW_ROAD,               &
                                  XABS_SW_WALL ,                                 &
                                  XABS_LW_WALL ,                                 &
                                  XABS_SW_GARDEN, XABS_LW_GARDEN  
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
CHARACTER(LEN=6),  INTENT(IN)  :: HPROGRAM ! program calling
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
CHARACTER(LEN=16) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=100):: YCOMMENT       ! Comment string
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!         Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_MISC_TEB_N',0,ZHOOK_HANDLE)
CALL INIT_IO_SURF_n(HPROGRAM,'TOWN  ','TEB   ','WRITE')
!
!-------------------------------------------------------------------------------
!
IF (LSURF_MISC_BUDGET) THEN
!
!*       Miscellaneous fields :
!        ----------------------
YRECFM='Z0_TOWN'
YCOMMENT='town roughness length'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XZ0_TOWN(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='XQF_BLD'
YCOMMENT='domestic heating'//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XQF_BLD(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='XQF_BLDWFR'
YCOMMENT='domestic heating'//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XQF_BLDWFR(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='XFLX_BLD'
YCOMMENT='heat flux from bld'//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XFLX_BLD(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='XQF_TOWN'
YCOMMENT='total anthropogenic heat'//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XQF_TOWN(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='XDQS_TOWN'
YCOMMENT='heat storage inside building'//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XDQS_TOWN(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='XTI_BLD_EQ'
YCOMMENT='building internal temperature'//' (oC)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XTI_BLD_EQ(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='XTI_BLDWFR'
YCOMMENT='building internal temperature'//' (oC)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XTI_BLDWFR(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='RN_ROAD'
YCOMMENT=' net radiation at road'//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XRN_ROAD(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='H_ROAD'
YCOMMENT='road sensible heat flux'//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XH_ROAD(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='LE_ROAD'
YCOMMENT='road latent heat flux'//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XLE_ROAD(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='GFLUX_ROAD'
YCOMMENT='net road conduction flux'//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XGFLUX_ROAD(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='RN_WALL'
YCOMMENT='net radiation for wall'//YRECFM//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XRN_WALL(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='H_WALL'
YCOMMENT='wall sensible heat flux'//YRECFM//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XH_WALL(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='GFLUX_WALL'
YCOMMENT='net wall conduction flux'//YRECFM//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XGFLUX_WALL(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='RN_ROOF'
YCOMMENT='net radiation for roof'//YRECFM//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XRN_ROOF(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='H_ROOF'
YCOMMENT='roof sensible heat flux'//YRECFM//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XH_ROOF(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='LE_ROOF'
YCOMMENT='roof latent heat flux'//YRECFM//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XLE_ROOF(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='GFLUX_ROOF'
YCOMMENT='net roof conduction flux'//YRECFM//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XGFLUX_ROOF(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='RUNOFF_TWN'
YCOMMENT='X_Y_'//YRECFM//' (kg/m2/s)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XRUNOFF(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='RN_GARDEN'
YCOMMENT='net radiation for GARDEN areas'//YRECFM//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XRN_GARDEN(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='H_GARDEN'
YCOMMENT='GARDEN area sensible heat flux'//YRECFM//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XH_GARDEN(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='LE_GARDEN'
YCOMMENT='GARDEN area latent heat flux'//YRECFM//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XLE_GARDEN(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='GFLUX_GARDEN'
YCOMMENT='net GARDEN area conduction flux'//YRECFM//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XGFLUX_GARDEN(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='RN_BLT'
YCOMMENT='net radiation for built surfaces'//YRECFM//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XRN_BLT(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='H_BLT'
YCOMMENT='built surface sensible heat flux'//YRECFM//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XH_BLT(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='LE_BLT'
YCOMMENT='built surface latent heat flux'//YRECFM//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XLE_BLT(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='GFLUX_BLT'
YCOMMENT='built surface conduction flux'//YRECFM//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XGFLUX_BLT(:),IRESP,HCOMMENT=YCOMMENT)
!
!
YRECFM='SWA_ROOF'
YCOMMENT='Sdown absorbed by roofs'//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XABS_SW_ROOF(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='SWA_SN_ROOF'
YCOMMENT='Sdown absorbed by snow on roofs'//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XABS_SW_SNOW_ROOF(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='LWA_ROOF'
YCOMMENT='Ldown absorbed by roofs'//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XABS_LW_ROOF(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='LWA_SN_ROOF'
YCOMMENT='Ldown absorbed by snow on roofs'//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XABS_LW_SNOW_ROOF(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='SWA_ROAD'
YCOMMENT='Sdown absorbed by roads'//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XABS_SW_ROAD(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='SWA_SN_ROAD'
YCOMMENT='Sdown absorbed by snow on roads'//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XABS_SW_SNOW_ROAD(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='LWA_ROAD'
YCOMMENT='Ldown absorbed by roads'//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XABS_LW_ROAD(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='LWA_SN_ROAD'
YCOMMENT='Ldown absorbed by snow on roads'//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XABS_LW_SNOW_ROAD(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='SWA_WALL'
YCOMMENT='Sdown absorbed by walls'//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XABS_SW_WALL(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='LWA_WALL'
YCOMMENT='Ldown absorbed by walls'//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XABS_LW_WALL(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='SWA_GARDEN'
YCOMMENT='Sdown absorbed by GARDEN areas'//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XABS_SW_GARDEN(:),IRESP,HCOMMENT=YCOMMENT)
!
YRECFM='LWA_GARDEN'
YCOMMENT='Ldown absorbed by GARDEN areas'//' (W/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XABS_LW_GARDEN(:),IRESP,HCOMMENT=YCOMMENT)
!
END IF
!
YRECFM='XTI_BLD'
YCOMMENT='building interior temperature'//YRECFM//' (oC)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XTI_BLD(:),IRESP,HCOMMENT=YCOMMENT)
!
!*       5.    Cumulated Energy fluxes 
!
IF (LSURF_BUDGETC) THEN
!
YRECFM='RUNOFFC_TWN'
YCOMMENT='X_Y_'//YRECFM//' (kg/m2)'
!
CALL WRITE_SURF(HPROGRAM,YRECFM,XRUNOFFC(:),IRESP,HCOMMENT=YCOMMENT)
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!         End of IO
!
CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_MISC_TEB_N',1,ZHOOK_HANDLE)
!
END SUBROUTINE WRITE_DIAG_MISC_TEB_n
