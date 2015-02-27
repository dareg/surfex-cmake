!     #########
      SUBROUTINE WRITE_DIAG_SEB_SEAICE_n(HPROGRAM)
!     #################################
!
!!****  *WRITE_DIAG_SEB_SEAICE_n* - write the seaice diagnostic fields
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
!!      S.Senesi                *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2014
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SFX_OASIS,      ONLY : LCPL_SEAICE
!
USE MODD_DIAG_SURF_ATM_n,ONLY : LRESET_BUDGETC
!
USE MODD_DIAG_SEAICE_n
!
USE MODD_SEAFLUX_n,  ONLY : XTICE, XICE_ALB, CSEAICE_SCHEME, LHANDLE_SIC
!
USE MODD_DIAG_SEAFLUX_n,ONLY : N2M, LRAD_BUDGET, LSURF_BUDGET,          &
                                 LCOEF, LSURF_VARS, LSURF_BUDGETC,      &
                                 XT2M_ICE, XQ2M_ICE, XHU2M_ICE,         &
                                 XZON10M_ICE, XMER10M_ICE, XWIND10M_ICE,&
                                 XRN_ICE, XH_ICE, XGFLUX_ICE, XRI_ICE,  &
                                 XCD_ICE, XCH_ICE,                      &
                                 XZ0_ICE, XZ0H_ICE, XQS_ICE, XSWU_ICE,  &
                                 XLWU_ICE, XSWBU_ICE, XFMU_ICE, XFMV_ICE,&
                                 XRNC_ICE, XHC_ICE, XGFLUXC_ICE,        &
                                 XSWUC_ICE, XLWUC_ICE, XFMUC_ICE,       &
                                 XFMVC_ICE, XLE_ICE, XLEC_ICE
USE MODI_INIT_IO_SURF_n
USE MODI_WRITE_SURF
USE MODI_END_IO_SURF_n
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
CHARACTER(LEN=12) :: YRECFM         ! Name of the article to be read
CHARACTER(LEN=100):: YCOMMENT       ! Comment string
CHARACTER(LEN=2)  :: YNUM
INTEGER           :: JSV, JSW
!
REAL(KIND=JPRB)   :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!         Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_SEB_SEAICE_N',0,ZHOOK_HANDLE)
!
CALL INIT_IO_SURF_n(HPROGRAM,'SEA   ','SEAFLX','WRITE')
!
IF(LCPL_SEAICE.OR.LHANDLE_SIC)THEN      
!
  YCOMMENT='Sea-ice temperature (K)'
  CALL WRITE_SURF(HPROGRAM,'TSICE',XTICE(:),IRESP,YCOMMENT)
!
  YCOMMENT='Sea-ice albedo (-)'
  CALL WRITE_SURF(HPROGRAM,'IALB',XICE_ALB(:),IRESP,YCOMMENT)
!
ENDIF
!
IF (TRIM(CSEAICE_SCHEME) == 'GELATO') THEN 
    YCOMMENT='Sea-ice thickness (m)'
    CALL WRITE_SURF(HPROGRAM,'SIT',XSIT(:),IRESP,YCOMMENT)
    !
    YCOMMENT='Sea-ice snow depth (m)'
    CALL WRITE_SURF(HPROGRAM,'SND',XSND(:),IRESP,YCOMMENT)
    !
    YCOMMENT='Sea mixed layer temp for Glt (K)'
    CALL WRITE_SURF(HPROGRAM,'SIMLT',XMLT(:),IRESP,YCOMMENT)
    !
ENDIF
!
!
!*       8.2.     Richardson number :
!               -----------------
IF (N2M>=1) THEN
   !
   YRECFM='RI_SEAICE'
   YCOMMENT='X_Y_'//YRECFM
   !
   CALL WRITE_SURF(HPROGRAM,YRECFM,XRI_ICE(:),IRESP,HCOMMENT=YCOMMENT)
   !
END IF
!
!*       8.3     Energy fluxes :
!               -------------
!
IF (LSURF_BUDGET) THEN

   YRECFM='RN_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
   !
   CALL WRITE_SURF(HPROGRAM,YRECFM,XRN_ICE(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='H_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
   !
   CALL WRITE_SURF(HPROGRAM,YRECFM,XH_ICE(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='LE_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
   !
   CALL WRITE_SURF(HPROGRAM,YRECFM,XLE_ICE(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='GFLX_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
   !
   CALL WRITE_SURF(HPROGRAM,YRECFM,XGFLUX_ICE(:),IRESP,HCOMMENT=YCOMMENT)
   !
   IF (LRAD_BUDGET) THEN
      !
      YRECFM='SWU_SEAICE'
      YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
      !
      CALL WRITE_SURF(HPROGRAM,YRECFM,XSWU_ICE(:),IRESP,HCOMMENT=YCOMMENT)
      !
      YRECFM='LWU_SEAICE'
      YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
      !
      CALL WRITE_SURF(HPROGRAM,YRECFM,XLWU_ICE(:),IRESP,HCOMMENT=YCOMMENT)
      !
      DO JSW=1, SIZE(XSWBU_ICE,2)
         YNUM=ACHAR(48+JSW)
         !
         YRECFM='SWU_SEAICE_'//YNUM
         YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
         !
         CALL WRITE_SURF(HPROGRAM,YRECFM,XSWBU_ICE(:,JSW),IRESP,HCOMMENT=YCOMMENT)
         !
      ENDDO
      !
   ENDIF
   !
   YRECFM='FMU_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (kg/ms2)'
   !
   CALL WRITE_SURF(HPROGRAM,YRECFM,XFMU_ICE(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='FMV_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (kg/ms2)'
   !
   CALL WRITE_SURF(HPROGRAM,YRECFM,XFMV_ICE(:),IRESP,HCOMMENT=YCOMMENT)
   !
END IF
!
IF (LSURF_BUDGETC) THEN
   !
   YRECFM='RNC_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
   !
   CALL WRITE_SURF(HPROGRAM,YRECFM,XRNC_ICE(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='HC_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
   !
   CALL WRITE_SURF(HPROGRAM,YRECFM,XHC_ICE(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='LEC_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
   !
   CALL WRITE_SURF(HPROGRAM,YRECFM,XLEC_ICE(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='GFLXC_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
   !
   CALL WRITE_SURF(HPROGRAM,YRECFM,XGFLUXC_ICE(:),IRESP,HCOMMENT=YCOMMENT)
   IF (LRAD_BUDGET .OR. (LSURF_BUDGETC .AND. .NOT.LRESET_BUDGETC)) THEN
      !
      YRECFM='SWUC_SEAICE'
      YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
      !
      CALL WRITE_SURF(HPROGRAM,YRECFM,XSWUC_ICE(:),IRESP,HCOMMENT=YCOMMENT)
      !
      YRECFM='LWUC_SEAICE'
      YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
      !
      CALL WRITE_SURF(HPROGRAM,YRECFM,XLWUC_ICE(:),IRESP,HCOMMENT=YCOMMENT)
      !
   ENDIF
   !
   YRECFM='FMUC_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (kg/ms)'
   !
   CALL WRITE_SURF(HPROGRAM,YRECFM,XFMUC_ICE(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='FMVC_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (kg/ms)'
   !
   CALL WRITE_SURF(HPROGRAM,YRECFM,XFMVC_ICE(:),IRESP,HCOMMENT=YCOMMENT)
   !
END IF
!
!*       8.4     transfer coefficients
!               ---------------------
!
IF (LCOEF) THEN
   !
   YRECFM='CD_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (W/s2)'
   !
   CALL WRITE_SURF(HPROGRAM,YRECFM,XCD_ICE(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='CH_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (W/s)'
   !
   CALL WRITE_SURF(HPROGRAM,YRECFM,XCH_ICE(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='Z0_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (M)'
   !
   CALL WRITE_SURF(HPROGRAM,YRECFM,XZ0_ICE(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='Z0H_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (M)'
   !
   CALL WRITE_SURF(HPROGRAM,YRECFM,XZ0H_ICE(:),IRESP,HCOMMENT=YCOMMENT)
   !
END IF
!
!
!*       8.5     Surface humidity
!               ----------------
!
IF (LSURF_VARS) THEN
   YRECFM='QS_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (KG/KG)'
   !
   CALL WRITE_SURF(HPROGRAM,YRECFM,XQS_ICE(:),IRESP,HCOMMENT=YCOMMENT)
   !
ENDIF
!

!
!*       8.6.     parameters at 2 and 10 meters :
!               -----------------------------
!
IF (N2M>=1) THEN
   !
   YRECFM='T2M_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (K)'
   !
   CALL WRITE_SURF(HPROGRAM,YRECFM,XT2M_ICE(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='Q2M_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (KG/KG)'
   !
   CALL WRITE_SURF(HPROGRAM,YRECFM,XQ2M_ICE(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='HU2M_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (-)'
   !
   CALL WRITE_SURF(HPROGRAM,YRECFM,XHU2M_ICE(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='ZON10M_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (M/S)'
   !
   CALL WRITE_SURF(HPROGRAM,YRECFM,XZON10M_ICE(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='MER10M_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (M/S)'
   !
   CALL WRITE_SURF(HPROGRAM,YRECFM,XMER10M_ICE(:),IRESP,HCOMMENT=YCOMMENT)
   !
   YRECFM='W10M_SEAICE'
   YCOMMENT='X_Y_'//YRECFM//' (M/S)'
   !
   CALL WRITE_SURF(HPROGRAM,YRECFM,XWIND10M_ICE(:),IRESP,HCOMMENT=YCOMMENT)
   !
END IF
!
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)

IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_SEB_SEAICE_N',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE WRITE_DIAG_SEB_SEAICE_n
