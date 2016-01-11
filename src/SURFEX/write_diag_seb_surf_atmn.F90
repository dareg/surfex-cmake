!     #########
      SUBROUTINE WRITE_DIAG_SEB_SURF_ATM_n (DTCO, DGU, DGUC, U, UG, &
                                            HPROGRAM)
!     #################################
!
!!****  *WRITE_DIAG_SEB_SURF_ATM_n* - writes surface diagnostics
!!
!!    PURPOSE
!!    -------
!!
!!
!!**  METHOD
!!    ------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      Modified    01/2006 : sea flux parameterization.
!!      Modified    08/2009 : cumulated diag
!!      Juan        6/12/2011: parallel bug , remove local ANY(XAVG_ZON10M) test
!!      B. Decharme  06/13   Add QS, evap and sublimation diags
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DIAG_n, ONLY : DIAG_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
!
USE MODD_SURF_PAR, ONLY : XUNDEF
!
USE MODI_INIT_IO_SURF_n
USE MODI_WRITE_SURF
USE MODI_END_IO_SURF_n
USE MODI_SUM_ON_ALL_PROCS
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DIAG_t), INTENT(INOUT) :: DGU
TYPE(DIAG_t), INTENT(INOUT) :: DGUC
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
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
!
INTEGER           :: JSW
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!         Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_SEB_SURF_ATM_N',0,ZHOOK_HANDLE)
CALL INIT_IO_SURF_n(DTCO, DGU, U, &
                     HPROGRAM,'FULL  ','SURF  ','WRITE')
!
!
!*       1.     Richardson number :
!               -----------------
!
IF (DGU%N2M>=1) THEN
  !        
  YRECFM='RI'
  YCOMMENT='X_Y_'//YRECFM
  !
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XRI(:),IRESP,HCOMMENT=YCOMMENT)
  !
ENDIF
!
!*       2.     parameters at surface, 2 and 10 meters :
!               ----------------------------------------
!
IF (DGU%N2M>=1.OR.DGU%LSURF_BUDGET.OR.DGU%LSURF_BUDGETC) THEN
  !
  YRECFM='TS'
  YCOMMENT='X_Y_'//YRECFM//' (K)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XTS(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='TSRAD'
  YCOMMENT='X_Y_'//YRECFM//' (K)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XTRAD(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='EMIS'
  YCOMMENT='X_Y_'//YRECFM//' (-)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XEMIS(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='SFCO2'
  YCOMMENT='X_Y_'//YRECFM//' (M.kgCO2.S-1.kgAIR-1)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XSFCO2(:),IRESP,HCOMMENT=YCOMMENT)
  !
ENDIF
!
IF (DGU%N2M>=1) THEN
  !
  YRECFM='T2M'
  YCOMMENT='X_Y_'//YRECFM//' (K)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XT2M(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='T2MMIN'
  YCOMMENT='X_Y_'//YRECFM//' (K)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XT2M_MIN(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='T2MMAX'
  YCOMMENT='X_Y_'//YRECFM//' (K)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XT2M_MAX(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='Q2M'
  YCOMMENT='X_Y_'//YRECFM//' (KG/KG)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XQ2M(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='HU2M'
  YCOMMENT='X_Y_'//YRECFM//' (-)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XHU2M(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='HU2MMIN'
  YCOMMENT='X_Y_'//YRECFM//' (-)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XHU2M_MIN(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='HU2MMAX'
  YCOMMENT='X_Y_'//YRECFM//' (-)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XHU2M_MAX(:),IRESP,HCOMMENT=YCOMMENT)
  !
  IF ( SUM_ON_ALL_PROCS(HPROGRAM,UG%G%CGRID,DGU%XZON10M(:)/= XUNDEF) > 0. ) THEN
    !
    YRECFM='ZON10M'
    YCOMMENT='X_Y_'//YRECFM//' (M/S)'
    CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XZON10M(:),IRESP,HCOMMENT=YCOMMENT)
    !
    YRECFM='MER10M'
    YCOMMENT='X_Y_'//YRECFM//' (M/S)'
    CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XMER10M(:),IRESP,HCOMMENT=YCOMMENT)
    !
    YRECFM='W10M'
    YCOMMENT='X_Y_'//YRECFM//' (M/S)'
    CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XWIND10M(:),IRESP,HCOMMENT=YCOMMENT)
    !
    YRECFM='W10MMAX'
    YCOMMENT='X_Y_'//YRECFM//' (M/S)'
    CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XWIND10M_MAX(:),IRESP,HCOMMENT=YCOMMENT)
    !
  ENDIF
  !
  IF (DGU%L2M_MIN_ZS) THEN
    !
    YRECFM='T2M_MIN_ZS'
    YCOMMENT='X_Y_'//YRECFM//' (K)'
    CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XT2M_MIN_ZS(:),IRESP,HCOMMENT=YCOMMENT)
    !
    YRECFM='Q2M_MIN_ZS'
    YCOMMENT='X_Y_'//YRECFM//' (KG/KG)'
    CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XQ2M_MIN_ZS(:),IRESP,HCOMMENT=YCOMMENT)
    !
    YRECFM='HU2M_MIN_ZS'
    YCOMMENT='X_Y_'//YRECFM//' (KG/KG)'
    CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XHU2M_MIN_ZS(:),IRESP,HCOMMENT=YCOMMENT)
    !
  END IF
  !
END IF
!
!*       3.     Energy fluxes :
!               -------------
!
IF (DGU%LSURF_BUDGET) THEN
  !
  YRECFM='RN'
  YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XRN(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='H'
  YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XH(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='LE'
  YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XLE(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='LEI'
  YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XLEI(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='GFLUX'
  YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XGFLUX(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='EVAP'
  YCOMMENT='X_Y_'//YRECFM//' (kg/m2/s)'
  !
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XEVAP(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='SUBL'
  YCOMMENT='X_Y_'//YRECFM//' (kg/m2/s)'
  !
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XSUBL(:),IRESP,HCOMMENT=YCOMMENT)
  !
  IF (DGU%LRAD_BUDGET) THEN
    !         
    YRECFM='SWD'
    YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
    CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XSWD(:),IRESP,HCOMMENT=YCOMMENT)
    !
    YRECFM='SWU'
    YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
    CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XSWU(:),IRESP,HCOMMENT=YCOMMENT)
    !
    YRECFM='LWD'
    YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
    CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XLWD(:),IRESP,HCOMMENT=YCOMMENT)
    !
    YRECFM='LWU'
    YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
    CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XLWU(:),IRESP,HCOMMENT=YCOMMENT)
    !
    DO JSW=1, SIZE(DGU%XSWBD,2)
      YNUM=ACHAR(48+JSW)
      !
      YRECFM='SWD_'//YNUM
      YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
      CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XSWBD(:,JSW),IRESP,HCOMMENT=YCOMMENT)
      !
      YRECFM='SWU_'//YNUM
      YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
      CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XSWBU(:,JSW),IRESP,HCOMMENT=YCOMMENT)
      !
    ENDDO
    !
  ENDIF
  !
  YRECFM='FMUNOSSO'
  YCOMMENT='X_Y_'//YRECFM//' (kg/ms2)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XFMU(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='FMVNOSSO'
  YCOMMENT='X_Y_'//YRECFM//' (kg/ms2)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XFMV(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='FMU'
  YCOMMENT='X_Y_'//YRECFM//' (kg/ms2)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XSSO_FMU(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='FMV'
  YCOMMENT='X_Y_'//YRECFM//' (kg/ms2)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XSSO_FMV(:),IRESP,HCOMMENT=YCOMMENT)
  !
END IF
!
! * Cumulated diag
!
IF (DGU%LSURF_BUDGETC) THEN
  !
  YRECFM='RNC'
  YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGUC%XRN(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='HC'
  YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGUC%XH(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='LEC'
  YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGUC%XLE(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='LEIC'
  YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGUC%XLEI(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='GFLUXC'
  YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGUC%XGFLUX(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='EVAPC'
  YCOMMENT='X_Y_'//YRECFM//' (kg/m2)'
  !
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGUC%XEVAP(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='SUBLC'
  YCOMMENT='X_Y_'//YRECFM//' (kg/m2)'
  !
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGUC%XSUBL(:),IRESP,HCOMMENT=YCOMMENT)
  !
  IF (DGU%LRAD_BUDGET .OR. (DGU%LSURF_BUDGETC .AND. .NOT.DGU%LRESET_BUDGETC)) THEN
    !        
    YRECFM='SWDC'
    YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
    CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGUC%XSWD(:),IRESP,HCOMMENT=YCOMMENT)
    !
    YRECFM='SWUC'
    YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
    CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGUC%XSWU(:),IRESP,HCOMMENT=YCOMMENT)
    !
    YRECFM='LWDC'
    YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
    CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGUC%XLWD(:),IRESP,HCOMMENT=YCOMMENT)
    !
    YRECFM='LWUC'
    YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
    CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGUC%XLWU(:),IRESP,HCOMMENT=YCOMMENT)
    !
  ENDIF
  !
  YRECFM='FMUC'
  YCOMMENT='X_Y_'//YRECFM//' (kg/ms)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGUC%XFMU(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='FMVC'
  YCOMMENT='X_Y_'//YRECFM//' (kg/ms)'
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGUC%XFMV(:),IRESP,HCOMMENT=YCOMMENT)
  !
END IF
!
!
!*       4.     Transfer coefficients
!               ---------------------
!
IF (DGU%LCOEF) THEN
  !
  YRECFM='CD'
  YCOMMENT='X_Y_'//YRECFM
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XCD(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='CH'
  YCOMMENT='X_Y_'//YRECFM
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XCH(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='CE'
  YCOMMENT='X_Y_'//YRECFM
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XCE(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='Z0'
  YCOMMENT='X_Y_'//YRECFM
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XZ0(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='Z0H'
  YCOMMENT='X_Y_'//YRECFM
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XZ0H(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='UREF'
  YCOMMENT='X_Y_'//YRECFM
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XUREF(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='ZREF'
  YCOMMENT='X_Y_'//YRECFM
  CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XZREF(:),IRESP,HCOMMENT=YCOMMENT)
  !
END IF
!
!
!*       5.     Surface humidity
!               ----------------
!
IF (DGU%LSURF_VARS) THEN
!
YRECFM='QS'
YCOMMENT='X_Y_'//YRECFM//' (kg/kg)'
!
 CALL WRITE_SURF(DGU, U, &
                  HPROGRAM,YRECFM,DGU%XQS(:),IRESP,HCOMMENT=YCOMMENT)
!
ENDIF
!
!-------------------------------------------------------------------------------
!
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_SEB_SURF_ATM_N',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE WRITE_DIAG_SEB_SURF_ATM_n
