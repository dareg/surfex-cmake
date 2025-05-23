!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
!     #########
      SUBROUTINE WRITE_DIAG_SEB_FLAKE_n (DTCO, DUO, U, CHF, DFO, D, DC, HPROGRAM)
!     #################################
!
!!****  *WRITE_DIAG_SEB_FLAKE_n* - writes FLAKE diagnostics
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
!!      V. Masson   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      Modified    01/2006 : sea flux parameterization.
!!      P.LeMoigne    04/2013 : Add accumulated diagnostics
!!      Modified    04/2013, P. Le Moigne: FLake chemistry
!!      S. Belamari 06/2014 : Introduce NBLOCK to avoid errors due to NBLOCK=0
!!                            when coupled with ARPEGE/ALADIN/AROME
!!      B. Decharme 02/2016 : NBLOCK instead of LCOUNTW for compilation in AAA
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DIAG_n, ONLY : DIAG_t, DIAG_OPTIONS_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_CH_FLAKE_n, ONLY : CH_FLAKE_t
!
USE MODD_XIOS, ONLY : LALLOW_ADD_DIM, YSWBAND_DIM_NAME
!
USE MODD_SURF_PAR,      ONLY : XUNDEF
!
#ifdef SFX_OL
USE MODD_IO_SURF_OL, ONLY : LDEF
#endif
!
USE MODI_INIT_IO_SURF_n
USE MODI_WRITE_SURF
USE MODI_END_IO_SURF_n
USE MODD_SURFEX_HOST
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DUO
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(CH_FLAKE_t), INTENT(INOUT) :: CHF
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DFO
TYPE(DIAG_t), INTENT(INOUT) :: D
TYPE(DIAG_t), INTENT(INOUT) :: DC
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
LOGICAL           :: GRESET
INTEGER           :: JSV, JSW
REAL(KIND=JPHOOK)   :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!         Initialisation for IO
!
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_SEB_FLAKE_N',0,ZHOOK_HANDLE)
!
GRESET=.TRUE.

IF (ASSOCIATED (YRSURFEX_HOST)) THEN
  CALL YRSURFEX_HOST%SET_GRESET (DUO%LRESETMINMAX, GRESET)
ENDIF

#ifdef SFX_OL
IF (LDEF) GRESET = .FALSE.
#endif

!
 CALL INIT_IO_SURF_n(DTCO, U, HPROGRAM,'WATER ','FLAKE ','WRITE','FLAKE_DIAGNOSTICS.OUT.nc')
!
!
!*       2.     Richardson number :
!               -----------------
!
IF (DFO%N2M>=1) THEN

  YRECFM='RI_WAT'
  YCOMMENT='Bulk-Richardson number for water'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XRI(:),IRESP,HCOMMENT=YCOMMENT)
!
END IF
!
!*       3.     Energy fluxes :
!               -------------
!
IF (DFO%LSURF_BUDGET) THEN

  YRECFM='RN_WAT'
  YCOMMENT='net radiation for water'//' (W/m2)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XRN(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='H_WAT'
  YCOMMENT='sensible heat flux for water'//' (W/m2)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XH(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='LE_WAT'
  YCOMMENT='total latent heat flux for water'//' (W/m2)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XLE(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='LEI_WAT'
  YCOMMENT='sublimation latent heat flux for water-ice'//' (W/m2)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XLEI(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='GFLUX_WAT'
  YCOMMENT='conduction flux for water'//' (W/m2)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XGFLUX(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='EVAP_WAT'
  YCOMMENT='X_Y_'//YRECFM//' (kg/m2/s)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XEVAP(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='SUBL_WAT'
  YCOMMENT='X_Y_'//YRECFM//' (kg/m2/s)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XSUBL(:),IRESP,HCOMMENT=YCOMMENT)
  !
  IF (DFO%LRAD_BUDGET) THEN
    !
    YRECFM='SWD_WAT'
    YCOMMENT='short wave downward radiation for water'//' (W/m2)'
    CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XSWD(:),IRESP,HCOMMENT=YCOMMENT)
    !
    YRECFM='SWU_WAT'
    YCOMMENT='short wave upward radiation for water'//' (W/m2)'
    CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XSWU(:),IRESP,HCOMMENT=YCOMMENT)
    !
    YRECFM='LWD_WAT'
    YCOMMENT='downward long wave radiation'//' (W/m2)'
    CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XLWD(:),IRESP,HCOMMENT=YCOMMENT)
    !
    YRECFM='LWU_WAT'
    YCOMMENT='upward long wave radiation'//' (W/m2)'
    CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XLWU(:),IRESP,HCOMMENT=YCOMMENT)
    !  
    IF (LALLOW_ADD_DIM)  THEN
      !
      YRECFM='SWD_WAT'
      YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
      CALL WRITE_SURF(DUO%CSELECT,&
           HPROGRAM,YRECFM,D%XSWBD(:,:),IRESP,HCOMMENT=YCOMMENT, HNAM_DIM=YSWBAND_DIM_NAME)
      !
      YRECFM='SWU_WAT'
      YCOMMENT='X_Y_'//YRECFM//' (W/m2)'
      CALL WRITE_SURF(DUO%CSELECT,&
           HPROGRAM,YRECFM,D%XSWBD(:,:),IRESP,HCOMMENT=YCOMMENT, HNAM_DIM=YSWBAND_DIM_NAME)  
      !
    ELSE
      !    
      DO JSW=1, SIZE(D%XSWBD,2)
        YNUM=ACHAR(48+JSW)
        !
        YRECFM='SWD_WAT_'//YNUM
        YCOMMENT='downward short wave radiation by spectral band '//' (W/m2)'
        CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XSWBD(:,JSW),IRESP,HCOMMENT=YCOMMENT)
       !
         YRECFM='SWU_WAT_'//YNUM
        YCOMMENT='upward short wave radiation by spectral band'//' (W/m2)'
        CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XSWBU(:,JSW),IRESP,HCOMMENT=YCOMMENT)
        !
      ENDDO
      !
    ENDIF
    !
  ENDIF
  !
  YRECFM='FMU_WAT'
  YCOMMENT='u-component of momentum flux for water'//' (kg/ms2)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XFMU(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='FMV_WAT'
  YCOMMENT='v-component of momentum flux for water'//' (kg/ms2)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XFMV(:),IRESP,HCOMMENT=YCOMMENT)
  !
END IF
!
IF (DFO%LSURF_BUDGET.OR.DFO%LSURF_BUDGETC) THEN
!
  YRECFM='TALB_WAT'
  YCOMMENT='total albedo over tile water (-)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XALBT(:),IRESP,HCOMMENT=YCOMMENT)
!
  YRECFM='WSN_WAT'
  YCOMMENT='snow water equivalent over tile water (-)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XSWE(:),IRESP,HCOMMENT=YCOMMENT)
!        
ENDIF
!
!
!*       4.     Transfer coefficients
!               ---------------------
!
IF (DFO%LCOEF) THEN

  YRECFM='CD_WAT'
  YCOMMENT='drag coefficient for wind over water (W/s2)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XCD(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='CH_WAT'
  YCOMMENT='drag coefficient for heat (W/s)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XCH(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='CE_WAT'
  YCOMMENT='drag coefficient for vapor (W/s/K)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XCE(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='Z0_WAT'
  YCOMMENT='roughness length over water (m)'
  CALL WRITE_SURF(DUO%CSELECT, HPROGRAM,YRECFM,D%XZ0(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='Z0H_WAT'
  YCOMMENT='thermal roughness length over water (m)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XZ0H(:),IRESP,HCOMMENT=YCOMMENT)
  !
END IF
!
!
!*       5.     Surface humidity
!               ----------------
!
IF (DFO%LSURF_VARS) THEN

  YRECFM='QS_WAT'
  YCOMMENT='specific humidity over water'//' (KG/KG)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XQS(:),IRESP,HCOMMENT=YCOMMENT)
  !
ENDIF
!

!
!*       6.     parameters at 2 and 10 meters :
!               -----------------------------
!
IF (DFO%N2M>=1) THEN
  !
  YRECFM='T2M_WAT'
  YCOMMENT='2 meters temperature'//' (K)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XT2M(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='T2MMIN_WAT'
  YCOMMENT='X_Y_'//YRECFM//' (K)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XT2M_MIN(:),IRESP,HCOMMENT=YCOMMENT)
  IF(GRESET)D%XT2M_MIN(:)=XUNDEF
  !
  YRECFM='T2MMAX_WAT'
  YCOMMENT='X_Y_'//YRECFM//' (K)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XT2M_MAX(:),IRESP,HCOMMENT=YCOMMENT)
  IF(GRESET)D%XT2M_MAX(:)=0.0
  !
  YRECFM='Q2M_WAT'
  YCOMMENT='2 meters specific humidity'//' (KG/KG)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XQ2M(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='HU2M_WAT'
  YCOMMENT='2 meters relative humidity'//' (KG/KG)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XHU2M(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='HU2MMIN_WAT'
  YCOMMENT='X_Y_'//YRECFM//' (-)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XHU2M_MIN(:),IRESP,HCOMMENT=YCOMMENT)
  IF(GRESET)D%XHU2M_MIN(:)=XUNDEF
  !
  YRECFM='HU2MMAX_WAT'
  YCOMMENT='X_Y_'//YRECFM//' (-)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XHU2M_MAX(:),IRESP,HCOMMENT=YCOMMENT)
  IF(GRESET)D%XHU2M_MAX(:)=-XUNDEF
  !
  YRECFM='ZON10M_WAT'
  YCOMMENT='10 meters zonal wind'//' (M/S)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XZON10M(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='MER10M_WAT'
  YCOMMENT='10 meters meridian wind'//' (M/S)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XMER10M(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='W10M_WAT'
  YCOMMENT='X_Y_'//YRECFM//' (M/S)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XWIND10M(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='W10MMAX_WAT'
  YCOMMENT='X_Y_'//YRECFM//' (M/S)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XWIND10M_MAX(:),IRESP,HCOMMENT=YCOMMENT)
  IF(GRESET)D%XWIND10M_MAX(:)=0.0
  !
END IF
!
!
!*       7.     chemical diagnostics:
!               --------------------
!
IF (CHF%SVF%NBEQ>0 .AND. CHF%CCH_DRY_DEP=="WES89 ") THEN
  DO JSV = 1,SIZE(CHF%CCH_NAMES,1)
    YRECFM='DVWT'//TRIM(CHF%CCH_NAMES(JSV))
    WRITE(YCOMMENT,'(A13,I3.3)')'(m/s) DV_WAT_',JSV
    CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,CHF%XDEP(:,JSV),IRESP,HCOMMENT=YCOMMENT)
  END DO
ENDIF
!
!
!*       8.     prognostic variable diagnostics:
!               --------------------------------
!
IF(DUO%LPROVAR_TO_DIAG)THEN
!
  YRECFM='TS_WAT'
  YCOMMENT='TS_WATER (K)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,D%XTS(:),IRESP,HCOMMENT=YCOMMENT)
!
ENDIF
!
IF (DFO%LSURF_BUDGETC) THEN
  !
  CALL END_IO_SURF_n(HPROGRAM)
  CALL INIT_IO_SURF_n(DTCO, U, HPROGRAM,'WATER ','FLAKE ','WRITE','FLAKE_DIAGNOSTICS.OUT.nc')
  !
  YRECFM='RNC_WAT'
  YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,DC%XRN(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='HC_WAT'
  YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,DC%XH(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='LEC_WAT'
  YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,DC%XLE(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='LEIC_WAT'
  YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,DC%XLEI(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='GFLUXC_WAT'
  YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,DC%XGFLUX(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='EVAPC_WAT'
  YCOMMENT='X_Y_'//YRECFM//' (kg/m2)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,DC%XEVAP(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='SUBLC_WAT'
  YCOMMENT='X_Y_'//YRECFM//' (kg/m2)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,DC%XSUBL(:),IRESP,HCOMMENT=YCOMMENT)
  !
  IF (DFO%LRAD_BUDGET .OR. (DFO%LSURF_BUDGETC .AND. .NOT.DUO%LRESET_BUDGETC)) THEN
    !
    YRECFM='SWDC_WAT'
    YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
    CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,DC%XSWD(:),IRESP,HCOMMENT=YCOMMENT)
    !
    YRECFM='SWUC_WAT'
    YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
    CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,DC%XSWU(:),IRESP,HCOMMENT=YCOMMENT)
    !
    YRECFM='LWDC_WAT'
    YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
    CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,DC%XLWD(:),IRESP,HCOMMENT=YCOMMENT)
    !
    YRECFM='LWUC_WAT'
    YCOMMENT='X_Y_'//YRECFM//' (J/m2)'
    CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,DC%XLWU(:),IRESP,HCOMMENT=YCOMMENT)
    !
  ENDIF
  !
  YRECFM='FMUC_WAT'
  YCOMMENT='X_Y_'//YRECFM//' (kg/ms)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,DC%XFMU(:),IRESP,HCOMMENT=YCOMMENT)
  !
  YRECFM='FMVC_WAT'
  YCOMMENT='X_Y_'//YRECFM//' (kg/ms)'
  CALL WRITE_SURF(DUO%CSELECT,HPROGRAM,YRECFM,DC%XFMV(:),IRESP,HCOMMENT=YCOMMENT)
  !
END IF
!
!-------------------------------------------------------------------------------
!
!         End of IO
!
 CALL END_IO_SURF_n(HPROGRAM)
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_SEB_FLAKE_N',1,ZHOOK_HANDLE)
!
!
END SUBROUTINE WRITE_DIAG_SEB_FLAKE_n
