!     #########
SUBROUTINE DIAG_INLINE_SEAFLUX_n(PTSTEP, PTA, PSST, PQA, &
     PPA, PPS, PRHOA, PZONA,                             &
     PMERA, PHT, PHW, PCD, PCDN, PCH, PCE, PRI, PHU,     &
     PZ0, PZ0H, PQSAT, PSFTH, PSFTQ, PSFZON, PSFMER,     &
     PDIR_SW, PSCA_SW, PLW, PDIR_ALB, PSCA_ALB, PICE_ALB,&
     PEMIS, PTRAD, PRAIN, PSNOW,                         & 
     TPGLT, PSIC, OHANDLE_SIC, PTICE,                    &
     PCD_ICE, PCDN_ICE, PCH_ICE, PCE_ICE, PRI_ICE,       &
     PZ0_ICE, PZ0H_ICE, PQSAT_ICE, PSFTH_ICE, PSFTQ_ICE, &
     PSFZON_ICE, PSFMER_ICE )
                                          
!     #####################################################################################
!
!!****  *DIAG_INLINE_SEAFLUX_n * - computes diagnostics during SEAFLUX time-step
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
!!      Modified    01/2006 : sea flux parameterization.
!!      B. Decharme 08/2009 : Diag for Earth System Model Coupling
!!      S. Riette   06/2009 CLS_2M becomes CLS_TQ, CLS_TQ and CLS_WIND have one
!!                          more argument (height of diagnostic)
!!      B. Decharme 04/2013 : Add EVAP and SUBL diag
!!      S. Senesi   01/2014 ! introduce fractional seaice and sea-ice model 
!!      S. Belamari 01/2014 : Wind module=XUNDEF if one component is XUNDEF
!!------------------------------------------------------------------
!

!
!
USE MODD_CSTS,           ONLY : XTTS
USE MODD_SURF_PAR,       ONLY : XUNDEF
USE MODD_SFX_OASIS,      ONLY : LCPL_SEA
USE MODD_SEAFLUX_n,      ONLY : LSBL, CSEAICE_SCHEME, LHANDLE_SIC
USE MODD_DIAG_SEAFLUX_n, ONLY : N2M, LSURF_BUDGET, LCOEF, LSURF_VARS,       &
                                       XT2M, XQ2M, XHU2M, XZON10M, XMER10M,   &
                                       XRN, XH, XLE, XGFLUX, XRI, XCD, XCH,   &
                                       XCE, XZ0, XZ0H, XQS, XSWD, XSWU, XLWD, &
                                       XLWU, XSWBD, XSWBU, XFMU, XFMV, XLE_ICE,  &
                                       LSURF_BUDGETC, XT2M_MIN, XT2M_MAX,     &
                                       XTS, XTSRAD, XHU2M_MIN, XHU2M_MAX,     &
                                       XWIND10M, XWIND10M_MAX, XEVAP, XSUBL,  &  
                                       XT2M_ICE, XQ2M_ICE, XHU2M_ICE,         &
                                       XZON10M_ICE, XMER10M_ICE, XWIND10M_ICE,&
                                       XRN_ICE, XH_ICE, XGFLUX_ICE, XRI_ICE,  &
                                       XCD_ICE, XCH_ICE,                      &
                                       XZ0_ICE, XZ0H_ICE, XQS_ICE, XSWU_ICE,  &
                                       XLWU_ICE, XSWBU_ICE, XFMU_ICE, XFMV_ICE
!
USE MODD_DIAG_SEAICE_n, ONLY : LDIAG_SEAICE, XSIT, XSND, XMLT
USE MODD_TYPES_GLT,     ONLY : T_GLT
USE MODD_GLT_PARAM ,    ONLY : GELATO_DIM=>NX
USE MODE_GLT_STATS ,    ONLY : GLT_AVHICEM, GLT_AVHSNWM
USE MODI_PARAM_CLS
USE MODI_CLS_TQ
USE MODI_CLS_WIND
USE MODI_DIAG_SURF_BUDGET_SEA
USE MODI_DIAG_SURF_BUDGETC_SEA
USE MODI_DIAG_CPL_ESM_SEA
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,               INTENT(IN) :: PTSTEP ! atmospheric time-step                 (s)
REAL, DIMENSION(:), INTENT(IN) :: PTA    ! atmospheric temperature
REAL, DIMENSION(:), INTENT(IN) :: PSST   ! sea surface temperature
REAL, DIMENSION(:), INTENT(IN) :: PQA    ! atmospheric specific humidity
REAL, DIMENSION(:), INTENT(IN) :: PPA    ! atmospheric level pressure
REAL, DIMENSION(:), INTENT(IN) :: PPS    ! surface pressure
REAL, DIMENSION(:), INTENT(IN) :: PRHOA  ! air density
REAL, DIMENSION(:), INTENT(IN) :: PZONA  ! zonal wind
REAL, DIMENSION(:), INTENT(IN) :: PMERA  ! meridian wind
REAL, DIMENSION(:), INTENT(IN) :: PHT    ! atmospheric level height
REAL, DIMENSION(:), INTENT(IN) :: PHW    ! atmospheric level height for wind
REAL, DIMENSION(:), INTENT(IN) :: PCD    ! drag coefficient for momentum
REAL, DIMENSION(:), INTENT(IN) :: PCDN   ! neutral drag coefficient
REAL, DIMENSION(:), INTENT(IN) :: PCH    ! drag coefficient for heat
REAL, DIMENSION(:), INTENT(IN) :: PCE    ! drag coefficient for vapor
REAL, DIMENSION(:), INTENT(IN) :: PRI    ! Richardson number
REAL, DIMENSION(:), INTENT(IN) :: PHU    ! near-surface humidity
REAL, DIMENSION(:), INTENT(IN) :: PZ0    ! roughness length for momentum
REAL, DIMENSION(:), INTENT(IN) :: PZ0H   ! roughness length for heat
REAL, DIMENSION(:), INTENT(IN) :: PQSAT  ! humidity at saturation
REAL, DIMENSION(:), INTENT(IN) :: PSFZON ! zonal friction
REAL, DIMENSION(:), INTENT(IN) :: PSFMER ! meridian friction
REAL, DIMENSION(:), INTENT(IN) :: PSFTH  ! heat flux  (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PSFTQ  ! water flux (kg/m2/s)
REAL, DIMENSION(:,:),INTENT(IN):: PDIR_SW   ! direct  solar radiation (on horizontal surf.)
!                                           !                                       (W/m2)
REAL, DIMENSION(:,:),INTENT(IN):: PSCA_SW   ! diffuse solar radiation (on horizontal surf.)
!                                           !                                       (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PLW       ! longwave radiation (on horizontal surf.)
REAL, DIMENSION(:), INTENT(IN) :: PTRAD     ! radiative temperature                 (K)
REAL, DIMENSION(:,:),INTENT(IN):: PDIR_ALB  ! direct albedo for each spectral band  (-)
REAL, DIMENSION(:,:),INTENT(IN):: PSCA_ALB  ! diffuse albedo for each spectral band (-)
REAL, DIMENSION(:)  ,INTENT(IN):: PICE_ALB  ! seaice albedo 
REAL, DIMENSION(:), INTENT(IN) :: PEMIS     ! emissivity                            (-)
!
REAL, DIMENSION(:), INTENT(IN) :: PRAIN     ! Rainfall (kg/m2/s)
REAL, DIMENSION(:), INTENT(IN) :: PSNOW     ! Snowfall (kg/m2/s)
!
TYPE(T_GLT)                    :: TPGLT     ! Sea-ice state , diagnostics and auxilliaries
REAL, DIMENSION(:), INTENT(IN) :: PSIC      ! Sea-ice cover
!
LOGICAL, INTENT(IN)               :: OHANDLE_SIC ! Do we weight seaice and open sea fluxes
REAL, DIMENSION(:), INTENT(IN)    :: PTICE      ! Seaice Surface Temperature
REAL, DIMENSION(:), INTENT(IN)    :: PCD_ICE    ! drag coefficient for momentum
REAL, DIMENSION(:), INTENT(IN)    :: PCDN_ICE   ! neutral drag coefficient
REAL, DIMENSION(:), INTENT(IN)    :: PCH_ICE    ! drag coefficient for heat
REAL, DIMENSION(:), INTENT(IN)    :: PCE_ICE    ! drag coefficient for vapor
REAL, DIMENSION(:), INTENT(IN)    :: PRI_ICE    ! Richardson number
REAL, DIMENSION(:), INTENT(IN)    :: PZ0_ICE    ! roughness length for momentum
REAL, DIMENSION(:), INTENT(IN)    :: PZ0H_ICE   ! roughness length for heat
REAL, DIMENSION(:), INTENT(IN)    :: PQSAT_ICE  ! humidity at saturation
REAL, DIMENSION(:), INTENT(IN)    :: PSFTH_ICE  ! heat flux  (W/m2)
REAL, DIMENSION(:), INTENT(IN)    :: PSFTQ_ICE  ! water flux (kg/m2/s)
REAL, DIMENSION(:), INTENT(IN)    :: PSFZON_ICE ! zonal friction
REAL, DIMENSION(:), INTENT(IN)    :: PSFMER_ICE ! meridian friction
!
!*      0.2    declarations of local variables
!
LOGICAL                         :: GSIC
REAL, DIMENSION(SIZE(PTA))      :: ZZ0W
REAL, DIMENSION(SIZE(PTA))      :: ZH
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DIAG_INLINE_SEAFLUX_N',0,ZHOOK_HANDLE)
!
! * Mean surface temperature need to couple with AGCM
!
IF (OHANDLE_SIC) THEN
   XTS   (:) = (1 - PSIC(:)) * PSST(:) + PSIC(:) * PTICE(:)
   XTSRAD(:) = PTRAD(:)
ELSE
   XTS   (:) = PSST (:)
   XTSRAD(:) = PTRAD(:)
ENDIF
!
IF (.NOT. LSBL) THEN
!
  IF (N2M==1) THEN        
    CALL PARAM_CLS(PTA, PSST, PQA, PPA, PRHOA, PZONA, PMERA, PHT, PHW, &
         PSFTH, PSFTQ, PSFZON, PSFMER,                                 &
         XT2M, XQ2M, XHU2M, XZON10M, XMER10M )  
    IF (OHANDLE_SIC) THEN
       CALL PARAM_CLS(PTA, PTICE, PQA, PPA, PRHOA, PZONA, PMERA, PHT, PHW, &
            PSFTH_ICE, PSFTQ_ICE, PSFZON_ICE, PSFMER_ICE,                  &
            XT2M_ICE, XQ2M_ICE, XHU2M_ICE, XZON10M_ICE, XMER10M_ICE  )  
    ENDIF
  ELSE IF (N2M==2) THEN
    ZH(:)=2.          
    CALL CLS_TQ(PTA, PQA, PPA, PPS, PHT,          &
                  PCD, PCH, PRI,                  &
                  PSST, PHU, PZ0H, ZH,            &
                  XT2M, XQ2M, XHU2M)
    ZH(:)=10.                
    CALL CLS_WIND(PZONA, PMERA, PHW,              &
                    PCD, PCDN, PRI, ZH,           &
                    XZON10M, XMER10M)  
    IF (OHANDLE_SIC) THEN
       ZH(:)=2.          
       CALL CLS_TQ(PTA, PQA, PPA, PPS, PHT,  &
            PCD_ICE, PCH_ICE, PRI_ICE,       &
            PTICE, PHU, PZ0H_ICE, ZH,        &
            XT2M_ICE, XQ2M_ICE, XHU2M_ICE)  
       ZH(:)=10.                
       CALL CLS_WIND(PZONA, PMERA, PHW,      &
            PCD_ICE, PCDN_ICE, PRI_ICE, ZH,  &
            XZON10M_ICE, XMER10M_ICE  )  
    ENDIF 
  END IF
!
  IF (N2M>=1) THEN
     IF (OHANDLE_SIC) THEN
        !
        XT2M    = XT2M    * (1 - PSIC) + XT2M_ICE    * PSIC
        XQ2M    = XQ2M    * (1 - PSIC) + XQ2M_ICE    * PSIC
        XHU2M   = XHU2M   * (1 - PSIC) + XHU2M_ICE   * PSIC
        !
        WHERE(XZON10M(:)/=XUNDEF.AND.XZON10M_ICE(:)/=XUNDEF)
              XZON10M(:) = XZON10M(:) * (1 - PSIC(:)) + XZON10M_ICE(:) * PSIC(:)
        ELSEWHERE
              XZON10M(:) = XUNDEF
        ENDWHERE
        WHERE(XMER10M(:)/=XUNDEF.AND.XMER10M_ICE(:)/=XUNDEF)
              XMER10M(:) = XMER10M(:) * (1 - PSIC(:)) + XMER10M_ICE(:) * PSIC(:)
        ELSEWHERE
              XMER10M(:) = XUNDEF
        ENDWHERE
        WHERE(XZON10M_ICE(:) /= XUNDEF .AND. XMER10M_ICE(:) /= XUNDEF)
              XWIND10M_ICE(:) = SQRT(XZON10M_ICE(:)**2+XMER10M_ICE(:)**2)
        ELSEWHERE
              XWIND10M_ICE(:) = XUNDEF
        ENDWHERE
        !
        XRI    = PRI     * (1 - PSIC) + PRI_ICE     * PSIC
        XRI_ICE=PRI_ICE
     ELSE
        XRI    =PRI
     ENDIF
    !
    XT2M_MIN(:) = MIN(XT2M_MIN(:),XT2M(:))
    XT2M_MAX(:) = MAX(XT2M_MAX(:),XT2M(:))
    !
    XHU2M_MIN(:) = MIN(XHU2M_MIN(:),XHU2M(:))
    XHU2M_MAX(:) = MAX(XHU2M_MAX(:),XHU2M(:))
    !
    WHERE(XZON10M(:) /= XUNDEF .AND. XMER10M(:) /= XUNDEF)
      XWIND10M(:) = SQRT(XZON10M(:)**2+XMER10M(:)**2)
    ELSEWHERE
      XWIND10M(:) = XUNDEF
    ENDWHERE
    XWIND10M_MAX(:) = MAX(XWIND10M_MAX(:),XWIND10M(:))
    !
  ENDIF
!
ELSE
  IF (N2M>=1) THEN
    XT2M    = XUNDEF
    XQ2M    = XUNDEF
    XHU2M   = XUNDEF
    XZON10M = XUNDEF
    XMER10M = XUNDEF
    XRI     = PRI
  ENDIF
ENDIF
!
IF (LSURF_BUDGET.OR.LSURF_BUDGETC) THEN
!
  CALL  DIAG_SURF_BUDGET_SEA   (XTTS, PSST, PRHOA, PSFTH, PSFTH_ICE,    &
                                  PSFTQ, PSFTQ_ICE,                     &
                                  PDIR_SW, PSCA_SW, PLW, PDIR_ALB,      &
                                  PSCA_ALB,PICE_ALB, PEMIS, PTRAD,      &
                                  PSFZON, PSFZON_ICE, PSFMER,           &
                                  PSFMER_ICE, OHANDLE_SIC, PSIC, PTICE, &
                                  XRN, XH, XLE, XLE_ICE, XGFLUX,        &
                                  XSWD, XSWU, XSWBD, XSWBU, XLWD, XLWU, &
                                  XFMU, XFMV, XEVAP, XSUBL,             &
                                  XRN_ICE, XH_ICE, XGFLUX_ICE,          &
                                  XSWU_ICE, XSWBU_ICE, XLWU_ICE,        &
                                  XFMU_ICE, XFMV_ICE                    )  
!
END IF
!
IF(LSURF_BUDGETC)THEN
  CALL DIAG_SURF_BUDGETC_SEA(PTSTEP, XRN, XH, XLE, XLE_ICE, XGFLUX,  &
                               XSWD, XSWU, XLWD, XLWU, XFMU, XFMV,   &
                               XEVAP, XSUBL, OHANDLE_SIC,            &
                               XRN_ICE, XH_ICE, XGFLUX_ICE,          &
                               XSWU_ICE, XLWU_ICE, XFMU_ICE, XFMV_ICE)
ENDIF
!
IF (LCOEF) THEN
   IF (OHANDLE_SIC) THEN 
      !
      !* Transfer coefficients
      !
      XCD = (1 - PSIC) * PCD + PSIC * PCD_ICE
      XCH = (1 - PSIC) * PCH + PSIC * PCH_ICE
      XCE = (1 - PSIC) * PCE + PSIC * PCE_ICE
      !
      !* Roughness lengths
      !
      ZZ0W = ( 1 - PSIC ) * 1.0/(LOG(PHW/PZ0)    **2)  +  &
                   PSIC   * 1.0/(LOG(PHW/PZ0_ICE)**2)  
      XZ0  = PHW  * EXP ( - SQRT ( 1./  ZZ0W ))
      ZZ0W = ( 1 - PSIC ) * 1.0/(LOG(PHW/PZ0H)    **2)  +  &
                   PSIC   * 1.0/(LOG(PHW/PZ0H_ICE)**2)  
      XZ0H = PHW  * EXP ( - SQRT ( 1./  ZZ0W ))
      !
   ELSE
      !
      !* Transfer coefficients
      !
      XCD = PCD
      XCH = PCH
      XCE = PCE
      !
      !* Roughness lengths
      !
      XZ0  = PZ0
      XZ0H = PZ0H
   ENDIF
   !
ENDIF
!
IF (LSURF_VARS) THEN
  !
  !* Humidity at saturation
  !
   IF (OHANDLE_SIC) THEN 
      XQS = (1 - PSIC) * PQSAT + PSIC * PQSAT_ICE
   ELSE 
      XQS = PQSAT
   ENDIF
ENDIF
!
! Diags from embedded Seaice model
! CALL DIAG_INLINE_SEAICE() : simply  : 
!
IF (LDIAG_SEAICE) THEN
   IF (TRIM(CSEAICE_SCHEME) == 'GELATO') THEN 
      GELATO_DIM=SIZE(PTA)
      XSIT  = RESHAPE(glt_avhicem(TPGLT%dom,TPGLT%sit),(/GELATO_DIM/))
      XSND  = RESHAPE(glt_avhsnwm(TPGLT%dom,TPGLT%sit),(/GELATO_DIM/))
      XMLT  = TPGLT%oce_all(:,1)%tml
   ELSE
      ! Placeholder for an alternate seaice scheme
   ENDIF
ENDIF
!
! Diags for Earth System Model coupling or for embedded Seaice model
! (we are actually using XCPL_.. variables for feeding the seaice model)
!
GSIC=(OHANDLE_SIC.AND.(CSEAICE_SCHEME /= 'NONE  '))
!
IF (LCPL_SEA.OR.GSIC) THEN
!
  CALL DIAG_CPL_ESM_SEA(PTSTEP,XZON10M,XMER10M,XFMU,XFMV,  &
                          XSWD,XSWU,XGFLUX,PSFTQ,PRAIN,    &
                          PSNOW,PLW,PTICE,PSFTH_ICE,       &
                          PSFTQ_ICE,PDIR_SW,PSCA_SW,       &
                          XSWU_ICE,XLWU_ICE,GSIC           )
! 
ENDIF
IF (LHOOK) CALL DR_HOOK('DIAG_INLINE_SEAFLUX_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_INLINE_SEAFLUX_n
