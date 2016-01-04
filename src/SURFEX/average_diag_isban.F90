!     #########
      SUBROUTINE AVERAGE_DIAG_ISBA_n (DGEI, DGI, OCANOPY, PPATCH, PLE, &
                                      PHW,PHT,PSFCO2,PTRAD)
!     #######################################
!
!
!!****  *AVERAGE_DIAG_ISBA_n*  
!!
!!    PURPOSE
!!    -------
!      Average the diagnostics from all ISBA tiles
!     
!!**  METHOD
!!    ------
!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!      
!!    REFERENCE
!!    ---------
!!      
!!    AUTHOR
!!    ------
!!	S. Belair           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/03/95 
!!      V.Masson    20/03/96  remove abnormal averages and average TS**4 instead
!!                            of TS
!!      (J.Stein)   27/03/96 use only H and LE in the soil scheme
!!      A. Boone    27/11/02 revised to output ALMA variables, and general applications
!!      B. Decharme 17/08/09 cumulative radiatif budget
!!      V. Masson   10/2013  coherence between canopy and min/max T2M diagnostics
!!      B. Decharme    04/13 Averaged Trad already done in average_diag.F90
!!                           Good dimension for CO2 flux
!!      P. Samuelsson  10/13 Added min max for XT2M
!!      B. Decharme    02/15 No dependence on HW for 10M Wind diags
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!
!
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_t
USE MODD_DIAG_ISBA_n, ONLY : DIAG_ISBA_t
!
USE MODD_SURF_PAR,    ONLY : XUNDEF
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
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DGEI
TYPE(DIAG_ISBA_t), INTENT(INOUT) :: DGI
!
LOGICAL, INTENT(IN) :: OCANOPY
REAL, DIMENSION(:,:), INTENT(IN) :: PPATCH
REAL, DIMENSION(:,:), INTENT(IN) :: PLE
!
REAL, DIMENSION(:), INTENT(IN)       :: PHW    ! atmospheric level height for wind (m)
REAL, DIMENSION(:), INTENT(IN)       :: PHT    ! atmospheric level height (m)
REAL, DIMENSION(:), INTENT(IN)       :: PSFCO2 ! CO2 flux (m/s*kg_CO2/kg_air)
REAL, DIMENSION(:), INTENT(IN)       :: PTRAD  ! Radiative temperature (K)
!
!*      0.2    declarations of local variables
!
INTEGER                              :: JPATCH ! tile loop counter
INTEGER                              :: JSWB   ! band loop counter
REAL, DIMENSION(SIZE(PPATCH,1))      :: ZSUMPATCH
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!       0.     Initialization
!              --------------
!
IF (LHOOK) CALL DR_HOOK('AVERAGE_DIAG_ISBA_N',0,ZHOOK_HANDLE)
ZSUMPATCH(:) = 0.
DO JPATCH=1,SIZE(PPATCH,2)
  ZSUMPATCH(:) = ZSUMPATCH(:) + PPATCH(:,JPATCH)
END DO
!
!       1.     Energy fluxes
!              -------------
!
IF (DGI%LSURF_BUDGET) THEN
  DGI%XAVG_RN(:)     = 0.
  DGI%XAVG_H (:)     = 0.
  DGI%XAVG_LE(:)     = 0.
  DGI%XAVG_LEI(:)    = 0.
  DGI%XAVG_GFLUX(:)  = 0.
  DGI%XAVG_SWD(:)    = 0.
  DGI%XAVG_SWU(:)    = 0.
  DGI%XAVG_LWD(:)    = 0.
  DGI%XAVG_LWU(:)    = 0.
  DGI%XAVG_FMU(:)    = 0.
  DGI%XAVG_FMV(:)    = 0.
  DGI%XAVG_SWBD(:,:) = 0.
  DGI%XAVG_SWBU(:,:) = 0.
  !
  DO JPATCH=1,SIZE(PPATCH,2)
    WHERE (ZSUMPATCH(:) > 0.)
!
! Net radiation
!
      DGI%XAVG_RN(:)  = DGI%XAVG_RN(:) +PPATCH(:,JPATCH) * DGI%XRN(:,JPATCH)
!
! Sensible heat flux
!
      DGI%XAVG_H (:)  = DGI%XAVG_H (:) +PPATCH(:,JPATCH) * DGI%XH (:,JPATCH)
!
! Total latent heat flux
!
      DGI%XAVG_LE(:)  = DGI%XAVG_LE(:) +PPATCH(:,JPATCH) * PLE(:,JPATCH)
!
! Sublimation latent heat flux
!
      DGI%XAVG_LEI(:) = DGI%XAVG_LEI(:) +PPATCH(:,JPATCH) * DGI%XLEI(:,JPATCH)
!
! Storage flux
!
      DGI%XAVG_GFLUX(:)  = DGI%XAVG_GFLUX(:) +PPATCH(:,JPATCH) * DGI%XGFLUX(:,JPATCH)
!
! Downwards SW radiation
!
      DGI%XAVG_SWD(:)  = DGI%XAVG_SWD(:) +PPATCH(:,JPATCH) * DGI%XSWD(:,JPATCH)
!
! Upwards SW radiation
!
      DGI%XAVG_SWU(:)  = DGI%XAVG_SWU(:) +PPATCH(:,JPATCH) * DGI%XSWU(:,JPATCH)
!
! Downwards LW radiation
!
      DGI%XAVG_LWD(:)  = DGI%XAVG_LWD(:) +PPATCH(:,JPATCH) * DGI%XLWD(:,JPATCH)
!
! Upwards LW radiation
!
      DGI%XAVG_LWU(:)  = DGI%XAVG_LWU(:) +PPATCH(:,JPATCH) * DGI%XLWU(:,JPATCH)
!
! Zonal wind stress
!
      DGI%XAVG_FMU(:)  = DGI%XAVG_FMU(:) +PPATCH(:,JPATCH) * DGI%XFMU(:,JPATCH)
!
! Meridian wind stress
!
      DGI%XAVG_FMV(:)  = DGI%XAVG_FMV(:) +PPATCH(:,JPATCH) * DGI%XFMV(:,JPATCH)
!
    END WHERE
  END DO
!
  DO JPATCH=1,SIZE(PPATCH,2)
    DO JSWB=1,SIZE(DGI%XSWBD,2)
      WHERE (ZSUMPATCH(:) > 0.)
!
! Downwards SW radiation for each spectral band
!
        DGI%XAVG_SWBD(:,JSWB)  = DGI%XAVG_SWBD(:,JSWB) +PPATCH(:,JPATCH) * DGI%XSWBD(:,JSWB,JPATCH)
!
! Upwards SW radiation for each spectral band
!
        DGI%XAVG_SWBU(:,JSWB)  = DGI%XAVG_SWBU(:,JSWB) +PPATCH(:,JPATCH) * DGI%XSWBU(:,JSWB,JPATCH)
!
      END WHERE
    END DO
  END DO
END IF
!
IF (DGEI%LSURF_BUDGETC) THEN
   DGI%XAVG_SWDC(:) = 0.
   DGI%XAVG_SWUC(:) = 0.
   DGI%XAVG_LWDC(:) = 0.
   DGI%XAVG_LWUC(:) = 0.
   DGI%XAVG_FMUC(:) = 0.
   DGI%XAVG_FMVC(:) = 0.
   DO JPATCH=1,SIZE(PPATCH,2)
      WHERE (ZSUMPATCH(:) > 0.)
!
!        Downwards SW radiation
!
         DGI%XAVG_SWDC(:) = DGI%XAVG_SWDC(:) + PPATCH(:,JPATCH) * DGI%XSWDC(:,JPATCH)
!
!        Upwards SW radiation
!
         DGI%XAVG_SWUC(:) = DGI%XAVG_SWUC(:) + PPATCH(:,JPATCH) * DGI%XSWUC(:,JPATCH)
!
!        Downwards LW radiation
!
         DGI%XAVG_LWDC(:) = DGI%XAVG_LWDC(:) + PPATCH(:,JPATCH) * DGI%XLWDC(:,JPATCH)
!
!        Upwards LW radiation
!
         DGI%XAVG_LWUC(:) = DGI%XAVG_LWUC(:) + PPATCH(:,JPATCH) * DGI%XLWUC(:,JPATCH)
!
!        Zonal wind stress
!
         DGI%XAVG_FMUC(:) = DGI%XAVG_FMUC(:) + PPATCH(:,JPATCH) * DGI%XFMUC(:,JPATCH)
!
!        Meridian wind stress
!
         DGI%XAVG_FMVC(:) = DGI%XAVG_FMVC(:) + PPATCH(:,JPATCH) * DGI%XFMVC(:,JPATCH)
!
    END WHERE
  END DO
ENDIF    
!
!
!       2.     surface temperature and 2 meters parameters
!              -------------------------------------------
!
DGI%XAVG_TS(:) = 0.0
DO JPATCH=1,SIZE(PPATCH,2)
    WHERE (ZSUMPATCH(:) > 0.)
       DGI%XAVG_TS(:)  = DGI%XAVG_TS(:) + PPATCH(:,JPATCH) * DGI%XTS(:,JPATCH)
    END WHERE
END DO
!
IF (.NOT. OCANOPY .AND. DGI%N2M>=1) THEN

  DGI%XAVG_T2M(:)  = 0.
  DGI%XAVG_Q2M(:)  = 0.
  DGI%XAVG_HU2M(:)  = 0.
  !
  DO JPATCH=1,SIZE(PPATCH,2)
    WHERE (ZSUMPATCH(:) > 0.)
!
! 2 meters temperature
!
      DGI%XAVG_T2M(:)  = DGI%XAVG_T2M(:) + PPATCH(:,JPATCH) * DGI%XT2M(:,JPATCH)
!
! 2 meters humidity
!
      DGI%XAVG_Q2M(:)  = DGI%XAVG_Q2M(:) + PPATCH(:,JPATCH) * DGI%XQ2M(:,JPATCH)
!
! 2 meters relative humidity
!
      DGI%XAVG_HU2M(:)  = DGI%XAVG_HU2M(:) + PPATCH(:,JPATCH) * DGI%XHU2M(:,JPATCH)
!
    END WHERE
  END DO
!
! 10 meters wind
!
  DGI%XAVG_ZON10M (:)  = 0.
  DGI%XAVG_MER10M (:)  = 0.
  DGI%XAVG_WIND10M(:)  = 0.
  DO JPATCH=1,SIZE(PPATCH,2)
    WHERE (ZSUMPATCH(:) > 0.)
      DGI%XAVG_ZON10M(:)  = DGI%XAVG_ZON10M (:) + PPATCH(:,JPATCH) * DGI%XZON10M (:,JPATCH)
      DGI%XAVG_MER10M(:)  = DGI%XAVG_MER10M (:) + PPATCH(:,JPATCH) * DGI%XMER10M (:,JPATCH)
      DGI%XAVG_WIND10M(:) = DGI%XAVG_WIND10M(:) + PPATCH(:,JPATCH) * DGI%XWIND10M(:,JPATCH)
    END WHERE
  ENDDO
!
  DGI%XAVG_T2M_MIN(:) = MIN(DGI%XAVG_T2M_MIN(:),DGI%XAVG_T2M(:))
  DGI%XAVG_T2M_MAX(:) = MAX(DGI%XAVG_T2M_MAX(:),DGI%XAVG_T2M(:))
!
  DGI%XAVG_HU2M_MIN(:) = MIN(DGI%XAVG_HU2M_MIN(:),DGI%XAVG_HU2M(:))
  DGI%XAVG_HU2M_MAX(:) = MAX(DGI%XAVG_HU2M_MAX(:),DGI%XAVG_HU2M(:))
!
  DGI%XAVG_WIND10M_MAX(:) = MAX(DGI%XAVG_WIND10M_MAX(:),DGI%XAVG_WIND10M(:))
!
END IF
!
! Richardson number
!
IF (DGI%N2M>=1) THEN

  DGI%XAVG_RI(:)  = 0.
  !
  DGI%XAVG_SFCO2(:)  = PSFCO2(:)
  !
  DO JPATCH=1,SIZE(PPATCH,2)
    WHERE (ZSUMPATCH(:) > 0.)
      DGI%XAVG_RI(:)  = DGI%XAVG_RI(:) + PPATCH(:,JPATCH) * DGI%XRI(:,JPATCH)
    END WHERE
  END DO
!
! min and max of XT2M
!
  DGI%XT2M_MIN(:,:) = MIN(DGI%XT2M_MIN(:,:),DGI%XT2M(:,:))
  DGI%XT2M_MAX(:,:) = MAX(DGI%XT2M_MAX(:,:),DGI%XT2M(:,:))
!
END IF
!
!       3.     Transfer coefficients
!              ---------------------
!
IF (DGI%LCOEF) THEN
  !
  DGI%XAVG_CD   (:) = 0.
  DGI%XAVG_CH   (:) = 0.
  DGI%XAVG_CE   (:) = 0.
  DGI%XAVG_Z0   (:) = 0.
  DGI%XAVG_Z0H  (:) = 0.
  DGI%XAVG_Z0EFF(:) = 0.
  !
  DO JPATCH=1,SIZE(PPATCH,2)
    WHERE (ZSUMPATCH(:) > 0.)
      !
      DGI%XAVG_CD(:)  = DGI%XAVG_CD(:) + PPATCH(:,JPATCH) * DGI%XCD(:,JPATCH)
      !
      DGI%XAVG_CH(:)  = DGI%XAVG_CH(:) + PPATCH(:,JPATCH) * DGI%XCH(:,JPATCH)
      !
      DGI%XAVG_CE(:)  = DGI%XAVG_CE(:) + PPATCH(:,JPATCH) * DGI%XCE(:,JPATCH)
      !
      !             
      DGI%XAVG_Z0(:)    = DGI%XAVG_Z0(:)    + PPATCH(:,JPATCH) * 1./(LOG(PHW(:)/DGI%XZ0_WITH_SNOW (:,JPATCH)))**2
      !      
      DGI%XAVG_Z0H(:)   = DGI%XAVG_Z0H(:)   + PPATCH(:,JPATCH) * 1./(LOG(PHT(:)/DGI%XZ0H_WITH_SNOW(:,JPATCH)))**2
      !      
      DGI%XAVG_Z0EFF(:) = DGI%XAVG_Z0EFF(:) + PPATCH(:,JPATCH) * 1./(LOG(PHW(:)/DGI%XZ0EFF        (:,JPATCH)))**2
      !      
    END WHERE
  END DO
  !
  DGI%XAVG_Z0(:)    = PHW(:) *  EXP( - SQRT(1./DGI%XAVG_Z0(:)) )
  !
  DGI%XAVG_Z0H(:)   = PHT(:) *  EXP( - SQRT(1./DGI%XAVG_Z0H(:)) )
  !
  DGI%XAVG_Z0EFF(:) = PHW(:) *  EXP( - SQRT(1./DGI%XAVG_Z0EFF(:)) )
  !
END IF
!
IF (DGI%LSURF_VARS) THEN
  DGI%XAVG_QS(:)  = 0.
  !
  DO JPATCH=1,SIZE(PPATCH,2)
    WHERE (ZSUMPATCH(:) > 0.)
!
! specific humidity at surface
!
      DGI%XAVG_QS(:)  = DGI%XAVG_QS(:) + PPATCH(:,JPATCH) * DGI%XQS(:,JPATCH)
!
    END WHERE
  END DO
END IF
!
IF (LHOOK) CALL DR_HOOK('AVERAGE_DIAG_ISBA_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE AVERAGE_DIAG_ISBA_n
