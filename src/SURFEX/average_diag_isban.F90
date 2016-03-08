!     #########
      SUBROUTINE AVERAGE_DIAG_ISBA_n (DIO, DGI, DGIC, DGIP, DGIPC, OSURF_BUDGETC, &
                                      OCANOPY, PPATCH, PLE, PHW, PHT ,PSFCO2, PTRAD)
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
USE MODD_DIAG_n, ONLY : DIAG_t, DIAG_PATCH_t, DIAG_OPTIONS_t
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
TYPE(DIAG_OPTIONS_t), INTENT(INOUT) :: DIO
TYPE(DIAG_t), INTENT(INOUT) :: DGI
TYPE(DIAG_t), INTENT(INOUT) :: DGIC
TYPE(DIAG_PATCH_t), INTENT(INOUT) :: DGIP
TYPE(DIAG_PATCH_t), INTENT(INOUT) :: DGIPC
!
LOGICAL, INTENT(IN) :: OSURF_BUDGETC
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
!
ZSUMPATCH(:) = 0.
DO JPATCH=1,SIZE(PPATCH,2)
  ZSUMPATCH(:) = ZSUMPATCH(:) + PPATCH(:,JPATCH)
END DO
!
!       1.     Energy fluxes
!              -------------
!
IF (DIO%LSURF_BUDGET) THEN
  DGI%XRN(:)     = 0.
  DGI%XH (:)     = 0.
  DGI%XLE(:)     = 0.
  DGI%XLEI(:)    = 0.
  DGI%XGFLUX(:)  = 0.
  DGI%XEVAP       (:) = 0.
  DGI%XSUBL       (:) = 0.  
  DGI%XSWD(:)    = 0.
  DGI%XSWU(:)    = 0.
  DGI%XLWD(:)    = 0.
  DGI%XLWU(:)    = 0.
  DGI%XFMU(:)    = 0.
  DGI%XFMV(:)    = 0.
  DGI%XSWBD(:,:) = 0.
  DGI%XSWBU(:,:) = 0.
  !
  DO JPATCH=1,SIZE(PPATCH,2)
    WHERE (ZSUMPATCH(:) > 0.)
!
! Net radiation
!
      DGI%XRN(:)  = DGI%XRN(:) +PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XRN(:)
!
! Sensible heat flux
!
      DGI%XH (:)  = DGI%XH (:) +PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XH (:)
!
! Total latent heat flux
!
      DGI%XLE(:)  = DGI%XLE(:) +PPATCH(:,JPATCH) * PLE(:,JPATCH)
!
! Sublimation latent heat flux
!
      DGI%XLEI(:) = DGI%XLEI(:) +PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XLEI(:)
!
! Storage flux
!
      DGI%XGFLUX(:)  = DGI%XGFLUX(:) +PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XGFLUX(:)
!
! Evapotranspiration
!
      DGI%XEVAP(:)  = DGI%XEVAP(:) + PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XEVAP(:)
!
! Sublimation
!
      DGI%XSUBL(:)  = DGI%XSUBL(:) + PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XSUBL(:)
      
!
! Downwards SW radiation
!
      DGI%XSWD(:)  = DGI%XSWD(:) +PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XSWD(:)
!
! Upwards SW radiation
!
      DGI%XSWU(:)  = DGI%XSWU(:) +PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XSWU(:)
!
! Downwards LW radiation
!
      DGI%XLWD(:)  = DGI%XLWD(:) +PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XLWD(:)
!
! Upwards LW radiation
!
      DGI%XLWU(:)  = DGI%XLWU(:) +PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XLWU(:)
!
! Zonal wind stress
!
      DGI%XFMU(:)  = DGI%XFMU(:) +PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XFMU(:)
!
! Meridian wind stress
!
      DGI%XFMV(:)  = DGI%XFMV(:) +PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XFMV(:)
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
        DGI%XSWBD(:,JSWB)  = DGI%XSWBD(:,JSWB) +PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XSWBD(:,JSWB)
!
! Upwards SW radiation for each spectral band
!
        DGI%XSWBU(:,JSWB)  = DGI%XSWBU(:,JSWB) +PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XSWBU(:,JSWB)
!
      END WHERE
    END DO
  END DO
END IF
!
IF (OSURF_BUDGETC) THEN
   DGIC%XRN        (:) = 0.
   DGIC%XH         (:) = 0.
   DGIC%XLE        (:) = 0.
   DGIC%XGFLUX     (:) = 0.
   DGIC%XLEI       (:) = 0.    
   DGIC%XEVAP      (:) = 0.
   DGIC%XSUBL      (:) = 0.
   DGIC%XSWD(:) = 0.
   DGIC%XSWU(:) = 0.
   DGIC%XLWD(:) = 0.
   DGIC%XLWU(:) = 0.
   DGIC%XFMU(:) = 0.
   DGIC%XFMV(:) = 0.
   DO JPATCH=1,SIZE(PPATCH,2)
      WHERE (ZSUMPATCH(:) > 0.)
!
! Net radiation
!
        DGIC%XRN(:)  = DGIC%XRN(:) + PPATCH(:,JPATCH) * DGIPC%AL(JPATCH)%XRN(:)
!
! Sensible heat flux
!
        DGIC%XH(:)  = DGIC%XH(:) + PPATCH(:,JPATCH) * DGIPC%AL(JPATCH)%XH(:)
!
! Total latent heat flux
!
        DGIC%XLE(:)  = DGIC%XLE(:) + PPATCH(:,JPATCH) * DGIPC%AL(JPATCH)%XLE(:)
!
! Storage flux
!
        DGIC%XGFLUX(:)  = DGIC%XGFLUX(:) + PPATCH(:,JPATCH) * DGIPC%AL(JPATCH)%XGFLUX(:)
!
! Total surface sublimation
!
        DGIC%XLEI(:)  = DGIC%XLEI(:) + PPATCH(:,JPATCH) * DGIPC%AL(JPATCH)%XLEI(:)              
!
! Evapotranspiration
!
        DGIC%XEVAP(:)  = DGIC%XEVAP(:) + PPATCH(:,JPATCH) * DGIPC%AL(JPATCH)%XEVAP(:)
!
! Sublimation
!
        DGIC%XSUBL(:)  = DGIC%XSUBL(:) + PPATCH(:,JPATCH) * DGIPC%AL(JPATCH)%XSUBL(:)
!
!        Downwards SW radiation
!
         DGIC%XSWD(:) = DGIC%XSWD(:) + PPATCH(:,JPATCH) * DGIPC%AL(JPATCH)%XSWD(:)
!
!        Upwards SW radiation
!
         DGIC%XSWU(:) = DGIC%XSWU(:) + PPATCH(:,JPATCH) * DGIPC%AL(JPATCH)%XSWU(:)
!
!        Downwards LW radiation
!
         DGIC%XLWD(:) = DGIC%XLWD(:) + PPATCH(:,JPATCH) * DGIPC%AL(JPATCH)%XLWD(:)
!
!        Upwards LW radiation
!
         DGIC%XLWU(:) = DGIC%XLWU(:) + PPATCH(:,JPATCH) * DGIPC%AL(JPATCH)%XLWU(:)
!
!        Zonal wind stress
!
         DGIC%XFMU(:) = DGIC%XFMU(:) + PPATCH(:,JPATCH) * DGIPC%AL(JPATCH)%XFMU(:)
!
!        Meridian wind stress
!
         DGIC%XFMV(:) = DGIC%XFMV(:) + PPATCH(:,JPATCH) * DGIPC%AL(JPATCH)%XFMV(:)
!
    END WHERE
  END DO
ENDIF    
!
!
!       2.     surface temperature and 2 meters parameters
!              -------------------------------------------
!
DGI%XTS(:) = 0.0
DO JPATCH=1,SIZE(PPATCH,2)
    WHERE (ZSUMPATCH(:) > 0.)
       DGI%XTS(:)  = DGI%XTS(:) + PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XTS(:)
    END WHERE
END DO
!
DGI%XALBT(:)   = 0.
DO JPATCH=1,SIZE(PPATCH,2)
  WHERE (ZSUMPATCH(:) > 0.)
!   Total albedo
    DGI%XALBT(:) = DGI%XALBT(:) + PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XALBT (:)       
  END WHERE
END DO
!
IF (.NOT. OCANOPY .AND. DIO%N2M>=1) THEN

  DGI%XT2M(:)  = 0.
  DGI%XQ2M(:)  = 0.
  DGI%XHU2M(:)  = 0.
  !
  DO JPATCH=1,SIZE(PPATCH,2)
    WHERE (ZSUMPATCH(:) > 0.)
!
! 2 meters temperature
!
      DGI%XT2M(:)  = DGI%XT2M(:) + PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XT2M(:)
!
! 2 meters humidity
!
      DGI%XQ2M(:)  = DGI%XQ2M(:) + PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XQ2M(:)
!
! 2 meters relative humidity
!
      DGI%XHU2M(:)  = DGI%XHU2M(:) + PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XHU2M(:)
!
    END WHERE
  END DO
!
! 10 meters wind
!
  DGI%XZON10M (:)  = 0.
  DGI%XMER10M (:)  = 0.
  DGI%XWIND10M(:)  = 0.
  DO JPATCH=1,SIZE(PPATCH,2)
    WHERE (ZSUMPATCH(:) > 0.)
      DGI%XZON10M(:)  = DGI%XZON10M (:) + PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XZON10M (:)
      DGI%XMER10M(:)  = DGI%XMER10M (:) + PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XMER10M (:)
      DGI%XWIND10M(:) = DGI%XWIND10M(:) + PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XWIND10M(:)
    END WHERE
  ENDDO
!
  DGI%XT2M_MIN(:) = MIN(DGI%XT2M_MIN(:),DGI%XT2M(:))
  DGI%XT2M_MAX(:) = MAX(DGI%XT2M_MAX(:),DGI%XT2M(:))
!
  DGI%XHU2M_MIN(:) = MIN(DGI%XHU2M_MIN(:),DGI%XHU2M(:))
  DGI%XHU2M_MAX(:) = MAX(DGI%XHU2M_MAX(:),DGI%XHU2M(:))
!
  DGI%XWIND10M_MAX(:) = MAX(DGI%XWIND10M_MAX(:),DGI%XWIND10M(:))
!
END IF
!
! Richardson number
!
IF (DIO%N2M>=1) THEN

  DGI%XRI(:)  = 0.
  !
  DGI%XSFCO2(:)  = PSFCO2(:)
  !
  DO JPATCH=1,SIZE(PPATCH,2)
    WHERE (ZSUMPATCH(:) > 0.)
      DGI%XRI(:)  = DGI%XRI(:) + PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XRI(:)
    END WHERE
  END DO
!
! min and max of XT2M
!
  DO JPATCH=1,SIZE(PPATCH,2)
    DGIP%AL(JPATCH)%XT2M_MIN(:) = MIN(DGIP%AL(JPATCH)%XT2M_MIN(:),DGIP%AL(JPATCH)%XT2M(:))
    DGIP%AL(JPATCH)%XT2M_MAX(:) = MAX(DGIP%AL(JPATCH)%XT2M_MAX(:),DGIP%AL(JPATCH)%XT2M(:))
  ENDDO
!
END IF
!
!       3.     Transfer coefficients
!              ---------------------
!
IF (DIO%LCOEF) THEN
  !
  DGI%XCD   (:) = 0.
  DGI%XCH   (:) = 0.
  DGI%XCE   (:) = 0.
  DGI%XZ0   (:) = 0.
  DGI%XZ0H  (:) = 0.
  DGI%XZ0EFF(:) = 0.
  !
  DO JPATCH=1,SIZE(PPATCH,2)
    WHERE (ZSUMPATCH(:) > 0.)
      !
      DGI%XCD(:)  = DGI%XCD(:) + PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XCD(:)
      !
      DGI%XCH(:)  = DGI%XCH(:) + PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XCH(:)
      !
      DGI%XCE(:)  = DGI%XCE(:) + PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XCE(:)
      !
      !             
      DGI%XZ0(:)    = DGI%XZ0(:)    + PPATCH(:,JPATCH) * 1./(LOG(PHW(:)/DGIP%AL(JPATCH)%XZ0 (:)))**2
      !      
      DGI%XZ0H(:)   = DGI%XZ0H(:)   + PPATCH(:,JPATCH) * 1./(LOG(PHT(:)/DGIP%AL(JPATCH)%XZ0H(:)))**2
      !      
      DGI%XZ0EFF(:) = DGI%XZ0EFF(:) + PPATCH(:,JPATCH) * 1./(LOG(PHW(:)/DGIP%AL(JPATCH)%XZ0EFF(:)))**2
      !      
    END WHERE
  END DO
  !
  DGI%XZ0(:)    = PHW(:) *  EXP( - SQRT(1./DGI%XZ0(:)) )
  !
  DGI%XZ0H(:)   = PHT(:) *  EXP( - SQRT(1./DGI%XZ0H(:)) )
  !
  DGI%XZ0EFF(:) = PHW(:) *  EXP( - SQRT(1./DGI%XZ0EFF(:)) )
  !
END IF
!
IF (DIO%LSURF_VARS) THEN
  DGI%XQS(:)  = 0.
  !
  DO JPATCH=1,SIZE(PPATCH,2)
    WHERE (ZSUMPATCH(:) > 0.)
!
! specific humidity at surface
!
      DGI%XQS(:)  = DGI%XQS(:) + PPATCH(:,JPATCH) * DGIP%AL(JPATCH)%XQS(:)
!
    END WHERE
  END DO
END IF
!
IF (LHOOK) CALL DR_HOOK('AVERAGE_DIAG_ISBA_N',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------
!
END SUBROUTINE AVERAGE_DIAG_ISBA_n
