!     #########
       SUBROUTINE DIAG_CPL_ESM_SEA (S, DGS, DGSI, PTSTEP, PSFTQ, PRAIN, PSNOW, &
                                    PLW, PSFTH_ICE, PSFTQ_ICE, PDIR_SW, PSCA_SW, OSIC)  
!     ###################################################################
!
!!****  *DIAG_CPL_ESM_SEA * - Computes diagnostics over sea for 
!!                Earth system model coupling or embedded seaice scheme
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
!!     B. Decharme 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    08/2009
!!      S.Senesi    01/2014  Adapt to embedded seaice scheme (SWU and LWU 
!!                           for seaice are provided as inputs)
!!      A.Voldoire  04/2015  Add LCPL_SEAICE test
!!------------------------------------------------------------------
!
USE MODD_DIAG_n, ONLY : DIAG_t
USE MODD_SEAFLUX_n, ONLY : SEAFLUX_t
!
USE MODD_CSTS,      ONLY : XSTEFAN, XLSTT
USE MODD_WATER_PAR, ONLY : XEMISWATICE
!
USE MODD_SFX_OASIS, ONLY : LCPL_SEAICE
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(SEAFLUX_t), INTENT(INOUT) :: S
TYPE(DIAG_t), INTENT(INOUT) :: DGS
TYPE(DIAG_t), INTENT(INOUT) :: DGSI
!
REAL,               INTENT(IN) :: PTSTEP    ! atmospheric time-step
REAL, DIMENSION(:), INTENT(IN) :: PSFTQ     ! water flux
REAL, DIMENSION(:), INTENT(IN) :: PRAIN     ! Rainfall
REAL, DIMENSION(:), INTENT(IN) :: PSNOW     ! Snowfall
REAL, DIMENSION(:), INTENT(IN) :: PLW       ! longwave radiation (on horizontal surf.)
REAL, DIMENSION(:), INTENT(IN) :: PSFTH_ICE ! heat flux  (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PSFTQ_ICE ! water flux (kg/m2/s)
REAL, DIMENSION(:,:),INTENT(IN):: PDIR_SW   ! direct  solar radiation (on horizontal surf.)
REAL, DIMENSION(:,:),INTENT(IN):: PSCA_SW   ! diffuse solar radiation (on horizontal surf.)
LOGICAL,            INTENT(IN) :: OSIC
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(S%XICE_ALB)) :: ZSWU, ZTICE4
!
INTEGER                      :: ISWB ! number of SW bands
INTEGER                      :: JSWB ! loop counter on number of SW bands
INTEGER                      :: INI  ! number of points
INTEGER                      :: JI   ! loop counter on number of points
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('DIAG_CPL_ESM_SEA',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
! Total or free-ice sea flux
!-------------------------------------------------------------------------------------
!
!* 10m wind speed (m)
!
S%XCPL_SEA_WIND(:) = S%XCPL_SEA_WIND(:) + PTSTEP * SQRT(DGS%XZON10M(:)**2+DGS%XMER10M(:)**2)
! 
!* wind stress (Pa.s)
!
S%XCPL_SEA_FWSU(:) = S%XCPL_SEA_FWSU(:) + PTSTEP * DGS%XFMU(:)
S%XCPL_SEA_FWSV(:) = S%XCPL_SEA_FWSV(:) + PTSTEP * DGS%XFMV(:)
S%XCPL_SEA_FWSM(:) = S%XCPL_SEA_FWSM(:) + PTSTEP * SQRT(DGS%XFMU(:)**2+DGS%XFMV(:)**2)
!
!* Solar net heat flux (J/m2)
!
S%XCPL_SEA_SNET(:) = S%XCPL_SEA_SNET(:) + PTSTEP * (DGS%XSWD(:) - DGS%XSWU(:))
!
!* Non solar heat flux (J/m2)
!
S%XCPL_SEA_HEAT(:) = S%XCPL_SEA_HEAT(:) + PTSTEP * (DGS%XGFLUX(:) + DGS%XSWU(:) - DGS%XSWD(:)) 
!
!* Evaporation (kg/m2)
!
S%XCPL_SEA_EVAP(:) = S%XCPL_SEA_EVAP(:) + PTSTEP * PSFTQ(:)
!
!* Precip (kg/m2)
! 
S%XCPL_SEA_RAIN(:) = S%XCPL_SEA_RAIN(:) + PTSTEP * PRAIN(:)
S%XCPL_SEA_SNOW(:) = S%XCPL_SEA_SNOW(:) + PTSTEP * PSNOW(:)
!
!-------------------------------------------------------------------------------------
! Ice flux
!-------------------------------------------------------------------------------------
IF (LCPL_SEAICE.OR.OSIC) THEN
!
  INI  = SIZE(PDIR_SW,1)
  ISWB = SIZE(PDIR_SW,2)
!
!* Solar net heat flux (J/m2)
!
  IF (OSIC) THEN
    ZSWU(:)=DGSI%XSWU(:)
  ELSE
    ZSWU(:)=0.0
    DO JSWB=1,ISWB
      DO JI=1,INI
         ZSWU(JI) = ZSWU(JI) + (PDIR_SW(JI,JSWB)+PSCA_SW(JI,JSWB)) * S%XICE_ALB(JI)
      ENDDO
    ENDDO
  ENDIF
!
  S%XCPL_SEAICE_SNET(:) = S%XCPL_SEAICE_SNET(:) + PTSTEP * (DGS%XSWD(:) - ZSWU(:))
!
!* Non solar heat flux (J/m2)
!
  IF (OSIC) THEN
    S%XCPL_SEAICE_HEAT(:) = S%XCPL_SEAICE_HEAT(:) + PTSTEP * &
              ( PLW(:) - DGSI%XLWU(:) - PSFTH_ICE(:) - XLSTT*PSFTQ_ICE(:) )
  ELSE
    ZTICE4(:)=S%XTICE(:)**4
    S%XCPL_SEAICE_HEAT(:) = S%XCPL_SEAICE_HEAT(:) + PTSTEP * ( XEMISWATICE*(PLW(:)-XSTEFAN*ZTICE4(:)) &
                                                         - PSFTH_ICE(:) - XLSTT*PSFTQ_ICE(:)      ) 
  ENDIF 
!
!* Sublimation (kg/m2)
!
  S%XCPL_SEAICE_EVAP(:) = S%XCPL_SEAICE_EVAP(:) + PTSTEP * PSFTQ_ICE(:)
!
ENDIF
!
IF (LHOOK) CALL DR_HOOK('DIAG_CPL_ESM_SEA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_CPL_ESM_SEA
