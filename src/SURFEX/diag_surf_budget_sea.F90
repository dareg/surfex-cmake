!     #########
SUBROUTINE  DIAG_SURF_BUDGET_SEA(DGS, DGSI, S, PTT, PRHOA, PSFTH, PSFTH_ICE, &
                                 PSFTQ, PSFTQ_ICE, PDIR_SW, PSCA_SW, PLW,    &
                                 PDIR_ALB, PSCA_ALB, PEMIS, PTRAD,           &
                                 PSFZON, PSFZON_ICE, PSFMER, PSFMER_ICE   ) 


!     ###############################################################################
!
!!****  *DIAG_SURF_BUDGET_WATER * - Computes diagnostics over water
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
!       B. decharme 04/2013 : Add EVAP and SUBL diag
!       S.Senesi    01/2014 : Handle fluxes on seaice
!!------------------------------------------------------------------
!
USE MODD_SEAFLUX_n, ONLY : SEAFLUX_t
USE MODD_DIAG_n, ONLY : DIAG_t
!
USE MODD_CSTS,           ONLY : XSTEFAN, XLSTT, XLVTT
USE MODD_WATER_PAR,      ONLY : XEMISWATICE, XALBSEAICE
! 
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(DIAG_t), INTENT(INOUT) :: DGS
TYPE(DIAG_t), INTENT(INOUT) :: DGSI
TYPE(SEAFLUX_t), INTENT(INOUT) :: S
!
REAL,               INTENT(IN) :: PTT       ! freezing temperature of water surface
REAL, DIMENSION(:), INTENT(IN) :: PRHOA     ! air density
REAL, DIMENSION(:), INTENT(IN) :: PSFTH     ! heat flux
REAL, DIMENSION(:), INTENT(IN) :: PSFTH_ICE ! heat flux on seaice
!
REAL, DIMENSION(:), INTENT(IN) :: PSFTQ     ! water flux
REAL, DIMENSION(:), INTENT(IN) :: PSFTQ_ICE ! water flux on seaice
!
REAL, DIMENSION(:,:),INTENT(IN):: PDIR_SW   ! direct  solar radiation (on horizontal surf.)
!                                           !                                       (W/m2)
REAL, DIMENSION(:,:),INTENT(IN):: PSCA_SW   ! diffuse solar radiation (on horizontal surf.)
!                                           !                                       (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PLW       ! longwave radiation (on horizontal surf.)
!
REAL, DIMENSION(:,:),INTENT(IN):: PDIR_ALB  ! direct albedo for each spectral band  (-)
REAL, DIMENSION(:,:),INTENT(IN):: PSCA_ALB  ! diffuse albedo for each spectral band (-)
REAL, DIMENSION(:), INTENT(IN) :: PEMIS     ! emissivity                            (-)
REAL, DIMENSION(:), INTENT(IN) :: PTRAD     ! radiative temperature                 (K)
!
REAL, DIMENSION(:), INTENT(IN) :: PSFZON    ! zonal friction
REAL, DIMENSION(:), INTENT(IN) :: PSFZON_ICE! zonal friction
REAL, DIMENSION(:), INTENT(IN) :: PSFMER    ! meridional friction
REAL, DIMENSION(:), INTENT(IN) :: PSFMER_ICE! meridional friction
!
!*      0.2    declarations of local variables
!
INTEGER                      :: I
INTEGER                      :: ISWB ! number of SW bands
INTEGER                      :: JSWB ! loop counter on number of SW bands
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_BUDGET_SEA',0,ZHOOK_HANDLE)
!
ISWB = SIZE(PDIR_SW,2)
! 
!* total incoming and outgoing SW
!
DO JSWB=1,ISWB
  DGS%XSWBD(:,JSWB) = PDIR_SW(:,JSWB)                    + PSCA_SW(:,JSWB)
  DGS%XSWBU(:,JSWB) = PDIR_SW(:,JSWB) * PDIR_ALB(:,JSWB) + PSCA_SW(:,JSWB) * PSCA_ALB(:,JSWB) 
ENDDO
!
DGS%XSWD(:) = 0.
DGS%XSWU(:) = 0.
DO JSWB=1,ISWB
   DGS%XSWD(:) = DGS%XSWD(:) + DGS%XSWBD(:,JSWB)
   DGS%XSWU(:) = DGS%XSWU(:) + DGS%XSWBU(:,JSWB)
ENDDO
!
!*incoming outgoing LW
!
DGS%XLWD(:)=PLW(:)
DGS%XLWU(:)=PEMIS(:)*XSTEFAN*PTRAD(:)**4 + (1.-PEMIS(:))*PLW(:)
!
!* net radiation
!
DGS%XRN(:)    =   DGS%XSWD(:) - DGS%XSWU(:)     + DGS%XLWD(:) - DGS%XLWU    (:)
!
IF (.NOT.S%LHANDLE_SIC) THEN
  !
  !* sensible heat flux
  !
  DGS%XH     = PSFTH
  !
  !* latent heat flux
  !
  WHERE (S%XSST<PTT  )
     DGS%XLE    = PSFTQ * XLSTT
     DGS%XLEI   = PSFTQ * XLSTT
     DGS%XEVAP  = PSFTQ
     DGS%XSUBL  = PSFTQ
  ELSEWHERE
     DGS%XLE    = PSFTQ * XLVTT
     DGS%XLEI   = 0.0
     DGS%XEVAP  = PSFTQ
     DGS%XSUBL  = 0.0
  END WHERE
  !
  !* wind stress
  !
  DGS%XFMU = PSFZON
  DGS%XFMV = PSFMER
  !
ELSE
  !
  !---------------------------------------------------------------------------- 
  ! Sea ice or mixed diag
  !---------------------------------------------------------------------------- 
  !
  !
  !* total incoming and outgoing SW
  !
  DO JSWB=1,ISWB
    DGSI%XSWBU(:,JSWB) = (PDIR_SW(:,JSWB) + PSCA_SW(:,JSWB)) * S%XICE_ALB(:) 
  ENDDO
  !
  DGSI%XSWU(:) = 0.
  DO JSWB=1,ISWB
     DGSI%XSWU(:) = DGSI%XSWU(:) + DGSI%XSWBU(:,JSWB)
  ENDDO
  !
  !*incoming outgoing LW
  !
  DGSI%XLWU(:)=XEMISWATICE*XSTEFAN*S%XTICE(:)**4 + (1.-XEMISWATICE)*PLW(:)
  !
  !* net radiation
  !
  DGSI%XRN(:) =   DGS%XSWD(:) - DGSI%XSWU(:) + DGS%XLWD(:) - DGSI%XLWU(:)
  !
  !* sensible heat flux
  !
  DGS%XH     = (1 - S%XSIC) * PSFTH + S%XSIC * PSFTH_ICE 
  DGSI%XH    =                                 PSFTH_ICE
  !
  !* latent heat flux
  !
  DGS%XLE     = (1 - S%XSIC) * PSFTQ * XLVTT + S%XSIC * PSFTQ_ICE * XLSTT
  DGS%XLEI    =                                         PSFTQ_ICE * XLSTT
  DGS%XEVAP   = (1 - S%XSIC) * PSFTQ         + S%XSIC * PSFTQ_ICE 
  DGS%XSUBL   =                                S%XSIC * PSFTQ_ICE 
  !
  !* ice storage flux
  !
  DGSI%XGFLUX = DGSI%XRN - DGSI%XH - DGS%XLEI
  !
  !* wind stress
  !
  DGS%XFMU  = (1 - S%XSIC) * PSFZON + S%XSIC * PSFZON_ICE
  DGSI%XFMU =                                  PSFZON_ICE
  DGS%XFMV  = (1 - S%XSIC) * PSFMER + S%XSIC * PSFMER_ICE
  DGSI%XFMV =                                  PSFMER_ICE
!  
ENDIF
!
!* total storage flux
!
DGS%XGFLUX = DGS%XRN - DGS%XH - DGS%XLE
!
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_BUDGET_SEA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_SURF_BUDGET_SEA
