!     #########
       SUBROUTINE DIAG_SURF_BUDGET_WATER (DGW, PTT, PTS, PRHOA, PSFTH, PSFTQ,  PDIR_SW, PSCA_SW, PLW, &
                                          PDIR_ALB, PSCA_ALB, PEMIS, PTRAD, PSFZON, PSFMER   )  
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
!                             Ts instead of Tsrad
!!------------------------------------------------------------------
!
USE MODD_DIAG_n, ONLY : DIAG_t
!
USE MODD_CSTS,           ONLY : XSTEFAN, XLSTT, XLVTT, XCPD
!
! 
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
TYPE(DIAG_t), INTENT(INOUT) :: DGW
!
REAL,               INTENT(IN) :: PTT       ! freezing temperature of water surface
REAL, DIMENSION(:), INTENT(IN) :: PTS       ! surface temperature (K)
REAL, DIMENSION(:), INTENT(IN) :: PRHOA     ! air density
REAL, DIMENSION(:), INTENT(IN) :: PSFTH     ! heat flux
REAL, DIMENSION(:), INTENT(IN) :: PSFTQ     ! water flux
REAL, DIMENSION(:,:),INTENT(IN):: PDIR_SW   ! direct  solar radiation (on horizontal surf.)
!                                           !                                       (W/m2)
REAL, DIMENSION(:,:),INTENT(IN):: PSCA_SW   ! diffuse solar radiation (on horizontal surf.)
!                                           !                                       (W/m2)
REAL, DIMENSION(:), INTENT(IN) :: PLW       ! longwave radiation (on horizontal surf.)
REAL, DIMENSION(:), INTENT(IN) :: PTRAD     ! radiative temperature                 (K)
REAL, DIMENSION(:,:),INTENT(IN):: PDIR_ALB  ! direct albedo for each spectral band  (-)
REAL, DIMENSION(:,:),INTENT(IN):: PSCA_ALB  ! diffuse albedo for each spectral band (-)
REAL, DIMENSION(:), INTENT(IN) :: PEMIS     ! emissivity                            (-)
REAL, DIMENSION(:), INTENT(IN) :: PSFZON    ! zonal friction
REAL, DIMENSION(:), INTENT(IN) :: PSFMER    ! meridional friction
!
!*      0.2    declarations of local variables
!
INTEGER                      :: ISWB ! number of SW bands
INTEGER                      :: JSWB ! loop counter on number of SW bands
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_BUDGET_WATER',0,ZHOOK_HANDLE)
ISWB = SIZE(PDIR_SW,2)
! 
!* total incoming and outgoing SW
!
DO JSWB=1,ISWB
  DGW%XSWBD(:,JSWB) = PDIR_SW(:,JSWB)                    + PSCA_SW(:,JSWB)
  DGW%XSWBU(:,JSWB) = PDIR_SW(:,JSWB) * PDIR_ALB(:,JSWB) + PSCA_SW(:,JSWB) * PSCA_ALB(:,JSWB) 
ENDDO
!
DGW%XSWD(:) = 0.
DGW%XSWU(:) = 0.
DO JSWB=1,ISWB
   DGW%XSWD(:)=DGW%XSWD(:)+DGW%XSWBD(:,JSWB)
   DGW%XSWU(:)=DGW%XSWU(:)+DGW%XSWBU(:,JSWB)
ENDDO
!
!*incoming outgoing LW
!
DGW%XLWD(:) = PLW(:)
DGW%XLWU(:) = PEMIS(:)*XSTEFAN*PTRAD(:)**4 + (1.-PEMIS(:))*PLW(:)
!
!* net radiation
!
DGW%XRN    =   DGW%XSWD(:) - DGW%XSWU(:) + DGW%XLWD(:) - DGW%XLWU(:)
!
!* sensible heat flux
!
DGW%XH     = PSFTH(:)
!
!* latent heat flux
!
WHERE (PTS<PTT  )
  DGW%XLE    = PSFTQ * XLSTT
  DGW%XLEI   = PSFTQ * XLSTT
  DGW%XEVAP  = PSFTQ
  DGW%XSUBL  = PSFTQ
ELSEWHERE
  DGW%XLE    = PSFTQ * XLVTT
  DGW%XLEI   = 0.0
  DGW%XEVAP  = PSFTQ
  DGW%XSUBL  = 0.0
END WHERE
!
!* storage flux
!
DGW%XGFLUX = DGW%XRN - DGW%XH - DGW%XLE
!
!* wind stress
!
DGW%XFMU = PSFZON
!
DGW%XFMV = PSFMER
IF (LHOOK) CALL DR_HOOK('DIAG_SURF_BUDGET_WATER',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_SURF_BUDGET_WATER
