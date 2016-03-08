!     ##################################################################
      SUBROUTINE DRY_WET_SOIL_ALBEDOS(PSAND, PCLAY, PVEGTYPE, IP  )  
!     ##################################################################
!
!!****  *DRY_WET_SOIL_ALBEDOS*  
!!
!!    PURPOSE
!!    -------
!       computes the albedo of bare soil, for dry or wet conditions
!
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
!!      
!!    AUTHOR
!!    ------
!!      V. Masson           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    17/12/99 
!       
!      (V. Masson)  16/02/01 Better fit with ISLSCP2; 
!                                            Ba et al 2001; 
!                                            Pinty et al 2000
!      (V. Masson) 01/2004  Add UV albedo
!      (R. Alkama) 05/2012  Add 7 new vegtype (19 rather than 12)
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ISBA_INIT_n, ONLY : ISBA_INIT_PGD_t
!
USE MODD_DATA_COVER_PAR, ONLY : NVT_PARK, NVT_TEBD, NVT_BONE, NVT_TRBE, NVT_TRBD, &
                                NVT_TEBE, NVT_TENE, NVT_BOBD, NVT_BOND, NVT_SHRB, &
                                NVT_C3, NVT_C4, NVT_IRR, NVT_GRAS, NVT_BOGR,      &
                                NVT_TROG                   
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!              -------------------------
!
REAL, DIMENSION(:), INTENT(IN)  :: PSAND       ! sand fraction
REAL, DIMENSION(:), INTENT(IN)  :: PCLAY       ! clay fraction
REAL, DIMENSION(:,:), INTENT(IN):: PVEGTYPE    ! vegetation type
!
TYPE(ISBA_INIT_PGD_t), INTENT(INOUT) :: IP
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DRY_WET_SOIL_ALBEDOS',0,ZHOOK_HANDLE)
IP%XALBVIS_DRY(:) = 0.05 +  (   0.05 + MAX(0.30 * PSAND(:), 0.10) )  &
                         * ( 1. - 0.9 * ( PVEGTYPE(:,NVT_C3)        &
                                        + PVEGTYPE(:,NVT_C4)        &
                                        + PVEGTYPE(:,NVT_IRR)       &
                                        + PVEGTYPE(:,NVT_GRAS)      &
                                        + PVEGTYPE(:,NVT_TROG)      &
                                        + PVEGTYPE(:,NVT_PARK)      &
                                        + PVEGTYPE(:,NVT_TRBE)      &
                                        + PVEGTYPE(:,NVT_BONE)      &
                                        + PVEGTYPE(:,NVT_TEBD)      &
                                        + PVEGTYPE(:,NVT_TRBD)      & 
                                        + PVEGTYPE(:,NVT_TEBE)      &
                                        + PVEGTYPE(:,NVT_TENE)      &
                                        + PVEGTYPE(:,NVT_BOBD)      &
                                        + PVEGTYPE(:,NVT_BOND)      &
                                        + PVEGTYPE(:,NVT_BOGR)      &
                                        + PVEGTYPE(:,NVT_SHRB))**2 )  
!
IP%XALBNIR_DRY(:) = IP%XALBVIS_DRY(:) + 0.10
!
IP%XALBUV_DRY (:) = 0.06 + 0.14 * PSAND(:)
!
IP%XALBVIS_WET(:) = IP%XALBVIS_DRY(:) / 2.
IP%XALBNIR_WET(:) = IP%XALBNIR_DRY(:) / 2.
IP%XALBUV_WET (:) = IP%XALBUV_DRY (:) / 2.
IF (LHOOK) CALL DR_HOOK('DRY_WET_SOIL_ALBEDOS',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DRY_WET_SOIL_ALBEDOS
