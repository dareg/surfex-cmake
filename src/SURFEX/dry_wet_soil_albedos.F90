!     ##################################################################
      SUBROUTINE DRY_WET_SOIL_ALBEDOS( KK  )  
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
USE MODD_ISBA_n, ONLY : ISBA_K_t
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
TYPE(ISBA_K_t), INTENT(INOUT) :: KK
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DRY_WET_SOIL_ALBEDOS',0,ZHOOK_HANDLE)
KK%XALBVIS_DRY(:) = 0.05 +  (   0.05 + MAX(0.30 * KK%XSAND(:,1), 0.10) )  &
                         * ( 1. - 0.9 * ( KK%XVEGTYPE(:,NVT_C3)        &
                                        + KK%XVEGTYPE(:,NVT_C4)        &
                                        + KK%XVEGTYPE(:,NVT_IRR)       &
                                        + KK%XVEGTYPE(:,NVT_GRAS)      &
                                        + KK%XVEGTYPE(:,NVT_TROG)      &
                                        + KK%XVEGTYPE(:,NVT_PARK)      &
                                        + KK%XVEGTYPE(:,NVT_TRBE)      &
                                        + KK%XVEGTYPE(:,NVT_BONE)      &
                                        + KK%XVEGTYPE(:,NVT_TEBD)      &
                                        + KK%XVEGTYPE(:,NVT_TRBD)      & 
                                        + KK%XVEGTYPE(:,NVT_TEBE)      &
                                        + KK%XVEGTYPE(:,NVT_TENE)      &
                                        + KK%XVEGTYPE(:,NVT_BOBD)      &
                                        + KK%XVEGTYPE(:,NVT_BOND)      &
                                        + KK%XVEGTYPE(:,NVT_BOGR)      &
                                        + KK%XVEGTYPE(:,NVT_SHRB))**2 )  
!
KK%XALBNIR_DRY(:) = KK%XALBVIS_DRY(:) + 0.10
!
KK%XALBUV_DRY (:) = 0.06 + 0.14 * KK%XSAND(:,1)
!
KK%XALBVIS_WET(:) = KK%XALBVIS_DRY(:) / 2.
KK%XALBNIR_WET(:) = KK%XALBNIR_DRY(:) / 2.
KK%XALBUV_WET (:) = KK%XALBUV_DRY (:) / 2.
!
IF (LHOOK) CALL DR_HOOK('DRY_WET_SOIL_ALBEDOS',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DRY_WET_SOIL_ALBEDOS
