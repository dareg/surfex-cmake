!   ############################################################################
SUBROUTINE PREPS_FOR_MEB_EBUD_RAD(PPS,                                         &
     PLAICV,PSNOWRHO,PSNOWSWE,PSNOWHEAT,                                       &
     PSNOWTEMP,PSNOWDZ,PSCOND,PHEATCAPS,PEMISNOW,PSIGMA_F,PCHIP                )
!   ############################################################################
!
!!****  *PREPS_FOR_MEB_EBUD_RAD*
!!
!!    PURPOSE
!!    -------
!
!     Get preliminary estimates of certain parameters needed for energy budget
!     solution of snowpack, and some other misc inputs needed by radiation 
!     routines for MEB.
!     
!!**  METHOD
!!    ------
!
!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!      
!!    REFERENCE
!!    ---------
!!
!!      
!!    AUTHOR
!!    ------
!!
!!	A. Boone                * CNRM-GAME, Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    02/2011
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!
USE MODD_CSTS,                ONLY : XTT
!
USE MODD_SNOW_PAR,            ONLY : XRHOSMAX_ES, XRHOSMIN_ES, XEMISSN
!
USE MODD_SURF_PAR,            ONLY : XUNDEF
!
USE MODE_SNOW3L,              ONLY : SNOW3LTHRMCOND, SNOW3LSCAP
!
USE MODE_MEB,                 ONLY : MEB_SHIELD_FACTOR
!
USE MODI_SNOW_HEAT_TO_T_WLIQ
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!
!*      0.1    Declaration of Arguments
!
REAL, DIMENSION(:),   INTENT(IN)  :: PLAICV
REAL, DIMENSION(:),   INTENT(IN)  :: PPS
REAL, DIMENSION(:,:), INTENT(IN)  :: PSNOWRHO, PSNOWSWE, PSNOWHEAT

REAL, DIMENSION(:),   INTENT(OUT) :: PSIGMA_F, PCHIP
REAL, DIMENSION(:),   INTENT(OUT) :: PEMISNOW
REAL, DIMENSION(:,:), INTENT(OUT) :: PSNOWDZ, PSCOND, PHEATCAPS, PSNOWTEMP
!
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PLAICV,1)) :: ZPSNA
REAL, DIMENSION(SIZE(PSNOWRHO,1),SIZE(PSNOWRHO,2)) :: ZSNOWLIQ, ZSNOWRHO
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!----------------------------------------------------
! 0) Initialization
!
IF (LHOOK) CALL DR_HOOK('PREPS_FOR_MEB_EBUD_RAD',0,ZHOOK_HANDLE)
!
! First, since snow might not exist, check using snow density as a
! proxy (since several computations below depend on this value). Note
! that if snow doesn't exist, the variables below depending on snow
! density are either zero (e.g. PSNOWDZ) or unused (e.g. PSCOND)
!
ZSNOWRHO(:,:)    = PSNOWRHO(:,:)
WHERE(PSNOWRHO(:,:)==XUNDEF)
   ZSNOWRHO(:,:) = XRHOSMIN_ES
ENDWHERE

! Snow layer thicknesses (m)

PSNOWDZ(:,:)     = PSNOWSWE(:,:)/ZSNOWRHO(:,:)             

! snow temperature (K) and liquid water content (kg m-3)

CALL SNOW_HEAT_TO_T_WLIQ(PSNOWHEAT,ZSNOWRHO,PSNOWTEMP,ZSNOWLIQ)

! Snow thermal conductivity:

PSCOND(:,:)      = SNOW3LTHRMCOND(ZSNOWRHO,PSNOWTEMP,PPS)  ! W m-1 K-1

! Snow heat capacity:

PHEATCAPS(:,:)   = SNOW3LSCAP(ZSNOWRHO)                    ! J m-3 K-1

! View factor: (1 - shielding factor)

ZPSNA(:)          = 0.
PCHIP(:)          = MEB_SHIELD_FACTOR(PLAICV,ZPSNA)
PSIGMA_F(:)       = 1.0 - PCHIP(:)

! snow emissivity

PEMISNOW(:)       = XEMISSN

IF (LHOOK) CALL DR_HOOK('PREPS_FOR_MEB_EBUD_RAD',1,ZHOOK_HANDLE)

END SUBROUTINE PREPS_FOR_MEB_EBUD_RAD
