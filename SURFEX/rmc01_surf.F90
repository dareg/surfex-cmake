!     #########
      SUBROUTINE RMC01_SURF(PZ, PLMO, PLK, PLEPS)
!     ##############################################################
!
!!****  *RMC01_SURF* -
!!
!!    PURPOSE
!!    -------
!!    This routine modifies the mixing and dissipative length near the SBL.
!!    (Redelsperger, Mahe and Carlotti, 2001)
!!
!!**  METHOD
!!    ------
!!
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
!!      Book 2
!!
!!    AUTHOR
!!    ------
!!
!!      V. Masson  - Meteo-France -
!!
!!    MODIFICATIONS
!!    -------------
!!     Original     07/2006
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_PAR,    ONLY : XUNDEF
USE MODD_CSTS,        ONLY : XKARMAN
USE MODD_CANOPY_TURB, ONLY : XALPSBL, XCMFS, XCED
!
USE MODE_SBLS
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declaration of arguments
!              ------------------------
!
REAL, DIMENSION(:,:),   INTENT(IN)  :: PZ   ! altitude of full levels
REAL, DIMENSION(:),     INTENT(IN)  :: PLMO ! Monin Obuhkov length
REAL, DIMENSION(:,:),   INTENT(OUT) :: PLK  ! Mixing length
REAL, DIMENSION(:,:),   INTENT(OUT) :: PLEPS! Dissipative length
!
!*       0.2   Declaration of local variables
!              ------------------------------
!
INTEGER :: JK        ! loop counter
!
REAL, DIMENSION(SIZE(PZ,1),SIZE(PZ,2)) :: ZZ_O_LMO ! height / LMO
REAL, DIMENSION(SIZE(PZ,1),SIZE(PZ,2)) :: ZPHIM    ! MO function
                                                   ! for stress
REAL, DIMENSION(SIZE(PZ,1),SIZE(PZ,2)) :: ZPHIE    ! MO function
                                                   ! for TKE
REAL, DIMENSION(SIZE(PZ,1),SIZE(PZ,2)) :: ZL       ! SBL length
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!*     1. MO quantities
!         -------------
!
! z/LMO
IF (LHOOK) CALL DR_HOOK('RMC01_SURF',0,ZHOOK_HANDLE)
DO JK=1,SIZE(PZ,2)
  WHERE (PLMO(:)==XUNDEF)
    ZZ_O_LMO(:,JK)=0.
  ELSEWHERE
    ZZ_O_LMO(:,JK)=PZ(:,JK)/PLMO(:)
  END WHERE
END DO
ZZ_O_LMO(:,:) = MAX(ZZ_O_LMO(:,:),-10.)
ZZ_O_LMO(:,:) = MIN(ZZ_O_LMO(:,:), 10.)
!
!
! MO function for stress
ZPHIM(:,:) = BUSINGER_PHIM(ZZ_O_LMO)
!
! MO function for TKE
ZPHIE(:,:) = BUSINGER_PHIE(ZZ_O_LMO)
!
!-------------------------------------------------------------------------------
!
!*     2. Modification of the mixing length
!         ---------------------------------
!
PLK(:,:) =  XKARMAN/SQRT(XALPSBL)/XCMFS &
            * PZ(:,:)/(ZPHIM(:,:)**2*SQRT(ZPHIE(:,:)))  
!
!-------------------------------------------------------------------------------
!
!*     3. Modification of the dissipative length
!         --------------------------------------
!
ZL = PLK * (XALPSBL**(3./2.)*XKARMAN*XCED) &
           / (XKARMAN/SQRT(XALPSBL)/XCMFS)  
!
WHERE (ZZ_O_LMO<0.)
  ZL = ZL/(1.-1.9*ZZ_O_LMO)
ELSEWHERE
  ZL = ZL/(1.-0.3*SQRT(ZZ_O_LMO))
ENDWHERE
!
PLEPS(:,:)=ZL(:,:)
IF (LHOOK) CALL DR_HOOK('RMC01_SURF',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE RMC01_SURF
