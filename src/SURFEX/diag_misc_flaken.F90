!     #########
       SUBROUTINE DIAG_MISC_FLAKE_n (DGMF, F)
!     ###############################################################################
!
!!****  *DIAG_MISC-FLAKE_n * - additional diagnostics for FLake
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
!!     P. Le Moigne 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/2005
!!------------------------------------------------------------------
!
!
USE MODD_FLAKE_n, ONLY : FLAKE_t
USE MODD_DIAG_MISC_FLAKE_n, ONLY : DIAG_MISC_FLAKE_t
!
USE MODD_SURF_PAR,           ONLY : XUNDEF
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(FLAKE_t), INTENT(INOUT) :: F
TYPE(DIAG_MISC_FLAKE_t), INTENT(INOUT) :: DGMF
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(DGMF%XZW_PROFILE),SIZE(F%XT_WML)) :: ZCSI      ! Vertical normalized coordinate
REAL, DIMENSION(SIZE(DGMF%XZW_PROFILE),SIZE(F%XT_WML)) :: ZSHAPE    ! Shape function
!
INTEGER         :: IZW
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DIAG_MISC_FLAKE_N',0,ZHOOK_HANDLE)
!
!* Flake temperature profile
!
DGMF%XTW_PROFILE(:,:) = XUNDEF
!
IF (DGMF%LWATER_PROFILE) THEN
!
   DO IZW=1,SIZE(DGMF%XZW_PROFILE)
      WHERE (F%XWATER_DEPTH(:)==F%XH_ML(:))
         ZCSI(IZW,:) = 0.
      ELSEWHERE
         ZCSI(IZW,:) = (DGMF%XZW_PROFILE(IZW) - F%XH_ML(:))/(F%XWATER_DEPTH(:) - F%XH_ML(:))
      END WHERE
      ZSHAPE(IZW,:) = (40./3.*F%XCT-20./3.)*ZCSI(IZW,:)   +     (18.-30.*F%XCT)*ZCSI(IZW,:)**2 &
                       + (20.*F%XCT-12.)   *ZCSI(IZW,:)**3+(5./3.-10./3.*F%XCT)*ZCSI(IZW,:)**4  
   END DO
!
   DO IZW=1,SIZE(DGMF%XZW_PROFILE)
      WHERE (F%XH_ML(:) >= DGMF%XZW_PROFILE(IZW))
         DGMF%XTW_PROFILE(IZW,:) =  F%XT_WML(:) 
      ELSEWHERE (F%XWATER_DEPTH(:) >= DGMF%XZW_PROFILE(IZW)) 
         DGMF%XTW_PROFILE(IZW,:) = F%XT_WML(:) - (F%XT_WML(:) - F%XT_BOT(:)) * ZSHAPE(IZW,:)
      END WHERE
   END DO
!
END IF
!
IF (LHOOK) CALL DR_HOOK('DIAG_MISC_FLAKE_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_MISC_FLAKE_n
