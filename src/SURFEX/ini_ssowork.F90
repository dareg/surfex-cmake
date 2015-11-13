!     #########
      SUBROUTINE INI_SSOWORK(PMESHLENGTH,PDLAT,PDLON)
!     ###############################################
!
!!**** *INI_SSOWORK* initializes and allocate work arrays for SSO reading
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!   
!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    10/12/97
!!
!
USE MODD_SURFEX_MPI, ONLY : NRANK, NPIO
USE MODD_PGDWORK,  ONLY : NSSO, N3D_ALL, X3D_ALL
USE MODD_SURF_PAR, ONLY : NUNDEF, XUNDEF
USE MODD_SURFEX_MPI, ONLY : NINDEX
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!----------------------------------------------------------------------------
!
!*    0.1    Declaration of arguments
!            ------------------------
!
REAL, OPTIONAL, INTENT(IN) :: PMESHLENGTH ! average mesh length in degrees
REAL, OPTIONAL, INTENT(IN) :: PDLAT       ! input file mesh size (in latitude,  degrees)
REAL, OPTIONAL, INTENT(IN) :: PDLON       ! input file mesh size (in longitude, degrees)

INTEGER :: IDIMF
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!----------------------------------------------------------------------------
!
!
!*    1.     Adapt subgrid mesh to input file resolution
!            -------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('INI_SSOWORK',0,ZHOOK_HANDLE)
!
IDIMF = SIZE(NINDEX)
!
IF (PRESENT(PMESHLENGTH) .AND. PRESENT(PDLAT) .AND. PRESENT(PDLON)) THEN
  IF (PDLAT/= XUNDEF .AND. PDLON /= XUNDEF) THEN
    NSSO = NINT( 2. * PMESHLENGTH / (PDLAT + PDLON) )
    NSSO = MAX(NSSO,3)
    NSSO = MIN(NSSO,10)
  ELSE
    NSSO = 10
  END IF
ELSE
  NSSO = 10
END IF
!
!----------------------------------------------------------------------------
!
!*    2.     Allocate subgrid arrays
!            -----------------------
!
IF (ALLOCATED(X3D_ALL)) DEALLOCATE(X3D_ALL)
!
ALLOCATE(X3D_ALL(IDIMF,NSSO,NSSO))
X3D_ALL(:,:,:) = -XUNDEF
!
IF (ALLOCATED(N3D_ALL)) DEALLOCATE(N3D_ALL)
!
ALLOCATE(N3D_ALL(IDIMF,NSSO,NSSO))
N3D_ALL(:,:,:) = 0
!
IF (LHOOK) CALL DR_HOOK('INI_SSOWORK',1,ZHOOK_HANDLE)
!
!----------------------------------------------------------------------------
!
END SUBROUTINE INI_SSOWORK
