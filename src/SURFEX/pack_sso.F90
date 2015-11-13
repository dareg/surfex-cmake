!     #########
      SUBROUTINE PACK_SSO(USS,HPROGRAM,KMASK,&
                         PAOSIP,PAOSIM,PAOSJP,PAOSJM,&
                         PHO2IP,PHO2IM,PHO2JP,PHO2JM,&
                         PSSO_SLOPE,PSSO_STDEV)
!     #########################################
!
!!****  *PACK_SSO* - routine to initialise the horizontal grid of a scheme
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!	V. Masson   *Meteo France*	
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2003 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_ATM_SSO_n, ONLY : SURF_ATM_SSO_t
!
USE MODI_PACK_SAME_RANK
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!              -------------------------
!
TYPE(SURF_ATM_SSO_t), INTENT(INOUT) :: USS
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM   ! calling program
INTEGER, DIMENSION(:), INTENT(IN) :: KMASK
REAL, DIMENSION(:), INTENT(OUT) :: PSSO_STDEV
REAL, DIMENSION(:), INTENT(OUT) :: PSSO_SLOPE
REAL, DIMENSION(:), INTENT(OUT) :: PAOSIP
REAL, DIMENSION(:), INTENT(OUT) :: PAOSIM
REAL, DIMENSION(:), INTENT(OUT) :: PAOSJP 
REAL, DIMENSION(:), INTENT(OUT) :: PAOSJM
REAL, DIMENSION(:), INTENT(OUT) :: PHO2IP
REAL, DIMENSION(:), INTENT(OUT) :: PHO2IM
REAL, DIMENSION(:), INTENT(OUT) :: PHO2JP
REAL, DIMENSION(:), INTENT(OUT) :: PHO2JM
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!---------------------------------------------------------------------------
!
!*       1.    Reading of type of grid
!              -----------------------
!
IF (LHOOK) CALL DR_HOOK('PACK_SSO',0,ZHOOK_HANDLE)
!
!*    1.      Number of points and packing
!             ----------------------------
!
CALL PACK_SAME_RANK(KMASK,USS%XAOSIP(:),PAOSIP(:))
CALL PACK_SAME_RANK(KMASK,USS%XAOSIM(:),PAOSIM(:))
CALL PACK_SAME_RANK(KMASK,USS%XAOSJP(:),PAOSJP(:))
CALL PACK_SAME_RANK(KMASK,USS%XAOSJM(:),PAOSJM(:))
!
CALL PACK_SAME_RANK(KMASK,USS%XHO2IP(:),PHO2IP(:))
CALL PACK_SAME_RANK(KMASK,USS%XHO2IM(:),PHO2IM(:))
CALL PACK_SAME_RANK(KMASK,USS%XHO2JP(:),PHO2JP(:))
CALL PACK_SAME_RANK(KMASK,USS%XHO2JM(:),PHO2JM(:))
!
CALL PACK_SAME_RANK(KMASK,USS%XSSO_STDEV(:),PSSO_STDEV(:))
CALL PACK_SAME_RANK(KMASK,USS%XSSO_SLOPE(:),PSSO_SLOPE(:))

!
IF (LHOOK) CALL DR_HOOK('PACK_SSO',1,ZHOOK_HANDLE)
!
!---------------------------------------------------------------------------
!
END SUBROUTINE PACK_SSO
