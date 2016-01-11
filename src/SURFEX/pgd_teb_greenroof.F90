!     #########
      SUBROUTINE PGD_TEB_GREENROOF (DTCO, UG, U, USS, GRM, TG, &
                                    HPROGRAM)
!     ##############################################################
!
!!**** *PGD_TEB_GREENROOF* monitor for averaging and interpolations of TEB GR physiographic fields
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
!!    C.de Munck & A. Lemonsu        Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original    07/2011
!!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
!
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SURF_ATM_SSO_n, ONLY : SURF_ATM_SSO_t
USE MODD_SURFEX_n, ONLY : TEB_GREENROOF_MODEL_t
USE MODD_GRID_n, ONLY : GRID_t
!
USE MODD_PGD_GRID,             ONLY : NL
USE MODD_DATA_COVER_PAR,       ONLY : NVEGTYPE
!
USE MODI_PGD_TEB_GREENROOF_PAR
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ABOR1_SFX
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SURF_ATM_SSO_t), INTENT(INOUT) :: USS
TYPE(TEB_GREENROOF_MODEL_t), INTENT(INOUT) :: GRM
TYPE(GRID_t), INTENT(INOUT) :: TG
!
 CHARACTER(LEN=6), INTENT(IN)  :: HPROGRAM   ! Type of program
!                                           ! F if all parameters must be specified
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
!*    0.3    Declaration of namelists
!            ------------------------
!
REAL(KIND=JPRB)          :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PGD_TEB_GREENROOF',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!*    1.      ISBA specific fields for green roofs
!             ------------------------------------
!
! for green roofs, CISBA = DIF / CSCOND = 'DEF '
GRM%TV%O%CISBA  = 'DIF'
GRM%TV%O%CSCOND = 'PL98 ' ! CSCOND_GR = 'PL98' !begin test 29092011 ! normalement pas besoin
GRM%TV%O%CHORT  = 'DEF '
GRM%TV%O%CKSAT  = 'DEF '
GRM%TV%O%LSOC   = .FALSE.
GRM%TV%O%LTR_ML = .FALSE.
!
ALLOCATE(GRM%TV%P%XRUNOFFB(TG%NDIM))
ALLOCATE(GRM%TV%P%XWDRAIN (TG%NDIM))
!
GRM%TV%P%XRUNOFFB(:) = 0.5 
GRM%TV%P%XWDRAIN (:) = 0.0
!
GRM%DTI%NTIME = 12
!
GRM%TV%O%LPAR = .TRUE.
!
CALL PGD_TEB_GREENROOF_PAR(DTCO, GRM%DTI, UG, U, USS, GRM%TV%O, GRM%TV%P, TG, &
                           HPROGRAM)
!
!
IF (LHOOK) CALL DR_HOOK('PGD_TEB_GREENROOF',1,ZHOOK_HANDLE)
!
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE PGD_TEB_GREENROOF
