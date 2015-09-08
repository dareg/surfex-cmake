!     #########
      SUBROUTINE PGD_TEB_GREENROOF (DTCO, DTGR, UG, U, USS, TGRO, TGRP, TG, &
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
USE MODD_DATA_TEB_GREENROOF_n, ONLY : DATA_TEB_GREENROOF_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SURF_ATM_SSO_n, ONLY : SURF_ATM_SSO_t
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TEB_GREENROOF_OPTIONS_t
USE MODD_TEB_GREENROOF_PGD_n, ONLY : TEB_GREENROOF_PGD_t
USE MODD_TEB_GRID_n, ONLY : TEB_GRID_t
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
TYPE(DATA_TEB_GREENROOF_t), INTENT(INOUT) :: DTGR
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SURF_ATM_SSO_t), INTENT(INOUT) :: USS
TYPE(TEB_GREENROOF_OPTIONS_t), INTENT(INOUT) :: TGRO
TYPE(TEB_GREENROOF_PGD_t), INTENT(INOUT) :: TGRP
TYPE(TEB_GRID_t), INTENT(INOUT) :: TG
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
TGRO%CISBA_GR  = 'DIF'
TGRO%CSCOND_GR = 'PL98 ' ! CSCOND_GR = 'PL98' !begin test 29092011 ! normalement pas besoin
TGRO%CHORT_GR  = 'DEF '
TGRO%CKSAT_GR  = 'DEF '
TGRO%LSOC_GR   = .FALSE.
TGRO%LTR_ML_GR = .FALSE.
!
ALLOCATE(TGRP%XRUNOFFB_GR(TG%NDIM))
ALLOCATE(TGRP%XWDRAIN_GR (TG%NDIM))
!
TGRP%XRUNOFFB_GR(:) = 0.5 
TGRP%XWDRAIN_GR (:) = 0.0
!
TGRO%NTIME_GR = 12
CALL PGD_TEB_GREENROOF_PAR(DTCO, DTGR, UG, U, USS, TGRO, TG, &
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
