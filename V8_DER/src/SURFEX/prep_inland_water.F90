!     #########
SUBROUTINE PREP_INLAND_WATER (DTCO, IOB, USS, &
                               FG, F, FSB, UG, U, WG, W, WSB, &
                              HPROGRAM,HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!     #################################################################################
!
!!****  *PREP_INLAND_WATER* - Chooses the surface schemes for lakes preparation
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
!!     S. Malardel
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!------------------------------------------------------------------
!
!
!
!
!
!
!
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
USE MODD_SURF_ATM_SSO_n, ONLY : SURF_ATM_SSO_t
!
USE MODD_FLAKE_GRID_n, ONLY : FLAKE_GRID_t
USE MODD_FLAKE_n, ONLY : FLAKE_t
USE MODD_FLAKE_SBL_n, ONLY : FLAKE_SBL_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_WATFLUX_GRID_n, ONLY : WATFLUX_GRID_t
USE MODD_WATFLUX_n, ONLY : WATFLUX_t
USE MODD_WATFLUX_SBL_n, ONLY : WATFLUX_SBL_t
!
USE MODI_PREP_WATFLUX
USE MODI_PREP_FLAKE
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
TYPE(SURF_ATM_SSO_t), INTENT(INOUT) :: USS
!
TYPE(FLAKE_GRID_t), INTENT(INOUT) :: FG
TYPE(FLAKE_t), INTENT(INOUT) :: F
TYPE(FLAKE_SBL_t), INTENT(INOUT) :: FSB
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(WATFLUX_GRID_t), INTENT(INOUT) :: WG
TYPE(WATFLUX_t), INTENT(INOUT) :: W
TYPE(WATFLUX_SBL_t), INTENT(INOUT) :: WSB
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM  ! program calling surf. schemes
 CHARACTER(LEN=28),  INTENT(IN)  :: HATMFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),   INTENT(IN)  :: HATMFILETYPE! type of the Atmospheric file
 CHARACTER(LEN=28),  INTENT(IN)  :: HPGDFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),   INTENT(IN)  :: HPGDFILETYPE! type of the Atmospheric file
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('PREP_INLAND_WATER',0,ZHOOK_HANDLE)
IF (U%CWATER=='WATFLX') THEN
  CALL PREP_WATFLUX(DTCO, IOB, &
                    UG, U, WG, W, WSB, &
                    HPROGRAM,HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
ELSE IF (U%CWATER=='FLAKE ') THEN
  CALL PREP_FLAKE(DTCO, IOB, USS, &
                  FG, F, FSB, UG, U, &
                  HPROGRAM,HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
END IF
IF (LHOOK) CALL DR_HOOK('PREP_INLAND_WATER',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE PREP_INLAND_WATER
