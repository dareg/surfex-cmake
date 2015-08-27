!     #################################################################################
SUBROUTINE PREP_SURF_ATM (B, BOP, DTCO, DTS, FG, F, FSB, IOB, ICP, IG, I, &
                          O, OR, SG, S, SSB, UG, U, USS, TCP, TGD, TGDO, &
                          TGDPE, TGDP, TGR, TGRO, TGRPE, TGRP, TG, T, TOP, TVG, WG, &
                          W, WSB, &
                          HPROGRAM,HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
!     #################################################################################
!
!!****  *PREP_SURF_ATM* - driver for surface fields preparation
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
!!     V. Masson 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!      P. Le Moigne 10/2005, Phasage Arome
!!------------------------------------------------------------------
!

!
!
!
!
!
!
!
!
USE MODD_BEM_n, ONLY : BEM_t
USE MODD_BEM_OPTION_n, ONLY : BEM_OPTIONS_t
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DATA_SEAFLUX_n, ONLY : DATA_SEAFLUX_t
USE MODD_FLAKE_GRID_n, ONLY : FLAKE_GRID_t
USE MODD_FLAKE_n, ONLY : FLAKE_t
USE MODD_FLAKE_SBL_n, ONLY : FLAKE_SBL_t
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
USE MODD_ISBA_CANOPY_n, ONLY : ISBA_CANOPY_t
USE MODD_ISBA_GRID_n, ONLY : ISBA_GRID_t
USE MODD_ISBA_n, ONLY : ISBA_t
USE MODD_OCEAN_n, ONLY : OCEAN_t
USE MODD_OCEAN_REL_n, ONLY : OCEAN_REL_t
USE MODD_SEAFLUX_GRID_n, ONLY : SEAFLUX_GRID_t
USE MODD_SEAFLUX_n, ONLY : SEAFLUX_t
USE MODD_SEAFLUX_SBL_n, ONLY : SEAFLUX_SBL_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SURF_ATM_SSO_n, ONLY : SURF_ATM_SSO_t
USE MODD_TEB_CANOPY_n, ONLY : TEB_CANOPY_t
USE MODD_TEB_GARDEN_n, ONLY : TEB_GARDEN_t
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TEB_GARDEN_OPTIONS_t
USE MODD_TEB_GARDEN_PGD_EVOL_n, ONLY : TEB_GARDEN_PGD_EVOL_t
USE MODD_TEB_GARDEN_PGD_n, ONLY : TEB_GARDEN_PGD_t
USE MODD_TEB_GREENROOF_n, ONLY : TEB_GREENROOF_t
USE MODD_TEB_GREENROOF_OPTION_n, ONLY : TEB_GREENROOF_OPTIONS_t
USE MODD_TEB_GREENROOF_PGD_EVOL_n, ONLY : TEB_GREENROOF_PGD_EVOL_t
USE MODD_TEB_GREENROOF_PGD_n, ONLY : TEB_GREENROOF_PGD_t
USE MODD_TEB_GRID_n, ONLY : TEB_GRID_t
USE MODD_TEB_n, ONLY : TEB_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_OPTIONS_t
USE MODD_WATFLUX_GRID_n, ONLY : WATFLUX_GRID_t
USE MODD_WATFLUX_n, ONLY : WATFLUX_t
USE MODD_WATFLUX_SBL_n, ONLY : WATFLUX_SBL_t
!
USE MODI_PREP_NATURE
USE MODI_PREP_SEA
USE MODI_PREP_INLAND_WATER
USE MODI_PREP_TOWN
!
USE MODE_READ_GRIB
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_SURF_VERSION
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(BEM_t), INTENT(INOUT) :: B
TYPE(BEM_OPTIONS_t), INTENT(INOUT) :: BOP
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DATA_SEAFLUX_t), INTENT(INOUT) :: DTS
TYPE(FLAKE_GRID_t), INTENT(INOUT) :: FG
TYPE(FLAKE_t), INTENT(INOUT) :: F
TYPE(FLAKE_SBL_t), INTENT(INOUT) :: FSB
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
TYPE(ISBA_CANOPY_t), INTENT(INOUT) :: ICP
TYPE(ISBA_GRID_t), INTENT(INOUT) :: IG
TYPE(ISBA_t), INTENT(INOUT) :: I
TYPE(OCEAN_t), INTENT(INOUT) :: O
TYPE(OCEAN_REL_t), INTENT(INOUT) :: OR
TYPE(SEAFLUX_GRID_t), INTENT(INOUT) :: SG
TYPE(SEAFLUX_t), INTENT(INOUT) :: S
TYPE(SEAFLUX_SBL_t), INTENT(INOUT) :: SSB
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(SURF_ATM_SSO_t), INTENT(INOUT) :: USS
TYPE(TEB_CANOPY_t), INTENT(INOUT) :: TCP
TYPE(TEB_GARDEN_t), INTENT(INOUT) :: TGD
TYPE(TEB_GARDEN_OPTIONS_t), INTENT(INOUT) :: TGDO
TYPE(TEB_GARDEN_PGD_EVOL_t), INTENT(INOUT) :: TGDPE
TYPE(TEB_GARDEN_PGD_t), INTENT(INOUT) :: TGDP
TYPE(TEB_GREENROOF_t), INTENT(INOUT) :: TGR
TYPE(TEB_GREENROOF_OPTIONS_t), INTENT(INOUT) :: TGRO
TYPE(TEB_GREENROOF_PGD_EVOL_t), INTENT(INOUT) :: TGRPE
TYPE(TEB_GREENROOF_PGD_t), INTENT(INOUT) :: TGRP
TYPE(TEB_GRID_t), INTENT(INOUT) :: TG
TYPE(TEB_t), INTENT(INOUT) :: T
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_VEG_OPTIONS_t), INTENT(INOUT) :: TVG
TYPE(WATFLUX_GRID_t), INTENT(INOUT) :: WG
TYPE(WATFLUX_t), INTENT(INOUT) :: W
TYPE(WATFLUX_SBL_t), INTENT(INOUT) :: WSB
!
 CHARACTER(LEN=6),  INTENT(IN) :: HPROGRAM    ! program calling surf. schemes
 CHARACTER(LEN=28), INTENT(IN) :: HATMFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),  INTENT(IN) :: HATMFILETYPE! type of the Atmospheric file
 CHARACTER(LEN=28), INTENT(IN) :: HPGDFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6),  INTENT(IN) :: HPGDFILETYPE! type of the Atmospheric file
!
!*      0.2    declarations of local variables
 CHARACTER(LEN=28)               :: YATMFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6)                :: YATMFILETYPE! type of the Atmospheric file
 CHARACTER(LEN=28)               :: YPGDFILE    ! name of the Atmospheric file
 CHARACTER(LEN=6)                :: YPGDFILETYPE! type of the Atmospheric file
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!
!-------------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('PREP_SURF_ATM',0,ZHOOK_HANDLE)
 CALL SURF_VERSION
!-------------------------------------------------------------------------------------
!
IF ( LEN_TRIM(HATMFILE)>0 ) THEN
  YATMFILE=HATMFILE
ELSE
  YATMFILE='                            '
ENDIF
!
IF ( LEN_TRIM(HPGDFILE)>0 ) THEN
  YPGDFILE=HPGDFILE
ELSE
  YPGDFILE='                            '
ENDIF
!
IF (  LEN_TRIM(HATMFILETYPE)>0 ) THEN
  YATMFILETYPE=HATMFILETYPE
ELSE
  YATMFILETYPE='      '
ENDIF
!
IF (  LEN_TRIM(HPGDFILETYPE)>0 ) THEN
  YPGDFILETYPE=HPGDFILETYPE
ELSE
  YPGDFILETYPE='      '
ENDIF
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! SEA Tile calculations:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
IF(U%NDIM_SEA>0) CALL PREP_SEA(IOB, &
   DTCO, DTS, O, OR, SG, S, SSB, UG, U, &
   HPROGRAM,YATMFILE,YATMFILETYPE,YPGDFILE,YPGDFILETYPE)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! INLAND WATER Tile calculations:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
IF(U%NDIM_WATER>0) CALL PREP_INLAND_WATER(DTCO, IOB, USS, &
   FG, F, FSB, UG, U, WG, W, WSB, &
   HPROGRAM,YATMFILE,YATMFILETYPE,YPGDFILE,YPGDFILETYPE)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! NATURAL SURFACE Tile calculations:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
IF(U%NDIM_NATURE>0) CALL PREP_NATURE(DTCO, IOB, ICP, IG, I, UG, U, USS, &
   HPROGRAM,YATMFILE,YATMFILETYPE,YPGDFILE,YPGDFILETYPE)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! URBAN Tile calculations:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
IF(U%NDIM_TOWN>0) CALL PREP_TOWN(B, BOP, DTCO, IOB, IG, I, UG, U, USS, TCP, TGD, &
                      TGDO, TGDPE, TGDP, TGR, TGRO, TGRPE, TGRP, TG, T, TOP, TVG, &
                      HPROGRAM,YATMFILE,YATMFILETYPE,YPGDFILE,YPGDFILETYPE)
!
 CALL CLEAR_GRIB_INDEX
!
IF (LHOOK) CALL DR_HOOK('PREP_SURF_ATM',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE PREP_SURF_ATM
