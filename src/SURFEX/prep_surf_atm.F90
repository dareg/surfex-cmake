!     #################################################################################
SUBROUTINE PREP_SURF_ATM(HPROGRAM,HATMFILE,HATMFILETYPE,HPGDFILE,HPGDFILETYPE)
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
USE MODD_IO_BUFF_n, ONLY : IOB => IO_BUFF
USE MODD_SURF_ATM_SSO_n, ONLY : USS => SURF_ATM_SSO
!
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
USE MODD_DATA_SEAFLUX_n, ONLY : DTS => DATA_SEAFLUX
USE MODD_OCEAN_n, ONLY : O => OCEAN
USE MODD_OCEAN_REL_n, ONLY : OR => OCEAN_REL
USE MODD_SEAFLUX_GRID_n, ONLY : SG => SEAFLUX_GRID
USE MODD_SEAFLUX_n, ONLY : S => SEAFLUX
USE MODD_SEAFLUX_SBL_n, ONLY : SSB => SEAFLUX_SBL
!
USE MODD_FLAKE_GRID_n, ONLY : FG => FLAKE_GRID
USE MODD_FLAKE_n, ONLY : F => FLAKE
USE MODD_FLAKE_SBL_n, ONLY : FSB => FLAKE_SBL
USE MODD_SURF_ATM_GRID_n, ONLY : UG => SURF_ATM_GRID
USE MODD_WATFLUX_GRID_n, ONLY : WG => WATFLUX_GRID
USE MODD_WATFLUX_n, ONLY : W => WATFLUX
USE MODD_WATFLUX_SBL_n, ONLY : WSB => WATFLUX_SBL
!
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
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
IF(U%NDIM_SEA>0) CALL PREP_SEA(DTCO, DTS, O, OR, SG, S, SSB, UG, U, &
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
IF(U%NDIM_NATURE>0) CALL PREP_NATURE(HPROGRAM,YATMFILE,YATMFILETYPE,YPGDFILE,YPGDFILETYPE)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
! URBAN Tile calculations:
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!
IF(U%NDIM_TOWN>0) CALL PREP_TOWN(HPROGRAM,YATMFILE,YATMFILETYPE,YPGDFILE,YPGDFILETYPE)
!
 CALL CLEAR_GRIB_INDEX
!
IF (LHOOK) CALL DR_HOOK('PREP_SURF_ATM',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE PREP_SURF_ATM
