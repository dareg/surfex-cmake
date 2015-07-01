!     ###########################################################
      SUBROUTINE ZOOM_PGD_INLAND_WATER(HPROGRAM,HINIFILE,HINIFILETYPE, &
                                       HFILE,HFILETYPE,OECOCLIMAP)
!     ###########################################################

!!
!!    PURPOSE
!!    -------
!!   This program prepares the physiographic data fields.
!!
!!    METHOD
!!    ------
!!   
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!    V. Masson                   Meteo-France
!!
!!    MODIFICATION
!!    ------------
!!
!!    Original     13/10/03
!!    B. Decharme  02/2014  Add LRM_RIVER
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
USE MODD_WATFLUX_GRID_n, ONLY : WG => WATFLUX_GRID
USE MODD_WATFLUX_n, ONLY : W => WATFLUX
!
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_PGD_FLAKE
USE MODI_PGD_WATFLUX
IMPLICIT NONE
!
!*    0.1    Declaration of dummy arguments
!            ------------------------------
!
CHARACTER(LEN=6),     INTENT(IN)  :: HPROGRAM    ! program calling
CHARACTER(LEN=28),    INTENT(IN)  :: HINIFILE    ! input atmospheric file name
CHARACTER(LEN=6),     INTENT(IN)  :: HINIFILETYPE! input atmospheric file type
CHARACTER(LEN=28),    INTENT(IN)  :: HFILE       ! output file name
CHARACTER(LEN=6),     INTENT(IN)  :: HFILETYPE   ! output file type
LOGICAL,              INTENT(IN)  :: OECOCLIMAP
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
LOGICAL :: LRM_RIVER ! dummy keys
!
!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ZOOM_PGD_INLAND_WATER',0,ZHOOK_HANDLE)
IF (U%CWATER=='NONE  ') THEN
  IF (LHOOK) CALL DR_HOOK('ZOOM_PGD_INLAND_WATER',1,ZHOOK_HANDLE)
  RETURN
ELSE IF (U%CWATER=='FLUX  ') THEN
  IF (LHOOK) CALL DR_HOOK('ZOOM_PGD_INLAND_WATER',1,ZHOOK_HANDLE)
  RETURN
ELSE IF (U%CWATER=='WATFLX') THEN
  CALL PGD_WATFLUX(DTCO, U, WG, W, &
                   HPROGRAM)
ELSE IF (U%CWATER=='FLAKE ') THEN
  LRM_RIVER=.TRUE.
  CALL PGD_FLAKE(HPROGRAM,OECOCLIMAP,LRM_RIVER)
END IF
IF (LHOOK) CALL DR_HOOK('ZOOM_PGD_INLAND_WATER',1,ZHOOK_HANDLE)
!
!_______________________________________________________________________________
!
END SUBROUTINE ZOOM_PGD_INLAND_WATER
