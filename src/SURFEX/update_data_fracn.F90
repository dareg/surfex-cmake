!     #########
      SUBROUTINE UPDATE_DATA_FRAC_n(PDATA_NATURE,PDATA_TOWN,PDATA_GARDEN,OGARDEN, &
                                    PDATA_BLD, PDATA_WALL_O_HOR                   )
!     #########################
!
!!**** *INI_DATA_FRAC* takes into account gardens into natural vegetation
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
!!
!!
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
!!    Original    09/2011
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_DATA_COVER_PAR, ONLY : JPCOVER
USE MODD_DATA_COVER_n, ONLY : XDATA_NATURE, XDATA_TOWN, XDATA_GARDEN, &
                              XDATA_BLD, XDATA_WALL_O_HOR, LGARDEN
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
REAL, DIMENSION(:), INTENT(IN)  :: PDATA_NATURE
REAL, DIMENSION(:), INTENT(IN)  :: PDATA_TOWN
REAL, DIMENSION(:), INTENT(IN)  :: PDATA_GARDEN
LOGICAL,            INTENT(IN)  :: OGARDEN
REAL, DIMENSION(:), INTENT(IN)  :: PDATA_BLD
REAL, DIMENSION(:), INTENT(IN)  :: PDATA_WALL_O_HOR
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('UPDATE_DATA_FRAC_n',0,ZHOOK_HANDLE)

IF (.NOT.ASSOCIATED(XDATA_NATURE)) THEN
  ALLOCATE(XDATA_NATURE    (JPCOVER))
  ALLOCATE(XDATA_TOWN      (JPCOVER))
  ALLOCATE(XDATA_GARDEN    (JPCOVER))
  ALLOCATE(XDATA_BLD       (JPCOVER))
  ALLOCATE(XDATA_WALL_O_HOR(JPCOVER))
ENDIF
!
LGARDEN = OGARDEN
!
XDATA_NATURE     = PDATA_NATURE
XDATA_TOWN       = PDATA_TOWN
XDATA_GARDEN     = PDATA_GARDEN
XDATA_BLD        = PDATA_BLD
XDATA_WALL_O_HOR = PDATA_WALL_O_HOR
!
IF (.NOT. OGARDEN) THEN
  XDATA_NATURE     = PDATA_NATURE + PDATA_GARDEN * PDATA_TOWN
  XDATA_TOWN       = PDATA_TOWN   * ( 1. - PDATA_GARDEN)
  XDATA_GARDEN     = 0.
  XDATA_BLD        = PDATA_BLD / (1. - PDATA_GARDEN)
  XDATA_WALL_O_HOR = PDATA_WALL_O_HOR / (1. - PDATA_GARDEN)
END IF
!
IF (LHOOK) CALL DR_HOOK('UPDATE_DATA_FRAC_n',1,ZHOOK_HANDLE)
!-------------------------------------------------------------------------------

END SUBROUTINE UPDATE_DATA_FRAC_n
