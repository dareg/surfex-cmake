!     ###########################################################
      SUBROUTINE ZOOM_PGD_TOWN (BOP, BDD, DTB, DTCO, DTT, IOB, UG, U, TGDO, TGDP, TG, &
                                TOP, TVG, &
                                HPROGRAM,HINIFILE,HINIFILETYPE,HFILE,HFILETYPE,OECOCLIMAP,OGARDEN)
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
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
!
!
!
!
USE MODD_BEM_OPTION_n, ONLY : BEM_OPTIONS_t
USE MODD_BLD_DESCRIPTION_n, ONLY : BLD_DESC_t
USE MODD_DATA_BEM_n, ONLY : DATA_BEM_t
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_DATA_TEB_n, ONLY : DATA_TEB_t
USE MODD_IO_BUFF_n, ONLY : IO_BUFF_t
USE MODD_SURF_ATM_GRID_n, ONLY : SURF_ATM_GRID_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TEB_GARDEN_OPTIONS_t
USE MODD_TEB_GARDEN_PGD_n, ONLY : TEB_GARDEN_PGD_t
USE MODD_TEB_GRID_n, ONLY : TEB_GRID_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
USE MODD_TEB_VEG_n, ONLY : TEB_VEG_OPTIONS_t
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
USE MODI_ZOOM_PGD_TEB
!
IMPLICIT NONE
!
!*    0.1    Declaration of dummy arguments
!            ------------------------------
!
!
TYPE(BEM_OPTIONS_t), INTENT(INOUT) :: BOP
TYPE(BLD_DESC_t), INTENT(INOUT) :: BDD
TYPE(DATA_BEM_t), INTENT(INOUT) :: DTB
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
TYPE(DATA_TEB_t), INTENT(INOUT) :: DTT
TYPE(IO_BUFF_t), INTENT(INOUT) :: IOB
TYPE(SURF_ATM_GRID_t), INTENT(INOUT) :: UG
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(TEB_GARDEN_OPTIONS_t), INTENT(INOUT) :: TGDO
TYPE(TEB_GARDEN_PGD_t), INTENT(INOUT) :: TGDP
TYPE(TEB_GRID_t), INTENT(INOUT) :: TG
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
TYPE(TEB_VEG_OPTIONS_t), INTENT(INOUT) :: TVG
!
 CHARACTER(LEN=6),     INTENT(IN)  :: HPROGRAM    ! program calling
 CHARACTER(LEN=28),    INTENT(IN)  :: HINIFILE    ! input atmospheric file name
 CHARACTER(LEN=6),     INTENT(IN)  :: HINIFILETYPE! input atmospheric file type
 CHARACTER(LEN=28),    INTENT(IN)  :: HFILE       ! output file name
 CHARACTER(LEN=6),     INTENT(IN)  :: HFILETYPE   ! output file type
LOGICAL,              INTENT(IN)  :: OECOCLIMAP  ! flag to use ecoclimap
LOGICAL,              INTENT(IN)  :: OGARDEN     ! flag to use garden
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
!------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ZOOM_PGD_TOWN',0,ZHOOK_HANDLE)
IF (U%CTOWN=='NONE  ') THEN
  IF (LHOOK) CALL DR_HOOK('ZOOM_PGD_TOWN',1,ZHOOK_HANDLE)
  RETURN
ELSE IF (U%CTOWN=='FLUX  ') THEN
  IF (LHOOK) CALL DR_HOOK('ZOOM_PGD_TOWN',1,ZHOOK_HANDLE)
  RETURN
ELSE IF (U%CTOWN=='TEB   ') THEN
  CALL ZOOM_PGD_TEB(BOP, BDD, DTB, DTCO, DTT, IOB, UG, U, TGDO, TGDP, TG, &
                               TOP, TVG, &
                    HPROGRAM,HINIFILE,HINIFILETYPE,OECOCLIMAP,OGARDEN)
END IF
IF (LHOOK) CALL DR_HOOK('ZOOM_PGD_TOWN',1,ZHOOK_HANDLE)
!
!_______________________________________________________________________________
!
END SUBROUTINE ZOOM_PGD_TOWN
