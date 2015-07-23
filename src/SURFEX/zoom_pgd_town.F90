!     ###########################################################
      SUBROUTINE ZOOM_PGD_TOWN(HPROGRAM,HINIFILE,HINIFILETYPE,HFILE,HFILETYPE,OECOCLIMAP,OGARDEN)
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
USE MODD_BEM_OPTION_n, ONLY : BOP => BEM_OPTIONS
USE MODD_BLD_DESCRIPTION_n, ONLY : BDD => BLD_DESC
USE MODD_DATA_BEM_n, ONLY : DTB => DATA_BEM
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
USE MODD_DATA_TEB_n, ONLY : DTT => DATA_TEB
USE MODD_IO_BUFF_n, ONLY : IOB => IO_BUFF
USE MODD_SURF_ATM_GRID_n, ONLY : UG => SURF_ATM_GRID
USE MODD_TEB_GARDEN_OPTION_n, ONLY : TGDO => TEB_GARDEN_OPTIONS
USE MODD_TEB_GARDEN_PGD_n, ONLY : TGDP => TEB_GARDEN_PGD
USE MODD_TEB_GRID_n, ONLY : TG => TEB_GRID
USE MODD_TEB_OPTION_n, ONLY : TOP => TEB_OPTIONS
USE MODD_TEB_VEG_n, ONLY : TVG => TEB_VEG_OPTIONS
!
USE MODD_SURF_ATM_n, ONLY : U => SURF_ATM
!
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
