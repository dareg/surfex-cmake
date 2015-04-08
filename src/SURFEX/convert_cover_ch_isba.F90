!     #########
      SUBROUTINE CONVERT_COVER_CH_ISBA   (PCOVER,OCOVER,PSOILRC_SO2,PSOILRC_O3)
!     ##############################################################
!
!!**** *CONVERT_COVER* convert surface cover classes into secondary 
!!                     physiographic variables for ISBA
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
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
!!    Original   01/2004
!     
!
!----------------------------------------------------------------------------
!
!*    0.     DECLARATION
!            -----------
!
USE MODD_DATA_COVER_n, ONLY : DTCO => DATA_COVER
!
USE MODD_DATA_COVER,     ONLY : XDATA_SOILRC_SO2, XDATA_SOILRC_O3 
!
USE MODD_DATA_COVER_PAR, ONLY : NVEGTYPE, JPCOVER
!
USE MODI_AV_PGD
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*    0.1    Declaration of arguments
!            ------------------------
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PCOVER
LOGICAL, DIMENSION(:), INTENT(IN)   :: OCOVER

REAL, DIMENSION(:,:),   INTENT(OUT)   :: PSOILRC_SO2
REAL, DIMENSION(:,:),   INTENT(OUT)   :: PSOILRC_O3
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!
!*    0.2    Declaration of local variables
!            ------------------------------
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('CONVERT_COVER_CH_ISBA',0,ZHOOK_HANDLE)
!
IF (ASSOCIATED(DTCO%XDATA_WEIGHT)) DEALLOCATE(DTCO%XDATA_WEIGHT)
!
 CALL AV_PGD (PSOILRC_SO2 ,PCOVER ,XDATA_SOILRC_SO2(:,:) ,'NAT','ARI',OCOVER,KDECADE=1)
 CALL AV_PGD (PSOILRC_O3  ,PCOVER ,XDATA_SOILRC_O3 (:,:) ,'NAT','ARI',OCOVER,KDECADE=1)
!
IF (ASSOCIATED(DTCO%XDATA_WEIGHT)) DEALLOCATE(DTCO%XDATA_WEIGHT)
!
IF (LHOOK) CALL DR_HOOK('CONVERT_COVER_CH_ISBA',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE CONVERT_COVER_CH_ISBA
