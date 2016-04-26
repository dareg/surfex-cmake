!##################################
SUBROUTINE CUMUL_DIAG_TEB_n (DMTC, DMT, TOP, PTSTEP)
!##################################
!
!
!!****  *CUMUL_DIAG_TEB_n*  
!!
!!    PURPOSE
!!    -------
!      Cumulates some diagnostics for TEB
!     
!!**  METHOD
!!    ------
!
!!    EXTERNAL
!!    --------
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
!!      C. de Munck       * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    02/2013
!!                  08/2013 (V. Masson) adds solar panels
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_DIAG_MISC_TEB_n, ONLY : DIAG_MISC_TEB_1P_t
USE MODD_TEB_OPTION_n, ONLY : TEB_OPTIONS_t
!
USE MODD_SURF_PAR,        ONLY :  XUNDEF
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(DIAG_MISC_TEB_1P_t), INTENT(INOUT) :: DMTC
TYPE(DIAG_MISC_TEB_1P_t), INTENT(INOUT) :: DMT
TYPE(TEB_OPTIONS_t), INTENT(INOUT) :: TOP
!
      REAL,               INTENT(IN) :: PTSTEP            ! time step
!
!*      0.2    declarations of local variables
!
INTEGER :: JI 

REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
!       0.     Initialization
!              --------------
IF (LHOOK) CALL DR_HOOK('CUMUL_DIAG_TEB_N',0,ZHOOK_HANDLE)
!
!       1.     Time-cumulated diagnostics for TEB
!              ----------------------------------
!
DO JI=1,SIZE(DMT%XRUNOFF_ROOF,1)
!
 IF (TOP%LSOLAR_PANEL) THEN
    IF (DMT%XTHER_PROD_BLD(JI) .NE. XUNDEF) THEN
     DMTC%XTHER_PROD_BLD(JI)     =  DMTC%XTHER_PROD_BLD(JI)     + DMT%XTHER_PROD_BLD(JI)  * PTSTEP
    ENDIF
    !
    IF (DMT%XPHOT_PROD_BLD(JI) .NE. XUNDEF) THEN
     DMTC%XPHOT_PROD_BLD(JI)     =  DMTC%XPHOT_PROD_BLD(JI)     + DMT%XPHOT_PROD_BLD(JI)  * PTSTEP
    ENDIF
 END IF

 IF (TOP%CBEM == 'BEM') THEN
    IF (DMT%XHVAC_COOL(JI) .NE. XUNDEF) THEN
     DMTC%XHVAC_COOL(JI)     =  DMTC%XHVAC_COOL(JI)       + DMT%XHVAC_COOL(JI)        * PTSTEP
    ENDIF
    !
    IF (DMT%XHVAC_HEAT(JI) .NE. XUNDEF) THEN
     DMTC%XHVAC_HEAT(JI)     =  DMTC%XHVAC_HEAT(JI)       + DMT%XHVAC_HEAT(JI)        * PTSTEP
    ENDIF
 ENDIF
 !
 IF (DMT%XRUNOFF_TOWN(JI) .NE. XUNDEF) THEN
  DMTC%XRUNOFF_TOWN(JI)      =  DMTC%XRUNOFF_TOWN(JI)     + DMT%XRUNOFF_TOWN(JI)      * PTSTEP
 ENDIF
 !
 IF (DMT%XRUNOFF_GARDEN(JI) .NE. XUNDEF) THEN
  DMTC%XRUNOFF_GARDEN(JI)    =  DMTC%XRUNOFF_GARDEN(JI)   + DMT%XRUNOFF_GARDEN(JI)    * PTSTEP
 ENDIF
 !
 IF (DMT%XRUNOFF_ROAD(JI) .NE. XUNDEF) THEN
  DMTC%XRUNOFF_ROAD(JI)      =  DMTC%XRUNOFF_ROAD(JI)     + DMT%XRUNOFF_ROAD(JI)      * PTSTEP
 ENDIF
 !
 IF (DMT%XRUNOFF_ROOF(JI) .NE. XUNDEF) THEN 
  DMTC%XRUNOFF_ROOF(JI)      =  DMTC%XRUNOFF_ROOF(JI)     + DMT%XRUNOFF_ROOF(JI)      * PTSTEP
 ENDIF
 !
 IF (DMT%XRUNOFF_STRLROOF(JI) .NE. XUNDEF) THEN
  DMTC%XRUNOFF_STRLROOF(JI)  =  DMTC%XRUNOFF_STRLROOF(JI) + DMT%XRUNOFF_STRLROOF(JI)  * PTSTEP
 ENDIF
 !
 IF (DMT%XDRAIN_GARDEN(JI) .NE. XUNDEF) THEN
   DMTC%XDRAIN_GARDEN(JI)    =  DMTC%XDRAIN_GARDEN(JI)    + DMT%XDRAIN_GARDEN(JI)     * PTSTEP
 ENDIF
 !
 IF (DMT%XIRRIG_GARDEN(JI) .NE. XUNDEF) THEN
   DMTC%XIRRIG_GARDEN(JI)    =  DMTC%XIRRIG_GARDEN(JI)    + DMT%XIRRIG_GARDEN(JI)     * PTSTEP
 ENDIF
 !
 IF (DMT%XIRRIG_ROAD(JI) .NE. XUNDEF) THEN
   DMTC%XIRRIG_ROAD(JI)      =  DMTC%XIRRIG_ROAD(JI)      + DMT%XIRRIG_ROAD(JI)       * PTSTEP
 ENDIF
 !
 IF (TOP%LGREENROOF) THEN 
    IF (DMT%XRUNOFF_GREENROOF(JI) .NE. XUNDEF) THEN
     DMTC%XRUNOFF_GREENROOF(JI) =  DMTC%XRUNOFF_GREENROOF(JI)+ DMT%XRUNOFF_GREENROOF(JI) * PTSTEP
    ENDIF
    !
    IF (DMT%XDRAIN_GREENROOF(JI) .NE. XUNDEF) THEN
     DMTC%XDRAIN_GREENROOF(JI)  =  DMTC%XDRAIN_GREENROOF(JI) + DMT%XDRAIN_GREENROOF(JI)  * PTSTEP
    ENDIF
    !
    IF (DMT%XIRRIG_GREENROOF(JI) .NE. XUNDEF) THEN
     DMTC%XIRRIG_GREENROOF(JI)  =  DMTC%XIRRIG_GREENROOF(JI) + DMT%XIRRIG_GREENROOF(JI)  * PTSTEP
    ENDIF
 ENDIF
 !
ENDDO
!
IF (LHOOK) CALL DR_HOOK('CUMUL_DIAG_TEB_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE CUMUL_DIAG_TEB_n
