!##################################
SUBROUTINE CUMUL_DIAG_TEB_n(PTSTEP)
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
USE MODD_SURF_PAR,        ONLY :  XUNDEF
USE MODD_TEB_OPTION_n,    ONLY :  CBEM, LGREENROOF, LSOLAR_PANEL
USE MODD_TEB_PANEL_n,     ONLY :  XFRAC_PANEL 
USE MODD_DIAG_MISC_TEB_n, ONLY :  XHVAC_COOL,        &
                                  XHVAC_HEAT,        &
                                  XRUNOFF_TOWN,      &
                                  XRUNOFF_GARDEN,    &
                                  XRUNOFF_ROAD,      &
                                  XRUNOFF_ROOF,      &
                                  XRUNOFF_STRLROOF,  &
                                  XRUNOFF_GREENROOF, &
                                  XDRAIN_GREENROOF,  &
                                  XDRAIN_GARDEN,     &
                                  XIRRIG_GREENROOF,  &
                                  XIRRIG_GARDEN,     &
                                  XIRRIG_ROAD,       &
                                  XTHER_PROD_BLD,    &
                                  XPHOT_PROD_BLD  
USE MODD_DIAG_CUMUL_TEB_n, ONLY : DGCT => DIAG_CUMUL_TEB
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
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
DO JI=1,SIZE(XRUNOFF_ROOF,1)
!
 IF (LSOLAR_PANEL) THEN
    IF (XTHER_PROD_BLD(JI) .NE. XUNDEF) THEN
     DGCT%XTHER_PROD_BLDC(JI)     =  DGCT%XTHER_PROD_BLDC(JI)     + XTHER_PROD_BLD(JI)  * PTSTEP
    ENDIF
    !
    IF (XPHOT_PROD_BLD(JI) .NE. XUNDEF) THEN
     DGCT%XPHOT_PROD_BLDC(JI)     =  DGCT%XPHOT_PROD_BLDC(JI)     + XPHOT_PROD_BLD(JI)  * PTSTEP
    ENDIF
 END IF

 IF (CBEM == 'BEM') THEN
    IF (XHVAC_COOL(JI) .NE. XUNDEF) THEN
     DGCT%XHVACC_COOL(JI)     =  DGCT%XHVACC_COOL(JI)       + XHVAC_COOL(JI)        * PTSTEP
    ENDIF
    !
    IF (XHVAC_HEAT(JI) .NE. XUNDEF) THEN
     DGCT%XHVACC_HEAT(JI)     =  DGCT%XHVACC_HEAT(JI)       + XHVAC_HEAT(JI)        * PTSTEP
    ENDIF
 ENDIF
 !
 IF (XRUNOFF_TOWN(JI) .NE. XUNDEF) THEN
  DGCT%XRUNOFFC_TOWN(JI)      =  DGCT%XRUNOFFC_TOWN(JI)     + XRUNOFF_TOWN(JI)      * PTSTEP
 ENDIF
 !
 IF (XRUNOFF_GARDEN(JI) .NE. XUNDEF) THEN
  DGCT%XRUNOFFC_GARDEN(JI)    =  DGCT%XRUNOFFC_GARDEN(JI)   + XRUNOFF_GARDEN(JI)    * PTSTEP
 ENDIF
 !
 IF (XRUNOFF_ROAD(JI) .NE. XUNDEF) THEN
  DGCT%XRUNOFFC_ROAD(JI)      =  DGCT%XRUNOFFC_ROAD(JI)     + XRUNOFF_ROAD(JI)      * PTSTEP
 ENDIF
 !
 IF (XRUNOFF_ROOF(JI) .NE. XUNDEF) THEN 
  DGCT%XRUNOFFC_ROOF(JI)      =  DGCT%XRUNOFFC_ROOF(JI)     + XRUNOFF_ROOF(JI)      * PTSTEP
 ENDIF
 !
 IF (XRUNOFF_STRLROOF(JI) .NE. XUNDEF) THEN
  DGCT%XRUNOFFC_STRLROOF(JI)  =  DGCT%XRUNOFFC_STRLROOF(JI) + XRUNOFF_STRLROOF(JI)  * PTSTEP
 ENDIF
 !
 IF (XDRAIN_GARDEN(JI) .NE. XUNDEF) THEN
   DGCT%XDRAINC_GARDEN(JI)    =  DGCT%XDRAINC_GARDEN(JI)    + XDRAIN_GARDEN(JI)     * PTSTEP
 ENDIF
 !
 IF (XIRRIG_GARDEN(JI) .NE. XUNDEF) THEN
   DGCT%XIRRIGC_GARDEN(JI)    =  DGCT%XIRRIGC_GARDEN(JI)    + XIRRIG_GARDEN(JI)     * PTSTEP
 ENDIF
 !
 IF (XIRRIG_ROAD(JI) .NE. XUNDEF) THEN
   DGCT%XIRRIGC_ROAD(JI)      =  DGCT%XIRRIGC_ROAD(JI)      + XIRRIG_ROAD(JI)       * PTSTEP
 ENDIF
 !
 IF (LGREENROOF) THEN 
    IF (XRUNOFF_GREENROOF(JI) .NE. XUNDEF) THEN
     DGCT%XRUNOFFC_GREENROOF(JI) =  DGCT%XRUNOFFC_GREENROOF(JI)+ XRUNOFF_GREENROOF(JI) * PTSTEP
    ENDIF
    !
    IF (XDRAIN_GREENROOF(JI) .NE. XUNDEF) THEN
     DGCT%XDRAINC_GREENROOF(JI)  =  DGCT%XDRAINC_GREENROOF(JI) + XDRAIN_GREENROOF(JI)  * PTSTEP
    ENDIF
    !
    IF (XIRRIG_GREENROOF(JI) .NE. XUNDEF) THEN
     DGCT%XIRRIGC_GREENROOF(JI)  =  DGCT%XIRRIGC_GREENROOF(JI) + XIRRIG_GREENROOF(JI)  * PTSTEP
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
