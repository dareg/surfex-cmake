!##################################
SUBROUTINE CUMUL_DIAG_TEB_n (DGCT, DGMT, TOP, &
                             PTSTEP)
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
USE MODD_DIAG_CUMUL_TEB_n, ONLY : DIAG_CUMUL_TEB_t
USE MODD_DIAG_MISC_TEB_n, ONLY : DIAG_MISC_TEB_t
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
TYPE(DIAG_CUMUL_TEB_t), INTENT(INOUT) :: DGCT
TYPE(DIAG_MISC_TEB_t), INTENT(INOUT) :: DGMT
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
DO JI=1,SIZE(DGMT%XRUNOFF_ROOF,1)
!
 IF (TOP%LSOLAR_PANEL) THEN
    IF (DGMT%XTHER_PROD_BLD(JI) .NE. XUNDEF) THEN
     DGCT%XTHER_PROD_BLDC(JI)     =  DGCT%XTHER_PROD_BLDC(JI)     + DGMT%XTHER_PROD_BLD(JI)  * PTSTEP
    ENDIF
    !
    IF (DGMT%XPHOT_PROD_BLD(JI) .NE. XUNDEF) THEN
     DGCT%XPHOT_PROD_BLDC(JI)     =  DGCT%XPHOT_PROD_BLDC(JI)     + DGMT%XPHOT_PROD_BLD(JI)  * PTSTEP
    ENDIF
 END IF

 IF (TOP%CBEM == 'BEM') THEN
    IF (DGMT%XHVAC_COOL(JI) .NE. XUNDEF) THEN
     DGCT%XHVACC_COOL(JI)     =  DGCT%XHVACC_COOL(JI)       + DGMT%XHVAC_COOL(JI)        * PTSTEP
    ENDIF
    !
    IF (DGMT%XHVAC_HEAT(JI) .NE. XUNDEF) THEN
     DGCT%XHVACC_HEAT(JI)     =  DGCT%XHVACC_HEAT(JI)       + DGMT%XHVAC_HEAT(JI)        * PTSTEP
    ENDIF
 ENDIF
 !
 IF (DGMT%XRUNOFF_TOWN(JI) .NE. XUNDEF) THEN
  DGCT%XRUNOFFC_TOWN(JI)      =  DGCT%XRUNOFFC_TOWN(JI)     + DGMT%XRUNOFF_TOWN(JI)      * PTSTEP
 ENDIF
 !
 IF (DGMT%XRUNOFF_GARDEN(JI) .NE. XUNDEF) THEN
  DGCT%XRUNOFFC_GARDEN(JI)    =  DGCT%XRUNOFFC_GARDEN(JI)   + DGMT%XRUNOFF_GARDEN(JI)    * PTSTEP
 ENDIF
 !
 IF (DGMT%XRUNOFF_ROAD(JI) .NE. XUNDEF) THEN
  DGCT%XRUNOFFC_ROAD(JI)      =  DGCT%XRUNOFFC_ROAD(JI)     + DGMT%XRUNOFF_ROAD(JI)      * PTSTEP
 ENDIF
 !
 IF (DGMT%XRUNOFF_ROOF(JI) .NE. XUNDEF) THEN 
  DGCT%XRUNOFFC_ROOF(JI)      =  DGCT%XRUNOFFC_ROOF(JI)     + DGMT%XRUNOFF_ROOF(JI)      * PTSTEP
 ENDIF
 !
 IF (DGMT%XRUNOFF_STRLROOF(JI) .NE. XUNDEF) THEN
  DGCT%XRUNOFFC_STRLROOF(JI)  =  DGCT%XRUNOFFC_STRLROOF(JI) + DGMT%XRUNOFF_STRLROOF(JI)  * PTSTEP
 ENDIF
 !
 IF (DGMT%XDRAIN_GARDEN(JI) .NE. XUNDEF) THEN
   DGCT%XDRAINC_GARDEN(JI)    =  DGCT%XDRAINC_GARDEN(JI)    + DGMT%XDRAIN_GARDEN(JI)     * PTSTEP
 ENDIF
 !
 IF (DGMT%XIRRIG_GARDEN(JI) .NE. XUNDEF) THEN
   DGCT%XIRRIGC_GARDEN(JI)    =  DGCT%XIRRIGC_GARDEN(JI)    + DGMT%XIRRIG_GARDEN(JI)     * PTSTEP
 ENDIF
 !
 IF (DGMT%XIRRIG_ROAD(JI) .NE. XUNDEF) THEN
   DGCT%XIRRIGC_ROAD(JI)      =  DGCT%XIRRIGC_ROAD(JI)      + DGMT%XIRRIG_ROAD(JI)       * PTSTEP
 ENDIF
 !
 IF (TOP%LGREENROOF) THEN 
    IF (DGMT%XRUNOFF_GREENROOF(JI) .NE. XUNDEF) THEN
     DGCT%XRUNOFFC_GREENROOF(JI) =  DGCT%XRUNOFFC_GREENROOF(JI)+ DGMT%XRUNOFF_GREENROOF(JI) * PTSTEP
    ENDIF
    !
    IF (DGMT%XDRAIN_GREENROOF(JI) .NE. XUNDEF) THEN
     DGCT%XDRAINC_GREENROOF(JI)  =  DGCT%XDRAINC_GREENROOF(JI) + DGMT%XDRAIN_GREENROOF(JI)  * PTSTEP
    ENDIF
    !
    IF (DGMT%XIRRIG_GREENROOF(JI) .NE. XUNDEF) THEN
     DGCT%XIRRIGC_GREENROOF(JI)  =  DGCT%XIRRIGC_GREENROOF(JI) + DGMT%XIRRIG_GREENROOF(JI)  * PTSTEP
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
