!     #########
SUBROUTINE DIAG_EVAP_CUMUL_ISBA_n (OSURF_BUDGETC, DE, DECK, DCK, DEK, DK, PEK, &
                             IO, PTSTEP, KSIZE, KPATCH, PRHOA)
!     ###############################################################################
!
!!****  *DIAG_EVAP-ISBA_n * - additional diagnostics for ISBA
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
!!     P. LeMoigne 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    01/2004
!!                     2008      New diag
!!      B. Decharme    2012      New snow diag LESL
!!                               Add carbon fluxes diag
!!                               Add isba water budget diag
!!      B. Decharme  04/2013     add Subsurface runoff if SGH (DIF option only) 
!!                               add sublimation
!!      P Samuelsson   04/2012   MEB
!!------------------------------------------------------------------
!
!
!
USE MODD_DIAG_EVAP_ISBA_n, ONLY : DIAG_EVAP_ISBA_PATCH_t, DIAG_EVAP_ISBA_t
USE MODD_DIAG_n, ONLY : DIAG_PATCH_t, DIAG_t
USE MODD_ISBA_n, ONLY : ISBA_PE_t
USE MODD_ISBA_OPTIONS_n, ONLY : ISBA_OPTIONS_t
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
LOGICAL, INTENT(IN) :: OSURF_BUDGETC
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DE
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DECK
TYPE(DIAG_t), INTENT(INOUT) :: DCK
TYPE(DIAG_EVAP_ISBA_t), INTENT(INOUT) :: DEK
TYPE(DIAG_t), INTENT(INOUT) :: DK
TYPE(ISBA_PE_t), INTENT(INOUT) :: PEK
TYPE(ISBA_OPTIONS_t), INTENT(INOUT) :: IO
!
REAL,    INTENT(IN)               :: PTSTEP        ! time step
INTEGER, INTENT(IN)               :: KSIZE  
INTEGER, INTENT(IN) :: KPATCH
!
REAL,    DIMENSION(:), INTENT(IN) :: PRHOA         ! air density for unit change
!
INTEGER :: JJ
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*      0.2    declarations of local variables
!
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DIAG_EVAP_CUMUL_ISBA_N',0,ZHOOK_HANDLE)
!
IF (OSURF_BUDGETC) THEN
!cdir nodep
  DO JJ=1,KSIZE
     !
     DCK%XRN   (JJ) = DCK%XRN   (JJ) + DK%XRN        (JJ) * PTSTEP
     DCK%XH    (JJ) = DCK%XH    (JJ) + DK%XH         (JJ) * PTSTEP
     DCK%XLE   (JJ) = DCK%XLE   (JJ) + PEK%XLE       (JJ) * PTSTEP
     DCK%XLEI  (JJ) = DCK%XLEI  (JJ) + DK%XLEI       (JJ) * PTSTEP
     DCK%XGFLUX(JJ) = DCK%XGFLUX(JJ) + DK%XGFLUX     (JJ) * PTSTEP
     !
     DECK%XLEG     (JJ) = DECK%XLEG     (JJ) + DEK%XLEG       (JJ) * PTSTEP
     DECK%XLEGI    (JJ) = DECK%XLEGI    (JJ) + DEK%XLEGI      (JJ) * PTSTEP
     DECK%XLEV     (JJ) = DECK%XLEV     (JJ) + DEK%XLEV       (JJ) * PTSTEP
     DECK%XLES     (JJ) = DECK%XLES     (JJ) + DEK%XLES       (JJ) * PTSTEP
     DECK%XLER     (JJ) = DECK%XLER     (JJ) + DEK%XLER       (JJ) * PTSTEP
     DECK%XLETR    (JJ) = DECK%XLETR    (JJ) + DEK%XLETR      (JJ) * PTSTEP
     DCK%XEVAP     (JJ) = DCK%XEVAP     (JJ) + DK%XEVAP       (JJ) * PTSTEP
     DCK%XSUBL     (JJ) = DCK%XSUBL     (JJ) + DK%XSUBL       (JJ) * PTSTEP
     DECK%XDRAIN   (JJ) = DECK%XDRAIN   (JJ) + DEK%XDRAIN     (JJ) * PTSTEP
     DECK%XQSB     (JJ) = DECK%XQSB     (JJ) + DEK%XQSB       (JJ) * PTSTEP
     DECK%XRUNOFF  (JJ) = DECK%XRUNOFF  (JJ) + DEK%XRUNOFF    (JJ) * PTSTEP
     DECK%XHORT    (JJ) = DECK%XHORT    (JJ) + DEK%XHORT      (JJ) * PTSTEP
     DECK%XDRIP    (JJ) = DECK%XDRIP    (JJ) + DEK%XDRIP      (JJ) * PTSTEP
     DECK%XRRVEG   (JJ) = DECK%XRRVEG   (JJ) + DEK%XRRVEG     (JJ) * PTSTEP
     DECK%XMELT    (JJ) = DECK%XMELT    (JJ) + DEK%XMELT      (JJ) * PTSTEP
     DECK%XIFLOOD  (JJ) = DECK%XIFLOOD  (JJ) + DEK%XIFLOOD    (JJ) * PTSTEP
     DECK%XPFLOOD  (JJ) = DECK%XPFLOOD  (JJ) + DEK%XPFLOOD    (JJ) * PTSTEP
     DECK%XLE_FLOOD(JJ) = DECK%XLE_FLOOD(JJ) + DEK%XLE_FLOOD  (JJ) * PTSTEP
     DECK%XLEI_FLOOD (JJ) = DECK%XLEI_FLOOD (JJ) + DEK%XLEI_FLOOD (JJ) * PTSTEP
     DECK%XIRRIG_FLUX(JJ) = DECK%XIRRIG_FLUX(JJ) + DEK%XIRRIG_FLUX(JJ) * PTSTEP
     !
     IF (IO%LMEB_PATCH(KPATCH)) THEN
       DECK%XLEVCV(JJ)  = DECK%XLEVCV (JJ) + DEK%XLEVCV  (JJ) * PTSTEP
       DECK%XLESC (JJ)  = DECK%XLESC  (JJ) + DEK%XLESC   (JJ) * PTSTEP
!       XLETRGVC (JJ)  = XLETRGVC (JJ) + XP_LETRGV   (JJ) * PTSTEP
       DECK%XLETRCV(JJ) = DECK%XLETRCV(JJ) + DEK%XLETRCV (JJ) * PTSTEP
!       XLERGVC  (JJ)  = XLERGVC  (JJ) + XP_LERGV    (JJ) * PTSTEP
       DECK%XLERCV (JJ)  = DECK%XLERCV(JJ) + DEK%XLERCV  (JJ) * PTSTEP
       DECK%XLE_C_A(JJ) = DECK%XLE_C_A(JJ) + DEK%XLE_C_A (JJ) * PTSTEP
       DECK%XLE_V_C(JJ) = DECK%XLE_V_C(JJ) + DEK%XLE_V_C (JJ) * PTSTEP
       DECK%XLE_G_C(JJ) = DECK%XLE_G_C(JJ) + DEK%XLE_G_C (JJ) * PTSTEP
       DECK%XLE_N_C(JJ) = DECK%XLE_N_C(JJ) + DEK%XLE_N_C (JJ) * PTSTEP
       !
       DECK%XSWNET_V  (JJ)  = DECK%XSWNET_V  (JJ)  + DEK%XSWNET_V    (JJ) * PTSTEP
       DECK%XSWNET_G  (JJ)  = DECK%XSWNET_G  (JJ)  + DEK%XSWNET_G    (JJ) * PTSTEP
       DECK%XSWNET_N  (JJ)  = DECK%XSWNET_N  (JJ)  + DEK%XSWNET_N    (JJ) * PTSTEP
       DECK%XSWNET_NS (JJ)  = DECK%XSWNET_NS (JJ)  + DEK%XSWNET_NS   (JJ) * PTSTEP
       DECK%XLWNET_V  (JJ)  = DECK%XLWNET_V  (JJ)  + DEK%XLWNET_V    (JJ) * PTSTEP
       DECK%XLWNET_G  (JJ)  = DECK%XLWNET_G  (JJ)  + DEK%XLWNET_G    (JJ) * PTSTEP
       DECK%XLWNET_N  (JJ)  = DECK%XLWNET_N  (JJ)  + DEK%XLWNET_N    (JJ) * PTSTEP
       DECK%XSWDOWN_GN(JJ)  = DECK%XSWDOWN_GN(JJ)  + DEK%XSWDOWN_GN  (JJ) * PTSTEP
       DECK%XLWDOWN_GN(JJ)  = DECK%XLWDOWN_GN(JJ)  + DEK%XLWDOWN_GN  (JJ) * PTSTEP
       DECK%XH_V_C    (JJ)  = DECK%XH_V_C    (JJ)  + DEK%XH_V_C      (JJ) * PTSTEP
       DECK%XH_G_C    (JJ)  = DECK%XH_G_C    (JJ)  + DEK%XH_G_C      (JJ) * PTSTEP
       DECK%XH_C_A    (JJ)  = DECK%XH_C_A    (JJ)  + DEK%XH_C_A      (JJ) * PTSTEP
       DECK%XH_N_C    (JJ)  = DECK%XH_N_C    (JJ)  + DEK%XH_N_C      (JJ) * PTSTEP
       DECK%XSR_GN    (JJ)  = DECK%XSR_GN    (JJ)  + DEK%XSR_GN      (JJ) * PTSTEP
       DECK%XMELTCV   (JJ)  = DECK%XMELTCV   (JJ)  + DEK%XMELTCV     (JJ) * PTSTEP
       DECK%XFRZCV    (JJ)  = DECK%XFRZCV    (JJ)  + DEK%XFRZCV      (JJ) * PTSTEP
     ENDIF
     !
     DCK%XSWD(JJ)  = DCK%XSWD(JJ) + DK%XSWD(JJ) * PTSTEP
     DCK%XSWU(JJ)  = DCK%XSWU(JJ) + DK%XSWU(JJ) * PTSTEP
     DCK%XLWD(JJ)  = DCK%XLWD(JJ) + DK%XLWD(JJ) * PTSTEP
     DCK%XLWU(JJ)  = DCK%XLWU(JJ) + DK%XLWU(JJ) * PTSTEP
     DCK%XFMU(JJ)  = DCK%XFMU(JJ) + DK%XFMU(JJ) * PTSTEP
     DCK%XFMV(JJ)  = DCK%XFMV(JJ) + DK%XFMV(JJ) * PTSTEP
     !
  END DO
  !
  IF (PEK%TSNOW%SCHEME=='3-L' .OR. PEK%TSNOW%SCHEME=='CRO') THEN
!cdir nodep
     DO JJ=1,KSIZE
        DECK%XLESL    (JJ) = DECK%XLESL    (JJ) + DEK%XLESL    (JJ) * PTSTEP
        DECK%XSNDRIFT (JJ) = DECK%XSNDRIFT (JJ) + DEK%XSNDRIFT (JJ) * PTSTEP
     END DO
  END IF
  !
  IF(IO%CPHOTO/='NON')THEN
!cdir nodep
     DO JJ=1,KSIZE
        !Transform units from kgCO2/kgair m/s to kgCO2/m2
        DECK%XGPP     (JJ)  =  DECK%XGPP(JJ) + DEK%XGPP(JJ) * PTSTEP
        DECK%XRESP_AUTO(JJ) =  DECK%XRESP_AUTO(JJ) + DEK%XRESP_AUTO(JJ) * PTSTEP
        DECK%XRESP_ECO (JJ) =  DECK%XRESP_ECO (JJ) +  DEK%XRESP_ECO(JJ) * PTSTEP
     END DO
  ELSE  
     DECK%XGPP      (:)=0.0
     DECK%XRESP_AUTO(:)=0.0
     DECK%XRESP_ECO (:)=0.0       
  ENDIF
  !
  IF(IO%LGLACIER)THEN
!cdir nodep
    DO JJ=1,KSIZE
       DECK%XICEFLUX(JJ)  = DECK%XICEFLUX(JJ) + DEK%XICEFLUX(JJ) * PTSTEP
    END DO  
  END IF
  !  
  IF(DE%LWATER_BUDGET)THEN
!cdir nodep
     DO JJ=1,KSIZE
        DECK%XDWG   (JJ)  =  DECK%XDWG   (JJ) + DEK%XDWG   (JJ) * PTSTEP
        DECK%XDWGI  (JJ)  =  DECK%XDWGI  (JJ) + DEK%XDWGI  (JJ) * PTSTEP
        DECK%XDWR   (JJ)  =  DECK%XDWR   (JJ) + DEK%XDWR   (JJ) * PTSTEP
        DECK%XDSWE  (JJ)  =  DECK%XDSWE  (JJ) + DEK%XDSWE  (JJ) * PTSTEP
        DECK%XWATBUD(JJ)  =  DECK%XWATBUD(JJ) + DEK%XWATBUD(JJ) * PTSTEP
     END DO
  ENDIF
  !  
END IF
IF (LHOOK) CALL DR_HOOK('DIAG_EVAP_CUMUL_ISBA_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE DIAG_EVAP_CUMUL_ISBA_n
