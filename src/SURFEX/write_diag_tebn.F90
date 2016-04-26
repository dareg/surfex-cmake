!     #########
SUBROUTINE WRITE_DIAG_TEB_n (DTCO, HSELECT, U, TM, GDM, GRM, HPROGRAM,HWRITE)
!     ###############################################################################
!
!!****  *WRITE_DIAG_TEB_n * - diagnostics for TEB
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
!!------------------------------------------------------------------
!
!
!
USE MODD_DATA_COVER_n, ONLY : DATA_COVER_t
USE MODD_SURF_ATM_n, ONLY : SURF_ATM_t
USE MODD_SURFEX_n, ONLY : TEB_MODEL_t
USE MODD_SURFEX_n, ONLY : TEB_GARDEN_MODEL_t
USE MODD_SURFEX_n, ONLY : TEB_GREENROOF_MODEL_t
!
USE MODD_SURF_PAR,   ONLY : XUNDEF
!
USE MODI_GOTO_WRAPPER_TEB_PATCH
USE MODI_WRITE_DIAG_SEB_TEB_n
USE MODI_WRITE_DIAG_MISC_TEB_n
USE MODI_WRITE_DIAG_PGD_TEB_n
USE MODI_WRITE_DIAG_PGD_GRDN_n
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
!
TYPE(DATA_COVER_t), INTENT(INOUT) :: DTCO
 CHARACTER(LEN=*), DIMENSION(:), INTENT(IN) :: HSELECT
TYPE(SURF_ATM_t), INTENT(INOUT) :: U
TYPE(TEB_MODEL_t), INTENT(INOUT) :: TM
TYPE(TEB_GARDEN_MODEL_t), INTENT(INOUT) :: GDM
TYPE(TEB_GREENROOF_MODEL_t), INTENT(INOUT) :: GRM
!
 CHARACTER(LEN=6),   INTENT(IN)  :: HPROGRAM ! program calling surf. schemes
 CHARACTER(LEN=3),   INTENT(IN)  :: HWRITE   ! 'PGD' : only physiographic fields are written
!                                           ! 'ALL' : all fields are written
!
!*      0.2    declarations of local variables
!
INTEGER         :: JTEB_PATCH
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_TEB_N',0,ZHOOK_HANDLE)
IF (HWRITE/='PGD') THEN
!        
   IF (TM%TD%O%XDIAG_TSTEP==XUNDEF .OR. &
         ABS(NINT(TM%TOP%TTIME%TIME/TM%TD%O%XDIAG_TSTEP)*TM%TD%O%XDIAG_TSTEP-TM%TOP%TTIME%TIME)<1.E-3 ) THEN
      CALL WRITE_DIAG_SEB_TEB_n(DTCO, HSELECT, U, TM%CHT, TM%TD%O, TM%TD%D, TM%TD%DUT, HPROGRAM)
      DO JTEB_PATCH=1,TM%TOP%NTEB_PATCH
        CALL GOTO_WRAPPER_TEB_PATCH(JTEB_PATCH, DMTC=TM%TD%DMTC, DMT=TM%TD%DMT, T=TM%T)
        CALL WRITE_DIAG_MISC_TEB_n(DTCO, HSELECT, U, TM%TD%DMTC%CUR, TM%TD%DMT%CUR, TM%TD%MTO, &
                                        TM%T%CUR, TM%TOP, HPROGRAM,JTEB_PATCH)
      END DO      
   END IF
!
ENDIF
!
IF (TM%TD%O%LPGD) THEN
  IF (TM%TD%O%XDIAG_TSTEP==XUNDEF .OR. &
          ABS(NINT(TM%TOP%TTIME%TIME/TM%TD%O%XDIAG_TSTEP)*TM%TD%O%XDIAG_TSTEP-TM%TOP%TTIME%TIME)<1.E-3 ) THEN
    IF (ASSOCIATED(TM%T%CUR%XBLD)) THEN
      !DO JTEB_PATCH=1,TM%TOP%NTEB_PATCH
      !  CALL GOTO_WRAPPER_TEB_PATCH(JTEB_PATCH, B=TM%B, T=TM%T)            
        CALL WRITE_DIAG_PGD_TEB_n(DTCO, HSELECT, U, TM%B%CUR, TM%BOP, TM%T%CUR, TM%TOP, TM%TPN, HPROGRAM)
      !ENDDO
      IF (TM%TOP%LGARDEN) &
        CALL WRITE_DIAG_PGD_GRDN_n(DTCO, HSELECT, U, TM%TD%MTO%LSURF_DIAG_ALBEDO, &
                        GDM%S, GDM%P, GDM%NPE%AL(1), GDM%O, HPROGRAM)
    ENDIF
  END IF
END IF
IF (LHOOK) CALL DR_HOOK('WRITE_DIAG_TEB_N',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------------
!
END SUBROUTINE WRITE_DIAG_TEB_n
