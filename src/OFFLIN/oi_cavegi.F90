SUBROUTINE OI_CAVEGI(VGAT1,VGAT2,VGAT3,VGBT1,VGBT2,VGBT3,VGCT1,VGCT2,       &
                      VGAH1,VGAH2,VGAH3,VGBH1,VGBH2,VGBH3,VGCH1,VGCH2,       &
                      SIGT2MP,SIGHP2,LSGOBS)   
!****------------------------------------------------------------------------
USE MODD_ASSIM, ONLY : SIGHP1, SIGT2MR, SIGH2MR, SIGT2MO, SIGH2MO, REPS3
!
!
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INTEGER :: J
! 
REAL ,INTENT(OUT) :: VGAT1(24),VGAT2(24),VGAT3(24)
REAL ,INTENT(OUT) :: VGBT1(24),VGBT2(24),VGBT3(24)
REAL ,INTENT(OUT) :: VGCT1(24),VGCT2(24)
REAL ,INTENT(OUT) :: VGAH1(24),VGAH2(24),VGAH3(24)
REAL ,INTENT(OUT) :: VGBH1(24),VGBH2(24),VGBH3(24)
REAL ,INTENT(OUT) :: VGCH1(24),VGCH2(24)
REAL ,INTENT(OUT) :: SIGT2MP(24),SIGHP2(24)
!
LOGICAL :: LSGOBS
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!**---------------------------------------------------------------------
!**  1. Initialisation des polynomes bruts et des champs de reference.
!**     -------------------------------------------------------------
!
! lecture des coefficients polynomiaux
!
IF (LHOOK) CALL DR_HOOK('OI_CAVEGI',0,ZHOOK_HANDLE)
DO J=1,24
  READ(61,'(6X,8F10.4)') VGAT1(J),VGAT2(J),VGAT3(J),&
      VGBT1(J),VGBT2(J),VGBT3(J),&
      VGCT1(J),VGCT2(J)    
ENDDO
!
DO J=1,24
  READ(61,'(6X,8F10.4)') VGAH1(J),VGAH2(J),VGAH3(J),&
      VGBH1(J),VGBH2(J),VGBH3(J),&
      VGCH1(J),VGCH2(J)    
ENDDO
!
DO J = 1,24   
  SIGT2MP(J)= MAX(.3 , 2.7*(1.0 - ((REAL(J)-15.)/9.)**2))   
  SIGHP2(J)= (REAL(J)*2.0/3. - 15.)*1.E-2
ENDDO
!
!**---------------------------------------------------------------------
!**  3. Initialisation des variables internes.
!**     -------------------------------------
!
LSGOBS = SIGT2MO > 0.0 .AND. SIGH2MO > 0.0 .AND. &
    (ABS(SIGH2MO-SIGH2MR) > REPS3 .OR. ABS(SIGT2MO-SIGT2MR) > REPS3)    
IF (LHOOK) CALL DR_HOOK('OI_CAVEGI',1,ZHOOK_HANDLE)
!
!**---------------------------------------------------------------------
END SUBROUTINE OI_CAVEGI


