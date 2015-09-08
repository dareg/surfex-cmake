! Oct-2012 P. Marguinaud 64b LFI
! Jan-2011 P. Marguinaud Thread-safe FA
SUBROUTINE FACTEC_MT64                                      &
&                     (FA, KREP, PA, KNBIT, KDEC, KE, KNUTIL)
USE FA_MOD, ONLY : FA_COM
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION
IMPLICIT NONE
!****
!      Sous-programme de calcul du FACTeur d'EChelle binaire associe
!     a un champ d'amplitude donnee, code sur KNBIT bits.
!**
!    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
!                PA     (Entree) ==> Amplitude du champ a compacter;
!                KNBIT  (Entree) ==> Nb de bits servant au compactage;
!                KDEC   (Entree) ==> Facteur d'echelle decimal;
!                KE     (Sortie) ==> Facteur d'echelle binaire;
!                KNUTIL (Sortie) ==> Nombre d'entiers utilises pour
!                                    representer le champ.
!
!     Modifications:
!
TYPE(FA_COM) :: FA
REAL (KIND=JPDBLR) PA
INTEGER (KIND=JPLIKB) KREP, KNBIT, KDEC, KE, KNUTIL
!
REAL (KIND=JPDBLR) ZTWO, ZHALF, ZTEN

!
!**
!     1.  -  CONTROLES ET INITIALISATIONS.
!-----------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('FACTEC_MT',0,ZHOOK_HANDLE)
ZTWO   = 2._JPDBLR
ZHALF  = 0.5_JPDBLR
ZTEN   = 10._JPDBLR
KREP   = 0
KE     = 0
KNUTIL = 0
IF (KNBIT.LE.0 .OR. KNBIT.GT.64) THEN
  KREP=-1
  WRITE (UNIT=FA%NULOUT,FMT=*)'****'
  WRITE (UNIT=FA%NULOUT,FMT=*)                               &
&         '**** FACTEC: ERROR, bits number out of range 1-64'
  WRITE (UNIT=FA%NULOUT,FMT=*)'****         KNBIT = ',KNBIT
  WRITE (UNIT=FA%NULOUT,FMT=*)                                 &
&         '****         Binary scale factor is not computed !!'
  WRITE (UNIT=FA%NULOUT,FMT=*)'****'
  IF (LHOOK) CALL DR_HOOK('FACTEC_MT',1,ZHOOK_HANDLE)
  RETURN
ENDIF
IF ( ABS(PA).LT.TINY(PA) ) THEN
  WRITE (UNIT=FA%NULOUT,FMT=*)'----'
  WRITE (UNIT=FA%NULOUT,FMT=*)                              &
&         '---- FACTEC: Warning, the range of the field is', &
&         ' considered as zero'
  WRITE (UNIT=FA%NULOUT,FMT=*)'----'
  KE = 0
  KNUTIL = 1
  IF (LHOOK) CALL DR_HOOK('FACTEC_MT',1,ZHOOK_HANDLE)
  RETURN
ENDIF
IF ( ABS(LOG10(ABS(PA))+REAL (KDEC, JPDBLR))  &
&     .GE. REAL (RANGE(PA), JPDBLR) ) THEN
  KREP=-1
  WRITE (UNIT=FA%NULOUT,FMT=*)'****'
  WRITE (UNIT=FA%NULOUT,FMT=*)                                  &
&                '**** FACTEC: ERROR, PA*10**KDEC exceeds real', &
&                'representation of KIND=',JPDBLR
  WRITE (UNIT=FA%NULOUT,FMT=*)                              &
&         '****         LOG10(ABS(PA)), KDEC, RANGE(PA) = ', &
&         LOG10(ABS(PA)), KDEC, RANGE(PA)
  WRITE (UNIT=FA%NULOUT,FMT=*)                                 &
&         '****         Binary scale factor is not computed !!'
  WRITE (UNIT=FA%NULOUT,FMT=*)'****'
  IF (LHOOK) CALL DR_HOOK('FACTEC_MT',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!
!**
!     2.  -  CALCUL DU FACTEUR D'ECHELLE BINAIRE
!-----------------------------------------------------------------------
!
! KE = FLOOR( LOG( (PA*10.**KDEC) / (2.**KNBIT-0.5) )/LOG(2.) ) + 1
KE = FLOOR( LOG( (PA*10._8**KDEC) /                &
&                 (2._8**KNBIT-0.5_8) )/LOG(2._8),  &
&                 KIND=JPLIKB ) + 1
! KE = FLOOR( LOG( (PA*ZTEN**KDEC) / (ZTWO**KNBIT-ZHALF) )/LOG(ZTWO) ) + 1
!
! KNUTIL = FLOOR( 0.5 + PA*(10.**KDEC)*(2.**(-KE)) )
KNUTIL = FLOOR( 0.5_8 + PA*(10._8**KDEC)*(2._8**(-KE)), &
&                KIND=JPLIKB )
! KNUTIL = FLOOR( ZHALF + PA*(ZTEN**KDEC)*(ZTWO**(-KE)) )
!
IF (LHOOK) CALL DR_HOOK('FACTEC_MT',1,ZHOOK_HANDLE)
END SUBROUTINE



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE FACTEC64                           &
&           (KREP, PA, KNBIT, KDEC, KE, KNUTIL)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKB)  KREP                                   !   OUT
REAL (KIND=JPDBLR)     PA                                     ! IN   
INTEGER (KIND=JPLIKB)  KNBIT                                  ! IN   
INTEGER (KIND=JPLIKB)  KDEC                                   ! IN   
INTEGER (KIND=JPLIKB)  KE                                     !   OUT
INTEGER (KIND=JPLIKB)  KNUTIL                                 !   OUT

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FACTEC_MT64                                  &
&           (FA, KREP, PA, KNBIT, KDEC, KE, KNUTIL)

END SUBROUTINE

SUBROUTINE FACTEC                             &
&           (KREP, PA, KNBIT, KDEC, KE, KNUTIL)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
REAL (KIND=JPDBLR)     PA                                     ! IN   
INTEGER (KIND=JPLIKM)  KNBIT                                  ! IN   
INTEGER (KIND=JPLIKM)  KDEC                                   ! IN   
INTEGER (KIND=JPLIKM)  KE                                     !   OUT
INTEGER (KIND=JPLIKM)  KNUTIL                                 !   OUT

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FACTEC_MT                                    &
&           (FA, KREP, PA, KNBIT, KDEC, KE, KNUTIL)

END SUBROUTINE

SUBROUTINE FACTEC_MT                              &
&           (FA, KREP, PA, KNBIT, KDEC, KE, KNUTIL)
USE FA_MOD, ONLY : FA_COM
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
TYPE (FA_COM)          FA                                     ! INOUT
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
REAL (KIND=JPDBLR)     PA                                     ! IN   
INTEGER (KIND=JPLIKM)  KNBIT                                  ! IN   
INTEGER (KIND=JPLIKM)  KDEC                                   ! IN   
INTEGER (KIND=JPLIKM)  KE                                     !   OUT
INTEGER (KIND=JPLIKM)  KNUTIL                                 !   OUT
! Local integers
INTEGER (KIND=JPLIKB)  IREP                                   !   OUT
INTEGER (KIND=JPLIKB)  INBIT                                  ! IN   
INTEGER (KIND=JPLIKB)  IDEC                                   ! IN   
INTEGER (KIND=JPLIKB)  IE                                     !   OUT
INTEGER (KIND=JPLIKB)  INUTIL                                 !   OUT
! Convert arguments

INBIT      = INT (     KNBIT, JPLIKB)
IDEC       = INT (      KDEC, JPLIKB)

CALL FACTEC_MT64                                  &
&           (FA, IREP, PA, INBIT, IDEC, IE, INUTIL)

KREP       = INT (      IREP, JPLIKM)
KE         = INT (        IE, JPLIKM)
KNUTIL     = INT (    INUTIL, JPLIKM)

END SUBROUTINE

!INTF KREP            OUT 
!INTF PA            IN    
!INTF KNBIT         IN    
!INTF KDEC          IN    
!INTF KE              OUT 
!INTF KNUTIL          OUT 
