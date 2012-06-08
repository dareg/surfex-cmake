C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FACTEC_MT (FA, KREP, PA, KNBIT, KDEC, KE, KNUTIL)
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
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
      IMPLICIT NONE
#include "precision.h"
      TYPE(FA_COM) :: FA
      REAL (KIND=JPDBLR) PA
      INTEGER KREP, KNBIT, KDEC, KE, KNUTIL
!
      REAL (KIND=JPDBLR) ZTWO, ZHALF, ZTEN
#include "facom_mt.h"
!
!**
!     1.  -  CONTROLES ET INITIALISATIONS.
!-----------------------------------------------------------------------
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FACTEC_MT',0,ZHOOK_HANDLE)
      ZTWO   = 2.
      ZHALF  = 0.5
      ZTEN   = 10.
      KREP   = 0
      KE     = 0
      KNUTIL = 0
      IF (KNBIT.LE.0 .OR. KNBIT.GT.64) THEN
        KREP=-1
        WRITE (UNIT=FA%NULOUT,FMT=*)'****'
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S         '**** FACTEC: ERROR, bits number out of range 1-64'
        WRITE (UNIT=FA%NULOUT,FMT=*)'****         KNBIT = ',KNBIT
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S         '****         Binary scale factor is not computed !!'
        WRITE (UNIT=FA%NULOUT,FMT=*)'****'
        IF (LHOOK) CALL DR_HOOK('FACTEC_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
      IF ( ABS(PA).LT.TINY(PA) ) THEN
        WRITE (UNIT=FA%NULOUT,FMT=*)'----'
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S         '---- FACTEC: Warning, the range of the field is',
     S         ' considered as zero'
        WRITE (UNIT=FA%NULOUT,FMT=*)'----'
        KE = 0
        KNUTIL = 1
        IF (LHOOK) CALL DR_HOOK('FACTEC_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
      IF ( ABS(LOG10(ABS(PA))+KDEC) .GE. RANGE(PA) ) THEN
        KREP=-1
        WRITE (UNIT=FA%NULOUT,FMT=*)'****'
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S                '**** FACTEC: ERROR, PA*10**KDEC exceeds real',
     S                'representation of KIND=',JPDBLR
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S         '****         LOG10(ABS(PA)), KDEC, RANGE(PA) = ',
     S         LOG10(ABS(PA)), KDEC, RANGE(PA)
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S         '****         Binary scale factor is not computed !!'
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
      KE = FLOOR( LOG( (PA*10._8**KDEC) /
     S                 (2._8**KNBIT-0.5_8) )/LOG(2._8) ) + 1
! KE = FLOOR( LOG( (PA*ZTEN**KDEC) / (ZTWO**KNBIT-ZHALF) )/LOG(ZTWO) ) + 1
!
! KNUTIL = FLOOR( 0.5 + PA*(10.**KDEC)*(2.**(-KE)) )
      KNUTIL = FLOOR( 0.5_8 + PA*(10._8**KDEC)*(2._8**(-KE)) )
! KNUTIL = FLOOR( ZHALF + PA*(ZTEN**KDEC)*(ZTWO**(-KE)) )
!
      IF (LHOOK) CALL DR_HOOK('FACTEC_MT',1,ZHOOK_HANDLE)
      END

