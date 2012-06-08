C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FALIMU_MT (FA,  KXNIVV, KXTRON, KXLATI, KXLONG )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-programme servant a obtenir les valeurs courantes
C     des LIMites Utilisateur en termes de Resolutions horizontale
C     et verticale, valides globalement pour tous les Fichiers et Cadres
C     ARPEGE manipules par le programme utilisateur.
C**
C        Arguments : KXNIVV ==> Nombre maximum de niveaux verticaux;
C  (tous de Sortie)  KXTRON ==> Troncature maximum;
C                    KXLATI ==> Nombre maximum de latitudes pole a pole;
C                    KXLONG ==> Nombre maxi de longitudes par parallele.
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KXNIVV, KXTRON, KXLATI, KXLONG
C
      INTEGER INUMER, INIMES, IREP
C
      LOGICAL LLVERG
C
C
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  SI PREMIERE UTILISATION, APPEL AU SOUS-PROGRAMME "FARINE".
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FALIMU_MT',0,ZHOOK_HANDLE)
      IF (FA%FALIMU_LLPREA) THEN
        CALL FARINE_MT (FA, 2)
        FA%FALIMU_LLPREA=.FALSE.
      ENDIF
C
C             Verrouillage global, si necessaire.
C
      IF (FA%LFAMUL) CALL LFIVER_MT (FA%LFI, FA%VRGLAS,'ON')
      LLVERG=FA%LFAMUL
C**
C     2.  -  RECOPIE DES VALEURS EN COMMON DANS LES ARGUMENTS.
C-----------------------------------------------------------------------
C
      KXNIVV=FA%NXNIVV
      KXTRON=FA%NXTRON
      KXLATI=FA%NXLATI
      KXLONG=FA%NXLONG
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE sous-programme "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
C
C          Deverrouillage global eventuel.
C
      IF (LLVERG) CALL LFIVER_MT (FA%LFI, FA%VRGLAS,'OFF')
C
      IF (FA%NIMSGA.EQ.2) THEN
        IREP=0
        INIMES=2
        CLNSPR='FALIMU'
        INUMER=FA%JPNIIL
        WRITE (UNIT=CLMESS,FMT='(''KXNIVV='',I4,'', KXTRON='',I4,
     S         '', KXLATI='',I4,'', KXLONG='',I4)')
     S     KXNIVV,KXTRON,KXLATI,KXLONG
        CALL FAIPAR_MT (FA, INUMER,INIMES,IREP,LLFATA,CLMESS,
     S               CLNSPR,CLACTI,.FALSE.)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FALIMU_MT',1,ZHOOK_HANDLE)
      END

