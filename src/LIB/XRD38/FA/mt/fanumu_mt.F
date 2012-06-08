C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FANUMU_MT (FA,  KNUMER, KRANG )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Ce sous-programme calcule le RANG du Numero d'Unite logique
C     *KNUMER* dans la table des unites logiques *FA%NULOGI*;
C     s'il n'y est pas trouve, le resultat est ZERO.
C        Ce sous-programme, appele par tous les sous-programmes
C     "fichier" du Logiciel de Fichiers ARPEGE, se charge lors de son
C     premier appel d'appeler le sous-programme preparatoire FARINE.
C**
C       Arguments : KNUMER (Entree) ==> Numero d'unite logique cherche;
C                   KRANG  (Sortie) ==> Rang dans la table des fichiers
C                                       du logiciel FA (0 si absent).
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KNUMER, KRANG
C
      INTEGER J, IRESUL
C
C
#include "facom2.h"
#include "facom_mt.h"
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FANUMU_MT',0,ZHOOK_HANDLE)
      IF (FA%FANUMU_LLPREA) THEN
        CALL FARINE_MT (FA, 2)
        FA%FANUMU_LLPREA=.FALSE.
      ENDIF
C
C          VERROUILLAGE GLOBAL (A CAUSE DE L'UTILISATION DE FA%NFIOUV )
C
      IF (FA%LFAMUL) CALL LFIVER_MT (FA%LFI, FA%VRGLAS,'ON')
C
      DO 10 J=1,FA%NFIOUV
C
      IF (KNUMER.EQ.FA%NULOGI(FA%NULIND(J))) THEN
        IRESUL=FA%NULIND(J)
        GOTO 20
      ENDIF
C
   10 CONTINUE
C
      IRESUL=0
C
   20 CONTINUE
C
C          DEVERROUILLAGE GLOBAL
C
      IF (FA%LFAMUL) CALL LFIVER_MT (FA%LFI, FA%VRGLAS,'OFF')
      KRANG=IRESUL
C
      IF (LHOOK) CALL DR_HOOK('FANUMU_MT',1,ZHOOK_HANDLE)
      END

