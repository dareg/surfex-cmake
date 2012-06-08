C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FANUCA_MT (FA,  CDNOMC, KRANGC, LDVERR )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Ce sous-programme calcule le RANG du Cadre de NOM
C     *CDNOMC* dans la table des noms de cadres *FA%CNOMCA*;
C     S'IL N'Y EST PAS TROUVE, LE RESULTAT EST ZERO.
C        Ce sous-programme, appele par plusieurs sous-programmes
C     du Logiciel de Fichiers ARPEGE, se charge lors de son
C     premier appel d'appeler le sous-programme preparatoire FARINE.
C**
C       Arguments : CDNOMC (Entree) ==> Nom de cadre cherche;
C                   KRANGC (Sortie) ==> Rang dans la table des cadres
C                                       du logiciel FA (0 si absent);
C                   LDVERR (Entree) ==> Sert, en mode multi-taches
C                                       seulement, a savoir si l'on doit
C                                       verrouiller globalement ou pas.
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KRANGC
C
      INTEGER J, IRESUL
C
      CHARACTER CDNOMC*(*)
C
      LOGICAL LDVERR, LLVERG
C
C
#include "facom2.h"
#include "facom_mt.h"
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FANUCA_MT',0,ZHOOK_HANDLE)
      IF (FA%FANUCA_LLPREA) THEN
        CALL FARINE_MT (FA, 2)
        FA%FANUCA_LLPREA=.FALSE.
      ENDIF
C
C          VERROUILLAGE GLOBAL (A CAUSE DE L'UTILISATION DE FA%NFIOUV )
C
      LLVERG=FA%LFAMUL.AND.LDVERR
      IF (LLVERG) CALL LFIVER_MT (FA%LFI, FA%VRGLAS,'ON')
C
      DO 10 J=1,FA%NCADEF
C
      IF (CDNOMC.EQ.FA%CNOMCA(FA%NCAIND(J))) THEN
        IRESUL=FA%NCAIND(J)
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
      IF (LLVERG) CALL LFIVER_MT (FA%LFI, FA%VRGLAS,'OFF')
      KRANGC=IRESUL
C
      IF (LHOOK) CALL DR_HOOK('FANUCA_MT',1,ZHOOK_HANDLE)
      END

