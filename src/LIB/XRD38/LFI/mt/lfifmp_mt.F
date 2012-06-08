C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFIFMP_MT (LFI, KNUMER, KRANFM )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-programme *INTERNE* du Logiciel de Fichiers Indexes LFI.
C     Calcule le rang de l'unite logique *KNUMER* dans la table des
C     Unites Logiques a Facteur Multiplicatif predefini;
C     si l'Unite Logique n'y est pas trouvee, le resultat est ZERO.
C
C        En mode Multi-Taches, il est necessaire de verrouiller
C     Globalement le code faisant appel a ce sous-programme.
C**
C       ARGUMENTS : KNUMER (Entree) ==> Numero d'unite logique cherche;
C                   KRANFM (Sortie) ==> Rang dans la table des unites
C                                       logiques a Facteur Multiplicatif
C                                       predefini (0 si pas trouve).
C
#ifndef f77
#include "precision.h"
#endif
C
      TYPE(LFICOM) :: LFI
      INTEGER KNUMER, KRANFM, J, IRANFM
#include "lficom_mt.h"
C**
C     1.  -  RECHERCHE DIRECTE DANS LA TABLE *LFI%MULOFM*.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFIFMP_MT',0,ZHOOK_HANDLE)
      DO 101 J=1,LFI%NULOFM
C
      IF (KNUMER.EQ.LFI%MULOFM(J)) THEN
        IRANFM=J
        GOTO 102
      ENDIF
C
  101 CONTINUE
C
      IRANFM=0
C
  102 CONTINUE
      KRANFM=IRANFM
C
      IF (LHOOK) CALL DR_HOOK('LFIFMP_MT',1,ZHOOK_HANDLE)
      END

