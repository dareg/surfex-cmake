C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFIOEF_MT (LFI, KREP, KNUMER, LDERFA )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Ce sous-programme permet d'obtenir l'option courante
C     de traitement des ERREURS detectees par le logiciel LFI,
C     pour un fichier particulier ouvert pour ce logiciel.
C
C        Noter que si le niveau global d'erreur fatale *LFI%NERFAG*
C     vaut 0 ou 2, l'option propre au fichier est inoperante.
C     *LFI%NERFAG* vaut par defaut 1, est reglable via le s/p "LFINEG",
C     et sa valeur peut etre obtenue par le s/p "LFIOEG".
C**
C     ARGUMENTS : KREP   (Sortie) ==> Code-REPonse du sous-programme;
C                 KNUMER (Entree) ==> NUMERo d'unite logique concernee;
C                 LDERFA (Sortie) ==> Option d'ERreur FAtale (vrai=oui).
C
#ifndef f77
#include "precision.h"
#endif
C
      TYPE(LFICOM) :: LFI
      INTEGER KREP, KNUMER, IRANG, IREP, INIMES
C
      LOGICAL LDERFA
C
#include "lficom2.h"
#include "lficom_mt.h"
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFIOEF_MT',0,ZHOOK_HANDLE)
      CALL LFINUM_MT (LFI, KNUMER,IRANG)
C
      IF (IRANG.EQ.0) THEN
        IREP=-1
      ELSE
        LDERFA=LFI%LERFAT(IRANG)
        LFI%NDEROP(IRANG)=19
        IREP=0
      ENDIF
C
      LLFATA=LLMOER (IREP,IRANG)
      KREP=IREP
C
      IF (LLFATA.OR.IXNIMS (IRANG).EQ.2) THEN
        INIMES=2
      ELSE
        IF (LHOOK) CALL DR_HOOK('LFIOEF_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='LFIOEF'
      WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,
     S       '', LDERFA= '',L1)') KREP,KNUMER,LDERFA
      CALL LFIEMS_MT (LFI, KNUMER,INIMES,IREP,LLFATA,
     S                CLMESS,CLNSPR,CLACTI)
C
      IF (LHOOK) CALL DR_HOOK('LFIOEF_MT',1,ZHOOK_HANDLE)
      END

