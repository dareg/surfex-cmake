C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFIOSF_MT (LFI, KREP, KNUMER, LDIMST )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Ce sous-programme permet d'obtenir l'option courante gouvernant
C     l'impression de STATISTIQUES a la fermeture d'un fichier
C     particulier, ouvert pour le logiciel LFI.
C
C        Noter que si le niveau global d'impression des statistiques
C     *LFI%NISTAG* vaut 0 ou 2, l'option propre au fichier est inoperante.
C     *LFI%NISTAG* vaut par defaut 1, est reglable via le s/p "LFINSG",
C     et sa valeur peut etre obtenue par le s/p "LFIOSG".
C**
C     ARGUMENTS : KREP   (Sortie) ==> Code-REPonse du sous-programme;
C                 KNUMER (Entree) ==> NUMERo d'unite logique concernee;
C                 LDIMST (Sortie) ==> Option d'IMpression des STatisti-
C                                     ques a la fermeture (vrai=oui).
C
#ifndef f77
#include "precision.h"
#endif
C
      TYPE(LFICOM) :: LFI
      INTEGER KREP, KNUMER, IRANG, IREP, INIMES
C
      LOGICAL LDIMST
C
#include "lficom2.h"
#include "lficom_mt.h"
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFIOSF_MT',0,ZHOOK_HANDLE)
      CALL LFINUM_MT (LFI, KNUMER,IRANG)
C
      IF (IRANG.NE.0) THEN
         IF (LFI%LMULTI) CALL LFIVER_MT (LFI, LFI%VERRUE(IRANG),'ON')
        LDIMST=LFI%LISTAT(IRANG)
        LFI%NDEROP(IRANG)=20
         IF (LFI%LMULTI) CALL LFIVER_MT (LFI, LFI%VERRUE(IRANG),'OFF')
        IREP=0
      ELSE
        IREP=-1
      ENDIF
C
      KREP=IREP
      LLFATA=LLMOER (IREP,IRANG)
C
      IF (LLFATA.OR.IXNIMS (IRANG).EQ.2) THEN
        INIMES=2
      ELSE
        IF (LHOOK) CALL DR_HOOK('LFIOSF_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='LFIOSF'
      WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,
     S       '', LDIMST= '',L1)') KREP,KNUMER,LDIMST
      CALL LFIEMS_MT (LFI, KNUMER,INIMES,IREP,LLFATA,
     S                CLMESS,CLNSPR,CLACTI)
C
      IF (LHOOK) CALL DR_HOOK('LFIOSF_MT',1,ZHOOK_HANDLE)
      END

