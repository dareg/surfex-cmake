C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFIOMF_MT (LFI, KREP, KNUMER, KNIMES )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Ce sous-programme permet d'obtenir la valeur courante du niveau
C     de messagerie propre aux actions faites via le logiciel LFI,
C     sur un fichier particulier ouvert pour ce logiciel.
C
C        Noter que si le niveau global de messagerie *LFI%NIMESG*
C     vaut 0 ou 2, le niveau propre au fichier est inoperant.
C     *LFI%NIMESG* vaut par defaut 1, est reglable via le s/p "LFINMG",
C     et sa valeur courante peut etre obtenue par le s/p "LFIOMG".
C**
C     ARGUMENTS : KREP   (Sortie) ==> Code-REPonse du sous-programme;
C                 KNUMER (Entree) ==> NUMERo d'unite logique concernee;
C                 KNIMES (Sortie) ==> NIveau courant de MESsagerie.
C
#ifndef f77
#include "precision.h"
#endif
C
      TYPE(LFICOM) :: LFI
      INTEGER KREP, KNUMER, KNIMES, IRANG, IREP, INIMES
C
#include "lficom2.h"
#include "lficom_mt.h"
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFIOMF_MT',0,ZHOOK_HANDLE)
      CALL LFINUM_MT (LFI, KNUMER,IRANG)
C
      IF (IRANG.EQ.0) THEN
        IREP=-1
      ELSE
        KNIMES=LFI%NIVMES(IRANG)
        LFI%NDEROP(IRANG)=21
        IREP=0
      ENDIF
C
      LLFATA=LLMOER (IREP,IRANG)
      KREP=IREP
C
      IF (LLFATA.OR.IXNIMS (IRANG).EQ.2) THEN
        INIMES=2
      ELSE
        IF (LHOOK) CALL DR_HOOK('LFIOMF_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='LFIOMF'
      WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,
     S       '', KNIMES='',I3)') KREP,KNUMER,KNIMES
      CALL LFIEMS_MT (LFI, KNUMER,INIMES,IREP,LLFATA,
     S                CLMESS,CLNSPR,CLACTI)
C
      IF (LHOOK) CALL DR_HOOK('LFIOMF_MT',1,ZHOOK_HANDLE)
      END

