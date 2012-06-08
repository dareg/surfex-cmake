C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFIOMG_MT (LFI, KNIVAU )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Ce sous-programme permet d'obtenir la valeur courante du NIVEAU
C      GLOBAL d'IMPRESSION des MESSAGES emis par le logiciel LFI.
C**
C        ARGUMENT : KNIVAU (Sortie) ==> Niveau global de messagerie.
C
C     Valeurs possibles (par defaut, 1):
C
C     0 : N'emettre que les messages d'erreurs reellement importants,
C         en pratique ceux relies a une erreur fatale. Le niveau indivi-
C         duel de messagerie des fichiers est alors inoperant.
C
C     1 : Seuls quelques messages "globaux" (c'est-a-dire non reliables
C         a un fichier deja ouvert) et les messages lies a un fichier
C         ouvert, messages dont le niveau est au plus egal au niveau de
C         messagerie individuelle du fichier.
C
C     2 : Emettre tous les messages possibles (donc jusqu'au niveau 2,
C         mais pas ceux lies au mode "mise au point"), meme si
C         ces messages concernent un fichier dont le niveau individuel
C         de messagerie est inferieur a 2.
C
#ifndef f77
#include "precision.h"
#endif
C
      TYPE(LFICOM) :: LFI
      INTEGER KNIVAU, INIMES, IREP, INUMER
#include "lficom_mt.h"
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFIOMG_MT',0,ZHOOK_HANDLE)
      IF (LFI%LFIOMG_LLPREA) THEN
        CALL LFIINI_MT (LFI, 2)
        LFI%LFIOMG_LLPREA=.FALSE.
      ENDIF
C
      KNIVAU=LFI%NIMESG
      INIMES=LFI%NIMESG
C
      IF (INIMES.EQ.2) THEN
        IREP=0
        INUMER=LFI%JPNIL
        CLNSPR='LFIOMG'
        WRITE (UNIT=CLMESS,FMT='(''KNIVAU='',I2)') KNIVAU
        CALL LFIEMS_MT (LFI, INUMER,INIMES,IREP,.FALSE.,CLMESS,
     S                  CLNSPR,CLACTI)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('LFIOMG_MT',1,ZHOOK_HANDLE)
      END

