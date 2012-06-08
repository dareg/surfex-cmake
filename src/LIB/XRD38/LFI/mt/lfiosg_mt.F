C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFIOSG_MT (LFI, KNIVAU )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Ce sous-programme permet d'obtenir la valeur courante du NIVEAU
C     GLOBAL d'IMPRESSION des STATISTIQUES concernant les fichiers
C     manipules par le logiciel LFI, lors de leur fermeture.
C**
C        ARGUMENT : KNIVAU (Sortie) ==> Niveau global d'impression des
C                                       statistiques a la fermeture.
C
C     Valeurs possibles (par defaut, 1):
C
C     0 : Dans ce cas, on n'imprime pas de statistiques lors de la fer-
C         meture d'un fichier, meme si l'option individuelle du fichier
C         concerne est a .TRUE. .
C
C     1 : Respect de l'option individuelle du fichier.
C
C     2 : Dans ce cas, force l'impression des statistiques a la ferme-
C         ture d'un fichier, meme si l'option individuelle du fichier
C         est a .FALSE. .
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
      IF (LHOOK) CALL DR_HOOK('LFIOSG_MT',0,ZHOOK_HANDLE)
      IF (LFI%LFIOSG_LLPREA) THEN
        CALL LFIINI_MT (LFI, 2)
        LFI%LFIOSG_LLPREA=.FALSE.
      ENDIF
C
      KNIVAU=LFI%NISTAG
      INIMES=LFI%NIMESG
C
      IF (INIMES.EQ.2) THEN
        IREP=0
        INUMER=LFI%JPNIL
        CLNSPR='LFIOSG'
        WRITE (UNIT=CLMESS,FMT='(''KNIVAU='',I2)') KNIVAU
        CALL LFIEMS_MT (LFI, INUMER,INIMES,IREP,.FALSE.,
     S                  CLMESS,CLNSPR,CLACTI)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('LFIOSG_MT',1,ZHOOK_HANDLE)
      END

