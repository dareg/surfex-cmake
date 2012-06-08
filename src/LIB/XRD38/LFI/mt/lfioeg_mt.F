C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFIOEG_MT (LFI, KNIVAU )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Ce sous-programme permet d'obtenir la valeur courante du NIVEAU
C      GLOBAL de traitement des ERREURS detectees par le logiciel LFI.
C**
C        ARGUMENT : KNIVAU (Sortie) ==> Niveau global de traitement des
C                                       erreurs detectees.
C
C     Valeurs possibles (par defaut, 1):
C
C     0 : Dans ce cas, toute erreur detectee sera fatale, meme si
C         l'erreur joue sur un fichier dont l'option individuelle de
C         traitement des erreurs est a .FALSE. .
C
C     1 : Seules les erreurs "globales" (c'est-a-dire non reliables a un
C         fichier deja ouvert) et les erreurs reliees a un fichier dont
C         l'option individuelle de traitement est .TRUE. seront fatales.
C
C     2 : Passer outre toute erreur detectee, meme si elle correspond
C         a un fichier dont l'option individuelle de traitement est
C         .TRUE. ; noter que dans ce dernier cas le code-reponse ne
C         sera pas nul. Par ailleurs, le code-reponse (-16) echappe a
C         ce mode de controle et est toujours fatal.
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
      IF (LHOOK) CALL DR_HOOK('LFIOEG_MT',0,ZHOOK_HANDLE)
      IF (LFI%LFIOEG_LLPREA) THEN
        CALL LFIINI_MT (LFI, 2)
        LFI%LFIOEG_LLPREA=.FALSE.
      ENDIF
C
      KNIVAU=LFI%NERFAG
      INIMES=LFI%NIMESG
C
      IF (INIMES.EQ.2) THEN
        IREP=0
        INUMER=LFI%JPNIL
        CLNSPR='LFIOEG'
        WRITE (UNIT=CLMESS,FMT='(''KNIVAU='',I2)') KNIVAU
        CALL LFIEMS_MT (LFI, INUMER,INIMES,IREP,.FALSE.,
     S                  CLMESS,CLNSPR,CLACTI)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('LFIOEG_MT',1,ZHOOK_HANDLE)
      END

