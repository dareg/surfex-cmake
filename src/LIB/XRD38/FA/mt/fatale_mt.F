C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FATALE_MT (FA,  KREP, KNUMER, LDERFA )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Ce sous-programme permet d'activer ou de desactiver l'option
C     rendant fatale toute erreur detectee sur un fichier particulier,
C     ouvert pour le Logiciel de Fichiers ARPEGE, de meme pour l'option
C     correspondante du logiciel LFI.
C        Cependant, tant que le niveau global d'erreur fatale *FA%NRFAGA*
C     vaut 0 ou 2, l'option propre au fichier est inoperante.
C     *FA%NRFAGA* vaut par defaut 1, et est reglable via le s/p "FANERG".
C**
C     Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                 KNUMER (Entree) ==> Numero d'Unite Logique concernee;
C                 LDERFA (Entree) ==> Option d'Erreur Fatale (Vrai=oui).
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KNUMER
C
      INTEGER IRANG, IREP, INIMES
C
      LOGICAL LDERFA, LLRLFI
C
#include "facom2.h"
#include "facom_mt.h"
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FATALE_MT',0,ZHOOK_HANDLE)
      CALL FANUMU_MT (FA, KNUMER,IRANG)
C
      IF (IRANG.NE.0) THEN
        FA%LERRFA(IRANG)=LDERFA
        CALL LFIERF_MT (FA%LFI, IREP,KNUMER,LDERFA)
        LLRLFI=IREP.NE.0
      ELSE
        IREP=-51
        LLRLFI=.FALSE.
      ENDIF
C
      LLFATA=LLMOER (IREP,IRANG)
      KREP=IREP
C
      IF (LLFATA.OR.IXNVMS (IRANG).EQ.2) THEN
        INIMES=2
      ELSE
        IF (LHOOK) CALL DR_HOOK('FATALE_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='FATALE'
      WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,
     S       '', LDERFA= '',L1)') KREP,KNUMER,LDERFA
      CALL FAIPAR_MT (FA, KNUMER,INIMES,IREP,LLFATA,CLMESS,
     S             CLNSPR,CLACTI,LLRLFI)
C
      IF (LHOOK) CALL DR_HOOK('FATALE_MT',1,ZHOOK_HANDLE)
      END

