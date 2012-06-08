C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FANIME_MT (FA,  KREP, KNUMER, KNIMES )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Ce sous-programme permet d'ajuster le Niveau de Messagerie
C     propre aux actions faites sur un fichier particulier, ouvert pour
C     le Logiciel de Fichiers ARPEGE, de meme que le Niveau correspon-
C     dant du logiciel LFI.
C        Cependant, tant que le Niveau de Messagerie Global *FA%NIMSGA*
C     vaut 0 ou 2, le niveau propre au fichier est inoperant.
C     *FA%NIMSGA* vaut par defaut 1, et est reglable via le s/p "FANMSG".
C**
C     Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                 KNUMER (Entree) ==> Numero d'Unite Logique concernee;
C                 KNIMES (Entree) ==> Niveau de Messagerie souhaite.
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KNUMER, KNIMES
C
      INTEGER IREP, IRANG, INIMEX, ILACTI, INIMES
C
      LOGICAL LLRLFI
C
#include "facom2.h"
#include "facom_mt.h"
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FANIME_MT',0,ZHOOK_HANDLE)
      CALL FANUMU_MT (FA, KNUMER,IRANG)
      INIMEX=0
C
      IF (IRANG.EQ.0) THEN
        IREP=-51
      ELSEIF (KNIMES.GE.0.AND.KNIMES.LE.2) THEN
        INIMEX=IXNVMS (IRANG)
        FA%NIVOMS(IRANG)=KNIMES
        CALL LFINIM_MT (FA%LFI, IREP,KNUMER,KNIMES)
        LLRLFI=IREP.NE.0
      ELSE
        IREP=-52
      ENDIF
C
      KREP=IREP
      LLFATA=LLMOER (IREP,IRANG)
C
      IF (LLFATA.OR.MAX (IXNVMS (IRANG),INIMEX).EQ.2) THEN
        INIMES=2
      ELSE
        IF (LHOOK) CALL DR_HOOK('FANIME_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='FANIME'
      WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,
     S       '', KNIMES='',I3)') KREP,KNUMER,KNIMES
      CALL FAIPAR_MT (FA, KNUMER,INIMES,IREP,LLFATA,CLMESS,
     S             CLNSPR,CLACTI,.FALSE.)
C
      IF (LHOOK) CALL DR_HOOK('FANIME_MT',1,ZHOOK_HANDLE)
      END

