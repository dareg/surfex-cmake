C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FAIRME_MT (FA,  KREP, KNUMER, CDSTTU )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C       Sous-programme de FERMETURE d'une unite logique "Fichier ARPEGE"
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KNUMER (Entree) ==> Numero de l'unite logique;
C                CDSTTU (Entree) ==> "STATUS" eventuel pour "CLOSE".
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KNUMER
C
      INTEGER IREP, IRANG, J, IPOSNU, IRANGC, INIMES, ILNOMC
C
      CHARACTER CDSTTU*(*), CLSTTU*7
C
      LOGICAL LLSTTU, LLVERF, LLRLFI
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  CONTROLES DES PARAMETRES D'APPEL, PUIS INITIALISATIONS.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FAIRME_MT',0,ZHOOK_HANDLE)
      IREP=0
      LLVERF=.FALSE.
      LLRLFI=.FALSE.
      CALL FANUMU_MT (FA, KNUMER,IRANG)
C
C         Verrouillage global eventuel.
C
      IF (FA%LFAMUL) CALL LFIVER_MT (FA%LFI, FA%VRGLAS,'ON')
C
      IF (IRANG.EQ.0) THEN
        IREP=-51
        GOTO 1001
      ELSEIF (LEN (CDSTTU).LE.0) THEN
        IREP=-65
        GOTO 1001
      ELSE
        LLSTTU=CDSTTU.EQ.'KEEP'.OR.CDSTTU.EQ.'DELETE'
C
        IF (LLSTTU) THEN
          CLSTTU=CDSTTU(1:MIN (LEN (CDSTTU),LEN (CLSTTU)))
        ELSE
          CLSTTU='DEFAUT'
        ENDIF
C
      ENDIF
C
      IF (FA%LFAMUL) CALL LFIVER_MT (FA%LFI, FA%VRFICH(IRANG),'ON')
      LLVERF=FA%LFAMUL
C
      IF (FA%LCREAF(IRANG).AND..NOT.LLSTTU) THEN
C
C         On force le relachement d'un fichier "parasite".
C
        LLSTTU=.TRUE.
        CLSTTU='DELETE'
      ENDIF
C**
C     2.  -  FERMETURE DU FICHIER, AU SENS DU LOGICIEL LFI.
C-----------------------------------------------------------------------
C
      CALL LFIFER_MT (FA%LFI, IREP,KNUMER,CLSTTU)
C
      IF (IREP.NE.0) THEN
        LLRLFI=.TRUE.
        GOTO 1001
      ENDIF
C**
C     3.  -  "NETTOYAGE" DES TABLES AYANT PERMIS DE GERER LE FICHIER.
C            ( au moins celles ayant un caractere "global" )
C-----------------------------------------------------------------------
C
      FA%NULOGI(IRANG)=FA%JPNIIL
C
      DO 301 J=1,FA%NFIOUV
C
      IF (FA%NULIND(J).EQ.IRANG) THEN
        IPOSNU=J
        GOTO 302
      ENDIF
C
  301 CONTINUE
C
      IREP=-66
      GOTO 1001
C
  302 CONTINUE
C
      FA%NFIOUV=FA%NFIOUV-1
C
      DO 303 J=IPOSNU,FA%NFIOUV
      FA%NULIND(J)=FA%NULIND(J+1)
  303 CONTINUE
C
      IF (FA%LFAMUL) THEN
        CALL LFIVER_MT (FA%LFI, FA%VRFICH(IRANG),'OFF')
        CALL LFIVER_MT (FA%LFI, FA%VRFICH(IRANG),'REL')
      ENDIF
C
      LLVERF=.FALSE.
      IRANGC=FA%NUCADR(IRANG)
      FA%NULCAD(IRANGC)=FA%NULCAD(IRANGC)-1
C
C        Si le cadre auquel etait rattache le fichier n'a plus d'autre
C     fichier rattache, et qu'on ne devait pas conserver ce cadre,
C     on le supprime.
C
      IF (FA%NULCAD(IRANGC).LE.0.AND.
     S   (FA%NGARDE(IRANGC).EQ.0.OR.(FA%NGARDE(IRANGC).EQ.1.AND.
     S   .NOT.FA%LIGARD)))
     S      CALL FACTUI_MT (FA, IREP,IRANGC)
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE SOUS-PROGRAMME "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
      KREP=IREP
      LLFATA=LLMOER (IREP,IRANG)
C
C        Deverrouillage(s) eventuel(s).
C
      IF (LLVERF) CALL LFIVER_MT (FA%LFI, FA%VRFICH(IRANG),'OFF')
      IF (FA%LFAMUL) CALL LFIVER_MT (FA%LFI, FA%VRGLAS,'OFF')
C
      IF (LLFATA) THEN
        INIMES=2
      ELSE
        INIMES=IXNVMS(IRANG)
      ENDIF
C
      IF (INIMES.EQ.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('FAIRME_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='FAIRME'
C
      IF (IREP.EQ.-65) THEN
        ILNOMC=8
        CLACTI(1:ILNOMC)=FA%CHAINC(:ILNOMC)
       ELSE
        ILNOMC=MIN ( LEN (CDSTTU), LEN (CLACTI) )
        CLACTI(1:ILNOMC)=CDSTTU(1:ILNOMC)
      ENDIF
C
      IF (INIMES.EQ.2) THEN
C
        ILNOMC=MIN (ILNOMC,FA%NCPCAD)
        WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,
     S         '', CDSTTU='''''',A,'''''''')') KREP,KNUMER,
     S       CLACTI(1:ILNOMC)
        CALL FAIPAR_MT (FA, KNUMER,INIMES,IREP,LLFATA,
     S               CLMESS,CLNSPR,
     S               CLACTI(1:ILNOMC),LLRLFI)
C
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FAIRME_MT',1,ZHOOK_HANDLE)
      END

