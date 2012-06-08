C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FANDAR_MT (FA,  KREP, KNUMER, KDATEF )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C      Sous-programme de definition d'une (Nouvelle) Date sur un fichier
C     ARpege.
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KNUMER (Entree) ==> Numero de l'unite logique;
C     (Tableau)  KDATEF (Entree) ==> Date elle-meme (FA%JPLDAT mots).
C*
C        En cas de modification effective (si le fichier etait deja muni
C     d'une date), il y a messagerie de niveau 1.
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KNUMER
      INTEGER KDATEF (FA%JPLDAT)
      INTEGER (KIND=JPDBLE) KLDATEF (FA%JPLDAT)
C
      INTEGER IRANG, IREP, INIMES, J
C
      LOGICAL LLVERF, LLRLFI, LLMODA
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  CONTROLES DES PARAMETRES D'APPEL, PUIS INITIALISATIONS.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FANDAR_MT',0,ZHOOK_HANDLE)
      LLVERF=.FALSE.
      LLRLFI=.FALSE.
      LLMODA=.FALSE.
      CALL FANUMU_MT (FA, KNUMER,IRANG)
C
      IF (IRANG.EQ.0) THEN
        IREP=-51
        GOTO 1001
      ENDIF
C
C         Verrouillage eventuel du fichier.
C
      IF (FA%LFAMUL) CALL LFIVER_MT (FA%LFI, FA%VRFICH(IRANG),'ON')
      LLVERF=FA%LFAMUL
C**
C     2.  -  DEFINITION PROPREMENT DITE VIA LE SOUS-PROGRAMME "FANDAI".
C            ( controles, puis mise a jour de FA%MADATE(.,IRANG) )
C-----------------------------------------------------------------------
C
      CALL FANDAI_MT (FA, IREP,IRANG,KDATEF,LLMODA)
C
      IF (IREP.EQ.0) THEN
C**
C     3.  -  ECRITURE DE LA DATE SUR LE FICHIER.
C-----------------------------------------------------------------------
C
        KLDATEF=KDATEF
        CALL LFIECR_MT (FA%LFI, IREP,KNUMER,FA%CPDATE,KLDATEF,FA%JPLDAT)
        LLRLFI=IREP.NE.0
        FA%LCREAF(IRANG)=FA%LCREAF(IRANG).AND.LLRLFI
      ENDIF
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE SOUS-PROGRAMME "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
      KREP=IREP
      LLFATA=LLMOER (IREP,IRANG)
C
C        Deverrouillage eventuel du fichier.
C
      IF (LLVERF) CALL LFIVER_MT (FA%LFI, FA%VRFICH(IRANG),'OFF')
C
      IF (LLFATA) THEN
        INIMES=2
      ELSEIF (IREP.NE.0) THEN
        INIMES=0
      ELSE
        INIMES=IXNVMS(IRANG)
      ENDIF
C
      IF (.NOT.LLFATA.AND.INIMES.EQ.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('FANDAR_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='FANDAR'
C
      IF (INIMES.GE.1.AND.LLMODA) THEN
        WRITE (UNIT=CLMESS,FMT=
     S         '(''MODIFICATION DE LA DATE, UNITE'',I3)') KNUMER
        CALL FAIPAR_MT (FA, KNUMER,INIMES,IREP,.FALSE.,CLMESS,
     S                  CLNSPR,CLACTI,.FALSE.)
      ENDIF
C
      IF (INIMES.EQ.2) THEN
C***** FAZZZZ - KREP=iiii, KNUMER=iii, KDATEF(1:5)=iiiii/ii/ii iii:ii, *****
C*****          KDATEF(7:8)=iiiiii-iiiiii                              *****
        WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,
     S       '', KDATEF(1:5)='',I5,2(''/'',I2),I3,'':'',I2.2,
     S       '', KDATEF(7:8)='',I6,''-'',I6)') KREP,KNUMER,
     S     (KDATEF(J),J=1,5),(KDATEF(J),J=7,8)
        CALL FAIPAR_MT (FA, KNUMER,INIMES,IREP,LLFATA,CLMESS,
     S               CLNSPR,CLACTI,LLRLFI)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FANDAR_MT',1,ZHOOK_HANDLE)
      END

