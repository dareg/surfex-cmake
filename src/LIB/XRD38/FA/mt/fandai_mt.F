C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FANDAI_MT (FA,  KREP, KRANG, KDATEF, LDMODA )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C      Sous-programme INTERNE du logiciel de Fichiers ARPEGE:
C     Definition d'une (Nouvelle) Date.
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KRANG  (Entree) ==> Rang de l'unite logique;
C     (Tableau)  KDATEF (Entree) ==> Date elle-meme (FA%JPLDAT mots).
C                LDMODA (Sortie) ==> Vrai s'il y a modification d'une
C                                    date deja definie.
C*
C       En mode multi-taches, il doit y avoir verrouillage du fichier
C     concerne avant l'appel au sous-programme.
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KRANG
      INTEGER KDATEF (FA%JPLDAT)
C
      INTEGER IMI123, IMAX69, IMINIM, J, ILMOIS, IDEBUT, INIMES, INUMER
C
      LOGICAL LDMODA
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  CONTROLES DES PARAMETRES D'APPEL.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FANDAI_MT',0,ZHOOK_HANDLE)
      LDMODA=.FALSE.
C
      IF (KRANG.LE.0.OR.KRANG.GT.FA%JPNXFA) THEN
        KREP=-66
        GOTO 1001
      ENDIF
C
C         Controle de la Date proprement dite.
C
      IMI123=MIN (KDATEF(1),KDATEF(2),KDATEF(3))
      IMAX69=MAX (KDATEF(6),KDATEF(9))
      IMINIM=KDATEF(1)
C
      DO 103 J=2,FA%JPLDAT
      IMINIM=MIN (IMINIM,KDATEF(J))
  103 CONTINUE
C
      IF (IMINIM.LT.0.OR.IMI123.LE.0.OR.KDATEF(2).GT.12.OR.
     S    KDATEF(3).GT.31.OR.KDATEF(4).GT.24.OR.KDATEF(5).GE.60.OR.
     S    IMAX69.GE.255.OR.
     S (KDATEF(10).LE.KDATEF(11).AND.(KDATEF(10)*KDATEF(11)).NE.0)) THEN
C
C        Erreur de syntaxe.
C
        KREP=-82
        GOTO 1001
      ELSEIF ((KDATEF(2).GT.7.OR.MOD (KDATEF(2),2).EQ.0).AND.
     S        (KDATEF(2).LE.7.OR.MOD (KDATEF(2),2).EQ.1)) THEN
C
C        Controle de coherence (annee,mois,jour).
C
        IF (KDATEF(2).EQ.2) THEN
          ILMOIS=28+MAX (0,1-MOD (KDATEF(1),4))
        ELSE
          ILMOIS=30
        ENDIF
C
        IF (KDATEF(3).GT.ILMOIS) THEN
          KREP=-82
          GOTO 1001
        ENDIF
C
      ENDIF
C
      KREP=0
C**
C     2.  -  SI DATE DEJA DEFINIE, COMPARAISON ANCIENNE/NOUVELLE.
C-----------------------------------------------------------------------
C
      IF (FA%LCREAF(KRANG)) THEN
        IDEBUT=1
      ELSE
C
        DO 201 J=1,FA%JPLDAT
C
        IF (FA%MADATE(J,KRANG).NE.KDATEF(J)) THEN
          LDMODA=.TRUE.
          IDEBUT=J
          GOTO 300
        ENDIF
C
  201   CONTINUE
C
C         Si on arrive ici, il y a redefinition a l'identique.
C
        GOTO 1001
      ENDIF
C**
C     3.  -  SI NECESSAIRE, MISE A JOUR DU TABLEAU "FA%MADATE".
C-----------------------------------------------------------------------
C
  300 CONTINUE
C
      DO 301 J=IDEBUT,FA%JPLDAT
      FA%MADATE(J,KRANG)=KDATEF(J)
  301 CONTINUE
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE EVENTUELLE,
C            VIA LE SOUS-PROGRAMME "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
      LLFATA=LLMOER (KREP,KRANG)
C
      IF (FA%LFAMOP.OR.LLFATA) THEN
        INIMES=2
        CLNSPR='FANDAI'
        INUMER=FA%JPNIIL
        WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KRANG='',I4,
     S       '', KDATEF(1:5)='',I5,2(''/'',I2),I3,'':'',I2.2,
     S       '', KDATEF(7:8)='',I6,''-'',I6,'', LDMODA= '',L1)')
     S     KREP,KRANG,(KDATEF(J),J=1,5),(KDATEF(J),J=7,8),LDMODA
        CALL FAIPAR_MT (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS,
     S                  CLNSPR,CLACTI, .FALSE.)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FANDAI_MT',1,ZHOOK_HANDLE)
      END

