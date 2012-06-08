C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFINFO_MT (LFI, KREP, KNUMER, CDNOMA, KLONG, KPOSEX )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        SOUS-PROGRAMME CHARGE DE RENSEIGNER SUR EXISTENCE ET CARACTERI-
C     STIQUES ( LONGUEUR, POSITION ) D'UN ARTICLE LOGIQUE, POUR UNE
C     UNITE LOGIQUE OUVERTE PAR LE LOGICIEL DE FICHIERS INDEXES *LFI*.
C**
C    ARGUMENTS : KREP   (SORTIE) ==> CODE-REPONSE DU SOUS-PROGRAMME;
C                KNUMER (ENTREE) ==> LFI%NUMERO DE L'UNITE LOGIQUE;
C                CDNOMA (ENTREE) ==> NOM DE L'ARTICLE A CHERCHER;
C                KLONG  (SORTIE) ==> LONGUEUR DE L'ARTICLE;
C                KPOSEX (SORTIE) ==> POSITION ( DANS LE FICHIER, DU PRE-
C                                    MIER MOT ) DE L'ARTICLE SUIVANT.
C
C       Noter que si l'unite logique est ouverte pour le logiciel LFI et
C     que l'article demande n'y est pas trouve, KREP, KLONG et KPOSEX
C     sont retournes a ZERO.
C
#ifndef f77
#include "precision.h"
#endif
C
      TYPE(LFICOM) :: LFI
      CHARACTER CDNOMA*(*), CLNOMA*(LFI%JPNCPN)
C
      INTEGER KREP, KNUMER, KLONG, KPOSEX, IREP, IRANG, ILCLNO, ILCDNO
      INTEGER IDECBL, IPOSBL, IARTEX, IRGPIM, INIMES, INBALO, IRETIN
C
      LOGICAL LLVERF
C
#include "lficom2.h"
#include "lficom_mt.h"
C**
C     1.  -  CONTROLES DES PARAMETRES D'APPEL, PUIS INITIALISATIONS.
C-----------------------------------------------------------------------
C
C        Appel legerement anticipe a LFINUM, garantissant l'initialisa-
C     tion des variables globales du logiciel a la 1ere utilisation.
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFINFO_MT',0,ZHOOK_HANDLE)
      CALL LFINUM_MT (LFI, KNUMER,IRANG)
      LLVERF=.FALSE.
      ILCDNO=LEN (CDNOMA)
      KLONG=0
      KPOSEX=0
C
      IF (ILCDNO.LE.0) THEN
        IREP=-15
        CLNOMA=LFI%CHINCO(:LFI%JPNCPN)
        ILCLNO=LFI%JPNCPN
        GOTO 1001
      ELSEIF (CDNOMA.EQ.' ') THEN
        IREP=-18
        CLNOMA=' '
        ILCLNO=1
        GOTO 1001
      ENDIF
C
C        Recherche de la longueur "utile" du nom d'article specifie.
C        (c'est-a-dire sans tenir compte des blancs terminaux eventuels)
C
      IDECBL=0
C
  101 CONTINUE
      IPOSBL=IDECBL+INDEX (CDNOMA(IDECBL+1:),' ')
C
      IF (IPOSBL.LE.IDECBL) THEN
        ILCLNO=ILCDNO
      ELSEIF (CDNOMA(IPOSBL:).EQ.' ') THEN
        ILCLNO=IPOSBL-1
      ELSE
        IDECBL=IPOSBL
        GOTO 101
      ENDIF
C
      IF (ILCLNO.LE.LFI%JPNCPN) THEN
        CLNOMA=CDNOMA(:ILCLNO)
      ELSE
        CLNOMA=CDNOMA(:LFI%JPNCPN)
        ILCLNO=LFI%JPNCPN
        IREP=-15
        GOTO 1001
      ENDIF
C
      IF (IRANG.EQ.0) THEN
        IREP=-1
        GOTO 1001
      ENDIF
C
       IF (LFI%LMULTI) CALL LFIVER_MT (LFI, LFI%VERRUE(IRANG),'ON')
      LLVERF=LFI%LMULTI
C
      INBALO=LFI%MDES1D(IXM(LFI%JPNALO,IRANG))
C
      IF (INBALO.NE.0) THEN
C**
C     2.  -  EXPLORATION DES (PAIRES DE) PAGES ET ARTICLES D'INDEX,
C            A LA RECHERCHE DE L'ARTICLE LOGIQUE DEMANDE.
C-----------------------------------------------------------------------
C
        CALL LFIRAN_MT (LFI, IREP,IRANG,CLNOMA(:ILCLNO),
     S                  IRGPIM,IARTEX,IRETIN)
C
        IF (IRETIN.EQ.1) THEN
          GOTO 903
        ELSEIF (IRETIN.EQ.2) THEN
          GOTO 904
        ELSEIF (IRETIN.NE.0) THEN
          GOTO 1001
        ENDIF
C
      ELSE
        IARTEX=0
        IREP=0
      ENDIF
C
      IF (IARTEX.EQ.0) THEN
        KLONG=0
        KPOSEX=0
      ELSE
C
C        ON COMPLETE LES CARACTERISTIQUES DE L'ARTICLE.
C
        IF (.NOT.LFI%LPHASP(IRGPIM)) THEN
C
          CALL LFIPHA_MT (LFI, IREP,IRANG,IRGPIM,IRETIN)
C
          IF (IRETIN.EQ.1) THEN
            GOTO 903
          ELSEIF (IRETIN.EQ.2) THEN
            GOTO 904
          ELSEIF (IRETIN.NE.0) THEN
            GOTO 1001
          ENDIF
C
        ENDIF
C
        KLONG=LFI%MLGPOS(IXM(2*IARTEX-1,IRGPIM))
        KPOSEX=LFI%MLGPOS(IXM(2*IARTEX,IRGPIM))
C
C        On met a jour ce qui a trait aux acces pseudo-sequentiels...
C     ceci surtout pour ne pas faire 2 recherches dans l'index lors
C     d'un appel a LFILEC qui suivrait l'appel a LFINFO.
C
        LFI%NDERGF(IRANG)=LFI%JPNAPP*LFI%MFACTM(IRANG)*
     S                    (LFI%MRGPIF(IRGPIM)-1)+IARTEX
        LFI%CNDERA(IRANG)=CLNOMA(:ILCLNO)
        LFI%NSUIVF(IRANG)=LFI%JPNIL
        LFI%NPRECF(IRANG)=LFI%JPNIL
      ENDIF
C
      GOTO 1001
C**
C     9.  - CI-DESSOUS, ETIQUETTES DE BRANCHEMENT EN CAS D'ERREUR E/S.
C-----------------------------------------------------------------------
C
  903 CONTINUE
      CLACTI='WRITE'
      GOTO 909
C
  904 CONTINUE
      CLACTI='READ'
C
  909 CONTINUE
C
C      AU CAS OU, ON FORCE LE CODE-REPONSE ENTREE/SORTIE A ETRE POSITIF.
C
      IREP=IABS (IREP)
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE SOUS-PROGRAMME "LFIEMS" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
      KREP=IREP
      LLFATA=LLMOER (IREP,IRANG)
C
      IF (IRANG.NE.0) THEN
        LFI%NDEROP(IRANG)=7
        LFI%NDERCO(IRANG)=IREP
         IF (LLVERF) CALL LFIVER_MT (LFI, LFI%VERRUE(IRANG),'OFF')
      ENDIF
C
      IF (LLFATA.OR.IXNIMS (IRANG).EQ.2) THEN
        INIMES=2
      ELSE
        IF (LHOOK) CALL DR_HOOK('LFINFO_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='LFINFO'
      WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,
     S       '', CDNOMA='''''',A,'''''', KLONG='',I7,'', KPOSEX='',I8)')
     S     KREP,KNUMER,CLNOMA(:ILCLNO),KLONG,KPOSEX
      CALL LFIEMS_MT (LFI, KNUMER,INIMES,IREP,LLFATA,
     S                CLMESS,CLNSPR,CLACTI)
C
      IF (LHOOK) CALL DR_HOOK('LFINFO_MT',1,ZHOOK_HANDLE)
      END

