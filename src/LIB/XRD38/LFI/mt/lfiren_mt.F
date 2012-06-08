C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFIREN_MT (LFI, KREP, KNUMER, CDNOM1, CDNOM2 )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        SOUS-PROGRAMME PERMETTANT DE RENOMMER UN ARTICLE (DE DONNEES)
C     SUR UNE UNITE LOGIQUE OUVERTE POUR LE LOGICIEL DE FICHIERS INDEXES
C     *LFI*. LE NOUVEAU NOM D'ARTICLE NE DOIT PAS Y ETRE DEJA UTILISE.
C**
C    ARGUMENTS : KREP   (SORTIE) ==> CODE-REPONSE DU SOUS-PROGRAMME;
C                KNUMER (ENTREE) ==> LFI%NUMERO DE L'UNITE LOGIQUE;
C                CDNOM1 (ENTREE) ==> NOM DE L'ARTICLE A RENOMMER;
C                CDNOM2 (ENTREE) ==> NOUVEAU NOM A DONNER A L'ARTICLE.
C
#ifndef f77
#include "precision.h"
#endif
C
      TYPE(LFICOM) :: LFI
      CHARACTER CDNOM1*(*), CDNOM2*(*), CLNOM1*(LFI%JPNCPN), 
     S          CLNOM2*(LFI%JPNCPN)
C
      INTEGER KREP, KNUMER, IRANG, IREP, ILCDN1, ILCLN1, ILCDN2, ILCLN2
      INTEGER IDECBL, IPOSBL, IARTEX, INBALO, IRGPIM, IRETIN, INIMES
C
      LOGICAL LLECR, LLVERF
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
      IF (LHOOK) CALL DR_HOOK('LFIREN_MT',0,ZHOOK_HANDLE)
      CALL LFINUM_MT (LFI, KNUMER,IRANG)
      LLVERF=.FALSE.
      IREP=0
      LLECR=.FALSE.
      ILCDN1=LEN (CDNOM1)
      ILCDN2=LEN (CDNOM2)
C
      IF (MIN0 (ILCDN1,ILCDN2).LE.0) THEN
C
        IREP=-15
C
        IF (ILCDN1.LE.0) THEN
          CLNOM1=LFI%CHINCO(:LFI%JPNCPN)
          ILCLN1=LFI%JPNCPN
        ELSE
          ILCLN1=MIN0 (ILCDN1,LFI%JPNCPN)
          CLNOM1=CDNOM1(:ILCLN1)
        ENDIF
C
        IF (ILCDN2.LE.0) THEN
          CLNOM2=LFI%CHINCO(:LFI%JPNCPN)
          ILCLN2=LFI%JPNCPN
        ELSE
          ILCLN2=MIN0 (ILCDN2,LFI%JPNCPN)
          CLNOM2=CDNOM2(:ILCLN2)
        ENDIF
C
        GOTO 1001
C
      ELSEIF (CDNOM1.EQ.' '.OR.CDNOM2.EQ.' ') THEN
C
        IREP=-18
C
        IF (CDNOM1.EQ.' ') THEN
          CLNOM1=' '
          ILCLN1=1
        ELSE
          ILCLN1=MIN0 (ILCDN1,LFI%JPNCPN)
          CLNOM1=CDNOM1(:ILCLN1)
        ENDIF
C
        IF (CDNOM2.EQ.' ') THEN
          CLNOM2=' '
          ILCLN2=1
        ELSE
          ILCLN2=MIN0 (ILCDN2,LFI%JPNCPN)
          CLNOM2=CDNOM2(:ILCLN2)
        ENDIF
C
        GOTO 1001
C
      ENDIF
C
C        Recherche de la longueur "utile" des noms d'article specifies.
C        (c'est-a-dire sans tenir compte des blancs terminaux eventuels)
C
      IDECBL=0
C
  101 CONTINUE
      IPOSBL=IDECBL+INDEX (CDNOM1(IDECBL+1:),' ')
C
      IF (IPOSBL.LE.IDECBL) THEN
        ILCLN1=ILCDN1
      ELSEIF (CDNOM1(IPOSBL:).EQ.' ') THEN
        ILCLN1=IPOSBL-1
      ELSE
        IDECBL=IPOSBL
        GOTO 101
      ENDIF
C
      IDECBL=0
C
  102 CONTINUE
      IPOSBL=IDECBL+INDEX (CDNOM2(IDECBL+1:),' ')
C
      IF (IPOSBL.LE.IDECBL) THEN
        ILCLN2=ILCDN2
      ELSEIF (CDNOM2(IPOSBL:).EQ.' ') THEN
        ILCLN2=IPOSBL-1
      ELSE
        IDECBL=IPOSBL
        GOTO 102
      ENDIF
C
      IF (ILCLN1.GT.LFI%JPNCPN) THEN
        ILCLN1=LFI%JPNCPN
        IREP=-15
      ENDIF
C
      IF (ILCLN2.GT.LFI%JPNCPN) THEN
        ILCLN2=LFI%JPNCPN
        IREP=-15
      ENDIF
C
      CLNOM1=CDNOM1(:ILCLN1)
      CLNOM2=CDNOM2(:ILCLN2)
      IF (IREP.NE.0) GOTO 1001
C
      IF (IRANG.EQ.0) THEN
        IREP=-1
        GOTO 1001
      ENDIF
C
       IF (LFI%LMULTI) CALL LFIVER_MT (LFI, LFI%VERRUE(IRANG),'ON')
      LLVERF=LFI%LMULTI
C
      IF (LFI%NEXPOR(IRANG).GT.0) THEN
C
C         Fichier en cours d'export... la seule modification acceptee
C         est l'ajout de nouveaux articles.
C
        IREP=-37
        GOTO 1001
      ENDIF
C
      IARTEX=0
      INBALO=LFI%MDES1D(IXM(LFI%JPNALO,IRANG))
C
      IF (INBALO.NE.0) THEN
C**
C     2.  -  EXPLORATION DES (PAIRES DE) PAGES ET ARTICLES D'INDEX,
C            A LA RECHERCHE DU NOUVEAU NOM D'ARTICLE, QUI NE DOIT
C            PAS ETRE LE NOM D'UN ARTICLE EXISTANT.
C-----------------------------------------------------------------------
C
        CALL LFIRAN_MT (LFI, IREP,IRANG,CLNOM2(:ILCLN2),
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
        IF (IARTEX.NE.0) THEN
          IREP=-25
          CLACTI=CLNOM2(:ILCLN2)
          GOTO 1001
        ENDIF
C**
C     3.  -  EXPLORATION DES (PAIRES DE) PAGES ET ARTICLES D'INDEX,
C            A LA RECHERCHE DE L'ARTICLE LOGIQUE A RENOMMER.
C-----------------------------------------------------------------------
C
        CALL LFIRAN_MT (LFI, IREP,IRANG,CLNOM1(:ILCLN1),
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
      ENDIF
C
      IF (IARTEX.EQ.0) THEN
        IREP=-20
        CLACTI=CLNOM1(:ILCLN1)
        GOTO 1001
      ENDIF
C**
C     4.  -  TOUT EST OK... ON EFFECTUE LE CHANGEMENT DE NOM.
C-----------------------------------------------------------------------
C
      LFI%CNOMAR(IXC(IARTEX,IRGPIM))=CLNOM2(:ILCLN2)
      LFI%LECRPI(IRGPIM,1)=.TRUE.
      LFI%NBRENO(IRANG)=LFI%NBRENO(IRANG)+1
C
C        On met a jour ce qui a trait aux acces pseudo-sequentiels...
C
      LFI%NDERGF(IRANG)=LFI%JPNAPP*LFI%MFACTM(IRANG)*
     S                  (LFI%MRGPIF(IRGPIM)-1)+IARTEX
      LFI%CNDERA(IRANG)=CLNOM2(:ILCLN2)
      LFI%NSUIVF(IRANG)=LFI%JPNIL
      LFI%NPRECF(IRANG)=LFI%JPNIL
C
      IF (.NOT.LFI%LMODIF(IRANG)) THEN
C
C         CAS DE LA PREMIERE MODIFICATION DEPUIS L'OUVERTURE DU FICHIER.
C
        LFI%LMODIF(IRANG)=.TRUE.
        CALL LFIMOE_MT (LFI, IREP,IRANG,IRETIN)
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
        LFI%NDEROP(IRANG)=13
        LFI%NDERCO(IRANG)=IREP
         IF (LLVERF) CALL LFIVER_MT (LFI, LFI%VERRUE(IRANG),'OFF')
      ENDIF
C
      IF (LLFATA.OR.IXNIMS (IRANG).EQ.2) THEN
        INIMES=2
      ELSE
        IF (LHOOK) CALL DR_HOOK('LFIREN_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='LFIREN'
      WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,
     S       '', CDNOM1='''''',A,'''''', CDNOM2='''''',A,'''''''')')
     S     KREP,KNUMER,CLNOM1(:ILCLN1),CLNOM2(:ILCLN2)
      CALL LFIEMS_MT (LFI, KNUMER,INIMES,IREP,LLFATA,
     S                CLMESS,CLNSPR,CLACTI)
C
      IF (LHOOK) CALL DR_HOOK('LFIREN_MT',1,ZHOOK_HANDLE)
      END

