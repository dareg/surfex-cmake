C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFIPXA_MT (LFI, KREP, KNUMER, CDNOMA, CDSTRU, CDSUIV,
     S                    KLSUIV )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-programme Preparatoire a l'eXport d'un Article d'un
C     fichier LFI vers un systeme a priori different.
C
C     Il s'agit, en l'occurrence, de decrire la structure interne
C     de cet article en termes de types de variables.  
C**
C    ARGUMENTS : KREP   (Sortie) ==> Code-Reponse du sous-programme;
C                KNUMER (Entree) ==> Numero d'Unite Logique associe;
C                CDNOMA (Entree) ==> Nom de l'article decrit;
C                CDSTRU (Entree) ==> Structure interne de cet article;
C                CDSUIV (Sortie) ==> Nom de l'article suivant sur le
C                                    fichier, s'il en existe;
C                KLSUIV (Sortie) ==> Longueur de cet article.
C
C     (s'il n'y a pas d'article suivant, on retourne CDSUIV=' ' et
C      KLSUIV=0)
C
C     Les syntaxes autorisees pour CDSTRU sont decrites dans le sous-
C     programmes *LFIDST*.
C     
#ifndef f77
#include "precision.h"
#endif
C
      TYPE(LFICOM) :: LFI
      CHARACTER CDNOMA*(*), CDSUIV*(*), CDSTRU*(*)
      CHARACTER*(LFI%JPNCPN) CLNOMA, CLSUIV, CLSTRU
C
      INTEGER KREP, KNUMER, KLSUIV
      INTEGER ILONEX, ILCLNO, ILCDNO, IRANMX, IDECBL, IPOSBL, ILCDST
      INTEGER IRANG, IREP, INBALO, J, IRANIE, INIMES, IARTEX
      INTEGER IRGPIM, IRGPIF, IARTIC, IRETIN, ILCDSU, ILCLSU, ILCLST
      INTEGER ILUSTR
C
      LOGICAL LLVERF, LLOUVR, LLEXUL
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
      IF (LHOOK) CALL DR_HOOK('LFIPXA_MT',0,ZHOOK_HANDLE)
      CALL LFINUM_MT (LFI, KNUMER,IRANG)
      IREP=0
      LLVERF=.FALSE.
      ILCDNO=LEN (CDNOMA)
      ILCDSU=LEN (CDSUIV)
      ILCDST=LEN (CDSTRU)
      CLNOMA=' '
      ILCLNO=1
      CLSTRU=' '
      ILCLST=1
      CLSUIV=' '
      ILCLSU=1
      KLSUIV=0
C
      IF (ILCDNO.LE.0) THEN
        IREP=-15
        CLNOMA=LFI%CHINCO(:LFI%JPNCPN)
        ILCLNO=LFI%JPNCPN
      ELSEIF (CDNOMA.EQ.' ') THEN
        IREP=-18
      ENDIF
C
      IF (ILCDSU.LE.0) THEN
        IREP=-15
        CLSUIV=LFI%CHINCO(:LFI%JPNCPN)
        ILCLSU=LFI%JPNCPN
      ENDIF
C
      IF (ILCDST.LE.0) THEN
        IREP=-15
        CLSTRU=LFI%CHINCO(:LFI%JPNCPN)
        ILCLST=LFI%JPNCPN
      ELSEIF (CDSTRU.EQ.' ') THEN
        IREP=-39
      ENDIF
C
      IF (IREP.NE.0) THEN
        GOTO 1001
      ELSE
        CDSUIV=' '
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
C        Recherche de la longueur "utile" de la structure specifiee.
C        (c'est-a-dire sans tenir compte des blancs terminaux eventuels)
C
      IDECBL=0
C
  102 CONTINUE
      IPOSBL=IDECBL+INDEX (CDSTRU(IDECBL+1:),' ')
C
      IF (IPOSBL.LE.IDECBL) THEN
        ILUSTR=ILCDST
      ELSEIF (CDSTRU(IPOSBL:).EQ.' ') THEN
        ILUSTR=IPOSBL-1
      ELSE
        IDECBL=IPOSBL
        GOTO 102
      ENDIF
C
      ILCLST=MIN0 (ILCLST,ILUSTR)
C
      IF (IRANG.EQ.0) THEN
        IREP=-1
        GOTO 1001
      ENDIF
C
       IF (LFI%LMULTI) CALL LFIVER_MT (LFI, LFI%VERRUE(IRANG),'ON')
      LLVERF=LFI%LMULTI
      IRANIE=LFI%NEXPOR(IRANG)
C
      IF (IRANIE.LE.0) THEN
        IREP=-38
        CLACTI='EXPORT'
        GOTO 1001
      ENDIF
C
      IRANMX=LFI%NRCFMX(IRANIE)
      IARTEX=0
      ILONEX=0
      INBALO=LFI%MDES1D(IXM(LFI%JPNALO,IRANG))

      IF (INBALO.NE.0) THEN
C**
C     2.  -  EXPLORATION DES (PAIRES DE) PAGES ET ARTICLES D'INDEX,
C            A LA RECHERCHE DE L'ARTICLE LOGIQUE DEMANDE.
C-----------------------------------------------------------------------
C
        CALL LFIRAN_MT (LFI, IREP,IRANG,CLNOMA(:ILCLNO),IRGPIM,
     S                  IARTEX,IRETIN)
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
        CLACTI=CLNOMA(:ILCLNO)
        GOTO 1001
      ENDIF
C
C        ON COMPLETE LES CARACTERISTIQUES DE L'ARTICLE.
C
      IRGPIF=LFI%MRGPIF(IRGPIM)
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
      ILONEX=LFI%MLGPOS(IXM(2*IARTEX-1,IRGPIM))
C**
C     8.  -  RECHERCHE DE L'ARTICLE LOGIQUE DE DONNEES SUIVANT.
C-----------------------------------------------------------------------
C
      CALL LFICAX_MT (LFI, IREP,IRANG,IRGPIM,IARTIC,IRETIN)
C
      IF (IRETIN.EQ.1) THEN
        GOTO 903
      ELSEIF (IRETIN.EQ.2) THEN
        GOTO 904
      ELSEIF (IRETIN.NE.0.OR.IARTIC.EQ.0) THEN
        GOTO 1001
      ENDIF
C
      IRGPIF=LFI%MRGPIF(IRGPIM)
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
      KLSUIV=LFI%MLGPOS(IXM(IARTIC,IRGPIM))
      CLSUIV=LFI%CNOMAR(IXC(IARTIC,IRGPIM))
C
C        Recherche de la longueur "utile" du nom d'article.
C        (c'est-a-dire sans tenir compte des blancs terminaux eventuels)
C
      IDECBL=0
C
  811 CONTINUE
      IPOSBL=IDECBL+INDEX (CLSUIV(IDECBL+1:),' ')
C
      IF (IPOSBL.LE.IDECBL) THEN
        ILCLSU=LFI%JPNCPN
      ELSEIF (CLSUIV(IPOSBL:).EQ.' ') THEN
        ILCLSU=IPOSBL-1
      ELSE
        IDECBL=IPOSBL
        GOTO 811
      ENDIF
C
      IF (ILCDSU.GE.ILCLSU) THEN
        CDSUIV=CLSUIV(:ILCLNO)
      ELSE
        IREP=-24
        CLACTI=CLSUIV
        GOTO 1001
      ENDIF
C
      GOTO 1001
C**
C     9.  - CI-DESSOUS, ETIQUETTES DE BRANCHEMENT EN CAS D'ERREUR E/S.
C-----------------------------------------------------------------------
C
  901 CONTINUE
      CLACTI='INQUIRE'
      GOTO 909
C
  902 CONTINUE
      CLACTI='OPEN'
      GOTO 909
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
        LFI%NDEROP(IRANG)=22
        LFI%NDERCO(IRANG)=IREP
         IF (LLVERF) CALL LFIVER_MT (LFI, LFI%VERRUE(IRANG),'OFF')
      ENDIF
C
      IF (LLFATA.OR.IXNIMS (IRANG).EQ.2) THEN
        INIMES=2
      ELSE
        IF (LHOOK) CALL DR_HOOK('LFIPXA_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='LFIPXA'
      WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,
     S       '', CDNOMA='''''',A,'''''', CDSTRU='''''',A,
     S       '''''', CDSUIV='''''',A,'''''', KLSUIV='',I7)')
     S     KREP,KNUMER,CLNOMA(:ILCLNO),CLSTRU(:ILCLST),
     S     CLSUIV(:ILCDSU),KLSUIV
      CALL LFIEMS_MT (LFI, KNUMER,INIMES,IREP,LLFATA,CLMESS,
     S                CLNSPR,CLACTI)
C
      IF (LHOOK) CALL DR_HOOK('LFIPXA_MT',1,ZHOOK_HANDLE)
      END

