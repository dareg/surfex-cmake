C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFILAF_MT (LFI, KREP, KNUMER, LDTOUT )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-programme donnant, pour une unite logique ouverte au sens
C     du logiciel de fichiers indexes *LFI*, la Liste des Articles logi-
C     ques de donnees presents dans le Fichier, liste donnee toutefois
C     dans l'ordre PHYSIQUE ou ceux-ci figurent dans le fichier.
C        Sur option on donne aussi des renseignements sur les articles
C     (physiques) de gestion propres au logiciel, ainsi que sur les
C     trous repertories dans l'index.
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KNUMER (Entree) ==> Numero de l'unite logique;
C                LDTOUT (Entree) ==> Vrai si on doit donner les rensei-
C                                    gnements optionnels (qui ne concer-
C                                    nent pas directement les articles
C                                    logiques de donnees).
C
#ifndef f77
#include "precision.h"
#endif
C
      TYPE(LFICOM) :: LFI
      INTEGER KREP, KNUMER, IMDESC, IREP, IRANG, INTROU, INBPIR, INBALO
      INTEGER INALDO, IFACTM, ILARPH, INALPP, INTPPI, INPPIM, INIMES, J
      INTEGER INAGES, IRESER, INUTIL, IPERTE, IPOSFI, IPOSDE, INEXCE
      INTEGER INABAL, INALDI, INTROI, INPIMD, INPIMF, INPILE, JRGPIF
      INTEGER IRGPFS, IRGPIM, IRANGM, IRPIMS, INALPI, ILONGA, IRECPI
      INTEGER IDERPU, IREC, IRETIN
C
      LOGICAL LDTOUT
C
#include "lficom2.h"
#include "lficom_mt.h"
C**
C     1.  -  CONTROLES DES PARAMETRES D'APPEL, PUIS INITIALISATIONS.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFILAF_MT',0,ZHOOK_HANDLE)
      IREP=0
      IRANG=0
      CLNSPR='LFILAF'
      CALL LFINUM_MT (LFI, KNUMER,IRANG)
C
      IF (IRANG.EQ.0) THEN
        IREP=-1
        GOTO 1001
      ENDIF
C
      IF (LFI%LMULTI) CALL LFIVER_MT (LFI, LFI%VERRUE(IRANG),'ON')
      INTROU=LFI%MDES1D(IXM(LFI%JPNTRU,IRANG))+LFI%NBTROU(IRANG)
      INBPIR=LFI%MDES1D(IXM(LFI%JPNPIR,IRANG))
      INBALO=LFI%MDES1D(IXM(LFI%JPNALO,IRANG))
      INALDO=INBALO-INTROU
      IFACTM=LFI%MFACTM(IRANG)
      ILARPH=LFI%JPLARD*IFACTM
      INALPP=LFI%JPNAPP*IFACTM
      INTPPI=(INBALO-1+INALPP)/INALPP
      INPPIM=LFI%NPPIMM(IRANG)
C
C         Envoi d'une banniere.
C
      WRITE (UNIT=LFI%NULOUT,FMT='(///)')
C
      IF (LFI%LFRANC) THEN
        WRITE (UNIT=CLMESS,
     S         FMT='(''Catalogue de l''''Unite Logique LFI''
     S ,I3,'' dans l''''ordre *PHYSIQUE* (sequentiel) des articles'')')
     S     KNUMER
      ELSE
        WRITE (UNIT=CLMESS,FMT='(''Catalog of LFI Logical Unit'',I3,
     S         '' in *PHYSICAL* (sequential) record order'')') KNUMER
      ENDIF
C
      INIMES=2
      LLFATA=.FALSE.
      CALL LFIEMS_MT (LFI, KNUMER,INIMES,IREP,LLFATA,
     S                CLMESS,CLNSPR,CLACTI)
C**
C     2.  -  SUR OPTION, RENSEIGNEMENTS SUR LES ARTICLES "DE GESTION".
C            (ARTICLE DOCUMENTAIRE, PAIRES D'ARTICLES D'INDEX)
C-----------------------------------------------------------------------
C
      IF (LDTOUT) THEN
        INAGES=1+2*INBPIR
        IRESER=ILARPH*INAGES
C
        IF (LFI%LFRANC) THEN
          WRITE (UNIT=LFI%NULOUT,FMT='(//,TR1,I6,
     S           '' article(s) "physique(s)" de gestion,'',I6,
     S           '' mots chacun, occupant donc'',I7,'' mots; detail:'',
     S /,TR10,''Article documentaire de la position 1 a'',I6,/,TR10,I6,
     S'' paire(s) d''''articles d''''index prereserves, de la position''
     S           ,I6,'' a'',I7)')
     S         INAGES,ILARPH,IRESER,ILARPH,INBPIR,ILARPH+1,IRESER
        ELSE
          WRITE (UNIT=LFI%NULOUT,FMT='(//,TR1,I6,
     S           '' "physical" records for file handling,'',I6,
     S           '' words each, occupying then'',I7,'' words; detail:'',
     S /,TR10,''Documentary record from position 1 to'',I6,/,TR10,I6,
     S'' pair(s) of pre-reserved index records, from position''
     S           ,I6,'' to'',I7)')
     S         INAGES,ILARPH,IRESER,ILARPH,INBPIR,ILARPH+1,IRESER
        ENDIF
C
        IF (INTPPI.LT.INBPIR) THEN
          INUTIL=INBPIR-INTPPI
          IPERTE=ILARPH*INUTIL*2
C
          IF (LFI%LFRANC) THEN
            WRITE (UNIT=LFI%NULOUT,
     S             FMT='(/,TR10,5(''=''),''> Il y a'',I3,
     S '' paire(s) d''''articles d''''index inutilises, representant'',
     S             I8,'' mots'')') INUTIL,IPERTE
          ELSE
            WRITE (UNIT=LFI%NULOUT,
     S             FMT='(/,TR10,5(''=''),''> There is (are)'',I3,
     S '' pair(s) of unused index records, leading to a loss of'',
     S             I8,'' words'')') INUTIL,IPERTE
          ENDIF
C
        ELSEIF (INTPPI.EQ.INBPIR) THEN
C
          IF (LFI%LFRANC) THEN
            WRITE (UNIT=LFI%NULOUT,
     S             FMT='(TR15,5(''-''),TR3,''pas de paire '',
     S        ''d''''articles d''''index inutilises ni excedentaires'',
     S          TR3,5(''-''))')
          ELSE
            WRITE (UNIT=LFI%NULOUT,
     S             FMT='(TR15,5(''-''),TR3,''no pair of '',
     S        ''unused or overflow pages'',
     S          TR3,5(''-''))')
          ENDIF
C
        ELSEIF (INTPPI.EQ.(INBPIR+1)) THEN
          IPOSFI=ILARPH*(LFI%MDES1D(IXM(ILARPH,IRANG))+1)
          IPOSDE=IPOSFI-2*ILARPH+1
C
          IF (LFI%LFRANC) THEN
            WRITE (UNIT=LFI%NULOUT,
     S             FMT='(TR10,''une paire d''''articles '',
     S             ''d''''index excedentaires, de la position'',
     S             I9,'' a'',I9)')
     S      IPOSDE,IPOSFI
          ELSE
            WRITE (UNIT=LFI%NULOUT,
     S             FMT='(TR10,''one pair of overflow index '',
     S             ''pages ,from position'',
     S             I9,'' to'',I9)')
     S      IPOSDE,IPOSFI
          ENDIF
C
        ELSE
          INEXCE=INTPPI-INBPIR
C
          IF (LFI%LFRANC) THEN
            WRITE (UNIT=LFI%NULOUT,
     S             FMT='(TR10,I6,'' paires d''''articles '',
     S           ''d''''index excedentaires, des positions:'')') INEXCE
C
            DO 201 J=1,INEXCE
            IPOSFI=ILARPH*(LFI%MDES1D(IXM(ILARPH+1-J,IRANG))+1)
            IPOSDE=IPOSFI-2*ILARPH+1
            WRITE (UNIT=LFI%NULOUT,FMT='(TR20,I9,'' a'',I9)') 
     S             IPOSDE,IPOSFI
  201       CONTINUE
C
          ELSE
            WRITE (UNIT=LFI%NULOUT,
     S             FMT='(TR10,I6,'' pairs of overflow index '',
     S           ''pages, from positions:'')') INEXCE
C
            DO 202 J=1,INEXCE
            IPOSFI=ILARPH*(LFI%MDES1D(IXM(ILARPH+1-J,IRANG))+1)
            IPOSDE=IPOSFI-2*ILARPH+1
            WRITE (UNIT=LFI%NULOUT,FMT='(TR20,I9,'' to'',I9)') 
     S             IPOSDE,IPOSFI
  202       CONTINUE
C
          ENDIF
C
        ENDIF
C
      ENDIF
C
      WRITE (UNIT=LFI%NULOUT,FMT='(//)')
C**
C     3.  -  RENSEIGNEMENTS INDIVIDUALISES SUR LES ARTICLES LOGIQUES.
C            (DONNEES, ET SUR OPTION TROUS REPERTORIES DANS L'INDEX)
C-----------------------------------------------------------------------
C
      IF (LFI%LFRANC) THEN
C
        IF (INBALO.EQ.0) THEN
          WRITE (UNIT=LFI%NULOUT,
     S           FMT='(/,TR10,5(''=''),''> L''''unite logique'',
     S I3,'' ne contient AUCUN ARTICLE LOGIQUE (ni donnees, ni trous)'',
     S           //)') KNUMER
          GOTO 1001
        ELSEIF (INBALO.EQ.INTROU) THEN
          WRITE (UNIT=LFI%NULOUT,
     S           FMT='(/,TR10,5(''=''),''> L''''unite logique'',
     S I3,'' ne contient QUE DES TROUS, pas de donnees)'',//)') KNUMER
          IF (.NOT.LDTOUT) GOTO 1001
        ENDIF
C
      ELSE
C
        IF (INBALO.EQ.0) THEN
          WRITE (UNIT=LFI%NULOUT,
     S           FMT='(/,TR10,5(''=''),''> The logical unit'',I3,
     S '' contains NO LOGICAL RECORD AT ALL (neither data, nor holes)'',
     S           //)') KNUMER
          GOTO 1001
        ELSEIF (INBALO.EQ.INTROU) THEN
          WRITE (UNIT=LFI%NULOUT,
     S           FMT='(/,TR10,5(''=''),''> The logical unit'',I3,
     S '' contains ONLY HOLES, no dat)'',//)') KNUMER
          IF (.NOT.LDTOUT) GOTO 1001
        ENDIF
C
      ENDIF
C*
C     3.1 -  BALAYAGE DES PAIRES D'ARTICLES D'INDEX, PAR ORDRE CROISSANT
C-----------------------------------------------------------------------
C
      INABAL=0
      INALDI=0
      INTROI=0
      INPIMD=2
      INPIMF=INPPIM
      IF (LFI%NPODPI(IRANG).EQ.2) INPIMD=3
      IF (LFI%NPODPI(IRANG).EQ.INPPIM) INPIMF=INPPIM-1
      INPILE=2
C
      DO 319 JRGPIF=1,INTPPI
      IRGPFS=JRGPIF+1
C
C        On fait en sorte que la P.A.I. concernee, ainsi que sa suivante
C     eventuelle, soient toutes les deux en memoire.
C
      IF (JRGPIF.EQ.INTPPI) THEN
        IRGPIM=LFI%MRGPIM(LFI%NPODPI(IRANG),IRANG)
        GOTO 314
C
      ELSEIF (JRGPIF.NE.1) THEN
C
C       Recherche de la P.A.I. dans les Paires de Pages d'Index memoire.
C
        DO 311 J=INPIMD,INPIMF
        IRGPIM=LFI%MRGPIM(J,IRANG)
C
        IF (LFI%MRGPIF(IRGPIM).EQ.JRGPIF) THEN
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
          GOTO 312
C
        ENDIF
C
  311   CONTINUE
C
C          Mise en memoire de la Paire d'Articles d'Index cherchee.
C
        CALL LFIPIM_MT (LFI, IREP,IRANG,IRANGM,IRGPIM,JRGPIF,IRGPFS,
     S                  INPILE,IRETIN)
C
        IF (IRETIN.EQ.1) THEN
          GOTO 903
        ELSEIF (IRETIN.EQ.2) THEN
          GOTO 904
        ELSEIF (IRETIN.NE.0) THEN
          GOTO 1001
        ELSEIF (IRANGM.GT.INPPIM) THEN
          INPPIM=IRANGM
          INPIMF=INPPIM
        ENDIF
C
      ELSE
        IRGPIM=LFI%MRGPIM(1,IRANG)
C
      ENDIF
C
  312 CONTINUE
C
      IF (IRGPFS.EQ.INTPPI) THEN
        IRPIMS=LFI%MRGPIM(LFI%NPODPI(IRANG),IRANG)
C
      ELSE
C
C       Recherche de la P.A.I. dans les Paires de Pages d'Index memoire.
C
        DO 313 J=INPIMD,INPIMF
        IRPIMS=LFI%MRGPIM(J,IRANG)
C
        IF (LFI%MRGPIF(IRPIMS).EQ.IRGPFS) THEN
C
          IF (.NOT.LFI%LPHASP(IRPIMS)) THEN
C
            CALL LFIPHA_MT (LFI, IREP,IRANG,IRPIMS,IRETIN)
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
          GOTO 314
C
        ENDIF
C
  313   CONTINUE
C
C          Mise en memoire de la Paire d'Articles d'Index cherchee.
C
        CALL LFIPIM_MT (LFI, IREP,IRANG,IRANGM,IRPIMS,IRGPFS,JRGPIF,
     S                  INPILE,IRETIN)
C
        IF (IRETIN.EQ.1) THEN
          GOTO 903
        ELSEIF (IRETIN.EQ.2) THEN
          GOTO 904
        ELSEIF (IRETIN.NE.0) THEN
          GOTO 1001
        ELSEIF (IRANGM.GT.INPPIM) THEN
          INPPIM=IRANGM
          INPIMF=INPPIM
        ENDIF
C
      ENDIF
C
  314 CONTINUE
      INALPI=MIN0 (INALPP,INBALO-INABAL)
C
C        Balayage de la Paire d'Article d'Index concernee.
C
      DO 318 J=1,INALPI
C
      IF (LFI%CNOMAR(IXC(J,IRGPIM)).NE.' ') THEN
C
C              Il s'agit d'un article logique de donnees; en plus de ses
C         caracteristiques tabulees, on verifie s'il n'y a pas de la
C         place "perdue" juste derriere les donnees, place recuperable
C         eventuellement en cas de reecriture plus longue de l'article
C         logique.
C
        INALDI=INALDI+1
        ILONGA=LFI%MLGPOS(IXM(2*J-1,IRGPIM))
        IPOSDE=LFI%MLGPOS(IXM(2*J  ,IRGPIM))
        IPOSFI=IPOSDE+ILONGA-1
C
        IF (J.EQ.1.AND.JRGPIF.GT.INBPIR) THEN
C
C          Cas du premier article logique d'une P.A.I. excedentaire;
C     dans ce cas, la P.A.I. est situee derriere l'article logique,
C     en occupant deux articles physiques.
C
          IRECPI=LFI%MDES1D(IXM(ILARPH+1-(JRGPIF-INBPIR),IRANG))
          IDERPU=ILARPH*(IRECPI-1)
C
        ELSEIF (J.EQ.INALPI.AND.JRGPIF.EQ.INTPPI) THEN
C
C          Cas du dernier article logique du fichier, sans P.A.I. situee
C     derriere: la derniere position utilisable sans modifier le nombre
C     d'articles physiques du fichier correspond a la fin du dernier
C     article physique contenant des donnees, ou a la fin du dernier
C     article physique ecrit sur le fichier.
C
          IMDESC=LFI%MDES1D(IXM(LFI%JPNAPH,IRANG))
          IREC=MAX0 (1+(IPOSFI-1)/ILARPH,IMDESC)
          IDERPU=ILARPH*IREC
C
C          Si on arrive au test ci-dessous, on est sur que l'article lo-
C     gique n'est pas le dernier du fichier.
C
        ELSEIF (J.NE.INALPP) THEN
C
C          Cas general, ou l'article logique n'est pas le dernier de sa
C     (Paire de) Page(s) d'Index.
C
          IDERPU=LFI%MLGPOS(IXM(2*J+2,IRGPIM))-1
C
        ELSE
C
C          Cas particulier ou l'article logique est le dernier de sa
C     (Paire de) Page(s) d'Index.
C
          IDERPU=LFI%MLGPOS(IXM(2,IRPIMS))-1
        ENDIF
C
        IF (IDERPU.EQ.IPOSFI) THEN
C
          IF (LFI%LFRANC) THEN
            WRITE (UNIT=LFI%NULOUT,
     S             FMT='(I7,''-eme article de donnees: "'',A,
     S             ''",'',I7,'' mots, position'',I9,'' a'',I9)')
     S       INALDI,LFI%CNOMAR(IXC(J,IRGPIM)),ILONGA,IPOSDE,IPOSFI
          ELSE
            WRITE (UNIT=LFI%NULOUT,
     S             FMT='(I7,''-th data record: "'',A,''",'',I7,
     S             '' words, position'',I9,'' to'',I9)')
     S       INALDI,LFI%CNOMAR(IXC(J,IRGPIM)),ILONGA,IPOSDE,IPOSFI
          ENDIF
C
        ELSE
C
C           On visualise en plus la place "perdue" derriere l'article.
C
          IF (LFI%LFRANC) THEN
            WRITE (UNIT=LFI%NULOUT,
     S             FMT='(I7,''-eme article de donnees: "'',A,
     S             ''",'',I7,'' mots, position'',I9,'' a'',I9,'' <'',SP,
     S             I8,'' >'')')
     S   INALDI,LFI%CNOMAR(IXC(J,IRGPIM)),ILONGA,IPOSDE,
     S   IPOSFI,IDERPU-IPOSFI
          ELSE
            WRITE (UNIT=LFI%NULOUT,
     S       FMT='(I7,''-th data record: '''''',A,
     S       '''''','',I7,'' words, position'',I9,'' to'',I9,'' <'',
     S       SP,I8,'' >'')')
     S   INALDI,LFI%CNOMAR(IXC(J,IRGPIM)),ILONGA,IPOSDE,
     S   IPOSFI,IDERPU-IPOSFI
          ENDIF
C
        ENDIF
C
      ELSEIF (LDTOUT) THEN
        INTROI=INTROI+1
        ILONGA=LFI%MLGPOS(IXM(2*J-1,IRGPIM))
        IPOSDE=LFI%MLGPOS(IXM(2*J  ,IRGPIM))
        IPOSFI=IPOSDE+ILONGA-1
C
        IF (LFI%LFRANC) THEN
          WRITE (UNIT=LFI%NULOUT,FMT='(TR1,5(''=''),''>'',T10,I6,
     S ''-eme TROU repertorie dans l''''index, longueur reutilisable:'',
     S         I7,'' mots, position'',I9,'' a'',I9)')
     S   INTROI,ILONGA,IPOSDE,IPOSFI
        ELSE
          WRITE (UNIT=LFI%NULOUT,FMT='(TR1,5(''=''),''>'',T10,I6,
     S ''-th HOLE cataloged within index, re-usable length:'',
     S         I7,'' words, position'',I9,'' to'',I9)')
     S   INTROI,ILONGA,IPOSDE,IPOSFI
        ENDIF
C
      ENDIF
C
  318 CONTINUE
C
      INABAL=INABAL+INALPI
  319 CONTINUE
C*
C     3.2 -  ENVOI DE MESSAGES RECAPITULATIFS.
C-----------------------------------------------------------------------
C
      IF (LFI%LFRANC) THEN
C
        IF (LDTOUT) THEN
          WRITE (UNIT=LFI%NULOUT,FMT='(//,T5,8(''-''),TR3,I7,
     S           '' articles logiques de donnees et'',I6,
     S           '' trous repertories listes'',TR3,8(''-''),//)')
     S    INALDI,INTROI
        ELSE
          WRITE (UNIT=LFI%NULOUT,FMT='(//,T5,8(''-''),TR3,I7,
     S       '' articles logiques de donnees listes'',TR3,8(''-''),//)')
     S    INALDI
        ENDIF
C
      ELSE
C
        IF (LDTOUT) THEN
          WRITE (UNIT=LFI%NULOUT,FMT='(//,T5,8(''-''),TR3,I7,
     S           '' logical records of data and'',I6,
     S           '' holes within index listed'',TR3,8(''-''),//)')
     S    INALDI,INTROI
        ELSE
          WRITE (UNIT=LFI%NULOUT,FMT='(//,T5,8(''-''),TR3,I7,
     S       '' logical records of data listed'',TR3,8(''-''),//)')
     S    INALDI
        ENDIF
C
      ENDIF
C
      IF (INALDI.EQ.INALDO.AND.(.NOT.LDTOUT.OR.INTROI.EQ.INTROU)) THEN
C
        IF (LFI%LFRANC) THEN
          WRITE (UNIT=CLMESS,FMT=
     S     '(''Fin du catalogue de l''''Unite Logique'',I3,'' ---'',I7,
     S       '' Articles logiques en tout'')') KNUMER,INBALO
        ELSE
          WRITE (UNIT=CLMESS,FMT=
     S     '(''End of catalog of Logical Unit'',I3,'' ---'',I7,
     S       '' logical Records for whole file'')') KNUMER,INBALO
        ENDIF
C
        CALL LFIEMS_MT (LFI, KNUMER,INIMES,IREP,LLFATA,
     S                  CLMESS,CLNSPR,CLACTI)
        WRITE (UNIT=LFI%NULOUT,FMT='(///)')
      ELSE
        IREP=-16
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
        LFI%NDEROP(IRANG)=18
        LFI%NDERCO(IRANG)=IREP
        IF (LFI%LMULTI) CALL LFIVER_MT (LFI, LFI%VERRUE(IRANG),'OFF')
      ENDIF
C
      IF (LLFATA.OR.IXNIMS (IRANG).EQ.2) THEN
        INIMES=2
      ELSE
        IF (LHOOK) CALL DR_HOOK('LFILAF_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,
     S    '', LDTOUT= '',L1)') KREP,KNUMER,LDTOUT
      CALL LFIEMS_MT (LFI, KNUMER,INIMES,IREP,LLFATA,
     S                CLMESS,CLNSPR,CLACTI)
C
      IF (LHOOK) CALL DR_HOOK('LFILAF_MT',1,ZHOOK_HANDLE)
      END

