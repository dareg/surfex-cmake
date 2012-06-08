C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFIPXF_MT (LFI, KREP, KNUMER, KNUMEX, CDCFGX, 
     S                      KLAREX, KXCNEX,
     S                      KFACEX, KNUTRA, CDNOMA, KLONG )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-programme Preparatoire a la realisation d'une
C        "version eXport" d'un Fichier LFI vers un systeme
C        a priori different. La methode utilisee suppose:
C
C     - que les fichiers a acces direct FORTRAN soient implantes ou
C       traitables comme des fichiers non bloques sur le systeme
C       destinataire;
C
C     - que l'on puisse ecrire par WRITE FORTRAN des fichiers
C       non bloques sur le systeme ou est fait la "version export",
C       aussi appelee "fichier export";
C
C     - que la conversion des variables, numeriques voire aussi
C       caracteres, soit faite au niveau des couches d'entrees/sorties
C       FORTRAN sur le systeme ou est fait le fichier export;
C
C     ( dans la pratique, les deux points qui precedent impliquent un
C       parametrage au niveau du langage de controle, a priori )
C
C     - que le programme utilisateur ait ouvert au prealable le fichier
C       LFI dont on veut realiser une "version export", et appelle le
C       sous-programme LFIPXF;
C
C     - que le programme utilisateur specifie le contenu des articles
C       a exporter en termes de types FORTRAN; ceci pouvant se faire
C       de deux manieres, eventuellement combinables:
C
C       1) Si le fichier contient (essentiellement) des donnees
C          utilisateur pouvant se decrire de maniere homogene,
C          par exemple rien que des variables reelles, et doit etre
C          exporte dans la totalite des articles, alors il suffira
C          d'appeler le sous-programme LFIXPH avec la description
C          correspondant a ces articles, que l'on peut aussi voir
C          comme une "description par defaut" (ou implicite);
C
C       2) Si ce n'est pas le cas, ou si une partie des articles ne
C          peut pas etre decrite de la meme maniere que les autres
C          articles, alors il faudra  que le programme utilisateur
C          specifie, pour chacun de ces articles,
C          le contenu en termes de types FORTRAN en appelant le sous-
C          programme LFIXPA: il s'agit la d'une description explicite,
C          ayant precedence sur une eventuelle description implicite;
C
C     - qu'en fin de compte le programme utilisateur appelle le sous-
C       programme LFIXPF qui fabriquera vraiment la version export,
C       a partir des specifications donnees via LFIPXF, LFIXPH, LFIXPA.
C**
C    ARGUMENTS : KREP   (Sortie) ==> Code-Reponse du sous-programme;
C                KNUMER (Entree) ==> Numero d'Unite Logique associe
C                                    au fichier LFI a exporter;
C                KNUMEX (Entree) ==> Numero d'Unite Logique associe
C                                    a la version export a realiser;
C                CDCFGX (Entree) ==> Configuration du systeme
C                                    destinataire du fichier export;
C                KLAREX (Entree) ==> Longueur d'ARticle Elementaire du
C                                    logiciel LFI du systeme destinatai-
C                                    re, exprimee en mots du systeme
C                                    destinataire;
C                                    (LFI%JPLARD du logiciel "distant")
C                KXCNEX (Entree) ==> Nombre maXimum de Caracteres par
C                                    Nom d'article du logiciel LFI du
C                                    systeme destinataire;
C                                    (LFI%JPNCPN du logiciel "distant")
C                KFACEX (Entree) ==> Facteur multiplicatif du fichier
C                                    export;
C                KNUTRA (Entree) ==> Numero d'Unite Logique utilisable
C                                    pour un fichier de travail eventuel
C                                    (si utilisation de LFIXPA), de type
C                                    LFI;
C                CDNOMA (Sortie) ==> Nom du premier article "candidat"
C                                    (potentiel) a l'export;
C                KLONG  (Sortie) ==> Longueur de cet article.
C
C     REMARQUE: Le fichier de travail n'est utilise que si l'on n'a pas
C               assez de place dans les tables pour stocker en memoire
C               les descripteurs. Mais si on en a besoin, il faut penser
C               que ce fichier occupera (temporairement, jusqu'a appel a
C               LFIXPF) une entree dans les tables LFI, et donc ne pas
C               avoir les tables saturees auparavant.
C
#ifndef f77
#include "precision.h"
#endif
C
      TYPE(LFICOM) :: LFI
      CHARACTER CDNOMA*(*), CLNOMA*(LFI%JPNCPN)
      CHARACTER CDCFGX*(*), CLCFGX*(LFI%JPXCCF)
C
      INTEGER KREP, KNUMER, KNUMEX, KLAREX, KXCNEX, KFACEX, KNUTRA
      INTEGER KLONG, ILCLNO, ILCDNO, IRANMX, ILCFGX, IDECBL, IPOSBL
      INTEGER IRANG, IREP, INUMER, INBALO, INTTRU, J, IRANIE, INIMES
      INTEGER IRGPIM, IRGPIF, IARTIC, IRETIN, ILCDCF
C
      LOGICAL LLVERG, LLVERF, LLEXUL, LLOUVR
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
      IF (LHOOK) CALL DR_HOOK('LFIPXF_MT',0,ZHOOK_HANDLE)
      CALL LFINUM_MT (LFI, KNUMER,IRANG)
      IREP=0
      INUMER=KNUMER
      LLVERF=.FALSE.
      LLVERG=.FALSE.
      ILCDNO=LEN (CDNOMA)
      ILCDCF=LEN (CDCFGX)
      CLNOMA=' '
      ILCLNO=1
      CLCFGX=' '
      ILCFGX=1
      KLONG=0
C
      IF (MIN0 (KLAREX,KXCNEX,KFACEX).LE.0) THEN
        IREP=-14
        GOTO 1001
      ELSEIF (ILCDNO.LE.0) THEN
        IREP=-15
        CLNOMA=LFI%CHINCO(:LFI%JPNCPN)
        ILCLNO=LFI%JPNCPN
      ENDIF
C
      IF (ILCDCF.LE.0) THEN
        IREP=-15
        CLCFGX=LFI%CHINCO(:LFI%JPNCPN)
        ILCFGX=LFI%JPNCPN
      ENDIF
C
      IF (IREP.NE.0) THEN
        GOTO 1001
      ELSE
        CDNOMA=' '
      ENDIF
C
C        Recherche de la longueur "utile" de la configuration specifiee.
C        (c'est-a-dire sans tenir compte des blancs terminaux eventuels)
C
      IDECBL=0
C
  101 CONTINUE
      IPOSBL=IDECBL+INDEX (CDCFGX(IDECBL+1:),' ')
C
      IF (IPOSBL.LE.IDECBL) THEN
        ILCFGX=ILCDCF
      ELSEIF (CDCFGX(IPOSBL:).EQ.' ') THEN
        ILCFGX=IPOSBL-1
      ELSE
        IDECBL=IPOSBL
        GOTO 101
      ENDIF
C
      IF (ILCFGX.LE.LFI%JPXCCF) THEN
        CLCFGX=CDCFGX(:ILCFGX)
      ELSE
        CLCFGX=CDCFGX(:LFI%JPXCCF)
        ILCFGX=LFI%JPXCCF
        IREP=-15
        GOTO 1001
      ENDIF
C
      DO 102 J=0,LFI%JPCFMX
C
      IF (CDCFGX.EQ.LFI%CFGMXD(J)) THEN
        IRANMX=J
        GOTO 103
      ENDIF      
C
  102 CONTINUE
C
C        Configuration du systeme destinataire inconnue ou non prevue.
C
      IREP=-32
      GOTO 1001
C
  103 CONTINUE
C
      IF (KXCNEX.GT.LFI%JPXCIE) THEN
        IREP=-33
        GOTO 1001
      ENDIF
C
C        Controle de validite FORTRAN et de non ouverture prealable
C        des Numeros d'Unite Logique KNUMEX et KNUTRA.
C
      INUMER=KNUMEX
      INQUIRE (UNIT=KNUMEX,EXIST=LLEXUL,OPENED=LLOUVR,ERR=901,
     S         IOSTAT=IREP)
      CLACTI='EXPORT'
C
      IF (.NOT.LLEXUL) THEN
        IREP=-30
        GOTO 1001
      ELSEIF (LLOUVR) THEN
        IREP=-34
        GOTO 1001
      ENDIF
C
      INUMER=KNUTRA
      INQUIRE (UNIT=KNUTRA,EXIST=LLEXUL,OPENED=LLOUVR,ERR=901,
     S         IOSTAT=IREP)
C
      IF (LFI%LFRANC) THEN
        CLACTI='DE TRAVAIL'
      ELSE
        CLACTI='WORK'
      ENDIF
C
      IF (.NOT.LLEXUL) THEN
        IREP=-30
        GOTO 1001
      ELSEIF (LLOUVR) THEN
        IREP=-34
        GOTO 1001
      ENDIF
C
      INUMER=KNUMER
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
        IREP=-35
        CLACTI='EXPORT'
        GOTO 1001
      ELSEIF (LFI%NIMPOR(IRANG).GT.0) THEN
        IREP=-35
        CLACTI='IMPORT'
        GOTO 1001
      ENDIF
C
      INBALO=LFI%MDES1D(IXM(LFI%JPNALO,IRANG))
      INTTRU=LFI%MDES1D(IXM(LFI%JPNTRU,IRANG))+LFI%NBTROU(IRANG)
C
      IF (INBALO.EQ.INTTRU) THEN
C
C         Fichier vide de donnees... inexportable.
C
        IREP=-36
        CLACTI='EXPORT'
        GOTO 1001
      ENDIF
C
C               Ouverture de l'unite logique KNUMEX.
C
      INUMER=KNUMEX
      OPEN (UNIT=KNUMEX,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
     S      FORM='UNFORMATTED',IOSTAT=IREP,ERR=902)
      REWIND (UNIT=KNUMEX,IOSTAT=IREP,ERR=906)
      INUMER=KNUMER
C**
C     2.  -  RECHERCHE DU PREMIER ARTICLE LOGIQUE DE DONNEES DU FICHIER.
C-----------------------------------------------------------------------
C
C       Reinitialisation des caracteristiques de type "pointeur".
C
      LFI%NDERGF(IRANG)=LFI%JPNIL
      LFI%CNDERA(IRANG)=' '
      LFI%NSUIVF(IRANG)=LFI%JPNIL
      LFI%NPRECF(IRANG)=LFI%JPNIL
C
      CALL LFICAX_MT (LFI, IREP,IRANG,IRGPIM,IARTIC,IRETIN)
C
      IF (IRETIN.EQ.1) THEN
        GOTO 903
      ELSEIF (IRETIN.EQ.2) THEN
        GOTO 904
      ELSEIF (IRETIN.NE.0) THEN
        GOTO 1001
      ELSEIF (IARTIC.EQ.0) THEN
        IREP=-16
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
      KLONG=LFI%MLGPOS(IXM(IARTIC,IRGPIM))
      CLNOMA=LFI%CNOMAR(IXC(IARTIC,IRGPIM))
C
C        Recherche de la longueur "utile" du nom d'article.
C        (c'est-a-dire sans tenir compte des blancs terminaux eventuels)
C
      IDECBL=0
C
  211 CONTINUE
      IPOSBL=IDECBL+INDEX (CLNOMA(IDECBL+1:),' ')
C
      IF (IPOSBL.LE.IDECBL) THEN
        ILCLNO=LFI%JPNCPN
      ELSEIF (CLNOMA(IPOSBL:).EQ.' ') THEN
        ILCLNO=IPOSBL-1
      ELSE
        IDECBL=IPOSBL
        GOTO 211
      ENDIF
C
      IF (ILCDNO.GE.ILCLNO) THEN
        CDNOMA=CLNOMA(:ILCLNO)
      ELSE
        IREP=-24
        CLACTI=CLNOMA
        GOTO 1001
      ENDIF
C**
C     3.  -  STOCKAGE DES PARAMETRES D'APPEL DANS LES TABLES.
C-----------------------------------------------------------------------
C
C           VERROUILLAGE GLOBAL EVENTUEL.
C
       IF (LFI%LMULTI) CALL LFIVER_MT (LFI, LFI%VERGLA,'ON')
      LLVERG=LFI%LMULTI
C
      IF (LFI%NUIMEX.LT.LFI%JPIMEX) THEN
C
        DO 301 J=1,LFI%JPIMEX
C
        IF (LFI%MNUIEX(J).EQ.LFI%JPNIL) THEN
          IRANIE=J
          LFI%NUIMEX=LFI%NUIMEX+1
          LFI%NINIEX(LFI%NUIMEX)=J
          LFI%MNUIEX(J)=KNUMER
          GOTO 302
        ENDIF
C
  301   CONTINUE
C
        IREP=-16
        GOTO 1001
C
      ELSE
C
C        Tables deja pleines...
C
        IREP=-37
        GOTO 1001
      ENDIF
C
  302   CONTINUE
C
C         Deverrouillage Global eventuel.
C
       IF (LFI%LMULTI) CALL LFIVER_MT (LFI, LFI%VERGLA,'OFF')
      LLVERG=.FALSE.
C
      LFI%NEXPOR(IRANG)=IRANIE
      LFI%NAEXPL(IRANIE)=0
      LFI%CNIMPL(IRANIE)=' '
      LFI%NIMPEX(IRANIE)=KNUMEX
      LFI%NUTRAV(IRANIE)=KNUTRA
      LFI%NLAPFD(IRANIE)=KLAREX*KFACEX
      LFI%NXCNLD(IRANIE)=KXCNEX
      LFI%NRCFMX(IRANIE)=IRANMX
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
      GOTO 909
C
  906 CONTINUE
      CLACTI='REWIND'
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
       IF (LLVERG) CALL LFIVER_MT (LFI, LFI%VERGLA,'OFF')
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
        IF (LHOOK) CALL DR_HOOK('LFIPXF_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='LFIPXF'
        WRITE (UNIT=CLMESS,FMT='(''ARGUMENTS='',I4,2('','',I3),A,
     S         '','',I5,2('','',I2),'','',I3,A,'','',I6)')
     S     KREP,KNUMER,KNUMEX,CLCFGX(:ILCFGX),KLAREX,KXCNEX,KFACEX,
     S     KNUTRA,CLNOMA(:ILCLNO),KLONG
      CALL LFIEMS_MT (LFI, INUMER,INIMES,IREP,LLFATA,
     S                CLMESS,CLNSPR,CLACTI)
C
      IF (LHOOK) CALL DR_HOOK('LFIPXF_MT',1,ZHOOK_HANDLE)
      END

