C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFIOUV_MT (LFI, KREP, KNUMER, LDNOMM, CDNOMF, 
     S                      CDSTTO, LDERFA,
     S                      LDIMST, KNIMES, KNBARP, KNBARI )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        SOUS-PROGRAMME D'OUVERTURE D'UNE UNITE LOGIQUE DEVANT ETRE
C     TRAITEE COMME UN FICHIER INDEXE, PAR LE LOGICIEL LFI.
C**
C     ARGUMENTS : KREP   (SORTIE) ==> CODE-REPONSE DU SOUS-PROGRAMME;
C                 KNUMER (ENTREE) ==> LFI%NUMERO DE L'UNITE LOGIQUE;
C                 LDNOMM (ENTREE) ==> VRAI SI L'UNITE LOGIQUE DOIT ETRE
C                                     ASSOCIEE A UN NOM DE FICHIER EXP-
C                                     LICITE LORS DE L'"OPEN" FORTRAN;
C                 CDNOMF (ENTREE) ==> NOM DE FICHIER EXPLICITE, SI
C                                     *LDNOMM* EST VRAI - MEME SI CE
C                                     N'EST PAS LE CAS, CE *DOIT* ETRE
C                                     UN OBJET DE TYPE "CHARACTER" .
C                 CDSTTO (ENTREE) ==> "STATUS" POUR L'"OPEN" FORTRAN
C                                     ('OLD','NEW','UNKNOWN','SCRATCH')
C                                     PAR DEFAUT, METTRE 'UNKNOWN';
C                 LDERFA (ENTREE) ==> OPTION D'ERREUR FATALE;
C                 LDIMST (ENTREE) ==> OPTION IMPRESSION DE STATISTIQUES
C                                     AU MOMENT DE LA FERMETURE;
C                 KNIMES (ENTREE) ==> NIVEAU DE LA MESSAGERIE (0,1 OU 2)
C                                     ( 0==>RIEN, 2==>TOUT )
C                 KNBARP (ENTREE) ==> NOMBRE D'ARTICLES LOGIQUES PREVUS,
C                                     CE QUI N'EST UTILISE QUE LORS DE
C                                     LA CREATION DU FICHIER,
C                                     ET QUI N'EMPECHE QUAND MEME PAS
C                                     D'AVOIR PLUS D'ARTICLES LOGIQUES;
C                 KNBARI (SORTIE) ==> NOMBRE D'ARTICLES LOGIQUES DE DON-
C                                     NEES SUR LE FICHIER, INITIALEMENT.
C                                     (ZERO SI CREATION)
      CHARACTER CPNOMD*(*)
      PARAMETER ( CPNOMD='%%%%% FICHIER SANS NOM %%%%%' )
C
C    Modifications:
C
C    02/06/97, Jean Clochard.
C
C              -Modification des impressions pour que l'annee puisse
C               etre imprimee avec 4 chiffres.
C
#ifndef f77
#include "precision.h"
#endif
C
      TYPE(LFICOM) :: LFI
      INTEGER KREP, KNUMER, KNIMES, KNBARP, KNBARI
#ifndef f77
      INTEGER (KIND=JPDBLE) IDATE, IHEURE
#else
      INTEGER IDATE, IHEURE
#endif
      INTEGER ILSTTU, IREPX, IRANG, IRANMS, INBARI, IDECBL, IPOSBL, J
      INTEGER ILNOMF, INLNOM, INIMES, IREP, ILNOMS, IFACTM, ILSTTO, IJ
      INTEGER IRANFM, ILACTI, ICOMPT, ITAILS, ICRITS, IPOFIN, ICRITG
      INTEGER ICRITD, ICRITR, IPOSCA, INREAD, INWRIT, IBASE, ILOREC
      INTEGER INAPHY, JREC, ILARPH, INALPP, IFACPH, IFACPP, INBPIR
      INTEGER IRANGD, IREC, INBALO, ILUTIL, IRGPIF, IRETIN
C
      LOGICAL LDNOMM, LDERFA, LDIMST, LLEXFI, LLNOUF, LLNOMS
      LOGICAL LLVERG, LLEXUL
C
      CHARACTER CDNOMF*(*), CDSTTO*(*)
      CHARACTER*(LFI%JPLSTX) CLSTTO
      CHARACTER*(LFI%JPLFTX) CLNOMF, CLNOMS
C
#include "lficom2.h"
#include "lficom_mt.h"
C
C     1.  -  CONTROLES DIVERS, ET INITIALISATIONS.
C-----------------------------------------------------------------------
C*
C     1.0 - PARTIE "ELEMENTAIRE".
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFIOUV_MT',0,ZHOOK_HANDLE)
      ILSTTU=MIN0 (LEN (CLSTTO), LEN (CDSTTO))
      IREPX=0
      IRANG=0
      IRANMS=0
      INBARI=0
      LLVERG=.FALSE.
C
C        Appel legerement anticipe a LFINUM, permettant une initialisa-
C     tion des variables globales du logiciel a la 1ere utilisation.
C
      CALL LFINUM_MT (LFI, KNUMER,IRANG)
C
      IF (LDNOMM) THEN
C
C        Recherche de la longueur "utile" du nom de fichier specifie.
C        (c'est-a-dire sans tenir compte des blancs terminaux eventuels)
C
        IDECBL=0
C
  101   CONTINUE
        IPOSBL=IDECBL+INDEX (CDNOMF(IDECBL+1:),' ')
C
        IF (IPOSBL.LE.IDECBL) THEN
          ILNOMF=LEN (CDNOMF)
        ELSEIF (CDNOMF(IPOSBL:).EQ.' ') THEN
          ILNOMF=MAX0 (IPOSBL-1,1)
        ELSE
          IDECBL=IPOSBL
          GOTO 101
        ENDIF
C
        IF (ILNOMF.GT.LFI%JPLFTX) THEN
          INLNOM=LFI%JPLFTX
          INIMES=LFI%NIMESG
C
          IF (INIMES.GE.1) THEN
C
C        Message preventif, car le controle de non ouverture d'un meme
C     fichier via deux unites logiques differentes risque de "sauter"
C     artificiellement... et pas forcement a cet appel.
C
C        Le code-reponse ci-dessous est bidon, mais permet de mettre
C     en relief le message via LFIEMS.
C
            IREP=LFI%JPNIL
            CLNSPR='LFIOUV'
C
            IF (LFI%LFRANC) THEN
              WRITE (UNIT=CLMESS,FMT=
     S               '(''ATTENTION: NOM DE FICHIER TRONQUE A'',I4,
     S                 '' CARACTERES...'')') LFI%JPLFTX
            ELSE
              WRITE (UNIT=CLMESS,FMT=
     S               '(''WARNING: FILE NAME TRUNCATED TO ONLY'',I4,
     S                 '' CHARACTERS...'')') LFI%JPLFTX
            ENDIF
C
            CALL LFIEMS_MT (LFI, KNUMER,INIMES,IREP,.FALSE.,
     S                      CLMESS,CLNSPR,
     S                      CLACTI)
          ENDIF
C
        ELSE
          INLNOM=ILNOMF
        ENDIF
C
        CLNOMF=CDNOMF(:INLNOM)
      ELSE
        ILNOMF=LEN (CPNOMD)
        CLNOMF=CPNOMD
        INLNOM=ILNOMF
      ENDIF
C
C       Ci-dessous, initialisations destinees a forcer l'impression
C     du nom du fichier en cas de problemes.
C
      CLNOMS=CLNOMF
      ILNOMS=INLNOM
      IFACTM=0
C
C        Controle de validite FORTRAN du Numero d'Unite Logique.
C
      INQUIRE (UNIT=KNUMER,EXIST=LLEXUL,ERR=901,IOSTAT=IREP)
C
      IF (.NOT.LLEXUL) THEN
        IREP=-30
        GOTO 1001
      ENDIF
C
C        CONTROLE DE L'ARGUMENT D'APPEL "KNIMES"
C
      IF (KNIMES.LT.0.OR.KNIMES.GT.2) THEN
        IREP=-2
        GOTO 1001
      ENDIF
C
C        CONTROLE DE L'ARGUMENT D'APPEL "CDSTTO"
C
      DO 103 J=1,LFI%JPNBST
      IF (CDSTTO.EQ.LFI%LFIOUV_CLSTEX(J)) GOTO 104
  103 CONTINUE
C
      ILACTI=MIN0 (LEN (CDSTTO),LEN (CLACTI))
      CLACTI=CDSTTO(:ILACTI)
      IREP=-7
      GOTO 1001
C
  104 CONTINUE
      ILSTTO=INDEX (CDSTTO,' ')-1
      IF (ILSTTO.GT.0) ILSTTU=ILSTTO
      CLSTTO=CDSTTO(:ILSTTU)
C
C               CONTROLE DE NON-OUVERTURE PREALABLE.
C
      IF (IRANG.NE.0) THEN
        IREP=-5
        GOTO 1001
      ENDIF
C
       IF (LFI%LMULTI) CALL LFIVER_MT (LFI, LFI%VERGLA,'ON')
      LLVERG=LFI%LMULTI
C
C        Recherche d'un eventuel facteur multiplicatif predefini pour
C     l'unite logique en question.
C
      CALL LFIFMP_MT (LFI, KNUMER,IRANFM)
      IFACTM=LFI%MFACTU(IRANFM)
C
      IF (LDNOMM) THEN
C
C        SI LE FICHIER EST NOMME, ON VERIFIE QU'IL N'A PAS ETE
C        DEJA OUVERT POUR UNE AUTRE UNITE LOGIQUE.
C
        DO 105 J=1,LFI%NBFIOU
        IJ=LFI%NUMIND(J)
C
        IF (CLNOMF.EQ.LFI%CNOMFI(IJ)(:MIN0 (LFI%JPLFTX,LFI%NLNOMF(IJ))))
     S  THEN
          ILACTI=MIN0(LEN (CLNOMF),LEN (CLACTI))
          CLACTI=CLNOMF(:ILACTI)
          IRANMS=IJ
          IREP=-13
          GOTO 1001
        ENDIF
C
  105   CONTINUE
C
      ENDIF
C
  110 CONTINUE
C*
C     1.1 - RECHERCHE D'UN EMPLACEMENT DISPONIBLE DANS LA TABLE DES
C           NUMEROS D'UNITES LOGIQUES *LFI%NUMERO* .
C           (Il faut IFACTM emplacements CONSECUTIFS)
C-----------------------------------------------------------------------
C
      IF ((LFI%NFACTM+IFACTM).GT.LFI%JPNXFI) THEN
        IREP=-6
        GOTO 1001
      ENDIF
C
      ICOMPT=0
      ITAILS=LFI%JPNXFI+1
      ICRITS=0
C
      DO 111 J=1,LFI%JPNXFI
C
      IF (LFI%NUMERO(J).EQ.LFI%JPNIL) THEN
        ICOMPT=ICOMPT+1
        IF (J.NE.LFI%JPNXFI.OR.ICOMPT.LT.IFACTM.OR.ICOMPT.GT.ITAILS)
     S    GOTO 111
        IPOFIN=LFI%JPNXFI
      ELSEIF (ICOMPT.LT.IFACTM.OR.ICOMPT.GT.ITAILS) THEN
        ICOMPT=0
C
        IF ((LFI%JPNXFI-J).LT.IFACTM) THEN
          GOTO 112
        ELSE
          GOTO 111
        ENDIF
C
      ELSE
        IPOFIN=J-1
      ENDIF
C
C       Les lignes qui suivent sont atteintes si on a trouve un espace
C     contigu suffisant dans la table LFI%NUMERO, et de taille inferieure
C     ou egale a ce qu'on aurait pu trouver precedemment.
C       On calcule alors un critere de cadrage (a gauche ou a droite)
C     dans cet espace, en privilegiant une occupation decentree.
C
      ICRITG=IABS (LFI%JPNXFI+1-2*(IPOFIN-ICOMPT+1))
      ICRITD=IABS (LFI%JPNXFI+1-2*IPOFIN)
C
      IF (ICRITG.GE.ICRITD) THEN
        ICRITR=ICRITG
        IPOSCA=IPOFIN-ICOMPT+1
      ELSE
        ICRITR=ICRITD
        IPOSCA=IPOFIN-IFACTM+1
      ENDIF
C
C       On retient l'espace trouve s'il est plus petit que ce qu'on
C     avait pu trouver precedemment, ou en cas d'egalite de taille
C     s'il est plus decentre.
C
      IF (ICOMPT.LT.ITAILS.OR.ICRITR.GT.ICRITS) THEN
        ITAILS=ICOMPT
        IRANG=IPOSCA
        ICRITS=ICRITR
      ENDIF
C
      ICOMPT=0
      IF ((LFI%JPNXFI-J).LT.IFACTM) GOTO 112
C
  111 CONTINUE
C
  112 CONTINUE
C
      IF (ITAILS.GT.LFI%JPNXFI) THEN
C
C         On n'a pas trouve d'espace ad hoc.
C
        IF (IFACTM.GT.1) THEN
          IREP=-27
        ELSE
          IREP=-16
        ENDIF
C
        GOTO 1001
C
      ENDIF
C
      IRANMS=IRANG
      IF (LFI%LMISOP) WRITE (UNIT=LFI%NULOUT,FMT=*)
     S  '====> LFIOUV - IRANG = ',IRANG, ', IFACTM = ',IFACTM
      LFI%LERFAT(IRANG)=LDERFA
      LFI%NIVMES(IRANG)=KNIMES
      INREAD=0
      INWRIT=0
C
C        CETTE INITIALISATION QUI PEUT PARAITRE BIEN COMPLIQUEE SERT
C     DE PARADE AU MAUVAIS COMPORTEMENT DU "READ" SUR UN FICHIER VIDE,
C     sur CRAY-2 sous UNICOS 4.0 et 5.0... ( Debut )
C
      CALL LFIDAH_MT (LFI, IDATE,IHEURE)
      IBASE=IHEURE+LFI%JPNIL
C
      DO 113 J=1,LFI%JPLDOC
      LFI%MDES1D(IXM(J,IRANG))=IBASE-J
  113 CONTINUE
C**
C     2.  -  OUVERTURE DU FICHIER AU SENS FORTRAN DU TERME (*OPEN*).
C-----------------------------------------------------------------------
C
      ILOREC=LFI%JPRECL*IFACTM
C
      IF (LDNOMM) THEN
C*
C     2.1 - CAS OU L'UNITE LOGIQUE DOIT ETRE ASSOCIEE A UN FICHIER
C           DONT LE NOM EST EXPLICITEMENT DONNE.
C-----------------------------------------------------------------------
C
        INQUIRE (FILE=CDNOMF,EXIST=LLEXFI,IOSTAT=IREP,ERR=901)
C
        IF (LLEXFI.AND.CLSTTO.EQ.'NEW'
     S      .OR..NOT.LLEXFI.AND.CLSTTO.EQ.'OLD') THEN
          CLACTI=CLSTTO
          IREP=-9
          IRANG=0
          IRANMS=0
          GOTO 1001
        ENDIF
C
        LLNOUF=CLSTTO.EQ.'NEW'.OR.CLSTTO.EQ.'SCRATCH'.OR..NOT.LLEXFI
C
C     APRES TOUS CES CONTROLES DE BASE, ON TENTE L'"OPEN" DU FICHIER .
C
        OPEN (UNIT=KNUMER,FILE=CDNOMF,STATUS=CLSTTO,ERR=902,
     S       FORM='UNFORMATTED',ACCESS='DIRECT',RECL=ILOREC,IOSTAT=IREP)
C
      ELSE
C*
C     2.2 - CAS OU L'UNITE LOGIQUE N'A PAS DE NOM DE FICHIER ASSOCIE
C           EXPLICITE; ON TENTE DIRECTEMENT L'"OPEN" .
C-----------------------------------------------------------------------
C
        IF (CLSTTO.NE.'OLD'.AND.CLSTTO.NE.'NEW') THEN
          OPEN (UNIT=KNUMER,STATUS=CLSTTO,FORM='UNFORMATTED',
     S          ACCESS='DIRECT',RECL=ILOREC,ERR=902,IOSTAT=IREP)
        ELSE
          OPEN (UNIT=KNUMER,FORM='UNFORMATTED',
     S          ACCESS='DIRECT',RECL=ILOREC,ERR=902,IOSTAT=IREP)
        ENDIF
C
        LLNOUF=CLSTTO.EQ.'SCRATCH'
C
      ENDIF
C*
C     2.3 - L'"OPEN" S'EST BIEN PASSE... ON ESSAIE DE RECUPERER LE NOM
C           *SYSTEME* EVENTUEL ASSOCIE A L'UNITE LOGIQUE.
C-----------------------------------------------------------------------
C
      INQUIRE (UNIT=KNUMER,NAMED=LLNOMS,NAME=CLNOMS,ERR=901,
     S         IOSTAT=IREP)
C
      IF (LLNOMS) THEN
C
C        Recherche de la longueur "utile" du nom systeme du fichier.
C        (c'est-a-dire sans tenir compte des blancs terminaux eventuels)
C
        IDECBL=0
C
  231   CONTINUE
        IPOSBL=IDECBL+INDEX (CLNOMS(IDECBL+1:),' ')
C
        IF (IPOSBL.LE.IDECBL) THEN
          ILNOMS=LEN (CLNOMS)
        ELSEIF (CLNOMS(IPOSBL:).EQ.' ') THEN
          ILNOMS=MAX0 (IPOSBL-1,1)
        ELSE
          IDECBL=IPOSBL
          GOTO 231
        ENDIF
C
        IF (.NOT.LDNOMM) THEN
          ILNOMF=ILNOMS
          INLNOM=ILNOMS
          CLNOMF=CLNOMS
        ENDIF
C
        DO 233 J=1,LFI%NBFIOU
        IJ=LFI%NUMIND(J)
C
        IF (CLNOMS.EQ.LFI%CNOMSY(IJ)(:LFI%NLNOMS(IJ))) THEN
          ILACTI=MIN0(LEN (CLNOMS),LEN (CLACTI))
          CLACTI=CLNOMS(:ILACTI)
          IREP=-13
          IRANG=0
          IRANMS=0
          GOTO 1001
        ENDIF
C
  233   CONTINUE
C
      ELSE
        ILNOMS=LEN (CPNOMD)
        CLNOMS=CPNOMD
      ENDIF
C
      IF (CLSTTO.EQ.'OLD'.OR..NOT.LLNOUF) THEN
C**
C     3.  -  DANS LE CAS OU LE FICHIER DEVAIT OU POUVAIT EXISTER AVANT
C            OUVERTURE, ON ESSAIE DE LIRE LES PREMIERS ARTICLES.
C-----------------------------------------------------------------------
C        ( L'ARTICLE DOCUMENTAIRE ET UNE PAIRE D'ARTICLES D'INDEX;
C          ON COMMENCE PAR L'ARTICLE NO. 3 CAR IL Y A PLUS DE CHANCES
C          D'AVOIR UNE MAUVAISE LECTURE POUR CELUI-CI )
C
C          DANS LE CAS DU "STATUS" 'UNKNOWN', IL S'AGIT DE LEVER
C       L'AMBIGUITE: FICHIER DEJA ECRIT PAR LE LOGICIEL, OU DEVANT ETRE
C       CREE PAR LUI ?
C
        DO 301 JREC=3,1,-2
        INAPHY=JREC
        CALL LFILDO_MT (LFI, IREP,KNUMER,JREC,LFI%MDES1D(IXM(1,IRANG)),
     S                  INREAD,IFACTM,IRETIN)
C
        IF (IRETIN.NE.0) THEN
          GOTO 302
        ENDIF
C
  301   CONTINUE
C
  302   CONTINUE
C
        IF (IREP.EQ.0) THEN
C
C           LECTURE OK... ON CONTROLE QUELQUES VALEURS "DOCUMENTAIRES"
C
C            Fin de la parade sur CRAY2, sous UNICOS 4.0 et 5.0 .
C
          DO 303 J=1,LFI%JPLDOC
          IF (LFI%MDES1D(IXM(J,IRANG)).NE.(IBASE-J)) GOTO 304
  303     CONTINUE
C
          LLNOUF=.TRUE.
          GOTO 390
C
  304     CONTINUE
          LLNOUF=.FALSE.
          ILARPH=LFI%MDES1D(IXM(LFI%JPLPAR,IRANG))
          INALPP=LFI%MDES1D(IXM(LFI%JPXAPI,IRANG))
          IFACPH=ILARPH/LFI%JPLARD
          IFACPP=INALPP/LFI%JPNAPP
C
          IF (MIN0 (ILARPH,INALPP).LE.0.OR.MOD (ILARPH,LFI%JPLARD).NE.0
     S        .OR.LFI%MDES1D(IXM(LFI%JPLMNA,IRANG)).NE.LFI%JPNCPN
     S        .OR.LFI%MDES1D(IXM(LFI%JPLLDO,IRANG)).NE.LFI%JPLDOC
     S        .OR.MOD (INALPP,LFI%JPNAPP).NE.0.OR.IFACPP.NE.IFACPH) THEN
            IREP=-10
            IRANG=0
            IRANMS=0
            GOTO 1001
          ELSEIF (LFI%MDES1D(IXM(LFI%JPFEAM,IRANG)).NE.0) THEN
            IREP=-11
            LLFATA=LLMOER (IREP,IRANG)
C
            IF (LLFATA) THEN
              IRANG=0
              IRANMS=0
              GOTO 1001
            ENDIF
C
C        SI L'ERREUR (-11) N'A PAS ETE FATALE, ON DONNE LA POSSIBILITE
C       DE TRAITER LE FICHIER DONT LA DERNIERE MODIFICATION N'A PAS ETE
C       "ENREGISTREE" . MAIS SANS AUCUNE GARANTIE ...
C
          ENDIF
C
          IF (IFACPH.NE.IFACTM) THEN
C
C     Messagerie de Niveau 1 pour prevenir de l'incident...
C
            INIMES=IXNIMS (IRANMS)
C
            IF (INIMES.GE.1) THEN
              CLNSPR='LFIOUV'
C
              IF (LFI%LFRANC) THEN
                WRITE (UNIT=CLMESS,FMT='(''Unite logique'',I3,
     S '', facteur multiplicatif lu sur fichier='',I3,'', attendu='',
     S                 I3)')KNUMER,IFACPH,IFACTM
              ELSE
                WRITE (UNIT=CLMESS,FMT='(''Logical Unit'',I3,
     S '', multiply factor read on file='',I3,'', expected='',I3)')
     S          KNUMER,IFACPH,IFACTM
              ENDIF
C
              IREPX=IREP
              IREP=0
              CALL LFIEMS_MT (LFI, KNUMER,INIMES,IREP,.FALSE.,
     S                        CLMESS,CLNSPR,CLACTI)
C
              IF (LFI%LFRANC) THEN
                ILUTIL=MIN0 (INLNOM,LFI%JPLFIX,LEN (CLMESS)-6)
                CLMESS='Nom='''//CLNOMF(:ILUTIL)//''''
              ELSE
                ILUTIL=MIN0 (INLNOM,LFI%JPLFIX,LEN (CLMESS)-7)
                CLMESS='Name='''//CLNOMF(:ILUTIL)//''''
              ENDIF
C
              CALL LFIEMS_MT (LFI, KNUMER,INIMES,IREP,.FALSE.,
     S                        CLMESS,CLNSPR,CLACTI)
C
              IF (LDNOMM.AND.CLNOMS.NE.CLNOMF) THEN
C
                IF (LFI%LFRANC) THEN
                  ILUTIL=MIN0 (ILNOMS,LFI%JPLFIX,LEN (CLMESS)-14)
                  CLMESS='Nom SYSTEME='''//CLNOMS(:ILUTIL)//''''
                ELSE
                  ILUTIL=MIN0 (ILNOMS,LFI%JPLFIX,LEN (CLMESS)-14)
                  CLMESS='SYSTEM Name='''//CLNOMS(:ILUTIL)//''''
                ENDIF
C
                CALL LFIEMS_MT (LFI, KNUMER,INIMES,IREP,.FALSE.,
     S                          CLMESS,CLNSPR,CLACTI)
              ENDIF
C
              IF (LFI%LFRANC) THEN
                CLMESS='On essaie de s''adapter au facteur '
     S               //'multiplicatif lu sur le fichier...'
              ELSE
                CLMESS='One tries to adapt to multiply '
     S               //'factor read on the file...'
              ENDIF
C
              CALL LFIEMS_MT (LFI, KNUMER,INIMES,IREP,.FALSE.,
     S                        CLMESS,CLNSPR,CLACTI)
              IREP=IREPX
            ENDIF
C
C        On va essayer de traiter le fichier avec la longueur d'Article
C     Physique lue sur le fichier. Pour cela, on doit d'abord le fermer,
C     puis on va recommencer le traitement depuis le paragraphe 1.1 .
C
            IRANG=0
            IRANMS=0
            CLOSE (UNIT=KNUMER,STATUS='KEEP',ERR=905,IOSTAT=IREP)
C
            IF (IFACPH.GT.LFI%JPFACX) THEN
              IREP=-28
              GOTO 1001
            ENDIF
C
            IFACTM=IFACPH
            GOTO 110
          ENDIF
C
        ELSEIF (CLSTTO.EQ.'OLD') THEN
          IREP=-12
          IRANG=0
          IRANMS=0
          GOTO 1001
        ELSE
          IREP=0
          LLNOUF=.TRUE.
        ENDIF
C
      ENDIF
C
  390 CONTINUE
C
C         Controle ultime avant le paragraphe suivant, dans la mesure
C     ou, contrairement au FORTRAN, on autorise les "STATUS" 'OLD'
C     et 'NEW' pour une unite logique sans nom de fichier explicite...
C     puisque le logiciel a sa propre mecanique de discrimination entre
C     un fichier "existant" ou "en mode creation".
C
      IF (LLNOUF.AND.CLSTTO.EQ.'OLD'
     S    .OR..NOT.LLNOUF.AND.CLSTTO.EQ.'NEW') THEN
        CLACTI=CLSTTO
        IREP=-9
        IRANG=0
        IRANMS=0
        GOTO 1001
      ENDIF
C**
C     4.  -  L'OUVERTURE FORTRAN EST OK, ON SAIT SI ON EST EN MODE
C            CREATION DU FICHIER INDEXE OU NON... ON COMMENCE A GARNIR
C            LES VARIABLES EN COMMON, MAIS SANS INCREMENTER *LFI%NBFIOU*
C            CAR ON PEUT ENCORE AVOIR DE (MAUVAISES) SURPRISES.
C-----------------------------------------------------------------------
C
      IREPX=IREP
      LFI%CNOMFI(IRANG)=CLNOMF
      LFI%NLNOMF(IRANG)=ILNOMF
      LFI%CNOMSY(IRANG)=CLNOMS
      LFI%NLNOMS(IRANG)=ILNOMS
      LFI%NDEROP(IRANG)=0
      LFI%CSTAOP(IRANG)=CLSTTO
      LFI%LNOUFI(IRANG)=LLNOUF
      LFI%LMODIF(IRANG)=.FALSE.
      LFI%NDERCO(IRANG)=IREP
      LFI%NTRULZ(IRANG)=0
      LFI%NRFPTZ(IRANG)=0
      LFI%NRFDTZ(IRANG)=0
      LFI%NBMOLU(IRANG)=0
      LFI%NBMOEC(IRANG)=0
      LFI%NDERGF(IRANG)=LFI%JPNIL
      LFI%CNDERA(IRANG)=' '
      LFI%MFACTM(IRANG)=IFACTM
      LFI%NSUIVF(IRANG)=LFI%JPNIL
      LFI%NPRECF(IRANG)=LFI%JPNIL
C
C     N.B.: LES PAGES D'INDEX DE RANG "IRANG" SONT AUTOMATIQUEMENT
C        "AFFECTEES" A L'UNITE LOGIQUE AYANT CE RANG, ET SERVENT
C        A Y STOCKER LA PREMIERE P.A.I. EN RANG DANS LE FICHIER.
C
C     ( LES PAGES D'INDEX DE RANG "IRANG+(J-1)*LFI%JPNXFI" OU J VARIE
C       DE 1 A LFI%JPNPIA, SONT AUTOMATIQUEMENT AFFECTEES A L'UNITE
C       LOGIQUE DE RANG "IRANG" )
C
      LFI%NBLECT(IRANG)=0
      LFI%NBNECR(IRANG)=0
      LFI%NREESP(IRANG)=0
      LFI%NREECO(IRANG)=0
      LFI%NREELO(IRANG)=0
      LFI%NBTROU(IRANG)=0
      LFI%NBRENO(IRANG)=0
      LFI%NBSUPP(IRANG)=0
      LFI%LISTAT(IRANG)=LDIMST
      LFI%LMIMAL(IRANG)=.FALSE.
       IF (LFI%LMULTI) CALL LFIVER_MT (LFI, LFI%VERRUE(IRANG),'ASGN')
C
      IF (LLNOUF) THEN
C*
C     4.1 - CAS DE CREATION DU FICHIER INDEXE - INITIALISATIONS DIVERSES
C-----------------------------------------------------------------------
C
        ILARPH=LFI%JPLARD*IFACTM
        INALPP=LFI%JPNAPP*IFACTM
C
        DO 412 J=1,ILARPH
        LFI%MLGPOS(IXM(J,IRANG))=0
  412   CONTINUE
C
        DO 413 J=1,LFI%JPNXNA
        LFI%CNOMAR(IXC(J,IRANG))=' '
  413   CONTINUE
C
        DO 414 J=1,ILARPH
        LFI%MDES1D(IXM(J,IRANG))=0
  414   CONTINUE
C
C         NOMBRE DE PAIRES D'ARTICLES D'INDEX RESERVES,
C         (ELLES OCCUPERONT LES ARTICLES 2 A (2*INBPIR+1) DU FICHIER)
C         ET REMPLISSAGE DE CERTAINS MOTS DE L'ARTICLE DOCUMENTAIRE.
C
        INBPIR=MAX0 (1,MIN0 (LFI%JPNXPR,1+(KNBARP-1)/INALPP))
        LFI%MDES1D(IXM(LFI%JPNPIR,IRANG))=INBPIR
        LFI%MDES1D(IXM(LFI%JPNAPH,IRANG))=1+2*INBPIR
        LFI%MDES1D(IXM(LFI%JPLPAR,IRANG))=ILARPH
        LFI%MDES1D(IXM(LFI%JPLMNA,IRANG))=LFI%JPNCPN
        LFI%MDES1D(IXM(LFI%JPLLDO,IRANG))=LFI%JPLDOC
        LFI%MDES1D(IXM(LFI%JPXAPI,IRANG))=INALPP
        LFI%MDES1D(IXM(LFI%JPFEAM,IRANG))=1
        LFI%NPODPI(IRANG)=1
        LFI%NALDPI(IRANG)=0
        LFI%NPPIMM(IRANG)=1
        IRANGD=IRANG
        CALL LFIDAH_MT (LFI, LFI%MDES1D(IXM(LFI%JPDCRE,IRANG)),
     S               LFI%MDES1D(IXM(LFI%JPHCRE,IRANG)))
C
C          ECRITURE DU PREMIER ARTICLE (DESCRIPTIF)
C
        IREC=1
        INAPHY=IREC
        CALL LFIEDO_MT (LFI, IREP,KNUMER,IREC,LFI%MDES1D(IXM(1,IRANG)),
     S                  INWRIT,IFACTM,IRETIN)
C
        IF (IRETIN.NE.0) THEN
          GOTO 904
        ENDIF
C
C
C     Remise a zero du descripteur en vue d'une fermeture normale.
C
        LFI%MDES1D(IXM(LFI%JPFEAM,IRANG))=0
C
C          ECRITURE DES ARTICLES CONTENANT LES PAIRES D'ARTICLES D'INDEX
C          "RESERVES".
C
        DO 415 J=1,INBPIR
        IREC=IREC+1
        INAPHY=IREC
        CALL LFIECC_MT (LFI, IREP,KNUMER,IREC,LFI%CNOMAR(IXC(1,IRANG)),
     S                  INWRIT,IFACTM,IRETIN)
C
        IF (IRETIN.NE.0) THEN
          GOTO 903
        ENDIF
C
        IREC=IREC+1
        INAPHY=IREC
        CALL LFIEDO_MT (LFI, IREP,KNUMER,IREC,LFI%MLGPOS(IXM(1,IRANG)),
     S                  INWRIT,IFACTM,IRETIN)
C
        IF (IRETIN.NE.0) THEN
          GOTO 904
        ENDIF
C
  415   CONTINUE
C
      ELSE
C*
C     4.2 - LE FICHIER EXISTAIT DEJA... ON LIT LA 1ERE PAIRE D'ARTICLES
C           D'INDEX ( + LA DERNIERE S'IL Y EN A AU MOINS 2 *UTILISEES* )
C-----------------------------------------------------------------------
C
        INBALO=LFI%MDES1D(IXM(LFI%JPNALO,IRANG))
        INBPIR=LFI%MDES1D(IXM(LFI%JPNPIR,IRANG))
        IREC=2
        INAPHY=IREC
        CALL LFILCC_MT (LFI, IREP,KNUMER,IREC,LFI%CNOMAR(IXC(1,IRANG)),
     S                  INREAD,IFACTM,IRETIN)
C
        IF (IRETIN.NE.0) THEN
          GOTO 904
        ENDIF
C
        IREC=3
        INAPHY=IREC
        CALL LFILDO_MT (LFI, IREP,KNUMER,IREC,LFI%MLGPOS(IXM(1,IRANG)),
     S                  INREAD,IFACTM,IRETIN)
C
        IF (IRETIN.NE.0) THEN
          GOTO 904
        ENDIF
C
        IF (INBALO.LE.INALPP) THEN
          LFI%NALDPI(IRANG)=INBALO
          LFI%NPODPI(IRANG)=1
          LFI%NPPIMM(IRANG)=1
          IRANGD=IRANG
        ELSE
C
C          CAS OU IL Y A AU MOINS 2 PAIRES D'ARTICLES D'INDEX UTILISEES.
C
          IRGPIF=1+(INBALO-1)/INALPP
          CALL LFIREC_MT (LFI, IRGPIF,IRANG,IREC)
          IRANGD=IRANG+LFI%JPNXFI
          INAPHY=IREC
          CALL LFILCC_MT (LFI, IREP,KNUMER,IREC,
     S                    LFI%CNOMAR(IXC(1,IRANGD)),
     S                    INREAD,IFACTM,IRETIN)
C
          IF (IRETIN.NE.0) THEN
            GOTO 904
          ENDIF
C
          IREC=IREC+1
          INAPHY=IREC
          CALL LFILDO_MT (LFI, IREP,KNUMER,IREC,
     S                    LFI%MLGPOS(IXM(1,IRANGD)),
     S                    INREAD,IFACTM,IRETIN)
C
          IF (IRETIN.NE.0) THEN
            GOTO 904
          ENDIF
C
          LFI%NALDPI(IRANG)=1+MOD (INBALO-1,INALPP)
          LFI%NPODPI(IRANG)=2
          LFI%NPPIMM(IRANG)=2
          LFI%MRGPIM(2,IRANG)=IRANGD
          LFI%MRGPIF(IRANGD)=IRGPIF
        ENDIF
C
      ENDIF
C**
C     5.  -  L'OUVERTURE AU SENS DU LOGICIEL DE FICHIERS INDEXES LFI
C            EST COMPLETE; ON MET DONC A JOUR LES DERNIERES VARIABLES
C            EN COMMON, DONT *LFI%NBFIOU*.
C-----------------------------------------------------------------------
C
C           REMARQUE: LA PREMIERE ET LA DERNIERE P.P.I. SONT TOUJOURS
C                     "PHASEES".
C
      DO 501 J=IRANG,IRANGD,LFI%JPNXFI
      LFI%LECRPI(J,1)=.FALSE.
      LFI%LECRPI(J,2)=.FALSE.
      LFI%LPHASP(J)=.TRUE.
  501 CONTINUE
C
      DO 502 J=0,LFI%JPNPDF-1
      LFI%NUMAPD(J,IRANG)=LFI%JPNIL
      LFI%NLONPD(J,IRANG)=0
      LFI%LECRPD(J,IRANG)=.FALSE.
  502 CONTINUE
C
      DO 503 J=1,IFACTM
      LFI%NUMERO(IRANG+J-1)=KNUMER
  503 CONTINUE
C
      LFI%NDERPD(IRANG)=LFI%JPNPDF-1
      LFI%NBFIOU=LFI%NBFIOU+1
      LFI%NFACTM=LFI%NFACTM+IFACTM
      LFI%NUMIND(LFI%NBFIOU)=IRANG
      INBARI=LFI%MDES1D(IXM(LFI%JPNALO,IRANG))-
     S       LFI%MDES1D(IXM(LFI%JPNTRU,IRANG))
      LFI%NBREAD(IRANG)=INREAD
      LFI%NBWRIT(IRANG)=INWRIT
      LFI%LTAMPL(IRANG)=LFI%LTAMLG
      LFI%LTAMPE(IRANG)=LFI%LTAMEG
      LFI%NEXPOR(IRANG)=LFI%JPNIL
      LFI%NIMPOR(IRANG)=LFI%JPNIL
C
      IREP=IREPX
      GOTO 1001
C**
C     9.  - CI-DESSOUS, ETIQUETTES DE BRANCHEMENT EN CAS D'ERREUR E/S.
C-----------------------------------------------------------------------
C
  901 CONTINUE
      CLACTI='INQUIRE'
      INAPHY=0
      GOTO 909
C
  902 CONTINUE
      CLACTI='OPEN'
      IRANG=0
      IRANMS=0
      INAPHY=0
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
  905 CONTINUE
      CLACTI='CLOSE'
      INAPHY=0
C
  909 CONTINUE
C
C      AU CAS OU, ON FORCE LE CODE-REPONSE ENTREE/SORTIE A ETRE POSITIF.
C
      IREP=IABS (IREP)
      LFI%NUMAPH(IRANG)=INAPHY
      IF (IRANG.EQ.0) LFI%MFACTM(0)=IFACTM
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE SOUS-PROGRAMME "LFIEMS" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
      KREP=IREP
      KNBARI=INBARI
      LLFATA=LLMOER (IREP,IRANG)
C
      IF (LLFATA) THEN
        INIMES=2
      ELSE
        INIMES=IXNIMS (IRANMS)
      ENDIF
C
       IF (LFI%LMULTI.AND.LLVERG) CALL LFIVER_MT (LFI, LFI%VERGLA,'OFF')
C
      IF (.NOT.LLFATA.AND.INIMES.EQ.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('LFIOUV_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='LFIOUV'
C
      IF (INIMES.GE.1) THEN
C
C           Impression du nom du fichier.
C
        IF (LFI%LFRANC) THEN
          ILUTIL=MIN0 (INLNOM,LFI%JPLFIX,LEN (CLMESS)-6)
          CLMESS='Nom='''//CLNOMF(:ILUTIL)//''''
        ELSE
          ILUTIL=MIN0 (INLNOM,LFI%JPLFIX,LEN (CLMESS)-7)
          CLMESS='Name='''//CLNOMF(:ILUTIL)//''''
        ENDIF
C
        CALL LFIEMS_MT (LFI, KNUMER,INIMES,IREP,.FALSE.,
     S                  CLMESS,CLNSPR,CLACTI)
C
        IF (LDNOMM.AND.CLNOMS(:ILNOMS).NE.CLNOMF(:INLNOM)) THEN
C
          IF (LFI%LFRANC) THEN
            ILUTIL=MIN0 (ILNOMS,LFI%JPLFIX,LEN (CLMESS)-14)
            CLMESS='Nom SYSTEME='''//CLNOMS(:ILUTIL)//''''
          ELSE
            ILUTIL=MIN0 (ILNOMS,LFI%JPLFIX,LEN (CLMESS)-14)
            CLMESS='SYSTEM Name='''//CLNOMS(:ILUTIL)//''''
          ENDIF
C
          CALL LFIEMS_MT (LFI, KNUMER,INIMES,IREP,.FALSE.,
     S                    CLMESS,CLNSPR,CLACTI)
        ENDIF
C
      ENDIF
C
      IF (INIMES.EQ.2) THEN
        WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,
     S     '', LDNOMM= '',L1,'', CDSTTO='''''',A,'''''', LDERFA= '',L1,
     S     '',  LDIMST= '',L1,
     S         '', KNIMES='',I2,'', KNBARP='',I6,'' KNBARI='',I6)')
     S   KREP,KNUMER,LDNOMM,CDSTTO(:ILSTTU),LDERFA,LDIMST,KNIMES,KNBARP,
     S   KNBARI
        CALL LFIEMS_MT (LFI, KNUMER,INIMES,IREP,LLFATA,
     S                  CLMESS,CLNSPR,CLACTI)
      ENDIF
C
C        LA MESSAGERIE QUI SUIT N'EST PAS EMISE EN CAS D'ERREUR FATALE.
C
      IF (INIMES.GE.1.AND.(IREP.EQ.0.OR.IREP.EQ.-11)) THEN
C
        IF (LLNOUF) THEN
C
          IF (LFI%LFRANC) THEN
            WRITE (UNIT=CLMESS,FMT='(''Unite'',I3,
     S   '' OUVERTE, CREATION de Fichier,'',I7,'' Articles prevus,'',I7,
     S   '' Articles gerables sans debordement'')')
     S      KNUMER,KNBARP,INALPP*INBPIR
          ELSE
            WRITE (UNIT=CLMESS,FMT='(''Unit'',I3,
     S   '' OPENED, File CREATION,'',I7,'' expected Records,'',I7,
     S   '' Records may be handled without overflow'')')
     S      KNUMER,KNBARP,INALPP*INBPIR
          ENDIF
C
        ELSE
C
          IF (LFI%LFRANC) THEN
            WRITE (UNIT=CLMESS,FMT='(''Unite'',I3,
     S       '' OUVERTE, derniere Modification OK a'',I9.6,''_'',I6.6,
     S       '','',I7,'' Articles de donnees,'',I9,'' mots en tout'')')
     S     KNUMER,LFI%MDES1D(IXM(LFI%JPDDMG,IRANG)),
     S     LFI%MDES1D(IXM(LFI%JPHDMG,IRANG)),
     S            KNBARI,ILARPH*LFI%MDES1D(IXM(LFI%JPNAPH,IRANG))
          ELSE
            WRITE (UNIT=CLMESS,FMT='(''Unit'',I3,
     S       '' OPENED, last Modification OK at'',I9.6,''_'',I6.6,
     S       '','',I7,'' data Records,'',I9,'' words in file'')')
     S     KNUMER,LFI%MDES1D(IXM(LFI%JPDDMG,IRANG)),
     S     LFI%MDES1D(IXM(LFI%JPHDMG,IRANG)),
     S            KNBARI,ILARPH*LFI%MDES1D(IXM(LFI%JPNAPH,IRANG))
          ENDIF
C
        ENDIF
C
        CALL LFIEMS_MT (LFI, KNUMER,INIMES,IREP,.FALSE.,
     S                  CLMESS,CLNSPR,CLACTI)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('LFIOUV_MT',1,ZHOOK_HANDLE)
      END

