C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFIEFR_MT (LFI, KNUMER, KNIMES, KCODE, LDFATA, 
     &                      CDMESS, CDNSPR, CDACTI )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        CE SOUS-PROGRAMME EST CHARGE DE FAIRE L'IMPRESSION DES MESSAGES
C     STANDARD EMIS PAR LE LOGICIEL DE FICHIERS INDEXES LFI, EN FAISANT
C     SI BESOIN EST L'"ABORT" DU PROGRAMME .
C        Les messages lies au mode "mise au point" sont emis directement
C     par les sous-programmes concernes.
C
C        Ce sous-programme est la V.O. (Version Originale, francaise),
C     et est appele par le sous-programme "chapeau" LFIEMS.
C     Pour la version anglaise, voir LFIENG.
C     ( For english version see subroutine LFIENG )
C**
C        ARGUMENTS : KNUMER ==> Numero eventuel de l'Unite Logique;
C        ( tous                 ( si LFI%JPNIL ==> pas d'Unite Logique )
C         d'Entree ) KNIMES ==> Niveau (0,1,2) du Message;
C                    KCODE  ==> Code correspondant a l'action en cause;
C                    LDFATA ==> Vrai si on doit avorter le programme;
C                    CDMESS ==> Si KNIMES#0, Message a emettre;
C                    CDNSPR ==> Nom du sous-programme appelant LFIEMS;
C                    CDACTI ==> Nom de l'action d'entree/sortie FORTRAN
C                               si KCODE >0), sinon fourre-tout (!) .
C*
C     !----------------------------------------------------------------!
C     ! TABLE DES VALEURS POSSIBLES DES CODES-REPONSES DU LOGICIEL LFI !
C     !----------------------------------------------------------------!
C
C-----------------------------------------------------------------------
C      0 ==> Aucune erreur n'a ete detectee, tout est OK.
C-----------------------------------------------------------------------
C valeur ==> Il s'agit (de la valeur absolue) du code-reponse FORTRAN
C positive   d'une instruction OPEN, READ, WRITE, CLOSE ou INQUIRE; pour
C            le sens exact voir le manuel de reference du constructeur.
C-----------------------------------------------------------------------
C     -1 ==> Unite Logique non ouverte pour le logiciel.
C-----------------------------------------------------------------------
C     -2 ==> Valeur d'un "NIVEAU" hors plage [0-2] .
C-----------------------------------------------------------------------
C     -3 ==> Option de verrou erronee (s/p a usage interne "LFIVER") .
C-----------------------------------------------------------------------
C     -4 ==> Changement explicite de mode Multi-Taches avec au moins une
C            unite logique ouverte-risque de problemes (s/p "LFIINI") .
C-----------------------------------------------------------------------
C     -5 ==> Unite Logique deja ouverte (LFIOUV, LFIAFM, LFISFM) .
C-----------------------------------------------------------------------
C     -6 ==> Pas assez de place dans les tables pour ouvrir l'Unite
C            Logique demandee (LFIOUV) .
C-----------------------------------------------------------------------
C     -7 ==> Argument illicite de "STATUS" pour l'instruction FORTRAN
C            "OPEN" (LFIOUV) .
C-----------------------------------------------------------------------
C     -8 ==> Incompatibilite entre "LDNOMM" et "CDSTTO" (LFIOUV) :
C            un fichier de "STATUS" 'OLD' ou 'NEW' doit etre nomme .
C            (CE CODE-REPONSE N'A PLUS DE SENS ACTUELLEMENT)
C-----------------------------------------------------------------------
C     -9 ==> Incompatibilite entre le "STATUS" 'NEW' ou 'OLD' et (respe-
C            ctivement) l'existence ou non du fichier (LFIOUV) .
C-----------------------------------------------------------------------
C    -10 ==> Le fichier considere n'est pas un fichier de type LFI, ou
C            ne peut pas etre traite par cette version du logiciel.
C            (LFIOUV)
C-----------------------------------------------------------------------
C    -11 ==> Fichier non ferme apres une modification (LFIOUV): cette
C            erreur n'est pas fatale si "LDERFA" est .FALSE., mais alors
C            integrite et coherence des donnees ne sont pas garanties.
C            Noter qu'une fois qu'un fichier a ce type de probleme, ce
C            code-reponse restera meme apres modification ulterieure.
C-----------------------------------------------------------------------
C    -12 ==> Fichier de "STATUS" 'OLD' mais erreur sur la lecture du
C            premier article physique du fichier (LFIOUV) .
C-----------------------------------------------------------------------
C    -13 ==> Fichier deja ouvert pour une autre unite logique LFI.
C            (LFIOUV)
C-----------------------------------------------------------------------
C    -14 ==> Argument d'appel de type ENTIER incorrect (souvent negatif)
C-----------------------------------------------------------------------
C    -15 ==> Argument d'appel de type CARACTERE incorrect (longueur).
C-----------------------------------------------------------------------
C    -16 ==> Incoherence Tables, Fichier, appels s/p internes, logiciel.
C            CETTE ERREUR NE PEUT PAS ETRE FILTREE. EST TOUJOURS FATALE.
C-----------------------------------------------------------------------
C    -17 ==> Trop d'articles logiques sur le fichier pour un de plus.
C            (par articles logiques on entend ceux lisibles par l'utili-
C             sateur, mais aussi les trous reperes dans l'index... qui
C             sont crees lors de reecritures d'articles de donnees ne
C             pouvant se faire sur place, et lors de suppression d'arti-
C             cles; ces trous peuvent etre "recycles" - LFIECR)
C-----------------------------------------------------------------------
C    -18 ==> Nom d'Article logique compose uniquement de BLANCS illicite
C            (pour le fonctionnement interne du logiciel LFI,
C             les trous d'index sont reperes par un nom d'article blanc)
C-----------------------------------------------------------------------
C    -19 ==> Un fichier ouvert avec le "STATUS" 'SCRATCH' ne peut pas
C            etre conserve: "CDSTTC" a 'KEEP' est illicite (LFIFER) .
C            si cette erreur n'est pas fatale, alors on execute un
C            "CLOSE" FORTRAN sans parametre "STATUS", de la meme maniere
C            que lorsque "CDSTTC" n'est ni a 'KEEP' ni a 'DELETE'.
C-----------------------------------------------------------------------
C    -20 ==> L'article logique demande n'existe pas dans le fichier.
C            (LFILEC, LFIREN, LFISUP)
C-----------------------------------------------------------------------
C    -21 ==> L'article logique demande est PLUS LONG sur le fichier;
C            si cette erreur n'est pas fatale, le resultat est une
C            lecture PARTIELLE de l'article, a la longueur demandee.
C            (LFILAP, LFILAS, LFILEC)
C-----------------------------------------------------------------------
C    -22 ==> L'article logique demande est PLUS COURT sur le fichier;
C            meme si cette erreur n'est pas fatale, AUCUNE LECTURE DE
C            DONNEES N'EST FAITE (LFILAP, LFILAS, LFILEC) .
C-----------------------------------------------------------------------
C    -23 ==> Il n'y a pas ou plus d'article "SUIVANT" a lire (LFILAS) .
C-----------------------------------------------------------------------
C    -24 ==> La variable caractere donnee en argument d'appel de sortie
C            est TROP COURTE pour y stocker le NOM de l'article, meme en
C            supprimant d'eventuels caracteres blancs en fin de nom.
C            (LFICAP, LFICAS, LFILAP, LFILAS)
C-----------------------------------------------------------------------
C    -25 ==> Le nouveau nom de l'article logique est (deja) celui d'un
C            autre article logique du fichier (LFIREN).
C-----------------------------------------------------------------------
C    -26 ==> Il n'y a pas ou plus d'article "PRECEDENT" a lire (LFILAP).
C-----------------------------------------------------------------------
C    -27 ==> Espace CONTIGU insuffisant dans les tables pour gerer le
C            fichier "multiple" demande (LFIOUV) .
C-----------------------------------------------------------------------
C    -28 ==> Facteur multiplicatif (de la longueur d'article physique
C            elementaire) trop grand pour la configuration du logiciel.
C            (LFIOUV, LFIAFM, LFIFMD)
C-----------------------------------------------------------------------
C    -29 ==> Pas assez de place dans les tables pour definir le facteur
C            multiplicatif a associer a l'Unite Logique (LFIAFM) .
C-----------------------------------------------------------------------
C    -30 ==> Numero d'Unite Logique FORTRAN illicite.
C-----------------------------------------------------------------------
C    -31 ==> Numero d'Unite Logique sans facteur multiplicatif predefini
C            (LFISFM)
C-----------------------------------------------------------------------
C
#ifndef f77
      USE SDL_MOD   , ONLY : SDL_SRLABORT
#endif
#ifndef f77
#include "precision.h"
#endif
C
      TYPE(LFICOM) :: LFI
      INTEGER KNUMER, KNIMES, KCODE, ILDMES, ILBLAN, INLNOM, INUMER
      INTEGER IDECBL, IPOSBL, ILACTI, ILACT2, ILNSPR, ILMESU, IJL, J, IJ
      INTEGER INBALO, ILMESA, INLIGN, IDECAL
C
      LOGICAL LDFATA
C
      CHARACTER  CDNSPR*(*), CLJOLI*6, CDMESS*(*), CLMESA*80, CDACTI*(*)
C
#include "lficom2.h"
#include "lficom_mt.h"
C**
C     1.  -  INITIALISATIONS.
C-----------------------------------------------------------------------
C
C        Recherche de la longueur "utile" de l'argument CDACTI.
C        (c'est-a-dire sans tenir compte des blancs terminaux eventuels)
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFIEFR_MT',0,ZHOOK_HANDLE)
      IDECBL=0
C
  101 CONTINUE
      IPOSBL=IDECBL+INDEX (CDACTI(IDECBL+1:),' ')
C
      IF (IPOSBL.LE.IDECBL) THEN
        ILACTI=LEN (CDACTI)
      ELSEIF (CDACTI(IPOSBL:).EQ.' ') THEN
        ILACTI=IPOSBL-1
      ELSE
        IDECBL=IPOSBL
        GOTO 101
      ENDIF
C
      ILACT2=MIN0 (ILACTI,LFI%JPNCPN)
      ILACTI=MIN0 (ILACT2,8)
      ILNSPR=MIN0 (LEN (CDNSPR),LFI%JPLSPX)
C
C        Prefixe (et eventuellement suffixe) pour le(s) message(s).
C
      IF (LDFATA) THEN
        CLJOLI=' *****'
      ELSEIF (KNIMES.EQ.0.OR.KCODE.NE.0) THEN
        CLJOLI=' */*/*'
      ELSE
        CLJOLI=' /////'
      ENDIF
C
      IF (KNIMES.NE.0) THEN
C**
C     2.  -  ON IMPRIME LE MESSAGE PREPARE PAR LE S/P APPELLANT LFIEMS.
C-----------------------------------------------------------------------
C
        ILMESU=MIN0 (LEN (CLMESS)-LEN (CLJOLI)-ILNSPR-4,
     S               LEN (CDMESS))
        CLMESS=CLJOLI//' '//CDNSPR(1:ILNSPR)//' - '//
     S             CDMESS(1:ILMESU)
        WRITE (UNIT=LFI%NULOUT,FMT='(A)') CLMESS
      ENDIF
C
      IF (KNIMES.EQ.0.OR.LDFATA) THEN
C**
C     3.  -  CONSTITUTION D'UN MESSAGE "AD HOC", EN FONCTION DE *KCODE*.
C-----------------------------------------------------------------------
C
C     En preambule, on cherche si l'unite logique concernee correspond
C     ou non a une unite logique ouverte pour le logiciel LFI.
C
        IF (KNUMER.EQ.LFI%JPNIL) THEN
          IJL=0
        ELSE
C
          DO 301 J=1,LFI%NBFIOU
          IJL=LFI%NUMIND(J)
          IF (KNUMER.EQ.LFI%NUMERO(IJL)) GOTO 302
  301     CONTINUE
C
          IJL=0
        ENDIF
C
  302   CONTINUE
C
        IF (KCODE.GT.0) THEN
C
          IF ((CDACTI.EQ.'READ'.OR.CDACTI.EQ.'WRITE')
     S        .AND.LFI%NUMAPH(IJL).GT.0) THEN
            WRITE (UNIT=CLMESS,FMT='(''ERREUR "'',A,''"'',I7,
     S             '',UNITE'',I3,'',NUM.ART'',I6,'',*'',I6,
     S             '' MOTS'')') CDACTI(1:ILACTI),KCODE,KNUMER,
     S                           LFI%NUMAPH(IJL),
     S                           LFI%JPLARD*LFI%MFACTM(IJL)
          ELSE
            WRITE (UNIT=CLMESS,
     S             FMT='(''ERREUR "'',A,''" FORTRAN, CODE=''
     S             ,I7,'', UNITE='',I3)') CDACTI(1:ILACTI),KCODE,KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-1) THEN
          WRITE (UNIT=CLMESS,FMT='(''UNITE LOGIQUE'',I3,
     S           '' NON OUVERTE POUR LE LOGICIEL LFI'')') KNUMER
C
        ELSEIF (KCODE.EQ.-2) THEN
C
          IF (KNUMER.EQ.LFI%JPNIL) THEN
            CLMESS='PARAMETRE DE NIVEAU "KNIVAU" HORS PLAGE [0-2]'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S   '(''NIVEAU DE MESSAGERIE HORS PLAGE [0-2], UNITE'',I3)') KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-3) THEN
          ILDMES=MIN0 (8,LEN (CDMESS))
          CLMESS='ACTION '''//CDMESS(1:ILDMES)
     S           //''' INCONNUE SUR LES VERROUS'
C
        ELSEIF (KCODE.EQ.-4) THEN
          CLMESS='CHANGEMENT MODE MULTI-TACHES AVEC '
     S               //'UNITE(S) OUVERTE(S)'
C
        ELSEIF (KCODE.EQ.-5) THEN
          WRITE (UNIT=CLMESS,FMT='(''UNITE LOGIQUE'',I3,
     S           '' DEJA OUVERTE POUR LFI - NE DEVRAIT PAS.'')') KNUMER
C
        ELSEIF (KCODE.EQ.-6) THEN
          WRITE (UNIT=CLMESS,FMT='(I3,'' ENTREES,'',
     S    '' PLUS ASSEZ DE PLACE DANS LES TABLES, UNITE'',I3)')
     S    LFI%JPNXFI,KNUMER
C
        ELSEIF (KCODE.EQ.-7) THEN
          WRITE (UNIT=CLMESS,FMT='(''STATUS FORTRAN '''''',A,
     S           '''''' INCONNU, UNITE'',I3)') CDACTI(1:ILACTI),KNUMER
C
        ELSEIF (KCODE.EQ.-8) THEN
          WRITE (UNIT=CLMESS,
     S           FMT='(''L''''UNITE'',I3,'' DE STATUS ''''''
     S,A,'''''' DOIT AVOIR UN NOM EXPLICITE'')') KNUMER,CDACTI(1:ILACTI)
C
        ELSEIF (KCODE.EQ.-9) THEN
C
        IF (CDACTI.EQ.'OLD') THEN
          WRITE (UNIT=CLMESS,FMT=
     S'(''STATUS ''''OLD'''' MAIS LE FICHIER N''''EXISTE PAS, UNITE'',
     S      I3)') KNUMER
        ELSE
          ILBLAN=INDEX (CDACTI(1:ILACTI),' ')
          IF (ILBLAN.GT.1) ILACTI=ILBLAN-1
          WRITE (UNIT=CLMESS,FMT=
     S'(''STATUS '''''',A,'''''' MAIS LE FICHIER EXISTE DEJA, UNITE'',
     S  I3)') CDACTI(1:ILACTI),KNUMER
        ENDIF
C
        ELSEIF (KCODE.EQ.-10) THEN
          WRITE (UNIT=CLMESS,FMT='(''INCOMPATIBILITE'',
     S           '' FICHIER / LOGICIEL, UNITE'',I3)') KNUMER
C
        ELSEIF (KCODE.EQ.-11) THEN
          WRITE (UNIT=CLMESS,
     S           FMT='(''UNITE'',I3,'' NON FERMEE APRES '',
     S           ''LA DERNIERE MODIFICATION'')') KNUMER
C
        ELSEIF (KCODE.EQ.-12) THEN
          WRITE (UNIT=CLMESS,FMT='(''UNITE'',I3,
     S  '' DE STATUS ''''OLD'''' - ERREUR LECTURE PREMIER ARTICLE'')')
     S       KNUMER
C
        ELSEIF (KCODE.EQ.-13) THEN
          INLNOM=1
          INUMER=LFI%JPNIL
C
          DO 131 J=1,LFI%NBFIOU
          IJ=LFI%NUMIND(J)
C
          IF (CDACTI.EQ.LFI%CNOMFI(IJ)) THEN
            INUMER=LFI%NUMERO(IJ)
            INLNOM=MIN0 (LFI%NLNOMF(IJ),LEN (CLMESS)-3)
            GOTO 132
          ENDIF
C
  131     CONTINUE
C
  132     CONTINUE
          CLMESS=' '''//CDACTI(1:INLNOM)//''''
          WRITE (UNIT=LFI%NULOUT,FMT='(A)') CLMESS
          WRITE (UNIT=CLMESS,FMT='(''UNITE'',I3,'' - FICHIER '',
     S           ''DEJA OUVERT POUR L''''UNITE'',I3)') KNUMER,INUMER
C
        ELSEIF (KCODE.EQ.-14) THEN
C
          IF (CDNSPR.EQ.'LFIECR'.OR.CDNSPR.EQ.'LFILEC'.OR.
     S        CDNSPR.EQ.'LFILAS'.OR.CDNSPR.EQ.'LFILAP') THEN
            WRITE (UNIT=CLMESS,FMT=
     S   '(''LONGUEUR D''''ARTICLE INCORRECTE, UNITE'',I3)') KNUMER
          ELSEIF (KNUMER.EQ.LFI%JPNIL) THEN
            CLMESS='RANG DANS LA TABLE *LFI%NUMERO* INCORRECT'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S   '(''ARGUMENT DE TYPE ENTIER INCORRECT, UNITE'',I3)') KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-15) THEN
          WRITE (UNIT=CLMESS,
     S           FMT='(''NOM D''''ARTICLE INCORRECT OU '',
     S           ''TROP LONG, UNITE'',I3)') KNUMER
C
        ELSEIF (KCODE.EQ.-16) THEN
          WRITE (UNIT=CLMESS,
     S           FMT='(''INCOHERENCE (TABLES, FICHIER, '',
     S           ''APPELS S/P INT, LOGICIEL), UNITE'',I3)') KNUMER
C
        ELSEIF (KCODE.EQ.-17) THEN
C
          IF (IJL.NE.0) THEN
            INBALO=LFI%MDES1D(IXM(LFI%JPNALO,IJL))
          ELSE
            INBALO=LFI%JPNIL
          ENDIF
C
          WRITE (UNIT=CLMESS,
     S           FMT='(I6,'' ARTICLES, INDEX PLEIN, UNITE'',
     S           I3)') INBALO,KNUMER
C
        ELSEIF (KCODE.EQ.-18) THEN
          WRITE (UNIT=CLMESS,
     S           FMT='(''ARTICLE DE NOM BLANC ILLICITE'',
     S           '', UNITE'',I3)') KNUMER
C
        ELSEIF (KCODE.EQ.-19) THEN
          WRITE (UNIT=CLMESS,FMT='(''UNITE'',I3,
     S           '' ''''SCRATCH'''', NE PEUT ETRE CONSERVEE'')') KNUMER
C
        ELSEIF (KCODE.EQ.-20) THEN
          WRITE (UNIT=CLMESS,FMT='(''ARTICLE "'',A,
     S           ''" NON TROUVE, UNITE'',I3)') CDACTI(1:ILACT2),KNUMER
C
        ELSEIF (KCODE.EQ.-21) THEN
          WRITE (UNIT=CLMESS,FMT='(''ARTICLE "'',A,
     S    ''" + *LONG* QUE DEMANDE, UNITE'',I3)')
     S      CDACTI(1:ILACT2),KNUMER
C
        ELSEIF (KCODE.EQ.-22) THEN
          WRITE (UNIT=CLMESS,FMT='(''ARTICLE "'',A,
     S    ''" + *COURT* QUE DEMANDE, UNITE'',I3)')
     S      CDACTI(1:ILACT2),KNUMER
C
        ELSEIF (KCODE.EQ.-23) THEN
          WRITE (UNIT=CLMESS,
     S           FMT='(''PAS OU PLUS D''''ARTICLE SUIVANT'',
     S    '' A LIRE, UNITE'',I3)') KNUMER
C
        ELSEIF (KCODE.EQ.-24) THEN
          WRITE (UNIT=CLMESS,FMT='(''VARIABLE CAR.TROP COURTE '',
     S    ''POUR "'',A,''", UNITE'',I3)')
     S      CDACTI(1:ILACT2),KNUMER
C
        ELSEIF (KCODE.EQ.-25) THEN
          WRITE (UNIT=CLMESS,
     S           FMT='(''NOUVEAU NOM D''''ARTICLE: "'',A,
     S    ''" DEJA UTILISE, UNITE'',I3)')
     S      CDACTI(1:ILACT2),KNUMER
C
        ELSEIF (KCODE.EQ.-26) THEN
          WRITE (UNIT=CLMESS,FMT='(''PAS OU PLUS D''''ARTICLE '',
     S    '' PRECEDENT A LIRE, UNITE'',I3)') KNUMER
C
        ELSEIF (KCODE.EQ.-27) THEN
          WRITE (UNIT=CLMESS,
     S           FMT='(''ESPACE CONTIGU INSUFFISANT DANS '',
     S    '' LES TABLES, UNITE'',I3)') KNUMER
C
        ELSEIF (KCODE.EQ.-28) THEN
C
          IF (KNUMER.EQ.LFI%JPNIL) THEN
            WRITE (UNIT=CLMESS,FMT='(''FACTEUR MULTIPLICATIF PAR '',
     S      ''DEFAUT SUPERIEUR AU MAXIMUM('',I3,'')'')') LFI%JPFACX
          ELSE
            WRITE (UNIT=CLMESS,FMT='(''FACTEUR MULTIPLICATIF '',
     S      ''DEMANDE SUPERIEUR AU MAXIMUM ('',I3,''), UNITE'',I3)')
     S        LFI%JPFACX,KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-29) THEN
          WRITE (UNIT=CLMESS,FMT='(I3,'' ENTREES,'',
     S    '' PAS DE PLACE POUR FACTEUR MULTIPLIC, UNITE'',I3)')
     S    LFI%JPXUFM,KNUMER
C
        ELSEIF (KCODE.EQ.-30) THEN
          WRITE (UNIT=CLMESS,
     S           FMT='(''LFI%NUMERO D''''UNITE LOGIQUE FORTRAN''
     S           ,I8,'' ILLICITE'')') KNUMER
C
        ELSEIF (KCODE.EQ.-31) THEN
          WRITE (UNIT=CLMESS,FMT='(''LFI%NUMERO UNITE LOGIQ'',I3,
     S       '' SANS FACTEUR MULTIPLICATIF PREDEFINI'')') KNUMER
C
C                  Pour les codes d'erreur non prevus...
C
        ELSEIF (KNUMER.EQ.LFI%JPNIL) THEN
          WRITE (UNIT=CLMESS,
     S           FMT='(''ERREUR GLOBALE *INCONNUE* LFI%NUMERO'',
     S                             I6)') KCODE
        ELSE
          WRITE (UNIT=CLMESS,
     S           FMT='(''ERREUR *INCONNUE* LFI%NUMERO'',I6,
     S           '' SUR UNITE LOGIQUE'',I3)') KCODE,KNUMER
        ENDIF
C
        ILMESA=LEN (CLMESA)
        ILMESU=ILMESA-1-2*LEN (CLJOLI)-ILNSPR-4
        CLMESA=CLJOLI//' '//CDNSPR(1:ILNSPR)//' - '
     S         //CLMESS(1:ILMESU)//CLJOLI
        WRITE (UNIT=LFI%NULOUT,FMT='(A)') CLMESA
C
C           Si l'unite logique correspond a une unite logique LFI
C     deja ouverte, on en imprime le nom.
C
        IF (IJL.NE.0) THEN
C
          IF (LFI%NLNOMF(IJL).LE.LFI%JPLFTX) THEN
            WRITE (UNIT=LFI%NULOUT,FMT='(A,/)') CLJOLI
     S             //' NOM - APPARENT MAIS'
     S             //' COMPLET - DE L''UNITE LOGIQUE LFI CONCERNEE:'
          ELSE
            WRITE (UNIT=CLMESS,FMT='(A,
     S             '' NOM - APPARENT, ET TRONQUE DE'',I4,
     S       '' CARACTERES - DE L''''UNITE LOGIQUE LFI CONCERNEE:'')')
     S      CLJOLI,LFI%NLNOMF(IJL)-LFI%JPLFTX
            WRITE (UNIT=LFI%NULOUT,FMT='(A,/)') CLMESS
          ENDIF
C
          INLIGN=(LFI%NLNOMF(IJL)-1)/LFI%JPLFIX
          IDECAL=0
C
          DO 801 J=1,INLIGN
          WRITE (UNIT=LFI%NULOUT,FMT='(A)')
     S           LFI%CNOMFI(IJL)(IDECAL+1:IDECAL+LFI%JPLFIX)//'...'
          IDECAL=IDECAL+LFI%JPLFIX
  801     CONTINUE
C
          IF (LFI%NLNOMF(IJL).LE.LFI%JPLFTX) THEN
            WRITE (UNIT=LFI%NULOUT,FMT='(A,/)') 
     S            LFI%CNOMFI(IJL)(IDECAL+1:LFI%NLNOMF(IJL))
          ELSE
            WRITE (UNIT=LFI%NULOUT,FMT='(A,/)') 
     S             LFI%CNOMFI(IJL)(IDECAL+1:LFI%JPLFTX)//'...'
          ENDIF
C
          IF (LFI%CNOMSY(IJL).NE.LFI%CNOMFI(IJL)) THEN
            WRITE (UNIT=LFI%NULOUT,FMT='(A,/)') CLJOLI//
     S ' NOM *SYSTEME* (APPARENT) DE L''UNITE LOGIQUE LFI CONCERNEE:'
            INLIGN=(LFI%NLNOMS(IJL)-1)/LFI%JPLFIX
            IDECAL=0
C
            DO 802 J=1,INLIGN
            WRITE (UNIT=LFI%NULOUT,FMT='(A)')
     S             LFI%CNOMSY(IJL)(IDECAL+1:IDECAL+LFI%JPLFIX)//'...'
            IDECAL=IDECAL+LFI%JPLFIX
  802       CONTINUE
C
            WRITE (UNIT=LFI%NULOUT,FMT='(A,/)') 
     S            LFI%CNOMSY(IJL)(IDECAL+1:LFI%NLNOMS(IJL))
          ENDIF
C
        ENDIF
C
        WRITE (UNIT=LFI%NULOUT,FMT='(A)') CLMESA
        IF (LDFATA.AND.KCODE.NE.0) THEN
C
C            Saborde le programme.
C
          CALL SDL_SRLABORT
        ENDIF
C
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('LFIEFR_MT',1,ZHOOK_HANDLE)
      END

