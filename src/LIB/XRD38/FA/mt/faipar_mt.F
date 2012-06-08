C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FAIPAR_MT (FA,  KNUMER, KNIMES, KCODE, LDFATA, CDMESS, 
     S                      CDNSPR, CDACTI, LDRLFI )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C**** 
C        Ce sous-programme est charge de FAIre PARt des Messages
C     emis par le logiciel de Fichiers ARPEGE, en faisant si besoin est
C     l'"ABORT" du programme .
C**
C        Arguments : KNUMER ==> Numero de l'unite logique eventuelle;
C                               ( FA%JPNIIL ==> pas d'unite logique )
C                    KNIMES ==> Niveau (0,1,2) du message;
C                    KCODE  ==> Code correspondant a l'action;
C                    LDFATA ==> Vrai si on doit avorter le programme;
C                    CDMESS ==> si KNIMES#0, Message a emettre;
C                    CDNSPR ==> Nom du Sous-Programme appelant;
C                    CDACTI ==> Nom de l'action d'Entree/Sortie FORTRAN
C                               (si KCODE >0), sinon fourre-tout !) .
C                    LDRLFI ==> Vrai si le code-reponse different de 0
C                               a ete detecte par le logiciel LFI.
C*
C     VALEURS POSSIBLES DES CODES-REPONSES DU LOGICIEL FICHIERS ARPEGE
C
C      0 ==> Aucune erreur detectee, message informatif.
C     >0 ==> Il s'agit du code-reponse d'une erreur detectee lors d'un
C            OPEN, READ, WRITE, CLOSE, ou INQUIRE FORTRAN; voir manuels.
C -(1a50)==> Code-reponse renvoye par le logiciel LFI.
C    -51 ==> Unite logique non ouverte pour le logiciel,
C             ou bien cadre non defini .
C    -52 ==> Valeur d'un "NIVEAU" hors plage [0-2] .
C    -53 ==> *** code-reponse inutilise pour le moment ***.
C    -54 ==> Changement explicite de mode Multi-Tasking avec au moins 1
C            unite logique ouverte-risque de problemes (s/p "FARINE") .
C    -55 ==> Unite logique deja ouverte ou cadre non redefinissable.
C    -56 ==> Trop d'unites logiques deja ouvertes ou de cadres definis.
C    -57 ==> Nouveau fichier, et cadre reference non defini au prealable
C    -58 ==> Fichier preexistant, et cadre reference non compatible .
C    -59 ==> Redefinition de cadre avec modification de parametre(s),
C            impossible car il y a au moins un fichier qui s'y rattache.
C    -60 ==> Un des articles definissant le cadre est manquant.
C    -61 ==>  "  "      "        "       "  " a une longueur inattendue.
C    -62 ==> L'article DATE est manquant.
C    -63 ==> L'article DATE a une longueur inattendue.
C    -64 ==> Argument d'appel de type ENTIER incorrect .
C    -65 ==> Argument d'appel de type CARACTERE incorrect .
C    -66 ==> INCOHERENCE tables, fichier, appels s/p internes, logiciel;
C            CE CODE EST LE SEUL QUE *FA%NRFAGA* VALANT 2 NE PEUT MASQUER,
C            ET PROVOQUE DONC UNE ERREUR FATALE DANS TOUS LES CAS.
C    -67 ==> Le cadre reference a au moins un fichier associe,
C            et ne peut donc pas etre supprime.
C    -68 ==> Nom d'article compose uniquement de BLANCS, non accepte
C            ( a cause du fonctionnement interne du logiciel LFI ),
C            ou cadre de nom blanc non accepte (les cadres non definis
C             sont reperes ainsi dans les tables du logiciel)
C    -69 ==> La variable CARACTERE donnee en argument est TROP COURTE
C            pour recevoir le NOM du cadre ou d'identificateur, meme en
C            supprimant les eventuels caracteres blancs en fin de nom .
C    -70 ==> Troncature hors plage [1-FA%NXTRON].
C    -71 ==> Nombre de latitudes de pole a pole hors plage [1-FA%NXLATI].
C    -72 ==> Nombre de niveaux verticaux hors plage [1-FA%NXNIVV].
C    -73 ==> Coefficient de dilatation inferieur a 1.
C    -74 ==> Un des nombres de points par parallele est <= 0 ou depasse
C            le nombre maximum annonce dans la description du cadre.
C    -75 ==> Un des nombres de points par parallele est inferieur
C            a celui qui le precede.
C    -76 ==> Un des nombres d'onde zonal maxi par parallele depasse
C            la troncature ou n'est pas positif.
C    -77 ==> Un des nombres d'onde zonal maxi par parallele est
C            inferieur a celui qui le precede.
C    -78 ==> Sur une latitude au moins, le nombre de points est insuffi-
C            sans vis-a-vis du nombre d'onde zonal maximum.
C    -79 ==> Inconsistance: le nombre de latitudes (de pole a pole)
C            est insuffisant vis-a-vis de la troncature.
C    -80 ==> Une valeur d'une des fonctions de la coordonnee hybride
C            est en-dehors de l'intervalle [0,1] .
C    -81 ==> Pour une couche, valeurs inadaptees des fonctions de la
C            coordonnee hybride.
C    -82 ==> L'article DATE a un contenu incorrect.
C    -83 ==> Le nombre maximum de points par parallele est en-dehors de
C            l'intervalle [1,FA%NXLONG]
C    -84 ==> Le nombre maximum de points par parallele est insuffisant
C            vis-a-vis de la troncature.
C    -85 ==> Fichier en mode creation, donc vide d'articles vis-a-vis
C            du logiciel (a la description du Cadre pres):
C            il faut le munir d'une date...
C    -86 ==> Prefixe et/ou suffixe "blanc" interdit.
C    -87 ==> Prefixe de champ trop long.
C    -88 ==> Troncature effective inferieure ou egale a
C            la sous-troncature "non compactee".
C    -89 ==> Article de type champ demande inexistant dans le fichier.
C    -90 ==> L'article de type champ existe, mais est trop long.
C    -91 ==> En-tete d'article champ incorrect.
C    -92 ==> L'article demande est en points-de-grille au lieu de
C            coefficients spectraux, ou vice-versa.
C    -93 ==> L'article demande est trop COURT sur le fichier.
C    -94 ==> L'article demande est trop LONG sur le fichier ( si cette
C            erreur n'est pas fatale, on traite le debut de l'article ).
C    -95 ==> Incoherence entre entete d'article champ et zone GRIB.
C    -96 ==> Niveau de codage GRIB trop grand.
C    -97 ==> Nombre de bits de codage trop grand.
C    -98 ==> Puissance de laplacien trop grande.
C    -99 ==> Sous-troncature superieure ou egale a la troncature
C   -100 ==> Une des lignes trigonometriques definissant le pole de
C            projection n'est pas dans [-1,1] .
C   -101 ==> Cosinus et Sinus de la longitude du pole de projection
C            incoherents.
C   -102 ==> Le Sinus d'une latitude de la grille n'est pas dans [-1,1].
C   -103 ==> Le Sinus d'une latitude est superieur ou egal a celui de
C            la latitude precedente.
C   -104 ==> Troncature maxi. incompatible avec cadre(s) deja defini(s).
C   -105 ==> Nombre maxi. de niveaux verticaux "  "       "      "     .
C   -106 ==> Nombre maxi. de latitudes pole a pole"       "      "     .
C   -107 ==> Nombre maxi. de longitudes "      "  "       "      "     .
C   -108 ==> Pression de reference aberrante.
C   -109 ==> Type de transformation horizontale hors plage [1-FA%NTYPTX]
C   -110 ==> Pas d'article "identificateur" (apres ceux du cadre)
C   -111 ==> Article de nom reserve, a ne pas ecrire ou lire ainsi
C   -112 ==> Caracteristiques des fichiers vraiment incompatibles,
C            recopie d'article "champ" impossible.
C   -113 ==> Sous-troncature implicite superieure ou egale a la
C            limite usager de troncature.
C   -114 ==> Ratio des troncatures (version LAM) superieur a 3.
C            Garde-fou, modifiable dans FARCIS, FACSIM et FAPULA.
C   -115 ==> (LAM) Le nombre maximum de points par parallele est insuffisant
C            vis-a-vis de la troncature.
C   -116 ==> (LAM) Le nombre maximum de points par meridien est insuffisant
C            vis-a-vis de la troncature.
C   -117 ==> (LAM) L'indicateur NDOMFP a une valeur hors de [-1,1]
C   -118 ==> (LAM) L'indice de depart (C+I) en longitude est hors plage [1-KNXLON].
C   -119 ==> (LAM) L'indice de fin (C+I) en longitude est errone.
C   -120 ==> (LAM) L'indice de depart (C+I) en latitude est hors plage [1-KNLATI].
C   -121 ==> (LAM) L'indice de fin (C+I) en latitude est errone.
C   -122 ==> (LAM) La largeur de la zone de relaxation est trop large en longitude.
C   -123 ==> (LAM) La largeur de la zone de relaxation est trop large en latitude.
C   -124 ==> Nombre de bits de codage nul alors qu'un codage est demande
C   -125 ==> Argument ayant une valeur incorrecte
C   -126 ==> Detection d'une incoherence dans un controle interne
C   -127 ==> (LAM) Description incoherente de l'ellipse contenant les troncatures
C   -128 ==> Message GRIB incorrect dans l'appel a DECF10
C <-1000 ==> Erreur rapportee par GRIBEX
C

      USE SDL_MOD   , ONLY : SDL_SRLABORT

#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KNUMER, KNIMES, KCODE
C
      INTEGER ILNSPR, ILACTI, ILACT2, ILMESU, ILMESA, INUME1, INUME2
C
      LOGICAL LDFATA, LDRLFI
C
      CHARACTER  CDNSPR*(*), CLJOLI*6, CDMESS*(*), CLMESA*80, CDACTI*(*)
#include "facom_mt.h"
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FAIPAR_MT',0,ZHOOK_HANDLE)
      ILNSPR=MIN (LEN (CDNSPR),6)
      ILACTI=MIN (LEN (CDACTI),8)
      ILACT2=MIN (LEN (CDACTI),FA%NCPCAD)
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
        ILMESU=MIN (LEN (CLMESS)-LEN (CLJOLI)-ILNSPR-4,LEN (CDMESS))
        CLMESS=CLJOLI//' '//CDNSPR(1:ILNSPR)//' - '//CDMESS(1:ILMESU)
        WRITE (UNIT=FA%NULOUT,FMT='(A)') CLMESS
      ENDIF
C
      IF (KNIMES.EQ.0.OR.LDFATA) THEN
C
C        CONSTITUTION D'UN MESSAGE D'ERREUR AD HOC, EN FONCTION DE KCODE
C
        IF (LDRLFI) THEN
          WRITE (UNIT=CLMESS,FMT='(''CODE-REPONSE DE *'',A6,''* ='',
     S           I4,'', UNITE ='',I3)') CDACTI(1:ILACTI),KCODE,KNUMER
C
        ELSEIF (KCODE.GT.0) THEN
          WRITE (UNIT=CLMESS,
     S           FMT='(''ERREUR "'',A,''" FORTRAN, CODE ='',
     S           I4,'', UNITE ='',I3)') CDACTI(1:ILACTI),KCODE,KNUMER
C
        ELSEIF (KCODE.EQ.-51) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='CADRE "'//CDACTI(1:ILACT2)//'" NON DEFINI'
          ELSE
            WRITE (UNIT=CLMESS,FMT='(''UNITE LOGIQUE'',I3,
     S             '' NON OUVERTE'')') KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-52) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='PARAMETRE DE NIVEAU "KNIVAU" HORS PLAGE [0-2]'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S   '(''NIVEAU DE MESSAGERIE HORS PLAGE [0-2], UNITE'',I3)') KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-53) THEN
          CLMESS='CADRE "'//CDACTI(1:ILACT2)//
     S           '", PARAMETRE(S) INCORRECT(S)'
C
        ELSEIF (KCODE.EQ.-54) THEN
          CLMESS=
     S      'CHANGEMENT EXPLICITE MODE MULTI AVEC UNITES OUVERTES'
C
        ELSEIF (KCODE.EQ.-55) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='CADRE "'//CDACTI(1:ILACT2)//
     S                '" NON REDEFINISSABLE'
          ELSE
            WRITE (UNIT=CLMESS,FMT='(''UNITE LOGIQUE'',I3,
     S             '' DEJA OUVERTE'')') KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-56) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            WRITE (UNIT=CLMESS,FMT=
     S         '(I2,'' CADRES DEFINIS, PAS DE PLACE POUR "'',A,''"'')')
     S      FA%JPNXCA,CDACTI(1:ILACT2)
          ELSE
            WRITE (UNIT=CLMESS,FMT='(I2,
     S           '' UNITES LOGIQUES OUVERTES, PAS DE PLACE POUR'',I3)')
     S      FA%JPNXFA,KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-57) THEN
          WRITE (UNIT=CLMESS,
     S           FMT='(''UNITE'',I3,'' A CREER, CADRE "'',A,
     S           ''" NON DEFINI'')') KNUMER,CDACTI(1:ILACT2)
C
        ELSEIF (KCODE.EQ.-58) THEN
          WRITE (UNIT=CLMESS,
     S           FMT='(''UNITE'',I3,'' /CADRE PREDEFINI "'',
     S           A,''" INCOMPATIBLES'')') KNUMER,CDACTI(1:ILACT2)
C
        ELSEIF (KCODE.EQ.-59) THEN
          CLMESS='CADRE "'//CDACTI(1:ILACT2)//
     S           '"  NON MODIFIABLE CAR UTILISE'
C
        ELSEIF (KCODE.EQ.-60) THEN
          WRITE (UNIT=CLMESS,FMT=
     S       '(''UN ARTICLE DU CADRE EST MANQUANT, UNITE'',I3)') KNUMER
C
        ELSEIF (KCODE.EQ.-61) THEN
C***** FAZZZZ - un article du cadre a une longueur inattendue, uniteiii *****
          WRITE (UNIT=CLMESS,FMT=
     S  '(''UN ARTICLE DU CADRE A UNE LONGUEUR INATTENDUE, UNITE'',I3)')
     S    KNUMER
C
        ELSEIF (KCODE.EQ.-62) THEN
          WRITE (UNIT=CLMESS,FMT=
     S       '(''PAS D''''ARTICLE DATE SUR L''''UNITE'',I3)') KNUMER
C
        ELSEIF (KCODE.EQ.-63) THEN
          WRITE (UNIT=CLMESS,FMT=
     S'(''L''''ARTICLE DATE A UNE LONGUEUR INATTENDUE SUR L'''' UNITE'',
     S    I3)')  KNUMER
C
        ELSEIF (KCODE.EQ.-64) THEN
C
          IF (CDNSPR.EQ.'FAIENC') THEN
            WRITE (UNIT=CLMESS,FMT=
     S   '(''NIVEAU VERTICAL HORS LIMITES, UNITE'',I3)') KNUMER
          ELSEIF (CDNSPR.EQ.'FAISAN'.OR.CDNSPR.EQ.'FALAIS') THEN
            WRITE (UNIT=CLMESS,FMT=
     S   '(''LONGUEUR D''''ARTICLE <=0, UNITE'',I3)') KNUMER
          ELSEIF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='RANG DANS UNE TABLE GLOBALE INCORRECT'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S   '(''ARGUMENT ENTIER INCORRECT, UNITE'',I3)') KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-65) THEN
C
          IF (CDNSPR.EQ.'FAIENC'.OR.CDNSPR.EQ.'FACILE'.OR.
     S        CDNSPR.EQ.'FAISAN'.OR.CDNSPR.EQ.'FALAIS') THEN
            WRITE (UNIT=CLMESS,
     S             FMT='(''NOM D''''ARTICLE INCORRECT OU '',
     S           ''TROP LONG, UNITE'',I3)') KNUMER
          ELSEIF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='NOM DE CADRE INCORRECT OU TROP LONG: "'
     S             //CDACTI(1:ILACT2)//'"'
          ELSE
            WRITE (UNIT=CLMESS,
     S             FMT='(''CHAINE DE CARACTERES INCORRECT'',
     S           ''E OU TROP LONGUE, UNITE'',I3)') KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-66) THEN
          WRITE (UNIT=CLMESS,
     S           FMT='(''INCOHERENCE (TABLES, FICHIER, '',
     S           ''APPELS S/P INT. OU FICHIER), UNITE'',I3)') KNUMER
C
        ELSEIF (KCODE.EQ.-67) THEN
          CLMESS='FICHIER ASSOCIE, CADRE"'//CDACTI(1:ILACT2)
     S           //'" DEVANT SUBSISTER'
C
        ELSEIF (KCODE.EQ.-68) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='CADRE DE NOM BLANC NON ACCEPTE'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S        '(''ARTICLE DE NOM BLANC NON ACCEPTE, UNITE'',I3)') KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-69) THEN
          WRITE (UNIT=CLMESS,FMT='(''VARIABLE CARACT. TROP COURTE '',
     S    ''POUR "'',A,''", UNITE'',I3)')
     S      CDACTI(1:ILACT2),KNUMER
C
        ELSEIF (KCODE.EQ.-70) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            WRITE (UNIT=CLMESS,FMT='(''TRONCATURE HORS PLAGE [1-'',
     S                               I3,''], CADRE "'',A,''"'')')
     S      FA%NXTRON,CDACTI(1:ILACT2)
          ELSE
            WRITE (UNIT=CLMESS,FMT='(''TRONCATURE HORS PLAGE [1-'',
     S                      I3,''], UNITE'',I3)') FA%NXTRON,KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-71) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            WRITE (UNIT=CLMESS,FMT='(''NOMBRE LATITUDES HORS [1-'',
     S             I4,''], CADRE "'',A,''"'')') 
     S             FA%NXLATI,CDACTI(1:ILACT2)
          ELSE
            WRITE (UNIT=CLMESS,FMT='(''NOMBRE LATITUDES HORS [1-'',
     S             I4,''], UNITE'',I3)') FA%NXLATI,KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-72) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            WRITE (UNIT=CLMESS,FMT='(''NB. NIV. VERTICAUX HORS [1-'',
     S             I3,''], CADRE "'',A,''"'')') 
     S             FA%NXNIVV,CDACTI(1:ILACT2)
          ELSE
            WRITE (UNIT=CLMESS,FMT='(''NB. NIV. VERTICAUX HORS [1-'',
     S             I3,''], UNITE'',I3)') FA%NXNIVV,KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-73) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='COEFF. DILATATION INFERIEUR A 1, CADRE "'
     S             //CDACTI(1:ILACT2)//'"'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S         '(''COEFF. DILATATION INFERIEUR A 1, UNITE'',I3)') KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-74) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='1 NB. POINTS/PARALL. HORS PLAGE, CADRE "'
     S             //CDACTI(1:ILACT2)//'"'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S         '(''1 NB. POINTS/PARALL. HORS PLAGE, UNITE'',I3)') KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-75) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='LISTE POINTS/PARALL. INCORRECTE, CADRE "'
     S             //CDACTI(1:ILACT2)//'"'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S         '(''LISTE POINTS/PARALL. INCORRECTE, UNITE'',I3)') KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-76) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='1 N.O.Z. MAX/PARALL. HORS PLAGE, CADRE "'
     S             //CDACTI(1:ILACT2)//'"'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S         '(''1 N.O.Z. MAX/PARALL. HORS PLAGE, UNITE'',I3)') KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-77) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='LISTE NOZMAX/PARALL. INCORRECTE, CADRE "'
     S             //CDACTI(1:ILACT2)//'"'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S         '(''LISTE NOZMAX/PARALL. INCORRECTE, UNITE'',I3)') KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-78) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='INCONSISTANCE NOZMAX/NB. POINTS, CADRE "'
     S             //CDACTI(1:ILACT2)//'"'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S         '(''INCONSISTANCE NOZMAX/NB. POINTS, UNITE'',I3)') KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-79) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='INCONSISTANCE NB.LAT/TRONCATURE, CADRE "'
     S             //CDACTI(1:ILACT2)//'"'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S         '(''INCONSISTANCE NB.LAT/TRONCATURE, UNITE'',I3)') KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-80) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='FONC. COORD. HYBRIDE HORS [0,1], CADRE "'
     S             //CDACTI(1:ILACT2)//'"'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S         '(''FONC. COORD. HYBRIDE HORS [0,1], UNITE'',I3)') KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-81) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
C***** FAZZZZ - fonc. coord. hybride inadaptees, CADRE "0123456789ABCDEF" *****
            CLMESS='FONC. COORD. HYBRIDE INADAPTEES, CADRE "'
     S             //CDACTI(1:ILACT2)//'"'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S         '(''FONC. COORD. HYBRIDE INADAPTEES, UNITE'',I3)') KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-82) THEN
C
          IF (CDNSPR.EQ.'FANDAR') THEN
C***** FAZZZZ - La DATE proposee a un contenu incorrect, UNITEiii *****
            WRITE (UNIT=CLMESS,FMT=
     S       '(''LA DATE PROPOSEE A UN CONTENU INCORRECT, UNITE'',I3)')
     S      KNUMER
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S       '(''L''''ARTICLE DATE A UN CONTENU INCORRECT, UNITE'',I3)')
     S      KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-83) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
C***** FAZZZZ - Nb. Max. longitu. hors [1-iiii], CADRE "0123456789abcdef" *****
            WRITE (UNIT=CLMESS,FMT='(''NB. MAX. LONGITU. HORS [1-'',
     S             I4,''], CADRE "'',A,''"'')') 
     S             FA%NXLONG,CDACTI(1:ILACT2)
          ELSE
            WRITE (UNIT=CLMESS,FMT='(''NB. MAX. LONGITU. HORS [1-'',
     S             I4,''], UNITE'',I3)') FA%NXLONG,KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-84) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='INCONSISTANCE KNXLON/TRONCATURE, CADRE "'
     S             //CDACTI(1:ILACT2)//'"'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S         '(''INCONSISTANCE KNXLON/TRONCATURE, UNITE'',I3)') KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-85) THEN
          WRITE (UNIT=CLMESS,FMT='(''UNITE'',I3,
     S           '' PAS ENCORE MUNIE D''''UNE DATE'')') KNUMER
C
        ELSEIF (KCODE.EQ.-86) THEN
          WRITE (UNIT=CLMESS,FMT=
     S         '(''PREFIXE/SUFFIXE BLANC INTERDIT, UNITE'',I3)') KNUMER
C
        ELSEIF (KCODE.EQ.-87) THEN
          WRITE (UNIT=CLMESS,FMT=
     S         '(''PREFIXE DE CHAMP TROP LONG, UNITE'',I3)') KNUMER
C
        ELSEIF (KCODE.EQ.-88) THEN
          WRITE (UNIT=CLMESS,FMT=
     S  '(''TRONCATURE <= SOUS-TRONCATURE "NON COMPACTEE", UNITE'',I3)')
     S       KNUMER
C
        ELSEIF (KCODE.EQ.-89) THEN
          WRITE (UNIT=CLMESS,FMT='(''ARTICLE-CHAMP "'',A,
     S           ''" INEXISTANT, UNITE'',I3)') CDACTI(1:ILACT2),KNUMER
C
        ELSEIF (KCODE.EQ.-90) THEN
          WRITE (UNIT=CLMESS,FMT='(''ARTICLE-CHAMP "'',A,
     S           ''" TROP LONG, UNITE'',I3)') CDACTI(1:ILACT2),KNUMER
C
        ELSEIF (KCODE.EQ.-91) THEN
C***** FAZZZZ - Entete d'article "0123456789abcdef" incorrect, uniteiii *****
          WRITE (UNIT=CLMESS,FMT='(''ENTETE D''''ARTICLE "'',A,
     S           ''" INCORRECT, UNITE'',I3)') CDACTI(1:ILACT2),KNUMER
C
        ELSEIF (KCODE.EQ.-92) THEN
C***** FAZZZZ - Desaccord Csp./pdg., Article "0123456789abcdef", uniteiii *****
          WRITE (UNIT=CLMESS,
     S           FMT='(''DESACCORD CSP./PDG., ARTICLE "'',A,
     S           ''" , UNITE'',I3)') CDACTI(1:ILACT2),KNUMER
C
        ELSEIF (KCODE.EQ.-93) THEN
          WRITE (UNIT=CLMESS,FMT='(''ARTICLE "'',A,
     S           ''" TROP *COURT* SUR L''''UNITE'',I3)')
     S        CDACTI(1:ILACT2),KNUMER
C
        ELSEIF (KCODE.EQ.-94) THEN
          WRITE (UNIT=CLMESS,FMT='(''ARTICLE "'',A,
     S           ''" TROP *LONG* SUR L''''UNITE'',I3)')
     S        CDACTI(1:ILACT2),KNUMER
C
        ELSEIF (KCODE.EQ.-95) THEN
C***** FAZZZZ - inconsis entete/GRIB article "0123456789abcdef", uniteiii *****
          WRITE (UNIT=CLMESS,
     S           FMT='(''INCONSIS ENTETE/GRIB ARTICLE "'',
     S           ''", UNITE'',I3)')
     S        CDACTI(1:ILACT2),KNUMER
C
        ELSEIF (KCODE.EQ.-96) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='NIVEAU IMPLICITE DE CODAGE GRIB TROP GRAND'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S       '(''NIVEAU DE CODAGE GRIB TROP GRAND, UNITE'',I3)') KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-97) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            WRITE (UNIT=CLMESS,FMT=
     S  '(''NOMBRE(S) IMPLICITE(S) DE BITS PAR VALEUR > MAXI ('',I2,
     S      '')'')')   FA%NBIMAX
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S  '(''NOMBRE(S) DE BITS PAR VALEUR > MAXI ('',I2,''), UNITE'',
     S      I3)')   FA%NBIMAX,KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-98) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS=
     S       'PUISSANCE IMPLICITE DE LAPLACIEN TROP GRANDE'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S     '(''PUISSANCE DE LAPLACIEN TROP GRANDE, UNITE'',I3)')
     S         KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-99) THEN
C
          IF (CDNSPR.EQ.'FARFLU') THEN
            CLMESS=
     S       'TRONCATURE MAX INFERIEURE OU = SOUS-TRONCATURE IMPLICITE'
          ELSEIF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS=
     S       'SS-TRONCATURE IMPLICITE SUPERIEURE OU = A UNE TRONCATURE'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S     '(''SOUS-TRONCATURE SUPERIEURE OU = A LA TRONCATURE, UNITE'',
     S       I3)')  KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-100) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='SIN/COS POLE PROJEC. HORS [-1,1], CADRE "'
     S             //CDACTI(1:ILACT2)//'"'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S     '(''SINUS/COSINUS DU POLE DE PROJECTION HORS [-1,1], UNITE'',
     S       I3)')  KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-101) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='INCOHERENCE COS/SIN LONGIT. POLE, CADRE "'
     S             //CDACTI(1:ILACT2)//'"'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S     '(''INCOHERENCE COS/SIN LONGITUDE POLE PROJECTION, UNITE'',
     S       I3)')  KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-102) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='SINUS D''UNE LATITUDE HORS [-1,1], CADRE "'
     S             //CDACTI(1:ILACT2)//'"'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S    '(''SINUS D''''UNE LATITUDE DE LA GRILLE HORS [-1,1], UNITE'',
     S       I3)')  KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-103) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
C***** FAZZZZ - sinus latitudes non decroissants, cadre '0123456789abcdef" *****
            CLMESS='SINUS LATITUDES NON DECROISSANTS, CADRE "'
     S             //CDACTI(1:ILACT2)//'"'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S    '(''SINUS DES LATITUDES NON DECROISSANTS, UNITE'',
     S       I3)')  KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-104) THEN
C
          CLMESS=
     S     'TRONCATURE MAXI INCOMPATIBLE AVEC CADRE(S) DEJA DEFINI(S)'
C
        ELSEIF (KCODE.EQ.-105) THEN
C
          CLMESS=
     S     'NB MAX NIVEAUX VERTICAUX/CADRE(S) DEFINI(S) INCOMPATIBLES'
C
        ELSEIF (KCODE.EQ.-106) THEN
C
          CLMESS=
     S     'NOMBRE MAXI DE LATITUDES/CADRE(S) DEFINI(S) INCOMPATIBLES'
C
        ELSEIF (KCODE.EQ.-107) THEN
C
          CLMESS=
     S     'NOMBRE MAX DE LONGITUDES/CADRE(S) DEFINI(S) INCOMPATIBLES'
C
        ELSEIF (KCODE.EQ.-108) THEN
C
C***** FAZZZZ - Pression de reference aberrante, cadre "0123456789abcdef" *****
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='PRESSION DE REFERENCE ABERRANTE, CADRE "'
     S             //CDACTI(1:ILACT2)//'"'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S             '(''PRESSION DE REFERENCE ABERRANTE, UNITE'',I3)')
     S          KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-109) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            WRITE (UNIT=CLMESS,
     S             FMT='(''TYPE TRANSFO HORIZ. HORS [1-'',
     S             I2,''], CADRE "'',A,''"'')') 
     S             FA%NTYPTX,CDACTI(1:ILACT2)
          ELSE
            WRITE (UNIT=CLMESS,
     S             FMT='(''TYPE TRANSFO HORIZ. HORS [1-'',
     S             I2,''], UNITE'',I3)') FA%NTYPTX,KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-110) THEN
          WRITE (UNIT=CLMESS,FMT=
     S    '(''PAS D''''ARTICLE DERRIERE CEUX DU CADRE, UNITE'',I3)')
     S           KNUMER
C
        ELSEIF (KCODE.EQ.-111) THEN
          WRITE (UNIT=CLMESS,FMT=
     S    '(''NOM RESERVE, NE PEUT ETRE UTILISE AINSI, UNITE'',I3)')
     S           KNUMER
C
        ELSEIF (KCODE.EQ.-112) THEN
          INUME1=KNUMER/1000
          INUME2=MOD (KNUMER,1000)
          WRITE (UNIT=CLMESS,
     S           FMT='(''COPIE DE '''''',A,'''''': UNITES'',
     S           I3,'' ET'',I3,'' INCOMPATIBLES'')') CDACTI(1:ILACT2),
     S      INUME1,INUME2
C
        ELSEIF (KCODE.EQ.-113) THEN
          CLMESS=
     S       'SS-TRONCATURE IMPLICITE SUPERIEURE OU = A TRONCATURE MAXI'
C
        ELSEIF (KCODE.EQ.-114) THEN
          CLMESS=
     S       'RATIO DES TRONCATURES HORIZ. (VERSION LAM) SUPERIEUR A 3:
     S        GARDE-FOU, MODIFIABLE DANS FARCIS+FACSIM'
C
        ELSEIF (KCODE.EQ.-115) THEN
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='KNXLON<(-2*KTYPTR+1), CADRE "'
     S             //CDACTI(1:ILACT2)//'"'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S         '(''KNXLON<(-2*KTYPTR+1), UNITE'',I3)') KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-116) THEN
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='KNXLON<(2*KTRONC+1), CADRE "'
     S             //CDACTI(1:ILACT2)//'"'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S         '(''KNXLON<(2*KTRONC+1), UNITE'',I3)') KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-117) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='INDICATEUR DE DOMAINE HORS [-1,1], CADRE "'
     S             //CDACTI(1:ILACT2)//'"'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S    '(''INDICATEUR DE DOMAINE HORS [-1,1], UNITE'',
     S       I3)')  KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-118) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS=
     S       'INDICE DE DEPART EN LONGITUDE HORS [1,KNXLON], CADRE "' 
     S             //CDACTI(1:ILACT2)//'"'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S    '(''INDICE DE DEPART EN LONGITUDE HORS [1,KNXLON], UNITE'',
     S       I3)')  KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-119) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='INDICE DE FIN EN LONGITUDE ERRONE, CADRE "'
     S             //CDACTI(1:ILACT2)//'"'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S    '(''INDICE DE FIN EN LONGITUDE ERRONE, UNITE'',
     S       I3)')  KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-120) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS=
     S       'INDICE DE DEPART EN LATITUDE HORS [1,KNLATI], CADRE "'
     S             //CDACTI(1:ILACT2)//'"'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S    '(''INDICE DE DEPART EN LATITUDE HORS [1,KNLATI], UNITE'',
     S       I3)')  KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-121) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='INDICE DE FIN EN LATITUDE ERRONE, CADRE "'
     S             //CDACTI(1:ILACT2)//'"'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S    '(''INDICE DE FIN EN LATITUDE ERRONE, UNITE'',
     S       I3)')  KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-122) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='ZONE DE RELAXATION TROP LARGE '//
     S                'EN LONGITUDE, CADRE "'
     S                //CDACTI(1:ILACT2)//'"'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S    '(''ZONE DE RELAXATION TROP LARGE EN LONGITUDE, UNITE'',
     S       I3)')  KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-123) THEN
C
          IF (KNUMER.EQ.FA%JPNIIL) THEN
            CLMESS='ZONE DE RELAXATION TROP LARGE '//
     S                'EN LATITUDE, CADRE "'
     S                //CDACTI(1:ILACT2)//'"'
          ELSE
            WRITE (UNIT=CLMESS,FMT=
     S    '(''ZONE DE RELAXATION TROP LARGE EN LATITUDE, UNITE'',
     S       I3)')  KNUMER
          ENDIF
C
        ELSEIF (KCODE.EQ.-124) THEN
          CLMESS=
     S       'NB DE BITS POUR CODER NUL'
C
        ELSEIF (KCODE.EQ.-125) THEN
          CLMESS=
     S       'ARGUMENT AYANT UNE VALEUR INCORRECTE'
C
        ELSEIF (KCODE.EQ.-126) THEN
          CLMESS=
     S       'INCOHERENCE DANS UN CONTROLE INTERNE'
C
        ELSEIF (KCODE.EQ.-127) THEN
          CLMESS=
     S       'INCOHERENCE DANS LE CONTROLE DE L''ELLIPSE'
C
        ELSEIF (KCODE.EQ.-128) THEN
          CLMESS=
     S       'MESSAGE GRIB INCORRECT POUR DECF10'
C
        ELSEIF (KCODE.LT.-1000) THEN
          CLMESS=
     S  'ERREUR RAPPORTEE PAR GRIBEX (cf manuel, avec kret=-krep-1000)'
C
        ELSEIF (KNUMER.EQ.FA%JPNIIL) THEN
          WRITE (UNIT=CLMESS,FMT='(''ERREUR GLOBALE NUMERO'',I6)') 
     S            KCODE
C
        ELSE
          WRITE (UNIT=CLMESS,FMT='(''ERREUR GRIB'',I6,
     S           '' SUR UNITE LOGIQUE'',I3)') 200+KCODE,KNUMER
        ENDIF
C
        ILMESA=LEN (CLMESA)
        ILMESU=ILMESA-1-2*LEN (CLJOLI)-ILNSPR-4
        CLMESA=CLJOLI//' '//CDNSPR(1:ILNSPR)//' - '//CLMESS(1:ILMESU)
     S         //CLJOLI
        WRITE (UNIT=FA%NULOUT,FMT='(A)') CLMESA
C
        WRITE (UNIT=FA%NULOUT,FMT=*) CLMESA
        IF (LDFATA.AND.KCODE.NE.0) THEN
          CALL SDL_SRLABORT
        ENDIF
C
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FAIPAR_MT',1,ZHOOK_HANDLE)
      END

