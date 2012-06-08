C Jan-2011 P. Marguinaud Interface to thread-safe FA
      INTEGER FUNCTION DECF10 (  KGRIB,  KLENG, KDECAL, KFAORI, KFAMOD,
     X                          KNBIMO )
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C     Fonction "DECF10"
C
C     But.
C     ---
C
C     Effectuer dans un message GRIB au niveau d'edition 1 le decalage
C     (en valeur) du descripteur "facteur d'echelle decimal".
C
C
C**   Interface.
C     ----------
C
C
C      Arguments d'appel:
C
C     KGRIB  (Entree+ ===> Tableau receptacle pour le message GRIB
C             Sortie)
C     KLENG  (Entree) ===> Longueur utile (mots) de KGRIB
C     KDECAL (Entree) ===> DECALage (additif) a effectuer sur la valeur
C                          courante du facteur d'echelle decimal.
C     KFAORI (Sortie) ===> FActeur decimal d'ORIgine
C     KFAMOD (Sortie) ===> FActeur decimal MODifie
C     KNBIMO (Entree) ===> Nombre de BIts par mot (entier)
C
C     Le code-retour de la fonction est:
C
C      0 si tout s'est bien passe
C     -1 si KGRIB ne commence pas par un entete GRIB
C     -2 si l'edition du GRIB ne convient pas
C     >0 code-retour (erreur) INXBIT
C
C
C     Methode.
C     --------
C
C     Verifie le niveau d'edition du message GRIB, decode le facteur
C     d'echelle decimal courant, y ajoute KDECAL, et remplace la valeur
C     par insertion dans le message.
C
C
C     Sous-programmes/fonctions externes utilises.
C     --------------------------------------------
C
C     INXBIT - INsertion d'une chaine binaire
C
C     Reference.
C     ----------
C
C     Aucune.
C
C
C     Commentaires.
C     -------------
C
C     Cettte fonction a ete developpee dans le cadre du pretraitement de
C     modeles etrangers arrivant sous forme de fichiers, pour DIAPASON
C     et SYNERGIE. Mais elle n'en depend pas vraiment, et peut donc
C     etre reutilisee dans un autre contexte.
C
C     Dans sa mouture initiale, seule l'edition 1 de GRIB est geree.
C
C     On utilise une primitive de bas niveau de la bibliotheque
C     "emos" du CEPMMT.
C
C     La convention de codage du signe du facteur d'echelle est
C     explicitement geree par la fonction courante.
C
C
C     Auteur.
C     -------
C
C     J. Clochard, Meteo France, DSI/OP/D - Juin 2001.
C
C
C     Modifications.
C     --------------
C
C     Aucune.
C
C
C
C     Constantes symboliques.
C
C     Signification des PARAMETER "locaux" au sous-programme :
C
C     JPDGRB => Nombre d'octets a decoder pour controler le message GRIB
C     JPGRIB => Nombre d'octets a decoder pour controler l'entete GRIB
C
      INTEGER, PARAMETER :: JPDGRB=8, JPGRIB=4
C
C     Arguments d'appel.
C
      INTEGER KLENG, KDECAL, KFAORI, KFAMOD, KNBIMO
C
      INTEGER KGRIB (KLENG)
C
C     Variables locales.
C
      INTEGER IAUXIL, IREPON, J, IDECAL, ILODES, IPIVOT, IFACOD, IFAC10
C
      INTEGER IDGRIB (JPDGRB)
C
C     Caracteres G, R, I, B en code CCITT IA-5
C
      INTEGER, PARAMETER :: IBLOCD (JPGRIB) = (/ 71, 82, 73, 66 /)
C
      CHARACTER CLOPER*1
C**
C     1.  -  INITIALISATIONS ET CONTROLES.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('DECF10',0,ZHOOK_HANDLE)
      IREPON = 0
      KFAORI = -9999
      KFAMOD = KFAORI
C
      IF ( KDECAL .EQ. 0 ) THEN
C
        GO TO 900
C
      ENDIF
C
C     Deocdage des premiers octets du message GRIB.
C
      CLOPER = 'D'
      IDECAL = 0
      IAUXIL = JPDGRB
      ILODES = 8
C
C      CALL INXBIT(  KGRIB,  KLENG, IDECAL, IDGRIB, IAUXIL, KNBIMO,
C     X             ILODES, CLOPER, IREPON )
C
      IF ( IREPON .NE. 0 ) THEN
C
        GO TO 900
C
      ENDIF
C
C     KGRIB commence-t-il bien par les 4 lettres GRIB ?
C
      DO 101 J = 1, JPGRIB
C
        IF( IDGRIB(J) .NE. IBLOCD(J) ) THEN
C
          IREPON = -1
C
          GO TO 900
C
        ENDIF
C
  101 CONTINUE
C
C     Controle du niveau d'edition de GRIB.
C
      IF( IDGRIB(JPDGRB) .NE. 1 ) THEN
C
        IREPON = -2
C
        GO TO 900
C
      ENDIF
C
      IDECAL = (8+26)*8
      IAUXIL = 1
      ILODES = 16
      IPIVOT = 2**15
C**
C     2.  -  DECODAGE DU FACTEUR D'ECHELLE DECIMAL COURANT.
C-----------------------------------------------------------------------
C
  200 CONTINUE
C
C      CALL INXBIT(  KGRIB,  KLENG, IDECAL, IFACOD, IAUXIL, KNBIMO,
C     X             ILODES, CLOPER, IREPON )
      IDECAL = IDECAL-ILODES*IAUXIL
C
      IF ( IREPON .NE. 0 ) THEN
C
        GO TO 900
C
      ELSEIF ( IFACOD .LE. IPIVOT ) THEN
C
        IFAC10 = IFACOD
C
      ELSE
C
        IFAC10 = IPIVOT-IFACOD
C
      ENDIF
C
      KFAORI = IFAC10
C**
C     3.  -  MODIFICATION DU FACTEUR D'ECHELLE DECIMAL COURANT.
C-----------------------------------------------------------------------
C
  300 CONTINUE
C
      IFAC10 = IFAC10+KDECAL
      KFAMOD = IFAC10
      CLOPER = 'C'
C
      IF ( IFAC10 .GE. 0 ) THEN
C
        IFACOD = IFAC10
C
      ELSE
C
        IFACOD = IPIVOT-IFAC10
C
      ENDIF
C
C      CALL INXBIT(  KGRIB,  KLENG, IDECAL, IFACOD, IAUXIL, KNBIMO,
C     X             ILODES, CLOPER, IREPON )
C**
C     9.  -  MESSAGERIE EVENTUELLE, RETOUR A L'APPLICATIF APPELANT.
C-----------------------------------------------------------------------
C
  900 CONTINUE
C
      DECF10 = IREPON
C
      IF (LHOOK) CALL DR_HOOK('DECF10',1,ZHOOK_HANDLE)
      END

