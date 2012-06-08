C**
C----- DESCRIPTION DES "PARAMETER" DU LOGICIEL DE "FICHIER ARPEGE" -----
C
C     JPNXFA = Nombre maximum de fichiers ouverts "simultanement"
C     JPNXCA =    "      "    "  cadres definissables "simultanement"
C     JPXNIV =    "      "    "  niveaux verticaux (champs d'altitude)
C     JPXTRO = Troncature maximum gerable
C     JPXLAT = Nombre maximum de latitudes de pole a pole
C     JPXLON = Nombre maximum de longitudes par parallele
C     JPLDAT = Longueur de l'article 'DATE', en mots
C     JPLB1P = Longueur du tableau "Bloc 1" pour sous-programmes GRIB
C     JPLB2P = Longueur du tableau "Bloc 2" pour sous-programmes GRIB
C     JPNIIL = Code "valeur absente" du logiciel pour les entiers
C     JPXCSP = Dimension maxi d'un champ en coefficients spectraux
C     JPXPDG = Dimension maxi d'un champ en points de grille
C     JPXCHA = Dimension maxi d'un champ ( maximum de JPXCSP et JPXPDG )
C     JPXPAH = Nombre maximum de latitudes par hemisphere
C     JPXIND = DIMENSIONING OF NOZPAR()
C     JPXGEO = DIMENSIONING OF SINLAT()
C     JPNVER = Numero de version du logiciel (qui est le contenu de
C              l'article dont le nom est l'identificateur du fichier)
C     JPUILA = Puissance de laplacien maximum pour laquelle les tableaux
C              servant a calculer laplacien et inverse sont precalcules
C     JPXNOM = Nombre maximum de caracteres par NOM d'article LFI.
C     JPXPRF =   "       "    "      "      par PReFixe de champ.
C     JPXSUF = JPXPRF+JPXNOM.
C     JPTNIV = Nombre de types de niveaux verticaux (re)connus.
C     CPDATE = Nom de l'article DATE
C
C         Noms des articles contenant les differentes parties du CADRE:
C
C     CPCADI = "Dimensions" (MTRONC, NNIVER, NLATIT, NXLOPA)
C     CPCAFS = Parametres de la transformation ARPEGE
C              (SSLAPO, SCLOPO, SSLOPO, SCODIL)
C     CPCARP = Tableaux lies a la reduction des points pres des poles
C     CPCASL = Tableau des sinus des latitudes
C     CPCACH = Valeurs des fonctions "A" et "B" de la coordonnee hybride
C     JPCADI et JPCAFS sont les longueurs des 2 premiers de ces articles
C
      INTEGER JPNXFA, JPNXCA, JPLDAT, JPNIIL, JPXNIV, JPXTRO, JPXLAT
      INTEGER JPUILA, JPXAU1, JPXLON, JPXAU2, JPXPAH, JPXIND, JPXGEO
      INTEGER JPXCSP, JPXCHA, JPLB1P, JPLB2P, JPCADI, JPCAFS, JPNVER
      INTEGER JPXPDG, JPXNOM, JPXPRF, JPXSUF, JPTNIV
C
      CHARACTER CPCADI*(*), CPCAFS*(*), CPCARP*(*), CPCACH*(*)
      CHARACTER CPCASL*(*), CPDATE*(*)
C
C
C Reglage de la troncature maximum gerable (JPXTRO)
C et du nombre maximum de niveaux verticaux (JPXNIV)
C
#if defined ( HIGHRES )
C
C     Setup high resolution parameters 
C
      PARAMETER ( JPXTRO=1280, JPXLAT=2560, JPXNIV=200 )
C
#else
C 
C     Setup low resolution parameters to save memory
C
      PARAMETER ( JPXTRO=599, JPXLAT=1200, JPXNIV=200 )
C
#endif
C
C
      PARAMETER ( JPNXFA=20, JPNXCA=20, JPLDAT=11, JPNIIL=-999 )
      PARAMETER ( JPUILA=3 )
      PARAMETER ( JPTNIV=122 )
      PARAMETER ( JPXAU1=(1+JPXLAT)/2, JPXLON=2*JPXLAT )
      PARAMETER ( JPXAU2=(2*JPXTRO)+4 )
      PARAMETER ( JPXPAH=(8*(8/JPXAU1)+JPXAU1*(JPXAU1/8))
     S                   /((8/JPXAU1)+(JPXAU1/8)) )
      PARAMETER ( JPXIND=(JPXAU1*(JPXAU1/JPXAU2)+JPXAU2*(JPXAU2/JPXAU1))
     S                   /((JPXAU1/JPXAU2)+(JPXAU2/JPXAU1)) )
      PARAMETER ( JPXGEO=(12*(12/JPXAU1)+JPXAU1*(JPXAU1/12))
     S                   /((12/JPXAU1)+(JPXAU1/12)) )
      PARAMETER ( JPXCSP=(1+JPXTRO)*(2+JPXTRO), JPXPDG=JPXLON*JPXLAT )
      PARAMETER ( JPXCHA=(JPXCSP*(JPXCSP/JPXPDG)+JPXPDG*(JPXPDG/JPXCSP))
     S                   /((JPXCSP/JPXPDG)+(JPXPDG/JPXCSP)) )
      PARAMETER ( JPLB1P=19, JPLB2P=17, JPXNOM=16, JPXPRF=8 )
      PARAMETER ( JPCADI=5, JPCAFS=4, JPNVER=1, JPXSUF=JPXNOM+JPXPRF )
      PARAMETER ( CPCADI='CADRE-DIMENSIONS', CPCAFS='CADRE-FRANKSCHMI' )
      PARAMETER ( CPCARP='CADRE-REDPOINPOL', CPCACH='CADRE-FOCOHYBRID' )
      PARAMETER ( CPCASL='CADRE-SINLATITUD', CPDATE='DATE-DES-DONNEES' )
C
