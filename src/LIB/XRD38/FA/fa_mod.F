C Jan-2011 P. Marguinaud Interface to thread-safe FA
      MODULE FA_MOD

      USE PARKIND1, ONLY : JPIM
      USE LFIMOD, ONLY : LFICOM

      IMPLICIT NONE
#include "precision.h"
      INTEGER, PARAMETER :: JPXNOM = 16
      INTEGER, PARAMETER :: JPXPRF = 8
      INTEGER, PARAMETER :: JPXSUF = JPXNOM+JPXPRF 
C****
C           !-----------------------------------------------!
C           ! Sous-programme du Logiciel de Fichiers ARPEGE !
C           !-----------------------------------------------!
C
C       - Version originale du logiciel: Mars 1990, auteur:
C           Jean CLOCHARD, Meteorologie Nationale, FRANCE.
C
C       - Juin 90: Ajout du "type de transformation horizontale"
C         aux elements definissant un cadre
C         (1=transformation de Frank-Schmidt "pure", 2=ARPEGE complete),
C         et ajout de la notion d'"Identificateur de Fichier".
C
C**
C---------------- VARIABLES "SIMPLES" GLOBALES -------------------------
C
C     NFIOUV = Nombre de fichiers ouverts simultanement
C     NCADEF =   "    "  cadres definis        "
C     LFAMUL = Option "utilisation du logiciel en mode multi-taches"
C     LFAMOP = Option "mode mise au point du logiciel"
C     VRGLAS = Verrou global du lociciel (en mode multi-taches)
C     NIMSGA = Niveau global de messagerie du logiciel
C     NRFAGA = Niveau de filtrage global des erreurs fatales detectees ?
C     NBIPDG = Nombre de bits IMPLICITE par valeur point-de-grille
C     NBICSP =    "   "   "   IMPLICITE par coefficient spectral
C     NSTROI = Sous-troncature non compactee IMPLICITE (coef. spectraux)
C     NPUILA = "Puissance de laplacien" IMPLICITE pour coeffi. spectraux
C     NMIDPL = Degre de modulation IMPLICITE de puissance de laplacien
C     NIGRIB = Niveau IMPLICITE de codage GRIB
C              (  0 ==> pas de codage, 1 ==> "standard", 2 => "Arpege" )
C              ( -1 ==> pas de codage (rangt coeff spec selon modele)  )
C              (  3 ==> "GRIBEX"      (rangt coeff spec selon modele)  )
C     NCPCAD = Nombre maximum de caracteres possibles par nom de cadre.
C     SPSMIN = Pression sol minimum   "     "    "     "      "    "
C     SPSMAX = Pression sol maximum   "     "    "     "      "    "
C     MPRESX = Pression maximum pour designer un niveau isobare
C     NBIMAC = Nombre de bits par mot machine
C     NBIMAX = Nombre maximum de bits par valeur elementaire compactee
C     NTYPTX = Type maximum reconnu de transformation horizontale
C     LIGARD = Option de conservation des cadres definis dynamiquement
C
C        Ci-dessous, limites en termes de dimensions reglables au niveau
C     usager, initialisees aux valeurs maximum du logiciel:
C
C     NXNIVV = Nombre maximum de niveaux verticaux au niveau usager
C     NXTRON = Troncature maximum                  "    "      "
C     NXLATI = Nombre maximum de latitudes de pole a pole "    "
C     NXLONG = Nombre maximum de longitudes par parallele "    "
C
C     CHAINC = Impression de substitut aux variables non reconnues
C**
C------------------------- TABLEAUX GLOBAUX ----------------------------
C
C     NULIND = Table d'indirection pour les autres tables "fichiers"
C     NCAIND =   "   "      "       "    "     "      "   "cadres"
C**
C--------- DESCRIPTIF DES ELEMENTS CONCERNANT UNE UNITE LOGIQUE --------
C
C     NULOGI = Numero de l'unite logique
C     NUCADR = Rang du cadre auquel le fichier est rattache
C     CIDENT = Identification (utilisateur) du fichier
C     LNOMME = Vrai si l'unite logique a un nom utilisateur (BOF...)
C     NIVOMS = Niveau de la messagerie defini par l'utilisateur
C     LERRFA = Vrai si toute erreur detectee doit etre fatale
C     LCREAF = Vrai si le fichier est en mode creation
C     NBFPDG = Nombre de bits EFFECTIF par valeur point-de-grille
C     NBFCSP =    "   "   "   EFFECTIF par coefficient spectral
C     NPUFLA = "Puissance de laplacien" EFFECTIVE pour coeffi. spectraux
C     NSTROF = Sous-troncature non compactee EFFECTIVE (coef. spectraux)
C     NMFDPL = Degre de modulation EFFECTIF de la puissance de laplacien
C     NFGRIB = Niveau EFFECTIF de codage GRIB
C     VRFICH = Verrou du fichier (en mode multi-taches)
C
C     MADATE = Contenu de l'article 'DATE'
C
C     CTNPRF = Prefixes de champ (re)connus.
C     NIVDSC = Descripteurs associes aux prefixes (re)connus.
C     NRASHO = Nombre de champs de coeff. spectraux ARPEGE ranges
C              horizontalement (donc a l'inverse du modele) detectes au
C              cours des lectures/ecritures
C     NRASVE = Nombre de champs de coeff. spectraux ARPEGE ranges
C              verticalement (donc comme le modele), detectes au cours
C              des lectures/ecritures
C
C**
C----------- DESCRIPTIF DES ELEMENTS CONSTITUANT UN CADRE --------------
C
C                     Elements DEFINISSANT un cadre:
C
C     CNOMCA = Nom du cadre (conventionnel; non ecrit sur le fichier)
C     NLCCAD = Longueur en caracteres du nom de cadre
C     NULCAD = Nombre d'unites logiques rattachees au cadre
C     MTRONC = Troncature associee aux donnees contenues dans le fichier
C     NNIVER = Nombre de niveaux verticaux
C     NLATIT =   "    "  latitudes de pole a pole
C     NXLOPA =   "    maximum de longitudes par cercle de latitude
C     NTYPTR = Type de transformation horizontale
C     SSLAPO = Sinus   de la latitude  du pole d'interet
C     SCLOPO = Cosinus  "   "     "     "    "      "
C     SSLOPO = Sinus    "   "     "     "    "      "
C     SCODIL = Coefficient de dilatation
C     SPREFE = Pression de reference (facteur multiplicatif de la premi-
C              ere fonction de la coordonnee hybride)
C     NLOPAR = Nombre de longitudes par parallele de l'hemisphere "nord"
C     NOZPAR =    "   d'onde zonal maximum "    "  "     "        "
C     SINLAT = Sinus des latitudes de la grille (*hemisphere nord seul*)
C     SFOHYB = Valeurs des fonctions "A" et "B" de la coordonnee hybride
C              (Definies aux LIMITES DE COUCHES)
C
C                Autres elements rattaches a un cadre:
C
C     NVAPDG = Nombre total de points de grille dans un champ horizontal
C     NGARDE = Option de conservation lors de la fermeture du dernier
C              fichier s'y rattachant. ( 0=Non, 1=selon LIGARD, 2=Oui )
C     LIMLAM = IMPLICIT SWITCH FOR SETTING OF GLOBAL/LAM CASE
C     NSFLAM = TOTAL NUMBER OF SPECTRAL COEFFICIENTS FOR LAM CASE
C
      TYPE FA_COM

      TYPE(LFICOM), POINTER :: LFI => NULL ()

      INTEGER NFIOUV, NCADEF, NIMSGA, NRFAGA, NBIPDG, NBICSP, NPUILA
      INTEGER NIGRIB, NCPCAD, NSTROI, NMIDPL, NBIMAC, NBIMAX, MPRESX
      INTEGER NXNIVV, NXTRON, NXLATI, NXLONG, NTYPTX
C
      INTEGER, POINTER :: NULIND (:), NCAIND (:)
      INTEGER, POINTER :: NULOGI (:), NUCADR (:), NIVOMS (:)
      INTEGER, POINTER :: NBFPDG (:), NBFCSP (:), NPUFLA (:)
      INTEGER, POINTER :: NFGRIB (:), NSTROF (:), NMFDPL (:)
      INTEGER, POINTER :: NRASHO (:), NRASVE (:)
      INTEGER, POINTER :: MADATE (:,:)
      INTEGER, POINTER :: MTRONC (:), NNIVER (:), NLATIT (:)
      INTEGER, POINTER :: NXLOPA (:), NULCAD (:), NLCCAD (:)
      INTEGER, POINTER :: NLOPAR (:,:), NOZPAR (:,:)
      INTEGER, POINTER :: NVAPDG (:), NTYPTR (:), NGARDE (:)
      INTEGER, POINTER :: NSFLAM (:), NIVDSC (:,:)
C
      REAL (KIND=JPDBLR) SPSMIN, SPSMAX, VRGLAS
C
      REAL (KIND=JPDBLR), POINTER :: SSLAPO (:), SCLOPO (:)
      REAL (KIND=JPDBLR), POINTER :: SCODIL (:), SINLAT (:,:)
      REAL (KIND=JPDBLR), POINTER :: SFOHYB (:,:,:), SPREFE (:)
      REAL (KIND=JPDBLR), POINTER :: VRFICH (:), SSLOPO (:)
C
      LOGICAL LFAMUL, LFAMOP, LIGARD
      LOGICAL, POINTER :: LNOMME (:), LERRFA (:), LCREAF (:)
      LOGICAL, POINTER :: LIMLAM (:)
C
      CHARACTER*(JPXNOM), POINTER :: CNOMCA (:), CIDENT (:)
      CHARACTER*(JPXNOM) CHAINC
      CHARACTER, POINTER :: CTNPRF (:)*8
C
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
      CHARACTER CPCADI*(16), CPCAFS*(16), CPCARP*(16), CPCACH*(16)
      CHARACTER CPCASL*(16), CPDATE*(16)
C


!*
!      FAMODU - MODULE POUR LE LOGICIEL FA

!        D. PARADIS       METEO FRANCE     21/7/00

!  ----------------------------------------------------------------------
!  1 - COEFFICIENTS MULTIPLICATIFS DES COEFF SPECTRAUX
!
!     XLAP1D  = Coefficients elementaires pour "laplacien" [n*(n+1)],
!              ainsi que leurs inverses.
!     XLAP2D  = Cf. XLAP1D, mais etale dans un champ complet, et calcules
!              pour les puissances 1 a JPUILA, ainsi que les inverses.
!     FLAP1D  = Tableau des "puissances de laplacien" pour l'ecriture
!                     (n*(n+1))**NPUFLA()
!              ( Ce dernier tableau n'est calcule que pour une
!                puissance de laplacien non nulle )
!     XLAP1DA = ALADIN version de XLAP1D
!     XLAP2DA = ALADIN version de XLAP2D
!     FLAP1DA = ALADIN version de FLAP1D
!     LIXLAP  = .T. s'il faut initialiser XLAPxDx
!     LIFLAP() = .T. s'il faut initialiser FLAP1Dx()
!
!  2 - EN-TETE POUR GRIBEX (SECTIONS 1 ET 2)
!
!     JPSEC1   = taille du tableau contenant les elements de la section 1
!     JPSEC2   = taille du tableau contenant les elements de la section 2
!     JPSEC4   = taille du tableau contenant les elements de la section 4
!     NSEC1    = tableau contenant les elts 2:21 de la section 1 de GRIBEX
!     LISEC1() = .T. s'il faut initialiser le tableau NSEC1 ci-dessus
!
!      Tableaux contenant des elements de la section 2 de GRIBEX:
!     NSEC2SP  = cas de la representation spectrale ARPEGE avec niveaux modele
!     NSEC2GG  = cas de la grille de Gauss
!     NSEC2LL  = cas de la grille latitude-longitude
!     NSEC2LA  = cas de la grille Lambert conforme (type general de Aladin)
!     NSEC2AL  = cas de la representation spectrale Aladin deguisee lat-lon
!     NSC2ALF  = supplement a NSEC2AL dependant du fichier
!     XSEC2    = coordonnee verticale
!     LISEC2() = .T. s'il faut initialiser les tableaux ci-dessus sauf NSC2ALF
!     LISC2F() = .T. s'il faut initialiser NSC2ALF
!
!  3 - PARAMETRES DEFINISSANT LE CODAGE GRIBEX
!
!     NCODGRI  = tableau contenant les parametres de codage implicites
!     NCOGRIF  = tableau contenant les parametres de codage pour chaque fichier
!     CIPREF   = tableau contenant les prefixes des noms des champs connus
!     CISUFF   = tableau contenant les suffixes des noms des champs connus
!     NCODPA   = tableau contenant les 6 descripteurs GRIB des champs connus
!                NCODPA(i,1)= numero de version de la table de code parametre,
!                             KSEC1(1)
!                NCODPA(i,2)= indicateur de parametre, KSEC1(6)
!                NCODPA(i,3)= indicateur de type de niveau, KSEC1(7)
!                NCODPA(i,4)= niveau (premier niveau de la couche), KSEC1(8)
!                NCODPA(i,5)= deuxieme niveau de la couche, KSEC1(9)
!                NCODPA(i,6)= indicateur de type de champ, KSEC1(18)
!                             (0 sauf si min/max dans le temps => 2,
!                                  ou si cumul dans le temps   => 4 )
!     JPXPAR   = nombre maximal de champs connus pour CISUFF, CIPREF et NCODPA
!     NBPARC   = nombre de champs connus pour CISUFF, CIPREF et NCODPA
!
!  4 - TABLEAU POUR CODER LES COEF SPECTRAUX ALADIN AVEC GRIBEX
!
!     NOMPAR   = tableau decrivant la position des coeff spectraux (CSP) ALADIN
!                en fonction du nombre d'onde zonal dans un tableau de donnees
!                stocke dans FA avec la methode de codage 3 ou -1 (rangt modele).
!                NOMPAR est l'equivalent de NOZPAR mais pour un rangement vertical
!                des CSP: NOMPAR(1)=NSMAX, NOMPAR(2)=NMSMAX (ces 2 valeurs de
!                NOMPAR sont speciales commme pour NOZPAR) et pour un JM donne,
!                compris entre 0 et NMSMAX, NOMPAR(2*JM+3) donne l'indice du
!                premier CSP associe a JM qui est contenu dans un champ spectral
!                et NOMPAR(2*JM+4) donne l'indice du dernier CSP associe a JM
!                qui est contenu dans un champ spectral.
!                
!  ----------------------------------------------------------------------

!  1 - COEFFICIENTS MULTIPLICATIFS DES COEFF SPECTRAUX

      REAL (KIND=JPDBLR), POINTER :: XLAP1D(:,:),   XLAP1DA(:,:)
      REAL (KIND=JPDBLR), POINTER :: XLAP2D(:,:,:), XLAP2DA(:,:,:)
      REAL (KIND=JPDBLR), POINTER :: FLAP1D(:,:),   FLAP1DA(:,:) 

      LOGICAL LIXLAP
      LOGICAL, POINTER :: LIFLAP(:)
 
!  2 - EN-TETE POUR GRIBEX (SECTIONS 1 ET 2)

      INTEGER :: JPSEC1, JPSEC2
      INTEGER :: JPSEC4
      INTEGER, POINTER :: NSEC1(:,:), NSEC2SP(:,:)
      INTEGER, POINTER :: NSEC2LL(:,:)
      INTEGER, POINTER :: NSEC2GG(:,:)
      INTEGER, POINTER :: NSEC2LA(:,:)
      INTEGER, POINTER :: NSEC2AL(:,:), NSC2ALF(:,:)

      REAL (KIND=JPDBLR), POINTER :: XSEC2(:,:)

      LOGICAL, POINTER :: LISEC1(:), LISEC2(:), LISC2F(:)
 
!  3 - PARAMETRES DEFINISSANT LE CODAGE GRIBEX

      INTEGER, POINTER :: NCODGRI(:), NCOGRIF(:,:)
      INTEGER JPXPAR
      INTEGER NBPARC
      CHARACTER(LEN=JPXPRF), POINTER :: CIPREF(:)
      CHARACTER(LEN=JPXSUF), POINTER :: CISUFF(:)
      INTEGER, POINTER :: NCODPA(:,:)

!  4 - TABLEAU POUR CODER LES COEF SPECTRAUX ALADIN AVEC GRIBEX

      INTEGER, POINTER :: NOMPAR(:,:)

      INTEGER (KIND=JPDBLE), POINTER :: ICHAMP(:), ICHAUX(:)

      LOGICAL :: FACADE_LLPREA = .TRUE.
      LOGICAL :: FACAGE_LLPREA = .TRUE.
      LOGICAL :: FACIES_LLPREA = .TRUE.
      LOGICAL :: FACTUM_LLPREA = .TRUE.
      LOGICAL :: FAGIOT_LLPREA = .TRUE.
      LOGICAL :: FAIFLA_LLPREA = .TRUE.
      LOGICAL :: FALIMU_LLPREA = .TRUE.
      LOGICAL :: FAMISO_LLPREA = .TRUE.
      LOGICAL :: FANERG_LLPREA = .TRUE.
      LOGICAL :: FANMSG_LLPREA = .TRUE.
      LOGICAL :: FANUCA_LLPREA = .TRUE.
      LOGICAL :: FANUMU_LLPREA = .TRUE.
      LOGICAL :: FAREGI_LLPREA = .TRUE.
      LOGICAL :: FARFLU_LLPREA = .TRUE.
      LOGICAL :: FARINE_LLPREA = .TRUE.
      LOGICAL :: FAVORI_LLPREA = .TRUE.
      LOGICAL :: FAXION_LLPREA = .TRUE.
      LOGICAL :: FARINE_LLDEFM = .FALSE.
      LOGICAL, POINTER :: FATRAN_LLINDIR(:)
      INTEGER, POINTER :: FATRAN_INDIRECT(:,:)
      INTEGER :: FAXION_ISCALX
      REAL (KIND=JPDBLR) FAXION_ZEPSIL

      INTEGER :: NULOUT = 0
      LOGICAL :: LOPENMP = .TRUE.

      END TYPE

      TYPE(FA_COM), SAVE, TARGET :: FA_COM_DEFAULT
      LOGICAL, SAVE :: FA_COM_DEFAULT_INIT = .FALSE.

      CONTAINS

      SUBROUTINE NEW_FA_DEFAULT ()
      USE LFIMOD, ONLY : LFICOM_DEFAULT, NEW_LFI_DEFAULT
      INTEGER(KIND=JPIM) :: IERR

      CALL NEW_LFI_DEFAULT
      IF (.NOT. FA_COM_DEFAULT_INIT) THEN
        CALL NEW_FA (FA_COM_DEFAULT, IERR)
        FA_COM_DEFAULT_INIT = .TRUE.
        FA_COM_DEFAULT%LFI => LFICOM_DEFAULT
      ENDIF

      END SUBROUTINE

      SUBROUTINE NEW_FA (FA, KERR, KPXTRO, KPXLAT, 
     S                   KPXNIV, KPNXFA, KPNXCA)
      TYPE(FA_COM) :: FA
      INTEGER(KIND=JPIM), INTENT(OUT) :: KERR
      INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KPXTRO
      INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KPXLAT
      INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KPXNIV
      INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KPNXFA
      INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN) :: KPNXCA


      FA%JPXNOM = JPXNOM
      FA%JPXPRF = JPXPRF
      FA%JPXSUF = JPXSUF

C
C Reglage de la troncature maximum gerable (JPXTRO)
C et du nombre maximum de niveaux verticaux (JPXNIV)
C
#if defined ( HIGHRES )
C
C     Setup high resolution parameters 
C
      FA%JPXTRO = 1280
      FA%JPXLAT = 2560 
      FA%JPXNIV = 200 
C
#else
C 
C     Setup low resolution parameters to save memory
C
      FA%JPXTRO = 599 
      FA%JPXLAT = 1200 
      FA%JPXNIV = 200 
C
#endif
C
C
      FA%JPNXFA=20 
      FA%JPNXCA=20 

      IF (PRESENT (KPXTRO)) FA%JPXTRO = KPXTRO
      IF (PRESENT (KPXLAT)) FA%JPXLAT = KPXLAT
      IF (PRESENT (KPXNIV)) FA%JPXNIV = KPXNIV
      IF (PRESENT (KPNXFA)) FA%JPNXFA = KPNXFA
      IF (PRESENT (KPNXCA)) FA%JPNXFA = KPNXCA

      FA%JPLDAT=11 
      FA%JPNIIL=-999 
      FA%JPUILA=3 
      FA%JPTNIV=12 
      FA%JPXAU1=(1+FA%JPXLAT)/2 
      FA%JPXLON=2*FA%JPXLAT 
      FA%JPXAU2=(2*FA%JPXTRO)+4 
      FA%JPXPAH=(8*(8/FA%JPXAU1)+FA%JPXAU1*(FA%JPXAU1/8))
     S          /((8/FA%JPXAU1)+(FA%JPXAU1/8)) 
      FA%JPXIND=(FA%JPXAU1*(FA%JPXAU1/FA%JPXAU2)+FA%JPXAU2*
     S          (FA%JPXAU2/FA%JPXAU1))
     S          /((FA%JPXAU1/FA%JPXAU2)+(FA%JPXAU2/FA%JPXAU1)) 
      FA%JPXGEO=(12*(12/FA%JPXAU1)+FA%JPXAU1*(FA%JPXAU1/12))
     S          /((12/FA%JPXAU1)+(FA%JPXAU1/12)) 
      FA%JPXCSP=(1+FA%JPXTRO)*(2+FA%JPXTRO) 
      FA%JPXPDG=FA%JPXLON*FA%JPXLAT 
      FA%JPXCHA=(FA%JPXCSP*(FA%JPXCSP/FA%JPXPDG)+
     S          FA%JPXPDG*(FA%JPXPDG/FA%JPXCSP))
     S          /((FA%JPXCSP/FA%JPXPDG)+(FA%JPXPDG/FA%JPXCSP)) 
      FA%JPLB1P=19 
      FA%JPLB2P=17 
      FA%JPCADI=5 
      FA%JPCAFS=4 
      FA%JPNVER=1 
      FA%CPCADI='CADRE-DIMENSIONS' 
      FA%CPCAFS='CADRE-FRANKSCHMI'
      FA%CPCARP='CADRE-REDPOINPOL' 
      FA%CPCACH='CADRE-FOCOHYBRID'
      FA%CPCASL='CADRE-SINLATITUD' 
      FA%CPDATE='DATE-DES-DONNEES' 

      FA%JPSEC1=37
      FA%JPSEC2=22+MAX(FA%JPXTRO-1,FA%JPXLAT)
      FA%JPSEC4=42
      FA%JPXPAR=500

      FA%XLAP1D  => NULL () 
      FA%XLAP1DA => NULL ()
      FA%XLAP2D  => NULL () 
      FA%XLAP2DA => NULL ()
      FA%FLAP1D  => NULL () 
      FA%FLAP1DA => NULL ()
      FA%ICHAMP  => NULL ()
      FA%ICHAUX  => NULL ()

      ALLOCATE (
     S  FA%NULIND (FA%JPNXFA),   FA%NCAIND (FA%JPNXFA), 
     S  FA%NULOGI (FA%JPNXFA),   FA%NUCADR (FA%JPNXFA), 
     S  FA%NIVOMS (0:FA%JPNXFA), FA%NBFPDG (FA%JPNXFA), 
     S  FA%NBFCSP (FA%JPNXFA),   FA%NPUFLA (FA%JPNXFA),
     S  FA%NFGRIB (FA%JPNXFA),   FA%NSTROF (FA%JPNXFA), 
     S  FA%NMFDPL (FA%JPNXFA),   FA%NRASHO (FA%JPNXFA), 
     S  FA%NRASVE (FA%JPNXFA),   FA%MADATE (FA%JPLDAT,FA%JPNXFA),
     S  FA%MTRONC (FA%JPNXCA),   FA%NNIVER (FA%JPNXCA), 
     S  FA%NLATIT (FA%JPNXCA),   FA%NXLOPA (FA%JPNXCA),   
     S  FA%NULCAD (FA%JPNXCA),   FA%NLCCAD (FA%JPNXCA),
     S  FA%NLOPAR (FA%JPXPAH,FA%JPNXCA), 
     S  FA%NOZPAR (FA%JPXIND,FA%JPNXCA),
     S  FA%NVAPDG (FA%JPNXCA), FA%NTYPTR (FA%JPNXCA), 
     S  FA%NGARDE (FA%JPNXCA),
     S  FA%NSFLAM (FA%JPNXCA), FA%NIVDSC (0:4,0:FA%JPTNIV),
     S  FA%SSLAPO (FA%JPNXCA), FA%SCLOPO (FA%JPNXCA),
     S  FA%SCODIL (FA%JPNXCA), FA%SINLAT (FA%JPXGEO,FA%JPNXCA),
     S  FA%SFOHYB (2,0:FA%JPXNIV,FA%JPNXCA), 
     S  FA%SPREFE (FA%JPNXCA),
     S  FA%VRFICH (FA%JPNXFA), FA%SSLOPO (FA%JPNXCA),
     S  FA%LNOMME (FA%JPNXFA), FA%LERRFA (0:FA%JPNXFA), 
     S  FA%LCREAF (FA%JPNXFA), FA%LIMLAM (FA%JPNXCA),
     S  FA%CNOMCA (FA%JPNXCA), FA%CIDENT (FA%JPNXFA),
     S  FA%CTNPRF (FA%JPTNIV), FA%LIFLAP (FA%JPNXFA),
     S  STAT = KERR )
      IF (KERR /= 0) RETURN

 
      ALLOCATE (
     S FA%NSEC1(2:21,FA%JPNXFA), 
     S FA%NSEC2SP(22,FA%JPNXCA),
     S FA%NSEC2LL(22,FA%JPNXCA),
     S FA%NSEC2GG(22+FA%JPXLAT,FA%JPNXCA),
     S FA%NSEC2LA(22,FA%JPNXCA),
     S FA%NSEC2AL(22,FA%JPNXCA), 
     S FA%NSC2ALF(FA%JPXTRO-1,FA%JPNXFA),
     S FA%XSEC2(10+2*(FA%JPXNIV+1),FA%JPNXCA),
     S FA%LISEC1(FA%JPNXFA), 
     S FA%LISEC2(FA%JPNXCA), 
     S FA%LISC2F(FA%JPNXFA),
     S STAT = KERR )
      IF (KERR /= 0) RETURN
 
      ALLOCATE (
     S FA%NCODGRI(10), 
     S FA%NCOGRIF(10,FA%JPNXFA),
     S FA%CIPREF(FA%JPXPAR),
     S FA%CISUFF(FA%JPXPAR),
     S FA%NCODPA(FA%JPXPAR,6),
     S FA%NOMPAR(2*FA%JPXTRO+4,FA%JPNXCA),
     S STAT = KERR )
      IF (KERR /= 0) RETURN

      ALLOCATE (
     S FA%FATRAN_LLINDIR(FA%JPNXCA), 
     S FA%FATRAN_INDIRECT((FA%JPXTRO+1)*(FA%JPXTRO+2),FA%JPNXCA), 
     S STAT = KERR)
      IF (KERR /= 0) RETURN

      FA%FATRAN_LLINDIR = .TRUE.

      END SUBROUTINE

      SUBROUTINE FREE_FA (FA, KERR)
      TYPE(FA_COM) :: FA
      INTEGER(KIND=JPIM), INTENT(OUT) :: KERR

      IF (ASSOCIATED (FA%XLAP1D  )) DEALLOCATE (FA%XLAP1D  )
      IF (ASSOCIATED (FA%XLAP1DA )) DEALLOCATE (FA%XLAP1DA )
      IF (ASSOCIATED (FA%XLAP2D  )) DEALLOCATE (FA%XLAP2D  )
      IF (ASSOCIATED (FA%XLAP2DA )) DEALLOCATE (FA%XLAP2DA )
      IF (ASSOCIATED (FA%FLAP1D  )) DEALLOCATE (FA%FLAP1D  )
      IF (ASSOCIATED (FA%FLAP1DA )) DEALLOCATE (FA%FLAP1DA )
      IF (ASSOCIATED (FA%ICHAMP  )) DEALLOCATE (FA%ICHAMP  )
      IF (ASSOCIATED (FA%ICHAUX  )) DEALLOCATE (FA%ICHAUX  )

      DEALLOCATE (
     S  FA%NULOGI, FA%NUCADR, 
     S  FA%NIVOMS, FA%NBFPDG, 
     S  FA%NBFCSP, FA%NPUFLA,
     S  FA%NFGRIB, FA%NSTROF, 
     S  FA%NMFDPL, FA%NRASHO, 
     S  FA%NRASVE, FA%MADATE,
     S  FA%MTRONC, FA%NNIVER, 
     S  FA%NLATIT, FA%NXLOPA,   
     S  FA%NULCAD, FA%NLCCAD,
     S  FA%NLOPAR, 
     S  FA%NOZPAR,
     S  FA%NVAPDG, FA%NTYPTR, 
     S  FA%NGARDE,
     S  FA%NSFLAM, FA%NIVDSC,
     S  FA%SSLAPO, FA%SCLOPO,
     S  FA%SCODIL, FA%SINLAT,
     S  FA%SFOHYB, 
     S  FA%SPREFE,
     S  FA%VRFICH, FA%SSLOPO,
     S  FA%LNOMME, FA%LERRFA, 
     S  FA%LCREAF, FA%LIMLAM,
     S  FA%CNOMCA, FA%CIDENT,
     S  FA%CTNPRF, FA%LIFLAP,
     S  STAT = KERR )
      IF (KERR /= 0) RETURN

 
      DEALLOCATE (
     S FA%NSEC1, 
     S FA%NSEC2SP,
     S FA%NSEC2LL,
     S FA%NSEC2GG,
     S FA%NSEC2LA,
     S FA%NSEC2AL, 
     S FA%NSC2ALF,
     S FA%XSEC2,
     S FA%LISEC1, 
     S FA%LISEC2, 
     S FA%LISC2F,
     S STAT = KERR )
      IF (KERR /= 0) RETURN
 
      DEALLOCATE (
     S FA%NCODGRI, 
     S FA%NCOGRIF,
     S FA%CIPREF,
     S FA%CISUFF,
     S FA%NCODPA,
     S FA%NOMPAR,
     S STAT = KERR )
      IF (KERR /= 0) RETURN

      DEALLOCATE (
     S FA%FATRAN_LLINDIR,
     S FA%FATRAN_INDIRECT,
     S STAT = KERR )
      IF (KERR /= 0) RETURN

      END SUBROUTINE

      END MODULE FA_MOD

