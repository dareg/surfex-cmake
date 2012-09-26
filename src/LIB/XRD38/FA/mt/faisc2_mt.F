C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FAISC2_MT (FA,  KREP, KRANGC )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C      Ce sous-programme initialise des tableaux "reference" de
C      l'en-tete GRIB, section 2: les differents types de grille
C      sont abordes (routine appelee une seule fois pour un cadre donne)
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KRANGC (Entree) ==> Rang dans la table des cadres;
C*
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KRANGC
C
      REAL ZPI, ZRAMDE, ZLATPRE, ZLATDER, ZLONPRE, ZLONDER
      REAL (KIND=JPDBLR), PARAMETER :: ONE = 1.0
C
      INTEGER INLAT, INIVAU, INUMER, INIMES
C
      LOGICAL LLMLAM
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     0.  -  INITIALISATIONS PREALABLES
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FAISC2_MT',0,ZHOOK_HANDLE)
      KREP=0
      IF (KRANGC.LE.0.OR.KRANGC.GT.FA%JPNXCA) THEN
        KREP=-66
        GOTO 1001
      ENDIF
      INLAT  = FA%NLATIT(KRANGC)
      INIVAU = FA%NNIVER(KRANGC)
      LLMLAM = FA%LIMLAM(KRANGC)
      ZPI = 2.*ASIN(1.)
C Conversion des radians en 1/1000 de degre
      ZRAMDE = 180000./ZPI
C
      IF (LLMLAM) GOTO 300
C**
C     1.  -  KSEC2 POUR LA REPRESENTATION SPECTRALE ARPEGE
C-----------------------------------------------------------------------
C
C Type de representation de donnees
C
C FA%SSLAPO=sinus latitude du pole d'interet
C (si=1, pole=pole N et pas de rotation)
C FA%SCODIL=coeff de dilation (si =1, pas de dilatation)
      IF ((1.-FA%SSLAPO(KRANGC)).LE.1.E-10) THEN
        FA%NSEC2SP(1,KRANGC)=70
        IF (ABS(FA%SCODIL(KRANGC)-1.).LE.1.E-10) THEN
          FA%NSEC2SP(1,KRANGC)=50
        ENDIF
      ELSE
        FA%NSEC2SP(1,KRANGC)=80
        IF (ABS(FA%SCODIL(KRANGC)-1.).LE.1.E-10) THEN
          FA%NSEC2SP(1,KRANGC)=60
        ENDIF
      ENDIF
C Troncature (3 fois la meme si triangulaire)
      FA%NSEC2SP(2,KRANGC) =FA%MTRONC(KRANGC)
      FA%NSEC2SP(3,KRANGC) =FA%MTRONC(KRANGC)
      FA%NSEC2SP(4,KRANGC) =FA%MTRONC(KRANGC)
C Type de representation
      FA%NSEC2SP(5,KRANGC) =1
C Mode de representation (2->complex packing)
      FA%NSEC2SP(6,KRANGC) =2
C Reserves
      FA%NSEC2SP(7:11,KRANGC)=0
C Nb de parametres pour la coord verticale
C On prend ici le cas de la coordonnee hybride
C mais le cas de la coord pression sera aisement
C pris en compte + tard (KSEC2(12)=0).
      FA%NSEC2SP(12,KRANGC)=2*(INIVAU+1)
C Latitude du pole sud de rotation
      FA%NSEC2SP(13,KRANGC)=0
C Longitude du pole sud de rotation
      FA%NSEC2SP(14,KRANGC)=0
C Lat et lon du pole d'etirement
      IF (FA%NTYPTR(KRANGC).GE.2) THEN
        FA%NSEC2SP(15,KRANGC)=ZRAMDE*ASIN(FA%SSLAPO(KRANGC))
        FA%NSEC2SP(16,KRANGC)=ZRAMDE*
     S    (SIGN(ONE,FA%SSLOPO(KRANGC))*ACOS(FA%SCLOPO(KRANGC)))
      ELSE
        FA%NSEC2SP(15,KRANGC)=0
        FA%NSEC2SP(16,KRANGC)=0
      ENDIF
C Reserves
      FA%NSEC2SP(17:22,KRANGC)=0
C**
C     2.  -  KSEC2 POUR LA GRILLE DE GAUSS (ARPEGE)
C-----------------------------------------------------------------------
C
C Type de representation de donnees
C
C FA%SSLAPO=sinus latitude du pole d'interet
C (si=1, pole=pole N et pas de rotation)
C FA%SCODIL=coeff de dilation (si =1, pas de dilatation)
      IF ((1.-FA%SSLAPO(KRANGC)).LE.1.E-10) THEN
        FA%NSEC2GG(1,KRANGC)=24
        IF (ABS(FA%SCODIL(KRANGC)-1.).LE.1.E-10) THEN
          FA%NSEC2GG(1,KRANGC)=4
        ENDIF
      ELSE
        FA%NSEC2GG(1,KRANGC)=34
        IF (ABS(FA%SCODIL(KRANGC)-1.).LE.1.E-10) THEN
          FA%NSEC2GG(1,KRANGC)=14
        ENDIF
      ENDIF
C Nb de pts sur un parallele
      FA%NSEC2GG(2,KRANGC)=FA%NXLOPA(KRANGC)
C Nb de pts sur une longitude
      FA%NSEC2GG(3,KRANGC)=INLAT
      ZLATPRE=ASIN(MAX(-ONE,MIN(ONE,FA%SINLAT(1,KRANGC))))
C Latitude (1/1000 degre) du premier pt de grille
      FA%NSEC2GG(4,KRANGC)=ZRAMDE*ZLATPRE
C Longitude (1/1000 degre) du premier pt de grille
      FA%NSEC2GG(5,KRANGC)=0
C Flag pour la resolution (0->on ne donne pas l'increment)
      FA%NSEC2GG(6,KRANGC)=0
C Latitude (1/1000 degre) du dernier pt de grille
      FA%NSEC2GG(7,KRANGC)=-FA%NSEC2GG(4,KRANGC)
C Longitude (1/1000 degre) du dernier pt de grille.
C (FA%NLOPAR(1,KRANGC)=nb de longitudes sur le 1er parallele)
      FA%NSEC2GG(8,KRANGC)=-360000/FA%NLOPAR(1,KRANGC)
C Increment zonal (1/1000 degre)
C Pas de sens ici.
      FA%NSEC2GG(9,KRANGC)=0
C Nb de paralleles entre le pole et l'equateur
      FA%NSEC2GG(10,KRANGC)=(INLAT+1)/2
C Flag pour le mode de balayage
      FA%NSEC2GG(11,KRANGC)=0
C Nombre de parametres pour la coord. verticale
      FA%NSEC2GG(12,KRANGC)=0
C Latitude du pole sud de rotation
      FA%NSEC2GG(13,KRANGC)=0
C Longitude du pole sud de rotation
      FA%NSEC2GG(14,KRANGC)=0
C Latitude et longitude du pole d'etirement
      IF (FA%NTYPTR(KRANGC).GE.2) THEN
        FA%NSEC2GG(15,KRANGC)=ZRAMDE*ASIN(FA%SSLAPO(KRANGC))
        FA%NSEC2GG(16,KRANGC)=ZRAMDE*
     S     (SIGN(ONE,FA%SSLOPO(KRANGC))*ACOS(FA%SCLOPO(KRANGC)))
      ELSE
        FA%NSEC2GG(15,KRANGC)=0
        FA%NSEC2GG(16,KRANGC)=0
      ENDIF
C Flag:      0 -> grille reguliere,  1 -> grille reduite
      IF (FA%NLOPAR(1,KRANGC)==FA%NLOPAR((1+INLAT)/2,KRANGC)) THEN
        FA%NSEC2GG(17,KRANGC)=0
      ELSE
        FA%NSEC2GG(17,KRANGC)=1
      ENDIF
C Flag:      0 -> Terre ronde     , 64 -> Terre ellipsoide
      FA%NSEC2GG(18,KRANGC)=0
C Flag sur les composantes des vecteurs (0->geographique, 8->grille)
      FA%NSEC2GG(19,KRANGC)=0
C Reserves
      FA%NSEC2GG(20:22,KRANGC)=0
C Pour les grilles reduites, nb de points sur chaque parallele
      FA%NSEC2GG(23:22+(1+INLAT)/2,KRANGC)=
     S                 FA%NLOPAR(1:(1+INLAT)/2,KRANGC)
      FA%NSEC2GG(23-MOD(INLAT,2)+(1+INLAT)/2:22+INLAT,KRANGC)=
     S                               FA%NLOPAR((1+INLAT)/2:1:-1,KRANGC)
      GOTO 600
C**
C     3.  -  KSEC2 POUR LA GRILLE LAT-LON (CAS FULL-POS, ARPEGE OU ALADIN)
C-------------------------------------------------------------------------
C
300   CONTINUE
C TEST POUR NEW EGGX
      IF (FA%SINLAT(1,KRANGC) .GE. 0) THEN
C OLD EGGX
      ZLATPRE=FA%SINLAT(7,KRANGC)
      ZLONPRE=FA%SINLAT(4,KRANGC)
      ZLATDER=FA%SINLAT(5,KRANGC)
      ZLONDER=FA%SINLAT(6,KRANGC)
C Type de representation de donnees
      ELSE
C NEW EGGX
      ZLATPRE=FA%SINLAT(16,KRANGC)
      ZLONPRE=FA%SINLAT(13,KRANGC)
      ZLATDER=FA%SINLAT(14,KRANGC)
      ZLONDER=FA%SINLAT(15,KRANGC)
      ENDIF
C
C Type de representation de donnees
      FA%NSEC2LL(1,KRANGC)=0
C Nb de pts sur un parallele
      FA%NSEC2LL(2,KRANGC)=FA%NXLOPA(KRANGC)
C Nb de pts sur une longitude
      FA%NSEC2LL(3,KRANGC)=INLAT
C Latitude (1/1000 degre) du premier pt de grille
      FA%NSEC2LL(4,KRANGC)=NINT(ZRAMDE*ZLATPRE)
C Longitude (1/1000 degre) du premier pt de grille
      FA%NSEC2LL(5,KRANGC)=NINT(ZRAMDE*ZLONPRE)
      IF (FA%NSEC2LL(5,KRANGC).GT.180000)
     S                      FA%NSEC2LL(5,KRANGC)=FA%NSEC2LL(5,KRANGC)-360000
C Flag pour la resolution (128->on donne l'increment: grille reguliere)
      FA%NSEC2LL(6,KRANGC)=128
C Latitude (1/1000 degre) du dernier pt de grille
      FA%NSEC2LL(7,KRANGC)=NINT(ZRAMDE*ZLATDER)
C Longitude (1/1000 degre) du dernier pt de grille
      FA%NSEC2LL(8,KRANGC)=NINT(ZRAMDE*ZLONDER)
      IF (FA%NSEC2LL(8,KRANGC).GT.180000)
     S                      FA%NSEC2LL(8,KRANGC)=FA%NSEC2LL(8,KRANGC)-360000
C Increment zonal (1/1000 degre)
      IF (ZLONPRE.GT.ZLONDER) THEN
        FA%NSEC2LL(9,KRANGC)=
     S    NINT((ZLONDER+2.*ZPI-ZLONPRE)*ZRAMDE/(FA%NXLOPA(KRANGC)-1))
      ELSE
        FA%NSEC2LL(9,KRANGC)=
     S        NINT((ZLONDER-ZLONPRE)*ZRAMDE/(FA%NXLOPA(KRANGC)-1))
      ENDIF
C Increment meridien (1/1000 degre)
      FA%NSEC2LL(10,KRANGC)=
     S          NINT((ZLATPRE-ZLATDER)*ZRAMDE/(INLAT-1))
C Flag pour le mode de balayage: W->E et S->N = 64; W->E et N->S = 0
C Full-Pos produit des champs lat-lon ranges S->N pour ARP et ALD.
C Or la BDAP attend un rangt N->S pour les grilles lat-lon.
C FA renverse donc les champs issus de Full-Pos avant codage GRIBEX.
C 
      FA%NSEC2LL(11,KRANGC)=0
C Nombre de parametres pour la coord. verticale
      FA%NSEC2LL(12,KRANGC)=0
C Latitude du pole sud de rotation
      FA%NSEC2LL(13,KRANGC)=0
C Longitude du pole sud de rotation
      FA%NSEC2LL(14,KRANGC)=0
C Latitude et longitude du pole d'etirement
      FA%NSEC2LL(15,KRANGC)=0
      FA%NSEC2LL(16,KRANGC)=0
C Flag:      0 -> grille reguliere,  1 -> grille reduite
      FA%NSEC2LL(17,KRANGC)=0
C Flag:      0 -> Terre ronde     , 64 -> Terre ellipsoide
      FA%NSEC2LL(18,KRANGC)=0
C Flag sur les composantes des vecteurs (0->geographique, 8->grille)
      FA%NSEC2LL(19,KRANGC)=0
C Reserves
      FA%NSEC2LL(20:22,KRANGC)=0
C**
C     4.  -  KSEC2 POUR LA GRILLE LAT-LON QUASI-REGULIERE ALADIN
C            (en fait, tenue de camouflage pour les coeff spectraux
C             que l'on va ranger en balayant le 1/4 de l'ellipse
C             verticalement: axe X=axe n (nb d'onde meridien) et
C             axe Y=axe m (nb d'onde zonal) afin de suivre le rangt
C             dans le modele. Seuls les coeff spectraux qui seront
C             compactes sont stockes sur la grille lat-lon, soit
C             tous sauf ceux des axes et ceux inclus dans le carre
C             delimite par la ss-tronc de non-compactage).
C-----------------------------------------------------------------------
C
C Type de representation de donnees
      FA%NSEC2AL(1,KRANGC)=0
C Nb de pts sur un parallele: valeur manquante
      FA%NSEC2AL(2,KRANGC)=2**16 -1
C Nb de pts sur une longitude: nombre d'onde zonal max -1
C associe au nombre d'onde meridien 1 (les CSP sur les axes sont
C extraits des champs de CSP puisque non compactes)
      FA%NSEC2AL(3,KRANGC)=(FA%NOZPAR(6,KRANGC)-
     S                      FA%NOZPAR(5,KRANGC)+1)/4 -1
C Latitude (1/1000 degre) du premier pt de grille: valeur bidon
      FA%NSEC2AL(4,KRANGC)=0
C Longitude (1/1000 degre) du premier pt de grille: valeur bidon
      FA%NSEC2AL(5,KRANGC)=0
C Flag pour la resolution (128->on donne l'increment: grille reguliere)
      FA%NSEC2AL(6,KRANGC)=0
C Latitude (1/1000 degre) du dernier pt de grille: valeur bidon
      FA%NSEC2AL(7,KRANGC)=40000
C Longitude (1/1000 degre) du dernier pt de grille: valeur bidon
      FA%NSEC2AL(8,KRANGC)=40000
C Increment zonal (1/1000 degre)
      FA%NSEC2AL(9,KRANGC)=2**16 -1
C Increment meridien (1/1000 degre): deduit des valeurs bidon
      FA%NSEC2AL(10,KRANGC)=(FA%NSEC2AL(7,KRANGC)-FA%NSEC2AL(4,KRANGC))/
     S                   (FA%NSEC2AL(3,KRANGC)-1)
C Flag pour le mode de balayage
      FA%NSEC2AL(11,KRANGC)=0
C Nombre de parametres pour la coord. verticale
      FA%NSEC2AL(12,KRANGC)=0
C Latitude du pole sud de rotation
      FA%NSEC2AL(13,KRANGC)=0
C Longitude du pole sud de rotation
      FA%NSEC2AL(14,KRANGC)=0
C Latitude et longitude du pole d'etirement
      FA%NSEC2AL(15,KRANGC)=0
      FA%NSEC2AL(16,KRANGC)=0
C Flag:      0 -> grille reguliere,  1 -> grille reduite
      FA%NSEC2AL(17,KRANGC)=1
C Flag:      0 -> Terre ronde     , 64 -> Terre ellipsoide
      FA%NSEC2AL(18,KRANGC)=0
C Flag sur les composantes des vecteurs (0->geographique, 8->grille)
      FA%NSEC2AL(19,KRANGC)=0
C Reserves
      FA%NSEC2AL(20:22,KRANGC)=0
C Les valeurs (22+1:22+FA%MTRONC(KRANGC)-1) representant les nb de pts
C le long de chaque parallele (ici, le nb de coeff spectraux
C pour un meme n (et -n), excepte le triangle et les axes non
C compactes) dependent de la ss-troncature qui depend du fichier
C et ne seront donc pas stockes dans le tableau FA%NSEC2AL qui
C depend du cadre. Le tableau FA%NSC2ALF(FA%JPXTRO-1,FA%JPNXFA) les
C contiendra.
C
C**
C     5.  -  KSEC2 POUR LA GRILLE LAMBERT CONFORME (CAS GENERAL ALADIN)
C-------------------------------------------------------------------------
C
C Type de representation de donnees
      FA%NSEC2LA(1,KRANGC)=3
C Nb de pts sur un parallele
      FA%NSEC2LA(2,KRANGC)=FA%NXLOPA(KRANGC)
C Nb de pts sur une longitude
      FA%NSEC2LA(3,KRANGC)=INLAT
C
C Les parametres communs sont regroupes
C Flag pour la resolution (128->on donne l'increment: grille reguliere)
      FA%NSEC2LA(6,KRANGC)=128
C Reserve
      FA%NSEC2LA(8,KRANGC)=0
C Flag pour le mode de balayage: W->E et S->N = 64; W->E et N->S = 0
      FA%NSEC2LA(11,KRANGC)=64
C Nombre de parametres pour la coord. verticale
      FA%NSEC2LA(12,KRANGC)=0
C Latitude (1/1000 degre) du premier pt de grille
      FA%NSEC2LA(4,KRANGC)=NINT(ZRAMDE*ZLATPRE)
C Longitude (1/1000 degre) du premier pt de grille
      FA%NSEC2LA(5,KRANGC)=NINT(ZRAMDE*ZLONPRE)
      IF (FA%NSEC2LA(5,KRANGC).GT.180000)
     S                      FA%NSEC2LA(5,KRANGC)=FA%NSEC2LA(5,KRANGC)-360000
C
C TEST POUR OLD/NEW EGGX
      IF (FA%SINLAT(1,KRANGC) .GE. 0) THEN
C Old EGGX
C Orientation de la grille
      FA%NSEC2LA(7,KRANGC)=NINT(ZRAMDE*FA%SINLAT(8,KRANGC))
      IF (FA%NSEC2LA(7,KRANGC).GT.180000)
     S                      FA%NSEC2LA(7,KRANGC)=FA%NSEC2LA(7,KRANGC)-360000
C Dimension de la maille dans la direction X
      FA%NSEC2LA(9,KRANGC)=NINT(FA%SINLAT(15,KRANGC))
C Dimension de la maille dans la direction Y
      FA%NSEC2LA(10,KRANGC)=NINT(FA%SINLAT(16,KRANGC))
C Flag pour le centre de projection
C (0: le pole Nord est sur le plan de projection
C  et 1 seul centre de projection est utilise;
C  128: idem sauf que c'est le pole Sud)
      IF (FA%SINLAT(9,KRANGC).GE.0) THEN
        FA%NSEC2LA(13,KRANGC)=0
      ELSE
        FA%NSEC2LA(13,KRANGC)=128
      ENDIF
C Premiere latitude depuis le pole ou le cone coupe la sphere
      FA%NSEC2LA(14,KRANGC)=NINT(ZRAMDE*FA%SINLAT(9,KRANGC))
C Deuxieme latitude depuis le pole ou le cone coupe la sphere
C Dans Aladin, le plan de projection est rarement secant (cela
C releve plus d'un domaine mal defini que d'un choix) et cette
C possibilite va disparaitre bientot. Comme le calcul de cette
C seconde latitude n'est pas aise (pb de convergence), on va
C declarer la grille tangente! mais avec un WARNING...
      FA%NSEC2LA(15,KRANGC)=FA%NSEC2LA(14,KRANGC)
      IF (ABS(FA%SINLAT(10,KRANGC)-SIN(FA%SINLAT(9,KRANGC))).GT.1.E-10 
     S    .AND. FA%LFAMOP) THEN
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S          ' FAISC2: WARNING !! La grille Lambert coupe en fait',
     S          ' la sphere, mais sera consideree comme tangente'
      ENDIF
      ELSE
C NEW EGGX
C Orientation de la grille
      FA%NSEC2LA(7,KRANGC)=NINT(ZRAMDE*FA%SINLAT(3,KRANGC))
      IF (FA%NSEC2LA(7,KRANGC).GT.180000)
     S                      FA%NSEC2LA(7,KRANGC)=FA%NSEC2LA(7,KRANGC)-360000
C Dimension de la maille dans la direction X
      FA%NSEC2LA(9,KRANGC)=NINT(FA%SINLAT(7,KRANGC))
C Dimension de la maille dans la direction Y
      FA%NSEC2LA(10,KRANGC)=NINT(FA%SINLAT(8,KRANGC))
C Flag pour le centre de projection
C (0: le pole Nord est sur le plan de projection
C  et 1 seul centre de projection est utilise;
C  128: idem sauf que c'est le pole Sud)
      IF (FA%SINLAT(4,KRANGC).GE.0) THEN
        FA%NSEC2LA(13,KRANGC)=0
      ELSE
        FA%NSEC2LA(13,KRANGC)=128
      ENDIF
C Premiere latitude depuis le pole ou le cone coupe la sphere
      FA%NSEC2LA(14,KRANGC)=NINT(ZRAMDE*FA%SINLAT(4,KRANGC))
C NEW EGGX toujours tangent
      FA%NSEC2LA(15,KRANGC)=FA%NSEC2LA(14,KRANGC)
      ENDIF

C Reserve
      FA%NSEC2LA(16,KRANGC)=0
C Flag:      0 -> grille reguliere
      FA%NSEC2LA(17,KRANGC)=0
C Flag:      0 -> Terre ronde     , 64 -> Terre ellipsoide
      FA%NSEC2LA(18,KRANGC)=0
C Flag sur les composantes des vecteurs (0->geographique, 8->grille)
      FA%NSEC2LA(19,KRANGC)=8
C Latitude du pole sud
      FA%NSEC2LA(20,KRANGC)=0
C Longitude du pole sud
      FA%NSEC2LA(21,KRANGC)=0
C Reserve
      FA%NSEC2LA(22,KRANGC)=0
C**
C     6.  -  PARTIE REELLE DE LA SECTION 2 DE GRIBEX
C-----------------------------------------------------------------------
C
600   CONTINUE
C Angle de rotation
      FA%XSEC2(1,KRANGC)=0.
C Coefficient d'etirement
      FA%XSEC2(2,KRANGC)=FA%SCODIL(KRANGC)
C Reserve
      FA%XSEC2(3:10,KRANGC)=0.
C Parametres pour la coordonnee verticale
      FA%XSEC2(11:11+INIVAU,KRANGC)=FA%SFOHYB(1,0:INIVAU,KRANGC)*
     S                           FA%SPREFE(KRANGC)
      FA%XSEC2(12+INIVAU:12+2*INIVAU,KRANGC)=
     S     FA%SFOHYB(2,0:INIVAU,KRANGC)
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE SOUS-PROGRAMME "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
C
      IF (FA%LFAMOP.OR.KREP.NE.0) THEN
        INIMES=2
        CLNSPR='FAISC2'
        INUMER=FA%JPNIIL
C
        WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KRANGC='',I4)')
     S     KREP, KRANGC
        CALL FAIPAR_MT (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS,
     S                  CLNSPR,CLNSPR,.FALSE.)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FAISC2_MT',1,ZHOOK_HANDLE)
      END

