C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FAINIG_MT (FA,  KREP,   KRANG,  CDPREF, KNIVAU, CDSUFF,
     S                    LDCOSP, KLCHAM, KSEC1, KSEC2, PSEC2, KSEC3,
     S                    PSEC3, KSEC4 )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C      Sous-programme INTERNE du logiciel de Fichiers ARPEGE:
C      INItialisation de l'entete Gribex d'un champ.
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KRANG  (Entree) ==> Rang de l'unite logique;
C                CDPREF (Entree) ==> Prefixe eventuel du nom d'article;
C                KNIVAU (Entree) ==> Niveau vertical eventuel;
C                CDSUFF (Entree) ==> Suffixe eventuel du nom d'article;
C                LDCOSP (Entree) ==> Vrai si le champ est represente
C                                    par des coefficients spectraux;
C                KLCHAM (Entree) ==> Longueur totale du champ;
C    ( Tableau ) KSEC1  (Sortie) ==> Image des parametres de la section 1
C                                    de GRIBEX;
C                KSEC2  (Sortie) ==> Image des parametres de la section 2
C                                    de GRIBEX, partie entiere;
C                PSEC2  (Sortie) ==> Image des parametres de la section 2
C                                    de GRIBEX, partie reelle;
C                KSEC3  (Sortie) ==> Image des parametres de la section 3
C                                    de GRIBEX, partie entiere;
C                PSEC3  (Sortie) ==> Image des parametres de la section 3
C                                    de GRIBEX, partie reelle;
C                KSEC4  (Sortie) ==> Image des parametres de la section 4
C                                    de GRIBEX, partie entiere;
C*
C     Modifications
C     -------------
C        R. El Ouaraini: 03-Oct-06 introduction du new EGGX pour tester ERPK
C        R. El Khatib : 11-Aug-2009 Bugfix for non-square geometries
C
C
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      REAL (KIND=JPDBLR) PSEC3(*), PSEC2(*)
C
      INTEGER KREP, KRANG, KNIVAU, KLCHAM, KSEC1(*), KSEC2(*), KSEC3(*)
      INTEGER KSEC4(*)
C
      CHARACTER CDPREF*(*), CDSUFF*(*)
C
      LOGICAL LDCOSP
C
      INTEGER IRANGC, INIMES, INUMER, INLAT, INIVAU, INBITS
      INTEGER INIPAR(6), ICPACK
C
      LOGICAL LLMLAM
C
      INTRINSIC LEN_TRIM
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     0.  -  CONTROLES ET INITIALISATIONS PREALABLES
C-----------------------------------------------------------------------
C
C  Controle de la bonne initialisation de la date
C
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FAINIG_MT',0,ZHOOK_HANDLE)
      IF (FA%LCREAF(KRANG)) THEN
        KREP=-85
        GOTO 1001
      ENDIF
      ICPACK=FA%NSTROF(KRANG)
      IRANGC=FA%NUCADR(KRANG)
      INLAT=FA%NLATIT(IRANGC)
      INIVAU=FA%NNIVER(IRANGC)
      LLMLAM=FA%LIMLAM(IRANGC)
      IF (LDCOSP) THEN
        INBITS=FA%NBFCSP(KRANG)
      ELSE
        INBITS=FA%NBFPDG(KRANG)
      ENDIF
C
C**
C     1.  -  SECTION 1: the product definition section
C-----------------------------------------------------------------------
C
C Appel a FAISC1 une seule fois pour un fichier: initialisation
C du tableau FA%NSEC1(2:21,KRANG) qui va servir comme base pour KSEC1:
C
      IF (FA%LISEC1(KRANG)) THEN
        CALL FAISC1_MT(FA, KREP,KRANG)
        IF (KREP.NE.0) GOTO 1001
        FA%LISEC1(KRANG)=.FALSE.
      ENDIF
      KSEC1(1:FA%JPSEC1)=0
      KSEC1(2:21)=FA%NSEC1(2:21,KRANG)
C
C  Initialisation de INIPAR (5 elts de KSEC1 (1 et 6:9) et un indicateur
C  de type de champ: 0->RAS; 2->min/max; 4->cumul)
      CALL FAIPAG_MT (FA,  KREP, KRANG, CDPREF, KNIVAU, CDSUFF, INIPAR )
      IF (KREP.NE.0) GOTO 1001
C  Element 1: version number of code table 2
      KSEC1(1) = INIPAR(1)
C  Element 6: parameter indicator
      KSEC1(6) = INIPAR(2)
      IF (INIPAR(2).LT.0.OR.INIPAR(2).GT.254.AND.FA%LFAMOP) THEN
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S         '----------------------------------------------------'
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S         '    FAINIG: warning, parameter indicator not defined'
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S         'for: ',CDPREF,'  ',CDSUFF,'. Set to 255, by default'
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S         '----------------------------------------------------'
        KSEC1(6) = 255
      ENDIF
C  Element 7: type of level indicator
      KSEC1(7) = INIPAR(3)
C  Element 8: height, pressure, etc of level or top of level
      KSEC1(8) = INIPAR(4)
C  Element 9: height, pressure, etc of level or bottom of level
      KSEC1(9) = INIPAR(5)
C Cas de la periode de reference
      IF (INIPAR(6)==2) THEN
C Convention dans FA (depuis fin 2000): l'echeance precedente
C est stockee dans FA%MADATE(10,KRANG).
        KSEC1(17)=KSEC1(16)
        KSEC1(16)=FA%MADATE(10,KRANG)
        KSEC1(18)=2
C Cas du cumul
      ELSEIF (INIPAR(6)==4) THEN
        KSEC1(17)=KSEC1(16)
        KSEC1(16)=FA%MADATE(10,KRANG)
        KSEC1(18)=4
C Nb de produits inclus dans le cumul: valeur bidon de 1
        KSEC1(19)=1
      ENDIF
C**
C     2.  -  SECTION 2: the grid definition section
C-----------------------------------------------------------------------
C
C Appel a FAISC2 une seule fois pour un cadre, pour initialiser
C les tableaux NSEC2xxx et FA%XSEC2.
C
      IF (FA%LISEC2(IRANGC)) THEN
        CALL FAISC2_MT(FA, KREP,IRANGC)
        IF (KREP.NE.0) GOTO 1001
        FA%LISEC2(IRANGC)=.FALSE.
      ENDIF
C
C Appel a FAIS2F une seule fois pour un fichier Aladin,
C pour initialiser le tableau FA%NSC2ALF (sauf redefinition
C de la ss-tronc dans FAGOTE).
C
      IF (LLMLAM.AND.FA%LISC2F(KRANG)) THEN
        CALL FAIS2F_MT(FA, KREP,KRANG)
        IF (KREP.NE.0) GOTO 1001
        FA%LISC2F(KRANG)=.FALSE.
      ENDIF
      KSEC2(1:FA%JPSEC2)=0
      IF (LLMLAM) THEN
        IF (LDCOSP) THEN
C  Le champ spectral que l'on doit coder va etre represente sur une
C  grille lat-lon quasi-reguliere puisque ce type de coeff. spectraux
C  n'est pas pris en compte dans GRIBEX.
          KSEC2(1:22)=FA%NSEC2AL(1:22,IRANGC)
          KSEC2(23:21+FA%NOMPAR(2,IRANGC))=
     S     FA%NSC2ALF(1:FA%NOMPAR(2,IRANGC)-1,KRANG)
        ELSE
        IF (FA%SINLAT(1,IRANGC) .GE. 0) THEN
C Old EGGX
          IF (FA%SINLAT(10,IRANGC).LT.0) THEN
C  Parametre de projection negatif, donc pas de projection
C  La grille de ce cadre est une grille lat-lon reguliere
C  du type Full-Pos (pour champ ARPEGE ou Aladin)
            KSEC2(1:22)=FA%NSEC2LL(1:22,IRANGC)
          ELSE
C  La grille de ce cadre est donc du type Lambert conforme
C  (cas general de la grille Aladin)
            KSEC2(1:22)=FA%NSEC2LA(1:22,IRANGC)
          ENDIF
      ELSE
C New EGGX
          IF (FA%SINLAT(2,IRANGC).LT.0) THEN
            KSEC2(1:22)=FA%NSEC2LL(1:22,IRANGC)
          ELSE
            KSEC2(1:22)=FA%NSEC2LA(1:22,IRANGC)
          ENDIF
        ENDIF
        ENDIF
      ELSE
        IF (LDCOSP) THEN
          KSEC2(1:22)=FA%NSEC2SP(1:22,IRANGC)
        ELSE
          KSEC2(1:22+INLAT)=FA%NSEC2GG(1:22+INLAT,IRANGC)
        ENDIF
      ENDIF
C
C Controle ultime: on regarde le prefixe pour s'assurer de la
C presence ou non d'une coordonnee hybride sur la verticale,
C seul cas qui impose une description dans la section 2 reelle.
C
      IF (CDPREF=='S') THEN
        KSEC2(12)=2*(INIVAU+1)
      ELSE
        KSEC2(12)=0
      ENDIF
C
      PSEC2(1:10+KSEC2(12))=FA%XSEC2(1:10+KSEC2(12),IRANGC)
C**
C     3.  -  SECTION 3: the bitmap section
C            As KSEC1(5)=128, the Section 3 is omitted => dummy values
C-----------------------------------------------------------------------
C
C     3.1  - INTEGER PART
C
C Flag: 0->bitmap included in the GRIB message, 1->not included
      KSEC3(1)=1
C Value used at missing data points in an INTEGER data field
      KSEC3(2)=0
C
C     3.2  - REAL PART
C
C Ignored
      PSEC3(1)=0.
C Value used at missing data points in an REAL data field
      PSEC3(2)=0.
C**
C     4.  -  SECTION 4: the binary data section (integer part only)
C-----------------------------------------------------------------------
C
C 1: Nb of data values in array PSEC4 to be encoded
      KSEC4(1)=KLCHAM
C 2: Nb of bits used for each encoded value
      KSEC4(2)=INBITS
C 3: Type of data (0:grid point; 128:spherical harmonic coeff)
      KSEC4(3)=0
C 4: Type of packing, only for spectral fields
C    but also to allow 2nd-order packing for grid points fields
C    and for Aladin spectral fields (seen as lat-lon grid points
C    by GRIBEX).
C    (0:simple packing; 64:complex packing and 2nd-order packing) 
      IF (FA%NCOGRIF(2,KRANG)==0) THEN
C If no Additional flags, then 2nd-order packing is not asked!
        KSEC4(4)=0
      ELSE
        KSEC4(4)=64
      ENDIF
      IF (LDCOSP.AND..NOT.LLMLAM) THEN
C For spherical harmonics coeff, complex packing is always done
        KSEC4(3)=128
        KSEC4(4)=64
      ENDIF
C 5: Data representation (0:float; 32:integer)
      KSEC4(5)=0
C 6: Additional flags indicator (0:no; 16:yes)
      KSEC4(6)=FA%NCOGRIF(2,KRANG)
      IF (LDCOSP.AND..NOT.LLMLAM) THEN
C For spherical harmonics coeff, additional flags indicator=0
        KSEC4(6)=0
      ENDIF
C 7: Reserved
      KSEC4(7)=FA%NCOGRIF(3,KRANG)
C 8: Nb of values indicator (0:single datum at each grid point; 64:matrix)
      KSEC4(8)=0
C 9: Secondary bitmaps indicator (0:no; 32:yes)
      KSEC4(9)=FA%NCOGRIF(4,KRANG)
C 10: Values width indicator (0:2nd order values have constant width; 16:not)
      KSEC4(10)=FA%NCOGRIF(5,KRANG)
C 11: Nb of bits for 2nd order values when these have constant width
      KSEC4(11)=FA%NCOGRIF(6,KRANG)
      IF (KSEC4(11).EQ.-99) KSEC4(11)=1-INBITS
C 12: General extended 2nd order packing (0:no; 8:yes)
C 13: Boustrophedonic ordering (0:no; 4:yes)
C 14,15: give the order of spatial differencing; if 0,0 then option rejected
      KSEC4(12:15)=FA%NCOGRIF(7:10,KRANG)
C 16: For complex packing, a pointer to the start of packed data values (octet nb)
      KSEC4(16)=0
C 17: For complex packing, the scaling factor factor P, stored as the INTEGER
C     value P*1000 (in the range -10000,+10000): defined later
      KSEC4(17)=0
C 18: For complex packing, the pentagonal resolution parameter J specifying
C     the truncation of the subset of the data not packed (32 bits)
      KSEC4(18)=0
C 19-20: Idem 18 for resolution parameters K and M
      KSEC4(19)=0
      KSEC4(20)=0
      IF (LDCOSP.AND..NOT.LLMLAM) THEN
C For spherical harmonics coeff (ARPEGE) only
        KSEC4(18)=ICPACK
        KSEC4(19)=ICPACK
        KSEC4(20)=ICPACK
      ENDIF
C 21-33: Reserved
C 34-42: 'X' decoding option
      KSEC4(21:FA%JPSEC4)=0
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE SOUS-PROGRAMME "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
      LLFATA=LLMOER (KREP,KRANG)
C
      IF (FA%LFAMOP.OR.LLFATA) THEN
        INIMES=2
        CLNSPR='FAINIG'
        INUMER=FA%JPNIIL
C
        WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KRANG='',I4,
     S       '', CDPREF='''''',A,'''''', KNIVAU='',I6,
     S       '', CDSUFF='''''',A,'''''', LDCOSP= '',L1)')
     S         KREP,KRANG,CDPREF(1:LEN_TRIM(CDPREF)),KNIVAU,
     S         CDSUFF(1:LEN_TRIM(CDSUFF)),LDCOSP
        CALL FAIPAR_MT (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS,
     S                  CLNSPR,
     S                  CLNSPR,.FALSE.)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FAINIG_MT',1,ZHOOK_HANDLE)
      END

