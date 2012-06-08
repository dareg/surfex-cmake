C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FACODX_MT (FA,  KREP,   KRANG,  CDPREF, KNIVAU, CDSUFF,
     S                    PSEC4, LDCOSP, KVALCO, KLONGD )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C      Sous-programme INTERNE du logiciel de Fichiers ARPEGE:
C      PREPARATION (codage GRIBEX) d'un CHAMP HORIZONTAL
C      destine a etre ecrit sur un fichier ARPEGE/ALADIN.
C       ( CODage d'un champ a l'aide de gribeX )
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KRANG  (Entree) ==> Rang de l'unite logique;
C                CDPREF (Entree) ==> Prefixe eventuel du nom d'article;
C                KNIVAU (Entree) ==> Niveau vertical eventuel;
C                CDSUFF (Entree) ==> Suffixe eventuel du nom d'article;
C    ( Tableau ) PSEC4  (Entree) ==> Valeurs REELLES du champ a ecrire;
C                LDCOSP (Entree) ==> Vrai si le champ est represente
C                                    par des coefficients spectraux;
C    ( Tableau ) KVALCO (Sortie) ==> Donnees destinees a l'ecriture;
C                KLONGD (Sortie) ==> Nombre de mots a ecrire;
C*
C       En mode multi-taches, il doit y avoir verrouillage du fichier
C     concerne avant l'appel au sous-programme.
C
C      Modifications
C      -------------
C         R. El Ouaraini : 03-Oct-06, introduire la nouvelle geometrie pour tester ERPK
C
C         JM AUDOIN  :  15 mai 2007 partie 5 changement unite 
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KRANG, KNIVAU, KLONGD
C
      INTEGER (KIND=JPDBLE) KVALCO(*)
      REAL (KIND=JPDBLR) PSEC4(*)
C
      LOGICAL LDCOSP
C
      CHARACTER CDPREF*(*), CDSUFF*(*)
C
      REAL (KIND=JPDBLR), ALLOCATABLE :: ZSEC4(:)
      INTEGER (KIND=JPDBLE), ALLOCATABLE :: IVALCO(:)
      REAL (KIND=JPDBLR) :: ZMIN, ZA
      REAL ZSEC2(10+2*(FA%JPXNIV+1)), ZSEC3(2), ZPULAP
C
      INTEGER ISEC0(2), ISEC1(FA%JPSEC1), ISEC2(FA%JPSEC2), ISEC3(2)
      INTEGER ISEC4(FA%JPSEC4), ILONSEC2
      INTEGER ILENG, IWORD, IRET, JM, IPULAP
      INTEGER ILCHAM, JN, IDECAL, ICPACK
      INTEGER ITRONC, ILOW, IHIGH, IDIMNC, INBITS
      INTEGER IL, IADD, IRANGC, IILCHAM, INIMES
      INTEGER INUMER,  INDEX, JLAT, JLON, IDECOPT
      INTEGER IFAORI, IFAMOD, INBIMO
C
      LOGICAL LLMLAM
C
      CHARACTER*1 CLOPER
C
      INTEGER DECF10
      EXTERNAL DECF10
      INTRINSIC LEN_TRIM
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  CONTROLES ET INITIALISATIONS.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FACODX_MT',0,ZHOOK_HANDLE)
      IF (KRANG.LE.0.OR.KRANG.GT.FA%JPNXFA) THEN
        KREP=-66
        GOTO 1001
      ENDIF
      ICPACK=FA%NSTROF(KRANG)
      IRANGC=FA%NUCADR(KRANG)
      LLMLAM=FA%LIMLAM(IRANGC)
      ITRONC=FA%MTRONC(IRANGC)
C
      IF (LLMLAM) THEN
        IF (LDCOSP) THEN
          ILONSEC2=21+FA%NOMPAR(2,IRANGC)
        ELSE
          ILONSEC2=22
        ENDIF
      ELSE
        IF (LDCOSP) THEN
          ILONSEC2=22
        ELSE
          ILONSEC2=22+FA%NLATIT(IRANGC)
        ENDIF
      ENDIF
C
      KVALCO(1)=FA%NFGRIB(KRANG)
      IDECAL=3
      IF (LDCOSP) THEN
        IF (LLMLAM) THEN
          ILCHAM=FA%NSFLAM(IRANGC)
        ELSE
          ILCHAM=(1+ITRONC)*(2+ITRONC)
        ENDIF
        KVALCO(2)=1
        INBITS=FA%NBFCSP(KRANG) 
        IDECAL=IDECAL+2
      ELSE
        ILCHAM=FA%NVAPDG(IRANGC)
        KVALCO(2)=0
        INBITS=FA%NBFPDG(KRANG)
      ENDIF
      KVALCO(3)=INBITS
      IILCHAM = ILCHAM
      IDECOPT = 0
C**
C     2.  -  PREPARATION DU TABLEAU DE DONNEES A ECRIRE SUR LE FICHIER.
C-----------------------------------------------------------------------
C
      ALLOCATE (ZSEC4 (ILCHAM))
C
      IF (LDCOSP .AND. LLMLAM) THEN
C
C       Champ ALADIN en coefficients spectraux... traitement particulier,
C     car non prevu dans GRIBEX (il y sera considere comme un champ lat-lon)
C     mais on a la possibilite de compacter une (pseudo-)puissance de
C     laplacien du champ a la place du champ, de maniere a augmenter
C     la precision du champ en "aplanissant" le spectre.
C
C     Determination de la puissance de Laplacien (en 1/1000 ieme)
C
        CALL FAPULA_MT (FA,  KREP, KRANG, PSEC4, IPULAP )
        ZPULAP=REAL(IPULAP,JPDBLR)/1000.
C       ZPULAP=0.
C       IPULAP=0
        IF (FA%LFAMOP) THEN
          print *,'FACODX: puissance de laplacien selectionee ',ZPULAP,
     S          ' pour une sous-tronc de ',ICPACK
        ENDIF
        IF (KREP.NE.0) GOTO 1001
C
C Transfert des coeff spectraux devant etre compactes de PSEC4 a ZSEC4
C avec prise en compte du coefficient (n**2+m**2)**zpulap. Les coefficients
C concernes sont ceux inclus dans le quart de l'ellipse, hors axes (coeff
C nuls), et hors du triangle non-compacte (sous-troncature).
        IILCHAM=0
C
        DO JM=1,FA%NOMPAR(2,IRANGC)
          ILOW=2+2*JM+1
          IADD=4* MAX(ICPACK+1-JM,1)
C
          DO INDEX=FA%NOMPAR(ILOW,IRANGC)+IADD,FA%NOMPAR(ILOW+1,IRANGC)
            IILCHAM=IILCHAM+1
            JN=(INDEX-FA%NOMPAR(ILOW,IRANGC))/4
            ZSEC4(IILCHAM)=PSEC4(INDEX) * ((REAL(JN**2+JM**2))**ZPULAP)
          ENDDO
        ENDDO
C Number of elements in sub-triangle+axes:IDIMNC
        IDIMNC=ILCHAM-IILCHAM
C Recherche de l'amplitude et du min du champ
        ZMIN=ZSEC4(1)
        ZA=0.
        ZMIN = MINVAL(ZSEC4(1:IILCHAM))
        ZA   = MAXVAL(ZSEC4(1:IILCHAM)) - ZMIN
C Recherche du facteur decimal optimal pour utiliser
C au mieux les INBITS dans le codage de ce champ
        IF (FA%LFAMOP) THEN
          WRITE (UNIT=FA%NULOUT,FMT=*)'FACODX: traitement du champ: ',
     S          CDPREF,KNIVAU,CDSUFF
        ENDIF
        CALL FACDEC_MT (FA, KREP, ZA, ZMIN, INBITS, IDECOPT)
        IF (KREP.NE.0) THEN
          IDECOPT = 0
          KREP = 0
        ENDIF
      ELSEIF(LDCOSP .AND. .NOT.LLMLAM) THEN
C
C          Transfert du tableau d'entree dans un tableau local
C     de maniere a eviter l'ecrasement du tableau d'entree par "GRIBEX".
C
        ZSEC4(1:IILCHAM) = PSEC4(1:IILCHAM)
        IDIMNC=(1+ICPACK)*(2+ICPACK)
      ELSE
C
C    CHAMPS NON SPECTRAUX: transfert du tableau d'entree dans un
C    tableau local de maniere a eviter son ecrasement par "GRIBEX".
C
C
        IDIMNC=0
C Tester si Nouvelle ou ancienne geometrie Aladin
      IF (FA%SINLAT(1,IRANGC) .GE. 0) THEN
        IF (LLMLAM .AND. FA%SINLAT(10,IRANGC).LT.0) THEN
C  Parametre de projection negatif, donc pas de projection:
C  Il s'agit d'une grille lat-lon reguliere du type Full-Pos
C  (pour champ ARPEGE ou Aladin). Il faut donc renverser
C  le champ afin de ranger Nord-Sud les valeurs plutot que Sud-Nord
C  (on conserve le rangt W-E consecutif).
C  Le but est de satisfaire la BDAP qui attend un rangt NW-->SE.
C
          IF (FA%LFAMOP) THEN
            WRITE (UNIT=FA%NULOUT,FMT=*)
     S              ' FACODX: Grille LAT-LON pour BDAP -> ',
     S              ' renversement des valeurs pour etre rangees NS'
          ENDIF
          DO JLAT=1,FA%NLATIT(IRANGC)
          DO JLON=1,FA%NXLOPA(IRANGC)
            JN=JLON+FA%NXLOPA(IRANGC)*(JLAT-1)
            INDEX=JLON+FA%NXLOPA(IRANGC)*(FA%NLATIT(IRANGC)-JLAT)
            ZSEC4(INDEX) = PSEC4(JN)
          ENDDO
          ENDDO
        ELSE
          ZSEC4(1:IILCHAM) = PSEC4(1:IILCHAM)
        ENDIF
      ELSE
        IF (LLMLAM .AND. FA%SINLAT(2,IRANGC).LT.0) THEN
          IF (FA%LFAMOP) THEN
            WRITE (UNIT=FA%NULOUT,FMT=*)
     S              ' FACODX: Grille LAT-LON pour BDAP -> ',
     S              ' renversement des valeurs pour etre rangees NS'
          ENDIF
          DO JLAT=1,FA%NLATIT(IRANGC)
          DO JLON=1,FA%NXLOPA(IRANGC)
            JN=JLON+FA%NXLOPA(IRANGC)*(JLAT-1)
            INDEX=JLON+FA%NXLOPA(IRANGC)*(FA%NLATIT(IRANGC)-JLAT)
            ZSEC4(INDEX) = PSEC4(JN)
          ENDDO
          ENDDO
        ELSE
          ZSEC4(1:IILCHAM) = PSEC4(1:IILCHAM)
        ENDIF
      ENDIF
C Recherche de l'amplitude et du min du champ
        ZMIN=ZSEC4(1)
        ZA=0.
        ZMIN = MINVAL(ZSEC4(1:IILCHAM))
        ZA   = MAXVAL(ZSEC4(1:IILCHAM)) - ZMIN
C Recherche du facteur decimal optimal pour utiliser
C au mieux les INBITS dans le codage de ce champ
        IF (FA%LFAMOP) THEN 
        WRITE (UNIT=FA%NULOUT,FMT=*)'FACODX: traitement du champ: ',
     S          CDPREF,KNIVAU,CDSUFF
        ENDIF
        CALL FACDEC_MT (FA, KREP, ZA, ZMIN, INBITS, IDECOPT)
        IF (KREP.NE.0) THEN
          IDECOPT = 0
          KREP = 0
        ENDIF
      ENDIF
C*
C     3.  -  INITIALISATION DE L'ENROBAGE GRIB
C-----------------------------------------------------------------------
C
C     3.1 -  Sections 1, 2, 3 et 4 (sf la partie reelle pour 4)
C
      CALL FAINIG_MT (FA,  KREP, KRANG, CDPREF, KNIVAU, CDSUFF, LDCOSP,
     S              IILCHAM, ISEC1, ISEC2, ZSEC2, ISEC3, ZSEC3, ISEC4 )
      IF (KREP.NE.0) THEN
        GOTO 1001
      ENDIF
C Prise en compte du facteur decimal
      ISEC1(23) = IDECOPT
C
C     3.2 -  Definition du type de codage
C
      CLOPER='C'
      IF (FA%NCOGRIF(1,KRANG)==1) CLOPER='K'
C*
C     4.  -  CODAGE GRIB PROPREMENT DIT
C-----------------------------------------------------------------------
C
      IRET=-1
C ILENG=longueur disponible en nb d'"entiers declares INTEGER" dans KVALCO.
C On part de l'hypothese ou le dimensionnement de KVALCO se fait
C dans la routine appelante a ILCHAM+2 (cas de l'absence de compactage).
      ILENG=(KIND(KVALCO)/KIND(ILENG))*(ILCHAM+2-IDECAL)
      IWORD=0
CDP
CDP  TEST AVEC UNE PUISSANCE DE LAPLACIEN IMPOSEE
CDP
CDP   CALL GRSMKP(0)
CDP ISEC4(17) = 2000
CDP
      IF (FA%LFAMOP) THEN
        WRITE (UNIT=FA%NULOUT,FMT=*)' FACODX: CLOPER = ',CLOPER
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S                    ' FACODX: IILCHAM, ILCHAM, IDECAL, ILENG = ',
     S                    IILCHAM, ILCHAM, IDECAL, ILENG
        WRITE (UNIT=FA%NULOUT,FMT=*)'       * ISEC1 = ',ISEC1
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S                    '       * ILONSEC2 ! ISEC2(1:ILONSEC2) = ',
     S                    ILONSEC2,' ! ', ISEC2(1:ILONSEC2)
        WRITE (UNIT=FA%NULOUT,FMT=*) '       * ZSEC2(1:2) = ',ZSEC2(1:2)
        IF (ISEC2(12).GT.0) WRITE (UNIT=FA%NULOUT,FMT=*)
     S          '       * ISEC2(12) ! ZSEC2(11:10+ISEC2(12)) = ',
     S                    ISEC2(12), ' ! ', ZSEC2(11:10+ISEC2(12))
        WRITE (UNIT=FA%NULOUT,FMT=*)'       * FA%JPSEC4 ! ISEC4 = ',
     S                              FA%JPSEC4,' ! ',ISEC4
        WRITE (UNIT=FA%NULOUT,FMT=*)'       * ZSEC4(1:20) = ',
     S                              ZSEC4(1:20)
      ENDIF
C     WARNING GRIBEX ENLEVE 
C      CALL GRSDBG(0)
C      CALL GRSVCK(0)

C      CALL GRIBEX(ISEC0,ISEC1,ISEC2,ZSEC2,ISEC3,ZSEC3,ISEC4,
C     S            ZSEC4,IILCHAM,KVALCO(IDECAL+1),ILENG,IWORD,
C     S            CLOPER,IRET)
C
      IF (IRET.GT.0) THEN
C Erreur rapportee par GRIBEX
        KREP=-1000-IRET
        GOTO 1001
      ELSEIF (IRET.LT.0) THEN
C Warning rapporte par GRIBEX
        WRITE (UNIT=FA%NULOUT,FMT=*)
        WRITE (UNIT=FA%NULOUT,FMT=*) 
     S               '!------------------------------------------'
        WRITE (UNIT=FA%NULOUT,FMT=*) 
     S               '!           FACODX:   WARNING !!!         !'
        WRITE (UNIT=FA%NULOUT,FMT=*) 
     S               '!------------------------------------------'
        WRITE (UNIT=FA%NULOUT,FMT=*) ' Code retour de GRIBEX = ',
     S        IRET,' pour le champ: ',CDPREF,KNIVAU,CDSUFF
        WRITE (UNIT=FA%NULOUT,FMT=*)
      ENDIF
C
C ISEC0(1) = nb d'octets dans le message GRIB
C IWORD    = nb de mots de JBDBLE octets (64 bits) du message GRIB
      IWORD=1+(ISEC0(1)-1)/JPDBLE
      KLONGD=IDECAL+IWORD+IDIMNC
      IF (FA%LFAMOP) THEN
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S         ' FACODX: longueur du GRIB en nb octets et en mots = ',
     S         ISEC0(1), IWORD
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S         ' FACODX: longueur de l''article FA en mots = ',
     S         KLONGD
        IF (ISEC4(4).EQ.64 .AND. ISEC4(3).EQ.128) THEN
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S           ' FACODX: complex packing with P=',ISEC4(17),
     S           ' and sub trunc = ',ISEC4(18)
        ENDIF
      ENDIF
C
C  CAS D'UN DEPASSEMENT DE LA TAILLE MAX DE L'ARTICLE FINAL
C  On ramene ce cas a celui d'un tableau trop petit dans GRIBEX.
C
      IF (KLONGD.GT.ILCHAM+2) THEN
        IF (FA%LFAMOP) THEN
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S            ' FACODX: article FA + long avec compactage',
     S            ' que sans => on le supprime'
        ENDIF
        IRET=710
        KREP=-1000-IRET
        GOTO 1001
      ENDIF
C
C
C*
C     5.  -  CHANGEMENT D'UNITE DE CERTAINS CHAMPS.
C            Il s'agit de champs dont les valeurs sont comprises
C            entre 0 et 1 dans le modele mais dont l'unite
C            conventionnelle dans le GRIB est le %.
C---------------------------------------------------------------
C
      INDEX = 0
      IF (CDPREF=='SURF') THEN
        IF (CDSUFF(1:6)=='NEBUL.'  .OR. CDSUFF(1:6)=='ALBEDO'  .OR.
     S      CDSUFF=='PROP.VEGETAT' ) THEN
          INDEX = 1
        ENDIF
      ELSEIF (CDPREF(1:LEN(TRIM(CDPREF)))=='CLS') THEN
        IF (CDSUFF=='HUMI.RELATIVE' .OR. CDSUFF=='MAXI.HUMI.REL' .OR.
     S      CDSUFF=='MINI.HUMI.REL' ) THEN
          INDEX = 1
        ENDIF
      ELSEIF (CDPREF(1:1)=='P'.OR.CDPREF(1:1)=='H') THEN
C         blindage 
        IF (CDSUFF(1:10)=='HUMI_RELAT') THEN
          INDEX = 1
        ENDIF
      ENDIF
      IF (INDEX==1) THEN
        IADD = -2
        INBIMO = 32  ! Nombre de BIts par mot (un mot=INTEGER)
        INDEX = DECF10 ( KVALCO(IDECAL+1), ILENG, IADD,
     S                   IFAORI, IFAMOD, INBIMO )
        IF (INDEX==-1) THEN
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S           'FACODX: pas d''entete GRIB au debut !'
          KREP=-128
          GOTO 1001
        ELSEIF (INDEX==-2) THEN
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S           'FACODX: edition du GRIB invalid pour DECF10 !'
          KREP=-128
          GOTO 1001
        ELSEIF (INDEX > 0) THEN
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S           'FACODX: ERREUR dans appel a INXBIT par DECF10 !'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S           '        avec code retour de INXBIT = ',INDEX
          KREP=-128
          GOTO 1001
        ELSEIF (INDEX < -2) THEN
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S           'FACODX: code retour inconnu de DECF10 : ',INDEX
          KREP=-126
          GOTO 1001
        ENDIF
      ENDIF
C
C*
C     6.  -  TRANSFERT DES COEFFICIENTS SPECTRAUX NON COMPACTES.
C-----------------------------------------------------------------------
C        (et non traites par GRIBEX) en fin d'article.
C
      IF (LDCOSP) THEN
        KVALCO(4)=ICPACK
        IF (LLMLAM) THEN
          KVALCO(5)=IPULAP
C Copy nonpacked part of PSEC4 (sub-triangle+axes) into KVALCO
          IILCHAM=0
          DO JM=0,FA%NOMPAR(2,IRANGC)
            IL=2+2*JM+1
            ILOW=FA%NOMPAR(IL,IRANGC)
C
            IF (JM.EQ.0) THEN
              IHIGH=FA%NOMPAR(IL+1,IRANGC)
            ELSE
              IHIGH=ILOW+4*(ICPACK+1-JM)-1
              IF (IHIGH.LE.ILOW) IHIGH=ILOW+3
            ENDIF     
C
            DO INDEX=ILOW,IHIGH
              IILCHAM=IILCHAM+1
              ZSEC4(IILCHAM)=PSEC4(INDEX)
            ENDDO
          ENDDO
          IF (IILCHAM.NE.IDIMNC) THEN
            WRITE (UNIT=FA%NULOUT,FMT='(A35,I10,A11,I10)')
     S            'FACODX: incoherence entre IILCHAM= ',IILCHAM,
     S            'et IDIMNC= ',IDIMNC
            KREP=-126
            GOTO 1001
          ENDIF
        ELSE
          KVALCO(5)=ISEC4(17)
C Recuperation des coeff spectraux non compactes sachant que le
C rangement est fait par colonnes de JM=cst juxtaposees
          ZSEC4(1:2*(ICPACK+1))=PSEC4(1:2*(ICPACK+1))
          IILCHAM=2*(ICPACK+1)-1
          INDEX=2*(ITRONC+1)-1
          DO JM=1,ICPACK
          DO JN=JM,ITRONC
            INDEX=INDEX+2
            IF (JN.LE.ICPACK) THEN
              IILCHAM=IILCHAM+2
              ZSEC4(IILCHAM) = PSEC4(INDEX)
              ZSEC4(IILCHAM+1) = PSEC4(INDEX+1)
            ENDIF
          ENDDO
          ENDDO
          IF (IILCHAM+1.NE.IDIMNC) THEN
            WRITE (UNIT=FA%NULOUT,FMT='(A35,I10,A11,I10)')
     S            'FACODX: incoherence entre IILCHAM+1= ',IILCHAM+1,
     S            'et IDIMNC= ',IDIMNC
            KREP=-126
            GOTO 1001
          ENDIF
        ENDIF
C Les IDIMNC coeff spectraux non compactes doivent etre transferes
C sur le tableau d'entiers KVALCO apres le IDECAL+IWORD ieme elt.
C
C       KVALCO(IDECAL+IWORD+1:KLONGD)=TRANSFER(ZSEC4,KVALCO,IDIMNC)
        ALLOCATE (IVALCO(IDIMNC))
        IVALCO(1:IDIMNC)=TRANSFER(ZSEC4,IVALCO,IDIMNC)
        KVALCO(IDECAL+IWORD+1:KLONGD)=IVALCO(1:IDIMNC)
        DEALLOCATE (IVALCO)
      ENDIF
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE SOUS-PROGRAMME "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
      IF (ALLOCATED(ZSEC4)) DEALLOCATE ( ZSEC4 )
C
C Cas particulier de l'erreur GRIBEX num 710: OUTPUT ARRAY TOO SMALL
C On s'en sert pour detecter un probleme de compactage lie a ce que
C le champ compacte+les descripteurs prennent plus de place que le
C champ non compacte...
C On sort donc du compactage (FACODX) pour demander un codage sans
C compactage (FACINE) avec rangement des valeurs selon le modele:
C FA%NFGRIB=-1.
C
      IF (IRET==710) THEN
        CLNSPR='FACODX'
        INIMES=2
        INUMER=FA%JPNIIL
C
        WRITE (UNIT=CLMESS,FMT='(''KREP='',I5,'', KRANG='',I4,
     S         '', CDPREF='''''',A,'''''', KNIVAU='',I6,
     S         '', CDSUFF='''''',A,'''''', LDCOSP= '',L1,
     S         '', KLONGD='',I6)')
     S     KREP, KRANG, CDPREF(1:LEN_TRIM(CDPREF)), KNIVAU,
     S     CDSUFF(1:LEN_TRIM(CDSUFF)), LDCOSP, KLONGD
        CALL FAIPAR_MT (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS,
     S                  CLNSPR,CLACTI,.FALSE.)
        CLMESS=
     S ' CAUTION: this field is not packed or it will occupy more space'
        CALL FAIPAR_MT (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS,
     S                  CLNSPR,CLACTI,.FALSE.)
        IF (LHOOK) CALL DR_HOOK('FACODX_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
C
C
      LLFATA=LLMOER (KREP,KRANG)
C
      IF (FA%LFAMOP.OR.LLFATA) THEN
        INIMES=2
        CLNSPR='FACODX'
        INUMER=FA%JPNIIL
C
        WRITE (UNIT=CLMESS,FMT='(''KREP='',I5,'', KRANG='',I4,
     S         '', CDPREF='''''',A,'''''', KNIVAU='',I6,
     S         '', CDSUFF='''''',A,'''''', LDCOSP= '',L1,
     S         '', KLONGD='',I6)')
     S     KREP, KRANG, CDPREF(1:LEN_TRIM(CDPREF)), KNIVAU,
     S     CDSUFF(1:LEN_TRIM(CDSUFF)), LDCOSP, KLONGD
        CALL FAIPAR_MT (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS,
     S                  CLNSPR,CLACTI,.FALSE.)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FACODX_MT',1,ZHOOK_HANDLE)
      END

