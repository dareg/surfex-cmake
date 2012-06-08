C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FADECX_MT (FA,  KREP,   KRANG,  CDNOMA, KVALCO, KLONGA,
     S                    KCHAMP, LDCOSP )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C      Sous-programme INTERNE du logiciel de Fichiers ARPEGE:
C      Controle de coherence et decodage (GRIBEX) d'un CHAMP
C      HORIZONTAL venant d'etre lu sur un fichier ARPEGE/ALADIN.
C       ( DECodage d'un champ gribeX )
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KRANG  (Entree) ==> Rang de l'unite logique;
C                CDNOMA (Entree) ==> Nom d'article (prefabrique);
C    ( Tableau ) KVALCO (Entree) ==> Donnees issues de la lecture;
C                KLONGA (Entree) ==> Nombre de mots lus;
C    ( Tableau ) KCHAMP (Sortie) ==> Valeurs REELLES du champ lu;
C                LDCOSP (Entree) ==> Vrai si le champ est represente
C                                    par des coefficients spectraux;
C*
C  MODIFICATION :
C         JM AUDOIN  15/05/2007  Partie 3.1  Blindage controle changement unite
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KRANG, KLONGA
C
      INTEGER (KIND=JPDBLE) KVALCO(*), KCHAMP(*)
C
      LOGICAL LDCOSP
C
      CHARACTER CDNOMA*(*)
C
      REAL (KIND=JPDBLR) ZSEC2(10+2*(FA%JPXNIV+1)), ZSEC3(2)
      REAL (KIND=JPDBLR), ALLOCATABLE ::  ZSEC4(:), ZCHAMP(:)
      REAL ZPULAP
C
      INTEGER (KIND=JPDBLE), ALLOCATABLE :: ICHAMP(:)
C
      INTEGER ISEC0(2), ISEC1(FA%JPSEC1), ISEC2(FA%JPSEC2), ISEC3(2)
      INTEGER ISEC4(FA%JPSEC4)
      INTEGER ILCHAM, ISTRIA, J, IDECAL, IPOFIN, ILONSEC2
      INTEGER ITRONC, IIND, ILOW, IHIGH
      INTEGER IL, IADD, IRANGC, IILCHAM, IERR, INIMES
      INTEGER IVALC3, IVALC4, IVALC5, IWORD
      INTEGER INUMER, ILENG, IRET, INDEX, JN, JM, JLAT, JLON
      INTEGER IFAORI, IFAMOD, INBIMO
      INTEGER I7,I10,I16
C
      LOGICAL LLMLAM, LLCOSP
C
      CHARACTER*1 CLOPER
C
      INTEGER DECF10
      EXTERNAL DECF10
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  CONTROLES ET INITIALISATIONS.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FADECX_MT',0,ZHOOK_HANDLE)
      KREP=0
      IF (KRANG.LE.0.OR.KRANG.GT.FA%JPNXFA) THEN
        KREP=-66
        GOTO 1001
      ENDIF
C
      CLOPER='D'
      ISTRIA=0
C**
C     2.  -  CONTROLE DES DONNEES DE L'ARTICLE
C-----------------------------------------------------------------------
C
      IF (KVALCO(1).NE.3.OR.
     S    KVALCO(2).LT.0.OR.KVALCO(2).GT.1.OR.
     S    (KVALCO(2).EQ.1.AND.KVALCO(4).LT.0)) THEN
        KREP=-91
        GOTO 1001
      ELSE
        LLCOSP=KVALCO(2).EQ.1
      ENDIF
C
      IF ((LLCOSP.AND..NOT.LDCOSP).OR.(.NOT.LLCOSP.AND.LDCOSP)) THEN
        KREP=-92
        GOTO 1001
      ENDIF
C
      IRANGC=FA%NUCADR(KRANG)
      LLMLAM=FA%LIMLAM(IRANGC)
      ITRONC=FA%MTRONC(IRANGC)
C
      IF (LDCOSP) THEN
        IF (LLMLAM) THEN
          ILCHAM=FA%NSFLAM(IRANGC)
          ILONSEC2=21+ITRONC
        ELSE    
          ILCHAM=(1+ITRONC)*(2+ITRONC)
          ILONSEC2=22
        ENDIF   
      ELSE
        ILCHAM=FA%NVAPDG(IRANGC)
        IF (LLMLAM) THEN
          ILONSEC2=22
        ELSE
          ILONSEC2=22+FA%NLATIT(IRANGC)
        ENDIF
      ENDIF
C
      ALLOCATE (ICHAMP(ILCHAM), ZCHAMP(ILCHAM), ZSEC4(ILCHAM))
C
C**
C     3.  -  DECODAGE GRIBEX DES DONNEES DE L'ARTICLE
C-----------------------------------------------------------------------
C
      IDECAL=3
      IVALC3=KVALCO(3)
      IF (LDCOSP) THEN
        IDECAL=IDECAL+2
C IVALC4=ss-tronc non compactee
C IVALC5=puissance de laplacien
        IVALC4=KVALCO(4)
        IVALC5=KVALCO(5)
      ENDIF
      IILCHAM=ILCHAM
C
C Pour Aladin, le calcul du nb de coeff spectraux qui ont
C ete compactes est plus complexe (certains ont ete retires
C pour ne pas etre compactes: ss-tronc triangulaire).
C
      IF (LDCOSP.AND.LLMLAM) THEN
        ISTRIA=4*(1+FA%NOZPAR(1,IRANGC)+FA%NOZPAR(2,IRANGC)+
     S            IVALC4*(IVALC4-1)/2)
        IILCHAM=ILCHAM-ISTRIA
      ENDIF
C ILENG=longueur disponible en entiers declares INTEGER dans KVALCO.
      ILENG=(KIND(KVALCO)/KIND(ILENG))*(KLONGA-IDECAL)
C
C     3.1 -  CHANGEMENT D'UNITE DE CERTAINS CHAMPS.
C            Il s'agit de champs dont les valeurs sont comprises
C            entre 0 et 1 dans le modele mais dont l'unite
C            conventionnelle dans le GRIB est le %.
C            Avant l'appel a GRIBEX, il faut leur redonner leurs
C            valeurs d'origine (comprises entre 0 et 1) en ajoutant
C            2 au facteur d'echelle decimal KSEC1(23) via DECF10.
C
      I7 =MIN(7,LEN(TRIM(CDNOMA)))
      I10=MIN(10,LEN(TRIM(CDNOMA)))
      I16=MIN(16,LEN(TRIM(CDNOMA)))
      IF (
     S CDNOMA(1:I10)=='SURFNEBUL.' .OR. 
     S CDNOMA(1:I10)=='SURFALBEDO' .OR.
     S CDNOMA=='SURFPROP.VEGETAT' .OR.
     S CDNOMA=='CLSHUMI.RELATIVE' .OR. CDNOMA=='CLSMAXI.HUMI.REL' .OR.
     S CDNOMA=='CLSMINI.HUMI.REL' .OR.
     S (CDNOMA(1:1)=='P'.AND.CDNOMA(I7:I16)=='HUMI_RELAT').OR.
     S (CDNOMA(1:1)=='H'.AND.CDNOMA(I7:I16)=='HUMI_RELAT')      ) THEN
        IADD = 2
        INBIMO = 32  ! Nombre de BIts par mot (un mot=INTEGER)
        INDEX = DECF10 ( KVALCO(IDECAL+1), ILENG, IADD,
     S                   IFAORI, IFAMOD, INBIMO )
        IF (INDEX==-1) THEN
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S           'FADECX: pas d''entete GRIB au debut !'
          KREP=-128
          GOTO 1001
        ELSEIF (INDEX==-2) THEN
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S           'FADECX: edition du GRIB invalid pour DECF10 !'
          KREP=-128
          GOTO 1001
        ELSEIF (INDEX > 0) THEN
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S           'FADECX: ERREUR dans appel a INXBIT par DECF10 !'
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S           'FADECX: avec code retour de INXBIT = ',INDEX
          KREP=-128
          GOTO 1001
        ELSEIF (INDEX < -2) THEN
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S           'FADECX: code retour inconnu de DECF10 : ',INDEX
          KREP=-126
          GOTO 1001
        ENDIF
      ENDIF
C
C     3.2 -  APPEL A GRIBEX
C
      IWORD=0
      IRET=-1
C      CALL GRIBEX(ISEC0,ISEC1,ISEC2,ZSEC2,ISEC3,ZSEC3,ISEC4,
C     S            KCHAMP,IILCHAM,KVALCO(IDECAL+1),ILENG,IWORD,
C     S            CLOPER,IRET)
      IF (FA%LFAMOP) THEN
        WRITE (UNIT=FA%NULOUT,FMT=*) 
     S         ' FADECX: KLONGA, IDECAL, ILENG, IILCHAM = ',
     S         KLONGA, IDECAL, ILENG, IILCHAM
        WRITE (UNIT=FA%NULOUT,FMT=*) '       * ISEC0 = ',ISEC0
        WRITE (UNIT=FA%NULOUT,FMT=*) '       * ISEC1 = ',ISEC1
        WRITE (UNIT=FA%NULOUT,FMT=*) 
     S                     '       * ILONSEC2 ! ISEC2(1:ILONSEC2) = ',
     S                     ILONSEC2, ' ! ', ISEC2(1:ILONSEC2)
        WRITE (UNIT=FA%NULOUT,FMT=*) '       * ZSEC2(1:2) = ',ZSEC2(1:2)
        IF (ISEC2(12).GT.0) WRITE (UNIT=FA%NULOUT,FMT=*)
     S          '       * ISEC2(12) ! ZSEC2(11:10+ISEC2(12)) = ',
     S                    ISEC2(12), ' ! ', ZSEC2(11:10+ISEC2(12))
        WRITE (UNIT=FA%NULOUT,FMT=*) '       * FA%JPSEC4 ! ISEC4 = ',
     S                               FA%JPSEC4,' ! ',ISEC4
      ENDIF
C*
C     3.1 -  CONTROLES DE COHERENCE
C-----------------------------------------------------------------------
C
      IF (IRET.GT.0) THEN
C Erreur rapportee par GRIBEX
        KREP=-1000-IRET
        WRITE (UNIT=FA%NULOUT,FMT=*) ' FADECX: IRET, KREP = ',IRET, KREP 
        GOTO 1001
      ELSEIF (IRET.LT.0) THEN
C Warning rapporte par GRIBEX
        WRITE (UNIT=FA%NULOUT,FMT=*)
        WRITE (UNIT=FA%NULOUT,FMT=*) 
     S         '!------------------------------------------'
        WRITE (UNIT=FA%NULOUT,FMT=*) 
     S         '!           FADECX:   WARNING !!!         !'
        WRITE (UNIT=FA%NULOUT,FMT=*) 
     S         '!------------------------------------------'
        WRITE (UNIT=FA%NULOUT,FMT=*) ' Code retour de GRIBEX = ',
     S        IRET,' pour le champ: ',CDNOMA
        WRITE (UNIT=FA%NULOUT,FMT=*)
      ENDIF
      IF (ISEC4(1).LT.IILCHAM) THEN
        KREP=-93
        IF (FA%LFAMOP) THEN
          WRITE (UNIT=FA%NULOUT,FMT=*) 
     S         'FADECX: ERREUR !!! Nbre de valeurs decodees = ',
     S            ISEC4(1),' et nbre de valeurs attendues = ',IILCHAM
        ENDIF
        GOTO 1001
      ELSEIF (ISEC4(1).GT.IILCHAM) THEN
        KREP=-94
        IF (FA%LFAMOP) THEN
          WRITE (UNIT=FA%NULOUT,FMT=*) 
     S         'FADECX: ERREUR !!! Nbre de valeurs decodees = ',
     S         ISEC4(1),' et nbre de valeurs attendues = ',IILCHAM
        ENDIF
        IF (LLMOER(KREP,KRANG)) GOTO 1001
      ENDIF
C
      IF (IVALC3.NE.ISEC4(2).AND.FA%LFAMOP) THEN
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S     ' FADECX: WARNING, le nb de bits de codage qui avait',
     S     ' ete demande ( ',IVALC3,' ) est different de celui qui a',
     S          ' ete finalement retenu ( ',ISEC4(2),' ) par GRIBEX.'
        WRITE (UNIT=FA%NULOUT,FMT=*)
     S         ' => Gain de place sans perte de precision'
      ENDIF
C
C  Dans le cas d'un champ spectral ARPEGE
C
      IF (LDCOSP.AND..NOT.LLMLAM.AND.(ISEC4(18).NE.IVALC4
     S    .OR.ISEC4(17).NE.IVALC5)) THEN                  
        IF (FA%LFAMOP) THEN
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S            'Ss-tronc non compactee dans GRIB = ',ISEC4(18),
     S            ' et on attend: ',IVALC4
          WRITE (UNIT=FA%NULOUT,FMT=*)
     S            'Puissance de laplacien dans GRIB = ',ISEC4(17),
     S            ' et on attend: ',IVALC5
        ENDIF
        KREP=-95
        GOTO 1001
      ENDIF
C
C Controle de l'adequation entre le nb de mots lus par LFI et le detail:
C ( enrobage FA + message GRIBEX + eventuelles valeurs non-compactees ).
C
      IWORD=1+(ISEC0(1)-1)/JPDBLE
      IF (FA%LFAMOP) THEN
        WRITE (UNIT=FA%NULOUT,FMT=*) ' FADECX: IWORD = ',IWORD
      ENDIF
      IPOFIN=IDECAL+IWORD
      IF (LDCOSP) THEN
        IF (LLMLAM) THEN
          IPOFIN=IPOFIN+ISTRIA
        ELSE
          IPOFIN=IPOFIN+(1+IVALC4)*(2+IVALC4)
        ENDIF
      ENDIF
C
      IF (KLONGA.LT.IPOFIN) THEN
        KREP=-93
        GOTO 1001
      ELSEIF (KLONGA.GT.IPOFIN) THEN
        KREP=-94
        IF (LLMOER(KREP,KRANG)) GOTO 1001
      ENDIF
C*
C     3.2 -  DEMODULATION DES COEFF. SPEC. ALADIN QUI ONT ETE COMPACTES
C-----------------------------------------------------------------------
C
      IF (LDCOSP.AND.LLMLAM) THEN
C  Transfert des donnees decodees et modulees entieres en nombres reels
C  pour les demoduler. Comme KCHAMP est a profil implicite, on ne peut
C  s'en servir pour la fonction TRANSFER => il faut passer par ICHAMP!
        ICHAMP(1:IILCHAM) = KCHAMP(1:IILCHAM)
        ZSEC4 (1:IILCHAM) = TRANSFER(ICHAMP,ZSEC4,IILCHAM)
        ZCHAMP=0.
        ZPULAP=REAL(IVALC5,JPDBLR) * (-0.001)
        IIND=0
        DO JM=1,FA%NOMPAR(2,IRANGC)
          ILOW=2+2*JM+1
          IADD=4*MAX(IVALC4+1-JM,1)
C
          DO INDEX=FA%NOMPAR(ILOW,IRANGC)+IADD,FA%NOMPAR(ILOW+1,IRANGC)
            IIND=IIND+1        
            JN=(INDEX-FA%NOMPAR(ILOW,IRANGC))/4
            ZCHAMP(INDEX)=ZSEC4(IIND) * ((REAL(JN**2+JM**2))**ZPULAP)
          ENDDO
        ENDDO
C  Transfert des donnees decodees et demodulees reelles en nombres entiers
C  disposes aux bons endroits du tableau definitif.
        ICHAMP(1:ILCHAM) = TRANSFER(ZCHAMP,ICHAMP,ILCHAM)
        KCHAMP(1:ILCHAM) = ICHAMP(1:ILCHAM)
      ENDIF
C*
C     3.3 -  TRANSFERT DES COEFFICIENTS SPECTRAUX NON COMPACTES.
C-----------------------------------------------------------------------
C        (et non fournis par GRIBEX) stockes en fin d'article.
C
      IF (LDCOSP) THEN
        IF (LLMLAM) THEN
          IIND=0
          DO JM=0,FA%NOMPAR(2,IRANGC)
            IL=2+2*JM+1
            ILOW=FA%NOMPAR(IL,IRANGC)
C
            IF (JM.EQ.0) THEN
              IHIGH=FA%NOMPAR(IL+1,IRANGC)
            ELSE
              IHIGH=ILOW+4*(IVALC4+1-JM)-1
              IF (IHIGH.LE.ILOW) IHIGH=ILOW+3
            ENDIF
C
            DO INDEX=ILOW,IHIGH
              IIND=IIND+1
              KCHAMP(INDEX)=KVALCO(IDECAL+IWORD+IIND)
            ENDDO
          ENDDO
        ELSE
C
C Cas ARPEGE
C
          KCHAMP(1:2*(IVALC4+1))=
     S        KVALCO(IDECAL+IWORD+1:IDECAL+IWORD+2*(IVALC4+1))
          IIND=2*(IVALC4+1)-1
          INDEX=2*(ITRONC+1)-1
          DO JM=1,IVALC4
          DO JN=JM,ITRONC
            INDEX=INDEX+2
            IF (JN.LE.IVALC4) THEN
              IIND=IIND+2
              KCHAMP(INDEX) = KVALCO(IDECAL+IWORD+IIND)
              KCHAMP(INDEX+1) = KVALCO(IDECAL+IWORD+IIND+1)
            ENDIF
          ENDDO
          ENDDO
C
        ENDIF
      ENDIF
C*
C     3.4 - Renversement des valeurs en pts de grille des champs
C            lat-lon afin de les ranger Sud-Nord plutot que Nord-Sud
C            (on conserve le rangt W-E consecutif) a l'image du rangt
C            initial effectue par FullPos.
C-----------------------------------------------------------------------
C
      IF ((ISEC2(1)==0.OR.ISEC2(1)==10.OR.ISEC2(1)==20.OR.
     S    ISEC2(1)==30) .AND. .NOT.LDCOSP) THEN
        IF (FA%LFAMOP) THEN
          WRITE (UNIT=FA%NULOUT,FMT=*) 
     S            ' FADECX: Grille LAT-LON issue BDAP -> ',
     S            ' renversement des valeurs pour etre rangees SN'
        ENDIF
        DO JLAT=1,FA%NLATIT(IRANGC)
        DO JLON=1,FA%NXLOPA(IRANGC)
          JN=JLON+FA%NXLOPA(IRANGC)*(JLAT-1)
          INDEX=JLON+FA%NXLOPA(IRANGC)*(FA%NLATIT(IRANGC)-JLAT)
          ICHAMP(INDEX) = KCHAMP(JN)
        ENDDO
        ENDDO
        KCHAMP(1:ILCHAM) = ICHAMP(1:ILCHAM)
      ENDIF
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE SOUS-PROGRAMME "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
      IF (ALLOCATED(ICHAMP)) DEALLOCATE ( ICHAMP, ZCHAMP, ZSEC4 )
      LLFATA=LLMOER (KREP,KRANG)
C
      IF (FA%LFAMOP.OR.LLFATA) THEN
        INIMES=2
        CLNSPR='FADECX'
        INUMER=FA%JPNIIL
C
        WRITE (UNIT=CLMESS,FMT='(''KREP='',I5,'', KRANG='',I4,
     S         '', CDNOMA='''''',A,'''''', KLONGA= '',I8,
     S         '', LDCOSP='',L1)')
     S     KREP, KRANG, CDNOMA, KLONGA, LDCOSP
        CALL FAIPAR_MT (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS,
     S                  CLNSPR,CDNOMA,.FALSE.)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FADECX_MT',1,ZHOOK_HANDLE)
      END

