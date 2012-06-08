C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FACINE_MT (FA,  KREP,   KRANG,  CDNOMA, KCHAMP, LDCOSP,
     S                      KVALCO, KLONGD, KB1PAR )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C      Sous-programme INTERNE du logiciel de Fichiers ARPEGE:
C      PREPARATION (codage GRIB ou non) d'un CHAMP HORIZONTAL
C      destine a etre ecrit sur un fichier ARPEGE/ALADIN.
C       ( Codage Interne d'un (Nouveau ?) champ a Ecrire )
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KRANG  (Entree) ==> Rang de l'unite logique;
C                CDNOMA (Entree) ==> Nom d'article (prefabrique);
C    ( Tableau ) KCHAMP (Entree) ==> Valeurs REELLES du champ a ecrire;
C                LDCOSP (Entree) ==> Vrai si le champ est represente
C                                    par des coefficients spectraux;
C    ( Tableau ) KVALCO (Sortie) ==> Donnees destinees a l'ecriture;
C                KLONGD (Sortie) ==> Nombre de mots a ecrire;
C    ( Tableau ) KB1PAR (Entree+ ==> Image des parametres de la section
C                        Sortie)     1 de GRIB.
C*
C       En mode multi-taches, il doit y avoir verrouillage du fichier
C     concerne avant l'appel au sous-programme.
C
C     Modifications
C     -------------
C  Juin 2001, D. Paradis, DT/DSI/DEV:
C
C    -retrait du compactage lorsqu'il conduit a un article de longueur
C     superieure a celle obtenue sans le compactage (permet aussi de
C     dimensionner KVALCO a ILCHAM+2 sans risquer un debordement)
C
C  Avril 2004, D. Paradis, DT/DSI/DEV:
C
C    -declaration de FA%ICHAMP et FA%ICHAUX en ALLOCATABLE (gain memoire)
C
C January 2010 Trygve Aspelien & Ryad El Khatib : 
C    - workaround against memory leak on IBM
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KRANG, KLONGD
C
      INTEGER (KIND=JPDBLE) KCHAMP(*), KVALCO(*)
      INTEGER KB1PAR (FA%JPLB1P)
C
      LOGICAL LDCOSP
C
      CHARACTER CDNOMA*(*)
C
      INTEGER ILCHAM, ISTRIA, IVALC1, IVALC2, J, IDECAL, ICPACK, IPUILA
      INTEGER ITRONC, IIND, ILOW, IHIGH, JTRON, IDIMNC, ILDISP, INBITS
      INTEGER IL, IADD, IRANGC, IARR, IILCHAM, INMOCC, IERR, INIMES
      INTEGER INUMER, IAUXIL, ITRONC2, ILONGFA, ILONGSEC, ILONGDATA
      INTEGER ILONGD
C
      INTEGER IB2PAR (FA%JPLB2P)
C
      LOGICAL LLARPE, LLMLAM
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  CONTROLES ET INITIALISATIONS.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FACINE_MT',0,ZHOOK_HANDLE)
      IF (KRANG.LE.0.OR.KRANG.GT.FA%JPNXFA) THEN
        KREP=-66
        GOTO 1001
      ENDIF
C**
C     2.  -  FABRICATION DE L'ARTICLE A ECRIRE SUR LE FICHIER.
C-----------------------------------------------------------------------
C
      IVALC1=FA%NFGRIB(KRANG)
      LLARPE=IVALC1.EQ.2
      IRANGC=FA%NUCADR(KRANG)
      LLMLAM=FA%LIMLAM(IRANGC)
C
C  Initialisation du nombre de valeurs du champ (ILCHAM)
C                 du type de champ (IVALC2): spectral/pdg
C                 du type de representation de donnees (IB2PAR(1))
C                 du nombre de bits utilises pour le compactage
C                 de la longueur (en mots de 64 bits) de l'enrobage de l'article
C                          FA (ILONGFA) + les donnees non compactees (en spectral)
C                 de la longueur (en octets) des sections 0, 1 et 5 du GRIB
C                          devant apparaitre dans l'article FA (ILONGSEC)
C                 de la longueur (en bits) de la section 4 du GRIB (ILONGDATA),
C                          devant etre un multiple de 16 bits.
C
C ILONGSEC= 4 (S0) + 24 (S1) + 4 (S5) pour le GRIB version 0
      ILONGSEC=32
      IF (LDCOSP) THEN
C
        IVALC2=1
        INBITS=FA%NBFCSP(KRANG)
C
        IF (LLMLAM) THEN
          IB2PAR(1)=34
          ILCHAM=FA%NSFLAM(IRANGC)
          IF (IVALC1.GT.0) THEN
C calcul du nombre de coefficients non compactes ISTRIA
            ICPACK=FA%NSTROF(KRANG)
            ITRONC=FA%MTRONC(IRANGC)
            ITRONC2=-FA%NTYPTR(IRANGC)
            ISTRIA=4*(1+ITRONC+ITRONC2+(ICPACK*(ICPACK-1)/2))
C
            ILONGFA=3+2*IVALC1+ISTRIA
C  Les 88 bits correspondent aux 11 octets d'enrobage GRIB V0 de la section 4
            ILONGDATA=(ILCHAM-ISTRIA)*INBITS + 88
          ENDIF
        ELSE
          IF (IVALC1.EQ.-1) THEN
            ILCHAM=(1+FA%MTRONC(IRANGC))*(2+FA%MTRONC(IRANGC))
          ELSE
            ILCHAM=(1+FA%MTRONC(IRANGC))**2
          ENDIF
          IB2PAR(1)=80
          IF (IVALC1.GT.0) THEN
            ICPACK=FA%NSTROF(KRANG)
C calcul du nombre de coefficients non compactes IDIMNC
            IDIMNC=(1+ICPACK)**2
            ILONGFA=3+2*IVALC1+IDIMNC
C  Les 144 bits correspondent aux 18 octets d'enrobage GRIB V0 de la section 4
            ILONGDATA=ILCHAM*INBITS + IDIMNC*(32-INBITS) + 144
          ENDIF
        ENDIF
C
      ELSE
C
        ILCHAM=FA%NVAPDG(IRANGC)
        IVALC2=0
        IB2PAR(1)=34
        INBITS=FA%NBFPDG(KRANG)
        IF (IVALC1.GT.0) THEN
          ILONGFA=1+2*IVALC1
C  Les 88 bits correspondent aux 11 octets d'enrobage GRIB V0 de la section 4
          ILONGDATA=ILCHAM*INBITS + 88
        ENDIF
C
      ENDIF
C
C  Retrait du compactage si celui-ci s'avere conduire
C  a un article plus long qu'en l'absence de compactage:
C
      IF (IVALC1.GT.0) THEN
C Arrondi de ILONGDATA au premier multiple de 16 superieur ou egal
        ILONGDATA=16*(1+(ILONGDATA-1)/16)
C Calcul du nombre de mots (64 bits) de la partie GRIB
        ILONGD=1+(ILONGDATA+8*ILONGSEC-1)/64
C On ajoute l'enrobage FA et les eventuelles donnees non compactees
        ILONGD=ILONGD+ILONGFA
        IF (ILONGD.GT.ILCHAM+2) THEN
          WRITE (FA%NULOUT,*)
     S    '///// FACINE:  the packing of article ',CDNOMA,
     S    ' is not done because it will generate'
          WRITE (FA%NULOUT,*)
     S    '         a size ( ',ILONGD,' words ) bigger than',
     S    ' the one ( ',ILCHAM+2,' words ) obtained without packing.'
          IVALC1=0
          WRITE (FA%NULOUT,*)
        ENDIF
      ENDIF
C
      ISTRIA = 0
C
      IF (IVALC1.EQ.-1.OR.IVALC1.EQ.0) THEN
C
C          Cas ou il n'y a aucun codage...
C     transfert du tableau d'entree a la suite des 2 valeurs
C     documentaires stockees ci-dessus dans KVALCO.
C
        DO 301 J=1,ILCHAM
        KVALCO(2+J)=KCHAMP(J)
  301   CONTINUE
C
        KLONGD=2+ILCHAM
C
      ELSE
C
C          Cas avec codage GRIB (standard ou non).
C
#ifndef RS6K
        IF ( ASSOCIATED(FA%ICHAMP) ) DEALLOCATE(FA%ICHAMP)
        ALLOCATE (FA%ICHAMP (ILCHAM))
#else
        IF (.NOT. ASSOCIATED(FA%ICHAMP)) THEN
          ALLOCATE (FA%ICHAMP (ILCHAM))
        ELSEIF ( ILCHAM .GT. SIZE(FA%ICHAMP) ) THEN
          DEALLOCATE(FA%ICHAMP)
          ALLOCATE (FA%ICHAMP (ILCHAM))
        ENDIF
#endif
C
        IDECAL=1+2*IVALC1
        KB1PAR(1)=98
        KB1PAR(2)=1
        KB1PAR(3)=254
        KB1PAR(4)=0
        KB1PAR(5)=255
        KB1PAR(9)=MOD (FA%MADATE(1,KRANG),100)
C
        DO 302 J=2,FA%JPLDAT
        KB1PAR(8+J)=FA%MADATE(J,KRANG)
  302   CONTINUE
C
        IB2PAR(6)=2
        IPUILA=FA%NPUFLA(KRANG)
        ITRONC=FA%MTRONC(IRANGC)
C
        IF (LDCOSP) THEN
C
C           Champ en coefficients spectraux... traitement particulier,
C     lie a la possibilite de compacter une (pseudo-)puissance de
C     laplacien du champ a la place du champ, de maniere a augmenter
C     la precision du champ en "aplanissant" le spectre.
C
          CALL FACSIM_MT (FA, KREP,KRANG,KCHAMP,FA%ICHAMP,IPUILA,ICPACK)
          IF (FA%LFAMOP) THEN
            print *,'FACINE: puissance Dolby selectionnee ',IPUILA
          ENDIF
          IF (KREP.NE.0) GOTO 1001
          IF (LLMLAM) THEN
C Copy the elements to be compacted from FA%ICHAMP to a work array
C This is that part of the quarter-ellipse which is out of the triangle of no compacting.
C In addition, the axes of ellipse are also excluded because of zero-coefficients
#ifndef RS6K
            ALLOCATE (FA%ICHAUX (ILCHAM))
#else
            IF (.NOT. ASSOCIATED(FA%ICHAUX)) THEN
              ALLOCATE (FA%ICHAUX (ILCHAM))
            ELSEIF ( ILCHAM .GT. SIZE(FA%ICHAUX) ) THEN
              DEALLOCATE(FA%ICHAUX)
              ALLOCATE (FA%ICHAUX (ILCHAM))
            ENDIF
#endif
            IIND=0
C
            DO 3021 JTRON=1,ITRONC
            ILOW=2+2*JTRON+1
            IADD=4* MAX(ICPACK+1-JTRON,1)
C
            DO 3021 J=FA%NOZPAR(ILOW,IRANGC)+IADD,
     S                FA%NOZPAR(ILOW+1,IRANGC)
            IIND=IIND+1
            FA%ICHAUX(IIND)=FA%ICHAMP(J)
 3021       CONTINUE
C Number of elements in sub-triangle+axes:ISTRIA
            ISTRIA=ILCHAM-IIND
            IDIMNC=0
          ELSE
            ISTRIA=IDIMNC
          ENDIF
C
          IDECAL=IDECAL+2
          ILDISP=ILCHAM+2-IDECAL-(IVALC1-1)*ISTRIA
C
          IF (.NOT.LLARPE) THEN
C
C            Recopie des coefficients spectraux "non compactes",
C     qui sont codes en fait sur 32 bits dans le cas standard de GRIB.
C
            DO 303 J=1,IDIMNC
            FA%ICHAMP(J)=KCHAMP(J)
  303       CONTINUE
C
          ENDIF
C
        ELSE
C
C          Transfert du tableau d'entree dans un tableau local
C     de maniere a eviter l'ecrasement du tableau d'entree par "CODEGA".
C
          DO 305 J=1,ILCHAM
          FA%ICHAMP(J)=KCHAMP(J)
  305     CONTINUE
C
          IDIMNC=0
          ILDISP=ILCHAM+2-IDECAL
        ENDIF
C*
C     3.1 -  CODAGE GRIB PROPREMENT DIT.
C-----------------------------------------------------------------------
C
        IARR=0
C
        IF (LDCOSP.AND.LLMLAM) THEN
          IILCHAM=ILCHAM-ISTRIA
          CALL CODEGA(FA%ICHAUX,IILCHAM,INBITS,FA%NBIMAC,KB1PAR,IB2PAR,
     S                FA%SFOHYB(1,0,IRANGC),2,KVALCO(IDECAL+1),ILDISP,
     S                INMOCC,IARR,     0,IPUILA,IERR,
     S                KVALCO(IDECAL-1),KVALCO(IDECAL),LLARPE)
        ELSE
          CALL CODEGA(FA%ICHAMP,ILCHAM,INBITS,FA%NBIMAC,KB1PAR,IB2PAR,
     S                FA%SFOHYB(1,0,IRANGC),2,KVALCO(IDECAL+1),ILDISP,
     S                INMOCC,IARR,ICPACK,IPUILA,IERR,
     S                KVALCO(IDECAL-1),KVALCO(IDECAL),LLARPE)
        ENDIF
#ifndef RS6K
        IF (ASSOCIATED( FA%ICHAMP )) DEALLOCATE ( FA%ICHAMP )
#endif
C
        IF (IERR.NE.0) THEN
          KREP=-200+IERR
          GOTO 1001
        ELSEIF (LDCOSP) THEN
          KVALCO(4)=ICPACK
          KVALCO(5)=IPUILA
C
          IF (LLARPE) THEN
C*
C     3.2 -  TRANSFERT DES COEFFICIENTS SPECTRAUX NON COMPACTES.
C-----------------------------------------------------------------------
C        (et non traites par CODEGA) en fin d'article.
C
           IF (LLMLAM) THEN
C Copy nonpacked part of kchamp (sub-triangle+axes) into ichaux
             IIND=0
C
             DO 3201 JTRON=0,ITRONC
             IL=2+2*JTRON+1
             ILOW=FA%NOZPAR(IL,IRANGC)
C
             IF (JTRON.EQ.0) THEN
               IHIGH=FA%NOZPAR(IL+1,IRANGC)
             ELSE
               IHIGH=ILOW+4*(ICPACK+1-JTRON)-1
               IF (IHIGH.LE.ILOW) IHIGH=ILOW+3
             ENDIF
C
             DO 3201 J=ILOW,IHIGH
             IIND=IIND+1
             FA%ICHAUX(IIND)=KCHAMP(J)
 3201        CONTINUE
C
             DO 3202 J=1,ISTRIA
             KVALCO(IDECAL+INMOCC+J)=FA%ICHAUX(J)
 3202        CONTINUE
C
#ifndef RS6K
             IF (ASSOCIATED( FA%ICHAUX )) DEALLOCATE ( FA%ICHAUX )
#endif
C
            ELSE
C
              DO 321 J=1,IDIMNC
              KVALCO(IDECAL+INMOCC+J)=KCHAMP(J)
  321         CONTINUE
C
            ENDIF
C
          ENDIF
C
        ENDIF
C
        KVALCO(3)=INBITS
C
        IF (LLMLAM) THEN
           KLONGD=IDECAL+INMOCC+ISTRIA
        ELSE
           KLONGD=IDECAL+INMOCC+IDIMNC
        ENDIF
C
      ENDIF
C
      KVALCO(1)=IVALC1
      KVALCO(2)=IVALC2
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
        CLNSPR='FACINE'
        INUMER=FA%JPNIIL
C
        WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KRANG='',I4,
     S         '', CDNOMA='''''',A,'''''', LDCOSP= '',L1,
     S         '', KLONGD='',I8)')
     S     KREP, KRANG, CDNOMA, LDCOSP, KLONGD
        CALL FAIPAR_MT (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS,
     S                  CLNSPR, CDNOMA,.FALSE.)
      ENDIF
C
#ifndef RS6K
      IF (ASSOCIATED( FA%ICHAUX )) DEALLOCATE ( FA%ICHAUX )
      IF (ASSOCIATED( FA%ICHAMP )) DEALLOCATE ( FA%ICHAMP )
#endif
      IF (LHOOK) CALL DR_HOOK('FACINE_MT',1,ZHOOK_HANDLE)
      END

