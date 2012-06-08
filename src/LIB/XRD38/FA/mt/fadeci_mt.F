C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FADECI_MT (FA,  KREP,   KRANG,  CDNOMA, KVALCO, KLONGA,
     S                    KCHAMP, LDCOSP )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C      Sous-programme INTERNE du logiciel de Fichiers ARPEGE:
C      Controle de coherence et decodage d'un CHAMP HORIZONTAL
C      venant d'etre lu sur un fichier ARPEGE/ALADIN.
C       ( DECodage Interne d'un champ lu )
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
C       En mode multi-taches, il doit y avoir verrouillage du fichier
C     concerne avant l'appel au sous-programme.
C
C     Modifications
C     -------------
C
C  Avril 2004, D. Paradis, DSI/DEV:
C
C    -Declaration ICHAUX en ALLOCATABLE, KCHAMP en profil implicite (gain mem.)
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
      REAL (KIND=JPDBLR) ZFOHYB (2)
C
      INTEGER ILCHAM, ISTRIA, J, IDECAL, ICPACK, IPUILA, IPOFIN
      INTEGER ITRONC, IIND, ILOW, IHIGH, JTRON, IDIMNC, INBITS
      INTEGER IL, IADD, IRANGC, IILCHAM, INDECO, IERR, INIMES
      INTEGER IVALC3, IVALC4, IVALC5, IJLENV, IJLENF, IDIZAI, IUNITE
      INTEGER INUMER
C
      INTEGER (KIND=JPDBLE), ALLOCATABLE :: ICHAUX(:)
      INTEGER IB1PAR (FA%JPLB1P), IB2PAR (FA%JPLB2P)
C
      LOGICAL LLARPE, LLMLAM, LLCOSP
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  CONTROLES ET INITIALISATIONS.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FADECI_MT',0,ZHOOK_HANDLE)
      IF (KRANG.LE.0.OR.KRANG.GT.FA%JPNXFA) THEN
        KREP=-66
        GOTO 1001
      ENDIF
C
      ISTRIA=0
C**
C     2.  -  CONTROLE DES DONNEES DE L'ARTICLE
C-----------------------------------------------------------------------
C
      IF (KVALCO(1).LT.-1.OR.KVALCO(1).GT.2.OR.
     S    KVALCO(2).LT.0.OR.KVALCO(2).GT.1.OR.
     S    (KVALCO(1).GT.0.AND.KVALCO(2).EQ.1.AND.KVALCO(4).LT.0)) THEN
        KREP=-91
        GOTO 1001
      ELSE
        LLARPE=KVALCO(1).EQ.2
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
C
      IF (LDCOSP) THEN
        IF (LLMLAM) THEN
           ILCHAM=FA%NSFLAM(IRANGC)
        ELSE    
          IF (KVALCO(1).EQ.-1) THEN
            ILCHAM=(1+FA%MTRONC(IRANGC))*(2+FA%MTRONC(IRANGC))
          ELSE
            ILCHAM=(1+FA%MTRONC(IRANGC))**2
          ENDIF
        ENDIF   
      ELSE
        ILCHAM=FA%NVAPDG(IRANGC)
      ENDIF
C
C**
C     3.  -  DECODAGE DES DONNEES DE L'ARTICLE
C-----------------------------------------------------------------------
C
      IF (KVALCO(1).EQ.-1 .OR. KVALCO(1).EQ.0) THEN
C
C          Cas ou il n'y a aucun codage... controle longueur d'article
C
        IF (KLONGA.LT.(ILCHAM+2)) THEN
          KREP=-93
          GOTO 1001
        ELSEIF (KLONGA.GT.(ILCHAM+2)) THEN
          KREP=-94
          IF (LLMOER(KREP,KRANG)) GOTO 1001
        ENDIF
C
C           Transfert du tableau d'entree a la suite des 2 valeurs
C     documentaires stockees en debut d'article.
C
        DO 301 J=1,ILCHAM
        KCHAMP(J)=KVALCO(2+J)
  301   CONTINUE
C
      ELSE
C*
C     3.1 -  DECODAGE GRIB PROPREMENT DIT (STANDARD OU NON).
C-----------------------------------------------------------------------
C
        IDECAL=1+2*KVALCO(1)
        IF (LDCOSP) IDECAL=IDECAL+2
        IVALC3=KVALCO(3)
        IVALC4=KVALCO(4)
        IVALC5=KVALCO(5)
        IF (LDCOSP.AND.LLMLAM) THEN
C
          ALLOCATE (ICHAUX (ILCHAM))
C
          ITRONC=FA%MTRONC(IRANGC)
          ISTRIA=FA%NOZPAR(4,IRANGC)-FA%NOZPAR(3,IRANGC)+1
          DO 310 JTRON=1,ITRONC 
            IADD=4*(IVALC4+1-JTRON)
            IF (IADD.LE.0) IADD=4
            ISTRIA=ISTRIA+IADD
 310      CONTINUE
          IILCHAM=ILCHAM-ISTRIA
          CALL DECOGA(ICHAUX,IILCHAM,INBITS,FA%NBIMAC,IB1PAR,IB2PAR,
     S               ZFOHYB(1),2,KVALCO(IDECAL+1),KLONGA-IDECAL,
     S               INDECO,IJLENV,IJLENF,ICPACK,IPUILA,IERR,
     S               KVALCO(IDECAL-1),KVALCO(IDECAL),LLARPE)
C
C  Controle de l'adequation entre nb de valeurs attendues/lues
C
          IF (IJLENF.LT.IILCHAM) THEN
            KREP=-93
            IF (FA%LFAMOP) THEN
              WRITE (UNIT=FA%NULOUT,FMT=*) 
     S                 'FADECI: erreur !!! Nbre de valeurs decodees = ',
     S                 IJLENF,' et nbre de valeurs attendues = ',IILCHAM
            ENDIF
            GOTO 1001
          ELSEIF (IJLENF.GT.IILCHAM) THEN
            KREP=-94
            IF (FA%LFAMOP) THEN
              WRITE (UNIT=FA%NULOUT,FMT=*) 
     S                 'FADECI: erreur !!! Nbre de valeurs decodees = ',
     S                 IJLENF,' et nbre de valeurs attendues = ',IILCHAM
            ENDIF
            IF (LLMOER(KREP,KRANG)) GOTO 1001
          ENDIF
          IIND=0
          DO 320 JTRON=1,ITRONC      
          ILOW=2+2*JTRON+1
          IADD=4*(IVALC4+1-JTRON)
          IF (IADD.LE.0) IADD=4  
          DO 320 J=FA%NOZPAR(ILOW,IRANGC)+IADD,FA%NOZPAR(ILOW+1,IRANGC)
          IIND=IIND+1        
          KCHAMP(J)=ICHAUX(IIND)
 320      CONTINUE          
C
          IF (ALLOCATED( ICHAUX )) DEALLOCATE ( ICHAUX )
C
        ELSE                         
          CALL DECOGA (KCHAMP,ILCHAM,INBITS,FA%NBIMAC,IB1PAR,IB2PAR,
     S               ZFOHYB(1),2,KVALCO(IDECAL+1),KLONGA-IDECAL,
     S               INDECO,IJLENV,IJLENF,ICPACK,IPUILA,IERR,
     S               KVALCO(IDECAL-1),KVALCO(IDECAL),LLARPE)
C
C  Controle de l'adequation entre nb de valeurs attendues/lues
C
          IF (IJLENF.LT.ILCHAM) THEN
            KREP=-93
            IF (FA%LFAMOP) THEN
              WRITE (UNIT=FA%NULOUT,FMT=*) 
     S                 'FADECI: erreur !!! Nbre de valeurs decodees = ',
     S                 IJLENF,' et nbre de valeurs attendues = ',ILCHAM
            ENDIF
            GOTO 1001
          ELSEIF (IJLENF.GT.ILCHAM) THEN
            KREP=-94
            IF (FA%LFAMOP) THEN
              WRITE (UNIT=FA%NULOUT,FMT=*) 
     S                 'FADECI: erreur !!! Nbre de valeurs decodees = ',
     S                 IJLENF,' et nbre de valeurs attendues = ',ILCHAM
            ENDIF
            IF (LLMOER(KREP,KRANG)) GOTO 1001
          ENDIF
        ENDIF
C
        IF (IERR.EQ.-2) THEN
          KREP=-93
          GOTO 1001
        ELSEIF (IERR.NE.0) THEN
          KREP=-200+IERR
          GOTO 1001
        ELSEIF (IVALC3.NE.INBITS.OR.(LDCOSP.AND.
     S          ((ICPACK.NE.IVALC4.AND..NOT.LLMLAM)
     S   .OR.(.NOT.LLMLAM.AND.IPUILA.NE.IVALC5)))) THEN                 
          KREP=-95
          GOTO 1001
        ELSEIF (IB1PAR(4).GT.64) THEN
C
C     Controle effectue s'il y a un bloc 2 en retour du decodage.
C
          IDIZAI=IB2PAR(1)/10
          IUNITE=IB2PAR(1)-IDIZAI*10
C
          IF ((LDCOSP.AND..NOT.LLMLAM.AND.
     S     (IUNITE.NE.0.OR.IDIZAI.LT.5.OR.IDIZAI.GT.8))
     S    .OR.(.NOT.LDCOSP.AND.
     S         (IUNITE.NE.4.OR.IDIZAI.LT.0.OR.IDIZAI.GT.3))) THEN
            KREP=-95
            GOTO 1001
          ENDIF
C
        ENDIF
        IF (LDCOSP.AND.LLMLAM) THEN
          ICPACK=IVALC4
          IPUILA=IVALC5
        ENDIF          
C
        IF (LDCOSP) THEN
C
          IF (LLARPE) THEN
            IF (LLMLAM) THEN
              IPOFIN=IDECAL+INDECO+ISTRIA
            ELSE   
              IDIMNC=(1+ICPACK)**2
              IPOFIN=IDECAL+INDECO+IDIMNC
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
C     3.2 -  TRANSFERT DES COEFFICIENTS SPECTRAUX NON COMPACTES.
C-----------------------------------------------------------------------
C        (et non fournis par DECOGA) stockes en fin d'article.
C
            IF (LLMLAM) THEN
       	       IIND=0
               DO 3201 JTRON=0,ITRONC
               IL=2+2*JTRON+1
               ILOW=FA%NOZPAR(IL,IRANGC)
               IF (JTRON.EQ.0) THEN
                  IHIGH=FA%NOZPAR(IL+1,IRANGC)
               ELSE
                  IHIGH=ILOW+4*(ICPACK+1-JTRON)-1
                  IF (IHIGH.LE.ILOW) IHIGH=ILOW+3
               ENDIF     
               DO 3201 J=ILOW,IHIGH
               IIND=IIND+1
               KCHAMP(J)=KVALCO(IDECAL+INDECO+IIND)
 3201          CONTINUE         
            ELSE
              DO 321 J=1,IDIMNC
              KCHAMP(J)=KVALCO(IDECAL+INDECO+J)
  321         CONTINUE
            ENDIF
C
          ENDIF
C*
C     3.3 -  SI NECESSAIRE, RECONSTITUTION DU SPECTRE.
C-----------------------------------------------------------------------
C
          IF (IPUILA.NE.0) THEN
            CALL FARCIS_MT (FA, KREP,KRANG,KCHAMP,ICPACK,IPUILA)
            IF (KREP.NE.0) GOTO 1001
          ENDIF
C
        ENDIF
C
      ENDIF
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
        CLNSPR='FADECI'
        INUMER=FA%JPNIIL
C
        WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KRANG='',I4,
     S         '', CDNOMA='''''',A,'''''', LDCOSP= '',L1,
     S         '', KLONGA='',I8)')
     S     KREP, KRANG, CDNOMA, LDCOSP, KLONGA
        CALL FAIPAR_MT (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS,
     S                  CLNSPR,CDNOMA,.FALSE.)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('FADECI_MT',1,ZHOOK_HANDLE)
      END

