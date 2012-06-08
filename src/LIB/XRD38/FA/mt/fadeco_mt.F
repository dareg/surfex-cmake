C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FADECO_MT (FA,  KREP,   KNUMER, CDPREF, KNIVAU, CDSUFF,
     S                    LDCOSP, CDNOMA, KLNOMA, KVALCO, KLONGD,
     S                    PCHAMP )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C      Sous-programme de controle et de DECODAGE d'un CHAMP HORIZONTAL
C      venant d'etre lu sur un fichier ARPEGE/ALADIN.
C       ( DECOdage de donnees )
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KNUMER (Entree) ==> Numero de l'unite logique;
C                CDPREF (Entree) ==> Prefixe eventuel du nom d'article;
C                KNIVAU (Entree) ==> Niveau vertical eventuel;
C                CDSUFF (Entree) ==> Suffixe eventuel du nom d'article;
C                LDCOSP (Entree) ==> Vrai si le champ est represente
C                                    par des coefficients spectraux;
C                CDNOMA (Sortie) ==> Nom de l'article-champ lu;
C                KLNOMA (Sortie) ==> Nombre de caracteres utiles dans
C                                    CDNOMA;
C    ( Tableau ) KVALCO (Entree) ==> Donnees issues de la lecture;
C                KLONGD (Entree) ==> Nombre de valeurs (mots de 64 bits
C                                    en principe) lues;
C    ( Tableau ) PCHAMP (Sortie) ==> Valeurs REELLES du champ lu.
C
C    Remarques:
C
C    - KVALCO est type entier, et doit avoir une longueur
C      suffisante pour stocker les donnees codees. Le dimensionnement
C      "tous terrains" est (2+ILCHAM), qui permet le cas echeant de
C      stocker un champ a pleine resolution sans codage effectif.
C      (ILCHAM est le nombre de valeurs du champ a decoder)
C
C    - CDNOMA doit avoir au moins FA%JPXNOM caracteres.
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KNUMER, KNIVAU, KLNOMA, KLONGD
C
C
      INTEGER IREP, ILPRFU, ILSUFU, ILNOMU, ILONGA, IRANG, INIMES
      INTEGER ILPREF, ILSUFF, ILCDNO, IRANGC, IVALC1
      INTEGER IB1PAR (FA%JPLB1P)
C
      REAL (KIND=JPDBLR) PCHAMP (*)
      INTEGER (KIND=JPDBLE) KVALCO(*)
C
      LOGICAL LLVERF, LLRLFI, LDCOSP, LLNOMU
C
      CHARACTER CDPREF*(*), CDSUFF*(*), CDNOMA*(*)
      CHARACTER CLPREF*(FA%JPXNOM), CLSUFF*(FA%JPXSUF)
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  CONTROLES ET INITIALISATIONS.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FADECO_MT',0,ZHOOK_HANDLE)
      LLVERF=.FALSE.
      LLRLFI=.FALSE.
      LLNOMU=.FALSE.
      ILPRFU=LEN (CDPREF)
      ILSUFU=LEN (CDSUFF)
      ILCDNO=LEN (CDNOMA)
      KLNOMA=0
      CALL FANUMU_MT (FA, KNUMER,IRANG)
C
      IF (IRANG.EQ.0) THEN
        IREP=-51
        GOTO 1001
      ELSEIF (ILCDNO.LT.FA%JPXNOM) THEN
        IREP=-65
        GOTO 1001
      ELSE
        CDNOMA=' '
      ENDIF
C
C         Verrouillage eventuel du fichier.
C
      IF (FA%LFAMUL) CALL LFIVER_MT (FA%LFI, FA%VRFICH(IRANG),'ON')
      LLVERF=FA%LFAMUL
C
      IF (FA%LCREAF(IRANG)) THEN
        IREP=-85
        GOTO 1001
      ENDIF
C**
C     2.  -  FABRICATION DU NOM D'ARTICLE VIA LE SOUS-PROGRAMME "FANFAR"
C            ( controles de CDPREF, KNIVAU, CDSUFF inclus )
C-----------------------------------------------------------------------
C
      CALL FANFAR_MT (FA, IREP,IRANG,CDPREF,KNIVAU,CDSUFF,CDNOMA,
     S                IB1PAR(6),ILPRFU,ILSUFU,ILNOMU)
      IF (IREP.NE.0) GOTO 1001
      LLNOMU=.TRUE.
      KLNOMA=ILNOMU
C**
C     3.  -  CONTROLE ET DECODAGE DE L'ARTICLE DEJA LU SUR LE FICHIER.
C-----------------------------------------------------------------------
C
C
C  Controle de l'homogeneite du type de rangement des coeff. spectraux
C  parmi les champs lus/ecrits: ces champs compactes avec
C  FA%NIGRIB=-1 ou 3 doivent etre ranges comme dans le modele ("verticalement"
C  soit selon des colonnes JM=cst consecutives) et contrairement si compactes
C  avec FA%NIGRIB= 0,1 ou 2.
C 
      IRANGC=FA%NUCADR(IRANG)
      IVALC1=KVALCO(1)
      IF (LDCOSP) THEN
        IF (IVALC1.EQ.-1.OR.IVALC1.EQ.3) THEN
          FA%NRASVE(IRANG)=FA%NRASVE(IRANG)+1
          IF (FA%NRASVE(IRANG).EQ.1.AND.FA%NRASHO(IRANG).GT.0) THEN
            WRITE(FA%NULOUT,*)
     S      '------------------------------------------------'
            WRITE(FA%NULOUT,*)' FADECO :  WARNING !!!!!           '
            WRITE(FA%NULOUT,*)' Un champ de coeff. spectraux avec'
            WRITE(FA%NULOUT,*)
     S      ' rangement type modele va etre lu alors que'
            WRITE(FA%NULOUT,*)
     S      ' d''autres champs spect. ont un rangt different.'
            WRITE(FA%NULOUT,*)
     S      ' ***  Prenez en compte cette heterogeneite!  ***'
            WRITE(FA%NULOUT,*)
     S      '------------------------------------------------'
          ENDIF
        ELSEIF (IVALC1.GE.0.AND.IVALC1.LE.2) THEN
          FA%NRASHO(IRANG)=FA%NRASHO(IRANG)+1
          IF (FA%NRASHO(IRANG).EQ.1.AND.FA%NRASVE(IRANG).GT.0) THEN
            WRITE(FA%NULOUT,*)
     S      '------------------------------------------------'
            WRITE(FA%NULOUT,*)' FADECO :  WARNING !!!!!           '
            WRITE(FA%NULOUT,*)' Un champ de coeff. spectraux avec'
            WRITE(FA%NULOUT,*)
     S      ' rangement autre que celui du modele va etre lu'
            WRITE(FA%NULOUT,*)
     S      ' alors que d''autres champs ont le rangt modele'
            WRITE(FA%NULOUT,*)
     S      ' ***  Prenez en compte cette heterogeneite!  ***'
            WRITE(FA%NULOUT,*)
     S      '------------------------------------------------'
          ENDIF
        ENDIF
      ENDIF
C
      IF (IVALC1.EQ.3) THEN
C Cas d'un champ gribe avec GRIBEX
        CALL FADECX_MT (FA,  IREP,   IRANG,  CDNOMA(1:ILNOMU), KVALCO, 
     S               KLONGD, PCHAMP, LDCOSP )
      ELSE
        CALL FADECI_MT (FA,  IREP,   IRANG,  CDNOMA(1:ILNOMU), KVALCO, 
     S               KLONGD, PCHAMP, LDCOSP )
      ENDIF
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE SOUS-PROGRAMME "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
      KREP=IREP
      LLFATA=LLMOER (IREP,IRANG)
C
C        Deverrouillage eventuel du fichier.
C
      IF (LLVERF) CALL LFIVER_MT (FA%LFI, FA%VRFICH(IRANG),'OFF')
C
      IF (LLFATA) THEN
        INIMES=2
      ELSE
        INIMES=IXNVMS(IRANG)
      ENDIF
C
      IF (.NOT.LLFATA.AND.INIMES.NE.2)  THEN 
        IF (LHOOK) CALL DR_HOOK('FADECO_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='FADECO'
C
      IF (ILPRFU.GE.1) THEN
        ILPREF=MIN (ILPRFU,LEN (CLPREF))
        CLPREF(1:ILPREF)=CDPREF(1:ILPREF)
      ELSE
        ILPREF=8
        CLPREF(1:ILPREF)=FA%CHAINC(:ILPREF)
      ENDIF
C
      IF (ILSUFU.GE.1) THEN
        ILSUFF=MIN (ILSUFU,LEN (CLSUFF))
        CLSUFF(1:ILSUFF)=CDSUFF(1:ILSUFF)
      ELSE
        ILSUFF=8
        CLSUFF(1:ILSUFF)=FA%CHAINC(:ILSUFF)
      ENDIF
C
      IF (.NOT.LLNOMU) THEN
        ILNOMU=MIN (ILPREF,FA%NCPCAD)
        CDNOMA(1:ILNOMU)=CLPREF(1:ILPREF)
      ENDIF
C
      WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,
     S       '', CDPREF='''''',A,'''''', KNIVAU='',I6,
     S       '', CDSUFF='''''',A,'''''', LDCOSP= '',L1)')
     S   KREP,KNUMER,CLPREF(1:ILPREF),KNIVAU,CLSUFF(1:ILSUFF),LDCOSP
      CALL FAIPAR_MT (FA, KNUMER,INIMES,IREP,LLFATA,CLMESS,
     S                CLNSPR,CDNOMA(1:ILNOMU),LLRLFI)
C
      IF (LHOOK) CALL DR_HOOK('FADECO_MT',1,ZHOOK_HANDLE)
      END

