C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FACOND_MT (FA,  KREP,   KNUMER, CDPREF, KNIVAU, CDSUFF,
     S                    PCHAMP, LDCOSP, CDNOMA, KLNOMA, PVALCO,
     S                    KLONGD )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C      Sous-programme de CODAGE d'un CHAMP HORIZONTAL destine a etre
C      ecrit sur un fichier ARPEGE/ALADIN.
C       ( COdage de (Nouvelles ?) Donnees )
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KNUMER (Entree) ==> Numero de l'unite logique;
C                CDPREF (Entree) ==> Prefixe eventuel du nom d'article;
C                KNIVAU (Entree) ==> Niveau vertical eventuel;
C                CDSUFF (Entree) ==> Suffixe eventuel du nom d'article;
C    ( Tableau ) PCHAMP (Entree) ==> Valeurs REELLES du champ a ecrire;
C                LDCOSP (Entree) ==> Vrai si le champ est represente
C                                    par des coefficients spectraux;
C                CDNOMA (Sortie) ==> Nom de l'article-champ a ecrire;
C                KLNOMA (Sortie) ==> Nombre de caracteres utiles dans
C                                    CDNOMA;
C    ( Tableau ) PVALCO (Sortie) ==> Donnees destinees a l'ecriture;
C                KLONGD (Sortie) ==> Nombre de valeurs (mots de 64 bits
C                                    en principe) a ecrire.
C
C    Remarques:
C
C    - PVALCO est type reel par commodite, et doit avoir une longueur
C      suffisante pour stocker les donnees codees. Le dimensionnement
C      "tous terrains" est (2+ILCHAM), qui permet le cas echeant de
C      stocker un champ a pleine resolution sans codage effectif.
C      (ILCHAM est le nombre de valeurs du champ a ecrire)
C
C    - CDNOMA doit avoir au moins FA%JPXNOM caracteres.
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KNUMER, KNIVAU, KLNOMA, KLONGD
C
      REAL (KIND=JPDBLR) PCHAMP (*), PVALCO (*)
C
      CHARACTER CDPREF*(*), CDSUFF*(*), CDNOMA*(*)
C
      INTEGER IREP, ILPRFU, ILSUFU, ILNOMU, ILONGA, IRANG, INIMES
      INTEGER ILPREF, ILSUFF, ILCDNO, IRANGC
      INTEGER IB1PAR (FA%JPLB1P)
C
      LOGICAL LLVERF, LLRLFI, LDCOSP, LLNOMU, LLNOPA
C
      CHARACTER CLPREF*(FA%JPXNOM), CLSUFF*(FA%JPXSUF)
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  CONTROLES ET INITIALISATIONS.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FACOND_MT',0,ZHOOK_HANDLE)
      LLVERF=.FALSE.
      LLRLFI=.FALSE.
      LLNOMU=.FALSE.
      LLNOPA=.FALSE.
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
C     3.  -  FABRICATION DE L'ARTICLE A ECRIRE SUR LE FICHIER.
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
      IF (LDCOSP) THEN
        IF (FA%NFGRIB(IRANG).EQ.-1.OR.FA%NFGRIB(IRANG).EQ.3) THEN
          FA%NRASVE(IRANG)=FA%NRASVE(IRANG)+1
          IF (FA%NRASVE(IRANG).EQ.1.AND.FA%NRASHO(IRANG).GT.0) THEN
            WRITE(FA%NULOUT,*)
     S      '------------------------------------------------'
            WRITE(FA%NULOUT,*)' FACOND :  WARNING !!!!!           '
            WRITE(FA%NULOUT,*)' Un champ de coeff. spectraux avec'
            WRITE(FA%NULOUT,*)
     S      ' rangement type modele va etre ecrit alors que'
            WRITE(FA%NULOUT,*)
     S      ' d''autres champs spec. ont un rangt different.'
            WRITE(FA%NULOUT,*)
     S      '------------------------------------------------'
          ENDIF
        ELSEIF (FA%NFGRIB(IRANG).GE.0.AND.FA%NFGRIB(IRANG).LE.2) THEN
          FA%NRASHO(IRANG)=FA%NRASHO(IRANG)+1
          IF (FA%NRASHO(IRANG).EQ.1.AND.FA%NRASVE(IRANG).GT.0) THEN
            WRITE(FA%NULOUT,*)
     S      '------------------------------------------------'
            WRITE(FA%NULOUT,*)
     S      ' FACOND :  WARNING !!!!!           '
            WRITE(FA%NULOUT,*)
     S      ' Un champ de coeff. spectraux avec'
            WRITE(FA%NULOUT,*)
     S      ' rangt autre que celui du modele va etre ecrit'
            WRITE(FA%NULOUT,*)
     S      ' alors que d''autres champs ont le rangt modele'
            WRITE(FA%NULOUT,*)
     S      '------------------------------------------------'
          ENDIF
        ENDIF
      ENDIF
C
500   CONTINUE
C
      IF (FA%NFGRIB(IRANG).EQ.3) THEN
C Cas d'un champ qu'il faut "griber" avec GRIBEX
        CALL FACODX_MT (FA,  IREP, IRANG, CDPREF, KNIVAU, CDSUFF, 
     S                  PCHAMP, LDCOSP, PVALCO, KLONGD)
C
C Cas particulier de l'erreur GRIBEX num 710: OUTPUT ARRAY TOO SMALL
C On s'en sert pour detecter un probleme de compactage lie a ce que
C le champ compacte + les descripteurs prennent plus de place que le
C champ non compacte...
C On sort donc du compactage (FACODX) pour demander un codage sans
C compactage (FACINE) avec rangement des valeurs selon le modele:
C FA%NFGRIB=-1.
C
        IF (IREP==-1710) THEN
          IREP = 0
          FA%NFGRIB(IRANG) = -1
          LLNOPA = .TRUE.
          GOTO 500
        ENDIF
      ELSE
        CALL FACINE_MT (FA,  IREP, IRANG, CDNOMA(1:ILNOMU), PCHAMP,
     S                LDCOSP, PVALCO, KLONGD, IB1PAR )
        IF (LLNOPA) FA%NFGRIB(IRANG) = 3
C  Le codage num 3 avait ete demande mais se revelait etre
C  plus gourmand en place que le num -1: on avait donc force
C  l'absence de compactage (-1). On revient maintenant au codage
C  num 3 pour ce cadre IRANG et les eventuels codages suivants.
C
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
        IF (LHOOK) CALL DR_HOOK('FACOND_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='FACOND'
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
      IF (LHOOK) CALL DR_HOOK('FACOND_MT',1,ZHOOK_HANDLE)
      END

