C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FAIENC_MT (FA,  KREP,   KNUMER, CDPREF, KNIVAU, CDSUFF,
     S                      PCHAMP, LDCOSP )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C      Sous-programme d'ECRITURE d'un CHAMP HORIZONTAL sur un fichier
C     ARPEGE.
C       ( Integration par Ecriture d'un (Nouveau ?) Champ )
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KNUMER (Entree) ==> Numero de l'unite logique;
C                CDPREF (Entree) ==> Prefixe eventuel du nom d'article;
C                KNIVAU (Entree) ==> Niveau vertical eventuel;
C                CDSUFF (Entree) ==> Suffixe eventuel du nom d'article;
C    ( Tableau ) PCHAMP (Entree) ==> Valeurs REELLES du champ a ecrire;
C                LDCOSP (Entree) ==> Vrai si le champ est represente
C                                    par des coefficients spectraux.
C
C     Modifications
C     -------------
C
C    Avril 1998: Partie "codage" (paragraphe 3 du sous-programme)
C                demenagee dans un sous-programme a usage interne au
C                logiciel (FACINE). Le but est de pouvoir, sur machine
C                a memoire distribuee, separer codage (via FACOND) et
C                ecriture (via FAISAN) afin de paralleliser le codage.
C
C  Avril 2004, D. Paradis, DSI/DEV:
C
C    -Declaration IVALCO en ALLOCATABLE (gain memoire)
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KNUMER, KNIVAU
C
      REAL (KIND=JPDBLR) PCHAMP (*)
C
      CHARACTER CDPREF*(*), CDSUFF*(*)
C
      INTEGER IREP, ILPRFU, ILSUFU, ILNOMU, ILONGA, IRANG, INIMES
      INTEGER ILPREF, ILSUFF
C
      INTEGER (KIND=JPDBLE), ALLOCATABLE :: IVALCO(:)
      INTEGER IB1PAR (FA%JPLB1P)
C
      INTEGER IVALC1, IRANGC, ILCHAM
C
      LOGICAL LLVERF, LLRLFI, LDCOSP, LLNOMU, LLMLAM, LLNOPA
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
      IF (LHOOK) CALL DR_HOOK('FAIENC_MT',0,ZHOOK_HANDLE)
      LLVERF=.FALSE.
      LLRLFI=.FALSE.
      LLNOMU=.FALSE.
      LLNOPA=.FALSE.
      ILPRFU=LEN (CDPREF)
      ILSUFU=LEN (CDSUFF)
      CALL FANUMU_MT (FA, KNUMER,IRANG)
C
      IF (IRANG.EQ.0) THEN
        IREP=-51
        GOTO 1001
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
      CALL FANFAR_MT (FA, IREP,IRANG,CDPREF,KNIVAU,CDSUFF,CLNOMA,
     S             IB1PAR(6), ILPRFU,ILSUFU,ILNOMU)
      IF (IREP.NE.0) GOTO 1001
      LLNOMU=.TRUE.
C**
C     3.  -  CALCUL D'UN MAJORANT POUR LA LONGUEUR DE L'ARTICLE (mots)
C            ( on va prendre le nombre de valeurs du champ +2 :
C              l'absence de compactage est un majorant et les 2 mots
C              correspondent a l'enrobage FA dans ce cas )
C-----------------------------------------------------------------------
C
      IVALC1=FA%NFGRIB(IRANG)
      IRANGC=FA%NUCADR(IRANG)
      LLMLAM=FA%LIMLAM(IRANGC)
      IF (LDCOSP) THEN
        IF (LLMLAM) THEN
          ILCHAM=FA%NSFLAM(IRANGC)
        ELSE
          IF (IVALC1.EQ.-1 .OR. IVALC1.EQ.3) THEN
            ILCHAM=(1+FA%MTRONC(IRANGC))*(2+FA%MTRONC(IRANGC))
          ELSE
            ILCHAM=(1+FA%MTRONC(IRANGC))**2
          ENDIF
        ENDIF
      ELSE
        ILCHAM=FA%NVAPDG(IRANGC)
      ENDIF
C
      ALLOCATE (IVALCO (ILCHAM+2))

C**
C     4.  -  FABRICATION DE L'ARTICLE A ECRIRE SUR LE FICHIER.
C-----------------------------------------------------------------------
C
C  Controle de l'homogeneite du type de rangement de coeff. spectraux
C  parmi les champs lus/ecrits: ces champs compactes avec
C  FA%NIGRIB=-1 ou 3 doivent etre ranges comme dans le modele ("verticalement"
C  soit selon des colonnes JM=cst consecutives) et contrairement si compactes
C  avec FA%NIGRIB= 0,1 ou 2.
C
      IRANGC=FA%NUCADR(IRANG)
      IF (LDCOSP) THEN
        IF (FA%NFGRIB(IRANG).EQ.-1 .OR. FA%NFGRIB(IRANG).EQ.3) THEN
          FA%NRASVE(IRANG)=FA%NRASVE(IRANG)+1
          IF (FA%NRASVE(IRANG).EQ.1 .AND. FA%NRASHO(IRANG).GT.0) THEN
            WRITE(FA%NULOUT,*)
     S      '------------------------------------------------'
            WRITE(FA%NULOUT,*)' FAIENC :  WARNING !!!!!           '
            WRITE(FA%NULOUT,*)' Un champ de coeff. spectraux avec'
            WRITE(FA%NULOUT,*)
     S      ' rangement type modele va etre ecrit alors que'
            WRITE(FA%NULOUT,*)
     S      ' les autres champs ont un rangement different.'
            WRITE(FA%NULOUT,*)
     S      '------------------------------------------------'
          ENDIF
        ELSEIF (FA%NFGRIB(IRANG).GE.0 .AND. FA%NFGRIB(IRANG).LE.2) THEN
          FA%NRASHO(IRANG)=FA%NRASHO(IRANG)+1
          IF (FA%NRASHO(IRANG).EQ.1 .AND. FA%NRASVE(IRANG).GT.0) THEN
            WRITE(FA%NULOUT,*)
     S      '------------------------------------------------'
            WRITE(FA%NULOUT,*)' FAIENC :  WARNING !!!!!           '
            WRITE(FA%NULOUT,*)' Un champ de coeff. spectraux avec'
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
     S               PCHAMP(1), LDCOSP, IVALCO, ILONGA)
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
        CALL FACINE_MT (FA,  IREP, IRANG, CLNOMA(1:ILNOMU), PCHAMP,
     S                LDCOSP, IVALCO, ILONGA, IB1PAR )
        IF (LLNOPA) FA%NFGRIB(IRANG) = 3
C  Le codage num 3 avait ete demande mais se revelait etre
C  plus gourmand en place que le num -1: on avait donc force
C  l'absence de compactage (-1). On revient maintenant au codage
C  num 3 pour ce cadre IRANG et les eventuels codages suivants.
C
      ENDIF
      IF (IREP.NE.0) GOTO 1001
C**
C     5.  -  ECRITURE DE L'ARTICLE "CHAMP" SUR LE FICHIER.
C-----------------------------------------------------------------------
C
      CALL LFIECR_MT (FA%LFI, IREP,KNUMER,CLNOMA(1:ILNOMU),
     S             IVALCO,ILONGA)
      LLRLFI=IREP.NE.0
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE SOUS-PROGRAMME "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
      IF (ALLOCATED( IVALCO )) DEALLOCATE ( IVALCO )
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
        IF (LHOOK) CALL DR_HOOK('FAIENC_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='FAIENC'
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
        CLNOMA(1:ILNOMU)=CLPREF(1:ILPREF)
      ENDIF
C
      WRITE (UNIT=CLMESS,FMT='(''KREP='',I5,'', KNUMER='',I3,
     S       '', CDPREF='''''',A,'''''', KNIVAU='',I6,
     S       '', CDSUFF='''''',A,'''''', LDCOSP= '',L1)')
     S   KREP,KNUMER,CLPREF(1:ILPREF),KNIVAU,CLSUFF(1:ILSUFF),LDCOSP
      CALL FAIPAR_MT (FA, KNUMER,INIMES,IREP,LLFATA,CLMESS,
     S                CLNSPR, CLNOMA(1:ILNOMU),LLRLFI)
C
      IF (LHOOK) CALL DR_HOOK('FAIENC_MT',1,ZHOOK_HANDLE)
      END

