C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FACOCH_MT (FA,  KREP, KNUME1, KNUME2, 
     S                      CDPREF, KNIVAU, CDSUFF )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C      Sous-programme de reCOpie d'un Champ Horizontal d'un fichier
C     ARPEGE sur un autre.
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KNUME1 (Entree) ==> Numero d'unite logique en entree;
C                KNUME2 (Entree) ==> Numero d'unite logique en sortie;
C                CDPREF (Entree) ==> Prefixe eventuel du nom d'article;
C                KNIVAU (Entree) ==> Niveau vertical eventuel;
C                CDSUFF (Entree) ==> Suffixe eventuel du nom d'article.
C
C     Modifications
C     -------------
C
C  Avril 2004, D. Paradis, DSI/DEV:
C    -Declaration IVALCO en ALLOCATABLE (gain memoire)
C  Juin  2004, D. Paradis, DSI/DEV:
C    -Prise en compte des codages type -1 et 3
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KNUME1, KNUME2, KNIVAU
C
      INTEGER IREP, ILPRFU, ILSUFU, ILNOMU, ILONGA, IRANC1, IRANC2
      INTEGER INIMES, J, INUMFI, IPOSEX, INPAHE, INPAHEL, JLAT, IZPAHEL
      INTEGER ISPAHEL, JNIV, ILPREF, ILSUFF, INUMRO, IRANG2, IGRIB
C
      INTEGER (KIND=JPDBLE), ALLOCATABLE :: IVALCO(:)
      INTEGER IRANG (2), INUMER (2), IB1PAR (3)
C
      LOGICAL LLVERF (2), LLRLFI, LLCOSP, LLMESS, LLNOMU
      LOGICAL LLMLAM1, LLMLAM2
C
      CHARACTER CDPREF*(*), CDSUFF*(*)
      CHARACTER CLAUXI*(FA%JPXNOM), CLPREF*(FA%JPXNOM), 
     S          CLSUFF*(FA%JPXSUF)
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  CONTROLES ET INITIALISATIONS.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FACOCH_MT',0,ZHOOK_HANDLE)
      LLRLFI=.FALSE.
      LLMESS=.FALSE.
      LLNOMU=.FALSE.
      ILPRFU=LEN (CDPREF)
      ILSUFU=LEN (CDSUFF)
      IRANC1=0
      IRANC2=0
      INUMER(1)=KNUME1
      INUMER(2)=KNUME2
      LLVERF(1)=.FALSE.
      LLVERF(2)=.FALSE.
      IRANG(2)=0
C
      DO 101 J=1,2
      INUMFI=J
      CALL FANUMU_MT (FA, INUMER(J),IRANG(J))
C
      IF (IRANG(J).EQ.0) THEN
        IREP=-51
        GOTO 1001
      ENDIF
C
C         Verrouillage eventuel du fichier.
C
      IF (FA%LFAMUL) CALL LFIVER_MT (FA%LFI, FA%VRFICH(IRANG(J)),'ON')
      LLVERF(J)=FA%LFAMUL
C
      IF (FA%LCREAF(IRANG(J))) THEN
        IREP=-85
        GOTO 1001
      ENDIF
C*
C       FABRICATION DU NOM D'ARTICLE VIA LE SOUS-PROGRAMME "FANFAR"
C            ( controles de CDPREF, KNIVAU, CDSUFF inclus )
C
      CALL FANFAR_MT (FA, IREP,IRANG(J),CDPREF,KNIVAU,
     S             CDSUFF,CLNOMA,IB1PAR,
     S             ILPRFU,ILSUFU,ILNOMU)
      IF (IREP.NE.0) GOTO 1001
  101 CONTINUE
C
      LLNOMU=.TRUE.
C**
C     2.  -  LECTURE DE L'ARTICLE SUR LE FICHIER, CONTROLES.
C-----------------------------------------------------------------------
C
      CALL LFINFO_MT (FA%LFI, IREP,KNUME1,CLNOMA(1:ILNOMU),
     S             ILONGA,IPOSEX)
C
      IF (IREP.NE.0) THEN
        LLRLFI=.TRUE.
        GOTO 1001
      ELSEIF (ILONGA.EQ.0) THEN
        IREP=-89
        GOTO 1001
      ELSEIF (ILONGA.GT.FA%JPXCHA+2) THEN
        IREP=-90
        GOTO 1001
      ENDIF
C
      ALLOCATE (IVALCO (ILONGA))
      CALL LFILEC_MT (FA%LFI, IREP,KNUME1,
     S             CLNOMA(1:ILNOMU),IVALCO,ILONGA)
      LLRLFI=IREP.NE.0
      IF (LLRLFI) GOTO 1001
C
      IF (IVALCO(1).LT.-1.OR.IVALCO(1).GT.3.OR.
     S    IVALCO(2).LT.0 .OR.IVALCO(2).GT.1.OR.
     S    (IVALCO(1).GT.0.AND.IVALCO(2).EQ.1.AND.IVALCO(4).LT.0)) THEN
        IREP=-91
        GOTO 1001
      ELSE
        LLCOSP=IVALCO(2).EQ.1
      ENDIF
C**
C     3.  -  CONTROLE DE COHERENCE ENTRE LES FICHIERS, VIS-A-VIS DU TYPE
C            DE DONNEES LUES (points de grille/coefficients spectraux).
C-----------------------------------------------------------------------
C
      IRANC1=FA%NUCADR(IRANG(1))
      IRANC2=FA%NUCADR(IRANG(2))
      INPAHE=(1+FA%NLATIT(IRANC1))/2
      LLMLAM1=FA%NTYPTR(IRANC1).LE. -1
      LLMLAM2=FA%NTYPTR(IRANC2).LE. -1
C
      IF (IRANC1.NE.IRANC2) THEN
C
C         On a pris ici une optique souple: n'est fatale qu'une erreur
C     vraiment grossiere. Toute autre discordance est signalee par un
C     message global de niveau 1.
C
        IF ( (LLMLAM1.AND..NOT.LLMLAM2).OR.
     S       (LLMLAM2.AND..NOT.LLMLAM1).OR.
     S       (LLCOSP.AND.((.NOT.LLMLAM1.AND..NOT.LLMLAM2.AND.
     S                     FA%MTRONC(IRANC1).NE.FA%MTRONC(IRANC2)) .OR.
     S                    (LLMLAM1.AND.LLMLAM2.AND.
     S                     FA%MTRONC(IRANC1).NE.FA%MTRONC(IRANC2).AND.
     S                     FA%NTYPTR(IRANC1).NE.FA%NTYPTR(IRANC2) ))   
     S                    ).OR.
     S       (.NOT.LLCOSP.AND.(FA%NLATIT(IRANC1).NE.
     S                         FA%NLATIT(IRANC2).OR.
     S                         FA%NVAPDG(IRANC1).NE.FA%NVAPDG(IRANC2)))
     S     ) THEN
          IREP=-112
          GOTO 1001
C
        ELSEIF (.NOT.LLCOSP) THEN
C
          IF (.NOT.LLMLAM1.AND..NOT.LLMLAM2) THEN
             INPAHEL=INPAHE
          ELSE
             INPAHEL=8
          ENDIF
          DO 301 JLAT=1,INPAHEL
          LLMESS=LLMESS.OR.FA%NLOPAR(JLAT,IRANC1).NE.
     S           FA%NLOPAR(JLAT,IRANC2)
  301     CONTINUE
C
          IF (LLMESS) THEN
            IREP=-112
            GOTO 1001
          ENDIF
C
        ENDIF
C
        LLMESS=FA%MTRONC(IRANC1).NE.FA%MTRONC(IRANC2).OR.
     S         FA%NTYPTR(IRANC1).NE.FA%NTYPTR(IRANC2).OR.
     S       (KNIVAU.GT.0.AND.(FA%NNIVER(IRANC1).NE.
     S                         FA%NNIVER(IRANC2)).OR.
     S                        (FA%SPREFE(IRANC1).NE.
     S                         FA%SPREFE(IRANC2))).OR.
     S         FA%NLATIT(IRANC1).NE.FA%NLATIT(IRANC2).OR.
     S         FA%SSLAPO(IRANC1).NE.FA%SSLAPO(IRANC2).OR.
     S         FA%SCLOPO(IRANC1).NE.FA%SCLOPO(IRANC2).OR.
     S         FA%SSLOPO(IRANC1).NE.FA%SSLOPO(IRANC2).OR.
     S         FA%SCODIL(IRANC1).NE.FA%SCODIL(IRANC2)
C
        IF (.NOT.LLMESS) THEN
C
          IF (.NOT.LLMLAM1.AND..NOT.LLMLAM2) THEN
             INPAHEL=INPAHE
             IZPAHEL=INPAHE
             ISPAHEL=INPAHE
          ELSE
             INPAHEL=8
             IZPAHEL=0
             ISPAHEL=18
          ENDIF
          DO 302 JLAT=1,INPAHEL
          LLMESS=FA%NLOPAR(JLAT,IRANC1).NE.FA%NLOPAR(JLAT,IRANC2)
     S           .OR.LLMESS
  302     CONTINUE
          DO 312 JLAT=1,IZPAHEL
          LLMESS=FA%NOZPAR(JLAT,IRANC1).NE.FA%NOZPAR(JLAT,IRANC2)
     S           .OR.LLMESS
  312     CONTINUE
          DO 322 JLAT=1,ISPAHEL
          LLMESS=FA%SINLAT(JLAT,IRANC1).NE.FA%SINLAT(JLAT,IRANC2)
     S           .OR.LLMESS
  322     CONTINUE
C
          IF (.NOT.LLMESS.AND.KNIVAU.GT.0) THEN
C
            DO 303 JNIV=0,FA%NNIVER(IRANC1)
            LLMESS=FA%SFOHYB(1,JNIV,IRANC1).NE.
     S             FA%SFOHYB(1,JNIV,IRANC2).OR.
     S      LLMESS.OR.FA%SFOHYB(2,JNIV,IRANC1).NE.
     S             FA%SFOHYB(2,JNIV,IRANC2)
  303       CONTINUE
C
          ENDIF
C
        ENDIF
C
      ENDIF
C**
C     4.  -  ECRITURE DE L'ARTICLE "CHAMP" SUR LE FICHIER.
C-----------------------------------------------------------------------
C
C        Deverrouillage eventuel de l'unite logique d'entree.
C
      IF (LLVERF(1)) CALL LFIVER_MT (FA%LFI, FA%VRFICH(IRANG(1)),'OFF')
      LLVERF(1)=.FALSE.
C
      CALL LFIECR_MT (FA%LFI, IREP,KNUME2,CLNOMA(1:ILNOMU),
     S             IVALCO,ILONGA)
      INUMFI=2
      LLRLFI=IREP.NE.0
      IF (LLRLFI) GOTO 1001
C
C  Controle de l'homogeneite du type de rangement de coeff. spectraux
C  parmi les champs lus/ecrits: ces champs compactes avec
C  IGRIB=-1 ou 3 doivent etre ranges comme dans le modele ("verticalement"
C  soit selon des colonnes JM=cst consecutives) et contrairement si compactes
C  avec IGRIB= 0,1 ou 2.
C
      IRANG2 = IRANG(2)
      IGRIB = IVALCO(1)
      IF (LLCOSP) THEN
        IF (IGRIB.EQ.-1 .OR. IGRIB.EQ.3) THEN
          FA%NRASVE(IRANG2)=FA%NRASVE(IRANG2)+1
          IF (FA%NRASVE(IRANG2).EQ.1 .AND. FA%NRASHO(IRANG2).GT.0) THEN
            WRITE(FA%NULOUT,*)
     S      '------------------------------------------------'
            WRITE(FA%NULOUT,*)' FACOCH : WARNING !!!!!           '
            WRITE(FA%NULOUT,*)
     S      ' Un champ de coef. spect. avec rangt type modele'
            WRITE(FA%NULOUT,*)' va etre ecrit sur l''unite ',KNUME2,
     S                ' alors que'
            WRITE(FA%NULOUT,*)
     S      ' d''autres champs y ont un rangement different.'
            WRITE(FA%NULOUT,*)
     S      '------------------------------------------------'
          ENDIF
        ELSEIF (IGRIB.GE.0 .AND. IGRIB.LE.2) THEN
          FA%NRASHO(IRANG2)=FA%NRASHO(IRANG2)+1
          IF (FA%NRASHO(IRANG2).EQ.1 .AND. FA%NRASVE(IRANG2).GT.0) THEN
            WRITE(FA%NULOUT,*)
     S      '------------------------------------------------'
            WRITE(FA%NULOUT,*)' FACOCH : WARNING !!!!!           '
            WRITE(FA%NULOUT,*)
     S      ' Un champ de coef. spect. avec rangt autre que'
            WRITE(FA%NULOUT,*)
     S      ' celui du modele va etre ecrit sur l''unite ', KNUME2
            WRITE(FA%NULOUT,*)
     S      ' alors que d''autres champs y ont le rangt modele'
            WRITE(FA%NULOUT,*)
     S      '------------------------------------------------'
          ENDIF
        ENDIF
      ENDIF
C
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE SOUS-PROGRAMME "FAIPAR" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
      IF (ALLOCATED( IVALCO )) DEALLOCATE ( IVALCO )
      KREP=IREP
      LLFATA=LLMOER (IREP,IRANG(INUMFI))
C
C        Deverrouillage eventuel des fichiers.
C
      IF (LLVERF(1)) CALL LFIVER_MT (FA%LFI, FA%VRFICH(IRANG(1)),'OFF')
      IF (LLVERF(2)) CALL LFIVER_MT (FA%LFI, FA%VRFICH(IRANG(2)),'OFF')
C
      CLNSPR='FACOCH'
C
C        Messages d'avertissement eventuels.
C
      IF (FA%NIMSGA.NE.0.AND.IREP.EQ.0) THEN
C
        IF (LLMESS) THEN
          INIMES=1
          WRITE (UNIT=CLMESS,FMT='(''*ATTENTION* - LES UNITES'',I3,
     S   '' ET'',I3,'' ONT DES CARACTERISTIQUES "CADRE" DIFFERENTES'')')
     S    KNUME1,KNUME2
          CALL FAIPAR_MT (FA, FA%JPNIIL,INIMES,IREP,.FALSE.,CLMESS,
     S                 CLNSPR,CLACTI,.FALSE.)
        ELSEIF (IRANC1.NE.IRANC2) THEN
          INIMES=1
          WRITE (UNIT=CLMESS,FMT='(''REMARQUE: CADRES '''''',A,
     S           '''''' ET '''''',A,
     S           '''''' DISTINCTS MAIS DE CONTENU IDENTIQUE (UNITES'',
     S           I3,'' ET'',I3,'' )'')')
     S      FA%CNOMCA(IRANC1)(1:FA%NLCCAD(IRANC1)),
     S      FA%CNOMCA(IRANC2)(1:FA%NLCCAD(IRANC2)),KNUME1,KNUME2
          CALL FAIPAR_MT (FA, FA%JPNIIL,INIMES,IREP,.FALSE.,CLMESS,
     S                 CLNSPR,CLACTI,.FALSE.)
        ENDIF
C
      ENDIF
C
      IF (LLFATA) THEN
        INIMES=2
      ELSE
        INIMES=IXNVMS(IRANG(INUMFI))
      ENDIF
C
      IF (.NOT.LLFATA.AND.INIMES.NE.2)  THEN 
        IF (LHOOK) CALL DR_HOOK('FACOCH_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
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
      WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUME1='',I3,
     S       '', KNUME2='',I3,'', CDPREF='''''',A,'''''', KNIVAU='',I6,
     S       '', CDSUFF='''''',A,'''''''')') KREP,KNUME1,KNUME2,
     S   CLPREF(1:ILPREF),KNIVAU,CLSUFF(1:ILSUFF)
C
      IF (IREP.EQ.-112) THEN
        INUMRO=1000*KNUME1+KNUME2
      ELSE
        INUMRO=INUMER(INUMFI)
      ENDIF
C
      CALL FAIPAR_MT (FA, INUMRO,INIMES,IREP,LLFATA,CLMESS,
     S                CLNSPR, CLNOMA(1:ILNOMU),LLRLFI)
C
      IF (LHOOK) CALL DR_HOOK('FACOCH_MT',1,ZHOOK_HANDLE)
      END

