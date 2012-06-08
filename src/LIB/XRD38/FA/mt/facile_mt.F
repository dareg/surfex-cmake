C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FACILE_MT (FA,  KREP, KNUMER, CDPREF, KNIVAU, CDSUFF, 
     S                      KCHAMP, LDCOSP )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C      Sous-programme de LECTURE d'un CHAMP HORIZONTAL sur un fichier
C     ARPEGE.
C       ( Champ d'Interet en LEcture )
C**
C    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
C                KNUMER (Entree) ==> Numero de l'unite logique;
C                CDPREF (Entree) ==> Prefixe eventuel du nom d'article;
C                KNIVAU (Entree) ==> Niveau vertical eventuel;
C                CDSUFF (Entree) ==> Suffixe eventuel du nom d'article;
C    ( Tableau ) KCHAMP (Sortie) ==> Valeurs REELLES du champ lu;
C                LDCOSP (Entree) ==> Vrai si le champ est represente
C                                    par des coefficients spectraux.
C     MODIF:
C     JM AUDOIN GMAP/EXT 10/05/95 intro de IVALC3 pour eviter ecrasement
C     D  PARADIS TTI/DEV 12/10/98 partie controle et decodage de l'article
C                                 demenagee dans un ss-prg a usage interne
C                                 du logiciel (FADECI).
C     D  PARADIS DSI/DEV 15/04/04 nettoyage code + declaration IVALCO en
C                                 ALLOCATABLE
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KREP, KNUMER, KNIVAU
C
      INTEGER IREP, ILPRFU, ILSUFU, ILNOMU, ILONGA, IRANG, INIMES
      INTEGER ILPREF, ILSUFF, IPOSEX, IRANGC
C
      INTEGER (KIND=JPDBLE) KCHAMP (*)
      INTEGER (KIND=JPDBLE), ALLOCATABLE :: IVALCO(:)
      INTEGER IB1PAR (FA%JPLB1P)
C
      LOGICAL LLVERF, LLRLFI, LDCOSP, LLNOMU
C
      CHARACTER CDPREF*(*), CDSUFF*(*)
      CHARACTER CLPREF*(FA%JPXNOM), CLSUFF*(FA%JPXSUF)
C
#include "facom2.h"
#include "facom_mt.h"
C**
C     1.  -  CONTROLES ET INITIALISATIONS.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FACILE_MT',0,ZHOOK_HANDLE)
      LLVERF=.FALSE.
      LLRLFI=.FALSE.
      LLNOMU=.FALSE.
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
     S             IB1PAR(6),ILPRFU,ILSUFU,ILNOMU)
      IF (IREP.NE.0) GOTO 1001
      LLNOMU=.TRUE.
C**
C     3.  -  LECTURE DE L'ARTICLE SUR LE FICHIER
C-----------------------------------------------------------------------
C
      CALL LFINFO_MT (FA%LFI, IREP,KNUMER,CLNOMA(1:ILNOMU),
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
C     ALLOCATE (IVALCO (2*ILONGA))
      ALLOCATE (IVALCO (ILONGA))
      CALL LFILEC_MT (FA%LFI, IREP,KNUMER,CLNOMA(1:ILNOMU),
     S             IVALCO,ILONGA)
      LLRLFI=IREP.NE.0
      IF (LLRLFI) GOTO 1001
C
C**
C     4.  -  CONTROLES ET DECODAGE DE L'ARTICLE
C----------------------------------------------
C
C  Controle de l'homogeneite du type de rangement de coeff. spectraux
C  parmi les champs lus/ecrits: ces champs compactes avec
C  FA%NIGRIB=-1 ou 3 doivent etre ranges comme dans le modele ("verticalement"
C  soit selon des colonnes JM=cst consecutives) et contrairement si compactes
C  avec FA%NIGRIB= 0,1 ou 2.
C 
      IRANGC=FA%NUCADR(IRANG)
      IF (LDCOSP) THEN
        IF (IVALCO(1).EQ.-1.OR.IVALCO(1).EQ.3) THEN
          FA%NRASVE(IRANG)=FA%NRASVE(IRANG)+1
          IF (FA%NRASVE(IRANG).EQ.1.AND.FA%NRASHO(IRANG).GT.0) THEN
            WRITE(FA%NULOUT,*)
     S      '------------------------------------------------'
            WRITE(FA%NULOUT,*)' FACILE :  WARNING !!!!!           '
            WRITE(FA%NULOUT,*)' Un champ de coeff. spectraux avec'
            WRITE(FA%NULOUT,*)
     S      ' rangement type modele va etre lu alors que'
            WRITE(FA%NULOUT,*)
     S      ' d''autres champs spec. ont un rangt different.'
            WRITE(FA%NULOUT,*)
     S      ' ***  Prenez en compte cette heterogeneite!  ***'
            WRITE(FA%NULOUT,*)
     S      '------------------------------------------------'
          ENDIF
        ELSEIF (IVALCO(1).GE.0.AND.IVALCO(1).LE.2) THEN
          FA%NRASHO(IRANG)=FA%NRASHO(IRANG)+1
          IF (FA%NRASHO(IRANG).EQ.1.AND.FA%NRASVE(IRANG).GT.0) THEN
            WRITE(FA%NULOUT,*)
     S      '------------------------------------------------'
            WRITE(FA%NULOUT,*)' FACILE :  WARNING !!!!!           '
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
      IF (IVALCO(1).EQ.3) THEN
C Cas d'un champ gribe avec GRIBEX
        CALL FADECX_MT (FA, IREP,IRANG,CLNOMA,
     S               IVALCO,ILONGA,KCHAMP,LDCOSP)
      ELSE
        CALL FADECI_MT (FA, IREP,IRANG,CLNOMA,
     S               IVALCO,ILONGA,KCHAMP,LDCOSP)
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
      LLFATA=LLMOER (IREP,IRANG)
C
C        Deverrouillage eventuel du fichier.
C
      IF (LLVERF) CALL LFIVER_MT (FA%LFI, FA%VRFICH(IRANG),'OFF')

      IF (LLFATA) THEN
        INIMES=2
      ELSE
        INIMES=IXNVMS(IRANG)
      ENDIF
C
      IF (.NOT.LLFATA.AND.INIMES.NE.2)  THEN 
        IF (LHOOK) CALL DR_HOOK('FACILE_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='FACILE'
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
      IF (LHOOK) CALL DR_HOOK('FACILE_MT',1,ZHOOK_HANDLE)
      END






