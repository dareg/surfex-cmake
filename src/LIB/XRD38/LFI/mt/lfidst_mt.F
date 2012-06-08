C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFIDST_MT (LFI, KREP, KRANG, KRANIE, CDSTRU, KLONG )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        SOUS-PROGRAMME *INTERNE* DU LOGICIEL DE FICHIERS INDEXES LFI;
C     Decorticage de la STructure decrivant un article logique
C     de donnees, en vue de son import ou export ulterieur.
C**
C    ARGUMENTS : KREP   (SORTIE) ==> CODE-REPONSE DU SOUS-PROGRAMME;
C                KRANG  (ENTREE) ==> RANG ( DANS LA TABLE *LFI%NUMERO* )
C                                    DE L'UNITE LOGIQUE CONCERNEE;
C                KRANIE (Entree) ==> Rang dans les tables d'import/exp;
C                CDSTRU (Entree) ==> Structure interne de cet article;
C                KLONG  (Entree) ==> Longueur (mots) de l'article.
C
C     Les syntaxes autorisees pour CDSTRU sont decrites ci-dessous.
C     La presence de crochets [ ] indique le cote optionnel, mais dans
C     tous les cas ce cote optionnel est reserve a la toute derniere
C     partie de la description. Si l'argument optionnel est present,
C     il doit etre coherent avec la longueur effective de l'article. 
C
C     'type [nbre]' ==> article homogene,ex: 'i', 'r 20'
C
C     'type_1 nbre_1 ... type_n [nbre_n]' ==> juxtaposition de types,ex:
C                                             'i 2 r', 'i 3 r 2 c 80000'
C
C     '(type_1 nbre_1 ... type_n nbre_n) [nbre]' ==> boucle,ex:
C                                                    '(i 1 r 1)'
C
C     ou une juxtaposition des possibilites ci-dessus.
C
C     Les blancs ne sont pas obligatoires, ils sont neutres et utilises
C     ci-dessus pour des questions de clarte.
C
#ifndef f77
#include "precision.h"
#endif
C
      TYPE(LFICOM) :: LFI
      CHARACTER CDSTRU*(*), CLSTRU*(LFI%JPNCPN)
C
      INTEGER KREP, KRANG, KRANIE, KLONG, IRANG, ILCLST, ILUSTR
      INTEGER J, INMOTS, INIPAR, IDEXPL, INGROU, IJ, IDECAL, IPOSTY
      INTEGER INIMES, INUMER, IDECGR, INTYPE
C
#include "lficom2.h"
#include "lficom_mt.h"
C**
C     1.  -  CONTROLES DES PARAMETRES D'APPEL, PUIS INITIALISATIONS.
C-----------------------------------------------------------------------
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFIDST_MT',0,ZHOOK_HANDLE)
      ILUSTR=LEN (CDSTRU)
      CLSTRU=' '
      ILCLST=1
C
      IF (KRANG.LE.0.OR.KRANG.GT.LFI%JPNXFI.OR.ILUSTR.LE.0.OR.
     S    KRANIE.LE.0.OR.KRANIE.GT.LFI%JPIMEX) THEN
        KREP=-16
        ILCLST=MIN0 (LFI%JPNCPN,LEN (CLSTRU))
        CLSTRU=LFI%CHINCO(:ILCLST)
        GOTO 1001
      ENDIF
C
      IRANG=KRANG
      KREP=0
      ILCLST=MIN0 (ILUSTR,LEN (CLSTRU))
      CLSTRU=CDSTRU(:ILCLST)
      INMOTS=0
      INIPAR=0
      INGROU=0
      IDECGR=0
      IDEXPL=LFI%NREXPL(LFI%NAEXPL(KRANG),KRANIE)
      IJ=0
C**
C     2.  -  DECORTICAGE PROPREMENT DIT.
C-----------------------------------------------------------------------
C*
C     2.1 -  DEBUT "BOUCLE" SUR LES GROUPES.
C-----------------------------------------------------------------------
C
  211 CONTINUE
C
      IDECAL=IJ
C
      IF (IDECAL.GE.ILUSTR) GOTO 301
C
      DO 212 J=IDECAL+1,ILUSTR
C
      IF (CDSTRU(J:J).EQ.'(') THEN
C
        INGROU=INGROU+1
        INIPAR=INIPAR+1
C
        IF (INIPAR.GT.1) THEN
          KREP=-40
          GOTO 1001
        ENDIF
C
        IJ=J
C
      ELSEIF (CDSTRU(J:J).NE.' ') THEN
C
        IPOSTY=INDEX (LFI%CTYPMX,CDSTRU(J:J))
C
        IF (IPOSTY.EQ.0) THEN
          KREP=-40
          GOTO 1001
        ENDIF
C
        INGROU=INGROU+1
        IJ=J-1
C
      ENDIF
C      
  212 CONTINUE
C*
C     2.2 -  DEBUT "BOUCLE" SUR LES TYPES.
C-----------------------------------------------------------------------
C
      INTYPE=0
C
  221 CONTINUE
C
      IDECAL=IJ
C
      DO 227 J=IDECAL+1,ILUSTR
C
      IF (CDSTRU(J:J).EQ.')') THEN
C
        INIPAR=INIPAR-1
C
        IF (INTYPE.EQ.0.OR.INIPAR.NE.0) THEN
          KREP=-40
          GOTO 1001
        ENDIF
C
        IJ=J
        
        GOTO 211
C
      ELSEIF (CDSTRU(J:J).NE.' ') THEN
C
        IPOSTY=INDEX (LFI%CTYPMX,CDSTRU(J:J))
C
        IF (IPOSTY.EQ.0) THEN
          KREP=-40
          GOTO 1001
        ENDIF
C
        INTYPE=INTYPE+1
C
        IF ((IDEXPL+2).GT.LFI%JPDEXP) THEN
          KREP=-42
          GOTO 1001
        ENDIF
C
        LFI%NDEXPL(IDEXPL+1,KRANIE)=IPOSTY
C
        IF (J.EQ.ILUSTR) THEN
C
          IF (INIPAR.EQ.0) THEN
            LFI%NDEXPL(IDEXPL+2,KRANIE)=LFI%JPNIL
            IDEXPL=IDEXPL+2
            GOTO 301
          ELSE
            KREP=-40
            GOTO 1001
          ENDIF
C
        ELSE
C
          CALL LFICHI_MT (LFI, KREP,CDSTRU(J+1:ILUSTR),
     S                    LFI%NDEXPL(IDEXPL+2,KRANIE),
     S                 IJ)
          IF (KREP.NE.0) GOTO 1001
C
          IDEXPL=IDEXPL+2
C
        ENDIF
C 
      ENDIF
C 
  227 CONTINUE
C


C
  301 CONTINUE
C


C**
C    10.  -  PHASE TERMINALE : MESSAGERIE INTERNE EVENTUELLE,
C            VIA LE SOUS-PROGRAMME "LFIEMS", PUIS RETOUR.
C-----------------------------------------------------------------------
C
 1001 CONTINUE
      LLFATA=LLMOER (KREP,KRANG)
C
      IF (LFI%LMISOP.OR.LLFATA) THEN
        INUMER=LFI%NUMERO(KRANG)
        INIMES=2
        CLNSPR='LFIDST'
        WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KRANG='',I3,
     S         '''''', CDSTRU='''''',A,'''''', KLONG='',I7)')
     S    KREP,KRANG,CLSTRU(:ILCLST),KLONG
        CALL LFIEMS_MT (LFI, INUMER,INIMES,KREP,.FALSE.,
     S                  CLMESS,CLNSPR,CLACTI)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('LFIDST_MT',1,ZHOOK_HANDLE)
      END

