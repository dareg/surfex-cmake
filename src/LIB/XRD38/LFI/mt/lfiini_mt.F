C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFIINI_MT (LFI, KOPTIO )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        CE SOUS-PROGRAMME EST CHARGE DES INITIALISATIONS DU LOGICIEL
C     DE FICHIERS INDEXES LFI .
C**
C        ARGUMENT : KOPTIO  ==> OPTION CONCERNANT LE MODE D'UTILISATION.
C                  (ENTREE)     (MULTI-TACHE OU NON)
C     VALEURS POSSIBLES : 0 ==> MODE MONO-TACHE PRESCRIT
C                         1 ==> MODE MULTI-TACHE PRESCRIT
C                         2 ==> UTILISATION DU MODE PAR DEFAUT SI C'EST
C                               LE PREMIER APPEL; SINON ON GARDE LE MODE
C                               PRESCRIT ANTERIEUREMENT .
C
#ifndef f77
#include "precision.h"
#endif
C
      TYPE(LFICOM) :: LFI
      INTEGER KOPTIO, JNPAGE, J, JRANG, IREP, INIMES, INUMER, IRGPIM
C
      LOGICAL LLNMUL, LLASGN, LLREL
#include "lficom_mt.h"
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFIINI_MT',0,ZHOOK_HANDLE)
      IF (LFI%LFIINI_LLPREA) THEN
C
C         C'EST LE PREMIER APPEL AU SOUS-PROGRAMME - INITIALISATIONS .
C
        LFI%NBFIOU=0
        LFI%NFACTM=0
        LFI%NULOFM=0
        LFI%NUIMEX=0
        LFI%NERFAG=1
        LFI%NISTAG=1
        LFI%NIMESG=1
        LFI%NPISAF=0
        LFI%LTAMLG=.FALSE.
        LFI%LTAMEG=.TRUE.
        LFI%LMISOP=.FALSE.
        LFI%LFRANC=.TRUE.
        LFI%LERFAT(0)=.TRUE.
        LFI%NIVMES(0)=0
C
C          L'indice zero dans LFI%MFACTU correspond au facteur multiplicatif
C     "par defaut" .
C
#ifdef HIGHRES
        LFI%MFACTU(0)=6
#else
        LFI%MFACTU(0)=1
#endif
        LFI%MFACTM(0)=0
        LFI%NUMAPH(0)=LFI%JPNIL
C
        LFI%CFGMXD(0)='locale'
        LFI%NBMOSD(0)=LFI%JPNBIM
        LFI%NBCASD(0)=LFI%JPNBIC
        LFI%CFGMXD(1)='ieee'
        LFI%NBMOSD(1)=32
        LFI%NBCASD(1)=8
        LFI%CFGMXD(2)='ultrix'
        LFI%NBMOSD(2)=32
        LFI%NBCASD(2)=8
        LFI%CFGMXD(3)='ieee_dp'
        LFI%NBMOSD(3)=64
        LFI%NBCASD(3)=8
        LFI%CFGMXD(4)='ultrix_dp'
        LFI%NBMOSD(4)=64
        LFI%NBCASD(4)=8
        LFI%CTYPMX='ircdl'
C
        DO 101 JNPAGE=1,LFI%JPNXPI
        LFI%MCOPIF(JNPAGE)=LFI%JPNIL
        LFI%MRGPIF(JNPAGE)=LFI%JPNIL
        LFI%LECRPI(JNPAGE,1)=.FALSE.
        LFI%LECRPI(JNPAGE,2)=.FALSE.
  101   CONTINUE
C
        DO 103 J=1,LFI%JPNPIA
C
        DO 102 JRANG=1,LFI%JPNXFI
        IRGPIM=JRANG+(J-1)*LFI%JPNXFI
        LFI%MCOPIF(IRGPIM)=JRANG
        LFI%MRGPIM(J,JRANG)=IRGPIM
  102   CONTINUE
C
  103   CONTINUE
C
        DO 104 JRANG=1,LFI%JPNXFI
        LFI%MRGPIF(JRANG)=1
        LFI%NUMERO(JRANG)=LFI%JPNIL
        LFI%NUMAPH(JRANG)=LFI%JPNIL
  104   CONTINUE
C
        DO 105 J=1,LFI%JPNCPN
        LFI%CHINCO(J:J)='?'
  105   CONTINUE
C
        DO 106 J=1,LFI%JPIMEX
        LFI%MNUIEX(J)=LFI%JPNIL
        LFI%NREXPL(0,J)=0
  106   CONTINUE
C
        LFI%LFIINI_LLPREA=.FALSE.
        LLNMUL=(KOPTIO.EQ.1).OR.(KOPTIO.EQ.2.AND.LFI%LFIINI_LLDEFM)
        LLASGN=LLNMUL
        LLREL=.FALSE.
C
      ELSEIF (KOPTIO.EQ.2) THEN
C
C         CE N'EST PAS LE PREMIER APPEL, MAIS COMME L'ARGUMENT VAUT 2,
C         ON LAISSE LES CHOSES EN PLACE .
C
        LLNMUL=LFI%LMULTI
        LLASGN=.FALSE.
        LLREL =.FALSE.
      ELSE
C
C     CE N'EST PAS LE PREMIER APPEL ET LE MODE EST PASSE 'EXPLICITEMENT'
C
        LLNMUL=KOPTIO.EQ.1
        LLASGN=LLNMUL.AND.(.NOT.LFI%LMULTI)
        LLREL =(.NOT.LLNMUL).AND.LFI%LMULTI
C
        IF ((LLASGN.OR.LLREL).AND.LFI%NBFIOU.NE.0) THEN
          IREP=-4
          GOTO 1001
        ENDIF
C
      ENDIF
C
      LFI%LMULTI=LLNMUL
C
C        Le controle de coherence ci-dessous ne peut etre fait qu'apres
C     l'initialisation eventuelle des variables globales du logiciel.
C
      IF (KOPTIO.LT.0.OR.KOPTIO.GT.2) THEN
        IREP=-2
        GOTO 1001
      ENDIF
C
      IREP=0
C
      IF (LLASGN) THEN
        CALL LFIVER_MT (LFI, LFI%VERGLA,'ASGN')
      ELSEIF (LLREL) THEN
        CALL LFIVER_MT (LFI, LFI%VERGLA,'REL')
      ENDIF
C
 1001 CONTINUE
C
C        MESSAGERIE EVENTUELLE, AVEC ABORT SI NECESSAIRE .
C
      LLFATA=IREP.NE.0.AND.LFI%NERFAG.NE.2
C
      IF (LLFATA) THEN
        INIMES=2
      ELSEIF (IREP.NE.0) THEN
        INIMES=0
      ELSEIF (LFI%NIMESG.EQ.2.OR.(LFI%NIMESG.EQ.1.AND.KOPTIO.NE.2)) THEN
        INIMES=LFI%NIMESG
      ELSE
        IF (LHOOK) CALL DR_HOOK('LFIINI_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      CLNSPR='LFIINI'
      INUMER=LFI%JPNIL
C
      IF (MAX0 (INIMES,LFI%NIMESG).EQ.2) THEN
C
        IF (LFI%LFRANC) THEN
          WRITE (UNIT=CLMESS,
     S           FMT='(''KOPTIO='',I5,'', CODE INTERNE='',
     S           I4)') KOPTIO,IREP
        ELSE
          WRITE (UNIT=CLMESS,
     S           FMT='(''KOPTIO='',I5,'', INTERNAL CODE='',
     S           I4)') KOPTIO,IREP
        ENDIF
C
        IF (INIMES.NE.2) CALL LFIEMS_MT (LFI, INUMER,LFI%NIMESG,IREP,
     S                                   .FALSE.,CLMESS,
     S                                   CLNSPR,CLACTI)
      ENDIF
C
      CALL LFIEMS_MT (LFI, INUMER,INIMES,IREP,LLFATA,
     S                CLMESS,CLNSPR,CLACTI)
C
      IF (INIMES.GE.1.AND.KOPTIO.NE.2) THEN
C
C        Cette messagerie de niveau 1 n'est pas emise en cas d'erreur,
C     meme non fatale.
C
        IF (LFI%LFRANC) THEN
C
          IF (KOPTIO.EQ.0) THEN
            CLMESS='$ Mode MONO-TACHE Prescrit explicitement... $'
          ELSE
            CLMESS='$ Mode MULTI-TACHE Prescrit explicitement... $'
          ENDIF
C
        ELSE
C
          IF (KOPTIO.EQ.0) THEN
            CLMESS='$ MONO-TASKING Mode explicitely Specified... $'
          ELSE
            CLMESS='$ MULTI-TASKING Mode explicitely Specified... $'
          ENDIF
C
        ENDIF
C
        CALL LFIEMS_MT (LFI, INUMER,INIMES,IREP,.FALSE.,
     S                  CLMESS,CLNSPR,CLACTI)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('LFIINI_MT',1,ZHOOK_HANDLE)
      END

