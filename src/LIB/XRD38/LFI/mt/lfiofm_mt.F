C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFIOFM_MT (LFI, KREP, KNUMER, KFACTM, LDOUVR )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-Programme permettant d'obtenir le Facteur Multiplicatif:
C
C     - effectif d'une Unite Logique FORTRAN deja ouverte pour le
C       logiciel de Fichiers Indexes *LFI*;
C     - prevu pour une Unite Logique FORTRAN, destinee a etre ouverte
C       ULTERIEUREMENT par le Logiciel de Fichiers Indexes *LFI*, en
C       supposant que l'on n'appelle pas ensuite LFIAFM ou LFIFMD
C       avant LFIOUV.
C
C       L'argument de sortie LDOUVR permet de savoir dans quel cas on se
C     trouve.
C**
C    ARGUMENTS : KREP   (Sortie) ==> Code-REPonse du sous-programme;
C                KNUMER (Entree) ==> NUMERo de l'unite logique;
C                KFACTM (Sortie) ==> FACteur Multiplicatif;
C                LDOUVR (Sortie) ==> Vrai si l'unite logique est deja
C                                    ouverte pour le logiciel LFI.
#ifndef f77
#include "precision.h"
#endif
C
      TYPE(LFICOM) :: LFI
      INTEGER KREP, KNUMER, KFACTM, IRANG, IREP, IRANFM, INIMES
C
      LOGICAL LLEXUL, LDOUVR
C
#include "lficom2.h"
#include "lficom_mt.h"
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFIOFM_MT',0,ZHOOK_HANDLE)
      CALL LFINUM_MT (LFI, KNUMER,IRANG)
      LDOUVR=IRANG.NE.0
C
      IF (LDOUVR) THEN
C
C        Unite logique deja ouverte pour le logiciel, on renvoie le
C     facteur multiplicatif effectif sous verrouillage eventuel.
C
        IF (LFI%LMULTI) CALL LFIVER_MT (LFI, LFI%VERRUE(IRANG),'ON')
        KFACTM=LFI%MFACTM(IRANG)
        IF (LFI%LMULTI) CALL LFIVER_MT (LFI, LFI%VERRUE(IRANG),'OFF')
      ELSE
C
C        Unite logique non (encore) ouverte pour le logiciel.
C
C        Controle de validite FORTRAN du Numero d'Unite Logique.
C
        INQUIRE (UNIT=KNUMER,EXIST=LLEXUL,ERR=901,IOSTAT=IREP)
C
        IF (.NOT.LLEXUL) THEN
          IREP=-30
          GOTO 1001
        ENDIF
C
C          On renvoie le facteur multiplicatif prevu,
C          sous verrouillage Global eventuel.
C
        IF (LFI%LMULTI) CALL LFIVER_MT (LFI, LFI%VERGLA,'ON')
        CALL LFIFMP_MT (LFI, KNUMER,IRANFM)
        KFACTM=LFI%MFACTU(IRANFM)
        IF (LFI%LMULTI) CALL LFIVER_MT (LFI, LFI%VERGLA,'OFF')
      ENDIF
C
      IREP=0
      GOTO 1001
C**
C     9.  - CI-DESSOUS, ETIQUETTE DE BRANCHEMENT EN CAS D'ERREUR INQUIRE
C-----------------------------------------------------------------------
C
  901 CONTINUE
      CLACTI='INQUIRE'
C
C      AU CAS OU, ON FORCE LE CODE-REPONSE ENTREE/SORTIE A ETRE POSITIF.
C
      IREP=IABS (IREP)
C**
C    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
C            VIA LE SOUS-PROGRAMME "LFIEMS" .
C-----------------------------------------------------------------------
C
 1001 CONTINUE
      KREP=IREP
      LLFATA=LLMOER (IREP,IRANG)
C
      IF (LLFATA) THEN
        INIMES=2
      ELSEIF (IRANG.EQ.0) THEN
        INIMES=LFI%NIMESG
      ELSE
        INIMES=IXNIMS (IRANG)
      ENDIF
C
      IF (INIMES.EQ.2) THEN
        CLNSPR='LFIOFM'
        WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,
     S         '', KFACTM='',I4,'', LDOUVR= '',L1)')
     S    KREP,KNUMER,KFACTM,LDOUVR
        CALL LFIEMS_MT (LFI, KNUMER,INIMES,IREP,LLFATA,
     S                  CLMESS,CLNSPR,CLACTI)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('LFIOFM_MT',1,ZHOOK_HANDLE)
      END

