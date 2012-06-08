C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFISFM_MT (LFI, KREP, KNUMER )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-Programme Suprimant un Facteur Multiplicatif
C     d'une Unite Logique FORTRAN, qui a ete fermee PRECEDEMMENT
C     par le Logiciel de Fichiers Indexes *LFI* .
C     (ou du moins, n'est pas ouverte pour ce logiciel)
C
C        Ce sous-programme permet de faire de la place dans les tables
C     decrivant les associations Unite Logique/facteur Multiplicatif.
C**
C    ARGUMENTS : KREP   (Sortie) ==> Code-REPonse du sous-programme;
C                KNUMER (Entree) ==> NUMERo de l'unite logique.
C
#ifndef f77
#include "precision.h"
#endif
C
      TYPE(LFICOM) :: LFI
      INTEGER KREP, KNUMER, IRANG, IREP, IRANFM, INIMES, IFACTM, J
C
      LOGICAL LLVERG, LLEXUL
C
#include "lficom2.h"
#include "lficom_mt.h"
C**
C     1.  -  CONTROLES DES PARAMETRES D'APPEL, INITIALISATIONS.
C-----------------------------------------------------------------------
C
C        Appel a LFINUM, permettant (le cas echeant) l'initialisation
C     variables globales du logiciel a la 1ere utilisation.
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFISFM_MT',0,ZHOOK_HANDLE)
      IFACTM=0
      LLVERG=.FALSE.
      CALL LFINUM_MT (LFI, KNUMER,IRANG)
C
      IF (IRANG.NE.0) THEN
        IREP=-5
        GOTO 1001
      ENDIF
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
C              Verrouillage Global eventuel.
C
       IF (LFI%LMULTI) CALL LFIVER_MT (LFI, LFI%VERGLA,'ON')
      LLVERG=LFI%LMULTI
C**
C     2.  -  TRAVAIL EFFECTIF SUR LES TABLES DECRIVANT LES ASSOCIATIONS
C            UNITES LOGIQUES/FACTEURS.
C-----------------------------------------------------------------------
C
      CALL LFIFMP_MT (LFI, KNUMER,IRANFM)
C
      IF (IRANFM.EQ.0) THEN
C
C          Unite logique non trouvee dans la table *LFI%MULOFM*.
C
        IREP=-31
        GOTO 1001
      ENDIF
C
      IFACTM=LFI%MFACTU(IRANFM)
      LFI%NULOFM=LFI%NULOFM-1
C
      DO 201 J=IRANFM,LFI%NULOFM
      LFI%MFACTU(J)=LFI%MFACTU(J+1)
      LFI%MULOFM(J)=LFI%MULOFM(J+1)
  201 CONTINUE
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
      IF (LLVERG) CALL LFIVER_MT (LFI, LFI%VERGLA,'OFF')
C
      IF (LLFATA) THEN
        INIMES=2
      ELSEIF (IRANG.EQ.0) THEN
        INIMES=LFI%NIMESG
      ELSE
        INIMES=IXNIMS (IRANG)
      ENDIF
C
      IF (INIMES.EQ.0)  THEN 
        IF (LHOOK) CALL DR_HOOK('LFISFM_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
      CLNSPR='LFISFM'
C
      IF (INIMES.EQ.2) THEN
        WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3)')
     S       KREP,KNUMER
        CALL LFIEMS_MT (LFI, KNUMER,INIMES,IREP,LLFATA,
     S                  CLMESS,CLNSPR,CLACTI)
      ENDIF
C
      IF (LFI%LFRANC) THEN
        WRITE (UNIT=CLMESS,FMT=
     S         '(''Suppression du Facteur Multiplicatif'',I3,
     S           '', Unite Logique'',I3)') IFACTM,KNUMER
      ELSE
        WRITE (UNIT=CLMESS,FMT=
     S         '(''Multiply Factor'',I3,
     S           '' suppressed, Logical Unit'',I3)') IFACTM,KNUMER
      ENDIF
C
      CALL LFIEMS_MT (LFI, KNUMER,INIMES,IREP,.FALSE.,
     S                CLMESS,CLNSPR,CLACTI)
C
      IF (LHOOK) CALL DR_HOOK('LFISFM_MT',1,ZHOOK_HANDLE)
      END

