C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFIAFM_MT (LFI, KREP, KNUMER, KFACTM )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-Programme permettant d'Attribuer un Facteur Multiplicatif
C     a une Unite Logique FORTRAN, destinee a etre ouverte
C     ULTERIEUREMENT par le Logiciel de Fichiers Indexes *LFI* .
C        Lors de cette ouverture ulterieure, LFIOUV essaiera de traiter
C     l'unite logique consideree comme un fichier a acces direct
C     non formatte de longueur d'article "Physique" LFI%JPLARD*KFACTM mots.
C**
C    ARGUMENTS : KREP   (Sortie) ==> Code-REPonse du sous-programme;
C                KNUMER (Entree) ==> NUMero de l'unite logique;
C                KFACTM (Entree) ==> FACteur Multiplicatif a attribuer.
C
#ifndef f77
#include "precision.h"
#endif
C
      TYPE(LFICOM) :: LFI
      INTEGER KREP, KNUMER, KFACTM, IRANG, IREP, IRANFM, INIMES, IFACTM
C
      LOGICAL LLVERG, LLEXUL
C
#include "lficom2.h"
#include "lficom_mt.h"
C**
C     1.  -  CONTROLES DES PARAMETRES D'APPEL, INITIALISATIONS.
C-----------------------------------------------------------------------
C
C        Appel legerement anticipe a LFINUM, permettant une initialisa-
C     tion des variables globales du logiciel a la 1ere utilisation.
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFIAFM_MT',0,ZHOOK_HANDLE)
      IFACTM=KFACTM
      LLVERG=.FALSE.
      CALL LFINUM_MT (LFI, KNUMER,IRANG)
C
      IF (KFACTM.LE.0) THEN
        IREP=-14
        GOTO 1001
      ELSEIF (KFACTM.GT.LFI%JPFACX) THEN
        IREP=-28
        GOTO 1001
      ELSEIF (IRANG.NE.0) THEN
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
      IF (IRANFM.NE.0) THEN
C
C          Redefinition du facteur multiplicatif.
C
        IFACTM=LFI%MFACTU(IRANFM)
      ELSEIF (LFI%NULOFM.GE.LFI%JPXUFM) THEN
C
C          Tables pleines...
C
        IREP=-29
        GOTO 1001
      ELSE
C
C          Cas standard.
C
        LFI%NULOFM=LFI%NULOFM+1
        IRANFM=LFI%NULOFM
        LFI%MULOFM(IRANFM)=KNUMER
        IFACTM=KFACTM
      ENDIF
C
      LFI%MFACTU(IRANFM)=KFACTM
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
        IF (LHOOK) CALL DR_HOOK('LFIAFM_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
      CLNSPR='LFIAFM'
C
      IF (INIMES.EQ.2) THEN
        WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3,
     S         '', KFACTM='',I4)') KREP,KNUMER,KFACTM
        CALL LFIEMS_MT (LFI, KNUMER,INIMES,IREP,LLFATA,
     S                  CLMESS,CLNSPR,CLACTI)
      ENDIF
C
      IF (IFACTM.EQ.KFACTM) THEN
C
        IF (LFI%LFRANC) THEN
          WRITE (UNIT=CLMESS,FMT=
     S           '(''Attribution du Facteur Multiplicatif'',I3,
     S             '' a l''''Unite Logique'',I3)') KFACTM,KNUMER
        ELSE
          WRITE (UNIT=CLMESS,FMT='(''Multiply Factor'',I3,
     S           '' specified for Logical Unit'',
     S             I3)') KFACTM,KNUMER
        ENDIF
C
        CALL LFIEMS_MT (LFI, KNUMER,INIMES,IREP,.FALSE.,CLMESS,
     S                  CLNSPR,CLACTI)
      ELSE
C
        IF (LFI%LFRANC) THEN
          WRITE (UNIT=CLMESS,FMT='(''Unite Logique'',I3,
     S           '': *NOUVEAU* Facteur Multiplicatif attribue='',I3)')
     S    KNUMER,KFACTM
        ELSE
          WRITE (UNIT=CLMESS,FMT='(''Logical Unit'',I3,
     S           '': *NEW* Multiply Factor specified='',I3)')
     S    KNUMER,KFACTM
        ENDIF
C
        CALL LFIEMS_MT (LFI, KNUMER,INIMES,IREP,.FALSE.,
     S                  CLMESS,CLNSPR,CLACTI)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('LFIAFM_MT',1,ZHOOK_HANDLE)
      END

