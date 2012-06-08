C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFILCC_MT (LFI, KREP, KNUMER, KREC, CDTAB, 
     S                      KNBLEC, KFACTM, KRETIN )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-programme charge des Lectures de Chaines de Caracteres
C     du logiciel LFI (articles d'index "noms").
C**
C     Arguments: KREP   (Sortie) ==> Code-reponse ( zero si OK; code-
C                                    reponse du "READ" FORTRAN sinon);
C                KNUMER (Entree) ==> NUMERo d'unite logique FORTRAN;
C                KREC   (Entree) ==> Numero d'enregistrement a lire;
C                KTAB   (Sortie) ==> Zone a lire, de Longueur
C                                    LFI%JPLARD*KFACTM *mots*;
C                KNBLEC (Entree) ==> Compteur de LECtures sur l'unite;
C                       +Sortie)
C                KFACTM (Entree) ==> FACteur Multiplicatif LFI de
C                                    l'unite logique;
C                KRETIN (Sortie) ==> Code-retour interne.
C
C
      TYPE(LFICOM) :: LFI
      INTEGER KREP, KNUMER, KREC, KNBLEC, KFACTM, KRETIN
C
      CHARACTER CDTAB (LFI%JPNXNA*KFACTM)*(LFI%JPNCPN)
#include "lficom_mt.h"
C
C        LECTURE .
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFILCC_MT',0,ZHOOK_HANDLE)
      READ (UNIT=KNUMER,REC=KREC,ERR=901,IOSTAT=KREP) CDTAB
C
      IF (LFI%LMISOP) THEN
        WRITE (UNIT=LFI%NULOUT,FMT=*)
     S          '+++++ LFILCC - READ / ',KNUMER,', REC = ',KREC,
     S          ' +++++'
      ENDIF
C
      KNBLEC=KNBLEC+1
      KRETIN=0
      GOTO 1001
C
  901 CONTINUE
      KRETIN=2
C
 1001 CONTINUE
C
      IF (LHOOK) CALL DR_HOOK('LFILCC_MT',1,ZHOOK_HANDLE)
      END

