C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFILDO_MT (LFI, KREP, KNUMER, KREC, KTAB, KNBLEC, 
     S                      KFACTM, KRETIN )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-programme charge des Lectures de DOnnees du logiciel LFI,
C     *SAUF* pour les articles d'index de type caractere.
C**
C     Arguments: KREP   (Sortie) ==> Code-reponse ( zero si OK; code-
C                                    reponse du "READ" FORTRAN sinon);
C                KNUMER (Entree) ==> NUMERo d'unite logique FORTRAN;
C                KREC   (Entree) ==> Numero d'enregistrement a lire;
C                KTAB   (Sortie) ==> Zone a lire, de Longueur
C                                    LFI%JPLARD*KFACTM *mots*;
C                KNBLEC (Entree  ==> Compteur de LECtures sur l'unite;
C                       +Sortie)
C                KFACTM (Entree) ==> FACteur Multiplicatif LFI de
C                                    l'unite logique;
C                KRETIN (Sortie) ==> Code-retour interne.
C
C
      TYPE(LFICOM) :: LFI
      INTEGER KREP, KNUMER, KREC, KNBLEC, KFACTM, KRETIN
#ifndef f77
#include "precision.h"
C
      INTEGER (KIND=JPDBLE)  KTAB (LFI%JPLARD*KFACTM)
#else
      INTEGER KTAB (LFI%JPLARD*KFACTM)
#endif
#include "lficom_mt.h"
C
C        LECTURE .
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFILDO_MT',0,ZHOOK_HANDLE)
      READ (UNIT=KNUMER,REC=KREC,ERR=901,IOSTAT=KREP) KTAB
C
      IF (LFI%LMISOP) THEN
        WRITE (UNIT=LFI%NULOUT,FMT=*)
     S          '+++++ LFILDO - READ / ',KNUMER,', REC = ',KREC,
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
      IF (LHOOK) CALL DR_HOOK('LFILDO_MT',1,ZHOOK_HANDLE)
      END

