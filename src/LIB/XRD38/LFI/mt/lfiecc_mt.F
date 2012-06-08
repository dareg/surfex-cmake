C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFIECC_MT (LFI, KREP, KNUMER, KREC, CDTAB, KNBECR, 
     S                      KFACTM, KRETIN )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-programme charge des Ecritures de Chaines de Caracteres
C     du logiciel LFI (articles d'index "noms").
C**
C     Arguments: KREP   (Sortie) ==> Code-reponse ( zero si OK; code-
C                                    reponse du "WRITE" FORTRAN sinon);
C                KNUMER (Entree) ==> NUMERo d'unite logique FORTRAN;
C                KREC   (Entree) ==> Numero d'enregistrement a ecrire;
C                KTAB   (Sortie) ==> Zone a ecrire, de Longueur
C                                    LFI%JPLARD*KFACTM *mots*;
C                KNBECR (Entree  ==> Compteur d'ECRitures sur l'unite;
C                       +Sortie)
C                KFACTM (Entree) ==> FACteur Multiplicatif LFI de
C                                    l'unite logique;
C                KRETIN (Sortie) ==> Code-retour interne.
C
C
      TYPE(LFICOM) :: LFI
      INTEGER KREP, KNUMER, KREC, KNBECR, KFACTM, KRETIN
C
      CHARACTER CDTAB (LFI%JPNXNA*KFACTM)*(LFI%JPNCPN)
#include "lficom_mt.h"
C
C        ECRITURE .
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFIECC_MT',0,ZHOOK_HANDLE)
      WRITE (UNIT=KNUMER,REC=KREC,ERR=901,IOSTAT=KREP) CDTAB
C
      IF (LFI%LMISOP) THEN
        WRITE (UNIT=LFI%NULOUT,FMT=*)
     S          '+++++ LFIECC - WRITE / ',KNUMER,', REC = ',KREC,
     S          ' +++++'
      ENDIF
C
      KNBECR=KNBECR+1
      KRETIN=0
      GOTO 1001
C
  901 CONTINUE
      KRETIN=1
C
 1001 CONTINUE
C
      IF (LHOOK) CALL DR_HOOK('LFIECC_MT',1,ZHOOK_HANDLE)
      END

