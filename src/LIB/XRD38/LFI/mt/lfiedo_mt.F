C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFIEDO_MT (LFI, KREP, KNUMER, KREC, KTAB, KNBECR, 
     S                      KFACTM, KRETIN )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Sous-programme charge des Ecritures de DOnnees du logiciel LFI,
C     *SAUF* pour les articles d'index de type caractere.
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
#ifndef f77
#include "precision.h"
C
      TYPE(LFICOM) :: LFI
      INTEGER KREP, KNUMER, KREC, KNBECR, KFACTM, KRETIN
C
      INTEGER (KIND=JPDBLE)  KTAB (LFI%JPLARD*KFACTM)
#else
      INTEGER KTAB (LFI%JPLARD*KFACTM)
#endif
#include "lficom_mt.h"
C
C        ECRITURE .
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFIEDO_MT',0,ZHOOK_HANDLE)
      WRITE (UNIT=KNUMER,REC=KREC,ERR=901,IOSTAT=KREP) KTAB
C
      IF (LFI%LMISOP) THEN
        WRITE (UNIT=LFI%NULOUT,FMT=*)
     S          '+++++ LFIEDO - WRITE / ',KNUMER,', REC = ',KREC,
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
      IF (LHOOK) CALL DR_HOOK('LFIEDO_MT',1,ZHOOK_HANDLE)
      END

