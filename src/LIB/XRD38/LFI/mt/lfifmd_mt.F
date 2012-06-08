C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFIFMD_MT (LFI, KFACMD )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Ce sous-programme permet de changer le Facteur Multiplicatif
C     par Defaut du logiciel de fichiers indexes LFI.
C        Apres appel reussi a ce sous-programme, toute ouverture d'unite
C     logique LFI pour laquelle il n'y a pas de facteur multiplicatif
C     predefini (via *LFIAFM*) se fera en traitant le fichier avec une
C     longueur PHYSIQUE d'article de LFI%JPLARD*KFACMD mots.
C
C        La valeur implicite de ce Facteur Multiplicatif par Defaut est
C     definie dans *LFIINI* ( en l'occurrence, il s'agit de 1 ) .
C**
C       ARGUMENT  : KFACMD (Entree) ==> Facteur Multiplicatif par Defaut
C
#ifndef f77
#include "precision.h"
#endif
C
      TYPE(LFICOM) :: LFI
      INTEGER KFACMD, INIMES, IREP, INUMER
#include "lficom_mt.h"
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFIFMD_MT',0,ZHOOK_HANDLE)
      IF (LFI%LFIFMD_LLPREA) THEN
        CALL LFIINI_MT (LFI, 2)
        LFI%LFIFMD_LLPREA=.FALSE.
      ENDIF
C
      IF (KFACMD.LE.0) THEN
        IREP=-14
      ELSEIF (KFACMD.GT.LFI%JPFACX) THEN
        IREP=-28
      ELSE
        IREP=0
C
C          Modification, sous Verrouillage Global eventuel.
C
        IF (LFI%LMULTI) CALL LFIVER_MT (LFI, LFI%VERGLA,'ON')
C
        LFI%MFACTU(0)=KFACMD
C
        IF (LFI%LMULTI) CALL LFIVER_MT (LFI, LFI%VERGLA,'OFF')
      ENDIF
C
C        MESSAGERIE EVENTUELLE, AVEC ABORT SI NECESSAIRE .
C
      LLFATA=IREP.NE.0.AND.LFI%NERFAG.NE.2
C
      IF (LLFATA) THEN
        INIMES=2
      ELSEIF (IREP.NE.0) THEN
        INIMES=0
      ELSEIF (LFI%NIMESG.EQ.0) THEN
        IF (LHOOK) CALL DR_HOOK('LFIFMD_MT',1,ZHOOK_HANDLE)
        RETURN
      ELSE
        INIMES=LFI%NIMESG
      ENDIF
C
      CLNSPR='LFIFMD'
      INUMER=LFI%JPNIL
C
      IF (INIMES.EQ.2) THEN
C
        IF (LFI%LFRANC) THEN
          WRITE (UNIT=CLMESS,
     S           FMT='(''KFACMD='',I5,'', CODE INTERNE='',
     S           I4)') KFACMD,IREP
        ELSE
          WRITE (UNIT=CLMESS,
     S           FMT='(''KFACMD='',I5,'', INTERNAL CODE='',
     S           I4)') KFACMD,IREP
        ENDIF
C
        CALL LFIEMS_MT (LFI, INUMER,INIMES,IREP,LLFATA,
     S                  CLMESS,CLNSPR,CLACTI)
      ENDIF
C
      IF (INIMES.GE.1) THEN
C
        IF (LFI%LFRANC) THEN
          WRITE (UNIT=CLMESS,FMT=
     S         '(''Reglage du Facteur Multiplicatif par Defaut a'',I3)')
     S         KFACMD
        ELSE
          WRITE (UNIT=CLMESS,FMT=
     S         '(''Default Multiply Factor set to'',I3)')
     S         KFACMD
        ENDIF
C
      ENDIF
C
      CALL LFIEMS_MT (LFI, INUMER,INIMES,IREP,LLFATA,
     S                CLMESS,CLNSPR,CLACTI)
C
      IF (LHOOK) CALL DR_HOOK('LFIFMD_MT',1,ZHOOK_HANDLE)
      END

