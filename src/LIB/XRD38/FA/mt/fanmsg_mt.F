C Jan-2011 P. Marguinaud Thread-safe FA
C R. El Khatib 30-Mar-2012 KULOUT
      SUBROUTINE FANMSG_MT (FA,  KNIVAU, KULOUT )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Ce sous-programme se charge de mettre le Niveau Global d'
C     impression des Messages du logiciel de Fichiers ARPEGE (*FA%NIMSGA*)
C     a la valeur KNIVAU, de meme que la variable correspondante du
C     du logiciel LFI. Par defaut, FA%NIMSGA vaut 1.
C**
C        Argument : KNIVAU (Entree) ==> Niveau Global d'Impression
C                                       des Messages.
C                                       Valeurs possibles:
C
C     0 : N'emettre que les messages d'erreurs reellement importants .
C     1 : N'emettre qu'un minimum de messages "globaux", et les messages
C         lies a un fichier ouvert qui sont de niveau au plus egal au
C         niveau de la messagerie pour ce fichier (Mode par defaut) .
C     2 : Emettre tous les messages possibles, meme s'ils ne correspon-
C         dent pas a un fichier ouvert avec le niveau de Messagerie 2 .
C
C                   KULOUT : logical unit number for printing
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KNIVAU
      INTEGER KULOUT
C
      INTEGER IREP, INIMES, INUMER
#include "facom_mt.h"
C
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE

      IF (LHOOK) CALL DR_HOOK('FANMSG_MT',0,ZHOOK_HANDLE)

      IF (FA%FANMSG_LLPREA) THEN
        CALL FARINE_MT (FA, 2)
        FA%FANMSG_LLPREA=.FALSE.
      ENDIF
C
      IF (KNIVAU.GE.0.AND.KNIVAU.LE.2) THEN
        INIMES=MAX (FA%NIMSGA,KNIVAU)
        FA%NIMSGA=KNIVAU
        CALL LFINMG_MT (FA%LFI, KNIVAU,KULOUT)
        IREP=0
      ELSE
        INIMES=FA%NIMSGA
        IREP=-52
      ENDIF
C
      LLFATA=IREP.NE.0.AND.FA%NRFAGA.NE.2
C
      IF (LLFATA) THEN
        INIMES=2
      ELSEIF (IREP.NE.0) THEN
        INIMES=0
      ELSEIF (INIMES.NE.2) THEN
        IF (LHOOK) CALL DR_HOOK('FANMSG_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      INUMER=FA%JPNIIL
      CLNSPR='FANMSG'
C
      IF (MAX (INIMES,FA%NIMSGA).EQ.2) THEN
        WRITE (UNIT=CLMESS,
     S         FMT='(''KNIVAU='',I5,'', CODE INTERNE='',I4)'
     S         ) KNIVAU,IREP
        IF (INIMES.NE.2) CALL FAIPAR_MT (FA, INUMER,FA%NIMSGA,IREP,
     S                                .FALSE.,CLMESS,
     S                                CLNSPR,CLACTI,.FALSE.)
      ENDIF
C
      CALL FAIPAR_MT (FA, INUMER,INIMES,IREP,LLFATA,CLMESS,
     S             CLNSPR,CLACTI,
     S             .FALSE.)
C
      IF (LHOOK) CALL DR_HOOK('FANMSG_MT',1,ZHOOK_HANDLE)
      END
