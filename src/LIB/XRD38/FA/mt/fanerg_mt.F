C Jan-2011 P. Marguinaud Thread-safe FA
      SUBROUTINE FANERG_MT (FA,  KNIVAU )
      USE FA_MOD, ONLY : FA_COM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        Ce sous-programme se charge de mettre le Niveau Global d'Erreur
C     Fatale du logiciel de Fichiers ARPEGE (*FA%NRFAGA*) a la valeur
C     KNIVAU, de meme que la variable correspondante du logiciel LFI.
C     Par defaut, FA%NRFAGA vaut 1.
C**
C        Argument : KNIVAU (Entree) ==> Niveau global d'erreur fatale.
C                                       Valeurs possibles:
C
C     0 : Rendre fatale toute erreur detectee, meme si elle correspond
C         a un fichier ouvert avec l'option "pas d'erreur fatale".
C     1 : Ne rend fatales que certaines erreurs "globales", c'est-a-dire
C         non reliables a un fichier ouvert, et les erreurs sur les fi-
C         chiers ouverts avec l'option "erreur fatale" (Mode par defaut)
C     2 : Passer outre toute erreur detectee, meme si elle correspond
C         a un fichier ouvert avec l'option "erreur fatale".
C         Neanmoins, le code-reponse eventuel ne sera pas zero.
C
#include "precision.h"
C
C
      TYPE(FA_COM) :: FA
      INTEGER KNIVAU
C
      INTEGER IREP, INIMES, INUMER
#include "facom_mt.h"
C
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FANERG_MT',0,ZHOOK_HANDLE)
      IF (FA%FANERG_LLPREA) THEN
        CALL FARINE_MT (FA, 2)
        FA%FANERG_LLPREA=.FALSE.
      ENDIF
C
      IF (KNIVAU.GE.0.AND.KNIVAU.LE.2) THEN
        FA%NRFAGA=KNIVAU
        CALL LFINEG_MT (FA%LFI, KNIVAU)
        IREP=0
      ELSE
        IREP=-52
      ENDIF
C
      LLFATA=IREP.NE.0.AND.FA%NRFAGA.NE.2
C
      IF (LLFATA) THEN
        INIMES=2
      ELSEIF (IREP.NE.0) THEN
        INIMES=0
      ELSEIF (FA%NIMSGA.EQ.2) THEN
        INIMES=2
      ELSE
        IF (LHOOK) CALL DR_HOOK('FANERG_MT',1,ZHOOK_HANDLE)
        RETURN
      ENDIF
C
      INUMER=FA%JPNIIL
      CLNSPR='FANERG'
C
      IF (MAX (INIMES,FA%NIMSGA).EQ.2) THEN
        WRITE (UNIT=CLMESS,
     S         FMT='(''KNIVAU='',I5,'', CODE INTERNE='',I4)'
     S         ) KNIVAU,IREP
      ENDIF
C
      CALL FAIPAR_MT (FA, INUMER,INIMES,IREP,LLFATA,CLMESS,
     S             CLNSPR,CLACTI,.FALSE.)
C
      IF (LHOOK) CALL DR_HOOK('FANERG_MT',1,ZHOOK_HANDLE)
      END

