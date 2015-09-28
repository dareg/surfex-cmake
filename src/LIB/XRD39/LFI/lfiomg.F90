! Oct-2012 P. Marguinaud 64b LFI
#include "lfisuffix.h"
! Jan-2011 P. Marguinaud Thread-safe LFI
SUBROUTINE LFIOMG_MT_BB            &
&                     (LFI, KNIVAU )
USE LFIMOD, ONLY : LFICOM
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
!****
!        Ce sous-programme permet d'obtenir la valeur courante du NIVEAU
!      GLOBAL d'IMPRESSION des MESSAGES emis par le logiciel LFI.
!**
!        ARGUMENT : KNIVAU (Sortie) ==> Niveau global de messagerie.
!
!     Valeurs possibles (par defaut, 1):
!
!     0 : N'emettre que les messages d'erreurs reellement importants,
!         en pratique ceux relies a une erreur fatale. Le niveau indivi-
!         duel de messagerie des fichiers est alors inoperant.
!
!     1 : Seuls quelques messages "globaux" (c'est-a-dire non reliables
!         a un fichier deja ouvert) et les messages lies a un fichier
!         ouvert, messages dont le niveau est au plus egal au niveau de
!         messagerie individuelle du fichier.
!
!     2 : Emettre tous les messages possibles (donc jusqu'au niveau 2,
!         mais pas ceux lies au mode "mise au point"), meme si
!         ces messages concernent un fichier dont le niveau individuel
!         de messagerie est inferieur a 2.
!
!
TYPE(LFICOM) :: LFI
INTEGER (KIND=JPLIKB) KNIVAU, INIMES, IREP, INUMER
CHARACTER(LEN=LFI%JPLSPX) CLNSPR
CHARACTER(LEN=LFI%JPLMES) CLMESS
CHARACTER(LEN=LFI%JPLFTX) CLACTI

!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('LFIOMG_MT',0,ZHOOK_HANDLE)
CLACTI=''
IF (LFI%LFIOMG_LLPREA) THEN
  CALL LFIINI_MT_BB              &
&                 (LFI, 2_JPLIKB )
  LFI%LFIOMG_LLPREA=.FALSE.
ENDIF
!
KNIVAU=LFI%NIMESG
INIMES=LFI%NIMESG
!
IF (INIMES.EQ.2) THEN
  IREP=0
  INUMER=LFI%JPNIL
  CLNSPR='LFIOMG'
  WRITE (UNIT=CLMESS,FMT='(''KNIVAU='',I2)') KNIVAU
  CALL LFIEMS_MT_BB                                      &
&                 (LFI, INUMER,INIMES,IREP,.FALSE.,CLMESS, &
&                  CLNSPR,CLACTI)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('LFIOMG_MT',1,ZHOOK_HANDLE)
END SUBROUTINE



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE LFIOMG_BB          &
&           (KNIVAU)
USE LFIMOD, ONLY : LFI => LFICOM_DEFAULT, &
&                   LFICOM_DEFAULT_INIT,   &
&                   NEW_LFI_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKB)  KNIVAU                                 !   OUT

IF (.NOT. LFICOM_DEFAULT_INIT) CALL NEW_LFI_DEFAULT ()

CALL LFIOMG_MT_BB          &
&           (LFI, KNIVAU)

END SUBROUTINE

SUBROUTINE LFIOMG_MM          &
&           (KNIVAU)
USE LFIMOD, ONLY : LFI => LFICOM_DEFAULT, &
&                   LFICOM_DEFAULT_INIT,   &
&                   NEW_LFI_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKM)  KNIVAU                                 !   OUT

IF (.NOT. LFICOM_DEFAULT_INIT) CALL NEW_LFI_DEFAULT ()

CALL LFIOMG_MT_MM          &
&           (LFI, KNIVAU)

END SUBROUTINE

SUBROUTINE LFIOMG_MT_MM          &
&           (LFI, KNIVAU)
USE LFIMOD, ONLY : LFICOM,              &
&                   LFICOM_DEFAULT_INIT, &
&                   NEW_LFI_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
TYPE (LFICOM)          LFI                                    ! INOUT
INTEGER (KIND=JPLIKM)  KNIVAU                                 !   OUT
! Local integers
INTEGER (KIND=JPLIKB)  INIVAU                                 !   OUT
! Convert arguments


CALL LFIOMG_MT_BB          &
&           (LFI, INIVAU)

KNIVAU     = INT (    INIVAU, JPLIKM)

END SUBROUTINE

!INTF KNIVAU          OUT 
