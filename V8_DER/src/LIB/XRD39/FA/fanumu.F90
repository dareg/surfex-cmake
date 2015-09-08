! Oct-2012 P. Marguinaud 64b LFI
#include "lfisuffix.h"
! Jan-2011 P. Marguinaud Thread-safe FA
SUBROUTINE FANUMU_MT_BB                   &
&                     (FA,  KNUMER, KRANG )
USE FA_MOD, ONLY : FA_COM
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
!****
!        Ce sous-programme calcule le RANG du Numero d'Unite logique
!     *KNUMER* dans la table des unites logiques *FA%NULOGI*;
!     s'il n'y est pas trouve, le resultat est ZERO.
!        Ce sous-programme, appele par tous les sous-programmes
!     "fichier" du Logiciel de Fichiers ARPEGE, se charge lors de son
!     premier appel d'appeler le sous-programme preparatoire FARINE.
!**
!       Arguments : KNUMER (Entree) ==> Numero d'unite logique cherche;
!                   KRANG  (Sortie) ==> Rang dans la table des fichiers
!                                       du logiciel FA (0 si absent).
!
!
!
TYPE(FA_COM) :: FA
INTEGER (KIND=JPLIKB) KNUMER, KRANG
!
INTEGER (KIND=JPLIKB) J, IRESUL
!
!

!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('FANUMU_MT',0,ZHOOK_HANDLE)

IF (FA%FANUMU_LLPREA) THEN
  CALL FARINE_MT_BB             &
&                 (FA, 2_JPLIKB )
  FA%FANUMU_LLPREA=.FALSE.
ENDIF
!
!          VERROUILLAGE GLOBAL (A CAUSE DE L'UTILISATION DE FA%NFIOUV )
!
IF (FA%LFAMUL) CALL LFIVER_MT_BB                      &
&                              (FA%LFI, FA%VRGLAS,'ON')
!
DO J=1,FA%NFIOUV
!
IF (KNUMER.EQ.FA%FICHIER(FA%NULIND(J))%NULOGI) THEN
  IRESUL=FA%NULIND(J)
  GOTO 20
ENDIF
!
ENDDO
!
IRESUL=0
!
20 CONTINUE
!
!          DEVERROUILLAGE GLOBAL
!
IF (FA%LFAMUL) CALL LFIVER_MT_BB                       &
&                              (FA%LFI, FA%VRGLAS,'OFF')
KRANG=IRESUL
!
IF (LHOOK) CALL DR_HOOK('FANUMU_MT',1,ZHOOK_HANDLE)
END SUBROUTINE



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE FANUMU_BB          &
&           (KNUMER, KRANG)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKB)  KNUMER                                 ! IN   
INTEGER (KIND=JPLIKB)  KRANG                                  !   OUT

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FANUMU_MT_BB             &
&           (FA, KNUMER, KRANG)

END SUBROUTINE

SUBROUTINE FANUMU_MM          &
&           (KNUMER, KRANG)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKM)  KNUMER                                 ! IN   
INTEGER (KIND=JPLIKM)  KRANG                                  !   OUT

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FANUMU_MT_MM             &
&           (FA, KNUMER, KRANG)

END SUBROUTINE

SUBROUTINE FANUMU_MT_MM          &
&           (FA, KNUMER, KRANG)
USE FA_MOD, ONLY : FA_COM,              &
&                   FA_COM_DEFAULT_INIT, &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
TYPE (FA_COM)          FA                                     ! INOUT
INTEGER (KIND=JPLIKM)  KNUMER                                 ! IN   
INTEGER (KIND=JPLIKM)  KRANG                                  !   OUT
! Local integers
INTEGER (KIND=JPLIKB)  INUMER                                 ! IN   
INTEGER (KIND=JPLIKB)  IRANG                                  !   OUT
! Convert arguments

INUMER     = INT (    KNUMER, JPLIKB)

CALL FANUMU_MT_BB             &
&           (FA, INUMER, IRANG)

KRANG      = INT (     IRANG, JPLIKM)

END SUBROUTINE

!INTF KNUMER        IN    
!INTF KRANG           OUT 
