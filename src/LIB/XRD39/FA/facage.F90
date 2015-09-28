! Oct-2012 P. Marguinaud 64b LFI
#include "lfisuffix.h"
! Jan-2011 P. Marguinaud Thread-safe FA
SUBROUTINE FACAGE_MT_BB                    &
&                     (FA,  CDNOMC, LDGARD )
USE FA_MOD, ONLY : FA_COM, JPNIIL
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
!****
!        Sous-programme servant a redefinir l'option de conservation
!     d'un cadre preexistant ( CAdre a Garder Eventuellement... )
!**
!        Arguments : CDNOMC ==> Nom symbolique du cadre;
!  (tous d'Entree)   LDGARD ==> Vrai si le cadre doit etre conserve meme
!                               apres la fermeture du dernier fichier
!                               qui s'y rattache.
!
!
!
TYPE(FA_COM) :: FA
INTEGER (KIND=JPLIKB) ILCDNO, ILNOMC, J, IRANGC
INTEGER (KIND=JPLIKB) IREP, INIMES, INUMER
!
LOGICAL LLVERG, LDGARD
!
CHARACTER CDNOMC*(*)
CHARACTER(LEN=FA%JPXNOM) CLACTI 
CHARACTER(LEN=FA%JPLMES) CLMESS 
CHARACTER(LEN=FA%JPLSPX) CLNSPR
LOGICAL                  LLFATA

!
!
!**
!     0.  -  SI PREMIERE UTILISATION, APPEL AU SOUS-PROGRAMME "FARINE".
!-----------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('FACAGE_MT',0,ZHOOK_HANDLE)
CLACTI=''
IF (FA%FACAGE_LLPREA) THEN
  CALL FARINE_MT_BB             &
&                 (FA, 2_JPLIKB )
  FA%FACAGE_LLPREA=.FALSE.
ENDIF
!**
!     1.  -  CONTROLE DE L'ARGUMENT "CDNOMC".
!-----------------------------------------------------------------------
!
LLVERG=.FALSE.
ILCDNO=INT (LEN (CDNOMC), JPLIKB)
ILNOMC=1
!
IF (ILCDNO.LE.0) THEN
  IREP=-65
  GOTO 1001
ELSEIF (CDNOMC.EQ.' ') THEN
  IREP=-68
  GOTO 1001
ENDIF
!
DO J=ILCDNO,1,-1
!
IF (CDNOMC(J:J).NE.' ') THEN
  ILNOMC=J
  GOTO 102
ENDIF
!
ENDDO
!
102 CONTINUE
!
IF (ILNOMC.GT.FA%NCPCAD) THEN
  IREP=-65
  GOTO 1001
ENDIF
!**
!     2.  -  RECHERCHE DU CADRE DANS LES TABLES.
!-----------------------------------------------------------------------
!
!             Verrouillage global prealable, si necessaire.
!
IF (FA%LFAMUL) CALL LFIVER_MT_BB                      &
&                              (FA%LFI, FA%VRGLAS,'ON')
LLVERG=FA%LFAMUL
!
CALL FANUCA_MT_BB                         &
&               (FA, CDNOMC,IRANGC,.FALSE.)
!
IF (IRANGC.EQ.0) THEN
  IREP=-51
  GOTO 1001
ENDIF
!**
!     3.  -  MISE A JOUR DU NIVEAU DE CONSERVATION.
!-----------------------------------------------------------------------
!
IF (LDGARD) THEN
  FA%CADRE(IRANGC)%NGARDE=2
ELSE
  FA%CADRE(IRANGC)%NGARDE=0
ENDIF
!
IREP=0
!**
!    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
!            VIA LE sous-programme "FAIPAR" .
!-----------------------------------------------------------------------
!
1001 CONTINUE
!
!          Deverrouillage global eventuel.
!
IF (LLVERG) CALL LFIVER_MT_BB                       &
&                           (FA%LFI, FA%VRGLAS,'OFF')
!
LLFATA=IREP.NE.0.AND.FA%NRFAGA.NE.2
!
IF (LLFATA.OR.FA%NIMSGA.EQ.2) THEN
  INIMES=2
  CLNSPR='FACAGE'
!
  IF (IREP.EQ.-65.AND.ILCDNO.LE.0) THEN
    ILNOMC=8
    CLACTI(1:ILNOMC)=FA%CHAINC(:ILNOMC)
  ELSE
    ILNOMC=MIN (INT (LEN (CLACTI), JPLIKB),ILNOMC)
    CLACTI(1:ILNOMC)=CDNOMC(1:ILNOMC)
  ENDIF
!
  ILNOMC=MIN (ILNOMC,FA%NCPCAD)
  WRITE (UNIT=CLMESS,                                 &
&         FMT='(''CDNOMC= '''''',A,'''''', LDGARD= '', &
&         L1,'', CODE INTERNE='',I4)')                 &
&         CLACTI(1:ILNOMC),LDGARD,IREP
  INUMER=JPNIIL
  CALL FAIPAR_MT_BB                                    &
&                 (FA, INUMER,INIMES,IREP,LLFATA,CLMESS, &
&               CLNSPR,CLACTI(1:ILNOMC),.FALSE.)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('FACAGE_MT',1,ZHOOK_HANDLE)
END SUBROUTINE



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE FACAGE_BB          &
&           (CDNOMC, LDGARD)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
CHARACTER (LEN=*)      CDNOMC                                 ! IN   
LOGICAL                LDGARD                                 ! IN   

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FACAGE_MT_BB              &
&           (FA, CDNOMC, LDGARD)

END SUBROUTINE

SUBROUTINE FACAGE_MM          &
&           (CDNOMC, LDGARD)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
CHARACTER (LEN=*)      CDNOMC                                 ! IN   
LOGICAL                LDGARD                                 ! IN   

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FACAGE_MT_MM              &
&           (FA, CDNOMC, LDGARD)

END SUBROUTINE

SUBROUTINE FACAGE_MT_MM          &
&           (FA, CDNOMC, LDGARD)
USE FA_MOD, ONLY : FA_COM,              &
&                   FA_COM_DEFAULT_INIT, &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
TYPE (FA_COM)          FA                                     ! INOUT
CHARACTER (LEN=*)      CDNOMC                                 ! IN   
LOGICAL                LDGARD                                 ! IN   
! Local integers
! Convert arguments


CALL FACAGE_MT_BB              &
&           (FA, CDNOMC, LDGARD)


END SUBROUTINE

!INTF CDNOMC        IN    
!INTF LDGARD        IN    
