! Oct-2012 P. Marguinaud 64b LFI
#include "lfisuffix.h"
! Jan-2011 P. Marguinaud Thread-safe FA
SUBROUTINE FANDAT_MT_BB                          &
&                     (FA,  KREP, KNUMER, KDATEF )
USE FA_MOD, ONLY : FA_COM
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
!****
!      Sous-programme de definition d'une (Nouvelle) Date sur un fichier
!     ARpege (sans ecriture de l'article correspondant)
!**
!    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
!                KNUMER (Entree) ==> Numero de l'unite logique;
!     (Tableau)  KDATEF (Entree) ==> Date elle-meme (FA%JPLDAT mots).
!*
!
!
!
TYPE(FA_COM) :: FA
INTEGER (KIND=JPLIKB) KREP, KNUMER
INTEGER (KIND=JPLIKB) KDATEF (FA%JPLDAT)
!
INTEGER (KIND=JPLIKB) IRANG, IREP, INIMES, J
!
LOGICAL LLVERF, LLRLFI, LLMODA
!
CHARACTER(LEN=FA%JPXNOM) CLACTI 
CHARACTER(LEN=FA%JPLMES) CLMESS 
CHARACTER(LEN=FA%JPLSPX) CLNSPR
LOGICAL                  LLFATA

!**
!     1.  -  CONTROLES DES PARAMETRES D'APPEL, PUIS INITIALISATIONS.
!-----------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('FANDAT_MT',0,ZHOOK_HANDLE)
CLACTI=''
LLVERF=.FALSE.
LLRLFI=.FALSE.
LLMODA=.FALSE.
CALL FANUMU_MT_BB                &
&               (FA, KNUMER,IRANG)
!
IF (IRANG.EQ.0) THEN
  IREP=-51
  GOTO 1001
ENDIF
!
!         Verrouillage eventuel du fichier.
!
IF (FA%LFAMUL) CALL LFIVER_MT_BB                             &
&                              (FA%LFI, FA%FICHIER(IRANG)%VRFICH,'ON')
LLVERF=FA%LFAMUL
!**
!     2.  -  DEFINITION PROPREMENT DITE VIA LE SOUS-PROGRAMME "FANDAI".
!            ( controles, puis mise a jour de FA%MADATE(.,IRANG) )
!-----------------------------------------------------------------------
!
CALL FANDAI_MT_BB                            &
&               (FA, IREP,IRANG,KDATEF,LLMODA)
!
IF (IREP.EQ.0) THEN
  LLRLFI=.FALSE.
  FA%FICHIER(IRANG)%LCREAF=.FALSE.
ENDIF
!**
!    10.  -  PHASE TERMINALE : MESSAGERIE, AVEC "ABORT" EVENTUEL,
!            VIA LE SOUS-PROGRAMME "FAIPAR" .
!-----------------------------------------------------------------------
!
1001 CONTINUE
KREP=IREP
LLFATA=LLMOER (IREP,IRANG)
!
!        Deverrouillage eventuel du fichier.
!
IF (LLVERF) CALL LFIVER_MT_BB                              &
&                           (FA%LFI, FA%FICHIER(IRANG)%VRFICH,'OFF')
!
IF (LLFATA) THEN
  INIMES=2
ELSEIF (IREP.NE.0) THEN
  INIMES=0
ELSE
  INIMES=IXNVMS(IRANG)
ENDIF
!
IF (.NOT.LLFATA.AND.INIMES.EQ.0)  THEN 
  IF (LHOOK) CALL DR_HOOK('FANDAT_MT',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!
CLNSPR='FANDAT'
!
IF (INIMES.GE.1.AND.LLMODA) THEN
  WRITE (UNIT=CLMESS,FMT=                                  &
&         '(''MODIFICATION DE LA DATE, UNITE'',I3)') KNUMER
  CALL FAIPAR_MT_BB                                     &
&                 (FA, KNUMER,INIMES,IREP,.FALSE.,CLMESS, &
&                  CLNSPR,CLACTI,.FALSE.)
ENDIF
!
IF (INIMES.EQ.2) THEN
!***** FAZZZZ - KREP=iiii, KNUMER=iii, KDATEF(1:5)=iiiii/ii/ii iii:ii, *****
!*****          KDATEF(7:8)=iiiiii-iiiiii                              *****
  WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3, &
&       '', KDATEF(1:5)='',I5,2(''/'',I2),I3,'':'',I2.2,   &
&       '', KDATEF(7:8)='',I6,''-'',I6)') KREP,KNUMER,     &
&     (KDATEF(J),J=1,5),(KDATEF(J),J=7,8)
  CALL FAIPAR_MT_BB                                    &
&                 (FA, KNUMER,INIMES,IREP,LLFATA,CLMESS, &
&               CLNSPR,CLACTI,LLRLFI)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('FANDAT_MT',1,ZHOOK_HANDLE)

CONTAINS

#include "facom2.llmoer.h"
#include "facom2.ixnvms.h"

END SUBROUTINE



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE FANDAT_BB             &
&           (KREP, KNUMER, KDATEF)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKB)  KREP                                   !   OUT
INTEGER (KIND=JPLIKB)  KNUMER                                 ! IN   
INTEGER (KIND=JPLIKB)  KDATEF     (*)                 ! IN   

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FANDAT_MT_BB                    &
&           (FA, KREP, KNUMER, KDATEF)

END SUBROUTINE

SUBROUTINE FANDAT_MM             &
&           (KREP, KNUMER, KDATEF)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KNUMER                                 ! IN   
INTEGER (KIND=JPLIKM)  KDATEF     (*)                 ! IN   

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FANDAT_MT_MM                    &
&           (FA, KREP, KNUMER, KDATEF)

END SUBROUTINE

SUBROUTINE FANDAT_MT_MM              &
&           (FA, KREP, KNUMER, KDATEF)
USE FA_MOD, ONLY : FA_COM,              &
&                   FA_COM_DEFAULT_INIT, &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
TYPE (FA_COM)          FA                                     ! INOUT
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KNUMER                                 ! IN   
INTEGER (KIND=JPLIKM)  KDATEF     (FA%JPLDAT)                 ! IN   
! Local integers
INTEGER (KIND=JPLIKB)  IREP                                   !   OUT
INTEGER (KIND=JPLIKB)  INUMER                                 ! IN   
INTEGER (KIND=JPLIKB)  IDATEF     (FA%JPLDAT)                 ! IN   
! Convert arguments

INUMER     = INT (    KNUMER, JPLIKB)
IDATEF     = INT (    KDATEF, JPLIKB)

CALL FANDAT_MT_BB                    &
&           (FA, IREP, INUMER, IDATEF)

KREP       = INT (      IREP, JPLIKM)

END SUBROUTINE

!INTF KREP            OUT                               
!INTF KNUMER        IN                                  
!INTF KDATEF        IN    DIMS=FA%JPLDAT                
