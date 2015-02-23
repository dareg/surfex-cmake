! Oct-2012 P. Marguinaud 64b LFI
#include "lfisuffix.h"
! Jan-2011 P. Marguinaud Thread-safe LFI
SUBROUTINE LFIOEF_MT_BB                          &
&                     (LFI, KREP, KNUMER, LDERFA )
USE LFIMOD, ONLY : LFICOM
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
!****
!        Ce sous-programme permet d'obtenir l'option courante
!     de traitement des ERREURS detectees par le logiciel LFI,
!     pour un fichier particulier ouvert pour ce logiciel.
!
!        Noter que si le niveau global d'erreur fatale *LFI%NERFAG*
!     vaut 0 ou 2, l'option propre au fichier est inoperante.
!     *LFI%NERFAG* vaut par defaut 1, est reglable via le s/p "LFINEG",
!     et sa valeur peut etre obtenue par le s/p "LFIOEG".
!**
!     ARGUMENTS : KREP   (Sortie) ==> Code-REPonse du sous-programme;
!                 KNUMER (Entree) ==> NUMERo d'unite logique concernee;
!                 LDERFA (Sortie) ==> Option d'ERreur FAtale (vrai=oui).
!
!
TYPE(LFICOM) :: LFI
INTEGER (KIND=JPLIKB) KREP, KNUMER, IRANG, IREP, INIMES
!
LOGICAL LDERFA
!
CHARACTER(LEN=LFI%JPLSPX) CLNSPR
CHARACTER(LEN=LFI%JPLMES) CLMESS
CHARACTER(LEN=LFI%JPLFTX) CLACTI
LOGICAL LLFATA

!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('LFIOEF_MT',0,ZHOOK_HANDLE)
CLACTI=''
CALL LFINUM_MT_BB                 &
&               (LFI, KNUMER,IRANG)
!
IF (IRANG.EQ.0) THEN
  IREP=-1
ELSE
  LDERFA=LFI%LERFAT(IRANG)
  LFI%NDEROP(IRANG)=19
  IREP=0
ENDIF
!
LLFATA=LLMOER (IREP,IRANG)
KREP=IREP
!
IF (LLFATA.OR.IXNIMS (IRANG).EQ.2) THEN
  INIMES=2
ELSE
  IF (LHOOK) CALL DR_HOOK('LFIOEF_MT',1,ZHOOK_HANDLE)
  RETURN
ENDIF
!
CLNSPR='LFIOEF'
WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KNUMER='',I3, &
&       '', LDERFA= '',L1)') KREP,KNUMER,LDERFA
CALL LFIEMS_MT_BB                              &
&               (LFI, KNUMER,INIMES,IREP,LLFATA, &
&                CLMESS,CLNSPR,CLACTI)
!
IF (LHOOK) CALL DR_HOOK('LFIOEF_MT',1,ZHOOK_HANDLE)

CONTAINS

#include "lficom2.ixnims.h"
#include "lficom2.llmoer.h"

END SUBROUTINE



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE LFIOEF_BB             &
&           (KREP, KNUMER, LDERFA)
USE LFIMOD, ONLY : LFI => LFICOM_DEFAULT, &
&                   LFICOM_DEFAULT_INIT,   &
&                   NEW_LFI_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKB)  KREP                                   !   OUT
INTEGER (KIND=JPLIKB)  KNUMER                                 ! IN   
LOGICAL                LDERFA                                 !   OUT

IF (.NOT. LFICOM_DEFAULT_INIT) CALL NEW_LFI_DEFAULT ()

CALL LFIOEF_MT_BB                     &
&           (LFI, KREP, KNUMER, LDERFA)

END SUBROUTINE

SUBROUTINE LFIOEF_MM             &
&           (KREP, KNUMER, LDERFA)
USE LFIMOD, ONLY : LFI => LFICOM_DEFAULT, &
&                   LFICOM_DEFAULT_INIT,   &
&                   NEW_LFI_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KNUMER                                 ! IN   
LOGICAL                LDERFA                                 !   OUT

IF (.NOT. LFICOM_DEFAULT_INIT) CALL NEW_LFI_DEFAULT ()

CALL LFIOEF_MT_MM                     &
&           (LFI, KREP, KNUMER, LDERFA)

END SUBROUTINE

SUBROUTINE LFIOEF_MT_MM               &
&           (LFI, KREP, KNUMER, LDERFA)
USE LFIMOD, ONLY : LFICOM,              &
&                   LFICOM_DEFAULT_INIT, &
&                   NEW_LFI_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
TYPE (LFICOM)          LFI                                    ! INOUT
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KNUMER                                 ! IN   
LOGICAL                LDERFA                                 !   OUT
! Local integers
INTEGER (KIND=JPLIKB)  IREP                                   !   OUT
INTEGER (KIND=JPLIKB)  INUMER                                 ! IN   
! Convert arguments

INUMER     = INT (    KNUMER, JPLIKB)

CALL LFIOEF_MT_BB                     &
&           (LFI, IREP, INUMER, LDERFA)

KREP       = INT (      IREP, JPLIKM)

END SUBROUTINE

!INTF KREP            OUT 
!INTF KNUMER        IN    
!INTF LDERFA          OUT 
