! Oct-2012 P. Marguinaud 64b LFI
#include "lfisuffix.h"
! Jan-2011 P. Marguinaud Thread-safe FA
SUBROUTINE FANDAI_MT_BB                                 &
&                     (FA,  KREP, KRANG, KDATEF, LDMODA )
USE FA_MOD, ONLY : FA_COM, JPNIIL
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
!****
!      Sous-programme INTERNE du logiciel de Fichiers ARPEGE:
!     Definition d'une (Nouvelle) Date.
!**
!    Arguments : KREP   (Sortie) ==> Code-reponse du sous-programme;
!                KRANG  (Entree) ==> Rang de l'unite logique;
!     (Tableau)  KDATEF (Entree) ==> Date elle-meme (FA%JPLDAT mots).
!                LDMODA (Sortie) ==> Vrai s'il y a modification d'une
!                                    date deja definie.
!*
!       En mode multi-taches, il doit y avoir verrouillage du fichier
!     concerne avant l'appel au sous-programme.
!
!
!
TYPE(FA_COM) :: FA
INTEGER (KIND=JPLIKB) KREP, KRANG
INTEGER (KIND=JPLIKB) KDATEF (FA%JPLDAT)
!
INTEGER (KIND=JPLIKB) IMI123, IMAX69, IMINIM
INTEGER (KIND=JPLIKB) J, ILMOIS, IDEBUT, INIMES, INUMER
!
LOGICAL LDMODA
!
CHARACTER(LEN=FA%JPXNOM) CLACTI 
CHARACTER(LEN=FA%JPLMES) CLMESS 
CHARACTER(LEN=FA%JPLSPX) CLNSPR
LOGICAL                  LLFATA

!**
!     1.  -  CONTROLES DES PARAMETRES D'APPEL.
!-----------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('FANDAI_MT',0,ZHOOK_HANDLE)
CLACTI=''
LDMODA=.FALSE.
!
IF (KRANG.LE.0.OR.KRANG.GT.FA%JPNXFA) THEN
  KREP=-66
  GOTO 1001
ENDIF
!
!         Controle de la Date proprement dite.
!
IMI123=MIN (KDATEF(1),KDATEF(2),KDATEF(3))
IMAX69=MAX (KDATEF(6),KDATEF(9))
IMINIM=KDATEF(1)
!
DO J=2,FA%JPLDAT
IMINIM=MIN (IMINIM,KDATEF(J))
ENDDO
!
IF (IMINIM.LT.0.OR.IMI123.LE.0.OR.KDATEF(2).GT.12.OR.              &
&    KDATEF(3).GT.31.OR.KDATEF(4).GE.24.OR.KDATEF(5).GE.60.OR.      &
&    IMAX69.GE.255.OR.                                              &
& (KDATEF(10).LE.KDATEF(11).AND.(KDATEF(10)*KDATEF(11)).NE.0)) THEN
!
!        Erreur de syntaxe.
!
  KREP=-82
  GOTO 1001
ELSEIF ((KDATEF(2).GT.7.OR.MOD (KDATEF(2),2_JPLIKB ).EQ.0).AND.  &
&        (KDATEF(2).LE.7.OR.MOD (KDATEF(2),2_JPLIKB ).EQ.1)) THEN
!
!        Controle de coherence (annee,mois,jour).
!
  IF (KDATEF(2).EQ.2) THEN
    ILMOIS=28+MAX (0_JPLIKB ,1-MOD (KDATEF(1),4_JPLIKB ))
  ELSE
    ILMOIS=30
  ENDIF
!
  IF (KDATEF(3).GT.ILMOIS) THEN
    KREP=-82
    GOTO 1001
  ENDIF
!
ENDIF
!
KREP=0
!**
!     2.  -  SI DATE DEJA DEFINIE, COMPARAISON ANCIENNE/NOUVELLE.
!-----------------------------------------------------------------------
!
IF (FA%FICHIER(KRANG)%LCREAF) THEN
  IDEBUT=1
ELSE
!
  DO J=1,FA%JPLDAT
!
  IF (FA%FICHIER(KRANG)%MADATE(J).NE.KDATEF(J)) THEN
    LDMODA=.TRUE.
    IDEBUT=J
    GOTO 300
  ENDIF
!
  ENDDO
!
!         Si on arrive ici, il y a redefinition a l'identique.
!
  GOTO 1001
ENDIF
!**
!     3.  -  SI NECESSAIRE, MISE A JOUR DU TABLEAU "FA%MADATE".
!-----------------------------------------------------------------------
!
300 CONTINUE
!
DO J=IDEBUT,FA%JPLDAT
FA%FICHIER(KRANG)%MADATE(J)=KDATEF(J)
ENDDO
!**
!    10.  -  PHASE TERMINALE : MESSAGERIE EVENTUELLE,
!            VIA LE SOUS-PROGRAMME "FAIPAR" .
!-----------------------------------------------------------------------
!
1001 CONTINUE
LLFATA=LLMOER (KREP,KRANG)
!
IF (FA%LFAMOP.OR.LLFATA) THEN
  INIMES=2
  CLNSPR='FANDAI'
  INUMER=JPNIIL
  WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KRANG='',I4,   &
&       '', KDATEF(1:5)='',I5,2(''/'',I2),I3,'':'',I2.2,    &
&       '', KDATEF(7:8)='',I6,''-'',I6,'', LDMODA= '',L1)') &
&     KREP,KRANG,(KDATEF(J),J=1,5),(KDATEF(J),J=7,8),LDMODA
  CALL FAIPAR_MT_BB                                     &
&                 (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS, &
&                  CLNSPR,CLACTI, .FALSE.)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('FANDAI_MT',1,ZHOOK_HANDLE)

CONTAINS

#include "facom2.llmoer.h"

END SUBROUTINE



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE FANDAI_BB                    &
&           (KREP, KRANG, KDATEF, LDMODA)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKB)  KREP                                   !   OUT
INTEGER (KIND=JPLIKB)  KRANG                                  ! IN   
INTEGER (KIND=JPLIKB)  KDATEF     (*)                 ! IN   
LOGICAL                LDMODA                                 !   OUT

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FANDAI_MT_BB                           &
&           (FA, KREP, KRANG, KDATEF, LDMODA)

END SUBROUTINE

SUBROUTINE FANDAI_MM                    &
&           (KREP, KRANG, KDATEF, LDMODA)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KRANG                                  ! IN   
INTEGER (KIND=JPLIKM)  KDATEF     (*)                 ! IN   
LOGICAL                LDMODA                                 !   OUT

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FANDAI_MT_MM                           &
&           (FA, KREP, KRANG, KDATEF, LDMODA)

END SUBROUTINE

SUBROUTINE FANDAI_MT_MM                     &
&           (FA, KREP, KRANG, KDATEF, LDMODA)
USE FA_MOD, ONLY : FA_COM,              &
&                   FA_COM_DEFAULT_INIT, &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
TYPE (FA_COM)          FA                                     ! INOUT
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KRANG                                  ! IN   
INTEGER (KIND=JPLIKM)  KDATEF     (FA%JPLDAT)                 ! IN   
LOGICAL                LDMODA                                 !   OUT
! Local integers
INTEGER (KIND=JPLIKB)  IREP                                   !   OUT
INTEGER (KIND=JPLIKB)  IRANG                                  ! IN   
INTEGER (KIND=JPLIKB)  IDATEF     (FA%JPLDAT)                 ! IN   
! Convert arguments

IRANG      = INT (     KRANG, JPLIKB)
IDATEF     = INT (    KDATEF, JPLIKB)

CALL FANDAI_MT_BB                           &
&           (FA, IREP, IRANG, IDATEF, LDMODA)

KREP       = INT (      IREP, JPLIKM)

END SUBROUTINE

!INTF KREP            OUT                               
!INTF KRANG         IN                                  
!INTF KDATEF        IN    DIMS=FA%JPLDAT                
!INTF LDMODA          OUT                               
