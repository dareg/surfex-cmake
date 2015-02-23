! Oct-2012 P. Marguinaud 64b LFI
! Jan-2011 P. Marguinaud Thread-safe FA
SUBROUTINE FANDAI_MT64                                  &
&                     (FA,  KREP, KRANG, KDATEF, LDMODA )
USE FA_MOD, ONLY : FA_COM, JPNIIL, &
                 & JD_YEA, JD_MON, JD_DAY, & 
                 & JD_HOU, JD_MIN, JD_TUN, & 
                 & JD_THO,         JD_IAN, &
                 & JD_CU1, JD_CU2,         &
                 & JD_DEX,         JD_SEM, &
                 & JD_SET, JD_CE1, JD_CE2, &
                 & JD_TST, JD_FMT

USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION
IMPLICIT NONE
!****
!      Sous-programme INTERNE du logiciel de Fichiers ARPEGE:
!     Definition d'une (Nouvelle) Date.
!**
!    Arguments : KREP   (Sortie)        ==> Code-reponse du sous-programme;
!                KRANG  (Entree)        ==> Rang de l'unite logique;
!     (Tableau)  KDATEF (Entree/Sortie) ==> Date elle-meme (FA%JPLDAT mots).
!                LDMODA (Sortie)        ==> Vrai s'il y a modification d'une
!                                           date deja definie.
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
INTEGER (KIND=JPLIKB) J, ILMOIS, INIMES, INUMER
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
KREP=0
!
IF (KRANG.LE.0.OR.KRANG.GT.FA%JPNXFA) THEN
  KREP=-66
  GOTO 1001
ENDIF
!
!         Controle de la Date proprement dite.
!
IMI123=MIN (KDATEF(JD_YEA),KDATEF(JD_MON),KDATEF(JD_DAY))
IMAX69=MAX (KDATEF(JD_TUN),KDATEF(JD_IAN))
IMINIM=KDATEF(JD_YEA)
!
DO J=2,FA%JPLDAT
IMINIM=MIN (IMINIM,KDATEF(J))
ENDDO
!
IF (IMINIM.LT.0.OR.IMI123.LE.0.OR.KDATEF(JD_MON).GT.12.OR.                    &
&   KDATEF(JD_DAY).GT.31.OR.KDATEF(JD_HOU).GE.24.OR.KDATEF(JD_MIN).GE.60.OR.  &
&   IMAX69.GE.255.OR.                                                         &
& (KDATEF(JD_CU1).LE.KDATEF(JD_CU2).AND.(KDATEF(JD_CU1)*KDATEF(JD_CU2)).NE.0)) THEN
!
!        Erreur de syntaxe.
!
  KREP=-82
  GOTO 1001
ELSEIF ((KDATEF(JD_MON).GT.7.OR.MOD (KDATEF(JD_MON),2_JPLIKB ).EQ.0).AND.  &
&       (KDATEF(JD_MON).LE.7.OR.MOD (KDATEF(JD_MON),2_JPLIKB ).EQ.1)) THEN
!
!        Controle de coherence (annee,mois,jour).
!
  IF (KDATEF(JD_MON).EQ.2) THEN
    ILMOIS=28+MAX (0_JPLIKB ,1-MOD (KDATEF(JD_YEA),4_JPLIKB ))
  ELSE
    ILMOIS=30
  ENDIF
!
  IF (KDATEF(JD_DAY).GT.ILMOIS) THEN
    KREP=-82
    GOTO 1001
  ENDIF
!
ENDIF


!**
!     2.  -  SI DATE DEJA DEFINIE, COMPARAISON ANCIENNE/NOUVELLE.
!-----------------------------------------------------------------------
!
IF (.NOT. FA%FICHIER(KRANG)%LCREAF) THEN
!
  DO J=1,FA%JPLDAT
!
  IF (FA%FICHIER(KRANG)%MADATE(J).NE.KDATEF(J)) THEN
    LDMODA=.TRUE.
    GOTO 300
  ENDIF
!
  ENDDO
!
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
FA%FICHIER(KRANG)%MADATE (:) = KDATEF (:)
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
  WRITE (UNIT=CLMESS,FMT='(''KREP='',I4,'', KRANG='',I4,    &
&       '', KDATEF(1:5)='',I5,2(''/'',I2),I3,'':'',I2.2,    &
&       '', KDATEF(7:8)='',I6,''-'',I6,                     &
&       '', LDMODA= '',L1)') &
&     KREP,KRANG,(KDATEF(J),J=1,5),(KDATEF(J),J=7,8),       &
&     LDMODA
  CALL FAIPAR_MT64                                      &
&               (FA, INUMER,INIMES,KREP,.FALSE.,CLMESS, &
&                CLNSPR,CLACTI, .FALSE.)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('FANDAI_MT',1,ZHOOK_HANDLE)

CONTAINS

#include "facom2.llmoer.h"

END SUBROUTINE



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE FANDAI64                     &
&           (KREP, KRANG, KDATEF, LDMODA)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKB)  KREP                                   !   OUT
INTEGER (KIND=JPLIKB)  KRANG                                  ! IN   
INTEGER (KIND=JPLIKB)  KDATEF     (FA%JPLDAT)                 ! IN   
LOGICAL                LDMODA                                 !   OUT

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FANDAI_MT64                            &
&           (FA, KREP, KRANG, KDATEF, LDMODA)

END SUBROUTINE

SUBROUTINE FANDAI                       &
&           (KREP, KRANG, KDATEF, LDMODA)
USE FA_MOD, ONLY : FA => FA_COM_DEFAULT, &
&                   FA_COM_DEFAULT_INIT,  &
&                   NEW_FA_DEFAULT
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KRANG                                  ! IN   
INTEGER (KIND=JPLIKM)  KDATEF     (FA%JPLDAT)                 ! IN   
LOGICAL                LDMODA                                 !   OUT

IF (.NOT. FA_COM_DEFAULT_INIT) CALL NEW_FA_DEFAULT ()

CALL FANDAI_MT                              &
&           (FA, KREP, KRANG, KDATEF, LDMODA)

END SUBROUTINE

SUBROUTINE FANDAI_MT                        &
&           (FA, KREP, KRANG, KDATEF, LDMODA)
USE FA_MOD, ONLY : FA_COM
USE LFI_PRECISION
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

CALL FANDAI_MT64                            &
&           (FA, IREP, IRANG, IDATEF, LDMODA)

KREP       = INT (      IREP, JPLIKM)

END SUBROUTINE

!INTF KREP            OUT                               
!INTF KRANG         IN                                  
!INTF KDATEF        IN    DIMS=FA%JPLDAT                
!INTF LDMODA          OUT                               
