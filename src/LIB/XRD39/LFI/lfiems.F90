! Oct-2012 P. Marguinaud 64b LFI
#include "lfisuffix.h"
! Jan-2011 P. Marguinaud Thread-safe LFI
SUBROUTINE LFIEMS_MT_BB                                   &
&                     (LFI, KNUMER, KNIMES, KCODE, LDFATA,  &
&                      CDMESS, CDNSPR, CDACTI )
USE LFIMOD, ONLY : LFICOM
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
!****
!        CE SOUS-PROGRAMME EST CHARGE DE FAIRE L'ECHO DES MESSAGES
!     EMIS PAR LE LOGICIEL DE FICHIERS INDEXES LFI, EN FAISANT SI
!     BESOIN EST L'"ABORT" DU PROGRAMME .
!        En l'occurrence, il s'agit d'un "chapeau" qui aiguille sur
!     LFIEFR ou LFIENG en fonction de la variable logique LFI%LFRANC.
!**
!        ARGUMENTS : KNUMER ==> Numero eventuel de l'Unite Logique;
!        ( tous                 ( si LFI%JPNIL ==> pas d'Unite Logique )
!         d'Entree ) KNIMES ==> Niveau (0,1,2) du Message;
!                    KCODE  ==> CODE CORRESPONDANT A L'ACTION;
!                    LDFATA ==> VRAI SI ON DOIT AVORTER LE PROGRAMME;
!                    CDMESS ==> SI KNIMES#0, MESSAGE A EMETTRE;
!                    CDNSPR ==> NOM DU SOUS-PROGRAMME APPELANT;
!                    CDACTI ==> NOM DE L'ACTION D'ENTREE/SORTIE FORTRAN
!                               (SI KCODE >0), SINON FOURRE-TOUT (!) .
!*
!        Pour la table des codes-reponses possibles, voir LFIEFR/LFIENG.
!
!
TYPE(LFICOM) :: LFI
INTEGER (KIND=JPLIKB) KNUMER, KNIMES, KCODE, ICODE, IREPON
!
LOGICAL LDFATA, LLEXUL
!
CHARACTER  CDNSPR*(*), CDMESS*(*), CDACTI*(*)

!**
!     1.  -  MODIFICATION EVENTUELLE DU CODE-REPONSE S'IL VAUT (-1).
!-----------------------------------------------------------------------
!*
!        Il s'agit en effet de discriminer entre un numero d'unite
!     logique licite pour le FORTRAN, mais effectivement non ouvert pour
!     le logiciel LFI, auquel cas le code-reponse est laisse a (-1),
!     et un numero d'unite logique FORTRAN carrement illicite, que l'on
!     traduit par le code-reponse (-30).
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('LFIEMS_MT',0,ZHOOK_HANDLE)
IF (KCODE.EQ.-1) THEN
  ICODE=-30
  INQUIRE (UNIT=KNUMER,EXIST=LLEXUL,ERR=101,IOSTAT=IREPON)
  IF (LLEXUL) ICODE=KCODE
ELSE
  ICODE=KCODE
ENDIF
!
101 CONTINUE
!**
!     2.  -  APPEL AU SOUS-PROGRAMME AD HOC EN FONCTION DE *LFI%LFRANC*.
!-----------------------------------------------------------------------
!
IF (LFI%LFRANC) THEN
  CALL LFIEFR_MT_BB                               &
&                 (LFI, KNUMER,KNIMES,ICODE,LDFATA, &
&                  CDMESS,CDNSPR,CDACTI)
ELSE
  CALL LFIENG_MT_BB                               &
&                 (LFI, KNUMER,KNIMES,ICODE,LDFATA, &
&                  CDMESS,CDNSPR,CDACTI)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('LFIEMS_MT',1,ZHOOK_HANDLE)
END SUBROUTINE



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE LFIEMS_BB                                      &
&           (KNUMER, KNIMES, KCODE, LDFATA, CDMESS, CDNSPR, &
&           CDACTI)
USE LFIMOD, ONLY : LFI => LFICOM_DEFAULT, &
&                   LFICOM_DEFAULT_INIT,   &
&                   NEW_LFI_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKB)  KNUMER                                 ! IN   
INTEGER (KIND=JPLIKB)  KNIMES                                 ! IN   
INTEGER (KIND=JPLIKB)  KCODE                                  ! IN   
LOGICAL                LDFATA                                 ! IN   
CHARACTER (LEN=*)      CDMESS                                 ! IN   
CHARACTER (LEN=*)      CDNSPR                                 ! IN   
CHARACTER (LEN=*)      CDACTI                                 ! IN   

IF (.NOT. LFICOM_DEFAULT_INIT) CALL NEW_LFI_DEFAULT ()

CALL LFIEMS_MT_BB                                              &
&           (LFI, KNUMER, KNIMES, KCODE, LDFATA, CDMESS, CDNSPR, &
&           CDACTI)

END SUBROUTINE

SUBROUTINE LFIEMS_MM                                      &
&           (KNUMER, KNIMES, KCODE, LDFATA, CDMESS, CDNSPR, &
&           CDACTI)
USE LFIMOD, ONLY : LFI => LFICOM_DEFAULT, &
&                   LFICOM_DEFAULT_INIT,   &
&                   NEW_LFI_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKM)  KNUMER                                 ! IN   
INTEGER (KIND=JPLIKM)  KNIMES                                 ! IN   
INTEGER (KIND=JPLIKM)  KCODE                                  ! IN   
LOGICAL                LDFATA                                 ! IN   
CHARACTER (LEN=*)      CDMESS                                 ! IN   
CHARACTER (LEN=*)      CDNSPR                                 ! IN   
CHARACTER (LEN=*)      CDACTI                                 ! IN   

IF (.NOT. LFICOM_DEFAULT_INIT) CALL NEW_LFI_DEFAULT ()

CALL LFIEMS_MT_MM                                              &
&           (LFI, KNUMER, KNIMES, KCODE, LDFATA, CDMESS, CDNSPR, &
&           CDACTI)

END SUBROUTINE

SUBROUTINE LFIEMS_MT_MM                                        &
&           (LFI, KNUMER, KNIMES, KCODE, LDFATA, CDMESS, CDNSPR, &
&           CDACTI)
USE LFIMOD, ONLY : LFICOM,              &
&                   LFICOM_DEFAULT_INIT, &
&                   NEW_LFI_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
TYPE (LFICOM)          LFI                                    ! INOUT
INTEGER (KIND=JPLIKM)  KNUMER                                 ! IN   
INTEGER (KIND=JPLIKM)  KNIMES                                 ! IN   
INTEGER (KIND=JPLIKM)  KCODE                                  ! IN   
LOGICAL                LDFATA                                 ! IN   
CHARACTER (LEN=*)      CDMESS                                 ! IN   
CHARACTER (LEN=*)      CDNSPR                                 ! IN   
CHARACTER (LEN=*)      CDACTI                                 ! IN   
! Local integers
INTEGER (KIND=JPLIKB)  INUMER                                 ! IN   
INTEGER (KIND=JPLIKB)  INIMES                                 ! IN   
INTEGER (KIND=JPLIKB)  ICODE                                  ! IN   
! Convert arguments

INUMER     = INT (    KNUMER, JPLIKB)
INIMES     = INT (    KNIMES, JPLIKB)
ICODE      = INT (     KCODE, JPLIKB)

CALL LFIEMS_MT_BB                                              &
&           (LFI, INUMER, INIMES, ICODE, LDFATA, CDMESS, CDNSPR, &
&           CDACTI)


END SUBROUTINE

!INTF KNUMER        IN    
!INTF KNIMES        IN    
!INTF KCODE         IN    
!INTF LDFATA        IN    
!INTF CDMESS        IN    
!INTF CDNSPR        IN    
!INTF CDACTI        IN    
