C Jan-2011 P. Marguinaud Thread-safe LFI
      SUBROUTINE LFIEMS_MT (LFI, KNUMER, KNIMES, KCODE, LDFATA, 
     S                      CDMESS, CDNSPR, CDACTI )
      USE LFIMOD, ONLY : LFICOM
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
C****
C        CE SOUS-PROGRAMME EST CHARGE DE FAIRE L'ECHO DES MESSAGES
C     EMIS PAR LE LOGICIEL DE FICHIERS INDEXES LFI, EN FAISANT SI
C     BESOIN EST L'"ABORT" DU PROGRAMME .
C        En l'occurrence, il s'agit d'un "chapeau" qui aiguille sur
C     LFIEFR ou LFIENG en fonction de la variable logique LFI%LFRANC.
C**
C        ARGUMENTS : KNUMER ==> Numero eventuel de l'Unite Logique;
C        ( tous                 ( si LFI%JPNIL ==> pas d'Unite Logique )
C         d'Entree ) KNIMES ==> Niveau (0,1,2) du Message;
C                    KCODE  ==> CODE CORRESPONDANT A L'ACTION;
C                    LDFATA ==> VRAI SI ON DOIT AVORTER LE PROGRAMME;
C                    CDMESS ==> SI KNIMES#0, MESSAGE A EMETTRE;
C                    CDNSPR ==> NOM DU SOUS-PROGRAMME APPELANT;
C                    CDACTI ==> NOM DE L'ACTION D'ENTREE/SORTIE FORTRAN
C                               (SI KCODE >0), SINON FOURRE-TOUT (!) .
C*
C        Pour la table des codes-reponses possibles, voir LFIEFR/LFIENG.
C
C
      TYPE(LFICOM) :: LFI
      INTEGER KNUMER, KNIMES, KCODE, ICODE, IREPON
C
      LOGICAL LDFATA, LLEXUL
C
      CHARACTER  CDNSPR*(*), CDMESS*(*), CDACTI*(*)
#include "lficom_mt.h"
C**
C     1.  -  MODIFICATION EVENTUELLE DU CODE-REPONSE S'IL VAUT (-1).
C-----------------------------------------------------------------------
C*
C        Il s'agit en effet de discriminer entre un numero d'unite
C     logique licite pour le FORTRAN, mais effectivement non ouvert pour
C     le logiciel LFI, auquel cas le code-reponse est laisse a (-1),
C     et un numero d'unite logique FORTRAN carrement illicite, que l'on
C     traduit par le code-reponse (-30).
C
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LFIEMS_MT',0,ZHOOK_HANDLE)
      IF (KCODE.EQ.-1) THEN
        ICODE=-30
        INQUIRE (UNIT=KNUMER,EXIST=LLEXUL,ERR=101,IOSTAT=IREPON)
        IF (LLEXUL) ICODE=KCODE
      ELSE
        ICODE=KCODE
      ENDIF
C
  101 CONTINUE
C**
C     2.  -  APPEL AU SOUS-PROGRAMME AD HOC EN FONCTION DE *LFI%LFRANC*.
C-----------------------------------------------------------------------
C
      IF (LFI%LFRANC) THEN
        CALL LFIEFR_MT (LFI, KNUMER,KNIMES,ICODE,LDFATA,
     S                  CDMESS,CDNSPR,CDACTI)
      ELSE
        CALL LFIENG_MT (LFI, KNUMER,KNIMES,ICODE,LDFATA,
     S                  CDMESS,CDNSPR,CDACTI)
      ENDIF
C
      IF (LHOOK) CALL DR_HOOK('LFIEMS_MT',1,ZHOOK_HANDLE)
      END

