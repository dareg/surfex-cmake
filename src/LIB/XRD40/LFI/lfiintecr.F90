! Oct-2012 P. Marguinaud 64b LFI
#include "lfisuffix.h"
! Jan-2011 P. Marguinaud Thread-safe LFI
SUBROUTINE LFIINTECR_MT_BB                                       &
&                        (LFI, KREP, KNUMER, CDNOMA, KTAB, KLONG )
USE LFIMOD, ONLY : LFICOM
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
!
!****
!        SOUS-PROGRAMME D'ECRITURE  D'UN TABLEAU D'ENTIERS A PARTIR 
!     D'UN ARTICLE (DE DONNEES) PAR *NOM*
!     SUR UNE UNITE LOGIQUE OUVERTE POUR LE LOGICIEL DE FICHIERS INDEXES
!     *LFI*; L'ARTICLE EN SORTIE EST UN "BLOC" DE DONNEES ADJACENTES.
!**
!    ARGUMENTS : KREP   (SORTIE) ==> CODE-REPONSE DU SOUS-PROGRAMME;
!                KNUMER (ENTREE) ==> LFI%NUMERO DE L'UNITE LOGIQUE;
!                CDNOMA (ENTREE) ==> NOM DE L'ARTICLE A RECHERCHER;
!                KTAB   (ENTREE) ==> PREMIER MOT A ECRIRE
!                KLONG  (ENTREE) ==> LONGUEUR DE L'ARTICLE A LIRE.
!
!
TYPE(LFICOM) :: LFI
CHARACTER CDNOMA*(*)
INTEGER (KIND=JPLIKB) KREP, KNUMER, KLONG, KTAB(KLONG)

INTEGER (KIND=JPLIKB) JI
INTEGER (KIND=JPDBLE)  ITAB (KLONG)

!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('LFIINTECR_MT',0,ZHOOK_HANDLE)
DO JI=1,KLONG
  ITAB(JI)=INT(KTAB(JI),KIND=JPDBLE)
ENDDO
CALL LFIECR_MT_BB                                        &
&               (LFI,  KREP, KNUMER, CDNOMA, ITAB, KLONG )
!
IF (LHOOK) CALL DR_HOOK('LFIINTECR_MT',1,ZHOOK_HANDLE)
END SUBROUTINE



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE LFIINTECR_BB                       &
&           (KREP, KNUMER, CDNOMA, KTAB, KLONG)
USE LFIMOD, ONLY : LFI => LFICOM_DEFAULT, &
&                   LFICOM_DEFAULT_INIT,   &
&                   NEW_LFI_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKB)  KREP                                   !   OUT
INTEGER (KIND=JPLIKB)  KNUMER                                 ! IN   
CHARACTER (LEN=*)      CDNOMA                                 ! IN   
INTEGER (KIND=JPLIKB)  KLONG                                  ! IN   
INTEGER (KIND=JPLIKB)  KTAB       (KLONG)                     ! IN   

IF (.NOT. LFICOM_DEFAULT_INIT) CALL NEW_LFI_DEFAULT ()

CALL LFIINTECR_MT_BB                               &
&           (LFI, KREP, KNUMER, CDNOMA, KTAB, KLONG)

END SUBROUTINE

SUBROUTINE LFIINTECR_MM                       &
&           (KREP, KNUMER, CDNOMA, KTAB, KLONG)
USE LFIMOD, ONLY : LFI => LFICOM_DEFAULT, &
&                   LFICOM_DEFAULT_INIT,   &
&                   NEW_LFI_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KNUMER                                 ! IN   
CHARACTER (LEN=*)      CDNOMA                                 ! IN   
INTEGER (KIND=JPLIKM)  KLONG                                  ! IN   
INTEGER (KIND=JPLIKM)  KTAB       (KLONG)                     ! IN   

IF (.NOT. LFICOM_DEFAULT_INIT) CALL NEW_LFI_DEFAULT ()

CALL LFIINTECR_MT_MM                               &
&           (LFI, KREP, KNUMER, CDNOMA, KTAB, KLONG)

END SUBROUTINE

SUBROUTINE LFIINTECR_MT_MM                         &
&           (LFI, KREP, KNUMER, CDNOMA, KTAB, KLONG)
USE LFIMOD, ONLY : LFICOM,              &
&                   LFICOM_DEFAULT_INIT, &
&                   NEW_LFI_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
TYPE (LFICOM)          LFI                                    ! INOUT
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KNUMER                                 ! IN   
CHARACTER (LEN=*)      CDNOMA                                 ! IN   
INTEGER (KIND=JPLIKM)  KLONG                                  ! IN   
INTEGER (KIND=JPLIKM)  KTAB       (KLONG)                     ! IN   
! Local integers
INTEGER (KIND=JPLIKB)  IREP                                   !   OUT
INTEGER (KIND=JPLIKB)  INUMER                                 ! IN   
INTEGER (KIND=JPLIKB)  ITAB       (KLONG)                     ! IN   
INTEGER (KIND=JPLIKB)  ILONG                                  ! IN   
! Convert arguments

INUMER     = INT (    KNUMER, JPLIKB)
ITAB       = INT (      KTAB, JPLIKB)
ILONG      = INT (     KLONG, JPLIKB)

CALL LFIINTECR_MT_BB                               &
&           (LFI, IREP, INUMER, CDNOMA, ITAB, ILONG)

KREP       = INT (      IREP, JPLIKM)

END SUBROUTINE

!INTF KREP            OUT                               
!INTF KNUMER        IN                                  
!INTF CDNOMA        IN                                  
!INTF KTAB          IN    DIMS=KLONG                    
!INTF KLONG         IN                                  
