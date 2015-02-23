! Oct-2012 P. Marguinaud 64b LFI
#include "lfisuffix.h"
! Jan-2011 P. Marguinaud Thread-safe LFI
SUBROUTINE LFIEDO_MT_BB                                      &
&                     (LFI, KREP, KNUMER, KREC, KTAB, KNBECR,  &
&                      KFACTM, KRETIN )
USE LFIMOD, ONLY : LFICOM
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
#ifdef USE_SAMIO
USE SAMIO_MOD
#endif
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
!****
!        Sous-programme charge des Ecritures de DOnnees du logiciel LFI,
!     *SAUF* pour les articles d'index de type caractere.
!**
!     Arguments: KREP   (Sortie) ==> Code-reponse ( zero si OK; code-
!                                    reponse du "WRITE" FORTRAN sinon);
!                KNUMER (Entree) ==> NUMERo d'unite logique FORTRAN;
!                KREC   (Entree) ==> Numero d'enregistrement a ecrire;
!                KTAB   (Sortie) ==> Zone a ecrire, de Longueur
!                                    LFI%JPLARD*KFACTM *mots*;
!                KNBECR (Entree  ==> Compteur d'ECRitures sur l'unite;
!                       +Sortie)
!                KFACTM (Entree) ==> FACteur Multiplicatif LFI de
!                                    l'unite logique;
!                KRETIN (Sortie) ==> Code-retour interne.
!
!
TYPE(LFICOM) :: LFI
INTEGER (KIND=JPLIKB) KREP, KNUMER, KREC
INTEGER (KIND=JPLIKB) KNBECR, KFACTM, KRETIN
!
INTEGER (KIND=JPDBLE)  KTAB (LFI%JPLARD*KFACTM)
LOGICAL :: LLSAMIO

!
!        ECRITURE .
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('LFIEDO_MT',0,ZHOOK_HANDLE)
LLSAMIO = .FALSE.
#ifdef USE_SAMIO
LLSAMIO = SAMIO_HAS_OPENED(KNUMER)
IF (LLSAMIO) THEN
   CALL SAMIO_WRITE (UNIT=KNUMER,REC=KREC,ERR=901,IOSTAT=KREP,  &
&        ARRAY=KTAB)
   IF (KREP /= 0) GOTO 901
ELSE
#endif
   WRITE (UNIT=KNUMER,REC=KREC,ERR=901,IOSTAT=KREP) KTAB
#ifdef USE_SAMIO
ENDIF
#endif
!
IF (LFI%LMISOP) THEN
  WRITE (UNIT=LFI%NULOUT,FMT=*)                             &
&          '+++++ LFIEDO - WRITE / ',KNUMER,', REC = ',KREC, &
&          ' +++++'
ENDIF
!
KNBECR=KNBECR+1
KRETIN=0
GOTO 1001
!
901 CONTINUE
KRETIN=1
!
1001 CONTINUE
!
IF (LHOOK) CALL DR_HOOK('LFIEDO_MT',1,ZHOOK_HANDLE)
END SUBROUTINE



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE LFIEDO_BB                                         &
&           (KREP, KNUMER, KREC, KTAB, KNBECR, KFACTM, KRETIN)
USE LFIMOD, ONLY : LFI => LFICOM_DEFAULT, &
&                   LFICOM_DEFAULT_INIT,   &
&                   NEW_LFI_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKB)  KREP                                   !   OUT
INTEGER (KIND=JPLIKB)  KNUMER                                 ! IN   
INTEGER (KIND=JPLIKB)  KREC                                   ! IN   
INTEGER (KIND=JPLIKB)  KFACTM                                 ! IN   
INTEGER (KIND=JPDBLE)  KTAB       (*)         !   OUT
INTEGER (KIND=JPLIKB)  KNBECR                                 ! INOUT
INTEGER (KIND=JPLIKB)  KRETIN                                 !   OUT

IF (.NOT. LFICOM_DEFAULT_INIT) CALL NEW_LFI_DEFAULT ()

CALL LFIEDO_MT_BB                                         &
&           (LFI, KREP, KNUMER, KREC, KTAB, KNBECR, KFACTM, &
&           KRETIN)

END SUBROUTINE

SUBROUTINE LFIEDO_MM                                         &
&           (KREP, KNUMER, KREC, KTAB, KNBECR, KFACTM, KRETIN)
USE LFIMOD, ONLY : LFI => LFICOM_DEFAULT, &
&                   LFICOM_DEFAULT_INIT,   &
&                   NEW_LFI_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KNUMER                                 ! IN   
INTEGER (KIND=JPLIKM)  KREC                                   ! IN   
INTEGER (KIND=JPLIKM)  KFACTM                                 ! IN   
INTEGER (KIND=JPDBLE)  KTAB       (*)         !   OUT
INTEGER (KIND=JPLIKM)  KNBECR                                 ! INOUT
INTEGER (KIND=JPLIKM)  KRETIN                                 !   OUT

IF (.NOT. LFICOM_DEFAULT_INIT) CALL NEW_LFI_DEFAULT ()

CALL LFIEDO_MT_MM                                         &
&           (LFI, KREP, KNUMER, KREC, KTAB, KNBECR, KFACTM, &
&           KRETIN)

END SUBROUTINE

SUBROUTINE LFIEDO_MT_MM                                   &
&           (LFI, KREP, KNUMER, KREC, KTAB, KNBECR, KFACTM, &
&           KRETIN)
USE LFIMOD, ONLY : LFICOM,              &
&                   LFICOM_DEFAULT_INIT, &
&                   NEW_LFI_DEFAULT
USE LFI_PRECISION, ONLY : JPDBLE, JPDBLR, JPLIKB, JPLIKM
IMPLICIT NONE
! Arguments
TYPE (LFICOM)          LFI                                    ! INOUT
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KNUMER                                 ! IN   
INTEGER (KIND=JPLIKM)  KREC                                   ! IN   
INTEGER (KIND=JPLIKM)  KFACTM                                 ! IN   
INTEGER (KIND=JPDBLE)  KTAB       (LFI%JPLARD*KFACTM)         !   OUT
INTEGER (KIND=JPLIKM)  KNBECR                                 ! INOUT
INTEGER (KIND=JPLIKM)  KRETIN                                 !   OUT
! Local integers
INTEGER (KIND=JPLIKB)  IREP                                   !   OUT
INTEGER (KIND=JPLIKB)  INUMER                                 ! IN   
INTEGER (KIND=JPLIKB)  IREC                                   ! IN   
INTEGER (KIND=JPLIKB)  INBECR                                 ! INOUT
INTEGER (KIND=JPLIKB)  IFACTM                                 ! IN   
INTEGER (KIND=JPLIKB)  IRETIN                                 !   OUT
! Convert arguments

INUMER     = INT (    KNUMER, JPLIKB)
IREC       = INT (      KREC, JPLIKB)
INBECR     = INT (    KNBECR, JPLIKB)
IFACTM     = INT (    KFACTM, JPLIKB)

CALL LFIEDO_MT_BB                                         &
&           (LFI, IREP, INUMER, IREC, KTAB, INBECR, IFACTM, &
&           IRETIN)

KREP       = INT (      IREP, JPLIKM)
KNBECR     = INT (    INBECR, JPLIKM)
KRETIN     = INT (    IRETIN, JPLIKM)

END SUBROUTINE

!INTF KREP            OUT                                                              
!INTF KNUMER        IN                                                                 
!INTF KREC          IN                                                                 
!INTF KTAB            OUT DIMS=LFI%JPLARD*KFACTM         KIND=JPDBLE                   
!INTF KNBECR        INOUT                                                              
!INTF KFACTM        IN                                                                 
!INTF KRETIN          OUT                                                              
