! Oct-2012 P. Marguinaud 64b LFI
! Jan-2011 P. Marguinaud Thread-safe LFI
SUBROUTINE LFIECC_MT64                                            &
&                     (LFI, KREP, KNUMER, KREC, CDTAB, KNBECR,  &
&                      KFACTM, KRETIN)
USE LFIMOD, ONLY : LFICOM
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION
IMPLICIT NONE
!****
!        Sous-programme charge des Ecritures de Chaines de Caracteres
!     du logiciel LFI (articles d'index "noms").
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
CHARACTER CDTAB (LFI%JPNXNA*KFACTM)*(LFI%JPNCPN)
!
!        ECRITURE .
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('LFIECC_MT64',0,ZHOOK_HANDLE)
!
   IF (KNUMER > 0) THEN
     WRITE (UNIT=KNUMER,REC=KREC,ERR=901,IOSTAT=KREP) CDTAB
   ELSE
!RJ: something fishy here
     CALL ABORT
   ENDIF
!
IF (LFI%LMISOP) THEN
  WRITE (UNIT=LFI%NULOUT,FMT=*)                             &
&          '+++++ LFIECC - WRITE / ',KNUMER,', REC = ',KREC, &
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
IF (LHOOK) CALL DR_HOOK('LFIECC_MT64',1,ZHOOK_HANDLE)
END SUBROUTINE



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE LFIECC64                                           &
&           (KREP, KNUMER, KREC, CDTAB, KNBECR, KFACTM, KRETIN)
USE LFIMOD, ONLY : LFI => LFICOM_DEFAULT, &
&                   LFICOM_DEFAULT_INIT,   &
&                   NEW_LFI_DEFAULT
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKB)  KREP                                   !   OUT
INTEGER (KIND=JPLIKB)  KNUMER                                 ! IN   
INTEGER (KIND=JPLIKB)  KREC                                   ! IN   
INTEGER (KIND=JPLIKB)  KFACTM                                 ! IN   
CHARACTER (LEN=*)      CDTAB      (LFI%JPNXNA*KFACTM)         !   OUT
INTEGER (KIND=JPLIKB)  KNBECR                                 ! INOUT
INTEGER (KIND=JPLIKB)  KRETIN                                 !   OUT

IF (.NOT. LFICOM_DEFAULT_INIT) CALL NEW_LFI_DEFAULT ()

CALL LFIECC_MT64                                           &
&           (LFI, KREP, KNUMER, KREC, CDTAB, KNBECR, KFACTM, &
&           KRETIN)

END SUBROUTINE

SUBROUTINE LFIECC                                             &
&           (KREP, KNUMER, KREC, CDTAB, KNBECR, KFACTM, KRETIN)
USE LFIMOD, ONLY : LFI => LFICOM_DEFAULT, &
&                   LFICOM_DEFAULT_INIT,   &
&                   NEW_LFI_DEFAULT
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KNUMER                                 ! IN   
INTEGER (KIND=JPLIKM)  KREC                                   ! IN   
INTEGER (KIND=JPLIKM)  KFACTM                                 ! IN   
CHARACTER (LEN=*)      CDTAB      (LFI%JPNXNA*KFACTM)         !   OUT
INTEGER (KIND=JPLIKM)  KNBECR                                 ! INOUT
INTEGER (KIND=JPLIKM)  KRETIN                                 !   OUT

IF (.NOT. LFICOM_DEFAULT_INIT) CALL NEW_LFI_DEFAULT ()

CALL LFIECC_MT                                             &
&           (LFI, KREP, KNUMER, KREC, CDTAB, KNBECR, KFACTM, &
&           KRETIN)

END SUBROUTINE

SUBROUTINE LFIECC_MT                                       &
&           (LFI, KREP, KNUMER, KREC, CDTAB, KNBECR, KFACTM, &
&           KRETIN)
USE LFIMOD, ONLY : LFICOM
USE LFI_PRECISION
IMPLICIT NONE
! Arguments
TYPE (LFICOM)          LFI                                    ! INOUT
INTEGER (KIND=JPLIKM)  KREP                                   !   OUT
INTEGER (KIND=JPLIKM)  KNUMER                                 ! IN   
INTEGER (KIND=JPLIKM)  KREC                                   ! IN   
INTEGER (KIND=JPLIKM)  KFACTM                                 ! IN   
CHARACTER (LEN=*)      CDTAB      (LFI%JPNXNA*KFACTM)         !   OUT
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

CALL LFIECC_MT64                                           &
&           (LFI, IREP, INUMER, IREC, CDTAB, INBECR, IFACTM, &
&           IRETIN)

KREP       = INT (      IREP, JPLIKM)
KNBECR     = INT (    INBECR, JPLIKM)
KRETIN     = INT (    IRETIN, JPLIKM)

END SUBROUTINE

!INTF KREP            OUT                                                              
!INTF KNUMER        IN                                                                 
!INTF KREC          IN                                                                 
!INTF CDTAB           OUT DIMS=LFI%JPNXNA*KFACTM         LEN=LFI%JPNCPN                
!INTF KNBECR        INOUT                                                              
!INTF KFACTM        IN                                                                 
!INTF KRETIN          OUT                                                              
