! Oct-2012 P. Marguinaud 64b LFI
! Jan-2011 P. Marguinaud Thread-safe LFI
SUBROUTINE LFILDO_MT64                                           &
&                     (LFI, KREP, KNUMER, KREC, KTAB, KNBLEC,  &
&                      KFACTM, KRETIN)
USE LFIMOD, ONLY : LFICOM
USE PARKIND1, ONLY : JPRB, JPIB, JPIA, JPIM
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION
IMPLICIT NONE
!****
!        Sous-programme charge des Lectures de DOnnees du logiciel LFI,
!     *SAUF* pour les articles d'index de type caractere.
!**
!     Arguments: KREP   (Sortie) ==> Code-reponse ( zero si OK; code-
!                                    reponse du "READ" FORTRAN sinon);
!                KNUMER (Entree) ==> NUMERo d'unite logique FORTRAN;
!                KREC   (Entree) ==> Numero d'enregistrement a lire;
!                KTAB   (Sortie) ==> Zone a lire, de Longueur
!                                    LFI%JPLARD*KFACTM *mots*;
!                KNBLEC (Entree  ==> Compteur de LECtures sur l'unite;
!                       +Sortie)
!                KFACTM (Entree) ==> FACteur Multiplicatif LFI de
!                                    l'unite logique;
!                KRETIN (Sortie) ==> Code-retour interne.
!
!
TYPE(LFICOM) :: LFI
INTEGER (KIND=JPLIKB) KREP, KNUMER, KREC, KNBLEC, KFACTM, KRETIN
!
INTEGER (KIND=JPDBLE)  KTAB (LFI%JPLARD*KFACTM)

!
!        LECTURE .
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('LFILDO_MT64',0,ZHOOK_HANDLE)
!
   IF (KNUMER > 0) THEN
     READ (UNIT=KNUMER,REC=KREC,ERR=901,IOSTAT=KREP) KTAB
   ELSE
!RJ: something fishy here
    CALL ABORT
   ENDIF
!
IF (LFI%LMISOP) THEN
  WRITE (UNIT=LFI%NULOUT,FMT=*)                            &
&          '+++++ LFILDO - READ / ',KNUMER,', REC = ',KREC, &
&          ' +++++'
ENDIF
!
KNBLEC=KNBLEC+1
KRETIN=0
GOTO 1001
!
901 CONTINUE
KRETIN=2
!
1001 CONTINUE
!
IF (LHOOK) CALL DR_HOOK('LFILDO_MT64',1,ZHOOK_HANDLE)
END SUBROUTINE



! Oct-2012 P. Marguinaud 64b LFI
SUBROUTINE LFILDO64                                          &
&           (KREP, KNUMER, KREC, KTAB, KNBLEC, KFACTM, KRETIN)
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
INTEGER (KIND=JPDBLE)  KTAB       (LFI%JPLARD*KFACTM)         !   OUT
INTEGER (KIND=JPLIKB)  KNBLEC                                 ! INOUT
INTEGER (KIND=JPLIKB)  KRETIN                                 !   OUT

IF (.NOT. LFICOM_DEFAULT_INIT) CALL NEW_LFI_DEFAULT ()

CALL LFILDO_MT64                                          &
&           (LFI, KREP, KNUMER, KREC, KTAB, KNBLEC, KFACTM, &
&           KRETIN)

END SUBROUTINE

SUBROUTINE LFILDO                                            &
&           (KREP, KNUMER, KREC, KTAB, KNBLEC, KFACTM, KRETIN)
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
INTEGER (KIND=JPDBLE)  KTAB       (LFI%JPLARD*KFACTM)         !   OUT
INTEGER (KIND=JPLIKM)  KNBLEC                                 ! INOUT
INTEGER (KIND=JPLIKM)  KRETIN                                 !   OUT

IF (.NOT. LFICOM_DEFAULT_INIT) CALL NEW_LFI_DEFAULT ()

CALL LFILDO_MT                                            &
&           (LFI, KREP, KNUMER, KREC, KTAB, KNBLEC, KFACTM, &
&           KRETIN)

END SUBROUTINE

SUBROUTINE LFILDO_MT                                      &
&           (LFI, KREP, KNUMER, KREC, KTAB, KNBLEC, KFACTM, &
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
INTEGER (KIND=JPDBLE)  KTAB       (LFI%JPLARD*KFACTM)         !   OUT
INTEGER (KIND=JPLIKM)  KNBLEC                                 ! INOUT
INTEGER (KIND=JPLIKM)  KRETIN                                 !   OUT
! Local integers
INTEGER (KIND=JPLIKB)  IREP                                   !   OUT
INTEGER (KIND=JPLIKB)  INUMER                                 ! IN   
INTEGER (KIND=JPLIKB)  IREC                                   ! IN   
INTEGER (KIND=JPLIKB)  INBLEC                                 ! INOUT
INTEGER (KIND=JPLIKB)  IFACTM                                 ! IN   
INTEGER (KIND=JPLIKB)  IRETIN                                 !   OUT
! Convert arguments

INUMER     = INT (    KNUMER, JPLIKB)
IREC       = INT (      KREC, JPLIKB)
INBLEC     = INT (    KNBLEC, JPLIKB)
IFACTM     = INT (    KFACTM, JPLIKB)

CALL LFILDO_MT64                                          &
&           (LFI, IREP, INUMER, IREC, KTAB, INBLEC, IFACTM, &
&           IRETIN)

KREP       = INT (      IREP, JPLIKM)
KNBLEC     = INT (    INBLEC, JPLIKM)
KRETIN     = INT (    IRETIN, JPLIKM)

END SUBROUTINE

!INTF KREP            OUT                                                              
!INTF KNUMER        IN                                                                 
!INTF KREC          IN                                                                 
!INTF KTAB            OUT DIMS=LFI%JPLARD*KFACTM         KIND=JPDBLE                   
!INTF KNBLEC        INOUT                                                              
!INTF KFACTM        IN                                                                 
!INTF KRETIN          OUT                                                              
